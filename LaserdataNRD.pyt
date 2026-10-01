# -*- coding: utf-8 -*-
"""
LaserdataNRD.pyt

Beräknar relativ punkttäthet i höjdskikt, NRD eller ORD, som raster för ett
intresseområde ur Lantmäteriets Laserdata Nedladdning, skog. Måtten följer
Campbell m.fl. (2018), Quantifying understory vegetation density using
small-footprint airborne lidar, Remote Sensing of Environment 215:330-342,
https://doi.org/10.1016/j.rse.2018.06.023:

    NRD(i, j) = antal punkter med i <= h < j / antal punkter med h < j
    ORD(i, j) = antal punkter med i <= h < j / antal punkter totalt

där h är punktens höjd över mark. NRD räknar bara de punkter som nått ned
till skiktet eller längre, och påverkas därför mindre av ett tätt krontak än
ORD. I studien förklarade NRD 0,15-1,85 m 44 % av variationen i fältmätt
buskskiktstäthet, ORD 14 %.

Datakälla
---------
STAC-katalog (öppen, ingen inloggning):
    https://api.lantmateriet.se/stac-hojd/v1
    Samling: dsm-skoglig-copc (den enda punktmolnssamlingen; övriga samlingar
             är markmodeller som raster)
    Sök:     POST /search  {"collections": [...], "bbox": [...], "limit": 100}
             Nästa sida via länken rel="next" (method + body).

Varje item är en ruta på 10 x 10 km (eller 5 x 10 km) i SWEREF 99 TM + RH 2000
(EPSG:5845), med en asset "data" som pekar på en COPC-fil (.copc.laz, ~1 GB).
properties: proj:bbox (rutans hörn i SWEREF 99 TM), pc:count, datetime,
start_datetime, end_datetime, skanningsomrade (t.ex. 26A024), punkttathet,
flyghojd. Många områden är skannade två gånger och båda finns kvar i
katalogen, med olika skanningsområde. Verktyget tar som standard den senaste
skanningen per ruta, men användaren kan välja skanning själv.

Nedladdning kräver OAuth2 (client credentials):
    POST https://apimanager.lantmateriet.se/oauth2/token
    Basic-auth med consumer key/secret, grant_type=client_credentials.
    Token gäller 3600 s. Nyckeln måste ha både API:et STAC-hojd och en
    beställning av Laserdata Nedladdning, skog, annars svarar dl1 med 403.

Hela rutor laddas aldrig ned. PDAL (ingår i ArcGIS Pro) läser COPC-filerna
direkt över HTTP och hämtar bara punkterna inom områdets utbredning. Området
delas i block som vart och ett ryms i minnesbudgeten. Blocken bearbetas i egna
processer (nrd_worker.py), flera åt gången, och varje process returnerar bara
antal punkter per cell och skikt. Rastren byggs sedan i ett svep med numpy,
eftersom cellerna är stora (standard 20 m) och räknetabellerna små.

Klasser: 1 oklassad, 2 mark, 7 lågt brus, 18 högt brus. Brus tas bort. Höjd
över mark räknas från en triangulering av markpunkterna, se nrd_worker.py.

Verktygstips (parameterförklaringar) skrivs till
LaserdataNRD.PunkttathetISkikt.pyt.xml från TOOLTIPS nedan när verktygslådan
laddas, så att texten bara finns på ett ställe.

Krav: ArcGIS Pro 3.x. arcpy, numpy, pdal och certifi ingår i arcgispro-py3.
Ingen licensnivå utöver Basic behövs.
"""

import base64
import datetime
import importlib
import importlib.util
import json
import math
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from xml.sax.saxutils import escape

import arcpy
import numpy as np

# PDAL-delen ligger i nrd_worker.py bredvid verktygslådan och körs i egna
# processer. reload: Pro håller modulen i minnet mellan körningar.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import nrd_worker as _nw  # noqa: E402
importlib.reload(_nw)

# ── Konstanter ────────────────────────────────────────────────────────────────

TOOL_VERSION = "1.0"
REPO_NAME = "arcgis-laserdata-nrd"

STAC_SEARCH_URL = "https://api.lantmateriet.se/stac-hojd/v1/search"
TOKEN_URL = "https://apimanager.lantmateriet.se/oauth2/token"
COLLECTION = "dsm-skoglig-copc"
GEOTORGET_URL = "https://geotorget.lantmateriet.se/geodataprodukter/laserdata-nedladdning-skog-api"
TERMS_URL = ("https://www.lantmateriet.se/globalassets/geodata/geodataprodukter/"
             "anvandningsvillkor-for-laserdata-nedladdning-skog.pdf")
PAPER_DOI = "https://doi.org/10.1016/j.rse.2018.06.023"

USER_AGENT = "arcgis-laserdata-nrd/" + TOOL_VERSION
HTTP_TIMEOUT = 60
HTTP_RETRIES = 3

# Token gäller 3600 s. Hämta en ny innan nästa block om den är äldre än så här.
TOKEN_MAX_AGE_S = 50 * 60

SWEREF99TM_WKID = 3006
RH2000_WKID = 5613

# Varje block läses med denna marginal, så att trianguleringen av markpunkterna
# har fullt underlag vid blockkanten.
READ_MARGIN_M = 20.0

# Minnestopp per läst punkt i en arbetsprocess: PDAL:s egen kopia och
# numpy-arrayen, arbetskopian (X, Y, Z, klass som float64) och höjd- och
# indexarrayerna. Uppmätt 213-220 byte med 0,9-5,4 miljoner punkter och 8
# attribut; Lantmäteriets filer har fler attribut, därav marginalen.
PEAK_BYTES_PER_POINT = 260
# Fast minne per arbetsprocess (Python, numpy, PDAL), uppmätt ~120 MB.
WORKER_BASE_BYTES = 150e6

MIN_BLOCK_M = 250
MAX_BLOCK_M = 5000

METRIC_NRD = "NRD (normaliserad relativ punkttäthet)"
METRIC_ORD = "ORD (total relativ punkttäthet)"
METRIC_CODE = {METRIC_NRD: "NRD", METRIC_ORD: "ORD"}

DEFAULT_BANDS = [[0.15, 1.85], [2.0, 5.0], [5.0, 10.0]]
MAX_BANDS = 20
MAX_BAND_TOP_M = 60.0
DEFAULT_CELL_SIZE = 20.0
DEFAULT_MIN_POINTS = 20
# Automatiskt läge (tomma fält): hälften av det lediga arbetsminnet, men minst
# 2 GB kvar åt Pro och Windows, och en process per processortråd utom en.
# Taket på antal processer skyddar Lantmäteriets server och nätverket; varje
# block öppnar dessutom fjärrfilen (~1 s), så små block lönar sig inte.
AUTO_MEMORY_SHARE = 0.5
AUTO_MEMORY_RESERVE_GB = 2.0
MIN_MEMORY_GB = 0.5
MIN_GB_PER_WORKER = 0.75
MAX_WORKERS = 16
DEFAULT_PREFIX = "laser"

# Varna över denna yta: körningen blir lång.
LARGE_AREA_KM2 = 500

AOI_POLYGONS = "Polygoner i ett lager"
AOI_EXTENT = "Utbredning (kartvy, lager eller koordinater)"

COUNT_NODATA = -1

MONTHS = ["januari", "februari", "mars", "april", "maj", "juni", "juli", "augusti",
          "september", "oktober", "november", "december"]

TOOL_SUMMARY = (
    "Beräknar relativ punkttäthet i valfria höjdskikt, NRD eller ORD, som raster för ett "
    "intresseområde ur Lantmäteriets Laserdata Nedladdning, skog. NRD i skiktet 0,15-1,85 m "
    "är ett mått på buskskiktets täthet (Campbell m.fl. 2018). Bara punkterna inom områdets "
    "utbredning hämtas, hela rutor laddas aldrig ned."
)

# Verktygstips per parameter, visas i verktygsdialogen. Se _write_tool_metadata.
TOOLTIPS = {
    "aoi_mode": (
        "Hur området avgränsas. 'Polygoner i ett lager' klipper resultatet till själva "
        "polygonerna. 'Utbredning' ger en rektangel: aktuell kartvy, utbredningen av ett "
        "lager, en ritad rektangel eller inskrivna koordinater."
    ),
    "aoi": (
        "Polygonlager som avgränsar området. Alla objekt i lagret slås ihop, eller bara de "
        "valda om det finns ett urval. Lagret kan ha vilket koordinatsystem som helst. "
        "Punkter hämtas inom polygonernas utbredning (bounding box), rastren klipps sedan "
        "till själva polygonerna."
    ),
    "aoi_extent": (
        "Rektangel som avgränsar området. I listan kan du välja kartvyns aktuella "
        "utbredning, eller ett lager för att använda dess utbredning. Du kan också rita en "
        "rektangel i kartan eller skriva in koordinater. Inskrivna koordinater tolkas i den "
        "aktiva kartans koordinatsystem; vilket som användes står i meddelandena."
    ),
    "metric": (
        "Vilket mått som beräknas i varje skikt.\n"
        "NRD = punkter i skiktet / punkter från marken upp till skiktets övre gräns. Räknar "
        "bara pulser som nått ned till skiktet, och påverkas därför mindre av ett tätt "
        "krontak ovanför. Rekommenderas för buskskikt och undervegetation.\n"
        "ORD = punkter i skiktet / alla punkter i cellen. Under ett tätt krontak blir ORD "
        "låg även där undervegetationen är tät, eftersom få pulser når ned.\n"
        "Utan något ovanför skiktet är de två måtten lika. I Campbell m.fl. (2018) förklarade "
        "NRD 0,15-1,85 m 44 % av fältmätt buskskiktstäthet, ORD 14 %."
    ),
    "bands": (
        "Höjdskikten som beräknas, i meter över mark, en rad per skikt. Ett raster skapas "
        "per skikt. Undre gränsen räknas in i skiktet, den övre inte, så en punkt exakt på "
        "gränsen hamnar i skiktet ovanför. Skikten får överlappa. Studien använde 0,15-1,85 m "
        "för buskskiktet. En undre gräns under ungefär 0,15 m fångar markpunkter som klassats "
        "fel eller ligger lite ovanför den triangulerade markytan, eftersom markmodellens "
        "höjdfel är några cm. Med NRD och undre gräns 0 blir värdet nära 1 överallt. Högst "
        "20 skikt, övre gräns högst 60 m."
    ),
    "cell_size": (
        "Rastrens cellstorlek i meter. Laserdata Skog har ungefär 2-5 punkter per m² räknat "
        "med alla returer (den nominella tätheten 1-2 per m² är lägre), en fjärdedel till en "
        "sjättedel av tätheten i studien (16 per m²). Varje cell behöver därför en stor yta "
        "för att få tillräckligt många punkter i nämnaren. Vid 20 m ligger ungefär "
        "1 000-2 000 punkter i en cell, men under tät skog når bara en del av dem ned under "
        "2 m. Studien karterade med 30 m. Under 10 m blir värdena brusiga. 1 till 200 m."
    ),
    "out_workspace": (
        "Geodatabas eller mapp där rastren sparas. I en mapp blir de GeoTIFF."
    ),
    "prefix": (
        "Början på utdatanamnen. Rastren heter <prefix>_<mått>_<cellstorlek>m_<undre>_<övre> "
        "med höjdgränserna i cm, till exempel laser_nrd_20m_015_185 för NRD 0,15-1,85 m i "
        "20 m celler (2,5 m skrivs 2p5m). Punktantalen heter "
        "<prefix>_nrd_<cellstorlek>m_antal_under_<övre> eller "
        "<prefix>_ord_<cellstorlek>m_antal_alla. Körningar med olika mått eller cellstorlek "
        "skriver alltså inte över varandra. Befintliga raster med samma namn skrivs över. "
        "Bara bokstäver, siffror och understreck, och första tecknet måste vara en bokstav."
    ),
    "add_group": (
        "Lägg resultatet i en grupp i den aktiva kartan, med samma namn som prefixet och "
        "måttet, till exempel 'laser - NRD 20 m'. Avmarkerat läggs varje raster till för "
        "sig. Utan aktiv karta läggs inget till."
    ),
    "pick_scans": (
        "Avmarkerat används den senaste skanningen för varje ruta, utan att katalogen "
        "behöver kontaktas medan dialogen är öppen. Markerat söker verktyget i Lantmäteriets "
        "öppna STAC-katalog efter alla skanningar som täcker området, och du väljer bland "
        "dem i listan nedanför. Sökningen görs igen när området ändras och tar ungefär en "
        "sekund. Den kräver ingen inloggning."
    ),
    "scans": (
        "Skanningarna som täcker området, en rad per skanningsområde med tid, lövläge, "
        "punkttäthet (alla returer, ur katalogens punktantal per ruta) och hur stor del av områdets utbredning den täcker. Många "
        "områden är skannade två gånger. En nyare skanning är ofta tätare, men kan vara "
        "gjord utan löv, och då ger lövträd och buskar få träffar: studien fann att "
        "punkttätheten är den viktigaste faktorn, men ett lövfritt buskskikt syns inte "
        "alls. Välj flera om området spänner över flera skanningsområden. Täcker två valda "
        "skanningar samma ruta används den senaste av dem. 'Med löv' betyder juni-augusti, "
        "'lövläge osäkert' maj, september och oktober, och 'utan löv' november-april."
    ),
    "min_points": (
        "Minsta antal punkter i nämnaren för att cellen ska få ett värde. Nämnaren är "
        "punkterna under skiktets övre gräns (NRD) eller alla punkter (ORD). Celler med "
        "färre blir NoData, eftersom ett förhållande mellan några få punkter mest är brus: "
        "med 20 punkter kan värdet bara ändras i steg om 0,05. Punktantalen kan sparas med "
        "rutan nedanför så att du kan välja en egen gräns efteråt."
    ),
    "save_counts": (
        "Spara nämnaren per cell som ett eget heltalsraster, så att du kan se hur säkert "
        "varje värde är. För NRD ett raster per unik övre gräns (punkter under gränsen), "
        "för ORD ett med alla punkter. Celler utan punkter har värdet 0."
    ),
    "consumer_key": (
        "Consumer key från Lantmäteriets API-portal. Nyckeln behöver både API:et STAC-hojd "
        "och en beställning av Laserdata Nedladdning, skog på Geotorget. Samma nyckel som "
        "till Laserdata Skog-verktyget."
    ),
    "consumer_secret": (
        "Consumer secret som hör till nyckeln. Visas dold i dialogen och skrivs aldrig till "
        "meddelandena."
    ),
    "memory_gb": (
        "Ungefär hur mycket arbetsminne verktyget får använda utöver ArcGIS Pro självt. "
        "Tomt = automatiskt: hälften av det lediga arbetsminnet när körningen startar, men "
        "minst 2 GB lämnas kvar åt Pro och Windows. Vilket värde som användes står i "
        "meddelandena. Området delas i block som var för sig ryms i minnet, så ytan "
        "begränsas inte av minnet. Mindre budget ger fler och mindre block, vilket blir något "
        "långsammare (ungefär 1 s extra per block). Varje process behöver ungefär 260 byte "
        "per läst punkt plus 150 MB. Minst 0,5 GB."
    ),
    "workers": (
        "Antal block som bearbetas samtidigt, vart och ett i en egen process. Tomt = "
        "automatiskt: en process per processortråd utom en, högst 16, och högst en per "
        "0,75 GB i minnesbudgeten. Aldrig fler än det finns block. Läsningen väntar mest på "
        "nätverket, så flera block samtidigt går fortare, men minnesbudgeten delas mellan "
        "dem och varje block blir mindre. Uppmätt: 4 samtidigt läste 2,5 gånger fortare än "
        "ett; över 4 är vinsten inte uppmätt mot Lantmäteriets server. 1 till 32."
    ),
    "out_rasters": "De skapade rastren, ett per skikt, plus punktantalen om de sparades.",
}


# =============================================================================
# Förlopp och format
# =============================================================================

def _fmt_duration(seconds):
    seconds = max(0, int(round(seconds)))
    if seconds < 60:
        return "{} s".format(seconds)
    if seconds < 3600:
        return "{} min".format(int(round(seconds / 60.0)))
    return "{} h {} min".format(seconds // 3600, (seconds % 3600) // 60)


def _fmt_count(n):
    return "{:,}".format(int(n)).replace(",", " ")


def _fmt_m(v):
    """Meter med decimalkomma och utan onödiga nollor: 0,15, 2, 1,85."""
    return "{:g}".format(round(float(v), 3)).replace(".", ",")


def _band_label(lo, hi):
    return "{}-{} m".format(_fmt_m(lo), _fmt_m(hi))


def _cm(v):
    """Gräns i hela cm för rasternamn: 0,15 -> 015, 10 -> 1000."""
    return "{:03d}".format(int(round(float(v) * 100)))


class _Steps:
    """Numrerade steg i förloppsindikatorn och i meddelandena."""

    def __init__(self, total, messages):
        self.total = total
        self.messages = messages
        self.k = 0
        self.t0 = time.time()
        self.t_step = self.t0

    def next(self, text):
        self.k += 1
        self.t_step = time.time()
        label = "Steg {} av {}: {}".format(self.k, self.total, text)
        arcpy.SetProgressor("default", label)
        self.messages.addMessage(label)

    def label(self, text):
        arcpy.SetProgressorLabel("Steg {} av {}: {}".format(self.k, self.total, text))

    def done(self, text=None):
        msg = "    klart på {}".format(_fmt_duration(time.time() - self.t_step))
        self.messages.addMessage(msg + (". " + text if text else "."))


# =============================================================================
# HTTP
# =============================================================================

def _http(req):
    """urlopen med omförsök vid tillfälliga fel. Returnerar (status, body)."""
    delay = 2
    for attempt in range(HTTP_RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as exc:
            # HTTPError först: den är också URLError och OSError.
            if exc.code in (408, 429) or exc.code >= 500:
                if attempt < HTTP_RETRIES - 1:
                    time.sleep(delay)
                    delay *= 2
                    continue
            raise
        except (urllib.error.URLError, OSError):
            if attempt < HTTP_RETRIES - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise


def _get_token(key, secret):
    auth = base64.b64encode("{}:{}".format(key, secret).encode("utf-8")).decode("ascii")
    req = urllib.request.Request(
        TOKEN_URL,
        data=urllib.parse.urlencode({"grant_type": "client_credentials"}).encode("ascii"),
        headers={"Authorization": "Basic " + auth, "User-Agent": USER_AGENT,
                 "Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    try:
        _status, body = _http(req)
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 401):
            raise ValueError(
                "Lantmäteriet godkände inte consumer key/secret (HTTP {}). Kontrollera "
                "nyckeln i Lantmäteriets API-portal.".format(exc.code)
            )
        raise
    return json.loads(body.decode("utf-8"))["access_token"]


def _check_access(url, token):
    """Hämta en byte av första filen, för ett begripligt fel i stället för PDAL:s."""
    req = urllib.request.Request(
        url, headers={"Authorization": "Bearer " + token, "Range": "bytes=0-0",
                      "User-Agent": USER_AGENT},
    )
    try:
        _http(req)
    except urllib.error.HTTPError as exc:
        if exc.code in (401, 403):
            raise ValueError(
                "Nyckeln saknar behörighet till punktmolnen (HTTP {}). Den behöver både "
                "API:et STAC-hojd i Lantmäteriets API-portal och en beställning av "
                "Laserdata Nedladdning, skog på Geotorget:\n{}".format(exc.code, GEOTORGET_URL)
            )
        raise


def _stac_search(bbox_wgs84):
    """Alla items i samlingen som skär bbox, över alla sidor."""
    body = {"collections": [COLLECTION], "bbox": list(bbox_wgs84), "limit": 100}
    url, method = STAC_SEARCH_URL, "POST"
    items = []
    seen_pages = set()
    while url:
        data = json.dumps(body).encode("utf-8") if method == "POST" else None
        req = urllib.request.Request(
            url, data=data, method=method,
            headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        )
        _status, raw = _http(req)
        page = json.loads(raw.decode("utf-8"))
        feats = page.get("features", [])
        # Skydd mot en pager som returnerar samma sida om och om igen.
        page_key = tuple(f.get("id") for f in feats)
        if page_key in seen_pages:
            break
        seen_pages.add(page_key)
        items.extend(feats)

        nxt = next((l for l in page.get("links", []) if l.get("rel") == "next"), None)
        if not nxt or not feats:
            break
        url = nxt["href"]
        method = nxt.get("method", "GET").upper()
        body = nxt.get("body", body)
    return items


class _Token:
    """Token som förnyas innan den går ut."""

    def __init__(self, key, secret):
        self._key, self._secret = key, secret
        self._token = _get_token(key, secret)
        self._time = time.time()

    def get(self):
        if time.time() - self._time > TOKEN_MAX_AGE_S:
            self._token = _get_token(self._key, self._secret)
            self._time = time.time()
        return self._token


# =============================================================================
# Geometri
# =============================================================================

def _tm():
    return arcpy.SpatialReference(SWEREF99TM_WKID)


def _rect_polygon(x0, y0, x1, y1):
    return arcpy.Polygon(arcpy.Array([arcpy.Point(x0, y0), arcpy.Point(x0, y1),
                                      arcpy.Point(x1, y1), arcpy.Point(x1, y0)]), _tm())


def _aoi_geometry(aoi_layer):
    """Alla (valda) polygoner i lagret, sammanslagna, i SWEREF 99 TM."""
    sr_in = arcpy.Describe(aoi_layer).spatialReference
    if sr_in is None or not (sr_in.factoryCode or sr_in.exportToString()):
        raise ValueError("Intresseområdet saknar koordinatsystem.")
    geom = None
    with arcpy.da.SearchCursor(aoi_layer, ["SHAPE@"], spatial_reference=_tm()) as cur:
        for (shape,) in cur:
            if shape is None or shape.area <= 0:
                continue
            geom = shape if geom is None else geom.union(shape)
    if geom is None:
        raise ValueError("Intresseområdet innehåller inga polygoner med yta.")
    return geom


def _active_map_sr():
    try:
        m = arcpy.mp.ArcGISProject("CURRENT").activeMap
        if m is not None and m.spatialReference is not None:
            return m.spatialReference, m.name
    except Exception:
        pass
    return None, None


def _densified_rect(xmin, ymin, xmax, ymax, sr):
    """Rektangel med förtätade kanter, omprojicerad till SWEREF 99 TM."""
    n = 16
    pts = ([(xmin + (xmax - xmin) * i / n, ymin) for i in range(n)]
           + [(xmax, ymin + (ymax - ymin) * i / n) for i in range(n)]
           + [(xmax - (xmax - xmin) * i / n, ymax) for i in range(n)]
           + [(xmin, ymax - (ymax - ymin) * i / n) for i in range(n)])
    poly = arcpy.Polygon(arcpy.Array([arcpy.Point(x, y) for x, y in pts]), sr)
    if sr.factoryCode != SWEREF99TM_WKID:
        poly = poly.projectAs(_tm())
    return poly


def _extent_polygon(value, text):
    """
    En GPExtent-parameter som polygon i SWEREF 99 TM. Returnerar (polygon, text om
    vilket koordinatsystem som användes).

    .value är ett geoprocessing-extentobjekt, inte arcpy.Extent. Hörnen är
    vanliga tal. Koordinatsystemet följer bara med när utbredningen kommer från
    ett lager eller en datakälla, och då bara som WKT2 i valueAsText efter de
    fyra talen. Inskrivna koordinater har inget; de tolkas i den aktiva kartans
    koordinatsystem, eftersom en utbredning vald i dialogen anges i kartans
    koordinater.
    """
    xmin, ymin, xmax, ymax = (float(value.XMin), float(value.YMin),
                              float(value.XMax), float(value.YMax))
    if not (xmax > xmin and ymax > ymin):
        raise ValueError("Utbredningen har ingen yta.")

    sr = None
    parts = (text or "").split(" ", 4)
    if len(parts) == 5 and parts[4].strip():
        sr = arcpy.SpatialReference()
        try:
            sr.loadFromString(parts[4].strip())
        except Exception:
            sr = None
    if sr is not None and (sr.factoryCode or sr.exportToString()):
        source = "utbredningens eget koordinatsystem"
    else:
        sr, map_name = _active_map_sr()
        if sr is not None and (sr.factoryCode or sr.exportToString()):
            source = "den aktiva kartans koordinatsystem ({})".format(map_name)
        else:
            sr = _tm()
            source = "ingen aktiv karta, så SWEREF 99 TM antas"
    poly = _densified_rect(xmin, ymin, xmax, ymax, sr)
    if poly is None or poly.area <= 0:
        raise ValueError("Utbredningen kunde inte omvandlas till SWEREF 99 TM.")
    return poly, "Utbredningen tolkas i {}: {}.".format(sr.name, source)


def _aoi_bbox_quick(p):
    """
    Områdets utbredning som polygon i SWEREF 99 TM, snabbt nog för
    updateParameters: lagrets utbredning i stället för dess polygoner.
    None om området inte är angivet eller inte kan läsas.
    """
    try:
        if p["aoi_mode"].valueAsText == AOI_EXTENT:
            if not p["aoi_extent"].valueAsText:
                return None
            return _extent_polygon(p["aoi_extent"].value, p["aoi_extent"].valueAsText)[0]
        if not p["aoi"].valueAsText:
            return None
        ext = arcpy.Describe(p["aoi"].value).extent
        sr = ext.spatialReference
        if sr is None or not (sr.factoryCode or sr.exportToString()):
            return None
        return _densified_rect(ext.XMin, ext.YMin, ext.XMax, ext.YMax, sr)
    except Exception:
        return None


def _grid(extent, cell):
    """Rutnät justerat till jämna multiplar av cellstorleken."""
    x0 = math.floor(extent.XMin / cell) * cell
    y0 = math.floor(extent.YMin / cell) * cell
    width = int(math.ceil((extent.XMax - x0) / cell))
    height = int(math.ceil((extent.YMax - y0) / cell))
    return {"resolution": cell, "origin_x": x0, "origin_y": y0,
            "width": width, "height": height}


# =============================================================================
# Skanningar i katalogen
# =============================================================================

def _actual_density(tile):
    """Punkter per m² med alla returer, ur pc:count och rutans yta."""
    x0, y0, x1, y1 = tile["bbox"]
    return tile["count"] / ((x1 - x0) * (y1 - y0)) if tile["count"] else None


def _tile_from_item(it):
    props = it.get("properties", {})
    asset = it.get("assets", {}).get("data")
    pb = props.get("proj:bbox") or (asset or {}).get("proj:bbox")
    if not asset or not pb:
        return None
    dt = props.get("datetime") or ""
    return {"id": it["id"], "href": asset["href"], "datetime": dt,
            "bbox": tuple(pb[:4]),
            "key": tuple(round(v) for v in pb[:4]),
            "count": props.get("pc:count") or 0,
            "start": props.get("start_datetime") or dt,
            "end": props.get("end_datetime") or dt,
            "area": props.get("skanningsomrade") or "",
            "flyghojd": props.get("flyghojd"),
            "punkttathet": props.get("punkttathet"),
            "modified": props.get("data_modified") or ""}


def _last_day(tile):
    """Sista flygdagen. end_datetime är midnatt efter sista dagen, så backa en dag då."""
    start, last = tile["start"][:10], tile["end"][:10]
    if last > start and tile["end"][11:19] == "00:00:00":
        try:
            last = (datetime.date.fromisoformat(last) - datetime.timedelta(days=1)).isoformat()
        except ValueError:
            pass
    return last or start


def _capture_period(tile):
    """'2021-03-07 - 2021-04-01', eller ett enda datum om start och slut är samma dag."""
    start = tile["start"][:10]
    if not start:
        return "okänt"
    last = _last_day(tile)
    return start if last in ("", start) else "{} - {}".format(start, last)


def _months(tiles):
    """Alla månader (år, månad) som rutornas insamling spänner över."""
    out = set()
    for t in tiles:
        try:
            d0 = datetime.date.fromisoformat(t["start"][:10])
            d1 = datetime.date.fromisoformat(_last_day(t))
        except ValueError:
            continue
        y, m = d0.year, d0.month
        while (y, m) <= (d1.year, d1.month):
            out.add((y, m))
            y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return sorted(out)


def _season(months):
    """Lövläge ur insamlingsmånaderna: 'med löv', 'utan löv' eller 'lövläge osäkert'."""
    kinds = set()
    for _y, m in months:
        kinds.add("med" if 6 <= m <= 8 else "utan" if m in (11, 12, 1, 2, 3, 4) else "osäkert")
    if kinds == {"med"}:
        return "med löv"
    if kinds == {"utan"}:
        return "utan löv"
    return "lövläge osäkert"


def _period_text(months):
    """'april 2026', 'mars-april 2026' eller 'december 2019 - januari 2020'."""
    if not months:
        return "okänd tid"
    (y0, m0), (y1, m1) = months[0], months[-1]
    if (y0, m0) == (y1, m1):
        return "{} {}".format(MONTHS[m0 - 1], y0)
    if y0 == y1:
        return "{}-{} {}".format(MONTHS[m0 - 1], MONTHS[m1 - 1], y0)
    return "{} {} - {} {}".format(MONTHS[m0 - 1], y0, MONTHS[m1 - 1], y1)


def _overlaps(aoi, bbox):
    """
    Sant om rutan och aoi har gemensam yta. disjoint() räcker inte: en ruta som
    bara nuddar områdets kant skulle då räknas med och visas med 0 % täckning.
    """
    rect = _rect_polygon(*bbox)
    return not aoi.disjoint(rect) and aoi.intersect(rect, 4).area > 1.0


def _coverage(tiles, aoi):
    """Andel av aoi:s yta som täcks av rutorna (0-1)."""
    if aoi is None or aoi.area <= 0 or not tiles:
        return 0.0
    union = None
    for t in tiles:
        r = _rect_polygon(*t["bbox"])
        union = r if union is None else union.union(r)
    return min(1.0, aoi.intersect(union, 4).area / aoi.area)


def _scan_catalog(items, aoi):
    """
    Skanningarna som täcker aoi, en per skanningsområde, senaste först. Varje
    post har code, tiles, months, season, density, coverage och label.
    """
    groups = {}
    for it in items:
        t = _tile_from_item(it)
        if t is None or not _overlaps(aoi, t["bbox"]):
            continue
        groups.setdefault(t["area"] or t["id"].split("-")[0].upper(), []).append(t)
    scans = []
    for code, tiles in groups.items():
        months = _months(tiles)
        dens = [d for d in (_actual_density(t) for t in tiles) if d]
        density = sum(dens) / len(dens) if dens else None
        cov = _coverage(tiles, aoi)
        label = "{} - {}, {}, {}täcker {} % av utbredningen".format(
            code, _period_text(months), _season(months),
            "{} p/m², ".format(_fmt_m(round(density, 1))) if density else "",
            int(round(100 * cov)))
        scans.append({"code": code, "tiles": tiles, "months": months,
                      "season": _season(months), "density": density, "coverage": cov,
                      "label": label, "latest": max(t["datetime"] for t in tiles)})
    scans.sort(key=lambda s: s["latest"], reverse=True)
    return scans


def _scan_code(label):
    """Skanningsområdets kod ur en rad i listan (eller en kod som redan är ren)."""
    return (label or "").split(" - ")[0].strip().upper()


def _pick_tiles(items, aoi, codes=None):
    """
    Rutor som skär själva polygonen, en insamling per ruta: den senaste, bland
    de valda skanningsområdena om codes är angivet.
    """
    newest = {}
    for it in items:
        t = _tile_from_item(it)
        if t is None:
            continue
        if codes is not None and (t["area"] or "").upper() not in codes:
            continue
        if not _overlaps(aoi, t["bbox"]):
            continue
        k = t["key"]
        if k not in newest or t["datetime"] > newest[k]["datetime"]:
            newest[k] = t
    return sorted(newest.values(), key=lambda t: t["id"])


# Katalogsökningar i dialogen, per områdesutbredning. Modulnivå, eftersom Pro
# inte alltid återanvänder samma verktygsobjekt mellan valideringsanropen.
_SCAN_CACHE = {}


def _scans_for_dialog(aoi):
    """(skanningar, felmeddelande) för aoi, sökt en gång per utbredning."""
    e = aoi.extent
    key = (round(e.XMin), round(e.YMin), round(e.XMax), round(e.YMax))
    if key not in _SCAN_CACHE:
        try:
            wgs = aoi.projectAs(arcpy.SpatialReference(4326)).extent
            items = _stac_search((wgs.XMin, wgs.YMin, wgs.XMax, wgs.YMax))
            _SCAN_CACHE[key] = (_scan_catalog(items, aoi), None)
        except Exception as exc:
            # Inte sparat vid fel, så att nästa validering försöker igen.
            return [], "Kunde inte söka i Lantmäteriets STAC-katalog: {}".format(exc)
    return _SCAN_CACHE[key]


# =============================================================================
# Block
# =============================================================================

def _available_memory_gb():
    """Ledigt arbetsminne i GB, eller None. psutil ingår i arcgispro-py3."""
    try:
        import psutil
        return psutil.virtual_memory().available / 1e9
    except Exception:
        return None


def _cpu_threads():
    """Processortrådar som processen får använda (logiska kärnor)."""
    n = getattr(os, "process_cpu_count", os.cpu_count)()
    return max(1, n or 1)


def _resolve_resources(memory_gb, workers):
    """
    Minnesbudget och antal processer, där None betyder automatiskt. Returnerar
    (minne i GB, processer, text om vad som valdes automatiskt eller "").
    """
    notes = []
    if memory_gb is None:
        avail = _available_memory_gb()
        if avail is None:
            memory_gb = 4.0
            notes.append("minnesbudget 4 GB (ledigt minne kunde inte läsas)")
        else:
            memory_gb = max(MIN_MEMORY_GB, min(AUTO_MEMORY_SHARE * avail,
                                               avail - AUTO_MEMORY_RESERVE_GB))
            notes.append("minnesbudget {} GB av {} GB ledigt".format(
                _fmt_m(round(memory_gb, 1)), _fmt_m(round(avail, 1))))
    memory_gb = max(MIN_MEMORY_GB, float(memory_gb))
    if workers is None:
        threads = _cpu_threads()
        workers = max(1, min(threads - 1, MAX_WORKERS, int(memory_gb / MIN_GB_PER_WORKER)))
        notes.append("{} processer ({} processortrådar)".format(workers, threads))
    workers = max(1, min(32, int(workers)))
    note = "Automatiskt: {}.".format(", ".join(notes)) if notes else ""
    # Även ett angivet antal måste rymmas: varje process behöver sitt fasta minne
    # och ett block av rimlig storlek.
    fit = max(1, int(memory_gb / MIN_GB_PER_WORKER))
    if workers > fit:
        note = (note + " " if note else "") + (
            "{} processer ryms inte i {} GB; använder {} (minst {} GB per process).".format(
                workers, _fmt_m(round(memory_gb, 1)), fit, _fmt_m(MIN_GB_PER_WORKER)))
        workers = fit
    return memory_gb, workers, note


def _plan_blocks(aoi, grid, tiles, memory_gb, workers):
    """
    Dela områdets rutnät i block som vart och ett ryms i minnesbudgeten delat
    med antalet samtidiga block. Blockstorleken räknas från den tätaste rutans
    punkttäthet (pc:count / rutans yta) och är en hel multipel av cellstorleken,
    så att varje cell ligger i exakt ett block. Blocken görs också så små att de
    blir minst lika många som processerna (men inte under MIN_BLOCK_M). Block
    som inte skär polygonen hoppas över.

    Returnerar (block, blocksida i m). Varje block har sina hörn, sitt delrutnät,
    läsområdet med marginal och de rutor det behöver läsa.
    """
    cell = grid["resolution"]
    margin = READ_MARGIN_M
    density = max(t["count"] / ((t["bbox"][2] - t["bbox"][0]) * (t["bbox"][3] - t["bbox"][1]))
                  for t in tiles) or 1.0
    per_worker = memory_gb * 1e9 / max(1, workers) - WORKER_BASE_BYTES
    points_per_block = max(1e5, per_worker / PEAK_BYTES_PER_POINT)
    side = math.sqrt(points_per_block / density) - 2 * margin
    # Minst lika många block som processer, annars står processer oanvända:
    # läsningen väntar mest på nätverket, så parallellitet väger tyngre än den
    # extra sekunden det kostar att öppna fjärrfilen för varje block.
    area = grid["width"] * grid["height"] * cell * cell
    side = min(side, math.sqrt(area / max(1, workers)))
    side = min(MAX_BLOCK_M, max(MIN_BLOCK_M, side))
    side_cells = max(1, int(side // cell))

    blocks = []
    n_cols = int(math.ceil(grid["width"] / side_cells))
    n_rows = int(math.ceil(grid["height"] / side_cells))
    for ry in range(n_rows):
        for cx in range(n_cols):
            c0, r0 = cx * side_cells, ry * side_cells
            w = min(side_cells, grid["width"] - c0)
            h = min(side_cells, grid["height"] - r0)
            x0 = grid["origin_x"] + c0 * cell
            y0 = grid["origin_y"] + r0 * cell
            x1, y1 = x0 + w * cell, y0 + h * cell
            if aoi.disjoint(_rect_polygon(x0, y0, x1, y1)):
                continue
            mx0, my0 = x0 - margin, y0 - margin
            mx1, my1 = x1 + margin, y1 + margin
            reads = []
            expected = 0.0
            for t in tiles:
                tx0, ty0, tx1, ty1 = t["bbox"]
                bx0, bx1 = max(tx0, mx0), min(tx1, mx1)
                by0, by1 = max(ty0, my0), min(ty1, my1)
                if bx1 <= bx0 or by1 <= by0:
                    continue
                bounds = "([{:.2f},{:.2f}],[{:.2f},{:.2f}])".format(bx0, bx1, by0, by1)
                reads.append((t, bounds))
                expected += t["count"] * (bx1 - bx0) * (by1 - by0) / ((tx1 - tx0) * (ty1 - ty0))
            if not reads:
                continue
            blocks.append({
                # Rad räknas från norr, som man läser en karta.
                "row": n_rows - ry, "col": cx + 1,
                "c0": c0, "r0": r0,
                "x0": x0, "y0": y0, "x1": x1, "y1": y1,
                "read_bbox": [mx0, my0, mx1, my1],
                "grid": {"resolution": cell, "origin_x": x0, "origin_y": y0,
                         "width": w, "height": h},
                "reads": reads, "expected": expected,
            })
    return blocks, side_cells * cell


# =============================================================================
# Arbetsprocesser
# =============================================================================

def _python_exe():
    """python.exe i den aktiva miljön. I Pro är sys.executable ArcGISPro.exe."""
    for cand in (os.path.join(sys.exec_prefix, "python.exe"), sys.executable):
        if cand and os.path.basename(cand).lower() == "python.exe" and os.path.isfile(cand):
            return cand
    raise ValueError("Hittar inte python.exe i den aktiva miljön ({}).".format(sys.exec_prefix))


class _WorkerPool:
    """
    Ett fast antal nrd_worker-processer som tar ett block i taget.

    Varför processer och inte trådar: flera PDAL-pipelines i trådar i samma
    process som arcpy kraschade processen (access violation i arcpy, uppmätt i
    Laserdata Skog). Processerna startas en gång (cirka 1,7 s för Python, numpy
    och PDAL) och återanvänds.

    En läsartråd per process väntar på svar och lägger dem i en kö; trådarna
    gör bara I/O mot rören, aldrig PDAL eller arcpy.
    """

    def __init__(self, n, workdir):
        self._results = queue.Queue()
        self._procs = []
        self._logs = []
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        exe = _python_exe()
        script = os.path.join(_HERE, "nrd_worker.py")
        env = dict(os.environ, PYTHONIOENCODING="utf-8")
        for i in range(n):
            log_path = os.path.join(workdir, "worker_{}.log".format(i))
            log = open(log_path, "w", encoding="utf-8", errors="replace")
            proc = subprocess.Popen(
                [exe, "-u", script], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=log, encoding="utf-8", errors="replace", creationflags=flags,
                env=env, cwd=workdir)
            self._procs.append(proc)
            self._logs.append((log, log_path))
            threading.Thread(target=self._reader, args=(i, proc), daemon=True).start()

    def _reader(self, i, proc):
        for line in proc.stdout:
            if line.startswith(_nw.RESULT_PREFIX):
                self._results.put((i, json.loads(line[len(_nw.RESULT_PREFIX):])))
        self._results.put((i, None))  # processen har avslutats

    def submit(self, i, task):
        self._procs[i].stdin.write(json.dumps(task) + "\n")
        self._procs[i].stdin.flush()

    def get(self, timeout):
        """(processindex, svar) eller None vid timeout. Svar None = processen dog."""
        try:
            return self._results.get(timeout=timeout)
        except queue.Empty:
            return None

    def log_tail(self, i, n=15):
        log, path = self._logs[i]
        log.flush()
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                return "".join(fh.readlines()[-n:])
        except OSError:
            return ""

    def close(self, kill=False):
        for proc in self._procs:
            try:
                if kill:
                    proc.kill()
                else:
                    proc.stdin.close()
            except OSError:
                pass
        for proc in self._procs:
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
        for log, _path in self._logs:
            log.close()


def _block_task(blk, token, bands, workdir):
    """Blocket som JSON-uppgift till en arbetsprocess (bara det den behöver)."""
    return {
        "block": {k: blk[k] for k in ("row", "col", "x0", "y0", "x1", "y1", "grid", "read_bbox")}
        | {"reads": [{"href": t["href"], "id": t["id"], "bounds": bounds}
                     for t, bounds in blk["reads"]]},
        "token": token, "bands": bands, "workdir": workdir,
    }


# =============================================================================
# Utdata och metadata
# =============================================================================

def _cell_tag(cell):
    """Cellstorlek i rasternamn: 20 -> 20m, 2,5 -> 2p5m (punkt är inte tillåten i namn)."""
    return "{:g}m".format(round(float(cell), 3)).replace(".", "p")


def _out_name(prefix, metric, cell, kind="ratio", lo=None, hi=None):
    """
    Rasternamn med prefix, mått, cellstorlek och höjdgränser i cm:
        <prefix>_nrd_20m_015_185         kvoten för skiktet 0,15-1,85 m
        <prefix>_nrd_20m_antal_under_185 nämnaren i NRD (punkter under 1,85 m)
        <prefix>_ord_20m_antal_alla      nämnaren i ORD (alla punkter)
    """
    base = "{}_{}_{}".format(prefix, metric.lower(), _cell_tag(cell))
    if kind == "under":
        return "{}_antal_under_{}".format(base, _cm(hi))
    if kind == "alla":
        return "{}_antal_alla".format(base)
    return "{}_{}_{}".format(base, _cm(lo), _cm(hi))


def _out_path(workspace, name):
    is_gdb = str(workspace).lower().endswith(".gdb")
    return os.path.join(workspace, name if is_gdb else name + ".tif")


def _planned_names(prefix, metric, cell, bands, save_counts):
    names = [_out_name(prefix, metric, cell, lo=lo, hi=hi) for lo, hi in bands]
    if save_counts:
        if metric == "NRD":
            for hi in sorted({hi for _lo, hi in bands}):
                names.append(_out_name(prefix, metric, cell, "under", hi=hi))
        else:
            names.append(_out_name(prefix, metric, cell, "alla"))
    return names


def _group_name(prefix, metric, cell):
    return "{} - {} {} m".format(prefix, metric, _fmt_m(cell))


def _default_workspace():
    try:
        gdb = arcpy.mp.ArcGISProject("CURRENT").defaultGeodatabase
        if gdb:
            return gdb
    except Exception:
        pass
    return None


def _norm_path(path):
    return os.path.normcase(os.path.normpath(str(path)))


def _release_from_maps(paths, messages, aprx=None):
    """
    Ta bort lager i projektets kartor som visar raster som ska skrivas över, och
    de grupper de låg i om grupperna blir tomma av det. Ett lager håller ett schemalås på rastret, och
    att skriva över det misslyckas då med ERROR 000464 (uppmätt vid en andra
    körning med samma prefix när den första körningens raster låg i kartan).
    """
    targets = {_norm_path(p) for p in paths}
    try:
        aprx = aprx or arcpy.mp.ArcGISProject("CURRENT")
    except Exception:
        return
    removed = 0
    for m in aprx.listMaps():
        parents = set()
        for lyr in m.listLayers():
            try:
                if (lyr.supports("DATASOURCE") and not lyr.isGroupLayer
                        and _norm_path(lyr.dataSource) in targets):
                    # longName är "grupp\lager" för ett lager i en grupp.
                    if "\\" in lyr.longName:
                        parents.add(lyr.longName.rsplit("\\", 1)[0])
                    m.removeLayer(lyr)
                    removed += 1
            except Exception:
                continue
        for lyr in m.listLayers():
            try:
                if lyr.isGroupLayer and lyr.longName in parents and not lyr.listLayers():
                    m.removeLayer(lyr)
            except Exception:
                continue
    if removed:
        messages.addMessage("    Tog bort {} lager ur kartan som visade raster som skrivs "
                            "över.".format(removed))


def _add_to_map(paths, group_name, messages, aprx=None):
    """
    Lägg rastren i den aktiva kartan, i en ny grupp om group_name är angivet.
    Rastren hamnar i gruppen i samma ordning som paths, uppifrån och ned.
    aprx kan anges vid test utanför Pro.
    """
    if aprx is None:
        try:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
        except Exception:
            return None
        m = aprx.activeMap
    else:
        # Ett projekt öppnat från fil har ingen aktiv karta; ta den första.
        m = aprx.activeMap or next(iter(aprx.listMaps()), None)
    if m is None:
        messages.addWarningMessage("Ingen aktiv karta - rastren läggs inte till.")
        return None
    if not group_name:
        for p in paths:
            m.addDataFromPath(p)
        return m
    group = m.createGroupLayer(group_name)
    # addLayerToGroup lägger lagret överst i gruppen, så lägg till bakifrån.
    for p in reversed(paths):
        lyr = m.addDataFromPath(p)
        m.addLayerToGroup(group, lyr, "TOP")
        m.removeLayer(lyr)
    messages.addMessage("    Gruppen '{}' i kartan {} har {} raster.".format(
        group_name, m.name, len(paths)))
    return m


def _write_raster_metadata(path, info, tiles, run):
    """
    Titel, sammanfattning och en beskrivning med alla inställningar, källrutor
    med insamlingsdatum, villkor och taggar.
    info: {"kind": "ratio"|"count", "lo", "hi", ...}; run: körningens inställningar.
    """
    h = escape
    metric = run["metric"]
    periods = sorted({_capture_period(t) for t in tiles})
    months = _months(tiles)
    band = _band_label(info["lo"], info["hi"]) if info.get("lo") is not None else None

    if info["kind"] == "ratio":
        if metric == "NRD":
            formula = ("NRD = antal punkter med {lo} m &lt;= h &lt; {hi} m / antal punkter med "
                       "h &lt; {hi} m").format(lo=_fmt_m(info["lo"]), hi=_fmt_m(info["hi"]))
            meaning = ("Andelen av de punkter som nått ned till skiktet eller längre som "
                       "träffade något i skiktet. Påverkas mindre av ett tätt krontak ovanför "
                       "än ORD.")
        else:
            formula = ("ORD = antal punkter med {lo} m &lt;= h &lt; {hi} m / antal punkter "
                       "totalt").format(lo=_fmt_m(info["lo"]), hi=_fmt_m(info["hi"]))
            meaning = ("Andelen av alla punkter i cellen som ligger i skiktet. Blir låg under "
                       "ett tätt krontak även där skiktet är tätt.")
        title = "{} {} ({:g} m celler)".format(metric, band, run["cell"])
        what = ("<p><b>Värde:</b> {} i höjdskiktet {} över mark, 0-1. {}</p>"
                "<p><b>Formel:</b> {}, där h är punktens höjd över mark.</p>").format(
            metric, h(band), meaning, formula)
        stat = info.get("stats")
    else:
        if metric == "NRD":
            title = "Punktantal under {} m ({:g} m celler)".format(_fmt_m(info["hi"]), run["cell"])
            what = ("<p><b>Värde:</b> antal punkter per cell med höjd över mark under {} m, "
                    "alltså nämnaren i NRD för skikten med den övre gränsen.</p>").format(
                _fmt_m(info["hi"]))
        else:
            title = "Punktantal, alla punkter ({:g} m celler)".format(run["cell"])
            what = "<p><b>Värde:</b> antal punkter per cell, nämnaren i ORD.</p>"
        stat = None

    settings = [
        ("Mått", metric),
        ("Höjdskikt i körningen", ", ".join(_band_label(lo, hi) for lo, hi in run["bands"])),
        ("Cellstorlek", "{} m".format(_fmt_m(run["cell"]))),
        ("Minsta antal punkter i nämnaren", "{} (celler med färre är NoData)".format(
            run["min_points"])),
        ("Val av skanning", run["scan_mode"]),
        ("Punkter som används", "Alla klasser utom brus (7 och 18), alla returer"),
        ("Höjd över mark", "Markpunkterna (klass 2) triangulerade (Delaunay), rastrerade "
                           "till {} m och bilinjärt interpolerade i varje punkts läge".format(
                               _fmt_m(_nw.DTM_RES))),
        ("Skiktgränser", "Undre gränsen ingår, övre inte. Punkter under markytan räknas "
                         "som punkter under skiktet"),
        ("Punkter lästa inom utbredningen", _fmt_count(run["points"])),
        ("varav markpunkter", _fmt_count(run["ground"])),
        ("Punkter utan markhöjd (ej räknade)", _fmt_count(run["no_height"])),
        ("Skapad", "{} med verktyget {} {} ({})".format(
            run["created"], h(run["tool_label"]), TOOL_VERSION, REPO_NAME)),
    ]
    if stat:
        settings.insert(4, ("Celler med värde / NoData p.g.a. för få punkter",
                            "{} / {} ({} %)".format(_fmt_count(stat["valid"]),
                                                    _fmt_count(stat["few"]),
                                                    _fmt_m(stat["few_pct"]))))
        if stat.get("mean") is not None:
            settings.insert(5, ("Medelvärde och median", "{} och {}".format(
                _fmt_m(stat["mean"]), _fmt_m(stat["median"]))))
    rows = "".join("<tr><td>{}</td><td>{}</td></tr>".format(h(k), h(str(v))) for k, v in settings)
    tile_rows = "".join(
        "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>"
        "</tr>".format(
            h(t["id"]), h(t["area"]), h(_capture_period(t)), h(_season(_months([t]))),
            "{} m".format(t["flyghojd"]) if t["flyghojd"] else "",
            "{} p/m²".format(_fmt_m(t["punkttathet"])) if t["punkttathet"] else "",
            _fmt_m(round(_actual_density(t), 1)) if _actual_density(t) else "")
        for t in tiles)
    ext = run["extent"]
    desc = (
        "{what}"
        "<p><b>Inställningar</b></p><table border='1' cellpadding='3'>{rows}</table>"
        "<p><b>Källa:</b> Lantmäteriet, Laserdata Nedladdning, skog (STAC-samling {coll}, "
        "{n} ruta/rutor). Flygburen laserskanning, klassificerat punktmoln.</p>"
        "<p><b>Insamlingsdatum:</b> {periods}. Utan löv ger lövträd och buskar få träffar, så "
        "täthet i låga skikt underskattas i lövskog.</p>"
        "<table border='1' cellpadding='3'><tr><th>Ruta</th><th>Skanningsområde</th>"
        "<th>Insamlad</th><th>Lövläge</th><th>Flyghöjd</th><th>Punkttäthet (nominell)</th>"
        "<th>Punkter per m² (alla returer)</th></tr>{tile_rows}</table>"
        "<p><b>Metod:</b> Campbell, M.J., Dennison, P.E., Hudak, A.T., Parham, L.M., Butler, "
        "B.W. (2018). Quantifying understory vegetation density using small-footprint "
        "airborne lidar. Remote Sensing of Environment 215, 330-342. {doi}</p>"
        "<p><b>Koordinatsystem:</b> SWEREF 99 TM (EPSG:3006). Utbredning X {x0:.0f} - "
        "{x1:.0f}, Y {y0:.0f} - {y1:.0f}. Klippt till intresseområdet.</p>"
    ).format(what=what, rows=rows, coll=COLLECTION, n=len(tiles),
             periods=h(", ".join(periods)), tile_rows=tile_rows, doi=PAPER_DOI,
             x0=ext.XMin, x1=ext.XMax, y0=ext.YMin, y1=ext.YMax)

    md = arcpy.metadata.Metadata(path)
    md.title = "{} från Laserdata Skog, {}".format(title, _period_text(months))
    md.summary = "{} ur Lantmäteriets Laserdata Nedladdning, skog, insamlad {} ({}).".format(
        title, ", ".join(periods), _season(months))
    md.description = desc
    md.tags = "Lantmäteriet, Laserdata Skog, laserskanning, punkttäthet, {}, {}".format(
        metric, "buskskikt, undervegetation" if info["kind"] == "ratio" else "punktantal")
    md.credits = "© Lantmäteriet, Laserdata Nedladdning, skog."
    md.accessConstraints = (
        "Användningsvillkor för Laserdata Nedladdning, skog: {}".format(TERMS_URL))
    md.save()


def _write_tool_metadata(tool_cls, toolbox_alias):
    """
    Skriv verktygets metadatafil med parameterförklaringar från TOOLTIPS.

    Pro läser verktygstipsen i dialogen från <verktygslåda>.<verktyg>.pyt.xml
    (elementet dialogReference per parameter). Det finns inget attribut på
    arcpy.Parameter för detta. Filen skrivs bara om innehållet har ändrats.
    """
    toolbox = os.path.splitext(os.path.basename(__file__))[0]
    path = os.path.join(_HERE, "{}.{}.pyt.xml".format(toolbox, tool_cls.__name__))

    def html(text):
        body = escape(text).replace("\n", "</SPAN></P><P><SPAN>")
        return escape('<DIV STYLE="text-align:Left;"><P><SPAN>{}</SPAN></P></DIV>'.format(body))

    tool = tool_cls()
    params = []
    for p in tool.getParameterInfo():
        tip = TOOLTIPS.get(p.name)
        if not tip:
            continue
        params.append(
            '<param name="{n}" displayname="{d}" type="{t}" direction="{r}">'
            "<dialogReference>{h}</dialogReference>"
            "<pythonReference>{h}</pythonReference></param>".format(
                n=p.name, d=escape(p.displayName, {'"': "&quot;"}),
                t=p.parameterType, r=p.direction, h=html(tip))
        )
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<metadata xml:lang="sv"><Esri><ArcGISFormat>1.0</ArcGISFormat></Esri>'
        '<tool name="{name}" displayname="{label}" toolboxalias="{alias}" xmlns="">'
        "<parameters>{params}</parameters><summary>{summary}</summary></tool>"
        "<dataIdInfo><idCitation><resTitle>{label}</resTitle></idCitation>"
        "<idAbs>{summary}</idAbs></dataIdInfo></metadata>\n"
    ).format(name=tool_cls.__name__, label=escape(tool.label), alias=toolbox_alias,
             params="".join(params), summary=html(TOOL_SUMMARY))

    try:
        with open(path, encoding="utf-8") as fh:
            if fh.read() == xml:
                return
    except OSError:
        pass
    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(xml)
    except OSError:
        # Skrivskyddad plats: verktyget fungerar ändå, bara utan verktygstips.
        pass


# =============================================================================
# Parametrar
# =============================================================================

def _parse_bands(value):
    """
    Skikten ur värdetabellen som [[undre, övre], ...] i meter. Raises ValueError
    med en läsbar text vid fel.
    """
    rows = value or []
    bands = []
    for i, row in enumerate(rows, 1):
        if row is None or len(row) < 2 or row[0] in (None, "") or row[1] in (None, ""):
            raise ValueError("Skikt {}: ange både undre och övre gräns.".format(i))
        try:
            lo, hi = float(str(row[0]).replace(",", ".")), float(str(row[1]).replace(",", "."))
        except ValueError:
            raise ValueError("Skikt {}: gränserna måste vara tal.".format(i))
        if lo < 0:
            raise ValueError("Skikt {}: undre gränsen kan inte vara negativ.".format(i))
        if hi <= lo:
            raise ValueError("Skikt {}: övre gränsen måste vara större än den undre.".format(i))
        if hi > MAX_BAND_TOP_M:
            raise ValueError("Skikt {}: övre gränsen får vara högst {:g} m.".format(
                i, MAX_BAND_TOP_M))
        bands.append([lo, hi])
    if not bands:
        raise ValueError("Ange minst ett höjdskikt.")
    if len(bands) > MAX_BANDS:
        raise ValueError("Högst {} skikt.".format(MAX_BANDS))
    names = [(_cm(lo), _cm(hi)) for lo, hi in bands]
    dup = sorted({"{}-{}".format(*n) for n in names if names.count(n) > 1})
    if dup:
        raise ValueError("Samma skikt finns flera gånger (gränser avrundade till hela cm): "
                         "{}.".format(", ".join(dup)))
    return bands


# =============================================================================
# Toolbox
# =============================================================================

class Toolbox:
    def __init__(self):
        self.label = "Lantmäteriet Laserdata NRD"
        self.alias = "laserdata_nrd"
        self.tools = [PunkttathetISkikt]
        _write_tool_metadata(PunkttathetISkikt, self.alias)


class PunkttathetISkikt:
    def __init__(self):
        self.label = "Punkttäthet i höjdskikt (NRD/ORD)"
        self.description = TOOL_SUMMARY + (
            " Kräver en consumer key och secret från Lantmäteriets API-portal med "
            "behörighet till STAC-hojd och Laserdata Nedladdning, skog."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_mode = arcpy.Parameter(
            displayName="Avgränsa området med", name="aoi_mode", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_mode.filter.type = "ValueList"
        p_mode.filter.list = [AOI_POLYGONS, AOI_EXTENT]
        p_mode.value = AOI_POLYGONS

        # Båda är Optional i ramverket; updateMessages kräver den som valts.
        p_aoi = arcpy.Parameter(
            displayName="Intresseområde (polygoner)", name="aoi", datatype="GPFeatureLayer",
            parameterType="Optional", direction="Input",
        )
        p_aoi.filter.list = ["Polygon"]

        p_extent = arcpy.Parameter(
            displayName="Utbredning", name="aoi_extent", datatype="GPExtent",
            parameterType="Optional", direction="Input",
        )
        p_extent.enabled = False

        p_metric = arcpy.Parameter(
            displayName="Mått", name="metric", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_metric.filter.type = "ValueList"
        p_metric.filter.list = [METRIC_NRD, METRIC_ORD]
        p_metric.value = METRIC_NRD

        p_bands = arcpy.Parameter(
            displayName="Höjdskikt (m över mark)", name="bands", datatype="GPValueTable",
            parameterType="Required", direction="Input",
        )
        p_bands.columns = [["GPDouble", "Undre gräns (m)"], ["GPDouble", "Övre gräns (m)"]]
        p_bands.value = [list(b) for b in DEFAULT_BANDS]

        p_cell = arcpy.Parameter(
            displayName="Cellstorlek (m)", name="cell_size", datatype="GPDouble",
            parameterType="Required", direction="Input",
        )
        p_cell.value = DEFAULT_CELL_SIZE

        p_ws = arcpy.Parameter(
            displayName="Utdata-arbetsyta", name="out_workspace", datatype="DEWorkspace",
            parameterType="Required", direction="Input",
        )
        default_ws = _default_workspace()
        if default_ws:
            p_ws.value = default_ws

        p_prefix = arcpy.Parameter(
            displayName="Namnprefix", name="prefix", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_prefix.value = DEFAULT_PREFIX

        p_group = arcpy.Parameter(
            displayName="Samla resultatet i en grupp i kartan", name="add_group",
            datatype="GPBoolean", parameterType="Optional", direction="Input",
        )
        p_group.value = True

        p_pick = arcpy.Parameter(
            displayName="Välj skanning själv", name="pick_scans", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category="Skanning",
        )
        p_pick.value = False

        p_scans = arcpy.Parameter(
            displayName="Skanningar", name="scans", datatype="GPString",
            parameterType="Optional", direction="Input", category="Skanning",
            multiValue=True,
        )
        p_scans.filter.type = "ValueList"
        p_scans.filter.list = []
        p_scans.enabled = False

        p_min = arcpy.Parameter(
            displayName="Minsta antal punkter i nämnaren", name="min_points", datatype="GPLong",
            parameterType="Optional", direction="Input", category="Tillförlitlighet",
        )
        p_min.value = DEFAULT_MIN_POINTS

        p_counts = arcpy.Parameter(
            displayName="Spara punktantal (nämnaren)", name="save_counts", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category="Tillförlitlighet",
        )
        p_counts.value = True

        p_key = arcpy.Parameter(
            displayName="Consumer key", name="consumer_key", datatype="GPString",
            parameterType="Required", direction="Input", category="Inloggning",
        )
        p_secret = arcpy.Parameter(
            displayName="Consumer secret", name="consumer_secret", datatype="GPStringHidden",
            parameterType="Required", direction="Input", category="Inloggning",
        )

        p_mem = arcpy.Parameter(
            displayName="Minnesbudget (GB)", name="memory_gb", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Avancerat",
        )
        # Tomt = automatiskt, se _resolve_resources.

        p_workers = arcpy.Parameter(
            displayName="Block samtidigt", name="workers", datatype="GPLong",
            parameterType="Optional", direction="Input", category="Avancerat",
        )
        p_workers.filter.type = "Range"
        p_workers.filter.list = [1, 32]

        p_out = arcpy.Parameter(
            displayName="Utdata-raster", name="out_rasters", datatype="DERasterDataset",
            parameterType="Derived", direction="Output", multiValue=True,
        )

        return [p_mode, p_aoi, p_extent, p_metric, p_bands, p_cell, p_ws, p_prefix, p_group,
                p_pick, p_scans, p_min, p_counts, p_key, p_secret, p_mem, p_workers, p_out]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        p = {q.name: q for q in parameters}
        by_extent = p["aoi_mode"].valueAsText == AOI_EXTENT
        p["aoi"].enabled = not by_extent
        p["aoi_extent"].enabled = by_extent

        pick = bool(p["pick_scans"].value)
        p["scans"].enabled = pick
        if not pick:
            return
        # Nätverk bara när användaren valt att välja skanning själv, och bara en
        # gång per områdesutbredning (_SCAN_CACHE).
        aoi = _aoi_bbox_quick(p)
        if aoi is None:
            return
        scans, _err = _scans_for_dialog(aoi)
        labels = [s["label"] for s in scans]
        # Från ett skript kan listan och värdet komma i samma anrop, och en kod
        # utan resten av raden ska matcha sin rad. Okända värden behålls, så att
        # updateMessages kan säga vilka skanningar som finns i stället.
        chosen = [v for v in (p["scans"].values or []) if v]
        by_code = {s["code"]: s["label"] for s in scans}
        kept = []
        for v in chosen:
            label = v if v in labels else by_code.get(_scan_code(v), v)
            if label not in kept:
                kept.append(label)
        p["scans"].filter.list = labels
        if kept != chosen:
            p["scans"].values = kept or None

    def updateMessages(self, parameters):
        p = {q.name: q for q in parameters}

        if p["aoi_mode"].valueAsText == AOI_EXTENT:
            if not p["aoi_extent"].valueAsText:
                p["aoi_extent"].setErrorMessage("Ange en utbredning.")
        elif not p["aoi"].valueAsText:
            p["aoi"].setErrorMessage("Ange ett polygonlager.")

        metric = METRIC_CODE.get(p["metric"].valueAsText, "NRD")
        bands = None
        try:
            bands = _parse_bands(p["bands"].value)
        except ValueError as exc:
            p["bands"].setErrorMessage(str(exc))
        if bands and metric == "NRD":
            zero = [_band_label(lo, hi) for lo, hi in bands if lo < 0.05]
            if zero:
                p["bands"].setWarningMessage(
                    "Med NRD och undre gräns nära 0 blir värdet nära 1 i varje cell, eftersom "
                    "nästan alla punkter under övre gränsen då ligger i skiktet: {}.".format(
                        ", ".join(zero)))

        cell = p["cell_size"].value
        if cell is not None and not (1 <= cell <= 200):
            p["cell_size"].setErrorMessage("Cellstorleken ska vara mellan 1 och 200 m.")
        elif cell is not None and cell < 10:
            p["cell_size"].setWarningMessage(
                "Laserdata Skog har ungefär 2-5 punkter per m². Vid {} m får en cell bara "
                "ungefär {}-{} punkter, och under tät skog når en liten del av dem ned i de låga "
                "skikten. Värdena blir brusiga.".format(
                    _fmt_m(cell), int(2 * cell * cell), int(5 * cell * cell)))

        mn = p["min_points"].value
        if mn is not None and mn < 1:
            p["min_points"].setErrorMessage("Ange minst 1.")

        mem = p["memory_gb"].value
        if mem is not None and mem < MIN_MEMORY_GB:
            p["memory_gb"].setErrorMessage("Ange minst 0,5 GB, eller lämna tomt för automatiskt.")
        elif mem is not None:
            avail = _available_memory_gb()
            if avail is not None and mem > avail:
                p["memory_gb"].setWarningMessage(
                    "Mer än det lediga arbetsminnet ({} GB). Windows börjar då skriva till disk "
                    "och körningen blir mycket långsam.".format(_fmt_m(round(avail, 1))))

        ws = p["out_workspace"].valueAsText
        prefix = (p["prefix"].valueAsText or "").strip()
        if prefix and not (prefix[0].isalpha() and all(c.isalnum() or c == "_" for c in prefix)):
            p["prefix"].setErrorMessage(
                "Prefixet får bara innehålla bokstäver, siffror och understreck, och "
                "måste börja med en bokstav."
            )
        elif prefix and ws and bands:
            names = _planned_names(prefix, metric, cell or DEFAULT_CELL_SIZE, bands,
                                   bool(p["save_counts"].value))
            existing = [n for n in names if arcpy.Exists(_out_path(ws, n))]
            if existing:
                p["prefix"].setWarningMessage("Skrivs över: " + ", ".join(existing))

        if p["pick_scans"].value:
            aoi = _aoi_bbox_quick(p)
            if aoi is None:
                p["scans"].setWarningMessage("Ange området först, så söks skanningarna fram.")
                return
            scans, err = _scans_for_dialog(aoi)
            if err:
                p["scans"].setErrorMessage(err)
                return
            if not scans:
                p["scans"].setErrorMessage(
                    "Inga skanningar täcker området. Laserdata Skog täcker ungefär 75 % av "
                    "Sverige, men inte fjällen.")
                return
            codes = {_scan_code(v) for v in (p["scans"].values or []) if v}
            chosen = [s for s in scans if s["code"] in codes]
            if not chosen:
                unknown = sorted(c for c in codes if c)
                if unknown:
                    p["scans"].setErrorMessage(
                        "{} täcker inte området. Skanningar som gör det: {}.".format(
                            ", ".join(unknown), ", ".join(s["code"] for s in scans)))
                else:
                    p["scans"].setErrorMessage("Välj minst en skanning, eller avmarkera "
                                               "'Välj skanning själv'.")
                return
            tiles = [t for s in chosen for t in s["tiles"]]
            cov = _coverage(tiles, aoi)
            notes = []
            if cov < 0.995:
                notes.append("De valda skanningarna täcker {} % av områdets utbredning; resten "
                             "blir NoData.".format(int(round(100 * cov))))
            leafless = [s["code"] for s in chosen if s["season"] != "med löv"]
            if leafless:
                notes.append("{} är inte skannad under sommaren. Lövträd och buskar ger då "
                             "få träffar och tätheten i låga skikt underskattas.".format(
                                 ", ".join(leafless)))
            if len({s["season"] for s in chosen}) > 1 or len(chosen) > 1:
                notes.append("Flera skanningar kan ge skarvar där de möts.")
            if notes:
                p["scans"].setWarningMessage(" ".join(notes))

    def execute(self, parameters, messages):
        p = {q.name: q for q in parameters}
        try:
            if p["aoi_mode"].valueAsText == AOI_EXTENT:
                if not p["aoi_extent"].valueAsText:
                    raise ValueError("Ange en utbredning.")
                aoi, note = _extent_polygon(p["aoi_extent"].value, p["aoi_extent"].valueAsText)
                messages.addMessage(note)
            else:
                if not p["aoi"].valueAsText:
                    raise ValueError("Ange ett polygonlager.")
                aoi = _aoi_geometry(p["aoi"].value)
            codes = None
            if p["pick_scans"].value:
                codes = [_scan_code(v) for v in (p["scans"].values or []) if v]
                if not codes:
                    raise ValueError("Välj minst en skanning, eller avmarkera 'Välj skanning "
                                     "själv'.")
            outputs = _run(
                aoi,
                (p["consumer_key"].valueAsText or "").strip(),
                (p["consumer_secret"].valueAsText or "").strip(),
                METRIC_CODE.get(p["metric"].valueAsText, "NRD"),
                _parse_bands(p["bands"].value),
                p["cell_size"].value or DEFAULT_CELL_SIZE,
                p["out_workspace"].valueAsText,
                p["prefix"].valueAsText.strip(),
                codes,
                p["min_points"].value or DEFAULT_MIN_POINTS,
                p["save_counts"].value is not False,
                p["add_group"].value is not False,
                p["memory_gb"].value or None,
                p["workers"].value or None,
                messages,
            )
        except ValueError as exc:
            error = str(exc)
        else:
            error = None
        # Utanför except-blocket, så att Pro inte skriver ut hela kedjan av
        # undantag under det läsbara felmeddelandet.
        if error:
            messages.addErrorMessage(error)
            raise arcpy.ExecuteError

        names = [q.name for q in parameters]
        arcpy.SetParameterAsText(names.index("out_rasters"), ";".join(outputs))

    def postExecute(self, parameters):
        return


# =============================================================================
# Körningens innehåll (separat funktion - går att testa utanför Pro)
# =============================================================================

def _run(aoi, key, secret, metric, bands, cell, workspace, prefix, codes, min_points,
         save_counts, add_group, memory_gb, workers, messages, aprx=None):
    """
    aoi: polygon i SWEREF 99 TM (från _aoi_geometry eller _extent_polygon).
    metric: "NRD" eller "ORD". bands: [[undre, övre], ...] i meter.
    codes: skanningsområden att välja bland, eller None för senaste per ruta.
    Returnerar sökvägarna till de skapade rastren, skikten först.
    aprx: projekt att lägga till rastren i vid test utanför Pro (annars CURRENT).
    """
    if not workspace:
        raise ValueError("Ange en utdata-arbetsyta.")
    if not arcpy.Exists(workspace):
        raise ValueError("Utdata-arbetsytan finns inte: {}".format(workspace))
    memory_gb, workers, auto_note = _resolve_resources(memory_gb, workers)
    if auto_note:
        messages.addMessage(auto_note)
    min_points = max(1, int(min_points))
    codes_set = {c.upper() for c in codes} if codes else None

    ext = aoi.extent
    grid = _grid(ext, cell)
    area_km2 = (ext.XMax - ext.XMin) * (ext.YMax - ext.YMin) / 1e6
    messages.addMessage("Intresseområdets utbredning: {:.2f} km² (SWEREF 99 TM), {} x {} celler "
                        "om {} m.".format(area_km2, grid["width"], grid["height"], _fmt_m(cell)))
    messages.addMessage("Mått: {}. Skikt: {}.".format(
        metric, ", ".join(_band_label(lo, hi) for lo, hi in bands)))
    if area_km2 > LARGE_AREA_KM2:
        messages.addWarningMessage(
            "Stort område ({:.0f} km²). Alla punkter inom utbredningen läses, ungefär 3 "
            "miljoner per km², så körningen kan ta lång tid.".format(area_km2))

    steps = _Steps(6, messages)
    workdir = os.path.join(arcpy.env.scratchFolder, "nrd_" + uuid.uuid4().hex[:8])
    os.makedirs(workdir)
    try:
        # ── 1. Katalogen ─────────────────────────────────────────────────────
        steps.next("söker rutor i Lantmäteriets STAC-katalog")
        wgs = aoi.projectAs(arcpy.SpatialReference(4326)).extent
        items = _stac_search((wgs.XMin, wgs.YMin, wgs.XMax, wgs.YMax))
        tiles = _pick_tiles(items, aoi, codes_set)
        if not tiles:
            if codes_set:
                raise ValueError("De valda skanningarna ({}) täcker inte området.".format(
                    ", ".join(sorted(codes_set))))
            raise ValueError(
                "Inga laserdata för området. Laserdata Skog täcker ungefär 75 % av Sverige, "
                "men inte fjällen.")
        scan_mode = ("Valda skanningar: {}".format(", ".join(sorted(codes_set))) if codes_set
                     else "Senaste skanningen per ruta")
        all_tiles = _pick_tiles(items, aoi)
        missing = {t["key"] for t in all_tiles} - {t["key"] for t in tiles}
        blocks, side = _plan_blocks(aoi, grid, tiles, memory_gb, workers)
        expected = sum(b["expected"] for b in blocks)
        steps.done("{} ruta/rutor. {}. Ungefär {} miljoner punkter väntas.".format(
            len(tiles), scan_mode, _fmt_m(round(expected / 1e6, 1))))
        for t in tiles:
            d = _actual_density(t)
            messages.addMessage("    {}: skanningsområde {}, insamlad {} ({}), {} punkter per m² "
                                "med alla returer.".format(
                t["id"], t["area"] or "okänt", _capture_period(t), _season(_months([t])),
                _fmt_m(round(d, 1)) if d else "?"))
        if missing:
            messages.addWarningMessage(
                "{} ruta/rutor i området saknas bland de valda skanningarna och blir "
                "NoData.".format(len(missing)))
        leafless = sorted({t["area"] for t in tiles if _season(_months([t])) != "med löv"})
        if leafless:
            messages.addWarningMessage(
                "Skanningsområde {} är inte skannat under sommaren (juni-augusti). Lövträd och "
                "buskar ger då få träffar, och tätheten i låga skikt underskattas i lövskog.{}".format(
                    ", ".join(leafless), "" if codes_set else
                    " Välj skanning själv under Skanning om det finns en sommarskanning."))
        if len({_capture_period(t) for t in tiles}) > 1:
            messages.addWarningMessage(
                "Rutorna är skannade vid olika tillfällen. Rastren kan ha en skarv vid "
                "rutgränsen, särskilt om årstid eller år skiljer sig.")
        workers = min(workers, len(blocks))
        messages.addMessage(
            "    {} block på upp till {:.0f} x {:.0f} m, {} åt gången, inom en minnesbudget "
            "på {} GB.".format(len(blocks), side, side, workers, _fmt_m(round(memory_gb, 1))))

        # ── 2. Behörighet ────────────────────────────────────────────────────
        steps.next("hämtar token och kontrollerar behörighet")
        token = _Token(key, secret)
        _check_access(tiles[0]["href"], token.get())
        if importlib.util.find_spec("pdal") is None:
            raise ValueError(
                "Python-paketet pdal saknas i den aktiva miljön ({}). Det ingår i "
                "standardmiljön arcgispro-py3.".format(sys.prefix))
        steps.done()

        # ── 3. Blocken ───────────────────────────────────────────────────────
        steps.next("läser punkter och räknar per cell, {} block, {} åt gången".format(
            len(blocks), workers))
        arcpy.SetProgressor("step", "", 0, 100, 1)
        H, W = grid["height"], grid["width"]
        n_all = np.zeros((H, W), np.int32)
        n_band = [np.zeros((H, W), np.int32) for _ in bands]
        n_below = [np.zeros((H, W), np.int32) for _ in bands]
        n_points = n_ground = n_nohgt = 0
        done_expected = 0.0
        t_blocks = time.time()
        pool = _WorkerPool(workers, workdir)
        failed = True
        try:
            pending = list(blocks)
            running = {}
            for i in range(workers):
                b = pending.pop(0)
                running[i] = b
                pool.submit(i, _block_task(b, token.get(), bands, workdir))
            k = 0
            while running:
                got = pool.get(timeout=1.0)
                if getattr(arcpy.env, "isCancelled", False):
                    raise ValueError("Avbrutet av användaren.")
                if got is None:
                    continue
                i, reply = got
                b = running.pop(i)
                if reply is None:
                    raise RuntimeError("Arbetsprocessen avslutades oväntat under block rad {}, "
                                       "kolumn {}:\n{}".format(b["row"], b["col"], pool.log_tail(i)))
                if not reply["ok"]:
                    raise RuntimeError("Fel i block rad {}, kolumn {}:\n{}".format(
                        b["row"], b["col"], reply["error"]))
                if pending:
                    nb = pending.pop(0)
                    running[i] = nb
                    pool.submit(i, _block_task(nb, token.get(), bands, workdir))

                res = reply["result"]
                k += 1
                n_points += res["points"]
                n_ground += res["ground"]
                n_nohgt += res["no_height"]
                if res["npz"]:
                    g = b["grid"]
                    sl = (slice(b["r0"], b["r0"] + g["height"]),
                          slice(b["c0"], b["c0"] + g["width"]))
                    with np.load(res["npz"]) as z:
                        n_all[sl] += z["n_all"]
                        for j in range(len(bands)):
                            n_band[j][sl] += z["band_{}".format(j)]
                            n_below[j][sl] += z["below_{}".format(j)]
                    os.remove(res["npz"])

                done_expected += b["expected"]
                elapsed = time.time() - t_blocks
                left = elapsed / done_expected * (expected - done_expected) if done_expected else 0
                messages.addMessage("    Block {} av {} (rad {}, kolumn {}): {} punkter.".format(
                    k, len(blocks), b["row"], b["col"], _fmt_count(res["points"])))
                steps.label("block {} av {} klara{}".format(
                    k, len(blocks), ", ca {} kvar".format(_fmt_duration(left)) if running else ""))
                arcpy.SetProgressorPosition(min(100, int(100 * done_expected / max(expected, 1))))
            failed = False
        finally:
            pool.close(kill=failed)

        if n_points == 0:
            raise ValueError("Inga punkter inom området.")
        note = "{} punkter, varav {} markpunkter.".format(_fmt_count(n_points),
                                                          _fmt_count(n_ground))
        if n_nohgt:
            note += " {} punkter ({} %) saknade markpunkter runt omkring och räknades inte.".format(
                _fmt_count(n_nohgt), _fmt_m(round(100.0 * n_nohgt / n_points, 1)))
        steps.done(note)

        # ── 4. Kvoter ────────────────────────────────────────────────────────
        steps.next("beräknar {} för {} skikt".format(metric, len(bands)))
        ratios, stats = [], []
        covered = n_all > 0
        for j, (lo, hi) in enumerate(bands):
            den = n_below[j] if metric == "NRD" else n_all
            ok = den >= min_points
            r = np.full((H, W), _nw.NODATA, np.float32)
            r[ok] = n_band[j][ok] / den[ok]
            vals = r[ok]
            few = int((covered & ~ok).sum())
            st = {"valid": int(ok.sum()), "few": few,
                  "few_pct": round(100.0 * few / max(1, int(covered.sum())), 1),
                  "mean": round(float(vals.mean()), 3) if vals.size else None,
                  "median": round(float(np.median(vals)), 3) if vals.size else None}
            ratios.append(r)
            stats.append(st)
            msg = "    {} {}: {} celler med värde".format(metric, _band_label(lo, hi),
                                                          _fmt_count(st["valid"]))
            if vals.size:
                msg += ", medel {}, median {}".format(_fmt_m(st["mean"]), _fmt_m(st["median"]))
            msg += "; {} celler ({} %) har färre än {} punkter i nämnaren och blir NoData.".format(
                _fmt_count(few), _fmt_m(st["few_pct"]), min_points)
            messages.addMessage(msg)
            if st["few_pct"] > 25:
                messages.addWarningMessage(
                    "{} {}: {} % av cellerna med punkter har för få punkter i nämnaren. Överväg "
                    "större celler eller en lägre gräns under Tillförlitlighet.".format(
                        metric, _band_label(lo, hi), _fmt_m(st["few_pct"])))
        steps.done()

        # ── 5. Skriv, klipp, metadata ────────────────────────────────────────
        jobs = []
        for j, (lo, hi) in enumerate(bands):
            jobs.append((_out_name(prefix, metric, cell, lo=lo, hi=hi), ratios[j],
                         {"kind": "ratio", "lo": lo, "hi": hi, "stats": stats[j]}))
        if save_counts:
            if metric == "NRD":
                for hi in sorted({hi for _lo, hi in bands}):
                    j = next(j for j, b in enumerate(bands) if b[1] == hi)
                    jobs.append((_out_name(prefix, metric, cell, "under", hi=hi), n_below[j],
                                 {"kind": "count", "lo": None, "hi": hi}))
            else:
                jobs.append((_out_name(prefix, metric, cell, "alla"), n_all,
                             {"kind": "count", "lo": None, "hi": None}))

        steps.next("skriver {} raster, klipper till området och skriver metadata".format(
            len(jobs)))
        arcpy.SetProgressor("step", "", 0, len(jobs), 1)
        sr = _tm()
        old_ocs = arcpy.env.outputCoordinateSystem
        old_overwrite = arcpy.env.overwriteOutput
        outputs = []
        run_info = {"metric": metric, "bands": bands, "cell": cell, "min_points": min_points,
                    "scan_mode": scan_mode, "points": n_points, "ground": n_ground,
                    "no_height": n_nohgt, "extent": ext,
                    "created": datetime.date.today().isoformat(),
                    "tool_label": PunkttathetISkikt().label}
        try:
            arcpy.env.outputCoordinateSystem = sr
            arcpy.env.overwriteOutput = True
            existing = [_out_path(workspace, name) for name, _a, _i in jobs
                        if arcpy.Exists(_out_path(workspace, name))]
            if existing:
                try:
                    _release_from_maps(existing, messages, aprx)
                except Exception as exc:
                    messages.addWarningMessage(
                        "Kunde inte ta bort gamla lager ur kartan: {}".format(exc))
            # Klippmallen i körningens temp-mapp, inte memory\ (se CLAUDE.md).
            clip_fc = arcpy.management.CopyFeatures([aoi], os.path.join(workdir, "aoi.shp"))[0]
            rect = "{} {} {} {}".format(ext.XMin, ext.YMin, ext.XMax, ext.YMax)
            ll = arcpy.Point(grid["origin_x"], grid["origin_y"])
            for k, (name, arr, info) in enumerate(jobs, 1):
                steps.label("skriver {} ({} av {})".format(name, k, len(jobs)))
                nodata = _nw.NODATA if info["kind"] == "ratio" else COUNT_NODATA
                tmp = os.path.join(workdir, "r{}.tif".format(k))
                # Rad 0 i räknetabellerna är i söder; rastret skrivs norr uppåt.
                arcpy.NumPyArrayToRaster(np.ascontiguousarray(arr[::-1]), ll, cell, cell,
                                         nodata).save(tmp)
                dst = _out_path(workspace, name)
                arcpy.management.Clip(tmp, rect, dst, clip_fc, str(nodata),
                                      "ClippingGeometry", "NO_MAINTAIN_EXTENT")
                try:
                    _write_raster_metadata(dst, info, tiles, run_info)
                except Exception as exc:
                    messages.addWarningMessage(
                        "Kunde inte skriva metadata för {}: {}".format(name, exc))
                outputs.append(dst)
                messages.addMessage("    Skapade {}.".format(dst))
                arcpy.SetProgressorPosition(k)
        finally:
            arcpy.env.outputCoordinateSystem = old_ocs
            arcpy.env.overwriteOutput = old_overwrite
        steps.done()

        # ── 6. Kartan ────────────────────────────────────────────────────────
        steps.next("lägger till rastren i kartan")
        try:
            _add_to_map(outputs, _group_name(prefix, metric, cell) if add_group else None,
                        messages, aprx)
        except Exception as exc:
            messages.addWarningMessage("Kunde inte lägga till rastren i kartan: {}".format(exc))
        steps.done()
        messages.addMessage("Klart på {}.".format(_fmt_duration(time.time() - steps.t0)))
        return outputs
    finally:
        arcpy.ResetProgressor()
        # arcpy kan hålla mappen öppen som arbetsyta. Släpp just den mappen
        # (utan argument skulle alla arbetsytor i Pro-sessionen släppas).
        try:
            arcpy.management.ClearWorkspaceCache(workdir)
        except Exception:
            pass
        shutil.rmtree(workdir, ignore_errors=True)
