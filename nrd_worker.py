# -*- coding: utf-8 -*-
"""
nrd_worker.py

PDAL-delen av LaserdataPunkttathetNRD.pyt: läser ett block ur Lantmäteriets COPC-filer,
räknar varje punkts höjd över mark och räknar punkter per cell och höjdskikt.

Körs som egen process, en per samtidigt block, och importerar aldrig arcpy.
Flera PDAL-pipelines i trådar i samma process som arcpy kraschade processen
(access violation i arcpy, uppmätt i Pro 3.6; båda laddar Esris gdal_e.dll).
I egna processer delar PDAL inget med arcpy.

Protokoll: en uppgift per rad som JSON på stdin, ett svar per rad på stdout
med prefixet RESULT_PREFIX (annan utskrift från PDAL ignoreras av anroparen).
Processen avslutas när stdin stängs.

Uppgift:  {"block": {...}, "token": "...", "bands": [[lo, hi], ...], "workdir": "..."}
Svar:     {"ok": true, "result": {...}} eller {"ok": false, "error": "traceback"}

Höjd över mark: markpunkterna (klass 2) i blocket plus marginal trianguleras
(filters.delaunay) och rastreras till ett fint rutnät (DTM_RES), som sedan
interpoleras bilinjärt i varje punkts läge. TIN:en är linjär inom varje
triangel, så felet mot TIN:en själv är högst några cm vid DTM_RES = 0,5 m.
Punkter utanför trianguleringen (inga markpunkter runt omkring, t.ex. mitt i
en sjö) får ingen höjd och räknas inte.

Räkning: för varje cell i blockets rutnät och varje skikt [lo, hi)
    n_band  = punkter med lo <= h < hi
    n_below = punkter med h < hi   (nämnaren i NRD, inklusive punkter under marken)
    n_all   = alla punkter med höjd (nämnaren i ORD, samma för alla skikt)
Skikten är halvöppna, så en punkt exakt på en gräns hamnar i det övre skiktet.
Blockets celler är halvöppna på samma sätt, så varje punkt räknas i exakt ett
block trots att blocken läses med marginal.
"""

import json
import os
import sys
import traceback

import numpy as np

RESULT_PREFIX = "@@NRD_RESULT "

NODATA = -9999.0
CLASS_GROUND = 2
NOISE_CLASSES = (7, 18)
DIMS = ["X", "Y", "Z", "Classification"]

# Upplösning (m) för det interna markrutnätet som punkthöjderna räknas från.
DTM_RES = 0.5


def import_pdal():
    # Pro:s PDAL-bygge har curl utan CA-certifikat. Utan detta misslyckas varje
    # HTTPS-anslutning och arbiter försöker om i all oändlighet.
    try:
        import certifi
        os.environ.setdefault("ARBITER_CA_INFO", certifi.where())
    except ImportError:
        pass
    import pdal
    return pdal


def ground_grid(pdal, ground, x0, y0, x1, y1, path):
    """
    TIN av markpunkterna, rastrerad till DTM_RES. Returnerar (array, x0, y_top)
    med NaN utanför trianguleringen. En enda array: filters.delaunay
    trianguerar varje array för sig.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    width = int(np.ceil((x1 - x0) / DTM_RES))
    height = int(np.ceil((y1 - y0) / DTM_RES))
    stages = [
        {"type": "filters.delaunay"},
        {"type": "filters.faceraster", "resolution": DTM_RES, "origin_x": x0, "origin_y": y0,
         "width": width, "height": height},
        {"type": "writers.raster", "filename": path, "data_type": "float32", "nodata": NODATA},
    ]
    pdal.Pipeline(json.dumps(stages), arrays=[ground]).execute()
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    arr = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    ds = None
    gdal.GetDriverByName("GTiff").Delete(path)
    arr[arr == NODATA] = np.nan
    # writers.raster skriver norr uppåt: gt[3] är övre kanten, gt[5] negativ.
    return arr, gt[0], gt[3], gt[1]


def sample_bilinear(grid, gx0, gytop, res, x, y):
    """Bilinjär interpolation i cellmittpunkterna. NaN där någon granne saknas."""
    fx = (x - gx0) / res - 0.5
    fy = (gytop - y) / res - 0.5
    h, w = grid.shape
    c0 = np.clip(np.floor(fx).astype(np.int64), 0, w - 2)
    r0 = np.clip(np.floor(fy).astype(np.int64), 0, h - 2)
    tx = np.clip(fx - c0, 0.0, 1.0)
    ty = np.clip(fy - r0, 0.0, 1.0)
    z00 = grid[r0, c0]
    z01 = grid[r0, c0 + 1]
    z10 = grid[r0 + 1, c0]
    z11 = grid[r0 + 1, c0 + 1]
    return (z00 * (1 - tx) * (1 - ty) + z01 * tx * (1 - ty)
            + z10 * (1 - tx) * ty + z11 * tx * ty)


def count_cells(col, row, h, w, bands, height):
    """
    Räkna punkter per cell. col/row är cellindex (rad 0 i söder), height är
    punkternas höjd över mark. Returnerar n_all (h, w) och per skikt
    (n_band, n_below), alla int32.
    """
    ncell = h * w
    idx = row * w + col
    n_all = np.bincount(idx, minlength=ncell).astype(np.int32).reshape(h, w)
    out = []
    for lo, hi in bands:
        below = height < hi
        band = below & (height >= lo)
        n_below = np.bincount(idx[below], minlength=ncell).astype(np.int32).reshape(h, w)
        n_band = np.bincount(idx[band], minlength=ncell).astype(np.int32).reshape(h, w)
        out.append((n_band, n_below))
    return n_all, out


def process_block(pdal, task):
    """Läs ett block, räkna höjd över mark och punkter per cell och skikt."""
    blk = task["block"]
    grid = blk["grid"]
    bands = task["bands"]
    cell = grid["resolution"]
    gx0, gy0, w, h = grid["origin_x"], grid["origin_y"], grid["width"], grid["height"]

    arrays = []
    for rd in blk["reads"]:
        stages = [{"type": "readers.copc", "bounds": rd["bounds"],
                   "filename": {"path": rd["href"],
                                "headers": {"Authorization": "Bearer " + task["token"]}}}]
        pipe = pdal.Pipeline(json.dumps(stages))
        pipe.execute()
        for a in pipe.arrays:
            keep = ~np.isin(a["Classification"], NOISE_CLASSES)
            arrays.append(np.stack([a["X"][keep], a["Y"][keep], a["Z"][keep],
                                    a["Classification"][keep].astype(np.float64)], axis=1))
        del pipe
    pts = np.concatenate(arrays) if arrays else np.zeros((0, 4))
    del arrays

    result = {"points": 0, "ground": 0, "no_height": 0, "npz": None}
    x, y = pts[:, 0], pts[:, 1]
    inside = (x >= blk["x0"]) & (x < blk["x1"]) & (y >= blk["y0"]) & (y < blk["y1"])
    is_ground = pts[:, 3] == CLASS_GROUND
    result["points"] = int(inside.sum())
    result["ground"] = int((inside & is_ground).sum())
    if not inside.any():
        return result

    ground = pts[is_ground]
    tag = "{:03d}_{:03d}".format(blk["row"], blk["col"])
    height = np.full(int(inside.sum()), np.nan)
    px, py, pz = x[inside], y[inside], pts[inside, 2]
    del pts, x, y
    if len(ground) >= 3:
        g = np.zeros(len(ground), dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8")])
        g["X"], g["Y"], g["Z"] = ground[:, 0], ground[:, 1], ground[:, 2]
        del ground
        rx0, ry0, rx1, ry1 = blk["read_bbox"]
        path = os.path.join(task["workdir"], "ground_{}.tif".format(tag)).replace("\\", "/")
        dtm, dx0, dytop, dres = ground_grid(pdal, g, rx0, ry0, rx1, ry1, path)
        height = pz - sample_bilinear(dtm, dx0, dytop, dres, px, py)
        del dtm, g

    ok = ~np.isnan(height)
    result["no_height"] = int((~ok).sum())
    col = np.floor((px[ok] - gx0) / cell).astype(np.int64)
    row = np.floor((py[ok] - gy0) / cell).astype(np.int64)
    # Punkter precis på blockets övre kant kan avrundas till cellen utanför.
    col = np.clip(col, 0, w - 1)
    row = np.clip(row, 0, h - 1)
    n_all, per_band = count_cells(col, row, h, w, bands, height[ok])

    npz = os.path.join(task["workdir"], "counts_{}.npz".format(tag))
    arrs = {"n_all": n_all}
    for k, (n_band, n_below) in enumerate(per_band):
        arrs["band_{}".format(k)] = n_band
        arrs["below_{}".format(k)] = n_below
    np.savez(npz, **arrs)
    result["npz"] = npz
    return result


def main():
    # Sökvägar kan innehålla å, ä och ö; läs och skriv alltid UTF-8.
    sys.stdin.reconfigure(encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    pdal = import_pdal()
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            out = {"ok": True, "result": process_block(pdal, json.loads(line))}
        except Exception:
            out = {"ok": False, "error": traceback.format_exc()}
        sys.stdout.write(RESULT_PREFIX + json.dumps(out) + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
