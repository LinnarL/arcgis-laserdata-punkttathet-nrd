# Laserdata NRD

ArcGIS Pro Python toolbox that maps relative point density in height bands above ground from
Lantmäteriet's Laserdata Nedladdning, skog. For each band and raster cell it computes either
metric from Campbell et al. (2018):

- **NRD** (normalised relative point density) = points in the band / points from the ground up to
  the band's upper limit
- **ORD** (overall relative point density) = points in the band / all points

NRD only counts the pulses that reached the band or went below it, so it is less affected by a
dense canopy above. In the paper, NRD 0.15-1.85 m explained 44 % of field-measured understory
density and ORD 14 %. With nothing above the band the two are equal.

Campbell, M.J., Dennison, P.E., Hudak, A.T., Parham, L.M., Butler, B.W. (2018). Quantifying
understory vegetation density using small-footprint airborne lidar. Remote Sensing of
Environment 215, 330-342. https://doi.org/10.1016/j.rse.2018.06.023

Only the points inside the area's bounding box are read, straight from Lantmäteriet's
cloud-optimised files. Whole tiles are never downloaded.

## Requirements

- ArcGIS Pro 3.x, Basic licence is enough. Developed and tested on 3.6 with Python 3.13.
- No extra packages. `pdal`, `numpy`, `certifi` and GDAL's Python bindings ship with the default
  `arcgispro-py3` environment.
- A consumer key and secret from Lantmäteriet's API portal with **both** the API `STAC-hojd` and
  an order of [Laserdata Nedladdning, skog](https://geotorget.lantmateriet.se/geodataprodukter/laserdata-nedladdning-skog-api)
  on Geotorget. The same key works for the Laserdata Skog toolbox. With only the API, the tool
  fails with HTTP 403 and says so.

## Install

1. Clone or download this repo. Keep `LaserdataNRD.pyt` and `nrd_worker.py` in the same folder;
   the toolbox runs the worker file in separate processes.
2. In ArcGIS Pro: Catalog, Toolboxes, Add Toolbox, select `LaserdataNRD.pyt`.
3. Open Lantmäteriet Laserdata NRD, Punkttäthet i höjdskikt (NRD/ORD).

## Parameters

| Parameter | Default | Notes |
|---|---|---|
| Avgränsa området med | Polygoner i ett lager | Polygon layer or extent |
| Intresseområde (polygoner) | - | Any coordinate system. All features, or the selection, are merged, and the rasters are clipped to the shapes |
| Utbredning | - | Display extent, a layer's extent, a drawn rectangle or typed coordinates (read in the active map's coordinate system) |
| Mått | NRD | NRD or ORD |
| Höjdskikt (m över mark) | 0.15-1.85, 2-5, 5-10 | Table typed in the dialog, one row per band. Lower limit included, upper excluded. Bands may overlap. Up to 20 bands, top at most 60 m |
| Cellstorlek (m) | 20 | 1-200 m. Warns below 10 m |
| Utdata-arbetsyta | project geodatabase | Geodatabase or folder (GeoTIFF) |
| Namnprefix | `laser` | Rasters are `<prefix>_nrd_015_185` etc., limits in cm. Existing ones are overwritten, with a warning in the dialog |
| Samla resultatet i en grupp i kartan | on | Group named `<prefix> - NRD 20 m` in the active map |
| Välj skanning själv | off | Under Skanning. Off uses the newest scan per tile |
| Skanningar | - | Under Skanning. Filled from the STAC catalogue when the box is ticked, see below |
| Minsta antal punkter i nämnaren | 20 | Under Tillförlitlighet. Cells with fewer points in the denominator are NoData |
| Spara punktantal (nämnaren) | on | Under Tillförlitlighet. Integer rasters `<prefix>_antal_under_185` (NRD, one per distinct upper limit) or `<prefix>_antal_alla` (ORD) |
| Consumer key / Consumer secret | - | Under Inloggning. The secret is a hidden field |
| Minnesbudget (GB) / Block samtidigt | 4 / 4 | Under Avancerat. As in Laserdata Skog |

Every parameter has a tooltip in the dialog. The text lives in `TOOLTIPS` in the `.pyt`, which
writes it to `LaserdataNRD.PunkttathetISkikt.pyt.xml` when the toolbox loads.

## Choosing a scan

Many areas have been scanned twice, and the catalogue keeps both scans, with different scanning
area codes (for example `19A024` from January 2020 and `26A024` from April 2026). By default the
tool uses the newest scan per tile. Tick **Välj skanning själv** and the dialog searches the open
STAC catalogue (no login, about a second, repeated only when the area changes) and lists each scan
covering the area:

```
26A024 - april 2026, utan löv, 3,8 p/m², täcker 64 % av utbredningen
22A025 - januari 2022, utan löv, 3 p/m², täcker 16 % av utbredningen
19A024 - december 2019 - januari 2020, utan löv, 3,4 p/m², täcker 64 % av utbredningen
```

Pick several when the area spans several scanning areas. Where two chosen scans cover the same
tile, the newer is used. The dialog warns when the selection leaves part of the area uncovered or
includes a scan outside June-August.

The choice matters for low bands. A newer scan is usually denser, and the paper found pulse
density the strongest factor for accuracy. But many scans are made without leaves (all
three in the example above, from southern Småland), and then deciduous shrubs and trees return few points and understory
density is underestimated. "Med löv" means June-August, "lövläge osäkert" May, September and
October, "utan löv" November-April.

## Output

One float raster per band (0-1, NoData where the denominator has too few points), plus the
optional count rasters, all in SWEREF 99 TM, snapped to whole multiples of the cell size and
clipped to the area. They are added to the active map, by default in one group. Layers that show
rasters about to be overwritten are removed from the map first; otherwise the overwrite fails on
a schema lock.

Each raster's metadata (Catalog, View Metadata) records the settings used: metric and formula,
the band and all bands in the run, cell size, minimum point count, how many cells fell below it,
mean and median, scan choice, which points were used and how height above ground was computed,
point counts, the tool version and date. It also has a table of the source tiles with capture
dates, leaf state, flying height and point density, the reference to the paper, credits to
Lantmäteriet and a link to the terms of use.

## How it is computed

1. Points in each block plus a 20 m margin are read with PDAL. Noise (classes 7 and 18) is dropped;
   every other class and every return is kept, as in the paper.
2. Ground points (class 2) are triangulated (Delaunay) and rasterised to 0.5 m. Each point's height
   above ground is its Z minus the bilinear interpolation of that grid at its position. Points
   with no ground triangle around them (typically inside a lake) get no height and are not
   counted; the log says how many.
3. Per cell and band: `n_band` = points with lower <= h < upper, `n_below` = points with h < upper,
   `n_all` = all points with a height. Points below the ground surface (h < 0) count as below the
   band. NRD = `n_band / n_below`, ORD = `n_band / n_all`.
4. Each block returns only these counts. The rasters are assembled in one pass with numpy, written,
   and clipped.

Blocked and unblocked runs give identical results, and the output matches a brute-force count of
the same points cell for cell (verified on a synthetic cloud with known answers, in 1, 4 and 9
blocks).

## About the data

- Lantmäteriet publishes a nominal density of 1-2 points per m². The actual number of points,
  from the catalogue's point count per tile, is 2.3-5 per m² (median about 3) because it includes
  all returns. That is what the tool counts, and what the scan list shows. The paper had 16.4 per
  m², so cells have to be much larger than its 11 x 2 m transects: at 20 m a cell holds about
  1 000-2 000 points. The paper mapped at 30 m.
- A lower limit below about 0.15 m picks up misclassified ground points, and the ground model's
  own error is a few cm. The paper names ground classification as a main error source.
- NRD also depends on canopy density and height; the paper found both reduce how well NRD
  predicts the understory. The count rasters show where the denominator is thin.
- Tiles scanned on different dates can leave a seam. The tool warns when this happens.
- Coverage is about 75 % of Sweden. The mountains are not included.
