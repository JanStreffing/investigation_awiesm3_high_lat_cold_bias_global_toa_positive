#!/usr/bin/env python3
"""Observed peatland areal fraction over the Siberian evaluation box.

Source: Hugelius et al. (2020) PNAS, "Large stocks of peatland carbon and
nitrogen are vulnerable to permafrost thaw", 10 km northern peatland maps
(Bolin Centre Database, https://bolin.su.se/data/hugelius-2020,
file hugelius-2020-2.zip).  Two grids are used:

  Histel_fraction.tif    percent of gridcell that is permafrost peatland
  Histosol_fraction.tif  percent of gridcell that is non-permafrost peatland

Their sum is total peatland areal coverage in percent.  Byte rasters, value
128 = nodata (ocean / outside domain), 0-100 = valid coverage percent.

The rasters are warped to a regular 0.05 deg lat/lon grid over the box and
averaged with cos(lat) area weights over land (= valid) pixels only.

Usage:
    python3 peat_fraction_siberia.py <dir_with_warped_tifs>

The warped files are produced by:
    gdalwarp -t_srs EPSG:4326 -te 60 55 180 75 -tr 0.05 0.05 -r near \
             -srcnodata none -dstnodata 200 -ot Byte in.tif out.tif
"""
import sys

import numpy as np
import rasterio

NODATA = 128

LAT0, LAT1 = 55.0, 75.0
LON0, LON1 = 60.0, 180.0
DX = 0.05


def load(path):
    with rasterio.open(path) as src:
        return src.read(1).astype(float)


def main(base):
    histel = load(f"{base}/Histel_fraction_v2.tif")
    histosol = load(f"{base}/Histosol_fraction_v2.tif")

    land = (histel != NODATA) & (histosol != NODATA) & (histel <= 100) & (histosol <= 100)
    peat = np.where(land, histel + histosol, 0.0)

    ny, nx = histel.shape
    lat = LAT1 - (np.arange(ny) + 0.5) * DX
    lon = LON0 + (np.arange(nx) + 0.5) * DX
    w = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nx))

    # 0.05 deg cell area in km2, R = 6371 km
    r = 6371.0
    cell = (np.deg2rad(DX) * r) ** 2 * np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nx))

    def report(tag, sel):
        m = land & sel
        if not m.any():
            print(f"{tag}: no land")
            return
        area_land = cell[m].sum()
        area_peat = (cell * peat / 100.0)[m].sum()
        frac = 100.0 * area_peat / area_land
        pf_permafrost = 100.0 * (cell * histel / 100.0)[m].sum() / area_land
        print(
            f"{tag:28s} land {area_land/1e6:6.3f}e6 km2   "
            f"peat {area_peat/1e6:6.3f}e6 km2   frac {frac:5.2f}%   "
            f"(of which Histel/permafrost {pf_permafrost:5.2f}%)"
        )

    lon2d = lon[None, :] * np.ones((ny, 1))
    lat2d = lat[:, None] * np.ones((1, nx))

    print("Hugelius et al. 2020 PNAS peatland fraction, box 55-75N 60-180E")
    print()
    report("FULL BOX 60-180E", np.ones_like(peat, bool))
    report("WEST 60-90E (WSL)", lon2d < 90)
    report("WEST-CENTRAL 60-105E", lon2d < 105)
    report("EAST 105-180E", lon2d >= 105)
    report("90-135E", (lon2d >= 90) & (lon2d < 135))
    report("135-180E", lon2d >= 135)
    print()
    for a, b in [(55, 60), (60, 65), (65, 70), (70, 75)]:
        report(f"lat {a}-{b}N", (lat2d >= a) & (lat2d < b))
    print()

    # How much would a DOMINANT-soil-type ancillary at ~TCO95 (~100 km) mark
    # as organic?  Aggregate to 1.0 deg lat x 2.0 deg lon (~110 x 95 km at 65N)
    # and count coarse cells whose peat fraction exceeds a threshold.
    fy, fx = int(round(1.0 / DX)), int(round(2.0 / DX))
    cy, cx = ny // fy, nx // fx
    pw = (peat * cell * land)[: cy * fy, : cx * fx].reshape(cy, fy, cx, fx).sum(axis=(1, 3))
    aw = (cell * land)[: cy * fy, : cx * fx].reshape(cy, fy, cx, fx).sum(axis=(1, 3))
    ac = cell[: cy * fy, : cx * fx].reshape(cy, fy, cx, fx).sum(axis=(1, 3))
    ok = aw > 0.5 * ac  # coarse cell is majority land
    coarse = np.where(ok, pw / np.maximum(aw, 1e-9), 0.0)
    print("Dominant-class emulation on ~1x2 deg (~100 km) cells, majority-land only:")
    for thr in (30, 40, 50, 60):
        sel = ok & (coarse > thr)
        frac_cells = 100.0 * aw[sel].sum() / aw[ok].sum()
        print(f"  land area in cells with peat > {thr}%:  {frac_cells:5.2f}%")
    print(f"  area-weighted mean peat fraction on coarse land cells: "
          f"{(coarse[ok] * aw[ok]).sum() / aw[ok].sum():.2f}%")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")

# ---------------------------------------------------------------------------
# MODEL SIDE, and a correction worth recording (2026-08-04).
#
# The model's own soil-type field (GRIB `slt`, ICMGGawi3INIT) over the Siberian
# box, masked with the IFS `lsm` from the SAME file so the grids are aligned:
#
#     medium 76.9 %, coarse 5.6 %, medium-fine 5.4 %, fine 4.0 %,
#     very fine 1.8 %, ORGANIC (type 6) 6.18 %,  water 0.2 %
#     organic by longitude: west 60-90E 7.1 %, mid 90-135E 8.8 %, east 135-180E 0.6 %
#
# An earlier pass reported organic as 0.10 % and "25.9 % unclassified". BOTH were
# artefacts of a LATITUDE-ORDERING MISMATCH: the GRIB-derived fields come out
# north->south, the campaign's remapped land mask (runs.LSMF) is south->north, so
# masking one with the other sampled the opposite hemisphere. Always force a
# common orientation before combining a GRIB-derived field with the campaign mask
# -- the campaign's own scripts are safe because mask and data come from the same
# remapped stream, but anything pulled straight from ICMGG*INIT is not.
#
# CONSEQUENCE. Against the observed areal fraction above (12.5 % box-wide) the
# model's 6.18 % is roughly half; against the fairer dominant-class target for a
# ~100 km grid (2.8 % of land in cells >50 % peat) it is about double. The field
# is therefore NOT badly wrong in the box mean -- the earlier "20-125x too low"
# framing rested on the bad number. What IS wrong is the PATTERN: the West
# Siberian Lowland should be the peat maximum (~26 %) and the model puts only
# 7.1 % there while over-weighting 90-135E.
#
# So reclassification is not a lever: the gap is a few percentage points, and the
# sensitivity is ~0.10-0.20 K for converting ALL land, i.e. order 0.01 K here --
# far below the +-0.242 K detection floor.
# ---------------------------------------------------------------------------
