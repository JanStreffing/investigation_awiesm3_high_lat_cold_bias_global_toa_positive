#!/usr/bin/env python3
"""Regenerate LPJ-GUESS peat-fraction forcing at the model grid's own precision.

WHY.  The canonical land-cover rewrite (a270270's fixed-peat work) requires every
coupled grid cell to find its own row in peat_frac.txt within 5e-5 degree per
axis, and explicitly forbids zero or nearest-neighbour substitution.  The files
shipped under input/lpj-guess/peat were written by create_peat_files.py with a
"{lon:.2f}\t{lat:.2f}" format, so their coordinates are rounded to two decimals
and miss that tolerance by up to ~2e-3.  Every coupled run then dies at LPJ-GUESS
initialisation with "Fixed peat map lookup failed".

COORDINATE SOURCE.  OASIS grids.nc, not the regridded land-use NetCDFs.  LPJ-GUESS
receives its grid through OASIS, so grids.nc *is* the grid the model asks about,
stored in float64.  The land-use NetCDFs carry the same grid rounded to ~4
decimals, which also misses 5e-5 (checked: dlon 2.3e-4 at 81.8168N).

PEAT VALUES.  Nearest-neighbour from the TL255 source file, exactly as before.
Only the coordinates change; the science does not.
"""

import argparse
import os
import numpy as np
import netCDF4
from scipy.spatial import cKDTree

SOURCE_PEAT = "/work/ab0246/a270092/input/lpj-guess/peat/TL255_peat_frac.txt"


def read_source(path):
    d = np.genfromtxt(path, names=True, dtype=None, encoding=None)
    return np.asarray(d["Lon"], float), np.asarray(d["Lat"], float), np.asarray(d["Peat_Frac"], float)


def read_oasis_grid(grids_nc, prefix):
    ds = netCDF4.Dataset(grids_nc)
    try:
        lon = np.asarray(ds.variables[f"{prefix}.lon"][:], dtype=float).ravel()
        lat = np.asarray(ds.variables[f"{prefix}.lat"][:], dtype=float).ravel()
    finally:
        ds.close()
    # LPJ-GUESS reports longitudes in [-180, 180); match that convention.
    # ">=", not ">": the dateline column sits at exactly 180 in grids.nc, and
    # LPJ-GUESS asks for it as -180.  Mapping only lon > 180 leaves it at +180,
    # which is 360 degrees from what the lookup wants and fails the 5e-5 test.
    lon = np.where(lon >= 180.0, lon - 360.0, lon)
    return lon, lat


def nearest_peat(src_lon, src_lat, src_val, dst_lon, dst_lat):
    """Plain 2-D nearest neighbour in degrees.

    Deliberately identical to create_peat_files.py's griddata(method="nearest"),
    which is Euclidean in (lon, lat) and handles neither the dateline nor the
    convergence of meridians.  Matching on the sphere instead would be more
    defensible, but it changes the peat fraction in 73 TCO95 cells by up to 0.85,
    and this script exists to fix coordinate precision, not to re-derive the
    forcing.  Change the interpolation as its own deliberate step if wanted.
    """
    tree = cKDTree(np.column_stack([src_lon, src_lat]))
    _, idx = tree.query(np.column_stack([dst_lon, dst_lat]), k=1)
    return src_val[idx]


def write_peat(path, lon, lat, val):
    # The dateline column is emitted under BOTH conventions, once at -180 and
    # once at +180.  The two files already shipped disagree with each other --
    # TCO95_peat_frac.txt spans [-176.9, 180.0] and TL255_peat_frac.txt spans
    # [-180.0, 179.56] -- so which one LPJ-GUESS asks with cannot be inferred,
    # and a miss here is a fatal "Fixed peat map lookup failed".  A duplicate row
    # is harmless to a lookup; a missing one stops the run.
    dl = np.isclose(lon, -180.0)
    if dl.any():
        lon = np.concatenate([lon, np.full(dl.sum(), 180.0)])
        lat = np.concatenate([lat, lat[dl]])
        val = np.concatenate([val, val[dl]])

    with open(path, "w") as f:
        f.write("Lon\tLat\tPeat_Frac\n")
        for x, y, v in zip(lon, lat, val):
            # 8 decimals: the model's coordinates arrive as float32 promoted to
            # double, so this reproduces them well inside the 5e-5 tolerance.
            f.write(f"{x:.8f}\t{y:.8f}\t{v:.8f}\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--grids-nc", required=True, help="OASIS grids.nc from a staged run")
    p.add_argument("--grid-prefix", required=True, help="atmosphere grid name, e.g. A096 (TCO95) or A128 (TL255)")
    p.add_argument("--out", required=True)
    p.add_argument("--source", default=SOURCE_PEAT)
    a = p.parse_args()

    slon, slat, sval = read_source(a.source)
    dlon, dlat = read_oasis_grid(a.grids_nc, a.grid_prefix)
    val = nearest_peat(slon, slat, sval, dlon, dlat)
    write_peat(a.out, dlon, dlat, val)
    print(f"{a.out}: {len(dlon)} cells from {a.grid_prefix}, "
          f"peat_frac min={val.min():.6f} max={val.max():.6f} mean={val.mean():.6f}")


if __name__ == "__main__":
    main()
