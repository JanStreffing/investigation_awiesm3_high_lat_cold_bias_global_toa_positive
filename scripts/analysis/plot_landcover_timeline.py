#!/usr/bin/env python3
"""
Global cropland (and pasture) evolution in an AWI-ESM3 TCO319 run.

Reads the yearly CMIP landcover fraction outputs written by LPJ-GUESS:
  cropFrac_yearly.out, pastureFrac_yearly.out   (column "Total", PERCENT of
  gridcell -- CMIPoutput declares these as "[%]", so they are divided by 100).

Areas are exact for the O320 octahedral grid that carries the Tco319
truncation: the Gauss-Legendre weight w_i of latitude row i is precisely that
row's band width in sin(latitude), so cell area = 2*pi*R^2 * w_i / n_i with
n_i = 20 + 4*(i-1) longitudes, mirrored about the equator.

This doubles as the check on whether land-use change is actually active: a run
with dynamic LUH3 land use shows cropland moving year on year, a run with it
disabled or pinned to a fixed year shows a flat line.
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_RUN = ("/work/bb1469/a270089/runtime/awiesm3-v3.4.2/"
               "AWI-ESM3-VEG-HR-CMIP7-1pctCO2_1949")
PLOTDIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "plots")

R_EARTH = 6371.0088   # km
N_OCT = 320           # O320: 640 Gaussian latitude rows, 421,120 points
PERCENT = 100.0       # cropFrac/pastureFrac are written in percent


def octahedral_cell_area():
    """(latitudes, per-cell area in km^2) for the O320 octahedral grid."""
    nlat = 2 * N_OCT
    x, w = np.polynomial.legendre.leggauss(nlat)
    order = np.argsort(-x)
    x, w = x[order], w[order]
    lat = np.degrees(np.arcsin(x))
    i = np.arange(nlat)
    nlon = 20 + 4 * np.minimum(i, nlat - 1 - i)
    return lat, 2.0 * np.pi * R_EARTH**2 * w / nlon


def area_for(lat_vals, grid_lat, grid_area, tol=1e-3):
    idx = np.abs(grid_lat[None, :] - lat_vals[:, None]).argmin(axis=1)
    resid = np.abs(grid_lat[idx] - lat_vals)
    if resid.max() > tol:
        raise SystemExit(f"latitude {lat_vals[resid.argmax()]} is not an O320 row "
                         f"(residual {resid.max():.4f} deg)")
    return grid_area[idx]


def read_var(lpj, chunk, fname, column):
    path = os.path.join(lpj, chunk, "run1", fname)
    header = open(path).readline().split()
    if column not in header:
        raise SystemExit(f"{path}: no column {column!r} (have {header})")
    cols = [header.index(c) for c in ("Lat", "Year", column)]
    # usecols is positional and pandas returns them in file order, so name by
    # sorted position rather than by the order requested above
    names = [n for _, n in sorted(zip(cols, ("Lat", "Year", "frac")))]
    return pd.read_csv(path, sep=r"\s+", skiprows=1, header=None,
                       usecols=cols, names=names)


def series(lpj, chunks, fname, column, grid_lat, grid_area, nyears):
    """{year: (area Mkm2, land-mean fraction)} for one landcover column."""
    res = {}
    for chunk in chunks:
        if len(res) >= nyears:
            break
        if not os.path.exists(os.path.join(lpj, chunk, "run1", fname)):
            continue
        df = read_var(lpj, chunk, fname, column)
        for year, g in df.groupby("Year"):
            a = area_for(g["Lat"].to_numpy(), grid_lat, grid_area)
            f = g["frac"].to_numpy() / PERCENT
            res[int(year)] = ((f * a).sum() * 1e-6, (f * a).sum() / a.sum())
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=DEFAULT_RUN)
    ap.add_argument("--nyears", type=int, default=20)
    ap.add_argument("--frac", action="store_true",
                    help="plot the land-mean fraction instead of absolute area")
    ap.add_argument("--fpc", action="store_true",
                    help="plot the FPC-based cropFrac/pastureFrac instead of the "
                         "prescribed fracLut land-use tiles")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    lpj = os.path.join(args.run, "outdata", "lpj_guess")
    expid = os.path.basename(args.run.rstrip("/"))
    out = args.out or os.path.join(
        PLOTDIR, f"cropland_timeline_{expid}_first{args.nyears}y.png")

    grid_lat, grid_area = octahedral_cell_area()
    chunks = sorted(d for d in os.listdir(lpj) if re.fullmatch(r"\d{8}-\d{8}", d))

    # fracLut is the prescribed land-use tile fraction (CMIPoutput.cpp:1394).
    # cropFrac/pastureFrac are NOT land use: they accumulate indiv_fpc
    # (CMIPoutput.cpp:1114/1126), i.e. foliar cover, which responds to climate.
    if args.fpc:
        want = [("cropFrac_yearly.out", "Total", "cropland (foliar cover)", "#b3591a"),
                ("pastureFrac_yearly.out", "Total", "pasture (foliar cover)", "#1a5c7a")]
    else:
        want = [("fracLut_yearly.out", "crp", "cropland (land-use tile)", "#b3591a"),
                ("fracLut_yearly.out", "pst", "pasture (land-use tile)", "#1a5c7a"),
                ("fracLut_yearly.out", "psl", "natural (land-use tile)", "#2e7d32")]
    data = {}
    for fname, column, label, colour in want:
        s = series(lpj, chunks, fname, column, grid_lat, grid_area, args.nyears)
        if not s:
            print(f"!! {fname} not found in any chunk - skipping")
            continue
        years = np.array(sorted(s))[:args.nyears]
        data[label] = (years,
                       np.array([s[y][0] for y in years]),
                       np.array([s[y][1] for y in years]),
                       colour)

    if not data:
        raise SystemExit(f"no landcover output found under {lpj}")

    for label, (years, area, frac, _) in data.items():
        print(f"{label:9s} {years[0]}-{years[-1]}: "
              f"{area[0]:7.3f} -> {area[-1]:7.3f} M km2  "
              f"(change {area[-1] - area[0]:+.3f}, range {area.ptp():.3f}); "
              f"land-mean {frac[0]:.5f} -> {frac[-1]:.5f}")
        d = np.diff(area)
        print(f"{'':9s} years with any change: {np.count_nonzero(np.abs(d) > 1e-9)}"
              f" of {len(d)};  largest single-year step {np.abs(d).max():.4f} M km2")

    fig, axes = plt.subplots(len(data), 1, figsize=(10, 4.0 * len(data)),
                             squeeze=False)
    for ax, (label, (years, area, frac, colour)) in zip(axes[:, 0], data.items()):
        y = frac if args.frac else area
        ax.plot(years, y, "-o", ms=4, lw=1.6, color=colour)
        if args.frac:
            ax.set_ylabel(f"{label}\nfraction of land")
            note = f"full range {years[0]}-{years[-1]}: {y.ptp():.2e}"
        else:
            ax.set_ylabel(f"{label} area  (million km$^2$)")
            note = f"full range {years[0]}-{years[-1]}: {y.ptp():.3f} M km$^2$"
        ax.grid(alpha=0.3)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.text(0.985, 0.06, note, transform=ax.transAxes, ha="right",
                fontsize=8.5, color="0.3")
    axes[-1, 0].set_xlabel("year")
    span = next(iter(data.values()))[0]
    axes[0, 0].set_title("Global land-use evolution (LPJ-GUESS yearly landcover output)\n"
                         f"{expid}  {span[0]}-{span[-1]}  ({len(span)} years)",
                         fontsize=11, fontweight="bold")

    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
