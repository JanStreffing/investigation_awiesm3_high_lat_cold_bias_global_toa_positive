#!/usr/bin/env python3
"""
Global-mean peatland stand fraction over the coupled TCO319 historical run.

Source: fpc.out, column PEAT_STANDFRAC = gridcell.landcover.frac[PEATLAND]
        (commonoutput.cpp:1750/1797), a 0-1 fraction of the gridcell.

Two means are reported because they answer different questions:

  * land mean      -- area-weighted mean over the land cells LPJ-GUESS ran,
                      i.e. "what fraction of the land surface is peatland"
  * absolute area  -- the same weighted by true cell area, in million km^2

Cell areas come from the exact TCO319 octahedral grid: the Gauss-Legendre
weight w_i of latitude row i is precisely the width of that row's band in
sin(latitude), so cell area = 2*pi*R^2 * w_i / n_i with n_i the row's
longitude count (20 + 4*(i-1), mirrored about the equator).
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
               "AWI-ESM3-VEG-HR-CMIP7-historical_1949")
PEAT_IN = "/work/ab0246/a270092/input/lpj-guess/peat/TCO319_peat_frac.txt"
PLOTDIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "plots")

R_EARTH = 6371.0088  # km
# The Tco319 spectral truncation is carried on the O320 octahedral grid:
# 640 Gaussian latitude rows, 20 + 4*(i-1) longitudes on row i, 421,120 points.
# Verified against the run's own latitudes to < 1e-5 deg; N=319 does not fit.
N_OCT = 320


def octahedral_cell_area():
    """(latitudes, per-cell area in km^2) for the O320 octahedral grid."""
    nlat = 2 * N_OCT
    x, w = np.polynomial.legendre.leggauss(nlat)
    order = np.argsort(-x)                      # north to south
    x, w = x[order], w[order]
    lat = np.degrees(np.arcsin(x))
    i = np.arange(nlat)
    row = np.minimum(i, nlat - 1 - i)           # mirror about the equator
    nlon = 20 + 4 * row
    area = 2.0 * np.pi * R_EARTH**2 * w / nlon  # w sums to 2 over the sphere
    return lat, area


def area_for(lat_vals, grid_lat, grid_area, tol=1e-3):
    """Map each cell latitude onto its octahedral row and return its area.

    tol must be relaxed to ~0.005 for the input peat file, whose latitudes are
    written with only two decimals.
    """
    idx = np.abs(grid_lat[None, :] - lat_vals[:, None]).argmin(axis=1)
    resid = np.abs(grid_lat[idx] - lat_vals)
    if resid.max() > tol:
        raise SystemExit(f"latitude {lat_vals[resid.argmax()]} is not an O320 row "
                         f"(residual {resid.max():.4f} deg)")
    return grid_area[idx]


def read_chunk(lpj, chunk):
    """(year, lat, PEAT_STANDFRAC) rows from one output chunk's fpc.out."""
    path = os.path.join(lpj, chunk, "run1", "fpc.out")
    header = open(path).readline().split()
    cols = [header.index(c) for c in ("Lat", "Year", "PEAT_STANDFRAC")]
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None,
                     usecols=cols, names=["Lat", "Year", "frac"])
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=DEFAULT_RUN, help="experiment runtime directory")
    ap.add_argument("--nyears", type=int, default=None,
                    help="use only the first N simulated years")
    ap.add_argument("--out", default=None, help="output png path")
    args = ap.parse_args()

    lpj = os.path.join(args.run, "outdata", "lpj_guess")
    expid = os.path.basename(args.run.rstrip("/"))
    out = args.out or os.path.join(
        PLOTDIR, "peat_standfrac_timeline_" + expid
        + (f"_first{args.nyears}y" if args.nyears else "") + ".png")

    grid_lat, grid_area = octahedral_cell_area()

    chunks = sorted(d for d in os.listdir(lpj)
                    if re.fullmatch(r"\d{8}-\d{8}", d)
                    and os.path.exists(os.path.join(lpj, d, "run1", "fpc.out")))
    if not chunks:
        raise SystemExit(f"no fpc.out under {lpj}")

    years, land_mean, peat_area, land_area_ts = [], [], [], []
    for chunk in chunks:
        if args.nyears and len(years) >= args.nyears:
            break
        df = read_chunk(lpj, chunk)
        for year, g in df.groupby("Year"):
            a = area_for(g["Lat"].to_numpy(), grid_lat, grid_area)
            f = g["frac"].to_numpy()
            years.append(int(year))
            land_mean.append(float((f * a).sum() / a.sum()))
            peat_area.append(float((f * a).sum() * 1e-6))   # million km^2
            land_area_ts.append(float(a.sum() * 1e-6))

    o = np.argsort(years)
    years = np.array(years)[o]
    land_mean = np.array(land_mean)[o]
    peat_area = np.array(peat_area)[o]
    land_area = np.array(land_area_ts)[o]

    if len(np.unique(years)) != len(years):
        raise SystemExit("duplicate years across chunks - overlapping output?")

    if args.nyears:
        k = args.nyears
        years, land_mean = years[:k], land_mean[:k]
        peat_area, land_area = peat_area[:k], land_area[:k]

    # what the raw forcing file implies, for reference
    pin = np.loadtxt(PEAT_IN, skiprows=1)
    a_in = area_for(pin[:, 1], grid_lat, grid_area, tol=6e-3)
    area_in = (pin[:, 2] * a_in).sum() * 1e-6
    mean_in = (pin[:, 2] * a_in).sum() / a_in.sum()

    print(f"years {years[0]}-{years[-1]}  ({len(years)} years)")
    print(f"LPJ-GUESS land area          : {land_area.mean():8.3f} M km2")
    print(f"peatland area  first / last  : {peat_area[0]:8.4f} / {peat_area[-1]:8.4f} M km2 "
          f"(drift {peat_area[-1] - peat_area[0]:+.2e})")
    print(f"land-mean frac first / last  : {land_mean[0]:8.6f} / {land_mean[-1]:8.6f} "
          f"(drift {land_mean[-1] - land_mean[0]:+.2e})")
    print(f"peak-to-peak land-mean frac  : {land_mean.ptp():.3e}")
    print(f"input forcing file           : {area_in:8.4f} M km2, land-mean {mean_in:.6f} "
          f"over {a_in.sum()*1e-6:.3f} M km2")

    fig, axes = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True)

    ax = axes[0]
    ax.plot(years, land_mean, "-o", ms=4, lw=1.6, color="#b3591a")
    ax.axhline(mean_in, ls="--", lw=1.2, color="0.45",
               label=f"input forcing file ({mean_in:.4f})")
    ax.set_ylabel("peatland fraction of land")
    ax.set_title("Global peatland stand fraction (PEAT_STANDFRAC, fpc.out)\n"
                 f"{expid}  (TCO319L137)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(years, peat_area, "-o", ms=4, lw=1.6, color="#1a5c7a")
    ax.axhline(area_in, ls="--", lw=1.2, color="0.45",
               label=f"input forcing file ({area_in:.3f} M km$^2$)")
    ax.set_ylabel("peatland area  (million km$^2$)")
    ax.set_xlabel("year")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.3)

    # a flat line on an auto-scaled axis is misleading; state the true span
    for ax, v, unit in ((axes[0], land_mean, ""), (axes[1], peat_area, " M km$^2$")):
        ax.text(0.985, 0.06,
                f"full range over {years[0]}-{years[-1]}: {v.ptp():.2e}{unit}",
                transform=ax.transAxes, ha="right", fontsize=8.5, color="0.3")

    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
