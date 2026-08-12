#!/usr/bin/env python3
"""
Plot the LPJ-GUESS peatland (wetland) fraction actually carried by the coupled
TCO319 historical run, first simulated year vs last available year.

Run:  AWI-ESM3-VEG-HR-CMIP7-historical_1949 (TCO319L137 + DARS2)
Var:  wetlandFrac_monthly.out  -- gridcell fraction occupied by the peatland
      stand.  Prescribed from file_peat (TCO319_peat_frac.txt) but rescaled
      each year by the LUH3 natural-land fraction, so it is not necessarily
      constant in time.

Also overlays the raw input peat file for reference, and reports how many of
the run's land points the input file fails to cover.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

RUN = "/work/bb1469/a270089/runtime/awiesm3-v3.4.2/AWI-ESM3-VEG-HR-CMIP7-historical_1949"
LPJ = os.path.join(RUN, "outdata", "lpj_guess")
PEAT_IN = "/work/ab0246/a270092/input/lpj-guess/peat/TCO319_peat_frac.txt"
GRIDLIST = os.path.join(RUN, "run_18800101-18811231", "work", "ece_gridlist_TCO319.txt")
OUT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "plots", "peat_hist_first_last.png",
)

MONTHS = slice(3, 15)  # columns Jan..Dec in wetlandFrac_monthly.out


def read_wetland(chunk, year):
    """Annual-mean wetland fraction for one year of one output chunk.

    wetlandFrac_monthly.out is written in PERCENT of gridcell; convert to a
    0-1 fraction so it is directly comparable with the input peat file.
    """
    path = os.path.join(LPJ, chunk, "run1", "wetlandFrac_monthly.out")
    a = np.loadtxt(path, skiprows=1)
    m = a[:, 2].astype(int) == year
    if not m.any():
        raise SystemExit(f"year {year} not in {path}")
    return a[m, 0], a[m, 1], a[m, MONTHS].mean(axis=1) / 100.0


def panel(ax, lon, lat, val, title, cmap, vmin, vmax, size, thresh=0.001):
    """Scatter only the cells that carry a signal, so the grey land shows through."""
    # cartopy 0.20 / matplotlib >=3.6 incompatibility in GeoAxes.autoscale_view
    ax._autoscaleXon = getattr(ax, "get_autoscalex_on", lambda: True)()
    ax._autoscaleYon = getattr(ax, "get_autoscaley_on", lambda: True)()
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor="0.90", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="#dceaf2", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=3)
    m = np.abs(val) > thresh
    sc = ax.scatter(lon[m], lat[m], c=val[m], s=size, cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree(), linewidths=0, zorder=2)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.text(0.01, 0.02, f"shown: {m.sum():,} of {len(val):,} cells",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.85))
    return sc


def main():
    chunks = sorted(d for d in os.listdir(LPJ)
                    if os.path.exists(os.path.join(LPJ, d, "run1", "wetlandFrac_monthly.out")))
    first_chunk, last_chunk = chunks[0], chunks[-1]
    first_year = int(first_chunk[:4])
    last_year = int(last_chunk[9:13])

    lon0, lat0, w0 = read_wetland(first_chunk, first_year)
    lon1, lat1, w1 = read_wetland(last_chunk, last_year)

    # cells are written in gridlist order and the gridlist is static, so the
    # two years share an index; assert it rather than assume it
    same = (lon0.shape == lon1.shape and np.allclose(lon0, lon1) and np.allclose(lat0, lat1))
    if same:
        dw = w1 - w0
    else:
        key0 = {(round(x, 3), round(y, 3)): i for i, (x, y) in enumerate(zip(lon0, lat0))}
        idx = np.array([key0.get((round(x, 3), round(y, 3)), -1)
                        for x, y in zip(lon1, lat1)])
        ok = idx >= 0
        lon1, lat1, w1 = lon1[ok], lat1[ok], w1[ok]
        dw = w1 - w0[idx[ok]]

    print(f"first: {first_year}  n={len(w0):7d}  mean={w0.mean():.5f}  "
          f"peat>0.01: {(w0 > 0.01).sum():6d}  max={w0.max():.3f}")
    print(f"last : {last_year}  n={len(w1):7d}  mean={w1.mean():.5f}  "
          f"peat>0.01: {(w1 > 0.01).sum():6d}  max={w1.max():.3f}")
    print(f"delta: mean={dw.mean():+.6f}  min={dw.min():+.4f}  max={dw.max():+.4f}  "
          f"cells changed >1e-6: {(np.abs(dw) > 1e-6).sum()}")

    # coverage of the input peat file against the run's actual land points
    pin = np.loadtxt(PEAT_IN, skiprows=1)
    pset = {(round(x, 2), round(y, 2)) for x, y in zip(pin[:, 0], pin[:, 1])}
    gl = np.loadtxt(GRIDLIST, usecols=(0, 1))
    gset = {(round(x, 2), round(y, 2)) for x, y in zip(gl[:, 0], gl[:, 1])}
    print(f"input peat file: {len(pset)} cells; run gridlist: {len(gset)} cells; "
          f"gridlist points absent from peat file: {len(gset - pset)}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 9),
                             subplot_kw={"projection": ccrs.Robinson()})
    s = 1.6
    vmax = 0.8
    sc0 = panel(axes[0, 0], pin[:, 0], pin[:, 1], pin[:, 2],
                f"Input forcing  TCO319_peat_frac.txt  ({len(pin):,} cells)",
                "YlOrBr", 0, vmax, s)
    panel(axes[0, 1], lon0, lat0, w0,
          f"Simulated peatland fraction  {first_year}", "YlOrBr", 0, vmax, s)
    panel(axes[1, 0], lon1, lat1, w1,
          f"Simulated peatland fraction  {last_year}", "YlOrBr", 0, vmax, s)
    lim = max(1e-6, np.abs(dw).max())
    sc1 = panel(axes[1, 1], lon1, lat1, dw,
                f"Difference  {last_year} - {first_year}   "
                f"(max |d| = {lim:.2e})", "RdBu_r", -lim, lim, 6, thresh=1e-6)

    cb0 = fig.colorbar(sc0, ax=axes[:, 0], orientation="horizontal",
                       fraction=0.04, pad=0.03, shrink=0.75)
    cb0.set_label("peatland fraction of gridcell")
    cb1 = fig.colorbar(sc1, ax=axes[:, 1], orientation="horizontal",
                       fraction=0.04, pad=0.03, shrink=0.75)
    cb1.set_label("fraction  /  difference")

    fig.suptitle("AWI-ESM3-VEG-HR-CMIP7-historical_1949 (TCO319L137, DARS2) - "
                 "LPJ-GUESS peatland fraction", fontsize=13, fontweight="bold")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=140, bbox_inches="tight")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
