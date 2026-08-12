#!/usr/bin/env python3
"""
Plot the CO2 concentration seen by LPJ-GUESS in the 1pctCO2 run.

The only CO2 concentration in the LPJ-GUESS output stream is the co2 column of
ifs_input.out (the rest are fluxes).  Each chunk is ~13 GB / 78 M lines, so a
full read is out of the question -- but the file is fixed-width and CO2 is
spatially uniform, so we seek straight to the first row of each day.

Layout, verified per chunk before use:
  * every line, header included, is exactly RECLEN bytes
  * rows are ordered (year, day, gridcell) with a constant cell count per day
  => byte offset of day d's first row = (1 + (d-1)*ncells) * RECLEN

Uniformity in space is not assumed; it is checked on sample days and the script
refuses to plot if CO2 varies across gridcells.
"""

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_RUN = ("/work/bb1469/a270089/runtime/awiesm3-v3.4.2/"
               "AWI-ESM3-VEG-HR-CMIP7-1pctCO2_1949")
PLOTDIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "plots")

RECLEN = 170       # bytes per line, header included
DAYS_PER_YEAR = 365
COL_LON, COL_LAT, COL_YEAR, COL_DAY, COL_CO2 = 0, 1, 2, 3, 4


def read_row(fh, index):
    """Parse data row `index` (0-based, excluding the header)."""
    fh.seek((1 + index) * RECLEN)
    return fh.read(RECLEN).split()


def chunk_layout(path):
    """(ndata, ncells, nyears) after verifying the fixed-width assumption."""
    size = os.path.getsize(path)
    if size % RECLEN:
        raise SystemExit(f"{path}: size {size} is not a multiple of {RECLEN} "
                         f"- not fixed-width, cannot seek")
    ndata = size // RECLEN - 1
    with open(path, "rb") as fh:
        header = fh.read(RECLEN).split()
        if header[COL_CO2] != b"co2":
            raise SystemExit(f"{path}: column {COL_CO2} is {header[COL_CO2]}, not co2")
        first = read_row(fh, 0)
        last = read_row(fh, ndata - 1)
    nyears = int(last[COL_YEAR]) - int(first[COL_YEAR]) + 1
    ndays = nyears * DAYS_PER_YEAR
    if ndata % ndays:
        raise SystemExit(f"{path}: {ndata} rows / {ndays} days is not an integer "
                         f"cell count")
    return ndata, ndata // ndays, nyears, int(first[COL_YEAR])


def chunk_co2(path, check_uniform=True):
    """(year, day, co2) for every day in one chunk, read by seeking."""
    ndata, ncells, nyears, year0 = chunk_layout(path)
    out = []
    with open(path, "rb") as fh:
        for d in range(nyears * DAYS_PER_YEAR):
            r = read_row(fh, d * ncells)
            out.append((int(r[COL_YEAR]), int(r[COL_DAY]), float(r[COL_CO2])))
            if check_uniform and d % 180 == 0:
                # same day, three widely separated gridcells
                for off in (ncells // 3, 2 * ncells // 3, ncells - 1):
                    s = read_row(fh, d * ncells + off)
                    if int(s[COL_DAY]) != int(r[COL_DAY]):
                        raise SystemExit(f"{path}: row ordering is not "
                                         f"(year, day, cell) as assumed")
                    if abs(float(s[COL_CO2]) - float(r[COL_CO2])) > 1e-6:
                        raise SystemExit(f"{path}: CO2 varies in space on day "
                                         f"{r[COL_DAY]} - cannot use one row per day")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default=DEFAULT_RUN)
    ap.add_argument("--nyears", type=int, default=20)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    lpj = os.path.join(args.run, "outdata", "lpj_guess")
    expid = os.path.basename(args.run.rstrip("/"))
    out = args.out or os.path.join(
        PLOTDIR, f"lpjg_co2_{expid}_first{args.nyears}y.png")

    chunks = sorted(d for d in os.listdir(lpj)
                    if re.fullmatch(r"\d{8}-\d{8}", d)
                    and os.path.exists(os.path.join(lpj, d, "run1", "ifs_input.out")))
    if not chunks:
        raise SystemExit(f"no ifs_input.out under {lpj}")

    rows = []
    for chunk in chunks:
        path = os.path.join(lpj, chunk, "run1", "ifs_input.out")
        ndata, ncells, nyears, year0 = chunk_layout(path)
        print(f"{chunk}: {ndata:,} rows, {ncells:,} cells/day, "
              f"{nyears} yr from {year0}")
        rows += chunk_co2(path)
        if len({r[0] for r in rows}) >= args.nyears:
            break

    rows.sort(key=lambda r: (r[0], r[1]))
    years = np.array([r[0] for r in rows])
    days = np.array([r[1] for r in rows])
    co2 = np.array([r[2] for r in rows])

    keep = years < years.min() + args.nyears
    years, days, co2 = years[keep], days[keep], co2[keep]
    t = years + (days - 1) / DAYS_PER_YEAR

    uy = np.unique(years)
    annual = np.array([co2[years == y].mean() for y in uy])
    y0 = uy[0]
    ideal = annual[0] * 1.01 ** (uy - y0)

    print(f"\nyears {uy[0]}-{uy[-1]}")
    print(f"CO2 first year mean : {annual[0]:.3f} ppm")
    print(f"CO2 last  year mean : {annual[-1]:.3f} ppm")
    print(f"ratio last/first    : {annual[-1] / annual[0]:.5f} "
          f"(1.01^{len(uy) - 1} = {1.01 ** (len(uy) - 1):.5f})")
    print(f"implied rate        : "
          f"{100 * ((annual[-1] / annual[0]) ** (1 / (len(uy) - 1)) - 1):.4f} %/yr")
    print(f"max |annual - 1%/yr|: {np.abs(annual - ideal).max():.3f} ppm")

    fig, axes = plt.subplots(2, 1, figsize=(10, 7.5),
                             gridspec_kw={"height_ratios": [2, 1]})

    ax = axes[0]
    ax.plot(t, co2, lw=0.8, color="#c0c0c0", label="daily (ifs_input.out)")
    ax.plot(uy + 0.5, annual, "-o", ms=4, lw=1.6, color="#b3591a",
            label="annual mean")
    ax.plot(uy + 0.5, ideal, "--", lw=1.3, color="#1a5c7a",
            label=f"{annual[0]:.1f} ppm $\\times$ 1.01$^{{y-{y0}}}$")
    ax.set_ylabel("CO$_2$  (ppm)")
    ax.set_title("CO$_2$ seen by LPJ-GUESS (ifs_input.out)\n"
                 f"{expid}  first {args.nyears} years",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(uy + 0.5, annual - ideal, "-o", ms=4, lw=1.4, color="#7a1a4a")
    ax.axhline(0, lw=1, color="0.5")
    ax.set_ylabel("annual $-$ 1%/yr  (ppm)")
    ax.set_xlabel("year")
    ax.grid(alpha=0.3)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
