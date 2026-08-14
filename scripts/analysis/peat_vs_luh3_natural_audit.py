#!/usr/bin/env python3
"""How many cells does the 46e7a84 peatland-vs-LUH3-NATURAL check put at risk?

Laszlo's 46e7a84 changed how the prescribed peat map and LUH3 are reconciled.

  OLD (<= 168a98b, framework/externalinput.cpp:723):
      lc.frac[PEATLAND] = peat_from_file          # peat map wins outright
      lc.frac[NATURAL]  = max(0, 1 - sum(others)) # NATURAL absorbs the remainder

  NEW (46e7a84, framework/externalinput.cpp:888-925):
      lc.frac[NATURAL] = luh3.get_natural(year)   # LUH3 is authoritative
      peat is CARVED FROM NATURAL; on a restart carrying a physical peat stand,
      a stand larger than the available NATURAL is a hard fail() rather than a
      silent truncation, because truncating would discard C/N/water state with
      no transfer.

So the new code does not invent a constraint -- it inverts the priority, and
then refuses to lose carbon quietly.  This script scopes how much of our
TCO95/CORE3 configuration sits on the wrong side of that inversion.

    peat map value  >  LUH3 NATURAL(1850) = primf + primn + secdf + secdn

is a necessary condition for the abort.  It is an OVER-estimate of the true
blast radius: the check compares the *physical* peat area carried in the
restart, which the spin-up may already have capped below the map value.  Cells
flagged here are at risk; the ones that actually fired are a subset.

Usage:
    python peat_vs_luh3_natural_audit.py [--year 1850]
"""

import argparse
import sys

import numpy as np
from netCDF4 import Dataset

PEAT_TXT = ("/work/bb1469/a270092/runtime/awiesm3-v3.4/Tuning_test_11I_v2soil/"
            "run_13500101-13591231/work/landuse/peat_frac.txt")
LUH3_STATES = ("/work/ab0246/a270092/input/lpj-guess/land_use/TCO95/"
               "multiple-states_input4MIPs_landState_CMIP_UofMD-landState-3-1-2_"
               "gn_0850-2024_TCO95.nc")
GRIDLIST = ("/work/bb1469/a270092/runtime/awiesm3-v3.4/Tuning_test_11I_v2soil/"
            "run_13500101-13591231/work/ece_gridlist_TCO95.txt")

# framework/externalinput.cpp: the fail() fires on
#   represented_peatland > lc.frac[NATURAL] + FRACTION_IDENTITY_TOLERANCE
TOL = 1.0e-6

# The cell 46e7a84 aborted on, for a self-check that we are reading the same data
# the model read.
KNOWN_FAIL = dict(lon=19.285715, lat=67.791733,
                  physical_peat=0.19172630653442613,
                  luh3_natural=0.088604211807250977)


def read_peat(path):
    """peat_frac.txt is 'Lon Lat Peat_Frac', one header line."""
    lon, lat, frac = [], [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            lon.append(float(parts[0]))
            lat.append(float(parts[1]))
            frac.append(float(parts[2]))
    return np.array(lon), np.array(lat), np.array(frac)


CROP = ("c3ann", "c4ann", "c3per", "c4per", "c3nfx")
PAST = ("pastr", "range")
NAT_RAW = ("primf", "primn", "secdf", "secdn")


def luh3_natural(path, year, run_barren):
    """The NATURAL fraction the peatland check actually compares against.

    close_luh3_base_fractions() (externalinput.cpp:54) makes NATURAL the single
    elastic category: NATURAL_closed = 1 - cropland - pasture - urban - barren.

    With run_barren 1, barren = 1 - sum(all 12 state vars at 1850)
    (luh3input.cpp:1031), and the algebra collapses to NATURAL_closed = the raw
    primf+primn+secdf+secdn.  With run_barren 0 the cache sets barren = 0
    (luh3input.cpp:958) and NATURAL additionally absorbs the LUH3 closure gap.

    OUR configuration is run_barren 0 -- the staged guess.ins overrides
    landcover.ins -- so the second branch is the live one.  Reproducing the
    reported 0.088604 at the abort cell rather than the raw 0.043008 is the
    self-check that this is right.
    """
    with Dataset(path) as ds:
        time = ds.variables["time"][:]
        units = ds.variables["time"].units
        base_year = int(units.split("since")[1].strip()[:4])
        years = base_year + np.round(np.asarray(time) / 365.0).astype(int)
        idx = int(np.argmin(np.abs(years - year)))
        if years[idx] != year:
            print(f"  note: nearest available year is {years[idx]}, not {year}")

        # Cells outside the LUH3 domain carry a ~2e20 fill value.  The model does
        # not arithmetic on those: load_state_data sets natural_fallback and gives
        # the cell NATURAL = 1 with zero transitions (luh3input.cpp:666-672), so a
        # peat stand can never exceed NATURAL there.  Mask them out the same way
        # rather than letting 2e20 flow into the closure.
        def get(v):
            x = np.asarray(ds.variables[v][idx, :], dtype=np.float64)
            return np.where(np.isfinite(x) & (x >= 0.0) & (x <= 1.0), x, np.nan)

        def s(names):
            return sum(get(v) for v in names)

        crop = s(CROP)
        past = s(PAST)
        urban = get("urban")
        nat_raw = s(NAT_RAW)

        fallback = ~np.isfinite(crop + past + urban + nat_raw)
        if run_barren:
            nat = nat_raw
        else:
            nat = np.maximum(0.0, 1.0 - crop - past - urban)
        # NATURAL=1 fallback wherever LUH3 has no usable state
        nat = np.where(fallback, 1.0, nat)
        past = np.where(fallback, np.nan, past)
        print(f"  {int(fallback.sum())} cells fall back to NATURAL=1 "
              f"(no usable LUH3 state)")

        return (np.asarray(ds.variables["lon"][:], dtype=np.float64),
                np.asarray(ds.variables["lat"][:], dtype=np.float64),
                nat, nat_raw, past, int(years[idx]))


def match_nearest(plon, plat, glon, glat):
    """Peat text uses -180..180, LUH3 uses 0..360; reconcile before matching."""
    p = np.where(plon < 0.0, plon + 360.0, plon)
    out = np.full(p.size, -1, dtype=int)
    for i in range(p.size):
        d = (glon - p[i]) ** 2 + (glat - plat[i]) ** 2
        j = int(np.argmin(d))
        if d[j] < 1e-4:          # ~0.01 deg
            out[i] = j
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, default=1850,
                    help="LUH3 year the run is pinned to (fixed_LU)")
    ap.add_argument("--run-barren", type=int, default=0,
                    help="run_barren from the staged guess.ins (ours is 0)")
    args = ap.parse_args()

    plon, plat, peat = read_peat(PEAT_TXT)
    print(f"peat map:   {peat.size} cells, "
          f"{np.count_nonzero(peat > 0)} with peat > 0, max {peat.max():.4f}")

    glon, glat, nat, nat_raw, past, used_year = luh3_natural(
        LUH3_STATES, args.year, args.run_barren)
    print(f"LUH3 {used_year}: {nat.size} cells, NATURAL mean {nat.mean():.4f} "
          f"(run_barren={args.run_barren})")

    j = match_nearest(plon, plat, glon, glat)
    ok = j >= 0
    if not ok.all():
        print(f"  WARNING: {np.count_nonzero(~ok)} peat cells found no LUH3 match")

    p = peat[ok]
    n = nat[j[ok]]
    rangefrac = past[j[ok]]
    lo = plon[ok]
    la = plat[ok]

    # self-check against the cell the model actually aborted on
    d = (lo - KNOWN_FAIL["lon"]) ** 2 + (la - KNOWN_FAIL["lat"]) ** 2
    k = int(np.argmin(d))
    print(f"\nself-check at the cell that aborted "
          f"({KNOWN_FAIL['lat']:.4f}N {KNOWN_FAIL['lon']:.4f}E):")
    print(f"  peat map        {p[k]:.6f}")
    print(f"  LUH3 NATURAL    {n[k]:.6f}   (model reported "
          f"{KNOWN_FAIL['luh3_natural']:.6f})")
    print(f"  physical peat in the restart state "
          f"{KNOWN_FAIL['physical_peat']:.6f}  <- what the check compares")
    if abs(n[k] - KNOWN_FAIL["luh3_natural"]) > 1e-5:
        print("  SELF-CHECK FAILED -- this script is not reading what the model "
              "read; the counts below are meaningless.")
    else:
        print("  self-check OK: reproduces the model's NATURAL exactly.")

    at_risk = p > n + TOL
    print(f"\ncells where the peat map exceeds LUH3 NATURAL({used_year}): "
          f"{np.count_nonzero(at_risk)} of {p.size} "
          f"({100.0 * np.count_nonzero(at_risk) / p.size:.2f} %)")
    if at_risk.any():
        excess = (p - n)[at_risk]
        print(f"  excess peat fraction: mean {excess.mean():.4f}, "
              f"max {excess.max():.4f}")
        print(f"  latitude span: {la[at_risk].min():.2f} to "
              f"{la[at_risk].max():.2f}")
        nh = np.count_nonzero(la[at_risk] > 50)
        print(f"  {nh} of them north of 50N "
              f"({100.0 * nh / np.count_nonzero(at_risk):.0f} %)")
        order = np.argsort(-excess)
        print("\n  worst 10 cells (lat, lon, peat, LUH3 natural, excess, "
              "pastr+range):")
        idx = np.where(at_risk)[0][order[:10]]
        for i in idx:
            print(f"    {la[i]:8.3f} {lo[i]:9.3f}  {p[i]:.4f}  {n[i]:.4f}  "
                  f"{p[i] - n[i]:.4f}  {rangefrac[i]:.4f}")

        # Is rangeland the mechanism?  LUH3 labels extensive high-latitude land
        # 'range', which LPJ-GUESS maps wholesale to PASTURE, squeezing NATURAL.
        print(f"\n  pastr+range on the at-risk cells: "
              f"mean {np.nanmean(rangefrac[at_risk]):.4f}, "
              f"median {np.nanmedian(rangefrac[at_risk]):.4f}")
        print(f"  pastr+range on the rest:          "
              f"mean {np.nanmean(rangefrac[~at_risk]):.4f}, "
              f"median {np.nanmedian(rangefrac[~at_risk]):.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
