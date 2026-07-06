#!/usr/bin/env python3
"""
Test for a veg-cooling -> sea-ice -> more-cooling positive feedback DRIFT.

Per-year (annual-mean) NH sea-ice area and Arctic T2m for both runs over
legs 1+2 (1350-1369). Spin-up drift is present in BOTH runs; the FEEDBACK
signature specific to the veg cooling is the WITH-minus-WITHOUT trend:
if the LPJG (veg-cooling) run grows ice / cools the Arctic faster than the
no-LPJG control, that is the amplifying feedback.

Usage: python seaice_t2m_trend.py [Y0 Y1]   (default 1350 1369)
Out:   /work/bb1469/a270092/eval/plots/seaice_t2m_trend.png
"""
import sys, glob, warnings
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1350, 1369)
years = np.arange(Y0, Y1 + 1)
WITHOUT = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WITH    = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"

def yr_file(run, var, y):
    pat = f"atm_remapped_1m_{var}_1m_{y}-{y}.nc" if run == "WITHOUT" else f"atm_remapped_1m_{var}_{y}-{y}.nc"
    g = glob.glob(f"{(WITHOUT if run=='WITHOUT' else WITH)}/{pat}")
    return g[0] if g else None

# grid cell area (km^2) from the model lat/lon
_probe = xr.open_dataset(yr_file("WITHOUT", "2t", Y0))
lat = _probe["lat"].values; lon = _probe["lon"].values
R = 6371.0
dlat = np.deg2rad(np.abs(np.gradient(lat))); dlon = np.deg2rad(np.abs(np.gradient(lon)))
area2d = (R**2) * np.outer(np.cos(np.deg2rad(lat)) * dlat, dlon)  # km^2, (lat,lon)
LAT2 = np.broadcast_to(lat[:, None], area2d.shape)

def annual(run, var, y):
    f = yr_file(run, var, y)
    d = xr.open_dataset(f, decode_times=False)[var]
    td = "time_counter" if "time_counter" in d.dims else "time"
    return d.mean(td).values  # (lat,lon) annual mean

def series(run):
    ice, arc, icw = [], [], []
    for y in years:
        ci = annual(run, "ci", y)       # sea-ice concentration 0..1
        t2 = annual(run, "2t", y) - 273.15
        nh = LAT2 > 0
        ice.append(np.nansum(np.where(nh, ci, 0) * area2d) / 1e6)          # 10^6 km^2 NH ice area
        arcm = LAT2 > 70
        arc.append(np.nansum(np.where(arcm, t2 * area2d, 0)) / np.nansum(np.where(arcm, area2d, 0)))  # >70N T2m
        # ci-weighted (ice-covered) T2m, NH
        w = np.where(nh, ci * area2d, 0)
        icw.append(np.nansum(w * t2) / np.nansum(w))
    return np.array(ice), np.array(arc), np.array(icw)

def trend(y, v):  # slope per decade + p-ish via stderr
    a, b = np.polyfit(y, v, 1)
    resid = v - (a * y + b); se = np.sqrt(np.sum(resid**2) / (len(y) - 2)) / (np.std(y) * np.sqrt(len(y)))
    return a * 10, se * 10  # per decade

WO = series("WITHOUT"); WI = series("WITH")
labels = ["NH sea-ice area [10^6 km^2]", "Arctic (>70N) T2m [C]", "ice-covered T2m [C]"]
print(f"=== Trends over {Y0}-{Y1} (per decade) ===\n")
print(f"{'metric':30s} {'WITHOUT':>18s} {'WITH':>18s} {'WITH-WITHOUT':>14s}")
for i, lab in enumerate(labels):
    two, twe = trend(years, WO[i]); tiw, tie = trend(years, WI[i])
    diff = tiw - two
    print(f"{lab:30s} {two:+7.3f}±{twe:5.3f}     {tiw:+7.3f}±{tie:5.3f}     {diff:+7.3f}")
print("\nmeans:")
for i, lab in enumerate(labels):
    print(f"  {lab:30s} WITHOUT={WO[i].mean():8.3f}   WITH={WI[i].mean():8.3f}   diff={WI[i].mean()-WO[i].mean():+.3f}")

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for i, lab in enumerate(labels):
    ax[i].plot(years, WO[i], "o-", color="C0", label="WITHOUT LPJG", ms=4)
    ax[i].plot(years, WI[i], "s-", color="C3", label="WITH LPJG", ms=4)
    for v, c in [(WO[i], "C0"), (WI[i], "C3")]:
        a, b = np.polyfit(years, v, 1); ax[i].plot(years, a*years+b, "--", color=c, lw=1)
    ax[i].set_title(lab); ax[i].set_xlabel("year"); ax[i].axvline(Y0+9.5, color="k", lw=.3, ls=":")
    ax[i].legend(fontsize=8)
fig.suptitle(f"Sea-ice & Arctic-T2m drift, legs 1-2 ({Y0}-{Y1}): is there a veg-cooling->ice feedback?", fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95])
OUT="/work/bb1469/a270092/eval/plots/seaice_t2m_trend.png"; fig.savefig(OUT, dpi=130, bbox_inches="tight"); print("\nwrote", OUT)
