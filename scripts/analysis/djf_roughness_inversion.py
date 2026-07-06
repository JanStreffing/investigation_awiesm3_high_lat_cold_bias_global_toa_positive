#!/usr/bin/env python3
"""
DJF roughness -> inversion confirmation: WITH vs WITHOUT LPJG.

Mechanism under test (report sec:lpjg): LPJG exports almost no boreal forest
cover (cvh~0.05 vs ~0.44 without LPJG). Forest z0m=2.0 m vs tundra/grass
0.034-0.10 m (susveg_mod.F90). In polar-night DJF a smoother surface -> weaker
turbulent coupling -> the stable boundary layer decouples and the surface
radiatively cools -> stronger low-level inversion -> colder screen T. This
script checks whether that atmospheric response actually shows up in the data,
now that we have a CLEAN same-ocean (CORE3) no-LPJG control.

WITH LPJG   = a270270 Tuning_test_06 Baseline (awiesm3)
WITHOUT LPJG= awicm3_noLPJG_CORE3_30y (this clean reference; same CORE3 ocean)

Usage: python djf_roughness_inversion.py [Y0 Y1]   (default 1350 1359 = leg 1)
Both runs are on the identical 192x400 remapped grid -> direct differencing.
"""
import sys, glob, warnings
import numpy as np, xarray as xr
warnings.filterwarnings('ignore')

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1350, 1359)

WITHOUT = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WITH    = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"

# per-run filename builders (two naming conventions)
def files(run, var, years):
    out = []
    for y in years:
        if run == "WITHOUT":
            stem = f"atm_remapped_1m_pl_t_1m_pl" if var == "pl_t" else f"atm_remapped_1m_{var}_1m"
        else:
            stem = f"atm_remapped_1m_{var}"
        f = f"{globals()[run]}/{stem}_{y}-{y}.nc"
        g = glob.glob(f)
        if g: out += g
    return sorted(out)

VARNAME = {"2t":"2t", "skt":"skt", "tcc":"tcc", "strd":"strd", "cvh":"cvh", "pl_t":"t"}

def load_djf(run, var, years, plev=None):
    fs = files(run, var, years)
    if not fs:
        return None
    ds = xr.open_mfdataset(fs, combine="by_coords", decode_times=True, use_cftime=True)
    da = ds[VARNAME[var]]
    if plev is not None:
        da = da.sel(pressure_levels=plev)
    tdim = "time_counter" if "time_counter" in da.dims else "time"
    mon = da[tdim].dt.month
    da = da.where((mon == 12) | (mon == 1) | (mon == 2), drop=True)
    return da.mean(tdim)  # DJF climatological mean map

# boreal boxes (lat S->N, lon 0..360)
BOXES = {
    "Siberia 55-70N,60-140E": dict(lat=(55, 70), lon=(60, 140)),
    "Pan-boreal 50-70N belt":  dict(lat=(50, 70), lon=(0, 360)),
}

def boxmean(da, box, landmask=None):
    lat = da['lat']; lon = da['lon']
    m = (lat >= box['lat'][0]) & (lat <= box['lat'][1]) & \
        (lon >= box['lon'][0]) & (lon <= box['lon'][1])
    sub = da.where(m)
    if landmask is not None:
        sub = sub.where(landmask > 0.5)
    w = np.cos(np.deg2rad(lat))
    return float(sub.weighted(w.broadcast_like(sub).fillna(0)).mean().values)

years = list(range(Y0, Y1 + 1))
print(f"=== DJF WITH vs WITHOUT LPJG, years {Y0}-{Y1} ===\n")

# land mask from WITH lsm (same grid applies to both); fall back to None
lm = None
lmf = glob.glob(f"{WITH}/atm_remapped_1m_lsm_{Y0}-{Y0}.nc")
if lmf:
    lm = xr.open_dataset(lmf[0])['lsm']
    if 'time_counter' in lm.dims: lm = lm.isel(time_counter=0)
    print("(land mask: WITH lsm)\n")

# ---- premise: cvh contrast ----
cvh_wo = load_djf("WITHOUT", "cvh", years)   # daily cvh? -> handled below
cvh_wi = load_djf("WITH",    "cvh", years)
# WITHOUT cvh is only daily; retry with daily builder if monthly missing
if cvh_wo is None:
    fs = sorted(sum([glob.glob(f"{WITHOUT}/atm_remapped_1d_cvh_1d_{y}-{y}.nc") for y in years], []))
    if fs:
        ds = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True)
        da = ds['cvh']; tdim = "time_counter" if "time_counter" in da.dims else "time"
        mon = da[tdim].dt.month
        cvh_wo = da.where((mon==12)|(mon==1)|(mon==2), drop=True).mean(tdim)

print("PREMISE  boreal high-veg cover cvh (DJF):")
for name, box in BOXES.items():
    a = boxmean(cvh_wo, box, lm) if cvh_wo is not None else float('nan')
    b = boxmean(cvh_wi, box, lm) if cvh_wi is not None else float('nan')
    print(f"  {name:26s}  WITHOUT={a:.3f}  WITH={b:.3f}   (expect ~0.44 vs ~0.05)")

# ---- fields ----
print("\nLoading DJF fields (2t, skt, T925, tcc, strd)...")
F = {}
for run in ("WITHOUT", "WITH"):
    F[run] = {
        "2t":   load_djf(run, "2t",  years),
        "skt":  load_djf(run, "skt", years),
        "t925": load_djf(run, "pl_t", years, plev=92500.0),
        "tcc":  load_djf(run, "tcc", years),
        "strd": load_djf(run, "strd", years),
    }

def strd_wm2(da):
    # strd may be accumulated J/m2 over the output interval; monthly-mean flux
    # is stored as W/m2 in these XIOS runs if values ~200-350. Detect & scale.
    v = abs(float(da.mean().values))
    return da if v < 1000 else da / 3600.0

print("\nRESULT  DJF boreal, WITHOUT vs WITH LPJG:")
hdr = f"  {'box':26s} {'run':8s} {'T2m':>7s} {'skt':>7s} {'T925':>7s} {'inv(T925-T2m)':>14s} {'sfcinv(T2m-skt)':>16s} {'cloud':>6s} {'LWdn':>7s}"
print(hdr)
for name, box in BOXES.items():
    for run in ("WITHOUT", "WITH"):
        f = F[run]
        t2  = boxmean(f["2t"],  box, lm) - 273.15
        sk  = boxmean(f["skt"], box, lm) - 273.15
        t9  = boxmean(f["t925"],box, lm) - 273.15
        cl  = boxmean(f["tcc"], box, lm)
        lw  = boxmean(strd_wm2(f["strd"]), box, lm)
        inv = t9 - t2
        sfc = t2 - sk
        print(f"  {name:26s} {run:8s} {t2:7.2f} {sk:7.2f} {t9:7.2f} {inv:14.2f} {sfc:16.2f} {cl:6.2f} {lw:7.1f}")
    print()
print("Interpretation: if WITH LPJG (low cvh, smooth) shows a LARGER low-level")
print("inversion (T925-T2m) and/or larger surface inversion (T2m-skt) with colder")
print("T2m, the roughness->decoupling->inversion mechanism is confirmed in data.")
