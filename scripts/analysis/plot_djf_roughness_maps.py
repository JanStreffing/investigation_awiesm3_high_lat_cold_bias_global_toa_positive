#!/usr/bin/env python3
"""
DJF roughness->inversion maps: WITH minus WITHOUT LPJG, boreal NH.
Shows the forest-cover collapse (cvh) spatially co-located with the colder
screen T, stronger low-level inversion, stronger surface decoupling, at
~unchanged cloud. Visual companion to djf_roughness_inversion.py (report sec:lpjg).

Usage: python plot_djf_roughness_maps.py [Y0 Y1]   (default 1350 1359)
Out:   /work/bb1469/a270092/eval/plots/djf_roughness_inversion_maps.png
"""
import sys, glob, warnings
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1350, 1359)
years = list(range(Y0, Y1 + 1))
OUT = "/work/bb1469/a270092/eval/plots/djf_roughness_inversion_maps.png"

WITHOUT = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WITH    = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"
VN = {"2t":"2t","skt":"skt","tcc":"tcc","strd":"strd","cvh":"cvh","pl_t":"t"}

def files(run, var, years):
    base = WITHOUT if run == "WITHOUT" else WITH
    out = []
    for y in years:
        if run == "WITHOUT":
            cands = ([f"atm_remapped_1m_pl_t_1m_pl_{y}-{y}.nc"] if var == "pl_t"
                     else [f"atm_remapped_1m_{var}_1m_{y}-{y}.nc",
                           f"atm_remapped_1d_{var}_1d_{y}-{y}.nc"])
        else:
            cands = [f"atm_remapped_1m_{var}_{y}-{y}.nc"]
        for c in cands:
            g = glob.glob(f"{base}/{c}")
            if g: out += g; break
    return sorted(out)

def djf_map(run, var, plev=None):
    fs = files(run, var, years)
    if not fs: return None
    ds = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True)
    da = ds[VN[var]]
    if plev is not None: da = da.sel(pressure_levels=plev)
    td = "time_counter" if "time_counter" in da.dims else "time"
    m = da[td].dt.month
    return da.where((m == 12) | (m == 1) | (m == 2), drop=True).mean(td).load()

print(f"Building DJF maps {Y0}-{Y1} ...")
D = {r: {v: djf_map(r, v, plev=(92500.0 if v == "pl_t" else None))
         for v in ("2t", "skt", "tcc", "strd", "cvh", "pl_t")} for r in ("WITHOUT", "WITH")}

def diff(v):  # WITH minus WITHOUT
    return (D["WITH"][v] - D["WITHOUT"][v])

lat = D["WITH"]["2t"]["lat"]; lon = D["WITH"]["2t"]["lon"]

# land mask (WITH lsm; same grid applies to both) -> focus on boreal land
lm = None
lmf = glob.glob(f"{WITH}/atm_remapped_1m_lsm_{Y0}-{Y0}.nc")
if lmf:
    lm = xr.open_dataset(lmf[0])["lsm"]
    if "time_counter" in lm.dims: lm = lm.isel(time_counter=0)
    lm = lm.reset_coords(drop=True)
def L(da):  # apply land mask if available
    return da.where(lm > 0.5) if lm is not None else da

def strd_wm2(da):  # detect accumulated J/m2 vs W/m2, return W/m2
    return da / 3600.0 if abs(float(da.mean().values)) > 1000 else da

dcvh = L(diff("cvh"))
dt2  = L(diff("2t"))
dinv = L((D["WITH"]["pl_t"] - D["WITH"]["2t"]) - (D["WITHOUT"]["pl_t"] - D["WITHOUT"]["2t"]))
dsfc = L((D["WITH"]["2t"] - D["WITH"]["skt"]) - (D["WITHOUT"]["2t"] - D["WITHOUT"]["skt"]))
dcld = L(diff("tcc"))
dlw  = L(strd_wm2(D["WITH"]["strd"]) - strd_wm2(D["WITHOUT"]["strd"]))

panels = [
    (dcvh, "Δ high-veg cover cvh", "BrBG",   0.6,  ""),
    (dt2,  "Δ 2 m T",             "RdBu_r", 6.0,  "K"),
    (dinv, "Δ low-level inversion (T925−T2m)", "RdBu_r", 3.0, "K"),
    (dsfc, "Δ surface inversion (T2m−skt)",    "RdBu_r", 2.0, "K"),
    (dcld, "Δ cloud cover",       "RdBu_r", 0.15, ""),
    (dlw,  "Δ LW down (strd)",    "RdBu_r", 8.0, "W/m²"),
]

fig = plt.figure(figsize=(15, 9))
proj = ccrs.NorthPolarStereo() if HAVE_CARTOPY else None
for i, (da, title, cmap, vlim, unit) in enumerate(panels, 1):
    ax = fig.add_subplot(2, 3, i, projection=proj) if HAVE_CARTOPY else fig.add_subplot(2, 3, i)
    if HAVE_CARTOPY:
        ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
        pc = ax.pcolormesh(lon, lat, da, transform=ccrs.PlateCarree(),
                           cmap=cmap, vmin=-vlim, vmax=vlim, shading="auto")
        ax.coastlines(linewidth=0.4); ax.add_feature(cfeature.LAND, facecolor="none")
        ax.gridlines(linewidth=0.2)
    else:
        sub = da.sel(lat=slice(45, 90))
        pc = ax.pcolormesh(sub["lon"], sub["lat"], sub, cmap=cmap, vmin=-vlim, vmax=vlim)
        ax.set_ylim(45, 90)
    ax.set_title(title + (f" [{unit}]" if unit else ""), fontsize=11)
    fig.colorbar(pc, ax=ax, shrink=0.7, pad=0.03)

fig.suptitle(f"DJF WITH − WITHOUT LPJG, boreal NH ({Y0}–{Y1}): forest loss → colder, "
             f"stronger inversion, decoupled surface — cloud ~unchanged", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print("wrote", OUT)
