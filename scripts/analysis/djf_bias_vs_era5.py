#!/usr/bin/env python3
"""
DJF 2m-T bias vs ERA5 over boreal land: current WITH-LPJG run, WITH + lever,
and WITHOUT-LPJG run (full cvh restoration = lever ceiling).

Answers: does the ~+1 K roughness lever close the boreal DJF cold bias, or
just dent it? Brackets the answer between 0 lever (WITH) and full (WITHOUT).

Usage: python djf_bias_vs_era5.py [Y0 Y1] [lever_K]   (default 1350 1369 1.0)
Out:   /work/bb1469/a270092/eval/plots/djf_bias_vs_era5.png
"""
import sys, glob, warnings
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    import cartopy.crs as ccrs, cartopy.feature as cfeature; HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1350, 1369)
LEVER = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
years = list(range(Y0, Y1 + 1))
WITHOUT = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WITH    = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"
ERA5    = "/work/ab0246/a270092/obs/era5/netcdf/T2M_DJF.nc"

def djf_2t(run):
    fs = []
    for y in years:
        pat = f"atm_remapped_1m_2t_1m_{y}-{y}.nc" if run == "WITHOUT" else f"atm_remapped_1m_2t_{y}-{y}.nc"
        g = glob.glob(f"{(WITHOUT if run=='WITHOUT' else WITH)}/{pat}");  fs += g
    ds = xr.open_mfdataset(sorted(fs), combine="by_coords", use_cftime=True)
    da = ds["2t"]; td = "time_counter" if "time_counter" in da.dims else "time"
    m = da[td].dt.month
    return da.where((m == 12) | (m == 1) | (m == 2), drop=True).mean(td).load()  # K

wo = djf_2t("WITHOUT"); wi = djf_2t("WITH")
lat = wi["lat"]; lon = wi["lon"]

# ERA5 -> model grid (periodic lon)
e = xr.open_dataset(ERA5, decode_times=False)["tas"].squeeze()
e = xr.concat([e, e.isel(lon=0).assign_coords(lon=360.0)], dim="lon")  # wrap
era = e.interp(lat=lat, lon=lon).drop_vars([c for c in ("time",) if c in e.coords], errors="ignore")

# boreal land mask
lm = None
lmf = glob.glob(f"{WITH}/atm_remapped_1m_lsm_{Y0}-{Y0}.nc")
if lmf:
    lm = xr.open_dataset(lmf[0])["lsm"]
    if "time_counter" in lm.dims: lm = lm.isel(time_counter=0)
    lm = lm.reset_coords(drop=True)
boreal = (lat >= 50) & (lat <= 80)
mask = boreal & ((lm > 0.5) if lm is not None else True)
W = np.cos(np.deg2rad(lat)).broadcast_like(wi)

def stat(field):  # area-weighted mean bias + rmse over boreal land
    b = (field - era).where(mask)
    w = W.where(np.isfinite(b))
    mean = float(b.weighted(w.fillna(0)).mean().values)
    rmse = float(np.sqrt((b**2).weighted(w.fillna(0)).mean().values))
    return mean, rmse

BOXES = {"Siberia 55-70N,60-140E": ((55,70),(60,140)), "Pan-boreal 50-70N land": ((50,80),(0,360))}
def boxbias(field, box):
    (la0,la1),(lo0,lo1)=box
    m = mask & (lat>=la0)&(lat<=la1)&(lon>=lo0)&(lon<=lo1)
    b=(field-era).where(m); w=W.where(np.isfinite(b))
    return float(b.weighted(w.fillna(0)).mean().values)

cases = {"WITH-LPJG (current)": wi, f"WITH + {LEVER:.1f}K lever": wi + LEVER, "WITHOUT-LPJG (full restore)": wo}

print(f"=== DJF 2m-T bias vs ERA5, boreal land, years {Y0}-{Y1} ({len(years)} yr) ===\n")
print(f"{'case':30s} {'Siberia':>9s} {'Pan-boreal':>11s} {'boreal RMSE':>12s}")
for name, f in cases.items():
    sib = boxbias(f, BOXES["Siberia 55-70N,60-140E"])
    pan = boxbias(f, BOXES["Pan-boreal 50-70N land"])
    _, rmse = stat(f)
    print(f"{name:30s} {sib:+9.2f} {pan:+11.2f} {rmse:12.2f}")
print("\n(bias = model - ERA5; negative = too cold)")

# maps
n=len(cases); fig=plt.figure(figsize=(6*n,5)); proj=ccrs.NorthPolarStereo() if HAVE_CARTOPY else None
for i,(name,f) in enumerate(cases.items(),1):
    b=(f-era).where(mask)
    ax=fig.add_subplot(1,n,i,projection=proj) if HAVE_CARTOPY else fig.add_subplot(1,n,i)
    if HAVE_CARTOPY:
        ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
        pc=ax.pcolormesh(lon,lat,b,transform=ccrs.PlateCarree(),cmap="RdBu_r",vmin=-12,vmax=12,shading="auto")
        ax.coastlines(linewidth=0.4); ax.gridlines(linewidth=0.2)
    else:
        pc=ax.pcolormesh(lon,lat.sel(lat=slice(45,90)),b.sel(lat=slice(45,90)),cmap="RdBu_r",vmin=-12,vmax=12)
    ax.set_title(name,fontsize=11); fig.colorbar(pc,ax=ax,shrink=0.7,label="2m-T bias vs ERA5 [K]")
fig.suptitle(f"DJF 2m-T bias vs ERA5, boreal NH ({Y0}-{Y1}): does the roughness lever close the cold bias?",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95])
OUT="/work/bb1469/a270092/eval/plots/djf_bias_vs_era5.png"; fig.savefig(OUT,dpi=130,bbox_inches="tight"); print("wrote",OUT)
