#!/usr/bin/env python3
"""
Point-specific test of the "cold climate suppresses LPJG establishment" hypothesis.

At the cells where the prescribed/obs vegetation HAS forest but LPJG does NOT
(cvh_noLPJG > 0.3 AND cvh_LPJG < 0.1, boreal land), check whether the model's
GROWING-SEASON climate falls below the actual establishment gates:
  gate 1: warmest-month mean T2m >= twmin_est = 5 C
  gate 2: growing-degree-days above 5C  >= gdd5min_est = 500
Compare WITH-LPJG (cold), WITHOUT-LPJG (warm, prescribed forest) and ERA5.
If WITH is below the gate / below ERA5 while WITHOUT & ERA5 clear it, the
establishment failure is climate-driven (fix the atmosphere, not the limits).

Usage: python establishment_climate_check.py [Y0 Y1]   (default 1370 1379)
"""
import sys, glob, warnings, calendar
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1370, 1379)
years = list(range(Y0, Y1 + 1))
WO = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
WI = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs"
ERA5MON = "/work/ab0246/a270092/obs/era5/netcdf/T2M_mon.nc"

def files(run, var, freq):  # freq 'm' or 'd'
    base = WO if run == "WO" else WI
    out = []
    for y in years:
        if run == "WO":
            pat = f"atm_remapped_1{freq}_{var}_1{freq}_{y}-{y}.nc"
        else:
            pat = f"atm_remapped_1{freq}_{var}_{y}-{y}.nc"
        out += glob.glob(f"{base}/{pat}")
    return sorted(out)

def clim_monthly(run):  # 12-month climatology of T2m in Celsius, (12,lat,lon)
    ds = xr.open_mfdataset(files(run, "2t", "m"), combine="by_coords", use_cftime=True)
    da = ds["2t"] - 273.15
    td = "time_counter" if "time_counter" in da.dims else "time"
    return da.groupby(da[td].dt.month).mean(td).load()

def warmest_month(run):
    return clim_monthly(run).max("month")

def gdd5(run):  # mean annual growing-degree-days above 5C, from daily T2m
    tot = None; n = 0
    for y in years:
        fs = files_year(run, "2t", "d", y)
        if not fs: continue
        ds = xr.open_dataset(fs[0]); da = ds["2t"] - 273.15
        td = "time_counter" if "time_counter" in da.dims else "time"
        g = (da.where(da > 5) - 5).fillna(0).sum(td)   # sum over days of max(0,T-5)
        tot = g if tot is None else tot + g; n += 1
    return (tot / n).load()

def files_year(run, var, freq, y):
    base = WO if run == "WO" else WI
    pat = f"atm_remapped_1{freq}_{var}_1{freq}_{y}-{y}.nc" if run == "WO" else f"atm_remapped_1{freq}_{var}_{y}-{y}.nc"
    return glob.glob(f"{base}/{pat}")

def cvh(run):
    if run == "WO":
        fs = sorted(sum([glob.glob(f"{WO}/atm_remapped_1d_cvh_1d_{y}-{y}.nc") for y in years], []))
        v = "cvh"
    else:
        fs = sorted(sum([glob.glob(f"{WI}/atm_remapped_1m_cvh_{y}-{y}.nc") for y in years], []))
        v = "cvh"
    ds = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True)
    td = "time_counter" if "time_counter" in ds[v].dims else "time"
    return ds[v].mean(td).load()

print(f"Loading fields {Y0}-{Y1} ...")
cvh_wo, cvh_wi = cvh("WO"), cvh("WI")
lat = cvh_wo["lat"]; lon = cvh_wo["lon"]
cm_wo, cm_wi = clim_monthly("WO"), clim_monthly("WI")     # (12,lat,lon) C
def seas(cm, months): return cm.sel(month=months).mean("month")
tw_wo, tw_wi = seas(cm_wo, [6, 7, 8]),  seas(cm_wi, [6, 7, 8])    # JJA (~warmest month)
tc_wo, tc_wi = seas(cm_wo, [12, 1, 2]), seas(cm_wi, [12, 1, 2])   # DJF (coldest season)
gdd_wo, gdd_wi = gdd5("WO"), gdd5("WI")

# ERA5 seasonal-mean climatologies (trusted files; siblings of T2M_DJF used earlier), regridded
def era_seas(fname):
    d = xr.open_dataset(f"/work/ab0246/a270092/obs/era5/netcdf/{fname}", decode_times=False)
    v = [x for x in d.data_vars if d[x].ndim >= 2][0]; a = d[v].squeeze()
    a = a.rename({k: vv for k, vv in {"longitude": "lon", "latitude": "lat"}.items() if k in a.coords or k in a.dims})
    if float(a.max()) > 100: a = a - 273.15
    if float(a.lat[0]) > float(a.lat[-1]): a = a.isel(lat=slice(None, None, -1))
    return xr.concat([a, a.isel(lon=0).assign_coords(lon=float(a.lon[0]) + 360)], dim="lon").interp(lat=lat, lon=lon)
era_tw, era_tc = era_seas("T2M_JJA.nc"), era_seas("T2M_DJF.nc")

# land mask + boreal
lm = None
lmf = glob.glob(f"{WI}/atm_remapped_1m_lsm_{Y0}-{Y0}.nc")
if lmf:
    lm = xr.open_dataset(lmf[0])["lsm"]; lm = lm.isel(time_counter=0) if "time_counter" in lm.dims else lm
    lm = lm.reset_coords(drop=True)
boreal_land = (lat >= 50) & (lat <= 75) & ((lm > 0.5) if lm is not None else True)

# THE mask: obs/prescribed forest present, LPJG failed
should = boreal_land & (cvh_wo > 0.3) & (cvh_wi < 0.1)
# control: LPJG-forest actually present
ok = boreal_land & (cvh_wi > 0.3)

W = np.cos(np.deg2rad(lat)).broadcast_like(cvh_wo)
def wm(field, m):
    b = field.where(m); w = W.where(np.isfinite(b))
    return float(b.weighted(w.fillna(0)).mean().values)
def frac_below(field, thr, m):
    b = field.where(m); tot = W.where(np.isfinite(b)).sum()
    return float((W.where(np.isfinite(b) & (b < thr)).sum() / tot).values)
def ncells(m):
    return int(m.where(boreal_land).sum().values)

print(f"\n=== 'trees SHOULD be (cvh_noLPJG>0.3) but are NOT (cvh_LPJG<0.1)' : {ncells(should)} boreal-land cells ===")
print(f"(control: LPJG forest present cvh_LPJG>0.3 : {ncells(ok)} cells)\n")

print("SUMMER (JJA) mean T2m  [proxy for warmest-month gate twmin_est = +5 C]")
print(f"  at SHOULD-but-NOT cells:  WITH-LPJG={wm(tw_wi,should):6.2f}   WITHOUT={wm(tw_wo,should):6.2f}   ERA5={wm(era_tw,should):6.2f} C")
print(f"     WITH-LPJG summer cold bias vs ERA5 = {wm(tw_wi,should)-wm(era_tw,should):+.2f} C")
print(f"     % of these cells BELOW the +5C gate:  WITH-LPJG={100*frac_below(tw_wi,5,should):5.1f}%   WITHOUT={100*frac_below(tw_wo,5,should):5.1f}%   ERA5={100*frac_below(era_tw,5,should):5.1f}%")
print(f"  (control, LPJG-forest cells):  WITH-LPJG={wm(tw_wi,ok):6.2f}   ERA5={wm(era_tw,ok):6.2f} C")

print("\nGDD5 (growing-degree-days >5C, daily)  [gate gdd5min_est = 500]")
print(f"  at SHOULD-but-NOT cells:  WITH-LPJG={wm(gdd_wi,should):6.0f}   WITHOUT={wm(gdd_wo,should):6.0f}")
print(f"     % of these cells BELOW the 500 gate:  WITH-LPJG={100*frac_below(gdd_wi,500,should):5.1f}%   WITHOUT={100*frac_below(gdd_wo,500,should):5.1f}%")
print(f"  (control, LPJG-forest cells):  WITH-LPJG={wm(gdd_wi,ok):6.0f}")

def frac_above(field, thr, m):
    b = field.where(m); tot = W.where(np.isfinite(b)).sum()
    return float((W.where(np.isfinite(b) & (b > thr)).sum() / tot).values)
print("\nWINTER (DJF) mean T2m  [proxy for coldest-month gate tcmin_est: BNE/BINE = -30 C; BNS/larch = none]")
print(f"  at SHOULD-but-NOT cells:  WITHOUT(forested)={wm(tc_wo,should):6.2f}   WITH-LPJG={wm(tc_wi,should):6.2f}   ERA5={wm(era_tc,should):6.2f} C")
print(f"     >> LPJG-induced winter cooling (WITH - WITHOUT) = {wm(tc_wi,should)-wm(tc_wo,should):+.2f} C   [the self-reinforcing feedback]")
print(f"     >> WITHOUT (forested) bias vs ERA5 = {wm(tc_wo,should)-wm(era_tc,should):+.2f} C   [the 'initial'/trigger climate]")
print(f"     % of cells with coldest-month BELOW -30C (evergreen EXCLUDED):")
print(f"        WITHOUT(initial/forested)={100*frac_below(tc_wo,-30,should):5.1f}%   ERA5={100*frac_below(era_tc,-30,should):5.1f}%   -->  WITH-LPJG(degraded)={100*frac_below(tc_wi,-30,should):5.1f}%")
print(f"  (control, LPJG-forest cells):  WITHOUT={wm(tc_wo,ok):6.2f}   WITH-LPJG={wm(tc_wi,ok):6.2f}   ERA5={wm(era_tc,ok):6.2f} C")

# map (North-Polar-Stereo + coastlines so the boreal cells read as geography)
try:
    import cartopy.crs as ccrs, cartopy.feature as cfeature; HAVE = True
except Exception:
    HAVE = False
panels = [
    (tw_wi.where(should),            "JJA-mean 2 m T [$^\\circ$C] @ should-forest cells\n(gate twmin_est $=+5$)", -2, 14, "viridis"),
    ((tw_wi - era_tw).where(should), "JJA-mean 2 m T bias: WITH-LPJG $-$ ERA5 [$^\\circ$C]",                    -6,  6, "RdBu_r"),
    (gdd_wi.where(should),           "GDD5 @ should-forest cells\n(gate gdd5min_est $=500$)",                        0, 1200, "viridis")]
proj = ccrs.NorthPolarStereo() if HAVE else None
fig = plt.figure(figsize=(16, 5.6))
for i, (fld, ttl, vmn, vmx, cm) in enumerate(panels, 1):
    a = fig.add_subplot(1, 3, i, projection=proj) if HAVE else fig.add_subplot(1, 3, i)
    if HAVE:
        a.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
        a.add_feature(cfeature.LAND, facecolor="0.90", zorder=0)
        a.add_feature(cfeature.OCEAN, facecolor="#eef3f7", zorder=0)
        p = a.pcolormesh(lon, lat, fld, transform=ccrs.PlateCarree(), cmap=cm, vmin=vmn, vmax=vmx, shading="auto", zorder=2)
        a.coastlines(linewidth=0.5, zorder=3); a.gridlines(linewidth=0.2, color="0.5", alpha=0.5)
    else:
        p = a.pcolormesh(lon, lat, fld, cmap=cm, vmin=vmn, vmax=vmx, shading="auto"); a.set_ylim(45, 78)
    a.set_title(ttl, fontsize=9.5); fig.colorbar(p, ax=a, shrink=0.72, pad=0.03)
fig.suptitle(f"Growing-season climate where LPJG fails to forest ({Y0}-{Y1}) vs establishment gates "
             f"(only the {ncells(should)} 'should-forest' cells are coloured)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
OUT = "/work/bb1469/a270092/eval/plots/establishment_climate_check.png"
fig.savefig(OUT, dpi=130, bbox_inches="tight"); print("\nwrote", OUT)
