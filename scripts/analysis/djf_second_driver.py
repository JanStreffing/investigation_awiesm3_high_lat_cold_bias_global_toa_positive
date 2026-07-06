#!/usr/bin/env python3
"""
Deep search for the SECOND (non-vegetation) boreal DJFM cold driver.

Diagnose the WITHOUT-LPJG run (full forest -> vegetation ruled out) against ERA5
over boreal land. Whatever bias survives here is the residual driver. Tests:
  1. Stable-BL over-decoupling  -> low-level inversion (T925-T2m) vs ERA5
  2. Circulation (Siberian High) -> MSLP vs ERA5
  3. Cloud / LW deficit          -> total cloud vs ERA5 (annual clim caveat)
Then spatially correlate the residual 2m-T cold bias with each candidate bias
to see which one explains where the model is too cold.

Season DJFM (matches ERA5 T_DJFM/T2M_DJFM/MSL_DJFM). Usage: python djf_second_driver.py [Y0 Y1]
Out: /work/bb1469/a270092/eval/plots/djf_second_driver.png
"""
import sys, glob, warnings
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    import cartopy.crs as ccrs; HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1360, 1369)
years = list(range(Y0, Y1 + 1))
RUN = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/oifs"
O = "/work/ab0246/a270092/obs/era5/netcdf"
MON = [12, 1, 2, 3]  # DJFM

def mdl(var, plev=None):
    fs = []
    for y in years:
        pat = f"atm_remapped_1m_pl_t_1m_pl_{y}-{y}.nc" if var == "pl_t" else f"atm_remapped_1m_{var}_1m_{y}-{y}.nc"
        fs += glob.glob(f"{RUN}/{pat}")
    ds = xr.open_mfdataset(sorted(fs), combine="by_coords", use_cftime=True)
    da = ds["t" if var == "pl_t" else var]
    if plev is not None: da = da.sel(pressure_levels=plev)
    td = "time_counter" if "time_counter" in da.dims else "time"
    return da.where(da[td].dt.month.isin(MON), drop=True).mean(td).load()

REF = mdl("2t"); LAT = REF["lat"]; LON = REF["lon"]

def regrid(da):
    da = da.rename({k: v for k, v in {"longitude": "lon", "latitude": "lat"}.items() if k in da.coords or k in da.dims})
    if float(da.lat[0]) > float(da.lat[-1]): da = da.isel(lat=slice(None, None, -1))
    da = xr.concat([da, da.isel(lon=0).assign_coords(lon=float(da.lon[0]) + 360.0)], dim="lon")
    return da.interp(lat=LAT, lon=LON)

def era(fname, var, plev=None, months=None):
    d = xr.open_dataset(f"{O}/{fname}", decode_times=False)
    da = d[var]
    if plev is not None:
        pc = [c for c in da.coords if "plev" in c.lower() or "lev" in c.lower()][0]
        da = da.sel({pc: plev})
    for sq in ("time", "valid_time"):
        if sq in da.dims and da.sizes.get(sq, 1) == 1: da = da.squeeze(sq)
    if months is not None and ("time" in da.dims or "valid_time" in da.dims):
        td = "time" if "time" in da.dims else "valid_time"
        da = da.isel({td: [i for i in range(da.sizes[td]) if (i % 12 + 1) in months]}).mean(td)
    return regrid(da)

# ---- model DJFM fields ----
M = dict(
    t2m = mdl("2t"), skt = mdl("skt"), t925 = mdl("pl_t", 92500.0), t1000 = mdl("pl_t", 100000.0),
    msl = mdl("msl"), tcc = mdl("tcc"), strr = mdl("str"),
)
# ---- ERA5 DJFM fields on model grid ----
E = dict(
    t2m  = era("T2M_DJFM.nc", "T2M"),
    t925 = era("T_DJFM.nc", "T", 92500.0), t1000 = era("T_DJFM.nc", "T", 100000.0),
    msl  = era("MSL_DJFM.nc", "MSL"), tcc = era("timmean_tcc.nc", "tcc"),
)
for k in ("t2m","t925","t1000","skt"):
    if k in E:  # to Celsius-consistent K already; pressure to hPa below
        pass

# boreal land mask
lm = None
lmf = glob.glob("/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/oifs/atm_remapped_1m_lsm_%d-%d.nc" % (Y0, Y0))
if lmf:
    lm = xr.open_dataset(lmf[0])["lsm"];  lm = lm.isel(time_counter=0) if "time_counter" in lm.dims else lm
    lm = lm.reset_coords(drop=True)
mask = (LAT >= 50) & (LAT <= 80) & ((lm > 0.5) if lm is not None else True)
Wt = np.cos(np.deg2rad(LAT)).broadcast_like(REF)

def wmean(da, extra=None):
    m = mask & (extra if extra is not None else True)
    b = da.where(m); w = Wt.where(np.isfinite(b))
    return float(b.weighted(w.fillna(0)).mean().values)

def wcorr(x, y):
    m = mask & np.isfinite(x) & np.isfinite(y)
    xv = x.where(m).values.ravel(); yv = y.where(m).values.ravel(); wv = Wt.where(m).values.ravel()
    g = np.isfinite(xv) & np.isfinite(yv) & np.isfinite(wv); xv, yv, wv = xv[g], yv[g], wv[g]
    mx = np.average(xv, weights=wv); my = np.average(yv, weights=wv)
    return np.average((xv-mx)*(yv-my), weights=wv) / (
        np.sqrt(np.average((xv-mx)**2, weights=wv)) * np.sqrt(np.average((yv-my)**2, weights=wv)))

# candidate bias fields (model - ERA5)
b_t2m = M["t2m"] - E["t2m"]                                   # the residual cold (K)
inv_m = M["t925"] - M["t2m"]; inv_e = E["t925"] - E["t2m"]
b_inv = inv_m - inv_e                                         # inversion too strong? (K)
b_msl = (M["msl"] - E["msl"]) / 100.0                          # hPa (Siberian High)
b_tcc = M["tcc"] - E["tcc"]                                    # cloud deficit
inv2_m = M["t925"] - M["t1000"]; inv2_e = E["t925"] - E["t1000"]
b_inv2 = inv2_m - inv2_e                                       # 925-1000 inversion bias

print(f"=== WITHOUT-LPJG (full forest) DJFM vs ERA5, boreal land, {Y0}-{Y1} ===")
print(f"  residual 2m-T bias        : {wmean(b_t2m):+6.2f} K   (still too cold => 2nd driver)")
print(f"  low-level inversion bias  : {wmean(b_inv):+6.2f} K   (T925-T2m, model minus ERA5; + = too strong)")
print(f"  925-1000hPa inversion bias: {wmean(b_inv2):+6.2f} K")
print(f"  ERA5 inversion (T925-T2m) : {wmean(inv_e):+6.2f} K    model: {wmean(inv_m):+6.2f} K")
print(f"  MSLP bias (Siberian High) : {wmean(b_msl):+6.2f} hPa")
print(f"  cloud-cover bias          : {wmean(b_tcc):+6.2f}     (annual ERA5 clim, caveat)")
print(f"  skt-T2m (model sfc inv)   : {wmean(M['skt']-M['t2m']):+6.2f} K")
print("\n--- which bias explains WHERE it's too cold? spatial r with residual 2m-T bias ---")
print(f"  r(coldbias, -inversion bias) = {wcorr(b_t2m, -b_inv):+.2f}   (+ => stronger inversion => colder)")
print(f"  r(coldbias, -MSLP bias)      = {wcorr(b_t2m, -b_msl):+.2f}   (+ => higher MSLP => colder)")
print(f"  r(coldbias, +cloud bias)     = {wcorr(b_t2m,  b_tcc):+.2f}   (+ => more cloud => warmer/less cold)")

# figure
panels = [(b_t2m,"residual 2m-T bias (model-ERA5)","RdBu_r",8,"K"),
          (b_inv,"inversion bias T925-T2m (+ = too strong)","RdBu_r",6,"K"),
          (b_msl,"MSLP bias (Siberian High)","RdBu_r",10,"hPa"),
          (b_tcc,"cloud-cover bias","RdBu_r",0.3,"")]
fig=plt.figure(figsize=(13,10)); proj=ccrs.NorthPolarStereo() if HAVE_CARTOPY else None
for i,(da,t,cm,vl,u) in enumerate(panels,1):
    ax=fig.add_subplot(2,2,i,projection=proj) if HAVE_CARTOPY else fig.add_subplot(2,2,i)
    dd=da.where(mask)
    if HAVE_CARTOPY:
        ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
        pc=ax.pcolormesh(LON,LAT,dd,transform=ccrs.PlateCarree(),cmap=cm,vmin=-vl,vmax=vl,shading="auto")
        ax.coastlines(linewidth=0.4); ax.gridlines(linewidth=0.2)
    else:
        pc=ax.pcolormesh(LON,LAT,dd,cmap=cm,vmin=-vl,vmax=vl)
    ax.set_title(t+(f" [{u}]" if u else ""),fontsize=11); fig.colorbar(pc,ax=ax,shrink=0.7)
fig.suptitle(f"Second-driver search: WITHOUT-LPJG (full forest) DJFM vs ERA5, boreal NH ({Y0}-{Y1})",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.96])
OUT="/work/bb1469/a270092/eval/plots/djf_second_driver.png"; fig.savefig(OUT,dpi=130,bbox_inches="tight"); print("\nwrote",OUT)
