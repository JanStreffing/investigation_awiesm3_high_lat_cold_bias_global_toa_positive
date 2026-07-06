#!/usr/bin/env python3
"""
NH sea-ice VOLUME trend, both runs, legs 1+2 (1350-1369). Volume is the
sensitive feedback indicator that area (concentration, winter-saturated) hides:
ice can thicken (storing the veg cooling) without area change. If WITH (LPJG,
veg-cooling) thickens relative to WITHOUT, that is the amplifying feedback.

m_ice = effective ice thickness [m] on native CORE3 mesh; volume = sum(m_ice*nod_area).
Usage: python seaice_volume_trend.py [Y0 Y1]
"""
import sys, glob, warnings
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

Y0, Y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1350, 1369)
years = np.arange(Y0, Y1 + 1)
FWO = "/work/bb1469/a270092/runtime/awicm3-develop/awicm3_noLPJG_CORE3_30y/outdata/fesom"
FWI = "/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/fesom"

md = xr.open_dataset("/work/ab0246/a270092/input/fesom2/core3/fesom.mesh.diag.nc")
narea = md["nod_area"].isel(nz=0).values            # m^2 per surface node
lat = xr.open_dataset(f"{FWO}/m_ice.fesom.{Y0}.nc")["lat"].values
nh = lat > 0; sh = lat < 0

def vol(base, y, hemi):                              # NH ice volume, 10^3 km^3
    f = glob.glob(f"{base}/m_ice.fesom.{y}.nc")
    if not f: return np.nan
    mi = xr.open_dataset(f[0])["m_ice"]
    td = "time" if "time" in mi.dims else "time_counter"
    m = mi.mean(td).values                           # annual-mean thickness [m]
    return np.nansum(np.where(hemi, m * narea, 0.0)) / 1e9 / 1e3   # m^3 -> km^3 -> 10^3 km^3

VWO = np.array([vol(FWO, y, nh) for y in years])
VWI = np.array([vol(FWI, y, nh) for y in years])

def trend(y, v):
    g = np.isfinite(v); y, v = y[g], v[g]
    a, b = np.polyfit(y, v, 1)
    se = np.sqrt(np.sum((v-(a*y+b))**2)/(len(y)-2))/(np.std(y)*np.sqrt(len(y)))
    return a*10, se*10

two, twe = trend(years, VWO); tiw, tie = trend(years, VWI)
print(f"=== NH sea-ice VOLUME [10^3 km^3], {Y0}-{Y1} ===\n")
print(f"  mean:  WITHOUT={np.nanmean(VWO):7.2f}   WITH={np.nanmean(VWI):7.2f}   diff={np.nanmean(VWI)-np.nanmean(VWO):+.2f}")
print(f"  trend/decade: WITHOUT={two:+.3f}±{twe:.3f}   WITH={tiw:+.3f}±{tie:.3f}   WITH-WITHOUT={tiw-two:+.3f}")
print("\n  (feedback signature = WITH thickening vs WITHOUT, i.e. WITH-WITHOUT trend > 0 & significant)")
print("\n year   WITHOUT    WITH")
for y, a, b in zip(years, VWO, VWI): print(f" {y}  {a:7.2f}  {b:7.2f}")

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(years, VWO, "o-", color="C0", label=f"WITHOUT LPJG ({two:+.2f}/dec)")
ax.plot(years, VWI, "s-", color="C3", label=f"WITH LPJG ({tiw:+.2f}/dec)")
for v,c in [(VWO,"C0"),(VWI,"C3")]:
    a,b=np.polyfit(years,v,1); ax.plot(years,a*years+b,"--",color=c,lw=1)
ax.axvline(Y0+9.5,color="k",lw=.3,ls=":"); ax.set_xlabel("year"); ax.set_ylabel("NH sea-ice volume [10^3 km^3]")
ax.set_title(f"NH sea-ice volume drift ({Y0}-{Y1}): veg-cooling feedback test"); ax.legend()
fig.tight_layout(); OUT="/work/bb1469/a270092/eval/plots/seaice_volume_trend.png"
fig.savefig(OUT,dpi=130,bbox_inches="tight"); print("\nwrote",OUT)
