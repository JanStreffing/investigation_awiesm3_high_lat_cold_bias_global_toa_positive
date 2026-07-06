#!/usr/bin/env python3
"""
Decisive climate-vs-competition figure for the v3 offline LPJG spin-up driven by
the corrected (flux-fixed) feedback-free OIFS-AMIP ~PI forcing.

Two questions, two rows of NH-polar maps:
  Row 1 (is the FORCING adequate / cold?):  AMIP warmest-month T2m [gate twmin_est=5],
        annual GDD5 [gate 500], coldest-month T2m [evergreen gate tcmin_est=-30].
  Row 2 (did VEG establish / close canopy?): v3 high-veg cvh(FRACH), TREEFPC, GRASSFPC.
If the forcing clears the gates (warm enough) yet cvh stays ~0.2 grass-dominated,
the failure is LPJG-internal competition (tune veg), not climate (tune atmosphere).
"""
import numpy as np, xarray as xr, pandas as pd, warnings, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    import cartopy.crs as ccrs, cartopy.feature as cfeature; HAVE=True
except Exception: HAVE=False

S="/tmp/claude-24456/-work-bb1469-a270092/e5d44e3f-0f81-4371-82f0-2dd6b5fb0d61/scratchpad"
FORC="/work/ab0246/a270092/input/lpj-guess/oifs_forcing/AMIP_noLPJG_1d_1870-1879_TCO95_PI.nc"
BOX=dict(la=(55,70), lo=(60,140))   # Siberian boreal (report box)

# ---------- AMIP forcing climate (per cell) ----------
d=xr.open_dataset(FORC, decode_times=False)
alat=d['lat'].values; alon=d['lon'].values
tas=d['tas'].values-273.15          # (3652, ncell) degC
import cftime
tc=d['time_counter']; dates=cftime.num2date(tc.values, tc.attrs['units'], tc.attrs.get('calendar','gregorian'))
mon=np.array([x.month for x in dates]); yr=np.array([x.year for x in dates])
# monthly climatology -> warmest / coldest month
mm=np.stack([tas[mon==m].mean(0) for m in range(1,13)])     # (12, ncell)
warmest=mm.max(0); coldest=mm.min(0)
# GDD5 = mean annual sum of max(0, T-5)
g=np.maximum(tas-5.0,0.0)
gdd5=np.stack([g[yr==y].sum(0) for y in range(1870,1880)]).mean(0)

# ---------- v3 spin-up equilibrium veg (per cell) ----------
fpc=pd.read_csv(f"{S}/v3_fpc_eq.txt", delim_whitespace=True)
vlat=fpc['Lat'].values; vlon=fpc['Lon'].values

def boxmask(la,lo): return (la>=BOX['la'][0])&(la<=BOX['la'][1])&(lo>=BOX['lo'][0])&(lo<=BOX['lo'][1])
def wbox(la,lo,v):
    m=boxmask(la,lo); w=np.cos(np.deg2rad(la[m]))
    return float(np.average(v[m], weights=w)) if m.sum() else np.nan

# ---------- print the conclusion stats ----------
am_bx=boxmask(alat,alon)
print(f"=== AMIP forcing, Siberian boreal box ({BOX['la'][0]}-{BOX['la'][1]}N, {BOX['lo'][0]}-{BOX['lo'][1]}E) ===")
print(f"  warmest-month T = {wbox(alat,alon,warmest):.1f} C   (gate twmin_est=+5; % cells >=5: {100*np.average((warmest[am_bx]>=5),weights=np.cos(np.deg2rad(alat[am_bx]))):.0f}%)")
print(f"  GDD5            = {wbox(alat,alon,gdd5):.0f}    (gate 500;      % cells >=500: {100*np.average((gdd5[am_bx]>=500),weights=np.cos(np.deg2rad(alat[am_bx]))):.0f}%)")
print(f"  coldest-month T = {wbox(alat,alon,coldest):.1f} C   (evergreen gate tcmin_est=-30; % cells < -30: {100*np.average((coldest[am_bx]<-30),weights=np.cos(np.deg2rad(alat[am_bx]))):.0f}%)")
print(f"\n=== v3 LPJG equilibrium veg, same box ===")
for c in ['FRACH','TREEFPC','GRASSFPC','FORESTFPC','BNE','BINE','BNS','C3G']:
    if c in fpc: print(f"  {c:9s} = {wbox(vlat,vlon,fpc[c].values):.3f}")

# ---------- figure (interpolated to regular grid + pcolormesh, per the eval template) ----------
if HAVE:
    sys.path.insert(0,'/work/ab0246/a270092/software/release_evaluation_tool2/scripts')
    from lpjg_helpers import interpolate_to_grid
    wrap=lambda lo: ((np.asarray(lo,dtype=float)+180)%360)-180
    fpc=fpc.copy(); fpc['Lon']=wrap(fpc['Lon'].values)
    amip=pd.DataFrame({'Lon':wrap(alon),'Lat':alat,'Year':3900,
                       'warmest':warmest,'gdd5':gdd5,'coldest':coldest})
    panels=[
      (amip,'warmest', "AMIP forcing: warmest-month T2m [$^\\circ$C]\n(gate twmin_est $=+5$)", 0,18,"YlOrRd"),
      (amip,'gdd5',    "AMIP forcing: GDD5\n(gate $=500$)",                                     0,1400,"YlGn"),
      (amip,'coldest', "AMIP forcing: coldest-month T2m [$^\\circ$C]\n(evergreen gate $=-30$)", -45,-10,"Blues_r"),
      (fpc,'FRACH',    "v3 spin-up: high-veg cover cvh (FRACH)\n(closed taiga $0.4$--$0.7$)",   0,0.7,"YlGn"),
      (fpc,'TREEFPC',  "v3 spin-up: TREEFPC",                                                   0,0.7,"YlGn"),
      (fpc,'GRASSFPC', "v3 spin-up: GRASSFPC",                                                  0,0.7,"YlOrBr"),
    ]
    fig=plt.figure(figsize=(17,11.8)); proj=ccrs.NorthPolarStereo()
    for i,(df,var,ttl,vmn,vmx,cm) in enumerate(panels,1):
        lon,lat,grid,_=interpolate_to_grid(df,var,3900,grid_res=1.0)
        lo,la=np.meshgrid(lon,lat)
        ax=fig.add_subplot(2,3,i,projection=proj); ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
        ax.add_feature(cfeature.OCEAN,facecolor="#eef3f7",zorder=0)
        im=ax.pcolormesh(lo,la,grid,cmap=cm,vmin=vmn,vmax=vmx,transform=ccrs.PlateCarree(),shading="auto",zorder=1)
        ax.coastlines(linewidth=0.4,zorder=3); ax.set_title(ttl,fontsize=9.5)
        fig.colorbar(im,ax=ax,shrink=0.7,pad=0.03)
    fig.suptitle("Decisive test: OIFS-AMIP (obs-SST, feedback-free ~PI) forcing vs the LPJG spin-up it drives  "
                 "(equilibrium yr 3900, NH boreal)", fontsize=13, fontweight='bold', y=1.0)
    fig.tight_layout(rect=[0,0,1,0.96], h_pad=4.0)
    OUT="/work/bb1469/a270092/eval/plots/v3_decisive_climate_vs_veg.png"
    fig.savefig(OUT,dpi=135,bbox_inches="tight"); print("\nwrote",OUT)
