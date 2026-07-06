#!/usr/bin/env python
# Test the causal direction: in the COUPLED baseline run, does boreal vegetation
# drift (needleleaf trees -> grass) DRIVE the cold bias? The atmosphere has no
# long-term memory, so if veg drives it, tree LAI and boreal 2m-T should fall
# together over the run (yr 1350->1379). Plot both.
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

RUN='/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata'
EVAL='/work/bb1469/a270092/eval'
LATMIN,LATMAX=55,75                         # boreal belt
TREE=['BNE','BINE','BNS']                   # boreal needleleaf
import subprocess
def cdo(args):
    return subprocess.run(args,capture_output=True,text=True)

# --- 1) boreal tree & grass LAI per year from lai.out (all chunks) ---
rows=[]
for f in sorted(glob.glob(f'{RUN}/lpj_guess/*/run1/lai.out')):
    df=pd.read_csv(f,sep=r'\s+')
    rows.append(df)
lai=pd.concat(rows).drop_duplicates(subset=['Lon','Lat','Year'])
bor=lai[(lai.Lat>=LATMIN)&(lai.Lat<=LATMAX)]
# area weight by cos(lat)
w=np.cos(np.deg2rad(bor.Lat))
g=bor.groupby('Year')
def wmean(col):
    return g.apply(lambda d: np.average(d[col], weights=np.cos(np.deg2rad(d.Lat))))
tree_lai=sum(wmean(c) for c in TREE)
grass_lai=wmean('C3G')
years=tree_lai.index.values

# --- 2) boreal 2m T, surface albedo, snow albedo per year (Siberia+Canada land) ---
import xarray as xr
def boxmean(vm, lat, lon):
    wlat=np.cos(np.deg2rad(lat))
    sib=vm.where((lat>=LATMIN)&(lat<=LATMAX)&(lon>=60)&(lon<=140))
    can=vm.where((lat>=LATMIN)&(lat<=LATMAX)&(lon>=235)&(lon<=300))
    aw=lambda x: float(x.weighted(wlat).mean(('lat','lon')))
    return 0.5*(aw(sib)+aw(can))
def openvar(y,var,fpvar=None):
    fp=f'{RUN}/oifs/atm_remapped_1m_{fpvar or var}_{y}-{y}.nc'
    if not os.path.exists(fp): return None
    ds=xr.open_dataset(fp,use_cftime=True)
    return ds, ds['lat'], ds['lon']%360
def annsum(ds,var):
    v=ds[var]; td='time_counter' if 'time_counter' in v.dims else 'time'
    return v.sum(td)
def annmean(ds,var):
    v=ds[var]; td='time_counter' if 'time_counter' in v.dims else 'time'
    return v.mean(td)

t2m=[]; salb=[]; snowalb=[]
for y in years:
    d,lat,lon=openvar(y,'2t'); t2m.append(boxmean(annmean(d,'2t'),lat,lon))
    # surface albedo = 1 - SSR/SSRD (net/down surface shortwave; accumulation cancels in ratio)
    try:
        dssr,lat,lon=openvar(y,'ssr'); dssrd,_,_=openvar(y,'ssrd')
        alb=1.0-(annsum(dssr,'ssr')/annsum(dssrd,'ssrd'))
        salb.append(boxmean(alb,lat,lon))
    except Exception as e:
        salb.append(np.nan)
    # snow albedo (asn) if present
    da=openvar(y,'asn')
    if da is not None:
        dd,lat,lon=da; snowalb.append(boxmean(annmean(dd,'asn'),lat,lon))
    else: snowalb.append(np.nan)
t2m=np.array(t2m); salb=np.array(salb); snowalb=np.array(snowalb)

def tr(v):
    m=np.isfinite(v); return np.polyfit(years[m],np.array(v)[m],1)[0]*10 if m.sum()>2 else np.nan

# --- plot: 3 stacked panels sharing year axis ---
fig,axs=plt.subplots(3,1,figsize=(9,9),sharex=True)
axs[0].plot(years,tree_lai.values,'-o',color='#0d47a1',label=f'needleleaf LAI (BNE+BINE+BNS)  [{tr(tree_lai.values):+.2f}/dec]')
axs[0].plot(years,grass_lai.values,'-s',color='#c9a800',mec='k',mew=0.3,label=f'C3 grass LAI  [{tr(grass_lai.values):+.2f}/dec]')
axs[0].set_ylabel('LAI (m²/m²)'); axs[0].legend(fontsize=8); axs[0].grid(alpha=0.25)
axs[0].set_title('Boreal 55–75°N (Siberia + N. America land): does veg drift drive the cooling?',fontweight='bold')
axs[1].plot(years,salb,'-D',color='#00897b',label=f'surface albedo 1−SSR/SSRD  [{tr(salb):+.4f}/dec]')
if np.isfinite(snowalb).any():
    axs[1].plot(years,snowalb,'-v',color='#5e35b1',label=f'snow albedo asn  [{tr(snowalb):+.4f}/dec]')
axs[1].set_ylabel('albedo'); axs[1].legend(fontsize=8); axs[1].grid(alpha=0.25)
axs[2].plot(years,t2m-273.15,'-^',color='#d62728',label=f'2 m temperature  [{tr(t2m):+.2f} K/dec]')
axs[2].set_ylabel('2 m T [°C]'); axs[2].set_xlabel('Model year'); axs[2].legend(fontsize=8); axs[2].grid(alpha=0.25)
plt.tight_layout(); plt.savefig(f'{EVAL}/plots/baseline_boreal_evolution.png',dpi=160,bbox_inches='tight')
import csv
with open(f'{EVAL}/data/baseline_boreal_evolution.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['year','tree_LAI','grass_LAI','sfc_albedo','snow_albedo','boreal_2mT_K'])
    for i,y in enumerate(years): w.writerow([y,f"{tree_lai.values[i]:.4f}",f"{grass_lai.values[i]:.4f}",f"{salb[i]:.4f}",f"{snowalb[i]:.4f}",f"{t2m[i]:.3f}"])
print("saved baseline_boreal_evolution.png")
print(f"{'yr':>6}{'treeLAI':>9}{'grassLAI':>9}{'sfcAlb':>8}{'2mT_C':>8}")
for i,y in enumerate(years): print(f"{y:>6}{tree_lai.values[i]:>9.3f}{grass_lai.values[i]:>9.3f}{salb[i]:>8.3f}{t2m[i]-273.15:>8.2f}")
print(f"\ntrends/decade: tree {tr(tree_lai.values):+.3f}  grass {tr(grass_lai.values):+.3f}  sfcAlb {tr(salb):+.4f}  2mT {tr(t2m):+.3f} K")
print(f"correlation tree-LAI vs 2mT: r={np.corrcoef(tree_lai.values,t2m)[0,1]:+.2f}; albedo vs 2mT: r={np.corrcoef(salb[np.isfinite(salb)],t2m[np.isfinite(salb)])[0,1]:+.2f}; tree vs albedo: r={np.corrcoef(tree_lai.values[np.isfinite(salb)],salb[np.isfinite(salb)])[0,1]:+.2f}")
