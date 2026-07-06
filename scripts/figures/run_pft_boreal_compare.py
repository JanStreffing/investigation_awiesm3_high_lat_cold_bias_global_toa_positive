#!/usr/bin/env python
# Boreal-belt dominant-PFT comparison: htessel spin-up (reference, no atm<->lpjg
# feedback) vs coupled Tuning_test_06 runs. Shows the boreal needleleaf -> C3 grass
# degradation over Siberia / N. Canada that co-locates with the cold bias.
import sys, os, numpy as np
sys.path.insert(0,'/work/ab0246/a270092/software/release_evaluation_tool2/scripts')
sys.path.insert(0,'/work/ab0246/a270092/software/release_evaluation_tool2')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.colors as mcolors, matplotlib.patches as mpatches
import cartopy.crs as ccrs, cartopy.feature as cfeature
from lpjg_helpers import read_lpjg_output, interpolate_dominant_pft, find_lpjg_latest_year

PFTS=['BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR','C3G','C4G']
PCOL={'BNE':'#0d47a1','BINE':'#1976d2','BNS':'#64b5f6','TeNE':'#4a148c','TeBS':'#7b1fa2',
      'IBS':'#ba68c8','TeBE':'#6a1b9a','TrBE':'#1b5e20','TrIBE':'#388e3c','TrBR':'#81c784',
      'C3G':'#fdd835','C4G':'#ff8f00','Barren':'#bdbdbd'}
TREES={'BNE','BINE','BNS','TeNE','TeBS','IBS','TeBE','TrBE','TrIBE','TrBR'}
GRASS={'C3G','C4G'}
EVAL='/work/bb1469/a270092/eval'
RUNS=[('/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_54_100YRES_2000YSPINUP_TCO95_CORE3/outdata/','Spin-up (HTESSEL, no feedback)'),
      ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/','06 Baseline (coupled)'),
      ('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1/outdata/','06V (coupled, aggressive)')]

fig=plt.figure(figsize=(18,7))
stats=[]
for i,(p,lab) in enumerate(RUNS):
    yr=find_lpjg_latest_year(p)
    df=read_lpjg_output(p,'lai.out',yr)
    lon,lat,grid,data,avail=interpolate_dominant_pft(df,yr,PFTS,1.0,0.2,5.0)
    ax=fig.add_subplot(1,3,i+1,projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    cols=[PCOL.get(p2,'#808080') for p2 in avail]
    cmap=mcolors.ListedColormap(cols); norm=mcolors.BoundaryNorm(np.arange(-0.5,len(avail),1),cmap.N)
    lo,la=np.meshgrid(lon,lat)
    ax.pcolormesh(lo,la,grid,cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),shading='auto',zorder=1)
    ax.add_feature(cfeature.OCEAN,facecolor='lightblue',alpha=0.4,zorder=2)
    ax.coastlines(resolution='110m',linewidth=0.4,zorder=3)
    ax.set_title(f'{lab}\n(yr {yr})',fontsize=11,fontweight='bold')
    # tree vs grass fraction of vegetated land north of 50N
    m=(la>=50)&np.isfinite(grid)
    names=np.array(avail)
    domnames=np.full(grid.shape,'',dtype=object)
    for k,nm in enumerate(avail): domnames[grid==k]=nm
    sel=domnames[m]
    ntot=np.sum(sel!='')
    ntree=np.sum(np.isin(sel,list(TREES))); ngrass=np.sum(np.isin(sel,list(GRASS)))
    stats.append((lab,yr,100*ntree/ntot,100*ngrass/ntot))
# legend
avail_all=[p for p in PFTS if p in PCOL]+['Barren']
patches=[mpatches.Patch(color=PCOL[p],label=p) for p in avail_all]
fig.legend(handles=patches,loc='lower center',ncol=7,fontsize=9,bbox_to_anchor=(0.5,-0.02))
fig.suptitle('Dominant PFT, boreal belt (>45°N): tree→grass loss in the coupled runs',
             fontsize=14,fontweight='bold',y=1.02)
plt.tight_layout()
plt.savefig(f'{EVAL}/plots/pft_boreal_compare.png',dpi=150,bbox_inches='tight')
print("saved pft_boreal_compare.png")
print(f"\n{'run':<38}{'year':>6}{'tree%>50N':>11}{'grass%>50N':>12}")
for lab,yr,tr,gr in stats: print(f"{lab:<38}{yr:>6}{tr:>11.1f}{gr:>12.1f}")
import csv
with open(f'{EVAL}/data/boreal_tree_grass_fraction.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['run','year','tree_pct_gt50N','grass_pct_gt50N'])
    for r in stats: w.writerow(r)
