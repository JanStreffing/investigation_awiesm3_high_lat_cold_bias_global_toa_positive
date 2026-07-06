#!/usr/bin/env python
# Boreal needleleaf & C3-grass LAI: spin-up end-state (the pool LPJG is started from,
# yr3850 state) vs coupled baseline END (yr1379) and their difference. Uses the eval
# tool's own LPJG gridding (lpjg_helpers.interpolate_to_grid, as in part24).
import sys, glob, numpy as np, pandas as pd
sys.path.insert(0,'/work/ab0246/a270092/software/release_evaluation_tool2/scripts')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
from lpjg_helpers import interpolate_to_grid

EVAL='/work/bb1469/a270092/eval'
SP=sorted(glob.glob('/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_54_100YRES_2000YSPINUP_TCO95_CORE3/outdata/lpj_guess/*/run1/lai.out'))[-1]
BL=sorted(glob.glob('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_06_Baseline/outdata/lpj_guess/*/run1/lai.out'))[-1]

def load(f):
    df=pd.read_csv(f, sep=r'\s+'); df['NLtree']=df['BNE']+df['BINE']+df['BNS']; return df
dsp=load(SP); dbl=load(BL)
spy=int(dsp.Year.max()); bly=1379

def grid(df,var,yr):
    lon,lat,g,_=interpolate_to_grid(df,var,yr,1.0); return lon,lat,g

rows=[('Needleleaf LAI (BNE+BINE+BNS)','NLtree',2.0,'YlGn','BrBG'),
      ('C3 grass LAI','C3G',2.0,'YlOrBr','RdBu_r')]
fig=plt.figure(figsize=(16,9)); proj=ccrs.NorthPolarStereo()
for r,(name,var,vmax,cmap,dcmap) in enumerate(rows):
    lo,la,gsp=grid(dsp,var,spy)
    _,_,gbl=grid(dbl,var,bly)
    gdiff=gbl-gsp
    sets=[(f'{name}\nspin-up end (init, yr{spy})',gsp,0,vmax,cmap,'LAI'),
          (f'{name}\nbaseline coupled (yr{bly})',gbl,0,vmax,cmap,'LAI'),
          (f'{name}\nbaseline − spin-up',gdiff,-1.0,1.0,dcmap,'Δ LAI')]
    LO,LA=np.meshgrid(lo,la)
    for c,(t,z,vm,vx,cm,cl) in enumerate(sets):
        ax=fig.add_subplot(2,3,r*3+c+1,projection=proj)
        ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
        im=ax.pcolormesh(LO,LA,np.ma.masked_invalid(z),cmap=cm,vmin=vm,vmax=vx,
                         transform=ccrs.PlateCarree(),shading='auto')
        ax.add_feature(cfeature.OCEAN,facecolor='lightblue',alpha=0.6,zorder=2)
        ax.coastlines(linewidth=0.4,zorder=3)
        ax.set_title(t,fontsize=9.5,fontweight='bold')
        cb=plt.colorbar(im,ax=ax,shrink=0.62,pad=0.03); cb.ax.tick_params(labelsize=7); cb.set_label(cl,fontsize=8)
fig.suptitle('Boreal LAI: spin-up start-state vs coupled baseline end (yr1379)',fontsize=13,fontweight='bold',y=1.01)
plt.tight_layout(); out=f'{EVAL}/plots/pft_lai_boreal_diff.png'
plt.savefig(out,dpi=150,bbox_inches='tight'); print('saved',out)
