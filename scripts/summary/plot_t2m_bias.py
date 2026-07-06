import sys, os, glob, math
sys.path.insert(0,'/work/ab0246/a270092/software/release_evaluation_tool2')
import numpy as np, xarray as xr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from bg_routines.ipcc_cmaps import get_bias_cmap

R='/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview'
T=R+'/t2m'; P=R+'/plots'
os.makedirs(P, exist_ok=True)
levels=[-8.0,-5.0,-3.0,-2.0,-1.0,-.6,-.2,.2,.6,1.0,2.0,3.0,5.0,8.0]
cmap=get_bias_cmap('2t')

era=xr.open_dataset(T+'/era5_clim.nc')
evar=[v for v in era.data_vars if 't2m' in v.lower() or v=='t2m'][0]
era5=np.squeeze(era[evar].values)
lat=era['lat'].values if 'lat' in era else era['latitude'].values
lon=era['lon'].values if 'lon' in era else era['longitude'].values

def wrmsd(diff):
    coslat=np.cos(np.deg2rad(lat)); w=np.sqrt(coslat)
    W=np.broadcast_to(w[:,None], diff.shape)
    return math.sqrt(np.sum(W*diff**2)/np.sum(W))
def wmean(diff):
    coslat=np.cos(np.deg2rad(lat)); w=np.sqrt(coslat)
    W=np.broadcast_to(w[:,None], diff.shape)
    return np.sum(W*diff)/np.sum(W)

# pretty labels
PRETTY={
 'Tuning_test_06_Baseline':'06 Baseline',
 'Tuning_test_06A_fesomA_albpnd028':'06A albpnd0.28',
 'Tuning_test_06D_HRlike':'06D HRlike (pond)',
 'Tuning_test_06H_fesomH_combo_g_rvice018':'06H meltpond+rvice',
 'Tuning_test_06O_1hcpl_mospp':'06O 1hcpl+mospp',
 'Tuning_test_06T_1hcpl_mospp_kpplow':'06T +kpplow',
 'Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1':'06V +entstpc3',
}

def load_bias(run):
    f=T+f'/{run}_clim.nc'
    ds=xr.open_dataset(f)
    mvar=[v for v in ds.data_vars if v=='2t' or '2t' in v][0]
    m=np.squeeze(ds[mvar].values)
    return m-era5

def plot_one(run, ax, title):
    diff=load_bias(run)
    r=wrmsd(diff); b=wmean(diff)
    d2, lon2 = add_cyclic_point(diff, coord=lon)
    ax.set_global(); ax.add_feature(cfeature.COASTLINE, zorder=3, linewidth=0.4)
    imf=ax.contourf(lon2, lat, d2, cmap=cmap, levels=levels, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax.set_title(title, fontweight='bold', fontsize=10)
    props=dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.6)
    ax.text(0.02,0.20,f'rmsd={r:.3f}',transform=ax.transAxes,fontsize=8,va='top',bbox=props,zorder=4)
    ax.text(0.02,0.10,f'bias={b:.3f}',transform=ax.transAxes,fontsize=8,va='top',bbox=props,zorder=4)
    return imf, r, b

cands=list(PRETTY.keys())
# individual maps
for run in cands:
    fig,ax=plt.subplots(figsize=(9,5),subplot_kw={'projection':ccrs.EqualEarth()})
    imf,r,b=plot_one(run,ax,PRETTY[run]+' − ERA5  (2m T, yr1350-1379)')
    cb=fig.colorbar(imf,ax=ax,orientation='horizontal',ticks=levels,pad=0.05,shrink=0.8)
    cb.set_label('K'); 
    for lab in cb.ax.xaxis.get_ticklabels()[::2]: lab.set_visible(False)
    plt.savefig(P+f'/t2m_bias_{PRETTY[run].split()[0]}.png',dpi=150,bbox_inches='tight'); plt.close()

# multipanel
fig,axes=plt.subplots(3,3,figsize=(16,11),subplot_kw={'projection':ccrs.EqualEarth()})
axes=axes.flatten()
imf=None
for i,run in enumerate(cands):
    imf,_,_=plot_one(run,axes[i],PRETTY[run])
for j in range(len(cands),len(axes)): axes[j].axis('off')
cax=fig.add_axes([0.25,0.06,0.5,0.02])
cb=fig.colorbar(imf,cax=cax,orientation='horizontal',ticks=levels); cb.set_label('2m air temperature bias vs ERA5 [K]')
for lab in cb.ax.xaxis.get_ticklabels()[::2]: lab.set_visible(False)
fig.suptitle('2m air temperature bias vs ERA5 (climatology yr 1350-1379)',fontsize=14,fontweight='bold',y=0.93)
plt.savefig(P+'/t2m_bias_multipanel.png',dpi=150,bbox_inches='tight'); plt.close()

# part8-style RMSE table for ALL 22 runs
print("=== part8-style (sqrt-coslat) RMSE all runs ===")
rows=[]
for f in sorted(glob.glob(T+'/Tuning_test_06*_clim.nc')):
    run=os.path.basename(f).replace('_clim.nc','')
    diff=load_bias(run)
    rows.append((run, wrmsd(diff), wmean(diff)))
for run,r,b in sorted(rows,key=lambda x:x[1]):
    print(f"{run:<46}{r:>8.3f}{b:>9.3f}")
import csv
with open(R+'/data/t2m_rmse_part8style.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['run','rmsd_sqrtcoslat','meanbias'])
    for row in sorted(rows,key=lambda x:x[1]): w.writerow([row[0],f"{row[1]:.4f}",f"{row[2]:.4f}"])
print("plots written to",P)
