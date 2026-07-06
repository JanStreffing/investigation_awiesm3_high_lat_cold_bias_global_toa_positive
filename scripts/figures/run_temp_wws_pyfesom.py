#!/usr/bin/env python
# Weddell-Sea potential-temperature depth-latitude section (pyfesom2 fallback for
# a270234's tripyview temp_fesom_core3_WWS.ipynb, which cannot decode these runs'
# pre-1582 'gregorian' calendar). Austral-winter (Sep) temperature averaged over the
# Weddell longitude band 50W-30W, vs depth & latitude: Baseline + anomalies.
import os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import pyfesom2 as pf

BASE='/work/bb1469/a270270/runtime/awiesm3-v3.4'; MESH='/work/ab0246/a270092/input/fesom2/core3/'
OUTP='/work/bb1469/a270092/eval/plots'
YEARS=range(1370,1380); WINTER=8
LONBAND=(-50,-30); LATR=(-78,-55)
RUNS=[('Tuning_test_06_Baseline','Baseline'),
      ('Tuning_test_06O_1hcpl_mospp','06O 1hcpl+mospp'),
      ('Tuning_test_06T_1hcpl_mospp_kpplow','06T +kpplow'),
      ('Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1','06V aggressive')]

mesh=pf.load_mesh(MESH, usepickle=False, usejoblib=False)
zlev=np.abs(np.array(mesh.zlev)); zmid=0.5*(zlev[:-1]+zlev[1:])      # layer mid-depths
x=mesh.x2; y=mesh.y2
inband=(x>=LONBAND[0])&(x<LONBAND[1])&(y>=LATR[0])&(y<LATR[1])
latbins=np.arange(LATR[0],LATR[1]+1,1.0); latc=0.5*(latbins[:-1]+latbins[1:])

def section(run):
    d=pf.get_data(f'{BASE}/{run}/outdata/fesom/','temp',YEARS,mesh,how=None,compute=False)
    da=d[WINTER::12].mean(dim='time').values         # (nod2, nz1)
    nz=da.shape[1]
    sec=np.full((len(latc),nz),np.nan)
    idx=np.where(inband)[0]
    ybin=np.digitize(y[idx],latbins)-1
    for b in range(len(latc)):
        sel=idx[ybin==b]
        if sel.size:
            with np.errstate(invalid='ignore'):
                sec[b,:]=np.nanmean(da[sel,:],axis=0)
    return sec  # (lat, depth)

print("computing sections ...")
secs={lab:section(run) for run,lab in RUNS}
nz=secs['Baseline'].shape[1]; depth=zmid[:nz]
base=secs['Baseline']
fig,axes=plt.subplots(2,2,figsize=(13,9),sharex=True,sharey=True)
axes=axes.flatten()
for i,(run,lab) in enumerate(RUNS):
    ax=axes[i]
    if lab=='Baseline':
        im=ax.contourf(latc,depth,base.T,levels=np.linspace(-2,2,21),cmap='RdYlBu_r',extend='both')
        cb=plt.colorbar(im,ax=ax); cb.set_label('°C')
        ax.set_title('Baseline (absolute)',fontweight='bold')
    else:
        an=secs[lab]-base
        im=ax.contourf(latc,depth,an.T,levels=np.linspace(-1.5,1.5,21),cmap='RdBu_r',extend='both')
        cb=plt.colorbar(im,ax=ax); cb.set_label('Δ°C vs Base')
        ax.set_title(f'{lab} − Baseline',fontweight='bold')
    ax.invert_yaxis(); ax.set_ylim(2000,0); ax.set_xlabel('Latitude'); ax.set_ylabel('Depth [m]')
fig.suptitle('West Weddell Sea (50°W–30°W) austral-winter potential temperature, yr1370-79',
             fontsize=13,fontweight='bold')
plt.tight_layout(); plt.savefig(f'{OUTP}/temp_wws_section_tuning.png',dpi=150,bbox_inches='tight')
print("saved temp_wws_section_tuning.png")
