#!/usr/bin/env python
# Adapted from a270234's mld_sic_sit_compare.ipynb to the coupled Tuning_test_06 runs.
# Antarctic austral-winter (Sep) MLD3 / sea-ice conc (a_ice) / thickness (m_ice):
#  - top-left panel = Baseline absolute
#  - remaining panels = run MINUS baseline (anomaly) -> shows WHERE each tuning acts.
import os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import pyfesom2 as pf

BASE='/work/bb1469/a270270/runtime/awiesm3-v3.4'
MESH='/work/ab0246/a270092/input/fesom2/core3/'
OUTP='/work/bb1469/a270092/eval/plots'
os.makedirs(OUTP, exist_ok=True)
YEARS=range(1370,1380); WINTER=8     # September (SH late winter / sea-ice max)

RUNS=[('Tuning_test_06_Baseline','Baseline'),
      ('Tuning_test_06A_fesomA_albpnd028','06A meltpond'),
      ('Tuning_test_06N_mospp','06N mospp'),
      ('Tuning_test_06O_1hcpl_mospp','06O 1hcpl+mospp'),
      ('Tuning_test_06T_1hcpl_mospp_kpplow','06T +kpplow'),
      ('Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1','06V aggressive')]

print("loading mesh ...")
mesh=pf.load_mesh(MESH, usepickle=False, usejoblib=False)
def wmean(run,var):
    d=pf.get_data(f'{BASE}/{run}/outdata/fesom/',var,YEARS,mesh,how=None,compute=False)
    return d[WINTER::12,:].mean(dim='time').values

# var, baseline abs levels, anomaly half-range, latmax, title, fname
specs=[('MLD3', (-2000,-100,21), 300, -50,'Antarctic winter MLD (Sep, yr1370-79)', 'mld_winter_antarctic.png','m'),
       ('a_ice',(0.0,1.0,21),    0.3, -50,'Antarctic winter sea-ice concentration (Sep)','sic_winter_antarctic.png','frac'),
       ('m_ice',(0.0,2.5,21),    1.0, -50,'Antarctic winter sea-ice thickness (Sep)','sit_winter_antarctic.png','m')]

for var,abslev,anom,latmax,title,fn,unit in specs:
    print(f"=== {var} ===")
    fields={}
    for run,lab in RUNS:
        try: fields[lab]=wmean(run,var); print(f"  {lab} ok")
        except Exception as e: print(f"  {lab} FAIL {e}")
    base=fields['Baseline']
    others=[lab for _,lab in RUNS if lab!='Baseline' and lab in fields]
    anom_arrs=[fields[l]-base for l in others]
    # auto-scale the anomaly colorbar to the data (98th pct of |anomaly|) so the
    # signal fills a symmetric diverging scale instead of saturating/vanishing.
    import math
    # per-variable percentile: sea-ice thickness has coastal fast-ice outliers that
    # inflate the 98th pct and wash out the broad pack-thickening signal, so use a
    # lower percentile for it; MLD/sic keep 98.
    pct=88 if var=='m_ice' else 98
    aflat=np.concatenate([a[np.isfinite(a)] for a in anom_arrs]) if anom_arrs else np.array([anom])
    arng=float(np.nanpercentile(np.abs(aflat),pct)) if aflat.size else anom
    if arng>0: arng=round(arng, -int(math.floor(math.log10(arng)))+1)   # 2 sig figs
    else: arng=anom
    arrs=[base]+anom_arrs
    titles=['Baseline (abs.)']+[f'{l} − Base' for l in others]
    levels=[abslev]+[(-arng,arng,21)]*len(others)
    cmaps=['Spectral_r']+['RdBu_r']*len(others)
    print(f"  anomaly colorbar half-range (98th pct) = {arng}")
    n=len(arrs); ncol=3; nrow=int(np.ceil(n/ncol))
    import cartopy.crs as ccrs, cartopy.feature as cfeature
    import matplotlib.tri as mtri
    x2,y2=mesh.x2, mesh.y2
    # native FESOM triangles, restricted to the Antarctic region (all vertices
    # south of -38) -> no interpolation, just far fewer triangles so tripcolor is fast
    elem_all=mesh.elem
    elem=elem_all[np.all(y2[elem_all] < -38.0, axis=1)]
    xt=x2[elem]; cyclic=(xt.max(axis=1)-xt.min(axis=1))>100.0   # dateline-spanning
    fig=plt.figure(figsize=(5*ncol,4.4*nrow))
    for i,(a,t,lv,cm) in enumerate(zip(arrs,titles,levels,cmaps)):
        # mask triangles touching a non-finite vertex -> no interpolation, native mesh
        bad=cyclic | np.any(~np.isfinite(a[elem]),axis=1)
        triang=mtri.Triangulation(x2,y2,triangles=elem,mask=bad)
        ax=fig.add_subplot(nrow,ncol,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180,180,-90,latmax],ccrs.PlateCarree())
        im=ax.tripcolor(triang, np.nan_to_num(a), cmap=cm, vmin=lv[0], vmax=lv[1],
                        shading='flat', transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND,facecolor='lightgray',zorder=2)
        ax.coastlines(linewidth=0.3,zorder=3)
        ax.set_title(t,fontsize=11)
        cb=plt.colorbar(im,ax=ax,shrink=0.7,pad=0.04); cb.ax.tick_params(labelsize=7)
        cb.set_label(unit,fontsize=8)
    fig.suptitle(title,y=1.0,fontsize=14,fontweight='bold')
    plt.tight_layout(); plt.savefig(f'{OUTP}/{fn}',dpi=135,bbox_inches='tight'); plt.close('all')
    print(f"  saved {fn}")
print("DONE")
