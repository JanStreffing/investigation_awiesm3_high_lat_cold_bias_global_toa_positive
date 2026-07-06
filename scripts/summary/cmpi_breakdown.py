#!/usr/bin/env python
# Per-variable and per-region CMPI breakdown vs Baseline.
# Goal (per J.S.): CMPI/RMSD rise vs present-day obs is partly expected for a PI run;
# what matters is WHERE the error changes -- a tuning should move the fields/regions it
# physically targets, and NOT degrade unrelated ones. This shows that map of change.
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

OUT='/work/ab0246/a270092/software/release_evaluation_tool2/output'
R='/work/bb1469/a270092/eval'
RUNS=[('Tuning_test_06_Baseline','Baseline'),
      ('Tuning_test_06A_fesomA_albpnd028','06A meltpond-alb'),
      ('Tuning_test_06D_HRlike','06D pond-combo'),
      ('Tuning_test_06H_fesomH_combo_g_rvice018','06H pond+rvice'),
      ('Tuning_test_06O_1hcpl_mospp','06O 1hcpl+mospp'),
      ('Tuning_test_06T_1hcpl_mospp_kpplow','06T +kpplow'),
      ('Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1','06V aggressive')]
VARORDER=['siconc','tas','clt','pr','rlut','uas','vas','ua','zg','zos','mlotst','thetao','so']

def load(run):
    f=f"{OUT}/{run}/cmpi/frac/{run}_fraction.csv"
    if not os.path.exists(f): return None
    df=pd.read_csv(f, sep=r'\s+')
    df=df[df['Variable']!='CMPI']
    df['FracMeanError']=pd.to_numeric(df['FracMeanError'],errors='coerce')
    return df

dfs={lab:load(run) for run,lab in RUNS if load(run) is not None}
labs=[l for l in [lab for _,lab in RUNS] if l in dfs]
print("loaded:",labs)

# --- per variable ---
pv={lab: dfs[lab].groupby('Variable')['FracMeanError'].mean() for lab in labs}
PV=pd.DataFrame(pv).reindex(VARORDER)
PV.to_csv(R+'/data/cmpi_by_variable.csv')
# --- per region ---
pr={lab: dfs[lab].groupby('Region')['FracMeanError'].mean() for lab in labs}
PR=pd.DataFrame(pr)
PR.to_csv(R+'/data/cmpi_by_region.csv')

base='Baseline'
def heat(M, fname, title, vmaxabs=None):
    others=[l for l in labs if l!=base]
    D=M[others].sub(M[base],axis=0)        # Δ vs baseline
    if vmaxabs is None: vmaxabs=np.nanpercentile(np.abs(D.values),98) or 0.05
    fig,ax=plt.subplots(figsize=(1.6+1.1*len(others), 0.5+0.42*len(D.index)))
    im=ax.imshow(D.values, cmap='RdBu_r', vmin=-vmaxabs, vmax=vmaxabs, aspect='auto')
    ax.set_xticks(range(len(others))); ax.set_xticklabels(others,rotation=30,ha='right',fontsize=8)
    ax.set_yticks(range(len(D.index))); ax.set_yticklabels(D.index,fontsize=8)
    for i in range(len(D.index)):
        for j in range(len(others)):
            v=D.values[i,j]
            if np.isfinite(v): ax.text(j,i,f"{v:+.02f}",ha='center',va='center',fontsize=6.5,
                                       color='white' if abs(v)>0.6*vmaxabs else 'black')
    ax.set_title(title,fontsize=10,fontweight='bold')
    cb=fig.colorbar(im,ax=ax,shrink=0.8); cb.set_label('Δ frac. error vs Baseline\n(red = worse, blue = better)',fontsize=8)
    plt.tight_layout(); plt.savefig(R+'/plots/'+fname,dpi=160,bbox_inches='tight'); plt.close()
    return D

Dv=heat(PV,'cmpi_breakdown_by_variable.png','CMPI fractional error change vs Baseline — by variable')
Dr=heat(PR,'cmpi_breakdown_by_region.png','CMPI fractional error change vs Baseline — by region')

print("\n=== Δ frac error vs baseline BY VARIABLE (red>0 worse) ===")
print(Dv.round(3).to_string())
print("\n=== Δ frac error vs baseline BY REGION ===")
print(Dr.round(3).to_string())
print("\n=== global CMPI ===")
for run,lab in RUNS:
    f=f"{OUT}/{run}/cmpi/frac/{run}_fraction.csv"
    if os.path.exists(f):
        for line in open(f):
            if line.startswith('CMPI'): print(f"  {lab:<22}{float(line.split()[-1]):.3f}")
