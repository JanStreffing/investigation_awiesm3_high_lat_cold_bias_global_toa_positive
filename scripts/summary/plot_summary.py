import os, glob, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview'
# load toa
toa={}
for f in glob.glob(R+'/toa/*_toa.txt'):
    run=os.path.basename(f).replace('_toa.txt',''); d=np.loadtxt(f); yr=d[:,0]; v=d[:,1]
    toa[run]=dict(last10=v[yr>=1370].mean(), mean=v.mean(), trend=np.polyfit(yr,v,1)[0]*10, yr=yr, v=v)
# load part8 rmse
rm={}
for line in open(R+'/data/t2m_rmse_part8style.csv'):
    p=line.strip().split(',')
    if p[0]=='run' or len(p)<3: continue
    rm[p[0]]=dict(rmse=float(p[1]), bias=float(p[2]))

short=lambda r: r.replace('Tuning_test_06','06').replace('_fesom','_').split('_')[0] if False else r
def lab(r):
    s=r.replace('Tuning_test_','')
    return s
# categories
cat={}
def setcat(keys,c):
    for k in keys: cat[k]=c
setcat(['Tuning_test_06_Baseline'],'Baseline')
setcat(['Tuning_test_06A_fesomA_albpnd028','Tuning_test_06B_fesomB_albpnd028_rfrac075',
        'Tuning_test_06C_fesomC_albpnd028_rfrac075_pndaspect13','Tuning_test_06D_HRlike',
        'Tuning_test_06G_fesomG_stronger_meltpond.','Tuning_test_06L_extra_strong_meltpond',
        'Tuning_test_06H_fesomH_combo_g_rvice018'],'Sea-ice meltpond')
setcat(['Tuning_test_06Q_ralbsead0045','Tuning_test_06R_openwater_albedo0075'],'Ocean/ice albedo')
setcat(['Tuning_test_06S_evp0'],'Sea-ice dynamics')
setcat(['Tuning_test_06E_fesomE_rvice018_ggaussb_m06','Tuning_test_06F_fesomF_F_lrdalb_true',
        'Tuning_test_06P_fesomP_entstpc3_1'],'Atmosphere physics')
setcat(['Tuning_test_06N_mospp','Tuning_test_06O_1hcpl_mospp','Tuning_test_06O4_1hcpl_mospp_kv0012',
        'Tuning_test_06O5_1hcpl_mospp_kv002','Tuning_test_06M_1hcpl','Tuning_test_06T_1hcpl_mospp_kpplow',
        'Tuning_test_06U_1hcpl_mospp_kpplow_openwater_albedo0075','Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1'],
        'Coupling/ocean mixing')
colors={'Baseline':'k','Sea-ice meltpond':'#1f77b4','Ocean/ice albedo':'#17becf','Sea-ice dynamics':'#9467bd',
        'Atmosphere physics':'#2ca02c','Coupling/ocean mixing':'#d62728'}

runs=sorted(toa)
# ---- scatter: TOA imbalance (last decade) vs T2m RMSE ----
fig,ax=plt.subplots(figsize=(10,7))
seen=set()
for r in runs:
    c=cat.get(r,'other'); col=colors.get(c,'gray')
    x=toa[r]['last10']; y=rm[r]['rmse']
    ax.scatter(x,y,c=col,s=90,edgecolor='k',linewidth=0.5,zorder=3,label=c if c not in seen else None)
    seen.add(c)
    tag=lab(r).replace('06','06').replace('_fesom','').replace('Tuning_test_','')
    # shorten label
    tag=r.replace('Tuning_test_06','06').replace('_fesomA','').replace('_fesomB','').replace('_fesomC','').replace('_fesomE','').replace('_fesomF','').replace('_fesomG','').replace('_fesomH','').replace('_fesomP','')
    ax.annotate(tag, (x,y), fontsize=6.5, xytext=(4,3), textcoords='offset points')
ax.axvline(toa['Tuning_test_06_Baseline']['last10'],color='gray',ls=':',lw=0.8)
ax.axhline(rm['Tuning_test_06_Baseline']['rmse'],color='gray',ls=':',lw=0.8)
ax.set_xlabel('Net TOA radiative imbalance, last decade 1370-1379 [W/m²]  (→ want smaller)')
ax.set_ylabel('2m T RMSE vs ERA5 [K]  (→ want smaller)')
ax.set_title('AWI-ESM3 Tuning_test_06: radiative imbalance vs near-surface temperature error',fontweight='bold')
ax.grid(alpha=0.25); ax.legend(fontsize=8,loc='upper left')
ax.text(0.99,0.01,'dotted lines = baseline; lower-left is better',transform=ax.transAxes,ha='right',fontsize=8,style='italic')
plt.tight_layout(); plt.savefig(R+'/plots/summary_scatter_TOA_vs_T2m.png',dpi=160); plt.close()

# ---- TOA timeseries all runs ----
fig,ax=plt.subplots(figsize=(11,6.5))
for r in runs:
    c=cat.get(r,'other'); col=colors.get(c,'gray')
    lw=2.2 if c=='Baseline' else 1.0
    ax.plot(toa[r]['yr'],toa[r]['v'],color=col,lw=lw,alpha=0.8)
ax.axhline(0,color='k',lw=0.8)
# legend by category
from matplotlib.lines import Line2D
handles=[Line2D([0],[0],color=colors[c],lw=2,label=c) for c in colors]
ax.legend(handles=handles,fontsize=9,loc='upper left')
ax.set_xlabel('Model year'); ax.set_ylabel('Net TOA imbalance [W/m²]')
ax.set_title('Net TOA radiative imbalance, annual global mean (all 22 Tuning_test_06 runs)',fontweight='bold')
ax.grid(alpha=0.25)
plt.tight_layout(); plt.savefig(R+'/plots/summary_TOA_timeseries.png',dpi=160); plt.close()
print("summary plots written")
