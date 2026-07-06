import csv,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R='/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview'
rows=[]
for d in csv.DictReader(open(R+'/data/MASTER_metrics.csv')):
    rows.append((d['run'],float(d['TOA_imbalance_Wm2']),float(d['t2m_RMSD_K']),float(d['t2m_bias_K'])))
cat={}
def sc(keys,c):
    for k in keys: cat[k]=c
sc(['Tuning_test_06_Baseline'],'Baseline')
sc(['Tuning_test_06A_fesomA_albpnd028','Tuning_test_06B_fesomB_albpnd028_rfrac075','Tuning_test_06C_fesomC_albpnd028_rfrac075_pndaspect13','Tuning_test_06D_HRlike','Tuning_test_06G_fesomG_stronger_meltpond.','Tuning_test_06L_extra_strong_meltpond','Tuning_test_06H_fesomH_combo_g_rvice018'],'Sea-ice meltpond')
sc(['Tuning_test_06Q_ralbsead0045','Tuning_test_06R_openwater_albedo0075'],'Ocean/ice albedo')
sc(['Tuning_test_06S_evp0'],'Sea-ice dynamics')
sc(['Tuning_test_06E_fesomE_rvice018_ggaussb_m06','Tuning_test_06F_fesomF_F_lrdalb_true','Tuning_test_06P_fesomP_entstpc3_1'],'Atmosphere physics')
sc(['Tuning_test_06N_mospp','Tuning_test_06O_1hcpl_mospp','Tuning_test_06O4_1hcpl_mospp_kv0012','Tuning_test_06O5_1hcpl_mospp_kv002','Tuning_test_06M_1hcpl','Tuning_test_06T_1hcpl_mospp_kpplow','Tuning_test_06U_1hcpl_mospp_kpplow_openwater_albedo0075','Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1'],'Coupling/ocean mixing')
col={'Baseline':'k','Sea-ice meltpond':'#1f77b4','Ocean/ice albedo':'#17becf','Sea-ice dynamics':'#9467bd','Atmosphere physics':'#2ca02c','Coupling/ocean mixing':'#d62728'}
fig,ax=plt.subplots(figsize=(11,7.5)); seen=set()
for run,toa,rmsd,bias in rows:
    c=cat.get(run,'other'); 
    ax.scatter(toa,rmsd,c=col.get(c,'gray'),s=110,edgecolor='k',lw=0.5,zorder=3,label=c if c not in seen else None); seen.add(c)
    tag=run.replace('Tuning_test_06','06')
    for s in ['_fesomA','_fesomB','_fesomC','_fesomE','_fesomF','_fesomG','_fesomH','_fesomP']: tag=tag.replace(s,'')
    ax.annotate(tag,(toa,rmsd),fontsize=6.8,xytext=(4,3),textcoords='offset points')
bl=[r for r in rows if r[0]=='Tuning_test_06_Baseline'][0]
ax.axvline(bl[1],color='gray',ls=':',lw=0.8); ax.axhline(bl[2],color='gray',ls=':',lw=0.8)
ax.set_xlabel('Net TOA radiative imbalance, 30-yr mean (part2_rad_balance) [W/m²]   → smaller = closer to equilibrium')
ax.set_ylabel('2m-T RMSD vs ERA5 (part8) [K]   → smaller = better')
ax.set_title('AWI-ESM3 Tuning_test_06: radiative imbalance vs near-surface temperature error',fontweight='bold')
ax.grid(alpha=0.25); ax.legend(fontsize=8.5,loc='upper right')
ax.annotate('BETTER\n(equilibrium + low bias)',(min(r[1] for r in rows),min(r[2] for r in rows)),
            fontsize=9,style='italic',color='green',ha='left',va='bottom')
plt.tight_layout(); plt.savefig(R+'/plots/summary_scatter_official.png',dpi=160); plt.close()
print("written summary_scatter_official.png")
