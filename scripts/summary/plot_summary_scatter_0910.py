"""Figure 2 extended: the 06 family plus rounds 09/10, on a PI-CORRECTED axis.

ERA5 is PRESENT DAY and the runs are PRE-INDUSTRIAL, about 1.1 K apart in global-mean
2m temperature.  Raw RMSD against ERA5 therefore contains an expected offset, and
ranking runs by it penalises exactly the runs that are correctly cold.  Two corrections:

  * y axis is now bias MINUS the -1.1 K PI expectation, so 0 = right for a PI run and
    the sign tells you which way the residual error goes.
  * centred RMSD, sqrt(RMSD^2 - bias^2), removes any uniform offset and isolates the
    spatial PATTERN error.  It is 1.41-1.64 across all 29 runs -- a +-8 % spread
    against a raw-RMSD spread of 64 % -- so the pattern is not what distinguishes
    them and plotting raw RMSD was measuring the mean, twice.


The original covered only Tuning_test_06*.  Rounds 09 and 10 branch from 06V and are
what the campaign now runs, so they belong here -- and putting them on shows the
trade-off the 06-only view could not: the newer runs move LEFT (better radiation) and
UP (worse near-surface temperature).  10A holds the best net TOA in the whole campaign
at +0.565 while its 2m RMSD is 47 % worse than the 06 baseline.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

BASE='/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive'
TOOL='/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview'
OFF=-1.1   # PI minus PD global-mean 2m T
rows=[]
for d in csv.DictReader(open(TOOL+'/data/MASTER_metrics.csv')):
    rows.append((d['run'],float(d['TOA_imbalance_Wm2']),float(d['t2m_bias_K'])-OFF))
new=[]
for d in csv.DictReader(open(BASE+'/data/MASTER_metrics_0910.csv')):
    if d['run'].endswith('_late'): continue          # same runs, longer window
    new.append((d['run'],float(d['TOA_imbalance_Wm2']),float(d['t2m_bias_K'])-OFF))

cat={}
def sc(keys,c):
    for k in keys: cat[k]=c
sc(['Tuning_test_06_Baseline'],'Baseline')
sc(['Tuning_test_06A_fesomA_albpnd028','Tuning_test_06B_fesomB_albpnd028_rfrac075','Tuning_test_06C_fesomC_albpnd028_rfrac075_pndaspect13','Tuning_test_06D_HRlike','Tuning_test_06G_fesomG_stronger_meltpond.','Tuning_test_06L_extra_strong_meltpond','Tuning_test_06H_fesomH_combo_g_rvice018'],'Sea-ice meltpond')
sc(['Tuning_test_06Q_ralbsead0045','Tuning_test_06R_openwater_albedo0075'],'Ocean/ice albedo')
sc(['Tuning_test_06S_evp0'],'Sea-ice dynamics')
sc(['Tuning_test_06E_fesomE_rvice018_ggaussb_m06','Tuning_test_06F_fesomF_F_lrdalb_true','Tuning_test_06P_fesomP_entstpc3_1'],'Atmosphere physics')
sc(['Tuning_test_06N_mospp','Tuning_test_06O_1hcpl_mospp','Tuning_test_06O4_1hcpl_mospp_kv0012','Tuning_test_06O5_1hcpl_mospp_kv002','Tuning_test_06M_1hcpl','Tuning_test_06T_1hcpl_mospp_kpplow','Tuning_test_06U_1hcpl_mospp_kpplow_openwater_albedo0075','Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1'],'Coupling/ocean mixing')
col={'Baseline':'k','Sea-ice meltpond':'#1f77b4','Ocean/ice albedo':'#17becf','Sea-ice dynamics':'#9467bd','Atmosphere physics':'#2ca02c','Coupling/ocean mixing':'#d62728','Rounds 09/10 (veg init + newSeaIce + G4)':'#ff7f0e'}

fig,ax=plt.subplots(figsize=(11.5,7.8)); seen=set()
for run,toa,rmsd in rows:
    c=cat.get(run,'other')
    ax.scatter(toa,rmsd,c=col.get(c,'gray'),s=110,edgecolor='k',lw=0.5,zorder=3,
               label=c if c not in seen else None); seen.add(c)
    tag=run.replace('Tuning_test_06','06')
    for s in ['_fesomA','_fesomB','_fesomC','_fesomE','_fesomF','_fesomG','_fesomH','_fesomP']: tag=tag.replace(s,'')
    ax.annotate(tag,(toa,rmsd),fontsize=6.8,xytext=(4,3),textcoords='offset points')
c='Rounds 09/10 (veg init + newSeaIce + G4)'
for run,toa,rmsd in new:
    ax.scatter(toa,rmsd,c=col[c],s=165,marker='D',edgecolor='k',lw=0.8,zorder=4,
               label=c if c not in seen else None); seen.add(c)
    ax.annotate(run,(toa,rmsd),fontsize=8.2,fontweight='bold',xytext=(5,4),
                textcoords='offset points',zorder=6)

bl=[r for r in rows if r[0]=='Tuning_test_06_Baseline'][0]
ax.axvline(bl[1],color='gray',ls=':',lw=0.8)
ax.axhline(0.0,color='gray',ls='-',lw=1.0)   # 0 = correct for a PI run
ax.axvline(0.0,color='#2a78d6',ls='--',lw=1.0)
ax.annotate('piControl goal 0',(0.0,ax.get_ylim()[1]),xytext=(3,-12),textcoords='offset points',
            fontsize=8,color='#2a78d6',va='top')
ax.set_xlabel('Net TOA radiative imbalance, 30-yr mean (part2) [W/m²]   → smaller = closer to equilibrium')
ax.set_ylabel('2m-T bias vs ERA5 minus the $-$1.1 K PI$-$PD offset [K]\n'
              '→ 0 = correct for a pre-industrial run')
ax.set_title('AWI-ESM3 coupled tuning: radiative imbalance vs PI-corrected 2m-T bias\n'
             'ERA5 is present-day; runs are pre-industrial, so $-$1.1 K of the raw bias is expected, not error',
             fontweight='bold',fontsize=11)
ax.grid(alpha=0.25); ax.legend(fontsize=8.5,loc='lower left',bbox_to_anchor=(0.005,0.09))
ax.scatter([0.0],[0.0],marker='*',s=420,c='none',edgecolor='green',lw=1.6,zorder=8)
ax.annotate('IDEAL\nequilibrium AND correct PI temperature',(0.0,0.0),
            xytext=(14,10),textcoords='offset points',fontsize=9,style='italic',color='green')
ax.annotate('', xy=(0.60,-0.90), xytext=(1.70,0.26),
            arrowprops=dict(arrowstyle='<->',color='#888888',lw=1.0,ls='--'))
ax.text(1.16,-0.34,'the trade-off: every W m$^{-2}$ of TOA\ncosts ~0.6 K of global cold',
        fontsize=8.5,style='italic',color='#666666',rotation=27,ha='center',va='center')
ax.text(0.99,0.02,'centred RMSD (pattern error) is 1.41-1.64 for ALL 29 runs:\n'
        'the spatial pattern does not distinguish them',transform=ax.transAxes,
        fontsize=8,style='italic',color='#555555',ha='right',va='bottom')
plt.tight_layout()
out=BASE+'/report/plots/summary_scatter_with0910_picorr.png'
plt.savefig(out,dpi=160); plt.close()
print('written',out)
