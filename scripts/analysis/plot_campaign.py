# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-27T14:39:21.

import json, csv, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

SP='/tmp/claude-24456/-work-ab0246-a270092-postprocessing-investigation-awiesm3-high-lat-cold-bias-global-toa-positive/c56eada6-56b5-4b15-83cd-c6ed69cf48d9/'
OUT='/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/plots/'
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTED='#8a8983'
C=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300']
plt.rcParams.update({'font.size':9,'axes.edgecolor':MUTED,'axes.linewidth':0.6,
    'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK2,
    'figure.facecolor':SURF,'axes.facecolor':SURF,'savefig.facecolor':SURF})
MON=['J','F','M','A','M','J','J','A','S','O','N','D']

# ---------------- FIGURE 1 : seasonal bias vs CRUNCEP3 ----------------
d=json.load(open(SP+'scratchpad/seasonal.json'))
ref=d['CRUNCEP3 reference']
series=['AMIP (atmosphere only, no ocean/no LPJG)','06 Baseline','080a (CRUNCEP veg init)',
        '08B = 06V','08F = 06V + LPJ levers']
boxes=['Boreal 55-70N','Siberia','E. Siberia']
fig,axes=plt.subplots(1,3,figsize=(11.6,3.9),sharey=True)
tbl={}
for j,(ax,b) in enumerate(zip(axes,boxes)):
    ax.axhspan(-0.5,0.5,color=MUTED,alpha=0.10,lw=0)
    ax.axhline(0,color=MUTED,lw=0.8,zorder=1)
    for i,s in enumerate(series):
        y=np.array(d[s][b])-np.array(ref[b])
        tbl.setdefault(s,{})[b]=y
        ax.plot(range(1,13),y,color=C[i],lw=2.0,solid_capstyle='round',zorder=3,
                marker='o',ms=3.2,mec=SURF,mew=0.8)
    ax.set_xticks(range(1,13)); ax.set_xticklabels(MON)
    ax.set_title(b,fontsize=9.5,color=INK,pad=6,loc='left',fontweight='bold')
    ax.grid(axis='y',color=MUTED,alpha=0.22,lw=0.5); ax.set_axisbelow(True)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    ax.set_xlim(0.6,12.4)
    if j==0: ax.set_ylabel('2 m T bias vs CRUNCEP3  [K]')
    ax.axvspan(5.5,8.5,color=MUTED,alpha=0.07,lw=0)
axes[1].text(7.0,axes[0].get_ylim()[1]*0.92,'JJA',ha='center',fontsize=8,color=MUTED)
fig.suptitle('Boreal 2 m-temperature bias by month — the atmosphere alone vs the coupled campaign',
             x=0.011,ha='left',fontsize=11.5,fontweight='bold',color=INK,y=1.015)
fig.text(0.011,0.925,'model $-$ CRUNCEP3 (1901–1910); land only, cos-lat weighted. '
         'Coupled runs: model yrs 1370–1379. AMIP: 1870–1879, prescribed SST, prescribed vegetation.',
         ha='left',fontsize=8.2,color=INK2)
h=[plt.Line2D([],[],color=C[i],lw=2.0) for i in range(len(series))]
fig.legend(h,series,loc='lower center',ncol=5,frameon=False,fontsize=8.2,
           bbox_to_anchor=(0.5,-0.085),labelcolor=INK2,handlelength=1.6,columnspacing=1.6)
fig.tight_layout(rect=[0,0.02,1,0.90])
fig.savefig(OUT+'campaign_boreal_seasonal_bias.png',dpi=170,bbox_inches='tight')
print('fig1 written')
for s in series:
    print(f"{s:46s} "+" | ".join(f"{b}: JJA {tbl[s][b][5:8].mean():+.2f} ANN {tbl[s][b].mean():+.2f}" for b in boxes))

# ---------------- FIGURE 2 : net TOA per model year ----------------
rows=defaultdict(list)
for r in csv.DictReader(open(SP+'toa_years.csv')):
    rows[r['label']].append((int(r['year']),float(r['netTOA'])))
order=['06 Baseline','06T','06V','07A (LPJ lever only)','080a (CRUNCEP veg init)','08B = 06V']
fig,ax=plt.subplots(figsize=(8.4,4.4))
ax.axhline(-0.16,color=MUTED,lw=1.2,ls=(0,(5,3)),zorder=2)
ax.text(1399.4,-0.16,' HR piControl goal',va='center',fontsize=8,color=INK2)
for i,lab in enumerate(order):
    xy=sorted(rows[lab]); x=[a for a,_ in xy]; y=[b for _,b in xy]
    ax.plot(x,y,color=C[i],lw=2.0,solid_capstyle='round',zorder=3)
    ax.plot(x[-1],y[-1],'o',ms=5,color=C[i],mec=SURF,mew=1.2,zorder=4)
    ax.text(x[-1]+0.6,y[-1],lab,fontsize=8.2,va='center',color=INK2)
ax.set_xlabel('model year'); ax.set_ylabel('global net TOA imbalance  [W m$^{-2}$]')
ax.grid(axis='y',color=MUTED,alpha=0.22,lw=0.5); ax.set_axisbelow(True)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.set_xlim(1349,1414)
fig.suptitle('The radiation objective across the campaign',x=0.011,ha='left',
             fontsize=11.5,fontweight='bold',color=INK,y=1.02)
fig.text(0.011,0.945,'Annual global-mean net TOA (tsr+ttr). The LPJ-GUESS levers are radiatively inert; '
         'only the 06T/06V ocean-mixing branch moves the budget — and it still stops ~1 W m$^{-2}$ short.',
         ha='left',fontsize=8.2,color=INK2)
fig.tight_layout(rect=[0,0,1,0.91])
fig.savefig(OUT+'campaign_net_toa_by_year.png',dpi=170,bbox_inches='tight')
print('fig2 written')
