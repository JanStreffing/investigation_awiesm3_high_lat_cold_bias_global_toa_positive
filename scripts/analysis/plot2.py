# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-28T07:39:33.

import json,csv,numpy as np,matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt; from collections import defaultdict
SP='/tmp/claude-24456/-work-ab0246-a270092-postprocessing-investigation-awiesm3-high-lat-cold-bias-global-toa-positive/c56eada6-56b5-4b15-83cd-c6ed69cf48d9/'
OUT='/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/plots/'
DAT='/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/data/'
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTED='#8a8983'
S={'AMIP baseline (new build)':'#2a78d6','080a (CRUNCEP veg init)':'#eb6834',
   '09A = baseline + newSeaIce':'#1baf7a','09B = 06T + newSeaIce':'#eda100',
   '09C = 06V + newSeaIce':'#e87ba4'}
LBL={'AMIP baseline (new build)':'AMIP  (atmosphere only)','080a (CRUNCEP veg init)':'080a  baseline, old sea ice',
     '09A = baseline + newSeaIce':'09A  baseline + newSeaIce','09B = 06T + newSeaIce':'09B  06T + newSeaIce',
     '09C = 06V + newSeaIce':'09C  06V + newSeaIce'}
plt.rcParams.update({'font.size':9,'axes.edgecolor':MUTED,'axes.linewidth':0.6,'xtick.color':INK2,
 'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK2,'figure.facecolor':SURF,
 'axes.facecolor':SURF,'savefig.facecolor':SURF})
MON=['J','F','M','A','M','J','J','A','S','O','N','D']
d=json.load(open(SP+'scratchpad/seasonal.json')); ref=d['CRUNCEP3 reference']
order=list(S.keys()); boxes=['Boreal 55-70N','Siberia','E. Siberia']
fig,axes=plt.subplots(1,3,figsize=(11.6,4.0),sharey=True)
for j,(ax,b) in enumerate(zip(axes,boxes)):
    ax.axvspan(5.5,8.5,color=MUTED,alpha=0.07,lw=0); ax.axhline(0,color=MUTED,lw=0.8,zorder=1)
    for s in order:
        y=np.array(d[s][b])-np.array(ref[b])
        ax.plot(range(1,13),y,color=S[s],lw=2.0,solid_capstyle='round',zorder=3,
                marker='o',ms=3.2,mec=SURF,mew=0.8)
    ax.set_xticks(range(1,13)); ax.set_xticklabels(MON)
    ax.set_title(b,fontsize=9.5,color=INK,pad=6,loc='left',fontweight='bold')
    ax.grid(axis='y',color=MUTED,alpha=0.22,lw=0.5); ax.set_axisbelow(True); ax.set_xlim(0.6,12.4)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    if j==0: ax.set_ylabel('2 m T bias vs CRUNCEP3  [K]')
fig.suptitle('Round 09 against the atmosphere-only floor — boreal 2 m-T bias by month',
             x=0.011,ha='left',fontsize=11.5,fontweight='bold',color=INK,y=1.015)
fig.text(0.011,0.93,'model $-$ CRUNCEP3 (1901–1910); land only, cos-lat weighted. Coupled: model yrs 1370–79. '
 'AMIP: 1872–79, prescribed SST + prescribed vegetation, same OIFS branch as round 09. Shaded band = JJA.',
 ha='left',fontsize=8.2,color=INK2)
h=[plt.Line2D([],[],color=S[s],lw=2.0) for s in order]
fig.legend(h,[LBL[s] for s in order],loc='lower center',ncol=5,frameon=False,fontsize=8.2,
           bbox_to_anchor=(0.5,-0.075),labelcolor=INK2,handlelength=1.6,columnspacing=1.4)
fig.tight_layout(rect=[0,0.02,1,0.90])
fig.savefig(OUT+'campaign_boreal_seasonal_bias.png',dpi=170,bbox_inches='tight')
with open(DAT+'campaign_boreal_seasonal_bias.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['series','box']+[f'm{i}' for i in range(1,13)])
    for s in d:
        for b in d[s]: w.writerow([s,b]+[round(v,3) for v in d[s][b]])
print('fig1 ok')

# ---- FIGURE 2 ----
rows=defaultdict(list)
for r in csv.DictReader(open(SP+'toa_years.csv')):
    try: rows[r['label']].append((int(r['year']),float(r['netTOA'])))
    except: pass
AMIP_FLOOR=0.6684
o2=['080a (CRUNCEP veg init)','09A = baseline + newSeaIce','09B = 06T + newSeaIce','09C = 06V + newSeaIce']
def sm(y,k=5):
    y=np.asarray(y); r=np.full(len(y),np.nan)
    for i in range(len(y)): a=max(0,i-k//2); b=min(len(y),i+k//2+1); r[i]=y[a:b].mean()
    return r
fig,ax=plt.subplots(figsize=(10.6,4.7))
fig.subplots_adjust(left=0.085,right=0.735,top=0.88,bottom=0.13)
ax.axhline(-0.16,color=INK2,lw=1.1,ls=(0,(5,3)),zorder=2)
ax.axhline(AMIP_FLOOR,color='#2a78d6',lw=1.4,ls=(0,(1.5,2)),zorder=2)
ax.axhspan(-0.16,AMIP_FLOOR,color='#2a78d6',alpha=0.06,lw=0,zorder=1)
ends=[]
for lab in o2:
    xy=[(a,b) for a,b in sorted(rows[lab]) if a<=1379]
    x=np.array([a for a,_ in xy]); y=np.array([b for _,b in xy]); c=S[lab]
    ax.plot(x,y,color=c,lw=0.9,alpha=0.30,zorder=3)
    ys=sm(y); ax.plot(x,ys,color=c,lw=2.2,solid_capstyle='round',zorder=4)
    ax.plot(x[-1],ys[-1],'o',ms=5.5,color=c,mec=SURF,mew=1.3,zorder=5)
    ends.append([ys[-1],LBL[lab],c,x[-1]])
ends.sort(key=lambda e:e[0]); MINS=0.16
for i in range(1,len(ends)):
    if ends[i][0]-ends[i-1][0]<MINS: ends[i][0]=ends[i-1][0]+MINS
for yv,lab,c,xe in ends:
    ax.annotate(lab,xy=(1.025,yv),xycoords=('axes fraction','data'),fontsize=8.4,va='center',
                color=c,annotation_clip=False)
ax.text(1352,AMIP_FLOOR+0.06,'atmosphere-only floor (AMIP, +0.67)',fontsize=8.2,color='#2a78d6')
ax.text(1352,-0.16+0.06,'HR piControl goal ($-$0.16)',fontsize=8.2,color=INK2)
ax.set_xlabel('model year'); ax.set_ylabel('global net TOA imbalance  [W m$^{-2}$]')
ax.grid(axis='y',color=MUTED,alpha=0.22,lw=0.5); ax.set_axisbelow(True)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.set_xlim(1349.5,1380.5); ax.set_ylim(-1.0,2.5)
fig.suptitle('The radiation objective — round 09 vs the atmosphere-only floor',x=0.008,ha='left',
             fontsize=11.5,fontweight='bold',color=INK,y=0.985)
fig.text(0.008,0.925,'Annual global-mean net TOA (thin) with 5-yr running mean (bold). The shaded band is the part of '
 'the gap ocean-side tuning could still close.',
 ha='left',fontsize=8.2,color=INK2)
fig.savefig(OUT+'campaign_net_toa_by_year.png',dpi=170)
import shutil; shutil.copy(SP+'toa_years.csv',DAT+'campaign_net_toa_by_year.csv')
print('fig2 ok')
for lab in o2:
    y=np.array([b for _,b in sorted(rows[lab])]); print(f"{lab:30s} last-decade {y[-10:].mean():.3f}")
