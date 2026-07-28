# Recovered from the Claude Code session transcript on 2026-07-28.
# Original location was the ephemeral session scratchpad, which was wiped.
# Session c56eada6-56b5-4b15-83cd-c6ed69cf48d9, written 07-27T14:42:12.

import csv, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
SP='/tmp/claude-24456/-work-ab0246-a270092-postprocessing-investigation-awiesm3-high-lat-cold-bias-global-toa-positive/c56eada6-56b5-4b15-83cd-c6ed69cf48d9/'
OUT='/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/plots/'
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTED='#8a8983'
plt.rcParams.update({'font.size':9,'axes.edgecolor':MUTED,'axes.linewidth':0.6,
 'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK2,
 'figure.facecolor':SURF,'axes.facecolor':SURF,'savefig.facecolor':SURF})
# fixed entity -> colour, shared with the seasonal figure
COL={'06 Baseline':'#eb6834','080a (CRUNCEP veg init)':'#1baf7a','06V':'#4a3aa7',
     '06T':'#008300','08B = 06V':'#eda100','08F = 06V + LPJ levers':'#e87ba4'}
ORDER=['06 Baseline','080a (CRUNCEP veg init)','06V','06T','08B = 06V','08F = 06V + LPJ levers']
rows=defaultdict(list)
for r in csv.DictReader(open(SP+'toa_years.csv')): rows[r['label']].append((int(r['year']),float(r['netTOA'])))
def smooth(y,k=5):
    y=np.asarray(y); out=np.full(len(y),np.nan)
    for i in range(len(y)):
        a=max(0,i-k//2); b=min(len(y),i+k//2+1); out[i]=y[a:b].mean()
    return out
fig,ax=plt.subplots(figsize=(9.0,4.6))
ax.axhspan(-0.4,0.1,color=MUTED,alpha=0.10,lw=0)
ax.axhline(-0.16,color=INK2,lw=1.1,ls=(0,(5,3)),zorder=2)
ends=[]
for lab in ORDER:
    xy=sorted(rows[lab]); x=np.array([a for a,_ in xy]); y=np.array([b for _,b in xy])
    c=COL[lab]
    ax.plot(x,y,color=c,lw=0.9,alpha=0.30,zorder=3)
    ys=smooth(y); ax.plot(x,ys,color=c,lw=2.2,solid_capstyle='round',zorder=4)
    ax.plot(x[-1],ys[-1],'o',ms=5.5,color=c,mec=SURF,mew=1.3,zorder=5)
    ends.append([ys[-1],lab,c,x[-1]])
# de-collide end labels
ends.sort(key=lambda e:e[0]); MIN=0.135
for i in range(1,len(ends)):
    if ends[i][0]-ends[i-1][0]<MIN: ends[i][0]=ends[i-1][0]+MIN
for yv,lab,c,xe in ends:
    ax.annotate(lab,xy=(xe,yv),xytext=(xe+1.2,yv),fontsize=8.4,va='center',color=INK2,
                annotation_clip=False)
ax.text(1401,-0.16,'HR piControl goal',fontsize=8.2,color=INK2,va='bottom')
ax.set_xlabel('model year'); ax.set_ylabel('global net TOA imbalance  [W m$^{-2}$]')
ax.grid(axis='y',color=MUTED,alpha=0.22,lw=0.5); ax.set_axisbelow(True)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
ax.set_xlim(1349.5,1400.5); ax.set_ylim(-1.0,2.5)
ax.set_xticks([1350,1360,1370,1380,1390,1400])
fig.suptitle('The radiation objective across the campaign',x=0.008,ha='left',
             fontsize=11.5,fontweight='bold',color=INK,y=1.03)
fig.text(0.008,0.945,'Annual global-mean net TOA (thin) with 5-yr running mean (bold). '
  'Adding an LPJ-GUESS lever moves nothing; only the 06T/06V ocean-mixing branch does — '
  'and it still stops ~1 W m$^{-2}$ above target.',ha='left',fontsize=8.2,color=INK2)
fig.tight_layout(rect=[0,0,0.80,0.90])
fig.savefig(OUT+'campaign_net_toa_by_year.png',dpi=170,bbox_inches='tight')
print('fig2 rewritten')
for lab in ORDER:
    xy=sorted(rows[lab]); y=np.array([b for _,b in xy])
    print(f"{lab:28s} last-decade mean {y[-10:].mean():5.2f}  final-5yr {y[-5:].mean():5.2f}  n={len(y)}")
