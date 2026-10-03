"""Annual global net TOA for the 15 series, with 11X, 11Y and 13A as references.

Net TOA = tsr + ttr (accumulated J/m2 per hourly step, /3600), cos-lat weighted global mean
of the monthly remapped output, annual mean.  Only years with all 12 months are used; 15B
and 15C read from their work dirs because their output was never moved to outdata.
Series colours: dataviz reference palette, validated (validate_palette.js, light).
Grey references: 11X (DP, old ice), 11Y (SP with the heat leak), 13A (DP, old ice).

Usage:  python3 scripts/figures/toa_timeseries_15series.py
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
R='/work/bb1469/a270092/runtime/awiesm3-v3.4'; ACC=3600.0
REPO=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LSMF=('/work/bb1469/a270270/runtime/awiesm3-v3.4/Tuning_test_08B_06V_06Tplus_ENTSTPC3_CRUNCEPinit/'
      'outdata/oifs/atm_remapped_1m_lsm_1350-1350.nc')
with xr.open_dataset(LSMF,decode_times=False) as d:
    lat=np.squeeze(d['lat'].values); m=np.squeeze(d['lsm'].values); m=m[0] if m.ndim==3 else m
W=np.cos(np.deg2rad(np.broadcast_to(lat[:,None],m.shape)))

def find(arm,v,y):
    for p in (f'{R}/{arm}/outdata/oifs/atm_remapped_1m_{v}_{y}-{y}.nc',
              *sorted(glob.glob(f'{R}/{arm}/run_*/work/atm_remapped_1m_{v}_{y}-{y}.nc'))):
        if os.path.isfile(p): return p
def ld(p):
    with xr.open_dataset(p,decode_times=False) as d:
        k=[c for c in d.data_vars if 'bnds' not in c and 'bounds' not in c][0]; return np.squeeze(d[k].values)
def series(arm,y0=1350,y1=1389):
    ys,vs=[],[]
    for y in range(y0,y1+1):
        a,b=find(arm,'tsr',y),find(arm,'ttr',y)
        if not (a and b): continue
        s,t=ld(a),ld(b)
        if s.shape[0]!=12 or t.shape[0]!=12: continue
        f=(s+t).mean(0)/ACC; ys.append(y); vs.append(float((f*W).sum()/W.sum()))
    return np.array(ys),np.array(vs)

MAIN=[('15F','15F  fixed SP ocean (baseline)','#2a78d6'),
      ('15C','15C  unfixed SP','#eb6834'),
      ('15B','15B  unfixed SP','#1baf7a'),
      ('15A','15A  unfixed SP','#eda100')]
REF=[('11X','11X  DP, old ice','#52514e','-'),
     ('11Y','11Y  SP leak, old ice','#8a8985','--'),
     ('13A','13A  DP, old ice','#52514e',':')]
data={a:series(a) for a,*_ in MAIN+REF}

print('annual global net TOA [W/m2]')
yrs=sorted({int(y) for a in data for y in data[a][0]})
arms=[a for a,*_ in MAIN+REF]
print('  year '+''.join(f'{a:>8}' for a in arms))
for y in yrs:
    row=''
    for a in arms:
        ys,vs=data[a]; i=np.where(ys==y)[0]
        row+=f'{vs[i[0]]:8.2f}' if len(i) else f'{"":>8}'
    print(f'  {y} '+row)
print('  mean '+''.join(f'{data[a][1].mean():8.2f}' for a in arms))

SURF,TXT1,TXT2,GRID='#fcfcfb','#0b0b0b','#52514e','#e4e3df'
fig,ax=plt.subplots(figsize=(11,5.2),dpi=150); fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
for a,lab,c,ls in REF:
    ys,vs=data[a]
    ax.plot(ys,vs,color=c,lw=1.2,ls=ls,label=lab,zorder=2)
    ax.annotate(a,(ys[-1],vs[-1]),xytext=(6,0),textcoords='offset points',va='center',fontsize=9,color=TXT2)
for a,lab,c in MAIN[::-1]:
    ys,vs=data[a]
    ax.plot(ys,vs,color=c,lw=2.2 if a=='15F' else 1.6,marker='o',ms=4.5 if a=='15F' else 3.5,
            mec=SURF,mew=1.0,label=lab,zorder=4 if a=='15F' else 3)
    ax.annotate(a,(ys[-1],vs[-1]),xytext=(6,0),textcoords='offset points',va='center',fontsize=9.5,
                color=TXT1,fontweight='bold' if a=='15F' else 'normal')
ax.axhline(0,color=TXT2,lw=0.9,zorder=1)
ax.grid(axis='y',color=GRID,lw=0.8); ax.set_axisbelow(True)
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_color(GRID)
ax.tick_params(colors=TXT2,labelsize=9)
ax.set_xlim(1349.5,1392); ax.set_xlabel('model year',color=TXT2,fontsize=10)
ax.set_ylabel('global net TOA  [W m$^{-2}$, positive = gaining]',color=TXT2,fontsize=10)
ax.set_title('Annual global net TOA: 15 series against 11X, 11Y and 13A',loc='left',color=TXT1,fontsize=12)
h,l=ax.get_legend_handles_labels(); order=[l.index(x[1]) for x in MAIN]+[l.index(x[1]) for x in REF]
ax.legend([h[i] for i in order],[l[i] for i in order],loc='upper left',fontsize=8.5,frameon=False,ncol=2)
fig.tight_layout()
out=os.path.join(REPO,'report','plots','toa_timeseries_15series.png'); fig.savefig(out,facecolor=SURF)
print('saved',out)
