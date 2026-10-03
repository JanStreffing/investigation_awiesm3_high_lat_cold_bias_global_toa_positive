"""Annual global net TOA for the 15F baseline and the 16 series (ocean heat-uptake levers).

Net TOA = tsr + ttr (accumulated J/m2 per hourly step, /3600), cos-lat weighted global mean
of the monthly remapped output, annual mean; only years with all 12 months are used.
Series colours: dataviz reference palette, validated (validate_palette.js, light).

Usage:  python3 scripts/figures/toa_timeseries_16series.py
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

ARMS=[('15F','15F  baseline (fixed SP ocean)','#2a78d6'),
      ('16A','16A  Redi 300','#eb6834'),
      ('16B','16B  Fox-Kemper MLE','#1baf7a'),
      ('16C','16C  GM 1500','#eda100')]
data={a:series(a) for a,*_ in ARMS}

print('annual global net TOA [W/m2]')
yrs=sorted({int(y) for a in data for y in data[a][0]})
print('  year '+''.join(f'{a:>8}' for a,*_ in ARMS))
for y in yrs:
    row=''
    for a,*_ in ARMS:
        ys,vs=data[a]; i=np.where(ys==y)[0]
        row+=f'{vs[i[0]]:8.2f}' if len(i) else f'{"":>8}'
    print(f'  {y} '+row)
for lo,hi in ((1350,1354),(1355,1359),(1360,1369)):
    row=''
    for a,*_ in ARMS:
        ys,vs=data[a]; s=(ys>=lo)&(ys<=hi)
        row+=f'{vs[s].mean():8.2f}' if s.any() else f'{"":>8}'
    print(f'  {lo}-{hi % 100:02d}'[:10].ljust(7)+row)

SURF,TXT1,TXT2,GRID='#fcfcfb','#0b0b0b','#52514e','#e4e3df'
fig,ax=plt.subplots(figsize=(10,5),dpi=150); fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
for a,lab,c in ARMS[::-1]:
    ys,vs=data[a]; base=(a=='15F')
    ax.plot(ys,vs,color=c,lw=2.2 if base else 1.6,marker='o',ms=5 if base else 4,
            mec=SURF,mew=1.0,label=lab,zorder=4 if base else 3)
    ax.annotate(a,(ys[-1],vs[-1]),xytext=(7,0),textcoords='offset points',va='center',
                fontsize=9.5,color=TXT1,fontweight='bold' if base else 'normal')
ax.axhline(0,color=TXT2,lw=0.9,zorder=1)
ax.grid(axis='y',color=GRID,lw=0.8); ax.set_axisbelow(True)
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_color(GRID)
ax.tick_params(colors=TXT2,labelsize=9)
ymax=max(int(data[x][0].max()) for x,*_ in ARMS)
ax.set_xlim(1349.5,ymax+2.5); ax.set_xticks(range(1350,ymax+2,2))
ax.set_xlabel('model year',color=TXT2,fontsize=10)
ax.set_ylabel('global net TOA  [W m$^{-2}$, positive = gaining]',color=TXT2,fontsize=10)
ax.set_title('Annual global net TOA: 15F baseline and the 16 ocean-mixing arms',loc='left',color=TXT1,fontsize=12)
h,l=ax.get_legend_handles_labels(); order=[l.index(x[1]) for x in ARMS]
ax.legend([h[i] for i in order],[l[i] for i in order],loc='upper left',fontsize=9,frameon=False)
fig.tight_layout()
out=os.path.join(REPO,'report','plots','toa_timeseries_16series.png'); fig.savefig(out,facecolor=SURF)
print('saved',out)
