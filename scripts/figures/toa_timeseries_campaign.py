"""Annual global energy balance: the PI200 control, the branches off it, and the clean PICAL run.

VAR=toa (default) draws net TOA, VAR=sfc the net surface energy (which includes the snow
enthalpy term); each writes its own <VAR>_timeseries_campaign.png.

Net TOA = tsr + ttr, cos-lat weighted global mean, annual mean of complete years, read from the
campaign store data/coupled_annual_diag.nc (scripts/analysis/coupled_annual_store.py), which also
keeps runs whose raw output has been deleted.  The 15 and 16 tuning series are deliberately out:
this is the production line only.  Context runs in grey: 11E and 11W (core3_beta mesh), 11X
(corrected CORE3, DP, old ice), 11Y (as 11X but single precision with the ocean heat leak), 13A
(DP, old ice).  NOTE the context set differs between the two variables: 11E, 11W and 11Y have
no surface flux output left on disk (their raw output was deleted and only the cached TOA
survives), so VAR=sfc draws 11X and 13A only.  That is missing data, not a selection.
PI200 continues from 16E at 1390; RCL_INPPMIN 70000 -> 50000 (S4) from 1500 and
RCL_INPSEA 0.2 -> 0.1 from 1510, both marked by vertical lines.

The production line runs PICAL -> PICAL_momixoff -> PICAL_ccnice.  PICAL ran 1850-1929 and
ended there; momix-off branched off it at 1920 and continues to 2050; PICAL_ccnice branched
off momix-off at 1940 and is the current line, carrying the sea-ice CCN weighting, S4 removed,
and albsn 0.83 for 1970-1979 then 0.80 from 1980.  The bright-albedo arms PICAL_v35def and
PICAL_v35def_darkice are drawn as context: they are the counter-example, an Arctic volume
runaway that fifty years of correction did not arrest.  PICAL's 1850-1919 are symlinked into
momix-off's outdata
(scripts/model/link_parent_years.sh) so the folder is one continuous record.  PICAL itself is
therefore NOT drawn separately - it would be the same curve up to 1919 and a ten-year overlap
after that.  Model years carry no forcing meaning anywhere here, so the line is drawn from 1350,
the year every other run in this campaign starts at, which puts its first year alongside the
first year of the others.  SHIFT below is that offset and nothing else depends on it.

Panel (a) is the absolute net TOA.  The 20-year branches off PI200 are three sets in three
windows, and at the scale of (a) they sit on top of the control, so (b)-(d) show each branch MINUS
PI200 over the same years instead: the control is then the zero line, which is also how the
branches were judged.  The shaded band is the detection threshold 2*sd*sqrt(2/N) from the
control's own year-to-year scatter in that window (no autocorrelation, so it is a floor).

Usage:  python3 scripts/figures/toa_timeseries_campaign.py
        VAR=sfc python3 scripts/figures/toa_timeseries_campaign.py
"""
import os
import numpy as np, xarray as xr, warnings
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
REPO=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DS=xr.open_dataset(os.path.join(REPO,'data','coupled_annual_diag.nc'))
STORE_RUNS=list(DS['run'].values.astype(str)); YRS=DS['year'].values.astype(int)
DIRS={'11E':'Tuning_test_11E_swemin15_K1'}
# All four PICAL-lineage runs carry the same -500: their own 1850 becomes 1350, which
# puts year 1 alongside year 1 of every other run.  Branch years therefore land in the
# right place relative to their parent without any further bookkeeping.
# Two ways to place the runs on the x-axis.
#
# COMPARISON (CRITICAL_PATH=0): every run shifted so its own year 1 lands on year 1 of
# every other, which is what you want when reading drift rates side by side.  -500 puts
# the PICAL lineage's 1850 at 1350.
#
# LINEAGE (CRITICAL_PATH=1, the default): the runs are laid END TO END in the order the
# campaign actually walked them, each starting where it branched off its parent, so the
# panel reads as one continuous path with no overlap.  The chain and its branch points:
#     PI200          own 1390-1619, unshifted
#     PICAL          took its LPJ-GUESS state from PI200 at 1619-12-31, own 1850 -> 1620
#     PICAL_momixoff branched off PICAL at PICAL's 1920; its record carries PICAL's
#                    1850-1919 symlinked in, so it draws the whole 1850-1940 stretch
#     PICAL_ccnice   branched off momixoff at momixoff's 1940
# All three PICAL-lineage runs share one calendar, so one offset (-230) places them all:
# own 1850 -> 1620, i.e. the year after PI200 ends.
SHIFT_COMPARE={'PICAL_momixoff': -500, 'PICAL_ccnice': -500,
               'PICAL_v35def': -500, 'PICAL_v35def_darkice': -500}
# PI200 is NOT on the path.  PICAL took only its LPJ-GUESS state from PI200 at
# 1619-12-31; FESOM came from core3_linfs at 1957 and the atmosphere was a cold start,
# which is what the -0.8 -> +1.0 excursion over PICAL's first fifteen years is.  A run
# that supplied one component's initial condition is not a predecessor, so the path
# starts at PICAL.  With PI200 gone nothing needs shifting: the lineage is one calendar
# already, and the axis can just be model years.
SHIFT_LINEAGE={}
VAR=os.environ.get('VAR','toa')
# CRITICAL_PATH=1 (the default) draws only the line the campaign actually took:
# PI200 -> PICAL/momixoff -> ccnice.  The grey context runs and the PI200 branch windows
# are exploration, not the path, and at this point they bury the run we care about.
# Set CRITICAL_PATH=0 to get the full figure back.
CRIT=os.environ.get('CRITICAL_PATH','1') not in ('0','false','False')
SHIFT = SHIFT_LINEAGE if CRIT else SHIFT_COMPARE
assert VAR in ('toa','sfc'), VAR
LABEL={'toa':'global net TOA','sfc':'global net surface energy'}[VAR]
TITLE={'toa':'net TOA','sfc':'net surface energy (snow enthalpy included)'}[VAR]

# In critical-path mode a run is drawn only as far as it WAS the path.  PI200 was
# extended past 1559 for one reason, stated in its own runscript: "control for the
# 1560-1579 ice, 1580-1599 GM and 1600-1619 basal-melt branches" - the branch windows
# this mode removes.  Those years are control for exploration, not the line the campaign
# took, and drawn in full they run past the current run and read as the ongoing control.
# Caveat worth keeping in mind: PICAL took its LPJ-GUESS state from PI200 at 1619-12-31,
# so the cut years did contribute the land initial condition even though nothing else
# on the path came through them.
# Each run is drawn only as far as it WAS the path: it stops at the branch point of
# whatever succeeded it.  In lineage coordinates PI200 ends at 1619 where PICAL takes
# over, and momix-off ends at 1710 (its own 1940) where ccnice branches off it; what
# momix-off did after that is a parallel continuation, not the line taken.
# momix-off's record carries PICAL's 1850-1919 symlinked in, so it draws the whole
# 1850-1940 stretch and stops at its own 1940 where ccnice branches off it.
CUT={'PICAL_momixoff': 1940} if CRIT else {}

def series(arm):
    name=DIRS.get(arm,arm)
    if name not in STORE_RUNS: return np.array([]),np.array([])
    v=DS[VAR].sel(run=name).values; k=np.isfinite(v)
    y=YRS[k]+SHIFT.get(arm,0); v=v[k]
    if arm in CUT:
        m=y<=CUT[arm]; y,v=y[m],v[m]
    return y,v

MAIN=[('PI200','PI200  v3.5 LR defaults, 200-yr PI control (S4 from 1500, INPSEA 0.1 from 1510)','#c0302f'),
      ('PICAL_momixoff','PICAL_momixoff  clean PI, calibrated basal melt, momix off from its year 1920 (own years 1850+)','#1f6fb4'),
      ('PICAL_ccnice','PICAL_ccnice  sea-ice CCN, S4 out, albsn 0.80 (own years 1940+); control for the mixing arms','#1baf7a'),
      # Branched off ccnice at 2020 with mix_scheme cvmix_TKE and nothing else changed, so
      # the two lines share every year before 2020 and the divergence after it is the
      # scheme.  ccnice is deliberately NOT cut at the branch point here: it continues as
      # the control and the comparison is the point of the figure.
      ('PICAL_cvTKE','PICAL_cvTKE  as ccnice but cvmix_TKE vertical mixing from 2020','#c07000')]
CTX=[('11E','11E  campaign reference (old mesh)','#b0afaa','-'),
     ('11W','11W  last on old mesh','#b0afaa','--'),
     ('11X','11X  DP, new mesh, old ice','#52514e','-'),
     ('11Y','11Y  SP heat leak, old ice','#52514e','--'),
     ('13A','13A  DP, old ice','#52514e',':')]
# One window per branch set.  Colours are a fixed order, reused across windows because the
# windows do not overlap in time; validated all-pairs against the panel surface (worst CVD
# pair dE 7.1, in the 6-8 band, which is legal only with the direct labels every line carries).
BRANCH_C=['#1f6fb4','#c07000','#00856b','#a03a8f']
WINDOWS=[(1560,1579,'(b) sea-ice and mixed-layer levers',
          [('PI200_spp','spp'),('PI200_mle','MLE'),('PI200_h0','h0min'),('PI200_all3','all three')]),
         (1580,1599,'(c) GM with the Rossby cutoff',
          [('PI200_gmR2500','GM 2500'),('PI200_gmR1500','GM 1500')]),
         (1600,1619,'(d) calibrated ice-shelf basal melt',
          [('PI200_cavg06','gamma 0.6')])]

if CRIT:
    CTX=[]; WINDOWS=[]
    MAIN=[m for m in MAIN if m[0]!='PI200']
data={a:series(a) for a,*_ in MAIN+CTX}
for _,_,_,arms in WINDOWS:
    for a,_ in arms: data[a]=series(a)

arms=[a for a,*_ in CTX+MAIN]
print(f'decadal mean {LABEL} [W/m2]   (n = complete years in the decade; PICAL on its shifted years)')
print('  decade   '+''.join(f'{a:>11}' for a in arms))
for lo in range(1350,int(max(data[a][0].max() for a in arms if len(data[a][0])))+1,10):
    row=''
    for a in arms:
        ys,vs=data[a]; s=(ys>=lo)&(ys<lo+10)
        row+=f'{vs[s].mean():7.2f} ({s.sum():>2})' if s.any() else f'{"":>11}'
    print(f'  {lo}-{(lo+9)%100:02d}  '+row)

# PI200 is the reference for the branch windows only; in critical-path mode it is not
# drawn and the windows are empty, so there is nothing to difference against.
py,pv=data.get('PI200',(np.array([]),np.array([])))
def diff(a,y0,y1):
    """branch minus PI200 on the years both have complete"""
    ys,vs=data[a]
    if not len(ys): return np.array([]),np.array([])
    k=(ys>=y0)&(ys<=y1); ys,vs=ys[k],vs[k]
    m={int(y):v for y,v in zip(py,pv)}
    ok=np.array([int(y) in m for y in ys])
    return ys[ok],vs[ok]-np.array([m[int(y)] for y in ys[ok]])

print(f'\nbranch minus PI200, window mean {LABEL} [W/m2]; threshold = 2*sd*sqrt(2/N) of the control')
for y0,y1,title,brs in WINDOWS:
    s=(py>=y0)&(py<=y1); thr=2*pv[s].std(ddof=1)*np.sqrt(2/max(s.sum(),1))
    print(f'  {y0}-{y1}  control {pv[s].mean():+.3f}, sd {pv[s].std(ddof=1):.3f}, threshold {thr:.3f}')
    for a,lab in brs:
        ys,dv=diff(a,y0,y1)
        if not len(ys): continue
        print(f'    {lab:<26}{dv.mean():+8.3f}  ({"detected" if abs(dv.mean())>thr else "within noise"}, n={len(ys)})')
cy,cv=data['PICAL_momixoff']
if len(cy):
    _sh=SHIFT.get("PICAL_momixoff",0)   # 0 in lineage mode: the axis is already its own calendar
    print(f'\nPICAL_momixoff: {len(cy)} complete years, own years {cy[0]-_sh}-{cy[-1]-_sh}, drawn at {cy[0]}-{cy[-1]}')
    print(f'  {LABEL} per year: '+' '.join(f'{v:+.2f}' for v in cv))


def end_labels(ax,items,pad=1.25):
    """Right-hand direct labels, pushed apart so they never overlap.

    items = [(x, y, text, colour, fontsize, weight)].  Labels are spread in axes
    fraction, nearest-first, keeping at least pad*label height between them; without
    this several runs end the window within a few hundredths of a W/m2 of each other
    and their labels land on top of one another.
    """
    fig=ax.figure; fig.canvas.draw()
    y0,y1=ax.get_ylim(); h=ax.get_window_extent().height
    ent=sorted(items,key=lambda e:e[1])
    gap=[pad*(e[4]*fig.dpi/72.0)/h*(y1-y0) for e in ent]
    ys=[e[1] for e in ent]
    for i in range(1,len(ys)):                       # push up from the bottom
        ys[i]=max(ys[i],ys[i-1]+gap[i])
    over=ys[-1]-(y1-0.02*(y1-y0))
    if over>0:                                       # then slide the stack back down
        for i in range(len(ys)-1,-1,-1):
            ys[i]-=over
            if i and ys[i]-ys[i-1]>gap[i]: break
            over=ys[i-1]-ys[i]+gap[i] if i else 0
    for (x,yv,txt,c,fs,fw),yl in zip(ent,ys):
        ax.annotate(txt,(x,yl),xytext=(6,0),textcoords='offset points',va='center',
                    fontsize=fs,color=c,fontweight=fw,
                    bbox=dict(fc=ax.get_facecolor(),ec='none',pad=1.0))


SURF,TXT1,TXT2,GRID='#fcfcfb','#0b0b0b','#52514e','#e4e3df'
# With no branch windows there is no bottom row to draw, so collapse to one panel
# instead of leaving three empty axes.
if WINDOWS:
    fig=plt.figure(figsize=(15,8.8),dpi=150); fig.patch.set_facecolor(SURF)
    gs=fig.add_gridspec(2,3,height_ratios=[2.05,1.0],hspace=0.34,wspace=0.13,
                        left=0.05,right=0.955,top=0.94,bottom=0.075)
    ax=fig.add_subplot(gs[0,:])
else:
    fig=plt.figure(figsize=(15,6.2),dpi=150); fig.patch.set_facecolor(SURF)
    gs=fig.add_gridspec(1,1,left=0.05,right=0.955,top=0.92,bottom=0.115)
    ax=fig.add_subplot(gs[0,0])
ax.set_facecolor(SURF)
for y0,y1,_,_ in WINDOWS:
    ax.axvspan(y0-0.5,y1+0.5,color='#efeee9',zorder=0)
    ax.axvline(y0-0.5,color=SURF,lw=1.4,zorder=1)
ends=[]
for a,lab,c,ls in CTX:
    ys,vs=data[a]
    if not len(ys): continue
    ax.plot(ys,vs,color=c,lw=1.1,ls=ls,label=lab,zorder=2)
    ends.append((ys[-1],vs[-1],a,TXT2,8.5,'normal'))
# MAIN is drawn in reverse so the FIRST entry ends up on top of the stacking order by
# default; z-order is set explicitly instead, rising with position in MAIN, so the run at
# the end of the list - the current one - is always the most visible line on the panel.
for i,(a,lab,c) in enumerate(MAIN[::-1]):
    ys,vs=data[a]
    if not len(ys): continue
    z=4+(len(MAIN)-1-i)
    last=(i==0)                      # i==0 is the final MAIN entry, i.e. the current run
    ax.plot(ys,vs,color=c,lw=3.0 if last else 2.0,marker='o',ms=5.2 if last else 4.2,
            mec=SURF,mew=1.0,label=lab,zorder=z)
    ends.append((ys[-1],vs[-1],a,TXT1,9.5,'bold'))
ax.axhline(0,color=TXT2,lw=0.9,zorder=1)
for xv,txt,dy in (() if CRIT else ((1499.5,'PI200: RCL_INPPMIN 50000 (S4)',30),(1509.5,'PI200: RCL_INPSEA 0.1',16))):
    ax.axvline(xv,color='#c0302f',lw=0.9,ls='--',zorder=1)
    ax.annotate(txt,(xv-1,ax.get_ylim()[0]),xytext=(0,dy),textcoords='offset points',fontsize=8.5,
                color='#c0302f',va='bottom',ha='right')
# where the production run switches to momix off (its year 1920)
_mx=1920+SHIFT.get('PICAL_momixoff',0)   # own-calendar axis in lineage mode, so no offset
if len(data['PICAL_momixoff'][0]) and data['PICAL_momixoff'][0].max()>=_mx:
    ax.axvline(_mx-0.5,color='#1f6fb4',lw=0.9,ls='--',zorder=1)
    ax.annotate('momix off from here',(_mx-0.5,ax.get_ylim()[1]),xytext=(4,-12),
                textcoords='offset points',fontsize=8.5,color='#1f6fb4',va='top',ha='left')
# where PICAL_ccnice changes FESOM vertical mixing scheme (namelist.oce, oce_dyn).
# KPP to 2059, cvmix_TKE 2060-2079, cvmix_TKE+cvmix_IDEMIX from 2080.  Both are on the
# production line itself, not on a branch, so they belong on this axis.
_cc = data['PICAL_ccnice'][0]
for _yr, _lab in ((2060, 'cvmix_TKE'), (2080, '+ cvmix_IDEMIX')):
    _x = _yr + SHIFT.get('PICAL_ccnice', 0)
    if len(_cc) and _cc.max() >= _x:
        ax.axvline(_x-0.5, color='#1a7f5a', lw=0.9, ls='--', zorder=1)
        ax.annotate(_lab, (_x-0.5, ax.get_ylim()[1]), xytext=(4,-12), textcoords='offset points',
                    fontsize=8.5, color='#1a7f5a', va='top', ha='left')
for y0,y1,title,_ in WINDOWS:
    ax.annotate(title[:3].strip(),(0.5*(y0+y1),ax.get_ylim()[0]),xytext=(0,4),textcoords='offset points',
                ha='center',va='bottom',fontsize=8.5,color=TXT2)
ax.grid(axis='y',color=GRID,lw=0.8); ax.set_axisbelow(True)
for s in ('top','right'): ax.spines[s].set_visible(False)
for s in ('left','bottom'): ax.spines[s].set_color(GRID)
ax.tick_params(colors=TXT2,labelsize=9)
ymax=int(max(data[x][0].max() for x in arms if len(data[x][0])))
_y0=int(min(ys[0] for ys,_ in data.values() if len(ys)))
_start=(_y0//10)*10
# more room on the right in lineage mode: the last run's name is set beyond its
# final point and was being clipped by the axes edge
ax.set_xlim(_start-0.5, ymax+(26 if CRIT else 4.5))
_step=20 if (ymax-_start)>320 else 10
ax.set_xticks(range(_start, ymax+2, _step))
ax.set_ylim(top=ax.get_ylim()[1]+0.75)
end_labels(ax,ends)
ax.set_xlabel(('model year  (the PICAL lineage on its own calendar; each run drawn to where its successor branched off)'
               if CRIT else 'model year  (the PICAL lineage is drawn shifted by -500; its own years are 1850+)'),color=TXT2,fontsize=10)
ax.set_ylabel(LABEL+'  [W m$^{-2}$, positive = gaining]',color=TXT2,fontsize=10)
ax.set_title((f'(a) Annual {TITLE}' if WINDOWS else f'Annual {TITLE}') +
             (': the critical path, momix-off then ccnice, each drawn to the branch point of its successor'
              if CRIT else ': the PI200 control and the PICAL production line, momix-off then ccnice, each from its first year'),
             loc='left',color=TXT1,fontsize=12)
# A run with no data for this VAR draws no line and so has no legend entry: filter MAIN
# the same way CTX already is, or the ordering raises on the first missing series.
h,l=ax.get_legend_handles_labels()
order=[l.index(x[1]) for x in MAIN if x[1] in l]+[l.index(x[1]) for x in CTX if x[1] in l]
ax.legend([h[i] for i in order],[l[i] for i in order],loc='upper left',fontsize=8.3,frameon=False,ncol=2)

axb=None
for j,(y0,y1,title,brs) in enumerate(WINDOWS):
    a2=fig.add_subplot(gs[1,j],sharey=axb) if axb is not None else fig.add_subplot(gs[1,j])
    axb=axb or a2; a2.set_facecolor(SURF)
    s=(py>=y0)&(py<=y1); thr=2*pv[s].std(ddof=1)*np.sqrt(2/max(s.sum(),1))
    a2.axhspan(-thr,thr,color='#efeee9',zorder=0)
    a2.axhline(0,color='#c0302f',lw=1.4,zorder=2)
    bends=[]
    for i,(a,lab) in enumerate(brs):
        ys,dv=diff(a,y0,y1)
        if not len(ys): continue
        c=BRANCH_C[i]
        a2.plot(ys,dv,color=c,lw=1.5,marker='o',ms=3.2,mec=SURF,mew=0.8,label=lab,zorder=3)
        bends.append((ys[-1],dv[-1],lab,c,8.2,'normal'))
    a2.grid(axis='y',color=GRID,lw=0.8); a2.set_axisbelow(True)
    for sp in ('top','right'): a2.spines[sp].set_visible(False)
    for sp in ('left','bottom'): a2.spines[sp].set_color(GRID)
    a2.tick_params(colors=TXT2,labelsize=8.5)
    a2.set_xlim(y0-0.5,y1+6.0); a2.set_xticks(range(y0,y1+1,5))
    end_labels(a2,bends)
    a2.set_title(title,loc='left',color=TXT1,fontsize=10)
    a2.set_xlabel('model year',color=TXT2,fontsize=9)
    if j==0:
        a2.set_ylabel(f'{"net TOA" if VAR=="toa" else "net surface"}, branch minus PI200  [W m$^{-2}$]',color=TXT2,fontsize=9)
        a2.annotate('zero line = PI200 control',(0.02,0.97),xycoords='axes fraction',
                    fontsize=7.8,color='#c0302f',va='top')
        a2.annotate('band = detection threshold',(0.02,0.90),xycoords='axes fraction',
                    fontsize=7.8,color=TXT2,va='top')
    else:
        a2.tick_params(labelleft=False)
out=os.path.join(REPO,'report','plots',f'{VAR}_timeseries_campaign.png'); fig.savefig(out,facecolor=SURF)
alt=os.path.join(REPO,'plots',f'{VAR}_timeseries_campaign.png'); fig.savefig(alt,facecolor=SURF)
print('saved',out); print('saved',alt)
