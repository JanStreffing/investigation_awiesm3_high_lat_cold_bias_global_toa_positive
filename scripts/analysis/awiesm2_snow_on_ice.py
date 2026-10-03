"""Snow on sea ice through the year: AWI-ESM2 PI against AWI-ESM3 11X and 13A.

WHY.  13A keeps 0.27-0.32 m of snow on Arctic ice through July-August, so its ice albedo
never leaves the snow branch. AWI-ESM2 runs the same FESOM thermodynamics but under ECHAM,
where FESOM compiles the `#else` branch of the snow-melt block: snow melts on energy alone,
with no `t > 273 K` skin gate. If AWI-ESM2 clears its snow in summer, the gate is implicated.

AWI-ESM2 reference: runtime/AWIESM2_pict, echam output linked from
/work/ba1066/a270107/esm_tools/EXP/PI_wisofix_c (run by a colleague, not in this campaign).

TRAPS.  m_snow is volume per grid area: divide by a_ice for snow depth on the floe.
a_ice is DAILY in the AWI-ESM3 output; aggregate to months before combining with m_snow.
AWI-ESM2 is on its own mesh; coordinates come from the file if present, else nod2d.out.

Usage:  python3 scripts/analysis/awiesm2_snow_on_ice.py
"""
import os, glob
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(v,'1')
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

A2   = '/work/ab0246/a270092/runtime/AWIESM2_pict/outdata/fesom'
R3   = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
_ML  = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
MN   = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

def to_monthly(a):
    nt=a.shape[0]
    if nt==12: return a
    ml=_ML.copy(); ml[1]=29 if nt==366 else 28
    o=np.empty((12,)+a.shape[1:]); i=0
    for m in range(12): o[m]=np.nanmean(a[i:i+ml[m]],axis=0); i+=ml[m]
    return o

def lat_of(ds, fallback):
    for k in ('lat','latitude','nod2_lat'):
        if k in ds.variables: return np.asarray(ds[k].values)
    return fallback

A2_LAT=None
with xr.open_dataset(f'{A2}/a_ice.fesom.6000.nc',decode_times=False) as d:
    print('AWI-ESM2 a_ice file:', dict(d.sizes), list(d.variables)[:10])
    A2_LAT=lat_of(d,None)
if A2_LAT is None:
    cands=glob.glob('/work/ab0246/a270092/runtime/AWIESM2_pict/**/namelist.config',recursive=True)
    mp=None
    for c in cands:
        for line in open(c):
            if 'meshpath' in line.lower(): mp=line.split('=')[1].strip().strip("'\",")
    print('mesh path from namelist.config:', mp)
    if mp and os.path.exists(os.path.join(mp,'nod2d.out')):
        with open(os.path.join(mp,'nod2d.out')) as f:
            n=int(f.readline())
            A2_LAT=np.loadtxt(f,max_rows=n,usecols=2)
        print('AWI-ESM2 nodes from nod2d.out:', A2_LAT.size)
if A2_LAT is None: raise SystemExit('no coordinates for the AWI-ESM2 mesh')

def cycle(path_fmt, years, latsel, lat_fallback=None):
    rows=[]
    for y in years:
        pa,ps=path_fmt.format(v='a_ice',y=y), path_fmt.format(v='m_snow',y=y)
        if not (os.path.exists(pa) and os.path.exists(ps)): continue
        with xr.open_dataset(pa,decode_times=False) as da, xr.open_dataset(ps,decode_times=False) as ds:
            ai=to_monthly(np.asarray(da['a_ice'].values,dtype=np.float64))
            ms=to_monthly(np.asarray(ds['m_snow'].values,dtype=np.float64))
            lat=lat_of(da,lat_fallback)
        m=latsel(lat); ai,ms=ai[:,m],ms[:,m]
        r=[]
        for k in range(12):
            s=(ai[k]>=0.15)&np.isfinite(ms[k])
            r.append(float((ms[k][s]/np.maximum(ai[k][s],1e-6)).mean()) if s.any() else np.nan)
        rows.append(r)
    return np.nanmean(np.array(rows),0) if rows else None

for hemi,sel in (('NH 60-90N',lambda l:l>=60),('SH 90-60S',lambda l:l<=-60)):
    print(f'\nsnow depth on ice (m_snow/a_ice, m), {hemi}, a_ice>=0.15   [FESOM albedo step at 1 mm]')
    print(f"{'run':>9} "+' '.join(f'{m:>6}' for m in MN))
    for name,fmt,yrs,fb in (('AWI-ESM2',A2+'/{v}.fesom.{y}.nc',(5998,5999,6000),A2_LAT),
                            ('11X',R3+'/11X/outdata/fesom/{v}.fesom.{y}.nc',(1357,1358,1359),None),
                            ('13A',R3+'/13A/outdata/fesom/{v}.fesom.{y}.nc',(1357,1358,1359),None)):
        c=cycle(fmt,yrs,sel,fb)
        print(f'{name:>9} '+(' '.join(f'{x:6.3f}' for x in c) if c is not None else 'missing'))
