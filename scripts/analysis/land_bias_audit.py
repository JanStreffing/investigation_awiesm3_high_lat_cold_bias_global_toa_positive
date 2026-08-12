"""How much land cold bias is actually left, per region, PERIOD-CORRECTED.

WHY.  Two levers are being considered: a strong Siberian summer one (+2 K) and a weak
global-land one.  The second is only worth building if the residual is genuinely global
rather than boreal, and every raw number in this campaign is contaminated by the epoch
mismatch -- the tuning runs are 1872-1915 on transient forcing, ERA5 is 1990-2014.

THE CORRECTION, measured not assumed.  amip_presentday is the SAME model and configuration
run over 1989-2015.  So the model's own epoch offset is

    offset(region, season) = presentday(1990-2014) - control(1872-1915)

and the period-clean bias is  raw_bias + offset.  Estimated corrections were used twice in
this campaign and were wrong both times (the boreal one by 0.7 K); this is the direct one.

A CAVEAT THAT MATTERS FOR THE TARGET.  Over Siberia the model warms only +0.42 K between
the epochs where observations imply ~+1.1 K, so the model also under-warms the historical
period.  Correcting with the MODEL's offset therefore gives the smaller (more conservative)
apparent bias; correcting with the OBSERVED epoch change gives a larger one.  Both are
printed, and the honest target is the range between them.

WHY AMIP MAKES THIS A CLEAN LAND TEST.  SST and sea ice are prescribed, so ocean 2 m
temperature is pinned to observations by construction.  The ocean row below is therefore a
METHOD CONTROL, not a result: if it is not near zero, the regridding or the masking is
wrong.  A land bias that survives next to a null ocean row is a land/atmosphere bias.

GRIDDING.  ERA5 0.25 deg is BIN-AVERAGED onto the model grid (area-weighted), not sampled
nearest-neighbour: at 0.25 -> ~1.9 deg a nearest-neighbour pick is one point out of ~55 and
over rough terrain that is a real sampling error, not just noise.
"""
import os, sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runs import RT, LSMF, ERA5_T2M, Y0, Y1

PRESENT = ('amip_presentday', 1990, 2014)
ARMS = [('control', 'amip_pi_base'),
        ('G4  tundra rsmin', 'amip_G4_tundra'),
        ('K1  +land albedo', 'amip_K1_landalb'),
        ('P5  +snow scheme', 'amip_P5_swemin15'),
        ('S4  +INPPMIN 50k', 'amip_S4_inppmin50000')]

# name, lat0, lat1, lon0, lon1, land?
REGIONS = [('Siberia box 55-75N', 55, 75, 60, 180, True),
           ('boreal land 55-70N', 55, 70, 0, 360, True),
           ('Arctic land 70-80N', 70, 80, 0, 360, True),
           ('NH midlat land 30-55N', 30, 55, 0, 360, True),
           ('tropical land 30S-30N', -30, 30, 0, 360, True),
           ('SH land 60-30S', -60, -30, 0, 360, True),
           ('GLOBAL land 60S-75N', -60, 75, 0, 360, True),
           ('[control] ocean 60S-60N', -60, 60, 0, 360, False)]
SEASONS = {'DJF': [11, 0, 1], 'MAM': [2, 3, 4], 'JJA': [5, 6, 7], 'SON': [8, 9, 10],
           'ANN': list(range(12))}
OBS_EPOCH_WARMING = 1.1   # Siberia, observed 1870s -> 1990s (report sub:perdirect)

print(__doc__)
print('=' * 100)

# ---------------------------------------------------------------- model grid + mask
with xr.open_dataset(LSMF, decode_times=False) as d:
    lsm = np.squeeze(d['lsm'].values)
    if lsm.ndim == 3:
        lsm = lsm[0]
    mlat, mlon = d['lat'].values, d['lon'].values
land = lsm > 0.5
# Greenland and Antarctica are permanent ice: their 2 m bias is an ice-sheet problem,
# not a vegetation/snow one, and they would dominate a "global land" mean.
gl = (mlat[:, None] > 59) & (mlat[:, None] < 84) & \
     (((mlon[None, :] > 300) & (mlon[None, :] < 350)))
land = land & ~gl & (mlat[:, None] > -62)
W2 = np.broadcast_to(np.cos(np.deg2rad(mlat))[:, None], lsm.shape)


def model_clim(run, y0, y1):
    acc = []
    for y in range(y0, y1 + 1):
        f = f'{RT}/{run}/outdata/oifs/atm_remapped_1m_2t_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d['2t'].values
        if a.shape[0] == 12:
            acc.append(a)
    if not acc:
        return None, 0
    return np.mean(acc, axis=0) - 273.15, len(acc)


# ---------------------------------------------------------------- ERA5, bin-averaged
with xr.open_dataset(ERA5_T2M, decode_times=False) as d:
    e = d['t2m'].values - 273.15                       # (300, 721, 1440) 1990-2014
    elat, elon = d['latitude'].values, d['longitude'].values
emon = np.stack([e[m::12].mean(axis=0) for m in range(12)])   # (12, 721, 1440)

# area-weighted bin average onto the model grid
def _edges(c):
    c = np.asarray(c, float)
    m = 0.5 * (c[1:] + c[:-1])
    return np.concatenate(([c[0] - (m[0] - c[0])], m, [c[-1] + (c[-1] - m[-1])]))


iy = np.clip(np.searchsorted(np.sort(_edges(mlat)), elat) - 1, 0, len(mlat) - 1)
if mlat[0] > mlat[-1]:
    iy = len(mlat) - 1 - iy
ix = np.clip(np.searchsorted(_edges(mlon), elon % 360) - 1, 0, len(mlon) - 1)
we = np.cos(np.deg2rad(elat))[:, None] * np.ones((1, len(elon)))
IJ = iy[:, None] * len(mlon) + ix[None, :]
den = np.zeros(len(mlat) * len(mlon))
np.add.at(den, IJ.ravel(), we.ravel())
ERA = np.empty((12, len(mlat), len(mlon)))
for m in range(12):
    num = np.zeros_like(den)
    np.add.at(num, IJ.ravel(), (emon[m] * we).ravel())
    ERA[m] = (num / np.where(den > 0, den, np.nan)).reshape(len(mlat), len(mlon))
print(f'ERA5 1990-2014 binned onto the model grid: '
      f'{100 * np.isfinite(ERA[0]).mean():.1f}% of cells filled\n')

# ---------------------------------------------------------------- climatologies
clim = {}
for lab, run in ARMS:
    c, n = model_clim(run, Y0, Y1)
    if c is None:
        print(f'  {lab}: no output, skipped')
        continue
    clim[lab] = c
    print(f'  {lab:20s} {n} years')
pres, npres = model_clim(PRESENT[0], PRESENT[1], PRESENT[2])
print(f'  {"presentday":20s} {npres} years  (epoch reference)\n')


def regmean(field, r):
    name, la0, la1, lo0, lo1, island = r
    sel = (mlat >= la0) & (mlat < la1)
    m = np.zeros_like(land)
    m[np.ix_(sel, (mlon >= lo0) & (mlon <= lo1))] = True
    m &= land if island else (lsm <= 0.5) & (np.abs(mlat)[:, None] < 62)
    m &= np.isfinite(field)
    return float(np.average(field[m], weights=W2[m])), int(m.sum())


# ---------------------------------------------------------------- the audit
for season, mons in SEASONS.items():
    era_s = ERA[mons].mean(axis=0)
    off = {}
    print(f'\n{"="*100}\n{season}   period-clean 2 m temperature bias vs ERA5 1990-2014 [K]'
          f'   (negative = model too cold)\n')
    hdr = f'  {"region":24s} {"cells":>6s} {"offset":>7s} '
    for lab, _ in ARMS:
        if lab in clim:
            hdr += f'{lab.split()[0]:>9s}'
    print(hdr + '   raw(control)')
    print('  ' + '-' * (len(hdr) + 12))
    for r in REGIONS:
        cs = clim['control'][mons].mean(axis=0)
        ps = pres[mons].mean(axis=0)
        eb, n = regmean(cs - era_s, r)
        o, _ = regmean(ps - cs, r)
        off[r[0]] = o
        line = f'  {r[0]:24s} {n:6d} {o:+7.2f} '
        for lab, _ in ARMS:
            if lab not in clim:
                continue
            b, _ = regmean(clim[lab][mons].mean(axis=0) - era_s, r)
            line += f'{b + o:+9.2f}'
        print(line + f'   {eb:+8.2f}')

    if season == 'JJA':
        sib = REGIONS[0]
        cs = clim['control'][mons].mean(axis=0)
        base, _ = regmean(cs - era_s, sib)
        best = max(regmean(clim[l][mons].mean(axis=0) - era_s, sib)[0]
                   for l, _ in ARMS if l in clim)
        print(f'\n  SIBERIA JJA TARGET RANGE, adopted stack:')
        print(f'    corrected with the MODEL epoch offset ({off[sib[0]]:+.2f}): '
              f'{best + off[sib[0]]:+.2f} K')
        print(f'    corrected with the OBSERVED epoch change '
              f'({OBS_EPOCH_WARMING:+.2f}):  {best + OBS_EPOCH_WARMING:+.2f} K')
        print(f'    -> a lever must deliver {abs(best + OBS_EPOCH_WARMING):.2f} to '
              f'{abs(best + off[sib[0]]):.2f} K of Siberian JJA warming.')

print('\n' + '=' * 100)
print('READ THE OCEAN ROW FIRST.  AMIP prescribes SST, so it is a method control and must')
print('be near zero.  If it is, every land number above is a land/atmosphere bias.')
