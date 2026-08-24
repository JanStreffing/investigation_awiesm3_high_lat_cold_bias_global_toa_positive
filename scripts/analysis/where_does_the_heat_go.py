"""Close the energy budget on the 1850 arms: net TOA sits near +1 W/m2, so where is it going?

THE PUZZLE.  The PI arms hold net TOA near +1 W/m2 for decades while global mean T2m barely
moves and the Gregory slope over the late window is indistinguishable from zero.  Energy
cannot simply vanish, so either it is being stored somewhere the atmosphere does not see, or
the model is not conserving it and one of the two numbers is wrong.

THE TEST.  FESOM writes thetaoga, the volume-mean ocean potential temperature, and volo, the
ocean volume.  Turn the trend in thetaoga into a heat uptake per square metre of the WHOLE
planet and compare it against the TOA imbalance directly:

    OHU  =  rho * cp * V * d(thetaoga)/dt  /  A_earth

If OHU accounts for the imbalance, the answer is that the heat is going into the deep ocean,
the model conserves energy, and the flat Gregory slope is expected rather than alarming: the
surface cannot equilibrate while the interior is still filling.  If OHU is much smaller than
the imbalance, there is a conservation problem and the TOA number cannot be trusted.

Sea-ice melt is included as a check but is expected to be tiny: latent heat of fusion times
the volume trend, over the same area.

CONSTANTS.  rho = 1027 kg/m3, cp = 3996 J/kg/K, L_f = 3.34e5 J/kg, rho_ice = 917 kg/m3,
A_earth = 5.101e14 m2.  Ocean volume is read from volo rather than assumed.

TRAP.  IFS accumulated fluxes divide by 3600.  Sea-ice volume is in 1e9 m3, not m3.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import glob
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
R270 = '/work/bb1469/a270270/runtime/awiesm3-v3.4'
ACC = 3600.0
RHO, CP = 1027.0, 3996.0
LF, RHO_ICE = 3.34e5, 917.0
A_EARTH = 5.101e14
SEC_YR = 365.25 * 86400.0
TRANSIENT_END = 1351


def res(root, pat):
    g = [x for x in glob.glob(f'{root}/{pat}') if os.path.isdir(x)]
    return g[0] if g else None


ARMS = [('06_Baseline untuned', res(R270, 'Tuning_test_06_Baseline')),
        ('11E campaign base',   res(R092, 'Tuning_test_11E_swemin15_K1')),
        ('11N LX4',             res(R092, '*11N*')),
        ('11Q LX4+RSBLB +S4',   res(R092, '11Q')),
        ('11W LX4+RSBLB -S4',   res(R092, '11W'))]


def oifs_year(root, var, y):
    f = f'{root}/outdata/oifs/atm_remapped_1m_{var}_{y}-{y}.nc'
    if not os.path.exists(f):
        return None, None
    with xr.open_dataset(f, decode_times=False) as d:
        n = var if var in d.data_vars else list(d.data_vars)[-1]
        a = d[n].values
        lat = d['lat'].values
    return (a, lat) if a.shape[0] == 12 else (None, None)


def fesom_year(root, var, y):
    g = glob.glob(f'{root}/outdata/fesom/{var}.fesom.{y}.nc')
    if not g:
        return np.nan
    with xr.open_dataset(g[0], decode_times=False) as d:
        n = var if var in d.data_vars else list(d.data_vars)[-1]
        return float(np.mean(d[n].values))


def trend(x, v):
    x = np.asarray(x, float); v = np.asarray(v, float)
    ok = np.isfinite(v); x, v = x[ok], v[ok]
    if len(x) < 6:
        return np.nan, np.nan
    b, a = np.polyfit(x, v, 1)
    r = v - (a + b * x)
    se = np.sqrt((r ** 2).sum() / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
    return b, se


print(__doc__)
print('=' * 104)
hdr = (f'  {"arm":22s} {"yr":>3s} {"netTOA":>8s} {"dT2m/dt":>10s} '
       f'{"dTocn/dt":>11s} {"OHU":>8s} {"ice":>7s} {"OHU+ice":>8s} {"closes?":>9s}')
print(hdr); print('  ' + '-' * (len(hdr) - 2))
print(f'  {"":22s} {"":>3s} {"W/m2":>8s} {"K/century":>10s} '
      f'{"mK/century":>11s} {"W/m2":>8s} {"W/m2":>7s} {"W/m2":>8s} {"":>9s}')

for label, root in ARMS:
    if root is None:
        print(f'  {label:22s} not found'); continue
    ys = sorted(int(f.split('_')[-1].split('-')[0])
                for f in os.listdir(f'{root}/outdata/oifs')
                if f.startswith('atm_remapped_1m_2t_') and f.endswith('.nc'))
    ys = [y for y in ys if y > TRANSIENT_END]
    Y, N, T, TO, IV, VO = [], [], [], [], [], []
    for y in ys:
        t2m, lat = oifs_year(root, '2t', y)
        tsr, _ = oifs_year(root, 'tsr', y)
        ttr, _ = oifs_year(root, 'ttr', y)
        if t2m is None or tsr is None or ttr is None:
            continue
        w = np.cos(np.deg2rad(lat))
        Y.append(y)
        N.append(float(np.average(((tsr + ttr) / ACC).mean(axis=0).mean(axis=1), weights=w)))
        T.append(float(np.average(t2m.mean(axis=0).mean(axis=1), weights=w)))
        TO.append(fesom_year(root, 'thetaoga', y))
        VO.append(fesom_year(root, 'volo', y))
        IV.append(fesom_year(root, 'sivoln', y) + fesom_year(root, 'sivols', y))
    if len(Y) < 8:
        print(f'  {label:22s} too short'); continue
    Y = np.array(Y)
    bT, _ = trend(Y, T)
    bO, seO = trend(Y, TO)
    bI, _ = trend(Y, IV)
    V = np.nanmean(VO)
    ohu = RHO * CP * V * bO / SEC_YR / A_EARTH if np.isfinite(bO) else np.nan
    # sivoln/sivols are 1e9 m3; losing ice RELEASES the latent heat into the system,
    # so a negative volume trend is a positive term on the uptake side
    ice = -bI * 1e9 * RHO_ICE * LF / SEC_YR / A_EARTH if np.isfinite(bI) else 0.0
    tot = ohu + ice
    nm = np.mean(N)
    closes = ('yes' if abs(tot - nm) < 0.25 * abs(nm) else
              'partly' if abs(tot - nm) < 0.5 * abs(nm) else 'NO')
    print(f'  {label:22s} {len(Y):3d} {nm:+8.3f} {bT*100:+10.2f} '
          f'{bO*1e3*100:+11.2f} {ohu:+8.3f} {ice:+7.3f} {tot:+8.3f} {closes:>9s}')

print('\n  OHU = rho*cp*V*d(thetaoga)/dt / A_earth, with V from volo.')
print('  "closes" compares OHU+ice against the mean net TOA over the same years.')
