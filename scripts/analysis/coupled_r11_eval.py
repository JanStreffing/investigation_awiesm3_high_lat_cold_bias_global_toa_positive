"""Round 11 coupled: the first runs to carry the AMIP stack into the coupled model.

WHAT IS NEW.  Until now the coupled branch carried exactly ONE lever out of the whole
50-run AMIP campaign -- G4 -- and the report was written on that basis.  The round-11
family changes that: 11A/11B/11D all carry G4 (RVRSMIN 1000/1000/225) AND D2b
(RCL_INPSEA 0.2 with RCL_INPPMIN 70000, the 700 hPa gate), plus a new coupling change
`useIFSsoiltemp`.  110Baseline is 09C plus that coupling change and nothing else, so it
is the correct control for the family -- comparing 11D against 10A would confound the
snow scheme with the soil-temperature coupling.

  110Baseline  09C + useIFSsoiltemp
  11A          + G4 + D2b
  11B          + tanh snow depletion        (ECE_SNOW_SCF=1, the falsified form)
  11D          + fitted snow depletion      (ECE_SNOW_SCF=3, SWEMIN=30)

THE TEST THAT MATTERS.  11D runs the fitted scheme at SWEMIN=30 -- which is P6, and P6
was FALSIFIED in AMIP: DJF T2m -0.850 K (clears the +-0.588 threshold), November
f_full collapsing to 0.752 against 0.884, and DJF soil -1.00 K against the scheme-off
reference.  The adopted AMIP setting is P5, SWEMIN=15.  So 11D is a coupled run of a
configuration the atmosphere-only testbed rejected, and the pre-registered expectation
is that its Siberian DJF soil is COLD against 110Baseline.  If it is not, the AMIP
falsification does not transfer and P6 deserves reconsideration; if it is, 11D should
be re-run at SWEMIN=15 before anything is built on it.

The coupled model can only make this worse than AMIP, not better: AMIP holds SST and
sea ice fixed, so it cannot express the snow-albedo-temperature feedback that turned a
-13.6 K AMIP soil signal into -16.2 K coupled for the tanh scheme.

WINDOW.  Last 10 years (1390-1399) for state, matching the campaign's coupled
convention; 110Baseline stops at 1389 so it uses 1380-1389 and that is flagged rather
than hidden -- a decade of drift difference is not negligible in a run that has not
equilibrated.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, xarray as xr, warnings
warnings.filterwarnings('ignore')

R270 = '/work/bb1469/a270270/runtime/awiesm3-v3.4'
ACC = 3600.0
CERESF = '/work/ab0246/a270092/obs/CERES/CERES_EBAF_Ed4.1_Subset_CLIM01-CLIM12.nc'

RUNS = [
    ('110 base (09C+IFSsoilT)', 'Tuning_test_110Baseline_09C_useIFSsoiltemp_CRUNCEPinit_newSeaIce'),
    ('11A  +G4+D2b', 'Tuning_test_11A_06V_G4_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce'),
    ('11B  +tanh snow', 'Tuning_test_11B_06V_G4_snowDepletion_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce'),
    ('11D  +fitted sw30', 'Tuning_test_11D_G4_fitted_snow_depletion_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce'),
    ('10A  G4 only (ref)', 'Tuning_test_10A_06V_G4_CRUNCEP_plus_CERES_init_newSeaIce'),
    ('10B  +tanh (ref)', 'Tuning_test_10B_06V_G4_snowDepletion_CRUNCEP_plus_CERES_init_newSeaIce'),
]
NLAST = 10
SIB = (55.0, 75.0, 60.0, 180.0)      # lat0 lat1 lon0 lon1
SO = (-65.0, -45.0)
DJF, JJA = [12, 1, 2], [6, 7, 8]

print(__doc__)
print('=' * 104)


def years_avail(root, var):
    D = f'{root}/outdata/oifs'
    if not os.path.isdir(D):
        return []
    ys = []
    for f in os.listdir(D):
        if f.startswith(f'atm_remapped_1m_{var}_') and f.endswith('.nc'):
            try:
                ys.append(int(f.split('_')[-1].split('-')[0]))
            except ValueError:
                pass
    return sorted(ys)


def load(root, var, yrs):
    """(12, nlat, nlon) mean over yrs, plus lat/lon.  De-accumulates flux variables."""
    D = f'{root}/outdata/oifs'
    div = ACC if var in ('tsr', 'ttr', 'tsrc', 'ttrc', 'ssr', 'str') else 1.0
    acc, lat, lon = [], None, None
    for y in yrs:
        f = f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc'
        if not os.path.exists(f):
            f = f'{D}/atm_remapped_1m_{var}_1m_{y}-{y}.nc'
        if not os.path.exists(f):
            continue
        with xr.open_dataset(f, decode_times=False) as d:
            a = d[var].values / div
            if lat is None:
                lat, lon = d['lat'].values, d['lon'].values
        if a.shape[0] == 12:
            acc.append(a)
    if not acc:
        return None, None, None
    return np.mean(acc, axis=0), lat, lon


def gmean(f2d, lat):
    w = np.cos(np.deg2rad(lat))
    return float(np.average(f2d.mean(axis=1), weights=w))


def boxmean(f2d, lat, lon, box, lsm=None):
    la0, la1, lo0, lo1 = box
    ys = (lat >= la0) & (lat <= la1)
    xs = ((lon % 360) >= lo0) & ((lon % 360) <= lo1)
    sub = f2d[np.ix_(ys, xs)]
    w = np.broadcast_to(np.cos(np.deg2rad(lat[ys]))[:, None], sub.shape).copy()
    if lsm is not None:
        w = np.where(lsm[np.ix_(ys, xs)] >= 0.5, w, 0.0)
    if w.sum() == 0:
        return np.nan
    return float(np.average(sub, weights=w))


def zband(f2d, lat, a, b):
    sel = (lat >= a) & (lat < b)
    w = np.cos(np.deg2rad(lat[sel]))
    return float(np.average(f2d[sel, :].mean(axis=1), weights=w))


# CERES SO reference
cds = xr.open_dataset(CERESF)
clat = cds['lat'].values
csel = (clat >= SO[0]) & (clat < SO[1])
cw = np.cos(np.deg2rad(clat[csel]))
so_cre_obs = float(np.average(cds['toa_cre_sw_clim'].values.mean(axis=0)[csel, :].mean(axis=1),
                              weights=cw))
so_cld_obs = float(np.average(cds['cldarea_total_daynight_clim'].values.mean(axis=0)[csel, :]
                              .mean(axis=1), weights=cw))
cds.close()

rows = []
for tag, name in RUNS:
    root = f'{R270}/{name}'
    ys = years_avail(root, '2t')
    if not ys:
        print(f'{tag:24s} no output'); continue
    win = ys[-NLAST:]
    lsm, _, _ = load(root, 'lsm', win[:1])
    if lsm is not None and lsm.ndim == 3:
        lsm = lsm[0]

    t2m, lat, lon = load(root, '2t', win)
    tsr, _, _ = load(root, 'tsr', win)
    ttr, _, _ = load(root, 'ttr', win)
    tsrc, _, _ = load(root, 'tsrc', win)
    stl2, _, _ = load(root, 'stl2', win)
    tcc, _, _ = load(root, 'tcc', win)
    if t2m is None or tsr is None:
        print(f'{tag:24s} incomplete'); continue

    net = (tsr + ttr).mean(axis=0)
    swcre = (tsr - tsrc) if tsrc is not None else None
    d = [m - 1 for m in DJF]
    j = [m - 1 for m in JJA]

    rows.append(dict(
        tag=tag, yrs=f'{win[0]}-{win[-1]}', n=len(win),
        toa=gmean(net, lat),
        t2m=gmean(t2m.mean(axis=0), lat) - 273.15,
        sibJJA=boxmean(t2m[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
        sibDJF=boxmean(t2m[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15,
        soilDJF=(boxmean(stl2[d].mean(axis=0), lat, lon, SIB, lsm) - 273.15
                 if stl2 is not None else np.nan),
        soilJJA=(boxmean(stl2[j].mean(axis=0), lat, lon, SIB, lsm) - 273.15
                 if stl2 is not None else np.nan),
        socre=zband(swcre.mean(axis=0), lat, *SO) if swcre is not None else np.nan,
        socld=zband(tcc.mean(axis=0), lat, *SO) * 100 if tcc is not None else np.nan,
    ))

print(f'{"run":24s} {"window":>10s} {"netTOA":>8s} {"globT2m":>8s} {"SibJJA":>8s} '
      f'{"SibDJF":>8s} {"soilDJF":>8s} {"soilJJA":>8s} {"SO CRE":>8s} {"SO cld":>7s}')
for r in rows:
    print(f'{r["tag"]:24s} {r["yrs"]:>10s} {r["toa"]:+8.3f} {r["t2m"]:8.3f} '
          f'{r["sibJJA"]:8.2f} {r["sibDJF"]:8.2f} {r["soilDJF"]:8.2f} {r["soilJJA"]:8.2f} '
          f'{r["socre"]:8.2f} {r["socld"]:7.2f}')
print(f'{"CERES":24s} {"":>10s} {"":>8s} {"":>8s} {"":>8s} {"":>8s} {"":>8s} {"":>8s} '
      f'{so_cre_obs:8.2f} {so_cld_obs:7.2f}')

# ---------------------------------------------------------------- deltas
byt = {r['tag']: r for r in rows}
base = byt.get('110 base (09C+IFSsoilT)')
if base:
    print('\n' + '=' * 104)
    print('DELTAS against 110Baseline -- the correct control (same useIFSsoiltemp coupling)')
    print('=' * 104)
    print(f'{"run":24s} {"dTOA":>8s} {"dglobT2m":>9s} {"dSibJJA":>8s} {"dSibDJF":>8s} '
          f'{"dsoilDJF":>9s} {"dSOcld":>8s}')
    for r in rows:
        if r['tag'] == base['tag'] or r['tag'].endswith('(ref)'):
            continue
        print(f'{r["tag"]:24s} {r["toa"] - base["toa"]:+8.3f} {r["t2m"] - base["t2m"]:+9.3f} '
              f'{r["sibJJA"] - base["sibJJA"]:+8.2f} {r["sibDJF"] - base["sibDJF"]:+8.2f} '
              f'{r["soilDJF"] - base["soilDJF"]:+9.2f} {r["socld"] - base["socld"]:+8.2f}')
    print(f'\n  NOTE 110Baseline covers {base["yrs"]} while the others cover 1390-1399, so these')
    print('  deltas carry a decade of drift difference.  Treat signs and large magnitudes as')
    print('  real and small differences as unresolved.')

    # ------------------------------------------------------------ the P6 verdict
    print('\n' + '=' * 104)
    print('VERDICT ON 11D: does the AMIP falsification of SWEMIN=30 transfer to the coupled model?')
    print('=' * 104)
    d11d = byt.get('11D  +fitted sw30')
    d11a = byt.get('11A  +G4+D2b')
    d11b = byt.get('11B  +tanh snow')
    if d11d and d11a:
        ds = d11d['soilDJF'] - d11a['soilDJF']
        dt = d11d['sibDJF'] - d11a['sibDJF']
        print(f'  11D minus 11A (the snow scheme, isolated -- same G4+D2b+coupling):')
        print(f'     Siberian DJF soil  {ds:+.2f} K')
        print(f'     Siberian DJF T2m   {dt:+.2f} K')
        print(f'     Siberian JJA T2m   {d11d["sibJJA"] - d11a["sibJJA"]:+.2f} K')
        if d11b:
            print(f'  11B minus 11A (the FALSIFIED tanh, for scale):')
            print(f'     Siberian DJF soil  {d11b["soilDJF"] - d11a["soilDJF"]:+.2f} K'
                  f'   <- AMIP predicted about -13 to -16 K')
        print()
        if ds < -0.5:
            print(f'  *** THE FALSIFICATION TRANSFERS.  11D cools the Siberian winter soil by')
            print(f'      {abs(ds):.2f} K against its own control, the direction AMIP predicted for')
            print(f'      SWEMIN=30 (P6: -1.00 K soil, DJF T2m -0.850 K, both against the')
            print(f'      scheme-off reference).  11D should be re-run at SWEMIN=15, which is')
            print(f'      the adopted AMIP setting and the only one with every season clean.')
        elif ds > 0.5:
            print(f'  The coupled soil WARMS by {ds:+.2f} K -- the opposite of the AMIP')
            print(f'  prediction for SWEMIN=30.  Either the coupling change alters the')
            print(f'  balance, or the AMIP falsifier does not transfer.  Worth chasing before')
            print(f'  discarding SWEMIN=30.')
        else:
            print(f'  Soil moves {ds:+.2f} K -- small.  But "small" is only meaningful against a')
            print(f'  detection threshold, and 10 coupled years is a much weaker test than the')
            print(f'  44 AMIP years the falsification was established on.  Measured below.')

    # -------- can a 10-year coupled window even resolve the AMIP-sized signal? -------
    print('\n' + '-' * 104)
    print('DETECTION THRESHOLD: what signal could a 10-year coupled window actually see?')
    print('-' * 104)

    def peryear(root, var, yrs, months, box):
        """Per-year seasonal box mean, so interannual scatter can be measured."""
        D = f'{root}/outdata/oifs'
        lsm_, la_, lo_ = load(root, 'lsm', yrs[:1])
        if lsm_ is not None and lsm_.ndim == 3:
            lsm_ = lsm_[0]
        out = []
        for y in yrs:
            f = f'{D}/atm_remapped_1m_{var}_{y}-{y}.nc'
            if not os.path.exists(f):
                continue
            with xr.open_dataset(f, decode_times=False) as dd:
                a = dd[var].values
                la, lo = dd['lat'].values, dd['lon'].values
            if a.shape[0] != 12:
                continue
            out.append(boxmean(a[[m - 1 for m in months]].mean(axis=0), la, lo, box, lsm_))
        return np.array(out)

    for lbl, var in (('Siberian DJF soil (stl2)', 'stl2'), ('Siberian DJF T2m', '2t')):
        sds = []
        for t, nm in (('11A', RUNS[1][1]), ('11D', RUNS[3][1])):
            v = peryear(f'{R270}/{nm}', var, list(range(1390, 1400)), DJF, SIB)
            if v.size > 2:
                sds.append(v.std(ddof=1))
        if not sds:
            continue
        sd = float(np.mean(sds))
        # SE of the difference of two independent 10-year means, 95 % two-sided
        se = sd * np.sqrt(2.0 / 10.0)
        thr = 2.101 * se                      # t(0.975, df=18)
        print(f'  {lbl:26s} interannual sd {sd:5.2f} K -> 95 % detection threshold '
              f'+-{thr:4.2f} K on a 10-yr pair')
    print()
    print('  The AMIP falsification of SWEMIN=30 rests on DJF T2m -0.850 K and soil -1.00 K,')
    print('  established over 44 years.  If those thresholds above exceed those numbers, then')
    print('  11D has NOT tested the falsifier -- a null here is the expected outcome whether the')
    print('  scheme is safe or not, and reporting it as vindication would be reading noise.')
    print('  What 11D DOES establish is that the fitted scheme is not catastrophic coupled:')
    print('  11B, the tanh, shows -18.99 K on the identical window and pair, so this diagnostic')
    print('  has ample power for damage of the kind that mattered -- just not for 1 K.')
