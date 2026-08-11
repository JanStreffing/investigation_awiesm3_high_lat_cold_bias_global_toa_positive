"""Is the boreal forest actually there?  PFTs and the growing season, coupled.

WHY THIS AND NOT A TEMPERATURE NUMBER.  The campaign's stated goal was never a T2m
target: it was a NH high-latitude boreal-summer cold bias "severe enough to lose the
boreal forest".  Round 27's coupled evaluation showed 11E is the first arm to warm
Siberian summer without a winter penalty, and I called Siberia "essentially solved".
That claim is about DIFFERENCES BETWEEN ARMS.  It says nothing about whether the
vegetation is right, which is the thing that actually has to be true.

THE TWO QUESTIONS, taken directly:
  1. Do the PFTs look normal?  Read fpc.out (foliar projective cover) at the Siberian
     boreal box and see what is growing.  The east-Siberian forest is dominated by LARCH,
     which in LPJ-GUESS is BNS (boreal needleleaf summergreen).  BNE/BINE are the
     evergreen conifers (spruce, pine), IBS the boreal broadleaf birch, C3G grass.  A
     boreal zone that comes out as grass, or as bare ground, is the failure mode this
     campaign exists to prevent.
  2. Is summer warm enough to let BNS grow where it should?  LPJ-GUESS gates
     establishment on bioclimatic limits, and the binding one in the boreal is GROWING
     DEGREE DAYS above 5 C.  GDD5 is computed here from the model's own daily 2 m
     temperature, so it answers "does the atmosphere supply the growing season the
     vegetation scheme requires", which a JJA mean cannot.

WHAT IS COMPARED.  11E against 11D and 110Baseline, all at their last available decade.
Note 110Baseline stops a decade earlier, so its vegetation has had less time -- flagged,
not hidden.  Vegetation is also SLOW: these runs start from an LPJ-GUESS spin-up state and
50 coupled years is not long for forest composition to re-equilibrate, so a PFT difference
between arms is a weaker signal than a PFT ABSENCE in all of them.

READING IT.  The diagnostic that matters is not whether 11E beats 11D by a few percent of
cover.  It is whether the boreal tree PFTs are present at all, at plausible cover, in the
band where observations put closed boreal forest.
"""
import os, sys, glob
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')

R092 = '/work/bb1469/a270092/runtime/awiesm3-v3.4'
R270 = '/work/bb1469/a270270/runtime/awiesm3-v3.4'
RSPIN = '/work/bb1469/a270270/runtime/lpjg-spinup'
ARMS = [('11E +sw15+K1', 'Tuning_test_11E_swemin15_K1'),
        ('11D +fitted sw30',
         'Tuning_test_11D_G4_fitted_snow_depletion_useIFSsoiltemp_CRUNCEP_plus_CERES_init_newSeaIce'),
        ('110 baseline',
         'Tuning_test_110Baseline_09C_useIFSsoiltemp_CRUNCEPinit_newSeaIce'),
        # THE REFERENCE.  The offline 2000-year LPJ-GUESS spin-up forced with observed
        # CRUNCEP+CERES -- i.e. the vegetation this scheme produces when the ATMOSPHERE IS
        # RIGHT.  It is also the state 11D/11E initialise from (ini_parent_date 3900-12-31),
        # so it is both the target and the starting point.  If it has forest and the
        # coupled runs do not, the coupled atmosphere is destroying it.  If it does NOT
        # have forest, the problem is in LPJ-GUESS or its parameters and no amount of
        # atmospheric tuning will produce a boreal forest.
        ('offline CRUNCEP 2000y',
         'LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPandCERES_daily_variability')]
# Siberian boreal box, the campaign's own definition for the forest question
BOX = dict(la=(55.0, 70.0), lo=(60.0, 140.0))
TREES = ['BNE', 'BINE', 'BNS', 'TeNE', 'TeBS', 'IBS', 'TeBE']
GRASS = ['C3G', 'C4G']

print(__doc__)
print('=' * 100)


def find(run, fn):
    for B in (R092, R270, RSPIN):
        hits = sorted(glob.glob(f'{B}/{run}/outdata/lpj_guess/*/run*/{fn}'))
        if hits:
            return hits
    return []


def _read_tail(f, nlines=60000):
    """Header + last nlines of a huge .out file.

    The 2000-year offline spin-up writes ONE 13 GB, 20.6-million-line fpc.out per rank
    (~10 300 cells x 2000 years).  Reading it whole killed the process with no output.
    LPJ-GUESS writes chronologically, so the final year is at the end: take the header
    plus a tail long enough to contain several years, then filter to the last one.
    """
    import subprocess, io
    head = subprocess.run(['head', '-1', f], capture_output=True, text=True).stdout
    tail = subprocess.run(['tail', '-n', str(nlines), f], capture_output=True, text=True).stdout
    return pd.read_csv(io.StringIO(head + tail), sep=r'\s+')


def last_year(run, fn):
    frames = []
    for f in find(run, fn):
        try:
            big = os.path.getsize(f) > 200 * 1024**2
            d = _read_tail(f) if big else pd.read_csv(f, sep=r'\s+')
        except Exception:
            continue
        if 'Lat' not in d or 'Year' not in d:
            continue
        d = d[(d.Lat >= BOX['la'][0]) & (d.Lat <= BOX['la'][1]) &
              (d.Lon >= BOX['lo'][0]) & (d.Lon <= BOX['lo'][1])]
        if len(d):
            frames.append(d[d.Year == d.Year.max()])
    return pd.concat(frames, ignore_index=True) if frames else None


print(f'Siberian boreal box {BOX["la"][0]:.0f}-{BOX["la"][1]:.0f}N, '
      f'{BOX["lo"][0]:.0f}-{BOX["lo"][1]:.0f}E; last year of the last decade on disk\n')
print('1. FOLIAR PROJECTIVE COVER BY PFT  (fpc.out, box mean)')
print('-' * 100)
cols = None
rows = {}
for lab, run in ARMS:
    d = last_year(run, 'fpc.out')
    if d is None:
        print(f'  {lab:18s} no fpc.out found')
        continue
    have = [c for c in TREES + GRASS if c in d.columns]
    if cols is None:
        cols = have
        print(f'  {"arm":18s} {"cells":>6s} {"year":>6s} ' +
              ' '.join(f'{c:>7s}' for c in cols) + f' {"TREE":>8s} {"GRASS":>7s}')
    tree = d[[c for c in TREES if c in d.columns]].sum(axis=1)
    gr = d[[c for c in GRASS if c in d.columns]].sum(axis=1)
    rows[lab] = (d, tree.mean(), gr.mean())
    print(f'  {lab:18s} {len(d):6d} {int(d.Year.max()):6d} ' +
          ' '.join(f'{d[c].mean():7.3f}' for c in cols) +
          f' {tree.mean():8.3f} {gr.mean():7.3f}')

print('\n  BNS = larch, the dominant east-Siberian boreal tree.  BNE/BINE = evergreen')
print('  conifer, IBS = boreal broadleaf (birch), C3G = grass.')

# ---------------------------------------------------------------- 2. is it forest?
print('\n2. IS IT FOREST?  fraction of Siberian cells by dominant cover')
print('-' * 100)
print(f'  {"arm":18s} {"tree>0.5":>9s} {"tree>0.2":>9s} {"grass>tree":>11s} '
      f'{"BNS>0.1":>8s} {"~bare<0.1":>10s}')
for lab, (d, tm, gm) in rows.items():
    tree = d[[c for c in TREES if c in d.columns]].sum(axis=1)
    gr = d[[c for c in GRASS if c in d.columns]].sum(axis=1)
    n = len(d)
    bns = d['BNS'] if 'BNS' in d.columns else pd.Series(np.zeros(n))
    print(f'  {lab:18s} {100*(tree>0.5).mean():8.1f}% {100*(tree>0.2).mean():8.1f}% '
          f'{100*(gr>tree).mean():10.1f}% {100*(bns>0.1).mean():7.1f}% '
          f'{100*((tree+gr)<0.1).mean():9.1f}%')

# ---------------------------------------------------------------- 3. growing season
print('\n3. DOES THE ATMOSPHERE SUPPLY THE GROWING SEASON?  GDD5 from daily 2m T')
print('-' * 100)
print('  LPJ-GUESS gates boreal establishment on growing degree days above 5 C.')
print('  Indicative GDD5 minima: BNS (larch) ~350, BNE ~600, IBS ~350.\n')
try:
    import xarray as xr
    for lab, run in ARMS:
        B = next((x for x in (R092, R270, RSPIN) if os.path.isdir(f'{x}/{run}')), R270)
        fs = sorted(glob.glob(f'{B}/{run}/outdata/oifs/atm_remapped_1d_2t_*.nc'))
        if not fs:
            fs = sorted(glob.glob(f'{B}/{run}/outdata/oifs/atmos_1h_sfc_2t_*.nc'))
        if not fs:
            print(f'  {lab:18s} no daily 2t on disk -- GDD5 not computable')
            continue
        with xr.open_dataset(fs[-1], decode_times=False) as ds:
            v = [k for k in ds.data_vars if k in ('2t', 'tas', 't2m')][0]
            a = ds[v].values
            lat, lon = ds['lat'].values, ds['lon'].values
        if a.ndim == 4:
            a = a[:, 0]
        t = a - 273.15
        gdd = np.clip(t - 5.0, 0, None).sum(axis=0)      # daily sum over the year
        sy = (lat >= BOX['la'][0]) & (lat <= BOX['la'][1])
        sx = (lon >= BOX['lo'][0]) & (lon <= BOX['lo'][1])
        g = gdd[np.ix_(sy, sx)]
        w = np.broadcast_to(np.cos(np.deg2rad(lat[sy]))[:, None], g.shape)
        print(f'  {lab:18s} GDD5 box mean {np.average(g, weights=w):7.0f} | '
              f'cells above 350: {100*(g > 350).mean():5.1f}% | above 600: '
              f'{100*(g > 600).mean():5.1f}%   ({os.path.basename(fs[-1])})')
except Exception as e:
    print(f'  GDD5 step failed: {e}')

print('\n  A box that clears 350 almost everywhere can support larch; one that does not')
print('  cannot, and no amount of tuning the SNOW will fix it -- that would be a')
print('  growing-season problem, not an albedo one.')
