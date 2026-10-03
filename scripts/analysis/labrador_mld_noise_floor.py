"""Is a 9-winter Labrador collapse outside what the 1200 s history does on its own?

March MLD2 in the convection box (56-62N 60-50W, area-weighted, from
labrador_ice_lineage.py var=MLD2). Null distribution: every 9-winter window inside the
continuous 1200 s old-build record of each mixing era. Arms: 2121-2129 of each 2120
branch (2120 is the OpenIFS cold start). Collapse winter: < 200 m.
"""
import glob, numpy as np, pandas as pd
df = pd.concat(pd.read_csv(f, names=['exp', 'year', 'month', 'region', 'v'])
               for f in glob.glob('data/clim/lab_ice/mld_*.csv'))
m = df[(df.region == 'interior') & (df.month == 3)].set_index(['exp', 'year']).v.sort_index()
def ser(e, a, b): return m.loc[e].loc[a:b]
LIN = [('PICAL', 1850, 1919), ('PICAL_momixoff', 1920, 1939), ('PICAL_ccnice', 1940, 2129)]
lin = pd.concat([ser(e, a, b) for e, a, b in LIN])
ERAS = [('KPP (1900-2059)', 1900, 2059), ('TKE (2060-2079)', 2060, 2079),
        ('TKE+IDEMIX, ccnice (2080-2129)', 2080, 2129)]
def stats(x):
    x = np.asarray(x); c = x < 200; run = best = 0
    for v in c: run = run + 1 if v else 0; best = max(best, run)
    return x.mean(), c.sum(), best
def windows(x, n=9):
    x = np.asarray(x); return np.array([x[i:i+n].mean() for i in range(len(x) - n + 1)])
print('NULL: 9-winter means of March MLD (m), 1200 s old build')
nulls = {}
for lab, a, b in ERAS:
    x = lin.loc[a:b]; w = windows(x); nulls[lab] = w
    mu, nc, best = stats(x)
    print(f'  {lab:32s} n={len(x):3d} yr  mean {mu:5.0f}  collapse winters {nc:2d} ({100*nc/len(x):.0f} %)  '
          f'longest run {best}  9-yr window min {w.min():4.0f}  5th pct {np.percentile(w,5):4.0f}  median {np.median(w):4.0f}')
cv = pd.concat([ser('PICAL_crunveg', 2101, 2119)])
w = windows(cv); mu, nc, best = stats(cv)
print(f'  {"TKE+IDEMIX, crunveg (2101-2119)":32s} n={len(cv):3d} yr  mean {mu:5.0f}  collapse winters {nc:2d}  longest run {best}  9-yr window min {w.min():4.0f}  median {np.median(w):4.0f}')
null = np.concatenate([nulls[ERAS[2][0]], w])
print(f'\nPooled IDEMIX-era null (ccnice 2080-2129 + crunveg 2101-2119): {len(null)} windows, '
      f'min {null.min():.0f}, 5th pct {np.percentile(null,5):.0f}, median {np.median(null):.0f}')
ARMS = [('PICAL_ccnice', 'ccnice cont., old 1200 s'), ('PICAL_crunveg_ob1200', 'old build 1200 s'),
        ('PICAL_crunveg_ts1200', 'new build 1200 s'), ('PICAL_crunveg_ob1800', 'old build 1800 s'),
        ('PICAL_crunveg_ihf0', 'new build 1800 s'), ('PICAL_crunveg_ihf1', 'new 1800 s + McPhee'),
        ('PICAL_crunveg_nx1800', 'new 1800 s + IDEMIX fix')]
print('\nARMS, 2121-2129 March MLD; percentile = share of null windows at or below the arm mean')
for e, lab in ARMS:
    x = ser(e, 2121, 2129); mu, nc, best = stats(x)
    print(f'  {lab:26s} {e:22s} n={len(x)}  mean {mu:5.0f}  pct {100*np.mean(null <= mu):5.1f}  collapses {nc}  longest {best}  '
          + ' '.join(f'{v:4.0f}' for v in x))
