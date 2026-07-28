"""Forcing-transfer test: same LPJ-GUESS, two forcings, 2000 yr to equilibrium.

How much of the boreal tree deficit is caused simply by moving LPJ-GUESS between the
forcing it was spun up under and the forcing our atmosphere delivers? Both arms already
existed:

  CRUNCEP arm  /work/bb1469/a270270/runtime/lpjg-spinup/
                 LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPcalibrated_v3
  AMIP arm     /work/bb1469/a270092/runtime/lpjg-spinup-develop/LR_2000y_PIforcing_v3

Controlled: byte-identical global.ins, firemodel "GLOBFIRM" in both .ins files, same
nyear_spinup=2000 / freenyears=100, same fixed 1850 CO2/ndep/LU, both reaching year 3900,
10043 gridcells in common. The material difference is the forcing file.

CAVEATS (both recorded in the report):
  * different binaries (a270092 vs a270270 build of the same LPJG version)
  * the AMIP forcing is the SUPERSEDED build, ~0.9 K warmer over Siberian JJA land than
    the corrected one, so the measured penalty is a LOWER BOUND
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = '/work/ab0246/a270092/postprocessing/investigation_awiesm3_high_lat_cold_bias_global_toa_positive/'
AMIP = ('/work/bb1469/a270092/runtime/lpjg-spinup-develop/LR_2000y_PIforcing_v3/'
        'outdata/lpj_guess/19000101-38991231/run1/fpc.out')
CRU = ('/work/bb1469/a270270/runtime/lpjg-spinup/LPJG-SPINUP_2000Y_TCO95_CORE3_CRUNCEPcalibrated_v3/'
       'outdata/lpj_guess/19000101-38991231/run1/fpc.out')
KEEP = ['TREEFPC', 'GRASSFPC', 'FORESTFPC', 'AGDD5']
YEAR = 3900


def load(f):
    hdr = open(f).readline().split()
    ix = {n: i for i, n in enumerate(hdr)}
    out = {}
    for L in open(f):
        p = L.split()
        if len(p) < len(hdr) or p[2] == 'Year':
            continue
        if int(p[2]) != YEAR:
            continue
        out[(round(float(p[0]), 4), round(float(p[1]), 4))] = [float(p[ix[k]]) for k in KEEP]
    return out


a, c = load(AMIP), load(CRU)
common = sorted(set(a) & set(c))
lat = np.array([k[1] for k in common])
A = np.array([a[k] for k in common])
C = np.array([c[k] for k in common])

BOX = {'NH 45N+': (45, 90, -180, 180),
       'Siberia': (55, 75, 60, 180),
       'E. Siberia': (55, 75, 90, 160)}
lon = np.array([k[0] for k in common])
print(f"{len(common)} common cells at equilibrium (yr {YEAR})\n")
print(f"  {'region':12s} {'var':10s} {'CRUNCEP':>9s} {'AMIP':>9s} {'diff':>9s} {'%':>7s}")
for bn, (la0, la1, lo0, lo1) in BOX.items():
    m = (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1)
    w = np.cos(np.deg2rad(lat[m]))
    for k in KEEP:
        i = KEEP.index(k)
        va, vc = np.average(A[m, i], weights=w), np.average(C[m, i], weights=w)
        pct = 100 * (va - vc) / vc if vc > 1e-9 else np.nan
        print(f"  {bn:12s} {k:10s} {vc:9.3f} {va:9.3f} {va-vc:+9.3f} {pct:+7.1f}")

# ---- zonal figure -----------------------------------------------------------
SURF, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8983'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.linewidth': 0.6,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.labelcolor': INK2, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF})
BLUE, ORANGE = '#2a78d6', '#eb6834'
edges = np.arange(40, 82, 2.0)
mid = 0.5 * (edges[:-1] + edges[1:])
fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.9))
for ax, (k, ttl) in zip(axes, [('TREEFPC', 'tree cover fraction'),
                               ('GRASSFPC', 'grass cover fraction'),
                               ('AGDD5', 'growing degree days above 5°C')]):
    i = KEEP.index(k)
    za, zc = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (lat >= lo) & (lat < hi)
        if m.sum() == 0:
            za.append(np.nan); zc.append(np.nan); continue
        w = np.cos(np.deg2rad(lat[m]))
        za.append(np.average(A[m, i], weights=w)); zc.append(np.average(C[m, i], weights=w))
    ax.plot(mid, zc, color=ORANGE, lw=2.2, marker='o', ms=3.2, mec=SURF, mew=0.8,
            label='CRUNCEP-forced (the spin-up in use)')
    ax.plot(mid, za, color=BLUE, lw=2.2, marker='o', ms=3.2, mec=SURF, mew=0.8,
            label='AMIP-forced (our atmosphere)')
    ax.axvspan(55, 75, color=MUTED, alpha=0.10, lw=0)
    ax.set_title(ttl, fontsize=9.5, color=INK, loc='left', fontweight='bold', pad=6)
    ax.set_xlabel('latitude'); ax.grid(axis='y', color=MUTED, alpha=0.22, lw=0.5)
    ax.set_axisbelow(True); ax.set_xlim(40, 80)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
axes[0].set_ylabel('fraction'); axes[2].set_ylabel('degree-days')
axes[0].legend(frameon=False, fontsize=8.2, labelcolor=INK2, loc='upper right')
fig.suptitle('Forcing transfer alone halves boreal tree cover — identical LPJ-GUESS, two forcings',
             x=0.008, ha='left', fontsize=11.5, fontweight='bold', color=INK, y=1.015)
fig.text(0.008, 0.925, 'Zonal means at equilibrium (yr 3900), 10043 common gridcells, byte-identical '
         'global.ins. Shaded band = the Siberian box.', ha='left', fontsize=8.2, color=INK2)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(REPO + 'plots/forcing_transfer_test.png', dpi=170, bbox_inches='tight')
print('\nwrote plots/forcing_transfer_test.png')
