"""Report figure: where the model's heat uptake actually lands, and how it gets there.

Panels (a) and (b) locate the heat in depth and in latitude; (c) and (d) split the
fixed-depth warming into the part that crossed density surfaces and the part that is
just the surfaces moving.  Inputs are produced by
scripts/analysis/ohc_drift_by_depth_band.py and scripts/analysis/heave_vs_watermass.py.
"""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, 'report', 'run_figures', 'ocean_heat_routing_11X.png')
ARM = '11X'

d = np.load(os.path.join(REPO, 'data', f'ohc_drift_{ARM}.npz'))
dT, ztop, zbot, lat, area, dzlay = d['dT'], d['ztop'], d['zbot'], d['lat'], d['area'], d['dz']
h = np.load(os.path.join(REPO, 'data', f'heave_watermass_{ARM}.npz'))
SIG, dth, hv, zbar = h['SIG'], h['dth'], h['dz'], h['zbar']

wet = np.isfinite(dT)
vol = area[:, None] * dzlay[None, :] * wet
zmid = 0.5 * (ztop + zbot)

BANDS = [('90-60S', -90, -60, '#1f3b73'), ('60-45S', -60, -45, '#4a90c0'),
         ('45-30S', -45, -30, '#7ec8a9'), ('30S-30N', -30, 30, '#e08214'),
         ('30-45N', 30, 45, '#c0392b'), ('45-60N', 45, 60, '#8e44ad'),
         ('60-90N', 60, 90, '#5d6d7e')]
SLABS = [('0-100 m', 0, 100), ('100-300 m', 100, 300), ('300-700 m', 300, 700),
         ('700-2000 m', 700, 2000), ('>2000 m', 2000, 1e9)]
SLABC = ['#f6d5a8', '#e8a55c', '#c0392b', '#7b2d26', '#3d1714']
RHO_CP = 1026.0 * 3996.0
COVER = 0.25          # a sigma level must be present in >=25% of the band's nodes
YRS = float(d['yrs'])
SEC = YRS * 365.0 * 86400.0
gocean = area[np.any(wet, axis=1)].sum()

fig, ax = plt.subplots(2, 2, figsize=(15.5, 11))

# ---- (a) temperature drift profile -------------------------------------------------
a = ax[0, 0]
for name, lo, hi, c in BANDS:
    sel = (lat >= lo) & (lat < hi)
    prof = []
    for k in range(len(ztop)):
        s = sel & wet[:, k]
        w = vol[s, k]
        prof.append(np.sum(dT[s, k] * w) / max(w.sum(), 1e-30) if w.sum() > 0 else np.nan)
    a.plot(prof, zmid, color=c, lw=1.6, label=name)
gprof = [np.sum(dT[wet[:, k], k] * vol[wet[:, k], k]) / max(vol[wet[:, k], k].sum(), 1e-30)
         for k in range(len(ztop))]
a.plot(gprof, zmid, 'k-', lw=3.0, label='GLOBAL', zorder=5)
a.axvline(0, color='0.5', lw=0.8)
a.set_ylim(2500, 0)
a.set_xlabel(f'temperature drift over {YRS:.0f} yr [K]')
a.set_ylabel('depth [m]')
a.set_title('(a) the heat is not in the top 100 m', fontsize=11, loc='left')
a.grid(alpha=0.3)
a.legend(fontsize=8, ncol=2)

# ---- (b) heat content by band and slab ---------------------------------------------
b = ax[0, 1]
ypos = np.arange(len(BANDS))
left_pos = np.zeros(len(BANDS))
left_neg = np.zeros(len(BANDS))
for j, (sn, z0, z1) in enumerate(SLABS):
    lev = (ztop >= z0) & (ztop < z1)
    vals = []
    for name, lo, hi, _ in BANDS:
        s = (lat >= lo) & (lat < hi)
        q = np.nansum(dT[s][:, lev] * vol[s][:, lev]) * RHO_CP
        vals.append(q / gocean / SEC)
    vals = np.array(vals)
    base = np.where(vals >= 0, left_pos, left_neg)
    b.barh(ypos, vals, left=base, color=SLABC[j], label=sn, height=0.65,
           edgecolor='white', linewidth=0.5)
    left_pos = left_pos + np.maximum(vals, 0)
    left_neg = left_neg + np.minimum(vals, 0)
b.axvline(0, color='0.3', lw=1.0)
b.set_yticks(ypos)
b.set_yticklabels([n for n, _, _, _ in BANDS], fontsize=9)
b.invert_yaxis()
b.set_xlabel('share of the global ocean heat gain [W m$^{-2}$ of global ocean]')
b.set_title(f'(b) total +{sum(left_pos)+sum(left_neg):.2f} W m$^{{-2}}$; '
            f'the Southern Ocean LOSES heat', fontsize=11, loc='left')
b.grid(alpha=0.3, axis='x')
b.legend(fontsize=8, loc='lower right')

# ---- (c) water-mass change, (d) heave ----------------------------------------------
SHOW = [('45-30S', -45, -30, '#7ec8a9'), ('30S-30N', -30, 30, '#e08214'),
        ('30-45N', 30, 45, '#c0392b'), ('45-60N', 45, 60, '#8e44ad')]
for panel, field, lab, ttl in ((ax[1, 0], dth, r'$\Delta\theta$ on $\sigma_0$ [K]',
                                '(c) heat that CROSSED density surfaces'),
                               (ax[1, 1], hv, r'$\Delta z$ of the $\sigma_0$ surface [m]',
                                '(d) heave: the surfaces moving')):
    # A sigma level is only plotted where it actually exists over a decent part of the
    # band.  Without this, levels that outcrop in a handful of nodes produce excursions
    # of hundreds of metres of "heave" that are sampling, not ocean.
    for name, lo, hi, c in SHOW:
        sel = (lat >= lo) & (lat < hi)
        nband = float(sel.sum())
        xs, zs = [], []
        for j in range(SIG.size):
            g = sel & np.isfinite(field[:, j]) & np.isfinite(zbar[:, j])
            if g.sum() / nband < COVER:
                continue
            w = area[g]
            zz = np.sum(zbar[g, j] * w) / w.sum()
            if zz < 80 or zz > 1600:
                continue
            zs.append(zz)
            xs.append(np.sum(field[g, j] * w) / w.sum())
        if len(zs) < 3:
            continue
        zs = np.array(zs)
        xs = np.array(xs)
        keep = np.concatenate([[True], np.diff(zs) > 0])   # sigma levels must deepen
        panel.plot(xs[keep], zs[keep], color=c, lw=1.8, marker='o', ms=3, label=name)
    panel.axvline(0, color='0.5', lw=0.8)
    panel.set_ylim(1600, 0)
    panel.set_xlabel(lab)
    panel.set_ylabel('mean depth of the isopycnal [m]')
    panel.set_title(ttl, fontsize=11, loc='left')
    panel.grid(alpha=0.3)
    panel.legend(fontsize=8)

fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=140, bbox_inches='tight')
print('wrote', OUT)
