"""Is the boreal tree collapse driven by the RADIATION deficit or the TEMPERATURE drop?

The forcing-transfer test showed that swapping CRUNCEP -> AMIP forcing alone costs -54 %
of Siberian TREEFPC, but could not say which variable carried it: both radiation and
temperature feed NPP, and therefore both feed the greff_min growth-efficiency mortality
that the report identifies as the actual killer (establishment gates are cleared).

A colleague's spin-ups separate them. Holding CRUNCEP temperature fixed and swapping only
the radiation to CERES isolates the radiation channel:

  CRUNCEPcalibrated_v3               CRUNCEP T, CRUNCEP radiation  (+21 W/m2 vs CERES)
  CRUNCEPandCERES                    CRUNCEP T, CERES radiation    (truth)
  CRUNCEPandCERES_daily_variability  as above, with daily variability retained
  LR_2000y_PIforcing_v3 (AMIP)       AMIP T,     AMIP radiation    (-11 W/m2 vs CERES)

So the radiation ladder is CRUNCEP > CERES > AMIP, with temperature constant across the
first three. Reading:

  * AGDD5 must be ~unchanged between arms 1-3. It is the control: it confirms only the
    radiation was swapped. If AGDD5 moves, the arms differ in more than radiation and the
    attribution below is void.
  * TREEFPC change from arm 1 -> 2 is the pure radiation effect, for a radiation change of
    known sign and roughly known size.
  * Comparing that with the full CRUNCEP -> AMIP drop apportions the -54 %.

All arms are 2000-yr spin-ups read at equilibrium (year 3900) on common gridcells.
"""
import numpy as np

SRC = {'CRUNCEP (T_cru, R_cru)': '/tmp/y3900/CRUNCEPcalibrated_v3.txt',
       'CRUNCEP+CERES (T_cru, R_ceres)': '/tmp/y3900/CRUNCEPandCERES.txt',
       'CRUNCEP+CERES dayvar': '/tmp/y3900/CRUNCEPandCERES_daily_variability.txt',
       'AMIP (T_amip, R_amip)': '/tmp/y3900/AMIP.txt'}
KEEP = ['TREEFPC', 'GRASSFPC', 'AGDD5']
BOX = {'NH 45N+': (45, 90, -180, 180), 'Siberia': (55, 75, 60, 180),
       'E. Siberia': (55, 75, 90, 160)}


def load(f):
    hdr = open(f).readline().split()
    ix = {n: i for i, n in enumerate(hdr)}
    out = {}
    for L in open(f):
        p = L.split()
        if len(p) < len(hdr) or p[2] == 'Year':
            continue
        out[(round(float(p[0]), 4), round(float(p[1]), 4))] = [float(p[ix[k]]) for k in KEEP]
    return out


D = {k: load(v) for k, v in SRC.items()}
common = sorted(set.intersection(*[set(d) for d in D.values()]))
lon = np.array([k[0] for k in common])
lat = np.array([k[1] for k in common])
M = {k: np.array([d[c] for c in common]) for k, d in D.items()}
print(f"{len(common)} gridcells common to all {len(D)} arms, equilibrium year 3900\n")

for bn, (la0, la1, lo0, lo1) in BOX.items():
    m = (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1)
    w = np.cos(np.deg2rad(lat[m]))
    print(f"--- {bn}  ({m.sum()} cells)")
    print(f"    {'arm':32s} {'TREEFPC':>9s} {'GRASSFPC':>9s} {'AGDD5':>9s}")
    base = None
    for k in SRC:
        v = [np.average(M[k][m, i], weights=w) for i in range(len(KEEP))]
        if base is None:
            base = v
            print(f"    {k:32s} {v[0]:9.3f} {v[1]:9.3f} {v[2]:9.1f}   (reference)")
        else:
            d = [v[i] - base[i] for i in range(3)]
            pc = 100 * d[0] / base[0] if base[0] > 1e-9 else np.nan
            print(f"    {k:32s} {v[0]:9.3f} {v[1]:9.3f} {v[2]:9.1f}   "
                  f"dTREE {d[0]:+.3f} ({pc:+.1f}%)  dAGDD5 {d[2]:+.1f}")
    print()

# ---- the apportionment ------------------------------------------------------
print("APPORTIONMENT of the CRUNCEP -> AMIP tree loss")
print(f"  {'region':12s} {'radiation only':>16s} {'full swap':>12s} {'radiation share':>16s}")
for bn, (la0, la1, lo0, lo1) in BOX.items():
    m = (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1)
    w = np.cos(np.deg2rad(lat[m]))
    t = {k: np.average(M[k][m, 0], weights=w) for k in SRC}
    rad = t['CRUNCEP+CERES (T_cru, R_ceres)'] - t['CRUNCEP (T_cru, R_cru)']
    full = t['AMIP (T_amip, R_amip)'] - t['CRUNCEP (T_cru, R_cru)']
    share = 100 * rad / full if abs(full) > 1e-9 else np.nan
    print(f"  {bn:12s} {rad:+16.3f} {full:+12.3f} {share:15.0f}%")
print("\nNOTE: the radiation-only arm swaps CRUNCEP->CERES (a REDUCTION of ~21 W/m2, i.e."
      "\ntowards truth), while the full swap goes CRUNCEP->AMIP (~32 W/m2, past truth)."
      "\nThe share is therefore a lower bound on the radiation contribution, not an exact split.")
