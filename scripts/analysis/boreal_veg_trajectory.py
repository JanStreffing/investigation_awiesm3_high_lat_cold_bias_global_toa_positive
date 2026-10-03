"""Year-by-year Siberian vegetation trajectory from the LPJ-GUESS annual tables.

WHY THIS EXISTS.  siberian_veg_matched.py compares one decade's final year across arms.
The questions after the 2100 LPJ-GUESS merge (canopy work: ltos, shrubs to low IFS
vegetation, iftreefracca 1) are about TRAJECTORIES: does the forest come back on the
production line after 2100, and does a forest started from the CRUNCEP state survive
longer than it did in 080a/11E before the merge.  So every year of every leg is read.

THE FRAME.  BNS (larch) is not gate-limited in this box; it is outcompeted by C3 grass
(skill bns-not-gate-limited).  So besides cover, the table carries the competition
terms: LAI and NPP of larch against grass, larch density (establishment success), and
grassPAR, the light that reaches the grass layer.  FRACH is what OpenIFS actually
receives as high-vegetation fraction, which is the coupling-relevant number.

OUTPUT.  One CSV per experiment: year, n cells, cos(lat)-weighted box means, and the
share of cells where grass FPC exceeds tree FPC.

Usage:  boreal_veg_trajectory.py <tag> <exp root> <out csv> [box]
        box = sib (55-75N, 60-180E, default) | boreal (50-70N, all longitudes)
"""
import glob
import sys
import numpy as np
import pandas as pd

BOXES = {'sib': (55.0, 75.0, 60.0, 180.0), 'boreal': (50.0, 70.0, 0.0, 360.0)}

WANT = {
    'fpc': ['BNE', 'BINE', 'BNS', 'IBS', 'C3G', 'HSE', 'HSS', 'LSE', 'LSS',
            'TREEFPC', 'GRASSFPC', 'FRACH', 'FRACL', 'AGDD5', 'PEATFPC', 'Total'],
    'lai': ['BNE', 'BNS', 'IBS', 'C3G', 'LSE', 'LSS', 'Total'],
    'anpp': ['BNE', 'BNS', 'IBS', 'C3G', 'Total'],
    'dens': ['BNE', 'BNS', 'IBS'],
    'cmass': ['BNS', 'IBS', 'C3G', 'Total'],
    'est_limits': ['mTmin20', 'agdd5', 'snowdepth', 'grassPAR', 'wcont_up'],
}


def read_table(fn, cols, box):
    hdr = open(fn).readline().split()
    use = ['Lon', 'Lat', 'Year'] + [c for c in cols if c in hdr]
    df = pd.read_csv(fn, sep=r'\s+', usecols=use)
    lo = df['Lon'] % 360
    la0, la1, lo0, lo1 = box
    return df[(df['Lat'] >= la0) & (df['Lat'] <= la1) & (lo >= lo0) & (lo <= lo1)]


def trajectory(root, box):
    legs = sorted(glob.glob(f'{root}/outdata/lpj_guess/*/run*/fpc.out'))
    frames = []
    for fpcfile in legs:
        d = fpcfile.rsplit('/', 1)[0]
        merged = None
        for tab, cols in WANT.items():
            fn = f'{d}/{tab}.out'
            try:
                t = read_table(fn, cols, box)
            except (OSError, ValueError):
                continue
            t = t.rename(columns={c: f'{tab}_{c}' for c in cols})
            merged = t if merged is None else merged.merge(t, on=['Lon', 'Lat', 'Year'], how='left')
        if merged is not None and len(merged):
            frames.append(merged)
    if not frames:
        return None
    df = pd.concat(frames)
    # overlapping legs (e.g. a 5-year and a 10-year folder for the same years): keep one
    df = df.drop_duplicates(subset=['Lon', 'Lat', 'Year'], keep='last')
    df['w'] = np.cos(np.deg2rad(df['Lat']))
    df['grass_gt_tree'] = (df['fpc_GRASSFPC'] > df['fpc_TREEFPC']).astype(float)
    vals = [c for c in df.columns if c not in ('Lon', 'Lat', 'Year', 'w')]
    out = []
    for yr, g in df.groupby('Year'):
        row = {'year': yr, 'ncell': len(g)}
        for c in vals:
            x = g[c]
            ok = x.notna()
            row[c] = np.average(x[ok], weights=g['w'][ok]) if ok.any() else np.nan
        out.append(row)
    return pd.DataFrame(out)


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(2)
    tag, root, out = sys.argv[1:4]
    box = BOXES[sys.argv[4] if len(sys.argv) > 4 else 'sib']
    t = trajectory(root, box)
    if t is None:
        print(f'{tag}: no fpc.out under {root}')
        sys.exit(1)
    t.insert(0, 'exp', tag)
    t.to_csv(out, index=False, float_format='%.5g')
    show = ['year', 'ncell', 'fpc_BNS', 'fpc_TREEFPC', 'fpc_GRASSFPC', 'fpc_FRACH',
            'grass_gt_tree', 'lai_BNS', 'lai_C3G', 'dens_BNS', 'est_limits_grassPAR']
    print(t[[c for c in show if c in t]].to_string(index=False, float_format='%.4f'))


if __name__ == '__main__':
    main()
