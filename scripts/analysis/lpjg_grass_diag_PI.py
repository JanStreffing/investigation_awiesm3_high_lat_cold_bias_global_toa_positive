#!/usr/bin/env python3
"""
Diagnose grass-vs-tree competition in the LPJG HR 2000yr spin-up, directly from
the raw .out text files (no netCDF). Reads last-year equilibrium at Siberian
boreal cells across all MPI ranks and characterizes the mechanism:
  - dens (tree density): ~0 => establishment failure; >0 => trees establish
  - cmass (biomass), anpp (NPP): tree vs grass productivity/standing stock
  - fpc + aggregates (TREEFPC/GRASSFPC/FRACH): who wins the cover
  - annual_burned_area: is fire resetting to grass?
Low memory: each rank file is read, filtered to the Siberian box + last year, discarded.
"""
import glob, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
import sys
R = sys.argv[1] if len(sys.argv) > 1 else "/work/bb1469/a270092/runtime/lpjg-spinup-develop/LR_2000y_PIforcing_1job/outdata/lpj_guess"
BOX = dict(la=(55, 70), lo=(60, 140))   # Siberian boreal

def gather(fn):
    frames = []
    for f in glob.glob(f"{R}/run*/{fn}"):
        try:
            d = pd.read_csv(f, delim_whitespace=True)
            d = d[(d.Lat >= BOX['la'][0]) & (d.Lat <= BOX['la'][1]) &
                  (d.Lon >= BOX['lo'][0]) & (d.Lon <= BOX['lo'][1])]
            if len(d):
                frames.append(d[d.Year == d.Year.max()])
        except Exception:
            pass
    return pd.concat(frames, ignore_index=True) if frames else None

def mean(df, c): return float(df[c].mean()) if (df is not None and c in df) else float('nan')

fpc = gather("fpc.out")
if fpc is None:
    print("NO Siberian cells found"); raise SystemExit
print(f"Siberian boreal cells (55-70N, 60-140E): {len(fpc)}  | equilibrium year {int(fpc.Year.max())}\n")

print("=== COVER: FPC fractional cover (who wins) ===")
print(f"  trees:  BNE={mean(fpc,'BNE'):.3f}  BINE={mean(fpc,'BINE'):.3f}  BNS(larch)={mean(fpc,'BNS'):.3f}  IBS={mean(fpc,'IBS'):.3f}")
print(f"  grass:  C3G={mean(fpc,'C3G'):.3f}")
print(f"  aggr:   TREEFPC={mean(fpc,'TREEFPC'):.3f}  GRASSFPC={mean(fpc,'GRASSFPC'):.3f}  FORESTFPC={mean(fpc,'FORESTFPC'):.3f}  FRACH(cvh)={mean(fpc,'FRACH'):.3f}")
print(f"  climate: AGDD5={mean(fpc,'AGDD5'):.0f}  AWCONT={mean(fpc,'AWCONT'):.3f}")

dens = gather("dens.out")
print("\n=== ESTABLISHMENT vs MORTALITY: tree density [indiv/m2] ===")
print(f"  BNE={mean(dens,'BNE'):.4f}  BINE={mean(dens,'BINE'):.4f}  BNS={mean(dens,'BNS'):.4f}   (near-0 => trees never establish; >0 => they establish but may die)")

cm = gather("cmass.out")
print("\n=== STANDING BIOMASS: cmass [kgC/m2] ===")
print(f"  trees:  BNE={mean(cm,'BNE'):.3f}  BINE={mean(cm,'BINE'):.3f}  BNS={mean(cm,'BNS'):.3f}   | grass C3G={mean(cm,'C3G'):.3f}")

an = gather("anpp.out")
print("\n=== PRODUCTIVITY: anpp [kgC/m2/yr]  (low tree NPP => greff-mortality risk) ===")
print(f"  trees:  BNE={mean(an,'BNE'):.4f}  BINE={mean(an,'BINE'):.4f}  BNS={mean(an,'BNS'):.4f}   | grass C3G={mean(an,'C3G'):.4f}")

fire = gather("annual_burned_area.out")
if fire is not None:
    cc = [c for c in fire.columns if c not in ('Lon', 'Lat', 'Year')]
    print("\n=== FIRE: annual_burned_area ===")
    for c in cc[:6]:
        print(f"  {c} = {mean(fire,c):.4f}")

print("\n--- verdict hints ---")
if dens is not None:
    td = mean(dens,'BNE')+mean(dens,'BINE')+mean(dens,'BNS')
    print(f"  total tree density = {td:.4f} indiv/m2  ({'NEAR-ZERO => establishment failure' if td<0.001 else 'nonzero => trees do establish -> mortality/competition'})")
