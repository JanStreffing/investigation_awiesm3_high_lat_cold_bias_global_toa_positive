import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
O="/work/bb1469/a270092/eval/lpjg_txt"
def load(v):
    try: return pd.read_csv(f"{O}/{v}_sib.txt", delim_whitespace=True)
    except: return None
def m(d,c): return float(d[c].mean()) if (d is not None and c in d) else float('nan')
fpc=load("fpc"); print(f"\nSiberian boreal cells: {0 if fpc is None else len(fpc)}\n")
print("COVER (fpc):  BNE=%.3f BINE=%.3f BNS=%.3f IBS=%.3f | C3G=%.3f"%(m(fpc,'BNE'),m(fpc,'BINE'),m(fpc,'BNS'),m(fpc,'IBS'),m(fpc,'C3G')))
print("  aggr: TREEFPC=%.3f GRASSFPC=%.3f FORESTFPC=%.3f FRACH=%.3f | AGDD5=%.0f AWCONT=%.3f"%(m(fpc,'TREEFPC'),m(fpc,'GRASSFPC'),m(fpc,'FORESTFPC'),m(fpc,'FRACH'),m(fpc,'AGDD5'),m(fpc,'AWCONT')))
d=load("dens"); print("\nDENSITY dens[indiv/m2]: BNE=%.4f BINE=%.4f BNS=%.4f  (total=%.4f)"%(m(d,'BNE'),m(d,'BINE'),m(d,'BNS'),m(d,'BNE')+m(d,'BINE')+m(d,'BNS')))
c=load("cmass"); print("BIOMASS cmass[kgC/m2]: BNE=%.3f BINE=%.3f BNS=%.3f | C3G=%.3f"%(m(c,'BNE'),m(c,'BINE'),m(c,'BNS'),m(c,'C3G')))
a=load("anpp"); print("NPP anpp[kgC/m2/yr]:   BNE=%.4f BINE=%.4f BNS=%.4f | C3G=%.4f"%(m(a,'BNE'),m(a,'BINE'),m(a,'BNS'),m(a,'C3G')))
f=load("annual_burned_area")
if f is not None:
    cc=[x for x in f.columns if x not in('Lon','Lat','Year')]
    print("FIRE annual_burned_area:", {x:round(m(f,x),4) for x in cc[:5]})
td=(m(d,'BNE')+m(d,'BINE')+m(d,'BNS')) if d is not None else float('nan')
print("\nVERDICT: total tree density=%.4f -> %s"%(td,"ESTABLISHMENT FAILURE (near-0)" if td<0.001 else "trees establish -> mortality/competition/fire"))
