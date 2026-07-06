#!/bin/bash
#SBATCH --job-name=lpjgLR
#SBATCH --account=bb1469
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/work/bb1469/a270092/eval/lpjg_txt_LR/lpjgLR_%j.log
set -e
W=/work/bb1469/a270092/runtime/lpjg-spinup-develop/LR_1000y_spinup/run_19000101-28991231/work
OUT=/work/bb1469/a270092/eval/lpjg_txt_LR; mkdir -p $OUT
for v in fpc cmass dens anpp annual_burned_area; do
  first=$(for f in $W/run*/output/$v.out; do [ $(wc -l < "$f" 2>/dev/null) -gt 1 ] && echo "$f" && break; done)
  [ -z "$first" ] && { echo "$v: none"; continue; }
  head -1 "$first" > $OUT/${v}_sib.txt
  # per box-cell last-year row, Siberian box (Lat55-70,Lon60-140), across all ranks
  awk '$2>=55&&$2<=70&&$1>=60&&$1<=140 && $3+0==$3 {k=$1"_"$2; if($3>=y[k]){y[k]=$3;r[k]=$0}} END{for(i in r)print r[i]}' \
      $W/run*/output/$v.out >> $OUT/${v}_sib.txt
  echo "$v: rows=$(($(wc -l < $OUT/${v}_sib.txt)-1))"
done
echo "=== summarize ==="
source ~/loadconda.sh; conda activate reval
python3 - <<'PYEOF'
import pandas as pd, warnings; warnings.filterwarnings('ignore')
O="/work/bb1469/a270092/eval/lpjg_txt_LR"
def load(v):
    try: return pd.read_csv(f"{O}/{v}_sib.txt", delim_whitespace=True)
    except Exception as e: print(v,"load err",e); return None
def m(d,c): return float(d[c].mean()) if (d is not None and c in d) else float('nan')
fpc=load("fpc"); print("\nLR Siberian boreal cells:", 0 if fpc is None else len(fpc), " yr~", int(fpc.Year.max()) if fpc is not None else "?","\n")
print("COVER(fpc): BNE=%.3f BINE=%.3f BNS=%.3f IBS=%.3f | C3G=%.3f"%(m(fpc,'BNE'),m(fpc,'BINE'),m(fpc,'BNS'),m(fpc,'IBS'),m(fpc,'C3G')))
print("  TREEFPC=%.3f GRASSFPC=%.3f FORESTFPC=%.3f FRACH=%.3f | AGDD5=%.0f"%(m(fpc,'TREEFPC'),m(fpc,'GRASSFPC'),m(fpc,'FORESTFPC'),m(fpc,'FRACH'),m(fpc,'AGDD5')))
d=load("dens"); print("DENSITY: BNE=%.4f BINE=%.4f BNS=%.4f (tot=%.4f)"%(m(d,'BNE'),m(d,'BINE'),m(d,'BNS'),(m(d,'BNE')+m(d,'BINE')+m(d,'BNS'))))
c=load("cmass"); print("BIOMASS: BNE=%.3f BINE=%.3f BNS=%.3f | C3G=%.3f"%(m(c,'BNE'),m(c,'BINE'),m(c,'BNS'),m(c,'C3G')))
a=load("anpp"); print("NPP:     BNE=%.4f BINE=%.4f BNS=%.4f | C3G=%.4f"%(m(a,'BNE'),m(a,'BINE'),m(a,'BNS'),m(a,'C3G')))
f=load("annual_burned_area")
if f is not None:
    cc=[x for x in f.columns if x not in('Lon','Lat','Year')]
    print("FIRE:", {x:round(m(f,x),4) for x in cc[:4]})
td=(m(d,'BNE')+m(d,'BINE')+m(d,'BNS')) if d is not None else float('nan')
tn=(m(a,'BNE')+m(a,'BINE')+m(a,'BNS')) if a is not None else float('nan')
print("\nVERDICT: tree density=%.4f, tree NPP=%.4f vs grass NPP=%.4f"%(td,tn,m(a,'C3G')))
PYEOF
