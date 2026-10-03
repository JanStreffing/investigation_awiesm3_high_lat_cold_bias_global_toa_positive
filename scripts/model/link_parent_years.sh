#!/bin/bash
# Link a parent run's output years into a branch's outdata, so the branch folder reads as one
# continuous record from the parent's start.
#
# Written for PICAL -> PICAL_momixoff on 2026-09-18: momix-off branched at 1920 and becomes the
# production line, so PICAL's 1850-1919 belong in front of it.  Symlinks, never copies: the years
# stay owned by the parent, and a reader can always tell which run actually produced a year by
# whether the file is a link.
#
# Usage: link_parent_years.sh <parent> <branch> <first year> <last year> [--go]
#        without --go it only reports what it would do.
set -u
ROOT=/work/bb1469/a270092/runtime/awiesm3-v3.4
P=$1; B=$2; Y0=$3; Y1=$4; GO=${5:-}
n=0; skip=0
for sub in oifs fesom; do
    src=$ROOT/$P/outdata/$sub; dst=$ROOT/$B/outdata/$sub
    [ -d "$src" ] || continue
    mkdir -p "$dst"
    for f in "$src"/*; do
        b=$(basename "$f")
        # year is either  ..._YYYY-YYYY.nc  (oifs)  or  ....YYYY.nc  (fesom, incl .gr.)
        y=$(echo "$b" | grep --color=never -oE '[0-9]{4}-[0-9]{4}\.nc$' | cut -d- -f1)
        [ -z "$y" ] && y=$(echo "$b" | grep --color=never -oE '\.[0-9]{4}\.nc$' | tr -d './nc')
        [ -z "$y" ] && continue
        [ "$y" -ge "$Y0" ] && [ "$y" -le "$Y1" ] || continue
        if [ -e "$dst/$b" ]; then skip=$((skip+1)); continue; fi
        n=$((n+1))
        [ "$GO" = "--go" ] && ln -s "$f" "$dst/$b"
    done
done
# LPJ-GUESS output is per-leg directories, not per-year files
for d in $ROOT/$P/outdata/lpj_guess/*/; do
    b=$(basename "$d"); y=${b:0:4}
    [ "$y" -ge "$Y0" ] && [ "$y" -le "$Y1" ] || continue
    if [ -e "$ROOT/$B/outdata/lpj_guess/$b" ]; then skip=$((skip+1)); continue; fi
    n=$((n+1))
    # strip the trailing slash the */ glob leaves: "ln -s dir/ target" does not make a
    # symlink named target, it tries to place the link INSIDE it.  This silently created
    # nothing on 2026-09-18 and the LPJ-GUESS years had to be linked by hand.
    [ "$GO" = "--go" ] && ln -s "${d%/}" "$ROOT/$B/outdata/lpj_guess/$b"
done
echo "$P $Y0-$Y1 -> $B : $n link(s) $([ "$GO" = "--go" ] && echo created || echo 'would be created'), $skip already present"
