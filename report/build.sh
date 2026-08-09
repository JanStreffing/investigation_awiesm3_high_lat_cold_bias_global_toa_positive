#!/bin/bash
# Build the campaign report and/or summary, then compress the PDF.
#
# WHY THE COMPRESSION STEP EXISTS.  The figures are matplotlib output at whatever dpi
# the plotting script happened to use, and several of the coupled NH maps are
# 4252x3767 px.  Printed at one text width (17.4 cm) that is ~620 dpi -- about three
# times more linear resolution than a PDF viewer or a printer can use.  Nothing is
# gained by carrying it, and report.pdf had reached 55 MB, over GitHub's 50 MB warning
# threshold and roughly ten times the size of everything else in the repository
# combined.
#
# Downsampling to 200 dpi takes report.pdf from 55 MB to ~5 MB with no visible loss:
# the raster figures are resampled bicubically, and all TEXT stays vector because
# pdfwrite re-embeds the fonts rather than rasterising them.  200 was chosen over
# ghostscript's /ebook preset (150 dpi) because the dense categorical maps -- the
# dominant-PFT panels especially -- show speckle at 150.
#
# THE SOURCE PNGs ARE LEFT ALONE deliberately.  Several are produced by a colleague's
# plotting scripts and get regenerated; downsampling them in place would be undone on
# the next run and would also degrade them for any other use.  Compressing at the PDF
# stage keeps the inputs authoritative and makes this idempotent.
#
# Usage:  ./build.sh            both documents
#         ./build.sh report     one of them
#         ./build.sh summary
set -eu
cd "$(dirname "$0")"

export PATH=/sw/spack-levante/texlive-live2025-3r2myy/bin/x86_64-linux:$PATH
GS=/sw/spack-levante/ghostscript-9.54.0-ed7q6u/bin/gs
DPI=200

DOCS=${*:-"report summary"}

for doc in $DOCS; do
  [ -f "$doc.tex" ] || { echo "no $doc.tex"; exit 1; }
  echo "=== $doc"

  # twice, so \ref and the table of contents resolve
  for pass in 1 2; do
    pdflatex -interaction=nonstopmode "$doc.tex" >/dev/null 2>&1 || true
  done

  err=$(grep -c '^!' "$doc.log" || true)
  und=$(grep -ci 'undefined' "$doc.log" || true)
  pages=$(grep -o 'Output written on.*' "$doc.log" | grep -o '[0-9]* pages' || echo '? pages')
  echo "    latex: $err errors, $und undefined refs, $pages"
  if [ "$err" != "0" ]; then
    echo "    *** LaTeX errors -- not compressing, fix them first"
    grep -n '^!' "$doc.log" | head -5
    continue
  fi

  before=$(stat -c%s "$doc.pdf")
  $GS -q -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite \
      -dCompatibilityLevel=1.5 -dDetectDuplicateImages=true \
      -dDownsampleColorImages=true -dColorImageDownsampleType=/Bicubic -dColorImageResolution=$DPI \
      -dDownsampleGrayImages=true  -dGrayImageDownsampleType=/Bicubic  -dGrayImageResolution=$DPI \
      -dDownsampleMonoImages=true  -dMonoImageDownsampleType=/Subsample -dMonoImageResolution=300 \
      -dAutoFilterColorImages=false -dColorImageFilter=/DCTEncode \
      -dEmbedAllFonts=true -dSubsetFonts=true \
      -sOutputFile="$doc.compressed.pdf" "$doc.pdf"

  # only accept the compressed file if it is smaller AND has the same page count --
  # a truncated ghostscript run must never silently replace a good PDF
  pb=$(pdfinfo "$doc.pdf" 2>/dev/null | awk '/^Pages/{print $2}' || echo '')
  pa=$(pdfinfo "$doc.compressed.pdf" 2>/dev/null | awk '/^Pages/{print $2}' || echo '')
  after=$(stat -c%s "$doc.compressed.pdf")
  if [ -n "$pb" ] && [ -n "$pa" ] && [ "$pb" != "$pa" ]; then
    echo "    *** page count changed $pb -> $pa, keeping the uncompressed PDF"
    rm -f "$doc.compressed.pdf"
  elif [ "$after" -ge "$before" ]; then
    echo "    already small ($((before/1048576)) MB); keeping as is"
    rm -f "$doc.compressed.pdf"
  else
    mv "$doc.compressed.pdf" "$doc.pdf"
    printf '    compressed: %.1f MB -> %.1f MB (%.1fx) at %s dpi\n' \
      "$(echo "$before/1048576" | bc -l)" "$(echo "$after/1048576" | bc -l)" \
      "$(echo "$before/$after" | bc -l)" "$DPI"
  fi

  rm -f "$doc.aux" "$doc.out" "$doc.toc" "$doc.log"
done
