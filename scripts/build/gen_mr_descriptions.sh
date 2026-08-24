#!/bin/bash
# Regenerate merge_request_descriptions.md from the pushed pr* branches.
#
# Lives in the repo rather than the scratchpad because the session scratchpad
# has been wiped mid-run twice, once truncating the generated file itself.
# Everything it needs is read from origin/pr* and upstream/main, so it can be
# re-run from scratch after a wipe.
#
# The stack sentence at the top of each description is written here and NOT in
# the commit message, because the merge request numbers are a prediction and
# should not be baked into git history.
set -u
REPO=${1:-/work/ab0246/a270092/model_codes/awiesm3-develop/oifs-48r1}
OUT=${2:-$(cd "$(dirname "$0")/../.." && pwd)/merge_request_descriptions.md}
BASE=upstream/main
FIRST_MR=91
BRANCHES=(pr1/landice-ism-coupling pr2/coupling-fixes-and-co2-spinup
          pr3/ocean-skin-and-surface-fixes pr4/snow-cover-depletion
          pr5/dms-marine-ccn pr6/namelist-tuning-exposures
          pr7/lpjg-raupach-roughness)

g () { git -C "$REPO" "$@"; }
stat_of () { g diff --numstat "$1" "$2" -- ifs-source | awk '{f++; a+=$1; d+=$2} END{printf "%d files, +%d / -%d", f, a, d}'; }

# does the branch's own commit apply cleanly straight onto vendor main?
WT=$(mktemp -d); g worktree add -q --detach "$WT" "$BASE" 2>/dev/null
declare -A SA
for b in "${BRANCHES[@]}"; do
  git -C "$WT" checkout -q --detach "$BASE"
  if git -C "$WT" cherry-pick -x "origin/$b" >/dev/null 2>&1; then SA[$b]=yes; else SA[$b]=no; fi
  git -C "$WT" cherry-pick --abort 2>/dev/null
done
g worktree remove --force "$WT" 2>/dev/null; g worktree prune

note_for () {   # $1 index, $2 standalone
  local i=$1 sa=$2 mr=$((FIRST_MR-1+i)) prev=$((FIRST_MR-2+i)) last=$((FIRST_MR-1+${#BRANCHES[@]}))
  if [ "$i" -eq 1 ]; then
    echo "Part of an ordered stack of ${#BRANCHES[@]} merge requests, !$FIRST_MR through !$last, splitting the AWI-ESM3 tuning campaign branch. This is the first of them, and it also applies cleanly to \`main\` on its own."
  elif [ "$sa" = yes ]; then
    echo "Part of an ordered stack of ${#BRANCHES[@]} merge requests, !$FIRST_MR through !$last. It is branched from !$prev, but it also applies cleanly to \`main\` on its own, so it does not have to wait for !$prev to land."
  else
    echo "Part of an ordered stack of ${#BRANCHES[@]} merge requests, !$FIRST_MR through !$last. It is branched from !$prev and overlaps its predecessors textually, so !$FIRST_MR through !$prev need to land first. Until they do, the diff shown here against \`main\` includes their changes as well."
  fi
}

{
cat <<'HDR'
# OpenIFS 48r1 upstreaming: merge request titles and descriptions

The campaign branch `movcav-landice+co2-concdriven` split into seven stacked branches on `git.smhi.se/jan.streffing/oifs48r1`, to be merged into `ec-earth/vendor/openifs/oifs48r1` `main` in the order below. Each branch is cut from the one before it, so a branch's own change is its diff against its predecessor, and its diff against `main` shrinks to that once its predecessors have landed.

The "standalone" column records whether the branch's own commit also applies cleanly straight onto `main`. The first two do, so they can be reviewed and merged in either order. The rest overlap their predecessors textually, almost always in `surfece.F90`, and need the stack order.

Merge request numbers assume the six still to be opened take !92 through !97 in order, which holds only if nobody else opens one on this project first. !91 is confirmed. Check the numbers before pasting, and fix the stack sentence at the top of any description whose number came out different.

Copy the fenced block for each merge request straight into the GitLab description field. The fences are not part of the text: everything inside one is the description, already in GitLab-flavoured markdown.

| # | MR | branch | own diff | standalone |
|---|----|--------|----------|------------|
HDR
prev=$BASE; i=0
for b in "${BRANCHES[@]}"; do
  i=$((i+1))
  printf '| %s | !%s | `%s` | %s | %s |\n' "$i" "$((FIRST_MR-1+i))" "$b" "$(stat_of "$prev" "origin/$b")" "${SA[$b]}"
  prev="origin/$b"
done
echo
i=0
for b in "${BRANCHES[@]}"; do
  i=$((i+1)); mr=$((FIRST_MR-1+i))
  echo; echo "---"; echo
  echo "## $i. \`$b\`  (!$mr)"; echo
  echo "Commit \`$(g rev-parse --short origin/$b)\`, author $(g log -1 --format='%an' origin/$b)."
  [ "$i" -eq 1 ] && { echo; echo "Already open as !$FIRST_MR. Replace its title and description with the text below, because GitLab does not re-read either when a branch is force-pushed."; }
  echo; echo "**Title**"; echo; echo '```'
  g log -1 --format=%s "origin/$b"
  echo '```'; echo; echo "**Description**"; echo; echo '```markdown'
  note_for $i "${SA[$b]}"
  echo
  # linkify each subsumed sha so a reviewer can click straight through to the
  # original on the fork; the commit message keeps them plain for `git log`
  g log -1 --format=%b "origin/$b" | awk 'BEGIN{RS="\0"} {sub(/\n+$/,"")} 1' \
    | sed -E 's|^- `([0-9a-f]{7})` |- [`\1`](https://git.smhi.se/jan.streffing/oifs48r1/-/commit/\1) |'
  echo '```'
done
} > "$OUT"
printf 'wrote %s: %s lines, %s fences, %s sections, %s em dashes\n' \
  "$OUT" "$(wc -l < "$OUT")" "$(grep -c '^```' "$OUT")" "$(grep -c '^## ' "$OUT")" "$(grep -c '—' "$OUT")"
