#!/usr/bin/env python3
"""Split the campaign's OpenIFS commits into an ordered stack of themed branches.

Replaying the 71 commits onto today's vendor main does not work: the branch has
merged vendor forward 18 times and those merges carry the conflict resolutions,
so 34 of 40 cherry-picks collide.  Instead we reconstruct each stage's tree *by
content*.

For every file changed between upstream/main and FINAL we walk the unified diff
and decide, line by line, whether that line exists yet at stage k:

  context line  -> always present
  '+' line      -> present iff its birth stage <= k   (from forward git-blame)
  '-' line      -> still present iff its death stage > k (from per-commit patches)

Stage 7 therefore reproduces FINAL exactly, which is asserted at the end.
"""
import subprocess, sys, os, collections

REPO = sys.argv[1]
ORDERED = sys.argv[2]
BASE = "upstream/main"
FINAL = "_final"

THEME_STAGE = {
    "LANDICE": 1, "MISC": 1,
    "COUPLE": 2, "XIOS": 2,
    "VOSKIN": 3, "SFCFIX": 3,
    "SNOWSCF": 4,
    "DMS": 5,
    "EXPOSE": 6,
    "LPJG": 7,
}

def git(*a, **kw):
    return subprocess.run(("git",) + a, cwd=REPO, capture_output=True, text=True,
                          errors="replace", **kw).stdout

# ---- commit -> stage, in chronological (ordered.txt) order -------------------
commits = []            # [(fullsha, stage, subject)]
for line in open(ORDERED):
    line = line.rstrip("\n")
    if not line:
        continue
    sha, theme, subj = line.split("|", 2)
    full = git("rev-parse", sha).strip()
    commits.append((full, THEME_STAGE[theme], subj))
sha2stage = {c[0]: c[1] for c in commits}
# the one merge commit that contributes surviving lines (conflict resolution);
# its 2 lines in surfece.F90 belong with the landice trunk they resolve.
sha2stage[git("rev-parse", "a242841").strip()] = 1

files = [f for f in git("diff", "--name-only", BASE, FINAL, "--", "ifs-source").split("\n") if f]

# ---- death stage of removed lines: last of our commits whose patch drops it --
def death_stages(path):
    """content -> list of stages that removed a line with that content, in order"""
    out = collections.defaultdict(list)
    for full, stage, _ in commits:
        patch = git("show", "--format=", "--unified=0", full, "--", path)
        for ln in patch.split("\n"):
            if ln.startswith("-") and not ln.startswith("---"):
                out[ln[1:]].append(stage)
    return out

# ---- birth stage of surviving lines: forward blame on FINAL -----------------
def birth_stages(path):
    """1-based new-file line number -> stage (0 = vendor)"""
    res = {}
    pc = git("blame", "--line-porcelain", FINAL, "--", path)
    cur = None
    for ln in pc.split("\n"):
        if len(ln) > 40 and ln[:40].isalnum() and " " in ln:
            parts = ln.split()
            if len(parts) >= 3 and len(parts[0]) == 40:
                cur = parts[0]
                res[int(parts[2])] = sha2stage.get(cur, 0)
    return res

STAGES = range(1, 8)
os.makedirs(f"{REPO}/.stagegen", exist_ok=True)
report = []

for path in files:
    old = git("show", f"{BASE}:{path}")
    new = git("show", f"{FINAL}:{path}")
    old_exists = bool(git("cat-file", "-e", f"{BASE}:{path}") or
                      subprocess.run(("git", "cat-file", "-e", f"{BASE}:{path}"),
                                     cwd=REPO, capture_output=True).returncode == 0)
    births = birth_stages(path)
    deaths = death_stages(path)
    deaths_used = collections.defaultdict(int)

    diff = git("diff", "-U1000000", BASE, FINAL, "--", path).split("\n")
    # skip header
    i = 0
    while i < len(diff) and not diff[i].startswith("@@"):
        i += 1
    if i >= len(diff):
        report.append((path, None, "no hunk header (binary/mode?)"))
        continue
    i += 1

    newno = 0
    rows = []   # (kind, text, stage)
    for ln in diff[i:]:
        if not ln:
            continue
        k, txt = ln[0], ln[1:]
        if k == "\\":
            continue
        if k == " ":
            newno += 1
            rows.append((" ", txt, 0))
        elif k == "+":
            newno += 1
            rows.append(("+", txt, births.get(newno, 0)))
        elif k == "-":
            st = deaths.get(txt)
            if st:
                idx = min(deaths_used[txt], len(st) - 1)
                deaths_used[txt] += 1
                d = st[idx]
            else:
                d = 1      # unattributable removal: fold into the trunk
            rows.append(("-", txt, d))

    # A '+' line blamed to a vendor commit is vendor text the diff chose to
    # re-emit (a duplicated/moved block).  It must exist in FINAL, so give it
    # the stage of the nearest neighbouring real addition.
    plus = [i for i, r in enumerate(rows) if r[0] == "+"]
    for pos, i in enumerate(plus):
        if rows[i][2]:
            continue
        st = 0
        for j in plus[pos::-1]:
            if rows[j][2]:
                st = rows[j][2]; break
        if not st:
            for j in plus[pos:]:
                if rows[j][2]:
                    st = rows[j][2]; break
        rows[i] = ("+", rows[i][1], st or 1)

    for k in STAGES:
        buf = []
        for kind, txt, st in rows:
            if kind == " ":
                buf.append(txt)
            elif kind == "+":
                if st and st <= k:
                    buf.append(txt)
            else:  # '-'
                if st > k:
                    buf.append(txt)
        d = f"{REPO}/.stagegen/{k}/{os.path.dirname(path)}"
        os.makedirs(d, exist_ok=True)
        with open(f"{REPO}/.stagegen/{k}/{path}", "w") as fh:
            fh.write("\n".join(buf) + ("\n" if new.endswith("\n") else ""))
    report.append((path, len(rows), None))

print(f"reconstructed {len([r for r in report if r[2] is None])} files")
for p, n, err in report:
    if err:
        print(f"  !! {p}: {err}")
