#!/usr/bin/env bash
# Statement comparison (a weak substitute for github.com/leanprover/comparator, which is not
# installed here).  Checks:
#  1. NahmStatement.lean imports only Mathlib; its definition block is byte for byte the block of
#     Nahm/Defs.lean (which imports only Mathlib); every theorem of NahmStatement.lean (docstring +
#     statement up to `:= by`) occurs byte for byte in NahmSolution.lean; the solution side
#     (NahmSolution.lean, Nahm.lean, Nahm/*.lean and the libraries it imports, ../lean-krtw/KRTW/*.lean,
#     ../rr-identities/lean/RRFusion/*.lean) has no `sorry`, `axiom`, `native_decide`, `unsafe`,
#     `implemented_by`, `extern`, and NahmStatement.lean has no `axiom`;
#  2. elaborated: `#check @name` (theorems) and `#print name` (definitions) give the same output from
#     both modules.
# Run from the project root after `lake build`.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.elan/bin:$PATH"
names=$(grep -oP '^theorem \K\w+' NahmStatement.lean)
defs=$(grep -oP '^(noncomputable )?def \K\w+' NahmStatement.lean)
python3 - <<'PY'
import re,sys,glob
st=open('NahmStatement.lean').read()
so=open('NahmSolution.lean').read()
df=open('Nahm/Defs.lean').read()
for name,src in (('NahmStatement.lean',st),('Nahm/Defs.lean',df)):
    imports=re.findall(r'^import (.*)$', src, re.M)
    if imports!=['Mathlib']: sys.exit(f"{name}: imports {imports}, expected only Mathlib")
blk=lambda s: s[s.index('-- BEGIN DEFINITIONS'):s.index('-- END DEFINITIONS')]
if blk(st)!=blk(df): sys.exit("definition blocks of NahmStatement.lean and Nahm/Defs.lean differ")
pat=re.compile(r"(/--(?:(?!-/)[\s\S])*-/\ntheorem (\w+)[\s\S]*?) := by\n")
thms=pat.findall(st)
if not thms: sys.exit("no theorems found in NahmStatement.lean")
bad=[n for t,n in thms if t not in so]
if bad: sys.exit("textual MISMATCH: "+", ".join(bad))
strip=lambda s: re.sub(r'--.*','',re.sub(r'/-(?:.|\n)*?-/','',s))
if re.search(r'\baxiom\b', strip(st)): sys.exit("NahmStatement.lean contains `axiom`")
for f in ['NahmSolution.lean','Nahm.lean']+sorted(glob.glob('Nahm/**/*.lean', recursive=True))+sorted(glob.glob('../lean-krtw/KRTW/**/*.lean', recursive=True))+sorted(glob.glob('../rr-identities/lean/RRFusion/**/*.lean', recursive=True)):
    src=strip(open(f).read())
    for b in ('sorry','native_decide','axiom ','unsafe ','implemented_by','extern'):
        if b in src: sys.exit(f"{f} contains `{b.strip()}`")
print(f"OK (textual): the definition block and the {len(thms)} statements of NahmStatement.lean occur verbatim in Nahm/Defs.lean and NahmSolution.lean; both import only Mathlib resp. the library; no sorry/axiom/native_decide on the solution side")
PY
tmp=$(mktemp -d)
for side in NahmStatement NahmSolution; do
  {
    echo "import $side"
    echo "open Nahm"
    for n in $names; do echo "#check @Nahm.$n"; done
    for d in $defs; do echo "#print Nahm.$d"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/NahmStatement.out" "$tmp/NahmSolution.out" > /dev/null; then
  echo "OK (elaborated): the $(echo $names | wc -w) statements and the $(echo $defs | wc -w) definitions are identical in NahmStatement and NahmSolution:"
  n=$(echo $names | wc -w)
  awk -v n="$n" '/^(@Nahm|Nahm)\./{c++} c<=n' "$tmp/NahmSolution.out"
else
  echo "MISMATCH between Statement and Solution:"; diff "$tmp/NahmStatement.out" "$tmp/NahmSolution.out"; exit 1
fi
