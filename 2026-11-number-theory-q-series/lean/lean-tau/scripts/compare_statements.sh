#!/usr/bin/env bash
# Statement comparison (a weak substitute for github.com/leanprover/comparator, which is not
# installed here).  Checks:
#  1. TauStatement.lean imports only Mathlib; its definition block is byte for byte the block of
#     Tau/Defs.lean (which imports only Mathlib); every theorem of TauStatement.lean (docstring +
#     statement up to `:= by`) occurs byte for byte in TauSolution.lean; the solution side
#     (TauSolution.lean, Tau.lean, Tau/*.lean) has no `sorry`, `axiom`, `native_decide`, `unsafe`,
#     `implemented_by`, `extern`, and TauStatement.lean has no `axiom`;
#  2. elaborated: `#check @name` (theorems) and `#print name` (definitions) give the same output from
#     both modules.
# Run from the project root after `lake build`.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.elan/bin:$PATH"
names=$(grep -oP '^theorem \K\w+' TauStatement.lean)
defs=$(grep -oP '^(noncomputable )?def \K\w+' TauStatement.lean)
python3 - <<'PY'
import re,sys,glob
st=open('TauStatement.lean').read()
so=open('TauSolution.lean').read()
df=open('Tau/Defs.lean').read()
for name,src in (('TauStatement.lean',st),('Tau/Defs.lean',df)):
    imports=re.findall(r'^import (.*)$', src, re.M)
    if imports!=['Mathlib']: sys.exit(f"{name}: imports {imports}, expected only Mathlib")
blk=lambda s: s[s.index('-- BEGIN DEFINITIONS'):s.index('-- END DEFINITIONS')]
if blk(st)!=blk(df): sys.exit("definition blocks of TauStatement.lean and Tau/Defs.lean differ")
pat=re.compile(r"(/--(?:(?!-/)[\s\S])*-/\ntheorem (\w+)[\s\S]*?) := by\n")
thms=pat.findall(st)
if not thms: sys.exit("no theorems found in TauStatement.lean")
bad=[n for t,n in thms if t not in so]
if bad: sys.exit("textual MISMATCH: "+", ".join(bad))
strip=lambda s: re.sub(r'--.*','',re.sub(r'/-(?:.|\n)*?-/','',s))
if re.search(r'\baxiom\b', strip(st)): sys.exit("TauStatement.lean contains `axiom`")
for f in ['TauSolution.lean','Tau.lean']+sorted(glob.glob('Tau/**/*.lean', recursive=True)):
    src=strip(open(f).read())
    for b in ('sorry','native_decide','axiom ','unsafe ','implemented_by','extern'):
        if b in src: sys.exit(f"{f} contains `{b.strip()}`")
print(f"OK (textual): the definition block and the {len(thms)} statements of TauStatement.lean occur verbatim in Tau/Defs.lean and TauSolution.lean; both import only Mathlib resp. the library; no sorry/axiom/native_decide on the solution side")
PY
tmp=$(mktemp -d)
for side in TauStatement TauSolution; do
  {
    echo "import $side"
    echo "open Tau"
    for n in $names; do echo "#check @Tau.$n"; done
    for d in $defs; do echo "#print Tau.$d"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/TauStatement.out" "$tmp/TauSolution.out" > /dev/null; then
  echo "OK (elaborated): the $(echo $names | wc -w) statements and the $(echo $defs | wc -w) definitions are identical in TauStatement and TauSolution:"
  n=$(echo $names | wc -w)
  awk -v n="$n" '/^(@Tau|Tau)\./{c++} c<=n' "$tmp/TauSolution.out"
else
  echo "MISMATCH between Statement and Solution:"; diff "$tmp/TauStatement.out" "$tmp/TauSolution.out"; exit 1
fi
