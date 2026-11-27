#!/usr/bin/env bash
# Statement comparison (a weak substitute for github.com/leanprover/comparator, which is not
# installed here).  Three checks:
#  1. textual: every theorem of PiSeriesStatement/*.lean (docstring + statement up to `:= by`)
#     occurs byte for byte in PiSeriesSolution.lean, and the solution side has no `sorry`,
#     `axiom` or `native_decide` anywhere (PiSeriesSolution.lean and PiSeries/);
#  2. the statement files import only Mathlib and contain no `def`, `axiom`, `notation`, `macro`;
#  3. elaborated: `#check @name` prints the same type from both modules.
# Run from the project root after `lake build`.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.elan/bin:$PATH"
names=$(cat PiSeriesStatement/*.lean | grep -oP '^theorem \K\w+')
python3 - <<'PY'
import re,sys,glob
so=open('PiSeriesSolution.lean').read()
pat=re.compile(r"(/--(?:.|\n)*?-/\ntheorem (\w+)(?:.|\n)*?) := by\n")
total=0
for f in sorted(glob.glob('PiSeriesStatement/*.lean')):
    st=open(f).read()
    imports=re.findall(r'^import (.*)$', st, re.M)
    if imports!=['Mathlib']: sys.exit(f"{f}: imports {imports}, expected only Mathlib")
    body=re.sub(r'/-(?:.|\n)*?-/','',st)
    for bad in ('def ','axiom ','notation','macro','instance','abbrev','opaque','local ','attribute'):
        if re.search(r'^\s*(noncomputable\s+)?'+bad, body, re.M): sys.exit(f"{f}: contains `{bad.strip()}`")
    thms=pat.findall(st)
    if not thms: sys.exit(f"{f}: no theorems found")
    bad=[n for t,n in thms if t not in so]
    if bad: sys.exit("textual MISMATCH: "+", ".join(bad))
    total+=len(thms)
for f in ['PiSeriesSolution.lean']+sorted(glob.glob('PiSeries/**/*.lean', recursive=True))+['PiSeries.lean']:
    src=re.sub(r'/-(?:.|\n)*?-/','',open(f).read())
    src=re.sub(r'--.*','',src)
    for bad in ('sorry','native_decide','axiom ','unsafe ','implemented_by','extern'):
        if bad in src: sys.exit(f"{f} contains `{bad.strip()}`")
print(f"OK (textual): {total} statements of PiSeriesStatement/*.lean occur verbatim in PiSeriesSolution.lean; statement files import only Mathlib and define nothing; no sorry/axiom/native_decide on the solution side")
PY
tmp=$(mktemp -d)
for side in PiSeriesStatement PiSeriesSolution; do
  {
    echo "import $side"
    for n in $names; do echo "#check @$n"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/PiSeriesStatement.out" "$tmp/PiSeriesSolution.out" > /dev/null; then
  echo "OK (elaborated): the $(echo $names | wc -w) statements have identical types in Statement and Solution:"
  cat "$tmp/PiSeriesSolution.out"
else
  echo "MISMATCH between Statement and Solution:"; diff "$tmp/PiSeriesStatement.out" "$tmp/PiSeriesSolution.out"; exit 1
fi
