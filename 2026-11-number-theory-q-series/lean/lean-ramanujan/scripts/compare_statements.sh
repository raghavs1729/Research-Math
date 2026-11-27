#!/usr/bin/env bash
# Statement comparison (a weak substitute for github.com/leanprover/comparator, which is not
# installed here).  Two checks:
#  1. textual: every theorem of RamanujanStatement.lean (docstring + statement up to `:= by`)
#     occurs byte for byte in RamanujanSolution.lean;
#  2. elaborated: `#check @name` prints the same type from both modules.
# Run from the project root after `lake build`.
set -euo pipefail
cd "$(dirname "$0")/.."
names=$(grep -oP '^theorem \K\w+' RamanujanStatement.lean)
python3 - <<'PY'
import re,sys
st=open('RamanujanStatement.lean').read(); so=open('RamanujanSolution.lean').read()
pat=re.compile(r"(/--(?:.|\n)*?-/\ntheorem (\w+)(?:.|\n)*?) := by\n")
thms=pat.findall(st)
if not thms: sys.exit("no theorems found")
bad=[n for t,n in thms if t not in so]
if bad: sys.exit("textual MISMATCH: "+", ".join(bad))
if "sorry" in so: sys.exit("RamanujanSolution.lean contains `sorry`")
print(f"OK (textual): {len(thms)} statements of RamanujanStatement.lean occur verbatim in RamanujanSolution.lean")
PY
tmp=$(mktemp -d)
for side in RamanujanStatement RamanujanSolution; do
  {
    echo "import $side"
    for n in $names; do echo "#check @$n"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/RamanujanStatement.out" "$tmp/RamanujanSolution.out" > /dev/null; then
  echo "OK (elaborated): the $(echo $names | wc -w) statements have identical types in Statement and Solution"
else
  echo "MISMATCH between Statement and Solution:"; diff "$tmp/RamanujanStatement.out" "$tmp/RamanujanSolution.out"; exit 1
fi
