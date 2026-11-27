#!/usr/bin/env bash
# Poor man's statement comparison (the real tool is github.com/leanprover/comparator, see
# comparator.json): print the elaborated type of every headline theorem, and the body of every
# definition used in the statements, once from the challenge module and once from the solution
# module, and diff the two listings.  Run from the project root after
#     lake build AllxChallenge AllxSolution
set -euo pipefail
cd "$(dirname "$0")/.."
names=$(python3 -c "import json;print(' '.join(json.load(open('comparator.json'))['theorem_names']))")
defs=$(python3 -c "import json;print(' '.join(json.load(open('comparator.json'))['definition_names']))")
tmp=$(mktemp -d)
for side in AllxChallenge AllxSolution; do
  {
    echo "import $side"
    echo "set_option pp.all false"
    for n in $names; do echo "#check @$n"; done
    for d in $defs; do echo "#print $d"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/AllxChallenge.out" "$tmp/AllxSolution.out" > /dev/null; then
  echo "OK: the $(echo $names | wc -w) headline statements and the $(echo $defs | wc -w) definitions elaborate identically in AllxChallenge and AllxSolution"
else
  echo "MISMATCH between challenge and solution:"; diff "$tmp/AllxChallenge.out" "$tmp/AllxSolution.out"; exit 1
fi
