#!/usr/bin/env bash
# Poor man's statement comparison (the real tool is github.com/leanprover/comparator, see
# comparator.json): print the elaborated type of every headline theorem once from the challenge
# module and once from the solution module, and diff the two listings.  Run from the project root
# after    lake build ChebyshevChallenge ChebyshevSolution
set -euo pipefail
cd "$(dirname "$0")/.."
names=$(python3 -c "import json;print(' '.join(json.load(open('comparator.json'))['theorem_names']))")
tmp=$(mktemp -d)
for side in ChebyshevChallenge ChebyshevSolution; do
  {
    echo "import $side"
    for n in $names; do echo "#check @$n"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/ChebyshevChallenge.out" "$tmp/ChebyshevSolution.out" > /dev/null; then
  echo "OK: the $(echo $names | wc -w) headline statements elaborate to identical types in ChebyshevChallenge and ChebyshevSolution"
else
  echo "MISMATCH between challenge and solution:"; diff "$tmp/ChebyshevChallenge.out" "$tmp/ChebyshevSolution.out"; exit 1
fi
