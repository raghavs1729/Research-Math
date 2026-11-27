#!/usr/bin/env bash
# Poor man's statement comparison (the real tool is github.com/leanprover/comparator, see
# comparator.json): print the elaborated type of every headline theorem once from Challenge and
# once from Solution, and diff the two listings.  Run from the project root after
#   lake build RaceChallengeDeps RaceChallenge RaceSolution
set -euo pipefail
cd "$(dirname "$0")/.."
names=$(python3 -c "import json;print(' '.join(json.load(open('comparator.json'))['theorem_names']))")
tmp=$(mktemp -d)
for side in RaceChallenge RaceSolution; do
  {
    echo "import $side"
    for n in $names; do echo "#check @$n"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/RaceChallenge.out" "$tmp/RaceSolution.out" > /dev/null; then
  echo "OK: the $(echo $names | wc -w) headline statements elaborate to identical types in Challenge and Solution"
else
  echo "MISMATCH between Challenge and Solution:"; diff "$tmp/RaceChallenge.out" "$tmp/RaceSolution.out"; exit 1
fi
