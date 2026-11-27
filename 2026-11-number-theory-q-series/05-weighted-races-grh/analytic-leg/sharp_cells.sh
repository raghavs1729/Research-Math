#!/bin/sh
# Sharp direct-certificate enclosures of delta'(M) on tiny cells at M = 300, 310, 320 (350),
# produced by the independent prover ../cert/certify_derivative.py (tolerance 0.003).
# Used by section [G] of certify_analytic_leg.py for the strong consistency check.
# About 25 s per cell on one core.
set -e
cd "$(dirname "$0")"
out=sharp_cells.txt
: > "$out"
for spec in "4 300 300.001" "4 310 310.001" "4 319.999 320" "3 300 300.001" "3 310 310.001" "3 320 320.001" "3 349.999 350"; do
  set -- $spec
  printf 'q=%s ' "$1" >> "$out"
  nice -n 10 python3 ../cert/certify_derivative.py cell --q "$1" --Ma "$2" --Mb "$3" --tol 0.003 >> "$out"
done
