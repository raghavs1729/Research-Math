#!/bin/sh
# Reproduces every result file of this directory (about 4 minutes of CPU, one process).
# Resumable: cells already PASS in a results file are skipped, provided the file header carries
# the SHA-256 of the current prover, shared module and zero file; otherwise the run is refused
# (delete the file to recompute).
set -e
cd "$(dirname "$0")"
nice -n 10 python3 certify_envelope.py selftest > selftest.out
nice -n 10 python3 certify_envelope.py range --q 4 --lo 1 --hi 2000 --ratio 1.01 --out results_q4.txt --workers 1
nice -n 10 python3 certify_envelope.py range --q 3 --lo 1 --hi 2000 --ratio 1.01 --out results_q3.txt --workers 1
nice -n 10 python3 certify_envelope.py point --q 4 --M 1.5 2 3 3.689 4.484 6 8 10 14 17.3 19.2 30 > points_q4.txt
nice -n 10 python3 certify_envelope.py point --q 3 --M 1.5 2 3 3.689 4.484 6 8 10 14 17.3 19.2 30 > points_q3.txt
nice -n 10 python3 consistency.py 4 > consistency_q4.txt
nice -n 10 python3 consistency.py 3 > consistency_q3.txt
nice -n 10 python3 analytic_leg.py > analytic_leg.out
