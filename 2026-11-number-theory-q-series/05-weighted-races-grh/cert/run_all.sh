#!/bin/sh
# Reproduces every result file in this directory (about 45 minutes on 3 cores).
# Usage: sh run_all.sh      (resumable: cells already certified in a result file are skipped)
set -e
cd "$(dirname "$0")"
python3 certify_derivative.py selftest
python3 certify_derivative.py aux --out aux_results.txt
# priority range (closes the gap between the coupling theorem and thm:monotone-large)
python3 certify_derivative.py range --q 4 --lo 3.6  --hi 17.4 --ratio 1.005 --tol 0.02 --workers 3 --out results_q4.txt
python3 certify_derivative.py range --q 3 --lo 4.4  --hi 19.2 --ratio 1.005 --tol 0.02 --workers 3 --out results_q3.txt
# stretch range (direct certificate on the range of prop:criterion-eval; coarser cells)
python3 certify_derivative.py range --q 4 --lo 17.4 --hi 320  --ratio 1.01 --tol 0.04 --workers 3 --out results_q4_stretch.txt
python3 certify_derivative.py range --q 3 --lo 19.2 --hi 350  --ratio 1.01 --tol 0.04 --workers 3 --out results_q3_stretch.txt
python3 crosscheck.py --q 4 > crosscheck_q4.txt
python3 crosscheck.py --q 3 > crosscheck_q3.txt
