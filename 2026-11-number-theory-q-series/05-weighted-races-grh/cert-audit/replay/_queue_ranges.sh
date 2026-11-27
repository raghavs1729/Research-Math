#!/bin/sh
# auditor's background queue: wait for plain_eval, then replay both ranges with ONE worker
cd "$(dirname "$0")"
while pgrep -f "plain_eval.py --q" >/dev/null; do sleep 5; done
date
python3 -B ../../cert/certify_derivative.py range --q 4 --lo 3.6 --hi 17.4 --ratio 1.005 --tol 0.02 --workers 1 --out results_q4.txt
date
python3 -B ../../cert/certify_derivative.py range --q 3 --lo 4.4 --hi 19.2 --ratio 1.005 --tol 0.02 --workers 1 --out results_q3.txt
date
