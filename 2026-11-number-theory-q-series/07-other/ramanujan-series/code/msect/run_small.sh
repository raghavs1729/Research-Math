#!/bin/sh
# multisection scan (q = 1, 2, 3) of Guillera's targets and Sun's A2, A3 over the seeds with at most 3 parameters
for s in Dixon3F2 Watson3F2 Whipple3F2 Seed1 Seed2 Seed3 Seed4 Seed9 Seed10 Kummer2F1 Gauss2nd Bailey2F1; do
  python3 multisect_scan.py $s 4 12 ../survey/verify_G.py msect/hits_G_small.jsonl
  python3 multisect_scan.py $s 4 12 ../survey/verify_B.py msect/hits_B_small.jsonl A2,A3
done
