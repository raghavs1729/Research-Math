#!/bin/sh
# multisection scan (blocks of 1, 2, 3 terms) of every surveyed identity over the main seed families
for t in A B C; do
  python3 multisect_scan.py $1 $2 12 ../survey/verify_$t.py msect/all_hits_$1.jsonl
done
