#!/bin/bash
# Second stage (after run_zeros_queue.sh was stopped at q = 56 for lack of time):
# the remaining moduli q <= 60 and a longer scan ("pass 2": |y| <= 40000, proposal step with the
# zeros below height 60, certificate with all located zeros) for near misses at q = 28, 36.
run() { while [ $(pgrep -fc "t2_zeros.py") -ge 4 ]; do sleep 5; done; "$@" & sleep 1; }
run bash -c 'python3 t2_zeros.py 57 100 2000 0.05 16,25 > raw/t2_zeros_57.jsonl 2> raw/t2_zeros_57.log'
run bash -c 'python3 t2_zeros.py 59 100 2000 0.05 26 > raw/t2_zeros_59.jsonl 2> raw/t2_zeros_59.log'
run bash -c 'python3 t2_zeros.py 28 300 40000 0.05 9,25 60 > raw/t2_zeros_28_pass2.jsonl 2> raw/t2_zeros_28_pass2.log'
run bash -c 'python3 t2_zeros.py 36 300 40000 0.05 13 60 > raw/t2_zeros_36_pass2.jsonl 2> raw/t2_zeros_36_pass2.log'
wait
# A third pass for (28,9) was launched by hand afterwards:
#   python3 t2_zeros.py 28 300 300000 0.05 9 40 > raw/t2_zeros_28_pass3.jsonl 2> raw/t2_zeros_28_pass3.log
