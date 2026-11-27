#!/bin/bash
cd /home/user/weighted-prime-races/monotonicity/akbary-francis
nice -n 10 ./scan scan_out/scan_in_q13_a12.txt 0 2e7 100 -3.5 -4.0 > scan_out/scan_q13_a12_2e7.out 2> scan_out/scan_q13_a12_2e7.err
nice -n 10 python3 verify_zeros.py 11 --count > logs/verify_q11.log 2>&1
nice -n 10 python3 fq_table.py 400 > logs/fq_table.log 2>&1
