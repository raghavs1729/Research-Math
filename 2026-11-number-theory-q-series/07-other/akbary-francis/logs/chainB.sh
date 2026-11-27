#!/bin/bash
cd /home/user/weighted-prime-races/monotonicity/akbary-francis
while pgrep -f "scan_in_q11_a1.txt" > /dev/null; do sleep 5; done
for a in 10 12; do
  nice -n 10 ./scan scan_out/scan_in_q13_a$a.txt 0 5e7 100 -3.5 -4.0 > scan_out/scan_q13_a${a}_5e7.out 2> scan_out/scan_q13_a${a}_5e7.err
done
nice -n 10 python3 verify_zeros.py 11 --count > logs/verify_q11.log 2>&1
