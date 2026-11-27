#!/bin/bash
cd /home/user/weighted-prime-races/monotonicity/akbary-francis
while pgrep -f "chainB.sh" > /dev/null; do sleep 5; done
nice -n 10 python3 fq_table.py 400 > logs/fq_table.log 2>&1
