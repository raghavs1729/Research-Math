#!/bin/bash
# Zero-side witness searches (t2_zeros.py) for the square classes with q <= 120 that have no
# prime-side witness (q = 2 mod 4 is covered by q/2; q = 7, 9, 16 are knife-edge moduli and are
# not attempted; q = 11, 19, 28, 32 and the candidate moduli 21, 40, 48, 60 were run separately
# with the argument "sq" = all square classes).  At most 4 jobs at a time.
# NOTE: this queue was stopped after q = 56 (time); the jobs for q >= 64 were NOT run.
# job = q:T:classes   ->   python3 t2_zeros.py q T 2000 0.05 classes
jobs="36:300:13,25 39:150:22 44:150:25 45:150:34 49:100:44 53:100:28 56:100:25 57:100:16,25 59:100:26 64:100:41 68:100:33 72:100:25,49 75:100:49 80:100:49 83:80:69 84:80:25,37 87:80:28 88:80:9,25 91:80:30 93:80:76 97:80:24 99:80:49,70 103:80:14,81 104:80:9,25 105:80:4,46 107:80:40 108:80:97 109:80:12,81 112:80:81 120:80:49"
for j in $jobs; do
  q=${j%%:*}; rest=${j#*:}; T=${rest%%:*}; cl=${rest##*:}
  while [ $(pgrep -fc "t2_zeros.py") -ge 4 ]; do sleep 5; done
  if [ ! -s raw/t2_zeros_${q}.jsonl ]; then
    python3 t2_zeros.py $q $T 2000 0.05 $cl > raw/t2_zeros_${q}.jsonl 2> raw/t2_zeros_${q}.log &
    sleep 1
  fi
done
wait
