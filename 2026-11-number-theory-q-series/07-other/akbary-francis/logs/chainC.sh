#!/bin/bash
# after chainA: prime sweep for (3), a = 1, prime moduli 11..149
cd /home/user/weighted-prime-races/monotonicity/akbary-francis
while pgrep -f "chainA.sh" > /dev/null; do sleep 5; done
mkdir -p sieve_out/sweep
: > sieve_out/sweep_a1.txt
for q in 11 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149; do
  lc=$(nice -n 10 python3 mertens_const.py $q 1 | head -n 1 | sed 's/.*\[\([-0-9.]*\) .*/\1/')
  X=1e9; if [ $q -eq 11 ] || [ $q -eq 19 ]; then X=1e10; fi
  nice -n 10 ./sieve $q 1 $X $lc sieve_out/sweep/samples_q${q}.txt sieve_out/sweep/records_q${q}.txt sieve_out/sweep/small_q${q}.txt
  first=$(grep -m1 FAIL sieve_out/sweep/records_q${q}.txt)
  echo "q=$q logC=$lc X=$X first: ${first:-none} | $(tail -n 1 sieve_out/sweep/records_q${q}.txt)" >> sieve_out/sweep_a1.txt
  rm -f sieve_out/sweep/small_q${q}.txt
done
