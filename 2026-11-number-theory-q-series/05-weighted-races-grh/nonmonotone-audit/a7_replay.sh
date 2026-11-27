#!/bin/sh
# a7_replay.sh -- provenance: re-run the audited code from a clean copy and diff.
#  (i) zeros_cert.py 232 100 from scratch -> zeros_q232.txt must be byte-identical
#  (ii) density.py q=163 M=3/5 6/5 -> DELTA lines identical up to the time= field
#  (iii) crosscheck.py on mutant m1 (gamma_1 deleted), Simpson only: does the
#        non-rigorous cross-check notice the incomplete list?
set -u
A=/home/user/weighted-prime-races/monotonicity/nonmonotone-audit
S=/home/user/weighted-prime-races/monotonicity/nonmonotone
R=$A/replay; rm -rf $R; mkdir -p $R; cp $S/*.py $S/zeros_q163.txt $R/
cd $R && nice -n 10 python3 zeros_cert.py 232 100 > zeros_cert_232.log 2>&1; echo "zeros_cert exit=$?"
cmp $R/zeros_q232.txt $S/zeros_q232.txt && echo "zeros_q232.txt BYTE-IDENTICAL"
diff <(grep -o "BLOCK.*W = [^;]*;" $R/zeros_q232.cert.txt) <(grep -o "BLOCK.*W = [^;]*;" $S/zeros_q232.cert.txt) && echo "cert block lines identical (up to timings)"
cd $R && nice -n 10 python3 density.py --q 163 --M 3/5 6/5 > density_163.log 2>&1; echo "density exit=$?"
diff <(grep DELTA $R/density_163.log | sed 's/time=.*//') <(grep -E "M=(3/5|6/5) " $S/raw_density_q163.txt | sed 's/time=.*//') && echo "DELTA lines identical (up to time=)"
cd $A/mutants/m1 && OMP_NUM_THREADS=1 nice -n 10 python3 crosscheck.py 163 3/5 6/5 --skip-mp --mc 0 --nroots 1 2>&1 | grep -E "Simpson|sigma"
