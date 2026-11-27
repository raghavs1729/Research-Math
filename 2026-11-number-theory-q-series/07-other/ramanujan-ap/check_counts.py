#!/usr/bin/env python3
"""Cross-checks of the pass-1 tables (via the CUM lines of apstep):
 (1) pi(x) at EVERY superchunk boundary x = K*10^10 against primecount (independent
     combinatorial algorithm, Deleglise-Rivat / Gourdon);
 (2) pi(10^n; q, a) for q = 3, 4, 5, 8, 12 against the OEIS values in oeis_reference.txt.
Usage: check_counts.py step_log [every]     (every: check primecount at every `every`-th boundary)
"""
import sys, subprocess
CQ = [1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12]
CA = [0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11]
every = int(sys.argv[2]) if len(sys.argv) > 2 else 1
cum = {}
for line in open(sys.argv[1]):
    w = line.split()
    if w and w[0] == "CUM": cum[int(w[1])] = list(map(int, w[2:]))
xs = sorted(cum)
print("boundaries:", len(xs), "from", xs[0], "to", xs[-1])
# internal consistency: class sums
for x in xs:
    n = cum[x]
    assert n[1] + n[2] + 1 == n[0], x            # mod 3: + prime 3
    assert n[3] + n[4] + 1 == n[0], x            # mod 4: + prime 2
    assert sum(n[5:9]) + 1 == n[0], x            # mod 5: + prime 5
    assert sum(n[9:13]) + 1 == n[0], x           # mod 8: + prime 2
    assert sum(n[13:17]) + 2 == n[0], x          # mod 12: + primes 2, 3
    assert n[13] + n[15] == n[1] and n[14] + n[16] + 1 == n[2]   # 12 -> 3 (prime 2 is 2 mod 3)
    assert n[13] + n[14] == n[3] and n[15] + n[16] + 1 == n[4]   # 12 -> 4 (prime 3 is 3 mod 4)
    assert n[9] + n[11] == n[3] and n[10] + n[12] == n[4]        # 8 -> 4
print("class sums consistent at all boundaries")
bad = 0; nchk = 0
for i, x in enumerate(xs):
    if i % every and x != xs[-1]: continue
    pc = int(subprocess.run(["primecount", str(x)], capture_output=True, text=True, check=True).stdout)
    nchk += 1
    if pc != cum[x][0]:
        bad += 1; print("MISMATCH pi(%d): table %d primecount %d" % (x, cum[x][0], pc))
print("primecount: %d boundaries checked, %d mismatches" % (nchk, bad))
# OEIS
ref = {}
for line in open("oeis_reference.txt"):
    if line.startswith("#"): continue
    A, name, vals = [t.strip() for t in line.split("|")]
    ref[A] = list(map(int, vals.split()))
# class index -> (A-number, offset to add to OEIS value to get our count of primes < 10^n in the class)
M = {3: ("A091098", 0), 4: ("A091099", 0),            # 4k+1, 4k+3
     1: ("A091115", 0), 2: ("A091119", 1),            # 6k+1 = 3k+1 ; 6k+5 plus the prime 2 = 3k+2
     5: ("A073505", 0), 6: ("A073507", 1), 7: ("A073506", 0), 8: ("A073508", 0),   # mod 10: 1; 7 (+2); 3; 9
     9: ("A091126", 0), 10: ("A091127", 0), 11: ("A091128", 0), 12: ("A091129", 0),
     13: ("A091161", 0), 14: ("A091162", 0), 15: ("A091163", 0), 16: ("A091164", 0)}
nbad = 0; nok = 0
for n in range(10, 20):
    x = 10 ** n
    if x not in cum: continue
    for c, (A, off) in M.items():
        want = ref[A][n - 1] + off
        if cum[x][c] != want:
            nbad += 1; print("OEIS MISMATCH n=%d class (%d,%d): table %d, %s+%d = %d" % (n, CQ[c], CA[c], cum[x][c], A, off, want))
        else: nok += 1
    print("10^%d: pi(x;4,1) = %d, pi(x;4,3) = %d, pi(x;3,1) = %d, pi(x;3,2) = %d" % (n, cum[x][3], cum[x][4], cum[x][1], cum[x][2]))
print("OEIS: %d values agree, %d mismatches" % (nok, nbad))
sys.exit(1 if bad or nbad else 0)
