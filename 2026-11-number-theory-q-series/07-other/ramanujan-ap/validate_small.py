#!/usr/bin/env python3
"""Compare apscan (C, floating screens) with smallx.py (Python, Arb) on x <= X:
failing primes and component starts per bin, for all 17 classes.
Usage: validate_small.py scan_output X"""
import sys, gzip, collections, math
import smallx
X = int(float(sys.argv[2]))
CQ = [1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12]
CA = [0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11]
bmax = int(math.floor(100 * math.log10(X))) - 1          # complete bins only
scan = collections.defaultdict(lambda: [0, 0])
op = gzip.open if sys.argv[1].endswith('.gz') else open
for line in op(sys.argv[1], 'rt'):
    w = line.split()
    if w and w[0] == 'B' and int(w[4]) <= bmax:
        r = scan[(int(w[3]), int(w[4]))]; r[0] += int(w[6]); r[1] += int(w[12])
P = smallx.primes_upto(X)
bad = 0; totf = tots = 0
for c in range(17):
    comps, nfail, nstart, failing, active = smallx.analyse(CQ[c], CA[c], X, P)
    for b in range(0, bmax + 1):
        # apscan covers primes >= 7 only
        f = sum(1 for p in failing if p >= 7 and int(math.floor(100 * math.log10(p) + 1e-12)) == b) if b < 100 else nfail[b]
        s = nstart[b] if b >= 85 else None
        got = scan.get((c, b), [0, 0])
        if got[0] != f or (s is not None and got[1] != s):
            bad += 1; print("MISMATCH class", c, "bin", b, "scan", got, "python", f, s)
        totf += f; tots += got[1]
    print(f"class ({CQ[c]},{CA[c]}): bins 0..{bmax} agree; failing primes >= 7: {sum(1 for p in failing if p >= 7 and p < 10 ** ((bmax + 1) / 100))}")
print("mismatches:", bad, " total failing primes compared:", totf, " component starts compared:", tots)
sys.exit(1 if bad else 0)
