#!/usr/bin/env python3
"""Re-decide in Arb the near-ties ("A" records) reported by apscan.
Each record: A p c n m decision rho-1 [start].  We certify the sign of
n^2 phi log p - e p m at 256 bits and compare with the decision taken in C."""
import sys, gzip
from flint import arb, ctx
ctx.prec = 256
CPHI = [1, 2,2, 2,2, 4,4,4,4, 4,4,4,4, 4,4,4,4]
E = arb(1).exp()
bad = tot = 0
for fn in sys.argv[1:]:
    op = gzip.open if fn.endswith('.gz') else open
    for line in op(fn, 'rt'):
        w = line.split()
        if not w or w[0] != 'A': continue
        p, c, n, m, dec = int(w[1]), int(w[2]), int(w[3]), int(w[4]), int(w[5])
        v = arb(n) ** 2 * CPHI[c] * arb(p).log() - E * p * m
        assert v > 0 or v < 0
        true = 1 if v > 0 else 0
        tot += 1
        rel = (v / (E * p * m)).str(8)
        ok = (true == dec)
        bad += (not ok)
        print(("OK   " if ok else "WRONG"), line.strip(), " true rho-1 =", rel)
print(f"# near-ties: {tot}, decided wrongly by the C code: {bad}")
sys.exit(1 if bad else 0)
