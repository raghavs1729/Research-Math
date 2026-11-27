#!/usr/bin/env python3
"""Consistency of the two certificates on their common range (a cross-check, not a proof).

Direct certificate (../cert/results_q*.txt, results_q*_stretch.txt):  lo <= delta'(M) <= hi < 0 on [c,d].
Envelope certificate (results_q*.txt here):                            -M delta'(M) >= LB > 0   on [a,b].
For every pair of cells whose intersection [x,y] has positive length both statements hold for all
M in [x,y], so necessarily
        LB <= M |delta'(M)| <= x * |lo|                      (taking M = x)
A violation would mean that one of the two certificates is wrong.  The ratio
        LB / (M |delta'(M)|)  lies in  [ LB/(y |lo|), LB/(x |hi|) ]
measures how much of the true derivative the envelope bound captures.

Scope (REFEREE-5 minor 7): on the common range TAILF and the far remainder are far below the
enclosure width of the direct certificate, so this comparison validates MAIN - NUM only; it would
not notice an error in TAILF (the referee's mutant "TAILF dropped" passes it).  TAILF is covered
by the referee's independent evaluation of the cell quantity (REFEREE-5.md section 3.3).
Usage: python3 consistency.py q
"""
import os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
q = int(sys.argv[1])

def parse(path, keys):
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if line.startswith("#") or len(p) < 4 or p[-1] != "PASS":
                continue
            d = dict(x.split("=", 1) for x in p[2:-1] if "=" in x)
            rows.append((Fraction(p[0]), Fraction(p[1])) + tuple(float(d[k]) for k in keys))
    return sorted(rows)

direct = []
for name in ("results_q%d.txt" % q, "results_q%d_stretch.txt" % q):
    direct += parse(os.path.join(HERE, "..", "cert", name), ("lo", "hi"))
direct.sort()
env = parse(os.path.join(HERE, "results_q%d.txt" % q), ("LB",))
npairs = nviol = 0
rmin, rmax = 1e9, 0.0
nec = 0.0
sample = []
j0 = 0
for a, b, LB in env:
    while j0 < len(direct) and direct[j0][1] <= a:
        j0 += 1
    j = j0
    while j < len(direct) and direct[j][0] < b:
        c, d, lo, hi = direct[j]
        x, y = max(a, c), min(b, d)
        if y > x:
            npairs += 1
            if not LB <= float(x) * abs(lo):
                nviol += 1
                print("VIOLATION", a, b, c, d, LB, lo)
            r1, r2 = LB / (float(y) * abs(lo)), LB / (float(x) * abs(hi))
            rmin, rmax = min(rmin, r1), max(rmax, r2)
            nec = max(nec, LB / (float(x) * abs(lo)))
            sample.append((float(x), float(y), LB, lo, hi, r1, r2))
        j += 1
print("# q=%d: direct cells %d on [%s, %s]; envelope cells %d on [%s, %s]"
      % (q, len(direct), float(direct[0][0]), float(direct[-1][1]), len(env), float(env[0][0]), float(env[-1][1])))
print("# overlapping cell pairs: %d ; violations of  LB <= x*|lo| : %d" % (npairs, nviol))
print("# necessary condition LB/(x*|lo|) <= 1:  largest value over all pairs = %.4f" % nec)
print("# LB / (M |delta'|) is enclosed, pair by pair, in [LB/(y|lo|), LB/(x|hi|)]; over all pairs: [%.4f, %.4f]"
      % (rmin, rmax))
print("#   (the upper end can exceed 1 only through the width of the direct enclosure, hi/lo >= 0.92)")
print("#   x         y        LB(env)      lo(direct)     hi(direct)    LB/(y|lo|)  LB/(x|hi|)")
step = max(1, len(sample) // 40)
for s in sample[::step]:
    print("%9.4f %9.4f  %.6e  %.6e  %.6e   %.4f     %.4f" % s)
