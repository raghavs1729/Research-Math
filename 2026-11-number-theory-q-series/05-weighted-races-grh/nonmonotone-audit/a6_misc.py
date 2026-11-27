#!/usr/bin/env python3
"""
a6_misc.py -- assorted audit checks.
 1. exact containment of the audited enclosures (raw_density_q*.txt) in my
    independent enclosures (a3_runs.log), and the comparisons on my enclosures.
 2. Gauss sums tau(chi) = i sqrt(q) (root number +1) for q = 163, 232; numerical
    functional-equation test Lambda(s) = Lambda(1-s).
 3. sup_x sqrt(x)|J0(x)| on a grid (vs sqrt(2/pi) and 1); -log J0(1) - 1/4 < 0.018.
 4. zeros_q*.txt == zeros list in zeros_q*.state.json; block log consistency.
"""
import re, json, math, cmath
from fractions import Fraction
from flint import arb, acb, ctx, dirichlet_char
import numpy as np
from scipy.special import j0

A = "/home/user/weighted-prime-races/monotonicity/nonmonotone-audit"
S = "/home/user/weighted-prime-races/monotonicity/nonmonotone"

# 1
theirs = {}
for q in (163, 232):
    for m in re.finditer(r"DELTA q=(\d+) M=(\S+) .*?lower=(\S+)\s+upper=(\S+)", open("%s/raw_density_q%d.txt" % (S, q)).read()):
        theirs[(int(m.group(1)), Fraction(m.group(2)))] = (Fraction(m.group(3)), Fraction(m.group(4)))
mine = {}
for line in open("%s/a3_runs.log" % A):
    m = re.match(r"q=(\d+) M=(\S+) N=(\d+) .*delta in \[(\S+), (\S+)\]", line)
    if m:
        key = (int(m.group(1)), Fraction(m.group(2)).limit_denominator(1000), int(m.group(3)))
        mine[key] = (Fraction(m.group(4)), Fraction(m.group(5)))
allok = True
for (q, M, N), (lo, hi) in sorted(mine.items()):
    t = theirs.get((q, M))
    if t:
        inside = lo <= t[0] and t[1] <= hi
        allok &= inside
        print("q=%d M=%s: audited [%s, %s] inside mine(N=%d) [%.12f, %.12f]: %s"
              % (q, M, float(t[0]), float(t[1]), N, float(lo), float(hi), inside))
print("ALL AUDITED ENCLOSURES INSIDE MINE:", allok)
for q, a, b in ((163, Fraction(3, 5), Fraction(6, 5)), (232, Fraction(18, 25), Fraction(1)), (163, Fraction(57, 100), Fraction(117, 100))):
    for N in (133, 114, 100, 80):
        if (q, a, N) in mine and (q, b, N) in mine:
            print("MINE q=%d N=%d: upper delta(%s) = %.12f < lower delta(%s) = %.12f : %s ; gap >= %.3e"
                  % (q, N, a, float(mine[(q, a, N)][1]), b, float(mine[(q, b, N)][0]),
                     mine[(q, a, N)][1] < mine[(q, b, N)][0], float(mine[(q, b, N)][0] - mine[(q, a, N)][1])))

# 2
ctx.prec = 100
for q, l in ((163, 162), (232, 115)):
    chi = dirichlet_char(q, l)
    tau = sum((chi(n) * acb.exp_pi_i(acb(2 * n) / q) for n in range(1, q)), acb(0))
    print("q=%d: tau(chi) = %s ; i sqrt(q) = %s" % (q, tau.str(15), acb(0, arb(q).sqrt()).str(15)))
    def Lam(s):
        s = acb(s)
        return (arb(q) / arb.pi()) ** ((s + 1) / 2) * ((s + 1) / 2).gamma() * chi.l(s)
    s = acb(0.3, 2.7)
    print("   Lambda(0.3+2.7i) = %s ; Lambda(0.7-2.7i) = %s" % (Lam(s).str(12), Lam(1 - s).str(12)))

# 3
x = np.linspace(1e-6, 2000, 4_000_001)
v = np.sqrt(x) * np.abs(j0(x))
print("max sqrt(x)|J0(x)| on (0,2000] grid = %.6f  (sqrt(2/pi) = %.6f)" % (v.max(), math.sqrt(2 / math.pi)))
ctx.prec = 64
print("-log J0(1) - 1/4 =", (-(arb(1).bessel_j(0)).log() - arb(1) / 4).str(12), "(< 0.018)")

# 4
for q in (163, 232):
    st = json.load(open("%s/zeros_q%d.state.json" % (S, q)))
    fz = [l.split()[0] for l in open("%s/zeros_q%d.txt" % (S, q)) if l.strip() and not l.startswith("#")]
    print("q=%d: file zeros == state zeros: %s ; n = %d ; blocks %s" % (q, fz == st["zeros"], len(fz), st["blocks"]))
