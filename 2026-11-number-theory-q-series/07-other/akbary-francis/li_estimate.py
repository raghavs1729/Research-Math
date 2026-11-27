#!/usr/bin/env python3
"""
li_estimate.py q [threshold] -- HEURISTIC size estimate (not used in any proof).
If the positive ordinates were linearly independent over Q (LI), S_1(u) would be
equidistributed like X = sum_j 2 w_j cos(theta_j), theta_j i.i.d. uniform.  Chernoff:
  P(X < -c) <= exp(-I(c)),   I(c) = sup_{lam>0} [ lam c - sum_j log I_0(2 w_j lam) ],
and the bound is sharp to leading exponential order.  Zeros above the file height
are replaced by their total weight (deterministic worst case: X - tail).
Prints I(c)/log(10): the expected log10 of the u-scale at which S(u) < -c first occurs.
"""
import sys, math
import numpy as np
from scipy.special import ive
from scipy.optimize import minimize_scalar
from afcommon import load_zeros
from closed_form import F_q
from flint import ctx

q = int(sys.argv[1])
c = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
ctx.prec = 100
F = float(F_q(q)[0].mid())
g = np.array([float(d) for _, d in load_zeros(q)])
w = 2.0 / (0.25 + g * g)
tail = F - w.sum()
print("q=%d: F=%.5f, listed weight %.5f, tail %.5f, %d ordinates" % (q, F, w.sum(), tail, len(g)))


def logI0(x):
    return np.log(ive(0, x)) + np.abs(x)


for cc in (c - 0.2, c - 0.1, c, c + 0.1, c + 0.2, c + 0.3):
    for label, extra in (("tail random (ignored)", 0.0), ("tail worst case", tail)):
        target = cc + extra         # need listed part < -(c + tail) in the worst case
        if target >= w.sum():
            print("  S < -%.2f [%s]: impossible with listed zeros" % (cc, label))
            continue
        f = lambda lam: -(lam * target - logI0(w * lam).sum())
        r = minimize_scalar(f, bounds=(0, 1e7), method="bounded", options={"xatol": 1e-6})
        print("  S < -%.2f [%s]: rate I = %.2f  => frequency ~ 10^-%.1f (lambda* = %.1f)"
              % (cc, label, -r.fun, -r.fun / math.log(10), r.x))
