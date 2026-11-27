#!/usr/bin/env python3
"""refine_u.py q a u [halfwidth] -- local minimiser of S_T(u) near u (double precision, all listed zeros);
prints u rounded to 4 decimals, S_T there, and S_T + tail estimate.  Exploratory; certify.py does the proof."""
import sys, math
import numpy as np
from flint import ctx
from afcommon import load_zeros, chi_arg_table
from closed_form import F_q
q, a, u0 = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
hw = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
ctx.prec = 100
F = float(F_q(q)[0].mid())
tab = chi_arg_table(q, a)
Z = load_zeros(q)
g = np.array([float(d) for _, d in Z]); ph = np.array([2 * math.pi * float(tab[l]) for l, _ in Z])
w = 2 / (0.25 + g * g)
S = lambda u: float((w * np.cos(g * u - ph)).sum())
us = u0 + np.arange(-int(hw / 1e-4), int(hw / 1e-4) + 1) * 1e-4
v, u = min((S(x), x) for x in us)
u = round(u, 4)
print("%.4f  S_T = %.6f  S_T + tail = %.6f  (tail = %.6f)" % (u, S(u), S(u) + F - w.sum(), F - w.sum()))
