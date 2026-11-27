#!/usr/bin/env python3
"""
f_two_ways.py q -- F_q two ways (prime q):
 (1) closed form (closed_form.py);
 (2) certified partial sum W over the pinned zero file (all zeros below T, with
     conjugates) + the smooth main-term ESTIMATE of the tail from the Riemann-von
     Mangoldt density of zeros (not rigorous):
        tail ~ (1/(pi T)) [ (q-2)(log(qT/2pi) + 1) + (log(T/2pi) + 1) ].
Agreement to ~4 digits confirms formulas (49), (50) of Akbary-Francis and the zero data.
"""
import sys, math
from flint import arb, ctx
from afcommon import load_zeros, EPS
from closed_form import F_q

q = int(sys.argv[1])
ctx.prec = 128
F, per = F_q(q)
zeros = load_zeros(q)
e = arb(EPS)
W = arb(0)
perW = {}
T = max(float(d) for _, d in zeros)
T = math.ceil(T)
for lab, d in zeros:
    g = arb(d) + arb(0, e.upper())
    w = 1 / (arb(1) / 4 + g * g)
    perW[lab] = perW.get(lab, arb(0)) + w
    W += 2 * w
est = ((q - 2) * (math.log(q * T / (2 * math.pi)) + 1) + (math.log(T / (2 * math.pi)) + 1)) / (math.pi * T)
print("q = %d, T = %d, %d windows" % (q, T, len(zeros)))
print("(1) closed form           F_q = %s" % F.str(20))
print("(2) certified partial sum W   = %s" % W.str(20))
print("    exact tail F_q - W        = %s" % (F - W).str(12))
print("    density estimate of tail  = %.6f" % est)
print("    W + estimate              = %.6f   (difference from closed form %.2e)" % (float(W.mid()) + est, float(W.mid()) + est - float(F.mid())))
# per conjugate pair: F(chi) + F(chibar) = 2 F(chi) versus sum over windows of chi and chibar
print("per character: closed form F(chi)   vs   partial sums over |gamma| < T (zeros of chi with gamma>0 and of chibar with gamma>0)")
for lab in sorted(perW, key=lambda s: (s != "zeta", int(s) if s != "zeta" else 0)):
    if lab == "zeta":
        part = 2 * perW[lab]
    else:
        part = perW[lab] + perW[str(pow(int(lab), -1, q))]
    print("   %-4s F(chi) = %s   partial = %s   tail = %.5f" % (lab, per[lab].str(12, radius=False), part.str(12, radius=False), float((per[lab] - part).mid())))
