#!/usr/bin/env python3
"""Zero ordinates used by the heuristic model (model.py).  NOT used by any theorem.

zeta: rigorous enclosures of the first NZ zeros (Arb, acb.zeta_zeros).
chi mod 5, 8, 12 (stretch moduli): sign changes of Hardy's Z on a grid, refined by
bisection -- uncertified (a close pair could be missed), adequate for a heuristic.
For complex characters mod 5 (chi and conj chi) zeros are found from sign changes of
the Hardy Z-function of chi on t in (-T, T) (zeros of conj chi are the negatives).
Output: zeros_zeta.txt, zeros_q{5,8,12}.txt  (lines: conductor index gamma)
"""
import sys
from flint import acb, arb, dirichlet_char, ctx
ctx.prec = 100
T = 400.0
def zeta_zeros(T):
    out = []; n = 1
    while True:
        zs = acb.zeta_zeros(n, 50)
        for z in zs:
            g = float(z.imag.mid())
            if g > T: return out
            out.append(g)
        n += 50
def hardy_zeros(q, idx, T, two_sided):
    chi = dirichlet_char(q, idx)
    f = lambda t: float(chi.hardy_z(acb(t)).real.mid())
    h = 0.02
    lo = -T if two_sided else 0.0
    out = []
    t = lo + 1e-9; fp = f(t)
    while t < T:
        t2 = t + h; f2 = f(t2)
        if fp == 0.0 or fp * f2 < 0:
            a, b, fa = t, t2, fp
            for _ in range(60):
                m = 0.5 * (a + b); fm = f(m)
                if fa * fm <= 0: b = m
                else: a, fa = m, fm
            out.append(0.5 * (a + b))
        t, fp = t2, f2
    return out
if __name__ == "__main__":
    with open("zeros_zeta.txt", "w") as f:
        for g in zeta_zeros(T): f.write("%.15f\n" % g)
    for q in (5, 8, 12):
        with open("zeros_q%d.txt" % q, "w") as f:
            for idx in range(2, q):
                from math import gcd
                if gcd(idx, q) != 1: continue
                chi = dirichlet_char(q, idx)
                if not chi.is_primitive(): 
                    f.write("# %d %d imprimitive (conductor %d): zeros of the inducing character are used by model.py\n" % (q, idx, chi.conductor()))
                    continue
                real = chi.is_real()
                for g in hardy_zeros(q, idx, T, not real):
                    f.write("%d %d %d %.12f\n" % (q, idx, 1 if real else 0, g))
