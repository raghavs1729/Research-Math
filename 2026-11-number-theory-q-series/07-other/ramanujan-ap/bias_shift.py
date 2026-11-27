#!/usr/bin/env python3
"""Heuristic size of the bias effect on the last counterexample (LI model, model.py).
For each q and c in {0, c_res}: X_alpha = the X with m_c(X) = alpha, where
m_c(X) = int_{log X}^infty P(Y >= W(e^u) + c B(u)) du  (expected log-measure of later failures).
Prints X_alpha and the ratio X_alpha(c = 0)/X_alpha(c = c_res), next to the observed ratios."""
import math, numpy as np
import model
def m_curve(q, c):
    r, tv, info = model.amps(q)
    us = np.arange(math.log(1e9), math.log(1e16), 0.02)
    lp = np.array([min(model.tail_prob(r, tv, model.W_of(math.exp(u)) + c * model.B_of(math.exp(u)))[1], 0.0) for u in us])
    cum = np.cumsum((10 ** lp)[::-1])[::-1] * 0.02
    return us, cum
def X_at(us, cum, alpha):
    i = int(np.argmax(cum < alpha))
    # log-linear interpolation
    t = (math.log(cum[i - 1]) - math.log(alpha)) / (math.log(cum[i - 1]) - math.log(cum[i]))
    return math.exp(us[i - 1] + t * 0.02)
for q, cres in ((3, 2), (4, 2), (5, 2), (8, 4), (12, 4)):
    u0, c0 = m_curve(q, 0); u1, c1 = m_curve(q, cres)
    out = []
    for alpha in (1e-2, 1e-3, 1e-4):
        a, b = X_at(u0, c0, alpha), X_at(u1, c1, alpha)
        out.append(f"alpha={alpha:g}: X(c=0) = {a:.3g}, X(c={cres}) = {b:.3g}, ratio {a/b:.2f}")
    print(f"q = {q}: " + "; ".join(out))
