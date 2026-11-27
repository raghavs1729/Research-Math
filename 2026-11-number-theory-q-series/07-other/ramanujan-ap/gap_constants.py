#!/usr/bin/env python3
"""Constants quoted in section 6 (the gap), certified in Arb where they enter an implication."""
from math import factorial
from flint import arb, ctx
ctx.prec = 200
for X in (10 ** 14,):
    L = arb(X).log()
    Wlow = (L / 2).exp() / L ** 4 * (1 + 14 / L + 145 / L ** 2)                 # Lemma A(ii) (L >= 31)
    aplus = L * (sum(arb(factorial(k)) / L ** (k + 1) for k in range(14)) + arb(5) * factorial(14) / L ** 15)
    two_lam = 2 * aplus; sqe_mu = arb(1).exp().sqrt() * L / (L - 1)
    print(f"x = 1e14: L = {L.str(8)}, W >= {Wlow.str(8)}, 2 lambda <= {two_lam.str(8)}, sqrt(e) mu = {sqe_mu.str(8)}")
    assert L > 31 and two_lam < arb("2.07") and sqe_mu < arb("1.71") and Wlow > arb("14.57")
    up, dn = arb("3"), arb("4.8")       # hypothesis (H): -4.8 < delta < 3 on [1e14/e, X0]
    tot = arb("2.07") * up + arb("1.71") * dn + dn * dn / arb(X).sqrt()
    print(f"   2.07 * 3 + 1.71 * 4.8 + 4.8^2/sqrt(x) = {tot.str(8)} < 14.57: {bool(tot < arb('14.57'))}")
    assert tot < arb("14.57")
# descriptive: how much of the worst-case budget the known zeros account for
import os
here = os.path.dirname(os.path.abspath(__file__))
for name, fn in (("chi4", "../../prototype/explore/chi4_zeros.txt"), ("chi3", "../../prototype/explore/chi3_zeros.txt"), ("zeta (to 400)", "zeros_zeta.txt")):
    zs = [float(l.split()[0]) for l in open(os.path.join(here, fn)) if l.strip()]
    print(f"{name}: {len(zs)} zeros up to {zs[-1]:.2f}, sum 2/gamma = {sum(2 / z for z in zs):.3f}")
import math
for X in (1e14, 2.04e18):
    L = math.log(X)
    print(f"x = {X:.3g}: phi L^2/(8 pi) = {2 * L * L / (8 * math.pi):.1f} (phi = 2)")
