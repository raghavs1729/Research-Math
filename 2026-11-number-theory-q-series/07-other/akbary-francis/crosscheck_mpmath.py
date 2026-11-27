#!/usr/bin/env python3
"""
crosscheck_mpmath.py q a u0 -- independent (non-rigorous, mpmath, no Arb in the numerics)
cross-check of the three numerical ingredients:
 (1) the first zeros of every character in zeros_q<q>.txt are zeros of mpmath's
     dirichlet L-function with character values exp(2 pi i e/ord) (and NOT of the
     conjugate character): labelling convention;
 (2) F_q from mpmath: log(q/pi) + 2 Re L'/L(1, chi) + psi((1+k)/2), L and L' at 1 from
     mpmath's Hurwitz zeta at 1 +- 1e-12 (central difference; ~20 correct digits);
 (3) S_T(u0) recomputed with mpmath at 60 digits from the decimal ordinates.
"""
import sys
import mpmath as mp
from afcommon import load_zeros, chi_arg_table
from flint import dirichlet_char   # used ONLY to read off the character table exponents

q, a, u0s = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
mp.mp.dps = 60
zeros = load_zeros(q)


def chi_values(n):
    chi = dirichlet_char(q, n)
    G = chi.group()
    ex = int(G.exponent())
    return [0] + [mp.exp(2j * mp.pi * int(chi.chi_exponent(b)) / ex) for b in range(1, q)], int(chi.parity())


# (1)
first = {}
for lab, d in zeros:
    first.setdefault(lab, []).append(d)
worst = 0
for lab, ds in first.items():
    for d in ds[:3]:
        s = mp.mpc(0.5, mp.mpf(d))
        if lab == "zeta":
            v = abs(mp.zeta(s))
            vc = None
        else:
            vals, _ = chi_values(int(lab))
            v = abs(mp.dirichlet(s, vals))
            vc = abs(mp.dirichlet(s, [mp.conj(z) for z in vals]))
            if not dirichlet_char(q, int(lab)).is_real():
                assert vc > 1e-3, (lab, d, vc)
        worst = max(worst, v)
        assert v < mp.mpf(10) ** (-25), (lab, d, v)
print("(1) first 3 ordinates of each of the %d labels: |L(1/2 + i d, chi)| < 1e-25 with mpmath (max %s); "
      "conjugate character does not vanish there" % (len(first), mp.nstr(worst, 3)))

# (2)  L(s, chi) = q^(-s) sum_a chi(a) zeta(s, a/q) (mpmath Hurwitz zeta), L'(1) by a central difference
mp.mp.dps = 45
F = 2 + mp.euler - mp.log(4 * mp.pi)
h = mp.mpf(10) ** (-12)
hz = {sg: [None] + [mp.zeta(1 + sg * h, mp.mpf(b) / q) for b in range(1, q)] for sg in (1, -1)}
for n in range(2, q):
    vals, k = chi_values(n)
    Lp = mp.power(q, -(1 + h)) * sum(vals[b] * hz[1][b] for b in range(1, q))
    Lm = mp.power(q, -(1 - h)) * sum(vals[b] * hz[-1][b] for b in range(1, q))
    L1 = (Lp + Lm) / 2
    dL = (Lp - Lm) / (2 * h)
    psi = mp.digamma(mp.mpf(1 + k) / 2)
    F += mp.log(q / mp.pi) + 2 * mp.re(dL / L1) + psi
print("(2) F_%d (mpmath, Hurwitz zeta + central difference) = %s" % (q, mp.nstr(F, 20)))

# (3)
mp.mp.dps = 60
tab = chi_arg_table(q, a)
u0 = mp.mpf(u0s)
S = mp.mpf(0)
W = mp.mpf(0)
for lab, d in zeros:
    g = mp.mpf(d)
    w = 2 / (mp.mpf(1) / 4 + g * g)
    ph = tab[lab]
    S += w * mp.cos(g * u0 - 2 * mp.pi * ph.numerator / ph.denominator)
    W += w
print("(3) S_T(u0) (mpmath) = %s ; W = %s ; F - W = %s ; S_T + F - W = %s"
      % (mp.nstr(S, 25), mp.nstr(W, 25), mp.nstr(F - W, 12), mp.nstr(S + F - W, 15)))
