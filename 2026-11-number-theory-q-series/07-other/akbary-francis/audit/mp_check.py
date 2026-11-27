#!/usr/bin/env python3
"""
mp_check.py q zerofile u0 [a]  -- independent (mpmath only) recomputation of
  F_q (Hadamard closed form, own characters, own L'/L(1) via Hurwitz zeta),
  W = sum over windows of 2/(1/4+d^2), S_T(u0) = sum 2 w cos(d u0 - arg chi(a)),
  the resulting bound S_T + (F_q - W), a rigorous a-priori error bound for using d in
  place of gamma (|gamma - d| < 1e-30), and the sensitivity of S_T to u0.
q prime.  u0 is parsed as an exact decimal (mpf at 80 digits).
Labels n are interpreted as Conrey labels chi_q(n, .) = e(ind(n) ind(m)/(q-1)), g = least
primitive root (verified separately by recert_zeros.py).
"""
import sys
import mpmath as mp

mp.mp.dps = 80
q = int(sys.argv[1]); path = sys.argv[2]; u0s = sys.argv[3]
a = int(sys.argv[4]) if len(sys.argv) > 4 else 1
u0 = mp.mpf(u0s)


def primroot(q):
    for g in range(2, q):
        if len({pow(g, k, q) for k in range(q - 1)}) == q - 1:
            return g


g = primroot(q); ind = {pow(g, k, q): k for k in range(q - 1)}
phi = q - 1


def chi(j, m):
    return mp.expjpi(mp.mpf(2 * j * ind[m % q]) / phi)


# F_q closed form: F(zeta) = 2 + gamma - log(4 pi); for primitive chi mod q, parity kappa:
# F(chi) = log(q/pi) + 2 Re L'/L(1, chi) + digamma((1+kappa)/2)
F = 2 + mp.euler - mp.log(4 * mp.pi)
Fparts = {"zeta": F}
for j in range(1, phi):
    kappa = j % 2  # chi_j(-1) = chi_j(g^{phi/2}) = (-1)^j
    L1 = sum(chi(j, m) * mp.zeta(1 + mp.mpf(0), mp.mpf(m) / q) for m in range(1, q)) if False else None
    # L(s) = q^{-s} sum chi(m) zeta(s, m/q); at s = 1 use L(1) = -(1/q) sum chi(m) digamma(m/q),
    # L'(1) from the derivative of q^{-s} zeta(s, m/q) with the poles cancelling:
    #   zeta(s,a) = 1/(s-1) - digamma(a) + O(s-1)  and  d/ds zeta(s,a) at s=1, finite part = -gamma_1(a).
    L1 = -sum(chi(j, m) * mp.digamma(mp.mpf(m) / q) for m in range(1, q)) / q
    # use a symmetric numerical derivative at high precision instead of Stieltjes constants
    def L(s):
        return q ** (-s) * sum(chi(j, m) * mp.zeta(s, mp.mpf(m) / q) for m in range(1, q))
    h = mp.mpf("1e-20")
    Lp = (L(1 + h) - L(1 - h)) / (2 * h)
    L1b = (L(1 + h) + L(1 - h)) / 2
    assert abs(L1 - L1b) < mp.mpf("1e-30"), (L1, L1b)
    Fj = mp.log(mp.mpf(q) / mp.pi) + 2 * mp.re(Lp / L1) + mp.digamma(mp.mpf(1 + kappa) / 2)
    Fparts[j] = Fj
    F += Fj
print("q = %d, a = %d" % (q, a))
print("F_q (closed form, mpmath) =", mp.nstr(F, 30))
# Euler-Kronecker constant gamma_q = gamma + sum_{chi != chi0} L'/L(1, chi) (paper Table: 2.61076 for 13)

rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
W = mp.mpf(0); S = mp.mpf(0); dS = mp.mpf(0); err = mp.mpf(0)
labs = {}
for lab, d in rows:
    gma = mp.mpf(d)
    w = 2 / (mp.mpf(1) / 4 + gma ** 2)
    if lab == "zeta":
        arg = 0
    else:
        n = int(lab)
        arg = 2 * mp.pi * mp.mpf(ind[n] * ind[a % q] % phi) / phi   # arg chi_n(a)
    W += w
    S += w * mp.cos(gma * u0 - arg)
    dS += -w * gma * mp.sin(gma * u0 - arg)
    # |d/dgamma (w cos(gamma u0 - arg))| <= w (u0 + 2 gamma w/2 ... ) ; crude: w*u0 + |dw/dgamma|
    err += (w * u0 + 2 * gma * w * w / 2) * mp.mpf("1e-30")
    labs[lab] = labs.get(lab, 0) + 1
print("windows", len(rows), labs)
print("W        =", mp.nstr(W, 30))
print("tail F-W =", mp.nstr(F - W, 25))
print("S_T(u0)  =", mp.nstr(S, 30), " (window error bound %s)" % mp.nstr(err, 3))
print("bound    =", mp.nstr(S + F - W, 25), " threshold -4 ; margin", mp.nstr(-4 - (S + F - W), 15))
print("dS_T/du at u0 =", mp.nstr(dS, 10))


def ST(u):
    s = mp.mpf(0)
    for lab, d in rows:
        gma = mp.mpf(d); w = 2 / (mp.mpf(1) / 4 + gma ** 2)
        arg = 0 if lab == "zeta" else 2 * mp.pi * mp.mpf(ind[int(lab)] * ind[a % q] % phi) / phi
        s += w * mp.cos(gma * u - arg)
    return s


if len(sys.argv) > 5:
    mp.mp.dps = 40
    for dl in ["1e-9", "1e-6", "1e-4", "1e-3", "1e-2", "0.05"]:
        for sg in (1, -1):
            v = ST(u0 + sg * mp.mpf(dl))
            print("S_T(u0 %s %s) = %s   (linear pred %s)" % ("+" if sg > 0 else "-", dl, mp.nstr(v, 15),
                  mp.nstr(S + sg * mp.mpf(dl) * dS, 15)))
