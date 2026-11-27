#!/usr/bin/env python3
"""Referee 3: the Bailey-basis facts (B-i)-(B-v) of F-all-k Sec. 3, the v_a rules, (3.1) and its T# analogue,
checked pointwise with exact rationals at rational (q, a) or (q, x); and where n >= 1 is needed.

(B-i)  sum_{r<=R} a^r q^{r^2} phi_e(r;a)/(q)_{R-r} = a^e q^{e^2} phi_e(R;a)               (exact, finite)
(B-ii) sum_{r>=0} a^r q^{r^2} phi_e(r;a) = a^e q^{e^2}/(aq;q)_inf                        (numerically, 120 digits)
(B-iii) h_a phi_e = (q^-e + a q^e) phi_e                                                   (exact)
(B-iv) phi_d(.;aq^c) = sum_e Gamma^{(c)}_{d,e}(a) phi_e(.;a), sum Gamma = (aq)_c            (exact, pointwise)
(B-v)  tau phi_n = alpha phi_{n-1} + beta phi_n + gamma phi_{n+1} (n >= 1)                (exact, pointwise)
v_a rules for F with kappa_0 = 0; failure of the pointwise h_a tau identity and of finiteness when kappa_0 != 0.
(3.1)  T Xi_{e,sigma} = w(e,sigma) sum_S (Q;Q)^{-1}_{sigma-S} sum_{e'} Gamma^{(2sigma-2S)}_{e,e'}(x q^{2S}) Xi_{e',S}
and the displayed formula for T# Xi#_{e,rho} (Sec. 5).
"""
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 120
q0, a0, x0 = Fr(2, 7), Fr(3, 11), Fr(3, 5)


def poch(a, p, n):
    r = Fr(1)
    for i in range(n):
        r *= 1 - a * p ** i
    return r


def qbin(n, k, q):
    if k < 0 or k > n:
        return Fr(0)
    return poch(q, q, n) / (poch(q, q, k) * poch(q, q, n - k))


def phi(e, R, a, q):
    if R < e:
        return Fr(0)
    return 1 / (poch(q, q, R - e) * poch(a * q, q, R + e))


def Gam(c, d, e, a, q):
    if not (d <= e <= d + c):
        return Fr(0)
    return ((1 - a * q ** (2 * e)) / (1 - a)) * (poch(a, q, d + e) / poch(a * q ** (c + 1), q, d + e)) * (-1) ** (e - d) \
        * q ** ((e - d) * (e - d - 1) // 2) * a ** (e - d) * q ** ((e + d) * (e - d)) * qbin(c, e - d, q)


def expand(f, a, q, Rb):
    c = []
    for R in range(Rb + 1):
        v = f(R) - sum(c[e] * phi(e, R, a, q) for e in range(R))
        c.append(v / phi(R, R, a, q))
    return c


def run():
    q, a = q0, a0
    ok = True
    Rb = 14
    # (B-i)
    okb1 = all(sum(a ** r * q ** (r * r) * phi(e, r, a, q) / poch(q, q, R - r) for r in range(R + 1))
               == a ** e * q ** (e * e) * phi(e, R, a, q) for e in range(5) for R in range(Rb))
    # (B-ii)
    okb2 = True
    for e in range(5):
        s = sum(a ** r * q ** (r * r) * phi(e, r, a, q) for r in range(0, 40))
        inf = mp.mpf(1)
        t = mp.mpf(a.numerator) / a.denominator * mp.mpf(q.numerator) / q.denominator
        qq = mp.mpf(q.numerator) / q.denominator
        while abs(t) > mp.mpf(10) ** -130:
            inf *= 1 - t
            t *= qq
        rhs = mp.mpf((a ** e * q ** (e * e)).numerator) / (a ** e * q ** (e * e)).denominator / inf
        okb2 &= abs(mp.mpf(s.numerator) / s.denominator - rhs) < mp.mpf(10) ** -100
    # (B-iii)
    okb3 = True
    for e in range(5):
        for R in range(Rb):
            lhs = q ** (-R) * (phi(e, R, a, q) - (phi(e, R - 1, a, q) if R >= 1 else 0)) + a * q ** R * phi(e, R, a, q)
            okb3 &= lhs == (q ** (-e) + a * q ** e) * phi(e, R, a, q)
    # (B-iv)
    okb4 = True
    for c in range(6):
        for d in range(5):
            okb4 &= sum(Gam(c, d, e, a, q) for e in range(d, d + c + 1)) == poch(a * q, q, c)
            for R in range(Rb):
                okb4 &= phi(d, R, a * q ** c, q) == sum(Gam(c, d, e, a, q) * phi(e, R, a, q) for e in range(d, d + c + 1))
    okb4 &= all(phi(d, R, a * q * q, q) == poch(a * q, q, 2) * phi(d + 1, R + 1, a, q) for d in range(5) for R in range(Rb))
    # (B-v)
    okb5 = True
    for n in range(1, 6):
        al = q ** (n - 1) / ((1 - a * q ** (2 * n - 1)) * (1 - a * q ** (2 * n)))
        be = -q ** (n - 1) * (1 + a * q ** (2 * n)) / ((1 - a * q ** (2 * n - 1)) * (1 - a * q ** (2 * n + 1)))
        ga = a * q ** (3 * n) / ((1 - a * q ** (2 * n)) * (1 - a * q ** (2 * n + 1)))
        okb5 &= al + be + ga == 0
        for R in range(Rb):
            okb5 &= q ** R * phi(n, R + 1, a, q) == al * phi(n - 1, R, a, q) + be * phi(n, R, a, q) + ga * phi(n + 1, R, a, q)
        c1 = expand(lambda R: phi(n, R + 1, a, q), a, q, Rb)
        c2 = expand(lambda R: q ** (2 * R) * phi(n, R + 1, a, q), a, q, Rb)
        okb5 &= all(c1[e] == 0 for e in range(Rb + 1) if abs(e - n) > 1) and sum(c1) == 1
        okb5 &= all(c2[e] == 0 for e in range(Rb + 1) if abs(e - n) > 1) and sum(c2) == 0
    print(f"(B-i) exact: {okb1}; (B-ii) 100 digits: {okb2}; (B-iii) exact: {okb3}; (B-iv) Gamma formula, sum, "
          f"phi_d(aq^2)=(aq)_2 phi_(d+1)(R+1): {okb4}; (B-v) incl. 'likewise' sums 1, 0: {okb5}")
    ok &= okb1 and okb2 and okb3 and okb4 and okb5
    # v_a rules for F = sum_{n>=1} kappa_n phi_n
    kap = {1: Fr(3, 2), 2: Fr(-5, 7), 4: Fr(2, 9)}
    F = lambda R: sum(k * phi(n, R, a, q) for n, k in kap.items())
    tauF = lambda R: q ** R * F(R + 1)
    hF = lambda G: (lambda R: q ** (-R) * (G(R) - (G(R - 1) if R >= 1 else 0)) + a * q ** R * G(R))
    vF = sum(kap.values())
    c_tau = expand(tauF, a, q, Rb)
    c_htau = expand(hF(tauF), a, q, Rb)
    c_sh = expand(lambda R: F(R + 1), a, q, Rb)
    c_q2 = expand(lambda R: q ** (2 * R) * F(R + 1), a, q, Rb)
    okv = sum(c_tau) == 0 and sum(c_htau) == (1 - 1 / q) * vF and sum(c_sh) == vF and sum(c_q2) == 0
    okv &= all(c_tau[e] == 0 for e in range(7, Rb + 1))
    pw = all(hF(tauF)(R) == F(R + 1) + a * q ** (2 * R) * F(R + 1) - F(R) / q for R in range(Rb))
    okv &= pw
    print(f"v_a rules (kappa_0 = 0): v(tau F) = 0, v(h tau F) = (1-1/q) v(F), v(F(.+1)) = v(F), v(q^2R F(.+1)) = 0, "
          f"pointwise h tau F identity: {okv}")
    ok &= okv
    # with kappa_0 != 0: tau phi_0 has an infinite expansion, and the pointwise identity fails at R = 0
    c0 = expand(lambda R: q ** R * phi(0, R + 1, a, q), a, q, Rb)
    nz = sum(1 for v in c0 if v != 0)
    G0 = lambda R: phi(0, R, a, q)
    t0 = lambda R: q ** R * G0(R + 1)
    fail0 = hF(t0)(0) != G0(1) + a * G0(1) - G0(0) / q
    print(f"kappa_0 != 0: tau phi_0 has {nz} nonzero coefficients among e <= {Rb} (not a finite combination: "
          f"v_a undefined) ; pointwise h tau identity fails at R = 0: {fail0}  -> n >= 1 is needed exactly there")
    # (3.1): T Xi_{e,sigma}, position space vs the displayed formula
    x = x0
    Q = q * q
    w = lambda r, s: x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s)
    ok31 = True
    for e in range(3):
        for sig in range(3):
            for S in range(sig + 1):
                for R in range(10):
                    lhs = sum(w(r, sig) * phi(e, r, x * q ** (2 * sig), q) / poch(q, q, R - r) for r in range(R + 1)) \
                        / poch(Q, Q, sig - S)
                    rhs = w(e, sig) / poch(Q, Q, sig - S) * sum(
                        Gam(2 * sig - 2 * S, e, ep, x * q ** (2 * S), q) * phi(ep, R, x * q ** (2 * S), q)
                        for ep in range(e, e + 2 * (sig - S) + 1))
                    ok31 &= lhs == rhs
    # T# Xi#_{e,rho}(r,s) = lam_{rho+e,e} sum_r (q)^{-1}_{rho-r} sum_e' Gamma^{Q,(rho-r)}_{e,e'}(x^2 q^{2r}) Xi#_{e',r}(r,s)
    oks = True
    for e in range(3):
        for rho in range(3):
            for r in range(rho + 1):
                for s in range(10):
                    lhs = sum(w(rho, S) * phi(e, S, x * x * q ** (2 * rho), Q) / poch(Q, Q, s - S) for S in range(s + 1)) \
                        / poch(q, q, rho - r)
                    lamv = x ** (rho + 2 * e) * q ** ((rho + e) ** 2 + e * e)
                    rhs = lamv / poch(q, q, rho - r) * sum(
                        Gam(rho - r, e, ep, x * x * q ** (2 * r), Q) * phi(ep, s, x * x * q ** (2 * r), Q)
                        for ep in range(e, e + rho - r + 1))
                    oks &= lhs == rhs
    print(f"(3.1) T Xi expansion: {ok31}; Sec. 5 T# Xi# expansion: {oks}")
    ok &= ok31 and oks
    print("ALL BAILEY-BASIS CHECKS " + ("PASSED" if ok else "FAILED"))


if __name__ == "__main__":
    run()
