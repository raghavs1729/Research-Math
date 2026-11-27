#!/usr/bin/env python3
"""
l1const.py -- zero-sum constants for HKMN Problem 23, with NO zero data.

For a primitive Dirichlet character chi of conductor f (f = 1: zeta), ASSUMING
GRH(chi), every nontrivial zero is rho = 1/2 + i*gamma (gamma real), and

    1/|rho(rho+1)| = 1/sqrt(A(A+2)),   A = 1/4 + gamma^2
                   = (1/pi) * int_0^pi  d(theta) / (A + 1 - cos(theta)).

(Elementary: int_0^pi d(theta)/(c - cos(theta)) = pi/sqrt(c^2-1), c > 1.)
Summing over rho (everything is positive, Tonelli):

    S(chi) := sum_rho 1/|rho(rho+1)| = (1/pi) int_0^pi  V_chi(alpha(theta))/alpha(theta)  d(theta),
    alpha(theta) = sqrt(5/4 - cos(theta))  in [1/2, 3/2],

where for real alpha > 0, with s = 1/2 + alpha  (Davenport Ch. 12 (17)-(18); under GRH
Re 1/(s-rho) = alpha/(alpha^2+gamma^2)):

    V_chi(alpha) := sum_rho alpha/(alpha^2+gamma^2)
                  = Re (L'/L)(s,chi) + (1/2) log(f/pi) + (1/2) psi((s+a)/2)        (f > 1, a = parity)
                  = (zeta'/zeta)(s) + 1/(s-1) + 1/s - (1/2) log(pi) + (1/2) psi(s/2)   (f = 1).

So S(chi) is an integral of a closed-form analytic function.  In the code
Re (L'/L) is computed by replacing Re (L'/L)(s,chi) with
(1/2)[(L'/L)(s,chi) + (L'/L)(s,chibar)].

NOTE (what is actually used): the quadrature is the trapezoidal rule with an exact,
elementary error formula (see S_integral); only real points s in [1,2] are evaluated.

Closed-form sandwich used as an independent check (no quadrature):
    1/(A+1)  <=  1/sqrt(A(A+2))  <=  (1/2)(1/A + 1/(A+2)),  hence
    V(sqrt5/2)/(sqrt5/2)  <=  S(chi)  <=  (1/2) [ V(1/2)/(1/2) + V(3/2)/(3/2) ].

Also computed, for the explicit formula of psi_1(x,chi) (HKMN23.md, Lemma E):
    odd  chi:  b  = (L'/L)(0,chi);              b0 = L''(-1)/(2 L'(-1))  [constant term at s=-1]
    even chi:  b  = L''(0)/(2 L'(0)) [constant term of L'/L at s=0];   cm1 = (L'/L)(-1,chi)
    zeta    :  (zeta'/zeta)(0) = log(2 pi),  (zeta'/zeta)(-1).

Output: JSON (decimal strings of ball midpoints and radii) + a readable table.
Usage:  python3 l1const.py [prec_bits]  > l1const.out
"""
import sys, json, math
from flint import acb, arb, acb_series, dirichlet_char, ctx

PREC = int(sys.argv[1]) if len(sys.argv) > 1 else 160
ctx.prec = PREC

CANDIDATE_MODULI = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24, 30, 40, 42, 48, 60]
EXTRA_MODULI = [11, 13]          # for comparison with published values only
ALL_MODULI = CANDIDATE_MODULI + EXTRA_MODULI


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def prim_chars(f):
    """Conrey numbers of the primitive characters mod f."""
    if f == 1:
        return [1]
    return [n for n in range(1, f) if math.gcd(n, f) == 1 and dirichlet_char(f, n).is_primitive()]


def conj_number(f, n):
    return pow(n, -1, f) if f > 1 else 1


def ljet(chi, s0, n):
    """Taylor coefficients c_0..c_{n-1} of L(s,chi) at s0."""
    ser = acb_series([acb(s0), acb(1)], prec=n).dirichlet_l(chi)
    return [ser[i] for i in range(n)]


def logder(chi, s):
    c = ljet(chi, s, 2)
    return c[1] / c[0]


def zeta_logder_deflated(s):
    """zeta'/zeta(s) + 1/(s-1), regular at s = 1."""
    ser = acb_series([acb(s), acb(1)], prec=2).zeta(1, True)   # D(s) = zeta(s) - 1/(s-1)
    D, D1 = ser[0], ser[1]
    w = acb(s) - 1
    return (D + w * D1) / (1 + w * D)


def Vfun(f, n, s):
    """Analytic continuation in s of V_chi(alpha), alpha = s - 1/2 (see docstring)."""
    s = acb(s)
    if f == 1:
        return zeta_logder_deflated(s) + 1 / s - arb.pi().log() / 2 + (s / 2).digamma() / 2
    chi = dirichlet_char(f, n)
    chib = dirichlet_char(f, conj_number(f, n))
    par = chi.parity()
    return (logder(chi, s) + logder(chib, s)) / 2 + (arb(f) / arb.pi()).log() / 2 + ((s + par) / 2).digamma() / 2


NTRAP = 64


def S_integral(f, n):
    """S(chi) = (1/pi) int_0^pi h(theta) dtheta,  h(theta) = W(5/4 - cos theta),
    W(beta) = sum_rho 1/(beta + gamma^2) = V(alpha)/alpha, alpha = sqrt(beta).

    Rigorous evaluation by the trapezoidal rule with an EXACT error formula.  For c > 1,
        1/(c - cos t) = (c^2-1)^(-1/2) * (1 + 2 sum_{k>=1} r^k cos(kt)),  r = c - sqrt(c^2-1).
    Here c = 5/4 + gamma^2 >= 5/4, so 0 < r <= 1/2, and h = a_0 + sum_{k>=1} a_k cos(k theta) with
    a_0 = S(chi) and 0 <= a_k <= 2 * 2^(-k) * S(chi).  The n-interval trapezoidal sum on [0, pi]
    of an even 2pi-periodic function is T_n = a_0 + a_{2n} + a_{4n} + ..., hence
        S <= T_n <= S * (1 + 2/(4^n - 1)),   i.e.   T_n (4^n-1)/(4^n+1) <= S(chi) <= T_n.
    Only values of V at real points s = 1/2 + alpha in [1, 2] are needed.
    """
    nt = NTRAP
    tot = arb(0)
    for i in range(nt + 1):
        th = arb.pi() * i / nt
        al = (arb(5) / 4 - th.cos()).sqrt()
        v = Vfun(f, n, acb(arb(1) / 2 + al)).real / al
        tot += v / 2 if i in (0, nt) else v
    T = tot / nt
    lo = T * (arb(4) ** nt - 1) / (arb(4) ** nt + 1)
    return lo.union(T)


def S_sandwich(f, n):
    r5 = arb(5).sqrt() / 2
    lo = Vfun(f, n, acb(1) / 2 + r5).real / r5
    V1 = Vfun(f, n, acb(1)).real * 2
    V3 = Vfun(f, n, acb(2)).real * 2 / 3
    return lo, (V1 + V3) / 2, V1


def lower_constants(f, n):
    if f == 1:
        ser = acb_series([acb(0), acb(1)], prec=2).zeta()
        z0 = ser[1] / ser[0]
        ser = acb_series([acb(-1), acb(1)], prec=2).zeta()
        zm1 = ser[1] / ser[0]
        return {"parity": 0, "zeta_logder_0": z0, "zeta_logder_m1": zm1}
    chi = dirichlet_char(f, n)
    par = chi.parity()
    c0 = ljet(chi, 0, 3)
    cm = ljet(chi, -1, 3)
    if par == 1:
        b = c0[1] / c0[0]
        b0 = cm[2] / cm[1]            # L''(-1)/(2L'(-1)) = c2/c1
        assert cm[0].contains(0)
        return {"parity": 1, "b": b, "b0": b0}
    else:
        assert c0[0].contains(0)
        b = c0[2] / c0[1]
        cm1 = cm[1] / cm[0]
        return {"parity": 0, "b": b, "cm1": cm1}


def bstr(x, d=30):
    x = arb(x) if not isinstance(x, (arb, acb)) else x
    return x.str(d)


def ball(x):
    """(mid, rad) as decimal strings, rad rounded up."""
    x = x if isinstance(x, arb) else arb(x)
    return [x.mid().str(40, radius=False), x.rad().str(5, radius=False)]


def cball(z):
    z = acb(z)
    return {"re": ball(z.real), "im": ball(z.imag)}


def main():
    conductors = sorted({f for q in ALL_MODULI for f in divisors(q) if f % 4 != 2})
    data = {"prec_bits": PREC, "chars": {}, "moduli": {}}
    Sball = {}
    print("# l1const.py  prec =", PREC, "bits")
    print("# S(chi) = sum_rho 1/|rho(rho+1)| under GRH(chi), by the trapezoidal rule (exact error formula) applied to a closed form")
    print("# f  n  parity order   S(chi) [Arb ball]   closed-form sandwich [lo, hi]   V=sum 1/|rho|^2")
    for f in conductors:
        for n in prim_chars(f):
            S = S_integral(f, n)
            lo, hi, V = S_sandwich(f, n)
            # consistency: the sandwich must contain the integral
            assert (S - lo) > -1e-25 and (hi - S) > -1e-25, (f, n, S, lo, hi)
            lc = lower_constants(f, n)
            Sball[(f, n)] = S
            if f == 1:
                par, order = 0, 1
            else:
                c = dirichlet_char(f, n)
                par, order = c.parity(), c.order()
            print(f"{f:3d} {n:3d} {par} {order:2d}  {S.str(25)}  [{lo.str(12, radius=False)}, {hi.str(12, radius=False)}]  {V.str(18)}")
            entry = {"conductor": f, "conrey": n, "parity": par, "order": order,
                     "S": ball(S), "S_lo_closed": ball(lo), "S_hi_closed": ball(hi), "V": ball(V)}
            for k, v in lc.items():
                if k != "parity":
                    entry[k] = cball(v)
            data["chars"][f"{f}.{n}"] = entry
            sys.stdout.flush()
    print()
    print("# q   phi  N(q,1)  L1(q) [Arb ball]   (2/3)N(q,1)   margin=(2/3)N-L1   rel.margin")
    for q in ALL_MODULI:
        L1 = arb(0)
        for f in divisors(q):
            if f % 4 == 2:
                continue
            for n in prim_chars(f):
                L1 += Sball[(f, n)]
        phi = sum(1 for k in range(1, q + 1) if math.gcd(k, q) == 1)
        N = sum(1 for b in range(q) if (b * b - 1) % q == 0) if q > 1 else 1
        thr = arb(2) * N / 3
        m = thr - L1
        print(f"{q:3d} {phi:3d} {N:2d}  {L1.str(22)}  {thr.str(10, radius=False)}  {m.str(12)}  {(m / thr).str(6, radius=False)}")
        data["moduli"][str(q)] = {"phi": phi, "N1": N, "L1": ball(L1), "margin": ball(m),
                                  "passes_sufficient_criterion": bool(m > 0), "fails_criterion": bool(m < 0)}
    with open("l1const.json", "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
