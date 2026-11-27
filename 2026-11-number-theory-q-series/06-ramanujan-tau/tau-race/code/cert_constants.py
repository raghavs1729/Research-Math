#!/usr/bin/env python3
"""Certified constants for the tau race (Arb ball arithmetic, python-flint).

  (C1) sigma0^2 = 2 * Lambda'/Lambda(1)  for  Lambda(s) = (2pi)^-s Gamma(s+11/2) L(s,Delta)
       (equal to sum_rho 1/|rho|^2 under GRH(Delta));  a0 = (L'/L)(0, Delta).
  (C2) (L_S'/L_S)(2) and Lambda_S'/Lambda_S(2) for S = Sym^2 Delta (Euler product, exact tau(p),
       Deligne tail), and the resulting bound C_S for sum_rho 2/|rho(rho+2)| under GRH(Sym^2).
  (C3) c_psi = sup_{0<y<=Y1} (y - psi(y))/sqrt(y)  and  theta(y) < y on (0, Y1].

Every number printed with "in [lo, hi]" is a rigorous enclosure given the correctness of
FLINT/Arb and of this program.  Usage: cert_constants.py [P] [Y1]   (defaults 10^7, 10^7)
"""
import os
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
import sys, json, hashlib, time
from flint import arb, acb, ctx, fmpz, fmpq
from tau_table import tau_list

ctx.prec = 200

from arbutil import encl
def lohi(x, digits=18):
    """outward-rounded decimal strings enclosing the arb x"""
    return encl(x, digits)

# ---------------------------------------------------------------- (C1)
def completed_logderiv(tau, NT=40):
    """Lambda*(w) = int_0^oo Delta(iy) y^(w-1) dy = (2pi)^-w Gamma(w) sum tau(n) n^-w,  w = s + 11/2.
    Lambda*(w) = sum_n tau(n) [E(w,n) + E(12-w,n)],  E(w,n) = int_1^oo e^{-2 pi n y} y^{w-1} dy
    (Delta(i/y) = y^12 Delta(iy)).  Derivative in w:
    Lambda*'(w) = sum_n tau(n) int_1^oo e^{-2 pi n y} log y (y^{w-1} - y^{11-w}) dy.
    Evaluated at w = 13/2 (s = 1).  Returns (Lambda*, Lambda*') as arbs."""
    two_pi = 2*arb.pi()
    w = arb(13)/2
    val = arb(0); der = arb(0)
    for n in range(1, NT+1):
        a = two_pi*n
        # E(w,n) = a^-w Gamma(w, a)
        E1 = a**(-w) * a.gamma_upper(w)                 # Gamma(w, a): z.gamma_upper(s)
        E2 = a**(-(12-w)) * a.gamma_upper(12-w)
        val += tau[n]*(E1+E2)
        # derivative term by rigorous integration on [1, Y] plus tail
        Y = arb(1 + (160 + n - 1)//n)      # exact integer endpoint; e^{-2 pi n (Y-1)} <= e^{-1005}
        def f(y, analytic, a=a):
            ly = y.log(analytic=analytic)
            return (-(a*y)).exp() * ly * ((ly*(arb(11)/2)).exp() - (ly*(arb(9)/2)).exp())
        I = acb.integral(f, acb(1), acb(Y), rel_tol=arb(2)**(-170), abs_tol=arb(2)**(-400))
        assert I.imag.contains(0)
        I = I.real
        # tail: int_Y^oo e^{-a y} log y (y^5.5 - y^4.5) dy <= int_Y^oo e^{-a y} y^6.5 dy
        #     <= e^{-aY} Y^6.5 int_0^oo e^{-a v} (1+v/Y)^6.5 dv <= e^{-aY} Y^6.5 / (a - 6.5/Y)
        tail = (-(a*Y)).exp() * Y**arb("6.5") / (a - arb("6.5")/Y)
        I = I + arb(0, tail.upper())
        der += tau[n]*I
    # tail in n > NT: |tau(n)| <= d(n) n^5.5 <= n^6.5;
    # |E1+E2| <= 2 int_1^oo e^{-a y} y^5.5 dy <= 2 e^{-a}/(a-5.5);  |B_n| <= e^{-a}/(a-6.5)   (n >= 2)
    t1 = arb(0); t2 = arb(0)
    # sum_{n>NT} n^6.5 e^{-2 pi n} * 2/(2 pi n - 6.5) <= sum n^6.5 e^{-2 pi n}  <= geometric bound
    # ratio of consecutive terms ((n+1)/n)^6.5 e^{-2pi} <= (1+1/NT)^6.5 e^{-2 pi} =: r < 1
    r = (1+arb(1)/NT)**arb("6.5") * (-two_pi).exp()
    first = arb(NT+1)**arb("6.5") * (-(two_pi*(NT+1))).exp()
    geo = first/(1-r)
    val += arb(0, geo.upper()); der += arb(0, geo.upper())
    return val, der

def c1(tau):
    val, der = completed_logderiv(tau)
    ld = der/val                               # Lambda'/Lambda at s = 1
    sigma0sq = 2*ld
    # L'/L(1) (analytic normalisation) = ld + log(2pi) - psi(13/2);  a0 = L'/L(0) = log 2pi - psi(11/2) - ld
    lp1 = ld + (2*arb.pi()).log() - (arb(13)/2).digamma()
    a0 = (2*arb.pi()).log() - (arb(11)/2).digamma() - ld
    L1 = val * (2*arb.pi())**(arb(13)/2) / (arb(13)/2).gamma()
    return dict(Lambda_at_1=val, sigma0sq=sigma0sq, LprimeL_1=lp1, a0=a0, L_at_1=L1)

# ---------------------------------------------------------------- (C2)
def primes_upto(N):
    import numpy as np
    s = np.ones(N+1, dtype=bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]

def c2(tau, primes, P):
    """-L_S'/L_S(2) = sum_p log p * sum_{k>=1} (alpha^{2k}+1+beta^{2k}) p^{-2k}
       = sum_p log p * [ (c2 x - 2 x^2)/(1 - c2 x + x^2) + x/(1-x) ],  x = p^-2, c2 = tau(p)^2/p^11 - 2 (exact rational)."""
    S = arb(0)
    for p in primes:
        if p > P: break
        x = fmpq(1, p*p)
        c2_ = fmpq(tau[p]*tau[p], p**11) - 2
        term = (c2_*x - 2*x*x)/(1 - c2_*x + x*x) + x/(1-x)
        S += arb(p).log() * arb(term)
    # tail: |sum_k (..) x^k| <= 3x/(1-x) = 3/(p^2-1);  sum_{p>P} log p/(p^2-1) <= 1.01 (log P + 1)/P   (P >= 10)
    tail = arb("3.03")*(arb(P).log()+1)/P
    S += arb(0, tail.upper())
    LS2 = -S
    pi = arb.pi()
    gS2 = -arb(3)/2*pi.log() + ((arb(3)/2).digamma() + (arb(13)/2).digamma() + arb(7).digamma())/2
    LamS2 = gS2 + LS2
    # C_S = sum_rho 2/|rho(rho+2)| <= (12/5) * Lambda_S'/Lambda_S(2)  under GRH(Sym^2) (see proof, Lemma S)
    CS = arb(12)/5 * LamS2
    gS0 = -arb(3)/2*pi.log() + ((arb(1)/2).digamma() + (arb(11)/2).digamma() + arb(6).digamma())/2
    # K_S >= (L_S'/L_S)(-2) = -gammaS'/gammaS(-2) - gammaS'/gammaS(3) - (L_S'/L_S)(3),  |(L_S'/L_S)(3)| <= -3 zeta'/zeta(3)
    def gS(s):
        s = arb(s)
        return -arb(3)/2*pi.log() + (((s+1)/2).digamma() + ((s+11)/2).digamma() + ((s+12)/2).digamma())/2
    from flint import acb_series
    z = acb_series([3, 1], prec=3).zeta()
    zlog3 = -(z[1]/z[0]).real                  # -zeta'/zeta(3) > 0
    KS = -gS(-2) - gS(3) + 3*zlog3
    # a_{-1} = (L'/L)(-1, Delta) = 2 log 2pi - psi(9/2) - psi(15/2) - (L'/L)(2, Delta);
    # -(L'/L)(2,Delta) = sum_p log p (lam x - 2 x^2)/(1 - lam x + x^2), x = p^-2, lam = tau(p)/p^5.5; tail <= 2.02 (log P+1)/P
    S1 = arb(0)
    for p in primes:
        if p > P: break
        ap = arb(p); x = 1/(ap*ap); lam = arb(tau[p])/(ap**5*ap.sqrt())
        S1 += ap.log()*(lam*x - 2*x*x)/(1 - lam*x + x*x)
    tail1 = arb("2.02")*(arb(P).log()+1)/P
    S1 += arb(0, tail1.upper())
    am1 = 2*(2*pi).log() - (arb(9)/2).digamma() - (arb(15)/2).digamma() + S1
    return dict(LS_logderiv_2=LS2, gammaS_logderiv_2=gS2, LambdaS_logderiv_2=LamS2, C_S=CS,
                gammaS_logderiv_0=gS0, three_LamS2=3*LamS2, aS0_lower=-gS0-3*LamS2,
                K_S=KS, a_minus1=am1, a_minus1_minus_2over9=am1-arb(2)/9, L_logderiv_2=-S1)

# ---------------------------------------------------------------- (C3)
def c3(primes, Y1):
    """sup over 0 < y <= Y1 of (y - psi(y))/sqrt(y), and check theta(y) < y.
    On [n, n') (consecutive prime powers) psi is constant and (y - psi)/sqrt(y) increases, so the sup
    over the interval is the limit at n'."""
    import math
    # list of prime powers <= Y1 with log p
    pp = []
    for p in primes:
        if p > Y1: break
        q = p
        while q <= Y1:
            pp.append((q, p, q == p)); q *= p
    pp.sort()
    ctx.prec = 100
    psi = arb(0); theta = arb(0)
    best = arb(2).sqrt()            # interval (0,2): sup of y/sqrt(y) = sqrt 2
    best_at = 2
    theta_ok = True
    logs = {}
    nxt = [q for (q, _, _) in pp[1:]] + [Y1]
    for (q, p, isprime), q2 in zip(pp, nxt):
        lp = logs.get(p)
        if lp is None:
            lp = arb(p).log()
            if p < 3200: logs[p] = lp
        psi += lp
        if isprime:
            theta += lp
            if not (theta < q): theta_ok = False
        v = (q2 - psi)/arb(q2).sqrt()
        if not (v < best):
            if v > best: best, best_at = v, q2
            else:        best = best.union(v)
    ctx.prec = 200
    return dict(c_psi=best, c_psi_at=best_at, theta_lt_y=theta_ok, psi_Y1=psi, theta_Y1=theta, Y1=Y1)

if __name__ == "__main__":
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 10**7
    Y1 = int(sys.argv[2]) if len(sys.argv) > 2 else 10**7
    t0 = time.time()
    tau = tau_list(max(P, 100)+1)
    primes = primes_upto(max(P, Y1))
    out = {}
    r1 = c1(tau)
    for k, v in r1.items(): print("C1 %-14s in [%s, %s]" % (k, *lohi(v, 30))); out[k] = lohi(v, 30)
    r2 = c2(tau, primes, P)
    for k, v in r2.items(): print("C2 %-22s in [%s, %s]" % (k, *lohi(v, 14))); out[k] = lohi(v, 14)
    r3 = c3(primes, Y1)
    print("C3 c_psi in [%s, %s] attained as y -> %d^-" % (*lohi(r3["c_psi"], 12), r3["c_psi_at"]))
    print("C3 theta(y) < y for all primes y <= Y1:", r3["theta_lt_y"], " psi(Y1) =", r3["psi_Y1"].str(15), " theta(Y1) =", r3["theta_Y1"].str(15))
    out.update(c_psi=lohi(r3["c_psi"], 12), c_psi_at=r3["c_psi_at"], theta_lt_y=r3["theta_lt_y"], P=P, Y1=Y1)
    assert r3["theta_lt_y"]
    assert r2["aS0_lower"] > 0 and r2["a_minus1_minus_2over9"] > 0
    js = json.dumps(out, indent=1, sort_keys=True)
    open(DATA + "/constants.json", "w").write(js + "\n")
    print("sha256(constants.json) =", hashlib.sha256((js+"\n").encode()).hexdigest())
    print("time %.1f s" % (time.time()-t0))
