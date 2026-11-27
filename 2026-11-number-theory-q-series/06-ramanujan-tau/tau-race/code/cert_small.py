#!/usr/bin/env python3
"""Certified finite range for the tau race (exact integers tau(p); Arb balls; every comparison certain).

With lam(p) = tau(p) p^(-11/2):
  Theta(t)   = sum_{p<=t} lam(p) log p,      Theta_1(t) = sum_{p<=t} lam(p) log p (t - p)   (piecewise linear)
  D_0, D_1   = sum lam(p), sum lam(p) p,     R(x) = D_0(x) - D_1(x)/x = sum_{p<=x} lam(p)(1 - p/x)
Checks, for N given (default 10^7):
  (A) Theta_1(p) > 0 at every prime 5 <= p <= N        [=> Theta_1 > 0 on (t*, N], t* the zero in (3,5)]
  (B) R(p) > 0 at every prime 7 <= p <= N              [=> R > 0 on (x*, N],  x* the zero in (5,7)]
  (C) lower bounds for G = Theta_1/t^{3/2} and N_R = R log x/sqrt x over all REAL arguments in [a, p_last], several a
      (minimum of the gap bounds over all gaps [q, p] of consecutive primes with p > a, q >= 11)
  (D) the constants A_1(T), A_2(T) of the kernel identity for several T
  (E) t*, x* enclosed.
Output: ../data/small_range_<N>.json (+ sha256).
"""
import os
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
import sys, json, hashlib, time
from flint import arb, ctx
from tau_table import tau_list
from arbutil import encl, dec_down, dec_up, primes_upto

def main(N):
    ctx.prec = 128
    t0 = time.time()
    tau = tau_list(N+1)
    primes = primes_upto(N)
    Th = arb(0); Mm = arb(0); D0 = arb(0); D1 = arb(0)
    failsA = []; failsB = []
    thresholds = [11, 100, 1000, 10**4, 10**5, 10**6]
    infG = {a: None for a in thresholds if a <= N}; infN = {a: None for a in thresholds if a <= N}
    Tlist = [p for p in [97, 199, 499, 997, 9973, 99991, 999983] if p <= N]
    AT = {}
    prev = None      # (p, Theta_1(p), R(p))
    first = {}
    absTh_max = arb(0)     # sup |Theta(p)|/sqrt(p)
    for p in primes:
        ap = arb(p)
        lam = arb(tau[p]) / (ap**5 * ap.sqrt())
        assert abs(lam) < 2                      # Deligne (sanity; certain comparison)
        lg = ap.log()
        Th += lam*lg; Mm += lam*lg*p; D0 += lam; D1 += lam*p
        T1 = ap*Th - Mm
        R = D0 - D1/ap
        if p <= 13: first[p] = dict(Theta=encl(Th, 12), Theta_1=encl(T1, 12), R=encl(R, 12), D0=encl(D0, 12), D1=encl(D1, 12))
        if p == 3: t_star = (lam*lg*3 + Mm - lam*lg*3) / Th        # zero of Theta(3) t - Mm(3) : t* = Mm(3)/Theta(3)
        if p == 5: x_star = D1/D0
        if p >= 5 and not (T1 > 0): failsA.append(p)
        if p >= 7 and not (R > 0): failsB.append(p)
        if prev is not None and p >= 11:
            q, T1q, Rq = prev
            # on [q, p]: Theta_1 linear => >= min(T1q, T1); t^{3/2} <= p^{3/2}
            gG = (T1q.min(T1)) / (ap*ap.sqrt()) if q >= 5 else None
            gN = (Rq.min(R)) * lg/ap.sqrt() if q >= 7 else None     # log x/sqrt x decreasing for x > e^2
            for a in infG:
                if p > a:            # the gap [q, p] meets [a, N]; its bound is valid on all of [q, p]
                    if gG is not None: infG[a] = gG if infG[a] is None else infG[a].min(gG)
                    if gN is not None: infN[a] = gN if infN[a] is None else infN[a].min(gN)
        if p in Tlist:
            A1 = D0 - Th/lg - T1/(ap*lg*lg)
            A2 = -D1 + (ap*Th - T1)/lg + T1/(lg*lg)
            AT[p] = dict(A1=encl(A1, 12), A2=encl(A2, 9), A1T_plus_A2=encl(A1*p + A2, 9),
                         Theta=encl(Th, 9), Theta_1=encl(T1, 6), D0=encl(D0, 9), D1=encl(D1, 6))
        prev = (p, T1, R)
    out = dict(N=N, n_primes=len(primes), failsA=failsA, failsB=failsB,
               t_star=encl(t_star, 15), x_star=encl(x_star, 15), first=first,
               infG={str(a): dec_down(v, 9) for a, v in infG.items() if v is not None},
               infN={str(a): dec_down(v, 9) for a, v in infN.items() if v is not None},
               kernel_constants={str(k): v for k, v in AT.items()},
               final=dict(Theta=encl(Th, 9), Theta_1=encl(T1, 4), D0=encl(D0, 9), D1=encl(D1, 4)))
    js = json.dumps(out, indent=1, sort_keys=True) + "\n"
    open(DATA + "/small_range_%d.json" % N, "w").write(js)
    print(js)
    print("PASS" if not failsA and not failsB else "FAIL", " sha256 =", hashlib.sha256(js.encode()).hexdigest(), " time %.1f s" % (time.time()-t0))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10**7)
