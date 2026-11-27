#!/usr/bin/env python3
"""AUDITOR: numerical sanity of (H2)+(H3) -- NOT part of any proof.

Under GRH + completeness of the list up to g_N, the closed form minus the
listed head is a sum over the unlisted zeros (all > g_N =: T):

    A(M) := V_T(M) / (8 M^2) = sum_{g > T} 1/(M^2 + g^2).

Necessary consequences, checked here in ball arithmetic for many M:
 (i)   A(M) > 0 and decreasing in M;
 (ii)  the divided difference  [A(M1)-A(M2)]/(M2^2-M1^2) = sum 1/((M1^2+g^2)(M2^2+g^2))
       must satisfy   0 < dd <= A(M1)/(M2^2+T^2)       (every g > T);
       a zero missing from the list at height g0 << T would violate it;
 (iii) A(M) and dd must agree with the Riemann-von Mangoldt main term
           sum_{g>T} g^-2 ~ (log(qT/2pi)+1)/(2 pi T),
           sum_{g>T} g^-4 ~ (log(qT/2pi)+1/3)/(6 pi T^3)
       to the accuracy expected from |S(T)| = O(1) (a few parts in 10^3);
 (iv)  the same for D_T:  D_T = d/dM V_T  =>  A2 := 2 V_T - M D_T = sum 16M^4/(M^2+g^2)^2 > 0,
       and  A2/(16 M^4) <= A(M)/(M^2+T^2).
The closed forms are the auditor's (Hurwitz-zeta series), the zero balls are
the auditor's.
"""
import math
import sys
from fractions import Fraction

from flint import arb, ctx

import indep_enclose as ie

ctx.prec = 200
bad = []
for q in (4, 3):
    zs = ie.zeros(q)
    T = zs[-1]
    Tf = float(T.mid())
    pts = ["2.5", "3.6", "4.4", "6", "8.5", "12", "17.4", "19.2", "40", "100", "320"]
    A = {}
    for Ms in pts:
        M = ie.qarb(Fraction(Ms))
        s2, D = ie.sigma2_D(q, M)
        v = arb(0)
        d = arb(0)
        for g in zs:
            den = M * M + g * g
            v += 8 * M * M / den
            d += 16 * M * g * g / (den * den)
        VT, DT = s2 - v, D - d
        A[Ms] = VT / (8 * M * M)
        A2 = 2 * VT - M * DT
        ok4 = bool(A2.lower() > 0) and bool((A2 / (16 * M ** 4)).upper() <= (A[Ms] / (M * M + T * T)).upper() * (1 + 1e-9))
        if not (VT.lower() > 0 and DT.lower() > 0 and ok4):
            bad.append((q, Ms, "positivity / (iv)"))
        print(f"q={q} M={Ms:>5}: V_T={VT.str(12)}  A(M)=V_T/(8M^2)={A[Ms].str(12)}  "
              f"A2/(16M^4)={(A2 / (16 * M ** 4)).str(6)}  (iv) {'ok' if ok4 else 'VIOLATED'}")
    rvm2 = (math.log(q * Tf / (2 * math.pi)) + 1) / (2 * math.pi * Tf)
    rvm4 = (math.log(q * Tf / (2 * math.pi)) + 1 / 3) / (6 * math.pi * Tf ** 3)
    print(f"   Riemann-von Mangoldt main terms: sum g^-2 ~ {rvm2:.6e}, sum g^-4 ~ {rvm4:.4e}")
    for M1s, M2s in zip(pts, pts[1:]):
        M1, M2 = ie.qarb(Fraction(M1s)), ie.qarb(Fraction(M2s))
        dd = (A[M1s] - A[M2s]) / (M2 * M2 - M1 * M1)
        cap = A[M1s] / (M2 * M2 + T * T)
        ok = bool(dd.lower() > 0) and bool(dd.upper() <= cap.upper())
        if not ok:
            bad.append((q, M1s, M2s, "(ii)"))
        print(f"   ({M1s:>5},{M2s:>5}): dd = {dd.str(8)}  cap A(M1)/(M2^2+T^2) = {cap.str(8)}  "
              f"dd/RvM4 = {float(dd.mid()) / rvm4:.4f}  {'ok' if ok else 'VIOLATED'}")
    r = float(A["2.5"].mid()) / rvm2
    print(f"   A(2.5)/RvM2 = {r:.5f}")
    if abs(r - 1) > 0.01:
        bad.append((q, "RvM"))
if bad:
    print("TAIL SANITY FAILURES:", bad)
    sys.exit(1)
print("tail sanity checks passed")
