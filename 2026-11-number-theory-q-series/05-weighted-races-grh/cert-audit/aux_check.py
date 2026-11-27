#!/usr/bin/env python3
"""AUDITOR: independent re-certification of aux_results.txt.

Uses the auditor's own closed forms (Hurwitz-zeta power series, with the
deflated zeta so that s = 1 is allowed) and exact-Fraction zero balls; nothing
is imported from the audited directory.  For every line of the audited
aux_results.txt the corresponding quantity is re-enclosed, the two enclosures
must intersect, and the sign / threshold statements are re-certified.

Inequalities (re-derived, AUDIT.md section 1.8-1.10):
  a_g = 4M^2/(M^2+g^2);  Q-term = a_g^{3/2}/2;
  A1 = sum_U a = V_T/2,  A2 = sum_U a^2 = 2 V_T - M D_T   (U = unlisted zeros)
  A2/(2 sqrt(amax)) <= Q_tail <= min( sqrt(amax) A1, sqrt(A1 A2) )/2,  amax = 4M^2/(M^2+g_N^2)
  W_M = sum_U 4/(M^2+g^2) = V_T(M)/(2M^2),   W_{M'} <= Lambda_tail <= sqrt(W_M W_{M'})  (M < M').
"""
import os
import re
import sys
from fractions import Fraction

from flint import arb, arb_series, ctx

import indep_enclose as ie

HERE = os.path.dirname(os.path.abspath(__file__))
AUX = os.path.join(HERE, "..", "cert", "aux_results.txt")
ctx.prec = 200
BAD = []


def check(msg, cond):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        BAD.append(msg)


def d1_deflated(q, s):
    """(log Lambda)'(s), (log Lambda)''(s) for real s >= 1 (deflated Hurwitz
    zeta: the poles of zeta(s,a1) and zeta(s,a2) cancel in the difference)."""
    x = arb_series([s, 1], prec=3)
    a1, a2 = (arb(1) / 4, arb(3) / 4) if q == 4 else (arb(1) / 3, arb(2) / 3)
    hz = x.zeta(a1, deflate=True) - x.zeta(a2, deflate=True)
    lg = hz.log()
    half = (x + 1) / 2
    L = half * (arb(q) / arb.pi()).log() + half.lgamma() + lg - x * arb(q).log()
    return L[1], 2 * L[2]


def sig2_D(q, M):
    d1, d2 = d1_deflated(q, M + arb(1) / 2)
    return 4 * M * d1, 4 * d1 + 4 * M * d2


def tail(q, M):
    zs = ie.zeros(q)
    s2, D = sig2_D(q, M)
    v = arb(0)
    d = arb(0)
    for g in zs:
        den = M * M + g * g
        v += 8 * M * M / den
        d += 16 * M * g * g / (den * den)
    return s2 - v, D - d, s2


def Q(q, Mf):
    M = ie.qarb(Mf)
    zs = ie.zeros(q)
    VT, DT, _ = tail(q, M)
    assert VT.lower() > 0 and DT.lower() > 0
    head = arb(0)
    for g in zs:
        den = M * M + g * g
        head += 4 * M ** 3 / (den * den.sqrt())
    A1 = VT / 2
    A2 = 2 * VT - M * DT
    assert A2.lower() > 0
    sa = (4 * M * M / (M * M + zs[-1] * zs[-1])).sqrt()
    hi1 = (sa * A1 / 2).upper()
    hi2 = ((A1 * A2).sqrt() / 2).upper()
    hi = hi1 if hi1 < hi2 else hi2
    lo = (A2 / (2 * sa)).lower()
    simple_hi = ((M / zs[-1]) * VT / 2).upper()        # the cruder bound (M/g_N) V_T/2
    return head + arb(lo).union(arb(hi)), head + arb(0).union(arb(simple_hi))


def Lam(q, Mf, Mpf):
    M, Mp = ie.qarb(Mf), ie.qarb(Mpf)
    zs = ie.zeros(q)
    head = arb(0)
    for g in zs:
        head += 4 / ((M * M + g * g) * (Mp * Mp + g * g)).sqrt()
    VT, _, s2 = tail(q, M)
    VTp, _, s2p = tail(q, Mp)
    assert VT.lower() > 0 and VTp.lower() > 0
    W, Wp = VT / (2 * M * M), VTp / (2 * Mp * Mp)
    t = arb(Wp.lower()).union(arb((W * Wp).sqrt().upper()))
    return 1 / (M * Mp) - head - t, (s2 * s2p).sqrt()


def iv(lo, hi):
    return arb(lo).union(arb(hi))


txt = open(AUX).read()
print("== deflated closed form agrees with the non-deflated one at s = 4.1 ==")
for q in (3, 4):
    a = d1_deflated(q, arb("4.1"))
    b = ie.loglam_d12(q, arb("4.1"))
    check(f"q={q}", a[0].overlaps(b[0]) and a[1].overlaps(b[1]) and float(a[0].rad()) < 1e-40)

print("== (a) Q ==")
for q in (4, 3):
    m = re.search(rf"q={q}: Q\(([\d.]+)\) in \[([\d.]+), ([\d.]+)\]", txt)
    Mreq, lo, hi = m.group(1), m.group(2), m.group(3)
    mine, crude = Q(q, Fraction(Mreq))
    print(f"   q={q} Q({Mreq}): audited [{lo}, {hi}]  auditor {mine.str(12)}  crude-tail {crude.str(8)}")
    check(f"q={q}: Q({Mreq}) < 1 re-certified (two-moment tail)", bool(mine.upper() < 1))
    check(f"q={q}: Q({Mreq}) < 1 re-certified with the crude tail (M/g_N) V_T/2 only", bool(crude.upper() < 1))
    check(f"q={q}: enclosures intersect", mine.overlaps(iv(lo, hi)))
    m = re.search(rf"q={q}: hence \(Q strictly increasing\)\s+([\d.]+) < M_0\(chi_{q}\) < ([\d.]+)", txt)
    below, above = m.group(1), m.group(2)
    qb, _ = Q(q, Fraction(below))
    qa, _ = Q(q, Fraction(above))
    check(f"q={q}: Q({below}) < 1 and Q({above}) > 1 re-certified  => {below} < M_0 < {above}",
          bool(qb.upper() < 1) and bool(qa.lower() > 1))

print("== (b1) sigma^2 = 2 ==")
for q in (4, 3):
    m = re.search(rf"q={q}: hence\s+([\d.]+) <= M\* < ([\d.]+)", txt)
    below, above = m.group(1), m.group(2)
    sb, _ = sig2_D(q, ie.qarb(Fraction(below)))
    sa, _ = sig2_D(q, ie.qarb(Fraction(above)))
    print(f"   q={q}: sigma^2({below}) = {sb.str(25)}, sigma^2({above}) = {sa.str(25)}")
    check(f"q={q}: sigma^2({below}) <= 2 < sigma^2({above}) re-certified", bool(sb.upper() <= 2) and bool(sa.lower() > 2))
    for Ms, s in ((below, sb), (above, sa)):
        m2 = re.search(rf"q={q}: sigma\^2\({re.escape(Ms)}\) in \[([\d.]+), ([\d.]+)\]", txt)
        check(f"q={q}: audited sigma^2({Ms}) interval contains the auditor's ball", iv(m2.group(1), m2.group(2)).contains(s))

print("== (b2) Lambda(M, M') ==")
for m in re.finditer(r"q=(\d) m=(\d) m'=(\d)\s+\[([-\d.]+), ([-\d.]+)\]\s+sign (.)\s+sigma\*sigma' <= ([\d.]+)", txt):
    q, a, b = int(m.group(1)), int(m.group(2)), int(m.group(3))
    lam, ss = Lam(q, Fraction(2 * a + 1, 2), Fraction(2 * b + 1, 2))
    sgn = "+" if lam.lower() > 0 else ("-" if lam.upper() < 0 else "?")
    check(f"q={q} m={a} m'={b}: auditor {lam.str(14)} sign {sgn}; audited [{m.group(4)}, {m.group(5)}] sign {m.group(6)}",
          lam.overlaps(iv(m.group(4), m.group(5))) and sgn == m.group(6) and bool(ss.upper() <= arb(m.group(7)) + arb("1e-9")))

print("== sigma_m^2, m = 0..3 ==")
for m in re.finditer(r"q=(\d) m=(\d): sigma\^2 in \[([\d.]+), ([\d.]+)\]", txt):
    q, a = int(m.group(1)), int(m.group(2))
    s, _ = sig2_D(q, ie.qarb(Fraction(2 * a + 1, 2)))
    check(f"q={q} m={a}: audited [{m.group(3)}, {m.group(4)}] contains auditor's {s.str(20)}", iv(m.group(3), m.group(4)).contains(s))

print()
if BAD:
    print("AUX CHECK FAILURES:", BAD)
    sys.exit(1)
print("aux_results.txt: every statement independently re-certified")
