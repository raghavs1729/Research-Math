#!/usr/bin/env python3
"""
nonresonance_certified.py -- certified version of the non-resonance bound
(GRH-ONLY.md, Theorems 6.1 and 6.2), added after Referee 4's item R5.

Two independent rigorous steps:

 (I)  (NR1)-(NR3) in EXACT INTEGER ARITHMETIC.  The stored ordinates (25
      significant digits) are multiplied by 10^24 and treated as integers; every
      linear form occurring in (NR1)-(NR3) has coefficient mass <= 4, so the
      certified enclosure radius 1e-20 of each ordinate moves it by at most
      4e-20.  A form whose exact value at the stored midpoints exceeds 4e-20 in
      absolute value is therefore certified nonzero.  (The lists are complete
      and simple up to T: certify_arb.py --K 0.)

 (II) J(kappa_0, b) of Theorem 6.1 in BALL ARITHMETIC (Arb, 200 bits) at an
      exactly rational feasible point: the float minimiser of nonresonance.py
      with every b_n rounded DOWN to 6 decimals (so b_n >= 0 and sum b_n <= 1
      hold exactly; this is asserted on exact rationals) and kappa_0 rounded
      to 6 decimals.  Ordinates enter as balls mid +/- 1e-20; the set K is
      decided on balls (an undecided zero is left in the tail, which only
      weakens the bound); V_t = V^+ + 1e-10 - V_head - V_K, where V^+ is the
      certified upper bound for V_M read from certified_bounds_results.txt
      (midpoint + radius, rounded UP to 10 decimals).

Output: for each (character, weight, tier) a decimal number d such that the
comparison  1 - J > d  is certified by Arb (outward-rounded lower bound for the
lower logarithmic density of {D_m > 0}).  Tier A uses only the in-repo
certified zero lists; tier B assumes in addition that all zeros up to 10^8/q
are simple (GRH-ONLY.md, Proposition 4.5 and Remark 4.6).

Also included: the grid-quadrature self-test of the closed formulas for
int Q^2 and int kappa^2 used in the proof of Theorem 6.1 (item 1), N = 1,2,3.

Runtime: about 15 s on one core.
"""
import bisect
import itertools
import math
import os
import re
import sys
from decimal import Decimal
from fractions import Fraction

import numpy as np
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from nonresonance import minimise            # float search only; nothing rigorous is imported

ctx.prec = 200
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")
SCALE = 10 ** 24
HEADS = (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20)


# ----------------------------------------------------------------------------
# (I) exact-integer verification of (NR1)-(NR3)
# ----------------------------------------------------------------------------

def load_int(fn):
    out = []
    for l in open(os.path.join(ZDIR, fn)):
        l = l.strip()
        if l:
            v = Decimal(l) * SCALE
            assert v == v.to_integral_value(), "ordinate has more than 24 decimals"
            out.append(int(v))
    return out


def vecs(N, order):
    seen = set()
    for k in range(1, order + 1):
        for idx in itertools.combinations_with_replacement(range(N), k):
            for sg in itertools.product((1, -1), repeat=k):
                v = [0] * N
                for i, s in zip(idx, sg):
                    v[i] += s
                if any(v):
                    seen.add(tuple(v))
    return seen


def nearest(L, s):
    j = bisect.bisect_left(L, s)
    return min(abs(L[jj] - s) for jj in (j - 1, j) if 0 <= jj < len(L))


def exact_nr(fn, log):
    """returns {N: True/False}; True iff all three margins exceed 4e-20 exactly"""
    g = load_int(fn)
    assert all(g[i] < g[i + 1] for i in range(len(g) - 1))
    T = g[-1] + SCALE // 4
    thr = 4 * SCALE // 10 ** 20                      # 4e-20 in units of 1e-24
    ok = {}
    for N in HEADS:
        head, rest = g[:N], g[N:]
        if not 3 * head[-1] + 3 * thr < T:
            ok[N] = False
            continue
        dot = lambda v: abs(sum(a * b for a, b in zip(v, head)))
        V4 = vecs(N, 4)
        m1 = min(dot(v) for v in V4)
        m2 = min(nearest(rest, dot(v)) for v in V4 if sum(map(abs, v)) <= 3)
        D = sorted([rest[i] + rest[j] for i in range(len(rest)) for j in range(i, len(rest))]
                   + [rest[j] - rest[i] for i in range(len(rest)) for j in range(i + 1, len(rest))])
        m3 = min(nearest(D, dot(v)) for v in V4 if sum(map(abs, v)) <= 2)
        ok[N] = min(m1, m2, m3) > thr
        log("    N = %2d: exact margins  NR1 = %.6e   NR2 = %.6e   NR3 = %.6e   (> 4e-20: %s)"
            % (N, m1 / SCALE, m2 / SCALE, m3 / SCALE, ok[N]))
    return ok


# ----------------------------------------------------------------------------
# (II) Arb evaluation of J at a rational feasible point
# ----------------------------------------------------------------------------

def load_balls(fn):
    eps = arb("1e-20")
    return [arb(l.strip()) + arb(0, 1) * eps
            for l in open(os.path.join(ZDIR, fn)) if l.strip()]


def qarb(fr):
    return arb(fr.numerator) / arb(fr.denominator)


def J_arb(k0, b, a, VK, Vt):
    """J(kappa_0, b) of Theorem 6.1; k0, b exact Fractions; a, VK, Vt arb"""
    k0 = qarb(k0)
    b = [qarb(x) for x in b]
    Sab = sum((x * y for x, y in zip(a, b)), arb(0))
    Saa = sum((x * x for x in a), arb(0))
    Sbb = sum((y * y for y in b), arb(0))
    Sa2b2 = sum((x * x * y * y for x, y in zip(a, b)), arb(0))
    Sd = sum(((x - y) ** 2 for x, y in zip(a, b)), arb(0))
    Sb = sum(b, arb(0))
    br = Sd / 2 + Saa * Sbb / 4 + Sab ** 2 / 2 - 3 * Sa2b2 / 8 \
        + (1 + Sbb / 2) * VK + (1 + Sb) ** 2 * Vt
    return (1 - k0) ** 2 + (1 - k0) * k0 * Sab + k0 ** 2 * br


def certified_floor(x, digits=9):
    """largest d (a decimal string with `digits` decimals, found by stepping
    down from the midpoint) such that the Arb comparison x > d is certified"""
    q = 10 ** digits
    d = math.floor(float(x.mid()) * q) + 1
    for _ in range(1000):
        if x > arb(d) / arb(q):
            return ("%%.%df" % digits) % (d / q), d
        d -= 1
    raise RuntimeError("no certified lower bound found")


def read_Vplus():
    """certified upper bounds for V_M from certified_bounds_results.txt, rounded
    UP to 10 decimals: {(char, m_string): (V_A, V_B)} as exact Fractions"""
    txt = open(os.path.join(HERE, "certified_bounds_results.txt")).read()
    out = {}
    for blk in txt.split("character ")[1:]:
        name = blk.split()[0]
        A = {}
        for m, mid, rad in re.findall(r"m = (\S+)\s+sigma\^2 = .*?V_M\^\+    = \[(\S+) \+/- (\S+)\]",
                                      blk, flags=re.S):
            up = Decimal(mid) + Decimal(rad)
            A[m] = Fraction(int((up * 10 ** 10).to_integral_value(rounding="ROUND_CEILING")), 10 ** 10)
        B = {}
        for m, v in re.findall(r"m = (\S+)\s+V_M <= (\S+)", blk):
            # printed midpoint to 12 significant digits of a ball of radius < 1e-12:
            # add 2e-10 before rounding up to 10 decimals
            up = Decimal(v) + Decimal("2e-10")
            B[m] = Fraction(int((up * 10 ** 10).to_integral_value(rounding="ROUND_CEILING")), 10 ** 10)
        for m in A:
            out[(name, m)] = (A[m], B[m])
    return out


def certify(name, zfile, Mstr, mkey, Vp, nr_ok, log):
    g = load_balls(zfile)
    T = g[-1] + arb("0.25")
    M = arb(Mstr)
    amp = [4 * M / (M * M + x * x).sqrt() for x in g]
    var = [x * x / 2 for x in amp]
    ampf = [float(x.mid()) for x in amp]
    log("  %s, M = %s (m = %s)" % (name, Mstr, mkey))
    res = {}
    for tier, V in (("A", Vp[0]), ("B", Vp[1])):
        best = None
        for N in HEADS:
            if not nr_ok.get(N, False):
                continue
            Tp = T - 2 * g[N - 1]
            if not (3 * g[N - 1] < T):
                continue
            inK = [i for i in range(N, len(g)) if g[i] < Tp]         # certified comparisons
            Vh = sum(var[:N], arb(0))
            VK = sum((var[i] for i in inK), arb(0))
            Vt = qarb(V) + arb("1e-10") - Vh - VK
            if not (Vt > 0):
                raise RuntimeError("V_t not certified positive")
            # float minimiser -> exact rational feasible point
            _, k0f, bf = minimise(ampf[:N], float(VK.mid()), float(Vt.mid()))
            b = [Fraction(math.floor(y * 10 ** 6), 10 ** 6) for y in bf]
            k0 = Fraction(round(k0f * 10 ** 6), 10 ** 6)
            assert k0 > 0 and all(y >= 0 for y in b) and sum(b) <= 1          # exact
            J = J_arb(k0, b, amp[:N], VK, Vt)
            dens = 1 - J
            s, d = certified_floor(dens)
            if best is None or d > best[1]:
                best = (s, d, N, J, k0, b)
        s, d, N, J, k0, b = best
        log("    tier %s: N = %2d, J in %s" % (tier, N, J.str(18, radius=True)))
        log("            certified:  lower log-density{D_m > 0}  >=  1 - J  >  %s" % s)
        log("            feasible point: kappa_0 = %s, b = [%s]"
            % (k0, ", ".join("%.6f" % float(y) for y in b)))
        res[tier] = (s, N)
    return res


def grid_selftest(log):
    """closed formulas of Theorem 6.1, item 1, against exact trapezoidal
    quadrature on a 16^N grid (exact for trigonometric degree < 16)"""
    rng = np.random.default_rng(1)
    worst = 0.0
    for N in (1, 2, 3):
        a = rng.uniform(0.05, 0.4, N)
        b = rng.uniform(0, 0.3, N)
        k0 = 0.93
        th = 2 * np.pi * np.arange(16) / 16
        G = np.meshgrid(*([th] * N), indexing="ij")
        H = sum(a[i] * np.cos(G[i]) for i in range(N))
        L = sum(b[i] * np.cos(G[i]) for i in range(N))
        q2 = np.mean((1 - k0 * (1 - L) * (1 + H)) ** 2)
        k2 = np.mean((k0 * (1 - L)) ** 2)
        Sab, Saa, Sbb = (a * b).sum(), (a * a).sum(), (b * b).sum()
        form = (1 - k0) ** 2 + (1 - k0) * k0 * Sab + k0 ** 2 * (
            ((a - b) ** 2).sum() / 2 + Saa * Sbb / 4 + Sab ** 2 / 2 - 3 * (a * a * b * b).sum() / 8)
        worst = max(worst, abs(q2 - form), abs(k2 - k0 ** 2 * (1 + Sbb / 2)))
    log("grid self-test of the closed formulas for int Q^2 and int kappa^2 (N = 1,2,3): "
        "max deviation %.2e" % worst)
    assert worst < 1e-14


def main():
    out = open(os.path.join(HERE, "nonresonance_certified_results.txt"), "w")

    def log(s=""):
        print(s)
        out.write(s + "\n")

    grid_selftest(log)
    Vp = read_Vplus()
    jobs = {
        "chi_4": ("chi4_zeros.txt", [("0.5", "0"), ("0.1", "-0.4"), ("0.25", "-0.25"),
                                     ("0.75", "0.25"), ("1", "0.5")]),
        "chi_3": ("chi3_zeros.txt", [("0.5", "0"), ("0.25", "-0.25"), ("1", "0.5")]),
    }
    summary = []
    for name, (zfile, weights) in jobs.items():
        log("=" * 78)
        log("%s: (NR1)-(NR3) in exact integer arithmetic (ordinates x 10^24)" % name)
        nr_ok = exact_nr(zfile, log)
        log("%s: Arb evaluation of J" % name)
        for Mstr, mkey in weights:
            VA, VB = Vp[(name, mkey)]
            log("  certified variance bounds used: V^+_A = %s, V^+_B = %s"
                % (float(VA), float(VB)))
            r = certify(name, zfile, Mstr, mkey, (VA, VB), nr_ok, log)
            summary.append((name, mkey, r["A"], r["B"]))
    log("=" * 78)
    log("SUMMARY (certified lower bounds for the lower logarithmic density of {D_m > 0}, GRH only)")
    for name, m, A, B in summary:
        log("  %s  m = %-6s  tier A: > %s (N = %d)    tier B: > %s (N = %d)"
            % (name, m, A[0], A[1], B[0], B[1]))
    out.close()


if __name__ == "__main__":
    main()
