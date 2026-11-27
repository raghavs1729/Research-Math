#!/usr/bin/env python3
"""
reflection_certified.py -- certified evaluation (Arb, 200 bits) of the
reflection bounds of DENSITY.md, Theorems 5.2 and 5.3 (GRH only, no LI).

Inputs (all from the repository, nothing recomputed here):
  * certified zero lists prototype/explore/chi{4,3}_zeros.txt (complete and
    simple up to T; enclosure radius 1e-20), as balls;
  * certified upper bounds V^+ for the variance V_M (tier A / tier B) and the
    closed-form balls for sigma_M^2, read from
    monotonicity/grhonly/certified_bounds_results.txt;
  * (NR1)-(NR3) in exact integer arithmetic: exact_nr() of
    monotonicity/grhonly/nonresonance_certified.py (imported unchanged).

Theorem 5.2 (no non-resonance needed), for 0 < c < 2 - sigma^2:
    nu_M((-inf,0]) <= (V + (1-c)^2) / (c^2 + (2 - sigma^2 - c)^2).
Theorem 5.3 (head of N zeros satisfying (NR1)-(NR3)):
    2 nu_M((-inf,0]) <= J_R,   see J_R_arb() below.

Output: reflection_certified_results.txt.  Every printed density is an
outward-rounded decimal d such that the Arb comparison 1 - bound > d is
certified.  Runtime: about 20 s.
"""
import math, os, re, sys
from fractions import Fraction
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
GRH = os.path.join(HERE, "..", "grhonly")
sys.path.insert(0, GRH)
import nonresonance_certified as nc          # exact_nr, read_Vplus, load_balls, qarb, certified_floor

ctx.prec = 200
HEADS = nc.HEADS


def read_sigma2():
    """closed-form balls sigma_M^2 from certified_bounds_results.txt: {(char, m): arb}"""
    txt = open(os.path.join(GRH, "certified_bounds_results.txt")).read()
    out = {}
    for blk in txt.split("character ")[1:]:
        name = blk.split()[0]
        for m, mid, rad in re.findall(r"m = (\S+)\s+sigma\^2 = \[(\S+) \+/- (\S+)\]", blk):
            out[(name, m)] = arb(mid) + arb(0, 1) * arb(rad) * 2     # radius doubled (decimal printing)
    return out


def thm52(V, s2):
    """min over a rational grid of c of B(c); returns (arb bound, c)"""
    best = None
    for i in range(900, 1101):
        c = Fraction(i, 1000)
        ca = nc.qarb(c)
        d2 = 2 - s2 - ca
        if not (d2 > 0):
            continue
        B = (V + (1 - ca) ** 2) / (ca * ca + d2 * d2)
        if best is None or float(B.mid()) < float(best[0].mid()):
            best = (B, c)
    return best


def J_R_arb(N, g, M, s2, Vplus, Ve_low, inK):
    """J_R of Theorem 5.3; returns arb upper bound for 2 nu((-inf,0]) or None"""
    w = [8 * M * M / (M * M + x * x) for x in g]
    s = [4 * M * x / (M * M + x * x) for x in g]
    wh = sum(w[:N], arb(0))
    rho = (s2 - wh) / 2                         # = (1/2) sum over zeros outside the head of m w
    rho_up = rho                                # ball enclosing the exact value
    if not (rho_up < 1):
        return None
    xmax = wh / (2 * (1 - rho_up))
    if not (xmax < 1):
        return None
    e = xmax * xmax / (1 - xmax)
    a0 = (1 + e) / (1 - rho_up)
    k = [x / (2 * (1 - rho_up) ** 2) for x in w[:N]]
    Ss = sum((x * x for x in s[:N]), arb(0))
    Sk = sum((x * x for x in k), arb(0))
    Ssk = sum((x * x * y * y for x, y in zip(s[:N], k)), arb(0))
    VKo = sum((s[i] * s[i] for i in inK), arb(0)) / 2
    Vto = Vplus - Ve_low - Ss / 2 - VKo
    if not (Vto > 0):
        raise RuntimeError("V_t^o not certified positive")
    return (a0 * a0 * Ss / 2 + Ss * Sk / 4 - Ssk / 8
            + (a0 * a0 + Sk / 2) * VKo + (a0 + sum(k, arb(0))) ** 2 * Vto)


def main():
    out = open(os.path.join(HERE, "reflection_certified_results.txt"), "w")

    def log(t=""):
        print(t)
        out.write(t + "\n")

    Vp = nc.read_Vplus()
    S2 = read_sigma2()
    jobs = {
        "chi_4": ("chi4_zeros.txt", [("0.5", "0"), ("0.1", "-0.4"), ("0.25", "-0.25"),
                                     ("0.75", "0.25"), ("1", "0.5")]),
        "chi_3": ("chi3_zeros.txt", [("0.5", "0"), ("0.25", "-0.25"), ("1", "0.5")]),
    }
    summary = []
    for name, (zfile, weights) in jobs.items():
        log("=" * 78)
        log("%s: (NR1)-(NR3) in exact integer arithmetic (imported check)" % name)
        nr_ok = nc.exact_nr(zfile, lambda *a: None)
        log("  head sizes passing: %s" % [N for N in HEADS if nr_ok.get(N, False)])
        g = nc.load_balls(zfile)
        T = g[-1] + arb("0.25")
        for Mstr, mkey in weights:
            M = arb(Mstr)
            s2 = S2[(name, mkey)]
            if not (s2 < 2):
                raise RuntimeError("sigma^2 < 2 not certified")
            w = [8 * M * M / (M * M + x * x) for x in g]
            Ve_low = sum((x * x for x in w), arb(0)) / 8       # V_e >= this (listed zeros are simple)
            log("  %s, M = %s (m = %s): sigma_M^2 in %s" % (name, Mstr, mkey, s2.str(15, radius=True)))
            row = [name, mkey]
            for tier, Vfr in (("A", Vp[(name, mkey)][0]), ("B", Vp[(name, mkey)][1])):
                V = nc.qarb(Vfr) + arb("1e-10")
                B, c = thm52(V, s2)
                d52, _ = nc.certified_floor(1 - B)
                log("    tier %s  Thm 5.2: nu((-inf,0]) <= %s  (c = %s);  density > %s"
                    % (tier, B.str(12, radius=True), c, d52))
                best = None
                for N in HEADS:
                    if not nr_ok.get(N, False) or not (3 * g[N - 1] < T):
                        continue
                    Tp = T - 2 * g[N - 1]
                    inK = [i for i in range(N, len(g)) if g[i] < Tp]
                    JR = J_R_arb(N, g, M, s2, V, Ve_low, inK)
                    if JR is None:
                        continue
                    dens = 1 - JR / 2
                    sdec, d = nc.certified_floor(dens)
                    if best is None or d > best[1]:
                        best = (sdec, d, N, JR / 2)
                if best:
                    log("    tier %s  Thm 5.3: N = %2d, nu((-inf,0]) <= %s;  density > %s"
                        % (tier, best[2], best[3].str(12, radius=True), best[0]))
                row.append((d52, best[0] if best else "-", best[2] if best else 0))
            summary.append(row)
    log("=" * 78)
    log("SUMMARY: certified lower bounds for the lower log-density of {D_m > 0} (GRH only)")
    log("  char   m       tier A: Thm 5.2 / Thm 5.3 (N)          tier B: Thm 5.2 / Thm 5.3 (N)")
    for name, m, A, Bt in summary:
        log("  %s  %-6s  %s / %s (%d)      %s / %s (%d)" % (name, m, A[0], A[1], A[2], Bt[0], Bt[1], Bt[2]))
    out.close()


if __name__ == "__main__":
    main()
