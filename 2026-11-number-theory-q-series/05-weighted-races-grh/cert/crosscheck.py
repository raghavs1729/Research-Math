#!/usr/bin/env python3
"""Independent, NON-rigorous cross-checks of the delta' certificate.

Nothing here is part of the proof.  It exists to catch a wrong formula or a
wrong normalisation in certify_derivative.py, using code paths that share
nothing with it:

 A. float64 evaluation (numpy + scipy.special) of
      delta(M)  = 1/2 + (1/pi) int (sin t/t) Phi_M(t) dt,
      delta'(M) = (1/pi) int (sin t/t) dPhi_M/dM dt,
    with the unseen zeros replaced by a Gaussian of the exact tail variance;
    the closed forms sigma^2, D come from mpmath (numerical differentiation of
    log Lambda built on mpmath.dirichlet), not from the Arb Dirichlet series.
 B. the float delta(M) is compared with the 30-point CERTIFIED density grid
    (prototype/explore/race_mono_grid_certified_results.txt);
 C. delta'(M) (formula) is compared with a central finite difference of the
    float delta, and both must lie inside the certified cell enclosure;
 D. integral consistency: for consecutive certified grid points m1 < m2 inside
    the certified range, sum_cells [lo,hi]*|cell cap [M1,M2]| must overlap the
    certified interval for delta(m2) - delta(m1).  This uses ONLY certified
    numbers on both sides.

Usage: crosscheck.py [--q 4] [--results results_q4.txt]
"""

from __future__ import annotations

import argparse
import os
import re
from fractions import Fraction

import numpy as np
from mpmath import mp, mpf
from scipy.special import j0, j1

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))


def zeros(q):
    with open(os.path.join(EXP, f"chi{q}_zeros.txt")) as fh:
        return np.array([float(x) for x in fh if x.strip()])


def closed_forms(q, M):
    """sigma^2(M), D(M) from mpmath (independent of the Arb code)."""
    mp.dps = 30
    chi = [0, 1, 0, -1] if q == 4 else [0, 1, -1]

    def logLam(x):
        return (mp.log(mpf(q) / mp.pi) * (x + 1) / 2 + mp.loggamma((x + 1) / 2)
                + mp.log(mp.dirichlet(x, chi)))

    s = mpf(M) + mpf(1) / 2
    d1 = mp.diff(logLam, s, 1)
    d2 = mp.diff(logLam, s, 2)
    return float(4 * mpf(M) * d1), float(4 * d1 + 4 * mpf(M) * d2)


def delta_and_derivative(q, M, g, tmax=None, n=None):
    sig2, D = closed_forms(q, M)
    c = 4 * M / np.sqrt(M * M + g * g)
    cp = 4 * g * g / (M * M + g * g) ** 1.5
    VT = sig2 - np.sum(c * c) / 2
    DT = D - np.sum(c * cp)
    if tmax is None:
        tmax = min(40.0, max(6.0, 60.0 / np.sqrt(sig2)))
    if n is None:
        n = int(tmax * 4000) // 2 * 2
    t = np.linspace(0.0, tmax, n + 1)
    w = np.ones(n + 1)
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    w *= (tmax / n) / 3                                  # Simpson
    P = np.ones_like(t)
    S = np.zeros_like(t)                                 # sum c' t J1/J0 handled safely below
    X = np.outer(t, c)                                   # (n+1) x N
    J0 = j0(X)
    J1 = j1(X)
    P = np.prod(J0, axis=1)
    # dP/dM = -sum_i c'_i t J1_i prod_{j != i} J0_j  via prefix/suffix products
    pre = np.cumprod(J0, axis=1)
    suf = np.cumprod(J0[:, ::-1], axis=1)[:, ::-1]
    pre = np.concatenate([np.ones((len(t), 1)), pre[:, :-1]], axis=1)
    suf = np.concatenate([suf[:, 1:], np.ones((len(t), 1))], axis=1)
    dP = -(t[:, None] * cp[None, :] * J1 * pre * suf).sum(axis=1)
    G = np.exp(-VT * t * t / 2)
    sinc = np.sinc(t / np.pi)
    Phi = P * G
    dPhi = G * (dP - P * (t * t / 2) * DT)
    delta = 0.5 + np.sum(w * sinc * Phi) / np.pi
    ddelta = np.sum(w * sinc * dPhi) / np.pi
    return delta, ddelta


def read_cells(path):
    cells = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if p[-1] != "PASS":
                continue
            cells.append((Fraction(p[0]), Fraction(p[1]),
                          float(p[3].split("=")[1]), float(p[2].split("=")[1])))
    cells.sort()
    return cells


def read_grid(q):
    path = os.path.join(EXP, "race_mono_grid_certified_results.txt")
    out = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            mm = re.match(r"== modulus (\d)", line)
            if mm:
                cur = int(mm.group(1))
                continue
            mm = re.match(r"\s*m=\s*([-\d.]+)\s+delta in \[([\d.]+), ([\d.]+)\]", line)
            if mm and cur == q:
                out[Fraction(mm.group(1))] = (float(mm.group(2)), float(mm.group(3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=4)
    ap.add_argument("--results", nargs="+")
    ap.add_argument("--points", default="")
    args = ap.parse_args()
    q = args.q
    res = args.results
    if not res:
        res = [os.path.join(HERE, f"results_q{q}.txt")]
        extra = os.path.join(HERE, f"results_q{q}_stretch.txt")
        if os.path.exists(extra):
            res.append(extra)
    g = zeros(q)
    cells = sorted(c for path in res for c in read_cells(path))
    lo_all, hi_all = cells[0][0], cells[-1][1]
    print(f"q={q}: {len(cells)} certified cells in {[os.path.basename(p) for p in res]}, "
          f"spanning [{float(lo_all)}, {float(hi_all)}]")
    grid = read_grid(q)

    print("\nB. float delta(M) vs certified density grid")
    worst = 0.0
    for m in sorted(grid):
        M = float(m) + 0.5
        if M < 2.0 or M > 120:
            continue
        d, _ = delta_and_derivative(q, M, g)
        lo, hi = grid[m]
        dev = max(lo - d, d - hi, 0.0)
        worst = max(worst, dev)
        print(f"   m={float(m):6.2f}  float delta={d:.10f}  certified [{lo:.10f}, {hi:.10f}]  "
              f"{'inside' if dev == 0 else 'outside by %.1e' % dev}")
    print(f"   worst excursion outside a certified density interval: {worst:.1e}")

    print("\nC. float delta'(M): formula vs central difference vs certified cell")
    if args.points:
        pts = [float(x) for x in args.points.split(",")]
    else:
        a, b = float(lo_all), float(hi_all)
        base = {4: [3.7, 5.0, 8.5, 12.0, 17.3, 17.397, 40.0, 100.5, 200.0],
                3: [4.5, 6.0, 10.0, 15.0, 19.1, 19.197, 40.0, 100.5, 200.0]}[q]
        pts = sorted(set([a + 0.003] + [x for x in base if a + 0.003 < x < b - 0.003] + [b - 0.003]))
    allin = True
    for M in pts:
        _, dd = delta_and_derivative(q, M, g)
        h = 1e-3
        dp, _ = delta_and_derivative(q, M + h, g)
        dm, _ = delta_and_derivative(q, M - h, g)
        fd = (dp - dm) / (2 * h)
        cell = [c for c in cells if c[0] <= Fraction(repr(M)) <= c[1]]
        if not cell:
            print(f"   M={M:.4f}: no certified cell")
            allin = False
            continue
        a, b, lo, hi = cell[0]
        ok = lo <= dd <= hi and lo <= fd <= hi
        allin &= ok
        print(f"   M={M:9.4f}  formula={dd:.8e}  FD={fd:.8e}  |diff|={abs(dd - fd):.1e}  "
              f"cell [{float(a)}, {float(b)}] enclosure [{lo:.6e}, {hi:.6e}]  "
              f"{'INSIDE' if ok else 'OUTSIDE'}")
    print(f"   all float values inside their certified cell enclosures: {allin}")

    print("\nD. integral of the certified delta' enclosures vs certified density differences")
    ms = sorted(m for m in grid if lo_all <= m + Fraction(1, 2) <= hi_all)
    allok = True
    for m1, m2 in zip(ms, ms[1:]):
        M1, M2 = m1 + Fraction(1, 2), m2 + Fraction(1, 2)
        ilo = ihi = 0.0
        cov = Fraction(0)
        for a, b, lo, hi in cells:
            ov = min(b, M2) - max(a, M1)
            if ov > 0:
                ilo += lo * float(ov)
                ihi += hi * float(ov)
                cov += ov
        if cov != M2 - M1:
            print(f"   m in [{float(m1)}, {float(m2)}]: not fully covered")
            allok = False
            continue
        dlo = grid[m2][0] - grid[m1][1]
        dhi = grid[m2][1] - grid[m1][0]
        ok = not (ihi < dlo or ilo > dhi)
        allok &= ok
        print(f"   m in [{float(m1):5.1f}, {float(m2):5.1f}]: int delta' in [{ilo:.8f}, {ihi:.8f}]   "
              f"delta(m2)-delta(m1) in [{dlo:.8f}, {dhi:.8f}]  "
              f"{'CONSISTENT' if ok else 'INCONSISTENT'}  (rel. width {abs(ihi - ilo) / abs(dlo):.2%})")
    print(f"   all consistent: {allok}")


if __name__ == "__main__":
    main()
