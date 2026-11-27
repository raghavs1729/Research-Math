#!/usr/bin/env python3
"""
crosscheck.py -- NON-RIGOROUS independent evaluations of delta_chi(M), to be
compared with the certified enclosures printed by density.py.

Usage:  nice -n 10 python3 crosscheck.py q M1 M2 ... [--mc N] [--out file]

Three independent computations (none of them uses Arb ball arithmetic for the
quadrature, none uses the Gauss-Legendre rule of density.py):

 (A) mpmath: sigma^2 from mpmath's own Dirichlet L-function (mp.dirichlet with
     derivative), the first ordinates re-located by mp.findroot on the
     mpmath Hardy function, and the Gil-Pelaez integral by tanh-sinh
     quadrature (mp.quad) with a Gaussian surrogate exp(-V t^2/2) for the
     unlisted zeros.
 (B) numpy/scipy: composite Simpson rule on a uniform grid, scipy.special.j0,
     same Gaussian surrogate.
 (C) Monte Carlo: random phases for the listed zeros + an independent
     N(0, V) variable for the unlisted ones, common random numbers for all M.

Also prints a consistency check of the variance identity: sigma^2 (closed
form) minus the listed zeros, against the Riemann-von Mangoldt estimate of
the unlisted tail.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.special import j0 as sj0
import mpmath as mpm
from mpmath import mp

from common import real_character, character_values

HERE = os.path.dirname(os.path.abspath(__file__))


def zero_path(q):
    if q in (3, 4):
        return os.path.join(HERE, "..", "..", "prototype", "explore", "chi%d_zeros.txt" % q)
    return os.path.join(HERE, "zeros_q%d.txt" % q)


def load_strings(q):
    with open(zero_path(q)) as fh:
        return [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("M", nargs="+")
    ap.add_argument("--mc", type=float, default=4e6)
    ap.add_argument("--nroots", type=int, default=12)
    ap.add_argument("--out", default=None)
    ap.add_argument("--skip-mp", action="store_true")
    args = ap.parse_args()
    q = args.q
    out = open(args.out, "a") if args.out else None

    def emit(s):
        print(s, flush=True)
        if out:
            out.write(s + "\n")
            out.flush()

    chi, d = real_character(q, 1)
    vals = character_values(chi, q)                 # exact chi(0..q-1)
    zs = load_strings(q)
    g = np.array([float(z) for z in zs])
    Ms = [Fraction(m) for m in args.M]
    emit("# crosscheck.py q=%d (Kronecker (%d/.)), %d listed zeros up to %.4f  [NON-RIGOROUS]"
         % (q, d, len(g), g[-1]))

    mp.dps = 25

    def Lmp(s, der=0):
        return mp.dirichlet(s, vals, der)

    def theta(t):
        return t / 2 * mp.log(q / mp.pi) + mp.im(mp.loggamma(mp.mpf(3) / 4 + 1j * t / 2))

    def Zmp(t):
        return mp.re(mp.exp(1j * theta(t)) * Lmp(mp.mpf(1) / 2 + 1j * t))

    # (A0) independent re-location of the first ordinates
    t0 = time.time()
    worst = 0
    for k in range(min(args.nroots, len(zs))):
        r = mp.findroot(Zmp, (mp.mpf(zs[k]) - mp.mpf("0.001"), mp.mpf(zs[k]) + mp.mpf("0.001")),
                        solver="secant", tol=1e-40, maxsteps=60)
        worst = max(worst, abs(r - mp.mpf(zs[k])))
    emit("(A0) mpmath re-location of the first %d ordinates: max |difference| = %s  (%.0fs)"
         % (min(args.nroots, len(zs)), mp.nstr(worst, 3), time.time() - t0))
    emit("     mpmath L(1/2,chi) = %s" % mp.nstr(Lmp(mp.mpf(1) / 2), 20))

    def sigma2_mp(M):
        s = M + mp.mpf(1) / 2
        return 4 * M * (mp.log(q / mp.pi) / 2 + mp.digamma((s + 1) / 2) / 2 + Lmp(s, 1) / Lmp(s))

    res = {}
    for fr in Ms:
        M = mp.mpf(fr.numerator) / fr.denominator
        Mf = float(M)
        s2 = sigma2_mp(M)
        c = 4 * Mf / np.sqrt(Mf * Mf + g * g)
        V = float(s2) - float(np.sum(c * c / 2))
        # RvM estimate of the unlisted tail: 8 M^2 int_G^inf (1/2pi) log(q t/2pi) / (M^2+t^2) dt,
        # G = midpoint between the last listed zero and the next expected one (count N(G) = #listed)
        G = g[-1]
        dens = lambda t: mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi)
        Gs = G + 0.5 / float(dens(G))
        Vest = float(8 * M * M * mp.quad(lambda t: dens(t) / (M * M + t * t), [Gs, 10 * Gs, 100 * Gs, mp.inf]))
        emit("M=%s: sigma^2 (mpmath closed form) = %s ; listed zeros give %.10f ; "
             "V_tail = %.8f  (RvM estimate of the tail: %.6f)"
             % (fr, mp.nstr(s2, 18), float(np.sum(c * c / 2)), V, Vest))
        # (B) Simpson
        Tend = min(60.0, 0.5 * math.sqrt(Mf * Mf + G * G) / (4 * Mf) * 2)
        n = int(Tend / 0.0005)
        n += n % 2
        t = np.linspace(0.0, Tend, n + 1)
        P = np.ones_like(t)
        for cj in c:
            P *= sj0(cj * t)
        f = P * np.exp(-V * t * t / 2) * np.sinc(t / np.pi)
        h = t[1] - t[0]
        simpson = h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
        dB = 0.5 + simpson / math.pi
        tailB = np.max(np.abs(f[-2000:]))
        emit("   (B) Simpson (h=%.4g, T=%.1f): delta = %.12f   (|integrand| at the end <= %.1e)"
             % (h, Tend, dB, tailB))
        res[fr] = {"B": dB, "V": V, "c": c}
        # (A) mpmath tanh-sinh
        if not args.skip_mp:
            t0 = time.time()
            mp.dps = 18
            cm = [mp.mpf(float(x)) for x in c]
            Vm = mp.mpf(V)

            def integrand(tt):
                if tt == 0:
                    return mp.mpf(1)
                p = mp.exp(-Vm * tt * tt / 2) * mp.sin(tt) / tt
                for cj in cm:
                    p *= mp.j0(cj * tt)
                    if abs(p) < 1e-30:
                        return mp.mpf(0)
                return p
            TA = 16.0
            pts = [mp.mpf(k) / 2 for k in range(0, int(2 * TA) + 1)]
            IA = mp.quad(integrand, pts, maxdegree=7)
            dA = mp.mpf(1) / 2 + IA / mp.pi
            mp.dps = 25
            emit("   (A) mpmath tanh-sinh on [0,%g]: delta = %s   (%.0fs)" % (TA, mp.nstr(dA, 13), time.time() - t0))
            res[fr]["A"] = float(dA)

    # (C) Monte Carlo with common random numbers
    N = int(args.mc)
    if N > 0:
        t0 = time.time()
        rng = np.random.default_rng(20261124)
        B = 100_000
        nb = N // B
        hits = {fr: [] for fr in Ms}
        for _ in range(nb):
            cs = np.cos(rng.uniform(0.0, 2 * np.pi, size=(B, len(g))))
            zt = rng.standard_normal(B)
            for fr in Ms:
                S = cs @ res[fr]["c"] + math.sqrt(res[fr]["V"]) * zt
                hits[fr].append((S > -1.0).mean())
        emit("   (C) Monte Carlo, %d samples in %d batches (common phases)  (%.0fs)" % (nb * B, nb, time.time() - t0))
        for fr in Ms:
            a = np.array(hits[fr])
            res[fr]["C"] = a.mean()
            emit("       M=%s: delta = %.6f +- %.6f (1 s.e.)" % (fr, a.mean(), a.std(ddof=1) / math.sqrt(nb)))
        for i in range(len(Ms)):
            for j in range(i + 1, len(Ms)):
                dd = np.array(hits[Ms[j]]) - np.array(hits[Ms[i]])
                emit("       delta(%s) - delta(%s) = %+.6f +- %.6f (1 s.e., paired)"
                     % (Ms[j], Ms[i], dd.mean(), dd.std(ddof=1) / math.sqrt(nb)))


if __name__ == "__main__":
    main()
