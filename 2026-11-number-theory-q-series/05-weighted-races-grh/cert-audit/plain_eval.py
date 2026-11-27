#!/usr/bin/env python3
"""AUDITOR: plain (NON-rigorous) high-precision evaluation of delta'(M), mpmath only.

No python-flint, no code from the audited directory.  30 significant digits.

  delta'(M) = (1/pi) int_0^inf (sin t/t) d/dM prod_g J0(c_g(M) t) dt
            ~ (1/pi) int_0^T  sinc(t) e^{-V_T t^2/2} ( dP(t) - P(t) (t^2/2) D_T ) dt

with P the product over ALL listed zeros, the unlisted zeros replaced by a
Gaussian of the exact variance V_T = sigma^2 - sum_listed c^2/2 (and
D_T = D - sum_listed c c'), T = 18/sigma, composite 24-point Gauss-Legendre on
panels of width 1/sigma.  The closed forms come from mpmath Hurwitz-zeta
s-derivatives.  Model error (Gaussian tail) ~ 1e-6 relative, quadrature
error far below that (checked by `--check`, which doubles T and the panel count).

Also evaluates delta(M) itself (same model), for the finite-difference check.

Usage: plain_eval.py --q 4 --M 3.6 5 8.5 [--check]
       plain_eval.py --q 4 --npoints 20 --lo 3.6 --hi 17.4 --out plain_q4.txt
"""
import argparse
import os
import sys
import time

from mpmath import mp, mpf

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))


def load(q):
    with open(os.path.join(ZDIR, f"chi{q}_zeros.txt")) as fh:
        return [mpf(x.strip()) for x in fh if x.strip()]


def closed(q, M):
    """sigma^2(M), D(M).  L(s,chi) = q^{-s} H(s), H = sum_a chi(a) zeta(s, a/q);
    Hurwitz-zeta s-derivatives from mpmath.zeta(s, a, n).  (mpmath.dirichlet(s,
    chi, derivative=2) is NOT used: in mpmath 1.4.1 it returns a wrong value;
    see AUDIT.md, finding M-9.)"""
    s = M + mpf(1) / 2
    if q == 4:
        ab = [(mpf(1) / 4, 1), (mpf(3) / 4, -1)]
    else:
        ab = [(mpf(1) / 3, 1), (mpf(2) / 3, -1)]
    H0, H1, H2 = (sum(e * mp.zeta(s, a, n) for a, e in ab) for n in (0, 1, 2))
    lp = -mp.log(q) + H1 / H0
    lpp = H2 / H0 - (H1 / H0) ** 2
    d1 = mp.log(mpf(q) / mp.pi) / 2 + mp.psi(0, (s + 1) / 2) / 2 + lp
    d2 = mp.psi(1, (s + 1) / 2) / 4 + lpp
    return 4 * M * d1, 4 * d1 + 4 * M * d2


_GL = {}


def gl(n):
    if n not in _GL:
        xs = []
        ws = []
        # roots of Legendre P_n by Newton
        for k in range(1, n + 1):
            x = mp.cos(mp.pi * (k - mpf(1) / 4) / (n + mpf(1) / 2))
            for _ in range(100):
                p0, p1 = mpf(1), x
                for j in range(2, n + 1):
                    p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
                dp = n * (x * p1 - p0) / (x * x - 1)
                dx = p1 / dp
                x -= dx
                if abs(dx) < mpf(10) ** (-mp.dps + 3):
                    break
            xs.append(x)
            ws.append(2 / ((1 - x * x) * dp * dp))
        _GL[n] = (xs, ws)
    return _GL[n]


def evaluate(q, M, zs=None, scale=1, npt=24):
    """Returns (delta, delta', sigma^2, D)."""
    M = mpf(M)
    if zs is None:
        zs = load(q)
    s2, D = closed(q, M)
    cs = [4 * M / mp.sqrt(M * M + g * g) for g in zs]
    cps = [4 * g * g / (M * M + g * g) ** mpf("1.5") for g in zs]
    VT = s2 - sum(c * c for c in cs) / 2
    DT = D - sum(c * cp for c, cp in zip(cs, cps))
    sig = mp.sqrt(s2)
    T = 18 * scale / sig
    npan = 18 * scale * scale
    h = T / npan
    xs, ws = gl(npt)
    acc0 = mpf(0)
    acc1 = mpf(0)
    n = len(zs)
    for k in range(npan):
        a = k * h
        for x, w in zip(xs, ws):
            t = a + h * (x + 1) / 2
            J0 = [mp.besselj(0, c * t) for c in cs]
            J1 = [mp.besselj(1, c * t) for c in cs]
            pre = [mpf(1)] * (n + 1)
            for i in range(n):
                pre[i + 1] = pre[i] * J0[i]
            P = pre[n]
            suf = mpf(1)
            a1 = mpf(0)
            for i in range(n - 1, -1, -1):
                a1 += cps[i] * J1[i] * pre[i] * suf
                suf *= J0[i]
            dP = -t * a1
            G = mp.exp(-VT * t * t / 2)
            sc = mp.sin(t) / t
            acc0 += w * h / 2 * sc * G * P
            acc1 += w * h / 2 * sc * G * (dP - P * t * t / 2 * DT)
    return mpf(1) / 2 + acc0 / mp.pi, acc1 / mp.pi, s2, D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--M", nargs="*", default=[])
    ap.add_argument("--npoints", type=int, default=0)
    ap.add_argument("--lo")
    ap.add_argument("--hi")
    ap.add_argument("--out")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--fd", nargs="*", default=[],
                    help="points at which the formula for delta' is compared with a central difference of delta (h = 1e-3, plus Richardson with 2h)")
    a = ap.parse_args()
    if a.fd:
        zs = load(a.q)
        out = open(a.out, "a") if a.out else sys.stdout
        h = mpf("0.001")
        for Ms in a.fd:
            M = mpf(Ms)
            _, dd, _, _ = evaluate(a.q, M, zs)
            d1p, _, _, _ = evaluate(a.q, M + h, zs)
            d1m, _, _, _ = evaluate(a.q, M - h, zs)
            d2p, _, _, _ = evaluate(a.q, M + 2 * h, zs)
            d2m, _, _, _ = evaluate(a.q, M - 2 * h, zs)
            fd1 = (d1p - d1m) / (2 * h)
            fd2 = (d2p - d2m) / (4 * h)
            rich = (4 * fd1 - fd2) / 3
            print(f"q={a.q} M={Ms} formula={mp.nstr(dd, 16)} central_fd(h=1e-3)={mp.nstr(fd1, 16)} "
                  f"richardson={mp.nstr(rich, 16)} |formula-richardson|={mp.nstr(abs(dd - rich), 3)}", file=out, flush=True)
        return
    pts = [mpf(x) for x in a.M]
    if a.npoints:
        lo, hi = mpf(a.lo), mpf(a.hi)
        # endpoints included; interior points deliberately irrational-looking offsets
        pts += [lo + (hi - lo) * k / (a.npoints - 1) for k in range(a.npoints)]
    zs = load(a.q)
    out = open(a.out, "a") if a.out else sys.stdout
    for M in pts:
        t0 = time.time()
        d, dd, s2, D = evaluate(a.q, M, zs)
        line = f"q={a.q} M={mp.nstr(M, 12)} delta={mp.nstr(d, 16)} ddelta={mp.nstr(dd, 16)} sigma2={mp.nstr(s2, 14)} D={mp.nstr(D, 14)}"
        if a.check:
            d2, dd2, _, _ = evaluate(a.q, M, zs, scale=2)
            line += f" |ddelta-ddelta(2T,4x panels)|={mp.nstr(abs(dd - dd2), 3)}"
        line += f" secs={time.time() - t0:.0f}"
        print(line, file=out, flush=True)


if __name__ == "__main__":
    main()
