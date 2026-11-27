#!/usr/bin/env python3
"""
lzeros_cert.py -- certified zeros of L(s, chi) for ALL nonprincipal characters
chi modulo a prime q (all primitive), including the complex ones, up to a
height T, with certified completeness.  Adapted from
monotonicity/nonmonotone/zeros_cert.py (odd real characters only).

Usage:  nice -n 10 python3 lzeros_cert.py q H [--pairs i,j,...] [--budget SEC]
        (GENERATE; resumable; state in state/zeros_q<q>_pair<n>.json, one per
         conjugate pair {n, nbar} of Conrey numbers, nbar = n^(-1) mod q)

Characters are Conrey characters chi_q(n, .) as implemented by Arb / python-flint
(dirichlet_char(q, n)); chi_q(nbar, .) is the complex conjugate (checked
numerically for every pair at start-up).

What is certified (Arb ball arithmetic; every inequality holds for all points
of the balls involved):

 [Z] For each character chi and each listed decimal d: the Hardy function
     Z(t, chi) (REAL-valued for real t, |Z(t, chi)| = |L(1/2 + it, chi)|,
     continuous; Arb acb_dirichlet_hardy_z) has certified opposite signs at
     d - EPS and d + EPS, EPS = 1e-30.  So each window holds at least one zero
     1/2 + i*gamma of L(s, chi) with gamma in (d - EPS, d + EPS), gamma > 0.
     Windows of one character are checked disjoint and increasing.

 [N] Completeness, for the pair {chi, chibar} at a common height T.  Let
       Lambda(s, chi) = (q/pi)^((s+k)/2) Gamma((s+k)/2) L(s, chi),
     k = 0 (even) or 1 (odd); entire; Lambda(s, chi) = eps(chi) Lambda(1-s, chibar),
     |eps| = 1; conj Lambda(conj s, chi) = Lambda(s, chibar).  Hence
       Lambda(1 - conj s, chi) = eps(chi) * conj Lambda(s, chi)
     (reflection in the critical line), so on the boundary of
     R_T = (-1, 2) x (-T, T) the change of argument on the left half equals the
     change on the right half; and the lower part of the right half for chi is
     the upper part for chibar.  Therefore
       W(T) := #{zeros of Lambda(., chi) in R_T, with multiplicity}
             = ( D(chi) + D(chibar) ) / pi,
     D(chi) = change of arg Lambda(s, chi) along C: 2 -> 2 + iT -> 1/2 + iT.
     D(chi) = theta_k(T) + [Arg L(2 + iT, chi) - Arg L(2, chi)] + hor(T, chi):
       * theta_k(T) = (T/2) log(q/pi) + Im lgamma((1/2 + k)/2 + iT/2)
         (lgamma principal branch, analytic for Re > 0, zero argument at s = 2);
       * on Re s = 2: |L - 1| <= zeta(2) - 1 < 1, so Re L > 0 and the continuous
         argument is the principal one at both ends;
       * hor: the segment [1/2, 2] + iT is covered by closed sub-segments; on
         each, L is evaluated on a ball containing it and must exclude 0, so the
         change of arg is the principal Arg of the ratio of the end values.
     The ball for W must contain a unique integer.  The zeros of Lambda(., chi)
     in R_T are exactly the nontrivial zeros of L(s, chi) with |Im| < T; those
     with Im < 0 are the conjugates of the zeros of L(s, chibar) with Im > 0.
     So if W = K(chi) + K(chibar) (numbers of windows in (0, T)), then every
     window holds EXACTLY ONE zero, it is SIMPLE, and there is no other zero in
     R_T: none off the critical line, none at s = 1/2, none between windows.
     (Real chi: chibar = chi, W = 2K.)

Classical inputs not re-proved: analytic continuation and functional equation
of Lambda; nontrivial zeros lie in 0 < Re s < 1 and L(s, chi) != 0 on Re s = 2
boundary pieces as argued above; correctness of Arb (acb_dirichlet_l,
acb_dirichlet_hardy_z, lgamma).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

from flint import arb, acb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
STATE = os.path.join(HERE, "state")
EPS = "1e-30"
PREC_SCAN = 96
PREC_REFINE = 224
PREC_ARG = 128
BLOCK = 20


def sign_of(x: arb) -> int:
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


class Char:
    def __init__(self, q: int, n: int):
        self.q, self.n = q, n
        self.chi = dirichlet_char(q, n)
        assert self.chi.conductor() == q and not self.chi.is_principal()
        self.k = int(self.chi.parity())

    def Z(self, t: arb) -> arb:
        return self.chi.hardy_z(acb(t)).real

    def L(self, s: acb) -> acb:
        return self.chi.l(s)

    def refine(self, a: arb, b: arb) -> str:
        old = ctx.prec
        ctx.prec = PREC_REFINE
        try:
            fa, fb = self.Z(a), self.Z(b)
            sa = sign_of(fa)
            lo, hi = a, b
            for _ in range(6):
                mid = arb(((lo + hi) / 2).mid())
                fm = self.Z(mid)
                sm = sign_of(fm)
                if sm == 0:
                    break
                if sm == sa:
                    lo, fa = mid, fm
                else:
                    hi, fb = mid, fm
            x0, f0, x1, f1 = lo, fa, hi, fb
            for _ in range(60):
                denom = (f1 - f0)
                x2 = arb((x1 - f1 * (x1 - x0) / denom).mid())
                if not (x2 > lo and x2 < hi):
                    x2 = arb(((lo + hi) / 2).mid())
                f2 = self.Z(x2)
                s2 = sign_of(f2)
                if s2 == 0:
                    x1 = x2
                    break
                if s2 == sa:
                    lo = x2
                else:
                    hi = x2
                x0, f0, x1, f1 = x1, f1, x2, f2
                if abs(float((x1 - x0).mid())) < 1e-45:
                    break
            return arb(x1.mid()).str(45, radius=False)
        finally:
            ctx.prec = old

    def certify_window(self, dstr: str) -> bool:
        old = ctx.prec
        ctx.prec = PREC_REFINE
        try:
            d = arb(dstr)
            e = arb(EPS)
            sm = sign_of(self.Z(d - e))
            sp = sign_of(self.Z(d + e))
            return sm != 0 and sp != 0 and sm == -sp
        finally:
            ctx.prec = old

    def theta(self, T: arb) -> arb:
        qpi = (arb(self.q) / arb.pi()).log()
        return T / 2 * qpi + acb((arb(1) / 2 + self.k) / 2, T / 2).lgamma().imag

    def horizontal_arg(self, T: arb):
        total = arb(0)
        nev = 0
        sig_hi = 2.0
        step = 1.0 / 16
        vcache = {}

        def val(sig: float) -> acb:
            nonlocal nev
            if sig not in vcache:
                vcache[sig] = self.L(acb(arb(sig), T))
                nev += 1
            return vcache[sig]

        while sig_hi > 0.5:
            sig_lo = max(0.5, sig_hi - step)
            hullre = arb(sig_lo).union(arb(sig_hi))
            B = self.L(acb(hullre, T))
            nev += 1
            ok = B.is_finite() and not B.contains(acb(0))
            if ok:
                delta = (val(sig_lo) / val(sig_hi)).arg()
                if delta.is_finite() and float(delta.rad()) < 1e-6:
                    total += delta
                    sig_hi = sig_lo
                    step = min(step * 2, 1.0 / 16)
                    continue
            step /= 2
            if step < 2.0 ** -40:
                raise ArithmeticError("horizontal segment near sigma=%.6f does not certify" % sig_hi)
        return total, nev

    def D(self, T: arb):
        """Ball for the change of arg Lambda(s, chi) along 2 -> 2+iT -> 1/2+iT."""
        old = ctx.prec
        ctx.prec = PREC_ARG
        try:
            L2T = self.L(acb(arb(2), T))
            L2 = self.L(acb(arb(2)))
            if not (L2T.real > 0 and L2.real > 0):
                raise ArithmeticError("Re L(2+it) not certified positive (impossible)")
            vert = L2T.arg() - L2.arg()
            hor, nev = self.horizontal_arg(T)
            return self.theta(T) + vert + hor, nev
        finally:
            ctx.prec = old


def check_conjugate(c1: Char, c2: Char):
    ctx.prec = 64
    for a in range(1, c1.q):
        v1 = c1.chi(a)
        v2 = c2.chi(a)
        if not abs(complex(v1) - complex(v2).conjugate()) < 1e-12:
            raise ArithmeticError("characters %d, %d are not conjugate" % (c1.n, c2.n))


def scan_char(C: Char, i0: int, i1: int, den: int):
    """Certified signs of Z at grid points i/den, i0 <= i <= i1."""
    vals = []
    for i in range(i0, i1 + 1):
        v = C.Z(arb(i) / den)
        if sign_of(v) == 0:
            old = ctx.prec
            ctx.prec = 4 * PREC_SCAN
            v = C.Z(arb(i) / den)
            ctx.prec = old
            if sign_of(v) == 0:
                raise ArithmeticError("scan value sign not certified: char %d, t = %d/%d" % (C.n, i, den))
        vals.append(v)
    return vals


def run_pair(q: int, n: int, H: int, budget: float, logf):
    nbar = pow(n, -1, q)
    assert n <= nbar
    chars = [Char(q, n)] if n == nbar else [Char(q, n), Char(q, nbar)]
    if n != nbar:
        check_conjugate(chars[0], chars[1])
    else:
        assert chars[0].chi.is_real()
    path = os.path.join(STATE, "zeros_q%d_pair%d.json" % (q, n))
    if os.path.exists(path):
        st = json.load(open(path))
    else:
        st = {"q": q, "pair": [c.n for c in chars], "parity": chars[0].k, "den": 32,
              "top_index": 0, "zeros": {str(c.n): [] for c in chars}, "blocks": [], "seconds": 0.0}
    DEN0 = st["den"]

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    t_start = time.time()
    ctx.prec = PREC_SCAN
    if not st["blocks"]:
        log("# q=%d pair %s parity k=%d (order %d)" % (q, st["pair"], chars[0].k, chars[0].chi.order()))
    while st["top_index"] < H * DEN0:
        t0 = time.time()
        i0 = st["top_index"]
        target = min(H, ((i0 + DEN0 - 1) // DEN0 // BLOCK + 1) * BLOCK)
        assert target * DEN0 > i0
        mult = 1
        while True:
            den = DEN0 * mult
            j0, j1 = i0 * mult, target * den
            allvals = [scan_char(C, j0, j1, den) for C in chars]
            # block top: grid point (on the coarse grid) in the last unit maximising min |Z|
            cand = [k for k in range(1, len(allvals[0])) if (j0 + k) % mult == 0 and (j0 + k) > (target - 1) * den]
            ktop = max(cand, key=lambda k: min(abs(float(v[k].mid())) for v in allvals))
            jtop = j0 + ktop
            new = {}
            for C, vals in zip(chars, allvals):
                sg = [sign_of(v) for v in vals]
                lst = []
                for k in range(ktop):
                    if sg[k] != sg[k + 1]:
                        a, b = arb(j0 + k) / den, arb(j0 + k + 1) / den
                        dstr = C.refine(a, b)
                        if not C.certify_window(dstr):
                            raise ArithmeticError("window around %s does not certify (char %d)" % (dstr, C.n))
                        d, e = arb(dstr), arb(EPS)
                        if not (d - e > a and d + e < b):
                            raise ArithmeticError("window around %s leaves its scan cell" % dstr)
                        lst.append(dstr)
                new[str(C.n)] = lst
            T = arb(jtop) / den
            K = sum(len(st["zeros"][str(C.n)]) + len(new[str(C.n)]) for C in chars)
            if len(chars) == 1:
                K *= 2
            Ds = [C.D(T) for C in chars]
            if len(chars) == 1:
                W = 2 * Ds[0][0] / arb.pi()
            else:
                W = (Ds[0][0] + Ds[1][0]) / arb.pi()
            try:
                Wint = int(W.unique_fmpz())
            except (TypeError, ValueError):
                raise ArithmeticError("count ball %s does not determine an integer" % W)
            if Wint == K:
                break
            log("  (count mismatch W=%d, windows=%d at T=%d/%d with grid 1/%d: refining grid)"
                % (Wint, K, jtop, den, den))
            mult *= 4
            if mult > 4096:
                raise ArithmeticError("count mismatch persists: multiple or off-line zero?")
        dt = time.time() - t0
        for C in chars:
            st["zeros"][str(C.n)].extend(new[str(C.n)])
        itop_num, itop_den = jtop, den
        # jtop is on the coarse grid by construction
        assert jtop % mult == 0
        st["top_index"] = jtop // mult
        st["blocks"].append({"T": "%d/%d" % (jtop // mult, DEN0), "K": K, "W": Wint})
        st["seconds"] += dt
        tmp = path + ".tmp"
        json.dump(st, open(tmp, "w"))
        os.replace(tmp, path)
        log("BLOCK pair %s top T = %d/%d = %.5f : W = %s -> %d = windows(chi) + windows(chibar) "
            "[%s] (grid 1/%d, %d L-evals) %.1fs"
            % (st["pair"], jtop // mult, DEN0, jtop / den, W.str(10, radius=True), Wint,
               ", ".join("%d:%d" % (C.n, len(st["zeros"][str(C.n)])) for C in chars), den,
               sum(d[1] for d in Ds), dt))
        if time.time() - t_start > budget:
            log("# budget reached for pair %s" % st["pair"])
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("H", type=int)
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--budget", type=float, default=1e9)
    args = ap.parse_args()
    q = args.q
    os.makedirs(STATE, exist_ok=True)
    reps = sorted({min(n, pow(n, -1, q)) for n in range(2, q)})
    if args.pairs:
        sel = [int(x) for x in args.pairs.split(",")]
        assert all(s in reps for s in sel), reps
        reps = sel
    logf = open(os.path.join(HERE, "zeros_q%d.cert.log" % q), "a")
    for n in reps:
        run_pair(q, n, args.H, args.budget, logf)


if __name__ == "__main__":
    main()
