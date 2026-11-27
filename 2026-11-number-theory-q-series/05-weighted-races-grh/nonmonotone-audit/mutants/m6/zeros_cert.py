#!/usr/bin/env python3
"""
zeros_cert.py -- certified zero data for L(s, chi), chi the ODD real primitive
character of conductor q (q = 163: the Kronecker symbol (-163/.)).

Usage:   nice -n 10 python3 zeros_cert.py q H [block]
         (resumable: finished blocks are stored in zeros_q<q>.state.json)

Output:  zeros_q<q>.txt        one decimal d per line; the certified window is
                               [d - 1e-20, d + 1e-20]
         zeros_q<q>.cert.txt   human-readable certificate log (appended)

What is certified (python-flint / Arb ball arithmetic; "certified" means the
stated inequality holds for every point of the balls involved):

 [Z]  Localisation.  The Hardy function Z(t) = Z(t, chi) (real for real t,
      |Z(t)| = |L(1/2 + it, chi)|, continuous) has certified opposite signs at
      d - 1e-20 and d + 1e-20 for every listed d.  Hence each window contains
      at least one zero of L on the critical line (an odd number, counted
      with multiplicity).  Windows are checked to be pairwise disjoint.

 [N]  Completeness and simplicity, by the argument principle for the
      completed function
          Lambda(s) = (q/pi)^((s+1)/2) Gamma((s+1)/2) L(s, chi)      (chi odd),
      which is entire and satisfies Lambda(s) = Lambda(1-s) (root number +1
      for a real primitive character: Gauss) and Lambda(conj s) = conj Lambda(s).
      Let R_T = (-1, 2) x (-T, T).  By the two symmetries the change of
      arg Lambda along the four quarters of the boundary of R_T is the same, so

          #{zeros of Lambda in R_T, with multiplicity} = (2/pi) * Delta_C arg Lambda,

      C = the path 2 -> 2 + iT -> 1/2 + iT.  Along C:
        * arg of the Gamma factor: Im[ ((s+1)/2) log(q/pi) + lgamma((s+1)/2) ]
          (lgamma = the principal branch, analytic for Re > 0, so its
          imaginary part IS the continuous argument); it vanishes at s = 2, so
          the change is theta(T) := (T/2) log(q/pi) + Im lgamma(3/4 + iT/2);
        * vertical part of arg L: |L(2+it) - 1| <= zeta(2) - 1 < 1, so
          Re L(2+it) > 0 and the continuous argument stays in (-pi/2, pi/2):
          the change is the principal Arg L(2 + iT)  (and L(2) > 0);
        * horizontal part: the segment [1/2, 2] + iT is covered by closed
          sub-segments; on each one L is evaluated on a ball containing the
          sub-segment and the result is required not to contain 0, so the
          change of arg along it is < pi in modulus and equals the principal
          Arg of the ratio of the endpoint values.
      The resulting ball for (2/pi) Delta_C arg Lambda must contain a unique
      integer W.  If W = 2K, where K is the number of certified windows in
      (0, T), then (each window and its mirror image hold >= 1 zero) every
      window contains EXACTLY ONE zero, that zero is SIMPLE, and Lambda has no
      other zero in R_T: none off the critical line, none on the real axis (so
      L(1/2, chi) != 0 and there is no real zero in (0, 1)), none between the
      windows.  Since all zeros of Lambda lie in 0 <= Re s <= 1 (classical),
      the list is the complete list of zeros with |Im s| < T.

 [0]  L(1/2, chi) > 0 is also certified directly (ball strictly positive).

Classical inputs NOT re-proved here: analytic continuation and the functional
equation of Lambda with root number +1; zeros of Lambda lie in the critical
strip; correctness of Arb's acb_dirichlet_l / acb_dirichlet_hardy_z / lgamma.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

from flint import arb, acb, ctx

from common import real_character

HERE = os.path.dirname(os.path.abspath(__file__))
DEN = 64                 # scan grid t = i / DEN (exact dyadic rationals)
EPS = "1e-20"
PREC_SCAN = 96
PREC_REFINE = 192
PREC_ARG = 128


def log(fh, msg):
    print(msg, flush=True)
    fh.write(msg + "\n")
    fh.flush()


def sign_of(x: arb) -> int:
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


class Cert:
    def __init__(self, q: int):
        self.q = q
        self.chi, self.d = real_character(q, 1)

    def Z(self, t: arb) -> arb:
        return self.chi.hardy_z(acb(t)).real

    def Zgrid(self, i: int, den: int = DEN) -> arb:
        return self.Z(arb(i) / den)

    # -- zero refinement ----------------------------------------------------
    def refine(self, a: arb, b: arb) -> str:
        """a < b exact points with certified opposite signs of Z.  Returns a
        30-digit decimal d (not itself certified here; see certify_window)."""
        old = ctx.prec
        ctx.prec = PREC_REFINE
        try:
            fa, fb = self.Z(a), self.Z(b)
            sa = sign_of(fa)
            lo, hi = a, b
            # bisection to width 2^-12 of the bracket, then safeguarded secant
            for _ in range(8):
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
            for _ in range(40):
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
                if abs(float((x1 - x0).mid())) < 1e-32:
                    break
            return arb(x1.mid()).str(32, radius=False)
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

    # -- argument principle ---------------------------------------------------
    def theta(self, T: arb) -> arb:
        qpi = (arb(self.q) / arb.pi()).log()
        return T / 2 * qpi - acb(arb(3) / 4, T / 2).lgamma().imag

    def horizontal_arg(self, T: arb):
        """Continuous change of arg L(sigma + iT) for sigma from 2 down to 1/2."""
        total = arb(0)
        nev = 0
        sig_hi = 2.0          # floats with exact binary values only (dyadic steps)
        step = 1.0 / 16
        vcache = {}

        def val(sig: float) -> acb:
            nonlocal nev
            if sig not in vcache:
                vcache[sig] = self.chi.l(acb(arb(sig), T))
                nev += 1
            return vcache[sig]

        minabs = None
        while sig_hi > 0.5:
            sig_lo = max(0.5, sig_hi - step)
            hullre = arb(sig_lo).union(arb(sig_hi))
            B = self.chi.l(acb(hullre, T))
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

    def count(self, T: arb):
        """Ball for W = #zeros of Lambda in (-1,2) x (-T,T), with multiplicity."""
        old = ctx.prec
        ctx.prec = PREC_ARG
        try:
            L2 = self.chi.l(acb(arb(2), T))
            if not L2.real > 0:
                raise ArithmeticError("Re L(2+iT) not certified positive (impossible)")
            vert = L2.arg()                       # principal Arg, in (-pi/2, pi/2)
            hor, nev = self.horizontal_arg(T)
            th = self.theta(T)
            W = 2 * (th + vert + hor) / arb.pi()
            return W, th, vert, hor, nev
        finally:
            ctx.prec = old


def main():
    q = int(sys.argv[1])
    H = int(sys.argv[2])
    block = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    state_path = os.path.join(HERE, "zeros_q%d.state.json" % q)
    zero_path = os.path.join(HERE, "zeros_q%d.txt" % q)
    cert_path = os.path.join(HERE, "zeros_q%d.cert.txt" % q)
    C = Cert(q)
    ctx.prec = PREC_SCAN
    if os.path.exists(state_path):
        with open(state_path) as fh:
            st = json.load(fh)
    else:
        st = {"q": q, "top_index": 0, "zeros": [], "blocks": [], "seconds": 0.0}
    fh = open(cert_path, "a")
    t_start = time.time()
    if not st["blocks"]:
        log(fh, "# zeros_cert.py  q=%d  character number %d = Kronecker symbol (%d/.), odd, "
                "real, primitive (checked for n <= 2q)" % (q, C.chi.number(), C.d))
        old = ctx.prec
        ctx.prec = 192
        L12 = C.chi.l(acb(arb(1) / 2)).real
        ctx.prec = old
        if not L12 > 0:
            raise ArithmeticError("L(1/2) not certified positive")
        log(fh, "# [0] L(1/2, chi) = %s  certified > 0" % L12.str(30, radius=True))
        z0 = C.Zgrid(0)
        if sign_of(z0) == 0:
            raise ArithmeticError("Z(0) sign not certified")

    while st["top_index"] <= (H - 1) * DEN:
        t0 = time.time()
        i0 = st["top_index"]
        # next block top: the first multiple of `block` at least one unit above i0
        target = min(H, ((i0 + DEN - 1) // DEN // block + 1) * block)
        # scan grid i0 .. target*DEN ; the block top is the grid point in the
        # last unit with the largest |Z| (far from zeros, cheap contour)
        idx = list(range(i0, target * DEN + 1))
        vals = [C.Zgrid(i) for i in idx]
        sg = [sign_of(v) for v in vals]
        if 0 in sg:
            raise ArithmeticError("a scan value has uncertified sign (index %d)" % idx[sg.index(0)])
        last = [k for k in range(len(idx)) if idx[k] > (target - 1) * DEN]
        ktop = max(last, key=lambda k: abs(float(vals[k].mid())))
        itop = idx[ktop]
        new = []
        for k in range(ktop):
            if sg[k] != sg[k + 1]:
                dstr = C.refine(arb(idx[k]) / DEN, arb(idx[k + 1]) / DEN)
                if not C.certify_window(dstr):
                    raise ArithmeticError("window around %s does not certify" % dstr)
                d = arb(dstr)
                e = arb(EPS)
                if not (d - e > arb(idx[k]) / DEN and d + e < arb(idx[k + 1]) / DEN):
                    raise ArithmeticError("window around %s leaves its scan cell" % dstr)
                new.append(dstr)
        T = arb(itop) / DEN
        K = len(st["zeros"]) + len(new)
        W, th, vert, hor, nev = C.count(T)
        try:
            Wint = int(W.unique_fmpz())
        except (TypeError, ValueError):
            Wint = None
        if Wint is None:
            raise ArithmeticError("count ball %s does not determine an integer" % W)
        dt = time.time() - t0
        log(fh, "BLOCK top T = %d/%d = %.6f : %d new windows, K = %d in (0,T); "
                "theta/pi = %s, ArgL(2+iT)/pi = %s, horizontal/pi = %s (%d L-evals); "
                "W = %s -> %d ; |Z(T)| = %.3f ; %.1fs"
            % (itop, DEN, itop / DEN, len(new), K,
               (th / arb.pi()).str(12, radius=False), (vert / arb.pi()).str(8, radius=False),
               (hor / arb.pi()).str(8, radius=False), nev, W.str(12, radius=True), Wint,
               abs(float(vals[ktop].mid())), dt))
        if Wint != 2 * K:
            log(fh, "*** COUNT MISMATCH: W = %d but 2K = %d: zeros missed by the scan "
                    "(or multiple / off-line zeros). Block NOT accepted." % (Wint, 2 * K))
            sys.exit(1)
        log(fh, "      => exactly %d zeros of Lambda in (-1,2)x(-T,T): the %d windows hold one "
                "SIMPLE zero each, no other zeros (none real, none off the line) up to T." % (Wint, K))
        st["zeros"].extend(new)
        st["top_index"] = itop
        st["blocks"].append({"itop": itop, "K": K, "W": Wint})
        st["seconds"] += dt
        with open(state_path, "w") as sfh:
            json.dump(st, sfh)
        with open(zero_path, "w") as zfh:
            zfh.write("# q=%d (Kronecker (%d/.)): %d ordinates, each window [d-1e-20,d+1e-20] holds exactly\n"
                      "# one simple zero; complete for 0 < gamma < T_cert = %d/%d (see zeros_q%d.cert.txt)\n"
                      % (q, C.d, K, itop, DEN, q))
            for z in st["zeros"]:
                zfh.write(z + "\n")
        if time.time() - t_start > 600:
            log(fh, "# time budget for this invocation reached; rerun to continue")
            break
    log(fh, "# state: %d zeros certified complete up to T_cert = %d/%d = %.6f ; cumulative %.1fs"
        % (len(st["zeros"]), st["top_index"], DEN, st["top_index"] / DEN, st["seconds"]))


if __name__ == "__main__":
    main()
