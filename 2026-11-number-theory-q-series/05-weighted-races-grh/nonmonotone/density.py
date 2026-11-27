#!/usr/bin/env python3
"""
density.py -- fail-closed Arb enclosures of the bias density

    delta_chi(M) = P( 1 + sum_{gamma>0} c_gamma(M) cos(theta_gamma) > 0 ),
    c_gamma(M) = 4M / sqrt(M^2 + gamma^2),   theta_gamma i.i.d. uniform,

for an ODD real primitive character chi of conductor q, at rational M.
This is the generalisation of prototype/explore/race_density_arb.py (written
for chi_3, chi_4 and integer m = M - 1/2) to an arbitrary such character and
arbitrary rational M.  Under GRH(chi) + LI(chi) the number delta_chi(M) is the
logarithmic density of the weighted race with weight exponent m = M - 1/2
(paper/weighted-races.tex); this program does NOT prove that statement, it
encloses the probability above.

Usage:  nice -n 10 python3 density.py --q 163 --M 3/5 6/5 [--out results_raw.txt]

Inputs taken from elsewhere
  * zeros_q<q>.txt + zeros_q<q>.cert.txt (zeros_cert.py): windows
    [d - 1e-20, d + 1e-20], each holding exactly one simple zero ordinate,
    and COMPLETENESS of the list for 0 < gamma < T_cert.  Only the following
    consequence is used: every zero ordinate gamma > 0 that is not in a
    listed window satisfies gamma > gamma_last (the last listed window).
    THE FILE IS BOUND TO ITS CERTIFICATE: the program refuses to start unless
    the sha256 of the zero file is recorded in a passing HASH-RECORD line of
    the certificate log (written by zeros_cert.py only after the windows and
    the count W = 2K were certified for exactly that file).  The hash is
    printed in every DELTA line.  (For q = 3, 4 the lists of
    prototype/explore/ are used, after `zeros_cert.py q --verify --file ...`.)
    Second, heuristic guard: the exact tail variance V must agree with the
    Riemann-von Mangoldt estimate of the unlisted zeros (rvm_guard); a list
    with a missing low zero fails it grossly.  The guard can only refuse.
    With --uncertified-zeros the hash requirement is waived and every output
    line is labelled UNCERTIFIED instead of DELTA (for experiments only).
  * GRH(chi): all zeros are 1/2 + i gamma, gamma real.  Used for the variance
    identity (below).  L(1/2, chi) != 0 is certified by zeros_cert.py.

Ingredients (each is an inequality proved in the comments / in NONMONOTONE.md)

 1. Gil-Pelaez.  Phi(t) = prod_gamma J_0(c_gamma t) is integrable (item 5)
    and  delta - 1/2 = (1/pi) int_0^inf Phi(t) sin(t)/t dt.

 2. Exact tail variance.  Under GRH and L(1/2,chi) != 0, summing over zeros
    WITH MULTIPLICITY,
        sum_{gamma>0} c_gamma^2/2 = sigma^2(M) := 4M (xi'/xi)(M + 1/2, chi)
                 = 4M [ (1/2) log(q/pi) + (1/2) psi((M + 3/2)/2) + (L'/L)(M + 1/2) ].
    The listed zeros are simple, so V := sigma^2 - sum_{listed} c^2/2 is
    EXACTLY sum over unlisted zeros (with multiplicity) of c^2/2 >= 0, and
    all unlisted zeros have gamma > gamma_last, i.e. c_gamma < c_top.
    (Under LI the unlisted zeros are simple too; the tail bounds below do not
    care: a zero of multiplicity k is k factors.)
    sigma^2 is enclosed in Arb: L' by the complex step Im L(s + ih)/h with a
    Cauchy/Abel remainder, cross-checked against Arb's power-series L-jet.

 3. Unseen product.  For 0 <= x <= 1:  exp(-x^2/4 - C4 x^4) <= J_0(x) <= exp(-x^2/4),
    C4 = 0.018 >= -log J_0(1) - 1/4 (checked in Arb at start-up).  Hence for
    0 <= t <= 1/c_top the product Psi over the unlisted zeros satisfies
        exp(-V t^2/2 - cb t^4) <= Psi(t) <= exp(-V t^2/2),
        cb := C4 * c_top^2 * 2V  >=  C4 * sum_unlisted c^4.

 4. On [0, T] (T <= 1/c_top): with P = finite product over the listed zeros,
        Phi sin t/t = P g m sin t/t  +  err,   g = exp(-V t^2/2),
        m = (1 + exp(-cb t^4))/2,   |err| <= |P| g (1 - exp(-cb t^4))/2 * |sin t/t|.
    The first term is integrated by composite 16-point Gauss-Legendre with
    certified nodes/weights and the Cauchy-estimate remainder; the second
    (E1) is bounded by interval subdivision (sup over cell * width).

 5. Tail in t.  |Phi(t)| <= prod_{j<=K} |J_0(c_j t)| <= prod_{j<=K} (c_j t)^(-1/2)
    (|J_0(x)| <= x^(-1/2), x > 0: Watson 13.74), so
        | int_T^inf Phi sin t/t dt | <= prod_{j<=K} c_j^(-1/2) * T^(-K/2) / (K/2) =: E2,
    for ANY K >= 1 (K is optimised; any K is rigorous).

Every printed enclosure is  1/2 + (I +- (rad(I) + Egl + E1 + E2))/pi.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from fractions import Fraction

from flint import acb, acb_series, arb, ctx
from mpmath import mp

from common import (real_character, character_values, max_partial_sum,
                    load_certified_zeros, HZ_SHA256)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CTX_BITS = 192
NODE_N = 16
C4 = arb(18) / 1000
TAIL_DENSITY_TARGET = 1e-15
T_DEN = 64


def hull(lo: arb, hi: arb) -> arb:
    if lo.lower() > hi.upper():
        raise ArithmeticError("invalid interval hull")
    return lo.union(hi)


def inflate(x: arb, radius: arb) -> arb:
    r = radius.upper()
    return x + hull(-r, r)


def require_finite(x, label: str) -> None:
    if not x.is_finite():
        raise ArithmeticError("non-finite Arb value: %s: %s" % (label, x))


def frac_to_arb(fr: Fraction) -> arb:
    return arb(fr.numerator) / arb(fr.denominator)


# ---------------------------------------------------------------------------
# sigma^2(M) = 4M (xi'/xi)(M + 1/2)
# ---------------------------------------------------------------------------

class Character:
    def __init__(self, q: int):
        self.q = q
        self.chi, self.d = real_character(q, 1)      # odd: a = 1
        self.smax = max_partial_sum(character_values(self.chi, q))

    def logderiv_L(self, s: arb) -> tuple[arb, arb, arb]:
        """Enclosures (L'/L)(s), L(s), L'(s) for real s > 1/4 + (a little)."""
        chi = self.chi
        h = arb(2) ** -60
        R = arb(1) / 4
        if not s - R > arb(1) / 10:
            raise ArithmeticError("s too small for the Abel bound")
        Ls = chi.l(acb(s)).real
        if not Ls > 0:
            raise ArithmeticError("L(s) not certified positive at s=%s" % s)
        Lh = chi.l(acb(s, h))
        dL = Lh.imag / h
        # |L(z)| <= smax |z| / Re z on Re z > 0 (Abel summation, partial sums of
        # chi bounded by smax, computed exactly).  L is real on the real axis,
        # so Im L(s+ih)/h = L'(s) - L'''(s) h^2/3! + ...; Cauchy on |z-s| <= R:
        ML = arb(self.smax) * (s + R) / (s - R)
        rem = ML * h ** 2 / (R ** 3 * (1 - (h / R) ** 2))
        dL = inflate(dL, rem)
        # independent enclosure: Arb's power-series jet of L
        jet = acb_series.dirichlet_l(acb_series([acb(s), acb(1)], prec=2), chi)
        L0j, L1j = jet[0].real, jet[1].real
        if not (L0j.overlaps(Ls) and L1j.overlaps(dL)):
            raise ArithmeticError("complex-step and jet enclosures of L, L' disagree")
        return dL / Ls, Ls, dL

    def sigma2(self, M: arb) -> arb:
        s = M + arb(1) / 2
        lpl, _, _ = self.logderiv_L(s)
        out = 4 * M * ((arb(self.q) / arb.pi()).log() / 2
                       + ((s + 1) / 2).digamma() / 2 + lpl)
        require_finite(out, "sigma2")
        if not out > 0:
            raise ArithmeticError("sigma2 not certified positive")
        return out


# ---------------------------------------------------------------------------
# certified Gauss-Legendre nodes (as in race_density_arb.py)
# ---------------------------------------------------------------------------

def legendre(n: int, x: arb) -> arb:
    p0 = arb(1)
    if n == 0:
        return p0
    p1 = x
    for k in range(1, n):
        p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
    return p1


def legendre_derivative(n: int, x: arb) -> arb:
    return n * (x * legendre(n, x) - legendre(n - 1, x)) / (x * x - 1)


def certified_gauss_legendre(n: int) -> list[tuple[arb, arb]]:
    mp.dps = 90
    xs, _ = mp.gauss_quadrature(n, "legendre")
    rad_mp = mp.mpf("1e-45")
    out = []
    for xm in xs:
        mid = mp.nstr(xm, 85)
        xb = arb(mid, "1e-45")
        left = arb(mp.nstr(xm - rad_mp, 85))
        right = arb(mp.nstr(xm + rad_mp, 85))
        pl, pr = legendre(n, left), legendre(n, right)
        opposite = (pl.upper() < 0 and pr.lower() > 0) or (pl.lower() > 0 and pr.upper() < 0)
        dp = legendre_derivative(n, xb)
        if (not xb.contains(left) or not xb.contains(right) or not opposite
                or dp.contains(0) or not legendre(n, xb).contains(0)):
            raise ArithmeticError("failed to isolate Gauss node near " + mid)
        w = 2 / ((1 - xb * xb) * dp * dp)
        if not w > 0:
            raise ArithmeticError("nonpositive Gauss weight enclosure")
        out.append((xb, w))
    if len(out) != n:
        raise ArithmeticError("wrong number of Gauss nodes")
    for (x, _), (y, _) in zip(out, out[1:]):
        if not (-1 < x.lower() and x.upper() < y.lower() and y.upper() < 1):
            raise ArithmeticError("Gauss node balls not strictly ordered")
    if not sum((w for _, w in out), arb(0)).contains(2):
        raise ArithmeticError("Gauss weights do not enclose total mass 2")
    return out


# ---------------------------------------------------------------------------
# the enclosure
# ---------------------------------------------------------------------------

def rvm_tail_estimate(q: int, M: float, gamma_last: float) -> float:
    """FLOAT (non-rigorous) Riemann-von Mangoldt estimate of the unlisted
    variance  sum_{gamma > gamma_last} 8M^2/(M^2+gamma^2): zero density
    (1/2pi) log(q t / 2pi), starting half a mean spacing above gamma_last."""
    with mp.workdps(20):
        dens = lambda t: mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi)
        G = gamma_last + 0.5 / float(dens(gamma_last))
        val = 8 * M * M * mp.quad(lambda t: dens(t) / (M * M + t * t), [G, 10 * G, 100 * G, mp.inf])
        return float(val)


def rvm_guard(q: int, M: arb, zeros: list[arb], V: arb) -> tuple[float, float, float]:
    """Plausibility guard (audit G1): refuse unless the exact tail variance V
    is within  3 * (variance of one zero at gamma_last) + 1%  of the RvM
    estimate.  It is a heuristic and can only REFUSE; it catches gross errors
    (a missing low zero: V = 7.21 instead of 0.0349 when gamma_1 of conductor
    163 is deleted), not the omission of a single high zero -- that is what
    the hash binding to the certificate is for."""
    Mf = float(M.mid())
    gl = float(zeros[-1].mid())
    est = rvm_tail_estimate(q, Mf, gl)
    tol = 3 * 8 * Mf * Mf / (Mf * Mf + gl * gl) + 0.01 * est
    Vf = float(V.mid())
    if not abs(Vf - est) <= tol:
        raise ArithmeticError(
            "REFUSED (RvM guard): tail variance V = %.6g but the Riemann-von Mangoldt estimate of "
            "the unlisted zeros is %.6g (tolerance %.3g): the zero list is incomplete or wrong"
            % (Vf, est, tol))
    return Vf, est, tol


def tail_bound(bs: list[arb], T: arb) -> tuple[arb, int]:
    """min over K of  prod_{j<=K} b_j^(-1/2) T^(-K/2) / (K/2)  (each K rigorous)."""
    best = None
    bestK = 0
    A = arb(1)
    for k, b in enumerate(bs, 1):
        A = A / (b * T).sqrt()
        if k >= 3:
            val = A / (arb(k) / 2)
            if best is None or float(val.upper()) < float(best.upper()):
                best, bestK = val, k
    require_finite(best, "E2")
    return best, bestK


def choose_T(bs: list[arb], tmax: arb) -> arb:
    """A dyadic T <= tmax, aiming at E2/pi <= TAIL_DENSITY_TARGET (heuristic
    choice only: the E2 actually used is recomputed in Arb for the chosen T)."""
    bm = [float(b.mid()) for b in bs]
    cap = math.floor(float(tmax.lower()) * T_DEN) / T_DEN
    cands = []
    logprod = 0.0
    for k, b in enumerate(bm, 1):
        logprod += math.log(b)
        if k >= 8 and k % 4 == 0:
            logT = (-0.5 * logprod - math.log(math.pi * (k / 2) * TAIL_DENSITY_TARGET)) * (2 / k)
            cands.append(math.exp(logT))
    Tf = min(cands) if cands else cap
    Tf = math.ceil(Tf * T_DEN) / T_DEN
    Tf = min(Tf, cap)
    T = arb(int(round(Tf * T_DEN))) / T_DEN
    if not (T > 0 and T <= tmax):
        raise ArithmeticError("no admissible truncation point")
    return T


class Model:
    def __init__(self, bs, V, cb, T, keep):
        self.bs = bs
        self.V = V
        self.cb = cb
        self.T = T
        self.keep = min(keep, len(bs))
        # suffix sums of b_j^2 (for the Gaussian majorant of dropped factors)
        suf = [arb(0)] * (len(bs) + 1)
        for j in range(len(bs) - 1, -1, -1):
            suf[j] = suf[j + 1] + bs[j] * bs[j]
        self.suf = suf
        self.bup = [float(b.upper()) * (1 + 1e-12) for b in bs]

    def main_integrand(self, t: arb) -> arb:
        p = arb(1)
        for b in self.bs:
            p *= (b * t).bessel_j(0)
        require_finite(p, "Bessel product")
        g = (-self.V * t * t / 2).exp()
        mid = (1 + (-self.cb * t ** 4).exp()) / 2
        return p * g * mid * t.sin() / t

    def e1_cell(self, a: arb, b: arb) -> tuple[arb, arb]:
        """(lower, upper) bounds of sup-type enclosure of
        |P(t)| g(t) (1 - exp(-cb t^4))/2 * |sin t / t|  on [a, b], times (b-a).
        Factors j >= keep: |J_0| <= 1 always, and 0 < J_0(x) <= exp(-x^2/4)
        when x = b_j t <= 1 (<= j_{0,1}); the latter is used for those j with
        b_j * b <= 1 (the b_j decrease, so these form a suffix)."""
        t = hull(a, b)
        p = arb(1)
        for bj in self.bs[:self.keep]:
            p *= (bj * t).bessel_j(0)
        y = abs(p)
        n = self.keep
        bu = float(b.upper())
        while n < len(self.bs) and self.bup[n] * bu > 1.0:
            n += 1
        if n < len(self.bs):
            if not (self.bs[n] * b <= 1):
                n += 1
            if n < len(self.bs):
                if not (self.bs[n] * b <= 1):
                    raise ArithmeticError("suffix selection failed")
                y = y * (-self.suf[n] * t * t / 4).exp()
        y = y * (-self.V * t * t / 2).exp() * (1 - (-self.cb * t ** 4).exp()) / 2
        if a >= 1:
            y = y / a                      # |sin t / t| <= 1/t <= 1/a
        require_finite(y, "E1 integrand")
        w = b - a
        up = y.upper() * w
        lo = y.lower()
        lo = (lo if lo > 0 else arb(0)) * w
        return lo, up


def composite_gauss(model: Model, nodes) -> tuple[arb, arb, int]:
    A = sum((b.upper() for b in model.bs), arb(0))
    r = min(arb(1) / 8, arb(4) / A)
    N = max(1, math.ceil(float((model.T / (arb(3) * r / 4)).upper())))
    h = model.T / N
    half = h / 2
    total = arb(0)
    for j in range(N):
        center = (arb(j) + arb(1) / 2) * h
        sub = arb(0)
        for x, w in nodes:
            sub += w * model.main_integrand(center + half * x)
        total += half * sub
    # Bound of |integrand| on every disc |z - xi| <= r, xi in [0, T], z = x+iy:
    #   |J_0(b z)| <= exp(b |y|);  |exp(-V z^2/2)| <= exp(V y^2/2);
    #   |sin z / z| <= exp(|y|);   Re z^4 >= -8 y^4  =>  |exp(-cb z^4)| <= exp(8 cb y^4).
    Vup = model.V.upper()
    cbup = model.cb.upper()
    Mb = (A * r + Vup * r * r / 2 + r).exp() * (1 + (8 * cbup * r ** 4).exp()) / 2
    n = len(nodes)
    Q = arb(math.factorial(n)) ** 4 / ((2 * n + 1) * arb(math.factorial(2 * n)) ** 2)
    if not h / r < 1:
        raise ArithmeticError("Gauss cell too long")
    err = N * Mb * h * Q * (h / r) ** (2 * n)
    require_finite(err, "Gauss remainder")
    return total, err, N


def interval_e1(model: Model, abs_tol: float, rel_tol: float, max_cells: int):
    n0 = max(256, int(math.ceil(float(model.T.upper()) * 32)))
    cells = []
    for j in range(n0):
        a = model.T * j / n0
        b = model.T * (j + 1) / n0
        lo, up = model.e1_cell(a, b)
        cells.append((float((up - lo).upper()), a, b, lo, up))
    while True:
        tot_lo = sum((c[3] for c in cells), arb(0))
        tot_up = sum((c[4] for c in cells), arb(0))
        gap = float((tot_up - tot_lo).upper())
        if gap <= max(abs_tol, rel_tol * float(tot_up.upper())) or len(cells) >= max_cells:
            return tot_up, len(cells), gap
        cells.sort(key=lambda c: -c[0])
        acc = 0.0
        k = 0
        while k < len(cells) and acc < 0.6 * gap:
            acc += cells[k][0]
            k += 1
        refine, cells = cells[:k], cells[k:]
        for _, a, b, _, _ in refine:
            m = (a + b) / 2
            for (u, v) in ((a, m), (m, b)):
                lo, up = model.e1_cell(u, v)
                cells.append((float((up - lo).upper()), u, v, lo, up))


def enclose(char: Character, M: arb, zeros: list[arb], nodes, keep: int = 40,
            e1_abs: float = 1e-10, e1_rel: float = 0.05, max_cells: int = 60000) -> dict:
    t0 = time.time()
    bs = [4 * M / (M * M + g * g).sqrt() for g in zeros]
    s2 = char.sigma2(M)
    finite_var = sum((b * b / 2 for b in bs), arb(0))
    V = s2 - finite_var
    if not V > 0:
        raise ArithmeticError("tail variance not certified positive: %s" % V)
    guard = rvm_guard(char.q, M, zeros, V)
    c_top = bs[-1].upper()                  # every unlisted zero has c < c_top
    tmax = 1 / c_top
    S4 = c_top ** 2 * 2 * V
    cb = (C4 * S4).upper()
    T = choose_T(bs, tmax)
    if not (T * c_top <= 1):
        raise ArithmeticError("T exceeds the range of the tail model")
    E2, K = tail_bound(bs, T)
    model = Model(bs, V, cb, T, keep)
    I, Egl, ncells = composite_gauss(model, nodes)
    E1, e1cells, e1gap = interval_e1(model, e1_abs, e1_rel, max_cells)
    raw_radius = E1.upper() + E2.upper() + Egl.upper() + I.rad()
    center = arb(1) / 2 + arb(I.mid()) / arb.pi()
    density = inflate(center, raw_radius / arb.pi())
    require_finite(density, "density")
    if density.lower() < 0 or density.upper() > 1:
        raise ArithmeticError("invalid probability enclosure: %s" % density)
    return dict(M=M, sigma2=s2, V=V, c1=bs[0], ctop=c_top, T=T, K=K, E1=E1, E2=E2, Egl=Egl,
                Irad=I.rad(), gauss_cells=ncells, e1_cells=e1cells, density=density,
                nzeros=len(zeros), seconds=time.time() - t0, cb=cb, guard=guard,
                n_gt1=sum(1 for b in bs[:K] if b * T > 1))


def fixed12(k: int) -> str:
    sign = "-" if k < 0 else ""
    whole, frac = divmod(abs(k), 10 ** 12)
    return "%s%d.%012d" % (sign, whole, frac)


def directed12(x: arb, lower: bool) -> str:
    scale = arb(10 ** 12)
    rounded = (x.lower() * scale).floor() if lower else (x.upper() * scale).ceil()
    k = rounded.unique_fmpz()
    if k is None:
        # the endpoint is within rounding error of a 1e-12 grid point: step outward
        k = (rounded.lower().floor() if lower else rounded.upper().ceil()).unique_fmpz()
        if k is None:
            raise ArithmeticError("directed decimal rounding failed")
    out = fixed12(int(k))
    parsed = arb(out)
    if lower:
        if not (parsed <= x.lower()):
            out = fixed12(int(k) - 1)
            if not (arb(out) <= x.lower()):
                raise ArithmeticError("printed lower endpoint is not outward")
    else:
        if not (parsed >= x.upper()):
            out = fixed12(int(k) + 1)
            if not (arb(out) >= x.upper()):
                raise ArithmeticError("printed upper endpoint is not outward")
    return out


def zero_file(q: int) -> str:
    if q == 4:
        return os.path.join(REPO, "prototype", "explore", "chi4_zeros.txt")
    if q == 3:
        return os.path.join(REPO, "prototype", "explore", "chi3_zeros.txt")
    return os.path.join(HERE, "zeros_q%d.txt" % q)


def setup():
    ctx.prec = CTX_BITS
    c4_exact = -(arb(1).bessel_j(0)).log() - arb(1) / 4
    if not c4_exact < C4:
        raise ArithmeticError("C4 is not an upper bound of -log J0(1) - 1/4")
    return certified_gauss_legendre(NODE_N)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--M", nargs="+", required=True, help="rationals p/q (M = m + 1/2)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--nzeros", type=int, default=0, help="use only the first n listed zeros")
    ap.add_argument("--keep", type=int, default=40)
    ap.add_argument("--e1-abs", type=float, default=1e-10)
    ap.add_argument("--e1-rel", type=float, default=0.05)
    ap.add_argument("--zeros", default=None, help="zero file (default: the certified file of q)")
    ap.add_argument("--uncertified-zeros", action="store_true",
                    help="waive the hash/certificate requirement; output is labelled UNCERTIFIED")
    args = ap.parse_args()

    nodes = setup()
    char = Character(args.q)
    zpath = args.zeros or zero_file(args.q)
    zeros, meta = load_certified_zeros(args.q, zpath, allow_uncertified=args.uncertified_zeros)
    tag = "DELTA" if meta["certified"] else "UNCERTIFIED"
    if args.nzeros:
        zeros = zeros[:args.nzeros]      # a prefix of a complete list is complete below the next zero
    out = open(args.out, "a") if args.out else None

    def emit(line: str) -> None:
        print(line, flush=True)
        if out:
            out.write(line + "\n")
            out.flush()

    emit("# density.py  q=%d  chi = Kronecker (%d/.)  zeros used: %d (gamma_1 = %s, gamma_last = %s)  "
         "max|partial sum| = %d" % (args.q, char.d, len(zeros), zeros[0].str(12, radius=False),
                                    zeros[-1].str(12, radius=False), char.smax))
    if meta["certified"]:
        emit("# zero file %s  sha256=%s  bound to a passing certificate record (K=%d, T_cert=%s, W=%d, %s)%s"
             % (os.path.basename(zpath), meta["sha256"], meta["K"], meta["T_cert"], meta["W"], meta["mode"],
                "  = hypothesis (H-Z)" if HZ_SHA256.get(args.q) == meta["sha256"] else ""))
    else:
        emit("# *** zero file %s  sha256=%s  has NO certificate record: output below is UNCERTIFIED ***"
             % (os.path.basename(zpath), meta["sha256"]))
    for Ms in args.M:
        fr = Fraction(Ms)
        r = enclose(char, frac_to_arb(fr), zeros, nodes, keep=args.keep,
                    e1_abs=args.e1_abs, e1_rel=args.e1_rel)
        d = r["density"]
        emit("%s q=%d M=%s (m=%s)  lower=%s  upper=%s  width=%.2e  | sigma2=%s  V_tail=%s  c_1=%.6f  "
             "T=%s  K=%d  E1=%.2e  E2=%.2e  Egl=%.2e  radI=%.2e  gauss_cells=%d  e1_cells=%d  time=%.0fs  "
             "zeros_sha256=%s"
             % (tag, args.q, fr, fr - Fraction(1, 2), directed12(d, True), directed12(d, False),
                float((d.upper() - d.lower()).upper()), r["sigma2"].str(15, radius=False),
                r["V"].str(10, radius=False), float(r["c1"].mid()), r["T"].str(8, radius=False), r["K"],
                float(r["E1"].upper()), float(r["E2"].upper()), float(r["Egl"].upper()),
                float(r["Irad"].upper()), r["gauss_cells"], r["e1_cells"], r["seconds"], meta["sha256"]))
        emit("#   guard M=%s: V_tail=%.8f vs RvM estimate %.8f (tolerance %.2e) OK;  "
             "%d of the K=%d envelope factors have b_j*T > 1"
             % (fr, r["guard"][0], r["guard"][1], r["guard"][2], r["n_gt1"], r["K"]))


if __name__ == "__main__":
    try:
        main()
    except ArithmeticError as exc:
        print("FAIL-CLOSED (density.py): %s" % exc, flush=True)
        sys.exit(2)
