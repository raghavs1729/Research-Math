#!/usr/bin/env python3
"""AUDITOR'S independent rigorous enclosure of delta'(M) on a cell [Ma, Mb].

Written from the mathematical definition, WITHOUT reusing certify_derivative.py.
It shares with the audited code only (i) the zero files, (ii) the hypotheses
GRH + variance identity + completeness of the zero lists, and (iii) FLINT/Arb.
Everything else follows a different route:

                         audited code                    this file
  quadrature             interval Riemann sums            Arb's rigorous Gauss-Legendre
                         on [0,T1] (first order)          integrator acb.integral (Petras),
                                                          integrand enclosed on COMPLEX balls
  zeros treated by       all listed zeros (511/537)       only the first N' (20..130) zeros;
  explicit Bessel                                         N' chosen per cell
  factors
  remaining zeros        unlisted zeros only:             listed zeros N'+1..N AND unlisted
                         Gaussian * [e^{-C4 S4 t^4}, 1],  zeros: power-sum expansion
                         first term of phi                -log J0(z) = sum lambda_n/n z^{2n}
                                                          through z^{2 NMAX} (NMAX = 6) with
                                                          exact rational Rayleigh sums
                                                          lambda_n, listed power sums exact,
                                                          monotone-ratio remainder at z0=3/2
  closed forms           Dirichlet series over prime      Hurwitz-zeta power series
  sigma^2, D             powers n <= 20000 + digamma      (arb_series.zeta) and lgamma series
  Bessel on balls        real 2nd-order Taylor            complex Lipschitz / Taylor bounds with
                                                          |J_n(z)| <= e^{|Im z|}
  coefficients on a      c = 4/sqrt(1+(g/M)^2) on the     monotone hulls of endpoint values
  cell                   M-ball
  large-t remainder      |J0(x)| <= x^{-1/2}              |J0(x)| <= sqrt(2/(pi x)), own T1
  zero parsing           arb(str, "1e-20")                exact Fraction -> arb, + 1e-20

MATHEMATICS (all re-derived by the auditor; see AUDIT.md section 3)
-------------------------------------------------------------------
c_g(M) = 4M/sqrt(M^2+g^2), c_g' = 4 g^2/(M^2+g^2)^{3/2}.
Split the zeros into H' (first N' listed), R (listed, N' < index <= N) and
U (unlisted, g > g_N).  For |t| c_{g_{N'+1}}(Mb) <= z0 < j_{0,1}, t COMPLEX:

  Phi_M(t) = P(t) * Psi(t),   P = prod_{H'} J0(c t),
  Psi(t)  = exp( - sum_{n=1}^{NMAX} (lam_n/n) S_2n t^{2n} - r1(t) ),
  d/dM Psi = - Psi * B(t),
  B(t)    = sum_{n=1}^{NMAX} 2 lam_n U_2n t^{2n} + r2(t),
  S_2n = sum_{R u U} c^{2n},   U_2n = sum_{R u U} c^{2n-1} c',
  |r1(t)| <= rho1 (|t|/z0)^{2NMAX+2} S_{2NMAX+2},
  |r2(t)| <= rho2 |t|^{2NMAX+2} z0^{-(2NMAX+1)} U_{2NMAX+2},
  rho1 = -log J0(z0) - sum_{n<=NMAX} (lam_n/n) z0^{2n} >= 0,
  rho2 = J1(z0)/J0(z0) - sum_{n<=NMAX} 2 lam_n z0^{2n-1} >= 0
(positive Taylor coefficients).  S_2 = 2 V', U_2 = D' with
  V'(M) = sigma^2(M) - sum_{H'} c^2/2,   D'(M) = D(M) - sum_{H'} c c'
exact under GRH + completeness; both nondecreasing in M (D' because
g_{N'+1}^2 > 3 Mb^2, checked).  For n >= 2 the R-part of S_2n, U_2n is
evaluated from the list and the U-part lies in [0, cG^{2n-2} 2 V_T(Mb)] resp.
[0, cG^{2n-2} D_T(Mb)], cG = c_{g_N}(Mb), V_T, D_T the exact sums over U.

  F(t; M) = sinc(t) * Psi * ( dP - P * B ),   dP = -t sum_{H'} c' J1(c t) prod_{others} J0,

is, for each fixed M, an ENTIRE function of t; `integrand` returns for every
complex ball in the disc |t| <= Rdom a ball containing F(ball; M) for every M
in the cell, and a non-finite ball outside the disc.  Hence acb.integral
returns a ball containing int_0^{T1} F(t; M) dt for every M in the cell.

Large t:  |d/dM Phi| <= (t^2/2) D(M) max_g prod_{g' in S\\{g}} |J0(c_{g'} t)|,
|J0(x)| <= sqrt(2/(pi x)) (Watson 13.74: x (J0^2+Y0^2) < 2/pi), S = zeros with
(pi/2) c_g(Ma) T1 >= 1, K = |S| >= 8:
  | int_{T1}^inf sinc * dPhi/dM | <= Dmax * prod_S (2/(pi c_g(Ma) T1))^{1/2}
                                     * (pi c_{g_1}(Ma) T1/2)^{1/2} * T1^2/(K-5).

Usage:  indep_enclose.py --q 4 --Ma 3.6 --Mb 3.618         (cell)
        indep_enclose.py --q 4 --Ma 8.5 --Mb 8.5           (single point)
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from decimal import Decimal
from fractions import Fraction

from flint import acb, arb, arb_series, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))
PREC = 110
NMAX = 6
Z0 = Fraction(3, 2)
MARGIN = Fraction(13, 10)        # Rdom = MARGIN * T1
ZERO_RAD = Fraction(1, 10 ** 20)
E2_GOAL = 1e-12                  # steering only


# ----------------------------------------------------------------------------
def qarb(x: Fraction) -> arb:
    return arb(x.numerator) / arb(x.denominator)


def ball0(r: arb) -> arb:
    """Real ball [-u, u], u >= r."""
    u = r.upper()
    if not u.is_finite():
        raise ArithmeticError("non-finite radius")
    return arb(0, u)


def cball0(r: arb) -> acb:
    """Complex ball (square) containing the disc |w| <= r."""
    b = ball0(r)
    return acb(b, b)


def hull(a: arb, b: arb) -> arb:
    return a.union(b)


_Z = {}


def zeros(q: int) -> list[arb]:
    """Zero balls [dec - 1e-20, dec + 1e-20], parsed as exact rationals."""
    if q in _Z:
        return _Z[q]
    fn = os.path.join(ZDIR, {4: "chi4_zeros.txt", 3: "chi3_zeros.txt"}[q])
    fr = []
    with open(fn, encoding="ascii") as fh:
        for line in fh:
            line = line.strip()
            if line:
                fr.append(Fraction(Decimal(line)))
    for a, b in zip(fr, fr[1:]):
        if not b - a > 2 * ZERO_RAD:
            raise ArithmeticError("zero windows not disjoint/ordered")
    if not fr[0] > 1:
        raise ArithmeticError("first zero suspicious")
    old = ctx.prec
    ctx.prec = 160
    rad = qarb(ZERO_RAD).upper()
    out = [qarb(f) + arb(0, rad) for f in fr]
    ctx.prec = old
    _Z[q] = out
    return out


# ---- closed forms through Hurwitz-zeta / lgamma power series -----------------
def loglam_d12(q: int, s: arb) -> tuple[arb, arb]:
    """((log Lambda)'(s), (log Lambda)''(s)), Lambda(s) = (q/pi)^{(s+1)/2}
    Gamma((s+1)/2) L(s,chi), chi the odd real primitive character mod q,
    L(s,chi) = q^{-s} sum_a chi(a) zeta(s, a/q)."""
    x = arb_series([s, 1], prec=3)
    if q == 4:
        hz = x.zeta(arb(1) / 4) - x.zeta(arb(3) / 4)
    elif q == 3:
        hz = x.zeta(arb(1) / 3) - x.zeta(arb(2) / 3)
    else:
        raise ValueError
    # log L = -s log q + log(hz)
    lg = hz.log()
    half = (x + 1) / 2
    logLam = half * (arb(q) / arb.pi()).log() + half.lgamma() + lg - x * arb(q).log()
    d1 = logLam[1]
    d2 = 2 * logLam[2]
    if not (d1.is_finite() and d2.is_finite()):
        raise ArithmeticError("closed form not finite")
    return d1, d2


def sigma2_D(q: int, M: arb) -> tuple[arb, arb]:
    d1, d2 = loglam_d12(q, M + arb(1) / 2)
    return 4 * M * d1, 4 * d1 + 4 * M * d2


# ---- Rayleigh sums lambda_n = sum_k j_{0,k}^{-2n} as exact rationals -----------
def rayleigh(nmax: int) -> list[Fraction]:
    """lam[n], n = 1..nmax, from J0 = sum a_k w^k (w = z^2), a_k = (-1/4)^k/k!^2
    and  -log J0 = sum (lam_n/n) w^n :  lam_m = -m a_m - sum_{n<m} lam_n a_{m-n}."""
    a = [Fraction(1)]
    for k in range(1, nmax + 1):
        a.append(a[-1] * Fraction(-1, 4) / (k * k))
    lam = [None]
    for m in range(1, nmax + 1):
        v = -m * a[m] - sum(lam[n] * a[m - n] for n in range(1, m))
        lam.append(v)
    return lam


LAM = rayleigh(NMAX + 1)
assert LAM[1] == Fraction(1, 4) and LAM[2] == Fraction(1, 32) and LAM[3] == Fraction(1, 192)
assert LAM[4] == Fraction(11, 12288) and all(v > 0 for v in LAM[1:])


# ---- Bessel J0, J1 on complex balls ---------------------------------------------
def bessel01(z: acb) -> tuple[acb, acb]:
    """Enclosures of J0, J1 on the complex ball z.
    |J_n(w)| <= e^{|Im w|} (Bessel's integral), |J1(w)| <= |w|/2 e^{|Im w|},
    J0' = -J1, J1' = (J0-J2)/2, J0'' = (J2-J0)/2, J1'' = (J3-3J1)/4, so with
    E = sup e^{|Im w|}, r = radius, z0 = centre:
      |J0(w)-J0(z0)| <= min( r E min(1,|w|max/2),  |J1(z0)| r + r^2 E/2 ),
      |J1(w)-J1(z0)| <= min( r E,                  |J0(z0)-J2(z0)| r/2 + r^2 E/2 )."""
    re, im = z.real, z.imag
    rr, ri = arb(re.rad()), arb(im.rad())
    z0 = acb(arb(re.mid()), arb(im.mid()))
    j0 = z0.bessel_j(0)
    j1 = z0.bessel_j(1)
    if rr == 0 and ri == 0:
        return j0, j1
    r = (rr * rr + ri * ri).sqrt()
    E = im.abs_upper().exp()
    az = abs(z).upper()
    half = az / 2
    lip0 = r * E * (half if half < 1 else arb(1))
    tay0 = abs(j1) * r + r * r * E / 2
    e0 = lip0 if lip0.upper() < tay0.upper() else tay0
    j2 = z0.bessel_j(2)
    lip1 = r * E
    tay1 = abs(j0 - j2) / 2 * r + r * r * E / 2
    e1 = lip1 if lip1.upper() < tay1.upper() else tay1
    return j0 + cball0(e0), j1 + cball0(e1)


NONFINITE = acb(arb("nan"), arb("nan"))


class Cell:
    def __init__(self, q: int, Ma: Fraction, Mb: Fraction, e2_goal: float = E2_GOAL,
                 z0: Fraction = Z0, verbose: bool = False):
        if not (0 < Ma <= Mb):
            raise ValueError("need 0 < Ma <= Mb")
        ctx.prec = PREC
        self.q, self.Ma, self.Mb = q, Ma, Mb
        zs = zeros(q)
        N = len(zs)
        a, b = qarb(Ma), qarb(Mb)
        gN = zs[-1]
        if not gN.lower() ** 2 > 3 * b.upper() ** 2:
            raise ArithmeticError("g_N^2 > 3 Mb^2 fails")

        def c_of(g, M):
            return 4 * M / (M * M + g * g).sqrt()

        def cp_of(g, M):
            d = M * M + g * g
            return 4 * g * g / (d * d.sqrt())

        ca = [c_of(g, a) for g in zs]           # c increasing in M
        cb = [c_of(g, b) for g in zs]
        # ---------------- large-t remainder / choice of T1 -----------------
        s2a, Da = sigma2_D(q, a)
        s2b, Db = sigma2_D(q, b)
        if not (s2a.lower() > 0 and Da.lower() > 0):
            raise ArithmeticError("sigma^2, D not positive")
        # D(M) <= sum_listed c c' (cell) + D_T(Mb)
        hd_b = arb(0)
        hv_b = arb(0)
        for g, c in zip(zs, cb):
            hd_b += c * cp_of(g, b)
            hv_b += c * c / 2
        VT_b = s2b - hv_b                         # unlisted variance at Mb
        DT_b = Db - hd_b                          # unlisted D at Mb (D_T nondecreasing)
        if not (VT_b.lower() > 0 and DT_b.lower() > 0):
            raise ArithmeticError("unlisted tail sums not positive")
        cc_hi = arb(0)
        for g, x, y in zip(zs, ca, cb):
            # c in [c(Ma), c(Mb)];  c' decreasing in M:  c' <= c'(Ma)
            cc_hi += arb(y.upper()) * arb(cp_of(g, a).upper())
        Dmax = arb(cc_hi.upper()) + arb(DT_b.upper())
        self.Dmax = Dmax

        def remainder(T1: Fraction):
            T = qarb(T1)
            prod = arb(1)
            K = 0
            for c in ca:
                x = arb.pi() * c * T / 2
                if not x.lower() >= 1:
                    break
                prod *= 1 / x.sqrt()
                K += 1
            if K < 8:
                return None, K
            f1 = (arb.pi() * ca[0] * T / 2).sqrt()
            E2 = Dmax * prod * f1 * T * T / (K - 5)
            if not E2.is_finite():
                return None, K
            return arb(E2.upper()), K

        sig = math.sqrt(float(s2a.mid()))
        T1 = None
        x = max(2.0 / sig, 0.2)
        while x < 60:
            cand = Fraction(math.ceil(x * 64), 64)
            E2, K = remainder(cand)
            if E2 is not None and float(E2.mid()) <= e2_goal:
                T1 = cand
                break
            x *= 1.1
        if T1 is None:
            raise ArithmeticError("no T1 found")
        self.T1, self.E2, self.K = T1, E2, K
        # ---------------- head size N' ---------------------------------------
        z0b = qarb(z0)
        j0z = z0b.bessel_j(0)
        if not (z0b < arb("2.4") and j0z.lower() > 0):
            raise ArithmeticError("z0 must be < j_{0,1}")
        Rdom = qarb(MARGIN * T1)
        Np = None
        for i in range(N):
            if (cb[i] * Rdom).upper() <= z0b.lower():
                Np = i                      # zeros[i:] are 'rest'; cb decreasing in i
                break
        if Np is None or Np < 8:
            raise ArithmeticError("no admissible head size")
        if not all((cb[i] * Rdom).upper() <= z0b.lower() for i in (Np, Np + 1, N - 1)):
            raise ArithmeticError("domain check")
        if not zs[Np].lower() ** 2 > 3 * b.upper() ** 2:
            raise ArithmeticError("g_{N'+1}^2 > 3 Mb^2 fails")
        self.Np, self.Rdom = Np, arb(Rdom.lower())
        head = zs[:Np]
        rest = zs[Np:]
        # ---------------- exact rest sums n = 1 --------------------------------
        def Vp_Dp(M, s2, D):
            v = arb(0)
            d = arb(0)
            for g in head:
                c = c_of(g, M)
                v += c * c / 2
                d += c * cp_of(g, M)
            return s2 - v, D - d
        Vpa, Dpa = Vp_Dp(a, s2a, Da)
        Vpb, Dpb = Vp_Dp(b, s2b, Db)
        if not (Vpa.lower() > 0 and Dpa.lower() > 0):
            raise ArithmeticError("rest sums not positive")
        Vp = hull(Vpa, Vpb)                 # V' nondecreasing in M
        Dp = hull(Dpa, Dpb)                 # D' nondecreasing in M (g^2 > 3M^2 on rest)
        # ---------------- higher power sums ------------------------------------
        cG = arb(cb[-1].upper())            # >= c_g(M) for every unlisted g, M in cell
        S = [None, 2 * Vp]
        U = [None, Dp]
        for n in range(2, NMAX + 2):
            sR = arb(0)
            uR = arb(0)
            for j in range(Np, N):
                g = zs[j]
                c = hull(ca[j], cb[j])
                cp = hull(cp_of(g, a), cp_of(g, b))      # c' monotone (decreasing) in M
                p = c ** (2 * n - 1)
                sR += p * c
                uR += p * cp
            sU = cG ** (2 * n - 2) * 2 * VT_b
            uU = cG ** (2 * n - 2) * DT_b
            S.append(sR + hull(arb(0), arb(sU.upper())))
            U.append(uR + hull(arb(0), arb(uU.upper())))
        self.S, self.U = S, U
        # coefficients of the exponent and of B
        self.ecoef = [qarb(LAM[n] / n) * S[n] for n in range(1, NMAX + 1)]
        self.bcoef = [qarb(2 * LAM[n]) * U[n] for n in range(1, NMAX + 1)]
        j1z = z0b.bessel_j(1)
        rho1 = -j0z.log() - sum((qarb(LAM[n] / n) * z0b ** (2 * n) for n in range(1, NMAX + 1)), arb(0))
        rho2 = j1z / j0z - sum((qarb(2 * LAM[n]) * z0b ** (2 * n - 1) for n in range(1, NMAX + 1)), arb(0))
        if not (rho1.lower() > 0 and rho2.lower() > 0):
            raise ArithmeticError("series remainders not positive")
        self.r1 = arb((rho1 / z0b ** (2 * NMAX + 2) * arb(S[NMAX + 1].upper())).upper())
        self.r2 = arb((rho2 / z0b ** (2 * NMAX + 1) * arb(U[NMAX + 1].upper())).upper())
        # head coefficients over the cell: monotone hulls
        self.c = [hull(ca[j], cb[j]) for j in range(Np)]
        self.cp = [hull(cp_of(zs[j], a), cp_of(zs[j], b)) for j in range(Np)]
        self.info = dict(T1=float(T1), K=K, Np=Np, E2=float(E2.mid()), sig2=float(s2a.mid()),
                         Vp=(float(Vpa.mid()), float(Vpb.mid())), Dmax=float(Dmax.mid()))
        self.nev = 0

    # ---------------------------------------------------------------------------
    def integrand(self, t: acb, analytic=False) -> acb:
        self.nev += 1
        at = abs(t)
        if not (at.upper() <= self.Rdom):
            return NONFINITE
        n = self.Np
        J0 = [None] * n
        J1 = [None] * n
        c = self.c
        for i in range(n):
            J0[i], J1[i] = bessel01(c[i] * t)
        pre = [None] * (n + 1)
        p = acb(1)
        pre[0] = p
        for i in range(n):
            p = p * J0[i]
            pre[i + 1] = p
        P = p
        suf = acb(1)
        acc = acb(0)
        cp = self.cp
        for i in range(n - 1, -1, -1):
            acc += cp[i] * J1[i] * pre[i] * suf
            suf = suf * J0[i]
        dP = -t * acc
        t2 = t * t
        expo = acb(0)
        B = acb(0)
        tp = acb(1)
        for k in range(NMAX):
            tp = tp * t2
            expo += self.ecoef[k] * tp
            B += self.bcoef[k] * tp
        au = arb(at.upper())
        tail_pow = au ** (2 * NMAX + 2)
        expo += cball0(self.r1 * tail_pow)
        B += cball0(self.r2 * tail_pow)
        F = t.sinc() * (-expo).exp() * (dP - P * B)
        return F

    def enclose(self, rel_tol=1e-11, abs_tol=1e-14, eval_limit=400000, depth_limit=400):
        ctx.prec = PREC
        I = acb.integral(self.integrand, 0, qarb(self.T1), rel_tol=arb(rel_tol), abs_tol=arb(abs_tol),
                         eval_limit=eval_limit, depth_limit=depth_limit)
        if not I.is_finite():
            raise ArithmeticError("integral not finite")
        val = (I.real + ball0(self.E2)) / arb.pi()
        return val


def enclose_cell(q: int, Ma, Mb, **kw) -> dict:
    t0 = time.time()
    Ma, Mb = Fraction(str(Ma)), Fraction(str(Mb))
    cell = Cell(q, Ma, Mb, **kw)
    v = cell.enclose()
    lo = math.nextafter(float(v.lower()), -math.inf)
    hi = math.nextafter(float(v.upper()), math.inf)
    return dict(q=q, Ma=Ma, Mb=Mb, lo=lo, hi=hi, neg=bool(v.upper() < 0), ball=v,
                nev=cell.nev, secs=time.time() - t0, **cell.info)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--Ma", required=True)
    ap.add_argument("--Mb", required=True)
    a = ap.parse_args()
    r = enclose_cell(a.q, a.Ma, a.Mb)
    print(f"q={r['q']} [{a.Ma}, {a.Mb}]  delta' in [{r['lo']:.10e}, {r['hi']:.10e}]  negative={r['neg']}  "
          f"T1={r['T1']:.4f} K={r['K']} N'={r['Np']} E2={r['E2']:.1e} evals={r['nev']} secs={r['secs']:.1f}")


if __name__ == "__main__":
    main()
