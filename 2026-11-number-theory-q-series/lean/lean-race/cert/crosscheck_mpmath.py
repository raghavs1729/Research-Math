#!/usr/bin/env python3
"""Independent, NON-RIGOROUS cross-check of cert/window_cert.py with mpmath (floating point).

No formula is shared with the certified script except the definition of the window:
  * the window is evaluated from its definition f(u) = S(20(27/50 - u)) * S(20(27/50 + u));
  * its transform is  fhat(r) = 2 [ sin(r b)/r + Int_b^h f(u) cos(r u) du ]  with the ramp integral
    done NUMERICALLY (96-point Gauss-Legendre on [b, h]; the ramp spans at most 6.4 wavelengths
    for r <= R2), not by the closed form and not by the moment series;
  * the weight uses mpmath's digamma;
  * the outer integral is mpmath's adaptive Gauss-Legendre on a partition of [0, R2] with step 1/2;
  * the tail beyond R2 = 800 is ignored (window_cert.py's rigorous bound for it is < 6e-7).

Run:  python3 cert/crosscheck_mpmath.py > cert/crosscheck_mpmath_output.txt
"""
from mpmath import mp, mpf, quad, cos, sin, digamma, log, pi, re, nstr
from mpmath.calculus.quadrature import GaussLegendre

mp.dps = 25
H = mpf(27) / 50
B = mpf(49) / 100
R2 = 800


def S(t):
    if t <= 0:
        return mpf(0)
    if t >= 1:
        return mpf(1)
    return 10 * t**3 - 15 * t**4 + 6 * t**5


def f(u):
    return S(20 * (H - u)) * S(20 * (H + u))


# 96-point Gauss-Legendre nodes on the ramp [B, H], with the window values precomputed
_nodes = GaussLegendre(mp).calc_nodes(6, mp.prec)          # nodes and weights on [-1, 1]
RAMP_NODES = [((B + H) / 2 + (H - B) / 2 * x, (H - B) / 2 * w) for x, w in _nodes]
RAMP_VALUES = [(u, w * f(u)) for u, w in RAMP_NODES]


def fhat(r):
    ramp = sum(wf * cos(r * u) for u, wf in RAMP_VALUES)
    plateau = sin(r * B) / r if r != 0 else B
    return 2 * (plateau + ramp)


def weight(r):
    return re(digamma(mpf(3) / 4 + 1j * r / 2)) + log(4 / pi)


def main():
    int_f = 2 * (B + sum(wf for _, wf in RAMP_VALUES))
    # self-test of the ramp quadrature against mpmath's adaptive quadrature at a few r
    for r in [mpf(3), mpf(111), mpf(799)]:
        direct = 2 * (quad(lambda u: cos(r * u), [0, B]) + quad(lambda u: f(u) * cos(r * u), [B, H]))
        assert abs(direct - fhat(r)) < mpf(10) ** -18, (r, direct, fhat(r))
    pts = [mpf(k) / 2 for k in range(0, 2 * R2 + 1)]
    gamma0 = quad(lambda r: fhat(r) ** 2 * weight(r), pts) / pi
    lo, hi = mpf("0.16116795432"), mpf("0.16116796545")
    print("crosscheck_mpmath.py (floating point, not a proof)")
    print(f"Integral f                 ~ {nstr(int_f, 15)}      (exact value 103/100)")
    print(f"(1/pi) Int_0^{R2} fhat^2 W   ~ {nstr(gamma0, 12)}")
    print(f"window_cert.py certifies chi4Gamma0 f in [{lo}, {hi}]")
    print(f"difference from the certified main part 0.16116795988 ~ {nstr(gamma0 - mpf('0.16116795988'), 3)}")
    print("agreement to 1e-6:", abs(gamma0 - mpf("0.16116795988")) < mpf(10) ** -6)


if __name__ == "__main__":
    main()
