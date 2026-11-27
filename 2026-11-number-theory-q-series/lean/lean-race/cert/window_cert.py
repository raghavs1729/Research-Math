#!/usr/bin/env python3
"""Certified upper bound for chi4Gamma0 of the concrete race window (Arb ball arithmetic).

THE WINDOW (all parameters rational; the Lean definition is `raceWindow` in RaceChallengeDeps.lean)

    S(t)  = 0                         (t <= 0)
          = 10 t^3 - 15 t^4 + 6 t^5   (0 <= t <= 1)        "quintic smoothstep", C^2 on R
          = 1                         (t >= 1)
    f(u)  = S(20 (27/50 - u)) * S(20 (27/50 + u))

i.e. with b = 49/100, d = 1/20, h = b + d = 27/50:  f is even, f = 1 on [-b, b],
f(u) = S((h - |u|)/d) on b <= |u| <= h, f = 0 for |u| >= h.  f is C^2 with support [-h, h],
2h = 27/25 < log 3.

THE CLAIM that this script certifies

    chi4Gamma0 f := (1/(2 pi)) * Integral_{-oo}^{oo} |fhat(r)|^2 * W(r) dr  <  GAMMA_BOUND = 81/500,
    fhat(r) = Integral f(u) e^{iru} du,    W(r) = Re psi(3/4 + i r/2) + log(4/pi),  psi = Gamma'/Gamma,

together with the exact rational value  Integral f = 2b + d = 103/100,  so that

    chi4Gamma0 f  <  (1/2) (Integral f)^2 - MARGIN,     MARGIN = 10609/20000 - 81/500 = 7369/20000.

HOW (every number on the certified path is an Arb ball or an exact rational; no floating point)

 1. fhat is real and even, W is even, so chi4Gamma0 f = (1/pi) * Integral_0^oo fhat(r)^2 W(r) dr.
 2. [0, R], R = 2000: rigorous quadrature with Arb's `acb.integral` (Petras' adaptive
    Gauss-Legendre with rigorous error bounds from complex-analytic enclosures).  The integrand is
    evaluated on complex balls as  fhat(r)^2 * ( (psi(3/4 + i r/2) + psi(3/4 - i r/2))/2 + log(4/pi) ),
    which is holomorphic for |Im r| < 3/2 and equals fhat(r)^2 W(r) on the real axis.  fhat(r):
      * at points with |r| < 40: the moment series  fhat(r) = sum_k (-1)^k m_{2k} r^{2k} / (2k)!,
        m_{2k} = Integral f(u) u^{2k} du computed EXACTLY in rational arithmetic, K = 150 terms,
        with the rigorous tail bound  |m_{2k}| <= 2 h^{2k+1}  (since |f| <= 1 on [-h, h]);
      * for Re r > 10: the closed form obtained by three integrations by parts
        (the second derivative of f is continuous, the third is piecewise polynomial),  s = r d,
            fhat(r) = (2/(r^3 d^2)) * [ (60/s - 720/s^3) (cos(r h) - cos(r b))
                                         - (360/s^2) (sin(r h) + sin(r b)) ];
        the two representations are checked against each other on the overlap;
      * on wide complex balls (which the quadrature uses for its error bounds) fhat and psi are
        enclosed in mean-value form, value at the midpoint + radius * (bound for the derivative),
        with |fhat'| <= 2 h^2 e^{h |Im r|} and |psi'(w)| <= 1/x^2 + 1/x for Re w >= x > 0.
 3. Tail r > R:  |fhat(r)| <= M3 / r^3 with M3 = L^1 norm of the third derivative of f
    = (2/d^2) * 40 sqrt(3)/3  (three integrations by parts; Integral_0^1 |S'''(t)| dt = 40 sqrt(3)/3),
    and for Re z = 3/4
        Re psi(z) <= log|z| + 4/27
    (Binet's second formula, DLMF 5.9.15:  psi(z) = log z - 1/(2z) - 2 Int_0^oo t dt/((t^2+z^2)(e^{2 pi t}-1));
     Re(1/(2z)) > 0, |t^2 + z^2| = |z - it| |z + it| >= (Re z)^2 = 9/16, Int_0^oo t dt/(e^{2 pi t}-1) = 1/24),
    and |3/4 + i r/2| <= r for r >= 1.  Hence, with c0 = 4/27 + log(4/pi) and log R + c0 > 0,
        (1/pi) Int_R^oo fhat^2 W  <=  (M3^2/pi) * ( (log R + c0)/(5 R^5) + 1/(25 R^5) ).
    (The LOWER bound printed for information uses, from the same formula,
     Re psi(z) >= log|z| - 1/(2|z|) - 4/27, hence W(r) > 0 and a nonnegative tail for r >= R.)

Run:  python3 cert/window_cert.py > cert/window_cert_output.txt     (python-flint 0.9.0)
The output is deterministic.  Exit status 0 iff the inequality is certified.
"""
import sys
from fractions import Fraction
from math import comb, factorial

import flint
from flint import acb, arb, ctx

ctx.prec = 300

# ----------------------------------------------------------------------------------------------
# parameters (exact rationals)
B = Fraction(49, 100)          # plateau half-width
D = Fraction(1, 20)            # ramp width
H = B + D                      # half-support, 27/50
R = 2000                       # quadrature range [0, R]
R_SPLIT = 20                   # the range is split as [0, R_SPLIT] + [R_SPLIT, R]
K = 150                        # moment series: terms k = 0..K
GAMMA_BOUND = Fraction(81, 500)
RAMP = {0: Fraction(1), 3: Fraction(-10), 4: Fraction(15), 5: Fraction(-6)}   # 1 - S(t) on [0,1]


def to_arb(q):
    q = Fraction(q)
    return arb(q.numerator) / arb(q.denominator)


# ----------------------------------------------------------------------------------------------
# exact moments m_{2k} = Integral f(u) u^{2k} du = 2 [ b^{2k+1}/(2k+1) + d * Int_0^1 ramp(t) (b + d t)^{2k} dt ]
def moment(n):
    total = B ** (n + 1) / (n + 1)
    ramp_part = Fraction(0)
    for i in range(n + 1):
        cb = comb(n, i) * B ** (n - i) * D ** i
        ramp_part += cb * sum(c / (i + j + 1) for j, c in RAMP.items())
    return 2 * (total + D * ramp_part)


MOMENTS = [moment(2 * k) for k in range(K + 1)]
INT_F = MOMENTS[0]
assert INT_F == 2 * B + D == Fraction(103, 100)
# sanity of the rigorous bound |m_{2k}| <= 2 h^{2k+1} used for the series tail
assert all(0 < MOMENTS[k] <= 2 * H ** (2 * k + 1) for k in range(K + 1))

SERIES_COEFF = [to_arb(Fraction((-1) ** k) * MOMENTS[k] / factorial(2 * k)) for k in range(K + 1)]
b_, d_, h_ = to_arb(B), to_arb(D), to_arb(H)
UNIT = acb(arb(0, 1), arb(0, 1))       # the complex ball {|re| <= 1, |im| <= 1}
LOG4PI = (arb(4) / arb.pi()).log()
INDET = acb(arb(0, float("inf")), arb(0, float("inf")))     # the indeterminate ball
assert not INDET.is_finite()


def split(z):
    """(midpoint as an exact point, upper bound for the distance of any point of the ball z to it)."""
    m = z.mid()
    return m, (z - m).abs_upper()


def fhat_series_point(r):
    """Moment series with rigorous tail (meant for narrow balls r, |r| < 40)."""
    zz = r * r
    acc = acb(0)
    for c in reversed(SERIES_COEFF):
        acc = acc * zz + c
    x = (r.abs_upper() * h_) ** 2                        # upper bound for (|r| h)^2
    q = x / ((2 * K + 3) * (2 * K + 4))
    assert q < 1
    tail = 2 * h_ * x ** (K + 1) / arb(factorial(2 * K + 2)) / (1 - q)
    return acc + tail * UNIT


def fhat_closed(r):
    """Closed form; holomorphic for r != 0, used only for Re r > 10."""
    s = r * d_
    y = 60 / s - 720 / s ** 3
    return 2 / (r ** 3 * d_ ** 2) * (y * ((r * h_).cos() - (r * b_).cos())
                                     - 360 / s ** 2 * ((r * h_).sin() + (r * b_).sin()))


def fhat(r):
    """Enclosure of fhat on the complex ball r (fhat is entire)."""
    if r.real > 30:
        return fhat_closed(r)                 # interval evaluation of the closed form
    # mean-value form: fhat(r) in fhat(m) + (r - m) * sup |fhat'|, and
    # |fhat'(xi)| = |Int i u f(u) e^{i xi u} du| <= h * 2h * exp(h |Im xi|)   (|f| <= 1 on [-h, h])
    m, rad = split(r)
    mid_val = fhat_series_point(m) if m.abs_upper() < 40 else fhat_closed(m)
    return mid_val + UNIT * (rad * 2 * h_ ** 2 * (h_ * r.imag.abs_upper()).exp())


def digamma_ball(z):
    """Enclosure of psi on the complex ball z with Re z > 1/5, in mean-value form:
    |psi'(w)| <= sum_n 1/|w+n|^2 <= 1/x^2 + 1/x for Re w >= x > 0."""
    x = z.real.lower()
    if not (x > arb(1) / 5):
        return INDET
    m, rad = split(z)
    return m.digamma() + UNIT * (rad * (1 / x ** 2 + 1 / x))


def weight(r):
    """(psi(3/4 + i r/2) + psi(3/4 - i r/2))/2 + log(4/pi); holomorphic for |Im r| < 3/2."""
    if not (r.imag.abs_upper() < 1):
        return INDET
    half_ir = acb(0, 1) * r / 2
    return (digamma_ball(arb(3) / 4 + half_ir) + digamma_ball(arb(3) / 4 - half_ir)) / 2 + LOG4PI


def integrand(r, _analytic):
    w = weight(r)
    if not w.is_finite():
        return INDET
    return fhat(r) ** 2 * w


def show(name, x, digits=30):
    print(f"{name} = {x.str(digits)}")


def main():
    print("window_cert.py: certified bound for chi4Gamma0 of the race window")
    print(f"python-flint {flint.__version__}, working precision {ctx.prec} bits")
    print(f"b = {B}, d = {D}, h = {H}; R = {R}, R_SPLIT = {R_SPLIT}, K = {K}")
    print(f"Integral f = {INT_F} (exact);  (1/2)(Integral f)^2 = {INT_F ** 2 / 2} (exact)")

    # --- consistency checks (not part of the bound, but they must pass) -----------------------
    for rr in ["10.5", "17", "25.25", "39"]:
        r = acb(arb(rr))
        a, c = fhat_series_point(r).real, fhat_closed(r).real
        assert a.overlaps(c) and (a - c).abs_upper() < arb("1e-40"), (rr, a, c)
    print("check: moment series and closed form agree on [10, 40] (4 sample points, |diff| < 1e-40)")
    m3 = 2 / d_ ** 2 * 40 * arb(3).sqrt() / 3
    for rr in ["50", "333.3", "1999"]:
        r = arb(rr)
        assert fhat_closed(acb(r)).real.abs_upper() < m3 / r ** 3
        z = acb(arb(3) / 4, r / 2)
        assert z.digamma().real < z.abs_upper().log() + arb(4) / 27
    print("check: |fhat(r)| <= M3/r^3 and Re psi(z) <= log|z| + 4/27 hold at 3 sample points")

    # --- main part: (1/pi) Int_0^R fhat^2 W ---------------------------------------------------
    opts = dict(rel_tol=arb("1e-20"), abs_tol=arb("1e-20"), eval_limit=10 ** 7, depth_limit=10 ** 5)
    i1 = acb.integral(integrand, arb(0), arb(R_SPLIT), **opts)
    i2 = acb.integral(integrand, arb(R_SPLIT), arb(R), **opts)
    assert i1.is_finite() and i2.is_finite()
    assert i1.imag.contains(0) and i2.imag.contains(0)
    main_part = (i1.real + i2.real) / arb.pi()
    show("(1/pi) Int_0^20   fhat^2 W", i1.real / arb.pi())
    show("(1/pi) Int_20^R   fhat^2 W", i2.real / arb.pi())
    show("main part  (1/pi) Int_0^R ", main_part)
    assert main_part.rad() < arb("1e-12")

    # --- tail ------------------------------------------------------------------------------------
    c0 = arb(4) / 27 + LOG4PI
    logr = arb(R).log()
    assert logr + c0 > 0
    tail = m3 ** 2 / arb.pi() * ((logr + c0) / (5 * arb(R) ** 5) + 1 / (25 * arb(R) ** 5))
    show("M3 (L1 norm of the third derivative of f)", m3)
    show("tail bound (1/pi) Int_R^oo", tail, 10)

    upper = main_part + tail
    lower = main_part - tail
    show("certified upper bound for chi4Gamma0 f", upper, 20)
    show("certified lower bound for chi4Gamma0 f", lower, 20)

    bound = to_arb(GAMMA_BOUND)
    ok = bool(upper < bound)
    half_sq = INT_F ** 2 / 2
    margin = half_sq - GAMMA_BOUND
    print(f"GAMMA_BOUND = {GAMMA_BOUND};  MARGIN = (1/2)(Integral f)^2 - GAMMA_BOUND = {margin}")
    if ok:
        print(f"CERTIFIED: chi4Gamma0 f < {GAMMA_BOUND}")
        print(f"CERTIFIED: chi4Gamma0 f < (1/2) * (Integral f)^2 - {margin},   Integral f = {INT_F}")
        return 0
    print("NOT CERTIFIED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
