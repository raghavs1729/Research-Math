#!/usr/bin/env python3
"""NON-CERTIFIED floating-point sanity check of the numerical hypothesis `hcert` of
`grh_chi4_smoothed_race_real`:

    chi4Gamma0 f < (1/2) (∫ f)^2 - margin,
    chi4Gamma0 f = (1/2π) ∫ |f^(r)|^2 [Re ψ(3/4 + i r/2) + log(4/π)] dr,   f^(r) = ∫ f(u) e^{iru} du,

for a C^2 window f supported in [-a/2, a/2], a < log 3: a plateau of half-width b with quintic
"smoothstep" ramps (1 - 10t^3 + 15t^4 - 6t^5) of width d, half-support a/2 = b + d.

This script proves nothing.  It only shows that the hypothesis is plausible, with room to spare,
for a concrete admissible window; a certified (interval-arithmetic) evaluation is what the Lean
theorem asks for.  Usage:  python3 scripts/numeric_sanity.py
"""
from mpmath import mp, mpf, log, pi, quad, cos, digamma, re

mp.dps = 20
half = mpf("0.549")          # a/2;  a = 1.098 < log 3 = 1.09861...
assert 2 * half < log(3)


def ramp(t):
    return 1 - (10 * t**3 - 15 * t**4 + 6 * t**5)


def report(d):
    b = half - d
    f = lambda u: mpf(1) if abs(u) <= b else (ramp((abs(u) - b) / d) if abs(u) < half else mpf(0))
    fh = lambda r: 2 * (quad(lambda u: cos(r * u), [0, b]) + quad(lambda u: f(u) * cos(r * u), [b, half]))
    weight = lambda r: re(digamma(mpf(3) / 4 + 1j * r / 2)) + log(4 / pi)
    pts = [0, 5, 10, 20, 40, 80, 160, 320, 640]
    gamma0 = 2 * quad(lambda r: fh(r) ** 2 * weight(r), pts) / (2 * pi)
    int_f = 2 * (b + quad(lambda u: f(u), [b, half]))
    half_g0 = int_f**2 / 2
    print(f"ramp width d = {mp.nstr(d, 3)}:  Gamma_0 ~ {mp.nstr(gamma0, 6)},  (1/2)(int f)^2 ~ "
          f"{mp.nstr(half_g0, 6)},  margin available ~ {mp.nstr(half_g0 - gamma0, 4)}")


for d in ["0.1", "0.2", "0.3"]:
    report(mpf(d))
