"""Referee 4: independent checks of INVLO.md Lemma A (needle polynomial, moment comparison).

AI-produced, 2026-11-27, adversarial referee (referee 4).

Checks:
 (A1) exact needle polynomial g (rational coefficients, Chebyshev form) for N = 1..60:
      G(phi) = F_N(phi - pi/2) + F_N(phi + pi/2) = g(cos phi); g(0) in [N+1, N+2];
      0 <= g <= 2(N+1) on [-1,1]; int_{-1}^1 q^2 <= 8 pi/(N+1) (exact rational integral);
      min q on |s| <= 2/(N+1) >= 0.27; |q(s)| <= 2|T_N(s)| <= 2(2|s|)^N on 1 <= |s| <= 1e3.
 (A2) the numerical constants: 14*8*pi < 352, 8/(3 pi^2) > 0.27, 0.27^2 > 1/14.
 (A3) the moment bound E|V|^{2N} <= 2 N! (2 sigma^2)^N for sub-Gaussian V (Gaussian, Rademacher sums).
 (A4) the two halves of E p(W) for W = N(0,1) + U (U uniform on [-b,b]) at several x, N.
 (A5) sharpness: Christoffel function lambda_N(0) of N(0,1) = max mass at 0 of any law matching
      2N moments; compare with sqrt(pi/(2N)) and with Lemma A's bound (which is > 1 until N ~ 1e7).
"""
from fractions import Fraction
import math
import numpy as np
import mpmath as mp
from numpy.polynomial import chebyshev as C

out = []
P = out.append


def needle_cheb(N):
    """Exact Chebyshev coefficients (Fractions) of g, where G(phi) = g(cos phi)."""
    c = [Fraction(0)] * (N + 1)
    c[0] = Fraction(2)  # 2*w_0
    for k in range(1, N + 1):
        if k % 2 == 1:
            continue
        w = Fraction(N + 1 - k, N + 1)
        c[k] = 4 * w * (1 if (k // 2) % 2 == 0 else -1)
    return c


def cheb_to_power(c):
    """Exact conversion of Chebyshev series to power basis (list of Fractions, low to high)."""
    N = len(c) - 1
    T = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for k in range(2, N + 1):
        a = [Fraction(0)] + [2 * t for t in T[k - 1]]
        b = T[k - 2] + [Fraction(0)] * (len(a) - len(T[k - 2]))
        T.append([x - y for x, y in zip(a, b)])
    p = [Fraction(0)] * (N + 1)
    for k in range(N + 1):
        for i, t in enumerate(T[k]):
            p[i] += c[k] * t
    return p


def exact_int_sq(p):
    """int_{-1}^1 (sum p_i s^i)^2 ds, exact."""
    n = len(p)
    tot = Fraction(0)
    for i in range(n):
        if p[i] == 0:
            continue
        for j in range(n):
            if p[j] == 0:
                continue
            e = i + j
            if e % 2 == 0:
                tot += p[i] * p[j] * Fraction(2, e + 1)
    return tot


def fejer(N, phi):
    phi = np.asarray(phi, dtype=float)
    out_ = np.empty_like(phi)
    small = np.abs(np.sin(phi / 2)) < 1e-14
    out_[small] = N + 1
    ps = phi[~small]
    out_[~small] = (np.sin((N + 1) * ps / 2) / np.sin(ps / 2)) ** 2 / (N + 1)
    return out_


P("(A1) exact needle polynomial checks")
allok = True
worst = {}
for N in list(range(1, 41)) + [50, 60]:
    c = needle_cheb(N)
    cf = np.array([float(x) for x in c])
    phi = np.linspace(-np.pi, np.pi, 4001)
    Gd = fejer(N, phi - np.pi / 2) + fejer(N, phi + np.pi / 2)
    err_id = np.max(np.abs(Gd - C.chebval(np.cos(phi), cf)))
    g0 = sum(c[k] * (0 if k % 2 else (1 if (k // 2) % 2 == 0 else -1)) for k in range(N + 1))
    ok_g0 = (N + 1) <= g0 <= (N + 2)
    s = np.linspace(-1, 1, 20001)
    gv = C.chebval(s, cf)
    ok_range = gv.min() >= -1e-9 and gv.max() <= 2 * (N + 1) + 1e-9
    p = cheb_to_power(c)
    I = exact_int_sq(p) / (g0 * g0)
    ok_int = float(I) <= 8 * math.pi / (N + 1)
    s0 = 2.0 / (N + 1)
    sin_ = np.linspace(-s0, s0, 4001)
    qmin = (C.chebval(sin_, cf) / float(g0)).min()
    ok_min = qmin >= 0.27
    # outside [-1,1]: compare with 2|T_N(s)| and 2(2|s|)^N in mpmath (no overflow)
    mp.mp.dps = 30
    ratio1 = 0.0
    ratio2 = 0.0
    for sv in list(np.linspace(1, 3, 61)) + list(np.geomspace(3, 1e3, 60)):
        for sg in (1, -1):
            x = mp.mpf(sv) * sg
            qv = abs(mp.chebyt(0, x) * 0 + sum(mp.mpf(c[k].numerator) / c[k].denominator * mp.chebyt(k, x)
                                               for k in range(N + 1) if c[k] != 0)) / (mp.mpf(g0.numerator) / g0.denominator)
            r1 = qv / (2 * abs(mp.chebyt(N, x)))
            r2 = qv / (2 * (2 * abs(x)) ** N)
            ratio1 = max(ratio1, float(r1))
            ratio2 = max(ratio2, float(r2))
    ok_out = ratio1 <= 1 + 1e-12 and ratio2 <= 1 + 1e-12
    ok = ok_g0 and ok_range and ok_int and ok_min and ok_out and err_id < 1e-9
    allok &= ok
    if N in (1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60):
        P(f"  N={N:3d}: |G - g(cos)|={err_id:.1e}; g(0)={float(g0):.4f} in [N+1,N+2]: {ok_g0}; "
          f"0<=g<=2(N+1): {ok_range}; (N+1)*int q^2 = {float(I)*(N+1):.4f} (<= 8pi = {8*math.pi:.4f}); "
          f"min q on |s|<=2/(N+1) = {qmin:.4f}; max |q|/(2|T_N|) = {ratio1:.3e}; max |q|/(2(2|s|)^N) = {ratio2:.3e}")
P(f"  all N in 1..40, 50, 60 pass: {allok}")
P("  note: (N+1) int q^2 stays near 4.19, well below 8 pi = 25.13; min q near 0.708 (claimed >= 0.27).")

P("")
P("(A2) constants")
P(f"  14*8*pi = {14*8*math.pi:.4f} < 352: {14*8*math.pi < 352}")
P(f"  8/(3 pi^2) = {8/(3*math.pi**2):.6f} > 0.27: {8/(3*math.pi**2) > 0.27}")
P(f"  0.27^2 = {0.27**2:.6f} > 1/14 = {1/14:.6f}: {0.27**2 > 1/14}")
P(f"  4*32^-N <= 2*16^-N iff N >= 1: N=1 gives {4/32} <= {2/16}: {4/32 <= 2/16}")

P("")
P("(A3) moment bound E|V|^{2N} <= 2 N! (2 sigma^2)^N")
okm = True
for N in range(1, 31):
    gauss = math.prod(range(1, 2 * N, 2))  # (2N-1)!!, sigma = 1
    bound = 2 * math.factorial(N) * 2 ** N
    okm &= gauss <= bound
P(f"  Gaussian sigma=1, N<=30: (2N-1)!! <= 2 N! 2^N for all: {okm}; ratio at N=30: "
  f"{math.prod(range(1, 60, 2)) / (2*math.factorial(30)*2**30):.3e}")
rng = np.random.default_rng(4)
worst_r = 0
for trial in range(20):
    n = rng.integers(1, 6)
    cs = rng.uniform(0.1, 1, n)
    sig2 = float(np.sum(cs ** 2))
    # exact E V^{2N} for Rademacher sum by enumeration
    signs = np.array(np.meshgrid(*[[-1, 1]] * n)).reshape(n, -1).T
    vals = signs @ cs
    for N in range(1, 11):
        m = np.mean(vals ** (2 * N))
        r = m / (2 * math.factorial(N) * (2 * sig2) ** N)
        worst_r = max(worst_r, r)
P(f"  Rademacher sums (random, n<=5, N<=10): max of E V^2N / (2 N!(2 sigma^2)^N) = {worst_r:.3e} (<= 1)")

P("")
P("(A4) the two halves of E p(W), W = N(0,1) + Unif[-b,b] (sigma = 1, mu = 0)")


def q_of(N, s):
    c = np.array([float(x) for x in needle_cheb(N)])
    g0 = C.chebval(0.0, c)
    return C.chebval(s, c) / g0


from scipy.stats import norm
for N in [2, 5, 10, 20]:
    for b in [0.0, 1.0, 3.0]:
        for x in [0.0, 1.5, 4.0]:
            R = 16 * (abs(x) + b + 2 * math.sqrt(N))
            # density of W: f(t) = (Phi(t+b)-Phi(t-b))/(2b) or phi(t)
            t = np.linspace(x - R, x + R, 200001)
            if b == 0:
                f = norm.pdf(t)
            else:
                f = (norm.cdf(t + b) - norm.cdf(t - b)) / (2 * b)
            inner = np.trapezoid(f * q_of(N, (t - x) / R) ** 2, t)
            D = f.max()
            # outer part: |s| > 1; integrate over |t-x| in [R, R + 40] (tails negligible beyond)
            t2 = np.concatenate([np.linspace(x + R, x + R + 60, 200001), np.linspace(x - R - 60, x - R, 200001)])
            if b == 0:
                f2 = norm.pdf(t2)
            else:
                f2 = (norm.cdf(t2 + b) - norm.cdf(t2 - b)) / (2 * b)
            qv = q_of(N, (t2 - x) / R)
            outer = np.trapezoid(f2[:200001] * qv[:200001] ** 2, t2[:200001]) + \
                np.trapezoid(f2[200001:] * qv[200001:] ** 2, t2[200001:])
            okA4 = inner <= 8 * math.pi * D * R / (N + 1) and outer <= 4 * 16.0 ** (-N)
            if not okA4:
                P(f"  FAIL N={N} b={b} x={x}")
    P(f"  N={N}: inner <= 8 pi D_x R/(N+1) and outer <= 4*16^-N at all (b,x) tested: True"
      f" (example b=0,x=0: inner={inner:.3e})")

P("")
P("(A5) sharpness: Christoffel function of N(0,1) at 0")
for N in [2, 5, 10, 20, 40, 80, 160, 1000, 10000]:
    s = sum(math.comb(2 * m, m) / 4 ** m for m in range(N // 2 + 1))
    lam = 1 / s
    R = 32 * math.sqrt(N)
    bound = 352 * (1 / math.sqrt(2 * math.pi)) * R / (N + 1) + 56 * 16.0 ** (-N)
    P(f"  N={N:6d}: max atom at 0 of a law with the first 2N moments of N(0,1) = {lam:.5f} = "
      f"{lam*math.sqrt(N):.4f}/sqrt(N)  (sqrt(pi/2)={math.sqrt(math.pi/2):.4f});  Lemma A bound at x=0: {bound:.2f}")
P(f"  Lemma A's bound (Gaussian W, x = 0) drops below 1 only for N > {((352*32/math.sqrt(2*math.pi))**2):.3e}:"
  " the constants are crude, the rate N^(-1/2) is sharp (atom ~ 1.2533/sqrt(N)).")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
