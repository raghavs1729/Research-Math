"""Referee 4: independent checks of INVLO.md Lemma D (densities of random cosine sums).

AI-produced, 2026-11-27, adversarial referee (referee 4).

 (D1) |J_0(y)| <= (2/(pi y))^{1/2}: y (J_0^2 + Y_0^2)(y) increases to 2/pi (Nicholson); sup ratio.
 (D2) 0 < J_0(y) <= exp(-y^2/4) on [0, j_{0,1}).
 (D3) |J_0(y)| < pi^{-1/2} for y > 2.
 (D4) int_R min(1, (2/(pi|xi|a))^{3/2}) dxi = 12/(pi a).
 (D5) Lemma D(a): sup density of sums with three coefficients >= a versus 6/(pi^2 a).
 (D6) Lemma D(b): pointwise factor bounds used in the proof, and sup density versus the bound,
      on adversarial configurations (#I = 4, b/a large; many terms; mixed J).
"""
import math
import numpy as np
import mpmath as mp
from scipy.special import j0, y0, jn_zeros

out = []
P = out.append

P("(D1) Nicholson monotonicity and |J_0| <= sqrt(2/(pi y))")
y = np.concatenate([np.geomspace(1e-6, 1, 20000), np.linspace(1, 5000, 2000001)])
h = y * (j0(y) ** 2 + y0(y) ** 2)
P(f"  y(J0^2+Y0^2): min diff of consecutive values = {np.min(np.diff(h)):.3e} (>= 0 means increasing);"
  f" max = {h.max():.12f}, 2/pi = {2/np.pi:.12f}")
r = np.abs(j0(y)) / np.sqrt(2 / (np.pi * y))
P(f"  sup |J0(y)|/sqrt(2/(pi y)) on (0,5000] = {r.max():.9f} (<= 1)")
mp.mp.dps = 40
vals = []
for Y in [10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6]:
    Ym = mp.mpf(Y)
    vals.append((Y, float(Ym * (mp.besselj(0, Ym) ** 2 + mp.bessely(0, Ym) ** 2) - 2 / mp.pi)))
P("  mpmath: y(J0^2+Y0^2) - 2/pi at y = 10..1e6: " + ", ".join(f"{a}: {b:.3e}" for a, b in vals) + " (all < 0)")

P("")
P("(D2) 0 < J_0(y) <= exp(-y^2/4) on [0, j01)")
j01 = jn_zeros(0, 1)[0]
yy = np.linspace(0, j01, 2000001)[:-1]
gap = np.exp(-yy ** 2 / 4) - j0(yy)
P(f"  j01 = {j01:.9f}; min(exp(-y^2/4) - J0) = {gap.min():.3e}; min J0 = {j0(yy).min():.3e} (> 0)")
mpgap = min(float(mp.e ** (-mp.mpf(t) ** 2 / 4) - mp.besselj(0, mp.mpf(t))) for t in np.linspace(0.001, float(j01) - 1e-6, 400))
P(f"  mpmath check on 400 points in (0, j01): min gap = {mpgap:.3e} (>= 0)")
P("  note: exp(-y^2/4) - J0(y) = y^4/32 - y^4/64 + O(y^6) = y^4/64 + O(y^6) near 0 (Taylor);"
  " the float value -2.2e-16 is rounding at y = 0, where the gap is exactly 0.")

P("")
P("(D3) |J_0(y)| < pi^{-1/2} for y > 2")
y3 = np.linspace(2, 2000, 4000001)
P(f"  max |J0| on [2,2000] = {np.abs(j0(y3)).max():.6f} < pi^-1/2 = {np.pi**-0.5:.6f}")

P("")
P("(D4) integral of min(1,(2/(pi |xi| a))^{3/2})")
for a in [0.1, 1.0, 7.0]:
    x0 = 2 / (np.pi * a)
    val = 2 * (x0 + 2 * x0)  # closed form: 2[x0 + x0^{3/2} * 2 x0^{-1/2}]
    from scipy.integrate import quad
    num = 2 * (quad(lambda t: 1.0, 0, x0)[0] + quad(lambda t: (x0 / t) ** 1.5, x0, np.inf)[0])
    P(f"  a={a}: numeric {num:.10f}, claimed 12/(pi a) = {12/(np.pi*a):.10f}")


def density_sup(coeffs, xmax=None, nxi=400001, ximax_factor=4000.0, nx=241):
    """sup_x of the density of sum c_k cos(theta_k) by Fourier inversion (trapezoid)."""
    coeffs = np.asarray(coeffs, dtype=float)
    a = coeffs.min()
    ximax = ximax_factor / a
    xi = np.linspace(0, ximax, nxi)
    phi = np.ones_like(xi)
    for c in coeffs:
        phi *= j0(xi * c)
    if xmax is None:
        xmax = coeffs.sum()
    xs = np.linspace(0, xmax, nx)
    best = 0.0
    bx = 0.0
    for i in range(0, nx, 20):
        X = xs[i:i + 20]
        f = np.trapezoid(phi[None, :] * np.cos(xi[None, :] * X[:, None]), xi, axis=1) / np.pi
        j = np.argmax(f)
        if f[j] > best:
            best, bx = f[j], X[j]
    tail = (2 / (np.pi * a)) ** 1.5 * 2 / np.sqrt(ximax) / np.pi if (coeffs >= a).sum() >= 3 else np.nan
    return best, bx, tail


P("")
P("(D5) Lemma D(a): sup density vs 6/(pi^2 a)")
for coeffs in [[1, 1, 1], [1, 1, 1, 0.3], [1, 1.5, 2], [1, 1, 1, 1], [1, 1, 1, 0.05, 0.05, 0.05]]:
    best, bx, tail = density_sup(coeffs)
    a = 1.0
    P(f"  coeffs {coeffs}: sup density ~ {best:.4f} at x ~ {bx:.3f} (truncation <= {tail:.1e});"
      f" bound 6/(pi^2 a) = {6/np.pi**2/a:.4f}: ok {best <= 6/np.pi**2/a}")

P("")
P("(D6) Lemma D(b)")
rng = np.random.default_rng(7)


def bound_b(cI, a, b):
    th = max(np.exp(-a * a / (b * b)), np.pi ** -0.5)
    S = np.sum(np.asarray(cI) ** 2)
    return (1 / (2 * np.pi)) * (2 * np.sqrt(np.pi) / np.sqrt(S) + 12 * th ** (len(cI) - 3) / (np.pi * a)), th


# pointwise factor bound for |xi| > 2/b
worst = 0.0
for trial in range(200):
    n = rng.integers(4, 30)
    a = rng.uniform(0.01, 1)
    b = a * rng.uniform(1, 50)
    cI = rng.uniform(a, b, n)
    th = max(np.exp(-a * a / (b * b)), np.pi ** -0.5)
    xi = np.linspace(2 / b * (1 + 1e-9), 2 / b + 200 / a, 20001)
    m = max(np.abs(j0(xi * c)).max() for c in cI)
    worst = max(worst, m / th)
P(f"  max over random (a,b,I) of max_{{|xi|>2/b}} |J0(xi c)|/theta = {worst:.6f} (<= 1)")

cases = [
    ("#I=4, b/a=1", [1, 1, 1, 1], 1, 1),
    ("#I=4, b/a=20", [0.05, 0.05, 0.05, 1.0], 0.05, 1.0),
    ("#I=5, b/a=3", [0.3, 0.5, 0.7, 0.9, 0.9], 0.3, 0.9),
    ("#I=20, b/a=4", list(np.linspace(0.05, 0.2, 20)), 0.05, 0.2),
    ("#I=60, b/a=10", list(np.linspace(0.01, 0.1, 60)), 0.01, 0.1),
]
for name, cI, a, b in cases:
    best, bx, tail = density_sup(cI, nxi=600001, ximax_factor=3000.0)
    bd, th = bound_b(cI, a, b)
    P(f"  {name}: sup density ~ {best:.4f} (at x~{bx:.3f}); D(b) bound = {bd:.4f} (theta={th:.3f}); ok {best <= bd}")

# a case where J is larger than I (extra small coefficients) -- bound unchanged, density smaller
cI = list(np.linspace(0.05, 0.2, 12))
cJ = cI + [0.01] * 5 + [0.3]
best, bx, _ = density_sup(cJ, nxi=600001, ximax_factor=3000.0)
bd, th = bound_b(cI, 0.05, 0.2)
P(f"  J = I + 6 extra terms: sup density ~ {best:.4f}; bound from I alone {bd:.4f}; ok {best <= bd}")
P("  remark: the first term (2pi)^-1 2 sqrt(pi) Sigma_I^-1/2 = 1/sqrt(pi Sigma_I) is exactly the Gaussian")
P("  density at 0 of a sum with variance Sigma_I/2, so the bound is essentially sharp for many equal terms.")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
