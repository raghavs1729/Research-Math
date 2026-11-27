"""Numerical checks for INVLO.md (moment method: short relations and atoms).

AI-produced (Anthropic's Claude), 2026-11-27.

(1) The needle polynomial q_N of Lemma A: q(0) = 1, 0 <= q <= 2 on [-1,1],
    int_{-1}^{1} q^2 <= 8 pi/(N+1), q >= 0.27 on |s| <= 2/(N+1), and
    |q(s)| <= 2 (2|s|)^N for |s| >= 1.
(2) |J_0(y)| <= sqrt(2/(pi |y|)) for y != 0 (used for the density bound).
(3) The moment identity behind Proposition B: for a finite model whose relations of
    length <= L are supported in S, the moments of order <= L under Haar measure on the
    closure torus equal those of X + Y (X: the S-part under Haar measure on its own torus,
    Y: the rest with independent uniform phases).  Checked exactly (rational arithmetic on
    the relation lattice) for a small example; and an example where the order-(L+1)
    moments differ.
(4) Lemma A itself, on an example: Z = a lattice-supported variable whose first 2N moments
    match those of a continuous W (Gauss quadrature nodes), compared with the bound.
"""
import itertools
from fractions import Fraction

import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.special import j0

out = []
P = out.append


def needle_coeffs(N):
    """Chebyshev coefficients (in s = cos phi) of G(phi) = F_N(phi-pi/2) + F_N(phi+pi/2).

    F_N(phi) = sum_{|k|<=N} (1-|k|/(N+1)) e^{ik phi}.  G(phi) = sum_k w_k e^{ik phi}
    (e^{-ik pi/2} + e^{ik pi/2}) = sum_k w_k 2 cos(k pi/2) e^{ik phi}; as a cosine series
    G = w_0*2 + sum_{k>=1} 2 w_k 2cos(k pi/2) cos(k phi).
    """
    c = np.zeros(N + 1)
    for k in range(N + 1):
        w = 1 - k / (N + 1)
        ck = 2 * np.cos(k * np.pi / 2)
        c[k] = w * ck * (1 if k == 0 else 2)
    return c


P("(1) needle polynomial checks")
for N in [1, 2, 3, 5, 10, 20, 50, 100, 300]:
    c = needle_coeffs(N)
    g0 = C.chebval(0.0, c)
    s = np.linspace(-1, 1, 200001)
    q = C.chebval(s, c) / g0
    # integral by Clenshaw-Curtis-free simple quadrature (fine grid, smooth integrand)
    I = np.trapezoid(q ** 2, s)
    s0 = 2 / (N + 1)
    inner = np.linspace(-s0, s0, 2001)
    qmin_inner = (C.chebval(inner, c) / g0).min()
    so = np.concatenate([np.linspace(1, 50, 20000), -np.linspace(1, 50, 20000)])
    qo = np.abs(C.chebval(so, c) / g0)
    ratio = (qo / (2 * (2 * np.abs(so)) ** N)).max()
    P(f"  N={N:4d}: q(0)={C.chebval(0.0, c)/g0:.6f}  min/max on [-1,1]=({q.min():.3e},{q.max():.4f})"
      f"  int q^2={I:.5f} <= 8pi/(N+1)={8*np.pi/(N+1):.5f}: {I <= 8*np.pi/(N+1)}"
      f"  min q on |s|<=2/(N+1)={qmin_inner:.4f} (>=0.27: {qmin_inner >= 0.27})"
      f"  max |q|/(2(2|s|)^N) on 1<=|s|<=50 = {ratio:.3e}")

P("")
P("(2) |J_0(y)| <= sqrt(2/(pi y)), y in (0, 2000]")
y = np.concatenate([np.linspace(1e-6, 10, 200001), np.linspace(10, 2000, 400001)])
r = np.abs(j0(y)) / np.sqrt(2 / (np.pi * y))
P(f"  max ratio = {r.max():.6f} (must be <= 1)")

P("")
P("(3) moment identity (exact)")


def relations_small(freqs, Lmax):
    """All nonzero k in Z^n with ||k||_1 <= Lmax and sum k_i f_i = 0 (f_i rational)."""
    n = len(freqs)
    rels = []
    rng = range(-Lmax, Lmax + 1)
    for k in itertools.product(rng, repeat=n):
        if 0 < sum(abs(t) for t in k) <= Lmax and sum(ki * fi for ki, fi in zip(k, freqs)) == 0:
            rels.append(k)
    return rels


def moments_haar(freqs, coeffs, jmax, Smask=None):
    """Exact moments E[(Re sum c_i e_{f_i})^j], j<=jmax, under Haar measure of the closure
    torus of the frequencies (Smask=None), or of the mixed model (S-part Haar, rest
    independent uniform phases) when Smask is given.  Coefficients real for simplicity:
    Re(c e_f) = c cos(f y)."""
    n = len(freqs)
    res = []
    # E prod_{i} (c_{a_i}/2)(e_{s_i f_{a_i}}) over sequences; group by net vector k.
    for j in range(jmax + 1):
        tot = Fraction(0)
        for seq in itertools.product(range(n), repeat=j):
            for signs in itertools.product((1, -1), repeat=j):
                k = [0] * n
                for a, sg in zip(seq, signs):
                    k[a] += sg
                if Smask is None:
                    ok = sum(ki * fi for ki, fi in zip(k, freqs)) == 0
                else:
                    ok = all(k[i] == 0 for i in range(n) if not Smask[i]) and \
                        sum(k[i] * freqs[i] for i in range(n) if Smask[i]) == 0
                if ok:
                    term = Fraction(1)
                    for a in seq:
                        term *= Fraction(coeffs[a]) / 2
                    tot += term
        res.append(tot)
    return res


# Example: frequencies 1, 2, 3 (relations 1+2-3=0 (length 3), 2*1-2=0 (length 3))
# together with sqrt(2)-like "free" ones modelled by large primes' reciprocals is not
# rational; instead take frequencies whose only short relations live in S = {1,2,3}:
# f = (1, 2, 3, 1000003/7, 1000033/11).  Relations involving the last two need
# coefficients divisible by 7 resp. 11 (by the argument with distinct primes), so length
# >= 7 + ... ; we use L = 4.
freqs = [Fraction(1), Fraction(2), Fraction(3), Fraction(1000003, 7), Fraction(1000033, 11)]
coeffs = [Fraction(1, 2), Fraction(1, 3), Fraction(1, 5), Fraction(1, 7), Fraction(1, 11)]
L = 4
rels = relations_small(freqs, L)
Smask = [any(r[i] != 0 for r in rels) for i in range(len(freqs))]
P(f"  frequencies {[str(f) for f in freqs]}")
P(f"  relations of length <= {L}: {len(rels)}; support S = {[str(freqs[i]) for i in range(len(freqs)) if Smask[i]]}")
mH = moments_haar(freqs, coeffs, L, None)
mM = moments_haar(freqs, coeffs, L, Smask)
P(f"  moments j<=L (Haar on closure torus): {[str(m) for m in mH]}")
P(f"  moments j<=L (mixed model):           {[str(m) for m in mM]}")
P(f"  equal: {mH == mM}")
# Order L+1..: the relation 7*(1000003/7) - 1000003*1 = 0 is long; within j <= 6 the
# mixed and true moments still agree here.  To see a difference, use a model where a
# length-(L+1) relation involves a free frequency.
freqs2 = [Fraction(1), Fraction(2), Fraction(3), Fraction(5, 2), Fraction(1000033, 11)]
rels2 = relations_small(freqs2, 2)
Smask2 = [any(r[i] != 0 for r in rels2) for i in range(len(freqs2))]
mH2 = moments_haar(freqs2, coeffs, 4, None)
mM2 = moments_haar(freqs2, coeffs, 4, Smask2)
P(f"  second example {[str(f) for f in freqs2]}, L=2: S = {[str(freqs2[i]) for i in range(len(freqs2)) if Smask2[i]]}")
P(f"    moments j<=4 (Haar): {[str(m) for m in mH2]}")
P(f"    moments j<=4 (mixed, S from L=2): {[str(m) for m in mM2]}")
P(f"    equal up to j=2: {mH2[:3] == mM2[:3]};  equal up to j=4: {mH2 == mM2} (expected False:"
  f" the length-3 relation 2*(5/2) - 2 - 3 = 0 involves 5/2, outside S for L=2)")

P("")
P("(4) Lemma A on an example: Z = Gauss-quadrature discretisation of W = N(0,1)")
# Z takes the N+1 Gauss-Hermite nodes with weights: its moments match W's up to 2N+1.
# Lemma A bounds P(|Z - x| <= 2R/(N+1)) by 352 D R/(N+1) + 56 16^{-N}, D = 1/sqrt(2pi),
# R = 16(|x| + 2 sqrt(N) sigma), sigma = 1 (sub-Gaussian proxy of N(0,1)).
from numpy.polynomial.hermite_e import hermegauss
for N in [2, 5, 10, 20, 40, 80]:
    nodes, w = hermegauss(N + 1)
    w = w / w.sum()
    x = 0.0
    R = 16 * (abs(x) + 2 * np.sqrt(N))
    rad = 2 * R / (N + 1)
    lhs = w[np.abs(nodes - x) <= rad].sum()
    maxatom = w.max()
    bound = 352 * (1 / np.sqrt(2 * np.pi)) * R / (N + 1) + 56 * 16.0 ** (-N)
    P(f"  N={N:3d}: max atom of Z = {maxatom:.4f} (~ c/sqrt(N): {maxatom*np.sqrt(N):.3f}/sqrt(N));"
      f"  P(|Z|<=2R/(N+1)) = {lhs:.4f} <= bound {bound:.2f}: {lhs <= bound}")

open(__file__.replace('.py', '_results.txt'), 'w').write("\n".join(out) + "\n")
print("\n".join(out))


# ---------------------------------------------------------------------------
# (5) J_0(y) <= exp(-y^2/4) on [0, j_{0,1}) (Weierstrass product + Rayleigh sum 1/4)
from scipy.special import jn_zeros
out2 = []
Q = out2.append
j01 = jn_zeros(0, 1)[0]
yy = np.linspace(0, j01 * (1 - 1e-12), 400001)
gap = np.exp(-yy ** 2 / 4) - j0(yy)
Q("(5) J_0(y) <= exp(-y^2/4) on [0, j_{0,1}):")
Q(f"  j_(0,1) = {j01:.6f}; min of exp(-y^2/4) - J_0(y) = {gap.min():.3e} (must be >= 0);"
  f" Rayleigh sum of j_(0,k)^-2 over first 10^5 zeros = {np.sum(jn_zeros(0, 100000) ** -2.0):.6f} (-> 1/4)")

# (6) Lemma D(b) on an example: Y = sum_{k=1}^n c_k cos(theta_k), c_k in [a,b]; compare the
# true density maximum (by numerical inversion of the characteristic function) with the bound.
Q("(6) Lemma D(b): density of sum c_k cos(theta_k) vs the bound")
rng = np.random.default_rng(1)
for n, a, b in [(10, 0.05, 0.1), (100, 0.01, 0.03), (1000, 0.002, 0.01)]:
    c = rng.uniform(a, b, n)
    xi = np.linspace(-60 / a, 60 / a, 2000001)
    phi = np.ones_like(xi)
    for ck in c:
        phi *= j0(xi * ck)
    dens0 = np.trapezoid(phi, xi) / (2 * np.pi)   # density at 0 (the maximum, by symmetry/unimodality)
    theta = max(np.exp(-a * a / (b * b)), np.pi ** -0.5)
    # Lemma D(b) as stated: first term 2 sqrt(pi) / sqrt(Sigma_I), Sigma_I = sum |c|^2
    bound = (1 / (2 * np.pi)) * (2 * np.sqrt(np.pi) / np.sqrt(np.sum(c ** 2)) + 12 / (np.pi * a) * theta ** (n - 3))
    gauss = 1 / np.sqrt(2 * np.pi * np.sum(c ** 2) / 2)
    Q(f"  n={n:5d}, c in [{a},{b}]: density at 0 = {dens0:.4f}; Gaussian value {gauss:.4f}; bound {bound:.4f}; ok: {dens0 <= bound}")
open(__file__.replace('.py', '_results.txt'), 'a').write("\n".join(out2) + "\n")
print("\n".join(out2))
