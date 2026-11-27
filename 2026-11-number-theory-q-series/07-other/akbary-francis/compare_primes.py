#!/usr/bin/env python3
"""
compare_primes.py q a samples.txt [out.txt] -- sanity check of the explicit formula
(Akbary-Francis Thm 1.18) against exact prime data (sieve.c), double precision.

Actual:     G(x) = phi(q) sqrt(x) log(x) * log f(x; q, a)        (from the sieve)
Leading:    G1(u) = -( S_a(u) + 2 R_{q,a} ),   u = log x,
            S_a(u) = sum over listed zeros 2 w cos(gamma u - arg chi(a)), w = 1/(1/4+gamma^2)
            (Thm 1.18, m = 2, GRH_q: the sum (7) is -S_a)
Refined:    G2(u): the same main terms WITHOUT expanding F_s(x) = int_x^oo t^s g(t) dt
            in 1/log x (Lemma 6.2 / Prop 3.5 before Lemma 3.4 is applied):
               F_s(x) = x^(s-1)/u + s E_1((1-s) u)        (exact)
               log f ~ -(1/phi) sum_chi conj(chi(a)) sum_rho F_rho(x)/rho
                        - sum_{k=2..6} (-mu(k)) (R_k/phi) F_{1/k}(x),   R_k = #{b : b^k = a mod q}
            (the k-sum models psi - theta = -sum_{k>=2} mu(k) sum_{b^k = a} psi(x^(1/k); q, b)
             with psi(y; q, b) ~ y/phi; neglected: the oscillation of those psi's,
             relative size x^(-1/4), and the O(1/x) terms J-hat, u(x), ...).
Not rigorous; zeros truncated at the height of the zero file.
"""
import sys, math
import numpy as np
from scipy.special import exp1
from afcommon import load_zeros, chi_arg_table

q, a = int(sys.argv[1]), int(sys.argv[2])
phi = q - 1
zeros = load_zeros(q)
tab = chi_arg_table(q, a)
g = np.array([float(d) for _, d in zeros])
ph = np.array([2 * math.pi * float(tab[lab]) for lab, _ in zeros])
w = 1.0 / (0.25 + g * g)
Rk = {k: sum(1 for b in range(1, q) if pow(b, k, q) == a % q) for k in range(2, 7)}
MU = {2: -1, 3: -1, 4: 0, 5: -1, 6: 1}


def F(s, u):
    return np.exp((s - 1) * u) / u + s * exp1((1 - s) * u)


rows = []
for line in open(sys.argv[3]):
    f = line.split()
    u, x, k, pk, logf = float(f[0]), int(f[1]), int(f[2]), int(f[3]), float(f[6])
    if u < 9.0 or abs(u * 10 - round(u * 10)) > 1e-9:   # every 0.1 in u, from x ~ 8100
        continue
    ux = math.log(x)
    G = phi * math.sqrt(x) * ux * logf
    S = float(np.sum(2 * w * np.cos(g * ux - ph)))
    G1 = -(S + 2 * Rk[2])
    rho = 0.5 + 1j * g
    Z = float(np.sum(2 * (np.exp(-1j * ph) * F(rho, ux) / rho).real))
    bias = sum(-MU[k] * Rk[k] * F(1.0 / k, ux).real for k in range(2, 7))
    G2 = math.sqrt(x) * ux * (-Z - bias)
    rows.append((u, x, k, G, G1, G2))
out = open(sys.argv[4], "w") if len(sys.argv) > 4 else sys.stdout
out.write("# u x k=pi(x;q,a) G_actual G1_leading G2_refined (G = phi sqrt(x) log x log f)\n")
for r in rows:
    out.write("%.1f %d %d %.5f %.5f %.5f\n" % r)
A = np.array(rows)
for lo in (9, 12, 15, 18, 20):
    sel = A[:, 0] >= lo
    if sel.sum() < 3:
        continue
    d1 = A[sel, 3] - A[sel, 4]
    d2 = A[sel, 3] - A[sel, 5]
    print("u >= %2d (%3d samples): G in [%.3f, %.3f], sd(G) = %.3f; G-G1: mean %+.4f rms %.4f max|.| %.4f; "
          "G-G2: mean %+.4f rms %.4f max|.| %.4f; corr(G,G1) = %.4f corr(G,G2) = %.5f"
          % (lo, sel.sum(), A[sel, 3].min(), A[sel, 3].max(), A[sel, 3].std(), d1.mean(), math.sqrt((d1 ** 2).mean()),
             abs(d1).max(), d2.mean(), math.sqrt((d2 ** 2).mean()), abs(d2).max(),
             np.corrcoef(A[sel, 3], A[sel, 4])[0, 1], np.corrcoef(A[sel, 3], A[sel, 5])[0, 1]))
