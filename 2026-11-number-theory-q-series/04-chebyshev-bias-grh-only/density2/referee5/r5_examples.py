"""Referee 5: comparison examples and remarks of INVLO.md / PRIMECORR.md.

AI-produced adversarial referee check, 2026-11-27.

(C1) PRIMECORR Sec. 2 ("What is new"): "a double zero anywhere in (T, pi T] violates
     H_corr(2,T)".  With a double zero gamma_0, Delta_2 >= 4 what(gamma_0/T)^2, and H_corr(2,T)
     asks Delta_2 < c_0^2, c_0 = min_{|xi|<=1} what.  So the claim needs 2 what(xi) >= c_0 on
     (1, pi).  Test two admissible windows v in C_c^infty(-1/2,1/2), v >= 0 even.
(C2) INVLO Sec. 3 item 2: x_n = alpha u_n / 2^{e_n} in (n, n+1), e_n - e_{n-1} >= log2 psi(n+1),
     u_n distinct odd primes; psi(t) = 2^t.  Exhaustive search of relations among x_1..x_4 with
     ||k||_1 <= 80; check ||k||_1 >= |k_{n*}| >= 2^{e_{n*}-e_{n*-1}} >= psi(x_{n*}); check that
     <x_m, x_m'> = alpha 2^{-e_m'} Z (index 1, TS_Z fails).
(C3) INVLO Sec. 4: {beta, 3 beta, 5 beta} has no relation of length <= 3 (shortest: 4).
(C4) PRIMECORR Remark P4(b): alpha = 0.0891, delta = 0.1, T = 1e3, C_2 = 1 -> L >= ?
(C5) Lemma D(b) on true chi_{-4} bands I = Gamma cap (T, 2T], M = 1/2 (Theorem C / P7 input):
     number of ordinates, the bound, and (T/log T)^{1/2}.
"""
import math
import os

import numpy as np
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
ZF = os.path.join(HERE, "..", "..", "..", "prototype", "explore", "chi4_zeros.txt")
out = []
P = out.append

# ---------------------------------------------------------------- (C1)
P("(C1) window dependence of 'a double zero anywhere in (T, pi T] violates H_corr(2,T)'")
U = np.linspace(-0.5, 0.5, 200001)[1:-1]
du = U[1] - U[0]


def bump(u, c, wdt):
    x = (u - c) / wdt
    o = np.zeros_like(u)
    m = np.abs(x) < 1
    o[m] = np.exp(-1 / (1 - x[m] ** 2))
    return o


windows = {
    "standard bump exp(-1/(1-4u^2))": bump(U, 0.0, 0.5),
    "two bumps at +-0.45, half-width 0.04": bump(U, 0.45, 0.04) + bump(U, -0.45, 0.04),
}
xi = np.linspace(0, np.pi, 3142)
for name, v in windows.items():
    v = v / (v.sum() * du)
    assert v[np.abs(U) > 0.4999].max() < 1e-12
    vh = np.array([np.sum(v * np.cos(x * U)) * du for x in xi])
    wh = vh ** 2
    c0 = wh[xi <= 1].min()
    bad = xi[(xi > 1) & (2 * wh < c0)]
    P(f"  {name}: c_0 = {c0:.4f}; what(pi) = {wh[-1]:.4f}; min_(1,pi] 2*what = {2 * wh[xi > 1].min():.4f};"
      f" xi in (1,pi] where a double zero at xi*T does NOT violate H_corr(2,T): "
      f"{'none' if len(bad) == 0 else f'[{bad.min():.3f}, {bad.max():.3f}]'}")

# ---------------------------------------------------------------- (C2)
P("")
P("(C2) INVLO Sec. 3 item 2, psi(t) = 2^t, alpha = 1")
psi = lambda t: 2.0 ** t
e, u = [], []
prev = None
for n in range(1, 5):
    en = 1 if prev is None else prev + max(1, math.ceil(math.log2(psi(n + 1))))
    while True:
        lo, hi = n * 2 ** en, (n + 1) * 2 ** en
        p = sympy.nextprime(lo)
        if p < hi and p % 2 == 1 and p not in u:
            break
        en += 1
    e.append(en)
    u.append(int(p))
    prev = en
x = [ui / 2 ** ei for ui, ei in zip(u, e)]
P(f"  e = {e}; u = {u}; x = {[round(t, 6) for t in x]}")
E = e[-1]
A = [ui * 2 ** (E - ei) for ui, ei in zip(u, e)]  # sum k_n A_n = 0  <=>  sum k_n x_n = 0
Lmax = 80
found = []
r = np.arange(-Lmax, Lmax + 1)
K1, K2, K3 = np.meshgrid(r, r, r, indexing="ij")
K1, K2, K3 = K1.ravel(), K2.ravel(), K3.ravel()
l1 = np.abs(K1) + np.abs(K2) + np.abs(K3)
keep = l1 <= Lmax
K1, K2, K3, l1 = K1[keep], K2[keep], K3[keep], l1[keep]
S = K1 * A[0] + K2 * A[1] + K3 * A[2]
ok = (S % A[3] == 0)
K4 = -S[ok] // A[3]
for k1, k2, k3, k4 in zip(K1[ok], K2[ok], K3[ok], K4):
    k = (int(k1), int(k2), int(k3), int(k4))
    if any(k) and sum(map(abs, k)) <= Lmax:
        found.append(k)
viol = []
for k in found:
    ns = max(i for i in range(4) if k[i] != 0)
    L1 = sum(map(abs, k))
    lowbd = 2 ** (e[ns] - (e[ns - 1] if ns > 0 else 0))
    if not (L1 >= psi(x[ns]) and abs(k[ns]) >= lowbd):
        viol.append(k)
short = sorted(found, key=lambda k: sum(map(abs, k)))[:3]
P(f"  relations with ||k||_1 <= {Lmax}: {len(found)}; shortest: {short}; GLI(psi)/divisibility violations: {len(viol)}")
for ns in range(1, 4):
    top = [k for k in found if max(i for i in range(4) if k[i] != 0) == ns]
    if top:
        kk = min(top, key=lambda k: sum(map(abs, k)))
        P(f"  shortest relation with top index n*={ns + 1}: k = {kk}, length {sum(map(abs, kk))};"
          f" |k_n*| = {abs(kk[ns])} >= 2^(e_n* - e_(n*-1)) = {2 ** (e[ns] - e[ns - 1])} >= psi(x_n*) = {psi(x[ns]):.2f}")
    else:
        P(f"  no relation with top index n*={ns + 1} and length <= {Lmax}"
          f" (consistent: 2^(e_n* - e_(n*-1)) = {2 ** (e[ns] - e[ns - 1])})")
idx_ok = all(math.gcd(u[m] * 2 ** (e[mp] - e[m]), u[mp]) == 1 for m in range(4) for mp in range(m + 1, 4))
P(f"  gcd(u_m 2^(e_m'-e_m), u_m') = 1 for all m < m': {idx_ok}  (so <x_m, x_m'> = 2^(-e_m') Z; index 1)")

# ---------------------------------------------------------------- (C3)
P("")
P("(C3) {beta, 3beta, 5beta}: shortest relation")
best = None
for k1 in range(-8, 9):
    for k3 in range(-8, 9):
        for k5 in range(-8, 9):
            if (k1, k3, k5) != (0, 0, 0) and k1 + 3 * k3 + 5 * k5 == 0:
                L1 = abs(k1) + abs(k3) + abs(k5)
                if best is None or L1 < best[0]:
                    best = (L1, (k1, k3, k5))
P(f"  shortest relation length {best[0]}, e.g. k = {best[1]} (claim: none of length <= 3)")

# ---------------------------------------------------------------- (C4)
P("")
P("(C4) Remark P4(b): L log L >= C_2 T / (alpha^2 delta)")
al, de, T, C2 = 0.0891, 0.1, 1e3, 1.0
target = C2 * T / (al ** 2 * de)
lo, hi = 2.0, 1e12
for _ in range(200):
    mid = math.sqrt(lo * hi)
    if mid * math.log(mid) < target:
        lo = mid
    else:
        hi = mid
P(f"  C_2 T/(alpha^2 delta) = {target:.4e}; smallest L with L log L >= that: {hi:.4e} (file: 1.1e5)")

# ---------------------------------------------------------------- (C5)
P("")
P("(C5) Lemma D(b) for the true bands Gamma cap (T, 2T] of L(s, chi_{-4}), M = 1/2, m_gamma = 1")
Z = np.array([float(l) for l in open(ZF) if l.strip()])
M = 0.5
for T in [20, 50, 100, 200, 300]:
    I = Z[(Z > T) & (Z <= 2 * T)]
    c = 4 * M / np.sqrt(M ** 2 + I ** 2)
    a, bb = c.min(), c.max()
    th = max(math.exp(-a * a / bb ** 2), math.pi ** -0.5)
    Sig = np.sum(c ** 2)
    t1 = 2 * math.sqrt(math.pi) / math.sqrt(Sig) / (2 * math.pi)
    t2 = 12 * th ** (len(I) - 3) / (math.pi * a) / (2 * math.pi)
    P(f"  T={T:4d}: #I = {len(I):4d} (#I/T = {len(I) / T:.3f}); a/b = {a / bb:.3f}; theta = {th:.4f};"
      f" D-bound = {t1:.3f} + {t2:.2e}; Gaussian peak 1/sqrt(pi Sigma) = {1 / math.sqrt(math.pi * Sig):.3f};"
      f" (T/log T)^(1/2)/M = {math.sqrt(T / math.log(T)) / M:.2f}")

open(os.path.join(HERE, "r5_examples_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
