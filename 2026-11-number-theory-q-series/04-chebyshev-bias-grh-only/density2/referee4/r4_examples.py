"""Referee 4: checks of the abstract examples of INVLO.md §3.

AI-produced, 2026-11-27, adversarial referee (referee 4).

Item 2 (rank one, GLI(psi) holds, TS/FD/TS_Z fail): psi(t) = 2t (so psi(t) log t / t -> infinity),
alpha = 1, e_n increasing with e_n - e_{n-1} >= log2 psi(n+1), u_n odd primes with
x_n = u_n / 2^{e_n} in (n, n+1).  We
 (i) build the first six x_n;
 (ii) for n* = 2..5 find, by exhaustive search (eliminating k_1), every relation
      sum_{i<=n*} k_i x_i = 0 with k_{n*} != 0 and ||k||_1 <= B, and compare the minimal length
      with 2^{e_{n*}-e_{n*-1}} and psi(x_{n*});
 (iii) check gcd(u_m 2^{e_m'-e_m}, u_m') = 1 (so <x_m, x_m'> = 2^{-e_m'} Z and the tail generates
      Z[1/2] alpha: TS_Z fails), and that X is inside Q alpha (TS and FD fail).
Item 3: X = {beta} u {n alpha}: the length-3 relations (n+1)alpha - n alpha - alpha = 0.
"""
from fractions import Fraction as Fr
from math import gcd, log2, ceil
from sympy import nextprime, isprime

out = []
P = out.append


def psi(t):
    return 2.0 * t


xs, es, us = [], [], []
e_prev = None
for n in range(1, 7):
    e = 0 if e_prev is None else e_prev + max(0, ceil(log2(psi(n + 1))))
    while True:
        lo = n * 2 ** e
        hi = (n + 1) * 2 ** e
        p = nextprime(lo)
        while p in us or p == 2:
            p = nextprime(p)
        if p < hi:
            break
        e += 1
    es.append(e)
    us.append(p)
    xs.append(Fr(p, 2 ** e))
    e_prev = e
P("Item 2 construction (alpha = 1, psi(t) = 2t):")
for n, (e, u, x) in enumerate(zip(es, us, xs), 1):
    P(f"  n={n}: e_n={e}, u_n={u} (prime: {isprime(u)}), x_n={float(x):.6f} in ({n},{n+1}): {n < x < n+1};"
      f" e_n - e_(n-1) = {e - es[n-2] if n > 1 else '-'} >= log2 psi(n+1) = {log2(psi(n+1)):.3f}")

P("")
P("  (ii) shortest relations with largest index n* (exhaustive, ||k||_1 <= B):")
for nstar, B in [(2, 80), (3, 60), (4, 40), (5, 26)]:
    X = xs[:nstar]
    best = None
    count = 0
    divis_ok = True

    def rec(i, budget, k):
        global best, count, divis_ok
        if i == nstar:
            if k[-1] == 0:
                return
            s = sum(ki * xi for ki, xi in zip(k, X[1:]))
            k1 = -s / X[0]
            if k1.denominator == 1 and abs(k1.numerator) <= budget:
                kk = [int(k1)] + list(k)
                L = sum(abs(t) for t in kk)
                count += 1
                Eprime = max(es[j] for j in range(nstar - 1) if kk[j] != 0)
                if kk[-1] % (2 ** (es[nstar - 1] - Eprime)) != 0:
                    divis_ok = False
                if best is None or L < best[0]:
                    best = (L, kk)
            return
        for t in range(-budget, budget + 1):
            k.append(t)
            rec(i + 1, budget - abs(t), k)
            k.pop()
    rec(1, B, [])
    need = 2 ** (es[nstar - 1] - es[nstar - 2])
    if best is None:
        P(f"    n*={nstar}: no relation with ||k||_1 <= {B}; claim: length >= 2^(e_n*-e_(n*-1)) = {need}"
          f" >= psi(x_n*) = {psi(float(X[-1])):.2f}  (consistent)")
    else:
        P(f"    n*={nstar}: {count} relations found, minimal length {best[0]} (k = {best[1]});"
          f" claim: >= {need} >= psi(x_n*) = {psi(float(X[-1])):.2f}: {best[0] >= need >= psi(float(X[-1]))};"
          f" 2^(e_n*-E') | k_n* in all found: {divis_ok}")

P("")
P("  (iii) TS_Z: gcd(u_m 2^(e_m'-e_m), u_m') for m < m':")
ok = True
for m in range(6):
    for mp_ in range(m + 1, 6):
        g = gcd(us[m] * 2 ** (es[mp_] - es[m]), us[mp_])
        ok &= (g == 1)
P(f"    all equal to 1: {ok}; so <x_m, x_m'> = 2^(-e_m') Z, and since e_n -> infinity the tail X_(>T)"
  f" generates Z[1/2] = <X> for every T (index 1): TS_Z fails.")
P("    X is contained in Q (= Q alpha): the Q-span of any infinite tail is Q, so TS fails; for FD,"
  " X_T minus any N_0 elements still spans Q once #X_T > N_0, so FD fails.")

P("")
P("Item 3: X = {beta} u {n alpha}, beta not in Q alpha.")
P("  relation (n+1)alpha - n alpha - alpha = 0: length 3, largest element (n+1)alpha -> infinity;"
  " so GLI(psi) fails as soon as psi > 3 eventually; TS holds with the form dual to beta.")
P("  Note: for psi <= 3 (e.g. psi = 0) GLI(psi) is vacuous, so 'incomparable' (item 4) needs psi > 3"
  " eventually; Theorem C's hypothesis psi(t) log t / t -> infinity guarantees this.")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
