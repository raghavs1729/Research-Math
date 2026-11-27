"""Referee 4: exact checks of the moment identities in INVLO.md Prop. B' and PRIMECORR.md Thm P1.

AI-produced, 2026-11-27, adversarial referee (referee 4).

Frequencies are integer vectors in a Z-basis of Q-independent reals (e.g. 1, sqrt2, sqrt3), so
Haar measure on the closure torus is the image of Haar measure on T^d and moments are constant
terms of Laurent polynomials (exact Fractions).

INVLO Prop. B':
 (B1) Gamma' = {1, 2, sqrt2, sqrt3, 2sqrt2+3sqrt3}: relations of length 3 (2*1 - 2) and 6
      (2 sqrt2 + 3 sqrt3 - (2sqrt2+3sqrt3)).  E Z^j = E W^j must hold for j <= L, and the first
      failure is at the length of the shortest relation not supported in S_L.
 (B2) sparse cyclic Gamma' = {1, 3, 5} (beta = 1): S_3 is EMPTY (no relation of length <= 3),
      contrary to "cyclic models are saturated as soon as L >= 3" (INVLO §4) taken literally.
 (B3) full cyclic {1,...,6}: S_3 = everything.
PRIMECORR Thm P1 (finite multiset models; Haar vs independent phases):
 (P1a) Delta_1 = 0, Delta_2 = sum_gamma 2 m(m-1) a^2.
 (P1b) odd k: mu_k = 0, E G^k = Delta_k >= 0; Delta_k > 0 iff an odd relation of length <= k exists
       (needs padding, not only P1(c)).
 (P1e) mu_4 = 12 s_2^2 - 6 s_4.
 (S3)  padding: a relation of length l contributes >= weight * (sum a^2)^{(k-l)/2} to Delta_k.
 (P4b) a relation of length 6 contributes 0 to Delta_2..Delta_5 and multinomial*weight to Delta_6.
"""
from fractions import Fraction as Fr
from itertools import product
from collections import defaultdict
import math

out = []
P = out.append


def mul(p, q):
    r = defaultdict(Fr)
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            r[e] += c1 * c2
    return {e: c for e, c in r.items() if c != 0}


def const_moments(terms, dim, jmax):
    """terms: list of (coef, exponent tuple). Returns [const term of f^j for j=0..jmax]."""
    f = defaultdict(Fr)
    for c, e in terms:
        f[tuple(e)] += c
    f = dict(f)
    zero = tuple([0] * dim)
    res = [Fr(1)]
    cur = {zero: Fr(1)}
    for j in range(1, jmax + 1):
        cur = mul(cur, f)
        res.append(cur.get(zero, Fr(0)))
    return res


def relations(vecs, Lmax):
    """All nonzero integer k with ||k||_1 <= Lmax and sum k_i vecs_i = 0."""
    n = len(vecs)
    d = len(vecs[0])
    found = []

    def rec(i, budget, k, acc):
        if i == n:
            if any(k) and all(a == 0 for a in acc):
                found.append(tuple(k))
            return
        for t in range(-budget, budget + 1):
            k.append(t)
            rec(i + 1, budget - abs(t), k, [a + t * b for a, b in zip(acc, vecs[i])])
            k.pop()
    rec(0, Lmax, [], [0] * d)
    return found


def B_prime_check(name, vecs, coeffs, Lvals, jmax):
    n = len(vecs)
    d = len(vecs[0])
    P(f"  {name}")
    # Z: Haar on closure torus
    termsZ = []
    for v, c in zip(vecs, coeffs):
        termsZ.append((c / 2, v))
        termsZ.append((c / 2, tuple(-x for x in v)))
    mZ = const_moments(termsZ, d, jmax)
    for L in Lvals:
        rels = relations(vecs, L)
        S = sorted({i for r in rels for i in range(n) if r[i] != 0})
        # W: S-part in basis coordinates, rest in fresh coordinates
        free = [i for i in range(n) if i not in S]
        D = d + len(free)
        termsW = []
        for i in range(n):
            if i in S:
                e = tuple(list(vecs[i]) + [0] * len(free))
            else:
                e = tuple([0] * d + [1 if f == i else 0 for f in free])
            termsW.append((coeffs[i] / 2, e))
            termsW.append((coeffs[i] / 2, tuple(-x for x in e)))
        mW = const_moments(termsW, D, jmax)
        agree = [j for j in range(jmax + 1) if mZ[j] == mW[j]]
        first_bad = next((j for j in range(jmax + 1) if mZ[j] != mW[j]), None)
        P(f"    L={L}: S_L = {S} (indices); moments agree for j <= L: {all(mZ[j] == mW[j] for j in range(L + 1))};"
          f" first disagreement at j = {first_bad}")
    return mZ


P("INVLO Prop B' (moment identity)")
# basis (1, sqrt2, sqrt3)
vecsB1 = [(1, 0, 0), (2, 0, 0), (0, 1, 0), (0, 0, 1), (0, 2, 3)]
coeffsB1 = [Fr(1, 2), Fr(1, 3), Fr(1, 5), Fr(1, 7), Fr(1, 11)]
B_prime_check("(B1) {1, 2, sqrt2, sqrt3, 2sqrt2+3sqrt3}", vecsB1, coeffsB1, [2, 3, 4, 5, 6], 8)
P("    expected: L=2: S empty, first failure at 3; L=3..5: S={0,1}, first failure at 6; L>=6: none")
vecsB2 = [(1,), (3,), (5,)]
coeffsB2 = [Fr(1), Fr(1), Fr(1)]
rels3 = relations(vecsB2, 3)
rels4 = relations(vecsB2, 4)
P(f"  (B2) sparse cyclic {{1,3,5}}: relations of length <= 3: {rels3}; of length <= 4: {len(rels4)}"
  f" (e.g. {rels4[:2]}): S_3 is empty, S_4 is everything")
B_prime_check("(B2) {1,3,5}", vecsB2, coeffsB2, [3, 4], 6)
vecsB3 = [(k,) for k in range(1, 7)]
rels = relations(vecsB3, 3)
S = sorted({i for r in rels for i in range(6) if r[i] != 0})
P(f"  (B3) full cyclic {{1..6}}: S_3 = {[i+1 for i in S]} (all)")
# (B4) a cyclic model inside the class C (M = 1/2, s = 0.1556): three frequencies n_i * beta,
# beta = 1e-4, close to the first three ordinates of L(s, chi_{-4}); shortest relation length.
nB4 = (60209, 102438, 129881)
beta = Fr(1, 10 ** 4)
best = None
Bmax = 400
for k2 in range(-Bmax, Bmax + 1):
    for k3 in range(-(Bmax - abs(k2)), Bmax - abs(k2) + 1):
        if k2 == 0 and k3 == 0:
            continue
        s_ = k2 * nB4[1] + k3 * nB4[2]
        if s_ % nB4[0] == 0:
            k1 = -s_ // nB4[0]
            L_ = abs(k1) + abs(k2) + abs(k3)
            if best is None or L_ < best[0]:
                best = (L_, (k1, k2, k3))
sig2 = sum(8 * 0.25 / (0.25 + float(n * beta) ** 2) for n in nB4)
P(f"  (B4) Gamma' = 1e-4 * {nB4} (= {[float(n*beta) for n in nB4]}), m = 1, M = 1/2: sigma^2 = {sig2:.4f} <= 0.1556;"
  f" shortest relation (search ||k||_1 <= {Bmax}): length {best[0]}, k = {best[1]}")
P(f"       so S_L = empty for all L < {best[0]}: a cyclic model in the class C that is NOT saturated at those lengths;"
  f" smaller beta pushes this further (not uniform over C).")
P("  (B5) a Sidon set need not be relation-free: {1,3,4} is Sidon (pairwise sums 2,4,5,6,7,8 distinct)"
  f" but has the length-3 relation 1 + 3 - 4 = 0: relations(len<=3) = {relations([(1,), (3,), (4,)], 3)}")
P("  Gauss-integer example (complex coefficients) is covered by linearity; real case suffices for the identity.")

P("")
P("PRIMECORR Thm P1 (duality on finite multiset models)")


def P1_model(name, zeros, weights, kmax, enum_k=5):
    """zeros: list of ordinate vectors (with repetition for multiplicity); weights a_rho (Fractions)."""
    d = len(zeros[0])
    n = len(zeros)
    terms = []
    for v, a in zip(zeros, weights):
        terms.append((a, v))
        terms.append((a, tuple(-x for x in v)))
    EH = const_moments(terms, d, kmax)
    termsI = []
    for i, a in enumerate(weights):
        e = tuple(1 if j == i else 0 for j in range(n))
        termsI.append((a, e))
        termsI.append((a, tuple(-x for x in e)))
    mu = const_moments(termsI, n, kmax)
    Delta = [EH[k] - mu[k] for k in range(kmax + 1)]
    P(f"  {name}")
    P(f"    E_Haar G^k  : {[str(x) for x in EH]}")
    P(f"    mu_k (ind.) : {[str(x) for x in mu]}")
    P(f"    Delta_k     : {[str(x) for x in Delta]}  (all >= 0: {all(x >= 0 for x in Delta)})")
    # brute-force Delta_k for small k by enumerating sequences
    for k in range(1, min(enum_k, kmax) + 1):
        tot = Fr(0)
        for seq in product(range(n), repeat=k):
            for sg in product((1, -1), repeat=k):
                kap = [0] * n
                for r, s in zip(seq, sg):
                    kap[r] += s
                if not any(kap):
                    continue
                acc = [0] * d
                for i in range(n):
                    if kap[i]:
                        acc = [a + kap[i] * b for a, b in zip(acc, zeros[i])]
                if all(a == 0 for a in acc):
                    w = Fr(1)
                    for r in seq:
                        w *= weights[r]
                    tot += w
        assert tot == Delta[k], (k, tot, Delta[k])
    P(f"    brute-force sum over sequences with net vector a multiset relation = Delta_k for k <= {min(enum_k, kmax)}: True")
    return EH, mu, Delta


# Model 1: ordinates 1 (double), 2, sqrt2 (simple); basis (1, sqrt2).
z1 = [(1, 0), (1, 0), (2, 0), (0, 1)]
a1 = [Fr(1, 2), Fr(1, 2), Fr(1, 3), Fr(1, 5)]
EH, mu, De = P1_model("model 1: ordinates {1 (m=2), 2, sqrt2}", z1, a1, 6, enum_k=5)
d2 = 2 * 2 * 1 * Fr(1, 2) ** 2
P(f"    (P1a) Delta_1 = {De[1]} (0); Delta_2 = {De[2]} vs sum 2m(m-1)a^2 = {d2}: {De[2] == d2}")
P(f"    (P1b) odd k: mu_3 = {mu[3]}, mu_5 = {mu[5]}; Delta_3 = {De[3]} > 0 (odd relations 1+1-2, 2*1-2 exist): M_3 = -Delta_3 < 0")
s2 = sum(a ** 2 for a in a1)
s4 = sum(a ** 4 for a in a1)
P(f"    (P1e) mu_4 = {mu[4]} vs 12 s2^2 - 6 s4 = {12*s2**2 - 6*s4}: {mu[4] == 12*s2**2 - 6*s4}")

# Model 2: one odd relation of length 3 only below length 5: ordinates 1, 2, sqrt2 (simple).
z2 = [(1, 0), (2, 0), (0, 1)]
a2 = [Fr(1, 2), Fr(1, 3), Fr(1, 5)]
EH2, mu2, De2 = P1_model("model 2: ordinates {1, 2, sqrt2} (only relation family 2*1 - 2 = 0)", z2, a2, 7, enum_k=5)
wk = a2[0] ** 2 * a2[1]
mult3 = 3  # orderings of the multiset {1,1,2} with signs fixed: 3!/2! = 3, times 2 for overall sign flip
P(f"    Delta_3 = {De2[3]}; single relation (2,-1) and its negative: 2 * 3!/(2!1!) * a1^2 a2 = {2*mult3*wk}")
pad = wk * sum(a ** 2 for a in a2)
P(f"    (S3) padding: Delta_5 = {De2[5]} >= weight * (sum a^2)^1 = {pad}: {De2[5] >= pad}"
  f" (ratio {float(De2[5]/pad):.2f}); Delta_7 = {De2[7]} >= weight*(sum a^2)^2 = {wk*sum(a**2 for a in a2)**2}")
P(f"    P1(b) 'iff' for k = 5 needs padding: Delta_5 > 0 although the only odd relations have length 3 (and 6+...)")

# Model 3: only a relation of length 6 (2 sqrt2 + 3 sqrt3 - (2sqrt2 + 3sqrt3) = 0), basis (sqrt2, sqrt3)
z3 = [(1, 0), (0, 1), (2, 3)]
a3 = [Fr(1, 2), Fr(1, 2), Fr(1, 2)]
EH3, mu3, De3 = P1_model("model 3: ordinates {sqrt2, sqrt3, 2sqrt2+3sqrt3} (shortest relation has length 6)", z3, a3, 6, enum_k=4)
multi = 2 * math.factorial(6) // (math.factorial(2) * math.factorial(3) * math.factorial(1))
w6 = a3[0] ** 2 * a3[1] ** 3 * a3[2]
P(f"    (P4b) Delta_2..Delta_5 = {[str(x) for x in De3[2:6]]} (all 0); Delta_6 = {De3[6]} = 2*6!/(2!3!1!) * weight = {multi*w6}")
P("    So a relation of length l contributes nothing to Delta_k for k < l, and contributes a multinomial")
P("    multiple of its weight at k = l: per-relation contributions are not 'of order 1' uniformly in l.")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
