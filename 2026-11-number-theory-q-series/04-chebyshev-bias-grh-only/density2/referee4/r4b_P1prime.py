"""Referee 4b: exact checks of PRIMECORR.md Lemma P1' (band duality) on finite multiset models.

AI-produced, 2026-11-27, adversarial referee (referee 4, second read).

Model: zeros rho with ordinate vectors in a Z-basis of Q-independent reals; R-zeros carry
weights beta_R(gamma) (1 at heights <= T*, in [0,1] above), band zeros carry weight 1.
  G_R = sum_R w_rho (z^v + z^-v),  G_I = sum_I (z^v + z^-v).
E_{a,b} := E_Haar[G_R^a G_I^b] - E_Haar[G_R^a] * mu_b(n_I),  mu_b(n) = E[(sum_{i<=n} 2cos)^b].
Checks:
 (Q1) E_{a,b} >= 0 and equals the brute-force sum over sequences whose net multiset vector is a
      relation with kappa_I != 0 (Lemma P1').
 (Q2) allowed relations among R-zeros contribute nothing, while the plain moment excess of
      G_R + G_I is positive (the point of band decoupling).
 (Q3) a band zero in a short relation with R-zeros of height <= T* gives E >= 1 (violates (D)).
 (Q4) a DOUBLE zero in the band gives E_{0,2} = 4 > 1: (D) forces simple band zeros, so GLI
      (which allows multiple zeros) does not imply H_dec.
 (Q5) an EMPTY band gives E_{a,b} = 0 for all a, b: (D) holds vacuously.
 (Q6) a relation through an R-zero above T* (weight 1/3) can stay below the threshold 1.
"""
from fractions import Fraction as Fr
from itertools import product
from collections import defaultdict

out = []
P = out.append


def mul(p, q):
    r = defaultdict(Fr)
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            r[tuple(a + b for a, b in zip(e1, e2))] += c1 * c2
    return {e: c for e, c in r.items() if c != 0}


def power_const(f, dim, k, extra=None):
    """constant terms of f^j * extra (dict) for j = 0..k; returns list of polys too."""
    zero = tuple([0] * dim)
    cur = {zero: Fr(1)}
    polys = [cur]
    for j in range(1, k + 1):
        cur = mul(cur, f)
        polys.append(cur)
    return polys


def laurent(zs, ws, dim):
    f = defaultdict(Fr)
    for v, w in zip(zs, ws):
        f[tuple(v)] += w
        f[tuple(-x for x in v)] += w
    return dict(f)


def mu_b(n, b):
    if n == 0:
        return Fr(1) if b == 0 else Fr(0)
    terms = []
    for i in range(n):
        e = tuple(1 if j == i else 0 for j in range(n))
        terms.append((e, Fr(1)))
        terms.append((tuple(-x for x in e), Fr(1)))
    f = defaultdict(Fr)
    for e, c in terms:
        f[e] += c
    polys = power_const(dict(f), n, b)
    return polys[b].get(tuple([0] * n), Fr(0))


def analyse(name, R, wR, I, Lmax=5, brute=True):
    dim = len((R + I)[0])
    fR = laurent(R, wR, dim) if R else {}
    fI = laurent(I, [Fr(1)] * len(I), dim) if I else {}
    zero = tuple([0] * dim)
    PR = power_const(fR, dim, Lmax) if R else [{zero: Fr(1)}] + [{}] * Lmax
    PI = power_const(fI, dim, Lmax) if I else [{zero: Fr(1)}] + [{}] * Lmax
    nI = len(I)
    P(f"  {name}: n_I = {nI}")
    res = {}
    for a in range(0, Lmax + 1):
        for b in range(1, Lmax + 1 - a):
            prod_ = mul(PR[a], PI[b]) if PR[a] and PI[b] else {}
            joint = prod_.get(zero, Fr(0))
            dec = PR[a].get(zero, Fr(0)) * mu_b(nI, b)
            E = joint - dec
            res[(a, b)] = E
            if brute and a + b <= 4:
                # brute force over sequences: a R-terms then b I-terms
                tot = Fr(0)
                for seqR in product(range(len(R)), repeat=a):
                    for sR in product((1, -1), repeat=a):
                        for seqI in product(range(nI), repeat=b):
                            for sI in product((1, -1), repeat=b):
                                kR = [0] * len(R)
                                kI = [0] * nI
                                for r_, s_ in zip(seqR, sR):
                                    kR[r_] += s_
                                for r_, s_ in zip(seqI, sI):
                                    kI[r_] += s_
                                if not any(kI):
                                    continue
                                acc = [0] * dim
                                for i, k in enumerate(kR):
                                    if k:
                                        acc = [x + k * y for x, y in zip(acc, R[i])]
                                for i, k in enumerate(kI):
                                    if k:
                                        acc = [x + k * y for x, y in zip(acc, I[i])]
                                if all(x == 0 for x in acc):
                                    w = Fr(1)
                                    for r_ in seqR:
                                        w *= wR[r_]
                                    tot += w
                assert tot == E, (name, a, b, tot, E)
    nonneg = all(v >= 0 for v in res.values())
    P(f"    E_(a,b) for a+b <= {Lmax}: " + ", ".join(f"({a},{b}):{str(v)}" for (a, b), v in sorted(res.items()) if v != 0)
      + ("" if any(v != 0 for v in res.values()) else "all zero"))
    bf = "brute-force sequence count = E for a+b <= 4: True" if brute else "brute force not run"
    P(f"    all E >= 0: {nonneg}; {bf}; (D) (all E < 1): {all(v < 1 for v in res.values())}")
    return res


P("Lemma P1' on finite multiset models (basis 1, sqrt2, sqrt3, sqrt5, sqrt7)")
e = lambda *c: tuple(c)
one, r2, r3, r5, r7 = e(1, 0, 0, 0, 0), e(0, 1, 0, 0, 0), e(0, 0, 1, 0, 0), e(0, 0, 0, 1, 0), e(0, 0, 0, 0, 1)
two = e(2, 0, 0, 0, 0)
# (Q1)/(Q2): R has the relation 2*1 - 2 = 0; band free
res = analyse("(Q2) R = {1, 2} (relation 2*1-2), I = {sqrt3, sqrt5, sqrt7}", [one, two], [Fr(1), Fr(1)], [r3, r5, r7])
# plain moment excess of G_R + G_I (one window, independent-phase prediction for all zeros)
dim = 5
fall = laurent([one, two, r3, r5, r7], [Fr(1)] * 5, dim)
polys = power_const(fall, dim, 5)
from math import comb
mu_all = [mu_b(5, k) for k in range(6)]
P(f"    plain moments of G_R+G_I minus independent prediction, k = 3..5: "
  f"{[str(polys[k].get((0,)*5, Fr(0)) - mu_all[k]) for k in range(3, 6)]} (positive: the R-relation is seen)")
analyse("(Q3) R = {1, sqrt2}, I = {1+sqrt2, sqrt3, sqrt5} (relation (1+sqrt2) - 1 - sqrt2)",
        [one, r2], [Fr(1), Fr(1)], [e(1, 1, 0, 0, 0), r3, r5])
analyse("(Q4) R = {1}, I = {sqrt2 (double), sqrt3, sqrt5}", [one], [Fr(1)], [r2, r2, r3, r5], Lmax=4)
P("    -> a double zero in the band violates (D) at (a,b) = (0,2): H_dec forces simple band zeros.")
analyse("(Q5) R = {1, 2}, I = empty", [one, two], [Fr(1), Fr(1)], [], Lmax=5, brute=False)
P("    -> with an empty band, (D) holds for every L: H_dec is vacuous (n_I = 0, mu_b(0) = 0 for b >= 1).")
analyse("(Q6) R = {1, sqrt2 [height > T*, weight 1/3]}, I = {1+sqrt2, sqrt3, sqrt5}",
        [one, r2], [Fr(1), Fr(1, 3)], [e(1, 1, 0, 0, 0), r3, r5], Lmax=4)
P("    -> the relation through the R-zero in (T*, T*+1] has weight 1/3 per sequence; (D) at (2,1) gives 4/3.")
P("       Such relations may or may not be excluded by (D) (here they are); Theorem P7 uses only the")
P("       exclusion of relations at heights <= T*, which Lemma P1' guarantees (weight exactly 1).")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
