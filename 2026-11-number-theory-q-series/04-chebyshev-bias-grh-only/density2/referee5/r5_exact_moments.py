"""Referee 5: exact (rational, group-ring) checks of the moment identities.

AI-produced adversarial referee check, 2026-11-27.

Frequencies are integer vectors in Z^d (coordinates in a Q-basis of formal independent
generators), so relations are decided exactly.  A trigonometric polynomial
sum_rho a_rho (e_{v_rho} + e_{-v_rho}) is a dict {exponent vector: coefficient}; the Haar
expectation on the closure torus is the coefficient of the zero vector.

(E1) INVLO Prop. B': moments of order <= L of Z (Haar on the closure torus) equal those of
     W = 1 + X + Y (X: S_L-part, Haar on its own torus; Y: rest, independent phases); they
     differ at order 6 because of a length-6 relation through the free part.
(E2) PRIMECORR Thm P1(b), padding: ordinates {1, 2, sqrt2}, weights 1: Delta_3, Delta_5.
(E3) PRIMECORR Lemma P1' and the closing Remark of Sec. 3: R-zeros {1,2} (relation of length
     3 among R-zeros), band zeros {x,y,z} independent: E_{a,b} = 0 for all a+b <= 6, while the
     plain-moment excess of G_R + G_I is 6, 0, 460 at orders 3, 4, 5.  Also a double band zero
     gives E_{0,2} = 4.  E_{a,b} is computed twice: from the definition (difference of
     expectations) and by enumerating sequences with kappa_I != 0 and total 0.
(E4) PRIMECORR Sec. 3, comparison, the spectrum "{n alpha} u U B_j u {z_x, x+z_x}" AS WRITTEN:
     the band (s1,s2) contains multiples of alpha, and (n+1)alpha - n alpha - alpha = 0 is a
     relation of length 3 at height <= T* with kappa_I != 0: E_{1,2} >= 1, so (D) FAILS.
(E5) The same spectrum with the multiples of alpha removed from [s1-1, s2+1]: E_{a,b} = 0
     for all a + b <= 5, while GLI(psi) still fails (24 alpha - 23 alpha - alpha = 0).
"""
import itertools
import os
from collections import defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
P = out.append


def poly_from(terms, dim):
    """terms: list of (vec, coeff); returns dict for sum coeff (e_vec + e_{-vec})."""
    G = defaultdict(Fraction)
    for v, c in terms:
        v = tuple(v)
        mv = tuple(-t for t in v)
        G[v] += c
        G[mv] += c
    return dict(G)


def mult(A, B):
    C = defaultdict(Fraction)
    for va, ca in A.items():
        for vb, cb in B.items():
            C[tuple(x + y for x, y in zip(va, vb))] += ca * cb
    return {k: v for k, v in C.items() if v != 0}


def power(G, k, dim):
    cur = {tuple([0] * dim): Fraction(1)}
    for _ in range(k):
        cur = mult(cur, G)
    return cur


def const(A, dim):
    return A.get(tuple([0] * dim), Fraction(0))


def l1_ball(n, L):
    """all k in Z^n with 0 < ||k||_1 <= L."""
    res = []

    def rec(i, rem, cur):
        if i == n:
            if any(cur):
                res.append(tuple(cur))
            return
        for t in range(-rem, rem + 1):
            cur.append(t)
            rec(i + 1, rem - abs(t), cur)
            cur.pop()
    rec(0, L, [])
    return res


def S_L(vecs, L):
    """indices occurring in some nonzero relation of length <= L."""
    d = len(vecs[0])
    S = set()
    for k in l1_ball(len(vecs), L):
        if all(sum(k[i] * vecs[i][j] for i in range(len(vecs))) == 0 for j in range(d)):
            S |= {i for i in range(len(vecs)) if k[i] != 0}
    return S


# ---------------------------------------------------------------------------
P("(E1) INVLO Proposition B' moment identity (exact)")
# basis (1, sqrt2, sqrt3); frequencies 1, 2, 3, sqrt2, 1+sqrt2, sqrt3, 5 sqrt3
names = ["1", "2", "3", "r2", "1+r2", "r3", "5r3"]
vecs = [(1, 0, 0), (2, 0, 0), (3, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (0, 0, 5)]
coef = [Fraction(1, 2), Fraction(1, 3), Fraction(1, 5), Fraction(2, 7), Fraction(1, 4), Fraction(3, 5), Fraction(1, 6)]
d = 3
for L in [3, 5]:
    S = S_L(vecs, L)
    P(f"  L={L}: S_L = {[names[i] for i in sorted(S)]}; free = {[names[i] for i in range(len(vecs)) if i not in S]}")
    Zm1 = poly_from([(vecs[i], coef[i] / 2) for i in range(len(vecs))], d)
    free = [i for i in range(len(vecs)) if i not in S]
    dimW = d + len(free)
    termsW = []
    for i in range(len(vecs)):
        if i in S:
            termsW.append((tuple(vecs[i]) + tuple([0] * len(free)), coef[i] / 2))
        else:
            e = [0] * len(free)
            e[free.index(i)] = 1
            termsW.append((tuple([0] * d) + tuple(e), coef[i] / 2))
    Wm1 = poly_from(termsW, dimW)
    mz, mw = [], []
    cz = {tuple([0] * d): Fraction(1)}
    cw = {tuple([0] * dimW): Fraction(1)}
    for j in range(1, 7):
        cz = mult(cz, Zm1)
        cw = mult(cw, Wm1)
        mz.append(const(cz, d))
        mw.append(const(cw, dimW))
    eq = [mz[j] == mw[j] for j in range(6)]
    P(f"    E(Z-1)^j == E(W-1)^j for j=1..6: {eq}  (claimed: equal for j <= L)")
    P(f"    order 6: Haar {mz[5]} vs mixed {mw[5]} (length-6 relation 5*r3 - (5r3) = 0 in the free part)")

# ---------------------------------------------------------------------------
P("")
P("(E2) PRIMECORR Theorem P1(b): padding, ordinates {1, 2, sqrt2}, weights 1")
vecs2 = [(1, 0), (2, 0), (0, 1)]
G = poly_from([(v, Fraction(1)) for v in vecs2], 2)
Gind = poly_from([(tuple(1 if j == i else 0 for j in range(3)), Fraction(1)) for i in range(3)], 3)
for k in [1, 2, 3, 4, 5, 6, 7]:
    D = const(power(G, k, 2), 2) - const(power(Gind, k, 3), 3)
    P(f"  Delta_{k} = {D}")
P("  odd relations of length <= 5 among {1,2,sqrt2}: lengths "
  f"{sorted(set(sum(map(abs, k)) for k in l1_ball(3, 5) if k[0] + 2 * k[1] == 0 and k[2] == 0 and sum(map(abs, k)) % 2 == 1))}")
P("  hand count: Delta_5 = 2*(20 + 30 + 60) = 220 (pads by rho = 1, 2, sqrt2)")


# ---------------------------------------------------------------------------
def E_ab_definition(R, I, a, b, d):
    """R, I: lists of (vec, weight) zeros (with multiplicity: repeat entries).
    Returns E_A[G_R^a G_I^b] - E_A[G_R^a] * E_ind[(G_I^ind)^b]."""
    GR = poly_from(R, d)
    GI = poly_from(I, d)
    nI = len(I)
    GIind = poly_from([(tuple(1 if j == i else 0 for j in range(nI)), w) for i, (v, w) in enumerate(I)], nI)
    joint = const(mult(power(GR, a, d), power(GI, b, d)), d)
    return joint - const(power(GR, a, d), d) * const(power(GIind, b, nI), nI)


def E_ab_sequences(R, I, a, b, d):
    """sum over sequences (a R-terms then b I-terms, signs) with kappa_I != 0 and total 0."""
    tot = Fraction(0)
    nR, nI = len(R), len(I)
    for rs in itertools.product(range(nR), repeat=a):
        for is_ in itertools.product(range(nI), repeat=b):
            for sg in itertools.product((1, -1), repeat=a + b):
                kI = [0] * nI
                for t, i in enumerate(is_):
                    kI[i] += sg[a + t]
                if not any(kI):
                    continue
                s = [0] * d
                w = Fraction(1)
                for t, r in enumerate(rs):
                    for j in range(d):
                        s[j] += sg[t] * R[r][0][j]
                    w *= R[r][1]
                for t, i in enumerate(is_):
                    for j in range(d):
                        s[j] += sg[a + t] * I[i][0][j]
                    w *= I[i][1]
                if not any(s):
                    tot += w
    return tot


P("")
P("(E3) PRIMECORR Lemma P1' (decoupled prediction), R = {1,2}, I = {x,y,z}, weights 1")
# basis (1, x, y, z)
R = [((1, 0, 0, 0), Fraction(1)), ((2, 0, 0, 0), Fraction(1))]
I = [((0, 1, 0, 0), Fraction(1)), ((0, 0, 1, 0), Fraction(1)), ((0, 0, 0, 1), Fraction(1))]
d = 4
Es = {}
for a in range(0, 6):
    for b in range(1, 7 - a):
        Es[(a, b)] = E_ab_definition(R, I, a, b, d)
P(f"  max_(a+b<=6, b>=1) |E_ab| = {max(abs(v) for v in Es.values())}  (claimed 0: relations among R-zeros only)")
chk = all(E_ab_sequences(R, I, a, b, d) == Es[(a, b)] for a in range(0, 3) for b in range(1, 4 - a))
P(f"  sequence enumeration agrees with the definition for a+b <= 3: {chk}")
Gall = poly_from(R + I, d)
Gind5 = poly_from([(tuple(1 if j == i else 0 for j in range(5)), Fraction(1)) for i in range(5)], 5)
P("  plain-moment excess Delta_k of G_R + G_I (one window), k = 3, 4, 5: " +
  ", ".join(str(const(power(Gall, k, d), d) - const(power(Gind5, k, 5), 5)) for k in [3, 4, 5]) +
  "  (file: 6, 0, 460)")
# double band zero
Idbl = [((0, 1, 0, 0), Fraction(1)), ((0, 1, 0, 0), Fraction(1)), ((0, 0, 1, 0), Fraction(1))]
P(f"  band with a double zero: E_02 = {E_ab_definition(R, Idbl, 0, 2, d)} (formula sum 2m(m-1) = 4);"
  f" by enumeration: {E_ab_sequences(R, Idbl, 0, 2, d)}")
# a relation through the band of length 3 with one R-zero: R={1}, I={x, 1+x}
R3 = [((1, 0), Fraction(1))]
I3 = [((0, 1), Fraction(1)), ((1, 1), Fraction(1))]
P(f"  R={{1}}, I={{x, 1+x}} (relation (1+x) - x - 1, height <= T*): E_12 = {E_ab_definition(R3, I3, 1, 2, 2)}"
  f" (>= 1, weight 1 per sequence; 4 sequences)")

# ---------------------------------------------------------------------------
P("")
P("(E4) PRIMECORR Sec. 3 comparison spectrum AS WRITTEN: Gamma = {n alpha} u B u {z_x, x+z_x}")
# alpha = 1; band (s1,s2) = (10.5, 21.5), g = 0.5; margins (10,11), (21,22) contain no integer.
# B = {x1, x2} at heights 12.3, 15.7 (independent generators); z_x at height T*+1.5 (weight 0).
s1, s2, g, Tstar = 10.5, 21.5, 0.5, 30.0
# basis: (alpha, x1, x2, z1, z2)
zeros = []
for n in range(1, 26):
    zeros.append(("%d*alpha" % n, float(n), (n, 0, 0, 0, 0)))
zeros += [("x1", 12.3, (0, 1, 0, 0, 0)), ("x2", 15.7, (0, 0, 1, 0, 0)),
          ("z1", Tstar + 1.5, (0, 0, 0, 1, 0)), ("z2", Tstar + 1.5, (0, 0, 0, 0, 1)),
          ("x1+z1", Tstar + 13.8, (0, 1, 0, 1, 0)), ("x2+z2", Tstar + 17.2, (0, 0, 1, 0, 1))]


def betaR(h):
    if h <= s1 - g or (s2 + g <= h <= Tstar):
        return Fraction(1)
    if h > Tstar + 1:
        return Fraction(0)
    raise ValueError("height in a margin or in (T*, T*+1]: %r" % h)


margin_hits = [nm for nm, h, v in zeros if (s1 - g < h < s1 + g) or (s2 - g < h < s2 + g)]
P(f"  (Z): zeros in the open margins: {margin_hits} (none, so (Z) holds)")
Iz = [(v, Fraction(1)) for nm, h, v in zeros if s1 < h < s2]
Rz = [(v, betaR(h)) for nm, h, v in zeros if not (s1 < h < s2)]
P(f"  band I contains: {[nm for nm, h, v in zeros if s1 < h < s2]}")
E12 = E_ab_definition(Rz, Iz, 1, 2, 5)
E04 = E_ab_definition(Rz, Iz, 0, 4, 5)
P(f"  E_12 = {E12}, E_04 = {E04}  -> (D) requires E < 1: FAILS (truncated model; the full spectrum"
  f" has E at least as large, all weights being >= 0)")

P("")
P("(E5) the same spectrum with the multiples of alpha in [s1-1, s2+1] removed")
zeros_fix = [z for z in zeros if not (z[0].endswith("alpha") and s1 - 1 <= z[1] <= s2 + 1)]
Iz = [(v, Fraction(1)) for nm, h, v in zeros_fix if s1 < h < s2]
Rz = [(v, betaR(h)) for nm, h, v in zeros_fix if not (s1 < h < s2)]
P(f"  band I contains: {[nm for nm, h, v in zeros_fix if s1 < h < s2]}")
mx = Fraction(0)
for a in range(0, 5):
    for b in range(1, 6 - a):
        if a + b <= 4 or (a + b == 5 and a <= 1):
            mx = max(mx, abs(E_ab_definition(Rz, Iz, a, b, 5)))
P(f"  max |E_ab| over a+b <= 4 and (a,b) in {{(0,5),(1,4)}}: {mx} (H_dec's (D) holds in this truncation)")
P("  GLI(psi), psi > 3 eventually, still fails: 24*alpha - 23*alpha - alpha = 0 (length 3, height 24),"
  " and in the infinite spectrum such relations occur in every gap between the sparse bands")

open(os.path.join(HERE, "r5_exact_moments_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
