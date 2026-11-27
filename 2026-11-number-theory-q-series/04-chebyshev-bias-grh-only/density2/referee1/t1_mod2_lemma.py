"""Referee test 1 (independent of the authors' scripts): Lemma 1.1 of TOY.md (mod-2 lemma).

Lemma 1.1: for G a subgroup of R and f = sum_i n_i 1_{g_i + L_i} (n_i in Z, L_i != {0}),
f is not congruent to 1_{0} mod 2.

Tests.
 A. G = Z, random configurations with subgroups qZ: check f mod 2 != 1_{0} on a window, and
    check the mechanism of the proof (super-increasing l_i in L_i, distinct subset sums,
    Delta_{l_1}...Delta_{l_r} f == 0 identically, Delta...1_{0} odd at 2^r points).
 B. G = Z + alpha Z = Z^2 (rank-2 hull): GF(2) linear algebra on a box B_R. Family = all
    cosets (meeting a small box) of a list of rank-1 subgroups Zv and of several finite-index
    sublattices. Compute the subspace of the GF(2)-span consisting of vectors supported in the
    small box B_s; Lemma 1.1 (and the stronger group-ring version: f = Q mod 2 with Q finitely
    supported forces Q = 0 mod 2) predicts it is {0}; in particular 1_{0} is not in the span.
    The margin R - s exceeds the sum of the chosen l_j, so a window solution would contradict
    the proof of the lemma, i.e. the test has teeth.
 C. Control: finite groups (Z/p)^2 with p odd, where the analogue of the lemma FAILS
    (1_{0} = 1_G + sum of the p+1 lines through 0, mod 2), showing that the test detects
    solutions when they exist and that torsion-freeness is essential.
"""
import itertools
import math
import random

import numpy as np

rng = random.Random(20261126)
out = []


def log(s=""):
    print(s)
    out.append(s)


# ---------------------------------------------------------------- A. G = Z
log("A. G = Z, random configurations f = sum n_i 1_{g_i + q_i Z}")
n_trials = 20000
viol = 0
mech_fail = 0
for trial in range(n_trials):
    r = rng.randint(1, 6)
    qs = [rng.randint(1, 12) for _ in range(r)]
    gs = [rng.randrange(q) for q in qs]
    ns = [rng.choice([-3, -2, -1, 1, 2, 3]) for _ in range(r)]
    Q = 1
    for q in qs:
        Q = Q * q // math.gcd(Q, q)
    W = 2 * Q + 5
    xs = np.arange(-W, W + 1)
    f = np.zeros_like(xs)
    for n, g, q in zip(ns, gs, qs):
        f += n * (((xs - g) % q) == 0)
    target = (xs == 0).astype(int)
    if np.all((f - target) % 2 == 0):
        viol += 1
    # mechanism of the proof
    ls = []
    for q in qs:
        bound = 2 * sum(abs(l) for l in ls)
        k = bound // q + 1
        ls.append(k * q * rng.choice([-1, 1]))
    sums = [sum(l for l, b in zip(ls, bits) if b) for bits in itertools.product([0, 1], repeat=r)]
    if len(set(sums)) != 2 ** r:
        mech_fail += 1
        continue

    def fval(x):
        return sum(n for n, g, q in zip(ns, gs, qs) if (x - g) % q == 0)

    def delta_all(func, x):
        tot = 0
        for bits in itertools.product([0, 1], repeat=r):
            s = sum(l for l, b in zip(ls, bits) if b)
            tot += (-1) ** (r - sum(bits)) * func(x + s)
        return tot

    for x in [rng.randint(-50, 50) for _ in range(3)]:
        if delta_all(fval, x) != 0:
            mech_fail += 1
    odd_pts = sum(1 for s in sums if delta_all(lambda y: int(y == 0), -s) % 2 != 0)
    if odd_pts != 2 ** r:
        mech_fail += 1
log(f"  trials = {n_trials}; configurations with f = 1_0 mod 2 on the window: {viol}")
log(f"  failures of the proof mechanism (distinct subset sums, Delta f = 0, Delta 1_0 odd at 2^r pts): {mech_fail}")


# ---------------------------------------------------------------- GF(2) helpers
def reduce_family(vectors, outside_mask):
    """Forward elimination over GF(2) with pivots only on 'outside' bits.
    Returns (number of pivot rows, list of nonzero inside-only remainders)."""
    pivots = {}  # leading outside bit -> row
    inside = []
    for v in vectors:
        while True:
            vo = v & outside_mask
            if vo == 0:
                break
            lead = vo.bit_length() - 1
            if lead in pivots:
                v ^= pivots[lead]
            else:
                pivots[lead] = v
                break
        if (v & outside_mask) == 0 and v != 0:
            # inside-only remainder; reduce it against previous inside remainders
            for u in inside:
                if v ^ u < v:
                    v ^= u
            if v:
                inside.append(v)
                inside.sort(reverse=True)
    return len(pivots), inside


# ---------------------------------------------------------------- B. G = Z^2
log()
log("B. G = Z + alpha Z (= Z^2), GF(2) span of coset indicators on a box")
R, s, rsmall = 34, 3, 3
pts = [(x, y) for x in range(-R, R + 1) for y in range(-R, R + 1)]
idx = {p: i for i, p in enumerate(pts)}
# inside = B_s occupies the LOWEST bits; outside = everything else
inside_pts = [p for p in pts if max(abs(p[0]), abs(p[1])) <= s]
order = inside_pts + [p for p in pts if max(abs(p[0]), abs(p[1])) > s]
bit = {p: i for i, p in enumerate(order)}
n_in = len(inside_pts)
outside_mask = ((1 << len(order)) - 1) ^ ((1 << n_in) - 1)


def vec(points):
    v = 0
    for p in points:
        if p in bit:
            v |= 1 << bit[p]
    return v


dirs = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1), (1, -2), (2, -1), (2, 0), (0, 2), (2, 2)]
lattices = {  # finite-index sublattices given by a membership test and coset representatives
    "2Z x Z": (lambda x, y: x % 2 == 0, [(0, 0), (1, 0)]),
    "Z x 2Z": (lambda x, y: y % 2 == 0, [(0, 0), (0, 1)]),
    "x=y mod 2": (lambda x, y: (x - y) % 2 == 0, [(0, 0), (1, 0)]),
    "2Z x 2Z": (lambda x, y: x % 2 == 0 and y % 2 == 0, [(0, 0), (1, 0), (0, 1), (1, 1)]),
    "3Z x Z": (lambda x, y: x % 3 == 0, [(0, 0), (1, 0), (2, 0)]),
    "x+y=0 mod 3": (lambda x, y: (x + y) % 3 == 0, [(0, 0), (1, 0), (2, 0)]),
    "Z^2": (lambda x, y: True, [(0, 0)]),
}
family = []
small = [(x, y) for x in range(-rsmall, rsmall + 1) for y in range(-rsmall, rsmall + 1)]
for v in dirs:
    seen = set()
    for p in small:
        # canonical label of the line p + Z v
        g = math.gcd(abs(v[0]), abs(v[1]))
        vp = (v[0] // g, v[1] // g)
        label = (p[0] * vp[1] - p[1] * vp[0], (p[0] * vp[0] + p[1] * vp[1]) % g if g > 1 else 0)
        # for non-primitive v the coset p + Zv is a sub-progression; label by residue along the line
        key = (label[0], (p[0] // vp[0] if vp[0] else p[1] // vp[1]) % g if g > 1 else 0)
        if key in seen:
            continue
        seen.add(key)
        line = [(p[0] + k * v[0], p[1] + k * v[1]) for k in range(-3 * R, 3 * R + 1)]
        family.append(vec(line))
for name, (test, reps) in lattices.items():
    for (a, b) in reps:
        family.append(vec([(x, y) for (x, y) in pts if test(x - a, y - b)]))
npiv, inside = reduce_family(family, outside_mask)
# one l_j per subgroup: l = v for the lines; for the lattices (0,1),(1,0),(1,1),(2,0),(0,1),(1,-1),(1,0)
sum_l = sum(max(abs(v[0]), abs(v[1])) for v in dirs) + (1 + 1 + 1 + 2 + 1 + 1 + 1)
log(f"  box B_R, R = {R} ({len(order)} points); inside box B_s, s = {s}; family size = {len(family)}")
log(f"  margin R - s = {R - s} >= sum over subgroups of |l_j|_inf = {sum_l} (so window solutions would be genuine)")
log(f"  rank of span restricted to outside coordinates = {npiv}")
log(f"  dimension of (span) intersect (functions supported in B_s) = {len(inside)}")
# direct check that 1_0 is not in the span (all coordinates)
target = vec([(0, 0)])
_, inside_with_target = reduce_family(family + [target], outside_mask)
log(f"  adding 1_0 to the family raises that dimension to {len(inside_with_target)} (=1 means 1_0 NOT in span)")

# ---------------------------------------------------------------- C. control: (Z/p)^2
log()
log("C. Control: finite groups (Z/p)^2 (torsion), all cosets of all nontrivial subgroups")
for p in [2, 3, 5]:
    elems = [(x, y) for x in range(p) for y in range(p)]
    bitp = {e: i for i, e in enumerate(elems)}

    def vp_(points):
        v = 0
        for e in points:
            v |= 1 << bitp[e]
        return v

    fam = []
    # subgroups of order p (lines through 0) and the whole group
    lines0 = []
    for d in [(1, 0)] + [(k, 1) for k in range(p)]:
        L = {((k * d[0]) % p, (k * d[1]) % p) for k in range(p)}
        lines0.append(L)
        for e in elems:
            fam.append(vp_({((a + e[0]) % p, (b + e[1]) % p) for (a, b) in L}))
    fam.append(vp_(elems))
    # rank test: is 1_0 in span(fam)?
    def rank(vs):
        piv = {}
        for v in vs:
            while v:
                lead = v.bit_length() - 1
                if lead in piv:
                    v ^= piv[lead]
                else:
                    piv[lead] = v
                    break
        return len(piv)
    r0 = rank(fam)
    r1 = rank(fam + [vp_([(0, 0)])])
    log(f"  p = {p}: 1_0 in GF(2)-span of cosets of nontrivial subgroups? {'YES' if r0 == r1 else 'no'}")
    if p % 2 == 1:
        # explicit identity: 1_G + sum of the p+1 lines through 0 == 1_0 mod 2
        tot = np.zeros((p, p), dtype=int)
        tot += 1
        for L in lines0:
            for (a, b) in L:
                tot[a, b] += 1
        ok = np.all((tot % 2) == np.array([[1 if (a, b) == (0, 0) else 0 for b in range(p)] for a in range(p)]))
        log(f"         explicit identity 1_G + sum_(p+1 lines through 0) = 1_0 (mod 2): {ok}")

with open(__file__.replace(".py", "_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
