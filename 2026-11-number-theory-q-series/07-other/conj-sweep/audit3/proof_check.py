# AUDIT-3: line-by-line mechanical check of the proof of Theorem 2 in R06-3-sumfree.md.
#   Part A: every inequality of the proof as a linear form L(f,d) >= 0; exact proof on the region
#           R = {f >= 1, d >= 2f+1, d >= 6} (vertices + recession directions), and which hypothesis it needs.
#   Part B: every interval claim of (I) and (II) as an equality/inclusion of explicit finite sets, f <= FM, d <= DM.
#   Part C: independent count of representations x = a+b+c (a<b<c in T) by inclusion-exclusion with FFT.
#   Part D: minimal (k, p, m) of T in the sense of the source (characteristic sequence).
import sys, itertools
from fractions import Fraction as Fr
import numpy as np

FM, DM = int(sys.argv[1]), int(sys.argv[2])

# ---------------------------------------------------------------- Part A
# L = (alpha, beta, gamma) means alpha*d + beta*f + gamma >= 0.  "a < b" for integers is b - a - 1 >= 0.
def lin(expr):  # expr: dict with keys 'd','f','1'
    return (expr.get('d', 0), expr.get('f', 0), expr.get('1', 0))
def L(d=0, f=0, c=0): return (d, f, c)
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def lt(a, b):  # a < b  (integers)
    return sub(sub(b, a), (0, 0, 1))
def le(a, b): return sub(b, a)
INEQ = [
 # name, linear form, where used
 ("seeds: f < g = d-f",                                   lt(L(0,1,0), L(1,-1,0))),
 ("seeds: g < h = 2d-f",                                  lt(L(1,-1,0), L(2,-1,0))),
 ("T ordered: max C = 3d-f-1 < E = 4d-f",                 lt(L(3,-1,-1), L(4,-1,0))),
 ("T ordered: E < min B_0 = 10d-3f-2",                    lt(L(4,-1,0), L(10,-3,-2))),
 ("T ordered: max B_j < min B_{j+1}, i.e. d-1 < m",       lt(L(1,0,-1), L(8,-2,-1))),
 ("(I) C: x <= 3d-f-1 < f+g+h = 3d-f",                    lt(L(3,-1,-1), L(3,-1,0))),
 ("(I) E: max(f+g+c) = 4d-f-1 < E",                       lt(L(4,-1,-1), L(4,-1,0))),
 ("(I) E: E < f + (2d-f) + (2d-f+1) = 4d-f+1",            lt(L(4,-1,0), L(4,-1,1))),
 ("(I) E: f <= g (least sum with two elements of C uses f)", le(L(0,1,0), L(1,-1,0))),
 ("(I) residues: e = 4d-f < m (all residues are the elements themselves)", lt(L(4,-1,0), L(8,-2,-1))),
 ("(I) residues: a = f < 2d-f-1 (f is not a U-residue)",  lt(L(0,1,0), L(2,-1,-1))),
 ("(I) residues: b = d-f < 2d-f-1",                       lt(L(1,-1,0), L(2,-1,-1))),
 ("(I) sigma < 2m: max sigma = 10d-3f-3 < 16d-4f-2",      lt(L(10,-3,-3), L(16,-4,-2))),
 ("(I) UUU lower: 3d-f-2 < 6d-3f-3",                      lt(L(3,-1,-2), L(6,-3,-3))),
 ("(I) UUU upper (printed 9d-3f-3) < 10d-3f-2",           lt(L(9,-3,-3), L(10,-3,-2))),
 ("(I) aUU lower: 3d-f-2 < 4d-f-2",                       lt(L(3,-1,-2), L(4,-1,-2))),
 ("(I) aUU upper (printed 6d-f-2) < 10d-3f-2",            lt(L(6,-1,-2), L(10,-3,-2))),
 ("(I) bUU lower: 3d-f-2 < 5d-3f-2",                      lt(L(3,-1,-2), L(5,-3,-2))),
 ("(I) bUU upper (printed 7d-3f-2) < 10d-3f-2",           lt(L(7,-3,-2), L(10,-3,-2))),
 ("(I) eUU lower: 3d-f-2 < 8d-3f-2",                      lt(L(3,-1,-2), L(8,-3,-2))),
 ("(I) eUU upper: 10d-3f-3 < 10d-3f-2 (needs distinctness)", lt(L(10,-3,-3), L(10,-3,-2))),
 ("(I) abU lower: 3d-f-2 < 3d-f-1 (tight)",               lt(L(3,-1,-2), L(3,-1,-1))),
 ("(I) abU upper: 4d-f-1 < 10d-3f-2",                     lt(L(4,-1,-1), L(10,-3,-2))),
 ("(I) aeU lower: 3d-f-2 < 6d-f-1",                       lt(L(3,-1,-2), L(6,-1,-1))),
 ("(I) aeU upper: 7d-f-1 < 10d-3f-2",                     lt(L(7,-1,-1), L(10,-3,-2))),
 ("(I) beU lower: 3d-f-2 < 7d-3f-1",                      lt(L(3,-1,-2), L(7,-3,-1))),
 ("(I) beU upper: 8d-3f-1 < 10d-3f-2",                    lt(L(8,-3,-1), L(10,-3,-2))),
 ("(I) abe: 3d-f-2 < 5d-f",                               lt(L(3,-1,-2), L(5,-1,0))),
 ("(I) abe: 5d-f < 10d-3f-2",                             lt(L(5,-1,0), L(10,-3,-2))),
 ("(II) C has >= 2 elements (two-sums c1+c2 exist)",      le(L(0,0,2), L(1,0,0))),
 ("(II) C has >= 3 elements (three-sums exist)",          le(L(0,0,3), L(1,0,0))),
 ("(II) x in [4d-f+1, .]: I2 starts by end(I1)+1: 5d-3f+1 <= 6d-f-2",  le(L(5,-3,1), L(6,-1,-2))),
 ("(II) I2 ends not before I1: 6d-f-3 <= 7d-3f-3",        le(L(6,-1,-3), L(7,-3,-3))),
 ("(II) I3 starts by end(I2)+1: 6d-3f+3 <= 7d-3f-2",      le(L(6,-3,3), L(7,-3,-2))),
 ("(II) I4 starts by end(I3)+1: 8d-3f+1 <= 9d-3f-5",      le(L(8,-3,1), L(9,-3,-5))),
 ("(II) summands <= E < x for x >= 4d-f+1",               lt(L(4,-1,0), L(4,-1,1))),
 ("(II) gap: J2 = [4d-f-1,..] starts at end(J1)+1 = 4d-f-1", le(L(4,-1,-1), L(4,-1,-1))),
 ("(II) gap: J3 starts by end(J2)+1: 5d-3f-1 <= 6d-f-2",  le(L(5,-3,-1), L(6,-1,-2))),
 ("(II) gap: J3 ends not before J2: 6d-f-3 <= 7d-3f-3",   le(L(6,-1,-3), L(7,-3,-3))),
 ("(II) gap: J4 starts by end(J3)+1: 6d-3f <= 7d-3f-2",   le(L(6,-3,0), L(7,-3,-2))),
 ("(II) gap: J5 starts by end(J4)+1: 8d-3f-1 <= 9d-3f-4", le(L(8,-3,-1), L(9,-3,-4))),
 ("(II) gap: small summands <= E < min B_0 <= beta",      lt(L(4,-1,0), L(10,-3,-2))),
 ("minimality: E+m = 12d-3f-1 > max B_0 = 11d-3f-3",      lt(L(11,-3,-3), L(12,-3,-1))),
 ("minimality: E+m < min B_1 = 18d-5f-3",                 lt(L(12,-3,-1), L(18,-5,-3))),
 ("minimality: 7d-2f >= 2 (period block has one entry != 1)", le(L(0,0,2), L(7,-2,0))),
 ("minimality: d-1 >= 1",                                 le(L(0,0,2), L(1,0,0))),
]
def val(Lf, f, d): return Lf[0]*d + Lf[1]*f + Lf[2]
def holds_on_region(Lf):
    # R = {f >= 1, d >= 6, d >= 2f+1}: vertices (f,d) = (1,6), (5/2,6); recession directions (0,1), (1,2)
    verts = [(Fr(1), Fr(6)), (Fr(5, 2), Fr(6))]
    rays = [(0, 1), (1, 2)]
    return all(val(Lf, f, d) >= 0 for f, d in verts) and all(Lf[0]*dd + Lf[1]*ff >= 0 for ff, dd in rays)
def hypothesis(Lf, FG=150, DG=500):
    # brute force on the grid 1 <= f <= FG, 1 <= d <= DG: a weakest hypothesis from a fixed ladder under which L >= 0
    pts = [(f, d) for f in range(1, FG+1) for d in range(1, DG+1)]
    for k in range(1, 9):
        if all(val(Lf, f, d) >= 0 for f, d in pts if d >= k): return "none" if k == 1 else f"d >= {k}"
    if all(val(Lf, f, d) >= 0 for f, d in pts if d >= f+1): return "d >= f+1"
    if all(val(Lf, f, d) >= 0 for f, d in pts if d >= 2*f): return "d >= 2f"
    if all(val(Lf, f, d) >= 0 for f, d in pts if d >= 2*f+1): return "d >= 2f+1"
    for k in range(2, 12):
        if all(val(Lf, f, d) >= 0 for f, d in pts if d >= max(k, 2*f+1)): return f"d >= 2f+1 and d >= {k}"
    return "FAILS"
def explicit(Lf):
    a, b, c = Lf
    if a == 0 and b == 0: return f"{c} >= 0 (identity)"
    return f"{a}d >= {-b}f + {-c}".replace("+ -", "- ")
print("== Part A: inequalities (alpha*d + beta*f + gamma >= 0) ==")
allA = True
for name, Lf in INEQ:
    ok = holds_on_region(Lf); allA &= ok
    fails = [(f, d) for f in range(1, 40) for d in range(max(6, 2*f+1), 200) if val(Lf, f, d) < 0]
    print(f"  [{'ok' if ok and not fails else 'FAIL'}] {name}:  <=> {explicit(Lf)};  implied by: {hypothesis(Lf)}")
print("  all inequalities hold on R (exact, vertices + rays):", allA)

# ---------------------------------------------------------------- Part B
def iv(a, b): return set(range(a, b+1))
def sumset(A, B): return {a + b for a in A for b in B}
def check_pair(f, d, deep=False):
    """returns list of failed claim names (empty = all claims of the proof hold as sets)"""
    bad = []
    g, h, m = d - f, 2*d - f, 8*d - 2*f - 1
    C = iv(2*d-f, 3*d-f-1); E = 4*d - f
    B0 = iv(10*d-3*f-2, 11*d-3*f-3)
    # (I) x in C, x = E
    small = sorted({f, g} | C)
    if len(small) != d + 2: bad.append("f, g, C not distinct")
    if f + g + h != 3*d - f or not all(x < f + g + h for x in C): bad.append("(I) C")
    sums_E = {a+b+c for a, b, c in itertools.combinations(small, 3)} if d <= 40 else None
    if sums_E is not None and E in sums_E: bad.append("(I) E")
    if {f+g+c for c in C} != iv(3*d-f, 4*d-f-1): bad.append("(I) E: f+g+c interval")
    two_with_C_min = min(small[0] + min(C) + (min(C)+1), 3*min(C)+3)
    if two_with_C_min != 4*d - f + 1: bad.append("(I) E: least sum with two elements of C")
    # (I) x in B_j: residues
    if {x % m for x in B0} != iv(2*d-f-1, 3*d-f-2): bad.append("(I) block residues")
    if not all(x < m for x in [f, g, E] + sorted(C)): bad.append("(I) residues not reduced")
    Uinf = iv(2*d-f-1, 3*d-f-2)          # residues carried by arbitrarily many elements (blocks; and C for r >= 2d-f)
    top = 3*d - f - 1                    # residue carried by exactly one element
    # exact sets of residue sums over triples of DISTINCT elements, by type
    U2 = sumset(Uinf, Uinf) | {top + r for r in Uinf}                 # two distinct U-elements
    U3 = sumset(sumset(Uinf, Uinf), Uinf) | {top + s for s in sumset(Uinf, Uinf)}
    U1 = Uinf | {top}
    a, b, e = f, g, E
    exact = {"UUU": U3, "aUU": {a+s for s in U2}, "bUU": {b+s for s in U2}, "eUU": {e+s for s in U2},
             "abU": {a+b+u for u in U1}, "aeU": {a+e+u for u in U1}, "beU": {b+e+u for u in U1}, "abe": {a+b+e}}
    printed = {"UUU": (6*d-3*f-3, 9*d-3*f-3), "aUU": (4*d-f-2, 6*d-f-2), "bUU": (5*d-3*f-2, 7*d-3*f-2),
               "eUU": (8*d-3*f-2, 10*d-3*f-3), "abU": (3*d-f-1, 4*d-f-1), "aeU": (6*d-f-1, 7*d-f-1),
               "beU": (7*d-3*f-1, 8*d-3*f-1), "abe": (5*d-f, 5*d-f)}
    targets = iv(2*d-f-1, 3*d-f-2) | iv(10*d-3*f-2, 11*d-3*f-3)
    if targets != {r for r in range(1, 2*m) if r % m in iv(2*d-f-1, 3*d-f-2)}: bad.append("(I) target set for sigma in (0,2m)")
    slack = {}
    for t in exact:
        lo, hi = printed[t]
        if not exact[t] <= iv(lo, hi): bad.append(f"(I) type {t}: exact set not inside printed interval")
        if min(exact[t]) != lo: bad.append(f"(I) type {t}: printed lower end not attained")
        slack[t] = hi - max(exact[t])
        if exact[t] & targets: bad.append(f"(I) type {t}: residue sum meets a target interval")
        if not (3*d-f-2 < lo and hi < 10*d-3*f-2): bad.append(f"(I) type {t}: printed interval not strictly between")
        if not (0 < min(exact[t]) and max(exact[t]) < 2*m): bad.append(f"(I) type {t}: sigma not in (0,2m)")
    # (II)
    C2 = {x + y for x, y in itertools.combinations(sorted(C), 2)}
    if C2 != iv(4*d-2*f+1, 6*d-2*f-3): bad.append("(II) two-sums of C")
    C3 = {s + c for s in C2 for c in C}  # superset of 3-sums of distinct elements; exact 3-sums computed for d <= 40
    if d <= 40:
        C3 = {x + y + z for x, y, z in itertools.combinations(sorted(C), 3)}
    else:  # c1 < c2 < c3: c3 + (two-sums of elements below c3)
        C3 = set()
        cs = sorted(C)
        for i, c3 in enumerate(cs):
            if i >= 2: C3 |= {c3 + s for s in range(2*cs[0]+1, 2*cs[i-1])}   # two-sums of the interval [cs[0], cs[i-1]]
    if C3 != iv(6*d-3*f+3, 9*d-3*f-6): bad.append("(II) three-sums of C")
    if {f + g + c for c in C} != iv(3*d-f, 4*d-f-1): bad.append("(II) [3d-f, 4d-f-1]")
    I1, I2, I3, I4 = {f+s for s in C2}, {g+s for s in C2}, C3, {E+s for s in C2}
    if I1 != iv(4*d-f+1, 6*d-f-3): bad.append("(II) f+c1+c2")
    if I2 != iv(5*d-3*f+1, 7*d-3*f-3): bad.append("(II) g+c1+c2")
    if I4 != iv(8*d-3*f+1, 10*d-3*f-3): bad.append("(II) E+c1+c2")
    if (I1 | I2 | I3 | I4) != iv(4*d-f+1, 10*d-3*f-3): bad.append("(II) covering of [4d-f+1, 10d-3f-3]")
    if not (min(I4) > E): bad.append("(II) E used for some x <= E")
    # gap after B_j (shift by (j+1)m): y = rho + s1 + s2
    R = iv(2*d-f-1, 3*d-f-2)
    if {x - m for x in B0} != R: bad.append("(II) gap: rho range")
    J1, J2, J3 = {r + f + g for r in R}, {r + f + c for r in R for c in C}, {r + g + c for r in R for c in C}
    J4, J5 = {r + s for r in R for s in C2}, {r + E + c for r in R for c in C}
    if J1 != iv(3*d-f-1, 4*d-f-2): bad.append("(II) gap f+g")
    if J2 != iv(4*d-f-1, 6*d-f-3): bad.append("(II) gap f+c")
    if J3 != iv(5*d-3*f-1, 7*d-3*f-3): bad.append("(II) gap g+c")
    if J4 != iv(6*d-3*f, 9*d-3*f-5): bad.append("(II) gap c1+c2")
    if J5 != iv(8*d-3*f-1, 10*d-3*f-3): bad.append("(II) gap E+c")
    gap0 = iv(11*d-3*f-2, 10*d-3*f-3+m)
    if {x - m for x in gap0} != iv(3*d-f-1, 10*d-3*f-3): bad.append("(II) gap: y range")
    if (J1 | J2 | J3 | J4 | J5) != iv(3*d-f-1, 10*d-3*f-3): bad.append("(II) gap covering")
    if not (E < min(B0) and max(B0) < min(gap0)): bad.append("(II) gap: summands not smaller / not distinct")
    # exhaustion of the complement of T above h, on (h, min B_2)
    X = 10*d-3*f-2 + 2*m - 1
    T = {f, g} | C | {E} | B0 | {x + m for x in B0}
    pieces = [C, iv(3*d-f, 4*d-f-1), {E}, iv(4*d-f+1, 10*d-3*f-3), B0, gap0, {x+m for x in B0}, {x+m for x in gap0}]
    union = set().union(*pieces)
    if union != iv(h, 10*d-3*f-3 + 2*m) or sum(len(p) for p in pieces) != len(union):   # consecutive, disjoint, from h to the end of the 2nd gap
        bad.append("exhaustion: pieces do not partition")
    if {x for x in union if x in T} != set().union(pieces[0], pieces[2], pieces[4], pieces[6]): bad.append("exhaustion: T-pieces")
    # deep: validate the multiplicity model of (I) against actual distinct elements of T (three blocks)
    if deep:
        Tel = [(f, 'a'), (g, 'b'), (E, 'e')] + [(c, 'U') for c in C] + [(x + j*m, 'U') for j in range(3) for x in B0]
        got = {}
        for (x, tx), (y, ty), (z, tz) in itertools.combinations(Tel, 3):
            key = ''.join(sorted(tx + ty + tz, key='abeU'.index))
            got.setdefault(key, set()).add(x % m + y % m + z % m)
        for t in exact:
            if got.get(t, set()) != exact[t]: bad.append(f"(I) type {t}: multiplicity model != brute force over elements")
        # and the statement itself: no element of B_0..B_2 is a sum of three distinct elements of T
        els = sorted(x for x, _ in Tel); es = set(els)
        s3 = {x + y + z for x, y, z in itertools.combinations(els, 3)}
        if s3 & {x + j*m for j in range(3) for x in B0}: bad.append("(I) blocks: brute force finds a representation")
    return bad, slack

print("== Part B: interval claims as finite sets ==")
nB = 0; failB = {}; slacks = {}
for f in range(1, FM+1):
    for d in range(max(6, 2*f+1), DM+1):
        bad, slack = check_pair(f, d, deep=(d <= 22))
        nB += 1
        if bad: failB[(f, d)] = bad
        for t, s in slack.items(): slacks.setdefault(t, set()).add(s)
print(f"  pairs f <= {FM}, max(6,2f+1) <= d <= {DM}: {nB}; pairs with a failed claim: {len(failB)}", dict(list(failB.items())[:5]))
print("  printed upper end minus true maximum of the residue sum, by type (0 = exact):", {t: sorted(s) for t, s in slacks.items()})
print("  below the range (2f+1 <= d < 6), failed claims:")
for f in range(1, 3):
    for d in range(2*f+1, 6):
        bad, _ = check_pair(f, d, deep=True)
        print(f"    (f,d)=({f},{d}):", bad)

# ---------------------------------------------------------------- Part C
def reps_fft(Tset, X):
    """r[x] = number of {a<b<c} in T with a+b+c = x, 0 <= x <= X, by inclusion-exclusion on ordered triples"""
    n = 1
    while n < 3*X + 4: n *= 2
    ind = np.zeros(n);
    for t in Tset:
        if t <= X: ind[t] = 1.0
    F = np.fft.rfft(ind)
    A3 = np.fft.irfft(F**3, n)
    ind2 = np.zeros(n)
    for t in Tset:
        if 2*t <= X: ind2[2*t] = 1.0
    A21 = np.fft.irfft(np.fft.rfft(ind2) * F, n)          # ordered (a,a,b), b arbitrary (includes b = a)
    A3, A21 = A3[:X+1], A21[:X+1]                         # only sums <= X are complete
    A3r, A21r = np.rint(A3), np.rint(A21)
    assert np.abs(A3 - A3r).max() < 0.2 and np.abs(A21 - A21r).max() < 0.2
    A111 = np.zeros(X+1)
    for t in Tset:
        if 3*t <= X: A111[3*t] = 1.0
    # ordered triples with all distinct = A3 - 3*A21 + 2*A111
    r = (A3r - 3*A21r + 2*A111) / 6.0
    rr = np.rint(r); assert np.abs(r - rr).max() < 1e-9
    assert rr.min() >= 0
    return rr.astype(np.int64)
print("== Part C: representation counts by FFT inclusion-exclusion (4 periods) ==")
nC = 0; badC = []; minrep = None; argmin = []
for f in range(1, FM+1):
    for d in range(max(6, 2*f+1), DM+1):
        m = 8*d-2*f-1; h = 2*d - f; X = 10*d-3*f-2 + 4*m - 1
        T = {f, d-f} | set(range(2*d-f, 3*d-f)) | {4*d-f}
        for j in range(4): T |= set(range(10*d-3*f-2+j*m, 11*d-3*f-2+j*m))
        r = reps_fft(T, X); nC += 1
        inT = np.zeros(X+1, dtype=bool); inT[sorted(T)] = True
        xs = np.arange(X+1)
        v1 = xs[(xs > h) & inT & (r > 0)]; v2 = xs[(xs > h) & ~inT & (r == 0)]
        if len(v1) or len(v2): badC.append((f, d, v1[:3].tolist(), v2[:3].tolist()))
        mr = int(r[(xs > h) & ~inT].min())
        if minrep is None or mr < minrep: minrep = mr
print(f"  pairs: {nC}; violations of (I) or (II): {len(badC)} {badC[:5]}; least number of representations of a non-element: {minrep}")
# sanity of the FFT counter against brute force
T = {2, 4} | set(range(10, 16)) | {22} | set(range(52, 58)) | set(range(95, 101))
r = reps_fft(T, 130); br = [0]*131
for a, b, c in itertools.combinations(sorted(T), 3):
    if a+b+c <= 130: br[a+b+c] += 1
print("  FFT counter equals brute force on T(2,6) up to 130:", list(r) == br)

# ---------------------------------------------------------------- Part D
def kpm(Tset, X):
    """minimal (k, p, m) of the source: m = least period of the tail of the characteristic sequence c on [1,X],
       preperiod = least N0 with c(n) = c(n+m) for all n > N0; k = #ones in [1,N0], p = #ones in a period"""
    c = np.zeros(X+1, dtype=np.int8); c[sorted(t for t in Tset if t <= X)] = 1
    tail0 = X // 2
    for m in range(1, X // 6):
        if np.array_equal(c[tail0:X+1-m], c[tail0+m:X+1]): break
    else: return None
    N0 = tail0
    while N0 >= 1 and c[N0] == c[N0+m]: N0 -= 1
    return int(c[1:N0+1].sum()), int(c[N0+1:N0+1+m].sum()), m
print("== Part D: minimal (k,p,m) of T ==")
badD = []; nD = 0
for f in range(1, FM+1):
    for d in range(max(6, 2*f+1), DM+1):
        m = 8*d-2*f-1; X = 10*d-3*f-2 + 12*m - 1
        T = {f, d-f} | set(range(2*d-f, 3*d-f)) | {4*d-f}
        for j in range(12): T |= set(range(10*d-3*f-2+j*m, 11*d-3*f-2+j*m))
        nD += 1
        if kpm(T, X) != (d+3, d, m): badD.append((f, d, kpm(T, X)))
print(f"  pairs: {nD}; (k,p,m) != (d+3, d, 8d-2f-1): {len(badD)} {badD[:5]}")
