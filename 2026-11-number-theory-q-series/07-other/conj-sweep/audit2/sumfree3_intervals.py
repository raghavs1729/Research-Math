"""Referee: mechanical check of every interval stated in the proof of the Theorem of R06, for d = 2..150.
(I)  the 8 residue types mod m = 8d-3 (subsets of {a,b,e} completed by elements of U), min and max of sigma with
     DISTINCT underlying elements, compared with the stated ranges; disjointness from [2d-2,3d-3] u [10d-5,11d-6].
(II) the sumsets used for the coverings, computed as sets, compared with the stated intervals; the unions."""
import itertools, json
def iv(lo, hi): return set(range(lo, hi + 1))
res = {'I_ok_from': None, 'II2_ok': [], 'II3_ok': []}
bad_I = []; bad_claims = []; cover2 = {}; cover3 = {}
for d in range(2, 151):
    m = 8 * d - 3
    C = list(range(2 * d - 1, 3 * d - 1)); E = 4 * d - 1
    assert len(C) == d and C[-1] == 3 * d - 2
    # --- (I): elements with their residues; block elements: residues [2d-2, 3d-3], arbitrarily many per residue
    # model: U-elements = C (residue = value) and two copies of each block residue (two blocks suffice for distinctness)
    Uel = [('C', c, c) for c in C] + [('B', k, r) for k in (0, 1, 2) for r in range(2 * d - 2, 3 * d - 2)]
    target = iv(2 * d - 2, 3 * d - 3) | iv(10 * d - 5, 11 * d - 6)
    stated = {(): (6 * d - 6, 9 * d - 6), ('a',): (4 * d - 3, 6 * d - 3), ('b',): (5 * d - 5, 7 * d - 5), ('e',): (8 * d - 5, 10 * d - 6),
              ('a', 'b'): (3 * d - 2, 4 * d - 2), ('a', 'e'): (6 * d - 2, 7 * d - 2), ('b', 'e'): (7 * d - 4, 8 * d - 4), ('a', 'b', 'e'): (5 * d - 1, 5 * d - 1)}
    val = {'a': 1, 'b': d - 1, 'e': 4 * d - 1}
    okI = True
    for typ, (lo, hi) in stated.items():
        base = sum(val[x] for x in typ); k = 3 - len(typ)
        sums = {base + sum(u[2] for u in comb) for comb in itertools.combinations(Uel, k)} if d <= 14 else None
        if sums is not None:
            if not (min(sums) >= lo and max(sums) <= hi): bad_claims.append((d, typ, min(sums), max(sums), lo, hi))
            if sums & target: okI = False
        if iv(lo, hi) & target: okI = False
        if lo <= 0 or hi >= 2 * m: okI = False
    if not okI: bad_I.append(d)
    # --- (II)
    two = {x + y for x, y in itertools.combinations(C, 2)}; three = {x + y + z for x, y, z in itertools.combinations(C, 3)} if d >= 3 else set()
    c_two = (two == iv(4 * d - 1, 6 * d - 5)); c_three = (three == iv(6 * d, 9 * d - 9)) if d >= 3 else None
    S1 = {1 + (d - 1) + c for c in C}
    cov = {1 + s for s in two} | {d - 1 + s for s in two} | three | {E + s for s in two}
    cover2[d] = (iv(3 * d - 1, 4 * d - 2) <= S1) and (iv(4 * d, 10 * d - 6) <= cov)
    rho = range(2 * d - 2, 3 * d - 2)
    pre = [1, d - 1] + C + [E]
    gap = {r + s1 + s2 for r in rho for s1, s2 in itertools.combinations(pre, 2)} if d >= 3 else set()
    cover3[d] = iv(3 * d - 2, 10 * d - 6) <= gap
    if d in (4, 5, 6): print(f'd={d}: [4d,10d-6] minus the four families of (II):', sorted(iv(4 * d, 10 * d - 6) - cov), '; two-sums = [4d-1,6d-5]:', c_two, '; three-sums = [6d,9d-9]:', c_three)
print('(I): stated ranges contain the true ranges (brute force, d <= 14): violations', bad_claims)
print('(I): d in 2..150 for which some type meets [2d-2,3d-3] u [10d-5,11d-6] or leaves (0,2m):', bad_I)
print('(II) coverings of [3d-1,4d-2] and [4d,10d-6] hold exactly for d in:', [d for d in cover2 if cover2[d]][:3], '...; fail for d =', [d for d in cover2 if not cover2[d]])
print('(II) covering of the gap residues [3d-2,10d-6] by rho + s1 + s2 fails for d =', [d for d in cover3 if not cover3[d]])
