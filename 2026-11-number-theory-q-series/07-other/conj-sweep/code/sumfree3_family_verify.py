# Theorem 2 of R06: for f >= 1, d >= max(6, 2f+1), g = d - f:
#   S_{f,g,g+d} = {f, d-f} U [2d-f, 3d-f-1] U {4d-f} U U_j [10d-3f-2+jm, 11d-3f-3+jm],  m = 8d-2f-1.
# (1) simulation; (2) mechanical check of the interval statements (I), (II) used in the proof.
import sys
from sumfree3 import seq
FM, DM = int(sys.argv[1]), int(sys.argv[2])
bad = []; cnt = 0
for f in range(1, FM+1):
    for d in range(max(6, 2*f+1), DM+1):
        m = 8*d - 2*f - 1; nt = 12*d + 40; cnt += 1
        S = seq((f, d-f, 2*d-f), nt)
        T = [f, d-f] + list(range(2*d-f, 3*d-f)) + [4*d-f]; j = 0
        while len(T) < nt: T += list(range(10*d-3*f-2 + j*m, 11*d-3*f-2 + j*m)); j += 1
        if S != T[:nt]: bad.append((f, d))
print(f'simulation: f <= {FM}, max(6,2f+1) <= d <= {DM}: {cnt} sequences, first 12d+40 terms; mismatches: {bad[:10]}')
# exceptions just below the range
exc = []
for f in range(1, FM+1):
    for d in range(2*f+1, max(6, 2*f+1)):
        m = 8*d - 2*f - 1; nt = 12*d + 40
        S = seq((f, d-f, 2*d-f), nt)
        T = [f, d-f] + list(range(2*d-f, 3*d-f)) + [4*d-f]; j = 0
        while len(T) < nt: T += list(range(10*d-3*f-2 + j*m, 11*d-3*f-2 + j*m)); j += 1
        if S != T[:nt]: exc.append((f, d))
print('pairs with 2f+1 <= d < 6 where the pattern fails:', exc)
# (2) interval statements
def check(f, d):
    m = 8*d - 2*f - 1; g = d - f
    C = set(range(2*d-f, 3*d-f)); E = 4*d - f; small = [f, g] + sorted(C) + [E]
    Bres = set(range(2*d-f-1, 3*d-f-1))                       # residues of the blocks
    # (I) residue types: all sums of three distinct elements of T (two periods of blocks suffice for residues), via residues with multiplicity rules
    res_elems = [(f, 'a'), (g, 'b'), (E, 'e')] + [(c, 'C') for c in C] + [(r, 'B') for r in Bres]
    ok = True
    n = len(res_elems)
    targets = Bres | {r + m for r in Bres}
    for i in range(n):
        for j in range(i, n):
            for k in range(j, n):
                x, y, z = res_elems[i], res_elems[j], res_elems[k]
                # equal indices allowed only for block residues (distinct blocks can share a residue)
                if (i == j and x[1] != 'B') or (j == k and y[1] != 'B'): continue
                if x[0] + y[0] + z[0] in targets: ok = False
    # (II) coverage of the gaps by sums of three distinct smaller elements
    T1 = small
    def is3sum(x, pool):
        pool = [p for p in pool if p < x]; s = set(pool)
        for i in range(len(pool)):
            for j in range(i+1, len(pool)):
                r = x - pool[i] - pool[j]
                if r > pool[j] and r in s: return True
        return False
    for x in list(range(3*d-f, 4*d-f)) + list(range(4*d-f+1, 10*d-3*f-2)):
        if not is3sum(x, T1): ok = False
    B0 = list(range(10*d-3*f-2, 11*d-3*f-2))
    for x in range(11*d-3*f-2, 10*d-3*f-2 + m):               # gap after B_0, using B_0 and two small elements
        if not any((x - b - s1) in set(small) and (x - b - s1) != s1 for b in B0 for s1 in small): ok = False
    # elements of C above h and E are not 3-sums of smaller elements
    for x in list(range(2*d-f+1, 3*d-f)) + [E]:
        if is3sum(x, small): ok = False
    return ok
fails = [(f, d) for f in range(1, 9) for d in range(max(6, 2*f+1), 41) if not check(f, d)]
print('mechanical check of (I) and (II) for f <= 8, d <= 40: failures', fails)
