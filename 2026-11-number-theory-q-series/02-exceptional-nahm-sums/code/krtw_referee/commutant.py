# Referee: commutant of T inside span{P^a Y^b Ssh^c Ush^d : |a|,|b|<=2, 0<=c,d<=2} at a rational point (q,x),
# acting on finite-support functions on Z^2 (exact arithmetic, flint fmpq_mat nullspace).
import sys, itertools, random
from fractions import Fraction as Fr
from flint import fmpq_mat, fmpq
exec(open('hamtest.py').read().split('random.seed')[0])
mons = [(a, b, c, d) for a in range(-2, 3) for b in range(-2, 3) for c in range(3) for d in range(3)]
def M_apply(m, G):
    a, b, c, d = m
    return lambda R, S: q ** (a * R) * Q ** (b * S) * G(R - c, S + d)
def M_dict(m, F):
    a, b, c, d = m
    return {(r + c, s - d): q ** (a * (r + c)) * Q ** (b * (s - d)) * v for (r, s), v in F.items()}
random.seed(2)
tests = [{(r0, s0): Fr(1)} for r0 in range(-3, 6) for s0 in range(-3, 6)]
pts = [(R, S) for R in range(-5, 9) for S in range(-5, 9)]
rows = []
for F in tests:
    TF = T(F)
    cols = []
    for m in mons:
        A = M_apply(m, TF); B = T(M_dict(m, F))
        cols.append([A(R, S) - B(R, S) for (R, S) in pts])
    for i in range(len(pts)):
        rows.append([cols[j][i] for j in range(len(mons))])
print("matrix", len(rows), "x", len(mons))
Mt = fmpq_mat(len(rows), len(mons), [fmpq(v.numerator, v.denominator) for row in rows for v in row])
Rr, rank = Mt.rref()
ncol = len(mons)
piv = []
for i in range(rank):
    for j in range(ncol):
        if Rr[i, j] != 0: piv.append(j); break
free = [j for j in range(ncol) if j not in piv]
nullity = len(free)
print("rank", rank, "nullity:", nullity)
class NS:
    def __init__(s): s.v = {}
    def __getitem__(s, ij): return s.v.get(ij, 0)
ns = NS()
for k, f in enumerate(free):
    ns.v[(f, k)] = 1
    for i, p in enumerate(piv): ns.v[(p, k)] = -Rr[i, f]
# also evaluate the document's H in this basis
H = {(-1, -1, 0, 0): 1, (-1, -1, 1, 0): -1 / q, (-1, 0, 0, 0): q ** 3 * x, (-1, 0, 1, 0): -q ** 2 * x,
     (1, 0, 0, 0): x / q, (1, 1, 0, 0): q ** 2 * x * x, (1, 1, 0, 1): -x * x}
vec = [H.get(m, Fr(0)) for m in mons]
res = [sum(r[j] * vec[j] for j in range(len(mons))) for r in rows]
print("document's H: number of nonzero commutator entries:", sum(1 for v in res if v != 0), "of", len(res))
for k in range(nullity):
    v = [ns[j, k] for j in range(len(mons))]
    print("null vector", k, {mons[j]: v[j] for j in range(len(mons)) if v[j] != 0})
