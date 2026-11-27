"""Exact certificate for the factorisation K_k(x) = B_{k+1}(x) J_k(x)  (conjecture (F_k) of proofs/KRTW-general.md, sect. 8)
for k = 2 and k = 3, by explicit chains of two elementary moves on generalized Nahm data (proofs/KRTW-F-moves.md).

A datum is (M, B, w, d):   F(M,B,w,d) = sum_{n in N^r} x^{w.n} q^{(1/2) n^T M n + B.n} / prod_i (q^{d_i};q^{d_i})_{n_i}.
All w_i >= 1, so for each power of x the sum is finite and F is an element of Q(q)[[x]].

Moves (each is an identity F(old) = F(new); proofs in the .md file, Lemma M1, M2):
  SPLIT(a)      [1/(q^d;q^d)_n = sum_{s+t=n} q^{dt}/((q^{2d};q^{2d})_s (q^{2d};q^{2d})_t)]:
                species a (d) -> a (2d) and a new last species t (2d) with the row of a, B_t = B_a + d, w_t = w_a.
  FUSE(a,b)     [Durfee: 1/((Q)_a (Q)_b) = sum_c Q^{(a-c)(b-c)}/((Q)_{a-c}(Q)_{b-c}(Q)_c), Q = q^d, d_a = d_b = d]:
                new last species c = "a+b": row_c = row_a + row_b, M_cc = M_aa + M_bb + 2 M_ab, B_c = B_a + B_b, w_c = w_a + w_b;
                afterwards M_ab += d.
  DEFUSE(a,b,c) the inverse of FUSE: allowed iff c has exactly the row/B/w of the fusion of a and b computed from
                the datum in which M_ab is lowered by d; then c is deleted and M_ab -= d.
The script (1) builds the data of K_k and of B_{k+1} (x) J_k from their definitions, (2) applies the chains, checking the
precondition of every DEFUSE exactly (Fractions), (3) checks that the two final data coincide after an explicit
permutation of the species, (4) independently evaluates all four sums (start and end of both chains) as power
series, for all powers x^0..x^D modulo q^N, by brute-force enumeration, and compares them with each other and with
the level-by-level evaluation of the original 2k-fold sum.   usage: python3 krtw_moves_certify.py [N] [D]
"""
import sys, itertools
from fractions import Fraction as Fr
from flint import fmpz_poly

# ---------------------------------------------------------------- data and moves
class Datum:
    def __init__(self, M, B, w, d, labs):
        self.M = [[Fr(v) for v in row] for row in M]; self.B = [Fr(v) for v in B]
        self.w = list(w); self.d = list(d); self.labs = list(labs)
    def copy(self):
        return Datum(self.M, self.B, self.w, self.d, self.labs)
    @property
    def r(self): return len(self.d)

def from_squares(vecs, lin, w, d, labs):
    """exponent = sum_v (v.n)^2 + lin.n"""
    n = len(d)
    M = [[Fr(0)] * n for _ in range(n)]
    for v in vecs:
        for p in range(n):
            for q in range(n):
                M[p][q] += 2 * v[p] * v[q]
    return Datum(M, lin, w, d, labs)

def split(D, a):
    D = D.copy(); r = D.r; da = D.d[a]
    row = [D.M[a][k] for k in range(r)]
    for k in range(r): D.M[k].append(row[k])
    D.M.append(row + [D.M[a][a]])
    D.B.append(D.B[a] + da); D.w.append(D.w[a]); D.d[a] = 2 * da; D.d.append(2 * da)
    D.labs.append(D.labs[a] + "''"); D.labs[a] = D.labs[a] + "'"
    return D

def fuse(D, a, b):
    assert a != b and D.d[a] == D.d[b], "FUSE needs two species of the same base"
    D = D.copy(); r = D.r; g = D.M[a][b]
    row = [D.M[a][k] + D.M[b][k] for k in range(r)]
    row[a] = D.M[a][a] + g; row[b] = D.M[b][b] + g
    cc = D.M[a][a] + D.M[b][b] + 2 * g
    D.M[a][b] += D.d[a]; D.M[b][a] += D.d[a]
    for k in range(r): D.M[k].append(row[k])
    D.M.append(row + [cc]); D.B.append(D.B[a] + D.B[b]); D.w.append(D.w[a] + D.w[b]); D.d.append(D.d[a])
    D.labs.append("[%s+%s]" % (D.labs[a], D.labs[b]))
    return D

def defuse(D, a, b, c):
    """exact check that c is the fusion of a and b, then remove it"""
    assert len({a, b, c}) == 3 and D.d[a] == D.d[b] == D.d[c], "DEFUSE: bases differ"
    r = D.r; dd = D.d[a]; g = D.M[a][b] - dd
    ok = (D.B[c] == D.B[a] + D.B[b] and D.w[c] == D.w[a] + D.w[b]
          and D.M[c][c] == D.M[a][a] + D.M[b][b] + 2 * g
          and D.M[a][c] == D.M[a][a] + g and D.M[b][c] == D.M[b][b] + g
          and all(D.M[c][k] == D.M[a][k] + D.M[b][k] for k in range(r) if k not in (a, b, c)))
    assert ok, "DEFUSE precondition fails for (%s,%s,%s)" % (D.labs[a], D.labs[b], D.labs[c])
    keep = [i for i in range(r) if i != c]
    E = Datum([[D.M[i][j] for j in keep] for i in keep], [D.B[i] for i in keep], [D.w[i] for i in keep],
              [D.d[i] for i in keep], [D.labs[i] for i in keep])
    ia, ib = keep.index(a), keep.index(b)
    E.M[ia][ib] -= dd; E.M[ib][ia] -= dd
    E.labs[ia] = "(%s|%s)" % (D.labs[a], D.labs[c]); E.labs[ib] = "(%s|%s)" % (D.labs[b], D.labs[c])
    return E

def tensor(D1, D2):
    n1, n2 = D1.r, D2.r
    M = [[Fr(0)] * (n1 + n2) for _ in range(n1 + n2)]
    for i in range(n1):
        for j in range(n1): M[i][j] = D1.M[i][j]
    for i in range(n2):
        for j in range(n2): M[n1 + i][n1 + j] = D2.M[i][j]
    return Datum(M, D1.B + D2.B, D1.w + D2.w, D1.d + D2.d, D1.labs + D2.labs)

def run(D, path, log):
    for mv in path:
        if mv[0] == 'p':
            log.append("   SPLIT  %s" % D.labs[mv[1]]); D = split(D, mv[1])
        elif mv[0] == 'f':
            log.append("   FUSE   %s , %s" % (D.labs[mv[1]], D.labs[mv[2]])); D = fuse(D, mv[1], mv[2])
        else:
            log.append("   DEFUSE %s  = fusion of  %s , %s   [precondition verified]" % (D.labs[mv[3]], D.labs[mv[1]], D.labs[mv[2]]))
            D = defuse(D, mv[1], mv[2], mv[3])
    return D

def same_up_to_permutation(D1, D2):
    if D1.r != D2.r: return None
    r = D1.r
    key = lambda D, i: (D.d[i], D.w[i], D.B[i], D.M[i][i], tuple(sorted(D.M[i])))
    if sorted(key(D1, i) for i in range(r)) != sorted(key(D2, i) for i in range(r)): return None
    def rec(i, perm, used):
        if i == r: return perm
        for j in range(r):
            if j in used or key(D1, i) != key(D2, j): continue
            if all(D1.M[i][t] == D2.M[j][perm[t]] for t in range(i)):
                res = rec(i + 1, perm + [j], used | {j})
                if res is not None: return res
        return None
    return rec(0, [], frozenset())

# ---------------------------------------------------------------- the four sums as data
def K_data(k):
    """K_k(x): species a_j = r_j - r_{j+1} (base q, j=1..k), b_l = s_l - s_{l-1} (base q^2, l=1..k)"""
    n = 2 * k
    def r(i): return [1 if j >= i else 0 for j in range(1, k + 1)] + [0] * k
    def s(i): return [0] * k + [1 if l <= i else 0 for l in range(1, k + 1)]
    vecs = []; w = [0] * n
    for i in range(1, k + 1):
        vecs.append([p + q for p, q in zip(r(i), s(i))]); vecs.append(s(i))
        w = [u + p + 2 * q for u, p, q in zip(w, r(i), s(i))]
    labs = ["a%d" % j for j in range(1, k + 1)] + ["b%d" % l for l in range(1, k + 1)]
    return from_squares(vecs, [0] * n, w, [1] * k + [2] * k, labs)

def B_data(m):
    """Bressoud chain B_m(x): species n_j = N_j - N_{j+1} (base q, j<m), n_m = N_m (base q^2)"""
    def Nv(i): return [1 if j >= i else 0 for j in range(1, m + 1)]
    vecs = [Nv(i) for i in range(1, m + 1)]
    w = [sum(v[j] for v in vecs) for j in range(m)]
    return from_squares(vecs, [0] * m, w, [1] * (m - 1) + [2], ["n%d" % j for j in range(1, m + 1)])

def C_data(k):
    """C_k(x) = sum_{0<=s_1<=...<=s_k} x^{sum s} q^{sum s^2}/prod (q^2;q^2)_{s_i-s_{i-1}}: species e_l = s_l - s_{l-1}"""
    def sg(i): return [1 if l <= i else 0 for l in range(1, k + 1)]
    vecs = [sg(i) for i in range(1, k + 1)]
    w = [sum(v[j] for v in vecs) for j in range(k)]
    return from_squares(vecs, [0] * k, w, [2] * k, ["e%d" % l for l in range(1, k + 1)])

def AGQ_data(j):
    """AG_j(x^2;q^2) = sum_{sg_1>=...>=sg_j>=0} x^{2 sum sg} q^{2 sum sg^2}/prod (q^2;q^2)_{sg_i-sg_{i+1}}"""
    def sg(i): return [1 if l >= i else 0 for l in range(1, j + 1)]
    vecs = []
    for i in range(1, j + 1): vecs += [sg(i), sg(i)]
    w = [sum(v[t] for v in vecs) for t in range(j)]
    return from_squares(vecs, [0] * j, w, [2] * j, ["g%d" % l for l in range(1, j + 1)])

# ---------------------------------------------------------------- exact evaluation modulo q^N, x-degrees <= Dmax
NQ = 40
def tr(p):
    return fmpz_poly(p.coeffs()[:NQ]) if p.degree() >= NQ else p
_ip = {}
def iqp(n, base):
    key = (n, base, NQ)
    if key not in _ip:
        c = [0] * NQ; c[0] = 1
        for i in range(1, n + 1):          # multiply by 1/(1-q^{base*i})
            st = base * i
            for e in range(st, NQ): c[e] += c[e - st]
        _ip[key] = fmpz_poly(c)
    return _ip[key]

def evaluate(D, Dmax):
    """coefficients of x^0..x^Dmax of F(D), each modulo q^NQ  (twice the exponent is tracked to allow half-integers)"""
    r = D.r; res = [fmpz_poly([0]) for _ in range(Dmax + 1)]
    n = [0] * r
    def rec(i, deg):
        if i == r:
            e2 = sum(n[p] * D.M[p][q] * n[q] for p in range(r) for q in range(r)) + 2 * sum(D.B[p] * n[p] for p in range(r))
            assert e2.denominator == 1 and e2 % 2 == 0 and e2 >= 0, "bad exponent"
            e = int(e2) // 2
            if e < NQ:
                c = fmpz_poly([0] * e + [1])
                for p in range(r): c = tr(c * iqp(n[p], D.d[p]))
                res[deg] = res[deg] + c
            return
        v = 0
        while deg + v * D.w[i] <= Dmax:
            n[i] = v; rec(i + 1, deg + v * D.w[i]); v += 1
        n[i] = 0
    rec(0, 0)
    return res

def K_direct(k, Dmax):
    """K_k(x) summed from its definition over (r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k) (independent of the data/moves code)"""
    res = [fmpz_poly([0]) for _ in range(Dmax + 1)]
    def rec(i, rprev, sprev, deg, e, c):
        if i == k + 1:
            res[deg] = res[deg] + tr(c * iqp(rprev, 1)) if e < NQ else res[deg]
            return
        s = sprev
        while deg + (k - i + 1) * 2 * s <= Dmax:
            r = 0
            while (rprev is None or r <= rprev) and deg + r + 2 * s <= Dmax:
                ee = e + (r + s) ** 2 + s * s
                if ee < NQ:
                    cc = tr(c * iqp(s - sprev, 2))
                    if rprev is not None: cc = tr(cc * iqp(rprev - r, 1))
                    rec(i + 1, r, s, deg + r + 2 * s, ee, cc)
                r += 1
            s += 1
    # the factor q^e is applied at the end
    def rec2(i, rprev, sprev, deg, e, c):
        if i == k + 1:
            if e < NQ:
                res[deg] = res[deg] + tr(fmpz_poly([0] * e + [1]) * tr(c * iqp(rprev, 1)))
            return
        s = sprev
        while deg + (k - i + 1) * 2 * s <= Dmax:
            r = 0
            while (rprev is None or r <= rprev) and deg + r + 2 * s <= Dmax:
                ee = e + (r + s) ** 2 + s * s
                if ee < NQ:
                    cc = tr(c * iqp(s - sprev, 2))
                    if rprev is not None: cc = tr(cc * iqp(rprev - r, 1))
                    rec2(i + 1, r, s, deg + r + 2 * s, ee, cc)
                r += 1
            s += 1
    rec2(1, None, 0, 0, 0, fmpz_poly([1]))
    return res

def mulx(A, Bs, Dmax):
    res = [fmpz_poly([0]) for _ in range(Dmax + 1)]
    for i, a in enumerate(A):
        if a == 0: continue
        for j, b in enumerate(Bs):
            if i + j > Dmax: break
            res[i + j] = res[i + j] + tr(a * b)
    return res

# ---------------------------------------------------------------- the chains (found by breadth-first search)
CHAINS = {
    # k: (moves on the fully split K_k,  moves on the fully split B_{k+1},  index of J-factor fusions)
    2: dict(K=[('f', 0, 3), ('F', 0, 4, 5), ('F', 0, 1, 5)],
            B=[('F', 0, 3, 4), ('F', 0, 1, 2)],
            BJ=[('f', 0, 2)]),
    3: dict(K=[('f', 0, 4), ('F', 0, 6, 7), ('f', 1, 5), ('F', 0, 1, 2), ('f', 1, 3), ('F', 4, 8, 9), ('f', 3, 5),
               ('F', 0, 3, 7), ('F', 4, 6, 8)],
            B=[('F', 0, 4, 5), ('F', 0, 1, 2)],
            BJ=[('f', 2, 5), ('f', 1, 5)]),
}

def main():
    global NQ
    NQ = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    allok = True
    for k in (2, 3):
        print("=" * 100); print("k = %d :   K_%d(x)  =  B_%d(x) * J_%d(x)" % (k, k, k + 1, k))
        ch = CHAINS[k]
        # --- K side
        D = K_data(k); K0 = D.copy(); log = []
        for a in range(k): log.append("   SPLIT  %s" % D.labs[a]); D = split(D, a)
        print(" K side: species", K0.labs, "bases", K0.d)
        DK = run(D, ch['K'], log)
        print("\n".join(log)); print("   -> rank %d, bases %s" % (DK.r, DK.d))
        # --- B side:  K_k = B_{k+1} * C_k / C_1 ;  B_{k+1} splits off a free species = C_1 (k=2), resp. J_3 = AG_1(x^2;q^2) (k=3)
        Bd = B_data(k + 1); B0 = Bd.copy(); log = []
        for a in range(k): log.append("   SPLIT  %s" % Bd.labs[a]); Bd = split(Bd, a)
        Bd = run(Bd, ch['B'], log)
        print(" B side: species", B0.labs, "bases", B0.d); print("\n".join(log))
        if k == 2:
            # after the two DEFUSEs the species n1' is free with exponent n^2, weight 1, base q^2: it is the factor C_1(x)
            i0 = 0
            free = (all(Bd.M[i0][j] == 0 for j in range(Bd.r) if j != i0) and Bd.M[i0][i0] == 2 and Bd.B[i0] == 0
                    and Bd.w[i0] == 1 and Bd.d[i0] == 2)
            print("   species %s is free, = C_1(x):" % Bd.labs[i0], free); allok &= free
            keep = [i for i in range(Bd.r) if i != i0]
            Bt = Datum([[Bd.M[i][j] for j in keep] for i in keep], [Bd.B[i] for i in keep], [Bd.w[i] for i in keep],
                       [Bd.d[i] for i in keep], [Bd.labs[i] for i in keep])       # = B_3/C_1
            T = tensor(Bt, C_data(2))                                             # (B_3/C_1) * C_2 = B_3 * J_2
            Jname = "C_2 (so that the product is B_3 C_2 / C_1 = B_3 J_2)"
        else:
            T = tensor(Bd, AGQ_data(1))                                           # B_4 * AG_1(x^2;q^2) = B_4 * J_3
            Jname = "AG_1(x^2;q^2) = J_3"
        log = []; print("   tensor with", Jname)
        DB = run(T, ch['BJ'], log); print("\n".join(log)); print("   -> rank %d, bases %s" % (DB.r, DB.d))
        perm = same_up_to_permutation(DK, DB)
        print(" final data of the two chains coincide after the permutation", perm, ":", perm is not None)
        allok &= perm is not None
        if perm is not None:
            print("   matching of species:");
            for i, j in enumerate(perm): print("      %-45s <->  %s" % (DK.labs[i], DB.labs[j]))
        # --- independent numerical confirmation
        kd = K_direct(k, Dmax)
        e0 = evaluate(K0, Dmax); e1 = evaluate(DK, Dmax); e2 = evaluate(DB, Dmax)
        bj = mulx(evaluate(B0, Dmax), evaluate(C_data(k), Dmax), Dmax)
        kc = mulx(kd, evaluate(C_data(1), Dmax), Dmax)
        t1 = (kd == e0); t2 = (e0 == e1); t3 = (e1 == e2); t4 = (bj == kc)
        print(" series check, x^0..x^%d modulo q^%d:  K_k (definition) = F(K datum): %s;  = F(end of K chain): %s;"
              "  = F(end of B chain): %s;  B_{k+1} C_k = K_k C_1: %s" % (Dmax, NQ, t1, t2, t3, t4))
        allok &= t1 and t2 and t3 and t4
    print("=" * 100)
    print("ALL CHECKS PASSED" if allok else "SOME CHECK FAILED")

if __name__ == "__main__":
    main()
