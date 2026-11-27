"""Series check (NOT a proof) of the statement of proofs/X-A-parafermion.md:

  for an affine Lie algebra g of type X^{(1)} (X simply laced) or of twisted type A_{2l-1}^{(2)}, D_{l+1}^{(2)},
  E_6^{(2)}, D_4^{(3)}, with horizontal subalgebra of type X (= C_l, B_l, F_4, G_2 in the twisted cases), level k >= 2,
  and every class c of the root lattice Qbar of X modulo k*Qbar:

      ( sum_{n>=0} mult_{k Lambda_0}(k Lambda_0 + mu - n delta) q^n ) * prod_{j>=1} (1-q^j)^{mult(j delta)}
            =  q^{(mu|mu)/(2k)} * F_c(q),                                                   mu in c,

  where F_c is the part of the generalized Nahm sum of the Sun--Wang pair (X, A_{k-1}),

      F(q) = sum_p q^{(1/2) sum_{i,j,s,t} G_ij (min(s,t) - st/k) p_i^(s) p_j^(t)} / prod_{i,s} (q^{d_i};q^{d_i})_{p_i^(s)},

  with sum_s s p_i^(s) = (coordinate of mu at alpha_i) mod k for all i.  Here G = C(X)D(X) is the Gram matrix of
  the simple roots of X with short roots of square length 2, d_i = G_ii/2.

The weight multiplicities are computed from the generalized Cartan matrix alone, by the Racah (Weyl--Kac)
recursion  sum_{w in W} eps(w) mult(lambda + rho - w rho) = [lambda = Lambda]  (lambda dominant),
with W the affine Weyl group, generated here by breadth-first search.  Nothing is taken from the papers cited
in the proof file except the statement being tested.

Type A_{2l}^{(2)} (pair (T_l, A_{k-1})): node 0 is the short end node (Kac's labelling, a_0 = 2, a_0^vee = 1), the
variable is q = e^{-delta} (so the alpha_0-depth beta_0 contributes q^{beta_0/2}), the lattice of finite parts has
basis alpha_1, ..., alpha_{l-1}, alpha_l/2 with Gram matrix C(T_l) (the tadpole Cartan matrix), all d_i = 1, and
mult(j delta) = l for all j.  For this type the statement tested is NOT covered by a cited theorem (see the proof file).

usage: python3 xa_parafermion_certify.py            (default list of cases)
"""
import sys, itertools, math
from fractions import Fraction as Fr
from flint import fmpz_poly


# ----------------------------------------------------------------------------------------------- affine data
def affine(name):
    """returns (norms, edges, r, label) ; nodes 0..l ; norms = square lengths of simple roots (short = 2)"""
    t, l = name
    if t == 'A1':   # A_l^{(1)}
        if l == 1:
            return [2, 2], [(0, 1, 2)], 1          # double edge handled separately (a_01 = a_10 = -2)
        return [2] * (l + 1), [(i, i + 1, 1) for i in range(1, l)] + [(0, 1, 1), (0, l, 1)], 1
    if t == 'D1':   # D_l^{(1)}, l >= 4 : 0,1 - 2 - ... - (l-2) - (l-1), l
        e = [(0, 2, 1), (1, 2, 1)] + [(i, i + 1, 1) for i in range(2, l - 2)] + [(l - 2, l - 1, 1), (l - 2, l, 1)]
        return [2] * (l + 1), e, 1
    if t == 'E1':   # E_6^{(1)} : 1-2-3-4-5, 3-6, 6-0
        assert l == 6
        return [2] * 7, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (3, 6, 1), (6, 0, 1)], 1
    if t == 'A2odd':   # A_{2l-1}^{(2)}, l >= 3 : 0,1 - 2 - ... - (l-1) => l (alpha_l long) ; l = 2: 0 - 2 - 1 with alpha_2 long
        if l == 2:
            return [2, 2, 4], [(0, 2, 1), (1, 2, 1)], 2
        e = [(0, 2, 1), (1, 2, 1)] + [(i, i + 1, 1) for i in range(2, l)]
        return [2] * l + [4], e, 2
    if t == 'D2':   # D_{l+1}^{(2)}, l >= 2 : 0 <= 1 - ... - (l-1) => l ; alpha_0, alpha_l short
        return [2] + [4] * (l - 1) + [2], [(i, i + 1, 1) for i in range(l)], 2
    if t == 'E2':   # E_6^{(2)} : 0 - 1 - 2 <= 3 - 4 ; alpha_3, alpha_4 long
        return [2, 2, 2, 4, 4], [(i, i + 1, 1) for i in range(4)], 2
    if t == 'D3':   # D_4^{(3)} : 0 - 1 <= 2 ; alpha_2 long (norm 6)
        return [2, 2, 6], [(0, 1, 1), (1, 2, 1)], 3
    if t == 'A2even':   # A_{2l}^{(2)} : 0 => 1 - ... - (l-1) => l ; norms 1, 2, ..., 2, 4 (alpha_0 short end)
        return [1] + [2] * (l - 1) + [4], [(i, i + 1, 1) for i in range(l)], 2
    raise ValueError(name)


def gram(norms, edges):
    n = len(norms)
    G = [[0] * n for _ in range(n)]
    for i in range(n): G[i][i] = norms[i]
    for (i, j, m) in edges:
        v = -max(norms[i], norms[j]) // 2 * m
        G[i][j] = G[j][i] = v
    return G


def cartan(G):
    n = len(G)
    A = [[Fr(2 * G[i][j], G[i][i]) for j in range(n)] for i in range(n)]      # a_ij = <alpha_i^vee, alpha_j>
    assert all(x.denominator == 1 for row in A for x in row)
    return [[int(x) for x in row] for row in A]


def null_vector(A):
    """primitive positive integer vector a with A a = 0"""
    n = len(A)
    for a0 in (1, 2):
        for a in itertools.product(range(1, 7), repeat=n - 1):
            v = (a0,) + a
            if all(sum(A[i][j] * v[j] for j in range(n)) == 0 for i in range(n)): return list(v)
    raise ValueError("no null vector")


# ----------------------------------------------------------------------------------------------- Weyl group data
def weyl_shifts(A, labels0, N):
    """all (gamma, sign) with gamma = mu - w mu in root coordinates, mu = weight with Dynkin labels labels0 (regular
    dominant), w in W, gamma_0 <= N.  BFS along length-increasing reflections."""
    n = len(A)
    start = (tuple([0] * n), tuple(labels0))
    seen = {start[0]: 1}
    frontier = [start]
    sign = 1
    while frontier:
        sign = -sign
        new = []
        for gam, lab in frontier:
            for j in range(n):
                if lab[j] > 0:
                    g2 = list(gam); g2[j] += lab[j]
                    if g2[0] > N: continue
                    g2 = tuple(g2)
                    if g2 in seen: continue
                    l2 = tuple(lab[i] - lab[j] * A[i][j] for i in range(n))
                    seen[g2] = sign
                    new.append((g2, l2))
        frontier = new
    return seen


class Mult:
    def __init__(self, A, k, N):
        self.A = A; self.n = len(A); self.k = k; self.N = N
        self.shifts = [(g, s) for g, s in weyl_shifts(A, [1] * self.n, N).items() if any(g)]
        self.memo = {}

    def labels(self, beta):
        A = self.A
        return [(self.k if j == 0 else 0) - sum(beta[i] * A[j][i] for i in range(self.n)) for j in range(self.n)]

    def dom(self, beta):
        """dominant W-conjugate of Lambda - beta, as beta' (None if not <= Lambda)"""
        beta = list(beta)
        if min(beta) < 0: return None
        lab = self.labels(beta)
        A = self.A; n = self.n
        while True:
            for j in range(n):
                if lab[j] < 0:
                    c = lab[j]
                    beta[j] += c
                    if beta[j] < 0: return None
                    for i in range(n): lab[i] -= c * A[i][j]
                    break
            else:
                return tuple(beta)

    def mult(self, beta):
        b = self.dom(beta)
        if b is None: return 0
        return self._m(b)

    def _m(self, b):
        if b in self.memo: return self.memo[b]
        if not any(b):
            self.memo[b] = 1; return 1
        assert b[0] <= self.N
        tot = 0
        n = self.n
        for g, s in self.shifts:
            ok = True
            for i in range(n):
                if g[i] > b[i]: ok = False; break
            if not ok: continue
            b2 = self.dom([b[i] - g[i] for i in range(n)])
            if b2 is None: continue
            tot -= s * self._m(b2)
        self.memo[b] = tot
        return tot


# ----------------------------------------------------------------------------------------------- Nahm side
def trunc(p, N):
    return fmpz_poly(p.coeffs()[:N + 1])


def inv_poch(d, n, N):
    """1/(q^d;q^d)_n mod q^{N+1}"""
    c = [0] * (N + 1); c[0] = 1
    for j in range(1, n + 1):
        st = d * j
        for e in range(st, N + 1): c[e] += c[e - st]
    return fmpz_poly(c)


def nahm_sectors(Gh, k, N, dlist=None):
    """dict: class (tuple in (Z/k)^l) -> list of coefficients of q^{Delta_c} F_c(q) up to q^N, Delta_c in [0,1) excluded:
    returns for each class a dict {exact exponent (Fraction): coefficient} up to exponent N."""
    l = len(Gh)
    var = [(i, s) for i in range(l) for s in range(1, k)]
    R = len(var)
    M = [[Fr(Gh[i][j]) * (min(s, t) - Fr(s * t, k)) for (j, t) in var] for (i, s) in var]
    # Cholesky (float) for Fincke--Pohst
    L = [[0.0] * R for _ in range(R)]
    for i in range(R):
        for j in range(i + 1):
            s_ = float(M[i][j]) - sum(L[i][t] * L[j][t] for t in range(j))
            L[i][j] = math.sqrt(s_) if i == j else s_ / L[j][j]
    # p^T M p = |L^T p|^2 ; (L^T p)_t = sum_{i>=t} L[i][t] p_i ; enumerate from i = R-1 down to 0
    out = {}
    cache = {}
    bound = 2.0 * N + 1e-6
    p = [0] * R
    pts = []

    def rec(i, acc):
        if i < 0:
            pts.append(tuple(p)); return
        # (L^T p)_i = L[i][i] p_i + sum_{j>i} L[j][i] p_j
        c = sum(L[j][i] * p[j] for j in range(i + 1, R))
        v = 0
        while True:
            y = L[i][i] * v + c
            # y increases with v ; but may start negative: need all v >= 0 with acc + y^2 <= bound
            if acc + y * y <= bound:
                p[i] = v
                rec(i - 1, acc + y * y)
            elif y > 0:
                break
            v += 1
        p[i] = 0

    rec(R - 1, 0.0)
    for pt in pts:
        E = sum(M[a][b] * pt[a] * pt[b] for a in range(R) for b in range(R)) / 2
        if E > N: continue
        cls = tuple(sum(s * pt[a] for a, (i, s) in enumerate(var) if i == ii) % k for ii in range(l))
        poly = fmpz_poly([1])
        Ei = int(math.floor(E))
        room = N - Ei
        for a, (i, s) in enumerate(var):
            if pt[a]:
                key = ((Gh[i][i] // 2) if dlist is None else dlist[i], pt[a])
                if key not in cache: cache[key] = inv_poch(key[0], key[1], N)
                poly = trunc(poly * cache[key], room)
        d = out.setdefault(cls, {})
        frac = E - Ei
        for e, cf in enumerate(poly.coeffs()):
            if cf: d[Ei + e + frac] = d.get(Ei + e + frac, 0) + int(cf)
    return out, R


# ----------------------------------------------------------------------------------------------- the check
def central_charge(Gh, A_h_coxeter, k, dlist=None):
    trD = sum(Gh[i][i] // 2 for i in range(len(Gh))) if dlist is None else sum(dlist)
    return Fr(trD * (k - 1) * A_h_coxeter, A_h_coxeter + k)


COXETER = {'A1': lambda l: l + 1, 'D1': lambda l: 2 * l - 2, 'E1': lambda l: 12,
           'A2odd': lambda l: 2 * l, 'D2': lambda l: 2 * l, 'E2': lambda l: 12, 'D3': lambda l: 6, 'A2even': lambda l: 2 * l + 1}
XNAME = {'A1': 'A_%d', 'D1': 'D_%d', 'E1': 'E_%d', 'A2odd': 'C_%d', 'D2': 'B_%d', 'E2': 'F_4', 'D3': 'G_2', 'A2even': 'T_%d'}
GNAME = {'A1': 'A_%d^(1)', 'D1': 'D_%d^(1)', 'E1': 'E_%d^(1)', 'A2odd': 'A_%d^(2)', 'D2': 'D_%d^(2)', 'E2': 'E_6^(2)', 'D3': 'D_4^(3)', 'A2even': 'A_%d^(2)'}


def check(name, k, N):
    t, l = name
    norms, edges, r = affine(name)
    G = gram(norms, edges)
    if name == ('A1', 1): G = [[2, -2], [-2, 2]]
    A = cartan(G)
    a = null_vector(A)
    n = len(A)
    lh = n - 1
    a0 = a[0]                      # 1, except 2 for A_{2l}^{(2)}; the variable is q = e^{-delta}: beta_0 contributes q^{beta_0/a0}
    even = (t == 'A2even')
    if even:
        Gh = [[0] * lh for _ in range(lh)]          # Gram of alpha_1, ..., alpha_{l-1}, alpha_l/2 = C(T_l)
        for i in range(lh):
            Gh[i][i] = 2 if i < lh - 1 else 1
            if i < lh - 1: Gh[i][i + 1] = Gh[i + 1][i] = -1
        dlist = [1] * lh
        mim = lambda j: lh
    else:
        Gh = [row[1:] for row in G[1:]]
        dlist = None
        # multiplicities of imaginary roots j*delta: lh if r | j, else (rank X_N - lh)/(r-1)
        rankXN = {'A1': lh, 'D1': lh, 'E1': lh, 'A2odd': 2 * lh - 1, 'D2': lh + 1, 'E2': 6, 'D3': 4}[t]
        mim = lambda j: lh if j % r == 0 else (rankXN - lh) // (r - 1)
    N0 = a0 * N                    # depth bound in units of alpha_0
    heis = fmpz_poly([1])          # in the variable u = q^{1/a0}
    for j in range(1, N + 1):
        f = [0] * (a0 * j + 1); f[0] = 1; f[a0 * j] = -1
        for _ in range(mim(j)): heis = trunc(heis * fmpz_poly(f), N0)
    MU = Mult(A, k, N0)
    sect, R = nahm_sectors(Gh, k, N, dlist)
    ok = True
    nclass = 0
    total = {}

    def beta_of(m, b0):
        """root coordinates of k Lambda_0 - (weight with finite part sum m_i (basis vector i) at alpha_0-depth b0), or None"""
        if not even:
            return [b0] + [b0 * a[i + 1] - m[i] for i in range(lh)]
        # finite part = sum_{i<l} (b0 - beta_i) alpha_i + (b0 - 2 beta_l) alpha_l/2
        if (b0 - m[lh - 1]) % 2: return None
        return [b0] + [b0 - m[i] for i in range(lh - 1)] + [(b0 - m[lh - 1]) // 2]

    for cls in itertools.product(range(k), repeat=lh):
        nclass += 1
        d = sect.get(cls, {})
        for e, cf in d.items(): total[e] = total.get(e, 0) + cf
        # representatives mu = sum m_i (basis vector i), m_i = cls_i - k*eps_i, eps in {0,1}: test all of them
        for eps in itertools.product((0, 1), repeat=lh):
            m = [cls[i] - k * eps[i] for i in range(lh)]
            Delta = Fr(sum(Gh[i][j] * m[i] * m[j] for i in range(lh) for j in range(lh)), 2 * k)
            ser = []
            for b0 in range(N0 + 1):
                b = beta_of(m, b0)
                ser.append(0 if b is None else MU.mult(b))
            lhs = trunc(fmpz_poly(ser) * heis, N0)
            lhs = [int(c) for c in lhs.coeffs()] + [0] * (N0 + 1)
            rhs = [0] * (N0 + 1)             # q^Delta F_c in the variable u = q^{1/a0}; its exponents must be integers
            for e, cf in d.items():
                ee = a0 * (e + Delta)
                if ee.denominator != 1:
                    ok = False; print("   non-integral exponent in class", cls); break
                if ee <= N0: rhs[int(ee)] += cf
            if lhs[:N0 + 1] != rhs[:N0 + 1]:
                ok = False
                print("   MISMATCH class", cls, "rep", m, lhs[:N0 + 1], rhs[:N0 + 1])
    c = central_charge(Gh, COXETER[t](lh), k, dlist)
    es = sorted(total)
    lead = ", ".join("%s q^%s" % (total[e], e) for e in es[:6])
    xn = XNAME[t] % lh if '%' in XNAME[t] else XNAME[t]
    gn = GNAME[t] % ({'A1': lh, 'D1': lh, 'E1': lh, 'A2odd': 2 * lh - 1, 'D2': lh + 1, 'A2even': 2 * lh}.get(t, 0)) if '%' in GNAME[t] else GNAME[t]
    print("%-10s level %d  <->  pair (%s, A_%d): rank %d, %d sectors x %d representatives, |W-shifts| = %d, depth <= %d : %s"
          % (gn, k, xn, k - 1, R, nclass, 2 ** lh, len(MU.shifts), N, "OK" if ok else "FAILED"))
    print("      c(X,A_{k-1}) = %s ; F = %s + ..." % (c, lead))
    return ok


CASES = [(('A1', 1), 2, 14), (('A1', 1), 3, 14), (('A1', 1), 4, 12), (('A1', 1), 5, 10),
         (('A1', 2), 2, 10), (('A1', 2), 3, 8), (('A1', 2), 4, 6), (('A1', 3), 2, 8), (('A1', 3), 3, 6), (('A1', 4), 2, 6),
         (('D1', 4), 2, 6), (('D1', 4), 3, 4), (('D1', 5), 2, 5), (('D1', 6), 2, 4), (('E1', 6), 2, 3),
         (('A2odd', 2), 2, 12), (('A2odd', 2), 3, 10), (('A2odd', 2), 4, 8),
         (('A2odd', 3), 2, 8), (('A2odd', 3), 3, 6), (('A2odd', 4), 2, 6),
         (('D2', 3), 2, 8), (('D2', 3), 3, 6), (('D2', 4), 2, 6),
         (('D3', 2), 2, 12), (('D3', 2), 3, 10), (('D3', 2), 4, 8),
         (('E2', 4), 2, 5),
         (('A2even', 1), 2, 10), (('A2even', 1), 3, 10), (('A2even', 1), 4, 9), (('A2even', 1), 5, 8), (('A2even', 1), 6, 8),
         (('A2even', 1), 7, 7), (('A2even', 1), 8, 6),
         (('A2even', 2), 2, 6), (('A2even', 2), 3, 5), (('A2even', 2), 4, 4), (('A2even', 3), 2, 5), (('A2even', 3), 3, 3), (('A2even', 4), 2, 4)]

if __name__ == "__main__":
    allok = True
    for name, k, N in CASES:
        allok &= check(name, k, N)
        sys.stdout.flush()
    print("ALL CHECKS PASSED" if allok else "SOME CHECK FAILED")
