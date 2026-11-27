"""sw_audit.py -- numerical audit of the explicit formulas of K. Sun, H. Wang, arXiv:2604.00847
(eqs (1.8)-(1.12), (2.1)-(2.8), (3.2), (3.5)-(3.41), Table 1).  Exact integer coefficients.

usage: python3 sw_audit.py SECTION [N]        (SECTION = name of a function below, or 'all')
Conventions (see SUN-WANG-AUDIT.md): characters include q^{h-c/24}; "effective" models = same characters,
labelled by h - h_min; R characters of SM(p,p') are normalised with leading coefficient 1
(q^{h-c/24} prod (1+q^n)/(1-q^n) * sum_k(...)).
"""
import sys, time
from fractions import Fraction as Fr
from sw_lib import Ctx, Ser, compare, report
from fusion import inv, cartanA, cartanT, cartanD, cartanE
from fusion2 import pair, cartanB, cartanC

F = Fr
LOG = []


def scal(M, a):
    return [[F(x) * a for x in row] for row in M]


def minmat(r):
    return [[min(i, j) + 1 for j in range(r)] for i in range(r)]


def rep(name, l, r):
    return report(name, l, r, LOG)


# ------------------------------------------------------------------ (1.8)-(1.10): (T2,A1), Zagier Ex. 7
def s_T2A1(N=60):
    c = Ctx(168, N)
    A = [[F(1, 2), F(-1, 2)], [F(-1, 2), 1]]
    NS = lambda w: c.sm_eff(2, 84, 'NS', w)
    R = lambda w: c.sm_eff(2, 84, 'R', w)
    W = c.sm_weights(2, 84)
    print("SM(84,2): #NS = %d, #R = %d" % (sum(1 for k in W if k[0] == 'NS'), sum(1 for k in W if k[0] == 'R')))
    f1 = c.nahm(A, [0, 0], C=F(-5, 84))
    f2 = c.nahm(A, [F(1, 2), F(-1, 2)], C=F(1, 21))
    f3 = c.nahm(A, [F(1, 2), 0], C=F(1, 84))
    rep("(1.8)", f1, NS(0) + NS(F(1, 2)) + NS(F(5, 2)) + NS(5) + 2 * R(F(1, 4)) + 2 * R(F(5, 4)))
    rep("(1.9)", f2, NS(F(5, 14)) + NS(F(6, 7)) + NS(F(13, 7)) - NS(F(20, 7)) + 2 * R(F(3, 28)) + 2 * R(F(87, 28)))
    rep("(1.10)", f3, NS(F(1, 14)) - NS(F(15, 14)) + NS(F(11, 7)) + NS(F(57, 14)) + 2 * R(F(23, 28)) - 2 * R(F(135, 28)))


# ------------------------------------------------------------------ (1.11), (1.12): (T1,D_r)
def s_T1D(N=40, rmax=8):
    for r in range(3, rmax + 1):
        den = 8 * (2 * r + 1)
        NN = N if r <= 6 else max(20, N // 2)
        c = Ctx(den, NN)
        A = inv(cartanD(r))
        t0 = time.time()
        f = c.nahm(A, [0] * r, C=F(-r, 8 * (2 * r + 1)))
        p, pp = 2, 8 * r + 4
        W = c.sm_weights(p, pp)
        ns = sorted(k[1] for k in W if k[0] == 'NS')
        rr = sorted(k[1] for k in W if k[0] == 'R')
        assert len(ns) == 2 * r + 1 and len(rr) == 2 * r + 1
        rhs = c.sm_eff(p, pp, 'NS', ns[0]) + c.sm_eff(p, pp, 'NS', ns[2 * r]) + 2 * c.sm_eff(p, pp, 'R', rr[r])
        rep("(1.11) r=%d [NS %s,%s; R %s]" % (r, ns[0], ns[2 * r], rr[r]), f, rhs)
        # (1.12): products
        m = 8 * r + 4
        p1 = c.prod(lambda n: -1 if (n % 4 != 2 and n % m not in (0, 4 * r + 1, m - 4 * r - 1)) else 0, F(1, 2))
        p2 = c.prod(lambda n: -1 if (n % 4 != 2 and n % m not in (0, 1, m - 1)) else 0, F(1, 2)).shift(F(r, 2))
        m2 = 4 * r + 2
        p3 = (c.prod(lambda n: -1 if n % 2 == 1 else 0) * c.prod(lambda n: -1 if n % m2 not in (0, r + 1, m2 - r - 1) else 0)).shift(F(r, 8))
        f0 = f.shift(F(r, 8 * (2 * r + 1)))
        rep("(1.12) r=%d" % r, f0, p1 + p2 + 2 * p3)
        print("     time %.1fs" % (time.time() - t0), flush=True)


# ------------------------------------------------------------------ section 2 identities quoted as known
def s_sec2(N=60):
    # (2.1) Andrews-Gordon, all s
    for r in range(1, 5):
        c = Ctx(1, N)
        for s in range(1, r + 2):
            B = [max(0, j - s + 2) for j in range(r)]      # N_s+...+N_r = sum_j max(0, j-s+1) n_j  (1-based j)
            B = [max(0, (j + 1) - s + 1) for j in range(r)]
            f = c.nahm(scal(minmat(r), 2), B)
            m = 2 * r + 3
            rhs = c.prod(lambda n: -1 if n % m not in (0, s, m - s) else 0)
            rep("(2.1) AG r=%d s=%d" % (r, s), f, rhs)
    # (2.2) (T1,T_2r)
    for r in range(1, 4):
        c = Ctx(2, N)
        f = c.nahm(minmat(2 * r), [0] * (2 * r))
        m = 4 * (r + 1)
        rhs = c.prod(lambda n: -1 if (n % 4 != 2 and n % m not in (0, 2 * r + 1, m - 2 * r - 1)) else 0, F(1, 2))
        rep("(2.2) (T1,T_%d)" % (2 * r), f, rhs)
    # (2.3) (T1,T_{2r-1})
    for r in range(1, 5):
        c = Ctx(2, N)
        f = c.nahm(minmat(2 * r - 1), [0] * (2 * r - 1))
        m = 2 * r + 1
        rhs = c.prod(lambda n: 1 if n % 4 == 2 else (-1 if n % 2 == 1 else 0), F(1, 2)) * c.prod(lambda n: -1 if n % m not in (0, r, m - r) else 0)
        rep("(2.3) (T1,T_%d)" % (2 * r - 1), f, rhs)
    # (2.4) (A1,B_r): exponent (1/2) sum M_i^2, M_i = 2n_i+...+2n_{r-1}+n_r, M_r = n_r ; d=(2,..,2,1)
    for r in range(2, 6):
        c = Ctx(2, N)
        M, d = pair(('A', 1), ('B', r))
        # check that the pair matrix reproduces the printed quadratic form
        Q = [[F(0)] * r for _ in range(r)]
        for i in range(r):
            v = [0] * r
            if i < r - 1:
                for j in range(i, r - 1):
                    v[j] = 2
                v[r - 1] = 1
            else:
                v[r - 1] = 1
            for a in range(r):
                for b in range(r):
                    Q[a][b] += v[a] * v[b]
        assert Q == [[F(x) for x in row] for row in M], "(2.4): printed quadratic form != A(A1,B_r)D"
        f = c.nahm(M, [0] * r, d)
        # (-q^{r/2}, -q^{r/2+1}, q^{r+1}; q^{r+1})_inf / (q^2;q^2)_inf  = sum_k q^{(r+1)k^2/2 + k/2}/(q^2;q^2)_inf
        th = c.theta(lambda k: F((r + 1) * k * k + k, 2))
        rhs = th * c.prod(lambda n: -1 if n % 2 == 0 else 0)
        rep("(2.4) (A1,B_%d)" % r, f, rhs)
    # (2.5) Bressoud
    for r in range(2, 5):
        c = Ctx(1, N)
        for s in range(1, r + 2):
            B = [max(0, (j + 1) - s + 1) for j in range(r)]
            f = c.nahm(scal(minmat(r), 2), B, [1] * (r - 1) + [2])
            m = 2 * r + 2
            if s == r + 1:
                th = c.theta(lambda k: (r + 1) * k * k, lambda k: (-1) ** abs(k))
            else:
                th = c.theta(lambda k: (r + 1) * k * k + (r + 1 - s) * k, lambda k: (-1) ** abs(k))
            rep("(2.5) Bressoud r=%d s=%d" % (r, s), f, th * c.euler_inv())
    # (2.6),(2.7): (A1,D_r) restricted sums: tested in total (sum of both)
    for r in range(3, 7):
        c = Ctx(4, N)
        f = c.nahm(scal(inv(cartanD(r)), 2), [0] * r)
        rhs = c.theta(lambda k: F(r * k * k, 4)) * c.euler_inv()
        rep("(2.6)+(2.7) (A1,D_%d)" % r, f, rhs)
    # (2.8) Warnaar (T1,C_r)
    for r in range(2, 8):
        c = Ctx(2, N)
        f = c.nahm(minmat(r), [0] * r, [1] * (r - 1) + [2])
        m = 2 * r + 3      # in t = q^{1/2}: (t^{r+1}, t^{r+2}, t^{2r+3}; t^{2r+3}) / ((-t^2;t^2)(t;t))
        num = c.prod(lambda n: 1 if n % m in (0, r + 1, r + 2) else 0, F(1, 2))
        den_ = c.prod(lambda n: (-1 if n % 2 == 0 else 0) + (1 if n % 4 == 0 else 0) - 1 + 0, F(1, 2))
        # 1/((-q;q)_inf (q^{1/2};q^{1/2})_inf): (-q;q) = (q^2;q^2)/(q;q) -> 1/(-q;q) = (q;q)/(q^2;q^2) = prod_{t-even}(1-t^n) / prod_{4|n}(1-t^n)
        den_ = c.prod(lambda n: (1 if n % 2 == 0 else 0) - (1 if n % 4 == 0 else 0) - 1, F(1, 2))
        rep("(2.8) (T1,C_%d)" % r, f, num * den_)
    # G2 identities at the end of section 2
    c = Ctx(1, N)
    f = c.nahm([[2, 3], [3, 6]], [0, 0], [1, 3])
    rep("sec.2 (T1,G2) Kanade-Russell", f, c.prod(lambda n: -1 if n % 9 in (1, 3, 6, 8) else 0))
    f = c.nahm([[2, -3], [-3, 6]], [0, 0], [1, 3])
    a = c.prod(lambda n: -1 if n % 9 in (1, 3, 6, 8) else 0)
    b1 = c.prod(lambda n: -1 if n % 9 in (3, 6) else 0)
    b2 = c.prod(lambda n: -1 if n % 9 in (2, 4, 5, 7) else 0)
    rep("sec.2 (G2,T1) Wang-Wang Conj 3.6", f, a * a + (b1 * b1 * b2).shift(1))
    # consistency of the printed quadratic forms with A(X,Y)D(X,Y)
    M, d = pair(('T', 1), ('G', 2)); print("   A(T1,G2)D =", [[str(x) for x in r_] for r_ in M], d)
    M, d = pair(('G', 2), ('T', 1)); print("   A(G2,T1)D =", [[str(x) for x in r_] for r_ in M], d)


# ------------------------------------------------------------------ (3.2), Table 1 rows (A1,E_n), (A1,A_r)
def s_A1E(N=60):
    c = Ctx(48, N)
    f = c.nahm(scal(inv(cartanE(8)), 2), [0] * 8, C=F(-1, 48))
    rep("(3.2)/Table (A1,E8) = chi_0 M(4,3)", f, c.vir(3, 4, 1, 1))
    # right side of (3.2): sum_m q^{2m^2}/(q)_{2m}
    from flint import fmpz_poly
    d = {}
    acc = None
    cc = Ctx(1, N)
    s = Ser(1, 0, fmpz_poly([0]), N)
    m = 0
    while 2 * m * m < N:
        s = s + cc.prod(lambda n: -1 if n <= 2 * m else 0).shift(2 * m * m)
        m += 1
    f1 = cc.nahm(scal(inv(cartanE(8)), 2), [0] * 8)
    rep("(3.2) as printed", f1, s)
    c = Ctx(240, N)
    f = c.nahm(scal(inv(cartanE(7)), 2), [0] * 7, C=F(-7, 240))
    rep("Table (A1,E7): chi_0+chi_{3/2} M(5,4)", f, c.vir(4, 5, 1, 1) + c.vir(4, 5, 3, 1))
    c = Ctx(84, N)
    f = c.nahm(scal(inv(cartanE(6)), 2), [0] * 6, C=F(-1, 28))
    # M(7,6): h_{r,s}, p=6, pp=7
    print("   M(7,6) weights:", {(r, s): str(c.vir_h(6, 7, r, s)) for r in range(1, 6) for s in range(1, 7) if (r, s) <= (6 - r, 7 - s)})
    rep("Table (A1,E6): chi_0+chi_5+2chi_{4/3} M(7,6)", f, c.vir(6, 7, 1, 1) + c.vir(6, 7, 5, 1) + 2 * c.vir(6, 7, 3, 1))


def s_A1A(N=60):
    for r in range(1, 7):
        k = r + 1
        cpf = F(2 * (k - 1), k + 2)
        den = 24 * (k + 2) * k
        c = Ctx(den, N if r <= 4 else N * 2 // 3)
        f = c.nahm(scal(inv(cartanA(r)), 2), [0] * r, C=-cpf / 24)
        rhs = 0
        for m in range(k):
            rhs = rhs + c.pf(k, 0, 2 * m)
        rep("Table (A1,A_%d) = sum_m PF Z_%d chi^0_{2m}" % (r, k), f, rhs)


# ------------------------------------------------------------------ (3.5) (A1,C_r)
def s_A1C(N=80):
    for r in range(2, 8):
        c = Ctx(24, N)
        f = c.nahm(scal(minmat(r), 2), [0] * r, [1] * (r - 1) + [2], C=F(-1, 24))
        k = 4 * (r + 1)
        rep("(3.5) (A1,C_%d): U(1)_%d chi_0-chi_{h=%d}" % (r, k, r + 1), f, c.u1(k, 0) - c.u1(k, k))
        M, d = pair(('A', 1), ('C', r))
        assert M == [[F(x) for x in row] for row in scal(minmat(r), 2)] and d == [1] * (r - 1) + [2]


# ------------------------------------------------------------------ (3.6)-(3.13) (T1,C_r)
def s_T1C(N=100, rmax=9):
    for r in range(2, rmax + 1):
        cc = F(3 * (r + 1), 2 * r + 3)
        den = 16 * (2 * r + 3)
        NN = N if r <= 6 else (N * 3 // 5 if r <= 8 else N * 2 // 5)
        c = Ctx(den, NN)
        M, d = pair(('T', 1), ('C', r))
        assert M == [[F(x) for x in row] for row in minmat(r)] and d == [1] * (r - 1) + [2]
        t0 = time.time()
        f = c.nahm(M, [0] * r, d, C=-cc / 24)
        p, pp = 4, 4 * r + 6
        W = c.sm_weights(p, pp)
        ns = sorted(k[1] for k in W if k[0] == 'NS')
        w = r // 2 + 1
        lab0 = W[('NS', F(0))]; lab1 = W[('NS', F(w))]
        rep("(3.6) r=%d: NS_0 %s - NS_%d %s" % (r, lab0, w, lab1), f, c.sm_eff(p, pp, 'NS', 0) - c.sm_eff(p, pp, 'NS', w))
        print("     SM(%d,4): %d NS weights; time %.1fs" % (pp, len(ns), time.time() - t0), flush=True)
    # r = 2 : (3.7)-(3.9)
    c = Ctx(112, N)
    M = minmat(2); d = [1, 2]
    W = c.sm_weights(4, 14)
    print("   SM_eff(14,4) NS weights:", [str(x) for x in sorted(k[1] for k in W if k[0] == 'NS')])
    NS = lambda w: c.sm_eff(4, 14, 'NS', w)
    rep("(3.7)", c.nahm(M, [0, 0], d, C=F(-3, 56)), NS(0) - NS(2))
    rep("(3.8)", c.nahm(M, [0, 1], d, C=F(1, 56)), NS(F(1, 14)) - NS(F(15, 14)))
    rep("(3.9)", c.nahm(M, [1, 1], d, C=F(9, 56)), NS(F(3, 14)) - NS(F(45, 14)))
    # r = 3 : (3.10)-(3.13)
    c = Ctx(144, N)
    M = minmat(3); d = [1, 1, 2]
    W = c.sm_weights(4, 18)
    print("   SM_eff(18,4) NS weights:", [str(x) for x in sorted(k[1] for k in W if k[0] == 'NS')])
    NS = lambda w: c.sm_eff(4, 18, 'NS', w)
    rep("(3.10)", c.nahm(M, [0, 0, 0], d, C=F(-1, 18)), NS(0) - NS(2))
    rep("(3.11)", c.nahm(M, [0, 0, 1], d, C=0), NS(F(1, 18)) - NS(F(55, 18)))
    rep("(3.12)", c.nahm(M, [0, 1, 1], d, C=F(1, 9)), NS(F(1, 6)) - NS(F(7, 6)))
    rep("(3.13)", c.nahm(M, [1, 1, 2], d, C=F(5, 18)), NS(F(1, 3)) - NS(F(13, 3)))


# ------------------------------------------------------------------ (3.14)-(3.18): F4, G2
def s_FG(N=100):
    c = Ctx(168, N)
    M, d = pair(('A', 1), ('F', 4))
    printed = [[8, 12, 8, 4], [12, 24, 16, 8], [8, 16, 12, 6], [4, 8, 6, 4]]
    print("   A(A1,F4)D == printed form of (3.14):", M == [[F(x) for x in row] for row in printed], d)
    f = c.nahm(printed, [0] * 4, [2, 2, 1, 1])
    rep("(3.14)", f, c.prod(lambda n: -1 if n % 14 in (2, 3, 4, 10, 11, 12) else 0))
    print("   M(7,6) eff weights 0 and 5:", c.vir_h(6, 7, 1, 1), c.vir_h(6, 7, 5, 1))
    rep("(3.15)", f, (c.vir(6, 7, 1, 1) - c.vir(6, 7, 5, 1)).shift(F(1, 28)))
    c = Ctx(120, N)
    M, d = pair(('T', 1), ('F', 4))
    printed = [[4, 6, 4, 2], [6, 12, 8, 4], [4, 8, 6, 3], [2, 4, 3, 2]]
    print("   A(T1,F4)D == printed form of (3.16):", M == [[F(x) for x in row] for row in printed], d)
    f = c.nahm(printed, [0] * 4, [2, 2, 1, 1])
    rep("(3.16)", f, c.prod(lambda n: -1 if n % 10 in (1, 2, 4, 6, 8, 9) else 0))
    rep("(3.17)", f, (c.vir_eff(2, 5, 0) * (c.vir(5, 6, 1, 1) - c.vir(5, 6, 4, 1))).shift(F(1, 20)))
    c = Ctx(144, N)
    M, d = pair(('A', 1), ('G', 2))
    print("   A(A1,G2)D =", [[str(x) for x in r_] for r_ in M], d)
    f = c.nahm(M, [0, 0], d, C=F(-1, 24))
    rep("(3.18)", f, c.u1(36, 0) - c.u1(36, 12) - c.u1(36, 24) + c.u1(36, 36))


# ------------------------------------------------------------------ (3.19)-(3.21): E6, E8
def s_E(N=60):
    c = Ctx(60, N)
    f = c.nahm(inv(cartanE(6)), [0] * 6, C=F(-1, 20))
    rep("(3.19)", f, c.vir_eff(2, 5, 0) * (c.vir(5, 6, 1, 1) + c.vir(5, 6, 4, 1) + 2 * c.vir(5, 6, 1, 3)))
    c = Ctx(33, N)
    f = c.nahm(inv(cartanE(8)), [0] * 8, C=F(-1, 33))
    rep("(3.20)", f, c.vir_eff(2, 11, 0))
    c = Ctx(4, N)
    f = c.nahm(inv(cartanE(7)), [0] * 7)
    th = c.theta(lambda k: F(3 * k * k, 4)) + c.theta(lambda k: F(7 * k * k, 4))
    two_f = f * 2
    rep("(T1,E7) [proofs/E7.md, not in S-W]", two_f, th * c.euler_inv())


def s_E8T1(N=24):
    c = Ctx(1, N)
    t0 = time.time()
    f = c.nahm(cartanE(8), [0] * 8, count=True)
    cs = f.coeffs()
    print("   (E8,T1) Nahm sum: ", [cs.get(i, 0) for i in range(min(N, 12))], " time %.0fs" % (time.time() - t0))
    ok = [cs.get(i, 0) for i in range(5)] == [1, 120, 1660, 12320, 68210]
    print("(3.21) printed coefficients 1+120q+1660q^2+12320q^3+68210q^4 :", "AGREE" if ok else "DISAGREE")


# ------------------------------------------------------------------ (3.22)-(3.25)
def s_rank12(N=100):
    c = Ctx(120, N)
    rep("(3.22)", c.nahm([[F(1, 2)]], [0], C=F(-1, 40)), c.vir_eff(3, 5, 0) + c.vir_eff(3, 5, F(1, 4)))
    rep("(3.23)", c.nahm([[F(1, 2)]], [F(1, 2)], C=F(1, 40)), c.vir_eff(3, 5, F(1, 20)) + c.vir_eff(3, 5, F(4, 5)))
    c = Ctx(24, N)
    rep("(3.24)", c.nahm(inv(cartanA(2)), [0, 0], C=F(-1, 24)), c.u1(3, 0) + 2 * c.u1(3, 2))
    rep("(3.25)", c.nahm(cartanA(2), [0, 0], C=F(-1, 24)), c.u1(4, 0) + c.u1(4, 4))
    rep("(3.25) (A1)_1 = eta^-1 sum q^{m^2}", c.nahm(cartanA(2), [0, 0], C=F(-1, 24)), (c.theta(lambda k: k * k) * c.euler_inv()).shift(F(-1, 24)))


# ------------------------------------------------------------------ (3.26)-(3.29): (T2,T1), Zagier Ex. 3
def s_T2T1(N=60):
    c = Ctx(96, N)
    A = [[1, -1], [-1, 2]]
    assert pair(('T', 2), ('T', 1))[0] in ([[2, -1], [-1, 1]], [[F(2), F(-1)], [F(-1), F(1)]])
    f = {1: c.nahm(A, [F(-3, 2), 2], C=F(25, 24)), 2: c.nahm(A, [0, 0], C=F(-5, 96)), 3: c.nahm(A, [F(-1, 2), 1], C=F(1, 6)),
         4: c.nahm(A, [F(1, 2), 0], C=F(1, 24)), 5: c.nahm(A, [0, 1], C=F(19, 96))}
    W = c.sm_weights(2, 8)
    print("   SM_eff(8,2) weights:", sorted((k[0], str(k[1])) for k in W))
    FNS = c.prod(lambda n: (-1 if n % 2 == 1 else 0) + (1 if n % 4 == 2 else 0), F(1, 2)).shift(F(-1, 48))   # q^{-1/48} prod (1+q^{n-1/2})
    FR = c.prod(lambda n: -1 if n % 2 == 1 else 0).shift(F(1, 24))                                        # q^{1/24} prod (1+q^n)
    NS = lambda w: c.sm_eff(2, 8, 'NS', w)
    R = lambda w: c.sm_eff(2, 8, 'R', w)
    rep("(3.26) f2", f[2], NS(0) * FNS)
    rep("(3.27) AS PRINTED: f3 = chi_{R,1/32} chi^F_{R,1/16}", f[3], R(F(1, 32)) * FR)
    rep("(3.28) AS PRINTED: f4 = chi_{R,5/32} chi^F_{R,1/16}", f[4], R(F(5, 32)) * FR)
    rep("(3.27) corrected: f4 = chi_{R,1/32} chi^F_{R,1/16}", f[4], R(F(1, 32)) * FR)
    rep("(3.28) corrected: f3 = 2 chi_{R,5/32} chi^F_{R,1/16}", f[3], 2 * R(F(5, 32)) * FR)
    rep("(3.29) f5", f[5], NS(F(1, 4)) * FNS)
    rep("f1 = 2 f4", f[1], 2 * f[4])
    rep("f5 with C = 19/24 [Zagier's misprint, corrected by Wang and by Sun-Wang]", c.nahm(A, [0, 1], C=F(19, 24)), NS(F(1, 4)) * FNS)


# ------------------------------------------------------------------ (3.30): (T1,A3), Zagier Ex. 9
def s_T1A3(N=60):
    c = Ctx(56, N)
    A = [[1, F(-1, 2)], [F(-1, 2), F(3, 4)]]
    f1 = c.nahm(A, [F(-1, 2), F(1, 4)], C=F(1, 28))
    f2 = c.nahm(A, [0, 0], C=F(-3, 56))
    f3 = c.nahm(A, [0, F(1, 2)], C=F(1, 56))
    W = c.sm_weights(2, 28)
    print("   SM_eff(28,2) NS:", [str(k[1]) for k in sorted(W) if k[0] == 'NS'], " R:", [str(k[1]) for k in sorted(W) if k[0] == 'R'])
    NS = lambda w: c.sm_eff(2, 28, 'NS', w)
    R = lambda w: c.sm_eff(2, 28, 'R', w)
    rep("p.13 f1", f1, NS(F(3, 14)) + NS(F(5, 7)) + 2 * R(F(5, 56)))
    rep("p.13 f2", f2, NS(0) + NS(F(3, 2)) + 2 * R(F(3, 8)))
    rep("p.13 f3", f3, NS(F(1, 14)) - NS(F(15, 14)) + 2 * R(F(53, 56)))
    g = c.nahm(inv(cartanA(3)), [0, 0, 0], C=F(-3, 56))
    rep("(3.30) (T1,A3) = RHS", g, NS(0) + NS(F(3, 2)) + 2 * R(F(3, 8)))
    rep("'(T1,A3) is identical to f1' [as printed]", g, f1)
    rep("(T1,A3) = f2", g, f2)


# ------------------------------------------------------------------ (3.31)-(3.36)
def s_A1A2(N=80):
    c = Ctx(30, N)
    A = [[F(4, 3), F(2, 3)], [F(2, 3), F(4, 3)]]
    f1 = c.nahm(A, [F(-2, 3), F(-1, 3)], C=F(1, 30))
    f2 = c.nahm(A, [F(-1, 3), F(-2, 3)], C=F(1, 30))
    f3 = c.nahm(A, [0, 0], C=F(-1, 30))
    rep("Ex.10 f1 = f2", f1, f2)
    rep("(3.31)", f1, 2 * c.vir(5, 6, 3, 3) + c.vir(5, 6, 2, 1) + c.vir(5, 6, 3, 1))
    rep("(3.32)", f3, c.vir(5, 6, 1, 1) + c.vir(5, 6, 4, 1) + 2 * c.vir(5, 6, 1, 3))
    # (3.33),(3.34): (A2,A1), D_2A = extension of Ising x M(5,4) by (1/2,3/2)
    c = Ctx(240, N)
    A = [[1, F(-1, 2)], [F(-1, 2), 1]]
    f1 = c.nahm(A, [F(-1, 2), 0], C=F(1, 20))
    f2 = c.nahm(A, [0, F(-1, 2)], C=F(1, 20))
    f3 = c.nahm(A, [0, 0], C=F(-1, 20))
    I = lambda w: c.vir_eff(3, 4, w)
    T = lambda w: c.vir_eff(4, 5, w)
    chi = {0: I(0) * T(0) + I(F(1, 2)) * T(F(3, 2)), F(1, 2): I(F(1, 16)) * T(F(7, 16)),
           F(1, 10): I(F(1, 16)) * T(F(3, 80)), F(3, 5): I(0) * T(F(3, 5)) + I(F(1, 2)) * T(F(1, 10)),
           'half': I(0) * T(F(3, 2)) + I(F(1, 2)) * T(0), 'tenth': I(0) * T(F(1, 10)) + I(F(1, 2)) * T(F(3, 5))}
    rep("Ex.11 f1 = f2", f1, f2)
    rep("(3.33)", f1, 3 * chi[F(1, 10)] + chi[F(3, 5)])
    rep("(3.34)", f3, chi[0] + 3 * chi[F(1, 2)])
    rep("  D2A: [0,3/2]+[1/2,0] = [1/16,7/16]", chi['half'], chi[F(1, 2)])
    rep("  D2A: [0,1/10]+[1/2,3/5] = [1/16,3/80]", chi['tenth'], chi[F(1, 10)])
    # (3.35),(3.36): (A1,A3)
    c = Ctx(48, N)
    f = c.nahm(scal(inv(cartanA(3)), 2), [0, 0, 0], C=F(-1, 24))
    rep("(3.35)", f, c.u1(3, 0) + c.u1(3, 3))
    rep("(3.36)", f, c.pf(4, 0, 0) + c.pf(4, 0, 4) + 2 * c.pf(4, 0, 2))


# ------------------------------------------------------------------ (3.37)-(3.41): (A1,B3), (A1,A5)
def s_A1B3(N=80):
    c = Ctx(96, N)
    M, d = pair(('A', 1), ('B', 3))
    printedA = [[2, 2, 2], [2, 4, 4], [1, 2, 3]]
    AD = [[F(printedA[i][j] * [2, 2, 1][j]) for j in range(3)] for i in range(3)]
    print("   A(A1,B3)D == printed A * diag(2,2,1):", M == AD, d)
    W = c.sm_weights(6, 8)
    print("   SM(8,6) NS:", [str(k[1]) for k in sorted(W) if k[0] == 'NS'], " R:", [str(k[1]) for k in sorted(W) if k[0] == 'R'])
    NS = lambda w: c.sm_eff(6, 8, 'NS', w)
    R = lambda w: c.sm_eff(6, 8, 'R', w)
    rep("(3.37)", c.nahm(AD, [0, 0, 0], d, C=F(-5, 96)), NS(0) - NS(3))
    rep("(3.38)", c.nahm(AD, [0, -1, -1], d, C=F(-1, 48)), NS(F(1, 32)) - NS(F(33, 32)))
    rep("(3.39)", c.nahm(AD, [-1, 0, F(1, 2)], d, C=F(1, 24)), R(F(3, 32)) - R(F(67, 32)))
    c = Ctx(96, min(N, 60))
    f = c.nahm(scal(inv(cartanA(5)), 2), [0] * 5, C=F(-5, 96))
    rep("(3.40)", f, c.pf(6, 0, 0) + c.pf(6, 0, 6) + 2 * c.pf(6, 0, 2) + 2 * c.pf(6, 0, 4))
    NS = lambda w: c.sm_eff(6, 8, 'NS', w)
    rep("(3.41)", f, NS(0) + NS(3) + 2 * NS(F(5, 6)))


# ------------------------------------------------------------------ Table 1: central charges
def s_table(N=0):
    hcox = {'A': lambda r: r + 1, 'B': lambda r: 2 * r, 'C': lambda r: 2 * r, 'D': lambda r: 2 * r - 2, 'T': lambda r: 2 * r + 1,
            'E': lambda r: {6: 12, 7: 18, 8: 30}[r], 'F': lambda r: 12, 'G': lambda r: 6}
    trD = {'A': lambda r: r, 'B': lambda r: 2 * r - 1, 'C': lambda r: r + 1, 'D': lambda r: r, 'T': lambda r: r,
           'E': lambda r: r, 'F': lambda r: 6, 'G': lambda r: 4}

    def cc(X, Y):
        return F(trD[X[0]](X[1]) * trD[Y[0]](Y[1]) * hcox[X[0]](X[1]), hcox[X[0]](X[1]) + hcox[Y[0]](Y[1]))
    rows = []
    for r in range(1, 8):
        rows += [(('A', 1), ('A', r), F(2 * r, r + 3)), (('A', 1), ('T', r), F(2 * r, 2 * r + 3)), (('T', 1), ('A', r), F(3 * r, r + 4)),
                 (('T', 1), ('T', 2 * r), F(3 * r, 2 * r + 2)), (('T', 1), ('T', 2 * r + 1), F(3 * (2 * r + 1), 4 * r + 6))]
        if r >= 2:
            rows += [(('A', 1), ('C', r), F(1)), (('T', 1), ('C', r), F(3 * (r + 1), 2 * r + 3)),
                     (('A', 1), ('B', r), F(2 * r - 1, r + 1)), (('T', 1), ('B', r), F(3 * (2 * r - 1), 2 * r + 3))]
        if r >= 3:
            rows += [(('A', 1), ('D', r), F(1)), (('T', 1), ('D', r), F(3 * r, 2 * r + 1))]
        for k in range(1, 4):
            rows += [(('A', 2 * r - 1), ('T', k), F(2 * k * r * (2 * r - 1), 2 * (k + r) + 1))]
    rows += [(('A', 1), ('E', 6), F(6, 7)), (('A', 1), ('E', 7), F(7, 10)), (('A', 1), ('E', 8), F(1, 2)), (('A', 1), ('F', 4), F(6, 7)),
             (('A', 1), ('G', 2), F(1)), (('T', 1), ('E', 6), F(6, 5)), (('T', 1), ('E', 7), F(1)), (('T', 1), ('E', 8), F(8, 11)),
             (('T', 1), ('F', 4), F(6, 5)), (('T', 1), ('G', 2), F(4, 3)), (('E', 8), ('T', 1), F(80, 11)), (('T', 1), ('A', 1), F(3, 5)),
             (('T', 1), ('A', 2), F(1)), (('T', 1), ('A', 3), F(9, 7)), (('T', 2), ('T', 1), F(5, 4)), (('T', 2), ('A', 1), F(10, 7)),
             (('A', 2), ('T', 1), F(1)), (('A', 1), ('A', 2), F(4, 5)), (('A', 2), ('A', 1), F(6, 5)), (('A', 1), ('A', 3), F(1)),
             (('A', 1), ('B', 3), F(5, 4))]
    bad = [(X, Y, c0, cc(X, Y)) for X, Y, c0 in rows if cc(X, Y) != c0]
    print("Table 1, column c: %d entries (families instantiated for r<=7, k<=3) compared with (1.6); mismatches: %s" % (len(rows), bad))


SECTIONS = ['s_table', 's_T2A1', 's_T1D', 's_sec2', 's_A1E', 's_A1A', 's_A1C', 's_T1C', 's_FG', 's_E', 's_E8T1', 's_rank12',
            's_T2T1', 's_T1A3', 's_A1A2', 's_A1B3']

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    args = [int(a) for a in sys.argv[2:]]
    todo = SECTIONS if which == 'all' else [which]
    for name in todo:
        print("==== %s %s" % (name, args), flush=True)
        t0 = time.time()
        globals()[name](*args)
        print("==== %s done in %.1fs" % (name, time.time() - t0), flush=True)
    nbad = sum(1 for x in LOG if not x[1])
    print("SUMMARY: %d comparisons, %d failed" % (len(LOG), nbad))
