"""Certificate for proofs/F-all-k.md  (KRTW (F_k) for all k, hence identities (ii), (iii)).

Every check prints OK/FAIL; the script ends with ALL CHECKS PASSED or a failure count.
Eigenfunctions are computed from the triangular systems of T and T# themselves (tower.py, sharp.py),
NOT from the H'/G recursions, so the eigen-relations checked in C3/C5 are independent confirmations.

  C1  [T,H'] = [T,G] = 0, [T#,H'#] = [T#,G#] = 0 (point masses, exact, two rational points) and
      H'# = W^-1 H'^t W, G# = W^-1 G^t W as kernels
  C2  Bailey-basis identities (symbolic in y = q^n, z = q^A): tau phi_n tridiagonal, column sums of tau, S, q^{2A}S
  C3  T phi = lam phi, H' phi = eps phi, G phi = g phi in position space (exact, rational point)
  C4  value at infinity of phi_{N,M}: product formula V_{N,M}; recursion denominator factorisation
  C5  T# phi# = lam phi#, H'# phi# = eps phi#, G# phi# = g phi#; product formula V#_{N,M}; Fricke factorisation
  C6  pairing <psi_a, phi_b> (full finite double sum) = delta_{ab} * single top term
  C7  product identity V V# (x^2q^{2L+2};q^2)_{2M} (xq;q)_{2N} = t_{N,M}  (symbolic)
  C8  K_k(x) (direct summation) = sum_{N>=M} lam^k t_{N,M} / ((xq)_inf (x^2q^2;q^2)_inf)  mod q^QMAX, all x-powers
  C9  u = sum_{N>=M} u~_{N,M} phi_{N,M} at points (R,S), mod x^XMAX, rational q
"""
import sys
from fractions import Fraction as Fr
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from tower import tower, lam
from sharp import eig_sharp
from check_II import target, betaB, betaJ
import verify_ops as VO

FAILS = 0


def report(label, ok, extra=''):
    global FAILS
    print('[%s] %s %s' % (label, 'OK' if ok else 'FAIL', extra))
    sys.stdout.flush()
    if not ok:
        FAILS += 1


def ev_poly(p, q0, x0):
    s = Fr(0)
    for (i, j), c in p.to_dict().items():
        s += Fr(int(c)) * x0 ** int(i) * q0 ** int(j)
    return s


def ev_rf(f, q0, x0):
    return ev_poly(f.n, q0, x0) / ev_poly(f.d, q0, x0)


def eps(N, M):
    return mono(0, -N - M) * RF(ONE + X * Q ** (2 * N)) * RF(ONE + X * Q ** (2 * M))


def gval(N, M):
    return mono(0, -2 * N) + mono(2, 2 * N) + mono(0, -2 * M) + mono(2, 2 * M)


def Vprod(N, M):
    return mono(2 * M, 2 * M * N + M * (2 * M - 1)) * poch(-X, M) / (
        qpoch(M, 2) * poch(Q ** (2 * (N - M + 1)), M, 2) * poch(X * Q ** M, M) * poch(X * X * Q ** (2 * N), M, 2))


def Vsprod(N, M):
    L = N - M
    if L == 0:
        return RF(1)
    s = -1 if L % 2 else 1
    return RF(s) * mono(L, L * (2 * L - 1) + 2 * M * L) * RF(ONE + Q ** L) * poch(X * X * Q ** 2, L, 2) / (
        qpoch(L, 2) * poch(X * Q ** (2 * M + 1), L) * poch(X * Q ** (2 * M + L), L))


# ---------------------------------------------------------------- eigenfunctions (genuine normalisation)
def phi_coeffs(N, M):
    """R-normalisation: slice S has coefficients c[S][e] relative to x q^{2S}; top c[M][L] = 1."""
    L = N - M
    out = {}
    for S in range(M + 1):
        g = tower(L, M - S)
        xs = X * Q ** (2 * S)
        fac = poch(X * X * Q ** (4 * S + 2), 2 * (M - S), 2)
        out[S] = {e: v.subs_x(xs) / fac for e, v in g.items()}
    return out


def phi_val(c, R, S, q0, x0):
    """phi(R,S) = sum_e c[S][e] / ((q)_{R-e} (x q^{2S+1})_{R+e})"""
    if S not in c:
        return Fr(0)
    tot = Fr(0)
    for e, v in c[S].items():
        if R < e:
            continue
        den = VO.qp(q0, R - e, q0) * VO.qp(x0 * q0 ** (2 * S + 1), R + e, q0)
        tot += ev_rf(v, q0, x0) / den
    return tot


def phis_val(cs, R, S, q0, x0):
    """phi#(R,S) = sum_e c[(e,R)] / ((Q;Q)_{S-e} (x^2 q^{2R} Q;Q)_{S+e})"""
    tot = Fr(0)
    Q0 = q0 * q0
    for (e, rho), v in cs.items():
        if rho != R or S < e:
            continue
        den = VO.qp(Q0, S - e, Q0) * VO.qp(x0 * x0 * q0 ** (2 * R) * Q0, S + e, Q0)
        tot += ev_rf(v, q0, x0) / den
    return tot


def main():
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    q0, x0 = Fr(3, 7), Fr(5, 11)

    # ---------------- C1
    tot = 0
    for (qq, xx) in [(Fr(3, 7), Fr(5, 11)), (Fr(-2, 9), Fr(13, 4))]:
        O = VO.Ops(qq, xx, RMAX=14, SMAX=9)
        for r0 in range(4):
            for s0 in range(4):
                d = {(r0, s0): Fr(1)}
                tot += VO.window_diff(O.T(O.Hp(d)), O.Hp(O.T(d)), 10, 5)
                tot += VO.window_diff(O.T(O.Gop(d)), O.Gop(O.T(d)), 10, 5)
                tot += VO.window_diff(O.Ts(O.Hps(d)), O.Hps(O.Ts(d)), r0 + 2, 5)
                tot += VO.window_diff(O.Ts(O.Gs(d)), O.Gs(O.Ts(d)), r0 + 2, 5)
    report('C1a', tot == 0, 'commutators vanish on 2*16 point masses (mismatches: %d)' % tot)
    # transposition: Hs(a,b) = H(b,a) w(b)/w(a), where K(a,b) := (K delta_b)(a)
    O = VO.Ops(q0, x0, RMAX=8, SMAX=8)
    bad = 0
    pts = [(r, s) for r in range(5) for s in range(5)]
    for b in pts:
        db = {b: Fr(1)}
        Hb = O.Hp(db); Gb = O.Gop(db)
        Hsb = O.Hps(db); Gsb = O.Gs(db)
        for a in pts:
            da = {a: Fr(1)}
            Ha = O.Hp(da); Ga = O.Gop(da)
            # H#(a,b) = H(b,a) w(b)/w(a)
            if Hsb.get(a, 0) != Ha.get(b, 0) * O.w(*b) / O.w(*a):
                bad += 1
            if Gsb.get(a, 0) != Ga.get(b, 0) * O.w(*b) / O.w(*a):
                bad += 1
    report('C1b', bad == 0, "H'# = W^-1 H'^t W and G# = W^-1 G^t W on 25x25 kernel entries (mismatches %d)" % bad)
    # T# = W^-1 T^t W
    bad = 0
    for b in pts:
        Tb = O.T({b: Fr(1)}); Tsb = O.Ts({b: Fr(1)})
        for a in pts:
            Ta = O.T({a: Fr(1)})
            if Tsb.get(a, 0) != Ta.get(b, 0) * O.w(*b) / O.w(*a):
                bad += 1
    report('C1c', bad == 0, 'T# = W^-1 T^t W on 25x25 kernel entries (mismatches %d)' % bad)

    # ---------------- C2 (symbolic, y = q^n, z = q^A): ratios phi_{n-1}(A)/phi_n(A+1) etc.
    import flint
    C4 = flint.fmpz_mpoly_ctx.get(('x', 'q', 'y', 'z'), 'lex')
    x4, q4, y4, z4 = C4.gens()
    one = C4.from_dict({(0, 0, 0, 0): 1})
    # identity  z * q*(1-x y^2/q)(1-x y^2)(1-x y^2 q) * y^? ... clear denominators:
    # alpha = y/(q (1-xy^2/q)(1-xy^2)), beta = -y(1+xy^2)/(q(1-xy^2/q)(1-xy^2 q)), gamma = x y^3/((1-xy^2)(1-xy^2q))
    # multiply by Dn := q*(q - x y^2)(1 - x y^2)(1 - x y^2 q) / q   (note q(1-xy^2/q) = q - x y^2)
    A1 = one - x4 * z4 * y4
    A2 = one - x4 * z4 * y4 * q4
    Bz = q4 * one * 0
    # basis numerators (times y^2 to clear z/y):  (1 - z q / y) -> (y - z q)/y ; (1 - z/y) -> (y - z)/y
    # identity: z = a A1 A2 + b (1-zq/y) A2 + c (1-z/y)(1-zq/y); multiply by y^2:
    # y^2 z = a y^2 A1 A2 + b y (y - zq) A2 + c (y - z)(y - zq)
    # a = y/((q - x y^2)(1 - x y^2)), b = -y(1+xy^2)/((q - x y^2)(1 - x y^2 q)), c = x y^3/((1-xy^2)(1-xy^2 q))
    D = (q4 - x4 * y4 ** 2) * (one - x4 * y4 ** 2) * (one - x4 * y4 ** 2 * q4)
    aN = y4 * (one - x4 * y4 ** 2 * q4)            # a * D
    bN = -y4 * (one + x4 * y4 ** 2) * (one - x4 * y4 ** 2)
    cN = x4 * y4 ** 3 * (q4 - x4 * y4 ** 2)
    lhs = y4 ** 2 * z4 * D
    rhs = aN * y4 ** 2 * A1 * A2 + bN * y4 * (y4 - z4 * q4) * A2 + cN * (y4 - z4) * (y4 - z4 * q4)
    report('C2a', (lhs - rhs).is_zero(), 'q^A phi_n(A+1) = alpha phi_{n-1} + beta phi_n + gamma phi_{n+1} (symbolic in q^n, q^A)')
    report('C2b', (aN + bN + cN).is_zero(), 'alpha_n + beta_n + gamma_n = 0 (symbolic)')
    # S phi_n and q^{2A} S phi_n: solve the 3x3 systems symbolically via sympy and check column sums
    import sympy as sp
    xs, qs, ys, zs = sp.symbols('x q y z')
    a, b, c = sp.symbols('a b c')
    B1 = (1 - xs * zs * ys) * (1 - xs * zs * ys * qs)
    B2 = (1 - zs * qs / ys) * (1 - xs * zs * ys * qs)
    B3 = (1 - zs / ys) * (1 - zs * qs / ys)
    for tgt, want, name in [(sp.Integer(1), 1, 'S'), (zs ** 2, 0, 'q^{2A} S')]:
        expr = sp.expand((a * B1 + b * B2 + c * B3 - tgt) * ys ** 2)
        eqs = sp.Poly(expr, zs).all_coeffs()
        sol = sp.solve(eqs, [a, b, c], dict=True)
        ok = len(sol) == 1 and sp.simplify(sol[0][a] + sol[0][b] + sol[0][c] - want) == 0
        report('C2c', ok, '%s phi_n is tridiagonal with column sum %d (symbolic)' % (name, want))

    # base-Q partial fractions used in section 5 (A := Q^S, symbols: a, Q, Qe := Q^e)
    aa, QQ, Qe, Zs = sp.symbols('a Q Qe Z')
    # phi^Q_e(S;a) ~ 1/((Q)_{S-e}(aQ)_{S+e}); common denominators: ratios as polynomials in Z = Q^S
    # (i) phi^Q_e(S;aQ) = c1 phi^Q_e(S;a) + c2 phi^Q_{e+1}(S;a), c1 + c2 = 1 - aQ
    c1, c2, c3 = sp.symbols('c1 c2 c3')
    # multiply by (Q)_{S-e}(aQ)_{S+e+1}:  (1-aQ) = c1 (1 - a Qe Z Q) + c2 (1 - Z/Qe)
    expr = sp.expand((1 - aa * QQ) - c1 * (1 - aa * Qe * Zs * QQ) - c2 * (1 - Zs / Qe))
    sol = sp.solve(sp.Poly(expr, Zs).all_coeffs(), [c1, c2], dict=True)
    ok_i = len(sol) == 1 and sp.simplify(sol[0][c1] + sol[0][c2] - (1 - aa * QQ)) == 0
    # (ii) Q^S phi^Q_e(S;aQ): (1-aQ) Z = b (1 - a Qe Z Q) - b (1 - Z/Qe), b = (1-aQ)/(Qe^{-1} - a Qe Q)
    bb = (1 - aa * QQ) / (1 / Qe - aa * Qe * QQ)
    ok_ii = sp.simplify((1 - aa * QQ) * Zs - bb * (1 - aa * Qe * Zs * QQ) + bb * (1 - Zs / Qe)) == 0
    # (iii) Q^S phi^Q_e(S;aQ^2) relative to a: (aQ;Q)_2 Z = c1 (1-aQe Q Z)(1-aQe Q^2 Z) + c2 (1-Z/Qe)(1-aQe Q^2 Z) + c3 (1-Z/Qe)(1-Z/(Qe Q))
    expr = sp.expand((1 - aa * QQ) * (1 - aa * QQ ** 2) * Zs - c1 * (1 - aa * Qe * QQ * Zs) * (1 - aa * Qe * QQ ** 2 * Zs)
                     - c2 * (1 - Zs / Qe) * (1 - aa * Qe * QQ ** 2 * Zs) - c3 * (1 - Zs / Qe) * (1 - Zs / (Qe * QQ)))
    sol = sp.solve(sp.Poly(expr, Zs).all_coeffs(), [c1, c2, c3], dict=True)
    ok_iii = len(sol) == 1 and sp.simplify(sol[0][c1] + sol[0][c2] + sol[0][c3]) == 0
    report('C2d', ok_i and ok_ii and ok_iii, 'base-q^2 partial fractions of section 5 (sums 1-aQ, 0, 0), symbolic')

    # ---------------- C3, C4
    O = VO.Ops(q0, x0, RMAX=16, SMAX=10)
    ok3 = True; ok4 = True
    for N in range(0, 5):
        for M in range(0, N + 1):
            c = phi_coeffs(N, M)
            L = N - M
            RW = L + 6
            F = {(R, S): phi_val(c, R, S, q0, x0) for R in range(0, RW + 3) for S in range(0, M + 1)}
            F = {k: v for k, v in F.items() if v != 0}
            lm = ev_rf(lam(N, M), q0, x0)
            TF = O.T(F); HF = O.Hp(F); GF = O.Gop(F)
            e0 = ev_rf(eps(N, M), q0, x0); g0 = ev_rf(gval(N, M), q0, x0)
            for R in range(0, RW + 1):
                for S in range(0, M + 2):
                    f = F.get((R, S), Fr(0))
                    if TF.get((R, S), 0) != lm * f or HF.get((R, S), 0) != e0 * f or GF.get((R, S), 0) != g0 * f:
                        ok3 = False
            # support: phi(R,S) = 0 for R < L
            for S in range(M + 1):
                for R in range(L):
                    if phi_val(c, R, S, q0, x0) != 0:
                        ok3 = False
            v = sum(c[0].values(), RF(0))
            if not (v == Vprod(N, M) * poch(X * Q, 2 * M)):
                ok4 = False
    report('C3', ok3, 'T phi = lam phi, H\' phi = eps phi, G phi = g phi, support R >= N-M (N <= 4, exact at q=3/7,x=5/11)')
    for N in range(5, 8):
        for M in range(0, N + 1):
            c = phi_coeffs(N, M)
            if not (sum(c[0].values(), RF(0)) == Vprod(N, M) * poch(X * Q, 2 * M)):
                ok4 = False
    report('C4a', ok4, 'v(phi_{N,M}(.,0)) = (xq;q)_{2M} V_{N,M} (product formula), N <= 7 (symbolic)')
    ok = True
    for N in range(1, 9):
        for M in range(1, N + 1):
            e = eps(N, M); g = gval(N, M)
            den = e * e + RF((ONE + X) ** 2) * (RF((ONE - X) ** 2) - g)
            fac = mono(0, -2 * (N + M)) * RF(ONE - Q ** (2 * N)) * RF(ONE - Q ** (2 * M)) * RF(ONE - X * X * Q ** (2 * N)) * RF(ONE - X * X * Q ** (2 * M))
            r = mono(2, 0) * RF(ONE + X) * (RF(X) + mono(0, -1)) / den
            if not (den == fac and Vprod(N, M) == r * Vprod(N - 1, M - 1).subs_x(X * Q ** 2)):
                ok = False
    report('C4b', ok, 'eps^2+(1+x)^2((1-x)^2-g) factorisation and V-recursion, N <= 8 (symbolic)')

    # ---------------- C5
    ok5 = True
    for N in range(0, 5):
        for M in range(0, N + 1):
            cs = eig_sharp(N, M)
            L = N - M
            SW = M + 5
            G = {(R, S): phis_val(cs, R, S, q0, x0) for R in range(0, L + 1) for S in range(0, SW + 4)}
            G = {k: v for k, v in G.items() if v != 0}
            lm = ev_rf(lam(N, M), q0, x0)
            TG = O.Ts(G); HG = O.Hps(G); GG = O.Gs(G)
            e0 = ev_rf(eps(N, M), q0, x0); g0 = ev_rf(gval(N, M), q0, x0)
            for R in range(0, L + 3):
                for S in range(0, SW + 1):
                    f = G.get((R, S), Fr(0))
                    if TG.get((R, S), 0) != lm * f or HG.get((R, S), 0) != e0 * f or GG.get((R, S), 0) != g0 * f:
                        ok5 = False
            for R in range(0, L + 1):
                for S in range(0, M):
                    if phis_val(cs, R, S, q0, x0) != 0:
                        ok5 = False
    report('C5a', ok5, "T# phi# = lam phi#, H'# phi# = eps phi#, G# phi# = g phi#, support S >= M (N <= 4, exact)")
    ok = True
    for N in range(0, 8):
        for M in range(0, N + 1):
            cs = eig_sharp(N, M)
            if not (sum((v for (e, rho), v in cs.items() if rho == 0), RF(0)) == Vsprod(N, M)):
                ok = False
    report('C5b', ok, 'v#(phi#_{N,M}(0,.)) = V#_{N,M} (product formula), N <= 7 (symbolic)')
    ok = True
    for N in range(1, 9):
        for M in range(0, N):
            L = N - M; Kp = N + M
            e = eps(N, M); g = gval(N, M)
            for rho in range(0, L):
                cr = mono(0, -rho) + mono(0, rho)
                den = g + RF(X) * cr * cr - e * cr
                fac = mono(0, -L - Kp) * RF(ONE - Q ** (L - rho)) * RF(ONE - Q ** (L + rho)) * RF(ONE - X * Q ** (Kp - rho)) * RF(ONE - X * Q ** (Kp + rho))
                if not (den == fac):
                    ok = False
    report('C5c', ok, 'Fricke factorisation g + x c^2 - eps c = q^{-L-N-M}(1-q^{L-r})(1-q^{L+r})(1-xq^{N+M-r})(1-xq^{N+M+r}), N <= 8')

    # ---------------- C6 pairing / biorthogonality (exact at rational point; sums are finite)
    ok = True
    modes = [(N, M) for N in range(0, 4) for M in range(0, N + 1)]
    data = {}
    for (N, M) in modes:
        data[(N, M)] = (phi_coeffs(N, M), eig_sharp(N, M))
    for a in modes:
        ca, csa = data[a]
        for b in modes:
            cb, csb = data[b]
            La, Ma = a[0] - a[1], a[1]
            Lb, Mb = b[0] - b[1], b[1]
            s = Fr(0)
            for r in range(Lb, La + 1):
                for S in range(Ma, Mb + 1):
                    s += O.w(r, S) * phis_val(csa, r, S, q0, x0) * phi_val(cb, r, S, q0, x0)
            if a == b:
                want = ev_rf(lam(*a), q0, x0) / (VO.qp(x0 * x0 * q0 ** (2 * La + 2), 2 * Ma, q0 * q0) * VO.qp(x0 * q0 ** (2 * Ma + 1), 2 * La, q0))
                if s != want:
                    ok = False
            elif s != 0:
                ok = False
    report('C6', ok, '<psi_a, phi_b> = delta_ab lam/((x^2q^{2L+2};q^2)_{2M}(xq^{2M+1};q)_{2L}) for all 10x10 pairs with N <= 3')

    # ---------------- C7 product identity
    ok = True
    for N in range(0, 9):
        for M in range(0, N + 1):
            L = N - M
            lhs = Vprod(N, M) * Vsprod(N, M) * poch(X * X * Q ** (2 * L + 2), 2 * M, 2) * poch(X * Q, 2 * N)
            if not (lhs == target(N, M)):
                ok = False
    report('C7', ok, 'V V# (x^2q^{2L+2};q^2)_{2M} (xq;q)_{2N} = t_{N,M}, N <= 8 (symbolic)')

    # ---------------- C8 (F_k) through the spectral formula, as series in q with polynomial coefficients in x
    from flint import fmpz_poly
    ok = check_C8(QMAX)
    report('C8', ok, 'K_k(x) = sum lam^k t/((xq)(x^2q^2;q^2)) for k <= 6, mod q^%d, all powers of x' % QMAX)

    # ---------------- C9 expansion of u
    ok = check_C9(q0)
    report('C9', ok, 'u(R,S) = 1/(q)_R = sum u~ phi_{N,M}(R,S) for R,S <= 2 mod x^7 at q = 3/7')

    print('ALL CHECKS PASSED' if FAILS == 0 else 'FAILURES: %d' % FAILS)


# ------------------------------------------------------------------ C8 machinery: bivariate truncated series
class BS:
    """dict (i,j) -> int : sum c x^i q^j, truncated q^j < QM"""
    QM = 30

    def __init__(self, d=None):
        self.d = d or {}

    def __add__(self, o):
        r = dict(self.d)
        for k, v in o.d.items():
            r[k] = r.get(k, 0) + v
        return BS({k: v for k, v in r.items() if v})

    def __mul__(self, o):
        r = {}
        for (i1, j1), v1 in self.d.items():
            for (i2, j2), v2 in o.d.items():
                j = j1 + j2
                if j < BS.QM:
                    k = (i1 + i2, j)
                    r[k] = r.get(k, 0) + v1 * v2
        return BS({k: v for k, v in r.items() if v})


def bs_from_rf(f):
    """expand a rational function n/d with d(0..,q) such that d has constant term 1 in q... : use q-adic expansion.
    Here: d is a product of factors (1 - x^a q^b) with b >= 1 or (1 - q^b) etc.  Implemented by series inversion in q
    over Z[x] (requires d(x, q=0) = +-1)."""
    num = {(int(i), int(j)): int(c) for (i, j), c in f.n.to_dict().items()}
    den = {(int(i), int(j)): int(c) for (i, j), c in f.d.to_dict().items()}
    # group by q-degree: den = sum_j D_j(x) q^j with D_0 = +-1
    D = {}
    for (i, j), c in den.items():
        D.setdefault(j, {})
        D[j][i] = D[j].get(i, 0) + c
    D0 = D.get(0, {})
    assert set(D0) == {0} and abs(D0[0]) == 1, 'denominator must be +-1 at q=0: %s' % f
    s0 = D0[0]
    Nn = {}
    for (i, j), c in num.items():
        Nn.setdefault(j, {})
        Nn[j][i] = Nn[j].get(i, 0) + c
    # solve F * D = N  for F = sum_j F_j(x) q^j
    Fq = {}
    for j in range(BS.QM):
        acc = dict(Nn.get(j, {}))
        for jj in range(1, j + 1):
            if jj in D and (j - jj) in Fq:
                for i1, c1 in D[jj].items():
                    for i2, c2 in Fq[j - jj].items():
                        acc[i1 + i2] = acc.get(i1 + i2, 0) - c1 * c2
        Fq[j] = {i: c * s0 for i, c in acc.items() if c}
    out = {}
    for j, row in Fq.items():
        for i, c in row.items():
            if c:
                out[(i, j)] = c
    return BS(out)


def check_C8(QM):
    BS.QM = QM
    import itertools
    ok = True
    # (xq;q)_inf (x^2 q^2; q^2)_inf
    P = BS({(0, 0): 1})
    for j in range(1, QM):
        P = P * BS({(0, 0): 1, (1, j): -1})
    for j in range(1, QM // 2 + 1):
        P = P * BS({(0, 0): 1, (2, 2 * j): -1})
    for k in range(1, 7):
        # spectral side: sum_{N>=M} lam^k t_{N,M}
        spec = BS()
        for N in range(0, QM):
            if (k + 2) * N * N - N >= QM:
                break
            for M in range(0, N + 1):
                if (N * N + M * M) * (k + 2) - N - M >= QM:
                    break
                tt = lam(N, M) ** k * target(N, M)
                spec = spec + bs_from_rf(tt)
        # direct side: (xq)(x^2q^2) K_k
        Kk = direct_K(k, QM)
        lhs = Kk * P
        if {kk: v for kk, v in lhs.d.items()} != spec.d:
            ok = False
            print('   C8 mismatch at k =', k)
    return ok


def direct_K(k, QM):
    """K_k(x) = sum x^{sum(r_i+2s_i)} prod q^{(r_i+s_i)^2+s_i^2}/((q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}}) mod q^QM,
    level-by-level transfer (exact integer arithmetic)."""
    # inverse Pochhammers as q-series
    def inv_poch(n, step):
        s = BS({(0, 0): 1})
        for i in range(1, n + 1):
            # 1/(1-q^{step i}) = sum q^{step i m}
            g = BS({(0, step * i * m): 1 for m in range((QM - 1) // (step * i) + 1)})
            s = s * g
        return s
    cacheq = {}; cacheQ = {}
    def iq(n):
        if n not in cacheq:
            cacheq[n] = inv_poch(n, 1)
        return cacheq[n]
    def iQ(n):
        if n not in cacheQ:
            cacheQ[n] = inv_poch(n, 2)
        return cacheQ[n]
    rmax = 0
    while (rmax + 1) ** 2 < QM:
        rmax += 1
    # F(R,S) := partial sum of the inner levels; start: u(r,s) = 1/(q)_r
    F = {(r, s): iq(r) for r in range(rmax + 1) for s in range(rmax + 1)}
    for lev in range(k):
        newF = {}
        last = (lev == k - 1)
        for R in range(rmax + 1 if not last else 1):
            for S in range(rmax + 1 if not last else 1):
                acc = BS()
                for r in range(0, (R if not last else rmax) + 1):
                    for s in range(S, rmax + 1):
                        e = (r + s) ** 2 + s * s
                        if e >= QM:
                            continue
                        val = F.get((r, s))
                        if val is None or not val.d:
                            continue
                        term = BS({(r + 2 * s, e): 1}) * val
                        if not last:
                            term = term * iq(R - r)
                        term = term * iQ(s - S)
                        acc = acc + term
                newF[(R, S)] = acc
        F = newF
    return F[(0, 0)]


def check_C9(q0):
    """expansion of u at points (R,S): u~_{N,M} = psi(u)/<psi,phi> = V# (x^2q^{2L+2};Q)_{2M} (xq^{2M+1};q)_{2L} / (x^2q^2;q^2)_inf"""
    XM = 7
    import flint
    def ser(f):
        a = {}
        for (i, j), c in f.n.to_dict().items():
            a[int(i)] = a.get(int(i), Fr(0)) + Fr(int(c)) * q0 ** int(j)
        b = {}
        for (i, j), c in f.d.to_dict().items():
            b[int(i)] = b.get(int(i), Fr(0)) + Fr(int(c)) * q0 ** int(j)
        A = [a.get(i, Fr(0)) for i in range(XM)]
        B = [b.get(i, Fr(0)) for i in range(XM)]
        inv = [Fr(0)] * XM
        inv[0] = 1 / B[0]
        for t in range(1, XM):
            inv[t] = -sum(B[j] * inv[t - j] for j in range(1, t + 1)) / B[0]
        return [sum(A[j] * inv[t - j] for j in range(t + 1)) for t in range(XM)]
    def mul(u, v):
        return [sum(u[j] * v[t - j] for j in range(t + 1)) for t in range(XM)]
    # 1/(x^2q^2;q^2)_inf as series in x at q0
    P = [Fr(0)] * XM
    for j in range(0, XM // 2 + 1):
        if 2 * j < XM:
            c = q0 ** (2 * j)
            for i in range(1, j + 1):
                c /= (1 - q0 ** (2 * i))
            P[2 * j] += c
    ok = True
    for R in range(0, 3):
        for S in range(0, 3):
            tot = [Fr(0)] * XM
            for M in range(S, XM):
                for L in range(0, R + 1):
                    N = L + M
                    if L + 2 * (M - S) >= XM:
                        continue
                    c = phi_coeffs(N, M)
                    if S not in c:
                        continue
                    ut = Vsprod(N, M) * poch(X * X * Q ** (2 * L + 2), 2 * M, 2) * poch(X * Q ** (2 * M + 1), 2 * L)
                    val = RF(0)
                    for e, v in c[S].items():
                        if R >= e:
                            val = val + v / (qpoch(R - e) * poch(X * Q ** (2 * S + 1), R + e))
                    tot = [tot[i] + w for i, w in enumerate(mul(ser(ut * val), P))]
            want = Fr(1) / VO.qp(q0, R, q0)
            if not (tot[0] == want and all(t == 0 for t in tot[1:])):
                ok = False
                print('   C9 mismatch at', R, S, tot[:4])
    return ok


if __name__ == '__main__':
    main()
