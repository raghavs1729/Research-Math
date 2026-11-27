"""Connection data for L4 = theta^4 - x(theta+1/2)^4 between the Frobenius basis at x=0 and holomorphic solutions at x=1,
and test of the pairing formulas for U2 and U3."""
from mpmath import mp, mpf, pi, zeta, gamma, log, matrix, lu_solve, pslq, taylor
import pickle
mp.dps = 110
K = 420           # series terms
x0 = mpf(1) / 2
# ---------- Frobenius at 0: Phi = sum_k A_k(eps) x^{k+eps}, truncated at eps^3
def pmul(a, b):
    r = [mpf(0)] * 4
    for i in range(4):
        for j in range(4 - i):
            r[i + j] += a[i] * b[j]
    return r
def pinv(a):      # 1/(a0 + a1 e + ...)
    r = [1 / a[0], 0, 0, 0]
    for n in range(1, 4):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r
half = mpf(1) / 2
A = [mpf(1), mpf(0), mpf(0), mpf(0)]
lx = log(x0)
xe = [lx**i / mp.factorial(i) for i in range(4)]          # x^eps
vals = [[mpf(0)] * 4 for _ in range(4)]                    # vals[m][j] = theta^m y_j (x0)
xp = mpf(1)
for k in range(K):
    # term: (k+eps)^m A_k(eps) x^k x^eps
    base = pmul(A, xe)
    cur = base
    for m in range(4):
        for j in range(4):
            vals[m][j] += cur[j] * xp
        cur = pmul(cur, [mpf(k), mpf(1), mpf(0), mpf(0)])
    f = pmul([k + half, mpf(1), 0, 0], pinv([mpf(k + 1), mpf(1), 0, 0]))
    f2 = pmul(f, f); A = pmul(A, pmul(f2, f2)); xp *= x0
# ---------- local holomorphic solutions at x=1 in s = 1-x: theta = -(1-s) d/ds
M = K + 8
def theta(c):
    return [(m * c[m] if m < len(c) else 0) - ((m + 1) * c[m + 1] if m + 1 < len(c) else 0) for m in range(len(c))]
def applyL(c):
    t1 = c; 
    th = [c]
    for _ in range(4): th.append(theta(th[-1]))
    # (theta+1/2)^4 c
    from math import comb
    sh = [sum(comb(4, j) * half**(4 - j) * th[j][m] for j in range(5)) for m in range(len(c))]
    # x = 1 - s
    xs = [sh[m] - (sh[m - 1] if m >= 1 else 0) for m in range(len(c))]
    return [th[4][m] - xs[m] for m in range(len(c))]
# build recurrence numerically: coefficient of s^p in L(sum c_m s^m) is linear in c_0..c_{p+3}; get the row by applying L to unit vectors lazily
def solve_series(c0, c1, c2):
    c = [mpf(c0), mpf(c1), mpf(c2)] + [mpf(0)] * (M - 3)
    for p in range(0, M - 4):
        # residual at power p with c_{p+3}=0, and derivative wrt c_{p+3}
        # only entries up to p+4 matter for power p: use a window
        w = c[:p + 5]
        w[p + 3] = mpf(0); w[p + 4] = mpf(0)
        r0 = applyL(w)[p]
        e = [mpf(0)] * (p + 5); e[p + 3] = mpf(1)
        d = applyL(e)[p]
        c[p + 3] = -r0 / d
    return c
def evalloc(c):
    s = 1 - x0
    out = []
    cur = c
    for m in range(4):
        out.append(sum(cur[i] * s**i for i in range(len(cur) - 6)))
        cur = theta(cur)
    return out
nu = solve_series(0, 1, mpf(3) / 4)
n3 = nu[3]
print('nu series: n2 = 3/4, n3 =', mp.nstr(n3, 30))
r2 = solve_series(0, 0, 1)
rstar = solve_series(1, 0, 3 * n3 - mpf(15) / 8)
r0a = solve_series(1, 0, 0)
Mat = matrix(4, 4)
for m in range(4):
    for j in range(4):
        Mat[m, j] = vals[m][j]
def coords(c):
    b = matrix(evalloc(c))
    return [x for x in lu_solve(Mat, b)]
c_nu, c_r2, c_rs, c_r0 = coords(nu), coords(r2), coords(rstar), coords(r0a)
print('coords nu  :', [mp.nstr(x, 25) for x in c_nu])
print('  expected nu from Mellin: coefficients of Res Gamma(-s)^4/Gamma(1/2-s)^4 x^s:  nu_3 =', mp.nstr(-6 / pi**2 * 0 + 1, 3), '(check below)')
print('coords r2  :', [mp.nstr(x, 25) for x in c_r2])
print('coords r*  :', [mp.nstr(x, 25) for x in c_rs])
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
print('eta_0 = c_r2[3]/2 =', mp.nstr(c_r2[3] / 2, 30), ' 16 L2/pi^2 =', mp.nstr(16 * L2 / pi**2, 30))
N = lambda e: (gamma(1 + e) / gamma(half + e))**4
Xi = lambda e: 2 * zeta(3, half + e) * N(e)
mp.dps = 140
Nt = taylor(N, 0, 3); Xt = taylor(Xi, 0, 3)
mp.dps = 100
def pairings(c):
    eta = [(-1)**i * c[3 - i] / 2 for i in range(4)]          # Y_rho(eps) = sum eta_i eps^i
    Ym = [eta[i] * (-1)**i for i in range(4)]                  # Y_rho(-eps)
    P1 = sum(Nt[1 - i] * Ym[i] for i in range(2))
    P3 = sum(Xt[3 - i] * Ym[i] for i in range(4))
    P3N = sum(Nt[3 - i] * Ym[i] for i in range(4))
    return P1, P3, P3N
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
U2, U3 = C['U2'][0], C['U3'][0]
for nm, c in (('r2', c_r2), ('r*', c_rs), ('nu', c_nu), ('r0a', c_r0)):
    P1, P3, P3N = pairings(c)
    print(nm, ' P1 =', mp.nstr(P1, 25), ' P3 =', mp.nstr(P3, 25), ' [e^3] N Y(-e) =', mp.nstr(P3N, 10), '  P3/(16 P1) =', mp.nstr(P3 / (16 * P1), 30))
print('1/U2     =', mp.nstr(1 / U2, 30))
print('1/U3 + 4 =', mp.nstr(1 / U3 + 4, 30))
pickle.dump({'nu': [str(x) for x in c_nu], 'r2': [str(x) for x in c_r2], 'rs': [str(x) for x in c_rs], 'r0a': [str(x) for x in c_r0], 'n3': str(n3)}, open('u3_coords.pkl', 'wb'))
