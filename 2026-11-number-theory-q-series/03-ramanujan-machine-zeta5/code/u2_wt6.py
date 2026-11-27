"""Test: are the U2/U3 constants critical L-values of weight-6 cusp forms of level 4 or 8 (times powers of pi and L(f8,j))?"""
from mpmath import mp, mpf, pi, gamma, zeta, pslq, gammainc, sqrt
import pickle, sympy as sp
mp.dps = 300
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
N = 900
def eta_prod(spec):          # spec: list of (m, e): prod eta(m tau)^e, returns q-exp coefficients list (index = power), needs integer leading power
    lead = sum(m * e for m, e in spec)
    assert lead % 24 == 0
    lead //= 24
    co = [0] * (N + 1); co[0] = 1
    for m, e in spec:
        for n in range(1, N // m + 1):
            st = m * n
            for _ in range(abs(e)):
                if e > 0:
                    for i in range(N, st - 1, -1):
                        co[i] -= co[i - st]
                else:
                    for i in range(st, N + 1):
                        co[i] += co[i - st]
    return [0] * lead + co[:N + 1 - lead]
f8 = eta_prod([(2, 4), (4, 4)])
sig = lambda n: sum(d for d in range(1, n + 1) if n % d == 0)
def E2diff(M):               # E2(tau) - M E2(M tau), q-exp, normalised: (M-1)/... keep integers: 24*sum sigma(n)(q^n - M q^{Mn}) + (M-1)
    r = [0] * (N + 1); r[0] = M - 1
    for n in range(1, N + 1):
        r[n] += 24 * sig(n)
        if n * M <= N: r[n * M] -= 24 * M * sig(n)
    return r
def mul(A, B):
    R = [0] * (N + 1)
    for i, x in enumerate(A):
        if x:
            for j in range(N + 1 - i):
                if B[j]: R[i + j] += x * B[j]
    return R
basis = [mul(f8, E2diff(M)) for M in (2, 4, 8)]
# Hecke T3 on weight 6, level 8: (T3 f)_n = a_{3n} + 3^5 a_{n/3}
def T3(f):
    return [f[3 * n] + (243 * f[n // 3] if n % 3 == 0 else 0) for n in range(0, N // 3)]
K = 12
M = sp.Matrix([[b[n] for n in range(1, K)] for b in basis])           # rows = basis, cols = coefficients 1..K-1
MT = sp.Matrix([[T3(b)[n] for n in range(1, K)] for b in basis])
X = MT * M.pinv()        # T3 in this basis (row convention): MT = X M
X = sp.simplify(X)
print('T3 eigenvalues:', X.eigenvals())
new = None
for val, mult, vecs in X.T.eigenvects():
    pass
ev = X.left_eigenvects() if hasattr(X, 'left_eigenvects') else None
for val, mult, vecs in X.T.eigenvects():
    if val != -12:
        v = vecs[0]
        g = [sum(v[i] * basis[i][n] for i in range(3)) for n in range(N + 1)]
        g = [x / g[1] for x in g]
        print('newform candidate eigenvalue', val, 'coeffs', g[1:12])
        new = [int(x) for x in g]
eta12 = eta_prod([(2, 12)])
def Lval(a, s, k, lev, eps):
    A = sqrt(mpf(lev)) / (2 * pi); tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        x = n / A
        tot += a[n] * (gammainc(s, x) / (mpf(n) / A)**s + eps * gammainc(k - s, x) / (mpf(n) / A)**(k - s))
    return tot / (A**s * gamma(s))
forms = {'eta(2t)^12 (N=4)': (eta12, 4), 'new N=8': (new, 8)}
for nm, (a, lev) in forms.items():
    # determine sign of functional equation numerically
    for eps in (1, -1):
        v1 = Lval(a[:400], mpf(7), 6, lev, eps); v2 = sum(mpf(a[n]) / mpf(n)**7 for n in range(1, N))
        if abs(v1 - v2) < mpf(10)**-6:
            good = eps
    Lg = {s: Lval(a, s, 6, lev, good) for s in (1, 2, 3, 4, 5)}
    print(nm, 'eps', good, [mp.nstr(Lg[s], 15) for s in (1, 2, 3, 4, 5)], flush=True)
    for name in ('U2', 'U3'):
        L, al1, be1, al2, be2 = C[name]
        targets = {'1/U': 1 / L, 'c_q': be1, 'c_q/L1': be1 / L1, 'c_q/L2': be1 / L2, 'c_q*L1': be1 * L1, 'c_q*L2': be1 * L2, 'al2': al2, 'be2': be2}
        for tn, T in targets.items():
            for a_ in range(-8, 9):
                for bas in ([Lg[5], Lg[4] * pi], [Lg[5], Lg[4] * pi, zeta(5), pi**5], [Lg[5], Lg[4] * pi, pi**3 * zeta(3) / pi**0, pi**5]):
                    r = pslq([T * pi**a_] + bas, maxcoeff=10**8, maxsteps=50000)
                    if r and r[0] != 0:
                        print('  HIT', nm, name, tn, 'a=', a_, len(bas), r, flush=True)
print('done')
