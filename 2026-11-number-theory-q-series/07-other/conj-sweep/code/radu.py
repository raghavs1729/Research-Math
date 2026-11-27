"""Radu's finite check for congruences of eta-quotient coefficients in arithmetic progressions.

Implements Lemma 4.5 of S. Radu, "An algorithmic approach to Ramanujan's congruences",
Ramanujan J. 20 (2009) 215-251, in the form quoted by Radu-Sellers (Int. J. Number Theory 7 (2011)
2249-2259, Lemma 2.4):

  f_r(q) = prod_{delta | M} (q^delta; q^delta)_inf^{r_delta} = sum c_r(n) q^n.
  If (m, M, N, t, r) is in Delta^*, a = (a_delta)_{delta | N} is an integer vector with
  p_{m,r}(gamma) + p^*_a(gamma) >= 0 for a complete set of representatives gamma of
  Gamma_0(N)\\SL2(Z)/Gamma_inf, and c_r(m n + t') = 0 (mod u) for all t' in P_{m,r}(t) and
  0 <= n <= floor(nu), then it holds for all n >= 0.

Everything here is exact (integers / Fractions). The coefficient check is done modulo u with FLINT
(nmod_poly); an independent check with PARI/GP is produced by radu_gp().
"""
from fractions import Fraction as Fr
from math import gcd, floor
import flint

def divisors(n): return [d for d in range(1, n+1) if n % d == 0]
def primes_of(n):
    ps = []; p = 2
    while p*p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps
def index_gamma0(N):
    x = Fr(N)
    for p in primes_of(N): x *= Fr(p+1, p)
    assert x.denominator == 1
    return int(x)
def kappa(m): return gcd(m*m - 1, 24)

def P_set(m, r, t):
    """P_{m,r}(t) = { t s + (s-1)/24 * sum(delta r_delta) mod m : s a square in (Z/24m)^* }."""
    S = sum(d*e for d, e in r.items())
    sq = set((x*x) % (24*m) for x in range(24*m) if gcd(x, 24*m) == 1)
    out = set()
    for s in sq:
        assert (s - 1) % 24 == 0
        out.add((t*s + (s-1)//24 * S) % m)
    return sorted(out)

def in_delta_star(m, M, N, t, r):
    k = kappa(m)
    if any(N % p for p in primes_of(m)): return False                      # (a)
    if any((m*N) % d for d, e in r.items() if e != 0): return False        # (b)
    if (k*N*sum(Fr(e*m*N, d) for d, e in r.items())) % 24 != 0: return False   # (c)
    if (k*N*sum(r.values())) % 8 != 0: return False                        # (d)
    S = sum(d*e for d, e in r.items())
    if N % (24*m // gcd(k*(-24*t - S), 24*m)) != 0: return False            # (e)
    if m % 2 == 0:                                                         # (f)
        prod = 1
        for d, e in r.items(): prod *= d**abs(e)
        s = 0
        while prod % 2 == 0: prod //= 2; s += 1
        j = prod
        if not (((k*N) % 4 == 0 and (N*s) % 8 == 0) or (s % 2 == 0 and (N*(1-j)) % 8 == 0)): return False
    return True

def cusp_reps(N):
    """A superset of representatives of Gamma_0(N)\\Gamma/Gamma_inf: all (a,c), c | N, 1<=a<=c, gcd(a,c)=1.
    (Every cusp of Gamma_0(N) is equivalent to some a/c with c | N, gcd(a,c)=1; listing more
    representatives than necessary only adds constraints.)"""
    return [(a, c) for c in divisors(N) for a in range(1, c+1) if gcd(a, c) == 1]

def p_mr(m, r, a, c):
    k = kappa(m)
    return min(sum(Fr(e * gcd(d*(a + k*lam*c), m*c)**2, 24*d*m) for d, e in r.items()) for lam in range(m))
def p_star(avec, c):
    return sum(Fr(e * gcd(d, c)**2, 24*d) for d, e in avec.items())

def best_a(m, N, r, bound=60):
    """Integer vector a on divisors of N with p_mr + p*_a >= 0 at all cusp representatives, with small nu.
    First the rational vector a* making all cusp inequalities tight is computed exactly (the matrix
    gcd(d,c)^2/(24 d) is invertible); then an ILP (scipy milp) searches integer vectors round(a*) + x,
    |x_i| <= bound, minimising sum a_d (index - d).  The result is re-verified exactly in certificate()."""
    import numpy as np, sympy
    from scipy.optimize import milp, LinearConstraint, Bounds
    ds = divisors(N); idx = index_gamma0(N)
    need = {}
    for a, c in cusp_reps(N):
        v = p_mr(m, r, a, c)
        need[c] = min(need.get(c, v), v)
    cs = list(need)
    Mx = sympy.Matrix([[sympy.Rational(gcd(d, c)**2, 24*d) for d in ds] for c in cs])
    astar = Mx.LUsolve(sympy.Matrix([-sympy.Rational(need[c].numerator, need[c].denominator) for c in cs]))
    base = [int(sympy.floor(x + sympy.Rational(1, 2))) for x in astar]
    A = [[float(Fr(gcd(d, c)**2, 24*d)) for d in ds] for c in cs]
    lb = [float(-need[c] - sum(Fr(gcd(d, c)**2, 24*d)*b0 for d, b0 in zip(ds, base))) - 1e-9 for c in cs]
    cost = np.array([idx - d for d in ds], dtype=float)
    res = milp(c=cost, constraints=LinearConstraint(np.array(A), lb=np.array(lb), ub=np.inf),
               integrality=np.ones(len(ds)), bounds=Bounds(-bound, bound), options={'time_limit': 120})
    assert res.x is not None, res.message
    return {d: b0 + int(round(x)) for d, b0, x in zip(ds, base, res.x)}

def nu_bound(m, M, N, t, r, avec):
    tmin = min(P_set(m, r, t))
    return (Fr(1, 24)*((sum(r.values()) + sum(avec.values()))*index_gamma0(N) - sum(d*e for d, e in avec.items()))
            - Fr(sum(d*e for d, e in r.items()), 24*m) - Fr(tmin, m))

def eta_series(r, D, u):
    """prod (q^d;q^d)^{r_d} mod u, truncated at q^D (nmod_poly)."""
    def f(d):
        c = [0]*D; k = 0
        while True:
            e1 = d*k*(3*k-1)//2; e2 = d*k*(3*k+1)//2
            if e1 >= D: break
            c[e1] = (c[e1] + (-1)**k) % u
            if k and e2 < D: c[e2] = (c[e2] + (-1)**k) % u
            k += 1
        return flint.nmod_poly(c, u)
    res = flint.nmod_poly([1], u)
    for d, e in r.items():
        if e == 0: continue
        p = f(d).pow_trunc(abs(e), D)
        if e < 0: p = p.inverse_series_trunc(D)
        res = res.mul_low(p, D)
    return res

def certificate(name, M, r, m, t, u, N=None, out=None, claimed=None):
    """Find N, a; verify hypotheses exactly; verify the congruence up to floor(nu). Returns True iff proved."""
    lines = []
    def P(*x):
        s = ' '.join(str(y) for y in x); lines.append(s); print(s, flush=True)
    r = {d: e for d, e in r.items() if e != 0}
    assert all(M % d == 0 for d in r)
    P(f'## {name}')
    P(f'eta-quotient exponents r = {r} (M = {M}); progression {m}n + {t}; modulus u = {u}')
    Pset = P_set(m, r, t)
    P(f'P_(m,r)(t) = {Pset}')
    if claimed is not None and not set(Pset) <= set(claimed):
        P(f'!! P contains residues outside the claimed set {claimed}'); 
    if N is None:
        N = next(n for n in range(1, 100000) if in_delta_star(m, M, n, t, r))
    assert in_delta_star(m, M, N, t, r)
    P(f'N = {N}: (m, M, N, t, r) in Delta^* verified; [SL2(Z):Gamma_0(N)] = {index_gamma0(N)}; kappa = {kappa(m)}')
    avec = best_a(m, N, r)
    reps = cusp_reps(N)
    worst = min(p_mr(m, r, a, c) + p_star(avec, c) for a, c in reps)
    assert worst >= 0
    P(f'a = { {d: e for d, e in avec.items() if e} }; min over {len(reps)} cusp representatives of p_mr + p*_a = {worst} >= 0 (exact)')
    nu = nu_bound(m, M, N, t, r, avec); B = floor(nu)
    P(f'nu = {nu} ; floor(nu) = {B}')
    D = m*(B+1) + 1
    ser = eta_series(r, D, u)
    ok = True; nz = 0
    for tp in Pset:
        for n in range(0, B+1):
            if int(ser[m*n + tp]) % u != 0:
                ok = False; P(f'!! c_r({m}*{n}+{tp}) = {int(ser[m*n+tp])} mod {u} != 0'); break
    P(f'checked c_r({m}n + t\') = 0 mod {u} for t\' in P and 0 <= n <= {B}: {"ALL ZERO -> congruence PROVED for all n" if ok else "FAILED"}')
    if out:
        with open(out, 'a') as fh: fh.write('\n'.join(lines) + '\n\n')
    return ok, (N, avec, B, Pset)
