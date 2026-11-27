"""Guillera's level-3 shifted series at z = -9/16 (G13 in survey/verify_G.py; arXiv:1802.02023, Tables 3-4):
       S := sum_{n>=0} (10n+7) (1)_n (5/6)_n (7/6)_n / (3/2)_n^3 (-9/16)^n   ?=   8 log 3 - (16/3) log 2.
PARTIAL RESULT (a reduction, not a proof).  The WZ pair with ratios
   sigma(n,k) = -2(k-2n)(k+n+1) / (3(k-n+1)(2k+2n+3)),
   rho(n,k)   = (k-n)(6n+5)(6n+7)(k+n+1) / (4(2n+1)(k-2n-2)(k-2n-1)(2k+2n+3)),
   R(n,k)     = 3(k-n)((2n+1)k - 10n^2 - 17n - 7) / (4(2n+1)(k-2n-2)(k-2n-1))
(from Ebisu's 2F1(a, 1-2a; 3/2-2a; -1/3), a -> -n, k -> k-n) is finite on the triangle 0 <= k < n and has poles on
n <= k <= 2n.  Telescoping on the triangle (F(1,0) = 1, T_N = sum_{k<N} F(N,k)):
   sum_{n=1}^{N-1} G(n,0) = T_1 - T_N + sum_{n=1}^{N-1} d_n,   d_n := F(n+1,n) + lim_{k->n} G(n,k),
with G(n,0) = -(9/35) t(n), t the summand of S, T_N -> 0, and d_n = -(4/35) u_n,
   u_n = (12n^2+18n+7) (1/2)_n (1)_n (5/6)_n (7/6)_n / ((3/2)_n (2)_n (5/4)_n (7/4)_n) (3/4)^n.
Hence  S = 28/9 + (4/9) U,  U = sum_{n>=1} u_n,  and Guillera's value is EQUIVALENT to
   U = sum_{n>=1} (6 + 2/(2n+1) - 1/(n+1)) (5/6)_n (7/6)_n / ((5/4)_n (7/4)_n) (3/4)^n  ?=  18 log 3 - 12 log 2 - 7.
This script checks the WZ equation exactly and everything else numerically (40 digits).  U is not evaluated."""
import sympy as sp, mpmath as mp
n, k = sp.symbols('n k')
sigma = -2*(k - 2*n)*(k + n + 1)/(3*(k - n + 1)*(2*k + 2*n + 3))
rho = (k - n)*(6*n + 5)*(6*n + 7)*(k + n + 1)/(4*(2*n + 1)*(k - 2*n - 2)*(k - 2*n - 1)*(2*k + 2*n + 3))
R = 3*(k - n)*((2*n + 1)*k - 10*n**2 - 17*n - 7)/(4*(2*n + 1)*(k - 2*n - 2)*(k - 2*n - 1))
assert sp.cancel(rho - 1 - (R.subs(k, k + 1)*sigma - R)) == 0
assert sp.cancel(rho*sigma.subs(n, n + 1) - sigma*rho.subs(k, k + 1)) == 0
print("WZ and compatibility: exact OK")
tau = sp.cancel(rho.subs(k, n - 1)*sigma.subs({n: n + 1, k: n - 1}, simultaneous=True))
gt = sp.cancel(sp.limit(sp.cancel(R*sigma.subs(k, k - 1)), k, n))
D = sp.factor(sp.cancel(tau + gt))
assert sp.cancel(D + 2*n*(12*n**2 + 18*n + 7)/(3*(2*n + 1)*(4*n + 1)*(4*n + 3))) == 0
print("d_n / F(n,n-1) =", D)
mp.mp.dps = 45
sf = sp.lambdify((n, k), sigma, 'mpmath'); rf_ = sp.lambdify((n, k), rho, 'mpmath'); Rf = sp.lambdify((n, k), R, 'mpmath'); Df = sp.lambdify(n, D, 'mpmath')
def F(a, b):
    v = mp.mpf(1)
    for j in range(1, a): v *= rf_(mp.mpf(j), mp.mpf(0))
    for j in range(0, b): v *= sf(mp.mpf(a), mp.mpf(j))
    return v
t = lambda j: mp.rf(1, j)*mp.rf(mp.mpf(5)/6, j)*mp.rf(mp.mpf(7)/6, j)/mp.rf(mp.mpf(3)/2, j)**3*(7 + 10*j)*(mp.mpf(-9)/16)**j
u = lambda j: (12*j*j + 18*j + 7)*mp.rf(mp.mpf(1)/2, j)*mp.rf(1, j)*mp.rf(mp.mpf(5)/6, j)*mp.rf(mp.mpf(7)/6, j)/(mp.rf(mp.mpf(3)/2, j)*mp.rf(2, j)*mp.rf(mp.mpf(5)/4, j)*mp.rf(mp.mpf(7)/4, j))*(mp.mpf(3)/4)**j
d = lambda j: F(j, j - 1)*Df(mp.mpf(j))
for j in range(1, 8):
    assert abs(Rf(mp.mpf(j), mp.mpf(0))*F(j, 0) + mp.mpf(9)/35*t(j)) < mp.mpf(10)**-40
    assert abs(d(j) + mp.mpf(4)/35*u(j)) < mp.mpf(10)**-40
for N in (2, 5, 12, 30):
    lhs = mp.fsum(Rf(mp.mpf(j), mp.mpf(0))*F(j, 0) for j in range(1, N))
    rhs = 1 - mp.fsum(F(N, j) for j in range(N)) + mp.fsum(d(j) for j in range(1, N))
    assert abs(lhs - rhs) < mp.mpf(10)**-38, N
print("triangle identity for N = 2, 5, 12, 30; G(n,0) = -(9/35) t(n); d_n = -(4/35) u_n: OK (numerically)")
S = mp.nsum(t, [0, mp.inf]); U = mp.nsum(u, [1, mp.inf])
print("S - (28/9 + 4U/9) =", mp.nstr(S - (mp.mpf(28)/9 + 4*U/9), 5), "  S - (8 log3 - 16/3 log2) =", mp.nstr(S - (8*mp.log(3) - mp.mpf(16)/3*mp.log(2)), 5),
      "  U - (18 log3 - 12 log2 - 7) =", mp.nstr(U - (18*mp.log(3) - 12*mp.log(2) - 7), 5))
