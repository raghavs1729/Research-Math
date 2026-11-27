#!/usr/bin/env python3
"""
THEOREM 12 (J. Guillera's conjecture, "Ramanujan series with a shift", arXiv:1802.02023, Tables 1-2,
level 2, q = -exp(-pi sqrt 9), z = -1/48: the half-shift of Ramanujan's series
sum (28n+3)(1/2)_n(1/4)_n(3/4)_n/(n!^3 (-48)^n) = 16/(pi sqrt 3)):

   sum_{n>=0} (28n+17) (1)_n (3/4)_n (5/4)_n / ((3/2)_n^3 (-48)^n)  =  32 log(27/16)  =  96 log 3 - 128 log 2,

equivalently  sum_{k>=1} (28k-11) binom(4k,2k) / (k (4k-1) (-48)^k binom(2k,k)^2) = (4/3) log 2 - log 3.

Verification script for the proof in RESULTS.md.
  psi(n,k) = (2k+2n+3/2) 3^(1-n) (-1)^n n! Gamma(2n+3/2) Gamma(k+n+1) Gamma(3k+2n+3/2)
             / ( Gamma(k+n+3/2) Gamma(n+1/2) Gamma(3k+4n+4) )
  (the coefficient of b+n in Seed4(1/2, b, k+1/2+n) at b = -n; Au's Seed4 = Dougall's 5F4 with (b/3, b/3+1/3, b/3+2/3)),
  gamma = R psi,   WZ:  psi(n+1,k) - psi(n,k) = gamma(n,k+1) - gamma(n,k).
"""
import sympy as sp
from sympy import Rational as Q
n, k, w = sp.symbols('n k w')
N_, X_ = sp.symbols('N X', nonnegative=True)

sigma = (k+n+1)*(4*k+4*n+7)*(6*k+4*n+3)*(6*k+4*n+5)*(6*k+4*n+7)/(4*(2*k+2*n+3)*(3*k+4*n+4)*(3*k+4*n+5)*(3*k+4*n+6)*(4*k+4*n+3))
rho = -(n+1)*(4*n+3)*(4*n+5)*(k+n+1)*(4*k+4*n+7)*(6*k+4*n+3)*(6*k+4*n+5)/(
    12*(2*n+1)*(2*k+2*n+3)*(3*k+4*n+4)*(3*k+4*n+5)*(3*k+4*n+6)*(3*k+4*n+7)*(4*k+4*n+3))
P = 36*k**3 + (156*n+120)*k**2 + (232*n**2+362*n+138)*k + (112*n**3+264*n**2+203*n+51)
cR = sp.symbols('cR')
R0 = P/((2*n+1)*(3*k+4*n+4)*(4*k+4*n+3))
num = sp.numer(sp.together(rho - 1 - cR*(R0.subs(k, k+1)*sigma - R0)))
sol = sp.solve(sp.Poly(sp.expand(num), n, k).coeffs(), cR)
print("constant of the certificate:", sol)
R = sol[cR]*R0
assert sp.cancel(rho - 1 - (R.subs(k, k+1)*sigma - R)) == 0
assert sp.cancel(rho*sigma.subs(n, n+1) - sigma*rho.subs(k, k+1)) == 0
print("(WZ) and compatibility: identically 0: OK;  R =", sp.factor(R))

# ---- closed form has these ratios (numerically, 40 digits)
import mpmath as mp
mp.mp.dps = 50
def psi(nn, kk):
    nn = mp.mpf(nn); kk = mp.mpf(kk)
    return ((2*kk+2*nn+mp.mpf(3)/2)*mp.mpf(3)**(1-nn)*(-1)**int(nn)*mp.factorial(nn)*mp.gamma(2*nn+mp.mpf(3)/2)*mp.gamma(kk+nn+1)
            * mp.gamma(3*kk+2*nn+mp.mpf(3)/2)/(mp.gamma(kk+nn+mp.mpf(3)/2)*mp.gamma(nn+mp.mpf(1)/2)*mp.gamma(3*kk+4*nn+4)))
_sf = sp.lambdify((n, k), sigma, 'mpmath'); _rf = sp.lambdify((n, k), rho, 'mpmath'); _Rf = sp.lambdify((n, k), R, 'mpmath')
sf = lambda a, b: _sf(mp.mpf(a), mp.mpf(b)); rf_ = lambda a, b: _rf(mp.mpf(a), mp.mpf(b)); Rf = lambda a, b: _Rf(mp.mpf(a), mp.mpf(b))
for (a, b) in [(0, 0), (0, 3), (2, 1), (5, 7), (3, 0)]:
    assert abs(psi(a, b+1)/psi(a, b) - sf(a, b)) < mp.mpf(10)**-40
    assert abs(psi(a+1, b)/psi(a, b) - rf_(a, b)) < mp.mpf(10)**-40
print("closed form of psi has the ratios sigma, rho: OK")

# ---- boundary data
t = lambda j: (28*j+17)*mp.rf(1, j)*mp.rf(mp.mpf(3)/4, j)*mp.rf(mp.mpf(5)/4, j)/(mp.rf(mp.mpf(3)/2, j)**3*mp.mpf(-48)**j)
assert sp.cancel(R.subs(k, 0) - (28*n+17)/(12*(2*n+1))) == 0
for j in range(6):
    assert abs(Rf(j, 0)*psi(j, 0) - t(j)/32) < mp.mpf(10)**-40
tr = (n+1)*(n+Q(3, 4))*(n+Q(5, 4))/((n+Q(3, 2))**3*(-48))*(28*n+45)/(28*n+17)
assert sp.cancel(R.subs({k: 0, n: n+1})/R.subs(k, 0)*rho.subs(k, 0) - tr) == 0
print("gamma(n,0) = t(n)/32 (n = 0 and equal ratios, exactly): OK")
for j in range(6):
    assert abs(psi(0, j) - mp.mpf(3)/4*(4*j+3)*mp.factorial(j)*mp.gamma(3*j+mp.mpf(3)/2)/(mp.gamma(j+mp.mpf(3)/2)*mp.factorial(3*j+3))) < mp.mpf(10)**-40
print("psi(0,k) = (3/4)(4k+3) k! Gamma(3k+3/2)/(Gamma(k+3/2)(3k+3)!): OK")

# ---- inequalities for the limit step (coefficientwise; n = N >= 0, k = X >= 0)
def nonneg(E):
    A, B = sp.fraction(sp.cancel(sp.together(E.subs({n: N_, k: X_}))))
    ca, cb = sp.Poly(A, N_, X_).coeffs(), sp.Poly(B, N_, X_).coeffs()
    return (all(c >= 0 for c in ca) and all(c > 0 for c in cb)) or (all(c <= 0 for c in ca) and all(c < 0 for c in cb))
assert nonneg(sigma)
assert nonneg((k + 3*n + 3)/(k + 3*n + 5) - sigma), "sigma <= (k+3n+3)/(k+3n+5)"
assert nonneg(-rho.subs(k, 0)) and nonneg(Q(1, 48) + rho.subs(k, 0)), "|rho(n,0)| <= 1/48"
assert nonneg(R) and nonneg(k + 2 - R), "0 <= R <= k+2"
print("0 < sigma <= (k+3n+3)/(k+3n+5), -1/48 <= rho(n,0) < 0, 0 <= R <= k+2: OK")

# ---- the boundary sum B = sum_k psi(0,k) = log(27/16)
# (b) generating function (trisection of the binomial series), checked as a power series to O(w^40)
# The coefficient of w^j in (1 + omega w) sqrt(1 - omega w) is omega^j c_j, c_j the coefficient in (1+w) sqrt(1-w);
# summing over the three cube roots of unity leaves 3 c_j for 3 | j and 0 otherwise.  So the identity is
#   c_0 = 1   and   (4k+3) Gamma(3k+3/2)/(3k+3)! = -(4 sqrt(pi)/3) c_{3k+3}   (k >= 0):
om = sp.exp(2*sp.pi*sp.I/3)
cser = sp.Poly(sp.series((1 + w)*sp.sqrt(1 - w), w, 0, 124).removeO(), w)
assert cser.coeff_monomial(1) == 1
for kk in range(0, 40):
    lhs = (4*kk + 3)*sp.gamma(3*kk + Q(3, 2))/sp.factorial(3*kk + 3)
    assert sp.simplify(lhs + sp.Rational(4, 3)*sp.sqrt(sp.pi)*cser.coeff_monomial(w**(3*kk + 3))) == 0
print("E(w) := sum (4k+3) Gamma(3k+3/2)/(3k+3)! w^(3k+3) = (4 sqrt(pi)/3) [1 - (1/3) sum_omega (1+omega w) sqrt(1-omega w)]  (coefficients k < 40 exactly): OK")
# (e) antiderivatives: d/dw log[(2 + om w + 2 s)(om + 2 conj(om) w - 2 om s)] = (1 - (1 + om w)/s)/w,  s^2 = 1 + om w + conj(om) w^2
s = sp.symbols('s')
for o in (sp.Integer(1), om, sp.conjugate(om)):
    oc = sp.conjugate(o)
    q = 1 + o*w + oc*w**2
    ds = (o + 2*oc*w)/(2*s)                     # s' from s^2 = q
    A = (2 + o*w + 2*s)*(o + 2*oc*w - 2*o*s)
    dA = sp.diff(A, w) + sp.diff(A, s)*ds
    expr = sp.together(dA/A - (1 - (1 + o*w)/s)/w)
    numer = sp.expand(sp.numer(expr))
    numer = sp.expand(sp.simplify(sp.expand(sp.rem(sp.Poly(numer, s), sp.Poly(s**2 - q, s)).as_expr())))
    assert sp.simplify(numer) == 0, (o, numer)
    A1 = (2 + o)*(o + 2*oc) if o != 1 else (3 + 2*sp.sqrt(3))*(3 - 2*sp.sqrt(3))
    assert sp.simplify(sp.Abs(sp.expand(A1)) - 3) == 0
    assert sp.simplify(sp.Abs(4*(o - 2*o)) - 4) == 0
print("antiderivative Lambda_omega = log[(2+omega w+2s)(omega+2 conj(omega) w-2 omega s)], |Lambda(1)| = 3, |Lambda(0)| = 4: OK")
# numerics of every step
mp.mp.dps = 45
B = mp.nsum(lambda j: psi(0, j), [0, mp.inf])
omn = mp.exp(2j*mp.pi/3)
J = [mp.quad(lambda x, o=o: (1 - (1+o*x)/mp.sqrt((1 + o*x + mp.conj(o)*x*x)))/x, [0, 0.5, 0.9, 1]) for o in (1, omn, mp.conj(omn))]
assert all(abs(j.real - mp.log(mp.mpf(3)/4)) < mp.mpf(10)**-18 for j in J)   # endpoint singularity limits the quadrature
I0 = mp.quad(lambda x: 3*((1-x**3)**-0.5 - 1)/x, [0, 0.5, 0.9, 0.99, 1])
assert abs(I0 - 2*mp.log(2)) < mp.mpf(10)**-18
assert abs(B - mp.log(mp.mpf(27)/16)) < mp.mpf(10)**-38
print("B = sum_k psi(0,k) = 2 log 2 + 3 log(3/4) = log(27/16): |B - log(27/16)| =", mp.nstr(abs(B - mp.log(mp.mpf(27)/16)), 5))
# telescoped identity at finite N
for NN in (1, 2, 4):
    lhs_ = mp.fsum(Rf(j, 0)*psi(j, 0) for j in range(NN))
    rhs_ = B - mp.nsum(lambda j: psi(NN, j), [0, mp.inf])
    assert abs(lhs_ - rhs_) < mp.mpf(10)**-35
print("sum_{n<N} gamma(n,0) = B - sum_k psi(N,k) for N = 1, 2, 4 (35 digits): OK")
# ---- final identities, 110 digits
mp.mp.dps = 130
S = mp.nsum(t, [0, mp.inf])
val = 32*mp.log(mp.mpf(27)/16)
print("final: |sum (28n+17)(1)_n(3/4)_n(5/4)_n/((3/2)_n^3 (-48)^n) - 32 log(27/16)| =", mp.nstr(abs(S - val), 5))
assert abs(S - val) < mp.mpf(10)**-110
S2 = mp.nsum(lambda j: (28*j-11)*mp.binomial(4*j, 2*j)/(j*(4*j-1)*mp.mpf(-48)**j*mp.binomial(2*j, j)**2), [1, mp.inf])
assert abs(S2 - (mp.mpf(4)/3*mp.log(2) - mp.log(3))) < mp.mpf(10)**-110
print("binomial form: sum_{k>=1} (28k-11) C(4k,2k)/(k(4k-1)(-48)^k C(2k,k)^2) = (4/3) log 2 - log 3: OK (110 digits)")
print("ALL CHECKS PASSED")
