#!/usr/bin/env python3
"""Independent checks for Theorem 3 (steps 2-5) and Theorem 4."""
import sympy as sp, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 40

print("=== Theorem 3")
# --- Psi_N closed form: exact symbolic derivative vs closed form, at exact points, N = 1..4
y, s = sp.symbols('y s', positive=True)
for N in (1, 2, 3, 4):
    direct = y**(2*N)*sp.diff(sp.log(y)/(1 - y), y, 2*N)
    # finite form from Leibniz: -(2N)!/(1-y)^{2N+1} y^{2N} [ log(1+u) - sum_{i<=2N} (-1)^{i-1} u^i/i ],  u = (1-y)/y
    u = (1 - y)/y
    leib = -sp.factorial(2*N)*y**(2*N)/(1 - y)**(2*N + 1)*(-sp.log(y) - sum((-1)**(i - 1)*u**i/i for i in range(1, 2*N + 1)))
    assert sp.simplify(direct - leib) == 0
    for yy in (sp.Rational(1, 20), sp.Rational(3, 10), sp.Rational(3, 4), sp.Rational(99, 100)):
        d = sp.N(direct.subs(y, yy), 40)
        ym = mp.mpf(yy.p)/yy.q
        integ = mp.quad(lambda t: t**(2*N)/(1 + t), mp.linspace(0, (1 - ym)/ym, 20))
        cf = -mp.factorial(2*N)*ym**(2*N)/(1 - ym)**(2*N + 1)*integ
        assert abs(mp.mpf(d) - cf) < mp.mpf(10)**-25*abs(cf), (N, yy, d, cf)
        assert abs(cf) <= mp.factorial(2*N)/((2*N + 1)*ym)
print("Psi_N: Leibniz form == y^{2N} d^{2N}/dy^{2N}[log y/(1-y)] symbolically (N<=4); integral closed form and bound (2N)!/((2N+1)y) at y = 1/20,3/10,3/4,99/100: OK")

def c(N): return (mp.mpf(1)/9)**N*mp.gamma(3*N + mp.mpf(3)/2)*mp.factorial(N)**2/(4*mp.factorial(2*N)*mp.factorial(3*N)*mp.gamma(N + mp.mpf(5)/6)*mp.gamma(N + mp.mpf(7)/6))
def F(N, X): return c(N)*(mp.mpf(3)/4)**X*mp.rf(-X, 2*N)*mp.gamma(X + N + 1)/mp.gamma(X + N + mp.mpf(3)/2)
def R(N, X): return (8*N**2 + 6*N - (2*N + 3)*X - X**2)/mp.mpf((3*N + 1)*(3*N + 2))
def G(N, X): return R(N, X)*F(N, X)
# exact derivative of F at integers via log-derivative pieces (avoid numerical differentiation):
def dF(N, k):
    # F = c (3/4)^x P(x) g(x), P = (-x)_{2N}, g = Gamma(x+N+1)/Gamma(x+N+3/2)
    g = mp.gamma(k + N + 1)/mp.gamma(k + N + mp.mpf(3)/2)
    P = mp.rf(-k, 2*N)
    dP = mp.mpf(0)
    for i in range(2*N):
        term = mp.mpf(-1)
        for j in range(2*N):
            if j != i: term *= (-k + j)
        dP += term
    lg = mp.log(mp.mpf(3)/4) + mp.digamma(k + N + 1) - mp.digamma(k + N + mp.mpf(3)/2)
    return c(N)*(mp.mpf(3)/4)**k*g*(dP + P*lg)
def D(N, kmax=700): return mp.fsum(dF(N, k) for k in range(kmax))
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
D0 = D(0)
print("D_0 + 3K/(2pi) =", mp.nstr(D0 + 3*K/(2*mp.pi), 5))
# growth claim of step 2: |F'(n,x)| / (x^{2n+3} (3/4)^x) bounded
for N in (0, 1, 3):
    print("  step 2: sup_{k<=2000} |F'(%d,k)|/(k^(2n+3)(3/4)^k) = %s; same for G' (num. diff): %s" % (N,
          mp.nstr(max(abs(dF(N, k))/(mp.mpf(k)**(2*N + 3)*(mp.mpf(3)/4)**k) for k in range(1, 2001, 7)), 5),
          mp.nstr(max(abs(mp.diff(lambda X: G(N, X), k))/(mp.mpf(k)**(2*N + 3)*(mp.mpf(3)/4)**k) for k in range(1, 400, 13)), 5)))
# G'(n,0) = -(3/(4pi)) t(n+1): exact formula of the text
t = lambda m: mp.mpf((4*m - 1)*3**m)/((2*m - 1)*m**2*mp.binomial(2*m, m)*mp.binomial(3*m, m))
for N in range(1, 12):
    gp = -R(N, 0)*c(N)*mp.factorial(2*N - 1)*mp.factorial(N)/mp.gamma(N + mp.mpf(3)/2)
    assert abs(gp + 3/(4*mp.pi)*t(N + 1)) < mp.mpf(10)**-35*abs(gp)
assert abs(-mp.mpf(3)/2*c(0)/mp.gamma(mp.mpf(3)/2) + 3/(4*mp.pi)*t(1)) < mp.mpf(10)**-38
print("step 3: G'(n,0) = -(3/(4 pi)) t(n+1), n = 0..11 (closed form of the text): OK   [n=1 base case by hand: G'(1,0) = -7/(160 pi), t(2) = 7/120]")
# telescoped: sum_{n<N} G'(n,0) = D_0 - D_N, and size of D_N
part = mp.mpf(0)
for N in range(1, 9):
    part += -3/(4*mp.pi)*t(N)
    DN = D(N, 900)
    bound = c(N)*mp.factorial(2*N)*4*mp.beta(N, mp.mpf(1)/2)/(3*mp.sqrt(mp.pi)*(2*N + 1))
    print("  N=%d: D_N=%s  (D_0 - sum_{n<N}G'(n,0)) - D_N = %s  bound=%s  bound*9^N*N=%s  |D_N|*9^N*N=%s" %
          (N, mp.nstr(DN, 8), mp.nstr(D0 - part - DN, 3), mp.nstr(bound, 6), mp.nstr(bound*9**N*N, 6), mp.nstr(abs(DN)*9**N*N, 6)))
for N in (50, 500, 5000):
    lb = (-N*mp.log(9) + mp.loggamma(3*N + 1.5) + 2*mp.loggamma(N + 1) - mp.log(4) - mp.loggamma(3*N + 1) - mp.loggamma(N + mp.mpf(5)/6) - mp.loggamma(N + mp.mpf(7)/6)
          + mp.log(4) + mp.loggamma(N) + mp.loggamma(0.5) - mp.loggamma(N + 0.5) - mp.log(3*mp.sqrt(mp.pi)*(2*N + 1)))
    print("  N=%d: bound*9^N*N = %s   (sqrt3/6 = %s)" % (N, mp.nstr(mp.exp(lb + N*mp.log(9))*N, 8), mp.nstr(mp.sqrt(3)/6, 8)))
# log-cosine evaluation
cl2 = lambda th: mp.clsin(2, th)
I = mp.quad(lambda tt: mp.log(3*tt/4)/((1 - 3*tt/4)*mp.sqrt(1 - tt)), mp.linspace(0, 1, 30))
I2 = 8/mp.sqrt(3)*mp.quad(lambda th: mp.log(mp.cos(th + mp.pi/6)) + mp.log(mp.cos(th - mp.pi/6)) - 2*mp.log(mp.cos(th)), mp.linspace(0, mp.pi/3, 30))
print("I (t-form) + 4K = %s ; I (theta-form) + 4K = %s" % (mp.nstr(I + 4*K, 3), mp.nstr(I2 + 4*K, 3)))
print("Cl2(pi/3) - 1.5 Cl2(2pi/3) = %s ; Cl2(2pi/3) - (sqrt3/2)K = %s ; 0.5Cl2(2pi/3) - Cl2(pi/3) + Cl2(2pi/3) = %s" %
      (mp.nstr(cl2(mp.pi/3) - 1.5*cl2(2*mp.pi/3), 3), mp.nstr(cl2(2*mp.pi/3) - mp.sqrt(3)/2*K, 3), mp.nstr(0.5*cl2(2*mp.pi/3) - cl2(mp.pi/3) + cl2(2*mp.pi/3), 3)))
for al in (0.2, 0.7, 1.3):
    assert abs(mp.quad(lambda p: mp.log(2*mp.cos(p)), [0, al]) - cl2(mp.pi - 2*al)/2) < mp.mpf(10)**-30
print("int_0^alpha log(2cos) = Cl2(pi-2alpha)/2: OK")
for (a_, b_, v) in ((mp.pi/6, mp.pi/2, -cl2(2*mp.pi/3)/2), (-mp.pi/6, mp.pi/6, cl2(2*mp.pi/3)), (0, mp.pi/3, cl2(mp.pi/3)/2)):
    assert abs(mp.quad(lambda p: mp.log(2*mp.cos(p)), mp.linspace(a_, b_, 12)) - v) < mp.mpf(10)**-25
print("three pieces of the log-cos integral: OK")

print("\n=== Theorem 4")
def sig(n, k): return Fr((4*k + 6*n + 1)*(4*k + 6*n + 3), 16*(k + n + 1)*(k + 3*n + 1))
def rho(n, k): return Fr(3*(6*n + 1)*(6*n + 5)*(4*k + 6*n + 1)*(4*k + 6*n + 3)*(4*k + 6*n + 5), 1024*(n + 1)*(k + n + 1)*(k + 3*n + 1)*(k + 3*n + 2)*(k + 3*n + 3))
def Q(n, k): return Fr(64*k**3 + (160*n + 144)*k**2 + (240*n**2 + 432*n + 92)*k + 360*n**3 + 612*n**2 + 230*n + 15, 64*(k + 3*n + 1)*(k + 3*n + 2))
# exact: 1 - sigma - (3/4)(n+1)/(k+3n+1) = (4nk+4k+16n+1)/(16(k+n+1)(k+3n+1))
for n in range(0, 30):
    for k in range(0, 60):
        assert 1 - sig(n, k) - Fr(3*(n + 1), 4*(k + 3*n + 1)) == Fr(4*n*k + 4*k + 16*n + 1, 16*(k + n + 1)*(k + 3*n + 1))
        assert Q(n, k) <= k + 3*n + 3
        assert (n + 1)*rho(n, k) - n == Q(n, k + 1)*sig(n, k) - Q(n, k)
        assert rho(n, k)*sig(n + 1, k) == sig(n, k)*rho(n, k + 1)
print("1-sigma-(3/4)(n+1)/(k+3n+1) = (4nk+4k+16n+1)/(16(k+n+1)(k+3n+1)) > 0; Q <= k+3n+3; WZ; compatibility (grid 30x60, exact): OK")
# bound h(n,k) <= h(n,0) ((3n+1)/(k+3n+1))^{3(n+1)/4}
worst = 0
for n in range(0, 12):
    p = mp.mpf(1)
    for k in range(0, 4000):
        b = (mp.mpf(3*n + 1)/(k + 3*n + 1))**(mp.mpf(3*(n + 1))/4)
        worst = max(worst, p/b)
        p *= mp.mpf(sig(n, k).numerator)/sig(n, k).denominator
print("max h(n,k)/[h(n,0)((3n+1)/(k+3n+1))^(3(n+1)/4)] over n<12, k<4000 =", mp.nstr(worst, 8), "(must be <= 1)")
# sum_k h(n,k)/h(n,0) = 3F2(3n/2+1/4, 3n/2+3/4, 1; n+1, 3n+1; 1); bound 1 + 4(3n+1)/(3n-1) <= 33/5 for n >= 2
for n in (1, 2, 3, 5, 10, 40):
    S = mp.hyp3f2(mp.mpf(3*n)/2 + 0.25, mp.mpf(3*n)/2 + 0.75, 1, n + 1, 3*n + 1, 1)
    bnd = mp.nsum(lambda kk: (mp.mpf(3*n + 1)/(kk + 3*n + 1))**(mp.mpf(3*(n + 1))/4), [0, mp.inf])
    print("  n=%2d: sum_k h(n,k)/h(n,0) = %s ; sum of the bounding series = %s ; integral bound 1+4(3n+1)/(3n-1) = %s" % (n, mp.nstr(S, 8), mp.nstr(bnd, 8), mp.nstr(1 + mp.mpf(4*(3*n + 1))/(3*n - 1), 6)))
# limit K (1/4)_K (3/4)_K / K!^2
for KK in (10**3, 10**6):
    v = KK*mp.exp(mp.loggamma(KK + 0.25) + mp.loggamma(KK + 0.75) - 2*mp.loggamma(KK + 1))/(mp.gamma(0.25)*mp.gamma(0.75))
    print("  K=%d: K (1/4)_K (3/4)_K/K!^2 = %s ; 1/(Gamma(1/4)Gamma(3/4)) = %s = 1/(sqrt2 pi) = %s" % (KK, mp.nstr(v, 12), mp.nstr(1/(mp.gamma(0.25)*mp.gamma(0.75)), 12), mp.nstr(1/(mp.sqrt(2)*mp.pi), 12)))
# telescoped identity for N = 1..5: sum_{n<N} Q(n,0)h(n,0) + N sum_k h(N,k) = 1/(sqrt2 pi)
hn0 = [mp.mpf(1)]
for n in range(0, 6): hn0.append(hn0[-1]*mp.mpf(rho(n, 0).numerator)/rho(n, 0).denominator)
part = mp.mpf(0)
for N in range(1, 6):
    part += mp.mpf(Q(N - 1, 0).numerator)/Q(N - 1, 0).denominator*hn0[N - 1]
    S = hn0[N]*mp.hyp3f2(mp.mpf(3*N)/2 + 0.25, mp.mpf(3*N)/2 + 0.75, 1, N + 1, 3*N + 1, 1)
    print("  N=%d: sum_{n<N} Q(n,0)h(n,0) + N sum_k h(N,k) - 1/(sqrt2 pi) = %s ;  N sum_k h(N,k) = %s <= 7N(27/32)^N = %s" %
          (N, mp.nstr(part + N*S - 1/(mp.sqrt(2)*mp.pi), 3), mp.nstr(N*S, 6), mp.nstr(7*N*(mp.mpf(27)/32)**N, 6)))
