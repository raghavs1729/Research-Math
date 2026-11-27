#!/usr/bin/env python3
"""Theorem 11 (B2)-(B5): symbolic checks + numerical O-constants, N up to 10^5."""
import sympy as sp, mpmath as mp, math, sys
mp.mp.dps = 40

# ---------- symbolic
x = sp.symbols('x', positive=True)
L = sp.log((1 - x)**2*(2*x + 1)**4/((1 + x)**2*(2*x - 1)**4))
xs = sp.sqrt(3)/2
print("L(x*)   =", sp.simplify(sp.expand_log(L.subs(x, xs), force=True)), " numeric:", sp.N(L.subs(x, xs), 30))
print("L'(x)   =", sp.simplify(sp.diff(L, x)), ";  L'(x*) =", sp.simplify(sp.diff(L, x).subs(x, xs)))
L2 = sp.lambdify(x, sp.diff(L, x, 2))
print("sup |L''| on [0.8,0.93] ~", max(abs(L2(0.8 + 0.13*i/1000)) for i in range(1001)))
phi = (-6*sp.log(2) + 2*(2*x + 1)*sp.log(2*x + 1) - 2*(1 - x)*sp.log(1 - x) - 2*(1 + x)*sp.log(1 + x) - 2*(2*x - 1)*sp.log(2*x - 1))
print("phi' - L =", sp.simplify(sp.expand_log(sp.diff(phi, x) - L, force=True)))
# phi(x*) = 0 exactly
e = sp.expand_log(phi.subs(x, xs), force=True)
print("phi(x*) numeric:", sp.N(e, 50), "; exp(phi(x*)) simplified:", sp.nsimplify(sp.simplify(sp.exp(phi.subs(x, xs)))))
# exact: exp(phi(x*)) = 2^-6 (1+s)^{2(1+s)} / ( ((2-s)/2)^{2-s} ((2+s)/2)^{2+s} (s-1)^{2(s-1)} ), s = sqrt3.
# with (1+s)^2 = 2(2+s), (s-1)^2 = 2(2-s):  = 2^-6 2^{1+s}(2+s)^{1+s} 2^{2-s} 2^{2+s} / ( (2-s)^{2-s} (2+s)^{2+s} 2^{s-1}(2-s)^{s-1} )
#   = 2^{-6+1+s+2-s+2+s-s+1} (2+s)^{-1} (2-s)^{-1} = 2^0 * 1 = 1.
psi0 = 2**sp.Rational(5, 2)*sp.pi**sp.Rational(-5, 2)*x*(1 - x)/((1 + x)*(2*x - 1)**4)
print("psi0(x*) - sqrt3/(sqrt2 pi^(5/2)) =", sp.simplify(psi0.subs(x, xs) - sp.sqrt(3)/(sp.sqrt(2)*sp.pi**sp.Rational(5, 2))),
      "; psi0(x*) sqrt(pi/12) - sqrt2/(4pi^2) =", sp.simplify(psi0.subs(x, xs)*sp.sqrt(sp.pi/12) - sp.sqrt(2)/(4*sp.pi**2)))

# ---------- symbolic Stirling: derive phi, -1/2 log N, psi0 from scratch
Ns, xx = sp.symbols('N x', positive=True)
def st(z, a):   # log Gamma(z + a) ~ (z + a - 1/2) log z - z + log(2pi)/2   (error O(1/z))
    return (z + a - sp.Rational(1, 2))*sp.log(z) - z + sp.log(2*sp.pi)/2
logC = ((6 - 8*Ns)*sp.log(2) - sp.log(sp.sqrt(2)*sp.pi) + 2*st(Ns, sp.Rational(3, 4)) + st(Ns, sp.Rational(1, 4)) - st(Ns, -sp.Rational(1, 4)) - st(2*Ns, -1))
loga = (logC + sp.log(2*xx*Ns) + st((4*xx + 2)*Ns, sp.Rational(1, 2)) - 2*st((1 - xx)*Ns, 0) - 2*st((1 + xx)*Ns, 1) - st((4*xx - 2)*Ns, sp.Rational(9, 2)))
claim = Ns*phi.subs(x, xx) - sp.log(Ns)/2 + sp.log(psi0.subs(x, xx))
d = sp.simplify(sp.expand(sp.expand_log(loga - claim, force=True)))
print("Stirling main terms minus claimed [N phi - (1/2) log N + log psi0] =", d)

# ---------- numerics
lg = mp.loggamma
def la(N, k):
    N = mp.mpf(N); k = mp.mpf(k)
    logC = (6 - 8*N)*mp.log(2) - mp.log(mp.sqrt(2)*mp.pi) + 2*lg(N + mp.mpf(3)/4) + lg(N + mp.mpf(1)/4) - lg(N - mp.mpf(1)/4) - lg(2*N - 1)
    return logC + mp.log(2*k + 1) + lg(4*k + 2*N + mp.mpf(1)/2) - 2*lg(N - k) - 2*lg(N + k + 1) - lg(4*k - 2*N + mp.mpf(9)/2)   # arg > 0 for k > N/2 - 9/8
def logr(N, k):   # float, exact formula for sigma(N,k); valid (positive factors) for k >= N/2 - 1, k <= N-2
    s = math.log1p(2.0/(2*k + 1)) + 2*math.log((N - k - 1)/(N + k + 1))
    z = 4*k - 2*N + 4.5
    for i in range(4): s += math.log1p((4*N - 4)/(z + i))
    return s
Lf = lambda t: math.log((1 - t)**2*(2*t + 1)**4/((1 + t)**2*(2*t - 1)**4))
phif = lambda t: (-6*mp.log(2) + 2*(2*t + 1)*mp.log(2*t + 1) - 2*(1 - t)*mp.log(1 - t) - 2*(1 + t)*mp.log(1 + t) - 2*(2*t - 1)*mp.log(2*t - 1))
psif = lambda t: mp.mpf(2)**2.5*mp.pi**-2.5*t*(1 - t)/((1 + t)*(2*t - 1)**4)
target = mp.sqrt(2)/(4*mp.pi**2)

print("\n(B2) uniform:  N*|log r(k) - L(k/N)| over k/N in [0.8,0.93];  ratio-expansion constant; summed-expansion constant")
for N in (100, 1000, 10**4, 10**5, 10**6):
    c1 = max(N*abs(logr(N, k) - Lf(k/N)) for k in range(math.ceil(0.8*N), math.floor(0.93*N) + 1))
    k0 = math.floor(math.sqrt(3)/2*N); J = math.floor(N**0.6)
    js = [j for j in range(-J - 1, J + 2) if 0.8*N <= k0 + j <= min(0.93*N, N - 2)]
    c2 = max(abs(logr(N, k0 + j) + 24*j/N)/(1/N + j*j/N**2) for j in js)
    # summed: log a(k0+j) - log a(k0) + 12 j^2/N
    cum = 0.0; c3 = 0.0; e3 = 0.0
    for j in range(1, max(js) + 1):
        cum += logr(N, k0 + j - 1); err = abs(cum + 12*j*j/N); c3 = max(c3, err/(j/N + j**3/N**2)); e3 = max(e3, err)
    cum = 0.0
    for j in range(-1, min(js) - 1, -1):
        cum -= logr(N, k0 + j); err = abs(cum + 12*j*j/N); c3 = max(c3, err/(abs(j)/N + abs(j)**3/N**2)); e3 = max(e3, err)
    print("  N=%8d  J=%6d  |j| range [%d,%d]  C(B2 uniform)=%.3f  C(ratio exp)=%.3f  C(summed)=%.3f  max err*N^(1/5)=%.3f"
          % (N, J, min(js), max(js), c1, c2, c3, e3*N**0.2))

print("\n(B4) tail-ratio claims and a(min I)")
def b4(N):
    k0 = math.floor(math.sqrt(3)/2*N); J = math.floor(N**0.6)
    if k0 + J > N - 2: return None
    up = logr(N, k0 + J) <= -12*J/N
    dn = logr(N, k0 - J - 1) >= 12*J/N
    return up, dn
firstok = None; lastbad = None
for N in range(4, 20001):
    r_ = b4(N)
    if r_ is None or not (r_[0] and r_[1]): lastbad = N
print("  r(k0+J) <= e^{-12J/N} and r(k0-J-1) >= e^{12J/N}: last N <= 20000 where it fails or k0+J > N-2:", lastbad)
for N in (10**5, 10**6, 10**7): print("  N=%d:" % N, b4(N))
lastbad = None
for N in list(range(4, 3001)) + [10**4, 10**5]:
    minI = max(0, -((-(N - 2))//2)); k0 = math.floor(math.sqrt(3)/2*N)
    dlt = la(N, minI) - la(N, k0)
    if not dlt <= -4*mp.mpf(N)**0.6: lastbad = N
    if N in (4, 10, 100, 1000, 10**4, 10**5): print("  N=%6d: log(a(min I)/a(k0)) = %s,  -4 N^(3/5) = %s,  N log(4/27) = %s" % (N, mp.nstr(dlt, 8), mp.nstr(-4*mp.mpf(N)**0.6, 8), mp.nstr(N*mp.log(mp.mpf(4)/27), 8)))
print("  a(min I) <= a(k0) e^{-4 N^(3/5)}: last failing N in 4..3000:", lastbad)

print("\n(B5) Stirling: N*(log a(N,k) - [N phi(x) - (1/2) log N + log psi0(x)]),  x = k/N")
for N in (100, 1000, 10**4, 10**5):
    row = []
    for xt in (0.8, 0.83, mp.sqrt(3)/2, 0.9, 0.93):
        k = int(mp.floor(xt*N)); t = mp.mpf(k)/N
        row.append(mp.nstr(N*(la(N, k) - (N*phif(t) - mp.log(N)/2 + mp.log(psif(t)))), 6))
    print("  N=%6d  x~0.8,0.83,x*,0.9,0.93:" % N, row)

print("\n(B3)-(B5) combined: T_N, centre, tails")
for N in (1000, 10**4, 10**5):
    vals = [mp.exp(la(N, k)) for k in range(max(0, N//2), N)]   # k >= N/2: argument of Gamma positive
    base = max(0, N//2)
    # contribution of k < N/2 (use log|Gamma|; sign positive by the parity argument)
    def la_low(k):
        Nm = mp.mpf(N); km = mp.mpf(k)
        logC = (6 - 8*Nm)*mp.log(2) - mp.log(mp.sqrt(2)*mp.pi) + 2*lg(Nm + mp.mpf(3)/4) + lg(Nm + mp.mpf(1)/4) - lg(Nm - mp.mpf(1)/4) - lg(2*Nm - 1)
        return logC + mp.log(2*km + 1) + lg(4*km + 2*Nm + mp.mpf(1)/2) - 2*lg(Nm - km) - 2*lg(Nm + km + 1) - mp.re(lg(4*km - 2*Nm + mp.mpf(9)/2))
    low = mp.fsum(mp.exp(la_low(k)) for k in range(0, base)) if N <= 10**4 else mp.mpf(0)
    T = mp.fsum(vals) + low
    k0 = int(mp.floor(mp.sqrt(3)/2*N)); J = int(mp.floor(mp.mpf(N)**0.6)); ak0 = vals[k0 - base]
    centre = mp.fsum(vals[k0 - J - base:k0 + J + 1 - base])
    up = mp.fsum(vals[k0 + J + 1 - base:]); dn = mp.fsum(vals[:k0 - J - base]) + low
    print("  N=%6d: T_N/target-1=%s  a(k0)sqrtN/psi0(x*)-1=%s  centre/(a(k0)sqrt(pi N/12))-1=%s  upper tail/(a(k0)sqrtN)=%s  lower=%s  (bound N^(2/5)e^(-11N^(1/5))=%s)"
          % (N, mp.nstr(T/target - 1, 3), mp.nstr(ak0*mp.sqrt(N)/psif(mp.sqrt(3)/2) - 1, 3), mp.nstr(centre/(ak0*mp.sqrt(mp.pi*N/12)) - 1, 3),
             mp.nstr(up/(ak0*mp.sqrt(N)), 3), mp.nstr(dn/(ak0*mp.sqrt(N)), 3), mp.nstr(mp.mpf(N)**0.4*mp.exp(-11*mp.mpf(N)**0.2), 3)))
