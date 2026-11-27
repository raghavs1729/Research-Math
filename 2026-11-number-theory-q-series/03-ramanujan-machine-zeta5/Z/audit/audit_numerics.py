# Independent audit script (does not reuse verify/z_checks.py).
from mpmath import mp, mpf, pi, log, gamma, loggamma, psi, exp, nsum, inf, quad, ellipk, gammainc, sqrt, harmonic, hyper, nstr, polyroots
mp.dps = 40
N = 400
# q-expansion of f = q prod (1-q^{2m})^4 (1-q^{4m})^4, coefficients c[1..N]
P = [0]*(N+1); P[0] = 1
for m in range(1, N+1):
    for step in (2*m, 4*m):
        if step > N: continue
        for _ in range(4):
            for i in range(N, step-1, -1):
                P[i] -= P[i-step]
c = [0]*(N+2)
for i in range(N): c[i+1] = P[i]
print("c_1..c_12:", c[1:13])
# Hecke sanity: multiplicativity c(3)c(5)=c(15), c(p^2)=c(p)^2-p^3 c(1)
assert c[3]*c[5] == c[15] and c[9] == c[3]**2 - 27 and c[25] == c[5]**2 - 125
def fit(t):
    q = exp(-2*pi*t); return sum(c[n]*q**n for n in range(1, N+1))
t = mpf('0.3'); print("Fricke check f(i/8t) - 64 t^4 f(it):", nstr(fit(1/(8*t)) - 64*t**4*fit(t), 5))
t0 = 1/sqrt(8)
def Lam(s):
    return sum(c[n]*((2*pi*n)**(-s)*gammainc(s, 2*pi*n*t0) + mpf(8)**(2-s)*(2*pi*n)**(s-4)*gammainc(4-s, 2*pi*n*t0)) for n in range(1, 120))
L = {s: Lam(s)*(2*pi)**s/gamma(s) for s in (1, 2, 3)}
# direct integration cross-check (no functional equation): L(f,s) = (2pi)^s/Gamma(s) int_0^inf f(it) t^{s-1} dt,
# with f(it) for small t taken from Fricke (independent check of the integral, not of the sign)
Ld = {s: (2*pi)**s/gamma(s)*quad(lambda t: fit(t)*t**(s-1), [t0, 1, 3, 10, inf]) + (2*pi)**s/gamma(s)*quad(lambda u: 64*u**4*fit(u)*(1/(8*u))**(s-1)*(1/(8*u**2)), [t0,1,3,10,inf]) for s in (1,2,3)}
for s in (1,2,3): print("L(f,%d) =" % s, nstr(L[s], 30), " direct:", nstr(L[s]-Ld[s], 3))
a = lambda n: (gamma(n+mpf(1)/2)/(sqrt(pi)*gamma(n+1)))**2
h = lambda n: psi(0, n+1)*0 + harmonic(2*n) - harmonic(n) - log(2)
S0 = hyper([0.5]*4, [1]*3, 1)
S0b = nsum(lambda n: a(n)**2, [0, inf])
S1 = nsum(lambda n: a(n)**2*h(n), [0, inf])
print("S0 =", nstr(S0, 30), " (nsum diff", nstr(S0-S0b, 3), ")")
print("S1 =", nstr(S1, 30))
# Taylor coefficients of F(eps) = (1/pi^2) sum_n Gamma(1/2+eps+n)^4/Gamma(1+eps+n)^4 via polygammas
def coeffs(n):
    g0 = 4*(loggamma(n+mpf(1)/2) - loggamma(n+1))
    g1 = 4*(psi(0, n+mpf(1)/2) - psi(0, n+1))
    g2 = 4*(psi(1, n+mpf(1)/2) - psi(1, n+1))
    g3 = 4*(psi(2, n+mpf(1)/2) - psi(2, n+1))
    e = exp(g0)
    return [e, e*g1, e*(g2/2 + g1**2/2), e*(g3/6 + g1*g2/2 + g1**3/6)]
F = [nsum(lambda n, j=j: coeffs(n)[j], [0, inf])/pi**2 for j in range(4)]
for j in range(4): print("F_%d(1) =" % j, nstr(F[j], 30))
# Independent cross-check of F_j: Cauchy integral of F(eps) on a circle |eps| = 1/8, F evaluated by nsum
import mpmath
def Fe(e): return nsum(lambda n: exp(4*(loggamma(n+mpf(1)/2+e) - loggamma(n+1+e))), [0, inf])/pi**2
mp.dps = 25
M = 24; r = mpf(1)/8
vals = [Fe(r*mpmath.expjpi(2*mpf(k)/M)) for k in range(M)]
Fc = [sum(vals[k]*mpmath.expjpi(-2*mpf(j*k)/M) for k in range(M))/M/r**j for j in range(4)]
mp.dps = 40
for j in range(4): print("  Cauchy F_%d diff:" % j, nstr(Fc[j].real - F[j], 3), nstr(Fc[j].imag, 3))
print("checks:")
print(" L1 + S1/2           =", nstr(L[1] + S1/2, 3))
print(" L2 - pi^2 S0/16     =", nstr(L[2] - pi**2*S0/16, 3))
print(" L3 - pi^2 L1/4      =", nstr(L[3] - pi**2*L[1]/4, 3))
print(" F0 - S0             =", nstr(F[0]-S0, 3))
print(" F1 - 8 S1           =", nstr(F[1]-8*S1, 3))
print(" F2 - 2pi^2 S0       =", nstr(F[2]-2*pi**2*S0, 3))
print(" F3 - (4pi^2/3)F1    =", nstr(F[3]-4*pi**2/3*F[1], 3))
print(" L1/F1 + 1/16        =", nstr(L[1]/F[1] + mpf(1)/16, 3))
print(" L2/F2 - 1/32        =", nstr(L[2]/F[2] - mpf(1)/32, 3))
print(" L3/F3 + 3/256       =", nstr(L[3]/F[3] + mpf(3)/256, 3))
# Lemma 1 / Lemma 2 integrals (mpmath ellipk takes parameter m = k^2)
def _safe(m):
    v = ellipk(m)
    return v if mpmath.isfinite(v) else mpf(0)   # endpoint node rounding to the log singularity; contribution negligible
K = lambda k: _safe(k**2); Kp = lambda k: _safe((1-k)*(1+k))
IKp2 = quad(lambda k: Kp(k)**2, [0, 1]); IKKp = quad(lambda k: K(k)*Kp(k), [0, 1])
print(" int K'^2 + pi^2 S1  =", nstr(IKp2 + pi**2*S1, 3))
print(" int K K' - pi^3S0/8 =", nstr(IKKp - pi**3*S0/8, 3))
print(" int K^2 - int K'^2/2 =", nstr(quad(lambda k: K(k)**2, [0,1]) - IKp2/2, 3))
for n in (0, 3):
    Mn = quad(lambda k: Kp(k)*k**(2*n), [0, 1]); Ml = quad(lambda k: Kp(k)*k**(2*n)*log(1/k), [0, 1])
    print(" M(%d) check:" % (2*n), nstr(Mn - pi**2/4*a(n), 3), " M' check:", nstr(Ml + pi**2/2*a(n)*h(n), 3))
# DLMF 19.12.1 check at k = 0.3
k = mpf('0.3'); ser = nsum(lambda n: a(n)*k**(2*n)*(log(1/k) - 2*h(n)), [0, inf])
print(" K'(k) series (Z3)   =", nstr(Kp(k) - ser, 3))
# dy/dk check at k=0.4
k = mpf('0.4'); dy = mpmath.diff(lambda x: Kp(x)/K(x), k)
print(" dy/dk check         =", nstr(dy + pi/(2*k*(1-k**2)*K(k)**2), 3))
# f(iy/2) = k^2 k' K^4/pi^4 at y = K'/K for k=0.4
y = Kp(k)/K(k); print(" f(iy/2) check       =", nstr(fit(y/2) - k**2*sqrt(1-k**2)*K(k)**4/pi**4, 3))
