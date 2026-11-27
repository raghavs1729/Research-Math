from mpmath import *
from case_half_quarter import Fj, eta_quot, NQ
mp.dps = 40
def Lfun(coef, lev, sign, k=4):
    A = sqrt(mpf(lev)) / (2 * pi)
    def L(s):
        tot = mpf(0)
        for n, c in sorted(coef.items()):
            x = n / A
            if x > 110: break
            tot += c * (gammainc(s, x) / (mpf(n) / A)**s + sign * gammainc(k - s, x) / (mpf(n) / A)**(k - s))
        return tot / (A**s * gamma(s))
    return L
def addq(*terms):
    coef = {}
    for mult, expo in terms:
        s, c = eta_quot(expo)
        for n, v in enumerate(c):
            if v and n + s <= NQ: coef[n + s] = coef.get(n + s, 0) + mult * v
    return coef
# (1/2,1/3)
g = addq((1, {2: 6, 6: 2}), (9, {2: 3, 6: 2, 18: 3}))
chi3 = lambda n: [0, 1, -1][n % 3]
f36 = {n: c * chi3(n) for n, c in g.items() if chi3(n)}
print("g~   :", [g.get(n, 0) for n in range(1, 26)])
print("f36  :", [f36.get(n, 0) for n in range(1, 26)])
F = Fj(mpf(1)/2, mpf(1)/3)
Lg = Lfun(g, 12, 1); Lg1, Lg2, Lg3 = Lg(1), Lg(2), Lg(3)
for sign in (1, -1):
    Lf = Lfun(f36, 36, sign); v = [Lf(1), Lf(2), Lf(3)]
    print("level 36 sign %+d: L(f,k) = %s ; L/F_k: %s" % (sign, [nstr(x, 15) for x in v], [nstr(v[i] / F[i + 1], 15) for i in range(3)]))
    print("   L(f,1) pi/(sqrt3 L(g,2)) = %s ; L(f,2)/(sqrt3 pi L(g,1)) = %s ; L(f,3)/(pi sqrt3 L(g,2)) = %s" % (nstr(v[0] * pi / (sqrt(3) * Lg2), 15), nstr(v[1] / (sqrt(3) * pi * Lg1), 15), nstr(v[2] / (pi * sqrt(3) * Lg2), 15)))
print("L(g,3)/(pi^2 L(g,1)) =", nstr(Lg3 / (pi**2 * Lg1), 15))
# direct Dirichlet-series sanity for L(f36,3)
print("direct sum L(f36,3) ~", nstr(sum(mpf(c) / n**3 for n, c in f36.items()), 8), "; L(g,3) direct ~", nstr(sum(mpf(c) / n**3 for n, c in g.items()), 8))
# multiplier nu(gamma) for g0 = a(q) eta^3 eta3^3 under gamma = [[1,0],[-3,1]]
def eta(t): 
    q = exp(2j * pi * t)
    return exp(2j * pi * t / 24) * nprod(lambda n: 1 - q**n, [1, inf])
def acub(t):
    q = exp(2j * pi * t)
    return nsum(lambda m, n: q**(m*m + m*n + n*n), [-25, 25], [-25, 25])
mp.dps = 20
g0 = lambda t: acub(t) * eta(t)**3 * eta(3 * t)**3
t0 = mpc('0.13', '0.71')
print("g0(gamma t)/((1-3t)^4 g0(t)) =", nstr(g0(t0 / (1 - 3 * t0)) / ((1 - 3 * t0)**4 * g0(t0)), 10), "; g0(t+1)/g0(t) =", nstr(g0(t0 + 1) / g0(t0), 10), "; g0(-1/(3t))/(9 t^4 g0(t)) =", nstr(g0(-1 / (3 * t0)) / (9 * t0**4 * g0(t0)), 10))
