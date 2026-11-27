from mpmath import *
from case_half_quarter import Fj, eta_quot, NQ
mp.dps = 40
F = Fj(mpf(1)/2, mpf(1)/3)
print("F_j(1), (1/2,1/3):", [nstr(x, 25) for x in F])
print("pslq [F3, pi^2 F1]:", pslq([F[3], pi**2 * F[1]], maxcoeff=10**4))
# candidate form: g = eta2^6 eta6^2 + 9 eta2^3 eta6^2 eta18^3
s1, c1 = eta_quot({2: 6, 6: 2}); s2, c2 = eta_quot({2: 3, 6: 2, 18: 3})
coef = {}
for n, c in enumerate(c1):
    if c and n + s1 <= NQ: coef[n + s1] = coef.get(n + s1, 0) + c
for n, c in enumerate(c2):
    if c and n + s2 <= NQ: coef[n + s2] = coef.get(n + s2, 0) + 9 * c
print("g coefficients:", [coef.get(n, 0) for n in range(1, 30)])
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
for lev in (36, 144, 12, 72):
    for sign in (1, -1):
        L = Lfun(coef, lev, sign)
        # consistency of the assumed functional equation: compare with a second evaluation using a different split point is not available here; test via L(s) at s=1,2,3 ratios
        v = [L(1), L(2), L(3)]
        print(lev, sign, [nstr(x, 15) for x in v], "ratios L1/F1, L2/F2, L3/F3:", [nstr(v[i] / F[i + 1], 12) for i in range(3)], " L2/(pi F1), L1 pi/F2:", nstr(v[1] / (pi * F[1]), 12), nstr(v[0] * pi / F[2], 12))
