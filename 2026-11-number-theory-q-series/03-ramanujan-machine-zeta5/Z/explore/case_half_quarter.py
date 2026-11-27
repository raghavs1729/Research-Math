"""Exploration (numerical): Zudilin's F_j(1) for (r,t) = (1/2,1/4) against the L-values of f1 = eta(4t)^16/(eta(2t)^4 eta(8t)^4)
(level 16, weight 4) and against S_0, S_1 of the (1/2,1/2) case."""
from mpmath import mp, mpf, pi, gamma, taylor, nsum, inf, gammainc, sqrt, hyper, pslq, nstr, log, binomial, harmonic
mp.dps = 40
def Fj(r, t, dps=40):
    r = mpf(r); t = mpf(t)
    pref = 1 / (gamma(r) * gamma(1 - r) * gamma(t) * gamma(1 - t))
    def G(e):
        # sum_n Gamma(r+e+n)Gamma(1-r+e+n)Gamma(t+e+n)Gamma(1-t+e+n)/Gamma(1+e+n)^4  via hypergeometric closed form
        c = gamma(r + e) * gamma(1 - r + e) * gamma(t + e) * gamma(1 - t + e) / gamma(1 + e)**4
        return pref * c * hyper([r + e, 1 - r + e, t + e, 1 - t + e, 1], [1 + e] * 4, 1)
    mp.dps = dps + 40
    v = taylor(G, 0, 3)
    mp.dps = dps
    return [+x for x in v]
NQ = 600
def eta_quot(expo):       # expo: dict d -> power; returns coefficients of prod eta(d t)^e as q-series, with leading power
    lead = sum(d * e for d, e in expo.items())
    assert lead % 24 == 0
    co = [0] * (NQ + 1); co[0] = 1
    for d, e in expo.items():
        for m in range(1, NQ // d + 1):
            step = d * m
            for _ in range(abs(e)):
                if e > 0:
                    for i in range(NQ, step - 1, -1): co[i] -= co[i - step]
                else:
                    for i in range(step, NQ + 1): co[i] += co[i - step]
    return lead // 24, co
def Lvals(expo, lev, k=4, sign=1):
    sh, co = eta_quot(expo)
    coef = {n + sh: co[n] for n in range(NQ + 1 - sh) if co[n]}
    A = sqrt(mpf(lev)) / (2 * pi)
    def L(s):
        tot = mpf(0)
        for n, c in coef.items():
            x = n / A
            if x > 110: break
            tot += c * (gammainc(s, x) / (mpf(n) / A)**s + sign * gammainc(k - s, x) / (mpf(n) / A)**(k - s))
        return tot / (A**s * gamma(s))
    return [L(1), L(2), L(3)], [coef.get(n, 0) for n in range(1, 12)]
if __name__ == '__main__':
    F = Fj(mpf(1) / 2, mpf(1) / 4)
    print("F_j(1), (1/2,1/4):", [nstr(x, 25) for x in F])
    (L1, L2, L3), c = Lvals({4: 16, 2: -4, 8: -4}, 16)
    print("f1 coefficients:", c)
    print("L(f1,k):", nstr(L1, 25), nstr(L2, 25), nstr(L3, 25))
    for kk, (L, Fv) in enumerate(zip((L1, L2, L3), F[1:]), 1):
        print("  L(f1,%d)/F_%d(1) = %s" % (kk, kk, nstr(L / Fv, 25)))
    # relations among F_0..F_3 with pi powers
    print("pslq [F3, pi^2 F1]:", pslq([F[3], pi**2 * F[1]], maxcoeff=10**4))
    print("pslq [F2, pi^2 F0, pi^2?]:", pslq([F[2], pi**2 * F[0], mpf(1), pi**2], maxcoeff=10**4))
    S0 = hyper([mpf(1) / 2] * 4, [1, 1, 1], 1)
    S1 = nsum(lambda n: (binomial(2 * n, n) / mpf(4)**n)**4 * (harmonic(2 * n) - harmonic(n) - log(2)), [0, inf])
    for name, v in (("F0", F[0]), ("F1", F[1]), ("F2/pi^2", F[2] / pi**2), ("F3/pi^2", F[3] / pi**2)):
        print("pslq [%s, S0, S1, pi S0, pi S1, S0/pi, S1/pi]:" % name, pslq([v, S0, S1, pi * S0, pi * S1, S0 / pi, S1 / pi], maxcoeff=10**4, maxsteps=10**5))
