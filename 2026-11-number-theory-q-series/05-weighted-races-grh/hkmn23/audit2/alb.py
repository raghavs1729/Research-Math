"""Second audit: the explicit bound LB(q) of section 4.4 (Theorem 5(b), q > 400), my own sieve.
LB(q) = phi/2 + 1 - 2 log(4 L+) - (log q)/3 - (phi/q^2)(zeta(2) log q + 1.81) - R(q),  L+ = phi(3 log q + 0.11),
R(q) = 2[(n_e + omega phi)(log X1 + 1) + |B|]/X1 + [n_o(log X1 + 1/2) + |C|]/X1^2 + 0.4(1+n_e)/X1^3 + 0.1 n_o/X1^4,
X1 = 16 (L+)^2, n_e = phi/2 - 1, n_o = phi/2, |B| <= phi(2 log q + 4.4), |C| <= phi(log q + 2.4).
Float64 (margins are >= 18 in absolute terms).  Also checks the q > 10^6 argument on extreme moduli.
"""
import numpy as np, math

Q = 10 ** 6
spf = np.zeros(Q + 1, dtype=np.int64)
for i in range(2, Q + 1):
    if spf[i] == 0:
        spf[i::i][spf[i::i] == 0] = i
phi = np.zeros(Q + 1); om = np.zeros(Q + 1, dtype=np.int64); N = np.zeros(Q + 1)
for q in range(1, Q + 1):
    n = q; ph = 1; w = 0; nn = 1
    while n > 1:
        p = int(spf[n]); e = 0
        while n % p == 0:
            n //= p; e += 1
        ph *= (p - 1) * p ** (e - 1); w += 1
        if p == 2:
            nn *= 1 if e == 1 else (2 if e == 2 else 4)
        else:
            nn *= 2
    phi[q] = ph; om[q] = w; N[q] = nn

z2 = math.pi ** 2 / 6


def LB(q, ph, w):
    lq = math.log(q)
    Lp = ph * (3 * lq + 0.11)
    X1 = 16 * Lp * Lp
    ne = ph / 2 - 1; no = ph / 2
    B = ph * (2 * lq + 4.4); C = ph * (lq + 2.4)
    R = (2 * ((ne + w * ph) * (math.log(X1) + 1) + B) / X1 + (no * (math.log(X1) + 0.5) + C) / X1 ** 2
         + 0.4 * (1 + ne) / X1 ** 3 + 0.1 * no / X1 ** 4)
    small = ph / q ** 2 * (z2 * lq + 1.81)
    return ph / 2 + 1 - 2 * math.log(4 * Lp) - lq / 3 - small - R, R, small


worst = (1e9, 0); Rmax = 0; smax = 0
for q in range(401, Q + 1):
    v, R, s = LB(q, phi[q], om[q])
    d = v - 2 / 3 * N[q]
    Rmax = max(Rmax, R); smax = max(smax, s)
    if d < worst[0]:
        worst = (d, q)
print(f'400<q<=1e6: min LB-(2/3)N = {worst[0]:.4f} at q={worst[1]} (phi={phi[worst[1]]:.0f}, N={N[worst[1]]:.0f}), max R = {Rmax:.3e}, max small = {smax:.3e}')
# exceptions to N <= phi/4
exc = [q for q in range(1, Q + 1) if N[q] > phi[q] / 4]
print('q with N(q,1) > phi/4:', exc)
# q > 1e6: worst case is large N/phi; test primorial-type moduli and 2^k*primorial
for q in [2 * 3 * 5 * 7 * 11 * 13 * 17 * 19, 8 * 3 * 5 * 7 * 11 * 13 * 17 * 19, 8 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37,
          10 ** 6 + 3, 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47]:
    n = q; ph = 1; w = 0; nn = 1; p = 2
    m = q; fac = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            fac[d] = fac.get(d, 0) + 1; m //= d
        d += 1
    if m > 1:
        fac[m] = fac.get(m, 0) + 1
    for p, e in fac.items():
        ph *= (p - 1) * p ** (e - 1); w += 1
        nn *= (2 if p > 2 else (1 if e == 1 else 2 if e == 2 else 4))
    v, R, s = LB(q, ph, w)
    print(f'q={q}: phi={ph} N={nn} LB={v:.4g} (2/3)N={2*nn/3:.4g}  phi/6={ph/6:.4g} R={R:.2e}')
# g(phi) of the note
g = lambda f: f / 3 + 1 - 2 * math.log(4 * f * (6 * math.log(f) + 2.2)) - (2 * math.log(f) + 0.7) / 3 - 0.04
print('g(707) =', g(707), ' min over 20..1e6 of g\' sign check:', all(g(f + 1) > g(f) for f in range(20, 100000)))
