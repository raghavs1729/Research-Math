"""Counterexamples to Conjecture 1 of arXiv:2412.12361v2 as printed (convergence of C[f_n] = 1 + f_1/(1 + f_2/(1 + ...))).
(A) bullet 2 with C = 1 + 4f < 0:  f_n = f constant < -1/4.  The N-th approximant is T^N(1), T(x) = 1 + f/x, an elliptic
    Moebius map (trace^2/det = -1/f in (0,4)); the orbit of the real point 1 does not converge.
    Exact: f = -1: approximants 0, oo, 1, 0, oo, 1, ...;  f = -1/2: 1/2, 0, oo, 1, 1/2, ...;  f = -1/3: period 6.
(B) bullet 4 with k = -2:  f_n = -1/4 - c/(n(n+1)),  c = 1/16 + beta^2/4 (beta > 0):  1 + 4 f_n = -4c/(n(n+1)) = O(n^-2).
    u_n = 2^-n Gamma(n+2+sigma)/Gamma(n+2), sigma = 1/2 + i beta, solves u_n = u_{n-1} + f_n u_{n-2} (checked numerically),
    P_n/Q_n = Re(a u_n)/Re(b u_n) oscillates like a function of beta*ln n: no limit.  Approximants printed for c = 1/4.
(C) control: c = 1/16 (beta = 0) and c = 1/32 converge."""
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, gamma, nstr
def approximants(f, N):
    P0, P1, Q0, Q1 = Fr(1), Fr(1), Fr(0), Fr(1); out = []
    for n in range(1, N + 1):
        P0, P1 = P1, P1 + f(n) * P0; Q0, Q1 = Q1, Q1 + f(n) * Q0
        out.append(None if Q1 == 0 else P1 / Q1)
    return out
for f in (Fr(-1), Fr(-1, 2), Fr(-1, 3)):
    a = approximants(lambda n: f, 24)
    per = {Fr(-1): 3, Fr(-1, 2): 4, Fr(-1, 3): 6}[f]
    assert all(a[i] == a[i + per] for i in range(24 - per)) and len(set(map(str, a[:per]))) > 1
    print("(A) f = %s: approximants %s ... (period %d): divergent" % (f, [str(x) if x is not None else 'oo' for x in a[:per + 1]], per))
mp.dps = 30
beta = mpf(3).sqrt() / 2; sigma = mpc(mpf(1) / 2, beta); c = mpf(1) / 16 + beta**2 / 4
assert abs(c - mpf(1) / 4) < mpf(10)**-25
u = lambda n: mpf(2)**(-n) * gamma(n + 2 + sigma) / gamma(n + 2)
res = max(abs(u(n) - u(n - 1) - (-mpf(1) / 4 - c / (n * (n + 1))) * u(n - 2)) / abs(u(n)) for n in range(1, 40))
print("(B) residual of the explicit complex solution (n < 40): %s" % nstr(res, 3)); assert res < mpf(10)**-25
def approx_float(f, Ns):
    P0, P1, Q0, Q1 = mpf(1), mpf(1), mpf(0), mpf(1); out = {}
    for n in range(1, max(Ns) + 1):
        fn = f(n); P0, P1 = P1, P1 + fn * P0; Q0, Q1 = Q1, Q1 + fn * Q0
        if n % 64 == 0:                       # rescale
            s = abs(Q1) + abs(P1); P0 /= s; P1 /= s; Q0 /= s; Q1 /= s
        if n in Ns: out[n] = P1 / Q1
    return out
Ns = [10**j for j in range(1, 7)] + [3 * 10**j for j in range(1, 6)]
for cc, label in ((mpf(1) / 4, "c = 1/4 (beta = sqrt(3)/2): DIVERGES"), (mpf(1) / 16, "c = 1/16: converges (slowly)"), (mpf(1) / 32, "c = 1/32: converges")):
    r = approx_float(lambda n: -mpf(1) / 4 - cc / (n * (n + 1)), Ns)
    print("    %-36s approximants at n = 10, 30, 100, ..., 10^6: %s" % (label, [nstr(r[n], 6) for n in sorted(Ns)]))
print("PASS")
