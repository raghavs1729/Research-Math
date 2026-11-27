"""Referee test 3 (independent): Proposition 2.1 of TOY.md (period doubling).

eps_n = 1 iff nu_2(n) is odd. Claims checked:
 (a) w^(2n) = 1 - w^(n) for all n (w^(0) = 1/2), w^(odd) = 0;
 (b) Fourier coefficients of Q_M = 2 pi M (R_M delta_0 + R_M delta_pi) are M(1+(-1)^n)/(M+in),
     with the closed form 2 pi M e^{-Mt}/(1-e^{-2 pi M}) for 2 pi M R_M delta_0 on (0,2pi)
     (delta normalised so that delta^(n) = 1/2pi);
 (c) coefficientwise: 1 + sum_k (-1)^k (Q_{M/2^k}(2^k t) - 2) has the coefficients 2M eps_|n|/(M+in);
 (d) pointwise at random non-dyadic t: the expansion agrees with Fejer means of the Fourier
     series of F_M, and F_M(t) + F_{M/2}(2t) = Q_M(t);
 (e) sup_t |Q_{M'} - 2| = O(M'), sup |F_M| finite, and ||F_M - 1||_2^2 from the expansion equals
     sum_{n != 0} eps_|n| 4M^2/(M^2+n^2) (the L^2 identity that fixes the constant 1).
"""
import math

import numpy as np

out = []


def log(s=""):
    print(s)
    out.append(s)


TWO_PI = 2 * math.pi


def nu2(n):
    n = np.asarray(n)
    v = np.zeros_like(n)
    m = n.copy()
    while True:
        even = (m % 2 == 0) & (m != 0)
        if not even.any():
            break
        v[even] += 1
        m[even] //= 2
    return v


n = np.arange(1, 10 ** 6 + 1)
eps = (nu2(n) % 2 == 1).astype(int)
eps_of = lambda k: (nu2(np.asarray(k)) % 2 == 1).astype(int)
ok_a1 = np.all(eps_of(2 * n[: 5 * 10 ** 5]) == 1 - eps[: 5 * 10 ** 5])
ok_a2 = np.all(eps[n % 2 == 1] == 0)
log(f"(a) eps(2n) = 1 - eps(n) for 1 <= n <= 5e5: {ok_a1};  eps(odd) = 0: {ok_a2};  (n = 0: 1/2 = 1 - 1/2)")


def Q(Mp, t):
    t = np.mod(t, TWO_PI)
    c = TWO_PI * Mp / (-math.expm1(-TWO_PI * Mp))
    return c * (np.exp(-Mp * t) + np.exp(-Mp * np.mod(t - math.pi, TWO_PI)))


# (b) Fourier coefficients of the closed form by exact integration on the pieces
def Q_coeff_exact(Mp, k):
    # (1/2pi) int_0^{2pi} Q(t) e^{-ikt} dt, computed piecewise in closed form
    c = TWO_PI * Mp / (-math.expm1(-TWO_PI * Mp))
    z = Mp + 1j * k
    I0 = (1 - np.exp(-TWO_PI * z)) / z  # int_0^{2pi} e^{-Mt} e^{-ikt}
    # int_0^{2pi} e^{-M((t-pi) mod 2pi)} e^{-ikt} dt = e^{-ik pi} * I0
    return c * (I0 + np.exp(-1j * k * math.pi) * I0) / TWO_PI


errb = 0.0
for Mp in [0.1, 0.5, 1.0, 3.0]:
    ks = np.arange(-50, 51)
    pred = Mp * (1 + (-1.0) ** ks) / (Mp + 1j * ks)
    got = np.array([Q_coeff_exact(Mp, k) for k in ks])
    errb = max(errb, np.max(np.abs(pred - got)))
log(f"(b) max |coefficients of the closed form of Q_M - M(1+(-1)^n)/(M+in)|, |n| <= 50, 4 values of M: {errb:.2e}")

# (c) coefficientwise identity
errc = 0.0
for M in [0.25, 0.5, 1.0, 2.0]:
    ns = np.arange(1, 200001)
    lhs = np.zeros(len(ns), dtype=complex)
    for k in range(0, 20):
        Mk = M / 2 ** k
        mask = ns % (2 ** k) == 0
        j = ns[mask] // 2 ** k
        lhs[mask] += (-1) ** k * Mk * (1 + (-1.0) ** j) / (Mk + 1j * j)
    rhs = 2 * M * eps_of(ns) / (M + 1j * ns)
    errc = max(errc, np.max(np.abs(lhs - rhs)))
log(f"(c) max over n <= 2e5 and M in {{.25,.5,1,2}} of |coefficient of the expansion - 2M eps_n/(M+in)|: {errc:.2e}")
log("    (constant term: Q^_{M'}(0) = 2 for every M', so the constant of the expansion is 1)")


def F_expansion(M, t, K=60):
    tot = np.ones_like(t, dtype=float)
    for k in range(K):
        Mk = M / 2 ** k
        tot += (-1) ** k * (Q(Mk, (2 ** k) * t) - 2)
    return tot


def F_fejer(M, t, N=200000):
    ns = np.arange(1, N + 1)
    c = eps_of(ns) * (1 - ns / (N + 1)) * 2 * M / (M + 1j * ns)
    return np.array([1 + 2 * np.real(np.sum(c * np.exp(1j * ns * tt))) for tt in np.atleast_1d(t)])


rng = np.random.default_rng(5)
ts = rng.uniform(0, TWO_PI, 8)
for M in [0.5, 2.0]:
    e = F_expansion(M, ts)
    f = F_fejer(M, ts)
    lhs = F_expansion(M, ts) + F_expansion(M / 2, 2 * ts)
    log(f"(d) M={M}: max |expansion - Fejer(N=2e5)| at 8 random t: {np.max(np.abs(e - f)):.2e};"
        f"  max |F_M(t) + F_(M/2)(2t) - Q_M(t)|: {np.max(np.abs(lhs - Q(M, ts))):.2e}")

# (e) uniform bounds and the L^2 identity
for Mp in [1e-1, 1e-2, 1e-3, 1e-4]:
    tt = np.linspace(0, TWO_PI, 200001)
    log(f"(e) M' = {Mp:g}: sup |Q_M' - 2| / M' = {np.max(np.abs(Q(Mp, tt) - 2)) / Mp:.4f}")
for M in [0.5, 2.0]:
    K = 2 ** 22
    tt = (np.arange(K) + 0.5) * TWO_PI / K
    Fv = F_expansion(M, tt, K=40)
    l2_num = np.mean((Fv - 1) ** 2)
    ns = np.arange(1, 10 ** 7 + 1)
    l2_pred = 2 * np.sum(eps_of(ns) * 4 * M ** 2 / (M ** 2 + ns ** 2.0))
    log(f"(e) M={M}: mean of expansion = {np.mean(Fv):.6f} (predicted 1); ||F-1||^2: grid {l2_num:.6f} vs "
        f"Parseval {l2_pred:.6f}; bound 8 M^2 zeta(2) = {8 * M * M * math.pi ** 2 / 6:.4f}; sup|F_M| = {np.max(np.abs(Fv)):.4f}")
# jump of F_M at t = 0: predicted 2M * (mass of w at 0) = 2M * sum_k (-1)^k pi 2^{-k} = 4 pi M / 3
for M in [0.5, 2.0]:
    h = 1e-9
    jump = F_expansion(M, np.array([h]))[0] - F_expansion(M, np.array([TWO_PI - h]))[0]
    # jump at pi j/2^a (j odd, a >= 1): sum_{k>=a} (-1)^k 2 pi M 2^{-k} = (-1)^a (4 pi M/3) 2^{-a}
    x = math.pi * 3 / 8  # a = 3
    jx = F_expansion(M, np.array([x + h]))[0] - F_expansion(M, np.array([x - h]))[0]
    log(f"    jump of F_M at 0 (M={M}): {jump:.6f}, predicted 4 pi M/3 = {4 * math.pi * M / 3:.6f};"
        f" at 3pi/8: {jx:.6f}, predicted -(4 pi M/3)/8 = {-(4 * math.pi * M / 3) / 8:.6f}")
log("    => the 2^a points of dyadic level a carry total jump 4 pi M/3 for every a: F_M is not BV, w is not a measure")

with open(__file__.replace(".py", "_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
