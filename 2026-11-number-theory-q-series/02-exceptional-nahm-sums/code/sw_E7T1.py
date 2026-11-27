"""sw_E7T1.py -- numerical test (NOT a proof) of a conjectural formula for the pair (E7,T1) of Sun-Wang
(arXiv:2604.00847: "for Nahm sums of type (E6,T1) and (E7,T1), we have not found good CFT interpretations"):

    sum_{n in N^7} q^{n^T C(E7) n / 2} / prod (q)_{n_i}  =  ( Theta_{E6}(q) + Theta_{A6}(q) ) / ( 2 (q;q)_inf^6 ),

Theta_L(q) = sum_{x in L} q^{|x|^2/2} for the root lattices E6, A6; i.e. q^{-1/4} * Nahm = (chi_0^{(E6)_1} + chi_0^{(A6)_1})/2.
(Compare (T1,E7) = (Theta_{<6>} + Theta_{<14>})/(2 eta), proofs/E7.md; <6>, <14> are the orthogonal complements of E6, A6 in E7.)
Also prints the first coefficients of the (E6,T1) sum (no formula found).
usage: python3 sw_E7T1.py N
"""
import sys, time
from fractions import Fraction as F
from sw_lib import Ctx, compare, report
from fusion import cartanE

N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
c = Ctx(3, N)


def theta_An(n):
    """sum over x in Z^{n+1}, sum x = 0, of q^{|x|^2/2}"""
    R = int((2 * N) ** 0.5) + 2
    P = 2 * N + 2
    dp = {0: [1] + [0] * (P - 1)}
    for i in range(n + 1):
        nd = {}
        for s, arr in dp.items():
            for x in range(-R, R + 1):
                e = x * x
                if e >= P:
                    continue
                t = nd.setdefault(s + x, [0] * P)
                for j in range(P - e):
                    if arr[j]:
                        t[j + e] += arr[j]
        dp = nd
    arr = dp[0]
    return c.from_dict({F(j, 2): v for j, v in enumerate(arr) if v and j % 2 == 0 and j // 2 < N})


# E6 = A2^3 with glue: Theta_E6 = a^3 + 2 c^3,  a = sum q^{m^2+mn+n^2},  c = sum q^{m^2+mn+n^2+m+n+1/3}
da, dc = {}, {}
R = int(N ** 0.5) * 2 + 4
for m in range(-R, R + 1):
    for n in range(-R, R + 1):
        e = m * m + m * n + n * n
        if e < N:
            da[F(e)] = da.get(F(e), 0) + 1
        e2 = F(3 * (m * m + m * n + n * n + m + n) + 1, 3)
        if e2 < N:
            dc[e2] = dc.get(e2, 0) + 1
A, Cc = c.from_dict(da), c.from_dict(dc)
ThE6 = A * A * A + 2 * (Cc * Cc * Cc)
ThA6 = theta_An(6)
print("Theta_E6:", ThE6.lead(5))
print("Theta_A6:", ThA6.lead(5))
e = c.euler_inv(); e6 = e * e * e; e6 = e6 * e6
t0 = time.time()
f = c.nahm(cartanE(7), [0] * 7)
print("(E7,T1) Nahm sum:", f.lead(8), "  [%.0fs]" % (time.time() - t0), flush=True)
report("2 * Nahm(C(E7)) = (Theta_E6 + Theta_A6)/(q)_inf^6", f * 2, (ThE6 + ThA6) * e6)
c1 = Ctx(1, min(N, 24))
t0 = time.time()
g = c1.nahm(cartanE(6), [0] * 6)
cs = g.coeffs()
print("(E6,T1) Nahm sum (no formula found):", [cs.get(i, 0) for i in range(min(N, 24))], "  [%.0fs]" % (time.time() - t0))
