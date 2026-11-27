# K_3(1), K_3(q) (and K_2) directly as q-series mod q^N vs KRTW product sides; also B,J theta evaluations.
import sys
from flint import fmpz_poly
from ser import qpoch_inv_table
N = int(sys.argv[1])
t1 = qpoch_inv_table(1, N, N); t2 = qpoch_inv_table(2, N, N)
def Ksum(k, lin):  # lin: 0 for (iii), 1 for (ii)  (extra exponent r_i+2s_i)
    out = [0] * N; cnt = [0]
    import itertools
    def rec(i, rs, ss, e0):
        # choose r_i<=r_{i-1}, s_i>=s_{i-1}; partial exponent lower bound e0
        if i == k:
            e = e0
            p = fmpz_poly([1])
            for j in range(k):
                rn = rs[j + 1] if j + 1 < k else 0
                sp = ss[j - 1] if j else 0
                p = p.mul_low(t1[rs[j] - rn], N).mul_low(t2[ss[j] - sp], N)
            c = p.coeffs()
            for a, v in enumerate(c):
                if a + e < N: out[a + e] += int(v)
            cnt[0] += 1
            return
        rtop = rs[-1] if rs else N
        slo = ss[-1] if ss else 0
        for r in range(0, rtop + 1):
            s = slo
            while True:
                ee = e0 + (r + s) ** 2 + s * s + lin * (r + 2 * s)
                # remaining indices i+1..k-1 have s>=s, so add at least (k-1-i)*s^2*2
                if ee + (k - 1 - i) * (2 * s * s + lin * 2 * s) >= N: break
                rec(i + 1, rs + [r], ss + [s], ee)
                s += 1
            if e0 + r * r + lin * r >= N: break
    rec(0, [], [], 0)
    return out, cnt[0]
def prod(exps, mod, den):
    p = fmpz_poly([1])
    for a in exps:
        j = a
        while j < N:
            p = p.mul_low(fmpz_poly([1] + [0] * (j - 1) + [-1]), N); j += mod
    for (a, st) in den:
        j = a
        while j < N:
            p = p.mul_low(fmpz_poly([1 if i % j == 0 else 0 for i in range(N)]), N); j += st
    return [int(v) for v in p.coeffs()] + [0] * (N - p.length())
for k in (1, 2, 3, 4):
    m = 2 * k + 4
    den = [(2, 2), (1, 1)]
    ii = prod([1, 2, 2 * k + 2, 2 * k + 3, m, m], m, den)
    iii = prod([k + 1, k + 2, k + 2, k + 3, m, m], m, den)
    Kq, c1 = Ksum(k, 1); K1, c2 = Ksum(k, 0)
    print("k=%d mod q^%d: (ii) %s  (iii) %s   [terms %d, %d]" % (k, N, Kq == ii[:N], K1 == iii[:N], c1, c2))
