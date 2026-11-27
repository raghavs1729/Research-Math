"""Referee check of (1.1): K_k from the defining 2k-fold sum (brute-force enumeration) vs omega(T^{k-1}u)
(ref_series.K_direct), as bivariate series mod x^DX, q^DQ."""
import sys, itertools
import ref_series as RS
from flint import fmpq_poly


def brute(k):
    DX, DQ = RS.DX, RS.DQ
    tot = RS.zero()
    # r_1 >= ... >= r_k >= 0, 0 <= s_1 <= ... <= s_k, sum (r_i + 2 s_i) <= DX
    def rec_r(i, prev, acc):
        if i == k:
            yield tuple(acc)
            return
        for r in range(prev, -1, -1):
            yield from rec_r(i + 1, r, acc + [r])
    def rec_s(i, prev, acc, budget):
        if i == k:
            yield tuple(acc)
            return
        for s in range(prev, budget + 1):
            yield from rec_s(i + 1, s, acc + [s], budget)
    for rs in rec_r(0, DX, []):
        if sum(rs) > DX:
            continue
        for ss in rec_s(0, 0, [], DX // 2):
            wx = sum(rs) + 2 * sum(ss)
            if wx > DX:
                continue
            eq = sum((rs[i] + ss[i])**2 + ss[i]**2 for i in range(k))
            if eq >= DQ:
                continue
            p = fmpq_poly([0] * eq + [1])
            for i in range(k):
                rn = rs[i + 1] if i + 1 < k else 0
                sp_ = ss[i - 1] if i > 0 else 0
                p = RS.trunc(p * RS.IQ[rs[i] - rn])
                p = RS.trunc(p * RS.IQQ[ss[i] - sp_])
            tot[wx] = RS.trunc(tot[wx] + p)
    return tot


if __name__ == "__main__":
    RS.DX, RS.DQ = 10, 45
    RS.IQ[:] = [RS.qpoch_inv(1, 1, n) for n in range(RS.DX + 2)]
    RS.IQQ[:] = [RS.qpoch_inv(2, 2, n) for n in range(RS.DX + 2)]
    for k in range(1, 6):
        print(f"k={k}: brute-force 2k-fold sum == omega(T^(k-1) u) mod x^{RS.DX} q^{RS.DQ}:",
              RS.equal(brute(k), RS.K_direct(k)))
