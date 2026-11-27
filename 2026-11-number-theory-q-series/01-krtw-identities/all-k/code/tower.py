"""Eigenfunction towers of the KRTW (iii)-chain transfer operator, in the Bailey basis relative to x (base q).

Phi_{N,M}(R;x) = sum_e c^{N,M}_e(x) / ((q)_{R-e} (xq)_{R+e})      (times 1/(x^2q^2;q^2)_inf)
eigen-equation (b = 0 term moved to the left):
 (lam_{N,M} - x^e q^{e^2}) c^{N,M}_e(x) = sum_{b=1}^{M} x^{2b} q^{2b^2} (x^2q^2;q^2)_{2b}/(q^2;q^2)_b
                                        * sum_d c^{N-b,M-b}_d(x q^{2b}) (x q^{2b})^d q^{d^2} Gamma^{(b)}_{d,e}(x)
seed: c^{L,0}_e = delta_{e,L}  (tower gamma^{L,M} = c^{L+M,M} with seed 1)
"""
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from functools import lru_cache


def qbin(n, k, step=1):
    if k < 0 or k > n:
        return RF(0)
    return qpoch(n, step) / (qpoch(k, step) * qpoch(n - k, step))


@lru_cache(maxsize=None)
def Gamma(b, d, e):
    """coefficient of 1/((q)_{R-e}(xq)_{R+e}) in 1/((q)_{R-d}(x q^{2b+1})_{R+d})"""
    n = e - d
    if n < 0 or n > 2 * b:
        return RF(0)
    if d + e == 0:
        return RF(1)
    num = RF(ONE - X * Q ** (2 * e)) * poch(X * Q, d + e - 1)          # (1-xq^{2e}) (x)_{d+e}/(1-x)
    den = poch(X * Q ** (2 * b + 1), d + e)
    sgn = -1 if n % 2 else 1
    return num / den * mono(n, n * (n - 1) // 2 + (e + d) * n) * qbin(2 * b, n) * sgn


def lam(N, M):
    return mono(N + M, N * N + M * M)


_tow = {}


def tower(L, M):
    """gamma^{L,M}: dict e -> RF, for the mode (N,M) = (L+M, M), seed 1"""
    key = (L, M)
    if key in _tow:
        return _tow[key]
    N = L + M
    if M == 0:
        res = {L: RF(1)}
        _tow[key] = res
        return res
    rhs = {}
    for b in range(1, M + 1):
        pref = mono(2 * b, 2 * b * b) * poch(X * X * Q ** 2, 2 * b, 2) / qpoch(b, 2)
        low = tower(L, M - b)          # mode (N-b, M-b) = (L + M - b, M - b)
        xs = X * Q ** (2 * b)
        for d, cd in low.items():
            cds = cd.subs_x(xs) * mono(d, 2 * b * d + d * d)
            for e in range(d, d + 2 * b + 1):
                g = Gamma(b, d, e)
                if g.iszero():
                    continue
                t = pref * cds * g
                rhs[e] = rhs.get(e, RF(0)) + t
    res = {}
    lm = lam(N, M)
    for e, v in rhs.items():
        if v.iszero():
            continue
        res[e] = v / (lm - mono(e, e * e))
    _tow[key] = res
    return res


if __name__ == '__main__':
    import sys
    for M in range(0, 3):
        for L in range(0, 3):
            t = tower(L, M)
            print('L=%d M=%d (mode N=%d,M=%d): support e in %s' % (L, M, L + M, M, sorted(t)))
            for e in sorted(t):
                print('   e=%d: %s' % (e, t[e].factor_str()))
            sys.stdout.flush()
