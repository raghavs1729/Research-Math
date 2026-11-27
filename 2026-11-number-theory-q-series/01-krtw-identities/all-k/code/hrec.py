"""Genuine eigenfunctions of T (top normalised) via the local H'-recursion, in A = R+S coordinates.

F(A,B) = sum_n kappa_n(B) phi_n(A;x),  phi_n(A;x) = 1/((q)_{A-n} (xq)_{A+n}).
kappa_n(M) = delta_{n,N};  for B < M:
  D_{n,B} kappa_n(B) = x^2 q^B [ alpha_{n+1} kappa_{n+1}(B+1) + beta_n kappa_n(B+1) + gamma_{n-1} kappa_{n-1}(B+1) ],
  D_{n,B} = (q^{-B}+xq^B)(q^{-n}+xq^n) - (q^{-N}+xq^N)(q^{-M}+xq^M),
  alpha_n = q^{n-1}/((1-xq^{2n-1})(1-xq^{2n})), beta_n = -q^{n-1}(1+xq^{2n})/((1-xq^{2n-1})(1-xq^{2n+1})),
  gamma_n = x q^{3n}/((1-xq^{2n})(1-xq^{2n+1})).
Path-sum decomposition by (#down, #zero, #up) steps.
"""
import sys
from ratfun import RF, X, Q, ONE, mono, poch, qpoch


def h(n):          # q^{-n} + x q^n  (as RF)
    return mono(0, -n) + mono(1, n)


def alpha(n):
    return mono(0, n - 1) / (RF(ONE - X * Q ** (2 * n - 1)) * RF(ONE - X * Q ** (2 * n)))


def beta(n):
    return RF(-1) * mono(0, n - 1) * RF(ONE + X * Q ** (2 * n)) / (RF(ONE - X * Q ** (2 * n - 1)) * RF(ONE - X * Q ** (2 * n + 1)))


def gamma(n):
    return mono(1, 3 * n) / (RF(ONE - X * Q ** (2 * n)) * RF(ONE - X * Q ** (2 * n + 1)))


def eig(N, M, classes=False):
    """returns dict (n,B) -> kappa; if classes, dict (n,B) -> {(nd,n0,nu): contribution}"""
    eps = h(N) * h(M)
    cur = {N: {(0, 0, 0): RF(1)}}
    out = {(N, M): dict(cur[N])}
    for B in range(M - 1, -1, -1):
        new = {}
        for n in range(N - (M - B), N + (M - B) + 1):
            if n < 0:
                continue
            D = h(B) * h(n) - eps
            pre = mono(2, B) / D
            acc = {}
            # from n+1 (down step in n): alpha_{n+1}
            if (n + 1) in cur:
                for (a, b, c), v in cur[n + 1].items():
                    key = (a + 1, b, c)
                    acc[key] = acc.get(key, RF(0)) + v * alpha(n + 1)
            if n in cur:
                for (a, b, c), v in cur[n].items():
                    key = (a, b + 1, c)
                    acc[key] = acc.get(key, RF(0)) + v * beta(n)
            if (n - 1) in cur:
                for (a, b, c), v in cur[n - 1].items():
                    key = (a, b, c + 1)
                    acc[key] = acc.get(key, RF(0)) + v * gamma(n - 1)
            res = {k: pre * v for k, v in acc.items() if not v.iszero()}
            if res:
                new[n] = res
                out[(n, B)] = res
        cur = new
    if classes:
        return out
    return {k: sum(v.values(), RF(0)) for k, v in out.items()}


if __name__ == '__main__':
    N = int(sys.argv[1]); M = int(sys.argv[2])
    out = eig(N, M, classes=True)
    for B in range(M, -1, -1):
        for n in range(N - (M - B), N + (M - B) + 1):
            if (n, B) not in out:
                continue
            cl = out[(n, B)]
            tot = sum(cl.values(), RF(0))
            print('B=%d n-N=%+d total: %s' % (B, n - N, tot.factor_str()))
            if len(cl) > 1:
                for k, v in sorted(cl.items()):
                    print('      class (down,zero,up)=%s: %s' % (k, v.factor_str()))
        print()
