"""WP-C numerical check (N2-PLAN.md, sketch (C)).

1. cf(N, k0) * sqrt(N) -> sqrt(3) / (sqrt(2) * pi^(5/2)),  k0 = floor(x* N), x* = sqrt(3)/2.
2. The bookkeeping of the eleven factorials m = c N + u (signed sums of c and u, I1, I2, the constant).
3. The exact "main term" identity used in Centre.lean:
     log(2^12 (4N-1)(2k0+1)) + sum(+-S(c_i,u_i,N)) - 22 N log 2 + (1/2) log N
       = log((4N-1)(2k0+1)/N^2) - (5/2) log(2 pi),
   S(c,u,N) = cN log N + cN(log c - 1) + u log N + u log c + (1/2) log(2 pi c N).
4. The log identities with a = log 2, b = log(sqrt 3 - 1).
"""
from mpmath import mp, mpf, sqrt, pi, log, loggamma, floor, exp

mp.dps = 40
x = sqrt(3) / 2
LIMIT = sqrt(3) / (sqrt(2) * pi ** (mpf(5) / 2))


def lf(m):
    assert m >= 0 and m == int(m)
    return loggamma(mpf(m) + 1)


def log_cf(N, k):
    assert 4 * k + 4 - 2 * N >= 0 and N - k - 1 >= 0 and 2 * N - 2 >= 0
    return (12 * log(2) + log(4 * N - 1) + log(2 * k + 1) + lf(4 * N) + lf(8 * k + 4 * N)
            + lf(4 * k + 4 - 2 * N)
            - (22 * N * log(2) + lf(2 * N) + lf(2 * N - 2) + lf(4 * k + 2 * N) + lf(8 * k + 8 - 4 * N)
               + 2 * lf(N - k - 1) + 2 * lf(N + k)))


def S(c, u, N):
    return c * N * log(N) + c * N * (log(c) - 1) + u * log(N) + u * log(c) + log(2 * pi * c * N) / 2


def table(th):
    """(sign, c, u) for the eleven factorials (squares listed with sign +-2)."""
    return [
        (+1, mpf(4), mpf(0)),
        (+1, 8 * x + 4, -8 * th),
        (+1, 4 * x - 2, 4 - 4 * th),
        (-1, mpf(2), mpf(0)),
        (-1, mpf(2), mpf(-2)),
        (-1, 4 * x + 2, -4 * th),
        (-1, 8 * x - 4, 8 - 8 * th),
        (-2, 1 - x, th - 1),
        (-2, 1 + x, -th),
    ]


print("1. cf(N,k0) sqrt(N) / limit")
for N in [8, 9, 10, 50, 100, 500, 1000, 3000, 5000, 20000]:
    k0 = int(floor(x * N))
    val = exp(log_cf(N, k0) + log(N) / 2)
    print(f"   N={N:6d} k0={k0:6d}  ratio={mp.nstr(val / LIMIT, 12)}")

print("2. bookkeeping")
I2 = (8 * log(8 * x + 4) + 4 * log(4 * x - 2) - 4 * log(4 * x + 2) - 8 * log(8 * x - 4)
      + 2 * log(1 - x) - 2 * log(1 + x))
I1 = (8 * log(2) + 4 * log(8 * x + 4) - 2 * log(4 * x - 2) - 4 * log(2) - 2 * log(4 * x + 2)
      + 4 * log(8 * x - 4) - 2 * log(1 - x) - 2 * log(1 + x))
print("   I2 =", mp.nstr(I2, 5), "  I1 - 22 log 2 =", mp.nstr(I1 - 22 * log(2), 5))
cst = 4 * log(4 * x - 2) + 2 * log(2) - 8 * log(8 * x - 4) + 2 * log(1 - x)
print("   const - log(1/16384) =", mp.nstr(cst - log(mpf(1) / 16384), 5))
print("   2^15 x* 4 (2pi)^(-5/2)/16384 / limit =",
      mp.nstr(2 ** 15 * x * 4 * (2 * pi) ** (-mpf(5) / 2) / 16384 / LIMIT, 12))
worst = 0
for N in [8, 9, 17, 100, 1001, 4567]:
    k0 = int(floor(x * N))
    th = x * N - k0
    tb = table(th)
    ms = [4 * N, 8 * k0 + 4 * N, 4 * k0 + 4 - 2 * N, 2 * N, 2 * N - 2, 4 * k0 + 2 * N,
          8 * k0 + 8 - 4 * N, N - k0 - 1, N + k0]
    for (s, c, u), m in zip(tb, ms):
        worst = max(worst, abs(c * N + u - m))
    sc = sum(s * c for s, c, u in tb)
    su = sum(s * u for s, c, u in tb)
    sclc = sum(s * c * log(c) for s, c, u in tb)
    main = (12 * log(2) + log(4 * N - 1) + log(2 * k0 + 1) + sum(s * S(c, u, N) for s, c, u in tb)
            - 22 * N * log(2) + log(N) / 2)
    target = log(mpf(4 * N - 1) * (2 * k0 + 1) / N ** 2) - mpf(5) / 2 * log(2 * pi)
    err = log_cf(N, k0) + log(N) / 2 - main      # sum of the Stirling errors, -> 0
    print(f"   N={N:5d}: sum c={mp.nstr(sc, 3)} sum u={mp.nstr(su, 3)} "
          f"sum c log c - 22 log 2={mp.nstr(sclc - 22 * log(2), 3)} "
          f"main-target={mp.nstr(main - target, 3)} stirling err={mp.nstr(err, 5)}")
print("   max |cN+u-m| =", mp.nstr(worst, 3))

print("4. log identities, a = log 2, b = log(sqrt3 - 1)")
a, b = log(2), log(sqrt(3) - 1)
for name, lhs, rhs in [
    ("log(8x+4)", log(8 * x + 4), 3 * a - b),
    ("log(4x-2)", log(4 * x - 2), a + b),
    ("log(4x+2)", log(4 * x + 2), 2 * a - b),
    ("log(8x-4)", log(8 * x - 4), 2 * a + b),
    ("log(1-x)", log(1 - x), 2 * b - 2 * a),
    ("log(1+x)", log(1 + x), -2 * b),
]:
    print(f"   {name}: diff = {mp.nstr(lhs - rhs, 3)}")
print("   limit / (4 sqrt3 (2pi)^(-5/2)) =", mp.nstr(LIMIT / (4 * sqrt(3) * (2 * pi) ** (-mpf(5) / 2)), 12))
