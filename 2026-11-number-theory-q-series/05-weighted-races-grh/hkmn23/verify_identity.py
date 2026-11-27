#!/usr/bin/env python3
"""
verify_identity.py -- independent numerical checks (mpmath, 30 digits; NOT part of a proof).

(1) Lemma A:  for x > 2,
      int_2^x (phi pi(t;q,a) - li(t)) dt = -li(2)(x-2) + W_1(x)/log x + int_2^x W_1(t) h''(t) dt,
      W_1(t) = phi * sum_{p<=t, p=a (q)} log p (t-p) - (t-2)^2/2,
      h''(t) = (x+t)/(t^2 log^2 t) + 2(x-t)/(t^2 log^3 t).
    Left side in closed form:  phi*sum_{p<=x,p=a}(x-p) - [x li(x) - li(x^2)] + [2 li(2) - li(4)].
    Right side: piecewise Gauss-Legendre quadrature between consecutive primes = a (q).
(2) The theorem itself on a grid: A(x) < 0 and A(x) <= -li(2)(x-2) + W_1(x)/log x  for p_1 <= x <= 2000 (p_1 = first prime = a mod q; before p_1 everything is explicit).
(3) HKMN's normalisation: A_0(x) = int_0^x (...) = A(x) + c_0 changes sign exactly once, at x_0 = 2.6111...
"""
import mpmath as mp
mp.mp.dps = 30

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\0\0"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

P = primes_upto(5000)
def phi_of(n):
    from math import gcd
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)

def W1(q, a, t):
    phi = phi_of(q)
    return phi * mp.fsum(mp.log(p) * (t - p) for p in P if p <= t and p % q == a % q) - (t - 2) ** 2 / 2

def A_closed(q, a, x):
    phi = phi_of(q)
    s = phi * mp.fsum((x - p) for p in P if p <= x and p % q == a % q)
    return s - (x * mp.li(x) - mp.li(x * x)) + (2 * mp.li(2) - mp.li(4))

def rhs(q, a, x):
    hpp = lambda t: (x + t) / (t ** 2 * mp.log(t) ** 2) + 2 * (x - t) / (t ** 2 * mp.log(t) ** 3)
    nodes = [mp.mpf(2)] + [mp.mpf(p) for p in P if 2 < p < x and p % q == a % q] + [mp.mpf(x)]
    integ = mp.fsum(mp.quad(lambda t: W1(q, a, t) * hpp(t), [nodes[i], nodes[i + 1]]) for i in range(len(nodes) - 1))
    return -mp.li(2) * (x - 2) + W1(q, a, x) / mp.log(x) + integ

print("(1) Lemma A, left - right:")
for (q, a) in [(2, 1), (5, 4), (20, 9), (24, 1), (1, 0), (7, 3)]:
    for x in (mp.mpf("50.5"), mp.mpf("200.25")):
        if q == 1:
            # q = 1: pi(t) = pi(t;2,1) + 1 for t >= 2;  identity for (1,·) holds with W_1 built from all primes
            l = A_closed(1, 0, x); r = rhs(1, 0, x)
        else:
            l = A_closed(q, a, x); r = rhs(q, a, x)
        print(f"   (q,a)=({q},{a})  x={x}:  left={mp.nstr(l, 15)}  right={mp.nstr(r, 15)}  diff={mp.nstr(l - r, 3)}")

print("(2) grid check p_1 <= x <= 2000 (step 0.37) for the 12 cases: max of A(x) and of A(x) - [-li(2)(x-2) + W_1(x)/log x]")
cases = [(2, 1), (3, 1), (4, 1), (5, 1), (5, 4), (8, 1), (12, 1), (15, 1), (15, 4), (20, 1), (20, 9), (24, 1)]
for (q, a) in cases:
    mA, mD, mW = -mp.inf, -mp.inf, -mp.inf
    p1 = min(p for p in P if p % q == a % q)
    x = mp.mpf(p1)
    while x <= 2000:
        A = A_closed(q, a, x); w = W1(q, a, x)
        mA = max(mA, A / (x - 2)); mD = max(mD, A - (-mp.li(2) * (x - 2) + w / mp.log(x))); mW = max(mW, w / x ** mp.mpf(1.5))
        x += mp.mpf("0.37")
    print(f"   ({q:2d},{a:2d}), x >= p_1 = {p1}: max A(x)/(x-2) = {mp.nstr(mA, 8)} (must be <= -li(2) = -1.0451...), max [A - bound] = {mp.nstr(mD, 5)} (must be <= 0), max W_1/x^1.5 = {mp.nstr(mW, 6)} (must be <= 0)")

print("(3) c_0 =", mp.nstr(mp.li(4) - 2 * mp.li(2), 20), "  root of li(x^2) - x li(x):", mp.nstr(mp.findroot(lambda x: mp.li(x * x) - x * mp.li(x), 2.6), 20))
print("    q=1: A(x)=int_2^x(pi-li) on a grid, max of A(x)/(x-2) for 2<x<=2000:",
      mp.nstr(max(A_closed(1, 0, mp.mpf(2) + mp.mpf(k) * mp.mpf("0.37")) / (mp.mpf(k) * mp.mpf("0.37")) for k in range(1, 5400)), 8),
      " (must be <= -(li(2)-1) = -0.04516...)")
