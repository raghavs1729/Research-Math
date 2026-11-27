#!/usr/bin/env python3
"""
x1.py -- constants for the lower limit 0 versus 2 (HKMN write the integral from 0, Johnston from 2).

  li(x) = PV int_0^x dt/log t,      int_0^x li(t) dt = x li(x) - li(x^2)     (x > 0)
  g(x)  := -int_0^x li(t) dt = li(x^2) - x li(x);   g'(x) = -li(x);  g(0+) = 0;
           g increases on (0, mu), decreases on (mu, oo), mu = 1.4513... the zero of li.
  c_0   := g(2) = li(4) - 2 li(2) = int_0^2 (0 - li(t)) dt  > 0.
  x_0   := the unique zero of g in (mu, oo).

For q >= 2 and a pair whose first prime p_1 = a (mod q) is >= 3:
  A_0(x) := int_0^x (phi(q) pi(t;q,a) - li(t)) dt = g(x) on (0, p_1],  and = c_0 + int_2^x(...) for x >= 2.
For q = 1:  A_0(x) = g(x) + sum_{p<=x} (x - p).

All roots are certified by sign changes in Arb (ball arithmetic), 20 digits.
"""
from flint import arb, ctx
ctx.prec = 200

def li(x):
    return arb(x).li()

def g(x):
    x = arb(x)
    return (x * x).li() - x * x.li()

def root(f, lo, hi, decreasing=True, bits=80):
    lo, hi = arb(lo), arb(hi)
    flo, fhi = f(lo), f(hi)
    assert (flo > 0 and fhi < 0) if decreasing else (flo < 0 and fhi > 0)
    for _ in range(bits):
        mid = (lo + hi) / 2
        fm = f(mid)
        if fm > 0:
            (lo, hi) = (mid, hi) if decreasing else (lo, mid)
        elif fm < 0:
            (lo, hi) = (lo, mid) if decreasing else (mid, hi)
        else:
            break
    return lo, hi

print("li(2)        =", li(2).str(30))
print("li(2) - 1    =", (li(2) - 1).str(30))
print("c_0 = li(4) - 2 li(2) = g(2) =", g(2).str(30))
mu = root(lambda x: x.li(), 1.2, 2, decreasing=False)
print("mu (zero of li) in [%s, %s]" % (mu[0].str(25, radius=False), mu[1].str(25, radius=False)))
x0 = root(g, 2, 3)
print("x_0 (zero of g = -int_0^x li) in [%s, %s]" % (x0[0].str(25, radius=False), x0[1].str(25, radius=False)))
print("g(3) =", g(3).str(20), "  c_0 - li(2)*(3-2) =", (g(2) - li(2)).str(20))

# q = 1:  A_0(x) = g(x) + pi(x) x - sum_{p<=x} p  between primes; decreasing on [2, 30] since pi(x) < li(x) there
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
print("check pi(p) < li(p) at primes <= 31 (so A_0' = pi - li < 0 on [2, 31]):",
      all(bool(li(p) > i + 1) for i, p in enumerate(primes)))
def A0_q1(x):
    x = arb(x)
    k = sum(1 for p in primes if p <= float(x.mid()))
    s = sum(p for p in primes if p <= float(x.mid()))
    return g(x) + k * x - s
for i in range(len(primes) - 1):
    p, pn = primes[i], primes[i + 1]
    vl, vr = A0_q1(p), g(pn) + (i + 1) * arb(pn) - sum(primes[:i + 1])
    print(f"  q=1: A_0({p}) = {vl.str(12, radius=False)},  A_0({pn}-) = {vr.str(12, radius=False)}")
    if vl > 0 and vr < 0:
        k, s = i + 1, sum(primes[:i + 1])
        r = root(lambda x: g(x) + k * x - s, p, pn)
        print("  x_1(1,1) in [%s, %s]" % (r[0].str(25, radius=False), r[1].str(25, radius=False)))
# threshold beyond which c_0 - (li(2)-1)(x-2) < 0
print("c_0/(li(2)-1) + 2 =", (g(2) / (li(2) - 1) + 2).str(12), " (for q = 1, A_0(x) <= c_0 - (li(2)-1)(x-2) once W_1^{(2,1)} <= 0)")
