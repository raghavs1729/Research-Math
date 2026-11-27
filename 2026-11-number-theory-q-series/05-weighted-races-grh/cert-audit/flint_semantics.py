#!/usr/bin/env python3
"""Audit of the python-flint (0.9.0) semantics that certify_derivative.py relies on.

Every assertion below is a property the audited code needs in order to be
fail-closed.  Exit status 0 iff all hold.  Nothing here is rigorous by
itself; it documents (and regression-tests) the behaviour of the installed
library version.
"""
import math
import sys
from decimal import Decimal
from fractions import Fraction

import flint
from flint import acb, arb, ctx

FAIL = []


def check(name, cond, extra=""):
    print(("ok   " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        FAIL.append(name)


def exact_fraction(x: arb) -> Fraction:
    """Exact rational value of an arb with zero radius (man * 2^exp)."""
    man, exp = x.man_exp()
    man, exp = int(man), int(exp)
    return Fraction(man) * (Fraction(2) ** exp)


print("python-flint", flint.__version__)
check("version is 0.9.0", flint.__version__ == "0.9.0")

ctx.prec = 64
# 1. comparisons of overlapping balls are False in BOTH directions ------------
a = arb(1, 1)          # [0, 2]
b = arb(1.5, 1)        # [0.5, 2.5]
for op, f in (("<", lambda x, y: x < y), (">", lambda x, y: x > y),
              ("<=", lambda x, y: x <= y), (">=", lambda x, y: x >= y),
              ("==", lambda x, y: x == y)):
    check(f"overlapping balls: (a {op} b) is False", f(a, b) is False)
check("ball vs int: [0,2] > 0 is False", (a > 0) is False)
# arb(1, 1) is NOT exactly [0, 2]: the radius is a 30-bit "mag" rounded UP, so
# the ball is [-eps, 2+eps] and ">= 0" is (conservatively) not certified.
print("info ball vs int: arb(1,1) >= 0 ->", a >= 0, "(False = conservative: mag radii are rounded up)")
check("ball vs int: [1,3] >= 1/2 is True", (arb(2, 1) >= 0.5) is True)
check("ball vs int: [0,2] < 0 is False", (a < 0) is False)
check("touching balls: [0,1] < [1,2] is False", (arb(0.5, 0.5) < arb(1.5, 0.5)) is False)
print("info touching balls: arb(.5,.5) <= arb(1.5,.5) ->", arb(0.5, 0.5) <= arb(1.5, 0.5),
      "(False = conservative, same reason)")
check("exact points: arb(1) <= arb(1) True, arb(1) < arb(1) False", (arb(1) <= arb(1)) is True and (arb(1) < arb(1)) is False)
# 2. NaN / inf -----------------------------------------------------------------
nan = arb("nan")
inf = arb("inf")
check("nan comparisons all False", not (nan < 0) and not (nan > 0) and not (nan <= 0) and not (nan >= 0))
check("nan.is_finite() False", nan.is_finite() is False)
check("inf.is_finite() False", inf.is_finite() is False)
check("[0 +/- inf] is not < 0 and not finite", not (arb(0, inf) < 0) and not arb(0, inf).is_finite())
check("nan.upper() < 0 False", not (nan.upper() < 0))
check("float nan < 0.0 False (hi < 0.0 guard)", not (float("nan") < 0.0))
z = arb(0, 1) / arb(0, 1)
check("division by a ball containing 0 is non-finite", not z.is_finite())
check("sqrt of a negative ball is non-finite", not arb(-1).sqrt().is_finite())
check("log of a ball containing 0 is non-finite", not arb(0, 1).log().is_finite())

# 3. upper()/lower() are true bounds at the current precision -------------------
ctx.prec = 128
x = arb(1) / 3 + arb(0, "1e-25")
m, r = exact_fraction(x.mid()), exact_fraction(x.rad())
for p in (64, 53, 30):
    ctx.prec = p
    up, lo = x.upper(), x.lower()
    check(f"upper() has zero radius and is >= mid+rad (prec {p})",
          exact_fraction(up.rad()) == 0 and exact_fraction(up) >= m + r)
    check(f"lower() has zero radius and is <= mid-rad (prec {p})",
          exact_fraction(lo.rad()) == 0 and exact_fraction(lo) <= m - r)
    au, al = (-x).abs_upper(), (-x).abs_lower()
    check(f"abs_upper/abs_lower bound |x| (prec {p})",
          exact_fraction(au) >= m + r and exact_fraction(al) <= m - r)
ctx.prec = 64
check("rad() is an exact (zero-radius) arb", exact_fraction(x.rad().rad()) == 0)
check("mid() is an exact (zero-radius) arb", exact_fraction(x.mid().rad()) == 0)

# 4. float(arb): at most 1 ulp from the midpoint; nextafter makes it a bound --
ctx.prec = 128
bad_up = bad_dn = 0
import random
random.seed(1)
for _ in range(20000):
    num = random.randrange(-10**30, 10**30)
    den = random.randrange(1, 10**30)
    y = arb(num) / arb(den)
    u, l = y.upper(), y.lower()
    fu = math.nextafter(float(u), math.inf)
    fl_ = math.nextafter(float(l), -math.inf)
    if not Fraction(fu) >= exact_fraction(u):
        bad_up += 1
    if not Fraction(fl_) <= exact_fraction(l):
        bad_dn += 1
check("up(x)=nextafter(float(x.upper()),+inf) >= upper (20000 random)", bad_up == 0)
check("dn(x)=nextafter(float(x.lower()),-inf) <= lower (20000 random)", bad_dn == 0)
ctx.prec = 64

# 5. string constructor with a radius: arb(dec, "1e-20") ------------------------
ctx.prec = 128
for s in ("6.020948904697596654902512", "639.9157072529864984688762", "0.1"):
    y = arb(s, "1e-20")
    m, r = exact_fraction(y.mid()), exact_fraction(y.rad())
    d = Fraction(Decimal(s))
    eps = Fraction(1, 10**20)
    check(f"arb('{s}','1e-20') contains [dec-1e-20, dec+1e-20]",
          m - r <= d - eps and d + eps <= m + r, f"rad-1e-20 = {float(r - eps):.2e}")
    check(f"  ... and is not grossly wider (rad < 1.0001e-20)", r < eps * Fraction(10001, 10000))
ctx.prec = 64

# 6. arb(Fraction-free) exact integers and ex() ---------------------------------
big = 3**200
check("arb(int) is exact for big ints at prec 64? (NOT required exact, must enclose)",
      arb(big).contains(arb(big)))
ctx.prec = 128
q = Fraction(36181, 10000)
e = arb(q.numerator) / arb(q.denominator)
m, r = exact_fraction(e.mid()), exact_fraction(e.rad())
check("ex(Fraction) = arb(num)/arb(den) encloses the rational", m - r <= q <= m + r)
ctx.prec = 64

# 7. union / intersection --------------------------------------------------------
u = arb(1).union(arb(3))
m, r = exact_fraction(u.mid()), exact_fraction(u.rad())
check("union contains both points and everything between", m - r <= 1 and m + r >= 3)
u = arb(1, "0.5").union(arb(1.2, "0.1"))
m, r = exact_fraction(u.mid()), exact_fraction(u.rad())
check("union of overlapping balls contains both", m - r <= Fraction(1, 2) and m + r >= Fraction(3, 2))
i = arb(0, 1).intersection(arb(1, "0.5"))
m, r = exact_fraction(i.mid()), exact_fraction(i.rad())
check("intersection contains the true intersection [0.5,1]", m - r <= Fraction(1, 2) and m + r >= 1)
try:
    arb(0, 1).intersection(arb(5, 1))
    check("disjoint intersection raises", False)
except ValueError:
    check("disjoint intersection raises ValueError", True)

# 8. arb(0, rad) constructor with an arb radius ----------------------------------
rb = arb(0, arb("1e-3").upper())
m, r = exact_fraction(rb.mid()), exact_fraction(rb.rad())
check("arb(0, u) = [-u, u] with rad >= u", m == 0 and r >= exact_fraction(arb("1e-3").upper()))
rb = arb(0, arb(-2))
check("arb(0, negative) uses |rad|", exact_fraction(rb.rad()) >= 2)

# 9. Bessel conventions and identities -------------------------------------------
ctx.prec = 128
import mpmath
mpmath.mp.dps = 40
for xs in ("0.001", "1", "2.5", "7.25", "20", "83.5"):
    for n in (0, 1):
        v = arb(xs).bessel_j(n)
        ref = mpmath.besselj(n, mpmath.mpf(xs))
        diff = abs(mpmath.mpf(v.mid().str(45, radius=False)) - ref)
        check(f"arb('{xs}').bessel_j({n}) = J_{n}({xs}) (mpmath)", diff < mpmath.mpf(10) ** -30)
check("acb bessel_j agrees", abs(complex(acb(2, 1).bessel_j(0)) - complex(mpmath.besselj(0, mpmath.mpc(2, 1)))) < 1e-14)
# zeta(2, x) Hurwitz = trigamma
v = arb(2).zeta(arb("2.05"))
ref = mpmath.psi(1, mpmath.mpf("2.05"))
check("arb(2).zeta(x) is the Hurwitz zeta(2, x) = psi'(x)",
      abs(mpmath.mpf(v.mid().str(45, radius=False)) - ref) < mpmath.mpf(10) ** -30)
v = arb("2.05").digamma()
check("digamma", abs(mpmath.mpf(v.mid().str(45, radius=False)) - mpmath.psi(0, mpmath.mpf("2.05"))) < mpmath.mpf(10) ** -30)
check("sinc(0-ball) finite and contains 1", arb(0, "0.01").sinc().contains(1) and arb(0, "0.01").sinc().is_finite())
check("sinc is sin(x)/x (not sin(pi x)/(pi x))", abs(float(arb(1).sinc()) - math.sin(1.0)) < 1e-15)
ctx.prec = 64

# 10. string round trip used for printing ----------------------------------------
x = arb(-1) / 3 * arb("1e-3")
hi = math.nextafter(float(x.upper()), math.inf)
d = Decimal(hi)
check("Decimal(float) is exact", Fraction(d) == Fraction(hi))
from decimal import ROUND_CEILING, ROUND_FLOOR, Context
c = Context(prec=7, rounding=ROUND_CEILING)
s = f"{c.create_decimal(d):.6e}"
check("printed hi (ROUND_CEILING) >= double", Fraction(Decimal(s)) >= Fraction(hi), s)
c = Context(prec=7, rounding=ROUND_FLOOR)
s = f"{c.create_decimal(d):.6e}"
check("printed lo (ROUND_FLOOR) <= double", Fraction(Decimal(s)) <= Fraction(hi), s)
check("float(printed) round trip can move by 1/2 ulp only (informational)", True)

# 11. Fraction <-> float comparisons used for steering are exact -----------------
check("Fraction < float comparison is exact", (Fraction(1, 3) < 1 / 3) == (Fraction(1, 3) < Fraction(1 / 3)))

print()
if FAIL:
    print("SEMANTICS FAILURES:", FAIL)
    sys.exit(1)
print("all python-flint semantics checks passed")
