#!/usr/bin/env python3
"""
Uniform-in-the-parameter bounds used in Theorems 6, 8, 9, 10 (RESULTS.md, "Uniformity lemma" and the
four proofs).  Every inequality is a statement  E >= 0  for a rational function E on a set
     n in {n0, n0+1, ...},  k (or x) >= x0,  parameter(s) in [-1/10, 1/10].
It is checked EXACTLY as follows: substitute n = n0 + N, x = x0 + X (N, X >= 0) and each parameter
p = (-1/10 + (1/10) T)/(1 + T), T >= 0 (a bijection of [0,oo) onto [-1/10, 1/10)); the endpoint p = 1/10 is
covered by continuity (all denominators stay positive: they are checked to have positive coefficients too).
Then E = A/B with polynomials A, B in the new variables; we assert that all coefficients of A are >= 0 and
all coefficients of B are > 0 (or all <= 0 / < 0).  That is a sufficient condition for E >= 0 on the set.
"""
import sympy as sp
from sympy import Rational as Q_
n, k, x, s, a, b, c = sp.symbols('n k x s a b c')
N, X, T1, T2 = sp.symbols('N X T1 T2', nonnegative=True)
TEN = Q_(1, 10)


def par(T):
    return (-TEN + TEN*T)/(1 + T)


def nonneg(E, sub):
    """E >= 0 after the substitution `sub`, by coefficient positivity of numerator and denominator."""
    E = sp.together(E.subs(sub, simultaneous=True))
    A, B = sp.fraction(sp.cancel(E))
    gens = sorted(A.free_symbols | B.free_symbols, key=str)
    if not gens:
        return A/B >= 0
    ca = sp.Poly(A, *gens).coeffs()
    cb = sp.Poly(B, *gens).coeffs()
    return (all(v >= 0 for v in ca) and all(v > 0 for v in cb)) or (all(v <= 0 for v in ca) and all(v < 0 for v in cb))


def check(name, E, sub):
    ok = nonneg(E, sub)
    print(("OK   " if ok else "FAIL ") + name)
    assert ok, name


# ------------------------------------------------------------------ Theorem 9 (pair of Theorem 2, rows x = k + c)
sig2 = (2*x+1)*(2*x+n+3)*(3*x+n+1)*(3*x+n+2)*(3*x+n+3)/((2*x+n+1)*(2*x+2*n+3)*(3*x+2*n+3)*(3*x+2*n+4)*(3*x+2*n+5))
rho2 = (n+1)*(4*n+3)*(4*n+5)*(2*x+n+2)*(3*x+n+1)/(9*(2*n+1)*(2*x+n+1)*(2*x+2*n+3)*(3*x+2*n+3)*(3*x+2*n+4))
R2 = (3*x+2*n+2)*(6*x+8*n+5)/(9*(2*n+1)*(2*x+n+1))
assert sp.cancel(rho2 - 1 - (R2.subs(x, x+1)*sig2 - R2)) == 0
print("Theorem 9 (x >= -1/10):")
check(" sigma(n,x) > 0, n >= 0", sig2, {n: N, x: -TEN + X})
check(" sigma(n,x) <= (2x+1)/(2x+5), n >= 1", (2*x+1)/(2*x+5) - sig2, {n: 1 + N, x: -TEN + X})
check(" rho(n,c) > 0, n >= 0, |c| <= 1/10", rho2, {n: N, x: par(T1)})
check(" rho(n,c) <= 1/5, n >= 0, |c| <= 1/10", Q_(1, 5) - rho2, {n: N, x: par(T1)})
check(" R(n,x) > 0", R2, {n: N, x: -TEN + X})
check(" R(n,x) <= x + 2", x + 2 - R2, {n: N, x: -TEN + X})

# ------------------------------------------------------------------ Theorem 10 (pair of Theorem 4, rows x = k + b)
sig4 = (4*x+6*n+1)*(4*x+6*n+3)/(16*(x+n+1)*(x+3*n+1))
rho4 = 3*(6*n+1)*(6*n+5)*(4*x+6*n+1)*(4*x+6*n+3)*(4*x+6*n+5)/(1024*(n+1)*(x+n+1)*(x+3*n+1)*(x+3*n+2)*(x+3*n+3))
Q4 = (64*x**3+(160*n+144)*x**2+(240*n**2+432*n+92)*x+360*n**3+612*n**2+230*n+15)/(64*(x+3*n+1)*(x+3*n+2))
assert sp.cancel((n+1)*rho4 - n - (Q4.subs(x, x+1)*sig4 - Q4)) == 0
print("Theorem 10 (x >= -1/10):")
check(" sigma(n,x) > 0", sig4, {n: N, x: -TEN + X})
check(" 1 - sigma(n,x) >= (3/4)(n+1)/(x+3n+1)", 1 - sig4 - Q_(3, 4)*(n+1)/(x+3*n+1), {n: N, x: -TEN + X})
check(" (3/4)(n+1) <= x+3n+1 (factors in [0,1))", x + 3*n + 1 - Q_(3, 4)*(n+1), {n: N, x: -TEN + X})
check(" rho(n,b) > 0, |b| <= 1/10", rho4, {n: N, x: par(T1)})
check(" rho(n,b) <= 9/10, n >= 0, |b| <= 1/10", Q_(9, 10) - rho4, {n: N, x: par(T1)})
check(" Q(n,x) >= 0", Q4, {n: N, x: -TEN + X})
check(" Q(n,x) <= x + 3n + 3", x + 3*n + 3 - Q4, {n: N, x: -TEN + X})
assert sp.limit(Q4.subs(n, 0)/x, x, sp.oo) == 1
print(" Q(0,x)/x -> 1: OK")

# ------------------------------------------------------------------ Theorem 8 (seed parameter s)
Ns = ((128*n + (128*s + 64))*k**3 + (320*n**2 + (384*s + 448)*n + (-192*s**2 + 192*s + 144))*k**2
      + (480*n**3 + (1056*s + 1104)*n**2 + (-352*s**2 + 736*s + 616)*n + (96*s**3 - 240*s**2 + 40*s + 92))*k
      + (720*n**4 + (1440*s + 1584)*n**3 + (-192*s**2 + 1440*s + 1072)*n**2 + (96*s**3 - 272*s**2 + 296*s + 260)*n
         + (-16*s**4 + 64*s**3 - 56*s**2 - 16*s + 15)))
Qs = Ns/(64*(2*n + 2*s + 1)*(k + 3*n + 1)*(k + 3*n + 2))
sigs = (4*k + 6*n - 2*s + 1)*(4*k + 6*n - 2*s + 3)/(16*(k + 3*n + 1)*(k + n - s + 1))
rhos = ((6*n + 2*s + 1)*(6*n + 2*s + 3)*(6*n + 2*s + 5)*(4*k + 6*n - 2*s + 1)*(4*k + 6*n - 2*s + 3)*(4*k + 6*n - 2*s + 5)
        / (1024*(n + 1)*(2*n + 2*s + 1)*(k + 3*n + 1)*(k + 3*n + 2)*(k + 3*n + 3)*(k + n - s + 1)))
assert sp.cancel((n + 1)*rhos - n - (Qs.subs(k, k + 1)*sigs - Qs)) == 0
print("Theorem 8 (|s| <= 1/10, k >= 0):")
S = {n: N, k: X, s: par(T1)}
check(" sigma_s > 0", sigs, S)
check(" 1 - sigma_s >= (7/10)(n+1)/(k+3n+1)", 1 - sigs - Q_(7, 10)*(n+1)/(k+3*n+1), S)
check(" rho_s(n,0) > 0", rhos.subs(k, 0), S)
check(" rho_s(n,0) <= 9/10", Q_(9, 10) - rhos.subs(k, 0), S)
check(" Q_s >= 0", Qs, S)
check(" Q_s <= 2(k+3n+3)", 2*(k + 3*n + 3) - Qs, S)
assert sp.simplify(sp.limit(Qs.subs(n, 0)/k, k, sp.oo) - 1) == 0
print(" Q_s(0,k)/k -> 1 for every s: OK")

# ------------------------------------------------------------------ Theorem 6 (parameters a, b)
u = k + 3*n + 2*a + b
sig6 = -Q_(1, 3)*(u+1)*(k+b-a+Q_(1, 2))/((u+Q_(3, 2))*(k+n+b+1))
rho6 = (-Q_(1, 3)*Q_(3, 4)*(u+1)*(u+2)*(u+3)/((u+Q_(3, 2))*(u+Q_(5, 2))*(u+Q_(7, 2)))
        * (2*n+2*a+Q_(3, 2))*(2*n+2*a+Q_(5, 2))*(n+a+1)**2/((2*n+2*a+1)*(2*n+2*a+2))
        * (-Q_(1, 2)-a-n)/((k+n+b+1)*(n+a+Q_(3, 4))*(n+a+Q_(5, 4))))
R6 = 3*(u+2)**2/(4*(u+2)**2 - 1)
assert sp.cancel(rho6 - 1 - (R6.subs(k, k+1)*sig6 - R6)) == 0
assert sp.cancel(rho6*sig6.subs(n, n+1) - sig6*rho6.subs(k, k+1)) == 0
print("Theorem 6 (|a|,|b| <= 1/10, k >= 0, n >= 0): WZ with (a,b) and compatibility exact: OK")
S = {n: N, k: X, a: par(T1), b: par(T2)}
check(" -sigma > 0", -sig6, S)
check(" |sigma| <= 1/3", Q_(1, 3) + sig6, S)
check(" rho(n,0) > 0", rho6.subs(k, 0), S)
check(" rho(n,0) <= 1/2", Q_(1, 2) - rho6.subs(k, 0), S)
check(" R > 0", R6, S)
check(" R <= 1", 1 - R6, S)
print("ALL CHECKS PASSED")
