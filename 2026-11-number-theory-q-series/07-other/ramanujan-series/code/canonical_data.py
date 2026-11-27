"""Canonical hypergeometric data of the eleven series (audit item R2): for the summand t(k),
t(k+1)/t(k) = z * prod (k + u_i) / prod (k + l_j) * P(k+1)/P(k); we print z, the multisets {u_i mod 1},
{l_j mod 1} after cancelling common classes, and the polynomial factor P.  Exact (sympy)."""
import sympy as sp
from collections import Counter
k = sp.symbols('k')
B = sp.binomial
T = {
 'T1':  (3*k+2)**2*16**k/((2*k+1)*(6*k+1)*(6*k+5)*B(6*k, 3*k)),
 'T2':  (8*k-3)*B(4*k, 2*k)/(k*(4*k-1)*9**k*B(2*k, k)**2),
 'T3':  (4*k-1)*3**k/((2*k-1)*k**2*B(2*k, k)*B(3*k, k)),
 'T4':  (360*k**3+612*k**2+230*k+15)*B(3*k, k)*B(6*k, 3*k)**2/((3*k+1)*(3*k+2)*2**(15*k)),
 'T5':  (77*k**3-86*k**2+29*k-3)*(-32)**k/(k**3*(3*k-1)*(3*k-2)*B(3*k, k)*B(4*k, k)**2),
 'T7':  B(2*k, k)*16**k*(1656*k**4-2604*k**3+1158*k**2-229*k+15)/(k**2*(2*k-1)**2*(6*k-1)*(6*k-5)*B(6*k, 3*k)*B(3*k, k)),
 'T11': (1344*k**3+944*k**2+156*k+9)*B(2*k, k)**2*B(4*k, 2*k)**3/((2*k+1)*2**(22*k)),
}
def frac1(x):
    r = x - sp.floor(x)
    return r if r != 0 else sp.Integer(1)
for name, t in T.items():
    r = sp.factor(sp.simplify(sp.combsimp((t.subs(k, k+1)/t).rewrite(sp.gamma))))
    num, den = sp.fraction(r)
    cnt = Counter(); z = sp.Integer(1); polys = []
    for poly, sign in ((num, 1), (den, -1)):
        c, facs = sp.factor_list(poly, k)
        z *= c**sign
        for f, m in facs:
            if sp.degree(f, k) == 1:
                a1, a0 = sp.Poly(f, k).all_coeffs()
                z *= a1**(sign*m)
                cnt[frac1(sp.Rational(a0, a1))] += sign*m
            else:
                lc = sp.Poly(f, k).LC(); z *= lc**(sign*m)
                polys.append((sign*m, sp.expand(f/lc)))
    up = sorted([(a, m) for a, m in cnt.items() if m > 0], key=lambda v: v[0])
    lo = sorted([(a, -m) for a, m in cnt.items() if m < 0], key=lambda v: v[0])
    f = lambda L: ' '.join(f"({a})^{m}" if m > 1 else f"({a})" for a, m in L)
    print(f"{name}: z = {z};  upper {f(up)}  /  lower {f(lo)}")
