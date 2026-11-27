"""AUDIT-4 part 5 (own code): symbolic / interval checks of the remaining analytic steps of R10."""
import sympy as sp
from flint import arb, ctx
ctx.prec = 200
m, n, u, v, a = sp.symbols('m n u v a', positive=True)
# (B1) monotonicity polynomials
print("(2m+1)^2(m+5/4) - (2m+2)^2(m+1/4) =", sp.expand((2*m+1)**2*(m+sp.Rational(5,4)) - (2*m+2)**2*(m+sp.Rational(1,4))), "  (>0: w^2(m+1/4) increasing)")
print("(2m+1)^2(m+3/2) - (2m+2)^2(m+1/2) =", sp.expand((2*m+1)**2*(m+sp.Rational(3,2)) - (2*m+2)**2*(m+sp.Rational(1,2))), "  (<0: w^2(m+1/2) decreasing)")
# Wallis limit: w_m^2 * m -> 1/pi
w2 = (sp.gamma(2*m+1)/sp.gamma(m+1)**2/4**m)**2
print("lim m*w_m^2 =", sp.limit(sp.simplify(sp.gammasimp(w2*m)), m, sp.oo))
# (2): constants
print("even: sqrt(pi(2m+1/2))/(pi(m+1/4)) - 2/sqrt(pi(2m+1/2)) =", sp.simplify(sp.sqrt(sp.pi*(2*m+sp.Rational(1,2)))/(sp.pi*(m+sp.Rational(1,4))) - 2/sp.sqrt(sp.pi*(2*m+sp.Rational(1,2)))))
odd = sp.sqrt(sp.pi*(2*m+sp.Rational(3,2)))/(sp.pi*(m+sp.Rational(5,4)))
nn = 2*m+1
print("odd : sqrt(pi(2m+3/2))/(pi(m+5/4)) - 2sqrt(n+1/2)/(sqrt(pi)(n+3/2)) at n=2m+1 =", sp.simplify(odd - 2*sp.sqrt(nn+sp.Rational(1,2))/(sp.sqrt(sp.pi)*(nn+sp.Rational(3,2)))))
t, j = sp.symbols('t j', positive=True)
print("even exponent: 2t^2/(m+t) - u^2/(n+u) at u=2t,n=2m =", sp.simplify(2*t**2/(m+t) - (2*t)**2/(2*m+2*t)))
print("odd exponent: 2j(j+1)/(m+1+j) - (u^2-1)/(n+u) at u=2j+1,n=2m+1 =", sp.simplify(2*j*(j+1)/(m+1+j) - ((2*j+1)**2-1)/(2*m+1+2*j+1)))
# Step 2 prefactor: 2n(n-1) - (n-u)(2n-1) = 2un - n - u
print("2n(n-1) - (n-u)(2n-1) =", sp.expand(2*n*(n-1) - (n-u)*(2*n-1)), " = n(2u-1)-u > 0 for 1<=u<n")
# chain
print("(u^2-1)/(n+u) - (v^2-2v)/(n+v-1) at u=v-1:", sp.simplify(((v-1)**2-1)/(n+v-1) - (v**2-2*v)/(n+v-1)))
d1 = sp.together((v**2-2*v)/(n+v-1) - (v**2/(n+v) - 2*v/n))
print("(v^2-2v)/(n+v-1) - [v^2/(n+v) - 2v/n] =", sp.factor(d1), " (numerator v(v^2+ ... ) -> must be >=0 for v>=1)")
d2 = sp.together(v**2/(n+v) - (v**2/n - v**3/n**2))
print("v^2/(n+v) - [v^2/n - v^3/n^2] =", sp.factor(d2))
# max of v exp(-v^2/n)
vv = sp.symbols('vv', positive=True)
crit = sp.solve(sp.diff(vv*sp.exp(-vv**2/n), vv), vv)
print("critical point of v e^{-v^2/n}:", crit, " value:", sp.simplify((vv*sp.exp(-vv**2/n)).subs(vv, crit[0])))
# derivative of log phi
ur = sp.symbols('ur', positive=True)
lp = sp.log(ur+1) - (ur**2-1)/(n+ur)
print("(log phi)' - [1/(u+1) - (u^2+2un+1)/(n+u)^2] =", sp.simplify(sp.diff(lp, ur) - (1/(ur+1) - (ur**2+2*ur*n+1)/(n+ur)**2)))
print("(u+1)(u^2+2un+1) - 2u^2 n =", sp.expand((ur+1)*(ur**2+2*ur*n+1) - 2*ur**2*n), " (>0)")
# threshold of (1.2 s - 1)(sqrt2 s - 1) >= s^2, s = sqrt n
s = sp.symbols('s', positive=True)
poly = sp.expand((sp.Rational(6,5)*s-1)*(sp.sqrt(2)*s-1) - s**2)
r = sp.solve(poly, s)
print("(1.2 s-1)(sqrt2 s-1)-s^2 =", poly, "; roots s =", [sp.N(x, 12) for x in r], "; n-threshold =", [sp.N(x**2, 12) for x in r])
print("coefficients:", sp.N(sp.Rational(6,5)*sp.sqrt(2)-1, 10), sp.N(sp.Rational(6,5)+sp.sqrt(2), 10), " (R10 prints 0.697 and 2.615: rounded in the safe direction)")
for N in (11, 12, 13):
    print("  n=%d: exact (1.2 sqrt n -1)(sqrt(2n)-1) - n = %s ; R10's 0.697n-2.615sqrt(n)+1 = %s" % (N, sp.N(poly.subs(s, sp.sqrt(N)), 8), sp.N(0.697*N-2.615*sp.sqrt(N)+1, 8)))
# brute interval check: (log phi)'(u) < 0 on u >= 1.2 sqrt n - 1, via the exact polynomial P(u) = (u+1)(u^2+2un+1) - (n+u)^2 > 0
bad = []
for N in range(3, 20001):
    sq = arb(N).sqrt()
    u0 = arb("1.2")*sq - 1
    P0 = (u0+1)*(u0*u0+2*u0*N+1) - (N+u0)**2
    # P'(u) = 3u^2 + 4un + 2u... check P increasing for u>=u0: P'(u)= 3u^2+2u(2n+1)+ (2n+1) - 2(n+u) = 3u^2 + 4un + 1 > 0
    if not bool(P0 > 0): bad.append(N)
print("n in 3..20000 where P(u0)>0 fails at u0=1.2sqrt(n)-1 (then phi need not decrease from u0):", bad, " [P'(u)=3u^2+4un+1>0, so P>0 on [u0,inf) otherwise]")
print("P'(u) =", sp.expand(sp.diff((ur+1)*(ur**2+2*ur*n+1) - (n+ur)**2, ur)))
# Step 3 constants
c = 4/(2*arb.pi()*arb(1).exp()).sqrt(); L = (1/c).log(); n0 = (arb("4.128")/L)**2
print("4/sqrt(2 pi e) =", c.str(20), " ln(sqrt(2 pi e)/4) =", L.str(20), " (4.128/L)^2 =", n0.str(20))
for N in (15990, 15991, 16100):
    print("  bound at n=%d: %s  <1: %s" % (N, (c*(arb("4.128")/arb(N).sqrt()).exp()).str(20), bool(c*(arb("4.128")/arb(N).sqrt()).exp() < 1)))
# is the case split exhaustive for integer u?  v<=1.2 sqrt n  or  u >= 1.2 sqrt n - 1  <=> v >= 1.2 sqrt n : yes (complementary, overlapping at equality)
