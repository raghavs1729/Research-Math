"""U3: the final coefficient extraction (exact, sympy).  Notation of ../THEOREM.md, section 7.
  (Y/N)(e) = tau0 + tau1 e + 2 pi^2 tau0 e^2 + (4 pi^2/3) tau1 e^3 + ...            (U2, Prop. 8)
  (Z/N)(e) = z0 + z1 e + z2 e^2 + z3 e^3 + ...,  z0 = e0, z1 = pi e1, z2 = pi^2 (e0+e2), z3 = (4 pi^3/3) e1 - 16 pi^2
  N(e)N(-e) = (e cot pi e)^4,   N(0) = pi^-2,   zeta(3,1/2+e) = 7 z3 - (pi^4/2) e + 186 z5 e^2 - (2 pi^6/3) e^3 + ...
  Y*(u) = u^3/2 - Y(u)/16 + Z(u)/32.
Claims:  192 mu4* = -32 [e^1] N(e) Y*(-e) = e1/pi^3 - 2 tau1/pi^4 =: E1,
         g*(0) = -2 [e^3] Xi(e) Y*(-e) = (93/4) z5 E1 + (pi^2/16)(e0 + e2 - 4 tau0),
         Legendre (g(0) mu4* - mu4 g*(0) = 1/6) with e2 = e0 + 8  <=>  E1 = 8 (1 + S1 K)/(pi^4 S0),
         1/U3 = g*(0)/(192 mu4*) - 4 = (93/4) z5 - 4 + (pi^6/64) S0 K/(1 + S1 K),  K = 4 + pi^2 (Z0 - 2 S0),
 where tau0 = pi^2 S0, tau1 = 8 pi^2 S1, e0 = pi^2 Z0."""
import sympy as sp
e, t0, t1, e0, e1, e2, z3s, z5 = sp.symbols('e tau0 tau1 e0 e1 e2 zeta3 zeta5')
S0, S1, Z0 = sp.symbols('S0 S1 Z0')
pi = sp.pi; R = sp.Rational
NN = (1 - R(4, 3) * pi**2 * e**2) / pi**4                      # (e cot pi e)^4 mod e^4
hz = 7 * z3s - pi**4 / 2 * e + 186 * z5 * e**2 - 2 * pi**6 / 3 * e**3
YN = lambda x: t0 + t1 * x + 2 * pi**2 * t0 * x**2 + R(4, 3) * pi**2 * t1 * x**3
ZN = lambda x: e0 + pi * e1 * x + pi**2 * (e0 + e2) * x**2 + (R(4, 3) * pi**3 * e1 - 16 * pi**2) * x**3
# N(e) * Y*(-e) = N(e) (-e)^3/2 + N(e)N(-e) [ -(Y/N)(-e)/16 + (Z/N)(-e)/32 ];  N(e) = pi^-2 + O(e)
core = sp.expand(NN * (-YN(-e) / 16 + ZN(-e) / 32))
NYs = lambda j: core.coeff(e, j) + (-R(1, 2) / pi**2 if j == 3 else 0)           # [e^j] N(e)Y*(-e), j <= 3
assert sp.simplify(NYs(3)) == 0          # <nu, r*> = 0: consistency
E1 = e1 / pi**3 - 2 * t1 / pi**4
assert sp.simplify(-32 * NYs(1) - E1) == 0
# [e^3] 2 zeta(3,1/2+e) N(e) Y*(-e): the e^3-part of Y* contributes 2*7 z3 * N(0) * (-1/2)
c3 = sp.expand(2 * hz * core).coeff(e, 3) + 2 * 7 * z3s * (-R(1, 2)) / pi**2
gs0 = sp.expand(-2 * c3)
target = R(93, 4) * z5 * E1 + pi**2 / 16 * (e0 + e2 - 4 * t0)
assert sp.simplify(gs0 - target) == 0
print("192 mu4* = E1 = e1/pi^3 - 2 tau1/pi^4;  g*(0) = (93/4) zeta5 E1 + (pi^2/16)(e0 + e2 - 4 tau0);  zeta(3) cancels: PASS")
# U2 quantities
mu4 = t1 / (6 * pi**4)                                                     # N(0) D/6, tau1 = pi^2 D
g0 = 4 * (pi**6 * t0 + 186 * z5 * t1) / pi**4
sub = {t0: pi**2 * S0, t1: 8 * pi**2 * S1, e0: pi**2 * Z0, e2: pi**2 * Z0 + 8}
K = 4 + pi**2 * (Z0 - 2 * S0)
leg = sp.expand((g0 * E1 / 192 - mu4 * target).subs(sub))
sol = sp.solve(sp.Eq(leg, R(1, 6)), e1)
assert len(sol) == 1
E1v = sp.simplify(E1.subs(sub).subs(e1, sol[0]))
assert sp.simplify(E1v - 8 * (1 + S1 * K) / (pi**4 * S0)) == 0
assert sp.simplify(sol[0] - 8 * (1 + 4 * S1 + pi**2 * S1 * Z0) / (pi * S0)) == 0
final = sp.simplify((target / E1).subs(sub).subs(e1, sol[0]) - 4)
assert sp.simplify(final - (R(93, 4) * z5 - 4 + pi**6 / 64 * S0 * K / (1 + S1 * K))) == 0
print("Legendre => E1 = 8(1 + S1 K)/(pi^4 S0), e1 = 8(1 + 4 S1 + pi^2 S1 Z0)/(pi S0);  1/U3 = (93/4) zeta5 - 4 + (pi^6/64) S0 K/(1 + S1 K): PASS")
