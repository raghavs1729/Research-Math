"""V5 (exact, sympy).  The final coefficient extraction of the proof of Theorem U2.
With  Y/N = tau0 + tau1 e + tau2 e^2 + tau3 e^3 + ...,  tau2 = 2 pi^2 tau0,  tau3 = (4 pi^2/3) tau1,
      N(e) N(-e) = (e cot(pi e))^4 = pi^-4 (1 - (4 pi^2/3) e^2 + O(e^4)),
      zeta(3,1/2+e) = 7 z3 - (pi^4/2) e + 186 z5 e^2 - (2 pi^6/3) e^3 + O(e^4):
  [e^1] N(e) Y(-e) = -tau1/pi^4,   [e^3] 2 zeta(3,1/2+e) N(e) Y(-e) = -2 (pi^6 tau0 + 186 z5 tau1)/pi^4,
hence 16 [e^1](N Y(-e)) / [e^3](Xi Y(-e)) = 8 tau1/(pi^6 tau0 + 186 z5 tau1).
Also the expansions of the three periodic functions at v = 1/2 - u -> 0 used in the bilateral-sum argument."""
import sympy as sp
e, t0, t1, z3, z5, v = sp.symbols('e tau0 tau1 zeta3 zeta5 v')
pi = sp.pi
NN = sp.series((e * sp.cot(pi * e))**4, e, 0, 4).removeO()
assert sp.expand(NN - (1 - sp.Rational(4, 3) * pi**2 * e**2) / pi**4) == 0
# Hurwitz zeta at 1/2: zeta(s,1/2) = (2^s - 1) zeta(s); Taylor: zeta(3,1/2+e) = sum_m binom(-3,m) zeta(3+m,1/2) e^m
hz = 7 * z3 - 3 * 15 * sp.zeta(4) * e + 6 * 31 * z5 * e**2 - 10 * 63 * sp.zeta(6) * e**3
assert sp.expand(hz - (7 * z3 - pi**4 / 2 * e + 186 * z5 * e**2 - 2 * pi**6 / 3 * e**3)) == 0
t2 = 2 * pi**2 * t0; t3 = sp.Rational(4, 3) * pi**2 * t1
YNm = t0 - t1 * e + t2 * e**2 - t3 * e**3
NY = sp.expand(NN * YNm)
c1 = NY.coeff(e, 1)
c3 = sp.expand(2 * hz * NN * YNm).coeff(e, 3)
assert sp.simplify(c1 + t1 / pi**4) == 0
assert sp.simplify(c3 + 2 * (pi**6 * t0 + 186 * z5 * t1) / pi**4) == 0
assert sp.simplify(16 * c1 / c3 - 8 * t1 / (pi**6 * t0 + 186 * z5 * t1)) == 0
print("coefficient extraction: PASS  (the zeta(3) terms cancel)")
# periodic functions near u = 1/2, v = 1/2 - u
S2 = sp.series(pi**2 / sp.sin(pi * v)**2, v, 0, 1).removeO()
S3 = sp.series(pi**3 * sp.cos(pi * v) / sp.sin(pi * v)**3, v, 0, 1).removeO()
S4 = sp.series(pi**4 / sp.sin(pi * v)**4, v, 0, 1).removeO()
cot4 = sp.series(sp.cot(pi * v)**4, v, 0, 0).removeO()
print("pi^2 sec^2(pi u) =", S2, "+ O(v);  pi^3 sec^2 tan =", S3, "+ O(v);  pi^4 sec^4 =", S4, "+ O(v);  cot^4(pi v) =", cot4, "+ O(1)")
assert sp.expand(S3 - 1 / v**3) == 0
assert sp.expand(S4 - (1 / v**4 + 2 * pi**2 / (3 * v**2) + sp.Rational(11, 45) * pi**4)) == 0
assert sp.expand(S2 - (1 / v**2 + pi**2 / 3)) == 0
assert sp.expand(cot4 - (1 / (pi * v)**4 - sp.Rational(4, 3) / (pi * v)**2)) == 0
print("V5 PASS")
