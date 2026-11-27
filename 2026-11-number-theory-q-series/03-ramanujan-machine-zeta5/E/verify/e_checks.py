"""Exact checks for ../THEOREM.md (three continued fractions for e, Table 3 of arXiv:1907.00205v4).
For each fraction PCF[a_n, b_n]:
 (1) U_n = (n+1)!^2 Pi(n) solves u_n = a_n u_{n-1} + b_n u_{n-2}   (polynomial identity, sympy);
 (2) p_n/q_n = a_0 - U_0/U_{-1} + 1/(U_{-1}^2 T_n),  T_n = sum_{j=0}^n prod_{i<=j}(-b_i)/(U_j U_{j-1})   (exact, n <= 40);
 (3) the telescoping identity for the summand (sympy), hence the value of T = lim T_n;
 (4) the resulting closed form equals the stated one (sympy), and the convergents approach it (60 digits)."""
import sympy as sp
from fractions import Fraction as Fr
from math import factorial
from mpmath import mp, mpf, e as E_
n, j = sp.symbols('n j')
R = sp.Rational
cases = {
 'e-i':   dict(a=n * (n + 6) + 7, b=-n**2 * (n + 3), Pi=n**3 + 6 * n**2 + 11 * n + 8,
               # summand t_j = (j+2)(j+3)/(6 (j+1)! Pi(j) Pi(j-1));  (j+2)(j+3)/(Pi(j-1)Pi(j)) = (j+1)N(j)/Pi(j-1) - N(j+1)/Pi(j) + 1/3
               core=(j + 2) * (j + 3), D0=lambda x: (x - 1)**3 + 6 * (x - 1)**2 + 11 * (x - 1) + 8, N=-(j**2 + 3 * j + 1) / 3, Q=R(1, 3), sign=1, pref=R(1, 6),
               value=6 * sp.E / (2 * sp.E - 3)),
 'e-ii':  dict(a=(n + 1) * (n + 3), b=n**2 * (n + 3), Pi=(n + 3) * (n**2 + 7 * n + 11),
               # summand t_j = (-1)^j/(6 (j+1)! c_{j+1} c_{j+2}), c_m = m^2+3m+1;  1/(c_{j+1}c_{j+2}) = (j+1)N(j)/c_{j+1} + N(j+1)/c_{j+2} + 1/2
               core=sp.Integer(1), D0=lambda x: (x + 1)**2 + 3 * (x + 1) + 1, N=-(j + 3) / 2, Q=R(1, 2), sign=-1, pref=R(1, 6),
               value=3 / (4 * sp.E - 10)),
 'e-iii': dict(a=n * (n + 3) + 3, b=-n**2 * (n + 2), Pi=sp.Integer(1),
               # summand t_j = (j+2)/(2 (j+1)!)
               core=None, value=4 * sp.E / (2 * sp.E - 1)),
}
mp.dps = 60
for name, c in cases.items():
    a, b, Pi = c['a'], c['b'], c['Pi']
    # (1)  U_n/U_{n-2} etc: divide the recurrence by (n-1)!^2:  (n+1)^2 n^2 Pi(n) = a_n n^2 Pi(n-1) + b_n Pi(n-2)
    res = sp.expand((n + 1)**2 * n**2 * Pi - a * n**2 * Pi.subs(n, n - 1) - b * Pi.subs(n, n - 2))
    assert res == 0, (name, res)
    # (2)
    af = sp.lambdify(n, a); bf = sp.lambdify(n, b); Pf = sp.lambdify(n, Pi)
    U = lambda m: Fr(factorial(m + 1)**2 * int(Pf(m)))
    p = [Fr(1), Fr(int(af(0)))]; q = [Fr(0), Fr(1)]
    T = Fr(0); prod = Fr(1)
    for m in range(0, 41):
        if m >= 1:
            p.append(int(af(m)) * p[-1] + int(bf(m)) * p[-2]); q.append(int(af(m)) * q[-1] + int(bf(m)) * q[-2])
            prod *= -int(bf(m))
        T += prod / (U(m) * U(m - 1))
        assert p[-1] / q[-1] == int(af(0)) - U(0) / U(-1) + 1 / (U(-1)**2 * T), (name, m)
    # (3), (4)
    if c['core'] is not None:
        D0 = c['D0'](j); D1 = c['D0'](j + 1); N = c['N']; Q = c['Q']; s = c['sign']
        assert sp.expand(c['core'] - ((j + 1) * N * D1 - s * N.subs(j, j + 1) * D0 + Q * D0 * D1)) == 0, name
        # sum_j s^j core/((j+1)! D0 D1) = N(0)/D0(0) + Q * sum_j s^j/(j+1)!
        tail = (sp.E - 1) if s == 1 else (1 - 1 / sp.E)
        Tval = c['pref'] * (N.subs(j, 0) / D0.subs(j, 0) + Q * tail)
        # check the summand really is pref * s^j core/((j+1)! D0 D1): compare with prod(-b_i)/(U_j U_{j-1}) for j <= 12
        prod = Fr(1)
        for m in range(0, 13):
            if m >= 1: prod *= -int(bf(m))
            lhs = prod / (U(m) * U(m - 1))
            rhs = sp.nsimplify(c['pref'] * s**m * c['core'].subs(j, m) / (sp.factorial(m + 1) * D0.subs(j, m) * D1.subs(j, m)))
            assert sp.Rational(lhs.numerator, lhs.denominator) == rhs, (name, m)
    else:
        Tval = sp.E - R(1, 2)
        prod = Fr(1)
        for m in range(0, 13):
            if m >= 1: prod *= -int(bf(m))
            assert prod / (U(m) * U(m - 1)) == Fr(m + 2, 2 * factorial(m + 1))
    U0 = sp.Integer(int(U(0))); Um1 = sp.Integer(int(U(-1)))
    closed = sp.simplify(a.subs(n, 0) - U0 / Um1 + 1 / (Um1**2 * Tval))
    assert sp.simplify(closed - c['value']) == 0, (name, closed)
    num = mpf(p[-1].numerator) / p[-1].denominator / (mpf(q[-1].numerator) / q[-1].denominator)
    print("%-6s U_n = (n+1)!^2 * (%s);  T = %s;  PCF = %s;  |p_40/q_40 - value| = %s" % (name, sp.factor(Pi), sp.simplify(Tval), c['value'], mp.nstr(abs(num - mpf(sp.N(c['value'], 70))), 3)))
print("PASS")
