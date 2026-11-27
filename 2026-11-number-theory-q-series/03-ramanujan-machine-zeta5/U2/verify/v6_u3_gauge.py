"""V6 (exact).  Gauge equivalence of the U2 and U3 recurrences (used in section 9; U3 itself is NOT proved).
R_j:  (n+1)^5 X_{n+1} + n^5 X_{n-1} = at_j(n) X_n,  at_2 = 2T(T^4+T^2+1/16), at_3 = 2T(T^4+5T^2+25/16), T = n+1/2.
Claim: if X satisfies R_2 for all n in a range, then  Z_n := P1(n) X_n + 32 (n+1)^5 X_{n+1},
       P1(n) = -(2n+1)(16n^4+56n^3+76n^2+46n+11),  satisfies R_3.   (Symbolic proof below.)
And with (w,v) the U2 sequences: G(w) = w^(3), G(v - 4w) = v^(3) (the U3 sequences: w3_0=1, v3_0=0, v3_1=1/8)."""
import sympy as sp
from fractions import Fraction as Fr
n = sp.symbols('n'); X0, X1 = sp.symbols('X0 X1')     # X_n, X_{n+1}; eliminate others with R_2
T = lambda m: m + sp.Rational(1, 2)
at2 = lambda m: 2 * T(m) * (T(m)**4 + T(m)**2 + sp.Rational(1, 16))
at3 = lambda m: 2 * T(m) * (T(m)**4 + 5 * T(m)**2 + sp.Rational(25, 16))
P1 = lambda m: -(2 * m + 1) * (16 * m**4 + 56 * m**3 + 76 * m**2 + 46 * m + 11)
Xm1 = (at2(n) * X0 - (n + 1)**5 * X1) / n**5                       # X_{n-1}
X2 = (at2(n + 1) * X1 - (n + 1)**5 * X0) / (n + 2)**5              # X_{n+2}
Z = lambda m, a, b: P1(m) * a + 32 * (m + 1)**5 * b
Zm1 = Z(n - 1, Xm1, X0); Z0 = Z(n, X0, X1); Zp1 = Z(n + 1, X1, X2)
res = sp.simplify(sp.together((n + 1)**5 * Zp1 + n**5 * Zm1 - at3(n) * Z0))
print("R_3 applied to G(X):", res); assert res == 0
def seq(at, x0, x1, N):
    s = [Fr(x0), Fr(x1)]
    for m in range(1, N):
        s.append((at(m) * s[m] - m**5 * s[m - 1]) / (m + 1)**5)
    return s
f2 = lambda m: Fr(8 * (m**5 + (m + 1)**5) - 12 * (m**3 + (m + 1)**3) + 7 * (2 * m + 1), 8)
f3 = lambda m: Fr(8 * (m**5 + (m + 1)**5) + 20 * (m**3 + (m + 1)**3) - 5 * (2 * m + 1), 8)
N = 40
w2 = seq(f2, 1, f2(0), N + 1); v2 = seq(f2, 0, Fr(1, 8), N + 1)
w3 = seq(f3, 1, f3(0), N); v3 = seq(f3, 0, Fr(1, 8), N)
P1i = lambda m: -(2 * m + 1) * (16 * m**4 + 56 * m**3 + 76 * m**2 + 46 * m + 11)
G = lambda s, m: P1i(m) * s[m] + 32 * (m + 1)**5 * s[m + 1]
for m in range(0, N):
    assert G(w2, m) == w3[m] and G(v2, m) - 4 * G(w2, m) == v3[m], m
print("G(w) = w^(3), G(v - 4w) = v^(3) for n < %d: PASS" % N)
print("V6 PASS")
