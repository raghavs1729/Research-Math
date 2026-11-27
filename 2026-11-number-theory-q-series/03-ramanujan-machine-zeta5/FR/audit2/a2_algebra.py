"""Audit 2: independent symbolic checks of Theorem C, Steps 6-7.

(1) Operator identity: sum_j [(j-x) g_j - (A(j-1)+k) g_{j-1} - B(j-1-x') g_{j-2}] z^j
    == Delta*theta(G) - N*G  for a generic polynomial G (coefficientwise).
(2) Partial fractions of N/Delta: alpha_1 - alpha_2 = omega/(B(l1-l2)), with x + x' = -C/B.
(3) Double-root case: N(1/lambda) = omega/(2 B lambda) when B = -A^2/4.
"""
import sympy as sp

A, k, B, C, x, z = sp.symbols('A k B C x z')
xp = -C / B - x  # other root
N = x + k * z + B * (1 - xp) * z**2
Delta = 1 - A * z - B * z**2
omega = 2 * k * B - A * (B + C)

# (1) generic polynomial of degree 6 with symbolic coefficients, starting at z^0
g = sp.symbols('g0:7')
G = sum(g[i] * z**i for i in range(7))
theta = lambda F: sp.expand(z * sp.diff(F, z))
lhs = sp.expand(Delta * theta(G) - N * G)
gg = lambda j: g[j] if 0 <= j < 7 else 0
rhs = sum(((j - x) * gg(j) - (A * (j - 1) + k) * gg(j - 1) - B * (j - 1 - xp) * gg(j - 2)) * z**j
          for j in range(0, 9))
print("(1) operator identity:", sp.simplify(lhs - sp.expand(rhs)) == 0)

# (2) partial fractions with l1 + l2 = A, l1 l2 = -B
l1, l2 = sp.symbols('l1 l2')
sub = {A: l1 + l2, B: -l1 * l2}
Nl = sp.expand(N.subs(sub))
# alpha_i = residue coefficient at 1/(1-l_i z): alpha_1 = [N/(1-l2 z)] at z = 1/l1
a1 = sp.simplify(Nl.subs(z, 1 / l1) / (1 - l2 / l1))
a2 = sp.simplify(Nl.subs(z, 1 / l2) / (1 - l1 / l2))
target = (omega / (B * (l1 - l2))).subs(sub)
print("(2) alpha1-alpha2 == omega/(B(l1-l2)):", sp.simplify(a1 - a2 - target) == 0)
# also alpha1 + alpha2 (should be in F_p automatically): print it
print("    alpha1+alpha2 =", sp.factor(sp.simplify((a1 + a2).subs({}))))

# (3) double root
lam = A / 2
val = sp.simplify(N.subs(B, -A**2 / 4).subs(z, 1 / lam) - (omega / (2 * B * lam)).subs(B, -A**2 / 4))
print("(3) N(1/lambda) == omega/(2 B lambda) when D = 0:", val == 0)
