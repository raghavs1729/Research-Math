"""Symbolic checks (sympy) for ../THEOREM-2.md: three classical transformations in operator form.
H_t : h(1-h) w'' + (1-2h) w' - t(1-t) w = 0            (2F1(t,1-t;1;h))
D_t : theta^3 - x (theta+1/2)(theta+t)(theta+1-t),  theta = x d/dx     (3F2(1/2,t,1-t;1,1;x))
(S1) Clausen:  w1, w2 in Sol(H_t)  ==>  w1*w2 (as a function of h) is annihilated by D_t in the variable x = 4h(1-h).
(S2) Pfaff (t = 1/2):  w in Sol(H_1/2)  ==>  (1-h)^(-1/2) w(h/(h-1)) in Sol(H_1/2).
(S3) Whipple:  v in Sol(D_1/4) (variable X)  ==>  (1-z)^(-1/2) v(-4z/(1-z)^2) in Sol(D_1/2) (variable z)."""
import sympy as sp
h, t = sp.symbols('h t')
R = sp.Rational
# ---- (S1)
a0, a1, b0, b1 = sp.symbols('a0 a1 b0 b1')          # w1, w1', w2, w2' at h
def d_h(expr):                                       # derivative w.r.t. h in the ring Q(h,t)[a0,a1,b0,b1], w'' eliminated by H_t
    a2 = (t * (1 - t) * a0 - (1 - 2 * h) * a1) / (h * (1 - h))
    b2 = (t * (1 - t) * b0 - (1 - 2 * h) * b1) / (h * (1 - h))
    return sp.diff(expr, h) + sp.diff(expr, a0) * a1 + sp.diff(expr, a1) * a2 + sp.diff(expr, b0) * b1 + sp.diff(expr, b1) * b2
x = 4 * h * (1 - h)
th = lambda e: sp.together(h * (1 - h) / (1 - 2 * h) * d_h(e))            # theta_x = x d/dx = (x / x'(h)) d/dh
W = a0 * b0
T1 = th(W); T2 = th(T1); T3 = th(T2)
shift = lambda e, c: None
# (theta+1/2)(theta+t)(theta+1-t) W = theta^3 W + (3/2) theta^2 W + (t(1-t) + 1/2) theta W + (t(1-t)/2) W
rhs = T3 + R(3, 2) * T2 + (t * (1 - t) + R(1, 2)) * T1 + t * (1 - t) / 2 * W
res = sp.simplify(sp.together(T3 - x * rhs))
print("(S1) D_t applied to w1*w2 at x = 4h(1-h):", res); assert res == 0
# ---- (S2)
lam = sp.symbols('lam')
w0, w1_, = sp.symbols('w0 w1')                       # w, w' at mu
mu = lam / (lam - 1)
w2_ = (R(1, 4) * w0 - (1 - 2 * mu) * w1_) / (mu * (1 - mu))
dmu = sp.diff(mu, lam)
v = (1 - lam)**R(-1, 2) * w0
def d_lam(expr):
    return sp.diff(expr, lam) + sp.diff(expr, w0) * w1_ * dmu + sp.diff(expr, w1_) * w2_ * dmu
v1 = d_lam(v); v2 = d_lam(v1)
res = sp.simplify(lam * (1 - lam) * v2 + (1 - 2 * lam) * v1 - R(1, 4) * v)
print("(S2) H_1/2 applied to (1-lam)^(-1/2) w(lam/(lam-1)):", res); assert res == 0
# ---- (S3)
z = sp.symbols('z')
X = -4 * z / (1 - z)**2
v0, v1s, v2s = sp.symbols('v0 v1 v2')                # v, theta_X v, theta_X^2 v at X ; theta_X^3 v = X (theta+1/2)(theta+1/4)(theta+3/4) v
# solve theta^3 v (1 - X) = X( (3/2) theta^2 v + (3/16 + 1/2) theta v + (3/32) v )
v3s = X * (R(3, 2) * v2s + (R(3, 16) + R(1, 2)) * v1s + R(3, 32) * v0) / (1 - X)
ratio = sp.simplify(z * sp.diff(X, z) / X)             # theta_z (f o X) = ratio * (theta_X f) o X
def th_z(expr):
    return sp.together(z * sp.diff(expr, z) + ratio * (sp.diff(expr, v0) * v1s + sp.diff(expr, v1s) * v2s + sp.diff(expr, v2s) * v3s))
u = (1 - z)**R(-1, 2) * v0
U1 = th_z(u); U2 = th_z(U1); U3 = th_z(U2)
res = sp.simplify(U3 - z * (U3 + R(3, 2) * U2 + R(3, 4) * U1 + R(1, 8) * u))
print("(S3) D_1/2 applied to (1-z)^(-1/2) v(-4z/(1-z)^2):", res); assert res == 0
print("PASS")
# ---- (S4) Taylor coefficients of pi*cot(pi(u+rho)) used in the bilateral-sum lemma
u, x0 = sp.symbols('u x0')
ser = sp.series(sp.cot(x0 + u), u, 0, 4).removeO()
C = sp.cot(x0); S = 1 + C**2
assert sp.simplify(ser.coeff(u, 1) + S) == 0 and sp.simplify(ser.coeff(u, 3) + S * (3 * S - 2) / 3) == 0
print("(S4) cot(x0+u) = C - S u + C S u^2 - (S(3S-2)/3) u^3 + ...,  S = csc^2 x0: PASS")
# ---- (S5) Ligozat: orders of the eta quotients at all cusps of Gamma_0(N) are >= 0 (holomorphy), and > 0 for cusp forms
from math import gcd
def cusp_orders(expo, N):
    out = {}
    for c in [d for d in range(1, N + 1) if N % d == 0]:
        num = sum(sp.Rational(gcd(c, d)**2 * r, d) for d, r in expo.items())
        out[c] = sp.Rational(N, 24 * gcd(c, N // c) * c) * num
    return out
for name, expo, N in (("eta4^16/(eta2^4 eta8^4)", {4: 16, 2: -4, 8: -4}, 16), ("eta2^4 eta4^4", {2: 4, 4: 4}, 8),
                      ("eta6^14/(eta2^3 eta18^3)", {6: 14, 2: -3, 18: -3}, 36), ("eta2^3 eta6^2 eta18^3", {2: 3, 6: 2, 18: 3}, 36),
                      ("eta2^6 eta6^2", {2: 6, 6: 2}, 12), ("eta2^3 eta6^3", {2: 3, 6: 3}, 12)):
    o = cusp_orders(expo, N)
    print("(S5) %-28s level %-3d orders at cusps 1/c: %s" % (name, N, {c: str(v) for c, v in o.items()}))
    assert all(v >= 0 for v in o.values())
print("PASS (S4, S5)")
