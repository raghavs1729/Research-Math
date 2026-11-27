#!/usr/bin/env python3
"""GRH tail: R_{q,a}(x) holds for all x >= X0(q).  Rigorous (exact rationals + Arb).

Notation: L = log x, t = 1/L, b(L) = e^{-L} li(e^L) = e^{-L} Ei(L),
S(L) = sum_{k<14} k!/L^{k+1},  M(x) = (e x/L) li(x/e) - li(x)^2 = x^2 [ b(L-1)/L - b(L)^2 ].

Step 1 (enclosure of li).  For L >= 30:  0 <= b(L) - S(L) <= KAPPA * 14!/L^15.
   Proof.  Integrating Ei'(u) = e^u/u by parts 14 times from u0 = 30,
     b(L) - S(L) = e^{-L} C0 + 14! e^{-L} int_30^L e^u u^{-15} du,  C0 = Ei(30) - e^30 S(30).
   The integrand g has (log g)' = 1 - 15/u >= 1/2 on [30, L], so the integral is <= 2 g(L);
   C0 > 0 (Arb), and e^{-L} C0 <= c0 14!/L^15 with c0 = C0 e^{-30} 30^15/14!, because
   L^15 e^{-L} decreases for L >= 15.  KAPPA := 2 + c0 (rounded up to a rational).
Step 2 (Lemma A).  For L >= 31:  M(x) >= x^2 L^{-6} (1 + 14/L + 145/L^2).
   By Step 1 (lower bound for b(L-1), upper for b(L)) it suffices that
     P(t) := (1-t)^14 { t S(1/t - 1) - (S(1/t) + KAPPA 14! t^15)^2 - t^6 (1 + 14 t + 145 t^2) } >= 0
   on (0, 1/31].  P is a polynomial with rational coefficients, P = t^9 Q, Q(0) = 1412;
   Q > 0 on [0, 1/31] is checked by Arb on a subdivision.
Step 3.  Let Pi(x) = phi(q) pi(x;q,a) = li(x) + D(x).  Grenie-Molteni (J. Number Theory 2019,
   arXiv:1709.07609, Cor. 1.2 with K = Q, L = Q(zeta_q), C = {a}; hypothesis: RH for zeta_L,
   i.e. RH + GRH for the characters mod q):  for x >= 2
     |D(x)| <= sqrt(x) A(L) + li(2),  A(L) = (1/(2 pi) + 3/L) log Delta + (L/(8 pi) + 1/(4 pi) + 6/L) phi.
   Then Pi(x)^2 - (e x/L) Pi(x/e) <= -M(x) + 2 li(x) |D(x)| + D(x)^2 + (e x/L) |D(x/e)|, and after
   division by x^{3/2} the claim follows from  LHS(L) > RHS(L)  where
     LHS(L) = e^{L/2} L^{-6} (1 + 14/L + 145/L^2)                      (increasing for L >= 31)
     RHS(L) = 2 a(L) A(L)/L + sqrt(e) A(L-1)/L + li2 (2 a(L) + e)/(L e^{L/2}) + (A(L) + li2 e^{-L/2})^2 e^{-L/2}
     a(L)   = L (S(L) + KAPPA 14!/L^15) >= L li(x)/x                    (RHS is decreasing, see RAMANUJAN-AP.md)
   so it is enough to check it at L0 = log X0.
"""
from fractions import Fraction as Fr
from math import factorial
from flint import arb, fmpq, fmpq_poly, ctx
ctx.prec = 300
N = 14
# ---- Step 1
u0 = arb(30)
S30 = sum(arb(factorial(k)) / u0 ** (k + 1) for k in range(N))
C0 = arb(30).ei() - u0.exp() * S30
c0 = C0 * (-u0).exp() * u0 ** 15 / factorial(14)
print("C0 =", C0.str(20), " c0 =", c0.str(20))
assert C0 > 0
KAPPA = Fr(5)
assert c0 + 2 < 5, "KAPPA too small"
print("KAPPA =", KAPPA, " (need > 2 + c0 =", (c0 + 2).str(12), ")")
# sanity: enclosure against Arb's Ei at a few L
for Lv in [30, 31, 35, 41, 60, 100]:
    Lb = arb(Lv)
    b = Lb.ei() * (-Lb).exp()
    S = sum(arb(factorial(k)) / Lb ** (k + 1) for k in range(N))
    up = S + arb(5) * factorial(14) / Lb ** 15
    assert b - S > 0 and up - b > 0, Lv
# ---- Step 2: exact polynomial P(t)
def poly(coeffs):  # list of Fractions, ascending
    return fmpq_poly([fmpq(c.numerator, c.denominator) for c in coeffs])
t = poly([Fr(0), Fr(1)])
one_minus_t = poly([Fr(1), Fr(-1)])
St = sum((poly([Fr(factorial(k))]) * t ** (k + 1) for k in range(N)), poly([Fr(0)]))        # S(1/t)
# (1-t)^14 * t * S(1/t - 1) = t * sum k! t^{k+1} (1-t)^{14-k-1}
T1 = sum((poly([Fr(factorial(k))]) * t ** (k + 2) * one_minus_t ** (N - k - 1) for k in range(N)), poly([Fr(0)]))
sq = (St + poly([KAPPA * factorial(14)]) * t ** 15) ** 2
main = t ** 6 * poly([Fr(1), Fr(14), Fr(145)])
P = T1 - one_minus_t ** N * (sq + main)
co = [P[i] for i in range(P.degree() + 1)]
assert all(c == 0 for c in co[:9]), co[:10]
Q = fmpq_poly(co[9:])
print("P = t^9 Q, Q(0) =", Q[0], " deg Q =", Q.degree())
assert Q[0] == 1412
# Q > 0 on [0, 1/31]: subdivision with Arb (Horner on balls)
def Qball(x):
    r = arb(0)
    for i in range(Q.degree(), -1, -1):
        r = r * x + arb(Q[i])
    return r
nsub = 4000
lo_min = None
for i in range(nsub):
    a_, b_ = arb(i) / (31 * nsub), arb(i + 1) / (31 * nsub)
    x = a_.union(b_)
    v = Qball(x)
    assert v > 0, (i, v)
    lo = v.lower()
    lo_min = lo if lo_min is None or lo < lo_min else lo_min
print("Lemma A verified: Q >= %s > 0 on [0, 1/31] (%d subintervals)" % (lo_min.str(8), nsub))
# independent sanity check of Lemma A against Arb's li
for Lv in [31, 33, 36, 40, 45, 50, 70, 100, 300]:
    Lb = arb(Lv)
    Mx2 = (Lb - 1).ei() * (-(Lb - 1)).exp() / Lb - (Lb.ei() * (-Lb).exp()) ** 2
    low = (1 + arb(14) / Lb + arb(145) / Lb ** 2) / Lb ** 6
    assert Mx2 > low
    print("  L=%4d  M L^6/x^2 = %s  >= 1+14/L+145/L^2 = %s" % (Lv, (Mx2 * Lb ** 6).str(10), (low * Lb ** 6).str(10)))
# ---- Step 3
li2 = arb(2).log().ei()
PI = arb.pi()
def A(L, phi, logD):
    return (1 / (2 * PI) + 3 / L) * logD + (L / (8 * PI) + 1 / (4 * PI) + 6 / L) * phi
def aplus(L):
    return L * (sum(arb(factorial(k)) / L ** (k + 1) for k in range(N)) + arb(5) * factorial(14) / L ** 15)
def margin(L, phi, logD):
    lhs = (L / 2).exp() / L ** 6 * (1 + 14 / L + 145 / L ** 2)
    rhs = (2 * aplus(L) * A(L, phi, logD) / L + arb(1).exp().sqrt() * A(L - 1, phi, logD) / L
           + li2 * (2 * aplus(L) + arb(1).exp()) / (L * (L / 2).exp())
           + (A(L, phi, logD) + li2 * (-L / 2).exp()) ** 2 * (-L / 2).exp())
    return lhs, rhs
FIELDS = {1: (1, 1), 3: (2, 3), 4: (2, 4), 5: (4, 125), 8: (4, 256), 12: (4, 144)}
print("\nq  phi  Delta_L   smallest X0 on a 1%-grid   log X0   LHS(L0)   RHS(L0)")
res = {}
for q, (phi, D) in FIELDS.items():
    logD = arb(D).log()
    # search X0 = c * 10^k on a multiplicative grid (ratio 1.01), smallest passing
    lo, hi = 31.0, 60.0
    # bisection on L in floats, then certify at a rounded-up X0
    for _ in range(60):
        mid = (lo + hi) / 2
        l, r = margin(arb(mid), phi, logD)
        if l > r: hi = mid
        else: lo = mid
    import math
    X0 = math.exp(hi)
    k = math.floor(math.log10(X0)); mant = math.ceil(X0 / 10 ** k * 100) / 100 + 0.01
    X0r = Fr(round(mant * 100), 100) * 10 ** k
    L0 = (arb(X0r.numerator) / X0r.denominator).log()
    l, r = margin(L0, phi, logD)
    assert l > r and L0 > 31
    res[q] = X0r
    print(f"{q:2d}  {phi}   {D:5d}   X0 = {float(X0r):.3e}   {float(L0.mid()):.4f}   {l.str(8)}   {r.str(8)}   certified: {bool(l > r)}")
