"""Checks for Z/THEOREM.md:  8 S_1 = -16 L(f,1),  S_1 = sum_k A_k (H_{2k} - H_k - log 2),  A_k = C(2k,k)^4/256^k,
f = eta(2 tau)^4 eta(4 tau)^4.
(1) exact q-series: eta(2t)^4 eta(4t)^4 = theta2^4 theta3^2 theta4^2 / 16  (thetas in the nome q = e^{2 pi i t}), to O(q^400);
(2) M(s) := int_0^1 K'(k) k^s dk = (pi/4) Gamma((s+1)/2)^2/Gamma(s/2+1)^2  (quadrature, s = 0, 1.3, 4);
    M(2n) = (pi^2/4) a_n,  M'(2n) = (pi^2/2) a_n (H_{2n} - H_n - log 2),  a_n = ((1/2)_n/n!)^2   (exact identity of digammas, checked numerically);
(3) K'(k) = sum_n a_n k^(2n) ( -log k - 2 (H_{2n} - H_n - log 2) )  on (0,1)  (numerically at several k);
(4) int_0^1 K'(k)^2 dk = -pi^2 S_1  (quadrature vs series);
(5) L(f,1) = 2 pi int_0^oo f(it) dt = (1/(2 pi^2)) int_0^1 K'(k)^2 dk  (L-value from the functional-equation series; quadrature);
(6) hence S_1 = -2 L(f,1); also S_0 = 16 L(f,2)/pi^2 by the same computation."""
from mpmath import mp, mpf, pi, quad, ellipk, gamma, log, sqrt, exp, gammainc, harmonic, psi, hyper, nsum, inf, binomial, agm
mp.dps = 40
# (1)
NQ = 400
def mul(A, B):
    R = [0] * (NQ + 1)
    for i, x in enumerate(A):
        if x:
            for j in range(NQ + 1 - i):
                if B[j]: R[i + j] += x * B[j]
    return R
th3 = [0] * (NQ + 1); th4 = [0] * (NQ + 1); th2h = [0] * (NQ + 1)   # theta2 = 2 q^(1/4) * th2h(q),  th2h = sum q^(n(n+1))
n = 0
while n * n <= NQ:
    th3[n * n] += 1 if n == 0 else 2; th4[n * n] += (1 if n == 0 else 2) * (-1)**n; n += 1
n = 0
while n * (n + 1) <= NQ:
    th2h[n * (n + 1)] += 1; n += 1
t2 = mul(th2h, th2h); t2_4 = mul(t2, t2)                 # theta2^4 = 16 q * t2_4
rhs = mul(mul(t2_4, mul(th3, th3)), mul(th4, th4))        # theta2^4 theta3^2 theta4^2/16 = q * rhs
co = [0] * (NQ + 1); co[0] = 1
for m in range(1, NQ + 1):
    for step in (2 * m, 4 * m):
        if step > NQ: continue
        for _ in range(4):
            for i in range(NQ, step - 1, -1): co[i] -= co[i - step]
assert co == rhs
print("(1) PASS: eta(2t)^4 eta(4t)^4 = theta2^4 theta3^2 theta4^2/16 to O(q^%d)" % NQ)
a = [1] + [0] * NQ; 
fcoef = [0] + co[:NQ]            # f = sum_{n>=1} fcoef[n] q^n
# (2)
Kp = lambda k: pi / (2 * agm(1, k))        # K'(k) = K(sqrt(1-k^2)) = pi/(2 AGM(1,k))
M = lambda s: pi / 4 * gamma((s + 1) / 2)**2 / gamma(s / 2 + 1)**2
for s in (mpf(0), mpf(13) / 10, mpf(4)):
    v = quad(lambda k: Kp(k) * k**s, [0, mpf(1) / 100, mpf(1) / 2, 1])
    assert abs(v - M(s)) < mpf(10)**-25, s
an = lambda n: (binomial(2 * n, n) / mpf(4)**n)**2
h = lambda n: harmonic(2 * n) - harmonic(n)
for n in (0, 1, 5, 12):
    assert abs(M(mpf(2 * n)) - pi**2 / 4 * an(n)) < mpf(10)**-30
    dM = M(mpf(2 * n)) * (psi(0, n + mpf(1) / 2) - psi(0, n + 1))
    assert abs(dM - pi**2 / 2 * an(n) * (h(n) - log(2))) < mpf(10)**-30
print("(2) PASS: moments of K'")
# (3)
for k in (mpf(1) / 10, mpf(1) / 2, mpf(9) / 10):
    ser = nsum(lambda n: an(int(n)) * k**(2 * int(n)) * (-log(k) - 2 * (h(int(n)) - log(2))), [0, inf])
    assert abs(ser - Kp(k)) < mpf(10)**-20, (k, ser, Kp(k))
print("(3) PASS: series of K'(k)")
# (4)
I2 = quad(lambda k: Kp(k)**2, [0, mpf(1) / 1000, mpf(1) / 100, mpf(1) / 2, 1])
mp.dps = 25
S1 = nsum(lambda n: an(int(n))**2 * (h(int(n)) - log(2)), [0, inf])
S0 = hyper([mpf(1) / 2] * 4, [1, 1, 1], 1)
print("(4) int K'^2 = %s ; -pi^2 S_1 = %s" % (mp.nstr(I2, 20), mp.nstr(-pi**2 * S1, 20)))
assert abs(I2 + pi**2 * S1) < mpf(10)**-12
# (5) L-values of f (level 8, weight 4, sign +1)
mp.dps = 40
def Lval(s, k=4, lev=8):
    A = sqrt(mpf(lev)) / (2 * pi); tot = mpf(0)
    for n in range(1, NQ + 1):
        if fcoef[n] == 0: continue
        x = n / A
        tot += fcoef[n] * (gammainc(s, x) / (mpf(n) / A)**s + gammainc(k - s, x) / (mpf(n) / A)**(k - s))
    return tot / (A**s * gamma(s))
L1, L2 = Lval(1), Lval(2)
F = lambda t: sum(fcoef[n] * exp(-2 * pi * n * t) for n in range(1, 200) if fcoef[n])
half = 1 / sqrt(mpf(8))
# int_0^oo f(it) dt, using f(i/(8t)) = 64 t^4 f(it) (Fricke, sign +1) to fold (0, 1/sqrt 8) onto (1/sqrt 8, oo):
# int_0^{1/sqrt 8} f(it) dt = int_{1/sqrt 8}^oo 8 u^2 f(iu) du
for t in (mpf(1) / 5, mpf(3) / 10):
    assert abs(F(1 / (8 * t)) - 64 * t**4 * F(t)) < mpf(10)**-30
I0 = quad(lambda t: F(t) * (1 + 8 * t * t), [half, 1, 2, 4, 8, 16])
I1 = quad(lambda t: F(t) * 2 * t, [half, 1, 2, 4, 8, 16])      # int_0^oo f(it) t dt  (t -> 1/(8t) maps t dt f to itself)
print("(5) L(f,1) = %s ; 2 pi int f = %s ; int K'^2/(2 pi^2) = %s" % (mp.nstr(L1, 25), mp.nstr(2 * pi * I0, 25), mp.nstr(I2 / (2 * pi**2), 25)))
assert abs(L1 - 2 * pi * I0) < mpf(10)**-25 and abs(L1 - I2 / (2 * pi**2)) < mpf(10)**-25
IKK = quad(lambda k: pi / (2 * agm(1, sqrt((1 - k) * (1 + k)))) * Kp(k), [0, mpf(1) / 1000, mpf(1) / 2, mpf(999) / 1000, 1])
print("(6) S_1 + 2 L(f,1) = %s ;  S_0 - 16 L(f,2)/pi^2 = %s ;  L(f,2) - int K K'/(2 pi) = %s" % (mp.nstr(S1 + 2 * L1, 5), mp.nstr(S0 - 16 * L2 / pi**2, 5), mp.nstr(L2 - IKK / (2 * pi), 5)))
assert abs(S1 + 2 * L1) < mpf(10)**-12 and abs(L2 - IKK / (2 * pi)) < mpf(10)**-25
assert abs(L2 - 4 * pi**2 * I1) < mpf(10)**-25 and abs(S0 - 16 * L2 / pi**2) < mpf(10)**-12
L3 = Lval(3)
assert abs(L3 - pi**2 / 4 * L1) < mpf(10)**-25
# (7) Zudilin's F_j(1) = [eps^j] (1/pi^2) sum_{n>=0} Gamma(1/2+eps+n)^4/Gamma(1+eps+n)^4, by Taylor (finite differences
#     are avoided: differentiate the summand analytically through polygamma is heavy, so use mpmath.taylor on nsum)
from mpmath import taylor
mp.dps = 30
G = lambda e: nsum(lambda n: (gamma(mpf(1) / 2 + e + n) / gamma(1 + e + n))**4, [0, inf]) / pi**2
Fj = taylor(G, 0, 3)
print("(7) F_j(1) =", [mp.nstr(x, 15) for x in Fj])
print("    L(f,1)/F_1(1) = %s (expect -1/16 = -0.0625), L(f,2)/F_2(1) = %s (1/32 = 0.03125), L(f,3)/F_3(1) = %s (-3/256 = -0.01171875)"
      % (mp.nstr(L1 / Fj[1], 12), mp.nstr(L2 / Fj[2], 12), mp.nstr(L3 / Fj[3], 12)))
assert abs(L1 / Fj[1] + mpf(1) / 16) < mpf(10)**-8 and abs(L2 / Fj[2] - mpf(1) / 32) < mpf(10)**-8 and abs(L3 / Fj[3] + mpf(3) / 256) < mpf(10)**-7
print("PASS")
