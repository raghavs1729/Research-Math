"""Checks for SR/THEOREM.md (families of arXiv:2412.16818, section 4.3 and eqs (8),(9)).
(1) Theorem SR: for z>0, alpha>=0, beta>-1, m = beta-alpha+1, z+m>0:  consecutive convergents of z + K_{n>=1} (n+alpha)(n+beta)/z
    bracket alpha*T_{alpha-1}/T_alpha (alpha>0), resp. 1/T_0 (alpha=0), T_c = int_0^oo e^{-zt} sech^m t tanh^c t dt.
(2) Family F1 (a_n = i+2j+1, b_n = n^2+(i+k)n): exact closed form of T_0 in Q + Q log 2 (k=0), Q + Q pi (k=1),
    computed symbolically, compared with quadrature; the coefficient of log 2 resp. pi is the nonzero number in the text.
(3) Family F2 (a_n = -2n+j-1, b_n = -n^2+jn+k): exact identity for the convergents and limit = smaller root.
(4) eqs (8),(9): PCF[4, 4n^2+2n], PCF[4, 4n^2+2n-2], PCF[4, 4n^2-2n] against the Gauss/lemniscate expressions."""
from fractions import Fraction as Fr
from mpmath import mp, mpf, quad, exp, cosh, tanh, inf, pi, log, sqrt, agm, binomial
mp.dps = 30
def convs(a, b, N):
    p0, p1, q0, q1 = mpf(1), mpf(a(0)), mpf(0), mpf(1); out = []
    for n in range(1, N + 1):
        p0, p1 = p1, a(n) * p1 + b(n) * p0; q0, q1 = q1, a(n) * q1 + b(n) * q0
        s = abs(q1) if q1 != 0 else (abs(q0) if q0 != 0 else mpf(1)); p0 /= s; p1 /= s; q0 /= s; q1 /= s
        out.append(p1 / q1 if q1 != 0 else None)
    return out
T = lambda z, m, c: quad(lambda t: exp(-z * t) / cosh(t)**m * tanh(t)**c, [0, 1, 4, 16, inf])
# (1)
for (z, al, be) in ((mpf(2), mpf(0), mpf(1) / 2), (mpf(2), mpf(1), mpf(-1) / 2), (mpf(3), mpf(0.7), mpf(2.2)), (mpf(1), mpf(0), mpf(0)), (mpf(5), mpf(2), mpf(0.25))):
    m = be - al + 1
    val = 1 / T(z, m, 0) if al == 0 else al * T(z, m, al - 1) / T(z, m, al)
    c = convs(lambda n: z, lambda n: (n + al) * (n + be), 4001)
    lo, hi = min(c[-1], c[-2]), max(c[-1], c[-2])
    assert lo <= val <= hi and hi - lo < mpf(10)**-2, (z, al, be)
    print("(1) z=%s alpha=%s beta=%s: value %s lies between convergents 4000, 4001 (gap %s)" % (mp.nstr(z, 3), mp.nstr(al, 3), mp.nstr(be, 3), mp.nstr(val, 15), mp.nstr(hi - lo, 3)))
# (2) exact evaluation of T_0 = 2^m int_0^1 u^(z+m-1) (1+u^2)^(-m) du as  rat + coef*const  (Fractions)
from math import comb
def T0_exact(i, j, k):
    if k == 0:      # 2^i int_0^1 w^(i+j) (1+w)^(-i-1) dw, w^(i+j) = sum C(i+j,l)(-1)^(i+j-l)(1+w)^l
        rat = Fr(0); coef = Fr(0)
        for l in range(i + j + 1):
            c = comb(i + j, l) * (-1)**(i + j - l)
            e = l - i - 1
            if e == -1: coef += c
            else: rat += c * (Fr(2)**(e + 1) - 1) / (e + 1)
        return 2**i * rat, 2**i * coef
    M, R = i + j + 1, i + 2     # 2^(i+2) int_0^1 u^(2M) (1+u^2)^(-R) du, u^(2M) = sum C(M,l)(-1)^(M-l)(1+u^2)^l
    # I_r = int_0^1 (1+u^2)^(-r) du = q_r + p_r*pi:  I_1 = pi/4, I_r = 1/((2r-2) 2^(r-1)) + (2r-3)/(2r-2) I_(r-1)
    I = {1: (Fr(0), Fr(1, 4))}
    for r in range(2, R + 1):
        q, pc = I[r - 1]; f = Fr(2 * r - 3, 2 * r - 2)
        I[r] = (Fr(1, (2 * r - 2) * 2**(r - 1)) + f * q, f * pc)
    rat = Fr(0); coef = Fr(0)
    for l in range(M + 1):
        c = comb(M, l) * (-1)**(M - l)
        e = l - R
        if e >= 0:      # int_0^1 (1+u^2)^e du
            rat += c * sum(Fr(comb(e, t), 2 * t + 1) for t in range(e + 1))
        else:
            q, pc = I[-e]; rat += c * q; coef += c * pc
    return 2**(i + 2) * rat, 2**(i + 2) * coef
def half_binom(a2, n):      # binomial(a2/2, n) exactly
    r = Fr(1)
    for t in range(n): r = r * (Fr(a2, 2) - t) / (t + 1)
    return r
for i in range(0, 5):
    for j in range(0, 5):
        for k in (0, 1):
            z, m = i + 2 * j + 1, i + k + 1
            rat, coef = T0_exact(i, j, k)
            pred = Fr((-1)**j * 2**i * comb(i + j, i)) if k == 0 else (-1)**j * 2**i * half_binom(2 * i + 2 * j + 1, i + 1)
            assert coef == pred and coef != 0, (i, j, k, coef, pred)
            const = log(2) if k == 0 else pi
            val = mpf(rat.numerator) / rat.denominator + mpf(coef.numerator) / coef.denominator * const
            assert abs(T(mpf(z), mpf(m), 0) - val) < mpf(10)**-20, (i, j, k)
            if (i, j) in ((0, 0), (1, 0), (2, 1), (3, 3)):
                print("(2) i=%d j=%d k=%d: 1/PCF[%d, n^2+%dn] = %s + %s*%s" % (i, j, k, z, i + k, rat, coef, "log2" if k == 0 else "pi"))
print("    all 50 cases (i,j <= 4): T_0 = rational + c*log 2 (k=0) or rational + c*pi (k=1), c as predicted and nonzero; agrees with quadrature")
# (3)
for j in range(-4, 5):
    for k in range(-4, 5):
        D = j * j + 4 * k
        if D <= 0: continue
        r1 = (j + sqrt(mpf(D))) / 2; r2 = (j - sqrt(mpf(D))) / 2
        if any(abs(r - round(r)) < mpf(10)**-12 and round(r) >= 1 for r in (r1, r2)): continue
        # exact formula: p_n/q_n = r1 + (r2 - 1)/S_n,  S_n = sum_{l=0}^{n} prod_{i=1}^{l} (i - r1)/(i + 1 - r2)
        c = convs(lambda n: -2 * n + j - 1, lambda n: -n * n + j * n + k, 60)
        S = mpf(0); term = mpf(1)
        for n in range(0, 61):
            S += term
            if n >= 1 and c[n - 1] is not None and abs(S) > mpf(10)**-10:
                assert abs(c[n - 1] - (r1 + (r2 - 1) / S)) < mpf(10)**-15 * max(1, abs(c[n - 1])), (j, k, n)
            term *= (n + 1 - r1) / (n + 2 - r2)
        cc = convs(lambda n: -2 * n + j - 1, lambda n: -n * n + j * n + k, 20000)[-1]
        assert abs(cc - r2) < mpf(10)**-2
print("(3) F2: convergent identity exact-to-precision for all |j|,|k| <= 4 with positive discriminant; limit = smaller root")
# (4)
G = 1 / agm(1, sqrt(mpf(2))); varpi = pi * G
for (b, expr, name) in ((lambda n: 4 * n * n + 2 * n, 2 * G / (4 * G - 3), "2G/(4G-3)"),
                        (lambda n: 4 * n * n + 2 * n - 2, (4 * G - 1) / (3 * G - 2), "(4G-1)/(3G-2)"),
                        (lambda n: 4 * n * n - 2 * n, 6 / (4 - varpi), "6/(4-varpi)")):
    c = convs(lambda n: mpf(4), b, 20001)
    lo, hi = min(c[-1], c[-2]), max(c[-1], c[-2])
    assert lo <= expr <= hi
    print("(4) %s = %s lies between convergents 20000, 20001 (gap %s)" % (name, mp.nstr(expr, 15), mp.nstr(hi - lo, 3)))
print("PASS")
