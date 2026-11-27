# Referee: numerical checks of Theorem B (B1)-(B5), the Bailey pair of (B1), and the limits, k <= 5.
import sys
from ser import S, monomial, poch, poch_inf, qpoch_inv_table
from flint import fmpz_poly
import direct
X = int(sys.argv[1]); N = int(sys.argv[2])
one = monomial(X, N, 0, 0)
def mono(a, b, c=1): return monomial(X, N, a, b, c)
def scal(s, c): return S(X, N, [p * c for p in s.c])
iQ = qpoch_inv_table(2, N, N); iq = qpoch_inv_table(1, N, N)
def invQ(n): s = S(X, N); s.c[0] = iQ[n]; return s
def invq(n): s = S(X, N); s.c[0] = iq[n]; return s
def P(c, a, b, step, n):  # (c x^a q^b ; q^step)_n
    return poch(X, N, [(c, a, b + step * j) for j in range(n) if b + step * j < N or a > 0])
def Pinf(c, a, b, step): return poch_inf(X, N, c, a, b, step)

xq_inf = Pinf(1, 1, 1, 1); x2Q_inf = Pinf(1, 2, 2, 2)
def bterm(k, n):
    # (-1)^n x^{(k+1)n} q^{(k+2)n^2} (x^2;Q)_n (1-xq^{2n}) / ((Q;Q)_n (1-x))
    if (k + 2) * n * n >= N: return None
    if n == 0: base = one
    else: base = P(-1, 1, 0, 0, 1) * P(1, 2, 2, 2, n - 1) * (one - mono(1, 2 * n))   # (1+x)(x^2Q;Q)_{n-1}(1-xq^{2n})
    return scal(mono((k + 1) * n, (k + 2) * n * n) * base * invQ(n), (-1) ** n)
def jterm(k, m):
    e = (k + 2) * m * m - m
    if e >= N: return None
    if m == 0: base = one
    else: base = P(1, 2, 2, 2, m - 1) * (one - mono(2, 4 * m))
    return scal(mono((k + 1) * m, e) * base * invQ(m), (-1) ** m)
def sumterms(f, k):
    s = S(X, N); n = 0
    while True:
        t = f(k, n)
        if t is None: break
        s = s + t; n += 1
    return s

C1inv = direct.C(1, X, N).inv()
res = []
# Bailey pair for (B1): beta_n = 1/(Q;Q)_n relative to a = x, base q
for n in range(0, 7):
    lhs = S(X, N)
    for r in range(0, n + 1):
        a = bterm(-1, r)  # with k=-1: x^0 q^{n^2}: alpha_r = (-1)^r q^{r^2}(x^2;Q)_r(1-xq^{2r})/((Q)_r(1-x))
        if a is None: continue
        lhs = lhs + a * invq(n - r) * P(1, 1, 1, 1, n + r).inv()
    res.append(("Bailey pair beta_n=1/(Q;Q)_n rel. x, n=%d" % n, lhs == invQ(n)))
for k in range(1, 6):
    Bk1 = direct.B(k + 1, X, N); Ck = direct.C(k, X, N); Jk = Ck * C1inv
    res.append(("(B1) k=%d" % k, xq_inf * Bk1 == sumterms(bterm, k)))
    res.append(("(B2) k=%d" % k, x2Q_inf * Jk == sumterms(jterm, k)))
    # (B3)
    rhs = S(X, N)
    def rec(i, ts, wl):
        global rhs
        if i == k:
            e = sum(t * t for t in ts)
            if e >= N: return
            tk = ts[-1]
            term = mono(sum(ts), e) * P(-1, 1, 1, 2, tk) * P(-1, 1, 2, 2, tk).inv() * invQ(tk)
            for j in range(k - 1): term = term * invQ(ts[j] - ts[j + 1])
            rhs = rhs + term; return
        top = ts[-1] if ts else X
        for t in range(0, top + 1):
            if t > wl: break
            rec(i + 1, ts + [t], wl - t)
    rec(0, [], X)
    res.append(("(B3) k=%d" % k, Bk1 == Pinf(-1, 1, 2, 2) * rhs))
    # (B4): 2(1-x)(1-x^2) (xq)(x^2Q) B J  ==  sum_{n,m} (-1)^{n+m} x^{..} q^{..} numerator/((Q)_n(Q)_m)
    lhs = scal(xq_inf * x2Q_inf * Bk1 * Jk * (one - mono(1, 0)) * (one - mono(2, 0)), 2)
    L = S(X, N)
    for n in range(0, N):
        if (k + 2) * n * n - n >= N: break
        for m in range(0, N):
            e = (k + 2) * (n * n + m * m) - n - m
            if e >= N: break
            num = P(1, 2, 0, 2, n) * P(1, 2, 0, 2, m) * (mono(0, n) + mono(0, m)) * (one + mono(1, n + m)) \
                  * (one - mono(1, 2 * n)) * (one - mono(1, 2 * m))
            L = L + scal(mono((k + 1) * (n + m), e) * num * invQ(n) * invQ(m), (-1) ** (n + m))
    res.append(("(B4) k=%d" % k, lhs == L))
    # (B5)
    if k % 2 == 0:
        mm = k // 2
        def frakJ(e_):
            out = S(X, N)
            def rec2(i, ss, wl):
                nonlocal out
                if i == mm:
                    ex = 2 * sum(s * s for s in ss)
                    if ex >= N: return
                    term = mono(2 * sum(ss), ex) * invQ(ss[-1]) * P(-1, 1, e_, 2, ss[-1]).inv()
                    for j in range(mm - 1): term = term * invQ(ss[j] - ss[j + 1])
                    out = out + term; return
                top = ss[-1] if ss else X
                for s in range(0, top + 1):
                    if 2 * s > wl: break
                    rec2(i + 1, ss + [s], wl - 2 * s)
            rec2(0, [], X)
            return out
        res.append(("(B5) B_{2m+1}=(-xq;q)J2 k=%d" % k, Bk1 == Pinf(-1, 1, 1, 1) * frakJ(2)))
        res.append(("(B5) J_{2m}=J1 k=%d" % k, Jk == frakJ(1)))
    # J_k = AG_{(k-1)/2}(x^2;Q) for odd k
    if k % 2 == 1:
        j = (k - 1) // 2
        if j == 0: ag = one
        else:
            ag = S(X, N)
            def rec3(i, ss, wl):
                global ag
                if i == j:
                    ex = 2 * sum(s * s for s in ss)
                    if ex >= N: return
                    term = mono(2 * sum(ss), ex) * invQ(ss[-1])
                    for t in range(j - 1): term = term * invQ(ss[t] - ss[t + 1])
                    ag = ag + term; return
                top = ss[-1] if ss else X
                for s in range(0, top + 1):
                    if 2 * s > wl: break
                    rec3(i + 1, ss + [s], wl - 2 * s)
            rec3(0, [], X)
        res.append(("J_k = AG_{(k-1)/2}(x^2;Q) k=%d" % k, Jk == ag))
for lab, ok in res: print("%-50s %s" % (lab, ok))
print("ALL:", all(o for _, o in res))
