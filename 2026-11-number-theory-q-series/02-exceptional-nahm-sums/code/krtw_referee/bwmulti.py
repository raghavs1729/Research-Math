# Referee: the explicit 2k-fold "Bartlett-Warnaar multisum" of section 5 at y=(sqrt x,-sqrt x),
# derived independently from BW Theorem 3.1 (r^(2l-1)=r^(2l) forced, f^(1)_{r,r} f^(1)_{r,s} = f^(2)_{r,s}),
# compared with K_k (= B_{k+1} J_k), k = 1,2,3. Also the k = 1 formula claimed in the last paragraph of section 5
# compared with q^N P'_{(2^N)}(sqrt x,-sqrt x;q) computed from Kostka-Foulkes polynomials.
import sys
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, fmpz_poly
import direct
X = int(sys.argv[1]); N = int(sys.argv[2])

def P(lst): return fmpq_poly(lst)
def qp_inv(sign, n):   # 1/(sign*q;q)_n for n>=0, i.e. 1/prod_{i=1}^n (1 - sign q^i)
    p = P([1])
    for i in range(1, n + 1):
        geo = [0] * N
        for j in range(0, N, i): geo[j] = sign ** (j // i)
        p = p.mul_low(P(geo), N)
    return p
def neg_q_poch(j):  # (-q^{1-j};q)_j = q^{-j(j-1)/2} * prod_{i=0}^{j-1}(1+q^i)  -> returns (shift, poly)
    p = P([1])
    for i in range(0, j): p = p * P([1] + [0] * (i - 1) + [1]) if i > 0 else p * 2
    return -j * (j - 1) // 2, p
def inv_mq(n):  # 1/(-q)_n for any integer n -> (shift, poly)
    if n >= 0: return 0, qp_inv(-1, n)
    return neg_q_poch(-n)
def omega(d):   # (shift, poly)
    d = abs(d)
    if d == 0: return 0, P([1])
    return d * (d - 1) // 2, P([1] + [0] * (d - 1) + [1]) * fmpq(1, 2)

def multisum(k):
    out = [P(0) for _ in range(X + 1)]
    minexp = [10 ** 9]
    def rec(l, A, Bv, wl):
        if l == k:
            A2 = A + [0]; B2 = Bv + [0]
            sh = 0; p = P([1])
            for i in range(k):
                a, b, an, bn = A2[i], B2[i], A2[i + 1], B2[i + 1]
                sh += a * a + b * b
                s1, p1 = omega(a - b); s2, p2 = inv_mq(a - bn); s3, p3 = inv_mq(b - an)
                sh += s1 + s2 + s3
                p = p.mul_low(p1, N + 200).mul_low(p2, N + 200).mul_low(p3, N + 200)
                p = p.mul_low(qp_inv(1, a - an), N + 200).mul_low(qp_inv(1, b - bn), N + 200)
            minexp[0] = min(minexp[0], sh)
            assert sh >= 0, (A, Bv, sh)
            if sh < N:
                out[sum(A) + sum(Bv)] += p.left_shift(sh).truncate(N)
            return
        ta = A[-1] if A else X; tb = Bv[-1] if Bv else X
        for a in range(0, ta + 1):
            for b in range(0, tb + 1):
                if a + b > wl: break
                rec(l + 1, A + [a], Bv + [b], wl - a - b)
    rec(0, [], [], X)
    return out, minexp[0]

for k in (1, 2, 3):
    ms, mn = multisum(k)
    K = direct.K(k, X, N)
    ok = all(ms[i] == fmpq_poly([int(v) for v in K.c[i].coeffs()]) for i in range(X + 1))
    print("k=%d: BW 2k-fold multisum at (sqrt x,-sqrt x) == K_k  mod x^%d q^%d : %s  (min q-exponent of a term: %d)" % (k, X + 1, N, ok, mn))

# k = 1 last paragraph of section 5: sum_{a+b=N} x^N q^{a^2+b^2-N} omega(a-b) 2^{[a!=b]} / ((Q)_a (Q)_b)
def Qinv(n):
    p = P([1])
    for i in range(1, n + 1):
        geo = [0] * N
        for j in range(0, N, 2 * i): geo[j] = 1
        p = p.mul_low(P(geo), N)
    return p
ms1, _ = multisum(1)
claimed_P = []   # claimed P'_{(2^n)} (no q^{|lambda|})
claimed_qP = []
for n in range(0, X + 1):
    tot = P(0)
    for a in range(0, n + 1):
        b = n - a
        s, w = omega(a - b)
        e = a * a + b * b - n + s
        t = (w * (2 if a != b else 1)).mul_low(Qinv(a), N).mul_low(Qinv(b), N)
        tot += t.left_shift(e).truncate(N) if e >= 0 else None
    claimed_P.append(tot)
    claimed_qP.append(tot.left_shift(n).truncate(N))
print("k=1 claimed formula (as P'_{2^N}) times q^N == multisum k=1 :", all(claimed_qP[i] == ms1[i] for i in range(X + 1)))
print("k=1 claimed formula itself == multisum k=1 (i.e. as q^N P') :", all(claimed_P[i] == ms1[i] for i in range(X + 1)))
for n in range(0, 4):
    print("  n=%d  multisum coeff: %s\n        claimed      : %s" % (n, ms1[n], claimed_P[n]))
