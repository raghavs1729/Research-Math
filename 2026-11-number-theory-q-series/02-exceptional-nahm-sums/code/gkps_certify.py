"""Series checks (NOT proofs; the proofs are in proofs/GKPS-signed.md) of the rank-three signed Nahm sum identities of
Gang--Kim--Park--Stubbs, arXiv:2411.06081, Section 4.1.3, in the form

    I(K, v, Q) := sum_{m in N^3} q^{m^T K m/2} (-q^{1/2})^{-v.m} q^{-Q.m} / ((q)_{m1}(q)_{m2}(q)_{m3}),      v = mu_0 - a,

and of the three lemmas used:  (L1) q-binomial theorem, (L2) q-Chu--Vandermonde in the form
sum_j (-1)^j q^{j(j-1)/2}/((q)_j (q)_{N-j} (q)_{M-j}) = q^{NM}/((q)_N (q)_M),  (L3) = Lemma 4 of proofs/E7.md.
usage: python3 gkps_certify.py [N]      (N = order in q, default 30; box bound adapts)
"""
import sys, itertools
from fractions import Fraction as Fr

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30

def mul(a, b):
    c = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j > N: break
                if y: c[i + j] += x * y
    return c

_IP = {}
def ip(n):
    """1/(q)_n mod q^{N+1}"""
    if n not in _IP:
        c = [0] * (N + 1); c[0] = 1
        for j in range(1, n + 1):
            for e in range(j, N + 1): c[e] += c[e - j]
        _IP[n] = c
    return _IP[n]

def signed_sum(K, v, Q=(0, 0, 0)):
    out = {}
    # positive definite K: exponent >= lam*|m|^2/2 - C|m| ; use a generous box and an exact cut
    B = 1
    while True:
        B += 1
        # smallest exponent on the faces m_i = B
        ok = True
        for i in range(3):
            for m in itertools.product(range(B + 1), repeat=2):
                mm = list(m[:i]) + [B] + list(m[i:])
                e = Fr(sum(K[a][b] * mm[a] * mm[b] for a in range(3) for b in range(3)), 2) - Fr(sum(v[a] * mm[a] for a in range(3)), 2) - sum(Q[a] * mm[a] for a in range(3))
                if e <= N: ok = False; break
            if not ok: break
        if ok and B > 6: break
        if B > 60: raise RuntimeError("box")
    for m in itertools.product(range(B + 1), repeat=3):
        e = Fr(sum(K[a][b] * m[a] * m[b] for a in range(3) for b in range(3)), 2) - Fr(sum(v[a] * m[a] for a in range(3)), 2) - sum(Q[a] * m[a] for a in range(3))
        if e > N: continue
        s = -1 if (sum(v[a] * m[a] for a in range(3)) % 2) else 1
        p = [1] + [0] * N
        for mi in m: p = mul(p, ip(mi))
        for k, c in enumerate(p):
            if c and e + k <= N: out[e + k] = out.get(e + k, 0) + s * c
    return {k: c for k, c in out.items() if c}

def prod_inv(res, mod):
    c = [0] * (N + 1); c[0] = 1
    for n in range(1, N + 1):
        if n % mod in res:
            for e in range(n, N + 1): c[e] += c[e - n]
    return {Fr(k): x for k, x in enumerate(c) if x}

G = prod_inv({1, 4}, 5); H = prod_inv({2, 3}, 5)
AG1 = prod_inv({2, 3, 4, 5}, 7); AG2 = prod_inv({1, 3, 4, 6}, 7); AG3 = prod_inv({1, 2, 5, 6}, 7)
def fermion(power):
    """(-q^{1/2};q)_inf^power as dict exponent -> coeff (exponents in (1/2)Z), up to q^N"""
    c = [0] * (2 * N + 1); c[0] = 1          # in t = q^{1/2}
    for _ in range(power):
        n = 1
        while n <= 2 * N:
            d = c[:]
            for e in range(n, 2 * N + 1): d[e] += c[e - n]
            c = d; n += 2
    return {Fr(k, 2): x for k, x in enumerate(c) if x}
def dmul(a, b):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 <= N: out[e1 + e2] = out.get(e1 + e2, 0) + c1 * c2
    return {k: c for k, c in out.items() if c}
def theta(lin):
    """sum_{u in Z} q^{u^2 + lin*u} / (q)_inf  (lin = 0 or 1)"""
    c = [0] * (N + 1)
    u = -N - 2
    while u <= N + 2:
        e = u * u + lin * u
        if 0 <= e <= N: c[e] += 1
        u += 1
    return dmul({Fr(k): x for k, x in enumerate(c) if x}, {Fr(k): x for k, x in enumerate(ip(N + 1)) if x})
def finv():
    """1/(-q^{1/2};q)_inf up to q^N"""
    c = [0] * (2 * N + 1); c[0] = 1
    n = 1
    while n <= 2 * N:
        d = c[:]
        for e in range(n, 2 * N + 1): d[e] -= d[e - n]
        c = d; n += 2
    return {Fr(k, 2): x for k, x in enumerate(c) if x}
def shift(a, e, sign=1):
    return {k + e: sign * c for k, c in a.items() if k + e <= N}
F1 = fermion(1); F2 = fermion(2)
FG = dmul(F1, G); FH = dmul(F1, H); FFG = dmul(F2, G); FFH = dmul(F2, H); ZERO = {}
ONE = {Fr(0): 1}; GG = dmul(G, G); HH = dmul(H, H)
T0G = dmul(theta(0), G); T0H = dmul(theta(0), H); T1G = dmul(theta(1), G); T1H = dmul(theta(1), H)
FI = finv(); HoF = dmul(FI, H); GoF = dmul(FI, G); mHoF = shift(HoF, Fr(-1, 2), -1)
AG9 = [prod_inv(set(range(1, 9)) - {i, 9 - i}, 9) for i in (1, 2, 3, 4)]
NAMES = {id(ONE): "1", id(GG): "G^2", id(HH): "H^2", id(T0G): "th0 G/(q)_inf", id(T0H): "th0 H/(q)_inf", id(T1G): "th1 G/(q)_inf",
         id(T1H): "th1 H/(q)_inf", id(HoF): "H/(-q^1/2)_inf", id(GoF): "G/(-q^1/2)_inf", id(mHoF): "-q^-1/2 H/(-q^1/2)_inf",
         id(AG9[0]): "AG mod 9, i=1", id(AG9[1]): "AG mod 9, i=2", id(AG9[2]): "AG mod 9, i=3", id(AG9[3]): "AG mod 9, i=4",
         id(FG): "(-q^1/2)_inf G", id(FH): "(-q^1/2)_inf H", id(FFG): "(-q^1/2)_inf^2 G", id(FFH): "(-q^1/2)_inf^2 H", id(ZERO): "0",
         id(G): "G", id(H): "H", id(AG1): "1/(q^2,q^3,q^4,q^5;q^7)", id(AG2): "1/(q,q^3,q^4,q^6;q^7)", id(AG3): "1/(q,q^2,q^5,q^6;q^7)"}

CASES = []
K34 = [[2, 1, 1], [1, 2, 1], [1, 1, 1]]
CASES += [("3-4 case 1", K34, (-2, -2, -3), (0, 0, 0), H), ("3-4 case 1, W_(1,1,1)", K34, (-2, -2, -3), (1, 1, 1), G),
          ("3-4 case 2", K34, (-2, -1, -2), (0, 0, 0), H), ("3-4 case 3", K34, (-1, -2, -2), (0, 0, 0), H)]
K32 = [[2, -1, -1], [-1, 2, 1], [-1, 1, 1]]
CASES += [("3-2 case 1", K32, (-2, 2, 1), (0, 0, 0), H), ("3-2 case 1, W_(1,-1,-1)", K32, (-2, 2, 1), (1, -1, -1), G),
          ("3-2 case 2", K32, (-2, 1, 0), (0, 0, 0), H)]
K310 = [[5, 4, -2], [4, 4, -2], [-2, -2, 2]]
CASES += [("3-10 version 1", K310, (-1, -2, 0), (0, 0, 0), G), ("3-10 version 1, W_(2,2,-1)", K310, (-1, -2, 0), (2, 2, -1), H),
          ("3-10 version 2", K310, (-2, -3, 0), (0, 0, 0), G), ("3-10 version 2, W_(2,2,-1)", K310, (-2, -3, 0), (2, 2, -1), H)]
K33 = [[2, 1, -1], [1, 2, 0], [-1, 0, 1]]
CASES += [("3-3", K33, (-2, -1, 0), (0, 0, 0), FH), ("3-3, W_(1,1,0)", K33, (-2, -1, 0), (1, 1, 0), FG)]
K39 = [[3, 2, 1], [2, 2, 1], [1, 1, 1]]
CASES += [("3-9", K39, (-2, -2, 0), (0, 0, 0), FH), ("3-9, W_(1,1,0)", K39, (-2, -2, 0), (1, 1, 0), FG)]
K36 = [[3, -1, -1], [-1, 1, 0], [-1, 0, 1]]
for mu0 in ((-4, 1, 1), (-2, 0, 1), (-2, 1, 0)):
    v = (mu0[0] - 1, mu0[1], mu0[2])
    CASES += [("3-6 mu0=%s" % (mu0,), K36, v, (0, 0, 0), ZERO), ("3-6 mu0=%s, W_(1,0,0)" % (mu0,), K36, v, (1, 0, 0), ZERO)]
K37 = [[3, -1, 1], [-1, 1, 0], [1, 0, 1]]
for mu0 in ((-2, 1, -1), (-4, 1, -2)):
    v = (mu0[0] - 1, mu0[1], mu0[2])
    CASES += [("3-7 mu0=%s" % (mu0,), K37, v, (0, 0, 0), ZERO), ("3-7 mu0=%s, W_(1,0,0)" % (mu0,), K37, v, (1, 0, 0), ZERO)]
CASES += [("3-7 (3)", K37, (-2, 0, 0), (0, 0, 0), FFH), ("3-7 (3), W_(1,0,0)", K37, (-2, 0, 0), (1, 0, 0), FFG)]
K38 = [[3, 2, -1], [2, 2, -1], [-1, -1, 1]]
CASES += [("3-8", K38, (-3, -3, 1), (0, 0, 0), ZERO), ("3-8, W_(1,1,0)", K38, (-3, -3, 1), (1, 1, 0), ZERO)]
K323 = [[5, -3, 2], [-3, 3, -1], [2, -1, 1]]
CASES += [("3-23", K323, (-5, 3, -3), (0, 0, 0), ZERO), ("3-23, W_(2,-1,1)", K323, (-5, 3, -3), (2, -1, 1), ONE),
          ("3-23, W_(3,-2,1)", K323, (-5, 3, -3), (3, -2, 1), ONE), ("3-23, W_(4,-3,2)", K323, (-5, 3, -3), (4, -3, 2), ZERO)]
K35 = [[3, -2, 2], [-2, 2, -1], [2, -1, 2]]
CASES += [("3-5", K35, (0, 0, -1), (0, 0, 0), FG), ("3-5, W_(1,-1,1)", K35, (0, 0, -1), (1, -1, 1), FH)]
K313 = [[4, 2, 2], [2, 2, 1], [2, 1, 2]]
CASES += [("3-13", K313, (-4, -2, -2), (0, 0, 0), HH), ("3-13, W_(2,1,1)", K313, (-4, -2, -2), (2, 1, 1), GG)]
K314 = [[4, -2, 2], [-2, 2, -1], [2, -1, 2]]
CASES += [("3-14 (1)", K314, (-4, 2, -2), (0, 0, 0), T1H), ("3-14 (1), W_(1,-1,1)", K314, (-4, 2, -2), (1, -1, 1), T0H),
          ("3-14 (1), W_(2,-1,1)", K314, (-4, 2, -2), (2, -1, 1), T0G), ("3-14 (1), W_(3,-2,2)", K314, (-4, 2, -2), (3, -2, 2), T1G),
          ("3-14 (2)", K314, (-5, 3, -3), (0, 0, 0), mHoF), ("3-14 (2), W_(1,-1,1)", K314, (-5, 3, -3), (1, -1, 1), HoF),
          ("3-14 (2), W_(2,-1,1)", K314, (-5, 3, -3), (2, -1, 1), GoF), ("3-14 (2), W_(3,-2,2)", K314, (-5, 3, -3), (3, -2, 2), GoF)]
K318 = [[6, 4, 2], [4, 4, 2], [2, 2, 2]]
CASES += [("3-18", K318, (-6, -4, -2), (0, 0, 0), AG9[0]), ("3-18, W_(1,1,1)", K318, (-6, -4, -2), (1, 1, 1), AG9[1]),
          ("3-18, W_(2,2,1)", K318, (-6, -4, -2), (2, 2, 1), AG9[2]), ("3-18, W_(3,2,1)", K318, (-6, -4, -2), (3, 2, 1), AG9[3])]
for a in (1, 2, 3, 4, -1, -2, -3, -4):
    K = [[a * a + 1, a * a, a], [a * a, a * a, a], [a, a, 2]]
    mu0 = (-1, -2, -1) if a > 0 else (1, 0, -1)
    v = tuple(mu0[i] - (a, a, 1)[i] for i in range(3))
    CASES += [("3-12 a=%d case 1" % a, K, v, (0, 0, 0), H), ("3-12 a=%d case 1, W_(a,a,1)" % a, K, v, (a, a, 1), G)]
    mu0 = (0, -1, -1); v = tuple(mu0[i] - (a, a, 1)[i] for i in range(3))
    CASES += [("3-12 a=%d case 2" % a, K, v, (0, 0, 0), H), ("3-12 a=%d case 2, W_(a,a,1)" % a, K, v, (a, a, 1), G)]
K316 = [[4, -2, -1], [-2, 2, 1], [-1, 1, 1]]
CASES += [("3-16", K316, (-4, 2, 1), (0, 0, 0), AG1), ("3-16, W_(1,0,0)", K316, (-4, 2, 1), (1, 0, 0), AG2), ("3-16, W_(2,-1,-1)", K316, (-4, 2, 1), (2, -1, -1), AG3)]
K317 = [[9, 5, 3], [5, 4, 1], [3, 1, 2]]
CASES += [("3-17", K317, (-5, -4, -2), (0, 0, 0), AG1), ("3-17, W_(2,1,1)", K317, (-5, -4, -2), (2, 1, 1), AG2), ("3-17, W_(3,2,1)", K317, (-5, -4, -2), (3, 2, 1), AG3)]

def lemmas():
    ok = True
    def poch(n):  # (q)_n as list
        c = [1] + [0] * N
        for j in range(1, n + 1):
            d = c[:]
            for e in range(j, N + 1): d[e] -= c[e - j]
            c = d
        return c
    def qpow(e):
        c = [0] * (N + 1)
        if 0 <= e <= N: c[e] = 1
        return c
    # (L2)
    for NN in range(0, 7):
        for M in range(0, 7):
            lhs = [0] * (N + 1)
            for j in range(0, min(NN, M) + 1):
                t = mul(mul(mul(ip(j), ip(NN - j)), ip(M - j)), qpow(j * (j - 1) // 2))
                lhs = [x + (-1) ** j * y for x, y in zip(lhs, t)]
            rhs = mul(mul(ip(NN), ip(M)), qpow(NN * M))
            ok &= (lhs == rhs)
    # (L1) with z = q^c: sum_j (-1)^j q^{j(j-1)/2 + c j}/((q)_j (q)_{M-j}) = (q^c;q)_M/(q)_M
    for c in range(0, 4):
        for M in range(0, 8):
            lhs = [0] * (N + 1)
            for j in range(0, M + 1):
                t = mul(mul(ip(j), ip(M - j)), qpow(j * (j - 1) // 2 + c * j))
                lhs = [x + (-1) ** j * y for x, y in zip(lhs, t)]
            num = [1] + [0] * N
            for i in range(M):
                d = num[:]
                for e in range(c + i, N + 1): d[e] -= num[e - c - i]
                num = d
            ok &= (lhs == mul(num, ip(M)))
    # (L3) Lemma 4 of E7.md, times (q)_inf: sum_m q^{m^2+Rm} (q)_inf/((q)_m (q)_{m+R+a}) = sum_j [a;j] q^{j^2+Rj}
    qinf = poch(N + 1)
    for R in range(0, 4):
        for a in range(0, 5):
            lhs = [0] * (N + 1)
            m = 0
            while m * m + R * m <= N:
                lhs = [x + y for x, y in zip(lhs, mul(mul(ip(m), ip(m + R + a)), qpow(m * m + R * m)))]
                m += 1
            lhs = mul(lhs, qinf)
            rhs = [0] * (N + 1)
            for j in range(0, a + 1):
                rhs = [x + y for x, y in zip(rhs, mul(mul(mul(poch(a), ip(j)), ip(a - j)), qpow(j * j + R * j)))]
            ok &= (lhs == rhs)
    # (L4) finite pentagonal number theorem: sum_{i=-n}^{n} (-1)^i q^{i(3i-1)/2}/((q)_{n+i}(q)_{n-i}) = 1/(q)_n
    for n in range(0, 7):
        lhs = [0] * (N + 1)
        for i in range(-n, n + 1):
            tt = mul(mul(ip(n + i), ip(n - i)), qpow(i * (3 * i - 1) // 2))
            lhs = [x + (-1) ** i * y for x, y in zip(lhs, tt)]
        ok &= (lhs == ip(n))
    # (L5) Durfee rectangle: sum_c q^{(A-c)(B-c)}/((q)_c (q)_{A-c} (q)_{B-c}) = 1/((q)_A (q)_B)
    for A in range(0, 6):
        for B in range(0, 6):
            lhs = [0] * (N + 1)
            for c in range(0, min(A, B) + 1):
                lhs = [x + y for x, y in zip(lhs, mul(mul(mul(ip(c), ip(A - c)), ip(B - c)), qpow((A - c) * (B - c))))]
            ok &= (lhs == mul(ip(A), ip(B)))
    return ok

if __name__ == "__main__":
    allok = lemmas()
    print("series compared up to q^%d" % (N - 1))
    print("lemmas (L1)-(L5):", "OK" if allok else "FAILED")
    for name, K, v, Q, rhs in CASES:
        lhs = signed_sum(K, v, Q)
        # compare up to q^{N-1} (a right-hand side with a prefactor q^{-1/2} is only known up to q^{N-1/2})
        good = ({k: c for k, c in lhs.items() if k <= N - 1} == {k: c for k, c in rhs.items() if k <= N - 1})
        allok &= good
        print("%-28s K=%s  mu0-a=%s  Q=%s  = %-26s : %s" % (name, K, v, Q, NAMES[id(rhs)], "OK" if good else "FAILED"))
        sys.stdout.flush()
    print("ALL CHECKS PASSED" if allok else "SOME CHECK FAILED")
