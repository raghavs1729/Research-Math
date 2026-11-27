# (HL_k) directly: sum_{lambda_1<=k} q^{|lambda|} P'_{2lambda}(sqrt x,-sqrt x;q) from Kostka-Foulkes (charge), vs K_k.
# Also the k=1 formula of the last paragraph of section 5, summed over ordered pairs and over a>=b.
# run from this directory: python3 hl_direct.py
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
Xmax, NN = 7, 30
sys.argv = ['x', '2', str(NN)]
exec(open('rs_generic.py').read().split('def schur2')[0])
from flint import fmpq_poly, fmpq
import direct
def qpq(n):
    p = fmpq_poly([1])
    for j in range(1, n + 1): p = p * fmpq_poly([1] + [0] * (j - 1) + [-1])
    return p
def parts(n, m):
    def rec(n, mx):
        if n == 0: yield (); return
        for p in range(min(n, mx), 0, -1):
            for r in rec(n - p, p): yield (p,) + r
    return list(rec(n, m))
def Pprime_special(mu):
    tot = fmpq_poly(0); n = sum(mu)
    for n2 in range(0, n // 2 + 1):
        n1 = n - n2; d = n1 - n2
        if d % 2: continue
        nu = (n1, n2) if n2 else (n1,)
        sv = (-1) ** n2
        for pw, cf in KF(nu, mu).items():
            tot += fmpq_poly([0] * pw + [cf * sv])
    b = fmpq_poly([1])
    for part in set(mu): b = b * qpq(mu.count(part))
    return tot.mul_low(qinv(b), NN)
for k in (1, 2, 3):
    K = direct.K(k, Xmax, NN)
    ok = True
    for s in range(0, Xmax + 1):
        tot = fmpq_poly(0) if s else fmpq_poly([1])
        if s:
            for lam in parts(s, k):
                tot += Pprime_special(tuple(2 * p for p in lam)).left_shift(s).truncate(NN)
        ok &= (tot == fmpq_poly([int(v) for v in K.c[s].coeffs()]))
    print("(HL_%d) via Kostka-Foulkes, x^0..x^%d mod q^%d:" % (k, Xmax, NN), ok)
def Qinv(n):
    p = fmpq_poly([1])
    for i in range(1, n + 1):
        p = p.mul_low(fmpq_poly([1 if j % (2 * i) == 0 else 0 for j in range(NN)]), NN)
    return p
for rng in ("ordered", "a>=b"):
    ok = True
    for n in range(1, Xmax + 1):
        tot = fmpq_poly(0)
        for a in range(0, n + 1):
            b = n - a
            if rng == "a>=b" and a < b: continue
            d = abs(a - b)
            w = fmpq_poly([1]) if d == 0 else fmpq_poly([1] + [0] * (d - 1) + [1]) * fmpq(1, 2)
            e = a * a + b * b - n + (d * (d - 1) // 2)
            t = (w * (2 if a != b else 1)).mul_low(Qinv(a), NN).mul_low(Qinv(b), NN)
            tot += t.left_shift(e).truncate(NN)
        ok &= (tot == Pprime_special((2,) * n))
    print("k=1 formula of sec.5 last paragraph (sum %s) == P'_{(2^N)}(sqrt x,-sqrt x):" % rng, ok)
