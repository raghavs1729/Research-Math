# arXiv:2611.05623 Conjectures 1 and 2 (2-adic): falsification test at large k.
# A_k: (f_6/(f_1 f_2 f_3))^k ;  B_k: (f_4 f_6^2/(f_1 f_3 f_12))^k  (designated summands, all parts odd)
import flint, sys
from fractions import Fraction
def v2(n):
    if n == 0: return 10**9
    c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c
def f(d, D, u):
    c = [0]*D; k = 0
    while True:
        e1 = d*k*(3*k-1)//2; e2 = d*k*(3*k+1)//2
        if e1 >= D: break
        c[e1] = (c[e1] + (-1)**k) % u
        if k and e2 < D: c[e2] = (c[e2] + (-1)**k) % u
        k += 1
    return flint.nmod_poly(c, u)
u = 2**62; NN = int(sys.argv[2]); D = 4*NN + 4; KMAX = int(sys.argv[1])
gA = f(6, D, u).mul_low((f(1, D, u).mul_low(f(2, D, u), D).mul_low(f(3, D, u), D)).inverse_series_trunc(D), D)
gB = f(4, D, u).mul_low(f(6, D, u).pow_trunc(2, D), D).mul_low((f(1, D, u).mul_low(f(3, D, u), D).mul_low(f(12, D, u), D)).inverse_series_trunc(D), D)
assert [int(gA[i]) for i in range(4)] == [1, 1, 3, 5]
print('B_1 first terms', [int(gB[i]) for i in range(8)])
badA = []; badB = []
pA = flint.nmod_poly([1], u); pB = flint.nmod_poly([1], u)
for k in range(1, KMAX+1):
    pA = pA.mul_low(gA, D); pB = pB.mul_low(gB, D)
    a3 = k*(k+1)*(k+14)//6; b3 = k*(k+2)*(k+7)//6
    assert int(pA[3]) == a3 % u and int(pB[3]) == b3 % u
    mA = min(min(v2(int(pA[4*n+3])), 62) for n in range(NN)); mB = min(min(v2(int(pB[4*n+3])), 62) for n in range(NN))
    eA = 3 if k % 16 == 2 else v2(a3)
    eB = 2 if k % 8 == 1 else (3 if k % 8 == 6 else v2(b3))
    if mA != eA: badA.append((k, eA, mA))
    if mB != eB: badB.append((k, eB, mB))
print(f'k <= {KMAX}, n < {NN}: Conjecture 1 (A) mismatches (k, predicted, observed min v2): {badA[:20]}')
print(f'k <= {KMAX}, n < {NN}: Conjecture 2 (B) mismatches: {badB[:20]}')
