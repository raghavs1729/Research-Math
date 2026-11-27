# arXiv:2611.05196, concluding conjecture: for p in {3,5,7}, alpha >= 1,
#   COS(p^(2a) n + (7 p^(2a) + 1)/8) = CES(same) = 0 mod p^a.
# COS/CES = sum of all odd/even parts in all cubic partitions of n (even parts in two colours):
#   sum COS(n) q^n = (1/(f_1 f_2)) * sum_n sigma_odd(n) q^n ;  sum CES(n) q^n = (1/(f_1 f_2)) * 2 sum_n (sigma(n) - sigma_odd(n)) q^n.
import flint, sys
def run(p, a, D):
    u = p**(a+2)
    def f(d):
        c = [0]*D; k = 0
        while True:
            e1 = d*k*(3*k-1)//2; e2 = d*k*(3*k+1)//2
            if e1 >= D: break
            c[e1] = (c[e1] + (-1)**k) % u
            if k and e2 < D: c[e2] = (c[e2] + (-1)**k) % u
            k += 1
        return flint.nmod_poly(c, u)
    C = (f(1).mul_low(f(2), D)).inverse_series_trunc(D)
    so = [0]*D; se = [0]*D
    for d in range(1, D):
        for m in range(d, D, d):
            if d % 2: so[m] += d
            else: se[m] += 2*d
    COS = C.mul_low(flint.nmod_poly([x % u for x in so], u), D); CES = C.mul_low(flint.nmod_poly([x % u for x in se], u), D)
    m = p**(2*a); t = (7*m + 1)//8
    ns = list(range(t, D, m))
    def v(x):
        x = int(x)
        if x == 0: return a+2
        k = 0
        while x % p == 0: x //= p; k += 1
        return k
    vc = [v(COS[n]) for n in ns]; ve = [v(CES[n]) for n in ns]
    print(f'p={p} alpha={a}: progression {m}n+{t}, {len(ns)} terms < {D}: min v_p(COS) = {min(vc)}, min v_p(CES) = {min(ve)}  (conjecture: >= {a})',
          '' if min(vc) >= a and min(ve) >= a else '  <-- FAILS at n index %d' % next(i for i in range(len(ns)) if vc[i] < a or ve[i] < a), flush=True)
    return COS, CES
# sanity: cubic partitions of 2: 2, 2', 1+1 -> COS(2) = 2, CES(2) = 4 ; COS(1) = 1
COS, CES = run(3, 1, 3000)
assert int(COS[1]) == 1 and int(COS[2]) == 2 and int(CES[2]) == 4
for p, a, D in [(5,1,20000),(7,1,40000),(3,2,60000),(3,3,400000),(5,2,400000),(7,2,1500000),(3,4,3000000),(5,3,6000000)]:
    run(p, a, D)
