"""AUDIT: direct tests for Theorems C and D from overpartition numbers (pbar = f_2/f_1^2 mod ell, pentagonal series).

 C:  pbar(16 n) = 0 mod 7 for n = 3 (8), (n/7) = +1, 16 n < NC.
 D:  pbar(1156 n) = 0 mod 19 for n = 3 (8), (n/17) = -1, (n/19) = +1, n <= nD;
     and, for ALL n <= nD, the modular form g of step 5 (q-expansion written by d.gp from its basis coordinates)
     against its direct definition  g(n) = h(4n),  h = (T(17^2) + 17^358) W,  W(m) = [(m/19)=1] (-1)^m pbar(m):
         g(n) = W(1156 n) + 17^358 ((4n/17) + 1) W(4n) + 17^717 W(4n/289)        (mod 19).
usage: direct_cd.py C NC   |   direct_cd.py D nD
"""
import sys, time, flint, resource
def leg(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1
def pent(d, D, ell):
    P = flint.nmod_poly([1], ell); k = 1
    while d*k*(3*k-1)//2 < D:
        sg = 1 if k % 2 == 0 else ell-1
        P[d*k*(3*k-1)//2] = sg
        if d*k*(3*k+1)//2 < D: P[d*k*(3*k+1)//2] = sg
        k += 1
    return P
def pbar(D, ell):
    t0 = time.time()
    f1 = pent(1, D, ell); g = f1.mul_low(f1, D); del f1
    g = g.inverse_series_trunc(D)
    pb = g.mul_low(pent(2, D, ell), D); del g
    assert [int(pb[i]) for i in range(8)] == [x % ell for x in (1, 2, 4, 8, 14, 24, 40, 64)]
    print(f'pbar mod {ell} for n < {D}: {time.time()-t0:.0f}s', flush=True)
    return pb
if sys.argv[1] == 'C':
    NC = int(sys.argv[2]); pb = pbar(NC, 7)
    ns = [n for n in range(3, (NC-1)//16 + 1, 8) if leg(n, 7) == 1]
    bad = [n for n in ns if int(pb[16*n]) != 0]
    print(f'C: pbar(16 n) mod 7, n = 3 (8), (n/7) = 1: {len(ns)} values of n <= {ns[-1]}; failures: {len(bad)} {bad[:5]}')
    ctrl = [n for n in range(3, 4000, 8) if leg(n, 7) == -1 and int(pb[16*n]) != 0]
    print(f'   (negative control: (n/7) = -1, n < 4000: {len(ctrl)} nonzero values)')
else:
    nD = int(sys.argv[2]); ell = 19; pb = pbar(1156*nD + 1, ell)
    W = lambda m: (int(pb[m]) * (1 if m % 2 == 0 else -1)) % ell if leg(m, ell) == 1 else 0
    c1 = pow(17, 358, ell); c2 = pow(17, 717, ell)
    gfile = [int(x) for x in open('d_g.txt').read().strip().strip('[]').split(',')]
    assert len(gfile) == 200001 and nD <= 200000
    mism = 0
    for n in range(nD + 1):
        v = W(1156*n) + c1*(leg(4*n, 17) + 1)*W(4*n)
        if n % 289 == 0: v += c2*W(4*n//289)
        if v % ell != gfile[n]: mism += 1
    ns = [n for n in range(3, nD + 1, 8) if leg(n, 17) == -1 and leg(n, 19) == 1]
    bad = [n for n in ns if int(pb[1156*n]) != 0]
    print(f'D: g(n) [basis expansion, d.gp] vs direct formula from pbar, all n <= {nD}: mismatches {mism}')
    print(f'D: pbar(1156 n) mod 19, n = 3 (8), (n/17) = -1, (n/19) = 1: {len(ns)} values of n <= {ns[-1]}; failures: {len(bad)} {bad[:5]}')
    ctrl = [n for n in range(3, 20000, 8) if leg(n, 17) == 1 and leg(n, 19) == 1 and int(pb[1156*n]) != 0]
    print(f'   (negative control: (n/17) = +1 instead, n < 20000: {len(ctrl)} nonzero values)')
print('peak RSS (GB):', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20, 2))
