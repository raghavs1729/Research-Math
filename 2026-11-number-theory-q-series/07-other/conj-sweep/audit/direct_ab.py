"""AUDIT: direct test of Theorems A/B and of the modular form h, from overpartition numbers only.

pbar mod ell is computed from the eta-quotient f_2/f_1^2 (Euler's pentagonal series), NOT from 1/phi(-q).
For every n <= Bd (Sturm bound) we form, straight from Shimura's formula applied to c0(n) := (-1)^n pbar(ell n),
    hd(n) := c0(p^2 n) + eps (n/p) p^(lam-1) c0(n) + p^(2lam-1) c0(n/p^2) - s eps p^(lam-1) c0(n)   (mod ell)
and check (i) hd(n) == h(n) for ALL n <= nmax, where h is the q-expansion written by ab.gp (basis expansion),
(ii) hd(n) == 0 for n = r mod 8, (iii) the congruence itself: pbar(ell p^2 n) == 0 for n = r (8), (n/p) = s.
usage: direct_ab.py ell p r s [nmax]
"""
import sys, time, flint
ell, p, r, s = map(int, sys.argv[1:5])
lam = (ell-3)//2
Bd = (lam+1)*4096
nmax = int(sys.argv[5]) if len(sys.argv) > 5 else Bd
def leg(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1
def pent(d, D):
    """prod_{n>=1} (1 - q^(dn)) mod ell to O(q^D): Euler's pentagonal number theorem (sparse construction)."""
    P = flint.nmod_poly([1], ell)
    k = 1
    while d*k*(3*k-1)//2 < D:
        sg = 1 if k % 2 == 0 else ell-1
        P[d*k*(3*k-1)//2] = sg
        if d*k*(3*k+1)//2 < D: P[d*k*(3*k+1)//2] = sg
        k += 1
    return P
t0 = time.time()
D = ell*p*p*nmax + 1
f1 = pent(1, D)
g = f1.mul_low(f1, D); del f1
g = g.inverse_series_trunc(D)
f2 = pent(2, D)
pb = g.mul_low(f2, D); del g, f2
print(f'pbar mod {ell} computed for n < {D} in {time.time()-t0:.0f}s; first values', [int(pb[i]) for i in range(8)], flush=True)
assert [int(pb[i]) for i in range(8)] == [x % ell for x in (1, 2, 4, 8, 14, 24, 40, 64)]
c0 = lambda n: (int(pb[ell*n]) * (1 if n % 2 == 0 else -1)) % ell
eps = leg(-1, p)**lam
c1 = pow(p, lam-1, ell); c2 = pow(p, 2*lam-1, ell)
hfile = [int(x) for x in open(f'h_{ell}.txt').read().strip().strip('[]').split(',')]
assert len(hfile) == Bd + 1
mism = 0; nz = 0; cnt = 0; congbad = 0; congcnt = 0
for n in range(nmax+1):
    v = c0(p*p*n) + eps*leg(n, p)*c1*c0(n) - s*eps*c1*c0(n)
    if n % (p*p) == 0: v += c2*c0(n//(p*p))
    v %= ell
    if v != hfile[n]: mism += 1
    if n % 8 == r:
        cnt += 1
        if v: nz += 1
        if leg(n, p) == s:
            congcnt += 1
            if int(pb[ell*p*p*n]) != 0: congbad += 1
import resource
print('peak RSS (GB):', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20, 2))
print(f'ell={ell}: n <= {nmax} (Sturm bound {Bd}): mismatches hd(n) vs basis-expansion h(n): {mism}; '
      f'n = {r} mod 8: {cnt} values, nonzero hd: {nz}; direct congruence pbar({ell}*{p}^2 n): {congcnt} admissible n, failures: {congbad}')
