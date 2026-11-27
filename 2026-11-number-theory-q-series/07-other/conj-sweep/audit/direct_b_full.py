"""AUDIT: direct full-range test for Theorem B (and for the form h of its proof) with low memory.

pbar(23 k) mod 23 for k <= KMAX = 169*45056 is obtained WITHOUT modular forms from the elementary identity
    1/G(q) = G(q)^22 / G(q)^23 = G(q)^22 / G(q^23)  (mod 23),   G(q) = f_1^2/f_2 = sum_n (-1)^n q^(n^2),
so   sum_k pbar(23 k) q^k = U_23(G^22) * (1/G)   (mod 23),
and U_23(G^22) = U_23(G^16 * G^6) = sum_{r mod 23} X_r * Y_{-r} (23-dissections X_r, Y_r of G^16, G^6;
G^16 and G^6 are computed to q^(23*KMAX) by chunked multiplication to keep the temporary memory small).
The result is cross-checked against a plain series inversion on an initial range, then used exactly as in
direct_ab.py:   c0(n) = (-1)^n pbar(23 n),
    hd(n) = c0(169 n) + (n/13) 13^9 c0(n) + 13^19 c0(n/169) - 13^9 c0(n)   (mod 23)       [(-1/13)^10 = +1, s = +1]
compared with the q-expansion h written by ab.gp for ALL n <= 45056, and the congruence pbar(23*169*n) = 0 for
n = 5 (8), (n/13) = +1.
"""
import sys, time, flint, resource
ell, p, r, s = 23, 13, 5, 1
lam = (ell-3)//2; Bd = (lam+1)*4096
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else Bd
KMAX = p*p*nmax; L = KMAX + 1; D = ell*KMAX + 1
def leg(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1
def theta_alt(D):
    G = flint.nmod_poly([1], ell); n = 1
    while n*n < D: G[n*n] = 2 if n % 2 == 0 else ell-2; n += 1
    return G
def pent(d, D):
    P = flint.nmod_poly([1], ell); k = 1
    while d*k*(3*k-1)//2 < D:
        sg = 1 if k % 2 == 0 else ell-1
        P[d*k*(3*k-1)//2] = sg
        if d*k*(3*k+1)//2 < D: P[d*k*(3*k+1)//2] = sg
        k += 1
    return P
t0 = time.time()
# G = f_1^2/f_2 = sum (-1)^n q^(n^2): checked on an initial range
D0 = 200000
f1 = pent(1, D0); assert f1.mul_low(f1, D0) == theta_alt(D0).mul_low(pent(2, D0), D0)
# Long products are done in NCH chunks (schoolbook on chunks) to keep FLINT's temporary memory small:
# a series of length D is stored as a list of NCH polynomials of length M.
NCH = 4; M = -(-D // NCH)
def chunks(A): return [A.right_shift(i*M).truncate(M) for i in range(NCH)]
def mulc(Ac, Bc):
    """chunked product modulo q^(NCH*M)"""
    Rc = []; carry = flint.nmod_poly([0], ell)
    for k in range(NCH):
        S = flint.nmod_poly([0], ell)
        for i in range(k + 1): S = S + Ac[i]*Bc[k-i]          # diagonal k: degree < 2M-1
        S = S + carry
        Rc.append(S.truncate(M)); carry = S.right_shift(M)
    return Rc
def dissect(Xc):
    out = []
    for rr in range(ell):
        out.append(flint.nmod_poly([int(Xc[i // M][i % M]) for i in range(rr, D, ell)], ell))
    return out
Gc = chunks(theta_alt(D))
G2c = mulc(Gc, Gc); del Gc
G4c = mulc(G2c, G2c)
G6c = mulc(G4c, G2c); del G2c
print(f'G^6: {time.time()-t0:.0f}s', flush=True)
Y = dissect(G6c); del G6c
print(f'G^6 dissected: {time.time()-t0:.0f}s', flush=True)
G8c = mulc(G4c, G4c); del G4c
G16c = mulc(G8c, G8c); del G8c
X = dissect(G16c); del G16c
print(f'G^16 dissected: {time.time()-t0:.0f}s', flush=True)
xpoly = flint.nmod_poly([0, 1], ell)
U = X[0].mul_low(Y[0], L)
for rr in range(1, ell):
    U = U + xpoly*X[rr].mul_low(Y[ell-rr], L)      # indices 23i+rr and 23j+(23-rr) add up to 23(i+j+1)
del X, Y
Hinv = theta_alt(L).inverse_series_trunc(L)
PB = U.mul_low(Hinv, L)                      # PB[k] = pbar(23 k) mod 23
print(f'pbar(23k) mod 23 for k <= {KMAX}: {time.time()-t0:.0f}s', flush=True)
# cross-check against plain inversion of f_1^2/f_2 on an initial range
D1 = 3000000
pb = pent(1, D1); pb = pb.mul_low(pb, D1).inverse_series_trunc(D1).mul_low(pent(2, D1), D1)
assert [int(pb[i]) for i in range(8)] == [x % ell for x in (1, 2, 4, 8, 14, 24, 40, 64)]
kchk = min((D1-1)//ell, KMAX)
assert all(int(PB[k]) == int(pb[ell*k]) for k in range(kchk + 1))
print(f'cross-check PB[k] = pbar(23k) for k <= {kchk} (plain inversion): OK', flush=True)
c0 = lambda n: (int(PB[n]) * (1 if n % 2 == 0 else -1)) % ell
eps = leg(-1, p)**lam; c1 = pow(p, lam-1, ell); c2 = pow(p, 2*lam-1, ell)
hfile = [int(x) for x in open(f'h_{ell}.txt').read().strip().strip('[]').split(',')]
assert len(hfile) == Bd + 1
mism = nz = cnt = congbad = congcnt = 0
for n in range(nmax + 1):
    v = c0(p*p*n) + eps*leg(n, p)*c1*c0(n) - s*eps*c1*c0(n)
    if n % (p*p) == 0: v += c2*c0(n//(p*p))
    v %= ell
    if v != hfile[n]: mism += 1
    if n % 8 == r:
        cnt += 1
        if v: nz += 1
        if leg(n, p) == s:
            congcnt += 1
            if int(PB[p*p*n]) != 0: congbad += 1
print('peak RSS (GB):', round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20, 2))
print(f'ell={ell}: n <= {nmax} (Sturm bound {Bd}): mismatches hd(n) vs basis-expansion h(n): {mism}; '
      f'n = {r} mod 8: {cnt} values, nonzero hd: {nz}; direct congruence pbar({ell}*{p}^2 n): {congcnt} admissible n, failures: {congbad}')
