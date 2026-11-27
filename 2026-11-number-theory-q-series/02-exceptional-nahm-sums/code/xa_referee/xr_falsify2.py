"""Correct Lie data, wrong Nahm data (or vice versa): sector identity must fail."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fractions import Fraction as Fr
from xr_lib import Affine, Mults, nahm_sectors, inv_poch_series, ip

def compare(Alie, k, N, G, D):
    M = Mults(Alie, k, N); H = Alie.heis_series(N); DU = Alie.DU; l = Alie.l
    secs = nahm_sectors(G, D, k, N, DU)
    bad = 0; tot = 0
    for c in itertools.product(range(k), repeat=l):
        mu = tuple(c)
        lhs = np.convolve(np.array(M.string(mu), dtype=object), H)[:N+1]
        rhs = np.zeros(N+1, dtype=object)
        sh = Fr(ip(G, mu, mu), 2*k)
        for p, E0, _ in secs.get(c, []):
            e = (sh+E0)*DU
            if e.denominator != 1: bad += 1; break
            e = int(e)
            if e > N: continue
            s = np.zeros(N+1, dtype=object); s[e] = 1
            for i in range(l):
                for t in range(1, k):
                    n = p[i*(k-1)+t-1]
                    if n: s = np.convolve(s, inv_poch_series(D[i], n, N, DU))[:N+1]
            rhs += s
        tot += 1
        if any(lhs[i] != rhs[i] for i in range(N+1)): bad += 1
    return bad, tot

C3 = Affine('C', 3); B3 = Affine('B', 3)
print("A_5^(2) strings vs (B_3,A_1) Nahm sectors (same basis order):", compare(C3, 2, 8, B3.G, B3.D))
print("D_4^(2) strings vs (C_3,A_1) Nahm sectors:", compare(B3, 2, 8, C3.G, C3.D))
print("A_5^(2) strings vs C_3 Gram with D=1 (q;q) denominators:", compare(C3, 2, 8, C3.G, [1,1,1]))
T1 = Affine('T', 1)
print("A_2^(2) k=4 strings vs Nahm with G=(2) [i.e. (A_1,A_3)]:", compare(T1, 4, 16, [[2]], [1]))
print("control: A_5^(2) strings vs (C_3,A_1):", compare(C3, 2, 8, C3.G, C3.D))
