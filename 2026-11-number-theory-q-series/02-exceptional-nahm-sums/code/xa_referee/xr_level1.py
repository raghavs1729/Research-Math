"""Sanity check of the Freudenthal engine: level 1, sum_n mult(L0+mu-n delta)q^n prod(1-q^j)^{mult} = q^{|mu|^2/2}
(one-dimensional vacuum spaces; Frenkel-Kac and twisted analogues) and A_1^(1) level 1: mult(L0 - n delta) = p(n)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from xr_lib import Affine, Mults, ip
for X, l, N in [('A',1,30),('A',2,20),('D',4,10),('E',6,6),('C',3,12),('B',3,12),('F',4,8),('G',2,16),('T',1,40),('T',2,24),('T',3,14)]:
    A = Affine(X, l); M = Mults(A, 1, N); H = A.heis_series(N)
    bad = 0; tot = 0
    for mu in itertools.product(range(-2, 3), repeat=l):
        e = A.DU * ip(A.G, mu, mu)
        if e % 2: bad += 1; continue
        e //= 2
        lhs = np.convolve(np.array(M.string(mu), dtype=object), H)[:N+1]
        rhs = [1 if i == e else 0 for i in range(N+1)]
        tot += 1
        if list(lhs) != rhs: bad += 1
    print(f"{X}_{l} level 1, depth {N}/{A.DU}: {tot} weights, failures {bad}")
M = Mults(Affine('A',1), 1, 30)
from sympy import npartitions
print("A_1^(1) level 1 mult(L0 - n delta) == p(n), n<=30:", all(M.mult((0,), n) == npartitions(n) for n in range(31)))
