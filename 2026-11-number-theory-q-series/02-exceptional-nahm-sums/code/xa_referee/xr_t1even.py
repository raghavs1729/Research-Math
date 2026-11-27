"""Consistency of Corollary 5 with T1-A-even.md Theorem A / (0.3) for r = 2, 4, 6:
   F_{T_1,r+1}(q) (direct)  ==  sum_N q^{-2N^2/(2k'+1)} * sum_n mult_{k' L0}^{A_1^(1)}(k'L0 + N alpha - n delta) q^n,  k' = r/2
and the total of the A_2^(2) side: sum over sectors m of q^{-m^2/(2k)} sum_n mult^{A_2^(2)}(kL0 + m eps - n delta) q^n (q)_inf."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
from xr_lib import Affine, Mults, nahm_sectors, inv_poch_series

def addser(d, e, coeffs):
    for i, c in enumerate(coeffs):
        if c: d[e + i] += c

def nahm_total(k, Nmax):
    secs = nahm_sectors([[1]], [1], k, 2 * Nmax, 2)
    F = defaultdict(int)
    for c, lst in secs.items():
        for p, E0, mu in lst:
            if E0 > Nmax: continue
            s = np.zeros(Nmax + 1, dtype=object); s[0] = 1
            for n in p:
                if n: s = np.convolve(s, inv_poch_series(1, n, Nmax, 1))[:Nmax + 1]
            addser(F, E0, s)
    return F

for r in (2, 4, 6):
    k = r + 1; kp = r // 2
    Nmax = {2: 14, 4: 12, 6: 10}[r]
    F = nahm_total(k, Nmax)
    A1 = Affine('A', 1)
    M = Mults(A1, kp, Nmax + 2 * 400 // (2 * kp + 1) + 2)
    R = defaultdict(int)
    Nb = 0
    for N in range(-20, 21):
        sh = Fr(-2 * N * N, 2 * kp + 1)
        st = M.string((N,))
        # mult(k'L0 + N alpha - n delta) nonzero only for n >= N^2/k'
        addser(R, sh, st)
    # A_2^(2) side
    T = Affine('T', 1); MT = Mults(T, k, 2 * Nmax + 2 * k + 2)
    H = T.heis_series(2 * Nmax + 2 * k + 2)
    S = defaultdict(int)
    for m in range(-3 * k, 3 * k + 1):
        if m % k != m % k: pass
    for m in range(k):  # one rep per class
        st = np.convolve(np.array(MT.string((m,)), dtype=object), H)[:2 * Nmax + 2 * k + 3]
        for i, c in enumerate(st):
            if c: S[Fr(i, 2) - Fr(m * m, 2 * k)] += c
    keys = sorted(e for e in set(F) | set(R) | set(S) if e <= Nmax)
    bad1 = [e for e in keys if F.get(e, 0) != R.get(e, 0)]
    bad2 = [e for e in keys if F.get(e, 0) != S.get(e, 0)]
    print(f"r={r} (k={k}, k'={kp}): exponents <= {Nmax}: {len(keys)} distinct; Nahm vs sl2-string form mismatches: {len(bad1)}; "
          f"Nahm vs A_2^(2) string form mismatches: {len(bad2)}; first coefficients {[(str(e), F[e]) for e in keys[:6]]}")
