"""Second audit, undecided classes of Theorem 6: L_1(q) (mpmath, al1) against (2/3)N, and a float scan of
the prime-side F_{q,a}(log x) on 1 <= x <= X (formula (F), my constants), to confirm that no prime-side
witness exists in that range.  Numerical (not a certificate; it only confirms 'no witness found')."""
import sys, math
import numpy as np
import mpmath as mp
from math import gcd
from sympy import totient, primefactors
from aw1 import all_chars, constants, sieve

X = float(sys.argv[1])
P = sieve(int(X) + 10)
for spec in sys.argv[2:]:
    q, a = map(int, spec.split(','))
    phi = int(totient(q))
    N = sum(1 for b in range(q) if gcd(b, q) == 1 and (b * b - a) % q == 0)
    chars = all_chars(q)
    L1 = float(sum(ch.S() for ch in chars))
    ne, no, B, C = constants(q, a, chars)
    B, C = float(B), float(C)
    pts = []; wts = []
    for p in P:
        p = int(p); lp = math.log(p); pk = p
        mpp = q
        while mpp % p == 0:
            mpp //= p
        while pk <= X:
            if q % p != 0 and (pk - a) % q == 0:
                pts.append(pk); wts.append(phi * lp)
            elif q % p == 0 and (pk - a) % mpp == 0:
                pts.append(pk); wts.append(int(totient(mpp)) * lp)
            pk *= p
    o = np.argsort(pts); pts = np.array(pts, float)[o]; wts = np.array(wts)[o]
    grid = np.unique(np.concatenate([np.geomspace(1, X, 300000), pts, pts - 1e-9]))
    grid = grid[grid >= 1]
    # S0(x) = sum w, S1(x) = sum w*n over n <= x
    cw = np.concatenate([[0], np.cumsum(wts)]); cwn = np.concatenate([[0], np.cumsum(wts * pts)])
    idx = np.searchsorted(pts, grid, side='right')
    lin = grid * cw[idx] - cwn[idx]
    x = grid; u = 1 / x
    with np.errstate(divide='ignore', invalid='ignore'):
        Te = np.where(u < 1, (x / 2) * ((1 + u) * np.log1p(u) + (1 - u) * np.log1p(-u)), math.log(2))
        To = np.where(u < 1, 1 - 0.5 * np.log1p(-u * u) - np.arctanh(np.minimum(u, 1 - 1e-16)) / u, 1 - math.log(2))
    ell = -ne * x * np.log(x) + B * x + no * np.log(x) + C - (1 + ne) * Te - no * To
    F = (lin - x * x / 2 - ell) / x ** 1.5
    i = int(np.argmax(F))
    print(f'({q},{a}) N={N} (2/3)N={2*N/3:.4f} L_1={L1:.6f} L_1-(2/3)N={L1-2*N/3:+.4f}  max F on [1,{X:.0e}] = {F[i]:.5f} at x={x[i]:.6g}; min F = {F.min():.4f}; max|F|<=L_1: {np.abs(F).max() <= L1}')
