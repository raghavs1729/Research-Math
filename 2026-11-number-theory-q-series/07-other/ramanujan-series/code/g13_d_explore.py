import sys, json, mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 30
eps = mp.mpf(10) ** -15
for (seed, sh, K, base, k0) in [('Eb076', [-1], -1, ['0'], '0'), ('Eb076', [-1], -1, ['1'], '1'), ('Eb070', [1], -1, ['1/2'], '1/2'), ('Eb011', [-1], 1, ['0'], '1'), ('Eb027', [-1], -2, ['0'], '0')]:
    full = [n for n in ALL if n.startswith(seed)][0]
    F = make_F(full, sh, K, base, k0)
    R = wz_mate(F); assert check_wz(F, R)
    W = WZNum(F, R)
    def safe(f):
        try: return f()
        except Exception as ex: return mp.nan
    print(seed, sh, K, base, k0)
    print('  F =', str(F)[:300])
    print('  sigma(0,k) =', pari.factor(pari.subst(F.ratio('k'), 'n', 0)), ' R(n,0) =', pari.factor(pari.subst(R, 'k', 0)))
    print('  F(0,0..3) =', [mp.nstr(safe(lambda j=j: W.Fv(eps, j)), 8) for j in range(4)])
    print('  G(0..3,0) =', [mp.nstr(safe(lambda j=j: W.Gv(j + eps, 0)), 8) for j in range(4)])
    print('  SG0 =', mp.nstr(safe(lambda: hsum_n(W, eps, 0)), 15), ' SF0 =', mp.nstr(safe(lambda: hsum_k(W, eps)), 15), ' SF8 =', mp.nstr(safe(lambda: hsum_k(W, 8 + eps)), 8))
