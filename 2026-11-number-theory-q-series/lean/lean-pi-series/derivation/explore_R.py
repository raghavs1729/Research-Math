"""Explore the calibration hits for Ramanujan's (6n+1)/2^8n [C1] and (42n+5)/2^12n [C2] series:
for each WZ pair found by the scan, print sigma, rho, R and the boundary behaviour."""
import sys, json
sys.path.insert(0, '../../ramanujan-series/code')
import mpmath as mp
from fractions import Fraction as Fr
from wz import *
from seeds_extra import SEEDS, load_ebisu
from numerics import WZNum
mp.mp.dps = 30
ALL = dict(SEEDS)
def make_F(name, shifts, K, base, k0):
    params, T = ALL[name]
    m = {}
    for i, p in enumerate(params):
        m[p] = lf_add({'1': Fr(base[i])}, {'n': Fr(shifts[i])})
    m['k'] = lf_add({'k': Fr(1), 'n': Fr(K)}, {'1': Fr(k0)})
    return T.subs(m)
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
T = {'C1': lambda k: (6*k+1)*mp.binomial(2*k,k)**3/mp.mpf(256)**k, 'C2': lambda k: (42*k+5)*mp.binomial(2*k,k)**3/mp.mpf(4096)**k}
want = sys.argv[1:]
for line in open('../../ramanujan-series/code/match_calib.out.txt'):
    if not line.startswith('{'): continue
    h = json.loads(line)
    if h['id'] not in want or h['j'] != 0: continue
    print('=' * 100); print(h)
    try:
        F = make_F(h['seed'], h['shifts'], h['K'], h['base'], h['k0'])
        R = wz_mate(F); assert check_wz(F, R)
        print('F =', F)
        print('sigma =', show(F.ratio('k'))); print('rho =', show(F.ratio('n'))); print('R =', show(R))
        W = WZNum(F, R); t = T[h['id']]
        def sf(f, *a):
            try: return mp.nstr(f(*a), 10)
            except Exception as e: return 'ERR'
        print('G(n,0)/t(n):', [sf(lambda n: W.Gv(n,0)/t(n), n) for n in range(1,4)])
        print('F(0,k):', [sf(W.Fv, 0, k) for k in range(4)], 'F(1,k):', [sf(W.Fv, 1, k) for k in range(4)])
        for n in (0,1,2): print(' G(%d,K): K=1e3,1e5,1e7:'%n, [sf(W.Gv, n, K) for K in (10**3,10**5,10**7)])
        for n in (1,2,5,10): print(' F(%d,k)/F(%d,0) k=1,10,100,1000:'%(n,n), [sf(lambda k: W.Fv(n,k)/W.Fv(n,0), k) for k in (1,10,100,1000)])
    except Exception as e:
        print('FAILED', repr(e))
