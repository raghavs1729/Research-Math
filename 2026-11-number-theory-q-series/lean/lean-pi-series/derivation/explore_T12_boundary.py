"""T12 boundary sum B = sum_k psi(0,k) = log(27/16): look for a second WZ direction (shift of a in Seed4(a,b,k+1/2),
coefficient of b at b = 0) whose G(m,0) is elementary (rational times geometric)."""
import sys
sys.path.insert(0, '../../ramanujan-series/code')
from fractions import Fraction as Fr
from wz import *
import mpmath as mp
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
params, T = SEEDS['Seed4']
print('Seed4 =', T)
for A in (1, -1, 2, -2, 0):
    for K in (0, 1, -1, 2, -2):
        if A == 0 and K == 0: continue
        try:
            F = T.subs({'a': lf_add({'1': Fr(1,2)}, {'n': Fr(A)}), 'b': 'b', 'k': lf_add({'k': Fr(1), 'n': Fr(K)}, {'1': Fr(1,2)})})
            R = wz_mate(F); assert check_wz(F, R)
            sig = pari.subst(F.ratio('k'), 'b', 0); rho = pari.subst(F.ratio('n'), 'b', 0); R0 = pari.subst(R, 'b', 0)
            g = pari.subst(R0, 'k', 0)
            gr = pari.subst(g, 'n', pari('n')+1)/g*pari.subst(rho, 'k', 0)
            print('A=%d K=%d' % (A, K)); print('   rho(n,0) =', show(pari.subst(rho,'k',0))); print('   R(n,0)  =', show(g)); print('   G ratio =', show(gr))
        except Exception as e:
            print('A=%d K=%d FAILED %r' % (A, K, e))
