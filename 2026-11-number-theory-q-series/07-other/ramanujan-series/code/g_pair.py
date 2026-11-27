"""print the data of one WZ pair: python3 g_pair.py SEED 'shifts' K 'base' k0"""
import sys, json
from explore import *
name, shifts, K, base, k0 = sys.argv[1], json.loads(sys.argv[2]), int(sys.argv[3]), json.loads(sys.argv[4]), sys.argv[5]
full = [n for n in ALL if n == name or n.startswith(name)][0]
F = make_F(full, shifts, K, base, k0)
print('F =', F)
R = wz_mate(F); assert check_wz(F, R)
sg, rh = F.ratio('k'), F.ratio('n')
fa = lambda x: str(pari.factor(pari.numerator(x))).replace('\n', ' ') + '  /  ' + str(pari.factor(pari.denominator(x))).replace('\n', ' ')
print('sigma =', fa(sg)); print('rho =', fa(rh)); print('R =', fa(R))
print('sigma(0,k) =', fa(pari.subst(sg, 'n', 0)))
print('rho(n,0) =', fa(pari.subst(rh, 'k', 0)))
print('R(n,0) =', fa(pari.subst(R, 'k', 0)))
print('R(0,k) =', fa(pari.subst(R, 'n', 0)))
