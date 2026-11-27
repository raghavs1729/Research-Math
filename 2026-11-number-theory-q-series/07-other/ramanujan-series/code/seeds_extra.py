"""Seed library: Au's seeds (wz.SEEDS) + classical summations at z != 1 + Ebisu-type one-parameter seeds."""
import json, os
from fractions import Fraction as Fr
from wz import *

# --- Dixon6F5 (Au, arXiv:2602.08721 Sec 2.1 (7))
SEEDS['Dixon6F5'] = (('a', 'b', 'c'), Term('2*a-b-c+2*k+1', [(2, '4*a-2*b-2*c'), (-1, 'b+c+k')], [
    ('a+1/2', 1), ('a-b+1', 1), ('a-c+1', 1), ('-b+k+1', 1), ('b+k', 1), ('-c+k+1', 1), ('c+k', 1),
    ('a-b-c+3/2', 1), ('2*a-b-c+1', 1), ('2*a-b-c+k+1', 1),
    ('k+1', -1), ('2*a-b+k+1', -1), ('2*a-c+k+1', -1), ('2*a-b-2*c+k+2', -1), ('2*a-2*b-c+k+2', -1)]))

# --- classical two/three-parameter summations at z = -1, 1/2 (not in Au's list)
# Kummer: 2F1(2a,b;1+2a-b;-1) = Gamma(1+2a-b)Gamma(1+a)/(Gamma(1+2a)Gamma(1+a-b))
SEEDS['Kummer2F1'] = (('a', 'b'), Term('1', [(-1, 'k')], [
    ('2*a+k', 1), ('b+k', 1), ('k+1', -1), ('1+2*a-b+k', -1),
    ('2*a', -1), ('b', -1), ('1+2*a', 1), ('1+a-b', 1), ('1+a', -1)]))
# Gauss second: 2F1(2a,2b;a+b+1/2;1/2) = Gamma(1/2)Gamma(a+b+1/2)/(Gamma(a+1/2)Gamma(b+1/2))
SEEDS['Gauss2nd'] = (('a', 'b'), Term('1', [(Fr(1, 2), 'k')], [
    ('2*a+k', 1), ('2*b+k', 1), ('k+1', -1), ('a+b+1/2+k', -1),
    ('2*a', -1), ('2*b', -1), ('a+1/2', 1), ('b+1/2', 1)]))
# Bailey: 2F1(2a,1-2a;2c;1/2) = Gamma(c)Gamma(c+1/2)/(Gamma(c+a)Gamma(c-a+1/2))
SEEDS['Bailey2F1'] = (('a', 'c'), Term('1', [(Fr(1, 2), 'k')], [
    ('2*a+k', 1), ('1-2*a+k', 1), ('k+1', -1), ('2*c+k', -1),
    ('2*a', -1), ('1-2*a', -1), ('2*c', 1), ('c', -1), ('c+1/2', -1), ('c+a', 1), ('c-a+1/2', 1)]))
# well-poised 4F3(-1): sum (a+2k)(a)_k(b)_k(c)_k (-1)^k/(k!(1+a-b)_k(1+a-c)_k) * 1/a = Gamma(1+a-b)Gamma(1+a-c)/(Gamma(1+a)Gamma(1+a-b-c))
SEEDS['WP4F3m1'] = (('a', 'b', 'c'), Term('a+2*k', [(-1, 'k')], [
    ('a+k', 1), ('b+k', 1), ('c+k', 1), ('k+1', -1), ('1+a-b+k', -1), ('1+a-c+k', -1),
    ('b', -1), ('c', -1), ('1+a-b-c', 1)]))
# Whipple 3F2(a,1-a,c;e,1+2c-e;1) = pi 2^{1-2c} Gamma(e)Gamma(1+2c-e)/(Gamma((a+e)/2)Gamma((a+1+2c-e)/2)Gamma((1-a+e)/2)Gamma((2+2c-a-e)/2))
# rescaled a->2a, e->2e:
SEEDS['Whipple3F2'] = (('a', 'c', 'e'), Term('1', [(2, '2*c')], [
    ('2*a+k', 1), ('1-2*a+k', 1), ('c+k', 1), ('k+1', -1), ('2*e+k', -1), ('1+2*c-2*e+k', -1),
    ('2*a', -1), ('1-2*a', -1), ('c', -1),
    ('a+e', 1), ('a+c-e+1/2', 1), ('e-a+1/2', 1), ('1+c-a-e', 1)]))


def load_ebisu(path=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ebisu_seeds.json')):
    out = {}
    for i, e in enumerate(json.load(open(path))):
        s, T = auto_seed_2F1(e['A'], e['B'], e['C'], e['x'])
        out['Eb%03d[%s|%s,%s;%s;%s]' % (i, e['label'], e['A'], e['B'], e['C'], e['x'])] = (('a',), T)
    return out


if __name__ == '__main__':
    # verify the seed criterion for the hand-coded extras
    for name in ['Dixon6F5', 'Kummer2F1', 'Gauss2nd', 'Bailey2F1', 'WP4F3m1', 'Whipple3F2']:
        params, T = SEEDS[name]
        res = []
        for p in params:
            F = T.subs({p: p + '+n'})
            try:
                R = wz_mate(F)
                res.append(check_wz(F, R))
            except NoMate as ex:
                res.append('NoMate ' + str(ex))
        print(name, res, flush=True)
