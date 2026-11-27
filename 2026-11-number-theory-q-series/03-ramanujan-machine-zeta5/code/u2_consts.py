from fractions import Fraction as Fr
from mpmath import mp, mpf, pi, gamma, zeta, pslq, exp, gammainc, sqrt, log, hyper
import asym_match as am
mp.dps = 320
cases = {'U1': (Fr(-15, 8), Fr(9, 8), Fr(-3, 2), Fr(-5, 2)),
         'U2': (Fr(-3, 2), Fr(7, 8), Fr(-1), Fr(-3)),
         'U3': (Fr(5, 2), Fr(-5, 8), Fr(1), Fr(-5))}
out = {}
for name, (c3, c1, d1, d2) in cases.items():
    a = lambda n, c3=c3, c1=c1: 8 * (n**5 + (n + 1)**5) + int(8 * c3) * (n**3 + (n + 1)**3) + int(8 * c1) * (2 * n + 1)
    L, (al1, be1, al2, be2), errs = am.limit(5, 8, am.at5(c3, c1), d1, d2, a, 170, 220)
    out[name] = (L, al1, be1, al2, be2)
    print(name, 'c_w =', mp.nstr(al1, 60), ' al2 =', mp.nstr(al2, 40), ' be2=', mp.nstr(be2, 40), flush=True)
import pickle
pickle.dump({k: [str(x) for x in v] for k, v in out.items()}, open('u_consts.pkl', 'wb'))
# U1 sanity: c_w = 4 A^3 / Gamma(1/4)^2 with A = sqrt(pi)/Gamma(3/4)^2  (w_n ~ c_w n^{-3/2}; series in t=n+1/2)
A = sqrt(pi) / gamma(mpf(3) / 4)**2
print('U1 predicted c_w', mp.nstr(4 * A**3 / gamma(mpf(1) / 4)**2, 60))
