import sympy as sp
from explore import *
import mpmath as mp
full = [n for n in ALL if n.startswith('Eb076')][0]
F = make_F(full, [-1], -1, ['0'], '0')
Rp = wz_mate(F); assert check_wz(F, Rp)
n, k = sp.symbols('n k')
sigma = sp.factor(sp.sympify(str(F.ratio('k')).replace('^', '**')))
rho = sp.factor(sp.sympify(str(F.ratio('n')).replace('^', '**')))
R = sp.factor(sp.sympify(str(Rp).replace('^', '**')))
print('sigma =', sigma); print('rho =', rho); print('R =', R)
assert sp.cancel(rho - 1 - (R.subs(k, k+1)*sigma - R)) == 0
print('R(n,0)*... ratio of G(n,0):', sp.factor(sp.cancel(R.subs({k: 0, n: n+1}, simultaneous=True)/R.subs(k, 0)*rho.subs(k, 0))))
tau = sp.cancel(rho.subs(k, n-1)*sigma.subs({n: n+1, k: n-1}, simultaneous=True))
gt = sp.cancel(sp.limit(sp.cancel(R*sigma.subs(k, k-1)), k, n))
D = sp.factor(sp.cancel(tau + gt))
print('tau =', sp.factor(tau)); print('G(n,n)/F(n,n-1) =', sp.factor(gt)); print('D =', D)
print('d_{n+1}/d_n =', sp.factor(sp.cancel(tau*D.subs(n, n+1)/D)))
