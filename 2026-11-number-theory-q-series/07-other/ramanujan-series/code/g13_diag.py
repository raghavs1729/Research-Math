import sympy as sp
n, k, e = sp.symbols('n k e')
sigma = (-2*k**2 + (2*n - 2)*k + (4*n**2 + 4*n))/(6*k**2 + 15*k + (-6*n**2 - 3*n + 9))
rho = ((36*n**2 + 72*n + 35)*k**2 + (36*n**2 + 72*n + 35)*k + (-36*n**4 - 108*n**3 - 107*n**2 - 35*n))/((16*n + 8)*k**3 + (-48*n**2 - 48*n - 12)*k**2 + (-48*n**2 - 64*n - 20)*k + (64*n**4 + 224*n**3 + 272*n**2 + 136*n + 24))
R = ((6*n + 3)*k**2 + (-36*n**2 - 54*n - 21)*k + (30*n**3 + 51*n**2 + 21*n))/((8*n + 4)*k**2 + (-32*n**2 - 40*n - 12)*k + (32*n**3 + 64*n**2 + 40*n + 8))
print('sigma =', sp.factor(sigma)); print('rho =', sp.factor(rho)); print('R =', sp.factor(R))
assert sp.cancel(rho - 1 - (R.subs(k, k+1)*sigma - R)) == 0
tau = sp.factor(sp.cancel(rho.subs(k, n-1)*sigma.subs({n: n+1, k: n-1}, simultaneous=True)))
gt = sp.factor(sp.cancel(sp.limit(sp.cancel(R*sigma.subs(k, k-1)), k, n)))
D = sp.factor(sp.cancel(tau + gt))
print('tau(n) = F(n+1,n)/F(n,n-1) =', tau)
print('G(n,n)/F(n,n-1) =', gt)
print('D(n) = d_n/F(n,n-1) =', D)
print('d_{n+1}/d_n =', sp.factor(sp.cancel(tau*D.subs(n, n+1)/D)))
