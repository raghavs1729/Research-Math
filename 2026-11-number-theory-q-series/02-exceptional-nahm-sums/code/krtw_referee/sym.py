import sympy as sp
u, v, x, q, al = sp.symbols('u v x q alpha')
lhs = u*(1-x*u**2)*(1-x**2*v**4) + v*(1-x*v**2)*(1-x**2*u**4)
rhs = (u+v)*(1+x*u*v)*(1-x*u**2)*(1-x*v**2)
print("E identity:", sp.expand(lhs-rhs) == 0)
# GOW L^(0)_k term at n=2, y=(alpha,-alpha) vs L_k term, symbolic in alpha,q for small (n,m), k = 1,2
def poch(a, n): return sp.prod([1 - a*q**i for i in range(n)])
def DeltaC(y1, y2): return (1-y1**2)*(1-y2**2)*(y1-y2)*(y1*y2-1)
ok = True
for k in (1, 2):
    for n in range(0, 3):
        for m in range(0, 3):
            y = [al, -al]; r = [n, m]
            t = DeltaC(y[0]*q**n, y[1]*q**m)/DeltaC(y[0], y[1])
            for i in range(2):
                t *= y[i]**(2*(k+1)*r[i]) * q**((k+1)*r[i]**2 + 2*sp.binomial(r[i], 2))
                for j in range(2):
                    t *= (-y[i]/y[j])**r[i] * poch(y[i]*y[j], r[i]) / poch(q*y[i]/y[j], r[i])
            X = al**2
            Lt = sp.Rational(1, 2)*(-1)**(n+m) * X**((k+1)*(n+m)) * q**((k+2)*(n*n+m*m)-n-m) * poch(X**2*0+X**2, 0)
            Q = q**2
            def pQ(a, nn): return sp.prod([1 - a*Q**i for i in range(nn)])
            Lt = sp.Rational(1, 2)*(-1)**(n+m) * X**((k+1)*(n+m)) * q**((k+2)*(n*n+m*m)-n-m) * pQ(X**2, n)*pQ(X**2, m) \
                 * (q**n+q**m)*(1+X*q**(n+m))*(1-X*q**(2*n))*(1-X*q**(2*m)) / (pQ(Q, n)*pQ(Q, m)*(1-X)*(1-X**2))
            ok &= sp.simplify(sp.factor(t - Lt)) == 0
print("L^(0)_k(alpha,-alpha) term == L_k term (k=1,2; n,m<=2):", ok)
