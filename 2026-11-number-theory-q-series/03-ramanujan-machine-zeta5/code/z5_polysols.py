import sympy as sp
n,c1=sp.symbols('n c1')
def allowed(d, s=5):
    # f polynomial of degree d, monic; solve (n+1)^s f(n+1) + n^s f(n-1) - a(n) f(n) = 0
    cs=sp.symbols('f0:%d'%d)
    f=n**d+sum(cs[i]*n**i for i in range(d))
    c3=sp.Rational(d*(d+s-1),2)
    a=n**s+(n+1)**s+c3*(n**3+(n+1)**3)+c1*(2*n+1)
    expr=sp.expand((n+1)**s*f.subs(n,n+1)+n**s*f.subs(n,n-1)-a*f)
    eqs=sp.Poly(expr,n).all_coeffs()
    eqs=[e for e in eqs if e!=0]
    G=sp.groebner(eqs,*cs,c1,order='lex')
    # last polynomial in c1 only
    polys=[g for g in G.exprs if g.free_symbols<= {c1}]
    return c3,[sp.factor(p) for p in polys]
for d in range(0,7):
    print(d,allowed(d))
