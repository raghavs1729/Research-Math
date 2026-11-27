from fractions import Fraction as Fr
from math import comb
N=14
def mul(a,b):
    c=[Fr(0)]*(N+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>N: break
            c[i+j]+=x*y
    return c
def theta(a): return [k*x for k,x in enumerate(a)]
def fit(w):
    # n^5 w_n + (n-1)^5 w_{n-2} = 2 t (t^4 + al t^2 + be) w_{n-1},  t = n-1/2 ; solve al,be from n=1,2 and test
    import sympy as sp
    al,be=sp.symbols('al be')
    eqs=[]
    for n in (1,2,3,4,5,6,7,8):
        t=sp.Rational(2*n-1,2)
        lhs=n**5*sp.Rational(w[n].numerator,w[n].denominator)+((n-1)**5*sp.Rational(w[n-2].numerator,w[n-2].denominator) if n>=2 else 0)
        eqs.append(sp.Eq(lhs, 2*t*(t**4+al*t**2+be)*sp.Rational(w[n-1].numerator,w[n-1].denominator)))
    sol=sp.solve(eqs[:2],[al,be],dict=True)
    ok=[all(sp.simplify(e.lhs-e.rhs).subs(s)==0 for e in eqs) for s in sol]
    return sol,ok
H=lambda n: sum(Fr(1,i) for i in range(1,n+1))
# U1 candidate: 3F2(1/2,1/2,1/2;1,1;x)^2
A3=[Fr(comb(2*k,k)**3,64**k) for k in range(N+1)]
print('3F2^2', fit(mul(A3,A3)))
# 4F3(1/2^4;1^3): y0, g ; psi = y0^2 + y0*theta g - g*theta y0
A4=[Fr(comb(2*k,k)**4,256**k) for k in range(N+1)]
g4=[A4[k]*4*(H(2*k)-H(k))*2 for k in range(N+1)]   # d/de of (1/2+e)_k^4/(1+e)_k^4 = 4*(sum 1/(1/2+j) - sum 1/(1+j)) = 4*(2H_{2k}-H_k - H_k)=8(H_2k-H_k)
psi=[x+y-z for x,y,z in zip(mul(A4,A4),mul(A4,theta(g4)),mul(g4,theta(A4)))]
print('wedge2 4F3', fit(psi))
# other candidates: 2F1(1/2,1/2;1;x)^2 * something? products
A2=[Fr(comb(2*k,k)**2,16**k) for k in range(N+1)]
print('2F1(1/2,1/2)^4?', fit(mul(mul(A2,A2),mul(A2,A2))))
print('A2*A4?', fit(mul(A2,A4)))
A5=[Fr(comb(2*k,k)**5,1024**k) for k in range(N+1)]
print('A5', fit(A5))
print('A2^2*A3', fit(mul(mul(A2,A2),A3)))
