"""Is the U3 recurrence gauge-equivalent to the U2 (or U1, or Euler zeta(5)) recurrence?
Search P0(n) w'_n = P1(n) w_n + P2(n) w_{n+1}, polynomials of degree <= D, by exact linear algebra."""
from fractions import Fraction as Fr
import sympy as sp
def seq(al, be, N):
    # n^5 w_n + (n-1)^5 w_{n-2} = 2t(t^4+al t^2+be) w_{n-1}, t = n-1/2
    w=[Fr(1)]
    for n in range(1,N+1):
        t=Fr(2*n-1,2)
        w.append((2*t*(t**4+al*t**2+be)*w[n-1]-((n-1)**5*w[n-2] if n>=2 else 0))/n**5)
    return w
par={'U1':(Fr(5,8),Fr(1,32)),'U2':(Fr(1),Fr(1,16)),'U3':(Fr(5),Fr(25,16)),'Z5':(Fr(5,2),Fr(5,16))}
def test(src,dst,D,N=80):
    w=seq(*par[src],N+2); v=seq(*par[dst],N+2)
    rows=[]
    for n in range(0,N):
        rows.append([v[n]*n**i for i in range(D+1)]+[-w[n]*n**i for i in range(D+1)]+[-w[n+1]*n**i for i in range(D+1)])
    M=sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in r] for r in rows[:3*(D+1)+6]])
    ns=M.nullspace()
    return ns
for src,dst in (('U2','U3'),('U1','U3'),('Z5','U3'),('U1','U2')):
    for D in (2,4,6):
        ns=test(src,dst,D)
        print(src,'->',dst,'D=',D,'nullity',len(ns))
        if ns:
            n=sp.symbols('n'); vec=ns[0]; vec=vec/vec[0] if vec[0]!=0 else vec
            P=[sp.factor(sum(vec[j*(D+1)+i]*n**i for i in range(D+1))) for j in range(3)]
            print('   P0,P1,P2 =',P); break
