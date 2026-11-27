from mpmath import mp, mpf, richardson, zeta, log, pi, pslq, catalan, identify, findpoly
mp.dps=120
def convergents(a,b,N):
    p0,p1=mpf(1),mpf(a(0)); q0,q1=mpf(0),mpf(1)
    out={}
    for n in range(1,N+1):
        p0,p1=p1,a(n)*p1+b(n)*p0; q0,q1=q1,a(n)*q1+b(n)*q0
        s=abs(q1); p0/=s;p1/=s;q0/=s;q1/=s
        out[n]=p1/q1
    return out
def rich(seq, ns, expo):
    # eliminate terms n^{-e} for e in expo successively (Neville-like using sequences at n, 2n, 4n...)
    vals=[seq[n] for n in ns]
    for e in expo:
        vals=[(2**e*vals[i+1]-vals[i])/(2**e-1) for i in range(len(vals)-1)]
    return vals
A=lambda c5,c3,c1: (lambda n: c5*(n**5+(n+1)**5)+c3*(n**3+(n+1)**3)+c1*(2*n+1))
cases={'U1':(8,-15,9),'U2':(8,-12,7),'U3':(8,20,-5)}
N=2**14
ns=[2**i for i in range(3,15)]
res={}
for name,(c5,c3,c1) in cases.items():
    cv=convergents(A(c5,c3,c1), lambda n:-64*n**10, N)
    print(name,[mp.nstr(cv[n],25) for n in (100,1000,N)])
    for expo in ([1,2,3,4,5,6,7,8,9,10],[2,3,4,5,6,7,8,9,10,11],[2,4,6,8,10,12,14,16,18],[1,3,5,7,9,11,13,15,17],[6,7,8,9,10,11,12,13,14],[6,8,10,12,14,16,18,20]):
        v=rich(cv,ns,expo)
        print('  ',expo[:3],mp.nstr(v[-1],45),mp.nstr(v[-1]-v[-2],5))
