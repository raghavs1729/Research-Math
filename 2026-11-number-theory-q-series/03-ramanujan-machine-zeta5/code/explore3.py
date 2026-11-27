from mpmath import mp, mpf, sqrt
mp.dps=30
def pcf(a,b,N):
    p0,p1=mpf(1),mpf(a(0)); q0,q1=mpf(0),mpf(1)
    out=[]
    for n in range(1,N+1):
        p0,p1=p1,a(n)*p1+b(n)*p0; q0,q1=q1,a(n)*q1+b(n)*q0
        s=abs(q1) if q1!=0 else 1; p0/=s;p1/=s;q0/=s;q1/=s
    return p1/q1 if q1!=0 else None
for j in range(-3,5):
    for k in range(-3,5):
        D=j*j+4*k
        if D<0: continue
        r1=(j+sqrt(D))/2; r2=(j-sqrt(D))/2
        v=[pcf(lambda n:-2*n+j-1, lambda n:-n*n+j*n+k, N) for N in (2000,8000)]
        print(j,k,'roots',mp.nstr(r1,8),mp.nstr(r2,8),'pcf',mp.nstr(v[0],12),mp.nstr(v[1],12))
