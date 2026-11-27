from mpmath import mp, mpf, quad, cosh, exp, inf, pi, log, sqrt, gamma, agm
mp.dps=30
def pcf(a,b,N):
    # a_0 + b_1/(a_1 + b_2/(a_2+...)) by forward recurrence
    p0,p1=mpf(1),mpf(a(0)); q0,q1=mpf(0),mpf(1)
    for n in range(1,N+1):
        p0,p1=p1,a(n)*p1+b(n)*p0; q0,q1=q1,a(n)*q1+b(n)*q0
        if n%64==0:
            s=abs(q1); p0/=s;p1/=s;q0/=s;q1/=s
    return p1/q1
def SR(z,m): return quad(lambda t: exp(-z*t)/cosh(t)**m,[0,inf])
for (i,j,k) in [(0,0,0),(0,0,1),(1,0,0),(2,1,1),(3,2,0)]:
    z=i+2*j+1; m=i+k+1
    v=[pcf(lambda n:z, lambda n:n*n+(i+k)*n, N) for N in (20000,20001,80000,80001)]
    print(i,j,k, v, 1/SR(z,m))
# eq (8): 4 + 6/(4+ ...), b_n = 4n^2+2n ; value 2G/(4G-3)?  Gauss constant
Gg=1/agm(1,sqrt(2))
L=pi*Gg   # lemniscate constant varpi = pi*G
print('Gauss',Gg,'Lemn',L)
v=[pcf(lambda n:4, lambda n:4*n*n+2*n, N) for N in (40000,40001)]
print(v, 2*Gg/(4*Gg-3), 1/SR(2,1.5)*2)
v=[pcf(lambda n:4, lambda n:4*n*n+2*n-2, N) for N in (40000,40001)]
print(v, (4*Gg-1)/(3*Gg-2))
v=[pcf(lambda n:4, lambda n:4*n*n-2*n, N) for N in (40000,40001)]
print(v, -6/(L-4), 7/(L-4)**2, 1/SR(2,0.5)*2)
