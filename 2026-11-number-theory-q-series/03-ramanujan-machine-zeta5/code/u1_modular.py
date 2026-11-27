from mpmath import mp, mpf, quad, pi, zeta, pslq, ellipk, exp, inf, gamma, qfrom, jtheta
mp.dps=60
exec(open('z5_values.py').read())
# check U1 = 3 / int_0^1 K^3 (2K' - K)(1-2z) dz
f=lambda z: ellipk(z)**3*(2*ellipk(1-z)-ellipk(z))*(1-2*z)
M=quad(f,[0,mpf(1)/100,mpf(1)/2,mpf(99)/100,1])
print('3/M',3/M); print('U1 ',U1)
# eta(2tau)^12 coefficients
N=400
# eta(2tau)^12 = q * prod (1-q^{2n})^12
coef=[0]*(N+1); coef[0]=1
for n in range(1,N//2+1):
    for _ in range(12):
        for k in range(N,2*n-1,-1):
            coef[k]-=coef[k-2*n]
a=[0]+coef[:N]   # a[n], n>=1 : coefficient of q^n in eta(2tau)^12
def Lf(s):
    # completed: Lambda(s) = (2pi/sqrt(4))^{-s} Gamma(s) L(s) = int_0^inf f(it/2) t^{s-1} dt,  f(i t/2)= sum a_n e^{-pi n t}
    # symmetric under t->1/t (level 4, weight 6): Lambda(s) = eps*Lambda(6-s).  compute directly by splitting at 1
    F=lambda t: sum(a[n]*exp(-pi*n*t) for n in range(1,N+1) if a[n])
    I1=quad(lambda t: F(t)*(t**(s-1)+t**(6-s-1)),[1,2,4,8,16])   # eps=+1 assumed, check below
    return I1*pi**s/gamma(s)
for s in (1,2,3,4,5):
    print(s, Lf(s))
L1,L2=Lf(1),Lf(2)
print(pslq([mpf(1), U1*zeta(5), U1*pi**5, U1*pi**3*L2, U1*pi**4*L1],maxcoeff=10**8,maxsteps=100000))
print(pslq([mpf(1), U1*zeta(5), U1*pi**5],maxcoeff=10**10,maxsteps=100000))
