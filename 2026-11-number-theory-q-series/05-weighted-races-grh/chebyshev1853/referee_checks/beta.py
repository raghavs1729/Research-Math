import mpmath as mp, sys
mp.mp.dps=25
q=int(sys.argv[1])
SG,c1,c3,p0=(3.9199328e-4,2.071193,2.153716634,3) if q==4 else (1.6531314e-5,2.081942,2.074662259,2)
SG,c1,c3=map(mp.mpf,(SG,c1,c3))
U1=mp.mpf(2)**22; ea=0.2/mp.log(U1)**2; ct=1+ea
hs=mp.sqrt(mp.pi)/2
F=lambda v: v*mp.exp(-v*v)+hs*mp.erfc(v)
YE,YI=mp.mpf(10)**12,mp.mpf(10)**16
def kflat(y,reg):
    Eb= ea*F(U1/mp.sqrt(YE)) if reg==1 else ea*hs
    v=hs-SG-Eb-(c1-1+ea)*mp.gamma(1.25)*y**-0.25-c1*mp.gamma(1.125)*y**-0.375
    t=0;k=3
    while p0**k<y:
        t+= (c3*mp.gamma(1+1/(2*mp.mpf(k)))*y**(1/(2*mp.mpf(k))) if reg<3 else ct*mp.gamma(1+1/mp.mpf(k))*y**(1/mp.mpf(k)))
        k+=2
    t+= (c3*mp.sqrt(p0)/(1-1/mp.mpf(p0))+1+1e-100) if reg<3 else ct*p0/(1-1/mp.mpf(p0))+1
    return v-t/mp.sqrt(y)
KE=kflat(YE,2); KI=kflat(YI,3)
print("Ebar_small",ea*F(U1/mp.sqrt(YE)),"K_E",KE,"K_I",KI)
def ks(y):
    r=kflat(y,1 if y<=YE else 2 if y<=YI else 3)
    if y<=YE: r=min(r,KE)
    if y<=YI: r=min(r,KI)
    return max(r,0)
def A(V): s=mp.sqrt(1+V); return -mp.log((s+1)/(s-1))
ca=2*mp.log((1+mp.sqrt(2))/2)
rho=mp.mpf(10)**0.01
def beta(a):
    v=ks(a/2)*(1-ca/mp.log(a)); Vp=1; Vn=1; i=1
    while True:
        Vc=2**i
        if a/(1+Vc)<500 or Vc>4096: break
        v+=ks(a/(1+Vc))*(A(Vc)-A(Vp))*mp.log(a)/mp.log(rho*a/Vp)**2
        Vp=Vc; Vn=Vc; i+=1
    An=mp.log(a/Vn)
    return v-mp.sqrt(mp.pi)*mp.sqrt(1+Vn)*(1+An)*mp.log(rho*a)/(An**2*a)
m=min((beta(mp.mpf(10)**(mp.mpf(j)/100)),j) for j in range(660,1401))
print("min beta on [10^6.6,10^14]:",m, " beta(X_A)=",beta(mp.mpf(10)**6.6))
print("A(64)-A(1)-ca",A(64)-A(1)-ca)
