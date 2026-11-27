import mpmath as mp
mp.mp.dps=40
N=700000
isp=bytearray([1])*(N+1); isp[0]=isp[1]=0
for i in range(2,int(N**.5)+1):
    if isp[i]: isp[i*i::i]=bytearray(len(isp[i*i::i]))
pr=[p for p in range(3,N+1) if isp[p]]
def C(alpha,x):
    x=mp.mpf(x); s=mp.mpf(0); a=mp.mpf(alpha)
    for p in pr:
        if p>1.2*x: break
        s+= (1 if p%4==1 else -1)*mp.exp(-(p/x)**a)
    return s
print("C_2048(623658)=",C(2048,623658))
print("C_4096(26862.745)=",C(4096,'26862.745'))
# D(x)=pi(x;4,1)-pi(x;4,3)
d=0
for p in pr:
    d+= 1 if p%4==1 else -1
    if d>0: print("first lead of 1 mod 4 at",p); break
