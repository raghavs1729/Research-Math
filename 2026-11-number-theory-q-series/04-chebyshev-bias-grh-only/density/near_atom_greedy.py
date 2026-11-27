# Start from the real-weight atom of Example 9.2 (M=1/2, F=0 on (4pi/3,2pi)), truncated at N.
# Compare its Q with the best integer sequence in {0,1} found by greedy bit-flips from roundings.
import numpy as np
L=8192; M=0.5
A0=2*np.pi/(1-np.exp(-2*np.pi*M)); A1=(2*np.pi/3)*np.exp(4*np.pi*M/3)/(1-np.exp(-2*np.pi*M/3))
b0=A1/(2*(A1-A0))
th=2*np.pi*np.arange(L)/L
def Fvals(m):
    N=len(m); c=np.zeros(L,complex); n=np.arange(1,N+1); c[1:N+1]=m*4*M/(M+1j*n)
    return 1+np.real(np.fft.ifft(c)*L)
def Qset(m,mask): f=Fvals(m)[mask]; return np.sqrt(np.mean(f**2))
def Qbest(m,eta): f=np.abs(Fvals(m)); k=int(eta*L); return np.sqrt(np.mean(np.partition(f,k)[:k]**2))
arc=(th>4*np.pi/3+0.05)&(th<2*np.pi-0.05)
for N in [30,60,120,240,480]:
    n=np.arange(1,N+1); mr=np.where(n%3==0,0.5,b0)
    best=None
    for start in [np.round(mr), np.floor(mr), np.ceil(mr)]:
        m=start.copy(); q=Qbest(m,0.25); improved=True
        while improved:
            improved=False
            for j in range(N):
                for v in (0.0,1.0):
                    if v==m[j]: continue
                    old=m[j]; m[j]=v; q2=Qbest(m,0.25)
                    if q2<q-1e-12: q=q2; improved=True
                    else: m[j]=old
        if best is None or q<best: best=q
    print(f"N={N:4d}  real Ex9.2 truncation: Q_arc={Qset(mr,arc):.4f}, Q_best(25%)={Qbest(mr,0.25):.4f} | best 0/1 integer (greedy): Q_best(25%)={best:.4f}",flush=True)
