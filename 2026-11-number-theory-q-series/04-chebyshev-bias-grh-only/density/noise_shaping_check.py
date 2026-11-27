"""Noise shaping of the integer rounding error on a target arc (heuristic for Remark 9.21).

One frequency band n in [N, 2N), random real targets t_n in [0.3, 0.8], cyclic model at M = 1/2.
k-th order sigma-delta: m_n = round(t_n + c0 - v_n), v_n = sum_{j>=1} h_j e_{n-j}, e_n = m_n - (t_n + c0 - v_n),
with H(z) = (1-z)^k (h_j = (-1)^j C(k,j)), |e_n| <= 1/2, 0 <= m_n <= 2^k + 2 (offset c0 = B/2).
The shaped part of the error, delta = h * (w e) (w a smooth taper of the band, to remove band-edge
artefacts that do not occur when the scheme runs continuously over all n), has
sum_n delta_n e^{in theta} = H(theta) * sum_n w_n e_n e^{in theta}, small where |H| is small.
Target set T = {theta : eta < dist(theta, 0) < 2 pi/q - eta}, q = 20, i.e. J' u -J' for the last
cell J of the comb (1, q).  |H| <= (2 sin(pi/q))^k on T.
Reported: rms of the F-error on T, and max over all frequencies of the coefficients of e^Psi * u in
multiplicity units, u smooth, = 1 on T (the quantity that the modified correction of Remark 9.21 pays).
"""
import numpy as np
from math import comb
rng=np.random.default_rng(3)
M=0.5; q=20; N=4096; L=1<<16
n=np.arange(N,2*N)
t=rng.uniform(0.3,0.8,len(n))
th=2*np.pi*np.arange(L)/L
eta=0.03; sig=0.01
d0=np.minimum(th,2*np.pi-th)
inT=(d0>eta)&(d0<2*np.pi/q-eta)
def smoothstep(x):
    x=np.clip(x,1e-12,1-1e-12); return 1/(1+np.exp(np.clip(1/x-1/(1-x),-700,700)))
u=smoothstep((d0-(eta-sig))/sig)*smoothstep(((2*np.pi/q-eta+sig)-d0)/sig)
u[(d0<eta-sig)|(d0>2*np.pi/q-eta+sig)]=0
def analyse(nn,dm,label):
    c=np.zeros(L,complex); c[nn]=2*M*dm/(M+1j*nn); c[-nn]=np.conj(c[nn])
    E=np.real(np.fft.ifft(c)*L)
    cP=np.zeros(L,complex); cP[nn]=2*M*M*dm/(M*M+nn*nn); cP[-nn]=cP[nn]
    ePsi=np.real(np.fft.ifft(cP)*L)
    pu=np.fft.fft(ePsi*u)/L; k=np.arange(L); kk=np.minimum(k,L-k)
    mult=np.abs(pu)*(M*M+kk**2)/(2*M*M)
    print(f"{label:40s} rms E on T = {np.sqrt(np.mean(E[inT]**2)):.3e}   max|mult(e^Psi u)| = {mult.max():.3e}",flush=True)
w=np.sin(np.pi*(n-N)/N)**2
analyse(n,w*(np.ceil(t)-t),"tapered round-up (Thm 9.6 rounding)")
for k in range(0,9):
    h=np.array([(-1)**i*comb(k,i) for i in range(k+1)],float)
    B=2**k+2; c0=B/2.0
    e=np.zeros(len(n)); m=np.zeros(len(n))
    for i in range(len(n)):
        v=sum(h[j]*e[i-j] for j in range(1,k+1) if i-j>=0)
        uu=t[i]+c0-v
        m[i]=np.round(uu); e[i]=m[i]-uu
    assert m.min()>=0 and m.max()<=B and np.abs(e).max()<=0.5+1e-12
    d=np.convolve(h,w*e); nn=np.arange(N,N+len(d))
    analyse(nn,d,f"k={k}: B={B}, sup_T|H|={(2*np.sin(np.pi/q))**k:.1e}")
