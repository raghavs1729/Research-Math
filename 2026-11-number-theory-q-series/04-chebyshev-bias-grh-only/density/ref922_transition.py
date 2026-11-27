#!/usr/bin/env python3
"""Referee toy run 2: low (round-up, Lemma D) bands j=0,1, offset ramp inside band 1, shaped bands j=2,3
(with spill 2->3). Checks integrality on [64, 16384), m_n >= 1, the (5.2)-type lower bound on targets,
physical vanishing of Psi - Psi_0 on S, and per-piece size of the spectral (1-D/M) residual (to attribute
any residual to discretisation rather than to the construction)."""
import math, numpy as np
from math import comb

def fexp(x):
    x = np.asarray(x, float); r = np.zeros_like(x); p = x > 0; r[p] = np.exp(-1.0/x[p]); return r
def sstep(x):
    x = np.clip(np.asarray(x, float), -1.0, 2.0); a, b = fexp(x), fexp(1-x); return a/(a+b)

q, M, K = 20, 0.5, 4
eta = math.pi/(4*q); sig = math.pi/(8*q); th1 = 2*math.pi/q - eta + sig; tau = 2*math.sin(th1/2)
L = 2**24
th = 2*np.pi*np.arange(L)/L; d0 = np.minimum(th, 2*np.pi - th)
u = sstep((d0-(eta-sig))/sig)*sstep((th1-d0)/sig)
T = (d0 > eta) & (d0 < 2*math.pi/q - eta)
Jp = T & (th > math.pi)
freq = np.fft.rfftfreq(L, 1.0/L).astype(np.int64)
a = lambda n: 2*M*M/(M*M + np.asarray(n, float)**2)
fwd = lambda f: np.fft.rfft(f).real/L
inv = lambda c: np.fft.irfft(c*L, n=L)
Dres = lambda c, mask: np.abs(inv(c)[mask] - np.fft.irfft(1j*freq*c*L, n=L)[mask]/M).max()

yb = np.linspace(-1, 1, 40001); dyb = yb[1]-yb[0]
cb = np.zeros_like(yb); inn = np.abs(yb) < 0.5; cb[inn] = np.exp(-1/(1-4*yb[inn]**2))
bt = np.convolve(cb, cb[::-1], mode="same")*dyb
def comb_g(Np, rho):
    w = 2*math.pi*rho/Np; hper = 2*math.pi/Np
    s = np.mod(th + hper/2, hper) - hper/2
    beta = np.interp(2*s/w, yb, bt, left=0, right=0)
    return beta/beta.mean(), beta > 0

v_q = 2*math.pi*M/q
A0 = 2*math.pi/(1-math.exp(-2*math.pi*M)); A1 = (2*math.pi/q)*math.exp(2*math.pi*M*(q-1)/q)/(1-math.exp(-v_q))
b0 = A1/(2*(A1-A0))
PsiC = a(freq)*np.where(freq % q == 0, 0.5, b0); PsiC[0] = 1.0
Psi0C = PsiC.copy()

N0 = 64; kk = 4; X = 1024
Coff = 2**(kk+1) + 1
cO = Coff*sstep(2*freq/X - 1); cO[0] = 0
O = inv(a(freq)*cO); muO = fwd(O*u)[0]
w0 = np.interp(4*d0/eta, yb, bt, left=0, right=0); w0 /= w0.mean()
offC = fwd(O*(1-u) + muO*w0)
lk = np.abs((offC[1:16384] - a(freq[1:16384])*cO[1:16384])/a(freq[1:16384]))
print(f"offset: mean={offC[0]:.1e}; max leak eps^O_n on n<16384: {lk.max():.3e} (n<1024: {lk[:1023].max():.3e}, n>=1024: {lk[1023:].max():.3e})")
PsiC = PsiC + offC
Emask = np.zeros(L, bool); pieces = [("offset", offC.copy())]
h = np.array([(-1)**i*comb(kk, i) for i in range(kk+1)], float)
spill_in = None
for j in range(4):
    Nj = N0*K**j; N1 = K*Nj; W = 2*N1; Np = 4*W; rho = 0.1
    g, Ej = comb_g(Np, rho); Emask |= Ej
    nB = np.arange(Nj, N1)
    targ = PsiC[Nj:N1]/a(nB)
    if j < 2:                                       # round-up step (Lemma D)
        dm = np.ceil(targ) - targ
        eC = np.zeros(L//2+1); eC[nB] = a(nB)*dm
        psi = inv(eC)*g
        kind = "round-up"
    else:                                           # shaped step (Lemma C)
        e = np.zeros(N1-Nj); m = np.zeros(N1-Nj)
        for i in range(N1-Nj):
            y = targ[i] + sum(h[r]*e[i-r] for r in range(1, kk+1) if i-r >= 0)
            m[i] = math.floor(y+0.5); e[i] = m[i]-y
        d = np.convolve(h, e); nn = np.arange(Nj, Nj+len(d))
        eC = np.zeros(L//2+1); eC[nn] = a(nn)*d
        eP = inv(eC); v = eP*u
        vC = fwd(v); vgC = fwd(v*g)
        P = inv(np.where(freq < W, vgC - vC, 0))
        psi = eP*(1-u) + (v-P)*g
        kind = f"shaped, spill={np.round(d[N1-Nj:],2)}"
    psiC = fwd(psi)
    ok_low = np.abs(psiC[:min(W, 7*N1) if j < 2 else W] - eC[:min(W, 7*N1) if j < 2 else W]).max()/a(Nj)
    print(f"band {j} [{Nj},{N1}) {kind}: min target {targ.min():.3f}, low-coef err {ok_low:.1e}, "
          f"max|psi| on (T or circle)\\E = {np.abs(psi[(T if j>=2 else np.ones(L,bool)) & ~Ej]).max():.1e}, "
          f"max copy mult {np.abs(psiC[W:])[:].max()/1:.2e}->units {np.max(np.abs(psiC[W:])/a(freq[W:])):.3e}")
    pieces.append((f"psi_{j}", psiC.copy()))
    PsiC = PsiC + psiC

nr = np.arange(N0, N0*K**4)
mf = PsiC[nr]/a(nr)
print(f"final on [{nr[0]},{nr[-1]}]: max dist to Z = {np.abs(mf - np.rint(mf)).max():.2e}, min m = {mf.min():.3f}, max m = {mf.max():.3f}")
S = Jp & ~Emask
dC = PsiC - Psi0C
print(f"physical max |Psi - Psi_0| on S: {np.abs(inv(dC)[S]).max():.2e};  fraction of J' kept: {S.sum()/Jp.sum():.3f}")
for name, c in pieces:
    print(f"   spectral (1-D/M) residual of {name} on S: {Dres(c, S):.2e}   (sup of piece {np.abs(inv(c)).max():.2e})")
