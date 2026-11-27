#!/usr/bin/env python3
"""Referee toy run of the FULL shaped cascade (Sec. 4.2-4.3 of the draft) over 3 consecutive shaped bands,
with base comb (1,q), offset O(1-u)+mu_O w_0, band-of-origin sigma-delta with spill, leak removal P_j.
Checks: integrality of the final multiplicities on all rounded bands (incl. spill zones), m_n >= 1,
the band-of-origin identity, vanishing of (1-D/M)(Psi - Psi_0) on T minus the combs, and the copy bound (3.3)."""
import math, numpy as np
from math import comb

def fexp(x):
    x = np.asarray(x, float); r = np.zeros_like(x); p = x > 0; r[p] = np.exp(-1.0/x[p]); return r
def sstep(x):
    x = np.clip(np.asarray(x, float), -1.0, 2.0); a, b = fexp(x), fexp(1-x); return a/(a+b)

q, M, K = 20, 0.5, 4
eta = math.pi/(4*q); sig = math.pi/(8*q); th1 = 2*math.pi/q - eta + sig; tau = 2*math.sin(th1/2)
L = 2**23
th = 2*np.pi*np.arange(L)/L; d0 = np.minimum(th, 2*np.pi - th)
u = sstep((d0-(eta-sig))/sig)*sstep((th1-d0)/sig)
T = (d0 > eta) & (d0 < 2*math.pi/q - eta)
Jp = T & (th > math.pi)                                   # J' (the arc near 2pi)
freq = np.fft.rfftfreq(L, 1.0/L).astype(np.int64)          # 0..L/2
a = lambda n: 2*M*M/(M*M + np.asarray(n, float)**2)
fwd = lambda f: np.fft.rfft(f)/L                           # coefficients hat f(n), n >= 0 (real even f)
def inv(c):                                                # even real function from hat f(n), n>=0
    return np.fft.irfft(c*L, n=L)

# bump b = c*c~ sampled
yb = np.linspace(-1, 1, 40001); dyb = yb[1]-yb[0]
cb = np.zeros_like(yb); inn = np.abs(yb) < 0.5; cb[inn] = np.exp(-1/(1-4*yb[inn]**2))
bt = np.convolve(cb, cb[::-1], mode="same")*dyb
def comb_g(Np, rho):
    w = 2*math.pi*rho/Np; hper = 2*math.pi/Np
    s = np.mod(th + hper/2, hper) - hper/2
    beta = np.interp(2*s/w, yb, bt, left=0, right=0)
    return beta/beta.mean(), beta > 0

# ---- base comb (1,q): coefficients a(n) m_base(n)
v_q = 2*math.pi*M/q
A0 = 2*math.pi/(1-math.exp(-2*math.pi*M)); A1 = (2*math.pi/q)*math.exp(2*math.pi*M*(q-1)/q)/(1-math.exp(-v_q))
b0 = A1/(2*(A1-A0))
mbase = np.where(freq % q == 0, 0.5, b0).astype(float)
PsiC = a(freq)*mbase; PsiC[0] = 1.0                       # coefficient array of Psi (n >= 0)
Psi0C = PsiC.copy()

# ---- offset O(1-u) + mu_O w_0 with c(n) = Coff * chi(n/X)
N0 = 64; X = N0; kk = 4; rho = 0.1
Coff = 2**(kk+1) + 1
chi = sstep(2*freq/X - 1)
cO = Coff*chi; cO[0] = 0
O = inv(a(freq)*cO)
OuC = fwd(O*u); muO = OuC[0].real
w0 = np.interp(4*d0/eta, yb, bt, left=0, right=0); w0 /= w0.mean()
offC = fwd(O*(1-u) + muO*w0)
print(f"offset: mean of O(1-u)+mu w0 = {offC[0].real:.2e}; max leak |eps^O_n| (n<=4096) = "
      f"{np.abs((offC[1:4097].real - a(freq[1:4097])*cO[1:4097])/a(freq[1:4097])).max():.3e}")
PsiC = PsiC + offC.real

# ---- three shaped bands
h = np.array([(-1)**i*comb(kk, i) for i in range(kk+1)], float)
rounded = np.zeros(L//2+1, bool); Emask = np.zeros(L, bool)
store = []
for j in range(3):
    Nj = N0*K**j; N1 = K*Nj; W = 2*N1; Np = 4*W
    g, Ej = comb_g(Np, rho); Emask |= Ej
    targ = PsiC[Nj:N1].real/a(np.arange(Nj, N1))
    e = np.zeros(N1-Nj); m = np.zeros(N1-Nj)
    for i in range(N1-Nj):
        y = targ[i] + sum(h[r]*e[i-r] for r in range(1, kk+1) if i-r >= 0)
        m[i] = math.floor(y+0.5); e[i] = m[i]-y
    d = np.convolve(h, e)                                  # supported on [Nj, N1+kk)
    nn = np.arange(Nj, Nj+len(d))
    eC = np.zeros(L//2+1); eC[nn] = a(nn)*d
    eP = inv(eC); v = eP*u
    vC = fwd(v); vgC = fwd(v*g)
    Lk = np.where(freq < W, vgC - vC, 0)                   # leak L_j(n), |n| < W
    P = inv(Lk)
    psi = eP*(1-u) + (v-P)*g
    psiC = fwd(psi)
    lowerr = np.abs(psiC[:W] - eC[:W]).max()
    vanish = np.abs(psi[T & ~Ej]).max()
    # copy bound (3.3) check on n >= W, up to the Nyquist
    V = np.abs(vC).max()
    lam_k = 1 + (math.pi/math.sqrt(3))*(kk/tau + 2.0/sig)
    hi = freq >= W
    mult = np.abs(psiC[hi].real)/a(freq[hi])
    ks = np.rint(freq[hi]/Np)
    ghat = np.abs(fwd(g).real)
    gk = ghat[np.minimum((ks*Np).astype(np.int64), L//2)]*(ks >= 1)
    main_bd = 258*K*K*tau**kk*lam_k*(ks**2*gk)              # 258 K^2 tau^k lam rho^-2 Xi(k* rho)
    print(f"band {j}: [{Nj},{N1}) W={W} N'={Np}: low-coef err/a(N)={lowerr/a(Nj):.1e}, max|psi| on T\\E={vanish:.1e}, "
          f"sup|v^|/a(N)={V/a(Nj):.2e}, max copy mult={mult.max():.3e}, max(copy/(3.3) main bound)={np.max(mult/np.maximum(main_bd,1e-300)*(main_bd>0)):.2e}, "
          f"max copy where k*=0: {mult[ks==0].max() if np.any(ks==0) else 0:.2e}")
    store.append((Nj, N1, e, d, targ, m))
    PsiC = PsiC + psiC.real
    rounded[Nj:N1] = True

# ---- final checks
n_r = np.nonzero(rounded)[0]
mfin = PsiC[n_r]/a(n_r)
print(f"final: max dist to Z on rounded bands [{n_r[0]},{n_r[-1]}] = {np.abs(mfin-np.rint(mfin)).max():.2e}; "
      f"min m_n = {mfin.min():.3f}; max m_n = {mfin.max():.3f}")
# band-of-origin identity: m_n - (real contributions) = sum_j (h * e1_Bj)_n ; check the spill zones
for (Nj, N1, e, d, targ, m) in store:
    print(f"   band [{Nj},{N1}): m - target = d on band: {np.abs((m - targ) - d[:N1-Nj]).max():.1e}; "
          f"spill d[N1..N1+k) = {np.round(d[N1-Nj:],3)}")
# vanishing of (1 - D/M)(Psi - Psi_0) on J' minus the combs
dC = PsiC - Psi0C
Dpart = np.fft.irfft((1j*freq*dC)*L, n=L)
Fc = inv(dC) - Dpart/M
S = Jp & ~Emask
print(f"(1-D/M)(Psi - Psi_0): max on S = J' minus combs: {np.abs(Fc[S]).max():.2e} (vs max overall {np.abs(Fc).max():.2e}); "
      f"measure fraction of S in J' = {S.sum()/Jp.sum():.3f}")
