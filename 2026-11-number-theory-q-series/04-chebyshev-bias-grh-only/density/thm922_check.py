#!/usr/bin/env python3
"""Numerical sanity checks for the draft of Theorem 9.22 (noise-shaped integer atoms).

Checks (all [num], none is part of the proof):
 (1) base comb (1,q): F_0 = 0 on the last cell J, b_0 > 1/2, for q = 20, 60 and M = 1/2, 0.02;
 (2) Lemma A (ell^1 bound): Lambda_s := ||((e^{i th}-1)^s u)^||_1 <= tau^s lambda(s);
 (3) Lemma B (shaped sup bound): V := sup_l |(e u)^(l)| <= a(N) sum_r C(k,r)(r+1)! N^-r Lambda_{k-r}
     <= (16/9) a(N) tau^k lambda(k), for band-of-origin sigma-delta errors;
 (4) Lemma C (one modified step): low coefficients exact, psi = 0 on T minus E, copy multiplicities
     vs. the bound 258 K^2 tau^k lambda(k) rho^-2 Xi(k_* rho) (+ tails);
 (5) Lemma O (offset leak): (M^2+l^2)|(O u)^(l)| vs Phi X^{-3/2} (M^2+4) U^O;
 (6) parameter table for q = 60 (theta' = 4/(1+beta)): Z, k_j, 2^{k_j}/(j+1)^{theta'}.
"""
import math
import numpy as np
from math import comb

out = []
def log(s=""):
    print(s, flush=True); out.append(s)

# ---------------------------------------------------------------- auxiliary profiles
def fexp(x):
    x = np.asarray(x, float)
    r = np.zeros_like(x); p = x > 0
    r[p] = np.exp(-1.0 / x[p]); return r
def sstep(x):                       # smooth step: 0 on (-inf,0], 1 on [1,inf)
    x = np.clip(np.asarray(x, float), -1.0, 2.0)
    a, b = fexp(x), fexp(1 - x); return a / (a + b)
xs = np.linspace(-0.5, 1.5, 400001)
ds = np.gradient(sstep(xs), xs[1] - xs[0])
S1 = np.abs(ds).max()
log(f"smooth step s: ||s'||_inf = {S1:.6f}")

# bump c(x) = exp(-1/(1-4x^2)) on (-1/2,1/2); b = c * c~ ; hat b(xi) = C(xi)^2/(2 pi), C = int c cos
yc = np.linspace(-0.5, 0.5, 200001)[1:-1]; dyc = yc[1] - yc[0]
cc = np.exp(-1 / (1 - 4 * yc ** 2))
def Cc(xi):
    return np.array([np.sum(cc * np.cos(z * yc)) * dyc for z in np.atleast_1d(xi)])
C0 = Cc(0.0)[0]
xg = np.linspace(1e-3, 120, 24000)
rat = (Cc(np.pi * xg) / C0) ** 2                  # hat b(pi x)/hat b(0)
C_Xi = max(1.0, (xg ** 3 * rat).max())
b0val = np.sum(cc * cc) * dyc                       # b(0) = int c^2 = ||b||_inf
m0b = 0.5 * (np.sum(cc) * dyc) ** 2                 # m_0(b) = (1/2) int b = (1/2)(int c)^2
log(f"bump b = c*c~: C_Xi = sup max(1, x^3 hb(pi x)/hb(0)) = {C_Xi:.4f}; b(0)/m_0(b) = {b0val/m0b:.4f}")
log(f"   hb(pi x)/hb(0) at x = 10, 20, 40: {rat[np.searchsorted(xg,10)]:.2e}, {rat[np.searchsorted(xg,20)]:.2e}, {rat[np.searchsorted(xg,40)]:.2e}")

def geometry(q):
    eta = math.pi / (4 * q); sig = math.pi / (8 * q)
    th1 = 2 * math.pi / q - eta + sig
    tau = 2 * math.sin(th1 / 2)
    return eta, sig, th1, tau

def make_u(th, q):
    eta, sig, th1, tau = geometry(q)
    d0 = np.minimum(th, 2 * np.pi - th)
    return sstep((d0 - (eta - sig)) / sig) * sstep((th1 - d0) / sig), d0

# ---------------------------------------------------------------- (1) base
log("\n(1) base comb (1,q): F_0 = 1 + 2M(b0 rho_0 + b1 rho_q) on the last cell J")
def saw(t, N, M):
    h = 2 * math.pi / N; s = np.mod(t, h)
    return h * np.exp(-M * s) / (1 - math.exp(-M * h)) - 1 / M
for q in (20, 60):
    for M in (0.5, 0.02):
        h1 = 2 * math.pi / q; x1 = 2 * math.pi * (q - 1) / q
        A0 = 2 * math.pi / (1 - math.exp(-2 * math.pi * M))
        A1 = h1 * math.exp(M * x1) / (1 - math.exp(-M * h1))
        b0 = A1 / (2 * (A1 - A0)); b1 = 0.5 - b0
        t = np.linspace(x1 + 1e-9, 2 * math.pi - 1e-9, 20001)
        F = 1 + 2 * M * (b0 * saw(t, 1, M) + b1 * saw(t, q, M))
        bnd = q / (2 * math.pi * (q - 1) * M) + 0.5
        log(f"  q={q:3d} M={M:5.2f}: b0={b0:.6f} (bound 1/((q-1)u)+1/2 = {bnd:.6f}), max|F_0| on J = {np.abs(F).max():.1e}")

# ---------------------------------------------------------------- (2) Lemma A
log("\n(2) Lemma A: Lambda_s <= tau^s lambda(s), lambda(s) = 1 + (pi/sqrt3)(s/tau + ||u'||)")
for q in (20, 60):
    L = 2 ** 20
    th = 2 * np.pi * np.arange(L) / L
    u, d0 = make_u(th, q)
    eta, sig, th1, tau = geometry(q)
    up = S1 / sig
    lam = lambda s: 1 + (math.pi / math.sqrt(3)) * (s / tau + up)
    z = np.exp(1j * th) - 1
    worst = 0
    for s in range(0, 13):
        Ls = np.abs(np.fft.fft(z ** s * u) / L).sum()
        worst = max(worst, Ls / (tau ** s * lam(s)))
        if s in (0, 2, 4, 8, 12):
            log(f"  q={q}: s={s:2d}  Lambda_s={Ls:.3e}  tau^s lambda(s)={tau**s*lam(s):.3e}  ratio={Ls/(tau**s*lam(s)):.3f}")
    log(f"  q={q}: tau={tau:.5f}, ||u'||={up:.2f}, max ratio over s<=12: {worst:.3f} (must be <= 1)")

# ---------------------------------------------------------------- (3),(4) one modified step
log("\n(3)-(4) one band [N,4N), band-of-origin sigma-delta of order k, modified correction")
M = 0.5; q = 20; K = 4
eta, sig, th1, tau = geometry(q)
up = S1 / sig
lam = lambda s: 1 + (math.pi / math.sqrt(3)) * (s / tau + up)
N = 256; N1 = K * N; W = 2 * N1; Np = 4 * W; rho = 0.1
L = 2 ** 21
th = 2 * np.pi * np.arange(L) / L
u, d0 = make_u(th, q)
T = (d0 > eta) & (d0 < 2 * math.pi / q - eta)
z = np.exp(1j * th) - 1
Lam = [np.abs(np.fft.fft(z ** s * u) / L).sum() for s in range(0, 20)]
a = lambda n: 2 * M * M / (M * M + np.asarray(n, float) ** 2)
freq = np.fft.fftfreq(L, d=1.0 / L).astype(np.int64)
rng = np.random.default_rng(7)
band = np.arange(N, N1)
targ = rng.uniform(0, 1, len(band)) + 3.0
# comb g = beta/hat beta(0), period 2 pi/Np, gap fraction rho
w = 2 * math.pi * rho / Np
yb = np.linspace(-1, 1, 40001); dyb = yb[1] - yb[0]
cb = np.zeros_like(yb); inn = np.abs(yb) < 0.5; cb[inn] = np.exp(-1 / (1 - 4 * yb[inn] ** 2))
bt = np.convolve(cb, cb[::-1], mode="same") * dyb
hper = 2 * math.pi / Np
sper = np.mod(th + hper / 2, hper) - hper / 2
beta = np.interp(2 * sper / w, yb, bt, left=0, right=0)
g = beta / beta.mean()
E = g > 0
ghat = np.fft.fft(g) / L
for k in (0, 2, 4, 6, 8):
    h = np.array([(-1) ** i * comb(k, i) for i in range(k + 1)], float)
    e = np.zeros(len(band)); m = np.zeros(len(band))
    for i in range(len(band)):
        fb = sum(h[r] * e[i - r] for r in range(1, k + 1) if i - r >= 0)
        uu = targ[i] + fb
        m[i] = math.floor(uu + 0.5); e[i] = m[i] - uu
    d = np.convolve(h, e)                       # band-of-origin piece, support [N, N1+k)
    nn = np.arange(N, N + len(d))
    assert np.allclose(d[:len(band)], m - targ)
    ehat = np.zeros(L, complex); ehat[nn] = a(nn) * d; ehat[-nn] = a(nn) * d
    eP = np.real(np.fft.ifft(ehat) * L)
    v = eP * u
    vhat = np.fft.fft(v) / L
    V = np.abs(vhat).max()
    Bfull = a(N) * sum(comb(k, r) * math.factorial(r + 1) * N ** (-r) * Lam[k - r] for r in range(k + 1))
    Bsimp = (16 / 9) * a(N) * tau ** k * lam(k)
    # leak and its removal with the same comb
    low = np.abs(freq) < W
    Lk = np.fft.fft(v * g) / L - vhat           # sum_{k != 0} vhat(n - k Np) ghat(k Np)
    Phat = np.where(low, Lk, 0)
    P = np.real(np.fft.ifft(Phat) * L)
    psi = eP * (1 - u) + (v - P) * g
    psihat = np.fft.fft(psi) / L
    exact = np.abs(psihat[low] - ehat[low]).max() / a(N)
    vanish = np.abs(psi[T & ~E]).max()
    hi = np.abs(freq) >= W
    mult = np.abs(psihat[hi]) * (M * M + freq[hi].astype(float) ** 2) / (2 * M * M)
    # pointwise main bound (16/9) tau^k lam(k) (M^2+n^2)/(M^2+N^2) |ghat(k_* Np)|  (Lemma C)
    fh = freq[hi]; ks = np.rint(fh / Np).astype(np.int64)
    gks = np.abs(ghat[np.mod(ks * Np, L)]) * (ks != 0)
    mainb = V * gks * (M * M + fh.astype(float) ** 2) / (2 * M * M)
    excess = (mult - mainb).max()
    crude = 258 * K * K * tau ** k * lam(k) * rho ** -2 * C_Xi
    log(f"  k={k}: V/a(N)={V/a(N):.2e}  Leibniz bd={Bfull/a(N):.2e}  simple bd={Bsimp/a(N):.2e}  "
        f"| low coeff err/a(N)={exact:.1e}  max|psi| on T\\E={vanish:.1e}  "
        f"| max copy mult={mult.max():.2e}  max(copy - sup-bound)={excess:.1e}  crude bd={crude:.2e}")

# ---------------------------------------------------------------- (5) Lemma O
log("\n(5) Lemma O: offset leak, f(x) = a(x) c(x), c(x) = A (log x)^th' chi(x/X)")
M = 0.5; q = 20; thp = 0.92; A = 1.0
L = 2 ** 20
th = 2 * np.pi * np.arange(L) / L
u, d0 = make_u(th, q)
uhat = np.fft.fft(u) / L
p = np.fft.fftfreq(L, d=1.0 / L)
Usharp = np.sum((1 + p ** 2) * (2 + np.abs(p)) * np.abs(uhat))
UO = 2 ** 3.5 * Usharp
x2 = np.linspace(-1, 3, 400001); chi2 = sstep(2 * x2 - 1)
d1 = np.gradient(chi2, x2[1] - x2[0]); d2 = np.gradient(d1, x2[1] - x2[0])
Xchi = 1 + np.abs(d1).max() + np.abs(d2).max()
log(f"  U^sharp = {Usharp:.3e}, U^O = {UO:.3e}, X_chi = 1+||chi'||+||chi''|| = {Xchi:.3f}")
for X in (2 ** 10, 2 ** 12, 2 ** 14):
    n = np.arange(1, L // 2)
    cn = A * np.log(n) ** thp * sstep(2 * n / X - 1)
    Ohat = np.zeros(L); Ohat[n] = a(n) * cn; Ohat[-n] = a(n) * cn
    O = np.real(np.fft.ifft(Ohat) * L)
    Ouhat = np.fft.fft(O * u) / L
    lhs = (np.abs(Ouhat) * (M * M + p ** 2))[np.abs(p) < L // 4].max()
    Phi = 44 * M * M * A * Xchi
    rhs = Phi * X ** -1.5 * (M * M + 4) * UO
    log(f"  X={X:6d}: max_l (M^2+l^2)|(Ou)^(l)| = {lhs:.3e}   Lemma O bound = {rhs:.3e}   (mean (Ou)^(0) = {Ouhat[0].real:.2e})")

# ---------------------------------------------------------------- (6) parameters for q = 60
log("\n(6) parameter table, q = 60, K = 4, delta = 1/2, rho_j = rho_0 (j+1)^{-3/2}")
q = 60; K = 4
eta, sig, th1, tau = geometry(q)
up = S1 / sig
lam = lambda s: 1 + (math.pi / math.sqrt(3)) * (s / tau + up)
beta = math.log2(1 / tau); thp = 4 / (1 + beta)
zeta32 = 2.6123753486854883
rho0 = 3 / (8 * q * zeta32)
SK = 4 / (1 - K ** -2) + 1.5 / (1 - 1 / K)
Z = 8 * 258 * C_Xi * SK
log(f"  tau={tau:.5f} beta={beta:.4f} theta'=4/(1+beta)={thp:.4f}  rho0={rho0:.5f}  S_K={SK:.4f}  Z={Z:.4e}")
def kj(j):
    Y = Z * K * K * (rho0 * (j + 1) ** -1.5) ** -2
    k = 1
    while Y * tau ** k * lam(k) > 2 ** k:
        k += 1
    return k
for j in (0, 1, 5, 10, 30, 100, 1000, 10 ** 4, 10 ** 6):
    k = kj(j)
    log(f"  j={j:8d}: k_j={k:3d}  2^k_j={2**k:.3e}  2^k_j/(j+1)^theta'={2**k/(j+1)**thp:.3e}  4k_j/tau={4*k/tau:.0f}")
Amax = max((2 ** (kj(j) + 1) + 1) / (j * math.log(4)) ** thp for j in list(range(1, 3000)) + [10 ** 4, 10 ** 5, 10 ** 6])
log(f"  A := sup_J (2^(k_J+1)+1)/(J log 4)^theta' ~ {Amax:.3e} (sampled)")

# ---------------------------------------------------------------- (7) the threshold j_0 (q = 60, M = 1/2)
log("\n(7) constants of (J1)-(J4) and X_low for q = 60, M = 1/2")
M = 0.5
L = 2 ** 22
th = 2 * np.pi * np.arange(L) / L
u, d0 = make_u(th, q)
uhat = np.fft.fft(u) / L
p = np.fft.fftfreq(L, d=1.0 / L)
Usharp = np.sum((1 + p ** 2) * (2 + np.abs(p)) * np.abs(uhat)); UO = 2 ** 3.5 * Usharp
U4 = np.mean(np.abs(np.real(np.fft.ifft(uhat * p ** 4) * L)))       # (1/2pi) int |u''''|
w0 = np.interp(4 * d0 / eta, yb, bt, left=0, right=0)                # b(4|theta|/eta), b even
w0 = w0 / w0.mean()                                                  # hat w_0(0) = 1
w0hat = np.abs(np.fft.fft(w0) / L)
W0 = np.max((M * M + p ** 2) * w0hat)
log(f"  U^sharp={Usharp:.3e}  U^O={UO:.3e}  U_4={U4:.3e}  W_0={W0:.3e}  ||w_0||_inf={w0.max():.3e}  X_chi={Xchi:.2f}")
A_ = Amax
need = 16 * 22 * Xchi * A_ * (M * M + 4) * UO * (1 + W0 / (M * M))
X_J3 = need ** (2 / 3)
j_J3 = math.ceil(math.log(X_J3, 4))
j_J1 = next(j for j in range(2, 200) if all(4 ** i >= 4 * kj(i) / tau and kj(i) <= 4 ** i for i in range(j, j + 200)))
def eps(j):
    return 2 ** kj(j) * U4 * ((1 + M * M) + 3 * C_Xi * (rho0 * (j + 1) ** -1.5) ** -2) * 4.0 ** (-3 * j)
j_J4 = next(j for j in range(2, 200) if sum(eps(i) for i in range(j, j + 120)) <= 1 / 16)
j0 = max(j_J1, j_J3, j_J4, 2)
Xlow = 82 * K * K * C_Xi * SK * (rho0 * j0 ** -1.5) ** -2
log(f"  (J1) holds from j = {j_J1};  (J3) needs X >= {X_J3:.3e}, i.e. j >= {j_J3};  (J4) holds from j = {j_J4}")
log(f"  => j_0 = {j0} (N_j0 = 4^{j0} = {4.0**j0:.3e});  X_low = 82 K^2 C_Xi S_K rho_(j0-1)^-2 = {Xlow:.3e}")
ser = sum((j + 1) ** 1.5 / (4 ** j - 1) for j in range(1, 200))
log(f"  sum_(j>=1) (j+1)^(3/2)/(4^j-1) = {ser:.6f};  1.3612 + it = {1.3612+ser:.4f};  8*(that)*8 zeta(3/2)/3 = {8*(1.3612+ser)*8*zeta32/3:.2f}")

open(__file__.replace(".py", "_results.txt"), "w").write("\n".join(out) + "\n")
