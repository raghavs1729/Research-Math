"""Numerical check of Proposition 9.17 (truncated to finitely many levels), M = 1/2.

Test A: 3 levels, ratio 1/10. delta m_n from the ANALYTIC formula
        delta m_n = eps sum_k (2/P_k) Re[ hatpsi_{l_k}(n l_k) sum_{x in X_k} e^{-inx} ],
        hatpsi_l(xi) = -i xi (xi^2 + M^2 l^2) hattheta0(xi)  (hattheta0 by quadrature + spline),
        then phi = synthesis of 2M dm_n/(M+in) by FFT, compared with the explicit bump solutions,
        and F = F_0 + phi on S;  also full F from all coefficients m0 + dm.
Test B: 3 levels with ratio ~1/2000 (other-level leakage negligible); m_{n_k+r} - m0_{n_k+r}
        from the analytic formula vs the predicted 2 eps q_k Re[hatpsi_0(1) hatnu_k(r)].
"""
import numpy as np
from scipy.interpolate import CubicSpline

M = 0.5
J0, J1 = 4*np.pi/3, 2*np.pi
LJ = J1 - J0

# ---------- bump theta0 and its Fourier transform ----------
uq, wq = np.polynomial.legendre.leggauss(6000)
uq = 0.5*uq; wq = 0.5*wq
def th0(u):
    u = np.asarray(u, float); out = np.zeros_like(u); m = np.abs(u) < 0.5
    out[m] = np.exp(-1/(1-4*u[m]**2)); return out
def _d1(u):
    q = 1-4*u**2; return (-8*u/q**2)*np.exp(-1/q)
def _d2(u):
    q = 1-4*u**2; p1 = -8*u/q**2; p2 = -8/q**2 - 128*u**2/q**3
    return (p2 + p1**2)*np.exp(-1/q)
_DER = {1: _d1, 2: _d2}
def th0_d(u, k):
    fk = _DER[k]
    u = np.asarray(u, float); out = np.zeros_like(u); m = np.abs(u) < 0.5
    out[m] = fk(u[m]); return out
thq = th0(uq)
def that_quad(xi):
    xi = np.atleast_1d(np.asarray(xi, float)); out = np.empty_like(xi)
    for i in range(0, len(xi), 4000):
        out[i:i+4000] = (np.cos(np.outer(xi[i:i+4000], uq)) @ (wq*thq))/(2*np.pi)
    return out
XI_MAX = 3000.0
xig = np.arange(0, XI_MAX+0.01, 0.01)
spl = CubicSpline(xig, that_quad(xig))
def that(xi):
    a = np.abs(xi); out = np.zeros_like(a, dtype=float); m = a <= XI_MAX
    out[m] = spl(a[m]); return out           # even function; beyond XI_MAX negligible (<1e-20)
def hpsi(xi, l):
    return -1j*xi*(xi**2 + (M*l)**2)*that(xi)
# C0 := sup_{l<=1} |hpsi_l(xi)| / min(|xi|, |xi|^-2)
xs = np.linspace(1e-6, XI_MAX, 600001)
T = that(xs)
C0 = max(np.max(((xs**2+M**2)*np.abs(T))[xs <= 1]), np.max((xs**3*(xs**2+M**2)*np.abs(T))[xs > 1]))
eps = 1/(32*C0)
hth1 = that_quad([1.0])[0]
print(f"hattheta0(1) = {hth1:.10f},  C0 = {C0:.6g},  eps = {eps:.6g}")

# ---------- Example 9.2 base ----------
A0 = 2*np.pi/(1-np.exp(-2*np.pi*M))
A1 = (2*np.pi/3)*np.exp(4*np.pi*M/3)/(1-np.exp(-2*np.pi*M/3))
b0 = A1/(2*(A1-A0)); b1 = 0.5 - b0
def rho(t, N):
    h = 2*np.pi/N; s = np.mod(t, h)
    return h*np.exp(-M*s)/(1-np.exp(-M*h)) - 1/M
def F0(t): return 1 + 2*M*(b0*rho(t, 1) + b1*rho(t, 3))
def m0(n): return np.where(n % 3 == 0, 0.5, b0)
print(f"b0 = {b0:.10f}")

# ---------- levels ----------
def build_levels(ells):
    levels = []; gaps = []
    for k, l in enumerate(ells, start=1):
        nk = int(np.ceil(1/l)); Pk = int(np.floor(LJ*2.0**(-k-2)/l))
        assert Pk >= 1, (k, l)
        y = J0 + (np.arange(1, Pk+1) - 0.5)*LJ/Pk
        j = np.rint(y*nk/(2*np.pi)).astype(np.int64)
        x = j*2*np.pi/nk
        keep = (x - 1.5*l > J0) & (x + 1.5*l < J1)
        if gaps:
            G = np.array(gaps)  # (a,b)
            lo = x - 1.5*l; hi = x + 1.5*l
            # interval intersection test, vectorised via sorting
            order = np.argsort(G[:, 0]); G = G[order]
            idx = np.searchsorted(G[:, 0], hi)  # gaps starting before hi
            # check the gap just before (max end among gaps starting before hi) -- gaps disjoint and sorted
            i = idx - 1
            bad = (i >= 0) & (G[np.maximum(i, 0), 1] > lo)
            keep &= ~bad
        xk = x[keep]; jk = j[keep]
        levels.append(dict(l=l, n=nk, P=Pk, x=xk, j=jk))
        gaps += [(a - l/2, a + l/2) for a in xk]
    return levels

def dm_formula(n, levels):
    n = np.asarray(n, np.int64)
    out = np.zeros(n.shape, float)
    for L in levels:
        # sum_x e^{-inx}, x = 2pi j/n_k  -> depends on n mod n_k
        if n.size > 200:
            ind = np.zeros(L['n']); np.add.at(ind, L['j'] % L['n'], 1.0)
            Sx = np.fft.fft(ind)[n % L['n']]          # sum_j e^{-2pi i r j/n_k}
        else:
            Sx = np.array([np.sum(np.exp(-2j*np.pi*(((nn % L['n'])*(L['j'] % L['n'])) % L['n'])/L['n'])) for nn in n.ravel()]).reshape(n.shape)
        out += eps*(2/L['P'])*np.real(hpsi(n*L['l'], L['l'])*Sx)
    return out

# =================== Test A ===================
print("\n=== Test A: ells = 0.08, 0.008, 0.0008 ===")
levA = build_levels([0.08, 0.008, 0.0008])
for k, L in enumerate(levA, 1):
    print(f" level {k}: l={L['l']}, n_k={L['n']}, P_k={L['P']}, kept={len(L['x'])}")
Lfft = 2**23
nmax = Lfft//2
n = np.arange(nmax+1)
dm = np.zeros(nmax+1); dm[1:] = dm_formula(n[1:], levA)
print(f" max|dm_n| = {np.max(np.abs(dm)):.4g}  (claimed bound 8 eps C0 = {8*eps*C0:.4g});  max|dm|/(eps C0) = {np.max(np.abs(dm))/(eps*C0):.4f}")
c = np.zeros(nmax+1, complex); c[1:] = 2*M*dm[1:]/(M + 1j*n[1:])
phi_c = Lfft*np.fft.irfft(c, Lfft)
t = 2*np.pi*np.arange(Lfft)/Lfft
# explicit physical solution
phi_p = np.zeros(Lfft); d2 = None
for L in levA:
    l = L['l']; a = eps/(L['P']*l)
    for xc in L['x']:
        # direct bump at x, u=(t-x)/l : phi = 2 M a l (th0'' - M l th0')(u)
        i0 = int(np.floor((xc - l/2)/(2*np.pi)*Lfft)); i1 = int(np.ceil((xc + l/2)/(2*np.pi)*Lfft))
        ii = np.arange(i0, i1+1) % Lfft; u = (t[ii] - xc)/l
        phi_p[ii] += 2*M*a*l*(th0_d(u, 2) - M*l*th0_d(u, 1))
        # reflected bump at -x, u=(-t-x)/l : phi = -2 M a l (th0'' + M l th0')(u)
        xr = (2*np.pi - xc)
        i0 = int(np.floor((xr - l/2)/(2*np.pi)*Lfft)); i1 = int(np.ceil((xr + l/2)/(2*np.pi)*Lfft))
        ii = np.arange(i0, i1+1) % Lfft; tt = t[ii]; u = (-(tt) - xc + 2*np.pi)/l
        phi_p[ii] += -2*M*a*l*(th0_d(u, 2) + M*l*th0_d(u, 1))
print(f" sup|phi| = {np.max(np.abs(phi_p)):.4g};  max|phi_coeff - phi_explicit| = {np.max(np.abs(phi_c-phi_p)):.3g}")
# S mask
inJ = (t > J0) & (t < J1)
inGap = np.zeros(Lfft, bool)
for L in levA:
    for xc in L['x']:
        inGap |= np.abs(t - xc) <= L['l']/2
S = inJ & ~inGap
lamS = S.mean()*2*np.pi
print(f" lambda(S)/|J| = {lamS/LJ:.4f} (>= 0.75 claimed)")
print(f" max_S |phi_coeff| = {np.max(np.abs(phi_c[S])):.3g}  vs sup|phi| = {np.max(np.abs(phi_c)):.3g}")
Fa = F0(t) + phi_c
print(f" max_S |F_0(closed form) + phi_coeff| = {np.max(np.abs(Fa[S])):.3g}")
# full F from all coefficients (F_0 part converges like 1/(N d): stay 0.01 away from jumps)
mall = np.zeros(nmax+1); mall[1:] = m0(n[1:]) + dm[1:]
call = np.zeros(nmax+1, complex); call[0] = 1.0; call[1:] = 2*M*mall[1:]/(M + 1j*n[1:])
Ffull = Lfft*np.fft.irfft(call, Lfft)
Sin = S & (t > J0 + 0.01) & (t < J1 - 0.01)
print(f" full coefficient synthesis (n<=4.2e6): max |F| on S (0.01 from ends of J) = {np.max(np.abs(Ffull[Sin])):.3g};"
      f" on gaps max|F| = {np.max(np.abs(Ffull[inGap & inJ])):.3g}")
print(f" weights: min m_n = {mall[1:].min():.6f}, max m_n = {mall[1:].max():.6f}")
# empty-interior proxy: largest run of S in J at this truncation vs spacing of next level
del Ffull, call, mall, phi_c, phi_p, c, dm

# =================== Test B ===================
print("\n=== Test B: ells = 0.08, 4e-5, 2e-8 (ratios ~1/2000) ===")
levB = build_levels([0.08, 4e-5, 2e-8])
for k, L in enumerate(levB, 1):
    print(f" level {k}: l={L['l']}, n_k={L['n']}, P_k={L['P']}, kept={len(L['x'])}, P'_k/P_k={len(L['x'])/L['P']:.5f}")
hp01 = -1j*hth1
for k in [2, 3]:
    L = levB[k-1]
    # nu_k-prediction: uniform measure on A_k = J minus gaps of earlier levels
    gaps = [(xc - LL['l']/2, xc + LL['l']/2) for LL in levB[:k-1] for xc in LL['x']]
    lamA = LJ - sum(b - a for a, b in gaps)
    def nuhat(r):
        if r == 0: return 1.0 + 0j
        tot = (np.exp(-1j*r*J1) - np.exp(-1j*r*J0))/(-1j*r)
        for a, b in gaps: tot -= (np.exp(-1j*r*b) - np.exp(-1j*r*a))/(-1j*r)
        return tot/lamA
    q = len(L['x'])/L['P']
    print(f" level {k}: q_k = P'_k/P_k = {q:.6f}, lambda(A_k)/|J| = {lamA/LJ:.6f}")
    print("   r   (m-m0)/eps [formula]   predicted 2q Re[hpsi0(1) nuhat(r)]   exact level-k term only   empirical nu_k")
    for r in range(-4, 5):
        nn = np.array([L['n'] + r])
        val = dm_formula(nn, levB)[0]/eps
        only = (2/L['P'])*np.real(hpsi(nn*L['l'], L['l'])*np.sum(np.exp(-2j*np.pi*(((nn[0] % L['n'])*(L['j'] % L['n'])) % L['n'])/L['n'])))[0]
        emp = np.mean(np.exp(-1j*r*L['x']))
        pred = 2*q*np.real(hp01*nuhat(r))
        pred_emp = 2*q*np.real(hp01*emp)
        print(f"  {r:+d}   {val:+.8e}        {pred:+.8e}                    {only:+.8e}        {pred_emp:+.8e}")
# limiting nu = 1_S/lambda(S): with all three levels removed
gaps = [(xc - LL['l']/2, xc + LL['l']/2) for LL in levB for xc in LL['x']]
lamS = LJ - sum(b - a for a, b in gaps)
print(f" lambda(S_3)/|J| = {lamS/LJ:.6f}")
