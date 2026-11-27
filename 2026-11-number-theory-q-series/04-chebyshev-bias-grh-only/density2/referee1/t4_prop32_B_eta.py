"""Referee test 4 (independent): Proposition 3.2 of TOY.md and the remarks after it.

 (a) the identity Psi = F/(2M) + I/2 + const (Psi = t/2 + 2 sum m_n sin(nt)/n, I = int_0^t F),
     exactly for random 0/1 truncations;
 (b) sum_{n>=1} sin(nt)/n = (pi - t)/2 on (0, 2pi), hence B_eta = Psi - pi/2 there;
 (c) on a flat arc {F = c} (the examples of test 2), B_eta is affine with slope c/2
     (Prop 3.2 is the case c = 0), consistent with the mechanism of the proof;
 (d) quadratic variation: int |B(t+h) - B(t)|^2 dt/2pi = pi h/2 - h^2/4 for every +-1 sequence;
 (e) the remark "For periodic eta the slope is +-1/2": true when eta_r = eta_{q-r} (symmetric),
     FALSE otherwise: for eta = (1,-1,-1) repeated (n = 1,2,0 mod 3), B_eta' = 1/2 + 1/(1+2cos t)
     off 2pi Z/3, and B_eta has logarithmic singularities at 2pi/3, 4pi/3.
"""
import math

import numpy as np
from scipy.integrate import quad

out = []


def log(s=""):
    print(s)
    out.append(s)


TWO_PI = 2 * math.pi
rng = np.random.default_rng(2026)

# (a)
N = 3000
for M in [0.5, 2.0]:
    m = rng.integers(0, 2, N).astype(float)
    n = np.arange(1, N + 1)
    a = 2 * M / (M + 1j * n)
    ts = rng.uniform(0, TWO_PI, 400)
    E = np.exp(1j * np.outer(ts, n))
    F = 1 + 2 * np.real(E @ (m * a))
    I = ts + 2 * np.real(((E - 1) / (1j * n)) @ (m * a))
    Psi = ts / 2 + 2 * (np.sin(np.outer(ts, n)) @ (m / n))
    d = Psi - F / (2 * M) - I / 2
    log(f"(a) M={M}: spread of Psi - F/(2M) - I/2 over 400 random t (N = {N}): {d.max() - d.min():.2e}; constant = {d.mean():.6f} (= -1/(2M) + sum m_n a_n/(in) term)")

# (b)
NN = 10 ** 6
ts = np.linspace(0.3, TWO_PI - 0.3, 7)
n = np.arange(1, NN + 1)
S = np.array([np.sum(np.sin(n * t) / n) for t in ts])
log(f"(b) max |sum_(n<=1e6) sin(nt)/n - (pi-t)/2| on [0.3, 2pi-0.3]: {np.max(np.abs(S - (math.pi - ts) / 2)):.2e}")


# (c) slope c/2 on flat arcs: closed-form examples from test 2 (verified there)
def beta_of(a):
    q = len(a)
    j = np.arange(q)
    B = (np.exp(2j * np.pi * np.outer(j, j) / q) @ np.array(a, dtype=float)) / q
    return B.real


def B_partial(eta_fun, t, N=10 ** 6):
    n = np.arange(1, N + 1)
    e = np.array([eta_fun(k) for k in range(1, N + 1)], dtype=float) if N <= 2000 else eta_fun(n)
    return np.array([np.sum(e * np.sin(n * tt) / n) for tt in np.atleast_1d(t)])


for (a, k, M, c) in [((0, 1, 0, 0, 0, 0, 0, 1), 3, 0.441271, 1), ((1, 2), 1, math.log(3) / math.pi, -1),
                     ((1, 0, 1, 1, 0, 1, 1, 0), 2, 0.902669, -1)]:
    q = len(a)
    if max(a) > 1:
        # B_eta is defined for 0/1 weights; for weights (1,2) use V_m = sum m_n sin(nt)/n directly:
        # Psi = t/2 + 2 V_m, and Psi' = w; on the flat arc Psi has slope c/2.
        n = np.arange(1, 10 ** 6 + 1)
        mn = np.array(a, dtype=float)[n % q]
        tA, tB = TWO_PI * (k + 0.2) / q, TWO_PI * (k + 0.8) / q
        PsiA = tA / 2 + 2 * np.sum(mn * np.sin(n * tA) / n)
        PsiB = tB / 2 + 2 * np.sum(mn * np.sin(n * tB) / n)
        log(f"(c) weights {a}, M={M:.6f}, flat arc {k} at level {c}: slope of Psi on the arc = {(PsiB - PsiA) / (tB - tA):.6f} (c/2 = {c / 2})")
    else:
        eta = lambda n_: 2 * np.array(a, dtype=float)[np.asarray(n_) % q] - 1
        tA, tB = TWO_PI * (k + 0.2) / q, TWO_PI * (k + 0.8) / q
        BA, BB = B_partial(eta, [tA, tB])
        log(f"(c) 0/1 weights {a}, M={M:.6f}, flat arc {k} at level {c}: slope of B_eta on the arc = {(BB - BA) / (tB - tA):.6f} (c/2 = {c / 2})")

# (d) quadratic variation via Parseval (exact for the truncation) and via a grid
Nq = 2 ** 15
eta = rng.choice([-1.0, 1.0], Nq)
n = np.arange(1, Nq + 1)
for h in [0.01, 0.3, 1.0, 3.0, 6.0]:
    pars = 2 * np.sum(np.sin(n * h / 2) ** 2 / n ** 2)
    log(f"(d) h={h}: 2 sum sin^2(nh/2)/n^2 (N=2^15) = {pars:.6f}; pi h/2 - h^2/4 = {math.pi * h / 2 - h * h / 4:.6f} (tail <= 2/N = {2 / Nq:.1e})")
# grid check of one h with an FFT-built B
K = 2 ** 18
coef = np.zeros(K, dtype=complex)
coef[1:Nq + 1] = eta / (2j * n)
coef[-Nq:] = -(eta / (2j * n))[::-1]
Bgrid = np.real(np.fft.ifft(coef) * K)
shift = 2 ** 13  # h = 2pi * shift / K
hh = TWO_PI * shift / K
qv = np.mean((np.roll(Bgrid, -shift) - Bgrid) ** 2)
log(f"(d) grid: h = {hh:.6f}: mean |B(t+h)-B(t)|^2 = {qv:.6f}; formula {math.pi * hh / 2 - hh * hh / 4:.6f}")

# (e) periodic eta: symmetric vs non-symmetric
log()
log("(e) slope of B_eta for periodic eta")
eta_sym = np.array([1.0, 1.0, -1.0, -1.0, 1.0])  # residues 0..4 mod 5: eta_1=eta_4=1, eta_2=eta_3=-1, eta_5(=eta_0)=1
eta_s = lambda n_: eta_sym[np.asarray(n_) % 5]
t1, t2 = 0.5, 1.0
B1, B2 = B_partial(eta_s, [t1, t2])
log(f"  symmetric q=5, eta_(1..5) = (1,-1,-1,1,1): slope on (0, 2pi/5) = {(B2 - B1) / (t2 - t1):.6f} (predicted -eta_q/2 = -0.5)")
eta_ns = np.array([-1.0, 1.0, -1.0])  # residues 0,1,2 mod 3: eta_1 = 1, eta_2 = -1, eta_3 = -1
eta_n = lambda n_: eta_ns[np.asarray(n_) % 3]
for (t1, t2) in [(0.3, 0.8), (1.0, 1.5), (2.3, 2.6), (3.0, 3.5)]:
    B1, B2 = B_partial(eta_n, [t1, t2])
    pred = quad(lambda s: 0.5 + 1 / (1 + 2 * math.cos(s)), t1, t2)[0]
    log(f"  non-symmetric q=3, eta_(1,2,3) = (1,-1,-1): [B(t2)-B(t1)]/(t2-t1) on ({t1},{t2}) = {(B2 - B1) / (t2 - t1):+.6f};"
        f" int(1/2 + 1/(1+2cos))/(t2-t1) = {pred / (t2 - t1):+.6f}")
t0 = TWO_PI / 3
for dlt in [1e-1, 1e-2, 1e-3]:
    Ba, Bb = B_partial(eta_n, [t0 - dlt, t0 - 10 * dlt], N=4 * 10 ** 6)
    log(f"  near 2pi/3: B(t0 - d) - B(t0 - 10 d) with d={dlt:g}: {Ba - Bb:.4f} (log singularity predicts ln(10)/sqrt(3) = {math.log(10) / math.sqrt(3):.4f})")

with open(__file__.replace(".py", "_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
