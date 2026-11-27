#!/usr/bin/env python3
"""
rank_one_deficit.py -- [num] checks for MAHLER.md Sections 4-5 (deficit criterion, the rank-one wall,
and the answers to questions (1)-(4) of the brief).

Cyclic (rank-one) finite models at M = 1/2:  F_N(t) = 1 + sum_{n<=N} (C_n e^{int} + conj),
C_n = w_n 2M/(M+in), with weights w_n of several kinds.  For each model we print
  delta_N   = -m(F_N)  (Mahler deficit; m = mean of log|F_N| over the circle),
  J_N       = log(1/|C_top|)  (the Jensen-Newton upper bound for delta_N, top = largest n with w_n != 0),
  m(2w_N), log||2w_N||_2  (w_N = 1/2 + sum w_n (e^{int}+e^{-int}); for integer w_n, 2w_N is an integer
               Laurent polynomial),
  rms of w_N on the near-zero set {|F_N| <= q} (q = a quantile of |F_N|), against ||w_N||_2.
Models: (a) w_n = 1 on 25N (the rescaled Thm 9.6 spectrum, M = 1/2, sigma^2 < 2), (b) random 0/1,
(c) the base-5 Toeplitz sequence of DENSITY.md Thm 9.15, (d) Ex 9.2 (real weights, F = 0 on an arc),
(e) beam-search 0/1 near-atoms (density2/toy_beam2.py, eta = 0.15, 0.30, N = 512; stored in
beam_best_eta*.npy), (f) Prop 7.3 (real weights, |F_100| < 2e-15 on [2, 2.3]; values on the arc in
50-digit arithmetic).
Also: the LI model on the true zeros of L(s, chi_{-4}) (prototype/explore/chi4_zeros.txt, 511 zeros
up to 640): Monte Carlo m(F_T - c) on T^n and rms of w_T on {|F_T| < q}.
Output: rank_one_deficit_results.txt
"""
import math, os
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
out = open(os.path.join(HERE, "rank_one_deficit_results.txt"), "w")
def log(s=""):
    print(s); out.write(s + "\n"); out.flush()

M = 0.5
def C_of(w):
    n = np.arange(len(w)); C = np.zeros(len(w), complex); C[1:] = w[1:] * 2 * M / (M + 1j * n[1:]); return C
def on_grid(coef, K, const, shift):
    """const + sum_{n>=1} (coef_n e^{in(t+shift)} + conj) on t = 2 pi k / K"""
    a = np.zeros(K, complex); n = np.arange(len(coef))
    a[:len(coef)] = coef * np.exp(1j * n * shift); a[0] = 0
    return const + 2 * (np.fft.ifft(a) * K).real

def report(name, w, q_eta=0.10, K=None, Fvals_override=None):
    N = len(w) - 1; top = int(np.nonzero(w)[0].max())
    K = K or max(2 ** 16, 64 * N)
    C = C_of(w); sh = 0.3183 * 2 * math.pi / K
    F = on_grid(C, K, 1.0, sh) if Fvals_override is None else Fvals_override
    W = on_grid(w.astype(complex), K, 0.5, sh)
    aF = np.abs(F)
    mF = float(np.mean(np.log(np.maximum(aF, 1e-300))))
    J = math.log(1 / abs(C[top]))
    q = float(np.quantile(aF, q_eta)); S = aF <= q
    m2w = float(np.mean(np.log(np.maximum(np.abs(2 * W), 1e-300))))
    l2w = 0.5 * math.log(float(np.mean((2 * W) ** 2)))
    log(f"  {name:28s} N={N:6d} delta_N={-mF:7.3f}  J_N={J:7.3f}  delta/logN={-mF / math.log(N):6.3f}  "
        f"m(2w)={m2w:6.3f} log||2w||={l2w:6.3f}  q_{q_eta:.2f}={q:.2e} rms(w|S)={math.sqrt(float(np.mean(W[S] ** 2))):7.3f} "
        f"||w||={math.sqrt(float(np.mean(W ** 2))):7.2f}")

log("Rank-one models at M = 1/2.  delta_N = -m(F_N) <= J_N = log(1/|C_top|) (Lemma 2.1).")
log("(a) w = 1 on 25N")
for N in [500, 5000, 50000]:
    w = np.zeros(N + 1); w[25::25] = 1; report("1 on 25N", w)
rng = np.random.default_rng(5)
log("(b) random 0/1")
for N in [128, 1024, 8192]:
    w = rng.integers(0, 2, N + 1).astype(float); w[0] = 0; w[N] = 1; report("random 0/1", w)
log("(c) base-5 Toeplitz: w_n = 1 iff last nonzero base-5 digit of n is 1 or 4")
def last_digit5(n):
    while n % 5 == 0: n //= 5
    return n % 5
for N in [125, 3125, 15625]:
    w = np.array([0.0] + [1.0 if last_digit5(n) in (1, 4) else 0.0 for n in range(1, N + 1)]); report("Toeplitz base 5", w)
log("(d) Example 9.2 (real weights b0 = 0.66722, 1/2; F = 0 on (4pi/3, 2pi))")
A0 = 2 * math.pi / (1 - math.exp(-2 * math.pi * M)); A1 = (2 * math.pi / 3) * math.exp(4 * math.pi * M / 3) / (1 - math.exp(-2 * math.pi * M / 3))
b0 = A1 / (2 * (A1 - A0)); b1 = 0.5 - b0
for N in [30, 300, 3000, 30000]:
    n = np.arange(N + 1); w = np.where(n % 3 == 0, b0 + b1, b0).astype(float); w[0] = 0
    report("Ex 9.2 (real)", w, q_eta=0.30)
log("(e) beam-search 0/1 near-atoms (density2 toy_beam2.py, width 150)")
for fn, eta in [("beam_best_eta0.15_N512.npy", 0.15), ("beam_best_eta0.3_N512.npy", 0.30)]:
    p = os.path.join(HERE, fn)
    if os.path.exists(p):
        m01 = np.load(p).astype(float); w = np.r_[0.0, m01]; report(f"beam eta={eta}", w, q_eta=eta)
log("(f) Prop 7.3 (real weights, N = 100): values on J = [2, 2.3] recomputed with 50 digits")
cn = {}
for L in open(os.path.join(REPO, "monotonicity", "density", "model_M0.5_coeffs.txt")):
    if L.strip():
        a, b = L.split(); cn[int(a)] = float.fromhex(b)
N = 100; w = np.zeros(N + 1)
for n_, c_ in cn.items(): w[n_] = c_ * (M * M + n_ * n_) / (4 * M * M)
K = 2 ** 16; t = 2 * math.pi * np.arange(K) / K + 0.3183 * 2 * math.pi / K
F = on_grid(C_of(w), K, 1.0, 0.3183 * 2 * math.pi / K)
J = np.nonzero((t >= 2) & (t <= 2.3))[0]
mp.mp.dps = 50
cnm = {n_: mp.mpf(cn[n_]) for n_ in cn}
for k in J:
    tt = mp.mpf(float(t[k]))
    F[k] = float(1 + mp.fsum(cnm[n_] * (mp.cos(n_ * tt) + 2 * n_ * mp.sin(n_ * tt)) for n_ in cnm))
log(f"     max|F_100| on J (50 digits) = {np.abs(F[J]).max():.3e}, |J|/2pi = {len(J) / K:.4f}")
report("Prop 7.3 (real)", w, q_eta=0.04, K=K, Fvals_override=F)
# The grid value of delta above is biased (|F| ~ 1e-15 on J, 13 sign changes there; the off-J values use
# weights rounded to double).  Exact value by Jensen's formula from certified roots (python-flint, Arb) of
# z^100 F_100(z) with the exact coefficients C_n = c_n (1/2 - i n) (c_n the stored doubles, read as exact).
try:
    from flint import acb_poly, acb, arb, ctx
    ctx.prec = 256
    co = [acb(0)] * 201; co[100] = acb(1)
    for n_, c_ in cn.items():
        Cn = acb(arb(c_)) * acb(arb(0.5), arb(-n_))
        co[100 + n_] = Cn; co[100 - n_] = Cn.conjugate()
    rts = acb_poly(co).roots(tol=1e-40, maxprec=4096)
    mex = math.log(float(abs(co[200]).mid())) + sum(max(0.0, math.log(float(abs(z).mid()))) for z in rts)
    log(f"     certified-roots value: delta_100 = -m(F_100) = {-mex:.6f}  (delta/logN = {-mex / math.log(100):.3f}; "
        f"use this, not the grid value)")
except ImportError:
    log("     (python-flint not installed: certified-roots value skipped)")
log()
log("LI model on the true zeros of L(s, chi_-4) (all zeros simple; under LI the hull is T^n).")
gam = np.array([float(x) for x in open(os.path.join(REPO, "prototype", "explore", "chi4_zeros.txt")) if x.strip()])
rngmc = np.random.default_rng(11)
log(f"  {len(gam)} zeros, gamma_1 = {gam[0]:.6f}; vertex bound (every ordinate is a vertex under LI): "
    f"log(2M/|M + i gamma_1|) = {math.log(2 * M / abs(M + 1j * gam[0])):+.4f}")
for T in [50, 100, 200, 400, 640]:
    g = gam[gam <= T]; n = len(g)
    amp = 2 * 2 * M / np.abs(M + 1j * g)            # F = 1 + sum amp_j cos(theta_j)
    acc = {0.0: 0.0, 1.0: 0.0, 2.0: 0.0}; tot = 0; wS = []; Fs = []
    for chunk in range(20):
        th = rngmc.uniform(0, 2 * math.pi, size=(50000, n))
        Fv = 1 + np.cos(th) @ amp
        Wv = 0.5 + 2 * np.cos(th + np.angle(M + 1j * g)) @ np.ones(n)    # w = 1/2 + sum 2 cos(theta_j + arg(M+i g))
        for cl in acc: acc[cl] += float(np.sum(np.log(np.abs(Fv - cl))))
        tot += len(Fv); q = 0.02
        sel = np.abs(Fv) < q; wS.append(Wv[sel]); Fs.append(Wv)
    wS = np.concatenate(wS); Fs = np.concatenate(Fs)
    V = float(np.sum(amp ** 2) / 2)
    log(f"  T={T:4d} n={n:4d} V={V:.4f}  m(F_T - c), c = 0, 1, 2: " + ", ".join(f"{acc[c] / tot:+.4f}" for c in acc) +
        f"   rms w_T on {{|F_T|<0.02}} = {math.sqrt(float(np.mean(wS ** 2))) if len(wS) else float('nan'):.2f} "
        f"(#={len(wS)}), ||w_T||_2 = {math.sqrt(float(np.mean(Fs ** 2))):.2f}")
out.close()
