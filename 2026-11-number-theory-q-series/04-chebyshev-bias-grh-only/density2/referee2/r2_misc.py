"""Referee 2: elementary inequalities of TOY.md Thm 2.3 (iii), (v), and a small-ball exploration.

(1) Counting in (iii): if lambda(E cap I) >= 0.99 lambda(I) and I is the union of q equal
    sub-arcs, the number b of sub-arcs with lambda(E cap J) < 0.9 lambda(J) satisfies
    b < 0.1 q, hence b <= floor(0.1 q); and 2 floor(0.1 q) < q - 1,
    (q-1) - 2 floor(0.1 q) >= 0.6 (q-1), for all q >= 2 (checked q <= 10^6).
    Sharpness of b <= floor(0.1 q): an extremal configuration is exhibited.
(2) rho(u) = eps e^{-eps u}/(1 - e^{-eps}) in [1 - eps, 1 + eps] for u in [0,1], eps in (0,1].
(3) (v): for kappa not odd, d = dist(kappa-1, 2Z): min over integers n of
    min_u |kappa - 1 - 2 n rho(u)| >= d/2 for eps <= eps_0(kappa); we compute the largest
    admissible eps on a grid of kappa.
(4) Exploration (not a proof): small-ball ratios lambda{|F_M - kappa| < delta}/(2 delta) for
    the period-doubling weights (covered by Thm 2.3) and for m_n = [nu_2(n) + nu_3(n) odd]
    (NOT covered), by a direct smoothed Fourier sum.  An atom of mass alpha at kappa would
    make the ratio grow like alpha/(2 delta).
"""
import numpy as np

out = ["Referee 2, r2_misc.py"]

# (1)
q = np.arange(2, 10 ** 6 + 1)
fl = np.floor(0.1 * q + 1e-12).astype(np.int64)
ok1 = np.all(2 * fl < q - 1)
ok2 = np.all((q - 1) - 2 * fl >= 0.6 * (q - 1) - 1e-9)
# Markov: b bad sub-arcs each missing > 0.1 of its mass, total missing <= 0.01 q (in sub-arc units)
# => 0.1 b < 0.01 q => b < 0.1 q => b <= ceil(0.1 q) - 1 <= floor(0.1 q).
ok3 = np.all(np.ceil(0.1 * q - 1e-12) - 1 <= fl)
out.append(f"(1) 2 floor(.1q) < q-1: {ok1};  (q-1)-2floor(.1q) >= .6(q-1): {ok2};  ceil(.1q)-1 <= floor(.1q): {ok3}  (q <= 1e6)")
# extremal: q = 25, b = 2 bad sub-arcs each with E-density just below 0.9, others full
qq = 25
b = int(np.ceil(0.1 * qq) - 1)
miss = b * 0.1000001
out.append(f"    extremal example q={qq}: b={b} bad sub-arcs missing {miss:.4f} sub-arc units <= 0.01 q = {0.01*qq:.4f}: {miss <= 0.01*qq}")
# (2)
eps = np.linspace(1e-6, 1, 4001)
u = np.linspace(0, 1, 401)
rho = eps[:, None] * np.exp(-eps[:, None] * u[None, :]) / (1 - np.exp(-eps[:, None]))
ok = np.all(rho <= 1 + eps[:, None] + 1e-12) and np.all(rho >= 1 - eps[:, None] - 1e-12)
out.append(f"(2) rho in [1-eps, 1+eps] on eps in (0,1], u in [0,1]: {ok};  max (rho-1)/eps = {np.max((rho-1)/eps[:,None]):.4f}, "
           f"min (rho-1)/eps = {np.min((rho-1)/eps[:,None]):.4f}")
# (3)
lines = []
for kappa in [0.0, 0.5, 2.0, -2.0, 1.9, 2.1, 0.999, -3.5]:
    d = abs((kappa - 1) - 2 * np.round((kappa - 1) / 2))
    best = None
    for e in np.geomspace(1e-5, 1, 400)[::-1]:
        uu = np.linspace(0, 1, 201)
        rr = e * np.exp(-e * uu) / (1 - np.exp(-e))
        nmax = int(abs(kappa) / (2 * (1 - e)) + 10 / e + 10) if e < 1 else int(abs(kappa) + 20)
        ns = np.arange(-nmax, nmax + 1)
        val = np.min(np.abs(kappa - 1 - 2 * ns[:, None] * rr[None, :]))
        if val >= d / 2:
            best = e
            break
    lines.append(f"kappa={kappa:+.3f}: d={d:.3f}, largest eps on grid with min_n,u|kappa-1-2n rho| >= d/2: {best:.3g}")
out.append("(3) " + "\n    ".join(lines))


# (4) exploration
def window(x, x0=0.5):
    out_ = np.zeros_like(x)
    out_[x <= x0] = 1.0
    mid = (x > x0) & (x < 1)
    y = (x[mid] - x0) / (1 - x0)
    g1, g2 = np.exp(-1.0 / (1 - y)), np.exp(-1.0 / y)
    out_[mid] = g1 / (g1 + g2)
    return out_


def nu(n, p):
    n = n.copy()
    v = np.zeros_like(n)
    mask = n % p == 0
    while mask.any():
        v[mask] += 1
        n[mask] //= p
        mask = (n % p == 0)
    return v


NC = 2 ** 21
L = 2 ** 23
n = np.arange(1, NC + 1)
v2, v3 = nu(n, 2), nu(n, 3)
seqs = {"period doubling [nu2 odd]": (v2 % 2 == 1).astype(float),
        "[nu2 + nu3 odd] (not covered)": ((v2 + v3) % 2 == 1).astype(float)}
out.append("(4) small-ball ratios lambda{|F_M - kappa| < delta}/(2 delta) (direct smoothed sums, N = 2^21, grid 2^23):")
for name, mseq in seqs.items():
    for M in [0.5, 2.0]:
        x = np.zeros(L, complex)
        x[1:NC + 1] = mseq * 2 * M / (M + 1j * n) * window(n / NC)
        F = 1 + 2 * (L * np.fft.ifft(x)).real
        row = f"    {name:32s} M={M}: "
        for kappa in [-1.0, 0.0, 1.0]:
            rs = [np.mean(np.abs(F - kappa) < dl) / (2 * dl) for dl in [1e-1, 1e-2, 3e-3]]
            row += f"kappa={kappa:+.0f}: " + ", ".join(f"{r_:.3f}" for r_ in rs) + " | "
        out.append(row)
out.append("    (columns: delta = 1e-1, 1e-2, 3e-3; growth like 1/delta would indicate an atom)")
out.append("    CAUTION: the apparent growth at kappa = +-1 for period doubling, M = 2, is NOT an atom: "
           "r2_pd_smallball.py evaluates F to sup-error < 1e-6 at level 2^22 and finds no mass within "
           "3e-4 of -1 and a bounded ratio at +1 (lumpy, cluster-like law).")

txt = "\n".join(out)
print(txt)
with open(__file__.replace(".py", "_results.txt"), "w") as f:
    f.write(txt + "\n")
