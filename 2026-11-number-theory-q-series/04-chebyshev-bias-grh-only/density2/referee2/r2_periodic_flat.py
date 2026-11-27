"""Referee 2: exhaustive search for flat pieces of F_M for PERIODIC 0/1 weights.

Independent derivation (referee's own): let m be Q-periodic on n >= 1, p(r) := m_r
(r = 1..Q-1), p(0) := m_Q.

* Symmetric case (p(-r) = p(r)): w = (1/2 - m_Q) + (2 pi/Q) sum_s c(s) delta_{2 pi s/Q},
  c(s) = sum_{r mod Q} p(r) e^{2 pi i r s/Q} (real).  With eps = 2 pi M/Q, on the open
  sub-arc (2 pi s/Q, 2 pi (s+1)/Q), t = 2 pi (s+u)/Q:
      F = (1 - 2 m_Q) + [2 eps/(1 - e^{-eps Q})] e^{-eps u} g_s(e^{-eps}),
      g_s(x) := sum_{j=0}^{Q-1} c(s-j) x^j.
  F is real-analytic on each sub-arc, so F has a level set of positive measure iff it is
  constant on some sub-arc, iff g_s(e^{-eps}) = 0 for some s: the level is then 1 - 2 m_Q.
  We find ALL roots x in (0,1) of all g_s, for all symmetric patterns with Q <= QMAX,
  i.e. all (pattern, M) with a flat piece, and record the levels.
* Non-symmetric case: F is constant on a sub-arc iff (D+M)F = M kappa there iff w is
  constant there (independent of M).  On open sub-arcs
      w = 1/2 - m_Q + sum_{r=1}^{Q-1} p(r) sin((Q/2 - r)t)/sin(Qt/2)
  (formula checked against direct Fourier sums in r2_subarcs.py).  We check that w is
  non-constant on every sub-arc for every non-symmetric pattern with Q <= QMAX_NS.
* A sample of the flat pieces found is re-verified by a direct smoothed Fourier sum of
  F_M (no closed form used), at the M found.
"""
import itertools
import numpy as np

QMAX = 22
QMAX_NS = 12
NC = 2 ** 19


def window(x, x0=0.5):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    out[x <= x0] = 1.0
    mid = (x > x0) & (x < 1)
    y = (x[mid] - x0) / (1 - x0)
    g1 = np.exp(-1.0 / (1 - y))
    g2 = np.exp(-1.0 / y)
    out[mid] = g1 / (g1 + g2)
    return out


def sym_patterns(Q):
    half = Q // 2
    for bits in itertools.product([0, 1], repeat=half + 1):
        p = np.zeros(Q, dtype=np.int64)
        p[0] = bits[0]
        for r in range(1, half + 1):
            p[r] = bits[r]
            p[Q - r] = bits[r]
        yield p


def minimal_period(p):
    Q = len(p)
    for d in range(1, Q + 1):
        if Q % d == 0 and np.all(p == np.roll(p, d)):
            return d
    return Q


def direct_F(p, M, L):
    """Smoothed direct Fourier sum of F_M for the Q-periodic weights p (p[0] = m_Q)."""
    Q = len(p)
    n = np.arange(1, NC + 1)
    coef = p[n % Q] * 2 * M / (M + 1j * n) * window(n / NC)
    x = np.zeros(L, complex)
    x[1:NC + 1] = coef
    return 1.0 + 2 * (L * np.fft.ifft(x)).real


def main():
    out = ["Referee 2, r2_periodic_flat.py: flat pieces of F_M for periodic 0/1 weights.",
           f"Symmetric patterns: all Q <= {QMAX} (primitive period = Q only); non-symmetric: all Q <= {QMAX_NS}."]
    flats = []
    npat = 0
    for Q in range(2, QMAX + 1):
        for p in sym_patterns(Q):
            if minimal_period(p) != Q or not p.any():
                continue
            npat += 1
            c = np.array([np.sum(p * np.cos(2 * np.pi * np.arange(Q) * s / Q)) for s in range(Q)])
            scale = np.sum(np.abs(c))
            for s in range(Q):
                g = np.array([c[(s - j) % Q] for j in range(Q)])  # coefficients of x^j
                gh = g[::-1].copy()  # highest degree first
                # g(1) = Q m_Q: if m_Q = 0, x = 1 (i.e. M = 0, excluded) is a root: deflate it
                h = gh.copy()
                while len(h) > 1 and abs(np.polyval(h, 1.0)) < 1e-9 * scale:
                    h, _ = np.polydiv(h, np.array([1.0, -1.0]))
                if len(h) <= 1:
                    continue
                roots = np.roots(h)
                for x in roots:
                    if abs(x.imag) < 1e-6 and 1e-12 < x.real < 1 - 1e-9:
                        xr = x.real
                        dg = np.polyder(gh)
                        for _ in range(6):  # Newton refinement on g itself
                            d = np.polyval(dg, xr)
                            if d != 0:
                                xr -= np.polyval(gh, xr) / d
                        if not (0 < xr < 1 - 1e-9):
                            continue
                        M = -Q * np.log(xr) / (2 * np.pi)
                        resid = abs(np.polyval(gh, xr)) / scale
                        flats.append((Q, tuple(int(v) for v in p), s, M, 1 - 2 * int(p[0]), resid))
    levels = sorted(set(f[4] for f in flats))
    out.append(f"symmetric primitive nonzero patterns examined: {npat}")
    out.append(f"(pattern, M, sub-arc) with a flat piece: {len(flats)}")
    out.append(f"levels of flat pieces found: {levels}   (prediction: 1 - 2 m_Q in {{1, -1}})")
    out.append(f"max relative residual |g_s(e^-eps)|/sum|c| at the (Newton-refined) roots: {max(f[5] for f in flats):.2e}")
    out.append(f"smallest M with a flat piece: {min(f[3] for f in flats):.4g}")
    for lev in levels:
        sub = [f for f in flats if f[4] == lev]
        out.append(f"  level {lev:+d}: {len(sub)} cases; M range [{min(f[3] for f in sub):.4g}, {max(f[3] for f in sub):.4g}]")
    # how often does a pattern have flat pieces for some M?
    pats = sorted(set((f[0], f[1]) for f in flats))
    out.append(f"distinct patterns having a flat piece for some M > 0: {len(pats)} of {npat}")
    out.append("first examples (Q, p(0..Q-1) with p(0)=m_Q, sub-arc s, M, level):")
    for f in flats[:12]:
        out.append(f"  Q={f[0]:2d} p={''.join(map(str, f[1]))} s={f[2]:2d} M={f[3]:.6f} level={f[4]:+d}")

    # direct-Fourier re-verification of a sample
    out.append("\nDirect smoothed Fourier sum of F_M at the M found (N = 2^19), on the flat sub-arc and")
    out.append("on a neighbouring sub-arc (u in [0.2,0.8]):")
    rng = np.random.default_rng(7)
    sample = [flats[k] for k in rng.choice(len(flats), size=min(10, len(flats)), replace=False)]
    for (Q, p, s, M, lev, _) in sample:
        p = np.array(p)
        L = Q * int(2 ** np.ceil(np.log2((NC + 2) / Q)))
        F = direct_F(p, M, L)
        per = L // Q
        js = np.arange(int(0.2 * per), int(0.8 * per))
        vals = F[s * per + js]
        nb = F[((s + 1) % Q) * per + js]
        out.append(f"  Q={Q:2d} p={''.join(map(str, p))} s={s:2d} M={M:.5f}: max|F - ({lev:+d})| on flat arc = "
                   f"{np.max(np.abs(vals - lev)):.2e};  spread on next arc = {np.ptp(nb):.3e}")

    # non-symmetric: w non-constant on every sub-arc
    out.append(f"\nNon-symmetric periodic patterns, Q <= {QMAX_NS}: min over patterns and sub-arcs of")
    out.append("the spread of w over u in [0.05,0.95] (must be > 0 for no flat piece at any M):")
    u = np.linspace(0.05, 0.95, 181)
    for Q in range(3, QMAX_NS + 1):
        worst = np.inf
        cnt = 0
        for bits in itertools.product([0, 1], repeat=Q):
            p = np.array(bits)
            if all(p[r] == p[(-r) % Q] for r in range(Q)):
                continue
            cnt += 1
            r = np.arange(1, Q)
            for s in range(Q):
                t = 2 * np.pi * (s + u) / Q
                w = 0.5 - p[0] + np.sum(p[None, 1:] * np.sin((Q / 2 - r)[None, :] * t[:, None]), axis=1) / np.sin(Q * t / 2)
                worst = min(worst, np.ptp(w))
        out.append(f"  Q={Q:2d}: {cnt:5d} non-symmetric patterns, min spread = {worst:.3e}")
    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace(".py", "_results.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
