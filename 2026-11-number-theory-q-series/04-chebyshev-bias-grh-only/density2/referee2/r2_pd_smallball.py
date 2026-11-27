"""Referee 2: does the period-doubling F_M have an atom at +-1?  (consistency check)

The proof of TOY.md Thm 2.3 with q_k = 2 (no bad sub-arcs since q < 10) gives at a density
point, for all large i: V_{i+1}(1) = -c_{i+1} = n = (kappa-1)/2.  For period doubling
c_k = 1_{k even} alternates, so there is NO level set of positive measure at all (any M).
A crude smoothed-sum exploration (r2_misc.py) showed small-ball ratios at kappa = +-1
growing as delta decreases (M = 2).  Here we resolve this.

Method (referee's own formula, validated in r2_subarcs.py against direct Fourier sums):
at level i (P = 2^{i+1}), F = F'_i + F''_i with F'_i = 1 + 2 B(s) e^{-eps u} on J_s,
B the periodic solution of B(s) = e^{-eps} B(s-1) + eps C(s), C = DFT of p.
F''_i -> 0 uniformly (checked below by comparing two levels), so the law of F is approximated
by the law of F'_i with an explicit sup-error, and small balls can be probed far below the
resolution of the smoothed sums.
"""
import numpy as np
from scipy.signal import lfilter


def Fprime_values(i, M, us=(0.25, 0.5, 0.75)):
    P = 2 ** (i + 1)
    n = np.arange(P)
    n1 = n.copy()
    n1[0] = 1
    v2 = np.zeros(P, dtype=np.int64)
    t = n1.copy()
    mask = t % 2 == 0
    while mask.any():
        v2[mask] += 1
        t[mask] //= 2
        mask = t % 2 == 0
    p = (v2 % 2 == 1).astype(float)
    p[0] = 0.0
    C = (P * np.fft.ifft(p)).real  # C(s) = sum_r p(r) e^{2 pi i r s/P}, real (symmetric p)
    eps = 2 * np.pi * M / P
    a = np.exp(-eps)
    Blin = lfilter([eps], [1, -a], C)
    BPm1 = Blin[-1] / (1 - a ** P)
    B = Blin + a ** (n + 1) * BPm1
    return {u: 1 + 2 * B * np.exp(-eps * u) for u in us}, P


def main():
    out = ["Referee 2, r2_pd_smallball.py: period doubling, law of F_M near +-1."]
    for M in [0.5, 2.0]:
        I = 21
        vals, P = Fprime_values(I, M)
        # uniform closeness of levels: compare F'_I with F'_{I+2} at the same points
        vals2, P2 = Fprime_values(I + 2, M, us=(0.5,))
        # point t = 2 pi (s + 1/2)/P equals sub-arc 4s+2 of level I+2 at u = 0 ... use u = 0.5 there:
        # t = 2 pi (4 s + 2 + 1/2)/P2 is NOT the same point; compare at matching points instead:
        # level-I point (s, u) <-> level-(I+2) sub-arc 4s + floor(4u), u' = 4u - floor(4u).
        s = np.arange(P)
        u = 0.5
        s2 = 4 * s + int(np.floor(4 * u))
        up = 4 * u - np.floor(4 * u)
        v_lo = vals[u]
        vals2b, _ = Fprime_values(I + 2, M, us=(up,))
        v_hi = vals2b[up][s2]
        diff = np.max(np.abs(v_hi - v_lo))
        out.append(f"M={M}: level I={I} (P=2^{I+1}): sup over sub-arcs |F'_I - F'_(I+2)| = {diff:.2e}  "
                   f"(bounds sup|F''_I| up to the tail beyond I+2)")
        x = np.concatenate([vals[k] for k in vals])
        out.append(f"   range of F: [{x.min():.4f}, {x.max():.4f}], mean {x.mean():.6f}")
        for kappa in [-1.0, 1.0, 0.0]:
            row = f"   kappa={kappa:+.0f}: ratio lambda(|F-kappa|<d)/(2d) for d = "
            ds = [1e-1, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]
            rs = [np.mean(np.abs(x - kappa) < d) / (2 * d) for d in ds]
            row += ", ".join(f"{d:.0e}:{r:.3f}" for d, r in zip(ds, rs))
            out.append(row)
        # where are the values near -1 / +1?  fraction of sub-arcs with |F'+1| < 1e-3 by position
        near = np.abs(vals[0.5] + 1) < 1e-3
        out.append(f"   #sub-arcs with |F'+1| < 1e-3: {near.sum()} of {P}; their t/2pi in "
                   f"[{(np.where(near)[0].min()/P) if near.any() else float('nan'):.4f}, "
                   f"{(np.where(near)[0].max()/P) if near.any() else float('nan'):.4f}]")
    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace(".py", "_results.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
