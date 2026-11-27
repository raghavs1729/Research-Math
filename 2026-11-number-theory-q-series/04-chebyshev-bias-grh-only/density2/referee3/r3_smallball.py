"""Referee 3: small-ball masses of F_M for Toeplitz weights (illustration of TOY.md Corollary 2.4).

AI-produced, 2026-11-27.

Not a proof of anything: a consistency test. For a single-hole Toeplitz sequence truncated at
depth K (period P = P_K), F = F' + F'' with F' = 1 + 2B(s)e^{-eps u} on J_s (eps = 2 pi M/P),
    B(s) = B(s-1)e^{-eps} + eps Wtot(s),   Wtot(s) = sum_{k : (P/P_k) | s} (P/P_k) V_k((s P_k/P) mod q_k)
(the general atom-mass recursion, checked independently in r3_toeplitz_forms.py), and
||F''||_2^2 <= (4 pi^2/3) M^2/P^2 (TOY (ii)). Hence, for every level kappa and radius eta,
    lam{|F'-kappa| < eta/2} - 4||F''||^2/eta^2  <=  lam{|F_M - kappa| < eta}  <=  lam{|F'-kappa| < 2 eta} + ||F''||^2/eta^2,
with lam normalised. Reported: the upper bound divided by 2 eta. A bounded ratio as eta -> 0 is what
"no atom at kappa" predicts; an atom makes it blow up like (atom mass)/(2 eta). The last case is a
periodic pattern with a flat sub-arc (found by r3_periodic_flat.py), where the ratio must blow up.
Also: the jumps of F' at 0, pi, pi/2, pi/4 for period doubling against TOY Proposition 2.1.
Output: r3_smallball_results.txt
"""
import numpy as np


def V_real(v):
    q = len(v)
    r = np.arange(q)
    return np.array([np.sum(v * np.cos(2 * np.pi * r * e / q)) for e in range(q)])


def B_array(qs, pats, M):
    P = int(np.prod(qs))
    Pk = np.cumprod([1] + list(qs))
    W = np.zeros(P)
    for k in range(1, len(qs) + 1):
        d = P // Pk[k]
        Vk = V_real(np.array(pats[k - 1], dtype=float))
        idx = np.arange(0, P, d)
        W[idx] += d * Vk[(idx // d) % qs[k - 1]]
    eps = 2 * np.pi * M / P
    j = np.arange(P)
    # periodic solution: B(0) = eps sum_j e^{-j eps} W(-j) / (1 - e^{-P eps})
    B0 = eps * np.sum(np.exp(-j * eps) * W[(-j) % P]) / (1 - np.exp(-P * eps))
    acc = np.cumsum(np.exp(j[1:] * eps) * W[1:])
    B = np.empty(P)
    B[0] = B0
    B[1:] = np.exp(-j[1:] * eps) * (B0 + eps * acc)
    return B, W, eps, P


def ball_measure(B, eps, kappa, eta):
    """Normalised measure of {u in [0,1): |1 - kappa + 2B e^{-eps u}| < eta}, averaged over s."""
    lo, hi = kappa - 1 - eta, kappa - 1 + eta          # need 2B y in (lo, hi), y in [e^{-eps}, 1]
    y0, y1 = np.exp(-eps), 1.0
    out = np.zeros_like(B)
    nz = np.abs(B) > 1e-300
    a = np.where(B[nz] > 0, lo / (2 * B[nz]), hi / (2 * B[nz]))
    b = np.where(B[nz] > 0, hi / (2 * B[nz]), lo / (2 * B[nz]))
    a = np.clip(a, y0, y1)
    b = np.clip(b, y0, y1)
    ln = np.where(b > a, (np.log(b) - np.log(a)) / eps, 0.0)
    out[nz] = ln
    out[~nz] = 1.0 if abs(1 - kappa) < eta else 0.0
    return out.mean()


def report(name, qs, pats, M, kappas, etas, lines, mult=1.0):
    B, W, eps, P = B_array(qs, pats, M)
    L2 = mult * (4 * np.pi ** 2 / 3) * M ** 2 / P ** 2
    lines.append(f"{name}: depth {len(qs)}, P = {P}, M = {M}, ||F''||^2 <= {L2:.1e}")
    for kappa in kappas:
        row = []
        for eta in etas:
            up = ball_measure(B, eps, kappa, 2 * eta) + L2 / eta ** 2
            row.append(up / (2 * eta))
        lines.append(f"    kappa = {kappa:+d}: upper bound of lam{{|F-kappa|<eta}}/(2 eta) for eta = "
                     + ", ".join(f"{e:.0e}" for e in etas) + ":  "
                     + "  ".join(f"{x:8.3f}" for x in row))
    return B, W, eps, P


def main():
    lines = []
    etas = [1e-1, 1e-2, 1e-3, 3e-4, 1e-4]
    K = 22
    pd_qs = [2] * K
    pd_pats = [[0, 1 if (k % 2 == 0) else 0] for k in range(1, K + 1)]   # c_k = 1_{k even}
    for M in (0.5, 1.0, 2.0):
        B, W, eps, P = report("period doubling", pd_qs, pd_pats, M, [-3, -1, 0, 1, 3], etas, lines)
        F = lambda s, u: 1 + 2 * B[s % P] * np.exp(-eps * u)
        jumps = []
        for a, sidx in ((0, 0), (0, P // 2), (1, P // 4), (2, P // 8), (3, P // 16)):
            jumps.append((sidx, F(sidx, 0.0) - F(sidx - 1, 1.0)))
        lines.append(f"    jumps of F' at t = 0, pi, pi/2, pi/4, pi/8: "
                     + ", ".join(f"{j:.6f}" for _, j in jumps)
                     + f"   (TOY Prop 2.1: 4 pi M/3 = {4*np.pi*M/3:.6f} at 0 and pi, "
                       f"+-(4 pi M/3) 2^-a at pi j/2^a)")
    rng = np.random.default_rng(5)
    for trial in range(3):
        qs, pats, P = [], [], 1
        while P * 5 < 3_000_000:
            q = int(rng.integers(2, 6))
            v = rng.integers(0, 2, size=q)
            v[0] = 0
            for r in range(1, q):
                v[q - r] = v[r]
            qs.append(q)
            pats.append(list(v))
            P *= q
        report(f"bounded q in {{2..5}}, random symmetric patterns (trial {trial})", qs, pats, 1.0,
               [-1, 0, 1], etas, lines)
    # adversarial: v_k = 1_{!=0} (q=2) except a non-constant symmetric q=4 pattern every 5th level
    qs, pats = [], []
    for k in range(1, 20):
        if k % 5 == 0:
            qs.append(4)
            pats.append([0, 1, 0, 1])
        else:
            qs.append(2)
            pats.append([0, 1])
    for M in (0.5, 1.7):
        report("near-periodic, bounded q (1_{!=0} except (0,1,0,1) at every 5th level)", qs, pats, M,
               [-1, 1], etas, lines)
    # control: periodic pattern with a flat sub-arc (Q = 5, m = 01101, M = 1.2455057521, level -1)
    m = [0, 1, 1, 0, 1]
    qs = [5] + [2] * 19
    pats = [[0] + m[:4]] + [[0, m[4]]] * 19
    report("control: periodic 01101 (flat sub-arc at level -1 for this M)", qs, pats, 1.2455057521,
           [-1, 1], etas, lines)
    txt = "\n".join(lines)
    print(txt)
    with open("r3_smallball_results.txt", "w") as f:
        f.write("r3_smallball.py results (AI-produced, 2026-11-27)\n" + txt + "\n")


if __name__ == "__main__":
    main()
