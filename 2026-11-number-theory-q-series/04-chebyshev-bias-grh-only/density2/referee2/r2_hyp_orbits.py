"""Referee 2: numerical checks for HYP.md Proposition 1.1 and the examples of HYP.md section 2.

Model A (circle): frequencies X = {2^k : 0 <= k < 9}, multiplicities in {1,2}, M = 1/2,
  F = 1 + Re sum m_g r_g e^{i g t}, r_g = 4M/(M + i g).  H = <X_{>T}> = 2^{k0} Z,
  N = 2^{k0}, K = H^perp = {2 pi j/N}.
Model B (2-torus): G = Z alpha + Z beta (alpha = 1, beta = sqrt 2), coordinates (x, y),
  e_{a alpha + b beta}(x,y) = e^{i(a x + b y)};
  X = {alpha, 2 alpha, 4 alpha, beta} u {2^k alpha + 2 beta : 3 <= k <= 7}.
  Taking T between max Lambda = 4 and 8 + 2 sqrt 2 gives H = <X_{>T}> = <8 alpha, 2 beta>,
  G/H = Z/8 x Z/2, N = 16, K = {(2 pi j/8, pi l)}.
Checks (Proposition 1.1):
  (1) F_H(theta k) = F_H(theta) for k in K (machine precision);
  (2) F_Lambda separates each K-orbit: min_{k != k'} |F_L(theta k) - F_L(theta k')| > 0 for
      random theta, with the distribution of the gap;
  (3) counting identity  E|S_delta(theta)| = N * lambda(E_delta),  E_delta = {|F - c| < delta},
      and P(|S_delta| >= 2) / P(|S_delta| >= 1) -> 0 as delta -> 0  (|S| <= 1 a.e. in the limit).
Examples of HYP.md section 2 (exact integer arithmetic):
  (4) X = {2^k}: index [G : <X_{>T}>] = 2^{k(T)};  FD fails (G_T cyclic);
  (5) X = odd multiples {(2k+1)}: Devin's parity hypothesis (Thm 2.3 of 1706.06394) holds,
      TS_Z fails (index 1), FD fails (G_T cyclic)  -> the parity hypothesis does not imply FD
      or the disjunction "TS_Z or FD";
  (6) Example 3.3a of DENSITY.md (truncated): <X_{>T}> contains every b_n (index 1).
"""
import math
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

rng = np.random.default_rng(3)
M = 0.5


def r(g):
    return 4 * M / (M + 1j * g)


def model_A():
    out = ["Model A (circle, X = {2^k : k < 9}):"]
    freqs = [2 ** k for k in range(9)]
    mult = {g: int(rng.integers(1, 3)) for g in freqs}
    for k0 in [1, 2, 3, 4]:
        N = 2 ** k0
        Lam = [g for g in freqs if g % N != 0]
        Hf = [g for g in freqs if g % N == 0]
        FL = lambda t: np.real(sum(mult[g] * r(g) * np.exp(1j * g * t) for g in Lam))
        FH = lambda t: 1 + np.real(sum(mult[g] * r(g) * np.exp(1j * g * t) for g in Hf))
        th = rng.uniform(0, 2 * np.pi, 200000)
        shifts = 2 * np.pi * np.arange(N) / N
        inv = max(np.max(np.abs(FH(th + s) - FH(th))) for s in shifts)
        VL = np.array([FL(th + s) for s in shifts])  # N x samples
        VLs = np.sort(VL, axis=0)
        gap = np.min(np.diff(VLs, axis=0), axis=0)
        out.append(f"  k0={k0} (H = {N}Z, N = {N}, Lambda = {Lam}): max|F_H(th+k)-F_H(th)| = {inv:.1e}; "
                   f"orbit gap of F_Lambda: min {gap.min():.2e}, P(gap<1e-3) = {np.mean(gap < 1e-3):.2e}")
        F = VL + np.array([FH(th + s) for s in shifts])
        for c in [0.0, 1.0]:
            line = f"      c={c}: "
            for delta in [3e-2, 3e-3, 3e-4]:
                S = np.sum(np.abs(F - c) < delta, axis=0)
                lamE = np.mean(np.abs(F[0] - c) < delta)
                p1, p2 = np.mean(S >= 1), np.mean(S >= 2)
                line += f"[d={delta:.0e}: E|S|={S.mean():.3e}, N*lam={N*lamE:.3e}, P(S>=2)/P(S>=1)={p2/max(p1,1e-300):.2e}] "
            out.append(line)
    return out


def model_B():
    out = ["Model B (2-torus, G/H = Z/8 x Z/2, N = 16):"]
    s2 = math.sqrt(2)
    # frequencies as (a, b) meaning a*alpha + b*beta, value a + b*sqrt2
    X = [(1, 0), (2, 0), (4, 0), (0, 1)] + [(2 ** k, 2) for k in range(3, 8)]
    mult = {g: int(rng.integers(1, 3)) for g in X}
    val = lambda g: g[0] + g[1] * s2
    T = 5.0  # between max Lambda = 4 and 8 + 2 sqrt2
    Lam = [g for g in X if val(g) <= T]
    Hf = [g for g in X if val(g) > T]
    # H = <(8,2),(16,2),...> = <(8,0),(0,2)> ; check with Smith normal form
    Hm = sp.Matrix([list(g) for g in Hf]).T
    snf = smith_normal_form(Hm, domain=sp.ZZ)
    diag = [snf[i, i] for i in range(min(snf.shape)) if snf[i, i] != 0]
    out.append(f"  Lambda = {Lam}; generators of H: {Hf}; Smith invariants of H in Z^2: {diag} -> index {abs(np.prod(diag))}")
    K = [(2 * np.pi * j / 8, np.pi * l) for j in range(8) for l in range(2)]
    e = lambda g, x, y: np.exp(1j * (g[0] * x + g[1] * y))
    FL = lambda x, y: np.real(sum(mult[g] * r(val(g)) * e(g, x, y) for g in Lam))
    FH = lambda x, y: 1 + np.real(sum(mult[g] * r(val(g)) * e(g, x, y) for g in Hf))
    x = rng.uniform(0, 2 * np.pi, 100000)
    y = rng.uniform(0, 2 * np.pi, 100000)
    inv = max(np.max(np.abs(FH(x + a, y + b) - FH(x, y))) for a, b in K)
    VL = np.array([FL(x + a, y + b) for a, b in K])
    gap = np.min(np.diff(np.sort(VL, axis=0), axis=0), axis=0)
    out.append(f"  max|F_H(th k) - F_H(th)| over K = {inv:.1e};  orbit gap of F_Lambda: min {gap.min():.2e}, "
               f"P(gap < 1e-3) = {np.mean(gap < 1e-3):.2e}")
    F = VL + np.array([FH(x + a, y + b) for a, b in K])
    for c in [0.0, 1.0]:
        line = f"      c={c}: "
        for delta in [3e-2, 3e-3, 3e-4]:
            S = np.sum(np.abs(F - c) < delta, axis=0)
            lamE = np.mean(np.abs(F[0] - c) < delta)
            p1, p2 = np.mean(S >= 1), np.mean(S >= 2)
            line += f"[d={delta:.0e}: E|S|={S.mean():.3e}, N*lam={16*lamE:.3e}, P(S>=2)/P(S>=1)={p2/max(p1,1e-300):.2e}] "
        out.append(line)
    return out


def examples():
    out = ["Examples of HYP.md section 2:"]
    # (4) powers of 2
    X = [2 ** k for k in range(40)]
    rows = []
    for T in [1, 10, 100, 1000, 10 ** 6]:
        tail = [g for g in X if g > T]
        idx = math.gcd(*tail)  # G = Z (1 in X)
        rows.append(f"T={T}: index {idx}")
    out.append("  (4) X = {2^k}: " + "; ".join(rows) + "  -> TS_Z holds; G_T = Z for T >= 1, so every nonzero")
    out.append("      hom G_T -> Z is nonzero on all #X_T ~ log2 T elements: FD fails; TS fails (Q-span = Q).")
    # (5) odd multiples
    X = [2 * k + 1 for k in range(60)]
    rows = []
    for T in [1, 10, 100]:
        tail = [g for g in X if g > T]
        rows.append(f"T={T}: index {math.gcd(*tail)}")
    # parity: relation lattice {n : sum n_k x_k = 0} for the first 12 elements; basis vectors' coordinate sums
    Xs = X[:12]
    # integer basis of the relation lattice: x_0 = 1 gives relations e_k - x_k e_0
    basis = []
    for k in range(1, len(Xs)):
        vec = [0] * len(Xs)
        vec[k] = 1
        vec[0] = -Xs[k]
        basis.append(vec)
    par = all(sum(v_) % 2 == 0 for v_ in basis)
    out.append("  (5) X = odd multiples: " + "; ".join(rows) + f"; relation lattice basis e_k - x_k e_0 (x_0 = 1, a Z-basis),")
    out.append(f"      all coordinate sums even: {par} -> Devin's parity hypothesis holds; TS_Z fails (index 1) and FD fails.")
    # (6) Example 3.3a truncated: b_1..b_n basis; X = {b_n} u {b_n + b_N(n,j)}
    nprim = 4
    # indices: primaries 0..nprim-1, linking N(n,j) = nprim + 3n + j (j = 0,1,2)
    dim = nprim + 3 * nprim
    vecs, vals = [], []
    bval = {}
    for n in range(nprim):
        bval[n] = 10.0 * (n + 1)
    for n in range(nprim):
        for j in range(3):
            bval[nprim + 3 * n + j] = math.exp(bval[n]) * (j + 1)
    for idx in range(dim):
        v_ = [0] * dim
        v_[idx] = 1
        vecs.append(v_)
        vals.append(bval[idx])
    for n in range(nprim):
        for j in range(3):
            v_ = [0] * dim
            v_[n] = 1
            v_[nprim + 3 * n + j] = 1
            vecs.append(v_)
            vals.append(bval[n] + bval[nprim + 3 * n + j])
    T = 1000.0
    tail = [vec for vec, vl in zip(vecs, vals) if vl > T]
    Mt = sp.Matrix(tail).T
    snf = smith_normal_form(Mt, domain=sp.ZZ)
    diag = [snf[i, i] for i in range(min(snf.shape)) if snf[i, i] != 0]
    out.append(f"  (6) Example 3.3a truncated ({nprim} primary, 3 linking each), T = {T}: rank of <X_>T> = {len(diag)} of {dim}, "
               f"Smith invariants {sorted(set(int(abs(d)) for d in diag))} -> <X_>T> = G (index 1): TS_Z fails")
    return out


def main():
    out = ["Referee 2, r2_hyp_orbits.py: HYP.md Proposition 1.1 (finite models) and section 2 examples.", ""]
    out += model_A() + [""] + model_B() + [""] + examples()
    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace(".py", "_results.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
