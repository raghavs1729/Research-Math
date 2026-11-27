"""Referee 2: independent check of TOY.md Theorem 2.3, step (i) and (ii).

Written from scratch (does not read or reuse toy_toeplitz*.py).

For random single-hole Toeplitz 0/1 sequences (several period structures, random
patterns, symmetric and not) we compute, by DIRECT smoothed Fourier summation
(flat-top C^infinity window, N = 2^20 coefficients, FFT evaluation):
  F'  (coefficients m'_n * 2M/(M+in), m' = m with entries on P Z set to 0),
  w'  (coefficients m'_n),
and compare at interior points of the level-i sub-arcs J_s (u in [0.25, 0.75])
with:
  (S1) symmetric case: F' = 1 + 2 B(s) e^{-eps u}, where B is the periodic solution of
       B(s) = e^{-eps} B(s-1) + eps C(s), C(s) := sum_{r mod P} p(r) e^{2 pi i r s/P}
       (atom masses computed by a direct DFT of p, NOT by the V-formula);
  (S2) recursion (R): C(s) = V_{i+1}(s mod q) whenever q does not divide s,
       V computed from the pattern v_{i+1};
  (N1) non-symmetric case: w' = 1/2 + sum_r p(r) sin((P/2-r)t)/sin(Pt/2)
       = 1/2 + Im[e^{i pi u} sum_r p_a(r) e^{-2 pi i r(s+u)/P}]/sin(pi u);
  (N2) the symmetry lemma: w'_{s-1}(u) - w'_s(u) = -Im[sum_r A_r e^{i w_r u}]/sin(pi u),
       A_r = p_a(r) e^{-2 pi i r s/P}(1 - e^{2 pi i r/P}), w_r = pi(1-2r/P);
       and this difference is not identically 0 when p_a != 0;
  (N3) in the non-symmetric case F' is NOT of the form 1 + 2B e^{-eps u} on sub-arcs;
  (P)  p_a == 0  iff  v_1..v_{i+1} symmetric;
  (ii) ||F''||^2 = sum_{n in PZ\\0} m 4M^2/(M^2+n^2) <= (4 pi^2/3) M^2/P^2.
"""
import numpy as np

rng = np.random.default_rng(20261127)
NC = 2 ** 20  # number of Fourier coefficients used


def window(x, x0=0.5):
    """Flat-top C^infinity window: 1 on [0, x0], 0 on [1, inf)."""
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    out[x <= x0] = 1.0
    mid = (x > x0) & (x < 1)
    y = (x[mid] - x0) / (1 - x0)
    g1 = np.exp(-1.0 / (1 - y))
    g2 = np.exp(-1.0 / y)
    out[mid] = g1 / (g1 + g2)
    return out


def toeplitz_seq(qs, vs, N):
    """m_n for n = 0..N (m_0 unused, set 0): m_n = v_{k(n)+1}(n/P_{k(n)} mod q_{k(n)+1})."""
    m = np.zeros(N + 1, dtype=np.int64)
    n = np.arange(1, N + 1)
    P = 1
    for q, v in zip(qs, vs):
        sel = (n % P == 0) & ((n // P) % q != 0)
        m[1:][sel] = np.asarray(v)[(n[sel] // P) % q]
        P *= q
        if P > N:
            break
    assert P > N, "not enough levels"
    return m


def smoothed_real_sum(coef, c0, L):
    """f(t_j) = c0 + 2 Re sum_{n=1}^{NC} coef[n] sigma(n/NC) e^{i n t_j}, t_j = 2 pi j/L."""
    x = np.zeros(L, complex)
    n = np.arange(1, NC + 1)
    x[1:NC + 1] = coef[1:NC + 1] * window(n / NC)
    S = L * np.fft.ifft(x)
    return c0 + 2 * S.real


def random_pattern(q, symmetric):
    v = rng.integers(0, 2, q)
    v[0] = 0
    if symmetric:
        for r in range(1, q):
            v[q - r] = v[r] if r <= q - r else v[q - r]
        v = np.array([0] + [v[min(r, q - r)] for r in range(1, q)])
    return v


def is_sym(v):
    q = len(v)
    return all(v[r] == v[(-r) % q] for r in range(q))


def run_case(name, qs, vs, i, M, out):
    """Level i (0-based as in TOY: P = P_{i+1}, q = q_{i+1})."""
    P = int(np.prod(qs[: i + 1]))
    q = qs[i]
    eps = 2 * np.pi * M / P
    N = NC
    m = toeplitz_seq(qs, vs, N + 1)
    n = np.arange(N + 2)
    # periodic part p (P-periodic extension of m'_n), check that m is P-periodic off PZ
    p = np.zeros(P, dtype=np.int64)
    p[1:] = m[1:P]
    offP = (n % P != 0) & (n >= 1)
    assert np.all(m[offP] == p[n[offP] % P]), "m not P-periodic off PZ"
    pa = (p - p[(-np.arange(P)) % P]) / 2.0
    sym_levels = all(is_sym(vs[k]) for k in range(i + 1))
    out.append(f"\n=== {name}: q_1.. = {[int(x) for x in qs[:i+3]]}..., level i = {i}, P = {P}, q = {q}, M = {M}, eps = {eps:.5f}")
    out.append(f"  v_1..v_{i+1} = {[list(map(int, vs[k])) for k in range(i + 1)]}")
    out.append(f"  patterns v_1..v_{i+1} symmetric: {sym_levels};  p_a == 0: {np.all(pa == 0)}   (P) {'OK' if sym_levels == bool(np.all(pa == 0)) else 'FAIL'}")

    L = P * int(2 ** np.ceil(np.log2((NC + 2) / P)))
    mp = np.where(n % P == 0, 0, m).astype(float)  # m'
    cF = np.zeros(NC + 1, complex)
    cw = np.zeros(NC + 1, complex)
    nn = np.arange(NC + 1)
    cF[1:] = mp[1:NC + 1] * 2 * M / (M + 1j * nn[1:])
    cw[1:] = mp[1:NC + 1]
    Fp = smoothed_real_sum(cF, 1.0, L)
    wp = smoothed_real_sum(cw, 0.5, L)
    per = L // P
    js = np.arange(per // 4, 3 * per // 4 + 1, max(1, per // 256))
    u = js / per
    # identical remainders: full F by direct sum, F'' := F - F' must be 2 pi/P-periodic
    cFull = np.zeros(NC + 1, complex)
    cFull[1:] = m[1:NC + 1] * 2 * M / (M + 1j * nn[1:])
    Ffull = smoothed_real_sum(cFull, 1.0, L)
    Fpp = Ffull - Fp
    res_rem = np.max(np.abs(np.roll(Fpp, -per) - Fpp))
    l2_num = np.mean(Fpp ** 2)

    # atom masses by direct DFT
    C = np.array([np.sum(p * np.exp(2j * np.pi * np.arange(P) * s / P)) for s in range(P)])
    V = np.array([np.sum(np.asarray(vs[i]) * np.exp(2j * np.pi * np.arange(q) * e / q)) for e in range(q)])
    res = {}
    if sym_levels:
        # periodic solution of B(s) = e^{-eps} B(s-1) + eps C(s)
        Brec = np.array([eps * np.sum(np.exp(-eps * np.arange(P)) * C[(s - np.arange(P)) % P])
                         for s in range(P)]) / (1 - np.exp(-eps * P))
        err_fit = 0.0
        for s in range(P):
            Bfit = (Fp[s * per + js] - 1) * np.exp(eps * u) / 2
            err_fit = max(err_fit, np.max(np.abs(Bfit - Brec[s].real)))
        res["S1 max|B_fit(s,u) - B_rec(s)|"] = err_fit
        res["S1 max|Im C(s)|"] = np.max(np.abs(C.imag))
        mask = np.array([s % q != 0 for s in range(P)])
        res["S2 max_{q∤s}|C(s) - V(s mod q)|"] = np.max(np.abs(C[mask] - V[np.arange(P)[mask] % q]))
        Rres = [Brec[s] - np.exp(-eps) * Brec[s - 1] - eps * V[s % q] for s in range(P) if s % q != 0]
        res["S2 max_{q∤s}|B(s)-e^-eps B(s-1)-eps V|"] = np.max(np.abs(Rres))
        # w' should be 1/2 on open sub-arcs
        err_w = max(np.max(np.abs(wp[s * per + js] - 0.5)) for s in range(P))
        res["S3 max|w' - 1/2| on sub-arcs"] = err_w
        # also: at q | s the jump involves lower levels (expected to differ from V)
        maskq = ~mask
        res["(info) max_{q|s}|C(s) - V(0)|"] = np.max(np.abs(C[maskq] - V[0]))
    else:
        errA = errB = 0.0
        for s in range(P):
            t = 2 * np.pi * (s + u) / P
            wA = 0.5 + np.sum(p[None, 1:] * np.sin((P / 2 - np.arange(1, P))[None, :] * t[:, None]), axis=1) / np.sin(P * t / 2)
            wB = 0.5 + np.imag(np.exp(1j * np.pi * u) * np.sum(pa[None, 1:] * np.exp(-2j * np.pi * np.arange(1, P)[None, :] * (s + u)[:, None] / P), axis=1)) / np.sin(np.pi * u)
            errA = max(errA, np.max(np.abs(wp[s * per + js] - wA)))
            errB = max(errB, np.max(np.abs(wp[s * per + js] - wB)))
        res["N1 max|w'_Fourier - formula sin/sin|"] = errA
        res["N1 max|w'_Fourier - formula Im[p_a]/sin(pi u)|"] = errB
        # symmetry lemma
        r = np.arange(1, P)
        om = np.pi * (1 - 2 * r / P)
        errL = 0.0
        mindiff = np.inf
        for s in range(1, P):
            A = pa[1:] * np.exp(-2j * np.pi * r * s / P) * (1 - np.exp(2j * np.pi * r / P))
            pred = -np.imag(np.sum(A[None, :] * np.exp(1j * om[None, :] * u[:, None]), axis=1)) / np.sin(np.pi * u)
            diff = wp[(s - 1) * per + js] - wp[s * per + js]
            errL = max(errL, np.max(np.abs(diff - pred)))
            mindiff = min(mindiff, np.max(np.abs(diff)))
        res["N2 max|w'_{s-1}-w'_s - (-Im sum A e^{iwu}/sin)|"] = errL
        res["N2 min_s max_u|w'_{s-1}-w'_s|  (must be > 0)"] = mindiff
        # F' not exponential: fit 1+2B e^{-eps u} at two points, measure misfit elsewhere
        mis = np.inf
        for s in range(P):
            Bf = (Fp[s * per + js] - 1) * np.exp(eps * u) / 2
            mis = min(mis, np.max(Bf) - np.min(Bf))
        res["N3 min_s (spread of (F'-1)e^{eps u}/2 over u)  (>0: not exponential)"] = mis
    # (ii)
    kk = np.arange(1, NC // P + 1)
    l2 = 2 * np.sum(m[P * kk] * 4 * M ** 2 / (M ** 2 + (P * kk) ** 2))
    res["(ii) ||F''||_2^2 (coefficient sum, truncated)"] = l2
    res["(ii) ||F - F'||_2^2 from the direct sums of F and F' (grid mean)"] = l2_num
    res["(ii) bound (4pi^2/3)M^2/P^2"] = 4 * np.pi ** 2 / 3 * M ** 2 / P ** 2
    res["(ii) max_t |F''(t + 2pi/P) - F''(t)|, F'' := F - F' (direct sums)"] = res_rem
    for k_, v_ in res.items():
        out.append(f"  {k_:62s} {v_:.3e}")
    return res


def main():
    out = ["Referee 2, r2_subarcs.py: direct smoothed Fourier sums vs the sub-arc formulas of TOY.md Thm 2.3.",
           f"N = {NC} coefficients, flat-top C^inf window, interior points u in [0.25,0.75]."]
    cases = []
    # symmetric examples
    for trial in range(4):
        qs = list(rng.choice([2, 3, 4, 5, 6, 7], size=30))
        vs = [random_pattern(q, True) for q in qs]
        i = int(rng.integers(1, 3))
        while int(np.prod(qs[: i + 1])) > 400:
            i -= 1
        M = float(rng.choice([0.5, 1.0, 3.0]))
        cases.append((f"symmetric #{trial}", qs, vs, i, M))
    # valuation family (period doubling)
    qs = [2] * 40
    vs = [np.array([0, 1 if (k + 1) % 2 == 0 else 0]) for k in range(40)]
    cases.append(("period doubling", qs, vs, 4, 0.5))
    # non-symmetric examples
    for trial in range(4):
        qs = list(rng.choice([3, 4, 5, 6, 7], size=30))
        vs = [random_pattern(q, rng.random() < 0.5) for q in qs]
        i = int(rng.integers(0, 2))
        k_bad = int(rng.integers(0, i + 1))
        while is_sym(vs[k_bad]):
            vs[k_bad] = random_pattern(qs[k_bad], False)
        M = float(rng.choice([0.5, 1.0, 3.0]))
        cases.append((f"non-symmetric #{trial} (asym at level {k_bad+1})", qs, vs, i, M))
    # asymmetry only at level 2 (v_1 symmetric), seen from level i = 1 and i = 2
    qs = [4, 5, 3] + [2] * 30
    vs = [np.array([0, 1, 0, 1]), np.array([0, 1, 1, 0, 0]), np.array([0, 1, 1])] + \
         [np.array([0, 1])] * 30
    cases.append(("asym only at level 2, seen at i=1", qs, vs, 1, 1.0))
    cases.append(("asym only at level 2, seen at i=2", qs, vs, 2, 1.0))
    cases.append(("asym only at level 2, seen at i=0 (symmetric there)", qs, vs, 0, 1.0))
    for c in cases:
        run_case(*c, out)
    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace(".py", "_results.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
