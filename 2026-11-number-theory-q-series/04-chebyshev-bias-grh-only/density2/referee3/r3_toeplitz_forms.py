"""Referee 3: independent checks of TOY.md Theorem 2.3, step (i) and the symmetry step.

AI-produced, 2026-11-27.

Independent method. For a P-periodic one-sided sequence p (p(0) = 0) the function
    F'(t) = 1 + sum_{n != 0} p(|n|) (2M/(M+in)) e^{int} = 1 + 4M Re S(t),
    S(t)  = sum_{n>=1} p(n) e^{int}/(M+in),
is computed from the radial-integral representation of the resolvent
    S(t) = -i int_0^inf e^{iMx} R(e^{it-x}) dx,     R(z) = sum_{n>=1} p(n) z^n = sum_{r<P} p(r) z^r/(1-z^P)
(valid off (2pi/P)Z: Abel summation, then dominated convergence), and its t-derivative
    S'(t) = int_0^inf e^{iMx} (zR'(z))|_{z=e^{it-x}} dx.
None of the trigonometric manipulations of the proof is used. Then:
  (A) the two closed forms of w' on open sub-arcs (TOY (i)) against w'_num := (DF' + M F')/(2M);
  (B) symmetric patterns: F' = 1 + 2B(s)e^{-eps u} on J_s (B(s) constant in u), and the
      recursion (R) B(s) = B(s-1)e^{-eps} + eps V_{i+1}(s mod q) for q not dividing s;
  (C) non-symmetric patterns: the difference formula of the symmetry step,
      w'_{s-1}(u) - w'_s(u) = -Im[sum_r A_r e^{i w_r u}]/sin(pi u),
      and the fact that it is not identically 0 (so adjacent sub-arcs never agree);
  (D) a sanity test of the integral method against a known closed form (all weights 1);
  (E) the multilevel closed form of the author's toy_toeplitz.py (imported, not modified) against
      the integral method, on the two cases of toy_toeplitz_results.txt (there it was compared
      only with slowly converging partial Fourier sums, errors ~1e-4).
Output: r3_toeplitz_forms_results.txt
"""
import warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning

warnings.simplefilter('ignore', IntegrationWarning)  # roundoff warnings at the 1e-14 level

RNG = np.random.default_rng(20261127)


def toeplitz_values(qs, pats, N):
    """m_1..m_N for the single-hole Toeplitz sequence of TOY 2.2 (pats[k-1] = v_k on Z/q_k)."""
    P = [1]
    for q in qs:
        P.append(P[-1] * q)
    m = np.zeros(N + 1, dtype=np.int64)
    for n in range(1, N + 1):
        k = 0
        while k + 1 < len(P) and n % P[k + 1] == 0:
            k += 1
        if k >= len(qs):
            raise ValueError("n divisible by the deepest period")
        m[n] = pats[k][(n // P[k]) % qs[k]]
    return m, P


def periodic_part(qs, pats, L):
    """p(r), r = 0..P_L - 1: the P_L-periodic sequence m' of TOY (i) at level L = i+1."""
    m, P = toeplitz_values(qs[:L], pats[:L], P_of(qs[:L]) - 1)
    p = np.zeros(P[L], dtype=float)
    p[1:] = m[1:P[L]]
    return p


def P_of(qs):
    out = 1
    for q in qs:
        out *= q
    return out


def make_R(p):
    P = len(p)
    r = np.arange(P)

    def R(z):
        return np.sum(p * z ** r) / (1 - z ** P)

    def zRp(z):  # z R'(z) = sum_r p(r) [r z^r + (P-r) z^{r+P}]/(1-z^P)^2
        return np.sum(p * (r * z ** r + (P - r) * z ** (r + P))) / (1 - z ** P) ** 2
    return R, zRp


def S_and_dS(p, M, t):
    R, zRp = make_R(p)
    P = len(p)

    def cquad(f):
        brk = [0.0, 1.0 / P, 10.0 / P, 1.0, 5.0, 40.0]
        tot = 0.0 + 0.0j
        for a, b in zip(brk[:-1], brk[1:]):
            re = quad(lambda x: f(x).real, a, b, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
            im = quad(lambda x: f(x).imag, a, b, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
            tot += re + 1j * im
        return tot
    S = -1j * cquad(lambda x: np.exp(1j * M * x) * R(np.exp(1j * t - x)))
    dS = cquad(lambda x: np.exp(1j * M * x) * zRp(np.exp(1j * t - x)))
    return S, dS


def Fprime_num(p, M, t):
    S, dS = S_and_dS(p, M, t)
    F = 1 + 4 * M * S.real
    dF = 4 * M * dS.real
    return F, dF


def w_first(p, t):
    P = len(p)
    r = np.arange(1, P)
    return 0.5 + np.sum(p[1:] * np.sin((P / 2 - r) * t)) / np.sin(P * t / 2)


def w_second(p, s, u):
    P = len(p)
    r = np.arange(1, P)
    pa = (p[1:] - p[(P - r) % P]) / 2
    val = np.exp(1j * np.pi * u) * np.sum(pa * np.exp(-2j * np.pi * r * (s + u) / P))
    return 0.5 + val.imag / np.sin(np.pi * u)


def diff_formula(p, s, u):
    P = len(p)
    r = np.arange(1, P)
    pa = (p[1:] - p[(P - r) % P]) / 2
    A = pa * np.exp(-2j * np.pi * r * s / P) * (1 - np.exp(2j * np.pi * r / P))
    om = np.pi * (1 - 2 * r / P)
    return -np.sum(A * np.exp(1j * om * u)).imag / np.sin(np.pi * u), np.max(np.abs(A))


def V_of(v):
    q = len(v)
    e = np.arange(q)
    return np.array([np.sum(v * np.exp(2j * np.pi * np.arange(q) * ee / q)) for ee in e])


def rand_pattern(q, symmetric):
    v = RNG.integers(0, 2, size=q)
    v[0] = 0
    if symmetric:
        for r in range(1, q):
            v[q - r] = v[r]
    return v


def main():
    out = []
    M = 0.7
    # (D) sanity: all weights 1  =>  F = -1 + 4 pi M e^{-Mt}/(1-e^{-2 pi M}) on (0, 2pi)
    errD = 0.0
    for t in [0.3, 1.1, 2.0, 3.7, 5.9]:
        # p == 1 for all n >= 1: R(z) = z/(1-z); use P = 1 trick via a direct R
        R = lambda z: z / (1 - z)
        re = quad(lambda x: (np.exp(1j * M * x) * R(np.exp(1j * t - x))).real, 0, 60, limit=400,
                  epsabs=1e-14, epsrel=1e-13)[0]
        im = quad(lambda x: (np.exp(1j * M * x) * R(np.exp(1j * t - x))).imag, 0, 60, limit=400,
                  epsabs=1e-14, epsrel=1e-13)[0]
        S = -1j * (re + 1j * im)
        F = 1 + 4 * M * S.real
        exact = -1 + 4 * np.pi * M * np.exp(-M * t) / (1 - np.exp(-2 * np.pi * M))
        errD = max(errD, abs(F - exact))
    out.append(f"(D) integral method vs exact closed form (all weights 1, M={M}): max error {errD:.2e}")

    cases = [([3, 2, 4], True), ([5, 3], True), ([2, 2, 2, 3], True), ([7, 3], True),
             ([3, 4, 2, 5], True), ([2, 3, 2, 2, 7], True),
             ([4, 3], False), ([3, 5], False), ([2, 3, 2], False), ([6, 2], False),
             ([5, 4, 3], False), ([2, 7, 2, 3], False)]
    worst = {"first": 0.0, "second": 0.0, "B_const": 0.0, "R": 0.0, "diff": 0.0}
    min_diff_sup = np.inf
    for qs, sym in cases:
        pats = [rand_pattern(q, sym) for q in qs]
        # make sure non-symmetric cases are really non-symmetric at the top level
        if not sym:
            # force asymmetry at the deepest level with q >= 3 (every q = 2 pattern is symmetric)
            kk = max(k for k in range(len(qs)) if qs[k] >= 3)
            v = pats[kk]
            if all(v[r] == v[len(v) - r] for r in range(1, len(v))):
                v[1], v[len(v) - 1] = 1, 0
        L = len(qs)
        p = periodic_part(qs, pats, L)
        P = len(p)
        q = qs[-1]
        eps = 2 * np.pi * M / P
        a = RNG.integers(0, P // q)           # a level-(L-1) arc
        ss = [a * q + e for e in range(q)]
        us = [0.17, 0.43, 0.71, 0.88]
        e1 = e2 = 0.0
        Bs = {}
        Fvals = {}
        for s in ss:
            for u in us:
                t = 2 * np.pi * (s + u) / P
                F, dF = Fprime_num(p, M, t)
                Fvals[(s, u)] = F
                wnum = (dF + M * F) / (2 * M)
                e1 = max(e1, abs(wnum - w_first(p, t)))
                e2 = max(e2, abs(wnum - w_second(p, s, u)))
            if sym:
                Bu = [(Fvals[(s, u)] - 1) * np.exp(eps * u) / 2 for u in us]
                Bs[s] = np.mean(Bu)
                worst["B_const"] = max(worst["B_const"], max(Bu) - min(Bu))
        worst["first"] = max(worst["first"], e1)
        worst["second"] = max(worst["second"], e2)
        line = f"qs={qs} sym={sym} pats={[list(map(int, v)) for v in pats]} P={P}: " \
               f"|w'_num - first form| <= {e1:.1e}, |w'_num - second form| <= {e2:.1e}"
        if sym:
            V = V_of(pats[-1])
            errR = max(abs(Bs[s] - (Bs[s - 1] * np.exp(-eps) + eps * V[s % q].real))
                       for s in ss[1:])
            worst["R"] = max(worst["R"], errR)
            # general mass formula: B(s) - B(s-1)e^{-eps} = eps * sum_{k: (P/P_k) | s} (P/P_k) V_k(s P_k/P mod q_k)
            Pk = [1]
            for qq in qs:
                Pk.append(Pk[-1] * qq)
            def Wtot(s_):
                tot = 0.0
                for k in range(1, L + 1):
                    d = P // Pk[k]
                    if s_ % d == 0:
                        tot += d * V_of(pats[k - 1])[(s_ // d) % qs[k - 1]].real
                return tot
            sall = list(range(0, P))
            Ball = {}
            for s2 in sall:
                uu = 0.5
                F2, _ = Fprime_num(p, M, 2 * np.pi * (s2 + uu) / P)
                Ball[s2] = (F2 - 1) * np.exp(eps * uu) / 2
            errW = max(abs(Ball[s2] - (Ball[(s2 - 1) % P] * np.exp(-eps) + eps * Wtot(s2))) for s2 in sall)
            worst["W"] = max(worst.get("W", 0.0), errW)
            line += f"; general mass recursion (all s, incl. q|s) error {errW:.1e}"
            line += f"; B(s) const to {worst['B_const']:.1e}; (R) error {errR:.1e}; " \
                    f"V_q={np.round(V.real, 6).tolist()}"
        else:
            # difference formula on an adjacent pair inside the arc, checked against the numbers
            s = ss[1]
            ed = 0.0
            for u in us:
                t0 = 2 * np.pi * (s - 1 + u) / P
                t1 = 2 * np.pi * (s + u) / P
                F0, dF0 = Fprime_num(p, M, t0)
                F1, dF1 = Fprime_num(p, M, t1)
                d_num = (dF0 + M * F0) / (2 * M) - (dF1 + M * F1) / (2 * M)
                d_for, Amax = diff_formula(p, s, u)
                ed = max(ed, abs(d_num - d_for))
            worst["diff"] = max(worst["diff"], ed)
            # sup over u of the difference for every adjacent pair (must be > 0 if p_a != 0)
            ug = np.linspace(0.01, 0.99, 197)
            sups = []
            for s2 in range(1, P):
                sups.append(max(abs(diff_formula(p, s2, u)[0]) for u in ug))
            min_diff_sup = min(min_diff_sup, min(sups))
            r_ = np.arange(1, P)
            pa_max = np.max(np.abs(p[1:] - p[(P - r_) % P])) / 2
            line += f"; max|p_a| = {pa_max:.1f}; difference formula error {ed:.1e}; min over all " \
                    f"adjacent pairs of sup_u|w'_(s-1) - w'_s| = {min(sups):.3f}"
        out.append(line)
    # (E) the author's closed form (toy_toeplitz.py) against the integral method
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from toy_toeplitz import closed_form, sym_pattern, sequence
    rngA = np.random.default_rng(1)          # same seed and order as toy_toeplitz.main()
    Mm = 0.5
    qsA = [2] * 8
    patsA = [np.array([0, c], float) for c in [0, 1, 1, 0, 1, 0, 0, 1]]
    qsB = [3, 5, 2, 3, 4]
    patsB = [sym_pattern(q, rngA) for q in qsB]
    for name, qsE, patsE in (("valuation q=2, c=(0,1,1,0,1,0,0,1)", qsA, patsA),
                             ("random patterns q=(3,5,2,3,4)", qsB, patsB)):
        PE = int(np.prod(qsE))
        pE = np.zeros(PE)
        mE = sequence(PE - 1, qsE, patsE, len(qsE))
        pE[1:] = mE[1:PE]
        errE = 0.0
        for _ in range(25):
            sE = int(RNG.integers(0, PE))
            uE = float(RNG.uniform(0.2, 0.8))
            tE = 2 * np.pi * (sE + uE) / PE
            Fc = closed_form(np.array([tE]), Mm, qsE, patsE, len(qsE))[0]
            Fi, _ = Fprime_num(pE, Mm, tE)
            errE = max(errE, abs(Fc - Fi))
        out.append(f"(E) author's closed form (toy_toeplitz.py), {name}, P={PE}: max |closed form - "
                   f"integral method| over 25 points = {errE:.1e}")
    out.append("")
    out.append(f"SUMMARY: max |w'_num - first form| = {worst['first']:.1e}; "
               f"max |w'_num - second form| = {worst['second']:.1e}; "
               f"B(s) constant in u to {worst['B_const']:.1e}; max (R) error = {worst['R']:.1e}; "
               f"max difference-formula error = {worst['diff']:.1e}; "
               f"max general-mass-recursion error = {worst.get('W', 0.0):.1e}; "
               f"min sup-difference over non-symmetric cases = {min_diff_sup:.3f}")
    txt = "\n".join(out)
    print(txt)
    with open("r3_toeplitz_forms_results.txt", "w") as f:
        f.write("r3_toeplitz_forms.py results (AI-produced, 2026-11-27)\n" + txt + "\n")


if __name__ == "__main__":
    main()
