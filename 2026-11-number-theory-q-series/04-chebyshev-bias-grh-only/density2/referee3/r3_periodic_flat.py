"""Referee 3: flat pieces of F_M for periodic weights (TOY.md, Remarks on Corollary 2.4 (1)).

AI-produced, 2026-11-27.

For a 0/1 sequence (m_n)_{n>=1} of minimal period Q, symmetric (m_r = m_{Q-r}), put p(r) = m_r
(p(0) = m_Q) and W(s) = sum_{r mod Q} p(r) cos(2 pi r s/Q). Between the atoms of w,
    F_M = (1 - 2 m_Q) + 2 B(s) e^{-eps u}  on J_s (eps = 2 pi M/Q),  B(s) = eps f_s(e^{-eps})/(1-e^{-Q eps}),
    f_s(x) = sum_{j=0}^{Q-1} x^j W(s-j),
so F_M is flat on J_s iff f_s has a root x = e^{-eps} in (0,1), and then the level is 1 - 2 m_Q.
(The derivation is the same as TOY (i)/(R); script r3_toeplitz_forms.py checks (R) independently.)

This script
 (a) counts primitive symmetric patterns with Q <= 22 and the triples (pattern, sub-arc, M)
     with a flat sub-arc (quoted in TOY: 11,956 patterns and 108,205 triples, from another
     referee's script which this referee did not read), and how many patterns have at least one;
 (b) re-checks several flat cases with an independent evaluation of F_M (radial-integral
     representation of the resolvent, as in r3_toeplitz_forms.py), confirming flatness and the
     level 1 - 2 m_Q;
 (c) for non-symmetric periodic patterns, checks on random examples that the analytic
     continuation of F_M from one sub-arc never matches the next one (no flat piece).
Output: r3_periodic_flat_results.txt
"""
import itertools
import warnings
import numpy as np
import mpmath as mp
from scipy.integrate import quad, IntegrationWarning

warnings.simplefilter("ignore", IntegrationWarning)


def primitive_period(seq):
    Q = len(seq)
    for d in range(1, Q + 1):
        if Q % d == 0 and all(seq[i] == seq[i % d] for i in range(Q)):
            return d
    return Q


def symmetric_patterns(Q):
    """All symmetric 0/1 sequences of period Q, as tuples (m_1, ..., m_Q)."""
    free = list(range(1, Q // 2 + 1))          # residues r in 1..Q-1 up to r <-> Q-r
    out = []
    for vals in itertools.product((0, 1), repeat=len(free)):
        m = [0] * (Q + 1)
        for r, x in zip(free, vals):
            m[r] = x
            m[Q - r] = x
        for mQ in (0, 1):
            m[Q] = mQ
            out.append(tuple(m[1:]))
    return out


n_multiple = [0]


def deflate_at_one(c, tol):
    """Divide sum_j c_j x^j by (x - 1) while x = 1 is a root (|value| <= tol). Returns (c, mult)."""
    mult = 0
    while len(c) > 1 and abs(sum(c)) <= tol * max(1, max(abs(x) for x in c)):
        # synthetic division of the polynomial (low -> high) by (x - 1)
        hi = list(c[::-1])                     # high -> low
        out = [hi[0]]
        for a in hi[1:-1]:
            out.append(a + out[-1])
        c = out[::-1]
        mult += 1
    return c, mult


def hp_roots_in_01(m, s, dps=60):
    """Distinct real roots of f_s in (0,1), in dps-digit arithmetic (x = 0, x = 1 removed)."""
    Q = len(m)
    with mp.workdps(dps):
        p = [mp.mpf(m[Q - 1])] + [mp.mpf(x) for x in m[:Q - 1]]
        W = [mp.fsum(p[r] * mp.cos(2 * mp.pi * r * ss / Q) for r in range(Q)) for ss in range(Q)]
        coef = [W[(s - j) % Q] for j in range(Q)]
        small = mp.mpf(10) ** (-(dps - 15))
        coef = [c if abs(c) > small else mp.mpf(0) for c in coef]
        while len(coef) > 1 and coef[-1] == 0:
            coef = coef[:-1]
        while len(coef) > 1 and coef[0] == 0:
            coef = coef[1:]
        coef, _ = deflate_at_one(coef, small)
        if len(coef) <= 1:
            return []
        rts = mp.polyroots(coef[::-1], maxsteps=2000, extraprec=10 * dps)
        out = []
        for z in rts:
            if abs(mp.im(z)) < mp.mpf(10) ** (-25) and 0 < mp.re(z) < 1:
                x = float(mp.re(z))
                if any(abs(x - y) < 1e-12 for y in out):
                    n_multiple[0] += 1
                else:
                    out.append(x)
        return sorted(out)


n_close = [0]


def flat_roots(m):
    """Distinct roots x in (0,1) of f_s for every sub-arc s, as (s, x); also the number of
    polynomials re-solved in high precision and the number of grid sign changes."""
    Q = len(m)
    p = np.array([m[Q - 1]] + list(m[:Q - 1]), dtype=float)   # p(0) = m_Q, p(r) = m_r
    r = np.arange(Q)
    W = np.array([np.sum(p * np.cos(2 * np.pi * r * s / Q)) for s in range(Q)])
    xg = np.linspace(1e-6, 1 - 1e-6, 4001)
    res, amb, signch = [], 0, 0
    for s in range(Q):
        c = [W[(s - j) % Q] for j in range(Q)]
        c = [x if abs(x) > 1e-12 else 0.0 for x in c]
        while len(c) > 1 and c[-1] == 0.0:
            c = c[:-1]
        while len(c) > 1 and c[0] == 0.0:
            c = c[1:]
        if len(c) <= 1:
            continue
        c, mult1 = deflate_at_one(c, 1e-11)
        vals = np.polyval(np.array(c[::-1]), xg)
        signch += int(np.sum(np.sign(vals[1:]) * np.sign(vals[:-1]) < 0))
        if len(c) <= 1:
            continue
        roots = np.roots(np.array(c[::-1]))
        good = [z for z in roots if abs(z.imag) < 1e-9 and 1e-7 < z.real < 1 - 1e-7]
        grey = [z for z in roots if (abs(z.imag) < 1e-3 and -1e-3 < z.real < 1 + 1e-3)
                and not (abs(z.imag) < 1e-9 and 1e-7 < z.real < 1 - 1e-7)
                and not (abs(z.imag) > 1e-9 and False)]
        # a genuinely complex root far from [0,1] is not grey; a real root clearly inside is good
        grey = [z for z in grey if not (abs(z.imag) >= 1e-3)]
        # close pairs of real roots may be one exact multiple root (seen: 2cos(4pi/9), Q = 18)
        gx = sorted(z.real for z in good)
        if any(b - a < 1e-5 for a, b in zip(gx[:-1], gx[1:])):
            grey = grey + [None]
            n_close[0] += 1
        if grey:
            amb += 1
            for x in hp_roots_in_01(m, s):
                res.append((s, x))
            continue
        for z in good:
            res.append((s, float(z.real)))
    merged = []
    for s, x in sorted(res):
        if merged and merged[-1][0] == s and abs(merged[-1][1] - x) < 1e-9:
            continue
        merged.append((s, x))
    return merged, amb, signch


def F_full(m, M, t):
    """Independent evaluation of F_M(t) for the periodic sequence m (period Q = len(m))."""
    Q = len(m)
    mm = np.array(m, dtype=float)
    rr = np.arange(1, Q + 1)

    def R(z):
        return np.sum(mm * z ** rr) / (1 - z ** Q)
    brk = [0.0, 1.0 / Q, 10.0 / Q, 1.0, 5.0, 40.0]
    tot = 0.0 + 0.0j
    for a, b in zip(brk[:-1], brk[1:]):
        re = quad(lambda x: (np.exp(1j * M * x) * R(np.exp(1j * t - x))).real, a, b, limit=400,
                  epsabs=1e-14, epsrel=1e-13)[0]
        im = quad(lambda x: (np.exp(1j * M * x) * R(np.exp(1j * t - x))).imag, a, b, limit=400,
                  epsabs=1e-14, epsrel=1e-13)[0]
        tot += re + 1j * im
    S = -1j * tot
    return 1 + 4 * M * S.real


def main():
    lines = []
    n_pat = 0
    n_trip = 0
    n_with = 0
    n_amb = 0
    n_sc = 0
    levels = {}
    per_Q = []
    examples = []
    seen_keys = set()
    for Q in range(1, 23):
        cnt_pat = cnt_trip = cnt_with = 0
        for m in symmetric_patterns(Q):
            if primitive_period(m) != Q:
                continue
            if Q == 1 and m[0] == 0:
                cnt_pat += 1     # m = 0: F = 1 identically (flat everywhere, every M); not counted
                continue
            cnt_pat += 1
            roots, amb, sc = flat_roots(m)
            n_amb += amb
            n_sc += sc
            if roots:
                cnt_with += 1
                lev = 1 - 2 * m[-1]
                levels[lev] = levels.get(lev, 0) + len(roots)
                key = (lev, Q)
                if Q in (5, 9, 14, 22) and key not in seen_keys:
                    seen_keys.add(key)
                    examples.append((m, roots[0]))
            cnt_trip += len(roots)
        per_Q.append((Q, cnt_pat, cnt_with, cnt_trip))
        n_pat += cnt_pat
        n_trip += cnt_trip
        n_with += cnt_with
    lines.append("(a) primitive symmetric 0/1 patterns, Q <= 22 (m = 0 at Q = 1 counted as a pattern, "
                 "excluded from the triples):")
    for Q, a, b, c in per_Q:
        lines.append(f"    Q={Q:2d}: patterns {a:5d}, with a flat sub-arc {b:5d}, triples {c:6d}")
    lines.append(f"    TOTAL: {n_pat} patterns (Q = 2..22: {n_pat - 2}); {n_with} have at least one flat "
                 f"sub-arc for some M; {n_trip} triples (pattern, sub-arc, M); levels {levels}; "
                 f"polynomials re-solved in 60-digit arithmetic (a root near the real segment [0,1] that double "
                 f"precision could not classify, or a close pair of real roots): {n_amb}, of which close pairs "
                 f"{n_close[0]}; exact multiple roots found (counted once): {n_multiple[0]}; sign changes of f_s "
                 f"on a 4001-point grid (odd-multiplicity roots only): {n_sc}")
    # (b) independent confirmation of some flat cases
    lines.append("(b) independent evaluation of F_M on flat sub-arcs (radial integral):")
    for m, (s, x) in examples:
        Q = len(m)
        M = -Q * np.log(x) / (2 * np.pi)
        lev = 1 - 2 * m[-1]
        devs = []
        for u in (0.15, 0.4, 0.6, 0.85):
            t = 2 * np.pi * (s + u) / Q
            devs.append(abs(F_full(m, M, t) - lev))
        # and on the next sub-arc (should not be flat in general)
        t2 = 2 * np.pi * ((s + 1) % Q + 0.5) / Q
        lines.append(f"    m={''.join(map(str, m))} (Q={Q}), sub-arc s={s}, M={M:.10f}: "
                     f"max |F_M - ({lev})| on J_s = {max(devs):.1e}; on J_(s+1) at u=1/2: "
                     f"F_M - ({lev}) = {F_full(m, M, t2) - lev:.3e}")
    # (c) non-symmetric periodic patterns: adjacent sub-arcs never agree (difference formula)
    rng = np.random.default_rng(3)
    worst = np.inf
    nchk = 0
    for _ in range(300):
        Q = int(rng.integers(3, 23))
        m = list(rng.integers(0, 2, size=Q))
        if all(m[r - 1] == m[Q - r - 1] for r in range(1, Q)):
            continue
        p = np.array([0.0] + [float(x) for x in m[:Q - 1]])     # p(0) = 0: the part off Q Z
        r = np.arange(1, Q)
        pa = (p[1:] - p[(Q - r) % Q]) / 2
        om = np.pi * (1 - 2 * r / Q)
        ug = np.linspace(0.02, 0.98, 97)
        for s in range(1, Q):
            A = pa * np.exp(-2j * np.pi * r * s / Q) * (1 - np.exp(2j * np.pi * r / Q))
            d = np.array([abs(np.sum(A * np.exp(1j * om * u)).imag / np.sin(np.pi * u)) for u in ug])
            worst = min(worst, d.max())
            nchk += 1
    lines.append(f"(c) non-symmetric periodic patterns: {nchk} adjacent pairs checked; "
                 f"min over pairs of sup_u |w'_(s-1) - w'_s| = {worst:.3f} (> 0, as Theorem 2.3 requires)")
    txt = "\n".join(lines)
    print(txt)
    with open("r3_periodic_flat_results.txt", "w") as f:
        f.write("r3_periodic_flat.py results (AI-produced, 2026-11-27)\n" + txt + "\n")


if __name__ == "__main__":
    main()
