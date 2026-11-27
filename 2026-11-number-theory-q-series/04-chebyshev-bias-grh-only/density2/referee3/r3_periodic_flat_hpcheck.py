"""Referee 3: validation of the root classifier of r3_periodic_flat.py in 60-digit arithmetic.

AI-produced, 2026-11-27.

For a random sample of (pattern, sub-arc) pairs with 12 <= Q <= 22, and for every pair on which a
naive double-precision classifier disagrees with flat_roots() of r3_periodic_flat.py, the distinct roots of
f_s in (0,1) are recomputed entirely in 60-digit arithmetic (exact zeros of W detected at that
precision, the exact root x = 1 removed when m_Q = 0, multiple roots counted once) and compared
with the output of flat_roots() of r3_periodic_flat.py.
Output: r3_periodic_flat_hpcheck_results.txt
"""
import random
import sys
from multiprocessing import Pool
import numpy as np
import mpmath as mp

sys.path.insert(0, ".")
import r3_periodic_flat as R  # noqa: E402


def hp_count(args):
    m, s = args
    Q = len(m)
    with mp.workdps(60):
        p = [mp.mpf(m[Q - 1])] + [mp.mpf(x) for x in m[:Q - 1]]
        W = [mp.fsum(p[r] * mp.cos(2 * mp.pi * r * ss / Q) for r in range(Q)) for ss in range(Q)]
        coef = [W[(s - j) % Q] for j in range(Q)]
        small = mp.mpf(10) ** (-45)
        coef = [c if abs(c) > small else mp.mpf(0) for c in coef]
        while len(coef) > 1 and coef[-1] == 0:
            coef = coef[:-1]
        while len(coef) > 1 and coef[0] == 0:
            coef = coef[1:]
        if len(coef) <= 1:
            return (m, s, [])
        coef, _ = R.deflate_at_one(coef, small)
        if len(coef) <= 1:
            return (m, s, [])
        rts = mp.polyroots(coef[::-1], maxsteps=3000, extraprec=600)
        out = []
        for z in rts:
            if abs(mp.im(z)) < mp.mpf(10) ** (-25) and 0 < mp.re(z) < 1:
                x = float(mp.re(z))
                if not any(abs(x - y) < 1e-12 for y in out):
                    out.append(x)
    return (m, s, sorted(out))


def naive(m, s):
    Q = len(m)
    p = np.array([m[Q - 1]] + list(m[:Q - 1]), dtype=float)
    r = np.arange(Q)
    W = np.array([np.sum(p * np.cos(2 * np.pi * r * ss / Q)) for ss in range(Q)])
    c = np.array([W[(s - j) % Q] for j in range(Q)])
    while len(c) > 1 and abs(c[-1]) < 1e-12:
        c = c[:-1]
    if len(c) <= 1:
        return []
    out = []
    for z in np.roots(c[::-1]):
        if abs(z.imag) < 1e-7 and 1e-9 < z.real < 1 - 1e-6:
            if not any(abs(z.real - y) < 1e-6 for y in out):
                out.append(z.real)
    return out


def main():
    rng = random.Random(11)
    pool_pats = {Q: [m for m in R.symmetric_patterns(Q) if R.primitive_period(m) == Q]
                 for Q in range(12, 23)}
    sample = []
    weights = [len(pool_pats[Q]) * Q for Q in range(12, 23)]
    Qs = list(range(12, 23))
    for _ in range(1600):
        Q = rng.choices(Qs, weights=weights)[0]
        m = rng.choice(pool_pats[Q])
        sample.append((m, rng.randrange(Q)))
    # pairs on which a naive double-precision classifier (plain numpy roots, |Im| < 1e-7,
    # 1e-9 < x < 1 - 1e-6, roots merged within 1e-6) disagrees with flat_roots()
    flagged = []
    naive_total = 0
    for Q in range(2, 23):
        for m in R.symmetric_patterns(Q):
            if R.primitive_period(m) != Q:
                continue
            cur = R.flat_roots(m)[0]
            for s in range(Q):
                a = naive(m, s)
                naive_total += len(a)
                if len(a) != len([x for (ss, x) in cur if ss == s]):
                    flagged.append((m, s))
    jobs = sample + flagged
    with Pool(4) as pool:
        res = pool.map(hp_count, jobs, chunksize=8)
    nd_s = nd_f = 0
    ex = []
    for k, (m, s, hp) in enumerate(res):
        dp = [x for (ss, x) in R.flat_roots(m)[0] if ss == s]
        same = len(dp) == len(hp) and all(abs(a - b) < 1e-7 for a, b in zip(dp, hp))
        if not same:
            if k < len(sample):
                nd_s += 1
            else:
                nd_f += 1
            if len(ex) < 5:
                ex.append((''.join(map(str, m)), s, dp, hp))
    lines = [f"random sample: {len(sample)} (pattern, sub-arc) pairs, 12 <= Q <= 22: "
             f"{nd_s} disagreements between flat_roots() and the 60-digit solve",
             f"pairs where the naive classifier disagrees with flat_roots(): {len(flagged)} "
             f"(naive total over Q = 2..22: {naive_total} triples); 60-digit solve vs flat_roots() "
             f"on these pairs: {nd_f} disagreements"]
    for e in ex:
        lines.append(f"   disagreement: m={e[0]} s={e[1]} double={e[2]} 60-digit={e[3]}")
    txt = "\n".join(lines)
    print(txt)
    with open("r3_periodic_flat_hpcheck_results.txt", "w") as f:
        f.write("r3_periodic_flat_hpcheck.py results (AI-produced, 2026-11-27)\n" + txt + "\n")


if __name__ == "__main__":
    main()
