"""Referee 3: exact check of the Galois step (iv) of TOY.md Theorem 2.3 and of Corollary 2.4,
Steps 1(b) and 2, plus the numbers quoted in the Remarks on Corollary 2.4 (3).

AI-produced, 2026-11-27.

Exact arithmetic: an element of Z[zeta_q] is stored as its integer coordinate vector in the
power basis 1, zeta, ..., zeta^{phi(q)-1} (reduction of x^k modulo the cyclotomic polynomial,
which is monic with integer coefficients). V(e) = sum_r v(r) zeta^{re}. Two values are equal iff
their vectors are equal; a value is rational iff all non-constant coordinates vanish.

Checks, for every pattern v : Z/q -> W with v(0) = 0 in the stated families:
 (iv-1) a value taken on >= 60% of e in {1..q-1} is a rational integer; also the largest
        proportion of {1..q-1} on which a NON-rational value is taken (the proof needs < 0.6;
        the orbit argument gives <= 1/2);
 (iv-2) on each orbit O of (Z/q)^x acting on {1..q-1}, every value taken on O is taken equally
        often (the fibre-counting claim of step (iv));
 (C1b)  q prime: if one V(e_0), e_0 != 0, is rational then V is constant on {1..q-1};
 (C2)   if V is constant = beta on {1..q-1} then v = -beta 1_{!=0} (Fourier inversion), so for
        0/1 patterns beta in {0,-1};
 (Rk3)  the Ramanujan-sum values for v = 1_{units} (q = 21 and q = p1 p2), and the smallest q
        for which some symmetric 0/1 pattern has an integer value n not in {0,-1} on >= 60% of
        {1..q-1} (the single-level arithmetic then allows a flat level kappa = 2n+1 with |kappa|>=3).
Output: r3_galois_results.txt
"""
import itertools
from math import gcd
import numpy as np
import sympy as sp

X = sp.symbols("x")


def reduction_table(q):
    phi = [int(c) for c in sp.Poly(sp.cyclotomic_poly(q, X), X).all_coeffs()][::-1]  # low->high
    d = len(phi) - 1
    tab = np.zeros((q, d), dtype=np.int64)
    cur = np.zeros(d + 1, dtype=np.int64)
    cur[0] = 1
    for k in range(q):
        tab[k] = cur[:d]
        nxt = np.zeros(d + 1, dtype=np.int64)
        nxt[1:] = cur[:d]
        if nxt[d] != 0:
            nxt = nxt - nxt[d] * np.array(phi, dtype=np.int64)
        cur = nxt
    return tab  # tab[k] = coordinates of zeta^k


def V_vectors(patterns, q, tab):
    """patterns: (n_pat, q) int array. Returns (n_pat, q, d) array of V(e), e = 0..q-1."""
    r = np.arange(q)
    out = np.zeros((patterns.shape[0], q, tab.shape[1]), dtype=np.int64)
    for e in range(q):
        out[:, e, :] = patterns @ tab[(r * e) % q]
    return out


def orbits(q):
    units = [a for a in range(1, q) if gcd(a, q) == 1]
    seen, orbs = set(), []
    for e in range(1, q):
        if e in seen:
            continue
        O = sorted({(a * e) % q for a in units})
        seen.update(O)
        orbs.append(O)
    return orbs


def is_prime(q):
    return sp.isprime(q)


def analyse(q, patterns, label, lines, stats):
    tab = reduction_table(q)
    Vv = V_vectors(patterns, q, tab)
    orbs = orbits(q)
    n_bad_iv1 = n_bad_iv2 = n_bad_c1b = n_bad_c2 = 0
    max_nonrat = 0.0
    big_int = []
    for ip in range(patterns.shape[0]):
        vals = [tuple(Vv[ip, e]) for e in range(1, q)]
        counts = {}
        for val in vals:
            counts[val] = counts.get(val, 0) + 1
        for val, c in counts.items():
            rational = all(x == 0 for x in val[1:])
            frac = c / (q - 1)
            if not rational:
                max_nonrat = max(max_nonrat, frac)
            if frac >= 0.6 and not rational:
                n_bad_iv1 += 1
            if frac >= 0.6 and rational and val[0] not in (0, -1):
                big_int.append((tuple(int(x) for x in patterns[ip]), int(val[0]), c))
        # (iv-2) equidistribution on orbits
        for O in orbs:
            oc = {}
            for e in O:
                oc[tuple(Vv[ip, e])] = oc.get(tuple(Vv[ip, e]), 0) + 1
            if len(set(oc.values())) != 1:
                n_bad_iv2 += 1
        # (C1b)
        if is_prime(q):
            rat_somewhere = any(all(x == 0 for x in val[1:]) for val in vals)
            if rat_somewhere and len(set(vals)) != 1:
                n_bad_c1b += 1
        # (C2)
        if len(set(vals)) == 1:
            beta = vals[0]
            if not all(x == 0 for x in beta[1:]):
                n_bad_c2 += 1
            else:
                b = beta[0]
                if not all(patterns[ip, r] == -b for r in range(1, q)):
                    n_bad_c2 += 1
    stats["patterns"] += patterns.shape[0]
    stats["viol"] += n_bad_iv1 + n_bad_iv2 + n_bad_c1b + n_bad_c2
    stats["max_nonrat"] = max(stats["max_nonrat"], max_nonrat)
    lines.append(f"{label} q={q}: {patterns.shape[0]} patterns; violations (iv-1) {n_bad_iv1}, "
                 f"(iv-2) {n_bad_iv2}, (C1b) {n_bad_c1b}, (C2) {n_bad_c2}; "
                 f"max share of a non-rational value {max_nonrat:.3f}; "
                 f"integer values n not in {{0,-1}} on >= 60%: {len(big_int)}"
                 + (f" (e.g. v={big_int[0][0]}, n={big_int[0][1]}, on {big_int[0][2]} of {q-1})"
                    if big_int else ""))
    return big_int


def symmetric_patterns(q, W=(0, 1)):
    free = list(range(1, q // 2 + 1))
    pats = []
    for vals in itertools.product(W, repeat=len(free)):
        v = [0] * q
        for r, x in zip(free, vals):
            v[r] = x
            v[q - r] = x
        pats.append(v)
    return np.array(pats, dtype=np.int64)


def all_patterns(q, W=(0, 1)):
    pats = []
    for vals in itertools.product(W, repeat=q - 1):
        pats.append([0] + list(vals))
    return np.array(pats, dtype=np.int64)


def main():
    lines = []
    stats = {"patterns": 0, "viol": 0, "max_nonrat": 0.0}
    first_big = None
    for q in range(2, 31):
        big = analyse(q, symmetric_patterns(q), "symmetric 0/1", lines, stats)
        if big and first_big is None:
            first_big = (q, big[0], len(big))
    for q in range(2, 15):
        analyse(q, all_patterns(q), "all 0/1 (non-symmetric included)", lines, stats)
    for q in range(2, 10):
        analyse(q, all_patterns(q, W=(0, 1, 2)), "all {0,1,2}-valued", lines, stats)
    lines.append("")
    lines.append(f"TOTAL: {stats['patterns']} patterns checked; total violations {stats['viol']}; "
                 f"largest share of {{1..q-1}} on which a non-rational value is taken: "
                 f"{stats['max_nonrat']:.3f} (proof needs < 0.6)")
    if first_big:
        q, (v, n, c), cnt = first_big
        lines.append(f"Smallest q with a symmetric 0/1 pattern having an integer value n not in "
                     f"{{0,-1}} on >= 60% of {{1..q-1}}: q={q} ({cnt} patterns; e.g. v={v}, "
                     f"n={n} on {c} of {q-1}, i.e. kappa = 2n+1 = {2*n+1})")
    # Remark (3): Ramanujan sums
    for q in (21, 15, 35, 77, 221):
        tab = reduction_table(q)
        v = np.array([[1 if (r > 0 and gcd(r, q) == 1) else 0 for r in range(q)]], dtype=np.int64)
        Vv = V_vectors(v, q, tab)[0]
        vals = [int(Vv[e][0]) if all(x == 0 for x in Vv[e][1:]) else None for e in range(1, q)]
        cnt1 = sum(1 for x in vals if x == 1)
        lines.append(f"Remark (3): v = 1_units mod {q}: V(e) values {sorted(set(vals))}; "
                     f"V = mu(q) = {sp.mobius(q)} on the units; V = 1 on {cnt1} of {q-1} "
                     f"(share {cnt1/(q-1):.3f})")
    txt = "\n".join(lines)
    print(txt)
    with open("r3_galois_results.txt", "w") as f:
        f.write("r3_galois.py results (AI-produced, 2026-11-27)\n" + txt + "\n")


if __name__ == "__main__":
    main()
