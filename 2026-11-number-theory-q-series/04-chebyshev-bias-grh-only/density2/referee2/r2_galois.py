"""Referee 2: check of TOY.md Theorem 2.3 step (iv) (Galois proportions), q = 2..30.

For patterns v : Z/q -> {0,1}, v(0) = 0, V(e) := sum_r v(r) zeta_q^{re}.
Values are represented EXACTLY in Z[zeta_q] = Z[x]/(Phi_q) (reduced integer vectors).
Claims tested:
  (a) V(a e) = sigma_a(V(e)) exactly, a in (Z/q)^x;
  (b) on each orbit O of (Z/q)^x acting on {1..q-1}, for each value beta taken on O,
      #{e in O : V(e) = beta}/|O| = 1/[Q(beta):Q].
      The degree is computed WITHOUT the Galois action: as the number of distinct
      eigenvalues of the multiplication-by-beta matrix on Q(zeta_q) (numerical, tol 1e-7),
      cross-checked on a sample by sympy's exact minimal polynomial.
  (c) consequence used in the proof: if some value is taken on >= 60% of {1..q-1},
      it is a rational integer.  We list which integers occur this way (all patterns,
      q <= 13; random patterns beyond).
"""
import itertools
import math
import numpy as np
import sympy as sp

rng = np.random.default_rng(11)
X = sp.symbols("x")


def cyclo(q):
    return [int(c) for c in sp.Poly(sp.cyclotomic_poly(q, X), X).all_coeffs()]  # highest first


def reduce_vec(a, phi_c):
    """a: integer coefficients a[k] of x^k (k = 0..len-1); reduce mod monic Phi (highest-first coeffs)."""
    a = list(map(int, a))
    d = len(phi_c) - 1
    for k in range(len(a) - 1, d - 1, -1):
        c = a[k]
        if c:
            # x^k = x^{k-d} * x^d, x^d = -(phi_c[1] x^{d-1} + ... + phi_c[d])
            for j in range(1, d + 1):
                a[k - j] -= c * phi_c[j]
            a[k] = 0
    out = a[:d] + [0] * max(0, d - len(a))
    return tuple(out[:d])


def V_exact(v, e, q, phi_c):
    a = [0] * q
    for r in range(q):
        if v[r]:
            a[(r * e) % q] += 1
    return reduce_vec(a, phi_c)


def sigma(vec, a_, q, phi_c):
    b = [0] * q
    for k, c in enumerate(vec):
        if c:
            b[(k * a_) % q] += c
    return reduce_vec(b, phi_c)


def mult_matrix(vec, phi_c):
    d = len(phi_c) - 1
    Mx = np.zeros((d, d))
    for j in range(d):
        # beta * x^j
        a = [0] * (2 * d)
        for k, c in enumerate(vec):
            a[k + j] += c
        Mx[:, j] = reduce_vec(a, phi_c)
    return Mx


def degree_numeric(vec, phi_c):
    ev = np.linalg.eigvals(mult_matrix(vec, phi_c))
    distinct = []
    for z in ev:
        if all(abs(z - w) > 1e-7 for w in distinct):
            distinct.append(z)
    return len(distinct)


def degree_exact(vec, q):
    z = sp.exp(2 * sp.pi * sp.I / q)
    beta = sum(c * z ** k for k, c in enumerate(vec))
    return sp.degree(sp.minimal_polynomial(beta, X), X)


def main():
    out = ["Referee 2, r2_galois.py: Galois proportions in TOY.md Thm 2.3 step (iv)."]
    total_checks = mism_b = mism_a = 0
    exact_checked = exact_mism = 0
    integer_majorities = {}
    for q in range(2, 31):
        phi_c = cyclo(q)
        units = [a for a in range(1, q) if math.gcd(a, q) == 1]
        orbits = []
        seen = set()
        for e0 in range(1, q):
            if e0 not in seen:
                O = sorted(set((a * e0) % q for a in units))
                orbits.append(O)
                seen.update(O)
        if q <= 13:
            pats = [np.array((0,) + b) for b in itertools.product([0, 1], repeat=q - 1)]
        else:
            pats = []
            for _ in range(150):
                v = rng.integers(0, 2, q)
                v[0] = 0
                pats.append(v)
                vs = v.copy()
                for r in range(1, q):
                    vs[r] = v[min(r, q - r)]
                pats.append(vs)
        deg_cache = {}
        for v in pats:
            vals = {e: V_exact(v, e, q, phi_c) for e in range(1, q)}
            # (a)
            for e in range(1, q):
                for a_ in units:
                    total_checks += 1
                    if vals[(a_ * e) % q] != sigma(vals[e], a_, q, phi_c):
                        mism_a += 1
            # (b)
            for O in orbits:
                counts = {}
                for e in O:
                    counts[vals[e]] = counts.get(vals[e], 0) + 1
                for beta, cnt in counts.items():
                    if beta not in deg_cache:
                        deg_cache[beta] = degree_numeric(beta, phi_c)
                    if cnt * deg_cache[beta] != len(O):
                        mism_b += 1
            # (c)
            counts = {}
            for e in range(1, q):
                counts[vals[e]] = counts.get(vals[e], 0) + 1
            for beta, cnt in counts.items():
                if cnt >= 0.6 * (q - 1):
                    is_int = all(c == 0 for c in beta[1:])
                    key = beta[0] if is_int else ("NON-INTEGER", beta)
                    integer_majorities.setdefault(q, set()).add(key)
        # exact cross-check of the numerical degree on a sample of values
        sample = list(deg_cache.items())
        rng.shuffle(sample)
        for beta, dnum in sample[:4 if q > 20 else 8]:
            exact_checked += 1
            if degree_exact(beta, q) != dnum:
                exact_mism += 1
        out.append(f"q={q:2d}: {len(pats):5d} patterns, {len(orbits)} orbits, distinct values {len(deg_cache):5d}; "
                   f"values taken on >=60% of {{1..q-1}}: {sorted(integer_majorities.get(q, []), key=str)}")
    out.append(f"\n(a) V(ae) = sigma_a V(e): {total_checks} checks, {mism_a} mismatches")
    out.append(f"(b) proportion on orbit = 1/deg: {mism_b} mismatches (all orbits, all values, all patterns above)")
    out.append(f"    numerical degree vs sympy exact minimal polynomial: {exact_checked} checked, {exact_mism} mismatches")
    nonint = [q for q, s in integer_majorities.items() if any(isinstance(k, tuple) for k in s)]
    out.append(f"(c) q with a NON-integer value taken on >= 60%: {nonint if nonint else 'none'}")
    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace(".py", "_results.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
