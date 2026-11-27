"""Referee test 2 (independent): Theorem 1.2 of TOY.md on explicit hulls where w IS a measure.

Circle (A = T, G = Z). For symmetric periodic integer weights m(n) = a[n mod q] (n != 0),
a[r] = a[q-r], one has w = (1/2 - a_0) + sum_j beta_j 2 pi delta_{2 pi j/q},
beta_j = (1/q) sum_r a_r e^{2 pi i r j/q}, and on the arc (t_k, t_{k+1}), t_j = 2 pi j/q,

    F(t) = (1 - 2 a_0) + 4 pi M/(1 - e^{-2 pi M}) e^{-M (t - t_k)} p_k(x),
    p_k(x) = sum_{d=0}^{q-1} beta_{(k-d) mod q} x^d,   x = e^{-2 pi M/q}.

So F is flat on arc k iff p_k(x) = 0 with x in (0,1). We enumerate all symmetric patterns
with q <= 8 and weights in {0,1} and in {0,1,2,3}, find every flat arc (all M > 0), and
record the levels (Theorem 1.2 predicts: odd integers only). Several flat arcs are checked
against Fejer means of the Fourier series of F (independent of the closed form).
Then: finite perturbations (w still a measure, a.c. part a non-constant polynomial):
no flat arcs survive.

Torus (A = T^2, G = Z + alpha Z, alpha = log 2/log 3). Weights
    m = 3*1_{Z\\0} - 2*1_{2Z\\0} + 2*1_{alpha Z\\0} - 1_{2 alpha Z\\0}   (values 1,3 on Z, 1,2 on alpha Z)
give w a measure (Haar measures on closed subgroups of T^2) and, at M = log 2/pi,
F = -3 on the quarter (pi,2pi)^2 of T^2: a positive-measure level set at an ODD level on a
rank-2 hull, consistent with Theorem 1.2. Along random orbits the a.e. derivative of F is
M(P - F) with P = -3, as in step (ii). A random search over such coset weights looks for
positive-measure level sets at EVEN levels.
"""
import itertools
import math

import numpy as np

out = []


def log(s=""):
    print(s)
    out.append(s)


TWO_PI = 2 * math.pi


def beta_of(a):
    q = len(a)
    j = np.arange(q)
    r = np.arange(q)
    B = (np.exp(2j * np.pi * np.outer(j, r) / q) @ np.array(a, dtype=float)) / q
    assert np.max(np.abs(B.imag)) < 1e-12
    return B.real


def flat_arcs(a):
    """All (k, M) with F flat on arc k (M > 0)."""
    q = len(a)
    b = beta_of(a)
    res = []
    for k in range(q):
        coeffs = np.array([b[(k - d) % q] for d in range(q)])  # coefficient of x^d
        if np.max(np.abs(coeffs)) < 1e-14:
            res.append((k, None))  # identically flat (m = 0 case)
            continue
        coeffs = coeffs.copy()
        coeffs[np.abs(coeffs) < 1e-12] = 0.0
        # factor out x^k (roots x = 0 are not admissible: they would mean M = infinity)
        nz = np.nonzero(coeffs)[0]
        coeffs = coeffs[nz[0]:nz[-1] + 1]
        poly = coeffs[::-1].copy()  # highest degree first
        # p_k(1) = a_0: when a_0 = 0, x = 1 (i.e. M = 0, the degenerate limit) is a root;
        # deflate it (possibly repeatedly) so that only genuine M > 0 remain
        while len(poly) > 1 and abs(np.polyval(poly, 1.0)) < 1e-12:
            poly, rem = np.polydiv(poly, np.array([1.0, -1.0]))
        if len(poly) <= 1:
            continue
        roots = np.roots(poly)
        for x in roots:
            if abs(x.imag) < 1e-9 and 1e-6 < x.real < 1 - 1e-7:
                res.append((k, -q * math.log(x.real) / TWO_PI))
    return res


def F_closed(a, M, t):
    q = len(a)
    b = beta_of(a)
    t = np.mod(t, TWO_PI)
    tot = (1 - 2 * a[0]) * np.ones_like(t)
    for j in range(q):
        s = np.mod(t - TWO_PI * j / q, TWO_PI)
        tot += 4 * math.pi * M * b[j] * np.exp(-M * s) / (1 - math.exp(-TWO_PI * M))
    return tot


def F_fejer(mfun, M, t, N=40000):
    """Fejer mean of the Fourier series 1 + sum_{n != 0} m(|n|) 2M/(M+in) e^{int}."""
    n = np.arange(1, N + 1)
    mn = np.array([mfun(k) for k in n], dtype=float)
    fej = 1 - n / (N + 1)
    c = mn * fej * 2 * M / (M + 1j * n)
    vals = []
    for tt in np.atleast_1d(t):
        vals.append(1 + 2 * np.real(np.sum(c * np.exp(1j * n * tt))))
    return np.array(vals)


def symmetric_patterns(q, W):
    free = list(range(0, q // 2 + 1))  # residues 0..floor(q/2) determine the pattern
    for vals in itertools.product(W, repeat=len(free)):
        a = [0] * q
        for r, v in zip(free, vals):
            a[r] = v
            a[(q - r) % q] = v
        yield tuple(a)


log("Circle: flat arcs of F for symmetric periodic integer weights (all M > 0)")
for W in [(0, 1), (0, 1, 2, 3)]:
    levels = {}
    n_pat = 0
    n_flat = 0
    examples = []
    for q in range(1, 9):
        for a in symmetric_patterns(q, W):
            n_pat += 1
            for k, M in flat_arcs(a):
                n_flat += 1
                c = 1 - 2 * a[0]
                levels[c] = levels.get(c, 0) + 1
                if M is not None and len(examples) < 400:
                    examples.append((a, k, M, c))
    log(f"  weights in {W}: {n_pat} symmetric patterns (q <= 8), {n_flat} flat arcs; "
        f"levels (level: count) = {dict(sorted(levels.items()))}")
    # independent check of a few flat arcs against Fejer means of the Fourier series
    rng = np.random.default_rng(1)
    picks = [examples[i] for i in rng.choice(len(examples), size=min(4, len(examples)), replace=False)]
    for a, k, M, c in picks:
        q = len(a)
        tmid = TWO_PI * (k + 0.5) / q
        ts = tmid + np.array([-0.25, 0.0, 0.25]) * TWO_PI / q
        fe = F_fejer(lambda n: a[n % q], M, ts)
        cl = F_closed(a, M, ts)
        log(f"    pattern a={a}, arc {k}, M={M:.6f}: level {c}; closed form {np.round(cl, 6)}, "
            f"Fejer(N=4e4) {np.round(fe, 4)}")

# the example used in the report
a = (1, 2)
M = math.log(3) / math.pi
ts = np.linspace(math.pi + 0.2, TWO_PI - 0.2, 5)
log()
log(f"Example: m(n) = 1 (n even), 2 (n odd), M = log3/pi = {M:.6f}: F on (pi, 2pi):")
log(f"  closed form {np.round(F_closed(a, M, ts), 10)}")
log(f"  Fejer mean  {np.round(F_fejer(lambda n: a[n % 2], M, ts), 4)}")
ts2 = np.linspace(0.2, math.pi - 0.2, 5)
log(f"  and on (0, pi): closed form {np.round(F_closed(a, M, ts2), 4)} (not flat)")

# finite perturbations of the flat examples: w is still a measure, P is a non-constant polynomial
log()
log("Finite perturbations (w stays a measure; P non-constant): spread of F on the formerly flat arc")
rng = np.random.default_rng(7)
min_spread = np.inf
for trial in range(200):
    a, k, M, c = examples[rng.integers(len(examples))]
    q = len(a)
    K = rng.integers(1, 6)
    pert = {int(n): int(rng.integers(1, 3)) for n in rng.choice(np.arange(1, 13), size=K, replace=False)}
    tt = np.linspace(TWO_PI * k / q + 1e-3, TWO_PI * (k + 1) / q - 1e-3, 2001)
    base = F_closed(a, M, tt)
    add = np.zeros_like(tt)
    for n, b in pert.items():
        add += b * 2 * np.real(2 * M / (M + 1j * n) * np.exp(1j * n * tt))
    vals = base + add
    min_spread = min(min_spread, vals.max() - vals.min())
log(f"  200 random perturbations: minimal oscillation (max - min) of F on the arc = {min_spread:.3e} (> 0: no flat arc)")

# ------------------------------------------------------------- torus
log()
log("Torus T^2, G = Z + alpha Z, alpha = log2/log3:")
alpha = math.log(2) / math.log(3)
M = math.log(2) / math.pi
M2 = M / alpha  # effective parameter in the alpha-direction
log(f"  M = log2/pi = {M:.6f}; effective parameter in direction alpha: M/alpha = {M2:.6f} (= log3/pi = {math.log(3)/math.pi:.6f})")


def F_torus(th1, th2):
    return F_closed((1, 3), M, th1) + F_closed((1, 2), M2, th2) - 1.0


g1 = np.linspace(0, TWO_PI, 801, endpoint=False) + 1e-4
T1, T2 = np.meshgrid(g1, g1, indexing="ij")
Fv = F_torus(T1, T2)
frac_m3 = np.mean(np.abs(Fv + 3) < 1e-9)
log(f"  fraction of a 801^2 grid with |F + 3| < 1e-9: {frac_m3:.4f} (quarter of T^2 expected: 0.25)")
# independent check by Fejer means in each direction at a few points of the flat quarter
pts = [(4.0, 4.5), (3.5, 5.9), (5.0, 3.3)]
for p in pts:
    v1 = F_fejer(lambda n: (1, 3)[n % 2], M, [p[0]])[0]
    v2 = F_fejer(lambda n: (1, 2)[n % 2], M2, [p[1]])[0]
    log(f"    theta={p}: Fejer value of F = {v1 + v2 - 1:.5f}")
# a.c. part of 2w: P = 2(1/2 - sum c_i) with c = (3,-2,2,-1)
P = 2 * (0.5 - (3 - 2 + 2 - 1))
log(f"  P (a.c. density of 2w) = {P}; level of the flat set = {P}  (odd, as Theorem 1.2 requires)")
# along random orbits: derivative of F equals M(P - F) away from the jump points
rng = np.random.default_rng(3)
maxerr = 0.0
for _ in range(20):
    th = rng.uniform(0, TWO_PI, 2)
    t = np.linspace(0, 40, 400001)
    f = F_torus(th[0] + t, th[1] + alpha * t)
    df = np.gradient(f, t)
    ph1 = np.mod(th[0] + t, math.pi)
    ph2 = np.mod(th[1] + alpha * t, math.pi)
    good = (np.minimum(ph1, math.pi - ph1) > 1e-3) & (np.minimum(ph2, math.pi - ph2) > 1e-3)
    err = np.max(np.abs(df[good] - M * (P - f[good])))
    maxerr = max(maxerr, err)
log(f"  along 20 random orbits (t in [0,40]): max |dF/dt - M(P - F)| off the jump points = {maxerr:.2e}")

# random search for positive-measure level sets at EVEN levels.
# F = F^(1)(theta_1) + F^(2)(theta_2) - 1 is flat exactly where both direction classes are flat
# at the same M; so in half of the draws alpha is CHOSEN as M_1/M_2 for two flat conditions
# (M_1 of class 1 at speed 1, M_2 of class 2 at speed alpha), giving genuine flat sets on T^2.
log()
log("Random search on T^2 = hull of G = Z + alpha Z (two direction classes, symmetric periodic weights in {0,1,2}):")
log("  type A draws: random M, alpha = log2/log3; type B draws: alpha = M_1/M_2 chosen so that both classes are flat")
rng = np.random.default_rng(11)
taus = [1e-2, 1e-4, 1e-6, 1e-8]
best = {(typ, par, tau): 0.0 for typ in "AB" for par in ("even", "odd") for tau in taus}
levels_B = {}
g = np.linspace(0, TWO_PI, 512, endpoint=False) + 1e-5
T1, T2 = np.meshgrid(g, g, indexing="ij")


def random_pattern():
    while True:
        q = int(rng.integers(1, 7))
        pat = tuple(next(itertools.islice(symmetric_patterns(q, (0, 1, 2)),
                                          int(rng.integers(0, 3 ** (q // 2 + 1))), None)))
        if any(pat):  # exclude the trivial pattern m = 0 in a direction
            return pat


n_draws = {"A": 0, "B": 0}
while n_draws["A"] < 200 or n_draws["B"] < 200:
    p1, p2 = random_pattern(), random_pattern()
    typ = "A" if n_draws["A"] < 200 else "B"
    if typ == "A":
        Mx, al = float(rng.uniform(0.05, 2.0)), alpha
    else:
        f1 = [x[1] for x in flat_arcs(p1) if x[1] is not None]
        f2 = [x[1] for x in flat_arcs(p2) if x[1] is not None]
        if not f1 or not f2:
            continue
        Mx = f1[int(rng.integers(len(f1)))]
        M2x = f2[int(rng.integers(len(f2)))]
        al = Mx / M2x
        if any(abs(al * qd - round(al * qd)) < 1e-9 for qd in range(1, 51)):
            continue  # alpha (numerically) rational: hull would not be T^2
    n_draws[typ] += 1
    Fv = F_closed(p1, Mx, T1) + F_closed(p2, Mx / al, T2) - 1.0
    de = np.abs(Fv / 2 - np.round(Fv / 2)) * 2
    do = np.abs((Fv - 1) / 2 - np.round((Fv - 1) / 2)) * 2
    for tau in taus:
        best[(typ, "even", tau)] = max(best[(typ, "even", tau)], float(np.mean(de < tau)))
        best[(typ, "odd", tau)] = max(best[(typ, "odd", tau)], float(np.mean(do < tau)))
    if typ == "B":
        flatvals = Fv[np.abs(Fv - np.round(Fv)) < 1e-9]
        for v in np.unique(np.round(flatvals)):
            levels_B[int(v)] = levels_B.get(int(v), 0) + 1
for typ in "AB":
    for tau in taus:
        log(f"  type {typ}, tau={tau:g}: max over 200 draws of the fraction of a 512^2 grid within tau of an EVEN integer = "
            f"{best[(typ, 'even', tau)]:.5f}; of an ODD integer = {best[(typ, 'odd', tau)]:.5f}")
log(f"  type B: integer levels of the flat sets found (level: number of draws) = {dict(sorted(levels_B.items()))}")

with open(__file__.replace(".py", "_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
