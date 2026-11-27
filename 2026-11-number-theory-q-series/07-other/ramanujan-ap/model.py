#!/usr/bin/env python3
"""Heuristic model (GRH + LI) for the failure set of R_{q,a}.  NOT part of any theorem.

Exact reformulation (no hypothesis).  With L = log x, phi = phi(q),
  eps(x) = (pi(x;q,a) - li(x)/phi) L / sqrt(x),  lam = L li(x)/x,  mu = L/(L-1),
  M(x) = (e x/L) li(x/e) - li(x)^2,  W(x) = M(x) L^2 / x^{3/2}:
  R_{q,a}(x) fails  <=>  T(x) := phi [2 lam eps(x) - sqrt(e) mu eps(x/e)] + phi^2 eps(x)^2/sqrt(x) >= W(x).
W does not depend on q.  Explicit formula (GRH), up to O(1/L):
  phi eps(x) ~ -c(q,a) beta(L) - sum_chi conj(chi)(a) sum_gamma x^{i gamma}/(1/2 + i gamma),
  c(q,a) = #{b mod q : b^2 = a},  beta(L) = L li(sqrt x)/(2 sqrt x) -> 1,
so   T(x) ~ -c(q,a) B(L) + Y_a(log x),   B(L) = 2 lam beta(L) - sqrt(e) mu beta(L-1) -> 2 - sqrt(e) = 0.3513,
     Y_a(u) = - sum_chi conj(chi)(a) sum_gamma (2 - sqrt(e) e^{-i gamma}) e^{i gamma u}/(1/2 + i gamma).
Under LI, Y_a(u) has the limiting distribution of  sum_j r_j cos(theta_j), theta_j iid uniform,
r_j = 2 |2 - sqrt(e) e^{-i gamma_j}| / |1/2 + i gamma_j| over the zero pairs of all L(s,chi), chi mod q
(the same for every a).  We print variances, Chernoff/saddle-point tail probabilities
P(Y >= w), the model failure probability P_fail(x) = P(Y >= W(x) + c B(L)) and the expected
logarithmic measure of {x > X : R fails},  m(X) = int_{log X}^infty P_fail(e^u) du.
"""
import sys, math, os
import numpy as np
from mpmath import mp, mpf, li as mpli, log as mplog, e as mpe, sqrt as mpsqrt
mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "prototype", "explore")
TCUT = 400.0
def load(fn):
    return np.array([float(l.split()[0]) for l in open(fn) if l.strip() and not l.startswith("#")])
import hashlib
PIN = {"chi4_zeros.txt": "0e73e00a6c1609bdb4eec62c13292435d0beee55750735a51f0fab84024ddfbc",
       "chi3_zeros.txt": "e3075235a847f43eb35057d3998a853f9383e27e958cc5fd3c4115dad7eb8aaa"}
for _fn, _h in PIN.items():       # the certified zero lists are pinned (see nonmonotone-audit/AUDIT.md, G1)
    assert hashlib.sha256(open(os.path.join(ROOT, _fn), "rb").read()).hexdigest() == _h, _fn
Z = {}
Z["zeta"] = load(os.path.join(HERE, "zeros_zeta.txt"))
Z["chi4"] = load(os.path.join(ROOT, "chi4_zeros.txt"))
Z["chi3"] = load(os.path.join(ROOT, "chi3_zeros.txt"))
def load_q(q):
    d = {}
    fn = os.path.join(HERE, "zeros_q%d.txt" % q)
    if not os.path.exists(fn): return d
    for l in open(fn):
        if l.startswith("#"): continue
        qq, idx, real, g = l.split()
        d.setdefault((int(idx), int(real)), []).append(float(g))
    return d
# characters mod q: list of (name, conductor, zero ordinates as "pairs")
def families(q):
    fam = [("zeta", 1, Z["zeta"])]
    if q == 3: fam.append(("chi3", 3, Z["chi3"]))
    if q == 4: fam.append(("chi4", 4, Z["chi4"]))
    if q in (5, 8, 12):
        d = load_q(q)
        if q == 8: fam.append(("chi4", 4, Z["chi4"]))
        if q == 12: fam += [("chi3", 3, Z["chi3"]), ("chi4", 4, Z["chi4"])]
        for (idx, real), zs in sorted(d.items()):
            zs = np.array(zs)
            if real: fam.append(("chi%d_%d" % (q, idx), q, zs[zs > 0]))
            elif idx == 2: fam.append(("chi%d_%d+conj" % (q, idx), q, zs))   # one entry per complex pair
    return fam
SQE = math.sqrt(math.e)
def amps(q):
    r = []; tailvar = 0.0; info = []
    for name, cond, zs in families(q):
        zs = zs[np.abs(zs) <= TCUT]
        rr = 2 * np.abs(2 - SQE * np.exp(-1j * zs)) / np.abs(0.5 + 1j * zs)
        r.append(rr)
        # tail: pairs with gamma > TCUT, density (1/2pi) log(cond t/2pi) per unit t (per character;
        # a complex pair has twice that), mean of |2 - sqrt(e) e^{-i g}|^2 = 4 + e
        mult = 2 if name.endswith("+conj") else 1
        tv = mult * 0.5 * 4 * (4 + math.e) * (1 / (2 * math.pi)) * (math.log(cond * TCUT / (2 * math.pi)) + 1) / TCUT
        tailvar += tv
        info.append((name, len(zs), float(np.sum(rr ** 2) / 2), tv))
    return np.concatenate(r), tailvar, info
def W_of(x):
    x = mpf(x); L = mplog(x)
    return float(((mpe * x / L) * mpli(x / mpe) - mpli(x) ** 2) * L ** 2 / x ** mpf(1.5))
def B_of(x):
    x = mpf(x); L = mplog(x)
    lam = L * mpli(x) / x; mu = L / (L - 1)
    beta = lambda y: mplog(y) * mpli(mpsqrt(y)) / (2 * mpsqrt(y))
    return float(2 * lam * beta(x) - mpsqrt(mpe) * mu * beta(x / mpe))
def tail_prob(r, tv, w):
    """saddle-point (Lugannani-Rice type leading term) and Chernoff bound for P(Y >= w)"""
    from scipy.special import ive
    def K(s):   # log prod I0(r s) + tv s^2/2
        z = r * s
        return float(np.sum(np.log(ive(0, z)) + z) + tv * s * s / 2)
    def K1(s):
        z = r * s
        return float(np.sum(r * ive(1, z) / ive(0, z)) + tv * s)
    def K2(s):
        z = r * s; q1 = ive(1, z) / ive(0, z)
        return float(np.sum(r * r * (1 - q1 / np.where(z > 0, z, 1) - q1 * q1)) + tv)
    lo, hi = 0.0, 1.0
    while K1(hi) < w: hi *= 2
    for _ in range(100):
        mid = (lo + hi) / 2
        if K1(mid) < w: lo = mid
        else: hi = mid
    s = hi
    logchern = K(s) - s * w
    sp = logchern - math.log(s * math.sqrt(2 * math.pi * K2(s)))
    return logchern / math.log(10), sp / math.log(10)     # log10 of bound, log10 of saddle-point estimate
CSQ = {(1, 0): 1, (3, 1): 2, (3, 2): 0, (4, 1): 2, (4, 3): 0, (5, 1): 2, (5, 4): 2, (5, 2): 0, (5, 3): 0,
       (8, 1): 4, (8, 3): 0, (8, 5): 0, (8, 7): 0, (12, 1): 4, (12, 5): 0, (12, 7): 0, (12, 11): 0}
if __name__ == "__main__":
    print("x          W(x)       B(L)")
    for k in [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 15, 16, 17, 18, 19]:
        x = 10 ** k
        print("1e%-5s  %10.4f  %8.5f" % (k, W_of(x), B_of(x)))
    for q in (1, 3, 4, 5, 8, 12):
        try: r, tv, info = amps(q)
        except Exception as ex:
            print("q =", q, "zeros missing:", ex); continue
        var = float(np.sum(r ** 2) / 2 + tv)
        print("\nq = %d: families %s" % (q, [(n, k, round(v, 4), round(t, 4)) for n, k, v, t in info]))
        print("   variance of Y = %.4f (sigma = %.4f), tail part %.4f" % (var, math.sqrt(var), tv))
        classes = sorted(set(c for (qq, a), c in CSQ.items() if qq == q))
        for c in classes:
            print("   c(q,a) = %d:   x: log10 P_fail (saddle point) [Chernoff bound]" % c)
            row = []
            us = np.arange(math.log(1e10), math.log(1e20), 0.05)
            lp = []
            for u in us:
                x = math.exp(u); w = W_of(x) + c * B_of(x)
                ch, sp = tail_prob(r, tv, w)
                lp.append(min(sp, 0.0))
            lp = np.array(lp)
            for k in [10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 15, 16, 17, 18]:
                i = int(round((math.log(10 ** k) - us[0]) / 0.05))
                ch, sp = tail_prob(r, tv, W_of(10 ** k) + c * B_of(10 ** k))
                row.append("1e%s: %.2f [%.2f]" % (k, sp, ch))
            print("      " + "; ".join(row))
            # expected log-measure beyond X
            dens = 10 ** lp
            cum = np.cumsum(dens[::-1])[::-1] * 0.05
            outs = []
            for k in [11, 12, 13, 14]:
                i = int(round((math.log(10 ** k) - us[0]) / 0.05))
                outs.append("m(1e%d) = %.2e" % (k, cum[i]))
            print("      expected log-measure of failures beyond X: " + ", ".join(outs))
