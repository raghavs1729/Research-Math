#!/usr/bin/env python3
"""
sanity_primes.py -- [num] sanity checks on the prime side (no rigor), q = 4.
 (1) reflection inequality of DENSITY.md Thm 5.1(b): empirical log-distribution of
     E^theta_{1/2}(x) = (theta(x;4,3)-theta(x;4,1))/sqrt(x) on [1e4, 1e8]:
     P(E <= t) <= P(E >= 2 - sigma^2 - t) for a grid of t, and near-symmetry about 1;
 (2) two-sided transversality of Thm 2.3 / (2.1)-(2.2):
     E^theta_M(x)/M - E^theta_M'(x)/M' in [kappa*lam, kappa*(2/(M M') - lam)] up to o(1).
Output: sanity_primes_results.txt.  Runtime ~40 s, ~1.5 GB.
"""
import os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
Z = [float(l) for l in open(os.path.join(HERE, "..", "..", "prototype", "explore", "chi4_zeros.txt")) if l.strip()]
X = 10 ** 8
out = open(os.path.join(HERE, "sanity_primes_results.txt"), "w")
def log(s=""):
    print(s); out.write(s + "\n")

sieve = np.ones(X // 2, dtype=bool)          # odd numbers 1,3,5,...
sieve[0] = False
for i in range(1, int(X ** 0.5) // 2 + 1):
    if sieve[i]:
        p = 2 * i + 1
        sieve[p * p // 2::p] = False
primes = 2 * np.nonzero(sieve)[0].astype(np.int64) + 1
del sieve
chi_star = np.where(primes % 4 == 3, 1.0, -1.0)
lp = np.log(primes.astype(np.float64))

ys = np.linspace(math.log(1e4), math.log(X), 2_000_001)
xs = np.exp(ys)
idx = np.searchsorted(primes, xs, side="right")

def Etheta(M):
    m = M - 0.5
    c = np.concatenate(([0.0], np.cumsum(chi_star * primes.astype(np.float64) ** m * lp)))
    return 2 * M * xs ** (-M) * c[idx]

s2 = 0.15556797992358593
E = Etheta(0.5)
log("q = 4, theta-form, x in [1e4, 1e8], %d log-spaced samples" % len(ys))
log("(1) mean %.4f  variance %.4f  (model: 1, sigma^2 = %.4f)  median %.4f  min %.4f  max %.4f"
    % (E.mean(), E.var(), s2, np.median(E), E.min(), E.max()))
log("    empirical P(E <= 0) = %.5f" % (E <= 0).mean())
ok = True
for t in (-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 0.9):
    L = (E <= t).mean(); Rr = (E >= 2 - s2 - t).mean(); Rs = (E >= 2 - t).mean()
    log("    t = %5.2f   P(E<=t) = %.5f   P(E >= 2-sigma^2-t) = %.5f   [P(E >= 2-t) = %.5f]" % (t, L, Rr, Rs))
    ok = ok and (L <= Rr + 0.01)
log("    reflection inequality (finite-range, tolerance 0.01): %s" % ("consistent" if ok else "VIOLATED"))

def lam(M, Mp):
    s = sum(4 / (math.hypot(M, g) * math.hypot(Mp, g)) for g in Z)
    T = Z[-1]                                  # crude tail: int_T^inf (4/t^2) (1/2pi) log(4t/2pi) dt
    tail = 4 * (math.log(4 * T / (2 * math.pi)) + 1) / (2 * math.pi * T)
    return 1 / (M * Mp) - s - tail
log("(2) transversality: D := E_M/M - E_M'/M' against [kappa*lam, kappa*(2/(MM') - lam)]")
for M, Mp in ((0.4, 0.5), (0.5, 0.6), (0.5, 1.0), (0.25, 0.5)):
    D = Etheta(M) / M - Etheta(Mp) / Mp
    k = Mp - M; l = lam(M, Mp)
    sel = xs >= 1e6
    log("    (M,M') = (%.2f,%.2f): lam = %.4f, interval [%.4f, %.4f];  observed on x >= 1e6: [%.4f, %.4f]"
        % (M, Mp, l, k * l, k * (2 / (M * Mp) - l), D[sel].min(), D[sel].max()))
out.close()
