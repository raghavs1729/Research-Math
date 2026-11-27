#!/usr/bin/env python3
"""Hassani's inequality as printed ([H21], Section 4, question 3):
      H_{q,a}(x):   pi(x;q,a)^2 < phi(q) (e x/log x) pi(x/e;q,a).
Claim (unconditional, q in {3,4,5,8,12}): H fails exactly for 1 < x < e*p1, p1 = least prime = a mod q
(where the right side is 0).  Three parts:
 (1) x <= 2e5: exact structure of the failure set (smallx.analyse with phi replaced by 1/phi; Arb).
 (2) 7 <= p <= X_s: at a prime p of the class H fails iff rho(p) >= phi^2 (rho as in apscan); the merged
     scan gives the exact maximum of rho over all scanned primes with pi(p/e) >= 1, and the stepping pass
     gives rho < 1 on unflagged rows.  Lemma 1 holds verbatim for H (same monotonicity).
 (3) x >= 3e13 (L >= 31, x/e >= 8e9): [BMOR] Thm 1.3, |pi(x;q,a) - Li(x)/phi| < x/(840 log^2 x) for
     x >= 8e9, 3 <= q <= 10^4, and Lemma A(i) give, after division by x^2/L^2, the sufficient condition
        (a+(L)/phi + k/L)^2 < 1 - phi k L/(L-1)^2 - e li(2) L e^{-L},  k = 1/840,
     whose left side decreases and right side increases in L; checked at L = 31 in Arb.
Usage: hassani_printed.py merged1 [merged2 ...]"""
import sys
from math import factorial
from flint import arb, ctx
import smallx
ctx.prec = 200
CQ = [1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12]
CA = [0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11]
PHI = [1, 2,2, 2,2, 4,4,4,4, 4,4,4,4, 4,4,4,4]
# (1)
P = smallx.primes_upto(200000)
orig_root = smallx.root
print("(1) exact failure set of H for x <= 2e5:")
for c in range(1, 17):
    q, a, phi = CQ[c], CA[c], PHI[c]
    save = smallx.PHI[q]
    smallx.PHI[q] = arb(1) / phi
    smallx.root = lambda n, m, ph, lo, phi=phi: orig_root(n, m, 1.0 / phi, lo)
    comps, nfail, nstart, failing, active = smallx.analyse(q, a, 200000, P)
    smallx.PHI[q] = save; smallx.root = orig_root
    p1 = next(p for p in P if p % q == a and q % p != 0)
    assert len(comps) == 1 and not active and comps[0][2] == "e*%d" % p1, (q, a, comps)
    print(f"   ({q},{a}): F_H ∩ (1, 2e5] = (1, e*{p1}) = (1, {comps[0][1]:.4f})")
# (2)
best = {}
for fn in sys.argv[1:]:
    for l in open(fn):
        w = l.split()
        if not w or w[0] != "BIN": continue
        c = int(w[1]); rho = float(w[6])
        if c not in best or rho > best[c][0]: best[c] = (rho, int(w[7]))
    for l in open(fn):
        if l.startswith("# coverage"): print("(2)", l.strip())
for c in range(1, 17):
    assert best[c][0] < PHI[c] ** 2 * 0.9
    print(f"   ({CQ[c]},{CA[c]}): max rho over scanned primes with pi(p/e) >= 1 is {best[c][0]:.4f} (at p = {best[c][1]}) < phi^2 = {PHI[c]**2}")
# (3)
L = arb(31); k = arb(1) / 840; li2 = arb(2).log().ei()
aplus = L * (sum(arb(factorial(j)) / L ** (j + 1) for j in range(14)) + arb(5) * factorial(14) / L ** 15)
print("(3) tail, L >= 31:")
for phi in (2, 4):
    lhs = (aplus / phi + k / L) ** 2
    rhs = 1 - phi * k * L / (L - 1) ** 2 - arb(1).exp() * li2 * L * (-L).exp()
    assert lhs < rhs
    print(f"   phi = {phi}: lhs = {lhs.str(8)} < rhs = {rhs.str(8)}")
