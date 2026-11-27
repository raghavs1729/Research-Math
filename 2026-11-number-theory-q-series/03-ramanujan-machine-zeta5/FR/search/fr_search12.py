"""Systematic search, degree (1,2):  a_n = A n + k,  b_n = B n^2 + C n + E  (integers), for factorial reduction (FR)
outside the families of Conjecture 1.4 of arXiv:2111.04468, and for irreducible b (Conjecture 1.1, "only if").
s(n) := max_{n-3<=m<=n} (ln GCD_m - ln m!)/m.   FR  <=>  s bounded below;  non-FR typically s(n) ~ -ln n + O(1).
Stage 1: n = 40, 80;  keep if s(80) > s(40) - 0.40.   Stage 2 (survivors): n = 160, 320, 640.
Classes:  DEG  (b_n = 0 for some n >= 1: terminating fraction),
          F6   (b reducible over Q and A^2 + 4B a nonzero square),   F7  (b reducible and 2kB = A(B + C)),
          IRR  (b irreducible over Q),  OUT (b reducible, not DEG, not F6, not F7)."""
import sys, math
from math import gcd, isqrt, lgamma, log
def issq(x): return x >= 0 and isqrt(x)**2 == x
def classify(A, k, B, C, E):
    disc = C * C - 4 * B * E
    red = issq(disc)
    cls = []
    if red:
        r = isqrt(disc)
        for num in (-C + r, -C - r):
            den = 2 * B
            if num % den == 0 and num // den >= 1: cls.append('DEG')
        if issq(A * A + 4 * B) and A * A + 4 * B != 0: cls.append('F6')
        if 2 * k * B == A * (B + C): cls.append('F7')
        if not cls: cls.append('OUT')
    else:
        cls.append('IRR')
        if issq(A * A + 4 * B) and A * A + 4 * B != 0: cls.append('i6')
        if 2 * k * B == A * (B + C): cls.append('i7')
    return cls
def svals(A, k, B, C, E, marks):
    p0, p1 = 1, k; q0, q1 = 0, 1
    out = []; best = None; mi = 0
    for n in range(1, marks[-1] + 1):
        a = A * n + k; b = (B * n + C) * n + E
        p0, p1 = p1, a * p1 + b * p0; q0, q1 = q1, a * q1 + b * q0
        if n > marks[mi] - 4:
            g = gcd(p1, q1)
            if g:
                v = (log(g) - lgamma(n + 1)) / n
                best = v if best is None or v > best else best
        if n == marks[mi]:
            out.append(best); best = None; mi += 1
    return out
if __name__ == '__main__':
    RA, RK, RB, RC = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    surv = []; cnt = {}
    for A in range(1, RA + 1):
        for k in range(-RK, RK + 1):
            for B in list(range(-RB, 0)) + list(range(1, RB + 1)):
                for C in range(-RC, RC + 1):
                    for E in range(-RC, RC + 1):
                        if gcd(gcd(A, k), gcd(gcd(B, C), E)) != 1 and gcd(A, k) ** 2 and all(x % 4 == 0 for x in (B, C, E)) and A % 2 == 0 and k % 2 == 0:
                            continue          # inflation by 2 of a smaller case (cheap partial dedup)
                        cls = classify(A, k, B, C, E)
                        key = cls[0]; cnt[key] = cnt.get(key, 0) + 1
                        s = svals(A, k, B, C, E, [40, 80])
                        if s[0] is None or s[1] is None: 
                            surv.append((A, k, B, C, E, cls, s)); continue
                        if s[1] > s[0] - 0.40: surv.append((A, k, B, C, E, cls, s))
    print("cases by class:", cnt, " stage-1 survivors:", len(surv))
    tab = {}
    for (A, k, B, C, E, cls, s) in surv:
        tab.setdefault(tuple(cls), []).append((A, k, B, C, E, s))
    for c, L in sorted(tab.items()): print(c, len(L))
    import pickle; pickle.dump(surv, open('stage1_%d_%d_%d_%d.pkl' % (RA, RK, RB, RC), 'wb'))
