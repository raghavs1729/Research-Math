#!/usr/bin/env python3
"""
t2_refine.py -- refine the certified zero brackets of raw/zeros_cache/ to high precision.

For the primitive character f.m and the cached brackets of zeros of Hardy's Z in (0, T]
(t2_zeros.py), each bracket [lo, hi] is replaced by a much narrower one INSIDE it:
a secant iteration on ball midpoints proposes x; the new bracket [x - d, x + d] is accepted only
if Z has certified opposite signs (Arb balls excluding 0) at its two ends and it lies inside the
old bracket.  So each new bracket again contains a zero on the critical line, and the brackets
stay pairwise disjoint.  Output: raw/zeros_hp/z_{f}_{m}_T{T}_P{P}.json  (list of [mid, rad] decimal
strings; the ball mid +- rad, rounded outward, contains the bracket).

Usage: python3 t2_refine.py f m T P        (P = target bits: bracket width about 2^-P)
"""
import sys, json, os
from flint import arb, acb, dirichlet_char, ctx


def main():
    f, m, T, P = (int(x) for x in sys.argv[1:5])
    os.makedirs("raw/zeros_hp", exist_ok=True)
    out_fn = f"raw/zeros_hp/z_{f}_{m}_T{T}_P{P}.json"
    if os.path.exists(out_fn):
        return
    ctx.prec = 128
    old = [arb(x) for x in json.load(open(f"raw/zeros_cache/z_{f}_{m}_T{T}.json"))]
    chi = dirichlet_char(f, m)
    ctx.prec = P + 96
    digits = int((P + 96) * 0.30103) + 5

    def Z(t):
        return chi.hardy_z(acb(t)).real

    def sgn(v):
        return 1 if v > 0 else (-1 if v < 0 else 0)

    res = []
    for b in old:
        lo, hi = arb(b.lower()), arb(b.upper())
        # secant iteration on midpoints, started from the ends of the old bracket
        x0, x1 = lo, hi
        f0, f1 = arb(Z(x0).mid()), arb(Z(x1).mid())
        target = arb(2) ** (-(P + 30))
        for _ in range(40):
            den = f1 - f0
            if den.contains(0):
                break
            x2 = arb((x1 - f1 * (x1 - x0) / den).mid())
            f2 = arb(Z(x2).mid())
            x0, f0, x1, f1 = x1, f1, x2, f2
            if abs(x1 - x0) < target:
                break
        x = x1
        d = arb(2) ** (-P)
        ok = False
        for _ in range(12):
            a, c = arb((x - d).mid()), arb((x + d).mid())
            s1, s2 = sgn(Z(a)), sgn(Z(c))
            if s1 * s2 < 0 and a > lo and c < hi:
                ok = True
                break
            d = d * 64
        assert ok, (f, m, b)
        nb = a.union(c)
        res.append([nb.mid().str(digits, radius=False), (nb.rad() * arb("1.0001")).upper().str(8, radius=False)])
    json.dump(res, open(out_fn + ".tmp", "w"))
    os.replace(out_fn + ".tmp", out_fn)
    sys.stderr.write(f"{f}.{m}: {len(res)} zeros refined to 2^-{P}\n")


if __name__ == "__main__":
    main()
