#!/usr/bin/env python3
"""
t2_lll_check.py -- independent re-evaluation of the lattice witnesses (raw/t2_zeros_lll_*.jsonl).

For each record (q, a, y_0 exact dyadic, claimed F_Z(y_0)):
 (1) the refined brackets of every character are re-read and checked to be positive and pairwise
     disjoint, and each is checked to lie inside a bracket of raw/zeros_cache/ (the first-stage
     certificate);
 (2) F_Z(y_0) is recomputed with mpmath (a different multiprecision library; 1100 bits) from the
     bracket midpoints, with phases reduced mod 2 pi explicitly and the characters' values taken
     from their exact exponents, and compared with the Arb ball of the record;
 (3) the remainder r = L_1(q) - sum of located weights is recomputed and the witness inequality
     F_Z(y_0) - r > (2/3) N re-checked.
This guards against a slip in the precision handling at |y_0| ~ 10^50..10^130; it does not
re-certify the sign changes (t2_refine.py does that in Arb).
"""
import sys, json, glob, math
import mpmath as mp
from fractions import Fraction
from flint import arb, acb, ctx, dirichlet_char, dirichlet_group
import t2_witness as TW
import l1const as LC

mp.mp.prec = 1100
ctx.prec = 128


def check(rec):
    q, a, T, P = rec["q"], rec["a"], int(rec["T"]), rec["P"]
    man, exp = rec["y0_exact"].split(" * 2^(")
    y0 = mp.mpf(int(man)) * mp.mpf(2) ** int(exp[:-1])
    table = TW.char_table(q)
    G = int(dirichlet_group(q).exponent()) if q > 1 else 1
    handled = set()
    tot = mp.mpf(0)
    found = mp.mpf(0)
    nbr = 0
    def hp(f, m):
        ctx.prec = 2200
        new = [(arb(x[0]), arb(x[1])) for x in json.load(open(f"raw/zeros_hp/z_{f}_{m}_T{T}_P{P}.json"))]
        old = [arb(x) for x in json.load(open(f"raw/zeros_cache/z_{f}_{m}_T{T}.json"))]
        assert len(new) == len(old)
        prev = None
        for (mid, rad), ob in zip(new, old):
            lo, hi = mid - rad, mid + rad
            assert lo > 0 and lo > ob.lower() and hi < ob.upper(), (f, m)
            if prev is not None:
                assert prev < lo
            prev = hi
        ctx.prec = 128
        return [mp.mpf(x[0]) for x in json.load(open(f"raw/zeros_hp/z_{f}_{m}_T{T}_P{P}.json"))]
    for chi, f, m, pr in table:
        key = (f, m); cj = (f, pow(m, -1, f)) if f > 1 else (1, 1)
        if key in handled:
            continue
        handled.add(key); handled.add(cj)
        e = Fraction(int(chi.chi_exponent(a)), G) if q > 1 else Fraction(0)      # chi(a) = exp(2 pi i e)
        wbar = mp.expjpi(-2 * mp.mpf(e.numerator) / e.denominator)
        if f == 1:
            zs = [mp.im(mp.zetazero(k)) for k in range(1, 400)]
            zs = [z for z in zs if z <= T]
        elif key == cj:
            zs = hp(f, m)
        else:
            zs = hp(f, m) + [-g for g in hp(*cj)]
        for g in zs:
            rho = mp.mpc(0.5, g)
            c = -wbar / (rho * (rho + 1))
            ph = mp.fmod(g * y0, 2 * mp.pi)
            tot += 2 * abs(c) * mp.cos(ph + mp.arg(c))
            found += 2 * abs(c)
            nbr += 1
    L1 = arb(0)
    for chi, f, m, pr in table:
        L1 += LC.S_integral(f, m)
    rem = mp.mpf(L1.mid().str(40, radius=False)) - found
    claimed = arb(rec["F_Z_at_y0"])
    diff = abs(tot - mp.mpf(claimed.mid().str(30, radius=False)))
    ok = (tot - rem) > mp.mpf(rec["threshold"]) and diff < mp.mpf(10) ** (-10) and nbr == rec["num_brackets"]
    print(f"({q},{a}): |y0| ~ 10^{mp.nstr(mp.log10(abs(y0)), 6)}; mpmath F_Z = {mp.nstr(tot, 15)}, Arb record {rec['F_Z_at_y0']}, |diff| = {mp.nstr(diff, 3)}; "
          f"r = {mp.nstr(rem, 8)}; F_Z - r = {mp.nstr(tot - rem, 10)} vs (2/3)N = {rec['threshold']:.6f}; cosines {nbr}; {'OK' if ok else 'MISMATCH'}")
    return ok


if __name__ == "__main__":
    allok = True
    for fn in sorted(glob.glob("raw/t2_zeros_lll_*.jsonl")):
        for l in open(fn):
            rec = json.loads(l)
            if rec.get("witness"):
                allok &= check(rec)
    print("ALL OK" if allok else "SOME MISMATCH")
