"""Second audit: central non-vanishing L(1/2,chi) != 0, certified in Arb, for every primitive
character of every conductor f | q, for a sample of moduli.  Own characters (achars.py);
L(1/2,chi) = f^{-1/2} sum_k chi(k) zeta(1/2, k/f) with flint's Hurwitz zeta (acb.zeta(s, a)).
(The author used flint's dirichlet L; this uses the Hurwitz route and an independent character table.)
"""
import sys
import flint
from flint import acb, arb, fmpq
from sympy import divisors
from achars import primitive_characters

flint.ctx.prec = 96


def run(q):
    worst = None; count = 0; bad = 0
    for f in divisors(q):
        if f < 3:
            continue
        half = acb(fmpq(1, 2))
        Z = [half.zeta(acb(fmpq(k, f))) for k in range(1, f + 1)]
        norm = arb(f).rsqrt()
        for vals in primitive_characters(f):
            s = acb(0)
            for k in range(1, f + 1):
                v = vals[k % f]
                if v is None:
                    continue
                s += acb.exp_pi_i(acb(fmpq(2 * v.numerator, v.denominator))) * Z[k - 1]
            L = s * norm
            count += 1
            absL = abs(L)
            if not (absL > 0):
                bad += 1
                print('  NOT certified nonzero:', f, vals[:6], L)
            lo = float(absL.lower())
            if worst is None or lo < worst[0]:
                worst = (lo, f)
    print(f'q={q}: {count} primitive characters of conductor f|q, f>=3; uncertified {bad}; min lower bound |L(1/2)| = {worst[0]:.4e} (f={worst[1]})')


if __name__ == '__main__':
    for q in [int(x) for x in sys.argv[1:]]:
        run(q)


def run_all(F):
    """every primitive character of every conductor 3 <= f <= F"""
    import math
    count = 0; bad = 0; worst = (1e9, 0)
    for f in range(3, F + 1):
        half = acb(fmpq(1, 2))
        P = primitive_characters(f)
        if not P:
            continue
        Z = [half.zeta(acb(fmpq(k, f))) if math.gcd(k, f) == 1 else None for k in range(1, f + 1)]
        norm = arb(f).rsqrt()
        for vals in P:
            s = acb(0)
            for k in range(1, f + 1):
                v = vals[k % f]
                if v is None:
                    continue
                s += acb.exp_pi_i(acb(fmpq(2 * v.numerator, v.denominator))) * Z[k - 1]
            absL = abs(s * norm)
            count += 1
            if not (absL > 0):
                bad += 1; print('  NOT certified nonzero:', f, absL)
            lo = float(absL.lower())
            if lo < worst[0]:
                worst = (lo, f)
    print(f'all conductors 3..{F}: {count} primitive characters; uncertified {bad}; min |L(1/2)| >= {worst[0]:.4e} at f={worst[1]}')
