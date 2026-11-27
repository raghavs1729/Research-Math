"""Referee check (i): exact commutators [T,H'], [T,G], [T#,H'#], [T#,G#] on point masses, and the
transposition formulas (5.1) H'# = W^-1 H'^t W, G# = W^-1 G^t W entrywise, at rational points (q,x)."""
import itertools, sys
from flint import fmpq
from ref_common import Point


def box(R0, R1, S0, S1):
    return [(R, S) for R in range(R0, R1 + 1) for S in range(S0, S1 + 1)]


def check_point(q, x, rmax=4, smax=3, Rbig=9, Sbig=7):
    P = Point(q, x)
    bad = 0
    n = 0
    for r0, s0 in itertools.product(range(rmax + 1), range(smax + 1)):
        d = {(r0, s0): fmpq(1)}
        # --- T with H', G
        Bev = box(0, Rbig, 0, s0 + 1)
        for op in (P.Hp, P.Gop):
            A = P.T(op(d, box(0, r0 + 3, 0, s0 + 2)), Bev)
            Td = P.T(d, box(0, Rbig + 1, 0, s0 + 2))
            B = op(Td, Bev)
            for k in Bev:
                n += 1
                if A[k] != B[k]:
                    bad += 1
        # --- T# with H'#, G#
        Bev = box(0, r0 + 1, 0, Sbig)
        for op in (P.Hps, P.Gs):
            A = P.Tsharp(op(d, box(0, r0 + 2, 0, s0 + 2)), Bev)
            Tsd = P.Tsharp(d, box(0, r0 + 3, 0, Sbig + 1))
            B = op(Tsd, Bev)
            for k in Bev:
                n += 1
                if A[k] != B[k]:
                    bad += 1
    # --- transposition: H'#(a,b) = w(b) H'(b,a) / w(a), same for G
    pts = box(0, rmax, 0, smax)
    for (opT, opS) in ((P.Hp, P.Hps), (P.Gop, P.Gs)):
        for a in pts:
            da = {a: fmpq(1)}
            colA = opT(da, pts)          # H'(b, a) for all b
            for b in pts:
                db = {b: fmpq(1)}
                lhs = opS(db, [a])[a]    # H'#(a, b)
                rhs = P.w(*b) * colA[b] / P.w(*a)
                n += 1
                if lhs != rhs:
                    bad += 1
    return n, bad


def negative_control(q, x):
    """Perturb H' (x^2 -> 2x^2 in the u-term) and check the commutator is detected."""
    P = Point(q, x)
    orig = P.Hp

    def Hbad(F, bx):
        out = orig(F, bx)
        for (R, S) in bx:
            out[(R, S)] -= P.x**2 * P.q**R * P.q**(2 * S) * P.get(F, R, S + 1)
        return out
    d = {(1, 1): fmpq(1)}
    Bev = box(0, 6, 0, 2)
    A = P.T(Hbad(d, box(0, 4, 0, 3)), Bev)
    B = Hbad(P.T(d, box(0, 7, 0, 3)), Bev)
    return sum(1 for k in Bev if A[k] != B[k])


if __name__ == "__main__":
    for (q, x) in ((fmpq(2, 7), fmpq(3, 11)), (fmpq(-3, 5), fmpq(5, 13))):
        n, bad = check_point(q, x)
        print(f"(q,x)=({q},{x}): {n} exact comparisons, {bad} mismatches")
        print(f"   negative control (perturbed H'): {negative_control(q, x)} mismatches (should be > 0)")
