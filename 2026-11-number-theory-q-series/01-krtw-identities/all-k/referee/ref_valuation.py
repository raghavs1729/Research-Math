"""Referee check of Lemma 6.2 (a), (b): x-adic valuations of the position-space eigenvectors.
ord_x f is estimated from exact values at x1 = 10^-25, x2 = 10^-50 (q fixed rational):
  ord ~ log(|f(x1)|/|f(x2)|) / log(x1/x2)   (exact rationals; the error is O(10^-25 * C)).
Checks: ord phi_{N,M}(R,S) >= 2(M-S) (Lemma 6.2(a)), ord phi#_{N,M}(R,S) >= 0 (6.2(b)), and reports whether the
bound 2(M-S) is attained (sharpness)."""
import math, sys
from flint import fmpq
from ref_common import Point
from ref_eigen import eig_T, eig_Tsharp


def lg(v):
    v = abs(v)
    return math.log10(int(v.p)) - math.log10(int(v.q))


def est_ord(a, b, step):
    if a == 0 and b == 0:
        return None
    return (lg(a) - lg(b)) / step


def main(qv, Nmax):
    x1, x2 = fmpq(1, 10**25), fmpq(1, 10**50)
    P1, P2 = Point(qv, x1), Point(qv, x2)
    worst_a, worst_b = 1e9, 1e9
    sharp = 0
    tot = 0
    for N in range(Nmax + 1):
        for M in range(N + 1):
            L = N - M
            Rmax = L + 2 * M + 4
            f1, f2 = eig_T(P1, N, M, Rmax), eig_T(P2, N, M, Rmax)
            for k in f1:
                o = est_ord(f1[k], f2[k], 25)
                if o is None:
                    continue
                tot += 1
                worst_a = min(worst_a, o - 2 * (M - k[1]))
                if abs(o - 2 * (M - k[1])) < 0.01:
                    sharp += 1
            Smax = M + L + 4
            g1, g2 = eig_Tsharp(P1, N, M, Smax), eig_Tsharp(P2, N, M, Smax)
            for k in g1:
                o = est_ord(g1[k], g2[k], 25)
                if o is None:
                    continue
                worst_b = min(worst_b, o)
    print(f"q={qv}, N<={Nmax}: min over (N,M,R,S) of [ord phi(R,S) - 2(M-S)] = {worst_a:.4f} "
          f"(bound attained at {sharp}/{tot} nonzero entries);  min ord phi#(R,S) = {worst_b:.4f}")


if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for qv in (fmpq(2, 7), fmpq(-3, 5)):
        main(qv, Nmax)
