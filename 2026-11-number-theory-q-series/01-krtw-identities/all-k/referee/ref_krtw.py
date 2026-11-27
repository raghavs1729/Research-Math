"""Referee check (v): KRTW identities (ii), (iii) as q-series, own code, exact integers mod q^D.

LHS(iii) = K_k(1), LHS(ii) = K_k(q), computed from the 2k-fold sum by the transfer recursion
  F_1(r,s) = 1/(q)_r,  F_{j+1}(R,S) = sum_{r<=R, s>=S} w(r,s) F_j(r,s) / ((q)_{R-r} (q^2;q^2)_{s-S}),
  LHS = sum_{r,s} w(r,s) F_k(r,s) / (q^2;q^2)_s,
  w(r,s) = q^{(r+s)^2+s^2} (iii)  or  q^{(r+s)^2+s^2+r+2s} (ii);
RHS = (q^a1,...,q^a6; q^{2k+4})_inf / ((q^2;q^2)_inf (q;q)_inf).
"""
import sys, math
from flint import fmpz_poly

D = int(sys.argv[1]) if len(sys.argv) > 1 else 200
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10


def tr(p):
    c = p.coeffs()
    return fmpz_poly(c[:D]) if len(c) > D else p


def partitions_series(step):
    """1/(q^step;q^step)_inf mod q^D via Euler's pentagonal recurrence on p(n)"""
    n_max = D // step + 1
    p = [0] * (n_max + 1)
    p[0] = 1
    for n in range(1, n_max + 1):
        s, k = 0, 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > n:
                break
            sg = 1 if k % 2 else -1
            s += sg * p[n - g1]
            if g2 <= n:
                s += sg * p[n - g2]
            k += 1
        p[n] = s
    c = [0] * D
    for n in range(n_max + 1):
        if n * step < D:
            c[n * step] = p[n]
    return fmpz_poly(c)


def inv_poch(step, n):
    """1/(q^step;q^step)_n mod q^D"""
    p = fmpz_poly([1])
    for i in range(1, n + 1):
        e = step * i
        if e >= D:
            break
        geo = [0] * D
        for j in range(0, D, e):
            geo[j] = 1
        p = tr(p * fmpz_poly(geo))
    return p


def mono(e):
    if e >= D:
        return fmpz_poly([])
    return fmpz_poly([0] * e + [1])


def lhs(k, which):
    Rm = int(math.isqrt(D)) + 1
    Sm = int(math.isqrt(D // 2)) + 1
    iq = [inv_poch(1, n) for n in range(Rm + 2)]
    iQ = [inv_poch(2, n) for n in range(Sm + 2)]

    def wexp(r, s):
        e = (r + s)**2 + s * s
        return e + (r + 2 * s if which == "ii" else 0)
    F = {(r, s): iq[r] for r in range(Rm + 1) for s in range(Sm + 1)}
    for _ in range(k - 1):
        WF = {(r, s): tr(F[(r, s)] * mono(wexp(r, s))) for (r, s) in F}
        ES = {}
        for (R, S) in F:
            acc = fmpz_poly([])
            for s in range(S, Sm + 1):
                acc += tr(WF[(R, s)] * iQ[s - S])
            ES[(R, S)] = tr(acc)
        G = {}
        for (R, S) in F:
            acc = fmpz_poly([])
            for r in range(R + 1):
                acc += tr(ES[(r, S)] * iq[R - r])
            G[(R, S)] = tr(acc)
        F = G
    tot = fmpz_poly([])
    for (r, s) in F:
        tot += tr(tr(F[(r, s)] * mono(wexp(r, s))) * iQ[s])
    return tr(tot)


def rhs(k, which):
    m = 2 * k + 4
    exps = [1, 2, 2 * k + 2, 2 * k + 3, 2 * k + 4, 2 * k + 4] if which == "ii" else \
           [k + 1, k + 2, k + 2, k + 3, 2 * k + 4, 2 * k + 4]
    num = fmpz_poly([1])
    for a in exps:
        j = 0
        while a + j * m < D:
            num = tr(num * (fmpz_poly([1]) - mono(a + j * m)))
            j += 1
    return tr(tr(num * partitions_series(1)) * partitions_series(2))


if __name__ == "__main__":
    print(f"KRTW (ii), (iii): series comparison mod q^{D}")
    allok = True
    for k in range(1, KMAX + 1):
        res = []
        for which in ("ii", "iii"):
            L, Rr = lhs(k, which), rhs(k, which)
            ok = (L == Rr)
            allok &= ok
            res.append(f"({which}) {'OK' if ok else 'FAIL'}")
        print(f"k={k}: " + "  ".join(res))
    # negative control: (iii) with k+1 on the right
    L = lhs(4, "iii")
    print("negative control (iii) k=4 vs RHS with k=5:", "detected" if L != rhs(5, "iii") else "NOT detected")
    print("ALL OK" if allok else "SOME FAILED")
