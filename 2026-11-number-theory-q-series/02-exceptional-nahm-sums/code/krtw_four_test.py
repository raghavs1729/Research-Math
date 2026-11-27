"""Numerical test of the four conjectural Andrews-Gordon-type identities of Kanade-Russell-Tsuchioka-Warnaar
(arXiv:2404.03851, display after Conjecture 4.9) at high order, for k = 1..KMAX.

    W(r,s) = prod_{i=1}^k q^{(r_i+s_i)^2+s_i^2} / ((q;q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}}),   r_{k+1}=s_0=0
    (i)   sum q^{sum s_i} W          = (q,q^{k+1},q^{k+2};q^{k+2})_inf / ((q;q^2)_inf (q)_inf)
    (ii)  sum q^{sum(r_i+2s_i)} W    = (q,q^2,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4};q^{2k+4})_inf / ((q^2;q^2)_inf (q)_inf)
    (iii) sum W                      = (q^{k+1},q^{k+2},q^{k+2},q^{k+3},q^{2k+4},q^{2k+4};q^{2k+4})_inf / (...)
    (iv)  sum Omega W                = (q^k,q^{k+1},q^{k+3},q^{k+4},q^{2k+4},q^{2k+4};q^{2k+4})_inf / (...)
          Omega = sum_{i<k} q^{r_i+2s_i}(1-q^{2s_{i+1}-2s_i}) + q^{r_k+2s_k}.

Method (exact, mod q^N): transfer recursion over the levels.  State after level i: (R, s) = (r_{i+1}, s_i);
    Phi_i(R,s) = sum_{r>=R, s'<=s} w_i(r,s) / ((q)_{r-R} (q^2;q^2)_{s-s'}) * Phi_{i-1}(r,s'),   Phi_0(r,s') = [s'=0],
and the sum is sum_s Phi_k(0,s).  This is just the definition, summed level by level (all terms with
exponent < N are kept: w_i(r,s) = q^{(r+s)^2+s^2+linear} bounds r and s).
For (iv), the i-th term of Omega multiplies level i by q^{r_i+2s_i} and replaces 1/(q^2;q^2)_{b} at level i+1
(b = s_{i+1}-s_i) by (1-q^{2b})/(q^2;q^2)_b = 1/(q^2;q^2)_{b-1} (0 for b = 0).

usage: python3 krtw_four_test.py [N] [KMAX]
"""
import sys
from flint import fmpz_poly
import krtw_ser as K
from krtw_ser import tr, mono, iqpoch

def multisum(k, lin, N, omega_term=None):
    """lin = (a,b): extra linear exponent a*r_i + b*s_i at every level.
    omega_term = i (1..k): insert the i-th term of Omega."""
    rmax = 2
    while rmax * rmax - 2 * rmax < N:      # margin: also valid for linear terms >= -(r+s) per level
        rmax += 1
    smax = rmax
    # Phi[(R,s)]
    Phi = {(r, 0): fmpz_poly([1]) for r in range(rmax + 1)}   # Phi_0(r, s') = [s'=0]
    for i in range(1, k + 1):
        new = {}
        for s in range(smax + 1):
            # g(r) = sum_{s'<=s} w(r,s) Phi_{i-1}(r,s') / (q^2;q^2)_{s-s'}
            g = {}
            for r in range(rmax + 1):
                e = (r + s) ** 2 + s * s + lin[0] * r + lin[1] * s
                if omega_term == i:
                    e += r + 2 * s
                if e >= N:
                    if r + s >= 1:
                        break          # e is increasing in r from here on
                    continue
                acc = fmpz_poly([0])
                for sp in range(s + 1):
                    p = Phi.get((r, sp))
                    if p is None:
                        continue
                    b = s - sp
                    if omega_term is not None and omega_term == i - 1:
                        if b == 0:
                            continue
                        den = iqpoch(b - 1, 2)
                    else:
                        den = iqpoch(b, 2)
                    acc = acc + tr(p * den)
                if acc != 0:
                    g[r] = tr(acc * mono(e))
            if not g:
                continue
            for R in range(rmax + 1):
                acc = fmpz_poly([0])
                for r, v in g.items():
                    if r >= R:
                        acc = acc + tr(v * iqpoch(r - R))
                if acc != 0:
                    new[(R, s)] = acc
        Phi = new
    tot = fmpz_poly([0])
    for (R, s), v in Phi.items():
        if R == 0:
            tot = tot + v
    return tot

def prod_side(factors_num, mod, N, den):
    """prod_{a in factors_num} (q^a; q^mod)_inf  /  prod_{(a,m) in den} (q^a;q^m)_inf   mod q^N (a>=1)"""
    c = [0] * N
    c[0] = 1
    for a in factors_num:
        n = a
        while n < N:
            for j in range(N - 1, n - 1, -1):
                c[j] -= c[j - n]
            n += mod
    for (a, m) in den:
        n = a
        while n < N:
            for j in range(n, N):
                c[j] += c[j - n]
            n += m
    return fmpz_poly(c)

def products(k, N):
    m = 2 * k + 4
    d2 = [(2, 2), (1, 1)]
    return {
        "i": prod_side([1, k + 1, k + 2], k + 2, N, [(1, 2), (1, 1)]),
        "ii": prod_side([1, 2, 2 * k + 2, 2 * k + 3, m, m], m, N, d2),
        "iii": prod_side([k + 1, k + 2, k + 2, k + 3, m, m], m, N, d2),
        "iv": prod_side([k, k + 1, k + 3, k + 4, m, m], m, N, d2),
    }

def first_diff(a, b, N):
    ca = a.coeffs() + [0] * N
    cb = b.coeffs() + [0] * N
    for n in range(N):
        if ca[n] != cb[n]:
            return n
    return None

def brute(k, lin, N, omega=False):
    """the same sums by direct enumeration of (r,s) (independent of multisum; small k and N only)"""
    import itertools
    tot = fmpz_poly([0])
    m = 2
    while m * m - 2 * m < N:
        m += 1
    for r in itertools.product(range(m + 1), repeat=k):
        if any(r[i] < r[i + 1] for i in range(k - 1)):
            continue
        for s in itertools.product(range(m + 1), repeat=k):
            if any(s[i] > s[i + 1] for i in range(k - 1)):
                continue
            e = sum((r[i] + s[i]) ** 2 + s[i] ** 2 + lin[0] * r[i] + lin[1] * s[i] for i in range(k))
            if e >= N:
                continue
            c = fmpz_poly([1])
            for i in range(k):
                c = tr(c * iqpoch(r[i] - (r[i + 1] if i < k - 1 else 0)))
                c = tr(c * iqpoch(s[i] - (s[i - 1] if i > 0 else 0), 2))
            if omega:
                om = fmpz_poly([0])
                for i in range(k - 1):
                    om = om + mono(r[i] + 2 * s[i]) - mono(r[i] + 2 * s[i] + 2 * s[i + 1] - 2 * s[i])
                om = om + mono(r[k - 1] + 2 * s[k - 1])
                c = tr(c * om)
            tot = tot + tr(c * mono(e))
    return tot

def selftest():
    """transfer recursion == brute force for k <= 3 mod q^30 (all four sums, and the x = 1/q linear terms)"""
    N0 = 30
    K.setN(N0)
    ok = True
    for k in (1, 2, 3):
        for lin in ((0, 1), (1, 2), (0, 0), (-1, -1)):
            ok &= (brute(k, lin, N0) == multisum(k, lin, N0))
        om = fmpz_poly([0])
        for i in range(1, k + 1):
            om = om + multisum(k, (0, 0), N0, omega_term=i)
        ok &= (brute(k, (0, 0), N0, omega=True) == om)
    # a deliberately wrong product must be detected
    P = products(3, N0)
    ok &= (first_diff(multisum(3, (0, 0), N0), P["iv"], N0) is not None)
    return ok

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    print("self-test (transfer recursion vs brute-force enumeration, k <= 3, mod q^30):", "OK" if selftest() else "FAIL")
    K.setN(N)
    print("KRTW four Andrews-Gordon-type identities, exact comparison mod q^%d" % N)
    allok = True
    for k in range(1, KMAX + 1):
        P = products(k, N)
        S1 = multisum(k, (0, 1), N)
        S2 = multisum(k, (1, 2), N)
        S3 = multisum(k, (0, 0), N)
        S4 = fmpz_poly([0])
        for i in range(1, k + 1):
            S4 = S4 + multisum(k, (0, 0), N, omega_term=i)
        res = []
        for name, Sx in (("i", S1), ("ii", S2), ("iii", S3), ("iv", S4)):
            d = first_diff(Sx, P[name], N)
            res.append("(%s) %s" % (name, "OK" if d is None else "FAILS at q^%d" % d))
            allok &= d is None
        print("k=%d: " % k + "   ".join(res), flush=True)
    print("ALL OK" if allok else "SOME FAIL")
