"""Check (I'):  sum_{L,M} t_{L,0}(x q^{2M}) gamma^{L,M}_d(x) == (x^2q^2;q^2)_inf alpha^{(0)}_d(x)
as power series in x (to x^XMAX), exactly, at a few rational values of q.
alpha^{(0)}_d = (-1)^d x^d q^{d^2 + d(d-1)/2} (x;q)_d (1 - x q^{2d}) / ((q;q)_d (1-x))."""
import sys
from fractions import Fraction
import flint
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from tower import tower
from check_II import target

XMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
DMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4


def poly_at_q(p, q0):
    """fmpz_mpoly in (x,q) -> fmpq_poly in x at q = q0"""
    coeffs = {}
    for (i, j), c in p.to_dict().items():
        i = int(i); j = int(j)
        coeffs[i] = coeffs.get(i, Fraction(0)) + Fraction(int(c)) * q0 ** j
    deg = max(coeffs) if coeffs else 0
    return flint.fmpq_poly([flint.fmpq(coeffs.get(i, Fraction(0)).numerator, coeffs.get(i, Fraction(0)).denominator) for i in range(deg + 1)])


def series(rf, q0, n):
    a = poly_at_q(rf.n, q0)
    b = poly_at_q(rf.d, q0)
    A = [a[i] for i in range(n)]
    B = [b[i] for i in range(n)]
    assert B[0] != 0
    inv = [flint.fmpq(0)] * n
    inv[0] = 1 / B[0]
    for k in range(1, n):
        s = flint.fmpq(0)
        for j in range(1, k + 1):
            if B[j] != 0:
                s += B[j] * inv[k - j]
        inv[k] = -s / B[0]
    return mulser(A, inv, n)


def alpha0(d):
    s = -1 if d % 2 else 1
    if d == 0:
        return RF(1)
    return RF(s) * mono(d, d * d + d * (d - 1) // 2) * poch(X, d) * RF(ONE - X * Q ** (2 * d)) / (qpoch(d) * RF(ONE - X))


def x2Qinf_series(q0, n):
    # (x^2 q^2; q^2)_inf = sum_j (-1)^j x^{2j} q^{j(j+1)} / (q^2;q^2)_j
    out = [flint.fmpq(0)] * n
    j = 0
    while 2 * j < n:
        c = Fraction((-1) ** j) * q0 ** (j * (j + 1))
        for i in range(1, j + 1):
            c /= (1 - q0 ** (2 * i))
        out[2 * j] += flint.fmpq(c.numerator, c.denominator)
        j += 1
    return out


def mulser(a, b, n):
    out = [flint.fmpq(0)] * n
    for i in range(n):
        if a[i] == 0:
            continue
        for j in range(n - i):
            out[i + j] += a[i] * b[j]
    return out


if __name__ == '__main__':
    for q0 in [Fraction(1, 3), Fraction(-2, 5), Fraction(3, 7)]:
        P = x2Qinf_series(q0, XMAX + 1)
        for d in range(0, DMAX + 1):
            rhs = mulser(P, series(alpha0(d), q0, XMAX + 1), XMAX + 1)
            lhs = [flint.fmpq(0)] * (XMAX + 1)
            for L in range(0, d + 1):
                M = 0
                while L + 2 * M <= XMAX:
                    if L + 2 * M >= d:
                        g = tower(L, M)
                        if d in g:
                            t = target(L, 0).subs_x(X * Q ** (2 * M)) * g[d]
                            s = series(t, q0, XMAX + 1)
                            lhs = [lhs[i] + s[i] for i in range(XMAX + 1)]
                    M += 1
            ok = all(lhs[i] == rhs[i] for i in range(XMAX + 1))
            print('q=%s d=%d: %s' % (q0, d, 'OK' if ok else 'FAIL'))
            if not ok:
                print('  lhs', lhs[:8]); print('  rhs', rhs[:8])
            sys.stdout.flush()
