"""High-precision limit of a PCF  a_n = lead*atil(n), b_n = -(lead n^s)^2  (double characteristic root),
via the two formal solutions  f(t) = t^d sum_j c_j t^(-2j),  t = n+1/2, of
     (n+1)^s f(n+1) + n^s f(n-1) = atil(n) f(n)
and matching with the exactly iterated numerator/denominator sequences at n = N, N+1.
Numerical tool only (used to obtain many digits of a limit); nothing here is part of a proof."""
from fractions import Fraction as Fr
from math import comb
from mpmath import mp, mpf

def series_coeffs(s, at, d, J):
    half = [mpf(comb(s, k)) / 2**k for k in range(s + 1)]
    d = mpf(d.numerator) / d.denominator
    def Trow(e, Rmax):
        # T(e,R) for even R = 0..Rmax : 2*sum_k half[k]*C(e,R-k) - at[s-R]
        bins = [mpf(1)]
        for i in range(1, Rmax + 1):
            bins.append(bins[-1] * (e - i + 1) / i)
        row = {}
        for R in range(0, Rmax + 1, 2):
            tot = mpf(0)
            for k in range(0, min(s, R) + 1):
                tot += half[k] * bins[R - k]
            tot = 2 * tot - at.get(s - R, 0)
            row[R] = tot
        return row
    rows = [Trow(d, 2 * J + 2)]
    assert abs(rows[0][0]) < mpf(10) ** (-mp.dps + 20) and abs(rows[0][2]) < mpf(10) ** (-mp.dps + 20), rows[0][2]
    c = [mpf(1)]
    for j in range(1, J + 1):
        rows.append(Trow(d - 2 * j, 2 * J + 2 - 2 * j))
        rhs = -sum(c[i] * rows[i][2 * j + 2 - 2 * i] for i in range(j))
        piv = rows[j][2]
        if abs(piv) < mpf(10) ** (-mp.dps + 20):
            if abs(rhs) > mpf(10) ** (-mp.dps + 40):
                raise ValueError("resonance with logarithm at j=%d (rhs=%s)" % (j, mp.nstr(rhs, 10)))
            c.append(mpf(0))      # harmless resonance: free coefficient, set to 0
            continue
        c.append(rhs / piv)
    return c, d

def evalser(c, d, t):
    tot = mpf(0); x = 1 / (t * t); pw = mpf(1); smallest = None
    for cj in c:
        term = cj * pw
        if smallest is None or abs(term) < smallest:
            smallest = abs(term)
        tot += term; pw *= x
    return tot * t ** d, abs(term)

def limit(s, lead, at, d1, d2, a, N, J):
    c1, e1 = series_coeffs(s, at, d1, J); c2, e2 = series_coeffs(s, at, d2, J)
    fp = [mpf(1), mpf(a(0)) / lead]   # f(n) = u_{n-1}/(lead^n n!^s):  p_{-1}=1, p_0=a_0
    fq = [mpf(0), mpf(1) / lead]
    for n in range(1, N + 1):
        fp.append((a(n) * fp[n] - lead * n**s * fp[n - 1]) / (lead * (n + 1) ** s))
        fq.append((a(n) * fq[n] - lead * n**s * fq[n - 1]) / (lead * (n + 1) ** s))
    t0, t1 = mpf(N) + mpf(1) / 2, mpf(N) + mpf(3) / 2
    F1a, er1 = evalser(c1, e1, t0); F1b, _ = evalser(c1, e1, t1)
    F2a, er2 = evalser(c2, e2, t0); F2b, _ = evalser(c2, e2, t1)
    det = F1a * F2b - F1b * F2a
    al1 = (fp[N] * F2b - fp[N + 1] * F2a) / det
    be1 = (fq[N] * F2b - fq[N + 1] * F2a) / det
    al2 = (F1a * fp[N + 1] - F1b * fp[N]) / det
    be2 = (F1a * fq[N + 1] - F1b * fq[N]) / det
    return al1 / be1, (al1, be1, al2, be2), (er1, er2)

def at5(c3, c1):
    # atil in t for s=5: n^5+(n+1)^5 = 2t^5+5t^3+(5/8)t ; n^3+(n+1)^3 = 2t^3+(3/2)t ; 2n+1 = 2t
    return {5: mpf(2), 3: mpf(5) + 2 * mpf(c3.numerator) / c3.denominator,
            1: mpf(5) / 8 + mpf(3) / 2 * mpf(c3.numerator) / c3.denominator + 2 * mpf(c1.numerator) / c1.denominator}

if __name__ == "__main__":
    mp.dps = 420
    cases = {'U2': (Fr(-3, 2), Fr(7, 8), Fr(-1), Fr(-3)),
             'U3': (Fr(5, 2), Fr(-5, 8), Fr(1), Fr(-5)),
             'Z5a': (Fr(6), Fr(0), Fr(2), Fr(-6))}
    for name, (c3, c1, d1, d2) in cases.items():
        a = lambda n, c3=c3, c1=c1: 8 * (n**5 + (n + 1)**5) + int(8 * c3) * (n**3 + (n + 1)**3) + int(8 * c1) * (2 * n + 1)
        try:
            for (N, J) in ((150, 200), (170, 220)):
                L, coef, errs = limit(5, 8, at5(c3, c1), d1, d2, a, N, J)
                print(name, N, J, mp.nstr(L, 300), 'trunc', mp.nstr(errs[0], 3), mp.nstr(errs[1], 3), flush=True)
        except ValueError as ex:
            print(name, 'ERR', ex, flush=True)
