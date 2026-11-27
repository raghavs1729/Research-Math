"""arXiv:2609.01562, Conjecture conj:496:  rho(n) < 1 for all n >= 496, where
     rho(n) = max_{1 <= a < n/2} R_n(a),   R_n(a) = G(n,a)/D(n),
     G(n,a) = [2 - (n-1)/(n-a)] C(n,a) C(n-2,a-1),   D(n) = 2 C(2n-2,n-1) - C(2n-1,n-1) = C(2n-2,n-1)/n.
Simplification used here (checked below against the definition):
     R_n(a) = (u+1) * a * 2(2n-1) * C(n,a)^2 / ( (n-1) * n * C(2n,n) ),   u = n - 2a.
Exact integer test of R_n(a) < 1:   (u+1) a 2(2n-1) C(n,a)^2 < (n-1) n C(2n,n).
Monotonicity (exact): R_n(a-1) < R_n(a)  <=>  a(a-1)(n-2a+3) < (n-a+1)^2 (n-2a+1).
For each n the maximum over 1 <= a < n/2 is therefore attained on [a_min - 1, n/2), a_min = least a >= 2 violating that.
"""
import sys
from fractions import Fraction as Fr
from math import comb
from flint import fmpz
def R_def(n, a): return Fr(2*(n-a) - (n-1), n-a) * comb(n, a) * comb(n-2, a-1) / Fr(comb(2*n-2, n-1), n)
def R_simpl(n, a): return Fr((n-2*a+1) * a * 2*(2*n-1) * comb(n, a)**2, (n-1) * n * comb(2*n, n))
assert all(R_def(n, a) == R_simpl(n, a) for n in range(3, 60) for a in range(1, (n+1)//2) if 2*a < n)
def anomaly(n, C2n, row_top):
    """returns True iff rho(n) > 1 ... exact; also returns whether equality occurs. row_top = C(n, floor((n-1)/2)) as fmpz."""
    amax = (n-1)//2
    # a_min: least a in [2, amax] with a(a-1)(n-2a+3) >= (n-a+1)^2 (n-2a+1)
    a_min = None
    for a in range(2, amax+1):
        if a*(a-1)*(n-2*a+3) >= (n-a+1)**2*(n-2*a+1): a_min = a; break
    lo = 1 if a_min is None else max(1, a_min - 1)
    if a_min is None: lo = amax
    rhs = fmpz((n-1)*n) * C2n
    c = row_top; a = amax; worst = None
    gt = False; eq = False
    while a >= lo:
        u = n - 2*a
        lhs = fmpz((u+1)*a*2*(2*n-1)) * c * c
        if lhs > rhs: gt = True
        if lhs == rhs: eq = True
        c = c * a // (n - a + 1)      # C(n, a-1)
        a -= 1
    return gt, eq
if __name__ == '__main__':
    NMAX = int(sys.argv[1])
    C2n = fmpz(comb(2*8, 8)); out = []
    ctop = fmpz(comb(8, 3))
    res = {}
    n = 8
    while n <= NMAX:
        gt, eq = anomaly(n, C2n, ctop); res[n] = (gt, eq)
        # update to n+1
        C2n = C2n * (2*n+1) * (2*n+2) // ((n+1)*(n+1))
        n += 1
        ctop = fmpz(comb(n, (n-1)//2)) if n < 200 else None
        if ctop is None:
            # C(n, floor((n-1)/2)) from the central binomial of 2*floor(n/2)
            pass
        if n >= 200:
            ctop = fmpz(comb(n, (n-1)//2))
    anom = [n for n in res if res[n][0]]
    print('rho(n) > 1 for all 8 <= n <= 375:', all(res[n][0] for n in range(8, 376)))
    print('n in 376..495 with rho(n) <= 1:', [n for n in range(376, 496) if not res[n][0]])
    print(f'largest n <= {NMAX} with rho(n) > 1:', max(anom))
    print(f'rho(n) = 1 occurs for n in:', [n for n in res if res[n][1]])
    print(f'rho(n) < 1 for all 496 <= n <= {NMAX}:', all((not res[n][0]) and (not res[n][1]) for n in range(496, NMAX+1)))
