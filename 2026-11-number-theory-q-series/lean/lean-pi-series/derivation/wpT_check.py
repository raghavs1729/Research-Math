"""WP-T sanity check (Tails.lean): the bounds on rho used in rho_le_below / rho_ge_above
and the tail bound for a(N,k).  Exact rational arithmetic for rho; floats for sqrt(3)/2."""
from fractions import Fraction as Fr
from math import sqrt, floor, ceil

XS = sqrt(3) / 2

def rho(n, k):
    n, k = Fr(n), Fr(k)
    return ((4*n+1)*(4*n+3)**2*(8*k+4*n+1)*(8*k+4*n+3)*(8*k-4*n+5)*(8*k-4*n+7)
            / (131072*n*(2*n-1)*(4*n-1)*(n-k)**2*(n+k+1)**2))

def sigma(n, k):
    n, k = Fr(n), Fr(k)
    return ((2*k+3)*(n-k-1)**2*(8*k+4*n+1)*(8*k+4*n+3)*(8*k+4*n+5)*(8*k+4*n+7)
            / ((2*k+1)*(n+k+1)**2*(8*k-4*n+9)*(8*k-4*n+11)*(8*k-4*n+13)*(8*k-4*n+15)))

def R(y):
    return (4*y*y-1)/(1-y*y)

ok = True
# (1) explicit bounds proved in Lean
#   below (k <= y n, y = x*-eta, n >= 3):  rho <= (1+4/n)(16R+408/n)^2/16384
#   above (y n <= k < n, y = x*+eta, n>=1): rho >= (16R)^2/(1+1/n)^2/16384
for eta in (0.001, 0.01, 0.05, 0.0999):
    yb, ya = XS - eta, XS + eta
    Rb, Ra = R(yb), R(ya)
    assert 1 <= Rb < 8 < Ra
    worst_b, worst_a = 0.0, 1e9
    for n in list(range(3, 400)) + [1000, 5000]:
        ub = (1+4/n)*(16*Rb+408/n)**2/16384
        lb = (16*Ra)**2/(1+1/n)**2/16384
        for k in range(0, floor(yb*n)+1):
            r = float(rho(n, k))
            if not (0 < r <= ub*(1+1e-12)):
                ok = False; print("BELOW FAIL", eta, n, k, r, ub)
            if n >= 2000: worst_b = max(worst_b, r)
        for k in range(ceil(ya*n), n):
            r = float(rho(n, k))
            if not (r >= lb*(1-1e-12)):
                ok = False; print("ABOVE FAIL", eta, n, k, r, lb)
            if n >= 2000: worst_a = min(worst_a, r)
    print(f"eta={eta}: rhobar(y-)={Rb*Rb/64:.5f} sup rho(n=5000)={worst_b:.5f};"
          f" rhobar(y+)={Ra*Ra/64:.5f} inf rho(n=5000)={worst_a:.5f}")

# (2) tail_bound for f = a (the WZ member): a(N,k) outside |k - x*N| <= delta N decays geometrically
def a_row(N):
    A = Fr(9, 256)
    for m in range(1, N):
        A *= rho(m, 0)
    row = [A]
    for k in range(N-1):
        row.append(row[-1]*sigma(N, k))
    return row
for delta in (0.02, 0.05, 0.09):
    prev = None
    for N in (50, 100, 200, 400):
        row = a_row(N)
        assert all(v > 0 for v in row)
        mx = max(float(row[k]) for k in range(N) if abs(k - XS*N) >= delta*N)
        root = mx ** (1.0/N)
        print(f"delta={delta} N={N}: max tail a = {mx:.3e}, N-th root = {root:.5f}")
print("ALL OK" if ok else "FAILURES")
