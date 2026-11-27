"""AUDIT-4 part 2 (own code): rigorous (arb ball arithmetic, 128 bits) tests of every displayed inequality of R10
for all 3 <= n <= NMAX and ALL 1 <= a < n/2.  A check 'x < y' is counted as passed only if the ball x-y is certainly negative.
"""
import sys
from flint import arb, fmpz, ctx
ctx.prec = 128
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
pi = arb.pi(); e1 = arb(1).exp()
def lt(x, y): return bool(x < y)
def le(x, y): return bool(x <= y)

# ---- (B1) for 1 <= m <= 2*NMAX+2 ; monotonicity of w_m^2(m+1/4), w_m^2(m+1/2)
c = fmpz(1); prev_lo = None; prev_hi = None; okB1 = True; margins = []
for m in range(0, 2*NMAX+3):
    if m > 0: c = c*(2*m)*(2*m-1)//(m*m)
    w = arb(c)/arb(4)**m
    lo = w*w*(arb(m)+arb(1)/4); hi = w*w*(arb(m)+arb(1)/2)
    if m >= 1:
        if not (lt(lo, 1/pi) and lt(1/pi, hi)): okB1 = False; print("B1 FAIL", m)
        if not (lt(prev_lo, lo) and lt(hi, prev_hi)): okB1 = False; print("B1 monotone FAIL", m)
    prev_lo, prev_hi = lo, hi
print("(B1) strict two-sided bound and monotonicity for 1<=m<=%d: %s ; pi*w^2(m+1/4), pi*w^2(m+1/2) at m=%d: %s %s"
      % (2*NMAX+2, okB1, 2*NMAX+2, (pi*prev_lo).str(12), (pi*prev_hi).str(12)))

# ---- per (n,a) checks
fails = {k: [] for k in ("B2/B3", "(2)", "pref<2", "(3)", "chain1", "chain2", "phi_le", "phi_final", "final", "final_n<12")}
worst_ratio = arb(0); worst_n = None       # max over n>=12 of rho(n)/bound
cnt = 0
C2n = fmpz(1)
import math
for n in range(1, NMAX+1):
    C2n = C2n*(2*n)*(2*n-1)//(n*n)          # C(2n,n)
    if n < 3: continue
    amax = (n-1)//2
    cen = fmpz(math.comb(n, n//2))
    sq = arb(n).sqrt()
    bound = 4/(2*pi*e1).sqrt()*(arb("4.128")/sq).exp()
    phi_cap = (arb(n)/(2*e1)).sqrt()*(arb("4.128")/sq).exp()
    pref2 = 2/(pi*(arb(n)+arb(1)/2)).sqrt()
    c = fmpz(n); rho = arb(0)
    for a in range(1, amax+1):
        if a > 1: c = c*(n-a+1)//a
        u = n-2*a; v = u+1; cnt += 1
        assert 1 <= u <= n-2 and v >= 2
        R = arb(fmpz((u+1)*a*2*(2*n-1))*c*c)/arb(fmpz((n-1)*n)*C2n)
        rho = rho.max(R) if hasattr(rho, "max") else (R if bool(R > rho) else rho)
        # (B2)/(B3): C(n,a)/C(n,floor(n/2)) <= exp(-E/2), E=(u^2 or u^2-1)/(n+u)
        E = arb(u*u if n % 2 == 0 else u*u-1)/(n+u)
        if not le(arb(c)/arb(cen), (-E/2).exp()): fails["B2/B3"].append((n, a))
        # (2)
        E2 = arb(u*u-1)/(n+u)
        if not lt(arb(c*c)/arb(C2n), pref2*(-E2).exp()): fails["(2)"].append((n, a))
        # prefactor
        if not (n-u)*(2*n-1) < 2*n*(n-1): fails["pref<2"].append((n, a))
        phi = v*(-E2).exp()
        if not lt(R, 2*pref2*phi): fails["(3)"].append((n, a))
        # chain: (u^2-1)/(n+u) >= v^2/(n+v) - 2v/n >= v^2/n - v^3/n^2 - 2v/n   (exact rationals)
        from fractions import Fraction as F
        e0 = F(u*u-1, n+u); e1_ = F(v*v, n+v)-F(2*v, n); e2_ = F(v*v, n)-F(v**3, n*n)-F(2*v, n)
        if not (e0 == F(v*v-2*v, n+v-1) and e0 >= e1_): fails["chain1"].append((n, a))
        if not e1_ >= e2_: fails["chain2"].append((n, a))
        if not le(phi, v*(-(arb(v*v)/n)).exp()*(arb(v**3)/(n*n)+arb(2*v)/n).exp()): fails["phi_le"].append((n, a))
        if n >= 12 and not le(phi, phi_cap): fails["phi_final"].append((n, a))
    if n >= 12:
        if not lt(rho, bound): fails["final"].append(n)
        q = rho/bound
        if bool(q > worst_ratio): worst_ratio = q; worst_n = n
    else:
        if not lt(rho, bound): fails["final_n<12"].append((n, rho.str(8), bound.str(8)))
    if n in (12, 100, 496, 1000, 2000, 3000, NMAX):
        print("  n=%5d  rho(n)=%s   R10 bound=%s" % (n, rho.str(15), bound.str(15)))
print("pairs (n,a) tested:", cnt, " for 3<=n<=%d" % NMAX)
for k, v_ in fails.items(): print("  failures of %-11s: %d %s" % (k, len(v_), v_[:8]))
print("max over 12<=n<=%d of rho(n)/bound = %s at n=%s" % (NMAX, worst_ratio.str(12), worst_n))
