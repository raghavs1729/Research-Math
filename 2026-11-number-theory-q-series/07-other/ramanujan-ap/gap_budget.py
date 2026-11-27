#!/usr/bin/env python3
"""What a zero list to height H buys in the gap (ESTIMATE, not a theorem; q = 3 or 4).

Under GRH, by monotonicity of psi(x;q,a) and the absolutely convergent explicit formula for
psi_1 = int psi, for every h = eta*sqrt(x) > 0
   Psi(x) := (phi psi(x;q,a) - x)/sqrt(x)  <=  eta/2 + sum_chi | sum_rho ((x+h)^{rho+1} - x^{rho+1})/(h rho (rho+1)) | / sqrt(x) + o(1)
(and symmetrically from below).  A zero pair at height gamma contributes at most
   min( 2/gamma , 4 sqrt(x)/(eta gamma^2) )      (crossover at T_c = 2 sqrt(x)/eta).
With the zeros known up to H (their sum S_low evaluated exactly; we allow S_low = 2.5 in absolute
value, its observed size) and the others taken in absolute value with density
(1/2pi) log(q gamma/2pi) per unit height and per character:
   B(x; H, eta) = eta/2 + S_low + sum_chi [ (l(T_c)^2 - l(H)^2)^+ / (2 pi) + (l(max(T_c,H)) + 1)/pi * min(1, T_c/H) ],
   l(T) = log(cond * T/(2 pi)).
No failure at x is guaranteed (for a class with c = 0; up to the O(1/log x) terms) when
   B(x) and B(x/e) <= W(x)/(2 lambda + sqrt(e) mu) =: allowed(x).
The first term eta/2 is the price of the window (it is the expected number of primes in it) and
does not depend on the kernel: eta <= 2*allowed forces T_c >= sqrt(x)/allowed(x) whatever the kernel.
Output: allowed(x), the kernel-independent minimal height sqrt(x)/allowed, and for several H the
best B over eta (crude kernel)."""
import math
from model import W_of
def lam(x):
    from mpmath import li, log
    return float(log(x) * li(x) / x)
def B(x, H, eta, conds=(1, 4), slow=2.5):
    Tc = 2 * math.sqrt(x) / eta
    tot = eta / 2 + (slow if H > 0 else 0.0)
    for cond in conds:
        l = lambda T: math.log(max(cond * T / (2 * math.pi), 1.0001))
        Hc = max(H, 2 * math.pi / cond * 1.0001)
        if Tc > Hc:
            tot += (l(Tc) ** 2 - l(Hc) ** 2) / (2 * math.pi) + (l(Tc) + 1) / math.pi
        else:   # all zeros below T_c known; unknown zeros lie above H >= T_c
            tot += (l(Hc) + 1) / math.pi * (Tc / Hc)
    return tot
print("x        W(x)   allowed   sqrt(x)/allowed |  best B over eta (crude kernel) for H = 640, 1e5, 1e6, 1e7, 1e8, 1e9")
for x in (1e14, 1e15, 1e16, 1e17, 3e17, 1e18, 2e18):
    L = math.log(x)
    allowed = W_of(x) / (2 * lam(x) + math.sqrt(math.e) * L / (L - 1))
    row = []
    for H in (640, 1e5, 1e6, 1e7, 1e8, 1e9):
        best = min((B(x, H, 10 ** (k / 50.0)), 10 ** (k / 50.0)) for k in range(-100, 200))
        row.append("%7.1f%s" % (best[0], "*" if best[0] <= allowed else " "))
    print("%-7.0e %7.1f %7.1f   %9.2e      | %s" % (x, W_of(x), allowed, math.sqrt(x) / allowed, " ".join(row)))
print("(* = within the allowed budget)")
print()
print("cost of certified zeros of L(s,chi_4) with python-flint/Arb on this machine (measured Hardy-Z times:")
print("0.07 s at t = 1e5, 0.66 s at t = 1e6, proportional to t), at 6 evaluations per zero:")
for H in (1e5, 1e6, 1e7, 1e8):
    N = H / math.pi * math.log(4 * H / (2 * math.pi * math.e))
    cost = N * 6 * 0.66e-6 * H / 2          # average height H/2
    print("   H = %.0e: %.2e zeros, about %.1e core-seconds = %.1f core-days" % (H, N, cost, cost / 86400))
