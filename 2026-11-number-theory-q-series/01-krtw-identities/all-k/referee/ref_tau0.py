"""Referee boundary check of (B-v): Bailey coefficients of tau phi_n (tau F(R) = q^R F(R+1)) relative to a,
exact at a rational point; n >= 1 should give support {n-1,n,n+1} with sum 0, n = 0 should NOT be finite."""
from flint import fmpq
from ref_common import poch
from ref_eigen import bailey_invert
q, a = fmpq(2, 7), fmpq(3, 11)
phi = lambda e, R: fmpq(0) if R < e else 1 / (poch(q, q, R - e) * poch(a * q, q, R + e))
for n in range(0, 4):
    vals = [q**R * phi(n, R + 1) for R in range(16)]
    kap = bailey_invert(vals, a, q, 15)
    supp = [e for e in range(16) if kap[e] != 0]
    print(f"n={n}: support of tau phi_n in phi_e basis (e<=15): {supp}; coefficient sum over e<=15 = {float(sum(kap)):.3e}")
