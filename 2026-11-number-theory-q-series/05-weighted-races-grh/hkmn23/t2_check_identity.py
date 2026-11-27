#!/usr/bin/env python3
"""
t2_check_identity.py -- numerical check (not part of a proof) of the identity used in Theorem 5(b):

   T_1 := - sum_{chi != chi_0 mod q} Re (L'/L)(1, chi*)  =  1 + gamma + I + sum_{p | q} phi(m_p) log p / (p^{o_p} - 1),
   I   := 1/2 + 2 int_1^oo (phi psi_1(x;q,1) - x^2/2) x^{-3} dx,
   m_p = q / p^{v_p(q)},  o_p = multiplicative order of p mod m_p (o_p = 1 if m_p = 1).

I is approximated by I_X = 1/2 + phi * sum_{n <= X, n = 1 (q)} Lambda(n) (1/n - 2/X + n/X^2) - log X  (exact value of the
truncated integral); the truncation error is O(L_1(q) X^{-1/2}) under GRH.
"""
import math, sys
import numpy as np
from flint import arb, acb, acb_series, dirichlet_char, ctx
import t2_witness as TW
ctx.prec = 80
X = 10 ** 7
pp = TW.PP(X)
def ld1(f, m):
    ser = acb_series([acb(1), acb(1)], prec=2).dirichlet_l(dirichlet_char(f, m))
    return float((ser[1] / ser[0]).real.mid())
for q in (3, 4, 5, 7, 12, 15, 16, 21, 40, 60, 97):
    phi = TW.phi_of(q)
    T1 = -sum(ld1(f, m) for chi, f, m, pr in TW.char_table(q) if not pr)
    mask = (pp.n % q) == 1
    n = pp.n[mask].astype(np.float64); lam = pp.lam[mask]
    IX = 0.5 + phi * float(np.sum(lam * (1 / n - 2 / X + n / X ** 2))) - math.log(X)
    corr = 0.0
    qq, p = q, 2
    while qq > 1:
        if qq % p == 0:
            pe = 1
            while qq % p == 0:
                qq //= p; pe *= p
            m = q // pe
            o = 1
            if m > 1:
                while pow(p, o, m) != 1:
                    o += 1
            corr += TW.phi_of(m) * math.log(p) / (p ** o - 1)
        p += 1
    rhs = 1 + 0.5772156649015329 + IX + corr
    print(f"q={q:3d}: T_1 = {T1:+.6f}   1+gamma+I_X+corr = {rhs:+.6f}   diff = {T1 - rhs:+.2e}   (I_X = {IX:+.6f}, corr = {corr:.6f}; allowed ~ 4 L_1 X^-1/2 ~ {4*phi*math.log(q+2)/X**0.5:.1e})")
