#!/usr/bin/env python3
"""
t2_closed.py -- closed-form witness at y_0 = 0 for the class a = 1, and central non-vanishing.

For a primitive character chi of conductor f > 1 and parity a, UNCONDITIONALLY
(zeros are symmetric under rho -> 1 - conj(rho); Davenport Ch. 12 (17)-(18) at s = 1 and s = 2):

    c(chi) := -Re sum_rho 1/(rho(rho+1))
            = Re[(L'/L)(2,chi) - (L'/L)(1,chi)] + (1/2)[psi(1 + a/2) - psi((1+a)/2)],

and for zeta:  c(zeta) = log(2 pi) + log 2 - 1/2 - (zeta'/zeta)(-1)   (from psi_1(1) = 0).
Under GRH(q):  F_{q,1}(0) = sum_{chi mod q} c(chi*) =: Phi(q)   (HKMN23.md, Theorem 5).

For every conductor f <= Q (f != 2 mod 4) this script computes, in Arb,
    c_f   = sum over primitive chi mod f of c(chi),
    minL  = min over primitive chi mod f of |L(1/2, chi)|, certified > 0  (central non-vanishing),
and writes one JSON line per conductor.   Usage: python3 t2_closed.py Q shard nshards
"""
import sys, json, math
from flint import acb, arb, acb_series, dirichlet_char, ctx
ctx.prec = 96

def ld(chi, s):
    ser = acb_series([acb(s), acb(1)], prec=2).dirichlet_l(chi)
    return ser[1] / ser[0]

def main():
    Q, shard, nsh = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    dg = {0: (arb(1).digamma() - (arb(1) / 2).digamma()) / 2, 1: ((arb(3) / 2).digamma() - arb(1).digamma()) / 2}
    fs = [f for f in range(3, Q + 1) if f % 4 != 2]
    # balance shards: interleave
    for idx, f in enumerate(fs):
        if idx % nsh != shard:
            continue
        cf = arb(0)
        minL = None
        cnt = 0
        nonvanish = True
        for n in range(1, f):
            if math.gcd(n, f) != 1:
                continue
            nb = pow(n, -1, f)
            if nb < n:
                continue          # conjugate pair handled once
            chi = dirichlet_char(f, n)
            if not chi.is_primitive():
                continue
            mult = 1 if nb == n else 2
            c = (ld(chi, 2) - ld(chi, 1)).real + dg[chi.parity()]
            cf += mult * c
            cnt += mult
            v = abs(chi.l(acb(arb(1) / 2)))
            if not (v > 0):
                nonvanish = False
            lo = float(v.lower())
            if minL is None or lo < minL:
                minL = lo
        print(json.dumps({"f": f, "nprim": cnt, "c_f": [cf.mid().str(30, radius=False), cf.rad().str(5, radius=False)],
                          "minL_half_lower": minL, "central_nonvanishing_certified": nonvanish}))
        sys.stdout.flush()

if __name__ == "__main__":
    main()
