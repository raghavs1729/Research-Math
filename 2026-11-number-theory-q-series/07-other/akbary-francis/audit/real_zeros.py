#!/usr/bin/env python3
"""real_zeros.py q -- rigorous (Arb ball) check that L(sigma, chi) != 0 for sigma in [1/2, 1] for every
non-principal chi mod prime q (own characters, Hurwitz zeta).  With the functional equation this
gives no real zero in (0, 1).  This is the hypothesis 'no singularity of script-L(s; q, 1) on (1-b, 1)'
of Akbary-Francis Thm 1.16(a) (zeta has no zero in (0,1): classical).
[1/2, 1-h]: subdivision into balls, Re or Im of the ball enclosure excludes 0.
[1-h, 1]: |L(s) - L(1)| <= h * max|L'| <= h * M/(r - h), M = max |L| on |z - 1| = r (maximum
principle; L entire), M bounded by balls covering the circle; L(1) = -(1/q) sum chi(m) digamma(m/q)."""
import sys
from flint import arb, acb, ctx
ctx.prec = 128
q = int(sys.argv[1]); phi = q - 1
g = next(g for g in range(2, q) if len({pow(g, k, q) for k in range(phi)}) == phi)
ind = {pow(g, k, q): k for k in range(phi)}


def L(s, j):
    tot = acb(0)
    for m in range(1, q):
        tot += (2 * acb.pi() * acb(0, 1) * (j * ind[m] % phi) / phi).exp() * acb.zeta(s, acb(m) / q)
    return acb(q) ** (-s) * tot


def nz(z):
    return not (z.real.contains(0) and z.imag.contains(0))


allok = True
h = arb(sys.argv[2] if len(sys.argv) > 2 else "0.05"); r = arb(sys.argv[3] if len(sys.argv) > 3 else "0.45")
for j in range(1, phi):
    ok = True; n = 0
    a = arb("0.5"); b = 1 - h
    N = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    step = (b - a) / N
    for i in range(N):
        s = acb(arb((a + step * (i + arb("0.5"))).mid(), (step / 2 * arb("1.0001")).mid()))
        if not nz(L(s, j)):
            ok = False; print("j=%d: cannot exclude zero near %s" % (j, s))
    L1 = -sum((2 * acb.pi() * acb(0, 1) * (j * ind[m] % phi) / phi).exp() * acb(arb(m) / q).digamma()
              for m in range(1, q)) / q
    # max |L| on circle |z-1| = r, 256 balls
    M = arb(0); K = 256
    for k in range(K):
        t = 2 * arb.pi() * (k + arb("0.5")) / K
        c = acb(1) + r * acb(t.cos(), t.sin())
        z = acb(c, None) if False else acb(arb(c.real.mid(), (r * arb.pi() / K * 1.01).mid()),
                                         arb(c.imag.mid(), (r * arb.pi() / K * 1.01).mid()))
        v = abs(L(z, j)).upper()
        M = M if M > v else arb(v)
    bound = h * M / (r - h)
    near1 = bool(abs(L1).lower() > bound.upper())
    ok = ok and near1
    allok = allok and ok
    print("j=%2d parity=%s  |L(1)| = %s  near-1 bound %s  [1/2,1]: %s" % (
        j, "odd" if j % 2 else "even", abs(L1).str(6), bound.str(4), "no zero" if ok else "NOT EXCLUDED"))
print("ALL: no real zero in [1/2,1]" if allok else "NOT ALL EXCLUDED")
