#!/usr/bin/env python3
"""
a1_character_zeros.py -- audit item 1 (independent of the audited code; nothing
is imported from ../nonmonotone).

 (a) chi = flint.dirichlet_char(q, l) used by the audited code (l = 162 for
     q = 163, l = 115 for q = 232): compare chi(n) (evaluated by flint as a
     complex number) with my own Kronecker symbol (-q/n) for 1 <= n <= 400.
     Also: L(2, chi) from Arb vs a direct mpmath sum with my Kronecker values,
     and L(1/2+it) at a few t vs mpmath's Hurwitz-based dirichlet().
 (b) semantics of acb.lgamma (continuous log-gamma, not log of Gamma).
 (c) independent float location of the first 20 ordinates with mpmath (own Z(t)),
     compared with zeros_q<q>.txt.
 (d) Riemann-von Mangoldt main term at T_cert vs the count K.
 (e) every listed window: certified sign change of Arb's Hardy Z at d -+ 1e-20,
     re-done here (own loop), and windows disjoint, inside (0, T_cert).
"""
import sys, math
from flint import arb, acb, ctx, dirichlet_char
import mpmath as mp

HERE = "/home/user/weighted-prime-races/monotonicity/nonmonotone"


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def kron(D, n):
    # Kronecker (D/n), D = discriminant, n >= 1
    r = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0:
            return 0
        r *= 1 if D % 8 in (1, 7) else -1
    return r * jacobi(D, n)


def main(q, l, nlist=20):
    D = -q
    ctx.prec = 128
    chi = dirichlet_char(q, l)
    print("q=%d l=%d: is_real=%s primitive=%s conductor=%d parity=%d" %
          (q, l, chi.is_real(), chi.is_primitive(), chi.conductor(), chi.parity()))
    bad = 0
    for n in range(1, 401):
        v = chi(n)                       # acb
        k = kron(D, n)
        if not acb(v).overlaps(acb(k)) or (abs(acb(v) - k).upper() > 1e-30):
            bad += 1
            print("  MISMATCH n=%d chi=%s kron=%d" % (n, v, k))
    print("  chi(n) == Kronecker(%d/n) for 1<=n<=400: %s" % (D, "YES" if bad == 0 else "NO (%d)" % bad))
    # chi(-1): parity
    print("  chi(q-1) =", chi(q - 1))
    # (a') L values via my own character list
    mp.mp.dps = 30
    clist = [kron(D, n) if n > 0 else 0 for n in range(q)]
    clist[0] = 0
    L2_mp = mp.dirichlet(2, clist)
    L2_arb = chi.l(acb(2))
    print("  L(2): arb %s   mpmath(own chi) %s   diff %.2e" % (L2_arb.real.str(20), mp.nstr(L2_mp, 20),
                                                              abs(float(L2_arb.real.mid()) - float(L2_mp))))
    for t in (0.0, 3.7, 15.25):
        a = chi.l(acb(arb(1) / 2, arb(t)))
        m = mp.dirichlet(mp.mpc(0.5, t), clist)
        print("  L(1/2+%gi): arb %s ; mpmath %s" % (t, a.str(15, radius=False), mp.nstr(m, 15)))
    # (b) lgamma semantics
    z = acb(arb(3) / 4, 60)
    lg = z.lgamma().imag
    argG = z.gamma().arg()
    print("  Im lgamma(3/4+60i) = %s ; Arg Gamma(3/4+60i) = %s  (continuous branch iff they differ by 2pi k, k != 0 here)"
          % (lg.str(12), argG.str(12)))
    ref = mp.im(mp.loggamma(mp.mpc(0.75, 60)))
    print("  mpmath Im loggamma(3/4+60i) = %s" % mp.nstr(ref, 15))
    # (c) independent zero location
    def Z(t):
        th = (t / 2) * mp.log(q / mp.pi) + mp.im(mp.loggamma(mp.mpc(0.75, t / 2)))
        return mp.re(mp.exp(1j * th) * mp.dirichlet(mp.mpc(0.5, t), clist))
    mp.mp.dps = 25
    found = []
    t, step = 0.0, 0.02
    zp = Z(mp.mpf(0))
    while len(found) < nlist:
        t2 = t + step
        z2 = Z(mp.mpf(t2))
        if z2 * zp < 0:
            found.append(mp.findroot(Z, (mp.mpf(t), mp.mpf(t2)), solver="anderson"))
        t, zp = t2, z2
    listed = [l.split()[0] for l in open("%s/zeros_q%d.txt" % (HERE, q)) if l.strip() and not l.startswith("#")]
    worst = 0
    for i, g in enumerate(found):
        d = abs(g - mp.mpf(listed[i]))
        worst = max(worst, d)
    print("  first %d ordinates, mpmath (own Z, step-0.02 scan) vs file: max |diff| = %s" % (nlist, mp.nstr(worst, 3)))
    print("   gamma_1..3 (mpmath):", [mp.nstr(g, 18) for g in found[:3]])
    # (d) RvM main term: N(T) ~ (T/2pi) log(qT/(2 pi e)) + chi-independent const (odd: -1/8 ... )
    Ttop = {163: 7617 / 64, 232: 6351 / 64}[q]
    main = Ttop / (2 * math.pi) * math.log(q * Ttop / (2 * math.pi * math.e))
    # theta(T)/pi is the smooth count incl. constants
    th = (Ttop / 2) * math.log(q / math.pi) + float(mp.im(mp.loggamma(mp.mpc(0.75, Ttop / 2))))
    print("  RvM at T_cert=%.6f: (T/2pi)log(qT/2pi e) = %.3f ; theta(T)/pi = %.3f ; listed K = %d"
          % (Ttop, main, th / math.pi, len(listed)))
    # (e) certified sign changes, own loop
    ctx.prec = 192
    eps = arb("1e-20")
    prev_hi = arb(0)
    ok = 0
    for s in listed:
        d = arb(s)
        a, b = d - eps, d + eps
        za = chi.hardy_z(acb(a)).real
        zb = chi.hardy_z(acb(b)).real
        if ((za > 0 and zb < 0) or (za < 0 and zb > 0)) and a > prev_hi:
            ok += 1
        else:
            print("  WINDOW FAIL", s)
        prev_hi = b
    print("  certified sign changes in disjoint increasing windows: %d / %d ; last window < T_cert: %s"
          % (ok, len(listed), prev_hi < arb(Ttop)))


if __name__ == "__main__":
    main(163, 162)
    main(232, 115)
