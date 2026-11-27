#!/usr/bin/env python3
"""
a3_indep_density.py -- independent rigorous enclosure of
    delta(M) = P(1 + sum_{gamma>0} c_gamma cos theta_gamma > 0),  c_gamma = 4M/sqrt(M^2+gamma^2),
written for the audit without using density.py.  Only the zero WINDOWS are taken
from ../nonmonotone/zeros_q<q>.txt (first N lines; sign changes re-verified in
a1_character_zeros.py), and COMPLETENESS below T0 is taken from MY OWN count
a2_count.py (W = 2N at T0, run separately; T0 passed on the command line).

Ingredients (re-derived in AUDIT.md):
  * Gil-Pelaez: delta = 1/2 + (1/pi) int_0^inf Phi(t) sin t/t dt.
  * sigma^2(M) = 4M[ (1/2)log(q/pi) + (1/2)psi((M+3/2)/2) + (L'/L)(M+1/2) ]  (GRH),
    L'/L from Arb's power-series jet ONLY (no complex step).
  * V = sigma^2 - sum_{j<=N} b_j^2/2 = sum over unlisted zeros (mult.) of c^2/2;
    all unlisted gamma >= T0, so c_gamma <= c_top := 4M/sqrt(M^2+T0^2).
  * Tail model (my form): for t <= T <= 1/c_top,
        g(t) e^{-C t^4 S4} <= Psi(t) <= g(t),  g = e^{-V t^2/2},  S4 = sum_unl c^4 <= 2 V c_top^2,
        C = -log J0(1) - 1/4 (enclosed in Arb);  hence  |Psi - g| <= g * C S4 t^4.
    Main term = acb.integral of P(t) g(t) sinc(t) on [0,T] (Arb's rigorous
    adaptive Gauss-Legendre integrator with analytic error bounds; integrand entire).
    Model error <= C*S4 * int_0^T |P| g t^3 |sin t| dt, bounded by an upper Riemann
    sum of Arb interval enclosures on cells of width h.
  * Tail in t: |J0(x)| <= sqrt(2/(pi x)) (x>0), so for any K,
        |int_T^inf Phi sinc| <= prod_{j<=K} (2/(pi b_j))^{1/2} T^{-K/2} / (K/2).
Usage: python3 a3_indep_density.py q l N T0num/T0den T M1 M2 ...
"""
import sys, time
from fractions import Fraction
from flint import arb, acb, acb_series, ctx, dirichlet_char

ctx.prec = 128
HERE = "/home/user/weighted-prime-races/monotonicity/nonmonotone"


def frac(s):
    f = Fraction(s)
    return arb(f.numerator) / f.denominator


def run(q, l, N, T0, T, M, cell=arb(1) / 64):
    t_start = time.time()
    chi = dirichlet_char(q, l)
    lines = [x.split()[0] for x in open("%s/zeros_q%d.txt" % (HERE, q)) if x.strip() and not x.startswith("#")]
    gam = [arb(s, "1e-20") for s in lines[:N]]
    assert len(gam) == N and gam[-1].upper() < T0
    b = [4 * M / (M * M + g * g).sqrt() for g in gam]
    # sigma^2 via Arb jet
    s = M + arb(1) / 2
    jet = acb_series.dirichlet_l(acb_series([acb(s), acb(1)], prec=2), chi)
    L0, L1 = jet[0].real, jet[1].real
    assert L0 > 0
    sig2 = 4 * M * ((arb(q) / arb.pi()).log() / 2 + ((M + arb(3) / 2) / 2).digamma() / 2 + L1 / L0)
    V = sig2 - sum((x * x / 2 for x in b), arb(0))
    assert V > 0, "V not positive"
    ctop = 4 * M / (M * M + T0 * T0).sqrt()
    if T is None:      # auto: dyadic T <= min(0.95/c_top, 40)
        T = arb(int(min(0.95 / float(ctop.upper()), 40.0) * 64)) / 64
    assert (T * ctop).upper() <= 1, "T beyond tail-model range"
    C = -(arb(1).bessel_j(0)).log() - arb(1) / 4
    S4 = 2 * V * ctop * ctop
    # main integral
    def f(z, analytic):
        p = acb(1)
        for x in b:
            p *= (x * z).bessel_j(0)
        return p * (-V * z * z / 2).exp() * z.sinc()
    I = acb.integral(f, 0, T, rel_tol=arb(2) ** -50, abs_tol=arb(2) ** -50, eval_limit=10 ** 7)
    Ire = I.real
    assert abs(I.imag).upper() < 1e-10
    # model error: upper Riemann sum
    n = int((T / cell).upper().ceil().unique_fmpz())
    emod = arb(0)
    for k in range(n):
        a, c = cell * k, cell * (k + 1)
        tt = a.union(c)
        p = arb(1)
        for x in b:
            p *= (x * tt).bessel_j(0)
        val = abs(p) * (-V * tt * tt / 2).exp() * tt ** 3 * abs(tt.sin())
        emod += val.upper() * (c - a)
    # (the last cell may overshoot T: harmless, integrand >= 0)
    emod = C * S4 * emod
    # tail in t
    best = None
    A = arb(1)
    for k, x in enumerate(b, 1):
        A = A * (2 / (arb.pi() * x * T)).sqrt()
        if k >= 2:
            e2 = A / (arb(k) / 2)
            if best is None or e2.upper() < best.upper():
                best, K = e2, k
    rad = emod.upper() + best.upper() + Ire.rad()
    dens = arb(1) / 2 + Ire / arb.pi()
    lo = (dens - rad / arb.pi()).lower()
    hi = (dens + rad / arb.pi()).upper()
    print("q=%d M=%s N=%d T0=%.4f T=%s: delta in [%s, %s]  width %.2e | sigma2=%s V=%s c_top=%.4e "
          "E_model=%.2e E_tail=%.2e (K=%d) rad(I)=%.1e  %.0fs"
          % (q, M.str(6, radius=False), N, float(T0.mid()), T.str(6, radius=False),
             lo.str(14, radius=False), hi.str(14, radius=False), float((hi - lo).mid()),
             sig2.str(15, radius=False), V.str(10, radius=False), float(ctop.mid()),
             float(emod.upper()), float(best.upper()), K, float(Ire.rad()), time.time() - t_start), flush=True)
    return lo, hi


if __name__ == "__main__":
    q, l, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    T0 = frac(sys.argv[4]); T = None if sys.argv[5] == 'auto' else frac(sys.argv[5])
    res = {}
    for Ms in sys.argv[6:]:
        res[Ms] = run(q, l, N, T0, T, frac(Ms))
    ks = list(res)
    for i in range(len(ks) - 1):
        a, b_ = res[ks[i]], res[ks[i + 1]]
        print("  upper(delta(%s)) < lower(delta(%s)) : %s ; lower(delta(%s)) > upper(delta(%s)) : %s"
              % (ks[i], ks[i + 1], a[1] < b_[0], ks[i], ks[i + 1], a[0] > b_[1]))
