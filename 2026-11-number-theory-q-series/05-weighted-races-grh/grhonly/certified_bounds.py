#!/usr/bin/env python3
"""
certified_bounds.py -- GRH-only (no LI) lower bounds for the logarithmic
density of the lead sets {x : D_m(x) > 0}, q = 4 and q = 3.

Companion to monotonicity/GRH-ONLY.md (Section 4: Lemmas 4.2, 4.3, Propositions 4.4,
4.5; Section 5: Theorems 5.1-5.3).  Ball arithmetic
(python-flint / Arb) throughout: every printed enclosure "[mid +/- rad]" is a
rigorous enclosure GIVEN the inputs listed under "Inputs" below.

What is computed, for chi = chi_4, chi_3 and a list of weights M = m + 1/2:

  sigma_M^2 = 4M (Lambda'/Lambda)(M + 1/2)        (closed form; = sum over zeros
                                                   with multiplicity m_gamma of
                                                   8M^2/(M^2+gamma^2), under GRH)
  tau_M(T)  = sigma_M^2 - sum_{gamma <= T} 8M^2/(M^2+gamma^2)   (exact tail,
                                                   multiplicity-weighted)
  eps_M(T)  >= sum_{gamma > T} (m_gamma - 1) * 8M^2/(M^2+gamma^2)
               (upper bound for the EXCESS of the true variance V_M of nu_M,
                which weights the zeros by m_gamma^2, over sigma_M^2)
  V_M^+     = sigma_M^2 + eps_M(T)  >=  V_M = Var(nu_M)
  Cantelli  : lower log-density{D_m > 0} >= 1/(1 + V_M^+)
  soft-threshold refinement (GRH-ONLY.md, Theorem 5.3)

Inputs (trusted, not re-proved here):
  (I1) prototype/explore/chi{4,3}_zeros.txt: every nontrivial zero of
       L(s, chi) with 0 < gamma <= T is simple, lies on the critical line, and
       is within 1e-20 of exactly one listed ordinate; T = (last listed
       ordinate) + 0.25, i.e. T_4 = 640.165707..., T_3 = 700.024747...
       (certified by prototype/explore/certify_arb.py --K 0, parts [C],[D]).
  (I2) GRH(chi) for the identities "zero sum = closed form" and for the
       multiplicity bound (Lemma 4.2 of GRH-ONLY.md).

The multiplicity bound used (GRH-ONLY.md Lemma 4.2): for every ordinate gamma
and every delta > 1/2, with x_0 = (1/2 + delta + a)/2 (a = 1 for odd chi),

  m_gamma <= mu_delta(gamma)
          := (delta/4) log(gamma^2 + 4 x_0^2) + delta * b_delta,
  b_delta := (1/2) log(q/(2 pi)) + 1/(24 x_0^2) + |zeta'/zeta|(1/2 + delta).

The log-weighted tail is evaluated EXACTLY through the closed form
(GRH-ONLY.md Lemma 4.3)

  J_A := sum_{gamma>0} gamma^{-2} log(1 + gamma^2/A^2)
       = int_A^infty (Lambda'/Lambda)(M + 1/2) M^{-2} dM
      <= J_A^+ (explicit closed form, see j_upper below).

Runtime: about 20 s on one core.
"""
import os
import sys
import time

from flint import arb, acb, ctx, dirichlet_char, arb_series

ctx.prec = 220
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")
EPS_ZERO = arb("1e-20")          # enclosure radius of each listed ordinate (I1)


def load_zeros(fn):
    out = []
    with open(os.path.join(ZDIR, fn)) as f:
        for line in f:
            line = line.strip()
            if line:
                mid = arb(line)
                out.append(mid + arb(0, 1) * EPS_ZERO)   # mid +/- 1e-20
    return out


def ball_pm(x):
    """the ball [-|x|, |x|]"""
    return arb(0, 1) * abs(x)


def l_derivative(chi, s0, k, radius):
    """Rigorous ball for L^(k)(s0, chi) by the Cauchy integral on
    |s - s0| = radius (L entire); same routine as certify_arb.py."""
    r = arb(radius)

    def integrand(th, analytic):
        z = (acb(0, 1) * th).exp()
        return chi.l(s0 + acb(r) * z) * z ** (-k)

    I = acb.integral(integrand, 0, 2 * arb.pi())
    fact = 1
    for j in range(2, k + 1):
        fact *= j
    return I * arb(fact) / (2 * arb.pi() * r ** k)


def cauchy_radius(s0f):
    """a radius whose circle stays away from s = 1 (Arb's Hurwitz
    representation of L is singular exactly there)"""
    for r in ("0.3", "0.2", "0.25", "0.15"):
        if abs(abs(s0f - 1.0) - float(r)) > 0.04:
            return arb(r)
    raise RuntimeError("no admissible radius")


def log_deriv_Lambda(chi, q, a, s0):
    """(Lambda'/Lambda)(s0, chi) for real s0 > 1/2, s0 an arb."""
    r = cauchy_radius(float(s0.mid()))
    L0 = l_derivative(chi, acb(s0), 0, r)
    L1 = l_derivative(chi, acb(s0), 1, r)
    lpl = (L1 / L0).real
    return (arb(q) / arb.pi()).log() / 2 + ((s0 + a) / 2).digamma() / 2 + lpl


def C_const(chi, a):
    """C(chi) = sum_{gamma>0} 2/gamma^2 = (log Lambda)''(1/2)
            = (1/4) psi'((1/2+a)/2) + (log L)''(1/2)."""
    half = acb(1) / 2
    vals = []
    for rad in ("0.25", "0.3"):
        L0 = l_derivative(chi, half, 0, arb(rad))
        L1 = l_derivative(chi, half, 1, arb(rad))
        L2 = l_derivative(chi, half, 2, arb(rad))
        vals.append((L2 / L0 - (L1 / L0) ** 2).real)
        if not L0.real > 0:
            raise RuntimeError("L(1/2) not certified positive")
    if not vals[0].overlaps(vals[1]):
        raise RuntimeError("Cauchy radii disagree")
    psi1 = acb((arb(1) / 2 + a) / 2).polygamma(1).real
    return psi1 / 4 + vals[0]


def zeta_logderiv_abs(sigma):
    """|zeta'/zeta|(sigma) for real sigma > 1 (arb)."""
    s = arb_series([sigma, 1], prec=2)
    z = s.zeta()
    return abs(z[1] / z[0])


def j_upper(q, a, A, Zval):
    """Closed-form upper bound J_A^+ >= J_A (Lemma 4.3):
       J_A <= log(q/(2 pi))/(2A)
              + (1/2)[ log(A+c)/A + log(1+c/A)/c ]
              - (1/2)[ 1/(cA) - log(1+c/A)/c^2 ]
              + Z(A+1/2)/A ,          c = 1/2 + a."""
    c = arb(1) / 2 + a
    A = arb(A)
    t1 = (arb(q) / (2 * arb.pi())).log() / (2 * A)
    t2 = ((A + c).log() / A + (1 + c / A).log() / c) / 2
    t3 = (1 / (c * A) - (1 + c / A).log() / c ** 2) / 2
    return t1 + t2 - t3 + Zval / A


def j_slack(a, A, Zval):
    """J_A^+ - J_A <= 1/(18 A^3) + 2 Z(A+1/2)/A   (how much the closed form
    can overshoot; printed for information only)."""
    A = arb(A)
    return 1 / (18 * A ** 3) + 2 * Zval / A


def upper(x):
    """a float >= every point of the ball x (for printing only)"""
    return float((x.mid() + x.rad()).upper()) if hasattr(x, "upper") else float(x.mid()) + float(x.rad())


def fl(x):
    return float(x.mid())


def run(name, q, a, chi, zeros, Mlist, out):
    def P(*args):
        s = " ".join(str(t) for t in args)
        print(s)
        out.write(s + "\n")

    P("=" * 78)
    P("character %s   (q = %d, a = %d),  %d certified simple zeros" % (name, q, a, len(zeros)))
    T = zeros[-1] + arb("0.25")
    P("  T = last ordinate + 0.25 = %s" % T.str(12, radius=False))
    g1 = zeros[0]

    # ---- C(chi) and the plain tail of sum gamma^-2
    C = C_const(chi, a)
    S2 = sum((1 / g ** 2 for g in zeros), arb(0))
    tail2 = C / 2 - S2                       # = sum_{gamma>T} m_gamma / gamma^2
    P("  C(chi) = sum 2/gamma^2            = %s" % C.str(25, radius=True))
    P("  sum_{gamma>T} m_gamma/gamma^2     = %s" % tail2.str(12, radius=True))
    if not tail2 > 0:
        raise RuntimeError("tail2 not certified positive")

    # ---- log-weighted tail through J_A
    A = 30
    ZA = zeta_logderiv_abs(arb(A) + arb(1) / 2)
    JA = j_upper(q, a, A, ZA)
    SlogA = sum(((1 + g ** 2 / arb(A) ** 2).log() / g ** 2 for g in zeros), arb(0))
    taillog = JA - SlogA                     # >= sum_{gamma>T} m_gamma log(1+gamma^2/A^2)/gamma^2
    P("  J_A^+ (A = %d)                    = %s   (overshoot <= %.2e)"
      % (A, JA.str(15, radius=True), upper(j_slack(a, A, ZA))))
    P("  sum_{gamma<=T} log(1+g^2/A^2)/g^2 = %s" % SlogA.str(15, radius=True))
    P("  => log-weighted tail  <=           %s" % taillog.str(10, radius=True))
    # sanity: the log-tail must exceed log(1+T^2/A^2) * tail2
    lowlog = (1 + T ** 2 / arb(A) ** 2).log() * tail2
    P("  (sanity: lower bound log(1+T^2/A^2)*tail2 = %s)" % lowlog.str(10, radius=True))
    if not taillog > lowlog:
        raise RuntimeError("log-tail sanity check failed")

    # ---- multiplicity bound: choose delta minimising the excess constant
    best = None
    for dnum in range(55, 201, 5):
        delta = arb(dnum) / 100
        x0 = (arb(1) / 2 + delta + a) / 2
        if not (arb(A) > 2 * x0):
            continue
        Z = zeta_logderiv_abs(arb(1) / 2 + delta)
        b = (arb(q) / (2 * arb.pi())).log() / 2 + 1 / (24 * x0 ** 2) + Z
        # mu_delta(gamma) - 1 <= (delta/4) log(1+gamma^2/A^2) + cc   (gamma real, A >= 2 x0)
        cc = delta * (arb(A).log() / 2 + b) - 1
        # excess of sum m_gamma^2/gamma^2 over sum m_gamma/gamma^2 beyond T:
        exc = (delta / 4) * taillog + cc * tail2
        if best is None or upper(exc) < upper(best[1]):
            best = (dnum, exc, delta, x0, b, cc)
    dnum, exc, delta, x0, b, cc = best
    P("  multiplicity bound: delta = %.2f  (sigma = %.2f),  b_delta = %s"
      % (dnum / 100, 0.5 + dnum / 100, b.str(10, radius=False)))
    muT = (delta / 4) * (T ** 2 + 4 * x0 ** 2).log() + delta * b
    P("    mu_delta(T) = %s  (so every zero just above T has multiplicity <= %d)"
      % (muT.str(8, radius=False), int(upper(muT))))
    P("  E2 := bound for sum_{gamma>T} (m_gamma^2 - m_gamma)/gamma^2  [distinct gamma]")
    P("      = (delta/4)*logtail + cc*tail2 = %s" % exc.str(10, radius=True))
    P("  ratio E2 / tail2 = %.4f  (average excess multiplicity allowed in the tail)"
      % (fl(exc) / fl(tail2)))
    Cplus = C + 2 * exc
    P("  C^+(chi) := C(chi) + 2*E2 = %s   [V_M <= C^+ (2M)^2 for all M > 0]"
      % Cplus.str(10, radius=True))

    # ---- per-weight results
    P("")
    P("  %7s %7s | %-22s %-12s %-12s | %-11s %-11s %-11s" %
      ("m", "M", "sigma_M^2 (closed form)", "tau_M(T)", "V_M^+", "1/(1+sig^2)", "1/(1+V^+)", "soft-thr"))
    rows = []
    for mstr in Mlist:
        m = arb(mstr)
        M = m + arb(1) / 2
        s0 = M + arb(1) / 2
        sig2 = 4 * M * log_deriv_Lambda(chi, q, a, s0)
        head = sum((8 * M ** 2 / (M ** 2 + g ** 2) for g in zeros), arb(0))
        tau = sig2 - head                     # exact multiplicity-weighted tail
        if not tau > 0:
            raise RuntimeError("tau not certified positive at m=%s" % mstr)
        # eps_M(T) <= 8 M^2 * [ (delta/4) logtail ] + cc * tau   (cc may have either sign;
        # the tau-term is an identity, the log-term uses 1/(M^2+g^2) <= 1/g^2)
        eps = 8 * M ** 2 * (delta / 4) * taillog + cc * tau
        Vplus = sig2 + eps
        cant_simple = 1 / (1 + sig2)
        cant = 1 / (1 + Vplus)

        # ---- soft-threshold refinement (Theorem 5.3).  The peeled set S consists of
        # certified simple zeros only; any threshold t > 0 is admissible.
        amps = [4 * M / (M ** 2 + g ** 2).sqrt() for g in zeros]     # decreasing
        tau_min = arb(0)
        Vtail = tau + eps                      # >= sum_{gamma>T} m_gamma^2 |r|^2/2
        bestst = None
        ampf = [fl(x) for x in amps]
        # candidate thresholds: the amplitudes themselves and a grid in between
        cands = sorted(set([t for t in ampf if t > upper(tau_min)]
                           + [ampf[0] * k / 400 for k in range(1, 401)
                              if ampf[0] * k / 400 > upper(tau_min)]))
        for tf in cands:
            # float pre-screen
            B = sum(x - tf for x in ampf if x > tf)
            if B >= 1:
                continue
            Vt = sum(min(x, tf) ** 2 / 2 for x in ampf) + upper(Vtail)
            val = Vt / (Vt + (1 - B) ** 2)
            if bestst is None or val < bestst[0]:
                bestst = (val, tf)
        # certify the best threshold in ball arithmetic
        tf = bestst[1]
        tb = arb(repr(tf))
        Bb = arb(0)
        Vb = Vtail
        for x in amps:
            if x > tb:
                Bb += x - tb
                Vb += tb ** 2 / 2
            elif x < tb:
                Vb += x ** 2 / 2
            else:                               # undecided comparison: worst case for both
                Bb += abs(x - tb)
                Vb += (abs(x) + abs(x - tb)) ** 2 / 2
        if not (Bb < 1):
            raise RuntimeError("soft threshold B >= 1")
        st = 1 - Vb / (Vb + (1 - Bb) ** 2)

        P("  %7s %7.4f | %-22s %-12s %-12s | %-11.6f %-11.6f %-11.6f  (tau=%.4f, B=%.4f)" %
          (mstr, fl(M), sig2.str(18, radius=False), tau.str(8, radius=False),
           Vplus.str(9, radius=False), fl(cant_simple) - fl(cant_simple.rad()),
           fl(cant) - float(cant.rad()), fl(st) - float(st.rad()), tf, fl(Bb)))
        rows.append((mstr, sig2, tau, eps, Vplus, cant_simple, cant, st, tf, Bb))
    P("")
    P("  full-precision enclosures:")
    for (mstr, sig2, tau, eps, Vplus, cs, cant, st, tf, Bb) in rows:
        P("   m = %-6s sigma^2 = %s" % (mstr, sig2.str(25, radius=True)))
        P("            tau_M(T) = %s" % tau.str(15, radius=True))
        P("            eps_M(T) <= %s" % eps.str(12, radius=True))
        P("            V_M^+    = %s" % Vplus.str(12, radius=True))
        P("            1/(1+sigma^2) = %s   [if all zeros are simple]" % cs.str(12, radius=True))
        P("            1/(1+V^+)     = %s   [GRH only, certified]" % cant.str(12, radius=True))
        P("            soft-threshold bound = %s  at tau = %r" % (st.str(12, radius=True), tf))

    # ---- Tier B: assume in addition that ALL zeros with 0 < gamma <= X are simple
    #      (Platt's verification height: X = 10^8/q; see GRH-ONLY.md, Remark 4.6).
    X = arb(10 ** 8) / q
    ZX = zeta_logderiv_abs(X + arb(1) / 2)
    c_a = arb(1) / 2 + a
    # U_X >= sum_{gamma>X} m_gamma/gamma^2 :  1/g^2 <= 2/(X^2+g^2) for g >= X, and
    # sum_gamma 2/(X^2+g^2) = (Lambda'/Lambda)(X+1/2)/X <= [log(q/pi)/2 + log((X+c_a)/2)/2 + Z]/X
    UX = ((arb(q) / arb.pi()).log() / 2 + ((X + c_a) / 2).log() / 2 + ZX) / X
    JX = j_upper(q, a, X, ZX)                 # >= sum_gamma log(1+g^2/X^2)/g^2
    # for g > X:  log(g^2+4x0^2) <= log(X^2) + log(1+g^2/X^2)
    logtailX = JX + 2 * X.log() * UX
    excB = (delta / 4) * logtailX + (delta * b - 1) * UX
    P("")
    P("  Tier B (all zeros up to X = 10^8/q = %.4g assumed simple):" % fl(X))
    P("    sum_{gamma>X} (m_gamma^2-m_gamma)/gamma^2 <= %s" % excB.str(6, radius=False))
    CplusB = C + 2 * excB
    P("    C^+_B(chi) = %s" % CplusB.str(14, radius=False))
    for (mstr, sig2, tau, eps, Vplus, cs, cant, st, tf, Bb) in rows:
        M = arb(mstr) + arb(1) / 2
        VB = sig2 + 8 * M ** 2 * excB
        cB = 1 / (1 + VB)
        P("    m = %-6s V_M <= %s   1/(1+V) >= %.8f   (simple-zero value %.8f)"
          % (mstr, VB.str(12, radius=False), fl(cB) - float(cB.rad()), fl(cs)))
    return dict(C=C, Cplus=Cplus, tail2=tail2, exc=exc, rows=rows, delta=delta, x0=x0, b=b,
                muT=muT, T=T)


def main():
    t0 = time.time()
    outfn = os.path.join(HERE, "certified_bounds_results.txt")
    with open(outfn, "w") as out:
        chi4 = dirichlet_char(4, 3)
        chi3 = dirichlet_char(3, 2)
        z4 = load_zeros("chi4_zeros.txt")
        z3 = load_zeros("chi3_zeros.txt")
        Ml = ["-0.45", "-0.4", "-0.3", "-0.25", "-0.1", "0", "0.25", "0.5", "1"]
        r4 = run("chi_4", 4, 1, chi4, z4, Ml, out)
        r3 = run("chi_3", 3, 1, chi3, z3, Ml, out)
        s = "\nelapsed %.1f s" % (time.time() - t0)
        print(s)
        out.write(s + "\n")


if __name__ == "__main__":
    main()
