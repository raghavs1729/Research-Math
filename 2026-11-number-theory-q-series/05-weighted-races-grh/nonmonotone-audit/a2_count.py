#!/usr/bin/env python3
"""
a2_count.py -- independent re-implementation of the zero count (audit item 1).

Contour (my variant): R = (-2, 3) x (-T, T), boundary split into four congruent
quarters by the symmetries Lambda(s) = Lambda(1-s) (root number +1) and
Lambda(conj s) = conj Lambda(s);  C : 3 -> 3+iT -> 1/2+iT.
Lambda(s) = (q/pi)^((s+1)/2) Gamma((s+1)/2) L(s,chi) (chi odd) has no zero on
Re s = 3 (Euler product) and none on Re s = -2 (mirror), and is nonzero at the
real points -2, 3.  So

   #zeros(R) = (2/pi) [ theta(T) + Arg L(3+iT) + Delta_{3 -> 1/2} arg L(sigma+iT) ],
   theta(T) = (T/2) log(q/pi) + Im loggamma(3/4 + iT/2),

since (s+1)/2 runs in Re > 0 (principal loggamma continuous there; vanishes at
s=3 real), |L(3+it) - 1| <= zeta(3)-1 < 1 keeps Re L(3+it) > 0 for all t, and
L(3) > 0.

Horizontal part: sub-segment [a,b]+iT, midpoint c.  I evaluate Arb's L on the
complex ball centred c+iT with radius (b-a)/2 (a DISC-type ball here: acb with
real radius (b-a)/2), F = L(ball), and the point value Lc = L(c+iT).  If
sup|F - Lc| < inf|Lc| then L(segment) lies in the open disc of centre Lc and
radius |Lc| (strictly), which lies in an open half-plane through 0, so the
continuous change of arg over [a,b] equals Arg(L(a)/Lc) + Arg(Lc/L(b))... with
both terms in (-pi/2, pi/2).  Otherwise bisect.  (Different criterion and
different bookkeeping from zeros_cert.py, which tests 0 notin rectangle.)

The result W must be an even integer 2N; N = #zeros with 0 < gamma < T (no
real zeros of Lambda then follows from comparing with the sign-change count).
"""
import sys, time
from flint import arb, acb, ctx, dirichlet_char

ctx.prec = 128
SIG_RIGHT = 3


def count(chi, q, T):
    T = arb(T)
    L3 = chi.l(acb(arb(SIG_RIGHT), T))
    if not L3.real > 0:
        raise ArithmeticError("Re L(3+iT) not > 0")
    vert = L3.arg()
    theta = T / 2 * (arb(q) / arb.pi()).log() + acb(arb(3) / 4, T / 2).lgamma().imag
    total = arb(0)
    nev = 0
    stack = [(arb(SIG_RIGHT), arb(1) / 2)]       # process from right to left
    while stack:
        a, b = stack.pop()                       # a > b, going from a to b
        c = (a + b) / 2
        rad = (a - b) / 2
        ball = acb(arb(c.mid(), rad.upper()), T)
        F = chi.l(ball)
        Lc = chi.l(acb(c, T))
        nev += 2
        dev = abs(F - Lc).upper()
        if F.is_finite() and dev < abs(Lc).lower():
            La = chi.l(acb(a, T))
            Lb = chi.l(acb(b, T))
            nev += 2
            d = (Lc / La).arg() + (Lb / Lc).arg()
            if not (d.rad() < 1e-10):
                raise ArithmeticError("arg not resolved")
            total += d
        else:
            if (a - b) < arb(2) ** -30:
                raise ArithmeticError("segment does not certify near sigma=%s" % c)
            stack.append((c, b))    # pushed first -> processed after (a,c)
            stack.append((a, c))
    W = 2 * (theta + vert + total) / arb.pi()
    return W, nev


if __name__ == "__main__":
    q = int(sys.argv[1]); l = int(sys.argv[2])
    chi = dirichlet_char(q, l)
    for Ts in sys.argv[3:]:
        t0 = time.time()
        num, den = (Ts.split("/") + ["1"])[:2]
        T = arb(int(num)) / int(den)
        W, nev = count(chi, q, T)
        n = W.unique_fmpz()
        print("q=%d T=%s (=%.6f): W = %s -> %s ; %d L-evals ; %.0fs" %
              (q, Ts, float(T.mid()), W.str(12, radius=True), n, nev, time.time() - t0), flush=True)
