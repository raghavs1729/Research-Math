#!/usr/bin/env python3
"""
certify.py q a u0 [tag] -- Arb certificate that

    S_a(u0) := sum_{chi mod q} conj(chi(a)) sum_{rho in Z(chi')} e^{i gamma u0} / (1/4 + gamma^2)
            <  -2 R_{q,a} - margin            (rho = 1/2 + i gamma, under GRH_q)

i.e. that the sum (7) of Akbary-Francis, which equals -S_a(log x) on GRH_q
(rho(rho-1) = -(1/4 + gamma^2)), exceeds 2 R_{q,a} at x = e^{u0}.

Inputs: zeros_q<q>.txt (sha256 pinned in SHA256SUMS; windows d -+ 1e-30, each
containing a zero 1/2 + i gamma, gamma > 0, of the labelled L-function; certified
by verify_zeros.py) and the closed form F_q (closed_form.py).

For a window (chi, gamma > 0): rho = 1/2 + i gamma is a zero of L(s, chi) and
conj(rho) is a zero of L(s, conj chi) (for real chi and for zeta: of the same
function).  Their joint contribution to S_a is
    2 w cos(gamma u0 - arg chi(a)),      w = 1/(1/4 + gamma^2),
and their joint contribution to F_q is 2 w.  All other zeros (GRH_q: on the line,
weights 1/(1/4+gamma^2) > 0) contribute to S_a at most their total weight
    tail = F_q - W,   W = sum over windows of 2 w
in absolute value.  Hence
    S_a(u0) <= S_T(u0) + (F_q - W).
(Only 'each window contains at least one zero, windows of one function are
disjoint' is used; completeness of the list is NOT needed for this inequality.)
"""
import json, os, sys
from fractions import Fraction
from flint import arb, acb, ctx
from afcommon import load_zeros, chi_arg_table, sha256_file, HERE, EPS
from closed_form import F_q


def R_qa(q, a, m=2):
    return sum(1 for b in range(1, q) if pow(b, m, q) == a % q)


def certified_sum(q, a, u0, zeros, Tmax=None):
    tab = chi_arg_table(q, a)
    e = arb(EPS)
    S = arb(0)
    W = arb(0)
    twopi = 2 * arb.pi()
    n = 0
    for lab, dstr in zeros:
        d = arb(dstr)
        if Tmax is not None and not d < Tmax:
            continue
        g = d + arb(0, e.upper())            # the ball [d - eps, d + eps]
        w = 2 / (arb(1) / 4 + g * g)
        ph = tab[lab]
        S += w * (g * u0 - twopi * ph.numerator / ph.denominator).cos()
        W += w
        n += 1
    return S, W, n


def main():
    q, a = int(sys.argv[1]), int(sys.argv[2])
    u0s = sys.argv[3]
    ctx.prec = 256
    assert pow(a, (q - 1) // 2, q) == 1, "a must be a square mod q (m = Ind_q(a) = 2)"
    R = R_qa(q, a)
    zeros = load_zeros(q)
    Fq, _ = F_q(q)
    u0 = arb(u0s)
    S, W, n = certified_sum(q, a, u0, zeros)
    tail = Fq - W
    assert tail > 0
    UB = S + tail
    thr = -2 * R
    ok = bool(UB < thr)
    margin = arb(thr) - UB
    print("q = %d, a = %d, R_{q,a} = %d, threshold 2R = %d" % (q, a, R, 2 * R))
    print("zero file sha256 = %s (%d windows)" % (sha256_file(os.path.join(HERE, "zeros_q%d.txt" % q)), n))
    print("u0   = %s      (x = e^u0 ~ 10^%.1f)" % (u0s, float(u0s) / 2.302585092994046))
    print("F_q  = %s   (closed form)" % Fq.str(25))
    print("W    = %s   (weight of the listed zeros and their conjugates)" % W.str(25))
    print("tail = F_q - W = %s" % tail.str(20))
    print("S_T(u0) = %s" % S.str(25))
    print("upper bound for S_a(u0):  S_T + tail = %s" % UB.str(25))
    print("CERTIFIED  S_a(u0) < -2R = %d :  %s ;  margin = -2R - (S_T + tail) = %s"
          % (thr, ok, margin.str(20)))
    print("alignment: -S_T/F_q = %.4f ; needed 2R/F_q = %.4f" % (-float(S.mid()) / float(Fq.mid()), 2 * R / float(Fq.mid())))
    cert = {"q": q, "a": a, "R": R, "u0": u0s, "F_q": Fq.str(30), "W": W.str(30), "tail": tail.str(20),
            "S_T": S.str(30), "upper_bound": UB.str(30), "certified": ok, "margin": margin.str(20),
            "zeros_sha256": sha256_file(os.path.join(HERE, "zeros_q%d.txt" % q)), "windows": n, "prec": ctx.prec}
    tag = ("_" + sys.argv[4]) if len(sys.argv) > 4 else ""
    json.dump(cert, open(os.path.join(HERE, "certificate_q%d_a%d%s.json" % (q, a, tag)), "w"), indent=1)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
