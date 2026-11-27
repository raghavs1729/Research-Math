#!/usr/bin/env python3
"""
closed_form.py -- F_q in closed form (Akbary-Francis (49), (50); Montgomery-Vaughan
Cor. 10.14 and 10.18), in Arb ball arithmetic, for a PRIME q.

  F(trivial) = sum_rho 1/(rho(1-rho)) over nontrivial zeros of zeta
             = 2 + gamma - log(4 pi)
  F(chi)     = sum_rho 1/(rho(1-rho)) over nontrivial zeros of L(s, chi), chi primitive mod q
             = log(q/pi) + 2 Re (L'/L)(1, chi) + psi((1+k)/2),   k = parity (0 even, 1 odd)
               [psi(1/2) = -gamma - 2 log 2, psi(1) = -gamma]
  F_q        = F(trivial) + sum_{chi != chi_0} F(chi)      (q prime: all chi != chi_0 primitive)

L(1, chi) and L'(1, chi) from the Laurent expansion of the Hurwitz zeta function,
  zeta(s, a) = 1/(s-1) + sum_n (-1)^n gamma_n(a) (s-1)^n / n!:
  L(s, chi) = q^(-s) sum_{a=1}^{q-1} chi(a) zeta(s, a/q), and sum_a chi(a) = 0, so
  L(1, chi)  = (1/q) sum_a chi(a) gamma_0(a/q)
  L'(1, chi) = -log(q) L(1, chi) - (1/q) sum_a chi(a) gamma_1(a/q).
Cross-check: L(1, chi) against Arb's acb_dirichlet_l, and L'(1,chi) against a
(non-rigorous) central difference.

Output: per-character F(chi) and F_q as balls; JSON to closed_form_q<q>.json
"""
import json, os, sys
from flint import arb, acb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))


def F_trivial():
    return 2 + arb.const_euler() - (4 * arb.pi()).log()


def L1_and_derivative(q, n):
    chi = dirichlet_char(q, n)
    s0 = acb(0)
    s1 = acb(0)
    for a in range(1, q):
        c = chi(a)
        s0 += c * acb.stieltjes(0, acb(arb(a) / q))
        s1 += c * acb.stieltjes(1, acb(arb(a) / q))
    L1 = s0 / q
    dL1 = -arb(q).log() * L1 - s1 / q
    return chi, L1, dL1


def F_chi(q, n):
    chi, L1, dL1 = L1_and_derivative(q, n)
    k = int(chi.parity())
    psi = -arb.const_euler() - (0 if k == 1 else 2 * arb(2).log())
    val = (arb(q) / arb.pi()).log() + 2 * (dL1 / L1).real + psi
    return val, L1, dL1, k


def F_q(q, verbose=False):
    tot = F_trivial()
    per = {"zeta": tot}
    for n in range(2, q):
        v, L1, dL1, k = F_chi(q, n)
        per[str(n)] = v
        tot += v
        if verbose:
            chi = dirichlet_char(q, n)
            Lref = chi.l(acb(1))
            assert (L1 - Lref).contains(acb(0)) or abs(complex(L1 - Lref)) < 1e-25
            old = ctx.prec
            ctx.prec = 400
            h = arb("1e-30")
            num = (chi.l(acb(1 + h)) - chi.l(acb(1 - h))) / (2 * h)
            ctx.prec = old
            assert abs(complex(num) - complex(dL1)) < 1e-20, (num, dL1)
            print("chi_%d(%2d,.) k=%d  L(1)=%s  L'/L(1)=%s  F(chi)=%s"
                  % (q, n, k, L1.str(15, radius=False), (dL1 / L1).str(15, radius=False), v.str(30)))
    return tot, per


if __name__ == "__main__":
    ctx.prec = 200
    for q in [int(a) for a in sys.argv[1:]] or [13, 11]:
        tot, per = F_q(q, verbose=True)
        print("F(zeta) = F_1 = %s" % per["zeta"].str(30))
        print("F_%d = %s" % (q, tot.str(40)))
        json.dump({"q": q, "F_q": tot.str(40), "per_character": {k: v.str(40) for k, v in per.items()}},
                  open(os.path.join(HERE, "closed_form_q%d.json" % q), "w"), indent=1)
