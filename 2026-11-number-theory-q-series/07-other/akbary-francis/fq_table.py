#!/usr/bin/env python3
"""
fq_table.py QMAX -- F_q (Akbary-Francis (9)) for ALL moduli q <= QMAX, composite included,
via the closed forms (49), (50):
   F_q = sum_{d | q} P(d),   P(1) = F_1 = 2 + gamma - log(4 pi),
   P(d) = sum over PRIMITIVE chi mod d of [ log(d/pi) + 2 Re L'/L(1, chi) + psi((1+k)/2) ]
(every character mod q is induced by exactly one primitive character of conductor d | q),
together with R_{q,1} = #{b mod q : b^2 = 1} and the comparison F_q < 2 R_{q,1}
(the hypothesis of Akbary-Francis Prop. 1.20 / 1.21, i.e. the cases in which their
method proves 'at most finitely many exceptions to (3)' for every square class).
Arb ball arithmetic, 100 bits.  Output: fq_table_<QMAX>.txt
"""
import sys, os, math
from flint import arb, acb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
ctx.prec = 100
QMAX = int(sys.argv[1])


def P(d):
    if d == 1:
        return 2 + arb.const_euler() - (4 * arb.pi()).log()
    g0 = [None] + [acb.stieltjes(0, acb(arb(a) / d)) for a in range(1, d)]
    g1 = [None] + [acb.stieltjes(1, acb(arb(a) / d)) for a in range(1, d)]
    units = [a for a in range(1, d) if math.gcd(a, d) == 1]
    tot = arb(0)
    logd = arb(d).log()
    base = (arb(d) / arb.pi()).log()
    for n in units:
        chi = dirichlet_char(d, n)
        if not chi.is_primitive():
            continue
        s0 = acb(0)
        s1 = acb(0)
        for a in units:
            c = chi(a)
            s0 += c * g0[a]
            s1 += c * g1[a]
        L1 = s0 / d
        dL1 = -logd * L1 - s1 / d
        k = int(chi.parity())
        psi = -arb.const_euler() - (0 if k == 1 else 2 * arb(2).log())
        tot += base + 2 * (dL1 / L1).real + psi
    return tot


def main():
    Pd = {}
    rows = []
    for q in range(1, QMAX + 1):
        F = arb(0)
        for d in range(1, q + 1):
            if q % d == 0:
                if d % 4 == 2:
                    continue          # no primitive characters of conductor 2 mod 4
                if d not in Pd:
                    Pd[d] = P(d)
                F += Pd[d]
        R = sum(1 for b in range(1, q + 1) if math.gcd(b, q) == 1 and (b * b - 1) % q == 0)
        if F < 2 * R:
            verdict = "F<2R"
        elif F > 2 * R:
            verdict = "F>2R"
        else:
            verdict = "undecided"
        rows.append((q, F, R, verdict))
        print("%4d  F_q = %-22s  2R = %3d  %s" % (q, F.str(12, radius=False), 2 * R, verdict), flush=True)
    with open(os.path.join(HERE, "fq_table_%d.txt" % QMAX), "w") as fh:
        fh.write("# q  F_q  2R_{q,1}  verdict   (closed forms (49),(50) of Akbary-Francis; Arb)\n")
        for q, F, R, v in rows:
            fh.write("%d %s %d %s\n" % (q, F.str(12, radius=False), 2 * R, v))
    good = [q for q, F, R, v in rows if v == "F<2R"]
    print("q <= %d with F_q < 2R_{q,1}: %s" % (QMAX, good))


if __name__ == "__main__":
    main()
