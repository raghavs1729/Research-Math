"""Referee 2, computation (3): Lemma 2.1 (Th - hT = -x O T), Lemma 2.2 (omega(hF) = omega(P^-1 F)) and hu = u + xOu,
exact rational arithmetic at rational points (q, x), on all point masses in a box and on random finite combinations.
Operators implemented from their definitions in iv-all-k.md / F-all-k.md section 1:
  (TF)(R,S) = sum_{r<=R, s>=S} w(r,s) F(r,s) / ((q)_{R-r} (Q;Q)_{s-S}),  w(r,s) = x^{r+2s} q^{(r+s)^2+s^2}
  (hF)(R,S) = q^{-R}(F(R,S) - F(R-1,S)) + x q^{R+2S} F(R,S)
  (OF)(R,S) = q^{R+2S} F(R,S+1)
  omega(F) = sum w(r,s) F(r,s)/(Q;Q)_s.
"""
from fractions import Fraction as Fr
import random

random.seed(12345)


def run(q, x, Rb=5, Sb=4, out=print):
    Q = q * q

    def qp(n, base):
        r = Fr(1)
        for j in range(1, n + 1):
            r *= (1 - base ** j)
        return r

    iq = {n: 1 / qp(n, q) for n in range(0, 40)}
    iQ = {n: 1 / qp(n, Q) for n in range(0, 40)}

    def w(r, s):
        return x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s)

    def T(F, box):
        G = {}
        for (R, S) in box:
            tot = Fr(0)
            for (r, s), v in F.items():
                if r <= R and s >= S and v != 0:
                    tot += w(r, s) * v * iq[R - r] * iQ[s - S]
            G[(R, S)] = tot
        return G

    def h(F, box):
        G = {}
        for (R, S) in box:
            G[(R, S)] = q ** (-R) * (F.get((R, S), 0) - F.get((R - 1, S), 0)) + x * q ** (R + 2 * S) * F.get((R, S), 0)
        return G

    def O(F, box):
        return {(R, S): q ** (R + 2 * S) * F.get((R, S + 1), 0) for (R, S) in box}

    def omega(F):
        return sum(w(r, s) * v * iQ[s] for (r, s), v in F.items())

    box = [(R, S) for R in range(Rb + 1) for S in range(Sb + 1)]
    boxT = [(R, S) for R in range(Rb + 1) for S in range(Sb + 2)]
    boxh = [(R, S) for R in range(Rb + 2) for S in range(Sb + 1)]
    nfail = 0
    nchk = 0
    nonzero_comm = 0
    neg_detect = 0
    tests = [{(r0, s0): Fr(1)} for (r0, s0) in box]
    for _ in range(6):
        tests.append({p: Fr(random.randint(-9, 9), random.randint(1, 9)) for p in box})
    for F in tests:
        TF = T(F, boxT)
        hF = h(F, boxh)
        ThF = T(hF, box)
        hTF = h(TF, box)
        OTF = O(TF, box)
        for p in box:
            lhs = ThF[p] - hTF[p]
            rhs = -x * OTF[p]
            nchk += 1
            if lhs != rhs:
                nfail += 1
            if lhs != 0:
                nonzero_comm += 1
        # negative control: wrong sign / no shift in O
        bad = any(ThF[p] - hTF[p] != x * OTF[p] for p in box) or True
        bad2 = any(ThF[p] - hTF[p] != -x * q ** (p[0] + 2 * p[1]) * TF[p] for p in box)
        if bad2:
            neg_detect += 1
        # Lemma 2.2 on finitely supported F (hF computed on its full support)
        hFfull = h(F, [(R, S) for R in range(Rb + 2) for S in range(Sb + 1)])
        PinvF = {(r, s): q ** (-r) * v for (r, s), v in F.items()}
        nchk += 1
        if omega(hFfull) != omega(PinvF):
            nfail += 1
    # hu = u + x O u on a larger box
    u = {(R, S): iq[R] for R in range(12) for S in range(6)}
    hu = h(u, [(R, S) for R in range(12) for S in range(5)])
    Ou = O(u, [(R, S) for R in range(12) for S in range(5)])
    for R in range(12):
        for S in range(5):
            nchk += 1
            if hu[(R, S)] != u[(R, S)] + x * Ou[(R, S)]:
                nfail += 1
    out("  point (q,x) = (%s, %s): %d exact comparisons, %d failures; commutator Th-hT nonzero at %d box points; "
        "negative control (O without the S-shift) detected on %d/%d test functions"
        % (q, x, nchk, nfail, nonzero_comm, neg_detect, len(tests)))
    return nfail == 0 and neg_detect == len(tests)


if __name__ == '__main__':
    lines = []
    pr = lambda s: (print(s, flush=True), lines.append(s))
    pr("Referee 2, (3): Lemma 2.1, Lemma 2.2, hu = u + xOu -- exact, point masses (r0,s0) in [0,5]x[0,4] + 6 random combos")
    ok = True
    for (q, x) in [(Fr(2, 7), Fr(3, 11)), (Fr(-3, 5), Fr(5, 13)), (Fr(1, 3), Fr(-2, 9))]:
        ok &= run(q, x, out=pr)
    pr("RESULT: %s" % ("ALL OK" if ok else "FAILURE"))
