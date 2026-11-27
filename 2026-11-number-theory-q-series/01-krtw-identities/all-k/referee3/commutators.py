#!/usr/bin/env python3
"""Referee 3, check (3): operator identities on point masses, exact rational arithmetic at rational (q, x).

Functions live on N^2 (value 0 outside).  T is implemented from its definition
    (TF)(R,S) = sum_{r<=R, s>=S} w(r,s) F(r,s) / ((q;q)_{R-r} (q^2;q^2)_{s-S}),   w(r,s) = x^{r+2s} q^{(r+s)^2+s^2},
H', G, h, O from their definitions in F-all-k Sec. 2 and iv-all-k Sec. 2 (as pointwise formulas, written
out by me), T# = reversed chain, H'#, G# from the explicit formula (5.1).

Checked on a box, for point masses delta_(R0,S0) in a box:
  [T,H'] = 0, [T,G] = 0                      (F-all-k Prop. 2.1)
  T h - h T + x O T = 0                      (iv-all-k Lemma 2.1)
  T# = W^{-1} T^t W, (5.1) = W^{-1} H'^t W and W^{-1} G^t W entrywise      ([C1b,c] of the author)
  [T#,H'#] = 0, [T#,G#] = 0 on point masses
  omega(h delta) = omega(P^{-1} delta)       (iv-all-k Lemma 2.2, point masses; finite sums)
Negative controls: H' with the coefficient of its U-term changed; h with xPY replaced by xP.
"""
from fractions import Fraction as Fr
import itertools

POINTS = [(Fr(2, 7), Fr(3, 5)), (Fr(-3, 4), Fr(5, 11)), (Fr(5, 3), Fr(-2, 9))]
BOX_R, BOX_S = 5, 4      # evaluation box
PM_R, PM_S = 3, 3        # point masses


def make(q, x, Hcoef=1, hvariant=False):
    Q = q * q

    def poch(a, p, n):
        r = Fr(1)
        for i in range(n):
            r *= (1 - a * p ** i)
        return r

    def qq(n):           # (q;q)_n, 1/(q;q)_n := 0 for n < 0 handled by callers
        return poch(q, q, n)

    def QQ(n):
        return poch(Q, Q, n)

    def w(r, s):
        return x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s)

    def T(F, Rb, Sb):
        """F: dict with finite support; returns TF on [0..Rb]x[0..Sb]"""
        out = {}
        for R in range(Rb + 1):
            for S in range(Sb + 1):
                acc = Fr(0)
                for (r, s), v in F.items():
                    if r <= R and s >= S and v != 0:
                        acc += w(r, s) * v / (qq(R - r) * QQ(s - S))
                out[(R, S)] = acc
        return out

    def Tsharp(G, Rb, Sb):
        out = {}
        for r in range(Rb + 1):
            for s in range(Sb + 1):
                acc = Fr(0)
                for (R, S), v in G.items():
                    if R >= r and S <= s and v != 0:
                        acc += w(R, S) * v / (qq(R - r) * QQ(s - S))
                out[(r, s)] = acc
        return out

    def get(F, R, S):
        if R < 0 or S < 0:
            return Fr(0)
        return F.get((R, S), Fr(0))

    # local operators: output on [0..Rb]x[0..Sb]; input F must be known on what they read
    def Hp(F, Rb, Sb):
        out = {}
        for R in range(Rb + 1):
            for S in range(Sb + 1):
                v = (q ** (-2 * S) + x) * q ** (-R) * (get(F, R, S) - get(F, R - 1, S))
                v += x * q ** R * (1 + x * q ** (2 * S)) * get(F, R, S)
                v -= Hcoef * x * x * q ** (R + 2 * S) * get(F, R, S + 1)
                out[(R, S)] = v
        return out

    def Gop(F, Rb, Sb):
        out = {}
        for R in range(Rb + 1):
            for S in range(Sb + 1):
                f0, f1, f2 = get(F, R, S), get(F, R - 1, S), get(F, R - 2, S)
                v = q ** (-2 * R - 2 * S) * (f0 - (1 + q) * f1 + q * f2)
                v += q ** (-2 * S) * f0 - x * (1 + 1 / q) * f1
                v += x * x * (1 + q ** (2 * R)) * q ** (2 * S) * (f0 - get(F, R, S + 1))
                out[(R, S)] = v
        return out

    def h(F, Rb, Sb):
        out = {}
        for R in range(Rb + 1):
            for S in range(Sb + 1):
                v = q ** (-R) * (get(F, R, S) - get(F, R - 1, S))
                v += x * q ** R * (1 if hvariant else q ** (2 * S)) * get(F, R, S)
                out[(R, S)] = v
        return out

    def O(F, Rb, Sb):
        return {(R, S): q ** (R + 2 * S) * get(F, R, S + 1) for R in range(Rb + 1) for S in range(Sb + 1)}

    # (5.1), with Sn F(R,S) = F(R+1,S), un F(R,S) = F(R,S-1), hQ_a = Y^{-1}(1-un) + aY
    def Hps(F, Rb, Sb):
        out = {}
        for R in range(Rb + 1):
            for S in range(Sb + 1):
                v = q ** (-R) * q ** (-2 * S) * (get(F, R, S) - get(F, R, S - 1))
                v += x * q ** (-R) * get(F, R, S)
                v += x * q ** R * (1 + x * q ** (2 * S)) * (get(F, R, S) - get(F, R + 1, S))
                out[(R, S)] = v
        return out

    def Gs(F, Rb, Sb):
        out = {}
        for R in range(Rb + 1):
            for S in range(Sb + 1):
                a = x * x * q ** (2 * R)
                hq = q ** (-2 * S) * (get(F, R, S) - get(F, R, S - 1)) + a * q ** (2 * S) * get(F, R, S)
                v = (q ** (-2 * R) + 1) * hq
                v -= x * (1 + 1 / q) * get(F, R + 1, S)
                v -= x * x * (1 + q) * q ** (2 * R) * q ** (2 * S) * get(F, R + 1, S)
                v += q * x * x * q ** (2 * R) * q ** (2 * S) * get(F, R + 2, S)
                out[(R, S)] = v
        return out

    return dict(T=T, Tsharp=Tsharp, Hp=Hp, Gop=Gop, h=h, O=O, Hps=Hps, Gs=Gs, w=w, qq=qq, QQ=QQ)


def maxabs(D1, D2, keys):
    return max(abs(D1.get(k, 0) - D2.get(k, 0)) for k in keys)


def kernel(op, Rb, Sb, need_R, need_S):
    """matrix of a local operator on the box: K[(a,b)] = (op delta_b)(a)"""
    K = {}
    for b in itertools.product(range(need_R + 1), range(need_S + 1)):
        out = op({b: Fr(1)}, Rb, Sb)
        for a, v in out.items():
            if v != 0:
                K[(a, b)] = v
    return K


def run():
    allok = True
    for (q, x) in POINTS:
        for label, kw in (("exact", {}), ("NEG-CONTROL H' U-coefficient 2", {"Hcoef": 2}),
                          ("NEG-CONTROL h with xP instead of xPY", {"hvariant": True})):
            op = make(q, x, **kw)
            T, Hp, Gop, h, O = op['T'], op['Hp'], op['Gop'], op['h'], op['O']
            keys = [(R, S) for R in range(BOX_R + 1) for S in range(BOX_S + 1)]
            bad = {"[T,H']": 0, "[T,G]": 0, "Th-hT+xOT": 0}
            for R0 in range(PM_R + 1):
                for S0 in range(PM_S + 1):
                    d = {(R0, S0): Fr(1)}
                    # T H' d : H'd has finite support; compute H'd on a box large enough to contain it
                    Hd = Hp(d, R0 + 2, S0 + 1)
                    THd = T(Hd, BOX_R, BOX_S)
                    Td = T(d, BOX_R, BOX_S + 1)
                    HTd = Hp(Td, BOX_R, BOX_S)
                    bad["[T,H']"] = max(bad["[T,H']"], maxabs(THd, HTd, keys))
                    Gd = Gop(d, R0 + 3, S0 + 1)
                    TGd = T(Gd, BOX_R, BOX_S)
                    GTd = Gop(Td, BOX_R, BOX_S)
                    bad["[T,G]"] = max(bad["[T,G]"], maxabs(TGd, GTd, keys))
                    hd = h(d, R0 + 2, S0 + 1)
                    Thd = T(hd, BOX_R, BOX_S)
                    hTd = h(Td, BOX_R, BOX_S)
                    OTd = O(Td, BOX_R, BOX_S)
                    diff = {k: Thd[k] - hTd[k] + x * OTd[k] for k in keys}
                    bad["Th-hT+xOT"] = max(bad["Th-hT+xOT"], max(abs(v) for v in diff.values()))
            ok = all(v == 0 for v in bad.values())
            if label == "exact":
                allok &= ok
            print(f"(q,x)=({q},{x}) [{label}]: " + ", ".join(f"{k}: {'0' if v == 0 else 'NONZERO'}" for k, v in bad.items()))
        # transposes and reversed chain
        op = make(q, x)
        w = op['w']
        Rb, Sb = 5, 5
        # kernels on the box (inputs supported in the box, outputs on the box); for local operators
        # the entries with both indices inside a smaller box are complete
        KH = kernel(op['Hp'], Rb + 1, Sb + 1, Rb + 1, Sb + 1)
        KG = kernel(op['Gop'], Rb + 2, Sb + 1, Rb + 2, Sb + 1)
        KHs = kernel(op['Hps'], Rb, Sb, Rb + 2, Sb + 1)
        KGs = kernel(op['Gs'], Rb, Sb, Rb + 2, Sb + 1)
        inner = [(R, S) for R in range(Rb + 1) for S in range(Sb + 1)]
        errH = errG = 0
        for a in inner:
            for b in inner:
                tH = w(b) * KH.get((b, a), 0) / w(a) if False else w(*b) * KH.get((b, a), 0) / w(*a)
                tG = w(*b) * KG.get((b, a), 0) / w(*a)
                errH = max(errH, abs(tH - KHs.get((a, b), 0)))
                errG = max(errG, abs(tG - KGs.get((a, b), 0)))
        # T# = W^{-1} T^t W  (kernel entries on the box; T kernel entries are finite)
        errT = 0
        for a in inner:
            for b in inner:
                Tab = op['T']({b: Fr(1)}, Rb, Sb).get(a, 0)        # T(a,b)
                Tsba = op['Tsharp']({a: Fr(1)}, Rb, Sb).get(b, 0)  # T#(b,a)
                errT = max(errT, abs(Tsba - w(*a) * Tab / w(*b)))
        ok2 = errH == 0 and errG == 0 and errT == 0
        allok &= ok2
        # commutators of T# with H'#, G# on point masses; T# delta has finite R-support (r <= R0) and
        # infinite S-support; H'#, G# read (R+1,S),(R+2,S),(R,S-1), so compute T# delta on a box
        bad = {"[T#,H'#]": 0, "[T#,G#]": 0}
        keys = [(R, S) for R in range(4) for S in range(4)]
        for R0 in range(3):
            for S0 in range(3):
                d = {(R0, S0): Fr(1)}
                Tsd = op['Tsharp'](d, 6, 4)
                HsTsd = op['Hps'](Tsd, 3, 3)
                Hsd = op['Hps'](d, R0 + 1, S0 + 2)
                TsHsd = op['Tsharp'](Hsd, 3, 3)
                bad["[T#,H'#]"] = max(bad["[T#,H'#]"], maxabs(HsTsd, TsHsd, keys))
                Gsd = op['Gs'](d, R0 + 1, S0 + 2)
                TsGsd = op['Tsharp'](Gsd, 3, 3)
                GsTsd = op['Gs'](Tsd, 3, 3)
                bad["[T#,G#]"] = max(bad["[T#,G#]"], maxabs(GsTsd, TsGsd, keys))
        ok3 = all(v == 0 for v in bad.values())
        allok &= ok3
        # Lemma 2.2 of iv-all-k on point masses
        Q = q * q
        errw = 0
        for R0 in range(4):
            for S0 in range(4):
                d = {(R0, S0): Fr(1)}
                hd = op['h'](d, R0 + 2, S0 + 1)
                om_h = sum(w(r, s) * v / op['QQ'](s) for (r, s), v in hd.items())
                om_P = w(R0, S0) * q ** (-R0) / op['QQ'](S0)
                errw = max(errw, abs(om_h - om_P))
        allok &= errw == 0
        print(f"(q,x)=({q},{x}): (5.1) H'# = W^-1 H'^t W: {'OK' if errH == 0 else 'FAIL'}; G# = W^-1 G^t W: "
              f"{'OK' if errG == 0 else 'FAIL'}; T# = W^-1 T^t W: {'OK' if errT == 0 else 'FAIL'}; "
              + ", ".join(f"{k}: {'0' if v == 0 else 'NONZERO'}" for k, v in bad.items())
              + f"; Lemma 2.2 on point masses: {'OK' if errw == 0 else 'FAIL'}")
    print("ALL EXACT COMMUTATOR CHECKS " + ("PASSED" if allok else "FAILED"))


if __name__ == "__main__":
    run()
