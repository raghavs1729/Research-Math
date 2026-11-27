"""Exact position-space checks (rational q, x) of the operator facts used in the proof.

Functions on N^2 are dicts {(R,S): value}.  Operators (R,S >= 0, F(-1,.) := 0):
  (T F)(R,S)  = sum_{r<=R, s>=S} w(r,s) F(r,s) / ((q)_{R-r} (q^2;q^2)_{s-S}),   w(r,s) = x^{r+2s} q^{(r+s)^2+s^2}
  (T# G)(r,s) = sum_{R>=r, S<=s} w(R,S) G(R,S) / ((q)_{R-r} (q^2;q^2)_{s-S})
  P = q^R, Y = q^{2S}, (sh F)(R,S) = F(R-1,S), (up F)(R,S) = F(R,S+1), (shi F)(R,S) = F(R+1,S), (dn F)(R,S) = F(R,S-1)
  H'  = (Y^-1 + x) P^-1 (1 - sh) + x P (1 + x Y) - x^2 P Y up
  G   = P^-2 Y^-1 (1 - sh)(1 - q sh) + Y^-1 - x(1+q^-1) sh + x^2 Y (1 - up) + x^2 P^2 Y (1 - up)
  H'# = P^-1 Y^-1 (1 - dn) + x P^-1 + x P (1 + x Y)(1 - shi)
  G#  = (P^-2 + 1) Y^-1 (1 - dn) + x^2 Y (1 + P^2) - x(1+q^-1) shi - x^2 (1+q) P^2 Y shi + q x^2 P^2 Y shi^2
Checks: [T,H'] = [T,G] = 0 and [T#,H'#] = [T#,G#] = 0 on point masses, evaluated on a window.
Sums over s >= S in T are truncated at s <= SMAX: exact as long as the input has support in s <= SMAX.
"""
import sys
from fractions import Fraction as Fr


def qp(a, n, q):
    p = Fr(1)
    for i in range(n):
        p *= (1 - a * q ** i)
    return p


class Ops:
    def __init__(self, q, x, RMAX, SMAX):
        self.q, self.x, self.RMAX, self.SMAX = q, x, RMAX, SMAX
        self.iq = [1 / qp(q, n, q) for n in range(RMAX + SMAX + 10)]
        self.iQ = [1 / qp(q * q, n, q * q) for n in range(RMAX + SMAX + 10)]

    def w(self, r, s):
        q, x = self.q, self.x
        return x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s)

    def T(self, F):
        out = {}
        for (r, s), v in F.items():
            if v == 0:
                continue
            c = self.w(r, s) * v
            for R in range(r, self.RMAX + 1):
                for S in range(0, s + 1):
                    out[(R, S)] = out.get((R, S), 0) + c * self.iq[R - r] * self.iQ[s - S]
        return out

    def Ts(self, G):
        out = {}
        for (R, S), v in G.items():
            if v == 0:
                continue
            c = self.w(R, S) * v
            for r in range(0, R + 1):
                for s in range(S, self.SMAX + 1):
                    out[(r, s)] = out.get((r, s), 0) + c * self.iq[R - r] * self.iQ[s - S]
        return out

    # elementary operators
    def mul(self, F, f):
        return {k: f(*k) * v for k, v in F.items()}

    def sh(self, F):    # F(R-1,S)
        return {(R + 1, S): v for (R, S), v in F.items()}

    def up(self, F):    # F(R,S+1)
        return {(R, S - 1): v for (R, S), v in F.items() if S >= 1}

    def shi(self, F):   # F(R+1,S)
        return {(R - 1, S): v for (R, S), v in F.items() if R >= 1}

    def dn(self, F):    # F(R,S-1)
        return {(R, S + 1): v for (R, S), v in F.items()}

    @staticmethod
    def add(*Fs):
        out = {}
        for c, F in Fs:
            for k, v in F.items():
                out[k] = out.get(k, 0) + c * v
        return out

    def Hp(self, F):
        q, x = self.q, self.x
        a = self.add((1, F), (-1, self.sh(F)))
        t1 = self.mul(a, lambda R, S: (q ** (-2 * S) + x) * q ** (-R))
        t2 = self.mul(F, lambda R, S: x * q ** R * (1 + x * q ** (2 * S)))
        t3 = self.mul(self.up(F), lambda R, S: -x * x * q ** (R + 2 * S))
        return self.add((1, t1), (1, t2), (1, t3))

    def Gop(self, F):
        q, x = self.q, self.x
        a = self.add((1, F), (-1, self.sh(F)))
        b = self.add((1, a), (-q, self.sh(a)))           # (1 - q sh)(1 - sh) F  (they commute)
        t1 = self.mul(b, lambda R, S: q ** (-2 * R) * q ** (-2 * S))
        t2 = self.mul(F, lambda R, S: q ** (-2 * S))
        t3 = self.mul(self.sh(F), lambda R, S: -x * (1 + 1 / q))
        c = self.add((1, F), (-1, self.up(F)))
        t4 = self.mul(c, lambda R, S: x * x * q ** (2 * S) * (1 + q ** (2 * R)))
        return self.add((1, t1), (1, t2), (1, t3), (1, t4))

    def Hps(self, F):
        q, x = self.q, self.x
        a = self.add((1, F), (-1, self.dn(F)))
        t1 = self.mul(a, lambda R, S: q ** (-R) * q ** (-2 * S))
        t2 = self.mul(F, lambda R, S: x * q ** (-R))
        b = self.add((1, F), (-1, self.shi(F)))
        t3 = self.mul(b, lambda R, S: x * q ** R * (1 + x * q ** (2 * S)))
        return self.add((1, t1), (1, t2), (1, t3))

    def Gs(self, F):
        q, x = self.q, self.x
        a = self.add((1, F), (-1, self.dn(F)))
        t1 = self.mul(a, lambda R, S: (q ** (-2 * R) + 1) * q ** (-2 * S))
        t2 = self.mul(F, lambda R, S: x * x * q ** (2 * S) * (1 + q ** (2 * R)))
        t3 = self.mul(self.shi(F), lambda R, S: -x * (1 + 1 / q))
        t4 = self.mul(self.shi(F), lambda R, S: -x * x * (1 + q) * q ** (2 * R + 2 * S))
        t5 = self.mul(self.shi(self.shi(F)), lambda R, S: q * x * x * q ** (2 * R + 2 * S))
        return self.add((1, t1), (1, t2), (1, t3), (1, t4), (1, t5))


def window_diff(A, B, Rw, Sw):
    bad = 0
    for R in range(Rw + 1):
        for S in range(Sw + 1):
            if A.get((R, S), 0) != B.get((R, S), 0):
                bad += 1
    return bad


if __name__ == '__main__':
    tot = 0
    for (q, x) in [(Fr(3, 7), Fr(5, 11)), (Fr(-2, 9), Fr(13, 4))]:
        O = Ops(q, x, RMAX=14, SMAX=9)
        for r0 in range(0, 4):
            for s0 in range(0, 4):
                d = {(r0, s0): Fr(1)}
                # T H' vs H' T : window R <= 10, S <= 5 (sums exact: inputs supported in s <= s0+1 <= SMAX)
                b1 = window_diff(O.T(O.Hp(d)), O.Hp(O.T(d)), 10, 5)
                b2 = window_diff(O.T(O.Gop(d)), O.Gop(O.T(d)), 10, 5)
                # T# H'# vs H'# T# : T# output supported on r <= r0+2, s >= s0; window s <= SMAX-3, r <= r0+2
                b3 = window_diff(O.Ts(O.Hps(d)), O.Hps(O.Ts(d)), r0 + 2, 5)
                b4 = window_diff(O.Ts(O.Gs(d)), O.Gs(O.Ts(d)), r0 + 2, 5)
                tot += b1 + b2 + b3 + b4
                print('q=%s x=%s delta(%d,%d): [T,H\']:%d [T,G]:%d [T#,H\'#]:%d [T#,G#]:%d mismatches'
                      % (q, x, r0, s0, b1, b2, b3, b4))
                sys.stdout.flush()
    print('TOTAL MISMATCHES:', tot)
