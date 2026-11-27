"""Referee's own helpers (no imports from the author's code).

Functions on N^2 are dicts {(R,S): value}; values are flint.fmpq (exact) at a rational point (q, x).
"""
from flint import fmpq


def poch(a, q, n):
    """(a;q)_n for n >= 0, exact; for n < 0 use (a;q)_n = 1/(a q^n; q)_{-n}."""
    if n >= 0:
        p = fmpq(1)
        for i in range(n):
            p *= (1 - a * q**i)
        return p
    return 1 / poch(a * q**n, q, -n)


def inv_qq(q, n):
    """1/(q;q)_n with 1/(q;q)_n := 0 for n < 0."""
    if n < 0:
        return fmpq(0)
    return 1 / poch(q, q, n)


class Point:
    def __init__(self, q, x):
        self.q = fmpq(q)
        self.x = fmpq(x)
        self.Q = self.q**2
        self._iq = {}
        self._iQ = {}

    def iq(self, n):  # 1/(q)_n
        if n not in self._iq:
            self._iq[n] = inv_qq(self.q, n)
        return self._iq[n]

    def iQ(self, n):  # 1/(Q;Q)_n
        if n not in self._iQ:
            self._iQ[n] = inv_qq(self.Q, n)
        return self._iQ[n]

    def w(self, r, s):
        q, x = self.q, self.x
        return x**(r + 2 * s) * q**((r + s)**2 + s**2)

    # ---------------- transfer operators ----------------
    def T(self, F, box):
        """(TF)(R,S) = sum_{r<=R, s>=S} w(r,s) F(r,s)/((q)_{R-r}(Q;Q)_{s-S}); F finite support dict."""
        out = {}
        for (R, S) in box:
            v = fmpq(0)
            for (r, s), f in F.items():
                if r <= R and s >= S and r >= 0 and s >= 0:
                    v += self.w(r, s) * f * self.iq(R - r) * self.iQ(s - S)
            out[(R, S)] = v
        return out

    def Tsharp(self, G, box):
        """(T#G)(r,s) = sum_{R>=r, S<=s} w(R,S) G(R,S)/((q)_{R-r}(Q;Q)_{s-S}); G finite support dict."""
        out = {}
        for (r, s) in box:
            v = fmpq(0)
            for (R, S), g in G.items():
                if R >= r and S <= s and R >= 0 and S >= 0:
                    v += self.w(R, S) * g * self.iq(R - r) * self.iQ(s - S)
            out[(r, s)] = v
        return out

    # ---------------- local operators (as printed) ----------------
    @staticmethod
    def get(F, R, S):
        if R < 0 or S < 0:
            return fmpq(0)
        return F.get((R, S), fmpq(0))

    def Hp(self, F, box):
        """H' = (Y^-1 + x) P^-1 (1 - Sft) + x P (1 + x Y) - x^2 P Y u ;  Sft F(R,S)=F(R-1,S), u F(R,S)=F(R,S+1)."""
        q, x = self.q, self.x
        g = self.get
        out = {}
        for (R, S) in box:
            P, Y = q**R, q**(2 * S)
            v = (1 / Y + x) / P * (g(F, R, S) - g(F, R - 1, S)) + x * P * (1 + x * Y) * g(F, R, S) \
                - x**2 * P * Y * g(F, R, S + 1)
            out[(R, S)] = v
        return out

    def Gop(self, F, box):
        """G = P^-2 Y^-1 (1-Sft)(1-q Sft) + Y^-1 - x(1+q^-1) Sft + x^2 (1+P^2) Y (1-u)."""
        q, x = self.q, self.x
        g = self.get
        out = {}
        for (R, S) in box:
            P, Y = q**R, q**(2 * S)
            f0, f1, f2 = g(F, R, S), g(F, R - 1, S), g(F, R - 2, S)
            # (1-Sft)(1-q Sft) = 1 - (1+q) Sft + q Sft^2
            v = (f0 - (1 + q) * f1 + q * f2) / (P**2 * Y) + f0 / Y - x * (1 + 1 / q) * f1 \
                + x**2 * (1 + P**2) * Y * (f0 - g(F, R, S + 1))
            out[(R, S)] = v
        return out

    def Hps(self, F, box):
        """H'# = P^-1 Y^-1 (1 - u#) + x P^-1 + x P (1 + x Y)(1 - S#);  S# F(R,S)=F(R+1,S), u# F(R,S)=F(R,S-1)."""
        q, x = self.q, self.x
        g = self.get
        out = {}
        for (R, S) in box:
            P, Y = q**R, q**(2 * S)
            v = (g(F, R, S) - g(F, R, S - 1)) / (P * Y) + x / P * g(F, R, S) \
                + x * P * (1 + x * Y) * (g(F, R, S) - g(F, R + 1, S))
            out[(R, S)] = v
        return out

    def Gs(self, F, box):
        """G# = (P^-2+1) h^Q_{x^2P^2} - x(1+q^-1) S# - x^2(1+q) P^2 Y S# + q x^2 P^2 Y S#^2,
        h^Q_a = Y^-1 (1 - u#) + a Y."""
        q, x = self.q, self.x
        g = self.get
        out = {}
        for (R, S) in box:
            P, Y = q**R, q**(2 * S)
            f0 = g(F, R, S)
            hQ = (f0 - g(F, R, S - 1)) / Y + x**2 * P**2 * Y * f0
            v = (1 / P**2 + 1) * hQ - x * (1 + 1 / q) * g(F, R + 1, S) \
                - x**2 * (1 + q) * P**2 * Y * g(F, R + 1, S) + q * x**2 * P**2 * Y * g(F, R + 2, S)
            out[(R, S)] = v
        return out
