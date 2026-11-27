#!/usr/bin/env python3
"""Referee 3: x-adic valuations (F-all-k Lemma 6.2, ord u~ in Thm 6.4), with x an indeterminate and q a fixed
rational number; exact arithmetic in Q(x) (numerator/denominator fmpq_poly, reduced by gcd).

phi_{N,M} and phi#_{N,M} are computed by the same position-space triangular solves as in eigen.py, now over
Q(x); their Bailey-basis coefficients c^{(S)}_e (base q, a = x q^{2S}) and c_{e,rho} (base Q, a = x^2 q^{2 rho})
are obtained by triangular inversion, and their x-adic orders are compared with Lemma 6.2:
   (a) ord c^{(S)}_e >= 2(M-S);   (b) ord c_{e,rho} >= 0;
   (c) ord psi_{N',M'}(Xi_{e,sigma}) >= e + 2 sigma, = for (e,sigma) = (L',M') (nonzero), and = 0 unless
       e <= L', sigma >= M'.
Also prints min ord over each slice, to see whether the bounds are attained.
"""
import sys
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3
QS = [Fr(2, 7), Fr(-3, 5)]


def lowdeg(p):
    cs = p.coeffs()
    for i, c in enumerate(cs):
        if c != 0:
            return i
    return None


class RF:
    __slots__ = ("n", "d")

    def __init__(self, n, d=None):
        if d is None:
            d = fmpq_poly([1])
        if n.is_zero():
            self.n, self.d = n, fmpq_poly([1])
            return
        g = n.gcd(d)
        n, d = n // g, d // g
        lc = d.leading_coefficient()
        self.n, self.d = n / lc, d / lc

    @staticmethod
    def c(v):
        return RF(fmpq_poly([fmpq(v.numerator, v.denominator)]))

    def __add__(self, o):
        return RF(self.n * o.d + o.n * self.d, self.d * o.d)

    def __sub__(self, o):
        return RF(self.n * o.d - o.n * self.d, self.d * o.d)

    def __mul__(self, o):
        return RF(self.n * o.n, self.d * o.d)

    def __truediv__(self, o):
        assert not o.n.is_zero()
        return RF(self.n * o.d, self.d * o.n)

    def iszero(self):
        return self.n.is_zero()

    def ord(self):
        if self.n.is_zero():
            return 10 ** 9
        return lowdeg(self.n) - lowdeg(self.d)


def xpow(k, coef=Fr(1)):
    return RF(fmpq_poly([0] * k + [fmpq(coef.numerator, coef.denominator)]))


ONE = RF.c(Fr(1))
ZERO = RF(fmpq_poly([]))


def lin(c0, c1, k):
    """c0 + c1 x^k"""
    cs = [fmpq(0)] * (k + 1)
    cs[0] += fmpq(c0.numerator, c0.denominator)
    cs[k] += fmpq(c1.numerator, c1.denominator)
    return RF(fmpq_poly(cs))


class Ctx:
    def __init__(self, q):
        self.q = q
        self.Q = q * q
        self._cache = {}

    def qq(self, n, p=None):
        p = self.q if p is None else p
        r = Fr(1)
        for i in range(1, n + 1):
            r *= 1 - p ** i
        return r

    def xpoch(self, xp, c, step, n):
        """prod_{i<n} (1 - x^xp c step^i)"""
        r = ONE
        for i in range(n):
            r = r * lin(Fr(1), -c * step ** i, xp)
        return r

    def w(self, r, s):
        return xpow(r + 2 * s, self.q ** ((r + s) ** 2 + s * s))

    def lam(self, N, M):
        return xpow(N + M, self.q ** (N * N + M * M))

    def phi_e(self, e, R, S):
        """phi_e(R; x q^{2S}) = 1/((q)_{R-e} (x q^{2S+1};q)_{R+e})"""
        if R < e:
            return ZERO
        key = ("p", e, R, S)
        if key not in self._cache:
            self._cache[key] = ONE / (RF.c(self.qq(R - e)) * self.xpoch(1, self.q ** (2 * S + 1), self.q, R + e))
        return self._cache[key]

    def phiQ_e(self, e, S, rho):
        """phi^Q_e(S; x^2 q^{2 rho}) = 1/((Q;Q)_{S-e} (x^2 q^{2rho} Q;Q)_{S+e})"""
        if S < e:
            return ZERO
        key = ("Q", e, S, rho)
        if key not in self._cache:
            self._cache[key] = ONE / (RF.c(self.qq(S - e, self.Q)) * self.xpoch(2, self.q ** (2 * rho + 2), self.Q, S + e))
        return self._cache[key]

    def eig(self, N, M, Rb):
        L = N - M
        lam = self.lam(N, M)
        ph = {}
        for S in range(M, -1, -1):
            for R in range(Rb + 1):
                rhs = ZERO
                for s in range(S, M + 1):
                    for r in range(R + 1):
                        if (r, s) == (R, S):
                            continue
                        v = ph.get((r, s))
                        if v is not None and not v.iszero():
                            rhs = rhs + self.w(r, s) * v / RF.c(self.qq(R - r) * self.qq(s - S, self.Q))
                if (R, S) == (L, M):
                    assert rhs.iszero()
                    ph[(R, S)] = ONE / self.xpoch(1, self.q ** (2 * M + 1), self.q, 2 * L)
                else:
                    ph[(R, S)] = rhs / (lam - self.w(R, S))
        return ph

    def eig_sharp(self, N, M, Sb):
        L = N - M
        lam = self.lam(N, M)
        G = {}
        for r in range(L, -1, -1):
            for s in range(Sb + 1):
                rhs = ZERO
                for R in range(r, L + 1):
                    for S in range(s + 1):
                        if (R, S) == (r, s):
                            continue
                        v = G.get((R, S))
                        if v is not None and not v.iszero():
                            rhs = rhs + self.w(R, S) * v / RF.c(self.qq(R - r) * self.qq(s - S, self.Q))
                if (r, s) == (L, M):
                    assert rhs.iszero()
                    G[(r, s)] = ONE / self.xpoch(2, self.q ** (2 * L + 2), self.Q, 2 * M)
                else:
                    G[(r, s)] = rhs / (lam - self.w(r, s))
        return G


def run(q):
    C = Ctx(q)
    ok = True
    lines = []
    modes = [(N, M) for N in range(NMAX + 1) for M in range(N + 1)]
    sharp = {}
    for (N, M) in modes:
        L = N - M
        Rb = L + 2 * M + 2
        ph = C.eig(N, M, Rb)
        slice_min = []
        for S in range(M + 1):
            c = []
            for R in range(Rb + 1):
                v = ph[(R, S)]
                for e in range(R):
                    v = v - c[e] * C.phi_e(e, R, S)
                c.append(v / C.phi_e(R, R, S))
            ords = [cc.ord() for cc in c if not cc.iszero()]
            mn = min(ords)
            slice_min.append(mn)
            if mn < 2 * (M - S):
                ok = False
                lines.append(f"  FAIL (a): mode {(N, M)} slice {S}: min ord {mn} < {2 * (M - S)}")
        Sb = M + L + 2
        G = C.eig_sharp(N, M, Sb)
        sharp[(N, M)] = G
        rho_min = []
        for rho in range(L + 1):
            c = []
            for S in range(Sb + 1):
                v = G[(rho, S)]
                for e in range(S):
                    v = v - c[e] * C.phiQ_e(e, S, rho)
                c.append(v / C.phiQ_e(S, S, rho))
            ords = [cc.ord() for cc in c if not cc.iszero()]
            mn = min(ords)
            rho_min.append(mn)
            if mn < 0:
                ok = False
                lines.append(f"  FAIL (b): mode {(N, M)} slice rho={rho}: min ord {mn}")
        lines.append(f"mode {(N, M)}: (a) min ord of c^(S), S=0..M: {slice_min} (bound 2(M-S) = "
                     f"{[2 * (M - S) for S in range(M + 1)]});  (b) min ord of c_(.,rho), rho=0..L: {rho_min}")
    # (c) psi_{N',M'}(Xi_{e,sigma})
    cok = True
    for (Np, Mp) in modes:
        Lp = Np - Mp
        G = sharp[(Np, Mp)]
        for e in range(0, Lp + 3):
            for sig in range(0, Mp + 3):
                val = ZERO
                for r in range(e, Lp + 1):
                    gv = G.get((r, sig))
                    if gv is None:
                        continue
                    val = val + C.w(r, sig) * gv * C.phi_e(e, r, sig)
                if e > Lp or sig < Mp:
                    if not val.iszero():
                        cok = False
                        lines.append(f"  FAIL (c) support: psi_{(Np, Mp)}(Xi_{(e, sig)}) != 0")
                    continue
                if val.ord() < e + 2 * sig:
                    cok = False
                    lines.append(f"  FAIL (c) ord: psi_{(Np, Mp)}(Xi_{(e, sig)}) ord {val.ord()} < {e + 2 * sig}")
                if (e, sig) == (Lp, Mp) and (val.iszero() or val.ord() != e + 2 * sig):
                    cok = False
                    lines.append(f"  FAIL (c) equality at {(Lp, Mp)}")
    ok &= cok
    lines.append(f"(c) for all modes N' <= {NMAX}, e <= L'+2, sigma <= M'+2: {'OK' if cok else 'FAIL'}")
    return ok, lines


if __name__ == "__main__":
    allok = True
    for q in QS:
        ok, lines = run(q)
        print(f"=== q = {q} (x indeterminate) ===")
        for l in lines:
            print(l)
        allok &= ok
    print("LEMMA 6.2 VALUATION CHECKS " + ("PASSED" if allok else "FAILED"))
