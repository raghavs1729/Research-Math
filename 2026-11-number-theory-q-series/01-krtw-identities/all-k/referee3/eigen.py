#!/usr/bin/env python3
"""Referee 3, check (4): eigenvectors of T and of the reversed chain T#, computed by a pointwise triangular
solve in POSITION space (not in the Bailey basis), exact rational arithmetic at rational (q, x).

phi_{N,M}: the solution of T phi = lam phi supported in S <= M, normalised by phi(L,M) = 1/(x q^{2M+1};q)_{2L}
(the value of Xi_{L,M} at (L,M)).  Order of the solve: S = M, M-1, ..., 0; R = 0, 1, 2, ...;
    (lam - w(R,S)) phi(R,S) = sum_{r<=R, S<=s<=M, (r,s) != (R,S)} w(r,s) phi(r,s) / ((q)_{R-r} (Q;Q)_{s-S}).
phi#_{N,M}: T# G = lam G supported in R <= L, normalised by G(L,M) = 1/(x^2 q^{2L+2};Q)_{2M}; order R = L..0, S = 0,1,..
Then: Bailey-basis coefficients (by my own triangular inversion), supports, values at infinity v_a
(coefficient sums), and the claims of F-all-k Sec. 3-6 and iv-all-k Lemma 3.1.
"""
import sys
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 150

POINTS = [(Fr(2, 7), Fr(3, 5)), (Fr(-1, 3), Fr(2, 11))]
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4


class Ctx:
    def __init__(self, q, x):
        self.q, self.x, self.Q = q, x, q * q

    def poch(self, a, p, n):
        r = Fr(1)
        for i in range(n):
            r *= 1 - a * p ** i
        return r

    def qq(self, n):
        return self.poch(self.q, self.q, n)

    def QQ(self, n):
        return self.poch(self.Q, self.Q, n)

    def w(self, r, s):
        q, x = self.q, self.x
        return x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s)

    def lam(self, N, M):
        return self.x ** (N + M) * self.q ** (N * N + M * M)

    def eps(self, N, M):
        q, x = self.q, self.x
        return q ** (-N - M) * (1 + x * q ** (2 * N)) * (1 + x * q ** (2 * M))

    def g(self, N, M):
        q, x = self.q, self.x
        return q ** (-2 * N) + x * x * q ** (2 * N) + q ** (-2 * M) + x * x * q ** (2 * M)

    def phi_e(self, e, R, a):
        if R < e:
            return Fr(0)
        return 1 / (self.qq(R - e) * self.poch(a * self.q, self.q, R + e))

    def phiQ_e(self, e, S, a):
        if S < e:
            return Fr(0)
        return 1 / (self.QQ(S - e) * self.poch(a * self.Q, self.Q, S + e))

    # ---------- eigenvectors ----------
    def eig(self, N, M, Rb):
        L = N - M
        lam = self.lam(N, M)
        ph = {}
        for S in range(M, -1, -1):
            for R in range(Rb + 1):
                rhs = Fr(0)
                for s in range(S, M + 1):
                    for r in range(R + 1):
                        if (r, s) == (R, S):
                            continue
                        v = ph.get((r, s), 0)
                        if v:
                            rhs += self.w(r, s) * v / (self.qq(R - r) * self.QQ(s - S))
                if (R, S) == (L, M):
                    assert rhs == 0, "inconsistent top equation"
                    ph[(R, S)] = 1 / self.poch(self.x * self.q ** (2 * M + 1), self.q, 2 * L)
                else:
                    d = lam - self.w(R, S)
                    assert d != 0
                    ph[(R, S)] = rhs / d
        return ph

    def eig_sharp(self, N, M, Sb):
        L = N - M
        lam = self.lam(N, M)
        G = {}
        for r in range(L, -1, -1):
            for s in range(Sb + 1):
                rhs = Fr(0)
                for R in range(r, L + 1):
                    for S in range(s + 1):
                        if (R, S) == (r, s):
                            continue
                        v = G.get((R, S), 0)
                        if v:
                            rhs += self.w(R, S) * v / (self.qq(R - r) * self.QQ(s - S))
                if (r, s) == (L, M):
                    assert rhs == 0
                    G[(r, s)] = 1 / self.poch(self.x * self.x * self.q ** (2 * L + 2), self.Q, 2 * M)
                else:
                    d = lam - self.w(r, s)
                    assert d != 0
                    G[(r, s)] = rhs / d
        return G

    # ---------- Bailey-basis expansions (my own triangular inversion) ----------
    def expand(self, f, a, Rb):
        """f: list f[0..Rb]; returns c[0..Rb] with f(R) = sum_e c_e phi_e(R;a) for R <= Rb"""
        c = []
        for R in range(Rb + 1):
            v = f[R] - sum(c[e] * self.phi_e(e, R, a) for e in range(R))
            c.append(v / self.phi_e(R, R, a))
        return c

    def expandQ(self, f, a, Sb):
        c = []
        for S in range(Sb + 1):
            v = f[S] - sum(c[e] * self.phiQ_e(e, S, a) for e in range(S))
            c.append(v / self.phiQ_e(S, S, a))
        return c

    # ---------- closed forms claimed in F-all-k ----------
    def V(self, N, M, x=None):
        q = self.q
        x = self.x if x is None else x
        Q = q * q
        num = x ** (2 * M) * q ** (2 * M * N + M * (2 * M - 1)) * self.poch(-x, q, M)
        den = self.poch(Q, Q, M) * self.poch(Q ** (N - M + 1), Q, M) * self.poch(x * q ** M, q, M) \
            * self.poch(x * x * q ** (2 * N), Q, M)
        return num / den

    def Vs(self, N, M):
        q, x, Q = self.q, self.x, self.Q
        L = N - M
        r = (-x) ** L * q ** (L * (2 * L - 1) + 2 * M * L)
        r *= self.poch(q, q * q, L) / (self.qq(L) * self.poch(q ** L, q, L))
        r *= self.poch(x * x * Q, Q, L) / (self.poch(x * q ** (2 * M + 1), q, L) * self.poch(x * q ** (2 * M + L), q, L))
        return r


def mpf(v):
    return mp.mpf(v.numerator) / v.denominator


def poch_inf_mp(a, p):
    a, p = mpf(a), mpf(p)
    r = mp.mpf(1)
    t = a
    while abs(t) > mp.mpf(10) ** (-170):
        r *= 1 - t
        t *= p
    return r


def run(q, x):
    C = Ctx(q, x)
    Q = C.Q
    out = []
    fails = []
    modes = [(N, M) for N in range(NMAX + 1) for M in range(N + 1)]
    Rb_big = 3 * NMAX + 6
    eigs, eigs_s = {}, {}
    for (N, M) in modes:
        L = N - M
        Rb = max(L + 2 * M + 3, 2 * NMAX + 3)
        ph = C.eig(N, M, Rb_big)
        eigs[(N, M)] = ph
        # support
        sup_ok = all(ph[(R, S)] == 0 for (R, S) in ph if R < L)
        # Bailey-basis coefficients per slice, support L <= e <= L + 2(M-S)
        coef_ok = True
        slices = {}
        for S in range(M + 1):
            f = [ph[(R, S)] for R in range(Rb_big + 1)]
            c = C.expand(f, x * q ** (2 * S), Rb_big)
            slices[S] = c
            for e, v in enumerate(c):
                inside = L <= e <= L + 2 * (M - S)
                if (v != 0) != inside and not (inside and v == 0):
                    coef_ok = False
                if not inside and v != 0:
                    coef_ok = False
            if S == M and (c[L] != 1 or any(c[e] for e in range(len(c)) if e != L)):
                coef_ok = False
        W = sum(slices[0])
        W_claim = C.poch(x * q, q, 2 * M) * C.V(N, M)
        # eigen-relations of H', G pointwise
        eps, g = C.eps(N, M), C.g(N, M)

        def get(R, S):
            return ph.get((R, S), Fr(0)) if R >= 0 and S >= 0 else Fr(0)
        hg_ok = True
        for R in range(Rb_big):
            for S in range(M + 1):
                Hv = (q ** (-2 * S) + x) * q ** (-R) * (get(R, S) - get(R - 1, S)) \
                    + x * q ** R * (1 + x * q ** (2 * S)) * get(R, S) - x * x * q ** (R + 2 * S) * get(R, S + 1)
                Gv = q ** (-2 * R - 2 * S) * (get(R, S) - (1 + q) * get(R - 1, S) + q * get(R - 2, S)) \
                    + q ** (-2 * S) * get(R, S) - x * (1 + 1 / q) * get(R - 1, S) \
                    + x * x * (1 + q ** (2 * R)) * q ** (2 * S) * (get(R, S) - get(R, S + 1))
                if Hv != eps * get(R, S) or Gv != g * get(R, S):
                    hg_ok = False
        # shift (3.2)
        shift_ok = None
        if M >= 1:
            C2 = Ctx(q, x * q * q)
            ph2 = C2.eig(N - 1, M - 1, Rb_big)
            shift_ok = all(ph[(R, S + 1)] == ph2[(R, S)] for R in range(Rb_big + 1) for S in range(M))
        # slice-0 equations at infinity (Prop 4.1 / iv Lemma 3.1): v_x(h_x F0), v_x(q^R F1)
        F0 = [ph[(R, 0)] for R in range(Rb_big + 1)]
        hF0 = [q ** (-R) * (F0[R] - (F0[R - 1] if R else 0)) + x * q ** R * F0[R] for R in range(Rb_big + 1)]
        mu0 = sum(C.expand(F0, x, Rb_big))
        mu1 = sum(C.expand(hF0, x, Rb_big))
        inf_ok = (1 + x) * mu1 == eps * mu0
        qF1_ok = None
        if M >= 1:
            F1 = [ph[(R, 1)] for R in range(Rb_big + 1)]
            qF1 = [q ** R * F1[R] for R in range(Rb_big + 1)]
            cq = C.expand(qF1, x, Rb_big)
            hqF1 = [q ** (-R) * (qF1[R] - (qF1[R - 1] if R else 0)) + x * q ** R * qF1[R] for R in range(Rb_big + 1)]
            Cs = Ctx(q, x * q * q)
            W_prev = sum(Cs.expand([Cs.eig(N - 1, M - 1, Rb_big)[(R, 0)] for R in range(Rb_big + 1)], x * q * q, Rb_big))
            v1 = C.poch(x * q, q, 2) * W_prev
            onepq = [(1 + q ** (2 * R)) * F1[R] for R in range(Rb_big + 1)]
            qF1_ok = (sum(cq) == 0 and sum(C.expand(hqF1, x, Rb_big)) == (1 - 1 / q) * v1
                      and sum(C.expand(onepq, x, Rb_big)) == v1)
            # Prop 4.1 recursion
            rec = x * x * (1 + x) * (x + 1 / q) * q ** (2 * N + 2 * M) * C.poch(x * q, q, 2) \
                / ((1 - q ** (2 * N)) * (1 - q ** (2 * M)) * (1 - x * x * q ** (2 * N)) * (1 - x * x * q ** (2 * M))) * W_prev
            if rec != W:
                fails.append(f"Prop 4.1 recursion fails at {(N, M)}")
        # Lemma 4.3 and iv Lemma 3.1 numerically: omega(phi) = lam W/(xq)_inf, omega(h phi) = eps/(1+x) omega(phi)
        om_ex = Fr(0)
        omh_ex = Fr(0)
        for (R, S), v in ph.items():
            if R == Rb_big:
                continue          # h phi at R needs phi(R-1): fine; but keep the same index set for both sums
            hv = q ** (-R) * (v - get(R - 1, S)) + x * q ** (R + 2 * S) * v
            om_ex += C.w(R, S) * v / C.QQ(S)
            omh_ex += C.w(R, S) * hv / C.QQ(S)
        om, omh = mpf(om_ex), mpf(omh_ex)
        om_claim = mpf(C.lam(N, M) * W) / poch_inf_mp(x * q, q)
        l43 = abs(om - om_claim) / max(abs(om_claim), mp.mpf(10) ** -140)
        l31 = abs(omh - mpf(eps / (1 + x)) * om) / max(abs(om), mp.mpf(10) ** -140)
        # reversed chain
        Sb = max(M + L + 3, 2 * NMAX + 3)
        Sb_big = 3 * NMAX + 6
        G = C.eig_sharp(N, M, Sb_big)
        eigs_s[(N, M)] = G
        supS_ok = all(G[(r, s)] == 0 for (r, s) in G if s < M)
        vs = {}
        coefS_ok = True
        for rho in range(L + 1):
            f = [G[(rho, s)] for s in range(Sb_big + 1)]
            c = C.expandQ(f, x * x * q ** (2 * rho), Sb_big)
            for e, v in enumerate(c):
                inside = M <= e <= M + L - rho
                if not inside and v != 0:
                    coefS_ok = False
            vs[rho] = sum(c)
        rec_ok = vs[L] == 1
        for rho in range(L):
            fac = -x * q ** (L + N + M - 1) * (1 - x * x * q ** (2 * rho + 2)) * (1 - q ** (2 * rho + 1)) \
                / ((1 - q ** (L - rho)) * (1 - q ** (L + rho)) * (1 - x * q ** (N + M - rho)) * (1 - x * q ** (N + M + rho)))
            if vs[rho] != fac * vs[rho + 1]:
                rec_ok = False
        Vs_ok = vs[0] == C.Vs(N, M)

        def gS(R, S):
            return G.get((R, S), Fr(0)) if R >= 0 and S >= 0 else Fr(0)
        hgS_ok = True
        for R in range(L + 1):
            for S in range(Sb_big):
                Hs = q ** (-R - 2 * S) * (gS(R, S) - gS(R, S - 1)) + x * q ** (-R) * gS(R, S) \
                    + x * q ** R * (1 + x * q ** (2 * S)) * (gS(R, S) - gS(R + 1, S))
                a = x * x * q ** (2 * R)
                hq = q ** (-2 * S) * (gS(R, S) - gS(R, S - 1)) + a * q ** (2 * S) * gS(R, S)
                Gsv = (q ** (-2 * R) + 1) * hq - x * (1 + 1 / q) * gS(R + 1, S) \
                    - x * x * (1 + q) * q ** (2 * R + 2 * S) * gS(R + 1, S) + q * x * x * q ** (2 * R + 2 * S) * gS(R + 2, S)
                if Hs != eps * gS(R, S) or Gsv != g * gS(R, S):
                    hgS_ok = False
        # Lemma 5.3 numerically
        psiu = mpf(sum(C.w(r, s) * v / C.qq(r) for (r, s), v in G.items()))
        psiu_claim = mpf(C.lam(N, M) * C.Vs(N, M)) / poch_inf_mp(x * x * Q, Q)
        l53 = abs(psiu - psiu_claim) / max(abs(psiu_claim), mp.mpf(10) ** -140)
        line = (f"(N,M)=({N},{M}): supp phi {sup_ok}, Bailey coeff. support {coef_ok}, H'phi=eps phi & G phi=g phi {hg_ok},"
                f" shift(3.2) {shift_ok}, W=(xq)_2M V [Cor4.2] {W == W_claim}, (1+x)v(h F0)=eps v(F0) {inf_ok},"
                f" v(q^R F1)=0 etc. {qF1_ok}, Lemma4.3 rel.err {mp.nstr(l43, 3)}, ivLemma3.1 rel.err {mp.nstr(l31, 3)};"
                f" supp phi# {supS_ok}, coeff. support {coefS_ok}, Prop5.1 rec {rec_ok}, v0=V# [Cor5.2] {Vs_ok},"
                f" H'#,G# eigen {hgS_ok}, Lemma5.3 rel.err {mp.nstr(l53, 3)}")
        out.append(line)
        good = sup_ok and coef_ok and hg_ok and (shift_ok is not False) and W == W_claim and inf_ok \
            and (qF1_ok is not False) and l43 < 1e-90 and l31 < 1e-90 and supS_ok and coefS_ok and rec_ok and Vs_ok \
            and hgS_ok and l53 < 1e-90
        if not good:
            fails.append(f"mode {(N, M)}: {line}")
    # Lemma 6.1 biorthogonality, all pairs
    pair_ok = True
    for a in modes:
        for b in modes:
            La, Ma = a[0] - a[1], a[1]
            G = eigs_s[a]
            ph = eigs[b]
            val = Fr(0)
            for (r, s), gv in G.items():
                pv = ph.get((r, s))
                if pv:
                    val += C.w(r, s) * gv * pv
            if a != b:
                if val != 0:
                    pair_ok = False
                    fails.append(f"psi_{a}(phi_{b}) = {val} != 0")
            else:
                N, M = a
                L = N - M
                claim = C.lam(N, M) / (C.poch(x * x * q ** (2 * L + 2), Q, 2 * M) * C.poch(x * q ** (2 * M + 1), q, 2 * L))
                if val != claim:
                    pair_ok = False
                    fails.append(f"psi_{a}(phi_{a}) = {val} != claim {claim}")
    out.append(f"Lemma 6.1 (all {len(modes)}^2 pairs, finite exact sums over the overlap of supports): {pair_ok}")
    return out, fails


if __name__ == "__main__":
    allfails = []
    for (q, x) in POINTS:
        print(f"=== (q, x) = ({q}, {x}), modes N <= {NMAX} ===", flush=True)
        out, fails = run(q, x)
        for l in out:
            print(l, flush=True)
        allfails += fails
    print("ALL EIGENVECTOR CHECKS " + ("PASSED" if not allfails else "FAILED:\n" + "\n".join(allfails)))
