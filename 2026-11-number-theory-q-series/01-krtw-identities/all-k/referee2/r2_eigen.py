"""Referee 2, computation (4): Lemma 3.1 on eigenvectors of T computed in POSITION space (no Bailey recursion).

At a rational point (q, x), phi = phi_{N,M} is obtained from T phi = lam phi, lam = x^{N+M} q^{N^2+M^2}, by the pointwise
triangular solve (slices S = M, M-1, ..., 0; R = 0, 1, ...):
    (lam - w(R,S)) phi(R,S) = sum_{(r,s) != (R,S), r<=R, S<=s<=M} w(r,s) phi(r,s)/((q)_{R-r}(Q;Q)_{s-S}),
phi = 0 for S > M, normalisation phi(L,M) = 1/(xq^{2M+1};q)_{2L} (top coefficient Xi_{L,M} = 1).  Exact Fractions.
Checks (exact unless marked num):
  E0  top slice = phi_L(R; xq^{2M});  H' phi = eps phi on the box (F-all-k Lemma 3.2, consistency of the eigenvectors)
  L31a  h Xi_{e,s} = (q^{-e} + x q^{2s+e}) Xi_{e,s}
  L21e  T(h phi) = lam (h phi - x O phi) on the box
  Bailey coefficients (relative to x) of F_0 = phi(.,0), h_x F_0, q^R F_1, (T h phi)(.,0) by a triangular solve on
  R = 0..E, E = L+2M, VERIFIED on R = E+1..Rmax (so the expansions are finite as claimed); v_x = coefficient sum:
  L31c1 v_x(q^R F_1) = 0 (M >= 1);  L31c2 (1+x) v_x(h_x F_0) = eps v_x(F_0);  L31c3 v_x((T h phi)(.,0)) = lam eps v_x(F_0)/(1+x)
  W     v_x(F_0) = W_{N,M} of F-all-k Cor 4.2 (consistency)
  num: Lemma 3.1  omega(h phi) = eps/(1+x) omega(phi); also omega(P^-1 phi); (b) (xq)_inf omega(h phi) = v_x((T h phi)(.,0));
       Lemma 4.3  (xq)_inf omega(phi) = lam v_x(F_0).  omega truncated at R <= Rmax (tail ~ |q|^{(Rmax+1)^2}).
  negative control: omega(h phi) vs eps_{other}/(1+x) omega(phi).
Usage: python3 r2_eigen.py [Nmax=4] [Rmax=22]
"""
import sys, time
from fractions import Fraction as Fr
import mpmath

mpmath.mp.dps = 420
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
Rmax = int(sys.argv[2]) if len(sys.argv) > 2 else 22


def run_point(q, x, out):
    Q = q * q
    allok = True

    def poch(a, base, n):
        r = Fr(1)
        for j in range(n):
            r *= (1 - a * base ** j)
        return r

    iq = {n: 1 / poch(q, q, n) for n in range(0, Rmax + 40)}
    iQ = {n: 1 / poch(Q, Q, n) for n in range(0, Rmax + 40)}

    def w(r, s):
        return x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s)

    def phiB(e, R, a):  # Bailey basis function phi_e(R; a)
        if R < e:
            return Fr(0)
        return iq[R - e] / poch(a * q, q, R + e)

    def eigvec(N, M):
        Lr = N - M
        lam = x ** (N + M) * q ** (N * N + M * M)
        phi = {}
        for S in range(M, -1, -1):
            for R in range(0, Rmax + 1):
                if (R, S) == (Lr, M):
                    phi[(R, S)] = 1 / poch(x * q ** (2 * M + 1), q, 2 * Lr)
                    continue
                tot = Fr(0)
                for s in range(S, M + 1):
                    for r in range(0, R + 1):
                        if (r, s) == (R, S):
                            continue
                        v = phi[(r, s)]
                        if v:
                            tot += w(r, s) * v * iq[R - r] * iQ[s - S]
                d = lam - w(R, S)
                assert d != 0
                phi[(R, S)] = tot / d
        return lam, phi

    def get(F, R, S):
        return F.get((R, S), Fr(0)) if R >= 0 and S >= 0 else Fr(0)

    def h(F, M):
        return {(R, S): q ** (-R) * (get(F, R, S) - get(F, R - 1, S)) + x * q ** (R + 2 * S) * get(F, R, S)
                for R in range(Rmax + 1) for S in range(M + 1)}

    def Oop(F, M):
        return {(R, S): q ** (R + 2 * S) * get(F, R, S + 1) for R in range(Rmax + 1) for S in range(M + 1)}

    def T(F, M):
        G = {}
        for R in range(Rmax + 1):
            for S in range(M + 1):
                tot = Fr(0)
                for s in range(S, M + 1):
                    for r in range(R + 1):
                        v = get(F, r, s)
                        if v:
                            tot += w(r, s) * v * iq[R - r] * iQ[s - S]
                G[(R, S)] = tot
        return G

    def Hprime(F, M):  # F-all-k: H' = (Y^-1 + x) h_{xY} - x^2 P Y u
        hF = h(F, M)
        return {(R, S): (q ** (-2 * S) + x) * hF[(R, S)] - x * x * q ** (R + 2 * S) * get(F, R, S + 1)
                for R in range(Rmax + 1) for S in range(M + 1)}

    def bailey_coeffs(f, E):
        """f(R) = sum_{e<=E} kappa_e phi_e(R;x): solve on R = 0..E, verify on R = E+1..Rmax."""
        kap = []
        for R in range(E + 1):
            acc = f[R] - sum(kap[e] * phiB(e, R, x) for e in range(R))
            kap.append(acc / phiB(R, R, x))
        ok = all(f[R] == sum(kap[e] * phiB(e, R, x) for e in range(E + 1)) for R in range(E + 1, Rmax + 1))
        return kap, ok

    def eps(N, M):
        return q ** (-N - M) * (1 + x * q ** (2 * N)) * (1 + x * q ** (2 * M))

    def Wformula(N, M):  # F-all-k Cor 4.2: W = (xq;q)_{2M} V_{N,M}
        V = x ** (2 * M) * q ** (2 * M * N + M * (2 * M - 1)) * poch(-x, q, M) / (
            poch(Q, Q, M) * poch(Q ** (N - M + 1), Q, M) * poch(x * q ** M, q, M) * poch(x * x * q ** (2 * N), Q, M))
        return poch(x * q, q, 2 * M) * V

    def tomp(fr):
        return mpmath.mpf(fr.numerator) / mpmath.mpf(fr.denominator)

    xq_inf = mpmath.mpf(1)
    xm, qm = tomp(x), tomp(q)
    j = 1
    while True:
        t = xm * qm ** j
        if abs(t) < mpmath.mpf(10) ** (-430):
            break
        xq_inf *= (1 - t)
        j += 1

    # L31a: h diagonal on Xi
    for e in range(4):
        for sg in range(3):
            Xi = {(R, S): (phiB(e, R, x * q ** (2 * sg)) if S == sg else Fr(0)) for R in range(Rmax + 1) for S in range(3)}
            hX = h(Xi, 2)
            c = q ** (-e) + x * q ** (2 * sg + e)
            if any(hX[p] != c * Xi[p] for p in hX):
                allok = False
                out("    L31a FAIL e=%d sigma=%d" % (e, sg))
    out("    L31a h Xi_{e,s} = (q^-e + x q^{2s+e}) Xi_{e,s}, e<4, s<3: %s" % ("OK" if allok else "FAIL"))

    worst = mpmath.mpf(0)
    worstb = mpmath.mpf(0)
    for N in range(Nmax + 1):
        for M in range(N + 1):
            t0 = time.time()
            Lr = N - M
            lam, phi = eigvec(N, M)
            ep = eps(N, M)
            msgs = []
            # E0
            top_ok = all(phi[(R, M)] == phiB(Lr, R, x * q ** (2 * M)) for R in range(Rmax + 1))
            Hp = Hprime(phi, M)
            # H'phi on S<=M; compare where all inputs are inside the box (R <= Rmax), all fine
            h_ok = all(Hp[(R, S)] == ep * phi[(R, S)] for R in range(Rmax + 1) for S in range(M + 1))
            hphi = h(phi, M)
            Ophi = Oop(phi, M)
            Thphi = T(hphi, M)
            l21_ok = all(Thphi[p] == lam * (hphi[p] - x * Ophi[p]) for p in Thphi)
            E = Lr + 2 * M
            F0 = [phi[(R, 0)] for R in range(Rmax + 1)]
            hF0 = [hphi[(R, 0)] for R in range(Rmax + 1)]
            qF1 = [q ** R * get(phi, R, 1) for R in range(Rmax + 1)]
            ThF = [Thphi[(R, 0)] for R in range(Rmax + 1)]
            k0, ok0 = bailey_coeffs(F0, E)
            k1, ok1 = bailey_coeffs(hF0, E)
            k2, ok2 = bailey_coeffs(qF1, E)
            k3, ok3 = bailey_coeffs(ThF, E)
            v0, v1, v2, v3 = sum(k0), sum(k1), sum(k2), sum(k3)
            supp0 = [e for e in range(E + 1) if k0[e] != 0]
            c1 = (v2 == 0)
            c2 = ((1 + x) * v1 == ep * v0)
            c3 = (v3 == lam * ep * v0 / (1 + x))
            cW = (v0 == Wformula(N, M))
            exact_ok = top_ok and h_ok and l21_ok and ok0 and ok1 and ok2 and ok3 and c1 and c2 and c3 and cW
            allok &= exact_ok
            # numerics
            om = sum(w(R, S) * phi[(R, S)] * iQ[S] for R in range(Rmax + 1) for S in range(M + 1))
            omh = sum(w(R, S) * hphi[(R, S)] * iQ[S] for R in range(Rmax + 1) for S in range(M + 1))
            omP = sum(w(R, S) * q ** (-R) * phi[(R, S)] * iQ[S] for R in range(Rmax + 1) for S in range(M + 1))
            omm, omhm, omPm = tomp(om), tomp(omh), tomp(omP)
            epm = tomp(ep)
            rel1 = abs(omhm - epm / (1 + xm) * omm) / abs(omm)
            rel2 = abs(omPm - epm / (1 + xm) * omm) / abs(omm)
            rel3 = abs(xq_inf * omhm - tomp(v3)) / abs(tomp(v3))
            rel4 = abs(xq_inf * omm - tomp(lam * v0)) / abs(tomp(lam * v0))
            Nn, Mn = N + 1, M  # a different mode
            relneg = abs(omhm - tomp(eps(Nn, Mn)) / (1 + xm) * omm) / abs(omm)
            worst = max(worst, rel1, rel2, rel3, rel4)
            if relneg < mpmath.mpf(10) ** (-5):
                allok = False
                msgs.append("NEG CONTROL NOT DETECTED")
            out("    (N,M)=(%d,%d): top/H'/L2.1/expansions %s; Bailey supp(F0)=%s; v_x(q^R F1)=0 %s; (1+x)v(h F0)=eps v(F0) %s; "
                "v(T h phi)=lam eps W/(1+x) %s; W=Cor4.2 %s | num rel.err L3.1 %.1e, P^-1 %.1e, (b) %.1e, L4.3 %.1e; neg %.1e  (%.1fs) %s"
                % (N, M, "OK" if (top_ok and h_ok and l21_ok and ok0 and ok1 and ok2 and ok3) else "FAIL",
                   supp0, "OK" if c1 else "FAIL", "OK" if c2 else "FAIL", "OK" if c3 else "FAIL", "OK" if cW else "FAIL",
                   rel1,
                   rel2, rel3, rel4, relneg, time.time() - t0, " ".join(msgs)))
    return allok, worst


if __name__ == '__main__':
    lines = []
    pr = lambda s: (print(s, flush=True), lines.append(s))
    pr("Referee 2, (4): Lemma 3.1 on position-space eigenvectors, N <= %d, R <= %d (exact Fractions; numerics mp.dps=420)"
       % (Nmax, Rmax))
    ok = True
    for (q, x) in [(Fr(2, 7), Fr(3, 11)), (Fr(-3, 5), Fr(5, 13)), (Fr(1, 3), Fr(-2, 9))]:
        pr("  point (q,x) = (%s, %s)" % (q, x))
        a, wst = run_point(q, x, pr)
        ok &= a
        pr("  worst relative error of the numerical checks at this point: %.2e" % wst)
    pr("RESULT: %s" % ("ALL EXACT CHECKS OK" if ok else "SOME CHECK FAILED"))
