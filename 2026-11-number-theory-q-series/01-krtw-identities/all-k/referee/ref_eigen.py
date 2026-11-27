"""Referee checks (ii), (iii): eigenvectors of T and T# computed in POSITION SPACE (triangular solve of
T phi = lambda phi pointwise, no Bailey-basis recursion), exact at rational (q, x).  Checks:
  * H' phi = eps phi, G phi = g phi;  H'# phi# = eps phi#, G# phi# = g phi#   (Lemma 3.2, Sec. 5)
  * the slice-0 Bailey coefficients of phi_{N,M} vanish beyond L+2M and sum to (xq;q)_{2M} V_{N,M}  (Cor 4.2)
  * the slice-R=0 base-Q Bailey coefficients of phi#_{N,M} vanish beyond M+L and sum to V#_{N,M}     (Cor 5.2)
  * shift (3.2): phi_{N,M}(R,S+1;x) = phi_{N-1,M-1}(R,S;xq^2)
  * biorthogonality psi_a(phi_b) (exact finite sums), Lemma 6.1
  * omega(phi) = lambda W/(xq)_inf (Lemma 4.3) and psi(u) = lambda V#/(x^2Q;Q)_inf (Lemma 5.3), numerically
  * value at infinity: (q)_inf (xq)_inf phi(R,0) -> W as R -> infinity (numerically)
"""
import sys, itertools
from flint import fmpq
import mpmath as mp
from ref_common import Point, poch

mp.mp.dps = 400


def eig_T(P, N, M, Rmax):
    """phi_{N,M} on [0,Rmax] x [0,M], normalised phi(L,M) = 1/(xq^{2M+1};q)_{2L}."""
    L = N - M
    q, x, Q = P.q, P.x, P.Q
    lam = x**(N + M) * q**(N * N + M * M)
    phi = {}
    for S in range(M, -1, -1):
        for R in range(0, Rmax + 1):
            if S == M and R < L:
                phi[(R, S)] = fmpq(0)
                continue
            if (R, S) == (L, M):
                phi[(R, S)] = 1 / poch(x * q**(2 * M + 1), q, 2 * L)
                continue
            rhs = fmpq(0)
            for r in range(R):
                rhs += P.w(r, S) * phi[(r, S)] * P.iq(R - r)
            for s in range(S + 1, M + 1):
                for r in range(R + 1):
                    rhs += P.w(r, s) * phi[(r, s)] * P.iq(R - r) * P.iQ(s - S)
            phi[(R, S)] = rhs / (lam - P.w(R, S))
    return phi


def eig_Tsharp(P, N, M, Smax):
    """phi#_{N,M} on [0,L] x [0,Smax], normalised phi#(L,M) = 1/(x^2 q^{2L+2};Q)_{2M}."""
    L = N - M
    q, x, Q = P.q, P.x, P.Q
    lam = x**(N + M) * q**(N * N + M * M)
    ph = {}
    for R in range(L, -1, -1):
        for S in range(0, Smax + 1):
            if R == L and S < M:
                ph[(R, S)] = fmpq(0)
                continue
            if (R, S) == (L, M):
                ph[(R, S)] = 1 / poch(x**2 * q**(2 * L + 2), Q, 2 * M)
                continue
            rhs = fmpq(0)
            for s in range(S):
                rhs += P.w(R, s) * ph[(R, s)] * P.iQ(S - s)
            for r in range(R + 1, L + 1):
                for s in range(S + 1):
                    rhs += P.w(r, s) * ph[(r, s)] * P.iq(r - R) * P.iQ(S - s)
            ph[(R, S)] = rhs / (lam - P.w(R, S))
    return ph


def bailey_invert(vals, a, p, n_max):
    """alpha_n = (1-a p^{2n})/(1-a) sum_{j<=n} (a;p)_{n+j} (-1)^{n-j} p^{C(n-j,2)} beta_j/(p;p)_{n-j}."""
    out = []
    for n in range(n_max + 1):
        s = fmpq(0)
        for j in range(n + 1):
            s += poch(a, p, n + j) * (-1)**(n - j) * p**((n - j) * (n - j - 1) // 2) * vals[j] / poch(p, p, n - j)
        out.append((1 - a * p**(2 * n)) / (1 - a) * s)
    return out


def V_closed(P, N, M):
    q, x, Q = P.q, P.x, P.Q
    return x**(2 * M) * q**(2 * M * N + M * (2 * M - 1)) * poch(-x, q, M) / (
        poch(Q, Q, M) * poch(Q**(N - M + 1), Q, M) * poch(x * q**M, q, M) * poch(x**2 * q**(2 * N), Q, M))


def Vs_closed(P, N, M):
    q, x, Q = P.q, P.x, P.Q
    L = N - M
    return (-x)**L * q**(L * (2 * L - 1) + 2 * M * L) * poch(q, Q, L) / (poch(q, q, L) * poch(q**L, q, L)) * \
        poch(x**2 * Q, Q, L) / (poch(x * q**(2 * M + 1), q, L) * poch(x * q**(2 * M + L), q, L))


def tomp(a):
    return mp.mpf(int(a.p)) / mp.mpf(int(a.q))


def qinf(a, q):
    return mp.qp(a, q)


def run(qv, xv, Nmax=5, Rext=10, Sext=10):
    P = Point(qv, xv)
    q, x, Q = P.q, P.x, P.Q
    eps = lambda N, M: q**(-N - M) * (1 + x * q**(2 * N)) * (1 + x * q**(2 * M))
    gg = lambda N, M: q**(-2 * N) + x**2 * q**(2 * N) + q**(-2 * M) + x**2 * q**(2 * M)
    modes = [(N, M) for N in range(Nmax + 1) for M in range(N + 1)]
    phis, phss = {}, {}
    nfail = 0
    ncheck = 0
    report = []
    for (N, M) in modes:
        L = N - M
        Rmax = L + 2 * M + Rext
        phi = eig_T(P, N, M, Rmax)
        phis[(N, M)] = phi
        # eigen-relations of H', G (on the interior of the box)
        bx = [(R, S) for R in range(Rmax) for S in range(M + 1)]
        full = dict(phi)
        Hphi = P.Hp(full, bx)
        Gphi = P.Gop(full, bx)
        for k in bx:
            ncheck += 2
            if Hphi[k] != eps(N, M) * full[k]:
                nfail += 1
            if Gphi[k] != gg(N, M) * full[k]:
                nfail += 1
        # slice-0 Bailey coefficients relative to a = x
        vals = [phi[(R, 0)] for R in range(Rmax + 1)]
        kap = bailey_invert(vals, x, q, Rmax)
        tail_ok = all(kap[e] == 0 for e in range(Rmax + 1) if e > L + 2 * M or e < L)
        Wnum = sum(kap)
        Wcl = poch(x * q, q, 2 * M) * V_closed(P, N, M)
        ncheck += 2
        if not tail_ok:
            nfail += 1
        if Wnum != Wcl:
            nfail += 1
        report.append(f"  phi_{{{N},{M}}}: Bailey support in [L,L+2M]: {tail_ok}; W == (xq)_2M V: {Wnum == Wcl}")
        # T#
        Smax = M + L + Sext
        ph = eig_Tsharp(P, N, M, Smax)
        phss[(N, M)] = ph
        bx = [(R, S) for R in range(L + 1) for S in range(Smax)]
        Hs = P.Hps(ph, bx)
        Gs = P.Gs(ph, bx)
        for k in bx:
            ncheck += 2
            if Hs[k] != eps(N, M) * ph[k]:
                nfail += 1
            if Gs[k] != gg(N, M) * ph[k]:
                nfail += 1
        vals = [ph[(0, S)] for S in range(Smax + 1)]
        kap = bailey_invert(vals, x**2, Q, Smax)
        tail_ok = all(kap[e] == 0 for e in range(Smax + 1) if e > M + L or e < M)
        Vs = sum(kap)
        ncheck += 2
        if not tail_ok:
            nfail += 1
        if Vs != Vs_closed(P, N, M):
            nfail += 1
        report.append(f"  phi#_{{{N},{M}}}: base-Q Bailey support in [M,M+L]: {tail_ok}; V# == closed form: {Vs == Vs_closed(P, N, M)}")
    # shift (3.2)
    P2 = Point(qv, fmpq(xv) * fmpq(qv)**2)
    for (N, M) in modes:
        if M == 0:
            continue
        L = N - M
        Rmax = L + 2 * M + Rext
        a = phis[(N, M)]
        b = eig_T(P2, N - 1, M - 1, Rmax)
        for R in range(Rmax + 1):
            for S in range(M):
                ncheck += 1
                if a[(R, S + 1)] != b[(R, S)]:
                    nfail += 1
    # biorthogonality (exact: finite sums)
    lam = lambda N, M: x**(N + M) * q**(N * N + M * M)
    nb = 0
    for A in modes:
        for B in modes:
            LA, MA = A[0] - A[1], A[1]
            LB, MB = B[0] - B[1], B[1]
            s = fmpq(0)
            for r in range(LB, LA + 1):
                for t in range(MA, MB + 1):
                    s += P.w(r, t) * phss[A][(r, t)] * phis[B][(r, t)]
            exp = fmpq(0)
            if A == B:
                exp = lam(*A) / (poch(x**2 * q**(2 * LA + 2), Q, 2 * MA) * poch(x * q**(2 * MA + 1), q, 2 * LA))
            ncheck += 1
            nb += 1
            if s != exp:
                nfail += 1
    report.append(f"  biorthogonality: {nb} exact pairings checked")
    # numerical: omega(phi) and psi(u), and the R->infinity value
    qm, xm = tomp(q), tomp(x)
    xqinf = qinf(xm * qm, qm)
    x2Qinf = qinf(xm**2 * qm**2, qm**2)
    worst = mp.mpf(0)
    worst_l = mp.mpf(0)
    for (N, M) in modes:
        L = N - M
        phi = phis[(N, M)]
        Rmax = L + 2 * M + Rext
        om = fmpq(0)
        for (r, s), v in phi.items():
            om += P.w(r, s) * v * P.iQ(s)
        om = tomp(om)
        pred = tomp(lam(N, M) * poch(x * q, q, 2 * M) * V_closed(P, N, M)) / xqinf
        worst = max(worst, abs(om - pred) / abs(pred))
        ph = phss[(N, M)]
        ps = fmpq(0)
        for (r, s), v in ph.items():
            ps += P.w(r, s) * v * P.iq(r)
        ps = tomp(ps)
        pred = tomp(lam(N, M) * Vs_closed(P, N, M)) / x2Qinf
        worst = max(worst, abs(ps - pred) / abs(pred))
        # value at infinity of slice 0
        lim = qinf(qm, qm) * xqinf * tomp(phi[(Rmax, 0)])
        W = tomp(poch(x * q, q, 2 * M) * V_closed(P, N, M))
        worst_lim = abs(lim - W) / abs(W)
        worst_l = max(worst_l, worst_lim)
    report.append(f"  (q)_inf (xq)_inf phi(Rmax,0) vs W (Rmax = L+2M+{Rext}): worst rel. err {mp.nstr(worst_l, 3)}")
    report.append(f"  omega(phi)=lam W/(xq)_inf and psi(u)=lam V#/(x^2Q;Q)_inf: worst rel. err {mp.nstr(worst, 3)}")
    return ncheck, nfail, report


if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for (qv, xv) in ((fmpq(2, 7), fmpq(3, 11)), (fmpq(-3, 5), fmpq(5, 13))):
        n, f, rep = run(qv, xv, Nmax=Nmax)
        print(f"(q,x)=({qv},{xv}), modes N<={Nmax}: {n} exact checks, {f} failures")
        for line in rep:
            print(line)
