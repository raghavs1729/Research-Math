"""Certificate for proofs/iv-all-k.md (KRTW identity (iv) for all k, and the refinement K^Om_k = B^(k)_{k+1} J_k).

Every check prints OK/FAIL; the script ends with ALL CHECKS PASSED or a failure count.
  D1  [T, h] = -x O T on point masses (exact, two rational points), h = P^-1(1-S) + x Y P, O = P Y up
  D2  omega(h F) = omega(P^-1 F) on point masses (exact)
  D3  K^Om_k(x) = x^-1 (K^b_k(x) - K_k(x))   (direct summations, all powers of x, mod q^QM, k <= 6);
      K^Om_k summed literally with the 2k-1 signed insertions of Omega, K^b_k = K_k with the extra factor q^{-r_1}
  D4  omega(h phi_{N,M}) = eps_{N,M}/(1+x) omega(phi_{N,M}) for N <= 5, via v_x((T h phi)(.,0)) computed from the
      eigenvectors of the triangular system (tower.py), exact in Q(x,q)
  D5  algebra: t (eps-1-x)/(x(1+x)) = lam_N lam_M (bO_N bJ_M + bO_M bJ_N)/(1+delta), N <= 12 (symbolic);
      bO_n := betaB_n (q^{-2n} + x - 1 + x^2 q^{2n})/(x(1+x));  and lam_n^2 bO_n = d_n (the alpha* of the change of a)
  D6  spectral side = direct side: (xq)_inf (x^2Q;Q)_inf K^Om_k = sum_{N>=M} lam^k t (eps-1-x)/(x(1+x)), k <= 6, mod q^QM
  D7  (xq)_inf B^(k)_{k+1}(x) = sum_n lam_n^{k+1} bO_n  (direct summation of B^(k)_{k+1}), k <= 6, mod q^QM
  D8  x = 1: sum_n lam_n^{k+1} bO_n |_{x=1} = (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_inf mod q^200, k <= 10 (integer series)
Negative controls: D3 with Omega's last term dropped; D6 with eps replaced by eps + x(1+x)q (must FAIL -> reported as detected).
"""
import sys
from fractions import Fraction as Fr
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from tower import lam
from check_II import target, betaB, betaJ
import verify_ops as VO
from certify_F_all_k import BS, bs_from_rf, phi_coeffs, ev_rf

FAILS = 0


def report(label, ok, extra=''):
    global FAILS
    print('[%s] %s %s' % (label, 'OK' if ok else 'FAIL', extra))
    sys.stdout.flush()
    if not ok:
        FAILS += 1


def eps(N, M):
    return mono(0, -N - M) * RF(ONE + X * Q ** (2 * N)) * RF(ONE + X * Q ** (2 * M))


def lamn(n):
    return mono(n, n * n)


def bO(n):
    """betaB_n (q^{-2n} + x - 1 + x^2 q^{2n}) / (x (1+x))"""
    return betaB(n) * (mono(0, -2 * n) + RF(X) - RF(1) + mono(2, 2 * n)) / (RF(X) * RF(ONE + X))


def dstar(e):
    """alpha*_e of the change of a (relative to x) applied to (L-inf)^2 of the pair beta_n = 1/(Q;Q)_n rel. to xq"""
    s = -1 if e % 2 else 1
    t = RF(s) * mono(2 * e, 3 * e * e + 2 * e) * poch(X * X * Q ** 2, e, 2) / qpoch(e, 2)
    if e >= 1:
        t = t + RF(s) * mono(2 * e - 1, 3 * e * e - 2 * e) * poch(X * X * Q ** 2, e - 1, 2) / qpoch(e - 1, 2)
    return t


# ------------------------------------------------------------------ direct summations (bivariate, mod q^QM)
def inv_poch_bs(n, step, QM):
    s = BS({(0, 0): 1})
    for i in range(1, n + 1):
        g = BS({(0, step * i * m): 1 for m in range((QM - 1) // (step * i) + 1)})
        s = s * g
    return s


class Chain:
    """sums over r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k of prod_i x^{r_i+2s_i} q^{(r_i+s_i)^2+s_i^2} / bonds,
    with an optional insertion: ('diag', i): q^{r_i+2s_i};  ('cross', i): q^{r_i+2s_{i+1}} (i<k);  ('rm1',): q^{-r_1}"""

    def __init__(self, QM):
        self.QM = QM
        self.cq = {}
        self.cQ = {}
        rmax = 0
        while (rmax + 1) ** 2 - (rmax + 1) < QM + 2:
            rmax += 1
        self.rmax = rmax

    def iq(self, n):
        if n not in self.cq:
            self.cq[n] = inv_poch_bs(n, 1, self.QM)
        return self.cq[n]

    def iQ(self, n):
        if n not in self.cQ:
            self.cQ[n] = inv_poch_bs(n, 2, self.QM)
        return self.cQ[n]

    def total(self, k, ins=None):
        QM, rmax = self.QM, self.rmax
        BS.QM = QM
        F = {}
        for lev in range(k, 0, -1):
            newF = {}
            for R in range(rmax + 1):
                for S in range(rmax + 1):
                    e = (R + S) ** 2 + S * S
                    if ins is not None:
                        if ins[0] == 'diag' and ins[1] == lev:
                            e += R + 2 * S
                        if ins[0] == 'cross' and ins[1] == lev:
                            e += R
                        if ins[0] == 'rm1' and lev == 1:
                            e -= R
                    if e >= QM:
                        continue
                    if lev == k:
                        acc = self.iq(R)
                    else:
                        acc = BS()
                        for (r2, s2), val in F.items():
                            if r2 > R or s2 < S:
                                continue
                            t = val * self.iq(R - r2) * self.iQ(s2 - S)
                            if ins is not None and ins[0] == 'cross' and ins[1] == lev:
                                t = t * BS({(0, 2 * s2): 1})
                            acc = acc + t
                    if acc.d:
                        newF[(R, S)] = acc * BS({(R + 2 * S, e): 1})
            F = newF
        tot = BS()
        for (R, S), val in F.items():
            tot = tot + val * self.iQ(S)
        return tot


def direct_Bk(k, QM):
    """B^(k)_{k+1}(x) = sum_{N_1>=...>=N_{k+1}>=0} x^{sum N} q^{sum N^2 + N_k + N_{k+1}} / ((q)_{N_1-N_2}...(q)_{N_k-N_{k+1}} (Q;Q)_{N_{k+1}})"""
    BS.QM = QM
    nmax = 0
    while (nmax + 1) ** 2 < QM:
        nmax += 1
    m = k + 1
    # level m (innermost): G(N) = x^N q^{N^2 + N} / (Q;Q)_N
    G = {}
    for N in range(nmax + 1):
        e = N * N + N
        if e < QM:
            G[N] = BS({(N, e): 1}) * inv_poch_bs(N, 2, QM)
    for lev in range(m - 1, 0, -1):
        newG = {}
        for N in range(nmax + 1):
            e = N * N + (N if lev >= k else 0)
            if e >= QM:
                continue
            acc = BS()
            for N2, val in G.items():
                if N2 <= N:
                    acc = acc + val * inv_poch_bs(N - N2, 1, QM)
            if acc.d:
                newG[N] = acc * BS({(N, e): 1})
        G = newG
    tot = BS()
    for N, val in G.items():
        tot = tot + val
    return tot


def xq_inf(QM):
    P = BS({(0, 0): 1})
    for j in range(1, QM):
        P = P * BS({(0, 0): 1, (1, j): -1})
    return P


def x2Q_inf(QM):
    P = BS({(0, 0): 1})
    for j in range(1, QM // 2 + 1):
        P = P * BS({(0, 0): 1, (2, 2 * j): -1})
    return P


def main():
    QM = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    KM = int(sys.argv[2]) if len(sys.argv) > 2 else 6

    # ---------------- D1, D2
    bad1 = bad2 = 0
    for (q, x) in [(Fr(3, 7), Fr(5, 11)), (Fr(-2, 9), Fr(13, 4))]:
        O = VO.Ops(q, x, RMAX=14, SMAX=9)

        def h(F):
            a = O.add((1, F), (-1, O.sh(F)))
            return O.add((1, O.mul(a, lambda R, S: q ** (-R))), (1, O.mul(F, lambda R, S: x * q ** (2 * S) * q ** R)))

        def Oop(F):
            return O.mul(O.up(F), lambda R, S: q ** (R + 2 * S))

        def omega(F):
            return sum((O.w(r, s) * v * O.iQ[s] for (r, s), v in F.items()), Fr(0))

        for r0 in range(4):
            for s0 in range(4):
                d = {(r0, s0): Fr(1)}
                lhs = O.add((1, O.T(h(d))), (-1, h(O.T(d))))
                rhs = O.mul(Oop(O.T(d)), lambda R, S: -x)
                bad1 += VO.window_diff(lhs, rhs, 10, 5)
                if omega(h(d)) != omega(O.mul(d, lambda R, S: q ** (-R))):
                    bad2 += 1
    report('D1', bad1 == 0, '[T,h] = -x O T on point masses (32 cases, window R<=10, S<=5)')
    report('D2', bad2 == 0, 'omega(h F) = omega(P^-1 F) on point masses (32 cases)')

    # ---------------- D3
    ch = Chain(QM)
    xinv_ok = True
    allok = True
    neg_detected = True
    for k in range(1, KM + 1):
        K = ch.total(k)
        Kb = ch.total(k, ('rm1',))
        KO = BS()
        for i in range(1, k + 1):
            KO = KO + ch.total(k, ('diag', i))
            if i < k:
                KO = KO + BS({(0, 0): -1}) * ch.total(k, ('cross', i))
        diff = Kb + BS({(0, 0): -1}) * K
        # x^{-1} diff
        if any(i == 0 for (i, j) in diff.d):
            xinv_ok = False
        sh = BS({(i - 1, j): v for (i, j), v in diff.d.items()})
        # compare only x-powers that are complete: all terms with q^j < QM are exact, so compare everything
        if sh.d != KO.d:
            allok = False
            print('   D3 mismatch at k =', k)
        if k == 3:
            KObad = KO + BS({(0, 0): -1}) * ch.total(k, ('diag', k))
            if KObad.d == sh.d:
                neg_detected = False
        main.cache = getattr(main, 'cache', {})
        main.cache[k] = KO
    report('D3', allok and xinv_ok, 'K^Om_k = x^-1 (K^b_k - K_k), k <= %d, mod q^%d, all powers of x' % (KM, QM))
    report('D3neg', neg_detected, 'negative control (Omega without its last term) detected')

    # ---------------- D4
    ok = True
    for N in range(0, 6):
        for M in range(0, N + 1):
            c = phi_coeffs(N, M)
            L = N - M
            # T h phi = h T phi - x O T phi = lam (h phi - x O phi); slice 0: lam (h_x F0 - x q^R F1)
            # v_x(h_x F0) = sum_e (q^-e + x q^e) c0_e ; v_x(q^R F1) = 0 ; so omega(h phi) (xq)_inf = lam sum_e (q^-e + x q^e) c0_e
            # check: sum_e (q^-e + x q^e) c0_e == eps/(1+x) * sum_e c0_e  (=: eps W/(1+x))
            c0 = c[0]
            mu1 = sum(((mono(0, -e) + mono(1, e)) * v for e, v in c0.items()), RF(0))
            mu0 = sum((v for e, v in c0.items()), RF(0))
            if not (mu1 * RF(ONE + X) == eps(N, M) * mu0):
                ok = False
                print('   D4 mismatch', N, M)
    report('D4', ok, 'v_x(h_x phi(.,0)) = eps/(1+x) v_x(phi(.,0)) for N <= 5 (eigenvectors from the triangular system)')
    # D4b: the composite statement omega(h phi) = eps/(1+x) omega(phi), directly: h is diagonal on Xi_{e,s} with
    # eigenvalue q^-e + x q^{2s+e}, and (xq)_inf omega(Xi_{e,s}) = w(e,s) (xq)_{2s} / (Q;Q)_s by Cauchy (B-ii),
    # so both sides are finite sums of rational functions (no use of v_x or of the commutator).
    ok = True
    for N in range(0, 6):
        for M in range(0, N + 1):
            c = phi_coeffs(N, M)
            om = RF(0)
            omh = RF(0)
            for s, row in c.items():
                for e, v in row.items():
                    base = v * mono(e + 2 * s, (e + s) ** 2 + s * s) * poch(X * Q, 2 * s) / qpoch(s, 2)
                    om = om + base
                    omh = omh + base * (mono(0, -e) + mono(1, 2 * s + e))
            if not (omh * RF(ONE + X) == eps(N, M) * om):
                ok = False
                print('   D4b mismatch', N, M)
    report('D4b', ok, 'omega(h phi_{N,M}) = eps/(1+x) omega(phi_{N,M}) for N <= 5 by finite Cauchy sums')

    # ---------------- D5
    ok = True
    for N in range(0, 13):
        for M in range(0, N + 1):
            lhs = target(N, M) * (eps(N, M) - RF(ONE + X)) / (RF(X) * RF(ONE + X))
            rhs = lamn(N) * lamn(M) * (bO(N) * betaJ(M) + bO(M) * betaJ(N))
            if N == M:
                rhs = rhs / 2
            if not (lhs == rhs):
                ok = False
                print('   D5 mismatch', N, M)
    for n in range(0, 13):
        if not (lamn(n) ** 2 * bO(n) == dstar(n)):
            ok = False
            print('   D5 d* mismatch', n)
    report('D5', ok, 'spectral weight identity (N <= 12) and lam_n^2 bO_n = alpha*_n (n <= 12), symbolic')

    # ---------------- D6, D7
    BS.QM = QM
    P1 = xq_inf(QM)
    P2 = x2Q_inf(QM)
    ok6 = True
    ok7 = True
    neg6 = True
    for k in range(1, KM + 1):
        spec = BS()
        specbad = BS()
        for N in range(0, QM):
            if (k + 2) * N * N - 3 * N >= QM + 2:
                break
            for M in range(0, N + 1):
                if (N * N + M * M) * (k + 2) - 3 * N - 3 * M >= QM + 2:
                    break
                wt = target(N, M) * (eps(N, M) - RF(ONE + X)) / (RF(X) * RF(ONE + X))
                spec = spec + bs_from_rf(lam(N, M) ** k * wt)
                if k == 2:
                    wtb = target(N, M) * (eps(N, M) + RF(X) * RF(ONE + X) * mono(0, 1) - RF(ONE + X)) / (RF(X) * RF(ONE + X))
                    specbad = specbad + bs_from_rf(lam(N, M) ** k * wtb)
        KO = main.cache[k]
        lhs = KO * P1 * P2
        if lhs.d != spec.d:
            ok6 = False
            print('   D6 mismatch at k =', k)
        if k == 2 and lhs.d == specbad.d:
            neg6 = False
        # D7
        Bk = direct_Bk(k, QM)
        lhs7 = Bk * P1
        rhs7 = BS()
        for n in range(0, QM):
            if (k + 2) * n * n - 2 * n >= QM + 2:
                break
            rhs7 = rhs7 + bs_from_rf(lamn(n) ** (k + 1) * bO(n))
        if lhs7.d != rhs7.d:
            ok7 = False
            print('   D7 mismatch at k =', k)
    report('D6', ok6, '(xq)(x^2Q;Q) K^Om_k = sum lam^k t (eps-1-x)/(x(1+x)), k <= %d, mod q^%d, all x' % (KM, QM))
    report('D6neg', neg6, 'negative control (eps -> eps + x(1+x)q) detected')
    report('D7', ok7, '(xq)_inf B^(k)_{k+1}(x) = sum_n lam_n^{k+1} bO_n, k <= %d, mod q^%d, all x' % (KM, QM))

    # ---------------- D8 (x = 1, integer series in q)
    QQ = 200
    ok = True
    for k in range(1, 11):
        lhs = [0] * QQ
        n = 0
        while (k + 2) * n * n - 2 * n < QQ:
            # x = 1: lam_n^{k+1} bO_n = (-1)^n (q^{(k+2)n^2-2n} + q^{(k+2)n^2+2n}) for n >= 1, 1 for n = 0
            # (computed here from the rational function, not from this formula)
            f = lamn(n) ** (k + 1) * bO(n)
            # evaluate at x = 1 as a Laurent polynomial / series in q: substitute x = 1 via subs
            f1 = f.subs_x(ONE)
            ser = bs_from_rf_q(f1, QQ)
            for j, c in ser.items():
                if j < QQ:
                    lhs[j] += c
            n += 1
        # (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_inf
        rhs = [0] * QQ
        rhs[0] = 1
        m = 2 * k + 4
        for a in (k, k + 4, m):
            e = a
            while e < QQ:
                new = rhs[:]
                for j in range(QQ - e):
                    new[j + e] -= rhs[j]
                rhs = new
                e += m
        if lhs != rhs:
            ok = False
            print('   D8 mismatch at k =', k)
    report('D8', ok, 'sum_n lam_n^{k+1} bO_n at x = 1 = (q^k,q^{k+4},q^{2k+4};q^{2k+4})_inf mod q^%d, k <= 10' % QQ)

    print('ALL CHECKS PASSED' if FAILS == 0 else 'FAILURES: %d' % FAILS)


def bs_from_rf_q(f, QQ):
    """univariate (x already substituted) rational function in q -> integer q-series dict mod q^QQ (denominator +-1 at q=0)"""
    num = {}
    for (i, j), c in f.n.to_dict().items():
        num[int(j)] = num.get(int(j), 0) + int(c)
    den = {}
    for (i, j), c in f.d.to_dict().items():
        den[int(j)] = den.get(int(j), 0) + int(c)
    lo = min(den)
    # normalise: divide num and den by q^lo
    den = {j - lo: c for j, c in den.items() if c}
    num = {j - lo: c for j, c in num.items() if c}
    assert abs(den[0]) == 1
    s0 = den[0]
    lo_n = min(num) if num else 0
    assert lo_n >= 0, 'negative q power'
    F = [0] * QQ
    for j in range(QQ):
        acc = num.get(j, 0)
        for jj in range(1, j + 1):
            if jj in den:
                acc -= den[jj] * F[j - jj]
        F[j] = acc * s0
    return {j: c for j, c in enumerate(F) if c}


if __name__ == '__main__':
    main()
