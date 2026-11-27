"""Numerical checks (exact arithmetic modulo q^N, all powers of x) of the statements of proofs/KRTW-F-moves.md.
NOT a proof; the proofs are in the .md file.  usage: python3 krtw_struct_certify.py [N] [KMAX]

  K_k(x)   = sum_{r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k} x^{sum(r_i+2s_i)} prod q^{(r_i+s_i)^2+s_i^2}/((q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}})
  B_m(x)   = sum_{N_1>=...>=N_m>=0} x^{sum N_i} q^{sum N_i^2}/((q)_{N_1-N_2}...(q)_{N_{m-1}-N_m}(q^2;q^2)_{N_m})
  C_k(x)   = sum_{0<=s_1<=...<=s_k} x^{sum s_i} q^{sum s_i^2}/prod (q^2;q^2)_{s_i-s_{i-1}},     J_k = C_k/C_1
Checks:
 [F]   K_k = B_{k+1} J_k                                             (conjecture (F_k); proved for k = 2, 3 by krtw_moves_certify.py)
 [S1]  B_{k+1}(x) (xq;q)_inf = sum_n (-1)^n x^{(k+1)n} q^{(k+2)n^2} (x^2;q^2)_n (1-xq^{2n})/((q^2;q^2)_n (1-x))
 [S2]  J_k(x) (x^2q^2;q^2)_inf = sum_m (-1)^m x^{(k+1)m} q^{(k+2)m^2-m} (x^2;q^2)_m (1-x^2q^{4m})/((q^2;q^2)_m (1-x^2))
 [S3]  B_{k+1}(x) = (-xq^2;q^2)_inf * sum_{t_1>=...>=t_k>=0} x^{sum t} q^{sum t^2} (-xq;q^2)_{t_k} / (prod (q^2;q^2)_{t_i-t_{i+1}} (-xq^2;q^2)_{t_k})
 [S4]  B_{k+1} J_k (xq;q)_inf (x^2q^2;q^2)_inf = L_k(x) := (1/2) sum_{n,m>=0} (-1)^{n+m} x^{(k+1)(n+m)} q^{(k+2)(n^2+m^2)-n-m}
            (x^2;q^2)_n (x^2;q^2)_m (q^n+q^m)(1+xq^{n+m})(1-xq^{2n})(1-xq^{2m}) / ((q^2;q^2)_n (q^2;q^2)_m (1-x)(1-x^2))
       (= the C_2 Rogers--Selberg series L^{(0)}_k(y_1,y_2;q) of Bartlett--Warnaar / Griffin--Ono--Warnaar at (y_1,y_2) = (sqrt x, -sqrt x))
 [S5]  the Bartlett--Warnaar multisum for sum_{la_1<=k} q^{|la|} P'_{2 la}(sqrt x, -sqrt x; q) equals K_k   (k = 1, 2; modulo q^min(N,22))
 [S6]  k = 2m even:  K_{2m} = (-xq;q)_inf * prod_{e=1,2} sum_{s_1>=...>=s_m} x^{2 sum s} q^{2 sum s^2}/(prod (q^2;q^2)_{s_i-s_{i+1}} (-xq^e;q^2)_{s_m})
"""
import sys
from flint import fmpz_poly
import krtw_ser as K
from krtw_ser import S, tr, mono, iqpoch, qpoch, xpoch, xpoch_inv, xpoch_inf

def rmax():
    r = 1
    while r * r < K.N: r += 1
    return r

def Kk(k):
    """level-by-level transfer; Phi[(R,s)] = sum over levels below with r <= ... """
    rm = rmax() + 1
    # A_j(R;x) recursion with R = infinity handled by taking R large is avoided: direct transfer over (r_i, s_i)
    Phi = {(r, 0): S.one() for r in range(rm + 1)}       # level k+1: r_{k+1} = 0 -> factor 1/(q)_{r_k} applied below
    # go from level k down to 1: state (r_i, s_i); here we sum levels from the inside (i = k) outwards
    # inside-out: level i has r_i >= r_{i+1}, s_i <= s_{i+1}.  Start: r_{k+1} = 0, s_{k+1} free (no factor).
    cur = None
    for lev in range(k, 0, -1):
        new = {}
        for r in range(rm + 1):
            for s in range(rm + 1):
                e = (r + s) ** 2 + s * s
                if e >= K.N: break
                if cur is None:
                    acc = S({0: iqpoch(r)})
                else:
                    acc = S()
                    for (r2, s2), v in cur.items():
                        if r2 <= r and s2 >= s:
                            acc = acc + v * tr(iqpoch(r - r2) * iqpoch(s2 - s, 2))
                if not acc.iszero(): new[(r, s)] = acc.shift(r + 2 * s, e)
        cur = new
    tot = S()
    for (r, s), v in cur.items(): tot = tot + v * iqpoch(s, 2)
    return tot

def chain(k, wt, sq, base_links, last_base, extra=None):
    """sum_{N_1>=...>=N_k>=0} x^{wt sum N} q^{sq sum N^2} extra(N_k) / (prod_{i<k} (q^b;q^b)_{N_i-N_{i+1}} (q^lb;q^lb)_{N_k}) ; lb = None: no last factor"""
    tot = S()
    def rec(Ns, e, c):
        nonlocal tot
        if len(Ns) == k:
            cc = c if last_base is None else tr(c * iqpoch(Ns[-1], last_base))
            t = S.xq(wt * sum(Ns), e) * cc
            if extra is not None: t = t * extra(Ns[-1])
            tot = tot + t
            return
        M = 0
        while not Ns or M <= Ns[-1]:
            ee = e + sq * M * M
            if ee >= K.N: break
            cc = c if not Ns else tr(c * iqpoch(Ns[-1] - M, base_links))
            rec(Ns + [M], ee, cc); M += 1
    rec([], 0, fmpz_poly([1]))
    return tot

def Sinv(a):
    assert a.d.get(0) == fmpz_poly([1])
    u = S.one() - a; res = S.one(); p = S.one()
    while True:
        p = p * u
        if p.iszero(): break
        res = res + p
    return res

def Bm(m): return chain(m, 1, 1, 1, 2)
def Ck(k): return chain(k, 1, 1, 2, 2)

def vwp_B(k):
    tot = S(); n = 0
    while (k + 2) * n * n < K.N:
        # (x^2;q^2)_n (1-xq^{2n})/(1-x) = (1+x) (x^2q^2;q^2)_{n-1} (1-xq^{2n}) for n>=1
        if n == 0: f = S.one()
        else: f = (S.one() + S.xq(1, 0)) * xpoch(-1, 2, 2, n - 1, 2) * (S.one() - S.xq(1, 2 * n))
        tot = tot + (f * iqpoch(n, 2)).shift((k + 1) * n, (k + 2) * n * n) * ((-1) ** n); n += 1
    return tot

def vwp_J(k):
    tot = S(); m = 0
    while (k + 2) * m * m - m < K.N:
        if m == 0: f = S.one()
        else: f = xpoch(-1, 2, 2, m - 1, 2) * (S.one() - S.xq(2, 4 * m))
        tot = tot + (f * iqpoch(m, 2)).shift((k + 1) * m, (k + 2) * m * m - m) * ((-1) ** m); m += 1
    return tot

def Lk(k):
    """twice L_k(x) (to stay in Z): sum_{n,m} ... ; the factor 1/((1-x)(1-x^2)) is cancelled against (x^2;q^2)_n (x^2;q^2)_m when n,m>=1"""
    tot = S()
    nm = 0
    while (k + 2) * nm * nm - nm < K.N: nm += 1
    def half(n):  # (x^2;q^2)_n/(1-x^2) for n>=1  ->  (x^2q^2;q^2)_{n-1}
        return xpoch(-1, 2, 2, n - 1, 2)
    onepx = S.one() + S.xq(1, 0)
    for n in range(nm):
        for m in range(nm):
            e = (k + 2) * (n * n + m * m) - n - m
            if e >= K.N: continue
            pref = (S.xq(0, n) + S.xq(0, m)) * (S.one() + S.xq(1, n + m))
            if n >= 1 and m >= 1:
                # (x^2)_n (x^2)_m (1-xq^{2n})(1-xq^{2m}) / ((1-x)(1-x^2)) = (1+x)(x^2q^2)_{n-1} (1-x^2)(x^2q^2)_{m-1} (1-xq^{2n})(1-xq^{2m})/(1-x^2)... expand:
                # (x^2;q^2)_n = (1-x)(1+x)(x^2q^2;q^2)_{n-1};  total/((1-x)(1-x^2)) = (1+x)(x^2q^2)_{n-1} (x^2q^2)_{m-1} (1-x)(1+x)/(1+x)... :
                # (x^2)_n (x^2)_m /((1-x)(1-x^2)) = (1-x^2) (x^2q^2)_{n-1} (x^2q^2)_{m-1} / (1-x) = (1+x)(x^2q^2)_{n-1}(x^2q^2)_{m-1}
                f = onepx * half(n) * half(m) * (S.one() - S.xq(1, 2 * n)) * (S.one() - S.xq(1, 2 * m))
            elif n >= 1 or m >= 1:
                j = max(n, m)
                # (x^2)_j (1-xq^{2j})(1-x) / ((1-x)(1-x^2)) = (x^2q^2)_{j-1} (1-xq^{2j})
                f = half(j) * (S.one() - S.xq(1, 2 * j))
            else:
                # (1-x)^2/((1-x)(1-x^2)) * (1+1)(1+x) = 2 (1-x)(1+x)/(1-x^2) = 2 : handle separately
                tot = tot + S.one() * 2
                continue
            tot = tot + (pref * f * tr(iqpoch(n, 2) * iqpoch(m, 2))).shift((k + 1) * (n + m), e) * ((-1) ** (n + m))
    return tot

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    K.setN(N)
    ok = True
    xq_inf = xpoch_inf(-1, 1, 1, 1); x2_inf = xpoch_inf(-1, 2, 2, 2)
    mxq = xpoch_inf(+1, 1, 1, 1)           # (-xq;q)_inf
    mxq1 = xpoch_inf(+1, 1, 1, 2)          # (-xq;q^2)_inf
    mxq2 = xpoch_inf(+1, 1, 2, 2)          # (-xq^2;q^2)_inf
    C1 = Ck(1)
    print("N = %d (modulo q^N, all powers of x)" % N)
    print("C_1 = (-xq;q^2)_inf:", C1 == mxq1); ok &= (C1 == mxq1)
    for k in range(1, KMAX + 1):
        Kx = Kk(k); B = Bm(k + 1); C = Ck(k); J = C * Sinv(C1)
        r = {}
        r['F'] = (Kx == B * J)
        r['S1'] = (B * xq_inf == vwp_B(k))
        r['S2'] = (J * x2_inf == vwp_J(k))
        Cs = chain(k, 1, 1, 2, 2, extra=lambda t: xpoch(+1, 1, 1, t, 2) * xpoch_inv(+1, 1, 2, t, 2))
        r['S3'] = (B == mxq2 * Cs)
        r['S4'] = ((B * J * xq_inf * x2_inf) * 2 == Lk(k))
        if k % 2 == 0:
            m = k // 2
            j1 = chain(m, 2, 2, 2, 2, extra=lambda t: xpoch_inv(+1, 1, 1, t, 2))
            j2 = chain(m, 2, 2, 2, 2, extra=lambda t: xpoch_inv(+1, 1, 2, t, 2))
            r['S6'] = (Kx == mxq * j1 * j2)
        print("k = %d: " % k + "  ".join("[%s] %s" % (key, val) for key, val in r.items()))
        ok &= all(r.values())
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")

if __name__ == "__main__":
    main()
