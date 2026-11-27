"""Numerical test of the x-refined FACTORISATION conjecture (F_k) for the KRTW identities (ii), (iii)
(proofs/KRTW-general.md, section 8).  Not a proof.

    K_k(x) := sum_{r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k} x^{sum (r_i+2s_i)} prod_i q^{(r_i+s_i)^2+s_i^2}
                                                     / ((q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}})
    (F_k):   K_k(x) = B_{k+1}(x) * J_k(x),
    B_{k+1}(x) = sum_{N_1>=...>=N_{k+1}>=0} x^{sum N_i} q^{sum N_i^2} / ((q)_{n_1}...(q)_{n_k} (q^2;q^2)_{n_{k+1}})   (Bressoud)
    J_k(x)     = (1/(x^2q^2;q^2)_inf) sum_{n>=0} (-1)^n x^{(k+1)n} q^{(k+1)n^2+n(n-1)} (x^2;q^2)_n (1-x^2q^{4n})
                                                                                      / ((q^2;q^2)_n (1-x^2)).
x = 1 gives (iii), x = q gives (ii) (Bressoud's even-modulus identities for B, Jacobi triple product for J).
Also tested: the same statement with Omega inserted (identity (iv)) for a few natural candidates.

Exact arithmetic modulo q^N, all powers of x.   usage: python3 krtw_Fk_test.py [N] [KMAX]
"""
import sys
from flint import fmpz_poly
import krtw_ser as K
from krtw_ser import S, tr, mono, iqpoch, xpoch, xpoch_inv

def transfer_x(k, N, omega_term=None):
    """K_k(x) (with the i-th term of Omega inserted if omega_term = i), summed level by level from the definition"""
    rmax = 2
    while rmax * rmax < N:
        rmax += 1
    Phi = {(r, 0): S.one() for r in range(rmax + 1)}
    for i in range(1, k + 1):
        new = {}
        for s in range(rmax + 1):
            g = {}
            for r in range(rmax + 1):
                e = (r + s) ** 2 + s * s + (r + 2 * s if omega_term == i else 0)
                if e >= N:
                    break
                acc = S()
                for sp in range(s + 1):
                    p = Phi.get((r, sp))
                    if p is None:
                        continue
                    b = s - sp
                    if omega_term is not None and omega_term == i - 1:
                        if b == 0:
                            continue
                        acc = acc + p * iqpoch(b - 1, 2)
                    else:
                        acc = acc + p * iqpoch(b, 2)
                if not acc.iszero():
                    g[r] = acc.shift(r + 2 * s, e)
            if not g:
                continue
            for R in range(rmax + 1):
                acc = S()
                for r, v in g.items():
                    if r >= R:
                        acc = acc + v * iqpoch(r - R)
                if not acc.iszero():
                    new[(R, s)] = acc
        Phi = new
    tot = S()
    for (R, s), v in Phi.items():
        if R == 0:
            tot = tot + v
    return tot

def bress(nvar, N, lin_from=None):
    """sum x^{sum N_i} q^{sum N_i^2 + N_{lin_from}+...+N_nvar} / ((q)_{n_1}...(q)_{n_{nvar-1}} (q^2;q^2)_{n_nvar})"""
    s = S()
    def rec(Ns, e, c):
        nonlocal s
        lev = len(Ns)
        if lev == nvar:
            s = s + S.xq(sum(Ns), e) * tr(c * iqpoch(Ns[-1], 2))
            return
        M = 0
        while lev == 0 or M <= Ns[-1]:
            ee = e + M * M + (M if (lin_from is not None and lev + 1 >= lin_from) else 0)
            if ee >= N:
                break
            cc = c if lev == 0 else tr(c * iqpoch(Ns[-1] - M))
            rec(Ns + [M], ee, cc)
            M += 1
    rec([], 0, fmpz_poly([1]))
    return s

def J(k, N):
    s = S.one()
    n = 1
    while True:
        e = (k + 1) * n * n + n * (n - 1)
        if e >= N:
            break
        a = xpoch(-1, 2, 2, n - 1, 2) * (S.one() - S.xq(2, 4 * n))      # (x^2;q^2)_n (1-x^2q^{4n})/(1-x^2)
        s = s + (a * iqpoch(n, 2)).shift((k + 1) * n, e) * ((-1) ** n)
        n += 1
    return s * xpoch_inv(-1, 2, 2, N, 2)

def Cchain(k, N):
    """C_k(x) = sum_{0<=s_1<=...<=s_k} x^{sum s_i} q^{sum s_i^2} / prod (q^2;q^2)_{s_i-s_{i-1}}"""
    s = S()
    def rec(ss, e, c):
        nonlocal s
        if len(ss) == k:
            s = s + S.xq(sum(ss), e) * c
            return
        v = ss[-1] if ss else 0
        while True:
            ee = e + v * v
            if ee >= N:
                break
            rec(ss + [v], ee, tr(c * iqpoch(v - (ss[-1] if ss else 0), 2)))
            v += 1
    rec([], 0, fmpz_poly([1]))
    return s

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    K.setN(N)
    print("x-refined factorisation (F_k): K_k(x) = B_{k+1}(x) J_k(x), all powers of x, mod q^%d" % N)
    allok = True
    for k in range(1, KMAX + 1):
        Jk = J(k, N)
        ok = (transfer_x(k, N) == bress(k + 1, N) * Jk)
        allok &= ok
        # Omega-weighted sum (identity (iv)) against candidates  B^{(i)}_{k+1}(x) J_k(x)
        Om = S()
        for i in range(1, k + 1):
            Om = Om + transfer_x(k, N, omega_term=i)
        cands = []
        for i in range(1, k + 2):
            if Om == bress(k + 1, N, lin_from=i) * Jk:
                cands.append(i)
        print("k=%d: (F_k) %s ;  Omega-sum(x) == B_{k+1}^{(linear terms from level i)}(x) J_k(x) for i in %s"
              % (k, "OK" if ok else "FAIL", cands), flush=True)
    print("ALL (F_k) OK" if allok else "SOME (F_k) FAIL")
    # proved statements used in section 8:  J_k(x) C_1(x) = C_k(x),  C_1(x) = (-xq;q^2)_inf;
    # and the specialisations x = 1, q of B_{k+1}(x) J_k(x) (resp. with Omega) are the products (iii), (ii), (iv)
    from krtw_ser import xpoch_inf
    from krtw_four_test import products, first_diff
    ok2 = (Cchain(1, N) == xpoch_inf(1, 1, 1, 2))
    ok3 = True
    for k in range(1, KMAX + 1):
        Jk = J(k, N)
        ok2 &= (Jk * Cchain(1, N) == Cchain(k, N))
        P = products(k, N)
        F = bress(k + 1, N) * Jk
        ok3 &= first_diff(F.at1(), P["iii"], N) is None
        ok3 &= first_diff(F.subx(1).at1(), P["ii"], N) is None
        ok3 &= first_diff((bress(k + 1, N, lin_from=k) * Jk).at1(), P["iv"], N) is None
    print("J_k(x) C_1(x) = C_k(x) and C_1(x) = (-xq;q^2)_inf, k <= %d: %s" % (KMAX, "OK" if ok2 else "FAIL"))
    print("B_{k+1}(x)J_k(x) at x=1, x=q and B^{(k)}_{k+1}(1)J_k(1) equal the products (iii), (ii), (iv), k <= %d: %s" % (KMAX, "OK" if ok3 else "FAIL"))
