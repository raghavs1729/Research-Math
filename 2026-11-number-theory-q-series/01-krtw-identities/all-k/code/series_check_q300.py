"""Independent series check of KRTW identities (ii), (iii), (iv) (arXiv:2404.03851v2, after Conj. 4.9) mod q^N.

Written from scratch for krtw4 (no code shared with ../rr-identities).  Sums are evaluated literally:
  sum over r_1..r_k, s_1..s_k >= 0 with 1/(q)_m = 0 for m < 0 (i.e. r decreasing, s increasing),
  r_{k+1} = s_0 = 0, by a level-by-level transfer from level k to level 1, in exact integer arithmetic (flint fmpz_poly).
(iv) uses Omega = sum_{i<k} q^{r_i+2s_i}(1 - q^{2s_{i+1}-2s_i}) + q^{r_k+2s_k}, expanded literally as the sum of 2k-1
signed monomial insertions (no rewriting of the bond factor).
usage: python3 series_check_q300.py [N=300] [KMAX=8]
"""
import sys
import flint
from flint import fmpz_poly

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8


def trunc(p):
    c = p.coeffs()[:N]
    return fmpz_poly(c)


def mul(a, b):
    return a.mul_low(b, N)


def mono(e):
    if e >= N:
        return fmpz_poly([])
    return fmpz_poly([0] * e + [1])


def inv_series(p):
    """1/p mod q^N for p with constant term 1 (Newton iteration with truncated products)"""
    assert p.coeffs() and p.coeffs()[0] == 1
    g = fmpz_poly([1])
    m = 1
    while m < N:
        m = min(2 * m, N)
        pg = p.mul_low(g, m)
        g = g.mul_low(fmpz_poly([2]) - pg, m)
    return g


def poch_inv(n, step):
    """1/(q^step; q^step)_n mod q^N"""
    p = fmpz_poly([1])
    for i in range(1, n + 1):
        p = mul(p, fmpz_poly([1] + [0] * (step * i - 1) + [-1]) if step * i < N else fmpz_poly([1]))
    return inv_series(p)


def prod_inf(residues, mod):
    """prod_{n >= 1, n mod 'mod' in residues (with multiplicity)} (1 - q^n) mod q^N; residues is a list"""
    p = fmpz_poly([1])
    for a in residues:
        n = a
        while n < N:
            p = mul(p, fmpz_poly([1] + [0] * (n - 1) + [-1]))
            n += mod
    return p


# caches
IQ = {}
IQQ = {}


def iq(n):
    if n not in IQ:
        IQ[n] = poch_inv(n, 1)
    return IQ[n]


def iQ(n):
    if n not in IQQ:
        IQQ[n] = poch_inv(n, 2)
    return IQQ[n]


def states():
    out = []
    for r in range(0, N):
        for s in range(0, N):
            if (r + s) ** 2 + s * s >= N:
                break
            out.append((r, s))
        if r * r >= N:
            break
    return out


ST = states()


def level_sum(k, lin, ins=None):
    """sum over chains of prod_i q^{(r_i+s_i)^2+s_i^2 + lin(r_i,s_i)} / bonds, with optional insertion
    ins = (i, a, b): extra factor q^{a r_i + b s_i}... general: ins is a list of (level i, function(r,s)->exponent, bond_flag)
    Here we only need: insertion of q^{r_i+2s_i} at level i, optionally times q^{2 s_{i+1}} (inner neighbour) with a sign."""
    raise NotImplementedError


def chain(k, lin, extra=None):
    """extra: None, or (i, mode) with mode 'diag' -> factor q^{r_i+2s_i},
       mode 'cross' -> factor q^{r_i + 2 s_{i+1}} (i < k).
    Transfer from level k (innermost, r_{k+1}=0) to level 1 (outermost, s_0=0)."""
    # F[(r,s)] for the current level = sum over inner levels (including this level's weight)
    F = {}
    for lev in range(k, 0, -1):
        newF = {}
        for (r, s) in ST:
            e = (r + s) ** 2 + s * s + lin(r, s)
            if extra is not None and extra[0] == lev:
                if extra[1] == 'diag':
                    e += r + 2 * s
                else:           # cross: q^{r_i} here, q^{2 s_{i+1}} attached to the inner level below
                    e += r
            if e >= N:
                continue
            if lev == k:
                acc = iq(r)     # 1/(q)_{r_k - r_{k+1}}, r_{k+1} = 0
            else:
                acc = fmpz_poly([])
                for (r2, s2), val in F.items():
                    if r2 > r or s2 < s:
                        continue
                    t = mul(val, iq(r - r2))
                    t = mul(t, iQ(s2 - s))
                    if extra is not None and extra[1] == 'cross' and extra[0] == lev:
                        t = mul(t, mono(2 * s2))
                    acc += t
                acc = trunc(acc)
            if acc.degree() >= 0:
                newF[(r, s)] = mul(acc, mono(e))
        F = newF
    tot = fmpz_poly([])
    for (r, s), val in F.items():
        tot += mul(val, iQ(s))       # 1/(q^2;q^2)_{s_1 - s_0}
    return trunc(tot)


def main():
    den = mul(prod_inf([2], 2), prod_inf([1], 1))         # (q^2;q^2)_inf (q;q)_inf
    deninv = inv_series(den)
    allok = True
    for k in range(1, KMAX + 1):
        m = 2 * k + 4
        # (ii)
        lhs = chain(k, lambda r, s: r + 2 * s)
        rhs = mul(prod_inf([1, 2, 2 * k + 2, 2 * k + 3, m, m], m), deninv)
        ok2 = (lhs == rhs)
        # (iii)
        lhs = chain(k, lambda r, s: 0)
        rhs = mul(prod_inf([k + 1, k + 2, k + 2, k + 3, m, m], m), deninv)
        ok3 = (lhs == rhs)
        # (iv): Omega = sum_{i<k} [q^{r_i+2s_i} - q^{r_i+2s_{i+1}}] + q^{r_k+2s_k}
        lhs = fmpz_poly([])
        for i in range(1, k + 1):
            lhs += chain(k, lambda r, s: 0, (i, 'diag'))
            if i < k:
                lhs -= chain(k, lambda r, s: 0, (i, 'cross'))
        lhs = trunc(lhs)
        rhs = mul(prod_inf([k, k + 1, k + 3, k + 4, m, m], m), deninv)
        ok4 = (lhs == rhs)
        print('k=%d  (ii): %s  (iii): %s  (iv): %s   [mod q^%d]' % (k, 'OK' if ok2 else 'FAIL', 'OK' if ok3 else 'FAIL',
                                                              'OK' if ok4 else 'FAIL', N))
        sys.stdout.flush()
        allok = allok and ok2 and ok3 and ok4
    # negative controls (k = 3): wrong linear term, wrong product
    k = 3; m = 10
    bad1 = chain(k, lambda r, s: r + s) == mul(prod_inf([1, 2, 8, 9, m, m], m), deninv)
    bad2 = chain(k, lambda r, s: 0) == mul(prod_inf([4, 5, 5, 7, m, m], m), deninv)
    print('negative controls detected:', (not bad1) and (not bad2))
    print('ALL OK' if allok and not bad1 and not bad2 else 'SOMETHING FAILED')


if __name__ == '__main__':
    main()
