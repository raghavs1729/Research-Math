"""Referee 2, computation (1): KRTW identity (iv) as a q-series, literally.

LHS(k) = sum_{r_1..r_k, s_1..s_k >= 0} Omega(r,s) prod_i q^{(r_i+s_i)^2+s_i^2} / ((q;q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}}),
r_{k+1} = s_0 = 0, 1/(q)_m = 0 for m < 0,
Omega = sum_{i=1}^{k-1} q^{r_i+2s_i} (1 - q^{2s_{i+1}-2s_i}) + q^{r_k+2s_k}       (KRTW 2404.03851v2, after Conj 4.9).
RHS(k) = (q^k,q^{k+1},q^{k+3},q^{k+4},q^{2k+4},q^{2k+4};q^{2k+4})_inf / ((q^2;q^2)_inf (q;q)_inf).

Method A (all k): level-by-level evaluation of the multisum.  G_i(R,S) := sum over levels i+1..k (weights and bonds),
carried with a flag "Omega term already inserted".  The Omega term i (< k) is inserted literally at the bond (i,i+1) as
the factor q^{r_i+2s_i}(1-q^{2(s_{i+1}-s_i)}); the term k as q^{r_k+2s_k} at level k.  No identity is used.
Method B (k <= 3): naive enumeration of all tuples.
Check: LHS * (q;q)_inf (q^2;q^2)_inf == numerator product, mod q^L.
Usage: python3 r2_iv_series.py [Nmax=300] [kmax=12]
"""
import sys, time, itertools
from flint import fmpz_poly
from r2_lib import mono, inv, poch, InvCache

Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 300
kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
L = Nmax + 1


def states_for(L):
    st = []
    R = 0
    while R * R < L:
        S = 0
        while (R + S) ** 2 + S * S < L:
            st.append((R, S))
            S += 1
        R += 1
    return st


def lhs_dp(k, L, cache, variant='iv'):
    """variant: 'iv' (KRTW Omega), 'iii' (no Omega), 'ii' (weight q^{r_i+2s_i} at every level),
    negative controls: 'nc_last' (last term q^{r_k+s_k}), 'nc_nofactor' (no (1-q^{..}) factors),
    'nc_half' ((1-q^{s_{i+1}-s_i}) instead of (1-q^{2s_{i+1}-2s_i}))."""
    st = states_for(L)
    wexp = {(r, s): (r + s) ** 2 + s * s for (r, s) in st}
    one = fmpz_poly([1])
    # kernels
    kern = {}
    kern_ins = {}

    def K(dr, ds):
        if (dr, ds) not in kern:
            kern[(dr, ds)] = cache.q(dr).mul_low(cache.Q(ds), L)
        return kern[(dr, ds)]

    def Kins(dr, ds):
        # literal insertion factor (1 - q^{2 ds}) (or (1-q^{ds}) for nc_half, 1 for nc_nofactor)
        if (dr, ds) not in kern_ins:
            if variant == 'nc_nofactor':
                f = one
            elif variant == 'nc_half':
                f = one - mono(ds, L)
            else:
                f = one - mono(2 * ds, L)
            kern_ins[(dr, ds)] = K(dr, ds).mul_low(f, L)
        return kern_ins[(dr, ds)]

    def level_factor(R, S):
        if variant == 'ii':
            return R + 2 * S
        return 0

    # level k
    G0 = {}
    G1 = {}
    for (R, S) in st:
        G0[(R, S)] = cache.q(R).mul_low(mono(level_factor(R, S), L), L)
        if variant == 'nc_last':
            G1[(R, S)] = cache.q(R).mul_low(mono(R + S, L), L)
        else:
            G1[(R, S)] = cache.q(R).mul_low(mono(R + 2 * S, L), L)
    for i in range(k - 1, 0, -1):
        N0, N1 = {}, {}
        for (R, S) in st:
            a0 = fmpz_poly()
            a1 = fmpz_poly()
            ains = fmpz_poly()
            for (r, s) in st:
                if r > R or s < S:
                    continue
                e = wexp[(r, s)]
                g0 = G0[(r, s)]
                g1 = G1[(r, s)]
                if not g0.is_zero():
                    t = K(R - r, s - S).mul_low(g0, L - e).left_shift(e)
                    a0 += t
                    if variant in ('iv', 'nc_last', 'nc_nofactor', 'nc_half'):
                        ains += Kins(R - r, s - S).mul_low(g0, L - e).left_shift(e)
                if not g1.is_zero():
                    a1 += K(R - r, s - S).mul_low(g1, L - e).left_shift(e)
            lf = level_factor(R, S)
            N0[(R, S)] = a0.mul_low(mono(lf, L), L) if lf else a0.truncate(L)
            N1[(R, S)] = (a1 + ains.mul_low(mono(R + 2 * S, L), L)).truncate(L)
        G0, G1 = N0, N1
    tot0 = fmpz_poly()
    tot1 = fmpz_poly()
    for (R, S) in st:
        e = wexp[(R, S)]
        tot0 += cache.Q(S).mul_low(G0[(R, S)], L - e).left_shift(e)
        tot1 += cache.Q(S).mul_low(G1[(R, S)], L - e).left_shift(e)
    tot0 = tot0.truncate(L)
    tot1 = tot1.truncate(L)
    if variant in ('iii', 'ii'):
        return tot0
    return tot1


def lhs_naive(k, L, cache):
    """naive enumeration over all tuples in a box (Pochhammer with negative index -> factor 0)."""
    Rmax = 0
    while Rmax * Rmax < L:
        Rmax += 1
    tot = fmpz_poly()
    rng = range(Rmax + 1)
    for rs in itertools.product(rng, repeat=k):
        if any(rs[i] < rs[i + 1] for i in range(k - 1)):
            continue  # 1/(q)_{neg} = 0
        base_r = sum(r * r for r in rs)
        if base_r >= L:
            continue
        for ss in itertools.product(rng, repeat=k):
            if any(ss[i] > ss[i + 1] for i in range(k - 1)):
                continue  # 1/(Q;Q)_{neg} = 0
            e = sum((rs[i] + ss[i]) ** 2 + ss[i] ** 2 for i in range(k))
            if e >= L:
                continue
            term = mono(e, L)
            for i in range(k):
                rn = rs[i + 1] if i + 1 < k else 0
                sp = ss[i - 1] if i > 0 else 0
                term = term.mul_low(cache.q(rs[i] - rn), L).mul_low(cache.Q(ss[i] - sp), L)
            # Omega, literally
            om = fmpz_poly()
            for i in range(k - 1):
                om += mono(rs[i] + 2 * ss[i], L) * (fmpz_poly([1]) - mono(2 * ss[i + 1] - 2 * ss[i], L))
            om += mono(rs[k - 1] + 2 * ss[k - 1], L)
            tot += term.mul_low(om, L)
    return tot.truncate(L)


def rhs_num(k, L, which='iv'):
    m = 2 * k + 4
    if which == 'iv':
        exps = [k, k + 1, k + 3, k + 4, m, m]
    elif which == 'iii':
        exps = [k + 1, k + 2, k + 2, k + 3, m, m]
    elif which == 'ii':
        exps = [1, 2, 2 * k + 2, 2 * k + 3, m, m]
    r = fmpz_poly([1])
    for a in exps:
        r = r.mul_low(poch(a, m, 'inf', L), L)
    return r


def main():
    t0 = time.time()
    out = []
    pr = lambda s: (print(s, flush=True), out.append(s))
    pr("Referee 2, (1): KRTW (iv) as q-series, literal multisum. Nmax = %d (mod q^%d), k = 1..%d" % (Nmax, L, kmax))
    nm = 4 * int(L ** 0.5) + 10
    cache = InvCache(L, nm)
    den = poch(1, 1, 'inf', L).mul_low(poch(2, 2, 'inf', L), L)  # (q;q)_inf (q^2;q^2)_inf
    allok = True
    # validate DP against naive enumeration, small L
    Ls = 61
    cs = InvCache(Ls, 4 * int(Ls ** 0.5) + 10)
    for k in (1, 2, 3):
        a = lhs_dp(k, Ls, cs)
        b = lhs_naive(k, Ls, cs)
        ok = (a - b).truncate(Ls).is_zero()
        allok &= ok
        pr("  DP vs naive enumeration, k=%d, mod q^%d: %s" % (k, Ls, "OK" if ok else "MISMATCH"))
    lhs = {}
    for k in range(1, kmax + 1):
        t1 = time.time()
        lhs[k] = lhs_dp(k, L, cache)
        lhs_den = lhs[k].mul_low(den, L)
        ok = (lhs_den - rhs_num(k, L)).truncate(L).is_zero()
        allok &= ok
        cs10 = [int(c) for c in lhs[k].coeffs()[:12]]
        pr("  (iv) k=%2d: LHS*(q)_inf(q^2;q^2)_inf == numerator mod q^%d: %s   [LHS = %s ...]  (%.1fs)"
           % (k, L, "OK" if ok else "FAIL", cs10, time.time() - t1))
    # positive controls (ii), (iii) for a few k
    for k in (1, 2, 3, 4):
        for which in ('ii', 'iii'):
            v = lhs_dp(k, L, cache, variant=which).mul_low(den, L)
            ok = (v - rhs_num(k, L, which)).truncate(L).is_zero()
            allok &= ok
            pr("  positive control (%s) k=%d mod q^%d: %s" % (which, k, L, "OK" if ok else "FAIL"))
    # negative controls
    def firstdiff(a, b):
        d = (a - b).truncate(L)
        cs = d.coeffs()
        for j, v in enumerate(cs):
            if v != 0:
                return j
        return None
    pr("  negative controls (each must FAIL; first differing power of q shown):")
    for k in (3, 5):
        j = firstdiff(lhs[k].mul_low(den, L), rhs_num(k + 1, L))
        pr("    NC1 LHS(iv,k=%d) vs RHS(iv,k=%d): %s" % (k, k + 1, "fails at q^%s (good)" % j if j is not None else "AGREES (bad)"))
        allok &= j is not None
        for v in ('nc_last', 'nc_nofactor', 'nc_half'):
            j = firstdiff(lhs_dp(k, L, cache, variant=v).mul_low(den, L), rhs_num(k, L))
            pr("    NC %-11s k=%d vs RHS(iv,k=%d): %s" % (v, k, k, "fails at q^%s (good)" % j if j is not None else "AGREES (bad)"))
            allok &= j is not None
    pr("RESULT: %s   (%.1fs)" % ("ALL OK" if allok else "SOME CHECK FAILED", time.time() - t0))
    return out


if __name__ == '__main__':
    main()
