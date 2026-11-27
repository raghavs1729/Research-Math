#!/usr/bin/env python3
"""
t2_lll.py -- zero-side witnesses by lattice reduction (simultaneous inhomogeneous approximation
on the heaviest ordinates), for the square classes that exhaustive scans cannot reach.

Under GRH(q),  F_{q,a}(y) = sum_j A_j cos(gamma_j y + theta_j) + (zeros not located),
one cosine per located zero (A_j = 2/|rho(rho+1)|, theta_j = arg(-conj(chi(a))/(rho(rho+1))); a zero
gamma of chi and the zero -gamma of chi-bar give the same cosine).  As in t2_zeros.py:
    F_{q,a}(y) >= F_Z(y) - r   for every real y,   r = L_1(q) - sum_j A_j   (closed form, Lemma S).

SEARCH (heuristic, floating/ball midpoints; nothing here is part of the certificate).  For the K
heaviest cosines we want  gamma_j y + theta_j = 0 (mod 2 pi)  simultaneously.  Lattice with rows
    e_j      :  ( 0, ..., round(s_j 2 pi 2^nu) at j, ..., 0 | 0 | 0 )                j = 1..K
    row_y    :  ( round(s_j gamma_j 2^(nu-ub)) )_j            | EPS | 0 )            y = z 2^-ub
    row_inh  :  ( round(s_j theta_j 2^nu) )_j                 | 0   | M )
with weights s_j = sqrt(A_j/A_1) and M huge.  After LLL the basis vector with last coordinate +-M
gives an integer z with all K phases small (embedding technique for the closest vector problem);
the other short basis vectors give approximate common periods, used to generate more candidates.
Candidates are ranked by the full sum F_Z and polished by a local search in y.

CERTIFICATE.  At the final y_0 (an exact dyadic rational, printed in full) F_Z(y_0) is evaluated in
ball arithmetic from the refined zero brackets (raw/zeros_hp/, t2_refine.py: each bracket has
certified opposite signs of Hardy's Z at its ends and lies inside a bracket of raw/zeros_cache/),
with working precision far beyond log2|y_0|, and the witness is  F_Z(y_0) - r > (2/3) N(q,a).

Usage:  python3 t2_lll.py q a T P [K list] [nu range]      e.g.  python3 t2_lll.py 28 9 300 480 40,60,80 80:400:8
Output: one JSON line (best certified value), log on stderr.
"""
import sys, json, math, time, random
from flint import arb, acb, fmpz, fmpz_mat, fmpq, dirichlet_char, ctx
import t2_witness as TW
import l1const as LC


def load(q, a, T, P):
    """cosines (gamma, A, theta) as Arb balls; remainder r as an Arb ball; L_1"""
    table = TW.char_table(q)
    handled = set()
    cos = []
    def hp(f, m):
        return [arb(x[0]) + arb(0, 1) * 0 + (arb(x[1]) * arb(-1)).union(arb(x[1])) for x in json.load(open(f"raw/zeros_hp/z_{f}_{m}_T{T}_P{P}.json"))]
    for chi, f, m, pr in table:
        key = (f, m)
        cj = (f, pow(m, -1, f)) if f > 1 else (1, 1)
        if key in handled:
            continue
        handled.add(key); handled.add(cj)
        wbar = acb(chi(a)).conjugate() if q > 1 else acb(1)
        if f == 1:
            zs = []
            k = 1
            while True:
                batch = acb.zeta_zeros(k, 50)
                stop = False
                for z in batch:
                    if z.imag > T:
                        stop = True
                        break
                    zs.append(z.imag)
                if stop:
                    break
                k += 50
        elif key == cj:
            zs = hp(f, m)
        else:
            zs = hp(f, m) + [-g for g in hp(*cj)]
        for g in zs:
            rho = acb(arb(1) / 2, g)
            c = -wbar / (rho * (rho + 1))
            cos.append((g, 2 * abs(c), c.arg()))
    # remainder
    ctx_prec = ctx.prec
    ctx.prec = 128
    L1 = arb(0)
    for chi, f, m, pr in table:
        L1 += LC.S_integral(f, m)
    ctx.prec = ctx_prec
    found = arb(0)
    for g, A, th in cos:
        found += A
    rem = L1 - found
    assert rem > -1e-20
    cos.sort(key=lambda t: -float(t[1].mid()))
    return cos, rem, L1


def FZ(cos, y):
    tot = arb(0)
    for g, A, th in cos:
        tot += A * (g * y + th).cos()
    return tot


def FZmid(cos, y):
    return float(FZ(cos, y).mid())


def lll_candidates(cos, K, nu, ub, log):
    """candidate y (arb exact dyadic) from one LLL run; also approximate periods"""
    top = cos[:K]
    A1 = top[0][1]
    two_nu = arb(2) ** nu
    def rnd(x):
        return fmpz(int((x.mid() + arb(1) / 2).floor().unique_fmpz())) if True else None
    n = K + 2
    rows = [[fmpz(0)] * n for _ in range(n)]
    EPS = fmpz(1) << 8
    M = fmpz(1) << (nu + 24)
    s = [(A / A1).sqrt() for g, A, th in top]
    for j, (g, A, th) in enumerate(top):
        rows[j][j] = rnd(s[j] * 2 * arb.pi() * two_nu)
        rows[K][j] = rnd(s[j] * g * arb(2) ** (nu - ub))
        rows[K + 1][j] = rnd(s[j] * th * two_nu)
    rows[K][K] = EPS
    rows[K + 1][K + 1] = M
    Mat = fmpz_mat(n, n, [x for r in rows for x in r])
    L = Mat.lll()
    cands, periods = [], []
    for i in range(n):
        last = L[i, n - 1]
        z = L[i, K]
        assert z % EPS == 0
        z = z // EPS
        if last == 0:
            if z != 0:
                periods.append(arb(z) / arb(2) ** ub)
        elif abs(last) == M:
            sign = 1 if last > 0 else -1
            cands.append(arb(sign * z) / arb(2) ** ub)
    return cands, periods


def main():
    q, a, T, P = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    Ks = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [40, 60, 80]
    lo, hi, st = (int(x) for x in sys.argv[6].split(":")) if len(sys.argv) > 6 else (80, 360, 8)
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    random.seed(seed)
    ctx.prec = P + 64
    t0 = time.time()
    N = sum(1 for b in range(q) if math.gcd(b, q) == 1 and (b * b - a) % q == 0)
    thr = arb(2) * N / 3
    cos, rem, L1 = load(q, a, T, P)
    log = sys.stderr.write
    log(f"({q},{a}): {len(cos)} cosines, L_1 = {L1.str(12)}, remainder r <= {float(rem.upper()):.6f}, threshold {float(thr.mid()):.4f}; need F_Z >= {float((thr + rem).upper()):.4f}  [{time.time()-t0:.0f}s]\n")
    best = (-1e9, None)
    pool = []
    for K in Ks:
        for nu in range(lo, hi + 1, st):
            if nu + 40 > P:
                continue
            for ub in (0, 6):
                try:
                    cands, periods = lll_candidates(cos, K, nu, ub, log)
                except Exception as e:
                    log(f"  LLL failed K={K} nu={nu}: {e}\n")
                    continue
                for y in cands:
                    sK = float(FZ(cos[:K], y).mid())
                    WK = sum(float(c[1].mid()) for c in cos[:K])
                    pool.append((sK, y, K, nu, ub, periods[:6]))
            log(f"  K={K} nu={nu}: pool {len(pool)}, best top-K alignment so far {max(p[0]/sum(float(c[1].mid()) for c in cos[:p[2]]) for p in pool) if pool else 0:.4f}  [{time.time()-t0:.0f}s]\n")
    # rank by the full sum; explore neighbours y + n * period
    pool.sort(key=lambda p: -p[0])
    evaluated = []
    for sK, y, K, nu, ub, periods in pool[:60]:
        f0 = FZmid(cos, y)
        evaluated.append((f0, y, K, nu))
        for per in periods[:4]:
            for nmul in (-3, -2, -1, 1, 2, 3):
                y2 = y + nmul * per
                f2 = FZmid(cos, y2)
                evaluated.append((f2, y2, K, nu))
    evaluated.sort(key=lambda p: -p[0])
    log(f"  best full sums before polish: {[round(e[0], 4) for e in evaluated[:8]]}  [{time.time()-t0:.0f}s]\n")
    # local polish of the best few
    finals = []
    for f0, y, K, nu in evaluated[:6]:
        step = arb(2) ** (-4)
        fy = f0
        for _ in range(36):
            improved = True
            while improved:
                improved = False
                for cand in (y - step, y + step):
                    fc = FZmid(cos, cand)
                    if fc > fy:
                        y, fy, improved = cand, fc, True
            step = step / 2
        # make y an exact dyadic rational with a short description
        ym = y.mid()
        finals.append((fy, ym, K, nu))
        log(f"    polished: F_Z = {fy:.5f} (K={K}, nu={nu})\n")
    finals.sort(key=lambda p: -p[0])
    fy, y0, K, nu = finals[0]
    y0 = arb(y0)                       # exact (a ball of radius 0: midpoint of a ball is a dyadic number)
    assert y0.rad() == 0
    val = FZ(cos, y0)
    lower = val - rem
    ok = bool(lower > thr)
    man, exp = y0.mid().man_exp()
    rec = {"q": q, "a": a, "N": N, "threshold": float(thr.mid()), "T": float(T), "method": "LLL",
           "num_brackets": len(cos), "L1": L1.str(15, radius=False), "remainder_upper": float(rem.upper()),
           "y0": float(y0.mid()), "y0_exact": f"{int(man)} * 2^({int(exp)})", "Ymax": abs(float(y0.mid())),
           "F_Z_at_y0": val.str(15), "F_lower_bound": float(lower.lower()), "witness": ok,
           "K": K, "nu": nu, "P": P, "seconds": round(time.time() - t0)}
    print(json.dumps(rec))
    log(f"({q},{a}): y0 ~ {float(y0.mid()):.6e}, F_Z = {val.str(10)}, lower = {float(lower.lower()):.5f}, threshold {float(thr.mid()):.5f}, witness = {ok}  [{time.time()-t0:.0f}s]\n")


if __name__ == "__main__":
    main()
