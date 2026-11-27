#!/usr/bin/env python3
"""
t2_zeros.py -- zero-side witnesses for HKMN23.md, Theorems 3 and 6, with a closed-form remainder.

Under GRH(q):
    F_{q,a}(y) = - sum_{chi mod q} conj(chi(a)) sum_{rho = 1/2 + i gamma of chi*} e^{i gamma y} / (rho (rho+1)).

A finite set Z of zeros is located by certified sign changes of Hardy's Z-function (Arb balls of
opposite certified signs at the two ends of brackets; each bracket then contains at least one zero
on the critical line, and one zero per bracket is counted; the brackets of one character are
checked to be pairwise disjoint when they are loaded, see _load).  NO completeness is claimed
or needed: the total weight S(chi*) = sum_rho 1/|rho(rho+1)| is known exactly (l1const.S_integral,
closed form), so the zeros not in Z have total weight
    r = sum_chi [ S(chi*) - sum_{rho in Z_chi} 1/|rho(rho+1)| ]   >= 0,
and for every real y
    F_{q,a}(y) >= F_Z(y) - r,       F_Z = the same sum restricted to Z.
A witness is a y_0 with F_Z(y_0) - r > (2/3) N(q,a), certified in ball arithmetic.

Zeros of chi-bar are the negatives of the zeros of chi, so only ordinates in (0, T] are searched,
for every primitive character involved; zeta zeros come from Arb (acb.zeta_zeros).

Usage:  python3 t2_zeros.py q T Ymax dy a1,a2,... [TSCAN]
        (classes a; "sq" = all square classes not in the LWS list; TSCAN: the floating-point scan
         uses only zeros with |gamma| <= TSCAN, default T; certification always uses all zeros)
Output: one JSON line per class, plus a log on stderr.
"""
import sys, json, math, time, os
import numpy as np
from flint import acb, arb, dirichlet_char, ctx
import t2_witness as TW
import l1const as LC

ctx.prec = 128
H = 0.04
TSCAN = 1e9


def hardy_sign(chi, t):
    v = chi.hardy_z(acb(t)).real
    if v > 0:
        return 1
    if v < 0:
        return -1
    return 0


def _load(fn):
    """read cached brackets and CHECK that they are pairwise disjoint and positive (so that the
    zeros counted, one per bracket, are distinct)"""
    zs = [arb(x) for x in json.load(open(fn))]
    for i in range(len(zs)):
        assert zs[i] > 0, (fn, i)
        if i:
            assert zs[i - 1] < zs[i], (fn, i, "overlapping brackets")
    return zs


def positive_zeros(f, m, T):
    """list of arb balls gamma in (0, T], each a certified bracket [lo, hi] containing a zero of L(s, chi_f(m,.)).
    Results are cached in raw/zeros_cache/ (decimal strings of outward-rounded balls)."""
    os.makedirs("raw/zeros_cache", exist_ok=True)
    fn = f"raw/zeros_cache/z_{f}_{m}_T{int(T)}.json"
    if os.path.exists(fn):
        return _load(fn)
    chi = dirichlet_char(f, m)
    out = []
    def Zf(t):
        ctx.prec = 64
        v = float(chi.hardy_z(acb(t)).real.mid())
        ctx.prec = 128
        return v
    t = 1e-9
    prev = Zf(t)
    n = int(T / H)
    for i in range(1, n + 1):
        t2 = i * H
        cur = Zf(t2)
        if prev * cur < 0:
            lo, hi, flo = t, t2, prev
            for _ in range(34):
                mid = 0.5 * (lo + hi)
                fm = Zf(mid)
                if fm == 0.0:
                    break
                if flo * fm < 0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            # widen slightly and certify in ball arithmetic (128 bits)
            w = max(hi - lo, 1e-11)
            lo2, hi2 = lo - w, hi + w
            ok = False
            for _ in range(6):
                s1, s2 = hardy_sign(chi, lo2), hardy_sign(chi, hi2)
                if s1 * s2 < 0:
                    ok = True
                    break
                w *= 8
                lo2, hi2 = lo - w, hi + w
            if ok and lo2 >= t and hi2 <= t2 and lo2 > 0:      # bracket inside its own scan cell
                out.append(arb(lo2).union(arb(hi2)))
        t, prev = t2, cur
    json.dump([x.str(40) for x in out], open(fn, "w"))
    return _load(fn)


def main():
    global TSCAN
    q, T, Ymax, dy = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    TSCAN = float(sys.argv[6]) if len(sys.argv) > 6 else T
    table = TW.char_table(q)
    phi = len(table)
    sqcount = {}
    for b in range(1, q):
        if math.gcd(b, q) == 1:
            sqcount[b * b % q] = sqcount.get(b * b % q, 0) + 1
    if sys.argv[5] == "sq":
        classes = [a for a in sorted(sqcount) if (q, a) not in TW.LWS18]
    else:
        classes = [int(x) for x in sys.argv[5].split(",")]
    t0 = time.time()
    prim = sorted({(f, m) for _, f, m, _ in table})
    pos = {}
    for (f, m) in prim:
        if f == 1:
            zs = []
            k = 1
            while True:
                batch = acb.zeta_zeros(k, 50)
                done = False
                for z in batch:
                    if float(z.imag.mid()) > T:
                        done = True
                        break
                    zs.append(z.imag)
                if done:
                    break
                k += 50
            pos[(f, m)] = zs
        else:
            pos[(f, m)] = positive_zeros(f, m, T)
    sys.stderr.write(f"q={q}: zeros located for {len(prim)} primitive characters, T={T}: {sum(len(v) for v in pos.values())} brackets, {time.time()-t0:.0f}s\n")

    def conj_of(f, m):
        return (f, pow(m, -1, f)) if f > 1 else (1, 1)

    def w_of(g):
        return 1 / ((arb(1) / 4 + g * g) * (arb(9) / 4 + g * g)).sqrt()

    # remainder
    rem = arb(0)
    Stot = arb(0)
    for chi, f, m, pr in table:
        S = LC.S_integral(f, m)
        found = arb(0)
        for g in pos[(f, m)]:
            found += w_of(g)
        for g in pos[conj_of(f, m)]:
            found += w_of(g)
        d = S - found
        assert d > -1e-20, (f, m, S, found)
        rem += d
        Stot += S
    rem_up = float(rem.upper())
    sys.stderr.write(f"q={q}: L_1 = {Stot.str(12)}, remainder r <= {rem_up:.6f}\n")

    for a in classes:
        N = sqcount.get(a, 0)
        thr = 2.0 * N / 3
        # float model: F_Z(y) = Re sum_j c_j exp(i g_j y)
        gs, cs = [], []
        for chi, f, m, pr in table:
            w = complex(acb(chi(a)).conjugate()) if q > 1 else 1.0
            for g in pos[(f, m)]:
                gg = float(g.mid())
                gs.append(gg); cs.append(-w / ((0.5 + 1j * gg) * (1.5 + 1j * gg)))
            for g in pos[conj_of(f, m)]:
                gg = -float(g.mid())
                gs.append(gg); cs.append(-w / ((0.5 + 1j * gg) * (1.5 + 1j * gg)))
        gs = np.array(gs); cs = np.array(cs)
        best = (-1e9, 0.0)
        ys = np.arange(-Ymax, Ymax + dy, dy)
        # the floating-point SCAN may use only the zeros below TSCAN (a proposal step; the
        # certificate below always uses all located zeros)
        sel = np.abs(gs) <= TSCAN
        gsel, csel = gs[sel], cs[sel]
        CH = max(1, int(4e7 // max(1, len(gsel))))
        cands = []
        for i in range(0, len(ys), CH):
            yy = ys[i:i + CH]
            vals = (np.exp(1j * np.outer(yy, gsel)) @ csel).real
            for j in np.argsort(vals)[-3:]:
                cands.append((float(vals[j]), float(yy[j])))
        cands.sort(reverse=True)
        for _, yc in cands[:40]:
            # local refinement with all zeros
            y1, st = yc, dy
            f1 = float((np.exp(1j * gs * y1) @ cs).real)
            for _ in range(40):
                for cand in (y1 - st, y1 + st):
                    fc = float((np.exp(1j * gs * cand) @ cs).real)
                    if fc > f1:
                        y1, f1 = cand, fc
                st /= 2
            if f1 > best[0]:
                best = (f1, y1)
        # local refinement
        def Ff(y):
            return float((np.exp(1j * gs * y) @ cs).real)
        y0 = best[1]
        step = dy
        for _ in range(60):
            for cand in (y0 - step, y0 + step):
                if Ff(cand) > Ff(y0):
                    y0 = cand
            step /= 2
        y0 = round(y0, 9)
        # certify at y0
        Y = arb(str(y0)) if False else arb(y0)
        tot = acb(0)
        for chi, f, m, pr in table:
            w = acb(chi(a)).conjugate() if q > 1 else acb(1)
            sub = acb(0)
            for g in pos[(f, m)]:
                sub += (acb(0, 1) * g * Y).exp() / ((acb(arb(1) / 2, g)) * (acb(arb(3) / 2, g)))
            for g in pos[conj_of(f, m)]:
                g2 = -g
                sub += (acb(0, 1) * g2 * Y).exp() / ((acb(arb(1) / 2, g2)) * (acb(arb(3) / 2, g2)))
            tot -= w * sub
        FZ = tot.real
        lower = FZ - rem
        ok = bool(lower > arb(2) * N / 3)
        rec = {"q": q, "a": a, "phi": phi, "N": N, "threshold": thr, "T": T, "Ymax": Ymax, "dy": dy,
               "num_brackets": int(len(gs)), "L1": Stot.str(15, radius=False), "remainder_upper": rem_up,
               "y0": y0, "F_Z_at_y0": FZ.str(15), "F_lower_bound": float(lower.lower()), "witness": ok,
               "imag_part_check": tot.imag.str(5)}
        print(json.dumps(rec))
        sys.stdout.flush()
        sys.stderr.write(f"  ({q},{a}): N={N} thr={thr:.4f} y0={y0} F_Z={float(FZ.mid()):.5f} lower={float(lower.lower()):.5f} witness={ok}  [{time.time()-t0:.0f}s]\n")


if __name__ == "__main__":
    main()
