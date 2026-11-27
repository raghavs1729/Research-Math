#!/usr/bin/env python3
"""Fail-closed Arb certificate of the ENVELOPE CRITERION (Theorem C / Corollary C'
of monotonicity/BRIDGE-ATTEMPT.md) on whole M-cells.

WHAT IS PROVED PER CELL
=======================
chi = real primitive character mod q, q in {3, 4};  M = m + 1/2;
c_j(M) = 4M/sqrt(M^2 + g_j^2) over the positive ordinates g_1 < g_2 < ... of L(s,chi);
delta(M) = 1/2 + (1/pi) int_0^oo Phi_M(t) sin t / t dt,  Phi_M = prod_j J_0(c_j t).

For a cell [M_a, M_b] the program prints a number LB and the word PASS only if it
has verified, in outward-rounded ball arithmetic, that LB > 0 and

        - M * delta'(M)  >=  MAIN - TAILF - NUM  >=  LB      for EVERY real M in [M_a, M_b],

hence delta'(M) <= -LB/M_b < 0 on the whole cell.  The three quantities are those
of Corollary C' (case (b), the n lowest ordinates):

    MAIN  = (1/pi) int_0^{t_1} exp(-sig2(M_b) t^2/2 - beta(z_1) S4(M_b) t^4) cos t dt
    TAILF = (1/pi) int_{t_1}^oo B(t) dt
    NUM   = S4(M_b)/(32 pi) int_0^oo t^2 B1(t) dt
    B(t)  = prod_{i<=n} b(c_i(M_a) t) * exp(-(v/2) min(t, z_*/c_n(M_a))^2),   B1 = B/b(c_1(M_a) t)
    v     = sig2(M_a) - sum_{i<=n} c_i(M_a)^2/2        (variance of the zeros above g_n at M_a)
    b(z)  = exp(-z^2/4) (z <= z_*), rho (z_* <= z <= z_far), (8/pi) z^{-1/2} (z >= z_far)
    rho = 0.446, z_* = 2 sqrt(log(1/rho)), z_far = (8/(pi rho))^2,
    beta(z_1) = (-log J_0(z_1) - z_1^2/4)/z_1^4,   t_1 <= min(pi/2, z_1/c_1(M_b)).

HYPOTHESES (not proved here; see METHOD.md)
    (H1) GRH(chi).  (H2) variance identity sig2(M) = 4M (log Lambda)'(M+1/2) and its
    M-derivative (paper prop:varid, lem:Dprop).  (H3) the zero files: each decimal is within
    1e-20 of an ordinate, zeros simple, list complete up to its last entry.  (H4) FLINT/Arb.
    (H5) Lemma D, Lemma J, Lemma L, Theorem C, Corollary C' of BRIDGE-ATTEMPT.md (pencil proofs,
    refereed in REFEREE-3.md).

NO FLOAT ON A CERTIFIED PATH.  Floats are used (i) to choose the t-grid, the cut T_end and the
candidate t_1 -- every such choice is an exact rational that is then re-verified in Arb -- and
(ii) to print (outward rounded).  Any undecided comparison, non-finite ball or failed hypothesis
check raises; a cell is PASS only if LB.lower() > 0.  "Fail-closed" refers to these checks; it
does NOT cover a wrong zero file: that is hypothesis (H3).  In series mode a cell compares the zero
list with the closed forms only through v > 0; `selftest` compares S4(series) with the list
enclosure at five values of M and detects gross errors only (REFEREE-5 minor 9).

USAGE
    certify_envelope.py cell  --q 4 --Ma 3.6 --Mb 3.636 [--n 12]
    certify_envelope.py range --q 4 --lo 1.5 --hi 400 --ratio 1.01 --out results_q4.txt
    certify_envelope.py point --q 4 --M 6          (degenerate cell M_a = M_b, for tables)
    certify_envelope.py selftest
"""
from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys
import time
from fractions import Fraction

from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "cert")))
# Re-used from the direct certificate: zero loading (ball = decimal +/- 1e-20, strict ordering
# checked), the Dirichlet-series closed forms of (log Lambda)', (log Lambda)'' for s >= 3, the
# complex-step (log Lambda)' for s >= 1, exact-rational helpers and outward-rounded printing.
from certify_derivative import (  # noqa: E402
    dec, e_dn, e_up, ex, fr, load_zeros, loglam_d1_complex_step, loglam_derivs_series, make_grid,
)

PREC = 128
RHO = Fraction(446, 1000)
Z1_CANDIDATES = [Fraction(k, 20) for k in (8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34,
                                           36, 38, 40, 42, 44, 46)]      # all < 2.4 < j_{0,1}
N_FINE = 3000            # fine t-grid: N_FINE steps on [0, FINE_SPAN/sigma_1]
FINE_SPAN = 12.0
GEO_RATIO = 1.01         # geometric t-grid beyond
DYAD = 1 << 24           # t-grid points are multiples of 2^-24 (exact in Arb)


class Fail(ArithmeticError):
    pass


def lo(x: arb) -> arb:
    v = x.lower()
    if not v.is_finite():
        raise Fail("non-finite lower endpoint")
    return v


def hi(x: arb) -> arb:
    v = x.upper()
    if not v.is_finite():
        raise Fail("non-finite upper endpoint")
    return v


_CONST = {}


def consts():
    """Balls for z_*, z_far and the exact rho (computed once per precision)."""
    if ctx.prec not in _CONST:
        rho = ex(RHO)
        zs = 2 * (-(rho.log())).sqrt()                 # z_* = 2 sqrt(log(1/rho))
        zfar = (8 / (arb.pi() * rho)) ** 2
        if not (zs > 1 and zs < arb(24) / 10):
            raise Fail("z_* out of range")
        _CONST[ctx.prec] = dict(rho=rho, zs=zs, zs_lo=lo(zs), zs_hi=hi(zs), zfar_hi=hi(zfar),
                                c8pi=8 / arb.pi())
    return _CONST[ctx.prec]


# --------------------------------------------------------------------------
# closed forms
# --------------------------------------------------------------------------
def c_ball(g: arb, M: arb) -> arb:
    """c(M) = 4/sqrt(1 + (g/M)^2); M occurs once."""
    r = g / M
    return 4 / (1 + r * r).sqrt()


def sig2_S4(q: int, M: Fraction, zeros_all: list[arb]) -> tuple[arb, arb, str]:
    """Enclosures of sig2(M) = sum c^2/2 and S4(M) = sum c^4 at the exact rational M.

    s = M + 1/2 >= 4:  sig2 = 4M h, S4 = 64 M h - 64 M^2 h',  h = (log Lambda)'(s),
                       h' = (log Lambda)''(s)  (Dirichlet series with rigorous tails).
    1 <= s < 4:        sig2 = 4M h by the complex step;  S4 is enclosed WITHOUT h' by
                       sum_{list} c^4 <= S4 <= sum_{list} c^4 + c_N^2 * 2 V_T,
                       V_T = sig2 - sum_{list} c^2/2 >= 0 (the zeros above the list have c <= c_N;
                       the list is complete, (H3)).
    """
    Mb = ex(M)
    s = Mb + arb(1) / 2
    if M + Fraction(1, 2) >= 4:
        d1, d2 = loglam_derivs_series(q, s)
        sig2 = 4 * Mb * d1
        S4 = 64 * Mb * d1 - 64 * Mb * Mb * d2
        how = "series"
    else:
        d1 = loglam_d1_complex_step(q, s)
        sig2 = 4 * Mb * d1
        s2 = arb(0)
        s4 = arb(0)
        for g in zeros_all:
            c2 = c_ball(g, Mb) ** 2
            s2 += c2
            s4 += c2 * c2
        VT = sig2 - s2 / 2
        if not hi(VT) > 0:
            raise Fail("tail variance not positive")
        cN2 = c_ball(zeros_all[-1], Mb) ** 2
        extra = hi(cN2) * 2 * hi(VT)
        S4 = lo(s4).union(hi(s4) + extra)
        how = "list"
    if not (lo(sig2) > 0 and lo(S4) > 0 and sig2.is_finite() and S4.is_finite()):
        raise Fail("closed forms not certified positive/finite")
    return sig2, S4, how


def beta_hi(z1: Fraction) -> arb:
    """Upper endpoint of beta(z1) = (-log J_0(z1) - z1^2/4)/z1^4, 0 < z1 < 2.4 < j_{0,1}."""
    if not (0 < z1 < Fraction(24, 10)):
        raise Fail("z1 out of range")
    z = ex(z1)
    J = z.bessel_j(0)
    if not lo(J) > 0:
        raise Fail("J_0(z1) not certified positive")
    b = (-(J.log()) - z * z / 4) / z ** 4
    if not lo(b) > 0:
        raise Fail("beta not positive")
    return hi(b)


# --------------------------------------------------------------------------
# the cell
# --------------------------------------------------------------------------
def dyad_floor(x: float) -> Fraction:
    return Fraction(math.floor(x * DYAD), DYAD)


def dyad_ceil(x: float) -> Fraction:
    return Fraction(math.ceil(x * DYAD), DYAD)


def certify_cell(q: int, Ma: Fraction, Mb: Fraction, n: int = 12) -> dict:
    t0 = time.time()
    ctx.prec = PREC
    K = consts()
    if not (0 < Ma <= Mb):
        raise Fail("bad cell")
    if n < 10:
        raise Fail("n >= 10 required (far remainder uses ordinates 2..10)")
    zeros_all = load_zeros(q)
    zeros = zeros_all[:n]
    A, Bm = ex(Ma), ex(Mb)

    # ---- amplitudes -------------------------------------------------------
    ca = [c_ball(g, A) for g in zeros]                 # c_i(M_a), i = 1..n
    ca_lo = [lo(c) for c in ca]
    for a, b in zip(ca_lo, ca_lo[1:]):
        if not (a > 0 and b > 0):
            raise Fail("amplitude not positive")
    cn_hi = hi(ca[-1])
    c1b_hi = hi(c_ball(zeros[0], Bm))                  # c_1(M_b) = max_j c_j(M) over the cell

    # ---- closed forms ------------------------------------------------------
    sig2_a, _S4a, how_a = sig2_S4(q, Ma, zeros_all)
    if Mb == Ma:
        sig2_b, S4_b, how_b = sig2_a, _S4a, how_a
    else:
        sig2_b, S4_b, how_b = sig2_S4(q, Mb, zeros_all)
    sig2b_hi, S4b_hi = hi(sig2_b), hi(S4_b)
    head = arb(0)
    for c in ca:
        head += c * c / 2
    v = sig2_a - head                                  # variance above g_n at M_a
    v_lo = lo(v)
    if not v_lo > 0:
        raise Fail("v_n(M_a) not certified positive")
    r_lo = lo(K["zs_lo"] / cn_hi)                      # lower bound of z_*/c_n(M_a)

    # ---- t-grid (exact dyadic rationals; float steering only) --------------
    sig1 = math.sqrt(max(float(sig2_a.mid()) - float(ca[0].mid()) ** 2 / 2, 1e-6))
    span = FINE_SPAN / sig1
    h = dyad_ceil(span / N_FINE)
    grid = [h * k for k in range(N_FINE + 1)]
    c10 = float(ca_lo[9].mid())
    Tend = dyad_ceil(float(K["zfar_hi"].mid()) / c10 * 1.0001)
    if Tend <= grid[-1]:
        Tend = grid[-1] * 2
    while grid[-1] < Tend:
        nxt = dyad_ceil(float(grid[-1]) * GEO_RATIO)
        grid.append(min(nxt, Tend) if nxt < Tend else Tend)
    for a, b in zip(grid, grid[1:]):
        if not a < b:
            raise Fail("t-grid not increasing")
    tb = [ex(t) for t in grid]
    for t in tb:
        if not t.is_exact():
            raise Fail("t-grid point not exact")
    TendB = tb[-1]
    if not lo(ca_lo[9] * TendB) > K["zfar_hi"]:
        raise Fail("T_end too small for the far remainder")

    # ---- envelope values at the grid points --------------------------------
    zs_lo, zs_hi, zfar_hi, rho = K["zs_lo"], K["zs_hi"], K["zfar_hi"], K["rho"]
    log_rho = rho.log()
    exp_amb = zs_lo * zs_lo / 4

    def log_b_upper(c_lo: arb, t: arb) -> arb:
        """Ball whose upper endpoint bounds log b(c t) for every c >= c_lo."""
        z = lo(c_lo * t)                               # exact lower bound of the argument
        if z < zs_lo:
            return -(z * z) / 4
        if z > zfar_hi:
            return (K["c8pi"]).log() - z.log() / 2
        if z > zs_hi:
            return log_rho
        return -exp_amb                                # z within the ball of z_*: b <= exp(-zs_lo^2/4)

    B_hi = []        # upper bounds of B(t_k)
    B1_hi = []       # upper bounds of B1(t_k)
    for t in tb:
        L1 = arb(0)
        for c in ca_lo[1:]:
            L1 += log_b_upper(c, t)
        m = t if t < r_lo else r_lo                    # exact min (both exact)
        L1 -= v_lo * m * m / 2
        L0 = L1 + log_b_upper(ca_lo[0], t)
        B1_hi.append(hi(L1.exp()))
        B_hi.append(hi(L0.exp()))

    # ---- remainders beyond T_end -------------------------------------------
    pref = K["c8pi"] ** 9
    for c in ca_lo[1:10]:
        pref = pref / c.sqrt()
    mT = TendB if TendB < r_lo else r_lo
    pref = pref * (-(v_lo * mT * mT / 2)).exp()
    rem2 = hi(pref * TendB ** arb("-1.5") / arb("1.5"))      # int_T^oo t^2 t^{-9/2} dt
    rem0 = hi(pref * TendB ** arb("-3.5") / arb("3.5"))      # int_T^oo     t^{-9/2} dt

    # ---- NUM (upper) ---------------------------------------------------------
    acc = arb(0)
    for k in range(len(tb) - 1):
        acc += tb[k + 1] * tb[k + 1] * B1_hi[k] * (tb[k + 1] - tb[k])
    num_int = hi(acc + rem2)
    NUM = hi(S4b_hi / (32 * arb.pi()) * num_int)

    # ---- suffix sums of B (upper), for TAILF from any grid point -------------
    suf = [arb(0)] * len(tb)
    run = arb(0) + rem0
    suf[-1] = run            # int_{T_end}^oo B <= rem0 (REFEREE-5 minor 3; index unreachable as t_1 < pi/2 < T_end)
    for k in range(len(tb) - 2, -1, -1):
        run = run + B_hi[k] * (tb[k + 1] - tb[k])
        suf[k] = run

    # ---- candidates for (z_1, t_1); MAIN (lower) -----------------------------
    half_pi_lo = lo(arb.pi() / 2)
    best = None
    c1f = float(c1b_hi.mid())
    for z1 in Z1_CANDIDATES:
        cap = min(math.pi / 2, float(z1) / c1f) * (1 - 1e-12)
        k1 = 0
        # largest grid index with t_k <= cap (float), then verified in Arb
        lo_i, hi_i = 0, len(grid) - 1
        while lo_i < hi_i:
            mid = (lo_i + hi_i + 1) // 2
            if float(grid[mid]) <= cap:
                lo_i = mid
            else:
                hi_i = mid - 1
        k1 = lo_i
        if k1 < 1:
            continue
        t1 = tb[k1]
        if not (t1 < half_pi_lo and hi(c1b_hi * t1) <= ex(z1)):
            raise Fail("t_1 not certified admissible")
        bh = beta_hi(z1)
        coef2 = sig2b_hi / 2
        coef4 = bh * S4b_hi
        s = arb(0)
        for k in range(1, k1 + 1):
            t = tb[k]
            t2 = t * t
            G = lo((-(coef2 * t2) - coef4 * t2 * t2).exp() * t.cos())
            if not G > 0:
                raise Fail("MAIN integrand not positive")
            s += G * (t - tb[k - 1])                 # lower Riemann sum (G decreasing on [0, pi/2])
            if float(G.mid()) < 1e-40:
                break
        main = lo(s / arb.pi())
        tailf = hi(suf[k1] / arb.pi())
        val = lo(main - tailf)
        if best is None or val > best[0]:
            best = (val, main, tailf, z1, grid[k1])
    if best is None:
        raise Fail("no admissible t_1")
    val, MAIN, TAILF, z1, t1 = best
    LB = lo(val - NUM)
    ok = bool(LB > 0)
    Kratio = hi(NUM / val) if val > 0 else None
    return dict(q=q, Ma=Ma, Mb=Mb, n=n, ok=ok, LB=LB, MAIN=MAIN, TAILF=TAILF, NUM=NUM, K=Kratio,
                dp_hi=hi(-LB / Bm) if ok else None, z1=z1, t1=t1, sig2_b=sig2_b, S4_b=S4_b,
                v_lo=v_lo, rem2=rem2, Tend=Tend, npts=len(tb), how=how_a + "/" + how_b,
                secs=time.time() - t0)


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------
def upf(x: arb) -> float:
    return math.nextafter(float(x.upper()), math.inf)


def dnf(x: arb) -> float:
    return math.nextafter(float(x.lower()), -math.inf)


def fmt_line(r: dict) -> str:
    """LB, MAIN rounded DOWN; TAILF, NUM, K, dprime_hi rounded UP."""
    k = "inf" if r["K"] is None else e_up(upf(r["K"]), 4)
    dp = "nan" if r["dp_hi"] is None else e_up(upf(r["dp_hi"]))
    return (f"{dec(r['Ma'])} {dec(r['Mb'])} LB={e_dn(dnf(r['LB']))} MAIN={e_dn(dnf(r['MAIN']))} "
            f"TAILF={e_up(upf(r['TAILF']), 3)} NUM={e_up(upf(r['NUM']))} K={k} dprime_hi={dp} "
            f"sig2b={e_up(upf(r['sig2_b']), 5)} S4b={e_up(upf(r['S4_b']), 5)} v={e_dn(dnf(r['v_lo']), 4)} "
            f"z1={float(r['z1']):.2f} t1={float(r['t1']):.5f} Tend={float(r['Tend']):.2f} "
            f"rem={e_up(upf(r['rem2']), 2)} npts={r['npts']} n={r['n']} cf={r['how']} "
            f"secs={r['secs']:.1f} {'PASS' if r['ok'] else 'FAIL'}")


def _work(job):
    q, Ma, Mb, n = job
    try:
        return job, certify_cell(q, Ma, Mb, n), None
    except Exception as e:      # fail closed
        return job, None, f"{type(e).__name__}: {e}"


def input_hashes(q: int) -> str:
    """SHA-256 of this prover, of the shared module and of the zero file (REFEREE-5 minor 8)."""
    import certify_derivative as cd
    files = (("prover", os.path.abspath(__file__)), ("shared", os.path.abspath(cd.__file__)),
             ("zeros", os.path.join(cd.ZERO_DIR, cd.ZERO_FILES[q])))
    parts = []
    for name, path in files:
        with open(path, "rb") as fh:
            parts.append(f"{name}={hashlib.sha256(fh.read()).hexdigest()}")
    return "# sha256 " + " ".join(parts)


def done_cells(path: str, hashline: str) -> set:
    """Cells already PASS in an existing results file.  The file is accepted for resuming only
    if its header carries exactly the current input hashes; otherwise the run is refused."""
    out = set()
    if os.path.exists(path):
        with open(path, encoding="ascii") as fh:
            lines = fh.readlines()
        if hashline + "\n" not in lines:
            raise SystemExit(f"{path}: produced by a different prover / shared module / zero file "
                             "(or has no hash header); delete it to recompute")
        for line in lines:
            p = line.split()
            if len(p) > 3 and not line.startswith("#") and p[-1] == "PASS":
                out.add((Fraction(p[0]), Fraction(p[1])))
    return out


def run_range(q, lo_, hi_, ratio, out, workers, n):
    import multiprocessing as mp
    grid = make_grid(fr(lo_), fr(hi_), ratio)
    hashline = input_hashes(q)
    done = done_cells(out, hashline)
    jobs = [(q, a, b, n) for a, b in zip(grid, grid[1:]) if (a, b) not in done]
    new = not os.path.exists(out)
    with open(out, "a", encoding="ascii") as fh:
        if new:
            fh.write(f"# certify_envelope.py range q={q} [{lo_}, {hi_}] ratio={ratio} n={n} prec={PREC}\n")
            fh.write(hashline + "\n")
            fh.write("# per cell: for ALL M in [M_a,M_b]:  -M delta'(M) >= MAIN - TAILF - NUM >= LB > 0, "
                     "hence delta'(M) <= dprime_hi < 0.  K = NUM/(MAIN-TAILF) (upper bound).\n")
        t0 = time.time()
        nfail = 0
        if workers > 1:
            pool = mp.Pool(workers)
            it = pool.imap(_work, jobs)
        else:
            it = map(_work, jobs)
        for job, r, err in it:
            if r is None:
                fh.write(f"{dec(job[1])} {dec(job[2])} ERROR {err} FAIL\n")
                nfail += 1
            else:
                fh.write(fmt_line(r) + "\n")
                nfail += 0 if r["ok"] else 1
            fh.flush()
        if workers > 1:
            pool.close()
            pool.join()
    summarize(out, fr(lo_), fr(hi_))
    print(f"range q={q}: {len(jobs)} cells run, {nfail} not certified, {time.time() - t0:.0f}s")


def summarize(out: str, lo_=None, hi_=None) -> None:
    rows = []
    with open(out, encoding="ascii") as fh:
        for line in fh:
            p = line.split()
            if line.startswith("#") or len(p) < 4:
                continue
            ok = p[-1] == "PASS"
            d = dict(x.split("=", 1) for x in p[2:-1] if "=" in x)
            rows.append((Fraction(p[0]), Fraction(p[1]), ok, d))
    rows.sort(key=lambda r: r[0])
    good = [r for r in rows if r[2]]
    bad = [r for r in rows if not r[2]]
    # contiguity of the certified cells
    runs = []
    for a, b, _, _ in good:
        if runs and runs[-1][1] == a:
            runs[-1][1] = b
        elif runs and runs[-1][1] > a:
            runs[-1][1] = max(runs[-1][1], b)
        else:
            runs.append([a, b])
    with open(out, "a", encoding="ascii") as fh:
        fh.write(f"# SUMMARY {time.strftime('%Y-%m-%d %H:%M:%S')}: {len(good)} PASS cells, {len(bad)} FAIL/ERROR lines\n")
        fh.write("# SUMMARY certified union: " + " U ".join(f"[{dec(a)}, {dec(b)}]" for a, b in runs) + "\n")
        if good:
            wk = max(good, key=lambda r: float(r[3]["K"]))
            wl = min(good, key=lambda r: float(r[3]["LB"]))
            wd = max(good, key=lambda r: float(r[3]["dprime_hi"]))
            fh.write(f"# SUMMARY worst K (upper bound) = {wk[3]['K']} on [{dec(wk[0])}, {dec(wk[1])}]\n")
            fh.write(f"# SUMMARY smallest LB = {wl[3]['LB']} on [{dec(wl[0])}, {dec(wl[1])}]\n")
            fh.write(f"# SUMMARY weakest certified bound delta' <= {wd[3]['dprime_hi']} on [{dec(wd[0])}, {dec(wd[1])}]\n")


def selftest() -> None:
    """Consistency checks of the closed forms (two independent evaluations each)."""
    ctx.prec = PREC
    for q in (4, 3):
        zs = load_zeros(q)
        for M in ("3.5", "5", "10", "40", "150"):
            Mf = fr(M)
            Mb = ex(Mf)
            s = Mb + arb(1) / 2
            d1, d2 = loglam_derivs_series(q, s)
            d1c = loglam_d1_complex_step(q, s)
            if not d1.overlaps(d1c):
                raise Fail("series vs complex step disagree")
            sig2 = 4 * Mb * d1
            S4 = 64 * Mb * d1 - 64 * Mb * Mb * d2
            s2 = arb(0)
            s4 = arb(0)
            for g in zs:
                c2 = c_ball(g, Mb) ** 2
                s2 += c2
                s4 += c2 * c2
            VT = sig2 - s2 / 2
            S4list = lo(s4).union(hi(s4) + hi(c_ball(zs[-1], Mb) ** 2) * 2 * hi(VT))
            if not (lo(VT) > 0 and S4.overlaps(S4list)):
                raise Fail("S4 closed form inconsistent with the list")
            print(f"q={q} M={M}: sig2={sig2.str(15)}  S4(series)={S4.str(15)}  "
                  f"S4(list) in [{float(lo(S4list)):.12f}, {float(hi(S4list)):.12f}]  "
                  f"V_T={VT.str(8)}  (series vs complex step: overlap)")
    K = consts()
    print("z_* =", K["zs"].str(20), " z_far <=", K["zfar_hi"].str(12))
    print("selftest OK")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cell"); c.add_argument("--q", type=int, required=True)
    c.add_argument("--Ma", required=True); c.add_argument("--Mb", required=True)
    c.add_argument("--n", type=int, default=12)
    p = sub.add_parser("point"); p.add_argument("--q", type=int, required=True)
    p.add_argument("--M", required=True, nargs="+"); p.add_argument("--n", type=int, default=12)
    r = sub.add_parser("range"); r.add_argument("--q", type=int, required=True)
    r.add_argument("--lo", required=True); r.add_argument("--hi", required=True)
    r.add_argument("--ratio", type=float, default=1.01); r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=1); r.add_argument("--n", type=int, default=12)
    s = sub.add_parser("summarize"); s.add_argument("--out", required=True)
    sub.add_parser("selftest")
    a = ap.parse_args()
    if a.cmd == "cell":
        print(fmt_line(certify_cell(a.q, fr(a.Ma), fr(a.Mb), a.n)))
    elif a.cmd == "point":
        for M in a.M:
            print(fmt_line(certify_cell(a.q, fr(M), fr(M), a.n)))
    elif a.cmd == "range":
        if a.workers > 2:
            raise SystemExit("at most 2 workers")
        run_range(a.q, a.lo, a.hi, a.ratio, a.out, a.workers, a.n)
    elif a.cmd == "summarize":
        summarize(a.out)
    elif a.cmd == "selftest":
        selftest()


if __name__ == "__main__":
    main()
