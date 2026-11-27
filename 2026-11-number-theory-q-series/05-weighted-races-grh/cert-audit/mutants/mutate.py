#!/usr/bin/env python3
"""AUDITOR: mutation tests of certify_derivative.py (on the auditor's own copy).

The copy `certify_derivative_copy.py` is byte-identical to the audited file
(checked against the sha256 passed on the command line by replay.sh).  Each
mutation is a textual patch of that copy (every pattern must occur exactly
once) or a patched zero file; the mutated program certifies test cells and the
outcome is classified:

  RAISES / FAIL      the mutated checker refuses to certify           (fail-closed)
  PASS, CONTAINS     prints PASS; its [lo,hi] still contains the true delta'
                     at M_a, (M_a+M_b)/2, M_b (the mutation is absorbed or
                     is numerically invisible at the resolution of the output)
  PASS, EXCLUDES     prints PASS with an enclosure that does NOT contain the
                     true value: a silently wrong certificate, visible only to
                     an independent evaluation ("killed by cross-check")

"True value" = the auditor's independent rigorous point enclosure
(indep_enclose.py, width ~1e-8), on the UNMUTATED zero lists.

Usage: mutate.py list | one NAME q Ma Mb | all [--out FILE] | art | driver
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
import types
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "certify_derivative_copy.py")
ZDIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "prototype", "explore"))
sys.path.insert(0, os.path.dirname(HERE))

CELLS = [(4, "3.6", "3.618"), (4, "17.3701", "17.4"), (3, "19.1282", "19.2")]

# name -> (description, [(old, new), ...])           code mutations
CODE = {
    "baseline": ("unmutated copy", []),
    "tail_sign_flip": ("flip the sign of the tail term: dP - P*BT -> dP + P*BT",
                       [("F = sc * Psi * (dP - P * BT)", "F = sc * Psi * (dP + P * BT)")]),
    "dP_sign_flip": ("flip the sign of the head derivative dP",
                     [("dP = -t * acc", "dP = t * acc")]),
    "sinc_sign_flip": ("truth known POSITIVE: delta -> 1 - delta (sinc -> -sinc)",
                       [("sc = t.sinc()", "sc = -t.sinc()")]),
    "J1_to_J2": ("replace J_1 by J_2 in bessel01", [("j1 = x0.bessel_j(1)", "j1 = x0.bessel_j(2)")]),
    "J1_to_x_over_2": ("replace J_1(x) by its first Taylor term x/2", [("j1 = x0.bessel_j(1)", "j1 = x0 / 2")]),
    "J1_to_minus_J1": ("replace J_1 by -J_1", [("j1 = x0.bessel_j(1)", "j1 = -x0.bessel_j(1)")]),
    "J0_to_cos": ("replace J_0 by cos in bessel01", [("j0 = x0.bessel_j(0)", "j0 = x0.cos()")]),
    "bessel_radius_div1000": ("shrink the Bessel Taylor radii by 1000",
                              [("e0 = abs(j1) * r + hr2", "e0 = (abs(j1) * r + hr2) / 1000"),
                               ("e1 = abs(d1) * r + hr2", "e1 = (abs(d1) * r + hr2) / 1000")]),
    "bessel_drop_second_order": ("drop the r^2/2 Lagrange term of the Bessel enclosures",
                                 [("hr2 = r * r / 2", "hr2 = arb(0)")]),
    "t_midpoint_sampling": ("evaluate the integrand at the t-midpoint instead of on the t-box",
                            [("t = hull(ta, tb)", "t = (ta + tb) / 2")]),
    "M_left_endpoint_only": ("use M = M_a instead of the M-ball for all head coefficients",
                             [("M = hull(a, b)\n        self.M = M", "M = a\n        self.M = M")]),
    "V_left_endpoint_only": ("tail variance/derivative taken at M_a only",
                             [("self.V = hull(VTa, VTb)", "self.V = VTa"),
                              ("self.DT = hull(DTa, DTb)", "self.DT = DTa")]),
    "zero_radius_removed": ("zero balls without the 1e-20 radius",
                            [('arb(line.strip(), "1e-20")', "arb(line.strip())")]),
    "E2_dropped": ("drop the large-t remainder E2", [("enc = (I1 + rad_ball(E2)) / pi", "enc = I1 / pi")]),
    "tail_model_radius_dropped": ("drop both hull factors of the tail model (Psi = G, BT = BT0)",
                                  [("Psi = G * hull(em, one)", "Psi = G"),
                                   ("BT = BT0 * hull(one, 1 + self.eta2 * t2)", "BT = BT0")]),
    "tail_model_wrong_direction": ("quartic correction with the wrong sign: Psi in G*[1, e^{+C4S4t^4}]",
                                   [("em = (-self.C4S4 * t2 * t2).exp()", "em = (self.C4S4 * t2 * t2).exp()")]),
    "gaussian_tail_dropped": ("drop the unlisted zeros altogether (G = 1)",
                              [("G = (-self.V * t2 / 2).exp()", "G = arb(1)")]),
    "tail_variance_x3": ("tail variance and its derivative multiplied by 3",
                         [("    VT = sig2 - hv\n    DT = D - hd", "    VT = (sig2 - hv) * 3\n    DT = (D - hd) * 3")]),
    "digamma_dropped": ("closed form without the digamma term",
                        [("d1 = (arb(q) / arb.pi()).log() / 2 + x.digamma() / 2 + lpl",
                          "d1 = (arb(q) / arb.pi()).log() / 2 + lpl")]),
    "chi_sign_error": ("chi_4(3) = +1 and chi_3(2) = +1 in the Dirichlet series",
                       [("return {1: 1, 3: -1}.get(n % 4, 0)", "return {1: 1, 3: 1}.get(n % 4, 0)"),
                        ("return {1: 1, 2: -1}.get(n % 3, 0)", "return {1: 1, 2: 1}.get(n % 3, 0)")]),
    "series_tail_dropped": ("drop the Dirichlet-series tail radii",
                            [("lpl = -s1 + rad_ball(tail1)", "lpl = -s1"), ("dlpl = s2 + rad_ball(tail2)", "dlpl = s2")]),
    "cprime_wrong_exponent": ("c' = (4/g)(1+(M/g)^2)^{-1} instead of ^{-3/2}",
                              [("cp = (4 / g) / ((1 + (M / g) ** 2) ** 3).sqrt()",
                                "cp = (4 / g) / ((1 + (M / g) ** 2) ** 2).sqrt()")]),
    "c_large_gamma_approx": ("c = 4M/g instead of 4M/sqrt(M^2+g^2) in the head",
                             [("c = 4 / (1 + (g / M) ** 2).sqrt()", "c = 4 * M / g")]),
    "watson_too_strong": ("large-t remainder with abs(J0(x)) <= 1/x instead of x^{-1/2}",
                          [("prod *= 1 / x.sqrt()            # (c_g(Ma) T1)^{-1/2} <= 1", "prod *= 1 / x")]),
    "z0_2p3": ("tail-model range z0 = 2.3 (< j_{0,1}: still a valid certificate)",
               [("Z0 = arb(1)            # tail-model range", 'Z0 = arb("2.3")   # tail-model range')]),
    "z0_2p5": ("tail-model range z0 = 2.5 (> j_{0,1} = 2.4048: INVALID)",
               [("Z0 = arb(1)            # tail-model range", 'Z0 = arb("2.5")   # tail-model range')]),
    "z0_6": ("tail-model range z0 = 6 (INVALID, but J0(6) > 0)",
             [("Z0 = arb(1)            # tail-model range", "Z0 = arb(6)        # tail-model range")]),
    "z0_range_check_removed_T1x4": ("remove the check c_{g_N} T1 <= z0 and integrate the tail model up to 4*T1max",
                                    [("if not (self.cG * T).upper() <= Z0.lower():", "if False:"),
                                     ("tmax = float(self.T1max.lower())", "tmax = 4 * float(self.T1max.lower())"),
                                     ("E2_TARGET = 1e-10 ", "E2_TARGET = 1e-300 ")]),
    "pass_test_inverted": ("PASS test inverted (hi > 0)",
                           [("ok = bool(enc.upper() < 0) and hi < 0.0", "ok = bool(enc.upper() > 0) and hi > 0.0")]),
}


def _z(q):
    with open(os.path.join(ZDIR, f"chi{q}_zeros.txt")) as fh:
        return [x.strip() for x in fh if x.strip()]


def _shift(s, d):
    from decimal import Decimal
    return str(Decimal(s) + Decimal(d))


# name -> (description, function(list_of_strings) -> list_of_strings)   zero-file mutations
ZERO = {
    "zero1_plus_1e-3": ("perturb gamma_1 by +1e-3", lambda z: [_shift(z[0], "0.001")] + z[1:]),
    "zero1_plus_0.05": ("perturb gamma_1 by +0.05", lambda z: [_shift(z[0], "0.05")] + z[1:]),
    "zero1_plus_0.5": ("perturb gamma_1 by +0.5", lambda z: [_shift(z[0], "0.5")] + z[1:]),
    "zero250_plus_0.05": ("perturb gamma_250 by +0.05", lambda z: z[:249] + [_shift(z[249], "0.05")] + z[250:]),
    "zero10_out_of_order": ("gamma_10 := gamma_11 + 0.01 (order broken)",
                            lambda z: z[:9] + [_shift(z[10], "0.01")] + z[10:]),
    "zero5_duplicated": ("gamma_5 listed twice", lambda z: z[:5] + [z[4]] + z[5:]),
    "drop_zero1": ("drop gamma_1 from the list", lambda z: z[1:]),
    "drop_zero3": ("drop gamma_3 from the list", lambda z: z[:2] + z[3:]),
    "drop_zero100": ("drop gamma_100 from the list", lambda z: z[:99] + z[100:]),
    "drop_zero400": ("drop gamma_400 from the list", lambda z: z[:399] + z[400:]),
    "drop_last_zero": ("drop the LAST listed zero (still a complete list: valid)", lambda z: z[:-1]),
    "drop_last_100": ("drop the last 100 listed zeros (still complete: valid)", lambda z: z[:-100]),
    "spurious_zero": ("insert a spurious zero at gamma_1 + 1", lambda z: [z[0], _shift(z[0], "1")] + z[1:]),
    "file_truncated_50": ("zero file truncated to 50 lines", lambda z: z[:50]),
}


def build(name: str):
    with open(SRC) as fh:
        src = fh.read()
    zdir = ZDIR
    if name in CODE:
        for old, new in CODE[name][1]:
            if src.count(old) != 1:
                raise RuntimeError(f"pattern for {name!r} occurs {src.count(old)} times: {old!r}")
            src = src.replace(old, new)
    elif name in ZERO:
        zdir = tempfile.mkdtemp(prefix="zeros_" + name + "_", dir=os.environ.get("MUT_TMP", None))
        for q in (3, 4):
            with open(os.path.join(zdir, f"chi{q}_zeros.txt"), "w") as fh:
                fh.write("\n".join(ZERO[name][1](_z(q))) + "\n")
    else:
        raise KeyError(name)
    mod = types.ModuleType("mutant_" + name)
    mod.__file__ = SRC
    exec(compile(src, "mutant:" + name, "exec"), mod.__dict__)
    mod.ZERO_DIR = zdir
    return mod


_TRUTH = {}


def truth(q, Ma, Mb):
    """Independent rigorous point enclosures of delta' at Ma, mid, Mb."""
    key = (q, Ma, Mb)
    if key not in _TRUTH:
        import indep_enclose as ie
        a, b = Fraction(Ma), Fraction(Mb)
        pts = [a, (a + b) / 2, b]
        _TRUTH[key] = [(float(p),) + (lambda r: (r["lo"], r["hi"]))(ie.enclose_cell(q, p, p)) for p in pts]
    return _TRUTH[key]


def run_one(name, q, Ma, Mb, tol=0.02):
    t0 = time.time()
    out = dict(name=name, q=q, Ma=Ma, Mb=Mb)
    try:
        mod = build(name)
        r = mod.certify_cell(q, Ma, Mb, None, tol, 6000)
        out.update(status="PASS" if r["ok"] else "FAIL", lo=r["lo"], hi=r["hi"], T1=r["T1"], K=r["K"],
                   N=r["N"], ncells=r["ncells"])
    except Exception as e:      # the driver's fail-closed path (_work) does the same
        out.update(status="RAISES", err=f"{type(e).__name__}: {e}")
    out["secs"] = round(time.time() - t0, 1)
    if out["status"] in ("PASS", "FAIL"):
        tr = truth(q, Ma, Mb)
        out["truth"] = tr
        out["contains_truth"] = all(out["lo"] <= tlo and thi <= out["hi"] for _, tlo, thi in tr)
    return out


def verdict(o):
    if o["status"] == "RAISES":
        return "RAISES (fail-closed)"
    if o["status"] == "FAIL":
        return "FAIL (fail-closed)"
    return "PASS, CONTAINS truth" if o["contains_truth"] else "PASS, EXCLUDES truth"


def cmd_all(outpath):
    names = list(CODE) + list(ZERO)
    res = []
    for name in names:
        for q, Ma, Mb in CELLS[:2] if name in CODE else CELLS:
            try:
                p = subprocess.run([sys.executable, "-B", os.path.abspath(__file__), "one", name, str(q), Ma, Mb],
                                   capture_output=True, text=True, timeout=900)
                o = json.loads(p.stdout.strip().splitlines()[-1])
            except subprocess.TimeoutExpired:
                o = dict(name=name, q=q, Ma=Ma, Mb=Mb, status="RAISES", err="timeout 900 s", secs=900)
            except Exception as e:
                o = dict(name=name, q=q, Ma=Ma, Mb=Mb, status="RAISES", err=f"harness: {e}", secs=0)
            o["verdict"] = verdict(o)
            desc = (CODE.get(name) or ZERO.get(name))[0]
            o["desc"] = desc
            res.append(o)
            enc = f"[{o['lo']:.6e}, {o['hi']:.6e}]" if "lo" in o else o.get("err", "")
            print(f"{name:32s} q={q} [{Ma},{Mb}]  {o['verdict']:24s} {enc}  ({o['secs']} s)", flush=True)
    with open(outpath, "w") as fh:
        json.dump(res, fh, indent=1)


# ------------------------------------------------------------------------------------
def cmd_art():
    """Artificial amplitude family with KNOWN POSITIVE derivative, fed to the
    unmutated checker's machinery (Cell.choose_T1 / Cell.integrate / the PASS test).

    Family: c_g(u) = (M0/u) c_g(M0), u in [ua, ub] around M0 (a pure dilation).
    Then X_u = 1 + (M0/u) S, S = sum c_g(M0) cos(theta_g), and
        delta~(u) = P(S > -u/M0),   d/du delta~ = f_S(-u/M0)/M0 > 0
    (f_S = density of S, strictly positive).  All fields the checker uses are
    filled in for this family with the same rigour as in Cell.__init__."""
    mod = build("baseline")
    arb, ctx = mod.arb, mod.ctx
    out = []
    for q, M0s, uas, ubs in ((4, "5", "5", "5.02"), (4, "17.4", "17.35", "17.4"), (3, "8", "8", "8.03")):
        mod.setup()
        M0, ua, ub = mod.fr(M0s), mod.fr(uas), mod.fr(ubs)
        zeros = mod.load_zeros(q)
        ctx.prec = mod.PREC_HI
        VT0, DT0, s20, D0 = mod.tail_VD(q, zeros, mod.ex(M0))
        c0 = [4 / (1 + (g / mod.ex(M0)) ** 2).sqrt() for g in zeros]
        u = mod.hull(mod.ex(ua), mod.ex(ub))
        lam = mod.ex(M0) / u                       # lambda(u) = M0/u, decreasing
        lamp = -mod.ex(M0) / (u * u)               # lambda'(u) < 0
        lam_a, lam_b = mod.ex(M0) / mod.ex(ua), mod.ex(M0) / mod.ex(ub)
        cell = object.__new__(mod.Cell)
        cell.q, cell.Ma, cell.Mb, cell.N, cell.zeros = q, ua, ub, len(zeros), zeros
        cell.cs = [lam * c for c in c0]
        cell.cps = [lamp * c for c in c0]
        cell.V = lam * lam * VT0                    # sum_tail c^2/2
        cell.DT = 2 * lam * lamp * VT0              # sum_tail c c'  (negative)
        cell.sig2 = lam * lam * s20
        cell.info = dict(VT=(mod.dn(cell.V), mod.up(cell.V)), DT=(mod.dn(cell.DT), mod.up(cell.DT)),
                         sig2=(mod.dn(cell.sig2), mod.up(cell.sig2)))
        cell.cs_a = [lam_b * c for c in c0]         # lower bounds of c over the cell
        cell.c1_b = lam_a * c0[0]                   # upper bound of the largest c
        cell.cG = arb((lam_a * c0[-1]).upper())
        # sum |c c'| over ALL zeros = |lam lam'| * 2 sigma^2(M0)
        cell.Dmax = arb((abs(lam * lamp) * 2 * s20).upper())
        cell.D_info = (0.0, 0.0, mod.up(cell.Dmax))
        j0z, j1z = mod.Z0.bessel_j(0), mod.Z0.bessel_j(1)
        C4 = (-(j0z.log()) - mod.Z0 ** 2 / 4) / mod.Z0 ** 4
        kap = (j1z / j0z - mod.Z0 / 2) / mod.Z0 ** 3
        cell.C4S4 = arb((C4 * cell.cG ** 2 * 2 * arb(cell.V.upper())).upper())
        cell.eta2 = arb((2 * kap * cell.cG ** 2).upper())
        cell.T1max = mod.Z0 / cell.cG
        ctx.prec = mod.PREC_LO
        T1, E2, K = cell.choose_T1()
        I1, model, ncells = cell.integrate(T1, 0.02, 6000)
        enc = (I1 + mod.rad_ball(E2)) / arb.pi()
        lo, hi = mod.dn(enc), mod.up(enc)
        ok = bool(enc.upper() < 0) and hi < 0.0     # the checker's PASS test, verbatim
        # non-rigorous value of f_S(-u/M0)/M0 at u = ua (numpy/scipy)
        import numpy as np
        from scipy.special import j0 as J0
        g = np.array([float(z.mid()) for z in zeros])
        m0 = float(M0)
        c = 4 * m0 / np.sqrt(m0 * m0 + g * g)
        vt = float(VT0.mid())
        tt = np.linspace(0, 40 / np.sqrt(float(s20.mid())), 40001)
        Phi = np.prod(J0(np.outer(tt, c)), axis=1) * np.exp(-vt * tt * tt / 2)
        x = -float(ua) / m0
        fS = np.trapezoid(np.cos(tt * x) * Phi, tt) / np.pi
        ref = fS / m0
        line = (f"artificial dilation family q={q} M0={M0s} u in [{uas},{ubs}]: checker enclosure "
                f"[{lo:.6e}, {hi:.6e}]  PASS={ok}  (float reference at u_a: {ref:.6e}; truth > 0)")
        print(line, flush=True)
        out.append(dict(q=q, M0=M0s, ua=uas, ub=ubs, lo=float(lo), hi=float(hi), ok=bool(ok), ref=float(ref),
                        fail_closed=bool((not ok) and lo > 0 and lo <= ref <= hi)))
    return out


def cmd_driver():
    """Driver-level behaviour: run_range/summarize on a failing mutant, and
    summarize() on hand-made result files."""
    out = {}
    tmp = tempfile.mkdtemp(prefix="driver_", dir=os.environ.get("MUT_TMP", None))
    # the range driver uses multiprocessing, so the mutant must be a real file
    with open(SRC) as fh:
        src = fh.read()
    for old_, new_ in CODE["sinc_sign_flip"][1] + [
            ("max_cells, max_depth=6):", "max_cells, max_depth=1):"),
            ('ZERO_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))',
             "ZERO_DIR = %r" % ZDIR)]:
        assert src.count(old_) == 1, old_
        src = src.replace(old_, new_)
    mfile = os.path.join(tmp, "mutant_sinc.py")
    with open(mfile, "w") as fh:
        fh.write(src)
    f = os.path.join(tmp, "r.txt")
    p = subprocess.run([sys.executable, "-B", mfile, "range", "--q", "4", "--lo", "17.3701", "--hi", "17.4",
                        "--ratio", "1.005", "--tol", "0.02", "--workers", "1", "--out", f],
                       capture_output=True, text=True, timeout=1200)
    txt = open(f).read()
    print(txt)
    base = build("baseline")
    # (i) gap in coverage
    f = os.path.join(tmp, "gap.txt")
    open(f, "w").write("3.6 3.7 hi=-1.0e-2 lo=-2.0e-2 PASS\n3.8 3.9 hi=-1.0e-2 lo=-2.0e-2 PASS\n")
    base.summarize(f, base.fr("3.6"), base.fr("3.9"))
    out["gap_detected"] = "fully certified: False" in open(f).read()
    # (ii) a FAIL line
    f = os.path.join(tmp, "fail.txt")
    open(f, "w").write("3.6 3.7 hi=-1.0e-2 lo=-2.0e-2 PASS\n3.7 3.9 hi=nan lo=nan FAIL\n")
    base.summarize(f, base.fr("3.6"), base.fr("3.9"))
    out["fail_line_detected"] = "fully certified: False" in open(f).read()
    # (iii) a line marked PASS whose hi is POSITIVE (tampered / inconsistent file)
    f = os.path.join(tmp, "pos.txt")
    open(f, "w").write("3.6 3.9 hi=+1.0e-2 lo=-2.0e-2 PASS\n")
    base.summarize(f, base.fr("3.6"), base.fr("3.9"))
    out["positive_hi_marked_PASS_rejected"] = "fully certified: False" in open(f).read()
    # (iv) overlapping cells
    f = os.path.join(tmp, "ovl.txt")
    open(f, "w").write("3.6 3.8 hi=-1.0e-2 lo=-2.0e-2 PASS\n3.7 3.9 hi=-1.0e-2 lo=-2.0e-2 PASS\n")
    base.summarize(f, base.fr("3.6"), base.fr("3.9"))
    out["overlap_accepted_as_cover"] = "fully certified: True" in open(f).read()
    # (v) float inputs are rejected
    try:
        base.certify_cell(4, 3.6, 3.618)
        out["float_endpoints_rejected"] = False
    except TypeError:
        out["float_endpoints_rejected"] = True
    print("driver/summarize checks:", out)
    return out


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "list":
        for k, v in list(CODE.items()) + list(ZERO.items()):
            print(f"{k:32s} {v[0]}")
    elif cmd == "one":
        o = run_one(sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
        print(json.dumps(o))
    elif cmd == "all":
        cmd_all(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "mutation_results.json"))
    elif cmd == "art":
        r = cmd_art()
        json.dump(r, open(os.path.join(HERE, "artificial_family.json"), "w"), indent=1)
    elif cmd == "driver":
        r = cmd_driver()
        json.dump(r, open(os.path.join(HERE, "driver_checks.json"), "w"), indent=1)
