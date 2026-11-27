#!/usr/bin/env python3
"""AUDITOR: enclosure-property tests of the audited code and of its mutants,
at two levels below the final output.

Output-level mutation tests (mutate.py) cannot see a mutation that makes the
computation non-rigorous but leaves it numerically accurate, because the final
enclosure is 4 % wide.  These tests check the defining property of an
enclosure directly, against independently computed true values:

 L1  bessel01(x-ball) must contain J_0(x), J_1(x) at sample points x of the
     ball (reference: Arb at 200 bits AND mpmath; random balls incl. balls
     touching 0, wide balls, balls with x_lo > 1 where the Watson cap acts).
 L2  Cell.integrand(t-box) must contain the true F(t, M) = sinc(t) dPhi_M/dM(t)
     at the corners/midpoints of the (t, M)-box (reference: the auditor's
     independent integrand indep_enclose.Cell.integrand at a thin (t, M), a
     rigorous point enclosure of relative width ~1e-6 or better).

For the unmutated code every sample must be contained (a single miss would be
a RIGOR BREAK).  For a mutant, misses show that the mutation is detectable.

Usage: proptest.py NAME [NAME ...]      (names as in mutate.py; 'baseline' = audited code)
"""
import json
import os
import random
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import mutate  # noqa: E402
import indep_enclose as ie  # noqa: E402
from flint import acb, arb, ctx  # noqa: E402
import mpmath  # noqa: E402

mpmath.mp.dps = 40


def level1(mod, nballs=1500, seed=7):
    rnd = random.Random(seed)
    mod.setup()
    miss0 = miss1 = n = 0
    worst = None
    for i in range(nballs):
        kind = i % 6
        x0 = rnd.choice([rnd.uniform(1e-4, 0.05), rnd.uniform(0.05, 3), rnd.uniform(3, 30)])
        if kind == 0:
            r = x0                      # ball [0, 2 x0]
        elif kind == 1:
            r = x0 * rnd.uniform(0.3, 0.99)
        elif kind == 2:
            r = rnd.uniform(0.2, 1.5)
            x0 = max(x0, r)
        elif kind == 3:
            r = 10 ** rnd.uniform(-6, -2)
        elif kind == 4:
            r = x0 * 10 ** rnd.uniform(-3, -1)
        else:
            x0 = rnd.uniform(1.2, 30)   # Watson cap active
            r = min(x0 - 1.01, rnd.uniform(0.01, 2.5))
        r = min(r, x0)
        X0, R = Fraction(x0), Fraction(r)
        ctx.prec = mod.PREC_LO
        xb = arb(float(X0), float(R))   # radius rounded up by Arb: ball contains [x0-r, x0+r]
        try:
            J0, J1 = mod.bessel01(xb)
        except Exception as e:           # baseline must never raise here
            miss0 += 1
            worst = f"raised {type(e).__name__}: {e} on x0={x0}, r={r}"
            continue
        ctx.prec = 200
        for f in (Fraction(-1), Fraction(1), Fraction(0), Fraction(rnd.uniform(-1, 1)), Fraction(rnd.uniform(-1, 1))):
            p = X0 + f * R
            if p < 0:
                continue
            pa = arb(p.numerator) / arb(p.denominator)
            t0, t1 = pa.bessel_j(0), pa.bessel_j(1)
            n += 1
            if not J0.overlaps(t0):
                miss0 += 1
                worst = worst or f"J0 miss: x0={x0}, r={r}, x={float(p)}"
            if not J1.overlaps(t1):
                miss1 += 1
                worst = worst or f"J1 miss: x0={x0}, r={r}, x={float(p)}"
            if n % 97 == 0:              # spot-check the reference itself against mpmath
                m0 = mpmath.besselj(0, mpmath.mpf(p.numerator) / p.denominator)
                assert abs(mpmath.mpf(t0.mid().str(45, radius=False)) - m0) < mpmath.mpf(10) ** -30
        ctx.prec = mod.PREC_LO
    return dict(samples=n, J0_misses=miss0, J1_misses=miss1, first=worst)


_REF = {}


def ref_F(q, M: Fraction, t: Fraction):
    """Rigorous point enclosure of F(t, M) from the independent implementation."""
    key = (q, M)
    if key not in _REF:
        _REF[key] = ie.Cell(q, M, M)
    cell = _REF[key]
    ctx.prec = ie.PREC
    v = cell.integrand(acb(ie.qarb(t)))
    if not v.is_finite():
        return None
    return v.real


def level2(mod, q, Ma, Mb, seed=11):
    rnd = random.Random(seed)
    mod.setup()
    a, b = Fraction(Ma), Fraction(Mb)
    cell = mod.Cell(q, mod.fr(Ma), mod.fr(Mb))
    T1, E2, K = cell.choose_T1()
    Tref = Fraction(float(ie.Cell(q, a, a).Rdom.mid())) * Fraction(9, 10)
    Tmax = min(Fraction(T1), Tref)
    boxes = []
    for width_div in (2048, 512, 128):
        h = Tmax / width_div
        for _ in range(25):
            k = rnd.randrange(0, width_div)
            boxes.append((k * h, (k + 1) * h))
    boxes.append((Fraction(0), Tmax / 2048))
    n = miss = 0
    first = None
    for ta, tb in boxes:
        ctx.prec = mod.PREC_LO
        try:
            F, _ = cell.integrand(mod.ex(ta), mod.ex(tb))
        except Exception as e:
            miss += 1
            first = first or f"raised {type(e).__name__}: {e}"
            continue
        for t in (ta, (ta + tb) / 2, tb):
            for M in (a, (a + b) / 2, b):
                if t == 0:
                    continue
                tv = ref_F(q, M, t)
                ctx.prec = mod.PREC_LO
                if tv is None:
                    continue
                n += 1
                if not F.overlaps(tv):
                    miss += 1
                    first = first or (f"t-box [{float(ta):.6f},{float(tb):.6f}], t={float(t):.6f}, M={float(M)}: "
                                      f"true F = {tv.str(12)}, box enclosure = {F.str(12)}")
    return dict(samples=n, misses=miss, first=first, T1=float(T1), tested_up_to=float(Tmax))


if __name__ == "__main__":
    names = sys.argv[1:] or ["baseline"]
    out = {}
    for name in names:
        try:
            mod = mutate.build(name)
        except Exception as e:
            print(name, "build failed:", e)
            continue
        res = {}
        try:
            res["L1"] = level1(mod)
        except Exception as e:
            res["L1"] = dict(error=f"{type(e).__name__}: {e}")
        for q, Ma, Mb in mutate.CELLS[:2]:
            try:
                res[f"L2 q={q} [{Ma},{Mb}]"] = level2(mod, q, Ma, Mb)
            except Exception as e:
                res[f"L2 q={q} [{Ma},{Mb}]"] = dict(error=f"{type(e).__name__}: {e}")
        out[name] = res
        print(name)
        for k, v in res.items():
            print("   ", k, json.dumps(v))
        sys.stdout.flush()
    path = os.path.join(HERE, "proptest_results.json")
    old = {}
    if os.path.exists(path):
        old = json.load(open(path))
    old.update(out)
    json.dump(old, open(path, "w"), indent=1)
