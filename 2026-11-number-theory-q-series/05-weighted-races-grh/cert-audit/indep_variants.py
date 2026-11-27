#!/usr/bin/env python3
"""AUDITOR: robustness of indep_enclose.py -- the same cells with different
internal parameters (series order NMAX, expansion radius z0, remainder goal,
domain margin, precision).  All enclosures of the same cell must have a common
point; reported: the intersection and the spread."""
import sys
from fractions import Fraction
import indep_enclose as ie

cells = [(4, "3.6", "3.618"), (4, "8.4754", "8.5177"), (4, "17.3701", "17.4"),
         (3, "4.4", "4.422"), (3, "19.1282", "19.2"), (4, "5", "5"), (3, "12", "12")]
variants = [("default (NMAX=6, z0=3/2, margin 1.3, E2 goal 1e-12, 110 bits)", {}),
            ("NMAX=3", dict(NMAX=3)), ("NMAX=9", dict(NMAX=9)),
            ("z0=1", dict(z0=Fraction(1))), ("z0=2.2", dict(z0=Fraction(11, 5))),
            ("E2 goal 1e-7", dict(e2=1e-7)), ("E2 goal 1e-20", dict(e2=1e-20)),
            ("margin 2.0", dict(MARGIN=Fraction(2))), ("precision 200 bits", dict(PREC=200))]
bad = 0
for q, a, b in cells:
    los, his = [], []
    print(f"q={q} [{a}, {b}]")
    for name, kw in variants:
        old = (ie.NMAX, ie.MARGIN, ie.PREC, ie.LAM)
        if "NMAX" in kw:
            ie.NMAX = kw["NMAX"]
            ie.LAM = ie.rayleigh(ie.NMAX + 1)
        if "MARGIN" in kw:
            ie.MARGIN = kw["MARGIN"]
        if "PREC" in kw:
            ie.PREC = kw["PREC"]
        ckw = {}
        if "z0" in kw:
            ckw["z0"] = kw["z0"]
        if "e2" in kw:
            ckw["e2_goal"] = kw["e2"]
        r = ie.enclose_cell(q, a, b, **ckw)
        ie.NMAX, ie.MARGIN, ie.PREC, ie.LAM = old
        los.append(r["lo"]); his.append(r["hi"])
        print(f"   {name:62s} [{r['lo']:.10e}, {r['hi']:.10e}]  T1={r['T1']:.3f} K={r['K']} N'={r['Np']} evals={r['nev']}")
    ok = max(los) <= min(his)
    bad += not ok
    print(f"   common intersection [{max(los):.10e}, {min(his):.10e}]  {'CONSISTENT' if ok else 'INCONSISTENT'}")
print("all variants mutually consistent:", bad == 0)
sys.exit(1 if bad else 0)
