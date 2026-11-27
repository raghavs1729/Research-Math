#!/usr/bin/env python3
"""AUDITOR: structural verification of the audited result files and comparison
with the auditor's from-scratch replay.

For each modulus:
 1. every data line of ../cert/results_q{q}.txt is parsed with exact rationals
    (own parser); it must end in PASS, have lo <= hi < 0 (as exact decimals),
    and M_a < M_b;
 2. the cells must tile the claimed range EXACTLY (no gap, no overlap);
 3. if replay/results_q{q}.txt exists: it must contain exactly the same cells
    with exactly the same printed hi/lo/E2/Emodel/rad/T1/K/N/tcells/VT/DT
    fields (everything except the timing field `secs` and the line order).

Exit status 0 iff all checks pass.
"""
import os
import sys
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "cert")
RANGES = {4: ("3.6", "17.4"), 3: ("4.4", "19.2")}
bad = []


def parse(path):
    cells = {}
    extra = []
    with open(path, encoding="ascii") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith("#"):
                if "not certified" in ln or "FAIL" in ln:
                    extra.append(ln)
                continue
            p = ln.split()
            a, b = Fraction(p[0]), Fraction(p[1])
            kv = dict(x.split("=", 1) for x in p[2:-1])
            key = (a, b)
            if key in cells:
                extra.append("duplicate cell " + ln)
            cells[key] = dict(status=p[-1], hi=kv.get("hi"), lo=kv.get("lo"),
                              fields={k: v for k, v in kv.items() if k != "secs"}, line=ln)
    return cells, extra


def check_file(q, path, label):
    lo, hi = (Fraction(x) for x in RANGES[q])
    cells, extra = parse(path)
    ok = True
    for e in extra:
        print(f"  [{label}] unexpected line: {e}")
        ok = False
    worst = None
    maxw = 0
    for (a, b), c in cells.items():
        if c["status"] != "PASS":
            print(f"  [{label}] not PASS: {c['line']}")
            ok = False
            continue
        try:
            h, l = Fraction(Decimal(c["hi"])), Fraction(Decimal(c["lo"]))
        except Exception:
            print(f"  [{label}] unparsable: {c['line']}")
            ok = False
            continue
        if not (a < b and l <= h < 0):
            print(f"  [{label}] bad enclosure: {c['line']}")
            ok = False
        if worst is None or h > worst[0]:
            worst = (h, a, b)
        maxw = max(maxw, float((h - l) / -l))
    pos = lo
    for (a, b) in sorted(cells):
        if a != pos:
            print(f"  [{label}] tiling broken at {float(pos)} (next cell starts at {float(a)})")
            ok = False
            break
        pos = b
    if pos != hi:
        print(f"  [{label}] cells end at {float(pos)}, not at {float(hi)}")
        ok = False
    print(f"  [{label}] q={q}: {len(cells)} cells, all PASS with hi < 0: {ok}; exact tiling of "
          f"[{RANGES[q][0]}, {RANGES[q][1]}]: {pos == hi}; worst hi = {float(worst[0]):.6e} on "
          f"[{float(worst[1])}, {float(worst[2])}]; max (hi-lo)/|lo| = {maxw:.4f}")
    return cells, ok


for q in (4, 3):
    print(f"q = {q}")
    aud, ok = check_file(q, os.path.join(CERT, f"results_q{q}.txt"), "audited")
    if not ok:
        bad.append(f"audited q={q}")
    rp = os.path.join(os.environ.get("AUDIT_REPLAY", os.path.join(HERE, "replay")), f"results_q{q}.txt")
    if os.path.exists(rp):
        rep, ok = check_file(q, rp, "replay ")
        if not ok:
            bad.append(f"replay q={q}")
        same = set(aud) == set(rep)
        ndiff = 0
        for k in set(aud) & set(rep):
            if aud[k]["fields"] != rep[k]["fields"]:
                ndiff += 1
                if ndiff <= 5:
                    print("  DIFFERENT:\n    audited:", aud[k]["line"], "\n    replay :", rep[k]["line"])
        print(f"  replay vs audited: same cells: {same}; cells with any differing field (except secs): {ndiff}")
        if not same or ndiff:
            bad.append(f"replay mismatch q={q}")
    else:
        print("  (no replay file yet)")
        if "--require-replay" in sys.argv:
            bad.append(f"no replay q={q}")

if bad:
    print("VERIFY FAILURES:", bad)
    sys.exit(1)
print("results files verified")
