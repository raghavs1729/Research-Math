#!/usr/bin/env python3
"""AUDITOR: independent enclosures vs the audited ones, plain values vs cells,
and certified-density differences vs the delta' enclosures.

Stages (all for one modulus, --q):
  indep   independent rigorous enclosure (indep_enclose.py) of delta' on EVERY
          audited cell (each cell split into --nsub sub-cells, union taken);
          writes indep_q{q}.txt.  Requirements checked: auditor's enclosure is
          negative and intersects the audited one.
  plain   every mpmath value of plain_q{q}.txt must lie in the audited [lo,hi]
          of each cell containing its M, and in the auditor's rigorous point
          enclosure at the same M.
  grid    certified densities of race_mono_grid_certified_results.txt:
          for consecutive grid points M1 < M2 inside the range,
            (a) mean value theorem: the certified difference quotient must meet
                the hull of the delta' enclosures of the cells meeting [M1,M2];
            (b) integral: sum_cells [lo,hi]*|cell cap [M1,M2]| must meet the
                certified interval for delta(M2)-delta(M1);
          both with the audited enclosures and with the auditor's.
"""
import argparse
import os
import re
import sys
import time
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.environ.get("AUDIT_OUT", HERE)      # where indep_q*.txt / plain_q*.txt are read and written
CERT = os.path.join(HERE, "..", "cert")
EXP = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))
RANGES = {4: ("3.6", "17.4"), 3: ("4.4", "19.2")}


def audited_cells(q):
    out = []
    with open(os.path.join(CERT, f"results_q{q}.txt")) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            kv = dict(x.split("=", 1) for x in p[2:-1])
            out.append((Fraction(p[0]), Fraction(p[1]), Fraction(Decimal(kv["lo"])), Fraction(Decimal(kv["hi"])),
                        p[-1], p[0], p[1]))
    out.sort()
    return out


def indep_cells(q):
    out = {}
    path = os.path.join(DIR, f"indep_q{q}.txt")
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            out[(Fraction(p[0]), Fraction(p[1]))] = (Fraction(Decimal(p[2])), Fraction(Decimal(p[3])))
    return out


def stage_indep(q, nsub, only=None):
    import indep_enclose as ie
    cells = audited_cells(q)
    path = os.path.join(DIR, f"indep_q{q}.txt")
    done = indep_cells(q)
    fh = open(path, "a")
    if not done:
        fh.write(f"# auditor's independent enclosures (indep_enclose.py, {nsub} sub-cells per cell), q={q}\n"
                 "# M_a M_b ind_lo ind_hi aud_lo aud_hi intersects ind_negative ind_inside_aud T1 K Nprime evals secs\n")
    nbad = 0
    t00 = time.time()
    for a, b, alo, ahi, st, sa, sb in cells:
        if (a, b) in done or (only and (sa, sb) not in only):
            continue
        t0 = time.time()
        lo = hi = None
        nev = 0
        info = None
        for k in range(nsub):
            x0 = a + (b - a) * k / nsub
            x1 = a + (b - a) * (k + 1) / nsub
            r = ie.enclose_cell(q, x0, x1)
            lo = r["lo"] if lo is None else min(lo, r["lo"])
            hi = r["hi"] if hi is None else max(hi, r["hi"])
            nev += r["nev"]
            info = r
        inter = not (Fraction(hi) < alo or Fraction(lo) > ahi)
        neg = hi < 0
        inside = alo <= Fraction(lo) and Fraction(hi) <= ahi
        if not (inter and neg):
            nbad += 1
        fh.write(f"{sa} {sb} {lo:.9e} {hi:.9e} {float(alo):.6e} {float(ahi):.6e} {inter} {neg} {inside} "
                 f"{info['T1']:.4f} {info['K']} {info['Np']} {nev} {time.time() - t0:.1f}\n")
        fh.flush()
    fh.close()
    print(f"q={q}: independent enclosures written to {path} ({time.time() - t00:.0f} s); failures: {nbad}")
    return nbad


def summary_indep(q):
    cells = audited_cells(q)
    ind = indep_cells(q)
    n = inter = neg = inside = 0
    worst = None
    maxrel = 0.0
    for a, b, alo, ahi, st, sa, sb in cells:
        if (a, b) not in ind:
            continue
        lo, hi = ind[(a, b)]
        n += 1
        inter += not (hi < alo or lo > ahi)
        neg += hi < 0
        inside += alo <= lo and hi <= ahi
        if worst is None or hi > worst[0]:
            worst = (hi, sa, sb)
        maxrel = max(maxrel, float((hi - lo) / -lo))
    print(f"q={q}: {n}/{len(cells)} audited cells independently enclosed; intersect: {inter}; "
          f"independent enclosure negative: {neg}; independent inside audited: {inside}; "
          f"worst independent upper bound {float(worst[0]):.6e} on [{worst[1]}, {worst[2]}]; "
          f"max relative width {maxrel:.4f}")
    return n == len(cells) and inter == n and neg == n


def stage_plain(q):
    import indep_enclose as ie
    cells = audited_cells(q)
    lo0, hi0 = (Fraction(x) for x in RANGES[q])
    ok = True
    rows = []
    with open(os.path.join(DIR, f"plain_q{q}.txt")) as fh:
        lines = [ln for ln in fh if ln.startswith("q=")]
    npts = len(lines)
    for k, ln in enumerate(lines):
        kv = dict(x.split("=", 1) for x in ln.split())
        M = lo0 + (hi0 - lo0) * k / (npts - 1)
        assert abs(float(M) - float(kv["M"])) < 1e-9, "point grid mismatch"
        v = Fraction(Decimal(kv["ddelta"]))
        cont = [(a, b, alo, ahi) for a, b, alo, ahi, *_ in cells if a <= M <= b]
        inside = bool(cont) and all(alo <= v <= ahi for _, _, alo, ahi in cont)
        r = ie.enclose_cell(q, M, M)
        ins2 = r["lo"] <= float(v) <= r["hi"]
        ok &= inside and ins2
        a, b, alo, ahi = cont[0]
        rows.append((float(M), float(v), float(a), float(b), float(alo), float(ahi), inside, r["lo"], r["hi"], ins2))
        print(f"  M={float(M):9.5f} plain={float(v):.10e} cell[{float(a)},{float(b)}] audited[{float(alo):.6e},{float(ahi):.6e}] "
              f"{'INSIDE' if inside else 'OUTSIDE'} ({len(cont)} cell(s)); rigorous point [{r['lo']:.10e},{r['hi']:.10e}] "
              f"{'INSIDE' if ins2 else 'OUTSIDE'}")
    print(f"q={q}: all {npts} plain values inside audited cell intervals and inside the auditor's point enclosures: {ok}")
    return ok


def read_grid(q):
    out = {}
    cur = None
    with open(os.path.join(EXP, "race_mono_grid_certified_results.txt")) as fh:
        for line in fh:
            mm = re.match(r"== modulus (\d)", line)
            if mm:
                cur = int(mm.group(1))
                continue
            mm = re.match(r"\s*m=\s*([-\d.]+)\s+delta in \[([\d.]+), ([\d.]+)\]", line)
            if mm and cur == q:
                out[Fraction(mm.group(1)) + Fraction(1, 2)] = (Fraction(Decimal(mm.group(2))), Fraction(Decimal(mm.group(3))))
    return out


def stage_grid(q):
    cells = audited_cells(q)
    ind = indep_cells(q)
    grid = read_grid(q)
    lo0, hi0 = (Fraction(x) for x in RANGES[q])
    Ms = sorted(M for M in grid if lo0 <= M <= hi0)
    ok = True
    for M1, M2 in zip(Ms, Ms[1:]):
        d_lo = grid[M2][0] - grid[M1][1]
        d_hi = grid[M2][1] - grid[M1][0]
        q_lo, q_hi = d_lo / (M2 - M1), d_hi / (M2 - M1)
        res = []
        for label, use_ind in (("audited", False), ("auditor", True)):
            hull_lo = hull_hi = None
            i_lo = i_hi = Fraction(0)
            cov = Fraction(0)
            for a, b, alo, ahi, *_ in cells:
                ov = min(b, M2) - max(a, M1)
                if ov <= 0:
                    continue
                if use_ind:
                    if (a, b) not in ind:
                        continue
                    alo, ahi = ind[(a, b)]
                hull_lo = alo if hull_lo is None else min(hull_lo, alo)
                hull_hi = ahi if hull_hi is None else max(hull_hi, ahi)
                i_lo += alo * ov
                i_hi += ahi * ov
                cov += ov
            full = cov == M2 - M1
            mvt = full and not (q_hi < hull_lo or q_lo > hull_hi)
            integ = full and not (i_hi < d_lo or i_lo > d_hi)
            ok &= mvt and integ
            res.append(f"{label}: int delta' in [{float(i_lo):.9f},{float(i_hi):.9f}] "
                       f"{'MVT ok' if mvt else 'MVT FAIL'}, {'integral ok' if integ else 'integral FAIL'}")
        print(f"  M in [{float(M1)}, {float(M2)}]: delta(M2)-delta(M1) in [{float(d_lo):.9f},{float(d_hi):.9f}]  " + "; ".join(res))
    print(f"q={q}: certified density differences consistent with the delta' enclosures: {ok}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["indep", "plain", "grid", "summary"])
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--nsub", type=int, default=4)
    a = ap.parse_args()
    if a.stage == "indep":
        bad = stage_indep(a.q, a.nsub)
        ok = summary_indep(a.q) and not bad
    elif a.stage == "summary":
        ok = summary_indep(a.q)
    elif a.stage == "plain":
        ok = stage_plain(a.q)
    else:
        ok = stage_grid(a.q)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
