#!/usr/bin/env python3
"""AUDITOR: turn the audit's raw outputs into the markdown tables of AUDIT.md."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def indep_table(q, nrows=14):
    rows = []
    with open(os.path.join(HERE, f"indep_q{q}.txt")) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.split()
            rows.append(p)
    rows.sort(key=lambda p: float(p[0]))
    n = len(rows)
    idx = sorted(set([round(i * (n - 1) / (nrows - 3)) for i in range(nrows - 2)] + [n - 3, n - 2, n - 1]))
    # cell with the worst RELATIVE audited width, and the cell whose audited hi is closest to 0 (= last)
    out = ["| cell [M_a, M_b] | audited lo | audited hi | independent lo | independent hi | meet | indep. < 0 | indep. inside audited | T1 / K / N′ (indep.) |",
           "|---|---|---|---|---|---|---|---|---|"]
    for i in idx:
        p = rows[i]
        out.append(f"| [{p[0]}, {p[1]}] | {p[4]} | {p[5]} | {float(p[2]):.6e} | {float(p[3]):.6e} | {p[6]} | {p[7]} | {p[8]} | {p[9]} / {p[10]} / {p[11]} |")
    return "\n".join(out), n


def plain_table(q):
    out = ["| M | plain mpmath δ′(M) | audited cell | audited [lo, hi] | in cell interval | auditor's rigorous point enclosure | in |",
           "|---|---|---|---|---|---|---|"]
    with open(os.path.join(HERE, f"compare_plain_q{q}.log")) as fh:
        for ln in fh:
            if not ln.startswith("  M="):
                continue
            a = ln.split()
            M = a[1]
            plain = a[2].split("=")[1]
            cell = a[3][4:]
            aud = a[4][7:]
            ins = a[5]
            rig = a[10]
            ins2 = a[11]
            out.append(f"| {M} | {plain} | {cell} | {aud} | {ins} | {rig} | {ins2} |")
    return "\n".join(out)


def mutation_table():
    res = json.load(open(os.path.join(HERE, "mutants", "mutation_results.json")))
    prop = json.load(open(os.path.join(HERE, "mutants", "proptest_results.json")))
    by = {}
    for o in res:
        by.setdefault(o["name"], []).append(o)
    out = ["| mutation | what is changed | outcome per test cell | enclosure-property test (L1 Bessel / L2 integrand misses) |",
           "|---|---|---|---|"]
    for name, lst in by.items():
        cells = []
        for o in lst:
            v = o["verdict"]
            if "lo" in o:
                v += f" [{o['lo']:.4e}, {o['hi']:.4e}]"
            elif "err" in o:
                v += f": {o['err'][:70]}"
            cells.append(f"q={o['q']} [{o['Ma']},{o['Mb']}]: {v}")
        pr = ""
        if name in prop:
            p = prop[name]
            l1 = p.get("L1", {})
            l2 = [v for k, v in p.items() if k.startswith("L2")]
            pr = (f"L1: {l1.get('J0_misses', '?')}+{l1.get('J1_misses', '?')} of {l1.get('samples', '?')}; "
                  + "L2: " + ", ".join(f"{v.get('misses', '?')}/{v.get('samples', '?')}" for v in l2))
        out.append(f"| `{name}` | {lst[0]['desc']} | " + "<br>".join(cells) + f" | {pr} |")
    return "\n".join(out)


if __name__ == "__main__":
    for q in (4, 3):
        t, n = indep_table(q)
        print(f"\n### independent vs audited, q={q} ({n} cells in total)\n")
        print(t)
    for q in (4, 3):
        print(f"\n### plain values, q={q}\n")
        print(plain_table(q))
    print("\n### mutations\n")
    print(mutation_table())
