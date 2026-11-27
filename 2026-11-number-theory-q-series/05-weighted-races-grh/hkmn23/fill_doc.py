#!/usr/bin/env python3
"""fill_doc.py -- (re)generate the result-dependent passages of HKMN23.md from t2_summary.json
and raw/t2_zeros_*.jsonl (run t2_assemble.py first).  The passages sit between HTML comment
markers <!--X-START--> ... <!--X-END-->; everything else in HKMN23.md is hand-written."""
import json, glob, re, math
S = json.load(open("t2_summary.json"))
und = {int(q): v for q, v in S["undecided"].items()}
u60 = [(q, a) for q in sorted(und) if q <= 60 for a in und[q]]
u120 = [(q, a) for q in sorted(und) if q > 60 for a in und[q]]
cand14 = {(7,1),(7,2),(7,4),(9,1),(9,4),(9,7),(14,1),(14,9),(14,11),(16,1),(16,9),(18,1),(18,7),(18,13)}
extra60 = [p for p in u60 if p not in cand14]
assert set(u60) >= cand14
fmt = lambda L: ", ".join(f"({q},{a})" for q, a in L) if L else "none"
ZW, ZB = {}, {}
for f in glob.glob("raw/t2_zeros_*.jsonl"):
    for l in open(f):
        d = json.loads(l); k = (d["q"], d["a"])
        if k not in ZB or d["F_lower_bound"] > ZB[k]["F_lower_bound"]:
            ZB[k] = d
        if d["witness"] and (k not in ZW or d["F_lower_bound"] > ZW[k]["F_lower_bound"]):
            ZW[k] = d
n60 = 1084
# pairs (2m, a), m odd, excluded only through the reduction (2m,a) = (m, a mod m) and a zero-side witness at m
PW = {}
for f in glob.glob("raw/t2_witness_120_*.jsonl"):
    for l in open(f):
        d = json.loads(l); PW[(d["q"], d["a"])] = d
red = sorted((q, a) for (q, a), d in PW.items()
             if q % 4 == 2 and d["N"] > 0 and not d["witness"] and (q // 2, a % (q // 2)) in ZW)
red60 = [p for p in red if p[0] <= 60]; red120 = [p for p in red if p[0] > 60]
b6 = f"""*All moduli q ≤ 60* (1084 pairs besides the 18 members; `t2_summary.out`). {n60 - len(u60)} pairs are
excluded: every non-square class (explicit prime-side witness, or Theorem 4 — NV(q) is
certified in this range), the class 1 at every modulus outside the 18 of Theorem 5(b), and
every other square class for which a prime-side or zero-side witness is certified, where for
`q = 2m`, m odd, the witness may be the one of `(m, a mod m)` (the reduction (i) of §0:
`F_{{2m,a}} = F_{{m,a}}`; used for {fmt(red60)}). **{len(u60)} square classes are
undecided**: the 14 of part (a)""" + (f", and {fmt(extra60)}" if extra60 else " and no others") + """.

(c) *Moduli 61 ≤ q ≤ 120.* Witnesses: prime side; Theorem 4 for non-square classes without
prime-side witness; zero side where it was run (`t2_summary.out` lists the pairs); and, for
{fmt(red120)}, the reduction to `q/2` and the zero-side witnesses at `q/2`. """ + (f"""Every pair is
excluded except {len(u120)} square classes `a ≠ 1` for which the zero-side search was not run or not
finished: {fmt(u120)}. None of these is claimed to be hard.""" if u120 else "Every pair is excluded.") + """"""
rows = []
lrows = []
for k in sorted(ZW):
    d = ZW[k]
    if d.get("method") == "LLL":
        lrows.append(f"    ({k[0]},{k[1]})".ljust(14) + f"{d['threshold']:.4f}    {d['remainder_upper']:.4f}        {d['y0']:+.4e}     {float(d['F_Z_at_y0'].split('[')[1].split(' ')[0]):9.4f}   {d['F_lower_bound']:9.4f}        {float(d['F_lower_bound'])/float(d['L1']):.3f}       K = {d['K']}, nu = {d['nu']}, T = {d['T']:.0f}, {d['num_brackets']} cosines")
        continue
    if k[0] in (21, 40, 48, 60) or k[1] == 1:
        continue
    rows.append(f"    ({k[0]},{k[1]})".ljust(14) + f"{d['threshold']:.4f}    {d['remainder_upper']:.4f}     {d['y0']:>16.6f}   {float(d['F_Z_at_y0'].split('[')[1].split(' ')[0]):9.4f}   {d['F_lower_bound']:9.4f}     T = {d['T']:.0f}, {d['num_brackets']} zeros")
nz = [k for k in sorted(ZB) if k not in ZW and k[0] not in (21, 40, 48, 60) and k[1] != 1]
ll = ("    (q,a)        (2/3)N    remainder r   y_0 (approx.)    F_Z(y_0)    certified F(y_0) ≥   / L_1(q)    parameters\n" + "\n".join(lrows))
q6 = ("Certified zero-side witnesses for square classes `a ≠ 1` (all have `F(y_0) ≥` the last column `> (2/3)N`):\n\n"
      "    (q,a)        (2/3)N    remainder r   y_0                F_Z(y_0)    certified F(y_0) ≥\n" + "\n".join(rows) + "\n\n"
      + ("Searched without success (best certified lower bound for `sup F`): "
         + "; ".join(f"({q},{a}): {ZB[(q,a)]['F_lower_bound']:.4f} against {ZB[(q,a)]['threshold']:.4f} (|y| ≤ {ZB[(q,a)]['Ymax']:.0f})" for q, a in nz) + "."
         if nz else "Every class searched on the zero side outside the knife-edge moduli received a witness (some only by the lattice method below).")
      + "\nFor `q = 11, 28, 32` the class 1 was searched as well: the maximum found is at `y ≈ 0` and equals `Φ(q)` up to the remainder (2.9147, 3.8555, 4.8940 against Φ = 2.9848, 3.9407, 5.0116; remainders 0.0745, 0.0914, 0.1297). For `q = 19` a larger value than `Φ(19) = 4.694` occurs at `y_0 = 2826.3` (`F ≥ 5.553`).")
s7 = (f"{n60 - len(u60)} of the 1084 pairs with `3 ≤ q ≤ 60` outside the 18 members are excluded under GRH(q) alone; "
      f"{len(u60)} square classes are undecided (the 14 of line 11" + (f", and {fmt(extra60)}" if extra60 else ", and no others") + "). "
      f"For `61 ≤ q ≤ 120`: all pairs excluded" + (f" except {len(u120)} square classes `a ≠ 1` whose zero-side search was not run or not finished." if u120 else "."))
doc = open("HKMN23.md").read()
for tag, txt in (("6B", b6), ("6Q", q6), ("S6", s7), ("LL", ll)):
    pat = re.compile(f"<!--{tag}-START-->.*?<!--{tag}-END-->", re.S)
    assert pat.search(doc), tag
    doc = pat.sub(lambda m: f"<!--{tag}-START-->{txt}<!--{tag}-END-->", doc)
open("HKMN23.md", "w").write(doc)
print("undecided q<=60:", len(u60), "extra:", extra60, "| 61..120:", len(u120))
