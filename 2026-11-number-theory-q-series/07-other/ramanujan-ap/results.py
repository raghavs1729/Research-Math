#!/usr/bin/env python3
"""Final table: for each class, the last failing prime found by the exact scan, the exact
end of the last component of the failure set (Arb, via last_interval.py), and the
statistic W(p*) + c B(L*) (the fluctuation level reached at the last failure, model.py).
Usage: results.py merged.txt [more merged files, e.g. from flagged-row scans]"""
import sys, re, subprocess, math
from model import W_of, B_of, CSQ
CQ = [1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12]
CA = [0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11]
last = {}; nfail = {}; ncomp = {}
for fn in sys.argv[1:]:
    for l in open(fn):
        if l.startswith("CLASS"):
            c = int(l.split()[1])
            m = re.search(r"failing_primes=(\d+)\s+components_started=(\d+).*last=\((\d+), (\d+), (\d+)\)", l)
            if not m: continue
            nfail[c] = nfail.get(c, 0) + int(m.group(1)); ncomp[c] = ncomp.get(c, 0) + int(m.group(2))
            t = (int(m.group(3)), int(m.group(4)), int(m.group(5)))
            if c not in last or t[0] > last[c][0]: last[c] = t
        elif l.startswith("# coverage"): print(l.strip())
print("| (q,a) | c | failing primes ≥ 7 | components starting at primes ≥ 7 | last failing prime p* | π(p*;q,a) | π(p*/e;q,a) | sup F (end of last component) | largest integer counterexample | W(p*) + c·B | start of last component |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for c in range(17):
    if c not in last: continue
    p, n, m = last[c]
    out = subprocess.run([sys.executable, "last_interval.py", str(CQ[c]), str(CA[c]), str(p), str(n), str(m)], capture_output=True, text=True)
    if out.returncode != 0:
        print("ERROR for class", c, out.stdout[-300:], out.stderr[-500:]); continue
    end = re.search(r"end of failure set: (.*)", out.stdout).group(1)
    N = int(re.search(r"largest integer counterexample: (\d+)", out.stdout).group(1))
    st = int(re.search(r"last component starts at the prime (\d+)", out.stdout).group(1))
    if end.startswith("root"):
        x = re.search(r"\[\[([0-9.]+)", end).group(1)
        endtxt = f"{x[:len(str(N)) + 7]}… (zero of G, attained)"
    else:
        l = re.search(r"e\*(\d+) = \[?([0-9.]+)", end)
        endtxt = f"e·{l.group(1)} = {l.group(2)[:len(str(N)) + 7]}… (not attained)"
    cs = CSQ[(CQ[c], CA[c])]
    lvl = W_of(p) + cs * B_of(p)
    sp = lambda v: f"{v:,}".replace(",", " ")
    print(f"| ({CQ[c]},{CA[c]}) | {cs} | {sp(nfail[c])} | {sp(ncomp[c])} | {sp(p)} | {sp(n)} | {sp(m)} | {endtxt} | {sp(N)} | {lvl:.3f} | {sp(st)} |")
