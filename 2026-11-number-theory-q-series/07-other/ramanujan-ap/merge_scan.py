#!/usr/bin/env python3
"""Merge apscan outputs.  Usage: merge_scan.py out1 [out2 ...]  > summary
Prints, per class: number of failing primes, first/last failing prime with (n, m),
and the per-bin table (bin = floor(100*log10 p)): primes, failing primes, max rho."""
import sys, gzip, collections
CQ = [1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12]
CA = [0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11]
bins = collections.defaultdict(lambda: [0, 0, 0, 0.0, 0, 0, 0, 0])   # (c,bin) -> nprimes,nfail,nfail0,maxrho,p,n,m,nstart
first = {}; last = {}; tot = collections.Counter(); nstart = collections.Counter(); amb = []; ends = {}
tasks = set()
for fn in sys.argv[1:]:
    op = gzip.open if fn.endswith('.gz') else open
    with op(fn, 'rt') as f:
        for line in f:
            w = line.split()
            if not w or w[0] == '#': continue
            if w[0] == 'B':
                a, b, c, bn = int(w[1]), int(w[2]), int(w[3]), int(w[4])
                tasks.add((a, b))
                r = bins[(c, bn)]
                r[0] += int(w[5]); r[1] += int(w[6]); r[2] += int(w[7])
                rho = float(w[8])
                if rho > r[3]: r[3:7] = [rho, int(w[9]), int(w[10]), int(w[11])]
                if len(w) > 12: r[7] += int(w[12]); nstart[c] += int(w[12])
            elif w[0] == 'F':
                c = int(w[3]); fp = tuple(map(int, w[4:7])); lp = tuple(map(int, w[7:10]))
                if c not in first or fp[0] < first[c][0]: first[c] = fp
                if c not in last or lp[0] > last[c][0]: last[c] = lp
                tot[c] += int(w[10])
            elif w[0] == 'A': amb.append(line.strip())
            elif w[0] == 'E':
                ends[(int(w[2]), int(w[3]))] = int(w[4])
# coverage
iv = sorted(tasks); cov = []
for a, b in iv:
    if cov and cov[-1][1] == a: cov[-1][1] = b
    else: cov.append([a, b])
print("# coverage (exact scan, all primes >= 7):", cov)
print("# ambiguous comparisons:", len(amb))
for c in range(17):
    if c not in tot and not any(k[0] == c for k in bins): continue
    print(f"CLASS {c} q={CQ[c]} a={CA[c]}  failing_primes={tot[c]}  components_started={nstart[c]}  first={first.get(c)}  last={last.get(c)}")
for c in range(17):
    for (cc, bn) in sorted(bins):
        if cc != c: continue
        r = bins[(cc, bn)]
        print(f"BIN {c} {bn} {r[0]} {r[1]} {r[2]} {r[3]:.15g} {r[4]} {r[5]} {r[6]} {r[7]}")
