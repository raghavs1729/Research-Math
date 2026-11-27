# AUDIT-3: the four pairs with 2f+1 <= d < 6, and spot checks of source table rows, from my own simulation (sf3 seq).
# (k, p, m) in the sense of the source: minimal period m of the characteristic sequence, k / p = number of entries
# in the minimal preperiod / in one period.  Periodicity here is "observed on the computed range" (the range is stated).
import subprocess, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def seq(f, g, h, X):
    out = subprocess.run([os.path.join(HERE, 'sf3'), 'seq', str(f), str(g), str(h), str(X)], capture_output=True, text=True, check=True).stdout
    return [int(t) for t in out.split()]
def kpm(S, X):
    c = np.zeros(X+1, dtype=np.int8); c[S] = 1
    tail0 = X // 3
    for m in range(1, X // 4):
        if np.array_equal(c[tail0:X+1-m], c[tail0+m:X+1]): break
    else: return None
    N0 = tail0
    while N0 >= 1 and c[N0] == c[N0+m]: N0 -= 1
    return int(c[1:N0+1].sum()), int(c[N0+1:N0+1+m].sum()), m, N0
def T(f, d, X):
    m = 8*d-2*f-1
    t = [f, d-f] + list(range(2*d-f, 3*d-f)) + [4*d-f]; j = 0
    while 10*d-3*f-2 + j*m <= X: t += list(range(10*d-3*f-2+j*m, 11*d-3*f-2+j*m)); j += 1
    return [x for x in t if x <= X]
print("pairs with 2f+1 <= d < 6:")
for f, d, X in [(1,3,3000), (1,4,6000), (1,5,20000), (2,5,400000)]:
    S = seq(f, d-f, 2*d-f, X); Tt = T(f, d, X)
    i = next(i for i in range(min(len(S), len(Tt))) if S[i] != Tt[i])
    r = kpm(S, X)
    print(f"  (f,d)=({f},{d}) S_{{{f},{d-f},{2*d-f}}} = {S[:16]}...; T = {Tt[:12]}...; first deviation at entry #{i+1}: S has {S[i]}, T has {Tt[i]};")
    print(f"     observed (k,p,m) = {r[:3]} on [1,{X}] (characteristic sequence periodic from {r[3]+1}); Theorem-2 pattern would be {(d+3, d, 8*d-2*f-1)}")
print("spot checks against the source's tables (rows with f+g = d and a few neighbours):")
rows = [((2,4,10),(9,6,43)), ((2,8,18),(13,10,75)), ((3,7,17),(13,10,73)), ((4,6,16),(13,10,71)), ((2,5,12),(10,7,51)),
        ((12,13,38),(28,25,175)), ((1,2,6),(6,6,28)), ((1,3,8),(8,6,38)), ((2,3,9),(11,7,44)), ((1,4,10),(9,6,48))]
for (f, g, h), src in rows:
    X = 60*src[2] + 2000
    r = kpm(seq(f, g, h, X), X)
    print(f"  ({f},{g},{h}): mine {r[:3]}, source {src}: {'AGREE' if r[:3] == src else 'DIFFER'}")
