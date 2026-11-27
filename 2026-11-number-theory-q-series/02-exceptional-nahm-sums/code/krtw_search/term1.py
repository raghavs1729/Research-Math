import sys
sys.path.insert(0, '/home/user/weighted-prime-races/monotonicity/rr-identities/code')
from fractions import Fraction as Fr
from classsearch import do
from moves import canon, terminals, defusions
from bfs1 import build, X
from bfs2 import fmtB
def Bdata(k):
    # vars n_1..n_k (d=1), n_{k+1} (d=2)
    def Nv(i): return [1 if j >= i else 0 for j in range(1, k+2)]
    vecs = [Nv(i) for i in range(1, k+2)]
    wt = [0]*(k+1)
    for v in vecs: wt = [w+x for w, x in zip(wt, v)]
    return build(vecs, wt, [1]*k+[2])
def Kdata(k):
    n = 2*k
    def r(i): return [1 if j >= i else 0 for j in range(1, k+1)] + [0]*k
    def s(i): return [0]*k + [1 if l <= i else 0 for l in range(1, k+1)]
    vecs = []; wt = [0]*n
    for i in range(1, k+1):
        vecs.append([x+y for x, y in zip(r(i), s(i))]); vecs.append(s(i))
        wt = [w+x+2*y for w, x, y in zip(wt, r(i), s(i))]
    return build(vecs, wt, [1]*k+[2]*k)
def Cdata(k):
    def sg(i): return [1 if l <= i else 0 for l in range(1, k+1)]
    vecs = [sg(i) for i in range(1, k+1)]
    wt = [0]*k
    for v in vecs: wt = [w+x for w, x in zip(wt, v)]
    return build(vecs, wt, [2]*k)
def splitall(M, B, d, labs):
    labs = labs[:]
    for a in range(len(d)):
        if d[a] == 1:
            M, B, d = do(M, B, d, ('p', a)); labs.append(labs[a]+"''"); labs[a] += "'"
    return M, B, d, labs
def showform(M, B, d, labs):
    print("   labels", labs, " d", d); print("   B =", [fmtB(b) for b in B])
    for row, l in zip(M, labs): print("     %-10s" % l, [str(x) for x in row])
if __name__ == "__main__":
    which = sys.argv[1]; k = int(sys.argv[2])
    if which == 'B':
        M, B, d = Bdata(k); labs = ["n%d" % i for i in range(1, k+2)]
    else:
        M, B, d = Kdata(k); labs = ["a%d" % i for i in range(1, k+1)] + ["b%d" % i for i in range(1, k+1)]
    M, B, d, labs = splitall(M, B, d, labs)
    print("split form rank", len(d)); showform(M, B, d, labs)
    T = terminals(M, B, d, labs)
    byrank = sorted(T.items(), key=lambda kv: len(kv[1][0]))
    print("terminals:", len(T), "ranks", sorted(set(len(v[0]) for v in T.values())))
    for key, (M2, B2, d2, l2, path) in byrank[:3]:
        print(" terminal rank", len(M2), "path", path); showform(M2, B2, d2, l2)
