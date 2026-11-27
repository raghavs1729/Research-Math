# base-Q-only bidirectional BFS (fusion / defusion) between two all-d=2 generalized Nahm data
import sys, itertools
sys.path.insert(0, '/home/user/weighted-prime-races/monotonicity/rr-identities/code')
from fractions import Fraction as Fr
from classsearch import do
from moves import canon, defusions
from bfs1 import build, X
from term1 import Bdata, Kdata, Cdata, splitall, showform
def tensor(D1, D2):
    (M1, B1, d1), (M2, B2, d2) = D1, D2
    n1, n2 = len(d1), len(d2)
    M = [[Fr(0)]*(n1+n2) for _ in range(n1+n2)]
    for i in range(n1):
        for j in range(n1): M[i][j] = M1[i][j]
    for i in range(n2):
        for j in range(n2): M[n1+i][n1+j] = M2[i][j]
    return M, list(B1)+list(B2), list(d1)+list(d2)
def exploreQ(M, B, d, maxrank, limit=300000, verbose=False):
    start = canon(M, B, d); seen = {start: (None, None)}; frontier = [(M, B, d, start)]; depth = 0
    while frontier and len(seen) < limit:
        nxt = []
        for (M, B, d, key) in frontier:
            mvs = defusions(M, B, d)
            if len(M) < maxrank:
                mvs = mvs + [('f', a, b) for a, b in itertools.combinations(range(len(M)), 2)]
            for mv in mvs:
                M2, B2, d2 = do(M, B, d, mv)
                k2 = canon(M2, B2, d2)
                if k2 not in seen:
                    seen[k2] = (key, mv); nxt.append((M2, B2, d2, k2))
        frontier = nxt; depth += 1
        if verbose: print("  depth", depth, "seen", len(seen), file=sys.stderr)
    return seen
def pathto(seen, key):
    p = []
    while seen[key][0] is not None:
        p.append(seen[key][1]); key = seen[key][0]
    return p[::-1]
def RR2data():
    # sum x^{2v} q^{2v^2}/(Q)_v : exponent 2 v^2 -> M=[[4]], weight 2
    return [[Fr(4)]], [2*X], [2]
if __name__ == "__main__":
    k = int(sys.argv[1]); capK = int(sys.argv[2]); capB = int(sys.argv[3])
    M, B, d = Kdata(k); labs = ["a%d" % i for i in range(1, k+1)] + ["b%d" % i for i in range(1, k+1)]
    MK, BK, dK, labsK = splitall(M, B, d, labs)
    sK = exploreQ(MK, BK, dK, capK, verbose=True); print("K side class:", len(sK), "min rank", min(len(key[0]) for key in sK))
    # target: (split B_{k+1}) x C_k / C_1 : use B_{k+1} x C_k and compare with K x C_1 instead
    MB, BB, dB = Bdata(k); labsB = ["n%d" % i for i in range(1, k+2)]
    MB, BB, dB, labsB = splitall(MB, BB, dB, labsB)
    # K x C_1
    MKC, BKC, dKC = tensor((MK, BK, dK), Cdata(1))
    MBC, BBC, dBC = tensor((MB, BB, dB), Cdata(k))
    s1 = exploreQ(MKC, BKC, dKC, capK+1, verbose=True); print("K*C1 class:", len(s1))
    s2 = exploreQ(MBC, BBC, dBC, capB, verbose=True); print("B*C class:", len(s2))
    common = set(s1) & set(s2); print("common:", len(common))
    if common:
        best = min(common, key=lambda key: len(pathto(s1, key))+len(pathto(s2, key)))
        print("meet rank", len(best[0])); print("path K*C1:", pathto(s1, best)); print("path B*C :", pathto(s2, best))
