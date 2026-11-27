import sys, itertools, pickle
sys.path.insert(0, '/home/user/weighted-prime-races/monotonicity/rr-identities/code')
from fractions import Fraction as Fr
from moves import canon, terminals
from bfs1 import X
from term1 import Bdata, Kdata, Cdata, splitall, showform
from bfsQ import tensor, exploreQ, pathto, RR2data
k = int(sys.argv[1]); capK = int(sys.argv[2]); capB = int(sys.argv[3]); lim = int(sys.argv[4]) if len(sys.argv) > 4 else 400000
M, B, d = Kdata(k); labs = ["a%d" % i for i in range(1, k+1)] + ["b%d" % i for i in range(1, k+1)]
MK, BK, dK, labsK = splitall(M, B, d, labs)
MB, BB, dB = Bdata(k); labsB = ["n%d" % i for i in range(1, k+2)]
MB, BB, dB, labsB = splitall(MB, BB, dB, labsB)
# J_k as base-Q data: k odd=2m-1: AG_{m-1}(x^2;Q) (rank m-1); here use C_k and compare with K*C_1 only if needed.
if k % 2 == 1:
    m = (k+1)//2
    # AG_{m-1}(x^2;q^2): sum x^{2 sum sig} q^{2 sum sig^2} / prod (Q)_{sig_i - sig_{i+1}}
    n = m-1
    def sg(i): return [1 if l >= i else 0 for l in range(1, n+1)]
    from bfs1 import build
    vecs = []
    wt = [0]*n
    for i in range(1, n+1):
        vecs.append(sg(i)); vecs.append(sg(i)); wt = [w+2*x for w, x in zip(wt, sg(i))]
    J = build(vecs, wt, [2]*n) if n > 0 else ([], [], [])
    # terminals of split B
    T = terminals(MB, BB, dB, labsB)
    targets = [v[:3] for v in T.values()]
    print("B terminals ranks:", [len(t[2]) for t in targets])
    sK = exploreQ(MK, BK, dK, capK, limit=lim, verbose=True); print("K class:", len(sK), "min rank", min(len(key[0]) for key in sK))
    for t in targets:
        D = tensor(t, J) if n > 0 else t
        sB = exploreQ(D[0], D[1], D[2], capB, limit=lim, verbose=True); print("B*J class:", len(sB))
        common = set(sK) & set(sB); print("common:", len(common))
        if common:
            best = min(common, key=lambda key: len(pathto(sK, key))+len(pathto(sB, key)))
            print("meet rank", len(best[0])); print("path K:", pathto(sK, best)); print("path B*J from terminal:", pathto(sB, best))
            print("K labels", labsK); 
            break
