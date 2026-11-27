import sys
sys.path.insert(0, '/home/user/weighted-prime-races/monotonicity/rr-identities/code')
from fractions import Fraction as Fr
from classsearch import do
from moves import canon, terminals
from bfs1 import X, build
from bfs2 import fmtB
from term1 import Bdata, Kdata, splitall
def run(name, data, path, labs, verbose=True):
    M, B, d = data; labs = labs[:]
    print("==", name)
    for mv in path:
        if mv[0] == 'f':
            new = "[%s+%s]" % (labs[mv[1]], labs[mv[2]])
            desc = "fuse %s , %s" % (labs[mv[1]], labs[mv[2]])
            M, B, d = do(M, B, d, mv); labs.append(new)
        elif mv[0] == 'F':
            c = mv[3]
            desc = "defuse %s  (as fusion of %s , %s)" % (labs[c], labs[mv[1]], labs[mv[2]])
            M, B, d = do(M, B, d, mv)
            la, lb, lc = labs[mv[1]], labs[mv[2]], labs[c]
            labs[mv[1]] = "(%s|%s)" % (la, lc); labs[mv[2]] = "(%s|%s)" % (lb, lc)
            labs = [l for i, l in enumerate(labs) if i != c]
        if verbose: print("  ", desc)
    print(" final labels:"); 
    for i, l in enumerate(labs): print("    v%d = %s   B=%s" % (i, l, fmtB(B[i])))
    for row in M: print("     ", [str(v) for v in row])
    return M, B, d, labs
k = 3
M, B, d = Kdata(k); labs = ["a%d" % i for i in range(1, k+1)] + ["b%d" % i for i in range(1, k+1)]
MK, BK, dK, labsK = splitall(M, B, d, labs)
pathK = [('f', 0, 4), ('F', 0, 6, 7), ('f', 1, 5), ('F', 0, 1, 2), ('f', 1, 3), ('F', 4, 8, 9), ('f', 3, 5), ('F', 0, 3, 7), ('F', 4, 6, 8)]
rK = run("K_3 split", (MK, BK, dK), pathK, labsK)

# B side
MB, BB, dB = Bdata(k); labsB = ["n%d" % i for i in range(1, k+2)]
MB, BB, dB, labsB = splitall(MB, BB, dB, labsB)
T = terminals(MB, BB, dB, labsB)
from bfsQ import tensor, RR2data
for key, (M2, B2, d2, l2, path) in T.items():
    print("B terminal labels", l2, "path", path)
    D = tensor((M2, B2, d2), RR2data())
    rB = run("B_4 terminal x RR2", D, [('f', 2, 5), ('f', 1, 5)], l2 + ["sg"])
    print("canon equal:", canon(*rB[:3]) == canon(*rK[:3]))
