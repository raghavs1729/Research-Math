import sys
sys.path.insert(0, '/home/user/weighted-prime-races/monotonicity/rr-identities/code')
from fractions import Fraction as Fr
from classsearch import do
from moves import canon
from bfs1 import KC, BC, X
def fmtB(b):
    # b = lin + X*w
    w = round(float(b/X)); lin = b - X*w
    return "%s+%sX" % (lin, w)
def run(name, data, path, labs):
    M, B, d = data
    print("==", name); labs = labs[:]
    for mv in path:
        if mv[0] == 'p':
            a = mv[1]; labs.append(labs[a]+"''"); labs[a] = labs[a]+"'"
        elif mv[0] == 'f':
            labs.append("(%s+%s)" % (labs[mv[1]], labs[mv[2]]))
        elif mv[0] == 'F':
            c = mv[3]; print("   defuse", labs[c], "= fusion of", labs[mv[1]], labs[mv[2]]); labs = [l for i, l in enumerate(labs) if i != c]
        M, B, d = do(M, B, d, mv)
        print(" after", mv, "labels", labs)
    print(" d =", d); print(" B =", [fmtB(b) for b in B])
    for row, l in zip(M, labs): print("   %-14s" % l, [str(x) for x in row])
    return M, B, d, labs
if __name__ == "__main__":
    r1 = run("K_2*C_1", KC(2), [('p', 0), ('f', 0, 3), ('p', 1), ('F', 0, 5, 7), ('F', 0, 1, 6)], ['a1', 'a2', 'b1', 'b2', 't'])
    r2 = run("B_3*C_2", BC(2), [('p', 0), ('p', 1), ('F', 0, 5, 6), ('F', 0, 1, 2), ('f', 1, 2)], ['n1', 'n2', 'n3', 'e1', 'e2'])
    print(canon(*r1[:3]) == canon(*r2[:3]))
