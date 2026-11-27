# Referee: replay the move chains of KRTW-F-moves.md section 3 with an independent implementation,
# and evaluate EVERY intermediate datum by brute force against direct evaluations.
import sys, itertools, time
from fractions import Fraction as Fr
from datum import D, split, fuse, defuse, direct_sum, evaluate
import direct
from ser import S

X = int(sys.argv[1]) if len(sys.argv) > 1 else 10
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40

def quad_datum(names, forms, base, weight):
    """exponent = sum_f (form . n)^2 ; forms: list of dict name->coef"""
    M = {}
    for i, a in enumerate(names):
        for b in names[i:]:
            M[(a, b)] = 2 * sum(f.get(a, 0) * f.get(b, 0) for f in forms)
    return D(names, M, {a: 0 for a in names}, weight, base)

def K_datum(k):
    A = ['a%d' % j for j in range(1, k + 1)]; Bn = ['b%d' % l for l in range(1, k + 1)]
    forms = []
    for i in range(1, k + 1):
        r = {('a%d' % j): 1 for j in range(i, k + 1)}
        s = {('b%d' % l): 1 for l in range(1, i + 1)}
        L = dict(r); L.update(s); forms.append(L); forms.append(s)
    w = {}
    for j in range(1, k + 1): w['a%d' % j] = j
    for l in range(1, k + 1): w['b%d' % l] = 2 * (k - l + 1)
    base = {a: 1 for a in A}; base.update({b: 2 for b in Bn})
    return quad_datum(A + Bn, forms, base, w)

def B_datum(m):
    names = ['n%d' % j for j in range(1, m + 1)]
    forms = [{('n%d' % j): 1 for j in range(i, m + 1)} for i in range(1, m + 1)]
    w = {('n%d' % j): j for j in range(1, m + 1)}
    base = {('n%d' % j): 1 for j in range(1, m)}; base['n%d' % m] = 2
    return quad_datum(names, forms, base, w)

def C_datum(k):
    names = ['e%d' % l for l in range(1, k + 1)]
    forms = [{('e%d' % l): 1 for l in range(1, i + 1)} for i in range(1, k + 1)]
    w = {('e%d' % l): k - l + 1 for l in range(1, k + 1)}
    return quad_datum(names, forms, {a: 2 for a in names}, w)

def AG1_datum():
    return D(['g1'], {('g1', 'g1'): 4}, {'g1': 0}, {'g1': 2}, {'g1': 2})

log = []
def check(Dt, target, label):
    t0 = time.time()
    F, mn = evaluate(Dt, X, N)
    ok = (F == target)
    msg = "%-60s rank %2d  minexp %s  %s  (%.1fs)" % (label, len(Dt.names), mn, "OK" if ok else "MISMATCH", time.time() - t0)
    print(msg); sys.stdout.flush(); log.append((label, ok))
    return ok

def DF(Dt, c, a, b):
    return defuse(Dt, a, b, c, '(%s|%s)' % (a, c), '(%s|%s)' % (b, c))
def FU(Dt, a, b):
    return fuse(Dt, a, b, '[%s+%s]' % (a, b))

def splitall(Dt, names):
    for a in names:
        Dt = split(Dt, a, a + "''")
        Dt.rename(a, a + "'")
    return Dt

print("X =", X, " N =", N)
# ---------- targets by direct evaluation
K2 = direct.K(2, X, N); K3 = direct.K(3, X, N)
B3 = direct.B(3, X, N); B4 = direct.B(4, X, N)
C1 = direct.C(1, X, N); C2 = direct.C(2, X, N); C3 = direct.C(3, X, N)
AG = direct.AG1sq(X, N)
print("direct: K3 == B4*C3/C1 :", K3 == B4 * C3 * C1.inv())
print("direct: J3 == AG1(x^2;Q):", C3 * C1.inv() == AG)
print("direct: K2 == B3*C2/C1 :", K2 == B3 * C2 * C1.inv())

# ---------- k = 2, K side
Dk = K_datum(2)
check(Dk, K2, "K2 start")
Dk = splitall(Dk, ['a1', 'a2'])
check(Dk, K2, "K2 after SPLIT")
Dk = FU(Dk, "a1'", "b2"); check(Dk, K2, "K2 FUSE a1',b2")
Dk = DF(Dk, "a2''", "a1'", "a1''"); check(Dk, K2, "K2 DEFUSE a2''=a1'+a1''")
Dk = DF(Dk, "[a1'+b2]", "(a1'|a2'')", "a2'"); check(Dk, K2, "K2 DEFUSE [a1'+b2]=(a1'|a2'')+a2'")
K2final = Dk
# ---------- k = 2, B side
Db = B_datum(3); check(Db, B3, "B3 start")
Db = splitall(Db, ['n1', 'n2']); check(Db, B3, "B3 after SPLIT")
Db = DF(Db, "n2''", "n1'", "n1''"); check(Db, B3, "B3 DEFUSE n2''=n1'+n1''")
Db = DF(Db, "n3", "(n1'|n2'')", "n2'"); check(Db, B3, "B3 DEFUSE n3=(n1'|n2'')+n2'")
free = "((n1'|n2'')|n3)"
print("  free species row:", [(k, Db.m(free, k)) for k in Db.names], "B", Db.B[free], "w", Db.w[free], "d", Db.d[free])
Db.remove(free)
check(Db, B3 * C1.inv(), "B3 datum minus free species  vs B3/C1")
Db = direct_sum(Db, C_datum(2)); check(Db, B3 * C2 * C1.inv(), "(B3/C1) (+) C2  vs B3 J2")
Db = FU(Db, "(n2'|n3)", "e1"); check(Db, K2, "FUSE (n2'|n3),e1  vs K2")
B2final = Db

def perm_match(D1, D2):
    if len(D1.names) != len(D2.names): return None
    c1 = D1.canon()
    for p in itertools.permutations(D2.names):
        if D2.canon(list(p)) == c1:
            return list(zip(D1.names, p))
    return None
print("k=2 final data match by permutation:", perm_match(K2final, B2final))

# ---------- k = 3, K side
Dk = K_datum(3); check(Dk, K3, "K3 start")
Dk = splitall(Dk, ['a1', 'a2', 'a3']); check(Dk, K3, "K3 after SPLIT")
steps = [
    ('F', "a1'", "b2"),
    ('D', "a2''", "a1'", "a1''"),
    ('F', "a2'", "b3"),
    ('D', "a3'", "(a1'|a2'')", "a2'"),
    ('F', "(a2'|a3')", "b2"),
    ('D', "[(a2'|a3')+b2]", "b3", "[a2'+b3]"),
    ('F', "b2", "(a1''|a2'')"),
    ('D', "[a1'+b2]", "((a1'|a2'')|a3')", "b2"),
    ('D', "[b2+(a1''|a2'')]", "(b3|[(a2'|a3')+b2])", "a3''"),
]
for st in steps:
    if st[0] == 'F':
        Dk = FU(Dk, st[1], st[2]); lab = "K3 FUSE %s,%s" % (st[1], st[2])
    else:
        Dk = DF(Dk, st[1], st[2], st[3]); lab = "K3 DEFUSE %s=%s+%s" % st[1:]
    check(Dk, K3, lab)
K3final = Dk
# ---------- k = 3, B side
Db = B_datum(4); check(Db, B4, "B4 start")
Db = splitall(Db, ['n1', 'n2', 'n3']); check(Db, B4, "B4 after SPLIT")
Db = DF(Db, "n2''", "n1'", "n1''"); check(Db, B4, "B4 DEFUSE n2''=n1'+n1''")
Db = DF(Db, "n3'", "(n1'|n2'')", "n2'"); check(Db, B4, "B4 DEFUSE n3'=(n1'|n2'')+n2'")
Db = direct_sum(Db, AG1_datum()); check(Db, B4 * AG, "B4 (+) AG1  vs B4*AG1")
Db = FU(Db, "n4", "g1"); check(Db, B4 * AG, "FUSE n4,g1")
Db = FU(Db, "(n2'|n3')", "g1"); check(Db, K3, "FUSE (n2'|n3'),g1  vs K3")
B3final = Db
print("k=3 final data match by permutation:", perm_match(K3final, B3final))
print("ALL OK:", all(ok for _, ok in log), " checks:", len(log))
