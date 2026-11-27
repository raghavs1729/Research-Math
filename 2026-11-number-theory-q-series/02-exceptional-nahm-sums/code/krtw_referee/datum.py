# Independent implementation of generalized Nahm data and the moves SPLIT / FUSE / DEFUSE,
# written from the statements of Lemmas M1, M2 (not from the project's code).
from fractions import Fraction as Fr
from ser import S, qpoch_inv_table
from flint import fmpz_poly
import itertools

class D:
    def __init__(self, names, M, B, w, d):
        self.names = list(names)
        self.M = {(a, b): Fr(v) for (a, b), v in M.items()}
        self.B = {a: Fr(v) for a, v in B.items()}
        self.w = dict(w); self.d = dict(d)
    def m(self, a, b):
        return self.M.get((a, b), self.M.get((b, a), Fr(0)))
    def copy(self):
        return D(self.names, self.M, self.B, self.w, self.d)
    def setm(self, a, b, v):
        self.M.pop((a, b), None); self.M.pop((b, a), None)
        self.M[(a, b)] = Fr(v)
    def rename(self, old, new):
        self.names = [new if n == old else n for n in self.names]
        M2 = {}
        for (a, b), v in self.M.items():
            M2[(new if a == old else a, new if b == old else b)] = v
        self.M = M2
        for dic in (self.B, self.w, self.d):
            dic[new] = dic.pop(old)
    def remove(self, c):
        self.names.remove(c)
        self.M = {(a, b): v for (a, b), v in self.M.items() if a != c and b != c}
        for dic in (self.B, self.w, self.d):
            dic.pop(c)
    def matrix(self):
        return [[self.m(a, b) for b in self.names] for a in self.names]
    def canon(self, order=None):
        o = order or self.names
        return (tuple(tuple(self.m(a, b) for b in o) for a in o),
                tuple(self.B[a] for a in o), tuple(self.w[a] for a in o), tuple(self.d[a] for a in o))

def split(D0, a, t):
    D1 = D0.copy()
    dd = D1.d[a]
    D1.names.insert(D1.names.index(a) + 1, t)
    for k in D0.names:
        if k != a:
            D1.setm(t, k, D0.m(a, k))
    D1.setm(t, t, D0.m(a, a)); D1.setm(a, t, D0.m(a, a))
    D1.d[a] = 2 * dd; D1.d[t] = 2 * dd
    D1.B[t] = D0.B[a] + dd; D1.w[t] = D0.w[a]
    return D1

def fuse(D0, a, b, c):
    assert D0.d[a] == D0.d[b]
    dd = D0.d[a]
    D1 = D0.copy()
    D1.names.append(c)
    D1.d[c] = dd; D1.w[c] = D0.w[a] + D0.w[b]; D1.B[c] = D0.B[a] + D0.B[b]
    D1.setm(c, c, D0.m(a, a) + D0.m(b, b) + 2 * D0.m(a, b))
    D1.setm(a, c, D0.m(a, a) + D0.m(a, b))
    D1.setm(b, c, D0.m(b, b) + D0.m(a, b))
    for k in D0.names:
        if k not in (a, b):
            D1.setm(c, k, D0.m(a, k) + D0.m(b, k))
    D1.setm(a, b, D0.m(a, b) + dd)
    return D1

def defuse(D1, a, b, c, newa, newb):
    """inverse of fuse: c removed, a,b renamed newa,newb. Checks the precondition."""
    dd = D1.d[a]
    assert D1.d[b] == dd and D1.d[c] == dd, "base mismatch"
    Mab = D1.m(a, b) - dd
    ok = (D1.w[c] == D1.w[a] + D1.w[b] and D1.B[c] == D1.B[a] + D1.B[b]
          and D1.m(c, c) == D1.m(a, a) + D1.m(b, b) + 2 * Mab
          and D1.m(a, c) == D1.m(a, a) + Mab and D1.m(b, c) == D1.m(b, b) + Mab
          and all(D1.m(c, k) == D1.m(a, k) + D1.m(b, k) for k in D1.names if k not in (a, b, c)))
    if not ok:
        raise AssertionError("DEFUSE precondition fails for %s,%s,%s" % (a, b, c))
    D0 = D1.copy()
    D0.remove(c)
    D0.setm(a, b, Mab)
    D0.rename(a, newa); D0.rename(b, newb)
    return D0

def direct_sum(D1, D2):
    names = D1.names + D2.names
    M = dict(D1.M); M.update(D2.M)
    B = dict(D1.B); B.update(D2.B)
    w = dict(D1.w); w.update(D2.w)
    d = dict(D1.d); d.update(D2.d)
    return D(names, M, B, w, d)

def evaluate(Dt, X, N):
    """F(D) mod x^{X+1}, q^N by brute-force enumeration; raises on negative q-exponent."""
    names = Dt.names
    r = len(names)
    Mm = Dt.matrix()
    Bv = [Dt.B[a] for a in names]; wv = [Dt.w[a] for a in names]; dv = [Dt.d[a] for a in names]
    tabs = {dd: qpoch_inv_table(dd, X, N) for dd in set(dv)}
    out = S(X, N)
    minexp = [None]
    n = [0] * r
    def rec(i, wleft, prod):
        if i == r:
            e = Fr(0)
            for a in range(r):
                if n[a]:
                    e += Bv[a] * n[a] + Mm[a][a] * n[a] * n[a] / 2
                    for b in range(a + 1, r):
                        if n[b]: e += Mm[a][b] * n[a] * n[b]
            assert e.denominator == 1, "non-integral exponent"
            e = int(e)
            if minexp[0] is None or e < minexp[0]: minexp[0] = e
            if e < 0: raise ValueError("negative exponent %d at %s" % (e, n))
            if e < N:
                xd = X - wleft
                out.c[xd] += prod.left_shift(e).truncate(N) if e else prod
            return
        for k in range(0, wleft // wv[i] + 1):
            n[i] = k
            p = prod if k == 0 else prod.mul_low(tabs[dv[i]][k], N)
            rec(i + 1, wleft - k * wv[i], p)
        n[i] = 0
    rec(0, X, fmpz_poly([1]))
    return out, minexp[0]
