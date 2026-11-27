import sys
sys.path.insert(0, '/home/user/weighted-prime-races/monotonicity/rr-identities/code')
from fractions import Fraction as Fr
from classsearch import explore, pathto
from moves import canon
X = Fr(1000003, 7)   # formal marker for the x-grading
def build(vecs_sq, wt, d):
    """exponent = sum_v (v.n)^2 ; x-weight = wt.n ; returns M,B,d"""
    n = len(d)
    M = [[Fr(0)]*n for _ in range(n)]
    for v in vecs_sq:
        for p in range(n):
            for q in range(n):
                M[p][q] += 2*v[p]*v[q]
    B = [X*w for w in wt]
    return M, B, list(d)
def KC(k):
    # vars: a_1..a_k (d=1), b_1..b_k (d=2), t (d=2)
    n = 2*k+1
    def r(i): return [1 if j >= i else 0 for j in range(1, k+1)] + [0]*k + [0]
    def s(i): return [0]*k + [1 if l <= i else 0 for l in range(1, k+1)] + [0]
    t = [0]*(2*k) + [1]
    vecs = []; wt = [0]*n
    for i in range(1, k+1):
        vecs.append([x+y for x, y in zip(r(i), s(i))]); vecs.append(s(i))
        wt = [w+x+2*y for w, x, y in zip(wt, r(i), s(i))]
    vecs.append(t); wt = [w+x for w, x in zip(wt, t)]
    return build(vecs, wt, [1]*k+[2]*k+[2])
def BC(k):
    # vars: n_1..n_k (d=1), n_{k+1} (d=2), e_1..e_k (d=2)
    n = 2*k+1
    def Nv(i): return [1 if j >= i else 0 for j in range(1, k+2)] + [0]*k
    def sg(i): return [0]*(k+1) + [1 if l <= i else 0 for l in range(1, k+1)]
    vecs = [Nv(i) for i in range(1, k+2)] + [sg(i) for i in range(1, k+1)]
    wt = [0]*n
    for v in vecs: wt = [w+x for w, x in zip(wt, v)]
    return build(vecs, wt, [1]*k+[2]+[2]*k)
if __name__ == '__main__':
    k = int(sys.argv[1]); maxrank = int(sys.argv[2]); limit = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
    M1, B1, d1 = KC(k); M2, B2, d2 = BC(k)
    print("start canon equal?", canon(M1, B1, d1) == canon(M2, B2, d2))
    s1 = explore(M1, B1, d1, maxrank, 4, limit); print("class of K*C1:", len(s1))
    s2 = explore(M2, B2, d2, maxrank, 4, limit); print("class of B*C :", len(s2))
    common = set(s1) & set(s2)
    print("common:", len(common))
    if common:
        best = min(common, key=lambda key: len(pathto(s1, key))+len(pathto(s2, key)))
        print("meeting point rank", len(best[0]), "d", best[0])
        print("path from K*C1:", pathto(s1, best))
        print("path from B*C :", pathto(s2, best))
