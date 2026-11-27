# two-parameter test: K_k(a,b) = sum a^{sum M_i} b^{sum s_i} q^{sum M_i^2+s_i^2}/..., a = x, b = w*x (w integer)
import sys
from lib2 import *
N = int(sys.argv[1]) if len(sys.argv)>1 else 26
K.setN(N)
w = 3
def K2ab(k):
    # brute force over r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k ; M_i=r_i+s_i
    tot = S()
    def rec(i, rs, ss, e, c):
        nonlocal tot
        if i == k:
            # close: factor 1/(q)_{r_k}
            c2 = tr(c*iqpoch(rs[-1]))
            Ms = sum(r+s for r, s in zip(rs, ss)); Ss = sum(ss)
            tot = tot + S.xq(Ms+Ss, e, w**Ss) * c2
            return
        rprev = rs[-1] if rs else None
        sprev = ss[-1] if ss else 0
        s = sprev
        while True:
            if e + (k-i)*2*s*s >= K.N: break
            r = 0
            while rprev is None or r <= rprev:
                ee = e + (r+s)**2 + s*s
                if ee >= K.N: break
                cc = tr(c*iqpoch(s-sprev, 2))
                if rprev is not None: cc = tr(cc*iqpoch(rprev-r))
                rec(i+1, rs+[r], ss+[s], ee, cc)
                r += 1
            s += 1
    rec(0, [], [], 0, fmpz_poly([1]))
    return tot
def bressw(nvar, wts):
    """sum prod (wts[i])^{N_i} x^{N_i} q^{N_i^2} /((q)...(q^2;q^2)_{N_nvar})"""
    s = S()
    def rec(Ns, e, c, coef):
        nonlocal s
        lev = len(Ns)
        if lev == nvar:
            s = s + S.xq(sum(Ns), e, coef) * tr(c*iqpoch(Ns[-1], 2)); return
        M = 0
        while lev == 0 or M <= Ns[-1]:
            ee = e + M*M
            if ee >= K.N: break
            cc = c if lev == 0 else tr(c*iqpoch(Ns[-1]-M))
            rec(Ns+[M], ee, cc, coef*wts[lev]**M)
            M += 1
    rec([], 0, fmpz_poly([1]), 1)
    return s
def scale(Sx, c):
    return S({k: v*(c**k) for k, v in Sx.d.items()})
k = 2
Kab = K2ab(k)
J = Cchain(k)*Sinv(Cchain(1))
print("sanity w->1 handled separately; testing candidates for K_2(x, %d x):"%w)
import itertools
for wts in itertools.product([1, w], repeat=k+1):
    B = bressw(k+1, list(wts))
    for nm, Jv in (("J(x)", J), ("J(wx)", scale(J, w))):
        print(wts, nm, Kab == B*Jv)
# also print low-order of Kab/ (B with wts (1,1,w)) to see
