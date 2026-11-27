import sys
from flint import fmpz_poly
import krtw_ser as K
from krtw_ser import S, tr, mono, iqpoch, qpoch, qbin, xpoch, xpoch_inv, xpoch_inf

def rmax_for(N):
    r = 1
    while r*r < N: r += 1
    return r

_A = {}
def A(j, R, lin=0):
    """A_j(R;x) = sum (q)_R/(q)_{R-r_1} prod x^{r_i+2s_i} q^{(r_i+s_i)^2+s_i^2 + lin*(r_i+2s_i)}/..., as S"""
    key = (j, R, lin, K.N)
    if key in _A: return _A[key]
    if j == 0:
        res = S.one()
    else:
        res = S()
        for r in range(R+1):
            b = 0
            while True:
                e = (r+b)**2 + b*b
                if e >= K.N: break
                # (x^{2b} q^{2b^2})^{j-1} A_{j-1}(r; x q^{2b})
                inner = A(j-1, r, lin).subx(2*b).shift(2*b*(j-1), 2*b*b*(j-1))
                c = tr(qbin(R, r) * iqpoch(b, 2))
                res = res + (inner * c).shift(r+2*b, e)
                b += 1
        if lin:
            pass
    _A[key] = res
    return res

def Kk(j):
    return A(j, rmax_for(K.N)+2)

def calK(j, i):
    """K_j(i;x) = sum_{r>=i,b} x^{r-i+2b} q^{(r+b)^2+b^2-ri} (x^{2b}q^{2b^2})^{j-1} A_{j-1}(r;xq^{2b}) /((q)_{r-i}(q^2;q^2)_b)"""
    res = S()
    r = i
    while True:
        if r*r - r*i >= K.N: break
        b = 0
        while True:
            e = (r+b)**2 + b*b - r*i
            if e >= K.N: break
            inner = A(j-1, r).subx(2*b).shift(2*b*(j-1), 2*b*b*(j-1))
            c = tr(iqpoch(r-i) * iqpoch(b, 2))
            res = res + (inner * c).shift(r-i+2*b, e)
            b += 1
        r += 1
    return res

def bress(nvar, lin=None):
    """sum x^{sum N_i} q^{sum N_i^2 + lin(N)} / ((q)_{n_1}...(q)_{n_{nvar-1}} (q^2;q^2)_{N_nvar}); lin: function of list N -> int"""
    s = S()
    N = K.N
    def rec(Ns, e, c):
        nonlocal s
        lev = len(Ns)
        if lev == nvar:
            ee = e + (lin(Ns) if lin else 0)
            if ee < N:
                s = s + S.xq(sum(Ns), ee) * tr(c * iqpoch(Ns[-1], 2))
            return
        M = 0
        while lev == 0 or M <= Ns[-1]:
            ee = e + M*M
            if ee >= N: break
            cc = c if lev == 0 else tr(c * iqpoch(Ns[-1]-M))
            rec(Ns+[M], ee, cc)
            M += 1
    rec([], 0, fmpz_poly([1]))
    return s

def Cchain(k, lin=None):
    """C_k(x) = sum_{0<=s_1<=...<=s_k} x^{sum s} q^{sum s^2 + lin(s)} / prod (q^2;q^2)_{s_i-s_{i-1}}"""
    s = S()
    N = K.N
    def rec(ss, e, c):
        nonlocal s
        lev = len(ss)
        if lev == k:
            ee = e + (lin(ss) if lin else 0)
            if ee < N:
                s = s + S.xq(sum(ss), ee) * c
            return
        prev = ss[-1] if ss else 0
        M = prev
        while True:
            # remaining levels all >= M
            ee = e + M*M
            if e + (k-lev)*M*M >= N: break
            rec(ss+[M], ee, tr(c*iqpoch(M-prev, 2)))
            M += 1
    rec([], 0, fmpz_poly([1]))
    return s

def Sinv(a):
    """1/a for S with constant term 1 (x^0 coefficient = 1 exactly) -- via geometric series in (1-a)"""
    one = S.one()
    u = one - a   # has x-degree>=1 terms only? require a.d[0]==1
    assert a.d.get(0) == fmpz_poly([1]), "need x^0 coeff 1"
    res = S.one(); p = S.one()
    # u has only x-degree >= 1 terms; each x-degree carries q-power >= 1 -> terminates
    while True:
        p = p * u
        if p.iszero(): break
        res = res + p
    return res
