# Compare A_k(R;x) with the C_2 Bailey-chain beta at x=(alpha,-alpha), N=(ceil(R/2), floor(R/2))
import sys
from lib2 import *
N = int(sys.argv[1]) if len(sys.argv)>1 else 24
K.setN(N)
def mqpoch(n):
    """(-q;q)_n for n>=0 as fmpz_poly"""
    p = fmpz_poly([1])
    for i in range(1, n+1):
        p = tr(p*(fmpz_poly([1]) + mono(i)))
    return p
_im = {}
def imq(n):
    """1/(-q;q)_n for n>=0 ; for n<0: (-q^{1+n};q)_{-n} = prod_{i=0}^{-n-1}(1+q^{1+n+i}) -> Laurent; return (poly, shift) with value = poly*q^{-shift}"""
    if n >= 0:
        if n not in _im: _im[n] = K.inv(mqpoch(n))
        return _im[n], 0
    k = -n
    # prod_{j=0}^{k-1} (1+q^{-j}) = q^{-k(k-1)/2} prod_{j=0}^{k-1}(1+q^j) = q^{-C(k,2)} * 2 * (-q;q)_{k-1}
    return 2*mqpoch(k-1), k*(k-1)//2
def omega(d):
    d = abs(d)
    if d == 0: return fmpz_poly([1]), 1
    return tr(mono(d*(d-1)//2)*(fmpz_poly([1])+mono(d))), 2   # numerator, denominator 2
def betaBW(k, N1, N2):
    """returns S with coefficients times 2^k (to stay integral): sum over chains (a_l,b_l), l=1..k"""
    tot = S()
    def rec(l, prev, e, c, sh, xdeg):
        nonlocal tot
        # prev=(A,B) are the upper indices (N for l=1, else r^{(l-1)}); choose (a,b)<= componentwise? need a<=A (q)_{A-a}, b<=B
        A, B = prev
        if l == k+1:
            # close with r^{(k+1)}=0: factors 1/((q)_A (q)_B (-q)_A (-q)_B)
            cc = tr(c*iqpoch(A)); cc = tr(cc*iqpoch(B))
            p1, s1 = imq(A); p2, s2 = imq(B)
            cc = tr(tr(cc*p1)*p2)
            ee = e - sh
            assert ee >= 0
            if ee < K.N:
                tot = tot + S.xq(xdeg, ee) * cc
            return
        for a in range(A+1):
            for b in range(B+1):
                ee = e + a*a + b*b
                if ee - sh - 50 >= K.N: continue
                cc = tr(c*iqpoch(A-a)); cc = tr(cc*iqpoch(B-b))
                p1, s1 = imq(A-b); p2, s2 = imq(B-a)
                cc = tr(tr(cc*p1)*p2)
                on, od = omega(a-b)
                cc = tr(cc*on)
                if od == 1: cc = cc*2
                rec(l+1, (a, b), ee, cc, sh+s1+s2, xdeg+a+b)
    # top level: links from N to r^{(1)}; the omega factor belongs to each r^{(l)}
    rec(1, (N1, N2), 0, fmpz_poly([1]), 0, 0)
    return tot
for k in (1, 2):
    for R in range(0, 5):
        N1, N2 = (R+1)//2, R//2
        b = betaBW(k, N1, N2)          # times 2^k
        a = A(k, R)
        # ratio at x^0
        c0 = b.d[0]
        ratio = a * K.inv(tr(c0*K.inv(fmpz_poly([c0.coeffs()[0]])))) if False else None
        # compare a*(b.d[0]) with b*(a.d[0])   (a.d[0]=1)
        lhs = a * c0
        print("k=%d R=%d N=(%d,%d): proportional:"%(k, R, N1, N2), lhs == b, " beta x^0 coeff:", c0.coeffs()[:8])
