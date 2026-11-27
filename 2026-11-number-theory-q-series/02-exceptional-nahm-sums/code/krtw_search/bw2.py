# Check: BW multisum at x=(alpha,-alpha) equals K_k(x).  BW_k = sum over (a_l,b_l), l=1..k of
#   prod_l x^{a_l+b_l} q^{a_l^2+b_l^2} omega(a_l-b_l) / ((q)_{a_l-a_{l+1}} (q)_{b_l-b_{l+1}} (-q)_{a_l-b_{l+1}} (-q)_{b_l-a_{l+1}})
import sys
from lib2 import *
from bw1 import mqpoch, imq, omega
N = int(sys.argv[1]) if len(sys.argv)>1 else 22
K.setN(N)
def BW(k):
    amax = rmax_for(K.N)+1
    tot = S()
    # recursion from the bottom: state (a,b) at level l, lower level (a',b') ; start with level k+1 = (0,0)
    def rec(l, low, e, c, sh, xdeg, twos):
        nonlocal tot
        ap, bp = low
        if l == 0:
            ee = e - sh
            if ee < K.N:
                assert ee >= 0
                tot = tot + S.xq(xdeg, ee) * c * (2**twos)
            return
        for a in range(ap, amax+1):
            if e + a*a - sh - 2*amax*amax >= K.N and False: break
            for b in range(bp, amax+1):
                ee = e + a*a + b*b
                # crude pruning: later shifts only reduce exponent by sh' ; estimate
                cc = tr(c*iqpoch(a-ap)); cc = tr(cc*iqpoch(b-bp))
                p1, s1 = imq(a-bp); p2, s2 = imq(b-ap)
                on, od = omega(a-b)
                cc0 = tr(tr(cc*p1)*p2)
                cc0 = tr(cc0*on)
                if ee - (sh+s1+s2) >= K.N + 4*amax: continue
                rec(l-1, (a, b), ee, cc0, sh+s1+s2, xdeg+a+b, twos + (0 if od == 2 else 0) - (1 if od == 2 else 0) )
    rec(k, (0, 0), 0, fmpz_poly([1]), 0, 0, 0)
    return tot
# note: factors: imq for negative index returns 2*(...) with shift; omega has denominator 2 when d != 0. Track powers of 2 separately:
def BW2(k):
    """exact with rational 2-powers handled: multiply all by 2^k then divide"""
    amax = rmax_for(K.N)+1
    tot = S()
    def rec(l, low, e, c, sh, xdeg):
        nonlocal tot
        ap, bp = low
        if l == 0:
            ee = e - sh
            if 0 <= ee < K.N: tot = tot + S.xq(xdeg, ee) * c
            return
        for a in range(ap, amax+1):
            for b in range(bp, amax+1):
                ee = e + a*a + b*b
                cc = tr(c*iqpoch(a-ap)); cc = tr(cc*iqpoch(b-bp))
                p1, s1 = imq(a-bp); p2, s2 = imq(b-ap)
                on, od = omega(a-b)
                cc = tr(tr(tr(cc*p1)*p2)*on)
                if od == 1: cc = cc*2
                if ee - (sh+s1+s2) >= K.N + 2*amax*amax: continue
                rec(l-1, (a, b), ee, cc, sh+s1+s2, xdeg+a+b)
    rec(k, (0, 0), 0, fmpz_poly([1]), 0, 0)
    return tot   # = 2^k * BW_k
for k in (1, 2):
    b = BW2(k); kk = calK(k, 0)
    print("k=%d: 2^k*BW_k == 2^k*K_k ?" % k, b == kk*(2**k))
    if not (b == kk*(2**k)):
        for d in sorted(b.d)[:4]: print("   BW x^%d" % d, b.d[d].coeffs()[:10], " K:", (kk.d.get(d) or fmpz_poly([0])).coeffs()[:10])
