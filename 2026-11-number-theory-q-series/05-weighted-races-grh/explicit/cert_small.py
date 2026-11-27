#!/usr/bin/env python3
"""
cert_small.py -- certified (exact integer / Arb ball arithmetic) versions of the finite checks
that the proofs in EXPLICIT.md actually use.  No floating-point decision anywhere in this file.

Usage: python3 cert_small.py [q ...]        (default: 4 3)

 A. EXACT.  Theta1(p) > 0 for every race prime p0 < p <= 10^4, as an inequality between integers:
        Theta1(p) = sum_{r<=p} s(r) log r (p-r) > 0   <=>   prod_{s(r)=+1} r^(p-r)  >  prod_{s(r)=-1} r^(p-r).
    (The two products are updated from one prime to the next: P(p') = P(p) * (prod_{r<=p} r)^(p'-p).)
    This is all that Theorem 1 needs from a computation (Step 7): Theta1 is piecewise linear with nodes
    at the race primes, so positivity at the nodes is positivity on (p0, 10^4].

 B. BALLS (python-flint / Arb, 256 bits).  For all race primes p <= NB = 2^22 = 4 194 304:
        Theta1(p) > 0 certified;  certified enclosures of min Theta1(p) and min Theta1(p)/p^{3/2}.

 C. BALLS.  Lemma 4(a) on (0, NB]:   y - c1 sqrt(y) <= vo(y) < y   with c1 = 2.071193 (q=4), 2.081942 (q=3):
        at every prime p <= NB:  (p - vo(p^-))/sqrt(p) <= c1  certified  (the sup over a gap is its right end),
                                 vo(p) < p certified               (the sup of vo(y) - y over a gap is its left end).
    Also certified enclosures of the sup (attained as y -> 1423^-) and of the sup over y >= 1427.

 D. BALLS + exact rationals.  Corollary 3 on the small range:  N(x) = R(x) log x/sqrt x > 1/2
        q=4:  on [23,25] and [25,29] (R >= R(left end), log x/sqrt x >= value at right end), and on every gap
              [p,p'] with 29 <= p < p' <= NB:   min(R(p),R(p')) * log p'/sqrt p' > 1/2;
        q=3:  on [6,7] (R >= R(6), log x/sqrt x >= log 6/sqrt 6 since x <= e^2), on [7,11], and on every gap
              with 11 <= p < p' <= NB.
    R(p) = D0 - D1/p is an exact rational.  Certified enclosure of the smallest gap bound.

 E. BALLS (64 bits; D0, D1 exact integers).  Corollary 3 on the whole range not covered by the analytic argument:
        q=4:  every gap [p,p'] with 29 <= p < p' <= 1.6e8:   R(p) w(p') > 1/2 and R(p') w(p') > 1/2,  w(x) = log x/sqrt x;
        q=3:  every gap [p,p'] with 11 <= p < p' <= 2.6e7.
    (On the gap R is monotone and w is decreasing, so N(x) >= min(R(p),R(p')) w(p').)
    Part E is the quick version of this check; cert_gaps.py carries it to 10^10 with a segmented sieve,
    and that is what the proof of Corollary 3 in EXPLICIT.md uses.
"""
import sys, time
from fractions import Fraction
from flint import arb, fmpz, ctx

ctx.prec = 256
NA = 10 ** 4
NB = 1 << 22
XE = {4: 160000000, 3: 26000000}     # part E: quick version; superseded by cert_gaps.py, which carries the same check to 1e10
C1 = {4: '2.071193', 3: '2.081942'}

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

ALLP = primes_upto(NB)

def sgn(p, q):
    if q == 4:
        return 0 if p == 2 else (1 if p % 4 == 3 else -1)
    return 0 if p == 3 else (1 if p % 3 == 2 else -1)

def lower(a):   # certified lower / upper bounds of an arb, as arb (exact endpoints)
    return a.lower()
def upper(a):
    return a.upper()

def part_A(q):
    t0 = time.time()
    Aplus = fmpz(1); Aminus = fmpz(1); Pplus = fmpz(1); Pminus = fmpz(1)
    prev = None; n = 0; bad = []
    for p in ALLP:
        if p > NA: break
        s = sgn(p, q)
        if s == 0: continue
        if prev is not None:
            d = p - prev
            Pplus *= Aplus ** d; Pminus *= Aminus ** d
            n += 1
            if not (Pplus > Pminus): bad.append(p)
        if s > 0: Aplus *= p
        else: Aminus *= p
        prev = p
    print('A. q=%d: exact integer test of Theta1(p) > 0 at the %d race primes p0 < p <= %d: failures %s   [largest integers ~ %d bits; %.1f s]'
          % (q, n, NA, bad if bad else 'none', Pplus.bit_length(), time.time() - t0))
    return not bad

def part_BCD(q):
    t0 = time.time()
    ell = 2 if q == 4 else 3
    c1 = arb(C1[q]); half = arb(1) / 2
    Th = arb(0); W = arb(0); theta = arb(0)
    D0 = 0; D1 = 0
    okB = okC = okD = True
    minT1 = None; minT1_at = 0; minG = None; minG_at = 0
    supc = None; supc_at = 0; supc2 = None; supc2_at = 0; supup = None; supup_at = 0
    minN = None; minN_at = 0
    first = True; prevp = None; Rprev = None
    nrace = 0; ngaps = 0
    startN = 29 if q == 4 else 11
    for p in ALLP:
        x = arb(p); lx = x.log(); sx = x.sqrt()
        # --- C: left limit at p
        vo = theta - (arb(ell).log() if p > ell else arb(0))
        v = (x - vo) / sx
        if not (v <= c1): okC = False
        if supc is None or upper(v) > upper(supc): supc = v; supc_at = p
        if p >= 1427 and (supc2 is None or upper(v) > upper(supc2)): supc2 = v; supc2_at = p
        theta += lx
        vo = theta - (arb(ell).log() if p >= ell else arb(0))
        u = (vo - x) / sx
        if not (u < 0): okC = False
        if supup is None or upper(u) > upper(supup): supup = u; supup_at = p
        s = sgn(p, q)
        if s == 0: continue
        # --- B
        Th += s * lx; W += s * x * lx
        D0 += s; D1 += s * p
        T1 = x * Th - W
        R = Fraction(D0) - Fraction(D1, p)
        if not first:
            nrace += 1
            if not (T1 > 0): okB = False
            if minT1 is None or lower(T1) < lower(minT1): minT1 = T1; minT1_at = p
            G = T1 / (x * sx)
            if minG is None or lower(G) < lower(minG): minG = G; minG_at = p
        # --- D
        if prevp is not None and prevp >= startN:
            r = min(R, Rprev)
            val = arb(r.numerator) / arb(r.denominator) * lx / sx
            ngaps += 1
            if not (val > half): okD = False
            if minN is None or lower(val) < lower(minN): minN = val; minN_at = prevp
        first = False; prevp = p; Rprev = R
    print('B. q=%d: Theta1(p) > 0 certified (Arb, %d bits) at all %d race primes p0 < p <= 2^22: %s' % (q, ctx.prec, nrace, okB))
    print('      min Theta1(p) = %s at p = %d ;  min Theta1(p)/p^1.5 = %s at p = %d' % (minT1.str(15), minT1_at, minG.str(15), minG_at))
    print('C. q=%d: Lemma 4(a) on (0, 2^22] with c1 = %s certified: %s' % (q, C1[q], okC))
    print('      sup (y - vo(y))/sqrt y = %s as y -> %d^- ;  sup over left limits at primes >= 1427: %s (y -> %d^-)' % (supc.str(15), supc_at, supc2.str(15), supc2_at))
    print('      sup (vo(y) - y)/sqrt y = %s at y = %d  (negative: vo(y) < y)' % (supup.str(15), supup_at))
    # --- D: the hand intervals
    hand = []
    if q == 4:
        # on [23,29): D0 = 2, D1 = 28, R = 2 - 28/x  (checked from the state below)
        def R4(xv): return Fraction(2) - Fraction(28, xv)
        for (a, b) in ((23, 25), (25, 29)):
            val = arb(R4(a).numerator) / arb(R4(a).denominator) * arb(b).log() / arb(b).sqrt()
            hand.append(('[%d,%d]' % (a, b), val, bool(val > half)))
    else:
        # on [5,7): D0 = 2, D1 = 7, R = 2 - 7/x ;  on [7,11): D0 = 1, D1 = 0, R = 1
        r6 = Fraction(2) - Fraction(7, 6)
        val = arb(r6.numerator) / arb(r6.denominator) * arb(6).log() / arb(6).sqrt()      # log x/sqrt x increasing on [6,7]
        hand.append(('[6,7]', val, bool(val > half)))
        val = arb(1) * min(arb(7).log() / arb(7).sqrt(), arb(11).log() / arb(11).sqrt(), key=lambda z: float(z.mid()))  # unimodal, max at e^2
        hand.append(('[7,11]', val, bool(val > half)))
    okD = okD and all(h[2] for h in hand)
    print('D. q=%d: N(x) > 1/2 certified on %s and on all %d gaps [p,p\'] with %d <= p < p\' <= 2^22: %s'
          % (q, ', '.join('%s (bound %s)' % (h[0], h[1].str(8)) for h in hand), ngaps, startN, okD))
    print('      smallest gap bound min(R(p),R(p\')) log p\'/sqrt p\' = %s (gap starting at p = %d)   [%.1f s]' % (minN.str(15), minN_at, time.time() - t0))
    return okB and okC and okD

def part_E(q):
    import numpy as np
    t0 = time.time()
    X = XE[q]
    # odd-only numpy sieve
    n = (X + 1) // 2
    sv = np.ones(n, dtype=bool); sv[0] = False            # index i <-> 2i+1
    for i in range(1, int(X ** 0.5) // 2 + 1):
        if sv[i]:
            pp = 2 * i + 1
            sv[(pp * pp) // 2::pp] = False
    odd = 2 * np.nonzero(sv)[0] + 1
    del sv
    old = ctx.prec; ctx.prec = 64
    half = arb(1) / 2
    startN = 29 if q == 4 else 11
    D0 = 0; D1 = 0; prevp = None; Rprev = None
    ok = True; ngaps = 0; minmid = 1e9; minat = 0
    def handle(p, s):
        nonlocal D0, D1, prevp, Rprev, ok, ngaps, minmid, minat
        D0 += s; D1 += s * p
        x = arb(p)
        R = arb(D0) - arb(D1) / x
        if prevp is not None and prevp >= startN:
            w = x.log() / x.sqrt()
            v1 = Rprev * w; v2 = R * w
            if not ((v1 > half) and (v2 > half)): ok = False
            ngaps += 1
            m = min(float(v1.mid()), float(v2.mid()))
            if m < minmid: minmid = m; minat = prevp
        prevp = p; Rprev = R
    if q == 3: handle(2, 1)
    CH = 1 << 20
    for a in range(0, len(odd), CH):
        for p in odd[a:a + CH].tolist():
            if q == 4:
                handle(p, 1 if (p & 3) == 3 else -1)
            else:
                if p == 3: continue
                handle(p, 1 if p % 3 == 2 else -1)
    ctx.prec = old
    print('E. q=%d: N(x) > 1/2 certified (Arb, 64 bits; D0, D1 exact) on all %d gaps [p,p\'] with %d <= p < p\' <= %d: %s'
          % (q, ngaps, startN, X, ok))
    print('      smallest gap bound ~ %.9f (gap starting at p = %d)   [%.1f s]' % (minmid, minat, time.time() - t0))
    return ok

if __name__ == '__main__':
    qs = [int(a) for a in sys.argv[1:]] or [4, 3]
    # sanity of the hard-wired states used in part D
    assert sum(sgn(p, 4) for p in ALLP if p <= 23) == 2 and sum(sgn(p, 4) * p for p in ALLP if p <= 23) == 28
    assert sum(sgn(p, 3) for p in ALLP if p <= 5) == 2 and sum(sgn(p, 3) * p for p in ALLP if p <= 5) == 7
    assert sum(sgn(p, 3) for p in ALLP if p <= 7) == 1 and sum(sgn(p, 3) * p for p in ALLP if p <= 7) == 0
    allok = True
    for q in qs:
        allok &= part_A(q)
        allok &= part_BCD(q)
        allok &= part_E(q)
        sys.stdout.flush()
    print('ALL CERTIFIED:', allok)
