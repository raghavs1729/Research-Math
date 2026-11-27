#!/usr/bin/env python3
"""
cert_gaps.py -- part E of cert_small.py carried to a large limit X with a segmented sieve:
certified (Arb ball arithmetic, exact integers D0, D1) proof that

        N(x) = R(x) log x / sqrt x  >  1/2      for every real x in [x_start, X],

x_start = 29 for q = 4 and 11 for q = 3 (the short ranges below are in cert_small.py, part D).

Method.  Between consecutive race primes p < p', R(x) = D0 - D1/x (D0, D1 the exact integer states
after p) is monotone and continuous, and w(x) = log x/sqrt x is decreasing (x > e^2).  Hence

        N(x) >= min(R(p), R(p')) * w(p')      on [p, p'].

For every such gap with p >= x_start the script checks, with Arb balls at 64 bits (so the comparison
is certain or the script reports a failure):   R(p) w(p') > 1/2   and   R(p') w(p') > 1/2,
where R(p) = D0 - D1/p is computed as a ball from the exact integers.  The last piece
[p_last, X] is treated in the same way with p' replaced by X (no prime at X is needed).

Nothing here is a floating-point decision; floats are used only to locate the minimum for the
report, and the reported minimum is then given as a certified ball.

Usage:  nohup nice -n 10 python3 cert_gaps.py q X > cert_gaps_q{q}.out &
        (X = 1e10: about 4.55e8 gaps; checkpoints are printed every 5e8)
"""
import sys, time, math
import numpy as np
from flint import arb, ctx

EXPECT = {  # exact states at 10^10 printed by the independent C program cesaro_scan (out_q*_1e10.txt)
    (4, 10 ** 10): (5960, 27201020910942),
    (3, 10 ** 10): (6264, 27948685258617),
}

def odd_prime_segments(X, seg_odd=1 << 24):
    """yield numpy arrays of the odd primes in consecutive segments of [3, X]"""
    r = math.isqrt(X) + 1
    sv = np.ones(r + 1, dtype=bool); sv[:2] = False
    for i in range(2, math.isqrt(r) + 1):
        if sv[i]: sv[i * i::i] = False
    base = [int(c) for c in np.nonzero(sv)[0] if c > 2]
    lo = 3
    while lo <= X:
        hi = min(lo + 2 * (seg_odd - 1), X)
        n = (hi - lo) // 2 + 1
        seg = np.ones(n, dtype=bool)
        for c in base:
            c2 = c * c
            if c2 > hi: break
            start = c2 if c2 >= lo else ((lo + c - 1) // c) * c
            if start % 2 == 0: start += c
            seg[(start - lo) // 2::c] = False
        yield lo + 2 * np.nonzero(seg)[0]
        lo += 2 * n

def main(q, X):
    t0 = time.time()
    ctx.prec = 64
    half = arb(1) / 2
    startN = 29 if q == 4 else 11
    D0 = 0; D1 = 0; prevp = 0; Rprev = None
    ok = True; ngaps = 0; nfail = 0
    minmid = 1e9; minat = 0; minv = None
    print('# cert_gaps q=%d X=%d  (Arb, %d bits; D0, D1 exact Python integers; segmented sieve)' % (q, X, ctx.prec)); sys.stdout.flush()
    next_chk = 5 * 10 ** 8
    if q == 3:      # the race prime 2
        D0 = 1; D1 = 2; prevp = 2; Rprev = arb(D0) - arb(D1) / arb(2)
    for parr in odd_prime_segments(X):
        if q == 4:
            sarr = np.where((parr & 3) == 3, 1, -1)
        else:
            keep = parr != 3
            parr = parr[keep]
            sarr = np.where(parr % 3 == 2, 1, -1)
        for p, s in zip(parr.tolist(), sarr.tolist()):
            D0 += s; D1 += s * p
            x = arb(p)
            R = arb(D0) - arb(D1) / x
            if prevp >= startN:
                w = x.log() / x.sqrt()
                v1 = Rprev * w; v2 = R * w
                if not ((v1 > half) and (v2 > half)):
                    ok = False; nfail += 1
                ngaps += 1
                m1 = float(v1.mid()); m2 = float(v2.mid())
                if m1 < minmid: minmid = m1; minat = prevp; minv = v1
                if m2 < minmid: minmid = m2; minat = prevp; minv = v2
            prevp = p; Rprev = R
        # checkpoint when the segment end passes a multiple of 5e8
        if len(parr) and int(parr[-1]) >= next_chk:
            seg_end = int(parr[-1])
            print('CHK  primes <= %d done: gaps %d, failures %d, smallest gap bound so far %s (gap starting at %d)  [%.0f s]'
                  % (seg_end, ngaps, nfail, minv.str(12), minat, time.time() - t0)); sys.stdout.flush()
            while next_chk <= seg_end: next_chk += 5 * 10 ** 8
    # the last piece [p_last, X]
    x = arb(X)
    RX = arb(D0) - arb(D1) / x
    w = x.log() / x.sqrt()
    v1 = Rprev * w; v2 = RX * w
    tail_ok = bool((v1 > half) and (v2 > half))
    ok = ok and tail_ok
    print('TAIL [%d, %d]: bounds %s, %s : %s' % (prevp, X, v1.str(12), v2.str(12), tail_ok))
    print('STATE at X: D0 = %d, D1 = %d' % (D0, D1), end='')
    if (q, X) in EXPECT:
        print('   (cesaro_scan: D0 = %d, D1 = %d : %s)' % (EXPECT[(q, X)] + (EXPECT[(q, X)] == (D0, D1),)))
        ok = ok and EXPECT[(q, X)] == (D0, D1)
    else:
        print()
    lowm = (minv - half).lower()
    print('RESULT q=%d: N(x) > 1/2 certified for every real x in [%d, %d]: %s' % (q, startN, X, ok))
    print('       gaps checked: %d (failures %d);  smallest gap bound = %s (gap starting at p = %d)' % (ngaps, nfail, minv.str(15), minat))
    print('       smallest certified margin over 1/2: >= %s' % lowm.str(12))
    print('       elapsed %.0f s' % (time.time() - t0))
    sys.stdout.flush()

if __name__ == '__main__':
    q = int(sys.argv[1]); X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 10
    assert q in (3, 4)
    main(q, X)
