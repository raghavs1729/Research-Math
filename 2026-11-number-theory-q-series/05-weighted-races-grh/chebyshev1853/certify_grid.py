#!/usr/bin/env python3
"""
certify_grid.py -- lower bounds for  N_C(x) = -C(x) log x / sqrt x  valid for ALL real x in a range (not only at
grid points), from Taylor expansions in h = log(x/x_j) around the grid points x_j.

  C(x_j e^h) = sum_p chi(p) phi_p(h),  phi_p(h) = exp(-z e^{-h}),  z = p/x_j.
  phi_p^{(k)}(h) = P_k(w) e^{-w},  w = z e^{-h},   P_0 = 1,  P_{k+1}(w) = w (P_k(w) - P_k'(w)).
  D_k := d^k/dh^k C(x_j e^h) at h=0  = sum_m c_{k,m} (M_m^+ - M_m^-),   M_m^{+-} = sum_{chi(p)=+-1} z^m e^{-z}.
  For |h| <= delta:
     -C(x_j e^h) >= -D_0 - sum_{k=1}^{K} |D_k| delta^k/k! - R,
     R = delta^{K+1}/(K+1)! * sum_m |c_{K+1,m}| e^{2 m delta} [ X m! + (m/e)^m ],   X = x_j e^delta
  (sum over ALL integers n >= 1 of v^m e^{-v}, v = n/X, is <= X m! + (m/e)^m for the unimodal function v^m e^{-v};
   and sup_{|h|<=delta} w^m e^{-w} <= e^{2 m delta} v^m e^{-v} with v = z e^{-delta} = p/X).
  On the cell, log x/sqrt x >= its value at the right end (cells with left end >= e^2).

Mode "scan":  uses the long double moments printed by abel_scan (K = 8, CUT = 50), with an a-priori error budget:
     |computed M_m - true M_m| <= ETA * (M_m^+ + M_m^-) + tail_m,
     ETA = nprimes * 2^-64 + 2^-45 + 2^-47
           (long double accumulation of nprimes nonnegative terms, each term carrying a relative error <= 2^-45
            from exp() and from the double z; 2^-47 for reading the printed values into doubles and for the double
            arithmetic of this script; see the note, section 5),
     tail_m = x Gamma(m+1, 50) + 50^m e^{-50} (m <= 50)   (primes p > 50 x, not summed by the C program).
Mode "arb":   recomputes the moments in Arb ball arithmetic on a fine grid; the prime tails, the Taylor remainder R
              and every comparison are also evaluated in Arb (certified, no floating-point decision).

Usage: python3 certify_grid.py scan out_q4_1e8.txt [target]
       python3 certify_grid.py arb q lo hi M [target]         grid 10^(j/M), lo <= x <= hi
"""
import sys, math, os
from fractions import Fraction

K = int(os.environ.get('TAYLOR_K', '8'))      # Taylor order (the scan files contain moments up to 8)


def stirling_polys(kmax):
    """coefficient lists of P_0..P_kmax (index = power of w)"""
    P = [[1]]
    for k in range(kmax):
        p = P[-1]
        d = [i * p[i] for i in range(1, len(p))] + [0]        # P'
        diff = [p[i] - (d[i] if i < len(d) else 0) for i in range(len(p))]
        P.append([0] + diff)
    return P


POL = stirling_polys(K + 1)


def remainder_bound(X, delta):
    tot = 0.0
    for m, c in enumerate(POL[K + 1]):
        if c:
            tot += abs(c) * math.exp(2 * m * delta) * (X * math.factorial(m) + (m / math.e) ** m if m else X + 1)
    return tot * delta ** (K + 1) / math.factorial(K + 1)


def upper_gamma(m, a):
    """Gamma(m+1, a) for integer m: m! e^{-a} sum_{i<=m} a^i/i!"""
    return math.factorial(m) * math.exp(-a) * sum(a ** i / math.factorial(i) for i in range(m + 1))


def scan_mode(fn, target):
    rows = []
    nprimes = None
    Mgrid = None
    for line in open(fn):
        if line.startswith('# abel_scan'):
            for tok in line.split():
                if tok.startswith('nprimes'):
                    nprimes = int(tok.split('=')[1])
                if tok.startswith('M='):
                    Mgrid = int(tok.split('=')[1])
            continue
        if line[0] == '#':
            continue
        f = line.split()
        j = int(f[0]); x = float(f[1])
        mom = [float(v) for v in f[5:]]
        Mp = mom[0::2]; Mm = mom[1::2]
        rows.append((j, x, Mp, Mm))
    delta = math.log(10) / (2 * Mgrid)
    ETA = nprimes * 2.0 ** -64 + 2.0 ** -45 + 2.0 ** -47     # last term: decimal->double conversion and the double arithmetic of this script
    print("# file %s: M = %d (delta = %.6f), nprimes = %d, ETA = %.3e, K = %d" % (fn, Mgrid, delta, nprimes, ETA, K))
    res = []
    worst_err = 0.0
    for j, x, Mp, Mm in rows:
        # derivatives and their error budgets
        D = []; E = []
        for k in range(K + 1):
            d = 0.0; e = 0.0
            for m, c in enumerate(POL[k]):
                if c:
                    d += c * (Mp[m] - Mm[m])
                    e += abs(c) * (ETA * (Mp[m] + Mm[m]) + x * upper_gamma(m, 50.0) + 50.0 ** m * math.exp(-50.0))
            D.append(d); E.append(e * (1 + 1e-12) + abs(d) * 1e-14)
        X = x * math.exp(delta)
        R = remainder_bound(X, delta)
        low = -D[0] - E[0]
        for k in range(1, K + 1):
            low -= (abs(D[k]) + E[k]) * delta ** k / math.factorial(k)
        low -= R
        low *= (1 - 1e-13)
        up = -D[0] + E[0] + sum((abs(D[k]) + E[k]) * delta ** k / math.factorial(k) for k in range(1, K + 1)) + R
        xl = x * math.exp(-delta)
        worst_err = max(worst_err, (E[0] + R) / abs(D[0]) if (D[0] and x >= 1) else 0)
        if xl >= math.e ** 2 and low > 0:
            nlow = low * math.log(X) / math.sqrt(X) * (1 - 1e-13)
        else:
            nlow = float('nan')
        res.append((j, x, xl, X, -D[0], low, up, nlow))
    return res, worst_err


def report(res, target, label):
    neg = [r for r in res if not (r[5] > 0)]
    pos_from = None
    for r in res:
        if not (r[5] > 0):
            pos_from = None
        elif pos_from is None:
            pos_from = r[2]
    print("%s: cells with certified -C > 0 on the whole cell: %d of %d ; all cells from x >= %.6g on (cells not certified: %d, all with x <= %.4g)"
          % (label, len(res) - len(neg), len(res), pos_from, len(neg), max([r[3] for r in neg]) if neg else 0))
    ok_from = None
    for r in res:
        if not (r[7] >= target):
            ok_from = None
        elif ok_from is None:
            ok_from = r[2]
    print("%s: N_C(x) >= %.4g certified for every x in [%.8g, %.8g]" % (label, target, ok_from, res[-1][3]))
    # minima of the certified lower bound per decade
    dec = {}
    for r in res:
        if r[7] == r[7]:
            d = math.floor(math.log10(r[1]) + 1e-9)
            if d not in dec or r[7] < dec[d][0]:
                dec[d] = (r[7], r[1], -r[4] and r[4] * math.log(r[1]) / math.sqrt(r[1]))
    for d in sorted(dec):
        print("   decade 10^%d: min certified lower bound of N_C over cells = %.6f (cell at x_j = %.6g, N_C(x_j) = %.6f)" % (d, dec[d][0], dec[d][1], dec[d][2]))


def arb_mode(q, lo, hi, M, target):
    """Everything that enters a decision is an Arb ball (80 bits), including the Taylor remainder R, the prime
    tails and the comparison with the target (revised after the referee report, item REP-2: these were IEEE
    doubles inflated by 1e-6 before).  Floats are used only (a) to choose the integer cut-off P_j ~ 60 x_j up to
    which primes are summed -- the tail bound is then computed in Arb from the ball c_j = P_j/x_j, whatever P_j
    is -- and (b) to print rounded decimals."""
    from flint import arb, ctx
    import numpy as np
    ctx.prec = 80
    CUT = 60
    tgt = arb(target)                      # target is passed as a decimal string
    jlo = math.ceil(math.log10(lo) * M - 1e-9); jhi = math.floor(math.log10(hi) * M + 1e-9)
    xs = [arb(10) ** (arb(j) / M) for j in range(jlo, jhi + 1)]
    P = int(CUT * 10 ** (jhi / M)) + 2
    bs = np.ones(P + 1, dtype=np.bool_); bs[:2] = False
    i = 2
    while i * i <= P:
        if bs[i]:
            bs[i * i::i] = False
        i += 1
    primes = np.nonzero(bs)[0]
    del bs
    delta = arb(10).log() / (2 * M)
    edelta = delta.exp()
    e2 = arb(2).exp()
    fact = [arb(math.factorial(m)) for m in range(K + 3)]
    # remainder: R = delta^(K+1)/(K+1)! * sum_m |c_{K+1,m}| e^{2 m delta} [X m! + (m/e)^m]
    rcoef = []
    for m, c in enumerate(POL[K + 1]):
        if c:
            mm = (arb(m) / arb(1).exp()) ** m if m else arb(1)
            rcoef.append((abs(c) * (2 * m * delta).exp(), fact[m], mm))
    rfac = delta ** (K + 1) / fact[K + 1]
    res = []
    for jj, x in enumerate(xs):
        xf = float(x.mid())
        Pj = int(CUT * xf)                 # primes p <= Pj are summed; the integers n > Pj are covered by the tail bound
        Mp = [arb(0) for _ in range(K + 1)]; Mm = [arb(0) for _ in range(K + 1)]
        for p in primes[:int(np.searchsorted(primes, Pj, side='right'))].tolist():
            if q == 4:
                if p == 2: continue
                s = 1 if p % 4 == 1 else -1
            else:
                if p == 3: continue
                s = 1 if p % 3 == 1 else -1
            z = p / x
            t = (-z).exp()
            S = Mp if s > 0 else Mm
            S[0] += t
            for k in range(1, K + 1):
                t = t * z
                S[k] += t
        # tails: f_m(v) = v^m e^{-v} is decreasing for v >= m, and every n > Pj has n/x > c := Pj/x, so
        #   sum_{n > Pj} f_m(n/x) <= f_m(c) + x int_c^oo f_m = c^m e^{-c} + x Gamma(m+1, c),
        #   Gamma(m+1, c) = m! e^{-c} sum_{i<=m} c^i/i!      (requires c >= m; checked in ball arithmetic)
        c = arb(Pj) / x
        if not (c > K + 2):
            raise RuntimeError("cut-off too small")
        ec = (-c).exp()
        tailm = []
        for m in range(K + 1):
            G = fact[m] * ec * sum((c ** i / fact[i] for i in range(m + 1)), arb(0))
            tailm.append(c ** m * ec + x * G)
        low = arb(0)
        D0 = None
        for k in range(K + 1):
            d = arb(0); tail = arb(0)
            for m, cf in enumerate(POL[k]):
                if cf:
                    d += cf * (Mp[m] - Mm[m])
                    tail += abs(cf) * tailm[m]
            # |true D_k - d| <= tail
            if k == 0:
                D0 = d
                low = -d - tail
            else:
                low -= (abs(d) + tail) * delta ** k / fact[k]
        X = x * edelta
        R = rfac * sum((a * (X * f + mm) for a, f, mm in rcoef), arb(0))
        low -= R
        pos = bool(low > 0)
        left = x / edelta
        left_ok = bool(left > e2)          # log x/sqrt x is decreasing for x >= e^2: its minimum on the cell is at the right end
        nl = low * X.log() / X.sqrt()
        ok = pos and left_ok and bool(nl > tgt)
        res.append(dict(j=jlo + jj, x=xf, left=float(left.mid()), right=float(X.mid()), negC=float((-D0).mid()),
                        low=float(low.lower()), pos=pos, nlow=float(nl.lower()) if (pos and left_ok) else float('nan'), ok=ok,
                        relR=float((R / abs(D0)).upper()) if D0 != 0 else 0.0, reltail=float((tailm[0] / abs(D0)).upper()) if D0 != 0 else 0.0))
    return res


def report_arb(res, target, M):
    # decisions (pos, ok) were taken in ball arithmetic; the decimals printed here are rounded for display,
    # interval ends are rounded inwards (left end up, right end down) at 6 decimals
    up6 = lambda v: math.ceil(v * 1e6 + 1e-3) / 1e6
    dn6 = lambda v: math.floor(v * 1e6 - 1e-3) / 1e6
    npos = sum(1 for r in res if r['pos'])
    print("arb: cells [x_j e^-delta, x_j e^delta] with -C > 0 certified on the whole cell: %d of %d" % (npos, len(res)))
    first = None
    for r in res:
        if not r['ok']:
            first = None
        elif first is None:
            first = r
    if first is None:
        print("arb: N_C(x) > %s not certified on the last cell" % target)
    else:
        print("arb: N_C(x) > %s certified (ball arithmetic) for every x in [%.6f, %.6f]   (%d consecutive cells, from j = %d)"
              % (target, up6(first['left']), dn6(res[-1]['right']), sum(1 for r in res if r['j'] >= first['j']), first['j']))
    print("arb: largest (Taylor remainder R)/|C(x_j)| = %.3e ;  largest (prime tail of M_0)/|C(x_j)| = %.3e"
          % (max(r['relR'] for r in res), max(r['reltail'] for r in res)))
    dec = {}
    for r in res:
        if r['nlow'] == r['nlow']:
            d = math.floor(math.log10(r['x']) + 1e-9)
            if d not in dec or r['nlow'] < dec[d][0]:
                dec[d] = (r['nlow'], r['x'], r['negC'] * math.log(r['x']) / math.sqrt(r['x']))
    for d in sorted(dec):
        print("   decade 10^%d: min lower bound of N_C over cells = %.6f (cell at x_j = %.6g, N_C(x_j) = %.6f)" % (d, dec[d][0], dec[d][1], dec[d][2]))
    tf = float(target)
    for r in res:
        if (r['nlow'] == r['nlow'] and abs(r['nlow'] - tf) < 0.004) or r['j'] % max(1, M // 2) == 0:
            print("   x_j = %.6f  cell [%.6f, %.6f]  -C(x_j) = %.12f  N_C(x_j) = %.8f  lower bound of N_C on the cell = %.8f  certified > %s: %s"
                  % (r['x'], r['left'], r['right'], r['negC'], r['negC'] * math.log(r['x']) / math.sqrt(r['x']), r['nlow'], target, r['ok']))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'scan':
        target = float(sys.argv[3]) if len(sys.argv) > 3 else 0.88
        res, werr = scan_mode(sys.argv[2], target)
        print("# largest relative error budget (rounding + truncation + Taylor remainder) of -C over the grid points x_j >= 1: %.3e" % werr)
        report(res, target, "scan")
        print("# selected cells:  x_j   -C(x_j)   certified lower bound of -C on cell   certified lower bound of N_C on cell")
        for r in res:
            if r[0] % 50 == 0 and r[1] >= 1:
                print("   %12.6g %18.10f %18.10f %10.6f" % (r[1], r[4], r[5], r[7]))
    else:
        q = int(sys.argv[2]); lo = float(sys.argv[3]); hi = float(sys.argv[4]); M = int(sys.argv[5])
        target = sys.argv[6] if len(sys.argv) > 6 else '0.88'
        res = arb_mode(q, lo, hi, M, target)
        print("# arb mode: q = %d, grid 10^(j/%d) in [%g, %g], Taylor order K = %d, 80-bit balls, primes up to about 60 x;" % (q, M, lo, hi, K))
        print("# moments, prime tails, Taylor remainder and all comparisons in Arb ball arithmetic (no floating-point decision)")
        report_arb(res, target, M)
