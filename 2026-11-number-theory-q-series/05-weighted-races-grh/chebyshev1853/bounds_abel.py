#!/usr/bin/env python3
"""
bounds_abel.py -- the closed-form lower bounds of CHEBYSHEV1853.md, their thresholds, the comparison with the
scan (abel_scan output), and the grid argument for the uniform constant.

Usage: python3 bounds_abel.py [q]         (q = 4 default, or 3)

Everything here is ordinary mpmath arithmetic at 30 digits on closed-form expressions (no interval arithmetic);
the constants S_Gamma, c1, c3 are the certified upper bounds from cert_abel.out, rounded UP.

Notation (q = 4:  ell = 2, p0 = 3;   q = 3:  ell = 3, p0 = 2):
  W(y)  = sum_p s(p) log p e^{-p/y},     s(p) = -chi(p)
  C(x)  = sum_p chi(p) e^{-p/x},         N_C(x) = -C(x) log x / sqrt x
  kappa(y):   W(y) >= kappa(y) sqrt(y)   for y >= 4          (Theorem 3 of the note, GRH(chi))
  kflat(y):   a minorant of kappa with  inf_{y' >= y} kappa(y') >= kstar(y) := min(kflat(y), later regime values)
"""
import sys, os, math
import mpmath as mp
mp.mp.dps = 30

q = int(sys.argv[1]) if len(sys.argv) > 1 else 4
HERE = os.path.dirname(os.path.abspath(__file__))

# ---- certified constants (cert_abel.out), rounded up -------------------------------------------------
if q == 4:
    SG = mp.mpf('0.00039199328')      # S_Gamma = sum_rho |Gamma(rho)|
    c1 = mp.mpf('2.071193')
    c3 = mp.mpf('2.153716634')
    p0 = 3
else:
    SG = mp.mpf('0.000016531314')
    c1 = mp.mpf('2.081942')
    c3 = mp.mpf('2.074662259')
    p0 = 2
U1 = mp.mpf(2) ** 22
eps_a = mp.mpf('0.2') / mp.log(U1) ** 2            # BKLNW Cor. 11.2, row k=2 (x >= 3 594 641)
cth = 1 + eps_a
# ladder for eps(u), u > U1:  0.2/log^2 on [2^22, 89967803),  1/log^3 U on [U, next) for U >= 89967803 (row k=3)
LAD = [(U1, eps_a)] + [(mp.mpf(U), 1 / mp.log(U) ** 3) for U in
                       [89967803, 10 ** 10, 10 ** 12, 10 ** 15, 10 ** 20, 10 ** 30, 10 ** 50, 10 ** 100, 10 ** 200, 10 ** 400]]
HALFSQRTPI = mp.sqrt(mp.pi) / 2
G54 = mp.gamma(mp.mpf(5) / 4)
G98 = mp.gamma(mp.mpf(9) / 8)


def F(v, up=True):
    """F(v) = v e^{-v^2} + (sqrt(pi)/2) erfc(v), decreasing from sqrt(pi)/2 to 0.  For v > 100 (F < 1e-4000) the
    value is replaced by the upper bound 1e-1000 (up=True) or the lower bound 0 (up=False)."""
    if v > 100:
        return mp.mpf(10) ** -1000 if up else mp.mpf(0)
    return v * mp.exp(-v * v) + HALFSQRTPI * mp.erfc(v)


def Etail_over_sqrt(y):
    """E_tail(y)/sqrt(y) = sum_i eps_i [F(U_i/sqrt y) - F(U_{i+1}/sqrt y)]."""
    r = mp.sqrt(y)
    tot = mp.mpf(0)
    for i, (U, e) in enumerate(LAD):
        a = F(U / r)
        b = F(LAD[i + 1][0] / r, up=False) if i + 1 < len(LAD) else mp.mpf(0)
        tot += e * (a - b)
    return tot


def tau(k, y):
    return U1 * mp.exp(-U1 ** k / y) * (1 + y / (k * U1 ** k))


def Bk(k, y):
    """bound for |W_k(y)|, k odd >= 3"""
    triv = cth * mp.gamma(1 + mp.mpf(1) / k) * y ** (mp.mpf(1) / k)
    sv = c3 * mp.gamma(1 + mp.mpf(1) / (2 * k)) * y ** (mp.mpf(1) / (2 * k)) + cth * tau(k, y)
    return min(triv, sv)


def kappa(y):
    """W(y) >= kappa(y) sqrt(y), y >= 4 (exact form of Theorem 3)."""
    y = mp.mpf(y)
    r = mp.sqrt(y)
    val = (HALFSQRTPI - SG) * r - Etail_over_sqrt(y) * r - c1 * G54 * y ** mp.mpf('0.25')
    val += max(0, (1 - eps_a) * G54 * y ** mp.mpf('0.25') - c1 * G98 * y ** mp.mpf('0.125'))
    k = 3
    while mp.mpf(p0) ** k < y:
        val -= Bk(k, y)
        k += 2
    val -= 1
    return val / r


def kappa_upper(y):
    """W(y) <= kappa_upper(y) sqrt(y), y >= 4."""
    y = mp.mpf(y)
    r = mp.sqrt(y)
    b = mp.mpf('0.7831887855') if q == 4 else mp.mpf('0.9481988267')
    val = (HALFSQRTPI + SG) * r + Etail_over_sqrt(y) * r + b + mp.mpf('0.2555') / y ** 2
    k = 3
    while mp.mpf(p0) ** k < y:
        val += Bk(k, y)
        k += 2
    val += 1
    # even k >= 4:  Q_k <= c_theta Gamma(1+1/k) y^{1/k} for p0^k < y, and the rest <= 1 (same argument as T_infinity)
    k = 4
    while mp.mpf(p0) ** k < y:
        val += cth * mp.gamma(1 + mp.mpf(1) / k) * y ** (mp.mpf(1) / k)
        k += 2
    val += 1
    return val / r


# ---- monotone minorant -------------------------------------------------------------------------------
Y_E = mp.mpf(10) ** 12       # below: E_tail/sqrt y <= EB_SMALL ; above: <= eps_a sqrt(pi)/2
Y_I = mp.mpf(10) ** 16       # below: sieve bound (c3) for all odd k ; above: trivial bound for all odd k
EB_SMALL = eps_a * F(U1 / mp.sqrt(Y_E))
EB_LARGE = eps_a * HALFSQRTPI
TAU_SMALL = cth * tau(3, Y_I)      # >= c_theta tau_k(y) for all odd k >= 3 and y <= Y_I (astronomically small)


def kflat(y, regime=None):
    """Minorant with the property:  kappa(y') >= kflat(y)  for all y' >= y in the same regime
    (regimes: [500, Y_E], (Y_E, Y_I], (Y_I, oo)).  Reason: every term of kappa below is monotone in the right
    direction for fixed k, and the odd powers k that are absent at y (p0^k >= y) but present at y' > p0^k
    contribute at y' at most  c3 (p0^k)^{1/(2k)-1/2} = c3 p0^{1/2} p0^{-k/2}  (regime I, II; sieve bound)
    resp.  c_theta (p0^k)^{1/k-1/2} = c_theta p0 p0^{-k/2}  (regime III; trivial bound), whose sum over the odd
    k with p0^k >= y is at most  c3 p0^{1/2} y^{-1/2}/(1-1/p0)  resp.  c_theta p0 y^{-1/2}/(1-1/p0)."""
    y = mp.mpf(y)
    r = mp.sqrt(y)
    if regime is None:
        regime = 1 if y <= Y_E else (2 if y <= Y_I else 3)
    val = HALFSQRTPI - SG
    val -= EB_SMALL if regime == 1 else EB_LARGE
    # squares and fourth powers (valid form for y >= 463: the max(0,.) is active)
    val -= (c1 - 1 + eps_a) * G54 * y ** mp.mpf('-0.25') + c1 * G98 * y ** mp.mpf('-0.375')
    k = 3
    tot = mp.mpf(0)
    if regime <= 2:
        while mp.mpf(p0) ** k < y:
            tot += c3 * mp.gamma(1 + mp.mpf(1) / (2 * k)) * y ** (mp.mpf(1) / (2 * k))
            k += 2
        tot += c3 * mp.sqrt(p0) / (1 - mp.mpf(1) / p0) + 1 + 40 * TAU_SMALL
    else:
        while mp.mpf(p0) ** k < y:
            tot += cth * mp.gamma(1 + mp.mpf(1) / k) * y ** (mp.mpf(1) / k)
            k += 2
        tot += cth * p0 / (1 - mp.mpf(1) / p0) + 1
    return val - tot / r


def kstar(y):
    """lower bound for inf_{y' >= y} kappa(y'), y >= 500."""
    y = mp.mpf(y)
    v = kflat(y)
    if y <= Y_E:
        v = min(v, KFE)
    if y <= Y_I:
        v = min(v, KFI)
    return v


KFE = kflat(Y_E, regime=2)      # lower bound for kappa on (Y_E, Y_I]
KFI = kflat(Y_I, regime=3)      # lower bound for kappa on (Y_I, oo)


# ---- the bound for N_C ------------------------------------------------------------------------------
CA = 2 * mp.log((1 + mp.sqrt(2)) / 2)


def A(V):
    s = mp.sqrt(1 + mp.mpf(V))
    return -mp.log((s + 1) / (s - 1))


VS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]


def NC_bound(x, kfun=kappa, x_hi=None, vs=VS):
    """Theorem 4: N_C(x') >= value, for x' = x (x_hi None) or for every x' in [x, x_hi] (then kfun must be a
    nondecreasing-safe minorant such as kstar).  Pieces v in (0,1], (V_{i-1}, V_i]."""
    x = mp.mpf(x)
    xh = mp.mpf(x_hi) if x_hi else x
    k0 = max(mp.mpf(0), kfun(x / 2))
    val = k0 * (1 - CA / mp.log(x))
    kmax = k0
    Vn = 1
    for i in range(1, len(vs)):
        Vp, Vc = vs[i - 1], vs[i]
        ymin = x / (1 + Vc)
        if ymin < 500:
            break
        ki = max(mp.mpf(0), kfun(ymin))
        val += ki * (A(Vc) - A(Vp)) * mp.log(x) / mp.log(xh / Vp) ** 2
        kmax = max(kmax, ki)
        Vn = Vc
    An = mp.log(x / Vn)
    val -= HALFSQRTPI * 2 * mp.sqrt(1 + Vn) * (1 + An) / An ** 2 * mp.log(xh) / x
    return val


def load(fn):
    rows = {}
    try:
        for line in open(os.path.join(HERE, fn)):
            if line[0] == '#':
                continue
            f = line.split()
            rows[float(f[1])] = (float(f[2]), float(f[3]), float(f[4]))
    except FileNotFoundError:
        pass
    return rows


if __name__ == '__main__':
    print("q = %d   S_Gamma <= %s   c1 = %s   c3 = %s   eps_a = %s" % (q, mp.nstr(SG, 10), c1, mp.nstr(c3, 10), mp.nstr(eps_a, 8)))
    print("sqrt(pi)/2 = %s ;  sqrt(pi)/2 - S_Gamma = %s ;  EB_SMALL = %s, EB_LARGE = %s, TAU_SMALL = %s"
          % (mp.nstr(HALFSQRTPI, 12), mp.nstr(HALFSQRTPI - SG, 12), mp.nstr(EB_SMALL, 3), mp.nstr(EB_LARGE, 6), mp.nstr(TAU_SMALL, 3)))
    # 1. threshold of positivity of kappa
    ys = [mp.mpf(10) ** (mp.mpf(j) / 200) for j in range(200, 200 * 6 + 1)]
    neg = [y for y in ys if kappa(y) <= 0]
    print("\n1. kappa(y) > 0 on the grid 10^(j/200) for all y >= %s   (largest grid y with kappa <= 0: %s)"
          % (mp.nstr(ys[ys.index(neg[-1]) + 1], 6), mp.nstr(neg[-1], 6)))
    bad = [y for y in ys if y >= 500 and kflat(y) > kappa(y) + mp.mpf(10) ** -25]
    print("   kflat <= kappa on the grid (y >= 500): %s" % (not bad))
    # numerical test of the defining property of kflat on the grid: kappa(y') >= kflat(y) for y' >= y (same regime)
    sel = [y for y in ys if y >= 500]
    kap = [kappa(y) for y in sel]
    sufmin = kap[:]
    for i in range(len(sel) - 2, -1, -1):
        sufmin[i] = min(sufmin[i], sufmin[i + 1])
    print("   min_{y' >= y on the grid} kappa(y') >= kflat(y) for all grid y in [500, 10^6]: %s" % all(sufmin[i] >= kflat(sel[i]) for i in range(len(sel))))
    for yy in (500, 69905, 10 ** 6, Y_E, Y_I):
        print("   kflat(%s) = %s" % (mp.nstr(mp.mpf(yy), 8), mp.nstr(kflat(yy), 10)))
    print("   KFE = kflat(Y_E, regime 2) = %s ;  KFI = kflat(Y_I, regime 3) = %s" % (mp.nstr(KFE, 10), mp.nstr(KFI, 10)))
    # T_infinity: (1/e) sum_{p >= p0, p not| q} log p (p0/p)^3 / (1 - p^-2)  <= 1 ?
    N = 10 ** 6
    sv = bytearray([1]) * (N + 1); sv[0] = sv[1] = 0
    i = 2
    while i * i <= N:
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
        i += 1
    tsum = mp.mpf(0)
    for p in range(2, N + 1):
        if sv[p] and (q % p):
            tsum += mp.log(p) * (mp.mpf(p0) / p) ** 3 / (1 - mp.mpf(p) ** -2)
    ttail = mp.mpf(p0) ** 3 * (mp.log(N) / (2 * mp.mpf(N) ** 2) + 1 / (4 * mp.mpf(N) ** 2)) * mp.mpf(9) / 8   # sum_{n>N} log n / n^3
    print("   T_infinity constant: (1/e) * sum_p log p (p0/p)^3/(1-p^-2) <= %s  (must be <= 1)" % mp.nstr((tsum + ttail) / mp.e, 8))
    # b - r^+(4) > 0
    bq = mp.mpf('0.7831887854') if q == 4 else mp.mpf('0.9481988266')
    b0q = mp.mpf('-0.0602669141') if q == 4 else mp.mpf('0.1438316446')
    s0q = mp.mpf('0.15556798') if q == 4 else mp.mpf('0.11322997')
    rp = lambda y: bq / (2 * y * y) + (1 + b0q + mp.log(1 + y)) / y + mp.mpf('0.4') * s0q / (y * y)
    print("   b - r^+(4) = %s  (must be > 0; r^+ is decreasing)" % mp.nstr(bq - rp(mp.mpf(4)), 8))
    # 2. comparison with the scan
    rows = {}
    for fn in ('out_q%d_1e8.txt' % q, 'out_q%d_1e9.txt' % q):
        rows.update(load(fn))
    print("\n2. W(y)/sqrt(y): lower bound kappa, scan value, upper bound;  N_C(x): bound of Theorem 4 (with kstar), scan value")
    print("   %10s %10s %10s %10s %10s | %10s %10s %7s" % ("x=y", "kappa", "W/sqrt", "kappa_up", "kstar", "NC_bound", "N_C scan", "ratio"))
    for e in [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]:
        x = 10.0 ** e
        key = min(rows, key=lambda t: abs(math.log(t) - math.log(x))) if rows else None
        if key is None or abs(math.log(key / x)) > 1e-6:
            sc = (float('nan'),) * 3
        else:
            sc = rows[key]
        kap = kappa(x); ku = kappa_upper(x)
        nb = NC_bound(x, kfun=kstar) if x >= 2000 else mp.nan
        nsc = -sc[0] * math.log(x) / math.sqrt(x)
        print("   %10.4g %10.6f %10.6f %10.6f %10.6f | %10.6f %10.6f %7.3f"
              % (x, kap, sc[1] / math.sqrt(x), ku, kstar(x) if x >= 500 else float('nan'), nb, nsc, float(nb) / nsc if nsc == nsc else float('nan')))
    viol = [(x, kappa(x), r[1] / math.sqrt(x)) for x, r in rows.items() if x >= 4 and kappa(x) > r[1] / math.sqrt(x)]
    viol2 = [(x) for x, r in rows.items() if x >= 4 and kappa_upper(x) < r[1] / math.sqrt(x)]
    print("   grid points (x >= 4) where kappa > W/sqrt (must be none): %d ;  where kappa_upper < W/sqrt: %d ;  points compared: %d"
          % (len(viol), len(viol2), len([x for x in rows if x >= 4])))
    viol3 = [x for x, r in rows.items() if x >= 2000 and NC_bound(x, kfun=kstar) > -r[0] * math.log(x) / math.sqrt(x)]
    print("   grid points (x >= 2000) where the N_C bound exceeds the scan value (must be none): %d" % len(viol3))
    # 3. grid argument for the uniform constant
    for target in (mp.mpf('0.88'), mp.mpf('0.885'), mp.mpf('0.875'), mp.mpf('0.5')):
        rho = mp.mpf(10) ** (mp.mpf(1) / 100)
        j = 100 * 14
        last_fail = None
        vals = []
        while j >= 100 * 3:
            a = mp.mpf(10) ** (mp.mpf(j) / 100)
            beta = NC_bound(a, kfun=kstar, x_hi=a * rho)
            vals.append((a, beta))
            if beta < target and last_fail is None:
                last_fail = a
                break
            j -= 1
        amin = vals[-2][0] if last_fail is not None else vals[-1][0]
        mn = min(v for a, v in vals if a >= amin)
        XA = mp.mpf(10) ** mp.mpf('6.6')
        if amin <= XA:
            print("\n   [target %s: min of beta(a) over the cells with 10^6.6 = %s <= a <= 10^14: %s]"
                  % (mp.nstr(target, 4), mp.nstr(XA, 9), mp.nstr(min(v for a, v in vals if a >= XA * (1 - mp.mpf(10) ** -12)), 8)))
        print("\n3. target %s: beta(a) >= target for every cell [a, 10^(1/100) a] with %s <= a <= 10^14 (min beta there = %s); fails at a = %s"
              % (mp.nstr(target, 4), mp.nstr(amin, 8), mp.nstr(mn, 8), mp.nstr(last_fail, 8) if last_fail else None))
        # x >= 10^14:
        xb = mp.mpf(10) ** 14
        kinf = min(KFE, KFI)
        # pieces V = 1..64: all y >= x/65 > 10^12, bracket = 1 + (A(64) - A(1) - c_a)/log x - corr > 1
        gain = A(64) - A(1) - CA
        corr14 = HALFSQRTPI * 2 * mp.sqrt(65) * (1 + mp.log(xb / 64)) / mp.log(xb / 64) ** 2 * mp.log(xb) / xb
        print("   for x >= 10^14 (pieces up to V = 64, all y >= x/65 > 10^12): kstar >= %s;  A(64) - A(1) - c_a = %s > 0, correction <= %s ;  => N_C >= %s"
              % (mp.nstr(kinf, 8), mp.nstr(gain, 6), mp.nstr(corr14, 3), mp.nstr(kinf - corr14, 8)))
    # 4. asymptotics
    print("\n4. kappa at large y:")
    for e in (12, 16, 17, 20, 30, 50, 100, 300):
        y = mp.mpf(10) ** e
        print("   y = 1e%d: kappa = %s, kflat = %s, E_tail/sqrt = %s" % (e, mp.nstr(kappa(y), 10), mp.nstr(kflat(y), 10), mp.nstr(Etail_over_sqrt(y), 4)))
