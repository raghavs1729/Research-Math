#!/usr/bin/env python3
"""
crosscheck.py -- component-wise numerical tests (double precision / mpmath; NOT certified) of the inequalities
behind Theorems 3 and 4 of CHEBYSHEV1853.md, so that a wrong constant in one component cannot hide behind slack
in another.

 A. Theorem 3, component by component, at y = 10^1 .. 10^6 (direct sums over primes / prime powers):
      squares        (sqrt(pi)/2) sqrt y - c1 Gamma(5/4) y^(1/4) - E(y)  <=  Q_2(y)  <=  (sqrt(pi)/2) sqrt y + E(y)
      fourth powers  Q_4(y) >= (1-eps_a) Gamma(5/4) y^(1/4) - c1 Gamma(9/8) y^(1/8)
      odd powers     |W_k(y)| <= B_k(y),  k = 3, 5, 7  (p0^k < y)
      high powers    sum over k with p0^k >= y (each parity)  <= 1
      zeros          -S_Gamma sqrt y <= -psi_e(y,chi) <= S_Gamma sqrt y + b + 0.2555/y^2     (y >= 4; uses GRH: this is a
                     test of Lemma 3 + Lemma 4 against the data, the zeros being on the line as far as computed)
    Also the squares inequality at every grid point of the scan (column Wsq).
 B. Theorem 4 with the TRUE function W in place of kappa (scan values): the right side, with kappa_* replaced by
    the minimum of W(y)/sqrt(y) over the grid points of each piece, must be <= N_C(x) (scan).  This tests
    Lemma 5 and the constants c_a, A(V) and the correction term independently of Theorem 3.
Usage: python3 crosscheck.py [q]
"""
import sys, os, math
import numpy as np
import mpmath as mp
sys.argv = sys.argv[:1] + (sys.argv[1:] or ['4'])
import bounds_abel as B          # reads q from sys.argv[1]

q = B.q
HERE = os.path.dirname(os.path.abspath(__file__))
ok_all = True


def flag(cond):
    global ok_all
    ok_all = ok_all and bool(cond)
    return "ok" if cond else "FAIL"


N = 60 * 10 ** 6
sv = np.ones(N + 1, dtype=np.bool_); sv[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = False
P = np.nonzero(sv)[0].astype(np.float64)
del sv
if q == 4:
    P = P[P != 2]; chi = np.where(P % 4 == 1, 1.0, -1.0)
else:
    P = P[P != 3]; chi = np.where(P % 3 == 1, 1.0, -1.0)
LP = np.log(P)
p0 = float(B.p0)
bq = 0.7831887854 if q == 4 else 0.9481988267

print("A. components of Theorem 3 (q = %d)" % q)
for e in range(1, 7):
    y = 10.0 ** e
    lim = 50 * y
    sel = P <= lim
    psi = 0.0; Q = {}; Wk = {}
    k = 1
    while p0 ** k <= lim:
        m = P ** k <= lim
        w = np.exp(-P[m] ** k / y)
        if k % 2 == 1:
            t = float(np.sum(chi[m] * LP[m] * w)); psi += t; Wk[k] = -t          # W_k = sum s log p e^{-p^k/y}
        else:
            t = float(np.sum(LP[m] * w)); psi += t; Q[k] = t
        k += 1
    W = Wk[1]
    r = math.sqrt(y)
    E = float(B.Etail_over_sqrt(y)) * r
    c1 = float(B.c1); epsa = float(B.eps_a); G54 = float(B.G54); G98 = float(B.G98)
    lo2 = float(B.HALFSQRTPI) * r - c1 * G54 * y ** 0.25 - E
    up2 = float(B.HALFSQRTPI) * r + E
    print("  y = 1e%d:  W = %.6f" % (e, W))
    print("     squares:  %.5f <= Q_2 = %.5f <= %.5f   %s" % (lo2, Q.get(2, 0.0), up2, flag(lo2 <= Q.get(2, 0.0) <= up2)))
    lo4 = (1 - epsa) * G54 * y ** 0.25 - c1 * G98 * y ** 0.125
    print("     fourth:   Q_4 = %.5f >= %.5f   %s" % (Q.get(4, 0.0), lo4, flag(Q.get(4, 0.0) >= lo4)))
    for kk in (3, 5, 7):
        if p0 ** kk < y:
            bk = float(B.Bk(kk, mp.mpf(y)))
            print("     k = %d:    |W_k| = %.5f <= B_k = %.5f   %s" % (kk, abs(Wk.get(kk, 0.0)), bk, flag(abs(Wk.get(kk, 0.0)) <= bk)))
    hi_odd = sum(abs(v) for kk, v in Wk.items() if kk >= 3 and p0 ** kk >= y)
    hi_even = sum(v for kk, v in Q.items() if kk >= 4 and p0 ** kk >= y)
    # the (few) powers with p0^k > 50 y contribute < e^{-50} each and are not summed
    print("     high powers (p0^k >= y): odd %.5f, even %.5f  (<= 1 each)   %s" % (hi_odd, hi_even, flag(hi_odd <= 1 and hi_even <= 1)))
    if y >= 4:
        SG = float(B.SG)
        print("     zeros:    %.6f <= -psi_e = %.6f <= %.6f   %s   [b = %.6f; -psi_e - b = %+.6f, S_Gamma sqrt y = %.6f]"
              % (-SG * r, -psi, SG * r + bq + 0.2555 / y ** 2, flag(-SG * r <= -psi <= SG * r + bq + 0.2555 / y ** 2), bq, -psi - bq, SG * r))
    # reassemble
    tot = -psi + Q.get(2, 0.0) - sum(v for kk, v in Wk.items() if kk >= 3) + sum(v for kk, v in Q.items() if kk >= 4)
    print("     decomposition: -psi_e + Q_2 - sum_odd W_k + sum_even Q_k = %.9f  vs W = %.9f   %s" % (tot, W, flag(abs(tot - W) < 1e-7 * max(1, abs(W)))))

rows = []
for fn in ('out_q%d_1e8.txt' % q, 'out_q%d_1e9.txt' % q):
    try:
        for line in open(os.path.join(HERE, fn)):
            if line[0] != '#':
                f = line.split(); rows.append((float(f[1]), float(f[2]), float(f[3]), float(f[4])))
    except FileNotFoundError:
        pass
rows.sort()
bad = 0; n = 0; worst = 1e9
for x, C, W, Wsq in rows:
    if x < 4:
        continue
    r = math.sqrt(x)
    E = float(B.Etail_over_sqrt(x)) * r
    lo2 = float(B.HALFSQRTPI) * r - float(B.c1) * float(B.G54) * x ** 0.25 - E
    up2 = float(B.HALFSQRTPI) * r + E
    n += 1
    if not (lo2 <= Wsq <= up2):
        bad += 1
    worst = min(worst, (up2 - Wsq) / r)
print("  squares inequality at the %d scan grid points y >= 4: violations %d ; min of (upper - Q_2)/sqrt y = %.3e   %s"
      % (n, bad, worst, flag(bad == 0)))

print("\nB. Theorem 4 with the true W (scan) in place of kappa")
xs = np.array([r[0] for r in rows]); Wn = np.array([r[2] / math.sqrt(r[0]) for r in rows])
CA = float(B.CA)
viol = 0; cnt = 0
print("   %10s %10s %10s %8s" % ("x", "bound", "N_C scan", "ratio"))
for x, C, W, Wsq in rows:
    if x < 2000 or x > 1.0000001e8:
        continue
    lx = math.log(x)

    def piece_min(ylo, yhi):
        m = (xs >= ylo * (1 - 1e-12)) & (xs <= yhi * (1 + 1e-12))
        # include the grid point just below ylo so that the minimum over the grid is a lower bound for an increasing-ish W/sqrt
        idx = np.nonzero(m)[0]
        i0 = max(idx[0] - 1, 0) if len(idx) else np.searchsorted(xs, ylo) - 1
        i1 = idx[-1] if len(idx) else i0
        return float(np.min(Wn[i0:i1 + 1]))
    k0 = piece_min(x / 2, x)
    val = k0 * (1 - CA / lx)
    Vn = 1
    vs = B.VS
    for i in range(1, len(vs)):
        Vp, Vc = vs[i - 1], vs[i]
        if x / (1 + Vc) < 500:
            break
        ki = piece_min(x / (1 + Vc), x / (1 + Vp))
        val += ki * float(B.A(Vc) - B.A(Vp)) * lx / math.log(x / Vp) ** 2
        Vn = Vc
    An = math.log(x / Vn)
    val -= math.sqrt(math.pi) * math.sqrt(1 + Vn) * (1 + An) / An ** 2 * lx / x
    NC = -C * lx / math.sqrt(x)
    cnt += 1
    if val > NC:
        viol += 1
    if abs(math.log10(x) - round(math.log10(x) * 2) / 2) < 1e-9:
        print("   %10.4g %10.6f %10.6f %8.4f" % (x, val, NC, val / NC))
print("   grid points tested: %d ; violations (bound > N_C): %d   %s" % (cnt, viol, flag(viol == 0)))
print("\nALL OK:", ok_all)
