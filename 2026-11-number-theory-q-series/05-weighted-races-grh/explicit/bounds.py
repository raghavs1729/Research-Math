#!/usr/bin/env python3
"""
bounds.py -- evaluates the explicit lower bounds of EXPLICIT.md and compares them with the sieve.

For chi = chi_4 (q=4) or chi_3 (q=3), under GRH(chi), Proposition 4 of EXPLICIT.md gives, for t >= p0:

   Theta1(t) >= t^{3/2} L(t),
   L(t) = 2/3 - sig0 - (4/5) c1 t^{-1/4} - c_theta * Eodd(t) - 1_{t > Y1^2} [ (2/3) eps_low + log(l) t^{-1/2} ],
   Eodd(t) = (3/4) t^{-1/6} + (5/6) t^{-3/10} + t^{-5/14} max(0, (log t/log p0 - 5)/2),

 with  sig0 = log(q/pi) - gamma + 2 L'/L(1,chi)           (closed form, = sum_rho |rho|^{-2} under GRH)
       c1   = sup_{0<y<=Y1} (y - vo(y))/sqrt(y)           (sieve, Y1 = 10^10),  vo(y) < y for y <= Y1 (sieve)
       eps_low, eps_up: |theta(y) - y| <= eps y for y >= Y1  (cited),  c_theta = 1 + eps_up
       l = the prime dividing q, p0 = least prime not dividing q.
 Variants:  'B'  : Broadbent-Kadiri-Lumley-Ng-Wilk, Math. Comp. 90 (2021), arXiv:2002.11068v2, Corollary 11.2
                   (k=3, M_3=1, X_0=89 967 803):  |theta(x)-x| < x/log^3 x,  so eps_low = eps_up = 1/log^3(Y1)
            'C'  : certified inputs only.  Y1 = 2^22 (Lemma 4(a) certified in ball arithmetic by cert_small.py, part C), and for
                   y >= 2^22 the row k=2, M_2=0.2, X_0=3 594 641 of the same Corollary 11.2:  eps_low = eps_up = 0.2/log^2(2^22)
            'RS' : Rosser-Schoenfeld 1962 (3.15), (3.16), in the form quoted by Johnston (Lemma 2.2):
                   x - x/log x < theta(x) < x + x/(2 log x) for x >= 41:  eps_up = 1/(2 log Y1), eps_low = 1/log Y1

 Corollary 3 (N >= 1/2): certified computation up to 10^10 (cert_gaps.py), and for x >= 10^10 the grid
 argument printed at the end of each report (margins of the double-precision evaluation are printed).

 Then (Theorem 2):  N(x) = R(x) log x / sqrt x >= L(x) + (Lstar(t1,x)/log x) (8/3 - 2 sqrt(t1/x) - (2/3)(t1/x)^{3/2})
 for every 10^4 <= t1 <= x, where Lstar(t1,x) = inf_{t1<=t<=x} L(t)  (= L(t1) if x <= Y1^2, else min(L(t1), L(Y1^2 +))).

 Usage: python3 bounds.py
"""
import os, sys, math
import mpmath as mp
import constants as K

HERE = os.path.dirname(os.path.abspath(__file__))

def parse_scan(q):
    d = dict(chk={}, bins={})
    for line in open(os.path.join(HERE, 'out_q%d_1e10.txt' % q)):
        f = line.split()
        if line.startswith('# cesaro_scan'):
            d['Y1'] = int(line.split('X=')[1].split()[0])
        elif line.startswith('CHK'):
            kv = dict(t.split('=') for t in f[1:]); d['chk'][int(kv['x'])] = kv
        elif line.startswith('C1 '):
            d['c1'] = float(line.split(':')[1].split()[0])
        elif line.startswith('CUP'):
            d['cup'] = float(line.split(':')[1].split()[0])
        elif line.startswith('EXACT'):
            d['badR'] = int(line.split(':')[1].split()[0])
        elif line.startswith('FLOAT'):
            d['badT'] = int(line.split(':')[1].split()[0])
        elif len(f) > 10 and f[0].isdigit():
            d['bins'][int(f[0])] = dict(Rmin=float(f[3]), Rat=int(f[4]), Rlow=float(f[5]), Nmin=float(f[7]), Nat=int(f[8]),
                                         Nlow=float(f[9]), Nmax=float(f[10]), Gmin=float(f[13]), Gat=int(f[14]), Glow=float(f[15]), Gmax=float(f[16]))
    return d

class Bound:
    def __init__(self, q, variant='B'):
        self.q = q; self.variant = variant
        self.sc = parse_scan(q)
        c = K.consts(q)
        self.sig0 = float(c['sig0']); self.b = float(c['b']); self.C = float(c['C'])
        self.Y1 = float(self.sc['Y1'])
        # round the sieve constant UP in the 6th decimal
        self.c1 = math.ceil(self.sc['c1'] * 1e6) / 1e6
        assert self.sc['cup'] < 0          # vo(y) < y on (0, Y1]
        if variant == 'C':
            self.Y1 = 2.0 ** 22
        LY = math.log(self.Y1)
        if variant == 'B':
            self.eps_low = self.eps_up = 1 / LY ** 3
        elif variant == 'C':
            self.eps_low = self.eps_up = 0.2 / LY ** 2
        else:
            self.eps_up = 1 / (2 * LY); self.eps_low = 1 / LY
        self.cth = 1 + self.eps_up
        self.l = 2 if q == 4 else 3
        self.p0 = 3 if q == 4 else 2
        self.c0 = 2.0 / 3 - self.sig0
        # threshold beyond which t^{-5/14}(log t/log p0 - 5) decreases
        self.tmon = math.exp(5 * math.log(self.p0) + 14.0 / 5)

    def Eodd(self, t):
        return 0.75 * t ** (-1 / 6) + (5 / 6) * t ** (-0.3) + t ** (-5 / 14) * max(0.0, (math.log(t) / math.log(self.p0) - 5) / 2)

    def terms(self, t, right=False):
        big = (t > self.Y1 ** 2) or right
        return dict(sq=0.8 * self.c1 * t ** (-0.25), odd=self.cth * self.Eodd(t),
                    far=((2 / 3) * self.eps_low + math.log(self.l) * t ** (-0.5)) if big else 0.0)

    def L(self, t, right=False):
        T = self.terms(t, right)
        return self.c0 - T['sq'] - T['odd'] - T['far']

    def Lstar(self, t1, x=float('inf')):
        """inf_{t1 <= t <= x} L(t), valid for t1 >= tmon: L is nondecreasing on [tmon, Y1^2] and on (Y1^2, oo)"""
        assert t1 >= self.tmon
        if t1 > self.Y1 ** 2 or x <= self.Y1 ** 2:
            return self.L(t1)
        return min(self.L(t1), self.L(self.Y1 ** 2, right=True))

    def t0(self):
        lo, hi = self.tmon, 1e12
        assert self.L(hi) > 0
        if self.L(lo) > 0: return lo
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if self.L(mid) > 0: hi = mid
            else: lo = mid
        return hi

    def phi(self, u):
        return 8 / 3 - 2 * math.sqrt(u) - (2 / 3) * u ** 1.5

    def B(self, x, t1=None, t1min=1e4):
        """refined lower bound for N(x); t1 optimised over a grid if not given"""
        if t1 is None:
            best = self.L(x)
            for j in range(0, 400):
                t = t1min * (x / t1min) ** (j / 400)
                best = max(best, self.L(x) + self.Lstar(t, x) / math.log(x) * self.phi(t / x))
            return best
        return self.L(x) + self.Lstar(t1, x) / math.log(x) * self.phi(t1 / x)

def fmt(x):
    return '%.0e' % x if x >= 1e7 else '%d' % x

def report(q, variant):
    b = Bound(q, variant)
    sc = b.sc
    print('=' * 110)
    print('q = %d   variant %s   (Lemma 4(a) used up to Y1 = %.6g,  Y1^2 = %.4g)' % (q, variant, b.Y1, b.Y1 ** 2))
    print('  sig0 = %.15f   c0 = 2/3 - sig0 = %.15f' % (b.sig0, b.c0))
    print('  c1 (rounded up) = %.6f   eps_low = %.3e  eps_up = %.3e  c_theta = %.6f   b = %.6f  C = 1+b0 = %.6f' % (b.c1, b.eps_low, b.eps_up, b.cth, b.b, b.C))
    print('  b t - log t - C at t = p0 = %d : %.4f  (>0 and increasing for t > 1/b)' % (b.p0, b.b * b.p0 - math.log(b.p0) - b.C))
    print('  L(t) is nondecreasing on [tmon, Y1^2] with tmon = %.1f' % b.tmon)
    t0 = b.t0()
    print('  analytic positivity threshold: L(t) > 0 for all t >= %.1f ;  L(10^4) = %.5f ;  limit c0 - (2/3)eps_low = %.6f' % (t0, b.L(1e4), b.c0 - (2 / 3) * b.eps_low))
    print('  sieve: R(p) <= 0 at %d primes, Theta1(p) <= 0 at %d primes (after the first race prime, p <= %.0e)' % (sc['badR'], sc['badT'], b.Y1))
    print('\n  table of L(t) and its error terms')
    print('  %-10s %-10s %-10s %-10s %-10s %-10s' % ('t', 'squares', 'odd k>=3', 't>Y1^2', 'L(t)', 'Lstar(t)'))
    Y2 = b.Y1 ** 2
    for t in [t for t in [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e12, 1e15, 1e20] if t < Y2] + [Y2]:
        T = b.terms(t)
        print('  %-10.3e %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f%s' % (t, T['sq'], T['odd'], T['far'], b.L(t), b.Lstar(t), '   (= Y1^2)' if t == Y2 else ''))
    T = b.terms(Y2, right=True)
    print('  %-10s %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f' % ('Y1^2 +', T['sq'], T['odd'], T['far'], b.L(Y2, True), b.L(Y2, True)))
    for t in [t for t in [1e15, 1e20, 1e30, 1e50, 1e100] if t > Y2]:
        T = b.terms(t)
        print('  %-10.3e %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f' % (t, T['sq'], T['odd'], T['far'], b.L(t), b.Lstar(t)))

    print('\n  consistency with the sieve, G(t) = Theta1(t)/t^{3/2}:  need L(t) <= G(t) pointwise;  rigorous test L(2^k) <= G_low(bin k);')
    print('  the last column is the heuristic comparison of the largest value of L on the bin with the smallest value of G on it')
    print('  %-4s %-12s %-10s %-10s %-10s %-8s %-12s' % ('k', '2^k', 'L(2^k)', 'G_low', 'G_max', 'ok', 'L(2^(k+1))<=G_low'))
    allok = True
    for k in sorted(sc['bins']):
        if (1 << k) < 1e4: continue
        B_ = sc['bins'][k]
        ok = b.L(float(1 << k)) <= B_['Glow']; allok &= ok
        print('  %-4d %-12d %-10.5f %-10.5f %-10.5f %-8s %-12s' % (k, 1 << k, b.L(float(1 << k)), B_['Glow'], B_['Gmax'], ok, b.L(float(1 << (k + 1))) <= B_['Glow']))
    print('  all rigorous consistency tests passed:', allok)

    print('\n  lower bounds for N(x) = R(x) log x/sqrt x  (simple: L(x);  refined: B(x) with t1 optimised) against exact sieve values')
    print('  %-12s %-12s %-12s %-12s %-12s %-10s' % ('x', 'L(x)', 'B(x)', 'N(x) exact', 'R(x) exact', 'bound*sqrt x/log x'))
    for X in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9, 10 ** 10):
        kv = sc['chk'][X]
        Bx = b.B(float(X))
        assert Bx <= float(kv['N'])
        print('  %-12.0e %-12.5f %-12.5f %-12.6f %-12.6f %-10.4f' % (X, b.L(float(X)), Bx, float(kv['N']), float(kv['R']), Bx * math.sqrt(X) / math.log(X)))
    print('  at the minimisers of N in each dyadic bin (x >= 2^14):')
    print('  %-4s %-12s %-10s %-10s %-10s' % ('k', 'x_min', 'B(x_min)', 'N_min', 'ratio'))
    worst = 0
    for k in sorted(sc['bins']):
        if (1 << k) < 1e4: continue
        B_ = sc['bins'][k]
        Bx = b.B(float(B_['Nat']))
        assert Bx <= B_['Nmin']
        worst = max(worst, Bx / B_['Nmin'])
        print('  %-4d %-12d %-10.5f %-10.5f %-10.3f' % (k, B_['Nat'], Bx, B_['Nmin'], Bx / B_['Nmin']))
    print('  largest ratio bound/true at bin minimisers: %.3f' % worst)

    # ---- Corollary 3:  N(x) >= 1/2.
    # For x in [a, rho a] (a >= 10^6), Theorem 2 with t1 = a/100 gives
    #     N(x) >= beta(a) := inf_{[a, rho a]} L  +  (inf_{t >= a/100} L) * phi(1/100) / log(rho a)
    # (t1/x <= 1/100, log x <= log(rho a); the infima are Lstar, which accounts for the jump of L at Y1^2).
    target = 0.5
    rho = 10 ** 0.02
    def beta(a):
        return b.Lstar(a, rho * a) + b.Lstar(a / 100) * b.phi(0.01) / math.log(rho * a)
    XE = 1e10          # limit of the certified computation (cert_gaps.py)
    AT = 1e30          # beyond AT the bound N(x) >= inf_{t >= AT} L(t) is used
    JT = int(round(math.log(AT / XE) / math.log(rho)))
    grid = [XE * rho ** j for j in range(0, JT)]
    margins = [beta(a) - target for a in grid]
    mg = min(margins); ag = grid[margins.index(mg)]
    mt = b.Lstar(AT) - target
    print('\n  Corollary 3 (N(x) >= %.2f), analytic part for x >= %.0e:' % (target, XE))
    if mg > 0 and mt > 0:
        print('    grid a = 10^10 * 10^(j/50), a < 10^30:  min (beta(a) - 1/2) = %.5f  (at a = %.3e);  beta(10^10) - 1/2 = %.5f' % (mg, ag, margins[0]))
        print('    tail x >= 10^30:  inf_{t >= 10^30} L(t) - 1/2 = %.5f   (limit of L - 1/2 = %.5f)' % (mt, b.c0 - (2 / 3) * b.eps_low - target))
        print('    => N(x) >= 1/2 for all x >= 10^10 with margin >= %.5f in the double-precision evaluation' % min(mg, mt))
    else:
        print('    NOT available in variant %s: grid margin %.5f, tail margin %.5f (limit of L = %.5f)' % (variant, mg, mt, b.c0 - (2 / 3) * b.eps_low))
    # for the record: how far down the same argument reaches, and how thin its margin is there
    j = 0; amin = XE
    while True:
        a = XE * rho ** (j - 1)
        if a < 1e6 or beta(a) < target: break
        amin = a; j -= 1
    firsts = [(XE * rho ** (j + i), beta(XE * rho ** (j + i)) - target) for i in range(0, 4)]
    print('    (record) the grid argument alone works down to a = %.4e; margins at the first grid points: %s'
          % (amin, ', '.join('%.3e: %.1e' % f for f in firsts)))
    return b

if __name__ == '__main__':
    for q in (4, 3):
        for v in ('B', 'C', 'RS'):
            report(q, v)
