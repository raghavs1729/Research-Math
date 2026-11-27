#!/usr/bin/env python3
"""
verify_small.py -- independent (Python) checks of everything the C scan and the proof rely on,
on the range x <= NMAX = 2^21.  Usage: python3 verify_small.py [q]   (q = 4 or 3)

 1. recomputes D0, D1, Theta, Theta1, R, N, G at primes and compares the dyadic-bin minima with
    the output of cesaro_scan (out_q{q}_1e10.txt);
 2. exact (rational) positivity of R at primes, and the tiny range below 11 treated by a dense grid;
 3. numerical check of the summation identities
       x R(x) = Theta1(x)/log x + int Theta1(t) h''(t) dt                       (Lemma 2)
       I_{m,m'}(x) = Phi(x)/log x + int Phi(t) (-g'(t)) dt                      (Lemma 6)
 4. symbolic check of h'' and g';
 5. the cone m <= 0 < m' <= 1, m < m': min of I_{m,m'} and of I - kappa x^{m-1} Theta1/log x at primes.
"""
import sys, os, math, re
from fractions import Fraction
import numpy as np
import mpmath as mp

q = int(sys.argv[1]) if len(sys.argv) > 1 else 4
NMAX = 1 << 21
HERE = os.path.dirname(os.path.abspath(__file__))
p0 = 3 if q == 4 else 2

sv = np.ones(NMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, int(NMAX ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = False
allp = np.nonzero(sv)[0]
if q == 4:
    P = allp[allp > 2]; s = np.where(P % 4 == 3, 1, -1)
else:
    P = allp[allp != 3]; s = np.where(P % 3 == 2, 1, -1)
P = P.astype(np.int64); s = s.astype(np.int64)
LD = np.longdouble
lp = np.log(P.astype(LD))
D0 = np.cumsum(s); D1 = np.cumsum(s * P)
Th = np.cumsum(s * lp); W = np.cumsum(s * P * lp)
x = P.astype(LD)
R = D0 - D1 / x
T1 = x * Th - W
N = R * np.log(x) / np.sqrt(x)
G = T1 / x ** 1.5
print('q=%d  primes in race up to 2^21: %d' % (q, len(P)))

# 1. compare with the C scan
fn = os.path.join(HERE, 'out_q%d_1e10.txt' % q)
ok = True
if os.path.exists(fn):
    for line in open(fn):
        f = line.split()
        if len(f) > 10 and f[0].isdigit():
            k = int(f[0])
            if k < 3 or k > 20: continue
            m = (P >= max(11, 1 << k)) & (P < (1 << (k + 1)))
            vals = dict(Rmin=float(f[3]), Rat=int(f[4]), Nmin=float(f[7]), Nat=int(f[8]), Nmax=float(f[10]), Gmin=float(f[13]), Gat=int(f[14]), Gmax=float(f[16]))
            mine = dict(Rmin=float(R[m].min()), Rat=int(P[m][R[m].argmin()]), Nmin=float(N[m].min()), Nat=int(P[m][N[m].argmin()]),
                        Nmax=float(N[m].max()), Gmin=float(G[m].min()), Gat=int(P[m][G[m].argmin()]), Gmax=float(G[m].max()))
            for key in vals:
                if abs(vals[key] - mine[key]) > 2e-6 * max(1, abs(vals[key])):
                    ok = False; print('  MISMATCH bin', k, key, vals[key], mine[key])
    print('1. dyadic bins k=3..20 of the C scan reproduced by numpy:', 'OK' if ok else 'FAILED')
    for line in open(fn):
        if line.startswith('CHK'):
            d = dict(t.split('=') for t in line.split()[1:])
            X = int(d['x'])
            if X <= NMAX:
                i = np.searchsorted(P, X, side='right') - 1
                Rx = Fraction(int(D0[i])) - Fraction(int(D1[i]), X)
                t1 = X * Th[i] - W[i]
                assert int(d['D0']) == D0[i] and int(d['D1']) == D1[i]
                assert abs(float(d['R']) - float(Rx)) < 1e-9 and abs(float(d['Theta1']) - float(t1)) < 1e-6 * abs(float(t1))
    print('   checkpoints x=10^2..10^6 (D0, D1 exact; R, Theta1) agree')

# 2. exact positivity of R at primes; small range
bad = [int(P[i]) for i in range(1, len(P)) if int(P[i]) * int(D0[i]) - int(D1[i]) <= 0]
print('2. exact: primes p (after the first) <= 2^21 with p*D0-D1 <= 0:', bad, ' ; primes with Theta1<=0:', int((T1[1:] <= 0).sum()))
def state(xx):
    i = np.searchsorted(P, xx, side='right') - 1
    if i < 0: return 0, 0, 0.0, 0.0
    return int(D0[i]), int(D1[i]), float(Th[i]), float(W[i])
grid = np.linspace(p0 + 1e-9, 64, 600001)
Rg = []; Ng = []; Gg = []
for xx in grid:
    a, b, th, w = state(xx)
    r = a - b / xx; Rg.append(r); Ng.append(r * math.log(xx) / math.sqrt(xx)); Gg.append((xx * th - w) / xx ** 1.5)
Rg = np.array(Rg); Ng = np.array(Ng); Gg = np.array(Gg)
print('   dense grid on (%d, 64]: min R = %.6g (at x=%.6f), all R>0: %s; all Theta1>0: %s' % (p0, Rg.min(), grid[Rg.argmin()], bool((Rg > 0).all()), bool((Gg > 0).all())))
for lo in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
    m = grid >= lo
    print('   inf over grid x in [%2d, 64] of N(x) = %.6f at x = %.5f ;  of G(x) = Theta1/x^1.5: %.6f at x = %.5f'
          % (lo, Ng[m].min(), grid[m][Ng[m].argmin()], Gg[m].min(), grid[m][Gg[m].argmin()]))
# last x with N(x) < 1/2 below 2^21 (N is evaluated on gaps through the lower bound min(R(p),R(p')) * w(p'))
w = np.log(x) / np.sqrt(x)
gap_low = np.minimum(R[:-1], R[1:]) * w[1:]
idx = np.nonzero((gap_low < 0.5) & (P[:-1] >= 11))[0]
print('   gaps [p,p\'] with p>=11 whose rigorous lower bound for N is < 1/2:', [(int(P[i]), int(P[i + 1])) for i in idx])

# 3. identities
mp.mp.dps = 25
def check_identity_R(X):
    X = mp.mpf(X)
    pts = [int(pp) for pp in P if pp <= X]
    lhs = sum(int(si) * (X - pp) for si, pp in zip(s, pts))          # x R(x)
    def h2(t):
        L = mp.log(t)
        return 2 / (t * L ** 2) + (X - t) * (1 / (t ** 2 * L ** 2) + 2 / (t ** 2 * L ** 3))
    tot = mp.mpf(0); th = mp.mpf(0); wv = mp.mpf(0)
    edges = pts + [X]
    for j, pp in enumerate(pts):
        th += int(s[j]) * mp.log(pp); wv += int(s[j]) * pp * mp.log(pp)
        a, b = edges[j], edges[j + 1]
        if b > a:
            tot += mp.quad(lambda t: (t * th - wv) * h2(t), [a, b])
    rhs = (X * th - wv) / mp.log(X) + tot
    return lhs, rhs
for X in ('200.5', '1000.25'):
    l, r = check_identity_R(X)
    print('3. Lemma 2 at x=%s:  x R(x) = %s,  Theta1/log x + int Theta1 h\'\' = %s,  diff = %s' % (X, mp.nstr(l, 15), mp.nstr(r, 15), mp.nstr(l - r, 3)))

def check_identity_pair(X, m, m2):
    X = mp.mpf(X); k = mp.mpf(m2) - m
    pts = [int(pp) for pp in P if pp <= X]
    lhs = sum(int(si) * mp.mpf(pp) ** m * (1 - (pp / X) ** k) for si, pp in zip(s, pts))
    def mg(t):
        L = mp.log(t)
        return (1 + (t / X) ** k) / (t * L ** 2) + 2 * (1 - (t / X) ** k) / (k * t * L ** 3)
    tot = mp.mpf(0); A = mp.mpf(0); Bv = mp.mpf(0)
    edges = pts + [X]
    Phi = lambda t: A - t ** (-k) * Bv
    for j, pp in enumerate(pts):
        A += int(s[j]) * mp.log(pp) * mp.mpf(pp) ** m; Bv += int(s[j]) * mp.log(pp) * mp.mpf(pp) ** m2
        a, b = edges[j], edges[j + 1]
        if b > a:
            tot += mp.quad(lambda t: Phi(t) * mg(t), [a, b])
    rhs = Phi(X) / mp.log(X) + tot
    return lhs, rhs
for (m, m2) in ((0, 1), (1, 2), (0.5, 3), (-1, 0.5)):
    l, r = check_identity_pair('700.3', m, m2)
    print('   Lemma 6 at x=700.3, (m,m\')=(%s,%s):  I = %s,  Phi/log x + int Phi(-g\') = %s,  diff = %s' % (m, m2, mp.nstr(l, 15), mp.nstr(r, 15), mp.nstr(l - r, 3)))

# 4. symbolic
try:
    import sympy as sp
    t, X, k = sp.symbols('t x kappa', positive=True)
    h = (X - t) / sp.log(t)
    h2 = 2 / (t * sp.log(t) ** 2) + (X - t) * (1 / (t ** 2 * sp.log(t) ** 2) + 2 / (t ** 2 * sp.log(t) ** 3))
    g = 1 / sp.log(t) + (1 - (t / X) ** k) / (k * sp.log(t) ** 2)
    mgp = (1 + (t / X) ** k) / (t * sp.log(t) ** 2) + 2 * (1 - (t / X) ** k) / (k * t * sp.log(t) ** 3)
    print('4. sympy: h\'\' formula residual =', sp.simplify(sp.diff(h, t, 2) - h2), ';  -g\' formula residual =', sp.simplify(-sp.diff(g, t) - mgp),
          ';  h\'(x) + 1/log x =', sp.simplify(sp.diff(h, t).subs(t, X) + 1 / sp.log(X)))
except ImportError:
    print('4. sympy not available')

# 5. the cone m <= 0, m < m' <= 1
xf = P.astype(np.float64); lf = np.log(xf); sf = s.astype(np.float64)
T1f = T1.astype(np.float64)
print('5. cone check at all race primes p0 < p <= 2^21  (I is continuous and monotone between primes, so primes suffice):')
for (m, m2) in ((0, 1), (0, 0.5), (0, 0.01), (-0.5, 0.5), (-0.25, 0.25), (-1, 0), (-1, 1), (-2, 1), (-3, -1), (-0.01, 0.0)):
    k = m2 - m
    Dm = np.cumsum(sf * xf ** m); Dm2 = np.cumsum(sf * xf ** m2)
    I = Dm - xf ** (-k) * Dm2
    lowb = k * xf ** (m - 1) * T1f / lf
    Mx = m + 0.5
    nrm = I * lf / (k * xf ** Mx)
    print('   (m,m\')=(%5s,%5s): min I = %.3e at p=%d ; min (I - kappa p^(m-1) Theta1/log p)/I = %.3e ; min of I log p/(kappa p^(m+1/2)) over p>=10^4: %.4f'
          % (m, m2, I[1:].min(), P[1:][I[1:].argmin()], ((I[1:] - lowb[1:]) / I[1:]).min(), nrm[xf >= 1e4].min()))
