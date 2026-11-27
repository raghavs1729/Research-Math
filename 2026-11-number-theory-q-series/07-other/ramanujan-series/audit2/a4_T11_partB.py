#!/usr/bin/env python3
"""
audit2 / a4: Theorem 11, Part B  (T_N := sum_{k<N} a(N,k) -> sqrt2/(4 pi^2)), step by step.

Everything is recomputed from the definition of a(N,k) in RESULTS.md:
  a(N,k) = C_N (2k+1) Gamma(4k+2N+1/2) / ((N-k-1)!^2 (N+k)!^2 Gamma(4k-2N+9/2)),
  C_N    = 2^(6-8N) Gamma(N+3/4)^2 Gamma(N+1/4) / (sqrt2 pi Gamma(N-1/4) Gamma(2N-1)).
(B1) is tested in exact rational arithmetic for all N <= 600 (every index), then in 50-digit floating
point for a few larger N; (B2)-(B5) are tested in 50-digit arithmetic for N up to 10^6, printing the
implied O-constants so that one can see that they stay bounded.
"""
import sys
from fractions import Fraction as Fr
from mpmath import mp, mpf, loggamma, log, exp, sqrt, pi, floor, nstr, diff

mp.dps = 50
OK = True
def check(cond, msg):
    global OK
    OK &= bool(cond)
    print(('  ok   ' if cond else '  FAIL ') + msg)

XS = sqrt(mpf(3))/2
def L(x):    return log((1-x)**2*(2*x+1)**4/((1+x)**2*(2*x-1)**4))
def phi(x):  return -6*log(2) + 2*(2*x+1)*log(2*x+1) - 2*(1-x)*log(1-x) - 2*(1+x)*log(1+x) - 2*(2*x-1)*log(2*x-1)
def psi0(x): return mpf(2)**2.5*pi**(-2.5)*x*(1-x)/((1+x)*(2*x-1)**4)
LIMIT = sqrt(mpf(2))/(4*pi**2)

def r_exact(N, k):          # sigma(N,k) as printed in RESULTS.md, exact
    num = Fr(2*k+3)*(N-k-1)**2; den = Fr(2*k+1)*(N+k+1)**2
    for i in range(4):
        num *= Fr(8*k+4*N+1, 2) + i; den *= Fr(8*k-4*N+9, 2) + i
    return num/den
def r_mp(N, k):
    num = mpf(2*k+3)*(N-k-1)**2; den = mpf(2*k+1)*(N+k+1)**2
    for i in range(4):
        num *= mpf(8*k+4*N+1)/2 + i; den *= mpf(8*k-4*N+9)/2 + i
    return num/den
def log_a(N, k):            # only used for 4k-2N+9/2 > 0
    logC = (6-8*N)*log(2) + 2*loggamma(N+mpf(3)/4) + loggamma(N+mpf(1)/4) - log(sqrt(mpf(2))*pi) - loggamma(N-mpf(1)/4) - loggamma(2*N-1)
    return logC + log(2*k+1) + loggamma(4*k+2*N+mpf(1)/2) - 2*loggamma(N-k) - 2*loggamma(N+k+1) - loggamma(4*k-2*N+mpf(9)/2)

# ------------------------------------------------------------------ algebraic facts
print('algebraic ingredients')
check(abs(phi(XS)) < mpf(10)**(-45) and abs(L(XS)) < mpf(10)**(-45), 'phi(x*) = 0 and L(x*) = 0 (50 digits; exact argument: 8x^3 = 6x and the log 2 / log(2+-sqrt3) bookkeeping)')
check(all(abs(diff(phi, x) - L(x)) < mpf(10)**(-30) for x in (mpf('0.81'), mpf('0.85'), XS, mpf('0.9'), mpf('0.92'))), "phi' = L at 5 points")
check(abs(diff(L, XS) + 24) < mpf(10)**(-30), "L'(x*) = -24;  L'(x) = -4/(1-x^2) - 16/(4x^2-1)")
check(abs(psi0(XS) - sqrt(mpf(3))/(sqrt(mpf(2))*pi**2.5)) < mpf(10)**(-45) and abs(psi0(XS)*sqrt(pi/12) - LIMIT) < mpf(10)**(-45),
      'psi0(x*) = sqrt3/(sqrt2 pi^(5/2)) and psi0(x*) sqrt(pi/12) = sqrt2/(4 pi^2)')
L2max = max(abs(diff(L, mpf(8)/10 + mpf(13)*i/10000, 2)) for i in range(101))
print(f'      sup |L\'\'| on [0.8,0.93] ~ {nstr(L2max, 6)}')

# ------------------------------------------------------------------ (B1)
print('(B1) shape, exact rational arithmetic, every index, N = 3..600')
bad = 0; minmargin = None
for N in range(3, 601):
    minI = (N-1)//2 if N % 2 else N//2 - 1          # smallest integer >= N/2 - 1
    assert minI >= N/2 - 1 > minI - 1
    rs = [r_exact(N, k) for k in range(0, N)]
    if any(x < 0 for x in rs): bad += 1             # positivity of a (a(N,0) > 0 is clear from the sign pattern)
    below = rs[:minI]
    if any(x < 1 for x in below): bad += 1
    if below:
        m = min(below) - 1
        minmargin = m if minmargin is None else min(minmargin, m)
    onI = rs[minI:N]
    if any(onI[i+1] >= onI[i] for i in range(len(onI)-1)): bad += 1
    # the three displayed inequalities of the text
    for k in range(1, N//2 + 1):
        if 6*N*N*k - 8*k**3 - 2*N*N - 8*k*k < 4*N*N*(k-1): bad += 1
    for k in range(0, N):
        if 2*k <= N - 4:        # k <= N/2 - 2
            lb = (Fr((N-k-1)*(N+2*k)**2, (N+k+1)*(N-2*k)**2))**2
            if rs[k] < lb: bad += 1
    if N % 2 and N >= 3:
        k = (N-3)//2
        if rs[k] < Fr(16, 81)*Fr(8*N-11, 2)**4 or Fr(16, 81)*Fr(8*N-11, 2)**4 < 1: bad += 1
check(bad == 0, f'a > 0, r(k) >= 1 for every k < min I, r strictly decreasing on I, and the three displayed bounds: no failure for 3 <= N <= 600 '
                f'(smallest r-1 below I: {float(minmargin):.4f})')
for N in (2000, 5001, 20000):
    minI = (N-1)//2 if N % 2 else N//2 - 1
    rs = [r_mp(N, k) for k in range(N)]
    good = all(x >= 1 for x in rs[:minI]) and all(rs[i+1] < rs[i] for i in range(minI, N-1))
    check(good, f'same two properties in 50-digit arithmetic for N = {N} (min over k < min I of r(k)-1 = {nstr(min(rs[:minI])-1, 5)})')

# ------------------------------------------------------------------ (B2)-(B5)
print('(B2)-(B5): N, then the measured constants')
hdr = ('N', 'B2a', 'B2b', 'B2c', 'B3', 'B4up', 'B4lo', 'B4min', 'B5a', 'B5b', 'T_N/lim-1')
rows = []
for N in (100, 1000, 10000, 100000, 1000000):
    k0 = int(floor(XS*N)); theta = XS*N - k0; J = int(floor(mpf(N)**0.6))
    minI = (N-1)//2 if N % 2 else N//2 - 1
    # B2a: sup_{k/N in [0.8,0.93]} N |log r(k) - L(k/N)|
    step = max(1, N//20000)
    ks = range(-(-8*N//10), 93*N//100 + 1, step)
    B2a = max(N*abs(log(r_mp(N, k)) - L(mpf(k)/N)) for k in ks)
    # B2b: sup_{|j|<=J+1} |log r(k0+j) + 24 j/N| / (1/N + j^2/N^2)
    B2b = max(abs(log(r_mp(N, k0+j)) + mpf(24*j)/N)/(mpf(1)/N + mpf(j*j)/N**2) for j in range(-J-1, J+2))
    # a(k0+j) for |j| <= J by the recursion, relative to a(k0)
    rel = {0: mpf(1)}
    for j in range(0, J+1):   rel[j+1] = rel[j]*r_mp(N, k0+j)
    for j in range(0, -J-1, -1): rel[j-1] = rel[j]/r_mp(N, k0+j-1)
    # B2c: sup_{|j|<=J} |log(a(k0+j)/a(k0)) + 12 j^2/N| * N^(1/5)
    B2c = max(abs(log(rel[j]) + mpf(12*j*j)/N) for j in range(-J, J+1))*mpf(N)**0.2
    centre = sum(rel[j] for j in range(-J, J+1))
    B3 = centre/sqrt(pi*N/12) - 1
    # tails, relative to a(k0)
    rJ, rmJ = r_mp(N, k0+J), r_mp(N, k0-J-1)
    c1 = rJ <= exp(-mpf(12*J)/N); c2 = rmJ >= exp(mpf(12*J)/N)
    up = mpf(0); v = rel[J]
    for k in range(k0+J, N-1):
        v *= r_mp(N, k); up += v                       # a(k0+J+1), ..., a(N-1)
    lo = mpf(0); v = rel[-J]; amin = None
    for k in range(k0-J-1, -1, -1):
        v /= r_mp(N, k); lo += v                       # a(k0-J-1), ..., a(0)
        if k == minI: amin = v
    B4up = up/sqrt(mpf(N)); B4lo = lo/sqrt(mpf(N))
    bound = mpf(N)**0.4*exp(-11*mpf(N)**0.2)
    B4min = log(amin)/(mpf(N)**0.6)                    # claim: <= -4
    # B5
    la0 = log_a(N, k0)
    B5a = max(abs(N*(log_a(N, k) - (N*phi(mpf(k)/N) - log(mpf(N))/2 + log(psi0(mpf(k)/N))))) for k in
              (int(0.8*N)+1, int(0.83*N), k0, int(0.9*N), int(0.93*N)))
    B5b = N*(exp(la0)*sqrt(mpf(N))/psi0(XS) - 1)
    TN = exp(la0)*(centre + up + lo + 0)               # centre includes |j|<=J, up/lo the rest
    rows.append((N, B2a, B2b, B2c, B3, B4up, B4lo, B4min, B5a, B5b, TN/LIMIT - 1))
    print(f'  N={N:8d} k0={k0} J={J} theta={nstr(theta,4)}')
    print(f'     (B2) sup N|log r - L| on [0.8,0.93] = {nstr(B2a,6)};  const in -24j/N+O(1/N+j^2/N^2) = {nstr(B2b,6)};  '
          f'sup|log a(k0+j)/a(k0) + 12j^2/N| N^(1/5) = {nstr(B2c,6)}')
    print(f'     (B3) centre/(a(k0) sqrt(pi N/12)) - 1 = {nstr(B3,6)}')
    print(f'     (B4) r(k0+J) <= e^(-12J/N): {c1};  r(k0-J-1) >= e^(12J/N): {c2};  upper tail/(a(k0) sqrt N) = {nstr(B4up,4)}, '
          f'lower = {nstr(B4lo,4)}  [claimed bound N^(2/5) e^(-11 N^(1/5)) = {nstr(bound,4)}];  log(a(min I)/a(k0))/N^(3/5) = {nstr(B4min,6)} (claim <= -4)')
    print(f'     (B5) sup N|log a - (N phi - (1/2)log N + log psi0)| at 5 points of the window = {nstr(B5a,6)};  N (a(k0) sqrt N/psi0(x*) - 1) = {nstr(B5b,6)}')
    print(f'     T_N / (sqrt2/(4 pi^2)) - 1 = {nstr(TN/LIMIT - 1, 5)}')
    if N >= 1000:
        check(c1 and c2 and B4up <= bound and B4lo <= bound and B4min <= -4 and abs(TN/LIMIT - 1) < mpf(10)**(-30),
              f'N = {N}: (B4) inequalities hold, tails below the claimed bound, T_N agrees with the limit to 30 digits')
check(all(rows[i][1] < 60 and rows[i][2] < 80 and rows[i][3] < 40 and abs(rows[i][8]) < 12 and abs(rows[i][9]) < 12 for i in range(1, len(rows))),
      'O-constants of (B2), (B5) stay bounded from N = 10^3 to 10^6')
check(all(abs(rows[i+1][4]) < abs(rows[i][4]) for i in range(len(rows)-1)), '(B3) error decreases with N')

# first N from which the two (B4) ratio inequalities hold without exception up to 30000
last_bad = 0
for N in range(20, 30001):
    k0 = int(floor(XS*N)); J = int(floor(mpf(N)**0.6))
    if k0 + J > N - 2: last_bad = N; continue
    if not (r_mp(N, k0+J) <= exp(-mpf(12*J)/N) and r_mp(N, k0-J-1) >= exp(mpf(12*J)/N)): last_bad = N
print(f'(B4) ratio inequalities hold for every N in ({last_bad}, 30000]')

print('ALL PART-B CHECKS OK' if OK else 'SOME PART-B CHECK FAILED')
sys.exit(0 if OK else 1)
