#!/usr/bin/env python3
"""
audit2 / a2: the WZ pairs of T4, T11, T5, T7 as actual numbers, in exact rational arithmetic.

Independent of ../proofs and of the ratio files: F(n,k) is computed here from a closed form
(Pochhammer products written out by the referee from the statements in RESULTS.md), G = R*F with R
typed in again below from RESULTS.md / the certificate text (NUM and P are read from the .gp files,
they are too long to retype; a1 checks them symbolically).  Checked on grids:

  * the closed form has the ratios sigma, rho printed in RESULTS.md,
  * the WZ equation  F(n+1,k) - F(n,k) = G(n,k+1) - G(n,k)  term by term (all cases of the triangle),
  * the boundary data (F(0,k), G(n,0) versus the summand of the conjecture),
  * the telescoped finite identity, including identity (A) of Theorem 11 and the finite identity of T7,
  * the rational inequalities used in the limit steps of T4 and T5 (as exact polynomial identities /
    coefficient positivity, via python-flint).
"""
import os, re, sys
from fractions import Fraction as Fr
from math import comb, factorial
import flint

HERE = os.path.dirname(os.path.abspath(__file__))
def gp_poly(fn, name):
    txt = open(os.path.join(HERE, '..', 'certificates', fn)).read()
    m = re.search(r'^' + name + r'\s*=\s*(.*?);', txt, flags=re.S | re.M)
    return m.group(1).replace('^', '**')

def poch(a, m):
    r = Fr(1)
    for i in range(m):
        r *= a + i
    return r

OK = True
def check(cond, msg):
    global OK
    OK &= bool(cond)
    print(('  ok   ' if cond else '  FAIL ') + msg)

# ============================================================================ T4
print('T4  (pair n*h, Q*h)')
def h4(n, k):   # closed form, h(0,0) = 1
    return Fr(1, 2**(2*k + 5*n)) * poch(Fr(1, 2), 3*n) * poch(Fr(1, 2), 2*k + 3*n) / (
        poch(Fr(1, 2), n) * factorial(k + n) * factorial(n) * factorial(k + 3*n))
def sg4(n, k):  return Fr((4*k+6*n+1)*(4*k+6*n+3), 16*(k+n+1)*(k+3*n+1))
def rho4(n, k): return Fr(3*(6*n+1)*(6*n+5)*(4*k+6*n+1)*(4*k+6*n+3)*(4*k+6*n+5),
                          1024*(n+1)*(k+n+1)*(k+3*n+1)*(k+3*n+2)*(k+3*n+3))
def Q4(n, k):   return Fr(64*k**3+(160*n+144)*k**2+(240*n**2+432*n+92)*k+360*n**3+612*n**2+230*n+15,
                          64*(k+3*n+1)*(k+3*n+2))
def t4(n):      return Fr((360*n**3+612*n**2+230*n+15)*comb(3*n, n)*comb(6*n, 3*n)**2, (3*n+1)*(3*n+2)*2**(15*n))
G = range(0, 26)
check(all(h4(n, k+1) == sg4(n, k)*h4(n, k) and h4(n+1, k) == rho4(n, k)*h4(n, k) for n in G for k in G),
      'closed form of h has the ratios sigma, rho of RESULTS.md (0<=n,k<=25)')
check(all((n+1)*h4(n+1, k) - n*h4(n, k) == Q4(n, k+1)*h4(n, k+1) - Q4(n, k)*h4(n, k) for n in G for k in G),
      'WZ equation for F = n h, G = Q h, term by term (0<=n,k<=25)')
check(all(h4(n, 0) == Fr(comb(3*n, n)*comb(6*n, 3*n)**2, 2**(15*n)) for n in range(40)), 'h(n,0) = C(3n,n) C(6n,3n)^2 / 2^(15n)')
check(all(h4(0, k) == poch(Fr(1, 4), k)*poch(Fr(3, 4), k)/factorial(k)**2 for k in range(40)), 'h(0,k) = (1/4)_k (3/4)_k / k!^2')
check(all(Q4(n, 0)*h4(n, 0) == t4(n)/64 for n in range(40)), 'Q(n,0) h(n,0) = t(n)/64')
check(all(sum(Q4(n, 0)*h4(n, 0) for n in range(N)) == sum(Q4(n, K)*h4(n, K) for n in range(N)) - N*sum(h4(N, k) for k in range(K))
          for N in range(0, 13) for K in range(0, 13)),
      'finite identity  sum_{n<N} Q(n,0)h(n,0) = sum_{n<N} Q(n,K)h(n,K) - N sum_{k<K} h(N,k)  (N,K<=12)')
# inequalities, symbolically
ctx = flint.fmpq_mpoly_ctx.get(['n', 'k'], 'lex'); n, k = ctx.gens()
def nonneg_coeffs(p): return all(c >= 0 for c in p.coeffs())
num = 16*(k+n+1)*(k+3*n+1) - (4*k+6*n+1)*(4*k+6*n+3) - 12*(n+1)*(k+n+1)
check(num == 4*n*k + 4*k + 16*n + 1, '1 - sigma - (3/4)(n+1)/(k+3n+1) = (4nk+4k+16n+1)/(16(k+n+1)(k+3n+1))  [>= 0]')
numQ = 64*k**3+(160*n+144)*k**2+(240*n**2+432*n+92)*k+360*n**3+612*n**2+230*n+15
d = (k+3*n+3)*64*(k+3*n+1)*(k+3*n+2) - numQ
check(nonneg_coeffs(d) and nonneg_coeffs(numQ), '0 <= Q(n,k) <= k+3n+3  (numerators have nonnegative coefficients)')
d = 27*1024*(n+1)**2*(3*n+1)*(3*n+2)*(3*n+3) - 32*3*(6*n+1)**2*(6*n+5)**2*(6*n+3)
check(nonneg_coeffs(d), 'rho(n,0) <= 27/32  (difference of cross products has nonnegative coefficients): %s' % d)
check(all(Fr(3*(nn+1), 4) <= 3*nn+1 for nn in range(0, 1000)), 'c = 3(n+1)/4 <= m = 3n+1, so every factor 1 - c/(j+m) is in [0,1)')
check(all(1 + Fr(4*(3*nn+1), 3*nn-1) <= Fr(33, 5) for nn in range(2, 2000)), '1 + m/(c-1) = 1 + 4(3n+1)/(3n-1) <= 33/5 for n >= 2')

# ============================================================================ T11
print('T11 (triangle pair a, gam)')
NUM_SRC = gp_poly('T11.gp', 'NUM')
def NUM(n, k): return Fr(eval(NUM_SRC, {'__builtins__': {}}, {'n': Fr(n), 'k': Fr(k)}))
def C11(n):    # = 2^(6-8n) Gamma(n+3/4)^2 Gamma(n+1/4) / (sqrt2 pi Gamma(n-1/4) Gamma(2n-1));
               # Gamma(3/4)^2 Gamma(1/4)/Gamma(-1/4) = -Gamma(1/4)Gamma(3/4)/4 = -pi sqrt2/4
    return Fr(2)**(6-8*n) * poch(Fr(3, 4), n)**2 * poch(Fr(1, 4), n) * Fr(-1, 4) / (poch(Fr(-1, 4), n) * factorial(2*n-2))
def a11(n, k):
    if k >= n: return Fr(0)
    # Gamma(4k+2n+1/2)/Gamma(4k-2n+9/2) = (4k-2n+9/2)_(4n-4)
    return C11(n)*(2*k+1)*poch(Fr(8*k-4*n+9, 2), 4*n-4) / (factorial(n-k-1)**2 * factorial(n+k)**2)
def R11(n, k): return NUM(n, k) / (131072*n*(2*n-1)*(2*n+1)*(4*n-1)*(k-n)**2*(2*k+1))
def gam11(n, k):
    if k < n: return R11(n, k)*a11(n, k)
    if k == n: return -a11(n+1, n)
    return Fr(0)
def sg11(n, k):
    return Fr(2*k+3, 2*k+1)*Fr(n-k-1, n+k+1)**2*poch(Fr(8*k+4*n+1, 2), 4)/poch(Fr(8*k-4*n+9, 2), 4)
def rho11(n, k):
    return (Fr(4*n+3, 4)**2*Fr(4*n+1, 4)*Fr(8*k+4*n+1, 2)*Fr(8*k+4*n+3, 2)*Fr(8*k-4*n+5, 2)*Fr(8*k-4*n+7, 2)
            / (512*n*Fr(4*n-1, 4)*(2*n-1)*(n-k)**2*(n+k+1)**2))
def t11(n): return Fr((1344*n**3+944*n**2+156*n+9)*comb(2*n, n)**2*comb(4*n, 2*n)**3, (2*n+1)*2**(22*n))
NMAX = 45
check(C11(1) == Fr(9, 256) and a11(1, 0) == t11(0)/256, 'C_1 = 9/256 and a(1,0) = t(0)/256 = 9/256')
check(all(a11(n, k) > 0 for n in range(1, NMAX) for k in range(n)), 'a(n,k) > 0 for 0 <= k < n <= 44')
check(all(a11(n, k+1) == sg11(n, k)*a11(n, k) for n in range(1, NMAX) for k in range(n-1))
      and all(a11(n+1, k) == rho11(n, k)*a11(n, k) for n in range(1, NMAX) for k in range(n)),
      'closed form of a has the ratios sigma, rho of RESULTS.md on the triangle (n <= 44)')
check(all(a11(n+1, k) - a11(n, k) == gam11(n, k+1) - gam11(n, k) for n in range(1, NMAX) for k in range(0, NMAX+5)),
      'WZ equation a(n+1,k)-a(n,k) = gam(n,k+1)-gam(n,k) for all 1<=n<=44, 0<=k<=49 (cases k<=n-2, k=n-1, k=n, k>n)')
check(all(gam11(n, 0) == -t11(n)/256 for n in range(1, NMAX)), 'gam(n,0) = -t(n)/256 for 1 <= n <= 44')
TN = lambda N: sum(a11(N, k) for k in range(N))
check(all(256*TN(N) == sum(t11(n) for n in range(N)) for N in range(1, NMAX)), 'identity (A): sum_{n<N} t(n) = 256 T_N for N = 1..44')
print('      256*T_N, N=1..6:', [float(256*TN(N)) for N in range(1, 7)])

# ============================================================================ T5
print('T5  (pair psi, gamma = R psi)')
def sg5(n, k):  return Fr((2*k+n+1)*(2*k+n+2), (2*k+2*n+3)*(2*k+4*n+3))
def rho5(n, k): return Fr(-(n+1)*(3*n+1)*(3*n+2)*(3*n+3)*(2*k+n+1), (4*n+1)*(4*n+3)*(2*k+2*n+3)*(2*k+4*n+3)*(2*k+4*n+5))
def R5(n, k):   return Fr((16*n+12)*k**2+(66*n**2+86*n+28)*k+77*n**3+145*n**2+88*n+17, (4*n+1)*(4*n+3)*(2*k+4*n+3))
psi_cache = {}
def psi5(n, k):
    if (n, k) not in psi_cache:
        if k > 0:   v = psi5(n, k-1)*sg5(n, k-1)
        elif n > 0: v = psi5(n-1, 0)*rho5(n-1, 0)
        else:       v = Fr(1)
        psi_cache[(n, k)] = v
    return psi_cache[(n, k)]
def t5(k): return Fr((77*k**3-86*k**2+29*k-3)*(-32)**k, k**3*(3*k-1)*(3*k-2)*comb(3*k, k)*comb(4*k, k)**2)
G = range(0, 26)
check(all(psi5(n+1, k) == rho5(n, k)*psi5(n, k) for n in G for k in G), 'psi is well defined by its ratios (rho consistent off the axis k = 0)')
check(all(psi5(n+1, k) - psi5(n, k) == R5(n, k+1)*psi5(n, k+1) - R5(n, k)*psi5(n, k) for n in G for k in G), 'WZ equation term by term (0<=n,k<=25)')
check(all(psi5(0, k) == Fr(4**k, (2*k+1)**2*comb(2*k, k)) for k in range(60)), 'psi(0,k) = 4^k / ((2k+1)^2 C(2k,k))')
check(all(R5(n, 0)*psi5(n, 0) == -t5(n+1)/3 for n in range(40)), 'gamma(n,0) = -t(n+1)/3')
num = (2*k+2*n+3)*(2*k+4*n+3) - (2*k+n+1)*(2*k+n+2)
check(2*num - (4*n+3)*(2*k+2*n+3) == 8*n*k + 6*k + 6*n**2 + 12*n + 5,
      '1 - sigma - (n+3/4)/(k+2n+3/2) = (8nk+6k+6n^2+12n+5)/(2(2k+2n+3)(2k+4n+3))  [>= 0]')
d = (4*n+1)*(4*n+3)*(2*n+3)*(4*n+3)*(4*n+5) - 8*(n+1)**2*(3*n+1)*(3*n+2)*(3*n+3)
check(nonneg_coeffs(d), '|rho(n,0)| <= 1/8')
numR = (16*n+12)*k**2+(66*n**2+86*n+28)*k+77*n**3+145*n**2+88*n+17
check(nonneg_coeffs(numR) and nonneg_coeffs((2*k+2*n+2)*(4*n+1)*(4*n+3)*(2*k+4*n+3) - numR), '0 <= R(n,k) <= 2k+2n+2')
check(all(1 + Fr(8*N+6, 4*N-1) <= 6 for N in range(1, 2000)), '1 + m/(c-1) = 1 + (2N+3/2)/(N-1/4) <= 6 for N >= 1')

# ============================================================================ T7
print('T7  (triangle pair Phi, Gam; values are multiplied by pi^2)')
P_SRC = gp_poly('T7.gp', 'P')
def P7(n, k): return Fr(eval(P_SRC, {'__builtins__': {}}, {'n': Fr(n), 'k': Fr(k)}))
def Phi7(n, k):     # pi^2 * Phi(n,k);  Gamma(m+1/2) = sqrt(pi) (1/2)_m
    if k >= n: return Fr(0)
    return (Fr((-1)**(k+1)*(4*k+4*n+3)*factorial(k)*factorial(4*n)*factorial(2*n), 2**(6*n-1))
            / (factorial(n-k-1)*(2*k+2*n+1)*(2*k+2*n+2)*poch(Fr(1, 2), k+2*n+1)*poch(Fr(1, 2), k+3*n+1)))
def R7(n, k): return P7(n, k) / (3*(6*n+1)*(6*n+5)*(k-n)*(2*k+4*n+3)*(2*k+6*n+3)*(2*k+6*n+5)*(4*k+4*n+3))
def Gam7(n, k):
    if k < n: return R7(n, k)*Phi7(n, k)
    if k == n: return -Phi7(n+1, n)
    return Fr(0)
def t7(k): return Fr(comb(2*k, k)*16**k*(1656*k**4-2604*k**3+1158*k**2-229*k+15), k**2*(2*k-1)**2*(6*k-1)*(6*k-5)*comb(6*k, 3*k)*comb(3*k, k))
NM = 40
check(all(Phi7(n+1, k) - Phi7(n, k) == Gam7(n, k+1) - Gam7(n, k) for n in range(0, NM) for k in range(0, NM+5)),
      'WZ equation for all 0<=n<=39, 0<=k<=44 (cases k<=n-2, k=n-1, k=n, k>n)')
check(all(Gam7(n, 0) == -t7(n+1)/6 for n in range(0, NM)), 'pi^2 Gam(n,0) = -t(n+1)/6 for 0 <= n <= 39 (n = 0 uses Gam(0,0) = -Phi(1,0))')
check(all(sum(t7(m) for m in range(1, N+1)) == 6*sum(Phi7(N, k) for k in range(N)) for N in range(1, NM)),
      'finite identity sum_{m<=N} t(m) = 6 pi^2 sum_{k<N} Phi(N,k), N = 1..39')
check(all(abs(Phi7(N, k+1)) <= abs(Phi7(N, k))/24 for N in range(2, 200) for k in range(N-1)), '|Phi(N,k+1)/Phi(N,k)| <= 1/24 on the support (N < 200)')
check(all(0 < Phi7(N+1, 0)/Phi7(N, 0) <= Fr(1, 6) for N in range(1, 400)), '0 < Phi(N+1,0)/Phi(N,0) <= 1/6 (N < 400)')
print('      partial sums of T7, N = 1,2,3,6,12:', [float(sum(t7(m) for m in range(1, N+1))) for N in (1, 2, 3, 6, 12)])


# ---- T7: the two bounds of the limit step, symbolically (they are only asserted in RESULTS.md)
ctx3 = flint.fmpq_mpoly_ctx.get(['k', 'v', 'm'], 'lex'); kk, vv, mm = ctx3.gens()
half = flint.fmpq(1, 2)
nn = kk + 2 + vv                       # support of the ratio: 0 <= k <= n-2
numS = (4*kk+4*nn+7)*(kk+1)*(nn-kk-1)*(2*kk+2*nn+1)*(2*kk+2*nn+2)
denS = (4*kk+4*nn+3)*(2*kk+2*nn+3)*(2*kk+2*nn+4)*(kk+2*nn+3*half)*(kk+3*nn+3*half)
check(nonneg_coeffs(denS - 24*numS), '|Phi(n,k+1)/Phi(n,k)| <= 1/24 for all 0 <= k <= n-2 (n = k+2+v: den - 24 num has nonnegative coefficients)')
nn = mm + 1                            # n >= 1
numR = (4*nn+7)*(4*nn+1)*(4*nn+2)*(4*nn+3)*(4*nn+4)*(2*nn+1)**2*(2*nn+2)**2
denR = 64*(4*nn+3)*nn*(2*nn+3)*(2*nn+4)*(2*nn+3*half)*(2*nn+5*half)*(3*nn+3*half)*(3*nn+5*half)*(3*nn+7*half)
check(nonneg_coeffs(denR - 6*numR), 'Phi(n+1,0)/Phi(n,0) = rho(n,0) <= 1/6 for all n >= 1 (n = m+1: den - 6 num has nonnegative coefficients)')
check(all(Phi7(n+1, 0)/Phi7(n, 0) == Fr((4*n+7)*(4*n+1)*(4*n+2)*(4*n+3)*(4*n+4)*(2*n+1)**2*(2*n+2)**2 * 2**5,
                                        64*(4*n+3)*n*(2*n+3)*(2*n+4)*(4*n+3)*(4*n+5)*(6*n+3)*(6*n+5)*(6*n+7)) for n in range(1, 60)),
      'rho(n,0) of the certificate file is the ratio Phi(n+1,0)/Phi(n,0) of the closed form (n < 60)')

# ============================================================================ T2 (inequalities only; WZ in a1)
print('T2  (inequalities of the limit step)')
lhs = lambda n_, k_: (2*k_+5)*(2*k_+n_+3)*(3*k_+n_+1)*(3*k_+n_+2)*(3*k_+n_+3)
rhs = lambda n_, k_: (2*k_+n_+1)*(2*k_+2*n_+3)*(3*k_+2*n_+3)*(3*k_+2*n_+4)*(3*k_+2*n_+5)
check(nonneg_coeffs(rhs(n+1, k) - lhs(n+1, k)), 'sigma(n,k) <= (2k+1)/(2k+5) for n >= 1')
check(nonneg_coeffs(9*(2*n+1)*(n+1)*(2*n+3)*(2*n+3)*(2*n+4) - 8*(n+1)*(4*n+3)*(4*n+5)*(n+2)*(n+1)), 'rho(n,0) <= 1/8')
check(nonneg_coeffs((k+2)*9*(2*n+1)*(2*k+n+1) - (3*k+2*n+2)*(6*k+8*n+5)), 'R(n,k) <= k+2')

print('ALL EXACT CHECKS OK' if OK else 'SOME EXACT CHECK FAILED')
sys.exit(0 if OK else 1)
