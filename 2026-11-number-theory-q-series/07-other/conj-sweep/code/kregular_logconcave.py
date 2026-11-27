# arXiv:2406.10987 (Craig? / Heim-Neuhauser et al.), Conjecture: N_k = N_infinity (= 26) for all k >= 30, where
# N_k = least N such that the k-regular partition function p_k(n) is log-concave at every n >= N
# (log-concave at n: p_k(n)^2 >= p_k(n-1) p_k(n+1)).   Test for k in [KMIN, KMAX], n <= NMAX.
import sys
from flint import fmpz_poly as Z
KMIN, KMAX, NMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
D = NMAX + 2
p = [0]*D; p[0] = 1
for n in range(1, D):                      # pentagonal recurrence
    s = 0; k = 1
    while True:
        g1 = k*(3*k-1)//2; g2 = k*(3*k+1)//2
        if g1 > n: break
        s += (-1)**(k+1) * p[n-g1]
        if g2 <= n: s += (-1)**(k+1) * p[n-g2]
        k += 1
    p[n] = s
P = Z(p)
def fk(k):
    c = [0]*D; j = 0
    while True:
        e1 = k*j*(3*j-1)//2; e2 = k*j*(3*j+1)//2
        if e1 >= D: break
        c[e1] += (-1)**j
        if j and e2 < D: c[e2] += (-1)**j
        j += 1
    return Z(c)
def Nk(seq, nmax):
    last_fail = 0
    for n in range(1, nmax):
        if seq[n]*seq[n] < seq[n-1]*seq[n+1]: last_fail = n
    return last_fail + 1
print('N_infinity on n <=', NMAX, ':', Nk(p, NMAX))
res = {}
for k in range(KMIN, KMAX+1):
    c = [int(v) for v in (P*fk(k)).coeffs()[:D]]
    res[k] = Nk(c, NMAX)
print('N_k (computed on n <= %d) for k = %d..%d:' % (NMAX, KMIN, min(KMAX, KMIN+25)), [res[k] for k in range(KMIN, min(KMAX, KMIN+25)+1)])
print('k >= 30 with N_k != 26:', [(k, res[k]) for k in res if k >= 30 and res[k] != 26])
