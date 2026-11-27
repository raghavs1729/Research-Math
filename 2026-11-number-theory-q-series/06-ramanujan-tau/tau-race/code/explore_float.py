import sys, math, time
import numpy as np
from tau_table import tau_list
N = int(sys.argv[1])
t = tau_list(N+1)
sieve = np.ones(N+1, dtype=bool); sieve[:2] = False
for i in range(2, int(N**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
primes = np.nonzero(sieve)[0]
lam = np.array([t[p]/p**5.5 for p in primes], dtype=np.longdouble)
P = primes.astype(np.longdouble); lg = np.log(P)
print("max |lam|", float(np.max(np.abs(lam))))
# Theta(p) = sum_{r<=p} lam log r ; Theta_1(p) = sum lam(r) log r (p - r) = p*Theta(p) - sum lam log r * r
Th = np.cumsum(lam*lg); T1 = P*Th - np.cumsum(lam*lg*P)
G = T1 / P**1.5
D0 = np.cumsum(lam); D1 = np.cumsum(lam*P); R = D0 - D1/P
print("first values p, Theta_1, R:")
for i in range(8): print(int(primes[i]), float(T1[i]), float(R[i]))
neg = primes[T1 <= 0]; print("primes with Theta_1(p)<=0:", neg[:20], len(neg))
negR = primes[R <= 0]; print("primes with R(p)<=0:", negR[:20], len(negR))
for lo in [10, 100, 10**3, 10**4, 10**5, 10**6]:
    m = primes >= lo
    if m.any():
        j = np.argmin(G[m]); print("inf G on p>=%d: %.6f at p=%d ; max %.6f" % (lo, float(G[m][j]), int(primes[m][j]), float(np.max(G[m]))))
        Nn = R[m]*lg[m]/np.sqrt(P[m]); j = np.argmin(Nn); print("   inf N=R log/sqrt: %.6f at %d" % (float(Nn[j]), int(primes[m][j])))
# sym2: Theta_S(y) = sum lam(p^2) log p = sum (lam^2-1) log p
ThS = np.cumsum((lam*lam-1)*lg)
for lo in [10,100,10**3,10**4,10**5,10**6]:
    m = primes>=lo
    if m.any(): print("sup |Theta_S(y)|/sqrt(y) y>=%d: %.4f ; sup|Theta(y)|/sqrt y %.4f; (y-theta)/sqrt y sup %.4f" % (lo, float(np.max(np.abs(ThS[m])/np.sqrt(P[m]))), float(np.max(np.abs(Th[m])/np.sqrt(P[m]))), float(np.max((P[m]-np.cumsum(lg)[m]+lg[m])/np.sqrt(P[m])))))
# unweighted D0 and KK sum
KK = np.cumsum(lam/np.sqrt(P))
print("D0<=0 primes:", primes[D0<=0][:30], int(np.sum(D0<=0)))
print("KK<=0 primes:", primes[KK<=0][:30], int(np.sum(KK<=0)), " KK at end", float(KK[-1]), 0.5*math.log(math.log(N)))
