"""Tests of Theorem 9.16 on finite combs (Lemma 9.0 closed form + coefficient synthesis)."""
import numpy as np
from mpmath import mp, mpf, findroot, exp, log, sqrt, pi as mpi

def rho(t, N, M):
    h = 2*np.pi/N
    s = np.mod(t, h)
    return h*np.exp(-M*s)/(1-np.exp(-M*h)) - 1/M

def F_comb(t, Ns, bs, M):
    return 1 + 2*M*sum(b*rho(t, N, M) for N, b in zip(Ns, bs))

def F_coeff(t, mfun, M, nmax):
    n = np.arange(1, nmax+1)
    m = mfun(n)
    c = 2*M*m/(M+1j*n)
    # Fejer-free direct partial sum (pointwise, away from jumps)
    out = np.ones_like(t)
    for chunk in range(0, nmax, 20000):
        nn = n[chunk:chunk+20000]; cc = c[chunk:chunk+20000]
        out += 2*np.real(np.exp(1j*np.outer(t, nn)) @ cc)
    return out

def g_inf(Ns, bs):
    # g = sum_k b_k z^N/(1-z^N); each term -> -1 at infinity
    return -sum(bs)

# (A) N=(1,4), b=(1,-1): M = (2/pi) log(1+sqrt2)
M = float(2/mpi*log(1+sqrt(2)))
print("M(1,4) =", repr(M))
t = np.linspace(np.pi/2+1e-3, np.pi-1e-3, 400)
F = F_comb(t, [1, 4], [1, -1], M)
print(" (1,4): max|F-1| on [pi/2,pi):", np.max(np.abs(F-1)), " g(inf)=", g_inf([1,4],[1,-1]), "=> c=", 1+2*g_inf([1,4],[1,-1]))
Fc = F_coeff(t[::20], lambda n: 1.0 - (n % 4 == 0), M, 400000)
print("   coefficient synthesis (n<=4e5) max|F-1|:", np.max(np.abs(Fc-1)))

# (B) N=(1,3), b=(2,-1): need A_1 = 2 A_0 on last cell -> (1+x+x^2)/3 = 2, x=e^u, u=2 pi M/3
x = (-1+np.sqrt(21))/2; u = np.log(x); M = 3*u/(2*np.pi)
t = np.linspace(4*np.pi/3+1e-3, 2*np.pi-1e-3, 400)
F = F_comb(t, [1, 3], [2, -1], M)
print("(1,3),(2,-1): M=", M, " max|F+1| on (4pi/3,2pi):", np.max(np.abs(F+1)), " g(inf)=", g_inf([1,3],[2,-1]),
      "=> c=", 1+2*g_inf([1,3],[2,-1]))
Fc = F_coeff(t[::20], lambda n: 2.0 - (n % 3 == 0), M, 400000)
print("   coefficient synthesis max|F+1|:", np.max(np.abs(Fc+1)))

# (C) N=(1,3), b=(3,-1): A_1 = 3/2 A_0 -> level c = 1-2*2 = -3
x = (-1+np.sqrt(1+4*8))/2; u = np.log(x); M = 3*u/(2*np.pi)   # A1/A0 = 3
F = F_comb(t, [1, 3], [3, -1], M)
print("(1,3),(3,-1): M=", M, " max|F-(-3)|:", np.max(np.abs(F+3)), " g(inf)=", g_inf([1,3],[3,-1]),
      "=> c=", 1+2*g_inf([1,3],[3,-1]))
# (D) pure comb q|n: w constant -1/2 off comb but F not constant
F = F_comb(t, [3], [1], 0.5)
print("comb 3|n, M=1/2: F on (4pi/3,2pi) range", F.min(), F.max(), "(not constant; kappa=-1 would be forced)")
