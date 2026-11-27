#!/usr/bin/env python3
"""Referee test of Lemma B with ADVERSARIAL errors e (worst case over e in [-1/2,1/2]^B, per l),
and of Lemma A, for q = 20, 60 and M = 1/2, 0.0025, at N = ceil(4k/tau) (commutator stressed)."""
import math, numpy as np
from math import comb

def fexp(x):
    x = np.asarray(x, float); r = np.zeros_like(x); p = x > 0; r[p] = np.exp(-1.0 / x[p]); return r
def sstep(x):
    x = np.clip(np.asarray(x, float), -1.0, 2.0); a, b = fexp(x), fexp(1 - x); return a / (a + b)

def setup(q, L=2**20):
    eta = math.pi/(4*q); sig = math.pi/(8*q); th1 = 2*math.pi/q - eta + sig; tau = 2*math.sin(th1/2)
    th = 2*np.pi*np.arange(L)/L; d0 = np.minimum(th, 2*np.pi - th)
    u = sstep((d0 - (eta - sig))/sig) * sstep((th1 - d0)/sig)
    uhat = np.real(np.fft.fft(u))/L
    up = 2.0/sig
    lam = lambda s: 1 + (math.pi/math.sqrt(3))*(s/tau + up)
    z = np.exp(1j*th) - 1
    Lam = [np.abs(np.fft.fft(z**s*u)/L).sum() for s in range(0, 16)]
    return tau, uhat, lam, Lam, L

for q in (20, 60):
    tau, uhat, lam, Lam, L = setup(q)
    uh = lambda p: uhat[np.mod(p, L)]
    print(f"q={q} tau={tau:.4f}")
    for s in range(0, 13, 3):
        print(f"   Lemma A s={s}: Lambda_s/(tau^s lam(s)) = {Lam[s]/(tau**s*lam(s)):.4f}")
    for M in (0.5, 0.0025):
        a = lambda x: 2*M*M/(M*M + np.asarray(x, float)**2)
        for k in (0, 2, 4, 6, 8, 10):
            N = max(2, math.ceil(4*k/tau))
            h = np.array([(-1)**i*comb(k, i) for i in range(k+1)], float)
            ms = np.arange(N, 4*N)
            ls = np.arange(-12*N, 12*N+1)
            C = np.zeros((len(ms), len(ls)))
            for i in range(k+1):
                C += h[i]*a(ms+i)[:, None]*(uh(ls[None, :] - ms[:, None] - i) + uh(ls[None, :] + ms[:, None] + i))
            worst = 0.5*np.abs(C).sum(axis=0).max()          # max_l max_e |hat v(l)|
            bd1 = a(N)*sum(comb(k, r)*math.factorial(r+1)*N**(-r)*Lam[k-r] for r in range(k+1))
            bd2 = (16/9)*a(N)*tau**k*lam(k)
            # pure-filter quantity (no commutator): a(N) * Lambda_k / 2 * 2
            print(f"   M={M:<6} k={k:2d} N={N:4d}: worst|v^|/a(N)={worst/a(N):.3e}  Leibniz bd={bd1/a(N):.3e}"
                  f"  (16/9)tau^k lam bd={bd2/a(N):.3e}  ratio worst/Leibniz={worst/bd1:.3f}")
