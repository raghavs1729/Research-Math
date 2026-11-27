#!/usr/bin/env python3
"""Float exploration (no rigor): general degree-one test function
Q = (1 + sum al_n cos th_n) - (B0 + sum be_n cos th_n) F  combined with the
reflection inequality (DENSITY.md, Theorem 4.6).  Any feasible point gives a
valid bound once certified; this script only searches."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")
def load(fn): return [float(l) for l in open(os.path.join(ZDIR, fn)) if l.strip()]

def Jgen(B0, al, be, a, VK, Vt):
    Sba = sum(x*y for x, y in zip(be, a)); Sbb = sum(x*x for x in be); Saa = sum(x*x for x in a)
    Sb2a2 = sum(x*x*y*y for x, y in zip(be, a))
    Sc = sum((p - q - B0*r)**2 for p, q, r in zip(al, be, a))
    return ((1-B0)**2 - (1-B0)*Sba + Sc/2 + Saa*Sbb/4 + Sba**2/2 - 3*Sb2a2/8
            + (B0*B0 + Sbb/2)*VK + (abs(B0) + sum(abs(x) for x in be))**2*Vt)

def bound(B0, al, be, a, VK, Vt, c0):
    if B0 - sum(abs(x) for x in be) < 0: return 9.0
    s1 = 1 - sum(abs(x) for x in al)
    if s1 <= 0: return 9.0
    s2 = max(0.0, 2*c0*B0 - 1 - sum(abs(2*c0*q - p) for p, q in zip(al, be)))
    return Jgen(B0, al, be, a, VK, Vt)/(s1*s1 + s2*s2)

def minimise(a, VK, Vt, c0, B0=1.0):
    N = len(a); al = [0.0]*N; be = [0.0]*N
    best = bound(B0, al, be, a, VK, Vt, c0); step = 0.05
    while step > 1e-7:
        imp = False
        for i in range(-1, 2*N):
            for sg in (1, -1):
                B1, al1, be1 = B0, list(al), list(be)
                if i == -1: B1 += sg*step
                elif i < N: al1[i] += sg*step
                else: be1[i-N] += sg*step
                v = bound(B1, al1, be1, a, VK, Vt, c0)
                if v < best - 1e-15: best, B0, al, be, imp = v, B1, al1, be1, True
            if 0 <= i < N:      # joint move al_i, be_i
                for sg in (1, -1):
                    for rr in (2*c0, 1.0, 0.5):
                        al1, be1 = list(al), list(be)
                        be1[i] += sg*step; al1[i] += sg*step*rr
                        v = bound(B0, al1, be1, a, VK, Vt, c0)
                        if v < best - 1e-15: best, al, be, imp = v, al1, be1, True
        if not imp: step /= 2
    return best, B0, al, be

if __name__ == "__main__":
    for name, zf, s2, VA, VB in (("chi_4", "chi4_zeros.txt", 0.15556797992358593, 0.1696308345, 0.1555736302),
                                  ("chi_3", "chi3_zeros.txt", 0.11322996985747234, 0.1254514113, 0.1132342076)):
        g = load(zf); M = 0.5; T = g[-1] + 0.25; c0 = 1 - s2/2
        amp = [4*M/math.sqrt(M*M + x*x) for x in g]; var = [x*x/2 for x in amp]
        for tier, V in (("A", VA), ("B", VB)):
            for N in (1, 2, 4, 8, 12, 16, 20):
                if 3*g[N-1] >= T: continue
                Tp = T - 2*g[N-1]
                VK = sum(var[i] for i in range(N, len(g)) if g[i] < Tp)
                Vt = V - sum(var[:N]) - VK
                val, B0, al, be = minimise(amp[:N], VK, Vt, c0)
                print("%s tier %s N=%2d bound %.6f dens>=%.6f B0=%.4f sum|al|=%.3f sum|be|=%.3f" %
                      (name, tier, N, val, 1-val, B0, sum(map(abs, al)), sum(map(abs, be))))
