#!/usr/bin/env python3
"""Float exploration (no rigor): reflection bound of DENSITY.md section 4bis,
alone and combined with the non-resonance test function of GRH-ONLY.md Thm 6.1."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "grhonly"))
from nonresonance import J_value
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")

def load(fn):
    return [float(l) for l in open(os.path.join(ZDIR, fn)) if l.strip()]

def refl_cheb(V, s2):
    best = (1, None)
    for i in range(0, 4001):
        c = 0.8 + i * 1e-4
        B = (V + (1 - c) ** 2) / (c * c + (2 - s2 - c) ** 2)
        if B < best[0]:
            best = (B, c)
    return best

def bound(k0, b, a, VK, Vt, s2):
    t = k0 * (1 - sum(b)) * (2 - s2) - 1
    J = J_value(k0, b, a, VK, Vt)
    return J / (1 + t * t) if t > 0 else J

def minimise(a, VK, Vt, s2, iters=600):
    N = len(a)
    k0, b = 1.0, [0.0] * N
    best = bound(k0, b, a, VK, Vt, s2)
    step = 0.1
    for it in range(iters):
        imp = False
        for i in range(-1, N):
            for sg in (1, -1):
                k1, b1 = k0, list(b)
                if i == -1:
                    k1 = k0 + sg * step
                    if k1 <= 0: continue
                else:
                    b1[i] += sg * step
                    if b1[i] < 0 or sum(b1) > 1 - 1e-9: continue
                v = bound(k1, b1, a, VK, Vt, s2)
                if v < best - 1e-15:
                    best, k0, b, imp = v, k1, b1, True
        if not imp:
            step /= 2
            if step < 1e-8: break
    return best, k0, b

for name, zf, s2, VA, VB in (("chi_4", "chi4_zeros.txt", 0.15556797992358593, 0.1696308345, 0.1555736302),
                              ("chi_3", "chi3_zeros.txt", 0.11322996985747234, 0.1254514113, 0.1132342076)):
    g = load(zf); M = 0.5
    T = g[-1] + 0.25
    amp = [4 * M / math.sqrt(M * M + x * x) for x in g]
    var = [x * x / 2 for x in amp]
    print(name, "zeros", len(g), "T", T)
    for tier, V in (("A", VA), ("B", VB)):
        B, c = refl_cheb(V, s2)
        print("  tier", tier, "Cantelli", V / (1 + V), " reflection-Chebyshev", B, "at c =", c)
        for N in (1, 2, 4, 6, 8, 10, 12, 16, 20):
            if 3 * g[N - 1] >= T: continue
            Tp = T - 2 * g[N - 1]
            VK = sum(var[i] for i in range(N, len(g)) if g[i] < Tp)
            Vt = V - sum(var[:N]) - VK
            val, k0, b = minimise(amp[:N], VK, Vt, s2)
            t = k0 * (1 - sum(b)) * (2 - s2) - 1
            print("    N=%2d  bound %.6f  density >= %.6f   k0=%.4f sum b=%.4f t=%.4f" % (N, val, 1 - val, k0, sum(b), t))
