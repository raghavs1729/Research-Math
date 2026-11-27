"""Sanity checks of the identities used in the proof of TOY.md Theorem 2.3.

AI-produced (Anthropic's Claude), 2026-11-26.

For a single-hole symmetric Toeplitz sequence (periods P_k = q_1...q_k, symmetric patterns
v_k with v_k(0)=0), fix a level i+1 and write F = F^{(i+1)} + F_{>=i+2}, where F^{(i+1)} is the
F-function of the P_{i+1}-periodic truncation (positions in P_{i+1}Z set to 0). Checks:
 (a) F_{>=i+2} is 2pi/P_{i+1}-periodic (identical on all level-(i+1) sub-arcs);
 (b) on each sub-arc J (t = 2pi(a'+u')/P_{i+1}), F^{(i+1)} = 1 + 2 B(a') e^{-eps' u'} exactly,
     and B(e) = B(e-1) e^{-eps'} + eps' V_{i+1}(e) for consecutive sub-arcs inside a level-i arc
     (e = a' mod q_{i+1} != 0);
 (c) ||F_{>=i+2}||_2^2 <= (8 pi^2/6) M^2 / P_{i+1}^2   (normalised L^2).
Uses the closed form of toy_toeplitz.py (checked there against direct Fourier sums).
Output: toy_toeplitz2_results.txt
"""
import numpy as np
from toy_toeplitz import closed_form, V_of, sym_pattern

def main():
    rng = np.random.default_rng(7)
    M = 0.5
    out = []
    for qs in ([2, 3, 2, 5, 3, 2, 2, 3], [3, 7, 2, 11, 2, 3], [2, 2, 31, 2, 3, 2]):
        pats = [sym_pattern(q, rng) for q in qs]
        P = [1]
        for q in qs:
            P.append(P[-1]*q)
        for i in range(0, len(qs) - 2):
            Pn = P[i+1]; q = qs[i]  # sub-arcs of level i+1 inside arcs of level i; q = q_{i+1}
            epsp = 2*np.pi*M/Pn
            nu = 512
            u = (np.arange(nu)+0.5)/nu
            # (a) fine part on several sub-arcs
            fines = []; Bs = []
            a = rng.integers(0, P[i])           # a level-i arc
            for e in range(q):
                ap = a*q + e                      # sub-arc index at level i+1
                t = 2*np.pi*(ap+u)/Pn
                F = closed_form(t, M, qs, pats)
                Fc = closed_form(t, M, qs[:i+1], pats[:i+1])
                fines.append(F - Fc)
                B = 0.5*(Fc - 1.0)*np.exp(epsp*u)
                Bs.append((B.mean(), B.max()-B.min()))
            fines = np.array(fines)
            dev_a = np.max(np.abs(fines - fines[0]))
            dev_b1 = max(b[1] for b in Bs)
            V = V_of(pats[i])                   # V_{i+1}
            rec = max(abs(Bs[e][0] - (Bs[e-1][0]*np.exp(-epsp) + epsp*V[e])) for e in range(1, q))
            L2 = np.mean(fines[0]**2)
            bound = (8*np.pi**2/6)*M**2/Pn**2
            out.append(f"qs={qs} level i+1={i+1} (q_{i+1}={q}, P={Pn}): (a) max diff of fine part "
                       f"across sub-arcs {dev_a:.1e}; (b) B constant on sub-arcs to {dev_b1:.1e}, "
                       f"recursion error {rec:.1e}; (c) ||fine||^2={L2:.2e} <= bound {bound:.2e}: {L2 <= bound}")
    txt = "\n".join(out)
    print(txt)
    with open("toy_toeplitz2_results.txt", "w") as f:
        f.write("toy_toeplitz2.py results (AI-produced, 2026-11-26)\n" + txt + "\n")

if __name__ == "__main__":
    main()
