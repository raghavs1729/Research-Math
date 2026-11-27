"""Numerical checks for TOY.md section 2.2 (single-hole Toeplitz sequences).

AI-produced (Anthropic's Claude), 2026-11-26.

Setting: q_1, q_2, ... >= 2, P_k = q_1...q_k, symmetric 0/1 patterns v_k on Z/q_k with
v_k(0) = 0; m_n = v_{k(n)+1}(n/P_{k(n)} mod q_{k(n)+1}), k(n) = max{j : P_j | n}.
F_M(t) = 1 + sum_{n != 0} m_|n| 2M/(M+in) e^{int}.

Checks:
 (1) the multilevel closed form
       F = 1 + sum_k (2/q_k) sum_sigma V_k(sigma) phi_{k-1}({y_{k-1} - sigma/q_k}),
       y_j = {P_j t/2pi}, eps_j = 2 pi M/P_j, phi_j(y) = eps_j e^{-eps_j y}/(1-e^{-eps_j}) - 1,
       V_k(sigma) = sum_r v_k(r) e^{2 pi i r sigma/q_k},
     against direct Fourier summation (K-level truncation = periodic sequence, compared
     with a long partial sum away from the atoms);
 (2) the local first-order structure on a level-i arc:
       F(u) = C + eps_i * X(u) + O(eps_i^2 ...),  X(u) = sum_l g_{i+l}(e_l)/Q_{i+l},
       g_k(e) = e - 2 T_k(e) + (1-2B-1)... (see TOY.md), increments g_k(e)-g_k(e-1) = 1+2V_k(e)
     (we check: (F(u) - F(0))/eps_i vs X(u) - X(0), after removing the exact coarse
      exponential, for the valuation family and for random patterns);
 (3) the valuation family: the local law of (F - F(left end))/eps_i on a level-i arc is
     close to uniform on an interval of length ~1 (so no concentration).
Outputs go to toy_toeplitz_results.txt.
"""
import numpy as np, sys

def phi(eps, y):
    return eps*np.exp(-eps*y)/(-np.expm1(-eps)) - 1.0

def V_of(pattern):
    q = len(pattern)
    r = np.arange(q)
    s = np.arange(q)
    return np.real(np.exp(2j*np.pi*np.outer(s, r)/q) @ np.asarray(pattern, float))

def closed_form(t, M, qs, pats, K=None):
    """F at points t (array) using levels 1..K."""
    K = len(qs) if K is None else K
    F = np.ones_like(t)
    P_prev = 1
    for k in range(K):
        q = qs[k]; V = V_of(pats[k])
        eps_prev = 2*np.pi*M/P_prev
        y = np.mod(P_prev*t/(2*np.pi), 1.0)
        acc = np.zeros_like(t)
        for sg in range(q):
            if V[sg] != 0:
                acc += V[sg]*phi(eps_prev, np.mod(y - sg/q, 1.0))
        F += (2.0/q)*acc
        P_prev *= q
    return F

def sequence(Nmax, qs, pats, K=None):
    K = len(qs) if K is None else K
    m = np.zeros(Nmax+1)
    P = [1]
    for q in qs[:K]:
        P.append(P[-1]*q)
    for n in range(1, Nmax+1):
        j = 0
        while j < K and n % P[j+1] == 0:
            j += 1
        if j < K:   # k(n) = j, level j+1
            r = (n // P[j]) % qs[j]
            m[n] = pats[j][r]
        else:
            m[n] = 0.0
    return m

def direct(t, M, m):
    n = np.arange(1, len(m))
    a = 2*M/(M+1j*n)
    # F = 1 + 2 Re sum m_n a_n e^{int}, chunked
    out = np.ones_like(t)
    for c0 in range(0, len(n), 20000):
        nn = n[c0:c0+20000]; aa = a[c0:c0+20000]*m[1+c0:1+c0+20000]
        out += 2*np.real(np.exp(1j*np.outer(t, nn)) @ aa)
    return out

def sym_pattern(q, rng, density=0.5):
    v = np.zeros(q)
    for r in range(1, q//2+1):
        b = 1.0 if rng.random() < density else 0.0
        v[r] = b; v[(q-r) % q] = b
    v[0] = 0.0
    return v

def main():
    rng = np.random.default_rng(1)
    out = []
    M = 0.5
    # ---- (1) closed form vs direct, K-level periodic truncation
    for name, qs, pats in [
        ("valuation q=2, c=(0,1,1,0,1,0,0,1)", [2]*8,
         [np.array([0, c], float) for c in [0, 1, 1, 0, 1, 0, 0, 1]]),
        ("random patterns q=(3,5,2,3,4)", [3, 5, 2, 3, 4],
         [sym_pattern(q, rng) for q in [3, 5, 2, 3, 4]]),
    ]:
        K = len(qs); P = int(np.prod(qs))
        m = sequence(400000, qs, pats, K)
        # sample points away from the atoms (grid 2pi j / P)
        tt = 2*np.pi*(np.arange(P) + rng.uniform(0.2, 0.8, P))/P
        tt = tt[:min(len(tt), 200)]
        Fc = closed_form(tt, M, qs, pats, K)
        errs = []
        for N in (100000, 200000, 400000):
            Fd = direct(tt, M, m[:N+1])
            errs.append(np.max(np.abs(Fc-Fd)))
        out.append(f"(1) {name}: P_K={P}, max|closed-direct| at N=1e5,2e5,4e5: "
                   + ", ".join(f"{e:.2e}" for e in errs)
                   + " (tail of the direct sum is O(1/(N dist)); should halve with N)")
    # ---- (2)+(3) local structure on level-i arcs
    def local_check(name, qs, pats, i, narcs=5, nu=4096):
        P = [1]
        for q in qs:
            P.append(P[-1]*q)
        eps = 2*np.pi*M/P[i]
        res = []
        for _ in range(narcs):
            a = rng.integers(0, P[i])
            u = (np.arange(nu)+0.5)/nu
            t = 2*np.pi*(a+u)/P[i]
            F = closed_form(t, M, qs, pats)
            # coarse part: levels 1..i give 1 + 2B e^{-eps u}; compute it exactly
            Fc = closed_form(t, M, qs[:i], pats[:i]) if i > 0 else np.ones_like(t)
            # predicted digit function X(u) = sum_l g_{i+l}(e_l)/Q_{i+l}
            X = np.zeros_like(u); Q = 1
            frac = u.copy()
            for l, k in enumerate(range(i, len(qs))):
                q = qs[k]; V = V_of(pats[k]); Q *= q
                e = np.floor(frac*q).astype(int); frac = frac*q - e
                T = np.array([V[d+1:].sum() for d in range(q)])
                g = np.arange(q) - 2*T
                X += g[e]/Q
            # model: F = Fc + eps*(X - u*(1+2B... )) ; we test the fine part only:
            fine = F - Fc
            # first-order prediction of the fine part: eps*(const - 2 sum T/Q)
            Xf = np.zeros_like(u); Q = 1; frac = u.copy()
            for k in range(i, len(qs)):
                q = qs[k]; V = V_of(pats[k]); Q *= q
                e = np.floor(frac*q).astype(int); frac = frac*q - e
                T = np.array([V[d+1:].sum() for d in range(q)])
                Xf += -2*T[e]/Q
            pred = eps*(Xf - Xf.mean()) + fine.mean()
            res.append(np.max(np.abs(fine - pred))/eps)
        return res
    for name, qs, pats, i in [
        ("valuation q=2 (16 levels), c random", [2]*16,
         [np.array([0, c], float) for c in rng.integers(0, 2, 16)], 6),
        ("valuation q=3 (10 levels), c random", [3]*10,
         [np.array([0, c, c], float) for c in rng.integers(0, 2, 10)], 4),
        ("random patterns q in {2,3,5} (9 levels)", [3, 5, 2, 3, 5, 2, 3, 2, 5],
         None, 3),
    ]:
        if pats is None:
            pats = [sym_pattern(q, rng) for q in qs]
        r = local_check(name, qs, pats, i)
        P_i = int(np.prod(qs[:i]))
        out.append(f"(2) {name}, level i={i} (eps_i={2*np.pi*M/P_i:.3e}): "
                   f"max |fine - first-order| / eps_i over 5 arcs = {max(r):.3e} (should be O(eps_i))")
    # ---- (3) law of the u-dependent part for the valuation family.
    # On a level-i arc, F = 1-2c_i + 2B e^{-eps u} + fine; first order:
    # (F - const)/eps = sum_l s_l e_l 2^{-l}, s_l = -2B + 2c_i - 2c_{i+l}.
    # On arcs where F vanishes somewhere, B ~ c_i - 1/2 and s_l = 1 - 2c_{i+l} = +-1: uniform law.
    qs = [2]*22; cs = rng.integers(0, 2, 22)
    pats = [np.array([0, c], float) for c in cs]
    i = 9; P_i = 2**i; eps = 2*np.pi*M/P_i
    c = np.concatenate([[0], np.cumsum(np.diff(np.concatenate([[0], cs])))])  # c_j, j=0..
    cvals = np.concatenate([[0], cs])   # c_0 = 0, c_j = cs[j-1]
    u = (np.arange(1 << 13)+0.5)/(1 << 13)
    worst_gen, worst_zero, nzero = 0.0, 0.0, 0
    for a in rng.integers(0, P_i, 400):
        t = 2*np.pi*(a+u)/P_i
        F = closed_form(t, M, qs, pats)
        Fc = closed_form(t, M, qs[:i], pats[:i])
        ci = cvals[i]
        B = 0.5*(Fc - 1 + 2*ci)*np.exp(eps*u)        # constant in u (check)
        Bm = B.mean()
        # predicted first-order u-dependence
        X = np.zeros_like(u); frac = u.copy()
        for l in range(1, len(qs) - i + 1):
            e = np.floor(frac*2).astype(int); frac = frac*2 - e
            s = -2*Bm + 2*ci - 2*cvals[i+l]
            X += s*e*2.0**(-l)
        Z = (F - F.mean())/eps
        dev = np.max(np.abs((Z - Z.mean()) - (X - X.mean())))
        worst_gen = max(worst_gen, dev)
        if F.min() < 0 < F.max():
            nzero += 1
            Zs = np.sort(Z)
            # compare with uniform on [Zs[0], Zs[-1]]: Kolmogorov distance
            Uref = Zs[0] + (Zs[-1]-Zs[0])*(np.arange(len(Zs))+0.5)/len(Zs)
            kd = np.max(np.abs(Zs - Uref))/(Zs[-1]-Zs[0])
            worst_zero = max(worst_zero, kd)
            if nzero == 1:
                out.append(f"(3) first arc with a sign change: B = {Bm:.4f}, c_i - 1/2 = {ci-0.5:.4f}, "
                           f"range of Z = {Zs[-1]-Zs[0]:.4f} (uniform law predicts 1)")
    out.append(f"(3) valuation q=2, 22 levels, 400 arcs of level 9 (eps={eps:.2e}): max deviation of "
               f"(F-mean)/eps from the first-order digit law = {worst_gen:.2e}; arcs with a sign change: "
               f"{nzero}, max normalised quantile distance to uniform there = {worst_zero:.3f}")
    txt = "\n".join(out)
    print(txt)
    with open("toy_toeplitz_results.txt", "w") as f:
        f.write("toy_toeplitz.py results (AI-produced, 2026-11-26)\n" + txt + "\n")

if __name__ == "__main__":
    main()
