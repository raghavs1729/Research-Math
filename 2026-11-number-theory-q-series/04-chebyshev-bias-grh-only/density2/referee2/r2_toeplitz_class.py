"""Referee 2: the class of sequences in TOY.md section 2.2 ("Setting").

Claims of TOY.md:
  (C1) "These are the regular Toeplitz sequences whose only hole is n = 0
        (for every k, m is P_k-periodic off P_k Z)."
  (C2) "every periodic sequence (period Q: take P_k = Q 2^{k-1}, with constant patterns
        from k = 2 on)".
Tests:
  (a) m_n = [nu_2(n) + nu_3(n) odd].  Every n != 0 is periodic (period 2^{nu_2+1} 3^{nu_3+1},
      also for the two-sided sequence m_|n|); the 6^k-skeleton has density
      (1-2^-k)(1-3^-k) -> 1 (regular); 0 is the only non-periodic position.  But m is NOT
      P-periodic off PZ for any P in [2, PMAX]  (so it is not in the class of Thm 2.3 for any
      period structure: the identification (C1) is false as stated).
  (b) a one-hole-per-period Toeplitz sequence whose hole is NOT at 0: P_1 = 2, q_k = 3
      (k >= 2), holes h_k = P_k/2 = 3^{k-1}: m_n = a for n even, m_n = c_j for n odd with
      nu_3(n) = j.  Every n >= 1 is filled (no hole in Z at all), but m is not P-periodic off
      PZ for any P in [2, PMAX]: "single-hole Toeplitz" (TOY.md intro) is wider than the class.
  (c) (C2): random periodic sequences are reproduced by the construction of TOY.md.
"""
import numpy as np

NMAX = 400000
PMAX = 3000
rng = np.random.default_rng(5)


def nu(n, p):
    n = np.asarray(n).copy()
    v = np.zeros_like(n)
    mask = n % p == 0
    while mask.any():
        v[mask] += 1
        n[mask] //= p
        mask = (n % p == 0) & (n != 0)
    return v


def periodic_off_PZ(m, P):
    """Is m (indexed 1..len-1) P-periodic off PZ on the available range?"""
    n = np.arange(1, len(m) - P)
    sel = n % P != 0
    return bool(np.all(m[n[sel]] == m[n[sel] + P]))


def toeplitz_seq(qs, vs, N):
    m = np.zeros(N + 1, dtype=np.int64)
    n = np.arange(1, N + 1)
    P = 1
    for q, v in zip(qs, vs):
        sel = (n % P == 0) & ((n // P) % q != 0)
        m[1:][sel] = np.asarray(v)[(n[sel] // P) % q]
        P *= q
        if P > N:
            break
    return m


def main():
    out = ["Referee 2, r2_toeplitz_class.py: the class of TOY.md section 2.2."]
    n = np.arange(NMAX + 1)
    n[0] = 1
    v2, v3 = nu(n, 2), nu(n, 3)
    m = ((v2 + v3) % 2).astype(np.int64)
    m[0] = 0
    # (a) periodicity of every position n <= 3000, two-sided
    bad = 0
    for k in range(1, 3001):
        p = 2 ** (int(v2[k]) + 1) * 3 ** (int(v3[k]) + 1)
        js = np.arange(k % p, NMAX + 1, p)
        js = js[js > 0]
        if not np.all(m[js] == m[k]):
            bad += 1
        # negative side: m_|k - j p| for j with k - j p < 0 equals m at |k - jp| = jp - k
        jn = np.arange(1, NMAX // p)
        pos = jn * p - k
        pos = pos[(pos > 0) & (pos <= NMAX)]
        if not np.all(m[pos] == m[k]):
            bad += 1
    out.append(f"(a) m_n = [nu2(n)+nu3(n) odd]: positions 1..3000 with a failure of their period (both sides): {bad}")
    for k in [1, 2, 4, 8]:
        dens = np.mean((v2[1:] < k) & (v3[1:] < k))
        out.append(f"    density of the 6^{k}-skeleton (nu2 < {k}, nu3 < {k}): {dens:.6f}  (predicted {(1-2.0**-k)*(1-3.0**-k):.6f})")
    good = [P for P in range(2, PMAX + 1) if periodic_off_PZ(m, P)]
    out.append(f"    P in [2,{PMAX}] with m P-periodic off PZ: {good if good else 'none'}  -> not in the class of Thm 2.3")
    # (b) hole at P_k/2
    a, c = 1, rng.integers(0, 2, 40)
    c[0] = 1 - a  # avoid the degenerate case c_0 = a (then m is 3-periodic off 3Z)
    mb = np.zeros(NMAX + 1, dtype=np.int64)
    even = (n % 2 == 0)
    mb[even] = a
    mb[~even] = c[v3[~even]]
    mb[0] = 0
    goodb = [P for P in range(2, PMAX + 1) if periodic_off_PZ(mb, P)]
    out.append(f"(b) hole class h_k = P_k/2 (P_k = 2*3^(k-1)), c_j = {list(map(int, c[:8]))}...: every n >= 1 filled;")
    out.append(f"    P in [2,{PMAX}] with m P-periodic off PZ: {goodb if goodb else 'none'}  -> not in the class of Thm 2.3")
    out.append("    (the class needs an infinite chain P_k | P_{k+1}, P_k -> oo, of such P)")
    # check it is P_k-periodic off h_k + P_k Z
    okb = True
    for k in range(1, 9):
        Pk = 2 * 3 ** (k - 1)
        hk = 3 ** (k - 1) % Pk if k > 1 else 1
        nn = np.arange(1, NMAX - Pk)
        sel = (nn % Pk) != hk
        okb &= bool(np.all(mb[nn[sel]] == mb[nn[sel] + Pk]))
    out.append(f"    ... but it IS P_k-periodic off h_k + P_k Z for k = 1..8: {okb}")
    # (c) periodic sequences
    fails = 0
    for _ in range(300):
        Q = int(rng.integers(2, 31))
        base = rng.integers(0, 2, Q)  # base[r] = m_r for r mod Q, base[0] = m_Q
        mm = base[np.arange(NMAX + 1) % Q]
        qs = [Q] + [2] * 40
        vs = [np.concatenate([[0], base[1:]])] + [np.array([0, base[0]])] * 40
        mt = toeplitz_seq(qs, vs, 100000)
        if not np.all(mt[1:] == mm[1:100001]):
            fails += 1
    out.append(f"(c) 300 random periodic sequences (Q in [2,30]) rebuilt by TOY's construction: {fails} failures")
    txt = "\n".join(out)
    print(txt)
    with open(__file__.replace(".py", "_results.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
