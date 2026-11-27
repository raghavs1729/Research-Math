"""Audit 2: independent brute-force test of Theorem C (no persistence of p in GCD_n, GCD_{n+1}).

For a_n = A n + k, b_n = B n^2 + C n + E (integers) and a prime p we compute V_n = (p_n, q_n) mod p
for n <= NMAX (vectorised over many (fraction, prime) pairs) and record whether some n has
V_n == V_{n+1} == 0 mod p (persistence; this state is absorbing).

Hypotheses of Theorem C: p odd, p does not divide B*omega (omega = 2kB - A(B+C)),
D = A^2 + 4B not a nonzero square mod p.
Reports: violations (hyp and persistence), and persistence counts when exactly one hypothesis is dropped.
"""
import numpy as np, random, sys
from collections import Counter

NMAX = 5000
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 12345)
PRIMES = [p for p in range(2, 400) if all(p % q for q in range(2, int(p**0.5) + 1))]
SMALL = [2, 3, 5, 7, 11, 13]


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def gen(kind):
    r = lambda m: random.randint(-m, m)
    nz = lambda m: random.choice([i for i in range(-m, m + 1) if i])
    if kind == 'generic':
        return nz(9), r(12), nz(9), r(15), r(15)
    if kind == 'irreducible':
        while True:
            B, C, E = nz(9), r(15), r(15)
            d = C * C - 4 * B * E
            if d < 0 or int(round(d ** 0.5)) ** 2 != d:
                return nz(9), r(12), B, C, E
    if kind == 'double_root':  # b = s (u n - v)^2
        s, u, v = nz(4), random.randint(1, 3), r(6)
        return nz(9), r(12), s * u * u, -2 * s * u * v, s * v * v
    if kind == 'D_zero':  # A = 2K, B = -K^2
        K = nz(4)
        return 2 * K, r(12), -K * K, r(15), r(15)
    if kind == 'b_vanishes':  # b = (n - N)(c1 n + c0), N positive integer
        N0, c1, c0 = random.randint(1, 30), nz(5), r(10)
        return nz(9), r(12), c1, c0 - c1 * N0, -c0 * N0
    if kind == 'omega_zero_Dnonsq':  # k = A(B+C)/(2B) integral
        while True:
            A, B, C, E = nz(9), nz(9), r(15), r(15)
            if (A * (B + C)) % (2 * B) == 0:
                return A, A * (B + C) // (2 * B), B, C, E
    raise ValueError


def run(pairs):
    P = np.array([q[1] for q in pairs], dtype=np.int64)
    A = np.array([q[0][0] for q in pairs], dtype=np.int64) % P
    K = np.array([q[0][1] for q in pairs], dtype=np.int64) % P
    B = np.array([q[0][2] for q in pairs], dtype=np.int64) % P
    C = np.array([q[0][3] for q in pairs], dtype=np.int64) % P
    E = np.array([q[0][4] for q in pairs], dtype=np.int64) % P
    pm, qm = np.ones_like(P), np.zeros_like(P)          # n = -1
    p0, q0 = K % P, np.ones_like(P)                     # n = 0
    pers = np.full(len(P), -1, dtype=np.int64)
    for n in range(1, NMAX + 2):
        nn = n % P
        an = (A * nn + K) % P
        bn = ((B * nn % P) * nn + C * nn + E) % P
        p1 = (an * p0 + bn * pm) % P
        q1 = (an * q0 + bn * qm) % P
        hit = (p0 == 0) & (q0 == 0) & (p1 == 0) & (q1 == 0) & (pers < 0)
        pers[hit] = n - 1
        pm, qm, p0, q0 = p0, q0, p1, q1
    return pers


def main():
    kinds = ['generic', 'irreducible', 'double_root', 'D_zero', 'b_vanishes', 'omega_zero_Dnonsq']
    pairs, info = [], []
    for kind in kinds:
        for _ in range(1500):
            f = gen(kind)
            for p in SMALL + random.sample(PRIMES[6:], 6):
                A, k, B, C, E = f
                om = 2 * k * B - A * (B + C)
                D = A * A + 4 * B
                if p == 2:
                    Dsq = (D % 2 == 1)
                else:
                    Dsq = legendre(D, p) == 1
                flags = dict(odd=p != 2, pB=B % p != 0, pom=om % p != 0, Dns=not Dsq)
                pairs.append((f, p)); info.append((kind, flags))
    pers = run(pairs)
    viol = []
    cnt = Counter()
    for (f, p), (kind, fl), n in zip(pairs, info, pers):
        hyp = all(fl.values())
        key = 'HYP' if hyp else 'drop:' + ','.join(sorted(k for k, v in fl.items() if not v))
        cnt[(kind, key, n >= 0)] += 1
        if hyp and n >= 0:
            viol.append((f, p, int(n)))
    print("pairs:", len(pairs), " NMAX:", NMAX)
    print("violations of Theorem C:", len(viol), viol[:10])
    agg = Counter()
    for (kind, key, ps), c in cnt.items():
        agg[(key, ps)] += c
    for key in sorted({k for k, _ in agg}):
        print(f"  {key:28s} persistent {agg[(key, True)]:6d}  not {agg[(key, False)]:6d}")
    print("per kind (HYP pairs tested):", {kd: cnt[(kd, 'HYP', False)] + cnt[(kd, 'HYP', True)] for kd in kinds})
    # examples of persistence with one hypothesis dropped
    ex = {}
    for (f, p), (kind, fl), n in zip(pairs, info, pers):
        bad = [k for k, v in fl.items() if not v]
        if n >= 0 and len(bad) == 1 and bad[0] not in ex:
            ex[bad[0]] = (f, p, int(n))
    print("examples (A,k,B,C,E), p, first n:", ex)


if __name__ == '__main__':
    main()
