#!/usr/bin/env python3
"""
Sieve proof that tau(p^4) = alpha has no solution, for alpha in {461, -599}.

tau(p^4) = t^4 - 3 p^11 t^2 + p^22 with t = tau(p).  So (U, V) = (p^11, t^2) is an integer solution of
        U^2 - 3UV + V^2 = alpha.                                                    (*)
T(U,V) = (V, 3V - U) preserves (*).  Every solution lies in the T-orbit of a solution with |U| <= B0
(B0 = sqrt(alpha/2) if alpha > 0, sqrt(|alpha|) if alpha < 0; see the proof file).  For each such seed
(Z_0, Z_1) = (U, V) the orbit is {(Z_n, Z_{n+1}) : n in Z}, Z_{n+1} = 3 Z_n - Z_{n-1}.
We show: for no n is Z_n an 11th power and Z_{n+1} a square, by exhibiting, for every residue class
n mod M, a prime q (with the period of Z mod q dividing M) such that Z_n mod q is not an 11th power
residue or Z_{n+1} mod q is not a quadratic residue (0 counts as a residue in both cases).
"""
import sys, json
from math import isqrt, gcd
import numpy as np

def seeds(alpha):
    B0 = isqrt(abs(alpha) // 2) if alpha > 0 else isqrt(abs(alpha))
    out = []
    for U in range(-B0, B0 + 1):
        D = 5 * U * U + 4 * alpha
        if D < 0: continue
        s = isqrt(D)
        if s * s != D: continue
        for sg in {s, -s}:
            if (3 * U + sg) % 2 == 0:
                V = (3 * U + sg) // 2
                assert U * U - 3 * U * V + V * V == alpha
                out.append((U, V))
    return out, B0

def period(q):
    # order of the matrix [[0,1],[-1,3]] mod q
    a, b, c, d = 0, 1, q - 1, 3
    A = (a, b, c, d); P = A; k = 1
    while P != (1, 0, 0, 1):
        P = ((P[0]*A[0] + P[1]*A[2]) % q, (P[0]*A[1] + P[1]*A[3]) % q,
             (P[2]*A[0] + P[3]*A[2]) % q, (P[2]*A[1] + P[3]*A[3]) % q)
        k += 1
    return k

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\0\0"
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def run(alpha, M, qmax, verbose=True):
    sd, B0 = seeds(alpha)
    qs = []
    for q in primes_upto(qmax):
        if q in (2, 5): continue
        # period divides q-1 or q+1 (or their halves); cheap pre-test before computing it
        if all(M % ((q + e)) and M % ((q + e) // 2) for e in (-1, 1)) and gcd(M, q * q - 1) < 3:
            continue
        T = period(q)
        if M % T == 0:
            qs.append((q, T))
    cert = {"alpha": alpha, "M": M, "B0": B0, "seeds": sd, "primes": qs, "kills": {}}
    ok_all = True
    for (U, V) in sd:
        alive = np.ones(M, dtype=bool)
        used = []
        for q, T in qs:
            if not alive.any(): break
            p11 = {pow(x, 11, q) for x in range(q)} if q % 11 == 1 else None   # None: every residue is an 11th power
            sq = {x * x % q for x in range(q)}
            z = [U % q, V % q]
            for n in range(2, T + 1): z.append((3 * z[-1] - z[-2]) % q)
            assert z[T] == z[0] and (3 * z[T] - z[T - 1]) % q == z[1]
            good = np.array([(p11 is None or z[n] in p11) and (z[n + 1] in sq) for n in range(T)], dtype=bool)
            before = int(alive.sum())
            alive &= np.tile(good, M // T)
            if int(alive.sum()) < before: used.append(q)
        left = np.nonzero(alive)[0].tolist()
        cert["kills"][f"{U},{V}"] = {"primes_used": used, "survivors": left[:20], "n_survivors": len(left)}
        if verbose:
            print(f"alpha={alpha} seed (Z0,Z1)=({U},{V}): survivors mod {M}: {len(left)} {left[:10]}  primes used: {used}")
        ok_all &= (len(left) == 0)
    cert["proved"] = ok_all
    return cert

if __name__ == "__main__":
    alpha = int(sys.argv[1]); M = int(sys.argv[2]); qmax = int(sys.argv[3])
    c = run(alpha, M, qmax)
    print("PROVED" if c["proved"] else "NOT PROVED")
    if len(sys.argv) > 4:
        json.dump(c, open(sys.argv[4], "w"))
