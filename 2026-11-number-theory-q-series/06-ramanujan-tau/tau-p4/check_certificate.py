#!/usr/bin/env python3
"""
Checker for certificate_tau_p4.txt (written by tau_p4.gp).  Python, exact integers, no shared code
with the GP generator or with sieve_tau_p4.py.

For each block (alpha, e, M) it checks, from scratch:
  1. the seed list is exactly the set of integer solutions (U, V) of U^2 - 3UV + V^2 = alpha with
     U^2 <= alpha/2 (alpha > 0) or U^2 <= |alpha| (alpha < 0)   [brute force over U and V];
  2. for every prime q used as a witness: q is prime and A^M = I mod q for A = [[0,1],[-1,3]]
     (so n -> (Z_n, Z_{n+1}) mod q has period dividing M for every seed, in both directions);
  3. every residue class n mod M of every seed has a line; on a witness line (n, q, code, z, w):
     z = Z_n mod q and w = Z_{n+1} mod q (recomputed), and
       code P: z is not in {x^e mod q : 0 <= x < q}     (brute-force enumeration)
       code S: w is not in {x^2 mod q : 0 <= x < q};
  4. SURVIVES lines are reported (expected: none, except n = 0 of seed (1, 9) for alpha = 55).
"""
import sys, re
from math import isqrt

def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))

def matpow_is_id(q, M):
    def mul(X, Y):
        return ((X[0]*Y[0] + X[1]*Y[2]) % q, (X[0]*Y[1] + X[1]*Y[3]) % q,
                (X[2]*Y[0] + X[3]*Y[2]) % q, (X[2]*Y[1] + X[3]*Y[3]) % q)
    R, B, k = (1, 0, 0, 1), (0, 1, q - 1, 3 % q), M
    while k:
        if k & 1: R = mul(R, B)
        B = mul(B, B); k >>= 1
    return R == (1, 0, 0, 1)

def all_seeds(alpha):
    B0 = isqrt(alpha // 2) if alpha > 0 else isqrt(-alpha)
    out = set()
    for U in range(-B0, B0 + 1):
        # V^2 - 3UV + U^2 - alpha = 0  =>  |V| <= 3|U| + sqrt(|alpha|) + 1
        lim = 3 * abs(U) + isqrt(abs(alpha)) + 2
        for V in range(-lim, lim + 1):
            if U * U - 3 * U * V + V * V == alpha:
                out.add((U, V))
    return B0, out

def main(path):
    lines = open(path).read().split("\n")
    i, ok, nblocks, nlines, survivors = 0, True, 0, 0, []
    while i < len(lines):
        m = re.match(r"# alpha = (-?\d+)  exponent e = (\d+)  modulus M = (\d+)  seed bound \|U\| <= (\d+)", lines[i])
        if not m:
            i += 1; continue
        alpha, e, M, B0 = map(int, m.groups()); nblocks += 1
        seeds = [tuple(map(int, s)) for s in re.findall(r"\[(-?\d+), (-?\d+)\]", lines[i + 1].split(":", 1)[1])]
        b0, true_seeds = all_seeds(alpha)
        if b0 != B0 or set(seeds) != true_seeds or len(seeds) != len(true_seeds):
            print("SEED MISMATCH", alpha, seeds, sorted(true_seeds)); ok = False
        i += 3
        pw_cache, sq_cache, per_ok = {}, {}, {}
        for s in seeds:
            hm = re.match(r"seed \[(-?\d+), (-?\d+)\]:", lines[i])
            if not hm or (int(hm.group(1)), int(hm.group(2))) != s:
                print("BAD SEED HEADER", alpha, s, lines[i][:40]); ok = False
            i += 1
            Z = [s[0], s[1]]
            for n in range(2, M + 2): Z.append(3 * Z[-1] - Z[-2])
            for n in range(M):
                parts = lines[i].split(); i += 1; nlines += 1
                if int(parts[0]) != n:
                    print("MISSING CLASS", alpha, s, n); ok = False; continue
                if parts[1] == "SURVIVES":
                    survivors.append((alpha, e, s, n)); continue
                q, code, z, w = int(parts[1]), parts[2], int(parts[3]), int(parts[4])
                if q not in per_ok:
                    per_ok[q] = is_prime(q) and matpow_is_id(q, M)
                    pw_cache[q] = {pow(x, e, q) for x in range(q)}
                    sq_cache[q] = {x * x % q for x in range(q)}
                good = per_ok[q] and z == Z[n] % q and w == Z[n + 1] % q
                if code == "P": good = good and z not in pw_cache[q]
                elif code == "S": good = good and w not in sq_cache[q]
                else: good = False
                if not good:
                    print("BAD LINE", alpha, s, lines[i - 1]); ok = False
        print(f"alpha={alpha:5d} e={e:2d} M={M:3d}: {len(seeds)} seeds, witness primes {sorted(per_ok)}")
    print(f"blocks: {nblocks}, class lines checked: {nlines}, survivors: {survivors}")
    print("CERTIFICATE OK" if ok else "CERTIFICATE FAILED")
    return ok

if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1] if len(sys.argv) > 1 else "certificate_tau_p4.txt") else 1)
