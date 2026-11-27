#!/usr/bin/env python3
"""Referee 3: brute-force enumeration of the (ii),(iii),(iv) multisums (all tuples r_i, s_i >= 0, the
printed convention 1/(q;q)_m = 1/(q^2;q^2)_m = 0 for m < 0, Omega literally), to cross-check series_q.py."""
import sys
import itertools
sys.argv = [sys.argv[0], sys.argv[1] if len(sys.argv) > 1 else "34", "4"]
import series_q as S
from flint import fmpz_poly

N = S.N


def brute(k, which):
    tot = fmpz_poly([])
    # bound: (r_i+s_i)^2 + s_i^2 < N for every i
    rng_r = range(0, int(N ** 0.5) + 1)
    rng_s = range(0, int((N / 2) ** 0.5) + 1)
    for rs in itertools.product(rng_r, repeat=k):
        for ss in itertools.product(rng_s, repeat=k):
            r = list(rs) + [0]
            s = [0] + list(ss)          # s[0] = s_0 = 0, s[i] = s_i
            if any(r[i] - r[i + 1] < 0 for i in range(k)):
                continue                # 1/(q;q)_m = 0, m < 0
            if any(s[i] - s[i - 1] < 0 for i in range(1, k + 1)):
                continue
            e = 0
            for i in range(k):
                ri, si = r[i], s[i + 1]
                e += (ri + si) ** 2 + si ** 2
                if which == "ii":
                    e += ri + 2 * si
            if e >= N:
                continue
            term = S.mono(e)
            for i in range(k):
                term = S.mul(term, S.A[r[i] - r[i + 1]])
                term = S.mul(term, S.B[s[i + 1] - s[i]])
            if which == "iv":
                # Omega = sum_{i=1}^{k-1} q^{r_i+2s_i}(1 - q^{2s_{i+1}-2s_i}) + q^{r_k+2s_k}
                om = fmpz_poly([])
                for i in range(1, k):
                    om += S.mono(r[i - 1] + 2 * s[i]) - S.mono(r[i - 1] + 2 * s[i + 1])
                om += S.mono(r[k - 1] + 2 * s[k])
                term = S.mul(term, om)
            tot += term
    return tot


for k in (1, 2, 3):
    for which, ex in (("ii", S.ex_ii), ("iii", S.ex_iii)):
        b = brute(k, which)
        d = S.multisum(k, ex)
        print(f"k={k} ({which}): brute == nested: {(b - d).is_zero()}")
    b = brute(k, "iv")
    d = S.multisum_omega(k, S.ex_iii)
    print(f"k={k} (iv):  brute == nested: {(b - d).is_zero()}")
