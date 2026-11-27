#!/usr/bin/env python3
"""Shared helpers: selection and identification of the real primitive character."""
from __future__ import annotations

import math

from flint import arb, acb, dirichlet_char


def kronecker(a: int, n: int) -> int:
    """Kronecker symbol (a/n) for n >= 1 (exact integer arithmetic)."""
    assert n >= 1
    if n == 1:
        return 1
    res = 1
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    if v:
        if a % 2 == 0:
            return 0
        if v % 2 == 1 and a % 8 in (3, 5):
            res = -res
    # Jacobi symbol (a/n), n odd positive
    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def real_character(q: int, parity: int = 1):
    """The real primitive Dirichlet character of conductor q and given parity
    (parity 1 = odd), as a python-flint dirichlet_char, together with the
    fundamental discriminant d such that chi(n) = Kronecker(d/n); the identity
    chi(n) = (d/n) is CHECKED for every 1 <= n <= 2q (it is periodic mod q)."""
    found = []
    for l in range(1, q):
        if math.gcd(l, q) != 1:
            continue
        chi = dirichlet_char(q, l)
        if (chi.is_real() and chi.is_primitive() and not chi.is_principal()
                and chi.conductor() == q and chi.parity() == parity):
            found.append(chi)
    if len(found) != 1:
        raise ArithmeticError("expected exactly one real primitive character "
                              "mod %d of parity %d, found %d" % (q, parity, len(found)))
    chi = found[0]
    d = -q if parity == 1 else q
    vals = character_values(chi, q)
    for n in range(1, 2 * q + 1):
        if vals[n % q] != kronecker(d, n):
            raise ArithmeticError("character is not the Kronecker symbol (%d/.) at n=%d" % (d, n))
    if vals[q - 1] != (-1 if parity == 1 else 1):
        raise ArithmeticError("parity check failed")
    return chi, d


def character_values(chi, q: int) -> list[int]:
    """Exact values chi(0..q-1) in {-1,0,1} of a real character."""
    vals = [0] * q
    exps = {}
    for n in range(1, q):
        if math.gcd(n, q) == 1:
            exps[n] = int(chi.chi_exponent(n))
    nz = sorted(set(exps.values()))
    # a real nonprincipal character takes exactly two values: exponent 0 <-> +1
    if len(nz) != 2 or nz[0] != 0:
        raise ArithmeticError("unexpected exponent set %r" % nz)
    for n, e in exps.items():
        vals[n] = 1 if e == 0 else -1
    # multiplicativity spot check (exact) to make sure 'nonzero exponent' means -1
    for n in range(1, q):
        for k in (2, 3, 5, 7):
            if vals[n] and vals[k % q]:
                if vals[(n * k) % q] != vals[n] * vals[k % q]:
                    raise ArithmeticError("multiplicativity check failed")
    if sum(vals) != 0:
        raise ArithmeticError("character sum over a period is not zero")
    return vals


def max_partial_sum(vals: list[int]) -> int:
    """max_x |sum_{n<=x} chi(n)| (exact; the sum over a full period vanishes)."""
    s = 0
    best = 0
    for v in vals:
        s += v
        best = max(best, abs(s))
    assert s == 0
    return best


def load_zeros(path: str) -> list[arb]:
    """Each listed decimal d stands for the certified window [d-1e-20, d+1e-20]."""
    with open(path, encoding="ascii") as fh:
        vals = [arb(line.split()[0], "1e-20") for line in fh
                if line.strip() and not line.startswith("#")]
    if not vals:
        raise ArithmeticError("empty zero file: " + path)
    for a, b in zip(vals, vals[1:]):
        if not a.upper() < b.lower():
            raise ArithmeticError("zero intervals are not strictly ordered")
    if not vals[0].lower() > 0:
        raise ArithmeticError("first zero window not positive")
    return vals
