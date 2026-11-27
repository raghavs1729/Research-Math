#!/usr/bin/env python3
"""Shared helpers: selection and identification of the real primitive character."""
from __future__ import annotations

import hashlib
import math
import os
import re
from fractions import Fraction

from flint import arb, acb, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))


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


def read_zero_windows(path: str) -> list[arb]:
    """LOW-LEVEL reader: each listed decimal d stands for the window
    [d-1e-20, d+1e-20].  Checks only positivity and strict ordering; it says
    NOTHING about completeness.  Certified programs must use
    load_certified_zeros instead."""
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


# ---------------------------------------------------------------------------
# binding of a zero file to its completeness certificate (audit finding G1)
# ---------------------------------------------------------------------------

# Hypothesis (H-Z) of NONMONOTONE.md: the files on which the theorems rest.
# These are the hashes pinned by the independent audit
# (monotonicity/nonmonotone-audit/AUDIT.md, section 5.2).
HZ_SHA256 = {
    163: "b9d38aba6df821dea401035ea6ae24aa045305e757f81ea1a6233b632ad25b9d",
    232: "42f40b29d79cc5228478acfd0480f23939c6f3b5d0479acaf178c661fa0af378",
}

RECORD_RE = re.compile(
    r"^HASH-RECORD sha256=([0-9a-f]{64}) file=(\S+) q=(\d+) K=(\d+) "
    r"T_cert=(\d+)/(\d+) W=(\d+) mode=(\S+) status=PASS\s*$")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def cert_log_path(q: int) -> str:
    return os.path.join(HERE, "zeros_q%d.cert.txt" % q)


def hash_record_line(sha: str, path: str, q: int, K: int, tnum: int, tden: int,
                     W: int, mode: str) -> str:
    return ("HASH-RECORD sha256=%s file=%s q=%d K=%d T_cert=%d/%d W=%d mode=%s status=PASS"
            % (sha, os.path.basename(path), q, K, tnum, tden, W, mode))


def find_hash_record(q: int, sha: str):
    """The LAST passing record for this sha256 in the certificate log of q
    (written by zeros_cert.py only after the windows and the count W = 2K have
    been certified for exactly the file with that hash), or None."""
    path = cert_log_path(q)
    if not os.path.exists(path):
        return None
    found = None
    with open(path) as fh:
        for line in fh:
            m = RECORD_RE.match(line)
            if m and m.group(1) == sha and int(m.group(3)) == q:
                found = dict(sha256=sha, file=m.group(2), K=int(m.group(4)),
                             T_cert=Fraction(int(m.group(5)), int(m.group(6))),
                             W=int(m.group(7)), mode=m.group(8))
    return found


def load_certified_zeros(q: int, path: str, allow_uncertified: bool = False):
    """Zero windows of the file, REFUSING unless the sha256 of the file is
    recorded in a passing HASH-RECORD of zeros_q<q>.cert.txt and the record is
    consistent with the file (number of windows K, all windows below T_cert,
    W = 2K).  Returns (windows, meta).  With allow_uncertified=True a file
    without record is accepted and meta['certified'] is False: the caller must
    then not label anything as certified."""
    sha = sha256_file(path)
    vals = read_zero_windows(path)
    rec = find_hash_record(q, sha)
    if rec is None:
        if not allow_uncertified:
            raise ArithmeticError(
                "REFUSED: zero file %s (sha256 %s) has no passing completeness record in %s. "
                "Run: python3 zeros_cert.py %d --verify [--file ... --T p/q]"
                % (path, sha, cert_log_path(q), q))
        return vals, dict(sha256=sha, certified=False, K=len(vals), T_cert=None)
    T = arb(rec["T_cert"].numerator) / arb(rec["T_cert"].denominator)
    if rec["K"] != len(vals) or rec["W"] != 2 * len(vals) or not (vals[-1].upper() < T):
        raise ArithmeticError("REFUSED: certificate record inconsistent with the zero file")
    rec["certified"] = True
    return vals, rec
