#!/usr/bin/env python3
"""
m0.py -- certified enclosure of M_0(chi), the root of

    Q(u) = sum_{gamma>0} 4 u^3 / (u^2 + gamma^2)^(3/2) = 1

(Theorem B of monotonicity/CORE-NOTES.md: under GRH + L(1/2,chi) != 0 + LI,
delta_chi is strictly decreasing for 0 < M <= M_0).

Usage:  nice -n 10 python3 m0.py q [--out file]

Q is strictly increasing in u (each term is).  With the certified complete
list of zeros below gamma_last and the exact tail variance
V(u) = sigma^2(u) - sum_listed c^2/2 = sum_unlisted 8u^2/(u^2+gamma^2)  (GRH),

    Q_list(u) <= Q(u) <= Q_list(u) + V(u) * u / (2 sqrt(u^2 + gamma_last^2)),

because 4u^3/(u^2+g^2)^(3/2) = [8u^2/(u^2+g^2)] * [u / (2 sqrt(u^2+g^2))] and
g > gamma_last for unlisted zeros.  So M_0 lies in (u_lo, u_hi) as soon as
the upper bound of Q(u_lo) is < 1 and Q_list(u_hi) > 1 (both Arb-certified).
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from fractions import Fraction

from flint import arb, ctx

from common import load_certified_zeros
from density import Character, zero_file, frac_to_arb, rvm_guard, CTX_BITS


def q_bounds(char: Character, zeros, u: arb):
    Ql = sum((4 * u ** 3 / ((u * u + g * g) ** 3).sqrt() for g in zeros), arb(0))
    V = char.sigma2(u) - sum((8 * u * u / (u * u + g * g) for g in zeros), arb(0))
    if not V > 0:
        raise ArithmeticError("tail variance not certified positive")
    rvm_guard(char.q, u, zeros, V)
    tail = V * u / (2 * (u * u + zeros[-1].lower() ** 2).sqrt())
    return Ql, Ql + tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    ctx.prec = CTX_BITS
    char = Character(args.q)
    zeros, meta = load_certified_zeros(args.q, zero_file(args.q))   # refuses uncertified files
    gf = [float(g.mid()) for g in zeros]
    Qf = lambda u: sum(4 * u ** 3 / (u * u + g * g) ** 1.5 for g in gf)
    lo, hi = 1e-3, 50.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if Qf(mid) < 1:
            lo = mid
        else:
            hi = mid
    root = 0.5 * (lo + hi)
    # certify a bracket on a decimal grid, widening until both sides certify
    for digits in (9, 8, 7, 6, 5, 4, 3):
        sc = 10 ** digits
        for k in (1, 2, 4, 8, 16, 32, 64):
            ulo = Fraction(math.floor(root * sc) - k, sc)
            uhi = Fraction(math.ceil(root * sc) + 1, sc)
            _, up_lo = q_bounds(char, zeros, frac_to_arb(ulo))
            l_hi, _ = q_bounds(char, zeros, frac_to_arb(uhi))
            if up_lo < 1 and l_hi > 1:
                line = ("M0 q=%d: Q(%s) <= %s < 1 < %s <= Q_list(%s)  =>  %s < M_0 < %s   "
                        "(%d listed zeros; float root of Q_list: %.10f)  zeros_sha256=%s"
                        % (args.q, ulo, up_lo.str(12, radius=False), l_hi.str(12, radius=False), uhi,
                           float(ulo), float(uhi), len(zeros), root, meta["sha256"]))
                print(line, flush=True)
                if args.out:
                    with open(args.out, "a") as fh:
                        fh.write(line + "\n")
                return
    raise ArithmeticError("could not certify a bracket for M_0")


if __name__ == "__main__":
    try:
        main()
    except ArithmeticError as exc:
        print("FAIL-CLOSED (m0.py): %s" % exc, flush=True)
        sys.exit(2)
