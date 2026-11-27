#!/usr/bin/env python3
"""
count_crosscheck.py -- test of the symmetric zero count of zeros_cert.py
(Cert.count) against the previously certified data for chi_4 and chi_3:
prototype/explore/certify_arb.py --K 0 proved, by a full-rectangle winding
number of L (a different method), that the 511 (resp. 537) listed ordinates
are ALL the zeros below T = 640.165707 (resp. 700.024747).  So Cert.count
must return exactly 1022 (resp. 1074) at those heights.

Usage:  nice -n 10 python3 count_crosscheck.py [--out file]
"""
import os
import sys
import time

from flint import arb, ctx

from zeros_cert import Cert
from common import read_zero_windows

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.join(HERE, "..", "..", "prototype", "explore")
out = open(sys.argv[2], "a") if len(sys.argv) > 2 and sys.argv[1] == "--out" else None
allok = True
for q, fn, Ts in ((4, "chi4_zeros.txt", "640.165707"), (3, "chi3_zeros.txt", "700.024747")):
    zeros = read_zero_windows(os.path.join(EXP, fn))
    C = Cert(q)
    for T in (arb(Ts), (zeros[49] + zeros[50]) / 2):
        t0 = time.time()
        T = arb(T.mid())
        K = sum(1 for g in zeros if g < T)
        amb = sum(1 for g in zeros if not (g < T) and not (g > T))
        W, th, vert, hor, nev = C.count(T)
        n = int(W.unique_fmpz())
        ok = (n == 2 * K and amb == 0)
        allok = allok and ok
        line = ("COUNT q=%d T=%s: symmetric count W = %s -> %d ; listed ordinates below T: %d ; "
                "2K = %d : %s  [%d L-evals, %.0fs]"
                % (q, T.str(12, radius=False), W.str(12, radius=True), n, K, 2 * K,
                   "AGREES" if ok else "*** DISAGREES ***", nev, time.time() - t0))
        print(line, flush=True)
        if out:
            out.write(line + "\n")
            out.flush()
sys.exit(0 if allok else 1)
