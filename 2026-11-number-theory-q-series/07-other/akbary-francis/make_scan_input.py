#!/usr/bin/env python3
"""make_scan_input.py q a [out] -- 'gamma phase' lines for scan.c from the pinned zero file."""
import sys, math
from afcommon import load_zeros, chi_arg_table
q, a = int(sys.argv[1]), int(sys.argv[2])
tab = chi_arg_table(q, a)
out = sys.argv[3] if len(sys.argv) > 3 else "scan_in_q%d_a%d.txt" % (q, a)
with open(out, "w") as fh:
    for lab, d in load_zeros(q):
        fh.write("%.17g %.17g\n" % (float(d), 2 * math.pi * float(tab[lab])))
