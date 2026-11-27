#!/usr/bin/env python3
"""
verify_zeros.py q [--count] -- re-certify an assembled zero file from scratch
(independent of the generation state).

 [H] the sha256 of zeros_q<q>.txt equals the pinned value in SHA256SUMS;
 [Z] for every line (label, d): the Hardy Z-function of the labelled character
     (label 'zeta': modulus 1) has certified opposite signs at d - 1e-30 and
     d + 1e-30; the windows of one label are disjoint, increasing and positive.
 [N] (--count) for every conjugate pair, at its T_cert from the header, the
     argument-principle count W equals the number of windows of chi and chibar
     below T_cert (see lzeros_cert.py); for zeta, arb.zeta_nzeros(T_cert).
"""
import re, sys, os, time
from flint import arb, acb, ctx, dirichlet_char
from afcommon import load_zeros, HERE, EPS
from lzeros_cert import Char, sign_of, PREC_REFINE


def main():
    q = int(sys.argv[1])
    t0 = time.time()
    zeros = load_zeros(q)
    ctx.prec = PREC_REFINE
    e = arb(EPS)
    by = {}
    for lab, d in zeros:
        by.setdefault(lab, []).append(d)
    for lab, ds in by.items():
        chi = dirichlet_char(1, 1) if lab == "zeta" else dirichlet_char(q, int(lab))
        prev = None
        for dstr in ds:
            d = arb(dstr)
            assert d - e > 0
            if prev is not None:
                assert prev + e < d - e, "windows not disjoint/increasing (%s, %s)" % (lab, dstr)
            sm = sign_of(chi.hardy_z(acb(d - e)).real)
            sp = sign_of(chi.hardy_z(acb(d + e)).real)
            assert sm != 0 and sm == -sp, "no certified sign change (%s, %s)" % (lab, dstr)
            prev = d
        print("[Z] label %-4s: %d windows certified (%.0fs)" % (lab, len(ds), time.time() - t0), flush=True)
    print("[Z] all %d windows of zeros_q%d.txt certified" % (len(zeros), q))
    if "--count" in sys.argv:
        hdr = [l for l in open(os.path.join(HERE, "zeros_q%d.txt" % q)) if l.startswith("#")]
        for l in hdr:
            m = re.search(r"zeta: T_cert = (\d+)", l)
            if m:
                T = int(m.group(1))
                N = int(arb(T).zeta_nzeros().unique_fmpz())
                K = sum(1 for d in by["zeta"] if arb(d) + e < T)
                assert N == K == len(by["zeta"]), (N, K)
                print("[N] zeta: N(%d) = %d = number of windows" % (T, N))
            m = re.search(r"pair \[([\d, ]+)\].*T_cert = (\d+)/(\d+)", l)
            if m:
                ns = [int(x) for x in m.group(1).split(",")]
                T = arb(int(m.group(2))) / int(m.group(3))
                chars = [Char(q, n) for n in ns]
                D = sum((c.D(T)[0] for c in chars), arb(0))
                if len(ns) == 1:
                    D = 2 * D
                ctx.prec = 128
                W = int((D / arb.pi()).unique_fmpz())
                ctx.prec = PREC_REFINE
                K = sum(1 for n in ns for d in by[str(n)] if arb(d) + e < T) * (2 if len(ns) == 1 else 1)
                Kall = sum(len(by[str(n)]) for n in ns) * (2 if len(ns) == 1 else 1)
                assert W == K == Kall, (ns, W, K, Kall)
                print("[N] pair %s: W(T = %s/%s) = %d = number of windows (with conjugates)" % (ns, m.group(2), m.group(3), W))
    print("verify_zeros: PASS (%.0fs)" % (time.time() - t0))


if __name__ == "__main__":
    main()
