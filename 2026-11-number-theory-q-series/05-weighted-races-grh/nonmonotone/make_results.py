#!/usr/bin/env python3
"""
make_results.py -- assemble results.txt from the raw outputs and CHECK the
claimed inequalities between the printed (outward-rounded, 12-decimal)
enclosure endpoints in exact rational arithmetic.  No floating point enters
any comparison that is labelled CERTIFIED.

Usage:  python3 make_results.py        (run after the other scripts)
"""
from __future__ import annotations

import os
import re
import sys
from fractions import Fraction

from common import HZ_SHA256, sha256_file, find_hash_record

HERE = os.path.dirname(os.path.abspath(__file__))


def read(name):
    p = os.path.join(HERE, name)
    return open(p).read() if os.path.exists(p) else ""


def deltas(text, hashes=None):
    """Only lines that START with DELTA count (UNCERTIFIED lines are ignored).  If `hashes`
    is a set, the zeros_sha256 of every line is added to it; a DELTA line without a hash
    (output of the pre-audit density.py) is rejected."""
    out = {}
    for m in re.finditer(r"^DELTA q=(\d+) M=(\S+) \(m=(\S+)\)\s+lower=(\S+)\s+upper=(\S+)\s+width=(\S+)(.*)$",
                         text, re.M):
        h = re.search(r"zeros_sha256=([0-9a-f]{64})", m.group(7))
        if hashes is not None:
            hashes.add(h.group(1) if h else "MISSING")
        out[(int(m.group(1)), Fraction(m.group(2)))] = (Fraction(m.group(4)), Fraction(m.group(5)),
                                                         m.group(4), m.group(5))
    return out


def crossvals(text):
    """{M: {'A': x, 'B': x, 'C': (mean, se)}} from crosscheck.py output."""
    out = {}
    cur = None
    for line in text.splitlines():
        m = re.match(r"M=(\S+): sigma\^2", line)
        if m:
            cur = Fraction(m.group(1))
            out.setdefault(cur, {})
            continue
        m = re.search(r"\((A|B)\) .*delta = ([0-9.]+)", line)
        if m and cur is not None:
            out[cur][m.group(1)] = Fraction(m.group(2))
            continue
        m = re.match(r"\s+M=(\S+): delta = ([0-9.]+) \+- ([0-9.]+)", line)
        if m:
            out.setdefault(Fraction(m.group(1)), {})["C"] = (Fraction(m.group(2)), Fraction(m.group(3)))
    return out


def main():
    L = []
    ok = True

    def say(s=""):
        L.append(s)

    def check(cond, text):
        nonlocal ok
        ok = ok and cond
        say("  [%s] %s" % ("OK" if cond else "FAILED", text))

    say("results.txt -- non-monotonicity of the weighted bias density for the real characters")
    say("of conductor 163 and 232.  Assembled by make_results.py from the raw outputs reproduced")
    say("verbatim below.  All enclosure endpoints are outward-rounded 12-decimal numbers; the")
    say("comparisons marked CERTIFIED are exact rational comparisons of those endpoints.")
    say("Every enclosure is tied to its zero file by sha256 (hypothesis (H-Z), audit finding G1).")
    say("M = m + 1/2;  delta(M) = P(1 + sum_gamma 4M/sqrt(M^2+gamma^2) cos(theta_gamma) > 0).")
    say("See NONMONOTONE.md for the theorem, the hypotheses and what is / is not rigorous.")
    say()
    for q, pairs, shape in (
            (163, [("3/5", "6/5"), ("57/100", "117/100"), ("57/100", "3/5")],
             [("3/10", "57/100"), ("117/100", "6/5"), ("6/5", "2")]),
            (232, [("18/25", "1"), ("71/100", "101/100"), ("71/100", "18/25"), ("1", "101/100")],
             [("2/5", "71/100"), ("101/100", "3/2")])):
        hs = set()
        D = deltas(read("raw_density_q%d.txt" % q), hs)
        if not D:
            say("q = %d: no density output found" % q)
            continue
        say("=" * 100)
        say("CONDUCTOR %d   (chi = Kronecker symbol (-%d/.), odd)" % (q, q))
        say("=" * 100)
        say("Zero file binding, hypothesis (H-Z):")
        zp = os.path.join(HERE, "zeros_q%d.txt" % q)
        actual = sha256_file(zp)
        rec = find_hash_record(q, actual)
        check(actual == HZ_SHA256[q], "sha256(zeros_q%d.txt) = %s equals the hash pinned by the audit" % (q, actual))
        check(rec is not None and rec["W"] == 2 * rec["K"],
              "passing HASH-RECORD for that hash in zeros_q%d.cert.txt: K = %s, T_cert = %s, W = %s"
              % ((q,) + ((rec["K"], rec["T_cert"], rec["W"]) if rec else ("-", "-", "-"))))
        check(hs == {HZ_SHA256[q]}, "every DELTA line of raw_density_q%d.txt carries zeros_sha256 = that hash" % q)
        for line in read("raw_m0.txt").splitlines():
            if line.startswith("M0 q=%d:" % q):
                check(("zeros_sha256=" + HZ_SHA256[q]) in line, "the M_0 line carries the same hash")
        say("CERTIFIED enclosures (GRH; zero data certified by zeros_cert.py):")
        for (qq, M) in sorted(D):
            lo, hi, slo, shi = D[(qq, M)]
            say("   M = %-8s (m = %-7s)   delta in [ %s , %s ]   width %.2e"
                % (M, M - Fraction(1, 2), slo, shi, float(hi - lo)))
        say("CERTIFIED comparisons:")
        for a, b in pairs:
            A, B = D.get((q, Fraction(a))), D.get((q, Fraction(b)))
            if A and B:
                check(A[1] < B[0], "INCREASE: delta(%s) <= %s < %s <= delta(%s);  delta(%s) - delta(%s) in [%s, %s]"
                      % (a, A[3], B[2], b, b, a, "%.12f" % float(B[0] - A[1]), "%.12f" % float(B[1] - A[0])))
        for a, b in shape:
            A, B = D.get((q, Fraction(a))), D.get((q, Fraction(b)))
            if A and B:
                check(A[0] > B[1], "decrease: delta(%s) >= %s > %s >= delta(%s)" % (a, A[2], B[3], b))
        for line in read("raw_m0.txt").splitlines():
            m = re.match(r"M0 q=%d:.*=>\s+(\S+) < M_0 < (\S+)" % q, line)
            if m:
                check(Fraction(m.group(2)) < min(Fraction(a) for a, _ in pairs),
                      "M_0 in (%s, %s): the increase lies above M_0 (consistent with Theorem B)"
                      % (m.group(1), m.group(2)))
        X = crossvals(read("raw_crosscheck_q%d.txt" % q))
        if X:
            say("NON-RIGOROUS cross-checks against the certified enclosures:")
            for M in sorted(X):
                if (q, M) not in D:
                    continue
                lo, hi = D[(q, M)][:2]
                for key, label in (("A", "mpmath tanh-sinh"), ("B", "Simpson")):
                    if key in X[M]:
                        v = X[M][key]
                        # these use the Gaussian surrogate for the unlisted zeros, which is
                        # the upper end of the certified two-sided tail model: they must lie
                        # in the enclosure
                        check(lo <= v <= hi, "M = %s: %s value %s lies in the enclosure" % (M, label, float(v)))
                if "C" in X[M]:
                    mean, se = X[M]["C"]
                    check(lo - 4 * se <= mean <= hi + 4 * se,
                          "M = %s: Monte Carlo %.6f +- %.6f is within 4 s.e. of the enclosure (%.2f s.e.)"
                          % (M, float(mean), float(se), float(abs(mean - (lo + hi) / 2) / se)))
        say()

    # sanity: chi mod 4
    hs4 = set()
    S = deltas(read("raw_sanity_chi4.txt"), hs4)
    if S:
        say("=" * 100)
        say("SANITY: the same code on chi mod 4 (zeros and published enclosures from prototype/explore)")
        say("=" * 100)
        p4 = os.path.join(HERE, "..", "..", "prototype", "explore", "chi4_zeros.txt")
        a4 = sha256_file(p4)
        r4 = find_hash_record(4, a4)
        check(r4 is not None and hs4 == {a4},
              "chi4_zeros.txt (sha256 %s) re-certified by zeros_cert.py --verify (K = %s, W = %s) and "
              "printed in every DELTA line" % (a4, r4["K"] if r4 else "-", r4["W"] if r4 else "-"))
        pub = {Fraction(1, 2): ("0.995927866782", "0.995927983908", "0.995928"),
               Fraction(3, 2): ("0.797630378395", "0.797630419780", "0.797630")}
        for M, (plo, phi, six) in pub.items():
            if (4, M) in S:
                lo, hi, slo, shi = S[(4, M)]
                check(Fraction(plo) <= lo and hi <= Fraction(phi),
                      "m = %s: [%s, %s] is contained in the published enclosure [%s, %s]"
                      % (M - Fraction(1, 2), slo, shi, plo, phi))
                r = Fraction(six)
                check(r - Fraction(5, 10 ** 7) < lo and hi < r + Fraction(5, 10 ** 7),
                      "m = %s: every point of the enclosure rounds to %s" % (M - Fraction(1, 2), six))
        say()
    hst = set()
    T = deltas(read("raw_stability_q163.txt"), hst)
    D = deltas(read("raw_density_q163.txt"))
    if T and D:
        say("=" * 100)
        say("STABILITY: conductor 163 with only the first 60 zeros (the other 73 treated as unlisted)")
        say("=" * 100)
        check(hst == {HZ_SHA256[163]}, "same certified zero file (hash in every DELTA line)")
        for key in sorted(T):
            if key in D:
                check(T[key][0] <= D[key][0] and D[key][1] <= T[key][1],
                      "M = %s: 60-zero enclosure [%s, %s] contains the 133-zero enclosure"
                      % (key[1], T[key][2], T[key][3]))
        say()
    W = read("raw_count_crosscheck.txt")
    if "COUNT" in W:
        say("=" * 100)
        say("ZERO-COUNT ROUTINE TESTED on chi_4, chi_3 (answers certified independently by certify_arb.py)")
        say("=" * 100)
        for line in W.splitlines():
            if line.startswith("COUNT"):
                check("AGREES" in line, line)
        say()
    X = read("raw_mutant_refusal.txt")
    if X:
        say("=" * 100)
        say("MUTATION DEMONSTRATION (audit G1): incomplete zero files are refused")
        say("=" * 100)
        blocks = [blk for blk in X.split("$ ")[1:]]
        for blk in blocks:
            cmd = blk.splitlines()[0]
            st = re.search(r"\[exit status (\d+)\]", blk)
            code = int(st.group(1)) if st else -1
            has_delta = any(l.startswith("DELTA") for l in blk.splitlines())
            has_rec = any(l.startswith("HASH-RECORD") for l in blk.splitlines())
            if "--uncertified-zeros" in cmd and code == 0:
                check(not has_delta and "UNCERTIFIED q=" in blk,
                      "`%s`: runs only with the waiver, and every output line is labelled UNCERTIFIED" % cmd)
            else:
                check(code != 0 and not has_delta and not has_rec and ("REFUSED" in blk or "VERIFY FAILED" in blk),
                      "`%s`: refused (exit status %d), no DELTA line, no hash record" % (cmd, code))
        say()
    say("OVERALL: %s" % ("ALL CHECKS PASSED" if ok else "*** SOME CHECK FAILED ***"))
    say()
    say("#" * 100)
    say("# RAW OUTPUTS (verbatim)")
    say("#" * 100)
    for name in ("zeros_q163.cert.txt", "raw_density_q163.txt", "raw_crosscheck_q163.txt",
                 "zeros_q232.cert.txt", "raw_density_q232.txt", "raw_crosscheck_q232.txt",
                 "raw_m0.txt", "zeros_q4.cert.txt", "zeros_q3.cert.txt", "raw_sanity_chi4.txt",
                 "raw_stability_q163.txt", "raw_count_crosscheck.txt", "raw_mutant_refusal.txt"):
        txt = read(name)
        if txt:
            say()
            say("----- %s " % name + "-" * (90 - len(name)))
            L.extend(txt.rstrip("\n").splitlines())
    with open(os.path.join(HERE, "results.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L[:L.index("#" * 100)]))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
