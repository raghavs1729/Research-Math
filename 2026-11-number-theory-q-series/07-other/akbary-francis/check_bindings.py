#!/usr/bin/env python3
"""
check_bindings.py -- confirm that every certificate_*.json is bound to its zero file.

How the binding works (nothing here is new; this script only re-checks it):
 * SHA256SUMS pins zeros_q<q>.txt;
 * afcommon.load_zeros(q), the only reader of zero files used by certify.py,
   verify_zeros.py, f_two_ways.py, compare_primes.py, li_estimate.py, refine_u.py,
   make_scan_input.py and crosscheck_mpmath.py, raises ValueError unless the sha256
   of the file on disk equals the pinned value;
 * certify.py writes that sha256 into the certificate as "zeros_sha256".
Checked here, for each certificate: zeros_sha256 == pinned hash == hash of the file now
on disk.  Exit status 0 iff all certificates pass.
"""
import glob, hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
pins = {}
for line in open(os.path.join(HERE, "SHA256SUMS")):
    parts = line.split()
    if len(parts) == 2:
        pins[parts[1]] = parts[0]
ok = True
for path in sorted(glob.glob(os.path.join(HERE, "certificate_*.json"))):
    c = json.load(open(path))
    name = "zeros_q%d.txt" % c["q"]
    cur = hashlib.sha256(open(os.path.join(HERE, name), "rb").read()).hexdigest()
    good = c.get("zeros_sha256") == pins.get(name) == cur
    ok = ok and good
    print("%-58s q=%-2d a=%-2d certified=%-5s %s  %s" % (os.path.basename(path), c["q"], c["a"], c["certified"],
                                                    name, "BOUND " + cur[:16] + "…" if good else "*** MISMATCH ***"))
print("check_bindings: %s" % ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
