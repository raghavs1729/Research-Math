#!/usr/bin/env python3
"""
assemble.py q -- merge state/zeros_q<q>_pair*.json and state/zeros_zeta.json into the
zero file zeros_q<q>.txt, and print its sha256 (to be pinned in SHA256SUMS and in
RESULT.md / certify.py; cf. audit finding G1 of monotonicity/nonmonotone-audit).

File format (ASCII): header lines '#...'; then one line per window
     <label> <d>
label = 'zeta' or the Conrey number n of chi_q(n, .); the window is (d - 1e-30, d + 1e-30)
and contains a zero 1/2 + i*gamma of zeta resp. L(s, chi_q(n, .)), gamma > 0.
(The zeros with negative ordinate are the conjugates of the zeros of the conjugate
character chi_q(n^-1 mod q, .).)
"""
import glob, hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
q = int(sys.argv[1])
lines = []
hdr = ["# zeros_q%d.txt: windows (d-1e-30, d+1e-30), each certified to contain exactly one simple zero" % q,
       "# 1/2 + i*gamma (gamma > 0) of the function named by the label; complete below T_cert (per pair):"]
zz = json.load(open(os.path.join(HERE, "state", "zeros_zeta.json")))
hdr.append("#   zeta: T_cert = %d, %d zeros" % (zz["T"], len(zz["zeros"])))
for d in zz["zeros"]:
    lines.append("zeta %s" % d)
Tmin = zz["T"]
for path in sorted(glob.glob(os.path.join(HERE, "state", "zeros_q%d_pair*.json" % q)),
                   key=lambda p: int(p.split("pair")[1].split(".")[0])):
    st = json.load(open(path))
    num, den = st["top_index"], st["den"]
    Tmin = min(Tmin, num / den)
    hdr.append("#   pair %s (parity %d): T_cert = %d/%d, zeros: %s; last count W = %d"
               % (st["pair"], st["parity"], num, den,
                  ", ".join("%s:%d" % (k, len(v)) for k, v in st["zeros"].items()), st["blocks"][-1]["W"]))
    for k, v in st["zeros"].items():
        for d in v:
            lines.append("%s %s" % (k, d))
hdr.append("# total %d windows; common completeness height %.5f" % (len(lines), Tmin))
out = os.path.join(HERE, "zeros_q%d.txt" % q)
with open(out, "w", encoding="ascii") as fh:
    fh.write("\n".join(hdr + lines) + "\n")
h = hashlib.sha256(open(out, "rb").read()).hexdigest()
print("%s  zeros_q%d.txt   (%d windows, common T = %.5f)" % (h, q, len(lines), Tmin))
