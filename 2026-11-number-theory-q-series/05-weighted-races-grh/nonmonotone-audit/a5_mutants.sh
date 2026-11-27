#!/bin/sh
# a5_mutants.sh -- mutation tests on COPIES of the audited code (audit item 4).
# Each mutant lives in mutants/<name>/ ; the originals are never touched.
set -u
A=/home/user/weighted-prime-races/monotonicity/nonmonotone-audit
S=/home/user/weighted-prime-races/monotonicity/nonmonotone
mk() { rm -rf "$A/mutants/$1"; mkdir -p "$A/mutants/$1"; cp "$S"/*.py "$S"/zeros_q163.txt "$S"/zeros_q163.cert.txt "$A/mutants/$1/"; }
run() { ( cd "$A/mutants/$1" && nice -n 10 python3 density.py --q 163 --M 3/5 6/5 2>&1 | tail -3 ); echo "exit=$?"; }

echo "== m1: delete gamma_1 from the zero file"
mk m1; sed -i '3d' "$A/mutants/m1/zeros_q163.txt"; run m1

echo "== m2: delete gamma_67 (mid-list) from the zero file"
mk m2; sed -i '69d' "$A/mutants/m2/zeros_q163.txt"; run m2

echo "== m3: tail model sign flip: exp(-V t^2/2) -> exp(+V t^2/2) in the main integrand"
mk m3; sed -i 's|g = (-self.V \* t \* t / 2).exp()|g = (self.V * t * t / 2).exp()|' "$A/mutants/m3/density.py"
grep -n "g = (self.V" "$A/mutants/m3/density.py"; run m3

echo "== m4: shrink a radius: drop E1 from the final radius"
mk m4; sed -i 's|raw_radius = E1.upper() + E2.upper()|raw_radius = 0 * E1.upper() + E2.upper()|' "$A/mutants/m4/density.py"
grep -n "raw_radius = 0" "$A/mutants/m4/density.py"; run m4

echo "== m5: zeros_cert.py with one detected sign change dropped (scan 'misses' a zero), H=20"
mk m5; rm -f "$A/mutants/m5/zeros_q163.txt" "$A/mutants/m5/zeros_q163.cert.txt"
sed -i 's|                new.append(dstr)|                if len(new) != 5: new.append(dstr)|' "$A/mutants/m5/zeros_cert.py"
grep -n "len(new) != 5" "$A/mutants/m5/zeros_cert.py"
( cd "$A/mutants/m5" && nice -n 10 python3 zeros_cert.py 163 20 2>&1 | tail -3 ); echo "exit=$?"

echo "== m6: theta(T) sign error in zeros_cert.py (Im lgamma subtracted), H=20"
mk m6; rm -f "$A/mutants/m6/zeros_q163.txt" "$A/mutants/m6/zeros_q163.cert.txt"
sed -i 's|return T / 2 \* qpi + acb(arb(3) / 4, T / 2).lgamma().imag|return T / 2 * qpi - acb(arb(3) / 4, T / 2).lgamma().imag|' "$A/mutants/m6/zeros_cert.py"
grep -n "qpi - acb" "$A/mutants/m6/zeros_cert.py"
( cd "$A/mutants/m6" && nice -n 10 python3 zeros_cert.py 163 20 2>&1 | tail -3 ); echo "exit=$?"
