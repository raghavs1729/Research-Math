#!/usr/bin/env bash
# Statement comparison (a weak substitute for github.com/leanprover/comparator, which is not
# installed here).  Checks:
#  1. KRTWStatement.lean imports only Mathlib; its definition block is byte for byte the block of
#     KRTW/Defs.lean (which imports only Mathlib); every theorem of KRTWStatement.lean (docstring +
#     statement up to `:= by`) occurs byte for byte in KRTWSolution.lean; the solution side
#     (KRTWSolution.lean, KRTW.lean, KRTW/*.lean) has no `sorry`, `axiom`, `native_decide`, `unsafe`,
#     `implemented_by`, `extern`, and KRTWStatement.lean has no `axiom`;
#  2. elaborated: `#check @name` (theorems) and `#print name` (definitions) give the same output from
#     both modules.
# Run from the project root after `lake build`.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.elan/bin:$PATH"
names=$(grep -oP '^theorem \K\w+' KRTWStatement.lean)
defs=$(grep -oP '^(noncomputable )?def \K\w+' KRTWStatement.lean)
python3 - <<'PY'
import re,sys,glob
st=open('KRTWStatement.lean').read()
so=open('KRTWSolution.lean').read()
df=open('KRTW/Defs.lean').read()
for name,src in (('KRTWStatement.lean',st),('KRTW/Defs.lean',df)):
    imports=re.findall(r'^import (.*)$', src, re.M)
    if imports!=['Mathlib']: sys.exit(f"{name}: imports {imports}, expected only Mathlib")
blk=lambda s: s[s.index('-- BEGIN DEFINITIONS'):s.index('-- END DEFINITIONS')]
if blk(st)!=blk(df): sys.exit("definition blocks of KRTWStatement.lean and KRTW/Defs.lean differ")
pat=re.compile(r"(/--(?:(?!-/)[\s\S])*-/\ntheorem (\w+)[\s\S]*?) := by\n")
thms=pat.findall(st)
if not thms: sys.exit("no theorems found in KRTWStatement.lean")
bad=[n for t,n in thms if t not in so]
if bad: sys.exit("textual MISMATCH: "+", ".join(bad))
strip=lambda s: re.sub(r'--.*','',re.sub(r'/-(?:.|\n)*?-/','',s))
if re.search(r'\baxiom\b', strip(st)): sys.exit("KRTWStatement.lean contains `axiom`")
for f in ['KRTWSolution.lean','KRTW.lean']+sorted(glob.glob('KRTW/**/*.lean', recursive=True)):
    src=strip(open(f).read())
    for b in ('sorry','native_decide','axiom ','unsafe ','implemented_by','extern'):
        if b in src: sys.exit(f"{f} contains `{b.strip()}`")
print(f"OK (textual): the definition block and the {len(thms)} statements of KRTWStatement.lean occur verbatim in KRTW/Defs.lean and KRTWSolution.lean; both import only Mathlib resp. the library; no sorry/axiom/native_decide on the solution side")
PY
tmp=$(mktemp -d)
for side in KRTWStatement KRTWSolution; do
  {
    echo "import $side"
    echo "open KRTW"
    for n in $names; do echo "#check @KRTW.$n"; done
    for d in $defs; do echo "#print KRTW.$d"; done
  } > "$tmp/$side.lean"
  lake env lean "$tmp/$side.lean" > "$tmp/$side.out"
done
if diff "$tmp/KRTWStatement.out" "$tmp/KRTWSolution.out" > /dev/null; then
  echo "OK (elaborated): the $(echo $names | wc -w) statements and the $(echo $defs | wc -w) definitions are identical in KRTWStatement and KRTWSolution:"
  n=$(echo $names | wc -w)
  awk -v n="$n" '/^(@KRTW|KRTW)\./{c++} c<=n' "$tmp/KRTWSolution.out"
else
  echo "MISMATCH between Statement and Solution:"; diff "$tmp/KRTWStatement.out" "$tmp/KRTWSolution.out"; exit 1
fi
