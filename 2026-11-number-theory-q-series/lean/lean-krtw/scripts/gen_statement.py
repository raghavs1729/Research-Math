#!/usr/bin/env python3
"""Generate KRTWStatement.lean (the trusted statement module, Mathlib only) from
  * the definitions of KRTW/Defs.lean (the block between `-- BEGIN DEFINITIONS` and
    `-- END DEFINITIONS`, inlined verbatim), and
  * the theorems of KRTWSolution.lean (docstring + statement copied byte for byte, every proof
    replaced by `sorry`).
Run from the project root:  python3 scripts/gen_statement.py"""
import re
defs = open('KRTW/Defs.lean').read()
block = defs[defs.index('-- BEGIN DEFINITIONS'):defs.index('-- END DEFINITIONS') + len('-- END DEFINITIONS')]
sol = open('KRTWSolution.lean').read()
pat = re.compile(r"(/--(?:(?!-/)[\s\S])*-/\ntheorem (\w+)[\s\S]*?) := by\n")
thms = pat.findall(sol)
out = ['''/-
KRTWStatement.lean — the TRUSTED statement module: WHAT IS CLAIMED.   (GENERATED FILE: see
scripts/gen_statement.py; the definitions are inlined verbatim from KRTW/Defs.lean and the statements
are copied byte for byte from KRTWSolution.lean, with every proof replaced by `sorry`.)

This module imports ONLY Mathlib.  The `sorry`s below are deliberate (this is the statement side);
KRTWSolution.lean (untrusted) proves exactly these statements from the library KRTW/.
-/
import Mathlib

open PowerSeries Finset

namespace KRTW

/-! ## Definitions (inlined verbatim from `KRTW/Defs.lean`) -/

''', block, '\n\n/-! ## Statements -/\n']
for t, n in thms:
    out.append('\n' + t + ' := by\n  sorry\n')
out.append('\nend KRTW\n')
open('KRTWStatement.lean', 'w').write(''.join(out))
print('KRTWStatement.lean:', len(thms), 'statements:', ', '.join(n for _, n in thms))
