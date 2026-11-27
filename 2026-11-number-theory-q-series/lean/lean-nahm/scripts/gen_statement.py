#!/usr/bin/env python3
"""Generate NahmStatement.lean (the trusted statement module, Mathlib only) from
  * the definitions of Nahm/Defs.lean (the block between `-- BEGIN DEFINITIONS` and
    `-- END DEFINITIONS`, inlined verbatim), and
  * the theorems of NahmSolution.lean (docstring + statement copied byte for byte, every proof
    replaced by `sorry`).
Run from the project root:  python3 scripts/gen_statement.py"""
import re
defs = open('Nahm/Defs.lean').read()
block = defs[defs.index('-- BEGIN DEFINITIONS'):defs.index('-- END DEFINITIONS') + len('-- END DEFINITIONS')]
sol = open('NahmSolution.lean').read()
pat = re.compile(r"(/--(?:(?!-/)[\s\S])*-/\ntheorem (\w+)[\s\S]*?) := by\n")
thms = pat.findall(sol)
out = ['''/-
NahmStatement.lean — the TRUSTED statement module: WHAT IS CLAIMED.   (GENERATED FILE: see
scripts/gen_statement.py; the definitions are inlined verbatim from Nahm/Defs.lean and the statements
are copied byte for byte from NahmSolution.lean, with every proof replaced by `sorry`.)

This module imports ONLY Mathlib.  The `sorry`s below are deliberate (this is the statement side);
NahmSolution.lean (untrusted) proves exactly these statements from the library Nahm/.
-/
import Mathlib

open PowerSeries Finset

namespace Nahm

/-! ## Definitions (inlined verbatim from `Nahm/Defs.lean`) -/

''', block, '\n\n/-! ## Statements -/\n']
for t, n in thms:
    out.append('\n' + t + ' := by\n  sorry\n')
out.append('\nend Nahm\n')
open('NahmStatement.lean', 'w').write(''.join(out))
print('NahmStatement.lean:', len(thms), 'statements:', ', '.join(n for _, n in thms))
