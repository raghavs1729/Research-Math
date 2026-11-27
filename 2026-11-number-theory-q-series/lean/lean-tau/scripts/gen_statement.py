#!/usr/bin/env python3
"""Generate TauStatement.lean (the trusted statement module, Mathlib only) from
  * the definitions of Tau/Defs.lean (the block between `-- BEGIN DEFINITIONS` and
    `-- END DEFINITIONS`, inlined verbatim), and
  * the theorems of TauSolution.lean (docstring + statement copied byte for byte, every proof
    replaced by `sorry`).
Run from the project root:  python3 scripts/gen_statement.py"""
import re
defs = open('Tau/Defs.lean').read()
block = defs[defs.index('-- BEGIN DEFINITIONS'):defs.index('-- END DEFINITIONS') + len('-- END DEFINITIONS')]
sol = open('TauSolution.lean').read()
pat = re.compile(r"(/--(?:(?!-/)[\s\S])*-/\ntheorem (\w+)[\s\S]*?) := by\n")
thms = pat.findall(sol)
out = ['''/-
TauStatement.lean — the TRUSTED statement module: WHAT IS CLAIMED.   (GENERATED FILE: see
scripts/gen_statement.py; the definitions are inlined verbatim from Tau/Defs.lean and the statements
are copied byte for byte from TauSolution.lean, with every proof replaced by `sorry`.)

This module imports ONLY Mathlib.  The `sorry`s below are deliberate (this is the statement side);
TauSolution.lean (untrusted) proves exactly these statements from the library Tau/.
-/
import Mathlib



namespace Tau

/-! ## Definitions (inlined verbatim from `Tau/Defs.lean`) -/

''', block, '\n\n/-! ## Statements -/\n']
for t, n in thms:
    out.append('\n' + t + ' := by\n  sorry\n')
out.append('\nend Tau\n')
open('TauStatement.lean', 'w').write(''.join(out))
print('TauStatement.lean:', len(thms), 'statements:', ', '.join(n for _, n in thms))
