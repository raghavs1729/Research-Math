/-
Tau/Defs.lean — the definitions used in the statements (imports only Mathlib).  The marked block
below is inlined verbatim in `TauStatement.lean`.
-/
import Mathlib

namespace Tau

-- BEGIN DEFINITIONS
/-- The sequence `a 0 = 1`, `a 1 = t`, `a (k + 2) = t · a (k + 1) − p ^ 11 · a k`.  For a prime `p`
and `t = τ(p)` (Ramanujan's tau function) this is `k ↦ τ(p ^ k)`, by the Hecke recurrence; here `t`
and `p` are arbitrary integers. -/
def heckeSeq (t p : ℤ) : ℕ → ℤ
  | 0 => 1
  | 1 => t
  | k + 2 => t * heckeSeq t p (k + 1) - p ^ 11 * heckeSeq t p k
-- END DEFINITIONS

end Tau
