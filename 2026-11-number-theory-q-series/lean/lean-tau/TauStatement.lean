/-
TauStatement.lean — the TRUSTED statement module: WHAT IS CLAIMED.   (GENERATED FILE: see
scripts/gen_statement.py; the definitions are inlined verbatim from Tau/Defs.lean and the statements
are copied byte for byte from TauSolution.lean, with every proof replaced by `sorry`.)

This module imports ONLY Mathlib.  The `sorry`s below are deliberate (this is the statement side);
TauSolution.lean (untrusted) proves exactly these statements from the library Tau/.
-/
import Mathlib

namespace Tau

/-! ## Definitions (inlined verbatim from `Tau/Defs.lean`) -/

-- BEGIN DEFINITIONS
/-- The sequence `a 0 = 1`, `a 1 = t`, `a (k + 2) = t · a (k + 1) − p ^ 11 · a k`.  For a prime `p`
and `t = τ(p)` (Ramanujan's tau function) this is `k ↦ τ(p ^ k)`, by the Hecke recurrence; here `t`
and `p` are arbitrary integers. -/
def heckeSeq (t p : ℤ) : ℕ → ℤ
  | 0 => 1
  | 1 => t
  | k + 2 => t * heckeSeq t p (k + 1) - p ^ 11 * heckeSeq t p k
-- END DEFINITIONS

/-! ## Statements -/

/-- **Theorem 1, first equation.** The equation `x ^ 22 − 3 x ^ 11 y ^ 2 + y ^ 4 = 461` has no
solution in integers `x, y`. -/
theorem no_solution_461 (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ 461 := by
  sorry

/-- **Theorem 1, second equation.** The equation `x ^ 22 − 3 x ^ 11 y ^ 2 + y ^ 4 = −599` has no
solution in integers `x, y`. -/
theorem no_solution_neg599 (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ -599 := by
  sorry

/-- **The fourth term of the Hecke sequence.** For all integers `t, p`, the term
`heckeSeq t p 4 = t ^ 4 − 3 p ^ 11 t ^ 2 + p ^ 22` is neither `461` nor `−599`. -/
theorem heckeSeq_four_ne (t p : ℤ) : heckeSeq t p 4 ≠ 461 ∧ heckeSeq t p 4 ≠ -599 := by
  sorry

/-- **Corollary (`τ(p ^ 4) ∉ {461, −599}`).** Let `τ : ℕ → ℤ` be any function with `τ 1 = 1` that
satisfies the Hecke recurrence `τ (p ^ (k + 2)) = τ p · τ (p ^ (k + 1)) − p ^ 11 · τ (p ^ k)` at the
powers of a prime `p` (Ramanujan's tau function is such a function, for every prime `p`; that fact is
NOT formalized here).  Then `τ (p ^ 4) ≠ 461` and `τ (p ^ 4) ≠ −599`.  (Primality of `p` is not
used.) -/
theorem tau_prime_pow_four_ne (τ : ℕ → ℤ) (h1 : τ 1 = 1) (p : ℕ) (_hp : p.Prime)
    (hrec : ∀ k : ℕ, τ (p ^ (k + 2)) = τ p * τ (p ^ (k + 1)) - (p : ℤ) ^ 11 * τ (p ^ k)) :
    τ (p ^ 4) ≠ 461 ∧ τ (p ^ 4) ≠ -599 := by
  sorry

end Tau
