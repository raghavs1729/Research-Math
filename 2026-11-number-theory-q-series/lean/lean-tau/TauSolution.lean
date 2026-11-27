/-
TauSolution.lean — the UNTRUSTED solution module: it proves, from the library `Tau/`, exactly the
statements of `TauStatement.lean` (same names, same docstrings, same statements; checked by
`scripts/compare_statements.sh`).  No `sorry`, no `axiom`, no `native_decide`.
The definition used in the statements (`heckeSeq`) is that of `Tau/Defs.lean` (Mathlib only),
inlined verbatim in `TauStatement.lean`.
-/
import Tau

namespace Tau

/-- **Theorem 1, first equation.** The equation `x ^ 22 − 3 x ^ 11 y ^ 2 + y ^ 4 = 461` has no
solution in integers `x, y`. -/
theorem no_solution_461 (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ 461 := by
  exact no_solution_461' x y

/-- **Theorem 1, second equation.** The equation `x ^ 22 − 3 x ^ 11 y ^ 2 + y ^ 4 = −599` has no
solution in integers `x, y`. -/
theorem no_solution_neg599 (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ -599 := by
  exact no_solution_M599' x y

/-- **The fourth term of the Hecke sequence.** For all integers `t, p`, the term
`heckeSeq t p 4 = t ^ 4 − 3 p ^ 11 t ^ 2 + p ^ 22` is neither `461` nor `−599`. -/
theorem heckeSeq_four_ne (t p : ℤ) : heckeSeq t p 4 ≠ 461 ∧ heckeSeq t p 4 ≠ -599 := by
  exact heckeSeq_four_ne' t p

/-- **Corollary (`τ(p ^ 4) ∉ {461, −599}`).** Let `τ : ℕ → ℤ` be any function with `τ 1 = 1` that
satisfies the Hecke recurrence `τ (p ^ (k + 2)) = τ p · τ (p ^ (k + 1)) − p ^ 11 · τ (p ^ k)` at the
powers of a prime `p` (Ramanujan's tau function is such a function, for every prime `p`; that fact is
NOT formalized here).  Then `τ (p ^ 4) ≠ 461` and `τ (p ^ 4) ≠ −599`.  (Primality of `p` is not
used.) -/
theorem tau_prime_pow_four_ne (τ : ℕ → ℤ) (h1 : τ 1 = 1) (p : ℕ) (_hp : p.Prime)
    (hrec : ∀ k : ℕ, τ (p ^ (k + 2)) = τ p * τ (p ^ (k + 1)) - (p : ℤ) ^ 11 * τ (p ^ k)) :
    τ (p ^ 4) ≠ 461 ∧ τ (p ^ 4) ≠ -599 := by
  -- primality of `p` is not needed
  rw [eq_heckeSeq τ p h1 hrec 4]
  exact heckeSeq_four_ne (τ p) p

end Tau
