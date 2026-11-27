/-
ChebLib/Cheb/Defs.lean — the Mathlib-only definitions used in the statements about Chebyshev's
assertion (in addition to those of ChebLib/Race/Defs.lean).  Nothing here depends on Zeta23.
-/
import ChebLib.Race.Defs

open scoped BigOperators
open Complex

noncomputable section

/-! ## 6. Chebyshev's sums (Abel weights `e^{−n/x}`) -/

/-- **Chebyshev's sum**  `Σ_p χ₋₄(p) e^{−p/x}`  over all primes `p` (`χ₋₄(2) = 0`,
`χ₋₄(p) = (−1)^{(p−1)/2}` for odd `p`; here `χ₋₄` is Mathlib's `ZMod.χ₄`). -/
def chebyshevSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)

/-- The logarithmically weighted sum  `Σ_p χ₋₄(p) (log p) e^{−p/x}`  over all primes `p`. -/
def chebyshevLogSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
    * Real.exp (-((p : ℕ) : ℝ) / x)

/-- The von Mangoldt weighted sum  `Σ_{n ≥ 1} Λ(n) χ₋₄(n) e^{−n/x}`  over prime powers. -/
def chebyshevPsi (x : ℝ) : ℝ :=
  ∑' n : ℕ, ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)
    * Real.exp (-(n : ℝ) / x)

/-- The sum  `Σ_ρ m_ρ |Γ(ρ)|`  over the distinct zeros `ρ` of `L(s,χ)` in the open critical strip
`0 < Re ρ < 1`, each weighted by its multiplicity `m_ρ` (the order of vanishing). -/
def gammaZeroSum {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : ℝ :=
  ∑' ρ : {ρ : ℂ | χ.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1},
    ((analyticOrderAt χ.LFunction (ρ : ℂ)).toNat : ℝ) * ‖Complex.Gamma (ρ : ℂ)‖

end
