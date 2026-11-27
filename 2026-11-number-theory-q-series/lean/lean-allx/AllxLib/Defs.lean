/-
AllxLib/Defs.lean — the Mathlib-only definitions used in the statements of AllxChallenge.lean in
addition to `chi4` and `GRHFor` (which are those of lean-chebyshev, `ChebLib/Race/Defs.lean`).
Nothing here depends on Zeta23 or ChebLib.
-/
import Mathlib

/-- `π(t; 4, a)`: the number of primes `p ≤ t` with `p ≡ a (mod 4)` (for a real cut-off `t`). -/
noncomputable def primeCountingMod4 (a : ℕ) (t : ℝ) : ℕ :=
  ((Finset.Icc 0 ⌊t⌋₊).filter (fun p => p.Prime ∧ p % 4 = a)).card
