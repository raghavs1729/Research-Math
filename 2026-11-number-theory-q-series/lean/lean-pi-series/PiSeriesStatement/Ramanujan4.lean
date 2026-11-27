/-
PiSeriesStatement/Ramanujan4.lean — TRUSTED statement file (imports only Mathlib, no definitions).
The `sorry` is deliberate: this is the statement side.  The proof is in `PiSeriesSolution.lean`.
-/
import Mathlib

/-- **Ramanujan's series for `4/π`** (S. Ramanujan, *Modular equations and approximations to π*,
Quart. J. Math. 45 (1914), 350–372, eq. (28)):

  `4/π = 1 + (7/4)(1/2)³ + (13/4²)(1·3/(2·4))³ + (19/4³)(1·3·5/(2·4·6))³ + ⋯`,

that is, since `1·3⋯(2n−1)/(2·4⋯(2n)) = C(2n,n)/4ⁿ`,

  `∑_{n ≥ 0} C(2n,n)³ (6n+1) / 2^{8n} = 4/π`.

In Lean: `Nat.choose (2 * n) n` is the binomial coefficient `C(2n,n)` (a natural number, cast to
`ℝ`), `2 ^ (8 * n)` is a real number raised to a natural exponent, and `HasSum f a` says that the
series `∑ f n` over `n : ℕ` converges (unconditionally) to `a`.  All terms are positive, so this
is the same as convergence of the partial sums `∑_{n<N}` to `4/π`. -/
theorem ramanujan_series_four_div_pi :
    HasSum (fun n : ℕ => (Nat.choose (2 * n) n : ℝ) ^ 3 * (6 * n + 1) / 2 ^ (8 * n))
      (4 / Real.pi) := by
  sorry
