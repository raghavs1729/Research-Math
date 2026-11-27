/-
PiSeriesStatement/SunPiSq.lean — TRUSTED statement file (imports only Mathlib, no definitions).
The `sorry` is deliberate: this is the statement side.  The proof is in `PiSeriesSolution.lean`.
-/
import Mathlib

/-- **A Ramanujan-type series for `1/π²` conjectured by Zhi-Wei Sun**
(Z.-W. Sun, *New series involving binomial coefficients*, arXiv:2307.03086,
Conjecture 5.33(i), eq. (5.31)):

  `∑_{k ≥ 0} (1344k³ + 944k² + 156k + 9) · C(2k,k)² · C(4k,2k)³ / ((2k+1) · 2^{22k}) = 64√2/π²`.

In Lean: `Nat.choose (2 * k) k` and `Nat.choose (4 * k) (2 * k)` are the binomial coefficients
`C(2k,k)` and `C(4k,2k)` (natural numbers, cast to `ℝ`), `Real.sqrt 2` is `√2`, `Real.pi ^ 2` is
`π²`, and `HasSum f a` says that the series `∑ f k` over `k : ℕ` converges (unconditionally) to
`a`.  All terms are positive, so this is the same as convergence of the partial sums to
`64√2/π²`. -/
theorem sun_series_sixty_four_sqrt_two_div_pi_sq :
    HasSum (fun k : ℕ => (1344 * (k : ℝ) ^ 3 + 944 * (k : ℝ) ^ 2 + 156 * k + 9)
        * (Nat.choose (2 * k) k : ℝ) ^ 2 * (Nat.choose (4 * k) (2 * k) : ℝ) ^ 3
        / ((2 * (k : ℝ) + 1) * 2 ^ (22 * k)))
      (64 * Real.sqrt 2 / Real.pi ^ 2) := by
  sorry
