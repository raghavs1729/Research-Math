/-
PiSeriesStatement/SunSqrt2.lean — TRUSTED statement file (imports only Mathlib, no definitions).
The `sorry` is deliberate: this is the statement side.  The proof is in `PiSeriesSolution.lean`.
-/
import Mathlib

/-- **A Ramanujan-type series for `1/π` conjectured by Zhi-Wei Sun**
(Z.-W. Sun, *New series involving binomial coefficients*, arXiv:2307.03086,
Conjecture 5.25(i), eq. (5.20)):

  `∑_{k ≥ 0} (360k³ + 612k² + 230k + 15) · C(3k,k) · C(6k,3k)² / ((3k+1)(3k+2) · 2^{15k}) = 32√2/π`.

In Lean: `Nat.choose (3 * k) k` and `Nat.choose (6 * k) (3 * k)` are the binomial coefficients
`C(3k,k)` and `C(6k,3k)` (natural numbers, cast to `ℝ`), `Real.sqrt 2` is `√2`, and `HasSum f a`
says that the series `∑ f k` over `k : ℕ` converges (unconditionally) to `a`.  All terms are
positive, so this is the same as convergence of the partial sums to `32√2/π`. -/
theorem sun_series_thirty_two_sqrt_two_div_pi :
    HasSum (fun k : ℕ => (360 * (k : ℝ) ^ 3 + 612 * (k : ℝ) ^ 2 + 230 * k + 15)
        * (Nat.choose (3 * k) k : ℝ) * (Nat.choose (6 * k) (3 * k) : ℝ) ^ 2
        / ((3 * (k : ℝ) + 1) * (3 * k + 2) * 2 ^ (15 * k)))
      (32 * Real.sqrt 2 / Real.pi) := by
  sorry
