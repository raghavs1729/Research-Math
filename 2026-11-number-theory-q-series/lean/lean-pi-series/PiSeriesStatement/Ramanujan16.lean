/-
PiSeriesStatement/Ramanujan16.lean — TRUSTED statement file (imports only Mathlib, no definitions).
The `sorry` is deliberate: this is the statement side.  The proof is in `PiSeriesSolution.lean`.
-/
import Mathlib

/-- **Ramanujan's series for `16/π`** (S. Ramanujan, *Modular equations and approximations to π*,
Quart. J. Math. 45 (1914), 350–372, eq. (29)):

  `16/π = 5 + (47/64)(1/2)³ + (89/64²)(1·3/(2·4))³ + (131/64³)(1·3·5/(2·4·6))³ + ⋯`,

that is, since `1·3⋯(2n−1)/(2·4⋯(2n)) = C(2n,n)/4ⁿ`,

  `∑_{n ≥ 0} C(2n,n)³ (42n+5) / 2^{12n} = 16/π`.

In Lean: `Nat.choose (2 * n) n` is the binomial coefficient `C(2n,n)` (a natural number, cast to
`ℝ`), `2 ^ (12 * n)` is a real number raised to a natural exponent, and `HasSum f a` says that the
series `∑ f n` over `n : ℕ` converges (unconditionally) to `a`.  All terms are positive, so this
is the same as convergence of the partial sums `∑_{n<N}` to `16/π`. -/
theorem ramanujan_series_sixteen_div_pi :
    HasSum (fun n : ℕ => (Nat.choose (2 * n) n : ℝ) ^ 3 * (42 * n + 5) / 2 ^ (12 * n))
      (16 / Real.pi) := by
  sorry
