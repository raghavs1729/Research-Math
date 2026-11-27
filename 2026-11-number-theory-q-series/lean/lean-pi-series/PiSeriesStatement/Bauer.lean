/-
PiSeriesStatement/Bauer.lean — TRUSTED statement file (imports only Mathlib, no definitions).
The `sorry` is deliberate: this is the statement side.  The proof is in `PiSeriesSolution.lean`.
-/
import Mathlib

/-- **Bauer's alternating series for `2/π`** (G. Bauer, *Von den Coefficienten der Reihen von
Kugelfunctionen einer Variablen*, J. reine angew. Math. 56 (1859), 101–121; it was rediscovered
by Ramanujan (first letter to Hardy, 1913) and is the subject of S. B. Ekhad and D. Zeilberger,
*A WZ proof of Ramanujan's formula for π*, arXiv:math/9306213):

  `2/π = 1 − 5 (1/2)³ + 9 (1·3/(2·4))³ − 13 (1·3·5/(2·4·6))³ + ⋯`,

that is, since `1·3⋯(2n−1)/(2·4⋯(2n)) = C(2n,n)/4ⁿ`,

  `∑_{n ≥ 0} (−1)ⁿ C(2n,n)³ (4n+1) / 2^{6n} = 2/π`.

This series converges only conditionally: its terms alternate in sign and have absolute value of
order `n^{-1/2}`, so the series of absolute values diverges and Mathlib's `HasSum` (unconditional
summability) would be false.  The statement is therefore the convergence of the partial sums:
`∑_{n<N} (−1)ⁿ C(2n,n)³ (4n+1) / 2^{6n} → 2/π` as `N → ∞`.

In Lean: `Finset.range N` is `{0, …, N−1}`, `Nat.choose (2 * n) n` is the binomial coefficient
`C(2n,n)` (a natural number, cast to `ℝ`), `2 ^ (6 * n)` is a real number raised to a natural
exponent, and `Filter.Tendsto u Filter.atTop (nhds a)` says that the sequence `u N` tends to `a`. -/
theorem bauer_series_two_div_pi :
    Filter.Tendsto (fun N : ℕ => ∑ n ∈ Finset.range N,
        (-1 : ℝ) ^ n * (Nat.choose (2 * n) n : ℝ) ^ 3 * (4 * n + 1) / 2 ^ (6 * n))
      Filter.atTop (nhds (2 / Real.pi)) := by
  sorry
