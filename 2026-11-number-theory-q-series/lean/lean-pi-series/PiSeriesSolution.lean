/-
PiSeriesSolution.lean — the UNTRUSTED solution module: it proves, from the library `PiSeries/`, exactly the
statements of `PiSeriesStatement/*.lean` (same names, same docstrings, same statements; checked by
`scripts/compare_statements.sh`).  No `sorry`, no `axiom`, no `native_decide`.
-/
import PiSeries

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
  exact PiSeries.Ramanujan6.ramanujan_six_n_plus_one

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
  exact PiSeries.Ramanujan42.ramanujan_forty_two_n_plus_five

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
  exact PiSeries.SunSqrt2.sun_series_32_sqrt2_div_pi

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
  exact PiSeries.Bauer.bauer_alternating

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
  exact PiSeries.SunPi2.sun_series_64_sqrt2_div_pi_sq
