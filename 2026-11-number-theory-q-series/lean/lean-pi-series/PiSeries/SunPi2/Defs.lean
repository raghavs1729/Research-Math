/-
# Sun's series  `∑_{k≥0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^{22k}) = 64√2/π²`
# (arXiv:2307.03086, Conjecture 5.33(i)).  Shared definitions.

The proof (RESULTS.md, Theorem 11, reorganised in `N2-PLAN.md`) uses a WZ pair supported on the
triangle `0 ≤ k < n`.  Its first member `a(n,k)` is defined here by its two ratios:

  `a(1,0) = 9/256`,  `a(n+1,0) = ρ(n,0) a(n,0)`,  `a(n,k+1) = σ(n,k) a(n,k)`,

so that `a(n,k) = 0` for `k ≥ n ≥ 1` (because `σ(n,n-1) = 0`).  `T N = ∑_{k<N} a(N,k)`.
-/
import Mathlib

open Finset

namespace PiSeries.SunPi2

/-- `a(n,k+1)/a(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ :=
  (2 * (k : ℝ) + 3) * ((n : ℝ) - k - 1) ^ 2
      * ((8 * (k : ℝ) + 4 * n + 1) * (8 * (k : ℝ) + 4 * n + 3) * (8 * (k : ℝ) + 4 * n + 5)
        * (8 * (k : ℝ) + 4 * n + 7)) /
    ((2 * (k : ℝ) + 1) * ((n : ℝ) + k + 1) ^ 2
      * ((8 * (k : ℝ) - 4 * n + 9) * (8 * (k : ℝ) - 4 * n + 11) * (8 * (k : ℝ) - 4 * n + 13)
        * (8 * (k : ℝ) - 4 * n + 15)))

/-- `a(n+1,k)/a(n,k)` (for `k < n`) -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  (4 * (n : ℝ) + 1) * (4 * (n : ℝ) + 3) ^ 2 * (8 * (k : ℝ) + 4 * n + 1) * (8 * (k : ℝ) + 4 * n + 3)
      * (8 * (k : ℝ) - 4 * n + 5) * (8 * (k : ℝ) - 4 * n + 7) /
    (131072 * (n : ℝ) * (2 * n - 1) * (4 * n - 1) * ((n : ℝ) - k) ^ 2 * ((n : ℝ) + k + 1) ^ 2)

/-- `A n = a(n,0)` for `n ≥ 1`: `A 1 = 9/256`, `A (n+1) = A n * ρ n 0`. -/
noncomputable def A (n : ℕ) : ℝ := 9 / 256 * ∏ m ∈ Ico 1 n, ρ m 0

/-- the first member of the WZ pair, `a(n,k)` (meaningful for `n ≥ 1`; zero for `k ≥ n ≥ 1`) -/
noncomputable def a (n k : ℕ) : ℝ := A n * ∏ j ∈ range k, σ n j

/-- the row sums `T_N = ∑_{k<N} a(N,k)` -/
noncomputable def T (N : ℕ) : ℝ := ∑ k ∈ range N, a N k

/-- Sun's summand -/
noncomputable def t (n : ℕ) : ℝ :=
  (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * (Nat.choose (2 * n) n : ℝ) ^ 2
      * (Nat.choose (4 * n) (2 * n) : ℝ) ^ 3 / ((2 * (n : ℝ) + 1) * 2 ^ (22 * n))

/-- the saddle point `x* = √3/2` -/
noncomputable def xstar : ℝ := Real.sqrt 3 / 2

/-- the closed form of `a(N,k)` in the region `k < N ≤ 2k+2` (natural-number subtractions:
`4k+4-2N`, `8k+8-4N`, `N-k-1`, `2N-2` are the true values there) -/
noncomputable def cf (N k : ℕ) : ℝ :=
  2 ^ 12 * (4 * (N : ℝ) - 1) * (2 * (k : ℝ) + 1) * ((4 * N).factorial : ℝ)
      * ((8 * k + 4 * N).factorial : ℝ) * ((4 * k + 4 - 2 * N).factorial : ℝ) /
    (2 ^ (22 * N) * ((2 * N).factorial : ℝ) * ((2 * N - 2).factorial : ℝ)
      * ((4 * k + 2 * N).factorial : ℝ) * ((8 * k + 8 - 4 * N).factorial : ℝ)
      * ((N - k - 1).factorial : ℝ) ^ 2 * ((N + k).factorial : ℝ) ^ 2)

lemma A_one : A 1 = 9 / 256 := by simp [A]

lemma A_succ {n : ℕ} (hn : 1 ≤ n) : A (n + 1) = A n * ρ n 0 := by
  unfold A
  rw [Finset.prod_Ico_succ_top hn, mul_assoc]

lemma a_zero_right (n : ℕ) : a n 0 = A n := by simp [a]

lemma a_succ_right (n k : ℕ) : a n (k + 1) = a n k * σ n k := by
  simp [a, prod_range_succ, mul_assoc]

end PiSeries.SunPi2
