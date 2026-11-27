/-
# WZ pairs supported on a triangle `0 ≤ k < n`

Several of the new series have a WZ pair `(a, R·a)` that lives on the triangle `k < n` (`n ≥ 1`):
`a(n,k+1) = σ(n,k) a(n,k)`, `a(n+1,k) = ρ(n,k) a(n,k)` (`k < n`), the WZ equation
`ρ − 1 = R(n,k+1) σ − R(n,k)` for `k + 1 < n`, and a second rational identity on the boundary row
`k = n − 1`.  Then the row sums `T n = ∑_{k<n} a(n,k)` satisfy `T(n+1) − T(n) = −R(n,0) a(n,0)`
(`tri_rowsum_succ`), i.e. `∑_{1≤n<N} (−R(n,0) a(n,0)) = T N − T 1` (`tri_sum`): the series
`∑ −R(n,0) a(n,0)` is reduced to the limit of the row sums.

Used for Sun's `64√2/π²` series (there with the copy in `SunPi2/PartA.lean`) and for the
zero-sum series of `SunZero.lean`.
-/
import Mathlib

open Finset

namespace PiSeries

/-- the term on the triangle with `a(1,0) = A1`, `a(n+1,0) = ρ(n,0) a(n,0)` (`n ≥ 1`),
`a(n,k+1) = σ(n,k) a(n,k)` -/
noncomputable def tri (A1 : ℝ) (σ ρ : ℕ → ℕ → ℝ) (n k : ℕ) : ℝ :=
  (A1 * ∏ m ∈ Ico 1 n, ρ m 0) * ∏ j ∈ range k, σ n j

section
variable {A1 : ℝ} {σ ρ R : ℕ → ℕ → ℝ}

lemma tri_one_zero : tri A1 σ ρ 1 0 = A1 := by simp [tri]

lemma tri_succ_right (n k : ℕ) : tri A1 σ ρ n (k + 1) = tri A1 σ ρ n k * σ n k := by
  simp [tri, prod_range_succ, mul_assoc]

lemma tri_succ_zero {n : ℕ} (hn : 1 ≤ n) : tri A1 σ ρ (n + 1) 0 = tri A1 σ ρ n 0 * ρ n 0 := by
  simp [tri, Finset.prod_Ico_succ_top hn, mul_assoc]

/-- `a(n+1,k) = ρ(n,k) a(n,k)` on the triangle -/
lemma tri_succ_left
    (compat : ∀ n k : ℕ, k + 1 < n → ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1))
    {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) :
    tri A1 σ ρ (n + 1) k = tri A1 σ ρ n k * ρ n k := by
  induction k with
  | zero => exact tri_succ_zero hn
  | succ k ih =>
    rw [tri_succ_right, ih (by omega), tri_succ_right, mul_assoc, compat n k hk, mul_assoc]

/-- **Row sums of a triangle WZ pair**: `T(n+1) = T(n) − R(n,0) a(n,0)`. -/
theorem tri_rowsum_succ
    (compat : ∀ n k : ℕ, k + 1 < n → ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1))
    (wz : ∀ n k : ℕ, k + 1 < n → ρ n k - 1 = R n (k + 1) * σ n k - R n k)
    (bdry : ∀ n : ℕ, 1 ≤ n →
      ρ n (n - 1) - 1 = -(ρ n (n - 1) * σ (n + 1) (n - 1)) - R n (n - 1))
    {n : ℕ} (hn : 1 ≤ n) :
    ∑ k ∈ range (n + 1), tri A1 σ ρ (n + 1) k
      = ∑ k ∈ range n, tri A1 σ ρ n k - R n 0 * tri A1 σ ρ n 0 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  set a := tri A1 σ ρ with ha
  have hleft : ∀ k, k < m + 1 → a (m + 1 + 1) k = a (m + 1) k * ρ (m + 1) k :=
    fun k hk => tri_succ_left compat (by omega) hk
  -- interior columns telescope
  have hint : ∑ k ∈ range m, (a (m + 1 + 1) k - a (m + 1) k)
      = R (m + 1) m * a (m + 1) m - R (m + 1) 0 * a (m + 1) 0 := by
    have h1 : ∀ k ∈ range m, a (m + 1 + 1) k - a (m + 1) k
        = R (m + 1) (k + 1) * a (m + 1) (k + 1) - R (m + 1) k * a (m + 1) k := by
      intro k hk
      have hk' := Finset.mem_range.mp hk
      rw [hleft k (by omega), ha, tri_succ_right]
      linear_combination (tri A1 σ ρ (m + 1) k) * wz (m + 1) k (by omega)
    rw [Finset.sum_congr rfl h1]
    exact Finset.sum_range_sub (fun k => R (m + 1) k * a (m + 1) k) m
  have hb := bdry (m + 1) (by omega)
  simp only [Nat.add_sub_cancel] at hb
  have hlast : a (m + 1 + 1) (m + 1) = a (m + 1) m * ρ (m + 1) m * σ (m + 1 + 1) m := by
    rw [ha, tri_succ_right, ← ha, hleft m (by omega)]
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_succ (fun k => a (m + 1) k),
    hlast, hleft m (by omega)]
  have hsplit : ∑ k ∈ range m, a (m + 1 + 1) k
      = ∑ k ∈ range m, a (m + 1) k + (R (m + 1) m * a (m + 1) m - R (m + 1) 0 * a (m + 1) 0) := by
    rw [← hint, Finset.sum_sub_distrib]; ring
  rw [hsplit]
  linear_combination (a (m + 1) m) * hb

/-- **The finite identity of a triangle WZ pair**:
`∑_{1≤n<N} (−R(n,0) a(n,0)) = ∑_{k<N} a(N,k) − a(1,0)`. -/
theorem tri_sum
    (compat : ∀ n k : ℕ, k + 1 < n → ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1))
    (wz : ∀ n k : ℕ, k + 1 < n → ρ n k - 1 = R n (k + 1) * σ n k - R n k)
    (bdry : ∀ n : ℕ, 1 ≤ n →
      ρ n (n - 1) - 1 = -(ρ n (n - 1) * σ (n + 1) (n - 1)) - R n (n - 1))
    {N : ℕ} (hN : 1 ≤ N) :
    ∑ n ∈ Ico 1 N, (-(R n 0 * tri A1 σ ρ n 0)) = ∑ k ∈ range N, tri A1 σ ρ N k - A1 := by
  induction N, hN using Nat.le_induction with
  | base => simp [tri_one_zero]
  | succ N hN ih =>
    rw [Finset.sum_Ico_succ_top hN, ih, tri_rowsum_succ compat wz bdry hN]
    ring

end

end PiSeries
