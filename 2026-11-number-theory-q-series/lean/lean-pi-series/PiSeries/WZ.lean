/-
# The Wilf–Zeilberger telescoping argument, in the form used for Ramanujan-type series for `1/π`

All the series proved in this project have a WZ pair of the special shape

  `F(n,k) = n · h(n,k)`,   `G(n,k) = Q(n,k) · h(n,k)`,
  `(n+1) h(n+1,k) − n h(n,k) = Q(n,k+1) h(n,k+1) − Q(n,k) h(n,k)`        (WZ)

with `h ≥ 0`, `G(n,0)` the summand of the series, `F(0,k) = 0`, and the whole value of the series
carried by one boundary term: `lim_K G(0,K) = L`, while `lim_K G(n,K) = 0` for `n ≥ 1` and
`N · ∑_k h(N,k) → 0`.  `hasSum_of_wz` is that argument, once and for all.

No Carlson theorem and no analytic continuation in `k` is used: `k` is an integer throughout.
-/
import Mathlib

open Finset Filter Topology

namespace PiSeries

/-- **Finite WZ identity.**  If `u(n+1,k) − u(n,k) = v(n,k+1) − v(n,k)` for all `n, k`, then for all
`N, K`: `∑_{k<K} u(N,k) − ∑_{k<K} u(0,k) = ∑_{n<N} v(n,K) − ∑_{n<N} v(n,0)`. -/
theorem wz_finite_identity (u v : ℕ → ℕ → ℝ)
    (h : ∀ n k, u (n + 1) k - u n k = v n (k + 1) - v n k) (N K : ℕ) :
    ∑ k ∈ range K, u N k - ∑ k ∈ range K, u 0 k
      = ∑ n ∈ range N, v n K - ∑ n ∈ range N, v n 0 := by
  have h1 : ∑ n ∈ range N, ∑ k ∈ range K, (u (n + 1) k - u n k)
      = ∑ n ∈ range N, (v n K - v n 0) := by
    refine Finset.sum_congr rfl (fun n _ => ?_)
    simp_rw [h]
    exact Finset.sum_range_sub (fun k => v n k) K
  rw [Finset.sum_comm] at h1
  have h2 : ∀ k, ∑ n ∈ range N, (u (n + 1) k - u n k) = u N k - u 0 k :=
    fun k => Finset.sum_range_sub (fun n => u n k) N
  simp_rw [h2, Finset.sum_sub_distrib] at h1
  exact h1

/-- `∑_{k<K} 1/(k+1)² ≤ 2`. -/
lemma sum_inv_sq_le (K : ℕ) : ∑ k ∈ range K, 1 / ((k : ℝ) + 1) ^ 2 ≤ 2 := by
  have key : ∀ K : ℕ, ∑ k ∈ range (K + 1), 1 / ((k : ℝ) + 1) ^ 2 ≤ 2 - 1 / ((K : ℝ) + 1) := by
    intro K
    induction K with
    | zero => norm_num
    | succ K ih =>
      rw [sum_range_succ]
      have h1 : (0 : ℝ) < (K : ℝ) + 1 := by positivity
      have h2 : (0 : ℝ) < (K : ℝ) + 1 + 1 := by positivity
      have h3 : 1 / (((K + 1 : ℕ) : ℝ) + 1) ^ 2
          ≤ 1 / ((K : ℝ) + 1) - 1 / (((K + 1 : ℕ) : ℝ) + 1) := by
        push_cast
        rw [div_sub_div _ _ h1.ne' h2.ne', div_le_div_iff₀ (by positivity) (by positivity)]
        nlinarith
      linarith
  cases K with
  | zero => simp
  | succ K =>
    have h0 := key K
    have h1 : (0 : ℝ) < 1 / ((K : ℝ) + 1) := by positivity
    linarith

/-- **The WZ limit argument.**  Let `h, Q : ℕ → ℕ → ℝ` satisfy the WZ equation
`(n+1) h(n+1,k) − n h(n,k) = Q(n,k+1) h(n,k+1) − Q(n,k) h(n,k)`, with `h ≥ 0`,
`h(n,k) ≤ C(n)/(k+1)²` and `|Q(n,k)| ≤ A(n)(k+1)` for `n ≥ 1`, `n · C(n) → 0`, and
`Q(0,K) h(0,K) → L` as `K → ∞`.  If the terms `Q(n,0) h(n,0)` are non-negative, then
`∑_n Q(n,0) h(n,0) = L`. -/
theorem hasSum_of_wz (h Q : ℕ → ℕ → ℝ) (C A : ℕ → ℝ) (L : ℝ)
    (hwz : ∀ n k : ℕ, ((n : ℝ) + 1) * h (n + 1) k - n * h n k
      = Q n (k + 1) * h n (k + 1) - Q n k * h n k)
    (hpos : ∀ n k, 0 ≤ h n k)
    (hC : ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ, h n k ≤ C n / ((k : ℝ) + 1) ^ 2)
    (hA : ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ, |Q n k| ≤ A n * ((k : ℝ) + 1))
    (hCl : Tendsto (fun n : ℕ => (n : ℝ) * C n) atTop (𝓝 0))
    (hL : Tendsto (fun K => Q 0 K * h 0 K) atTop (𝓝 L))
    (hg : ∀ n, 0 ≤ Q n 0 * h n 0) :
    HasSum (fun n => Q n 0 * h n 0) L := by
  -- `C n ≥ 0` for `n ≥ 1`
  have hC0 : ∀ n, 1 ≤ n → 0 ≤ C n := by
    intro n hn
    have := (hpos n 0).trans (hC n hn 0)
    simpa using this
  -- partial sums in `k` are bounded
  have hpart : ∀ n, 1 ≤ n → ∀ K, ∑ k ∈ range K, h n k ≤ 2 * C n := by
    intro n hn K
    calc ∑ k ∈ range K, h n k ≤ ∑ k ∈ range K, C n * (1 / ((k : ℝ) + 1) ^ 2) :=
          Finset.sum_le_sum (fun k _ => by rw [mul_one_div]; exact hC n hn k)
      _ = C n * ∑ k ∈ range K, 1 / ((k : ℝ) + 1) ^ 2 := by rw [Finset.mul_sum]
      _ ≤ C n * 2 := mul_le_mul_of_nonneg_left (sum_inv_sq_le K) (hC0 n hn)
      _ = 2 * C n := mul_comm _ _
  have hsum : ∀ n, 1 ≤ n → Summable (h n) :=
    fun n hn => summable_of_sum_range_le (hpos n) (hpart n hn)
  have htsum : ∀ n, 1 ≤ n → ∑' k, h n k ≤ 2 * C n :=
    fun n hn => Real.tsum_le_of_sum_range_le (hpos n) (hpart n hn)
  -- `G(n,K) → 0` for `n ≥ 1`
  have hgK : ∀ n, 1 ≤ n → Tendsto (fun K => Q n K * h n K) atTop (𝓝 0) := by
    intro n hn
    have hA0 : 0 ≤ A n := by
      have := (abs_nonneg _).trans (hA n hn 0)
      simpa using this
    have hlim : Tendsto (fun K : ℕ => A n * C n * (1 / ((K : ℝ) + 1))) atTop (𝓝 0) := by
      have := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul (A n * C n)
      simpa using this
    refine squeeze_zero_norm (fun K => ?_) hlim
    rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (hpos n K)]
    have hK : (0 : ℝ) < (K : ℝ) + 1 := by positivity
    calc |Q n K| * h n K ≤ (A n * ((K : ℝ) + 1)) * (C n / ((K : ℝ) + 1) ^ 2) :=
          mul_le_mul (hA n hn K) (hC n hn K) (hpos n K) (mul_nonneg hA0 hK.le)
      _ = A n * C n * (1 / ((K : ℝ) + 1)) := by field_simp
  -- the finite identity, and `K → ∞`
  have hfin := wz_finite_identity (fun n k => (n : ℝ) * h n k) (fun n k => Q n k * h n k)
    (fun n k => by simpa using hwz n k)
  have hpartial : ∀ N, 1 ≤ N →
      ∑ n ∈ range N, Q n 0 * h n 0 = L - (N : ℝ) * ∑' k, h N k := by
    intro N hN
    have hL1 : Tendsto (fun K => ∑ k ∈ range K, (N : ℝ) * h N k
        - ∑ k ∈ range K, ((0 : ℕ) : ℝ) * h 0 k) atTop (𝓝 ((N : ℝ) * ∑' k, h N k)) := by
      have h1 : Tendsto (fun K => ∑ k ∈ range K, (N : ℝ) * h N k) atTop
          (𝓝 ((N : ℝ) * ∑' k, h N k)) := by
        simp_rw [← Finset.mul_sum]
        exact ((hsum N hN).hasSum.tendsto_sum_nat).const_mul _
      simpa using h1
    have hR1 : Tendsto (fun K => ∑ n ∈ range N, Q n K * h n K - ∑ n ∈ range N, Q n 0 * h n 0)
        atTop (𝓝 (∑ n ∈ range N, (if n = 0 then L else 0) - ∑ n ∈ range N, Q n 0 * h n 0)) := by
      refine (tendsto_finsetSum _ (fun n _ => ?_)).sub tendsto_const_nhds
      by_cases hn : n = 0
      · subst hn; simpa using hL
      · simpa [hn] using hgK n (Nat.one_le_iff_ne_zero.mpr hn)
    have heq : (fun K => ∑ k ∈ range K, (N : ℝ) * h N k - ∑ k ∈ range K, ((0 : ℕ) : ℝ) * h 0 k)
        = (fun K => ∑ n ∈ range N, Q n K * h n K - ∑ n ∈ range N, Q n 0 * h n 0) :=
      funext (fun K => hfin N K)
    rw [heq] at hL1
    have huniq := tendsto_nhds_unique hL1 hR1
    have hite : ∑ n ∈ range N, (if n = 0 then L else 0) = L := by
      rw [Finset.sum_ite_eq' (range N) 0 (fun _ => L)]
      simp [Nat.pos_of_ne_zero (Nat.one_le_iff_ne_zero.mp hN)]
    rw [hite] at huniq
    linarith
  -- `N ∑_k h(N,k) → 0`
  have hlim : Tendsto (fun N : ℕ => (N : ℝ) * ∑' k, h N k) atTop (𝓝 0) := by
    have h2 : Tendsto (fun N : ℕ => 2 * ((N : ℝ) * C N)) atTop (𝓝 0) := by
      simpa using hCl.const_mul 2
    refine squeeze_zero' (Eventually.of_forall fun N =>
      mul_nonneg (Nat.cast_nonneg N) (tsum_nonneg (hpos N))) ?_ h2
    filter_upwards [eventually_ge_atTop 1] with N hN
    calc (N : ℝ) * ∑' k, h N k ≤ (N : ℝ) * (2 * C N) :=
          mul_le_mul_of_nonneg_left (htsum N hN) (Nat.cast_nonneg N)
      _ = 2 * ((N : ℝ) * C N) := by ring
  rw [hasSum_iff_tendsto_nat_of_nonneg hg]
  have hfinal : Tendsto (fun N : ℕ => L - (N : ℝ) * ∑' k, h N k) atTop (𝓝 (L - 0)) :=
    tendsto_const_nhds.sub hlim
  rw [sub_zero] at hfinal
  refine hfinal.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  exact (hpartial N hN).symm

end PiSeries
