/-
# The WZ limit argument for pairs whose terms change sign

Same argument as `PiSeries.hasSum_of_wz`, without the hypothesis `h ≥ 0`: the bounds are on `|h|`,
and the conclusion is the convergence of the partial sums `∑_{n<N} Q(n,0) h(n,0)` (which is all one
can say when the series converges only conditionally, as Bauer's series does).
-/
import PiSeries.WZ
import PiSeries.Term

open Finset Filter Topology

namespace PiSeries

/-- **The WZ limit argument, signed version.**  Let `h, Q : ℕ → ℕ → ℝ` satisfy the WZ equation
`(n+1) h(n+1,k) − n h(n,k) = Q(n,k+1) h(n,k+1) − Q(n,k) h(n,k)`, with
`|h(n,k)| ≤ C(n)/(k+1)²` and `|Q(n,k)| ≤ A(n)(k+1)` for `n ≥ 1`, `n · C(n) → 0`, and
`Q(0,K) h(0,K) → L` as `K → ∞`.  Then `∑_{n<N} Q(n,0) h(n,0) → L` as `N → ∞`. -/
theorem tendsto_of_wz_signed (h Q : ℕ → ℕ → ℝ) (C A : ℕ → ℝ) (L : ℝ)
    (hwz : ∀ n k : ℕ, ((n : ℝ) + 1) * h (n + 1) k - n * h n k
      = Q n (k + 1) * h n (k + 1) - Q n k * h n k)
    (hC : ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ, |h n k| ≤ C n / ((k : ℝ) + 1) ^ 2)
    (hA : ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ, |Q n k| ≤ A n * ((k : ℝ) + 1))
    (hCl : Tendsto (fun n : ℕ => (n : ℝ) * C n) atTop (𝓝 0))
    (hL : Tendsto (fun K => Q 0 K * h 0 K) atTop (𝓝 L)) :
    Tendsto (fun N => ∑ n ∈ range N, Q n 0 * h n 0) atTop (𝓝 L) := by
  have hC0 : ∀ n, 1 ≤ n → 0 ≤ C n := by
    intro n hn
    have := (abs_nonneg _).trans (hC n hn 0)
    simpa using this
  have hpart : ∀ n, 1 ≤ n → ∀ K, ∑ k ∈ range K, |h n k| ≤ 2 * C n := by
    intro n hn K
    calc ∑ k ∈ range K, |h n k| ≤ ∑ k ∈ range K, C n * (1 / ((k : ℝ) + 1) ^ 2) :=
          Finset.sum_le_sum (fun k _ => by rw [mul_one_div]; exact hC n hn k)
      _ = C n * ∑ k ∈ range K, 1 / ((k : ℝ) + 1) ^ 2 := by rw [Finset.mul_sum]
      _ ≤ C n * 2 := mul_le_mul_of_nonneg_left (sum_inv_sq_le K) (hC0 n hn)
      _ = 2 * C n := mul_comm _ _
  have hsumabs : ∀ n, 1 ≤ n → Summable (fun k => |h n k|) :=
    fun n hn => summable_of_sum_range_le (fun k => abs_nonneg _) (hpart n hn)
  have hsum : ∀ n, 1 ≤ n → Summable (h n) := fun n hn => Summable.of_abs (hsumabs n hn)
  have htsum : ∀ n, 1 ≤ n → |∑' k, h n k| ≤ 2 * C n := by
    intro n hn
    have h1 : |∑' k, h n k| ≤ ∑' k, |h n k| := by
      have := norm_tsum_le_tsum_norm (f := h n) (hsum n hn).norm
      simpa only [Real.norm_eq_abs] using this
    exact h1.trans (Real.tsum_le_of_sum_range_le (fun k => abs_nonneg _) (hpart n hn))
  have hgK : ∀ n, 1 ≤ n → Tendsto (fun K => Q n K * h n K) atTop (𝓝 0) := by
    intro n hn
    have hA0 : 0 ≤ A n := by
      have := (abs_nonneg _).trans (hA n hn 0)
      simpa using this
    have hlim : Tendsto (fun K : ℕ => A n * C n * (1 / ((K : ℝ) + 1))) atTop (𝓝 0) := by
      have := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul (A n * C n)
      simpa using this
    refine squeeze_zero_norm (fun K => ?_) hlim
    rw [Real.norm_eq_abs, abs_mul]
    have hK : (0 : ℝ) < (K : ℝ) + 1 := by positivity
    calc |Q n K| * |h n K| ≤ (A n * ((K : ℝ) + 1)) * (C n / ((K : ℝ) + 1) ^ 2) :=
          mul_le_mul (hA n hn K) (hC n hn K) (abs_nonneg _) (mul_nonneg hA0 hK.le)
      _ = A n * C n * (1 / ((K : ℝ) + 1)) := by field_simp
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
  have hlim : Tendsto (fun N : ℕ => (N : ℝ) * ∑' k, h N k) atTop (𝓝 0) := by
    have h2 : Tendsto (fun N : ℕ => 2 * ((N : ℝ) * C N)) atTop (𝓝 0) := by
      simpa using hCl.const_mul 2
    refine squeeze_zero_norm' ?_ h2
    filter_upwards [eventually_ge_atTop 1] with N hN
    rw [Real.norm_eq_abs, abs_mul, abs_of_nonneg (Nat.cast_nonneg N)]
    calc (N : ℝ) * |∑' k, h N k| ≤ (N : ℝ) * (2 * C N) :=
          mul_le_mul_of_nonneg_left (htsum N hN) (Nat.cast_nonneg N)
      _ = 2 * ((N : ℝ) * C N) := by ring
  have hfinal : Tendsto (fun N : ℕ => L - (N : ℝ) * ∑' k, h N k) atTop (𝓝 (L - 0)) :=
    tendsto_const_nhds.sub hlim
  rw [sub_zero] at hfinal
  refine hfinal.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  exact (hpartial N hN).symm

/-- `hterm σ ρ n k = hterm σ ρ n 0 · ∏_{j<k} σ(n,j)`, the product written as an `hterm` with `ρ = 1`. -/
lemma hterm_eq_mul (σ ρ : ℕ → ℕ → ℝ) (n k : ℕ) :
    hterm σ ρ n k = hterm σ ρ n 0 * hterm σ (fun _ _ => (1 : ℝ)) n k := by
  simp [hterm]

/-- Signed version of `hterm_le_of_sigma_le` (with `a = 1`): if `0 < σ(n,k) ≤ ((k+1)/(k+2))²` for
all `k`, then `|h(n,k)| ≤ |h(n,0)|/(k+1)²`, whatever the sign of `ρ`. -/
lemma abs_hterm_le_of_sigma_le {σ ρ : ℕ → ℕ → ℝ} (n : ℕ) (hσ0 : ∀ k, 0 < σ n k)
    (hσ : ∀ k : ℕ, σ n k ≤ (((k : ℝ) + 1) / ((k : ℝ) + 1 + 1)) ^ 2) (k : ℕ) :
    |hterm σ ρ n k| ≤ |hterm σ ρ n 0| / ((k : ℝ) + 1) ^ 2 := by
  have h1 := hterm_le_of_sigma_le (σ := σ) (ρ := fun _ _ => (1 : ℝ)) (a := 1) le_rfl n hσ0
    (by simp [hterm]) hσ k
  have h0 : 0 ≤ hterm σ (fun _ _ => (1 : ℝ)) n k := by
    unfold hterm
    exact mul_nonneg (by simp) (prod_nonneg fun j _ => (hσ0 j).le)
  have h2 : hterm σ (fun _ _ => (1 : ℝ)) n k ≤ 1 / ((k : ℝ) + 1) ^ 2 := by
    simpa [hterm] using h1
  rw [hterm_eq_mul σ ρ n k, abs_mul, abs_of_nonneg h0, div_eq_mul_one_div]
  exact mul_le_mul_of_nonneg_left h2 (abs_nonneg _)

end PiSeries
