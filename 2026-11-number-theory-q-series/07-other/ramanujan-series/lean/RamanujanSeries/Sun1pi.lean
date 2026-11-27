/-
# A Ramanujan-type series for 1/π conjectured by Z.-W. Sun
(arXiv:2307.03086v10, Conjecture 5.25(i), eq. (5.20))

  ∑_{k ≥ 0} (360k³+612k²+230k+15) C(3k,k) C(6k,3k)² / ((3k+1)(3k+2) 2^(15k)) = 32 √2 / π.

Formal proof by the Wilf–Zeilberger method (Theorem 4 of RESULTS.md).  With

  h(n,k) = c(n) ∏_{j<k} σ(n,j),   u(n,k) = n h(n,k),   v(n,k) = Q(n,k) h(n,k)

(`σ`, `ρ`, `Q` the explicit rational functions below; `h(n,k+1) = σ(n,k) h(n,k)`,
`h(n+1,k) = ρ(n,k) h(n,k)`), the pair `(u, v)` is a WZ pair:
`u(n+1,k) - u(n,k) = v(n,k+1) - v(n,k)`.  It comes from Gauss's ₂F₁(1) summation theorem used as a
WZ seed (Au's `Seed1(1/2+n, 1+2n, k+n)`).  Everything is checked here:
the WZ equation, the telescoping, the limit `v(0,K) → 1/(Γ(1/4)Γ(3/4)) = √2/(2π)` (through
Euler's limit formula for `Γ` and the reflection formula), the vanishing of the other boundary
terms, and the identification of `64 v(n,0)` with Sun's summand.
-/
import Mathlib

open Finset Filter Topology

namespace SunSeries1pi

/-- `h(n,k+1)/h(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ :=
  ((4 * (k : ℝ) + 6 * n + 1) * (4 * (k : ℝ) + 6 * n + 3)) /
    (16 * ((k : ℝ) + n + 1) * ((k : ℝ) + 3 * n + 1))

/-- `h(n+1,k)/h(n,k)` -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  3 * (6 * (n : ℝ) + 1) * (6 * (n : ℝ) + 5) * (4 * (k : ℝ) + 6 * n + 1) * (4 * (k : ℝ) + 6 * n + 3) *
      (4 * (k : ℝ) + 6 * n + 5) /
    (1024 * ((n : ℝ) + 1) * ((k : ℝ) + n + 1) * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2) *
      ((k : ℝ) + 3 * n + 3))

/-- the WZ certificate -/
noncomputable def Q (n k : ℕ) : ℝ :=
  (64 * (k : ℝ) ^ 3 + (160 * (n : ℝ) + 144) * (k : ℝ) ^ 2 +
      (240 * (n : ℝ) ^ 2 + 432 * (n : ℝ) + 92) * (k : ℝ) +
      360 * (n : ℝ) ^ 3 + 612 * (n : ℝ) ^ 2 + 230 * (n : ℝ) + 15) /
    (64 * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2))

/-- `c(n) = h(n,0)` -/
noncomputable def c : ℕ → ℝ
  | 0 => 1
  | n + 1 => c n * ρ n 0

/-- the hypergeometric term -/
noncomputable def h (n : ℕ) : ℕ → ℝ
  | 0 => c n
  | k + 1 => h n k * σ n k

/-- first member of the WZ pair -/
noncomputable def u (n k : ℕ) : ℝ := (n : ℝ) * h n k

/-- second member of the WZ pair -/
noncomputable def v (n k : ℕ) : ℝ := Q n k * h n k

@[simp] lemma h_zero (n : ℕ) : h n 0 = c n := rfl
lemma h_succ (n k : ℕ) : h n (k + 1) = h n k * σ n k := rfl
lemma c_succ (n : ℕ) : c (n + 1) = c n * ρ n 0 := rfl

/-- compatibility of the two ratios -/
lemma compat (n k : ℕ) : ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1) := by
  unfold ρ σ
  push_cast
  have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h3 : ((k : ℝ) + 3 * n + 1) ≠ 0 := by positivity
  have h4 : ((k : ℝ) + 3 * n + 2) ≠ 0 := by positivity
  have h5 : ((k : ℝ) + 3 * n + 3) ≠ 0 := by positivity
  have h6 : ((k : ℝ) + (n + 1) + 1) ≠ 0 := by positivity
  have h7 : ((k : ℝ) + 3 * (n + 1) + 1) ≠ 0 := by positivity
  have h8 : ((k : ℝ) + 1 + n + 1) ≠ 0 := by positivity
  have h9 : ((k : ℝ) + 1 + 3 * n + 1) ≠ 0 := by positivity
  have h10 : ((k : ℝ) + 1 + 3 * n + 2) ≠ 0 := by positivity
  have h11 : ((k : ℝ) + 1 + 3 * n + 3) ≠ 0 := by positivity
  field_simp
  ring

lemma h_succ_left (n k : ℕ) : h (n + 1) k = h n k * ρ n k := by
  induction k with
  | zero => simp [c_succ]
  | succ k ih =>
    rw [h_succ, ih, h_succ, mul_assoc, compat, mul_assoc]

/-- the WZ equation, as an identity between the rational functions -/
lemma wz_rat (n k : ℕ) :
    ((n : ℝ) + 1) * ρ n k - n = Q n (k + 1) * σ n k - Q n k := by
  unfold ρ σ Q
  push_cast
  have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h3 : ((k : ℝ) + 3 * n + 1) ≠ 0 := by positivity
  have h4 : ((k : ℝ) + 3 * n + 2) ≠ 0 := by positivity
  have h5 : ((k : ℝ) + 3 * n + 3) ≠ 0 := by positivity
  have h9 : ((k : ℝ) + 1 + 3 * n + 1) ≠ 0 := by positivity
  have h10 : ((k : ℝ) + 1 + 3 * n + 2) ≠ 0 := by positivity
  field_simp
  ring

/-- **The WZ equation** `u(n+1,k) - u(n,k) = v(n,k+1) - v(n,k)`. -/
theorem wz (n k : ℕ) : u (n + 1) k - u n k = v n (k + 1) - v n k := by
  unfold u v
  rw [h_succ_left, h_succ]
  push_cast
  linear_combination (h n k) * wz_rat n k

/-- telescoping in `k` -/
lemma sum_k (n K : ℕ) : ∑ k ∈ range K, (u (n + 1) k - u n k) = v n K - v n 0 := by
  simp_rw [wz]
  exact Finset.sum_range_sub (fun k => v n k) K

/-- **The finite WZ identity**, valid for all `N`, `K`. -/
theorem finite_identity (N K : ℕ) :
    ∑ k ∈ range K, u N k - ∑ k ∈ range K, u 0 k
      = ∑ n ∈ range N, v n K - ∑ n ∈ range N, v n 0 := by
  have h1 : ∑ n ∈ range N, ∑ k ∈ range K, (u (n + 1) k - u n k)
      = ∑ n ∈ range N, (v n K - v n 0) := Finset.sum_congr rfl (fun n _ => sum_k n K)
  rw [Finset.sum_comm] at h1
  have h2 : ∀ k, ∑ n ∈ range N, (u (n + 1) k - u n k) = u N k - u 0 k :=
    fun k => Finset.sum_range_sub (fun n => u n k) N
  simp_rw [h2, Finset.sum_sub_distrib] at h1
  exact h1

lemma u_zero (k : ℕ) : u 0 k = 0 := by simp [u]

/-! ### Bounds -/

lemma σ_pos (n k : ℕ) : 0 < σ n k := by
  unfold σ; positivity

lemma ρ0_pos (n : ℕ) : 0 < ρ n 0 := by
  unfold ρ; push_cast; positivity

lemma c_pos (n : ℕ) : 0 < c n := by
  induction n with
  | zero => simp [c]
  | succ n ih => rw [c_succ]; exact mul_pos ih (ρ0_pos n)

lemma h_pos (n k : ℕ) : 0 < h n k := by
  induction k with
  | zero => simpa using c_pos n
  | succ k ih => rw [h_succ]; exact mul_pos ih (σ_pos n k)

/-- for `n ≥ 1` the ratio in `k` is at most `(k+3n)/(k+3n+2)` -/
lemma σ_le (n k : ℕ) :
    σ (n + 1) k ≤ ((k : ℝ) + 3 * n + 3) / ((k : ℝ) + 3 * n + 5) := by
  unfold σ
  push_cast
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have key : ((k : ℝ) + 3 * n + 3) * (16 * ((k : ℝ) + (n + 1) + 1) * ((k : ℝ) + 3 * (n + 1) + 1))
      - (4 * (k : ℝ) + 6 * (n + 1) + 1) * (4 * (k : ℝ) + 6 * (n + 1) + 3) * ((k : ℝ) + 3 * n + 5)
      = 16 * (k : ℝ) ^ 2 * n + 60 * k * (n : ℝ) ^ 2 + 112 * k * n + 33 * k + 36 * (n : ℝ) ^ 3
        + 156 * (n : ℝ) ^ 2 + 195 * n + 69 := by ring
  have hpos : (0 : ℝ) ≤ 16 * (k : ℝ) ^ 2 * n + 60 * k * (n : ℝ) ^ 2 + 112 * k * n + 33 * k
      + 36 * (n : ℝ) ^ 3 + 156 * (n : ℝ) ^ 2 + 195 * n + 69 := by positivity
  linarith

lemma ρ0_le (n : ℕ) : ρ n 0 ≤ 27 / 32 := by
  unfold ρ
  push_cast
  rw [div_le_iff₀ (by positivity)]
  have key : (27 / 32 : ℝ) * (1024 * ((n : ℝ) + 1) * ((0 : ℝ) + n + 1) * ((0 : ℝ) + 3 * n + 1) *
        ((0 : ℝ) + 3 * n + 2) * ((0 : ℝ) + 3 * n + 3))
      - 3 * (6 * (n : ℝ) + 1) * (6 * (n : ℝ) + 5) * (4 * (0 : ℝ) + 6 * n + 1) *
        (4 * (0 : ℝ) + 6 * n + 3) * (4 * (0 : ℝ) + 6 * n + 5)
      = 34992 * (n : ℝ) ^ 4 + 92016 * (n : ℝ) ^ 3 + 87480 * (n : ℝ) ^ 2 + 35190 * n + 4959 := by
    ring
  have hpos : (0 : ℝ) ≤ 34992 * (n : ℝ) ^ 4 + 92016 * (n : ℝ) ^ 3 + 87480 * (n : ℝ) ^ 2
      + 35190 * n + 4959 := by positivity
  linarith

lemma c_le (n : ℕ) : c n ≤ (27 / 32) ^ n := by
  induction n with
  | zero => simp [c]
  | succ n ih =>
    rw [c_succ, pow_succ]
    exact mul_le_mul ih (ρ0_le n) (ρ0_pos n).le (by positivity)

/-- decay in `k` of the rows `n ≥ 1` -/
lemma h_le (n k : ℕ) :
    h (n + 1) k ≤ c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4)) /
      (((k : ℝ) + 3 * n + 3) * ((k : ℝ) + 3 * n + 4)) := by
  induction k with
  | zero =>
    rw [h_zero]
    apply le_of_eq
    push_cast
    have h1 : (3 * (n : ℝ) + 3) ≠ 0 := by positivity
    have h2 : (3 * (n : ℝ) + 4) ≠ 0 := by positivity
    field_simp
    ring
  | succ k ih =>
    rw [h_succ]
    have hc := (c_pos (n + 1)).le
    calc h (n + 1) k * σ (n + 1) k
        ≤ (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4)) /
            (((k : ℝ) + 3 * n + 3) * ((k : ℝ) + 3 * n + 4))) *
          (((k : ℝ) + 3 * n + 3) / ((k : ℝ) + 3 * n + 5)) :=
          mul_le_mul ih (σ_le n k) (σ_pos _ _).le (by positivity)
      _ = c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4)) /
            ((((k + 1 : ℕ) : ℝ) + 3 * n + 3) * (((k + 1 : ℕ) : ℝ) + 3 * n + 4)) := by
          push_cast
          have h1 : ((k : ℝ) + 3 * n + 3) ≠ 0 := by positivity
          have h2 : ((k : ℝ) + 3 * n + 4) ≠ 0 := by positivity
          have h3 : ((k : ℝ) + 3 * n + 5) ≠ 0 := by positivity
          have h4 : ((k : ℝ) + 1 + 3 * n + 3) ≠ 0 := by positivity
          have h5 : ((k : ℝ) + 1 + 3 * n + 4) ≠ 0 := by positivity
          field_simp
          ring

lemma sum_h_le_aux (n K : ℕ) :
    ∑ k ∈ range K, h (n + 1) k ≤
      c (n + 1) * (3 * (n : ℝ) + 4) * (1 - (3 * (n : ℝ) + 3) / ((K : ℝ) + 3 * n + 3)) := by
  induction K with
  | zero =>
    have h1 : (3 * (n : ℝ) + 3) ≠ 0 := by positivity
    simp [h1]
  | succ K ih =>
    rw [Finset.sum_range_succ]
    have hK := h_le n K
    have heq : c (n + 1) * (3 * (n : ℝ) + 4) * (1 - (3 * (n : ℝ) + 3) / ((K : ℝ) + 3 * n + 3))
        + c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4)) /
            (((K : ℝ) + 3 * n + 3) * ((K : ℝ) + 3 * n + 4))
        = c (n + 1) * (3 * (n : ℝ) + 4) *
            (1 - (3 * (n : ℝ) + 3) / (((K + 1 : ℕ) : ℝ) + 3 * n + 3)) := by
      push_cast
      have h1 : ((K : ℝ) + 3 * n + 3) ≠ 0 := by positivity
      have h2 : ((K : ℝ) + 3 * n + 4) ≠ 0 := by positivity
      have h3 : ((K : ℝ) + 1 + 3 * n + 3) ≠ 0 := by positivity
      field_simp
      ring
    linarith

lemma sum_h_le (n K : ℕ) :
    ∑ k ∈ range K, h (n + 1) k ≤ c (n + 1) * (3 * (n : ℝ) + 4) := by
  refine (sum_h_le_aux n K).trans ?_
  have hc : 0 ≤ c (n + 1) * (3 * (n : ℝ) + 4) := mul_nonneg (c_pos _).le (by positivity)
  have hfrac : 0 ≤ (3 * (n : ℝ) + 3) / ((K : ℝ) + 3 * n + 3) := by positivity
  nlinarith

lemma summable_h (n : ℕ) : Summable (h (n + 1)) :=
  summable_of_sum_range_le (fun k => (h_pos _ k).le) (sum_h_le n)

lemma tsum_h_le (n : ℕ) : ∑' k, h (n + 1) k ≤ c (n + 1) * (3 * (n : ℝ) + 4) :=
  Real.tsum_le_of_sum_range_le (fun k => (h_pos _ k).le) (sum_h_le n)

lemma tsum_h_nonneg (n : ℕ) : 0 ≤ ∑' k, h (n + 1) k :=
  tsum_nonneg (fun k => (h_pos _ k).le)

lemma Q_nonneg (n k : ℕ) : 0 ≤ Q n k := by
  unfold Q; positivity

lemma Q_le (n k : ℕ) : Q n k ≤ (k : ℝ) + 3 * n + 3 := by
  unfold Q
  rw [div_le_iff₀ (by positivity)]
  have key : ((k : ℝ) + 3 * n + 3) * (64 * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2))
      - (64 * (k : ℝ) ^ 3 + (160 * (n : ℝ) + 144) * (k : ℝ) ^ 2 +
          (240 * (n : ℝ) ^ 2 + 432 * (n : ℝ) + 92) * (k : ℝ) +
          360 * (n : ℝ) ^ 3 + 612 * (n : ℝ) ^ 2 + 230 * (n : ℝ) + 15)
      = 416 * (k : ℝ) ^ 2 * n + 240 * (k : ℝ) ^ 2 + 1488 * k * (n : ℝ) ^ 2 + 1872 * k * n + 612 * k
        + 1368 * (n : ℝ) ^ 3 + 2844 * (n : ℝ) ^ 2 + 1882 * n + 369 := by ring
  have hpos : (0 : ℝ) ≤ 416 * (k : ℝ) ^ 2 * n + 240 * (k : ℝ) ^ 2 + 1488 * k * (n : ℝ) ^ 2
      + 1872 * k * n + 612 * k + 1368 * (n : ℝ) ^ 3 + 2844 * (n : ℝ) ^ 2 + 1882 * n + 369 := by
    positivity
  linarith

lemma v_nonneg (n k : ℕ) : 0 ≤ v n k := mul_nonneg (Q_nonneg n k) (h_pos n k).le

/-- the boundary terms of the rows `n ≥ 1` are `O(1/K)` -/
lemma v_le (n K : ℕ) :
    v (n + 1) K ≤ 2 * (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4))) * (1 / ((K : ℝ) + 1)) := by
  unfold v
  have hB : 0 ≤ c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4)) :=
    mul_nonneg (c_pos _).le (by positivity)
  have hQ : Q (n + 1) K ≤ (K : ℝ) + 3 * n + 6 := by
    have := Q_le (n + 1) K
    push_cast at this
    linarith
  have hfrac : ((K : ℝ) + 3 * n + 6) / (((K : ℝ) + 3 * n + 3) * ((K : ℝ) + 3 * n + 4))
      ≤ 2 * (1 / ((K : ℝ) + 1)) := by
    rw [mul_one_div, div_le_div_iff₀ (by positivity) (by positivity)]
    have key : 2 * (((K : ℝ) + 3 * n + 3) * ((K : ℝ) + 3 * n + 4))
        - ((K : ℝ) + 3 * n + 6) * ((K : ℝ) + 1)
        = (K : ℝ) ^ 2 + 9 * K * n + 7 * K + 18 * (n : ℝ) ^ 2 + 39 * n + 18 := by ring
    have hpos : (0 : ℝ) ≤ (K : ℝ) ^ 2 + 9 * K * n + 7 * K + 18 * (n : ℝ) ^ 2 + 39 * n + 18 := by
      positivity
    linarith
  calc Q (n + 1) K * h (n + 1) K
      ≤ ((K : ℝ) + 3 * n + 6) * (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4)) /
          (((K : ℝ) + 3 * n + 3) * ((K : ℝ) + 3 * n + 4))) :=
        mul_le_mul hQ (h_le n K) (h_pos _ _).le (by positivity)
    _ = (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4))) *
          (((K : ℝ) + 3 * n + 6) / (((K : ℝ) + 3 * n + 3) * ((K : ℝ) + 3 * n + 4))) := by ring
    _ ≤ (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4))) * (2 * (1 / ((K : ℝ) + 1))) :=
        mul_le_mul_of_nonneg_left hfrac hB
    _ = _ := by ring

lemma tendsto_v_succ (n : ℕ) : Tendsto (fun K => v (n + 1) K) atTop (𝓝 0) := by
  have hlim : Tendsto (fun K : ℕ =>
      2 * (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4))) * (1 / ((K : ℝ) + 1))) atTop
      (𝓝 (2 * (c (n + 1) * ((3 * (n : ℝ) + 3) * (3 * (n : ℝ) + 4))) * 0)) :=
    (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul _
  rw [mul_zero] at hlim
  exact squeeze_zero (fun K => v_nonneg _ K) (fun K => v_le n K) hlim

/-! ### The row `n = 0`: Euler's limit formula -/

lemma h0_eq (K : ℕ) :
    h 0 K = (∏ j ∈ range K, ((1 / 4 : ℝ) + j)) * (∏ j ∈ range K, ((3 / 4 : ℝ) + j)) /
      ((K.factorial : ℝ)) ^ 2 := by
  induction K with
  | zero => simp [c]
  | succ K ih =>
    rw [h_succ, ih, Finset.prod_range_succ, Finset.prod_range_succ, Nat.factorial_succ]
    unfold σ
    push_cast
    have h1 : ((K.factorial : ℝ)) ≠ 0 := by positivity
    have h2 : ((K : ℝ) + 1) ≠ 0 := by positivity
    field_simp
    ring

lemma prod_pos_aux (s : ℝ) (hs : 0 < s) (K : ℕ) : 0 < ∏ j ∈ range K, (s + (j : ℝ)) :=
  Finset.prod_pos (fun j _ => by positivity)

/-- `v(0,K)` in terms of Euler's sequences for `Γ(1/4)` and `Γ(3/4)` -/
lemma v0_mul (K : ℕ) :
    v 0 K * (Real.GammaSeq (1 / 4) K * Real.GammaSeq (3 / 4) K)
      = (K : ℝ) * (4 * (K : ℝ) + 5) / (4 * ((K : ℝ) + 1) * ((K : ℝ) + 2)) := by
  unfold v Q Real.GammaSeq
  rw [h0_eq, Finset.prod_range_succ, Finset.prod_range_succ]
  have hP1 := (prod_pos_aux (1 / 4) (by norm_num) K).ne'
  have hP3 := (prod_pos_aux (3 / 4) (by norm_num) K).ne'
  have hf : ((K.factorial : ℝ)) ≠ 0 := by positivity
  have hK : (0 : ℝ) ≤ (K : ℝ) := Nat.cast_nonneg K
  have hpow : (K : ℝ) ^ (1 / 4 : ℝ) * (K : ℝ) ^ (3 / 4 : ℝ) = K := by
    rw [← Real.rpow_add_of_nonneg hK (by norm_num) (by norm_num)]
    norm_num
  have h1 : ((1 / 4 : ℝ) + K) ≠ 0 := by positivity
  have h2 : ((3 / 4 : ℝ) + K) ≠ 0 := by positivity
  have h3 : ((K : ℝ) + 1) ≠ 0 := by positivity
  have h4 : ((K : ℝ) + 2) ≠ 0 := by positivity
  have h5 : ((K : ℝ) + 3 * (0 : ℕ) + 1) ≠ 0 := by positivity
  have h6 : ((K : ℝ) + 3 * (0 : ℕ) + 2) ≠ 0 := by positivity
  push_cast
  have e : ∀ A B : ℝ,
      (64 * (K : ℝ) ^ 3 + (160 * (0 : ℝ) + 144) * (K : ℝ) ^ 2 +
            (240 * (0 : ℝ) ^ 2 + 432 * (0 : ℝ) + 92) * (K : ℝ) +
            360 * (0 : ℝ) ^ 3 + 612 * (0 : ℝ) ^ 2 + 230 * (0 : ℝ) + 15) /
          (64 * ((K : ℝ) + 3 * 0 + 1) * ((K : ℝ) + 3 * 0 + 2)) *
        ((∏ j ∈ range K, ((1 / 4 : ℝ) + j)) * (∏ j ∈ range K, ((3 / 4 : ℝ) + j)) /
          ((K.factorial : ℝ)) ^ 2) *
        (A * (K.factorial : ℝ) / ((∏ j ∈ range K, ((1 / 4 : ℝ) + j)) * ((1 / 4 : ℝ) + K)) *
          (B * (K.factorial : ℝ) / ((∏ j ∈ range K, ((3 / 4 : ℝ) + j)) * ((3 / 4 : ℝ) + K))))
      = (A * B) * (4 * (K : ℝ) + 5) / (4 * ((K : ℝ) + 1) * ((K : ℝ) + 2)) := by
    intro A B
    field_simp
    ring
  rw [e, hpow]

/-- the limit of the boundary term of the row `n = 0` -/
lemma tendsto_v_zero : Tendsto (fun K => v 0 K) atTop (𝓝 (Real.sqrt 2 / (2 * Real.pi))) := by
  -- the Gamma product
  have hG : Tendsto (fun K => Real.GammaSeq (1 / 4) K * Real.GammaSeq (3 / 4) K) atTop
      (𝓝 (Real.Gamma (1 / 4) * Real.Gamma (3 / 4))) :=
    (Real.GammaSeq_tendsto_Gamma _).mul (Real.GammaSeq_tendsto_Gamma _)
  have hval : Real.Gamma (1 / 4) * Real.Gamma (3 / 4) = Real.pi / (Real.sqrt 2 / 2) := by
    have := Real.Gamma_mul_Gamma_one_sub (1 / 4)
    rw [show (1 : ℝ) - 1 / 4 = 3 / 4 by norm_num,
      show Real.pi * (1 / 4) = Real.pi / 4 by ring, Real.sin_pi_div_four] at this
    exact this
  have hs2 : (0 : ℝ) < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  have hγ : Real.pi / (Real.sqrt 2 / 2) ≠ 0 := by
    have := Real.pi_pos
    positivity
  rw [hval] at hG
  -- the rational factor
  have hr : Tendsto (fun K : ℕ => (K : ℝ) * (4 * (K : ℝ) + 5) / (4 * ((K : ℝ) + 1) * ((K : ℝ) + 2)))
      atTop (𝓝 1) := by
    have he : Tendsto (fun K : ℕ => (7 * (K : ℝ) + 8) / (4 * ((K : ℝ) + 1) * ((K : ℝ) + 2))) atTop
        (𝓝 0) := by
      have hlim : Tendsto (fun K : ℕ => 2 * (1 / ((K : ℝ) + 1))) atTop (𝓝 (2 * 0)) :=
        (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul _
      rw [mul_zero] at hlim
      refine squeeze_zero (fun K => by positivity) (fun K => ?_) hlim
      rw [mul_one_div, div_le_div_iff₀ (by positivity) (by positivity)]
      have hK : (0 : ℝ) ≤ (K : ℝ) := Nat.cast_nonneg K
      nlinarith [mul_nonneg hK hK]
    have h1 : Tendsto (fun K : ℕ => 1 - (7 * (K : ℝ) + 8) / (4 * ((K : ℝ) + 1) * ((K : ℝ) + 2)))
        atTop (𝓝 (1 - 0)) := tendsto_const_nhds.sub he
    rw [sub_zero] at h1
    refine h1.congr (fun K => ?_)
    have h3 : ((K : ℝ) + 1) ≠ 0 := by positivity
    have h4 : ((K : ℝ) + 2) ≠ 0 := by positivity
    field_simp
    ring
  have hdiv := hr.div hG hγ
  have hfinal : (1 : ℝ) / (Real.pi / (Real.sqrt 2 / 2)) = Real.sqrt 2 / (2 * Real.pi) := by
    have := Real.pi_pos
    field_simp
  rw [hfinal] at hdiv
  refine hdiv.congr' ?_
  have hne : ∀ᶠ K : ℕ in atTop,
      Real.GammaSeq (1 / 4) K * Real.GammaSeq (3 / 4) K ≠ 0 := hG.eventually_ne hγ
  filter_upwards [hne] with K hK
  rw [Pi.div_apply, ← v0_mul K, mul_div_assoc, div_self hK, mul_one]

/-! ### The limit argument -/

/-- letting `K → ∞` in the finite identity -/
lemma partial_sum_v (n : ℕ) :
    ∑ j ∈ range (n + 1), v j 0
      = Real.sqrt 2 / (2 * Real.pi) - ((n : ℝ) + 1) * ∑' k, h (n + 1) k := by
  have hL : Tendsto (fun K => ∑ k ∈ range K, u (n + 1) k - ∑ k ∈ range K, u 0 k) atTop
      (𝓝 (((n : ℝ) + 1) * ∑' k, h (n + 1) k)) := by
    have h1 : Tendsto (fun K => ∑ k ∈ range K, h (n + 1) k) atTop (𝓝 (∑' k, h (n + 1) k)) :=
      (summable_h n).hasSum.tendsto_sum_nat
    have h2 := h1.const_mul ((n : ℝ) + 1)
    refine h2.congr (fun K => ?_)
    simp [u, Finset.mul_sum]
  have hR : Tendsto (fun K => ∑ j ∈ range (n + 1), v j K - ∑ j ∈ range (n + 1), v j 0) atTop
      (𝓝 ((∑ j ∈ range n, (0 : ℝ) + Real.sqrt 2 / (2 * Real.pi)) - ∑ j ∈ range (n + 1), v j 0)) := by
    refine Tendsto.sub ?_ tendsto_const_nhds
    have h1 : Tendsto (fun K => ∑ j ∈ range n, v (j + 1) K + v 0 K) atTop
        (𝓝 (∑ j ∈ range n, (0 : ℝ) + Real.sqrt 2 / (2 * Real.pi))) :=
      (tendsto_finsetSum _ (fun j _ => tendsto_v_succ j)).add tendsto_v_zero
    refine h1.congr (fun K => ?_)
    rw [Finset.sum_range_succ']
  have heq : (fun K => ∑ k ∈ range K, u (n + 1) k - ∑ k ∈ range K, u 0 k)
      = (fun K => ∑ j ∈ range (n + 1), v j K - ∑ j ∈ range (n + 1), v j 0) :=
    funext (fun K => finite_identity (n + 1) K)
  rw [heq] at hL
  have := tendsto_nhds_unique hL hR
  simp only [Finset.sum_const_zero, zero_add] at this
  linarith

/-- `∑_n v(n,0) = √2/(2π)` -/
theorem hasSum_v0 : HasSum (fun n => v n 0) (Real.sqrt 2 / (2 * Real.pi)) := by
  rw [hasSum_iff_tendsto_nat_of_nonneg (fun n => v_nonneg n 0)]
  rw [← tendsto_add_atTop_iff_nat 1]
  -- the remainder `(n+1) ∑_k h(n+1,k)` tends to `0`
  have hrem : Tendsto (fun n : ℕ => ((n : ℝ) + 1) * ∑' k, h (n + 1) k) atTop (𝓝 0) := by
    have hgeom : Tendsto (fun m : ℕ => 4 * ((m : ℝ) ^ 2 * (27 / 32 : ℝ) ^ m)) atTop (𝓝 (4 * 0)) :=
      (tendsto_pow_const_mul_const_pow_of_abs_lt_one 2
        (by rw [abs_of_pos] <;> norm_num : |(27 / 32 : ℝ)| < 1)).const_mul 4
    rw [mul_zero] at hgeom
    have hgeom' : Tendsto (fun n : ℕ => 4 * ((((n + 1 : ℕ) : ℝ)) ^ 2 * (27 / 32 : ℝ) ^ (n + 1)))
        atTop (𝓝 0) := (tendsto_add_atTop_iff_nat 1).mpr hgeom
    refine squeeze_zero
      (fun n => mul_nonneg (by positivity) (tsum_h_nonneg n)) (fun n => ?_) hgeom'
    have h1 := tsum_h_le n
    have h2 := c_le (n + 1)
    have hn : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
    calc ((n : ℝ) + 1) * ∑' k, h (n + 1) k
        ≤ ((n : ℝ) + 1) * (c (n + 1) * (3 * (n : ℝ) + 4)) :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = (((n : ℝ) + 1) * (3 * (n : ℝ) + 4)) * c (n + 1) := by ring
      _ ≤ (((n : ℝ) + 1) * (3 * (n : ℝ) + 4)) * (27 / 32) ^ (n + 1) :=
          mul_le_mul_of_nonneg_left h2 (by positivity)
      _ ≤ (4 * ((n : ℝ) + 1) ^ 2) * (27 / 32) ^ (n + 1) := by
          apply mul_le_mul_of_nonneg_right _ (by positivity)
          nlinarith
      _ = 4 * ((((n + 1 : ℕ) : ℝ)) ^ 2 * (27 / 32 : ℝ) ^ (n + 1)) := by
          push_cast; ring
  have : Tendsto (fun n : ℕ => Real.sqrt 2 / (2 * Real.pi) - ((n : ℝ) + 1) * ∑' k, h (n + 1) k)
      atTop (𝓝 (Real.sqrt 2 / (2 * Real.pi) - 0)) := tendsto_const_nhds.sub hrem
  rw [sub_zero] at this
  refine this.congr (fun n => ?_)
  exact (partial_sum_v n).symm

/-! ### Identification of `64 v(n,0)` with Sun's summand -/

lemma centralBinom_step (n : ℕ) :
    ((Nat.centralBinom (3 * n + 3) : ℕ) : ℝ) * ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3))
      = 8 * ((6 * (n : ℝ) + 1) * (6 * n + 3) * (6 * n + 5)) * (Nat.centralBinom (3 * n) : ℝ) := by
  have h1 := Nat.succ_mul_centralBinom_succ (3 * n)
  have h2 := Nat.succ_mul_centralBinom_succ (3 * n + 1)
  have h3 := Nat.succ_mul_centralBinom_succ (3 * n + 2)
  have e1 : ((3 * (n : ℝ) + 1)) * (Nat.centralBinom (3 * n + 1) : ℝ)
      = 2 * (2 * (3 * (n : ℝ)) + 1) * (Nat.centralBinom (3 * n) : ℝ) := by exact_mod_cast h1
  have e2 : ((3 * (n : ℝ) + 1 + 1)) * (Nat.centralBinom (3 * n + 1 + 1) : ℝ)
      = 2 * (2 * (3 * (n : ℝ) + 1) + 1) * (Nat.centralBinom (3 * n + 1) : ℝ) := by exact_mod_cast h2
  have e3 : ((3 * (n : ℝ) + 2 + 1)) * (Nat.centralBinom (3 * n + 2 + 1) : ℝ)
      = 2 * (2 * (3 * (n : ℝ) + 2) + 1) * (Nat.centralBinom (3 * n + 2) : ℝ) := by exact_mod_cast h3
  have e2' : ((3 * (n : ℝ) + 2)) * (Nat.centralBinom (3 * n + 2) : ℝ)
      = 2 * (6 * (n : ℝ) + 3) * (Nat.centralBinom (3 * n + 1) : ℝ) := by
    rw [show (3 * n + 1 + 1 : ℕ) = 3 * n + 2 from rfl] at e2
    linarith
  have e3' : ((3 * (n : ℝ) + 3)) * (Nat.centralBinom (3 * n + 3) : ℝ)
      = 2 * (6 * (n : ℝ) + 5) * (Nat.centralBinom (3 * n + 2) : ℝ) := by
    rw [show (3 * n + 2 + 1 : ℕ) = 3 * n + 3 from rfl] at e3
    linarith
  calc ((Nat.centralBinom (3 * n + 3) : ℕ) : ℝ) * ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3))
      = ((3 * (n : ℝ) + 3) * (Nat.centralBinom (3 * n + 3) : ℝ)) * ((3 * (n : ℝ) + 1) * (3 * n + 2)) := by ring
    _ = 2 * (6 * (n : ℝ) + 5) * ((3 * (n : ℝ) + 2) * (Nat.centralBinom (3 * n + 2) : ℝ)) * (3 * (n : ℝ) + 1) := by
        rw [e3']; ring
    _ = 2 * (6 * (n : ℝ) + 5) * (2 * (6 * (n : ℝ) + 3)) * ((3 * (n : ℝ) + 1) * (Nat.centralBinom (3 * n + 1) : ℝ)) := by
        rw [e2']; ring
    _ = _ := by rw [e1]; ring

lemma choose_step (n : ℕ) :
    ((Nat.choose (3 * n + 3) (n + 1) : ℕ) : ℝ) * (((n : ℝ) + 1) * (2 * n + 1) * (2 * n + 2))
      = ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3)) * (Nat.choose (3 * n) n : ℝ) := by
  have hA := Nat.add_one_mul_choose_eq (3 * n + 2) n
  have hB := Nat.choose_mul_succ_eq (3 * n + 1) n
  have hC := Nat.choose_mul_succ_eq (3 * n) n
  rw [show 3 * n + 1 + 1 - n = 2 * n + 2 by omega] at hB
  rw [show 3 * n + 1 - n = 2 * n + 1 by omega] at hC
  have eA : (3 * (n : ℝ) + 2 + 1) * (Nat.choose (3 * n + 2) n : ℝ)
      = (Nat.choose (3 * n + 2 + 1) (n + 1) : ℝ) * ((n : ℝ) + 1) := by exact_mod_cast hA
  have eB : (Nat.choose (3 * n + 1) n : ℝ) * (3 * (n : ℝ) + 1 + 1)
      = (Nat.choose (3 * n + 1 + 1) n : ℝ) * (2 * (n : ℝ) + 2) := by exact_mod_cast hB
  have eC : (Nat.choose (3 * n) n : ℝ) * (3 * (n : ℝ) + 1)
      = (Nat.choose (3 * n + 1) n : ℝ) * (2 * (n : ℝ) + 1) := by exact_mod_cast hC
  rw [show (3 * n + 2 + 1 : ℕ) = 3 * n + 3 from rfl] at eA
  rw [show (3 * n + 1 + 1 : ℕ) = 3 * n + 2 from rfl] at eB
  linear_combination (-(2 * (n : ℝ) + 2) * (2 * n + 1)) * eA
    - ((3 * (n : ℝ) + 3) * (2 * n + 1)) * eB - ((3 * (n : ℝ) + 3) * (3 * n + 2)) * eC

lemma c_eq (n : ℕ) :
    c n = (Nat.choose (3 * n) n : ℝ) * (Nat.centralBinom (3 * n) : ℝ) ^ 2 / 2 ^ (15 * n) := by
  induction n with
  | zero => simp [c]
  | succ n ih =>
    rw [c_succ, ih]
    unfold ρ
    have hcb := centralBinom_step n
    have hch := choose_step n
    rw [show 3 * (n + 1) = 3 * n + 3 by ring, show 15 * (n + 1) = 15 * n + 15 by ring, pow_add]
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
    have h3 : (2 * (n : ℝ) + 2) ≠ 0 := by positivity
    have h4 : (3 * (n : ℝ) + 1) ≠ 0 := by positivity
    have h5 : (3 * (n : ℝ) + 2) ≠ 0 := by positivity
    have h6 : (3 * (n : ℝ) + 3) ≠ 0 := by positivity
    have eX : ((Nat.choose (3 * n + 3) (n + 1) : ℕ) : ℝ)
        = ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3)) * (Nat.choose (3 * n) n : ℝ) /
          (((n : ℝ) + 1) * (2 * n + 1) * (2 * n + 2)) := by
      rw [eq_div_iff (by positivity)]
      exact hch
    have eV : ((Nat.centralBinom (3 * n + 3) : ℕ) : ℝ)
        = 8 * ((6 * (n : ℝ) + 1) * (6 * n + 3) * (6 * n + 5)) * (Nat.centralBinom (3 * n) : ℝ) /
          ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3)) := by
      rw [eq_div_iff (by positivity)]
      exact hcb
    rw [eX, eV]
    push_cast
    field_simp
    ring

/-- Sun's summand -/
noncomputable def sunTerm (k : ℕ) : ℝ :=
  (360 * (k : ℝ) ^ 3 + 612 * (k : ℝ) ^ 2 + 230 * (k : ℝ) + 15) * (Nat.choose (3 * k) k : ℝ) *
      (Nat.centralBinom (3 * k) : ℝ) ^ 2 /
    ((3 * (k : ℝ) + 1) * (3 * (k : ℝ) + 2) * 2 ^ (15 * k))

lemma v0_eq (n : ℕ) : 64 * v n 0 = sunTerm n := by
  unfold v sunTerm Q
  rw [h_zero, c_eq]
  push_cast
  have h1 : (3 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : (3 * (n : ℝ) + 2) ≠ 0 := by positivity
  have h3 : ((0 : ℝ) + 3 * n + 1) ≠ 0 := by positivity
  have h4 : ((0 : ℝ) + 3 * n + 2) ≠ 0 := by positivity
  field_simp
  ring

/-- **Sun's conjecture** (arXiv:2307.03086v10, Conj. 5.25(i)):
`∑_{k≥0} (360k³+612k²+230k+15) C(3k,k) C(6k,3k)² / ((3k+1)(3k+2) 2^(15k)) = 32√2/π`. -/
theorem sun_series_1pi :
    HasSum (fun k : ℕ =>
      (360 * (k : ℝ) ^ 3 + 612 * (k : ℝ) ^ 2 + 230 * (k : ℝ) + 15) * (Nat.choose (3 * k) k : ℝ) *
          (Nat.centralBinom (3 * k) : ℝ) ^ 2 /
        ((3 * (k : ℝ) + 1) * (3 * (k : ℝ) + 2) * 2 ^ (15 * k)))
      (32 * Real.sqrt 2 / Real.pi) := by
  have h := hasSum_v0.mul_left 64
  simp_rw [v0_eq] at h
  have hval : 64 * (Real.sqrt 2 / (2 * Real.pi)) = 32 * Real.sqrt 2 / Real.pi := by
    have := Real.pi_pos
    field_simp
    ring
  rw [hval] at h
  exact h

end SunSeries1pi
