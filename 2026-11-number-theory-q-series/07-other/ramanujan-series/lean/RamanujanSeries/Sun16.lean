/-
# A conjecture of Z.-W. Sun (arXiv:2506.01870v6, Conjecture 3.2, first identity)

  ∑_{k ≥ 0} (3k+2)² 16^k / ((2k+1)(6k+1)(6k+5) C(6k,3k)) = π / (2 √3).

Formal proof by the Wilf–Zeilberger method.  The WZ pair is

  f(n,k) = c(n) ∏_{j<k} σ(n,j),   g(n,k) = R(n,k) f(n,k),

with the explicit rational functions `σ`, `ρ`, `R` below (so that `f(n,k+1) = σ(n,k) f(n,k)`,
`f(n+1,k) = ρ(n,k) f(n,k)`), coming from the strange evaluation
`₂F₁(2a, 1-a; 2a+1/2; -1/3)` (Ebisu) used as a WZ seed with `a ↦ 1/2 + n`, `k ↦ k + n`.
Everything (the WZ equation, the telescoping, both limits and the evaluation of
`∑ f(0,k) = ∑ (-1/3)^k/(2k+1) = π/(2√3)`) is checked here.
-/
import Mathlib

open Finset Filter Topology

namespace SunSeries16

/-- `f(n,k+1)/f(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ :=
  -(((k : ℝ) + 3 * n + 1) * (2 * k + 1)) / (3 * (2 * (k : ℝ) + 6 * n + 3) * ((k : ℝ) + n + 1))

/-- `f(n+1,k)/f(n,k)` -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  2 * ((n : ℝ) + 1) * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2) * ((k : ℝ) + 3 * n + 3) /
    ((2 * (k : ℝ) + 6 * n + 3) * (2 * (k : ℝ) + 6 * n + 5) * (2 * (k : ℝ) + 6 * n + 7) *
      ((k : ℝ) + n + 1))

/-- the WZ certificate -/
noncomputable def R (n k : ℕ) : ℝ :=
  3 * ((k : ℝ) + 3 * n + 2) ^ 2 / ((2 * (k : ℝ) + 6 * n + 3) * (2 * (k : ℝ) + 6 * n + 5))

/-- `c(n) = f(n,0)` -/
noncomputable def c : ℕ → ℝ
  | 0 => 1
  | n + 1 => c n * ρ n 0

/-- the first member of the WZ pair -/
noncomputable def f (n : ℕ) : ℕ → ℝ
  | 0 => c n
  | k + 1 => f n k * σ n k

/-- the second member of the WZ pair -/
noncomputable def g (n k : ℕ) : ℝ := R n k * f n k

@[simp] lemma f_zero (n : ℕ) : f n 0 = c n := rfl
lemma f_succ (n k : ℕ) : f n (k + 1) = f n k * σ n k := rfl
lemma c_succ (n : ℕ) : c (n + 1) = c n * ρ n 0 := rfl

/-- compatibility of the two ratios -/
lemma compat (n k : ℕ) : ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1) := by
  unfold ρ σ
  push_cast
  have h1 : (2 * (k : ℝ) + 6 * n + 3) ≠ 0 := by positivity
  have h2 : (2 * (k : ℝ) + 6 * n + 5) ≠ 0 := by positivity
  have h3 : (2 * (k : ℝ) + 6 * n + 7) ≠ 0 := by positivity
  have h4 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h5 : (2 * (k : ℝ) + 6 * (n + 1) + 3) ≠ 0 := by positivity
  have h6 : ((k : ℝ) + (n + 1) + 1) ≠ 0 := by positivity
  have h7 : (2 * ((k : ℝ) + 1) + 6 * n + 3) ≠ 0 := by positivity
  have h8 : (2 * ((k : ℝ) + 1) + 6 * n + 5) ≠ 0 := by positivity
  have h9 : (2 * ((k : ℝ) + 1) + 6 * n + 7) ≠ 0 := by positivity
  have h10 : ((k : ℝ) + 1 + n + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma f_succ_left (n k : ℕ) : f (n + 1) k = f n k * ρ n k := by
  induction k with
  | zero => simp [c_succ]
  | succ k ih =>
    rw [f_succ, ih, f_succ, mul_assoc, compat, mul_assoc]

/-- the WZ equation, as an identity between the rational functions -/
lemma wz_rat (n k : ℕ) : ρ n k - 1 = R n (k + 1) * σ n k - R n k := by
  unfold ρ σ R
  push_cast
  have h1 : (2 * (k : ℝ) + 6 * n + 3) ≠ 0 := by positivity
  have h2 : (2 * (k : ℝ) + 6 * n + 5) ≠ 0 := by positivity
  have h3 : (2 * (k : ℝ) + 6 * n + 7) ≠ 0 := by positivity
  have h4 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h7 : (2 * ((k : ℝ) + 1) + 6 * n + 3) ≠ 0 := by positivity
  have h8 : (2 * ((k : ℝ) + 1) + 6 * n + 5) ≠ 0 := by positivity
  field_simp
  ring

/-- **The WZ equation** `f(n+1,k) - f(n,k) = g(n,k+1) - g(n,k)`. -/
theorem wz (n k : ℕ) : f (n + 1) k - f n k = g n (k + 1) - g n k := by
  unfold g
  rw [f_succ_left, f_succ]
  linear_combination (f n k) * wz_rat n k

/-- telescoping in `k` -/
lemma sum_k (n K : ℕ) : ∑ k ∈ range K, (f (n + 1) k - f n k) = g n K - g n 0 := by
  simp_rw [wz]
  exact Finset.sum_range_sub (fun k => g n k) K

/-- **The finite WZ identity**, valid for all `N`, `K`. -/
theorem finite_identity (N K : ℕ) :
    ∑ k ∈ range K, f N k - ∑ k ∈ range K, f 0 k
      = ∑ n ∈ range N, g n K - ∑ n ∈ range N, g n 0 := by
  have h1 : ∑ n ∈ range N, ∑ k ∈ range K, (f (n + 1) k - f n k)
      = ∑ n ∈ range N, (g n K - g n 0) := Finset.sum_congr rfl (fun n _ => sum_k n K)
  rw [Finset.sum_comm] at h1
  have h2 : ∀ k, ∑ n ∈ range N, (f (n + 1) k - f n k) = f N k - f 0 k :=
    fun k => Finset.sum_range_sub (fun n => f n k) N
  simp_rw [h2, Finset.sum_sub_distrib] at h1
  exact h1

/-! ### Bounds -/

lemma σ_abs_le (n k : ℕ) : |σ n k| ≤ 1 / 3 := by
  unfold σ
  have hd : (0 : ℝ) < 3 * (2 * (k : ℝ) + 6 * n + 3) * ((k : ℝ) + n + 1) := by positivity
  rw [neg_div, abs_neg, abs_of_nonneg (by positivity), div_le_iff₀ hd]
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  nlinarith [mul_nonneg hk hn, mul_nonneg hk hk, mul_nonneg hn hn]

lemma ρ0_pos (n : ℕ) : 0 < ρ n 0 := by
  unfold ρ
  push_cast
  positivity

lemma ρ0_le (n : ℕ) : ρ n 0 ≤ 1 / 3 := by
  unfold ρ
  push_cast
  have hd : (0 : ℝ) < (2 * (0 : ℝ) + 6 * n + 3) * (2 * (0 : ℝ) + 6 * n + 5) * (2 * (0 : ℝ) + 6 * n + 7) *
      ((0 : ℝ) + n + 1) := by positivity
  rw [div_le_iff₀ hd]
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have h2 : (0 : ℝ) ≤ (n : ℝ) ^ 2 := by positivity
  have h3 : (0 : ℝ) ≤ (n : ℝ) ^ 3 := by positivity
  have h4 : (0 : ℝ) ≤ (n : ℝ) ^ 4 := by positivity
  nlinarith [h2, h3, h4, hn]

lemma c_pos (n : ℕ) : 0 < c n := by
  induction n with
  | zero => simp [c]
  | succ n ih => rw [c_succ]; exact mul_pos ih (ρ0_pos n)

lemma c_le (n : ℕ) : c n ≤ (1 / 3) ^ n := by
  induction n with
  | zero => simp [c]
  | succ n ih =>
    rw [c_succ, pow_succ]
    exact mul_le_mul ih (ρ0_le n) (ρ0_pos n).le (by positivity)

lemma f_abs_le (n k : ℕ) : |f n k| ≤ (1 / 3) ^ n * (1 / 3) ^ k := by
  induction k with
  | zero => simpa [abs_of_pos (c_pos n)] using c_le n
  | succ k ih =>
    rw [f_succ, abs_mul, pow_succ, ← mul_assoc]
    exact mul_le_mul ih (σ_abs_le n k) (abs_nonneg _) (by positivity)

lemma R_nonneg (n k : ℕ) : 0 ≤ R n k := by
  unfold R; positivity

lemma R_le_one (n k : ℕ) : R n k ≤ 1 := by
  unfold R
  have hd : (0 : ℝ) < (2 * (k : ℝ) + 6 * n + 3) * (2 * (k : ℝ) + 6 * n + 5) := by positivity
  rw [div_le_one hd]
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  nlinarith [mul_nonneg hk hn, mul_nonneg hk hk, mul_nonneg hn hn]

lemma g_abs_le (n k : ℕ) : |g n k| ≤ (1 / 3) ^ n * (1 / 3) ^ k := by
  unfold g
  rw [abs_mul, abs_of_nonneg (R_nonneg n k)]
  calc R n k * |f n k| ≤ 1 * |f n k| := mul_le_mul_of_nonneg_right (R_le_one n k) (abs_nonneg _)
    _ = |f n k| := one_mul _
    _ ≤ _ := f_abs_le n k

lemma summable_f (n : ℕ) : Summable (f n) := by
  refine Summable.of_abs ?_
  refine Summable.of_nonneg_of_le (fun k => abs_nonneg _) (f_abs_le n) ?_
  exact (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left _

lemma tsum_f_abs_le (n : ℕ) : |∑' k, f n k| ≤ (1 / 3) ^ n * (3 / 2) := by
  have hs : Summable (fun k : ℕ => (1 / 3 : ℝ) ^ n * (1 / 3) ^ k) :=
    (summable_geometric_of_lt_one (by norm_num) (by norm_num)).mul_left _
  have h1 : |∑' k, f n k| ≤ ∑' k, |f n k| := by
    have h := norm_tsum_le_tsum_norm (f := f n) (summable_f n).norm
    simp only [Real.norm_eq_abs] at h
    exact h
  have h2 : ∑' k, |f n k| ≤ ∑' k, (1 / 3 : ℝ) ^ n * (1 / 3) ^ k :=
    (summable_f n).abs.tsum_le_tsum (f_abs_le n) hs
  have h3 : ∑' k, (1 / 3 : ℝ) ^ n * (1 / 3) ^ k = (1 / 3) ^ n * (3 / 2) := by
    rw [tsum_mul_left, tsum_geometric_of_lt_one (by norm_num) (by norm_num)]
    norm_num
  linarith

/-! ### The boundary sum `∑ f(0,k)` -/

lemma f0_eq (k : ℕ) : f 0 k = (-1 / 3) ^ k / (2 * (k : ℝ) + 1) := by
  induction k with
  | zero => simp [c]
  | succ k ih =>
    rw [f_succ, ih]
    unfold σ
    push_cast
    have h1 : (2 * (k : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (2 * ((k : ℝ) + 1) + 1) ≠ 0 := by positivity
    have h3 : ((k : ℝ) + 1) ≠ 0 := by positivity
    have h4 : (2 * (k : ℝ) + 3) ≠ 0 := by positivity
    rw [pow_succ]
    field_simp
    ring

lemma hasSum_f0 : HasSum (f 0) (Real.pi / (2 * Real.sqrt 3)) := by
  have hs3 : (0 : ℝ) < Real.sqrt 3 := Real.sqrt_pos.mpr (by norm_num)
  have hx : ‖(Real.sqrt 3)⁻¹‖ < 1 := by
    rw [Real.norm_eq_abs, abs_of_pos (inv_pos.mpr hs3)]
    rw [inv_lt_one_iff₀]
    right
    rw [show (1 : ℝ) = Real.sqrt 1 by simp]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  have h := Real.hasSum_arctan hx
  rw [Real.arctan_inv_sqrt_three] at h
  have h' := h.mul_left (Real.sqrt 3)
  have hsq : Real.sqrt 3 * Real.sqrt 3 = 3 := Real.mul_self_sqrt (by norm_num)
  have hfun : f 0 = fun k : ℕ =>
      Real.sqrt 3 * ((-1) ^ k * (Real.sqrt 3)⁻¹ ^ (2 * k + 1) / ((2 * k + 1 : ℕ) : ℝ)) := by
    funext k
    rw [f0_eq]
    have hpow : ((Real.sqrt 3)⁻¹) ^ (2 * k + 1) = (1 / 3) ^ k * (Real.sqrt 3)⁻¹ := by
      rw [pow_succ, pow_mul]
      congr 1
      rw [inv_pow, Real.sq_sqrt (by norm_num)]
      simp
    rw [hpow]
    have h1 : (2 * (k : ℝ) + 1) ≠ 0 := by positivity
    have h3 : (-1 / 3 : ℝ) ^ k = (-1) ^ k * (1 / 3) ^ k := by
      rw [← mul_pow]; norm_num
    rw [h3]
    push_cast
    have hs3' : Real.sqrt 3 ≠ 0 := hs3.ne'
    field_simp
  have hval : Real.pi / (2 * Real.sqrt 3) = Real.sqrt 3 * (Real.pi / 6) := by
    rw [div_eq_iff (by positivity)]
    linear_combination (-(Real.pi / 3)) * hsq
  rw [hfun, hval]
  exact h'

/-! ### The limit argument -/

/-- letting `K → ∞` in the finite identity -/
lemma partial_sum_g (N : ℕ) :
    ∑ n ∈ range N, g n 0 = ∑' k, f 0 k - ∑' k, f N k := by
  have hL : Tendsto (fun K => ∑ k ∈ range K, f N k - ∑ k ∈ range K, f 0 k) atTop
      (𝓝 (∑' k, f N k - ∑' k, f 0 k)) :=
    ((summable_f N).hasSum.tendsto_sum_nat).sub ((summable_f 0).hasSum.tendsto_sum_nat)
  have hgK : ∀ n, Tendsto (fun K => g n K) atTop (𝓝 0) := by
    intro n
    refine squeeze_zero_norm (a := fun K : ℕ => (1 / 3 : ℝ) ^ n * (1 / 3) ^ K) (fun K => ?_) ?_
    · rw [Real.norm_eq_abs]; exact (g_abs_le n K)
    · have : Tendsto (fun K : ℕ => (1 / 3 : ℝ) ^ n * (1 / 3) ^ K) atTop (𝓝 ((1 / 3) ^ n * 0)) :=
        (tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)).const_mul _
      simpa using this
  have hR : Tendsto (fun K => ∑ n ∈ range N, g n K - ∑ n ∈ range N, g n 0) atTop
      (𝓝 (∑ n ∈ range N, (0 : ℝ) - ∑ n ∈ range N, g n 0)) :=
    (tendsto_finsetSum _ (fun n _ => hgK n)).sub tendsto_const_nhds
  have heq : (fun K => ∑ k ∈ range K, f N k - ∑ k ∈ range K, f 0 k)
      = (fun K => ∑ n ∈ range N, g n K - ∑ n ∈ range N, g n 0) :=
    funext (fun K => finite_identity N K)
  rw [heq] at hL
  have := tendsto_nhds_unique hL hR
  simp at this
  linarith

lemma g0_nonneg (n : ℕ) : 0 ≤ g n 0 := by
  unfold g
  exact mul_nonneg (R_nonneg n 0) (by simpa using (c_pos n).le)

/-- `∑_n g(n,0) = π/(2√3)` -/
theorem hasSum_g0 : HasSum (fun n => g n 0) (Real.pi / (2 * Real.sqrt 3)) := by
  rw [hasSum_iff_tendsto_nat_of_nonneg g0_nonneg]
  have h0 : ∑' k, f 0 k = Real.pi / (2 * Real.sqrt 3) := hasSum_f0.tsum_eq
  have hlim : Tendsto (fun N => ∑' k, f N k) atTop (𝓝 0) := by
    refine squeeze_zero_norm (a := fun N : ℕ => (1 / 3 : ℝ) ^ N * (3 / 2)) (fun N => ?_) ?_
    · rw [Real.norm_eq_abs]; exact tsum_f_abs_le N
    · have : Tendsto (fun N : ℕ => (1 / 3 : ℝ) ^ N * (3 / 2)) atTop (𝓝 (0 * (3 / 2))) :=
        (tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)).mul_const _
      simpa using this
  have : Tendsto (fun N => ∑' k, f 0 k - ∑' k, f N k) atTop
      (𝓝 (Real.pi / (2 * Real.sqrt 3) - 0)) := by
    rw [← h0]
    exact tendsto_const_nhds.sub hlim
  simpa [partial_sum_g] using this

/-! ### Identification of `g(n,0)` with Sun's summand -/

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
    have : (3 * (n : ℝ) + 1 + 1) = 3 * (n : ℝ) + 2 := by ring
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

lemma c_eq (n : ℕ) : c n = 16 ^ n / ((6 * (n : ℝ) + 1) * (Nat.centralBinom (3 * n) : ℝ)) := by
  induction n with
  | zero => simp [c]
  | succ n ih =>
    rw [c_succ, ih]
    unfold ρ
    have hcb : (0 : ℝ) < (Nat.centralBinom (3 * n) : ℝ) := by exact_mod_cast Nat.centralBinom_pos _
    have hcb' : (0 : ℝ) < (Nat.centralBinom (3 * n + 3) : ℝ) := by exact_mod_cast Nat.centralBinom_pos _
    have hstep := centralBinom_step n
    rw [show 3 * (n + 1) = 3 * n + 3 by ring]
    push_cast
    have h1 : (6 * (n : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (6 * ((n : ℝ) + 1) + 1) ≠ 0 := by positivity
    have h3 : (2 * (0 : ℝ) + 6 * n + 3) ≠ 0 := by positivity
    have h4 : (2 * (0 : ℝ) + 6 * n + 5) ≠ 0 := by positivity
    have h5 : (2 * (0 : ℝ) + 6 * n + 7) ≠ 0 := by positivity
    have h6 : ((0 : ℝ) + n + 1) ≠ 0 := by positivity
    rw [div_mul_div_comm, div_eq_div_iff (by positivity) (by positivity)]
    have key : ((Nat.centralBinom (3 * n + 3) : ℕ) : ℝ) * ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3))
        = 8 * ((6 * (n : ℝ) + 1) * (6 * n + 3) * (6 * n + 5)) * (Nat.centralBinom (3 * n) : ℝ) := hstep
    rw [pow_succ]
    linear_combination (2 * 16 ^ n * ((n : ℝ) + 1) * (6 * n + 7)) * key

/-- Sun's summand -/
noncomputable def sunTerm (k : ℕ) : ℝ :=
  (3 * (k : ℝ) + 2) ^ 2 * 16 ^ k /
    ((2 * k + 1) * (6 * k + 1) * (6 * k + 5) * (Nat.centralBinom (3 * k) : ℝ))

lemma g0_eq (n : ℕ) : g n 0 = sunTerm n := by
  unfold g sunTerm R
  rw [f_zero, c_eq]
  have hcb : (0 : ℝ) < (Nat.centralBinom (3 * n) : ℝ) := by exact_mod_cast Nat.centralBinom_pos _
  push_cast
  have h1 : (6 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h3 : (6 * (n : ℝ) + 5) ≠ 0 := by positivity
  have h4 : (2 * (0 : ℝ) + 6 * n + 3) ≠ 0 := by positivity
  have h5 : (2 * (0 : ℝ) + 6 * n + 5) ≠ 0 := by positivity
  field_simp
  ring

/-- **Sun's conjecture** (arXiv:2506.01870v6, Conj. 3.2):
`∑_{k≥0} (3k+2)² 16^k / ((2k+1)(6k+1)(6k+5) C(6k,3k)) = π/(2√3)`. -/
theorem sun_series_16 :
    HasSum (fun k : ℕ => (3 * (k : ℝ) + 2) ^ 2 * 16 ^ k /
      ((2 * k + 1) * (6 * k + 1) * (6 * k + 5) * (Nat.centralBinom (3 * k) : ℝ)))
      (Real.pi / (2 * Real.sqrt 3)) := by
  have h := hasSum_g0
  simp_rw [g0_eq] at h
  exact h

end SunSeries16
