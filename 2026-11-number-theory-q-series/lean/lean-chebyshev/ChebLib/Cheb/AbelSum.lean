/-
ChebLib/Cheb/AbelSum.lean — an Abelian theorem for the weight `e^{−n/x}`.

If the partial sums `A(y) = Σ_{n ≤ y} a_n` satisfy `A(y)/√y → 1`, then
`(Σ_n a_n e^{−n/x}) / √x → Γ(3/2) = √π/2`.

Proof idea: `Σ_n a_n e^{−n/x} = (1/x) ∫_0^∞ A(y) e^{−y/x} dy` (as `e^{−n/x} = (1/x) ∫_n^∞ e^{−y/x} dy`
and Fubini), `(1/x) ∫_0^∞ √y e^{−y/x} dy = Γ(3/2) √x`, and `|A(y) − √y| ≤ ε √y` for `y ≥ Y`,
`≤ C_Y` for `y ≤ Y`.
-/
import Mathlib

open MeasureTheory Real Filter Topology

noncomputable section

namespace ChebLib

lemma integral_exp_neg_div_Ioi {x : ℝ} (hx : 0 < x) (c : ℝ) :
    ∫ y in Set.Ioi c, Real.exp (-(y / x)) = x * Real.exp (-(c / x)) := by
  have hneg : -(1 / x) < 0 := by
    have : 0 < 1 / x := by positivity
    linarith
  have h := integral_exp_mul_Ioi hneg c
  have e : (fun y : ℝ => Real.exp (-(y / x))) = fun y => Real.exp (-(1 / x) * y) := by
    funext y
    congr 1
    ring
  rw [e, h]
  have e2 : -(1 / x) * c = -(c / x) := by ring
  rw [e2]
  field_simp

lemma integrableOn_exp_neg_div_Ioi {x : ℝ} (hx : 0 < x) (c : ℝ) :
    IntegrableOn (fun y : ℝ => Real.exp (-(y / x))) (Set.Ioi c) := by
  have hneg : -(1 / x) < 0 := by
    have : 0 < 1 / x := by positivity
    linarith
  have h := integrableOn_exp_mul_Ioi hneg c
  have e : (fun y : ℝ => Real.exp (-(y / x))) = fun y => Real.exp (-(1 / x) * y) := by
    funext y
    congr 1
    ring
  rw [e]
  exact h

lemma integral_sqrt_mul_exp_neg_div {x : ℝ} (hx : 0 < x) :
    ∫ y in Set.Ioi (0 : ℝ), Real.sqrt y * Real.exp (-(y / x)) =
      Real.sqrt π / 2 * (x * Real.sqrt x) := by
  have hr : 0 < 1 / x := by positivity
  have h := Real.integral_rpow_mul_exp_neg_mul_Ioi (a := 3 / 2) (r := 1 / x) (by norm_num) hr
  have e : ∀ y : ℝ, Real.sqrt y * Real.exp (-(y / x)) =
      y ^ ((3 : ℝ) / 2 - 1) * Real.exp (-(1 / x * y)) := by
    intro y
    rw [Real.sqrt_eq_rpow, show (3 : ℝ) / 2 - 1 = 1 / 2 by norm_num,
      show -(y / x) = -(1 / x * y) by ring]
  simp_rw [e]
  rw [h]
  have hG : Real.Gamma (3 / 2) = Real.sqrt π / 2 := by
    rw [show (3 : ℝ) / 2 = 1 / 2 + 1 by norm_num, Real.Gamma_add_one (by norm_num),
      Real.Gamma_one_half_eq]
    ring
  have hx' : (1 / (1 / x)) ^ ((3 : ℝ) / 2) = x * Real.sqrt x := by
    rw [one_div_one_div, show (3 : ℝ) / 2 = 1 + 1 / 2 by norm_num, Real.rpow_add hx,
      Real.rpow_one, Real.sqrt_eq_rpow]
  rw [hG, hx']
  ring

lemma integrableOn_sqrt_mul_exp_neg_div {x : ℝ} (hx : 0 < x) :
    IntegrableOn (fun y : ℝ => Real.sqrt y * Real.exp (-(y / x))) (Set.Ioi 0) := by
  refine Integrable.of_integral_ne_zero ?_
  rw [integral_sqrt_mul_exp_neg_div hx]
  positivity

/-- Absolute convergence of the Abel sums of a linearly bounded sequence. -/
theorem summable_mul_exp_neg_div (a : ℕ → ℝ) (C : ℝ) (ha : ∀ n, |a n| ≤ C * n) {x : ℝ}
    (hx : 0 < x) : Summable (fun n : ℕ => a n * Real.exp (-(n / x))) := by
  have hr : ‖Real.exp (-(1 / x))‖ < 1 := by
    rw [Real.norm_eq_abs, Real.abs_exp, Real.exp_lt_one_iff]
    have : 0 < 1 / x := by positivity
    linarith
  have hs := (summable_pow_mul_geometric_of_norm_lt_one 1 hr).mul_left C
  refine Summable.of_norm_bounded hs fun n => ?_
  have e : Real.exp (-(n / x)) = Real.exp (-(1 / x)) ^ n := by
    rw [← Real.exp_nat_mul]
    congr 1
    ring
  rw [norm_mul, Real.norm_eq_abs, Real.norm_eq_abs, Real.abs_exp, e, pow_one, ← mul_assoc]
  exact mul_le_mul_of_nonneg_right (ha n) (by positivity)

lemma abel_sum_eq_integral (a : ℕ → ℝ) (C : ℝ) (ha : ∀ n, |a n| ≤ C * n) {x : ℝ}
    (hx : 0 < x) :
    x * ∑' n : ℕ, a n * Real.exp (-(n / x)) =
      ∫ y in Set.Ioi (0 : ℝ), (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) * Real.exp (-(y / x)) := by
  set F : ℕ → ℝ → ℝ := fun n =>
    (Set.Ici (n : ℝ)).indicator (fun y => a n * Real.exp (-(y / x))) with hF
  have hInt : ∀ n, Integrable (F n) := by
    intro n
    rw [hF, integrable_indicator_iff measurableSet_Ici, integrableOn_Ici_iff_integrableOn_Ioi]
    exact (integrableOn_exp_neg_div_Ioi hx _).const_mul _
  have hval : ∀ n, ∫ y, F n y = a n * (x * Real.exp (-(n / x))) := by
    intro n
    rw [hF, integral_indicator measurableSet_Ici, integral_Ici_eq_integral_Ioi,
      integral_const_mul, integral_exp_neg_div_Ioi hx]
  have hnorm : ∀ n, ∫ y, ‖F n y‖ = x * (|a n| * Real.exp (-(n / x))) := by
    intro n
    have : (fun y => ‖F n y‖) =
        (Set.Ici (n : ℝ)).indicator (fun y => |a n| * Real.exp (-(y / x))) := by
      funext y
      rw [hF, norm_indicator_eq_indicator_norm]
      congr 1
      funext z
      rw [norm_mul, Real.norm_eq_abs, Real.norm_eq_abs, Real.abs_exp]
    rw [this, integral_indicator measurableSet_Ici, integral_Ici_eq_integral_Ioi,
      integral_const_mul, integral_exp_neg_div_Ioi hx]
    ring
  have habs : ∀ n, |(fun n => |a n|) n| ≤ C * n := fun n => by
    simpa using ha n
  have hsum : Summable fun n => ∫ y, ‖F n y‖ := by
    simp_rw [hnorm]
    exact (summable_mul_exp_neg_div (fun n => |a n|) C habs hx).mul_left x
  have hS := hasSum_integral_of_summable_integral_norm hInt hsum
  have htsum : ∀ y : ℝ, ∑' n, F n y =
      (Set.Ici (0 : ℝ)).indicator
        (fun y => (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) * Real.exp (-(y / x))) y := by
    intro y
    by_cases hy : 0 ≤ y
    · rw [Set.indicator_of_mem (Set.mem_Ici.mpr hy)]
      rw [tsum_eq_sum (s := Finset.Icc 0 ⌊y⌋₊), Finset.sum_mul]
      · refine Finset.sum_congr rfl fun n hn => ?_
        have : (n : ℝ) ≤ y := (Nat.le_floor_iff hy).mp (Finset.mem_Icc.mp hn).2
        simp [hF, this]
      · intro n hn
        have : ¬ (n : ℝ) ≤ y := by
          intro h
          exact hn (Finset.mem_Icc.mpr ⟨Nat.zero_le _, (Nat.le_floor_iff hy).mpr h⟩)
        simp [hF, this]
    · rw [Set.indicator_of_notMem (by simpa using hy)]
      have : ∀ n, F n y = 0 := by
        intro n
        have : ¬ (n : ℝ) ≤ y := by
          intro h
          exact hy ((Nat.cast_nonneg n).trans h)
        simp [hF, this]
      simp [this]
  have h1 : ∑' n, a n * (x * Real.exp (-(n / x))) =
      ∫ y in Set.Ioi (0 : ℝ), (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) * Real.exp (-(y / x)) := by
    rw [← integral_Ici_eq_integral_Ioi, ← integral_indicator measurableSet_Ici]
    simp_rw [← htsum, ← hval]
    exact hS.tsum_eq
  rw [← h1, ← tsum_mul_left]
  congr 1
  funext n
  ring

lemma measurable_partial_sum_mul_exp (a : ℕ → ℝ) (x : ℝ) :
    Measurable (fun y : ℝ => (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) * Real.exp (-(y / x))) := by
  have h1 : Measurable (fun y : ℝ => ∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) :=
    (measurable_from_nat (f := fun m : ℕ => ∑ n ∈ Finset.Icc 0 m, a n)).comp Nat.measurable_floor
  exact h1.mul (by fun_prop)

lemma abel_sum_bound (a : ℕ → ℝ) (C : ℝ) (ha : ∀ n, |a n| ≤ C * n) {x : ℝ}
    (hx : 0 < x) (M ε : ℝ)
    (hb : ∀ y : ℝ, 0 < y → |(∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) - Real.sqrt y| ≤ M + ε * Real.sqrt y) :
    |x * (∑' n : ℕ, a n * Real.exp (-(n / x))) - Real.sqrt π / 2 * (x * Real.sqrt x)| ≤
      M * x + ε * (Real.sqrt π / 2 * (x * Real.sqrt x)) := by
  have hI1 := integrableOn_exp_neg_div_Ioi hx 0
  have hI2 := integrableOn_sqrt_mul_exp_neg_div hx
  have hG : IntegrableOn (fun y : ℝ => (M + ε * Real.sqrt y) * Real.exp (-(y / x)))
      (Set.Ioi 0) := by
    have : (fun y : ℝ => (M + ε * Real.sqrt y) * Real.exp (-(y / x))) =
        fun y => M * Real.exp (-(y / x)) + ε * (Real.sqrt y * Real.exp (-(y / x))) := by
      funext y; ring
    rw [this]
    exact (hI1.const_mul M).add (hI2.const_mul ε)
  have hGval : ∫ y in Set.Ioi (0 : ℝ), (M + ε * Real.sqrt y) * Real.exp (-(y / x)) =
      M * x + ε * (Real.sqrt π / 2 * (x * Real.sqrt x)) := by
    have : (fun y : ℝ => (M + ε * Real.sqrt y) * Real.exp (-(y / x))) =
        fun y => M * Real.exp (-(y / x)) + ε * (Real.sqrt y * Real.exp (-(y / x))) := by
      funext y; ring
    rw [this, integral_add (hI1.const_mul M) (hI2.const_mul ε), integral_const_mul,
      integral_const_mul, integral_exp_neg_div_Ioi hx, integral_sqrt_mul_exp_neg_div hx]
    simp
  have hA : IntegrableOn
      (fun y : ℝ => (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) * Real.exp (-(y / x))) (Set.Ioi 0) := by
    refine Integrable.mono' (hG.add hI2) (measurable_partial_sum_mul_exp a x).aestronglyMeasurable ?_
    rw [ae_restrict_iff' measurableSet_Ioi]
    refine Eventually.of_forall fun y hy => ?_
    have h1 := hb y hy
    have h2 : |∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n| ≤ (M + ε * Real.sqrt y) + Real.sqrt y := by
      have := abs_le.mp h1
      have hs := Real.sqrt_nonneg y
      rw [abs_le]
      constructor <;> linarith [this.1, this.2]
    have he := Real.exp_pos (-(y / x))
    rw [norm_mul, Real.norm_eq_abs, Real.norm_eq_abs, Real.abs_exp]
    calc |∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n| * Real.exp (-(y / x))
        ≤ ((M + ε * Real.sqrt y) + Real.sqrt y) * Real.exp (-(y / x)) :=
          mul_le_mul_of_nonneg_right h2 he.le
      _ = _ := by simp only [Pi.add_apply]; ring
  have hdiff : x * (∑' n : ℕ, a n * Real.exp (-(n / x))) - Real.sqrt π / 2 * (x * Real.sqrt x) =
      ∫ y in Set.Ioi (0 : ℝ), ((∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) * Real.exp (-(y / x)) -
        Real.sqrt y * Real.exp (-(y / x))) := by
    rw [abel_sum_eq_integral a C ha hx, ← integral_sqrt_mul_exp_neg_div hx, ← integral_sub hA hI2]
  rw [hdiff, ← hGval, ← Real.norm_eq_abs]
  refine norm_integral_le_of_norm_le hG ?_
  rw [ae_restrict_iff' measurableSet_Ioi]
  refine Eventually.of_forall fun y hy => ?_
  rw [← sub_mul, norm_mul, Real.norm_eq_abs, Real.norm_eq_abs, Real.abs_exp]
  exact mul_le_mul_of_nonneg_right (hb y hy) (Real.exp_pos _).le

/-- **Abelian theorem.**  If `|a_n| ≤ C n` and `(Σ_{n ≤ y} a_n)/√y → 1` as `y → ∞`, then
`(Σ_n a_n e^{−n/x})/√x → √π/2` as `x → ∞`. -/
theorem tendsto_abel_sum_div_sqrt (a : ℕ → ℝ) (C : ℝ) (ha : ∀ n, |a n| ≤ C * n)
    (h : Tendsto (fun y : ℝ => (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) / Real.sqrt y) atTop (𝓝 1)) :
    Tendsto (fun x : ℝ => (∑' n : ℕ, a n * Real.exp (-(n / x))) / Real.sqrt x) atTop
      (𝓝 (Real.sqrt π / 2)) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  have hδ : 0 < ε / 2 := by positivity
  obtain ⟨Y₀, hY₀⟩ := (Metric.tendsto_atTop.mp h) (ε / 2) hδ
  set Y : ℝ := max Y₀ 1 with hYdef
  have hY1 : 1 ≤ Y := le_max_right _ _
  set M : ℝ := ∑ n ∈ Finset.Icc 0 ⌊Y⌋₊, |a n| + Real.sqrt Y with hMdef
  have hM0 : 0 ≤ M := by
    have : 0 ≤ ∑ n ∈ Finset.Icc 0 ⌊Y⌋₊, |a n| := Finset.sum_nonneg fun _ _ => abs_nonneg _
    have := Real.sqrt_nonneg Y
    linarith
  have hb : ∀ y : ℝ, 0 < y →
      |(∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) - Real.sqrt y| ≤ M + ε / 2 * Real.sqrt y := by
    intro y hy
    have hsy : 0 < Real.sqrt y := Real.sqrt_pos.mpr hy
    by_cases hyY : Y ≤ y
    · have h1 := hY₀ y ((le_max_left _ _).trans hyY)
      rw [Real.dist_eq] at h1
      have h2 : (∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) - Real.sqrt y =
          ((∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) / Real.sqrt y - 1) * Real.sqrt y := by
        field_simp
      rw [h2, abs_mul, abs_of_pos hsy]
      have := mul_le_mul_of_nonneg_right h1.le hsy.le
      linarith
    · have hyY' : y ≤ Y := (not_le.mp hyY).le
      have h1 : |∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n| ≤ ∑ n ∈ Finset.Icc 0 ⌊Y⌋₊, |a n| := by
        refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
        refine Finset.sum_le_sum_of_subset_of_nonneg ?_ fun _ _ _ => abs_nonneg _
        exact Finset.Icc_subset_Icc le_rfl (Nat.floor_le_floor hyY')
      have h2 : Real.sqrt y ≤ Real.sqrt Y := Real.sqrt_le_sqrt hyY'
      have h3 : |(∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n) - Real.sqrt y| ≤
          |∑ n ∈ Finset.Icc 0 ⌊y⌋₊, a n| + |Real.sqrt y| := abs_sub _ _
      rw [abs_of_pos hsy] at h3
      have : 0 ≤ ε / 2 * Real.sqrt y := by positivity
      linarith
  have hc0 : 0 ≤ Real.sqrt π / 2 := by positivity
  have hc1 : Real.sqrt π / 2 ≤ 1 := by
    have : Real.sqrt π ≤ 2 := by
      rw [Real.sqrt_le_iff]
      refine ⟨by norm_num, ?_⟩
      have := Real.pi_le_four
      linarith
    linarith
  refine ⟨max 1 ((2 * M / ε) ^ 2 + 1), fun x hxN => ?_⟩
  have hx1 : 1 ≤ x := (le_max_left _ _).trans hxN
  have hx : 0 < x := by linarith
  have hsx : 0 < Real.sqrt x := Real.sqrt_pos.mpr hx
  have hxl : (2 * M / ε) ^ 2 < x := by
    have := (le_max_right _ _).trans hxN
    linarith
  have hsl : 2 * M / ε < Real.sqrt x := (Real.lt_sqrt (by positivity)).mpr hxl
  have hMs : M < ε / 2 * Real.sqrt x := by
    rw [div_lt_iff₀ hε] at hsl
    linarith
  have key := abel_sum_bound a C ha hx M (ε / 2) hb
  have hpos : 0 < x * Real.sqrt x := by positivity
  rw [Real.dist_eq]
  have e : (∑' n : ℕ, a n * Real.exp (-(n / x))) / Real.sqrt x - Real.sqrt π / 2 =
      (x * (∑' n : ℕ, a n * Real.exp (-(n / x))) - Real.sqrt π / 2 * (x * Real.sqrt x)) /
        (x * Real.sqrt x) := by
    field_simp
  rw [e, abs_div, abs_of_pos hpos, div_lt_iff₀ hpos]
  calc _ ≤ M * x + ε / 2 * (Real.sqrt π / 2 * (x * Real.sqrt x)) := key
    _ ≤ M * x + ε / 2 * (1 * (x * Real.sqrt x)) := by gcongr
    _ < ε / 2 * Real.sqrt x * x + ε / 2 * (1 * (x * Real.sqrt x)) := by gcongr
    _ = ε * (x * Real.sqrt x) := by ring

end ChebLib
