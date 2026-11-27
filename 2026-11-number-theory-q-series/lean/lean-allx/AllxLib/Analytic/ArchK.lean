/-
AllxLib/Analytic/ArchK.lean — the archimedean term of the explicit formula for the test function
`Kt t` is at most `log t + 3` (for `t ≥ 2`).

By `arch_weil` the term is
`k(0)(log(4/π) − γ) − I(3/2) + Σ_n (k(0)/(n+1) − I(2n+7/2))`,  `I(b) = ∫ Kt t u e^{−b|u|} du`,
`k(0) = t − 1`.  Lower bounds for `I(b)`: the part `u ≤ 0` is an exact sum of exponential integrals,
the part `u ≥ 0` is bounded below by `∫₀^{log t} e^{u/2}(t − e^u) e^{−bu} du`.
-/
import AllxLib.Analytic.ArchWeil
import AllxLib.Analytic.KernelK

open MeasureTheory Set Filter Topology

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace Allx

open Real in
lemma integral_exp_neg_mul (c L : ℝ) (hc : c ≠ 0) :
    ∫ u in (0 : ℝ)..L, exp (-c * u) = (1 - exp (-c * L)) / c := by
  have h : ∀ u ∈ uIcc (0 : ℝ) L, HasDerivAt (fun u => -exp (-c * u) / c) (exp (-c * u)) u := by
    intro u _
    have := ((hasDerivAt_exp_mul (-c) u).neg).div_const c
    refine this.congr_deriv ?_
    field_simp
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt h
    (Continuous.intervalIntegrable (by fun_prop) _ _)]
  simp only [mul_zero, exp_zero]
  ring

open Real in
/-- The part `u ≤ 0` of `∫ Kt t u e^{−b|u|} du`. -/
theorem integral_Iic_Kt {t : ℝ} (ht : 2 ≤ t) {b : ℝ} (hb : 0 < b) :
    ∫ u in Iic (0 : ℝ), Kt t u * exp (-b * |u|)
      = t * (3 / (b + 3 / 2) - 3 / (b + 5 / 2) + 1 / (b + 7 / 2))
        - (3 / (b + 5 / 2) - 3 / (b + 7 / 2) + 1 / (b + 9 / 2)) := by
  have hcong : ∀ u ∈ Iic (0 : ℝ), Kt t u * exp (-b * |u|)
      = 3 * t * exp ((b + 3 / 2) * u) + (-(3 * t + 3)) * exp ((b + 5 / 2) * u)
        + (t + 3) * exp ((b + 7 / 2) * u) + (-1) * exp ((b + 9 / 2) * u) := by
    intro u hu
    have hu' : u ≤ 0 := hu
    rw [Kt_of_nonpos ht hu', abs_of_nonpos hu',
      show -b * -u = b * u by ring,
      show (b + 3 / 2) * u = b * u + (u / 2 + u) by ring,
      show (b + 5 / 2) * u = b * u + (u / 2 + (u + u)) by ring,
      show (b + 7 / 2) * u = b * u + (u / 2 + (u + (u + u))) by ring,
      show (b + 9 / 2) * u = b * u + (u / 2 + (u + (u + (u + u)))) by ring]
    simp only [exp_add]
    ring
  rw [setIntegral_congr_fun measurableSet_Iic hcong]
  have i1 : IntegrableOn (fun u : ℝ => 3 * t * exp ((b + 3 / 2) * u)) (Iic 0) :=
    (integrableOn_exp_mul_Iic (by linarith) 0).const_mul _
  have i2 : IntegrableOn (fun u : ℝ => (-(3 * t + 3)) * exp ((b + 5 / 2) * u)) (Iic 0) :=
    (integrableOn_exp_mul_Iic (by linarith) 0).const_mul _
  have i3 : IntegrableOn (fun u : ℝ => (t + 3) * exp ((b + 7 / 2) * u)) (Iic 0) :=
    (integrableOn_exp_mul_Iic (by linarith) 0).const_mul _
  have i4 : IntegrableOn (fun u : ℝ => (-1) * exp ((b + 9 / 2) * u)) (Iic 0) :=
    (integrableOn_exp_mul_Iic (by linarith) 0).const_mul _
  have i12 : IntegrableOn (fun u : ℝ => 3 * t * exp ((b + 3 / 2) * u)
      + (-(3 * t + 3)) * exp ((b + 5 / 2) * u)) (Iic 0) := i1.add i2
  have i123 : IntegrableOn (fun u : ℝ => 3 * t * exp ((b + 3 / 2) * u)
      + (-(3 * t + 3)) * exp ((b + 5 / 2) * u) + (t + 3) * exp ((b + 7 / 2) * u)) (Iic 0) :=
    i12.add i3
  rw [integral_add i123 i4, integral_add i12 i3, integral_add i1 i2, integral_const_mul,
    integral_const_mul, integral_const_mul, integral_const_mul,
    integral_exp_mul_Iic (by linarith) 0, integral_exp_mul_Iic (by linarith) 0,
    integral_exp_mul_Iic (by linarith) 0, integral_exp_mul_Iic (by linarith) 0]
  simp only [mul_zero, exp_zero]
  ring

open Real in
theorem integrable_Kt_mul_exp {t : ℝ} (ht : 2 ≤ t) (b : ℝ) (hb : 0 ≤ b) :
    Integrable (fun u : ℝ => Kt t u * exp (-b * |u|)) := by
  refine (integrable_Kt (by linarith)).mul_bdd (c := 1) (by fun_prop) ?_
  refine Eventually.of_forall fun u => ?_
  rw [Real.norm_eq_abs, abs_of_pos (exp_pos _)]
  apply exp_le_one_iff.mpr
  nlinarith [abs_nonneg u]

open Real in
/-- The part `u ≥ 0` of `∫ Kt t u e^{−b|u|} du`, from below. -/
theorem integral_Ioi_Kt_ge {t : ℝ} (ht : 2 ≤ t) {b : ℝ} (hb : 0 ≤ b) :
    ∫ u in (0 : ℝ)..log t, (t * exp (-(b - 1 / 2) * u) - exp (-(b - 3 / 2) * u))
      ≤ ∫ u in Ioi (0 : ℝ), Kt t u * exp (-b * |u|) := by
  have hL : 0 ≤ log t := log_nonneg (by linarith)
  have hfi := integrable_Kt_mul_exp ht b hb
  rw [intervalIntegral.integral_of_le hL]
  have hg : IntegrableOn (fun u : ℝ => t * exp (-(b - 1 / 2) * u) - exp (-(b - 3 / 2) * u))
      (Ioc 0 (log t)) :=
    (Continuous.integrableOn_Icc (by fun_prop)).mono_set Ioc_subset_Icc_self
  calc ∫ u in Ioc 0 (log t), (t * exp (-(b - 1 / 2) * u) - exp (-(b - 3 / 2) * u))
      ≤ ∫ u in Ioc 0 (log t), Kt t u * exp (-b * |u|) := by
        refine setIntegral_mono_on hg hfi.integrableOn measurableSet_Ioc fun u hu => ?_
        have hu0 : 0 ≤ u := hu.1.le
        rw [abs_of_nonneg hu0]
        have h1 := le_Kt_of_nonneg t hu0
        have h2 := mul_le_mul_of_nonneg_right h1 (exp_pos (-b * u)).le
        refine le_trans (le_of_eq ?_) h2
        rw [show -(b - 1 / 2) * u = u / 2 + -b * u by ring,
          show -(b - 3 / 2) * u = u / 2 + u + -b * u by ring]
        simp only [exp_add]
        ring
    _ ≤ ∫ u in Ioi 0, Kt t u * exp (-b * |u|) := by
        refine setIntegral_mono_set hfi.integrableOn ?_ ?_
        · exact Eventually.of_forall fun u => mul_nonneg (Kt_nonneg ht u) (exp_pos _).le
        · exact Eventually.of_forall fun u hu => hu.1

open Real in
theorem integral_Kt_split {t : ℝ} (ht : 2 ≤ t) {b : ℝ} (hb : 0 ≤ b) :
    ∫ u : ℝ, Kt t u * exp (-b * |u|)
      = (∫ u in Iic (0 : ℝ), Kt t u * exp (-b * |u|))
        + ∫ u in Ioi (0 : ℝ), Kt t u * exp (-b * |u|) := by
  have hfi := integrable_Kt_mul_exp ht b hb
  rw [← intervalIntegral.integral_Iic_add_Ioi hfi.integrableOn hfi.integrableOn]

open Real in
/-- `I(3/2) ≥ (29/20) t − 79/60 − log t`. -/
theorem integral_Kt_three_halves_ge {t : ℝ} (ht : 2 ≤ t) :
    29 / 20 * t - 79 / 60 - log t ≤ ∫ u : ℝ, Kt t u * exp (-(3 / 2) * |u|) := by
  have ht0 : 0 < t := by linarith
  rw [integral_Kt_split ht (by norm_num), integral_Iic_Kt ht (by norm_num)]
  have h := integral_Ioi_Kt_ge ht (b := 3 / 2) (by norm_num)
  have e : ∫ u in (0 : ℝ)..log t, (t * exp (-(3 / 2 - 1 / 2) * u) - exp (-(3 / 2 - 3 / 2) * u))
      = t - 1 - log t := by
    have e1 : ∀ u : ℝ, t * exp (-(3 / 2 - 1 / 2) * u) - exp (-(3 / 2 - 3 / 2) * u)
        = t * exp (-1 * u) - 1 := by
      intro u
      rw [show -(3 / 2 - 1 / 2 : ℝ) * u = -1 * u by ring, show -(3 / 2 - 3 / 2 : ℝ) * u = 0 by ring,
        exp_zero]
    simp_rw [e1]
    rw [intervalIntegral.integral_sub (Continuous.intervalIntegrable (by fun_prop) _ _)
      (Continuous.intervalIntegrable (by fun_prop) _ _), intervalIntegral.integral_const_mul,
      integral_exp_neg_mul 1 (log t) one_ne_zero, intervalIntegral.integral_const]
    rw [show -1 * log t = -log t by ring, exp_neg, exp_log ht0]
    simp only [sub_zero, smul_eq_mul, mul_one, div_one]
    field_simp
  rw [e] at h
  have : t * (3 / (3 / 2 + 3 / 2) - 3 / (3 / 2 + 5 / 2) + 1 / (3 / 2 + 7 / 2))
      - (3 / (3 / 2 + 5 / 2) - 3 / (3 / 2 + 7 / 2) + 1 / (3 / 2 + 9 / 2))
      = 9 / 20 * t - 19 / 60 := by norm_num; ring
  rw [this]
  linarith

open Real in
/-- `I(2n+7/2)` from below. -/
theorem integral_Kt_ge {t : ℝ} (ht : 2 ≤ t) (n : ℕ) :
    t * (3 / (2 * (n : ℝ) + 5) - 3 / (2 * (n : ℝ) + 6) + 1 / (2 * (n : ℝ) + 7))
        - (3 / (2 * (n : ℝ) + 6) - 3 / (2 * (n : ℝ) + 7) + 1 / (2 * (n : ℝ) + 8))
        + t / (2 * (n : ℝ) + 3) - 1 / (((n : ℝ) + 2) * (2 * (n : ℝ) + 3)) - 1 / (2 * (n : ℝ) + 2)
      ≤ ∫ u : ℝ, Kt t u * exp (-(2 * (n : ℝ) + 7 / 2) * |u|) := by
  have ht0 : 0 < t := by linarith
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hb : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by linarith
  have hL : 1 / 2 ≤ log t := by
    have h2 : log 2 ≤ log t := log_le_log (by norm_num) ht
    have := Real.log_two_gt_d9
    linarith
  rw [integral_Kt_split ht hb.le, integral_Iic_Kt ht hb]
  have h := integral_Ioi_Kt_ge ht (b := 2 * (n : ℝ) + 7 / 2) hb.le
  have e : ∫ u in (0 : ℝ)..log t, (t * exp (-(2 * (n : ℝ) + 7 / 2 - 1 / 2) * u)
        - exp (-(2 * (n : ℝ) + 7 / 2 - 3 / 2) * u))
      = t * ((1 - exp (-(2 * (n : ℝ) + 3) * log t)) / (2 * (n : ℝ) + 3))
        - (1 - exp (-(2 * (n : ℝ) + 2) * log t)) / (2 * (n : ℝ) + 2) := by
    rw [show (2 * (n : ℝ) + 7 / 2 - 1 / 2) = 2 * (n : ℝ) + 3 by ring,
      show (2 * (n : ℝ) + 7 / 2 - 3 / 2) = 2 * (n : ℝ) + 2 by ring,
      intervalIntegral.integral_sub (Continuous.intervalIntegrable (by fun_prop) _ _)
      (Continuous.intervalIntegrable (by fun_prop) _ _), intervalIntegral.integral_const_mul,
      integral_exp_neg_mul _ (log t) (by linarith), integral_exp_neg_mul _ (log t) (by linarith)]
  rw [e] at h
  -- `t e^{−(2n+3) log t} = e^{−(2n+2) log t} ≤ 1/(n+2)`
  have hE : t * exp (-(2 * (n : ℝ) + 3) * log t) = exp (-(2 * (n : ℝ) + 2) * log t) := by
    nth_rewrite 1 [← exp_log ht0]
    rw [← exp_add]
    congr 1
    ring
  have hE1 : exp (-(2 * (n : ℝ) + 2) * log t) ≤ 1 / ((n : ℝ) + 2) := by
    have h1 := Real.add_one_le_exp ((2 * (n : ℝ) + 2) * log t)
    have h2 : (n : ℝ) + 2 ≤ exp ((2 * (n : ℝ) + 2) * log t) := by nlinarith
    rw [show -(2 * (n : ℝ) + 2) * log t = -((2 * (n : ℝ) + 2) * log t) by ring, exp_neg,
      one_div]
    exact inv_anti₀ (by linarith) h2
  have hE0 := exp_pos (-(2 * (n : ℝ) + 2) * log t)
  have h3 : t / (2 * (n : ℝ) + 3) - 1 / (((n : ℝ) + 2) * (2 * (n : ℝ) + 3)) - 1 / (2 * (n : ℝ) + 2)
      ≤ t * ((1 - exp (-(2 * (n : ℝ) + 3) * log t)) / (2 * (n : ℝ) + 3))
        - (1 - exp (-(2 * (n : ℝ) + 2) * log t)) / (2 * (n : ℝ) + 2) := by
    have a1 : t * ((1 - exp (-(2 * (n : ℝ) + 3) * log t)) / (2 * (n : ℝ) + 3))
        = t / (2 * (n : ℝ) + 3) - exp (-(2 * (n : ℝ) + 2) * log t) / (2 * (n : ℝ) + 3) := by
      rw [← hE]; ring
    have a2 : exp (-(2 * (n : ℝ) + 2) * log t) / (2 * (n : ℝ) + 3)
        ≤ 1 / (((n : ℝ) + 2) * (2 * (n : ℝ) + 3)) := by
      rw [← div_div]
      exact div_le_div_of_nonneg_right hE1 (by linarith)
    have a3 : (1 - exp (-(2 * (n : ℝ) + 2) * log t)) / (2 * (n : ℝ) + 2) ≤ 1 / (2 * (n : ℝ) + 2) :=
      div_le_div_of_nonneg_right (by linarith) (by linarith)
    rw [a1]
    linarith
  have e5 : (2 * (n : ℝ) + 7 / 2 + 3 / 2) = 2 * (n : ℝ) + 5 := by ring
  have e6 : (2 * (n : ℝ) + 7 / 2 + 5 / 2) = 2 * (n : ℝ) + 6 := by ring
  have e7 : (2 * (n : ℝ) + 7 / 2 + 7 / 2) = 2 * (n : ℝ) + 7 := by ring
  have e8 : (2 * (n : ℝ) + 7 / 2 + 9 / 2) = 2 * (n : ℝ) + 8 := by ring
  rw [e5, e6, e7, e8]
  linarith

open Real in
/-- The `n`-th term of the series in `arch_weil` for `Kt t`. -/
theorem arch_term_Kt_le {t : ℝ} (ht : 2 ≤ t) (n : ℕ) :
    Kt t 0 / ((n : ℝ) + 1) - ∫ u : ℝ, Kt t u * exp (-(2 * (n : ℝ) + 7 / 2) * |u|)
      ≤ (t + 1) * (1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2)) := by
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have h := integral_Kt_ge ht n
  rw [Kt_zero ht]
  -- rational inequalities
  have hc : 1 / ((n : ℝ) + 1) - 1 / (2 * (n : ℝ) + 3) - 3 / (2 * (n : ℝ) + 5)
      + 3 / (2 * (n : ℝ) + 6) - 1 / (2 * (n : ℝ) + 7) ≤ 1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2) := by
    rw [← sub_nonneg]
    have e : 1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2) - (1 / ((n : ℝ) + 1) - 1 / (2 * (n : ℝ) + 3)
        - 3 / (2 * (n : ℝ) + 5) + 3 / (2 * (n : ℝ) + 6) - 1 / (2 * (n : ℝ) + 7))
        = (4 * (n : ℝ) ^ 4 + 40 * (n : ℝ) ^ 3 + 137 * (n : ℝ) ^ 2 + 197 * (n : ℝ) + 96)
          / (((n : ℝ) + 1) * ((n : ℝ) + 2) * (2 * (n : ℝ) + 3) * (2 * (n : ℝ) + 5)
            * (2 * (n : ℝ) + 6) * (2 * (n : ℝ) + 7)) := by
      field_simp
      ring
    rw [e]
    positivity
  have hd : -1 / ((n : ℝ) + 1) + (3 / (2 * (n : ℝ) + 6) - 3 / (2 * (n : ℝ) + 7)
      + 1 / (2 * (n : ℝ) + 8)) + 1 / (2 * (n : ℝ) + 2) ≤ 0 := by
    have e : -1 / ((n : ℝ) + 1) + (3 / (2 * (n : ℝ) + 6) - 3 / (2 * (n : ℝ) + 7)
        + 1 / (2 * (n : ℝ) + 8)) + 1 / (2 * (n : ℝ) + 2)
        = -(12 * (n : ℝ) ^ 2 + 96 * (n : ℝ) + 204)
          / ((2 * (n : ℝ) + 2) * (2 * (n : ℝ) + 6) * (2 * (n : ℝ) + 7) * (2 * (n : ℝ) + 8)) := by
      field_simp
      ring
    rw [e]
    apply div_nonpos_of_nonpos_of_nonneg
    · nlinarith
    · positivity
  have hq : 1 / (((n : ℝ) + 2) * (2 * (n : ℝ) + 3)) ≤ 1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2) := by
    have e : 1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2) = 1 / (((n : ℝ) + 2) * ((n : ℝ) + 1)) := by
      field_simp
      ring
    rw [e]
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith
  have ht0 : (0 : ℝ) ≤ t := by linarith
  have key : (t - 1) / ((n : ℝ) + 1)
      - (t * (3 / (2 * (n : ℝ) + 5) - 3 / (2 * (n : ℝ) + 6) + 1 / (2 * (n : ℝ) + 7))
        - (3 / (2 * (n : ℝ) + 6) - 3 / (2 * (n : ℝ) + 7) + 1 / (2 * (n : ℝ) + 8))
        + t / (2 * (n : ℝ) + 3) - 1 / (((n : ℝ) + 2) * (2 * (n : ℝ) + 3)) - 1 / (2 * (n : ℝ) + 2))
      = t * (1 / ((n : ℝ) + 1) - 1 / (2 * (n : ℝ) + 3) - 3 / (2 * (n : ℝ) + 5)
          + 3 / (2 * (n : ℝ) + 6) - 1 / (2 * (n : ℝ) + 7))
        + (-1 / ((n : ℝ) + 1) + (3 / (2 * (n : ℝ) + 6) - 3 / (2 * (n : ℝ) + 7)
          + 1 / (2 * (n : ℝ) + 8)) + 1 / (2 * (n : ℝ) + 2))
        + 1 / (((n : ℝ) + 2) * (2 * (n : ℝ) + 3)) := by ring
  have h1 := mul_le_mul_of_nonneg_left hc ht0
  nlinarith

open ArchWeilAux in
/-- **The archimedean term of `Kt t` is at most `log t + 3`.** -/
theorem arch_Kt_le {t : ℝ} (ht : 2 ≤ t) :
    ∃ a : ℝ, a ≤ Real.log t + 3 ∧
      (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, Zeta23.paperFT (fun u => ((Kt t u : ℝ) : ℂ)) (r : ℂ)
        * ((W4 r : ℝ) : ℂ) = (a : ℂ) := by
  have hM : ∀ r : ℝ, ‖Zeta23.paperFT (fun u => ((Kt t u : ℝ) : ℂ)) (r : ℂ)‖ * (1 + r ^ 2)
      ≤ 3 * (t + 1) ^ (3 / 2 : ℝ) + 24 / 5 * (3 / 7 * t + 1 / 5) :=
    fun r => norm_integral_Kt_mul_cexp_mul_le ht r
  obtain ⟨hsum, -, heq⟩ := arch_weil (continuous_Kt t) (integrable_Kt (by linarith)) hM
  refine ⟨_, ?_, heq⟩
  -- the series
  have hser : ∑' n : ℕ, (Kt t 0 / ((n : ℝ) + 1)
      - ∫ u : ℝ, Kt t u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|)) ≤ t + 1 := by
    refine hsum.tsum_le_of_sum_range_le fun N => ?_
    have hpart : ∀ N : ℕ, ∑ n ∈ Finset.range N, (1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2))
        = 1 - 1 / ((N : ℝ) + 1) := by
      intro N
      induction N with
      | zero => simp
      | succ N ih =>
        rw [Finset.sum_range_succ, ih]
        push_cast
        ring
    calc ∑ n ∈ Finset.range N, (Kt t 0 / ((n : ℝ) + 1)
          - ∫ u : ℝ, Kt t u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|))
        ≤ ∑ n ∈ Finset.range N, (t + 1) * (1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 2)) :=
          Finset.sum_le_sum fun n _ => arch_term_Kt_le ht n
      _ = (t + 1) * (1 - 1 / ((N : ℝ) + 1)) := by rw [← Finset.mul_sum, hpart]
      _ ≤ t + 1 := by
          have : (0 : ℝ) ≤ 1 / ((N : ℝ) + 1) := by positivity
          nlinarith
  have h32 := integral_Kt_three_halves_ge ht
  rw [Kt_zero ht] at hser ⊢
  have hlog : Real.log (4 / Real.pi) ≤ 1 / 3 := by
    have h1 : Real.log (4 / Real.pi) ≤ Real.log (4 / 3) :=
      Real.log_le_log (by positivity) (by
        rw [div_le_div_iff₀ Real.pi_pos (by norm_num)]
        linarith [Real.pi_gt_three])
    have h2 : Real.log (4 / 3) ≤ 4 / 3 - 1 := Real.log_le_sub_one_of_pos (by norm_num)
    linarith
  have hγ := Real.one_half_lt_eulerMascheroniConstant
  have hneg : (t - 1) * (Real.log (4 / Real.pi) - Real.eulerMascheroniConstant) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (by linarith) (by linarith)
  linarith

end Allx
