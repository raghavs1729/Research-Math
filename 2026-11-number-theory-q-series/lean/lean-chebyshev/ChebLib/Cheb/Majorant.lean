/-
ChebLib/Cheb/Majorant.lean — the majorant test function
`g(u) = (5/2) e^{−(3/2)|u|} − (3/2) e^{−(5/2)|u|}`.

It is `C²` (the kinks of the two exponentials cancel), even, `g(0) = 1`, and its Fourier transform
(convention `∫ g(u) e^{iru} du`) is the product of two Cauchy kernels
`30 / ((9/4 + r²)(25/4 + r²)) > 0`.
-/
import Mathlib

open Complex MeasureTheory Real

noncomputable section

namespace ChebLib

/-- `g(u) = (5/2) e^{−(3/2)|u|} − (3/2) e^{−(5/2)|u|}`. -/
def maj (u : ℝ) : ℝ := 5 / 2 * Real.exp (-(3 / 2) * |u|) - 3 / 2 * Real.exp (-(5 / 2) * |u|)

theorem maj_neg (u : ℝ) : maj (-u) = maj u := by
  simp [maj]

/-- `g(log n) = (5/2) n^{-3/2} − (3/2) n^{-5/2}` for `n ≥ 1`. -/
theorem maj_log {n : ℕ} (hn : 0 < n) :
    maj (Real.log n) = 5 / 2 * (n : ℝ) ^ (-(3 / 2) : ℝ) - 3 / 2 * (n : ℝ) ^ (-(5 / 2) : ℝ) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hlog : 0 ≤ Real.log n := Real.log_nonneg h1
  unfold maj
  rw [abs_of_nonneg hlog, Real.rpow_def_of_pos hn', Real.rpow_def_of_pos hn']
  congr 3 <;> ring

/-- Gluing two differentiable functions at `0` when their derivatives agree there. -/
private lemma hasDerivAt_glue {F A B A' B' : ℝ → ℝ} (hA : ∀ t, HasDerivAt A (A' t) t)
    (hB : ∀ t, HasDerivAt B (B' t) t) (hFA : ∀ u, 0 ≤ u → F u = A u)
    (hFB : ∀ u, u ≤ 0 → F u = B u) (h0 : A' 0 = B' 0) (u : ℝ) :
    HasDerivAt F (if 0 ≤ u then A' u else B' u) u := by
  rcases lt_trichotomy u 0 with hu | hu | hu
  · rw [if_neg (not_le.mpr hu)]
    refine (hB u).congr_of_eventuallyEq ?_
    filter_upwards [Iio_mem_nhds hu] with x hx using hFB x (le_of_lt hx)
  · subst hu
    rw [if_pos le_rfl]
    have h1 : HasDerivWithinAt F (A' 0) (Set.Ici 0) 0 :=
      (hA 0).hasDerivWithinAt.congr (fun x hx => hFA x hx) (hFA 0 le_rfl)
    have h2 : HasDerivWithinAt F (A' 0) (Set.Iic 0) 0 := by
      rw [h0]
      exact (hB 0).hasDerivWithinAt.congr (fun x hx => hFB x hx) (hFB 0 le_rfl)
    have h3 := h2.union h1
    rwa [Set.Iic_union_Ici, hasDerivWithinAt_univ] at h3
  · rw [if_pos hu.le]
    refine (hA u).congr_of_eventuallyEq ?_
    filter_upwards [Ioi_mem_nhds hu] with x hx using hFA x (le_of_lt hx)

private lemma hasDerivAt_exp_mul (c a t : ℝ) :
    HasDerivAt (fun x : ℝ => c * Real.exp (a * x)) (c * a * Real.exp (a * t)) t := by
  have h := (((hasDerivAt_id' t).const_mul a).exp).const_mul c
  have e : c * a * Real.exp (a * t) = c * (Real.exp (a * t) * (a * 1)) := by ring
  rw [e]
  exact h

/-- Integrability of the Fourier kernel `e^{-b|u|} e^{iru}`. -/
private lemma integrableOn_kernel_Ioi {b : ℝ} (hb : 0 < b) (r : ℝ) :
    IntegrableOn (fun u : ℝ => (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)) (Set.Ioi 0) := by
  have h := integrableOn_exp_mul_complex_Ioi (a := -(b : ℂ) + I * r) (by simp [hb]) 0
  refine h.congr_fun ?_ measurableSet_Ioi
  intro x hx
  have hx' : (0 : ℝ) < x := hx
  simp only [abs_of_pos hx']
  push_cast
  rw [← Complex.exp_add]
  congr 1
  ring

private lemma integrableOn_kernel_Iic {b : ℝ} (hb : 0 < b) (r : ℝ) :
    IntegrableOn (fun u : ℝ => (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)) (Set.Iic 0) := by
  have h := integrableOn_exp_mul_complex_Iic (a := (b : ℂ) + I * r) (by simp [hb]) 0
  refine h.congr_fun ?_ measurableSet_Iic
  intro x hx
  have hx' : x ≤ 0 := hx
  simp only [abs_of_nonpos hx']
  push_cast
  rw [← Complex.exp_add]
  congr 1
  ring

private lemma integrable_kernel {b : ℝ} (hb : 0 < b) (r : ℝ) :
    Integrable (fun u : ℝ => (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)) := by
  have h := (integrableOn_kernel_Iic hb r).union (integrableOn_kernel_Ioi hb r)
  rwa [Set.Iic_union_Ioi, integrableOn_univ] at h

private lemma integrable_exp_neg_mul_abs {b : ℝ} (hb : 0 < b) :
    Integrable (fun u : ℝ => Real.exp (-b * |u|)) := by
  have h := (integrable_kernel hb 0).norm
  simpa only [Complex.ofReal_zero, mul_zero, zero_mul, Complex.exp_zero, mul_one,
    Complex.norm_real, Real.norm_eq_abs, Real.abs_exp] using h

private lemma integral_kernel {b : ℝ} (hb : 0 < b) (r : ℝ) :
    ∫ u : ℝ, (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)
      = ((2 * b / (b ^ 2 + r ^ 2) : ℝ) : ℂ) := by
  rw [← intervalIntegral.integral_Iic_add_Ioi (integrableOn_kernel_Iic hb r)
    (integrableOn_kernel_Ioi hb r)]
  have h1 : ∫ u : ℝ in Set.Iic 0, (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)
      = ∫ u : ℝ in Set.Iic 0, cexp (((b : ℂ) + I * r) * u) := by
    refine setIntegral_congr_fun measurableSet_Iic ?_
    intro x hx
    have hx' : x ≤ 0 := hx
    simp only [abs_of_nonpos hx']
    push_cast
    rw [← Complex.exp_add]
    congr 1
    ring
  have h2 : ∫ u : ℝ in Set.Ioi 0, (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)
      = ∫ u : ℝ in Set.Ioi 0, cexp ((-(b : ℂ) + I * r) * u) := by
    refine setIntegral_congr_fun measurableSet_Ioi ?_
    intro x hx
    have hx' : (0 : ℝ) < x := hx
    simp only [abs_of_pos hx']
    push_cast
    rw [← Complex.exp_add]
    congr 1
    ring
  rw [h1, h2, integral_exp_mul_complex_Iic (by simp [hb]),
    integral_exp_mul_complex_Ioi (by simp [hb])]
  simp only [Complex.ofReal_zero, mul_zero, Complex.exp_zero]
  have hne1 : (b : ℂ) + I * r ≠ 0 := by
    intro h
    have := congrArg Complex.re h
    simp at this
    exact hb.ne' this
  have hne2 : (b : ℂ) - I * r ≠ 0 := by
    intro h
    have := congrArg Complex.re h
    simp at this
    exact hb.ne' this
  have hne3 : ((b : ℂ) ^ 2 + (r : ℂ) ^ 2) ≠ 0 := by
    exact_mod_cast (by positivity : (b ^ 2 + r ^ 2 : ℝ) ≠ 0)
  rw [show -(b : ℂ) + I * r = -((b : ℂ) - I * r) by ring, neg_div_neg_eq]
  push_cast
  rw [div_add_div _ _ hne1 hne2, div_eq_div_iff (mul_ne_zero hne1 hne2) hne3]
  linear_combination (2 * (b : ℂ) * (r : ℂ) ^ 2) * Complex.I_sq

/-- First derivative of `maj`. -/
private def maj1 (u : ℝ) : ℝ :=
  if 0 ≤ u then 5 / 2 * (-(3 / 2)) * Real.exp (-(3 / 2) * u) - 3 / 2 * (-(5 / 2)) * Real.exp (-(5 / 2) * u)
  else 5 / 2 * (3 / 2) * Real.exp (3 / 2 * u) - 3 / 2 * (5 / 2) * Real.exp (5 / 2 * u)

/-- Second derivative of `maj`. -/
private def maj2 (u : ℝ) : ℝ :=
  45 / 8 * Real.exp (-(3 / 2) * |u|) - 75 / 8 * Real.exp (-(5 / 2) * |u|)

private lemma hasDerivAt_maj (u : ℝ) : HasDerivAt maj (maj1 u) u := by
  refine hasDerivAt_glue
    (A := fun x => 5 / 2 * Real.exp (-(3 / 2) * x) - 3 / 2 * Real.exp (-(5 / 2) * x))
    (B := fun x => 5 / 2 * Real.exp (3 / 2 * x) - 3 / 2 * Real.exp (5 / 2 * x))
    (fun t => (hasDerivAt_exp_mul (5 / 2) (-(3 / 2)) t).sub (hasDerivAt_exp_mul (3 / 2) (-(5 / 2)) t))
    (fun t => (hasDerivAt_exp_mul (5 / 2) (3 / 2) t).sub (hasDerivAt_exp_mul (3 / 2) (5 / 2) t))
    ?_ ?_ ?_ u
  · intro x hx
    simp only [maj, abs_of_nonneg hx]
  · intro x hx
    simp only [maj, abs_of_nonpos hx]
    congr 3 <;> ring
  · norm_num

private lemma hasDerivAt_maj1 (u : ℝ) : HasDerivAt maj1 (maj2 u) u := by
  have h := hasDerivAt_glue (F := maj1)
    (A := fun x => 5 / 2 * (-(3 / 2)) * Real.exp (-(3 / 2) * x)
      - 3 / 2 * (-(5 / 2)) * Real.exp (-(5 / 2) * x))
    (B := fun x => 5 / 2 * (3 / 2) * Real.exp (3 / 2 * x) - 3 / 2 * (5 / 2) * Real.exp (5 / 2 * x))
    (fun t => (hasDerivAt_exp_mul (5 / 2 * (-(3 / 2))) (-(3 / 2)) t).sub
      (hasDerivAt_exp_mul (3 / 2 * (-(5 / 2))) (-(5 / 2)) t))
    (fun t => (hasDerivAt_exp_mul (5 / 2 * (3 / 2)) (3 / 2) t).sub
      (hasDerivAt_exp_mul (3 / 2 * (5 / 2)) (5 / 2) t))
    ?_ ?_ ?_ u
  · convert h using 1
    unfold maj2
    split_ifs with hu
    · rw [abs_of_nonneg hu]; ring
    · rw [abs_of_neg (not_le.mp hu)]
      rw [show -(3 / 2 : ℝ) * -u = 3 / 2 * u by ring, show -(5 / 2 : ℝ) * -u = 5 / 2 * u by ring]
      ring
  · intro x hx
    simp only [maj1, if_pos hx]
  · intro x hx
    rcases hx.eq_or_lt with rfl | hx'
    · simp [maj1]; norm_num
    · simp only [maj1, if_neg (not_le.mpr hx')]
  · norm_num

private lemma abs_maj1_le (u : ℝ) :
    |maj1 u| ≤ 15 / 4 * Real.exp (-(3 / 2) * |u|) + 15 / 4 * Real.exp (-(5 / 2) * |u|) := by
  unfold maj1
  split_ifs with hu
  · rw [abs_of_nonneg hu]
    have h1 := Real.exp_pos (-(3 / 2) * u)
    have h2 := Real.exp_pos (-(5 / 2) * u)
    rw [abs_le]; constructor <;> nlinarith
  · rw [abs_of_neg (not_le.mp hu)]
    rw [show -(3 / 2 : ℝ) * -u = 3 / 2 * u by ring, show -(5 / 2 : ℝ) * -u = 5 / 2 * u by ring]
    have h1 := Real.exp_pos (3 / 2 * u)
    have h2 := Real.exp_pos (5 / 2 * u)
    rw [abs_le]; constructor <;> nlinarith

private lemma integrable_maj : Integrable maj :=
  ((integrable_exp_neg_mul_abs (b := 3 / 2) (by norm_num)).const_mul (5 / 2)).sub
    ((integrable_exp_neg_mul_abs (b := 5 / 2) (by norm_num)).const_mul (3 / 2))

/-- `g` is twice differentiable with continuous second derivative, and `g, g', g''` are
integrable. -/
theorem maj_data : ∃ g1 g2 : ℝ → ℝ, (∀ u, HasDerivAt maj (g1 u) u) ∧
    (∀ u, HasDerivAt g1 (g2 u) u) ∧ Continuous g2 ∧ Integrable maj ∧ Integrable g1 ∧
    Integrable g2 := by
  refine ⟨maj1, maj2, hasDerivAt_maj, hasDerivAt_maj1, ?_, integrable_maj, ?_, ?_⟩
  · unfold maj2
    fun_prop
  · have hcont : Continuous maj1 :=
      continuous_iff_continuousAt.mpr fun u => (hasDerivAt_maj1 u).continuousAt
    refine Integrable.mono'
      (((integrable_exp_neg_mul_abs (b := 3 / 2) (by norm_num)).const_mul (15 / 4)).add
        ((integrable_exp_neg_mul_abs (b := 5 / 2) (by norm_num)).const_mul (15 / 4)))
      hcont.aestronglyMeasurable (Filter.Eventually.of_forall fun u => ?_)
    rw [Real.norm_eq_abs]
    exact abs_maj1_le u
  · exact ((integrable_exp_neg_mul_abs (b := 3 / 2) (by norm_num)).const_mul (45 / 8)).sub
      ((integrable_exp_neg_mul_abs (b := 5 / 2) (by norm_num)).const_mul (75 / 8))

/-- The Fourier transform of `g`: `∫ g(u) e^{iru} du = 30 / ((9/4 + r²)(25/4 + r²))`. -/
theorem integral_maj_mul_cexp (r : ℝ) :
    ∫ u : ℝ, (maj u : ℂ) * cexp (I * r * u)
      = ((30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2)) : ℝ) : ℂ) := by
  have e : ∀ u : ℝ, (maj u : ℂ) * cexp (I * r * u)
      = (5 / 2 : ℂ) * ((Real.exp (-(3 / 2) * |u|) : ℂ) * cexp (I * r * u))
        - (3 / 2 : ℂ) * ((Real.exp (-(5 / 2) * |u|) : ℂ) * cexp (I * r * u)) := by
    intro u
    simp only [maj]
    push_cast
    ring
  simp_rw [e]
  rw [integral_sub ((integrable_kernel (b := 3 / 2) (by norm_num) r).const_mul _)
    ((integrable_kernel (b := 5 / 2) (by norm_num) r).const_mul _),
    integral_const_mul, integral_const_mul, integral_kernel (by norm_num), integral_kernel (by norm_num)]
  have hR : 5 / 2 * (2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2)) - 3 / 2 * (2 * (5 / 2) / ((5 / 2) ^ 2 + r ^ 2))
      = 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2)) := by
    have h5 : (9 / 4 + r ^ 2 : ℝ) ≠ 0 := by positivity
    have h6 : (25 / 4 + r ^ 2 : ℝ) ≠ 0 := by positivity
    rw [show ((3 : ℝ) / 2) ^ 2 = 9 / 4 by norm_num, show ((5 : ℝ) / 2) ^ 2 = 25 / 4 by norm_num]
    field_simp
    ring
  rw [← hR]
  push_cast
  ring

end ChebLib
