/-
ChebLib/Cheb/Kernel.lean — the Abel kernel `h_x(u) = e^{u/2} exp(−e^u/x)` in logarithmic variable.

`h_x(log n) / √n = e^{−n/x}`, and its Fourier transform (convention `∫ h(u) e^{iru} du`) is the
Mellin transform of `e^{−t/x}` on the critical line:  `x^{1/2+ir} Γ(1/2+ir)`.
-/
import Mathlib

open Complex MeasureTheory Real

noncomputable section

namespace ChebLib

/-- `h_x(u) = e^{u/2} exp(−e^u/x)`. -/
def abelKernel (x u : ℝ) : ℝ := Real.exp (u / 2) * Real.exp (-(Real.exp u / x))

/-- The first derivative of `h_x`. -/
def abelKernel1 (x u : ℝ) : ℝ := (1 / 2 - Real.exp u / x) * abelKernel x u

/-- The second derivative of `h_x`. -/
def abelKernel2 (x u : ℝ) : ℝ :=
  ((1 / 2 - Real.exp u / x) ^ 2 - Real.exp u / x) * abelKernel x u

theorem hasDerivAt_abelKernel (x u : ℝ) : HasDerivAt (abelKernel x) (abelKernel1 x u) u := by
  have h1 : HasDerivAt (fun u : ℝ => Real.exp (u / 2)) (Real.exp (u / 2) * (1 / 2)) u :=
    ((hasDerivAt_id u).div_const 2).exp
  have h2 : HasDerivAt (fun u : ℝ => Real.exp (-(Real.exp u / x)))
      (Real.exp (-(Real.exp u / x)) * (-(Real.exp u / x))) u :=
    (((Real.hasDerivAt_exp u).div_const x).neg).exp
  have h := h1.mul h2
  refine h.congr_deriv ?_
  simp only [abelKernel1, abelKernel]
  ring

theorem hasDerivAt_abelKernel1 (x u : ℝ) : HasDerivAt (abelKernel1 x) (abelKernel2 x u) u := by
  have h1 : HasDerivAt (fun u : ℝ => 1 / 2 - Real.exp u / x) (-(Real.exp u / x)) u :=
    ((Real.hasDerivAt_exp u).div_const x).const_sub (1 / 2)
  have h := h1.mul (hasDerivAt_abelKernel x u)
  refine h.congr_deriv ?_
  simp only [abelKernel2, abelKernel1]
  ring

theorem continuous_abelKernel2 (x : ℝ) : Continuous (abelKernel2 x) := by
  show Continuous (fun u : ℝ => ((1 / 2 - Real.exp u / x) ^ 2 - Real.exp u / x) *
    (Real.exp (u / 2) * Real.exp (-(Real.exp u / x))))
  fun_prop

/-- `e^{au} exp(−e^u/x)` is integrable on `ℝ` for `a, x > 0` (substitute `t = e^u`). -/
theorem integrable_exp_mul_exp_neg_exp {x a : ℝ} (hx : 0 < x) (ha : 0 < a) :
    Integrable (fun u : ℝ => Real.exp (a * u) * Real.exp (-(Real.exp u / x))) := by
  have h0 : IntegrableOn (fun t : ℝ => t ^ (a - 1) * Real.exp (-(1 / x) * t ^ (1 : ℝ)))
      (Set.Ioi 0) :=
    integrableOn_rpow_mul_exp_neg_mul_rpow (by linarith) one_pos (by positivity)
  have h1 := (integrableOn_image_iff_integrableOn_abs_deriv_smul (s := Set.univ)
    (f := Real.exp) (f' := Real.exp) MeasurableSet.univ
    (fun u _ => (Real.hasDerivAt_exp u).hasDerivWithinAt) Real.exp_injective.injOn
    (fun t : ℝ => t ^ (a - 1) * Real.exp (-(1 / x) * t ^ (1 : ℝ)))).1
    (by rwa [Set.image_univ, Real.range_exp])
  rw [integrableOn_univ] at h1
  refine h1.congr (Filter.Eventually.of_forall fun u => ?_)
  simp only [smul_eq_mul, abs_of_pos (Real.exp_pos u), Real.rpow_one]
  rw [← Real.exp_mul, ← mul_assoc, ← Real.exp_add]
  have e1 : u + u * (a - 1) = a * u := by ring
  have e2 : -(1 / x) * Real.exp u = -(Real.exp u / x) := by ring
  rw [e1, e2]

theorem integrable_abelKernel {x : ℝ} (hx : 0 < x) : Integrable (abelKernel x) := by
  have h := integrable_exp_mul_exp_neg_exp hx (a := 1 / 2) (by norm_num)
  refine h.congr (Filter.Eventually.of_forall fun u => ?_)
  simp only [abelKernel]
  rw [show (1 / 2 : ℝ) * u = u / 2 by ring]

theorem integrable_abelKernel1 {x : ℝ} (hx : 0 < x) : Integrable (abelKernel1 x) := by
  have h1 := integrable_exp_mul_exp_neg_exp hx (a := 1 / 2) (by norm_num)
  have h2 := integrable_exp_mul_exp_neg_exp hx (a := 3 / 2) (by norm_num)
  refine ((h1.const_mul (1 / 2)).sub (h2.const_mul (1 / x))).congr
    (Filter.Eventually.of_forall fun u => ?_)
  simp only [abelKernel1, abelKernel, Pi.sub_apply]
  have e : Real.exp (3 / 2 * u) = Real.exp u * Real.exp (u / 2) := by
    rw [← Real.exp_add]; congr 1; ring
  rw [e, show (1 / 2 : ℝ) * u = u / 2 by ring]
  ring

theorem integrable_abelKernel2 {x : ℝ} (hx : 0 < x) : Integrable (abelKernel2 x) := by
  have h1 := integrable_exp_mul_exp_neg_exp hx (a := 1 / 2) (by norm_num)
  have h2 := integrable_exp_mul_exp_neg_exp hx (a := 3 / 2) (by norm_num)
  have h3 := integrable_exp_mul_exp_neg_exp hx (a := 5 / 2) (by norm_num)
  refine (((h1.const_mul (1 / 4)).sub (h2.const_mul (2 / x))).add
    (h3.const_mul (1 / x ^ 2))).congr (Filter.Eventually.of_forall fun u => ?_)
  simp only [abelKernel2, abelKernel, Pi.sub_apply, Pi.add_apply]
  have e : Real.exp (3 / 2 * u) = Real.exp u * Real.exp (u / 2) := by
    rw [← Real.exp_add]; congr 1; ring
  have e' : Real.exp (5 / 2 * u) = Real.exp u * Real.exp u * Real.exp (u / 2) := by
    rw [← Real.exp_add, ← Real.exp_add]; congr 1; ring
  rw [e, e', show (1 / 2 : ℝ) * u = u / 2 by ring]
  ring

/-- `h_x(log n) = √n e^{−n/x}`. -/
theorem abelKernel_log (x : ℝ) {n : ℕ} (hn : 0 < n) :
    abelKernel x (Real.log n) = Real.sqrt n * Real.exp (-(n / x)) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [abelKernel]
  rw [Real.exp_half, Real.exp_log hn']

/-- `h_x(−log n) = n^{-1/2} e^{−1/(n x)}`. -/
theorem abelKernel_neg_log (x : ℝ) {n : ℕ} (hn : 0 < n) :
    abelKernel x (-Real.log n) = (Real.sqrt n)⁻¹ * Real.exp (-(1 / (n * x))) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  simp only [abelKernel]
  rw [Real.exp_half, Real.exp_neg, Real.exp_log hn', Real.sqrt_inv]
  have e : (n : ℝ)⁻¹ / x = 1 / (n * x) := by
    rw [one_div, mul_inv, div_eq_mul_inv]
  rw [e]

/-- The Fourier transform of the Abel kernel: `∫ h_x(u) e^{iru} du = x^{1/2+ir} Γ(1/2+ir)`. -/
theorem integral_abelKernel_mul_cexp {x : ℝ} (hx : 0 < x) (r : ℝ) :
    ∫ u : ℝ, (abelKernel x u : ℂ) * cexp (I * r * u)
      = (x : ℂ) ^ ((1 / 2 : ℂ) + I * r) * Complex.Gamma (1 / 2 + I * r) := by
  have hare : 0 < ((1 / 2 : ℂ) + I * r).re := by simp
  have hx' : (0 : ℝ) < 1 / x := by positivity
  have key := integral_cpow_mul_exp_neg_mul_Ioi hare hx'
  have cv := integral_image_eq_integral_abs_deriv_smul (s := Set.univ)
    (f := Real.exp) (f' := Real.exp) MeasurableSet.univ
    (fun u _ => (Real.hasDerivAt_exp u).hasDerivWithinAt) Real.exp_injective.injOn
    (fun t : ℝ => (t : ℂ) ^ ((1 / 2 : ℂ) + I * r - 1) * cexp (-(((1 / x : ℝ) : ℂ) * t)))
  rw [Set.image_univ, Real.range_exp, Measure.restrict_univ] at cv
  rw [key] at cv
  have hxx : (1 / ((1 / x : ℝ) : ℂ)) = (x : ℂ) := by
    push_cast; simp
  rw [hxx] at cv
  rw [cv]
  refine integral_congr_ae (Filter.Eventually.of_forall fun u => ?_)
  simp only [abelKernel, abs_of_pos (Real.exp_pos u)]
  have hcp : ((Real.exp u : ℝ) : ℂ) ^ ((1 / 2 : ℂ) + I * r - 1)
      = cexp (u * ((1 / 2 : ℂ) + I * r - 1)) := by
    rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast (Real.exp_pos u).ne'),
      ← Complex.ofReal_log (Real.exp_pos u).le, Real.log_exp]
  rw [hcp, Complex.real_smul]
  push_cast
  simp only [← Complex.exp_add]
  congr 1
  ring

/-- `|Γ(1/2 + ir)|² = π / cosh(π r)`. -/
theorem norm_Gamma_half_add_I_mul_sq (r : ℝ) :
    ‖Complex.Gamma (1 / 2 + I * r)‖ ^ 2 = π / Real.cosh (π * r) := by
  have h := Complex.Gamma_mul_Gamma_one_sub (1 / 2 + I * r)
  have hc : (1 : ℂ) - (1 / 2 + I * r) = (starRingEnd ℂ) (1 / 2 + I * r) := by
    apply Complex.ext
    · simp; norm_num
    · simp
  rw [hc, Complex.Gamma_conj, Complex.mul_conj'] at h
  have hs : Complex.sin (π * (1 / 2 + I * r)) = Complex.cosh (π * r) := by
    rw [show (π : ℂ) * (1 / 2 + I * r) = (π * r) * I + π / 2 by ring,
      Complex.sin_add_pi_div_two, Complex.cos_mul_I]
  rw [hs] at h
  exact_mod_cast h

/-- Partial sum lower bound for `cosh`. -/
theorem cosh_ge_poly8 (y : ℝ) :
    1 + y ^ 2 / 2 + y ^ 4 / 24 + y ^ 6 / 720 + y ^ 8 / 40320 ≤ Real.cosh y := by
  have h := sum_le_hasSum (Finset.range 5)
    (fun n _ => div_nonneg (by rw [pow_mul]; positivity) (Nat.cast_nonneg _))
    (Real.hasSum_cosh y)
  norm_num [Finset.sum_range_succ, Nat.factorial] at h
  linarith

/-- The Gamma factor on the critical line is dominated by a product of two Cauchy kernels:
`(9/4 + r²)(25/4 + r²) |Γ(1/2 + ir)| ≤ (225/16) √π`  (equality at `r = 0`). -/
theorem norm_Gamma_half_mul_le (r : ℝ) :
    (9 / 4 + r ^ 2) * (25 / 4 + r ^ 2) * ‖Complex.Gamma (1 / 2 + I * r)‖
      ≤ 225 / 16 * Real.sqrt π := by
  have hG2 := norm_Gamma_half_add_I_mul_sq r
  have hcosh : 0 < Real.cosh (π * r) := Real.cosh_pos _
  have hpoly : ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2)) ^ 2
      ≤ (225 / 16) ^ 2 * Real.cosh (π * r) := by
    have h1 := cosh_ge_poly8 (π * r)
    have hs : 0 ≤ r ^ 2 := sq_nonneg r
    have hpi : 9 ≤ π ^ 2 := by nlinarith [Real.pi_gt_three]
    have hY : 9 * r ^ 2 ≤ (π * r) ^ 2 := by rw [mul_pow]; nlinarith
    have hY0 : 0 ≤ 9 * r ^ 2 := by positivity
    have hY2 : (9 * r ^ 2) ^ 2 ≤ (π * r) ^ 4 := by
      calc (9 * r ^ 2) ^ 2 ≤ ((π * r) ^ 2) ^ 2 := pow_le_pow_left₀ hY0 hY 2
        _ = (π * r) ^ 4 := by ring
    have hY3 : (9 * r ^ 2) ^ 3 ≤ (π * r) ^ 6 := by
      calc (9 * r ^ 2) ^ 3 ≤ ((π * r) ^ 2) ^ 3 := pow_le_pow_left₀ hY0 hY 3
        _ = (π * r) ^ 6 := by ring
    have hY4 : (9 * r ^ 2) ^ 4 ≤ (π * r) ^ 8 := by
      calc (9 * r ^ 2) ^ 4 ≤ ((π * r) ^ 2) ^ 4 := pow_le_pow_left₀ hY0 hY 4
        _ = (π * r) ^ 8 := by ring
    have hp : ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2)) ^ 2 ≤ (225 / 16) ^ 2 *
        (1 + 9 * r ^ 2 / 2 + (9 * r ^ 2) ^ 2 / 24 + (9 * r ^ 2) ^ 3 / 720
          + (9 * r ^ 2) ^ 4 / 40320) := by
      nlinarith [hs, pow_nonneg hs 2, pow_nonneg hs 3, pow_nonneg hs 4]
    linarith
  apply le_of_pow_le_pow_left₀ two_ne_zero (by positivity)
  rw [mul_pow, mul_pow (225 / 16 : ℝ), Real.sq_sqrt Real.pi_pos.le, hG2, mul_div_assoc',
    div_le_iff₀ hcosh]
  nlinarith [mul_le_mul_of_nonneg_right hpoly Real.pi_pos.le]

end ChebLib
