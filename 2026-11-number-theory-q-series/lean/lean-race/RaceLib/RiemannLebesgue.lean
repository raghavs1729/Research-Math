/-
RaceLib/RiemannLebesgue.lean — the Riemann–Lebesgue lemma in the form
`∫ e^{irL} F(r) dr → 0` as `L → +∞` (Mathlib only; no integrability hypothesis is needed, because
Bochner integrals of non-integrable functions are `0`).
-/
import Mathlib

open Complex MeasureTheory Filter Topology
open scoped FourierTransform

noncomputable section

namespace RaceLib

theorem tendsto_integral_cexp_mul (F : ℝ → ℂ) :
    Tendsto (fun L : ℝ => ∫ r : ℝ, cexp (I * r * L) * F r) atTop (𝓝 0) := by
  have h := Real.tendsto_integral_exp_smul_cocompact F
  have hφ : Tendsto (fun L : ℝ => -L / (2 * Real.pi)) atTop (cocompact ℝ) := by
    have h1 : Tendsto (fun L : ℝ => -L / (2 * Real.pi)) atTop atBot :=
      tendsto_neg_atTop_atBot.atBot_div_const (by positivity)
    exact h1.mono_right atBot_le_cocompact
  refine (h.comp hφ).congr fun L => ?_
  simp only [Function.comp_apply]
  refine integral_congr_ae (Eventually.of_forall fun r => ?_)
  simp only [Circle.smul_def, Real.fourierChar_apply, smul_eq_mul]
  congr 2
  have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  push_cast
  field_simp

end RaceLib
