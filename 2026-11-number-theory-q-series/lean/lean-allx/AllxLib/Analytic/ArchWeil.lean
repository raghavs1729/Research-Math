/-
AllxLib/Analytic/ArchWeil.lean — the archimedean term of the explicit formula of `L(s, χ₋₄)` in
Weil's form, for a general test function.

For `k : ℝ → ℝ` continuous and integrable whose Fourier transform `k̂ = paperFT k` satisfies
`|k̂(r)| (1 + r²) ≤ M`:

  `(1/2π) ∫ k̂(r) (Re ψ(3/4 + ir/2) + log(4/π)) dr
      = k(0) (log(4/π) − γ) − ∫ k(u) e^{−(3/2)|u|} du
        + Σ_n (k(0)/(n+1) − ∫ k(u) e^{−(2n+7/2)|u|} du)`.

Proof: the partial-fraction series `Re ψ(3/4 + ir/2) = −γ − P_{3/2}(r) + Σ_n (1/(n+1) − P_{2n+7/2}(r))`,
`P_b(r) = 2b/(b²+r²)` (`ChebLib.PoissonAux.re_digamma_arch`), Fourier inversion
`∫ k̂ = 2π k(0)` and Parseval `∫ k̂ P_b = 2π ∫ k(u) e^{−b|u|} du` (`P_b` is the Fourier transform of
`e^{−b|u|}`), and termwise integration.
-/
import ChebLib.Cheb.Poisson
import ChebLib.Cheb.EFExt

open Complex MeasureTheory Real Filter Topology
open scoped FourierTransform
open Zeta23

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace Allx

/-- The archimedean weight of `χ₋₄`: `W(r) = Re ψ(3/4 + ir/2) + log(4/π)`. -/
def W4 (r : ℝ) : ℝ := (Complex.digamma (3 / 4 + Complex.I * r / 2)).re + Real.log (4 / Real.pi)

namespace ArchWeilAux

open ChebLib.PoissonAux

/-! ### The Fourier transform of `e^{−b|u|}` (copied from `ChebLib/Cheb/Majorant.lean`) -/

lemma integrableOn_kernel_Ioi {b : ℝ} (hb : 0 < b) (r : ℝ) :
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

lemma integrableOn_kernel_Iic {b : ℝ} (hb : 0 < b) (r : ℝ) :
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

lemma integrable_kernel {b : ℝ} (hb : 0 < b) (r : ℝ) :
    Integrable (fun u : ℝ => (Real.exp (-b * |u|) : ℂ) * cexp (I * r * u)) := by
  have h := (integrableOn_kernel_Iic hb r).union (integrableOn_kernel_Ioi hb r)
  rwa [Set.Iic_union_Ioi, integrableOn_univ] at h

/-- `u ↦ e^{−b|u|}` is integrable (`b > 0`). -/
lemma integrable_exp_neg_mul_abs {b : ℝ} (hb : 0 < b) :
    Integrable (fun u : ℝ => Real.exp (-b * |u|)) := by
  have h := (integrable_kernel hb 0).norm
  simpa only [Complex.ofReal_zero, mul_zero, zero_mul, Complex.exp_zero, mul_one,
    Complex.norm_real, Real.norm_eq_abs, Real.abs_exp] using h

/-- `∫ e^{−b|u|} e^{iru} du = 2b/(b²+r²)`. -/
lemma integral_kernel {b : ℝ} (hb : 0 < b) (r : ℝ) :
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

/-! ### The Poisson kernel `2b/(b²+r²)` -/

lemma continuous_poisson {b : ℝ} (hb : 0 < b) :
    Continuous (fun r : ℝ => 2 * b / (b ^ 2 + r ^ 2)) :=
  continuous_const.div (by fun_prop) (fun r => by positivity)

lemma poisson_nonneg {b : ℝ} (hb : 0 < b) (r : ℝ) : 0 ≤ 2 * b / (b ^ 2 + r ^ 2) := by
  positivity

lemma poisson_le {b : ℝ} (hb : 0 < b) (r : ℝ) : 2 * b / (b ^ 2 + r ^ 2) ≤ 2 * b / b ^ 2 := by
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_nonneg hb.le (sq_nonneg r), sq_nonneg b]

/-! ### Fourier facts for a test function with `|k̂(r)| (1 + r²) ≤ M` -/

section Fourier

variable {k : ℝ → ℝ} {M : ℝ}

lemma continuous_paperFT (hki : Integrable k) :
    Continuous (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ)) := by
  unfold paperFT
  refine continuous_of_dominated (bound := fun u => ‖(k u : ℂ)‖) (fun r => ?_) (fun r => ?_)
    hki.ofReal.norm ?_
  · exact hki.ofReal.aestronglyMeasurable.mul (Continuous.aestronglyMeasurable (by fun_prop))
  · refine Eventually.of_forall fun u => ?_
    rw [norm_mul, Complex.norm_exp]
    simp
  · exact Eventually.of_forall fun u => by fun_prop

lemma norm_paperFT_le
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) (r : ℝ) :
    ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ ≤ M * (1 / (1 ^ 2 + r ^ 2)) := by
  have h : (0 : ℝ) < 1 ^ 2 + r ^ 2 := by positivity
  rw [mul_one_div, le_div_iff₀ h, one_pow]
  exact hM r

lemma nonneg_of_bound
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) : 0 ≤ M := by
  have h := hM 0
  have h0 : 0 ≤ ‖paperFT (fun u => (k u : ℂ)) ((0 : ℝ) : ℂ)‖ * (1 + (0 : ℝ) ^ 2) := by positivity
  linarith

/-- `k̂` is integrable. -/
lemma integrable_paperFT (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ)) :=
  Integrable.mono' ((integrable_inv_sq_add one_pos).const_mul M)
    (continuous_paperFT hki).aestronglyMeasurable
    (Eventually.of_forall fun r => norm_paperFT_le hM r)

/-- Mathlib's Fourier transform of `k` is integrable. -/
lemma integrable_fourier (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    Integrable (𝓕 (fun u : ℝ => (k u : ℂ))) := by
  have h := (integrable_paperFT hki hM).comp_mul_left' (R := -(2 * π))
    (neg_ne_zero.mpr (by positivity))
  refine h.congr (Eventually.of_forall fun w => ?_)
  simp only
  rw [EF.paperFT_ofReal_eq_fourier]
  field_simp

/-- **(F1)** Fourier inversion at `0`: `∫ k̂(r) dr = 2π k(0)`. -/
theorem integral_paperFT (hk : Continuous k) (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    ∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) = ((2 * π * k 0 : ℝ) : ℂ) := by
  have h := EF.paper_inversion (k := fun u => (k u : ℂ)) (continuous_ofReal.comp hk) hki.ofReal
    (integrable_fourier hki hM) 0
  simp only [Complex.ofReal_zero, mul_zero, Complex.exp_zero, mul_one] at h
  have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  push_cast
  rw [h]
  field_simp

/-- **(F2)** Parseval against the Poisson kernel:
`∫ k̂(r) · 2b/(b²+r²) dr = 2π ∫ k(u) e^{−b|u|} du`. -/
theorem integral_paperFT_mul_poisson (hk : Continuous k) (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M)
    {b : ℝ} (hb : 0 < b) :
    ∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((2 * b / (b ^ 2 + r ^ 2) : ℝ) : ℂ)
      = ((2 * π * ∫ u : ℝ, k u * Real.exp (-b * |u|) : ℝ) : ℂ) := by
  have hE : Integrable (fun u : ℝ => ((Real.exp (-b * |u|) : ℝ) : ℂ)) :=
    (integrable_exp_neg_mul_abs hb).ofReal
  have h := EF.integral_k_mul_weight (k := fun u => (k u : ℂ)) (continuous_ofReal.comp hk)
    hki.ofReal (integrable_fourier hki hM) hE
  have hin : ∀ τ : ℝ, ∫ u : ℝ, ((Real.exp (-b * |u|) : ℝ) : ℂ) * cexp (-I * τ * u)
      = ((2 * b / (b ^ 2 + τ ^ 2) : ℝ) : ℂ) := by
    intro τ
    have h1 := integral_kernel hb (-τ)
    rw [neg_sq] at h1
    rw [← h1]
    congr 1
    ext u
    push_cast
    ring_nf
  simp_rw [hin] at h
  have h2 : ∫ u : ℝ, (k u : ℂ) * ((Real.exp (-b * |u|) : ℝ) : ℂ)
      = ((∫ u : ℝ, k u * Real.exp (-b * |u|) : ℝ) : ℂ) := by
    rw [← integral_complex_ofReal]
    congr 1
    ext u
    push_cast
    rfl
  have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have h3 : ((2 * π : ℝ) : ℂ) * (1 / (2 * π) : ℂ) = 1 := by
    push_cast
    field_simp
  rw [Complex.ofReal_mul, ← h2, h, ← mul_assoc, h3, one_mul]

/-- `k̂ · g` is integrable for `g` continuous and bounded. -/
lemma integrable_paperFT_mul_bdd (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M)
    {g : ℝ → ℝ} (hg : Continuous g) {C : ℝ} (hC : ∀ r, |g r| ≤ C) :
    Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((g r : ℝ) : ℂ)) := by
  refine (integrable_paperFT hki hM).mul_bdd (c := C)
    (continuous_ofReal.comp hg).aestronglyMeasurable (Eventually.of_forall fun r => ?_)
  rw [Complex.norm_real, Real.norm_eq_abs]
  exact hC r

lemma integrable_paperFT_mul_poisson (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M)
    {b : ℝ} (hb : 0 < b) :
    Integrable (fun r : ℝ =>
      paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((2 * b / (b ^ 2 + r ^ 2) : ℝ) : ℂ)) :=
  integrable_paperFT_mul_bdd hki hM (continuous_poisson hb) (C := 2 * b / b ^ 2) fun r => by
    rw [abs_of_nonneg (poisson_nonneg hb r)]
    exact poisson_le hb r

end Fourier

/-! ### The series terms against `1/(1+r²)` -/

/-- `q n = 1/(n+1) − 2/(2n+9/2) = (1/π) ∫ psiTerm n r/(1+r²) dr`. -/
def q (n : ℕ) : ℝ := 1 / ((n : ℝ) + 1) - 2 / (2 * (n : ℝ) + 9 / 2)

lemma summable_q : Summable q := by
  have h1 : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ 2) := by
    have h := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
    refine h.congr fun n => ?_
    push_cast
    rfl
  refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_) (h1.mul_left (5 / 4))
  · unfold q
    rw [sub_nonneg, div_le_div_iff₀ (by positivity) (by positivity)]
    linarith
  · unfold q
    have hn : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
    rw [div_sub_div _ _ (by positivity) (by positivity), mul_one_div,
      div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [mul_nonneg hn hn, mul_nonneg hn (mul_nonneg hn hn)]

lemma inv_mul_psiTerm_eq (n : ℕ) (r : ℝ) :
    1 / (1 ^ 2 + r ^ 2) * psiTerm n r
      = 1 / ((n : ℝ) + 1) * (1 / (1 ^ 2 + r ^ 2))
        - 2 * (2 * (n : ℝ) + 7 / 2)
          * (1 / ((1 ^ 2 + r ^ 2) * ((2 * (n : ℝ) + 7 / 2) ^ 2 + r ^ 2))) := by
  rw [psiTerm_eq]
  have h1 : ((1 : ℝ) ^ 2 + r ^ 2) ≠ 0 := by positivity
  have h2 : ((2 * (n : ℝ) + 7 / 2) ^ 2 + r ^ 2) ≠ 0 := by positivity
  field_simp

lemma integrable_inv_mul_psiTerm (n : ℕ) :
    Integrable (fun r : ℝ => 1 / (1 ^ 2 + r ^ 2) * psiTerm n r) := by
  have hc : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by positivity
  have h := ((integrable_inv_sq_add one_pos).const_mul (1 / ((n : ℝ) + 1))).sub
    ((integrable_inv_mul_inv one_pos hc).const_mul (2 * (2 * (n : ℝ) + 7 / 2)))
  exact h.congr (Eventually.of_forall fun r => (inv_mul_psiTerm_eq n r).symm)

lemma integral_inv_mul_psiTerm (n : ℕ) :
    ∫ r : ℝ, 1 / (1 ^ 2 + r ^ 2) * psiTerm n r = π * q n := by
  have hc : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by positivity
  simp_rw [inv_mul_psiTerm_eq]
  rw [integral_sub ((integrable_inv_sq_add one_pos).const_mul _)
    ((integrable_inv_mul_inv one_pos hc).const_mul _), integral_const_mul, integral_const_mul,
    integral_inv_sq_add one_pos, integral_inv_mul_inv one_pos hc]
  unfold q
  have h3 : (2 * (n : ℝ) + 9 / 2) ≠ 0 := by positivity
  field_simp
  ring

/-- `r ↦ (1/(1+r²)) Σ_n psiTerm n r` is integrable. -/
lemma integrable_inv_mul_tsum_psiTerm :
    Integrable (fun r : ℝ => 1 / (1 ^ 2 + r ^ 2) * ∑' n : ℕ, psiTerm n r) := by
  have h := integrable_tsum_of_nonneg integrable_inv_mul_psiTerm
    (fun n r => mul_nonneg (by positivity) (psiTerm_nonneg n r))
    (fun r => (summable_psiTerm r).mul_left _)
    (by simp_rw [integral_inv_mul_psiTerm]; exact summable_q.mul_left π)
  refine h.congr (Eventually.of_forall fun r => ?_)
  simp only
  rw [tsum_mul_left]

/-! ### Termwise integration -/

section Main

variable {k : ℝ → ℝ} {M : ℝ}

/-- The `n`-th term of the series in `arch_weil`. -/
def aTerm (k : ℝ → ℝ) (n : ℕ) : ℝ :=
  k 0 / ((n : ℝ) + 1) - ∫ u : ℝ, k u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|)

lemma paperFT_mul_psiTerm_eq (n : ℕ) (r : ℝ) :
    paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ)
      = paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((1 / ((n : ℝ) + 1) : ℝ) : ℂ)
        - paperFT (fun u => (k u : ℂ)) (r : ℂ)
          * ((2 * (2 * (n : ℝ) + 7 / 2) / ((2 * (n : ℝ) + 7 / 2) ^ 2 + r ^ 2) : ℝ) : ℂ) := by
  rw [psiTerm_eq, ← mul_sub, ← Complex.ofReal_sub]

lemma integrable_paperFT_mul_psiTerm (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) (n : ℕ) :
    Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ)) := by
  have hc : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by positivity
  have h := ((integrable_paperFT hki hM).mul_const (((1 / ((n : ℝ) + 1) : ℝ)) : ℂ)).sub
    (integrable_paperFT_mul_poisson hki hM hc)
  exact h.congr (Eventually.of_forall fun r => (paperFT_mul_psiTerm_eq n r).symm)

lemma integral_paperFT_mul_psiTerm (hk : Continuous k) (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) (n : ℕ) :
    ∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ)
      = ((2 * π * aTerm k n : ℝ) : ℂ) := by
  have hc : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by positivity
  simp_rw [paperFT_mul_psiTerm_eq]
  rw [integral_sub ((integrable_paperFT hki hM).mul_const _)
    (integrable_paperFT_mul_poisson hki hM hc), EF.cintegral_mul_const,
    integral_paperFT hk hki hM, integral_paperFT_mul_poisson hk hki hM hc]
  unfold aTerm
  push_cast
  ring

lemma summable_integral_norm (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    Summable (fun n : ℕ =>
      ∫ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ)‖) := by
  have hM0 : 0 ≤ M := nonneg_of_bound hM
  refine Summable.of_nonneg_of_le (fun n => integral_nonneg fun r => norm_nonneg _) (fun n => ?_)
    (summable_q.mul_left (M * π))
  have hb : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ)‖
      ≤ M * (1 / (1 ^ 2 + r ^ 2) * psiTerm n r) := by
    intro r
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (psiTerm_nonneg n r),
      ← mul_assoc]
    exact mul_le_mul_of_nonneg_right (norm_paperFT_le hM r) (psiTerm_nonneg n r)
  calc ∫ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ)‖
      ≤ ∫ r : ℝ, M * (1 / (1 ^ 2 + r ^ 2) * psiTerm n r) :=
        integral_mono (integrable_paperFT_mul_psiTerm hki hM n).norm
          ((integrable_inv_mul_psiTerm n).const_mul M) hb
    _ = M * π * q n := by
        rw [integral_const_mul, integral_inv_mul_psiTerm, mul_assoc]

/-- `k̂ · Σ_n psiTerm n` is integrable. -/
lemma integrable_paperFT_mul_tsum (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    Integrable (fun r : ℝ =>
      paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ)) := by
  have hcont : Continuous (fun r : ℝ =>
      paperFT (fun u => (k u : ℂ)) (r : ℂ) * (((1 + r ^ 2 : ℝ)) : ℂ)) :=
    (continuous_paperFT hki).mul (continuous_ofReal.comp (by fun_prop))
  have hG : Integrable (fun r : ℝ =>
      (((1 / (1 ^ 2 + r ^ 2) * ∑' n : ℕ, psiTerm n r : ℝ)) : ℂ)) :=
    integrable_inv_mul_tsum_psiTerm.ofReal
  have h := hG.bdd_mul (c := M) hcont.aestronglyMeasurable
    (Eventually.of_forall fun r => by
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
      exact hM r)
  refine h.congr (Eventually.of_forall fun r => ?_)
  simp only
  have h2 : ((1 + r ^ 2 : ℝ) : ℂ) * ((1 / (1 ^ 2 + r ^ 2) : ℝ) : ℂ) = 1 := by
    rw [← Complex.ofReal_mul, one_pow, mul_one_div_cancel (ne_of_gt (by positivity)),
      Complex.ofReal_one]
  rw [Complex.ofReal_mul]
  calc paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((1 + r ^ 2 : ℝ) : ℂ)
        * (((1 / (1 ^ 2 + r ^ 2) : ℝ) : ℂ) * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ))
      = paperFT (fun u => (k u : ℂ)) (r : ℂ)
        * (((1 + r ^ 2 : ℝ) : ℂ) * ((1 / (1 ^ 2 + r ^ 2) : ℝ) : ℂ))
        * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ) := by ring
    _ = _ := by rw [h2, mul_one]

/-- Termwise integration of `k̂ · Σ_n psiTerm n`. -/
lemma hasSum_aTerm (hk : Continuous k) (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    HasSum (fun n : ℕ => ((2 * π * aTerm k n : ℝ) : ℂ))
      (∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ)) := by
  have h := hasSum_integral_of_summable_integral_norm (integrable_paperFT_mul_psiTerm hki hM)
    (summable_integral_norm hki hM)
  simp_rw [integral_paperFT_mul_psiTerm hk hki hM] at h
  have e : ∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ)
      = ∫ r : ℝ, ∑' n : ℕ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((psiTerm n r : ℝ) : ℂ) := by
    congr 1
    ext r
    rw [Complex.ofReal_tsum, tsum_mul_left]
  rw [e]
  exact h

lemma summable_aTerm (hk : Continuous k) (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    Summable (aTerm k) := by
  have h := Complex.summable_ofReal.mp (hasSum_aTerm hk hki hM).summable
  exact (summable_mul_left_iff (by positivity : (2 * π : ℝ) ≠ 0)).mp h

lemma integral_paperFT_mul_tsum (hk : Continuous k) (hki : Integrable k)
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    ∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ)
      = ((2 * π * ∑' n : ℕ, aTerm k n : ℝ) : ℂ) := by
  rw [← (hasSum_aTerm hk hki hM).tsum_eq, ← tsum_mul_left, Complex.ofReal_tsum]

lemma paperFT_mul_W4_eq (r : ℝ) :
    paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((W4 r : ℝ) : ℂ)
      = paperFT (fun u => (k u : ℂ)) (r : ℂ)
          * ((Real.log (4 / π) - Real.eulerMascheroniConstant : ℝ) : ℂ)
        - paperFT (fun u => (k u : ℂ)) (r : ℂ)
          * ((2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2) : ℝ) : ℂ)
        + paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((∑' n : ℕ, psiTerm n r : ℝ) : ℂ) := by
  unfold W4
  rw [re_digamma_arch, ← mul_sub, ← mul_add, ← Complex.ofReal_sub, ← Complex.ofReal_add]
  congr 2
  ring

end Main

end ArchWeilAux

open ArchWeilAux ChebLib.PoissonAux in
/-- **The archimedean term in Weil's form.** -/
theorem arch_weil {k : ℝ → ℝ} (hk : Continuous k) (hki : Integrable k) {M : ℝ}
    (hM : ∀ r : ℝ, ‖paperFT (fun u => (k u : ℂ)) (r : ℂ)‖ * (1 + r ^ 2) ≤ M) :
    Summable (fun n : ℕ => k 0 / ((n : ℝ) + 1)
        - ∫ u : ℝ, k u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|)) ∧
    Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((W4 r : ℝ) : ℂ)) ∧
    (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, paperFT (fun u => (k u : ℂ)) (r : ℂ) * ((W4 r : ℝ) : ℂ)
      = ((k 0 * (Real.log (4 / Real.pi) - Real.eulerMascheroniConstant)
          - (∫ u : ℝ, k u * Real.exp (-(3 / 2) * |u|))
          + ∑' n : ℕ, (k 0 / ((n : ℝ) + 1)
              - ∫ u : ℝ, k u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|)) : ℝ) : ℂ) := by
  have h32 : (0 : ℝ) < 3 / 2 := by norm_num
  have h1 : Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ)
      * ((Real.log (4 / π) - Real.eulerMascheroniConstant : ℝ) : ℂ)) :=
    (integrable_paperFT hki hM).mul_const _
  have h2 : Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ)
      * ((2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2) : ℝ) : ℂ)) :=
    integrable_paperFT_mul_poisson hki hM h32
  have h3 := integrable_paperFT_mul_tsum hki hM
  have h12 : Integrable (fun r : ℝ => paperFT (fun u => (k u : ℂ)) (r : ℂ)
        * ((Real.log (4 / π) - Real.eulerMascheroniConstant : ℝ) : ℂ)
      - paperFT (fun u => (k u : ℂ)) (r : ℂ)
        * ((2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2) : ℝ) : ℂ)) := h1.sub h2
  refine ⟨summable_aTerm hk hki hM, ?_, ?_⟩
  · exact (h12.add h3).congr (Eventually.of_forall fun r => (paperFT_mul_W4_eq r).symm)
  · simp_rw [paperFT_mul_W4_eq]
    rw [integral_add h12 h3, integral_sub h1 h2, EF.cintegral_mul_const,
      integral_paperFT hk hki hM, integral_paperFT_mul_poisson hk hki hM h32,
      integral_paperFT_mul_tsum hk hki hM]
    have hπ : (π : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
    have ha : (fun n : ℕ => k 0 / ((n : ℝ) + 1)
        - ∫ u : ℝ, k u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|)) = aTerm k := rfl
    rw [ha]
    push_cast
    field_simp

end Allx
