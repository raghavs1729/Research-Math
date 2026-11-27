/-
ChebLib/Cheb/Poisson.lean — the archimedean term of the explicit formula for the majorant, in closed
form.

For `b > 0`:  `(1/2π) ∫ 2b/(b² + r²) · Re ψ(3/4 + ir/2) dr = ψ(3/4 + b/2)`  (Poisson integral of the
harmonic function `Re ψ`), proved termwise from the partial-fraction series
`Re ψ(a + it) = −γ − a/(a²+t²) + Σ_n (1/(n+1) − (n+1+a)/((n+1+a)²+t²))`
(`Zeta23.MuFields.re_digamma_vertical`) and `∫ dr/((b²+r²)(c²+r²)) = π/(bc(b+c))`.
With `30/((9/4+r²)(25/4+r²)) = (5/2)·3/(9/4+r²) − (3/2)·5/(25/4+r²)`,
`ψ(3/2) = 2 − γ − 2 log 2`, `ψ(2) = 1 − γ`:

  `(1/2π) ∫ 30/((9/4+r²)(25/4+r²)) (Re ψ(3/4 + ir/2) + log(4/π)) dr = 7/2 − γ − 3 log 2 − log π`.
-/
import Mathlib
import Zeta23.GammaFacts.Mu

open Complex MeasureTheory Real Filter Topology

noncomputable section

namespace ChebLib

namespace PoissonAux

/-! ### Elementary integrals -/

lemma integrable_inv_sq_add {b : ℝ} (hb : 0 < b) :
    Integrable (fun r : ℝ => 1 / (b ^ 2 + r ^ 2)) := by
  have h := (integrable_inv_one_add_sq.comp_div hb.ne').const_mul (b ^ 2)⁻¹
  refine h.congr (Eventually.of_forall fun r => ?_)
  simp only
  field_simp

lemma integral_inv_sq_add {b : ℝ} (hb : 0 < b) :
    ∫ r : ℝ, 1 / (b ^ 2 + r ^ 2) = π / b := by
  have h : (fun r : ℝ => 1 / (b ^ 2 + r ^ 2))
      = fun r : ℝ => (b ^ 2)⁻¹ * (fun x : ℝ => (1 + x ^ 2)⁻¹) (r / b) := by
    funext r
    simp only
    field_simp
  rw [h, integral_const_mul, Measure.integral_comp_div (fun x : ℝ => (1 + x ^ 2)⁻¹) b,
    integral_univ_inv_one_add_sq, abs_of_pos hb, smul_eq_mul]
  field_simp

lemma integrable_inv_mul_inv {b c : ℝ} (hb : 0 < b) (hc : 0 < c) :
    Integrable (fun r : ℝ => 1 / ((b ^ 2 + r ^ 2) * (c ^ 2 + r ^ 2))) := by
  refine Integrable.mono' ((integrable_inv_sq_add hb).const_mul (1 / c ^ 2)) ?_
    (Eventually.of_forall fun r => ?_)
  · refine Continuous.aestronglyMeasurable ?_
    exact continuous_const.div (by fun_prop) (fun r => by positivity)
  · rw [Real.norm_eq_abs, abs_of_pos (by positivity)]
    rw [div_mul_div_comm, one_mul, mul_comm (c ^ 2)]
    apply one_div_le_one_div_of_le (by positivity)
    have : 0 < b ^ 2 + r ^ 2 := by positivity
    nlinarith [sq_nonneg r]

lemma integral_inv_mul_inv_of_ne {b c : ℝ} (hb : 0 < b) (hc : 0 < c) (hbc : b ≠ c) :
    ∫ r : ℝ, 1 / ((b ^ 2 + r ^ 2) * (c ^ 2 + r ^ 2)) = π / (b * c * (b + c)) := by
  have hne : c - b ≠ 0 := sub_ne_zero.mpr (Ne.symm hbc)
  have hne2 : c ^ 2 - b ^ 2 ≠ 0 := by
    have : c ^ 2 - b ^ 2 = (c - b) * (c + b) := by ring
    rw [this]
    exact mul_ne_zero hne (by positivity)
  have h : (fun r : ℝ => 1 / ((b ^ 2 + r ^ 2) * (c ^ 2 + r ^ 2)))
      = fun r : ℝ => (1 / (c ^ 2 - b ^ 2)) * (1 / (b ^ 2 + r ^ 2) - 1 / (c ^ 2 + r ^ 2)) := by
    funext r
    have h1 : b ^ 2 + r ^ 2 ≠ 0 := by positivity
    have h2 : c ^ 2 + r ^ 2 ≠ 0 := by positivity
    field_simp
    ring
  rw [h, integral_const_mul, integral_sub (integrable_inv_sq_add hb) (integrable_inv_sq_add hc),
    integral_inv_sq_add hb, integral_inv_sq_add hc]
  have : c ^ 2 - b ^ 2 = (c - b) * (c + b) := by ring
  rw [this]
  field_simp
  ring

lemma tendsto_div_sq_add_atTop (b : ℝ) (hb : 0 < b) :
    Tendsto (fun r : ℝ => r / (b ^ 2 + r ^ 2)) atTop (𝓝 0) ∧
    Tendsto (fun r : ℝ => r / (b ^ 2 + r ^ 2)) atBot (𝓝 0) := by
  have key : ∀ r : ℝ, r ≠ 0 → ‖r / (b ^ 2 + r ^ 2)‖ ≤ |r|⁻¹ := by
    intro r hr
    have hr' : 0 < |r| := abs_pos.mpr hr
    rw [Real.norm_eq_abs, abs_div, abs_of_pos (show 0 < b ^ 2 + r ^ 2 by positivity),
      div_le_iff₀ (by positivity), inv_mul_eq_div, le_div_iff₀ hr']
    nlinarith [sq_abs r, sq_nonneg b]
  constructor
  · refine squeeze_zero_norm' ?_ (tendsto_inv_atTop_zero.comp tendsto_abs_atTop_atTop)
    filter_upwards [eventually_ne_atTop 0] with r hr using key r hr
  · refine squeeze_zero_norm' ?_ (tendsto_inv_atTop_zero.comp tendsto_abs_atBot_atTop)
    filter_upwards [eventually_ne_atBot 0] with r hr using key r hr

lemma integral_inv_sq_sq {b : ℝ} (hb : 0 < b) :
    ∫ r : ℝ, 1 / ((b ^ 2 + r ^ 2) * (b ^ 2 + r ^ 2)) = π / (b * b * (b + b)) := by
  have hderiv : ∀ r : ℝ, HasDerivAt (fun r : ℝ => r / (b ^ 2 + r ^ 2))
      (2 * b ^ 2 * (1 / ((b ^ 2 + r ^ 2) * (b ^ 2 + r ^ 2))) - 1 / (b ^ 2 + r ^ 2)) r := by
    intro r
    have h1 : b ^ 2 + r ^ 2 ≠ 0 := by positivity
    have hd : HasDerivAt (fun r : ℝ => b ^ 2 + r ^ 2) (2 * r) r := by
      simpa using ((hasDerivAt_id r).pow 2).const_add (b ^ 2)
    have := (hasDerivAt_id' r).fun_div hd h1
    refine this.congr_deriv ?_
    field_simp
    ring
  have hint : Integrable (fun r : ℝ =>
      2 * b ^ 2 * (1 / ((b ^ 2 + r ^ 2) * (b ^ 2 + r ^ 2))) - 1 / (b ^ 2 + r ^ 2)) :=
    ((integrable_inv_mul_inv hb hb).const_mul _).sub (integrable_inv_sq_add hb)
  have h0 := integral_of_hasDerivAt_of_tendsto hderiv hint
    (tendsto_div_sq_add_atTop b hb).2 (tendsto_div_sq_add_atTop b hb).1
  rw [integral_sub ((integrable_inv_mul_inv hb hb).const_mul _) (integrable_inv_sq_add hb),
    integral_const_mul, integral_inv_sq_add hb] at h0
  have h2 : ∫ r : ℝ, 1 / ((b ^ 2 + r ^ 2) * (b ^ 2 + r ^ 2)) = (π / b) / (2 * b ^ 2) := by
    rw [eq_div_iff (by positivity)]
    linarith
  rw [h2]
  field_simp
  ring

lemma integral_inv_mul_inv {b c : ℝ} (hb : 0 < b) (hc : 0 < c) :
    ∫ r : ℝ, 1 / ((b ^ 2 + r ^ 2) * (c ^ 2 + r ^ 2)) = π / (b * c * (b + c)) := by
  by_cases hbc : b = c
  · subst hbc
    exact integral_inv_sq_sq hb
  · exact integral_inv_mul_inv_of_ne hb hc hbc

/-! ### The kernel `30/((9/4+r²)(25/4+r²))` -/

/-- The majorant's Fourier transform `30/((9/4+r²)(25/4+r²))`. -/
def majK (r : ℝ) : ℝ := 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))

lemma majK_eq (r : ℝ) :
    majK r = 15 / 2 * (1 / ((3 / 2) ^ 2 + r ^ 2) - 1 / ((5 / 2) ^ 2 + r ^ 2)) := by
  unfold majK
  have h1 : (9 / 4 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have h2 : (25 / 4 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have h3 : ((3 / 2) ^ 2 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have h4 : ((5 / 2) ^ 2 + r ^ 2 : ℝ) ≠ 0 := by positivity
  field_simp
  ring

lemma majK_nonneg (r : ℝ) : 0 ≤ majK r := by
  unfold majK
  positivity

lemma integrable_majK : Integrable majK := by
  have h := ((integrable_inv_sq_add (b := 3 / 2) (by norm_num)).sub
    (integrable_inv_sq_add (b := 5 / 2) (by norm_num))).const_mul (15 / 2)
  exact h.congr (Eventually.of_forall fun r => (majK_eq r).symm)

lemma integral_majK : ∫ r : ℝ, majK r = 2 * π := by
  simp_rw [majK_eq]
  rw [integral_const_mul, integral_sub (integrable_inv_sq_add (by norm_num))
    (integrable_inv_sq_add (by norm_num)), integral_inv_sq_add (by norm_num),
    integral_inv_sq_add (by norm_num)]
  ring

lemma majK_mul_eq {c : ℝ} (hc : 0 < c) (r : ℝ) :
    majK r * (2 * c / (c ^ 2 + r ^ 2))
      = 15 * c * (1 / (((3 / 2) ^ 2 + r ^ 2) * (c ^ 2 + r ^ 2))
        - 1 / (((5 / 2) ^ 2 + r ^ 2) * (c ^ 2 + r ^ 2))) := by
  rw [majK_eq]
  have h3 : ((3 / 2) ^ 2 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have h4 : ((5 / 2) ^ 2 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have h5 : (c ^ 2 + r ^ 2 : ℝ) ≠ 0 := by positivity
  field_simp

lemma integrable_majK_mul {c : ℝ} (hc : 0 < c) :
    Integrable (fun r : ℝ => majK r * (2 * c / (c ^ 2 + r ^ 2))) := by
  have h := ((integrable_inv_mul_inv (b := 3 / 2) (by norm_num) hc).sub
    (integrable_inv_mul_inv (b := 5 / 2) (by norm_num) hc)).const_mul (15 * c)
  exact h.congr (Eventually.of_forall fun r => (majK_mul_eq hc r).symm)

lemma integral_majK_mul {c : ℝ} (hc : 0 < c) :
    ∫ r : ℝ, majK r * (2 * c / (c ^ 2 + r ^ 2))
      = 2 * π * (5 / (3 / 2 + c) - 3 / (5 / 2 + c)) := by
  simp_rw [majK_mul_eq hc]
  rw [integral_const_mul, integral_sub (integrable_inv_mul_inv (by norm_num) hc)
    (integrable_inv_mul_inv (by norm_num) hc), integral_inv_mul_inv (by norm_num) hc,
    integral_inv_mul_inv (by norm_num) hc]
  field_simp
  ring

/-! ### The series terms -/

/-- The `n`-th term of the partial-fraction series of `Re ψ(3/4 + ir/2)`. -/
def psiTerm (n : ℕ) (r : ℝ) : ℝ :=
  1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 3 / 4) / (((n : ℝ) + 1 + 3 / 4) ^ 2 + (r / 2) ^ 2)

lemma psiTerm_eq (n : ℕ) (r : ℝ) :
    psiTerm n r = 1 / ((n : ℝ) + 1)
      - 2 * (2 * (n : ℝ) + 7 / 2) / ((2 * (n : ℝ) + 7 / 2) ^ 2 + r ^ 2) := by
  unfold psiTerm
  have h1 : (((n : ℝ) + 1 + 3 / 4) ^ 2 + (r / 2) ^ 2) ≠ 0 := by positivity
  have h2 : ((2 * (n : ℝ) + 7 / 2) ^ 2 + r ^ 2) ≠ 0 := by positivity
  congr 1
  field_simp
  ring

lemma psiTerm_nonneg (n : ℕ) (r : ℝ) : 0 ≤ psiTerm n r := by
  unfold psiTerm
  rw [sub_nonneg, div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [sq_nonneg (r / 2), Nat.cast_nonneg (α := ℝ) n]

lemma integrable_majK_mul_psiTerm (n : ℕ) :
    Integrable (fun r : ℝ => majK r * psiTerm n r) := by
  have hc : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by positivity
  have h := (integrable_majK.mul_const (1 / ((n : ℝ) + 1))).sub (integrable_majK_mul hc)
  refine h.congr (Eventually.of_forall fun r => ?_)
  simp only [Pi.sub_apply]
  rw [psiTerm_eq]
  ring

lemma integral_majK_mul_psiTerm (n : ℕ) :
    ∫ r : ℝ, majK r * psiTerm n r
      = 2 * π * (1 / ((n : ℝ) + 1) - 5 / (2 * (n : ℝ) + 5) + 3 / (2 * (n : ℝ) + 6)) := by
  have hc : (0 : ℝ) < 2 * (n : ℝ) + 7 / 2 := by positivity
  have h : (fun r : ℝ => majK r * psiTerm n r)
      = fun r : ℝ => majK r * (1 / ((n : ℝ) + 1))
        - majK r * (2 * (2 * (n : ℝ) + 7 / 2) / ((2 * (n : ℝ) + 7 / 2) ^ 2 + r ^ 2)) := by
    funext r
    rw [psiTerm_eq]
    ring
  rw [h, integral_sub (integrable_majK.mul_const _) (integrable_majK_mul hc),
    integral_mul_const, integral_majK, integral_majK_mul hc]
  rw [show (3 / 2 + (2 * (n : ℝ) + 7 / 2)) = 2 * (n : ℝ) + 5 by ring,
    show (5 / 2 + (2 * (n : ℝ) + 7 / 2)) = 2 * (n : ℝ) + 6 by ring]
  ring

/-! ### The numerical series -/

lemma three_halves_mem : (3 / 2 : ℂ) ∈ Complex.integerComplement := by
  rintro ⟨k, hk⟩
  have h2 : ((2 * k : ℤ) : ℂ) = ((3 : ℤ) : ℂ) := by
    push_cast
    rw [hk]
    norm_num
  have h3 : 2 * k = 3 := by exact_mod_cast h2
  omega

lemma digamma_three_halves :
    Complex.digamma (3 / 2) = ((2 - Real.eulerMascheroniConstant - 2 * Real.log 2 : ℝ) : ℂ) := by
  have hs : ∀ m : ℕ, (1 / 2 : ℂ) ≠ -m := by
    intro m h
    have h1 := congrArg Complex.re h
    simp at h1
    have : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
    linarith
  have h := Complex.digamma_apply_add_one (1 / 2) hs
  rw [show (1 / 2 : ℂ) + 1 = 3 / 2 by norm_num, Complex.digamma_one_half] at h
  rw [h]
  push_cast
  rw [Complex.ofReal_log (by norm_num)]
  push_cast
  ring

lemma hasSum_three_halves :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) - 1 / (3 / 2 + (n : ℝ) + 1))
      (8 / 3 - 2 * Real.log 2) := by
  have h := Zeta23.DigammaSeries.hasSum_digamma_series three_halves_mem
  rw [digamma_three_halves] at h
  rw [← Complex.hasSum_ofReal]
  convert h using 1
  · funext n
    push_cast
    ring
  · push_cast
    ring

lemma hasSum_telescope :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 3)) (3 / 2) := by
  have hpart : ∀ N : ℕ, ∑ n ∈ Finset.range N, (1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 3))
      = 3 / 2 - 1 / ((N : ℝ) + 1) - 1 / ((N : ℝ) + 2) := by
    intro N
    induction N with
    | zero => norm_num
    | succ N ih =>
      rw [Finset.sum_range_succ, ih]
      push_cast
      field_simp
      ring
  rw [hasSum_iff_tendsto_nat_of_nonneg]
  · simp_rw [hpart]
    have h1 : Tendsto (fun N : ℕ => 1 / ((N : ℝ) + 1)) atTop (𝓝 0) :=
      tendsto_one_div_add_atTop_nhds_zero_nat
    have h2 : Tendsto (fun N : ℕ => 1 / ((N : ℝ) + 2)) atTop (𝓝 0) := by
      have := h1.comp (tendsto_add_atTop_nat 1)
      refine this.congr fun N => ?_
      simp only [Function.comp_apply]
      push_cast
      ring_nf
    have := (tendsto_const_nhds (x := (3 / 2 : ℝ))).sub h1 |>.sub h2
    simpa using this
  · intro n
    rw [sub_nonneg]
    apply one_div_le_one_div_of_le (by positivity)
    linarith

lemma hasSum_arch_coeffs :
    HasSum (fun n : ℕ => 1 / ((n : ℝ) + 1) - 5 / (2 * (n : ℝ) + 5) + 3 / (2 * (n : ℝ) + 6))
      (53 / 12 - 5 * Real.log 2) := by
  have h := (hasSum_three_halves.mul_left (5 / 2)).sub (hasSum_telescope.mul_left (3 / 2))
  have hf : (fun n : ℕ => 1 / ((n : ℝ) + 1) - 5 / (2 * (n : ℝ) + 5) + 3 / (2 * (n : ℝ) + 6))
      = fun n : ℕ => 5 / 2 * (1 / ((n : ℝ) + 1) - 1 / (3 / 2 + (n : ℝ) + 1))
        - 3 / 2 * (1 / ((n : ℝ) + 1) - 1 / ((n : ℝ) + 3)) := by
    funext n
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (2 * (n : ℝ) + 5) ≠ 0 := by positivity
    have h3 : (2 * (n : ℝ) + 6) ≠ 0 := by positivity
    have h4 : (3 / 2 + (n : ℝ) + 1) ≠ 0 := by positivity
    have h5 : ((n : ℝ) + 3) ≠ 0 := by positivity
    field_simp
    ring
  have hv : (53 / 12 - 5 * Real.log 2 : ℝ)
      = 5 / 2 * (8 / 3 - 2 * Real.log 2) - 3 / 2 * (3 / 2) := by ring
  rw [hf, hv]
  exact h

/-! ### Termwise integration -/

lemma integrable_tsum_of_nonneg {F : ℕ → ℝ → ℝ} (hint : ∀ n, Integrable (F n))
    (hnn : ∀ n r, 0 ≤ F n r) (hpt : ∀ r, Summable (fun n => F n r))
    (hsum : Summable (fun n => ∫ r, F n r)) :
    Integrable (fun r => ∑' n, F n r) := by
  refine ⟨?_, ?_⟩
  · refine aestronglyMeasurable_of_tendsto_ae atTop
      (f := fun N r => ∑ n ∈ Finset.range N, F n r) (fun N => ?_)
      (Eventually.of_forall fun r => ?_)
    · exact (integrable_finsetSum _ fun n _ => hint n).aestronglyMeasurable
    · exact (hpt r).hasSum.tendsto_sum_nat
  · rw [hasFiniteIntegral_iff_ofReal (Eventually.of_forall fun r => tsum_nonneg fun n => hnn n r)]
    have h1 : ∀ r, ENNReal.ofReal (∑' n, F n r) = ∑' n, ENNReal.ofReal (F n r) := fun r =>
      ENNReal.ofReal_tsum_of_nonneg (fun n => hnn n r) (hpt r)
    simp_rw [h1]
    rw [lintegral_tsum (fun n => (hint n).aemeasurable.ennreal_ofReal)]
    have h2 : ∀ n, ∫⁻ r, ENNReal.ofReal (F n r) = ENNReal.ofReal (∫ r, F n r) := fun n =>
      (ofReal_integral_eq_lintegral_ofReal (hint n) (Eventually.of_forall (hnn n))).symm
    simp_rw [h2]
    rw [← ENNReal.ofReal_tsum_of_nonneg (fun n => integral_nonneg (hnn n)) hsum]
    exact ENNReal.ofReal_lt_top

lemma summable_psiTerm (r : ℝ) : Summable (fun n : ℕ => psiTerm n r) :=
  Zeta23.MuFields.summable_re_terms (a := 3 / 4) (by norm_num) (by norm_num) (r / 2)

lemma summable_integral_majK_mul_psiTerm :
    Summable (fun n : ℕ => ∫ r : ℝ, majK r * psiTerm n r) := by
  simp_rw [integral_majK_mul_psiTerm]
  exact (hasSum_arch_coeffs.mul_left (2 * π)).summable

lemma integrable_tsum_majK_mul_psiTerm :
    Integrable (fun r : ℝ => ∑' n : ℕ, majK r * psiTerm n r) :=
  integrable_tsum_of_nonneg integrable_majK_mul_psiTerm
    (fun n r => mul_nonneg (majK_nonneg r) (psiTerm_nonneg n r))
    (fun r => (summable_psiTerm r).mul_left _) summable_integral_majK_mul_psiTerm

lemma integral_tsum_majK_mul_psiTerm :
    ∫ r : ℝ, ∑' n : ℕ, majK r * psiTerm n r = 2 * π * (53 / 12 - 5 * Real.log 2) := by
  have hnorm : Summable (fun n : ℕ => ∫ r : ℝ, ‖majK r * psiTerm n r‖) := by
    have : ∀ n : ℕ, ∀ r : ℝ, ‖majK r * psiTerm n r‖ = majK r * psiTerm n r := fun n r =>
      by rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (majK_nonneg r) (psiTerm_nonneg n r))]
    simp_rw [this]
    exact summable_integral_majK_mul_psiTerm
  have h1 := hasSum_integral_of_summable_integral_norm integrable_majK_mul_psiTerm hnorm
  have h2 := hasSum_arch_coeffs.mul_left (2 * π)
  simp_rw [← integral_majK_mul_psiTerm] at h2
  exact h1.unique h2

/-! ### Assembly -/

lemma re_digamma_arch (r : ℝ) :
    (Complex.digamma (3 / 4 + I * r / 2)).re
      = -Real.eulerMascheroniConstant - 2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2)
        + ∑' n : ℕ, psiTerm n r := by
  have h := Zeta23.MuFields.re_digamma_vertical (a := 3 / 4) (by norm_num) (by norm_num) (r / 2)
  have hz : (3 / 4 + I * r / 2 : ℂ) = (((3 / 4 : ℝ)) : ℂ) + I * ((r / 2 : ℝ) : ℂ) := by
    push_cast
    ring
  rw [hz, h]
  have h1 : ((3 / 4 : ℝ) ^ 2 + (r / 2) ^ 2) ≠ 0 := by positivity
  have h2 : ((3 / 2 : ℝ) ^ 2 + r ^ 2) ≠ 0 := by positivity
  have h3 : (3 / 4 : ℝ) / ((3 / 4) ^ 2 + (r / 2) ^ 2) = 2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2) := by
    field_simp
    ring
  rw [h3]
  rfl

lemma archIntegrand_eq (r : ℝ) :
    30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))
        * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / π))
      = majK r * (-Real.eulerMascheroniConstant + Real.log (4 / π))
        - majK r * (2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2))
        + ∑' n : ℕ, majK r * psiTerm n r := by
  rw [re_digamma_arch, tsum_mul_left]
  unfold majK
  ring

end PoissonAux

/-- The archimedean integrand for the majorant is integrable. -/
theorem integrable_majFT_mul_archWeight :
    Integrable (fun r : ℝ => 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))
      * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / π))) := by
  have h := ((PoissonAux.integrable_majK.mul_const
      (-Real.eulerMascheroniConstant + Real.log (4 / π))).sub
    (PoissonAux.integrable_majK_mul (c := 3 / 2) (by norm_num))).add
    PoissonAux.integrable_tsum_majK_mul_psiTerm
  exact h.congr (Eventually.of_forall fun r => (PoissonAux.archIntegrand_eq r).symm)

/-- The archimedean term for the majorant, in closed form. -/
theorem integral_majFT_mul_archWeight :
    ∫ r : ℝ, 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))
        * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / π))
      = 2 * π * (7 / 2 - Real.eulerMascheroniConstant - 3 * Real.log 2 - Real.log π) := by
  simp_rw [PoissonAux.archIntegrand_eq]
  have h1 : Integrable (fun r : ℝ =>
      PoissonAux.majK r * (-Real.eulerMascheroniConstant + Real.log (4 / π))) :=
    PoissonAux.integrable_majK.mul_const _
  have h2 : Integrable (fun r : ℝ =>
      PoissonAux.majK r * (2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2))) :=
    PoissonAux.integrable_majK_mul (c := 3 / 2) (by norm_num)
  have h12 : Integrable (fun r : ℝ =>
      PoissonAux.majK r * (-Real.eulerMascheroniConstant + Real.log (4 / π))
        - PoissonAux.majK r * (2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2))) := h1.sub h2
  rw [integral_add h12 PoissonAux.integrable_tsum_majK_mul_psiTerm, integral_sub h1 h2,
    integral_mul_const, PoissonAux.integral_majK,
    PoissonAux.integral_majK_mul (c := 3 / 2) (by norm_num),
    PoissonAux.integral_tsum_majK_mul_psiTerm]
  have hlog : Real.log (4 / π) = 2 * Real.log 2 - Real.log π := by
    rw [Real.log_div (by norm_num) Real.pi_ne_zero, show (4 : ℝ) = 2 ^ 2 by norm_num,
      Real.log_pow]
    push_cast
    ring
  rw [hlog]
  ring

end ChebLib
