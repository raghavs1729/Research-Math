/-
ChebLib/Cheb/PsiAbel.lean — the Abel-summed von Mangoldt sum of `χ₋₄` under GRH.

With `h_x(u) = e^{u/2} exp(−e^u/x)` the test function `k = h_x − h_1` satisfies the hypotheses of
the extended explicit formula `EF_ext`; its transform is `(x^{1/2+ir} − 1) Γ(1/2+ir)` and its prime
side is `ψ(x) − ψ(1) + D(x)`, where `ψ(x) = Σ_n Λ(n) χ₋₄(n) e^{−n/x}` and `D(x)` is bounded.  Hence,
under GRH for `L(s,χ₋₄)`,

    |ψ(x)| ≤ (Σ_ρ m_ρ |Γ(1/2 + iγ_ρ)| + o(1)) √x      (x → ∞),

the archimedean term being `o(√x)` by the Riemann–Lebesgue lemma.
-/
import ChebLib.Cheb.EFExt
import ChebLib.Cheb.Kernel
import ChebLib.Cheb.AbelSum
import ChebLib.Race.S1

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace ChebLib

/-! ### Notation -/

/-- `χ₋₄(n)` as a real number. -/
def X4 (n : ℕ) : ℝ := ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)

/-- `Λ(n)`. -/
abbrev vM (n : ℕ) : ℝ := ArithmeticFunction.vonMangoldt n

theorem one_lt_four : 1 < 4 := by norm_num

/-- The zero configuration of `L(s,χ₋₄)`. -/
abbrev Z4 : ZeroConfig := LZeros (LSeam_of (χ := chi4) one_lt_four RaceLib.chi4_isPrimitive)

/-- `Σ_ρ m_ρ |Γ(1/2 + i Im ρ)|` over the nontrivial zeros of `L(s,χ₋₄)`. -/
def gammaSum : ℝ :=
  ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖

/-- `ψ(x) = Σ_n Λ(n) χ₋₄(n) e^{−n/x}`. -/
def psiAbel (x : ℝ) : ℝ := ∑' n : ℕ, vM n * X4 n * Real.exp (-(n / x))

theorem abs_X4_le_one (n : ℕ) : |X4 n| ≤ 1 := RaceLib.chi4_abs_le_one n

theorem vM_nonneg (n : ℕ) : 0 ≤ vM n := ArithmeticFunction.vonMangoldt_nonneg

theorem vM_le_self (n : ℕ) : vM n ≤ n :=
  ArithmeticFunction.vonMangoldt_le_log.trans (Real.log_le_self (Nat.cast_nonneg n))

theorem abs_vM_mul_X4_le (n : ℕ) : |vM n * X4 n| ≤ 1 * (n : ℝ) := by
  rw [abs_mul, abs_of_nonneg (vM_nonneg n), one_mul]
  calc vM n * |X4 n| ≤ vM n * 1 := mul_le_mul_of_nonneg_left (abs_X4_le_one n) (vM_nonneg n)
    _ ≤ n := by rw [mul_one]; exact vM_le_self n

theorem summable_psiAbel {x : ℝ} (hx : 0 < x) :
    Summable (fun n : ℕ => vM n * X4 n * Real.exp (-(n / x))) :=
  summable_mul_exp_neg_div _ 1 abs_vM_mul_X4_le hx

theorem summable_vM_div_sq : Summable (fun n : ℕ => vM n / (n : ℝ) ^ 2) := by
  have h := ArithmeticFunction.LSeriesSummable_vonMangoldt (s := 2) (by norm_num)
  refine Summable.of_norm_bounded h.norm ?_
  intro n
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [vM]
  · rw [LSeries.norm_term_eq]
    have hne : n ≠ 0 := Nat.pos_iff_ne_zero.mp hn
    simp only [hne, if_false]
    have e1 : ‖(ArithmeticFunction.vonMangoldt n : ℂ)‖ = ArithmeticFunction.vonMangoldt n := by
      rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (vM_nonneg n)]
    have e2 : ((n : ℝ)) ^ (2 : ℂ).re = (n : ℝ) ^ 2 := by
      have : (2 : ℂ).re = ((2 : ℕ) : ℝ) := by simp
      rw [this, Real.rpow_natCast]
    rw [e1, e2, Real.norm_eq_abs, abs_of_nonneg (div_nonneg (vM_nonneg n) (by positivity))]

/-! ### The archimedean weight of `χ₋₄` -/

theorem archWeight_eq_muq (r : ℝ) :
    archWeight chi4 r = 2 * Real.pi * muq (parity chi4) 4 r := by
  unfold archWeight muq
  rw [RaceLib.archParity_eq_parity]
  have hpi : (Real.pi : ℝ) ≠ 0 := Real.pi_ne_zero
  push_cast
  field_simp
  ring

theorem continuous_archWeight : Continuous (archWeight chi4) := by
  have h : archWeight chi4 = fun r => 2 * Real.pi * muq (parity chi4) 4 r :=
    funext archWeight_eq_muq
  rw [h]
  exact continuous_const.mul (GammaChi.muq_smooth (parity chi4) 4).continuous

theorem exists_archWeight_bound :
    ∃ K : ℝ, 0 ≤ K ∧ ∀ r : ℝ, |archWeight chi4 r| ≤ K * (1 + |r|) ^ (1 / 2 : ℝ) := by
  obtain ⟨K₂, hK₂, hmu⟩ := abs_muq_le_of_gammaFactsChi
    (GammaChi.gammaFactsChi (parity_le_one' chi4) (by norm_num : 1 ≤ 4)) (by norm_num : 1 ≤ 4)
  refine ⟨2 * Real.pi * K₂, by positivity, fun r => ?_⟩
  rw [archWeight_eq_muq, abs_mul, abs_of_pos (by positivity : (0:ℝ) < 2 * Real.pi)]
  calc 2 * Real.pi * |muq (parity chi4) 4 r|
      ≤ 2 * Real.pi * (K₂ * (1 + |r|) ^ (1 / 2 : ℝ)) :=
        mul_le_mul_of_nonneg_left (hmu r) (by positivity)
    _ = _ := by ring

/-! ### The Gamma factor on the critical line -/

theorem continuous_Gamma_half : Continuous (fun r : ℝ => Complex.Gamma (1 / 2 + I * (r : ℂ))) := by
  rw [continuous_iff_continuousAt]
  intro r
  have hs : ∀ m : ℕ, (1 / 2 + I * (r : ℂ)) ≠ -(m : ℂ) := by
    intro m h
    have := congrArg Complex.re h
    simp at this
    have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
    linarith
  have hc : ContinuousAt (fun r : ℝ => (1 / 2 : ℂ) + I * (r : ℂ)) r := by fun_prop
  exact ContinuousAt.comp (g := Complex.Gamma) (Complex.continuousAt_Gamma _ hs) hc

theorem norm_Gamma_half_le (r : ℝ) :
    ‖Complex.Gamma (1 / 2 + I * (r : ℂ))‖ ≤ 225 / 16 * Real.sqrt Real.pi / (1 + r ^ 2) := by
  have h := norm_Gamma_half_mul_le r
  have hG : 0 ≤ ‖Complex.Gamma (1 / 2 + I * (r : ℂ))‖ := norm_nonneg _
  rw [le_div_iff₀ (by positivity)]
  refine le_trans ?_ h
  rw [mul_comm]
  refine mul_le_mul_of_nonneg_right ?_ hG
  nlinarith [sq_nonneg r, sq_nonneg (r ^ 2)]

/-- The function `F(r) = Γ(1/2 + ir) W(r)`, `W` the archimedean weight of `χ₋₄`. -/
def archF (r : ℝ) : ℂ := Complex.Gamma (1 / 2 + I * (r : ℂ)) * ((archWeight chi4 r : ℝ) : ℂ)

theorem integrable_archF : Integrable archF := by
  obtain ⟨K, hK, hW⟩ := exists_archWeight_bound
  have hmeas : AEStronglyMeasurable archF :=
    (continuous_Gamma_half.mul (Complex.continuous_ofReal.comp continuous_archWeight)).aestronglyMeasurable
  have hdom : Integrable (fun r : ℝ =>
      2 * (225 / 16 * Real.sqrt Real.pi) * K * (1 + ‖r‖) ^ (-(3 / 2) : ℝ)) :=
    (integrable_one_add_norm (by rw [Module.finrank_self]; norm_num)).const_mul _
  refine hdom.mono' hmeas (Eventually.of_forall fun r => ?_)
  unfold archF
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs]
  have hb : 0 < 1 + |r| := by linarith [abs_nonneg r]
  have hsq : (1 + |r|) ^ (2 : ℝ) ≤ 2 * (1 + r ^ 2) := by
    rw [Real.rpow_two]
    nlinarith [sq_abs r, sq_nonneg (1 - |r|)]
  have key : (1 + |r|) ^ (1 / 2 : ℝ) ≤ 2 * (1 + r ^ 2) * (1 + |r|) ^ (-(3 / 2) : ℝ) := by
    rw [show (1 / 2 : ℝ) = 2 + (-(3 / 2)) by norm_num, Real.rpow_add hb]
    exact mul_le_mul_of_nonneg_right hsq (Real.rpow_nonneg hb.le _)
  have hrp : 0 ≤ (1 + |r|) ^ (-(3 / 2) : ℝ) := Real.rpow_nonneg hb.le _
  have hc : (0:ℝ) ≤ 225 / 16 * Real.sqrt Real.pi := by positivity
  have h1 : (0:ℝ) < 1 + r ^ 2 := by positivity
  calc ‖Complex.Gamma (1 / 2 + I * (r : ℂ))‖ * |archWeight chi4 r|
      ≤ (225 / 16 * Real.sqrt Real.pi / (1 + r ^ 2)) * (K * (1 + |r|) ^ (1 / 2 : ℝ)) :=
        mul_le_mul (norm_Gamma_half_le r) (hW r) (abs_nonneg _) (by positivity)
    _ ≤ (225 / 16 * Real.sqrt Real.pi / (1 + r ^ 2))
          * (K * (2 * (1 + r ^ 2) * (1 + |r|) ^ (-(3 / 2) : ℝ))) := by gcongr
    _ = 2 * (225 / 16 * Real.sqrt Real.pi) * K * (1 + |r|) ^ (-(3 / 2) : ℝ) := by
        field_simp

/-! ### The zero sum of the Gamma factor -/

theorem Z4_re (hGRH : GRHFor chi4) (ρ : Z4.carrier) : (ρ : ℂ).re = 1 / 2 :=
  hGRH ρ ρ.2.2.1 ρ.2.2.2 ρ.2.1

theorem Z4_gammaOf (hGRH : GRHFor chi4) (ρ : Z4.carrier) :
    gammaOf (ρ : ℂ) = (((ρ : ℂ).im : ℝ) : ℂ) :=
  RaceLib.gammaOf_eq_im (Z4_re hGRH ρ)

theorem summable_Z4_inv_sq (hGRH : GRHFor chi4) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) / (1 + (ρ : ℂ).im ^ 2)) := by
  obtain ⟨A₀, hA₀, hloc⟩ := localCountChi (χ := chi4) one_lt_four RaceLib.chi4_isPrimitive
  have hloc' : ∀ t : ℝ, (Z4.N t (t + 1) : ℝ) ≤ A₀ * Real.log (|t| + 3) := by
    intro t; rw [LZeros_N]; exact hloc t
  refine (Zeta23.WeilEF.zero_sum_inv_sq_gen Z4 hA₀ hloc').congr fun ρ => ?_
  rw [Z4_gammaOf hGRH ρ, Complex.normSq_ofReal, sq]

theorem summable_gammaSum (hGRH : GRHFor chi4) :
    Summable (fun ρ : Z4.carrier =>
      (Z4.mult ρ : ℝ) * ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖) := by
  refine Summable.of_nonneg_of_le (fun ρ => mul_nonneg (Nat.cast_nonneg _) (norm_nonneg _))
    (fun ρ => ?_) ((summable_Z4_inv_sq hGRH).mul_left (225 / 16 * Real.sqrt Real.pi))
  calc (Z4.mult ρ : ℝ) * ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖
      ≤ (Z4.mult ρ : ℝ) * (225 / 16 * Real.sqrt Real.pi / (1 + (ρ : ℂ).im ^ 2)) :=
        mul_le_mul_of_nonneg_left (norm_Gamma_half_le _) (Nat.cast_nonneg _)
    _ = 225 / 16 * Real.sqrt Real.pi * ((Z4.mult ρ : ℝ) / (1 + (ρ : ℂ).im ^ 2)) := by ring

/-! ### The test function `h_x − h_1` -/

/-- `k_x = h_x − h_1`, complex valued. -/
def kAbel (x : ℝ) (u : ℝ) : ℂ := ((abelKernel x u - abelKernel 1 u : ℝ) : ℂ)

def kAbel1 (x : ℝ) (u : ℝ) : ℂ := ((abelKernel1 x u - abelKernel1 1 u : ℝ) : ℂ)

def kAbel2 (x : ℝ) (u : ℝ) : ℂ := ((abelKernel2 x u - abelKernel2 1 u : ℝ) : ℂ)

theorem hasDerivAt_kAbel (x u : ℝ) : HasDerivAt (kAbel x) (kAbel1 x u) u :=
  ((hasDerivAt_abelKernel x u).sub (hasDerivAt_abelKernel 1 u)).ofReal_comp

theorem hasDerivAt_kAbel1 (x u : ℝ) : HasDerivAt (kAbel1 x) (kAbel2 x u) u :=
  ((hasDerivAt_abelKernel1 x u).sub (hasDerivAt_abelKernel1 1 u)).ofReal_comp

theorem continuous_kAbel2 (x : ℝ) : Continuous (kAbel2 x) :=
  Complex.continuous_ofReal.comp ((continuous_abelKernel2 x).sub (continuous_abelKernel2 1))

theorem integrable_kAbel {x : ℝ} (hx : 0 < x) : Integrable (kAbel x) :=
  ((integrable_abelKernel hx).sub (integrable_abelKernel one_pos)).ofReal

theorem integrable_kAbel1 {x : ℝ} (hx : 0 < x) : Integrable (kAbel1 x) :=
  ((integrable_abelKernel1 hx).sub (integrable_abelKernel1 one_pos)).ofReal

theorem integrable_kAbel2 {x : ℝ} (hx : 0 < x) : Integrable (kAbel2 x) :=
  ((integrable_abelKernel2 hx).sub (integrable_abelKernel2 one_pos)).ofReal

theorem integrable_abelKernel_mul_cexp {x : ℝ} (hx : 0 < x) (τ : ℝ) :
    Integrable (fun u : ℝ => (abelKernel x u : ℂ) * cexp (I * τ * u)) := by
  refine (integrable_abelKernel hx).ofReal.mul_bdd (c := 1) (by fun_prop) ?_
  filter_upwards with u
  rw [RaceLib.norm_cexp_I_mul_ofReal_mul_ofReal]

/-- The transform of `k_x`: `(x^{1/2+iτ} − 1) Γ(1/2+iτ)`. -/
theorem paperFT_kAbel {x : ℝ} (hx : 0 < x) (τ : ℝ) :
    paperFT (kAbel x) τ
      = ((x : ℂ) ^ ((1 / 2 : ℂ) + I * τ) - 1) * Complex.Gamma (1 / 2 + I * τ) := by
  unfold paperFT kAbel
  have e : ∀ u : ℝ, ((abelKernel x u - abelKernel 1 u : ℝ) : ℂ) * cexp (I * τ * u)
      = (abelKernel x u : ℂ) * cexp (I * τ * u) - (abelKernel 1 u : ℂ) * cexp (I * τ * u) := by
    intro u; push_cast; ring
  simp_rw [e]
  rw [integral_sub (integrable_abelKernel_mul_cexp hx τ) (integrable_abelKernel_mul_cexp one_pos τ),
    integral_abelKernel_mul_cexp hx, integral_abelKernel_mul_cexp one_pos]
  simp only [Complex.ofReal_one, Complex.one_cpow]
  ring

theorem cpow_half_add {x : ℝ} (hx : 0 < x) (r : ℝ) :
    (x : ℂ) ^ ((1 / 2 : ℂ) + I * r) = (Real.sqrt x : ℂ) * cexp (I * r * (Real.log x : ℂ)) := by
  have hx0 : (x : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hx.ne'
  rw [Complex.cpow_add _ _ hx0]
  congr 1
  · rw [Real.sqrt_eq_rpow, Complex.ofReal_cpow hx.le]
    norm_num
  · rw [Complex.cpow_def_of_ne_zero hx0, Complex.ofReal_log hx.le]
    congr 1
    ring

theorem norm_cpow_half_sub_one_le {x : ℝ} (hx : 0 < x) (r : ℝ) :
    ‖(x : ℂ) ^ ((1 / 2 : ℂ) + I * r) - 1‖ ≤ Real.sqrt x + 1 := by
  refine (norm_sub_le _ _).trans ?_
  rw [norm_one, cpow_half_add hx, norm_mul, RaceLib.norm_cexp_I_mul_ofReal_mul_ofReal,
    Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _), mul_one]

/-! ### The prime side of `k_x` -/

/-- The bounded dual term `D(x) = Σ_n Λ(n) χ₋₄(n) n^{-1} (e^{−1/(nx)} − e^{−1/n})`. -/
def dualTerm (x : ℝ) (n : ℕ) : ℝ :=
  vM n * X4 n / n * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1))))

theorem exp_diff_bounds {x : ℝ} (hx : 1 ≤ x) {n : ℕ} (hn : 0 < n) :
    0 ≤ Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1)))
      ∧ Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1))) ≤ 1 / n := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hx0 : 0 < x := by linarith
  have hab : 1 / ((n : ℝ) * x) ≤ 1 / ((n : ℝ) * 1) := by
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith
  constructor
  · have := Real.exp_le_exp.mpr (neg_le_neg hab)
    linarith
  · have h1 : Real.exp (-(1 / ((n : ℝ) * x))) ≤ 1 := by
      rw [Real.exp_le_one_iff]
      have : 0 < 1 / ((n : ℝ) * x) := by positivity
      linarith
    have h2 : 1 - 1 / ((n : ℝ) * 1) ≤ Real.exp (-(1 / ((n : ℝ) * 1))) := by
      have := Real.add_one_le_exp (-(1 / ((n : ℝ) * 1)))
      linarith
    rw [mul_one] at h2 ⊢
    linarith

theorem abs_dualTerm_le {x : ℝ} (hx : 1 ≤ x) (n : ℕ) : |dualTerm x n| ≤ vM n / (n : ℝ) ^ 2 := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [dualTerm, vM]
  · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    obtain ⟨h0, h1⟩ := exp_diff_bounds hx hn
    unfold dualTerm
    rw [abs_mul, abs_div, abs_mul, abs_of_nonneg (vM_nonneg n), abs_of_pos hn', abs_of_nonneg h0]
    calc vM n * |X4 n| / n * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1))))
        ≤ vM n * 1 / n * (1 / n) :=
          mul_le_mul (div_le_div_of_nonneg_right
            (mul_le_mul_of_nonneg_left (abs_X4_le_one n) (vM_nonneg n)) hn'.le) h1 h0
            (div_nonneg (mul_nonneg (vM_nonneg n) zero_le_one) hn'.le)
      _ = vM n / (n : ℝ) ^ 2 := by field_simp

theorem summable_dualTerm {x : ℝ} (hx : 1 ≤ x) : Summable (dualTerm x) :=
  Summable.of_norm_bounded summable_vM_div_sq fun n => by
    rw [Real.norm_eq_abs]; exact abs_dualTerm_le hx n

theorem abs_tsum_dualTerm_le {x : ℝ} (hx : 1 ≤ x) :
    |∑' n : ℕ, dualTerm x n| ≤ ∑' n : ℕ, vM n / (n : ℝ) ^ 2 := by
  rw [← Real.norm_eq_abs]
  refine (norm_tsum_le_tsum_norm (summable_dualTerm hx).norm).trans ?_
  exact (summable_dualTerm hx).norm.tsum_le_tsum (fun n => by
    rw [Real.norm_eq_abs]; exact abs_dualTerm_le hx n) summable_vM_div_sq

theorem kAbel_log (x : ℝ) {n : ℕ} (hn : 0 < n) :
    kAbel x (Real.log n)
      = ((Real.sqrt n * (Real.exp (-(n / x)) - Real.exp (-(n / 1))) : ℝ) : ℂ) := by
  unfold kAbel
  rw [abelKernel_log x hn, abelKernel_log 1 hn]
  congr 1
  ring

theorem kAbel_neg_log (x : ℝ) {n : ℕ} (hn : 0 < n) :
    kAbel x (-Real.log n)
      = (((Real.sqrt n)⁻¹ * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1)))) : ℝ) : ℂ) := by
  unfold kAbel
  rw [abelKernel_neg_log x hn, abelKernel_neg_log 1 hn]
  congr 1
  ring

theorem coeff_chi4 (n : ℕ) : coeff chi4 n = ((X4 n : ℝ) : ℂ) := RaceLib.chi4_natCast n

/-- The `n`-th term of the prime side of the explicit formula for `k_x`. -/
theorem primeTerm_kAbel (x : ℝ) (n : ℕ) :
    ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
        * (coeff chi4 n * kAbel x (Real.log n)
          + (starRingEnd ℂ) (coeff chi4 n) * kAbel x (-Real.log n))
      = ((vM n * X4 n * Real.exp (-(n / x)) - vM n * X4 n * Real.exp (-(n / 1))
          + dualTerm x n : ℝ) : ℂ) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [dualTerm, vM]
  · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
    have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
    rw [coeff_chi4, Complex.conj_ofReal, kAbel_log x hn, kAbel_neg_log x hn]
    unfold dualTerm vM
    generalize Real.exp (-(n / x)) = a
    generalize Real.exp (-(n / 1)) = b
    generalize Real.exp (-(1 / (n * x))) = c
    generalize Real.exp (-(1 / (n * 1))) = d
    norm_cast
    generalize Real.sqrt n = s at hs hss
    rw [hss]
    field_simp

theorem summable_primeSide_kAbel {x : ℝ} (hx : 1 ≤ x) :
    Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖kAbel x (Real.log n)‖ + ‖kAbel x (-Real.log n)‖)) := by
  have hx0 : 0 < x := by linarith
  have hb1 : Summable (fun n : ℕ => (n : ℝ) * Real.exp (-(n / x))) :=
    summable_mul_exp_neg_div (fun n => (n : ℝ)) 1 (fun n => by simp) hx0
  have hb2 : Summable (fun n : ℕ => (n : ℝ) * Real.exp (-(n / 1))) :=
    summable_mul_exp_neg_div (fun n => (n : ℝ)) 1 (fun n => by simp) one_pos
  refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_)
    ((hb1.add hb2).add summable_vM_div_sq)
  · exact mul_nonneg (div_nonneg (vM_nonneg n) (Real.sqrt_nonneg _)) (by positivity)
  · rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp
    · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
      have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
      obtain ⟨h0, h1⟩ := exp_diff_bounds hx hn
      rw [kAbel_log x hn, kAbel_neg_log x hn, Complex.norm_real, Complex.norm_real,
        Real.norm_eq_abs, Real.norm_eq_abs, abs_mul, abs_mul, abs_of_pos hs,
        abs_of_pos (inv_pos.mpr hs), abs_of_nonneg h0]
      have hE : |Real.exp (-(n / x)) - Real.exp (-(n / 1))|
          ≤ Real.exp (-(n / x)) + Real.exp (-(n / 1)) := by
        rw [abs_le]
        constructor <;> linarith [Real.exp_pos (-((n : ℝ) / x)), Real.exp_pos (-((n : ℝ) / 1))]
      have hΛ := vM_nonneg n
      have hΛn := vM_le_self n
      have e : ArithmeticFunction.vonMangoldt n / Real.sqrt n
          * (Real.sqrt n * |Real.exp (-(n / x)) - Real.exp (-(n / 1))|
            + (Real.sqrt n)⁻¹ * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1)))))
          = vM n * |Real.exp (-(n / x)) - Real.exp (-(n / 1))|
            + vM n / n * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1)))) := by
        have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
        unfold vM
        generalize |Real.exp (-(n / x)) - Real.exp (-(n / 1))| = a
        generalize Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1))) = c
        generalize Real.sqrt n = s at hs hss
        rw [hss]
        field_simp
      rw [e]
      have t1 : vM n * |Real.exp (-(n / x)) - Real.exp (-(n / 1))|
          ≤ (n : ℝ) * Real.exp (-(n / x)) + (n : ℝ) * Real.exp (-(n / 1)) := by
        calc vM n * |Real.exp (-(n / x)) - Real.exp (-(n / 1))|
            ≤ (n : ℝ) * (Real.exp (-(n / x)) + Real.exp (-(n / 1))) :=
              mul_le_mul hΛn hE (abs_nonneg _) hn'.le
          _ = _ := by ring
      have t2 : vM n / n * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1))))
          ≤ vM n / (n : ℝ) ^ 2 := by
        calc vM n / n * (Real.exp (-(1 / (n * x))) - Real.exp (-(1 / (n * 1))))
            ≤ vM n / n * (1 / n) := mul_le_mul_of_nonneg_left h1 (div_nonneg hΛ hn'.le)
          _ = _ := by field_simp
      linarith

/-! ### The explicit formula for `k_x` -/

/-- The explicit formula of `L(s,χ₋₄)` for the Abel kernel, under GRH: for `x ≥ 1`,
`Σ_ρ m_ρ (x^{1/2+iγ_ρ} − 1) Γ(1/2+iγ_ρ) = −(ψ(x) − ψ(1) + D(x))
  + (1/2π) ∫ (x^{1/2+ir} − 1) Γ(1/2+ir) W(r) dr`. -/
theorem EF_kAbel (hGRH : GRHFor chi4) {x : ℝ} (hx : 1 ≤ x) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℂ)
        * (((x : ℂ) ^ ((1 / 2 : ℂ) + I * (((ρ : ℂ).im : ℝ) : ℂ)) - 1)
          * Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ)))) ∧
      ∑' ρ : Z4.carrier, (Z4.mult ρ : ℂ)
        * (((x : ℂ) ^ ((1 / 2 : ℂ) + I * (((ρ : ℂ).im : ℝ) : ℂ)) - 1)
          * Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ)))
      = -((psiAbel x - psiAbel 1 + ∑' n : ℕ, dualTerm x n : ℝ) : ℂ)
        + (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, ((x : ℂ) ^ ((1 / 2 : ℂ) + I * (r : ℂ)) - 1) * archF r := by
  have hx0 : 0 < x := by linarith
  obtain ⟨hsum, hEF⟩ := EF_ext (χ := chi4) one_lt_four RaceLib.chi4_isPrimitive hGRH
    (hasDerivAt_kAbel x) (hasDerivAt_kAbel1 x) (continuous_kAbel2 x) (integrable_kAbel hx0)
    (integrable_kAbel1 hx0) (integrable_kAbel2 hx0) (summable_primeSide_kAbel hx)
  have hterm : ∀ ρ : Z4.carrier, (Z4.mult ρ : ℂ) * paperFT (kAbel x) (gammaOf ρ)
      = (Z4.mult ρ : ℂ) * (((x : ℂ) ^ ((1 / 2 : ℂ) + I * (((ρ : ℂ).im : ℝ) : ℂ)) - 1)
          * Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))) := by
    intro ρ
    rw [Z4_gammaOf hGRH ρ, paperFT_kAbel hx0]
  refine ⟨hsum.congr hterm, ?_⟩
  rw [← tsum_congr hterm, hEF]
  unfold literatureRHSchi
  congr 1
  · congr 1
    rw [tsum_congr (primeTerm_kAbel x), ← Complex.ofReal_tsum]
    congr 1
    have s1 := summable_psiAbel hx0
    have s2 := summable_psiAbel one_pos
    have s3 := summable_dualTerm hx
    rw [(s1.sub s2).tsum_add s3, s1.tsum_sub s2]
    rfl
  · congr 1
    refine integral_congr_ae (Eventually.of_forall fun r => ?_)
    simp only
    rw [paperFT_kAbel hx0, mul_assoc]
    rfl

/-! ### The bound for `ψ(x)` -/

theorem tendsto_archJ :
    Tendsto (fun x : ℝ => ∫ r : ℝ, cexp (I * r * (Real.log x : ℂ)) * archF r) atTop (𝓝 0) :=
  (RaceLib.tendsto_integral_cexp_mul archF).comp Real.tendsto_log_atTop

theorem integral_cpow_mul_archF {x : ℝ} (hx : 0 < x) :
    ∫ r : ℝ, ((x : ℂ) ^ ((1 / 2 : ℂ) + I * (r : ℂ)) - 1) * archF r
      = (Real.sqrt x : ℂ) * (∫ r : ℝ, cexp (I * r * (Real.log x : ℂ)) * archF r)
        - ∫ r : ℝ, archF r := by
  have hI : Integrable (fun r : ℝ => cexp (I * r * (Real.log x : ℂ)) * archF r) := by
    refine integrable_archF.bdd_mul (c := 1) (by fun_prop) ?_
    filter_upwards with r
    rw [RaceLib.norm_cexp_I_mul_ofReal_mul_ofReal]
  rw [← integral_const_mul, ← integral_sub (hI.const_mul _) integrable_archF]
  refine integral_congr_ae (Eventually.of_forall fun r => ?_)
  simp only
  rw [cpow_half_add hx]
  ring

/-- **The Abel-summed von Mangoldt sum of `χ₋₄` under GRH**: for every `ε > 0`,
`|Σ_n Λ(n) χ₋₄(n) e^{−n/x}| ≤ (Σ_ρ m_ρ |Γ(1/2+iγ_ρ)| + ε) √x` for all large `x`. -/
theorem abs_psiAbel_le (hGRH : GRHFor chi4) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ x : ℝ in atTop, |psiAbel x| ≤ (gammaSum + ε) * Real.sqrt x := by
  set D₀ : ℝ := ∑' n : ℕ, vM n / (n : ℝ) ^ 2 with hD₀
  set C₁ : ℝ := |psiAbel 1| + D₀ + gammaSum + ‖∫ r : ℝ, archF r‖ with hC₁
  have hgs0 : 0 ≤ gammaSum :=
    tsum_nonneg fun ρ => mul_nonneg (Nat.cast_nonneg _) (norm_nonneg _)
  have hJ := tendsto_archJ
  rw [Metric.tendsto_nhds] at hJ
  have hJε := hJ ε hε
  have hsqrt : ∀ᶠ x : ℝ in atTop, C₁ ≤ ε / 2 * Real.sqrt x := by
    have h := (Real.tendsto_sqrt_atTop.const_mul_atTop (by positivity : (0:ℝ) < ε / 2))
    exact h.eventually_ge_atTop C₁
  filter_upwards [hJε, hsqrt, eventually_ge_atTop (1 : ℝ)] with x hJx hCx hx
  have hx0 : 0 < x := by linarith
  obtain ⟨hsum, hEF⟩ := EF_kAbel hGRH hx
  rw [dist_zero_right] at hJx
  -- bound for the zero side
  have hzero : ‖∑' ρ : Z4.carrier, (Z4.mult ρ : ℂ)
        * (((x : ℂ) ^ ((1 / 2 : ℂ) + I * (((ρ : ℂ).im : ℝ) : ℂ)) - 1)
          * Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ)))‖
      ≤ (Real.sqrt x + 1) * gammaSum := by
    refine (norm_tsum_le_tsum_norm hsum.norm).trans ?_
    unfold gammaSum
    rw [← tsum_mul_left]
    refine hsum.norm.tsum_le_tsum (fun ρ => ?_) ((summable_gammaSum hGRH).mul_left _)
    rw [norm_mul, norm_mul, Complex.norm_natCast]
    calc (Z4.mult ρ : ℝ) * (‖(x : ℂ) ^ ((1 / 2 : ℂ) + I * (((ρ : ℂ).im : ℝ) : ℂ)) - 1‖
          * ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖)
        ≤ (Z4.mult ρ : ℝ) * ((Real.sqrt x + 1)
          * ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖) := by
          gcongr
          exact norm_cpow_half_sub_one_le hx0 _
      _ = _ := by ring
  -- bound for the archimedean side
  have harch : ‖(1 / (2 * Real.pi) : ℂ)
        * ∫ r : ℝ, ((x : ℂ) ^ ((1 / 2 : ℂ) + I * (r : ℂ)) - 1) * archF r‖
      ≤ Real.sqrt x * ε / 2 + ‖∫ r : ℝ, archF r‖ := by
    rw [integral_cpow_mul_archF hx0, norm_mul]
    have hpi : ‖(1 / (2 * Real.pi) : ℂ)‖ ≤ 1 / 2 := by
      have : (1 / (2 * Real.pi) : ℂ) = ((1 / (2 * Real.pi) : ℝ) : ℂ) := by push_cast; ring
      rw [this, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
      rw [div_le_div_iff₀ (by positivity) (by norm_num)]
      linarith [Real.pi_gt_three]
    have h2 : ‖(Real.sqrt x : ℂ) * (∫ r : ℝ, cexp (I * r * (Real.log x : ℂ)) * archF r)
        - ∫ r : ℝ, archF r‖ ≤ Real.sqrt x * ε + ‖∫ r : ℝ, archF r‖ := by
      refine (norm_sub_le _ _).trans ?_
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Real.sqrt_nonneg _)]
      gcongr
    have h3 : (0:ℝ) ≤ Real.sqrt x * ε + ‖∫ r : ℝ, archF r‖ := by positivity
    calc ‖(1 / (2 * Real.pi) : ℂ)‖ * _ ≤ 1 / 2 * (Real.sqrt x * ε + ‖∫ r : ℝ, archF r‖) :=
          mul_le_mul hpi h2 (norm_nonneg _) (by norm_num)
      _ ≤ Real.sqrt x * ε / 2 + ‖∫ r : ℝ, archF r‖ := by
          linarith [norm_nonneg (∫ r : ℝ, archF r)]
  -- assemble
  have hD : |∑' n : ℕ, dualTerm x n| ≤ D₀ := abs_tsum_dualTerm_le hx
  have hkey : ((psiAbel x - psiAbel 1 + ∑' n : ℕ, dualTerm x n : ℝ) : ℂ)
      = -(∑' ρ : Z4.carrier, (Z4.mult ρ : ℂ)
        * (((x : ℂ) ^ ((1 / 2 : ℂ) + I * (((ρ : ℂ).im : ℝ) : ℂ)) - 1)
          * Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))))
        + (1 / (2 * Real.pi) : ℂ)
          * ∫ r : ℝ, ((x : ℂ) ^ ((1 / 2 : ℂ) + I * (r : ℂ)) - 1) * archF r := by
    rw [hEF]; ring
  have hnorm : |psiAbel x - psiAbel 1 + ∑' n : ℕ, dualTerm x n|
      ≤ (Real.sqrt x + 1) * gammaSum + (Real.sqrt x * ε / 2 + ‖∫ r : ℝ, archF r‖) := by
    rw [← Real.norm_eq_abs, ← Complex.norm_real, hkey]
    refine (norm_add_le _ _).trans (add_le_add ?_ harch)
    rw [norm_neg]
    exact hzero
  have h4 : |psiAbel x| ≤ |psiAbel x - psiAbel 1 + ∑' n : ℕ, dualTerm x n| + |psiAbel 1|
      + |∑' n : ℕ, dualTerm x n| := by
    have e : psiAbel x = (psiAbel x - psiAbel 1 + ∑' n : ℕ, dualTerm x n) + psiAbel 1
        - ∑' n : ℕ, dualTerm x n := by ring
    calc |psiAbel x| = |(psiAbel x - psiAbel 1 + ∑' n : ℕ, dualTerm x n) + psiAbel 1
          - ∑' n : ℕ, dualTerm x n| := by rw [← e]
      _ ≤ _ := (abs_sub _ _).trans (add_le_add (abs_add_le _ _) le_rfl)
  have hsx : 0 ≤ Real.sqrt x := Real.sqrt_nonneg x
  calc |psiAbel x| ≤ (Real.sqrt x + 1) * gammaSum + (Real.sqrt x * ε / 2 + ‖∫ r : ℝ, archF r‖)
        + |psiAbel 1| + D₀ := by linarith
    _ = gammaSum * Real.sqrt x + ε / 2 * Real.sqrt x + C₁ := by rw [hC₁]; ring
    _ ≤ gammaSum * Real.sqrt x + ε / 2 * Real.sqrt x + ε / 2 * Real.sqrt x := by linarith
    _ = (gammaSum + ε) * Real.sqrt x := by ring

end ChebLib
