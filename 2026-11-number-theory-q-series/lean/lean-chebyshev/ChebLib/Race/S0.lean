/-
RaceLib/S0.lean — the GRH bound for the smoothed twisted prime-power sum (stage S0).

For a primitive Dirichlet character χ mod q > 1, a test function g = f ⋆ f̃ (autocorrelation of
f ∈ C²_c) and a shift L beyond the support of g, Zeta23's literature-form explicit formula
`Zeta23.ThmE.EF_lit_chi_L`, applied to u ↦ g(u − L) and to g, gives under GRH for L(s,χ)

    ‖S(L) − Γ_L‖ ≤ Σ_ρ m_ρ ĝ(γ_ρ) = Γ_0 − Σ_n Λ(n) n^{-1/2} (χ(n) g(log n) + conj χ(n) g(−log n)).
-/
import ChebLib.Race.Defs
import Zeta23.ThmE.MainChi
import Zeta23.ThmE.GammaFactsChiProof

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace RaceLib

/-! ### Bridges between the challenge definitions and Zeta23 -/

theorem raceFT_eq_paperFT : raceFT = paperFT := rfl

theorem autocorr_eq_weilTest (f : ℝ → ℂ) : autocorr f = weilTest f f := by
  funext x
  simp only [autocorr, weilTest, convolution_def, ContinuousLinearMap.mul_apply', tilde, neg_sub]

theorem archParity_eq_parity {q : ℕ} (χ : DirichletCharacter ℂ q) :
    archParity χ = parity χ := rfl

/-! ### Shifted test functions -/

/-- The shifted test function `u ↦ g(u − L)`. -/
def shift (g : ℝ → ℂ) (L : ℝ) : ℝ → ℂ := fun u => g (u - L)

theorem shift_contDiff {g : ℝ → ℂ} (hg : ContDiff ℝ 2 g) (L : ℝ) : ContDiff ℝ 2 (shift g L) :=
  hg.comp (contDiff_id.sub contDiff_const)

theorem shift_hasCompactSupport {g : ℝ → ℂ} (hg : HasCompactSupport g) (L : ℝ) :
    HasCompactSupport (shift g L) :=
  hg.comp_homeomorph (Homeomorph.subRight L)

/-- `(g(· − L))^(z) = e^{izL} ĝ(z)`. -/
theorem paperFT_shift (g : ℝ → ℂ) (L : ℝ) (z : ℂ) :
    paperFT (shift g L) z = cexp (I * z * L) * paperFT g z := by
  unfold paperFT shift
  rw [← cintegral_const_mul,
    ← integral_sub_right_eq_self (fun u : ℝ => cexp (I * z * L) * (g u * cexp (I * z * u))) L]
  congr 1; ext u
  rw [mul_left_comm, ← Complex.exp_add]
  congr 2
  push_cast
  ring

/-! ### GRH makes the ordinates real, and autocorrelations have nonnegative transform -/

theorem gammaOf_eq_im {ρ : ℂ} (h : ρ.re = 1 / 2) : gammaOf ρ = (ρ.im : ℂ) := by
  unfold gammaOf
  rw [div_eq_iff Complex.I_ne_zero]
  apply Complex.ext <;> simp [h]

theorem paperFT_weilTest_self_ofReal {f : ℝ → ℂ} (hf : Continuous f) (hfs : HasCompactSupport f)
    (t : ℝ) : paperFT (weilTest f f) t = (Complex.normSq (paperFT f t) : ℂ) := by
  rw [paperFT_weilTest hf hf hfs hfs, Complex.conj_ofReal, Complex.mul_conj]

theorem norm_cexp_I_mul_ofReal_mul_ofReal (t L : ℝ) : ‖cexp (I * (t : ℂ) * (L : ℂ))‖ = 1 := by
  rw [show I * (t : ℂ) * (L : ℂ) = ((t * L : ℝ) : ℂ) * I by push_cast; ring]
  exact Complex.norm_exp_ofReal_mul_I _

/-! ### The core inequality on the explicit-formula right-hand sides -/

variable {q : ℕ} [NeZero q] {χ : DirichletCharacter ℂ q}

/-- Under GRH for `L(s,χ)`: the zero side at shift `L` is dominated by the zero side at shift `0`
(which is a sum of nonnegative terms).  Stated on the right-hand sides of the explicit formula. -/
theorem norm_literatureRHSchi_shift_le (hq : 1 < q) (hprim : χ.IsPrimitive) (hGRH : GRHFor χ)
    {f : ℝ → ℂ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (L : ℝ) :
    ‖literatureRHSchi (parity χ) q (coeff χ) (shift (weilTest f f) L)‖
      ≤ (literatureRHSchi (parity χ) q (coeff χ) (weilTest f f)).re := by
  have hg2 : ContDiff ℝ 2 (weilTest f f) := weilTest_contDiff hf hf.continuous hfc
  have hgc : HasCompactSupport (weilTest f f) := weilTest_hasCompactSupport hfc hfc
  obtain ⟨hsumL, hEFL⟩ := EF_lit_chi_L hq hprim (shift (weilTest f f) L)
    (shift_contDiff hg2 L) (shift_hasCompactSupport hgc L)
  obtain ⟨-, hEF0⟩ := EF_lit_chi_L hq hprim (weilTest f f) hg2 hgc
  set Z := LZeros (LSeam_of hq hprim) with hZ
  have hre : ∀ ρ : Z.carrier, (ρ : ℂ).re = 1 / 2 := fun ρ =>
    hGRH ρ ρ.2.2.1 ρ.2.2.2 ρ.2.1
  set b : Z.carrier → ℝ := fun ρ => (Z.mult ρ : ℝ) * Complex.normSq (paperFT f ((ρ : ℂ).im : ℝ))
    with hb
  have h0 : ∀ ρ : Z.carrier,
      (Z.mult ρ : ℂ) * paperFT (weilTest f f) (gammaOf ρ) = ((b ρ : ℝ) : ℂ) := by
    intro ρ
    rw [gammaOf_eq_im (hre ρ), paperFT_weilTest_self_ofReal hf.continuous hfc, hb]
    push_cast
    rfl
  have hL : ∀ ρ : Z.carrier,
      ‖(Z.mult ρ : ℂ) * paperFT (shift (weilTest f f) L) (gammaOf ρ)‖ = b ρ := by
    intro ρ
    rw [gammaOf_eq_im (hre ρ), paperFT_shift, paperFT_weilTest_self_ofReal hf.continuous hfc,
      norm_mul, norm_mul, norm_cexp_I_mul_ofReal_mul_ofReal, Complex.norm_natCast,
      Complex.norm_real, Real.norm_of_nonneg (Complex.normSq_nonneg _), one_mul]
  rw [← hEFL, ← hEF0]
  calc ‖∑' ρ : Z.carrier, (Z.mult ρ : ℂ) * paperFT (shift (weilTest f f) L) (gammaOf ρ)‖
      ≤ ∑' ρ : Z.carrier, ‖(Z.mult ρ : ℂ) * paperFT (shift (weilTest f f) L) (gammaOf ρ)‖ :=
        norm_tsum_le_tsum_norm hsumL.norm
    _ = ∑' ρ : Z.carrier, b ρ := tsum_congr hL
    _ = (∑' ρ : Z.carrier, ((b ρ : ℝ) : ℂ)).re := by
        rw [← Complex.ofReal_tsum, Complex.ofReal_re]
    _ = (∑' ρ : Z.carrier, (Z.mult ρ : ℂ) * paperFT (weilTest f f) (gammaOf ρ)).re := by
        congr 1
        exact tsum_congr fun ρ => (h0 ρ).symm

/-! ### The right-hand sides in the challenge's notation -/

omit [NeZero q] in
/-- For a shift beyond the support, the `conj χ` leg of the prime side vanishes and the explicit
formula's right-hand side is `−S(L) + Γ_L`. -/
theorem literatureRHSchi_shift_eq (χ : DirichletCharacter ℂ q) {g : ℝ → ℂ} {A L : ℝ}
    (hgs : tsupport g ⊆ Icc (-A) A) (hL : A < L) :
    literatureRHSchi (parity χ) q (coeff χ) (shift g L)
      = -primePowerSum χ g L + archTerm χ g L := by
  unfold literatureRHSchi primePowerSum archTerm
  congr 1
  · congr 1
    refine tsum_congr fun n => ?_
    have hz : g (-Real.log n - L) = 0 := by
      apply image_eq_zero_of_notMem_tsupport
      intro hmem
      have h1 := (hgs hmem).1
      have h2 : 0 ≤ Real.log (n : ℝ) := Real.log_natCast_nonneg n
      linarith
    simp only [shift, coeff]
    rw [hz, mul_zero, add_zero, mul_assoc]
  · congr 1
    refine integral_congr_ae (Filter.Eventually.of_forall fun r => ?_)
    simp only [paperFT_shift, archWeight, archParity_eq_parity, raceFT_eq_paperFT]

omit [NeZero q] in
/-- At shift `0` the explicit formula's right-hand side is `Γ_0` minus the two-sided prime sum. -/
theorem literatureRHSchi_eq (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) :
    literatureRHSchi (parity χ) q (coeff χ) g = archTerm χ g 0 - primePowerSumSym χ g := by
  unfold literatureRHSchi primePowerSumSym archTerm
  rw [neg_add_eq_sub]
  congr 2
  refine integral_congr_ae (Filter.Eventually.of_forall fun r => ?_)
  simp only [archWeight, archParity_eq_parity, raceFT_eq_paperFT, Complex.ofReal_zero,
    mul_zero, Complex.exp_zero, one_mul]

/-! ### Stage S0 -/

/-- **S0, general primitive χ.**  Under GRH for `L(s,χ)`, for `g = f ⋆ f̃` with `f ∈ C²_c` supported
in `[−A/2, A/2]` and every shift `L > A`:
`‖S(L) − Γ_L‖ ≤ Re (Γ_0 − Σ_n Λ(n) n^{-1/2} (χ(n) g(log n) + conj χ(n) g(−log n)))`. -/
theorem norm_primePowerSum_sub_archTerm_le (hq : 1 < q) (hprim : χ.IsPrimitive) (hGRH : GRHFor χ)
    {f : ℝ → ℂ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {A : ℝ}
    (hA : tsupport f ⊆ Icc (-(A / 2)) (A / 2)) {L : ℝ} (hL : A < L) :
    ‖primePowerSum χ (autocorr f) L - archTerm χ (autocorr f) L‖
      ≤ (archTerm χ (autocorr f) 0 - primePowerSumSym χ (autocorr f)).re := by
  have h := norm_literatureRHSchi_shift_le hq hprim hGRH hf hfc L
  rw [literatureRHSchi_shift_eq χ (tsupport_weilTest_subset hA hA) hL, literatureRHSchi_eq] at h
  rw [autocorr_eq_weilTest, ← norm_neg, neg_sub, sub_eq_neg_add]
  exact h

/-- **The explicit formula at shift 0, in the challenge's notation** (unconditional): for
`g ∈ C²_c`, the zero sum converges absolutely and `Σ_ρ m_ρ ĝ(γ_ρ) = Γ_0 − (two-sided prime sum)`. -/
theorem zeroSum_eq (hq : 1 < q) (hprim : χ.IsPrimitive) {g : ℝ → ℂ} (hg : ContDiff ℝ 2 g)
    (hgc : HasCompactSupport g) :
    zeroSum χ g = archTerm χ g 0 - primePowerSumSym χ g := by
  rw [← literatureRHSchi_eq, ← (EF_lit_chi_L hq hprim g hg hgc).2]
  rfl

/-- The zero sum of the explicit formula converges absolutely (part of `EF_lit_chi_L`). -/
theorem summable_zeroSum (hq : 1 < q) (hprim : χ.IsPrimitive) {g : ℝ → ℂ} (hg : ContDiff ℝ 2 g)
    (hgc : HasCompactSupport g) :
    Summable (fun ρ : {ρ : ℂ | χ.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1} =>
      ((analyticOrderAt χ.LFunction (ρ : ℂ)).toNat : ℂ) * raceFT g (((ρ : ℂ) - 1 / 2) / I)) :=
  (EF_lit_chi_L hq hprim g hg hgc).1

/-- The archimedean integrand `ĝ(r) [Re ψ(1/4 + κ/2 + ir/2) + log(q/π)]` is integrable for
`g ∈ C²_c` (so the integrals `Γ_L` are genuine Bochner integrals). -/
theorem integrable_raceFT_mul_archWeight (χ : DirichletCharacter ℂ q) {g : ℝ → ℂ}
    (hg : ContDiff ℝ 2 g) (hgc : HasCompactSupport g) :
    Integrable (fun r : ℝ => raceFT g r * (archWeight χ r : ℂ)) := by
  have hq : 1 ≤ q := Nat.one_le_iff_ne_zero.mpr (NeZero.ne q)
  have h := integrable_paperFT_mul_muq hg hgc
    (GammaChi.gammaFactsChi (parity_le_one' χ) hq) hq
  have h2 := h.const_mul ((2 * Real.pi : ℝ) : ℂ)
  refine h2.congr (Filter.Eventually.of_forall fun r => ?_)
  have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  simp only [muq, archWeight, raceFT_eq_paperFT, archParity_eq_parity]
  push_cast
  field_simp
  ring

/-- **S0, zero-sum form.**  `‖S(L) − Γ_L‖ ≤ Re Σ_ρ m_ρ ĝ(γ_ρ)`. -/
theorem norm_primePowerSum_sub_archTerm_le_zeroSum (hq : 1 < q) (hprim : χ.IsPrimitive)
    (hGRH : GRHFor χ) {f : ℝ → ℂ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {A : ℝ}
    (hA : tsupport f ⊆ Icc (-(A / 2)) (A / 2)) {L : ℝ} (hL : A < L) :
    ‖primePowerSum χ (autocorr f) L - archTerm χ (autocorr f) L‖
      ≤ (zeroSum χ (autocorr f)).re := by
  have hg2 : ContDiff ℝ 2 (autocorr f) := by
    rw [autocorr_eq_weilTest]; exact weilTest_contDiff hf hf.continuous hfc
  have hgc : HasCompactSupport (autocorr f) := by
    rw [autocorr_eq_weilTest]; exact weilTest_hasCompactSupport hfc hfc
  rw [zeroSum_eq hq hprim hg2 hgc]
  exact norm_primePowerSum_sub_archTerm_le hq hprim hGRH hf hfc hA hL

end RaceLib
