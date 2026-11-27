/-
RaceLib/Abel.lean — a smoothed Abel-summation lemma (Mathlib only).

If the partial sums `C(x) = Σ_{n ≤ x} c(n)` of a complex sequence satisfy `C(x)/√x → 1`, then for
`g ∈ C¹_c(ℝ)`:   `Σ_n c(n) n^{-1/2} g(log n − L) → (1/2) ∫ g`  as `L → +∞`.
-/
import Mathlib

open Complex MeasureTheory Set Filter Topology
open scoped BigOperators

noncomputable section

namespace RaceLib

/-- The weight `t ↦ g(log t − L)/√t`. -/
def wfun (g : ℝ → ℂ) (L : ℝ) (t : ℝ) : ℂ := (Real.sqrt t)⁻¹ • g (Real.log t - L)

/-- The derivative of `wfun g L` at `t > 0`. -/
def wderiv (g : ℝ → ℂ) (L : ℝ) (t : ℝ) : ℂ :=
  (Real.sqrt t)⁻¹ • (t⁻¹ • deriv g (Real.log t - L))
    + (-(1 / (2 * Real.sqrt t)) / (Real.sqrt t) ^ 2) • g (Real.log t - L)

theorem hasDerivAt_comp_log_sub {g : ℝ → ℂ} (hg : Differentiable ℝ g) (L : ℝ) {t : ℝ}
    (ht : 0 < t) :
    HasDerivAt (fun t => g (Real.log t - L)) (t⁻¹ • deriv g (Real.log t - L)) t := by
  have h1 : HasDerivAt (fun t => Real.log t - L) t⁻¹ t :=
    (Real.hasDerivAt_log ht.ne').sub_const L
  exact (hg (Real.log t - L)).hasDerivAt.scomp t h1

theorem hasDerivAt_wfun {g : ℝ → ℂ} (hg : Differentiable ℝ g) (L : ℝ) {t : ℝ} (ht : 0 < t) :
    HasDerivAt (wfun g L) (wderiv g L t) t := by
  have hs : Real.sqrt t ≠ 0 := (Real.sqrt_pos.mpr ht).ne'
  have h3 : HasDerivAt (fun t => (Real.sqrt t)⁻¹)
      (-(1 / (2 * Real.sqrt t)) / (Real.sqrt t) ^ 2) t :=
    (Real.hasDerivAt_sqrt ht.ne').inv hs
  exact h3.smul (hasDerivAt_comp_log_sub hg L ht)

theorem wderiv_mul_sqrt (g : ℝ → ℂ) (L : ℝ) {t : ℝ} (ht : 0 < t) :
    wderiv g L t * (Real.sqrt t : ℂ)
      = t⁻¹ • deriv g (Real.log t - L) - (1 / 2 : ℂ) * (t⁻¹ • g (Real.log t - L)) := by
  have hs : Real.sqrt t ≠ 0 := (Real.sqrt_pos.mpr ht).ne'
  have hsC : ((Real.sqrt t : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hs
  have hsqC : (t : ℂ) = ((Real.sqrt t : ℝ) : ℂ) ^ 2 := by exact_mod_cast (Real.sq_sqrt ht.le).symm
  unfold wderiv
  simp only [Complex.real_smul]
  push_cast
  rw [hsqC]
  field_simp
  ring

theorem continuousOn_comp_log_sub {h : ℝ → ℂ} (hh : Continuous h) (L : ℝ) {a b : ℝ} (ha : 0 < a) :
    ContinuousOn (fun t => h (Real.log t - L)) (Icc a b) :=
  hh.comp_continuousOn ((Real.continuousOn_log.mono
    (fun _ ht => (lt_of_lt_of_le ha ht.1).ne')).sub continuousOn_const)

theorem continuousOn_inv_Icc {a b : ℝ} (ha : 0 < a) : ContinuousOn (fun t : ℝ => t⁻¹) (Icc a b) :=
  continuousOn_inv₀.mono (fun _ ht => (lt_of_lt_of_le ha ht.1).ne')

theorem continuousOn_wderiv {g : ℝ → ℂ} (hg : ContDiff ℝ 1 g) (L : ℝ) {a b : ℝ} (ha : 0 < a) :
    ContinuousOn (wderiv g L) (Icc a b) := by
  have hsq : ContinuousOn (fun t : ℝ => Real.sqrt t) (Icc a b) := Real.continuous_sqrt.continuousOn
  have hne : ∀ t ∈ Icc a b, Real.sqrt t ≠ 0 := fun t ht =>
    (Real.sqrt_pos.mpr (lt_of_lt_of_le ha ht.1)).ne'
  have hg' : Continuous (deriv g) := hg.continuous_deriv le_rfl
  unfold wderiv
  refine ContinuousOn.add ?_ ?_
  · have h1 : ContinuousOn (fun t : ℝ => (Real.sqrt t)⁻¹) (Icc a b) := hsq.inv₀ hne
    have h2 : ContinuousOn (fun t : ℝ => t⁻¹ • deriv g (Real.log t - L)) (Icc a b) :=
      (continuousOn_inv_Icc ha).smul (continuousOn_comp_log_sub hg' L ha)
    exact h1.smul h2
  · have h1 : ContinuousOn (fun t : ℝ => -(1 / (2 * Real.sqrt t)) / (Real.sqrt t) ^ 2) (Icc a b) := by
      refine ContinuousOn.div ?_ (hsq.pow 2) (fun t ht => pow_ne_zero 2 (hne t ht))
      refine ContinuousOn.neg (ContinuousOn.div continuousOn_const (continuousOn_const.mul hsq) ?_)
      exact fun t ht => mul_ne_zero two_ne_zero (hne t ht)
    exact h1.smul (continuousOn_comp_log_sub hg.continuous L ha)

section
variable {g : ℝ → ℂ} {A : ℝ}

theorem wfun_eq_zero_of_le (hgs : tsupport g ⊆ Icc (-A) A) {L t : ℝ} (ht0 : 0 ≤ t)
    (ht : t ≤ Real.exp (L - A - 1)) : wfun g L t = 0 := by
  unfold wfun
  rcases ht0.eq_or_lt with h0 | hpos
  · rw [← h0]; simp
  · have hlog : Real.log t ≤ L - A - 1 := by
      rw [← Real.log_exp (L - A - 1)]
      exact Real.log_le_log hpos ht
    have : g (Real.log t - L) = 0 := by
      apply image_eq_zero_of_notMem_tsupport
      intro hmem
      have := (hgs hmem).1
      linarith
    rw [this, smul_zero]

theorem wfun_eq_zero_of_ge (hgs : tsupport g ⊆ Icc (-A) A) {L t : ℝ}
    (ht : Real.exp (L + A + 1) ≤ t) : wfun g L t = 0 := by
  unfold wfun
  have hlog : L + A + 1 ≤ Real.log t := by
    rw [← Real.log_exp (L + A + 1)]
    exact Real.log_le_log (Real.exp_pos _) ht
  have : g (Real.log t - L) = 0 := by
    apply image_eq_zero_of_notMem_tsupport
    intro hmem
    have := (hgs hmem).2
    linarith
  rw [this, smul_zero]

/-- The main term: `∫_a^b w'(t) √t dt = −(1/2) ∫ g`, `a = e^{L−A−1}`, `b = e^{L+A+1}`. -/
theorem integral_wderiv_mul_sqrt (hg : ContDiff ℝ 1 g) (hgs : tsupport g ⊆ Icc (-A) A)
    (hA : 0 ≤ A) (L : ℝ) :
    ∫ t in Ioc (Real.exp (L - A - 1)) (Real.exp (L + A + 1)), wderiv g L t * (Real.sqrt t : ℂ)
      = -((1 / 2 : ℂ) * ∫ u, g u) := by
  set a := Real.exp (L - A - 1) with ha_def
  set b := Real.exp (L + A + 1) with hb_def
  have ha : 0 < a := Real.exp_pos _
  have hab : a ≤ b := Real.exp_le_exp.mpr (by linarith)
  have hpos : ∀ t ∈ uIcc a b, 0 < t := fun t ht => by
    rw [uIcc_of_le hab] at ht; exact lt_of_lt_of_le ha ht.1
  have hg' : Continuous (deriv g) := hg.continuous_deriv le_rfl
  have hgd : Differentiable ℝ g := hg.differentiable one_ne_zero
  have hla : Real.log a - L = -A - 1 := by rw [ha_def, Real.log_exp]; ring
  have hlb : Real.log b - L = A + 1 := by rw [hb_def, Real.log_exp]; ring
  have hgz : ∀ x, x ∉ Icc (-A) A → g x = 0 := fun x hx =>
    image_eq_zero_of_notMem_tsupport (fun h => hx (hgs h))
  -- integrability
  have hc1 : ContinuousOn (fun t : ℝ => t⁻¹ • deriv g (Real.log t - L)) (uIcc a b) := by
    rw [uIcc_of_le hab]
    exact (continuousOn_inv_Icc ha).smul (continuousOn_comp_log_sub hg' L ha)
  have hc2 : ContinuousOn (fun t : ℝ => t⁻¹ • g (Real.log t - L)) (uIcc a b) := by
    rw [uIcc_of_le hab]
    exact (continuousOn_inv_Icc ha).smul (continuousOn_comp_log_sub hg.continuous L ha)
  rw [← intervalIntegral.integral_of_le hab]
  have hcongr : ∫ t in a..b, wderiv g L t * (Real.sqrt t : ℂ)
      = ∫ t in a..b, (t⁻¹ • deriv g (Real.log t - L) - (1 / 2 : ℂ) * (t⁻¹ • g (Real.log t - L))) :=
    intervalIntegral.integral_congr fun t ht => wderiv_mul_sqrt g L (hpos t ht)
  rw [hcongr, intervalIntegral.integral_sub hc1.intervalIntegrable
    ((hc2.intervalIntegrable).const_mul _), intervalIntegral.integral_const_mul]
  -- first piece: an exact derivative
  have h1 : ∫ t in a..b, t⁻¹ • deriv g (Real.log t - L) = 0 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt
      (fun t ht => hasDerivAt_comp_log_sub hgd L (hpos t ht)) hc1.intervalIntegrable]
    rw [hla, hlb, hgz (A + 1) (fun h => by linarith [h.2]),
      hgz (-A - 1) (fun h => by linarith [h.1]), sub_zero]
  -- second piece: substitution u = log t − L
  have h2 : ∫ t in a..b, t⁻¹ • g (Real.log t - L) = ∫ u, g u := by
    have hsub := intervalIntegral.integral_deriv_smul_comp (a := a) (b := b)
      (f := fun t => Real.log t - L) (f' := fun t => t⁻¹) (g := g)
      (fun t ht => (Real.hasDerivAt_log (hpos t ht).ne').sub_const L)
      (by rw [uIcc_of_le hab]; exact continuousOn_inv_Icc ha) hg.continuous
    simp only [Function.comp_apply] at hsub
    rw [hsub, hla, hlb]
    apply intervalIntegral.integral_eq_integral_of_support_subset
    intro x hx
    by_contra hmem
    exact hx (hgz x (fun h => hmem ⟨by linarith [h.1], by linarith [h.2]⟩))
  rw [h1, h2, zero_sub]

/-- Abel summation: the smoothed sum equals `(1/2)∫g` minus an error integral against
`C(t) − √t`. -/
theorem tsum_wfun_mul_eq (c : ℕ → ℂ) (hg : ContDiff ℝ 1 g) (hgs : tsupport g ⊆ Icc (-A) A)
    (hA : 0 ≤ A) (L : ℝ) :
    ∑' n : ℕ, wfun g L n * c n
      = (1 / 2 : ℂ) * (∫ u, g u)
        - ∫ t in Ioc (Real.exp (L - A - 1)) (Real.exp (L + A + 1)),
            wderiv g L t * ((∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ)) := by
  set a := Real.exp (L - A - 1) with ha_def
  set b := Real.exp (L + A + 1) with hb_def
  have ha : 0 < a := Real.exp_pos _
  have hab : a ≤ b := Real.exp_le_exp.mpr (by linarith)
  have hgd : Differentiable ℝ g := hg.differentiable one_ne_zero
  have hDc : ContinuousOn (wderiv g L) (Icc a b) := continuousOn_wderiv hg L ha
  have hderiv : ∀ t ∈ Icc a b, deriv (wfun g L) t = wderiv g L t := fun t ht =>
    (hasDerivAt_wfun hgd L (lt_of_lt_of_le ha ht.1)).deriv
  -- step 1: the series is a finite sum
  have hfin : ∑' n : ℕ, wfun g L n * c n = ∑ k ∈ Finset.Ioc ⌊a⌋₊ ⌊b⌋₊, wfun g L k * c k := by
    apply tsum_eq_sum
    intro n hn
    rw [Finset.mem_Ioc, not_and_or, not_lt, not_le] at hn
    rcases hn with hn | hn
    · rw [wfun_eq_zero_of_le hgs (Nat.cast_nonneg n)
        ((Nat.le_floor_iff ha.le).mp hn), zero_mul]
    · rw [wfun_eq_zero_of_ge hgs (Nat.lt_of_floor_lt hn).le, zero_mul]
  -- step 2: Abel summation
  have hint : IntegrableOn (deriv (wfun g L)) (Icc a b) :=
    (hDc.integrableOn_Icc).congr_fun (fun t ht => (hderiv t ht).symm) measurableSet_Icc
  have habel := sum_mul_eq_sub_sub_integral_mul c ha.le hab
    (fun t ht => (hasDerivAt_wfun hgd L (lt_of_lt_of_le ha ht.1)).differentiableAt) hint
  rw [wfun_eq_zero_of_le hgs ha.le le_rfl, wfun_eq_zero_of_ge hgs le_rfl, zero_mul, zero_mul,
    sub_zero, zero_sub] at habel
  -- step 3: replace deriv by wderiv and split C = √ + (C − √)
  have hcongr : ∫ t in Ioc a b, deriv (wfun g L) t * ∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k
      = ∫ t in Ioc a b, wderiv g L t * ∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k :=
    setIntegral_congr_fun measurableSet_Ioc (fun t ht => by
      simp only [hderiv t (Ioc_subset_Icc_self ht)])
  have hI1 : IntegrableOn (fun t => wderiv g L t * ∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) (Ioc a b) :=
    (integrableOn_mul_sum_Icc c ha.le hDc.integrableOn_Icc).mono_set Ioc_subset_Icc_self
  have hI2 : IntegrableOn (fun t => wderiv g L t * (Real.sqrt t : ℂ)) (Ioc a b) := by
    refine (ContinuousOn.integrableOn_Icc ?_).mono_set Ioc_subset_Icc_self
    exact hDc.mul (Complex.continuous_ofReal.comp Real.continuous_sqrt).continuousOn
  have hsplit : ∫ t in Ioc a b, wderiv g L t * ((∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ))
      = (∫ t in Ioc a b, wderiv g L t * ∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k)
        - ∫ t in Ioc a b, wderiv g L t * (Real.sqrt t : ℂ) := by
    rw [← integral_sub hI1 hI2]
    congr 1; ext t; ring
  rw [hfin, habel, hcongr, hsplit, integral_wderiv_mul_sqrt hg hgs hA L]
  ring

/-- **Smoothed Abel summation.**  If `C(x)/√x → 1` then
`Σ_n c(n) n^{-1/2} g(log n − L) → (1/2)∫g` as `L → +∞`. -/
theorem tendsto_tsum_wfun_mul {c : ℕ → ℂ}
    (hc : Tendsto (fun x : ℝ => (∑ k ∈ Finset.Icc 0 ⌊x⌋₊, c k) / (Real.sqrt x : ℂ)) atTop (𝓝 1))
    (hg : ContDiff ℝ 1 g) (hgs : tsupport g ⊆ Icc (-A) A) (hA : 0 ≤ A) :
    Tendsto (fun L : ℝ => ∑' n : ℕ, wfun g L n * c n) atTop (𝓝 ((1 / 2 : ℂ) * ∫ u, g u)) := by
  have hgc : HasCompactSupport g := IsCompact.of_isClosed_subset isCompact_Icc
    (isClosed_tsupport g) hgs
  obtain ⟨Mg, hMg⟩ := hgc.exists_bound_of_continuous hg.continuous
  obtain ⟨Mg', hMg'⟩ := hgc.deriv.exists_bound_of_continuous (hg.continuous_deriv le_rfl)
  have hMg0 : 0 ≤ Mg := (norm_nonneg _).trans (hMg 0)
  have hMg'0 : 0 ≤ Mg' := (norm_nonneg _).trans (hMg' 0)
  rw [Metric.tendsto_atTop]
  intro ε hε
  set K := (Mg' + Mg) * (Real.exp (2 * A + 2) - 1) with hK
  have hexp : 0 ≤ Real.exp (2 * A + 2) - 1 := by
    have := Real.add_one_le_exp (2 * A + 2); linarith
  have hK0 : 0 ≤ K := mul_nonneg (add_nonneg hMg'0 hMg0) hexp
  set δ := ε / (K + 1) with hδ
  have hδ0 : 0 < δ := div_pos hε (by linarith)
  have hδK : δ * K < ε := by
    rw [hδ, div_mul_eq_mul_div, div_lt_iff₀ (by linarith)]
    nlinarith
  obtain ⟨x₀, hx₀⟩ := (Metric.tendsto_atTop.mp hc) δ hδ0
  refine ⟨Real.log (max x₀ 1) + A + 1, fun L hL => ?_⟩
  set a := Real.exp (L - A - 1) with ha_def
  set b := Real.exp (L + A + 1) with hb_def
  have ha : 0 < a := Real.exp_pos _
  have hab : a ≤ b := Real.exp_le_exp.mpr (by linarith)
  have hax : max x₀ 1 ≤ a := by
    rw [ha_def, ← Real.exp_log (lt_of_lt_of_le one_pos (le_max_right x₀ 1))]
    exact Real.exp_le_exp.mpr (by linarith)
  rw [tsum_wfun_mul_eq c hg hgs hA L, dist_eq_norm, sub_sub_cancel_left, norm_neg]
  have hbound : ∀ t ∈ Ioc a b,
      ‖wderiv g L t * ((∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ))‖
        ≤ δ * (Mg' + Mg) / a := by
    intro t ht
    have ht0 : 0 < t := lt_trans ha ht.1
    have hs : 0 < Real.sqrt t := Real.sqrt_pos.mpr ht0
    have hsC : ((Real.sqrt t : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hs.ne'
    have hE : ‖(∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ)‖ ≤ δ * Real.sqrt t := by
      have h := hx₀ t ((le_max_left x₀ 1).trans (hax.trans ht.1.le))
      rw [dist_eq_norm] at h
      have heq : (∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ)
          = ((∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) / (Real.sqrt t : ℂ) - 1) * (Real.sqrt t : ℂ) := by
        field_simp
      rw [heq, norm_mul, Complex.norm_real, Real.norm_of_nonneg hs.le]
      exact mul_le_mul_of_nonneg_right h.le hs.le
    have hD : ‖wderiv g L t‖ * Real.sqrt t ≤ (Mg' + Mg) / t := by
      have h1 : ‖wderiv g L t‖ * Real.sqrt t = ‖wderiv g L t * (Real.sqrt t : ℂ)‖ := by
        rw [norm_mul, Complex.norm_real, Real.norm_of_nonneg hs.le]
      rw [h1, wderiv_mul_sqrt g L ht0]
      calc ‖t⁻¹ • deriv g (Real.log t - L) - (1 / 2 : ℂ) * (t⁻¹ • g (Real.log t - L))‖
          ≤ ‖t⁻¹ • deriv g (Real.log t - L)‖ + ‖(1 / 2 : ℂ) * (t⁻¹ • g (Real.log t - L))‖ :=
            norm_sub_le _ _
        _ ≤ t⁻¹ * Mg' + t⁻¹ * Mg := by
            rw [norm_smul, norm_mul, norm_smul, Real.norm_of_nonneg (inv_pos.mpr ht0).le]
            have hhalf : ‖(1 / 2 : ℂ)‖ ≤ 1 := by norm_num
            have hi : 0 ≤ t⁻¹ := (inv_pos.mpr ht0).le
            have e1 := mul_le_mul_of_nonneg_left (hMg' (Real.log t - L)) hi
            have e2 := mul_le_mul_of_nonneg_left (hMg (Real.log t - L)) hi
            have e3 : ‖(1 / 2 : ℂ)‖ * (t⁻¹ * ‖g (Real.log t - L)‖)
                ≤ 1 * (t⁻¹ * ‖g (Real.log t - L)‖) :=
              mul_le_mul_of_nonneg_right hhalf (mul_nonneg hi (norm_nonneg _))
            linarith
        _ = (Mg' + Mg) / t := by field_simp
    calc ‖wderiv g L t * ((∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ))‖
        = ‖wderiv g L t‖ * ‖(∑ k ∈ Finset.Icc 0 ⌊t⌋₊, c k) - (Real.sqrt t : ℂ)‖ := norm_mul _ _
      _ ≤ ‖wderiv g L t‖ * (δ * Real.sqrt t) := mul_le_mul_of_nonneg_left hE (norm_nonneg _)
      _ = δ * (‖wderiv g L t‖ * Real.sqrt t) := by ring
      _ ≤ δ * ((Mg' + Mg) / t) := mul_le_mul_of_nonneg_left hD hδ0.le
      _ ≤ δ * ((Mg' + Mg) / a) := by
          apply mul_le_mul_of_nonneg_left _ hδ0.le
          exact div_le_div_of_nonneg_left (add_nonneg hMg'0 hMg0) ha ht.1.le
      _ = δ * (Mg' + Mg) / a := by ring
  have hnorm := norm_setIntegral_le_of_norm_le_const (μ := volume) (s := Ioc a b)
    (by simp) hbound
  rw [Real.volume_real_Ioc_of_le hab] at hnorm
  refine lt_of_le_of_lt hnorm ?_
  have hba : (b - a) / a = Real.exp (2 * A + 2) - 1 := by
    rw [sub_div, div_self ha.ne', hb_def, ha_def, ← Real.exp_sub]
    congr 2; ring
  calc δ * (Mg' + Mg) / a * (b - a) = δ * ((Mg' + Mg) * ((b - a) / a)) := by ring
    _ = δ * K := by rw [hba]
    _ < ε := hδK

end

end RaceLib
