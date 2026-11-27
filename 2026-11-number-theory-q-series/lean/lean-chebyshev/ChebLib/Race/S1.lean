/-
RaceLib/S1.lean — the smoothed race for χ₋₄ (stage S1).

Under GRH for `L(s,χ₋₄)` and the prime number theorem `θ(x) ~ x` (an explicit hypothesis), for
`g = f ⋆ f̃`, `f ∈ C²_c` supported in `[−a/2, a/2]` with `a < log 3`:

    limsup_{L→∞} Re Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L)  ≤  Re Γ_0 − (1/2) Re ĝ(0).

Ingredients: S0 (`‖S(L) − Γ_L‖ ≤ Re Γ_0`), the Riemann–Lebesgue lemma (`Γ_L → 0`), and smoothed
Abel summation against `Σ_{n ≤ x, n not prime} Λ(n) χ₋₄(n) = √x + o(√x)` (prime squares).
-/
import ChebLib.Race.Chi4
import ChebLib.Race.Abel
import ChebLib.Race.NonPrime
import ChebLib.Race.RiemannLebesgue

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

namespace RaceLib

/-- `Γ_L → 0` as `L → +∞` (Riemann–Lebesgue). -/
theorem tendsto_archTerm {q : ℕ} (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) :
    Tendsto (fun L : ℝ => archTerm χ g L) atTop (𝓝 0) := by
  have h := (tendsto_integral_cexp_mul (fun r => raceFT g r * (archWeight χ r : ℂ))).const_mul
    ((1 / (2 * Real.pi) : ℂ))
  rw [mul_zero] at h
  refine h.congr fun L => ?_
  unfold archTerm
  congr 1
  refine integral_congr_ae (Eventually.of_forall fun r => ?_)
  simp only [mul_assoc]

theorem chi4_natCast (n : ℕ) :
    chi4 (n : ZMod 4) = (((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) : ℂ) := by
  rw [Complex.ofReal_intCast]; rfl

/-- The prime powers that are not primes: `S(L) − P(L) = Σ_n c(n) n^{-1/2} g(log n − L)` with
`c = npCoeff`. -/
theorem primePowerSum_sub_primeSum_chi4 {g : ℝ → ℂ} {A : ℝ} (hgs : tsupport g ⊆ Icc (-A) A)
    (L : ℝ) :
    primePowerSum chi4 g L - primeSum chi4 g L = ∑' n : ℕ, wfun g L n * (npCoeff n : ℂ) := by
  set F : ℕ → ℂ := fun n => ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
    * chi4 (n : ZMod 4) * g (Real.log n - L) with hF
  set P : ℕ → ℂ := fun n => chi4 (n : ZMod 4) * ((Real.log n / Real.sqrt n : ℝ) : ℂ)
    * g (Real.log n - L) with hP
  have hzero : ∀ n ∉ Finset.range (⌊Real.exp (L + A)⌋₊ + 1), g (Real.log n - L) = 0 := by
    intro n hn
    rw [Finset.mem_range, not_lt] at hn
    have h1 : Real.exp (L + A) < n := Nat.lt_of_floor_lt (Nat.lt_of_succ_le hn)
    have h2 : L + A < Real.log n := by
      rw [← Real.log_exp (L + A)]
      exact Real.log_lt_log (Real.exp_pos _) h1
    apply image_eq_zero_of_notMem_tsupport
    intro hmem
    have := (hgs hmem).2
    linarith
  have hFs : Summable F := summable_of_ne_finset_zero (s := Finset.range (⌊Real.exp (L + A)⌋₊ + 1))
    (fun n hn => by simp only [hF, hzero n hn, mul_zero])
  have hPs : Summable ({n : ℕ | n.Prime}.indicator P) :=
    summable_of_ne_finset_zero (s := Finset.range (⌊Real.exp (L + A)⌋₊ + 1)) (fun n hn => by
      by_cases hp : n ∈ {n : ℕ | n.Prime}
      · rw [Set.indicator_of_mem hp]; simp only [hP, hzero n hn, mul_zero]
      · rw [Set.indicator_of_notMem hp])
  have hprime : primeSum chi4 g L = ∑' n : ℕ, {n : ℕ | n.Prime}.indicator P n :=
    tsum_subtype {n : ℕ | n.Prime} P
  have hpp : primePowerSum chi4 g L = ∑' n : ℕ, F n := rfl
  rw [hpp, hprime, ← hFs.tsum_sub hPs]
  refine tsum_congr fun n => ?_
  by_cases hp : n.Prime
  · rw [Set.indicator_of_mem (show n ∈ {n : ℕ | n.Prime} from hp)]
    have hc : npCoeff n = 0 := by unfold npCoeff; rw [if_pos hp]
    simp only [hF, hP, hc, ArithmeticFunction.vonMangoldt_apply_prime hp]
    push_cast
    ring
  · rw [Set.indicator_of_notMem (show n ∉ {n : ℕ | n.Prime} from hp), sub_zero]
    have hc : npCoeff n
        = ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) := by
      unfold npCoeff; rw [if_neg hp]
    simp only [hF, hc, wfun, chi4_natCast, Complex.real_smul]
    push_cast
    ring

/-- Prime squares: under PNT, `S(L) − P(L) → (1/2) ĝ(0)` as `L → +∞`. -/
theorem tendsto_primePowerSum_sub_primeSum_chi4 (hPNT : PNTTheta) {g : ℝ → ℂ}
    (hg : ContDiff ℝ 1 g) {A : ℝ} (hA : 0 ≤ A) (hgs : tsupport g ⊆ Icc (-A) A) :
    Tendsto (fun L : ℝ => primePowerSum chi4 g L - primeSum chi4 g L) atTop
      (𝓝 ((1 / 2 : ℂ) * raceFT g 0)) := by
  have hc : Tendsto (fun x : ℝ => (∑ k ∈ Finset.Icc 0 ⌊x⌋₊, (npCoeff k : ℂ)) / (Real.sqrt x : ℂ))
      atTop (𝓝 1) := by
    have h := (Complex.continuous_ofReal.tendsto 1).comp (tendsto_sum_npCoeff_div_sqrt hPNT)
    rw [Complex.ofReal_one] at h
    refine h.congr fun x => ?_
    simp only [Function.comp_apply]
    push_cast
    rfl
  have h := tendsto_tsum_wfun_mul hc hg hgs hA
  have hft : raceFT g 0 = ∫ u, g u := by simp [raceFT]
  rw [hft]
  exact h.congr fun L => (primePowerSum_sub_primeSum_chi4 hgs L).symm

/-- **S1, two-sided form.**  Under GRH for `L(s,χ₋₄)` and PNT, for every `ε > 0` and all large `L`:
`‖P(L) + (1/2) ĝ(0)‖ ≤ Re Γ_0 + ε`. -/
theorem chi4_norm_primeSum_add_le (hGRH : GRHFor chi4) (hPNT : PNTTheta)
    {f : ℝ → ℂ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {a : ℝ}
    (ha : a < Real.log 3) (hsupp : tsupport f ⊆ Icc (-(a / 2)) (a / 2)) {ε : ℝ} (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀,
      ‖primeSum chi4 (autocorr f) L + (1 / 2 : ℂ) * raceFT (autocorr f) 0‖
        ≤ (archTerm chi4 (autocorr f) 0).re + ε := by
  have hg2 : ContDiff ℝ 2 (autocorr f) := by
    rw [autocorr_eq_weilTest]; exact weilTest_contDiff hf hf.continuous hfc
  have hgs : tsupport (autocorr f) ⊆ Icc (-a) a := by
    rw [autocorr_eq_weilTest]; exact tsupport_weilTest_subset hsupp hsupp
  have hgs' : tsupport (autocorr f) ⊆ Icc (-(max a 0)) (max a 0) :=
    hgs.trans (Icc_subset_Icc (neg_le_neg (le_max_left _ _)) (le_max_left _ _))
  have hΓ := tendsto_archTerm chi4 (autocorr f)
  have hN := tendsto_primePowerSum_sub_primeSum_chi4 hPNT (hg2.of_le (by norm_num))
    (le_max_right a 0) hgs'
  obtain ⟨L₁, h1⟩ := (Metric.tendsto_atTop.mp hΓ) (ε / 2) (half_pos hε)
  obtain ⟨L₂, h2⟩ := (Metric.tendsto_atTop.mp hN) (ε / 2) (half_pos hε)
  refine ⟨max (max L₁ L₂) (a + 1), fun L hL => ?_⟩
  have hL1 : L ≥ L₁ := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hL
  have hL2 : L ≥ L₂ := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hL
  have hLa : a < L := lt_of_lt_of_le (by linarith) (le_trans (le_max_right _ _) hL)
  have b0 := chi4_norm_primePowerSum_sub_archTerm_le hGRH hf hfc ha hsupp hLa
  have b1 := h1 L hL1
  rw [dist_zero_right] at b1
  have b2 := h2 L hL2
  rw [dist_eq_norm] at b2
  set S := primePowerSum chi4 (autocorr f) L
  set P := primeSum chi4 (autocorr f) L
  set Γ := archTerm chi4 (autocorr f) L
  set h0 := raceFT (autocorr f) 0
  have hid : P + (1 / 2 : ℂ) * h0 = (S - Γ) + Γ - ((S - P) - (1 / 2 : ℂ) * h0) := by ring
  rw [hid]
  have t1 : ‖(S - Γ) + Γ - ((S - P) - (1 / 2 : ℂ) * h0)‖
      ≤ ‖(S - Γ) + Γ‖ + ‖(S - P) - (1 / 2 : ℂ) * h0‖ := norm_sub_le _ _
  have t2 : ‖(S - Γ) + Γ‖ ≤ ‖S - Γ‖ + ‖Γ‖ := norm_add_le _ _
  linarith

/-- **S1, race form.**  If moreover `Re Γ_0 < (1/2) Re ĝ(0) − margin`, then for all large `L`:
`Re P(L) < −margin`. -/
theorem chi4_primeSum_re_lt (hGRH : GRHFor chi4) (hPNT : PNTTheta)
    {f : ℝ → ℂ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {a : ℝ}
    (ha : a < Real.log 3) (hsupp : tsupport f ⊆ Icc (-(a / 2)) (a / 2)) {margin : ℝ}
    (hcert : (archTerm chi4 (autocorr f) 0).re
      < (1 / 2) * (raceFT (autocorr f) 0).re - margin) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, (primeSum chi4 (autocorr f) L).re < -margin := by
  obtain ⟨L₀, h⟩ := chi4_norm_primeSum_add_le hGRH hPNT hf hfc ha hsupp
    (ε := ((1 / 2) * (raceFT (autocorr f) 0).re - margin - (archTerm chi4 (autocorr f) 0).re) / 2)
    (by linarith)
  refine ⟨L₀, fun L hL => ?_⟩
  have h1 := h L hL
  have h2 : (primeSum chi4 (autocorr f) L + (1 / 2 : ℂ) * raceFT (autocorr f) 0).re
      ≤ ‖primeSum chi4 (autocorr f) L + (1 / 2 : ℂ) * raceFT (autocorr f) 0‖ :=
    Complex.re_le_norm _
  have h3 : ((1 / 2 : ℂ) * raceFT (autocorr f) 0).re = (1 / 2) * (raceFT (autocorr f) 0).re := by
    simp
  rw [Complex.add_re, h3] at h2
  linarith

/-! ### Real-valued windows -/

theorem autocorr_ofReal (f : ℝ → ℝ) (x : ℝ) :
    autocorr (fun x => (f x : ℂ)) x = ((∫ t, f t * f (t - x) : ℝ) : ℂ) := by
  unfold autocorr
  rw [← integral_complex_ofReal]
  congr 1; ext t
  simp

/-- For a real-valued window the smoothed prime sum for `χ₋₄` is the real number
`chi4RaceSum f L`. -/
theorem primeSum_chi4_ofReal (f : ℝ → ℝ) (L : ℝ) :
    primeSum chi4 (autocorr fun x => (f x : ℂ)) L = ((chi4RaceSum f L : ℝ) : ℂ) := by
  unfold primeSum chi4RaceSum
  have h : ∀ p : Nat.Primes,
      chi4 ((p : ℕ) : ZMod 4) * ((Real.log (p : ℕ) / Real.sqrt (p : ℕ) : ℝ) : ℂ)
        * autocorr (fun x => (f x : ℂ)) (Real.log (p : ℕ) - L)
      = ((((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * (Real.log (p : ℕ) / Real.sqrt (p : ℕ))
          * ∫ t, f t * f (t - (Real.log (p : ℕ) - L)) : ℝ) : ℂ) := by
    intro p
    rw [chi4_natCast, autocorr_ofReal]
    push_cast
    ring
  rw [tsum_congr h, ← Complex.ofReal_tsum]

theorem primeSum_chi4_im_eq_zero (f : ℝ → ℝ) (L : ℝ) :
    (primeSum chi4 (autocorr fun x => (f x : ℂ)) L).im = 0 := by
  rw [primeSum_chi4_ofReal, Complex.ofReal_im]

theorem contDiff_ofReal_comp {f : ℝ → ℝ} (hf : ContDiff ℝ 2 f) :
    ContDiff ℝ 2 (fun x => (f x : ℂ)) :=
  Complex.ofRealCLM.contDiff.comp hf

theorem hasCompactSupport_ofReal_comp {f : ℝ → ℝ} (hfc : HasCompactSupport f) :
    HasCompactSupport (fun x => (f x : ℂ)) :=
  hfc.comp_left Complex.ofReal_zero

theorem tsupport_ofReal_comp_subset (f : ℝ → ℝ) :
    tsupport (fun x => (f x : ℂ)) ⊆ tsupport f :=
  tsupport_comp_subset Complex.ofReal_zero f

/-- `ĝ(r) = |f̂(r)|²` for real `r` and `g = f ⋆ f̃`. -/
theorem raceFT_autocorr_ofReal {F : ℝ → ℂ} (hF : Continuous F) (hFc : HasCompactSupport F)
    (r : ℝ) : raceFT (autocorr F) r = ((‖raceFT F r‖ ^ 2 : ℝ) : ℂ) := by
  rw [autocorr_eq_weilTest, raceFT_eq_paperFT, paperFT_weilTest_self_ofReal hF hFc,
    Complex.normSq_eq_norm_sq]

/-- `Re Γ_0 = chi4Gamma0 f` for a real-valued window. -/
theorem archTerm_chi4_zero_ofReal {f : ℝ → ℝ} (hf : Continuous f) (hfc : HasCompactSupport f) :
    archTerm chi4 (autocorr fun x => (f x : ℂ)) 0 = ((chi4Gamma0 f : ℝ) : ℂ) := by
  have hF : Continuous (fun x => (f x : ℂ)) := Complex.continuous_ofReal.comp hf
  have hFc := hasCompactSupport_ofReal_comp hfc
  unfold archTerm chi4Gamma0
  have h : ∀ r : ℝ, cexp (I * r * ((0 : ℝ) : ℂ)) * raceFT (autocorr fun x => (f x : ℂ)) r
      * (archWeight chi4 r : ℂ)
      = ((‖raceFT (fun x => (f x : ℂ)) r‖ ^ 2
          * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi)) : ℝ) : ℂ) := by
    intro r
    rw [raceFT_autocorr_ofReal hF hFc, archWeight_chi4]
    simp
  simp only [h]
  rw [integral_complex_ofReal]
  push_cast
  ring

/-- The integrand of `chi4Gamma0 f` is integrable. -/
theorem integrable_chi4Gamma0_integrand {f : ℝ → ℝ} (hf : ContDiff ℝ 2 f)
    (hfc : HasCompactSupport f) :
    Integrable (fun r : ℝ => ‖raceFT (fun x => (f x : ℂ)) r‖ ^ 2
      * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi))) := by
  have hF2 := contDiff_ofReal_comp hf
  have hFc := hasCompactSupport_ofReal_comp hfc
  have hg2 : ContDiff ℝ 2 (autocorr fun x => (f x : ℂ)) := by
    rw [autocorr_eq_weilTest]; exact weilTest_contDiff hF2 hF2.continuous hFc
  have hgc : HasCompactSupport (autocorr fun x => (f x : ℂ)) := by
    rw [autocorr_eq_weilTest]; exact weilTest_hasCompactSupport hFc hFc
  have h := (integrable_raceFT_mul_archWeight chi4 hg2 hgc).re
  refine h.congr (Eventually.of_forall fun r => ?_)
  simp only [raceFT_autocorr_ofReal hF2.continuous hFc, archWeight_chi4, ← Complex.ofReal_mul,
    RCLike.re_to_complex, Complex.ofReal_re]

/-- `ĝ(0) = (∫ f)²` for a real-valued window. -/
theorem raceFT_autocorr_zero_ofReal {f : ℝ → ℝ} (hf : Continuous f) (hfc : HasCompactSupport f) :
    raceFT (autocorr fun x => (f x : ℂ)) 0 = (((∫ x, f x) ^ 2 : ℝ) : ℂ) := by
  have hF : Continuous (fun x => (f x : ℂ)) := Complex.continuous_ofReal.comp hf
  have hFc := hasCompactSupport_ofReal_comp hfc
  have h0 : raceFT (fun x => (f x : ℂ)) ((0 : ℝ) : ℂ) = ((∫ x, f x : ℝ) : ℂ) := by
    unfold raceFT
    rw [← integral_complex_ofReal]
    congr 1; ext u
    simp
  have h := raceFT_autocorr_ofReal hF hFc 0
  rw [h0, Complex.norm_real, Real.norm_eq_abs, sq_abs] at h
  rw [← Complex.ofReal_zero]
  exact h

/-- **S1 for a real-valued window, in real terms.** -/
theorem chi4RaceSum_lt (hGRH : GRHFor chi4) (hPNT : PNTTheta)
    {f : ℝ → ℝ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {a : ℝ}
    (ha : a < Real.log 3) (hsupp : tsupport f ⊆ Icc (-(a / 2)) (a / 2)) {margin : ℝ}
    (hcert : chi4Gamma0 f < (1 / 2) * (∫ x, f x) ^ 2 - margin) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum f L < -margin := by
  have h := chi4_primeSum_re_lt hGRH hPNT (contDiff_ofReal_comp hf)
    (hasCompactSupport_ofReal_comp hfc) ha ((tsupport_ofReal_comp_subset f).trans hsupp)
    (margin := margin) (by
      rw [archTerm_chi4_zero_ofReal hf.continuous hfc,
        raceFT_autocorr_zero_ofReal hf.continuous hfc, Complex.ofReal_re, Complex.ofReal_re]
      exact hcert)
  obtain ⟨L₀, hL₀⟩ := h
  refine ⟨L₀, fun L hL => ?_⟩
  have := hL₀ L hL
  rwa [primeSum_chi4_ofReal, Complex.ofReal_re] at this

/-- **Two-sided S1 for a real-valued window, in real terms.** -/
theorem abs_chi4RaceSum_add_le (hGRH : GRHFor chi4) (hPNT : PNTTheta)
    {f : ℝ → ℝ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {a : ℝ}
    (ha : a < Real.log 3) (hsupp : tsupport f ⊆ Icc (-(a / 2)) (a / 2)) {ε : ℝ} (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, |chi4RaceSum f L + (1 / 2) * (∫ x, f x) ^ 2| ≤ chi4Gamma0 f + ε := by
  obtain ⟨L₀, hL₀⟩ := chi4_norm_primeSum_add_le hGRH hPNT (contDiff_ofReal_comp hf)
    (hasCompactSupport_ofReal_comp hfc) ha ((tsupport_ofReal_comp_subset f).trans hsupp) hε
  refine ⟨L₀, fun L hL => ?_⟩
  have h := hL₀ L hL
  rw [archTerm_chi4_zero_ofReal hf.continuous hfc, raceFT_autocorr_zero_ofReal hf.continuous hfc,
    primeSum_chi4_ofReal, Complex.ofReal_re] at h
  have hc : ((chi4RaceSum f L : ℝ) : ℂ) + (1 / 2 : ℂ) * (((∫ x, f x) ^ 2 : ℝ) : ℂ)
      = ((chi4RaceSum f L + (1 / 2) * (∫ x, f x) ^ 2 : ℝ) : ℂ) := by
    push_cast; ring
  rwa [hc, Complex.norm_real, Real.norm_eq_abs] at h

end RaceLib
