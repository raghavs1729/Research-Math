/-
ChebAllxLib/Integral.lean — second-order partial summation over the INFINITE range, in integral
form (Mathlib only; no number theory): Lemma 1 of `monotonicity/chebyshev1853/CHEBYSHEV1853.md`
for convex test functions.

Let `a ≥ 0` and `F : ℝ → ℝ` be twice differentiable on `[a, ∞)` with `F'' ≥ 0` continuous there and
`F(t) → 0` as `t → ∞`.  Then

* `F' ≤ 0`, `F ≥ 0` on `[a, ∞)` and `(T − p) F'(T) → 0` (`deriv_nonpos_of_tendsto_zero`,
  `tendsto_mul_deriv`): the boundary terms of the partial summation vanish by convexity alone;
* `F(p) = ∫_p^∞ (t − p) F''(t) dt` for every `p ≥ a` (`integral_Ioi_sub_mul_deriv2`);
* for every real sequence `c` with `c(n) = 0` for `n < a` and `Σ_n |c(n)| F(n) < ∞`:

      Σ_n c(n) F(n) = ∫_a^∞ C₁(t) F''(t) dt,     C₁(t) = Σ_{n ≤ t} c(n) (t − n),

  the integrand being integrable (`tsum_mul_eq_integral`).  The exchange of sum and integral is
  Tonelli's theorem for the kernels `c(n) (t − n)₊ F''(t)`; no growth hypothesis on `c` or `C₁`
  is needed.

Also: the closed form of `F''` for `F(t) = e^{−t/x}/log t` and its positivity on `(1, ∞)`.
-/
import Mathlib

open Filter Topology MeasureTheory Set

namespace ChebAllx

section Kernel

variable {a : ℝ} {F F' F'' : ℝ → ℝ}

/-- Taylor's formula with integral remainder:
`∫_p^T (t − p) F''(t) dt = F(p) − F(T) + F'(T) (T − p)`. -/
lemma integral_sub_mul_deriv2 (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF'' : ContinuousOn F'' (Ici a))
    {p T : ℝ} (hap : a ≤ p) (hpT : p ≤ T) :
    ∫ t in p..T, (t - p) * F'' t = F p - F T + F' T * (T - p) := by
  have hg : ∀ t ∈ Set.uIcc p T,
      HasDerivAt (fun t => F t + F' t * (p - t)) ((p - t) * F'' t) t := by
    intro t ht
    rw [Set.uIcc_of_le hpT] at ht
    have hat : a ≤ t := hap.trans ht.1
    have h2 := (hF' t hat).mul ((hasDerivAt_const t p).sub (hasDerivAt_id t))
    refine ((hF t hat).add h2).congr_deriv ?_
    simp only [Pi.sub_apply, id]
    ring
  have hcont : ContinuousOn (fun t => (p - t) * F'' t) (Set.uIcc p T) := by
    rw [Set.uIcc_of_le hpT]
    exact (continuousOn_const.sub continuousOn_id).mul (hF''.mono fun t ht => hap.trans ht.1)
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt hg hcont.intervalIntegrable
  have e : ∫ t in p..T, (t - p) * F'' t = -∫ t in p..T, (p - t) * F'' t := by
    rw [← intervalIntegral.integral_neg]
    refine intervalIntegral.integral_congr fun t _ => ?_
    ring
  rw [e, hFTC]
  ring

lemma deriv_monotoneOn (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t)
    (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t) : MonotoneOn F' (Ici a) := by
  refine monotoneOn_of_deriv_nonneg (convex_Ici a)
    (fun t ht => (hF' t ht).continuousAt.continuousWithinAt) (fun t ht => ?_) (fun t ht => ?_)
  · rw [interior_Ici] at ht
    exact (hF' t ht.le).differentiableAt.differentiableWithinAt
  · rw [interior_Ici] at ht
    rw [(hF' t ht.le).deriv]
    exact hF''0 t ht.le

/-- A convex function tending to `0` at infinity has nonpositive derivative. -/
lemma deriv_nonpos_of_tendsto_zero (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Tendsto F atTop (𝓝 0)) {t₀ : ℝ} (ht₀ : a ≤ t₀) : F' t₀ ≤ 0 := by
  by_contra hpos
  push Not at hpos
  have hmono := deriv_monotoneOn hF' hF''0
  have hge : ∀ t, t₀ ≤ t → F t₀ + F' t₀ * (t - t₀) ≤ F t := by
    intro t ht
    have h := (convex_Ici t₀).mul_sub_le_image_sub_of_le_deriv (f := F) (C := F' t₀)
      (fun s hs => (hF s (ht₀.trans hs)).continuousAt.continuousWithinAt)
      (fun s hs => by
        rw [interior_Ici] at hs
        exact (hF s (ht₀.trans hs.le)).differentiableAt.differentiableWithinAt)
      (fun s hs => by
        rw [interior_Ici] at hs
        rw [(hF s (ht₀.trans hs.le)).deriv]
        exact hmono ht₀ (ht₀.trans hs.le) hs.le)
      t₀ self_mem_Ici t ht ht
    linarith
  have htop : Tendsto (fun t => F t₀ + F' t₀ * (t - t₀)) atTop atTop :=
    tendsto_atTop_add_const_left _ _
      ((tendsto_id.atTop_add tendsto_const_nhds).const_mul_atTop hpos)
  have hFtop : Tendsto F atTop atTop :=
    tendsto_atTop_mono' _ ((eventually_ge_atTop t₀).mono hge) htop
  exact not_tendsto_nhds_of_tendsto_atTop hFtop 0 hlim

lemma antitoneOn_of_tendsto_zero (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Tendsto F atTop (𝓝 0)) : AntitoneOn F (Ici a) := by
  refine antitoneOn_of_deriv_nonpos (convex_Ici a)
    (fun t ht => (hF t ht).continuousAt.continuousWithinAt) (fun t ht => ?_) (fun t ht => ?_)
  · rw [interior_Ici] at ht
    exact (hF t ht.le).differentiableAt.differentiableWithinAt
  · rw [interior_Ici] at ht
    rw [(hF t ht.le).deriv]
    exact deriv_nonpos_of_tendsto_zero hF hF' hF''0 hlim ht.le

lemma nonneg_of_tendsto_zero (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Tendsto F atTop (𝓝 0)) {t : ℝ} (ht : a ≤ t) : 0 ≤ F t := by
  have hA := antitoneOn_of_tendsto_zero hF hF' hF''0 hlim
  exact le_of_tendsto hlim ((eventually_ge_atTop t).mono fun s hs => hA ht (ht.trans hs) hs)

/-- For a convex function tending to `0`: `0 ≤ −F'(T) (T − p) ≤ 2 F((T+p)/2)`. -/
lemma neg_deriv_mul_le (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Tendsto F atTop (𝓝 0)) {p T : ℝ} (hap : a ≤ p) (hpT : p ≤ T) :
    -(F' T * (T - p)) ≤ 2 * F ((T + p) / 2) := by
  have hmono := deriv_monotoneOn hF' hF''0
  set m := (T + p) / 2 with hm
  have hpm : p ≤ m := by rw [hm]; linarith
  have hmT : m ≤ T := by rw [hm]; linarith
  have h := (convex_Icc m T).image_sub_le_mul_sub_of_deriv_le (f := F) (C := F' T)
    (fun s hs => (hF s (hap.trans (hpm.trans hs.1))).continuousAt.continuousWithinAt)
    (fun s hs => by
      rw [interior_Icc] at hs
      exact (hF s (hap.trans (hpm.trans hs.1.le))).differentiableAt.differentiableWithinAt)
    (fun s hs => by
      rw [interior_Icc] at hs
      rw [(hF s (hap.trans (hpm.trans hs.1.le))).deriv]
      exact hmono (hap.trans (hpm.trans hs.1.le)) (hap.trans hpT) hs.2.le)
    m (left_mem_Icc.2 hmT) T (right_mem_Icc.2 hmT) hmT
  have hFT := nonneg_of_tendsto_zero hF hF' hF''0 hlim (hap.trans hpT)
  have e : T - m = (T - p) / 2 := by rw [hm]; ring
  rw [e] at h
  linarith

/-- For a convex function tending to `0`: `F'(T) (T − p) → 0` as `T → ∞`. -/
lemma tendsto_mul_deriv (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Tendsto F atTop (𝓝 0)) {p : ℝ} (hap : a ≤ p) :
    Tendsto (fun T => F' T * (T - p)) atTop (𝓝 0) := by
  have hhalf : Tendsto (fun T : ℝ => (T + p) / 2) atTop atTop :=
    (tendsto_id.atTop_add tendsto_const_nhds).atTop_div_const (by norm_num)
  have hup : Tendsto (fun T : ℝ => -(2 * F ((T + p) / 2))) atTop (𝓝 0) := by
    simpa using ((hlim.comp hhalf).const_mul 2).neg
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hup tendsto_const_nhds ?_ ?_
  · filter_upwards [eventually_ge_atTop p] with T hT
    have := neg_deriv_mul_le hF hF' hF''0 hlim hap hT
    linarith
  · filter_upwards [eventually_ge_atTop p] with T hT
    exact mul_nonpos_of_nonpos_of_nonneg
      (deriv_nonpos_of_tendsto_zero hF hF' hF''0 hlim (hap.trans hT)) (sub_nonneg.2 hT)

/-- **The kernel representation.**  For a convex `F` on `[a, ∞)` tending to `0` and every `p ≥ a`:
`F(p) = ∫_p^∞ (t − p) F''(t) dt`. -/
theorem integral_Ioi_sub_mul_deriv2 (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF'' : ContinuousOn F'' (Ici a))
    (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t) (hlim : Tendsto F atTop (𝓝 0)) {p : ℝ} (hap : a ≤ p) :
    IntegrableOn (fun t => (t - p) * F'' t) (Ioi p) ∧ ∫ t in Ioi p, (t - p) * F'' t = F p := by
  have hT : Tendsto (fun T => ∫ t in p..T, (t - p) * F'' t) atTop (𝓝 (F p)) := by
    have h1 : Tendsto (fun T => F p - F T + F' T * (T - p)) atTop (𝓝 (F p)) := by
      simpa using (tendsto_const_nhds.sub hlim).add (tendsto_mul_deriv hF hF' hF''0 hlim hap)
    refine h1.congr' ?_
    filter_upwards [eventually_ge_atTop p] with T hT
    exact (integral_sub_mul_deriv2 hF hF' hF'' hap hT).symm
  have hcont : ContinuousOn (fun t => (t - p) * F'' t) (Ici p) :=
    (continuousOn_id.sub continuousOn_const).mul (hF''.mono fun t ht => hap.trans ht)
  have hint : IntegrableOn (fun t => (t - p) * F'' t) (Ioi p) := by
    refine integrableOn_Ioi_of_intervalIntegral_norm_tendsto (F p) p (b := id) (l := atTop)
      (fun T => ?_) tendsto_id ?_
    · exact (hcont.mono fun t ht => ht.1).integrableOn_Icc.mono_set Ioc_subset_Icc_self
    · refine hT.congr' ?_
      filter_upwards [eventually_ge_atTop p] with T hT
      refine intervalIntegral.integral_congr_ae ?_
      refine Filter.Eventually.of_forall fun t ht => ?_
      rw [Set.uIoc_of_le hT] at ht
      rw [Real.norm_of_nonneg (mul_nonneg (sub_nonneg.2 ht.1.le) (hF''0 t (hap.trans ht.1.le)))]
  exact ⟨hint, tendsto_nhds_unique (intervalIntegral_tendsto_integral_Ioi p hint tendsto_id) hT⟩

/-- **Second-order partial summation over the infinite range (Lemma 1 of the note, for convex
test functions).**  Let `a ≥ 0`, `F` twice differentiable on `[a, ∞)` with `F'' ≥ 0` continuous
and `F(t) → 0`.  For every real sequence `c` vanishing below `a` with `Σ |c(n)| F(n) < ∞`:

    Σ_n c(n) F(n) = ∫_a^∞ (Σ_{n ≤ t} c(n) (t − n)) F''(t) dt,

and the integrand is integrable on `(a, ∞)`. -/
theorem tsum_mul_eq_integral (ha : 0 ≤ a) {c : ℕ → ℝ} (hc : ∀ n : ℕ, (n : ℝ) < a → c n = 0)
    (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t)
    (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t) (hF'' : ContinuousOn F'' (Ici a))
    (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t) (hlim : Tendsto F atTop (𝓝 0))
    (hs : Summable fun n : ℕ => |c n| * F n) :
    IntegrableOn (fun t => (∑ n ∈ Finset.Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t) (Ioi a) ∧
      ∑' n, c n * F n = ∫ t in Ioi a, (∑ n ∈ Finset.Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t := by
  -- the kernels
  set G : ℕ → ℝ → ℝ := fun n t => c n * (Ioi (n : ℝ)).indicator (fun t => (t - n) * F'' t) t
    with hG
  have hKint : ∀ n : ℕ, a ≤ n →
      Integrable ((Ioi (n : ℝ)).indicator (fun t => (t - n) * F'' t)) := fun n hn =>
    (integral_Ioi_sub_mul_deriv2 hF hF' hF'' hF''0 hlim hn).1.integrable_indicator
      measurableSet_Ioi
  have hKval : ∀ n : ℕ, a ≤ n →
      ∫ t, (Ioi (n : ℝ)).indicator (fun t => (t - n) * F'' t) t = F n := fun n hn => by
    rw [integral_indicator measurableSet_Ioi]
    exact (integral_Ioi_sub_mul_deriv2 hF hF' hF'' hF''0 hlim hn).2
  have hKnn : ∀ n : ℕ, a ≤ n → ∀ t, 0 ≤ (Ioi (n : ℝ)).indicator (fun t => (t - n) * F'' t) t :=
    fun n hn t => Set.indicator_nonneg
      (fun s hs => mul_nonneg (sub_nonneg.2 (le_of_lt (mem_Ioi.1 hs)))
        (hF''0 s (hn.trans (le_of_lt (mem_Ioi.1 hs))))) t
  have hGint : ∀ n, Integrable (G n) := by
    intro n
    by_cases hn : a ≤ n
    · exact (hKint n hn).const_mul (c n)
    · have : G n = fun _ => 0 := by
        funext t
        simp only [hG, hc n (not_le.1 hn), zero_mul]
      rw [this]
      exact integrable_zero _ _ _
  have hGnorm : ∀ n, ∫ t, ‖G n t‖ = |c n| * F n := by
    intro n
    by_cases hn : a ≤ n
    · have : (fun t => ‖G n t‖) = fun t =>
          |c n| * (Ioi (n : ℝ)).indicator (fun t => (t - n) * F'' t) t := by
        funext t
        simp only [hG, norm_mul, Real.norm_eq_abs, abs_of_nonneg (hKnn n hn t)]
      rw [this, integral_const_mul, hKval n hn]
    · simp only [hG, hc n (not_le.1 hn), zero_mul, norm_zero, integral_zero, abs_zero]
  have hGval : ∀ n, ∫ t, G n t = c n * F n := by
    intro n
    by_cases hn : a ≤ n
    · simp only [hG]
      rw [integral_const_mul, hKval n hn]
    · simp only [hG, hc n (not_le.1 hn), zero_mul, integral_zero]
  have hsum : Summable fun n => ∫ t, ‖G n t‖ := by
    simpa only [hGnorm] using hs
  -- the pointwise sum of the kernels
  have hpt : ∀ t, ∑' n, G n t
      = (Ioi a).indicator (fun t => (∑ n ∈ Finset.Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t) t := by
    intro t
    by_cases hta : a < t
    · rw [Set.indicator_of_mem (mem_Ioi.2 hta)]
      have ht0 : 0 ≤ t := ha.trans hta.le
      rw [tsum_eq_sum (s := Finset.Icc 0 ⌊t⌋₊), Finset.sum_mul]
      · refine Finset.sum_congr rfl fun n hn => ?_
        have hnt : (n : ℝ) ≤ t := (Nat.le_floor_iff ht0).1 (Finset.mem_Icc.1 hn).2
        simp only [hG]
        rcases hnt.lt_or_eq with h | h
        · rw [Set.indicator_of_mem (mem_Ioi.2 h)]
          ring
        · rw [Set.indicator_of_notMem (by rw [mem_Ioi, h]; exact lt_irrefl _), h]
          ring
      · intro n hn
        have hnt : t < n := by
          by_contra hle
          push Not at hle
          exact hn (Finset.mem_Icc.2 ⟨Nat.zero_le _, Nat.le_floor hle⟩)
        simp only [hG]
        rw [Set.indicator_of_notMem (by rw [mem_Ioi]; exact not_lt.2 hnt.le), mul_zero]
    · rw [Set.indicator_of_notMem (by rwa [mem_Ioi])]
      have : ∀ n, G n t = 0 := by
        intro n
        simp only [hG]
        by_cases hn : a ≤ n
        · rw [Set.indicator_of_notMem (by
            rw [mem_Ioi]; exact not_lt.2 ((not_lt.1 hta).trans hn)), mul_zero]
        · rw [hc n (not_le.1 hn), zero_mul]
      simp only [this, tsum_zero]
  have hHS := hasSum_integral_of_summable_integral_norm hGint hsum
  have hval : ∑' n, c n * F n
      = ∫ t in Ioi a, (∑ n ∈ Finset.Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t := by
    have h := hHS.tsum_eq
    simp only [hGval] at h
    rw [h, ← integral_indicator measurableSet_Ioi]
    exact integral_congr_ae (Filter.Eventually.of_forall hpt)
  refine ⟨?_, hval⟩
  -- integrability: the sum of the kernels is integrable
  have hzero : ∀ t, ∀ n ∉ Finset.Icc 0 ⌊t⌋₊, G n t = 0 := by
    intro t n hn
    simp only [hG]
    by_cases hna : a ≤ n
    · have hnt : t ≤ n := by
        by_contra hlt
        push Not at hlt
        exact hn (Finset.mem_Icc.2 ⟨Nat.zero_le _, Nat.le_floor hlt.le⟩)
      rw [Set.indicator_of_notMem (by rw [mem_Ioi]; exact not_lt.2 hnt), mul_zero]
    · rw [hc n (not_le.1 hna), zero_mul]
  have hI : Integrable (fun t => ∑' n, G n t) := by
    have hmeas : AEStronglyMeasurable (fun t => ∑' n, G n t) := by
      refine aestronglyMeasurable_of_tendsto_ae atTop
        (f := fun N t => ∑ n ∈ Finset.range N, G n t) (fun N => ?_) ?_
      · exact (integrable_finsetSum _ fun n _ => hGint n).aestronglyMeasurable
      · exact Filter.Eventually.of_forall fun t =>
          (summable_of_ne_finset_zero (hzero t)).hasSum.tendsto_sum_nat
    refine ⟨hmeas, ?_⟩
    rw [hasFiniteIntegral_iff_enorm]
    calc ∫⁻ t, ‖∑' n, G n t‖ₑ ≤ ∫⁻ t, ∑' n, ‖G n t‖ₑ := by
          refine lintegral_mono fun t => ?_
          rw [tsum_eq_sum (hzero t)]
          exact (enorm_sum_le _ _).trans (ENNReal.sum_le_tsum _)
      _ = ∑' n, ∫⁻ t, ‖G n t‖ₑ :=
          lintegral_tsum fun n => (hGint n).aestronglyMeasurable.enorm
      _ = ∑' n, ENNReal.ofReal (|c n| * F n) := by
          refine tsum_congr fun n => ?_
          rw [← hGnorm n, ofReal_integral_norm_eq_lintegral_enorm (hGint n)]
      _ = ENNReal.ofReal (∑' n, |c n| * F n) :=
          (ENNReal.ofReal_tsum_of_nonneg (fun n => by
            rw [← hGnorm n]; exact integral_nonneg fun _ => norm_nonneg _) hs).symm
      _ < ⊤ := ENNReal.ofReal_lt_top
  rw [← integrable_indicator_iff measurableSet_Ioi]
  exact hI.congr (Filter.Eventually.of_forall hpt)

end Kernel

/-! ### The test function `F(t) = e^{−t/x}/log t` -/

section ExpLog

/-- `F(t) = e^{−t/x}/log t`. -/
noncomputable def chebF (x t : ℝ) : ℝ := Real.exp (-t / x) / Real.log t

/-- `F'(t) = −e^{−t/x} (1/(x log t) + 1/(t log² t))`. -/
noncomputable def chebF' (x t : ℝ) : ℝ :=
  -Real.exp (-t / x) * (1 / (x * Real.log t) + 1 / (t * Real.log t ^ 2))

/-- `F''(t) = e^{−t/x} (1/(x² log t) + 2/(x t log² t) + 1/(t² log² t) + 2/(t² log³ t))`. -/
noncomputable def chebF'' (x t : ℝ) : ℝ :=
  Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
    + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3))

variable {x t : ℝ}

lemma hasDerivAt_exp_neg_div (x t : ℝ) :
    HasDerivAt (fun t : ℝ => Real.exp (-t / x)) (Real.exp (-t / x) * (-1 / x)) t := by
  have h : HasDerivAt (fun t : ℝ => -t / x) (-1 / x) t := (hasDerivAt_id t).neg.div_const x
  exact h.exp

lemma hasDerivAt_chebF (hx : x ≠ 0) (ht : 1 < t) : HasDerivAt (chebF x) (chebF' x t) t := by
  have ht0 : t ≠ 0 := by linarith
  have hl : Real.log t ≠ 0 := (Real.log_pos ht).ne'
  have h := (hasDerivAt_exp_neg_div x t).div (Real.hasDerivAt_log ht0) hl
  refine h.congr_deriv ?_
  unfold chebF'
  field_simp
  ring

lemma hasDerivAt_chebF' (hx : x ≠ 0) (ht : 1 < t) : HasDerivAt (chebF' x) (chebF'' x t) t := by
  have ht0 : t ≠ 0 := by linarith
  have hl : Real.log t ≠ 0 := (Real.log_pos ht).ne'
  have hL := Real.hasDerivAt_log ht0
  have h1 : HasDerivAt (fun t : ℝ => 1 / (x * Real.log t))
      ((0 * (x * Real.log t) - 1 * (x * t⁻¹)) / (x * Real.log t) ^ 2) t :=
    (hasDerivAt_const t (1 : ℝ)).div (hL.const_mul x) (mul_ne_zero hx hl)
  have h2 : HasDerivAt (fun t : ℝ => 1 / (t * Real.log t ^ 2))
      ((0 * (t * Real.log t ^ 2) - 1 * (1 * Real.log t ^ 2 + t * (2 * Real.log t ^ 1 * t⁻¹)))
        / (t * Real.log t ^ 2) ^ 2) t :=
    (hasDerivAt_const t (1 : ℝ)).div ((hasDerivAt_id t).mul (hL.pow 2))
      (mul_ne_zero ht0 (pow_ne_zero 2 hl))
  have h := (hasDerivAt_exp_neg_div x t).neg.mul (h1.add h2)
  refine h.congr_deriv ?_
  unfold chebF''
  simp only [Pi.neg_apply, Pi.add_apply]
  field_simp
  ring

lemma chebF''_pos (hx : 0 < x) (ht : 1 < t) : 0 < chebF'' x t := by
  have ht0 : 0 < t := by linarith
  have hl : 0 < Real.log t := Real.log_pos ht
  unfold chebF''
  positivity

lemma continuousOn_chebF'' (hx : x ≠ 0) : ContinuousOn (chebF'' x) (Ici 2) := by
  have hl : ∀ t ∈ Ici (2 : ℝ), Real.log t ≠ 0 := fun t ht =>
    (Real.log_pos (by linarith [mem_Ici.1 ht])).ne'
  have ht0 : ∀ t ∈ Ici (2 : ℝ), t ≠ 0 := fun t ht => by linarith [mem_Ici.1 ht]
  have hlog : ContinuousOn Real.log (Ici (2 : ℝ)) :=
    Real.continuousOn_log.mono fun t ht => ht0 t ht
  have hid : ContinuousOn (fun t : ℝ => t) (Ici (2 : ℝ)) := continuousOn_id
  unfold chebF''
  refine ContinuousOn.mul (Real.continuous_exp.comp_continuousOn
    ((continuousOn_id.neg).div_const x)) ?_
  refine ContinuousOn.add (ContinuousOn.add (ContinuousOn.add ?_ ?_) ?_) ?_
  · exact continuousOn_const.div (continuousOn_const.mul hlog)
      fun t ht => mul_ne_zero (pow_ne_zero 2 hx) (hl t ht)
  · exact continuousOn_const.div ((continuousOn_const.mul hid).mul (hlog.pow 2))
      fun t ht => mul_ne_zero (mul_ne_zero hx (ht0 t ht)) (pow_ne_zero 2 (hl t ht))
  · exact continuousOn_const.div ((hid.pow 2).mul (hlog.pow 2))
      fun t ht => mul_ne_zero (pow_ne_zero 2 (ht0 t ht)) (pow_ne_zero 2 (hl t ht))
  · exact continuousOn_const.div ((hid.pow 2).mul (hlog.pow 3))
      fun t ht => mul_ne_zero (pow_ne_zero 2 (ht0 t ht)) (pow_ne_zero 3 (hl t ht))

lemma tendsto_chebF (hx : 0 < x) : Tendsto (chebF x) atTop (𝓝 0) := by
  have h1 : Tendsto (fun t : ℝ => Real.exp (-t / x)) atTop (𝓝 0) :=
    Real.tendsto_exp_atBot.comp ((tendsto_neg_atTop_atBot).atBot_div_const hx)
  exact h1.div_atTop Real.tendsto_log_atTop

end ExpLog

/-- A function `e^{−g}` with `g` concave is convex. -/
lemma convexOn_exp_neg_of_concaveOn {s : Set ℝ} {g : ℝ → ℝ} (hg : ConcaveOn ℝ s g) :
    ConvexOn ℝ s (fun t => Real.exp (-g t)) := by
  refine ⟨hg.1, fun u hu v hv a b ha hb hab => ?_⟩
  have h1 := hg.2 hu hv ha hb hab
  have h2 := convexOn_exp.2 (mem_univ (-g u)) (mem_univ (-g v)) ha hb hab
  simp only [smul_eq_mul] at h1 h2 ⊢
  refine le_trans (Real.exp_le_exp.2 ?_) h2
  linarith

end ChebAllx
