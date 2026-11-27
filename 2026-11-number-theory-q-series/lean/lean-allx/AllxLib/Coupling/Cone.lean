/-
Copyright (c) 2026. Released under the repository's license.

# The cone `m ≤ 0, m < m' ≤ 1` (`EXPLICIT.md`, Theorem 5, deterministic part)

For a real weight `w`, `m ≤ 0`, `m < m' ≤ 1`, `κ = m' − m` and `x ≥ 2`: if
`Θ₁(t) = ∑_{p ≤ t} w(p) log p (t − p) ≥ 0` on `[2, x]`, then

  `D_m(x) − x^{−κ} D_{m'}(x) ≥ κ x^{m−1} Θ₁(x)/log x`        (`Coupling.cone_lower_bound`)

and hence (`Coupling.cone_lead_inheritance`) `Θ₁(x) > 0` and `D_{m'}(x) ≥ 0` give `D_m(x) > 0`.
In the note the hypothesis `Θ₁ > 0` on `(p₀, ∞)` is Theorem 1 (proved there under GRH);
here it is an explicit hypothesis.

Proof: Lemma 2 (`second_order_partial_summation`) with `F = f/log t`,
`f(t) = t^m − x^{−κ} t^{m'}`; `F(x) = 0`, `F'(x) = −κ x^{m−1}/log x`, and
`F'' = f''/L − 2f'/(tL²) + f(L + 2)/(t²L³) ≥ 0` (`L = log t`) because `f ≥ 0`, `f' ≤ 0`,
`f'' ≥ 0` on `[2, x]` in the cone.
-/
import AllxLib.Coupling.Convexity

open Real Finset MeasureTheory intervalIntegral

namespace Coupling

section Cone

variable (m m' x : ℝ)

/-- `f(t) = t^m − x^{−κ} t^{m'}`. -/
noncomputable def coneF0 (t : ℝ) : ℝ := t ^ m - x ^ (-(m' - m)) * t ^ m'
/-- `f'(t)`. -/
noncomputable def coneF1 (t : ℝ) : ℝ := m * t ^ (m - 1) - x ^ (-(m' - m)) * (m' * t ^ (m' - 1))
/-- `f''(t)`. -/
noncomputable def coneF2 (t : ℝ) : ℝ :=
  m * (m - 1) * t ^ (m - 2) - x ^ (-(m' - m)) * (m' * (m' - 1) * t ^ (m' - 2))

/-- `F = f/log t`. -/
noncomputable def coneG0 (t : ℝ) : ℝ := coneF0 m m' x t / log t
/-- `F' = f'/L − f/(t L²)`. -/
noncomputable def coneG1 (t : ℝ) : ℝ :=
  coneF1 m m' x t / log t - coneF0 m m' x t / (t * log t ^ 2)
/-- `F'' = f''/L − 2f'/(tL²) + f(L + 2)/(t²L³)`. -/
noncomputable def coneG2 (t : ℝ) : ℝ :=
  coneF2 m m' x t / log t - 2 * coneF1 m m' x t / (t * log t ^ 2) +
    coneF0 m m' x t * (log t + 2) / (t ^ 2 * log t ^ 3)

variable {m m' x}

lemma hasDerivAt_coneF0 {t : ℝ} (ht : 0 < t) :
    HasDerivAt (coneF0 m m' x) (coneF1 m m' x t) t := by
  unfold coneF0 coneF1
  exact (Real.hasDerivAt_rpow_const (Or.inl ht.ne')).sub
    ((Real.hasDerivAt_rpow_const (Or.inl ht.ne')).const_mul _)

lemma hasDerivAt_coneF1 {t : ℝ} (ht : 0 < t) :
    HasDerivAt (coneF1 m m' x) (coneF2 m m' x t) t := by
  unfold coneF1 coneF2
  have h1 := (Real.hasDerivAt_rpow_const (p := m - 1) (Or.inl ht.ne')).const_mul m
  have h2 := ((Real.hasDerivAt_rpow_const (p := m' - 1) (Or.inl ht.ne')).const_mul m').const_mul
    (x ^ (-(m' - m)))
  refine (h1.sub h2).congr_deriv ?_
  rw [show m - 1 - 1 = m - 2 by ring, show m' - 1 - 1 = m' - 2 by ring]
  ring

lemma hasDerivAt_coneG0 {t : ℝ} (ht : 1 < t) :
    HasDerivAt (coneG0 m m' x) (coneG1 m m' x t) t := by
  have ht0 : 0 < t := by linarith
  have hl : log t ≠ 0 := (Real.log_pos ht).ne'
  have h := (hasDerivAt_coneF0 (m := m) (m' := m') (x := x) ht0).div
    (Real.hasDerivAt_log ht0.ne') hl
  unfold coneG0 coneG1
  refine h.congr_deriv ?_
  field_simp

lemma hasDerivAt_coneG1 {t : ℝ} (ht : 1 < t) :
    HasDerivAt (coneG1 m m' x) (coneG2 m m' x t) t := by
  have ht0 : 0 < t := by linarith
  have hl : log t ≠ 0 := (Real.log_pos ht).ne'
  have hlog := Real.hasDerivAt_log ht0.ne'
  have ha := (hasDerivAt_coneF1 (m := m) (m' := m') (x := x) ht0).div hlog hl
  have hden := (hasDerivAt_id t).mul (hlog.pow 2)
  have hb := (hasDerivAt_coneF0 (m := m) (m' := m') (x := x) ht0).div hden
    (by simp only [Pi.mul_apply, Pi.pow_apply, id]; exact mul_ne_zero ht0.ne' (pow_ne_zero _ hl))
  unfold coneG1 coneG2
  refine (ha.sub hb).congr_deriv ?_
  simp only [Pi.mul_apply, Pi.pow_apply, id, Nat.cast_ofNat]
  field_simp
  ring

/-- Sign facts: on `0 < t ≤ x` in the cone, `f ≥ 0`, `f' ≤ 0`, `f'' ≥ 0`. -/
lemma cone_signs (hm : m ≤ 0) (hmm' : m < m') (hm'1 : m' ≤ 1) {t : ℝ} (ht : 0 < t)
    (htx : t ≤ x) :
    0 ≤ coneF0 m m' x t ∧ coneF1 m m' x t ≤ 0 ∧ 0 ≤ coneF2 m m' x t := by
  have hx : 0 < x := ht.trans_le htx
  set κ := m' - m with hκ
  have hκ0 : 0 < κ := by linarith
  set r := (t / x) ^ κ with hr
  have hr0 : 0 < r := Real.rpow_pos_of_pos (div_pos ht hx) _
  have hr1 : r ≤ 1 := Real.rpow_le_one (div_pos ht hx).le ((div_le_one hx).2 htx) hκ0.le
  -- `x^{−κ} t^{a + κ} = r t^a`
  have key : ∀ a : ℝ, x ^ (-κ) * t ^ (a + κ) = r * t ^ a := by
    intro a
    rw [Real.rpow_add ht, hr, Real.div_rpow ht.le hx.le, Real.rpow_neg hx.le]
    field_simp
  have e0 : x ^ (-κ) * t ^ m' = r * t ^ m := by rw [← key]; congr 2; rw [hκ]; ring
  have e1 : x ^ (-κ) * t ^ (m' - 1) = r * t ^ (m - 1) := by
    rw [← key]; congr 2; rw [hκ]; ring
  have e2 : x ^ (-κ) * t ^ (m' - 2) = r * t ^ (m - 2) := by
    rw [← key]; congr 2; rw [hκ]; ring
  have p0 : 0 < t ^ m := Real.rpow_pos_of_pos ht _
  have p1 : 0 < t ^ (m - 1) := Real.rpow_pos_of_pos ht _
  have p2 : 0 < t ^ (m - 2) := Real.rpow_pos_of_pos ht _
  refine ⟨?_, ?_, ?_⟩
  · unfold coneF0
    rw [e0]; nlinarith
  · unfold coneF1
    have : x ^ (-κ) * (m' * t ^ (m' - 1)) = m' * r * t ^ (m - 1) := by
      rw [mul_left_comm, e1]; ring
    rw [this]
    have hsign : m - m' * r ≤ 0 := by
      rcases le_or_gt 0 m' with h | h
      · nlinarith
      · nlinarith
    nlinarith
  · unfold coneF2
    have : x ^ (-κ) * (m' * (m' - 1) * t ^ (m' - 2)) = m' * (m' - 1) * r * t ^ (m - 2) := by
      rw [mul_left_comm, e2]; ring
    rw [this]
    have hsign : 0 ≤ m * (m - 1) - m' * (m' - 1) * r := by
      rcases le_or_gt 0 m' with h | h
      · have : m' * (m' - 1) ≤ 0 := by nlinarith
        nlinarith
      · have h1 : m' * (m' - 1) < m * (m - 1) := by nlinarith
        have h2 : 0 < m' * (m' - 1) := by nlinarith
        nlinarith
    nlinarith

lemma coneG2_nonneg (hm : m ≤ 0) (hmm' : m < m') (hm'1 : m' ≤ 1) {t : ℝ} (ht : 1 < t)
    (htx : t ≤ x) : 0 ≤ coneG2 m m' x t := by
  obtain ⟨h0, h1, h2⟩ := cone_signs hm hmm' hm'1 (by linarith) htx
  have hl : 0 < log t := Real.log_pos ht
  have ht0 : 0 < t := by linarith
  unfold coneG2
  have a : 0 ≤ coneF2 m m' x t / log t := div_nonneg h2 hl.le
  have b : 0 ≤ -(2 * coneF1 m m' x t / (t * log t ^ 2)) := by
    rw [← neg_div]; exact div_nonneg (by linarith) (by positivity)
  have c : 0 ≤ coneF0 m m' x t * (log t + 2) / (t ^ 2 * log t ^ 3) := by positivity
  linarith

lemma continuousOn_coneG2 : ContinuousOn (coneG2 m m' x) (Set.Icc 2 x) := by
  intro t ht
  have ht0 : t ≠ 0 := by linarith [ht.1]
  have hl : log t ≠ 0 := (Real.log_pos (by linarith [ht.1])).ne'
  apply ContinuousAt.continuousWithinAt
  unfold coneG2 coneF0 coneF1 coneF2
  have c1 : ∀ p : ℝ, ContinuousAt (fun t : ℝ => t ^ p) t :=
    fun p => Real.continuousAt_rpow_const _ _ (Or.inl ht0)
  have cl : ContinuousAt Real.log t := Real.continuousAt_log ht0
  fun_prop (disch := simp [*])

end Cone

/-- **`EXPLICIT.md`, Theorem 5 (deterministic part).** For a real weight `w`, `m ≤ 0`,
`m < m' ≤ 1`, `x ≥ 2`: if `Θ₁(t) = ∑_{p ≤ t} w(p) log p (t − p) ≥ 0` for `t ∈ [2, x]`, then
`D_m(x) − x^{−(m'−m)} D_{m'}(x) ≥ (m' − m) x^{m−1} Θ₁(x)/log x`. -/
theorem cone_lower_bound (w : ℕ → ℝ) {m m' x : ℝ} (hm : m ≤ 0) (hmm' : m < m') (hm'1 : m' ≤ 1)
    (hx : 2 ≤ x) (hnn : ∀ t ∈ Set.Icc 2 x, 0 ≤ Theta1 (fun p => w p * log p) t) :
    (m' - m) * x ^ (m - 1) * Theta1 (fun p => w p * log p) x / log x ≤
      primeSum w m x - x ^ (-(m' - m)) * primeSum w m' x := by
  have hF : ∀ t ∈ Set.Icc 2 x, HasDerivAt (coneG0 m m' x) (coneG1 m m' x t) t :=
    fun t ht => hasDerivAt_coneG0 (by linarith [ht.1])
  have hF' : ∀ t ∈ Set.Icc 2 x, HasDerivAt (coneG1 m m' x) (coneG2 m m' x t) t :=
    fun t ht => hasDerivAt_coneG1 (by linarith [ht.1])
  have key := second_order_partial_summation (fun p => w p * log p) hx hF hF'
    continuousOn_coneG2
  have hx0 : 0 < x := by linarith
  have hlx : 0 < log x := Real.log_pos (by linarith)
  -- left side
  have hlhs : ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, w p * log p * coneG0 m m' x p =
      primeSum w m x - x ^ (-(m' - m)) * primeSum w m' x := by
    rw [primeSum, primeSum, mul_sum, ← sum_sub_distrib]
    refine sum_congr rfl fun p hp => ?_
    have hp2 : (2 : ℝ) ≤ p := by exact_mod_cast (mem_filter.1 hp).2.two_le
    have hlp : log p ≠ 0 := (Real.log_pos (by linarith)).ne'
    unfold coneG0 coneF0
    field_simp
  -- boundary terms
  have hf0 : coneF0 m m' x x = 0 := by
    unfold coneF0
    rw [← Real.rpow_add hx0]; ring_nf
  have hG0 : coneG0 m m' x x = 0 := by unfold coneG0; rw [hf0, zero_div]
  have hG1 : coneG1 m m' x x = -((m' - m) * x ^ (m - 1)) / log x := by
    unfold coneG1 coneF1
    rw [hf0, zero_div, sub_zero]
    have : x ^ (-(m' - m)) * (m' * x ^ (m' - 1)) = m' * x ^ (m - 1) := by
      rw [mul_left_comm, ← Real.rpow_add hx0]; ring_nf
    rw [this]
    congr 1; ring
  rw [hlhs] at key
  rw [key, hG0, hG1, zero_mul, zero_sub]
  have hI : 0 ≤ ∫ t in (2 : ℝ)..x, Theta1 (fun p => w p * log p) t * coneG2 m m' x t :=
    intervalIntegral.integral_nonneg hx fun t ht =>
      mul_nonneg (hnn t ht) (coneG2_nonneg hm hmm' hm'1 (by linarith [ht.1]) ht.2)
  have : -(-((m' - m) * x ^ (m - 1)) / log x * Theta1 (fun p => w p * log p) x) =
      (m' - m) * x ^ (m - 1) * Theta1 (fun p => w p * log p) x / log x := by ring
  linarith

/-- **Lead inheritance for every `x` (`EXPLICIT.md`, Theorem 5 / Corollary 4 shape).** In the
cone `m ≤ 0, m < m' ≤ 1`, if `Θ₁ ≥ 0` on `[2, x]`, `Θ₁(x) > 0` and `D_{m'}(x) ≥ 0`, then
`D_m(x) > 0`. -/
theorem cone_lead_inheritance (w : ℕ → ℝ) {m m' x : ℝ} (hm : m ≤ 0) (hmm' : m < m')
    (hm'1 : m' ≤ 1) (hx : 2 ≤ x)
    (hnn : ∀ t ∈ Set.Icc 2 x, 0 ≤ Theta1 (fun p => w p * log p) t)
    (hpos : 0 < Theta1 (fun p => w p * log p) x) (hD' : 0 ≤ primeSum w m' x) :
    0 < primeSum w m x := by
  have h := cone_lower_bound w hm hmm' hm'1 hx hnn
  have hx0 : 0 < x := by linarith
  have hlx : 0 < log x := Real.log_pos (by linarith)
  have hb : 0 < (m' - m) * x ^ (m - 1) * Theta1 (fun p => w p * log p) x / log x := by
    have : 0 < m' - m := by linarith
    have := Real.rpow_pos_of_pos hx0 (m - 1)
    positivity
  have := mul_nonneg (Real.rpow_pos_of_pos hx0 (-(m' - m))).le hD'
  linarith

end Coupling
