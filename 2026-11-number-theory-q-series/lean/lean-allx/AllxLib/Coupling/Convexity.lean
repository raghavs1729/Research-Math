/-
Copyright (c) 2026. Released under the repository's license.

# The exact convexity identity (Target B; `explicit/EXPLICIT.md`, Lemma 2 and (2.1))

For a real weight `c : ℕ → ℝ` on the primes put (sums over primes `p`)

* `Θ c t  = ∑_{p ≤ t} c(p)`,
* `Θ₁ c t = ∑_{p ≤ t} c(p) (t − p)`  (`= ∫_0^t Θ c`; continuous, piecewise linear).

**Lemma 2 (second-order partial summation; no hypothesis).** If `x ≥ 2` and `F` is twice
differentiable on `[2, x]` with continuous `F''`, then

  `∑_{p ≤ x} c(p) F(p) = F(x) Θ(x) − F'(x) Θ₁(x) + ∫_2^x Θ₁(t) F''(t) dt`.

(The note integrates from `p₀`, the least prime not dividing `q`; since `Θ₁ = 0` on
`[2, p₀]` the two integrals agree.) The proof here is Taylor's formula with integral
remainder applied to each prime, followed by an exchange of a finite sum and an integral;
it uses no Stieltjes integration.

With `c(p) = w(p) log p` and `h(t) = (x − t)/log t` (convex on `(1, x]`):

  `x R(x) = Θ₁(x)/log x + ∫_2^x Θ₁(t) h''(t) dt`,   `h'' > 0` on `(1, x]`        (2.1)

(`Coupling.cesaro_convexity_identity`), hence **if `Θ₁ ≥ 0` on `[2, x]` and `Θ₁(x) > 0`
then `R(x) > 0`**, indeed `R(x) ≥ Θ₁(x)/(x log x)` (`Coupling.cesaro_pos_of_theta1`).
-/
import AllxLib.Coupling.Abel

open Real Finset MeasureTheory intervalIntegral

namespace Coupling

/-- `Θ c t = ∑_{p ≤ t, p prime} c(p)`. -/
noncomputable def Theta (c : ℕ → ℝ) (t : ℝ) : ℝ :=
  ∑ p ∈ Icc 0 ⌊t⌋₊ with p.Prime, c p

/-- `Θ₁ c t = ∑_{p ≤ t, p prime} c(p) (t − p)`. -/
noncomputable def Theta1 (c : ℕ → ℝ) (t : ℝ) : ℝ :=
  ∑ p ∈ Icc 0 ⌊t⌋₊ with p.Prime, c p * (t - p)

/-- For `0 ≤ t ≤ x`, `Θ₁ c t = ∑_{p ≤ x} c(p) max(t − p, 0)`. -/
lemma Theta1_eq_sum_max (c : ℕ → ℝ) {t x : ℝ} (ht : 0 ≤ t) (htx : t ≤ x) :
    Theta1 c t = ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, c p * max (t - p) 0 := by
  unfold Theta1
  have hsub : (Icc 0 ⌊t⌋₊).filter Nat.Prime ⊆ (Icc 0 ⌊x⌋₊).filter Nat.Prime :=
    filter_subset_filter _ (Icc_subset_Icc le_rfl (Nat.floor_le_floor htx))
  rw [← sum_subset hsub]
  · refine sum_congr rfl fun p hp => ?_
    have hp' : (p : ℝ) ≤ t := by
      have := (mem_Icc.1 (mem_filter.1 hp).1).2
      exact (Nat.le_floor_iff ht).1 this
    rw [max_eq_left (by linarith)]
  · intro p hp hpn
    have hpt : t < p := by
      by_contra h
      push Not at h
      exact hpn (mem_filter.2 ⟨mem_Icc.2 ⟨Nat.zero_le _, Nat.le_floor h⟩, (mem_filter.1 hp).2⟩)
    rw [max_eq_right (by linarith), mul_zero]

/-- Taylor's formula with integral remainder, in the form used for each prime:
`F p = F x + F'(x)(p − x) + ∫_p^x (t − p) F''(t) dt`. -/
lemma taylor_integral {F F' F'' : ℝ → ℝ} {a x p : ℝ} (hap : a ≤ p) (hpx : p ≤ x)
    (hF : ∀ t ∈ Set.Icc a x, HasDerivAt F (F' t) t)
    (hF' : ∀ t ∈ Set.Icc a x, HasDerivAt F' (F'' t) t)
    (hF'' : ContinuousOn F'' (Set.Icc a x)) :
    F p = F x + F' x * (p - x) + ∫ t in p..x, (t - p) * F'' t := by
  have hsub : Set.uIcc p x ⊆ Set.Icc a x := by
    rw [Set.uIcc_of_le hpx]; exact Set.Icc_subset_Icc hap le_rfl
  have hg : ∀ t ∈ Set.uIcc p x,
      HasDerivAt (fun t => F t + F' t * (p - t)) ((p - t) * F'' t) t := by
    intro t ht
    have h1 := hF t (hsub ht)
    have h2 := (hF' t (hsub ht)).mul ((hasDerivAt_const t p).sub (hasDerivAt_id t))
    refine (h1.add h2).congr_deriv ?_
    simp only [Pi.sub_apply, id]
    ring
  have hcont : ContinuousOn (fun t => (p - t) * F'' t) (Set.uIcc p x) :=
    (continuousOn_const.sub continuousOn_id).mul (hF''.mono hsub)
  have hFTC := integral_eq_sub_of_hasDerivAt hg (hcont.intervalIntegrable)
  have hneg : ∫ t in p..x, (t - p) * F'' t = -∫ t in p..x, (p - t) * F'' t := by
    rw [← intervalIntegral.integral_neg]; congr 1; ext t; ring
  rw [hneg, hFTC]; ring

/-- `∫_p^x (t − p) F''(t) dt = ∫_a^x max(t − p, 0) F''(t) dt` for `a ≤ p ≤ x`. -/
lemma integral_max_eq {F'' : ℝ → ℝ} {a x p : ℝ} (hap : a ≤ p) (hpx : p ≤ x)
    (hF'' : ContinuousOn F'' (Set.Icc a x)) :
    ∫ t in p..x, (t - p) * F'' t = ∫ t in a..x, max (t - p) 0 * F'' t := by
  have hcont : ContinuousOn (fun t => max (t - p) 0 * F'' t) (Set.Icc a x) :=
    ((continuous_id.sub continuous_const).max continuous_const).continuousOn.mul hF''
  have hi1 : IntervalIntegrable (fun t => max (t - p) 0 * F'' t) volume a p :=
    (hcont.mono (by rw [Set.uIcc_of_le hap]; exact Set.Icc_subset_Icc le_rfl hpx)).intervalIntegrable
  have hi2 : IntervalIntegrable (fun t => max (t - p) 0 * F'' t) volume p x :=
    (hcont.mono (by rw [Set.uIcc_of_le hpx]; exact Set.Icc_subset_Icc hap le_rfl)).intervalIntegrable
  rw [← integral_add_adjacent_intervals hi1 hi2]
  have hz : ∫ t in a..p, max (t - p) 0 * F'' t = 0 := by
    rw [integral_congr (g := fun _ => (0 : ℝ))]
    · simp
    · intro t ht
      rw [Set.uIcc_of_le hap] at ht
      simp [max_eq_right (by linarith [ht.2] : t - p ≤ 0)]
  rw [hz, zero_add]
  refine integral_congr fun t ht => ?_
  rw [Set.uIcc_of_le hpx] at ht
  simp [max_eq_left (by linarith [ht.1] : 0 ≤ t - p)]

/-- **Lemma 2 of `EXPLICIT.md` (second-order partial summation).** For any real weight `c`
on the primes, `x ≥ 2`, and `F` twice differentiable on `[2, x]` with continuous `F''`,
`∑_{p ≤ x} c(p) F(p) = F(x) Θ(x) − F'(x) Θ₁(x) + ∫_2^x Θ₁(t) F''(t) dt`. -/
theorem second_order_partial_summation (c : ℕ → ℝ) {F F' F'' : ℝ → ℝ} {x : ℝ} (hx : 2 ≤ x)
    (hF : ∀ t ∈ Set.Icc 2 x, HasDerivAt F (F' t) t)
    (hF' : ∀ t ∈ Set.Icc 2 x, HasDerivAt F' (F'' t) t)
    (hF'' : ContinuousOn F'' (Set.Icc 2 x)) :
    ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, c p * F p =
      F x * Theta c x - F' x * Theta1 c x + ∫ t in (2 : ℝ)..x, Theta1 c t * F'' t := by
  set S := (Icc 0 ⌊x⌋₊).filter Nat.Prime with hS
  have hmem : ∀ p ∈ S, (2 : ℝ) ≤ p ∧ (p : ℝ) ≤ x := by
    intro p hp
    obtain ⟨hp1, hp2⟩ := mem_filter.1 hp
    exact ⟨by exact_mod_cast hp2.two_le,
      (Nat.le_floor_iff (by linarith)).1 (mem_Icc.1 hp1).2⟩
  -- Taylor at each prime
  have hT : ∀ p ∈ S, c p * F p = c p * F x - F' x * (c p * (x - p)) +
      ∫ t in (2 : ℝ)..x, c p * (max (t - p) 0 * F'' t) := by
    intro p hp
    obtain ⟨h2, hx'⟩ := hmem p hp
    rw [taylor_integral h2 hx' hF hF' hF'', integral_max_eq h2 hx' hF'',
      intervalIntegral.integral_const_mul]
    ring
  rw [sum_congr rfl hT, sum_add_distrib, sum_sub_distrib, ← sum_mul, ← mul_sum]
  have hint : ∀ p ∈ S, IntervalIntegrable (fun t => c p * (max (t - p) 0 * F'' t)) volume 2 x :=
    fun p _ => (continuousOn_const.mul
      (((continuous_id.sub continuous_const).max continuous_const).continuousOn.mul
        (hF''.mono (by rw [Set.uIcc_of_le hx])))).intervalIntegrable
  rw [← intervalIntegral.integral_finsetSum hint]
  have hΘ1 : Theta1 c x = ∑ p ∈ S, c p * (x - p) := rfl
  have hΘ : Theta c x = ∑ p ∈ S, c p := rfl
  rw [hΘ1, hΘ, mul_comm (F x)]
  congr 1
  refine integral_congr fun t ht => ?_
  rw [Set.uIcc_of_le hx] at ht
  rw [Theta1_eq_sum_max c (by linarith [ht.1]) ht.2, sum_mul]
  refine sum_congr rfl fun p _ => ?_
  ring

/-! ### The kernel `h(t) = (x − t)/log t` -/

/-- `h'(t)` for `h(t) = (x − t)/log t`. -/
noncomputable def hD1 (x t : ℝ) : ℝ := -(1 / log t) - (x - t) / (t * log t ^ 2)

/-- `h''(t)` for `h(t) = (x − t)/log t`. -/
noncomputable def hD2 (x t : ℝ) : ℝ :=
  2 / (t * log t ^ 2) + (x - t) / (t ^ 2 * log t ^ 2) + 2 * (x - t) / (t ^ 2 * log t ^ 3)

lemma hasDerivAt_h {x t : ℝ} (ht : 1 < t) :
    HasDerivAt (fun t => (x - t) / log t) (hD1 x t) t := by
  have ht0 : t ≠ 0 := by positivity
  have hl : log t ≠ 0 := (Real.log_pos ht).ne'
  have := ((hasDerivAt_const t x).sub (hasDerivAt_id t)).div (Real.hasDerivAt_log ht0) hl
  refine this.congr_deriv ?_
  unfold hD1
  simp only [Pi.sub_apply, id]
  field_simp
  ring

lemma hasDerivAt_hD1 {x t : ℝ} (ht : 1 < t) :
    HasDerivAt (fun t => hD1 x t) (hD2 x t) t := by
  have ht0 : t ≠ 0 := by positivity
  have hl : log t ≠ 0 := (Real.log_pos ht).ne'
  have hlog := Real.hasDerivAt_log ht0
  have h1 := (hasDerivAt_const t (1 : ℝ)).div hlog hl
  have hden := (hasDerivAt_id t).mul (hlog.pow 2)
  have h2 := ((hasDerivAt_const t x).sub (hasDerivAt_id t)).div hden
    (by simp only [Pi.mul_apply, Pi.pow_apply, id]; exact mul_ne_zero ht0 (pow_ne_zero _ hl))
  have := h1.neg.sub h2
  show HasDerivAt (fun t => -(1 / log t) - (x - t) / (t * log t ^ 2)) _ t
  refine this.congr_deriv ?_
  unfold hD2
  simp only [Pi.mul_apply, Pi.pow_apply, Pi.sub_apply, id, Nat.cast_ofNat]
  field_simp
  ring

lemma hD2_pos {x t : ℝ} (ht : 1 < t) (htx : t ≤ x) : 0 < hD2 x t := by
  have hl : 0 < log t := Real.log_pos ht
  have ht0 : 0 < t := by linarith
  unfold hD2
  have : 0 ≤ x - t := by linarith
  positivity

lemma continuousOn_hD2 (x : ℝ) : ContinuousOn (hD2 x) (Set.Icc 2 x) := by
  intro t ht
  have ht0 : t ≠ 0 := by linarith [ht.1]
  have hl : log t ≠ 0 := (Real.log_pos (by linarith [ht.1])).ne'
  apply ContinuousAt.continuousWithinAt
  unfold hD2
  fun_prop (disch := simp [*])

/-- **Identity (2.1) of `EXPLICIT.md`.** For any real weight `w` and `x ≥ 2`, with
`Θ₁(t) = ∑_{p ≤ t} w(p) log p (t − p)`,
`x R(x) = Θ₁(x)/log x + ∫_2^x Θ₁(t) h''(t) dt`, where `h(t) = (x − t)/log t`. -/
theorem cesaro_convexity_identity (w : ℕ → ℝ) {x : ℝ} (hx : 2 ≤ x) :
    x * cesaroSum w x =
      Theta1 (fun p => w p * log p) x / log x +
        ∫ t in (2 : ℝ)..x, Theta1 (fun p => w p * log p) t * hD2 x t := by
  have hF : ∀ t ∈ Set.Icc 2 x, HasDerivAt (fun t => (x - t) / log t) (hD1 x t) t :=
    fun t ht => hasDerivAt_h (by linarith [ht.1])
  have hF' : ∀ t ∈ Set.Icc 2 x, HasDerivAt (fun t => hD1 x t) (hD2 x t) t :=
    fun t ht => hasDerivAt_hD1 (by linarith [ht.1])
  have key := second_order_partial_summation (fun p => w p * log p) hx hF hF'
    (continuousOn_hD2 x)
  have hlx : log x ≠ 0 := (Real.log_pos (by linarith)).ne'
  have hx0 : x ≠ 0 := by linarith
  have hlhs : ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, w p * log p * ((x - p) / log p) =
      x * cesaroSum w x := by
    rw [cesaroSum, mul_sum]
    refine sum_congr rfl fun p hp => ?_
    have hp2 : (2 : ℝ) ≤ p := by exact_mod_cast (mem_filter.1 hp).2.two_le
    have hlp : log p ≠ 0 := (Real.log_pos (by linarith)).ne'
    field_simp
  rw [hlhs] at key
  rw [key]
  simp only [sub_self, zero_div, zero_mul, zero_sub, hD1]
  field_simp
  ring

/-- **Target B (positivity transfer).** If `x ≥ 2`, `Θ₁(t) ≥ 0` for all `t ∈ [2, x]` and
`Θ₁(x) > 0`, then `R(x) ≥ Θ₁(x)/(x log x) > 0`. -/
theorem cesaro_pos_of_theta1 (w : ℕ → ℝ) {x : ℝ} (hx : 2 ≤ x)
    (hnn : ∀ t ∈ Set.Icc 2 x, 0 ≤ Theta1 (fun p => w p * log p) t)
    (hpos : 0 < Theta1 (fun p => w p * log p) x) :
    Theta1 (fun p => w p * log p) x / (x * log x) ≤ cesaroSum w x ∧ 0 < cesaroSum w x := by
  have hlx : 0 < log x := Real.log_pos (by linarith)
  have hx0 : 0 < x := by linarith
  have hI : 0 ≤ ∫ t in (2 : ℝ)..x, Theta1 (fun p => w p * log p) t * hD2 x t :=
    intervalIntegral.integral_nonneg hx fun t ht =>
      mul_nonneg (hnn t ht) (hD2_pos (by linarith [ht.1]) ht.2).le
  have key := cesaro_convexity_identity w hx
  have hle : Theta1 (fun p => w p * log p) x / (x * log x) ≤ cesaroSum w x := by
    rw [div_le_iff₀ (by positivity)]
    have : Theta1 (fun p => w p * log p) x / log x ≤ x * cesaroSum w x := by linarith
    rw [div_le_iff₀ hlx] at this
    linarith
  exact ⟨hle, lt_of_lt_of_le (by positivity) hle⟩

end Coupling
