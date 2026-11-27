/-
AllxLib/Final.lean — assembly: under GRH for `L(s, χ₋₄)` the Cesàro mean of Chebyshev's race
`R(x) = Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x)` is positive for EVERY real `x > 3`.

* finite range `3 < t ≤ 10⁷`: `theta1_wlog_pos_of_le` (kernel-checked, no hypothesis);
* analytic range `t ≥ 10⁷`: `theta1_wlog_pos_of_ge` (GRH for `L(s, χ₋₄)`);
* `Coupling.cesaro_pos_of_theta1`: `x R(x) = Θ₁(x)/log x + ∫₂ˣ Θ₁(t) h''(t) dt` with `h'' > 0`.
-/
import AllxLib.Analytic.Main
import AllxLib.Finite.Final
import AllxLib.Coupling.Cone
import AllxLib.Defs

open Real Finset MeasureTheory

noncomputable section

namespace Allx

/-- `Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p) > 0` for every `t > 3`, under GRH for `L(s,χ₋₄)`. -/
theorem theta1_wlog_pos (hGRH : GRHFor chi4) {t : ℝ} (h3 : 3 < t) :
    0 < Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t := by
  rcases le_total t 10000000 with h | h
  · exact theta1_wlog_pos_of_le h3 h
  · exact theta1_wlog_pos_of_ge hGRH h

theorem theta1_wlog_nonneg (hGRH : GRHFor chi4) {t : ℝ} (h2 : 2 ≤ t) :
    0 ≤ Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t := by
  rcases le_total t 10000000 with h | h
  · exact theta1_wlog_nonneg_of_le h2 h
  · exact (theta1_wlog_pos_of_ge hGRH h).le

/-- **The Cesàro mean of Chebyshev's race is positive for every `x > 3`** (GRH for `L(s,χ₋₄)`),
with the lower bound `Θ₁(x)/(x log x)`. -/
theorem cesaro_w4_ge_and_pos (hGRH : GRHFor chi4) {x : ℝ} (hx : 3 < x) :
    Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) x / (x * Real.log x)
        ≤ Coupling.cesaroSum Coupling.w4 x ∧ 0 < Coupling.cesaroSum Coupling.w4 x :=
  Coupling.cesaro_pos_of_theta1 Coupling.w4 (by linarith)
    (fun _ ht => theta1_wlog_nonneg hGRH ht.1) (theta1_wlog_pos hGRH hx)

/-- The finite range, with no hypothesis: `R(x) > 0` for `3 < x ≤ 10⁷`. -/
theorem cesaro_w4_pos_of_le {x : ℝ} (hx : 3 < x) (hN : x ≤ 10000000) :
    0 < Coupling.cesaroSum Coupling.w4 x :=
  (Coupling.cesaro_pos_of_theta1 Coupling.w4 (by linarith)
    (fun _ ht => theta1_wlog_nonneg_of_le ht.1 (ht.2.trans hN)) (theta1_wlog_pos_of_le hx hN)).2

/-! ### The Cesàro mean as an integral of `π(t;4,3) − π(t;4,1)` -/

theorem primeSum_w4_zero (t : ℝ) :
    Coupling.primeSum Coupling.w4 0 t
      = (primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ) := by
  unfold Coupling.primeSum primeCountingMod4
  have e : ∀ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime, Coupling.w4 p * (p : ℝ) ^ (0 : ℝ)
      = (if p % 4 = 3 then (1 : ℝ) else 0) - (if p % 4 = 1 then (1 : ℝ) else 0) := by
    intro p _
    rw [Real.rpow_zero, mul_one, w4_eq]
  rw [sum_congr rfl e, sum_sub_distrib, sum_boole, sum_boole, filter_filter, filter_filter]

/-- `R(x) = (1/x) ∫₀ˣ (π(t;4,3) − π(t;4,1)) dt` for `x > 0` (no hypothesis). -/
theorem cesaro_w4_eq_integral {x : ℝ} (hx : 0 < x) :
    Coupling.cesaroSum Coupling.w4 x
      = (1 / x) * ∫ t in (0 : ℝ)..x, ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ)) := by
  rw [Coupling.cesaro_eq_integral Coupling.w4 hx]
  congr 1
  refine intervalIntegral.integral_congr fun t _ => ?_
  exact primeSum_w4_zero t

/-- `∫₀ˣ (π(t;4,3) − π(t;4,1)) dt > 0` for every `x > 3`, under GRH for `L(s,χ₋₄)`. -/
theorem integral_race_pos (hGRH : GRHFor chi4) {x : ℝ} (hx : 3 < x) :
    0 < ∫ t in (0 : ℝ)..x, ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ)) := by
  have hx0 : 0 < x := by linarith
  have h := (cesaro_w4_ge_and_pos hGRH hx).2
  rw [cesaro_w4_eq_integral hx0] at h
  have h1 : 0 < (1 / x) := by positivity
  exact (mul_pos_iff_of_pos_left h1).mp h

/-! ### Inheritance of leads -/

/-- `D₀(x) = D₁(x)/x + R(x)` where `D_m(x) = Σ_{p ≤ x} (−χ₋₄(p)) p^m`. -/
theorem sum_w4_eq (x : ℝ) (hx : 0 < x) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, Coupling.w4 p
      = (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, Coupling.w4 p * (p : ℝ)) / x
        + Coupling.cesaroSum Coupling.w4 x := by
  unfold Coupling.cesaroSum
  rw [sum_div, ← sum_add_distrib]
  refine sum_congr rfl fun p _ => ?_
  field_simp
  ring

/-- **Inheritance of leads** (`m = 0` against `m' = 1`): for every `x > 3`, under GRH for
`L(s,χ₋₄)`, if `Σ_{p ≤ x} (−χ₋₄(p)) p ≥ 0` then `Σ_{p ≤ x} (−χ₋₄(p)) > 0`. -/
theorem lead_inheritance_w4 (hGRH : GRHFor chi4) {x : ℝ} (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, Coupling.w4 p * (p : ℝ)) :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, Coupling.w4 p := by
  have hx0 : 0 < x := by linarith
  rw [sum_w4_eq x hx0]
  have h1 := (cesaro_w4_ge_and_pos hGRH hx).2
  have h2 : 0 ≤ (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, Coupling.w4 p * (p : ℝ)) / x :=
    div_nonneg h hx0.le
  linarith

/-- `Σ_{p ≤ x} (−χ₋₄(p)) = π(x;4,3) − π(x;4,1)`. -/
theorem sum_w4_eq_counts (x : ℝ) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, Coupling.w4 p
      = (primeCountingMod4 3 x : ℝ) - (primeCountingMod4 1 x : ℝ) := by
  rw [← primeSum_w4_zero]
  unfold Coupling.primeSum
  refine sum_congr rfl fun p _ => ?_
  rw [Real.rpow_zero, mul_one]

/-! ### The cone of weights `m ≤ 0`, `m < m' ≤ 1` (`EXPLICIT.md`, Theorem 5) -/

/-- For `m ≤ 0`, `m < m' ≤ 1` and every `x > 3`, under GRH for `L(s,χ₋₄)`:
`D_m(x) − x^{−(m'−m)} D_{m'}(x) ≥ (m'−m) x^{m−1} Θ₁(x)/log x`. -/
theorem cone_w4_lower_bound (hGRH : GRHFor chi4) {m m' x : ℝ} (hm : m ≤ 0) (hmm' : m < m')
    (hm'1 : m' ≤ 1) (hx : 3 < x) :
    (m' - m) * x ^ (m - 1) * Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) x / Real.log x
      ≤ Coupling.primeSum Coupling.w4 m x - x ^ (-(m' - m)) * Coupling.primeSum Coupling.w4 m' x :=
  Coupling.cone_lower_bound Coupling.w4 hm hmm' hm'1 (by linarith)
    (fun _ ht => theta1_wlog_nonneg hGRH ht.1)

/-- Inheritance of leads in the cone: `D_{m'}(x) ≥ 0 ⟹ D_m(x) > 0` for every `x > 3`. -/
theorem cone_w4_lead_inheritance (hGRH : GRHFor chi4) {m m' x : ℝ} (hm : m ≤ 0) (hmm' : m < m')
    (hm'1 : m' ≤ 1) (hx : 3 < x) (h : 0 ≤ Coupling.primeSum Coupling.w4 m' x) :
    0 < Coupling.primeSum Coupling.w4 m x :=
  Coupling.cone_lead_inheritance Coupling.w4 hm hmm' hm'1 (by linarith)
    (fun _ ht => theta1_wlog_nonneg hGRH ht.1) (theta1_wlog_pos hGRH hx) h

end Allx
