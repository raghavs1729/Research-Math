/-
AllxLib/Analytic/ZeroSum.lean — the two positive zero sums, in closed form, under GRH for `L(s,χ₋₄)`.

* `sumMaj_eq`:   `Σ_ρ m_ρ · 30/((9/4+γ²)(25/4+γ²)) = 7/2 − γ − 3 log 2 − log π + 5A − 3B`
  (this is `ChebLib.EF_maj`, rewritten with `A = Aser`, `B = Bser`);
* `sumOne_le`:   `Σ_ρ m_ρ · 3/(9/4+γ²) ≤ 2 − γ − log π + 2A`.

The second kernel `3/(9/4+r²)` is the transform of `e^{−(3/2)|u|}`, which has a corner; the
explicit formula `EF_ext` needs `C²`.  We apply it to
`gN c u = (c e^{−(3/2)|u|} − (3/2) e^{−c|u|})/(c − 3/2)`, whose transform
`3c(c+3/2)/((9/4+r²)(c²+r²))` tends to `3/(9/4+r²)` as `c → ∞`, bound the right-hand side by
`c/(c−3/2) · (2 − γ − log π + 2A) + (3/(c−3/2)) Σ Λ(n)/n²` (archimedean side: `arch_weil`), and let
`c → ∞` (Fatou for sums of nonnegative terms).
-/
import ChebLib.Cheb.ZeroSum
import AllxLib.Analytic.ArchWeil
import AllxLib.Analytic.KernelK
import AllxLib.Analytic.Numerics

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE ChebLib

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace Allx

/-! ### The `C²` approximation of `e^{−(3/2)|u|}` -/

/-- `gN c u = (c e^{−(3/2)|u|} − (3/2) e^{−c|u|})/(c − 3/2)`. -/
def gN (c u : ℝ) : ℝ :=
  (c * Real.exp (-(3 / 2) * |u|) - 3 / 2 * Real.exp (-c * |u|)) / (c - 3 / 2)

/-- First derivative of `gN c`. -/
def gN1 (c u : ℝ) : ℝ :=
  if 0 ≤ u then 3 / 2 * c / (c - 3 / 2) * (-Real.exp (-(3 / 2) * u) + Real.exp (-c * u))
  else 3 / 2 * c / (c - 3 / 2) * (Real.exp (3 / 2 * u) - Real.exp (c * u))

/-- Second derivative of `gN c`. -/
def gN2 (c u : ℝ) : ℝ :=
  3 / 2 * c / (c - 3 / 2) * (3 / 2 * Real.exp (-(3 / 2) * |u|) - c * Real.exp (-c * |u|))

theorem gN_neg (c u : ℝ) : gN c (-u) = gN c u := by simp [gN]

theorem gN_zero {c : ℝ} (hc : 3 / 2 < c) : gN c 0 = 1 := by
  have : c - 3 / 2 ≠ 0 := by linarith
  unfold gN
  simp only [abs_zero, mul_zero, Real.exp_zero, mul_one]
  exact div_self this

theorem hasDerivAt_gN (c u : ℝ) : HasDerivAt (gN c) (gN1 c u) u := by
  refine hasDerivAt_glue
    (A := fun u => (c * Real.exp (-(3 / 2) * u) - 3 / 2 * Real.exp (-c * u)) / (c - 3 / 2))
    (B := fun u => (c * Real.exp (3 / 2 * u) - 3 / 2 * Real.exp (c * u)) / (c - 3 / 2))
    (A' := fun u => 3 / 2 * c / (c - 3 / 2) * (-Real.exp (-(3 / 2) * u) + Real.exp (-c * u)))
    (B' := fun u => 3 / 2 * c / (c - 3 / 2) * (Real.exp (3 / 2 * u) - Real.exp (c * u)))
    (fun t => ?_) (fun t => ?_) (fun u hu => ?_) (fun u hu => ?_) (by simp) u
  · have h := (((hasDerivAt_exp_mul (-(3 / 2)) t).const_mul c).sub
      ((hasDerivAt_exp_mul (-c) t).const_mul (3 / 2))).div_const (c - 3 / 2)
    refine h.congr_deriv ?_
    ring
  · have h := (((hasDerivAt_exp_mul (3 / 2) t).const_mul c).sub
      ((hasDerivAt_exp_mul c t).const_mul (3 / 2))).div_const (c - 3 / 2)
    refine h.congr_deriv ?_
    ring
  · simp only [gN, abs_of_nonneg hu]
  · simp only [gN, abs_of_nonpos hu]
    rw [show -(3 / 2 : ℝ) * -u = 3 / 2 * u by ring, show -c * -u = c * u by ring]

theorem gN2_eq (c u : ℝ) : gN2 c u =
    if 0 ≤ u then 3 / 2 * c / (c - 3 / 2)
        * (3 / 2 * Real.exp (-(3 / 2) * u) - c * Real.exp (-c * u))
    else 3 / 2 * c / (c - 3 / 2) * (3 / 2 * Real.exp (3 / 2 * u) - c * Real.exp (c * u)) := by
  unfold gN2
  split_ifs with hu
  · rw [abs_of_nonneg hu]
  · rw [abs_of_nonpos (le_of_not_ge hu),
      show -(3 / 2 : ℝ) * -u = 3 / 2 * u by ring, show -c * -u = c * u by ring]

theorem hasDerivAt_gN1 (c u : ℝ) : HasDerivAt (gN1 c) (gN2 c u) u := by
  rw [gN2_eq]
  refine hasDerivAt_glue
    (A := fun u => 3 / 2 * c / (c - 3 / 2) * (-Real.exp (-(3 / 2) * u) + Real.exp (-c * u)))
    (B := fun u => 3 / 2 * c / (c - 3 / 2) * (Real.exp (3 / 2 * u) - Real.exp (c * u)))
    (A' := fun u => 3 / 2 * c / (c - 3 / 2)
      * (3 / 2 * Real.exp (-(3 / 2) * u) - c * Real.exp (-c * u)))
    (B' := fun u => 3 / 2 * c / (c - 3 / 2) * (3 / 2 * Real.exp (3 / 2 * u) - c * Real.exp (c * u)))
    (fun t => ?_) (fun t => ?_) (fun u hu => ?_) (fun u hu => ?_) (by simp) u
  · have h := (((hasDerivAt_exp_mul (-(3 / 2)) t).neg).add
      (hasDerivAt_exp_mul (-c) t)).const_mul (3 / 2 * c / (c - 3 / 2))
    refine h.congr_deriv ?_
    ring
  · have h := ((hasDerivAt_exp_mul (3 / 2) t).sub
      (hasDerivAt_exp_mul c t)).const_mul (3 / 2 * c / (c - 3 / 2))
    exact h
  · simp only [gN1, if_pos hu]
  · rcases eq_or_lt_of_le hu with h0 | hlt
    · subst h0; simp [gN1]
    · simp only [gN1, if_neg (not_le.mpr hlt)]

theorem continuous_gN (c : ℝ) : Continuous (gN c) := by unfold gN; fun_prop

theorem continuous_gN2 (c : ℝ) : Continuous (gN2 c) := by unfold gN2; fun_prop

open ArchWeilAux in
theorem integrable_gN {c : ℝ} (hc : 3 / 2 < c) : Integrable (gN c) :=
  (((integrable_exp_neg_mul_abs (b := 3 / 2) (by norm_num)).const_mul c).sub
    ((integrable_exp_neg_mul_abs (b := c) (by linarith)).const_mul (3 / 2))).div_const _

open ArchWeilAux in
theorem integrable_gN2 {c : ℝ} (hc : 3 / 2 < c) : Integrable (gN2 c) :=
  (((integrable_exp_neg_mul_abs (b := 3 / 2) (by norm_num)).const_mul (3 / 2)).sub
    ((integrable_exp_neg_mul_abs (b := c) (by linarith)).const_mul c)).const_mul _

open ArchWeilAux in
theorem integrable_gN1 {c : ℝ} (hc : 3 / 2 < c) : Integrable (gN1 c) := by
  have hcont : Continuous (gN1 c) :=
    continuous_iff_continuousAt.mpr fun u => (hasDerivAt_gN1 c u).continuousAt
  have hK : 0 ≤ 3 / 2 * c / (c - 3 / 2) := by
    apply div_nonneg <;> linarith
  refine Integrable.mono'
    (((integrable_exp_neg_mul_abs (b := 3 / 2) (by norm_num)).add
      (integrable_exp_neg_mul_abs (b := c) (by linarith))).const_mul (3 / 2 * c / (c - 3 / 2)))
    hcont.aestronglyMeasurable (Eventually.of_forall fun u => ?_)
  rw [Real.norm_eq_abs]
  simp only [Pi.add_apply]
  unfold gN1
  split_ifs with hu
  · rw [abs_mul, abs_of_nonneg hK, abs_of_nonneg hu]
    refine mul_le_mul_of_nonneg_left ?_ hK
    have h1 := Real.exp_pos (-(3 / 2) * u)
    have h2 := Real.exp_pos (-c * u)
    rw [abs_le]; constructor <;> linarith
  · rw [abs_mul, abs_of_nonneg hK, abs_of_nonpos (le_of_not_ge hu),
      show -(3 / 2 : ℝ) * -u = 3 / 2 * u by ring, show -c * -u = c * u by ring]
    refine mul_le_mul_of_nonneg_left ?_ hK
    have h1 := Real.exp_pos (3 / 2 * u)
    have h2 := Real.exp_pos (c * u)
    rw [abs_le]; constructor <;> linarith

/-- The transform of `gN c`. -/
def gNhat (c r : ℝ) : ℝ := 3 * c * (3 / 2 + c) / ((9 / 4 + r ^ 2) * (c ^ 2 + r ^ 2))

open ArchWeilAux in
theorem integral_gN_mul_cexp {c : ℝ} (hc : 3 / 2 < c) (r : ℝ) :
    ∫ u : ℝ, ((gN c u : ℝ) : ℂ) * cexp (I * r * u) = ((gNhat c r : ℝ) : ℂ) := by
  have hc0 : 0 < c := by linarith
  have e : ∀ u : ℝ, ((gN c u : ℝ) : ℂ) * cexp (I * r * u)
      = ((c / (c - 3 / 2) : ℝ) : ℂ) * ((Real.exp (-(3 / 2) * |u|) : ℂ) * cexp (I * r * u))
        - ((3 / 2 / (c - 3 / 2) : ℝ) : ℂ) * ((Real.exp (-c * |u|) : ℂ) * cexp (I * r * u)) := by
    intro u
    simp only [gN]
    push_cast
    ring
  simp_rw [e]
  rw [integral_sub ((integrable_kernel (b := 3 / 2) (by norm_num) r).const_mul _)
    ((integrable_kernel (b := c) hc0 r).const_mul _),
    integral_const_mul, integral_const_mul, integral_kernel (by norm_num), integral_kernel hc0]
  have hne : c - 3 / 2 ≠ 0 := by linarith
  have h5 : (9 / 4 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have h6 : (c ^ 2 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have hR : c / (c - 3 / 2) * (2 * (3 / 2) / ((3 / 2) ^ 2 + r ^ 2))
      - 3 / 2 / (c - 3 / 2) * (2 * c / (c ^ 2 + r ^ 2)) = gNhat c r := by
    unfold gNhat
    rw [show ((3 : ℝ) / 2) ^ 2 = 9 / 4 by norm_num, div_mul_div_comm, div_mul_div_comm,
      div_sub_div _ _ (mul_ne_zero hne h5) (mul_ne_zero hne h6),
      div_eq_div_iff (mul_ne_zero (mul_ne_zero hne h5) (mul_ne_zero hne h6)) (mul_ne_zero h5 h6)]
    ring
  rw [← hR]
  push_cast
  ring

theorem gNhat_nonneg {c : ℝ} (hc : 3 / 2 < c) (r : ℝ) : 0 ≤ gNhat c r := by
  unfold gNhat
  have : 0 < c := by linarith
  positivity

open ArchWeilAux in
/-- `∫ e^{−a|u|} du = 2/a`. -/
theorem integral_exp_neg_mul_abs {a : ℝ} (ha : 0 < a) :
    ∫ u : ℝ, Real.exp (-a * |u|) = 2 / a := by
  have h := integral_kernel ha 0
  simp only [Complex.ofReal_zero, mul_zero, zero_mul, Complex.exp_zero, mul_one] at h
  have h2 : ∫ u : ℝ, ((Real.exp (-a * |u|) : ℝ) : ℂ) = ((∫ u : ℝ, Real.exp (-a * |u|) : ℝ) : ℂ) :=
    integral_ofReal
  rw [h2] at h
  have h3 := Complex.ofReal_injective h
  rw [h3]
  field_simp
  ring

open ArchWeilAux in
theorem integral_gN_mul_exp {c : ℝ} (hc : 3 / 2 < c) {β : ℝ} (hβ : 0 < β) :
    ∫ u : ℝ, gN c u * Real.exp (-β * |u|)
      = c / (c - 3 / 2) * (2 / (3 / 2 + β)) - 3 / 2 / (c - 3 / 2) * (2 / (c + β)) := by
  have e : ∀ u : ℝ, gN c u * Real.exp (-β * |u|)
      = c / (c - 3 / 2) * Real.exp (-(3 / 2 + β) * |u|)
        - 3 / 2 / (c - 3 / 2) * Real.exp (-(c + β) * |u|) := by
    intro u
    unfold gN
    rw [show -(3 / 2 + β) * |u| = -(3 / 2) * |u| + -β * |u| by ring,
      show -(c + β) * |u| = -c * |u| + -β * |u| by ring, Real.exp_add, Real.exp_add]
    ring
  simp_rw [e]
  rw [integral_sub ((integrable_exp_neg_mul_abs (b := 3 / 2 + β) (by linarith)).const_mul _)
    ((integrable_exp_neg_mul_abs (b := c + β) (by linarith)).const_mul _),
    integral_const_mul, integral_const_mul, integral_exp_neg_mul_abs (by linarith),
    integral_exp_neg_mul_abs (by linarith)]

/-! ### The archimedean side for `gN c` -/

open ArchWeilAux ChebLib.PoissonAux in
/-- The archimedean term of `gN c` is at most `c/(c−3/2) · (2 − γ − log π)` for `c ≥ 20`. -/
theorem arch_gN_le {c : ℝ} (hc : 20 ≤ c) :
    Integrable (fun r : ℝ => paperFT (fun u => ((gN c u : ℝ) : ℂ)) (r : ℂ) * ((W4 r : ℝ) : ℂ)) ∧
    ∃ a : ℝ, a ≤ c / (c - 3 / 2)
        * (2 - Real.eulerMascheroniConstant - Real.log Real.pi) ∧
      (1 / (2 * Real.pi) : ℂ)
        * ∫ r : ℝ, paperFT (fun u => ((gN c u : ℝ) : ℂ)) (r : ℂ) * ((W4 r : ℝ) : ℂ) = (a : ℂ) := by
  have hc' : 3 / 2 < c := by linarith
  have hc0 : 0 < c := by linarith
  have hcd : 0 < c - 3 / 2 := by linarith
  have hM : ∀ r : ℝ, ‖paperFT (fun u => ((gN c u : ℝ) : ℂ)) (r : ℂ)‖ * (1 + r ^ 2)
      ≤ 3 * (3 / 2 + c) / c := by
    intro r
    have : paperFT (fun u => ((gN c u : ℝ) : ℂ)) (r : ℂ) = ((gNhat c r : ℝ) : ℂ) :=
      integral_gN_mul_cexp hc' r
    rw [this, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (gNhat_nonneg hc' r)]
    unfold gNhat
    have x0 : 0 ≤ r ^ 2 := sq_nonneg r
    rw [div_mul_eq_mul_div, div_le_div_iff₀ (by positivity) hc0]
    have h1 : 1 + r ^ 2 ≤ 9 / 4 + r ^ 2 := by linarith
    have h2 : c * c ≤ c ^ 2 + r ^ 2 := by nlinarith
    have h3 : (0 : ℝ) ≤ 3 * (3 / 2 + c) := by linarith
    calc 3 * c * (3 / 2 + c) * (1 + r ^ 2) * c = 3 * (3 / 2 + c) * ((1 + r ^ 2) * (c * c)) := by ring
      _ ≤ 3 * (3 / 2 + c) * ((9 / 4 + r ^ 2) * (c ^ 2 + r ^ 2)) := by
          apply mul_le_mul_of_nonneg_left _ h3
          exact mul_le_mul h1 h2 (by positivity) (by positivity)
  obtain ⟨hsum, hint, heq⟩ := arch_weil (continuous_gN c) (integrable_gN hc') hM
  refine ⟨hint, _, ?_, heq⟩
  rw [gN_zero hc']
  -- the two ingredients
  set a : ℝ := c / (c - 3 / 2) with ha
  set d : ℝ := 3 / 2 / (c - 3 / 2) with hd
  have had : a - d = 1 := by rw [ha, hd, ← sub_div, div_self hcd.ne']
  have ha0 : 0 ≤ a := by rw [ha]; positivity
  have hd0 : 0 < d := by rw [hd]; positivity
  have hI : ∀ {β : ℝ}, 0 < β → ∫ u : ℝ, gN c u * Real.exp (-β * |u|)
      = a * (2 / (3 / 2 + β)) - d * (2 / (c + β)) := fun hβ => integral_gN_mul_exp hc' hβ
  have hterm : ∀ n : ℕ, (1 / ((n : ℝ) + 1)
      - ∫ u : ℝ, gN c u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|))
      = a * (1 / ((n : ℝ) + 1) - 1 / (3 / 2 + (n : ℝ) + 1))
        - d * (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2))) := by
    intro n
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    rw [hI (by linarith)]
    have e1 : (2 : ℝ) / (3 / 2 + (2 * (n : ℝ) + 7 / 2)) = 1 / (3 / 2 + (n : ℝ) + 1) := by
      rw [div_eq_div_iff (by linarith) (by linarith)]; ring
    rw [e1]
    linear_combination (-(1 / ((n : ℝ) + 1))) * had
  have hS1 := hasSum_three_halves
  have hsum' : Summable (fun n : ℕ => a * (1 / ((n : ℝ) + 1) - 1 / (3 / 2 + (n : ℝ) + 1))
        - d * (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2)))) := by
    refine hsum.congr fun n => ?_
    rw [gN_zero hc']
    exact hterm n
  -- the second series is summable and nonnegative
  have hS2 : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2))) := by
    have h := ((hS1.summable.mul_left a).sub hsum').mul_left (1 / d)
    refine h.congr fun n => ?_
    field_simp
    ring
  have hpos : ∀ n : ℕ, 0 ≤ 1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2)) := by
    intro n
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    rw [sub_nonneg, div_le_div_iff₀ (by linarith) (by linarith)]
    linarith
  have hS2ge : 1 - 2 / (c + 7 / 2) ≤ ∑' n : ℕ, (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2))) := by
    have := hS2.le_tsum 0 (fun n _ => hpos n)
    simpa using this
  have htsum : ∑' n : ℕ, (1 / ((n : ℝ) + 1)
      - ∫ u : ℝ, gN c u * Real.exp (-(2 * (n : ℝ) + 7 / 2) * |u|))
      = a * (8 / 3 - 2 * Real.log 2)
        - d * ∑' n : ℕ, (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2))) := by
    rw [tsum_congr hterm, (hS1.summable.mul_left a).tsum_sub (hS2.mul_left d), tsum_mul_left,
      tsum_mul_left, hS1.tsum_eq]
  rw [htsum, hI (by norm_num : (0 : ℝ) < 3 / 2)]
  have hlog : Real.log (4 / Real.pi) = 2 * Real.log 2 - Real.log Real.pi := by
    rw [Real.log_div (by norm_num) Real.pi_ne_zero, show (4 : ℝ) = 2 ^ 2 by norm_num,
      Real.log_pow]
    push_cast
    ring
  have hL0 : 0 ≤ Real.log (4 / Real.pi) :=
    Real.log_nonneg (by rw [le_div_iff₀ Real.pi_pos]; linarith [Real.pi_le_four])
  have hγ := Real.eulerMascheroniConstant_lt_two_thirds
  -- the `d`-bracket is nonnegative
  have h2a : 2 / (c + 3 / 2) ≤ 2 / 21 := by
    rw [div_le_div_iff₀ (by linarith) (by norm_num)]; linarith
  have h2b : 2 / (c + 7 / 2) ≤ 2 / 23 := by
    rw [div_le_div_iff₀ (by linarith) (by norm_num)]; linarith
  have hbr : 0 ≤ (Real.log (4 / Real.pi) - Real.eulerMascheroniConstant) - 2 / (c + 3 / 2)
      + ∑' n : ℕ, (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2))) := by
    linarith
  have key : 1 * (Real.log (4 / Real.pi) - Real.eulerMascheroniConstant)
      - (a * (2 / (3 / 2 + 3 / 2)) - d * (2 / (c + 3 / 2)))
      + (a * (8 / 3 - 2 * Real.log 2)
        - d * ∑' n : ℕ, (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2))))
      = a * (2 - Real.eulerMascheroniConstant - Real.log Real.pi)
        - d * ((Real.log (4 / Real.pi) - Real.eulerMascheroniConstant) - 2 / (c + 3 / 2)
          + ∑' n : ℕ, (1 / ((n : ℝ) + 1) - 2 / (c + (2 * (n : ℝ) + 7 / 2)))) := by
    have : (1 : ℝ) = a - d := had.symm
    rw [this, hlog]
    ring
  rw [key]
  have := mul_nonneg hd0.le hbr
  linarith

/-! ### The explicit formula for `gN c` and the limit -/

theorem exp_neg_three_halves_log {n : ℕ} (hn : 0 < n) :
    Real.exp (-(3 / 2) * Real.log n) = ((n : ℝ) * Real.sqrt n)⁻¹ := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  rw [← rpow_neg_three_halves hn', Real.rpow_def_of_pos hn']
  congr 1
  ring

/-- The real prime-side term of the explicit formula for `gN c`. -/
def pterm (c : ℝ) (n : ℕ) : ℝ := vM n * (-X4 n) * (2 * gN c (Real.log n) / Real.sqrt n)

theorem pterm_le {c : ℝ} (hc : 3 / 2 < c) (n : ℕ) :
    |pterm c n| ≤ (2 * c / (c - 3 / 2) + 3 / (c - 3 / 2)) * (vM n / (n : ℝ) ^ 2) ∧
    pterm c n ≤ 2 * c / (c - 3 / 2) * (vM n * (-X4 n) / (n : ℝ) ^ 2)
      + 3 / (c - 3 / 2) * (vM n / (n : ℝ) ^ 2) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [pterm, vM]
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
  have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
  have hlog : 0 ≤ Real.log n := Real.log_nonneg h1
  have hcd : 0 < c - 3 / 2 := by linarith
  have hΛ := vM_nonneg n
  have hX := abs_X4_le_one n
  have hX1 : -X4 n ≤ 1 := by linarith [neg_abs_le (X4 n)]
  have hX2 : -1 ≤ -X4 n := by linarith [le_abs_self (X4 n)]
  -- `e = e^{-c log n}/√n ∈ [0, 1/n²]`
  set e : ℝ := Real.exp (-c * Real.log n) / Real.sqrt n with he
  have he0 : 0 ≤ e := by rw [he]; positivity
  have he1 : e ≤ 1 / (n : ℝ) ^ 2 := by
    have : Real.exp (-c * Real.log n) ≤ Real.exp (-(3 / 2) * Real.log n) :=
      Real.exp_le_exp.mpr (by nlinarith)
    rw [exp_neg_three_halves_log hn] at this
    rw [he, div_le_iff₀ hs]
    calc Real.exp (-c * Real.log n) ≤ ((n : ℝ) * Real.sqrt n)⁻¹ := this
      _ = 1 / (n : ℝ) ^ 2 * Real.sqrt n := by
          generalize Real.sqrt n = s at hs hss
          rw [hss]; field_simp
  have hg : 2 * gN c (Real.log n) / Real.sqrt n
      = 2 * c / (c - 3 / 2) * (1 / (n : ℝ) ^ 2) - 3 / (c - 3 / 2) * e := by
    rw [he]
    unfold gN
    rw [abs_of_nonneg hlog, exp_neg_three_halves_log hn]
    generalize Real.exp (-c * Real.log n) = E
    generalize Real.sqrt n = s at hs hss
    rw [hss]
    ring
  have hK1 : 0 ≤ 2 * c / (c - 3 / 2) := by apply div_nonneg <;> linarith
  have hK2 : 0 ≤ 3 / (c - 3 / 2) := by positivity
  have hq : 0 ≤ 1 / (n : ℝ) ^ 2 := by positivity
  have e1 : pterm c n = 2 * c / (c - 3 / 2) * (vM n * (-X4 n) / (n : ℝ) ^ 2)
      - 3 / (c - 3 / 2) * (vM n * (-X4 n) * e) := by
    unfold pterm
    rw [hg]
    ring
  have b1 : |vM n * (-X4 n) / (n : ℝ) ^ 2| ≤ vM n / (n : ℝ) ^ 2 := by
    rw [abs_div, abs_mul, abs_of_nonneg hΛ, abs_neg, abs_of_pos (by positivity : (0:ℝ) < (n:ℝ) ^ 2)]
    exact div_le_div_of_nonneg_right (mul_le_of_le_one_right hΛ hX) (by positivity)
  have b2 : |vM n * (-X4 n) * e| ≤ vM n / (n : ℝ) ^ 2 := by
    rw [abs_mul, abs_mul, abs_of_nonneg hΛ, abs_neg, abs_of_nonneg he0]
    calc vM n * |X4 n| * e ≤ vM n * 1 * (1 / (n : ℝ) ^ 2) :=
          mul_le_mul (mul_le_mul_of_nonneg_left hX hΛ) he1 he0 (by positivity)
      _ = vM n / (n : ℝ) ^ 2 := by ring
  constructor
  · rw [e1]
    calc |2 * c / (c - 3 / 2) * (vM n * (-X4 n) / (n : ℝ) ^ 2)
          - 3 / (c - 3 / 2) * (vM n * (-X4 n) * e)|
        ≤ |2 * c / (c - 3 / 2) * (vM n * (-X4 n) / (n : ℝ) ^ 2)|
          + |3 / (c - 3 / 2) * (vM n * (-X4 n) * e)| := abs_sub _ _
      _ = 2 * c / (c - 3 / 2) * |vM n * (-X4 n) / (n : ℝ) ^ 2|
          + 3 / (c - 3 / 2) * |vM n * (-X4 n) * e| := by
          rw [abs_mul, abs_mul, abs_of_nonneg hK1, abs_of_nonneg hK2]
      _ ≤ 2 * c / (c - 3 / 2) * (vM n / (n : ℝ) ^ 2) + 3 / (c - 3 / 2) * (vM n / (n : ℝ) ^ 2) :=
          add_le_add (mul_le_mul_of_nonneg_left b1 hK1) (mul_le_mul_of_nonneg_left b2 hK2)
      _ = _ := by ring
  · rw [e1]
    have := mul_le_mul_of_nonneg_left (neg_le_of_abs_le b2) hK2
    linarith

theorem summable_pterm {c : ℝ} (hc : 3 / 2 < c) : Summable (pterm c) :=
  Summable.of_norm_bounded (summable_vM_div_sq.mul_left _) fun n => by
    rw [Real.norm_eq_abs]; exact (pterm_le hc n).1

theorem abs_gN_log_div_le {c : ℝ} (hc : 3 / 2 < c) {n : ℕ} (hn : 0 < n) :
    |gN c (Real.log n)| / Real.sqrt n ≤ (c + 3 / 2) / (c - 3 / 2) * (1 / (n : ℝ) ^ 2) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
  have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
  have hlog : 0 ≤ Real.log n := Real.log_nonneg h1
  have hcd : 0 < c - 3 / 2 := by linarith
  have hE2 : Real.exp (-c * Real.log n) ≤ Real.exp (-(3 / 2) * Real.log n) :=
    Real.exp_le_exp.mpr (by nlinarith)
  have hE20 := Real.exp_pos (-c * Real.log n)
  have hE10 := Real.exp_pos (-(3 / 2) * Real.log n)
  have hg : |gN c (Real.log n)| ≤ (c + 3 / 2) / (c - 3 / 2) * Real.exp (-(3 / 2) * Real.log n) := by
    unfold gN
    rw [abs_of_nonneg hlog, abs_div, abs_of_pos hcd, div_le_iff₀ hcd]
    have : (c + 3 / 2) / (c - 3 / 2) * Real.exp (-(3 / 2) * Real.log n) * (c - 3 / 2)
        = (c + 3 / 2) * Real.exp (-(3 / 2) * Real.log n) := by
      rw [mul_assoc, mul_comm (Real.exp _), ← mul_assoc, div_mul_cancel₀ _ hcd.ne']
    rw [this, abs_le]
    constructor <;> nlinarith
  rw [div_le_iff₀ hs]
  refine hg.trans (le_of_eq ?_)
  rw [exp_neg_three_halves_log hn]
  generalize Real.sqrt n = s at hs hss
  rw [hss]
  field_simp

/-- The explicit formula for `gN c` under GRH, as an upper bound. -/
theorem EF_gN_le (hGRH : GRHFor chi4) {c : ℝ} (hc : 20 ≤ c) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * gNhat c (ρ : ℂ).im) ∧
      ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * gNhat c (ρ : ℂ).im
        ≤ c / (c - 3 / 2) * (2 - Real.eulerMascheroniConstant - Real.log Real.pi + 2 * Aser)
          + 3 / (c - 3 / 2) * (7 / 6) := by
  have hc' : 3 / 2 < c := by linarith
  have hcd : 0 < c - 3 / 2 := by linarith
  -- hypotheses of `EF_ext`
  have hps : Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖((gN c (Real.log n) : ℝ) : ℂ)‖ + ‖((gN c (-Real.log n) : ℝ) : ℂ)‖)) := by
    refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_)
      (summable_vM_div_sq.mul_left (2 * ((c + 3 / 2) / (c - 3 / 2))))
    · exact mul_nonneg (div_nonneg (vM_nonneg n) (Real.sqrt_nonneg _)) (by positivity)
    · rcases Nat.eq_zero_or_pos n with rfl | hn
      · simp
      · rw [gN_neg, Complex.norm_real, Real.norm_eq_abs]
        have h := abs_gN_log_div_le hc' hn
        have hΛ := vM_nonneg n
        calc ArithmeticFunction.vonMangoldt n / Real.sqrt n
              * (|gN c (Real.log n)| + |gN c (Real.log n)|)
            = 2 * vM n * (|gN c (Real.log n)| / Real.sqrt n) := by unfold vM; ring
          _ ≤ 2 * vM n * ((c + 3 / 2) / (c - 3 / 2) * (1 / (n : ℝ) ^ 2)) :=
              mul_le_mul_of_nonneg_left h (by positivity)
          _ = _ := by ring
  obtain ⟨hsum, hEF⟩ := EF_ext (χ := chi4) one_lt_four RaceLib.chi4_isPrimitive hGRH
    (k := fun u => ((gN c u : ℝ) : ℂ)) (k1 := fun u => ((gN1 c u : ℝ) : ℂ))
    (k2 := fun u => ((gN2 c u : ℝ) : ℂ)) (fun u => (hasDerivAt_gN c u).ofReal_comp)
    (fun u => (hasDerivAt_gN1 c u).ofReal_comp)
    (Complex.continuous_ofReal.comp (continuous_gN2 c)) (integrable_gN hc').ofReal
    (integrable_gN1 hc').ofReal (integrable_gN2 hc').ofReal hps
  have hFT : ∀ τ : ℝ, paperFT (fun u => ((gN c u : ℝ) : ℂ)) τ = ((gNhat c τ : ℝ) : ℂ) := fun τ =>
    integral_gN_mul_cexp hc' τ
  have hterm : ∀ ρ : Z4.carrier,
      (Z4.mult ρ : ℂ) * paperFT (fun u => ((gN c u : ℝ) : ℂ)) (gammaOf ρ)
        = (((Z4.mult ρ : ℝ) * gNhat c (ρ : ℂ).im : ℝ) : ℂ) := by
    intro ρ
    rw [Z4_gammaOf hGRH ρ, hFT]
    push_cast
    ring
  have hsumR : Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * gNhat c (ρ : ℂ).im) :=
    Complex.summable_ofReal.mp (hsum.congr hterm)
  refine ⟨hsumR, ?_⟩
  -- the archimedean side
  obtain ⟨-, a, ha, harch⟩ := arch_gN_le hc
  have harch' : (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, paperFT (fun u => ((gN c u : ℝ) : ℂ)) r
      * ((((Complex.digamma (1 / 4 + (parity chi4 : ℂ) / 2 + I * r / 2)).re
        + Real.log ((4 : ℕ) / Real.pi) : ℝ)) : ℂ) = (a : ℂ) := by
    rw [← harch]
    congr 2
    funext r
    have hw : (Complex.digamma (1 / 4 + (parity chi4 : ℂ) / 2 + I * r / 2)).re
        + Real.log ((4 : ℕ) / Real.pi)
        = (Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi) :=
      RaceLib.archWeight_chi4 r
    rw [hw]
    rfl
  -- the prime side
  have hprime : ∀ n : ℕ, ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
        * (coeff chi4 n * ((gN c (Real.log n) : ℝ) : ℂ)
          + (starRingEnd ℂ) (coeff chi4 n) * ((gN c (-Real.log n) : ℝ) : ℂ))
      = ((-(pterm c n) : ℝ) : ℂ) := by
    intro n
    rw [coeff_chi4, Complex.conj_ofReal, gN_neg]
    unfold pterm vM
    push_cast
    ring
  have hreal : ((∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * gNhat c (ρ : ℂ).im : ℝ) : ℂ)
      = ((∑' n : ℕ, pterm c n + a : ℝ) : ℂ) := by
    rw [Complex.ofReal_tsum, ← tsum_congr hterm, hEF]
    unfold literatureRHSchi
    rw [harch', tsum_congr hprime, ← Complex.ofReal_tsum, tsum_neg]
    push_cast
    ring
  rw [Complex.ofReal_injective hreal]
  -- bound the prime side
  have hA : Summable (fun n : ℕ => vM n * (-X4 n) / (n : ℝ) ^ 2) := summable_Aser
  have hB := summable_vM_div_sq
  have hp : ∑' n : ℕ, pterm c n ≤ 2 * c / (c - 3 / 2) * Aser + 3 / (c - 3 / 2) * (7 / 6) := by
    have h1 : ∑' n : ℕ, pterm c n
        ≤ ∑' n : ℕ, (2 * c / (c - 3 / 2) * (vM n * (-X4 n) / (n : ℝ) ^ 2)
          + 3 / (c - 3 / 2) * (vM n / (n : ℝ) ^ 2)) :=
      (summable_pterm hc').tsum_le_tsum (fun n => (pterm_le hc' n).2)
        ((hA.mul_left _).add (hB.mul_left _))
    rw [(hA.mul_left _).tsum_add (hB.mul_left _), tsum_mul_left, tsum_mul_left] at h1
    have h2 : ∑' n : ℕ, vM n / (n : ℝ) ^ 2 ≤ 7 / 6 := tsum_vonMangoldt_div_sq_le
    have h3 : (0 : ℝ) ≤ 3 / (c - 3 / 2) := by positivity
    have h4 : ∑' n : ℕ, vM n * (-X4 n) / (n : ℝ) ^ 2 = Aser := rfl
    rw [h4] at h1
    nlinarith [mul_le_mul_of_nonneg_left h2 h3]
  have e : c / (c - 3 / 2) * (2 - Real.eulerMascheroniConstant - Real.log Real.pi + 2 * Aser)
      = 2 * c / (c - 3 / 2) * Aser
        + c / (c - 3 / 2) * (2 - Real.eulerMascheroniConstant - Real.log Real.pi) := by ring
  rw [e]
  linarith

/-- `gNhat N r → 3/(9/4+r²)` as `N → ∞`. -/
theorem tendsto_gNhat (r : ℝ) :
    Tendsto (fun N : ℕ => gNhat (N : ℝ) r) atTop (𝓝 (3 / (9 / 4 + r ^ 2))) := by
  have h1 : Tendsto (fun N : ℕ => 1 + 3 / 2 / (N : ℝ)) atTop (𝓝 (1 + 0)) :=
    tendsto_const_nhds.add (tendsto_const_div_atTop_nhds_zero_nat _)
  have h2 : Tendsto (fun N : ℕ => 1 + r ^ 2 / (N : ℝ) * (1 / (N : ℝ))) atTop (𝓝 (1 + 0 * 0)) :=
    tendsto_const_nhds.add ((tendsto_const_div_atTop_nhds_zero_nat _).mul
      (tendsto_const_div_atTop_nhds_zero_nat _))
  have h3 := (h1.div h2 (by norm_num)).const_mul (3 / (9 / 4 + r ^ 2))
  rw [show 3 / (9 / 4 + r ^ 2) * ((1 + 0) / (1 + 0 * 0)) = 3 / (9 / 4 + r ^ 2) by norm_num] at h3
  refine h3.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have x0 : 0 ≤ r ^ 2 := sq_nonneg r
  unfold gNhat
  simp only [Pi.div_apply]
  have hNn : (N : ℝ) ≠ 0 := hN'.ne'
  have hA : (9 / 4 + r ^ 2 : ℝ) ≠ 0 := by positivity
  have hB : ((N : ℝ) ^ 2 + r ^ 2) ≠ 0 := by positivity
  have hC : (1 + r ^ 2 / (N : ℝ) * (1 / (N : ℝ))) ≠ 0 := by positivity
  rw [div_mul_div_comm, div_eq_div_iff (mul_ne_zero hA hC) (mul_ne_zero hA hB)]
  field_simp
  ring

/-- **`Σ_ρ m_ρ · 3/(9/4+γ_ρ²) ≤ 2 − γ − log π + 2A`** under GRH for `L(s,χ₋₄)`. -/
theorem sumOne_le (hGRH : GRHFor chi4) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2))) ∧
      ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2))
        ≤ 2 - Real.eulerMascheroniConstant - Real.log Real.pi + 2 * Aser := by
  set L : ℝ := 2 - Real.eulerMascheroniConstant - Real.log Real.pi + 2 * Aser with hL
  have hnn : ∀ ρ : Z4.carrier, 0 ≤ (Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2)) := fun ρ =>
    mul_nonneg (Nat.cast_nonneg _) (by positivity)
  have hfin : ∀ F : Finset Z4.carrier,
      ∑ ρ ∈ F, (Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2)) ≤ L := by
    intro F
    have hf : Tendsto (fun N : ℕ => ∑ ρ ∈ F, (Z4.mult ρ : ℝ) * gNhat (N : ℝ) (ρ : ℂ).im) atTop
        (𝓝 (∑ ρ ∈ F, (Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2)))) :=
      tendsto_finsetSum _ fun ρ _ => (tendsto_gNhat _).const_mul _
    have hg : Tendsto (fun N : ℕ => (1 / (1 - 3 / 2 / (N : ℝ))) * L
        + (3 / (N : ℝ)) / (1 - 3 / 2 / (N : ℝ)) * (7 / 6)) atTop
        (𝓝 ((1 / (1 - 0)) * L + 0 / (1 - 0) * (7 / 6))) := by
      have hd : Tendsto (fun N : ℕ => 1 - 3 / 2 / (N : ℝ)) atTop (𝓝 (1 - 0)) :=
        tendsto_const_nhds.sub (tendsto_const_div_atTop_nhds_zero_nat _)
      exact ((tendsto_const_nhds.div hd (by norm_num)).mul_const L).add
        (((tendsto_const_div_atTop_nhds_zero_nat 3).div hd (by norm_num)).mul_const _)
    rw [show (1 / (1 - 0)) * L + 0 / (1 - 0) * (7 / 6) = L by norm_num] at hg
    refine le_of_tendsto_of_tendsto hf hg ?_
    filter_upwards [eventually_ge_atTop 20] with N hN
    have hN' : (20 : ℝ) ≤ N := by exact_mod_cast hN
    obtain ⟨hs, hle⟩ := EF_gN_le hGRH hN'
    have h1 : ∑ ρ ∈ F, (Z4.mult ρ : ℝ) * gNhat (N : ℝ) (ρ : ℂ).im
        ≤ ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * gNhat (N : ℝ) (ρ : ℂ).im :=
      hs.sum_le_tsum F fun ρ _ =>
        mul_nonneg (Nat.cast_nonneg _) (gNhat_nonneg (by linarith) _)
    refine h1.trans (hle.trans (le_of_eq ?_))
    have hN0 : (N : ℝ) ≠ 0 := by linarith
    have hN1 : (N : ℝ) - 3 / 2 ≠ 0 := by linarith
    have hN2 : 1 - 3 / 2 / (N : ℝ) ≠ 0 := by
      rw [sub_ne_zero]; intro h; rw [eq_div_iff hN0] at h; linarith
    rw [← hL]
    field_simp
  exact ⟨summable_of_sum_le hnn hfin, Real.tsum_le_of_sum_le hnn hfin⟩

/-- **`Σ_ρ m_ρ · 30/((9/4+γ_ρ²)(25/4+γ_ρ²)) = 7/2 − γ − 3 log 2 − log π + 5A − 3B`** under GRH. -/
theorem sumMaj_eq (hGRH : GRHFor chi4) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * majF (ρ : ℂ).im) ∧
      ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * majF (ρ : ℂ).im
        = 7 / 2 - Real.eulerMascheroniConstant - 3 * Real.log 2 - Real.log Real.pi
          + 5 * Aser - 3 * Bser := by
  obtain ⟨hs, he⟩ := EF_maj hGRH
  refine ⟨hs, ?_⟩
  have hA := summable_Aser
  have hB := summable_Bser
  have ht : ∀ n : ℕ, ArithmeticFunction.vonMangoldt n * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ))
        * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)
      = 5 * (ArithmeticFunction.vonMangoldt n * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) / (n : ℝ) ^ 2)
        - 3 * (ArithmeticFunction.vonMangoldt n * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) / (n : ℝ) ^ 3) := by
    intro n
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp
    · have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
      field_simp
      ring
  have : ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * majF (ρ : ℂ).im
      = ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im := rfl
  rw [this, he, tsum_congr ht, (hA.mul_left 5).tsum_sub (hB.mul_left 3), tsum_mul_left,
    tsum_mul_left]
  unfold Aser Bser
  ring

end Allx
