/-
ChebLib/Cheb/Headline.lean — the results of ChebLib/Cheb in the notation of the statement layer
(ChebLib/Race/Defs.lean, ChebLib/Cheb/Defs.lean): bridges from the Zeta23 objects to Mathlib-only
definitions, the window form of the logarithmic statement, and Chebyshev's literal form.
-/
import ChebLib.Cheb.Defs
import ChebLib.Cheb.Final

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace ChebLib

/-! ### `χ₋₄` -/

theorem chi4_apply_natCast (n : ℕ) :
    chi4 (n : ZMod 4) = (((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) : ℂ) := RaceLib.chi4_natCast n

theorem chi4_natCast_two : ((ZMod.χ₄ ((2 : ℕ) : ZMod 4) : ℤ) : ℝ) = 0 := RaceLib.chi4_natCast_two

theorem chi4_natCast_of_odd {n : ℕ} (hn : n % 2 = 1) :
    ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = (-1) ^ ((n - 1) / 2) := by
  rw [ZMod.χ₄_eq_neg_one_pow hn]
  have h : n / 2 = (n - 1) / 2 := by omega
  rw [h]
  push_cast
  rfl

/-! ### The explicit formula for integrable test functions -/

/-- The extended explicit formula in the notation of the statement layer. -/
theorem zeroSum_eq_of_integrable {q : ℕ} [NeZero q] {χ : DirichletCharacter ℂ q} (hq : 1 < q)
    (hprim : χ.IsPrimitive) (hGRH : GRHFor χ)
    {k k1 k2 : ℝ → ℂ} (hk1 : ∀ u, HasDerivAt k (k1 u) u) (hk2 : ∀ u, HasDerivAt k1 (k2 u) u)
    (hk2c : Continuous k2) (hki : Integrable k) (hk1i : Integrable k1) (hk2i : Integrable k2)
    (hps : Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖k (Real.log n)‖ + ‖k (-Real.log n)‖))) :
    Summable (fun ρ : {ρ : ℂ | χ.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1} =>
        ((analyticOrderAt χ.LFunction (ρ : ℂ)).toNat : ℂ) * raceFT k (((ρ : ℂ) - 1 / 2) / I)) ∧
      zeroSum χ k = archTerm χ k 0 - primePowerSumSym χ k := by
  obtain ⟨hsum, hEF⟩ := EF_ext hq hprim hGRH hk1 hk2 hk2c hki hk1i hk2i hps
  refine ⟨hsum, ?_⟩
  rw [← RaceLib.literatureRHSchi_eq, ← hEF]
  rfl

/-! ### The zero sum -/

theorem Z4_eq (hGRH : GRHFor chi4) (ρ : Z4.carrier) :
    (ρ : ℂ) = 1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ) := by
  apply Complex.ext
  · simp [Z4_re hGRH ρ]
  · simp

theorem gammaZeroSum_chi4_eq (hGRH : GRHFor chi4) : gammaZeroSum chi4 = gammaSum := by
  change (∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * ‖Complex.Gamma (ρ : ℂ)‖) = gammaSum
  unfold gammaSum
  refine tsum_congr fun ρ => ?_
  rw [← Z4_eq hGRH ρ]

theorem summable_gammaZeroSum_chi4 (hGRH : GRHFor chi4) :
    Summable (fun ρ : {ρ : ℂ | chi4.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1} =>
      ((analyticOrderAt chi4.LFunction (ρ : ℂ)).toNat : ℝ) * ‖Complex.Gamma (ρ : ℂ)‖) := by
  change Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * ‖Complex.Gamma (ρ : ℂ)‖)
  refine (summable_gammaSum hGRH).congr fun ρ => ?_
  rw [← Z4_eq hGRH ρ]

theorem gammaZeroSum_chi4_lt (hGRH : GRHFor chi4) :
    gammaZeroSum chi4 < Real.sqrt Real.pi / 2 := by
  rw [gammaZeroSum_chi4_eq hGRH]
  exact gammaSum_lt hGRH

/-! ### The three Abel sums -/

theorem chebyshevPsi_eq (x : ℝ) : chebyshevPsi x = psiAbel x := by
  unfold chebyshevPsi psiAbel
  refine tsum_congr fun n => ?_
  rw [neg_div]
  rfl

theorem chebyshevLogSum_eq (x : ℝ) : chebyshevLogSum x = phiAbel x := tsum_primes_chi4_log_exp x

theorem chebyshevSum_eq (x : ℝ) : chebyshevSum x = ∑' n : ℕ, chebA n * Real.exp (-(n / x)) :=
  tsum_primes_chi4_exp x

theorem summable_chebyshevPsi {x : ℝ} (hx : 0 < x) :
    Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(n : ℝ) / x)) := by
  refine (summable_psiAbel hx).congr fun n => ?_
  rw [neg_div]
  rfl

theorem summable_chebyshevLogSum {x : ℝ} (hx : 0 < x) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  have h : Summable (fun n : ℕ => ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * Real.log n
      * Real.exp (-(n : ℝ) / x)) := by
    have := summable_mul_exp_neg_div (fun n : ℕ => X4 n * Real.log n) 1 (fun n => by
      rw [one_mul, abs_mul, abs_of_nonneg (Real.log_natCast_nonneg n)]
      calc |X4 n| * Real.log n ≤ 1 * Real.log n :=
            mul_le_mul_of_nonneg_right (abs_X4_le_one n) (Real.log_natCast_nonneg n)
        _ ≤ n := by rw [one_mul]; exact Real.log_le_self (Nat.cast_nonneg n)) hx
    refine this.congr fun n => ?_
    rw [neg_div]
    rfl
  exact h.subtype {n : ℕ | n.Prime}

theorem abs_chebyshevPsi_le (hGRH : GRHFor chi4) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ x : ℝ in atTop, |chebyshevPsi x| ≤ (gammaZeroSum chi4 + ε) * Real.sqrt x := by
  rw [gammaZeroSum_chi4_eq hGRH]
  filter_upwards [abs_psiAbel_le hGRH hε] with x hx
  rw [chebyshevPsi_eq]
  exact hx

theorem abs_chebyshevLogSum_add_le (hGRH : GRHFor chi4) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ x : ℝ in atTop, |chebyshevLogSum x + Real.sqrt Real.pi / 2 * Real.sqrt x|
      ≤ (gammaZeroSum chi4 + ε) * Real.sqrt x := by
  rw [gammaZeroSum_chi4_eq hGRH]
  filter_upwards [abs_phiAbel_add_le hGRH hε] with x hx
  rw [chebyshevLogSum_eq]
  exact hx

theorem chebyshevLogSum_le (hGRH : GRHFor chi4) :
    ∃ c : ℝ, 0 < c ∧ ∀ᶠ x : ℝ in atTop, chebyshevLogSum x ≤ -c * Real.sqrt x := by
  obtain ⟨c, hc, x₀, h⟩ := phiAbel_le hGRH
  refine ⟨c, hc, ?_⟩
  filter_upwards [eventually_ge_atTop x₀] with x hx
  rw [chebyshevLogSum_eq]
  exact h x hx

/-! ### The window form: a smoothed race with a concrete window -/

/-- With the window `g(u) = e^{u/2} exp(−e^u)`:
`Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L) = e^{−L/2} Σ_p χ₋₄(p) (log p) e^{−p/e^L}`. -/
theorem tsum_window_eq (L : ℝ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
        * (Real.log (p : ℕ) / Real.sqrt (p : ℕ))
        * (Real.exp ((Real.log (p : ℕ) - L) / 2) * Real.exp (-Real.exp (Real.log (p : ℕ) - L)))
      = chebyshevLogSum (Real.exp L) / Real.sqrt (Real.exp L) := by
  unfold chebyshevLogSum
  rw [← tsum_div_const]
  refine tsum_congr fun p => ?_
  have hp : (0 : ℝ) < ((p : ℕ) : ℝ) := by exact_mod_cast p.2.pos
  have hs : 0 < Real.sqrt ((p : ℕ) : ℝ) := Real.sqrt_pos.mpr hp
  have hE : 0 < Real.sqrt (Real.exp L) := Real.sqrt_pos.mpr (Real.exp_pos L)
  have e1 : Real.exp ((Real.log (p : ℕ) - L) / 2)
      = Real.sqrt ((p : ℕ) : ℝ) / Real.sqrt (Real.exp L) := by
    rw [show (Real.log (p : ℕ) - L) / 2 = Real.log (p : ℕ) / 2 - L / 2 by ring, Real.exp_sub,
      Real.exp_half, Real.exp_half, Real.exp_log hp]
  have e2 : Real.exp (Real.log (p : ℕ) - L) = ((p : ℕ) : ℝ) / Real.exp L := by
    rw [Real.exp_sub, Real.exp_log hp]
  rw [e1, e2, neg_div]
  field_simp

theorem window_sum_le (hGRH : GRHFor chi4) :
    ∃ c : ℝ, 0 < c ∧ ∀ᶠ L : ℝ in atTop,
      ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
        * (Real.log (p : ℕ) / Real.sqrt (p : ℕ))
        * (Real.exp ((Real.log (p : ℕ) - L) / 2) * Real.exp (-Real.exp (Real.log (p : ℕ) - L)))
      ≤ -c := by
  obtain ⟨c, hc, h⟩ := chebyshevLogSum_le hGRH
  refine ⟨c, hc, ?_⟩
  filter_upwards [Real.tendsto_exp_atTop.eventually h] with L hL
  rw [tsum_window_eq]
  have hE : 0 < Real.sqrt (Real.exp L) := Real.sqrt_pos.mpr (Real.exp_pos L)
  rw [div_le_iff₀ hE]
  exact hL

/-! ### Chebyshev's assertion -/

theorem tendsto_chebyshevSum (hGRH : GRHFor chi4) : Tendsto chebyshevSum atTop atBot :=
  chebyshev_bias hGRH

/-- Chebyshev's literal form: `e^{−3c} − e^{−5c} + e^{−7c} + e^{−11c} − e^{−13c} − ⋯ → +∞` as
`c → 0⁺` (sign `(−1)^{(p+1)/2}` at the odd prime `p`). -/
theorem tendsto_chebyshev_literal (hGRH : GRHFor chi4) :
    Tendsto (fun c : ℝ => ∑' p : Nat.Primes,
      if (p : ℕ) = 2 then 0 else (-1 : ℝ) ^ (((p : ℕ) + 1) / 2) * Real.exp (-((p : ℕ) : ℝ) * c))
      (𝓝[>] 0) atTop := by
  have h := (tendsto_neg_atBot_atTop.comp (tendsto_chebyshevSum hGRH)).comp
    tendsto_inv_nhdsGT_zero
  refine h.congr fun c => ?_
  simp only [Function.comp_apply]
  unfold chebyshevSum
  rw [← tsum_neg]
  refine tsum_congr fun p => ?_
  by_cases h2 : (p : ℕ) = 2
  · rw [if_pos h2, h2, chi4_natCast_two, zero_mul, neg_zero]
  · rw [if_neg h2]
    have hodd : (p : ℕ) % 2 = 1 := (p.2.eq_two_or_odd).resolve_left h2
    rw [chi4_natCast_of_odd hodd, div_inv_eq_mul, neg_mul]
    have e : ((p : ℕ) + 1) / 2 = ((p : ℕ) - 1) / 2 + 1 := by omega
    rw [e, pow_succ]
    ring

end ChebLib
