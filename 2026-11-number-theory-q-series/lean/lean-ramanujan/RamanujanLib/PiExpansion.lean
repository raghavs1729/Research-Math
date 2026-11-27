/-
RamanujanLib/PiExpansion.lean — Ramanujan's five-term expansion of `π(x)`:

  `|π(x) − x (1/L + 1/L² + 2/L³ + 6/L⁴ + 24/L⁵)| ≤ C · x / L⁶`,  `L = log x`,  for all `x ≥ 2`.

Partial summation `π(x) = θ(x)/log x + ∫_2^x θ(t) / (t log² t) dt` (Mathlib), the bound
`|θ(t) − t| ≤ K t / log⁵ t` (`ThetaBound`), and the expansion of `∫_2^x dt / log² t`
(`LogIntegral`).
-/
import RamanujanLib.LogIntegral
import RamanujanLib.ThetaBound

open Real Set MeasureTheory

noncomputable section

namespace RamanujanLib

theorem exists_primeCounting_expansion :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ) - x * P (log x)⁻¹| ≤ C * (x / log x ^ 6) := by
  obtain ⟨K, hK0, hK⟩ := exists_theta_sub_le
  obtain ⟨C₁, hC₁⟩ := exists_li_expansion
  obtain ⟨C₇, hC₇0, hC₇⟩ := exists_J_le 7
  refine ⟨K + C₁ + K * C₇ / log 2, fun x hx => ?_⟩
  have hl : 0 < log x := log_pos_of_two_le hx
  have hl2 : 0 < log (2 : ℝ) := log_pos (by norm_num)
  have hx0 : 0 < x := by linarith
  have hq : 0 ≤ x / log x ^ 6 := by positivity
  -- the error integral
  set R : ℝ := ∫ t in (2 : ℝ)..x, (Chebyshev.theta t / (t * log t ^ 2) - 1 / log t ^ 2) with hR
  have hθint : IntervalIntegrable (fun t : ℝ => Chebyshev.theta t / (t * log t ^ 2)) volume 2 x := by
    rw [intervalIntegrable_iff, uIoc_of_le hx, ← integrableOn_Icc_iff_integrableOn_Ioc]
    exact Chebyshev.integrableOn_theta_div_id_mul_log_sq x
  have h1int : IntervalIntegrable (fun t : ℝ => 1 / log t ^ 2) volume 2 x :=
    intervalIntegrable_inv_log_pow 2 (by norm_num) (by linarith)
  have hRsplit : R = (∫ t in (2 : ℝ)..x, Chebyshev.theta t / (t * log t ^ 2)) - J 2 x := by
    rw [hR, intervalIntegral.integral_sub hθint h1int, J]
  have hRbound : |R| ≤ K * J 7 x := by
    have h7int : IntervalIntegrable (fun t : ℝ => K * (1 / log t ^ 7)) volume 2 x :=
      (intervalIntegrable_inv_log_pow 7 (by norm_num) (by linarith)).const_mul K
    have := intervalIntegral.norm_integral_le_of_norm_le (μ := volume) hx
      (f := fun t : ℝ => Chebyshev.theta t / (t * log t ^ 2) - 1 / log t ^ 2)
      (g := fun t : ℝ => K * (1 / log t ^ 7))
      (Filter.Eventually.of_forall fun t ht => ?_) h7int
    · rw [Real.norm_eq_abs, intervalIntegral.integral_const_mul] at this
      exact this
    · have ht2 : 2 ≤ t := ht.1.le
      have htl : 0 < log t := log_pos_of_two_le ht2
      have ht0 : 0 < t := by linarith
      have hE := hK t ht2
      have heq : Chebyshev.theta t / (t * log t ^ 2) - 1 / log t ^ 2
          = (Chebyshev.theta t - t) / (t * log t ^ 2) := by
        field_simp
      rw [Real.norm_eq_abs, heq, abs_div, abs_of_pos (by positivity : 0 < t * log t ^ 2),
        div_le_iff₀ (by positivity)]
      calc |Chebyshev.theta t - t| ≤ K * (t / log t ^ 5) := hE
        _ = K * (1 / log t ^ 7) * (t * log t ^ 2) := by
            field_simp
  -- assemble
  have hπ := Chebyshev.primeCounting_eq_theta_div_log_add_integral hx
  have hli := hC₁ x hx
  have hEx := hK x hx
  have hJ7 := hC₇ x hx
  have hmain : (Nat.primeCounting ⌊x⌋₊ : ℝ) - x * P (log x)⁻¹
      = (Chebyshev.theta x - x) / log x + (x / log x + J 2 x - x * P (log x)⁻¹) + R := by
    rw [hπ, hRsplit]
    field_simp
    ring
  have hterm1 : |(Chebyshev.theta x - x) / log x| ≤ K * (x / log x ^ 6) := by
    rw [abs_div, abs_of_pos hl, div_le_iff₀ hl]
    calc |Chebyshev.theta x - x| ≤ K * (x / log x ^ 5) := hEx
      _ = K * (x / log x ^ 6) * log x := by
          field_simp
  have hterm3 : |R| ≤ K * C₇ / log 2 * (x / log x ^ 6) := by
    have hlog : log 2 ≤ log x := log_two_le_log hx
    have h76 : x / log x ^ 7 ≤ (x / log x ^ 6) / log 2 := by
      rw [div_div, pow_succ]
      gcongr
    calc |R| ≤ K * J 7 x := hRbound
      _ ≤ K * (C₇ * (x / log x ^ 7)) := by gcongr
      _ ≤ K * (C₇ * ((x / log x ^ 6) / log 2)) := by gcongr
      _ = K * C₇ / log 2 * (x / log x ^ 6) := by ring
  rw [hmain, add_mul, add_mul]
  exact (abs_add_three _ _ _).trans (add_le_add (add_le_add hterm1 hli) hterm3)

/-- The same expansion with the main term written out. -/
theorem exists_primeCounting_expansion' :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ)
          - x * (1 / log x + 1 / log x ^ 2 + 2 / log x ^ 3 + 6 / log x ^ 4 + 24 / log x ^ 5)|
        ≤ C * (x / log x ^ 6) := by
  obtain ⟨C, hC⟩ := exists_primeCounting_expansion
  refine ⟨C, fun x hx => ?_⟩
  have hP : P (log x)⁻¹
      = 1 / log x + 1 / log x ^ 2 + 2 / log x ^ 3 + 6 / log x ^ 4 + 24 / log x ^ 5 := by
    unfold P
    ring
  rw [← hP]
  exact hC x hx

end RamanujanLib
