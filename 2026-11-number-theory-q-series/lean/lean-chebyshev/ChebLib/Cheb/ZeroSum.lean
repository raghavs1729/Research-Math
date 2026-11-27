/-
ChebLib/Cheb/ZeroSum.lean — the zero sum of the Gamma factor is small:
under GRH for `L(s,χ₋₄)`,   `Σ_ρ m_ρ |Γ(1/2 + iγ_ρ)| < √π / 2`.

Proof.  The explicit formula (`EF_ext`) for the majorant
`g(u) = (5/2) e^{−(3/2)|u|} − (3/2) e^{−(5/2)|u|}`, whose transform is
`ĝ(r) = 30/((9/4+r²)(25/4+r²)) > 0`, reads

    Σ_ρ m_ρ ĝ(γ_ρ) = (7/2 − γ − 3 log 2 − log π) + Σ_n Λ(n) (−χ₋₄(n)) (5n − 3)/n³  < 16/15

(archimedean side: `Poisson.lean`; numerical inequality: `Numeric.lean`), and
`|Γ(1/2 + ir)| ≤ (15/32) √π ĝ(r)` (`Kernel.lean`).  No information on individual zeros is used.
-/
import ChebLib.Cheb.PsiAbel
import ChebLib.Cheb.Majorant
import ChebLib.Cheb.Poisson
import ChebLib.Cheb.Numeric

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace ChebLib

/-- `ĝ(r) = 30/((9/4+r²)(25/4+r²))`. -/
def majFT (r : ℝ) : ℝ := 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))

theorem majFT_nonneg (r : ℝ) : 0 ≤ majFT r := by unfold majFT; positivity

theorem rpow_neg_three_halves {t : ℝ} (ht : 0 < t) :
    t ^ (-(3 / 2) : ℝ) = (t * Real.sqrt t)⁻¹ := by
  rw [Real.rpow_neg ht.le, show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.rpow_add ht,
    Real.rpow_one, ← Real.sqrt_eq_rpow]

theorem rpow_neg_five_halves {t : ℝ} (ht : 0 < t) :
    t ^ (-(5 / 2) : ℝ) = (t ^ 2 * Real.sqrt t)⁻¹ := by
  rw [Real.rpow_neg ht.le, show (5 / 2 : ℝ) = 2 + 1 / 2 by norm_num, Real.rpow_add ht,
    Real.rpow_two, ← Real.sqrt_eq_rpow]

/-- The `n`-th term of the prime side of the explicit formula for the majorant. -/
theorem primeTerm_maj (n : ℕ) :
    ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
        * (coeff chi4 n * ((maj (Real.log n) : ℝ) : ℂ)
          + (starRingEnd ℂ) (coeff chi4 n) * ((maj (-Real.log n) : ℝ) : ℂ))
      = ((-(ArithmeticFunction.vonMangoldt n * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ))
          * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)) : ℝ) : ℂ) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
    have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
    have key : ArithmeticFunction.vonMangoldt n / Real.sqrt n
          * (X4 n * maj (Real.log n) + X4 n * maj (Real.log n))
        = -(ArithmeticFunction.vonMangoldt n * (-X4 n) * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)) := by
      rw [maj_log hn, rpow_neg_three_halves hn', rpow_neg_five_halves hn']
      generalize Real.sqrt n = s at hs hss
      rw [hss]
      field_simp
      ring
    rw [coeff_chi4, Complex.conj_ofReal, maj_neg]
    change _ = ((-(ArithmeticFunction.vonMangoldt n * (-X4 n)
      * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)) : ℝ) : ℂ)
    rw [← key]
    push_cast
    ring

theorem abs_maj_log_le {n : ℕ} (hn : 0 < n) :
    |maj (Real.log n)| ≤ 5 / 2 * ((n : ℝ) * Real.sqrt n)⁻¹ := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
  rw [maj_log hn, rpow_neg_three_halves hn', rpow_neg_five_halves hn']
  have h1 : ((n : ℝ) ^ 2 * Real.sqrt n)⁻¹ ≤ ((n : ℝ) * Real.sqrt n)⁻¹ := by
    apply inv_anti₀ (by positivity)
    have : (n : ℝ) ≤ (n : ℝ) ^ 2 := by nlinarith
    exact mul_le_mul_of_nonneg_right this hs.le
  have h2 : 0 < ((n : ℝ) ^ 2 * Real.sqrt n)⁻¹ := by positivity
  rw [abs_le]
  constructor <;> linarith

theorem summable_primeSide_maj :
    Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖((maj (Real.log n) : ℝ) : ℂ)‖ + ‖((maj (-Real.log n) : ℝ) : ℂ)‖)) := by
  refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_) (summable_vM_div_sq.mul_left 5)
  · exact mul_nonneg (div_nonneg (vM_nonneg n) (Real.sqrt_nonneg _)) (by positivity)
  · rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp
    · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
      have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
      have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
      rw [maj_neg, Complex.norm_real, Real.norm_eq_abs]
      have h := abs_maj_log_le hn
      calc ArithmeticFunction.vonMangoldt n / Real.sqrt n
            * (|maj (Real.log n)| + |maj (Real.log n)|)
          ≤ ArithmeticFunction.vonMangoldt n / Real.sqrt n
            * (5 / 2 * ((n : ℝ) * Real.sqrt n)⁻¹ + 5 / 2 * ((n : ℝ) * Real.sqrt n)⁻¹) :=
            mul_le_mul_of_nonneg_left (add_le_add h h)
              (div_nonneg (vM_nonneg n) (Real.sqrt_nonneg _))
        _ = 5 * (vM n / (n : ℝ) ^ 2) := by
            unfold vM
            generalize Real.sqrt n = s at hs hss
            rw [hss]
            field_simp
            ring

/-- **The explicit formula for the majorant, under GRH**: the positive zero sum
`Σ_ρ m_ρ ĝ(γ_ρ)` converges and equals
`(7/2 − γ − 3 log 2 − log π) + Σ_n Λ(n) (−χ₋₄(n)) (5n − 3)/n³`. -/
theorem EF_maj (hGRH : GRHFor chi4) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im) ∧
      ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im
        = 7 / 2 - Real.eulerMascheroniConstant - 3 * Real.log 2 - Real.log Real.pi
          + ∑' n : ℕ, ArithmeticFunction.vonMangoldt n
              * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3) := by
  obtain ⟨g1, g2, hd1, hd2, hc2, hi0, hi1, hi2⟩ := maj_data
  obtain ⟨hsum, hEF⟩ := EF_ext (χ := chi4) one_lt_four RaceLib.chi4_isPrimitive hGRH
    (k := fun u => ((maj u : ℝ) : ℂ)) (k1 := fun u => ((g1 u : ℝ) : ℂ))
    (k2 := fun u => ((g2 u : ℝ) : ℂ)) (fun u => (hd1 u).ofReal_comp) (fun u => (hd2 u).ofReal_comp)
    (Complex.continuous_ofReal.comp hc2) hi0.ofReal hi1.ofReal hi2.ofReal summable_primeSide_maj
  have hFT : ∀ τ : ℝ, paperFT (fun u => ((maj u : ℝ) : ℂ)) τ = ((majFT τ : ℝ) : ℂ) := fun τ =>
    integral_maj_mul_cexp τ
  have hterm : ∀ ρ : Z4.carrier,
      (Z4.mult ρ : ℂ) * paperFT (fun u => ((maj u : ℝ) : ℂ)) (gammaOf ρ)
        = (((Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im : ℝ) : ℂ) := by
    intro ρ
    rw [Z4_gammaOf hGRH ρ, hFT]
    push_cast
    ring
  have hsumR : Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im) :=
    Complex.summable_ofReal.mp (hsum.congr hterm)
  refine ⟨hsumR, ?_⟩
  apply Complex.ofReal_injective
  rw [Complex.ofReal_tsum, ← tsum_congr hterm, hEF]
  unfold literatureRHSchi
  have harch : (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, paperFT (fun u => ((maj u : ℝ) : ℂ)) r
      * ((((Complex.digamma (1 / 4 + (parity chi4 : ℂ) / 2 + I * r / 2)).re
        + Real.log ((4 : ℕ) / Real.pi) : ℝ)) : ℂ)
      = ((7 / 2 - Real.eulerMascheroniConstant - 3 * Real.log 2 - Real.log Real.pi : ℝ) : ℂ) := by
    have e : ∀ r : ℝ, paperFT (fun u => ((maj u : ℝ) : ℂ)) r
        * ((((Complex.digamma (1 / 4 + (parity chi4 : ℂ) / 2 + I * r / 2)).re
          + Real.log ((4 : ℕ) / Real.pi) : ℝ)) : ℂ)
        = ((30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))
          * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi)) : ℝ) : ℂ) := by
      intro r
      rw [hFT]
      have hw : (Complex.digamma (1 / 4 + (parity chi4 : ℂ) / 2 + I * r / 2)).re
          + Real.log ((4 : ℕ) / Real.pi)
          = (Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi) :=
        RaceLib.archWeight_chi4 r
      rw [hw]
      unfold majFT
      push_cast
      ring
    simp_rw [e]
    rw [integral_complex_ofReal, integral_majFT_mul_archWeight]
    have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
    push_cast
    field_simp
  rw [harch, tsum_congr primeTerm_maj, ← Complex.ofReal_tsum, tsum_neg]
  push_cast
  ring

/-- **The zero sum of the Gamma factor is small**: under GRH for `L(s,χ₋₄)`,
`Σ_ρ m_ρ |Γ(1/2 + iγ_ρ)| < √π / 2`. -/
theorem gammaSum_lt (hGRH : GRHFor chi4) : gammaSum < Real.sqrt Real.pi / 2 := by
  obtain ⟨hsum, hEF⟩ := EF_maj hGRH
  have hnum := maj_numeric
  rw [← hEF] at hnum
  have hle : gammaSum ≤ 15 / 32 * Real.sqrt Real.pi
      * ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im := by
    unfold gammaSum
    rw [← tsum_mul_left]
    refine (summable_gammaSum hGRH).tsum_le_tsum (fun ρ => ?_) (hsum.mul_left _)
    have h := norm_Gamma_half_mul_le (ρ : ℂ).im
    have hpos : (0:ℝ) < (9 / 4 + (ρ : ℂ).im ^ 2) * (25 / 4 + (ρ : ℂ).im ^ 2) := by positivity
    have hG : ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖
        ≤ 15 / 32 * Real.sqrt Real.pi * majFT (ρ : ℂ).im := by
      unfold majFT
      rw [show 15 / 32 * Real.sqrt Real.pi
          * (30 / ((9 / 4 + (ρ : ℂ).im ^ 2) * (25 / 4 + (ρ : ℂ).im ^ 2)))
          = 225 / 16 * Real.sqrt Real.pi
            / ((9 / 4 + (ρ : ℂ).im ^ 2) * (25 / 4 + (ρ : ℂ).im ^ 2)) by ring,
        le_div_iff₀ hpos]
      linarith
    calc (Z4.mult ρ : ℝ) * ‖Complex.Gamma (1 / 2 + I * (((ρ : ℂ).im : ℝ) : ℂ))‖
        ≤ (Z4.mult ρ : ℝ) * (15 / 32 * Real.sqrt Real.pi * majFT (ρ : ℂ).im) :=
          mul_le_mul_of_nonneg_left hG (Nat.cast_nonneg _)
      _ = _ := by ring
  have hsp : 0 < Real.sqrt Real.pi := Real.sqrt_pos.mpr Real.pi_pos
  calc gammaSum ≤ 15 / 32 * Real.sqrt Real.pi
        * ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * majFT (ρ : ℂ).im := hle
    _ < 15 / 32 * Real.sqrt Real.pi * (16 / 15) :=
        mul_lt_mul_of_pos_left hnum (by positivity)
    _ = Real.sqrt Real.pi / 2 := by ring

end ChebLib
