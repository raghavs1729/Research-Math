/-
RaceLib/Chi4.lean — the character χ₋₄: primitive, odd, and the S0 bound in closed form.
-/
import ChebLib.Race.S0

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

namespace RaceLib

theorem chi4_apply (a : ZMod 4) : chi4 a = ((ZMod.χ₄ a : ℤ) : ℂ) := rfl

theorem chi4_three : chi4 (3 : ZMod 4) = -1 := by
  rw [chi4_apply]; have : ZMod.χ₄ (3 : ZMod 4) = -1 := by decide
  rw [this]; simp

theorem chi4_two : chi4 (2 : ZMod 4) = 0 := by
  rw [chi4_apply]; have : ZMod.χ₄ (2 : ZMod 4) = 0 := by decide
  rw [this]; simp

theorem chi4_neg_one : chi4 (-1 : ZMod 4) = -1 := by
  have : (-1 : ZMod 4) = 3 := by decide
  rw [this, chi4_three]

theorem chi4_not_even : ¬ chi4.Even := by
  intro h
  have h1 : chi4 (-1 : ZMod 4) = 1 := h
  rw [chi4_neg_one] at h1
  norm_num at h1

theorem archParity_chi4 : archParity chi4 = 1 := by
  unfold archParity
  rw [if_neg chi4_not_even]

theorem chi4_ne_one : chi4 ≠ 1 := by
  intro h
  have h3 := chi4_three
  rw [h] at h3
  have hu : IsUnit (3 : ZMod 4) := by decide
  rw [MulChar.one_apply hu] at h3
  norm_num at h3

/-- `χ₋₄` is primitive (its conductor is `4`). -/
theorem chi4_isPrimitive : chi4.IsPrimitive := by
  rw [DirichletCharacter.isPrimitive_def]
  have hd : DirichletCharacter.conductor chi4 ∣ 4 := DirichletCharacter.conductor_dvd_level _
  have hmem : DirichletCharacter.conductor chi4 ∈ Nat.divisors 4 :=
    Nat.mem_divisors.mpr ⟨hd, by norm_num⟩
  have hdiv : Nat.divisors 4 = {1, 2, 4} := by decide
  rw [hdiv] at hmem
  simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
  rcases hmem with h | h | h
  · exact absurd (DirichletCharacter.eq_one_iff_conductor_eq_one.mpr h) chi4_ne_one
  · exfalso
    have hf : chi4.FactorsThrough 2 := h ▸ DirichletCharacter.factorsThrough_conductor chi4
    have heq := hf.eq_changeLevel
    have hu : IsUnit (3 : ZMod 4) := by decide
    have h3 : chi4 ((hu.unit : (ZMod 4)ˣ) : ZMod 4)
        = hf.χ₀ (ZMod.cast ((hu.unit : (ZMod 4)ˣ) : ZMod 4)) := by
      conv_lhs => rw [heq]
      exact DirichletCharacter.changeLevel_eq_cast_of_dvd _ hf.dvd _
    rw [IsUnit.unit_spec, chi4_three] at h3
    have hc : (ZMod.cast (3 : ZMod 4) : ZMod 2) = 1 := by decide
    rw [hc, map_one] at h3
    norm_num at h3
  · exact h

/-- The archimedean weight of `χ₋₄` (odd, modulus 4): `Re ψ(3/4 + ir/2) + log(4/π)`. -/
theorem archWeight_chi4 (r : ℝ) :
    archWeight chi4 r = (Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi) := by
  unfold archWeight
  rw [archParity_chi4]
  norm_num

/-- If `g` is supported in `[−a, a]` with `a < log 3`, the two-sided prime sum of the explicit
formula for `χ₋₄` at shift `0` is empty (`Λ(1) = 0`, `χ₋₄(2) = 0`, and `log n ≥ log 3` for
`n ≥ 3`). -/
theorem primePowerSumSym_chi4_eq_zero {g : ℝ → ℂ} {a : ℝ} (hgs : tsupport g ⊆ Icc (-a) a)
    (ha : a < Real.log 3) : primePowerSumSym chi4 g = 0 := by
  unfold primePowerSumSym
  have hterm : ∀ n : ℕ, ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
      * (chi4 (n : ZMod 4) * g (Real.log n) + conj (chi4 (n : ZMod 4)) * g (-Real.log n)) = 0 := by
    intro n
    by_cases hn : n < 3
    · interval_cases n
      · simp
      · simp
      · have h2 : chi4 ((2 : ℕ) : ZMod 4) = 0 := by
          rw [Nat.cast_ofNat]; exact chi4_two
        rw [h2]; simp
    · have h3 : (3 : ℝ) ≤ n := by exact_mod_cast not_lt.mp hn
      have hlog : Real.log 3 ≤ Real.log n := Real.log_le_log (by norm_num) h3
      have h1 : g (Real.log n) = 0 := by
        apply image_eq_zero_of_notMem_tsupport
        intro hmem
        have := (hgs hmem).2
        linarith
      have h2 : g (-Real.log n) = 0 := by
        apply image_eq_zero_of_notMem_tsupport
        intro hmem
        have := (hgs hmem).1
        linarith
      rw [h1, h2]; simp
  simp only [hterm, tsum_zero]

/-- **S0 for χ₋₄, closed form.**  Under GRH for `L(s,χ₋₄)`, for `g = f ⋆ f̃` with `f ∈ C²_c`
supported in `[−a/2, a/2]`, `a < log 3`, and every shift `L > a`:  `‖S(L) − Γ_L‖ ≤ Re Γ_0`. -/
theorem chi4_norm_primePowerSum_sub_archTerm_le (hGRH : GRHFor chi4)
    {f : ℝ → ℂ} (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) {a : ℝ}
    (ha : a < Real.log 3) (hsupp : tsupport f ⊆ Icc (-(a / 2)) (a / 2)) {L : ℝ} (hL : a < L) :
    ‖primePowerSum chi4 (autocorr f) L - archTerm chi4 (autocorr f) L‖
      ≤ (archTerm chi4 (autocorr f) 0).re := by
  have h := norm_primePowerSum_sub_archTerm_le (by norm_num : 1 < 4) chi4_isPrimitive hGRH hf hfc
    hsupp hL
  have hgs : tsupport (autocorr f) ⊆ Icc (-a) a := by
    rw [autocorr_eq_weilTest]; exact tsupport_weilTest_subset hsupp hsupp
  rwa [primePowerSumSym_chi4_eq_zero hgs ha, sub_zero] at h

end RaceLib
