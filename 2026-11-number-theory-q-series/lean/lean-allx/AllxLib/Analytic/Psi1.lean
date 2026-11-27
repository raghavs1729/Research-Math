/-
AllxLib/Analytic/Psi1.lean — the bound for `ψ₁(t, χ₋₄) = Σ_{n ≤ t} Λ(n) χ₋₄(n) (t − n)` under GRH.

The explicit formula `EF_ext` for the `C²` test function `Kt t` reads

    Σ_ρ m_ρ K̂(γ_ρ) = −(Σ_n Λ(n)χ(n) H(t − n) + Σ_n Λ(n)χ(n) (t − 1/n)(1 − (1 − 1/n)³)/n) + arch,

with `|Σ_ρ| ≤ (27/100)(t+1)^{3/2} + O(t)` (`ZeroSum.lean`, `Numerics.lean`), `arch ≤ log t + 3`
(`ArchK.lean`), the second prime sum is `≤ (7/2) t` in absolute value, and
`Σ_n Λ(n)χ(n) H(t − n)` differs from `ψ₁(t,χ)` by at most `(1/3) log(t+1)`.  Hence

    ψ₁(t, χ₋₄) ≤ (27/100)(t+1)^{3/2} + 4t + log t + 4 + (1/3) log(t+1)      (t ≥ 2).
-/
import AllxLib.Analytic.ZeroSum
import AllxLib.Analytic.ArchK
import AllxLib.Analytic.PrimePowers

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE ChebLib

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace Allx

/-! ### The prime side -/

/-- The weight of the dual prime sum: `(t − 1/n)(1 − (1 − 1/n)³)/n`. -/
def plw (t : ℝ) (n : ℕ) : ℝ := (t - (n : ℝ)⁻¹) * (1 - (1 - (n : ℝ)⁻¹) ^ 3) / n

theorem plw_bounds {t : ℝ} (ht : 2 ≤ t) {n : ℕ} (hn : 0 < n) :
    0 ≤ plw t n ∧ plw t n ≤ 3 * t / (n : ℝ) ^ 2 := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hx0 : 0 < (n : ℝ)⁻¹ := inv_pos.mpr hn'
  have hx1 : (n : ℝ)⁻¹ ≤ 1 := inv_le_one_of_one_le₀ h1
  have e : plw t n = (t - (n : ℝ)⁻¹) * (1 - (1 - (n : ℝ)⁻¹) ^ 3) * (n : ℝ)⁻¹ := by
    unfold plw; rw [div_eq_mul_inv]
  have e2 : 3 * t / (n : ℝ) ^ 2 = t * (3 * (n : ℝ)⁻¹) * (n : ℝ)⁻¹ := by
    field_simp
  rw [e, e2]
  generalize (n : ℝ)⁻¹ = x at hx0 hx1
  have hA : 0 ≤ t - x := by linarith
  have hB : 0 ≤ 1 - (1 - x) ^ 3 := by nlinarith [mul_pos hx0 hx0]
  have hB' : 1 - (1 - x) ^ 3 ≤ 3 * x := by nlinarith [mul_pos hx0 hx0]
  constructor
  · positivity
  · apply mul_le_mul_of_nonneg_right _ hx0.le
    exact mul_le_mul (by linarith) hB' hB (by linarith)

theorem Hhat_le_div_sq {t : ℝ} (ht : 2 ≤ t) {n : ℕ} (hn : 0 < n) :
    Hhat (t - n) ≤ (t + 1) ^ 3 / (n : ℝ) ^ 2 := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  rcases le_or_gt (n : ℝ) (t + 1) with h | h
  · have h1 := (Hhat_sub_max (t - n)).2
    have h2 : max (t - n) 0 ≤ t := max_le (by linarith) (by linarith)
    rw [le_div_iff₀ (by positivity)]
    have h3 : Hhat (t - n) ≤ t + 1 := by linarith
    have h4 : 0 ≤ Hhat (t - n) := Hhat_nonneg _
    have h5 : (n : ℝ) ^ 2 ≤ (t + 1) ^ 2 := pow_le_pow_left₀ hn'.le h 2
    calc Hhat (t - n) * (n : ℝ) ^ 2 ≤ (t + 1) * (t + 1) ^ 2 :=
          mul_le_mul h3 h5 (by positivity) (by linarith)
      _ = (t + 1) ^ 3 := by ring
  · rw [Hhat_of_le_neg_one (by linarith)]
    positivity

/-- The `n`-th term of the prime side of the explicit formula for `Kt t`. -/
theorem primeTerm_Kt {t : ℝ} (ht : 2 ≤ t) (n : ℕ) :
    ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
        * (coeff chi4 n * ((Kt t (Real.log n) : ℝ) : ℂ)
          + (starRingEnd ℂ) (coeff chi4 n) * ((Kt t (-Real.log n) : ℝ) : ℂ))
      = ((vM n * X4 n * Hhat (t - n) + vM n * X4 n * plw t n : ℝ) : ℂ) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [vM]
  · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
    have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
    rw [coeff_chi4, Complex.conj_ofReal, Kt_log t hn, Kt_neg_log ht hn]
    unfold plw vM
    generalize (t - (n : ℝ)⁻¹) * (1 - (1 - (n : ℝ)⁻¹) ^ 3) = P
    generalize Hhat (t - n) = H
    norm_cast
    generalize Real.sqrt n = s at hs hss
    rw [hss]
    field_simp

theorem summable_primeSide_Kt {t : ℝ} (ht : 2 ≤ t) :
    Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖((Kt t (Real.log n) : ℝ) : ℂ)‖ + ‖((Kt t (-Real.log n) : ℝ) : ℂ)‖)) := by
  refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_)
    (summable_vM_div_sq.mul_left ((t + 1) ^ 3 + 3 * t))
  · exact mul_nonneg (div_nonneg (vM_nonneg n) (Real.sqrt_nonneg _)) (by positivity)
  · rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp
    · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
      have hs : 0 < Real.sqrt n := Real.sqrt_pos.mpr hn'
      have hss : (n : ℝ) = Real.sqrt n ^ 2 := (Real.sq_sqrt hn'.le).symm
      have hΛ := vM_nonneg n
      rw [Complex.norm_real, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs,
        abs_of_nonneg (Kt_nonneg ht _), abs_of_nonneg (Kt_nonneg ht _), Kt_log t hn,
        Kt_neg_log ht hn]
      obtain ⟨hp0, hp1⟩ := plw_bounds ht hn
      have hH := Hhat_le_div_sq ht hn
      have e : ArithmeticFunction.vonMangoldt n / Real.sqrt n
            * (Real.sqrt n * Hhat (t - n)
              + (Real.sqrt n)⁻¹ * ((t - (n : ℝ)⁻¹) * (1 - (1 - (n : ℝ)⁻¹) ^ 3)))
          = vM n * Hhat (t - n) + vM n * plw t n := by
        unfold plw vM
        generalize (t - (n : ℝ)⁻¹) * (1 - (1 - (n : ℝ)⁻¹) ^ 3) = P
        generalize Hhat (t - n) = H
        generalize Real.sqrt n = s at hs hss
        rw [hss]
        field_simp
      rw [e]
      calc vM n * Hhat (t - n) + vM n * plw t n
          ≤ vM n * ((t + 1) ^ 3 / (n : ℝ) ^ 2) + vM n * (3 * t / (n : ℝ) ^ 2) :=
            add_le_add (mul_le_mul_of_nonneg_left hH hΛ) (mul_le_mul_of_nonneg_left hp1 hΛ)
        _ = ((t + 1) ^ 3 + 3 * t) * (vM n / (n : ℝ) ^ 2) := by ring

theorem summable_H {t : ℝ} (ht : 2 ≤ t) : Summable (fun n : ℕ => vM n * X4 n * Hhat (t - n)) := by
  refine Summable.of_norm_bounded (summable_vM_div_sq.mul_left ((t + 1) ^ 3)) fun n => ?_
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [vM]
  · rw [Real.norm_eq_abs, abs_mul, abs_mul, abs_of_nonneg (vM_nonneg n),
      abs_of_nonneg (Hhat_nonneg _)]
    calc vM n * |X4 n| * Hhat (t - n) ≤ vM n * 1 * ((t + 1) ^ 3 / (n : ℝ) ^ 2) :=
          mul_le_mul (mul_le_mul_of_nonneg_left (abs_X4_le_one n) (vM_nonneg n))
            (Hhat_le_div_sq ht hn) (Hhat_nonneg _) (by have := vM_nonneg n; positivity)
      _ = (t + 1) ^ 3 * (vM n / (n : ℝ) ^ 2) := by ring

theorem abs_P_le {t : ℝ} (ht : 2 ≤ t) (n : ℕ) :
    |vM n * X4 n * plw t n| ≤ 3 * t * (vM n / (n : ℝ) ^ 2) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [vM]
  · obtain ⟨hp0, hp1⟩ := plw_bounds ht hn
    rw [abs_mul, abs_mul, abs_of_nonneg (vM_nonneg n), abs_of_nonneg hp0]
    calc vM n * |X4 n| * plw t n ≤ vM n * 1 * (3 * t / (n : ℝ) ^ 2) :=
          mul_le_mul (mul_le_mul_of_nonneg_left (abs_X4_le_one n) (vM_nonneg n)) hp1 hp0
            (by have := vM_nonneg n; positivity)
      _ = 3 * t * (vM n / (n : ℝ) ^ 2) := by ring

theorem summable_P {t : ℝ} (ht : 2 ≤ t) : Summable (fun n : ℕ => vM n * X4 n * plw t n) :=
  Summable.of_norm_bounded (summable_vM_div_sq.mul_left (3 * t)) fun n => by
    rw [Real.norm_eq_abs]; exact abs_P_le ht n

theorem abs_tsum_P_le {t : ℝ} (ht : 2 ≤ t) : |∑' n : ℕ, vM n * X4 n * plw t n| ≤ 7 / 2 * t := by
  have h1 : |∑' n : ℕ, vM n * X4 n * plw t n| ≤ ∑' n : ℕ, 3 * t * (vM n / (n : ℝ) ^ 2) := by
    rw [← Real.norm_eq_abs]
    refine (norm_tsum_le_tsum_norm (summable_P ht).norm).trans ?_
    exact (summable_P ht).norm.tsum_le_tsum (fun n => by rw [Real.norm_eq_abs]; exact abs_P_le ht n)
      (summable_vM_div_sq.mul_left _)
  rw [tsum_mul_left] at h1
  have h2 : ∑' n : ℕ, vM n / (n : ℝ) ^ 2 ≤ 7 / 6 := tsum_vonMangoldt_div_sq_le
  nlinarith

/-! ### `Σ Λ χ H(t − n)` against `ψ₁(t,χ)` -/

theorem psi1chi_le_tsum_H {t : ℝ} (ht : 2 ≤ t) :
    psi1chi t - 1 / 3 * Real.log (t + 1) ≤ ∑' n : ℕ, vM n * X4 n * Hhat (t - n) := by
  have ht0 : 0 ≤ t := by linarith
  set F : ℕ := ⌊t⌋₊ with hF
  have hFt : (F : ℝ) ≤ t := Nat.floor_le ht0
  have htF : t < (F : ℝ) + 1 := Nat.lt_floor_add_one t
  have hF2 : 2 ≤ F := Nat.le_floor (by exact_mod_cast ht)
  -- `ψ₁` as a series with the weight `max (t − n) 0`
  have hpsi : psi1chi t = ∑' n : ℕ, vM n * X4 n * max (t - n) 0 := by
    rw [tsum_eq_sum (s := Finset.Ioc 0 F)]
    · unfold psi1chi
      refine Finset.sum_congr rfl fun n hn => ?_
      have hnF : (n : ℝ) ≤ F := by exact_mod_cast (Finset.mem_Ioc.mp hn).2
      rw [max_eq_left (by linarith)]
      rfl
    · intro n hn
      rcases Nat.eq_zero_or_pos n with rfl | hn0
      · simp [vM]
      · have : F < n := by
          by_contra h
          exact hn (Finset.mem_Ioc.mpr ⟨hn0, not_lt.mp h⟩)
        have hnF : (F : ℝ) + 1 ≤ n := by exact_mod_cast this
        rw [max_eq_right (by linarith)]
        ring
  have hsM : Summable (fun n : ℕ => vM n * X4 n * max (t - n) 0) := by
    refine summable_of_ne_finset_zero (s := Finset.Ioc 0 F) fun n hn => ?_
    rcases Nat.eq_zero_or_pos n with rfl | hn0
    · simp [vM]
    · have : F < n := by
        by_contra h
        exact hn (Finset.mem_Ioc.mpr ⟨hn0, not_lt.mp h⟩)
      have hnF : (F : ℝ) + 1 ≤ n := by exact_mod_cast this
      rw [max_eq_right (by linarith)]
      ring
  -- the difference is supported on `{F, F+1}`
  have hd0 : ∀ n : ℕ, n ∉ ({F, F + 1} : Finset ℕ) → Hhat (t - n) - max (t - n) 0 = 0 := by
    intro n hn
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hn
    rcases Nat.lt_or_ge n F with h | h
    · have : (n : ℝ) + 1 ≤ F := by exact_mod_cast h
      rw [Hhat_of_one_le (by linarith), max_eq_left (by linarith)]
      ring
    · have h2 : F + 2 ≤ n := by omega
      have : (F : ℝ) + 2 ≤ n := by exact_mod_cast h2
      rw [Hhat_of_le_neg_one (by linarith), max_eq_right (by linarith)]
      ring
  have hdiff : ∑' n : ℕ, vM n * X4 n * Hhat (t - n) - ∑' n : ℕ, vM n * X4 n * max (t - n) 0
      = ∑ n ∈ ({F, F + 1} : Finset ℕ), vM n * X4 n * (Hhat (t - n) - max (t - n) 0) := by
    rw [← (summable_H ht).tsum_sub hsM, tsum_eq_sum (s := ({F, F + 1} : Finset ℕ))]
    · refine Finset.sum_congr rfl fun n _ => ?_
      ring
    · intro n hn
      rw [← mul_sub, hd0 n hn, mul_zero]
  have hterm : ∀ n : ℕ, 0 < n → (n : ℝ) ≤ t + 1 →
      -(1 / 6 * Real.log (t + 1)) ≤ vM n * X4 n * (Hhat (t - n) - max (t - n) 0) := by
    intro n hn hnt
    have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    obtain ⟨d0, d1⟩ := Hhat_sub_max (t - n)
    have hΛ0 := vM_nonneg n
    have hΛ1 : vM n ≤ Real.log (t + 1) :=
      (ArithmeticFunction.vonMangoldt_le_log).trans (Real.log_le_log hn' hnt)
    have hX := abs_X4_le_one n
    have hX2 : -1 ≤ X4 n := by linarith [neg_abs_le (X4 n)]
    have hL : 0 ≤ Real.log (t + 1) := Real.log_nonneg (by linarith)
    have h1 : -(vM n * (Hhat (t - n) - max (t - n) 0))
        ≤ vM n * X4 n * (Hhat (t - n) - max (t - n) 0) := by
      have := mul_nonneg hΛ0 d0
      nlinarith
    have h2 : vM n * (Hhat (t - n) - max (t - n) 0) ≤ Real.log (t + 1) * (1 / 6) :=
      mul_le_mul hΛ1 d1 d0 hL
    linarith
  have hsum2 : -(1 / 3 * Real.log (t + 1))
      ≤ ∑ n ∈ ({F, F + 1} : Finset ℕ), vM n * X4 n * (Hhat (t - n) - max (t - n) 0) := by
    rw [Finset.sum_pair (by omega)]
    have a1 := hterm F (by omega) (by linarith)
    have a2 := hterm (F + 1) (by omega) (by push_cast; linarith)
    linarith
  rw [hpsi]
  linarith

/-! ### The zero side -/

/-- The bound for `‖K̂(r)‖`. -/
def zbound (t r : ℝ) : ℝ :=
  (t + 1) ^ (3 / 2 : ℝ) * (1 / (9 / 4 + r ^ 2) + 5 / 12 * majF r) + (3 / 7 * t + 1 / 5) * majF r

theorem tsum_zbound_le (hGRH : GRHFor chi4) {t : ℝ} (ht : 2 ≤ t) :
    Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * zbound t (ρ : ℂ).im) ∧
      ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * zbound t (ρ : ℂ).im
        ≤ 27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 81 / 125 * (3 / 7 * t + 1 / 5) := by
  obtain ⟨hs1, hle1⟩ := sumOne_le hGRH
  obtain ⟨hsM, heM⟩ := sumMaj_eq hGRH
  set T : ℝ := (t + 1) ^ (3 / 2 : ℝ) with hT
  have hT0 : 0 ≤ T := Real.rpow_nonneg (by linarith) _
  have e : ∀ ρ : Z4.carrier, (Z4.mult ρ : ℝ) * zbound t (ρ : ℂ).im
      = T / 3 * ((Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2)))
        + (5 / 12 * T + (3 / 7 * t + 1 / 5)) * ((Z4.mult ρ : ℝ) * majF (ρ : ℂ).im) := by
    intro ρ
    unfold zbound
    rw [← hT]
    ring
  have hsum : Summable (fun ρ : Z4.carrier => (Z4.mult ρ : ℝ) * zbound t (ρ : ℂ).im) :=
    ((hs1.mul_left _).add (hsM.mul_left _)).congr fun ρ => (e ρ).symm
  refine ⟨hsum, ?_⟩
  rw [tsum_congr e, (hs1.mul_left _).tsum_add (hsM.mul_left _), tsum_mul_left, tsum_mul_left, heM]
  set L1 : ℝ := 2 - Real.eulerMascheroniConstant - Real.log Real.pi + 2 * Aser with hL1
  set L2 : ℝ := 7 / 2 - Real.eulerMascheroniConstant - 3 * Real.log 2 - Real.log Real.pi
    + 5 * Aser - 3 * Bser with hL2
  set S1 : ℝ := ∑' ρ : Z4.carrier, (Z4.mult ρ : ℝ) * (3 / (9 / 4 + (ρ : ℂ).im ^ 2)) with hS1
  have hS10 : 0 ≤ S1 := tsum_nonneg fun ρ => mul_nonneg (Nat.cast_nonneg _) (by positivity)
  have hnum : L1 / 3 + 5 / 12 * L2 ≤ 27 / 100 := by
    have := Sbar_numeric
    rw [hL1, hL2]
    linarith
  have hL2nn : 0 ≤ L2 := by
    rw [← heM]
    exact tsum_nonneg fun ρ => mul_nonneg (Nat.cast_nonneg _) (majF_nonneg _)
  have hL2le : L2 ≤ 81 / 125 := by linarith
  have h1 : T / 3 * S1 ≤ T / 3 * L1 := mul_le_mul_of_nonneg_left hle1 (by positivity)
  have h2 : (T / 3 * L1 + 5 / 12 * T * L2) ≤ 27 / 100 * T := by
    have := mul_le_mul_of_nonneg_left hnum hT0
    linarith
  have h3 : (3 / 7 * t + 1 / 5) * L2 ≤ (3 / 7 * t + 1 / 5) * (81 / 125) :=
    mul_le_mul_of_nonneg_left hL2le (by linarith)
  linarith

/-! ### The bound -/

/-- **The bound for `ψ₁(t, χ₋₄)` under GRH for `L(s, χ₋₄)`.** -/
theorem psi1chi_le (hGRH : GRHFor chi4) {t : ℝ} (ht : 2 ≤ t) :
    psi1chi t ≤ 27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 4 * t + Real.log t + 4
      + 1 / 3 * Real.log (t + 1) := by
  have ht1 : 1 < t := by linarith
  obtain ⟨hsum, hEF⟩ := EF_ext (χ := chi4) one_lt_four RaceLib.chi4_isPrimitive hGRH
    (k := fun u => ((Kt t u : ℝ) : ℂ)) (k1 := fun u => ((Kt1 t u : ℝ) : ℂ))
    (k2 := fun u => ((Kt2 t u : ℝ) : ℂ)) (fun u => (hasDerivAt_Kt t u).ofReal_comp)
    (fun u => (hasDerivAt_Kt1 t u).ofReal_comp)
    (Complex.continuous_ofReal.comp (continuous_Kt2 t)) (integrable_Kt ht1).ofReal
    (integrable_Kt1 ht1).ofReal (integrable_Kt2 ht1).ofReal (summable_primeSide_Kt ht)
  -- the zero side
  obtain ⟨hzs, hzle⟩ := tsum_zbound_le hGRH ht
  have hzterm : ∀ ρ : Z4.carrier,
      ‖(Z4.mult ρ : ℂ) * paperFT (fun u => ((Kt t u : ℝ) : ℂ)) (gammaOf ρ)‖
        ≤ (Z4.mult ρ : ℝ) * zbound t (ρ : ℂ).im := by
    intro ρ
    rw [Z4_gammaOf hGRH ρ, norm_mul, Complex.norm_natCast]
    exact mul_le_mul_of_nonneg_left (norm_integral_Kt_mul_cexp_le ht _) (Nat.cast_nonneg _)
  have hZ : ‖∑' ρ : Z4.carrier,
      (Z4.mult ρ : ℂ) * paperFT (fun u => ((Kt t u : ℝ) : ℂ)) (gammaOf ρ)‖
      ≤ 27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 81 / 125 * (3 / 7 * t + 1 / 5) := by
    refine (norm_tsum_le_tsum_norm hsum.norm).trans ?_
    exact (hsum.norm.tsum_le_tsum hzterm hzs).trans hzle
  -- the archimedean side
  obtain ⟨a, ha, harch⟩ := arch_Kt_le ht
  have harch' : (1 / (2 * Real.pi) : ℂ) * ∫ r : ℝ, paperFT (fun u => ((Kt t u : ℝ) : ℂ)) r
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
  -- the explicit formula, in real terms
  have hH := summable_H ht
  have hP := summable_P ht
  have hEF' : (∑' ρ : Z4.carrier,
      (Z4.mult ρ : ℂ) * paperFT (fun u => ((Kt t u : ℝ) : ℂ)) (gammaOf ρ))
      = ((-(∑' n : ℕ, vM n * X4 n * Hhat (t - n) + ∑' n : ℕ, vM n * X4 n * plw t n) + a : ℝ) : ℂ) := by
    rw [hEF]
    unfold literatureRHSchi
    rw [harch', tsum_congr (primeTerm_Kt ht), ← Complex.ofReal_tsum, hH.tsum_add hP]
    push_cast
    ring
  have hre : -(∑' n : ℕ, vM n * X4 n * Hhat (t - n) + ∑' n : ℕ, vM n * X4 n * plw t n) + a
      ≥ -(27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 81 / 125 * (3 / 7 * t + 1 / 5)) := by
    rw [hEF', Complex.norm_real, Real.norm_eq_abs] at hZ
    linarith [neg_abs_le (-(∑' n : ℕ, vM n * X4 n * Hhat (t - n)
      + ∑' n : ℕ, vM n * X4 n * plw t n) + a)]
  have hPb := abs_tsum_P_le ht
  have hPb' := neg_abs_le (∑' n : ℕ, vM n * X4 n * plw t n)
  have hpsi := psi1chi_le_tsum_H ht
  linarith

end Allx
