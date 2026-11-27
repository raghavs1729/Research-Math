/-
RaceLib/Window.lean — the concrete race window `raceWindow`: it is `C²`, compactly supported in
`[−27/50, 27/50]`, `27/25 < log 3`, `∫ raceWindow = 103/100`; and the race theorem for it, whose
only inputs are GRH for `L(s,χ₋₄)` and the externally certified inequality `RaceWindowCert`.
-/
import RaceLib.S1
import RaceLib.PNT

open MeasureTheory Set Filter Topology

noncomputable section

namespace RaceLib

/-! ### Powers of the positive part are `C²` -/

theorem hasDerivAt_posPow (m : ℕ) (t : ℝ) :
    HasDerivAt (fun x : ℝ => (max x 0) ^ (m + 2)) ((m + 2 : ℕ) * (max t 0) ^ (m + 1)) t := by
  rcases lt_trichotomy t 0 with ht | rfl | ht
  · have heq : (fun x : ℝ => (max x 0) ^ (m + 2)) =ᶠ[𝓝 t] fun _ => (0 : ℝ) := by
      filter_upwards [Iio_mem_nhds ht] with x hx
      rw [max_eq_right (le_of_lt hx)]
      exact zero_pow (by omega)
    have h := (hasDerivAt_const t (0 : ℝ)).congr_of_eventuallyEq heq
    exact h.congr_deriv (by rw [max_eq_right ht.le, zero_pow (by omega), mul_zero])
  · have h1 : HasDerivWithinAt (fun x : ℝ => (max x 0) ^ (m + 2)) 0 (Iic 0) 0 := by
      refine (hasDerivWithinAt_const (0 : ℝ) (Iic 0) (0 : ℝ)).congr (fun x hx => ?_) ?_
      · rw [max_eq_right hx]; exact zero_pow (by omega)
      · rw [max_self]; exact zero_pow (by omega)
    have h2 : HasDerivWithinAt (fun x : ℝ => (max x 0) ^ (m + 2)) 0 (Ici 0) 0 := by
      have h := (hasDerivAt_pow (m + 2) (0 : ℝ)).hasDerivWithinAt (s := Ici 0)
      have h' : HasDerivWithinAt (fun x : ℝ => x ^ (m + 2)) 0 (Ici 0) 0 :=
        h.congr_deriv (by
          rw [show m + 2 - 1 = m + 1 from rfl, zero_pow (by omega), mul_zero])
      refine h'.congr (fun x hx => ?_) ?_
      · rw [max_eq_left hx]
      · rw [max_self]
    have h := h1.union h2
    rw [Iic_union_Ici, hasDerivWithinAt_univ] at h
    exact h.congr_deriv (by rw [max_self, zero_pow (by omega), mul_zero])
  · have heq : (fun x : ℝ => (max x 0) ^ (m + 2)) =ᶠ[𝓝 t] fun x => x ^ (m + 2) := by
      filter_upwards [Ioi_mem_nhds ht] with x hx
      rw [max_eq_left (le_of_lt hx)]
    have h := (hasDerivAt_pow (m + 2) t).congr_of_eventuallyEq heq
    exact h.congr_deriv (by rw [max_eq_left ht.le]; rfl)

theorem continuous_posPow (k : ℕ) : Continuous (fun x : ℝ => (max x 0) ^ k) :=
  (continuous_id.max continuous_const).pow k

/-- `x ↦ (max x 0)^(m+3)` is `C²`. -/
theorem contDiff_two_posPow (m : ℕ) : ContDiff ℝ 2 (fun x : ℝ => (max x 0) ^ (m + 3)) := by
  have hd1 : deriv (fun x : ℝ => (max x 0) ^ (m + 3))
      = fun x => ((m + 1 + 2 : ℕ) : ℝ) * (max x 0) ^ (m + 2) :=
    funext fun t => (hasDerivAt_posPow (m + 1) t).deriv
  have hd2 : deriv (fun x : ℝ => ((m + 1 + 2 : ℕ) : ℝ) * (max x 0) ^ (m + 2))
      = fun x => ((m + 1 + 2 : ℕ) : ℝ) * (((m + 2 : ℕ) : ℝ) * (max x 0) ^ (m + 1)) :=
    funext fun t => ((hasDerivAt_posPow m t).const_mul _).deriv
  rw [show (2 : WithTop ℕ∞) = 1 + 1 from rfl, contDiff_succ_iff_deriv]
  refine ⟨fun t => (hasDerivAt_posPow (m + 1) t).differentiableAt, fun h => ?_, ?_⟩
  · exact absurd h (by decide)
  · rw [hd1, contDiff_one_iff_deriv]
    refine ⟨fun t => ((hasDerivAt_posPow m t).const_mul _).differentiableAt, ?_⟩
    rw [hd2]
    exact continuous_const.mul (continuous_const.mul (continuous_posPow _))

/-! ### The smoothstep -/

theorem smoothstep_of_nonpos {t : ℝ} (h : t ≤ 0) : smoothstep t = 0 := by
  unfold smoothstep; rw [if_pos h]

theorem smoothstep_of_mem {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1) :
    smoothstep t = 10 * t ^ 3 - 15 * t ^ 4 + 6 * t ^ 5 := by
  unfold smoothstep
  rcases h0.eq_or_lt with h | h
  · rw [← h]; norm_num
  · rw [if_neg (not_le.mpr h), if_pos h1]

theorem smoothstep_of_one_le {t : ℝ} (h : 1 ≤ t) : smoothstep t = 1 := by
  unfold smoothstep
  rw [if_neg (not_le.mpr (lt_of_lt_of_le one_pos h))]
  split_ifs with h1
  · have : t = 1 := le_antisymm h1 h
    rw [this]; norm_num
  · rfl

theorem smoothstep_eq (t : ℝ) :
    smoothstep t = (10 * (max t 0) ^ 3 - 15 * (max t 0) ^ 4 + 6 * (max t 0) ^ 5)
      - (10 * (max (t - 1) 0) ^ 3 + 15 * (max (t - 1) 0) ^ 4 + 6 * (max (t - 1) 0) ^ 5) := by
  rcases le_or_gt t 0 with h0 | h0
  · rw [smoothstep_of_nonpos h0, max_eq_right h0, max_eq_right (by linarith)]
    norm_num
  rcases le_or_gt t 1 with h1 | h1
  · rw [smoothstep_of_mem h0.le h1, max_eq_left h0.le, max_eq_right (by linarith)]
    norm_num
  · rw [smoothstep_of_one_le h1.le, max_eq_left h0.le, max_eq_left (by linarith)]
    ring

theorem contDiff_smoothstep : ContDiff ℝ 2 smoothstep := by
  have h3 := contDiff_two_posPow 0
  have h4 := contDiff_two_posPow 1
  have h5 := contDiff_two_posPow 2
  have hs : ContDiff ℝ 2 (fun t : ℝ => t - 1) := contDiff_id.sub contDiff_const
  have key : ContDiff ℝ 2 (fun t : ℝ =>
      (10 * (max t 0) ^ 3 - 15 * (max t 0) ^ 4 + 6 * (max t 0) ^ 5)
      - (10 * (max (t - 1) 0) ^ 3 + 15 * (max (t - 1) 0) ^ 4 + 6 * (max (t - 1) 0) ^ 5)) :=
    (((contDiff_const.mul h3).sub (contDiff_const.mul h4)).add (contDiff_const.mul h5)).sub
      (((contDiff_const.mul (h3.comp hs)).add (contDiff_const.mul (h4.comp hs))).add
        (contDiff_const.mul (h5.comp hs)))
  have : smoothstep = fun t : ℝ =>
      (10 * (max t 0) ^ 3 - 15 * (max t 0) ^ 4 + 6 * (max t 0) ^ 5)
      - (10 * (max (t - 1) 0) ^ 3 + 15 * (max (t - 1) 0) ^ 4 + 6 * (max (t - 1) 0) ^ 5) :=
    funext smoothstep_eq
  rw [this]
  exact key

/-! ### The window -/

theorem contDiff_raceWindow : ContDiff ℝ 2 raceWindow := by
  have h1 : ContDiff ℝ 2 (fun u : ℝ => 20 * (27 / 50 - u)) :=
    contDiff_const.mul (contDiff_const.sub contDiff_id)
  have h2 : ContDiff ℝ 2 (fun u : ℝ => 20 * (27 / 50 + u)) :=
    contDiff_const.mul (contDiff_const.add contDiff_id)
  exact (contDiff_smoothstep.comp h1).mul (contDiff_smoothstep.comp h2)

theorem raceWindow_of_ge {u : ℝ} (h : 27 / 50 ≤ u) : raceWindow u = 0 := by
  unfold raceWindow
  rw [smoothstep_of_nonpos (by linarith), zero_mul]

theorem raceWindow_of_le {u : ℝ} (h : u ≤ -(27 / 50)) : raceWindow u = 0 := by
  unfold raceWindow
  rw [smoothstep_of_nonpos (t := 20 * (27 / 50 + u)) (by linarith), mul_zero]

theorem raceWindow_of_abs_le {u : ℝ} (h : |u| ≤ 49 / 100) : raceWindow u = 1 := by
  have h' := abs_le.mp h
  unfold raceWindow
  rw [smoothstep_of_one_le (by linarith [h'.2]), smoothstep_of_one_le (by linarith [h'.1]), mul_one]

theorem raceWindow_right {u : ℝ} (h1 : 49 / 100 ≤ u) (h2 : u ≤ 27 / 50) :
    raceWindow u = 10 * (20 * (27 / 50 - u)) ^ 3 - 15 * (20 * (27 / 50 - u)) ^ 4
      + 6 * (20 * (27 / 50 - u)) ^ 5 := by
  unfold raceWindow
  rw [smoothstep_of_mem (by linarith) (by linarith), smoothstep_of_one_le (by linarith), mul_one]

theorem raceWindow_left {u : ℝ} (h1 : -(27 / 50) ≤ u) (h2 : u ≤ -(49 / 100)) :
    raceWindow u = 10 * (20 * (27 / 50 + u)) ^ 3 - 15 * (20 * (27 / 50 + u)) ^ 4
      + 6 * (20 * (27 / 50 + u)) ^ 5 := by
  unfold raceWindow
  rw [smoothstep_of_one_le (by linarith), smoothstep_of_mem (t := 20 * (27 / 50 + u))
    (by linarith) (by linarith), one_mul]

theorem raceWindow_neg (u : ℝ) : raceWindow (-u) = raceWindow u := by
  unfold raceWindow
  rw [sub_neg_eq_add, ← sub_eq_add_neg, mul_comm]

theorem raceWindow_piecewise (u : ℝ) :
    raceWindow (-u) = raceWindow u
      ∧ (|u| ≤ 49 / 100 → raceWindow u = 1)
      ∧ (49 / 100 ≤ u → u ≤ 27 / 50 →
          raceWindow u = 10 * (20 * (27 / 50 - u)) ^ 3 - 15 * (20 * (27 / 50 - u)) ^ 4
            + 6 * (20 * (27 / 50 - u)) ^ 5)
      ∧ (27 / 50 ≤ |u| → raceWindow u = 0) := by
  refine ⟨raceWindow_neg u, raceWindow_of_abs_le, raceWindow_right, fun h => ?_⟩
  rcases le_abs.mp h with h | h
  · exact raceWindow_of_ge h
  · exact raceWindow_of_le (by linarith)

theorem support_raceWindow_subset : Function.support raceWindow ⊆ Ioo (-(27 / 50)) (27 / 50) := by
  intro u hu
  rw [Function.mem_support] at hu
  constructor
  · by_contra h
    exact hu (raceWindow_of_le (not_lt.mp h))
  · by_contra h
    exact hu (raceWindow_of_ge (not_lt.mp h))

theorem tsupport_raceWindow_subset : tsupport raceWindow ⊆ Icc (-(27 / 50)) (27 / 50) :=
  closure_minimal (support_raceWindow_subset.trans Ioo_subset_Icc_self) isClosed_Icc

theorem hasCompactSupport_raceWindow : HasCompactSupport raceWindow :=
  IsCompact.of_isClosed_subset isCompact_Icc (isClosed_tsupport _) tsupport_raceWindow_subset

/-- `27/25 < log 3`, i.e. `e^{27/25} < 3`. -/
theorem log_three_gt : (27 / 25 : ℝ) < Real.log 3 := by
  rw [Real.lt_log_iff_exp_lt (by norm_num)]
  have h1 : Real.exp (27 / 25 : ℝ) = Real.exp 1 * Real.exp (2 / 25) := by
    rw [← Real.exp_add]; norm_num
  have h2 : Real.exp (2 / 25 : ℝ) ≤ 25 / 23 := by
    have h := Real.add_one_le_exp (-(2 / 25 : ℝ))
    rw [Real.exp_neg] at h
    have hpos : 0 < Real.exp (2 / 25 : ℝ) := Real.exp_pos _
    have h3 : (23 / 25 : ℝ) ≤ (Real.exp (2 / 25))⁻¹ := by linarith
    rw [le_inv_comm₀ (by norm_num) hpos] at h3
    linarith [h3, show ((23 : ℝ) / 25)⁻¹ = 25 / 23 by norm_num]
  have h3 := Real.exp_one_lt_d9
  have hpos : 0 < Real.exp (1 : ℝ) := Real.exp_pos _
  rw [h1]
  nlinarith

/-! ### The integral of the window -/

theorem integral_quintic :
    ∫ x in (0 : ℝ)..1, (10 * x ^ 3 - 15 * x ^ 4 + 6 * x ^ 5) = 1 / 2 := by
  have hderiv : ∀ x ∈ uIcc (0 : ℝ) 1,
      HasDerivAt (fun x : ℝ => 5 / 2 * x ^ 4 - 3 * x ^ 5 + x ^ 6)
        (10 * x ^ 3 - 15 * x ^ 4 + 6 * x ^ 5) x := by
    intro x _
    have h := (((hasDerivAt_pow 4 x).const_mul (5 / 2 : ℝ)).sub
      ((hasDerivAt_pow 5 x).const_mul (3 : ℝ))).add (hasDerivAt_pow 6 x)
    exact h.congr_deriv (by
      simp only [Nat.cast_ofNat, show (4 - 1 : ℕ) = 3 from rfl, show (5 - 1 : ℕ) = 4 from rfl,
        show (6 - 1 : ℕ) = 5 from rfl]
      ring)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv
    (Continuous.intervalIntegrable (by fun_prop) _ _)]
  norm_num

theorem integral_raceWindow : ∫ x, raceWindow x = 103 / 100 := by
  have hcont : Continuous raceWindow := contDiff_raceWindow.continuous
  have hsupp : Function.support raceWindow ⊆ Ioc (-(27 / 50)) (27 / 50) :=
    support_raceWindow_subset.trans Ioo_subset_Ioc_self
  rw [← intervalIntegral.integral_eq_integral_of_support_subset hsupp,
    ← intervalIntegral.integral_add_adjacent_intervals (b := -(49 / 100))
      (hcont.intervalIntegrable _ _) (hcont.intervalIntegrable _ _),
    ← intervalIntegral.integral_add_adjacent_intervals (a := -(49 / 100)) (b := 49 / 100)
      (hcont.intervalIntegrable _ _) (hcont.intervalIntegrable _ _)]
  -- left ramp
  have hL : ∫ x in (-(27 / 50 : ℝ))..(-(49 / 100)), raceWindow x = 1 / 40 := by
    have hcongr : ∫ x in (-(27 / 50 : ℝ))..(-(49 / 100)), raceWindow x
        = ∫ x in (-(27 / 50 : ℝ))..(-(49 / 100)),
            (fun t : ℝ => 10 * t ^ 3 - 15 * t ^ 4 + 6 * t ^ 5) (54 / 5 + 20 * x) := by
      refine intervalIntegral.integral_congr fun x hx => ?_
      rw [uIcc_of_le (by norm_num)] at hx
      simp only
      rw [raceWindow_left hx.1 hx.2, show (20 : ℝ) * (27 / 50 + x) = 54 / 5 + 20 * x by ring]
    have h := intervalIntegral.integral_comp_add_mul (a := -(27 / 50 : ℝ)) (b := -(49 / 100))
      (fun t : ℝ => 10 * t ^ 3 - 15 * t ^ 4 + 6 * t ^ 5) (by norm_num : (20 : ℝ) ≠ 0) (54 / 5)
    have e1 : (54 / 5 + 20 * -(27 / 50) : ℝ) = 0 := by norm_num
    have e2 : (54 / 5 + 20 * -(49 / 100) : ℝ) = 1 := by norm_num
    rw [e1, e2, integral_quintic] at h
    rw [hcongr]
    refine h.trans ?_
    norm_num
  -- plateau
  have hM : ∫ x in (-(49 / 100 : ℝ))..(49 / 100), raceWindow x = 49 / 50 := by
    have hcongr : ∫ x in (-(49 / 100 : ℝ))..(49 / 100), raceWindow x
        = ∫ _x in (-(49 / 100 : ℝ))..(49 / 100), (1 : ℝ) := by
      refine intervalIntegral.integral_congr fun x hx => ?_
      rw [uIcc_of_le (by norm_num)] at hx
      exact raceWindow_of_abs_le (abs_le.mpr hx)
    rw [hcongr, intervalIntegral.integral_const]
    norm_num
  -- right ramp
  have hR : ∫ x in (49 / 100 : ℝ)..(27 / 50), raceWindow x = 1 / 40 := by
    have hcongr : ∫ x in (49 / 100 : ℝ)..(27 / 50), raceWindow x
        = ∫ x in (49 / 100 : ℝ)..(27 / 50),
            (fun t : ℝ => 10 * t ^ 3 - 15 * t ^ 4 + 6 * t ^ 5) (54 / 5 - 20 * x) := by
      refine intervalIntegral.integral_congr fun x hx => ?_
      rw [uIcc_of_le (by norm_num)] at hx
      simp only
      rw [raceWindow_right hx.1 hx.2, show (20 : ℝ) * (27 / 50 - x) = 54 / 5 - 20 * x by ring]
    have h := intervalIntegral.integral_comp_sub_mul (a := (49 / 100 : ℝ)) (b := 27 / 50)
      (fun t : ℝ => 10 * t ^ 3 - 15 * t ^ 4 + 6 * t ^ 5) (by norm_num : (20 : ℝ) ≠ 0) (54 / 5)
    have e1 : (54 / 5 - 20 * (27 / 50) : ℝ) = 0 := by norm_num
    have e2 : (54 / 5 - 20 * (49 / 100) : ℝ) = 1 := by norm_num
    rw [e1, e2, integral_quintic] at h
    rw [hcongr]
    refine h.trans ?_
    norm_num
  rw [hL, hM, hR]
  norm_num

/-! ### The race for the concrete window -/

theorem tsupport_raceWindow_subset' :
    tsupport raceWindow ⊆ Icc (-((27 / 25 : ℝ) / 2)) ((27 / 25) / 2) := by
  have h : ((27 / 25 : ℝ) / 2) = 27 / 50 := by norm_num
  rw [h]
  exact tsupport_raceWindow_subset

theorem raceWindowCert_iff :
    RaceWindowCert ↔ chi4Gamma0 raceWindow < (1 / 2) * (∫ x, raceWindow x) ^ 2 - 7369 / 20000 := by
  unfold RaceWindowCert
  rw [integral_raceWindow]
  norm_num

/-- **The race for the concrete window.**  Inputs: GRH for `L(s,χ₋₄)` and the externally certified
inequality `RaceWindowCert`. -/
theorem raceWindow_race (hGRH : GRHFor chi4) (hcert : RaceWindowCert) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum raceWindow L < -(7369 / 20000) :=
  chi4RaceSum_lt hGRH pntTheta contDiff_raceWindow hasCompactSupport_raceWindow log_three_gt
    tsupport_raceWindow_subset' (raceWindowCert_iff.mp hcert)

/-- **Two-sided form for the concrete window** (GRH only, no numerical input). -/
theorem raceWindow_two_sided (hGRH : GRHFor chi4) {ε : ℝ} (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀,
      |chi4RaceSum raceWindow L + 10609 / 20000| ≤ chi4Gamma0 raceWindow + ε := by
  obtain ⟨L₀, h⟩ := abs_chi4RaceSum_add_le hGRH pntTheta contDiff_raceWindow
    hasCompactSupport_raceWindow log_three_gt tsupport_raceWindow_subset' hε
  refine ⟨L₀, fun L hL => ?_⟩
  have := h L hL
  rw [integral_raceWindow] at this
  norm_num at this ⊢
  exact this

end RaceLib
