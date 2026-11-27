/-
RamanujanLib/Main.lean — Ramanujan's inequality `π(x)² < (e x / log x) π(x / e)` for all
sufficiently large real `x`.

With `u = 1 / log x` one has `1 / log (x / e) = u / (1 − u)`, and the five-term expansion
`π(y) = y · P(1 / log y) + O(y / log⁶ y)` reduces the inequality to
`(P u + C u⁶)² < u (P (u/(1−u)) − C (u/(1−u))⁶)` for small `u > 0`.  The difference of the two
sides is `u⁶ · G C u` with `G C` continuous at `0` and `G C 0 = 1` (a polynomial identity).
-/
import RamanujanLib.PiExpansion

open Real Filter Topology

noncomputable section

namespace RamanujanLib

/-- The polynomial with `u · P(u/(1−u)) − P(u)² = u⁶ · Q(u) / (1−u)⁵`.  Its constant term `1` is
the coefficient of `x² / log⁶ x` in Ramanujan's computation. -/
def Q (u : ℝ) : ℝ :=
  1 + 129 * u - 275 * u ^ 2 + 217 * u ^ 3 - 40 * u ^ 4 + 1024 * u ^ 5 - 3468 * u ^ 6
    + 4452 * u ^ 7 - 2592 * u ^ 8 + 576 * u ^ 9

/-- The normalised gap, including the `O`-terms with constant `C`. -/
def G (C u : ℝ) : ℝ :=
  Q u / (1 - u) ^ 5 - C * u / (1 - u) ^ 6 - 2 * C * P u - C ^ 2 * u ^ 6

theorem G_zero (C : ℝ) : G C 0 = 1 := by simp [G, Q, P]

theorem continuousAt_G (C : ℝ) : ContinuousAt (G C) 0 := by
  have h5 : (1 - (0 : ℝ)) ^ 5 ≠ 0 := by norm_num
  have h6 : (1 - (0 : ℝ)) ^ 6 ≠ 0 := by norm_num
  unfold G Q P
  fun_prop (disch := assumption)

/-- The algebraic identity behind Ramanujan's inequality. -/
theorem key_identity (C u : ℝ) (hu : 1 - u ≠ 0) :
    u * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6) - (P u + C * u ^ 6) ^ 2 = u ^ 6 * G C u := by
  unfold G Q P
  field_simp
  ring

/-- Pointwise step: the expansion with constant `C` gives the inequality at every `x ≥ 2e` at
which the normalised gap is positive. -/
theorem ramanujan_of_expansion (C : ℝ)
    (hC : ∀ y : ℝ, 2 ≤ y →
      |(Nat.primeCounting ⌊y⌋₊ : ℝ) - y * P (log y)⁻¹| ≤ C * (y / log y ^ 6))
    (x : ℝ) (hx : 2 * exp 1 ≤ x) (hG : 0 < G C (log x)⁻¹) :
    (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
  have he : 0 < exp 1 := exp_pos 1
  have he2 : (2 : ℝ) ≤ exp 1 := by
    have := Real.add_one_le_exp (1 : ℝ)
    linarith
  have hx2 : 2 ≤ x := by nlinarith
  have hx0 : 0 < x := by linarith
  have hxe : 2 ≤ x / exp 1 := by
    rw [le_div_iff₀ he]
    exact hx
  have hlog : log (x / exp 1) = log x - 1 := by
    rw [log_div hx0.ne' he.ne', log_exp]
  have hL1 : 0 < log x - 1 := by
    rw [← hlog]
    exact log_pos (by linarith)
  have hL : 0 < log x := by linarith
  set L : ℝ := log x with hLdef
  set u : ℝ := L⁻¹ with hu
  have hu0 : 0 < u := inv_pos.mpr hL
  have h1u : 1 - u ≠ 0 := by
    have h : 1 - u = (L - 1) / L := by
      rw [hu]
      field_simp
    rw [h]
    exact div_ne_zero hL1.ne' hL.ne'
  have hv : (L - 1)⁻¹ = u / (1 - u) := by
    rw [hu]
    field_simp
  -- upper bound for `π(x)`
  have h1 : (Nat.primeCounting ⌊x⌋₊ : ℝ) ≤ x * (P u + C * u ^ 6) := by
    have h := (abs_le.mp (hC x hx2)).2
    have : x / L ^ 6 = x * u ^ 6 := by rw [hu, inv_pow, div_eq_mul_inv]
    rw [this] at h
    linarith
  -- lower bound for `π(x / e)`
  have h2 : x / exp 1 * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6)
      ≤ (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
    have h := (abs_le.mp (hC (x / exp 1) hxe)).1
    rw [hlog, hv] at h
    have : x / exp 1 / (L - 1) ^ 6 = x / exp 1 * (u / (1 - u)) ^ 6 := by
      rw [← hv, inv_pow, div_eq_mul_inv]
    rw [this] at h
    linarith
  -- the inequality between the main terms
  have hkey : (P u + C * u ^ 6) ^ 2 < u * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6) := by
    have hid := key_identity C u h1u
    have hpos : 0 < u ^ 6 * G C u := mul_pos (pow_pos hu0 6) hG
    linarith
  have hπ0 : (0 : ℝ) ≤ (Nat.primeCounting ⌊x⌋₊ : ℝ) := Nat.cast_nonneg _
  calc (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      ≤ (x * (P u + C * u ^ 6)) ^ 2 := pow_le_pow_left₀ hπ0 h1 2
    _ = x ^ 2 * (P u + C * u ^ 6) ^ 2 := by ring
    _ < x ^ 2 * (u * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6)) := by
        gcongr
    _ = exp 1 * x / L * (x / exp 1 * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6)) := by
        rw [hu]
        field_simp
    _ ≤ exp 1 * x / L * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
        gcongr

/-- **Ramanujan's inequality for the prime counting function**, for all sufficiently large real
`x`:  `π(x)² < (e x / log x) · π(x / e)`. -/
theorem ramanujan_eventually :
    ∃ x₀ : ℝ, ∀ x : ℝ, x₀ ≤ x →
      (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        < exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
  obtain ⟨C, hC⟩ := exists_primeCounting_expansion
  have hev : ∀ᶠ u : ℝ in 𝓝 0, 0 < G C u :=
    (continuousAt_G C).eventually (lt_mem_nhds (by rw [G_zero]; norm_num))
  have hinv : Tendsto (fun x : ℝ => (log x)⁻¹) atTop (𝓝 0) :=
    tendsto_log_atTop.inv_tendsto_atTop
  obtain ⟨x₀, hx₀⟩ :=
    eventually_atTop.mp ((hinv.eventually hev).and (eventually_ge_atTop (2 * exp 1)))
  exact ⟨x₀, fun x hx => ramanujan_of_expansion C hC x (hx₀ x hx).2 (hx₀ x hx).1⟩

/-! ### The size of the gap: Ramanujan's asymptotic `x² / log⁶ x` -/

/-- Two-sided bounds for the gap `(e x / log x) π(x/e) − π(x)²` in terms of `G (±C)`. -/
theorem gap_bounds (C : ℝ)
    (hC : ∀ y : ℝ, 2 ≤ y →
      |(Nat.primeCounting ⌊y⌋₊ : ℝ) - y * P (log y)⁻¹| ≤ C * (y / log y ^ 6))
    (x : ℝ) (hx : 2 * exp 1 ≤ x) (hP : 0 ≤ P (log x)⁻¹ - C * (log x)⁻¹ ^ 6) :
    x ^ 2 * (log x)⁻¹ ^ 6 * G C (log x)⁻¹
        ≤ exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ)
          - (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2 ∧
      exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ)
          - (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        ≤ x ^ 2 * (log x)⁻¹ ^ 6 * G (-C) (log x)⁻¹ := by
  have he : 0 < exp 1 := exp_pos 1
  have he2 : (2 : ℝ) ≤ exp 1 := by
    have := Real.add_one_le_exp (1 : ℝ)
    linarith
  have hx2 : 2 ≤ x := by nlinarith
  have hx0 : 0 < x := by linarith
  have hxe : 2 ≤ x / exp 1 := by
    rw [le_div_iff₀ he]
    exact hx
  have hlog : log (x / exp 1) = log x - 1 := by
    rw [log_div hx0.ne' he.ne', log_exp]
  have hL1 : 0 < log x - 1 := by
    rw [← hlog]
    exact log_pos (by linarith)
  have hL : 0 < log x := by linarith
  set L : ℝ := log x with hLdef
  set u : ℝ := L⁻¹ with hu
  have hu0 : 0 < u := inv_pos.mpr hL
  have h1u : 1 - u ≠ 0 := by
    have h : 1 - u = (L - 1) / L := by
      rw [hu]
      field_simp
    rw [h]
    exact div_ne_zero hL1.ne' hL.ne'
  have hv : (L - 1)⁻¹ = u / (1 - u) := by
    rw [hu]
    field_simp
  have hxu : x / L ^ 6 = x * u ^ 6 := by rw [hu, inv_pow, div_eq_mul_inv]
  have hxv : x / exp 1 / (L - 1) ^ 6 = x / exp 1 * (u / (1 - u)) ^ 6 := by
    rw [← hv, inv_pow, div_eq_mul_inv]
  have hA := abs_le.mp (hC x hx2)
  have hB := abs_le.mp (hC (x / exp 1) hxe)
  rw [hxu] at hA
  rw [hlog, hv, hxv] at hB
  have h1 : (Nat.primeCounting ⌊x⌋₊ : ℝ) ≤ x * (P u + C * u ^ 6) := by linarith [hA.2]
  have h1' : x * (P u - C * u ^ 6) ≤ (Nat.primeCounting ⌊x⌋₊ : ℝ) := by linarith [hA.1]
  have h2 : x / exp 1 * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6)
      ≤ (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by linarith [hB.1]
  have h2' : (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ)
      ≤ x / exp 1 * (P (u / (1 - u)) + C * (u / (1 - u)) ^ 6) := by linarith [hB.2]
  have hπ0 : (0 : ℝ) ≤ (Nat.primeCounting ⌊x⌋₊ : ℝ) := Nat.cast_nonneg _
  have hcoef : 0 ≤ exp 1 * x / L := by positivity
  have hsq : (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2 ≤ (x * (P u + C * u ^ 6)) ^ 2 :=
    pow_le_pow_left₀ hπ0 h1 2
  have hsq' : (x * (P u - C * u ^ 6)) ^ 2 ≤ (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2 :=
    pow_le_pow_left₀ (mul_nonneg hx0.le hP) h1' 2
  have hm := mul_le_mul_of_nonneg_left h2 hcoef
  have hm' := mul_le_mul_of_nonneg_left h2' hcoef
  have e1 : exp 1 * x / L * (x / exp 1 * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6))
      = x ^ 2 * (u * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6)) := by
    rw [hu]
    field_simp
  have e2 : exp 1 * x / L * (x / exp 1 * (P (u / (1 - u)) + C * (u / (1 - u)) ^ 6))
      = x ^ 2 * (u * (P (u / (1 - u)) + C * (u / (1 - u)) ^ 6)) := by
    rw [hu]
    field_simp
  have id1 := key_identity C u h1u
  have id2 := key_identity (-C) u h1u
  constructor
  · have : x ^ 2 * u ^ 6 * G C u
        = x ^ 2 * (u * (P (u / (1 - u)) - C * (u / (1 - u)) ^ 6))
          - (x * (P u + C * u ^ 6)) ^ 2 := by
      rw [mul_assoc, ← id1]
      ring
    rw [this, ← e1]
    linarith
  · have : x ^ 2 * u ^ 6 * G (-C) u
        = x ^ 2 * (u * (P (u / (1 - u)) + C * (u / (1 - u)) ^ 6))
          - (x * (P u - C * u ^ 6)) ^ 2 := by
      rw [mul_assoc, ← id2]
      ring
    rw [this, ← e2]
    linarith

/-- **Ramanujan's asymptotic for the gap:**
`(e x / log x) π(x/e) − π(x)² = (1 + o(1)) · x² / log⁶ x`. -/
theorem ramanujan_gap_tendsto :
    Tendsto (fun x : ℝ =>
      (exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ)
        - (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2) / (x ^ 2 / log x ^ 6)) atTop (𝓝 1) := by
  obtain ⟨C, hC⟩ := exists_primeCounting_expansion
  have hinv : Tendsto (fun x : ℝ => (log x)⁻¹) atTop (𝓝 0) :=
    tendsto_log_atTop.inv_tendsto_atTop
  have hlo : Tendsto (fun x : ℝ => G C (log x)⁻¹) atTop (𝓝 1) := by
    have := (continuousAt_G C).tendsto.comp hinv
    rwa [G_zero] at this
  have hhi : Tendsto (fun x : ℝ => G (-C) (log x)⁻¹) atTop (𝓝 1) := by
    have := (continuousAt_G (-C)).tendsto.comp hinv
    rwa [G_zero] at this
  -- eventually `P u − C u⁶ ≥ 0`
  have hH : ContinuousAt (fun u : ℝ => 1 + u + 2 * u ^ 2 + 6 * u ^ 3 + 24 * u ^ 4 - C * u ^ 5) 0 := by
    fun_prop
  have hHev : ∀ᶠ u : ℝ in 𝓝 0, 0 < 1 + u + 2 * u ^ 2 + 6 * u ^ 3 + 24 * u ^ 4 - C * u ^ 5 :=
    hH.eventually (lt_mem_nhds (by norm_num))
  have hev := (hinv.eventually hHev).and (eventually_ge_atTop (2 * exp 1))
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlo hhi ?_ ?_
  · filter_upwards [hev] with x ⟨hH, hx⟩
    have he2 : (2 : ℝ) ≤ exp 1 := by
      have := Real.add_one_le_exp (1 : ℝ)
      linarith
    have hx1 : 1 < x := by nlinarith
    have hL : 0 < log x := log_pos hx1
    have hP : 0 ≤ P (log x)⁻¹ - C * (log x)⁻¹ ^ 6 := by
      have : P (log x)⁻¹ - C * (log x)⁻¹ ^ 6 = (log x)⁻¹ *
          (1 + (log x)⁻¹ + 2 * (log x)⁻¹ ^ 2 + 6 * (log x)⁻¹ ^ 3 + 24 * (log x)⁻¹ ^ 4
            - C * (log x)⁻¹ ^ 5) := by
        unfold P
        ring
      rw [this]
      exact mul_nonneg (inv_pos.mpr hL).le hH.le
    have hb := (gap_bounds C hC x hx hP).1
    have hden : 0 < x ^ 2 / log x ^ 6 := by positivity
    rw [le_div_iff₀ hden]
    have : G C (log x)⁻¹ * (x ^ 2 / log x ^ 6) = x ^ 2 * (log x)⁻¹ ^ 6 * G C (log x)⁻¹ := by
      rw [inv_pow]
      ring
    rw [this]
    exact hb
  · filter_upwards [hev] with x ⟨hH, hx⟩
    have he2 : (2 : ℝ) ≤ exp 1 := by
      have := Real.add_one_le_exp (1 : ℝ)
      linarith
    have hx1 : 1 < x := by nlinarith
    have hL : 0 < log x := log_pos hx1
    have hP : 0 ≤ P (log x)⁻¹ - C * (log x)⁻¹ ^ 6 := by
      have : P (log x)⁻¹ - C * (log x)⁻¹ ^ 6 = (log x)⁻¹ *
          (1 + (log x)⁻¹ + 2 * (log x)⁻¹ ^ 2 + 6 * (log x)⁻¹ ^ 3 + 24 * (log x)⁻¹ ^ 4
            - C * (log x)⁻¹ ^ 5) := by
        unfold P
        ring
      rw [this]
      exact mul_nonneg (inv_pos.mpr hL).le hH.le
    have hb := (gap_bounds C hC x hx hP).2
    have hden : 0 < x ^ 2 / log x ^ 6 := by positivity
    rw [div_le_iff₀ hden]
    have : G (-C) (log x)⁻¹ * (x ^ 2 / log x ^ 6)
        = x ^ 2 * (log x)⁻¹ ^ 6 * G (-C) (log x)⁻¹ := by
      rw [inv_pow]
      ring
    rw [this]
    exact hb

end RamanujanLib
