/-
AllxLib/Analytic/KernelK.lean — the test function `Kt t` of the explicit formula (Mathlib only).

    Kt t = (Bk ½ (t+1) − 2 Bk ½ t + Bk ½ (t−1))/6 − t · Bk ½ 1 + Bk (3/2) 1,

i.e. `Kt t u = e^{u/2} (H(t − e^u) − (t − e^u)(1 − e^u)₊³)`, `H = Hhat`.  For `t ≥ 2`:
`Kt t (log n) = √n · H(t − n)`,  `Kt t (−log n) = (t − 1/n)(1 − (1 − 1/n)³)/√n`, and with
`s = 1/2 + ir`

    ∫ Kt t u e^{iru} du = ((t+1)^{s+3} − 2 t^{s+3} + (t−1)^{s+3} − 6t)/(s(s+1)(s+2)(s+3))
                           + 6/((s+1)(s+2)(s+3)(s+4)),

whose norm is at most
`(t+1)^{3/2} (1/(9/4+r²) + (5/12) m(r)) + (3t/7 + 1/5) m(r)`,  `m(r) = 30/((9/4+r²)(25/4+r²))`.
-/
import AllxLib.Analytic.Kernel

open Real MeasureTheory Set

noncomputable section

namespace Allx

/-! ### The hat average of the positive part -/

/-- `H(y) = ((y+1)₊³ − 2 y₊³ + (y−1)₊³)/6`. -/
def Hhat (y : ℝ) : ℝ := ((max (y + 1) 0) ^ 3 - 2 * (max y 0) ^ 3 + (max (y - 1) 0) ^ 3) / 6

lemma Hhat_of_one_le {y : ℝ} (h : 1 ≤ y) : Hhat y = y := by
  unfold Hhat
  rw [max_eq_left (by linarith), max_eq_left (by linarith), max_eq_left (by linarith)]
  ring

lemma Hhat_of_le_neg_one {y : ℝ} (h : y ≤ -1) : Hhat y = 0 := by
  unfold Hhat
  rw [max_eq_right (by linarith), max_eq_right (by linarith), max_eq_right (by linarith)]
  ring

/-- `0 ≤ H(y) − y₊ ≤ 1/6`. -/
lemma Hhat_sub_max (y : ℝ) : 0 ≤ Hhat y - max y 0 ∧ Hhat y - max y 0 ≤ 1 / 6 := by
  rcases le_total y (-1) with h1 | h1
  · rw [Hhat_of_le_neg_one h1, max_eq_right (by linarith)]; norm_num
  rcases le_total y 0 with h2 | h2
  · have e : Hhat y - max y 0 = (y + 1) ^ 3 / 6 := by
      unfold Hhat
      rw [max_eq_left (by linarith : (0 : ℝ) ≤ y + 1), max_eq_right h2,
        max_eq_right (by linarith : y - 1 ≤ 0)]
      ring
    rw [e]
    have a1 : (0 : ℝ) ≤ (y + 1) ^ 3 := pow_nonneg (by linarith) 3
    have a2 : (y + 1) ^ 3 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    constructor <;> linarith
  rcases le_total y 1 with h3 | h3
  · have e : Hhat y - max y 0 = (1 - y) ^ 3 / 6 := by
      unfold Hhat
      rw [max_eq_left (by linarith : (0 : ℝ) ≤ y + 1), max_eq_left h2,
        max_eq_right (by linarith : y - 1 ≤ 0)]
      ring
    rw [e]
    have a1 : (0 : ℝ) ≤ (1 - y) ^ 3 := pow_nonneg (by linarith) 3
    have a2 : (1 - y) ^ 3 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    constructor <;> linarith
  · rw [Hhat_of_one_le h3, max_eq_left h2]; norm_num

lemma max_le_Hhat (y : ℝ) : max y 0 ≤ Hhat y := by linarith [(Hhat_sub_max y).1]

lemma Hhat_nonneg (y : ℝ) : 0 ≤ Hhat y := (le_max_right y 0).trans (max_le_Hhat y)

/-! ### The test function -/

/-- The test function. -/
def Kt (t u : ℝ) : ℝ :=
  (Bk (1 / 2) (t + 1) u - 2 * Bk (1 / 2) t u + Bk (1 / 2) (t - 1) u) / 6
    - t * Bk (1 / 2) 1 u + Bk (3 / 2) 1 u

/-- Its first derivative. -/
def Kt1 (t u : ℝ) : ℝ :=
  (Bk1 (1 / 2) (t + 1) u - 2 * Bk1 (1 / 2) t u + Bk1 (1 / 2) (t - 1) u) / 6
    - t * Bk1 (1 / 2) 1 u + Bk1 (3 / 2) 1 u

/-- Its second derivative. -/
def Kt2 (t u : ℝ) : ℝ :=
  (Bk2 (1 / 2) (t + 1) u - 2 * Bk2 (1 / 2) t u + Bk2 (1 / 2) (t - 1) u) / 6
    - t * Bk2 (1 / 2) 1 u + Bk2 (3 / 2) 1 u

theorem hasDerivAt_Kt (t u : ℝ) : HasDerivAt (Kt t) (Kt1 t u) u := by
  have h := (((((hasDerivAt_Bk (1 / 2) (t + 1) u).sub
    ((hasDerivAt_Bk (1 / 2) t u).const_mul 2)).add (hasDerivAt_Bk (1 / 2) (t - 1) u)).div_const 6).sub
    ((hasDerivAt_Bk (1 / 2) 1 u).const_mul t)).add (hasDerivAt_Bk (3 / 2) 1 u)
  exact h

theorem hasDerivAt_Kt1 (t u : ℝ) : HasDerivAt (Kt1 t) (Kt2 t u) u := by
  have h := (((((hasDerivAt_Bk1 (1 / 2) (t + 1) u).sub
    ((hasDerivAt_Bk1 (1 / 2) t u).const_mul 2)).add (hasDerivAt_Bk1 (1 / 2) (t - 1) u)).div_const 6).sub
    ((hasDerivAt_Bk1 (1 / 2) 1 u).const_mul t)).add (hasDerivAt_Bk1 (3 / 2) 1 u)
  exact h

theorem continuous_Kt (t : ℝ) : Continuous (Kt t) := by
  have := continuous_Bk (1 / 2) (t + 1)
  have := continuous_Bk (1 / 2) t
  have := continuous_Bk (1 / 2) (t - 1)
  have := continuous_Bk (1 / 2) 1
  have := continuous_Bk (3 / 2) 1
  unfold Kt; fun_prop

theorem continuous_Kt2 (t : ℝ) : Continuous (Kt2 t) := by
  have := continuous_Bk2 (1 / 2) (t + 1)
  have := continuous_Bk2 (1 / 2) t
  have := continuous_Bk2 (1 / 2) (t - 1)
  have := continuous_Bk2 (1 / 2) 1
  have := continuous_Bk2 (3 / 2) 1
  unfold Kt2; fun_prop

theorem integrable_Kt {t : ℝ} (ht : 1 < t) : Integrable (Kt t) := by
  have h := (((((integrable_Bk (σ := 1 / 2) (a := t + 1) (by norm_num) (by linarith)).sub
    ((integrable_Bk (σ := 1 / 2) (a := t) (by norm_num) (by linarith)).const_mul 2)).add
    (integrable_Bk (σ := 1 / 2) (a := t - 1) (by norm_num) (by linarith))).div_const 6).sub
    ((integrable_Bk (σ := 1 / 2) (a := 1) (by norm_num) (by norm_num)).const_mul t)).add
    (integrable_Bk (σ := 3 / 2) (a := 1) (by norm_num) (by norm_num))
  exact h

theorem integrable_Kt1 {t : ℝ} (ht : 1 < t) : Integrable (Kt1 t) := by
  have h := (((((integrable_Bk1 (σ := 1 / 2) (a := t + 1) (by norm_num) (by linarith)).sub
    ((integrable_Bk1 (σ := 1 / 2) (a := t) (by norm_num) (by linarith)).const_mul 2)).add
    (integrable_Bk1 (σ := 1 / 2) (a := t - 1) (by norm_num) (by linarith))).div_const 6).sub
    ((integrable_Bk1 (σ := 1 / 2) (a := 1) (by norm_num) (by norm_num)).const_mul t)).add
    (integrable_Bk1 (σ := 3 / 2) (a := 1) (by norm_num) (by norm_num))
  exact h

theorem integrable_Kt2 {t : ℝ} (ht : 1 < t) : Integrable (Kt2 t) := by
  have h := (((((integrable_Bk2 (σ := 1 / 2) (a := t + 1) (by norm_num) (by linarith)).sub
    ((integrable_Bk2 (σ := 1 / 2) (a := t) (by norm_num) (by linarith)).const_mul 2)).add
    (integrable_Bk2 (σ := 1 / 2) (a := t - 1) (by norm_num) (by linarith))).div_const 6).sub
    ((integrable_Bk2 (σ := 1 / 2) (a := 1) (by norm_num) (by norm_num)).const_mul t)).add
    (integrable_Bk2 (σ := 3 / 2) (a := 1) (by norm_num) (by norm_num))
  exact h

/-! ### Pointwise formulas -/

/-- `Kt t u = e^{u/2} (H(t − e^u) − (t − e^u)(1 − e^u)₊³)`. -/
theorem Kt_eq (t u : ℝ) :
    Kt t u = exp (u / 2) * (Hhat (t - exp u) - (t - exp u) * (mk 1 u) ^ 3) := by
  unfold Kt Bk Hhat mk
  rw [show (3 / 2 : ℝ) * u = u / 2 + u by ring, exp_add, show (1 / 2 : ℝ) * u = u / 2 by ring,
    show t + 1 - exp u = t - exp u + 1 by ring, show t - 1 - exp u = t - exp u - 1 by ring]
  ring

theorem Kt_of_nonneg (t : ℝ) {u : ℝ} (hu : 0 ≤ u) : Kt t u = exp (u / 2) * Hhat (t - exp u) := by
  rw [Kt_eq, mk_of_ge (by simpa using Real.one_le_exp hu)]
  ring

theorem Kt_of_nonpos {t : ℝ} (ht : 2 ≤ t) {u : ℝ} (hu : u ≤ 0) :
    Kt t u = exp (u / 2) * ((t - exp u) * (1 - (1 - exp u) ^ 3)) := by
  have he : exp u ≤ 1 := by simpa using Real.exp_le_one_iff.mpr hu
  rw [Kt_eq, mk_of_le he, Hhat_of_one_le (by linarith)]
  ring

theorem Kt_nonneg {t : ℝ} (ht : 2 ≤ t) (u : ℝ) : 0 ≤ Kt t u := by
  rcases le_total 0 u with hu | hu
  · rw [Kt_of_nonneg t hu]
    exact mul_nonneg (exp_pos _).le (Hhat_nonneg _)
  · have he : exp u ≤ 1 := by simpa using Real.exp_le_one_iff.mpr hu
    have he0 := exp_pos u
    rw [Kt_of_nonpos ht hu]
    have h1 : (1 - exp u) ^ 3 ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
    exact mul_nonneg (exp_pos _).le (mul_nonneg (by linarith) (by linarith))

theorem Kt_zero {t : ℝ} (ht : 2 ≤ t) : Kt t 0 = t - 1 := by
  rw [Kt_of_nonneg t le_rfl]
  simp only [zero_div, exp_zero, one_mul]
  exact Hhat_of_one_le (by linarith)

/-- For `u ≥ 0`:  `e^{u/2}(t − e^u) ≤ Kt t u`. -/
theorem le_Kt_of_nonneg (t : ℝ) {u : ℝ} (hu : 0 ≤ u) : exp (u / 2) * (t - exp u) ≤ Kt t u := by
  rw [Kt_of_nonneg t hu]
  exact mul_le_mul_of_nonneg_left ((le_max_left _ _).trans (max_le_Hhat _)) (exp_pos _).le

theorem Kt_log (t : ℝ) {n : ℕ} (hn : 0 < n) : Kt t (log n) = √n * Hhat (t - n) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  rw [Kt_of_nonneg t (Real.log_nonneg h1), Real.exp_half, exp_log hn']

theorem Kt_neg_log {t : ℝ} (ht : 2 ≤ t) {n : ℕ} (hn : 0 < n) :
    Kt t (-log n) = (√n)⁻¹ * ((t - (n : ℝ)⁻¹) * (1 - (1 - (n : ℝ)⁻¹) ^ 3)) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  rw [Kt_of_nonpos ht (by linarith [Real.log_nonneg h1]), Real.exp_half, Real.exp_neg, exp_log hn',
    Real.sqrt_inv]

/-! ### The Fourier transform -/

section FT

open Complex

theorem integrable_Bk_mul_cexp {σ a : ℝ} (hσ : 0 < σ) (ha : 0 < a) (r : ℝ) :
    Integrable (fun u : ℝ => ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u)) := by
  have h := (integrable_Bk hσ ha).ofReal (𝕜 := ℂ)
  refine h.mul_bdd (c := 1) (by fun_prop) ?_
  refine Filter.Eventually.of_forall fun u => ?_
  rw [show I * (r : ℂ) * (u : ℂ) = ((r * u : ℝ) : ℂ) * I by push_cast; ring,
    Complex.norm_exp_ofReal_mul_I]

/-- The Fourier transform of `Kt t`, with `s = 1/2 + ir`. -/
theorem integral_Kt_mul_cexp {t : ℝ} (ht : 1 < t) (r : ℝ) :
    ∫ u : ℝ, ((Kt t u : ℝ) : ℂ) * cexp (I * r * u)
      = (((t + 1 : ℝ) : ℂ) ^ ((1 / 2 : ℂ) + I * r + 3) - 2 * (t : ℂ) ^ ((1 / 2 : ℂ) + I * r + 3)
          + ((t - 1 : ℝ) : ℂ) ^ ((1 / 2 : ℂ) + I * r + 3) - 6 * t)
          / (((1 / 2 : ℂ) + I * r) * ((1 / 2 : ℂ) + I * r + 1) * ((1 / 2 : ℂ) + I * r + 2)
            * ((1 / 2 : ℂ) + I * r + 3))
        + 6 / (((1 / 2 : ℂ) + I * r + 1) * ((1 / 2 : ℂ) + I * r + 2) * ((1 / 2 : ℂ) + I * r + 3)
            * ((1 / 2 : ℂ) + I * r + 4)) := by
  set s : ℂ := (1 / 2 : ℂ) + I * r with hs
  have i1 := integrable_Bk_mul_cexp (σ := 1 / 2) (a := t + 1) (by norm_num) (by linarith) r
  have i2 := integrable_Bk_mul_cexp (σ := 1 / 2) (a := t) (by norm_num) (by linarith) r
  have i3 := integrable_Bk_mul_cexp (σ := 1 / 2) (a := t - 1) (by norm_num) (by linarith) r
  have i4 := integrable_Bk_mul_cexp (σ := 1 / 2) (a := 1) (by norm_num) (by norm_num) r
  have i5 := integrable_Bk_mul_cexp (σ := 3 / 2) (a := 1) (by norm_num) (by norm_num) r
  have v1 := integral_Bk_mul_cexp (σ := 1 / 2) (a := t + 1) (by norm_num) (by linarith) r
  have v2 := integral_Bk_mul_cexp (σ := 1 / 2) (a := t) (by norm_num) (by linarith) r
  have v3 := integral_Bk_mul_cexp (σ := 1 / 2) (a := t - 1) (by norm_num) (by linarith) r
  have v4 := integral_Bk_mul_cexp (σ := 1 / 2) (a := 1) (by norm_num) (by norm_num) r
  have v5 := integral_Bk_mul_cexp (σ := 3 / 2) (a := 1) (by norm_num) (by norm_num) r
  have e : ∀ u : ℝ, ((Kt t u : ℝ) : ℂ) * cexp (I * r * u)
      = (1 / 6 : ℂ) * (((Bk (1 / 2) (t + 1) u : ℝ) : ℂ) * cexp (I * r * u))
        - (1 / 3 : ℂ) * (((Bk (1 / 2) t u : ℝ) : ℂ) * cexp (I * r * u))
        + (1 / 6 : ℂ) * (((Bk (1 / 2) (t - 1) u : ℝ) : ℂ) * cexp (I * r * u))
        - (t : ℂ) * (((Bk (1 / 2) 1 u : ℝ) : ℂ) * cexp (I * r * u))
        + ((Bk (3 / 2) 1 u : ℝ) : ℂ) * cexp (I * r * u) := by
    intro u
    unfold Kt
    push_cast
    ring
  simp_rw [e]
  have j1 := i1.const_mul (1 / 6 : ℂ)
  have j2 := i2.const_mul (1 / 3 : ℂ)
  have j3 := i3.const_mul (1 / 6 : ℂ)
  have j4 := i4.const_mul (t : ℂ)
  have j12 : Integrable (fun u : ℝ => (1 / 6 : ℂ) * (((Bk (1 / 2) (t + 1) u : ℝ) : ℂ) * cexp (I * r * u))
        - (1 / 3 : ℂ) * (((Bk (1 / 2) t u : ℝ) : ℂ) * cexp (I * r * u))) := j1.sub j2
  have j123 : Integrable (fun u : ℝ => (1 / 6 : ℂ) * (((Bk (1 / 2) (t + 1) u : ℝ) : ℂ) * cexp (I * r * u))
        - (1 / 3 : ℂ) * (((Bk (1 / 2) t u : ℝ) : ℂ) * cexp (I * r * u))
        + (1 / 6 : ℂ) * (((Bk (1 / 2) (t - 1) u : ℝ) : ℂ) * cexp (I * r * u))) := j12.add j3
  have j1234 : Integrable (fun u : ℝ => (1 / 6 : ℂ) * (((Bk (1 / 2) (t + 1) u : ℝ) : ℂ) * cexp (I * r * u))
        - (1 / 3 : ℂ) * (((Bk (1 / 2) t u : ℝ) : ℂ) * cexp (I * r * u))
        + (1 / 6 : ℂ) * (((Bk (1 / 2) (t - 1) u : ℝ) : ℂ) * cexp (I * r * u))
        - (t : ℂ) * (((Bk (1 / 2) 1 u : ℝ) : ℂ) * cexp (I * r * u))) := j123.sub j4
  rw [integral_add j1234 i5, integral_sub j123 j4, integral_add j12 j3, integral_sub j1 j2,
    integral_const_mul, integral_const_mul, integral_const_mul, integral_const_mul,
    v1, v2, v3, v4, v5]
  have hre : s.re = 1 / 2 := by simp [hs]
  have hne : ∀ j : ℕ, s + j ≠ 0 := by
    intro j h
    have := congrArg Complex.re h
    simp [hre] at this
    have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    linarith
  have h0 := hne 0
  have h1 := hne 1
  have h2 := hne 2
  have h3 := hne 3
  have h4 := hne 4
  push_cast at h0 h1 h2 h3 h4
  rw [add_zero] at h0
  have c1 : ((1 / 2 : ℝ) : ℂ) = 1 / 2 := by push_cast; ring
  have c3 : ((3 / 2 : ℝ) : ℂ) = 3 / 2 := by push_cast; ring
  rw [c1, c3, ← hs, Complex.ofReal_one, Complex.one_cpow, Complex.one_cpow,
    show (1 / 2 : ℂ) + 3 + I * r = s + 3 by rw [hs]; ring,
    show (1 / 2 : ℂ) + 1 + I * r = s + 1 by rw [hs]; ring,
    show (1 / 2 : ℂ) + 2 + I * r = s + 2 by rw [hs]; ring,
    show (3 / 2 : ℂ) + I * r = s + 1 by rw [hs]; ring,
    show (3 / 2 : ℂ) + 1 + I * r = s + 2 by rw [hs]; ring,
    show (3 / 2 : ℂ) + 2 + I * r = s + 3 by rw [hs]; ring,
    show (3 / 2 : ℂ) + 3 + I * r = s + 4 by rw [hs]; ring]
  field_simp
  ring

/-- `‖c + ir‖² = c² + r²`. -/
lemma sq_norm_add_I_mul (c r : ℝ) : ‖(c : ℂ) + I * r‖ ^ 2 = c ^ 2 + r ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp
  ring

/-- The second difference of `a ↦ a^{s+3}` at `t`, `s = 1/2 + ir`. -/
theorem norm_second_diff_cpow_le {t : ℝ} (ht : 2 ≤ t) (r : ℝ) :
    ‖((t + 1 : ℝ) : ℂ) ^ ((1 / 2 : ℂ) + I * r + 3) - 2 * (t : ℂ) ^ ((1 / 2 : ℂ) + I * r + 3)
        + ((t - 1 : ℝ) : ℂ) ^ ((1 / 2 : ℂ) + I * r + 3)‖
      ≤ ‖((1 / 2 : ℂ) + I * r + 2) * ((1 / 2 : ℂ) + I * r + 3)‖ * (t + 1) ^ (3 / 2 : ℝ) := by
  set s : ℂ := (1 / 2 : ℂ) + I * r with hs
  have hre : s.re = 1 / 2 := by simp [hs]
  have hne : ∀ j : ℕ, s + j ≠ 0 := by
    intro j h
    have := congrArg Complex.re h
    simp [hre] at this
    have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    linarith
  have h2 := hne 2
  have h3 := hne 3
  push_cast at h2 h3
  have key := norm_second_diff_le (t := t) (C := (t + 1) ^ (3 / 2 : ℝ))
    (f := fun a : ℝ => (a : ℂ) ^ (s + 3) / (s + 3) / (s + 2))
    (f' := fun a : ℝ => (a : ℂ) ^ (s + 2) / (s + 2))
    (f'' := fun a : ℝ => (a : ℂ) ^ (s + 1)) ?_ ?_ ?_
  · have e : ((t + 1 : ℝ) : ℂ) ^ (s + 3) - 2 * (t : ℂ) ^ (s + 3) + ((t - 1 : ℝ) : ℂ) ^ (s + 3)
        = ((s + 2) * (s + 3)) * (((t + 1 : ℝ) : ℂ) ^ (s + 3) / (s + 3) / (s + 2)
          - 2 * ((t : ℂ) ^ (s + 3) / (s + 3) / (s + 2))
          + ((t - 1 : ℝ) : ℂ) ^ (s + 3) / (s + 3) / (s + 2)) := by
      field_simp
    rw [e, norm_mul]
    exact mul_le_mul_of_nonneg_left key (norm_nonneg _)
  · intro x hx
    have hx0 : x ≠ 0 := by linarith [hx.1]
    have h := (hasDerivAt_ofReal_cpow_const' (x := x) hx0 (r := s + 2) (by
      intro h
      have := congrArg Complex.re h
      simp [hre] at this
      linarith)).div_const (s + 2)
    rw [show s + 2 + 1 = s + 3 by ring] at h
    exact h
  · intro x hx
    have hx0 : x ≠ 0 := by linarith [hx.1]
    have h := hasDerivAt_ofReal_cpow_const' (x := x) hx0 (r := s + 1) (by
      intro h
      have := congrArg Complex.re h
      simp [hre] at this
      linarith)
    rw [show s + 1 + 1 = s + 2 by ring] at h
    exact h
  · intro x hx
    have hx0 : 0 < x := by linarith [hx.1]
    rw [Complex.norm_cpow_eq_rpow_re_of_pos hx0]
    have : (s + 1).re = 3 / 2 := by simp [hre]; norm_num
    rw [this]
    exact Real.rpow_le_rpow hx0.le hx.2 (by norm_num)

end FT

/-! ### Real inequalities between the kernels -/

/-- `m(r) = 30/((9/4+r²)(25/4+r²))`, the transform of the `C²` majorant of lean-chebyshev. -/
def majF (r : ℝ) : ℝ := 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2))

lemma majF_nonneg (r : ℝ) : 0 ≤ majF r := by unfold majF; positivity

lemma kernel_ineq_i {p0 p1 x : ℝ} (h0 : 0 < p0) (h1 : 0 < p1) (hx : 0 ≤ x)
    (e0 : p0 ^ 2 = 1 / 4 + x) (e1 : p1 ^ 2 = 9 / 4 + x) :
    1 / (p0 * p1) ≤ 1 / (9 / 4 + x) + 5 / 12 * (30 / ((9 / 4 + x) * (25 / 4 + x))) := by
  have hA : (0 : ℝ) < 9 / 4 + x := by linarith
  have hB : (0 : ℝ) < 25 / 4 + x := by linarith
  have e : 1 / (9 / 4 + x) + 5 / 12 * (30 / ((9 / 4 + x) * (25 / 4 + x)))
      = (75 / 4 + x) / ((9 / 4 + x) * (25 / 4 + x)) := by
    field_simp
    ring
  rw [e, div_le_div_iff₀ (mul_pos h0 h1) (mul_pos hA hB), one_mul]
  -- `(9/4+x)(25/4+x) ≤ (75/4+x) p0 p1`, from the squares
  have hsq : ((9 / 4 + x) * (25 / 4 + x)) ^ 2 ≤ ((75 / 4 + x) * (p0 * p1)) ^ 2 := by
    have : ((75 / 4 + x) * (p0 * p1)) ^ 2 = (75 / 4 + x) ^ 2 * ((1 / 4 + x) * (9 / 4 + x)) := by
      rw [← e0, ← e1]; ring
    rw [this]
    nlinarith [mul_nonneg hx hx, mul_nonneg (mul_nonneg hx hx) hx,
      mul_nonneg (mul_nonneg (mul_nonneg hx hx) hx) hx]
  exact (sq_le_sq₀ (by positivity) (by positivity)).mp hsq

lemma kernel_ineq_ii {p0 p1 p2 p3 x : ℝ} (h0 : 0 < p0) (h1 : 0 < p1) (h2 : 0 < p2) (h3 : 0 < p3)
    (hx : 0 ≤ x) (e0 : p0 ^ 2 = 1 / 4 + x) (e1 : p1 ^ 2 = 9 / 4 + x) (e2 : p2 ^ 2 = 25 / 4 + x)
    (e3 : p3 ^ 2 = 49 / 4 + x) :
    1 / (p0 * p1 * p2 * p3) ≤ 1 / 14 * (30 / ((9 / 4 + x) * (25 / 4 + x))) := by
  have hA : (0 : ℝ) < 9 / 4 + x := by linarith
  have hB : (0 : ℝ) < 25 / 4 + x := by linarith
  -- `p0 p3 ≥ (7/15) p1 p2`
  have hsq : (7 / 15 * (p1 * p2)) ^ 2 ≤ (p0 * p3) ^ 2 := by
    have a : (7 / 15 * (p1 * p2)) ^ 2 = 49 / 225 * ((9 / 4 + x) * (25 / 4 + x)) := by
      rw [← e1, ← e2]; ring
    have b : (p0 * p3) ^ 2 = (1 / 4 + x) * (49 / 4 + x) := by rw [← e0, ← e3]; ring
    rw [a, b]
    nlinarith [mul_nonneg hx hx]
  have hle : 7 / 15 * (p1 * p2) ≤ p0 * p3 := (sq_le_sq₀ (by positivity) (by positivity)).mp hsq
  have e : 1 / 14 * (30 / ((9 / 4 + x) * (25 / 4 + x))) = 1 / (7 / 15 * (p1 * p2) * (p1 * p2)) := by
    rw [← e1, ← e2]
    field_simp
    ring
  rw [e]
  apply one_div_le_one_div_of_le (by positivity)
  calc 7 / 15 * (p1 * p2) * (p1 * p2) ≤ p0 * p3 * (p1 * p2) :=
        mul_le_mul_of_nonneg_right hle (by positivity)
    _ = p0 * p1 * p2 * p3 := by ring

lemma kernel_ineq_iii {p1 p2 p3 p4 x : ℝ} (h1 : 0 < p1) (h2 : 0 < p2) (h3 : 0 < p3) (h4 : 0 < p4)
    (e1 : p1 ^ 2 = 9 / 4 + x) (e2 : p2 ^ 2 = 25 / 4 + x)
    (e3 : p3 ^ 2 = 49 / 4 + x) (e4 : p4 ^ 2 = 81 / 4 + x) :
    1 / (p1 * p2 * p3 * p4) ≤ 1 / 30 * (30 / ((9 / 4 + x) * (25 / 4 + x))) := by
  have h13 : p1 ≤ p3 := (sq_le_sq₀ h1.le h3.le).mp (by rw [e1, e3]; linarith)
  have h24 : p2 ≤ p4 := (sq_le_sq₀ h2.le h4.le).mp (by rw [e2, e4]; linarith)
  have e : 1 / 30 * (30 / ((9 / 4 + x) * (25 / 4 + x))) = 1 / (p1 * p2 * p1 * p2) := by
    rw [← e1, ← e2]
    field_simp
  rw [e]
  apply one_div_le_one_div_of_le (by positivity)
  calc p1 * p2 * p1 * p2 ≤ p1 * p2 * p3 * p4 :=
    mul_le_mul (mul_le_mul_of_nonneg_left h13 (by positivity)) h24 h2.le (by positivity)

section FTbound

open Complex

/-- **The norm of the transform of `Kt t`.** -/
theorem norm_integral_Kt_mul_cexp_le {t : ℝ} (ht : 2 ≤ t) (r : ℝ) :
    ‖∫ u : ℝ, ((Kt t u : ℝ) : ℂ) * cexp (I * r * u)‖
      ≤ (t + 1) ^ (3 / 2 : ℝ) * (1 / (9 / 4 + r ^ 2) + 5 / 12 * majF r)
        + (3 / 7 * t + 1 / 5) * majF r := by
  rw [integral_Kt_mul_cexp (by linarith) r]
  set s : ℂ := (1 / 2 : ℂ) + I * r with hs
  have x0 : 0 ≤ r ^ 2 := sq_nonneg r
  have n0 : ‖s‖ ^ 2 = 1 / 4 + r ^ 2 := by
    have := sq_norm_add_I_mul (1 / 2) r
    rw [show ((1 / 2 : ℝ) : ℂ) = 1 / 2 by push_cast; ring] at this
    rw [hs, this]; norm_num
  have nj : ∀ c : ℝ, ‖s + c‖ ^ 2 = (1 / 2 + c) ^ 2 + r ^ 2 := by
    intro c
    have := sq_norm_add_I_mul (1 / 2 + c) r
    rw [show (((1 / 2 + c : ℝ)) : ℂ) + I * r = s + c by rw [hs]; push_cast; ring] at this
    exact this
  have n1 : ‖s + 1‖ ^ 2 = 9 / 4 + r ^ 2 := by
    have := nj 1; push_cast at this; rw [this]; norm_num
  have n2 : ‖s + 2‖ ^ 2 = 25 / 4 + r ^ 2 := by
    have := nj 2; push_cast at this; rw [this]; norm_num
  have n3 : ‖s + 3‖ ^ 2 = 49 / 4 + r ^ 2 := by
    have := nj 3; push_cast at this; rw [this]; norm_num
  have n4 : ‖s + 4‖ ^ 2 = 81 / 4 + r ^ 2 := by
    have := nj 4; push_cast at this; rw [this]; norm_num
  have pos : ∀ {z : ℂ} {y : ℝ}, ‖z‖ ^ 2 = y → 0 < y → 0 < ‖z‖ := by
    intro z y h hy
    have : 0 < ‖z‖ ^ 2 := by rw [h]; exact hy
    refine lt_of_le_of_ne (norm_nonneg z) fun h0 => ?_
    rw [← h0] at this
    norm_num at this
  have p0 : 0 < ‖s‖ := pos n0 (by linarith)
  have p1 : 0 < ‖s + 1‖ := pos n1 (by linarith)
  have p2 : 0 < ‖s + 2‖ := pos n2 (by linarith)
  have p3 : 0 < ‖s + 3‖ := pos n3 (by linarith)
  have p4 : 0 < ‖s + 4‖ := pos n4 (by linarith)
  have hd := norm_second_diff_cpow_le ht r
  rw [← hs] at hd
  have ht0 : (0 : ℝ) ≤ t := by linarith
  have hT : (0 : ℝ) ≤ (t + 1) ^ (3 / 2 : ℝ) := Real.rpow_nonneg (by linarith) _
  -- split
  have e : (((t + 1 : ℝ) : ℂ) ^ (s + 3) - 2 * (t : ℂ) ^ (s + 3) + ((t - 1 : ℝ) : ℂ) ^ (s + 3) - 6 * t)
        / (s * (s + 1) * (s + 2) * (s + 3)) + 6 / ((s + 1) * (s + 2) * (s + 3) * (s + 4))
      = (((t + 1 : ℝ) : ℂ) ^ (s + 3) - 2 * (t : ℂ) ^ (s + 3) + ((t - 1 : ℝ) : ℂ) ^ (s + 3))
          / (s * (s + 1) * (s + 2) * (s + 3))
        - 6 * t / (s * (s + 1) * (s + 2) * (s + 3)) + 6 / ((s + 1) * (s + 2) * (s + 3) * (s + 4)) := by
    ring
  rw [e]
  have b1 : ‖(((t + 1 : ℝ) : ℂ) ^ (s + 3) - 2 * (t : ℂ) ^ (s + 3) + ((t - 1 : ℝ) : ℂ) ^ (s + 3))
        / (s * (s + 1) * (s + 2) * (s + 3))‖ ≤ (t + 1) ^ (3 / 2 : ℝ) * (1 / (‖s‖ * ‖s + 1‖)) := by
    rw [norm_div, norm_mul, norm_mul, norm_mul]
    rw [norm_mul] at hd
    rw [div_le_iff₀ (by positivity)]
    calc _ ≤ ‖s + 2‖ * ‖s + 3‖ * (t + 1) ^ (3 / 2 : ℝ) := hd
      _ = (t + 1) ^ (3 / 2 : ℝ) * (1 / (‖s‖ * ‖s + 1‖)) * (‖s‖ * ‖s + 1‖ * ‖s + 2‖ * ‖s + 3‖) := by
          field_simp
  have b2 : ‖6 * (t : ℂ) / (s * (s + 1) * (s + 2) * (s + 3))‖
      = 6 * t * (1 / (‖s‖ * ‖s + 1‖ * ‖s + 2‖ * ‖s + 3‖)) := by
    rw [norm_div, norm_mul, norm_mul, norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg ht0]
    simp
    ring
  have b3 : ‖(6 : ℂ) / ((s + 1) * (s + 2) * (s + 3) * (s + 4))‖
      = 6 * (1 / (‖s + 1‖ * ‖s + 2‖ * ‖s + 3‖ * ‖s + 4‖)) := by
    rw [norm_div, norm_mul, norm_mul, norm_mul]
    simp
    ring
  have k1 := kernel_ineq_i p0 p1 x0 n0 n1
  have k2 := kernel_ineq_ii p0 p1 p2 p3 x0 n0 n1 n2 n3
  have k3 := kernel_ineq_iii p1 p2 p3 p4 n1 n2 n3 n4
  have hm : majF r = 30 / ((9 / 4 + r ^ 2) * (25 / 4 + r ^ 2)) := rfl
  rw [← hm] at k1 k2 k3
  calc ‖_ - 6 * (t : ℂ) / (s * (s + 1) * (s + 2) * (s + 3)) + 6 / ((s + 1) * (s + 2) * (s + 3) * (s + 4))‖
      ≤ ‖(((t + 1 : ℝ) : ℂ) ^ (s + 3) - 2 * (t : ℂ) ^ (s + 3) + ((t - 1 : ℝ) : ℂ) ^ (s + 3))
          / (s * (s + 1) * (s + 2) * (s + 3))‖
        + ‖6 * (t : ℂ) / (s * (s + 1) * (s + 2) * (s + 3))‖
        + ‖(6 : ℂ) / ((s + 1) * (s + 2) * (s + 3) * (s + 4))‖ :=
        (norm_add_le _ _).trans (add_le_add_left (norm_sub_le _ _) _)
    _ ≤ (t + 1) ^ (3 / 2 : ℝ) * (1 / (9 / 4 + r ^ 2) + 5 / 12 * majF r)
          + 6 * t * (1 / 14 * majF r) + 6 * (1 / 30 * majF r) := by
        rw [b2, b3]
        have c1 := b1.trans (mul_le_mul_of_nonneg_left k1 hT)
        have c2 : 6 * t * (1 / (‖s‖ * ‖s + 1‖ * ‖s + 2‖ * ‖s + 3‖)) ≤ 6 * t * (1 / 14 * majF r) :=
          mul_le_mul_of_nonneg_left k2 (by positivity)
        have c3 : 6 * (1 / (‖s + 1‖ * ‖s + 2‖ * ‖s + 3‖ * ‖s + 4‖)) ≤ 6 * (1 / 30 * majF r) :=
          mul_le_mul_of_nonneg_left k3 (by norm_num)
        linarith
    _ = _ := by ring

/-- A crude decay bound for the transform of `Kt t` (hypothesis of `arch_weil`). -/
theorem norm_integral_Kt_mul_cexp_mul_le {t : ℝ} (ht : 2 ≤ t) (r : ℝ) :
    ‖∫ u : ℝ, ((Kt t u : ℝ) : ℂ) * cexp (I * r * u)‖ * (1 + r ^ 2)
      ≤ 3 * (t + 1) ^ (3 / 2 : ℝ) + 24 / 5 * (3 / 7 * t + 1 / 5) := by
  have h := norm_integral_Kt_mul_cexp_le ht r
  have hT : (0 : ℝ) ≤ (t + 1) ^ (3 / 2 : ℝ) := Real.rpow_nonneg (by linarith) _
  have x0 : 0 ≤ r ^ 2 := sq_nonneg r
  have a1 : 1 / (9 / 4 + r ^ 2) * (1 + r ^ 2) ≤ 1 := by
    rw [one_div, inv_mul_le_iff₀ (by positivity)]; linarith
  have a2 : majF r * (1 + r ^ 2) ≤ 24 / 5 := by
    unfold majF
    rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    nlinarith [mul_nonneg x0 x0]
  have hmn := majF_nonneg r
  calc _ ≤ ((t + 1) ^ (3 / 2 : ℝ) * (1 / (9 / 4 + r ^ 2) + 5 / 12 * majF r)
          + (3 / 7 * t + 1 / 5) * majF r) * (1 + r ^ 2) :=
        mul_le_mul_of_nonneg_right h (by positivity)
    _ = (t + 1) ^ (3 / 2 : ℝ) * (1 / (9 / 4 + r ^ 2) * (1 + r ^ 2) + 5 / 12 * (majF r * (1 + r ^ 2)))
          + (3 / 7 * t + 1 / 5) * (majF r * (1 + r ^ 2)) := by ring
    _ ≤ (t + 1) ^ (3 / 2 : ℝ) * (1 + 5 / 12 * (24 / 5)) + (3 / 7 * t + 1 / 5) * (24 / 5) := by
        gcongr
    _ = _ := by ring

end FTbound

end Allx
