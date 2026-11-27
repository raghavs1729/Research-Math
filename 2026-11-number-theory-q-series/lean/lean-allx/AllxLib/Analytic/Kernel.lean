/-
AllxLib/Analytic/Kernel.lean — the smoothed Riesz kernel (Mathlib only).

Building block:  `Bk σ a u = e^{σu} (a − e^u)₊³`.  It is `C²` in `u`, decays like `e^{σu}` at `−∞`,
vanishes for `u ≥ log a`, and

    ∫ Bk σ a u · e^{iru} du = 6 a^{σ+3+ir} / ((σ+ir)(σ+1+ir)(σ+2+ir)(σ+3+ir)).

The test function of the explicit formula is

    Kt t = (Bk ½ (t+1) − 2 Bk ½ t + Bk ½ (t−1))/6 − t · Bk ½ 1 + Bk (3/2) 1,

i.e. `Kt t u = e^{u/2} H(t − e^u) (1 − (1 − e^u)₊³)` with `H(y) = ((y+1)₊³ − 2y₊³ + (y−1)₊³)/6`
(the average of `y₊` against the hat function on `[−1,1]`; `H(y) = y` for `y ≥ 1`).
-/
import Mathlib

open Real MeasureTheory Set

noncomputable section

namespace Allx

/-! ### Powers of the positive part -/

/-- Gluing two differentiable functions at `0` when their derivatives agree there. -/
lemma hasDerivAt_glue {F A B A' B' : ℝ → ℝ} (hA : ∀ t, HasDerivAt A (A' t) t)
    (hB : ∀ t, HasDerivAt B (B' t) t) (hFA : ∀ u, 0 ≤ u → F u = A u)
    (hFB : ∀ u, u ≤ 0 → F u = B u) (h0 : A' 0 = B' 0) (u : ℝ) :
    HasDerivAt F (if 0 ≤ u then A' u else B' u) u := by
  rcases lt_trichotomy u 0 with hu | hu | hu
  · rw [if_neg (not_le.mpr hu)]
    refine (hB u).congr_of_eventuallyEq ?_
    filter_upwards [Iio_mem_nhds hu] with x hx using hFB x (le_of_lt hx)
  · subst hu
    rw [if_pos le_rfl]
    have h1 : HasDerivWithinAt F (A' 0) (Set.Ici 0) 0 :=
      (hA 0).hasDerivWithinAt.congr (fun x hx => hFA x hx) (hFA 0 le_rfl)
    have h2 : HasDerivWithinAt F (A' 0) (Set.Iic 0) 0 := by
      rw [h0]
      exact (hB 0).hasDerivWithinAt.congr (fun x hx => hFB x hx) (hFB 0 le_rfl)
    have h3 := h2.union h1
    rwa [Set.Iic_union_Ici, hasDerivWithinAt_univ] at h3
  · rw [if_pos hu.le]
    refine (hA u).congr_of_eventuallyEq ?_
    filter_upwards [Ioi_mem_nhds hu] with x hx using hFA x (le_of_lt hx)

lemma hasDerivAt_pos3 (y : ℝ) :
    HasDerivAt (fun y : ℝ => (max y 0) ^ 3) (3 * (max y 0) ^ 2) y := by
  have h := hasDerivAt_glue (F := fun y : ℝ => (max y 0) ^ 3) (A := fun y => y ^ 3)
    (B := fun _ => 0) (A' := fun y => 3 * y ^ 2) (B' := fun _ => 0)
    (fun t => by simpa using hasDerivAt_pow 3 t) (fun t => hasDerivAt_const t 0)
    (fun u hu => by simp [max_eq_left hu]) (fun u hu => by simp [max_eq_right hu]) (by norm_num) y
  convert h using 1
  split_ifs with hy
  · rw [max_eq_left hy]
  · rw [max_eq_right (le_of_not_ge hy)]; norm_num

lemma hasDerivAt_pos2 (y : ℝ) :
    HasDerivAt (fun y : ℝ => (max y 0) ^ 2) (2 * max y 0) y := by
  have h := hasDerivAt_glue (F := fun y : ℝ => (max y 0) ^ 2) (A := fun y => y ^ 2)
    (B := fun _ => 0) (A' := fun y => 2 * y) (B' := fun _ => 0)
    (fun t => by simpa using hasDerivAt_pow 2 t) (fun t => hasDerivAt_const t 0)
    (fun u hu => by simp [max_eq_left hu]) (fun u hu => by simp [max_eq_right hu]) (by norm_num) y
  convert h using 1
  split_ifs with hy
  · rw [max_eq_left hy]
  · rw [max_eq_right (le_of_not_ge hy)]; norm_num

lemma hasDerivAt_pos1_of_pos2 (y : ℝ) :
    HasDerivAt (fun y : ℝ => (max y 0) ^ 2 * y) (3 * (max y 0) ^ 2) y := by
  have h := (hasDerivAt_pos3 y)
  have e : (fun y : ℝ => (max y 0) ^ 2 * y) = fun y : ℝ => (max y 0) ^ 3 := by
    funext z
    rcases le_total 0 z with hz | hz
    · rw [max_eq_left hz]; ring
    · rw [max_eq_right hz]; ring
  rw [e]; exact h

/-! ### The building block -/

/-- `(a − e^u)₊`. -/
def mk (a u : ℝ) : ℝ := max (a - exp u) 0

/-- `Bk σ a u = e^{σu} (a − e^u)₊³`. -/
def Bk (σ a u : ℝ) : ℝ := exp (σ * u) * (mk a u) ^ 3

/-- First derivative of `Bk σ a`. -/
def Bk1 (σ a u : ℝ) : ℝ := exp (σ * u) * (σ * (mk a u) ^ 3 - 3 * exp u * (mk a u) ^ 2)

/-- Second derivative of `Bk σ a`. -/
def Bk2 (σ a u : ℝ) : ℝ :=
  exp (σ * u) * (σ ^ 2 * (mk a u) ^ 3 - 3 * (2 * σ + 1) * exp u * (mk a u) ^ 2
    + 6 * exp u ^ 2 * mk a u)

lemma mk_nonneg (a u : ℝ) : 0 ≤ mk a u := le_max_right _ _

lemma mk_le {a : ℝ} (ha : 0 ≤ a) (u : ℝ) : mk a u ≤ a :=
  max_le (by linarith [exp_pos u]) ha

lemma mk_of_le {a u : ℝ} (h : exp u ≤ a) : mk a u = a - exp u := max_eq_left (by linarith)

lemma mk_of_ge {a u : ℝ} (h : a ≤ exp u) : mk a u = 0 := max_eq_right (by linarith)

lemma mk_of_log_lt {a u : ℝ} (ha : 0 < a) (h : log a < u) : mk a u = 0 := by
  apply mk_of_ge
  have := exp_lt_exp.mpr h
  rw [exp_log ha] at this
  exact this.le

lemma exp_le_of_le_log {a u : ℝ} (ha : 0 < a) (h : u ≤ log a) : exp u ≤ a := by
  have := exp_le_exp.mpr h
  rwa [exp_log ha] at this

lemma hasDerivAt_mk3 (a u : ℝ) :
    HasDerivAt (fun u => (mk a u) ^ 3) (-(3 * exp u * (mk a u) ^ 2)) u := by
  have hv : HasDerivAt (fun u : ℝ => a - exp u) (-exp u) u := by
    simpa using (hasDerivAt_exp u).const_sub a
  have h := (hasDerivAt_pos3 (a - exp u)).comp u hv
  refine HasDerivAt.congr_deriv (f := fun u => (mk a u) ^ 3) h ?_
  unfold mk; ring

lemma hasDerivAt_mk2 (a u : ℝ) :
    HasDerivAt (fun u => (mk a u) ^ 2) (-(2 * exp u * mk a u)) u := by
  have hv : HasDerivAt (fun u : ℝ => a - exp u) (-exp u) u := by
    simpa using (hasDerivAt_exp u).const_sub a
  have h := (hasDerivAt_pos2 (a - exp u)).comp u hv
  refine HasDerivAt.congr_deriv (f := fun u => (mk a u) ^ 2) h ?_
  unfold mk; ring

lemma hasDerivAt_exp_mul (σ u : ℝ) : HasDerivAt (fun u : ℝ => exp (σ * u)) (σ * exp (σ * u)) u := by
  have h := ((hasDerivAt_id' u).const_mul σ).exp
  refine h.congr_deriv ?_
  ring

theorem hasDerivAt_Bk (σ a u : ℝ) : HasDerivAt (Bk σ a) (Bk1 σ a u) u := by
  have h := (hasDerivAt_exp_mul σ u).mul (hasDerivAt_mk3 a u)
  refine HasDerivAt.congr_deriv (f := Bk σ a) h ?_
  unfold Bk1; ring

theorem hasDerivAt_Bk1 (σ a u : ℝ) : HasDerivAt (Bk1 σ a) (Bk2 σ a u) u := by
  have h1 : HasDerivAt (fun u => σ * (mk a u) ^ 3 - 3 * exp u * (mk a u) ^ 2)
      (σ * (-(3 * exp u * (mk a u) ^ 2))
        - (3 * exp u * (mk a u) ^ 2 + 3 * exp u * (-(2 * exp u * mk a u)))) u := by
    have ha := (hasDerivAt_mk3 a u).const_mul σ
    have hb := ((hasDerivAt_exp u).const_mul 3).mul (hasDerivAt_mk2 a u)
    exact ha.sub hb
  have h := (hasDerivAt_exp_mul σ u).mul h1
  refine HasDerivAt.congr_deriv (f := Bk1 σ a) h ?_
  unfold Bk2; ring

theorem continuous_mk (a : ℝ) : Continuous (mk a) := by
  unfold mk; fun_prop

theorem continuous_Bk (σ a : ℝ) : Continuous (Bk σ a) := by
  have := continuous_mk a
  unfold Bk; fun_prop

theorem continuous_Bk1 (σ a : ℝ) : Continuous (Bk1 σ a) := by
  have := continuous_mk a
  unfold Bk1; fun_prop

theorem continuous_Bk2 (σ a : ℝ) : Continuous (Bk2 σ a) := by
  have := continuous_mk a
  unfold Bk2; fun_prop

/-! ### Integrability -/

lemma integrable_of_exp_bound {f : ℝ → ℝ} (hf : Continuous f) {L C σ : ℝ} (hσ : 0 < σ)
    (h1 : ∀ u ≤ L, |f u| ≤ C * exp (σ * u)) (h2 : ∀ u, L < u → f u = 0) : Integrable f := by
  have hI : IntegrableOn f (Iic L) := by
    have hexp : IntegrableOn (fun u => C * exp (σ * u)) (Iic L) :=
      (integrableOn_exp_mul_Iic hσ L).const_mul C
    refine Integrable.mono' hexp hf.aestronglyMeasurable.restrict ?_
    refine (ae_restrict_mem measurableSet_Iic).mono fun u hu => ?_
    rw [Real.norm_eq_abs]
    exact h1 u hu
  have hJ : IntegrableOn f (Ioi L) := by
    have h0 : IntegrableOn (fun _ : ℝ => (0 : ℝ)) (Ioi L) := integrableOn_zero
    exact h0.congr_fun (fun u hu => (h2 u hu).symm) measurableSet_Ioi
  have := hI.union hJ
  rwa [Iic_union_Ioi, integrableOn_univ] at this

lemma abs_Bk_le {a : ℝ} (ha : 0 ≤ a) (σ u : ℝ) : |Bk σ a u| ≤ a ^ 3 * exp (σ * u) := by
  unfold Bk
  rw [abs_mul, abs_of_pos (exp_pos _), abs_of_nonneg (pow_nonneg (mk_nonneg a u) 3), mul_comm]
  exact mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (mk_nonneg a u) (mk_le ha u) 3)
    (exp_pos _).le

lemma Bk_of_log_lt {a u : ℝ} (ha : 0 < a) (σ : ℝ) (h : log a < u) : Bk σ a u = 0 := by
  unfold Bk; rw [mk_of_log_lt ha h]; ring

lemma Bk1_of_log_lt {a u : ℝ} (ha : 0 < a) (σ : ℝ) (h : log a < u) : Bk1 σ a u = 0 := by
  unfold Bk1; rw [mk_of_log_lt ha h]; ring

lemma Bk2_of_log_lt {a u : ℝ} (ha : 0 < a) (σ : ℝ) (h : log a < u) : Bk2 σ a u = 0 := by
  unfold Bk2; rw [mk_of_log_lt ha h]; ring

theorem integrable_Bk {σ a : ℝ} (hσ : 0 < σ) (ha : 0 < a) : Integrable (Bk σ a) :=
  integrable_of_exp_bound (continuous_Bk σ a) hσ (L := log a) (C := a ^ 3)
    (fun u _ => abs_Bk_le ha.le σ u) (fun _ hu => Bk_of_log_lt ha σ hu)

theorem integrable_Bk1 {σ a : ℝ} (hσ : 0 < σ) (ha : 0 < a) : Integrable (Bk1 σ a) := by
  refine integrable_of_exp_bound (continuous_Bk1 σ a) hσ (L := log a)
    (C := σ * a ^ 3 + 3 * a * a ^ 2) (fun u hu => ?_) (fun u hu => Bk1_of_log_lt ha σ hu)
  have he : exp u ≤ a := exp_le_of_le_log ha hu
  have hm0 := mk_nonneg a u
  have hm := mk_le ha.le u
  have h3 : (mk a u) ^ 3 ≤ a ^ 3 := pow_le_pow_left₀ hm0 hm 3
  have h2 : (mk a u) ^ 2 ≤ a ^ 2 := pow_le_pow_left₀ hm0 hm 2
  unfold Bk1
  rw [abs_mul, abs_of_pos (exp_pos _), mul_comm]
  refine mul_le_mul_of_nonneg_right ?_ (exp_pos _).le
  have t1 : |σ * (mk a u) ^ 3| ≤ σ * a ^ 3 := by
    rw [abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left h3 hσ.le
  have t2 : |3 * exp u * (mk a u) ^ 2| ≤ 3 * a * a ^ 2 := by
    rw [abs_of_nonneg (by positivity)]
    exact mul_le_mul (by linarith) h2 (by positivity) (by positivity)
  exact (abs_sub _ _).trans (add_le_add t1 t2)

theorem integrable_Bk2 {σ a : ℝ} (hσ : 0 < σ) (ha : 0 < a) : Integrable (Bk2 σ a) := by
  refine integrable_of_exp_bound (continuous_Bk2 σ a) hσ (L := log a)
    (C := σ ^ 2 * a ^ 3 + 3 * (2 * σ + 1) * a * a ^ 2 + 6 * a ^ 2 * a) (fun u hu => ?_)
    (fun u hu => Bk2_of_log_lt ha σ hu)
  have he : exp u ≤ a := exp_le_of_le_log ha hu
  have he0 := (exp_pos u).le
  have hm0 := mk_nonneg a u
  have hm := mk_le ha.le u
  have h3 : (mk a u) ^ 3 ≤ a ^ 3 := pow_le_pow_left₀ hm0 hm 3
  have h2 : (mk a u) ^ 2 ≤ a ^ 2 := pow_le_pow_left₀ hm0 hm 2
  have he2 : exp u ^ 2 ≤ a ^ 2 := pow_le_pow_left₀ he0 he 2
  unfold Bk2
  rw [abs_mul, abs_of_pos (exp_pos _), mul_comm]
  refine mul_le_mul_of_nonneg_right ?_ (exp_pos _).le
  have t1 : |σ ^ 2 * (mk a u) ^ 3| ≤ σ ^ 2 * a ^ 3 := by
    rw [abs_of_nonneg (by positivity)]
    exact mul_le_mul_of_nonneg_left h3 (by positivity)
  have t2 : |3 * (2 * σ + 1) * exp u * (mk a u) ^ 2| ≤ 3 * (2 * σ + 1) * a * a ^ 2 := by
    rw [abs_of_nonneg (by positivity)]
    exact mul_le_mul (mul_le_mul_of_nonneg_left he (by positivity)) h2 (by positivity)
      (by positivity)
  have t3 : |6 * exp u ^ 2 * mk a u| ≤ 6 * a ^ 2 * a := by
    rw [abs_of_nonneg (by positivity)]
    exact mul_le_mul (mul_le_mul_of_nonneg_left he2 (by norm_num)) hm hm0 (by positivity)
  calc |σ ^ 2 * (mk a u) ^ 3 - 3 * (2 * σ + 1) * exp u * (mk a u) ^ 2 + 6 * exp u ^ 2 * mk a u|
      ≤ |σ ^ 2 * (mk a u) ^ 3 - 3 * (2 * σ + 1) * exp u * (mk a u) ^ 2|
        + |6 * exp u ^ 2 * mk a u| := abs_add_le _ _
    _ ≤ (|σ ^ 2 * (mk a u) ^ 3| + |3 * (2 * σ + 1) * exp u * (mk a u) ^ 2|)
        + |6 * exp u ^ 2 * mk a u| := by gcongr; exact abs_sub _ _
    _ ≤ _ := by linarith

/-! ### The Fourier transform of the building block -/

open Complex in
/-- `∫ Bk σ a u e^{iru} du = 6 a^{σ+3+ir} / ((σ+ir)(σ+1+ir)(σ+2+ir)(σ+3+ir))`. -/
theorem integral_Bk_mul_cexp {σ a : ℝ} (hσ : 0 < σ) (ha : 0 < a) (r : ℝ) :
    ∫ u : ℝ, ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u)
      = 6 * (a : ℂ) ^ ((σ : ℂ) + 3 + I * r)
        / (((σ : ℂ) + I * r) * ((σ : ℂ) + 1 + I * r) * ((σ : ℂ) + 2 + I * r)
          * ((σ : ℂ) + 3 + I * r)) := by
  set s : ℂ := (σ : ℂ) + I * r with hs
  have hsre : s.re = σ := by simp [hs]
  -- integrability of the pieces on `Iic (Real.log a)`
  have hint : ∀ j : ℕ, IntegrableOn (fun u : ℝ => cexp ((s + j) * u)) (Iic (Real.log a)) := fun j =>
    integrableOn_exp_mul_complex_Iic (by simp [hsre]; positivity) _
  have hval : ∀ j : ℕ, ∫ u : ℝ in Iic (Real.log a), cexp ((s + j) * u)
      = (a : ℂ) ^ (s + j) / (s + j) := by
    intro j
    rw [integral_exp_mul_complex_Iic (by simp [hsre]; positivity)]
    congr 1
    rw [Complex.cpow_def_of_ne_zero (by exact_mod_cast ha.ne'), Complex.ofReal_log ha.le, mul_comm]
  -- the integrand vanishes on `Ioi (Real.log a)`
  have hzero : ∫ u : ℝ in Ioi (Real.log a), ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u) = 0 := by
    rw [setIntegral_congr_fun (g := fun _ => (0 : ℂ)) measurableSet_Ioi]
    · simp
    · intro u hu
      show ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u) = 0
      rw [Bk_of_log_lt ha σ hu]; simp
  -- on `Iic (Real.log a)` the integrand is a sum of four exponentials
  have hexpand : ∀ u ∈ Iic (Real.log a), ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u)
      = (a : ℂ) ^ 3 * cexp ((s + (0 : ℕ)) * u) - 3 * (a : ℂ) ^ 2 * cexp ((s + (1 : ℕ)) * u)
        + 3 * (a : ℂ) * cexp ((s + (2 : ℕ)) * u) - cexp ((s + (3 : ℕ)) * u) := by
    intro u hu
    have he : Real.exp u ≤ a := exp_le_of_le_log ha hu
    unfold Bk
    rw [mk_of_le he]
    have e0 : cexp ((s + (0 : ℕ)) * u) = cexp (σ * u) * cexp (I * r * u) := by
      rw [← Complex.exp_add]; congr 1; rw [hs]; push_cast; ring
    have e1 : cexp ((s + (1 : ℕ)) * u) = cexp (σ * u) * cexp (I * r * u) * cexp u := by
      rw [← Complex.exp_add, ← Complex.exp_add]; congr 1; rw [hs]; push_cast; ring
    have e2 : cexp ((s + (2 : ℕ)) * u) = cexp (σ * u) * cexp (I * r * u) * cexp u ^ 2 := by
      rw [← Complex.exp_nat_mul, ← Complex.exp_add, ← Complex.exp_add]
      congr 1; rw [hs]; push_cast; ring
    have e3 : cexp ((s + (3 : ℕ)) * u) = cexp (σ * u) * cexp (I * r * u) * cexp u ^ 3 := by
      rw [← Complex.exp_nat_mul, ← Complex.exp_add, ← Complex.exp_add]
      congr 1; rw [hs]; push_cast; ring
    rw [e0, e1, e2, e3]
    push_cast
    ring
  have hsplit : ∫ u : ℝ, ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u)
      = ∫ u : ℝ in Iic (Real.log a), ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u) := by
    have hfi : Integrable (fun u : ℝ => ((Bk σ a u : ℝ) : ℂ) * cexp (I * r * u)) := by
      have h := (integrable_Bk hσ ha).ofReal (𝕜 := ℂ)
      refine h.mul_bdd (c := 1) (by fun_prop) ?_
      refine Filter.Eventually.of_forall fun u => ?_
      rw [show I * (r : ℂ) * (u : ℂ) = ((r * u : ℝ) : ℂ) * I by push_cast; ring,
        Complex.norm_exp_ofReal_mul_I]
    rw [← intervalIntegral.integral_Iic_add_Ioi hfi.integrableOn hfi.integrableOn, hzero, add_zero]
  rw [hsplit, setIntegral_congr_fun measurableSet_Iic hexpand]
  have i0 := (hint 0).const_mul ((a : ℂ) ^ 3)
  have i1 := (hint 1).const_mul (3 * (a : ℂ) ^ 2)
  have i2 := (hint 2).const_mul (3 * (a : ℂ))
  have i3 := hint 3
  have j01 : IntegrableOn (fun u : ℝ => (a : ℂ) ^ 3 * cexp ((s + (0 : ℕ)) * u)
      - 3 * (a : ℂ) ^ 2 * cexp ((s + (1 : ℕ)) * u)) (Iic (Real.log a)) := i0.sub i1
  have j012 : IntegrableOn (fun u : ℝ => (a : ℂ) ^ 3 * cexp ((s + (0 : ℕ)) * u)
      - 3 * (a : ℂ) ^ 2 * cexp ((s + (1 : ℕ)) * u) + 3 * (a : ℂ) * cexp ((s + (2 : ℕ)) * u))
      (Iic (Real.log a)) := j01.add i2
  rw [integral_sub j012 i3, integral_add j01 i2, integral_sub i0 i1,
    integral_const_mul, integral_const_mul, integral_const_mul, hval 0, hval 1, hval 2, hval 3]
  have ha0 : (a : ℂ) ≠ 0 := by exact_mod_cast ha.ne'
  have hne : ∀ j : ℕ, s + j ≠ 0 := by
    intro j h
    have := congrArg Complex.re h
    simp [hsre] at this
    have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    linarith
  have p1 : (a : ℂ) ^ (s + (1 : ℕ)) = (a : ℂ) ^ s * a := by
    rw [Complex.cpow_add _ _ ha0]; simp
  have p2 : (a : ℂ) ^ (s + (2 : ℕ)) = (a : ℂ) ^ s * a ^ 2 := by
    rw [Complex.cpow_add _ _ ha0]; simp
  have p3 : (a : ℂ) ^ (s + (3 : ℕ)) = (a : ℂ) ^ s * a ^ 3 := by
    rw [Complex.cpow_add _ _ ha0]; simp
  have p0 : (a : ℂ) ^ (s + (0 : ℕ)) = (a : ℂ) ^ s := by simp
  have p3' : (a : ℂ) ^ (s + 3) = (a : ℂ) ^ s * a ^ 3 := by
    have := p3; push_cast at this; exact this
  have h0 := hne 0
  have h1 := hne 1
  have h2 := hne 2
  have h3 := hne 3
  push_cast at h0 h1 h2 h3
  rw [p0, p1, p2, p3, show (σ : ℂ) + 3 + I * r = s + 3 by rw [hs]; ring,
    show (σ : ℂ) + 1 + I * r = s + 1 by rw [hs]; ring,
    show (σ : ℂ) + 2 + I * r = s + 2 by rw [hs]; ring, p3']
  push_cast
  rw [add_zero] at h0 ⊢
  field_simp
  ring

/-! ### A second-difference bound -/

/-- If `f'' ` is bounded by `C` on `[t−1, t+1]` then `‖f(t+1) − 2 f(t) + f(t−1)‖ ≤ C`. -/
theorem norm_second_diff_le {f f' f'' : ℝ → ℂ} {t C : ℝ}
    (hf : ∀ x ∈ Icc (t - 1) (t + 1), HasDerivAt f (f' x) x)
    (hf' : ∀ x ∈ Icc (t - 1) (t + 1), HasDerivAt f' (f'' x) x)
    (hC : ∀ x ∈ Icc (t - 1) (t + 1), ‖f'' x‖ ≤ C) :
    ‖f (t + 1) - 2 * f t + f (t - 1)‖ ≤ C := by
  have hmem : ∀ σ ∈ Icc (0 : ℝ) 1, t + σ ∈ Icc (t - 1) (t + 1) ∧ t - σ ∈ Icc (t - 1) (t + 1) := by
    intro σ hσ
    exact ⟨⟨by linarith [hσ.1], by linarith [hσ.2]⟩, ⟨by linarith [hσ.2], by linarith [hσ.1]⟩⟩
  -- `g σ = f (t+σ) + f (t−σ)`
  have hg : ∀ σ ∈ uIcc (0 : ℝ) 1, HasDerivAt (fun σ => f (t + σ) + f (t - σ))
      (f' (t + σ) - f' (t - σ)) σ := by
    intro σ hσ
    rw [uIcc_of_le zero_le_one] at hσ
    obtain ⟨h1, h2⟩ := hmem σ hσ
    have a1 : HasDerivAt (fun σ : ℝ => t + σ) 1 σ := by simpa using (hasDerivAt_id σ).const_add t
    have a2 : HasDerivAt (fun σ : ℝ => t - σ) (-1) σ := by simpa using (hasDerivAt_id σ).const_sub t
    have b1 := (hf _ h1).scomp σ a1
    have b2 := (hf _ h2).scomp σ a2
    have := b1.add b2
    refine this.congr_deriv ?_
    simp
    ring
  have hcont : ContinuousOn (fun σ => f' (t + σ) - f' (t - σ)) (uIcc (0 : ℝ) 1) := by
    rw [uIcc_of_le zero_le_one]
    intro σ hσ
    obtain ⟨h1, h2⟩ := hmem σ hσ
    have c1 : ContinuousAt (fun σ : ℝ => f' (t + σ)) σ :=
      ContinuousAt.comp (f := fun σ : ℝ => t + σ) (hf' _ h1).continuousAt (by fun_prop)
    have c2 : ContinuousAt (fun σ : ℝ => f' (t - σ)) σ :=
      ContinuousAt.comp (f := fun σ : ℝ => t - σ) (hf' _ h2).continuousAt (by fun_prop)
    exact (c1.sub c2).continuousWithinAt
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt hg hcont.intervalIntegrable
  have hbound : ∀ σ ∈ Icc (0 : ℝ) 1, ‖f' (t + σ) - f' (t - σ)‖ ≤ 2 * C * σ := by
    intro σ hσ
    obtain ⟨h1, h2⟩ := hmem σ hσ
    have := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le (f := f') (f' := f'')
      (s := Icc (t - 1) (t + 1)) (fun x hx => (hf' x hx).hasDerivWithinAt) hC (convex_Icc _ _) h2 h1
    rw [show t + σ - (t - σ) = 2 * σ by ring, Real.norm_eq_abs,
      abs_of_nonneg (by linarith [hσ.1])] at this
    linarith
  have e : f (t + 1) - 2 * f t + f (t - 1)
      = (fun σ => f (t + σ) + f (t - σ)) 1 - (fun σ => f (t + σ) + f (t - σ)) 0 := by
    simp; ring
  rw [e, ← hFTC]
  calc ‖∫ σ in (0 : ℝ)..1, (f' (t + σ) - f' (t - σ))‖
      ≤ ∫ σ in (0 : ℝ)..1, ‖f' (t + σ) - f' (t - σ)‖ :=
        intervalIntegral.norm_integral_le_integral_norm zero_le_one
    _ ≤ ∫ σ in (0 : ℝ)..1, 2 * C * σ := by
        refine intervalIntegral.integral_mono_on zero_le_one ?_ ?_ hbound
        · exact (hcont.norm).intervalIntegrable
        · exact (by fun_prop : Continuous fun σ : ℝ => 2 * C * σ).intervalIntegrable _ _
    _ = C := by
        rw [intervalIntegral.integral_const_mul, integral_id]
        ring

end Allx
