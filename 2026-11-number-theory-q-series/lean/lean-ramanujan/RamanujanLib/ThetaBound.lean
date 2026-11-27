/-
RamanujanLib/ThetaBound.lean — `|θ(x) − x| ≤ K · x / log⁵ x` for all `x ≥ 2`.

Input: `MediumPNT` (zeta23's port of PrimeNumberTheoremAnd),
`ψ(x) − x = O(x · exp(−c (log x)^{1/10}))`, and Mathlib's `ψ − θ ≤ C √x`.
-/
import Mathlib
import Zeta23.FromPNTPlus.MediumPNT

open Real Filter Topology Asymptotics

noncomputable section

namespace RamanujanLib

/-- Eventually `exp(−c (log x)^{1/10}) ≤ c^{-50} / log⁵ x`. -/
theorem eventually_exp_le (c : ℝ) (hc : 0 < c) :
    ∀ᶠ x : ℝ in atTop,
      exp (-c * (log x) ^ ((1 : ℝ) / 10)) ≤ (1 / c ^ 50) * (1 / log x ^ 5) := by
  have h1 : Tendsto (fun x : ℝ => c * (log x) ^ ((1 : ℝ) / 10)) atTop atTop :=
    ((tendsto_rpow_atTop (by norm_num)).comp tendsto_log_atTop).const_mul_atTop hc
  have h2 := (Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 50).comp h1
  have h3 : ∀ᶠ x : ℝ in atTop,
      (c * (log x) ^ ((1 : ℝ) / 10)) ^ 50 * exp (-(c * (log x) ^ ((1 : ℝ) / 10))) ≤ 1 :=
    (h2.eventually (ge_mem_nhds (zero_lt_one' ℝ))).mono fun x hx => hx
  filter_upwards [h3, eventually_gt_atTop 1] with x hx hx1
  have hl : 0 < log x := log_pos hx1
  have hpow : ((log x) ^ ((1 : ℝ) / 10)) ^ 50 = log x ^ 5 := by
    rw [← rpow_natCast, ← rpow_mul hl.le]
    norm_num
  rw [mul_pow, hpow, neg_mul_eq_neg_mul] at hx
  have hpos : 0 < c ^ 50 * log x ^ 5 := by positivity
  calc exp (-c * (log x) ^ ((1 : ℝ) / 10))
      = (c ^ 50 * log x ^ 5 * exp (-c * (log x) ^ ((1 : ℝ) / 10))) / (c ^ 50 * log x ^ 5) := by
        field_simp
    _ ≤ 1 / (c ^ 50 * log x ^ 5) := by gcongr
    _ = (1 / c ^ 50) * (1 / log x ^ 5) := by rw [one_div_mul_one_div]

/-- Eventually `√x ≤ x / log⁵ x`. -/
theorem eventually_sqrt_le : ∀ᶠ x : ℝ in atTop, √x ≤ x / log x ^ 5 := by
  have h := (isLittleO_log_rpow_rpow_atTop (5 : ℝ) (show (0 : ℝ) < 1 / 2 by norm_num)).bound
    (zero_lt_one' ℝ)
  filter_upwards [h, eventually_gt_atTop 1] with x hx hx1
  have hl : 0 < log x := log_pos hx1
  have hx0 : 0 < x := by linarith
  rw [Real.norm_of_nonneg (rpow_nonneg hl.le _), Real.norm_of_nonneg (rpow_nonneg hx0.le _),
    one_mul, ← sqrt_eq_rpow] at hx
  have h5 : log x ^ (5 : ℝ) = log x ^ 5 := by exact_mod_cast rpow_natCast (log x) 5
  rw [h5] at hx
  rw [le_div_iff₀ (pow_pos hl 5)]
  calc √x * log x ^ 5 ≤ √x * √x := by gcongr
    _ = x := mul_self_sqrt hx0.le

/-- `|θ(x) − x| ≤ K · x / log⁵ x` for all sufficiently large `x`. -/
theorem eventually_theta_sub_le :
    ∃ K : ℝ, ∀ᶠ x : ℝ in atTop, |Chebyshev.theta x - x| ≤ K * (x / log x ^ 5) := by
  obtain ⟨c, hc, hO⟩ := MediumPNT
  obtain ⟨K₁, hK₁⟩ := hO.bound
  obtain ⟨C₂, hC₂⟩ := Chebyshev.psi_sub_theta_le_mul_sqrt
  refine ⟨|K₁| * (1 / c ^ 50) + |C₂|, ?_⟩
  filter_upwards [hK₁, eventually_exp_le c hc, eventually_sqrt_le, eventually_gt_atTop 1]
    with x h1 h2 h3 hx1
  have hl : 0 < log x := log_pos hx1
  have hx0 : 0 < x := by linarith
  simp only [Pi.sub_apply, id, Real.norm_eq_abs] at h1
  rw [abs_mul, abs_of_pos hx0, abs_of_pos (exp_pos _)] at h1
  have hq : 0 ≤ x / log x ^ 5 := by positivity
  have hψ : |Chebyshev.psi x - x| ≤ |K₁| * (1 / c ^ 50) * (x / log x ^ 5) := by
    calc |Chebyshev.psi x - x| ≤ K₁ * (x * exp (-c * (log x) ^ ((1 : ℝ) / 10))) := h1
      _ ≤ |K₁| * (x * ((1 / c ^ 50) * (1 / log x ^ 5))) := by
          gcongr
          exact le_abs_self _
      _ = |K₁| * (1 / c ^ 50) * (x / log x ^ 5) := by ring
  have hd1 : 0 ≤ Chebyshev.psi x - Chebyshev.theta x := sub_nonneg.mpr (Chebyshev.theta_le_psi x)
  have hd2 : Chebyshev.psi x - Chebyshev.theta x ≤ |C₂| * (x / log x ^ 5) := by
    calc Chebyshev.psi x - Chebyshev.theta x ≤ C₂ * √x := hC₂ x
      _ ≤ |C₂| * (x / log x ^ 5) := by
          gcongr
          exact le_abs_self _
  have key : Chebyshev.theta x - x
      = (Chebyshev.psi x - x) - (Chebyshev.psi x - Chebyshev.theta x) := by ring
  rw [key, add_mul]
  calc |Chebyshev.psi x - x - (Chebyshev.psi x - Chebyshev.theta x)|
      ≤ |Chebyshev.psi x - x| + |Chebyshev.psi x - Chebyshev.theta x| := abs_sub _ _
    _ ≤ _ := by
        rw [abs_of_nonneg hd1]
        linarith

/-- **Prime number theorem with a `log⁻⁵` saving, on all of `[2, ∞)`**:
`|θ(x) − x| ≤ K · x / log⁵ x` for every `x ≥ 2`. -/
theorem exists_theta_sub_le :
    ∃ K : ℝ, 0 ≤ K ∧ ∀ x : ℝ, 2 ≤ x → |Chebyshev.theta x - x| ≤ K * (x / log x ^ 5) := by
  obtain ⟨K, hK⟩ := eventually_theta_sub_le
  obtain ⟨T₀, hT₀⟩ := eventually_atTop.mp hK
  set T : ℝ := max 2 T₀ with hT
  have hT2 : (2 : ℝ) ≤ T := le_max_left _ _
  have hlT : 0 < log T := log_pos (by linarith)
  refine ⟨max |K| (3 * log T ^ 5), le_max_of_le_left (abs_nonneg _), fun x hx => ?_⟩
  have hl : 0 < log x := log_pos (by linarith)
  have hq : 0 ≤ x / log x ^ 5 := by positivity
  rcases le_total T x with hxT | hxT
  · calc |Chebyshev.theta x - x| ≤ K * (x / log x ^ 5) := hT₀ x ((le_max_right _ _).trans hxT)
      _ ≤ |K| * (x / log x ^ 5) := by gcongr; exact le_abs_self _
      _ ≤ _ := by gcongr; exact le_max_left _ _
  · have h4 : Chebyshev.theta x ≤ log 4 * x := Chebyshev.theta_le_log4_mul_x (by linarith)
    have hlog4 : log 4 ≤ 2 := by
      have : log 4 = 2 * log 2 := by
        rw [show (4 : ℝ) = 2 ^ 2 by norm_num, log_pow]; norm_num
      have := Real.log_two_lt_d9
      linarith
    have hθ := Chebyshev.theta_nonneg x
    have habs : |Chebyshev.theta x - x| ≤ 3 * x := by
      rw [abs_le]
      constructor <;> nlinarith
    have hlx : log x ≤ log T := log_le_log (by linarith) hxT
    have : 3 * x ≤ 3 * log T ^ 5 * (x / log x ^ 5) := by
      rw [mul_div_assoc', le_div_iff₀ (pow_pos hl 5)]
      have : log x ^ 5 ≤ log T ^ 5 := pow_le_pow_left₀ hl.le hlx 5
      nlinarith
    calc |Chebyshev.theta x - x| ≤ 3 * x := habs
      _ ≤ 3 * log T ^ 5 * (x / log x ^ 5) := this
      _ ≤ _ := by gcongr; exact le_max_right _ _

end RamanujanLib
