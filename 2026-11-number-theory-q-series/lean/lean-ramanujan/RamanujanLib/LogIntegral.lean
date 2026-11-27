/-
RamanujanLib/LogIntegral.lean — the integrals `J k x = ∫_2^x dt / log^k t` (Mathlib only).

* integration by parts: `J k x = x / log^k x − 2 / log^k 2 + k · J (k+1) x`;
* the bound `J k x ≤ C_k · x / log^k x` for all `x ≥ 2`;
* the five-term expansion
  `x / log x + J 2 x = x (1/L + 1/L² + 2/L³ + 6/L⁴ + 24/L⁵) + O(x / L⁶)`, `L = log x`.
-/
import Mathlib

open Real Set MeasureTheory

noncomputable section

namespace RamanujanLib

/-- `J k x = ∫_2^x dt / log^k t`. -/
def J (k : ℕ) (x : ℝ) : ℝ := ∫ t in (2 : ℝ)..x, 1 / log t ^ k

/-- Ramanujan's main term as a polynomial in `u = 1 / log x`:
`P u = u + u² + 2u³ + 6u⁴ + 24u⁵ = Σ_{k=0}^{4} k! u^{k+1}`. -/
def P (u : ℝ) : ℝ := u + u ^ 2 + 2 * u ^ 3 + 6 * u ^ 4 + 24 * u ^ 5

theorem log_pos_of_two_le {t : ℝ} (ht : 2 ≤ t) : 0 < log t := log_pos (by linarith)

theorem log_two_le_log {t : ℝ} (ht : 2 ≤ t) : log 2 ≤ log t := log_le_log (by norm_num) ht

theorem hasDerivAt_div_log_pow (k : ℕ) {t : ℝ} (ht : 1 < t) :
    HasDerivAt (fun t : ℝ => t / log t ^ k) (1 / log t ^ k - k / log t ^ (k + 1)) t := by
  have ht0 : t ≠ 0 := by linarith
  have hl : log t ≠ 0 := (log_pos ht).ne'
  have h : HasDerivAt (fun t : ℝ => t / log t ^ k)
      ((1 * log t ^ k - t * (k * log t ^ (k - 1) * t⁻¹)) / (log t ^ k) ^ 2) t :=
    (hasDerivAt_id' t).fun_div ((hasDerivAt_log ht0).fun_pow k) (pow_ne_zero k hl)
  convert h using 1
  rcases k with _ | k
  · simp
  · simp only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one]
    field_simp
    ring

theorem continuousOn_inv_log_pow (k : ℕ) {a b : ℝ} (ha : 1 < a) (hb : 1 < b) :
    ContinuousOn (fun t : ℝ => 1 / log t ^ k) (uIcc a b) := by
  intro t ht
  rw [mem_uIcc] at ht
  have h1 : 1 < t := by rcases ht with h | h <;> linarith [h.1]
  have ht0 : t ≠ 0 := by linarith
  have : log t ^ k ≠ 0 := pow_ne_zero _ (log_pos h1).ne'
  exact ContinuousAt.continuousWithinAt (by fun_prop)

theorem intervalIntegrable_inv_log_pow (k : ℕ) {a b : ℝ} (ha : 1 < a) (hb : 1 < b) :
    IntervalIntegrable (fun t : ℝ => 1 / log t ^ k) volume a b :=
  (continuousOn_inv_log_pow k ha hb).intervalIntegrable

/-- Integration by parts on `[a, b]`, `1 < a`, `1 < b`. -/
theorem integral_inv_log_pow_sub (k : ℕ) {a b : ℝ} (ha : 1 < a) (hb : 1 < b) :
    ∫ t in a..b, (1 / log t ^ k - k / log t ^ (k + 1)) = b / log b ^ k - a / log a ^ k := by
  have hint : IntervalIntegrable (fun t : ℝ => 1 / log t ^ k - k / log t ^ (k + 1)) volume a b := by
    have h2 : IntervalIntegrable (fun t : ℝ => (k : ℝ) / log t ^ (k + 1)) volume a b := by
      simpa [div_eq_mul_inv] using (intervalIntegrable_inv_log_pow (k + 1) ha hb).const_mul (k : ℝ)
    exact (intervalIntegrable_inv_log_pow k ha hb).sub h2
  refine intervalIntegral.integral_eq_sub_of_hasDerivAt (f := fun t : ℝ => t / log t ^ k) ?_ hint
  intro t ht
  rw [mem_uIcc] at ht
  exact hasDerivAt_div_log_pow k (by rcases ht with h | h <;> linarith [h.1])

/-- Integration by parts on `[a, b]`:
`∫_a^b dt/log^k t = b / log^k b − a / log^k a + k ∫_a^b dt/log^{k+1} t`. -/
theorem integral_inv_log_pow_ibp (k : ℕ) {a b : ℝ} (ha : 1 < a) (hb : 1 < b) :
    ∫ t in a..b, 1 / log t ^ k
      = b / log b ^ k - a / log a ^ k + k * ∫ t in a..b, 1 / log t ^ (k + 1) := by
  have h := integral_inv_log_pow_sub k ha hb
  have h2 : IntervalIntegrable (fun t : ℝ => (k : ℝ) / log t ^ (k + 1)) volume a b := by
    simpa [div_eq_mul_inv] using (intervalIntegrable_inv_log_pow (k + 1) ha hb).const_mul (k : ℝ)
  rw [intervalIntegral.integral_sub (intervalIntegrable_inv_log_pow k ha hb) h2] at h
  have h3 : ∫ t in a..b, (k : ℝ) / log t ^ (k + 1) = k * ∫ t in a..b, 1 / log t ^ (k + 1) := by
    rw [← intervalIntegral.integral_const_mul]
    exact intervalIntegral.integral_congr fun _ _ => by ring
  rw [h3] at h
  linarith

/-- Integration by parts: `J k x = x / log^k x − 2 / log^k 2 + k · J (k+1) x`. -/
theorem J_ibp (k : ℕ) {x : ℝ} (hx : 2 ≤ x) :
    J k x = x / log x ^ k - 2 / log 2 ^ k + k * J (k + 1) x :=
  integral_inv_log_pow_ibp k (by norm_num) (by linarith)

theorem J_nonneg (k : ℕ) {x : ℝ} (hx : 2 ≤ x) : 0 ≤ J k x :=
  intervalIntegral.integral_nonneg hx fun _ ht =>
    div_nonneg zero_le_one (pow_nonneg (log_pos_of_two_le ht.1).le _)

/-- `x / log^k x` is bounded below by a positive constant on `[2, ∞)`. -/
theorem exists_pos_le_div_log_pow (k : ℕ) : ∃ c : ℝ, 0 < c ∧ ∀ x : ℝ, 2 ≤ x → c ≤ x / log x ^ k := by
  refine ⟨1 / (k.factorial : ℝ), by positivity, fun x hx => ?_⟩
  have hl : 0 < log x := log_pos_of_two_le hx
  have h := Real.pow_div_factorial_le_exp (log x) hl.le k
  rw [exp_log (by linarith)] at h
  have hf : (0 : ℝ) < k.factorial := by positivity
  rw [div_le_iff₀ hf] at h
  rw [le_div_iff₀ (pow_pos hl k), div_mul_eq_mul_div, div_le_iff₀ hf]
  linarith

/-- `∫_2^x dt / log^k t ≤ C · x / log^k x` for all `x ≥ 2`. -/
theorem exists_J_le (k : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ x : ℝ, 2 ≤ x → J k x ≤ C * (x / log x ^ k) := by
  obtain ⟨c, hc, hcx⟩ := exists_pos_le_div_log_pow k
  -- split point: beyond `T`, `log t ≥ 2k`
  set T : ℝ := max 2 (exp (2 * k)) with hT
  have hT2 : (2 : ℝ) ≤ T := le_max_left _ _
  have hTe : exp (2 * k) ≤ T := le_max_right _ _
  have hlogT : ∀ t : ℝ, T ≤ t → 2 * (k : ℝ) ≤ log t := fun t ht => by
    have := log_le_log (exp_pos _) (hTe.trans ht)
    rwa [log_exp] at this
  -- the part up to `T` is a constant `A`
  set A : ℝ := J k T with hA
  have hA0 : 0 ≤ A := J_nonneg k hT2
  -- the part beyond `T`
  have htail : ∀ x : ℝ, T ≤ x → ∫ t in T..x, 1 / log t ^ k ≤ 2 * (x / log x ^ k) := by
    intro x hx
    have hT1 : (1 : ℝ) < T := by linarith
    have hx1 : (1 : ℝ) < x := by linarith
    have hint1 := intervalIntegrable_inv_log_pow k hT1 hx1
    have hint2 : IntervalIntegrable (fun t : ℝ => (k : ℝ) / log t ^ (k + 1)) volume T x := by
      simpa [div_eq_mul_inv] using (intervalIntegrable_inv_log_pow (k + 1) hT1 hx1).const_mul (k : ℝ)
    have hmono : ∫ t in T..x, 1 / log t ^ k
        ≤ ∫ t in T..x, 2 * (1 / log t ^ k - k / log t ^ (k + 1)) := by
      refine intervalIntegral.integral_mono_on hx hint1 ((hint1.sub hint2).const_mul 2) ?_
      intro t ht
      have hl : 0 < log t := log_pos (by linarith [ht.1])
      have h2k := hlogT t ht.1
      have hk : (k : ℝ) / log t ^ (k + 1) ≤ (1 / 2) * (1 / log t ^ k) := by
        rw [pow_succ, div_le_iff₀ (by positivity)]
        have : (1 / 2) * (1 / log t ^ k) * (log t ^ k * log t) = log t / 2 := by
          field_simp
        rw [this]
        linarith
      linarith
    rw [intervalIntegral.integral_const_mul, integral_inv_log_pow_sub k hT1 hx1] at hmono
    have : 0 ≤ T / log T ^ k := div_nonneg (by linarith) (pow_nonneg (log_pos hT1).le _)
    linarith
  refine ⟨2 + A / c, by positivity, fun x hx => ?_⟩
  have hcx' := hcx x hx
  have hAc : A ≤ A / c * (x / log x ^ k) := by
    calc A = A / c * c := by field_simp
      _ ≤ A / c * (x / log x ^ k) := by gcongr
  rcases le_total x T with hxT | hxT
  · -- `x ≤ T`: the integral is at most `A`
    have hle : J k x ≤ A := by
      have hsplit := intervalIntegral.integral_add_adjacent_intervals
        (intervalIntegrable_inv_log_pow k (a := 2) (b := x) (by norm_num) (by linarith))
        (intervalIntegrable_inv_log_pow k (a := x) (b := T) (by linarith) (by linarith))
      have hnn : 0 ≤ ∫ t in x..T, 1 / log t ^ k :=
        intervalIntegral.integral_nonneg hxT fun t ht =>
          div_nonneg zero_le_one (pow_nonneg (log_pos (by linarith [ht.1])).le _)
      unfold J at hA ⊢
      linarith
    have : 0 ≤ 2 * (x / log x ^ k) := by linarith
    linarith
  · have hsplit := intervalIntegral.integral_add_adjacent_intervals
      (intervalIntegrable_inv_log_pow k (a := 2) (b := T) (by norm_num) (by linarith))
      (intervalIntegrable_inv_log_pow k (a := T) (b := x) (by linarith) (by linarith))
    have := htail x hxT
    unfold J at hA ⊢
    linarith

/-- Five-term expansion of `x / log x + ∫_2^x dt / log² t` (which is `li x − li 2 + 2 / log 2`):
it is `x · P(1 / log x) + O(x / log⁶ x)` on `[2, ∞)`. -/
theorem exists_li_expansion :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x →
      |x / log x + J 2 x - x * P (log x)⁻¹| ≤ C * (x / log x ^ 6) := by
  obtain ⟨C6, hC6, hJ6⟩ := exists_J_le 6
  obtain ⟨c, hc, hcx⟩ := exists_pos_le_div_log_pow 6
  set c0 : ℝ := 2 / log 2 ^ 2 + 2 * (2 / log 2 ^ 3) + 6 * (2 / log 2 ^ 4) + 24 * (2 / log 2 ^ 5)
    with hc0
  have hc0nn : 0 ≤ c0 := by
    have : 0 < log (2 : ℝ) := log_pos (by norm_num)
    positivity
  refine ⟨120 * C6 + c0 / c, fun x hx => ?_⟩
  have hl : 0 < log x := log_pos_of_two_le hx
  have h2 := J_ibp 2 hx
  have h3 := J_ibp 3 hx
  have h4 := J_ibp 4 hx
  have h5 := J_ibp 5 hx
  push_cast at h2 h3 h4 h5
  have hmain : x / log x + J 2 x - x * P (log x)⁻¹ = 120 * J 6 x - c0 := by
    have hP : x * P (log x)⁻¹
        = x / log x + x / log x ^ 2 + 2 * (x / log x ^ 3) + 6 * (x / log x ^ 4)
          + 24 * (x / log x ^ 5) := by
      unfold P
      field_simp
    rw [hP, hc0]
    linarith
  rw [hmain]
  have h6 := hJ6 x hx
  have h6nn := J_nonneg 6 hx
  have hcx' := hcx x hx
  have hc0le : c0 ≤ c0 / c * (x / log x ^ 6) := by
    calc c0 = c0 / c * c := by field_simp
      _ ≤ c0 / c * (x / log x ^ 6) := by gcongr
  have hnn : 0 ≤ C6 * (x / log x ^ 6) := by positivity
  rw [abs_le]
  constructor <;> nlinarith
end RamanujanLib
