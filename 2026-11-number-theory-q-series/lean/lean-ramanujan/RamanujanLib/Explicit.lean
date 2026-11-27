/-
RamanujanLib/Explicit.lean — an explicit criterion (Mathlib only, no prime number theorem):

  if `|π(y) − (c + ∫_2^y dt/log t)| ≤ A √y log y` for all `y ≥ e³⁹`, with `0 ≤ A ≤ 1/25` and
  `|c| ≤ 2`, then `π(x)² < (e x / log x) π(x/e)` for every real `x ≥ e⁴⁰`.

Schoenfeld's bound under the Riemann Hypothesis has this shape with `A = 1/(8π)` and
`c = li 2`.  The proof is Ramanujan's computation with every constant explicit:

* eight integrations by parts from the base point `a = e¹⁶`:
  `∫_a^x dt/log t = x P₈(u) − a P₈(1/16) + 8! ∫_a^x dt/log⁹ t`, `u = 1/log x`, and
  `0 ≤ ∫_a^x dt/log⁹ t ≤ (16/7) x u⁹`;
* the identity `u P₈(u/(1−u)) − S₈(u)² − u⁶ (1 + 14u) = u⁸ M₈(u)/(1−u)⁸` with `M₈ > 0` on
  `[0, 1/40]`, where `S₈(u) = P₈(u) + 92160 u⁹`;
* `e^{L/2} (L + 14) ≥ 0.159 L⁷` for `L ≥ 40`, which beats the total error `≤ 0.152 · x^{3/2}`.
-/
import RamanujanLib.LogIntegral

open Real Set MeasureTheory

noncomputable section

namespace RamanujanLib

/-- `P₈ u = Σ_{k=1}^{8} (k−1)! u^k`. -/
def P8 (u : ℝ) : ℝ :=
  u + u ^ 2 + 2 * u ^ 3 + 6 * u ^ 4 + 24 * u ^ 5 + 120 * u ^ 6 + 720 * u ^ 7 + 5040 * u ^ 8

/-- `S₈ u = P₈ u + 92160 u⁹` (`92160 = 8! · 16/7`). -/
def S8 (u : ℝ) : ℝ := P8 u + 92160 * u ^ 9

/-- The numerator polynomial of the explicit gap. -/
def M8 (u : ℝ) : ℝ :=
  145 + 252 * u - 137251 * u ^ 2 + 1113306 * u ^ 3 - 3729481 * u ^ 4 + 6568624 * u ^ 5
    - 8426305 * u ^ 6 - 458918 * u ^ 7 - 36887316 * u ^ 8 - 226329600 * u ^ 9
    - 4201059456 * u ^ 10 + 49039960320 * u ^ 11 - 195377904000 * u ^ 12
    + 418685414400 * u ^ 13 - 546718394880 * u ^ 14 + 450858355200 * u ^ 15
    - 230543366400 * u ^ 16 + 67018752000 * u ^ 17 - 8493465600 * u ^ 18

/-! ### Numerical facts -/

theorem exp_sixteen_le : exp 16 ≤ 8900000 := by
  have h := Real.exp_one_lt_d9
  calc exp 16 = exp 1 ^ 16 := by rw [← Real.exp_nat_mul]; norm_num
    _ ≤ (2.7182818286 : ℝ) ^ 16 := pow_le_pow_left₀ (exp_pos 1).le h.le 16
    _ ≤ 8900000 := by norm_num

theorem exp_twenty_ge : 485000000 ≤ exp 20 := by
  have h := Real.exp_one_gt_d9
  calc (485000000 : ℝ) ≤ (2.7182818283 : ℝ) ^ 20 := by norm_num
    _ ≤ exp 1 ^ 20 := pow_le_pow_left₀ (by norm_num) h.le 20
    _ = exp 20 := by rw [← Real.exp_nat_mul]; norm_num

theorem sqrt_exp_eq (t : ℝ) : √(exp t) = exp (t / 2) := by
  have h : exp t = exp (t / 2) * exp (t / 2) := by
    rw [← exp_add]
    congr 1
    ring
  rw [h]
  exact sqrt_mul_self (exp_pos _).le

/-- `e^{L/2} (L + 14) ≥ 0.159 L⁷` for `L ≥ 40`. -/
theorem growth {L : ℝ} (hL : 40 ≤ L) : 0.159 * L ^ 7 ≤ exp (L / 2) * (L + 14) := by
  have h1 : L / 40 ≤ exp (L / 40 - 1) := by
    have := Real.add_one_le_exp (L / 40 - 1)
    linarith
  have h2 : (L / 40) ^ 7 ≤ exp ((L - 40) / 2) := by
    calc (L / 40) ^ 7 ≤ exp (L / 40 - 1) ^ 7 := pow_le_pow_left₀ (by positivity) h1 7
      _ = exp (7 * (L / 40 - 1)) := by rw [← Real.exp_nat_mul]; norm_num
      _ ≤ exp ((L - 40) / 2) := by
          apply Real.exp_le_exp.mpr
          linarith
  have h3 : exp (L / 2) = exp 20 * exp ((L - 40) / 2) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [h3]
  have h4 : 485000000 * (L / 40) ^ 7 * 54 ≤ exp 20 * exp ((L - 40) / 2) * (L + 14) := by
    have hp : 0 ≤ (L / 40) ^ 7 := by positivity
    have h := mul_le_mul exp_twenty_ge h2 hp (exp_pos 20).le
    have h54 : (54 : ℝ) ≤ L + 14 := by linarith
    exact mul_le_mul h h54 (by norm_num) (by positivity)
  have h5 : 0.159 * L ^ 7 ≤ 485000000 * (L / 40) ^ 7 * 54 := by
    have e : (485000000 : ℝ) * (L / 40) ^ 7 * 54 = (485000000 * 54 / 40 ^ 7) * L ^ 7 := by ring
    rw [e]
    have hL7 : (0 : ℝ) ≤ L ^ 7 := by positivity
    exact mul_le_mul_of_nonneg_right (by norm_num) hL7
  linarith

/-! ### Polynomial facts on `[0, 1/40]` -/

theorem S8_bounds {u : ℝ} (hu0 : 0 ≤ u) (hu : u ≤ 1 / 40) : 0 ≤ S8 u ∧ S8 u ≤ 1.03 * u := by
  have e : S8 u = u * (1 + u + 2 * u ^ 2 + 6 * u ^ 3 + 24 * u ^ 4 + 120 * u ^ 5 + 720 * u ^ 6
      + 5040 * u ^ 7 + 92160 * u ^ 8) := by
    unfold S8 P8
    ring
  have hk : ∀ k : ℕ, u ^ k ≤ (1 / 40) ^ k := fun k => pow_le_pow_left₀ hu0 hu k
  have h2 := hk 2
  have h3 := hk 3
  have h4 := hk 4
  have h5 := hk 5
  have h6 := hk 6
  have h7 := hk 7
  have h8 := hk 8
  norm_num at h2 h3 h4 h5 h6 h7 h8
  have hb : 1 + u + 2 * u ^ 2 + 6 * u ^ 3 + 24 * u ^ 4 + 120 * u ^ 5 + 720 * u ^ 6
      + 5040 * u ^ 7 + 92160 * u ^ 8 ≤ 1.03 := by linarith
  have hb0 : 0 ≤ 1 + u + 2 * u ^ 2 + 6 * u ^ 3 + 24 * u ^ 4 + 120 * u ^ 5 + 720 * u ^ 6
      + 5040 * u ^ 7 + 92160 * u ^ 8 := by positivity
  rw [e]
  constructor
  · exact mul_nonneg hu0 hb0
  · calc u * (1 + u + 2 * u ^ 2 + 6 * u ^ 3 + 24 * u ^ 4 + 120 * u ^ 5 + 720 * u ^ 6
          + 5040 * u ^ 7 + 92160 * u ^ 8) ≤ u * 1.03 := mul_le_mul_of_nonneg_left hb hu0
      _ = 1.03 * u := by ring

theorem M8_pos {u : ℝ} (hu0 : 0 ≤ u) (hu : u ≤ 1 / 40) : 0 < M8 u := by
  have hk : ∀ k : ℕ, u ^ k ≤ (1 / 40) ^ k := fun k => pow_le_pow_left₀ hu0 hu k
  have hn : ∀ k : ℕ, 0 ≤ u ^ k := fun k => pow_nonneg hu0 k
  have h2 := hk 2
  have h4 := hk 4
  have h6 := hk 6
  have h7 := hk 7
  have h8 := hk 8
  have h9 := hk 9
  have h10 := hk 10
  have h12 := hk 12
  have h14 := hk 14
  have h16 := hk 16
  have h18 := hk 18
  norm_num at h2 h4 h6 h7 h8 h9 h10 h12 h14 h16 h18
  unfold M8
  linarith [hn 3, hn 5, hn 11, hn 13, hn 15, hn 17]

theorem gap_identity8 (u : ℝ) (hu : 1 - u ≠ 0) :
    u * P8 (u / (1 - u)) - S8 u ^ 2 - u ^ 6 * (1 + 14 * u) = u ^ 8 * M8 u / (1 - u) ^ 8 := by
  unfold S8 P8 M8
  field_simp
  ring

theorem gap_lower8 {u : ℝ} (hu0 : 0 < u) (hu : u ≤ 1 / 40) :
    u ^ 6 * (1 + 14 * u) ≤ u * P8 (u / (1 - u)) - S8 u ^ 2 := by
  have h1u : 0 < 1 - u := by linarith
  have hid := gap_identity8 u h1u.ne'
  have hpos : 0 ≤ u ^ 8 * M8 u / (1 - u) ^ 8 :=
    div_nonneg (mul_nonneg (pow_nonneg hu0.le 8) (M8_pos hu0.le hu).le) (pow_nonneg h1u.le 8)
  linarith

/-! ### The logarithmic integral from the base point `e¹⁶` -/

theorem intervalIntegrable_inv_log {a b : ℝ} (ha : 1 < a) (hb : 1 < b) :
    IntervalIntegrable (fun t : ℝ => 1 / log t) volume a b := by
  simpa using intervalIntegrable_inv_log_pow 1 ha hb

/-- Eight integrations by parts. -/
theorem integral_inv_log_eq_P8 {a x : ℝ} (ha : 1 < a) (hx : 1 < x) :
    ∫ t in a..x, 1 / log t
      = x * P8 (log x)⁻¹ - a * P8 (log a)⁻¹ + 40320 * ∫ t in a..x, 1 / log t ^ 9 := by
  have h1 : ∫ t in a..x, 1 / log t
      = x / log x - a / log a + ∫ t in a..x, 1 / log t ^ 2 := by
    simpa using integral_inv_log_pow_ibp 1 ha hx
  have h2 : ∫ t in a..x, 1 / log t ^ 2
      = x / log x ^ 2 - a / log a ^ 2 + 2 * ∫ t in a..x, 1 / log t ^ 3 := by
    simpa using integral_inv_log_pow_ibp 2 ha hx
  have h3 : ∫ t in a..x, 1 / log t ^ 3
      = x / log x ^ 3 - a / log a ^ 3 + 3 * ∫ t in a..x, 1 / log t ^ 4 := by
    simpa using integral_inv_log_pow_ibp 3 ha hx
  have h4 : ∫ t in a..x, 1 / log t ^ 4
      = x / log x ^ 4 - a / log a ^ 4 + 4 * ∫ t in a..x, 1 / log t ^ 5 := by
    simpa using integral_inv_log_pow_ibp 4 ha hx
  have h5 : ∫ t in a..x, 1 / log t ^ 5
      = x / log x ^ 5 - a / log a ^ 5 + 5 * ∫ t in a..x, 1 / log t ^ 6 := by
    simpa using integral_inv_log_pow_ibp 5 ha hx
  have h6 : ∫ t in a..x, 1 / log t ^ 6
      = x / log x ^ 6 - a / log a ^ 6 + 6 * ∫ t in a..x, 1 / log t ^ 7 := by
    simpa using integral_inv_log_pow_ibp 6 ha hx
  have h7 : ∫ t in a..x, 1 / log t ^ 7
      = x / log x ^ 7 - a / log a ^ 7 + 7 * ∫ t in a..x, 1 / log t ^ 8 := by
    simpa using integral_inv_log_pow_ibp 7 ha hx
  have h8 : ∫ t in a..x, 1 / log t ^ 8
      = x / log x ^ 8 - a / log a ^ 8 + 8 * ∫ t in a..x, 1 / log t ^ 9 := by
    simpa using integral_inv_log_pow_ibp 8 ha hx
  have ex : x * P8 (log x)⁻¹ = x / log x + x / log x ^ 2 + 2 * (x / log x ^ 3)
      + 6 * (x / log x ^ 4) + 24 * (x / log x ^ 5) + 120 * (x / log x ^ 6)
      + 720 * (x / log x ^ 7) + 5040 * (x / log x ^ 8) := by
    unfold P8
    ring
  have ea : a * P8 (log a)⁻¹ = a / log a + a / log a ^ 2 + 2 * (a / log a ^ 3)
      + 6 * (a / log a ^ 4) + 24 * (a / log a ^ 5) + 120 * (a / log a ^ 6)
      + 720 * (a / log a ^ 7) + 5040 * (a / log a ^ 8) := by
    unfold P8
    ring
  rw [ex, ea]
  linarith

/-- The remainder integral beyond `e¹⁶`. -/
theorem integral_inv_log_pow_nine_bounds {a x : ℝ} (ha : exp 16 ≤ a) (hax : a ≤ x) :
    0 ≤ ∫ t in a..x, 1 / log t ^ 9 ∧
      ∫ t in a..x, 1 / log t ^ 9 ≤ 16 / 7 * (x / log x ^ 9) := by
  have he : (17 : ℝ) ≤ exp 16 := by
    have := Real.add_one_le_exp (16 : ℝ)
    linarith
  have ha1 : 1 < a := by linarith
  have hx1 : 1 < x := by linarith
  have hlog : ∀ t : ℝ, a ≤ t → 16 ≤ log t := fun t ht => by
    have := log_le_log (exp_pos 16) (ha.trans ht)
    rwa [log_exp] at this
  constructor
  · exact intervalIntegral.integral_nonneg hax fun t ht =>
      div_nonneg zero_le_one (pow_nonneg (by linarith [hlog t ht.1]) _)
  · have hint1 := intervalIntegrable_inv_log_pow 9 ha1 hx1
    have hint2 : IntervalIntegrable (fun t : ℝ => ((9 : ℕ) : ℝ) / log t ^ (9 + 1)) volume a x := by
      simpa [div_eq_mul_inv] using
        (intervalIntegrable_inv_log_pow (9 + 1) ha1 hx1).const_mul ((9 : ℕ) : ℝ)
    have hmono : ∫ t in a..x, 1 / log t ^ 9
        ≤ ∫ t in a..x, 16 / 7 * (1 / log t ^ 9 - ((9 : ℕ) : ℝ) / log t ^ (9 + 1)) := by
      refine intervalIntegral.integral_mono_on hax hint1 ((hint1.sub hint2).const_mul (16 / 7)) ?_
      intro t ht
      have h16 := hlog t ht.1
      have hl : 0 < log t := by linarith
      have hk : ((9 : ℕ) : ℝ) / log t ^ (9 + 1) ≤ 9 / 16 * (1 / log t ^ 9) := by
        rw [pow_succ, div_le_iff₀ (by positivity)]
        have : 9 / 16 * (1 / log t ^ 9) * (log t ^ 9 * log t) = 9 / 16 * log t := by
          field_simp
        rw [this]
        push_cast
        linarith
      linarith
    rw [intervalIntegral.integral_const_mul, integral_inv_log_pow_sub 9 ha1 hx1] at hmono
    have : 0 ≤ a / log a ^ 9 := div_nonneg (by linarith) (pow_nonneg (by linarith [hlog a le_rfl]) _)
    linarith

/-- Explicit two-sided bounds for `∫_2^x dt/log t`, `x ≥ e¹⁶`. -/
theorem Li_bounds {x : ℝ} (hx : exp 16 ≤ x) :
    x * P8 (log x)⁻¹ - 12999998 ≤ ∫ t in (2 : ℝ)..x, 1 / log t ∧
      ∫ t in (2 : ℝ)..x, 1 / log t ≤ x * S8 (log x)⁻¹ + 12999998 := by
  have he : (17 : ℝ) ≤ exp 16 := by
    have := Real.add_one_le_exp (16 : ℝ)
    linarith
  have hub := exp_sixteen_le
  have ha1 : (1 : ℝ) < exp 16 := by linarith
  have hx1 : 1 < x := by linarith
  have hsplit := intervalIntegral.integral_add_adjacent_intervals
    (intervalIntegrable_inv_log (a := 2) (b := exp 16) (by norm_num) ha1)
    (intervalIntegrable_inv_log (a := exp 16) (b := x) ha1 hx1)
  -- the head `∫_2^{e¹⁶}`
  have hl2 : 0 < log (2 : ℝ) := log_pos (by norm_num)
  have hhead0 : 0 ≤ ∫ t in (2 : ℝ)..exp 16, 1 / log t :=
    intervalIntegral.integral_nonneg (by linarith) fun t ht =>
      div_nonneg zero_le_one (log_pos (by linarith [ht.1])).le
  have hhead : ∫ t in (2 : ℝ)..exp 16, 1 / log t ≤ 1.45 * exp 16 := by
    have hmono : ∫ t in (2 : ℝ)..exp 16, 1 / log t ≤ ∫ _t in (2 : ℝ)..exp 16, 1 / log 2 := by
      refine intervalIntegral.integral_mono_on (by linarith)
        (intervalIntegrable_inv_log (by norm_num) ha1) intervalIntegrable_const ?_
      intro t ht
      exact one_div_le_one_div_of_le hl2 (log_le_log (by norm_num) ht.1)
    rw [intervalIntegral.integral_const, smul_eq_mul] at hmono
    have hinv : 1 / log 2 ≤ 1.45 := by
      rw [div_le_iff₀ hl2]
      have := Real.log_two_gt_d9
      linarith
    have : (exp 16 - 2) * (1 / log 2) ≤ (exp 16 - 2) * 1.45 :=
      mul_le_mul_of_nonneg_left hinv (by linarith)
    linarith
  -- the tail `∫_{e¹⁶}^x`
  have htail := integral_inv_log_eq_P8 (a := exp 16) (x := x) ha1 hx1
  rw [log_exp] at htail
  obtain ⟨h90, h9⟩ := integral_inv_log_pow_nine_bounds (a := exp 16) (x := x) le_rfl hx
  have hP16 : 0 ≤ P8 (16 : ℝ)⁻¹ ∧ P8 (16 : ℝ)⁻¹ ≤ 0.07 := by
    unfold P8
    constructor <;> norm_num
  have hS : x * S8 (log x)⁻¹ = x * P8 (log x)⁻¹ + 92160 * (x / log x ^ 9) := by
    unfold S8
    rw [inv_pow]
    ring
  have haP : 0 ≤ exp 16 * P8 (16 : ℝ)⁻¹ := mul_nonneg (by linarith) hP16.1
  have haP' : exp 16 * P8 (16 : ℝ)⁻¹ ≤ exp 16 * 0.07 :=
    mul_le_mul_of_nonneg_left hP16.2 (by linarith)
  constructor
  · rw [← hsplit, htail]
    linarith
  · rw [← hsplit, htail, hS]
    linarith

/-! ### The explicit inequality -/

/-- The real-variable core: with `π(x) ≤ x S + K + A r L` and
`π(x/e) ≥ (x/E) Pv − K − A q (L − 1)`, the square of the first bound is below
`E x u` times the second. -/
theorem explicit_core (x L u r q E A K S Pv : ℝ)
    (hx : x = r * r) (hr : 485000000 ≤ r) (hL : 40 ≤ L) (hu : u = L⁻¹)
    (hE0 : 0 < E) (hE1 : E ≤ 2.72) (hq0 : 0 ≤ q) (hEq : E * q ≤ 1.649 * r)
    (hA0 : 0 ≤ A) (hA : A ≤ 1 / 25) (hK0 : 0 ≤ K) (hK : K ≤ 13000000)
    (hS0 : 0 ≤ S) (hS : S ≤ 1.03 * u)
    (hgap : u ^ 6 * (1 + 14 * u) ≤ u * Pv - S ^ 2)
    (hgrow : 0.159 * L ^ 7 ≤ r * (L + 14)) :
    (x * S + K + A * (r * L)) ^ 2 < E * x * u * (x / E * Pv - K - A * (q * (L - 1))) := by
  have hL0 : 0 < L := by linarith
  have hr0 : 0 < r := by linarith
  have hx0 : 0 < x := by rw [hx]; positivity
  have hu0 : 0 < u := by rw [hu]; positivity
  have huL : u * L = 1 := by rw [hu]; field_simp
  have hu40 : u ≤ 1 / 40 := by
    rw [hu, inv_eq_one_div]
    exact one_div_le_one_div_of_le (by norm_num) hL
  have hxr0 : 0 < x * r := by positivity
  -- `L² ≤ r / 100000`
  have hL2 : L ^ 2 ≤ r / 100000 := by
    have h4 : (2560000 : ℝ) ≤ L ^ 4 := by
      have := pow_le_pow_left₀ (by norm_num : (0 : ℝ) ≤ 40) hL 4
      norm_num at this
      exact this
    have h5 : r * (L + 14) ≤ r * (1.35 * L) := mul_le_mul_of_nonneg_left (by linarith) hr0.le
    have h6 : 0.159 * 2560000 * (L ^ 2 * L) ≤ 0.159 * L ^ 7 := by
      have e : 0.159 * L ^ 7 = 0.159 * L ^ 4 * (L ^ 2 * L) := by ring
      rw [e]
      have : (0 : ℝ) ≤ L ^ 2 * L := by positivity
      nlinarith
    have h7 : (407040 * L ^ 2) * L ≤ (1.35 * r) * L := by nlinarith
    have h8 := le_of_mul_le_mul_right h7 hL0
    linarith
  -- the gap
  have g1 : 0.159 * (x * r) ≤ x ^ 2 * (u ^ 6 * (1 + 14 * u)) := by
    have e : u ^ 6 * (1 + 14 * u) = (L + 14) / L ^ 7 := by
      rw [hu]
      field_simp
    rw [e, ← mul_div_assoc, le_div_iff₀ (by positivity)]
    have h := mul_le_mul_of_nonneg_left hgrow hxr0.le
    have e2 : x * r * (r * (L + 14)) = x ^ 2 * (L + 14) := by
      rw [hx]
      ring
    linarith
  have g2 : 0.159 * (x * r) ≤ x ^ 2 * (u * Pv) - x ^ 2 * S ^ 2 := by
    have h := mul_le_mul_of_nonneg_left hgap (sq_nonneg x)
    have e : x ^ 2 * (u * Pv - S ^ 2) = x ^ 2 * (u * Pv) - x ^ 2 * S ^ 2 := by ring
    linarith
  -- the five error terms
  have hS' : S ≤ 0.02575 := by linarith
  have b1 : 2 * x * S * K ≤ 0.0015 * (x * r) := by
    have hSK : S * K ≤ 0.02575 * 13000000 := mul_le_mul hS' hK hK0 (by norm_num)
    have h : 2 * (S * K) ≤ 0.0015 * r := by linarith
    calc 2 * x * S * K = x * (2 * (S * K)) := by ring
      _ ≤ x * (0.0015 * r) := mul_le_mul_of_nonneg_left h hx0.le
      _ = 0.0015 * (x * r) := by ring
  have b2 : 2 * A * x * r * (S * L) ≤ 0.0824 * (x * r) := by
    have hSL : S * L ≤ 1.03 := by
      calc S * L ≤ 1.03 * u * L := mul_le_mul_of_nonneg_right hS hL0.le
        _ = 1.03 := by rw [mul_assoc, huL, mul_one]
    have hSL0 : 0 ≤ S * L := mul_nonneg hS0 hL0.le
    have h : A * (S * L) ≤ 1 / 25 * 1.03 := mul_le_mul hA hSL hSL0 (by norm_num)
    calc 2 * A * x * r * (S * L) = (x * r) * (2 * (A * (S * L))) := by ring
      _ ≤ (x * r) * (2 * (1 / 25 * 1.03)) := mul_le_mul_of_nonneg_left (by linarith) hxr0.le
      _ = 0.0824 * (x * r) := by ring
  have b3 : (K + A * (r * L)) ^ 2 ≤ 0.0001 * (x * r) := by
    have h1 : (K + A * (r * L)) ^ 2 ≤ 2 * K ^ 2 + 2 * (A * (r * L)) ^ 2 := by
      nlinarith [sq_nonneg (K - A * (r * L))]
    have h2 : 2 * K ^ 2 ≤ 0.00005 * (x * r) := by
      have hK2 : K ^ 2 ≤ 13000000 ^ 2 := pow_le_pow_left₀ hK0 hK 2
      have hrr : (485000000 : ℝ) * 485000000 ≤ r * r := mul_le_mul hr hr (by norm_num) hr0.le
      have hxr : (485000000 : ℝ) * 485000000 * 485000000 ≤ x * r := by
        rw [hx]
        exact mul_le_mul hrr hr (by norm_num) (by positivity)
      norm_num at hK2 hxr
      linarith
    have h3 : 2 * (A * (r * L)) ^ 2 ≤ 0.00005 * (x * r) := by
      have hA2 : A ^ 2 ≤ (1 / 25) ^ 2 := pow_le_pow_left₀ hA0 hA 2
      have e : (A * (r * L)) ^ 2 = A ^ 2 * (x * L ^ 2) := by
        rw [hx]
        ring
      rw [e]
      have hxL : x * L ^ 2 ≤ x * (r / 100000) := mul_le_mul_of_nonneg_left hL2 hx0.le
      have h : A ^ 2 * (x * L ^ 2) ≤ (1 / 25) ^ 2 * (x * (r / 100000)) :=
        mul_le_mul hA2 hxL (by positivity) (by positivity)
      have e2 : ((1 : ℝ) / 25) ^ 2 * (x * (r / 100000)) = (1 / 62500000) * (x * r) := by ring
      rw [e2] at h
      linarith
    linarith
  have b4 : E * u * x * K ≤ 0.002 * (x * r) := by
    have hEu : E * u ≤ 2.72 * (1 / 40) := mul_le_mul hE1 hu40 hu0.le (by norm_num)
    have hEuK : E * u * K ≤ 2.72 * (1 / 40) * 13000000 := mul_le_mul hEu hK hK0 (by norm_num)
    have h : E * u * K ≤ 0.002 * r := by linarith
    calc E * u * x * K = x * (E * u * K) := by ring
      _ ≤ x * (0.002 * r) := mul_le_mul_of_nonneg_left h hx0.le
      _ = 0.002 * (x * r) := by ring
  have b5 : A * x * (E * q) * (u * (L - 1)) ≤ 0.066 * (x * r) := by
    have h1 : u * (L - 1) ≤ 1 := by
      have : u * (L - 1) = 1 - u := by rw [mul_sub, huL, mul_one]
      linarith
    have h1' : 0 ≤ u * (L - 1) := mul_nonneg hu0.le (by linarith)
    have hEq0 : 0 ≤ E * q := mul_nonneg hE0.le hq0
    have h2 : A * (E * q) ≤ 1 / 25 * (1.649 * r) := mul_le_mul hA hEq hEq0 (by norm_num)
    have h3 : A * (E * q) * (u * (L - 1)) ≤ 1 / 25 * (1.649 * r) * 1 :=
      mul_le_mul h2 h1 h1' (by positivity)
    calc A * x * (E * q) * (u * (L - 1)) = x * (A * (E * q) * (u * (L - 1))) := by ring
      _ ≤ x * (1 / 25 * (1.649 * r) * 1) := mul_le_mul_of_nonneg_left h3 hx0.le
      _ = 0.06596 * (x * r) := by ring
      _ ≤ 0.066 * (x * r) := mul_le_mul_of_nonneg_right (by norm_num) hxr0.le
  have eU : (x * S + K + A * (r * L)) ^ 2
      = x ^ 2 * S ^ 2 + 2 * x * S * K + 2 * A * x * r * (S * L) + (K + A * (r * L)) ^ 2 := by
    ring
  have eW : E * x * u * (x / E * Pv - K - A * (q * (L - 1)))
      = x ^ 2 * (u * Pv) - E * u * x * K - A * x * (E * q) * (u * (L - 1)) := by
    field_simp
  rw [eU, eW]
  linarith

/-- **Explicit criterion.**  If `|π(y) − (c + ∫_2^y dt/log t)| ≤ A √y log y` for all `y ≥ e³⁹`,
with `0 ≤ A ≤ 1/25` and `|c| ≤ 2`, then Ramanujan's inequality holds for every real
`x ≥ e⁴⁰`. -/
theorem ramanujan_of_sqrt_log_bound (A c : ℝ) (hA0 : 0 ≤ A) (hA : A ≤ 1 / 25) (hc : |c| ≤ 2)
    (h : ∀ y : ℝ, exp 39 ≤ y →
      |(Nat.primeCounting ⌊y⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..y, 1 / log t)| ≤ A * (√y * log y))
    (x : ℝ) (hx : exp 40 ≤ x) :
    (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
  have hE0 : 0 < exp 1 := exp_pos 1
  have hE1 : exp 1 ≤ 2.72 := by
    have := Real.exp_one_lt_d9
    linarith
  have hx0 : 0 < x := lt_of_lt_of_le (exp_pos 40) hx
  have hL : 40 ≤ log x := by rwa [le_log_iff_exp_le hx0]
  have hy39 : exp 39 ≤ x / exp 1 := by
    rw [le_div_iff₀ hE0, ← exp_add]
    norm_num
    exact hx
  have hlogy : log (x / exp 1) = log x - 1 := by rw [log_div hx0.ne' hE0.ne', log_exp]
  have h39x : exp 39 ≤ x := le_trans (exp_le_exp.mpr (by norm_num)) hx
  have h16x : exp 16 ≤ x := le_trans (exp_le_exp.mpr (by norm_num)) hx
  have h16y : exp 16 ≤ x / exp 1 := le_trans (exp_le_exp.mpr (by norm_num)) hy39
  have hLiU := (Li_bounds h16x).2
  have hLiL := (Li_bounds h16y).1
  have hπx := (abs_le.mp (h x h39x)).2
  have hπy := (abs_le.mp (h (x / exp 1) hy39)).1
  rw [hlogy] at hLiL hπy
  have hc' := abs_le.mp hc
  have hsx : √x = exp (log x / 2) := by
    conv_lhs => rw [← exp_log hx0]
    exact sqrt_exp_eq _
  have hr : 485000000 ≤ √x := by
    rw [hsx]
    exact exp_twenty_ge.trans (exp_le_exp.mpr (by linarith))
  have hgrow : 0.159 * log x ^ 7 ≤ √x * (log x + 14) := by
    rw [hsx]
    exact growth hL
  have hEq : exp 1 * √(x / exp 1) ≤ 1.649 * √x := by
    have h1 : (exp 1 * √(x / exp 1)) ^ 2 = exp 1 * x := by
      rw [mul_pow, sq_sqrt (by positivity)]
      field_simp
    have h2 : (1.649 * √x) ^ 2 = 2.719201 * x := by
      rw [mul_pow, sq_sqrt hx0.le]
      norm_num
    have h3 : (exp 1 * √(x / exp 1)) ^ 2 ≤ (1.649 * √x) ^ 2 := by
      rw [h1, h2]
      exact mul_le_mul_of_nonneg_right (by linarith [Real.exp_one_lt_d9]) hx0.le
    exact (pow_le_pow_iff_left₀ (by positivity) (by positivity) (by norm_num)).mp h3
  have hL0 : 0 < log x := by linarith
  have hL1 : 0 < log x - 1 := by linarith
  have hu0 : 0 < (log x)⁻¹ := inv_pos.mpr hL0
  have hu40 : (log x)⁻¹ ≤ 1 / 40 := by
    rw [inv_eq_one_div]
    exact one_div_le_one_div_of_le (by norm_num) hL
  have hv : (log x - 1)⁻¹ = (log x)⁻¹ / (1 - (log x)⁻¹) := by
    field_simp
  rw [hv] at hLiL
  obtain ⟨hS0, hS1⟩ := S8_bounds hu0.le hu40
  have hU : (Nat.primeCounting ⌊x⌋₊ : ℝ)
      ≤ x * S8 (log x)⁻¹ + 13000000 + A * (√x * log x) := by linarith [hc'.2]
  have hW : x / exp 1 * P8 ((log x)⁻¹ / (1 - (log x)⁻¹)) - 13000000
      - A * (√(x / exp 1) * (log x - 1)) ≤ (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
    linarith [hc'.1]
  have hcore := explicit_core x (log x) (log x)⁻¹ (√x) (√(x / exp 1)) (exp 1) A 13000000
    (S8 (log x)⁻¹) (P8 ((log x)⁻¹ / (1 - (log x)⁻¹))) (mul_self_sqrt hx0.le).symm hr hL rfl hE0 hE1
    (sqrt_nonneg _) hEq hA0 hA (by norm_num) le_rfl hS0 hS1 (gap_lower8 hu0 hu40) hgrow
  have hπ0 : (0 : ℝ) ≤ (Nat.primeCounting ⌊x⌋₊ : ℝ) := Nat.cast_nonneg _
  calc (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      ≤ (x * S8 (log x)⁻¹ + 13000000 + A * (√x * log x)) ^ 2 := pow_le_pow_left₀ hπ0 hU 2
    _ < exp 1 * x * (log x)⁻¹ * (x / exp 1 * P8 ((log x)⁻¹ / (1 - (log x)⁻¹)) - 13000000
          - A * (√(x / exp 1) * (log x - 1))) := hcore
    _ ≤ exp 1 * x * (log x)⁻¹ * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) :=
        mul_le_mul_of_nonneg_left hW (by positivity)
    _ = exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by ring

end RamanujanLib
