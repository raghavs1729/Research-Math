/-
Copyright (c) 2026. Released under the repository's license.

# The final elementary inequality

For `t ≥ 10^7` all the lower-order terms produced by the elementary estimates are
dominated by the main term `(2/3) log 2 · t^{3/2}` (`Allx.final_ineq`).

Method: every term is bounded by an explicit constant times `X = t^{3/2}`, using
`t^{b-a} ≤ T₀^{-a} t^b` and `log t ≤ log T₀ · T₀^{-a} · t^a` for `t ≥ T₀ = 10^7`
(the latter needs `a log T₀ ≥ 1`), together with rational bounds for `T₀^{-a}`,
`log T₀` and `log 2`.
-/
import Mathlib

open Real

noncomputable section

namespace Allx

/-- `t^{b-a} ≤ c · t^b` for `t ≥ T > 0`, `a ≥ 0`, if `T^{-a} ≤ c`. -/
lemma rpow_le_mul_rpow {T t a b e c : ℝ} (hT : 0 < T) (hTt : T ≤ t) (ha : 0 ≤ a)
    (he : e = b - a) (hc : T ^ (-a) ≤ c) : t ^ e ≤ c * t ^ b := by
  have ht : 0 < t := lt_of_lt_of_le hT hTt
  have h1 : t ^ e = t ^ b * t ^ (-a) := by
    rw [he, sub_eq_add_neg, Real.rpow_add ht]
  have h2 : t ^ (-a) ≤ T ^ (-a) := Real.rpow_le_rpow_of_nonpos hT hTt (by linarith)
  have h3 : 0 ≤ t ^ b := Real.rpow_nonneg ht.le b
  rw [h1, mul_comm c]
  exact mul_le_mul_of_nonneg_left (h2.trans hc) h3

/-- `log t ≤ log T · T^{-a} · t^a` for `t ≥ T > 0` provided `a log T ≥ 1`. -/
lemma log_le_mul_rpow {T t a : ℝ} (hT : 0 < T) (hTt : T ≤ t) (ha : 0 < a)
    (h1 : 1 ≤ a * log T) : log t ≤ log T * T ^ (-a) * t ^ a := by
  have ht : 0 < t := lt_of_lt_of_le hT hTt
  have hy : (t / T) ^ a = T ^ (-a) * t ^ a := by
    rw [Real.div_rpow ht.le hT.le, Real.rpow_neg hT.le, div_eq_mul_inv, mul_comm]
  have hq : 1 ≤ t / T := by rw [le_div_iff₀ hT]; linarith
  have hy1 : 1 ≤ (t / T) ^ a := Real.one_le_rpow hq ha.le
  have hlog : log t = log T + log (t / T) := by
    rw [Real.log_div ht.ne' hT.ne']; ring
  have h2 : a * log (t / T) ≤ (t / T) ^ a - 1 := by
    rw [← Real.log_rpow (div_pos ht hT)]
    exact Real.log_le_sub_one_of_pos (Real.rpow_pos_of_pos (div_pos ht hT) a)
  rw [mul_assoc, ← hy, hlog]
  set y := (t / T) ^ a
  set Lr := log (t / T)
  have h3 : (y - 1) * 1 ≤ (y - 1) * (a * log T) :=
    mul_le_mul_of_nonneg_left h1 (by linarith)
  have h4 : a * Lr ≤ a * ((y - 1) * log T) := by nlinarith
  have h5 : Lr ≤ (y - 1) * log T := le_of_mul_le_mul_left h4 ha
  linarith

/-- `t^{b-a} log t ≤ L · c · t^b` for `t ≥ T`, if `log T ≤ L`, `T^{-a} ≤ c`, `a log T ≥ 1`. -/
lemma rpow_mul_log_le {T t a b e c L : ℝ} (hT : 0 < T) (hTt : T ≤ t) (ha : 0 < a)
    (he : e = b - a) (h1 : 1 ≤ a * log T) (hL : log T ≤ L) (hc : T ^ (-a) ≤ c) :
    t ^ e * log t ≤ L * c * t ^ b := by
  have ht : 0 < t := lt_of_lt_of_le hT hTt
  have hlog := log_le_mul_rpow hT hTt ha h1
  have hlT : 0 ≤ log T := by nlinarith
  have hTa : 0 ≤ T ^ (-a) := Real.rpow_nonneg hT.le _
  have hb : t ^ e * t ^ a = t ^ b := by
    rw [← Real.rpow_add ht, he]; congr 1; ring
  have hte : 0 ≤ t ^ e := Real.rpow_nonneg ht.le e
  have htb : 0 ≤ t ^ b := Real.rpow_nonneg ht.le b
  calc t ^ e * log t ≤ t ^ e * (log T * T ^ (-a) * t ^ a) :=
        mul_le_mul_of_nonneg_left hlog hte
    _ = log T * T ^ (-a) * t ^ b := by rw [← hb]; ring
    _ ≤ L * c * t ^ b := by
        apply mul_le_mul_of_nonneg_right _ htb
        exact mul_le_mul hL hc hTa (hlT.trans hL)

/-- `T^{-a} ≤ c` follows from `1 ≤ c^n T^k` when `a n = k`. -/
lemma rpow_neg_le_of_pow {T c a : ℝ} (n k : ℕ) (hn : n ≠ 0) (hT : 0 < T) (hc : 0 ≤ c)
    (hak : a * n = k) (h : 1 ≤ c ^ n * T ^ k) : T ^ (-a) ≤ c := by
  have h0 : 0 ≤ T ^ (-a) := Real.rpow_nonneg hT.le _
  rw [← pow_le_pow_iff_left₀ h0 hc hn]
  have h1 : (T ^ (-a)) ^ n = (T ^ k)⁻¹ := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hT.le, neg_mul, hak, Real.rpow_neg hT.le,
      Real.rpow_natCast]
  rw [h1]
  have hTk : 0 < T ^ k := pow_pos hT k
  rw [inv_le_iff_one_le_mul₀ hTk]
  exact h

/-- **The final elementary inequality.** -/
theorem final_ineq {t : ℝ} (ht : 10000000 ≤ t) :
    27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 4 * t + log t + 4 + 1 / 3 * log (t + 1)
      + t * (log 2 + log (√t + 2)) + t / 2 * log t
      + 3 / 4 * log 4 * t ^ (4 / 3 : ℝ) + 4 / 7 * t ^ (7 / 6 : ℝ) * log t
      + 5 / 6 * log 4 * t ^ (6 / 5 : ℝ) + 4 / 11 * t ^ (11 / 10 : ℝ) * log t
    < 2 / 3 * log 2 * t ^ (3 / 2 : ℝ) := by
  have hT : (0 : ℝ) < 10000000 := by norm_num
  have ht0 : 0 < t := by linarith
  set X := t ^ (3 / 2 : ℝ) with hX
  have hXpos : 0 < X := Real.rpow_pos_of_pos ht0 _
  -- logarithms of constants
  have hl2lo : (0.6931 : ℝ) < log 2 := lt_trans (by norm_num) Real.log_two_gt_d9
  have hl2hi : log 2 < (0.6932 : ℝ) := lt_trans Real.log_two_lt_d9 (by norm_num)
  have hl4 : log 4 = 2 * log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  have hl10 : 3 * log 10 ≤ 10 * log 2 := by
    have h : log ((10 : ℝ) ^ 3) ≤ log ((2 : ℝ) ^ 10) :=
      Real.log_le_log (by norm_num) (by norm_num)
    rw [Real.log_pow, Real.log_pow] at h
    norm_num at h
    linarith
  have hl10' : log 2 ≤ log 10 := Real.log_le_log (by norm_num) (by norm_num)
  have hL0 : log 10000000 = 7 * log 10 := by
    rw [show (10000000 : ℝ) = 10 ^ 7 by norm_num, Real.log_pow]; norm_num
  have hL0hi : log 10000000 ≤ (16.18 : ℝ) := by rw [hL0]; linarith
  have hL0lo : (3 : ℝ) ≤ log 10000000 := by rw [hL0]; linarith
  -- powers of `T₀`
  have c12 : (10000000 : ℝ) ^ (-(1 / 2 : ℝ)) ≤ 0.00032 :=
    rpow_neg_le_of_pow 2 1 (by norm_num) hT (by norm_num) (by norm_num) (by norm_num)
  have c16 : (10000000 : ℝ) ^ (-(1 / 6 : ℝ)) ≤ 0.0682 :=
    rpow_neg_le_of_pow 6 1 (by norm_num) hT (by norm_num) (by norm_num) (by norm_num)
  have c13 : (10000000 : ℝ) ^ (-(1 / 3 : ℝ)) ≤ 0.00466 :=
    rpow_neg_le_of_pow 3 1 (by norm_num) hT (by norm_num) (by norm_num) (by norm_num)
  have c310 : (10000000 : ℝ) ^ (-(3 / 10 : ℝ)) ≤ 0.008 :=
    rpow_neg_le_of_pow 10 3 (by norm_num) hT (by norm_num) (by norm_num) (by norm_num)
  have c25 : (10000000 : ℝ) ^ (-(2 / 5 : ℝ)) ≤ 0.00159 :=
    rpow_neg_le_of_pow 5 2 (by norm_num) hT (by norm_num) (by norm_num) (by norm_num)
  -- pure powers
  have p12 : t ≤ 0.00032 * X := by
    have h := rpow_le_mul_rpow (e := 1) (b := 3 / 2) hT ht (by norm_num) (by norm_num) c12
    rwa [Real.rpow_one] at h
  have p16 : t ^ (4 / 3 : ℝ) ≤ 0.0682 * X :=
    rpow_le_mul_rpow (b := 3 / 2) hT ht (by norm_num) (by norm_num) c16
  have p310 : t ^ (6 / 5 : ℝ) ≤ 0.008 * X :=
    rpow_le_mul_rpow (b := 3 / 2) hT ht (by norm_num) (by norm_num) c310
  -- powers times `log t`
  have q12 : t * log t ≤ 16.18 * 0.00032 * X := by
    have h := rpow_mul_log_le (e := 1) (b := 3 / 2) hT ht (by norm_num) (by norm_num)
      (by linarith) hL0hi c12
    rwa [Real.rpow_one] at h
  have q13 : t ^ (7 / 6 : ℝ) * log t ≤ 16.18 * 0.00466 * X :=
    rpow_mul_log_le (b := 3 / 2) hT ht (by norm_num) (by norm_num) (by linarith) hL0hi c13
  have q25 : t ^ (11 / 10 : ℝ) * log t ≤ 16.18 * 0.00159 * X :=
    rpow_mul_log_le (b := 3 / 2) hT ht (by norm_num) (by norm_num) (by linarith) hL0hi c25
  -- the term `(t + 1)^{3/2}`
  have hA : (t + 1) ^ (3 / 2 : ℝ) ≤ 1.0001 * X := by
    have h1 : t + 1 ≤ t * (1 + 1 / 10000000) := by nlinarith
    have h2 : (t + 1) ^ (3 / 2 : ℝ) ≤ (t * (1 + 1 / 10000000)) ^ (3 / 2 : ℝ) :=
      Real.rpow_le_rpow (by linarith) h1 (by norm_num)
    rw [Real.mul_rpow ht0.le (by norm_num)] at h2
    have h3 : ((1 : ℝ) + 1 / 10000000) ^ (3 / 2 : ℝ) ≤ ((1 : ℝ) + 1 / 10000000) ^ (2 : ℝ) :=
      Real.rpow_le_rpow_of_exponent_le (by norm_num) (by norm_num)
    rw [Real.rpow_two] at h3
    have h4 : ((1 : ℝ) + 1 / 10000000) ^ (3 / 2 : ℝ) ≤ 1.0001 := h3.trans (by norm_num)
    calc (t + 1) ^ (3 / 2 : ℝ) ≤ X * (1 + 1 / 10000000) ^ (3 / 2 : ℝ) := h2
      _ ≤ X * 1.0001 := mul_le_mul_of_nonneg_left h4 hXpos.le
      _ = 1.0001 * X := mul_comm _ _
  -- small logarithms
  have hlogt : log t ≤ t := by linarith [Real.log_le_sub_one_of_pos ht0]
  have hlogt1 : log (t + 1) ≤ t := by
    linarith [Real.log_le_sub_one_of_pos (show (0 : ℝ) < t + 1 by linarith)]
  have hlogt0 : 0 ≤ log t := Real.log_nonneg (by linarith)
  -- `log (√t + 2) ≤ log 2 + log t / 2`
  have hsqrt2 : (2 : ℝ) ≤ √t := by
    rw [Real.le_sqrt' (by norm_num)]; linarith
  have hlogs : log (√t + 2) ≤ log 2 + log t / 2 := by
    have h : log (√t + 2) ≤ log (2 * √t) :=
      Real.log_le_log (by positivity) (by linarith)
    rw [Real.log_mul (by norm_num) (by positivity), Real.log_sqrt ht0.le] at h
    exact h
  have hD : t * (log 2 + log (√t + 2)) ≤ 2 * 0.6932 * t + 1 / 2 * (t * log t) := by
    have h1 : t * (log 2 + log (√t + 2)) ≤ t * (2 * 0.6932 + log t / 2) :=
      mul_le_mul_of_nonneg_left (by linarith) ht0.le
    linarith
  -- products with `log 4`, `log 2`
  have hl4hi : log 4 ≤ (1.3864 : ℝ) := by rw [hl4]; linarith
  have hl4nn : 0 ≤ log (4 : ℝ) := by rw [hl4]; linarith
  have hF : log 4 * t ^ (4 / 3 : ℝ) ≤ 1.3864 * (0.0682 * X) :=
    mul_le_mul hl4hi p16 (Real.rpow_nonneg ht0.le _) (by norm_num)
  have hH : log 4 * t ^ (6 / 5 : ℝ) ≤ 1.3864 * (0.008 * X) :=
    mul_le_mul hl4hi p310 (Real.rpow_nonneg ht0.le _) (by norm_num)
  have hR : 0.6931 * X ≤ log 2 * X := mul_le_mul_of_nonneg_right hl2lo.le hXpos.le
  linarith [hA, p12, q12, q13, q25, hlogt, hlogt1, hD, hF, hH, hR, hXpos]

end Allx
