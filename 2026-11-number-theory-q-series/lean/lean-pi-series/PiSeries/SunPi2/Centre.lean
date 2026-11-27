/-
# WP-C: the centre term.  `cf N ⌊x* N⌋₊ · √N → √3 / (√2 π^{5/2})`, by Stirling's formula for the
# eleven factorials of `cf` (N2-PLAN.md, sketch (C); numerics in `derivation/wpC_check.py`).

With `k0 = ⌊x* N⌋₊ = x* N − θ_N`, `θ_N ∈ [0,1)`, every factorial argument is `m = c N + (r − q θ_N)`.
The Stirling main terms `St c u N` add up (with the signs of `cf`) to
`22 N log 2 − 12 log 2 − (5/2) log N − (5/2) log (2π)` identically (`main_term`), so
`cf N k0 · √N = Rf N · (4N−1)(2k0+1)/N²` with `log (Rf N) → −(5/2) log (2π)`.
-/
import PiSeries.SunPi2.Defs
import PiSeries.SunPi2.Stirling

open Filter Topology

namespace PiSeries.SunPi2

/-! ### `√3`, `x*`, the fractional part -/

lemma sqrt3_sq : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)

lemma sqrt3_gt : (17 / 10 : ℝ) < Real.sqrt 3 := by
  nlinarith [sqrt3_sq, Real.sqrt_nonneg 3]

lemma sqrt3_lt : Real.sqrt 3 < (18 / 10 : ℝ) := by
  nlinarith [sqrt3_sq, Real.sqrt_nonneg 3]

lemma xstar_gt_C : (17 / 20 : ℝ) < xstar := by unfold xstar; linarith [sqrt3_gt]

lemma xstar_lt_C : xstar < (9 / 10 : ℝ) := by unfold xstar; linarith [sqrt3_lt]

/-- `k0 N = ⌊x* N⌋₊` -/
noncomputable def k0 (N : ℕ) : ℕ := ⌊xstar * (N : ℝ)⌋₊

/-- `θ_N = x* N − ⌊x* N⌋₊ ∈ [0,1)` -/
noncomputable def θc (N : ℕ) : ℝ := xstar * N - k0 N

lemma θc_nonneg (N : ℕ) : 0 ≤ θc N := by
  have : (0 : ℝ) ≤ xstar * N := mul_nonneg (by linarith [xstar_gt_C]) (Nat.cast_nonneg N)
  exact sub_nonneg.2 (Nat.floor_le this)

lemma θc_lt_one (N : ℕ) : θc N < 1 := by
  have := Nat.lt_floor_add_one (xstar * (N : ℝ))
  unfold θc k0; linarith

lemma k0_cast (N : ℕ) : (k0 N : ℝ) = xstar * N - θc N := by unfold θc; ring

lemma k0_lt (N : ℕ) (hN : 1 ≤ N) : k0 N + 1 ≤ N := by
  have hN' : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have h1 : (k0 N : ℝ) < N := by
    have := θc_nonneg N
    rw [k0_cast]; nlinarith [xstar_lt_C]
  have : k0 N < N := by exact_mod_cast h1
  omega

lemma two_N_le (N : ℕ) : 2 * N ≤ 4 * k0 N + 4 := by
  have h1 : (2 * N : ℝ) < 4 * k0 N + 4 := by
    have := θc_lt_one N
    have hN' : (0 : ℝ) ≤ N := Nat.cast_nonneg N
    rw [k0_cast]; nlinarith [xstar_gt_C]
  have : 2 * N < 4 * k0 N + 4 := by exact_mod_cast h1
  omega

/-! ### The logarithms of the constants `c`, in terms of `log 2` and `log (√3 − 1)` -/

lemma sqrt3_sub_one_pos : 0 < Real.sqrt 3 - 1 := by linarith [sqrt3_gt]

lemma log_c2 : Real.log (8 * xstar + 4) = 3 * Real.log 2 - Real.log (Real.sqrt 3 - 1) := by
  have hs := sqrt3_sq
  have hpos := sqrt3_sub_one_pos
  have : 8 * xstar + 4 = 2 ^ 3 / (Real.sqrt 3 - 1) := by
    unfold xstar; rw [eq_div_iff hpos.ne']; linear_combination 4 * hs
  rw [this, Real.log_div (by positivity) hpos.ne', Real.log_pow]; push_cast; ring

lemma log_c3 : Real.log (4 * xstar - 2) = Real.log 2 + Real.log (Real.sqrt 3 - 1) := by
  have hpos := sqrt3_sub_one_pos
  have : 4 * xstar - 2 = 2 * (Real.sqrt 3 - 1) := by unfold xstar; ring
  rw [this, Real.log_mul (by norm_num) hpos.ne']

lemma log_c6 : Real.log (4 * xstar + 2) = 2 * Real.log 2 - Real.log (Real.sqrt 3 - 1) := by
  have hs := sqrt3_sq
  have hpos := sqrt3_sub_one_pos
  have : 4 * xstar + 2 = 2 ^ 2 / (Real.sqrt 3 - 1) := by
    unfold xstar; rw [eq_div_iff hpos.ne']; linear_combination 2 * hs
  rw [this, Real.log_div (by positivity) hpos.ne', Real.log_pow]; push_cast; ring

lemma log_c7 : Real.log (8 * xstar - 4) = 2 * Real.log 2 + Real.log (Real.sqrt 3 - 1) := by
  have hpos := sqrt3_sub_one_pos
  have : 8 * xstar - 4 = 2 ^ 2 * (Real.sqrt 3 - 1) := by unfold xstar; ring
  rw [this, Real.log_mul (by norm_num) hpos.ne', Real.log_pow]; push_cast; ring

lemma log_c8 : Real.log (1 - xstar) = 2 * Real.log (Real.sqrt 3 - 1) - 2 * Real.log 2 := by
  have hs := sqrt3_sq
  have hpos := sqrt3_sub_one_pos
  have : 1 - xstar = (Real.sqrt 3 - 1) ^ 2 / 2 ^ 2 := by
    unfold xstar; rw [eq_div_iff (by norm_num)]; linear_combination -hs
  rw [this, Real.log_div (by positivity) (by norm_num), Real.log_pow, Real.log_pow]; push_cast; ring

lemma log_c9 : Real.log (1 + xstar) = -(2 * Real.log (Real.sqrt 3 - 1)) := by
  have hs := sqrt3_sq
  have hpos := sqrt3_sub_one_pos
  have : 1 + xstar = ((Real.sqrt 3 - 1) ^ 2)⁻¹ := by
    unfold xstar; apply eq_inv_of_mul_eq_one_left
    linear_combination (Real.sqrt 3 / 2) * hs
  rw [this, Real.log_inv, Real.log_pow]; push_cast; ring

lemma log_four : Real.log 4 = 2 * Real.log 2 := by
  rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; push_cast; ring

/-! ### Stirling main terms -/

/-- the Stirling main term of `log (cN+u)!` -/
noncomputable def St (c u : ℝ) (N : ℕ) : ℝ :=
  c * N * Real.log N + c * N * (Real.log c - 1) + u * Real.log N + u * Real.log c
    + 1 / 2 * Real.log (2 * Real.pi * c * N)

lemma St_eq (c u : ℝ) (N : ℕ) (hc : 0 < c) (hN : 1 ≤ N) :
    St c u N = c * N * Real.log N + c * N * (Real.log c - 1) + u * Real.log N + u * Real.log c
      + 1 / 2 * (Real.log 2 + Real.log Real.pi + Real.log c + Real.log N) := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hpi := Real.pi_pos
  unfold St
  rw [Real.log_mul (by positivity) hN'.ne', Real.log_mul (by positivity) hc.ne',
    Real.log_mul (by norm_num) hpi.ne']

/-- Stirling for `m = c N + (r − q θ_N)`. -/
lemma stir (c : ℝ) (hc : 0 < c) (m : ℕ → ℕ) (r q : ℝ)
    (hm : ∀ᶠ N in atTop, (m N : ℝ) = c * N + (r - q * θc N)) :
    Tendsto (fun N : ℕ => Real.log ((m N).factorial : ℝ) - St c (r - q * θc N) N) atTop (𝓝 0) := by
  refine tendsto_log_factorial_linear c hc m (fun N => r - q * θc N) (|r| + |q|) ?_ hm
  intro N
  have h0 := θc_nonneg N
  have h1 := θc_lt_one N
  calc |r - q * θc N| ≤ |r| + |q * θc N| := abs_sub _ _
    _ ≤ |r| + |q| := by
      rw [abs_mul, abs_of_nonneg h0]; nlinarith [abs_nonneg q]

/-- The signed sum of the eleven Stirling main terms. -/
lemma main_term (N : ℕ) (hN : 1 ≤ N) :
    St 4 (0 - 0 * θc N) N + St (8 * xstar + 4) (0 - 8 * θc N) N
        + St (4 * xstar - 2) (4 - 4 * θc N) N
        - St 2 (0 - 0 * θc N) N - St 2 (-2 - 0 * θc N) N - St (4 * xstar + 2) (0 - 4 * θc N) N
        - St (8 * xstar - 4) (8 - 8 * θc N) N - 2 * St (1 - xstar) (-1 - (-1) * θc N) N
        - 2 * St (1 + xstar) (0 - 1 * θc N) N
      = 22 * N * Real.log 2 - 12 * Real.log 2 - 5 / 2 * Real.log N
        - 5 / 2 * Real.log (2 * Real.pi) := by
  have hx1 := xstar_gt_C
  have hx2 := xstar_lt_C
  rw [St_eq 4 _ N (by norm_num) hN, St_eq (8 * xstar + 4) _ N (by linarith) hN,
    St_eq (4 * xstar - 2) _ N (by linarith) hN, St_eq 2 (0 - 0 * θc N) N (by norm_num) hN,
    St_eq 2 (-2 - 0 * θc N) N (by norm_num) hN, St_eq (4 * xstar + 2) _ N (by linarith) hN,
    St_eq (8 * xstar - 4) _ N (by linarith) hN, St_eq (1 - xstar) _ N (by linarith) hN,
    St_eq (1 + xstar) _ N (by linarith) hN, log_c2, log_c3, log_c6, log_c7, log_c8, log_c9,
    log_four, Real.log_mul two_ne_zero Real.pi_pos.ne']
  ring

/-! ### The factorial part of `cf N k0 · √N` -/

/-- `cf N k0 · √N` with the prefactor `(4N−1)(2k0+1)` replaced by `N²`. -/
noncomputable def Rf (N : ℕ) : ℝ :=
  2 ^ 12 * (N : ℝ) ^ 2 * Real.sqrt N * ((4 * N).factorial : ℝ)
      * ((8 * k0 N + 4 * N).factorial : ℝ) * ((4 * k0 N + 4 - 2 * N).factorial : ℝ) /
    (2 ^ (22 * N) * ((2 * N).factorial : ℝ) * ((2 * N - 2).factorial : ℝ)
      * ((4 * k0 N + 2 * N).factorial : ℝ) * ((8 * k0 N + 8 - 4 * N).factorial : ℝ)
      * ((N - k0 N - 1).factorial : ℝ) ^ 2 * ((N + k0 N).factorial : ℝ) ^ 2)

lemma cf_alg (n s A B C D k : ℝ) (hn : n ≠ 0) :
    2 ^ 12 * (4 * n - 1) * (2 * k + 1) * A * B * C / D * s
      = 2 ^ 12 * n ^ 2 * s * A * B * C / D * ((4 * n - 1) * (2 * k + 1) / n ^ 2) := by
  field_simp

lemma cf_mul_sqrt (N : ℕ) (hN : 1 ≤ N) :
    cf N (k0 N) * Real.sqrt N
      = Rf N * ((4 * (N : ℝ) - 1) * (2 * (k0 N : ℝ) + 1) / (N : ℝ) ^ 2) := by
  have hN' : (N : ℝ) ≠ 0 := by positivity
  unfold cf Rf
  exact cf_alg _ _ _ _ _ _ _ hN'

lemma log_Rf_alg (n f1 f2 f3 f4 f5 f6 f7 f8 f9 : ℝ) (M : ℕ) (hn : 0 < n) (h1 : 0 < f1)
    (h2 : 0 < f2) (h3 : 0 < f3) (h4 : 0 < f4) (h5 : 0 < f5) (h6 : 0 < f6) (h7 : 0 < f7)
    (h8 : 0 < f8) (h9 : 0 < f9) :
    Real.log (2 ^ 12 * n ^ 2 * Real.sqrt n * f1 * f2 * f3 /
        (2 ^ (22 * M) * f4 * f5 * f6 * f7 * f8 ^ 2 * f9 ^ 2))
      = 12 * Real.log 2 + 5 / 2 * Real.log n + Real.log f1 + Real.log f2 + Real.log f3
        - 22 * M * Real.log 2 - Real.log f4 - Real.log f5 - Real.log f6 - Real.log f7
        - 2 * Real.log f8 - 2 * Real.log f9 := by
  have hs : 0 < Real.sqrt n := Real.sqrt_pos.2 hn
  rw [Real.log_div (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_sqrt hn.le]
  simp only [Real.log_pow]
  push_cast
  ring

lemma Rf_pos (N : ℕ) (hN : 1 ≤ N) : 0 < Rf N := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hs : 0 < Real.sqrt N := Real.sqrt_pos.2 hN'
  unfold Rf
  positivity

lemma fact_pos (n : ℕ) : (0 : ℝ) < (n.factorial : ℝ) := by exact_mod_cast n.factorial_pos

lemma tendsto_log_Rf :
    Tendsto (fun N : ℕ => Real.log (Rf N)) atTop (𝓝 (-(5 / 2 * Real.log (2 * Real.pi)))) := by
  have hx1 := xstar_gt_C
  have hx2 := xstar_lt_C
  have e1 := stir 4 (by norm_num) (fun N => 4 * N) 0 0
    (Eventually.of_forall fun N => by push_cast; ring)
  have e2 := stir (8 * xstar + 4) (by linarith) (fun N => 8 * k0 N + 4 * N) 0 8
    (Eventually.of_forall fun N => by push_cast; rw [k0_cast]; ring)
  have e3 := stir (4 * xstar - 2) (by linarith) (fun N => 4 * k0 N + 4 - 2 * N) 4 4
    (Eventually.of_forall fun N => by
      rw [Nat.cast_sub (two_N_le N)]; push_cast; rw [k0_cast]; ring)
  have e4 := stir 2 (by norm_num) (fun N => 2 * N) 0 0
    (Eventually.of_forall fun N => by push_cast; ring)
  have e5 := stir 2 (by norm_num) (fun N => 2 * N - 2) (-2) 0 (by
    filter_upwards [eventually_ge_atTop 1] with N hN
    rw [Nat.cast_sub (by omega)]; push_cast; ring)
  have e6 := stir (4 * xstar + 2) (by linarith) (fun N => 4 * k0 N + 2 * N) 0 4
    (Eventually.of_forall fun N => by push_cast; rw [k0_cast]; ring)
  have e7 := stir (8 * xstar - 4) (by linarith) (fun N => 8 * k0 N + 8 - 4 * N) 8 8
    (Eventually.of_forall fun N => by
      rw [Nat.cast_sub (by have := two_N_le N; omega)]; push_cast; rw [k0_cast]; ring)
  have e8 := stir (1 - xstar) (by linarith) (fun N => N - k0 N - 1) (-1) (-1) (by
    filter_upwards [eventually_ge_atTop 1] with N hN
    have h := k0_lt N hN
    rw [Nat.cast_sub (by omega), Nat.cast_sub (by omega)]; push_cast; rw [k0_cast]; ring)
  have e9 := stir (1 + xstar) (by linarith) (fun N => N + k0 N) 0 1
    (Eventually.of_forall fun N => by push_cast; rw [k0_cast]; ring)
  have hErr := ((((((((e1.add e2).add e3).sub e4).sub e5).sub e6).sub e7).sub
    (e8.const_mul 2)).sub (e9.const_mul 2)).sub_const (5 / 2 * Real.log (2 * Real.pi))
  have e0 : ((0 : ℝ) + 0 + 0 - 0 - 0 - 0 - 0 - 2 * 0 - 2 * 0) - 5 / 2 * Real.log (2 * Real.pi)
      = -(5 / 2 * Real.log (2 * Real.pi)) := by ring
  rw [e0] at hErr
  refine hErr.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have hmain := main_term N hN
  rw [Rf, log_Rf_alg _ _ _ _ _ _ _ _ _ _ _ hN' (fact_pos _) (fact_pos _) (fact_pos _)
    (fact_pos _) (fact_pos _) (fact_pos _) (fact_pos _) (fact_pos _) (fact_pos _)]
  linarith

lemma tendsto_Rf :
    Tendsto Rf atTop (𝓝 (Real.exp (-(5 / 2 * Real.log (2 * Real.pi))))) := by
  have h := (Real.continuous_exp.tendsto _).comp tendsto_log_Rf
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  exact Real.exp_log (Rf_pos N hN)

/-! ### The prefactor -/

lemma tendsto_prefactor :
    Tendsto (fun N : ℕ => (4 * (N : ℝ) - 1) * (2 * (k0 N : ℝ) + 1) / (N : ℝ) ^ 2) atTop
      (𝓝 (4 * (2 * xstar))) := by
  have h1 : Tendsto (fun N : ℕ => 1 / (N : ℝ)) atTop (𝓝 0) := tendsto_one_div_atTop_nhds_zero_nat
  have hθ : Tendsto (fun N : ℕ => θc N / (N : ℝ)) atTop (𝓝 0) :=
    squeeze_zero (fun N => div_nonneg (θc_nonneg N) (Nat.cast_nonneg N))
      (fun N => div_le_div_of_nonneg_right (θc_lt_one N).le (Nat.cast_nonneg N)) h1
  have h := ((tendsto_const_nhds (x := (4 : ℝ))).sub h1).mul
    (((tendsto_const_nhds (x := 2 * xstar)).sub (hθ.const_mul 2)).add h1)
  have e0 : ((4 : ℝ) - 0) * (2 * xstar - 2 * 0 + 0) = 4 * (2 * xstar) := by ring
  rw [e0] at h
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  rw [k0_cast]
  field_simp

/-! ### The constant -/

lemma centre_const :
    Real.exp (-(5 / 2 * Real.log (2 * Real.pi))) * (4 * (2 * xstar))
      = Real.sqrt 3 / (Real.sqrt 2 * Real.pi ^ (5 / 2 : ℝ)) := by
  have hpi := Real.pi_pos
  have h2 : Real.exp (5 / 2 * Real.log 2) = 4 * Real.sqrt 2 := by
    have : 5 / 2 * Real.log 2 = Real.log (4 * Real.sqrt 2) := by
      rw [Real.log_mul (by norm_num) (by positivity), Real.log_sqrt (by norm_num), log_four]
      ring
    rw [this, Real.exp_log (by positivity)]
  have h3 : Real.exp (-(5 / 2 * Real.log (2 * Real.pi)))
      = (4 * Real.sqrt 2)⁻¹ * (Real.pi ^ (5 / 2 : ℝ))⁻¹ := by
    rw [← h2, Real.rpow_def_of_pos hpi, ← Real.exp_neg, ← Real.exp_neg, ← Real.exp_add,
      Real.log_mul two_ne_zero hpi.ne']
    congr 1; ring
  have hpow : 0 < Real.pi ^ (5 / 2 : ℝ) := Real.rpow_pos_of_pos hpi _
  have hs2 : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)
  rw [h3]; unfold xstar
  field_simp

/-! ### The centre term -/

theorem tendsto_cf_centre :
    Tendsto (fun N : ℕ => cf N ⌊xstar * N⌋₊ * Real.sqrt N) atTop
      (𝓝 (Real.sqrt 3 / (Real.sqrt 2 * Real.pi ^ (5 / 2 : ℝ)))) := by
  have h := tendsto_Rf.mul tendsto_prefactor
  rw [centre_const] at h
  refine h.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  exact (cf_mul_sqrt N hN).symm

end PiSeries.SunPi2
