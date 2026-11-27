/-
# WP-S: Stirling's formula along linear sequences `m N = c N + u N` (`u` bounded), in logarithmic form.
-/
import Mathlib

open Filter Topology

namespace PiSeries.SunPi2

/-- `log n! = log (stirlingSeq n) + ½ log (2n) + n (log n − 1)`. -/
lemma log_factorial_eq_stirling (n : ℕ) :
    Real.log (n.factorial : ℝ)
      = Real.log (Stirling.stirlingSeq n) + 1 / 2 * Real.log (2 * n) + n * (Real.log n - 1) := by
  rw [Stirling.log_stirlingSeq_formula]
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  · have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    rw [Real.log_div hn'.ne' (Real.exp_pos 1).ne', Real.log_exp]; ring

theorem tendsto_log_factorial_linear (c : ℝ) (hc : 0 < c) (m : ℕ → ℕ) (u : ℕ → ℝ) (B : ℝ)
    (hu : ∀ N, |u N| ≤ B) (hm : ∀ᶠ N in atTop, (m N : ℝ) = c * N + u N) :
    Tendsto (fun N : ℕ => Real.log ((m N).factorial : ℝ)
        - (c * N * Real.log N + c * N * (Real.log c - 1) + u N * Real.log N + u N * Real.log c
           + 1 / 2 * Real.log (2 * Real.pi * c * N)))
      atTop (𝓝 0) := by
  have hB : 0 ≤ B := (abs_nonneg _).trans (hu 0)
  have hNtop : Tendsto (fun N : ℕ => (N : ℝ)) atTop atTop := tendsto_natCast_atTop_atTop
  have hcN : Tendsto (fun N : ℕ => c * (N : ℝ)) atTop atTop := hNtop.const_mul_atTop hc
  -- `m N → ∞`
  have hmR : Tendsto (fun N : ℕ => (m N : ℝ)) atTop atTop := by
    have h1 : Tendsto (fun N : ℕ => c * (N : ℝ) + (-B)) atTop atTop :=
      tendsto_atTop_add_const_right _ (-B) hcN
    refine tendsto_atTop_mono' _ ?_ h1
    filter_upwards [hm] with N hN
    rw [hN]; linarith [neg_abs_le (u N), hu N]
  have hmN : Tendsto m atTop atTop := tendsto_natCast_atTop_iff.mp hmR
  have h1 : Tendsto (fun N => Real.log (Stirling.stirlingSeq (m N))) atTop
      (𝓝 (Real.log (Real.sqrt Real.pi))) :=
    ((Real.continuousAt_log (by positivity)).tendsto.comp
      Stirling.tendsto_stirlingSeq_sqrt_pi).comp hmN
  -- `m / (cN) → 1`
  have hratio : Tendsto (fun N : ℕ => (m N : ℝ) / (c * N)) atTop (𝓝 1) := by
    have h0 : Tendsto (fun N : ℕ => u N / (c * N)) atTop (𝓝 0) := by
      refine squeeze_zero_norm' (a := fun N : ℕ => B / (c * N)) ?_
        (tendsto_const_nhds.div_atTop hcN)
      filter_upwards [eventually_gt_atTop 0] with N hN
      have hN' : (0 : ℝ) < N := by exact_mod_cast hN
      rw [Real.norm_eq_abs, abs_div, abs_of_pos (by positivity : (0 : ℝ) < c * N)]
      gcongr; exact hu N
    have : Tendsto (fun N : ℕ => 1 + u N / (c * N)) atTop (𝓝 (1 + 0)) :=
      tendsto_const_nhds.add h0
    rw [add_zero] at this
    refine this.congr' ?_
    filter_upwards [hm, eventually_gt_atTop 0] with N hN hpos
    have hN' : (0 : ℝ) < N := by exact_mod_cast hpos
    have : (0 : ℝ) < c * N := by positivity
    rw [hN]; field_simp
  have h2 : Tendsto (fun N : ℕ => Real.log ((m N : ℝ) / (c * N))) atTop (𝓝 0) := by
    have := (Real.continuousAt_log one_ne_zero).tendsto.comp hratio
    simpa [Function.comp_def] using this
  -- `m log (m/(cN)) − u → 0`
  have h3 : Tendsto (fun N : ℕ => (m N : ℝ) * Real.log ((m N : ℝ) / (c * N)) - u N) atTop
      (𝓝 0) := by
    have hup : Tendsto (fun N : ℕ => B ^ 2 / (c * N)) atTop (𝓝 0) :=
      tendsto_const_nhds.div_atTop hcN
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hup ?_ ?_
    · filter_upwards [hm, eventually_gt_atTop 0, hmR.eventually_gt_atTop 0] with N hN hpos hmpos
      have hN' : (0 : ℝ) < N := by exact_mod_cast hpos
      have hy : (0 : ℝ) < c * N := by positivity
      have h := Real.one_sub_inv_le_log_of_pos (div_pos hmpos hy)
      have h' : (m N : ℝ) * (1 - ((m N : ℝ) / (c * N))⁻¹) ≤
          (m N : ℝ) * Real.log ((m N : ℝ) / (c * N)) :=
        mul_le_mul_of_nonneg_left h hmpos.le
      have e : (m N : ℝ) * (1 - ((m N : ℝ) / (c * N))⁻¹) = u N := by
        rw [inv_div]; field_simp; linarith
      linarith
    · filter_upwards [hm, eventually_gt_atTop 0, hmR.eventually_gt_atTop 0] with N hN hpos hmpos
      have hN' : (0 : ℝ) < N := by exact_mod_cast hpos
      have hy : (0 : ℝ) < c * N := by positivity
      have h := Real.log_le_sub_one_of_pos (div_pos hmpos hy)
      have h' : (m N : ℝ) * Real.log ((m N : ℝ) / (c * N)) ≤
          (m N : ℝ) * ((m N : ℝ) / (c * N) - 1) :=
        mul_le_mul_of_nonneg_left h hmpos.le
      have e : (m N : ℝ) * ((m N : ℝ) / (c * N) - 1) - u N = (u N) ^ 2 / (c * N) := by
        rw [hN]; field_simp; ring
      have hsq : (u N) ^ 2 ≤ B ^ 2 := by
        have := hu N
        rw [← sq_abs]; gcongr
      have : (u N) ^ 2 / (c * N) ≤ B ^ 2 / (c * N) := by gcongr
      linarith
  have hsum := (h1.sub (tendsto_const_nhds (x := Real.log (Real.sqrt Real.pi)))).add
    ((h2.const_mul (1 / 2)).add h3)
  simp only [sub_self, mul_zero, add_zero] at hsum
  refine hsum.congr' ?_
  filter_upwards [hm, eventually_gt_atTop 0, hmR.eventually_gt_atTop 0] with N hN hpos hmpos
  have hN' : (0 : ℝ) < N := by exact_mod_cast hpos
  have hy : (0 : ℝ) < c * N := by positivity
  have hpi := Real.pi_pos
  rw [log_factorial_eq_stirling, Real.log_sqrt hpi.le,
    Real.log_div hmpos.ne' hy.ne', Real.log_mul hc.ne' hN'.ne',
    Real.log_mul (by positivity) hmpos.ne',
    show 2 * Real.pi * c * (N : ℝ) = 2 * (Real.pi * (c * N)) by ring,
    Real.log_mul (by positivity) (by positivity), Real.log_mul hpi.ne' hy.ne',
    Real.log_mul hc.ne' hN'.ne']
  have hu' : u N = (m N : ℝ) - c * N := by linarith
  rw [hu']
  ring

end PiSeries.SunPi2
