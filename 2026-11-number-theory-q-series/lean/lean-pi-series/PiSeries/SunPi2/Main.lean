/-
# Sun's series  `∑_{k≥0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^{22k}) = 64√2/π²`

Conjectured by Zhi-Wei Sun (arXiv:2307.03086, Conjecture 5.33(i)).  Proof: Theorem 11 of
`monotonicity/ramanujan-series/RESULTS.md`, reorganised as in `N2-PLAN.md`:

* `PartA.lean`   : the WZ pair on the triangle `k < n` and the exact identity `∑_{n<N} t n = 256 T N`,
                   `T N = ∑_{k<N} a N k`; the factorial closed form of `a N k` near the diagonal;
* `Tails.lean`   : `a N k ≤ B qᴺ` for `|k − x*N| ≥ δN` (`x* = √3/2`), from the ratio in `n`;
* `Ratio.lean`   : the Gaussian profile of `k ↦ a N k` near `x*N`, from the ratio in `k`;
* `Gaussian.lean`: `(1/√N) ∑_{|j|≤δN} exp(−c j²/N) → √(π/c)`;
* `Stirling.lean`, `Centre.lean`: `a N ⌊x*N⌋ √N → √3/(√2 π^{5/2})`;
* `Assembly.lean`: these imply `T N → √3/(√2 π^{5/2}) · √(π/12)`.
-/
import PiSeries.SunPi2.Assembly
import PiSeries.SunPi2.PartA
import PiSeries.SunPi2.Gaussian
import PiSeries.SunPi2.Centre
import PiSeries.SunPi2.Ratio
import PiSeries.SunPi2.Tails

open Finset Filter Topology

namespace PiSeries.SunPi2

/-- `√3/(√2 π^{5/2}) · √(π/12) = √2/(4π²)` -/
lemma const_eq :
    Real.sqrt 3 / (Real.sqrt 2 * Real.pi ^ (5 / 2 : ℝ)) * Real.sqrt (Real.pi / 12)
      = Real.sqrt 2 / (4 * Real.pi ^ 2) := by
  have hpi : 0 < Real.pi := Real.pi_pos
  have h1 : Real.pi ^ (5 / 2 : ℝ) = Real.pi ^ 2 * Real.sqrt Real.pi := by
    rw [Real.sqrt_eq_rpow, ← Real.rpow_natCast, ← Real.rpow_add hpi]
    norm_num
  have h2 : Real.sqrt (Real.pi / 12) = Real.sqrt Real.pi / (2 * Real.sqrt 3) := by
    rw [Real.sqrt_div hpi.le, show (12 : ℝ) = 2 ^ 2 * 3 by norm_num,
      Real.sqrt_mul (by norm_num), Real.sqrt_sq (by norm_num)]
  rw [h1, h2]
  have hs2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
  have h3 : 0 < Real.sqrt 3 := Real.sqrt_pos.mpr (by norm_num)
  have h4 : 0 < Real.sqrt Real.pi := Real.sqrt_pos.mpr hpi
  have h5 : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  rw [div_mul_div_comm, div_eq_div_iff (by positivity) (by positivity)]
  linear_combination (-(2 * Real.sqrt 3 * Real.sqrt Real.pi * Real.pi ^ 2)) * hs2

/-- **`T N → √2/(4π²)`** -/
theorem tendsto_T : Tendsto T atTop (𝓝 (Real.sqrt 2 / (4 * Real.pi ^ 2))) := by
  have hc : Tendsto (fun N : ℕ => a N ⌊xstar * N⌋₊ * Real.sqrt N) atTop
      (𝓝 (Real.sqrt 3 / (Real.sqrt 2 * Real.pi ^ (5 / 2 : ℝ)))) := by
    refine tendsto_cf_centre.congr' ?_
    filter_upwards [eventually_ge_atTop 8] with N hN
    have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
    have hxg := xstar_gt_086
    have hxl := xstar_lt_087
    have hk0le : (⌊xstar * N⌋₊ : ℝ) ≤ xstar * N := Nat.floor_le (by nlinarith)
    have hk0lt : xstar * N < (⌊xstar * N⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one _
    have h1 : ⌊xstar * N⌋₊ < N := by
      have : (⌊xstar * N⌋₊ : ℝ) < N := by nlinarith
      exact_mod_cast this
    have h2 : N ≤ 2 * ⌊xstar * N⌋₊ + 2 := by
      have : (N : ℝ) ≤ 2 * (⌊xstar * N⌋₊ : ℝ) + 2 := by nlinarith
      exact_mod_cast this
    rw [a_eq_cf h1 h2]
  have hB : 0 ≤ (∑' n, t n) / 256 :=
    div_nonneg (tsum_nonneg (fun n => (t_pos n).le)) (by norm_num)
  have h := tendsto_T_of ((∑' n, t n) / 256) _ hB (fun n k hn hk => a_pos hn hk)
    (fun δ hδ hδ' => tail_bound a _ (fun n k hn hk => a_pos hn hk)
      (fun n k hn hk => (a_le_T hn hk).trans (T_le n hn))
      (fun n k hn hk => a_succ_left hn hk) δ hδ hδ')
    prod_sigma_window tendsto_gaussian_sum hc
  rw [const_eq] at h
  exact h

/-- **Sun's conjecture 5.33(i) of arXiv:2307.03086:**
`∑_{k≥0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^{22k}) = 64√2/π²`. -/
theorem sun_series_64_sqrt2_div_pi_sq :
    HasSum (fun k : ℕ => (1344 * (k : ℝ) ^ 3 + 944 * (k : ℝ) ^ 2 + 156 * k + 9)
        * (Nat.choose (2 * k) k : ℝ) ^ 2 * (Nat.choose (4 * k) (2 * k) : ℝ) ^ 3
        / ((2 * (k : ℝ) + 1) * 2 ^ (22 * k)))
      (64 * Real.sqrt 2 / Real.pi ^ 2) := by
  have h := hasSum_t_of_tendsto_T tendsto_T
  have e : 256 * (Real.sqrt 2 / (4 * Real.pi ^ 2)) = 64 * Real.sqrt 2 / Real.pi ^ 2 := by
    have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
    field_simp
    ring
  rw [e] at h
  exact h

end PiSeries.SunPi2
