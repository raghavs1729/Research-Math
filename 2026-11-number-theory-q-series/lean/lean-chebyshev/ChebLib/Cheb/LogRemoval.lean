/-
ChebLib/Cheb/LogRemoval.lean — removing the logarithmic weight from an Abel-summed sequence.

If `Φ(y) = Σ_n a_n (log n) e^{−n/y} ≤ −c √y` for all large `y` (with `|a_n| ≤ 1`, `a_n = 0` for
`n < 3`), then `G(x) = Σ_n a_n e^{−n/x} → −∞`.

Proof idea (positive kernel).  For `n ≥ 3`, `x > 0`:
  `e^{−n/x} = (log n) ∫_0^∞ n^{−s} e^{−n/x} ds`,  `n^{−s} = Γ(s)^{-1} ∫_0^∞ v^{s−1} e^{−n v} dv`,
so with `Φ_s(x) := Σ_n a_n (log n) n^{−s} e^{−n/x}`:
  `G(x) = ∫_0^∞ Φ_s(x) ds`,   `Φ_s(x) = Γ(s)^{-1} ∫_0^∞ v^{s−1} Φ(1/(v + 1/x)) dv`.
For all `y > 0`:  `Φ(y) ≤ C₀ e^{−2/y} − c √y · 1[y ≥ x₀]`  with
`C₀ = Σ_{n ≥ 3} (log n) e^{−(n−2)/x₀}`.  Hence for `x ≥ 2 x₀`:
  `Φ_s(x) ≤ C₀ 2^{−s} − c √(x/2) x^{−s} / Γ(s+1)`   (restrict the negative part to `v ≤ 1/x`),
  `G(x) ≤ C₀ / log 2 − c √(x/2) ∫_0^1 x^{−s} ds = C₀ / log 2 − c √(x/2) (1 − 1/x) / log x → −∞`
(using `Γ(s+1) ≤ 1` on `[0,1]`).

Formalisation notes.  We work with `w = 1/x`, `b_n = a_n log n` (so `|b_n| ≤ log n`), and replace the
exact evaluation of `∫_0^1 x^{−s} ds` by the cruder `∫_0^{1/4} x^{−s} ds ≥ x^{−1/4}/4`, and
`√y ≥ √x / 2` for `y ≥ x/2`; this gives `G(x) ≤ C₀ ∫_0^∞ 2^{−s} ds − (c/8) x^{1/4} → −∞`.
The helper lemmas live in the namespace `ChebLib.LogRemoval`:
* `integrable_tsum_of_summable_integral_norm` : integrability of a countable sum;
* `integral_rpow_mul_tsum` (interchange in `v`), `integral_tsum_rpow` (interchange in `s`);
* `tsum_le_C0` (bound `Φ(y) ≤ C₀ e^{−2/y}` for `y ≤ x₀`), `tsum_rpow_le` (bound for `Φ_s`),
  `tsum_div_log_le` (bound for `G`).
-/
import Mathlib

open MeasureTheory Real Filter Topology Set

noncomputable section

namespace ChebLib

namespace LogRemoval

/-- Integrability of a countable sum, given summability of the `L¹` norms. -/
lemma integrable_tsum_of_summable_integral_norm {α : Type*} [MeasurableSpace α] {μ : Measure α}
    {F : ℕ → α → ℝ} (hF_int : ∀ i, Integrable (F i) μ)
    (hF_sum : Summable fun i ↦ ∫ a, ‖F i a‖ ∂μ)
    (hsum : ∀ᵐ a ∂μ, Summable (fun i => F i a)) :
    Integrable (fun a => ∑' i, F i a) μ := by
  have hmeas : AEStronglyMeasurable (fun a => ∑' i, F i a) μ := by
    refine aestronglyMeasurable_of_tendsto_ae atTop
      (f := fun N a => ∑ i ∈ Finset.range N, F i a) ?_ ?_
    · intro N
      exact (integrable_finsetSum _ fun i _ => hF_int i).aestronglyMeasurable
    · filter_upwards [hsum] with a ha using ha.hasSum.tendsto_sum_nat
  refine ⟨hmeas, ?_⟩
  have h1 : ∫⁻ a, ‖∑' i, F i a‖ₑ ∂μ ≤ ∫⁻ a, ∑' i, ‖F i a‖ₑ ∂μ :=
    lintegral_mono fun a => enorm_tsum_le_tsum_enorm
  rw [lintegral_tsum (fun i => (hF_int i).aestronglyMeasurable.enorm)] at h1
  refine lt_of_le_of_lt h1 (lt_top_iff_ne_top.2 ?_)
  have (i : ℕ) : ∫⁻ a, ‖F i a‖ₑ ∂μ = ENNReal.ofReal (∫ a, ‖F i a‖ ∂μ) :=
    (ofReal_integral_norm_eq_lintegral_enorm (hF_int i)).symm
  rw [funext this]
  exact hF_sum.tsum_ofReal_ne_top

lemma summable_log_mul_exp {w : ℝ} (hw : 0 < w) :
    Summable (fun n : ℕ => Real.log n * Real.exp (-(n * w))) := by
  have h := summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) 1 (r := Real.exp (-w))
    (by rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]; exact Real.exp_lt_one_iff.2 (by linarith))
  refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_) h
  · exact mul_nonneg (Real.log_natCast_nonneg n) (Real.exp_pos _).le
  · rw [pow_one, ← Real.exp_nat_mul, mul_neg]
    exact mul_le_mul_of_nonneg_right (Real.log_le_self (Nat.cast_nonneg n)) (Real.exp_pos _).le

lemma summable_exp_neg_nat_mul {w : ℝ} (hw : 0 < w) :
    Summable (fun n : ℕ => Real.exp (-(n * w))) := by
  have h := summable_geometric_of_lt_one (r := Real.exp (-w)) (Real.exp_pos _).le
    (Real.exp_lt_one_iff.2 (by linarith))
  refine h.congr fun n => ?_
  rw [← Real.exp_nat_mul, mul_neg]

lemma summable_mul_exp {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {w : ℝ} (hw : 0 < w) :
    Summable (fun n : ℕ => b n * Real.exp (-(n * w))) := by
  refine Summable.of_norm_bounded (summable_log_mul_exp hw) fun n => ?_
  rw [Real.norm_eq_abs, abs_mul, abs_of_pos (Real.exp_pos _)]
  exact mul_le_mul_of_nonneg_right (hb n) (Real.exp_pos _).le

lemma summable_log_mul_exp_sub_two {w : ℝ} (hw : 0 < w) :
    Summable (fun n : ℕ => Real.log n * Real.exp (-(((n : ℝ) - 2) * w))) := by
  refine ((summable_log_mul_exp hw).mul_left (Real.exp (2 * w))).congr fun n => ?_
  rw [mul_left_comm, ← Real.exp_add]
  congr 2; ring

/-- The constant `C₀`. -/
def C0 (w : ℝ) : ℝ := ∑' n : ℕ, Real.log n * Real.exp (-(((n : ℝ) - 2) * w))

lemma C0_nonneg (w : ℝ) : 0 ≤ C0 w :=
  tsum_nonneg fun n => mul_nonneg (Real.log_natCast_nonneg n) (Real.exp_pos _).le

lemma tsum_le_C0 {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {w₀ u : ℝ} (hw₀ : 0 < w₀)
    (hu : w₀ ≤ u) :
    ∑' n : ℕ, b n * Real.exp (-(n * u)) ≤ C0 w₀ * Real.exp (-(2 * u)) := by
  have hupos : 0 < u := lt_of_lt_of_le hw₀ hu
  rw [C0, ← tsum_mul_right]
  refine (summable_mul_exp hb hupos).tsum_le_tsum (fun n => ?_)
    ((summable_log_mul_exp_sub_two hw₀).mul_right _)
  have h1 : b n ≤ Real.log n := (le_abs_self _).trans (hb n)
  have hlog := Real.log_natCast_nonneg n
  rw [mul_assoc, ← Real.exp_add]
  refine (mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le).trans ?_
  rcases Nat.lt_or_ge n 2 with hn | hn
  · have : Real.log (n : ℝ) = 0 := by
      interval_cases n <;> simp
    rw [this]; simp
  · refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 ?_) hlog
    have h2 : (2 : ℝ) ≤ n := by exact_mod_cast hn
    nlinarith [mul_le_mul_of_nonneg_left hu (sub_nonneg.2 h2)]

lemma eq_zero_of_abs_le_log {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {n : ℕ} (hn : n < 2) :
    b n = 0 := by
  have h := hb n
  have : Real.log (n : ℝ) = 0 := by interval_cases n <;> simp
  rw [this] at h
  exact abs_nonpos_iff.1 h

lemma one_div_natCast_rpow_le_one (n : ℕ) {s : ℝ} (hs : 0 ≤ s) : (1 / (n : ℝ)) ^ s ≤ 1 := by
  refine Real.rpow_le_one (by positivity) ?_ hs
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  · rw [div_le_one (by exact_mod_cast hn)]; exact_mod_cast hn

/-- Step A: the Mellin-type identity in the variable `v`. -/
lemma integral_rpow_mul_tsum {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {s w : ℝ} (hs : 0 < s)
    (hw : 0 < w) :
    IntegrableOn (fun v : ℝ => v ^ (s - 1) * ∑' n : ℕ, b n * Real.exp (-(n * (v + w)))) (Ioi 0) ∧
    ∫ v in Ioi 0, v ^ (s - 1) * ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) =
      Real.Gamma s * ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) := by
  set F : ℕ → ℝ → ℝ :=
    fun n v => b n * Real.exp (-(n * w)) * (v ^ (s - 1) * Real.exp (-(n * v))) with hF
  have hb0 : b 0 = 0 := eq_zero_of_abs_le_log hb (by norm_num)
  have hG : 0 < Real.Gamma s := Real.Gamma_pos_of_pos hs
  have hint_base : ∀ n : ℕ, 0 < n →
      IntegrableOn (fun v : ℝ => v ^ (s - 1) * Real.exp (-(n * v))) (Ioi 0) := by
    intro n hn
    apply Integrable.of_integral_ne_zero
    rw [Real.integral_rpow_mul_exp_neg_mul_Ioi hs (by exact_mod_cast hn)]
    positivity
  have hF_int : ∀ n, Integrable (F n) (volume.restrict (Ioi 0)) := by
    intro n
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp [hF, hb0]
    · exact (hint_base n hn).const_mul _
  have hF_norm : ∀ n, ∫ v in Ioi 0, ‖F n v‖ =
      |b n| * Real.exp (-(n * w)) * ((1 / (n : ℝ)) ^ s * Real.Gamma s) := by
    intro n
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp [hF, hb0]
    · rw [← Real.integral_rpow_mul_exp_neg_mul_Ioi hs (by exact_mod_cast hn : (0 : ℝ) < n),
        ← integral_const_mul]
      refine setIntegral_congr_fun measurableSet_Ioi fun v hv => ?_
      simp only [hF, Real.norm_eq_abs, abs_mul, abs_of_pos (Real.exp_pos _)]
      rw [abs_of_nonneg (Real.rpow_nonneg (le_of_lt hv) _)]
  have hF_val : ∀ n, ∫ v in Ioi 0, F n v =
      b n * Real.exp (-(n * w)) * ((1 / (n : ℝ)) ^ s * Real.Gamma s) := by
    intro n
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp [hF, hb0]
    · rw [← Real.integral_rpow_mul_exp_neg_mul_Ioi hs (by exact_mod_cast hn : (0 : ℝ) < n),
        ← integral_const_mul]
  have hF_sum : Summable fun n => ∫ v in Ioi 0, ‖F n v‖ := by
    simp_rw [hF_norm]
    refine Summable.of_nonneg_of_le (fun n => by positivity) (fun n => ?_)
      ((summable_log_mul_exp hw).mul_right (Real.Gamma s))
    have h1 : (1 / (n : ℝ)) ^ s ≤ 1 := one_div_natCast_rpow_le_one n hs.le
    have h2 : 0 ≤ (1 / (n : ℝ)) ^ s := by positivity
    have h3 : 0 < Real.exp (-(n * w)) := Real.exp_pos _
    have h4 := hb n
    have h5 := abs_nonneg (b n)
    calc |b n| * Real.exp (-(n * w)) * ((1 / (n : ℝ)) ^ s * Real.Gamma s)
        ≤ Real.log n * Real.exp (-(n * w)) * (1 * Real.Gamma s) := by gcongr
      _ = Real.log n * Real.exp (-(n * w)) * Real.Gamma s := by ring
  have hterm : ∀ (v : ℝ) (n : ℕ), F n v = v ^ (s - 1) * (b n * Real.exp (-(n * (v + w)))) := by
    intro v n
    simp only [hF]
    rw [mul_add, neg_add, Real.exp_add]; ring
  have hsum : ∀ v ∈ Ioi (0 : ℝ), Summable (fun n => F n v) := by
    intro v hv
    simp_rw [hterm v]
    exact (summable_mul_exp hb (add_pos hv hw)).mul_left _
  have hEq : ∀ v ∈ Ioi (0 : ℝ), ∑' n, F n v =
      v ^ (s - 1) * ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) := by
    intro v _
    simp_rw [hterm v]
    rw [tsum_mul_left]
  have hI := integrable_tsum_of_summable_integral_norm hF_int hF_sum
    (ae_restrict_of_forall_mem measurableSet_Ioi hsum)
  have hS := (hasSum_integral_of_summable_integral_norm hF_int hF_sum).tsum_eq
  constructor
  · exact hI.congr (ae_restrict_of_forall_mem measurableSet_Ioi hEq)
  · rw [← setIntegral_congr_fun measurableSet_Ioi hEq, ← hS]
    simp_rw [hF_val]
    rw [← tsum_mul_left]
    congr 1; ext n; ring

/-- Step B: integrating out `s`. -/
lemma integral_tsum_rpow {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {w : ℝ} (hw : 0 < w) :
    IntegrableOn (fun s : ℝ => ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w))) (Ioi 0) ∧
    ∫ s : ℝ in Ioi 0, ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) =
      ∑' n : ℕ, b n / Real.log n * Real.exp (-(n * w)) := by
  set H : ℕ → ℝ → ℝ := fun n s => b n * Real.exp (-(n * w)) * (1 / (n : ℝ)) ^ s with hH
  have hb01 : ∀ n, n < 2 → b n = 0 := fun n hn => eq_zero_of_abs_le_log hb hn
  have hlog : ∀ n : ℕ, 2 ≤ n → 0 < Real.log n := fun n hn => Real.log_pos (by exact_mod_cast hn)
  have hpow : ∀ n : ℕ, 2 ≤ n → ∀ s : ℝ, (1 / (n : ℝ)) ^ s = Real.exp (-Real.log n * s) := by
    intro n hn s
    have : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
    rw [Real.rpow_def_of_pos (by positivity), one_div, Real.log_inv]
  have hbase_int : ∀ n : ℕ, 2 ≤ n → IntegrableOn (fun s : ℝ => (1 / (n : ℝ)) ^ s) (Ioi 0) := by
    intro n hn
    have := integrableOn_exp_mul_Ioi (a := -Real.log n) (neg_lt_zero.2 (hlog n hn)) 0
    exact this.congr_fun (fun s _ => (hpow n hn s).symm) measurableSet_Ioi
  have hbase_val : ∀ n : ℕ, 2 ≤ n → ∫ s : ℝ in Ioi 0, (1 / (n : ℝ)) ^ s = 1 / Real.log n := by
    intro n hn
    rw [setIntegral_congr_fun measurableSet_Ioi (fun s _ => hpow n hn s),
      integral_exp_mul_Ioi (neg_lt_zero.2 (hlog n hn)) 0]
    simp [neg_div_neg_eq]
  have hH_int : ∀ n, Integrable (H n) (volume.restrict (Ioi 0)) := by
    intro n
    rcases Nat.lt_or_ge n 2 with hn | hn
    · simp [hH, hb01 n hn]
    · exact (hbase_int n hn).const_mul _
  have hH_norm : ∀ n, ∫ s in Ioi 0, ‖H n s‖ =
      |b n| * Real.exp (-(n * w)) * (1 / Real.log n) := by
    intro n
    rcases Nat.lt_or_ge n 2 with hn | hn
    · simp [hH, hb01 n hn]
    · rw [← hbase_val n hn, ← integral_const_mul]
      refine setIntegral_congr_fun measurableSet_Ioi fun s _ => ?_
      simp only [hH, Real.norm_eq_abs, abs_mul, abs_of_pos (Real.exp_pos _)]
      rw [abs_of_nonneg (Real.rpow_nonneg (by positivity) _)]
  have hH_val : ∀ n, ∫ s in Ioi 0, H n s =
      b n * Real.exp (-(n * w)) * (1 / Real.log n) := by
    intro n
    rcases Nat.lt_or_ge n 2 with hn | hn
    · simp [hH, hb01 n hn]
    · rw [← hbase_val n hn, ← integral_const_mul]
  have hH_sum : Summable fun n => ∫ s in Ioi 0, ‖H n s‖ := by
    simp_rw [hH_norm]
    refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_) (summable_exp_neg_nat_mul hw)
    · have := Real.log_natCast_nonneg n
      positivity
    · rcases Nat.lt_or_ge n 2 with hn | hn
      · simp [hb01 n hn, (Real.exp_pos _).le]
      · have h3 : 0 < Real.exp (-(n * w)) := Real.exp_pos _
        have h4 := hb n
        have h5 := hlog n hn
        rw [mul_one_div, div_le_iff₀ h5]
        nlinarith [mul_le_mul_of_nonneg_right h4 h3.le]
  have hterm : ∀ (s : ℝ) (n : ℕ), H n s = b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) := by
    intro s n
    simp only [hH]; ring
  have hsum : ∀ s ∈ Ioi (0 : ℝ), Summable (fun n => H n s) := by
    intro s hs
    refine Summable.of_norm_bounded (summable_log_mul_exp hw) fun n => ?_
    have h1 : (1 / (n : ℝ)) ^ s ≤ 1 := one_div_natCast_rpow_le_one n (le_of_lt hs)
    have h2 : 0 ≤ (1 / (n : ℝ)) ^ s := by positivity
    have h3 : 0 < Real.exp (-(n * w)) := Real.exp_pos _
    have h4 := hb n
    have h5 := abs_nonneg (b n)
    simp only [hH, Real.norm_eq_abs, abs_mul, abs_of_pos h3, abs_of_nonneg h2]
    calc |b n| * Real.exp (-(n * w)) * (1 / (n : ℝ)) ^ s
        ≤ Real.log n * Real.exp (-(n * w)) * 1 := by gcongr
      _ = Real.log n * Real.exp (-(n * w)) := by ring
  have hI := integrable_tsum_of_summable_integral_norm hH_int hH_sum
    (ae_restrict_of_forall_mem measurableSet_Ioi hsum)
  have hS := (hasSum_integral_of_summable_integral_norm hH_int hH_sum).tsum_eq
  constructor
  · exact hI.congr (ae_restrict_of_forall_mem measurableSet_Ioi fun s _ => by simp_rw [hterm s])
  · have : ∀ s ∈ Ioi (0 : ℝ), ∑' n, H n s =
        ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) :=
      fun s _ => by simp_rw [hterm s]
    rw [← setIntegral_congr_fun measurableSet_Ioi this, ← hS]
    simp_rw [hH_val]
    congr 1; ext n; ring

lemma Gamma_add_one_le_one {s : ℝ} (hs0 : 0 ≤ s) (hs1 : s ≤ 1) : Real.Gamma (s + 1) ≤ 1 := by
  have h := Real.convexOn_Gamma.le_max_of_mem_Icc (x := 1) (y := 2) (z := s + 1)
    (by simp) (by simp) ⟨by linarith, by linarith⟩
  simpa [Real.Gamma_one, Real.Gamma_two] using h

/-- Step (iv): bound for `Φ_s`. -/
lemma tsum_rpow_le {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {s w C M : ℝ} (hs : 0 < s)
    (hw : 0 < w) (hC : 0 ≤ C)
    (h1 : ∀ v : ℝ, 0 < v →
      ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) ≤ C * Real.exp (-(2 * v)))
    (h2 : ∀ v : ℝ, 0 < v → v ≤ w → ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) ≤ -M) :
    ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) ≤
      C * (1 / 2 : ℝ) ^ s - M * w ^ s / Real.Gamma (s + 1) := by
  obtain ⟨hI, hval⟩ := integral_rpow_mul_tsum hb hs hw
  have hG : 0 < Real.Gamma s := Real.Gamma_pos_of_pos hs
  have hg1 : IntegrableOn (fun v : ℝ => C * (v ^ (s - 1) * Real.exp (-(2 * v)))) (Ioi 0) := by
    have : IntegrableOn (fun v : ℝ => v ^ (s - 1) * Real.exp (-(2 * v))) (Ioi 0) := by
      apply Integrable.of_integral_ne_zero
      rw [Real.integral_rpow_mul_exp_neg_mul_Ioi hs two_pos]; positivity
    exact this.const_mul C
  have hg2 : IntegrableOn (fun v : ℝ => M * v ^ (s - 1)) (Ioc 0 w) := by
    have := intervalIntegral.intervalIntegrable_rpow' (a := 0) (b := w) (r := s - 1) (by linarith)
    exact ((intervalIntegrable_iff_integrableOn_Ioc_of_le hw.le).1 this).const_mul M
  have hg2' : IntegrableOn ((Ioc 0 w).indicator (fun v : ℝ => M * v ^ (s - 1))) (Ioi 0) :=
    (hg2.integrable_indicator measurableSet_Ioc).integrableOn
  have hval2 : ∫ v in Ioi 0, (Ioc 0 w).indicator (fun v : ℝ => M * v ^ (s - 1)) v =
      M * (w ^ s / s) := by
    rw [setIntegral_indicator measurableSet_Ioc,
      Set.inter_eq_right.2 (Ioc_subset_Ioi_self), integral_const_mul,
      ← intervalIntegral.integral_of_le hw.le, integral_rpow (Or.inl (by linarith))]
    simp [Real.zero_rpow hs.ne']
  have key : ∫ v in Ioi 0, v ^ (s - 1) * ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) ≤
      ∫ v in Ioi 0, (C * (v ^ (s - 1) * Real.exp (-(2 * v))) -
        (Ioc 0 w).indicator (fun v : ℝ => M * v ^ (s - 1)) v) := by
    refine setIntegral_mono_on hI (hg1.sub hg2') measurableSet_Ioi fun v hv => ?_
    have hv0 : (0 : ℝ) < v := hv
    have hvpow : 0 ≤ v ^ (s - 1) := Real.rpow_nonneg hv0.le _
    have hCe : 0 ≤ C * Real.exp (-(2 * v)) := mul_nonneg hC (Real.exp_pos _).le
    by_cases hvw : v ≤ w
    · rw [Set.indicator_of_mem (show v ∈ Ioc 0 w from ⟨hv0, hvw⟩)]
      have h3 : ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) ≤ C * Real.exp (-(2 * v)) - M := by
        linarith [h2 v hv0 hvw]
      calc v ^ (s - 1) * ∑' n : ℕ, b n * Real.exp (-(n * (v + w)))
          ≤ v ^ (s - 1) * (C * Real.exp (-(2 * v)) - M) := mul_le_mul_of_nonneg_left h3 hvpow
        _ = _ := by ring
    · rw [Set.indicator_of_notMem (fun h : v ∈ Ioc 0 w => hvw h.2)]
      calc v ^ (s - 1) * ∑' n : ℕ, b n * Real.exp (-(n * (v + w)))
          ≤ v ^ (s - 1) * (C * Real.exp (-(2 * v))) := mul_le_mul_of_nonneg_left (h1 v hv0) hvpow
        _ = _ := by ring
  rw [hval, integral_sub hg1 hg2', hval2, integral_const_mul,
    Real.integral_rpow_mul_exp_neg_mul_Ioi hs two_pos] at key
  have e : (C * (1 / 2 : ℝ) ^ s - M * w ^ s / (s * Real.Gamma s)) * Real.Gamma s =
      C * ((1 / 2 : ℝ) ^ s * Real.Gamma s) - M * (w ^ s / s) := by
    field_simp
  rw [Real.Gamma_add_one hs.ne']
  refine le_of_mul_le_mul_right ?_ hG
  rw [e]
  linarith [key]

/-- Step (v): integrating the bound for `Φ_s` over `s`. -/
lemma tsum_div_log_le {b : ℕ → ℝ} (hb : ∀ n, |b n| ≤ Real.log n) {w C M : ℝ} (hw : 0 < w)
    (hw1 : w ≤ 1) (hM : 0 ≤ M)
    (h : ∀ s : ℝ, 0 < s → ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) ≤
      C * (1 / 2 : ℝ) ^ s - M * w ^ s / Real.Gamma (s + 1)) :
    ∑' n : ℕ, b n / Real.log n * Real.exp (-(n * w)) ≤
      C * (∫ s : ℝ in Ioi 0, (1 / 2 : ℝ) ^ s) - M * w ^ (1 / 4 : ℝ) / 4 := by
  obtain ⟨hI, hval⟩ := integral_tsum_rpow hb hw
  rw [← hval]
  have hg1 : IntegrableOn (fun s : ℝ => C * (1 / 2 : ℝ) ^ s) (Ioi 0) := by
    have h0 : Real.log (1 / 2 : ℝ) < 0 := Real.log_neg (by norm_num) (by norm_num)
    have := integrableOn_exp_mul_Ioi h0 0
    have h2 : IntegrableOn (fun s : ℝ => (1 / 2 : ℝ) ^ s) (Ioi 0) :=
      this.congr_fun (fun s _ => (Real.rpow_def_of_pos (by norm_num) s).symm) measurableSet_Ioi
    exact h2.const_mul C
  have hg2 : IntegrableOn ((Ioc 0 (1 / 4 : ℝ)).indicator (fun _ : ℝ => M * w ^ (1 / 4 : ℝ)))
      (Ioi 0) := by
    have : IntegrableOn (fun _ : ℝ => M * w ^ (1 / 4 : ℝ)) (Ioc 0 (1 / 4 : ℝ)) :=
      integrableOn_const (by simp)
    exact (this.integrable_indicator measurableSet_Ioc).integrableOn
  have hval2 : ∫ s in Ioi 0, (Ioc 0 (1 / 4 : ℝ)).indicator (fun _ : ℝ => M * w ^ (1 / 4 : ℝ)) s =
      M * w ^ (1 / 4 : ℝ) / 4 := by
    rw [setIntegral_indicator measurableSet_Ioc,
      Set.inter_eq_right.2 (Ioc_subset_Ioi_self), setIntegral_const,
      Real.volume_real_Ioc_of_le (by norm_num)]
    simp; ring
  have key : ∫ s : ℝ in Ioi 0, ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) ≤
      ∫ s in Ioi 0, (C * (1 / 2 : ℝ) ^ s -
        (Ioc 0 (1 / 4 : ℝ)).indicator (fun _ : ℝ => M * w ^ (1 / 4 : ℝ)) s) := by
    refine setIntegral_mono_on hI (hg1.sub hg2) measurableSet_Ioi fun s hs => ?_
    have hs0 : (0 : ℝ) < s := hs
    refine (h s hs0).trans ?_
    have hGpos : 0 < Real.Gamma (s + 1) := Real.Gamma_pos_of_pos (by linarith)
    have hws : 0 ≤ w ^ s := Real.rpow_nonneg hw.le _
    by_cases hs4 : s ≤ 1 / 4
    · rw [Set.indicator_of_mem (show s ∈ Ioc 0 (1 / 4 : ℝ) from ⟨hs0, hs4⟩)]
      have hG1 : Real.Gamma (s + 1) ≤ 1 := Gamma_add_one_le_one hs0.le (by linarith)
      have h3 : w ^ (1 / 4 : ℝ) ≤ w ^ s := Real.rpow_le_rpow_of_exponent_ge hw hw1 hs4
      have h4 : w ^ s ≤ w ^ s / Real.Gamma (s + 1) := by
        rw [le_div_iff₀ hGpos]; nlinarith
      have h5 : M * w ^ (1 / 4 : ℝ) ≤ M * (w ^ s / Real.Gamma (s + 1)) :=
        mul_le_mul_of_nonneg_left (h3.trans h4) hM
      rw [mul_div_assoc]; linarith
    · rw [Set.indicator_of_notMem (fun h : s ∈ Ioc 0 (1 / 4 : ℝ) => hs4 h.2)]
      have : 0 ≤ M * w ^ s / Real.Gamma (s + 1) := div_nonneg (mul_nonneg hM hws) hGpos.le
      linarith
  rw [integral_sub hg1 hg2, hval2, integral_const_mul] at key
  exact key

end LogRemoval

open LogRemoval

/-- **Removing the logarithmic weight.**  Let `a : ℕ → ℝ` be bounded by `1` and vanish for `n < 3`.
If the log-weighted Abel sums satisfy `Σ_n a_n (log n) e^{−n/y} ≤ −c √y` for all `y ≥ x₀` (some
`c > 0`), then the unweighted Abel sums `Σ_n a_n e^{−n/x}` tend to `−∞`. -/
theorem tendsto_atBot_of_log_weighted_le (a : ℕ → ℝ) (ha : ∀ n, |a n| ≤ 1)
    (ha3 : ∀ n, n < 3 → a n = 0) (c x₀ : ℝ) (hc : 0 < c)
    (hΦ : ∀ y : ℝ, x₀ ≤ y →
      ∑' n : ℕ, a n * Real.log n * Real.exp (-(n / y)) ≤ -c * Real.sqrt y) :
    Tendsto (fun x : ℝ => ∑' n : ℕ, a n * Real.exp (-(n / x))) atTop atBot := by
  set b : ℕ → ℝ := fun n => a n * Real.log n with hbdef
  have hb : ∀ n, |b n| ≤ Real.log n := by
    intro n
    have := Real.log_natCast_nonneg n
    simp only [hbdef, abs_mul, abs_of_nonneg this]
    exact mul_le_of_le_one_left this (ha n)
  set x₁ : ℝ := max x₀ 1 with hx₁
  have hx₁1 : 1 ≤ x₁ := le_max_right _ _
  have hx₁0 : x₀ ≤ x₁ := le_max_left _ _
  have hx₁pos : 0 < x₁ := by linarith
  set w₀ : ℝ := x₁⁻¹ with hw₀def
  have hw₀ : 0 < w₀ := inv_pos.2 hx₁pos
  set K : ℝ := ∫ s : ℝ in Ioi 0, (1 / 2 : ℝ) ^ s with hK
  have hΦ' : ∀ u : ℝ, 0 < u → u ≤ w₀ →
      ∑' n : ℕ, b n * Real.exp (-(n * u)) ≤ -c * Real.sqrt u⁻¹ := by
    intro u hu huw
    have h1 : x₀ ≤ u⁻¹ := by
      have : x₁ ≤ u⁻¹ := by rw [le_inv_comm₀ hx₁pos hu]; exact huw
      linarith
    have := hΦ u⁻¹ h1
    simpa only [div_inv_eq_mul, hbdef] using this
  have main : ∀ x : ℝ, max (2 * x₁) 1 ≤ x →
      ∑' n : ℕ, a n * Real.exp (-(n / x)) ≤ C0 w₀ * K - c / 8 * x ^ (1 / 4 : ℝ) := by
    intro x hx
    have hx2 : 2 * x₁ ≤ x := le_trans (le_max_left _ _) hx
    have hx1 : 1 ≤ x := le_trans (le_max_right _ _) hx
    have hxpos : 0 < x := by linarith
    set w : ℝ := x⁻¹ with hwdef
    have hw : 0 < w := inv_pos.2 hxpos
    have hw1 : w ≤ 1 := inv_le_one_of_one_le₀ hx1
    have hwx : w * x = 1 := inv_mul_cancel₀ hxpos.ne'
    have h2w : 2 * w ≤ w₀ := by
      have h1 : w₀ * x₁ = 1 := inv_mul_cancel₀ hx₁pos.ne'
      have h2 : 2 * w * x₁ ≤ w * x := by nlinarith
      have h3 : (2 * w) * x₁ ≤ w₀ * x₁ := by linarith
      exact le_of_mul_le_mul_right h3 hx₁pos
    set M : ℝ := c * (Real.sqrt x / 2) with hM
    have hMnn : 0 ≤ M := by positivity
    have hP1 : ∀ v : ℝ, 0 < v →
        ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) ≤ C0 w₀ * Real.exp (-(2 * v)) := by
      intro v hv
      rcases le_total (v + w) w₀ with h | h
      · have := hΦ' (v + w) (by linarith) h
        have h3 : 0 ≤ c * Real.sqrt (v + w)⁻¹ := by positivity
        have h4 : 0 ≤ C0 w₀ * Real.exp (-(2 * v)) :=
          mul_nonneg (C0_nonneg _) (Real.exp_pos _).le
        linarith
      · refine (tsum_le_C0 hb hw₀ h).trans ?_
        exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.2 (by linarith)) (C0_nonneg _)
    have hP2 : ∀ v : ℝ, 0 < v → v ≤ w →
        ∑' n : ℕ, b n * Real.exp (-(n * (v + w))) ≤ -M := by
      intro v hv hvw
      have := hΦ' (v + w) (by linarith) (by linarith)
      have h3 : Real.sqrt x / 2 ≤ Real.sqrt (v + w)⁻¹ := by
        rw [Real.le_sqrt (by positivity) (by positivity), div_pow, Real.sq_sqrt hxpos.le,
          le_inv_comm₀ (by positivity) (by linarith)]
        have : (x / 2 ^ 2)⁻¹ = 4 * w := by rw [hwdef]; field_simp; norm_num
        rw [this]; linarith
      have := mul_le_mul_of_nonneg_left h3 hc.le
      rw [hM]; linarith
    have hΦs : ∀ s : ℝ, 0 < s → ∑' n : ℕ, b n * (1 / (n : ℝ)) ^ s * Real.exp (-(n * w)) ≤
        C0 w₀ * (1 / 2 : ℝ) ^ s - M * w ^ s / Real.Gamma (s + 1) :=
      fun s hs => tsum_rpow_le hb hs hw (C0_nonneg _) hP1 hP2
    have hfin := tsum_div_log_le hb hw hw1 hMnn hΦs
    have hG : ∑' n : ℕ, a n * Real.exp (-(n / x)) =
        ∑' n : ℕ, b n / Real.log n * Real.exp (-(n * w)) := by
      refine tsum_congr fun n => ?_
      rw [hwdef, ← div_eq_mul_inv]
      rcases Nat.lt_or_ge n 3 with hn | hn
      · simp [hbdef, ha3 n hn]
      · have : Real.log n ≠ 0 :=
          (Real.log_pos (by exact_mod_cast (by omega : 1 < n))).ne'
        simp only [hbdef]; rw [mul_div_assoc, div_self this, mul_one]
    rw [hG]
    refine hfin.trans (le_of_eq ?_)
    have e1 : Real.sqrt x = x ^ (1 / 4 : ℝ) * x ^ (1 / 4 : ℝ) := by
      rw [Real.sqrt_eq_rpow, ← Real.rpow_add hxpos]; norm_num
    have e2 : w ^ (1 / 4 : ℝ) = (x ^ (1 / 4 : ℝ))⁻¹ := by rw [hwdef, Real.inv_rpow hxpos.le]
    have e3 : 0 < x ^ (1 / 4 : ℝ) := Real.rpow_pos_of_pos hxpos _
    rw [hM, e1, e2]; field_simp; ring
  have h1 : Tendsto (fun x : ℝ => c / 8 * x ^ (1 / 4 : ℝ)) atTop atTop :=
    (tendsto_rpow_atTop (by norm_num)).const_mul_atTop (by positivity)
  rw [tendsto_atBot]
  intro B
  filter_upwards [eventually_ge_atTop (max (2 * x₁) 1), tendsto_atTop.1 h1 (C0 w₀ * K - B)]
    with x hx hx'
  linarith [main x hx]

end ChebLib
