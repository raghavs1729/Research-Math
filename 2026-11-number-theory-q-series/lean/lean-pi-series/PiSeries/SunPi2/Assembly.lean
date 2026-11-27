/-
# Sun's `64√2/π²` series: assembly of the Laplace argument

`tendsto_T_of` derives `T N → ψ √(π/12)` from the four analytic inputs (tails, window profile,
Gaussian sum, value at the centre), stated here as hypotheses; `Main.lean` supplies them.
-/
import PiSeries.SunPi2.Defs

open Finset Filter Topology

namespace PiSeries.SunPi2

lemma xstar_gt_086 : (0.86 : ℝ) < xstar := by
  unfold xstar
  have h : (1.72 : ℝ) < Real.sqrt 3 := by
    rw [show (1.72 : ℝ) = Real.sqrt (1.72 ^ 2) by rw [Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  linarith

lemma xstar_lt_087 : xstar < (0.87 : ℝ) := by
  unfold xstar
  have h : Real.sqrt 3 < (1.74 : ℝ) := by
    rw [show (1.74 : ℝ) = Real.sqrt (1.74 ^ 2) by rw [Real.sqrt_sq (by norm_num)]]
    exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)
  linarith

lemma a_eq_mul_prod_Ico (N : ℕ) {k0 k : ℕ} (h : k0 ≤ k) :
    a N k = a N k0 * ∏ i ∈ Ico k0 k, σ N i := by
  unfold a
  rw [mul_assoc, Finset.prod_range_mul_prod_Ico _ h]

/-- from the two-sided bounds on the products of `σ` to two-sided bounds on `a N k / a N k0` -/
lemma window_pointwise {N k0 k : ℕ} {ε C : ℝ} (hak : 0 < a N k) (hak0 : 0 < a N k0)
    (h1 : k0 ≤ k →
      Real.exp (-(12 + ε) * ((k : ℝ) - k0) ^ 2 / N - C / N) ≤ ∏ i ∈ Ico k0 k, σ N i ∧
      ∏ i ∈ Ico k0 k, σ N i ≤ Real.exp (-(12 - ε) * ((k : ℝ) - k0) ^ 2 / N + C / N))
    (h2 : k ≤ k0 →
      Real.exp ((12 - ε) * ((k : ℝ) - k0) ^ 2 / N - C / N) ≤ ∏ i ∈ Ico k k0, σ N i ∧
      ∏ i ∈ Ico k k0, σ N i ≤ Real.exp ((12 + ε) * ((k : ℝ) - k0) ^ 2 / N + C / N)) :
    a N k0 * Real.exp (-(12 + ε) * ((k : ℝ) - k0) ^ 2 / N - C / N) ≤ a N k ∧
    a N k ≤ a N k0 * Real.exp (-(12 - ε) * ((k : ℝ) - k0) ^ 2 / N + C / N) := by
  rcases le_total k0 k with h | h
  · obtain ⟨l, u⟩ := h1 h
    rw [a_eq_mul_prod_Ico N h]
    exact ⟨mul_le_mul_of_nonneg_left l hak0.le, mul_le_mul_of_nonneg_left u hak0.le⟩
  · obtain ⟨l, u⟩ := h2 h
    have e := a_eq_mul_prod_Ico N h
    constructor
    · have h3 : a N k0 ≤ a N k * Real.exp ((12 + ε) * ((k : ℝ) - k0) ^ 2 / N + C / N) := by
        rw [e]; exact mul_le_mul_of_nonneg_left u hak.le
      have h4 : -(12 + ε) * ((k : ℝ) - k0) ^ 2 / N - C / N
          = -((12 + ε) * ((k : ℝ) - k0) ^ 2 / N + C / N) := by ring
      rw [h4, Real.exp_neg]
      rw [← div_eq_mul_inv, div_le_iff₀ (Real.exp_pos _)]
      exact h3
    · have h3 : a N k * Real.exp ((12 - ε) * ((k : ℝ) - k0) ^ 2 / N - C / N) ≤ a N k0 := by
        rw [e]; exact mul_le_mul_of_nonneg_left l hak.le
      have h4 : -(12 - ε) * ((k : ℝ) - k0) ^ 2 / N + C / N
          = -((12 - ε) * ((k : ℝ) - k0) ^ 2 / N - C / N) := by ring
      rw [h4, Real.exp_neg]
      rw [← div_eq_mul_inv, le_div_iff₀ (Real.exp_pos _)]
      exact h3

/-- the splitting of `T N` into window and tails, for one `N` -/
lemma T_bounds {N k0 M : ℕ} {ε C B q : ℝ} (hM : M ≤ k0) (hk0 : k0 + M + 1 ≤ N)
    (hposN : ∀ k, k < N → 0 < a N k)
    (hw : ∀ k : ℕ, k0 - M ≤ k → k ≤ k0 + M →
      a N k0 * Real.exp (-(12 + ε) * ((k : ℝ) - k0) ^ 2 / N - C / N) ≤ a N k ∧
      a N k ≤ a N k0 * Real.exp (-(12 - ε) * ((k : ℝ) - k0) ^ 2 / N + C / N))
    (ht : ∀ k : ℕ, k < N → (k < k0 - M ∨ k0 + M < k) → a N k ≤ B * q ^ N)
    (hBq : 0 ≤ B * q ^ N) :
    a N k0 * Real.exp (-C / N) * ∑ i ∈ range (2 * M + 1), Real.exp (-(12 + ε) * ((i : ℝ) - M) ^ 2 / N)
      ≤ T N ∧
    T N ≤ a N k0 * Real.exp (C / N) * ∑ i ∈ range (2 * M + 1), Real.exp (-(12 - ε) * ((i : ℝ) - M) ^ 2 / N)
      + N * (B * q ^ N) := by
  have hsplit : T N = ∑ k ∈ range (k0 - M), a N k + ∑ i ∈ range (2 * M + 1), a N (k0 - M + i)
      + ∑ k ∈ Ico (k0 + M + 1) N, a N k := by
    unfold T
    rw [← Finset.sum_range_add_sum_Ico _ (show k0 - M ≤ N by omega),
      ← Finset.sum_Ico_consecutive _ (show k0 - M ≤ k0 + M + 1 by omega) hk0,
      Finset.sum_Ico_eq_sum_range, show k0 + M + 1 - (k0 - M) = 2 * M + 1 by omega]
    ring
  have hcast : ∀ i : ℕ, (((k0 - M + i : ℕ) : ℝ)) - k0 = (i : ℝ) - M := by
    intro i
    push_cast [Nat.cast_sub hM]
    ring
  have hlow1 : 0 ≤ ∑ k ∈ range (k0 - M), a N k :=
    Finset.sum_nonneg (fun k hk => (hposN k (by have := Finset.mem_range.mp hk; omega)).le)
  have hlow2 : 0 ≤ ∑ k ∈ Ico (k0 + M + 1) N, a N k :=
    Finset.sum_nonneg (fun k hk => (hposN k (Finset.mem_Ico.mp hk).2).le)
  have hup1 : ∑ k ∈ range (k0 - M), a N k ≤ ((k0 - M : ℕ) : ℝ) * (B * q ^ N) := by
    have := Finset.sum_le_card_nsmul (range (k0 - M)) (fun k => a N k) (B * q ^ N)
      (fun k hk => ht k (by have := Finset.mem_range.mp hk; omega)
        (Or.inl (Finset.mem_range.mp hk)))
    simpa using this
  have hup2 : ∑ k ∈ Ico (k0 + M + 1) N, a N k ≤ ((N - (k0 + M + 1) : ℕ) : ℝ) * (B * q ^ N) := by
    have := Finset.sum_le_card_nsmul (Ico (k0 + M + 1) N) (fun k => a N k) (B * q ^ N)
      (fun k hk => ht k (Finset.mem_Ico.mp hk).2
        (Or.inr (by have := (Finset.mem_Ico.mp hk).1; omega)))
    simpa using this
  have hcard : ((k0 - M : ℕ) : ℝ) + ((N - (k0 + M + 1) : ℕ) : ℝ) ≤ (N : ℝ) := by
    have : (k0 - M) + (N - (k0 + M + 1)) ≤ N := by omega
    exact_mod_cast this
  constructor
  · rw [hsplit, Finset.mul_sum]
    have hmid : ∑ i ∈ range (2 * M + 1),
        a N k0 * Real.exp (-C / N) * Real.exp (-(12 + ε) * ((i : ℝ) - M) ^ 2 / N)
        ≤ ∑ i ∈ range (2 * M + 1), a N (k0 - M + i) := by
      refine Finset.sum_le_sum (fun i hi => ?_)
      have hi' := Finset.mem_range.mp hi
      have h := (hw (k0 - M + i) (by omega) (by omega)).1
      rw [hcast] at h
      have e : -(12 + ε) * ((i : ℝ) - M) ^ 2 / N - C / N
          = -C / N + -(12 + ε) * ((i : ℝ) - M) ^ 2 / N := by ring
      rw [e, Real.exp_add, ← mul_assoc] at h
      exact h
    linarith
  · rw [hsplit]
    have hmid : ∑ i ∈ range (2 * M + 1), a N (k0 - M + i)
        ≤ a N k0 * Real.exp (C / N) *
          ∑ i ∈ range (2 * M + 1), Real.exp (-(12 - ε) * ((i : ℝ) - M) ^ 2 / N) := by
      rw [Finset.mul_sum]
      refine Finset.sum_le_sum (fun i hi => ?_)
      have hi' := Finset.mem_range.mp hi
      have h := (hw (k0 - M + i) (by omega) (by omega)).2
      rw [hcast] at h
      have e : -(12 - ε) * ((i : ℝ) - M) ^ 2 / N + C / N
          = C / N + -(12 - ε) * ((i : ℝ) - M) ^ 2 / N := by ring
      rw [e, Real.exp_add, ← mul_assoc] at h
      exact h
    have hsum : ((k0 - M : ℕ) : ℝ) * (B * q ^ N) + ((N - (k0 + M + 1) : ℕ) : ℝ) * (B * q ^ N)
        ≤ (N : ℝ) * (B * q ^ N) := by
      rw [← add_mul]
      exact mul_le_mul_of_nonneg_right hcard hBq
    linarith

/-- for one `ε`: eventual upper and lower bounds for `T N` with explicit limits -/
lemma T_squeeze (B ψ : ℝ) (hB : 0 ≤ B)
    (hpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < a n k)
    (htail : ∀ δ : ℝ, 0 < δ → δ < 1 / 10 → ∃ q : ℝ, 0 < q ∧ q < 1 ∧ ∃ N0 : ℕ, ∀ N ≥ N0,
      ∀ k : ℕ, k < N → δ * N ≤ |(k : ℝ) - xstar * N| → a N k ≤ B * q ^ N)
    (hwin : ∀ ε : ℝ, 0 < ε → ∃ δ > 0, ∃ C : ℝ, ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ,
      |(k : ℝ) - ⌊xstar * N⌋₊| ≤ δ * N →
      (⌊xstar * N⌋₊ ≤ k →
        Real.exp (-(12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N)
          ≤ ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i ∧
        ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i
          ≤ Real.exp (-(12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)) ∧
      (k ≤ ⌊xstar * N⌋₊ →
        Real.exp ((12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N)
          ≤ ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i ∧
        ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i
          ≤ Real.exp ((12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)))
    (hgauss : ∀ c δ : ℝ, 0 < c → 0 < δ → Tendsto (fun N : ℕ => (1 / Real.sqrt N) *
        ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-c * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N))
      atTop (𝓝 (Real.sqrt (Real.pi / c))))
    (hcentre : Tendsto (fun N : ℕ => a N ⌊xstar * N⌋₊ * Real.sqrt N) atTop (𝓝 ψ))
    (ε : ℝ) (hε : 0 < ε) (hε1 : ε < 1) :
    (∃ U : ℕ → ℝ, Tendsto U atTop (𝓝 (ψ * Real.sqrt (Real.pi / (12 - ε)))) ∧
      ∀ᶠ N in atTop, T N ≤ U N) ∧
    (∃ V : ℕ → ℝ, Tendsto V atTop (𝓝 (ψ * Real.sqrt (Real.pi / (12 + ε)))) ∧
      ∀ᶠ N in atTop, V N ≤ T N) := by
  obtain ⟨δ1, hδ1, C, N1, hW⟩ := hwin ε hε
  set δ : ℝ := min δ1 (1 / 20) with hδdef
  have hδ : 0 < δ := lt_min hδ1 (by norm_num)
  have hδle1 : δ ≤ δ1 := min_le_left _ _
  have hδle2 : δ ≤ 1 / 20 := min_le_right _ _
  obtain ⟨q, hq0, hq1, N2, hT⟩ := htail (δ / 2) (by positivity) (by linarith)
  have hxg := xstar_gt_086
  have hxl := xstar_lt_087
  -- the bounds for all large `N`
  have hev : ∀ᶠ N : ℕ in atTop,
      a N ⌊xstar * N⌋₊ * Real.exp (-C / N) *
          ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-(12 + ε) * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N)
        ≤ T N ∧
      T N ≤ a N ⌊xstar * N⌋₊ * Real.exp (C / N) *
          ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-(12 - ε) * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N)
        + N * (B * q ^ N) := by
    filter_upwards [eventually_ge_atTop N1, eventually_ge_atTop N2, eventually_ge_atTop 13,
      eventually_ge_atTop ⌈2 / δ⌉₊] with N hN1 hN2 hN13 hNδ
    have hNpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
    have hN13' : (13 : ℝ) ≤ N := by exact_mod_cast hN13
    have hNδ' : 2 / δ ≤ N := (Nat.le_ceil _).trans (by exact_mod_cast hNδ)
    have hδN : 2 ≤ δ * N := by
      rw [div_le_iff₀ hδ] at hNδ'; linarith
    set M := ⌊δ * N⌋₊ with hMdef
    set k0 := ⌊xstar * N⌋₊ with hk0def
    have hMle : (M : ℝ) ≤ δ * N := Nat.floor_le (by positivity)
    have hMlt : δ * N < (M : ℝ) + 1 := Nat.lt_floor_add_one _
    have hk0le : (k0 : ℝ) ≤ xstar * N := Nat.floor_le (by nlinarith)
    have hk0lt : xstar * N < (k0 : ℝ) + 1 := Nat.lt_floor_add_one _
    have hM : M ≤ k0 := Nat.floor_le_floor (by nlinarith)
    have hk0N : k0 + M + 1 ≤ N := by
      have : (k0 : ℝ) + M + 1 ≤ N := by nlinarith
      exact_mod_cast this
    have hposN : ∀ k, k < N → 0 < a N k := fun k hk => hpos N k (by omega) hk
    refine T_bounds hM hk0N hposN ?_ ?_ (by positivity)
    · intro k hk1 hk2
      have habs : |(k : ℝ) - k0| ≤ δ1 * N := by
        have h1 : (k0 : ℝ) ≤ k + M := by exact_mod_cast (show k0 ≤ k + M by omega)
        have h2 : (k : ℝ) ≤ k0 + M := by exact_mod_cast hk2
        have h3 : δ * N ≤ δ1 * N := mul_le_mul_of_nonneg_right hδle1 hNpos.le
        rw [abs_le]; constructor <;> linarith
      obtain ⟨h1, h2⟩ := hW N hN1 k habs
      exact window_pointwise (hposN k (by omega)) (hposN k0 (by omega)) h1 h2
    · intro k hk hcase
      refine hT N hN2 k hk ?_
      rcases hcase with h | h
      · have h1 : (k : ℝ) + M + 1 ≤ k0 := by exact_mod_cast (show k + M + 1 ≤ k0 by omega)
        have h2 : δ / 2 * N ≤ -((k : ℝ) - xstar * N) := by nlinarith
        exact h2.trans (neg_le_abs _)
      · have h1 : (k0 : ℝ) + M + 1 ≤ k := by exact_mod_cast (show k0 + M + 1 ≤ k by omega)
        have h2 : δ / 2 * N ≤ (k : ℝ) - xstar * N := by nlinarith
        exact h2.trans (le_abs_self _)
  have hexp : ∀ c : ℝ, Tendsto (fun N : ℕ => Real.exp (c / N)) atTop (𝓝 1) := by
    intro c
    have h := (Real.continuous_exp.tendsto 0).comp (tendsto_const_div_atTop_nhds_zero_nat c)
    rw [Real.exp_zero] at h
    exact h
  have hNq : Tendsto (fun N : ℕ => (N : ℝ) * (B * q ^ N)) atTop (𝓝 0) := by
    have := (tendsto_self_mul_const_pow_of_lt_one hq0.le hq1).const_mul B
    rw [mul_zero] at this
    refine this.congr (fun N => ?_)
    ring
  constructor
  · refine ⟨fun N : ℕ => (a N ⌊xstar * N⌋₊ * Real.sqrt N) * Real.exp (C / N) *
        ((1 / Real.sqrt N) * ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1),
          Real.exp (-(12 - ε) * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N)) + N * (B * q ^ N), ?_, ?_⟩
    · have h := ((hcentre.mul (hexp C)).mul (hgauss (12 - ε) δ (by linarith) hδ)).add hNq
      simpa using h
    · filter_upwards [hev, eventually_ge_atTop 1] with N hN hN1
      have hs : Real.sqrt N ≠ 0 := by
        have : (0 : ℝ) < N := by exact_mod_cast hN1
        positivity
      refine hN.2.trans (le_of_eq ?_)
      field_simp
  · refine ⟨fun N : ℕ => (a N ⌊xstar * N⌋₊ * Real.sqrt N) * Real.exp (-C / N) *
        ((1 / Real.sqrt N) * ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1),
          Real.exp (-(12 + ε) * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N)), ?_, ?_⟩
    · have h := (hcentre.mul (hexp (-C))).mul (hgauss (12 + ε) δ (by linarith) hδ)
      simpa using h
    · filter_upwards [hev, eventually_ge_atTop 1] with N hN hN1
      have hs : Real.sqrt N ≠ 0 := by
        have : (0 : ℝ) < N := by exact_mod_cast hN1
        positivity
      refine le_trans (le_of_eq ?_) hN.1
      field_simp

/-- **Assembly of the Laplace argument**: from the tail bound, the window profile, the Gaussian
sum and the value at the centre, `T N → ψ √(π/12)`. -/
theorem tendsto_T_of (B ψ : ℝ) (hB : 0 ≤ B)
    (hpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < a n k)
    (htail : ∀ δ : ℝ, 0 < δ → δ < 1 / 10 → ∃ q : ℝ, 0 < q ∧ q < 1 ∧ ∃ N0 : ℕ, ∀ N ≥ N0,
      ∀ k : ℕ, k < N → δ * N ≤ |(k : ℝ) - xstar * N| → a N k ≤ B * q ^ N)
    (hwin : ∀ ε : ℝ, 0 < ε → ∃ δ > 0, ∃ C : ℝ, ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ,
      |(k : ℝ) - ⌊xstar * N⌋₊| ≤ δ * N →
      (⌊xstar * N⌋₊ ≤ k →
        Real.exp (-(12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N)
          ≤ ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i ∧
        ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i
          ≤ Real.exp (-(12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)) ∧
      (k ≤ ⌊xstar * N⌋₊ →
        Real.exp ((12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N)
          ≤ ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i ∧
        ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i
          ≤ Real.exp ((12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)))
    (hgauss : ∀ c δ : ℝ, 0 < c → 0 < δ → Tendsto (fun N : ℕ => (1 / Real.sqrt N) *
        ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-c * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N))
      atTop (𝓝 (Real.sqrt (Real.pi / c))))
    (hcentre : Tendsto (fun N : ℕ => a N ⌊xstar * N⌋₊ * Real.sqrt N) atTop (𝓝 ψ)) :
    Tendsto T atTop (𝓝 (ψ * Real.sqrt (Real.pi / 12))) := by
  have key := T_squeeze B ψ hB hpos htail hwin hgauss hcentre
  have hcont : ∀ s : ℝ, ContinuousAt (fun e : ℝ => ψ * Real.sqrt (Real.pi / (12 + s * e))) 0 := by
    intro s
    refine continuousAt_const.mul (Real.continuous_sqrt.continuousAt.comp ?_)
    refine continuousAt_const.div (continuousAt_const.add (continuousAt_const.mul continuousAt_id)) ?_
    norm_num
  rw [tendsto_order]
  constructor
  · intro b hb
    have h0 : b < (fun e : ℝ => ψ * Real.sqrt (Real.pi / (12 + 1 * e))) 0 := by simpa using hb
    obtain ⟨r, hr, hball⟩ := Metric.eventually_nhds_iff.mp ((hcont 1).eventually (lt_mem_nhds h0))
    have hε0 : 0 < min (r / 2) (1 / 2) := lt_min (by positivity) (by norm_num)
    have hεr : dist (min (r / 2) (1 / 2)) 0 < r := by
      rw [Real.dist_eq, sub_zero, abs_of_pos hε0]
      exact (min_le_left _ _).trans_lt (by linarith)
    have hb' := hball hεr
    simp only [one_mul] at hb'
    obtain ⟨V, hV, hVT⟩ := (key _ hε0 ((min_le_right _ _).trans_lt (by norm_num))).2
    filter_upwards [hV.eventually (lt_mem_nhds hb'), hVT] with N h1 h2
    exact h1.trans_le h2
  · intro b hb
    have h0 : (fun e : ℝ => ψ * Real.sqrt (Real.pi / (12 + (-1) * e))) 0 < b := by simpa using hb
    obtain ⟨r, hr, hball⟩ := Metric.eventually_nhds_iff.mp ((hcont (-1)).eventually (gt_mem_nhds h0))
    have hε0 : 0 < min (r / 2) (1 / 2) := lt_min (by positivity) (by norm_num)
    have hεr : dist (min (r / 2) (1 / 2)) 0 < r := by
      rw [Real.dist_eq, sub_zero, abs_of_pos hε0]
      exact (min_le_left _ _).trans_lt (by linarith)
    have hb' := hball hεr
    have e : (12 : ℝ) + -1 * min (r / 2) (1 / 2) = 12 - min (r / 2) (1 / 2) := by ring
    rw [e] at hb'
    obtain ⟨U, hU, hTU⟩ := (key _ hε0 ((min_le_right _ _).trans_lt (by norm_num))).1
    filter_upwards [hU.eventually (gt_mem_nhds hb'), hTU] with N h1 h2
    exact h2.trans_lt h1

end PiSeries.SunPi2
