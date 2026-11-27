import Mathlib

/-!
# Riemann-sum form of the Gaussian integral (package WP-G)

`(1/√N) ∑_{|j| ≤ δN} exp(-c j²/N) → √(π/c)`.
-/

open Finset Filter Topology

namespace PiSeries.SunPi2

/-- Symmetric reindexing of a sum of an even function over `|j| ≤ M`. -/
private lemma sum_symm_range (h : ℝ → ℝ) (heven : ∀ x, h (-x) = h x) (M : ℕ) :
    ∑ i ∈ range (2 * M + 1), h ((i : ℝ) - M) =
      h 0 + 2 * ∑ j ∈ range M, h ((j : ℝ) + 1) := by
  have h1 : 2 * M + 1 = M + (M + 1) := by ring
  rw [h1, Finset.sum_range_add, Finset.sum_range_succ']
  have h2 : ∑ i ∈ range M, h ((i : ℝ) - M) = ∑ j ∈ range M, h ((j : ℝ) + 1) := by
    rw [← Finset.sum_range_reflect]
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [Finset.mem_range] at hj
    rw [← heven]
    congr 1
    have : ((M - 1 - j : ℕ) : ℝ) = (M : ℝ) - 1 - j := by
      rw [Nat.cast_sub (by omega), Nat.cast_sub (by omega)]; simp
    rw [this]; ring
  have h3 : ∑ k ∈ range M, h (((M + (k + 1) : ℕ) : ℝ) - M) =
      ∑ j ∈ range M, h ((j : ℝ) + 1) := by
    refine Finset.sum_congr rfl fun j _ => ?_
    congr 1; push_cast; ring
  rw [h2, h3]
  simp
  ring

/-- Explicit error bound for the Gaussian Riemann sum. -/
private lemma gaussian_sum_bound (c δ : ℝ) (hc : 0 < c) (hδ : 0 < δ) (N : ℕ) (hN : 1 ≤ N) :
    |(1 / Real.sqrt N) *
        ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-c * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N)
      - Real.sqrt (Real.pi / c)| ≤ (1 + 2 / (c * δ)) / Real.sqrt N := by
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hsqN : 0 < Real.sqrt N := Real.sqrt_pos.mpr hNpos
  set M : ℕ := ⌊δ * N⌋₊ with hM
  set b : ℝ := c / N with hb
  have hbpos : 0 < b := div_pos hc hNpos
  set g : ℝ → ℝ := fun s => Real.exp (-b * s ^ 2) with hg
  have hcδ : 0 < c * δ := mul_pos hc hδ
  -- rewrite the sum in terms of `g`
  have hsum : ∑ i ∈ range (2 * M + 1), Real.exp (-c * ((i : ℝ) - M) ^ 2 / N) =
      ∑ i ∈ range (2 * M + 1), g ((i : ℝ) - M) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [hg, hb]
    congr 1; ring
  have heven : ∀ x, g (-x) = g x := fun x => by simp [hg]
  have hg0 : g 0 = 1 := by simp [hg]
  have hsymm := sum_symm_range g heven M
  rw [hg0] at hsymm
  set A : ℝ := ∑ j ∈ range M, g ((j : ℝ) + 1) with hA
  -- basic facts about `g`
  have hanti : AntitoneOn g (Set.Ici 0) := by
    intro x hx y _ hxy
    simp only [hg]
    apply Real.exp_le_exp.mpr
    have hx0 : (0 : ℝ) ≤ x := hx
    have : x ^ 2 ≤ y ^ 2 := pow_le_pow_left₀ hx0 hxy 2
    nlinarith
  have hint : MeasureTheory.Integrable g := integrable_exp_neg_mul_sq hbpos
  have hI : ∫ x in Set.Ioi (0 : ℝ), g x = Real.sqrt (Real.pi / b) / 2 :=
    integral_gaussian_Ioi b
  have hsqrt : Real.sqrt (Real.pi / b) = Real.sqrt (Real.pi / c) * Real.sqrt N := by
    rw [← Real.sqrt_mul (div_pos Real.pi_pos hc).le]
    congr 1
    rw [hb]; field_simp
  -- upper bound
  have hup : A ≤ Real.sqrt (Real.pi / b) / 2 := by
    have := AntitoneOn.sum_range_le_integral (f := g) (N := M)
      (hanti.mono Set.Icc_subset_Ici_self) hint.integrableOn
      (fun t _ => (Real.exp_pos _).le)
    rw [hI] at this
    refine le_trans (le_of_eq ?_) this
    refine Finset.sum_congr rfl fun j _ => ?_
    push_cast; rfl
  -- lower bound
  have hlow : Real.sqrt (Real.pi / b) / 2 - 1 / (c * δ) ≤ 1 + A := by
    have h1 := AntitoneOn.integral_le_sum (x₀ := 0) (a := M + 1) (f := g)
      (hanti.mono Set.Icc_subset_Ici_self)
    have h2 : ∑ i ∈ range (M + 1), g (0 + (i : ℝ)) = 1 + A := by
      rw [Finset.sum_range_succ']
      simp only [zero_add, Nat.cast_zero, hg0, hA]
      rw [add_comm]
      congr 1
      refine Finset.sum_congr rfl fun j _ => ?_
      push_cast; rfl
    rw [h2, zero_add] at h1
    have hM1 : (0 : ℝ) ≤ ((M + 1 : ℕ) : ℝ) := Nat.cast_nonneg _
    have h3 := intervalIntegral.integral_Ioi_sub_Ioi (f := g) (μ := MeasureTheory.volume)
      hint.integrableOn hM1
    rw [hI] at h3
    -- tail bound
    have htail : ∫ x in Set.Ioi (((M + 1 : ℕ) : ℝ)), g x ≤ 1 / (c * δ) := by
      have hlt : δ * N < ((M + 1 : ℕ) : ℝ) := by
        push_cast; exact Nat.lt_floor_add_one _
      have hneg : -(c * δ) < 0 := by linarith
      have hle : ∫ x in Set.Ioi (((M + 1 : ℕ) : ℝ)), g x ≤
          ∫ x in Set.Ioi (((M + 1 : ℕ) : ℝ)), Real.exp (-(c * δ) * x) := by
        refine MeasureTheory.setIntegral_mono_on hint.integrableOn
          (integrableOn_exp_mul_Ioi hneg _) measurableSet_Ioi fun x hx => ?_
        have hx' : δ * N < x := lt_trans hlt hx
        have hx0 : 0 < x := lt_trans (mul_pos hδ hNpos) hx'
        simp only [hg]
        apply Real.exp_le_exp.mpr
        have : c * δ * x ≤ b * x ^ 2 := by
          rw [hb, div_mul_eq_mul_div, le_div_iff₀ hNpos]
          nlinarith [mul_pos hc hx0, mul_nonneg (mul_pos hc hx0).le (sub_nonneg.mpr hx'.le)]
        linarith
      rw [integral_exp_mul_Ioi hneg] at hle
      refine le_trans hle ?_
      rw [neg_div_neg_eq]
      apply div_le_div_of_nonneg_right _ hcδ.le
      apply Real.exp_le_one_iff.mpr
      have : 0 ≤ c * δ * ((M + 1 : ℕ) : ℝ) := mul_nonneg hcδ.le hM1
      linarith
    linarith
  -- combine
  rw [hsum, hsymm, abs_le]
  have hrw : 1 / Real.sqrt N * (1 + 2 * A) - Real.sqrt (Real.pi / c) =
      (1 + 2 * A - Real.sqrt (Real.pi / b)) / Real.sqrt N := by
    rw [hsqrt]; field_simp
  rw [hrw, ← neg_div]
  constructor
  · apply div_le_div_of_nonneg_right _ hsqN.le
    have : 2 / (c * δ) = 2 * (1 / (c * δ)) := by ring
    linarith
  · apply div_le_div_of_nonneg_right _ hsqN.le
    have : 0 < 2 / (c * δ) := by positivity
    linarith

theorem tendsto_gaussian_sum (c δ : ℝ) (hc : 0 < c) (hδ : 0 < δ) :
    Tendsto (fun N : ℕ => (1 / Real.sqrt N) *
        ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-c * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N))
      atTop (𝓝 (Real.sqrt (Real.pi / c))) := by
  rw [tendsto_iff_norm_sub_tendsto_zero]
  have hlim : Tendsto (fun N : ℕ => (1 + 2 / (c * δ)) / Real.sqrt N) atTop (𝓝 0) :=
    tendsto_const_nhds.div_atTop (Real.tendsto_sqrt_atTop.comp tendsto_natCast_atTop_atTop)
  refine squeeze_zero' (Eventually.of_forall fun N => norm_nonneg _) ?_ hlim
  filter_upwards [eventually_ge_atTop 1] with N hN
  rw [Real.norm_eq_abs]
  exact gaussian_sum_bound c δ hc hδ N hN

end PiSeries.SunPi2
