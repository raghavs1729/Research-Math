/-
# WP-R (N2-PLAN.md): the ratio `σ(N,k)` in the window `|k − x*N| ≤ δN`.

`σ(N,k) = Λ(k/N)·(1+O(1/N))`, `Λ(x) = (1−x)²(2x+1)⁴/((1+x)²(2x−1)⁴)`, `Λ(x*) = 1`,
`(log Λ)'(x*) = −24`.  Hence `log σ(N,k) = −24(k−x*N)/N + o(|k−x*N|/N) + O(1/N)` and, multiplying up,
`∏ σ(N,i)` over `i` between `k0 = ⌊x*N⌋` and `k` is `exp(∓(12±ε)(k−k0)²/N + O(1/N))`.
-/
import PiSeries.SunPi2.Defs

open Finset Filter Topology

namespace PiSeries.SunPi2

/-! ### elementary logarithm bounds -/

lemma log_sub_log_nonneg {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) : 0 ≤ Real.log b - Real.log a :=
  sub_nonneg.mpr (Real.log_le_log ha hab)

lemma log_sub_log_le {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) :
    Real.log b - Real.log a ≤ (b - a) / a := by
  have hb : 0 < b := lt_of_lt_of_le ha hab
  rw [← Real.log_div hb.ne' ha.ne']
  have h := Real.log_le_sub_one_of_pos (div_pos hb ha)
  calc Real.log (b / a) ≤ b / a - 1 := h
    _ = (b - a) / a := by field_simp

lemma log_sub_log_le' {a b n c : ℝ} (hn : 0 < n) (ha : 0 < a) (hab : a ≤ b)
    (h : (b - a) * n ≤ c * a) : Real.log b - Real.log a ≤ c / n :=
  (log_sub_log_le ha hab).trans (by rw [div_le_div_iff₀ ha hn]; exact h)

/-! ### the saddle point -/

lemma xstar_sq_R : xstar ^ 2 = 3 / 4 := by
  unfold xstar
  rw [div_pow, Real.sq_sqrt (by norm_num)]
  norm_num

lemma xstar_pos_R : 0 < xstar := by unfold xstar; positivity

lemma xstar_gt_R : (0.86 : ℝ) < xstar := by
  by_contra h
  rw [not_lt] at h
  nlinarith [xstar_sq_R, xstar_pos_R]

lemma xstar_lt_R : xstar < (0.87 : ℝ) := by
  by_contra h
  rw [not_lt] at h
  nlinarith [xstar_sq_R, xstar_pos_R]

/-- `g = log Λ` -/
noncomputable def g (x : ℝ) : ℝ :=
  2 * Real.log (1 - x) + 4 * Real.log (2 * x + 1) - 2 * Real.log (1 + x) - 4 * Real.log (2 * x - 1)

lemma g_xstar : g xstar = 0 := by
  have h1 := xstar_gt_R
  have h2 := xstar_lt_R
  have key : (1 - xstar) * (2 * xstar + 1) ^ 2 = (1 + xstar) * (2 * xstar - 1) ^ 2 := by
    linear_combination (-8 * xstar) * xstar_sq_R
  have h := congrArg Real.log key
  rw [Real.log_mul (by linarith) (by nlinarith), Real.log_mul (by linarith) (by nlinarith),
    Real.log_pow, Real.log_pow] at h
  push_cast at h
  unfold g
  linarith

lemma hasDerivAt_g : HasDerivAt g (-24) xstar := by
  have hx1 := xstar_gt_R
  have hx2 := xstar_lt_R
  have ne1 : (1 - xstar) ≠ 0 := by linarith
  have ne2 : (2 * xstar + 1) ≠ 0 := by linarith
  have ne3 : (1 + xstar) ≠ 0 := by linarith
  have ne4 : (2 * xstar - 1) ≠ 0 := by linarith
  have h1 : HasDerivAt (fun x : ℝ => 1 - x) (-1) xstar := by
    simpa using (hasDerivAt_id xstar).const_sub 1
  have h2 : HasDerivAt (fun x : ℝ => 2 * x + 1) 2 xstar := by
    simpa using ((hasDerivAt_id xstar).const_mul 2).add_const 1
  have h3 : HasDerivAt (fun x : ℝ => 1 + x) 1 xstar := by
    simpa using (hasDerivAt_id xstar).const_add 1
  have h4 : HasDerivAt (fun x : ℝ => 2 * x - 1) 2 xstar := by
    simpa using ((hasDerivAt_id xstar).const_mul 2).sub_const 1
  have h := ((((h1.log ne1).const_mul 2).add ((h2.log ne2).const_mul 4)).sub
    ((h3.log ne3).const_mul 2)).sub ((h4.log ne4).const_mul 4)
  have e : (2 * (-1 / (1 - xstar)) + 4 * (2 / (2 * xstar + 1)) - 2 * (1 / (1 + xstar))
      - 4 * (2 / (2 * xstar - 1)) : ℝ) = -24 := by
    have e1 : (1 - xstar) * (1 + xstar) = 1 / 4 := by linear_combination (-1 : ℝ) * xstar_sq_R
    have e2 : (2 * xstar + 1) * (2 * xstar - 1) = 2 := by linear_combination (4 : ℝ) * xstar_sq_R
    have a1 : (2 * (-1 / (1 - xstar)) - 2 * (1 / (1 + xstar)) : ℝ) = -16 := by
      have : (2 * (-1 / (1 - xstar)) - 2 * (1 / (1 + xstar)) : ℝ)
          = -4 / ((1 - xstar) * (1 + xstar)) := by
        field_simp
        ring
      rw [this, e1]; norm_num
    have a2 : (4 * (2 / (2 * xstar + 1)) - 4 * (2 / (2 * xstar - 1)) : ℝ) = -8 := by
      have : (4 * (2 / (2 * xstar + 1)) - 4 * (2 / (2 * xstar - 1)) : ℝ)
          = -16 / ((2 * xstar + 1) * (2 * xstar - 1)) := by
        field_simp
        ring
      rw [this, e2]; norm_num
    linarith
  rw [e] at h
  exact h

lemma g_window (ε : ℝ) (hε : 0 < ε) :
    ∃ δ > 0, ∀ x : ℝ, |x - xstar| < δ → |g x + 24 * (x - xstar)| ≤ ε * |x - xstar| := by
  have h := hasDerivAt_iff_isLittleO.mp hasDerivAt_g
  have h2 := h.def hε
  rw [Metric.eventually_nhds_iff] at h2
  obtain ⟨δ, hδ, H⟩ := h2
  refine ⟨δ, hδ, fun x hx => ?_⟩
  have h3 := H (y := x) (by rwa [Real.dist_eq])
  simp only [g_xstar, Real.norm_eq_abs, smul_eq_mul, sub_zero] at h3
  have e : g x + 24 * (x - xstar) = g x - (x - xstar) * -24 := by ring
  rw [e]
  exact h3

/-- `g(k/n)` in terms of the linear forms occurring in `σ` -/
lemma g_div (n k : ℝ) (hn : 0 < n) (h1 : 0 < n - k) (h2 : 0 < 8 * k + 4 * n) (h3 : 0 < n + k)
    (h4 : 0 < 8 * k - 4 * n) :
    g (k / n) = 2 * Real.log (n - k) + 4 * Real.log (8 * k + 4 * n) - 2 * Real.log (n + k)
      - 4 * Real.log (8 * k - 4 * n) := by
  have e1 : 1 - k / n = (n - k) / n := by field_simp
  have e2 : 2 * (k / n) + 1 = (8 * k + 4 * n) / (4 * n) := by field_simp; ring
  have e3 : 1 + k / n = (n + k) / n := by field_simp
  have e4 : 2 * (k / n) - 1 = (8 * k - 4 * n) / (4 * n) := by field_simp; ring
  have h4n : (4 * n) ≠ 0 := by positivity
  unfold g
  rw [e1, e2, e3, e4, Real.log_div h1.ne' hn.ne', Real.log_div h2.ne' h4n,
    Real.log_div h3.ne' hn.ne', Real.log_div h4.ne' h4n]
  ring

/-! ### `log σ` versus `g(k/N)` -/

lemma log_quot_form (a b c1 c2 c3 c4 d e f1 f2 f3 f4 : ℝ) (ha : 0 < a) (hb : 0 < b)
    (hc1 : 0 < c1) (hc2 : 0 < c2) (hc3 : 0 < c3) (hc4 : 0 < c4) (hd : 0 < d) (he : 0 < e)
    (hf1 : 0 < f1) (hf2 : 0 < f2) (hf3 : 0 < f3) (hf4 : 0 < f4) :
    Real.log (a * b ^ 2 * (c1 * c2 * c3 * c4) / (d * e ^ 2 * (f1 * f2 * f3 * f4)))
      = Real.log a + 2 * Real.log b + Real.log c1 + Real.log c2 + Real.log c3 + Real.log c4
        - Real.log d - 2 * Real.log e - Real.log f1 - Real.log f2 - Real.log f3 - Real.log f4 := by
  rw [Real.log_div (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow]
  push_cast
  ring

/-- the real-variable core of `log_sigma_window` -/
lemma log_sigma_real (n k : ℝ) (hn : 20 ≤ n) (hk1 : 0.83 * n ≤ k) (hk2 : k ≤ 0.9 * n) :
    0 < (2 * k + 3) * (n - k - 1) ^ 2
          * ((8 * k + 4 * n + 1) * (8 * k + 4 * n + 3) * (8 * k + 4 * n + 5) * (8 * k + 4 * n + 7)) /
        ((2 * k + 1) * (n + k + 1) ^ 2
          * ((8 * k - 4 * n + 9) * (8 * k - 4 * n + 11) * (8 * k - 4 * n + 13)
            * (8 * k - 4 * n + 15))) ∧
    |Real.log ((2 * k + 3) * (n - k - 1) ^ 2
          * ((8 * k + 4 * n + 1) * (8 * k + 4 * n + 3) * (8 * k + 4 * n + 5) * (8 * k + 4 * n + 7)) /
        ((2 * k + 1) * (n + k + 1) ^ 2
          * ((8 * k - 4 * n + 9) * (8 * k - 4 * n + 11) * (8 * k - 4 * n + 13)
            * (8 * k - 4 * n + 15)))) - g (k / n)| ≤ 66 / n := by
  have hn0 : 0 < n := by linarith
  have p1 : 0 < 2 * k + 1 := by linarith
  have p2 : 0 < n - k - 1 := by linarith
  have p3 : 0 < n + k := by linarith
  have p4 : 0 < 8 * k + 4 * n := by linarith
  have p5 : 0 < 8 * k - 4 * n := by linarith
  have q1 : 0 < 2 * k + 3 := by linarith
  have q2 : 0 < n - k := by linarith
  have q3 : 0 < n + k + 1 := by linarith
  have c1 : 0 < 8 * k + 4 * n + 1 := by linarith
  have c2 : 0 < 8 * k + 4 * n + 3 := by linarith
  have c3 : 0 < 8 * k + 4 * n + 5 := by linarith
  have c4 : 0 < 8 * k + 4 * n + 7 := by linarith
  have f1 : 0 < 8 * k - 4 * n + 9 := by linarith
  have f2 : 0 < 8 * k - 4 * n + 11 := by linarith
  have f3 : 0 < 8 * k - 4 * n + 13 := by linarith
  have f4 : 0 < 8 * k - 4 * n + 15 := by linarith
  refine ⟨by positivity, ?_⟩
  rw [log_quot_form _ _ _ _ _ _ _ _ _ _ _ _ q1 p2 c1 c2 c3 c4 p1 q3 f1 f2 f3 f4,
    g_div n k hn0 q2 p4 p3 p5]
  -- the nine brackets
  have A0 := log_sub_log_nonneg p1 (show 2 * k + 1 ≤ 2 * k + 3 by linarith)
  have A1 := log_sub_log_le' (c := 2) hn0 p1 (show 2 * k + 1 ≤ 2 * k + 3 by linarith) (by linarith)
  have B0 := log_sub_log_nonneg p2 (show n - k - 1 ≤ n - k by linarith)
  have B1 := log_sub_log_le' (c := 20) hn0 p2 (show n - k - 1 ≤ n - k by linarith) (by linarith)
  have D0 := log_sub_log_nonneg p3 (show n + k ≤ n + k + 1 by linarith)
  have D1 := log_sub_log_le' (c := 1) hn0 p3 (show n + k ≤ n + k + 1 by linarith) (by linarith)
  have E10 := log_sub_log_nonneg p4 (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 1 by linarith)
  have E11 := log_sub_log_le' (c := 1 / 4) hn0 p4
    (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 1 by linarith) (by linarith)
  have E20 := log_sub_log_nonneg p4 (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 3 by linarith)
  have E21 := log_sub_log_le' (c := 3 / 4) hn0 p4
    (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 3 by linarith) (by linarith)
  have E30 := log_sub_log_nonneg p4 (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 5 by linarith)
  have E31 := log_sub_log_le' (c := 5 / 4) hn0 p4
    (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 5 by linarith) (by linarith)
  have E40 := log_sub_log_nonneg p4 (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 7 by linarith)
  have E41 := log_sub_log_le' (c := 7 / 4) hn0 p4
    (show 8 * k + 4 * n ≤ 8 * k + 4 * n + 7 by linarith) (by linarith)
  have F10 := log_sub_log_nonneg p5 (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 9 by linarith)
  have F11 := log_sub_log_le' (c := 9 / 2) hn0 p5
    (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 9 by linarith) (by linarith)
  have F20 := log_sub_log_nonneg p5 (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 11 by linarith)
  have F21 := log_sub_log_le' (c := 11 / 2) hn0 p5
    (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 11 by linarith) (by linarith)
  have F30 := log_sub_log_nonneg p5 (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 13 by linarith)
  have F31 := log_sub_log_le' (c := 13 / 2) hn0 p5
    (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 13 by linarith) (by linarith)
  have F40 := log_sub_log_nonneg p5 (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 15 by linarith)
  have F41 := log_sub_log_le' (c := 15 / 2) hn0 p5
    (show 8 * k - 4 * n ≤ 8 * k - 4 * n + 15 by linarith) (by linarith)
  -- express everything through `n⁻¹`
  simp only [div_eq_mul_inv] at A1 B1 D1 E11 E21 E31 E41 F11 F21 F31 F41 ⊢
  generalize n⁻¹ = r at A1 B1 D1 E11 E21 E31 E41 F11 F21 F31 F41 ⊢
  rw [abs_le]
  constructor <;> linarith

theorem log_sigma_window (ε : ℝ) (hε : 0 < ε) :
    ∃ δ > 0, ∃ C : ℝ, ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ, |(k : ℝ) - xstar * N| ≤ δ * N →
      0 < σ N k ∧
      |Real.log (σ N k) + 24 * ((k : ℝ) - xstar * N) / N| ≤ ε * |(k : ℝ) - xstar * N| / N + C / N := by
  obtain ⟨δ0, hδ0, H⟩ := g_window ε hε
  refine ⟨min (δ0 / 2) 0.03, by positivity, 66, 20, fun N hN k hk => ?_⟩
  have hNr : (20 : ℝ) ≤ N := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < N := by linarith
  have hδ1 : min (δ0 / 2) 0.03 ≤ 0.03 := min_le_right _ _
  have hδ2 : min (δ0 / 2) 0.03 ≤ δ0 / 2 := min_le_left _ _
  have hk' := abs_le.mp hk
  have m1 : min (δ0 / 2) 0.03 * (N : ℝ) ≤ 0.03 * N := mul_le_mul_of_nonneg_right hδ1 hNpos.le
  have m2 : 0.86 * (N : ℝ) ≤ xstar * N := mul_le_mul_of_nonneg_right xstar_gt_R.le hNpos.le
  have m3 : xstar * (N : ℝ) ≤ 0.87 * N := mul_le_mul_of_nonneg_right xstar_lt_R.le hNpos.le
  have hk1 : 0.83 * (N : ℝ) ≤ k := by linarith [hk'.1]
  have hk2 : (k : ℝ) ≤ 0.9 * N := by linarith [hk'.2]
  obtain ⟨hpos, hlog⟩ := log_sigma_real N k hNr hk1 hk2
  have e1 : (k : ℝ) / N - xstar = ((k : ℝ) - xstar * N) / N := by field_simp
  have hx : |(k : ℝ) / N - xstar| < δ0 := by
    rw [e1, abs_div, abs_of_pos hNpos, div_lt_iff₀ hNpos]
    have : min (δ0 / 2) 0.03 * (N : ℝ) ≤ δ0 / 2 * N := mul_le_mul_of_nonneg_right hδ2 hNpos.le
    nlinarith
  have hg := H _ hx
  rw [e1, abs_div, abs_of_pos hNpos] at hg
  refine ⟨hpos, ?_⟩
  have e2 : Real.log (σ N k) + 24 * ((k : ℝ) - xstar * N) / N
      = (Real.log (σ N k) - g (k / N)) + (g (k / N) + 24 * (((k : ℝ) - xstar * N) / N)) := by ring
  rw [e2]
  have hlog' : |Real.log (σ N k) - g (k / N)| ≤ 66 / N := hlog
  calc |(Real.log (σ N k) - g (k / N)) + (g (k / N) + 24 * (((k : ℝ) - xstar * N) / N))|
      ≤ |Real.log (σ N k) - g (k / N)| + |g (k / N) + 24 * (((k : ℝ) - xstar * N) / N)| :=
        abs_add_le _ _
    _ ≤ 66 / N + ε * (|(k : ℝ) - xstar * N| / N) := add_le_add hlog' hg
    _ = ε * |(k : ℝ) - xstar * N| / N + 66 / N := by ring

/-! ### summing up -/

lemma sum_window_aux (u : ℕ → ℝ) (c n W e C : ℝ) (hn : 0 < n) (he : 0 ≤ e) (a b : ℕ) (hab : a ≤ b)
    (h : ∀ i ∈ Ico a b, |u i + 24 * ((i : ℝ) - c) / n| ≤ e * |(i : ℝ) - c| / n + C / n)
    (hW : ∀ i ∈ Ico a b, |(i : ℝ) - c| ≤ W) :
    |∑ i ∈ Ico a b, u i + 12 * (((b : ℝ) - c) ^ 2 - ((a : ℝ) - c) ^ 2) / n|
      ≤ ((b : ℝ) - a) * (e * W + C + 12) / n := by
  induction b, hab using Nat.le_induction with
  | base => simp
  | succ b hab ih =>
    have ih' := ih (fun i hi => h i (by rw [mem_Ico] at hi ⊢; omega))
      (fun i hi => hW i (by rw [mem_Ico] at hi ⊢; omega))
    have hb := h b (by rw [mem_Ico]; omega)
    have hWb := hW b (by rw [mem_Ico]; omega)
    rw [Finset.sum_Ico_succ_top hab]
    have e1 : ∑ i ∈ Ico a b, u i + u b + 12 * ((((b + 1 : ℕ) : ℝ) - c) ^ 2 - ((a : ℝ) - c) ^ 2) / n
        = (∑ i ∈ Ico a b, u i + 12 * (((b : ℝ) - c) ^ 2 - ((a : ℝ) - c) ^ 2) / n)
          + (u b + 24 * ((b : ℝ) - c) / n) + 12 / n := by
      push_cast; ring
    rw [e1]
    have hW' : e * |(b : ℝ) - c| / n ≤ e * W / n := by gcongr
    have e2 : (((b + 1 : ℕ) : ℝ) - a) * (e * W + C + 12) / n
        = ((b : ℝ) - a) * (e * W + C + 12) / n + (e * W / n + C / n + 12 / n) := by
      push_cast; ring
    rw [e2]
    have h12 : 0 ≤ 12 / n := by positivity
    calc |(∑ i ∈ Ico a b, u i + 12 * (((b : ℝ) - c) ^ 2 - ((a : ℝ) - c) ^ 2) / n)
          + (u b + 24 * ((b : ℝ) - c) / n) + 12 / n|
        ≤ |∑ i ∈ Ico a b, u i + 12 * (((b : ℝ) - c) ^ 2 - ((a : ℝ) - c) ^ 2) / n|
          + |u b + 24 * ((b : ℝ) - c) / n| + |12 / n| := abs_add_three _ _ _
      _ ≤ ((b : ℝ) - a) * (e * W + C + 12) / n + (e * W / n + C / n + 12 / n) := by
          rw [abs_of_nonneg h12]; linarith

lemma window_core (N : ℕ) (hN : (0 : ℝ) < N) (e C₁ δ₁ : ℝ) (he : 0 ≤ e)
    (H : ∀ i : ℕ, |(i : ℝ) - xstar * N| ≤ δ₁ * N → 0 < σ N i ∧
      |Real.log (σ N i) + 24 * ((i : ℝ) - xstar * N) / N| ≤ e * |(i : ℝ) - xstar * N| / N + C₁ / N)
    (a b : ℕ) (hab : a ≤ b) (Y : ℝ) (hY : ∀ i ∈ Ico a b, |(i : ℝ) - xstar * N| ≤ Y)
    (hYδ : Y ≤ δ₁ * N) :
    0 < ∏ i ∈ Ico a b, σ N i ∧
    |Real.log (∏ i ∈ Ico a b, σ N i)
        + 12 * (((b : ℝ) - xstar * N) ^ 2 - ((a : ℝ) - xstar * N) ^ 2) / N|
      ≤ ((b : ℝ) - a) * (e * Y + C₁ + 12) / N := by
  have hpos : ∀ i ∈ Ico a b, 0 < σ N i := fun i hi => (H i ((hY i hi).trans hYδ)).1
  refine ⟨Finset.prod_pos hpos, ?_⟩
  rw [Real.log_prod (fun i hi => (hpos i hi).ne')]
  exact sum_window_aux (fun i => Real.log (σ N i)) (xstar * N) N Y e C₁ hN he a b hab
    (fun i hi => (H i ((hY i hi).trans hYδ)).2) hY

lemma absorb (ε D y n : ℝ) (hε : 0 < ε) (hn : 0 < n) :
    (ε / 2 * y ^ 2 + D * y) / n ≤ ε * y ^ 2 / n + D ^ 2 / (2 * ε) / n := by
  have h : D * y ≤ ε / 2 * y ^ 2 + D ^ 2 / (2 * ε) := by
    rw [← sub_nonneg]
    have e : ε / 2 * y ^ 2 + D ^ 2 / (2 * ε) - D * y = (ε * y - D) ^ 2 / (2 * ε) := by
      field_simp; ring
    rw [e]; positivity
  rw [← add_div]
  apply div_le_div_of_nonneg_right _ hn.le
  linarith

theorem prod_sigma_window (ε : ℝ) (hε : 0 < ε) :
    ∃ δ > 0, ∃ C : ℝ, ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ, |(k : ℝ) - ⌊xstar * N⌋₊| ≤ δ * N →
      (⌊xstar * N⌋₊ ≤ k →
        Real.exp (-(12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N) ≤ ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i ∧
        ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i ≤ Real.exp (-(12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)) ∧
      (k ≤ ⌊xstar * N⌋₊ →
        Real.exp ((12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N) ≤ ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i ∧
        ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i ≤ Real.exp ((12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)) := by
  obtain ⟨δ₁, hδ₁, C₀, N₁, H₀⟩ := log_sigma_window (ε / 2) (by positivity)
  set C₁ : ℝ := max C₀ 0 with hC₁
  have hC₁0 : 0 ≤ C₁ := le_max_right _ _
  set D : ℝ := ε / 2 + C₁ + 36 with hD
  refine ⟨δ₁ / 2, by positivity, D ^ 2 / (2 * ε), max N₁ (⌈2 / δ₁⌉₊ + 1), fun N hN k hk => ?_⟩
  have hN1 : N₁ ≤ N := le_trans (le_max_left _ _) hN
  have hN2 : ⌈2 / δ₁⌉₊ + 1 ≤ N := le_trans (le_max_right _ _) hN
  have hNpos : (0 : ℝ) < N := by
    have : 1 ≤ N := by omega
    exact_mod_cast this
  have hN3 : 2 ≤ δ₁ * N := by
    have h1 : 2 / δ₁ ≤ (⌈2 / δ₁⌉₊ : ℝ) := Nat.le_ceil _
    have h2 : (⌈2 / δ₁⌉₊ : ℝ) ≤ N := by
      have : ⌈2 / δ₁⌉₊ ≤ N := by omega
      exact_mod_cast this
    have h3 : 2 / δ₁ ≤ (N : ℝ) := h1.trans h2
    rw [div_le_iff₀ hδ₁] at h3
    linarith
  -- the pointwise estimate with the nonnegative constant `C₁`
  have H : ∀ i : ℕ, |(i : ℝ) - xstar * N| ≤ δ₁ * N → 0 < σ N i ∧
      |Real.log (σ N i) + 24 * ((i : ℝ) - xstar * N) / N|
        ≤ ε / 2 * |(i : ℝ) - xstar * N| / N + C₁ / N := by
    intro i hi
    obtain ⟨h1, h2⟩ := H₀ N hN1 i hi
    refine ⟨h1, h2.trans ?_⟩
    have : C₀ / N ≤ C₁ / N := div_le_div_of_nonneg_right (le_max_left _ _) hNpos.le
    linarith
  have hθ0 : ((⌊xstar * N⌋₊ : ℕ) : ℝ) ≤ xstar * N :=
    Nat.floor_le (mul_nonneg xstar_pos_R.le hNpos.le)
  have hθ1 : xstar * N < ((⌊xstar * N⌋₊ : ℕ) : ℝ) + 1 := Nat.lt_floor_add_one _
  set k0 : ℕ := ⌊xstar * N⌋₊ with hk0
  set θ : ℝ := xstar * N - k0 with hθ
  have hθa : 0 ≤ θ := by rw [hθ]; linarith
  have hθb : θ < 1 := by rw [hθ]; linarith
  have hk' := abs_le.mp hk
  constructor
  · intro hle
    have hleR : (k0 : ℝ) ≤ k := by exact_mod_cast hle
    set y : ℝ := (k : ℝ) - k0 with hy
    have hy0 : 0 ≤ y := by rw [hy]; linarith
    have hY : ∀ i ∈ Ico k0 k, |(i : ℝ) - xstar * N| ≤ y + 1 := by
      intro i hi
      rw [mem_Ico] at hi
      have h1 : (k0 : ℝ) ≤ i := by exact_mod_cast hi.1
      have h2 : (i : ℝ) ≤ k := by exact_mod_cast hi.2.le
      rw [abs_le]
      constructor <;> linarith
    obtain ⟨hP, hL⟩ := window_core N hNpos (ε / 2) C₁ δ₁ (by positivity) H k0 k hle (y + 1) hY
      (by linarith [hk'.2])
    set L : ℝ := Real.log (∏ i ∈ Ico k0 k, σ N i) with hLdef
    have e1 : L + 12 * y ^ 2 / N
        = (L + 12 * (((k : ℝ) - xstar * N) ^ 2 - ((k0 : ℝ) - xstar * N) ^ 2) / N)
          + 24 * y * θ / N := by
      rw [hy, hθ]; ring
    have hyθ : 24 * y * θ / N ≤ 24 * y / N := by
      apply div_le_div_of_nonneg_right _ hNpos.le
      nlinarith
    have hyθ0 : 0 ≤ 24 * y * θ / N := by positivity
    have hmain : |L + 12 * y ^ 2 / N| ≤ ε * y ^ 2 / N + D ^ 2 / (2 * ε) / N := by
      rw [e1]
      refine (abs_add_le _ _).trans ?_
      rw [abs_of_nonneg hyθ0]
      refine le_trans ?_ (absorb ε D y N hε hNpos)
      have e2 : (ε / 2 * y ^ 2 + D * y) / N
          = y * (ε / 2 * (y + 1) + C₁ + 12) / N + 24 * y / N := by
        rw [hD]; ring
      rw [e2]
      exact add_le_add hL hyθ
    have hmain' := abs_le.mp hmain
    rw [← Real.exp_log hP]
    constructor
    · apply Real.exp_le_exp.mpr
      have e3 : -(12 + ε) * y ^ 2 / N - D ^ 2 / (2 * ε) / N
          = -(12 * y ^ 2 / N) - (ε * y ^ 2 / N + D ^ 2 / (2 * ε) / N) := by ring
      rw [e3]
      linarith [hmain'.1]
    · apply Real.exp_le_exp.mpr
      have e3 : -(12 - ε) * y ^ 2 / N + D ^ 2 / (2 * ε) / N
          = -(12 * y ^ 2 / N) + (ε * y ^ 2 / N + D ^ 2 / (2 * ε) / N) := by ring
      rw [e3]
      linarith [hmain'.2]
  · intro hle
    have hleR : (k : ℝ) ≤ k0 := by exact_mod_cast hle
    set y : ℝ := (k0 : ℝ) - k with hy
    have hy0 : 0 ≤ y := by rw [hy]; linarith
    have hsq : ((k : ℝ) - k0) ^ 2 = y ^ 2 := by rw [hy]; ring
    have hY : ∀ i ∈ Ico k k0, |(i : ℝ) - xstar * N| ≤ y + 1 := by
      intro i hi
      rw [mem_Ico] at hi
      have h1 : (k : ℝ) ≤ i := by exact_mod_cast hi.1
      have h2 : (i : ℝ) ≤ k0 := by exact_mod_cast hi.2.le
      rw [abs_le]
      constructor <;> linarith
    obtain ⟨hP, hL⟩ := window_core N hNpos (ε / 2) C₁ δ₁ (by positivity) H k k0 hle (y + 1) hY
      (by linarith [hk'.1])
    set L : ℝ := Real.log (∏ i ∈ Ico k k0, σ N i) with hLdef
    have e1 : L - 12 * y ^ 2 / N
        = (L + 12 * (((k0 : ℝ) - xstar * N) ^ 2 - ((k : ℝ) - xstar * N) ^ 2) / N)
          + 24 * y * θ / N := by
      rw [hy, hθ]; ring
    have hyθ : 24 * y * θ / N ≤ 24 * y / N := by
      apply div_le_div_of_nonneg_right _ hNpos.le
      nlinarith
    have hyθ0 : 0 ≤ 24 * y * θ / N := by positivity
    have hmain : |L - 12 * y ^ 2 / N| ≤ ε * y ^ 2 / N + D ^ 2 / (2 * ε) / N := by
      rw [e1]
      refine (abs_add_le _ _).trans ?_
      rw [abs_of_nonneg hyθ0]
      refine le_trans ?_ (absorb ε D y N hε hNpos)
      have e2 : (ε / 2 * y ^ 2 + D * y) / N
          = (y * (ε / 2 * (y + 1) + C₁ + 12) / N) + 24 * y / N := by
        rw [hD]; ring
      rw [e2]
      exact add_le_add hL hyθ
    have hmain' := abs_le.mp hmain
    rw [← Real.exp_log hP, hsq]
    constructor
    · apply Real.exp_le_exp.mpr
      have e3 : (12 - ε) * y ^ 2 / N - D ^ 2 / (2 * ε) / N
          = 12 * y ^ 2 / N - (ε * y ^ 2 / N + D ^ 2 / (2 * ε) / N) := by ring
      rw [e3]
      linarith [hmain'.1]
    · apply Real.exp_le_exp.mpr
      have e3 : (12 + ε) * y ^ 2 / N + D ^ 2 / (2 * ε) / N
          = 12 * y ^ 2 / N + (ε * y ^ 2 / N + D ^ 2 / (2 * ε) / N) := by ring
      rw [e3]
      linarith [hmain'.2]

end PiSeries.SunPi2
