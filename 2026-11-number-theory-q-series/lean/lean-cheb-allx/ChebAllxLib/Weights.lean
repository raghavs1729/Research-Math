/-
ChebAllxLib/Weights.lean — (1) the integral form of the identity for Chebyshev's sum,
`−C(x) = ∫_2^∞ Θ₁(t) F''(t) dt` with `F(t) = e^{−t/x}/log t` (Lemma 1(ii) of
`monotonicity/chebyshev1853/CHEBYSHEV1853.md`; no hypothesis); (2) Fujii's weights
`e^{−(p/x)^α}`, `0 < α ≤ 1`, and the weights `p^{−σ}`, `σ > 1`, as instances of the theorem for
convex weights.
-/
import ChebAllxLib.Main
import ChebAllxLib.Integral

open Filter Topology MeasureTheory Set

noncomputable section

namespace ChebAllx

/-! ### The integral identity for Chebyshev's sum -/

/-- `c(n) = [n prime] (−χ₋₄(n)) log n`. -/
def primeCoeff (n : ℕ) : ℝ := if n.Prime then Coupling.w4 n * Real.log n else 0

lemma theta1_eq_sum_primeCoeff (t : ℝ) :
    Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t
      = ∑ n ∈ Finset.Icc 0 ⌊t⌋₊, primeCoeff n * (t - n) := by
  unfold Coupling.Theta1 primeCoeff
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl fun n _ => ?_
  split_ifs
  · rfl
  · rw [zero_mul]

/-- **The integral identity (no hypothesis).**  For every `x > 0`, with
`Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p)` and `F(t) = e^{−t/x}/log t`:
`Σ_p χ₋₄(p) e^{−p/x} = −∫_2^∞ Θ₁(t) F''(t) dt`, the integrand being integrable. -/
theorem chebyshevSum_eq_neg_integral {x : ℝ} (hx : 0 < x) :
    IntegrableOn (fun t => Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t * chebF'' x t)
        (Ioi 2) ∧
      chebyshevSum x = -∫ t in Ioi (2 : ℝ),
        Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t * chebF'' x t := by
  have hc : ∀ n : ℕ, (n : ℝ) < 2 → primeCoeff n = 0 := by
    intro n hn
    have hn2 : n < 2 := by exact_mod_cast hn
    unfold primeCoeff
    rw [if_neg]
    intro hp
    exact absurd hp.two_le (by omega)
  have hs : Summable fun n : ℕ => |primeCoeff n| * chebF x n := by
    have h : Summable fun n : ℕ => Real.exp (n * (-1 / x)) :=
      Real.summable_exp_nat_mul_iff.2 (by
        have : 0 < 1 / x := one_div_pos.2 hx
        rw [neg_div]; linarith)
    refine Summable.of_nonneg_of_le (fun n => ?_) (fun n => ?_) h
    · exact mul_nonneg (abs_nonneg _)
        (div_nonneg (Real.exp_pos _).le (Real.log_natCast_nonneg n))
    · unfold primeCoeff chebF
      have he : Real.exp (-(n : ℝ) / x) = Real.exp (n * (-1 / x)) := by congr 1; ring
      split_ifs with hp
      · have hlog : 0 < Real.log n := Real.log_pos (by exact_mod_cast hp.one_lt)
        rw [abs_mul, abs_of_pos hlog, ← he]
        have h1 := abs_w4_le_one n
        calc |Coupling.w4 n| * Real.log n * (Real.exp (-(n : ℝ) / x) / Real.log n)
            = |Coupling.w4 n| * Real.exp (-(n : ℝ) / x) := by field_simp
          _ ≤ 1 * Real.exp (-(n : ℝ) / x) := mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le
          _ = Real.exp (-(n : ℝ) / x) := one_mul _
      · rw [abs_zero, zero_mul]
        exact (Real.exp_pos _).le
  have key := tsum_mul_eq_integral (a := 2) (by norm_num) hc
    (fun t ht => hasDerivAt_chebF hx.ne' (by linarith))
    (fun t ht => hasDerivAt_chebF' hx.ne' (by linarith))
    (continuousOn_chebF'' hx.ne') (fun t ht => (chebF''_pos hx (by linarith)).le)
    (tendsto_chebF hx) hs
  simp only [← theta1_eq_sum_primeCoeff] at key
  refine ⟨key.1, ?_⟩
  rw [← key.2]
  unfold chebyshevSum
  rw [ChebLib.tsum_primes_eq (fun n => ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)
    * Real.exp (-(n : ℝ) / x)), ← tsum_neg]
  refine tsum_congr fun n => ?_
  unfold primeCoeff chebF Coupling.w4
  split_ifs with hp
  · have hlog : Real.log n ≠ 0 := (Real.log_pos (by exact_mod_cast hp.one_lt)).ne'
    field_simp
  · rw [zero_mul, neg_zero]

/-! ### Fujii's weights and powers -/

lemma convexOn_fujii {x α : ℝ} (hx : 0 < x) (hα0 : 0 ≤ α) (hα1 : α ≤ 1) :
    ConvexOn ℝ (Ici 3) (fun t : ℝ => Real.exp (-(t / x) ^ α)) := by
  have hg : ConcaveOn ℝ (Ici (3 : ℝ)) (fun t : ℝ => (t / x) ^ α) := by
    have h1 : ConcaveOn ℝ (Ici (0 : ℝ)) (fun t : ℝ => t ^ α) := Real.concaveOn_rpow hα0 hα1
    have h2 : ConcaveOn ℝ (Ici (0 : ℝ)) (fun t : ℝ => (x ^ α)⁻¹ • t ^ α) :=
      h1.smul (inv_nonneg.2 (Real.rpow_nonneg hx.le α))
    refine (h2.subset (Ici_subset_Ici.2 (by norm_num)) (convex_Ici 3)).congr ?_
    intro t ht
    have ht0 : 0 ≤ t := by linarith [mem_Ici.1 ht]
    simp only [smul_eq_mul]
    rw [Real.div_rpow ht0 hx.le, div_eq_inv_mul]
  exact convexOn_exp_neg_of_concaveOn hg

lemma summable_fujii {x α : ℝ} (hx : 0 < x) (hα0 : 0 < α) :
    Summable fun n : ℕ => Real.exp (-((n : ℝ) / x) ^ α) := by
  have h1 := tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero (2 / α) 1 one_pos
  have hy : Tendsto (fun n : ℕ => ((n : ℝ) / x) ^ α) atTop atTop :=
    (tendsto_rpow_atTop hα0).comp ((tendsto_natCast_atTop_atTop).atTop_div_const hx)
  have h2 := h1.comp hy
  refine Summable.of_norm_bounded_eventually_nat
    (g := fun n : ℕ => x ^ 2 * (1 / (n : ℝ) ^ 2))
    ((Real.summable_one_div_nat_pow.2 one_lt_two).mul_left _) ?_
  filter_upwards [(tendsto_order.1 h2).2 1 zero_lt_one, eventually_gt_atTop 0] with n hn hn0
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn0
  have hnx : 0 ≤ (n : ℝ) / x := by positivity
  have e : (((n : ℝ) / x) ^ α) ^ (2 / α) = ((n : ℝ) / x) ^ (2 : ℕ) := by
    rw [← Real.rpow_mul hnx, mul_div_cancel₀ _ hα0.ne', Real.rpow_two]
  simp only [Function.comp_apply, e, neg_mul, one_mul] at hn
  rw [Real.norm_of_nonneg (Real.exp_pos _).le, mul_one_div, le_div_iff₀ (by positivity)]
  rw [div_pow, div_mul_eq_mul_div, div_lt_one (by positivity)] at hn
  linarith

lemma convexOn_rpow_neg {σ : ℝ} (hσ : 0 ≤ σ) : ConvexOn ℝ (Ici 3) (fun t : ℝ => t ^ (-σ)) := by
  have hg : ConcaveOn ℝ (Ici (3 : ℝ)) (fun t : ℝ => σ • Real.log t) :=
    ((strictConcaveOn_log_Ioi.concaveOn).subset
      (fun t ht => lt_of_lt_of_le (by norm_num) (mem_Ici.1 ht)) (convex_Ici 3)).smul hσ
  refine (convexOn_exp_neg_of_concaveOn hg).congr fun t ht => ?_
  have ht0 : 0 < t := lt_of_lt_of_le (by norm_num) (mem_Ici.1 ht)
  simp only [smul_eq_mul]
  rw [Real.rpow_def_of_pos ht0]
  congr 1
  ring

lemma summable_chi4_fujii {x α : ℝ} (hx : 0 < x) (hα0 : 0 < α) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(((p : ℕ) : ℝ) / x) ^ α)) := by
  refine Summable.of_norm_bounded ((summable_fujii hx hα0).subtype {n : ℕ | n.Prime}) fun p => ?_
  rw [norm_mul, Real.norm_of_nonneg (Real.exp_pos _).le, Real.norm_eq_abs]
  have h := abs_w4_le_one (p : ℕ)
  unfold Coupling.w4 at h
  rw [abs_neg] at h
  exact (mul_le_mul_of_nonneg_right h (Real.exp_pos _).le).trans (by rw [one_mul]; rfl)

/-- **Fujii's weights.**  Under GRH for `L(s, χ₋₄)`: `Σ_p χ₋₄(p) e^{−(p/x)^α} < 0` for every `x > 0`
and every `0 < α ≤ 1`. -/
theorem fujii_sum_neg (hGRH : GRHFor chi4) {x α : ℝ} (hx : 0 < x) (hα0 : 0 < α) (hα1 : α ≤ 1) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(((p : ℕ) : ℝ) / x) ^ α) < 0 :=
  tsum_mul_neg_of_convexOn hGRH (convexOn_fujii hx hα0.le hα1)
    ((summable_fujii hx hα0).subtype {n : ℕ | n.Prime}) (Real.exp_pos _)

/-- Under GRH for `L(s, χ₋₄)`: `Σ_p χ₋₄(p) p^{−σ} < 0` for every `σ > 1`. -/
theorem prime_dirichlet_sum_neg (hGRH : GRHFor chi4) {σ : ℝ} (hσ : 1 < σ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * ((p : ℕ) : ℝ) ^ (-σ) < 0 :=
  tsum_mul_neg_of_convexOn hGRH (convexOn_rpow_neg (by linarith))
    ((Real.summable_nat_rpow.2 (by linarith : -σ < -1)).subtype {n : ℕ | n.Prime})
    (Real.rpow_pos_of_pos (by norm_num) _)

/-! ### From the sum without logarithms to the sums with logarithms

`1/log p = ∫_0^∞ p^{−u} du`, hence `C(x) = ∫_0^∞ (Σ_p χ₋₄(p) log p · p^{−u} e^{−p/x}) du`, and each
integrand is a log-weighted sum with the convex weight `t^{−u} e^{−t/x}`. -/

lemma integral_rpow_neg_Ioi {n : ℕ} (hn : 2 ≤ n) :
    IntegrableOn (fun u : ℝ => (n : ℝ) ^ (-u)) (Ioi 0) ∧
      ∫ u in Ioi (0 : ℝ), (n : ℝ) ^ (-u) = 1 / Real.log n := by
  have hn0 : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
  have hL : 0 < Real.log n := Real.log_pos (by exact_mod_cast (by omega : 1 < n))
  have e : (fun u : ℝ => (n : ℝ) ^ (-u)) = fun u => Real.exp (-Real.log n * u) := by
    funext u
    rw [Real.rpow_def_of_pos hn0]
    congr 1
    ring
  rw [e]
  refine ⟨integrableOn_exp_mul_Ioi (neg_lt_zero.2 hL) 0, ?_⟩
  rw [integral_exp_mul_Ioi (neg_lt_zero.2 hL) 0, mul_zero, Real.exp_zero]
  field_simp

/-- **From `C` to the log-weighted sums (no hypothesis).**  For every `x > 0`:
`Σ_p χ₋₄(p) e^{−p/x} = ∫_0^∞ (Σ_p χ₋₄(p) log p · p^{−u} e^{−p/x}) du`. -/
theorem chebyshevSum_eq_integral_log_sums {x : ℝ} (hx : 0 < x) :
    chebyshevSum x = ∫ u in Ioi (0 : ℝ), ∑' p : Nat.Primes,
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
        * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  set G : Nat.Primes → ℝ → ℝ := fun p u =>
    ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x)) with hG
  have hlog : ∀ p : Nat.Primes, 0 < Real.log (p : ℕ) := fun p =>
    Real.log_pos (by exact_mod_cast p.2.one_lt)
  have hG' : ∀ (p : Nat.Primes) (u : ℝ), G p u
      = (((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
          * Real.exp (-((p : ℕ) : ℝ) / x)) * ((p : ℕ) : ℝ) ^ (-u) := by
    intro p u
    simp only [hG]
    ring
  have hint : ∀ p, Integrable (G p) (volume.restrict (Ioi 0)) := by
    intro p
    have h := (integral_rpow_neg_Ioi p.2.two_le).1.const_mul
      (((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * Real.exp (-((p : ℕ) : ℝ) / x))
    exact h.congr (Filter.Eventually.of_forall fun u => (hG' p u).symm)
  have hval : ∀ p : Nat.Primes, ∫ u in Ioi (0 : ℝ), G p u
      = ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x) := by
    intro p
    simp_rw [hG']
    rw [integral_const_mul, (integral_rpow_neg_Ioi p.2.two_le).2]
    have := (hlog p).ne'
    field_simp
  have hnorm : ∀ p : Nat.Primes, ∫ u in Ioi (0 : ℝ), ‖G p u‖
      = |((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)| * Real.exp (-((p : ℕ) : ℝ) / x) := by
    intro p
    have hp0 : (0 : ℝ) < (p : ℕ) := by exact_mod_cast p.2.pos
    have h1 : ∀ u : ℝ, ‖G p u‖ = (|((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)| * Real.log (p : ℕ)
        * Real.exp (-((p : ℕ) : ℝ) / x)) * ((p : ℕ) : ℝ) ^ (-u) := by
      intro u
      rw [hG', Real.norm_eq_abs, abs_mul, abs_mul, abs_mul, abs_of_pos (hlog p),
        abs_of_pos (Real.exp_pos _), abs_of_pos (Real.rpow_pos_of_pos hp0 _)]
    simp_rw [h1]
    rw [integral_const_mul, (integral_rpow_neg_Ioi p.2.two_le).2]
    have := (hlog p).ne'
    field_simp
  have hsum : Summable fun p : Nat.Primes => ∫ u in Ioi (0 : ℝ), ‖G p u‖ := by
    simp_rw [hnorm]
    refine Summable.of_nonneg_of_le (fun p => by positivity) (fun p => ?_)
      (summable_primes_exp_neg_div hx)
    have h := abs_w4_le_one (p : ℕ)
    unfold Coupling.w4 at h
    rw [abs_neg] at h
    exact (mul_le_mul_of_nonneg_right h (Real.exp_pos _).le).trans (by rw [one_mul])
  have h := (hasSum_integral_of_summable_integral_norm hint hsum).tsum_eq
  simp_rw [hval] at h
  unfold chebyshevSum
  rw [h]

lemma convexOn_rpow_neg_mul_exp {x u : ℝ} (hx : 0 < x) (hu : 0 ≤ u) :
    ConvexOn ℝ (Ici 3) (fun t : ℝ => t ^ (-u) * Real.exp (-t / x)) := by
  have hlogc : ConcaveOn ℝ (Ici (3 : ℝ)) Real.log :=
    (strictConcaveOn_log_Ioi.concaveOn).subset
      (fun t ht => lt_of_lt_of_le (by norm_num) (mem_Ici.1 ht)) (convex_Ici 3)
  have hg : ConcaveOn ℝ (Ici (3 : ℝ)) (fun t : ℝ => u • Real.log t + (1 / x) • t) :=
    (hlogc.smul hu).add ((concaveOn_id (convex_Ici 3)).smul (one_div_pos.2 hx).le)
  refine (convexOn_exp_neg_of_concaveOn hg).congr fun t ht => ?_
  have ht0 : 0 < t := lt_of_lt_of_le (by norm_num) (mem_Ici.1 ht)
  simp only [smul_eq_mul]
  rw [Real.rpow_def_of_pos ht0, ← Real.exp_add]
  congr 1
  ring

/-- Under GRH for `L(s, χ₋₄)`: for every `x > 0` and `u ≥ 0`,
`Σ_p χ₋₄(p) log p · p^{−u} e^{−p/x} < 0`. -/
theorem log_sum_rpow_exp_neg (hGRH : GRHFor chi4) {x u : ℝ} (hx : 0 < x) (hu : 0 ≤ u) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x)) < 0 := by
  refine tsum_log_mul_neg_of_convexOn hGRH (convexOn_rpow_neg_mul_exp hx hu) ?_
    (mul_pos (Real.rpow_pos_of_pos (by norm_num) _) (Real.exp_pos _))
  refine Summable.of_nonneg_of_le (fun p => ?_) (fun p => ?_)
    (summable_primes_log_mul_exp_neg_div hx)
  · have hp0 : (0 : ℝ) < (p : ℕ) := by exact_mod_cast p.2.pos
    exact mul_nonneg (Real.log_natCast_nonneg _)
      (mul_nonneg (Real.rpow_pos_of_pos hp0 _).le (Real.exp_pos _).le)
  · have hp1 : (1 : ℝ) ≤ (p : ℕ) := by exact_mod_cast p.2.one_lt.le
    have h1 : ((p : ℕ) : ℝ) ^ (-u) ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos hp1 (by linarith)
    exact mul_le_mul_of_nonneg_left
      ((mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le).trans (by rw [one_mul]))
      (Real.log_natCast_nonneg _)

end ChebAllx
