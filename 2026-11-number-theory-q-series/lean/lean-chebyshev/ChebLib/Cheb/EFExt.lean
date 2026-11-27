/-
ChebLib/Cheb/EFExt.lean — the explicit formula for test functions without compact support.

Zeta23's `EF_lit_chi_L` is the explicit formula of a primitive Dirichlet L-function for test
functions `k ∈ C²_c(ℝ)`.  Under the Riemann hypothesis for `L(s,χ)` (all ordinates `γ_ρ` real) we
extend it, by truncation `k_j = η_j k` and dominated convergence on each of its three sides, to
every `k` that is twice differentiable with continuous `k''`, with `k, k', k''` integrable, and
whose prime side converges absolutely.
-/
import ChebLib.Race.S0

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory Set Filter Topology
open Zeta23 Zeta23.EF Zeta23.ThmE

noncomputable section

set_option backward.isDefEq.respectTransparency false

namespace ChebLib

/-! ### The cut-off functions -/

/-- A fixed smooth bump: `1` on `[−1,1]`, `0` outside `(−2,2)`, values in `[0,1]`. -/
def bump0 : ContDiffBump (0 : ℝ) := ⟨1, 2, zero_lt_one, one_lt_two⟩

/-- `η_j(u) = η(u/(j+1))`. -/
def cutoff (j : ℕ) (u : ℝ) : ℝ := bump0 (u / ((j : ℝ) + 1))

theorem bump0_contDiff : ContDiff ℝ 2 (bump0 : ℝ → ℝ) := bump0.contDiff

theorem cutoff_nonneg (j : ℕ) (u : ℝ) : 0 ≤ cutoff j u := bump0.nonneg

theorem cutoff_le_one (j : ℕ) (u : ℝ) : cutoff j u ≤ 1 := bump0.le_one

theorem cutoff_eq_one {j : ℕ} {u : ℝ} (h : |u| ≤ (j : ℝ) + 1) : cutoff j u = 1 := by
  apply bump0.one_of_mem_closedBall
  have hj : (0 : ℝ) < (j : ℝ) + 1 := by positivity
  rw [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs, abs_div, abs_of_pos hj,
    div_le_iff₀ hj]
  simpa [bump0] using h

theorem tendsto_cutoff (u : ℝ) : Tendsto (fun j : ℕ => cutoff j u) atTop (𝓝 1) := by
  apply tendsto_const_nhds.congr'
  filter_upwards [eventually_ge_atTop ⌈|u|⌉₊] with j hj
  refine (cutoff_eq_one ?_).symm
  have h1 : |u| ≤ (⌈|u|⌉₊ : ℝ) := Nat.le_ceil _
  have h2 : (⌈|u|⌉₊ : ℝ) ≤ (j : ℝ) := by exact_mod_cast hj
  linarith

theorem cutoff_contDiff (j : ℕ) : ContDiff ℝ 2 (cutoff j) :=
  bump0_contDiff.comp (contDiff_id.div_const _)

theorem cutoff_hasCompactSupport (j : ℕ) : HasCompactSupport (cutoff j) := by
  have hj : (0 : ℝ) < (j : ℝ) + 1 := by positivity
  apply HasCompactSupport.of_support_subset_isCompact (isCompact_closedBall (0 : ℝ) (2 * (j + 1)))
  intro u hu
  rw [Function.mem_support] at hu
  by_contra hmem
  apply hu
  apply bump0.zero_of_le_dist
  rw [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs, not_le] at hmem
  rw [dist_zero_right, Real.norm_eq_abs, abs_div, abs_of_pos hj, le_div_iff₀ hj]
  simpa [bump0] using hmem.le

/-- Bounds for the first two derivatives of the bump. -/
theorem exists_bump0_deriv_bounds : ∃ B₁ B₂ : ℝ, 0 ≤ B₁ ∧ 0 ≤ B₂ ∧
    (∀ v, |deriv (bump0 : ℝ → ℝ) v| ≤ B₁) ∧ (∀ v, |deriv (deriv (bump0 : ℝ → ℝ)) v| ≤ B₂) := by
  have hc1 : Continuous (deriv (bump0 : ℝ → ℝ)) := bump0_contDiff.continuous_deriv (by norm_num)
  have hd1 : ContDiff ℝ 1 (deriv (bump0 : ℝ → ℝ)) := by
    have := bump0_contDiff
    rw [show (2 : WithTop ℕ∞) = 1 + 1 by norm_num, contDiff_succ_iff_deriv] at this
    exact this.2.2
  have hc2 : Continuous (deriv (deriv (bump0 : ℝ → ℝ))) := hd1.continuous_deriv le_rfl
  have hs0 : HasCompactSupport (bump0 : ℝ → ℝ) := bump0.hasCompactSupport
  obtain ⟨B₁, hB₁⟩ := hc1.bounded_above_of_compact_support hs0.deriv
  obtain ⟨B₂, hB₂⟩ := hc2.bounded_above_of_compact_support hs0.deriv.deriv
  refine ⟨B₁, B₂, (norm_nonneg _).trans (hB₁ 0), (norm_nonneg _).trans (hB₂ 0), hB₁, hB₂⟩

theorem hasDerivAt_cutoff (j : ℕ) (u : ℝ) :
    HasDerivAt (cutoff j) (deriv (bump0 : ℝ → ℝ) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1)) u := by
  have h0 : HasDerivAt (bump0 : ℝ → ℝ) (deriv (bump0 : ℝ → ℝ) (u / ((j : ℝ) + 1)))
      (u / ((j : ℝ) + 1)) :=
    (bump0_contDiff.differentiable (by norm_num)).differentiableAt.hasDerivAt
  have h1 : HasDerivAt (fun u : ℝ => u / ((j : ℝ) + 1)) (1 / ((j : ℝ) + 1)) u :=
    (hasDerivAt_id u).div_const _
  have := h0.comp u h1
  have e : (bump0 : ℝ → ℝ) ∘ (fun u : ℝ => u / ((j : ℝ) + 1)) = cutoff j := rfl
  rw [e, mul_one_div] at this
  exact this

theorem hasDerivAt_cutoff' (j : ℕ) (u : ℝ) :
    HasDerivAt (fun u : ℝ => deriv (bump0 : ℝ → ℝ) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1))
      (deriv (deriv (bump0 : ℝ → ℝ)) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1) ^ 2) u := by
  have hd1 : ContDiff ℝ 1 (deriv (bump0 : ℝ → ℝ)) := by
    have := bump0_contDiff
    rw [show (2 : WithTop ℕ∞) = 1 + 1 by norm_num, contDiff_succ_iff_deriv] at this
    exact this.2.2
  have h0 : HasDerivAt (deriv (bump0 : ℝ → ℝ))
      (deriv (deriv (bump0 : ℝ → ℝ)) (u / ((j : ℝ) + 1))) (u / ((j : ℝ) + 1)) :=
    (hd1.differentiable one_ne_zero).differentiableAt.hasDerivAt
  have h1 : HasDerivAt (fun u : ℝ => u / ((j : ℝ) + 1)) (1 / ((j : ℝ) + 1)) u :=
    (hasDerivAt_id u).div_const _
  have := (h0.comp u h1).div_const ((j : ℝ) + 1)
  refine HasDerivAt.congr_deriv this ?_
  rw [mul_one_div, div_div, ← sq]

/-! ### The truncations of a test function -/

section Trunc

variable {k k1 k2 : ℝ → ℂ}

/-- The truncation `k_j = η_j k`. -/
def trunc (k : ℝ → ℂ) (j : ℕ) : ℝ → ℂ := fun u => (cutoff j u : ℂ) * k u

theorem contDiff_two_of_hasDerivAt (hk1 : ∀ u, HasDerivAt k (k1 u) u)
    (hk2 : ∀ u, HasDerivAt k1 (k2 u) u) (hk2c : Continuous k2) : ContDiff ℝ 2 k := by
  have e1 : deriv k = k1 := funext fun u => (hk1 u).deriv
  have e2 : deriv k1 = k2 := funext fun u => (hk2 u).deriv
  rw [show (2 : WithTop ℕ∞) = 1 + 1 by norm_num, contDiff_succ_iff_deriv]
  refine ⟨fun u => (hk1 u).differentiableAt, by simp, ?_⟩
  rw [e1, show (1 : WithTop ℕ∞) = 0 + 1 by norm_num, contDiff_succ_iff_deriv]
  refine ⟨fun u => (hk2 u).differentiableAt, by simp, ?_⟩
  rw [e2]
  exact contDiff_zero.mpr hk2c

theorem trunc_contDiff (hk : ContDiff ℝ 2 k) (j : ℕ) : ContDiff ℝ 2 (trunc k j) :=
  (Complex.ofRealCLM.contDiff.comp (cutoff_contDiff j)).mul hk

theorem trunc_hasCompactSupport (k : ℝ → ℂ) (j : ℕ) : HasCompactSupport (trunc k j) := by
  have h : HasCompactSupport (fun u : ℝ => (cutoff j u : ℂ)) :=
    (cutoff_hasCompactSupport j).comp_left Complex.ofReal_zero
  exact h.mul_right

theorem norm_trunc_le (k : ℝ → ℂ) (j : ℕ) (u : ℝ) : ‖trunc k j u‖ ≤ ‖k u‖ := by
  rw [trunc, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (cutoff_nonneg j u)]
  exact mul_le_of_le_one_left (norm_nonneg _) (cutoff_le_one j u)

theorem tendsto_trunc (k : ℝ → ℂ) (u : ℝ) : Tendsto (fun j : ℕ => trunc k j u) atTop (𝓝 (k u)) := by
  have h := ((Complex.continuous_ofReal.tendsto 1).comp (tendsto_cutoff u)).mul_const (k u)
  simpa [trunc] using h

theorem deriv_deriv_trunc (hk1 : ∀ u, HasDerivAt k (k1 u) u)
    (hk2 : ∀ u, HasDerivAt k1 (k2 u) u) (j : ℕ) (u : ℝ) :
    deriv (deriv (trunc k j)) u
      = ((deriv (deriv (bump0 : ℝ → ℝ)) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1) ^ 2 : ℝ) : ℂ) * k u
        + 2 * ((deriv (bump0 : ℝ → ℝ) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1) : ℝ) : ℂ) * k1 u
        + (cutoff j u : ℂ) * k2 u := by
  have hd : deriv (trunc k j) = fun u =>
      ((deriv (bump0 : ℝ → ℝ) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1) : ℝ) : ℂ) * k u
        + (cutoff j u : ℂ) * k1 u := by
    funext v
    exact ((hasDerivAt_cutoff j v).ofReal_comp.mul (hk1 v)).deriv
  rw [hd]
  have h := (((hasDerivAt_cutoff' j u).ofReal_comp.mul (hk1 u)).add
    ((hasDerivAt_cutoff j u).ofReal_comp.mul (hk2 u)))
  refine Eq.trans (show _ = _ from h.deriv) ?_
  ring

end Trunc

/-! ### The extension -/

variable {q : ℕ} [NeZero q] {χ : DirichletCharacter ℂ q}

/-- **The explicit formula for non-compactly supported test functions, under GRH.**
Let `χ` be a primitive character mod `q > 1` and assume the Riemann hypothesis for `L(s,χ)`.  Let
`k` be twice differentiable (`k' = k1`, `k'' = k2`, `k2` continuous) with `k, k', k''` integrable
and `Σ_n Λ(n) n^{-1/2} (|k(log n)| + |k(−log n)|) < ∞`.  Then the zero sum converges absolutely,
the transform obeys `‖k̂(τ)‖ (1+τ²) ≤ M` on the real line, and

    Σ_ρ m_ρ k̂(γ_ρ) = −Σ_n Λ(n) n^{-1/2} (χ(n) k(log n) + conj χ(n) k(−log n))
                      + (1/2π) ∫ k̂(r) [Re ψ(1/4 + κ/2 + ir/2) + log(q/π)] dr. -/
theorem EF_ext (hq : 1 < q) (hprim : χ.IsPrimitive) (hGRH : GRHFor χ)
    {k k1 k2 : ℝ → ℂ} (hk1 : ∀ u, HasDerivAt k (k1 u) u) (hk2 : ∀ u, HasDerivAt k1 (k2 u) u)
    (hk2c : Continuous k2) (hki : Integrable k) (hk1i : Integrable k1) (hk2i : Integrable k2)
    (hps : Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖k (Real.log n)‖ + ‖k (-Real.log n)‖))) :
    Summable (fun ρ : (LZeros (LSeam_of hq hprim)).carrier =>
        ((LZeros (LSeam_of hq hprim)).mult ρ : ℂ) * paperFT k (gammaOf ρ)) ∧
      ∑' ρ : (LZeros (LSeam_of hq hprim)).carrier,
          ((LZeros (LSeam_of hq hprim)).mult ρ : ℂ) * paperFT k (gammaOf ρ)
        = literatureRHSchi (parity χ) q (coeff χ) k := by
  classical
  set Z := LZeros (LSeam_of hq hprim) with hZ
  have hq1 : 1 ≤ q := hq.le
  have hkC2 : ContDiff ℝ 2 k := contDiff_two_of_hasDerivAt hk1 hk2 hk2c
  obtain ⟨B₁, B₂, hB₁0, hB₂0, hB₁, hB₂⟩ := exists_bump0_deriv_bounds
  set M : ℝ := (∫ u, ‖k u‖) + (B₂ * (∫ u, ‖k u‖) + 2 * B₁ * (∫ u, ‖k1 u‖) + ∫ u, ‖k2 u‖) with hM
  have hI0 : 0 ≤ ∫ u, ‖k u‖ := integral_nonneg fun _ => norm_nonneg _
  have hI1 : 0 ≤ ∫ u, ‖k1 u‖ := integral_nonneg fun _ => norm_nonneg _
  have hI2 : 0 ≤ ∫ u, ‖k2 u‖ := integral_nonneg fun _ => norm_nonneg _
  have hM0 : 0 ≤ M := by positivity
  -- the truncations
  have hjC2 : ∀ j, ContDiff ℝ 2 (trunc k j) := fun j => trunc_contDiff hkC2 j
  have hjc : ∀ j, HasCompactSupport (trunc k j) := fun j => trunc_hasCompactSupport k j
  -- uniform bound for the transforms of the truncations
  have hbound : ∀ (j : ℕ) (τ : ℝ), ‖paperFT (trunc k j) τ‖ * (1 + τ ^ 2) ≤ M := by
    intro j τ
    obtain ⟨Lam, hLam⟩ := exists_abs_le_of_hasCompactSupport (hjc j)
    refine (norm_paperFT_mul_one_add_sq_le (hjC2 j) hLam τ).trans ?_
    have hj : (0 : ℝ) < (j : ℝ) + 1 := by positivity
    have hj1 : (1 : ℝ) ≤ (j : ℝ) + 1 := by linarith [Nat.cast_nonneg (α := ℝ) j]
    have h1 : ∫ u, ‖trunc k j u‖ ≤ ∫ u, ‖k u‖ :=
      integral_mono_of_nonneg (Eventually.of_forall fun _ => norm_nonneg _) hki.norm
        (Eventually.of_forall (norm_trunc_le k j))
    have hdom : Integrable (fun u => B₂ * ‖k u‖ + 2 * B₁ * ‖k1 u‖ + ‖k2 u‖) :=
      ((hki.norm.const_mul B₂).add (hk1i.norm.const_mul (2 * B₁))).add hk2i.norm
    have h2 : ∫ u, ‖deriv (deriv (trunc k j)) u‖
        ≤ ∫ u, (B₂ * ‖k u‖ + 2 * B₁ * ‖k1 u‖ + ‖k2 u‖) := by
      refine integral_mono_of_nonneg (Eventually.of_forall fun _ => norm_nonneg _) hdom
        (Eventually.of_forall fun u => ?_)
      simp only
      rw [deriv_deriv_trunc hk1 hk2]
      have e1 : ‖((deriv (deriv (bump0 : ℝ → ℝ)) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1) ^ 2 : ℝ) : ℂ)‖
          ≤ B₂ := by
        rw [Complex.norm_real, Real.norm_eq_abs, abs_div, abs_of_pos (by positivity : (0:ℝ) < ((j : ℝ) + 1) ^ 2)]
        refine (div_le_self (abs_nonneg _) ?_).trans (hB₂ _)
        nlinarith
      have e2 : ‖((deriv (bump0 : ℝ → ℝ) (u / ((j : ℝ) + 1)) / ((j : ℝ) + 1) : ℝ) : ℂ)‖ ≤ B₁ := by
        rw [Complex.norm_real, Real.norm_eq_abs, abs_div, abs_of_pos hj]
        exact (div_le_self (abs_nonneg _) hj1).trans (hB₁ _)
      have e3 : ‖(cutoff j u : ℂ)‖ ≤ 1 := by
        rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (cutoff_nonneg j u)]
        exact cutoff_le_one j u
      calc ‖_ * k u + 2 * _ * k1 u + (cutoff j u : ℂ) * k2 u‖
          ≤ ‖_ * k u‖ + ‖2 * _ * k1 u‖ + ‖(cutoff j u : ℂ) * k2 u‖ := norm_add₃_le
        _ ≤ B₂ * ‖k u‖ + 2 * B₁ * ‖k1 u‖ + ‖k2 u‖ := by
            rw [norm_mul, norm_mul, norm_mul, norm_mul]
            have n2 : ‖(2 : ℂ)‖ = 2 := by norm_num
            rw [n2]
            gcongr
            exact mul_le_of_le_one_left (norm_nonneg _) e3
    have h3 : ∫ u, (B₂ * ‖k u‖ + 2 * B₁ * ‖k1 u‖ + ‖k2 u‖)
        = B₂ * (∫ u, ‖k u‖) + 2 * B₁ * (∫ u, ‖k1 u‖) + ∫ u, ‖k2 u‖ := by
      have e := integral_add ((hki.norm.const_mul B₂).add (hk1i.norm.const_mul (2 * B₁))) hk2i.norm
      have e' := integral_add (hki.norm.const_mul B₂) (hk1i.norm.const_mul (2 * B₁))
      simp only [Pi.add_apply] at e e'
      rw [e, e', integral_const_mul, integral_const_mul]
    rw [hM]
    linarith
  -- pointwise convergence of the transforms on the real line
  have hFTlim : ∀ τ : ℝ, Tendsto (fun j : ℕ => paperFT (trunc k j) τ) atTop (𝓝 (paperFT k τ)) := by
    intro τ
    unfold paperFT
    refine tendsto_integral_of_dominated_convergence (fun u => ‖k u‖) ?_ hki.norm ?_ ?_
    · intro j
      exact ((hjC2 j).continuous.mul (by fun_prop)).aestronglyMeasurable
    · intro j
      refine Eventually.of_forall fun u => ?_
      rw [norm_mul, RaceLib.norm_cexp_I_mul_ofReal_mul_ofReal τ u, mul_one]
      exact norm_trunc_le k j u
    · refine Eventually.of_forall fun u => ?_
      exact (tendsto_trunc k u).mul_const _
  have hbound_lim : ∀ τ : ℝ, ‖paperFT k τ‖ * (1 + τ ^ 2) ≤ M := fun τ =>
    le_of_tendsto' (((hFTlim τ).norm).mul_const _) fun j => hbound j τ
  have hbd : ∀ (f : ℝ → ℂ), (∀ τ : ℝ, ‖f τ‖ * (1 + τ ^ 2) ≤ M) → ∀ τ : ℝ, ‖f τ‖ ≤ M / (1 + τ ^ 2) := by
    intro f hf τ
    rw [le_div_iff₀ (by positivity)]
    exact hf τ
  -- GRH: ordinates are real
  have hre : ∀ ρ : Z.carrier, (ρ : ℂ).re = 1 / 2 := fun ρ => hGRH ρ ρ.2.2.1 ρ.2.2.2 ρ.2.1
  have hgam : ∀ ρ : Z.carrier, gammaOf (ρ : ℂ) = (((ρ : ℂ).im : ℝ) : ℂ) := fun ρ =>
    RaceLib.gammaOf_eq_im (hre ρ)
  have hnsq : ∀ ρ : Z.carrier, Complex.normSq (gammaOf (ρ : ℂ)) = (ρ : ℂ).im ^ 2 := by
    intro ρ
    rw [hgam ρ, Complex.normSq_ofReal, sq]
  -- the summable majorant on the zero side
  obtain ⟨A₀, hA₀, hloc⟩ := localCountChi hq hprim
  have hloc' : ∀ t : ℝ, (Z.N t (t + 1) : ℝ) ≤ A₀ * Real.log (|t| + 3) := by
    intro t; rw [hZ, LZeros_N]; exact hloc t
  have hsumZ : Summable (fun ρ : Z.carrier =>
      M * ((Z.mult ρ : ℝ) / (1 + Complex.normSq (gammaOf (ρ : ℂ))))) :=
    (Zeta23.WeilEF.zero_sum_inv_sq_gen Z hA₀ hloc').mul_left M
  have hzb : ∀ (f : ℝ → ℂ), (∀ τ : ℝ, ‖f τ‖ * (1 + τ ^ 2) ≤ M) → ∀ ρ : Z.carrier,
      ‖(Z.mult ρ : ℂ) * f ((ρ : ℂ).im)‖
        ≤ M * ((Z.mult ρ : ℝ) / (1 + Complex.normSq (gammaOf (ρ : ℂ)))) := by
    intro f hf ρ
    rw [norm_mul, Complex.norm_natCast, hnsq ρ]
    calc (Z.mult ρ : ℝ) * ‖f ((ρ : ℂ).im)‖ ≤ (Z.mult ρ : ℝ) * (M / (1 + (ρ : ℂ).im ^ 2)) := by
          gcongr
          exact hbd f hf _
      _ = M * ((Z.mult ρ : ℝ) / (1 + (ρ : ℂ).im ^ 2)) := by ring
  have hsum : Summable (fun ρ : Z.carrier => (Z.mult ρ : ℂ) * paperFT k (gammaOf ρ)) := by
    refine Summable.of_norm_bounded hsumZ fun ρ => ?_
    have h := hzb (fun τ => paperFT k τ) hbound_lim ρ
    rw [← hgam ρ] at h
    exact h
  refine ⟨hsum, ?_⟩
  -- the explicit formula for the truncations
  have hEF : ∀ j : ℕ, ∑' ρ : Z.carrier, (Z.mult ρ : ℂ) * paperFT (trunc k j) (gammaOf ρ)
      = literatureRHSchi (parity χ) q (coeff χ) (trunc k j) := fun j =>
    (EF_lit_chi_L hq hprim (trunc k j) (hjC2 j) (hjc j)).2
  -- zero side converges
  have hzero : Tendsto (fun j : ℕ => ∑' ρ : Z.carrier,
      (Z.mult ρ : ℂ) * paperFT (trunc k j) (gammaOf ρ)) atTop
      (𝓝 (∑' ρ : Z.carrier, (Z.mult ρ : ℂ) * paperFT k (gammaOf ρ))) := by
    refine tendsto_tsum_of_dominated_convergence hsumZ ?_ ?_
    · intro ρ
      rw [hgam ρ]
      exact (hFTlim _).const_mul _
    · refine Eventually.of_forall fun j ρ => ?_
      have h := hzb (fun τ => paperFT (trunc k j) τ) (hbound j) ρ
      rw [← hgam ρ] at h
      exact h
  -- prime side converges
  have hcoef : ∀ n : ℕ, ‖coeff χ n‖ ≤ 1 := fun n => DirichletCharacter.norm_le_one χ _
  have hprime : Tendsto (fun j : ℕ => ∑' n : ℕ,
      ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
        * (coeff χ n * trunc k j (Real.log n)
          + (starRingEnd ℂ) (coeff χ n) * trunc k j (-Real.log n))) atTop
      (𝓝 (∑' n : ℕ, ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
        * (coeff χ n * k (Real.log n) + (starRingEnd ℂ) (coeff χ n) * k (-Real.log n)))) := by
    refine tendsto_tsum_of_dominated_convergence hps ?_ ?_
    · intro n
      exact (((tendsto_trunc k _).const_mul _).add ((tendsto_trunc k _).const_mul _)).const_mul _
    · refine Eventually.of_forall fun j n => ?_
      have hΛ : 0 ≤ ArithmeticFunction.vonMangoldt n / Real.sqrt n :=
        div_nonneg ArithmeticFunction.vonMangoldt_nonneg (Real.sqrt_nonneg _)
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hΛ]
      refine mul_le_mul_of_nonneg_left ?_ hΛ
      refine (norm_add_le _ _).trans (add_le_add ?_ ?_)
      · rw [norm_mul]
        calc ‖coeff χ n‖ * ‖trunc k j (Real.log n)‖ ≤ 1 * ‖k (Real.log n)‖ :=
              mul_le_mul (hcoef n) (norm_trunc_le k j _) (norm_nonneg _) zero_le_one
          _ = _ := one_mul _
      · rw [norm_mul, RCLike.norm_conj]
        calc ‖coeff χ n‖ * ‖trunc k j (-Real.log n)‖ ≤ 1 * ‖k (-Real.log n)‖ :=
              mul_le_mul (hcoef n) (norm_trunc_le k j _) (norm_nonneg _) zero_le_one
          _ = _ := one_mul _
  -- archimedean side converges
  obtain ⟨K₂, hK₂, hmu⟩ := abs_muq_le_of_gammaFactsChi
    (GammaChi.gammaFactsChi (parity_le_one' χ) hq1) hq1
  have hW : ∀ r : ℝ, |(Complex.digamma (1 / 4 + (parity χ : ℂ) / 2 + I * r / 2)).re
      + Real.log (q / Real.pi)| ≤ 2 * Real.pi * (K₂ * (1 + |r|) ^ (1 / 2 : ℝ)) := by
    intro r
    have h := hmu r
    unfold muq at h
    rw [abs_mul, abs_of_pos (by positivity : (0:ℝ) < 1 / (2 * Real.pi)), add_comm] at h
    have hpi : (0:ℝ) < 2 * Real.pi := by positivity
    calc |_| = 2 * Real.pi * (1 / (2 * Real.pi) * |_|) := by field_simp
      _ ≤ _ := mul_le_mul_of_nonneg_left h hpi.le
  have harch : Tendsto (fun j : ℕ => ∫ r : ℝ, paperFT (trunc k j) r
      * ((((Complex.digamma (1 / 4 + (parity χ : ℂ) / 2 + I * r / 2)).re
        + Real.log (q / Real.pi) : ℝ)) : ℂ)) atTop
      (𝓝 (∫ r : ℝ, paperFT k r
      * ((((Complex.digamma (1 / 4 + (parity χ : ℂ) / 2 + I * r / 2)).re
        + Real.log (q / Real.pi) : ℝ)) : ℂ))) := by
    refine tendsto_integral_of_dominated_convergence
      (fun r : ℝ => 4 * Real.pi * M * K₂ * (1 + ‖r‖) ^ (-(3 / 2) : ℝ)) ?_ ?_ ?_ ?_
    · intro j
      exact (RaceLib.integrable_raceFT_mul_archWeight χ (hjC2 j) (hjc j)).aestronglyMeasurable
    · exact (integrable_one_add_norm (by rw [Module.finrank_self]; norm_num)).const_mul _
    · intro j
      refine Eventually.of_forall fun r => ?_
      rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs]
      have hb : 0 < 1 + |r| := by linarith [abs_nonneg r]
      have h1 := hbound j r
      have hsq : (1 + |r|) ^ (2 : ℝ) ≤ 2 * (1 + r ^ 2) := by
        rw [Real.rpow_two]
        nlinarith [sq_abs r, sq_nonneg (1 - |r|)]
      have key : (1 + |r|) ^ (1 / 2 : ℝ) ≤ 2 * (1 + r ^ 2) * (1 + |r|) ^ (-(3 / 2) : ℝ) := by
        rw [show (1 / 2 : ℝ) = 2 + (-(3 / 2)) by norm_num, Real.rpow_add hb]
        exact mul_le_mul_of_nonneg_right hsq (Real.rpow_nonneg hb.le _)
      have hrp : 0 ≤ (1 + |r|) ^ (-(3 / 2) : ℝ) := Real.rpow_nonneg hb.le _
      calc ‖paperFT (trunc k j) r‖ * |_|
          ≤ ‖paperFT (trunc k j) r‖ * (2 * Real.pi * (K₂ * (1 + |r|) ^ (1 / 2 : ℝ))) := by
            gcongr
            exact hW r
        _ ≤ ‖paperFT (trunc k j) r‖
              * (2 * Real.pi * (K₂ * (2 * (1 + r ^ 2) * (1 + |r|) ^ (-(3 / 2) : ℝ)))) := by
            gcongr
        _ = (‖paperFT (trunc k j) r‖ * (1 + r ^ 2))
              * (4 * Real.pi * K₂ * (1 + |r|) ^ (-(3 / 2) : ℝ)) := by ring
        _ ≤ M * (4 * Real.pi * K₂ * (1 + |r|) ^ (-(3 / 2) : ℝ)) := by
            gcongr
        _ = 4 * Real.pi * M * K₂ * (1 + |r|) ^ (-(3 / 2) : ℝ) := by ring
    · refine Eventually.of_forall fun r => ?_
      exact (hFTlim r).mul_const _
  have hRHS : Tendsto (fun j : ℕ => literatureRHSchi (parity χ) q (coeff χ) (trunc k j)) atTop
      (𝓝 (literatureRHSchi (parity χ) q (coeff χ) k)) := by
    unfold literatureRHSchi
    exact hprime.neg.add (harch.const_mul _)
  exact tendsto_nhds_unique hzero (hRHS.congr fun j => (hEF j).symm)

end ChebLib
