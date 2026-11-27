/-
ChebLib/Cheb/Final.lean — Chebyshev's assertion under GRH for `L(s,χ₋₄)`.

  * `abs_psiAbel_le`  (PsiAbel.lean):  `|Σ_n Λ(n) χ₋₄(n) e^{−n/x}| ≤ (G + ε) √x`,
    `G = Σ_ρ m_ρ |Γ(1/2+iγ_ρ)|`;
  * `gammaSum_lt`     (ZeroSum.lean):  `G < √π/2`;
  * prime powers: `Σ_n c(n) e^{−n/x} ~ (√π/2) √x` for `c(n) = Λ(n) χ₋₄(n) [n not prime]`
    (prime squares, by the prime number theorem and the Abelian theorem of AbelSum.lean);
  * hence `Σ_p χ₋₄(p) (log p) e^{−p/x} ≤ −c √x` for large `x`;
  * hence (LogRemoval.lean) `Σ_p χ₋₄(p) e^{−p/x} → −∞`.
-/
import ChebLib.Cheb.ZeroSum
import ChebLib.Cheb.LogRemoval
import ChebLib.Race.PNT

open scoped BigOperators
open Complex MeasureTheory Set Filter Topology

noncomputable section

namespace ChebLib

/-! ### Primes, prime powers -/

/-- `a_n = χ₋₄(n)` if `n` is prime, `0` otherwise. -/
def chebA (n : ℕ) : ℝ := if n.Prime then X4 n else 0

theorem abs_chebA_le_one (n : ℕ) : |chebA n| ≤ 1 := by
  unfold chebA
  split_ifs
  · exact abs_X4_le_one n
  · simp

theorem chebA_eq_zero_of_lt_three {n : ℕ} (h : n < 3) : chebA n = 0 := by
  unfold chebA
  interval_cases n
  · simp [Nat.not_prime_zero]
  · simp [Nat.not_prime_one]
  · rw [if_pos Nat.prime_two]
    exact RaceLib.chi4_natCast_two

theorem vM_mul_X4_eq (n : ℕ) : vM n * X4 n = chebA n * Real.log n + RaceLib.npCoeff n := by
  unfold chebA RaceLib.npCoeff
  by_cases hp : n.Prime
  · rw [if_pos hp, if_pos hp, add_zero]
    unfold vM
    rw [ArithmeticFunction.vonMangoldt_apply_prime hp]
    ring
  · rw [if_neg hp, if_neg hp, zero_mul, zero_add]
    rfl

theorem abs_npCoeff_le (n : ℕ) : |RaceLib.npCoeff n| ≤ 1 * (n : ℝ) := by
  unfold RaceLib.npCoeff
  split_ifs
  · simp
  · exact abs_vM_mul_X4_le n

theorem abs_chebA_mul_log_le (n : ℕ) : |chebA n * Real.log n| ≤ 1 * (n : ℝ) := by
  rw [abs_mul, one_mul]
  have h0 : 0 ≤ Real.log n := Real.log_natCast_nonneg n
  rw [abs_of_nonneg h0]
  calc |chebA n| * Real.log n ≤ 1 * Real.log n :=
        mul_le_mul_of_nonneg_right (abs_chebA_le_one n) h0
    _ ≤ n := by rw [one_mul]; exact Real.log_le_self (Nat.cast_nonneg n)

/-- `Φ(x) = Σ_p χ₋₄(p) (log p) e^{−p/x}`, as a sum over `n`. -/
def phiAbel (x : ℝ) : ℝ := ∑' n : ℕ, chebA n * Real.log n * Real.exp (-(n / x))

/-- The contribution of the proper prime powers, `Σ_n c(n) e^{−n/x}`. -/
def npAbel (x : ℝ) : ℝ := ∑' n : ℕ, RaceLib.npCoeff n * Real.exp (-(n / x))

theorem psiAbel_eq {x : ℝ} (hx : 0 < x) : psiAbel x = phiAbel x + npAbel x := by
  unfold psiAbel phiAbel npAbel
  rw [← (summable_mul_exp_neg_div _ 1 abs_chebA_mul_log_le hx).tsum_add
    (summable_mul_exp_neg_div _ 1 abs_npCoeff_le hx)]
  refine tsum_congr fun n => ?_
  rw [vM_mul_X4_eq]
  ring

/-- Prime squares: `Σ_n c(n) e^{−n/x} / √x → √π/2` (prime number theorem + Abelian theorem). -/
theorem tendsto_npAbel :
    Tendsto (fun x : ℝ => npAbel x / Real.sqrt x) atTop (𝓝 (Real.sqrt Real.pi / 2)) :=
  tendsto_abel_sum_div_sqrt RaceLib.npCoeff 1 abs_npCoeff_le
    (RaceLib.tendsto_sum_npCoeff_div_sqrt RaceLib.pntTheta)

/-! ### The logarithmically weighted sum over primes -/

/-- **Two-sided form.**  Under GRH for `L(s,χ₋₄)`, for every `ε > 0` and all large `x`:
`|Σ_p χ₋₄(p) (log p) e^{−p/x} + (√π/2) √x| ≤ (G + ε) √x`, `G = Σ_ρ m_ρ |Γ(1/2+iγ_ρ)|`. -/
theorem abs_phiAbel_add_le (hGRH : GRHFor chi4) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ x : ℝ in atTop,
      |phiAbel x + Real.sqrt Real.pi / 2 * Real.sqrt x| ≤ (gammaSum + ε) * Real.sqrt x := by
  have h1 := abs_psiAbel_le hGRH (half_pos hε)
  have h2 := tendsto_npAbel
  rw [Metric.tendsto_nhds] at h2
  filter_upwards [h1, h2 (ε / 2) (half_pos hε), eventually_gt_atTop (0 : ℝ)] with x hx1 hx2 hx0
  have hs : 0 < Real.sqrt x := Real.sqrt_pos.mpr hx0
  rw [Real.dist_eq] at hx2
  have h3 : |npAbel x - Real.sqrt Real.pi / 2 * Real.sqrt x| ≤ ε / 2 * Real.sqrt x := by
    have e : npAbel x - Real.sqrt Real.pi / 2 * Real.sqrt x
        = (npAbel x / Real.sqrt x - Real.sqrt Real.pi / 2) * Real.sqrt x := by
      field_simp
    rw [e, abs_mul, abs_of_pos hs]
    exact mul_le_mul_of_nonneg_right hx2.le hs.le
  have e : phiAbel x + Real.sqrt Real.pi / 2 * Real.sqrt x
      = psiAbel x - (npAbel x - Real.sqrt Real.pi / 2 * Real.sqrt x) := by
    rw [psiAbel_eq hx0]; ring
  rw [e]
  calc |psiAbel x - (npAbel x - Real.sqrt Real.pi / 2 * Real.sqrt x)|
      ≤ |psiAbel x| + |npAbel x - Real.sqrt Real.pi / 2 * Real.sqrt x| := abs_sub _ _
    _ ≤ (gammaSum + ε / 2) * Real.sqrt x + ε / 2 * Real.sqrt x := add_le_add hx1 h3
    _ = (gammaSum + ε) * Real.sqrt x := by ring

/-- **One-sided form.**  Under GRH for `L(s,χ₋₄)` there are `c > 0` and `x₀` with
`Σ_p χ₋₄(p) (log p) e^{−p/x} ≤ −c √x` for all `x ≥ x₀`. -/
theorem phiAbel_le (hGRH : GRHFor chi4) :
    ∃ c : ℝ, 0 < c ∧ ∃ x₀ : ℝ, ∀ y : ℝ, x₀ ≤ y → phiAbel y ≤ -c * Real.sqrt y := by
  have hlt := gammaSum_lt hGRH
  set c : ℝ := (Real.sqrt Real.pi / 2 - gammaSum) / 2 with hc
  have hc0 : 0 < c := by rw [hc]; linarith
  obtain ⟨x₀, hx₀⟩ := eventually_atTop.mp (abs_phiAbel_add_le hGRH hc0)
  refine ⟨c, hc0, x₀, fun y hy => ?_⟩
  have h := (abs_le.mp (hx₀ y hy)).2
  have hs : 0 ≤ Real.sqrt y := Real.sqrt_nonneg y
  have e : (gammaSum + c) * Real.sqrt y - Real.sqrt Real.pi / 2 * Real.sqrt y
      = -c * Real.sqrt y := by rw [hc]; ring
  linarith

/-! ### Chebyshev's assertion -/

/-- Sum over `n` with the prime indicator. -/
theorem tendsto_chebA (hGRH : GRHFor chi4) :
    Tendsto (fun x : ℝ => ∑' n : ℕ, chebA n * Real.exp (-(n / x))) atTop atBot := by
  obtain ⟨c, hc0, x₀, hx₀⟩ := phiAbel_le hGRH
  exact tendsto_atBot_of_log_weighted_le chebA abs_chebA_le_one
    (fun n hn => chebA_eq_zero_of_lt_three hn) c x₀ hc0 hx₀

theorem tsum_primes_eq (f : ℕ → ℝ) :
    ∑' p : Nat.Primes, f p = ∑' n : ℕ, if n.Prime then f n else 0 := by
  have h : ∑' p : Nat.Primes, f p = ∑' n : ℕ, {n : ℕ | n.Prime}.indicator f n :=
    tsum_subtype {n : ℕ | n.Prime} f
  rw [h]
  refine tsum_congr fun n => ?_
  by_cases hp : n.Prime
  · rw [Set.indicator_of_mem (show n ∈ {n : ℕ | n.Prime} from hp), if_pos hp]
  · rw [Set.indicator_of_notMem (show n ∉ {n : ℕ | n.Prime} from hp), if_neg hp]

theorem tsum_primes_chi4_exp (x : ℝ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)
      = ∑' n : ℕ, chebA n * Real.exp (-(n / x)) := by
  rw [tsum_primes_eq (fun n => ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * Real.exp (-(n : ℝ) / x))]
  refine tsum_congr fun n => ?_
  unfold chebA X4
  split_ifs
  · rw [neg_div]
  · rw [zero_mul]

theorem tsum_primes_chi4_log_exp (x : ℝ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
        * Real.exp (-((p : ℕ) : ℝ) / x)
      = phiAbel x := by
  rw [tsum_primes_eq (fun n => ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * Real.log n
    * Real.exp (-(n : ℝ) / x))]
  unfold phiAbel
  refine tsum_congr fun n => ?_
  unfold chebA X4
  split_ifs
  · rw [neg_div]
  · rw [zero_mul, zero_mul]

/-- The series in Chebyshev's assertion converges absolutely for every `x > 0`. -/
theorem summable_primes_chi4_exp {x : ℝ} (hx : 0 < x) :
    Summable (fun p : Nat.Primes =>
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  have h : Summable (fun n : ℕ => ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * Real.exp (-(n : ℝ) / x)) := by
    have := summable_mul_exp_neg_div (fun n : ℕ => X4 n) 1 (fun n => by
      rw [one_mul]
      rcases Nat.eq_zero_or_pos n with rfl | hn
      · simp [X4]
      · exact (abs_X4_le_one n).trans (by exact_mod_cast hn)) hx
    refine this.congr fun n => ?_
    rw [neg_div]
    rfl
  exact h.subtype {n : ℕ | n.Prime}

/-- **Chebyshev's assertion under GRH for `L(s,χ₋₄)`**:
`Σ_p χ₋₄(p) e^{−p/x} → −∞` as `x → +∞`. -/
theorem chebyshev_bias (hGRH : GRHFor chi4) :
    Tendsto (fun x : ℝ => ∑' p : Nat.Primes,
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)) atTop atBot := by
  refine (tendsto_chebA hGRH).congr fun x => ?_
  exact (tsum_primes_chi4_exp x).symm

end ChebLib
