/-
ChebLib/Cheb/Numeric.lean — the numerical inequality behind the zero-sum bound.

With `g(u) = (5/2) e^{−(3/2)|u|} − (3/2) e^{−(5/2)|u|}`, the explicit formula for `L(s,χ₋₄)` reads
  `Σ_ρ m_ρ ĝ(γ_ρ) = (7/2 − γ − 3 log 2 − log π) + Σ_n Λ(n) (−χ₋₄(n)) (5n − 3)/n³`,
and the right-hand side (numerically `≈ 0.0462`) must be shown to be `< 16/15`.
-/
import Mathlib

open Real

noncomputable section

namespace ChebLib

/-! ### Elementary bounds on the summand -/

/-- The weight `(5n − 3)/n³` is nonnegative (also for `n = 0`, where it is `0`). -/
private lemma majCoef_nonneg (n : ℕ) : 0 ≤ (5 * (n : ℝ) - 3) / (n : ℝ) ^ 3 := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp
  · have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
    apply div_nonneg
    · linarith
    · positivity

/-- The weight `(5n − 3)/n³` is at most `5/n²`. -/
private lemma majCoef_le (n : ℕ) (hn : 1 ≤ n) :
    (5 * (n : ℝ) - 3) / (n : ℝ) ^ 3 ≤ 5 / (n : ℝ) ^ 2 := by
  have h1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have h0 : (0 : ℝ) < n := by linarith
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [sq_nonneg (n : ℝ)]

/-- The values of `χ₄` on naturals, as real numbers. -/
private lemma chi4_cases (n : ℕ) :
    (n % 4 = 3 ∧ ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = -1) ∨
    (n % 4 ≠ 3 ∧ (((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = 0 ∨ ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = 1)) := by
  rw [ZMod.χ₄_nat_eq_if_mod_four]
  have h4 : n % 4 < 4 := Nat.mod_lt _ (by norm_num)
  by_cases h2 : n % 2 = 0
  · right
    refine ⟨by omega, Or.inl ?_⟩
    simp [h2]
  · by_cases h1 : n % 4 = 1
    · right
      refine ⟨by omega, Or.inr ?_⟩
      simp [h2, h1]
    · left
      refine ⟨by omega, ?_⟩
      simp [h2, h1]

/-- The summand. -/
private def majTerm (n : ℕ) : ℝ :=
  ArithmeticFunction.vonMangoldt n
    * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)

/-- Terms outside the progression `3 mod 4` are nonpositive. -/
private lemma majTerm_nonpos (n : ℕ) (hn : n % 4 ≠ 3) : majTerm n ≤ 0 := by
  have hΛ : 0 ≤ ArithmeticFunction.vonMangoldt n := ArithmeticFunction.vonMangoldt_nonneg
  have hc := majCoef_nonneg n
  unfold majTerm
  rcases chi4_cases n with ⟨h, _⟩ | ⟨_, h | h⟩
  · exact absurd h hn
  · rw [h]; simp
  · rw [h]
    have := mul_nonneg hΛ hc
    nlinarith

/-- Terms in the progression `3 mod 4` are at most `5 log n / n²`. -/
private lemma majTerm_le (n : ℕ) (hn : 1 ≤ n) :
    majTerm n ≤ 5 * Real.log n / (n : ℝ) ^ 2 := by
  have hΛ : 0 ≤ ArithmeticFunction.vonMangoldt n := ArithmeticFunction.vonMangoldt_nonneg
  have hΛ' : ArithmeticFunction.vonMangoldt n ≤ Real.log n :=
    ArithmeticFunction.vonMangoldt_le_log
  have hc := majCoef_nonneg n
  have hc' := majCoef_le n hn
  have hlog : 0 ≤ Real.log n := Real.log_natCast_nonneg n
  have h5 : (0 : ℝ) ≤ 5 / (n : ℝ) ^ 2 := by positivity
  have hmain : ArithmeticFunction.vonMangoldt n * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)
      ≤ Real.log n * (5 / (n : ℝ) ^ 2) := mul_le_mul hΛ' hc' hc hlog
  have hnn : 0 ≤ Real.log n * (5 / (n : ℝ) ^ 2) := mul_nonneg hlog h5
  have hmn : 0 ≤ ArithmeticFunction.vonMangoldt n * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3) :=
    mul_nonneg hΛ hc
  have e : 5 * Real.log n / (n : ℝ) ^ 2 = Real.log n * (5 / (n : ℝ) ^ 2) := by ring
  rw [e]
  unfold majTerm
  rcases chi4_cases n with ⟨_, h⟩ | ⟨_, h | h⟩
  · rw [h]; linarith
  · rw [h]; simpa using hnn
  · rw [h]; linarith

/-- Absolute value bound for the summand. -/
private lemma abs_majTerm_le (n : ℕ) (hn : 1 ≤ n) :
    |majTerm n| ≤ 5 * (ArithmeticFunction.vonMangoldt n / (n : ℝ) ^ 2) := by
  have hΛ : 0 ≤ ArithmeticFunction.vonMangoldt n := ArithmeticFunction.vonMangoldt_nonneg
  have hc := majCoef_nonneg n
  have hc' := majCoef_le n hn
  have hmain : ArithmeticFunction.vonMangoldt n * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)
      ≤ ArithmeticFunction.vonMangoldt n * (5 / (n : ℝ) ^ 2) :=
    mul_le_mul_of_nonneg_left hc' hΛ
  have hmn : 0 ≤ ArithmeticFunction.vonMangoldt n * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3) :=
    mul_nonneg hΛ hc
  have e : 5 * (ArithmeticFunction.vonMangoldt n / (n : ℝ) ^ 2)
      = ArithmeticFunction.vonMangoldt n * (5 / (n : ℝ) ^ 2) := by ring
  rw [e, abs_le]
  unfold majTerm
  rcases chi4_cases n with ⟨_, h⟩ | ⟨_, h | h⟩
  · rw [h]; constructor <;> linarith
  · rw [h]; constructor <;> simp <;> linarith
  · rw [h]; constructor <;> linarith

private lemma summable_majTerm : Summable majTerm := by
  have h := ArithmeticFunction.LSeriesSummable_vonMangoldt (s := 2) (by norm_num)
  have h' := (h.norm).mul_left 5
  refine Summable.of_norm_bounded h' ?_
  intro n
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [majTerm]
  · rw [Real.norm_eq_abs, LSeries.norm_term_eq]
    have hne : n ≠ 0 := Nat.pos_iff_ne_zero.mp hn
    have hΛ : 0 ≤ ArithmeticFunction.vonMangoldt n := ArithmeticFunction.vonMangoldt_nonneg
    simp only [hne, if_false]
    have e1 : ‖(ArithmeticFunction.vonMangoldt n : ℂ)‖ = ArithmeticFunction.vonMangoldt n := by
      rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hΛ]
    have e2 : ((n : ℝ)) ^ (2 : ℂ).re = (n : ℝ) ^ 2 := by
      have : (2 : ℂ).re = ((2 : ℕ) : ℝ) := by simp
      rw [this, Real.rpow_natCast]
    rw [e1, e2]
    exact abs_majTerm_le n hn

/-- The prime side of the explicit formula for the majorant converges absolutely. -/
theorem summable_majPrime :
    Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n
      * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)) :=
  summable_majTerm

/-! ### The telescoping tail bound -/

/-- A general telescoping bound for a summable series. -/
private lemma tsum_le_of_telescope {g P : ℕ → ℝ} (hs : Summable g) (hP : ∀ N, 0 ≤ P N)
    (h : ∀ N, g N ≤ P N - P (N + 1)) : ∑' n, g n ≤ P 0 := by
  have key : ∀ N, ∑ i ∈ Finset.range N, g i ≤ P 0 - P N := by
    intro N
    induction N with
    | zero => simp
    | succ N ih =>
      rw [Finset.sum_range_succ]
      linarith [h N]
  refine le_of_tendsto' hs.hasSum.tendsto_sum_nat (fun N => ?_)
  linarith [key N, hP N]

/-- The potential `D(t) = 5 (log t + 1) / (4 t)`, a primitive-type majorant for
`t ↦ 5 log t / t²` along steps of length `4`. -/
private def majPot (t : ℝ) : ℝ := 5 * (Real.log t + 1) / (4 * t)

private lemma majPot_nonneg (t : ℝ) (ht : 1 ≤ t) : 0 ≤ majPot t := by
  unfold majPot
  have := Real.log_nonneg ht
  positivity

/-- Key step inequality: `5 log(x+4)/(x+4)² ≤ D(x) − D(x+4)` for `x ≥ 7`. -/
private lemma majPot_step (x : ℝ) (hx : 7 ≤ x) :
    5 * Real.log (x + 4) / (x + 4) ^ 2 ≤ majPot x - majPot (x + 4) := by
  have hx0 : 0 < x := by linarith
  have hx4 : 0 < x + 4 := by linarith
  -- `L = log (x+4) ≥ 3 log 2 ≥ 11/7`
  have hL : 11 / 7 ≤ Real.log (x + 4) := by
    have h8 : Real.log 8 ≤ Real.log (x + 4) := Real.log_le_log (by norm_num) (by linarith)
    have e8 : Real.log 8 = 3 * Real.log 2 := by
      rw [show (8 : ℝ) = 2 ^ 3 by norm_num, Real.log_pow]; norm_num
    have := Real.log_two_gt_d9
    linarith
  -- `δ = log (x+4) − log x ≤ 4/x`
  obtain ⟨δ, hδ⟩ : ∃ δ, δ = Real.log (x + 4) - Real.log x := ⟨_, rfl⟩
  have hδx : δ * x ≤ 4 := by
    have h1 : Real.log ((x + 4) / x) ≤ (x + 4) / x - 1 :=
      Real.log_le_sub_one_of_pos (by positivity)
    rw [Real.log_div hx4.ne' hx0.ne'] at h1
    have h2 : (x + 4) / x - 1 = 4 / x := by field_simp; ring
    rw [h2, ← hδ] at h1
    calc δ * x ≤ 4 / x * x := mul_le_mul_of_nonneg_right h1 hx0.le
      _ = 4 := by field_simp
  have hδ7 : δ * 7 ≤ 4 := by
    by_cases hδ0 : 0 ≤ δ
    · nlinarith
    · linarith
  have hδxx : δ * x * x ≤ 4 * x := mul_le_mul_of_nonneg_right hδx hx0.le
  have hl : Real.log x = Real.log (x + 4) - δ := by linarith
  unfold majPot
  rw [hl, ← sub_nonneg]
  have e : 5 * (Real.log (x + 4) - δ + 1) / (4 * x) - 5 * (Real.log (x + 4 + 4 - 4) + 1) / (4 * (x + 4))
        - 5 * Real.log (x + 4) / (x + 4) ^ 2
      = 5 * (16 * Real.log (x + 4) + 4 * (x + 4) - δ * (x + 4) ^ 2) / (4 * x * (x + 4) ^ 2) := by
    rw [show x + 4 + 4 - 4 = x + 4 by ring]
    field_simp
    ring
  rw [show x + 4 + 4 - 4 = x + 4 by ring] at e
  rw [e]
  apply div_nonneg
  · nlinarith
  · positivity

/-- The potential along the tail: `P(N) = D(m)` with `m` the largest integer `≡ 3 mod 4`
that is `≤ N + 10`. -/
private def majP (N : ℕ) : ℝ := majPot ((4 * ((N + 7) / 4) + 3 : ℕ) : ℝ)

private lemma majTerm_tail_le (N : ℕ) : majTerm (N + 11) ≤ majP N - majP (N + 1) := by
  by_cases h : (N + 11) % 4 = 3
  · have e1 : 4 * ((N + 7) / 4) + 3 = N + 7 := by omega
    have e2 : 4 * ((N + 1 + 7) / 4) + 3 = N + 11 := by omega
    unfold majP
    rw [e1, e2]
    have h1 := majTerm_le (N + 11) (by omega)
    have h2 := majPot_step ((N : ℝ) + 7) (by have : (0 : ℝ) ≤ N := Nat.cast_nonneg N; linarith)
    have e3 : ((N + 11 : ℕ) : ℝ) = (N : ℝ) + 7 + 4 := by push_cast; ring
    have e4 : ((N + 7 : ℕ) : ℝ) = (N : ℝ) + 7 := by push_cast; ring
    rw [e3] at h1 ⊢
    rw [e4]
    linarith
  · have e : 4 * ((N + 1 + 7) / 4) + 3 = 4 * ((N + 7) / 4) + 3 := by omega
    unfold majP
    rw [e, sub_self]
    exact majTerm_nonpos _ h

private lemma majTerm_tail : ∑' n, majTerm (n + 11) ≤ majPot 7 := by
  have hs : Summable (fun n => majTerm (n + 11)) := (summable_nat_add_iff 11).2 summable_majTerm
  have h := tsum_le_of_telescope (P := majP) hs
    (fun N => majPot_nonneg _ (by
      have : 1 ≤ 4 * ((N + 7) / 4) + 3 := by omega
      exact_mod_cast this))
    majTerm_tail_le
  have e : majP 0 = majPot 7 := by
    unfold majP
    norm_num
  rw [e] at h
  exact h

/-! ### The head of the series -/

private lemma majTerm_three : majTerm 3 = Real.log 3 * (12 / 27) := by
  unfold majTerm
  rw [ArithmeticFunction.vonMangoldt_apply_prime (by norm_num : Nat.Prime 3),
    ZMod.χ₄_nat_three_mod_four (by norm_num : 3 % 4 = 3)]
  norm_num

private lemma majTerm_five : majTerm 5 = -(Real.log 5 * (22 / 125)) := by
  unfold majTerm
  rw [ArithmeticFunction.vonMangoldt_apply_prime (by norm_num : Nat.Prime 5),
    ZMod.χ₄_nat_one_mod_four (by norm_num : 5 % 4 = 1)]
  norm_num

private lemma majTerm_seven : majTerm 7 = Real.log 7 * (32 / 343) := by
  unfold majTerm
  rw [ArithmeticFunction.vonMangoldt_apply_prime (by norm_num : Nat.Prime 7),
    ZMod.χ₄_nat_three_mod_four (by norm_num : 7 % 4 = 3)]
  norm_num

private lemma majTerm_head :
    ∑ n ∈ Finset.range 11, majTerm n
      ≤ Real.log 3 * (12 / 27) - Real.log 5 * (22 / 125) + Real.log 7 * (32 / 343) := by
  simp only [Finset.sum_range_succ, Finset.sum_range_zero]
  have h0 := majTerm_nonpos 0 (by norm_num)
  have h1 := majTerm_nonpos 1 (by norm_num)
  have h2 := majTerm_nonpos 2 (by norm_num)
  have h4 := majTerm_nonpos 4 (by norm_num)
  have h6 := majTerm_nonpos 6 (by norm_num)
  have h8 := majTerm_nonpos 8 (by norm_num)
  have h9 := majTerm_nonpos 9 (by norm_num)
  have h10 := majTerm_nonpos 10 (by norm_num)
  rw [majTerm_three, majTerm_five, majTerm_seven]
  linarith

/-! ### Logarithms of small numbers in terms of `log 2` -/

private lemma log_three_ge : 2 * Real.log 2 - 1 / 3 ≤ Real.log 3 := by
  have h : Real.log (4 / 3) ≤ 4 / 3 - 1 := Real.log_le_sub_one_of_pos (by norm_num)
  rw [Real.log_div (by norm_num) (by norm_num)] at h
  have e4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  linarith

private lemma log_three_le : Real.log 3 ≤ 2 * Real.log 2 - 1 / 4 := by
  have h : Real.log (3 / 4) ≤ 3 / 4 - 1 := Real.log_le_sub_one_of_pos (by norm_num)
  rw [Real.log_div (by norm_num) (by norm_num)] at h
  have e4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  linarith

private lemma log_five_ge : 2 * Real.log 2 ≤ Real.log 5 := by
  have h : Real.log 4 ≤ Real.log 5 := Real.log_le_log (by norm_num) (by norm_num)
  have e4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
  linarith

private lemma log_seven_le : Real.log 7 ≤ 3 * Real.log 2 := by
  have h : Real.log 7 ≤ Real.log 8 := Real.log_le_log (by norm_num) (by norm_num)
  have e8 : Real.log 8 = 3 * Real.log 2 := by
    rw [show (8 : ℝ) = 2 ^ 3 by norm_num, Real.log_pow]; norm_num
  linarith

/-- **The numerical inequality.**  `7/2 − γ − 3 log 2 − log π + Σ_n Λ(n) (−χ₋₄(n)) (5n − 3)/n³
< 16/15`  (the left side is `≈ −0.3014 + 0.3476 = 0.0462`). -/
theorem maj_numeric :
    7 / 2 - Real.eulerMascheroniConstant - 3 * Real.log 2 - Real.log π
      + ∑' n : ℕ, ArithmeticFunction.vonMangoldt n
          * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)
      < 16 / 15 := by
  have hsplit : ∑' n : ℕ, ArithmeticFunction.vonMangoldt n
          * (-((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)) * ((5 * (n : ℝ) - 3) / (n : ℝ) ^ 3)
      = ∑ n ∈ Finset.range 11, majTerm n + ∑' n, majTerm (n + 11) :=
    (summable_majTerm.sum_add_tsum_nat_add 11).symm
  rw [hsplit]
  have hhead := majTerm_head
  have htail := majTerm_tail
  have hpot : majPot 7 = 5 * (Real.log 7 + 1) / 28 := by
    unfold majPot; norm_num
  rw [hpot] at htail
  have hγ := Real.one_half_lt_eulerMascheroniConstant
  have hlogπ : Real.log 3 ≤ Real.log π := Real.log_le_log (by norm_num) Real.pi_gt_three.le
  have h2lo := Real.log_two_gt_d9
  have h2hi := Real.log_two_lt_d9
  have h3lo := log_three_ge
  have h3hi := log_three_le
  have h5 := log_five_ge
  have h7 := log_seven_le
  linarith

end ChebLib
