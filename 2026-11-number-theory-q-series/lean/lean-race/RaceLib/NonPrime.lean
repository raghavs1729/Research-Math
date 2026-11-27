/-
RaceLib/NonPrime.lean — prime squares from the prime number theorem (Mathlib only).

Let `c(n) = Λ(n) χ₋₄(n)` for `n` not prime and `c(n) = 0` for `n` prime.  Then, ASSUMING the prime
number theorem in the form `θ(x)/x → 1`,   `Σ_{n ≤ x} c(n) = √x + o(√x)`:
the squares of odd primes contribute `θ(√x) − O(1)` and the higher powers `O(x^{1/3} log x)`.
-/
import Mathlib

open Filter Topology Finset

noncomputable section

namespace RaceLib

/-- `c(n) = Λ(n) χ₋₄(n)` for `n` not prime, and `0` for `n` prime. -/
def npCoeff (n : ℕ) : ℝ :=
  if n.Prime then 0 else ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)

theorem chi4_abs_le_one (n : ℕ) : |((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)| ≤ 1 := by
  rw [ZMod.χ₄_nat_eq_if_mod_four]
  split_ifs <;> simp

theorem chi4_sq_of_odd {n : ℕ} (hn : n % 2 = 1) : ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) ^ 2 = 1 := by
  rw [ZMod.χ₄_nat_eq_if_mod_four]
  split_ifs with h1 h2
  · omega
  · simp
  · simp

theorem chi4_natCast_two : ((ZMod.χ₄ ((2 : ℕ) : ZMod 4) : ℤ) : ℝ) = 0 := by
  have : ZMod.χ₄ ((2 : ℕ) : ZMod 4) = 0 := by decide
  rw [this]; simp

theorem npCoeff_zero : npCoeff 0 = 0 := by simp [npCoeff]

theorem npCoeff_prime_pow {p k : ℕ} (hp : p.Prime) (hk : k ≠ 0) :
    npCoeff (p ^ k)
      = if k = 1 then 0 else Real.log p * ((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) ^ k := by
  unfold npCoeff
  have hprime : (p ^ k).Prime ↔ k = 1 := by
    constructor
    · intro h; exact ((h.pow_eq_iff).mp rfl).2
    · rintro rfl; simpa using hp
  by_cases h1 : k = 1
  · rw [if_pos (hprime.mpr h1), if_pos h1]
  · rw [if_neg (fun h => h1 (hprime.mp h)), if_neg h1,
      ArithmeticFunction.vonMangoldt_apply_pow hk, ArithmeticFunction.vonMangoldt_apply_prime hp]
    simp only [Nat.cast_pow, map_pow, Int.cast_pow]

/-- The contribution of the `k`-th powers of primes to `Σ_{n ≤ x} c(n)`. -/
def npInner (x : ℝ) (k : ℕ) : ℝ :=
  ∑ p ∈ Ioc 0 ⌊x ^ ((1 : ℝ) / k)⌋₊ with p.Prime, npCoeff (p ^ k)

theorem sum_npCoeff_eq {x : ℝ} (hx : 0 ≤ x) :
    ∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n
      = ∑ k ∈ Icc 1 (⌊Real.log x / Real.log 2⌋₊ + 2), npInner x k := by
  rw [← Finset.add_sum_Ioc_eq_sum_Icc (Nat.zero_le _), npCoeff_zero, zero_add]
  unfold npInner
  rw [← Chebyshev.sum_PrimePow_eq_sum_sum' npCoeff hx (Nat.le_add_right _ 2)]
  symm
  apply Finset.sum_filter_of_ne
  intro n _ hn
  apply ArithmeticFunction.vonMangoldt_ne_zero_iff.mp
  intro h0
  apply hn
  unfold npCoeff
  split_ifs
  · rfl
  · rw [h0, zero_mul]

theorem abs_npInner_sub_le {x : ℝ} (hx : 1 ≤ x) {k : ℕ} (hk : 1 ≤ k) :
    |npInner x k - (if k = 2 then Chebyshev.theta (x ^ ((1 : ℝ) / (2 : ℕ))) else 0)|
      ≤ Real.log 2 + Chebyshev.theta (x ^ ((1 : ℝ) / 3)) := by
  have hlog2 : 0 ≤ Real.log 2 := Real.log_nonneg (by norm_num)
  have hθ3 : 0 ≤ Chebyshev.theta (x ^ ((1 : ℝ) / 3)) := Chebyshev.theta_nonneg _
  have hterm : ∀ p ∈ (Ioc 0 ⌊x ^ ((1 : ℝ) / k)⌋₊).filter Nat.Prime,
      npCoeff (p ^ k)
        = if k = 1 then 0 else Real.log p * ((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) ^ k := fun p hp =>
    npCoeff_prime_pow (Finset.mem_filter.mp hp).2 (by omega)
  unfold npInner
  rw [Finset.sum_congr rfl hterm]
  rcases Nat.lt_or_ge k 2 with h1 | h2
  · -- k = 1
    have hk1 : k = 1 := by omega
    subst hk1
    simp only [if_true, Finset.sum_const_zero, OfNat.one_ne_ofNat, if_false, sub_zero, abs_zero]
    linarith
  rcases h2.eq_or_lt with h2 | h3
  · -- k = 2
    subst h2
    simp only [OfNat.ofNat_ne_one, if_false, if_true]
    refine le_trans ?_ (le_add_of_nonneg_right hθ3)
    unfold Chebyshev.theta
    rw [← Finset.sum_sub_distrib]
    have hpt : ∀ p ∈ (Ioc 0 ⌊x ^ ((1 : ℝ) / (2 : ℕ))⌋₊).filter Nat.Prime,
        Real.log p * ((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) ^ 2 - Real.log p
          = if p = 2 then -Real.log 2 else 0 := by
      intro p hp
      have hpp : p.Prime := (Finset.mem_filter.mp hp).2
      rcases hpp.eq_two_or_odd with h | h
      · subst h
        rw [chi4_natCast_two, if_pos rfl]
        norm_num
      · have hne : p ≠ 2 := by omega
        rw [chi4_sq_of_odd h, if_neg hne]
        ring
    rw [Finset.sum_congr rfl hpt, Finset.sum_ite_eq']
    split_ifs
    · rw [abs_neg, abs_of_nonneg hlog2]
    · rw [abs_zero]; exact hlog2
  · -- k ≥ 3
    have hk1 : k ≠ 1 := by omega
    have hk2 : k ≠ 2 := by omega
    simp only [hk1, hk2, if_false, sub_zero]
    calc |∑ p ∈ (Ioc 0 ⌊x ^ ((1 : ℝ) / k)⌋₊).filter Nat.Prime,
            Real.log p * ((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) ^ k|
        ≤ ∑ p ∈ (Ioc 0 ⌊x ^ ((1 : ℝ) / k)⌋₊).filter Nat.Prime,
            |Real.log p * ((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) ^ k| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ p ∈ (Ioc 0 ⌊x ^ ((1 : ℝ) / k)⌋₊).filter Nat.Prime, Real.log p := by
          apply Finset.sum_le_sum
          intro p hp
          have hp1 : (1 : ℝ) ≤ p := by
            exact_mod_cast (Finset.mem_filter.mp hp).2.one_le
          have hlp : 0 ≤ Real.log p := Real.log_nonneg hp1
          rw [abs_mul, abs_pow, abs_of_nonneg hlp]
          calc Real.log p * |((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ)| ^ k ≤ Real.log p * 1 := by
                apply mul_le_mul_of_nonneg_left _ hlp
                exact pow_le_one₀ (abs_nonneg _) (chi4_abs_le_one p)
            _ = Real.log p := mul_one _
      _ = Chebyshev.theta (x ^ ((1 : ℝ) / k)) := rfl
      _ ≤ Chebyshev.theta (x ^ ((1 : ℝ) / 3)) := by
          apply Chebyshev.theta_mono
          apply Real.rpow_le_rpow_of_exponent_le hx
          have : (3 : ℝ) ≤ k := by exact_mod_cast h3
          rw [div_le_div_iff₀ (by linarith) (by norm_num)]
          linarith
      _ ≤ Real.log 2 + Chebyshev.theta (x ^ ((1 : ℝ) / 3)) := by linarith

theorem abs_sum_npCoeff_sub_theta_le {x : ℝ} (hx : 1 ≤ x) :
    |∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n - Chebyshev.theta (Real.sqrt x)|
      ≤ (Real.log x / Real.log 2 + 2) * (Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3)) := by
  have hx0 : 0 ≤ x := by linarith
  set N := ⌊Real.log x / Real.log 2⌋₊ + 2 with hN
  have hsqrt : x ^ ((1 : ℝ) / (2 : ℕ)) = Real.sqrt x := by
    rw [Real.sqrt_eq_rpow]; norm_num
  have h2mem : 2 ∈ Icc 1 N := Finset.mem_Icc.mpr ⟨by norm_num, by omega⟩
  have hsum : ∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n - Chebyshev.theta (Real.sqrt x)
      = ∑ k ∈ Icc 1 N,
          (npInner x k - (if k = 2 then Chebyshev.theta (x ^ ((1 : ℝ) / (2 : ℕ))) else 0)) := by
    rw [Finset.sum_sub_distrib, Finset.sum_ite_eq' (Icc 1 N) 2, if_pos h2mem, hsqrt,
      sum_npCoeff_eq hx0]
  have hθ : Chebyshev.theta (x ^ ((1 : ℝ) / 3)) ≤ Real.log 4 * x ^ ((1 : ℝ) / 3) :=
    Chebyshev.theta_le_log4_mul_x (Real.rpow_nonneg hx0 _)
  have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hlogx : 0 ≤ Real.log x := Real.log_nonneg hx
  have hNle : (N : ℝ) ≤ Real.log x / Real.log 2 + 2 := by
    rw [hN]; push_cast
    have := Nat.floor_le (div_nonneg hlogx hlog2.le)
    linarith
  have hB0 : 0 ≤ Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3) := by
    have : 0 ≤ Real.log 4 * x ^ ((1 : ℝ) / 3) :=
      mul_nonneg (Real.log_nonneg (by norm_num)) (Real.rpow_nonneg hx0 _)
    linarith
  rw [hsum]
  calc |∑ k ∈ Icc 1 N,
          (npInner x k - (if k = 2 then Chebyshev.theta (x ^ ((1 : ℝ) / (2 : ℕ))) else 0))|
      ≤ ∑ k ∈ Icc 1 N,
          |npInner x k - (if k = 2 then Chebyshev.theta (x ^ ((1 : ℝ) / (2 : ℕ))) else 0)| :=
        Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _k ∈ Icc 1 N, (Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3)) := by
        apply Finset.sum_le_sum
        intro k hk
        exact (abs_npInner_sub_le hx (Finset.mem_Icc.mp hk).1).trans (by linarith)
    _ = (N : ℝ) * (Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3)) := by
        rw [Finset.sum_const, Nat.card_Icc, nsmul_eq_mul]
        simp
    _ ≤ (Real.log x / Real.log 2 + 2) * (Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3)) :=
        mul_le_mul_of_nonneg_right hNle hB0

theorem abs_sum_npCoeff_sub_theta_div_le {x : ℝ} (hx : 1 ≤ x) :
    |(∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n - Chebyshev.theta (Real.sqrt x)) / Real.sqrt x|
      ≤ 78 * x ^ (-(1 / 12 : ℝ)) := by
  have hx0 : 0 < x := by linarith
  have hs : 0 < Real.sqrt x := Real.sqrt_pos.mpr hx0
  have h12 : (1 : ℝ) ≤ x ^ ((1 : ℝ) / 12) := Real.one_le_rpow hx (by norm_num)
  have h3 : (1 : ℝ) ≤ x ^ ((1 : ℝ) / 3) := Real.one_le_rpow hx (by norm_num)
  have hlogx : Real.log x ≤ 12 * x ^ ((1 : ℝ) / 12) := by
    have := Real.log_le_rpow_div hx0.le (show (0 : ℝ) < 1 / 12 by norm_num)
    calc Real.log x ≤ x ^ ((1 : ℝ) / 12) / (1 / 12) := this
      _ = 12 * x ^ ((1 : ℝ) / 12) := by ring
  have hlog2 : (1 / 2 : ℝ) < Real.log 2 := by
    have := Real.log_two_gt_d9; linarith
  have hlog2' : Real.log 2 < 1 := by
    have := Real.log_two_lt_d9; linarith
  have hlog4 : Real.log 4 < 2 := by
    have h : Real.log 4 = 2 * Real.log 2 := by
      rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]; norm_num
    linarith
  have hlogx0 : 0 ≤ Real.log x := Real.log_nonneg hx
  have hA : Real.log x / Real.log 2 + 2 ≤ 26 * x ^ ((1 : ℝ) / 12) := by
    have h1 : Real.log x / Real.log 2 ≤ Real.log x / (1 / 2) :=
      div_le_div_of_nonneg_left hlogx0 (by norm_num) hlog2.le
    have h2 : Real.log x / (1 / 2) = 2 * Real.log x := by ring
    linarith
  have hB : Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3) ≤ 3 * x ^ ((1 : ℝ) / 3) := by
    nlinarith
  have hB0 : 0 ≤ Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3) := by
    have : 0 ≤ Real.log 4 * x ^ ((1 : ℝ) / 3) :=
      mul_nonneg (Real.log_nonneg (by norm_num)) (by linarith)
    linarith
  rw [abs_div, abs_of_pos hs, div_le_iff₀ hs]
  have hpow : x ^ ((1 : ℝ) / 12) * x ^ ((1 : ℝ) / 3) = x ^ (-(1 / 12 : ℝ)) * Real.sqrt x := by
    rw [Real.sqrt_eq_rpow, ← Real.rpow_add hx0, ← Real.rpow_add hx0]
    norm_num
  calc |∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n - Chebyshev.theta (Real.sqrt x)|
      ≤ (Real.log x / Real.log 2 + 2) * (Real.log 2 + Real.log 4 * x ^ ((1 : ℝ) / 3)) :=
        abs_sum_npCoeff_sub_theta_le hx
    _ ≤ (26 * x ^ ((1 : ℝ) / 12)) * (3 * x ^ ((1 : ℝ) / 3)) :=
        mul_le_mul hA hB hB0 (by linarith)
    _ = 78 * (x ^ ((1 : ℝ) / 12) * x ^ ((1 : ℝ) / 3)) := by ring
    _ = 78 * x ^ (-(1 / 12 : ℝ)) * Real.sqrt x := by rw [hpow]; ring

/-- **Prime squares from PNT.**  If `θ(x)/x → 1` then `(Σ_{n ≤ x} c(n))/√x → 1`. -/
theorem tendsto_sum_npCoeff_div_sqrt
    (hPNT : Tendsto (fun x : ℝ => Chebyshev.theta x / x) atTop (𝓝 1)) :
    Tendsto (fun x : ℝ => (∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n) / Real.sqrt x) atTop (𝓝 1) := by
  have h1 : Tendsto (fun x : ℝ =>
      (∑ n ∈ Icc 0 ⌊x⌋₊, npCoeff n - Chebyshev.theta (Real.sqrt x)) / Real.sqrt x)
      atTop (𝓝 0) := by
    apply squeeze_zero_norm' (a := fun x : ℝ => 78 * x ^ (-(1 / 12 : ℝ)))
    · filter_upwards [eventually_ge_atTop (1 : ℝ)] with x hx
      rw [Real.norm_eq_abs]
      exact abs_sum_npCoeff_sub_theta_div_le hx
    · simpa using (tendsto_rpow_neg_atTop (show (0 : ℝ) < 1 / 12 by norm_num)).const_mul (78 : ℝ)
  have h2 : Tendsto (fun x : ℝ => Chebyshev.theta (Real.sqrt x) / Real.sqrt x) atTop (𝓝 1) :=
    hPNT.comp Real.tendsto_sqrt_atTop
  have h3 := h1.add h2
  rw [zero_add] at h3
  refine h3.congr fun x => ?_
  rw [← add_div, sub_add_cancel]

end RaceLib
