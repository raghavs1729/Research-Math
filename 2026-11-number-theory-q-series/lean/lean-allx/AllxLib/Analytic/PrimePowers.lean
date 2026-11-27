/-
Copyright (c) 2026. Released under the repository's license.

# Prime powers: from `ψ₁(t, χ₋₄)` to `Θ₁(t)` (elementary part)

With `c(p) = −χ₋₄(p) log p` and `Θ₁(t) = ∑_{p ≤ t} c(p) (t − p)` we prove, for every `t ≥ 1`,

  `Θ₁(t) ≥ −ψ₁(t, χ₋₄) + (2/3) log 2 · t^{3/2} − (explicit lower-order terms)`

(`Allx.theta1_wlog_ge`). The mechanism: the squares of odd prime powers all contribute with
the sign `+1` to `Ψ(u, χ₋₄)`, which gives the main term `ψ(√u)`; the remaining prime powers
are controlled by the Costa-Pereira inequality and Chebyshev's bounds from Mathlib.
-/
import AllxLib.Coupling.Convexity

open Real Finset MeasureTheory
open ArithmeticFunction hiding log map_zero

noncomputable section

namespace Allx

open Chebyshev

/-- `χ₋₄(n)` as a real number. -/
def X4r (n : ℕ) : ℝ := ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)

/-- `ψ₁(t, χ₋₄) = Σ_{n ≤ t} Λ(n) χ₋₄(n) (t − n)`. -/
def psi1chi (t : ℝ) : ℝ := ∑ n ∈ Ioc 0 ⌊t⌋₊, vonMangoldt n * X4r n * (t - n)

/-! ### The character -/

lemma X4r_eq (n : ℕ) :
    X4r n = if n % 2 = 0 then 0 else if n % 4 = 1 then 1 else -1 := by
  unfold X4r
  rw [ZMod.χ₄_nat_eq_if_mod_four]
  split_ifs <;> simp

lemma neg_one_le_X4r (n : ℕ) : -1 ≤ X4r n := by
  rw [X4r_eq]; split_ifs <;> norm_num

lemma X4r_sq (m : ℕ) : X4r (m ^ 2) = if 2 ∣ m then 0 else 1 := by
  have h : X4r (m ^ 2) = X4r m ^ 2 := by
    unfold X4r
    rw [Nat.cast_pow, map_pow]
    push_cast
    rfl
  rw [h, X4r_eq]
  have h2 : (2 ∣ m) ↔ m % 2 = 0 := Nat.dvd_iff_mod_eq_zero
  split_ifs <;> first | omega | norm_num

/-! ### The even prime powers -/

/-- `E2(y) = Σ_{m ≤ y, 2 ∣ m} Λ(m)`. -/
def E2 (y : ℝ) : ℝ := ∑ m ∈ Ioc 0 ⌊y⌋₊ with 2 ∣ m, vonMangoldt m

lemma E2_mono : Monotone E2 := by
  intro x y hxy
  exact sum_le_sum_of_subset_of_nonneg
    (filter_subset_filter _ (Ioc_subset_Ioc le_rfl (Nat.floor_le_floor hxy)))
    (fun _ _ _ => vonMangoldt_nonneg)

lemma E2_nonneg (y : ℝ) : 0 ≤ E2 y := sum_nonneg fun _ _ => vonMangoldt_nonneg

lemma E2_le_log {y : ℝ} (hy : 1 ≤ y) : E2 y ≤ log y := by
  have hN : 0 < ⌊y⌋₊ := Nat.floor_pos.2 hy
  set N := ⌊y⌋₊ with hNdef
  set J := Nat.log 2 N with hJ
  have hpow : 2 ^ J ≤ N := Nat.pow_log_le_self 2 hN.ne'
  have h1 : E2 y ≤ ∑ d ∈ (2 ^ J).divisors, vonMangoldt d := by
    unfold E2
    rw [← sum_filter_ne_zero]
    refine sum_le_sum_of_subset_of_nonneg ?_ (fun _ _ _ => vonMangoldt_nonneg)
    intro m hm
    simp only [mem_filter, mem_Ioc] at hm
    obtain ⟨⟨⟨-, hmN⟩, h2m⟩, hΛ⟩ := hm
    rw [vonMangoldt_ne_zero_iff, isPrimePow_nat_iff] at hΛ
    obtain ⟨p, k, hp, -, rfl⟩ := hΛ
    have hp2 : p = 2 :=
      ((Nat.prime_dvd_prime_iff_eq Nat.prime_two hp).1
        (Nat.prime_two.dvd_of_dvd_pow h2m)).symm
    subst hp2
    have hk : k ≤ J := Nat.le_log_of_pow_le (by norm_num) hmN
    exact Nat.mem_divisors.2 ⟨pow_dvd_pow 2 hk, by positivity⟩
  rw [vonMangoldt_sum] at h1
  have h2 : log ((2 ^ J : ℕ) : ℝ) ≤ log y := by
    apply log_le_log (by positivity)
    calc ((2 ^ J : ℕ) : ℝ) ≤ (N : ℝ) := by exact_mod_cast hpow
      _ ≤ y := Nat.floor_le (by linarith)
  linarith

/-! ### Abel summation with `f(u) = u` -/

lemma integral_sum_Icc (a : ℕ → ℝ) {t : ℝ} (ht : 0 ≤ t) :
    ∫ u in (0 : ℝ)..t, ∑ n ∈ Icc 0 ⌊u⌋₊, a n = ∑ n ∈ Icc 0 ⌊t⌋₊, a n * (t - n) := by
  have hdiff : ∀ u ∈ Set.Icc (0 : ℝ) t, DifferentiableAt ℝ (fun u : ℝ => u) u :=
    fun u _ => differentiableAt_id
  have hint : IntegrableOn (deriv fun u : ℝ => u) (Set.Icc 0 t) := by
    simp only [deriv_id'']; exact integrableOn_const (by simp)
  have key := sum_mul_eq_sub_integral_mul a ht hdiff hint
  simp only [deriv_id'', one_mul] at key
  rw [intervalIntegral.integral_of_le ht]
  have h2 : ∫ u in Set.Ioc (0 : ℝ) t, ∑ k ∈ Icc 0 ⌊u⌋₊, a k =
      t * ∑ k ∈ Icc 0 ⌊t⌋₊, a k - ∑ k ∈ Icc 0 ⌊t⌋₊, (k : ℝ) * a k := by linarith
  rw [h2, mul_sum, ← sum_sub_distrib]
  refine sum_congr rfl fun k _ => ?_
  ring

lemma intervalIntegrable_sum_Icc (a : ℕ → ℝ) {t : ℝ} (ht : 0 ≤ t) :
    IntervalIntegrable (fun u : ℝ => ∑ n ∈ Icc 0 ⌊u⌋₊, a n) volume 0 t := by
  have h := integrableOn_mul_sum_Icc (m := 0) a (le_refl (0 : ℝ)) (b := t)
    (g := fun _ => (1 : ℝ)) (integrableOn_const (by simp))
  simp only [one_mul] at h
  rw [intervalIntegrable_iff_integrableOn_Icc_of_le ht]
  exact h

/-! ### The non-prime part of `Ψ(u, χ₋₄)` -/

/-- Coefficients of `Θ(u) + Ψ(u, χ₋₄)`: only non-primes contribute. -/
def gco (n : ℕ) : ℝ := if n.Prime then 0 else vonMangoldt n * X4r n

/-- `G(u) = Θ(u) + Ψ(u, χ₋₄)`. -/
def Gsum (u : ℝ) : ℝ := ∑ n ∈ Icc 0 ⌊u⌋₊, gco n

lemma sum_gco_mul (t : ℝ) :
    ∑ n ∈ Icc 0 ⌊t⌋₊, gco n * (t - n) =
      Coupling.Theta1 (fun p => Coupling.w4 p * log p) t + psi1chi t := by
  unfold Coupling.Theta1 psi1chi
  rw [sum_filter, ← add_sum_Ioc_eq_sum_Icc (Nat.zero_le _),
    ← add_sum_Ioc_eq_sum_Icc (Nat.zero_le _)]
  simp only [gco, Nat.not_prime_zero, if_false, ArithmeticFunction.map_zero, zero_mul, zero_add]
  rw [← sum_add_distrib]
  refine sum_congr rfl fun n _ => ?_
  split_ifs with hn
  · rw [vonMangoldt_apply_prime hn]
    simp only [Coupling.w4, X4r]
    ring
  · ring

/-- **Step (A).** -/
lemma Gsum_ge {u : ℝ} (hu : 0 ≤ u) :
    ψ (√u) - E2 (√u) - ψ (u ^ (3 : ℝ)⁻¹) - ψ (u ^ (5 : ℝ)⁻¹) ≤ Gsum u := by
  set N := ⌊u⌋₊ with hN
  set S := (Ioc 0 N).filter (fun n => ¬ n.Prime) with hS
  have hG : Gsum u = ∑ n ∈ S, vonMangoldt n * X4r n := by
    unfold Gsum
    rw [← add_sum_Ioc_eq_sum_Icc (Nat.zero_le _), hS, sum_filter]
    simp only [gco, Nat.not_prime_zero, if_false, ArithmeticFunction.map_zero, zero_mul, zero_add]
    refine sum_congr rfl fun n _ => ?_
    split_ifs <;> rfl
  have hsplit : ∑ n ∈ S, vonMangoldt n * X4r n =
      ∑ n ∈ S, vonMangoldt n * (X4r n + 1) - ∑ n ∈ S, vonMangoldt n := by
    rw [← sum_sub_distrib]
    exact sum_congr rfl fun n _ => by ring
  have hnp : ∑ n ∈ S, vonMangoldt n = ψ u - θ u := (psi_sub_theta_eq_sum_not_prime u).symm
  have hCP := psi_sub_theta_le_psi_add_psi_add_psi u
  have hsqrt : u ^ (2 : ℝ)⁻¹ = √u := by rw [Real.sqrt_eq_rpow]; norm_num
  rw [hsqrt] at hCP
  -- the squares
  have hsq : 2 * ψ (√u) - E2 (√u) ≤ ∑ n ∈ S, vonMangoldt n * (X4r n + 1) := by
    have himg : ∑ m ∈ Ioc 0 ⌊√u⌋₊, vonMangoldt (m ^ 2) * (X4r (m ^ 2) + 1) =
        ∑ n ∈ (Ioc 0 ⌊√u⌋₊).image (fun m => m ^ 2), vonMangoldt n * (X4r n + 1) := by
      rw [sum_image]
      intro a _ b _ hab
      exact Nat.pow_left_injective (by norm_num : (2 : ℕ) ≠ 0) hab
    have hsub : (Ioc 0 ⌊√u⌋₊).image (fun m => m ^ 2) ⊆ S := by
      intro n hn
      obtain ⟨m, hm, rfl⟩ := mem_image.1 hn
      rw [mem_Ioc] at hm
      rw [hS, mem_filter, mem_Ioc]
      refine ⟨⟨pow_pos hm.1 2, ?_⟩, Nat.Prime.not_prime_pow le_rfl⟩
      have h1 : (m : ℝ) ≤ √u := (Nat.le_floor_iff (Real.sqrt_nonneg u)).1 hm.2
      have h2 : (m : ℝ) ^ 2 ≤ u := by
        calc (m : ℝ) ^ 2 ≤ (√u) ^ 2 := by gcongr
          _ = u := Real.sq_sqrt hu
      exact Nat.le_floor (by exact_mod_cast h2)
    have hle := sum_le_sum_of_subset_of_nonneg hsub
      (f := fun n => vonMangoldt n * (X4r n + 1))
      (fun n _ _ => mul_nonneg vonMangoldt_nonneg (by linarith [neg_one_le_X4r n]))
    have heval : ∑ m ∈ Ioc 0 ⌊√u⌋₊, vonMangoldt (m ^ 2) * (X4r (m ^ 2) + 1) =
        2 * ψ (√u) - E2 (√u) := by
      unfold Chebyshev.psi E2
      rw [sum_filter, mul_sum, ← sum_sub_distrib]
      refine sum_congr rfl fun m _ => ?_
      rw [vonMangoldt_apply_pow (by norm_num), X4r_sq]
      split_ifs <;> ring
    linarith
  linarith

lemma integral_Gsum {t : ℝ} (ht : 0 ≤ t) :
    ∫ u in (0 : ℝ)..t, Gsum u =
      Coupling.Theta1 (fun p => Coupling.w4 p * log p) t + psi1chi t := by
  unfold Gsum
  rw [integral_sum_Icc gco ht, sum_gco_mul]

/-! ### Chebyshev bounds along `u ↦ u ^ a` -/

lemma psi_rpow_le {a u t : ℝ} (ha : 0 < a) (hu : 0 ≤ u) (hut : u ≤ t) (ht : 1 ≤ t) :
    ψ (u ^ a) ≤ log 4 * u ^ a + 2 * a * log t * u ^ (a / 2) := by
  have hlogt : 0 ≤ log t := log_nonneg ht
  have hlog4 : 0 ≤ log (4 : ℝ) := log_nonneg (by norm_num)
  rcases lt_or_ge u 1 with h1 | h1
  · have : u ^ a < 2 := by
      have : u ^ a < 1 := Real.rpow_lt_one hu h1 ha
      linarith
    rw [psi_eq_zero_of_lt_two this]
    have h2 := mul_nonneg hlog4 (Real.rpow_nonneg hu a)
    have h3 := mul_nonneg (mul_nonneg (mul_nonneg zero_le_two ha.le) hlogt)
      (Real.rpow_nonneg hu (a / 2))
    linarith
  · have hx : 1 ≤ u ^ a := Real.one_le_rpow h1 ha.le
    have h := psi_le hx
    have hs : √(u ^ a) = u ^ (a / 2) := by
      rw [Real.sqrt_eq_rpow, ← Real.rpow_mul hu]; congr 1; ring
    have hl : log (u ^ a) = a * log u := Real.log_rpow (by linarith) a
    rw [hs, hl] at h
    have hlu : log u ≤ log t := log_le_log (by linarith) hut
    have hpos : 0 ≤ u ^ (a / 2) := Real.rpow_nonneg hu _
    have : 2 * u ^ (a / 2) * (a * log u) ≤ 2 * a * log t * u ^ (a / 2) := by
      have := mul_le_mul_of_nonneg_left hlu (mul_nonneg hpos ha.le)
      nlinarith
    linarith

lemma integral_rpow_zero {r t : ℝ} (hr : -1 < r) :
    ∫ u in (0 : ℝ)..t, u ^ r = t ^ (r + 1) / (r + 1) := by
  rw [integral_rpow (Or.inl hr), Real.zero_rpow (by linarith), sub_zero]

lemma integral_psi_rpow_le {a t : ℝ} (ha : 0 < a) (ht : 1 ≤ t) :
    ∫ u in (0 : ℝ)..t, ψ (u ^ a) ≤
      log 4 * (t ^ (a + 1) / (a + 1)) + 2 * a * log t * (t ^ (a / 2 + 1) / (a / 2 + 1)) := by
  have ht0 : (0 : ℝ) ≤ t := by linarith
  have hmono : MonotoneOn (fun u : ℝ => ψ (u ^ a)) (Set.uIcc 0 t) := by
    intro x hx y _ hxy
    rw [Set.uIcc_of_le ht0] at hx
    exact psi_mono (Real.rpow_le_rpow hx.1 hxy ha.le)
  have hI1 : IntervalIntegrable (fun u : ℝ => u ^ a) volume 0 t :=
    intervalIntegral.intervalIntegrable_rpow' (by linarith)
  have hI2 : IntervalIntegrable (fun u : ℝ => u ^ (a / 2)) volume 0 t :=
    intervalIntegral.intervalIntegrable_rpow' (by linarith)
  have hle := intervalIntegral.integral_mono_on ht0 hmono.intervalIntegrable
    ((hI1.const_mul (log 4)).add (hI2.const_mul (2 * a * log t)))
    (fun u hu => psi_rpow_le ha hu.1 hu.2 ht)
  rw [intervalIntegral.integral_add (hI1.const_mul _) (hI2.const_mul _),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_rpow_zero (by linarith), integral_rpow_zero (by linarith)] at hle
  exact hle

/-! ### The main statement -/

/-- **Main statement of this file.** -/
theorem theta1_wlog_ge {t : ℝ} (ht : 1 ≤ t) :
    - psi1chi t
      + 2 / 3 * log 2 * t ^ (3 / 2 : ℝ)
      - t * (log 2 + log (√t + 2)) - t / 2 * log t
      - 3 / 4 * log 4 * t ^ (4 / 3 : ℝ) - 4 / 7 * t ^ (7 / 6 : ℝ) * log t
      - 5 / 6 * log 4 * t ^ (6 / 5 : ℝ) - 4 / 11 * t ^ (11 / 10 : ℝ) * log t
    ≤ Coupling.Theta1 (fun p => Coupling.w4 p * log p) t := by
  have ht0 : (0 : ℝ) ≤ t := by linarith
  have huIcc : Set.uIcc (0 : ℝ) t = Set.Icc 0 t := Set.uIcc_of_le ht0
  -- integrability
  have hm1 : MonotoneOn (fun u : ℝ => ψ (√u)) (Set.uIcc 0 t) :=
    fun x _ y _ hxy => psi_mono (Real.sqrt_le_sqrt hxy)
  have hm2 : MonotoneOn (fun u : ℝ => E2 (√u)) (Set.uIcc 0 t) :=
    fun x _ y _ hxy => E2_mono (Real.sqrt_le_sqrt hxy)
  have hm3 : MonotoneOn (fun u : ℝ => ψ (u ^ (3 : ℝ)⁻¹)) (Set.uIcc 0 t) := by
    intro x hx y _ hxy
    rw [huIcc] at hx
    exact psi_mono (Real.rpow_le_rpow hx.1 hxy (by norm_num))
  have hm5 : MonotoneOn (fun u : ℝ => ψ (u ^ (5 : ℝ)⁻¹)) (Set.uIcc 0 t) := by
    intro x hx y _ hxy
    rw [huIcc] at hx
    exact psi_mono (Real.rpow_le_rpow hx.1 hxy (by norm_num))
  have hi1 := hm1.intervalIntegrable (μ := volume)
  have hi2 := hm2.intervalIntegrable (μ := volume)
  have hi3 := hm3.intervalIntegrable (μ := volume)
  have hi5 := hm5.intervalIntegrable (μ := volume)
  have hiG : IntervalIntegrable Gsum volume 0 t := intervalIntegrable_sum_Icc gco ht0
  -- integrate (A)
  have hA := intervalIntegral.integral_mono_on ht0 (((hi1.sub hi2).sub hi3).sub hi5) hiG
    (fun u hu => Gsum_ge hu.1)
  rw [integral_Gsum ht0, intervalIntegral.integral_sub ((hi1.sub hi2).sub hi3) hi5,
    intervalIntegral.integral_sub (hi1.sub hi2) hi3, intervalIntegral.integral_sub hi1 hi2] at hA
  -- the main term
  have hsqrtI : IntervalIntegrable (fun u : ℝ => √u) volume 0 t :=
    (Real.continuous_sqrt.intervalIntegrable _ _)
  have h1 : 2 / 3 * log 2 * t ^ (3 / 2 : ℝ) - t * (log 2 + log (√t + 2)) ≤
      ∫ u in (0 : ℝ)..t, ψ (√u) := by
    have hlow : IntervalIntegrable
        (fun u : ℝ => log 2 * √u - (log 2 + log (√t + 2))) volume 0 t :=
      (hsqrtI.const_mul _).sub intervalIntegrable_const
    have hle := intervalIntegral.integral_mono_on ht0 hlow hi1 (fun u hu => by
      have h := psi_ge' (Real.sqrt_nonneg u)
      have hl : log (√u + 2) ≤ log (√t + 2) :=
        log_le_log (by positivity) (by linarith [Real.sqrt_le_sqrt hu.2])
      show log 2 * √u - (log 2 + log (√t + 2)) ≤ ψ (√u)
      linarith)
    have hsq : ∫ u in (0 : ℝ)..t, √u = t ^ (3 / 2 : ℝ) / (3 / 2) := by
      have : (fun u : ℝ => √u) = fun u : ℝ => u ^ (1 / 2 : ℝ) := by
        funext u; exact Real.sqrt_eq_rpow u
      rw [this, integral_rpow_zero (by norm_num)]
      norm_num
    rw [intervalIntegral.integral_sub (hsqrtI.const_mul _) intervalIntegrable_const,
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const, hsq] at hle
    simp only [sub_zero, smul_eq_mul] at hle
    linarith
  -- the even prime powers
  have h2 : ∫ u in (0 : ℝ)..t, E2 (√u) ≤ t / 2 * log t := by
    have hle := intervalIntegral.integral_mono_on ht0 hi2
      (intervalIntegrable_const (c := log t / 2)) (fun u hu => by
      have h := E2_mono (Real.sqrt_le_sqrt hu.2)
      have h' := E2_le_log (Real.one_le_sqrt.2 ht)
      rw [Real.log_sqrt ht0] at h'
      show E2 (√u) ≤ log t / 2
      linarith)
    rw [intervalIntegral.integral_const] at hle
    simp only [sub_zero, smul_eq_mul] at hle
    linarith
  -- cube roots and fifth roots
  have h3 := integral_psi_rpow_le (a := (3 : ℝ)⁻¹) (by norm_num) ht
  have h5 := integral_psi_rpow_le (a := (5 : ℝ)⁻¹) (by norm_num) ht
  have e43 : ((3 : ℝ)⁻¹ + 1) = 4 / 3 := by norm_num
  have e76 : ((3 : ℝ)⁻¹ / 2 + 1) = 7 / 6 := by norm_num
  have e65 : ((5 : ℝ)⁻¹ + 1) = 6 / 5 := by norm_num
  have e1110 : ((5 : ℝ)⁻¹ / 2 + 1) = 11 / 10 := by norm_num
  rw [e43, e76] at h3
  rw [e65, e1110] at h5
  linarith

end Allx
