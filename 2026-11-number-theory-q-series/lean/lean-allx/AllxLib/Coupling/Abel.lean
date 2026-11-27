/-
Copyright (c) 2026. Released under the repository's license.

# The exact Abel identity for weighted prime sums (Target A)

For an arbitrary real weight `w : ℕ → ℝ` on the primes (in the application
`w p = −χ(p)` for the nontrivial character `χ` mod 4 or mod 3, or `w ≡ 1`) and a real
exponent `m`, put

  `D w m x = ∑_{p ≤ x, p prime} w(p) p^m`.

This file proves, for all real `m, m'` and `x ≥ 2` (with `κ = m' − m`),

  `D w m x − x^(−κ) D w m' x = κ ∫_2^x t^(−κ−1) D w m' t dt`        (`Coupling.abel_identity`)

which is the identity `R_m(x) = x^{−κ} R_{m'}(x) + κ ∫_2^x t^{−κ−1} R_{m'}(t) dt` of
`GRH-ONLY.md`, Theorem 2.1(a), and the "exact identity" behind Theorem A of
`monotone-coupling.tex` / `CORE-NOTES.md`. It holds for actual primes, with no hypothesis.

Corollaries:
* `Coupling.lead_inheritance`: if `κ > 0` and the integral is positive, then
  `D w m' x ≥ 0 ⟹ D w m x > 0`.
* `Coupling.cesaro_eq_integral`: the Cesàro form
  `R w x := ∑_{p ≤ x} w(p) (1 − p/x) = (1/x) ∫_0^x D w 0 t dt` for `x > 0`.
-/
import Mathlib

open Real Finset MeasureTheory intervalIntegral

namespace Coupling

/-- `D w m x = ∑_{p ≤ x, p prime} w(p) p^m` (real exponent `m`, real cut-off `x`). -/
noncomputable def primeSum (w : ℕ → ℝ) (m x : ℝ) : ℝ :=
  ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, w p * (p : ℝ) ^ m

/-- The Cesàro (Riesz) mean `R w x = ∑_{p ≤ x, p prime} w(p) (1 − p/x)`. -/
noncomputable def cesaroSum (w : ℕ → ℝ) (x : ℝ) : ℝ :=
  ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, w p * (1 - p / x)

/-- The coefficient sequence `c k = [k prime] w(k) k^m`. -/
noncomputable def primeCoeff (w : ℕ → ℝ) (m : ℝ) (k : ℕ) : ℝ :=
  if k.Prime then w k * (k : ℝ) ^ m else 0

lemma primeSum_eq_sum_coeff (w : ℕ → ℝ) (m x : ℝ) :
    primeSum w m x = ∑ k ∈ Icc 0 ⌊x⌋₊, primeCoeff w m k := by
  rw [primeSum, sum_filter]; rfl

/-- **The exact Abel identity (Target A).** For every real weight `w`, all real `m, m'` and
every `x ≥ 2`,
`D w m x − x^(−(m'−m)) D w m' x = (m'−m) ∫_2^x t^(−(m'−m)−1) D w m' t dt`. -/
theorem abel_identity (w : ℕ → ℝ) (m m' : ℝ) {x : ℝ} (hx : 2 ≤ x) :
    primeSum w m x - x ^ (-(m' - m)) * primeSum w m' x =
      (m' - m) * ∫ t in (2 : ℝ)..x, t ^ (-(m' - m) - 1) * primeSum w m' t := by
  set κ := m' - m with hκ
  have hc0 : primeCoeff w m' 0 = 0 := by simp [primeCoeff, Nat.not_prime_zero]
  have hc1 : primeCoeff w m' 1 = 0 := by simp [primeCoeff, Nat.not_prime_one]
  have hdiff : ∀ t ∈ Set.Icc (2 : ℝ) x, DifferentiableAt ℝ (fun t : ℝ => t ^ (-κ)) t := by
    intro t ht
    exact (Real.hasDerivAt_rpow_const (Or.inl (by linarith [ht.1]))).differentiableAt
  have hderiv : ∀ t, deriv (fun t : ℝ => t ^ (-κ)) t = -κ * t ^ (-κ - 1) :=
    fun t => Real.deriv_rpow_const t (-κ)
  have hint : IntegrableOn (deriv fun t : ℝ => t ^ (-κ)) (Set.Icc 2 x) := by
    rw [show (deriv fun t : ℝ => t ^ (-κ)) = fun t => -κ * t ^ (-κ - 1) from funext hderiv]
    refine ContinuousOn.integrableOn_Icc ?_
    refine continuousOn_const.mul ?_
    exact fun t ht => (Real.continuousAt_rpow_const _ _
      (Or.inl (by linarith [ht.1] : t ≠ 0))).continuousWithinAt
  have key := sum_mul_eq_sub_integral_mul₁ (primeCoeff w m') hc0 hc1 x hdiff hint
  -- left side: termwise `k^(-κ) * [k prime] w k k^m' = [k prime] w k k^m`
  have hlhs : ∑ k ∈ Icc 0 ⌊x⌋₊, (k : ℝ) ^ (-κ) * primeCoeff w m' k = primeSum w m x := by
    rw [primeSum_eq_sum_coeff]
    refine sum_congr rfl fun k _ => ?_
    unfold primeCoeff
    split_ifs with hk
    · have hk0 : (0 : ℝ) < k := by exact_mod_cast hk.pos
      rw [mul_left_comm, ← Real.rpow_add hk0]
      congr 2; rw [hκ]; ring
    · simp
  rw [hlhs, ← primeSum_eq_sum_coeff] at key
  rw [key]
  have hI : ∫ t in Set.Ioc (2 : ℝ) x, deriv (fun t : ℝ => t ^ (-κ)) t *
      ∑ k ∈ Icc 0 ⌊t⌋₊, primeCoeff w m' k =
      -κ * ∫ t in (2 : ℝ)..x, t ^ (-κ - 1) * primeSum w m' t := by
    rw [intervalIntegral.integral_of_le hx, ← MeasureTheory.integral_const_mul]
    refine setIntegral_congr_fun measurableSet_Ioc fun t _ => ?_
    simp only [hderiv, primeSum_eq_sum_coeff]
    ring
  rw [hI]; ring

/-- **Lead inheritance.** If `m < m'`, `x ≥ 2` and
`∫_2^x t^(−(m'−m)−1) D w m' t dt > 0`, then `D w m' x ≥ 0` implies `D w m x > 0`. -/
theorem lead_inheritance (w : ℕ → ℝ) {m m' x : ℝ} (hmm' : m < m') (hx : 2 ≤ x)
    (hI : 0 < ∫ t in (2 : ℝ)..x, t ^ (-(m' - m) - 1) * primeSum w m' t)
    (hD' : 0 ≤ primeSum w m' x) : 0 < primeSum w m x := by
  have h := abel_identity w m m' hx
  have hx0 : 0 < x ^ (-(m' - m)) := Real.rpow_pos_of_pos (by linarith) _
  have : 0 < (m' - m) * ∫ t in (2 : ℝ)..x, t ^ (-(m' - m) - 1) * primeSum w m' t :=
    mul_pos (by linarith) hI
  nlinarith [mul_nonneg hx0.le hD']

/-- `∫_0^x D w 0 t dt = ∑_{p ≤ x} w(p) (x − p)` for `x ≥ 0` (Abel summation with `f(t) = t`). -/
theorem integral_primeSum_zero (w : ℕ → ℝ) {x : ℝ} (hx : 0 ≤ x) :
    ∫ t in (0 : ℝ)..x, primeSum w 0 t =
      ∑ p ∈ Icc 0 ⌊x⌋₊ with p.Prime, w p * (x - p) := by
  have hdiff : ∀ t ∈ Set.Icc (0 : ℝ) x, DifferentiableAt ℝ (fun t : ℝ => t) t :=
    fun t _ => differentiableAt_id
  have hint : IntegrableOn (deriv fun t : ℝ => t) (Set.Icc 0 x) := by
    simp only [deriv_id'']; exact integrableOn_const (by simp)
  have key := sum_mul_eq_sub_integral_mul (primeCoeff w 0) hx hdiff hint
  simp only [deriv_id'', one_mul] at key
  rw [intervalIntegral.integral_of_le hx]
  have h1 : ∫ t in Set.Ioc (0 : ℝ) x, primeSum w 0 t =
      ∫ t in Set.Ioc (0 : ℝ) x, ∑ k ∈ Icc 0 ⌊t⌋₊, primeCoeff w 0 k := by
    simp_rw [primeSum_eq_sum_coeff]
  rw [h1]
  have h2 : ∫ t in Set.Ioc (0 : ℝ) x, ∑ k ∈ Icc 0 ⌊t⌋₊, primeCoeff w 0 k =
      x * ∑ k ∈ Icc 0 ⌊x⌋₊, primeCoeff w 0 k -
        ∑ k ∈ Icc 0 ⌊x⌋₊, (k : ℝ) * primeCoeff w 0 k := by linarith
  rw [h2, mul_sum, ← sum_sub_distrib, sum_filter]
  refine sum_congr rfl fun k _ => ?_
  unfold primeCoeff
  split_ifs <;> simp; ring

/-- **The Cesàro form.** For `x > 0`,
`R w x = ∑_{p ≤ x} w(p)(1 − p/x) = (1/x) ∫_0^x D w 0 t dt`. -/
theorem cesaro_eq_integral (w : ℕ → ℝ) {x : ℝ} (hx : 0 < x) :
    cesaroSum w x = (1 / x) * ∫ t in (0 : ℝ)..x, primeSum w 0 t := by
  rw [integral_primeSum_zero w hx.le, cesaroSum, mul_sum]
  refine sum_congr rfl fun k _ => ?_
  field_simp

/-- The weight `s(p) = −χ₄(p)` of the race mod 4 (`+1` on `p ≡ 3`, `−1` on `p ≡ 1 mod 4`). -/
noncomputable def w4 (n : ℕ) : ℝ := -((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)

/-- Target A specialised to the race mod 4. -/
theorem abel_identity_chi4 (m m' : ℝ) {x : ℝ} (hx : 2 ≤ x) :
    primeSum w4 m x - x ^ (-(m' - m)) * primeSum w4 m' x =
      (m' - m) * ∫ t in (2 : ℝ)..x, t ^ (-(m' - m) - 1) * primeSum w4 m' t :=
  abel_identity w4 m m' hx

end Coupling
