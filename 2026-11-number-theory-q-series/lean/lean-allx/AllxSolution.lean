/-
AllxSolution.lean — the UNTRUSTED solution module: the statements of AllxChallenge.lean,
byte-identical, each PROVED by delegating to AllxLib (built on lean-chebyshev's ChebLib, on the
Zeta23 library — anthropics/formal-math, folder zeta23/ — and on Mathlib).

  chebyshev_race_cesaro_pos            →  Allx.cesaro_w4_ge_and_pos        (the headline)
  chebyshev_race_cesaro_eq_integral    →  Allx.cesaro_w4_eq_integral
  chebyshev_race_integral_pos          →  Allx.integral_race_pos
  chebyshev_race_theta1_pos            →  Allx.theta1_wlog_pos
  chebyshev_race_cesaro_lower_bound    →  Allx.cesaro_w4_ge_and_pos
  chebyshev_race_lead_inheritance      →  Allx.lead_inheritance_w4
  chebyshev_race_cone_lower_bound,
  chebyshev_race_cone_lead_inheritance →  Allx.cone_w4_lower_bound, Allx.cone_w4_lead_inheritance
  chebyshev_race_cesaro_pos_of_le      →  Allx.cesaro_w4_pos_of_le         (kernel-checked, no GRH)
  grh_chi4_psi1_le                     →  Allx.psi1chi_le
  chebyshev_race_sum_eq_counts, neg_chi4_natCast, grhFor_of_dirichlet_zeros_on_line  (conventions)

The definitions used in the statements are `chi4`, `GRHFor` (lean-chebyshev, ChebLib/Race/Defs.lean)
and `primeCountingMod4` (AllxLib/Defs.lean); all three are Mathlib-only and are inlined verbatim in
AllxChallenge.lean.  Nothing in this file is part of the trusted base.
-/
import AllxLib

open Finset MeasureTheory

/-- **The Cesàro mean of Chebyshev's race is positive at EVERY `x > 3`, under GRH for one
L-function.**  Assume that every zero of the Dirichlet L-function `L(s, χ₋₄)` in the strip
`0 < Re s < 1` has real part `1/2` (`GRHFor chi4`).  Then for every real `x > 3`

    `R(x) = Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) > 0`,

the sum over the primes `p ≤ x`; `−χ₋₄(p)` is `+1` for `p ≡ 3 (mod 4)`, `−1` for `p ≡ 1 (mod 4)` and
`0` for `p = 2` (`neg_chi4_natCast`).  By `chebyshev_race_cesaro_eq_integral`,
`R(x) = (1/x) ∫₀ˣ (π(t;4,3) − π(t;4,1)) dt`: Chebyshev's bias holds on average at every `x`, not
only for all sufficiently large `x`. -/
theorem chebyshev_race_cesaro_pos (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x) := by
  exact (Allx.cesaro_w4_ge_and_pos hGRH hx).2

/-- **The Cesàro mean is the mean of `π(t;4,3) − π(t;4,1)`** (no hypothesis).  For every `x > 0`,
`Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) = (1/x) ∫₀ˣ (π(t;4,3) − π(t;4,1)) dt`, where
`π(t;4,a) = primeCountingMod4 a t` is the number of primes `p ≤ t` with `p ≡ a (mod 4)`. -/
theorem chebyshev_race_cesaro_eq_integral (x : ℝ) (hx : 0 < x) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x)
      = (1 / x) * ∫ t in (0 : ℝ)..x,
          ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ)) := by
  exact Allx.cesaro_w4_eq_integral hx

/-- **The integral form.**  Under GRH for `L(s, χ₋₄)`:
`∫₀ˣ (π(t;4,3) − π(t;4,1)) dt > 0` for every real `x > 3`. -/
theorem chebyshev_race_integral_pos (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    0 < ∫ t in (0 : ℝ)..x, ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ)) := by
  exact Allx.integral_race_pos hGRH hx

/-- **The integrated θ-race is positive at every `t > 3`** (Theorem 1 of
`monotonicity/explicit/EXPLICIT.md`, positivity part).  Under GRH for `L(s, χ₋₄)`:
`Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p) > 0` for every real `t > 3`. -/
theorem chebyshev_race_theta1_pos (hGRH : GRHFor chi4) (t : ℝ) (ht : 3 < t) :
    0 < ∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p) := by
  exact Allx.theta1_wlog_pos hGRH ht

/-- **The lower bound `R(x) ≥ Θ₁(x)/(x log x)`** (Theorem 2 of `EXPLICIT.md`, first inequality).
Under GRH for `L(s, χ₋₄)`, for every real `x > 3`:
`Σ_{p ≤ x} (−χ₋₄(p)) log p · (x − p) / (x log x) ≤ Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x)`. -/
theorem chebyshev_race_cesaro_lower_bound (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (x - p)) / (x * Real.log x)
      ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x) := by
  exact (Allx.cesaro_w4_ge_and_pos hGRH hx).1

/-- **Inheritance of leads.**  Under GRH for `L(s, χ₋₄)`, for every real `x > 3`: if the class
`3 mod 4` leads or ties in the race weighted by `p`, i.e. `Σ_{p ≤ x} (−χ₋₄(p)) p ≥ 0`, then it leads
strictly in Chebyshev's race: `π(x;4,1) < π(x;4,3)`. -/
theorem chebyshev_race_lead_inheritance (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ)) :
    primeCountingMod4 1 x < primeCountingMod4 3 x := by
  have h1 := Allx.lead_inheritance_w4 hGRH hx h
  rw [Allx.sum_w4_eq_counts] at h1
  exact_mod_cast sub_pos.mp h1

/-- **The cone of weights** (Theorem 5 of `EXPLICIT.md`).  Under GRH for `L(s, χ₋₄)`, for all real
`m ≤ 0`, `m < m' ≤ 1` and every real `x > 3`, with `D_m(x) = Σ_{p ≤ x} (−χ₋₄(p)) p^m`:
`D_m(x) − x^{−(m'−m)} D_{m'}(x) ≥ (m'−m) x^{m−1} Θ₁(x)/log x`  (the right side is positive by
`chebyshev_race_theta1_pos`). -/
theorem chebyshev_race_cone_lower_bound (hGRH : GRHFor chi4) (m m' x : ℝ) (hm : m ≤ 0)
    (hmm' : m < m') (hm'1 : m' ≤ 1) (hx : 3 < x) :
    (m' - m) * x ^ (m - 1) * (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (x - p)) / Real.log x
      ≤ (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m)
        - x ^ (-(m' - m))
          * ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m' := by
  exact Allx.cone_w4_lower_bound hGRH hm hmm' hm'1 hx

/-- **Inheritance of leads in the whole cone.**  Under GRH for `L(s, χ₋₄)`, for all real `m ≤ 0`,
`m < m' ≤ 1` and every real `x > 3`: if `Σ_{p ≤ x} (−χ₋₄(p)) p^{m'} ≥ 0` then
`Σ_{p ≤ x} (−χ₋₄(p)) p^m > 0`. -/
theorem chebyshev_race_cone_lead_inheritance (hGRH : GRHFor chi4) (m m' x : ℝ) (hm : m ≤ 0)
    (hmm' : m < m') (hm'1 : m' ≤ 1) (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m') :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m := by
  exact Allx.cone_w4_lead_inheritance hGRH hm hmm' hm'1 hx h

/-- **The finite range, with no hypothesis.**  For every real `x` with `3 < x ≤ 10⁷`,
`Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) > 0`.  (Proved through `Θ₁ > 0` on `(3, 10⁷]`, which is checked by
kernel reduction, `decide +kernel`, on 1000 chunks of `10⁴` integers: no `native_decide`.) -/
theorem chebyshev_race_cesaro_pos_of_le (x : ℝ) (hx : 3 < x) (hN : x ≤ 10000000) :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x) := by
  exact Allx.cesaro_w4_pos_of_le hx hN

/-- **The explicit bound behind the analytic range.**  Under GRH for `L(s, χ₋₄)`, for every real
`t ≥ 2`:  `ψ₁(t, χ₋₄) = Σ_{n ≤ t} Λ(n) χ₋₄(n) (t − n) ≤ (27/100)(t+1)^{3/2} + 4t + log t + 4 +
(1/3) log(t+1)`.  (The constant `27/100` bounds `Σ_ρ m_ρ/|ρ(ρ+1)|`, whose value is about `0.154`;
it is obtained in closed form, with no numerical information on the zeros.) -/
theorem grh_chi4_psi1_le (hGRH : GRHFor chi4) (t : ℝ) (ht : 2 ≤ t) :
    ∑ n ∈ Ioc 0 ⌊t⌋₊, ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * (t - n)
      ≤ 27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 4 * t + Real.log t + 4 + 1 / 3 * Real.log (t + 1) := by
  exact Allx.psi1chi_le hGRH ht

/-- `Σ_{p ≤ x} (−χ₋₄(p)) = π(x;4,3) − π(x;4,1)`: the unweighted sum is Chebyshev's race. -/
theorem chebyshev_race_sum_eq_counts (x : ℝ) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ)
      = (primeCountingMod4 3 x : ℝ) - (primeCountingMod4 1 x : ℝ) := by
  exact Allx.sum_w4_eq_counts x

/-- The sign convention: `−χ₋₄(n)` is `+1` for `n ≡ 3`, `−1` for `n ≡ 1 (mod 4)`, `0` for even `n`
(`χ₋₄` is Mathlib's `ZMod.χ₄`). -/
theorem neg_chi4_natCast (n : ℕ) :
    -((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = if n % 4 = 3 then 1 else if n % 4 = 1 then -1 else 0 := by
  rw [Allx.chi4_nat]
  split_ifs <;> simp_all

/-- The hypothesis `GRHFor χ` is the instance at `χ` of the statement "all zeros of all Dirichlet
L-functions in the open critical strip lie on the line `Re s = 1/2`". -/
theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ := by
  exact fun s h1 h2 h3 => h q χ s h1 h2 h3
