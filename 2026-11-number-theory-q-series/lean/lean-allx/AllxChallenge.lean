/-
AllxChallenge.lean — the TRUSTED challenge module: WHAT IS CLAIMED.   (GENERATED FILE: see
scripts/gen_challenge.py; the definitions `GRHFor`, `chi4` are inlined verbatim from
lean-chebyshev's ChebLib/Race/Defs.lean, `primeCountingMod4` from AllxLib/Defs.lean, and the
statements are copied byte for byte from AllxSolution.lean, with every proof replaced by `sorry`.)

THE CESÀRO MEAN OF CHEBYSHEV'S RACE IS POSITIVE AT EVERY x > 3, UNDER GRH FOR ONE L-FUNCTION.

  * `chebyshev_race_cesaro_pos` — the headline: if every zero of L(s,χ₋₄) in 0 < Re s < 1 has real
    part 1/2 (`GRHFor chi4`), then for EVERY real x > 3
        Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) > 0        (−χ₋₄(p) = +1, −1, 0 for p ≡ 3, 1 mod 4, p = 2);
  * `chebyshev_race_cesaro_eq_integral` (no hypothesis) — that sum is (1/x) ∫₀ˣ (π(t;4,3) − π(t;4,1)) dt;
    `chebyshev_race_integral_pos` — hence the integral is positive for every x > 3 under GRH;
  * `chebyshev_race_theta1_pos` — Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p) > 0 for every t > 3;
    `chebyshev_race_cesaro_lower_bound` — R(x) ≥ Θ₁(x)/(x log x);
  * `chebyshev_race_lead_inheritance` — for every x > 3: Σ_{p ≤ x} (−χ₋₄(p)) p ≥ 0 ⟹ π(x;4,1) < π(x;4,3);
  * `chebyshev_race_cone_lower_bound`, `chebyshev_race_cone_lead_inheritance` — the cone of weights
    m ≤ 0, m < m' ≤ 1: D_m(x) − x^{−(m'−m)} D_{m'}(x) ≥ (m'−m) x^{m−1} Θ₁(x)/log x, hence
    D_{m'}(x) ≥ 0 ⟹ D_m(x) > 0, for every x > 3  (D_m(x) = Σ_{p ≤ x} (−χ₋₄(p)) p^m);
  * `chebyshev_race_cesaro_pos_of_le` (no hypothesis) — the headline inequality for 3 < x ≤ 10⁷
    (kernel-checked);
  * `grh_chi4_psi1_le` — ψ₁(t,χ₋₄) ≤ (27/100)(t+1)^{3/2} + 4t + log t + 4 + (1/3) log(t+1), t ≥ 2;
  * `chebyshev_race_sum_eq_counts`, `neg_chi4_natCast`, `grhFor_of_dirichlet_zeros_on_line` — the
    conventions: Σ_{p ≤ x} (−χ₋₄(p)) = π(x;4,3) − π(x;4,1); the sign of −χ₋₄; the shape of GRH.

This module imports ONLY Mathlib.  Every hypothesis is an explicit binder: there is no `axiom`.
AllxSolution.lean (untrusted) proves exactly these statements from AllxLib/, lean-chebyshev's
ChebLib/ and the Zeta23 library.  The `sorry`s below are deliberate (this is the challenge side).
-/
import Mathlib

/-! ## Definitions (inlined verbatim) -/

/-- The Riemann hypothesis for the single Dirichlet L-function `L(s,χ)`: every zero in the open
critical strip `0 < Re s < 1` has real part `1/2`.  This is a hypothesis of the theorems, passed as
an explicit argument; it is never an axiom.  (Its binder structure is that of the statement
`∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ), 0 < s.re → s.re < 1 →
DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2`, instantiated at one character.) -/
def GRHFor {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → χ.LFunction s = 0 → s.re = 1 / 2

/-- The nontrivial Dirichlet character modulo 4 (`χ₋₄(n) = 1, −1, 0` for `n ≡ 1, 3 mod 4` and `n`
even), as a `ℂ`-valued Dirichlet character, from Mathlib's `ZMod.χ₄ : MulChar (ZMod 4) ℤ`. -/
def chi4 : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)

/-- `π(t; 4, a)`: the number of primes `p ≤ t` with `p ≡ a (mod 4)` (for a real cut-off `t`). -/
noncomputable def primeCountingMod4 (a : ℕ) (t : ℝ) : ℕ :=
  ((Finset.Icc 0 ⌊t⌋₊).filter (fun p => p.Prime ∧ p % 4 = a)).card

open Finset MeasureTheory

/-! ## Statements -/

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
  sorry

/-- **The Cesàro mean is the mean of `π(t;4,3) − π(t;4,1)`** (no hypothesis).  For every `x > 0`,
`Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) = (1/x) ∫₀ˣ (π(t;4,3) − π(t;4,1)) dt`, where
`π(t;4,a) = primeCountingMod4 a t` is the number of primes `p ≤ t` with `p ≡ a (mod 4)`. -/
theorem chebyshev_race_cesaro_eq_integral (x : ℝ) (hx : 0 < x) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x)
      = (1 / x) * ∫ t in (0 : ℝ)..x,
          ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ)) := by
  sorry

/-- **The integral form.**  Under GRH for `L(s, χ₋₄)`:
`∫₀ˣ (π(t;4,3) − π(t;4,1)) dt > 0` for every real `x > 3`. -/
theorem chebyshev_race_integral_pos (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    0 < ∫ t in (0 : ℝ)..x, ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ)) := by
  sorry

/-- **The integrated θ-race is positive at every `t > 3`** (Theorem 1 of
`monotonicity/explicit/EXPLICIT.md`, positivity part).  Under GRH for `L(s, χ₋₄)`:
`Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p) > 0` for every real `t > 3`. -/
theorem chebyshev_race_theta1_pos (hGRH : GRHFor chi4) (t : ℝ) (ht : 3 < t) :
    0 < ∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p) := by
  sorry

/-- **The lower bound `R(x) ≥ Θ₁(x)/(x log x)`** (Theorem 2 of `EXPLICIT.md`, first inequality).
Under GRH for `L(s, χ₋₄)`, for every real `x > 3`:
`Σ_{p ≤ x} (−χ₋₄(p)) log p · (x − p) / (x log x) ≤ Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x)`. -/
theorem chebyshev_race_cesaro_lower_bound (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (x - p)) / (x * Real.log x)
      ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x) := by
  sorry

/-- **Inheritance of leads.**  Under GRH for `L(s, χ₋₄)`, for every real `x > 3`: if the class
`3 mod 4` leads or ties in the race weighted by `p`, i.e. `Σ_{p ≤ x} (−χ₋₄(p)) p ≥ 0`, then it leads
strictly in Chebyshev's race: `π(x;4,1) < π(x;4,3)`. -/
theorem chebyshev_race_lead_inheritance (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ)) :
    primeCountingMod4 1 x < primeCountingMod4 3 x := by
  sorry

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
  sorry

/-- **Inheritance of leads in the whole cone.**  Under GRH for `L(s, χ₋₄)`, for all real `m ≤ 0`,
`m < m' ≤ 1` and every real `x > 3`: if `Σ_{p ≤ x} (−χ₋₄(p)) p^{m'} ≥ 0` then
`Σ_{p ≤ x} (−χ₋₄(p)) p^m > 0`. -/
theorem chebyshev_race_cone_lead_inheritance (hGRH : GRHFor chi4) (m m' x : ℝ) (hm : m ≤ 0)
    (hmm' : m < m') (hm'1 : m' ≤ 1) (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m') :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m := by
  sorry

/-- **The finite range, with no hypothesis.**  For every real `x` with `3 < x ≤ 10⁷`,
`Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) > 0`.  (Proved through `Θ₁ > 0` on `(3, 10⁷]`, which is checked by
kernel reduction, `decide +kernel`, on 1000 chunks of `10⁴` integers: no `native_decide`.) -/
theorem chebyshev_race_cesaro_pos_of_le (x : ℝ) (hx : 3 < x) (hN : x ≤ 10000000) :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x) := by
  sorry

/-- **The explicit bound behind the analytic range.**  Under GRH for `L(s, χ₋₄)`, for every real
`t ≥ 2`:  `ψ₁(t, χ₋₄) = Σ_{n ≤ t} Λ(n) χ₋₄(n) (t − n) ≤ (27/100)(t+1)^{3/2} + 4t + log t + 4 +
(1/3) log(t+1)`.  (The constant `27/100` bounds `Σ_ρ m_ρ/|ρ(ρ+1)|`, whose value is about `0.154`;
it is obtained in closed form, with no numerical information on the zeros.) -/
theorem grh_chi4_psi1_le (hGRH : GRHFor chi4) (t : ℝ) (ht : 2 ≤ t) :
    ∑ n ∈ Ioc 0 ⌊t⌋₊, ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * (t - n)
      ≤ 27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 4 * t + Real.log t + 4 + 1 / 3 * Real.log (t + 1) := by
  sorry

/-- `Σ_{p ≤ x} (−χ₋₄(p)) = π(x;4,3) − π(x;4,1)`: the unweighted sum is Chebyshev's race. -/
theorem chebyshev_race_sum_eq_counts (x : ℝ) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ)
      = (primeCountingMod4 3 x : ℝ) - (primeCountingMod4 1 x : ℝ) := by
  sorry

/-- The sign convention: `−χ₋₄(n)` is `+1` for `n ≡ 3`, `−1` for `n ≡ 1 (mod 4)`, `0` for even `n`
(`χ₋₄` is Mathlib's `ZMod.χ₄`). -/
theorem neg_chi4_natCast (n : ℕ) :
    -((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = if n % 4 = 3 then 1 else if n % 4 = 1 then -1 else 0 := by
  sorry

/-- The hypothesis `GRHFor χ` is the instance at `χ` of the statement "all zeros of all Dirichlet
L-functions in the open critical strip lie on the line `Re s = 1/2`". -/
theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ := by
  sorry
