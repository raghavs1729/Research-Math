/-
ChebAllxChallenge.lean — the TRUSTED challenge module: WHAT IS CLAIMED.   (GENERATED FILE: see
scripts/gen_challenge.py; the definitions `GRHFor`, `chi4` are inlined verbatim from
lean-chebyshev's ChebLib/Race/Defs.lean, `chebyshevSum`, `chebyshevLogSum` from
ChebLib/Cheb/Defs.lean, and the statements are copied byte for byte from ChebAllxSolution.lean,
with every proof replaced by `sorry`.)

CHEBYSHEV'S 1853 ASSERTION AT EVERY x > 0, AND FOR EVERY CONVEX WEIGHT, UNDER GRH FOR ONE
L-FUNCTION.

  * `grh_chebyshev_sum_neg` — the headline: if every zero of L(s,χ₋₄) in 0 < Re s < 1 has real part
    1/2 (`GRHFor chi4`), then for EVERY real x > 0
        Σ_p χ₋₄(p) e^{−p/x} < 0        (χ₋₄(p) = (−1)^{(p−1)/2} for odd p, χ₋₄(2) = 0);
    `grh_chebyshev_log_sum_neg` — the same for Σ_p χ₋₄(p) log p · e^{−p/x};
  * `chebyshev_sum_neg_of_le` (NO hypothesis) — the same inequality for every 0 < x ≤ 250000,
    unconditionally;
  * `grh_chi4_convex_weight_sum_neg`, `grh_chi4_convexOn_weight_sum_neg`,
    `grh_chi4_convex_log_weight_sum_neg` — every convex weight: Σ_p χ₋₄(p) w(p) < 0 for every w convex
    on [3, ∞) (as a sequence or as a function) with Σ_p w(p) < ∞ and w(3) > 0;
  * `grh_fujii_sum_neg` — Fujii's weights: Σ_p χ₋₄(p) e^{−(p/x)^α} < 0 for every x > 0, 0 < α ≤ 1;
    `grh_chi4_prime_dirichlet_sum_neg` — Σ_p χ₋₄(p) p^{−σ} < 0 for every σ > 1;
  * `chi4_log_sum_eq_theta1_series` (no hypothesis) — the identity behind these:
        −Σ_p χ₋₄(p) log p · f(p) = Σ_{k ≥ 0} Θ₁(k+4) (f(k+3) − 2 f(k+4) + f(k+5)),
    Θ₁(m) = Σ_{p ≤ m} (−χ₋₄(p)) log p · (m − p);
  * `chebyshevSum_eq_neg_integral`, `chebyshevSum_integrand_integrableOn` (no hypothesis) — the
    integral identity of the note:  Σ_p χ₋₄(p) e^{−p/x} = −∫_2^∞ Θ₁(t) F''(t) dt,  F = e^{−t/x}/log t;
    `chebyshev_kernel_hasDerivAt`, `chebyshev_kernel_second_deriv_pos` — F'' in closed form, F'' > 0;
  * `chebyshevSum_eq_integral_log_sums` (no hypothesis), `grh_chi4_log_sum_rpow_exp_neg` — from the sum
    without logarithms to the sums with logarithms:  Σ_p χ₋₄(p) e^{−p/x} =
    ∫_0^∞ (Σ_p χ₋₄(p) log p · p^{−u} e^{−p/x}) du,  each integrand negative under GRH;
  * `second_order_summation_by_parts_infinite`, `second_order_partial_summation_infinite` — the two
    general summation-by-parts identities over the infinite range (discrete and integral form);
  * `chebyshevSum_summable`, `chebyshevLogSum_summable`, `chi4_natCast_two`, `chi4_natCast_of_odd`,
    `fujii_summable`, `grhFor_of_dirichlet_zeros_on_line` — the conventions: convergence of the series, the values of
    χ₋₄, the shape of GRH.

This module imports ONLY Mathlib.  Every hypothesis is an explicit binder: there is no `axiom`.
ChebAllxSolution.lean (untrusted) proves exactly these statements from ChebAllxLib/, lean-allx's
AllxLib/, lean-chebyshev's ChebLib/ and the Zeta23 library.  The `sorry`s below are deliberate
(this is the challenge side).
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

noncomputable section

/-- **Chebyshev's sum**  `Σ_p χ₋₄(p) e^{−p/x}`  over all primes `p` (`χ₋₄(2) = 0`,
`χ₋₄(p) = (−1)^{(p−1)/2}` for odd `p`; here `χ₋₄` is Mathlib's `ZMod.χ₄`). -/
def chebyshevSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)

/-- The logarithmically weighted sum  `Σ_p χ₋₄(p) (log p) e^{−p/x}`  over all primes `p`. -/
def chebyshevLogSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
    * Real.exp (-((p : ℕ) : ℝ) / x)

end

open Finset MeasureTheory

/-! ## Statements -/

/-- **Chebyshev's 1853 assertion, pointwise: the Abel sum is negative at EVERY `x > 0`, under GRH
for one L-function.**  Assume that every zero of the Dirichlet L-function `L(s, χ₋₄)` in the strip
`0 < Re s < 1` has real part `1/2` (`GRHFor chi4`).  Then for every real `x > 0`

    `C(x) = Σ_p χ₋₄(p) e^{−p/x} < 0`,

the sum over all primes `p` (`chebyshevSum x`), with `χ₋₄(2) = 0` and `χ₋₄(p) = (−1)^{(p−1)/2}` for
odd `p` (`chi4_natCast_two`, `chi4_natCast_of_odd`); the series converges absolutely
(`chebyshevSum_summable`).  The classical theorem (Hardy–Littlewood, Landau; `grh_chebyshev_bias`
in `monotonicity/lean-chebyshev`) is `C(x) → −∞` as `x → ∞`; this is the statement for every
single `x`. -/
theorem grh_chebyshev_sum_neg (hGRH : GRHFor chi4) (x : ℝ) (hx : 0 < x) :
    chebyshevSum x < 0 := by
  sorry

/-- **The logarithmically weighted sum is negative at EVERY `x > 0`.**  Under GRH for `L(s, χ₋₄)`:
`Σ_p χ₋₄(p) log p · e^{−p/x} < 0` for every real `x > 0` (`chebyshevLogSum x`; the series converges
absolutely, `chebyshevLogSum_summable`). -/
theorem grh_chebyshev_log_sum_neg (hGRH : GRHFor chi4) (x : ℝ) (hx : 0 < x) :
    chebyshevLogSum x < 0 := by
  sorry

/-- **A finite range with NO hypothesis.**  For every real `0 < x ≤ 250000`:
`C(x) = Σ_p χ₋₄(p) e^{−p/x} < 0`, unconditionally.  (The proof uses only the kernel-checked fact
`Θ₁(t) ≥ 0` for `2 ≤ t ≤ 10⁷` of `monotonicity/lean-allx`: the second-order summation by parts is
applied to the weight minus its secant at `10⁷`, cut off there, and the remainder is bounded
trivially by `((10⁷+1)² + x) e^{−10⁷/x}`.  Compare Theorem 5 of
`monotonicity/chebyshev1853/CHEBYSHEV1853.md`, which certifies `0 < x ≤ 69 905` numerically.) -/
theorem chebyshev_sum_neg_of_le (x : ℝ) (hx : 0 < x) (hN : x ≤ 250000) :
    chebyshevSum x < 0 := by
  sorry

/-- **Every convex weight.**  Under GRH for `L(s, χ₋₄)`: let `w : ℕ → ℝ` be convex on `n ≥ 3`
(`2 w(n+1) ≤ w(n) + w(n+2)`), with `Σ_p w(p)` convergent (`Summable`, i.e. absolutely) and
`w(3) > 0`.  Then

    `Σ_p χ₋₄(p) w(p) < 0`.

No smoothness, monotonicity or sign condition is assumed: a convex sequence which is summable over
the primes is automatically nonincreasing and nonnegative on `n ≥ 3`.  (If `w(3) = 0` then `w`
vanishes on `n ≥ 3` and the sum is `0`.)  The values `w(0), w(1), w(2)` play no role
(`χ₋₄(2) = 0`). -/
theorem grh_chi4_convex_weight_sum_neg (hGRH : GRHFor chi4) (w : ℕ → ℝ)
    (hconv : ∀ n, 3 ≤ n → 2 * w (n + 1) ≤ w n + w (n + 2))
    (hsum : Summable fun p : Nat.Primes => w p) (hw3 : 0 < w 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w p < 0 := by
  sorry

/-- **Every convex weight, real-variable form.**  Under GRH for `L(s, χ₋₄)`: let `w : ℝ → ℝ` be
convex on `[3, ∞)` with `Σ_p w(p)` convergent and `w(3) > 0`.  Then `Σ_p χ₋₄(p) w(p) < 0`.
(Theorem 1 of `monotonicity/chebyshev1853/CHEBYSHEV1853.md`, without its hypotheses `w ∈ C²`,
`w ≥ 0`, `w` nonincreasing, `t·w(t) → 0`.) -/
theorem grh_chi4_convexOn_weight_sum_neg (hGRH : GRHFor chi4) (w : ℝ → ℝ)
    (hconv : ConvexOn ℝ (Set.Ici 3) w)
    (hsum : Summable fun p : Nat.Primes => w (p : ℕ)) (hw3 : 0 < w 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w (p : ℕ) < 0 := by
  sorry

/-- **Every convex weight, with the logarithm.**  Under GRH for `L(s, χ₋₄)`: let `f : ℕ → ℝ` be
convex on `n ≥ 3` with `Σ_p log p · f(p)` (absolutely) convergent and `f(3) > 0`.  Then
`Σ_p χ₋₄(p) log p · f(p) < 0`. -/
theorem grh_chi4_convex_log_weight_sum_neg (hGRH : GRHFor chi4) (f : ℕ → ℝ)
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) (hf3 : 0 < f 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p < 0 := by
  sorry

/-- **Fujii's weights.**  Under GRH for `L(s, χ₋₄)`: for every real `x > 0` and every exponent
`0 < α ≤ 1`,  `Σ_p χ₋₄(p) e^{−(p/x)^α} < 0`  (the series converges absolutely, `fujii_summable`;
`α = 1` is `grh_chebyshev_sum_neg`).  The weight `e^{−(t/x)^α}` is convex exactly because
`α ≤ 1`; for `α > 1` the statement is not claimed. -/
theorem grh_fujii_sum_neg (hGRH : GRHFor chi4) (x α : ℝ) (hx : 0 < x) (hα0 : 0 < α)
    (hα1 : α ≤ 1) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(((p : ℕ) : ℝ) / x) ^ α) < 0 := by
  sorry

/-- The series of `grh_fujii_sum_neg` converges absolutely for `x > 0`, `α > 0` (unconditional). -/
theorem fujii_summable (x α : ℝ) (hx : 0 < x) (hα0 : 0 < α) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(((p : ℕ) : ℝ) / x) ^ α)) := by
  sorry

/-- **The prime Dirichlet series on the real axis.**  Under GRH for `L(s, χ₋₄)`:
`Σ_p χ₋₄(p) p^{−σ} < 0` for every real `σ > 1` (`(p : ℝ) ^ (-σ)` is the real power). -/
theorem grh_chi4_prime_dirichlet_sum_neg (hGRH : GRHFor chi4) (σ : ℝ) (hσ : 1 < σ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * ((p : ℕ) : ℝ) ^ (-σ) < 0 := by
  sorry

/-- **The series identity behind the theorems (no hypothesis on zeros): second-order summation by
parts over the infinite range.**  For every `f : ℕ → ℝ` convex on `n ≥ 3` with
`Σ_p log p · |f(p)| < ∞`:

    `−Σ_p χ₋₄(p) log p · f(p) = Σ_{k ≥ 0} Θ₁(k+4) · (f(k+3) − 2 f(k+4) + f(k+5))`,

where `Θ₁(m) = Σ_{p ≤ m} (−χ₋₄(p)) log p · (m − p)` is the integrated θ-race of
`chebyshev_race_theta1_pos` (`monotonicity/lean-allx`), positive at every `m ≥ 4` under GRH; the
second differences are `≥ 0` by convexity.  This is the discrete form of
`−Σ_p χ(p) log p · F(p) = ∫ Θ₁(t) F''(t) dt` (Lemma 1 of
`monotonicity/chebyshev1853/CHEBYSHEV1853.md`); there are no boundary terms. -/
theorem chi4_log_sum_eq_theta1_series (f : ℕ → ℝ)
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) :
    HasSum (fun k : ℕ =>
        (∑ p ∈ (Icc 0 (k + 4)).filter Nat.Prime,
          -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (((k + 4 : ℕ) : ℝ) - p))
          * (f (k + 3) - 2 * f (k + 4) + f (k + 5)))
      (-∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p) := by
  sorry

/-- **The integral identity of the note (no hypothesis).**  For every real `x > 0`

    `Σ_p χ₋₄(p) e^{−p/x} = −∫_2^∞ Θ₁(t) F''(t) dt`,

where `Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p)` (the function of `chebyshev_race_theta1_pos`,
`monotonicity/lean-allx`: positive for every `t > 3` under GRH, and `0` on `[2, 3]`) and
`F''(t) = e^{−t/x} (1/(x² log t) + 2/(x t log² t) + 1/(t² log² t) + 2/(t² log³ t))` is the second
derivative of `F(t) = e^{−t/x}/log t` (`chebyshev_kernel_hasDerivAt`), positive on `(1, ∞)`
(`chebyshev_kernel_second_deriv_pos`).  The integral is a Lebesgue integral over `(2, ∞)` and the
integrand is integrable (`chebyshevSum_integrand_integrableOn`).  This is Lemma 1(ii) of
`monotonicity/chebyshev1853/CHEBYSHEV1853.md`. -/
theorem chebyshevSum_eq_neg_integral (x : ℝ) (hx : 0 < x) :
    chebyshevSum x = -∫ t in Set.Ioi (2 : ℝ),
      (∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p))
        * (Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
          + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3))) := by
  sorry

/-- The integrand of `chebyshevSum_eq_neg_integral` is Lebesgue integrable on `(2, ∞)`. -/
theorem chebyshevSum_integrand_integrableOn (x : ℝ) (hx : 0 < x) :
    IntegrableOn (fun t : ℝ =>
      (∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p))
        * (Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
          + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3)))) (Set.Ioi 2) := by
  sorry

/-- **The kernel.**  For `x ≠ 0` and `t > 1`: `F(t) = e^{−t/x}/log t` has derivative
`F'(t) = −e^{−t/x} (1/(x log t) + 1/(t log² t))`, and `F'` has derivative
`F''(t) = e^{−t/x} (1/(x² log t) + 2/(x t log² t) + 1/(t² log² t) + 2/(t² log³ t))`. -/
theorem chebyshev_kernel_hasDerivAt (x t : ℝ) (hx : x ≠ 0) (ht : 1 < t) :
    HasDerivAt (fun t : ℝ => Real.exp (-t / x) / Real.log t)
        (-Real.exp (-t / x) * (1 / (x * Real.log t) + 1 / (t * Real.log t ^ 2))) t ∧
      HasDerivAt
        (fun t : ℝ => -Real.exp (-t / x) * (1 / (x * Real.log t) + 1 / (t * Real.log t ^ 2)))
        (Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
          + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3))) t := by
  sorry

/-- `F''(t) > 0` for every `x > 0` and `t > 1`. -/
theorem chebyshev_kernel_second_deriv_pos (x t : ℝ) (hx : 0 < x) (ht : 1 < t) :
    0 < Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
      + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3)) := by
  sorry

/-- **From the sum without logarithms to the sums with logarithms (no hypothesis).**  Since
`1/log p = ∫_0^∞ p^{−u} du`, for every real `x > 0`

    `Σ_p χ₋₄(p) e^{−p/x} = ∫_0^∞ ( Σ_p χ₋₄(p) log p · p^{−u} e^{−p/x} ) du`

(Lebesgue integral over `u ∈ (0, ∞)`; `(p : ℝ) ^ (-u)` is the real power).  By
`grh_chi4_log_sum_rpow_exp_neg` every integrand is negative under GRH: Chebyshev's sum is an
average, with a positive measure, of logarithmically weighted sums.  (A variant of Lemma 5(iii) of
`monotonicity/chebyshev1853/CHEBYSHEV1853.md`, which uses Volterra's kernel instead.) -/
theorem chebyshevSum_eq_integral_log_sums (x : ℝ) (hx : 0 < x) :
    chebyshevSum x = ∫ u in Set.Ioi (0 : ℝ), ∑' p : Nat.Primes,
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
        * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  sorry

/-- Under GRH for `L(s, χ₋₄)`: for every real `x > 0` and `u ≥ 0`,
`Σ_p χ₋₄(p) log p · p^{−u} e^{−p/x} < 0`  (`u = 0` is `grh_chebyshev_log_sum_neg`). -/
theorem grh_chi4_log_sum_rpow_exp_neg (hGRH : GRHFor chi4) (x u : ℝ) (hx : 0 < x) (hu : 0 ≤ u) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x)) < 0 := by
  sorry

/-- **Second-order summation by parts over the infinite range, discrete form** (pure analysis; no
number theory).  Let `f : ℕ → ℝ` be a convex sequence tending to `0` and `a : ℕ → ℝ` any sequence
with `Σ_n |a(n)| f(n) < ∞`.  Then

    `Σ_n a(n) f(n) = Σ_{k ≥ 0} A₁(k+1) · (f(k) − 2 f(k+1) + f(k+2))`,   `A₁(k+1) = Σ_{n ≤ k} a(n) (k+1−n)`.

No boundary terms and no growth condition on `a`. -/
theorem second_order_summation_by_parts_infinite (a f : ℕ → ℝ)
    (hconv : ∀ k, 2 * f (k + 1) ≤ f k + f (k + 2))
    (h0 : Filter.Tendsto f Filter.atTop (nhds 0)) (hs : Summable fun n => |a n| * f n) :
    HasSum (fun k : ℕ => (∑ n ∈ range (k + 1), a n * ((k : ℝ) + 1 - n))
        * (f k - 2 * f (k + 1) + f (k + 2))) (∑' n, a n * f n) := by
  sorry

/-- **Second-order partial summation over the infinite range, integral form** (pure analysis;
Lemma 1 of `monotonicity/chebyshev1853/CHEBYSHEV1853.md` for convex test functions).  Let `a ≥ 0`,
let `F` be twice differentiable on `[a, ∞)` (`F' `, `F''` its derivatives) with `F''` continuous and
`≥ 0` there and `F(t) → 0` as `t → ∞`, and let `c : ℕ → ℝ` vanish below `a` with
`Σ_n |c(n)| F(n) < ∞`.  Then, with `C₁(t) = Σ_{n ≤ t} c(n) (t − n)`:

    `C₁ · F''` is Lebesgue integrable on `(a, ∞)`   and   `Σ_n c(n) F(n) = ∫_a^∞ C₁(t) F''(t) dt`.

No hypothesis on the size of `c`, `C₁` or `F'` is needed: the boundary terms `F(T)·Σ_{n≤T} c(n)` and
`F'(T)·C₁(T)` never appear (Tonelli for the kernels `(t − n)₊ F''(t)`). -/
theorem second_order_partial_summation_infinite (a : ℝ) (ha : 0 ≤ a) (c : ℕ → ℝ)
    (hc : ∀ n : ℕ, (n : ℝ) < a → c n = 0) (F F' F'' : ℝ → ℝ)
    (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t) (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t)
    (hF'' : ContinuousOn F'' (Set.Ici a)) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Filter.Tendsto F Filter.atTop (nhds 0)) (hs : Summable fun n : ℕ => |c n| * F n) :
    IntegrableOn (fun t => (∑ n ∈ Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t) (Set.Ioi a) ∧
      ∑' n, c n * F n = ∫ t in Set.Ioi a, (∑ n ∈ Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t := by
  sorry

/-- The series defining Chebyshev's sum `chebyshevSum x` converges absolutely for `x > 0`
(unconditional; from `monotonicity/lean-chebyshev`).  (Lean's `∑'` returns `0` for a divergent
series; this excludes that reading of `grh_chebyshev_sum_neg`.) -/
theorem chebyshevSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes =>
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  sorry

/-- The series defining `chebyshevLogSum x` converges absolutely for `x > 0` (unconditional; from
`monotonicity/lean-chebyshev`). -/
theorem chebyshevLogSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  sorry

/-- `χ₋₄(2) = 0`. -/
theorem chi4_natCast_two : ((ZMod.χ₄ ((2 : ℕ) : ZMod 4) : ℤ) : ℝ) = 0 := by
  sorry

/-- `χ₋₄(n) = (−1)^{(n−1)/2}` for odd `n` (Chebyshev's sign: `+1` for `n ≡ 1`, `−1` for `n ≡ 3`
mod 4). -/
theorem chi4_natCast_of_odd (n : ℕ) (hn : n % 2 = 1) :
    ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = (-1) ^ ((n - 1) / 2) := by
  sorry

/-- **The shape of the hypothesis.**  `GRHFor χ` is the instance at the single character `χ` of the
statement "for every modulus `N ≥ 1` and every Dirichlet character `χ` mod `N`, every zero of
`L(s,χ)` with `0 < Re s < 1` has `Re s = 1/2`": a proof of the latter (with exactly this binder
structure) discharges the hypothesis `GRHFor χ` of every theorem above, for every `χ`. -/
theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ := by
  sorry
