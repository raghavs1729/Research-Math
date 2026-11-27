/-
RamanujanStatement.lean — the TRUSTED statement module: WHAT IS CLAIMED.

This module imports ONLY Mathlib and contains no definitions of its own: the statements are
written with Mathlib's `Nat.primeCounting` (the number of primes `≤ n`), `Nat.floor` (`⌊x⌋₊`),
`Real.log`, `Real.exp`, `Chebyshev.theta` and (last statement only) `RiemannHypothesis`.  The
first four statements have no hypotheses.  The last two are conditional; their hypotheses are
explicit binders, described in the docstrings.  There are no `axiom`s.
`RamanujanSolution.lean` (untrusted) proves exactly these statements from `RamanujanLib/` and the
Zeta23 library.  The `sorry`s below are deliberate (this is the statement side).
-/
import Mathlib

/-- **Ramanujan's inequality for the prime counting function.**

Write `π(x)` for the number of primes `p ≤ x` (`x` real).  Ramanujan (second notebook; Berndt,
*Ramanujan's Notebooks*, Part IV, Chapter 24) asserted that

  `π(x)² < (e · x / log x) · π(x / e)`   for all sufficiently large `x`.

This is that statement, for real `x`: there is a real `x₀` such that the strict inequality holds
for every real `x ≥ x₀`.  In Lean, `π(x)` is `Nat.primeCounting ⌊x⌋₊` (primes `≤ ⌊x⌋`, which are the
primes `≤ x` for `x ≥ 0`), `e` is `Real.exp 1` and `log` is the natural logarithm.  The statement
is unconditional; `x₀` is not explicit. -/
theorem ramanujan_prime_counting_inequality :
    ∃ x₀ : ℝ, ∀ x : ℝ, x₀ ≤ x →
      (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        < Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ) := by
  sorry

/-- **Ramanujan's asymptotic for the gap** (the computation behind the inequality):

  `(e x / log x) · π(x / e) − π(x)² = (1 + o(1)) · x² / log⁶ x`   as `x → ∞`,

stated as: the quotient of the left side by `x² / log⁶ x` tends to `1`.  Since the limit is
positive, this implies `ramanujan_prime_counting_inequality`, and it shows that the inequality
holds with a margin of exact order `x² / log⁶ x`. -/
theorem ramanujan_gap_asymptotic :
    Filter.Tendsto (fun x : ℝ =>
      (Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ)
        - (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2) / (x ^ 2 / Real.log x ^ 6))
      Filter.atTop (nhds 1) := by
  sorry

/-- **Five-term expansion of the prime counting function** (the input of Ramanujan's argument):
with `L = log x`,

  `|π(x) − x (1/L + 1/L² + 2/L³ + 6/L⁴ + 24/L⁵)| ≤ C · x / L⁶`   for every real `x ≥ 2`,

for an absolute (non-explicit) constant `C`. -/
theorem primeCounting_five_term_expansion :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ)
          - x * (1 / Real.log x + 1 / Real.log x ^ 2 + 2 / Real.log x ^ 3 + 6 / Real.log x ^ 4
            + 24 / Real.log x ^ 5)|
        ≤ C * (x / Real.log x ^ 6) := by
  sorry

/-- **Prime number theorem with a `log⁻⁵` saving, on all of `[2, ∞)`:**
`|θ(x) − x| ≤ K · x / log⁵ x` for every real `x ≥ 2`, where `θ(x) = Σ_{p ≤ x} log p` and `K` is an
absolute (non-explicit) constant. -/
theorem chebyshev_theta_sub_le :
    ∃ K : ℝ, ∀ x : ℝ, 2 ≤ x → |Chebyshev.theta x - x| ≤ K * (x / Real.log x ^ 5) := by
  sorry

/-- **Explicit criterion: a Schoenfeld-type bound gives Ramanujan's inequality for `x ≥ e⁴⁰`.**

Hypothesis `hS`: for every real `x ≥ 2657`,
`|π(x) − (c + ∫_2^x dt / log t)| < √x · log x / (8π)`, where `c` is any constant with `|c| ≤ 2`.
With `c = li 2 = 1.04516…` the function `c + ∫_2^x dt / log t` is the logarithmic integral `li x`
and `hS` is Schoenfeld's bound (Math. Comp. 30 (1976), Corollary 1), which he proves assuming the
Riemann Hypothesis.  (The principal-value integral `li` is not in Mathlib, which is why the constant
is a parameter.)

Conclusion: `π(x)² < (e x / log x) · π(x / e)` for every real `x ≥ e⁴⁰ ≈ 2.35 · 10¹⁷`.

The implication itself is unconditional and fully proved: it uses neither the Riemann Hypothesis
nor the prime number theorem.  For comparison, the threshold of Hassani under the Riemann
Hypothesis is `138 766 146 692 471 228 ≈ 1.39 · 10¹⁷`. -/
theorem ramanujan_of_schoenfeld_bound (c : ℝ) (hc : |c| ≤ 2)
    (hS : ∀ x : ℝ, 2657 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..x, 1 / Real.log t)|
        < √x * Real.log x / (8 * Real.pi)) :
    ∀ x : ℝ, Real.exp 40 ≤ x →
      (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        < Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ) := by
  sorry

/-- **Reduction of the sharp statement under the Riemann Hypothesis.**

Three hypotheses, all explicit binders (none is proved here, none is an axiom):

* `hRH`: the Riemann Hypothesis, in Mathlib's form `RiemannHypothesis`;
* `hSchoenfeld`: the Riemann Hypothesis implies Schoenfeld's bound
  `|π(x) − (c + ∫_2^x dt / log t)| < √x log x / (8π)` for `x ≥ 2657` (a literature theorem:
  Schoenfeld 1976, Corollary 1, with `c = li 2`);
* `hFinite`: Ramanujan's inequality at every integer `n` with `38 358 837 683 ≤ n < e⁴⁰`
  (a statement about finitely many integers).

Conclusion: `π(x)² < (e x / log x) · π(x / e)` for every real `x ≥ 38 358 837 683`, the sharp
threshold of Dudek and Platt (Exp. Math. 24 (2015)), who prove it under the Riemann Hypothesis.

Status of `hFinite` in the literature (as recorded in `../RAMANUJAN-TRIAGE.md`, not checked
here): Dudek and Platt verified the integers up to `1.15 · 10¹⁶` by computation, and Axler
(arXiv:2207.02486) proves the inequality unconditionally for `38 358 837 683 ≤ x ≤ 10¹⁹`, which
contains the range `[38 358 837 683, e⁴⁰)`, `e⁴⁰ ≈ 2.35 · 10¹⁷`. -/
theorem ramanujan_sharp_under_RH_reduction (c : ℝ) (hc : |c| ≤ 2)
    (hRH : RiemannHypothesis)
    (hSchoenfeld : RiemannHypothesis → ∀ x : ℝ, 2657 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..x, 1 / Real.log t)|
        < √x * Real.log x / (8 * Real.pi))
    (hFinite : ∀ n : ℕ, 38358837683 ≤ n → (n : ℝ) < Real.exp 40 →
      (Nat.primeCounting n : ℝ) ^ 2
        < Real.exp 1 * (n : ℝ) / Real.log (n : ℝ)
          * (Nat.primeCounting ⌊(n : ℝ) / Real.exp 1⌋₊ : ℝ)) :
    ∀ x : ℝ, 38358837683 ≤ x →
      (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        < Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ) := by
  sorry
