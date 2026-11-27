/-
ChallengeDeps.lean — the TRUSTED definition layer of the smoothed prime race statements.

Everything the statements in Challenge.lean mention is defined HERE, from Mathlib alone.  This
module imports nothing from the Zeta23 library.  A reader auditing WHAT is claimed needs to read
only this file and Challenge.lean (trusting Mathlib and the Lean kernel); the proofs live in
RaceLib/ and in the Zeta23 library (anthropics/formal-math, folder zeta23/).

Conventions (they are those of Zeta23's literature-form explicit formula `EF_lit_chi`, which
transcribes [Iwaniec–Kowalski, Thm 5.12] for a primitive Dirichlet character χ mod q > 1):
  * Fourier transform  ĝ(z) = ∫ g(u) e^{izu} du  (no 2π in the exponent; complex argument z);
  * a nontrivial zero ρ of L(s,χ) is written ρ = 1/2 + iγ_ρ, i.e. γ_ρ = (ρ − 1/2)/i (complex in
    general; real exactly when Re ρ = 1/2), and is counted with its multiplicity m_ρ;
  * the archimedean weight is  Re ψ(1/4 + κ/2 + ir/2) + log(q/π),  ψ = Γ'/Γ, κ = 0 for even χ and
    κ = 1 for odd χ.
-/
import Mathlib

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory

noncomputable section

/-! ## 1. Test functions -/

/-- Fourier transform with the convention `ĝ(z) = ∫ g(u) e^{izu} du`, for complex `z`. -/
def raceFT (g : ℝ → ℂ) (z : ℂ) : ℂ := ∫ u : ℝ, g u * Complex.exp (Complex.I * z * (u : ℂ))

/-- The autocorrelation `g = f ⋆ f̃` of `f`, where `f̃(u) = conj (f (−u))`:
`g(x) = ∫ f(t) conj(f(t − x)) dt`.  Its Fourier transform is `|f̂|² ≥ 0` on the real line. -/
def autocorr (f : ℝ → ℂ) : ℝ → ℂ := fun x => ∫ t : ℝ, f t * conj (f (t - x))

/-! ## 2. The hypothesis -/

/-- The Riemann hypothesis for the single Dirichlet L-function `L(s,χ)`: every zero in the open
critical strip `0 < Re s < 1` has real part `1/2`.  This is a hypothesis of the theorems, passed as
an explicit argument; it is never an axiom.  (Its binder structure is that of the statement
`∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ), 0 < s.re → s.re < 1 →
DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2`, instantiated at one character.) -/
def GRHFor {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → χ.LFunction s = 0 → s.re = 1 / 2

/-! ## 3. The three sides of the explicit formula -/

open scoped Classical in
/-- The parity `κ` of `χ` as a natural number: `0` if `χ` is even (`χ(−1) = 1`), else `1`. -/
def archParity {q : ℕ} (χ : DirichletCharacter ℂ q) : ℕ := if χ.Even then 0 else 1

/-- The archimedean weight `Re ψ(1/4 + κ/2 + ir/2) + log(q/π)` of the explicit formula for
`L(s,χ)`, where `ψ = Γ'/Γ` is Mathlib's `Complex.digamma`. -/
def archWeight {q : ℕ} (χ : DirichletCharacter ℂ q) (r : ℝ) : ℝ :=
  (Complex.digamma (1 / 4 + (archParity χ : ℂ) / 2 + I * r / 2)).re + Real.log (q / Real.pi)

/-- The archimedean term `Γ_L` of the explicit formula for the shifted test function
`u ↦ g(u − L)`:  `Γ_L = (1/2π) ∫ e^{irL} ĝ(r) [Re ψ(1/4 + κ/2 + ir/2) + log(q/π)] dr`. -/
def archTerm {q : ℕ} (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) (L : ℝ) : ℂ :=
  (1 / (2 * Real.pi) : ℂ) *
    ∫ r : ℝ, Complex.exp (I * r * L) * raceFT g r * (archWeight χ r : ℂ)

/-- The smoothed twisted sum over prime powers
`S(L) = Σ_{n ≥ 1} Λ(n) χ(n) n^{-1/2} g(log n − L)`  (`Λ(0) = 0`, so `n = 0` contributes nothing). -/
def primePowerSum {q : ℕ} (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) (L : ℝ) : ℂ :=
  ∑' n : ℕ, ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ) * χ (n : ZMod q)
    * g (Real.log n - L)

/-- The two-sided prime-power sum of the explicit formula at shift `0`:
`Σ_{n ≥ 1} Λ(n) n^{-1/2} (χ(n) g(log n) + conj(χ(n)) g(−log n))`. -/
def primePowerSumSym {q : ℕ} (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) : ℂ :=
  ∑' n : ℕ, ((ArithmeticFunction.vonMangoldt n / Real.sqrt n : ℝ) : ℂ)
    * (χ (n : ZMod q) * g (Real.log n) + conj (χ (n : ZMod q)) * g (-Real.log n))

/-- The sum over the distinct nontrivial zeros `ρ` of `L(s,χ)` (zeros with `0 < Re ρ < 1`), each
weighted by its multiplicity `m_ρ` (the analytic order of vanishing), of `ĝ(γ_ρ)`,
`γ_ρ = (ρ − 1/2)/i`:   `Σ_ρ m_ρ ĝ(γ_ρ)`. -/
def zeroSum {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) : ℂ :=
  ∑' ρ : {ρ : ℂ | χ.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1},
    ((analyticOrderAt χ.LFunction (ρ : ℂ)).toNat : ℂ) * raceFT g (((ρ : ℂ) - 1 / 2) / I)

/-- The smoothed twisted sum over primes only:
`P(L) = Σ_p χ(p) (log p) p^{-1/2} g(log p − L)`. -/
def primeSum {q : ℕ} (χ : DirichletCharacter ℂ q) (g : ℝ → ℂ) (L : ℝ) : ℂ :=
  ∑' p : Nat.Primes, χ ((p : ℕ) : ZMod q) * ((Real.log (p : ℕ) / Real.sqrt (p : ℕ) : ℝ) : ℂ)
    * g (Real.log (p : ℕ) - L)

/-! ## 4. The character χ₋₄ and the prime number theorem -/

/-- The nontrivial Dirichlet character modulo 4 (`χ₋₄(n) = 1, −1, 0` for `n ≡ 1, 3 mod 4` and `n`
even), as a `ℂ`-valued Dirichlet character, from Mathlib's `ZMod.χ₄ : MulChar (ZMod 4) ℤ`. -/
def chi4 : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)

/-- The prime number theorem in the form `θ(x) ~ x`, i.e. `θ(x)/x → 1`, where
`θ(x) = Σ_{p ≤ x} log p` is Mathlib's `Chebyshev.theta`.  (It is PROVED in this project, as the
statement `prime_number_theorem_theta`, from the Zeta23 library's port of the PNT+ project; the
Mathlib revision pinned here does not contain it.) -/
def PNTTheta : Prop :=
  Filter.Tendsto (fun x : ℝ => Chebyshev.theta x / x) Filter.atTop (nhds 1)

/-! ## 5. The race for a real-valued window, in real terms -/

/-- For a real-valued window `f`, with `g = f ⋆ f̃` (so `g(x) = ∫ f(t) f(t − x) dt`), the smoothed
`χ₋₄`-twisted sum over primes
`Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L)`,  a real number.
It is `(sum over p ≡ 1 mod 4) − (sum over p ≡ 3 mod 4)` of the positive weights
`(log p) p^{-1/2} g(log p − L)` when `g ≥ 0`; "negative" means: the class `3 mod 4` is ahead. -/
def chi4RaceSum (f : ℝ → ℝ) (L : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
    * (Real.log (p : ℕ) / Real.sqrt (p : ℕ)) * ∫ t : ℝ, f t * f (t - (Real.log (p : ℕ) - L))

/-- For a real-valued window `f`: the archimedean term at shift `0` for `χ₋₄` and `g = f ⋆ f̃`,
`Γ_0 = (1/2π) ∫ |f̂(r)|² [Re ψ(3/4 + ir/2) + log(4/π)] dr`,  a real number. -/
def chi4Gamma0 (f : ℝ → ℝ) : ℝ :=
  (1 / (2 * Real.pi)) * ∫ r : ℝ, ‖raceFT (fun x => (f x : ℂ)) r‖ ^ 2
    * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi))

end
