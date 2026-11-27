/-
RaceSolution.lean — the UNTRUSTED solution module: the statements of RaceChallenge.lean, byte-identical,
each PROVED by delegating to RaceLib (which is built on the Zeta23 library, anthropics/formal-math,
folder zeta23/, and on Mathlib).

  explicit_formula_zero_sum                →  RaceLib.zeroSum_eq  (= Zeta23.ThmE.EF_lit_chi_L)
  zeroSum_summable_and_archTerm_integrable →  RaceLib.summable_zeroSum, .integrable_raceFT_mul_archWeight
  grh_smoothed_prime_power_sum_bound       →  RaceLib.norm_primePowerSum_sub_archTerm_le
  grh_smoothed_prime_power_sum_le_zero_sum →  RaceLib.norm_primePowerSum_sub_archTerm_le_zeroSum
  grhFor_of_dirichlet_zeros_on_line        →  (definitional)
  chi4_isPrimitive, chi4_archWeight        →  RaceLib.chi4_isPrimitive, RaceLib.archWeight_chi4
  grh_chi4_smoothed_prime_power_sum_bound  →  RaceLib.chi4_norm_primePowerSum_sub_archTerm_le
  prime_number_theorem_theta               →  RaceLib.pntTheta  (from Zeta23's `MediumPNT`)
  grh_chi4_smoothed_race_two_sided         →  RaceLib.chi4_norm_primeSum_add_le
  grh_chi4_smoothed_race                   →  RaceLib.chi4_primeSum_re_lt
  chi4_primeSum_eq_chi4RaceSum,
  chi4_archTerm_zero_eq_chi4Gamma0         →  RaceLib.primeSum_chi4_ofReal, .archTerm_chi4_zero_ofReal
  chi4Gamma0_integrand_integrable          →  RaceLib.integrable_chi4Gamma0_integrand
  grh_chi4_smoothed_race_real_two_sided    →  RaceLib.abs_chi4RaceSum_add_le
  grh_chi4_smoothed_race_real              →  RaceLib.chi4RaceSum_lt
  raceWindow_piecewise, raceWindow_admissible,
  raceWindow_integral, raceWindowCert_iff  →  RaceLib/Window.lean
  grh_chi4_race_concrete_window            →  RaceLib.raceWindow_race
  grh_chi4_race_concrete_window_two_sided  →  RaceLib.raceWindow_two_sided

The definitions used in the statements are those of RaceChallengeDeps.lean (Mathlib only; inlined
verbatim in RaceChallenge.lean).  Nothing in this file is part of the trusted base.
-/
import RaceChallengeDeps
import RaceLib

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory

/-- **The explicit formula at shift 0** (unconditional; a restatement, in the notation of this
project, of Zeta23's `EF_lit_chi_L`).  For a primitive Dirichlet character `χ` modulo `q > 1` and a
test function `g ∈ C²_c(ℝ)`:

    Σ_ρ m_ρ ĝ(γ_ρ) = Γ_0 − Σ_{n ≥ 1} Λ(n) n^{-1/2} (χ(n) g(log n) + conj(χ(n)) g(−log n)),

the sum on the left running over the distinct zeros `ρ = 1/2 + iγ_ρ` of `L(s,χ)` in the open
critical strip, with multiplicity `m_ρ`, and `Γ_0 = (1/2π) ∫ ĝ(r) [Re ψ(1/4 + κ/2 + ir/2) +
log(q/π)] dr`. -/
theorem explicit_formula_zero_sum {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) (hq : 1 < q)
    (hχ : χ.IsPrimitive) (g : ℝ → ℂ) (hg : ContDiff ℝ 2 g) (hgc : HasCompactSupport g) :
    zeroSum χ g = archTerm χ g 0 - primePowerSumSym χ g := by
  exact RaceLib.zeroSum_eq hq hχ hg hgc

/-- **The sums and integrals above are genuine.**  For a primitive `χ` modulo `q > 1` and
`g ∈ C²_c(ℝ)`: the zero sum `Σ_ρ m_ρ ĝ(γ_ρ)` defining `zeroSum χ g` converges absolutely, and the
archimedean integrand `ĝ(r) [Re ψ(1/4 + κ/2 + ir/2) + log(q/π)]` defining `archTerm χ g L` is
Lebesgue integrable (multiplying it by the unimodular factor `e^{irL}` does not change that).
(Lean's `∑'` and `∫` return `0` for divergent sums and non-integrable functions; this statement
excludes that reading.  The prime sums have finitely many nonzero terms.) -/
theorem zeroSum_summable_and_archTerm_integrable {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q)
    (hq : 1 < q) (hχ : χ.IsPrimitive) (g : ℝ → ℂ) (hg : ContDiff ℝ 2 g)
    (hgc : HasCompactSupport g) :
    Summable (fun ρ : {ρ : ℂ | χ.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1} =>
        ((analyticOrderAt χ.LFunction (ρ : ℂ)).toNat : ℂ) * raceFT g (((ρ : ℂ) - 1 / 2) / I))
      ∧ Integrable (fun r : ℝ => raceFT g r * (archWeight χ r : ℂ)) := by
  exact ⟨RaceLib.summable_zeroSum hq hχ hg hgc, RaceLib.integrable_raceFT_mul_archWeight χ hg hgc⟩

/-- **S0: the GRH bound for a smoothed twisted prime-power sum.**  Let `χ` be a primitive Dirichlet
character modulo `q > 1` and ASSUME the Riemann hypothesis for `L(s,χ)` (`GRHFor χ`).  Let
`f ∈ C²_c(ℝ)` be supported in `[−A/2, A/2]` and let `g = f ⋆ f̃` be its autocorrelation (so `g` is
supported in `[−A, A]` and `ĝ = |f̂|² ≥ 0` on the real line).  Then for every shift `L > A`:

    ‖S(L) − Γ_L‖ ≤ Re ( Γ_0 − Σ_{n ≥ 1} Λ(n) n^{-1/2} (χ(n) g(log n) + conj(χ(n)) g(−log n)) ),

where `S(L) = Σ_n Λ(n) χ(n) n^{-1/2} g(log n − L)` and
`Γ_L = (1/2π) ∫ e^{irL} ĝ(r) [Re ψ(1/4 + κ/2 + ir/2) + log(q/π)] dr`.
The right-hand side is the zero sum `Σ_ρ m_ρ ĝ(γ_ρ)` at shift `0` (see
`explicit_formula_zero_sum`); no information about individual zeros, no simplicity and no linear
independence of the ordinates is used. -/
theorem grh_smoothed_prime_power_sum_bound {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q)
    (hq : 1 < q) (hχ : χ.IsPrimitive) (hGRH : GRHFor χ) (f : ℝ → ℂ) (hf : ContDiff ℝ 2 f)
    (hfc : HasCompactSupport f) (A : ℝ) (hA : tsupport f ⊆ Set.Icc (-(A / 2)) (A / 2)) (L : ℝ)
    (hL : A < L) :
    ‖primePowerSum χ (autocorr f) L - archTerm χ (autocorr f) L‖
      ≤ (archTerm χ (autocorr f) 0 - primePowerSumSym χ (autocorr f)).re := by
  exact RaceLib.norm_primePowerSum_sub_archTerm_le hq hχ hGRH hf hfc hA hL

/-- **S0, zero-sum form.**  Same hypotheses as `grh_smoothed_prime_power_sum_bound`:
`‖S(L) − Γ_L‖ ≤ Re Σ_ρ m_ρ ĝ(γ_ρ)`  (multiplicities included). -/
theorem grh_smoothed_prime_power_sum_le_zero_sum {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q)
    (hq : 1 < q) (hχ : χ.IsPrimitive) (hGRH : GRHFor χ) (f : ℝ → ℂ) (hf : ContDiff ℝ 2 f)
    (hfc : HasCompactSupport f) (A : ℝ) (hA : tsupport f ⊆ Set.Icc (-(A / 2)) (A / 2)) (L : ℝ)
    (hL : A < L) :
    ‖primePowerSum χ (autocorr f) L - archTerm χ (autocorr f) L‖
      ≤ (zeroSum χ (autocorr f)).re := by
  exact RaceLib.norm_primePowerSum_sub_archTerm_le_zeroSum hq hχ hGRH hf hfc hA hL

/-- **The shape of the hypothesis.**  `GRHFor χ` is the instance at the single character `χ` of the
statement "for every modulus `N ≥ 1` and every Dirichlet character `χ` mod `N`, every zero of
`L(s,χ)` with `0 < Re s < 1` has `Re s = 1/2`": a proof of the latter (with exactly this binder
structure) discharges the hypothesis `GRHFor χ` of every theorem below, for every `χ`. -/
theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ := by
  exact fun s h0 h1 hz => h q χ s h0 h1 hz

/-- `χ₋₄` is a primitive Dirichlet character modulo 4. -/
theorem chi4_isPrimitive : chi4.IsPrimitive := by
  exact RaceLib.chi4_isPrimitive

/-- The archimedean weight of `χ₋₄` (an odd character of modulus 4) is
`Re ψ(3/4 + ir/2) + log(4/π)`. -/
theorem chi4_archWeight (r : ℝ) :
    archWeight chi4 r = (Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi) := by
  exact RaceLib.archWeight_chi4 r

/-- **S0 for χ₋₄, closed form.**  ASSUME the Riemann hypothesis for `L(s,χ₋₄)`.  Let `f ∈ C²_c(ℝ)`
be supported in `[−a/2, a/2]` with `a < log 3`, and `g = f ⋆ f̃`.  Then for every shift `L > a`:

    ‖ Σ_n Λ(n) χ₋₄(n) n^{-1/2} g(log n − L) − Γ_L ‖ ≤ Re Γ_0,

with `Γ_L = (1/2π) ∫ e^{irL} ĝ(r) [Re ψ(3/4 + ir/2) + log(4/π)] dr` (see `chi4_archWeight`).
(The prime side of the explicit formula at shift `0` is empty because `g` vanishes at `±log n` for
`n ≥ 3`, `Λ(1) = 0` and `χ₋₄(2) = 0`.) -/
theorem grh_chi4_smoothed_prime_power_sum_bound (hGRH : GRHFor chi4) (f : ℝ → ℂ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (L : ℝ) (hL : a < L) :
    ‖primePowerSum chi4 (autocorr f) L - archTerm chi4 (autocorr f) L‖
      ≤ (archTerm chi4 (autocorr f) 0).re := by
  exact RaceLib.chi4_norm_primePowerSum_sub_archTerm_le hGRH hf hfc ha hsupp hL

/-- **The prime number theorem**, in the form `θ(x)/x → 1` (`θ(x) = Σ_{p ≤ x} log p`).  Proved here
from `MediumPNT` (`ψ(x) − x = O(x exp(−c (log x)^{1/10}))`), which the Zeta23 library ports from
the PrimeNumberTheoremAnd project, and Mathlib's `ψ − θ = O(√x)`.  It is the input for the prime
squares in the race statements below. -/
theorem prime_number_theorem_theta : PNTTheta := by
  exact RaceLib.pntTheta

/-- **S1: the smoothed race, two-sided form.**  ASSUME the Riemann hypothesis for `L(s,χ₋₄)`
(`GRHFor chi4`).  Let `f ∈ C²_c(ℝ)` be supported in `[−a/2, a/2]` with `a < log 3`, and
`g = f ⋆ f̃`.  Then for every `ε > 0` and all sufficiently large `L`:

    ‖ Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L) + (1/2) ĝ(0) ‖ ≤ Re Γ_0 + ε.

(`ĝ(0) = ∫ g = |∫ f|²`; the term `−(1/2) ĝ(0)` is the contribution of the squares of primes, the
source of Chebyshev's bias; `Γ_0`, see `grh_chi4_smoothed_prime_power_sum_bound`, bounds the
contribution of ALL the zeros of `L(s,χ₋₄)`, with multiplicity.) -/
theorem grh_chi4_smoothed_race_two_sided (hGRH : GRHFor chi4) (f : ℝ → ℂ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (ε : ℝ) (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀,
      ‖primeSum chi4 (autocorr f) L + (1 / 2 : ℂ) * raceFT (autocorr f) 0‖
        ≤ (archTerm chi4 (autocorr f) 0).re + ε := by
  exact RaceLib.chi4_norm_primeSum_add_le hGRH RaceLib.pntTheta hf hfc ha hsupp hε

/-- **S1: the smoothed race** (complex-valued window).  ASSUME the Riemann hypothesis for
`L(s,χ₋₄)` and — for the chosen window `f ∈ C²_c(ℝ)`, supported in `[−a/2, a/2]` with `a < log 3`,
`g = f ⋆ f̃` — the numerical inequality

    hcert :  Re Γ_0 < (1/2) Re ĝ(0) − margin

(a statement about explicit integrals, to be supplied by a certified computation; it is an explicit
hypothesis, not an axiom).  Then for all sufficiently large `L`:

    Re Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L) < −margin. -/
theorem grh_chi4_smoothed_race (hGRH : GRHFor chi4) (f : ℝ → ℂ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (margin : ℝ)
    (hcert : (archTerm chi4 (autocorr f) 0).re
      < (1 / 2) * (raceFT (autocorr f) 0).re - margin) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, (primeSum chi4 (autocorr f) L).re < -margin := by
  exact RaceLib.chi4_primeSum_re_lt hGRH RaceLib.pntTheta hf hfc ha hsupp hcert

/-- For a real-valued window `f`, the complex smoothed prime sum for `χ₋₄` and `g = f ⋆ f̃` is the
real number `chi4RaceSum f L`. -/
theorem chi4_primeSum_eq_chi4RaceSum (f : ℝ → ℝ) (L : ℝ) :
    primeSum chi4 (autocorr fun x => (f x : ℂ)) L = ((chi4RaceSum f L : ℝ) : ℂ) := by
  exact RaceLib.primeSum_chi4_ofReal f L

/-- For a real-valued continuous compactly supported window `f`, the archimedean term `Γ_0` for
`χ₋₄` and `g = f ⋆ f̃` is the real number `chi4Gamma0 f`. -/
theorem chi4_archTerm_zero_eq_chi4Gamma0 (f : ℝ → ℝ) (hf : Continuous f)
    (hfc : HasCompactSupport f) :
    archTerm chi4 (autocorr fun x => (f x : ℂ)) 0 = ((chi4Gamma0 f : ℝ) : ℂ) := by
  exact RaceLib.archTerm_chi4_zero_ofReal hf hfc

/-- For a real-valued `C²_c` window `f`, the integrand defining `chi4Gamma0 f` is Lebesgue
integrable (so `chi4Gamma0 f` is a genuine integral). -/
theorem chi4Gamma0_integrand_integrable (f : ℝ → ℝ) (hf : ContDiff ℝ 2 f)
    (hfc : HasCompactSupport f) :
    Integrable (fun r : ℝ => ‖raceFT (fun x => (f x : ℂ)) r‖ ^ 2
      * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi))) := by
  exact RaceLib.integrable_chi4Gamma0_integrand hf hfc

/-- **S1 for a real-valued window, two-sided form.**  ASSUME the Riemann hypothesis for
`L(s,χ₋₄)`.  Let `f : ℝ → ℝ` be `C²`, supported in `[−a/2, a/2]` with `a < log 3`, and
`g(x) = ∫ f(t) f(t − x) dt`.  Then for every `ε > 0` and all sufficiently large `L`:

    | Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L) + (1/2) (∫ f)² |  ≤  Γ_0 + ε,

`Γ_0 = (1/2π) ∫ |f̂(r)|² [Re ψ(3/4 + ir/2) + log(4/π)] dr`. -/
theorem grh_chi4_smoothed_race_real_two_sided (hGRH : GRHFor chi4) (f : ℝ → ℝ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (ε : ℝ) (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, |chi4RaceSum f L + (1 / 2) * (∫ x, f x) ^ 2| ≤ chi4Gamma0 f + ε := by
  exact RaceLib.abs_chi4RaceSum_add_le hGRH RaceLib.pntTheta hf hfc ha hsupp hε

/-- **S1: the smoothed Chebyshev race modulo 4** (real-valued window; the headline).
ASSUME the Riemann hypothesis for the single L-function `L(s,χ₋₄)` (`hGRH`).  Let `f : ℝ → ℝ` be a
`C²` window supported in `[−a/2, a/2]` with `a < log 3`, let `g(x) = ∫ f(t) f(t − x) dt`, and ASSUME
the numerical inequality for this window

    hcert :  Γ_0 < (1/2) (∫ f)² − margin,
    Γ_0 = (1/2π) ∫ |f̂(r)|² [Re ψ(3/4 + ir/2) + log(4/π)] dr

(to be supplied by a certified computation; an explicit hypothesis, not an axiom).  Then there is
`L₀` such that for every `L ≥ L₀`:

    Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L)  <  −margin.

With `margin ≥ 0`: for all large `L` the smoothed sum is strictly negative.  When `f ≥ 0` the
weights `(log p) p^{-1/2} g(log p − L)` are nonnegative, and this says that the primes `≡ 3 mod 4`
are ahead of the primes `≡ 1 mod 4` in the smoothed race, for EVERY large `L` (no exceptional set,
no linear-independence or simplicity hypothesis on the zeros). -/
theorem grh_chi4_smoothed_race_real (hGRH : GRHFor chi4) (f : ℝ → ℝ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (margin : ℝ)
    (hcert : chi4Gamma0 f < (1 / 2) * (∫ x, f x) ^ 2 - margin) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum f L < -margin := by
  exact RaceLib.chi4RaceSum_lt hGRH RaceLib.pntTheta hf hfc ha hsupp hcert

/-- **The concrete window is what it should be.**  `raceWindow` is even, equals `1` on
`[−49/100, 49/100]`, equals the quintic ramp `10t³ − 15t⁴ + 6t⁵`, `t = 20(27/50 − u)`, on
`[49/100, 27/50]`, and vanishes for `|u| ≥ 27/50`.  (This is the function whose `Γ_0` the script
`cert/window_cert.py` bounds.) -/
theorem raceWindow_piecewise (u : ℝ) :
    raceWindow (-u) = raceWindow u
      ∧ (|u| ≤ 49 / 100 → raceWindow u = 1)
      ∧ (49 / 100 ≤ u → u ≤ 27 / 50 →
          raceWindow u = 10 * (20 * (27 / 50 - u)) ^ 3 - 15 * (20 * (27 / 50 - u)) ^ 4
            + 6 * (20 * (27 / 50 - u)) ^ 5)
      ∧ (27 / 50 ≤ |u| → raceWindow u = 0) := by
  exact RaceLib.raceWindow_piecewise u

/-- **The concrete window is admissible**: it is `C²`, compactly supported, supported in
`[−a/2, a/2]` with `a = 27/25`, and `27/25 < log 3`. -/
theorem raceWindow_admissible :
    ContDiff ℝ 2 raceWindow ∧ HasCompactSupport raceWindow
      ∧ tsupport raceWindow ⊆ Set.Icc (-((27 / 25 : ℝ) / 2)) ((27 / 25) / 2)
      ∧ (27 / 25 : ℝ) < Real.log 3 := by
  exact ⟨RaceLib.contDiff_raceWindow, RaceLib.hasCompactSupport_raceWindow,
    RaceLib.tsupport_raceWindow_subset', RaceLib.log_three_gt⟩

/-- `∫ raceWindow = 103/100` (exact). -/
theorem raceWindow_integral : ∫ x, raceWindow x = 103 / 100 := by
  exact RaceLib.integral_raceWindow

/-- The certificate `RaceWindowCert` (`Γ_0 < 81/500`) is literally the numerical hypothesis of
`grh_chi4_smoothed_race_real` for `f = raceWindow` and `margin = 7369/20000`. -/
theorem raceWindowCert_iff :
    RaceWindowCert ↔ chi4Gamma0 raceWindow < (1 / 2) * (∫ x, raceWindow x) ^ 2 - 7369 / 20000 := by
  exact RaceLib.raceWindowCert_iff

/-- **The smoothed Chebyshev race modulo 4 for a concrete window.**  There are exactly two inputs:

  * `hGRH : GRHFor chi4` — the Riemann hypothesis for the single L-function `L(s,χ₋₄)`;
  * `hcert : RaceWindowCert` — the inequality `chi4Gamma0 raceWindow < 81/500` about one explicit
    integral, certified outside Lean in Arb ball arithmetic (see the docstring of
    `RaceWindowCert` for the script and its sha256).

Conclusion: there is `L₀` such that for every `L ≥ L₀`

    Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L)  <  −7369/20000 = −0.36845,

where `g(x) = ∫ f(t) f(t − x) dt` for the explicit `C²` window `f = raceWindow ≥ 0`: for all large
`L` the primes `≡ 3 mod 4` are ahead of the primes `≡ 1 mod 4` in this smoothed race. -/
theorem grh_chi4_race_concrete_window (hGRH : GRHFor chi4) (hcert : RaceWindowCert) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum raceWindow L < -(7369 / 20000) := by
  exact RaceLib.raceWindow_race hGRH hcert

/-- **Two-sided form for the concrete window** (GRH only; no numerical input).  For every `ε > 0`
and all sufficiently large `L`:  `|chi4RaceSum raceWindow L + 10609/20000| ≤ Γ_0 + ε`,
`Γ_0 = chi4Gamma0 raceWindow`, `10609/20000 = (1/2)(∫ raceWindow)²`. -/
theorem grh_chi4_race_concrete_window_two_sided (hGRH : GRHFor chi4) (ε : ℝ) (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀,
      |chi4RaceSum raceWindow L + 10609 / 20000| ≤ chi4Gamma0 raceWindow + ε := by
  exact RaceLib.raceWindow_two_sided hGRH hε
