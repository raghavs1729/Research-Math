/-
ChebyshevSolution.lean — the UNTRUSTED solution module: the statements of ChebyshevChallenge.lean,
byte-identical, each PROVED by delegating to ChebLib (built on the Zeta23 library,
anthropics/formal-math, folder zeta23/, and on Mathlib).

Part I  (compactly supported windows; proofs in ChebLib/Race, taken over from monotonicity/lean-race):
  explicit_formula_zero_sum, zeroSum_summable_and_archTerm_integrable,
  grh_smoothed_prime_power_sum_bound, grh_smoothed_prime_power_sum_le_zero_sum,
  grhFor_of_dirichlet_zeros_on_line, chi4_isPrimitive, chi4_archWeight,
  grh_chi4_smoothed_prime_power_sum_bound, prime_number_theorem_theta,
  grh_chi4_smoothed_race_two_sided, grh_chi4_smoothed_race, chi4_primeSum_eq_chi4RaceSum,
  chi4_archTerm_zero_eq_chi4Gamma0, chi4Gamma0_integrand_integrable,
  grh_chi4_smoothed_race_real_two_sided, grh_chi4_smoothed_race_real.

Part II (Abel weights e^{-n/x}; proofs in ChebLib/Cheb):
  chi4_apply_natCast, chi4_natCast_two, chi4_natCast_of_odd      →  ChebLib.chi4_*
  grh_explicit_formula_integrable_test                           →  ChebLib.zeroSum_eq_of_integrable
  grh_chi4_gammaZeroSum_summable, grh_chi4_gammaZeroSum_lt       →  ChebLib.summable_gammaZeroSum_chi4,
                                                                    ChebLib.gammaZeroSum_chi4_lt
  chebyshevPsi_summable, chebyshevLogSum_summable,
  chebyshevSum_summable                                          →  ChebLib.summable_*
  grh_chi4_chebyshevPsi_bound                                    →  ChebLib.abs_chebyshevPsi_le
  grh_chi4_chebyshevLogSum_two_sided, grh_chi4_chebyshevLogSum_neg → ChebLib.abs_chebyshevLogSum_add_le,
                                                                    ChebLib.chebyshevLogSum_le
  grh_chi4_smoothed_race_abel_window                             →  ChebLib.window_sum_le
  grh_chebyshev_bias                                             →  ChebLib.tendsto_chebyshevSum
  grh_chebyshev_assertion                                        →  ChebLib.tendsto_chebyshev_literal

The definitions used in the statements are those of ChebLib/Race/Defs.lean and ChebLib/Cheb/Defs.lean
(Mathlib only; inlined verbatim in ChebyshevChallenge.lean).  Nothing in this file is part of the
trusted base.
-/
import ChebLib

open scoped BigOperators ComplexConjugate
open Complex MeasureTheory

/-! ## Part I: compactly supported windows -/

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

/-! ## Part II: Abel weights and Chebyshev's assertion -/

/-- The Dirichlet character `chi4` (`χ₋₄` as a `ℂ`-valued Dirichlet character mod 4) takes at a
natural number `n` the value of Mathlib's quadratic character `ZMod.χ₄`. -/
theorem chi4_apply_natCast (n : ℕ) :
    chi4 (n : ZMod 4) = (((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) : ℂ) := by
  exact ChebLib.chi4_apply_natCast n

/-- `χ₋₄(2) = 0`. -/
theorem chi4_natCast_two : ((ZMod.χ₄ ((2 : ℕ) : ZMod 4) : ℤ) : ℝ) = 0 := by
  exact ChebLib.chi4_natCast_two

/-- `χ₋₄(n) = (−1)^{(n−1)/2}` for odd `n` (Chebyshev's sign: `+1` for `n ≡ 1`, `−1` for `n ≡ 3`
mod 4). -/
theorem chi4_natCast_of_odd (n : ℕ) (hn : n % 2 = 1) :
    ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = (-1) ^ ((n - 1) / 2) := by
  exact ChebLib.chi4_natCast_of_odd hn

/-- **The explicit formula for test functions without compact support, under GRH.**  Let `χ` be a
primitive Dirichlet character modulo `q > 1` and ASSUME the Riemann hypothesis for `L(s,χ)`.  Let
`k : ℝ → ℂ` be twice differentiable, `k' = k1`, `k'' = k2` with `k2` continuous, with `k`, `k'`,
`k''` Lebesgue integrable and `Σ_n Λ(n) n^{-1/2} (|k(log n)| + |k(−log n)|) < ∞`.  Then the zero sum
converges absolutely and

    Σ_ρ m_ρ k̂(γ_ρ) = Γ_0 − Σ_{n ≥ 1} Λ(n) n^{-1/2} (χ(n) k(log n) + conj(χ(n)) k(−log n)),

exactly as in `explicit_formula_zero_sum` (which requires compact support but no hypothesis on the
zeros).  Proof: truncation `k_j = η_j k` and dominated convergence on the three sides; GRH is used
to have real ordinates `γ_ρ`, where `|k̂_j(γ)| ≤ M/(1+γ²)` uniformly in `j`. -/
theorem grh_explicit_formula_integrable_test {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q)
    (hq : 1 < q) (hχ : χ.IsPrimitive) (hGRH : GRHFor χ) (k k1 k2 : ℝ → ℂ)
    (hk1 : ∀ u, HasDerivAt k (k1 u) u) (hk2 : ∀ u, HasDerivAt k1 (k2 u) u) (hk2c : Continuous k2)
    (hki : Integrable k) (hk1i : Integrable k1) (hk2i : Integrable k2)
    (hps : Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n / Real.sqrt n
      * (‖k (Real.log n)‖ + ‖k (-Real.log n)‖))) :
    Summable (fun ρ : {ρ : ℂ | χ.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1} =>
        ((analyticOrderAt χ.LFunction (ρ : ℂ)).toNat : ℂ) * raceFT k (((ρ : ℂ) - 1 / 2) / I))
      ∧ zeroSum χ k = archTerm χ k 0 - primePowerSumSym χ k := by
  exact ChebLib.zeroSum_eq_of_integrable hq hχ hGRH hk1 hk2 hk2c hki hk1i hk2i hps

/-- Under the Riemann hypothesis for `L(s,χ₋₄)`, the sum `Σ_ρ m_ρ |Γ(ρ)|` over the zeros in the
critical strip converges. -/
theorem grh_chi4_gammaZeroSum_summable (hGRH : GRHFor chi4) :
    Summable (fun ρ : {ρ : ℂ | chi4.LFunction ρ = 0 ∧ 0 < ρ.re ∧ ρ.re < 1} =>
      ((analyticOrderAt chi4.LFunction (ρ : ℂ)).toNat : ℝ) * ‖Complex.Gamma (ρ : ℂ)‖) := by
  exact ChebLib.summable_gammaZeroSum_chi4 hGRH

/-- **The zero sum is smaller than the bias.**  ASSUME the Riemann hypothesis for `L(s,χ₋₄)`.  Then

    Σ_ρ m_ρ |Γ(ρ)|  <  √π / 2 = Γ(1/2)/2,

the sum running over all zeros of `L(s,χ₋₄)` in the critical strip, with multiplicity.  (The true
value of the left side is about `4·10⁻⁴`; what is proved is the inequality.)  No numerical
information about the zeros is used: the explicit formula for the test function
`g(u) = (5/2) e^{−(3/2)|u|} − (3/2) e^{−(5/2)|u|}`, whose transform `30/((9/4+r²)(25/4+r²))`
dominates `(32/15) |Γ(1/2+ir)|/√π`, evaluates the positive sum `Σ_ρ m_ρ ĝ(γ_ρ)` in closed form as
`7/2 − γ − 3 log 2 − log π + Σ_n Λ(n) (−χ₋₄(n)) (5n − 3)/n³`, which is shown to be `< 16/15`. -/
theorem grh_chi4_gammaZeroSum_lt (hGRH : GRHFor chi4) :
    gammaZeroSum chi4 < Real.sqrt Real.pi / 2 := by
  exact ChebLib.gammaZeroSum_chi4_lt hGRH

/-- The series defining `chebyshevPsi x` converges absolutely for `x > 0` (unconditional). -/
theorem chebyshevPsi_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun n : ℕ => ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(n : ℝ) / x)) := by
  exact ChebLib.summable_chebyshevPsi hx

/-- The series defining `chebyshevLogSum x` converges absolutely for `x > 0` (unconditional). -/
theorem chebyshevLogSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  exact ChebLib.summable_chebyshevLogSum hx

/-- The series defining Chebyshev's sum `chebyshevSum x` converges absolutely for `x > 0`
(unconditional).  (Lean's `∑'` returns `0` for a divergent series; this excludes that reading of
the statements below.) -/
theorem chebyshevSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes =>
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)) := by
  exact ChebLib.summable_primes_chi4_exp hx

/-- **The Abel-summed von Mangoldt sum under GRH.**  ASSUME the Riemann hypothesis for `L(s,χ₋₄)`.
Then for every `ε > 0` and all sufficiently large `x`:

    | Σ_{n ≥ 1} Λ(n) χ₋₄(n) e^{−n/x} |  ≤  ( Σ_ρ m_ρ |Γ(ρ)| + ε ) √x.

(From the explicit formula `Σ_n Λ(n) χ(n) e^{−n/x} = −Σ_ρ Γ(ρ) x^ρ + …`, obtained here from
`grh_explicit_formula_integrable_test` with `k(u) = e^{u/2} (exp(−e^u/x) − exp(−e^u))`; the
archimedean term is `o(√x)` by the Riemann–Lebesgue lemma.) -/
theorem grh_chi4_chebyshevPsi_bound (hGRH : GRHFor chi4) (ε : ℝ) (hε : 0 < ε) :
    ∀ᶠ x : ℝ in Filter.atTop, |chebyshevPsi x| ≤ (gammaZeroSum chi4 + ε) * Real.sqrt x := by
  exact ChebLib.abs_chebyshevPsi_le hGRH hε

/-- **The logarithmically weighted Chebyshev sum, two-sided.**  ASSUME the Riemann hypothesis for
`L(s,χ₋₄)`.  Then for every `ε > 0` and all sufficiently large `x`:

    | Σ_p χ₋₄(p) (log p) e^{−p/x} + (√π/2) √x |  ≤  ( Σ_ρ m_ρ |Γ(ρ)| + ε ) √x.

The main term `−(√π/2) √x` is the contribution of the squares of primes (`χ₋₄(p²) = 1`), evaluated
with the prime number theorem (`prime_number_theorem_theta`); by `grh_chi4_gammaZeroSum_lt` the
right-hand side is smaller than the main term. -/
theorem grh_chi4_chebyshevLogSum_two_sided (hGRH : GRHFor chi4) (ε : ℝ) (hε : 0 < ε) :
    ∀ᶠ x : ℝ in Filter.atTop, |chebyshevLogSum x + Real.sqrt Real.pi / 2 * Real.sqrt x|
      ≤ (gammaZeroSum chi4 + ε) * Real.sqrt x := by
  exact ChebLib.abs_chebyshevLogSum_add_le hGRH hε

/-- **The logarithmically weighted Chebyshev sum is eventually negative, of order `√x`.**  ASSUME
the Riemann hypothesis for `L(s,χ₋₄)`.  Then there is `c > 0` such that for all sufficiently
large `x`:

    Σ_p χ₋₄(p) (log p) e^{−p/x}  ≤  −c √x.

No hypothesis other than GRH for this one L-function: no linear independence, no simplicity of
zeros, no numerical certificate. -/
theorem grh_chi4_chebyshevLogSum_neg (hGRH : GRHFor chi4) :
    ∃ c : ℝ, 0 < c ∧ ∀ᶠ x : ℝ in Filter.atTop, chebyshevLogSum x ≤ -c * Real.sqrt x := by
  exact ChebLib.chebyshevLogSum_le hGRH

/-- **A smoothed Chebyshev race for primes with a concrete window, no numerical hypothesis.**
ASSUME the Riemann hypothesis for `L(s,χ₋₄)`.  For the window `g(u) = e^{u/2} exp(−e^u)` (positive,
integrable, of rapid decay; not compactly supported) there is `c > 0` such that for all sufficiently
large `L`:

    Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L)  ≤  −c  <  0.

(This is `grh_chi4_chebyshevLogSum_neg` at `x = e^L`, divided by `√x`.  Compare
`grh_chi4_smoothed_race_real`, for compactly supported windows, which needs a numerical hypothesis
on the window.) -/
theorem grh_chi4_smoothed_race_abel_window (hGRH : GRHFor chi4) :
    ∃ c : ℝ, 0 < c ∧ ∀ᶠ L : ℝ in Filter.atTop,
      ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
        * (Real.log (p : ℕ) / Real.sqrt (p : ℕ))
        * (Real.exp ((Real.log (p : ℕ) - L) / 2) * Real.exp (-Real.exp (Real.log (p : ℕ) - L)))
      ≤ -c := by
  exact ChebLib.window_sum_le hGRH

/-- **Chebyshev's assertion (1853) under GRH for `L(s,χ₋₄)`.**  ASSUME that every zero of the
Dirichlet L-function `L(s,χ₋₄)` in the critical strip `0 < Re s < 1` has real part `1/2`
(`hGRH : GRHFor chi4`; an explicit hypothesis, not an axiom).  Then

    Σ_p χ₋₄(p) e^{−p/x}  →  −∞      as x → +∞,

the sum over all primes `p`, with `χ₋₄(2) = 0` and `χ₋₄(p) = (−1)^{(p−1)/2}` for odd `p`
(`chi4_natCast_two`, `chi4_natCast_of_odd`); the series converges absolutely for every `x > 0`
(`chebyshevSum_summable`).  In words: with the weights `e^{−p/x}`, the primes `≡ 3 mod 4` outweigh
the primes `≡ 1 mod 4` by an unbounded amount.  (Hardy–Littlewood and Landau showed that the
assertion is equivalent to the absence of zeros of `L(s,χ₋₄)` with real part `> 1/2`; this is the
implication from GRH.)  The proof uses, besides `hGRH`, only Mathlib, the Zeta23 library's explicit
formula for primitive Dirichlet L-functions, and its port of the prime number theorem. -/
theorem grh_chebyshev_bias (hGRH : GRHFor chi4) :
    Filter.Tendsto chebyshevSum Filter.atTop Filter.atBot := by
  exact ChebLib.tendsto_chebyshevSum hGRH

/-- **Chebyshev's assertion in his own form.**  ASSUME the Riemann hypothesis for `L(s,χ₋₄)`.  Then

    e^{−3c} − e^{−5c} + e^{−7c} + e^{−11c} − e^{−13c} − e^{−17c} + e^{−19c} + ⋯  →  +∞    as c → 0⁺,

the sum over the odd primes `p` with sign `(−1)^{(p+1)/2}` (`+` for `p ≡ 3`, `−` for `p ≡ 1` mod 4).
(`𝓝[>] 0` is the filter of right neighbourhoods of `0`.) -/
theorem grh_chebyshev_assertion (hGRH : GRHFor chi4) :
    Filter.Tendsto (fun c : ℝ => ∑' p : Nat.Primes,
      if (p : ℕ) = 2 then 0 else (-1 : ℝ) ^ (((p : ℕ) + 1) / 2) * Real.exp (-((p : ℕ) : ℝ) * c))
      (nhdsWithin 0 (Set.Ioi 0)) Filter.atTop := by
  exact ChebLib.tendsto_chebyshev_literal hGRH
