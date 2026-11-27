# lean-chebyshev: GRH for L(s, χ₋₄) implies Chebyshev's assertion (Lean 4)

A machine-checked proof, in Lean 4 with Mathlib and the public library `zeta23`
(`anthropics/formal-math`), of the following implication.

> **Theorem (`grh_chebyshev_bias`).** Assume that every zero of the Dirichlet L-function L(s, χ₋₄)
> in the strip 0 < Re s < 1 has real part 1/2. Then
> Σ_p χ₋₄(p) e^{−p/x} → −∞ as x → +∞,
> the sum over all primes, χ₋₄(2) = 0 and χ₋₄(p) = (−1)^{(p−1)/2} for odd p.

This is the assertion Chebyshev made in 1853 (in his own form, `grh_chebyshev_assertion`:
e^{−3c} − e^{−5c} + e^{−7c} + e^{−11c} − e^{−13c} − ⋯ → +∞ as c → 0⁺). Hardy–Littlewood (1917) and
Landau (1918) showed it is equivalent to L(s, χ₋₄) having no zero with real part > 1/2; what is
formalized here is the implication from GRH, not the converse.

The Riemann hypothesis for L(s, χ₋₄) is an **explicit hypothesis** of the theorems (`GRHFor chi4`),
never an axiom. There is **no other hypothesis**: no numerical certificate, no data on zeros, no
linear independence, no simplicity of zeros. The prime number theorem used for the prime squares is
proved (from zeta23's port of `MediumPNT`).

Status: `lake build` succeeds; no `sorry` outside the challenge file; every headline theorem
depends only on `propext`, `Classical.choice`, `Quot.sound`; no `native_decide`; no `axiom`.
Nobody has reviewed this development other than the Lean kernel and the agent that wrote it; the
statements below are short and are what a reader should check.

## The exact Lean statements (headline)

From `ChebyshevChallenge.lean` (which imports only Mathlib):

```lean
def chi4 : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)

def GRHFor {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → χ.LFunction s = 0 → s.re = 1 / 2

def chebyshevSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)

theorem grh_chebyshev_bias (hGRH : GRHFor chi4) :
    Filter.Tendsto chebyshevSum Filter.atTop Filter.atBot

theorem grh_chebyshev_assertion (hGRH : GRHFor chi4) :
    Filter.Tendsto (fun c : ℝ => ∑' p : Nat.Primes,
      if (p : ℕ) = 2 then 0 else (-1 : ℝ) ^ (((p : ℕ) + 1) / 2) * Real.exp (-((p : ℕ) : ℝ) * c))
      (nhdsWithin 0 (Set.Ioi 0)) Filter.atTop

theorem chebyshevSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes =>
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x))

theorem chi4_natCast_two : ((ZMod.χ₄ ((2 : ℕ) : ZMod 4) : ℤ) : ℝ) = 0
theorem chi4_natCast_of_odd (n : ℕ) (hn : n % 2 = 1) :
    ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = (-1) ^ ((n - 1) / 2)
theorem chi4_apply_natCast (n : ℕ) :
    chi4 (n : ZMod 4) = (((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) : ℂ)

theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ
```

Reading notes (what a reader must check to trust that this is Chebyshev's assertion):

- `chi4` is Mathlib's quadratic character `ZMod.χ₄` (values 0, 1, 0, −1 on 0, 1, 2, 3 mod 4) as a
  `ℂ`-valued Dirichlet character; `χ.LFunction` is Mathlib's `DirichletCharacter.LFunction` (the
  analytic continuation of the L-series).
- Lean's `∑'` is `0` for a non-summable family; `chebyshevSum_summable` shows the series is
  absolutely convergent for every x > 0, so `chebyshevSum x` is the genuine sum. (For x ≤ 0 the value
  is irrelevant to a limit at +∞.)
- `Filter.Tendsto f atTop atBot` is "f(x) → −∞ as x → +∞"; `nhdsWithin 0 (Set.Ioi 0)` is c → 0⁺.
- The hypothesis `GRHFor chi4` has exactly the shape that a theorem of the form
  `∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ), 0 < s.re → s.re < 1 →
  DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2` would discharge
  (`grhFor_of_dirichlet_zeros_on_line`). Nothing here depends on, or says anything about, any
  announced proof of such a theorem.

## Other statements proved (all in `ChebyshevChallenge.lean`, 31 in total)

Under `GRHFor chi4` only (targets T2 and intermediate results of T3):

| name | statement |
|---|---|
| `grh_chi4_chebyshevLogSum_two_sided` | for every ε > 0, eventually \|Σ_p χ₋₄(p) (log p) e^{−p/x} + (√π/2)√x\| ≤ (G + ε)√x, with G = Σ_ρ m_ρ \|Γ(ρ)\| (`gammaZeroSum chi4`) |
| `grh_chi4_chebyshevLogSum_neg` | there is c > 0 with Σ_p χ₋₄(p) (log p) e^{−p/x} ≤ −c√x eventually |
| `grh_chi4_smoothed_race_abel_window` | **T2 without numerical hypothesis**: for the concrete window g(u) = e^{u/2} exp(−e^u), Σ_p χ₋₄(p) (log p) p^{−1/2} g(log p − L) ≤ −c < 0 for all large L |
| `grh_chi4_gammaZeroSum_lt`, `grh_chi4_gammaZeroSum_summable` | G = Σ_ρ m_ρ \|Γ(ρ)\| converges and G < √π/2 (multiplicities included; no zero numerics) |
| `grh_chi4_chebyshevPsi_bound` | for every ε > 0, eventually \|Σ_n Λ(n) χ₋₄(n) e^{−n/x}\| ≤ (G + ε)√x |
| `grh_explicit_formula_integrable_test` | for every primitive χ mod q > 1, under `GRHFor χ`: the explicit formula for test functions k with k, k′, k″ integrable (k″ continuous) and absolutely convergent prime side; no compact support |

Unconditional: `chebyshevPsi_summable`, `chebyshevLogSum_summable`, `chebyshevSum_summable`,
`chi4_isPrimitive`, `prime_number_theorem_theta` (θ(x)/x → 1), `explicit_formula_zero_sum`.

Part I of the challenge file (16 statements, **taken over unchanged from `monotonicity/lean-race`**,
whose proof files are copied to `ChebLib/Race/`): the explicit-formula bound for smoothed
prime-power sums of any primitive character under GRH (target **T1**:
`grh_smoothed_prime_power_sum_bound`, `grh_smoothed_prime_power_sum_le_zero_sum`), and the smoothed
race for compactly supported windows f ⋆ f̃ (`grh_chi4_smoothed_race_real_two_sided`, and
`grh_chi4_smoothed_race_real`, which carries the explicit numerical hypothesis `hcert` on the
window). Part II (15 statements) is new here.

## Hypotheses that remain

- `GRHFor chi4` (or `GRHFor χ` in the general statements). That is all for Part II.
- Part I only: `grh_chi4_smoothed_race` and `grh_chi4_smoothed_race_real` keep the numerical
  hypothesis `hcert` for a compactly supported window (as in lean-race). The statement
  `grh_chi4_smoothed_race_abel_window` is the version with a concrete window and no such hypothesis.

## What is not formalized

- The converse (Landau: if the assertion holds then L(s, χ₋₄) has no zero with Re s > 1/2), and
  hence the equivalence.
- Any modulus other than 4; any rate for the divergence of `chebyshevSum` beyond what the proof
  gives internally (it shows `chebyshevSum x ≤ C − c′ x^{1/4}`; the true order is √x / log x).
- The Hardy–Littlewood formula Σ_n Λ(n)χ(n)e^{−n/x} = −Σ_ρ Γ(ρ)x^ρ + ⋯ with its lower-order terms
  made explicit: only the bound `grh_chi4_chebyshevPsi_bound` that follows from it is stated.
- The official Comparator run (`leanprover/comparator` with `lean4export`, `nanoda`, `landrun`) was
  not performed; `comparator.json` is provided, and `scripts/compare_statements.sh` checks that the
  31 statements elaborate to identical types in the challenge and solution modules.
- No human has read the proofs. The trusted base is: the Lean kernel, Mathlib, the definitions and
  statements in `ChebyshevChallenge.lean`.

## Pins and build

| | |
|---|---|
| Lean | `leanprover/lean4:v4.33.0-rc2` (`lean-toolchain`) |
| Mathlib | `51e6992efd06126df61a496bebf8f49482a4e129` |
| zeta23 | `https://github.com/anthropics/formal-math`, commit `11b6c7e3c927df2658f513b6b68ea7f8929c62c5`, `subDir = "zeta23"` (Apache-2.0) |

```
cd monotonicity/lean-chebyshev
lake exe cache get                                   # Mathlib oleans; never build Mathlib from source
lake build ChebyshevChallenge ChebyshevSolution      # builds the needed zeta23 modules, then this project
lake env lean scripts/PrintAxioms.lean               # axiom audit, expected output in logs/print-axioms.txt
./scripts/compare_statements.sh                      # challenge vs solution statements
python3 scripts/gen_challenge.py --check             # the challenge file is the generated one
```

Do not run `lake clean` (it removes the Mathlib cache). On 4 cores the zeta23 modules that are
imported (`Zeta23.ThmE.MainChi`, `Zeta23.ThmE.GammaFactsChiProof`, `Zeta23.FromPNTPlus.MediumPNT`
and their dependencies, 110 modules) build in about 9 minutes; this project's own modules in a few
minutes more (`logs/build.txt`).

Module names: the Palomar template calls the statement modules `Challenge` / `Solution`. Here they
are `ChebyshevChallenge` / `ChebyshevSolution`, because the zeta23 package already ships root
modules named `Challenge`, `ChallengeDeps`, `Solution`, and Lake refuses to disambiguate imports
between two packages.

`#print axioms` (full output in `logs/print-axioms.txt`; all 31 lines have the same axiom list):

```
'grh_chi4_gammaZeroSum_lt' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_chebyshevPsi_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_chebyshevLogSum_two_sided' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_chebyshevLogSum_neg' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_smoothed_race_abel_window' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chebyshev_bias' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chebyshev_assertion' depends on axioms: [propext, Classical.choice, Quot.sound]
```

## The proof, and where each step is

Write ψ(x) = Σ_n Λ(n)χ₋₄(n)e^{−n/x}, Φ(x) = Σ_p χ₋₄(p)(log p)e^{−p/x}, C(x) = Σ_p χ₋₄(p)e^{−p/x}.

1. **Explicit formula without compact support** (`ChebLib/Cheb/EFExt.lean`, `EF_ext`).
   zeta23 proves Weil's explicit formula for primitive L(s, χ) and test functions in C²_c
   (`Zeta23.ThmE.EF_lit_chi_L`). Under GRH the ordinates γ_ρ are real, so for the truncations
   k_j = η_j k one has |k̂_j(γ)| ≤ M/(1+γ²) uniformly; dominated convergence on the zero side
   (Σ_ρ m_ρ/(1+γ_ρ²) < ∞ from zeta23), the archimedean side and the prime side extends the formula
   to every k with k, k′, k″ integrable.
2. **The Abel kernel** (`Kernel.lean`). h_x(u) = e^{u/2}exp(−e^u/x) has h_x(log n)/√n = e^{−n/x} and
   transform x^{1/2+ir}Γ(1/2+ir); |Γ(1/2+ir)|² = π/cosh(πr).
3. **ψ(x)** (`PsiAbel.lean`). The explicit formula for k = h_x − h_1 (the difference makes the dual
   prime side absolutely convergent) gives
   ψ(x) = ψ(1) − D(x) − Σ_ρ m_ρ(x^{1/2+iγ_ρ} − 1)Γ(1/2+iγ_ρ) + (1/2π)∫(x^{1/2+ir} − 1)Γ(1/2+ir)W(r)dr,
   with D bounded and W the archimedean weight; the last integral is o(√x) by the Riemann–Lebesgue
   lemma. Hence |ψ(x)| ≤ (G + ε)√x, G = Σ_ρ m_ρ|Γ(1/2+iγ_ρ)|.
4. **G < √π/2 without zero numerics** (`Majorant.lean`, `Poisson.lean`, `Numeric.lean`,
   `ZeroSum.lean`). The C² function g(u) = (5/2)e^{−(3/2)|u|} − (3/2)e^{−(5/2)|u|} has transform
   ĝ(r) = 30/((9/4+r²)(25/4+r²)) > 0, and |Γ(1/2+ir)| ≤ (15/32)√π ĝ(r) (from cosh(πr) ≥ a degree-8
   Taylor polynomial). The explicit formula for g is a closed form:
   Σ_ρ m_ρ ĝ(γ_ρ) = (7/2 − γ − 3 log 2 − log π) + Σ_n Λ(n)(−χ₋₄(n))(5n − 3)/n³
   (the archimedean side is a Poisson integral of Re ψ, computed from the partial-fraction series of
   the digamma function: (5/2)ψ(3/2) − (3/2)ψ(2) + log(4/π)). The right side is about 0.0462
   (floating point); `maj_numeric` proves it is < 16/15, which is what is needed, with three explicit
   terms and a telescoping bound for the tail.
5. **Prime powers** (`AbelSum.lean`, `ChebLib/Race/NonPrime.lean`, `ChebLib/Race/PNT.lean`).
   ψ(x) − Φ(x) = Σ_n c(n)e^{−n/x} with c supported on proper prime powers; its partial sums are
   ~ √y by the prime number theorem (prime squares, χ₋₄(p²) = 1), and an Abelian theorem gives
   Σ_n c(n)e^{−n/x} ~ Γ(3/2)√x = (√π/2)√x. So Φ(x) = −(√π/2)√x + O((G+ε)√x) ≤ −c√x (`Final.lean`).
6. **Removing log p** (`LogRemoval.lean`). From e^{−n/x} = (log n)∫₀^∞ n^{−s}e^{−n/x}ds and
   n^{−s} = Γ(s)^{−1}∫₀^∞ v^{s−1}e^{−nv}dv, the sum C(x) is an average of Φ(1/(v+1/x)) against a
   positive kernel, so the one-sided bound Φ(y) ≤ −c√y for large y (and |Φ(y)| ≤ C₀e^{−2/y} for
   small y) gives C(x) → −∞.

File map: `ChebLib/Race/` (from lean-race: definitions, explicit-formula bound, compact windows,
PNT), `ChebLib/Cheb/` (new: steps 1–6, `Defs.lean` for the Mathlib-only definitions, `Headline.lean`
for the statements in final form), `ChebyshevSolution.lean` (statements + one-line delegations),
`ChebyshevChallenge.lean` (generated: definitions inlined, proofs replaced by `sorry`),
`comparator.json`, `scripts/`, `logs/`, `PROGRESS.md`.

Sizes: `ChebLib/Cheb/` is about 3,500 lines, `ChebLib/Race/` about 1,400.

## Numerical orientation (floating point, not part of the proof)

Zeros of L(s, χ₋₄): first ordinates 6.0209, 10.2438, 12.9881, …; G = Σ_ρ |Γ(ρ)| ≈ 3.9·10⁻⁴ against
√π/2 ≈ 0.886; Σ_ρ ĝ(γ_ρ) ≈ 0.0462 against the proved bound 16/15.

## Provenance

Written by an AI agent (Claude Code) in one session on 2026-11-24, with sub-agents proving
individual lemma files against fixed statements; all files were compiled together and the axiom
audit was run on the final build. The proofs in `ChebLib/Race/` come from `monotonicity/lean-race`
in this repository. zeta23 is Apache-2.0; its `MediumPNT` is a port from the
PrimeNumberTheoremAnd project (see zeta23's `NOTICE`).
