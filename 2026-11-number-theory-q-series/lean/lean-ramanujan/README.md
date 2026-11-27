# lean-ramanujan: Ramanujan's prime-counting inequality in Lean 4

Status (2026-11-24): **all six statements of section 1 compile with no `sorry`, no `axiom`, no
`native_decide`; `#print axioms` gives `[propext, Classical.choice, Quot.sound]` for each.**
Four are unconditional (the headline among them). Two are conditional, with their hypotheses as
explicit binders. The prime number theorem with error term is proved, not assumed (zeta23's
`MediumPNT`).

## 1. What is proved

`π(x)` is the number of primes `≤ x`; in Lean, `Nat.primeCounting ⌊x⌋₊` for real `x`.

### Headline: Ramanujan's inequality for all sufficiently large real x

```lean
theorem ramanujan_prime_counting_inequality :
    ∃ x₀ : ℝ, ∀ x : ℝ, x₀ ≤ x →
      (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        < Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ)
```

In words: there is `x₀` such that `π(x)² < (e x / log x) · π(x / e)` for every real `x ≥ x₀`.
No hypotheses.

### The quantitative form (Ramanujan's computation)

```lean
theorem ramanujan_gap_asymptotic :
    Filter.Tendsto (fun x : ℝ =>
      (Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ)
        - (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2) / (x ^ 2 / Real.log x ^ 6))
      Filter.atTop (nhds 1)
```

`(e x / log x) π(x/e) − π(x)² = (1 + o(1)) x² / log⁶ x`. This implies the headline and gives the
exact order of the margin.

### The two inputs, as statements of their own

```lean
theorem primeCounting_five_term_expansion :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ)
          - x * (1 / Real.log x + 1 / Real.log x ^ 2 + 2 / Real.log x ^ 3 + 6 / Real.log x ^ 4
            + 24 / Real.log x ^ 5)|
        ≤ C * (x / Real.log x ^ 6)

theorem chebyshev_theta_sub_le :
    ∃ K : ℝ, ∀ x : ℝ, 2 ≤ x → |Chebyshev.theta x - x| ≤ K * (x / Real.log x ^ 5)
```

### Conditional statements (stretch S2; hypotheses are binders, see section 6)

```lean
theorem ramanujan_of_schoenfeld_bound (c : ℝ) (hc : |c| ≤ 2)
    (hS : ∀ x : ℝ, 2657 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..x, 1 / Real.log t)|
        < √x * Real.log x / (8 * Real.pi)) :
    ∀ x : ℝ, Real.exp 40 ≤ x →
      (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
        < Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ)

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
        < Real.exp 1 * x / Real.log x * (Nat.primeCounting ⌊x / Real.exp 1⌋₊ : ℝ)
```

## 2. Relation to Ramanujan's claim

Ramanujan (second notebook; Berndt, *Ramanujan's Notebooks* Part IV, ch. 24) derived from
`π(x) = x Σ_{k=0}^{4} k!/log^{k+1} x + O(x/log⁶ x)` that
`π(x)² − (e x/log x) π(x/e) = −x²/log⁶ x + …` and concluded that `π(x)² < (e x/log x) π(x/e)` for
`x` sufficiently large. `ramanujan_prime_counting_inequality` is that conclusion and
`ramanujan_gap_asymptotic` is the displayed computation with the lower-order terms as `o(1)`.

- **Real `x`.** The statement is for every real `x ≥ x₀`, not only integers.
- **`π`.** `Nat.primeCounting n` counts primes `≤ n`; primes `≤ x` are the primes `≤ ⌊x⌋` for
  `x ≥ 0`. `⌊x⌋₊` is the natural floor, `0` for negative `x`; only large `x` matter.
- **Parsing.** `Real.exp 1 * x / Real.log x * π` is `((e·x)/log x)·π`.
- **Not explicit.** `x₀` is existential. A threshold is necessary: the inequality fails at
  `x = 38 358 837 682` (the largest integer failure under RH, Dudek–Platt 2015). See section 5.

## 3. Proof

The unconditional part is about 700 lines in `RamanujanLib/` (`LogIntegral`, `ThetaBound`,
`PiExpansion`, `Main`) on top of Mathlib and one zeta23 theorem.

1. `ThetaBound.lean`. zeta23's `MediumPNT`: `ψ(x) − x = O(x exp(−c (log x)^{1/10}))`. Since
   `y⁵⁰ e^{−y} → 0`, `exp(−c (log x)^{1/10}) ≤ c⁻⁵⁰/log⁵ x` eventually; with Mathlib's
   `0 ≤ ψ − θ ≤ C√x` this gives `|θ(x) − x| ≤ K x/log⁵ x` eventually, and on the remaining
   compact range `θ(x) ≤ x log 4` suffices. Result: the bound for every `x ≥ 2`.
2. `LogIntegral.lean` (Mathlib only). `J k x = ∫_2^x dt/log^k t`. Integration by parts
   `J k x = x/log^k x − 2/log^k 2 + k·J (k+1) x` from the derivative of `t/log^k t`;
   `J k x ≤ C_k x/log^k x` on `[2,∞)` (beyond `T = max(2, e^{2k})` the integrand is at most twice
   that derivative); hence `x/log x + J 2 x = x·P(1/log x) + O(x/log⁶ x)` with
   `P(u) = u + u² + 2u³ + 6u⁴ + 24u⁵`.
3. `PiExpansion.lean`. Mathlib's partial summation
   `π(x) = θ(x)/log x + ∫_2^x θ(t)/(t log² t) dt`, with `θ(t) = t + (θ(t) − t)`: the main part is
   step 2, the error is at most `K x/log⁶ x + K·J 7 x`.
4. `Main.lean`. Put `u = 1/log x`; then `1/log(x/e) = u/(1−u)`. With the expansion constant `C`,
   `π(x) ≤ x(P(u) + C u⁶)` and `π(x/e) ≥ (x/e)(P(u/(1−u)) − C (u/(1−u))⁶)`, so it suffices that
   `(P(u) + C u⁶)² < u (P(u/(1−u)) − C (u/(1−u))⁶)`. The identity (`key_identity`, by
   `field_simp; ring`)

   `u (P(u/(1−u)) − C (u/(1−u))⁶) − (P(u) + C u⁶)² = u⁶ · G C u`,
   `G C u = Q(u)/(1−u)⁵ − C u/(1−u)⁶ − 2 C P(u) − C² u⁶`,
   `Q(u) = 1 + 129u − 275u² + 217u³ − 40u⁴ + 1024u⁵ − 3468u⁶ + 4452u⁷ − 2592u⁸ + 576u⁹`

   reduces everything to `G C` being continuous at `0` with `G C 0 = 1 > 0`. The `1` is Ramanujan's
   `65 − 64`. For the asymptotic, the same identity with `−C` gives the upper bound, and
   `G C (1/log x) ≤ gap/(x²/log⁶ x) ≤ G (−C) (1/log x)` squeezes to `1`.

Nothing was imported or copied from PrimeNumberTheoremAnd's `Ramanujan.lean`; its `criterion`
(Dudek–Platt Lemma 2.1) is replaced by the identity above.

## 4. Layout

| file | role |
|---|---|
| `RamanujanStatement.lean` | trusted statements (6): `import Mathlib` only, no definitions, proofs `sorry` |
| `RamanujanSolution.lean` | the same statements verbatim, each proved by one delegation to `RamanujanLib` |
| `RamanujanLib/LogIntegral.lean` | `∫_2^x dt/log^k t`: integration by parts, bounds, expansion (Mathlib only) |
| `RamanujanLib/ThetaBound.lean` | `|θ(x) − x| ≤ K x/log⁵ x` from `MediumPNT` |
| `RamanujanLib/PiExpansion.lean` | five-term expansion of `π(x)` |
| `RamanujanLib/Main.lean` | the identity, the inequality, the asymptotic |
| `RamanujanLib/Explicit.lean` | explicit criterion: a `√x log x` bound for `π − li` gives the inequality for `x ≥ e⁴⁰` (Mathlib only) |
| `RamanujanLib/Conditional.lean` | Schoenfeld form, integers-to-reals step, sharp reduction (Mathlib only) |
| `scripts/PrintAxioms.lean` | axiom audit |
| `scripts/compare_statements.sh` | statements of the two modules agree textually and after elaboration |
| `comparator.json`, `scripts/verify-comparator.sh` | `leanprover/comparator` run with NanoDa kernel replay |
| `PROGRESS.md` | log |

## 5. What is not formalized

- **An explicit `x₀` (stretch S1).** Not available from this route: `MediumPNT` is stated as
  `∃ c > 0, (ψ − id) =O[atTop] fun x ↦ x * exp (−c * log x ^ (1/10))`. Both `c` and the `O`
  constant and threshold are existential in the statement, so nothing downstream can extract
  values from it. An explicit `x₀` needs an
  explicit PNT error term, which is not in zeta23 or Mathlib at these pins. For comparison
  (`../RAMANUJAN-TRIAGE.md`, Topic 1): the best refereed unconditional tail is `x ≥ exp(3158.442)`
  (Axler, arXiv:2207.02486), from explicit zero-free regions.
- **Schoenfeld's bound under RH**, i.e. `RiemannHypothesis → |π(x) − li x| < √x log x/(8π)` for
  `x ≥ 2657`. It is a hypothesis (`hSchoenfeld`) of the conditional statement. Proving it needs a bound for
  `ψ(x) − x` under RH with explicit constants. zeta23 has Weil's explicit formula and the
  Riemann–von Mangoldt formula, but I found no such explicit bound in it or in Mathlib, and
  deriving one was not attempted.
- **The finite verification** (`hFinite`): no computation of `π` up to `e⁴⁰` is done or imported.
- **The principal-value integral `li`**: not in Mathlib; the conditional statements use
  `c + ∫_2^x dt/log t` with a free constant `|c| ≤ 2` instead (`c = li 2 = 1.045…` gives `li`).

## 6. Under the Riemann Hypothesis (stretch S2)

What the literature says (`../RAMANUJAN-TRIAGE.md`, Topic 1): under RH the inequality holds for
all `x ≥ 38 358 837 683` and fails at `38 358 837 682` (Dudek–Platt 2015).

What is proved here is a **reduction**, in two steps.

1. `ramanujan_of_schoenfeld_bound` (fully proved, no RH, no PNT): a bound of Schoenfeld's shape,
   `|π(x) − (c + ∫_2^x dt/log t)| < √x log x/(8π)` for `x ≥ 2657`, implies the inequality for every
   real `x ≥ e⁴⁰ ≈ 2.35·10¹⁷`. The library version `RamanujanLib.ramanujan_of_sqrt_log_bound` takes
   any constant `A ≤ 1/25` in place of `1/(8π)` and needs the bound only for `y ≥ e³⁹`.
   For comparison, Hassani's RH threshold is `138 766 146 692 471 228 ≈ 1.39·10¹⁷`; a back-of-envelope
   calculation (not formalized) says no argument from Schoenfeld's bound alone can go much below
   that, so `e⁴⁰` is within a factor `1.7` of what this input allows.
2. `ramanujan_sharp_under_RH_reduction`: `RiemannHypothesis`, plus `hSchoenfeld`
   (RH ⇒ Schoenfeld's bound), plus `hFinite` (the inequality at the integers of
   `[38 358 837 683, e⁴⁰)`), imply the inequality for every real `x ≥ 38 358 837 683`. The step from
   integers to reals is proved (`ramanujan_real_of_nat`: `t/log t` and `π(t/e)` are nondecreasing
   for `t ≥ e`, and `π(x) = π(⌊x⌋)`).

So the sharp RH statement is reduced to two named literature inputs. `hSchoenfeld` is Schoenfeld
1976, Corollary 1. `hFinite` is contained in Axler's unconditional result for
`38 358 837 683 ≤ x ≤ 10¹⁹` (arXiv:2207.02486, as recorded in the triage; not checked here);
Dudek–Platt's own computation reaches `1.15·10¹⁶`.

How the criterion is proved (`Explicit.lean`, 480 lines, every constant explicit):

- From the base point `a = e¹⁶`, eight integrations by parts:
  `∫_a^x dt/log t = x P₈(u) − a P₈(1/16) + 8!·∫_a^x dt/log⁹ t`, `u = 1/log x`,
  `P₈(u) = Σ_{k=1}^{8} (k−1)! u^k`, and `0 ≤ ∫_a^x dt/log⁹ t ≤ (16/7) x u⁹` because `log t ≥ 16`.
  With `∫_2^a dt/log t ≤ 1.45 a` and `e¹⁶ ≤ 8.9·10⁶`:
  `x P₈(u) − 1.3·10⁷ ≤ c + ∫_2^x dt/log t ≤ x S₈(u) + 1.3·10⁷`, `S₈ = P₈ + 92160 u⁹`.
- `u P₈(u/(1−u)) − S₈(u)² − u⁶(1 + 14u) = u⁸ M₈(u)/(1−u)⁸` with an explicit polynomial `M₈` of
  degree 18, positive on `[0, 1/40]` (`field_simp; ring` and a term-by-term bound).
- `e^{L/2}(L + 14) ≥ 0.159 L⁷` for `L ≥ 40` (from `e ≥ 2.7182818283` and `1 + t ≤ e^t`), so the gap
  is at least `0.159·x^{3/2}`, against a total error of at most `0.152·x^{3/2}`
  (`explicit_core`: five error terms, bounded one by one).

## 7. Pins and build

- Lean `leanprover/lean4:v4.33.0-rc2` (`lean-toolchain`)
- Mathlib `51e6992efd06126df61a496bebf8f49482a4e129`
- zeta23: `https://github.com/anthropics/formal-math` at
  `11b6c7e3c927df2658f513b6b68ea7f8929c62c5`, `subDir = "zeta23"` (Lake git require). Only
  `Zeta23.FromPNTPlus.MediumPNT` and its imports are used.

```bash
export PATH="$HOME/.elan/bin:$PATH"
cd lean/lean-ramanujan
lake exe cache get     # Mathlib oleans; do not build Mathlib from source
lake build             # RamanujanLib, RamanujanStatement, RamanujanSolution
lake env lean scripts/PrintAxioms.lean
./scripts/compare_statements.sh
```

**Do not run `lake clean`** (it deletes the Mathlib cache).

Measured in a fresh container (4 cores, 16 GB): Mathlib cache download 1 min 38 s; first build of
`Zeta23.FromPNTPlus.MediumPNT` and its zeta23 imports a few minutes; this project alone 35 s.

`lake build` ends with

```
warning: RamanujanStatement.lean:25:8: declaration uses `sorry`
warning: RamanujanStatement.lean:38:8: declaration uses `sorry`
warning: RamanujanStatement.lean:51:8: declaration uses `sorry`
warning: RamanujanStatement.lean:62:8: declaration uses `sorry`
warning: RamanujanStatement.lean:80:8: declaration uses `sorry`
warning: RamanujanStatement.lean:107:8: declaration uses `sorry`
✔ [8717/8723] Built RamanujanLib.Conditional (7.1s)
✔ [8718/8723] Built RamanujanLib.PiExpansion (7.9s)
✔ [8719/8723] Built RamanujanLib.Main (11s)
✔ [8720/8723] Built RamanujanLib (4.0s)
✔ [8722/8723] Built RamanujanSolution (5.4s)
Build completed successfully (8723 jobs).
```

The six `sorry` warnings are the statement file; `RamanujanLib/` and `RamanujanSolution.lean`
contain none. The other warnings in the log come from zeta23 (deprecations).

`lake env lean scripts/PrintAxioms.lean`:

```
'ramanujan_prime_counting_inequality' depends on axioms: [propext, Classical.choice, Quot.sound]
'ramanujan_gap_asymptotic' depends on axioms: [propext, Classical.choice, Quot.sound]
'primeCounting_five_term_expansion' depends on axioms: [propext, Classical.choice, Quot.sound]
'chebyshev_theta_sub_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'ramanujan_of_schoenfeld_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'ramanujan_sharp_under_RH_reduction' depends on axioms: [propext, Classical.choice, Quot.sound]
```

`./scripts/compare_statements.sh`:

```
OK (textual): 6 statements of RamanujanStatement.lean occur verbatim in RamanujanSolution.lean
OK (elaborated): the 6 statements have identical types in Statement and Solution
```

### Comparator and independent kernel

`./scripts/verify-comparator.sh` runs `leanprover/comparator` on `comparator.json`
(`RamanujanStatement` against `RamanujanSolution`, six theorem names, permitted axioms `propext`,
`Quot.sound`, `Classical.choice`, NanoDa replay enabled), with the tool revisions zeta23 pins for
this toolchain (comparator `75c730e`, lean4export `9fb131b`, nanoda_lib `68d5ca9`, landrun
`811cfff`). It needs `git`, `go`, `cargo` and `lake`. Result here (5 min 2 s):

```
Running nanoda kernel on solution
Nanoda kernel accepts the solution
Running Lean default kernel on solution.
Lean default kernel accepts the solution
Your solution is okay!
```

Two caveats. Comparator calls Landrun with `--best-effort`, so I cannot say whether the Landlock
sandbox was actually enforced in this container. And the script points comparator at the Landrun
binary directly rather than through zeta23's `landrun-wrapper.sh`, which rejects the `--`
delimiter that this comparator revision already passes.

## 8. Licence and attribution

zeta23 is Apache-2.0 (Anthropic, PBC); its `FromPNTPlus/` modules are ported from the
PrimeNumberTheoremAnd project of Kontorovich, Tao et al. This project depends on zeta23 through
Lake and copies no code from it or from PrimeNumberTheoremAnd.
