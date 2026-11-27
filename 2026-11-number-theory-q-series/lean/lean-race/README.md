# lean-race: a smoothed Chebyshev prime race modulo 4 under an explicit GRH hypothesis (Lean 4)

Status (2026-11-24): **everything below compiles, with no `sorry`, no `axiom`, no `native_decide`,
and `#print axioms` = `[propext, Classical.choice, Quot.sound]` for all 22 headline theorems.**
Stage S0 and stage S1 of `../advances/FORMAL-MATH-REPO.md` are both done, the prime number theorem
needed in S1 is proved rather than assumed, and the race theorem is instantiated at a concrete
window. **The final theorem has exactly two inputs: GRH for `L(s,χ₋₄)` and one externally certified
numerical inequality with a published hash.**

```lean
theorem grh_chi4_race_concrete_window (hGRH : GRHFor chi4) (hcert : RaceWindowCert) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum raceWindow L < -(7369 / 20000)
```

The project is a small layer (about 1,560 lines in `RaceLib/`) on top of the public library
`anthropics/formal-math`, folder `zeta23/`, whose literature-form explicit formula for a primitive
Dirichlet L-function (`Zeta23.ThmE.EF_lit_chi_L`) carries the analytic weight.

## 1. What is proved

Notation. `χ₋₄` is the nontrivial character mod 4 (`chi4`). For a window `f` put
`g = f ⋆ f̃`, `g(x) = ∫ f(t) conj(f(t − x)) dt` (`autocorr f`), so `ĝ = |f̂|² ≥ 0` on the real line,
with `ĝ(z) = ∫ g(u) e^{izu} du` (`raceFT`). `Λ` is von Mangoldt, `ψ = Γ'/Γ` is `Complex.digamma`.

- `S(L) = Σ_n Λ(n) χ(n) n^{-1/2} g(log n − L)` (`primePowerSum χ g L`)
- `P(L) = Σ_p χ(p) (log p) p^{-1/2} g(log p − L)` (`primeSum χ g L`; primes only)
- `Γ_L = (1/2π) ∫ e^{irL} ĝ(r) [Re ψ(1/4 + κ/2 + ir/2) + log(q/π)] dr` (`archTerm χ g L`; `κ` = parity of `χ`)
- `Σ_ρ m_ρ ĝ(γ_ρ)` over the distinct zeros `ρ = 1/2 + iγ_ρ` with `0 < Re ρ < 1`, with multiplicity (`zeroSum χ g`)
- `GRHFor χ := ∀ s : ℂ, 0 < s.re → s.re < 1 → χ.LFunction s = 0 → s.re = 1 / 2`

### S1 for a general real window

```lean
theorem grh_chi4_smoothed_race_real (hGRH : GRHFor chi4) (f : ℝ → ℝ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (margin : ℝ)
    (hcert : chi4Gamma0 f < (1 / 2) * (∫ x, f x) ^ 2 - margin) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum f L < -margin
```

with (from `RaceChallengeDeps.lean`, Mathlib only)

```lean
def chi4RaceSum (f : ℝ → ℝ) (L : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
    * (Real.log (p : ℕ) / Real.sqrt (p : ℕ)) * ∫ t : ℝ, f t * f (t - (Real.log (p : ℕ) - L))

def chi4Gamma0 (f : ℝ → ℝ) : ℝ :=
  (1 / (2 * Real.pi)) * ∫ r : ℝ, ‖raceFT (fun x => (f x : ℂ)) r‖ ^ 2
    * ((Complex.digamma (3 / 4 + I * r / 2)).re + Real.log (4 / Real.pi))
```

In words: assume the Riemann hypothesis for the single L-function `L(s,χ₋₄)`. Take a `C²` window
`f` supported in `[−a/2, a/2]`, `a < log 3`, for which the explicit number
`Γ_0 = (1/2π) ∫ |f̂(r)|² [Re ψ(3/4 + ir/2) + log(4/π)] dr` is smaller than `(1/2)(∫f)² − margin`.
Then for **every** sufficiently large `L`,
`Σ_p χ₋₄(p) (log p) p^{-1/2} g(log p − L) < −margin`. With `margin ≥ 0` and `f ≥ 0` this says the
primes `≡ 3 mod 4` are ahead in the smoothed race for all large `L`: no exceptional set, no
linear-independence hypothesis, no simplicity hypothesis, multiplicities of zeros included.

### All 22 statements (`RaceChallenge.lean` / `RaceSolution.lean`)

| name | content | hypotheses beyond the window |
|---|---|---|
| `explicit_formula_zero_sum` | `zeroSum χ g = archTerm χ g 0 − primePowerSumSym χ g` for primitive `χ` mod `q > 1`, `g ∈ C²_c` | none |
| `zeroSum_summable_and_archTerm_integrable` | the zero sum converges absolutely; the archimedean integrand is integrable | none |
| `grh_smoothed_prime_power_sum_bound` | **S0, general `χ`**: `‖S(L) − Γ_L‖ ≤ Re(Γ_0 − Σ_n Λ(n) n^{-1/2}(χ(n) g(log n) + conj χ(n) g(−log n)))` for `L > A` | `GRHFor χ` |
| `grh_smoothed_prime_power_sum_le_zero_sum` | **S0, zero-sum form**: `‖S(L) − Γ_L‖ ≤ Re Σ_ρ m_ρ ĝ(γ_ρ)` | `GRHFor χ` |
| `grhFor_of_dirichlet_zeros_on_line` | `GRHFor χ` follows from the all-moduli statement with the same binders | none |
| `chi4_isPrimitive` | `chi4.IsPrimitive` | none |
| `chi4_archWeight` | the weight for `χ₋₄` is `Re ψ(3/4 + ir/2) + log(4/π)` | none |
| `grh_chi4_smoothed_prime_power_sum_bound` | **S0 for `χ₋₄`, closed form**: `‖S(L) − Γ_L‖ ≤ Re Γ_0` when `a < log 3`, `L > a` | `GRHFor chi4` |
| `prime_number_theorem_theta` | `θ(x)/x → 1` (`PNTTheta`) | none |
| `grh_chi4_smoothed_race_two_sided` | **S1**: for all `ε > 0`, eventually `‖P(L) + ĝ(0)/2‖ ≤ Re Γ_0 + ε` | `GRHFor chi4` |
| `grh_chi4_smoothed_race` | **S1**: eventually `Re P(L) < −margin` | `GRHFor chi4`, `hcert` |
| `chi4_primeSum_eq_chi4RaceSum` | for real `f`, `P(L)` is the real number `chi4RaceSum f L` | none |
| `chi4_archTerm_zero_eq_chi4Gamma0` | for real `f`, `Γ_0` is the real number `chi4Gamma0 f` | none |
| `chi4Gamma0_integrand_integrable` | the integrand of `chi4Gamma0 f` is integrable | none |
| `grh_chi4_smoothed_race_real_two_sided` | **S1, real**: eventually `|chi4RaceSum f L + (∫f)²/2| ≤ chi4Gamma0 f + ε` | `GRHFor chi4` |
| `grh_chi4_smoothed_race_real` | **S1, real, general window** (above) | `GRHFor chi4`, `hcert` |
| `raceWindow_piecewise` | `raceWindow` is even, `1` on `[−49/100, 49/100]`, the quintic ramp on `[49/100, 27/50]`, `0` beyond | none |
| `raceWindow_admissible` | `raceWindow` is `C²`, compactly supported in `[−27/50, 27/50]`, and `27/25 < log 3` | none |
| `raceWindow_integral` | `∫ raceWindow = 103/100` | none |
| `raceWindowCert_iff` | `RaceWindowCert ↔ chi4Gamma0 raceWindow < (1/2)(∫ raceWindow)² − 7369/20000` | none |
| `grh_chi4_race_concrete_window` | **the final theorem**: eventually `chi4RaceSum raceWindow L < −7369/20000` | `GRHFor chi4`, `RaceWindowCert` |
| `grh_chi4_race_concrete_window_two_sided` | eventually `|chi4RaceSum raceWindow L + 10609/20000| ≤ chi4Gamma0 raceWindow + ε` | `GRHFor chi4` |

S0 verbatim:

```lean
theorem grh_smoothed_prime_power_sum_bound {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q)
    (hq : 1 < q) (hχ : χ.IsPrimitive) (hGRH : GRHFor χ) (f : ℝ → ℂ) (hf : ContDiff ℝ 2 f)
    (hfc : HasCompactSupport f) (A : ℝ) (hA : tsupport f ⊆ Set.Icc (-(A / 2)) (A / 2)) (L : ℝ)
    (hL : A < L) :
    ‖primePowerSum χ (autocorr f) L - archTerm χ (autocorr f) L‖
      ≤ (archTerm χ (autocorr f) 0 - primePowerSumSym χ (autocorr f)).re

theorem grh_chi4_smoothed_prime_power_sum_bound (hGRH : GRHFor chi4) (f : ℝ → ℂ)
    (hf : ContDiff ℝ 2 f) (hfc : HasCompactSupport f) (a : ℝ) (ha : a < Real.log 3)
    (hsupp : tsupport f ⊆ Set.Icc (-(a / 2)) (a / 2)) (L : ℝ) (hL : a < L) :
    ‖primePowerSum chi4 (autocorr f) L - archTerm chi4 (autocorr f) L‖
      ≤ (archTerm chi4 (autocorr f) 0).re
```

## 2. The concrete window, its certificate, and the inputs that remain

### The window

```lean
def smoothstep (t : ℝ) : ℝ :=
  if t ≤ 0 then 0 else if t ≤ 1 then 10 * t ^ 3 - 15 * t ^ 4 + 6 * t ^ 5 else 1

def raceWindow (u : ℝ) : ℝ := smoothstep (20 * (27 / 50 - u)) * smoothstep (20 * (27 / 50 + u))

def RaceWindowCert : Prop := chi4Gamma0 raceWindow < 81 / 500
```

`raceWindow` is an even plateau: `1` on `[−49/100, 49/100]`, quintic smoothstep ramps of width
`1/20`, `0` outside `[−27/50, 27/50]`. Proved in Lean (`RaceLib/Window.lean`): this piecewise
description, `ContDiff ℝ 2 raceWindow`, compact support, `tsupport ⊆ [−a/2, a/2]` with `a = 27/25`,
`27/25 < log 3`, and `∫ raceWindow = 103/100`.

```lean
theorem raceWindow_admissible :
    ContDiff ℝ 2 raceWindow ∧ HasCompactSupport raceWindow
      ∧ tsupport raceWindow ⊆ Set.Icc (-((27 / 25 : ℝ) / 2)) ((27 / 25) / 2)
      ∧ (27 / 25 : ℝ) < Real.log 3

theorem raceWindow_integral : ∫ x, raceWindow x = 103 / 100

theorem raceWindowCert_iff :
    RaceWindowCert ↔ chi4Gamma0 raceWindow < (1 / 2) * (∫ x, raceWindow x) ^ 2 - 7369 / 20000

theorem grh_chi4_race_concrete_window (hGRH : GRHFor chi4) (hcert : RaceWindowCert) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀, chi4RaceSum raceWindow L < -(7369 / 20000)

theorem grh_chi4_race_concrete_window_two_sided (hGRH : GRHFor chi4) (ε : ℝ) (hε : 0 < ε) :
    ∃ L₀ : ℝ, ∀ L ≥ L₀,
      |chi4RaceSum raceWindow L + 10609 / 20000| ≤ chi4Gamma0 raceWindow + ε
```

### The certificate (outside Lean)

`RaceWindowCert` is a `Prop`-valued hypothesis, in the pattern of zeta23's `EnclOK`: Lean does not
prove it and it is not an axiom. It is the statement `Γ_0 < 81/500` for
`Γ_0 = chi4Gamma0 raceWindow = (1/2π) ∫ |f̂(r)|² [Re ψ(3/4 + ir/2) + log(4/π)] dr`.
`cert/window_cert.py` certifies it in Arb ball arithmetic (python-flint 0.9.0, 300 bits, outward
rounding, no floating point on the certified path):

| quantity | value |
|---|---|
| `(1/π) ∫_0^2000 f̂² W` (rigorous quadrature, `acb.integral`) | `0.16116795988245043 ± 1.67e-18` |
| tail bound `(1/π) ∫_2000^∞` | `≤ 5.561928010e-9` |
| **certified enclosure of `Γ_0`** | **`0.16116795432 < Γ_0 < 0.16116796545`** |
| claimed bound | `81/500 = 0.162` |
| `∫ f` (exact, also proved in Lean) | `103/100` |
| `(1/2)(∫f)²` | `10609/20000 = 0.53045` |
| **margin** `(1/2)(∫f)² − 81/500` | **`7369/20000 = 0.36845`** |

How the script works:
- `f̂` at points with `|r| < 40`: the moment series `Σ (−1)^k m_{2k} r^{2k}/(2k)!` with the moments
  computed exactly in rational arithmetic (151 terms, rigorous tail).
- `f̂` for `Re r > 10`: a closed form from three integrations by parts. The two representations are
  checked against each other on the overlap.
- Wide complex balls (needed by the quadrature's error bounds): mean-value enclosures, using
  `|f̂'| ≤ 2h² e^{h|Im r|}` and `|ψ'(w)| ≤ 1/x² + 1/x` for `Re w ≥ x > 0`.
- Tail: `|f̂(r)| ≤ M3/r³` with `M3` the L¹ norm of the third derivative of `f`, `800·40√3/3`, and
  `Re ψ(z) ≤ log|z| + 4/27` on `Re z = 3/4` from Binet's second formula (DLMF 5.9.15).

The closed form, the decay of `f̂` and the Binet bound are derived in the script's header. They are
paper arguments, not machine-checked; the script spot-checks each of them at sample points.

```
sha256  361b1a54173be69d2efb3931ce11b3dd86db01560206b1be3ebaf6963d6a62d6  cert/window_cert.py
sha256  33487c8ed1389a14c74695e7a7b75f09e808a88529ead0dc420b45c845763ce3  cert/window_cert_output.txt
sha256  990b20686590c2463986b802515a826a358138d44e036f66606d70e016f881bf  cert/crosscheck_mpmath.py
sha256  b63a66d87b12feb1312b0fb41b190fd6cd28321d04f2c6c3ade8018f7cdb9a83  cert/crosscheck_mpmath_output.txt
```

The first two hashes are also in the docstring of `RaceWindowCert`. Reproduce with
`python3 cert/window_cert.py | cmp - cert/window_cert_output.txt` (2 seconds; the output is
deterministic) and `cd cert && sha256sum -c SHA256SUMS`.

Independent cross-check (`cert/crosscheck_mpmath.py`, floating point, not a proof; it shares no
formula with the certified script: `f̂` by numerical quadrature of the ramp, mpmath's digamma,
range `[0, 800]`): `0.161167958443`, which differs from the certified main part by `−1.44e-9`.
That is consistent: the certified bound for everything beyond `r = 800` is below `6e-7`.

### Inputs that remain

Exactly two, both explicit binders of `grh_chi4_race_concrete_window` (never `axiom`s):

1. **`hGRH : GRHFor chi4`**: the zeros of `L(s,χ₋₄)` in the open strip lie on the line.
2. **`hcert : RaceWindowCert`**: `chi4Gamma0 raceWindow < 81/500`, certified outside Lean as above.

The two-sided statements (`…_two_sided`) need no numerical input at all. The general-window theorem
`grh_chi4_smoothed_race_real` keeps its abstract `hcert`; `scripts/numeric_sanity.py` (floating
point, proves nothing) explores other windows.

**Not done: a proof of `RaceWindowCert` inside Lean.** I looked for a cheap route and did not find
one. The inequality has a factor 3.3 of slack (`0.161` against `0.530`), but the crude bounds that
are easy to formalize lose more than that. A quadratic majorant `W(r) ≤ A + Br²` with Plancherel,
or Jensen applied to `Re ψ(z) ≤ log|z|`, both give bounds above `(1/2)(∫f)²`, because the weight is
negative near `r = 0` (`W(0) = −0.844`), exactly where `|f̂|²` has its mass, and grows only
logarithmically. A real proof would need, in Lean, enclosures of `Re ψ(3/4 + ir/2)` on a grid
(Mathlib has almost no quantitative digamma API), a quadrature error bound, the decay of `f̂` and a
Binet-type bound for the tail. That is a separate project (`../advances/FORMAL-MATH-REPO.md` §5
estimates 500 to 1,500 lines for a u-side variant) and was not attempted.

The prime number theorem is **not** a remaining hypothesis: `prime_number_theorem_theta` is proved
in `RaceLib/PNT.lean` from `MediumPNT` (`ψ(x) − x = O(x exp(−c (log x)^{1/10}))`), which zeta23
ports from the PrimeNumberTheoremAnd project in `Zeta23/FromPNTPlus/MediumPNT.lean`, plus Mathlib's
`Chebyshev.isBigO_psi_sub_theta_sqrt`. That zeta23 module is outside zeta23's default build target,
so zeta23's CI does not compile it; it was compiled here (3 min 10 s) and its own
`#print axioms MediumPNT` reports the three standard axioms.

## 3. About GRH

The only claim made is about shape. `GRHFor χ` has the binder structure of the archive statement
`DirichletLZerosOnLine` described in `../GRH-PROOF-DIGEST.md` §1 and §5, instantiated at one
character, and `grhFor_of_dirichlet_zeros_on_line` checks that a term of that type discharges it
(`fun s h0 h1 hz => h q χ s h0 h1 hz`). Nothing here imports, checks or depends on the announced
proof. Note also that the archive is pinned to Lean `v4.33.1` / Mathlib `0df444a3…`, and this
project to Lean `v4.33.0-rc2` / Mathlib `51e6992e…`, so the two cannot be linked without moving
one of them.

## 4. How the true statements differ from the informal plan

- **S0 bound.** The right-hand side is the zero sum at shift 0, `Σ_ρ m_ρ ĝ(γ_ρ)`. By the explicit
  formula it equals `Γ_0` *minus* the two-sided prime sum at shift 0, so in general it is not the
  archimedean term alone. For `χ₋₄` and `supp g ⊂ (−log 3, log 3)` that prime sum is empty
  (`Λ(1) = 0`, `χ₋₄(2) = 0`) and the bound is `Γ_0`.
- **Real part.** The right-hand sides are complex expressions in Lean, so the statements take
  `.re`. Under GRH the zero sum is a sum of nonnegative reals; the real-window versions are stated
  with real numbers throughout.
- **Complex windows.** `f : ℝ → ℂ` and complex `χ` are allowed in S0, so the left side is a complex
  modulus. No conjugate appears on the `S(L)` side: for `L > A` the `conj χ` leg of the explicit
  formula vanishes.
- **Shift condition.** "Beyond the support" is `A < L` where `supp f ⊆ [−A/2, A/2]`
  (so `supp g ⊆ [−A, A]`).
- **Smoothness.** `f ∈ C²` is required because `EF_lit_chi` is stated for `C²_c` test functions.
  The `cos²` ramps in the feasibility table are only `C¹`; use e.g. quintic ramps.
- **S1 route.** `Γ_L → 0` is the Riemann–Lebesgue lemma (`Real.tendsto_integral_exp_smul_cocompact`).
  The prime squares come from one smoothed Abel summation: with `c(n) = Λ(n) χ₋₄(n)` for `n` not
  prime, `Σ_{n ≤ x} c(n) = θ(√x) + O(x^{1/3} log x) = √x + o(√x)`, hence
  `Σ_n c(n) n^{-1/2} g(log n − L) → (1/2)∫g`. Cubes and higher powers are inside that one estimate.
- **Two-sided form.** S1 is proved as `‖P(L) + ĝ(0)/2‖ ≤ Re Γ_0 + ε`, which is stronger than the
  one-sided race statement and needs no numerical input.
- **Module names.** The statement layer is `RaceChallengeDeps` / `RaceChallenge` / `RaceSolution`,
  not `ChallengeDeps` / `Challenge` / `Solution`: zeta23 ships root modules with exactly those
  names, and Lake refuses `import ChallengeDeps` when two packages provide it
  ("could not disambiguate the module").

## 5. What is not formalized

- A proof of `RaceWindowCert` inside Lean (it is certified in Arb, outside Lean; see §2).
- GRH itself, in any form.
- The unsmoothed race, the first Riesz mean, Theorem A of `../CORE-NOTES.md`, anything with LI.
- S1 for characters other than `χ₋₄` (S0 is general; S1 is written for `chi4` only).
- An effective `L₀`: the theorems give existence only.
- `leanprover/comparator` was not run (it is not installed here). `comparator.json` is provided.
  As a weaker substitute, `scripts/gen_challenge.py --check` verifies that the challenge statements
  are textual copies of the solution statements, and `scripts/compare_statements.sh` verifies that
  the 22 statements elaborate to identical types on both sides.

## 6. Layout

| file | role |
|---|---|
| `RaceChallengeDeps.lean` | trusted definitions, `import Mathlib` only |
| `RaceChallenge.lean` | trusted statements, `import Mathlib` only, definitions inlined, proofs `sorry`. **Generated** by `scripts/gen_challenge.py` |
| `RaceSolution.lean` | the same statements, each proved by one delegation to `RaceLib` |
| `RaceLib/S0.lean` | bridges to zeta23; shifted test functions; the S0 inequality |
| `RaceLib/Chi4.lean` | `χ₋₄` is primitive and odd; S0 in closed form |
| `RaceLib/Abel.lean` | smoothed Abel summation (Mathlib only) |
| `RaceLib/NonPrime.lean` | `Σ_{n ≤ x, n not prime} Λ(n)χ₋₄(n) ~ √x` from `θ(x) ~ x` (Mathlib only) |
| `RaceLib/RiemannLebesgue.lean` | `∫ e^{irL} F(r) dr → 0` (Mathlib only) |
| `RaceLib/PNT.lean` | `θ(x)/x → 1` from zeta23's `MediumPNT` |
| `RaceLib/S1.lean` | the race statements, complex and real |
| `RaceLib/Window.lean` | the concrete window: `C²`, support, `∫ = 103/100`, the final theorem |
| `cert/window_cert.py`, `cert/window_cert_output.txt`, `cert/SHA256SUMS` | the Arb certificate for `RaceWindowCert` |
| `cert/crosscheck_mpmath.py` and its output | independent floating-point cross-check |
| `scripts/PrintAxioms.lean`, `comparator.json` | generated from the theorem list |
| `scripts/numeric_sanity.py` | non-certified numerics for other windows |

To read **what is claimed** you need `RaceChallengeDeps.lean` (139 lines) and the statements in
`RaceChallenge.lean`, and you trust Mathlib's `DirichletCharacter.LFunction`, `analyticOrderAt`,
`Complex.digamma`, `ArithmeticFunction.vonMangoldt`, `Chebyshev.theta`, `ZMod.χ₄`, `tsum` and the
Bochner integral. Lean's `∑'` and `∫` are `0` for divergent sums and non-integrable functions;
`zeroSum_summable_and_archTerm_integrable` and `chi4Gamma0_integrand_integrable` rule that reading
out, and the prime sums have finitely many nonzero terms.

## 7. Pins and build

- Lean `leanprover/lean4:v4.33.0-rc2` (`lean-toolchain`)
- Mathlib `51e6992efd06126df61a496bebf8f49482a4e129`
- zeta23: `https://github.com/anthropics/formal-math` at
  `11b6c7e3c927df2658f513b6b68ea7f8929c62c5`, `subDir = "zeta23"`, as a Lake **git** require
  (it worked through the proxy; no path dependency). All revisions are in `lake-manifest.json`.

```bash
export PATH="$HOME/.elan/bin:$PATH"
cd lean/lean-race
lake exe cache get     # Mathlib oleans; do not build Mathlib from source
lake build             # RaceLib, RaceChallengeDeps, RaceChallenge, RaceSolution
lake env lean scripts/PrintAxioms.lean
python3 scripts/gen_challenge.py --check && ./scripts/compare_statements.sh
```

**Do not run `lake clean`.** It deletes the build products of every dependency, Mathlib's cache
included, and the next `lake build` then compiles Mathlib from source. (This happened once during
this work; the build was stopped after 7 minutes and the cache restored with `lake exe cache get`.)

Measured on this machine (4 cores, 16 GB, other jobs running, `nice -n 5`):

| step | wall time |
|---|---|
| `lake update` + Mathlib cache download and unpack (8,681 files) | 2 min 53 s |
| zeta23 cone of `Zeta23.ThmE.MainChi` (first build) | 24 min 23 s |
| `Zeta23.FromPNTPlus.MediumPNT` on top of that | 3 min 10 s |
| second, complete build from the Mathlib cache: 110 zeta23 modules + `RaceLib` | 21 min 7 s |
| this project alone, zeta23 already built (`lake build`) | 39 s |

The largest resident size seen for one `lean` process was about 7 GB
(`Zeta23.ThmE.PrimeSideChi`; most of it memory-mapped Mathlib `.olean` files, system "used" memory
stayed near 3 GB). Output of the last `lake build` ends with

```
✔ [8819/8823] Built RaceLib.Window (13s)
✔ [8820/8823] Built RaceLib (5.6s)
✔ [8822/8823] Built RaceSolution (7.1s)
Build completed successfully (8823 jobs).
```

The only warnings from this project's modules are the 22 deliberate `declaration uses 'sorry'` in
`RaceChallenge.lean`. Build logs are in `logs/` (git-ignored).

`lake env lean scripts/PrintAxioms.lean`:

```
'explicit_formula_zero_sum' depends on axioms: [propext, Classical.choice, Quot.sound]
'zeroSum_summable_and_archTerm_integrable' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_smoothed_prime_power_sum_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_smoothed_prime_power_sum_le_zero_sum' depends on axioms: [propext, Classical.choice, Quot.sound]
'grhFor_of_dirichlet_zeros_on_line' depends on axioms: [propext, Classical.choice, Quot.sound]
'chi4_isPrimitive' depends on axioms: [propext, Classical.choice, Quot.sound]
'chi4_archWeight' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_smoothed_prime_power_sum_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'prime_number_theorem_theta' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_smoothed_race_two_sided' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_smoothed_race' depends on axioms: [propext, Classical.choice, Quot.sound]
'chi4_primeSum_eq_chi4RaceSum' depends on axioms: [propext, Classical.choice, Quot.sound]
'chi4_archTerm_zero_eq_chi4Gamma0' depends on axioms: [propext, Classical.choice, Quot.sound]
'chi4Gamma0_integrand_integrable' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_smoothed_race_real_two_sided' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_smoothed_race_real' depends on axioms: [propext, Classical.choice, Quot.sound]
'raceWindow_piecewise' depends on axioms: [propext, Classical.choice, Quot.sound]
'raceWindow_admissible' depends on axioms: [propext, Classical.choice, Quot.sound]
'raceWindow_integral' depends on axioms: [propext, Classical.choice, Quot.sound]
'raceWindowCert_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_race_concrete_window' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_race_concrete_window_two_sided' depends on axioms: [propext, Classical.choice, Quot.sound]
```

## 8. Proof in one paragraph

Apply `EF_lit_chi_L` to `k(u) = g(u − L)` and to `k = g`. Since `k̂(z) = e^{izL} ĝ(z)` and
`ĝ(γ) = |f̂(γ)|²` for real `γ`, GRH makes every term of the zero sum at shift `L` have modulus equal
to the corresponding (nonnegative) term at shift `0`; the triangle inequality for the absolutely
convergent sum gives S0. For `L` beyond the support the `conj χ` leg of the prime side vanishes,
leaving `−S(L) + Γ_L`. For S1 write `P = (S − Γ_L) + Γ_L − (S − P)`: the first term is bounded by
`Γ_0` (S0), the second tends to `0` (Riemann–Lebesgue), and `S − P`, the sum over proper prime
powers, tends to `ĝ(0)/2` by Abel summation against `Σ_{n ≤ x, n not prime} Λ(n)χ₋₄(n) ~ √x`,
which is the prime number theorem applied at `√x`.

## 9. Licence and attribution

zeta23 is Apache-2.0 (Anthropic, PBC); its `FromPNTPlus/` modules are ported from the
PrimeNumberTheoremAnd project (see zeta23's `NOTICE`). This project depends on zeta23 through Lake
and copies no code from it. The layout follows the Palomar template used in that repository.
