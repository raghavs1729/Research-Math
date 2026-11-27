# lean-allx: the Cesàro mean of Chebyshev's race is positive at every x > 3, under GRH for L(s, χ₋₄) (Lean 4)

A machine-checked proof, in Lean 4 with Mathlib, the public library `zeta23`
(`anthropics/formal-math`) and `monotonicity/lean-chebyshev` of this repository, of:

> **Theorem (`chebyshev_race_cesaro_pos`).** Assume that every zero of the Dirichlet L-function
> L(s, χ₋₄) in the strip 0 < Re s < 1 has real part 1/2. Then for **every** real x > 3
>
> Σ_{p ≤ x} (−χ₋₄(p)) (1 − p/x) > 0,
>
> the sum over the primes p ≤ x, with −χ₋₄(p) = +1 for p ≡ 3 (mod 4), −1 for p ≡ 1 (mod 4), 0 for p = 2.
> Equivalently (`chebyshev_race_cesaro_eq_integral`, `chebyshev_race_integral_pos`):
> ∫₀ˣ (π(t;4,3) − π(t;4,1)) dt > 0 for every x > 3.

So Chebyshev's bias holds on average at every x, not only for all sufficiently large x. The paper
proof is `monotonicity/explicit/EXPLICIT.md` (Theorems 1 and 2); the Lean proof follows its
structure (finite range by computation, analytic range from the explicit formula, Ingham's
second-order partial summation with a positive kernel) but with different, cruder constants; see
"How the Lean proof differs from the paper proof".

The Riemann hypothesis for L(s, χ₋₄) is an **explicit hypothesis** of the theorems (`GRHFor chi4`),
never an axiom. There is **no other hypothesis**: no numerical certificate, no data on zeros, no
linear independence, no simplicity of zeros, no hypothesis on ζ.

Status: `lake build` succeeds; no `sorry` outside the challenge file; all 13 headline theorems
depend only on `propext`, `Classical.choice`, `Quot.sound`; no `native_decide` (hence no
`Lean.ofReduceBool`); no `axiom`. Nobody has reviewed this development other than the Lean kernel
and the agents that wrote it; the statements below are short and are what a reader should check.

## The exact Lean statements

From `AllxChallenge.lean` (which imports only Mathlib; 3 definitions, 13 statements):

```lean
def GRHFor {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → χ.LFunction s = 0 → s.re = 1 / 2

def chi4 : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)

noncomputable def primeCountingMod4 (a : ℕ) (t : ℝ) : ℕ :=
  ((Finset.Icc 0 ⌊t⌋₊).filter (fun p => p.Prime ∧ p % 4 = a)).card

open Finset MeasureTheory

theorem chebyshev_race_cesaro_pos (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x)

theorem chebyshev_race_cesaro_eq_integral (x : ℝ) (hx : 0 < x) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x)
      = (1 / x) * ∫ t in (0 : ℝ)..x,
          ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ))

theorem chebyshev_race_integral_pos (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    0 < ∫ t in (0 : ℝ)..x, ((primeCountingMod4 3 t : ℝ) - (primeCountingMod4 1 t : ℝ))

theorem chebyshev_race_theta1_pos (hGRH : GRHFor chi4) (t : ℝ) (ht : 3 < t) :
    0 < ∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p)

theorem chebyshev_race_cesaro_lower_bound (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x) :
    (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (x - p)) / (x * Real.log x)
      ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x)

theorem chebyshev_race_lead_inheritance (hGRH : GRHFor chi4) (x : ℝ) (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ)) :
    primeCountingMod4 1 x < primeCountingMod4 3 x

theorem chebyshev_race_cone_lower_bound (hGRH : GRHFor chi4) (m m' x : ℝ) (hm : m ≤ 0)
    (hmm' : m < m') (hm'1 : m' ≤ 1) (hx : 3 < x) :
    (m' - m) * x ^ (m - 1) * (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (x - p)) / Real.log x
      ≤ (∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m)
        - x ^ (-(m' - m))
          * ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m'

theorem chebyshev_race_cone_lead_inheritance (hGRH : GRHFor chi4) (m m' x : ℝ) (hm : m ≤ 0)
    (hmm' : m < m') (hm'1 : m' ≤ 1) (hx : 3 < x)
    (h : 0 ≤ ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m') :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (p : ℝ) ^ m

theorem chebyshev_race_cesaro_pos_of_le (x : ℝ) (hx : 3 < x) (hN : x ≤ 10000000) :
    0 < ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime,
      -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * (1 - (p : ℝ) / x)

theorem grh_chi4_psi1_le (hGRH : GRHFor chi4) (t : ℝ) (ht : 2 ≤ t) :
    ∑ n ∈ Ioc 0 ⌊t⌋₊, ArithmeticFunction.vonMangoldt n * ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * (t - n)
      ≤ 27 / 100 * (t + 1) ^ (3 / 2 : ℝ) + 4 * t + Real.log t + 4 + 1 / 3 * Real.log (t + 1)

theorem chebyshev_race_sum_eq_counts (x : ℝ) :
    ∑ p ∈ (Icc 0 ⌊x⌋₊).filter Nat.Prime, -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ)
      = (primeCountingMod4 3 x : ℝ) - (primeCountingMod4 1 x : ℝ)

theorem neg_chi4_natCast (n : ℕ) :
    -((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = if n % 4 = 3 then 1 else if n % 4 = 1 then -1 else 0

theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ
```

Reading notes (what a reader must check to trust that these say what the prose says):

- `chi4` is Mathlib's quadratic character `ZMod.χ₄` (values 0, 1, 0, −1 on 0, 1, 2, 3 mod 4) as a
  `ℂ`-valued Dirichlet character; `χ.LFunction` is Mathlib's `DirichletCharacter.LFunction` (the
  analytic continuation of the L-series). `GRHFor` and `chi4` are character for character those of
  `monotonicity/lean-chebyshev`.
- `(Icc 0 ⌊x⌋₊).filter Nat.Prime` is the set of primes p ≤ x (`⌊x⌋₊` is the natural floor);
  `neg_chi4_natCast` states the sign of each term; `chebyshev_race_sum_eq_counts` identifies the
  unweighted sum with π(x;4,3) − π(x;4,1).
- `∫ t in (0:ℝ)..x, f t` is Mathlib's interval integral (Lebesgue); the integrand is a bounded
  step function.
- `chebyshev_race_lead_inheritance` is "D₁(x) ≥ 0 ⟹ D₀(x) > 0" of `EXPLICIT.md`, Corollary 4, with
  D₀ expressed through the two counting functions.
- `chebyshev_race_cone_lower_bound` and `chebyshev_race_cone_lead_inheritance` are Theorem 5 of
  `EXPLICIT.md` (first inequality and its consequence) for the whole cone of real exponents
  m ≤ 0, m < m′ ≤ 1; `(p : ℝ) ^ m` is the real power `Real.rpow`.
- `chebyshev_race_cesaro_pos_of_le` has no hypothesis at all: it is the finite part of the proof.
- The hypothesis `GRHFor chi4` has exactly the shape that a theorem "all zeros of all Dirichlet
  L-functions in the open critical strip are on the line" would discharge
  (`grhFor_of_dirichlet_zeros_on_line`). Nothing here depends on, or says anything about, any
  announced proof of such a theorem.

## Structure of the proof

Write c(p) = −χ₋₄(p) log p and Θ₁(t) = Σ_{p ≤ t} c(p)(t − p) (`Coupling.Theta1`).

1. **Partial summation** (`AllxLib/Coupling/`, ported unchanged from `monotonicity/lean`):
   x R(x) = Θ₁(x)/log x + ∫₂ˣ Θ₁(t) h″(t) dt with h(t) = (x − t)/log t, h″ > 0
   (`cesaro_convexity_identity`), hence Θ₁ ≥ 0 on [2, x] and Θ₁(x) > 0 imply
   R(x) ≥ Θ₁(x)/(x log x) > 0 (`cesaro_pos_of_theta1`).
2. **Finite range 3 < t ≤ 10⁷** (`AllxLib/Finite/`, no hypothesis): Θ₁(k) > 0 for every integer
   4 ≤ k ≤ 10⁷, by **kernel reduction** (`decide +kernel`) of an integer scan in 1000 chunks of 10⁴
   integers; Θ₁ is affine between consecutive integers. For each prime p the scan finds a with
   2^a ≤ p^256 < 2^(a+1) and *checks both inequalities*, which gives
   a log 2 ≤ 256 log p < (a+1) log 2 and an integer lower bound V(k) for 256 Θ₁(k)/log 2 with
   V(k+1) = V(k) + U(k); the scan checks V(k) > 0. Primality inside a chunk [n₀, n₀+K) is decided by
   a sieve made of big-number operations: the bit mask of the proper multiples of q is
   ((2^(qJ) − 1)/(2^q − 1)) shifted, OR-ed over 2 ≤ q < B with B² ≥ n₀ + K (`sieve_bit_iff`).
   The chunk boundary states are numerals produced by `scripts/gen_chunks.py`; they are not
   trusted: each chunk theorem is checked by the kernel and the semantic invariant is transported
   by `Allx.Good.step` (`scan_sound`).
3. **Analytic range t ≥ 10⁷** (`AllxLib/Analytic/`, under `GRHFor chi4`):
   - `PrimePowers.lean` (`theta1_wlog_ge`): Θ₁(t) ≥ −ψ₁(t,χ) + ∫₀ᵗ [ψ(√u) − (powers of 2) −
     ψ(u^{1/3}) − ψ(u^{1/5})] du — the squares of odd prime powers carry the bias, the rest is
     bounded by Mathlib's Costa-Pereira inequality — and then ψ(y) ≥ (y−1) log 2 − log(y+2),
     ψ(y) ≤ y log 4 + 2√y log y (Mathlib's Chebyshev bounds).
   - `Kernel.lean`, `KernelK.lean`: the C² test function
     K_t(u) = e^{u/2} (H(t − e^u) − (t − e^u)(1 − e^u)₊³), H(y) = ((y+1)₊³ − 2y₊³ + (y−1)₊³)/6
     (H is y₊ averaged against the hat function, so H(y) = y for y ≥ 1 and 0 ≤ H − y₊ ≤ 1/6):
     K_t(log n) = √n H(t − n), its Fourier transform is
     ((t+1)^{s+3} − 2t^{s+3} + (t−1)^{s+3} − 6t)/(s(s+1)(s+2)(s+3)) + 6/((s+1)(s+2)(s+3)(s+4)),
     s = 1/2 + ir, and ‖K̂_t(r)‖ ≤ (t+1)^{3/2}(1/(9/4+r²) + (5/12)m(r)) + (3t/7 + 1/5)m(r) with
     m(r) = 30/((9/4+r²)(25/4+r²)). (A second-difference bound gives the factor (t+1)^{3/2}/|s(s+1)|
     with constant 1.)
   - `ArchWeil.lean` (`arch_weil`): the archimedean term of the explicit formula in Weil's form,
     for any continuous integrable k with |k̂(r)| ≪ 1/(1+r²); `ArchK.lean`: for K_t it is ≤ log t + 3.
   - `ZeroSum.lean`: under GRH, Σ_ρ m_ρ·3/(9/4+γ_ρ²) ≤ 2 − γ − log π + 2A (`sumOne_le`; the kernel
     has a corner in u-space, so the explicit formula is applied to a C² family and a Fatou limit
     is taken) and Σ_ρ m_ρ m(γ_ρ) = 7/2 − γ − 3 log 2 − log π + 5A − 3B (`sumMaj_eq`), where
     A = (L′/L)(2,χ₋₄), B = (L′/L)(3,χ₋₄) as Dirichlet series.
   - `Numerics.lean` (`Sbar_numeric`): (2A + 2 − γ − log π)/3 + (5/12)(7/2 − γ − 3 log 2 − log π +
     5A − 3B) ≤ 27/100 (certified value 0.2548; the true value is 0.1716). Head of A up to n = 47,
     telescoping tails, γ ≥ 0.5752 from H₂₅₅ − 8 log 2, log p in units of log 2 from integer
     inequalities p^k ≤ 2^m. This number bounds Σ_ρ m_ρ/|ρ(ρ+1)| (true value ≈ 0.154).
   - `Psi1.lean` (`psi1chi_le`): the explicit formula `ChebLib.EF_ext` for K_t gives
     ψ₁(t,χ₋₄) ≤ (27/100)(t+1)^{3/2} + 4t + log t + 4 + (1/3) log(t+1) for t ≥ 2.
   - `Elementary.lean` (`final_ineq`), `Main.lean`: (2/3) log 2 = 0.4621 beats 0.27 plus the
     lower-order terms from t = 10⁷ on (total 0.4105 at 10⁷).
4. **Assembly** (`AllxLib/Final.lean`): `theta1_wlog_pos`, `cesaro_w4_ge_and_pos`, the integral
   identity, the inheritance of leads, and the cone of weights (`AllxLib/Coupling/Cone.lean`, the
   deterministic part of Theorem 5, ported from `monotonicity/lean` with one line changed).

## How the Lean proof differs from the paper proof

`EXPLICIT.md` proves Θ₁(t) ≥ L_χ(t) t^{3/2} with sharp constants (2/3 − 0.1556 − …), using the
explicit formula for ψ₁ with all its terms, a computed bound θ(y) ≥ y − 2.07√y up to 10¹⁰ and the
published bound of Broadbent–Kadiri–Lumley–Ng–Wilk beyond; its analytic part works from t = 4892.
None of those inputs is available in Lean, so the formal proof uses

- Chebyshev's elementary ψ(y) ≥ (y−1) log 2 − log(y+2) (bias constant (2/3) log 2 instead of 2/3),
- the closed-form majorant 27/100 for Σ_ρ m_ρ/|ρ(ρ+1)| instead of 0.1556,
- a C² smoothing of the Riesz kernel instead of the exact formula for ψ₁ (only an upper bound
  for ψ₁ is proved),

and consequently an analytic threshold of 10⁷ instead of 4892, which is why the kernel check goes
up to 10⁷. The conclusion "R(x) > 0 for every x > 3" is the same.

## What is not formalized

- The quantitative statements of `EXPLICIT.md`: Θ₁(t) ≥ L_χ(t) t^{3/2} (Theorem 1), the refined
  lower bound of Theorem 2, Corollary 3 (R(x) > √x/(2 log x) for x ≥ 23). Only positivity and the
  inequality R(x) ≥ Θ₁(x)/(x log x) are formalized.
- Of Theorem 5 (the cone of weights) the first lower bound and the inheritance of leads are
  formalized; the second lower bound, in terms of `L_χ`, is not (it needs Theorem 1's `L_χ`).
- Theorem 6 (the converse, by Landau's theorem; hence the equivalence with GRH) and Theorem 7 (the
  character mod 3). The mod-3 analogue would need the whole analytic chain for χ₋₃ (lean-chebyshev
  only treats χ₋₄).
- The exact explicit formula for ψ₁(t,χ) as an identity; a two-sided bound for ψ₁.
- The official Comparator run (`leanprover/comparator` with `lean4export`, `nanoda`, `landrun`) was
  not performed; `comparator.json` is provided, and `scripts/compare_statements.sh` checks that the
  13 statements and the 3 definitions elaborate identically in the challenge and solution modules.
- No human has read the proofs. The trusted base is: the Lean kernel (including its GMP-backed
  evaluation of natural-number arithmetic, which `decide +kernel` relies on), Mathlib, and the
  definitions and statements in `AllxChallenge.lean`.

## Files

| file | lines | content |
|---|---|---|
| `AllxChallenge.lean` | 168 | **trusted**: definitions and the 13 statements with `sorry`; imports only Mathlib (generated) |
| `AllxSolution.lean` | 143 | the same statements, proved by one-line delegation |
| `AllxLib/Final.lean` | 129 | assembly |
| `AllxLib/Coupling/Abel.lean`, `Convexity.lean`, `Cone.lean` | 613 | second-order partial summation, the cone (from `monotonicity/lean`) |
| `AllxLib/Finite/Defs.lean`, `Sieve.lean`, `Scan.lean` | 665 | kernel programs, sieve correctness, soundness of the scan |
| `AllxLib/Finite/Chunks01–08.lean`, `Chain*.lean`, `Final.lean` | generated | 1000 `decide +kernel` chunk theorems and their chain |
| `AllxLib/Analytic/Kernel.lean`, `KernelK.lean` | 929 | the test function and its transform |
| `AllxLib/Analytic/ArchWeil.lean`, `ArchK.lean` | 785 | archimedean term |
| `AllxLib/Analytic/ZeroSum.lean`, `Numerics.lean` | 1331 | zero sums in closed form and their numerical bound |
| `AllxLib/Analytic/Psi1.lean`, `PrimePowers.lean`, `Elementary.lean`, `Main.lean` | 883 | ψ₁ bound, prime powers, final inequality, analytic range |
| `scripts/gen_chunks.py`, `gen_numerics.py`, `gen_challenge.py` | | generators of the chunk files, of `Numerics.lean`, of the challenge file |
| `scripts/sim_scan.py`, `budget.py` | | planning aids (floating point; not part of the proof) |
| `scripts/PrintAxioms.lean`, `compare_statements.sh`, `comparator.json` | | audit |
| `logs/` | | build output, axiom audit, statement comparison |

## Pins and build

| | |
|---|---|
| Lean | `leanprover/lean4:v4.33.0-rc2` (`lean-toolchain`) |
| Mathlib | `51e6992efd06126df61a496bebf8f49482a4e129` |
| zeta23 | `https://github.com/anthropics/formal-math`, commit `11b6c7e3c927df2658f513b6b68ea7f8929c62c5`, `subDir = "zeta23"` (Apache-2.0) |
| lean-chebyshev | `../lean-chebyshev` in this repository (path dependency; `ChebLib`) |

```
cd monotonicity/lean-allx
lake exe cache get                               # Mathlib oleans; never build Mathlib from source
lake build AllxChallenge AllxSolution            # zeta23 modules, ChebLib, then this project
lake env lean scripts/PrintAxioms.lean           # axiom audit, expected output in logs/print-axioms.txt
./scripts/compare_statements.sh                  # challenge vs solution statements
python3 scripts/gen_challenge.py --check         # the challenge file is the generated one
python3 scripts/gen_chunks.py 1e7                # regenerates AllxLib/Finite/Chunks*.lean, Chain*.lean (optional)
```

Do not run `lake clean` (it removes the Mathlib cache). On 4 cores: the needed zeta23 modules and
ChebLib take about 25 minutes; this project's own modules about 7 minutes, of which the 1000 kernel
checks are about 6 (each chunk 1.0–1.4 s of kernel time, about 620 MB per chunk file; the chunk
files set `Elab.async false` because running the 125 checks of a file concurrently exhausts the
memory). See `logs/build.txt`.

## Provenance

Written by an AI agent (Claude Code) in one session on 2026-11-24, with sub-agents proving
individual lemma files against fixed statements (finite check, Weil form of the archimedean term,
prime-power bookkeeping and the final inequality, numerical bounds); all files were compiled
together and the axiom audit was run on the final build. zeta23 is Apache-2.0. The prime number theorem (zeta23's
`MediumPNT`, used by lean-chebyshev) plays no role in the argument here: the bias comes from
Chebyshev's elementary lower bound for ψ.
