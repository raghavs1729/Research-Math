AI-produced (Anthropic's Claude), adversarial referee report, 2026-11-26.

# Referee report 1 on `monotonicity/density2/TOY.md`

**Version reviewed.** I reviewed `TOY.md` at commit `6e037df` (2026-11-26 23:18 UTC, 205 lines). That
is the file I read at the start.

While this review was in progress, the file was changed by commits `b4424eb`, `6cf6671` and
`47afd65` (the last at 2026-11-27 00:01 UTC); it now has 315 lines. The changes are:
- a new §2.2 with **Theorem 2.3** (single-hole Toeplitz weights), which is **not covered by this
  report**;
- a new summary and section headings: §2 is now "Beyond bounded variation: the Grandi mechanism and
  Toeplitz weights", which disposes of my objection to the old §2 title;
- a rewording of the last sentence of §2.1 (see §5 and R3 below);
- a new paragraph after Conjecture 5.1.

Every other passage quoted below is unchanged in the latest version. Section numbers 1, 3, 4 and 5 are
the same in both versions.

**Scope.** I read all of `TOY.md`. For context only, I also read `DENSITY.md` Lemma 1.1 and Remark 7.5.
I checked every proof line by line.

I recomputed everything numerical with new scripts in `monotonicity/density2/referee1/`; each
script's output is in the matching `*_output.txt`. I did not read or reuse the authors' scripts.

I checked the Breuer–Simon citation of §4 against the text of arXiv:1002.0823v2.

**Bottom line.**
- I found **no error in any labelled result**. Lemma 1.1, Theorem 1.2, Lemma 1.3 and Propositions
  2.1, 3.1 and 3.2 are correct.
- Three proofs have small gaps, each fixable in a few lines: Theorem 1.2 step (iii), Lemma 1.3,
  and the displayed estimate in the proof of Proposition 3.1.
- Two side statements are **false as written**:
  - "the a.c. part of `w` is an integer trigonometric polynomial" (Remark (3), and again in §2);
  - "For periodic `η` the slope is `±1/2`" (§3).
- One strategic conclusion in §2 is stated as fact but not proved.
- Conjecture 5.1 has a trivial counterexample as stated.

I also give two short supplements:
- **O1** strengthens the parity conclusion of Theorem 1.2.
- **O2** shows that, under the hypotheses of Theorem 1.2, the weights live on finitely many arithmetic
  progressions. So `F` factors through a finite-dimensional torus. In particular this **answers the
  question left open in Remark (2)**: `w` is never a finite measure for the true race.

## Verdicts

| Item | Verdict |
|---|---|
| Lemma 1.1 (mod-2 lemma) | CORRECT |
| Theorem 1.2 | CORRECT WITH GAP (step (iii); fix in §2 below) |
| Lemma 1.3 | CORRECT WITH GAP (complex `u`; extension of (1.1) to Borel `f`) |
| Remarks (1)–(3) after Lemma 1.3 | (3) WRONG as written (factor 2); (1) overstated; (2) question answerable (O2) |
| Proposition 2.1 | CORRECT |
| §2 prose after Proposition 2.1 | OVERSTATED (unlabelled heuristic conclusion; misquotes Theorem 1.2) |
| Proposition 3.1 | CORRECT WITH GAP (displayed `L²` estimate wrong as written; conclusion right) |
| Proposition 3.2 | CORRECT |
| §3 remark "For periodic `η` the slope is `±1/2`" | WRONG (counterexample in §7) |
| §4 | CORRECT (citations verified) |
| §5 numerics [num] | REPRODUCED independently (with a grid caveat) |
| Conjecture 5.1 | needs `|W| ≥ 2` and `M`-dependent constants |

---------------------------------------------------------------------------

## 1. Lemma 1.1 — CORRECT

**Proof check.**
- With `|ℓ_{i+1}| > Σ_{j≤i}|ℓ_j|` (the factor 2 is not needed), the `2^r` subset sums are distinct.
  To see this, look at the largest index in the symmetric difference of two subsets.
- `Δ_{ℓ_i} 1_{g_i+L_i} = 0` pointwise, and the operators `Δ` commute and are `ℤ`-linear.
- `Δ_{ℓ_1}⋯Δ_{ℓ_r}1_{{0}} = Σ_T (−1)^{r−|T|} 1_{{−ℓ_T}}`.
- If `f − 1_{{0}} = 2h`, then `Δ⋯1_{{0}} = −2Δ⋯h` would be even everywhere. That is a contradiction.
- The case `r = 0` is trivial.
- Implicit hypothesis: `L_i ⊂ G` (or work on `ℝ`).

Torsion-freeness is essential, and the statement does use `G ⊂ ℝ`. In `(ℤ/p)²` with `p` odd,
`1_{{0}} ≡ 1_G + Σ_{p+1 lines through 0} 1_line (mod 2)`.

**Numerics (`referee1/t1_mod2_lemma.py`).**
- (A) `G = ℤ`, subgroups `qℤ`: 20000 random configurations, 0 with `f ≡ 1_{{0}}`. The proof mechanism
  was verified each time: distinct subset sums, `Δ⋯f ≡ 0`, and `Δ⋯1_{{0}}` odd at `2^r` points.
- (B) `G = ℤ + αℤ ≅ ℤ²`, GF(2) linear algebra on the box of radius 34. The family has 185 cosets:
  lines in 11 directions and 7 finite-index sublattices. The margin 31 exceeds `Σ_j|ℓ_j|_∞ = 26`, so a
  window solution would contradict the proof, which gives the test teeth.
  - The span contains **no** nonzero function supported in the box of radius 3.
  - In particular `1_{{0}}` is not in the span.
  - This also tests the stronger statement O1.
- (C) Control: in `(ℤ/3)²` and `(ℤ/5)²`, `1_{{0}}` **is** in the span (explicit identity above
  verified); in `(ℤ/2)²` it is not.

## 2. Theorem 1.2 — CORRECT WITH GAP (minor)

**(i) Cohen and the decomposition `2w = Pλ_A + σ`: correct.**
- `|2ŵ| ≤ ‖2w‖`. An integer-valued bounded function takes finitely many values `V`.
- The Lagrange polynomials, taken in the convolution algebra `M(A)` with unit `δ_e`, give idempotents
  with transform `1_{2ŵ=v}`.
- Cohen's theorem for compact `A` with discrete dual `G` puts each level set in the coset ring `R(G)`.
  `G` is discrete, so every subgroup counts.
- The `ℤ`-span of coset indicators is closed under products, because the intersection of two cosets is
  empty or a coset. Hence every `1_E`, `E ∈ R(G)`, is a `ℤ`-combination of coset indicators, and
  `2ŵ = Σ n_i 1_{g_i+L_i}` follows.
- The transform of `e_g m_{L^⊥}` is `1_{g+L}`, with the convention `μ̂(γ) = ∫ ē_γ dμ`.
- Singularity for `L ≠ {0}`:
  - `L` is torsion-free, hence infinite, and `A/L^⊥` is dual to `L`, so `[A : L^⊥] = ∞`.
  - A closed subgroup of positive Haar measure is open (Steinhaus), hence of finite index. So
    `λ_A(L^⊥) = 0`, and `e_g m_{L^⊥} ⊥ λ_A`.
- Uniqueness of Fourier–Stieltjes transforms gives `2w = Σ n_i e_{g_i} m_{L_i^⊥}` as measures.
- `P` is well defined: it is the density of the a.c. part.

Minor: "(Rudin, *Fourier analysis on groups*, §3.1; ...)" should cite the theorem number of Cohen's
theorem, not the section.

**(ii) Correct.** `(D+M)h = Mσ`, so `Dh = Mσ − Mhλ_A`. Its a.c. part is `−Mh`, and Lemma 1.3 applies to
`h ∈ L² ⊂ L¹`. `h` is real, since `2w` is a real measure and `R_M` preserves real functions.

**(iii) Gap.** I object to:

> "for a.e. `θ ∈ S`, `t = 0` is a Lebesgue density point of `{t : τ_tθ ∈ S}` and a point of
> differentiability of `F_θ`, so the derivative there is `0`."

- Lemma 1.3 only provides *some* BV function that agrees with `F_θ` a.e. in `t`. "Differentiable at
  `t = 0`" and "`F_θ(0) = c`" refer to one particular representative.
- To pass from "a.e. `t`, for each `θ`" to "`t = 0`, for a.e. `θ`", that representative must be
  covariant (`F̃_{τ_rθ} = F̃_θ(·+r)`) and must satisfy `F̃_θ(0) = F(θ)` for a.e. `θ`. Neither is arranged.

The fix does not need `t = 0` at all.
- For a.e. `θ`, take the right-continuous representative `F̃_θ`. It is canonical, hence covariant.
- `F̃_θ` equals `F(τ_tθ)` for a.e. `t`, and is differentiable at a.e. `t` with derivative
  `M(P − F)(τ_tθ)`.
- Let `Z_θ := {t : τ_tθ ∈ S}`. Almost every `t ∈ Z_θ` is a density point of `Z_θ` at which all of this
  holds. Approaching `t` through points of `Z_θ` where `F̃_θ = c` gives derivative `0`.
- Hence `(P − c)(τ_tθ) = 0` for a.e. `t ∈ Z_θ`. Fubini and invariance give `P = c` a.e. on `S`.

**Why `P ≡ c` follows.**
- `P − c` is a trigonometric polynomial vanishing a.e. on a set of positive measure. By `DENSITY.md`
  Lemma 1.1, all its coefficients vanish.
- That lemma is stated for the race hull, but its proof uses only that a finitely generated subgroup of
  `G` is `≅ ℤ^n` and duality. It applies verbatim to any `A` of §0; say so.
- "`c = n_0`": `n_0` is undefined. Write `c = P̂(0)`.

**Parity conclusion: correct.**
- `2ŵ(0) = 1` and `2ŵ(g) = 2m(g)` is even for `g ≠ 0`, so `σ̂ = 2ŵ − c1_{{0}} ≡ (1−c)1_{{0}} (mod 2)`.
- `σ̂ = Σ_{L_i≠0} n_i 1_{g_i+L_i}`, so Lemma 1.1 forbids `1 − c` odd.
- (O1 gives a stronger and cleaner form.)

**Numerics (`referee1/t2_theorem12_examples.py`): no counterexample, consistent with the theorem.**

*Circle.* For symmetric periodic weights `m(n) = a[n mod q]` I used the exact closed form of `F` on each
arc: `F` is flat on arc `k` iff a polynomial `p_k(x)`, `x = e^{−2πM/q}`, vanishes. I enumerated all
symmetric patterns with `q ≤ 8`:

| Weights | Patterns | Flat arcs (over all `M > 0`) | Levels |
|---|---|---|---|
| `{0,1}` | 90 | 135 | `±1` only |
| `{0,1,2,3}` | 1700 | 4501 | `−5, −3, −1, 1` only |

- Samples were confirmed by Fejér means of the Fourier series, independently of the closed form.
- Example: `m(n) = 1` (`n` even), `2` (`n` odd), `M = log3/π` gives `F = −1` on `(π, 2π)`.
- Adding random finite perturbations (`w` is still a measure, `P` non-constant) destroys flatness: the
  minimal oscillation on the former flat arc is 0.25.

*Rank-2 hull `A = 𝕋²`, `G = ℤ + αℤ`.* Explicit example:
- `α = log2/log3` (irrational).
- Weights `1,3` on even/odd elements of `ℤ`, and `1,2` on even/odd elements of `αℤ`. Then `w` is a
  measure (Haar measures on closed subgroups).
- `M = log2/π`.

Results:
- `F = −3` on the quarter `(π,2π)²` of `𝕋²`. The grid fraction is 0.2494, and Fejér means give −2.9998.
- Along 20 random orbits, `dF/dt = M(P − F)` with `P = −3`, to 2e-5. This is exactly step (ii).
- Random search over 400 rank-2 examples:
  - In 200 of them, `α` is chosen so that both directions are flat at the same `M`. Genuine flat sets
    appear (grid fraction 0.167 for all `τ ≤ 1e-4`), and their levels are `{−7,−5,−3,−1,1}`.
  - In every draw, the fraction of points within `τ` of an **even** integer goes to 0 as `τ → 0`.

**Why no rank-2 counterexample can exist** (beyond the proof itself): see O2.
- Under the hypotheses, `F` is a finite sum of one-dimensional sawtooth-exponential functions of the
  phases of finitely many characters, plus a trigonometric polynomial.
- Within each cell, `F − c` is a combination of exponentials of distinct linear forms. So flatness
  forces each direction to be flat separately, and the level is `P̂(0)`, which is odd.

**Remark (3) is WRONG as written.**

> "(3) The proof shows exactly where integrality enters: the a.c. part of `w` is an integer
> trigonometric polynomial, while `ŵ(0) = 1/2`."

- The a.c. part of `w` is `P/2`.
- In the flat case `P = c` is odd, so the a.c. part is `c/2 ∉ ℤ`.
- By O1, `P̂(0)` is **always** odd, so the a.c. part of `w` is *never* an integer polynomial.
- The intended statement is about `2w`. The same slip recurs in §2 (see §5 below).

**Remark (1) is overstated.**

> "On the circle this is Helson's theorem plus the parity argument of `DENSITY.md` Remark 7.5, now for
> sets of positive measure and every hull."

- Theorem 1.2 adds a hypothesis that Remark 7.5 does not make: `w` a finite measure.
- Remark 7.5 allows integer weights `O(n^k)` on an arc and gets eventual quasi-polynomiality from
  Pólya–Carlson. Polynomially growing weights are not covered by Theorem 1.2.
- So Theorem 1.2 is complementary to Remark 7.5, not an extension of it.

**Remark (2).** "whether that is excluded was not checked" — it can be settled; see O2.

## 3. Lemma 1.3 — CORRECT WITH GAP (minor)

What I checked:

| Point | Finding |
|---|---|
| (a) `Du_ε = φ_ε⋆μ` | Correct. Computed against trigonometric polynomials, which are dense in `C(A)`, so `‖Du_ε‖_{L¹} ≤ ‖μ‖`. |
| Variation bound | Correct. `(u_ε)_θ = φ_ε * u_θ`; then Fubini and invariance; lower semicontinuity of the essential variation (a countable supremum) plus Fatou. |
| Flow-invariance of `A_0` | Correct. |
| (b) Vague convergence and domination | Correct. `|∫ψ(t)f(τ_tθ) d(φ_ε*ν_θ)| ≤ ‖ψ‖_∞‖f‖_∞|ν_θ|(supp ψ+[−ε,ε])`, integrable by (a). This needs `φ_ε ≥ 0` with support in `[−ε,ε]` (say so). |
| (c) Lebesgue decomposition | Correct, and it really gives `ρ_θ(t) = ρ(τ_tθ)`. The singular part of `ν_θ` is carried by `{t : τ_tθ ∈ Z}`. `Z` is `λ_A`-null by Fubini. On `{τ_tθ ∉ Z}`, `ν_θ = ρ̃(τ_tθ)dt`. So `μ = ρ̃λ_A + (measure on Z)`. Uniqueness of the decomposition gives `ρ̃ = ρ` a.e., and Fubini transfers this to a.e. orbit. |
| Fubini/invariance | No misuse found. Joint measurability of `(θ,t) ↦ u(τ_tθ)` uses metrizability of `A` (Borel(A×ℝ) is the product σ-algebra), which is in §0. |

Gaps to close:
1. > "`ρ̃(φ) := limsup_{δ→0, δ∈ℚ} ν_φ((−δ,δ))/2δ`"

   This is undefined for complex `ν_φ`, yet the lemma allows complex `u`. Treat `Re u` and `Im u`
   separately (the lemma is linear), or restrict to real `u`: that suffices, since `h` is real.
2. > "Both sides of (1.1) are finite measures in `f`, so (1.1) holds for bounded Borel `f`."

   The left side is not yet defined for Borel `f`: measurability of `θ ↦ ∫ψ(t)f(τ_tθ)dν_θ(t)` is needed.
   The fix is a monotone-class argument.
   - Consider the class of bounded Borel `f` for which this map is measurable and (1.1) holds.
   - It contains `C(A)` and is a vector space.
   - It is closed under bounded pointwise limits by dominated convergence: in `t` because
     `|ν_θ|(supp ψ) < ∞`, and in `θ` by the bound from (a).
3. `u_ε` must be defined pointwise by the integral (for `θ` with `u_θ ∈ L¹_loc`), not as an `L¹` class.
   Otherwise "`(u_ε)_θ = φ_ε * u_θ`" is only an a.e. statement.
4. Notation: `φ` denotes a point of `A` in (c), but a test function in the statement and in (a).
5. Theorem 1.2 uses that the pointwise a.e. derivative of a good representative equals the density of
   the a.c. part of `ν_θ`. This is standard, but say it.

(The Euclidean analogue is the slicing theorem for BV functions in the Ambrosio–Fusco–Pallara monograph;
citing it would help the reader.)

## 4. Proposition 2.1 — CORRECT

**Checked.**
- `ν₂(2n) = ν₂(n)+1` gives `ŵ(2n) = 1 − ŵ(n)` (also at `n = 0`) and `ŵ(odd) = 0`.
- `w(t)+w(2t) = π(δ_0+δ_π)` holds with the normalisation `δ̂_0(n) = 1/2π`. State this, since the
  formula for `2πM'R_{M'}δ_0` depends on it.
- `R_M[u(2·)] = ½(R_{M/2}u)(2·)` holds.
- The iteration, and the bound `‖F_{M'}−1‖_2 ≤ 2√2 M'(Σn^{−2})^{1/2}`, are correct.
- `sup|Q_{M'}−2| = O(M')`; numerically `sup|Q_{M'} − 2|/M' → π`.
- The constant: `2Σ_{k<K}(−1)^k + (−1)^K F_{M/2^K}(2^K·) → 1` in `L²`. Uniform convergence of the
  series, hence a.e. equality and boundedness.

**Numerics (`referee1/t3_prop21_period_doubling.py`).**
- The coefficient identity holds exactly for `n ≤ 2·10^5` and 4 values of `M`.
- The expansion agrees with Fejér means of the Fourier series (`N = 2·10^5`) at random points to
  ≤ 3.4e-4, which is the expected Fejér error.
- `F_M(t) + F_{M/2}(2t) = Q_M(t)` to 1e-15.
- `‖F_M − 1‖²` from the expansion equals the Parseval value; e.g. 0.627982 for `M = 1/2`.
- `sup|F_{1/2}| ≈ 2.36`.

Optional addition, making the point of §2 precise: `w` is **not** a measure here.
- By Helson, because `{n : ν₂(n) odd}` is not eventually periodic.
- Equivalently, `F_M` is not of bounded variation. It jumps at every dyadic point:
  - at `0` and `π` by `4πM/3` (verified numerically at 0);
  - at each of the `2^a` points `πj/2^a` (`j` odd, `a ≥ 1`) by `(4πM/3)2^{−a}`.

  So the total jump mass is `Σ_a 4πM/3 = ∞`.

## 5. §2 prose after Proposition 2.1 — OVERSTATED

> "and the obstruction "a.c. part of `w` is an integer polynomial" of Theorem 1.2 has no analogue."

This misquotes Theorem 1.2 (it is `2w`; see Remark (3) above). "Has no analogue" is not a statement
that anything proves.

> "So integrality-by-parity is not a mechanism that can reach the true race."

- Nothing in the note proves this. Proposition 2.1 is one example in which `w` is not a measure.
- It does not show that no parity-type argument can go beyond bounded variation.
- Under the labelling convention of the header this sentence, and the paragraph containing it, must be
  marked **Heuristic**.
- The section title "why the parity argument stops at bounded variation" claims the same and should
  be softened.
- What *can* be proved is narrower: Theorem 1.2 cannot be applied to the true race (O2).
- In the latest version (`47afd65`; §2 retitled) the sentence reads "So the parity argument of
  Theorem 1.2, *as it stands*, cannot reach the true race." This narrower claim is true; O2 proves it. But the note still
  derives it from the Grandi discussion, which is not a proof. Add the argument of O2, or label the
  sentence.

"(and, on the hull, by the prime squares)" is unclear. The prime squares create the positive mean density
`½` of the prime atoms through the explicit formula; say so or drop it.

The displayed formula `(D+M)E^θ_M = 2MΣ_p χ*(p)(log p)p^{−1/2}δ_{log p}` is consistent with
`GRH-ONLY.md` §0 (`M = m + ½`).

## 6. Proposition 3.1 — CORRECT WITH GAP

> "`‖F_{m^{(j)}} − F_{m^∞}‖_2² ≤ 2(max W)² Σ_{n>N_j}|a_n|² + (finitely many coordinates) → 0`"

As written this does not follow.
- The coordinates `n ≤ N_j` where `m^{(j)}` and `m^∞` differ are not "finitely many" uniformly in `j`.
- Example: `m^{(j)} = 1` on `(j, N_j]` and 0 elsewhere, with `m^∞ = 0`. The left side is
  `2Σ_{j<n≤N_j}|a_n|²`, which neither term controls.

Correct version: for fixed `K`,
`‖F_{m^{(j)}} − F_{m^∞}‖² ≤ 2Σ_{n≤K}|m^{(j)}_n − m^∞_n|²|a_n|² + 2(max W)²Σ_{n>K}|a_n|²`.
The first term vanishes for `j ≥ j(K)` along the subsequence. Then let `K → ∞`.

The rest is correct:
- Weak convergence of the laws comes from `L²` convergence.
- Portmanteau applied to the closed set `{|F| ≤ τ}`.
- `τ → 0` by continuity from above.
- The Chebyshev converse.

**Quantifiers on `η`: correct.**
- The forward direction gives `λ{F_m=0} ≥ 2πη` with the same `η`. It only needs `liminf τ*_η(N) = 0`.
- The converse needs `η' < η`, because Chebyshev loses an arbitrarily small amount. This is the right
  form.
- "Equivalently ..." is correct: `s(τ)` is monotone, and with infinite sequences dominated convergence
  gives the `L²` limit.
- Integrality of `W` is never used; any finite `W ⊂ [0,∞)` works.

## 7. Proposition 3.2 — CORRECT; the remark after it is WRONG

**The proof is correct.**
- `Ψ' = w` and `(F/2M + I/2)' = w` as distributions on `ℝ`, so they differ by a constant.
- `2m_n = 1 + η_n`.
- `Σ_{n≥1} sin(nt)/n = (π−t)/2` on `(0,2π)`.

Numerics (`referee1/t4_prop32_B_eta.py`):
- The identity `Ψ − F/2M − I/2 = const` holds exactly for random 0/1 truncations (spread 1e-14).
- `Σ_{n≤10^6} sin(nt)/n = (π−t)/2` to 3e-6.
- The mechanism generalises: on a flat set `{F = c}`, `B_η` (or `Ψ`) has slope `c/2`. I checked this on
  three flat-arc examples, with `c = ±1`.

> "For periodic `η` the slope is `±1/2`."

**This is false** unless `η_r = η_{q−r}` for all `r`. That condition is equivalent to the two-sided
weights `m_{|n|}` being periodic, i.e. to `w` being a measure.

For `|z| = 1` off the `q`-th roots of unity, with `N(z) = Σ_{r=1}^q η_r z^r`:

    B_η' = Re N(z)/(1−z^q) = −η_q/2 + Σ_{1≤r<q/2} (η_r − η_{q−r}) sin((q/2−r)t)/(2 sin(qt/2)).

Counterexample: `η = (1,−1,−1)` repeated (`n ≡ 1,2,0 mod 3`), i.e. the toy with `m_n = 1_{n≡1 (3)}`.
- Then `B_η' = ½ + 1/(1+2cos t)`.
- `B_η` has **logarithmic singularities** at `2π/3` and `4π/3`, so it is not even bounded.
- Numerically `B(t_0−δ) − B(t_0−10δ) → ln10/√3`: 1.3324 at `δ = 10^{−3}`, against 1.3294 predicted.
- Difference quotients match `∫(½ + 1/(1+2cos))` to 1e-5.

The quadratic-variation identity is correct, and was checked by Parseval and on an FFT grid. It holds for
`0 ≤ h ≤ 2π`; state the range.

## 8. §4 — CORRECT

- `F = 1 + 2Re g` on `𝕋`, so the pseudocontinuation equivalence is right.
- Breuer–Simon, arXiv:1002.0823v2, checked against the text:
  - Thm 1.3: continuation across an open interval `I` ⇒ right limits reflectionless on `I`.
  - Thm 1.4: the `H¹`-type bound (1.8) on `I` suffices.
  - Thm 5.1: finitely many values ⇒ strong natural boundary, or eventually periodic.
  - All are for intervals (arcs), as the note says.
- Minor: "`H²` norm `≍ √N`" needs the ones of `m` to have positive density near `N`.

## 9. §5 numerics and Conjecture 5.1

**Reproduction with my own code (`referee1/t5_numerics_check.py`, `t5b_exhaustive_N22.py`).**

*Exhaustive search, `N = 22`, grid `16N`:*

| `η` | .05 | .10 | .15 | .20 | .30 |
|---|---|---|---|---|---|
| `τ*` (mine) | 5.98e-3 | 1.37e-2 | 3.21e-2 | 5.96e-2 | 0.165 |

These match the exhaustive values recorded in `PROGRESS.md` §6.3; `TOY.md` itself prints none.

*Beam search (my own, width 200):*

| `η` | 0.40 | 0.25 | 0.05 |
|---|---|---|---|
| `TOY.md` table, `N = 1024` | 0.282 | 0.0615 | 0.0017 |
| mine, `N = 1024` | 0.2820 | 0.0620 | 1.65e-3 |
| mine, `N = 256` | 0.284 | 0.0647 | 2.8e-3 |

**Grid caveat.** "(`M = 1/2`, grid `8N`–`16N` points)" should say which grid each number used.
- The beam outputs are robust: re-evaluated on a 64× finer grid, they change by at most 7%
  (1.65e-3 → 1.76e-3 for `η = 0.05`; under 2% otherwise).
- An unstructured optimiser *can* exploit an `8N` grid. At `N = 22`, `η = 0.05`, exhaustive search on an
  `8N` grid returns 3.8e-3, but the true quantile of that sequence is 1.7e-2, a factor 4.4 larger.
- The tabulated `τ*` are optimiser values, i.e. upper bounds for `τ*_η(N)`. "The values have stopped
  moving (plateau)" is a statement about the optimiser, not a lower bound.

Not checked: the `N = 2048` values, the `{0,1,2}` values, and the real-weight value 2e-4.

> "moreover `sup_m λ{|F_m| ≤ τ}/2π ≍ log|W| / log(1/τ)`."

This is false for `|W| = 1`, `W ≠ {0}`.
- `W = {1}` gives the single function `F = 4πMe^{−Mt}/(1−e^{−2πM}) − 1` on `(0,2π)`. For `M = ½` it has a
  simple zero, so the left side is `≍ τ` while the right side is 0.
- The constants must also depend on `M`. `PROGRESS.md` reports that `M = 0.125` plateaus far from 0.
  More decisively, `‖F_M − 1‖_2² ≤ 8ζ(2)(max W)²M²` uniformly in `m`. By Chebyshev,
  `sup_m λ{|F_m| ≤ τ} → 0` as `M → 0` for each fixed `τ < 1`, while the right side does not depend on `M`.

---------------------------------------------------------------------------

## Required fixes

- **R1** (Remark (3); §2). Replace "the a.c. part of `w` is an integer trigonometric polynomial, while
  `ŵ(0) = 1/2`" with: "the a.c. part of `2w` is an integer trigonometric polynomial `P`, and parity
  forces `P̂(0)` odd (while `2ŵ(0) = 1`)". Correct the quotation of Theorem 1.2 in §2 the same way.
- **R2** (§3). Replace "For periodic `η` the slope is `±1/2`" with: "for periodic `η` with
  `η_r = η_{q−r}` (equivalently `w` a measure) `B_η` is piecewise affine with slope `−η_q/2 = ±1/2`;
  otherwise `B_η` has logarithmic singularities". Give the formula above, or the `q = 3` counterexample.
- **R3** (§2). Label the paragraph after Proposition 2.1 as **Heuristic**, in particular "So
  integrality-by-parity is not a mechanism that can reach the true race". Soften the section title.
  Optionally replace the sentence with the provable statement of O2.
  - The latest version (`47afd65`) has already narrowed the sentence to "... *as it stands*, cannot
    reach the true race", and has retitled §2. The narrowed sentence is true, but it needs the proof of
    O2 (or a label).
  - The "has no analogue" sentence is unchanged and still needs fixing.
- **R4** (Proposition 3.1). Replace the displayed `L²` estimate by the split at a fixed `K` given in §6.
- **R5** (Theorem 1.2 (iii)). Replace the "`t = 0` / point of differentiability" sentence by the
  a.e.-`t` argument along orbits in §2, or fix a covariant representative explicitly. Define `n_0`
  (it is `P̂(0)`).
- **R6** (Lemma 1.3). Handle complex `u`, or state the lemma for real `u`. Justify the extension of
  (1.1) to bounded Borel `f` by the monotone-class argument. State `φ_ε ≥ 0`.
- **R7** (Conjecture 5.1). Add `|W| ≥ 2`, and say that the implied constants depend on `M` (fixed `M`).

## Optional improvements

**O1 (sharper parity statement).**
The proof of Lemma 1.1 gives more. If `Σ n_i 1_{g_i+L_i}` (all `L_i ≠ {0}`) is `≡ Q (mod 2)` with `Q`
finitely supported, then `Q ≡ 0 (mod 2)`.

*Proof.* In the group ring `𝔽_2[G]`, which is an integral domain because `G` is torsion-free abelian,
`Δ_{ℓ_1}⋯Δ_{ℓ_r}Q` corresponds to `Q·Π_i([−ℓ_i] + 1)`. This is nonzero unless `Q ≡ 0`.

*Consequence.* Apply this to `σ̂ ≡ 1_{{0}} − P̂ (mod 2)`. Then for **every** integer-weight `w` that is a
measure, with no flatness assumption: `P̂ ≡ 1_{{0}} (mod 2)`. That is, `P̂(0)` is odd and `P̂(g)` is even
for `g ≠ 0`.

Theorem 1.2 then reads: `F = c` on positive measure ⇒ `P ≡ c` ⇒ `c = P̂(0)` is odd. This is also the
correct form of Remark (3). (Test (B) of `t1` checks this stronger statement for one family on `ℤ²`.)

**O2 (structure under the hypotheses of Theorem 1.2; answers Remark (2)).**

*Claim.* Assume integer weights, `w` a finite measure, and `Σ m(g)²/(1+g²) < ∞`. Then
`2ŵ = b + Σ_i n_i 1_{g_i+ℤℓ_i}`, with `b` finitely supported and finitely many `ℓ_i ≠ 0`.

*Proof.*
1. Cohen gives `2ŵ = Σ n_i 1_{g_i+L_i}`.
2. Every subgroup of `ℝ` is either cyclic (possibly `{0}`) or dense. Let `D` be the sum of the terms
   with dense `L_i`.
3. `supp D` is locally finite in `ℝ`:
   - `m` is integer-valued and `Σ m²/(1+g²) < ∞`, so `#{g ∈ supp m : |g| ≤ T} ≤ (1+T²)Σ m²/(1+g²)`.
     Hence `supp m` is locally finite.
   - Cyclic cosets are locally finite.
   - `supp D ⊂ supp(2ŵ) ∪ ∪_{cyclic}(g_i+L_i)`.
4. If `D ≠ 0`, pick an isolated point `x` of `supp D`, isolated at radius `δ`.
5. Choose `ℓ_i ∈ L_i∖{0}` for the dense `L_i`, with `|ℓ_{i+1}| > Σ_{j≤i}|ℓ_j|` and `Σ|ℓ_i| < δ`. Choose
   them from the largest down, each in a prescribed small interval; this is possible because a dense
   subgroup meets every interval.
6. Then `0 = (ΠΔ_{ℓ_i}D)(x) = ±D(x)`, a contradiction. So `D = 0`. ∎

*Consequences.*
- (a) `m` is bounded and supported in finitely many arithmetic progressions plus a finite set. `F`
  depends only on finitely many characters, so it factors through a finite-dimensional torus with a
  linear flow. "On every hull" in Theorem 1.2 is therefore a finite-dimensional statement. This also
  explains why every rank-2 example in `t2` is a sum of one-direction functions.
- (b) For the true race, `Σ_{0<γ≤T} m_γ = N(T,χ) ~ (T/2π) log T` (Riemann–von Mangoldt), whereas (a)
  forces `O(T)`.

  Hence **`w` is not a finite measure on the true hull**. Theorem 1.2 provably says nothing about the
  true race. `PROGRESS.md` §3.3 asserts this without proof; Remark (2) says it was not checked.

**O3.** Remark (1): say "complementary to `DENSITY.md` Remark 7.5" (see §2) rather than presenting
Theorem 1.2 as its extension.

**O4.** Proposition 2.1: state the normalisation `δ̂_0 = 1/2π`. Add that `w` is not a measure, and give
the jump structure of `F_M`: jump `4πM/3` at 0, a jump at every dyadic point, continuous elsewhere. This
is what §2 needs.

**O5.** Cite Cohen's theorem by theorem number. Note that `DENSITY.md` Lemma 1.1 applies to any `A` of §0.

**O6.** §5: give the grid for each number, and re-evaluate the final sequences on a fine grid (they
survive). Point out that the table values are upper bounds for `τ*_η(N)`.

**O7.** Lemma 1.1: the factor 2 is unnecessary; add `L_i ⊂ G`; note that torsion-freeness is essential,
with the `(ℤ/3)²` example.

**O8.** §3: give the range `0 ≤ h ≤ 2π` of the quadratic-variation identity. §4: add the density
condition behind "`≍ √N`".

---------------------------------------------------------------------------

## Files (all in `monotonicity/density2/referee1/`)

| Script | Output | What it tests |
|---|---|---|
| `t1_mod2_lemma.py` | `t1_mod2_lemma_output.txt` | Lemma 1.1 on `ℤ`; GF(2) test on `ℤ²` (also O1); torsion control |
| `t2_theorem12_examples.py` | `t2_theorem12_examples_output.txt` | Theorem 1.2: flat arcs on the circle; the `𝕋²` example; orbit derivative; random rank-2 search |
| `t3_prop21_period_doubling.py` | `t3_prop21_period_doubling_output.txt` | Proposition 2.1 |
| `t4_prop32_B_eta.py` | `t4_prop32_B_eta_output.txt` | Proposition 3.2; quadratic-variation identity; the periodic-`η` counterexample |
| `t5_numerics_check.py` | `t5_numerics_check_part1_output.txt` | §5: exhaustive search with grid checks |
| `t5_numerics_check.py` | `t5_numerics_check_part2_output.txt` | §5: beam search with grid checks |
| `t5b_exhaustive_N22.py` | `t5b_exhaustive_N22_output.txt` | §5: exhaustive search at `N = 22` |
