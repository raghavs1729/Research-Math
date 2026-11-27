# Referee 2 report on `DENSITY.md`: headlines R (time reversal) and A (atom = jump)

Independent adversarial referee, 2026-11-26. Scope as briefed: Theorems 5.1–5.4 with
`reflection_certified.py` (R), and Theorem 2.1, Corollary 2.2, Theorem 2.3 and Proposition 2.4
(A). I also checked the labelling of Proposition 7.3 and the wording of `RESULTS.md`. I
re-derived everything below by hand, or with my own scripts. I read `PROGRESS.md` §5 only
after I had reached the verdicts.

## Verdicts

| | verdict |
|---|---|
| **R** (Theorems 5.1–5.4) | **PROOF COMPLETE**, given the repository's quoted inputs (I1)–(I4) and (Q1)–(Q2), Proposition 4.4 of `GRH-ONLY.md`, and the in-repo certified zero lists. One minor slip in a side statement (finding 1): the slip is true as stated, but its proof as written covers only part of the range. |
| **A** (Theorem 2.1, Cor. 2.2, Thm 2.3, Prop. 2.4) | **PROOF COMPLETE**, given (I1)–(I4). |
| Certified numbers reproduced | **YES.** The rerun output is byte-identical. My independent recomputation agrees to at least 12 digits. |

No error was found that affects any stated theorem or certified constant.

## What I recomputed

1. **Byte-identical rerun.** I reran `reflection_certified.py` from a scratch copy, with only
   the output path changed, using python-flint 0.9.0 at 200 bits. The output is
   **byte-identical** to `reflection_certified_results.txt`. I did the same for
   `model_certified.py`: also identical. The repository was left untouched (`git status`
   clean).
2. **Independent mpmath recomputation** (`indep_JR.py`, 40 digits). It imports no repository
   module and reads only the zero lists.
   - σ²_M from the closed form `4M(Λ'/Λ)(M+½)`, with `Λ'/Λ = ½log(q/π) + ½ψ((s+1)/2) + L'/L`:
     - χ₋₄: `σ²_{1/2} = 0.155567979923586`.
     - χ₋₃: `σ²_{1/2} = 0.113229969857472`.
   - M_σ: `1.82391477877` for χ₋₄ and `2.12775095918` for χ₋₃. Both match (I4).
   - V⁺ rebuilt from the formula of `GRH-ONLY.md` Prop. 4.4, using Lemma 4.3's `J_A^+` and
     Lemma 4.2's `b_δ`, `c'`. For χ₋₄ at m = 0 this gives `0.1696308344383`; the repository
     has `0.169630834438`.
   - Theorem 5.2 (c = 1.019): `ν ≤ 0.09884975954862`, so the bound is 0.901150240.
   - Theorem 5.3 (N = 20), written from the formula in the text:
     - χ₋₄: `ν ≤ 0.08907513807014`, so the bound is **0.91092486193**.
     - χ₋₃: `ν ≤ 0.06523860156028`, so the bound is **0.93476139844**.

   All other weights in the table of §5 also agree.
3. **Closed form of ∫κ²Y_h² and ∫κ².** I checked the closed forms used in Theorem 5.3
   (items 1 and 3) by exact grid quadrature: a 16^N grid for N = 1, 2, 3, five random
   parameter sets each. The maximum deviation is 2.2·10⁻¹⁶. The repository has a grid
   self-test only for Theorem 6.1, not for `J_R`, so this check is new.
4. **(NR1)–(NR3) for N = 20, independently.** My own exact-integer code (ordinates ×10²⁴,
   coordinatewise enumeration of the 118 720 vectors with `|n|₁ ≤ 4`) gives these margins:
   - χ₋₄: NR1 2.765974e-3, NR2 2.375512e-4, NR3 1.476929e-6. These are identical to
     `nonresonance_certified_results.txt`.
   - χ₋₃: smallest margin NR3 = 7.941163e-7. This matches `GRH-ONLY.md` §6.

   All margins are far above the threshold 4·10⁻²⁰.
5. **Non-LI consistency test** (`checks.py`). I used the Proposition 7.3 model: integer
   frequencies, a cyclic hull, LI as false as it can be, and a real atom at 0.
   - The identity `F(t) + F(−t) − 2 − R(t)` vanishes to 1.4·10⁻¹⁴, with `max|R| = σ²`.
   - `P(F ≤ 0) = 0.110`. This is below the Theorem 5.2 bound with optimised c (0.375), below
     the c = 1 form (0.485), and below 1/2.
   - The one-sided reflection inequalities of Theorem 5.1(b) hold on 81 levels.

All scripts and outputs are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/densref/`:
`indep_JR.py`, `indep_JR_out.txt`, `checks.py`, `checks_out.txt`, `rerun_reflection.py`,
`rerun_reflection_results.txt`, `rerun_model.py`, `rerun_model_results.txt`.

## Hand checks (all passed)

### R(a): the identity and the bound |R| ≤ σ²

- **Identity.** `F_M = 1 + Re Σ_{γ∈Γ} m_γ r_γ e_γ`, with `r_γ = 4M/(M+iγ)`. This gives
  `Re(r e^{iφ}) = (w/2)cos φ + s sin φ`, where `w = 2Re r = 8M²/(M²+γ²)` and
  `s = −Im r = 4Mγ/(M²+γ²)`.
- Since `e_γ∘inv = ē_γ`, we get `F∘inv = 1 + Re Σ m_γ r̄_γ e_γ`. Hence
  `F + F∘inv = 2 + Σ m_γ w_γ cos φ_γ = 2 + R`. This is an identity in L²(A), so it holds
  a.e.
- **Haar measure.** `inv` preserves λ_A: A is a compact abelian group, the closure of the
  subgroup ι(ℝ), and Haar measure is inversion-invariant.
- **Multiplicities.** The coefficient of R carries `m_γ` to the **first** power. So
  `|R| ≤ Σ m_γ w_γ = σ²_M`, which is exactly the Hadamard closed form, counted with
  multiplicity.
- **Convergence.** The series for R converges absolutely and uniformly
  (`Σ m_γ/γ² < ∞`), so R is continuous and the bound holds everywhere on A. The variance
  correctly uses `m_γ²`: `V = V^e + V^o`, with `V^e = ⅛Σ m²w²` and `V^o = ½Σ m²s²`.
- **Pointwise statement.** F_M itself is only an L² class, because `Σ m_γ|r_γ| = ∞`.
  Theorem 5.1(a) correctly says "a.e."; see finding 6.

### R(b): from the hull to the race

- **Limiting law.** ν_M is the law of F_M under λ_A by `GRH-ONLY.md` Lemma 1.2 and
  Theorem 1.3. Kronecker–Weyl on the orbit closure needs no LI: the proof uses only that a
  character `ψ_n` is trivial on A iff `ω_n = 0`. Together with the explicit formula (Q1) and
  the mean-value Lemma 1.1, this gives W_p convergence.
- **Direction of the portmanteau step.** Correct. The repository uses `limsup` over the
  closed set `(−∞,0]` (GRH-ONLY Thm 5.1). Equivalently,
  `\underline d{D_m>0} ≥ ν((0,∞)) = 1 − ν((−∞,0])` (DENSITY (I3)/(1.1)). Moving between
  `dt/t` on [2,X] and `dy` on [0, log X] costs O(1/Y).
- **Citations.** The Rubinstein–Sarnak, Devin and Bailleul references are background only;
  the proof is self-contained given (Q1).

### R(c): the inequalities

- **5.1(b).** `X' ≤ t ⇒ X ≥ 2 − σ² − t`. The W_∞ reformulation is a correct quantile-coupling
  equivalence.
- **5.1(c).** `{F≤0} = {Y ≤ −(1+R/2)}` and `{F∘inv≤0} = {Y ≥ 1+R/2}`. These are disjoint
  because `1 + R/2 ≥ 1 − σ²/2 > 0`, and they have equal measure.
- **5.2.** This is a Chebyshev bound on the two disjoint sets N and N′, with
  `λ(N′) ≥ λ(N)`. It is correct for `0 < c < 2 − σ²`; see finding 1.
- **5.3, the multiplier.** `ρ = ½Σ_{∉head} m w` uses simple head zeros. Then
  `1+R/2 ≥ (1−ρ)(1+x)`, `1/(1+x) ≤ 1 − x + e`, `κ ≥ 1/(1+R/2) > 0`, and Markov's inequality
  gives `λ{κ|Y| ≥ 1} ≤ ∫κ²Y²`.
- **5.3, the integrals.** Items 1–5 check out.
  - Items 1–3 are justified by (NR1)–(NR3) as written in `GRH-ONLY.md` §6. In item 2, the
    monomials of `κ²Y_h` have odd degree, so n ≠ 0. In item 3, the diagonal terms
    `cos 2φ_γ` are excluded by size.
  - Item 4 uses the buffer `T − T' = 2γ_N`.
  - Item 5 uses `‖κ‖_∞ ≤ a₀ + Σk_n` and `∫Y_t² = ½Σ_{γ>T'} m²s²`.
- **5.3, monotonicity.** `J_R` is increasing in V and decreasing in V^e, as stated.

### R(d): `reflection_certified.py`, line by line

- **No floating point on the certified path.**
  - Floats are used only to choose c (`thm52`) and as the starting guess in
    `certified_floor`. In both places the bound itself, and the final comparison, are done
    in Arb.
  - The set K is decided by certified Arb comparisons. An undecided zero would be left in
    the tail, which only weakens the bound and stays covered by NR3, because NR3 covers all
    pairs in (γ_N, T].
- **Zeros.** The zeros are balls of radius 10⁻²⁰ around the stored decimals, certified by
  `certify_arb.py` [C]/[D]. These are sign changes of Hardy's Z together with a winding
  number count. The count covers heights up to `γ_last + 0.25`, which is exactly the
  `T = g[-1] + 0.25` the script uses (640.1657 and 700.0247). `cert-audit/AUDIT.md` reports
  a rerun of this count (511/511, 537/537); I did not rerun it.
- **Tail beyond the last zero.** It enters only through ρ, via the closed-form ball for σ²
  (radius doubled after decimal parsing), and through V⁺ from Prop. 4.4, read from the
  printed ball and rounded **up**, plus 10⁻¹⁰.
- **Lower bound for V^e.** `V^e_low` is the sum over the 511 (537) listed simple zeros,
  which is a valid lower bound.
- **Non-resonance inputs.** These are the repository's certified ones: `exact_nr()`,
  imported unchanged from `grhonly/nonresonance_certified.py`. They are checked for the
  same head N and the same T, and gated by `3γ_N < T`.

### A: the atom equals the jump

- **(2.1).** This is the triangle inequality plus the definition of λ_χ.
- **(2.2).** This is a comparison of Fourier coefficients:
  `4/(M+iγ) − 4/(M'+iγ) = 4κ/((M+iγ)(M'+iγ))`, and constants `κ/(MM')`.
- **Positivity on I = (0, M_σ].** `λ_χ > 0` on I by Cauchy–Schwarz:
  `Σ4/(|M+iγ||M'+iγ|) ≤ σ_Mσ_{M'}/(2MM')`, and σ_M is strictly increasing in M.
- **Citation.** The prime-side citation (`monotone-coupling.tex` Theorem A(ii), with range
  `λ_χ > 0`, and (iii) `σσ' < 2`) is accurate, as is `GRH-ONLY.md` Cor. 2.4 / Thm 7.1(a).
  However, the proofs do not use it: everything goes through the hull form (2.1)–(2.2)
  together with (I1)–(I4).
- **Theorem 2.1(a).** `{F_{M'} ≥ 0} ⊂ {F_M > 0}` a.e. Hence `\bar d(M') ≤ d(M)`.
- **Theorem 2.1(b).**
  - d is lower semicontinuous and `\bar d` is upper semicontinuous, by L²-continuity of M ↦
    F_M and the portmanteau theorem.
  - Squeezing gives `d(M₀⁺) = \bar d(M₀⁺) = d(M₀)` and `\bar d(M₀⁻) = d(M₀⁻) = \bar d(M₀)`.
  - Therefore `ν_{M₀}({0}) = \bar d(M₀) − d(M₀)` equals the jump. Both inequalities hold,
    because this is an equality.
  - Each of the four true densities lies between d and `\bar d` by (2.3), so each has the
    same one-sided limits and the same jump.
- **Theorem 2.1(c).** The tie sets `{F_M = 0}`, M ∈ [M₁,M₂], are pairwise disjoint modulo
  null sets. They all lie inside `{F_{M₁} ≥ 0} ∖ {F_{M₂} > 0}`, which has measure
  `\bar d(M₁) − d(M₂)`. So the atoms are countable and their total mass is at most 1, and
  in fact `≤ 1 − d(M_σ) < 1`.
- **Theorem 2.1(d).** No atom gives `d = \bar d`, and the densities exist by (2.3). This is
  the standard direction. Near-ties are handled by (1.1).
- **Direction of the stated equivalences.** The text claims "no atom ⟺ δ continuous at m",
  which is proved because it is an equality. It claims "δ continuous ⇒ density exists",
  which is proved. The converse "density exists ⇒ no atom" is explicitly **not** claimed:
  see after Cor. 2.2, in §8, and in RESULTS item 3. This is correct.
- **Range.** `m ≤ 1.3239` (χ₋₄) and `m ≤ 1.6277` (χ₋₃) lie inside `M_σ − ½`, by my M_σ
  values above.
- **Theorem 2.3 and Proposition 2.4.** I re-derived the constants a_i, b_i. The cocycle
  version of Prop. 2.4 is right, including the symmetric `G_{M,M'}` on both sides of M₁.

### Labelling of Proposition 7.3

Proposition 7.3 is stated for "the exact analogue of F_{1/2} with the ordinates replaced by
the positive integers". §8 lists it under "statements about model functions on the circle".
It is correctly labelled as a statement about **models**, not about the true race.

## Findings, by severity (none fatal)

1. **(Minor; mathematical, repairable in one line.)** The "in particular" bound of
   Theorem 5.2 (DENSITY.md l.706) has a gap in its proof.
   - **The gap.** The bound `ν_M((−∞,0]) ≤ V_M/(1+(1−σ_M²)²)` is obtained by taking c = 1.
     That choice is admissible only when `1 < 2 − σ²`, i.e. **σ² < 1**. The theorem,
     however, assumes only σ² < 2. For `1 ≤ σ² < 2` the bound is still true, but for a
     reason not given: `V ≥ σ² ≥ 1 > (1−σ²)²` makes it weaker than Cantelli's `V/(1+V)`.
   - **Unqualified restatements.** The same bound is repeated with no range at all in §0
     row R (l.26), the RESULTS.md table (l.50) and Remark 7.4 (l.964). For σ² ≥ 2 it is not
     proved, and for large σ² it is not implied by Cantelli.
   - **Fix.** Add "(for σ² < 1; for 1 ≤ σ² < 2 it follows from Cantelli)" and state σ² < 2
     wherever the bound is quoted.
   - **Impact.** None on the certified numbers: they use c = 1.019 at σ² ≈ 0.156.
2. **(Minor; overclaim in RESULTS.md, item 6, ll.38–40.)** RESULTS.md says: "A proof must
   use integrality of the multiplicities or arithmetic of the dense orbit."
   - DENSITY Remark 7.4 and §8 item 5 say less. A proof must use one of five missing
     features (α)–(ε): integer weights, the explicit formula, weights ≥ 1 on every
     frequency, a non-cyclic G, and the growth of N(T) (unbounded F). DENSITY only *guesses*
     (α)/(β).
   - This is exactly the overclaim the second internal read flagged (PROGRESS §5, second
     table, item 1). It was fixed in DENSITY but not carried over to RESULTS.md.
   - Also, "Such a model has the hull" should read "a hull (the circle)".
3. **(Cosmetic.)** In the table after Theorem 5.4, at χ₋₄, m = −0.25, tier A, Thm 5.3's
   0.97853 beats 6.2's 0.9784, but neither entry is bold ("bold = best known").
4. **(Cosmetic.)** In RESULTS item 2, "(was 0.1157)" is the previous value for q = 4 only.
   For q = 3 the previous value was 0.0888.
5. **(Remark, not an error.)** Two of the certified consequences are immediate rather than
   separate results.
   - The bound `ν_{1/2}({0}) < 0.0891` (resp. 0.0653) follows at once from the lower-density
     bound.
   - The same goes for "upper minus lower density < 0.0891", since
     `\bar d − \underline d ≤ 1 − \underline d`.
6. **(Cosmetic.)** §0 row R and RESULTS item 2 write `F(θ)+F(θ⁻¹) = 2 + R(θ)` without "a.e.".
   F_M is only an L² class. Theorem 5.1(a) itself is correct.
7. **(Wording.)** RESULTS item 1 says the problem "is Problem 13 … in the form 'ties have
   density 0'". Per LITERATURE.md, Problem 13 asks for an *unconditional* proof that ties
   have density 0. The note's target is GRH-conditional and stronger (no atom ⇒ near-ties
   density 0). This is loose wording, not a claim of progress. "Nothing here concerns
   Problem 12" is consistent throughout.
8. **(Outside the headline, noted for completeness.)** For tier B, `read_Vplus` adds 2·10⁻¹⁰
   to a *printed* 12-digit value rather than parsing a ball. This relies on a comment that
   the radius is below 10⁻¹². Tier-B numbers are not headline results and do not improve
   the repository's values anyway.

## What I could not check

- **Inputs I trusted without re-deriving.**
  - The explicit formula with error term (Q1) from `weighted-races.tex`, and the mean-value
    Lemma 1.1 of `GRH-ONLY.md`.
  - The Hadamard/multiplicity Lemma 4.2. I checked only that my implementation of
    Prop. 4.4 and Lemma 4.3 reproduces V⁺, and re-derived Lemma 4.3's integrals.
- **Not rerun.**
  - The zero localisation and completeness run of `certify_arb.py` (winding number). I rely
    on its recorded audit rerun in `cert-audit/AUDIT.md`.
  - Arb itself is trusted.
- **Novelty claims.** I did not check them, including that the time-reversal identity is
  absent from the literature and the Kaczorowski–Ramaré attribution.
- **Sections outside the brief.** Not reviewed: §3 (TS/FD), §4, §6, Propositions 7.1–7.2,
  and the proof details of Proposition 7.3 beyond the byte-identical rerun of
  `model_certified.py`.
