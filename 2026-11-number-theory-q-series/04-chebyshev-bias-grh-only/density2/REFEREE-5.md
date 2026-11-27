AI-produced (Anthropic's Claude), adversarial referee report, 2026-11-27.

# Referee report 5: `INVLO.md` and `PRIMECORR.md`

**Scope.** I read `monotonicity/density2/INVLO.md` and `monotonicity/density2/PRIMECORR.md` in full: all labelled
results, remarks, comparisons and summaries. Context was read only for definitions and cited inputs:
- `DENSITY.md` §1 and §3;
- `HYP.md`;
- the predecessor's Lemma 9.13, Remark 9.11 and Conjecture U (scratchpad copy `DENSITY_pred2.md`);
- `THM922_DRAFT.md` §0, §6 and §8, only to check what Corollary P6 assumes;
- `LINE1.md` Proposition 4.1, which `INVLO.md` §4 cites.

**Independence.** Before this report was written I did not read any `REFEREE-*.md`, the directories
`referee1/`–`referee4/`, or `PROGRESS.md` and `RESULTS.md`. Every attribution to "referee 4" in the two files is
treated as an unverified claim. Where such a claim mattered, I redid the computation in `referee5/`. I read the
author's scripts and their saved outputs, but did not rerun them: they overwrite their own result files. None of
their code is reused.

**Bottom line.**
- Every labelled mathematical result in both files is correct as stated: Lemmas A, D, P0, P0′ and P1′,
  Propositions B′, B, E and P5, Theorems C, P1, P2 and P7, Corollaries C′ and P3, and Corollary P6 relative to
  the draft.
- In particular, I found no vacuous or degenerate case of the hypotheses of Theorem P7, and no way to meet them
  without the conclusion (§3.3).
- One worked example is wrong. In the comparison after Theorem P7, the spectrum offered to show that "H_dec is
  strictly weaker than GLI" fails condition (D), and the next bullet reuses it (§2.1). A one-line repair restores
  both claims.
- The remaining problems are smaller:
  - a window-dependent claim in a discussion paragraph (§2.2);
  - a heuristic sentence labelled rigorous (§2.7);
  - status lines and summaries that say more, or other, than what is proved (§2.3–2.6).

## 1. Verdict table

| # | File: item | Verdict | Note |
|---|---|---|---|
| A1 | INVLO: Lemma A (needle polynomial, constants) | CORRECT | all constants re-derived; §3.1 |
| A2 | INVLO: sharpness remark (Gauss–Hermite) | CORRECT | central weight `·√N → √(π/2) = 1.2533` |
| A3 | INVLO: Lemma D (a), (b) | CORRECT | Bessel facts verified in multiprecision |
| A4 | INVLO: Proposition B′ | CORRECT | exact check: equal to order `5`, differ at `6` |
| A5 | INVLO: Proposition B | CORRECT | |
| A6 | INVLO: Hypothesis GLI(ψ); "LI ⇒ GLI" | CORRECT | |
| A7 | INVLO: Theorem C | CORRECT | trivial slip "`⌊L/2⌋ ≥ L/3`", not needed (§5 a) |
| A8 | INVLO: Corollary C′ | CORRECT | |
| A9 | INVLO: Proposition E | CORRECT | |
| A10 | INVLO: Remarks (1), (2) | CORRECT | "`N > 2·10⁷`" reproduced |
| A11 | INVLO: §3, items 1–5 | CORRECT | item 2 constructed and checked |
| A12 | INVLO: §4 (assessment of Conjecture U) | CORRECT | except one sentence of the verdict, row A13 and §2.6 |
| A13 | INVLO: header, §0, statements about numerical checks, §4 verdict | WRONG AS STATED (minor) | status line §2.3; `N ≤ 300` claim §2.4; range §2.6 |
| B1 | PRIMECORR: §1 window, explicit formula (1.1), its numerical check | CORRECT | reproduced to `2·10⁻⁸` with another window |
| B2 | PRIMECORR: Lemma P0 | CORRECT | |
| B3 | PRIMECORR: Theorem P1 (a)–(e) | CORRECT | padding confirmed exactly |
| B4 | PRIMECORR: H_corr; Theorem P2 | CORRECT | |
| B5 | PRIMECORR: Corollary P3 | CORRECT | |
| B6 | PRIMECORR: "What is new and what is not" (after P3) | WRONG AS STATED | "`(T, πT]`" depends on the window; §2.2 |
| B7 | PRIMECORR: §3 status and numerical-check paragraph | CORRECT | numbers reproduce; the window is not of P7 type (§5 h) |
| B8 | PRIMECORR: Lemma P0′ | CORRECT | range of validity and majorants checked; minor gap only in the "Alternatively" route (§5 c) |
| B9 | PRIMECORR: Setting; H_dec (Z), (N), (D); "On the prime side" | CORRECT | |
| B10 | PRIMECORR: Lemma P1′ | CORRECT | sign, positivity and prediction checked exactly |
| B11 | PRIMECORR: Theorem P7 (Cases A and B) | CORRECT | degenerate cases in §3.3 |
| B12 | PRIMECORR: comparison paragraph after P7 | WRONG AS STATED | the spectrum taken from referee 4's item A13 fails (D); §2.1 |
| B13 | PRIMECORR: remark "why plain moments do not suffice" | CORRECT | `6, 0, 460` reproduced |
| B14 | PRIMECORR: Remark P4(a) | CORRECT (as a labelled heuristic) | arithmetic of the exponents checked |
| B15 | PRIMECORR: Remark P4(b) | CORRECT WITH GAP | "cannot exclude …" is not proved; range; §2.7 |
| B16 | PRIMECORR: Proposition P5 | CORRECT | |
| B17 | PRIMECORR: Corollary P6 (assumptions from `THM922_DRAFT.md`) | CORRECT | states the draft correctly; cite 9.23(i)–(iii) (§5 d) |
| B18 | PRIMECORR: "What this shows, and what it does not"; closing Remark | CORRECT | |
| B19 | PRIMECORR: header, §0 summary | WRONG AS STATED (minor) | status line §2.3; Theorem P7 summary §2.5; Proposition E range §2.6 |
| B20 | PRIMECORR: §6 literature | not checked | §7 |

## 2. Items that are not CORRECT

### 2.1 The comparison spectrum after Theorem P7 fails (D) (row B12)

> - **For simple spectra, H_dec is strictly weaker than GLI.** Take the spectrum of `REFEREE-4b.md` A13:
>   - `Γ = {nα} ∪ ⋃_jB_j ∪ {z_x, x+z_x : x ∈ ⋃_jB_j}`; […]
>
>   Every relation through `x ∈ B_j` uses `z_x` and `x + z_x`, which lie outside `supp β_R`. So
>   `E_{a,b} = 0` and H_dec holds along `L_j`. But GLI(ψ) fails for every `ψ > 3`: the set contains
>   `(n+1)α − nα − α = 0`.

**Why it fails.**
- In Lemma P1′ and in (D) the band is `I = Γ ∩ (s_1, s_2)`. That is every element of the spectrum in the band,
  not only `B_j`.
- `{nα}` must mean the whole progression `{nα : n ≥ 1}`. Otherwise GLI could not fail through
  `(n+1)α − nα − α = 0` at arbitrarily large height.
- The band has length `s_{2,j} − s_{1,j} ≥ s_{1,j} → ∞`. So for large `j` it contains consecutive multiples
  `nα, (n+1)α`. If one of them lies in an open margin, (Z) fails. Otherwise both lie in
  `[s_{1,j}+g, s_{2,j}−g]`.
- Once `s_{1,j} − g ≥ α`, the window gives `β_R(α) = 1`.
- So `(n+1)α − nα − α = 0` is a relation of length `3`, supported at height `≤ s_{2,j} < T*_j`, with
  `κ_I = e_{(n+1)α} − e_{nα} ≠ 0`. It is realised by `4` sequences of weight `1` in `G_R G_I²`.
- Hence `E_{1,2} ≥ 4`, and (D) fails at `(a,b) = (1,2)` for every `L_j ≥ 3`.
- Relations inside the band also break (D). For instance `(n+1)α + (n−1)α − 2nα = 0` gives `E_{0,4} ≥ 1`.

**Exact check.** `referee5/r5_exact_moments.py`, part (E4), takes a truncation with `α = 1`, band
`(10.5, 21.5)` and `g = ½`. Then (Z) holds, and the result is `E_{1,2} = 240` and `E_{0,4} = 3960`. Every weight is
`≥ 0`, so the full spectrum has `E` at least this large.

**Consequences.**
- The bullet "For simple spectra, H_dec is strictly weaker than GLI" has no valid example.
- The bullet "*H_dec does not imply 5.4(1).* The spectrum of the previous item satisfies H_dec" is false as
  written, for the same reason.
- Both conclusions are still true, after the following repair.

**Repair.** Delete from the progression the multiples that lie in `⋃_j [s_{1,j} − 1, s_{2,j} + 1]`. Keep at least
`α` and `2α`, by taking `α < (s_{1,1} − 1)/2`. Then:
- `I_j = B_j`, and (Z) holds for `g_j < ½`, with `B_j` kept away from the margins as before;
- every relation with `κ_I ≠ 0` contains some `z_x` with `x ∈ B_j`, and `β_R(z_x) = 0` because
  `z_x > T*_j + 1`. So `E ≡ 0`. Part (E5) of the same script checks this on a truncation: `E_{a,b} = 0`;
- GLI(ψ) still fails for every `ψ > 3`. The gap `(s_{2,j}+1, s_{1,j+1}−1)` has length `> T*_j → ∞`, so it contains
  consecutive multiples, and `(n+1)α − nα − α = 0` occurs at arbitrarily large height;
- in the 5.4(1) bullet, replace "`nα` via `2nα`" (`2nα` may have been deleted) by "`nα` via any other multiple
  of `α` in `Γ_T`";
- TS still holds, because `x ∈ B_j` is not in the `ℚ`-span of the tail.

### 2.2 "A double zero anywhere in `(T, πT]` violates H_corr(2,T)" (row B6)

> Already at `ℓ = 2`, a double zero anywhere in `(T, πT]` violates H_corr(2,T) (referee 4,
> `referee4/r4_P5_P4.py`).

**The claim depends on the window.**
- If `γ_0` is the only multiple zero, it gives exactly `Δ_2 = 4ŵ(γ_0/T)²`, by P1(a).
- H_corr(2,T) asks for `Δ_2 < c_0²`, where `c_0 = min_{|ξ|≤1} ŵ`.
- So the claim holds if and only if `2ŵ(ξ) ≥ c_0` on `(1, π]`.
- The paper fixes only that `v ∈ C_c^∞(−½,½)` is even, nonnegative and has `∫v = 1`. That does not imply the
  inequality.

**Counterexample** (`referee5/r5_examples.py`, part (C1)). Let `v` be two bumps of half-width `0.04` centred at
`±0.45`.
- Then `c_0 = 0.811` and `ŵ(π) = 0.024`.
- `2ŵ(ξ) < c_0` for all `ξ ∈ [1.956, π]`.
- So a lone double zero in `[1.96T, πT]` does not violate H_corr(2,T).

For the standard bump `exp(−1/(1−4u²))` the claim is true: there `min_{(1,π]} 2ŵ = 1.338 > c_0 = 0.961`.

**What survives.** The conclusion drawn from it ("the effective tolerance is much smaller than `c_0^ℓ`") still
holds. For even `ℓ ≥ 4`, padding gives `Δ_ℓ ≥ ŵ(γ_0/T)²(2s_2)^{(ℓ−2)/2}` with `s_2 ≍ T log T`. That exceeds `c_0^ℓ`
unless `ŵ(γ_0/T)` is very small.

### 2.3 Status lines (rows A13, B19)

Both files say, on line 7: "Status: **not yet refereed**."

The bodies say otherwise:
- `INVLO.md` cites referee 4 (lines 327 and 399–400);
- `PRIMECORR.md` §3 says "This section was read by referee 4 (`REFEREE-4b.md`), and the corrections from that read
  are included", and it credits several items to referee 4.

The headers contradict the bodies.

### 2.4 Numerical claim for Lemma A (row A13)

> `invlo_checks.py` (1) confirms these properties numerically for `N ≤ 300`.

The script's own output for `N = 300` ends with
"`max |q|/(2(2|s|)^N) on 1<=|s|<=50 = nan`": double precision overflows. So the bound `|q(s)| ≤ 2(2|s|)^N` for
`|s| ≥ 1` is not confirmed at `N = 300`.

The property itself is true: it is Chebyshev's inequality. `referee5/r5_lemmaA.py`, part (A4), confirms it in
60-digit arithmetic for `N ≤ 400` and `1 ≤ |s| ≤ 10³`. The worst ratio is `0.25`, at `N = 1`.

### 2.5 The §0 summary of Theorem P7 claims more than is proved (row B19)

> Suppose that for infinitely many `L` some **nonempty** band at height `s`, with `s/log s = o(L)`, is free in this
> sense. Then `ν_M` has no atoms.

Theorem P7 also needs three things this sentence leaves out:
- the band shape `2s_1 ≤ s_2 ≤ 3s_1`;
- zero-free margins, condition (Z);
- freeness only relative to the zeros up to `T* = ⌊L/2⌋⁴`.

The shape is essential to the proof. Case B needs `≫ s_1` free ordinates, to make the Lemma D(b) density
`≪ (s_1/log s_1)^{1/2}`. A narrow nonempty band with one or two zeros gives no bounded density at all, since
Lemma D needs three large free coefficients. So the summary, as worded, claims more than is proved.

The parenthesis "For bounded `s` the hypothesis already implies TS" is correct even under the looser wording,
because Case A does not use the shape.

### 2.6 The range of Proposition E in two summaries (rows A13, B19, B15)

`INVLO.md` lines 408–410:

> Proposition E: for every `δ`, all but a proportion `O(δ)` of every dyadic block up to height
> `C_2^{−1}α²δ·L log L` is covered by relations of length `≤ L`.

`PRIMECORR.md` lines 506–508 say the same. Proposition E proves this only for blocks `(T, 2T]` with `T ≥ log³L`.
`INVLO.md` §0 states the range correctly. The fix is to add "from height `log³L`".

### 2.7 Remark P4(b): a consistency claim labelled rigorous (row B15)

> Unlike (a), this part is not heuristic, except its last sentence. It is a rigorous consequence of Theorem P1
> and `INVLO.md` Proposition E. […] A hypothesis "the moments of order `≤ K` are Gaussian up to a relative error
> `η`", with `K` fixed or `K = o(T/log T)`, cannot exclude the structure an atom needs, however small `η` is.

**What is rigorous.**
- By P1, a relation of length `ℓ'` contributes nothing to `Δ_ℓ` when `ℓ < ℓ'` or `ℓ ≢ ℓ' (mod 2)`.
- The weights quoted for `ℓ = ℓ'` are right. They are the multinomial `ℓ'!/Π|κ_ρ|!` times `Π a_ρ^{|κ_ρ|}`, which
  lies in `[c_0^{ℓ'}, ℓ'!]`.
- Proposition E allows relation lengths up to `L ≳ C_2T/(α²δ log T)`.

**What is not proved.** "Cannot exclude the structure an atom needs" is a consistency statement. It needs an
abstract spectrum that has no relation of length `≤ K` and does have that structure. By `INVLO.md` §4 the
structure is saturation in the range of Proposition E together with the failure of TS, TS_ℤ and FD. No such
spectrum is given.

**Fix.** Either supply such a spectrum, for instance a variant of `DENSITY.md` Example 3.3a in which every relation
has length `> K`, or move the sentence into the heuristic part. Two smaller points: add the lower end `log³L`
(§2.6), and note that the value `C_2 = 1` in the example is only illustrative, since Proposition E gives no value
for `C_2`.

## 3. What was checked in the items found correct

### 3.1 INVLO.md

- **Lemma A.**
  - Needle polynomial:
    - the Fejér identities;
    - `g(0) = N+1+F_N(π) ∈ [N+1, N+2]`;
    - `∫g² ≤ 4∫F_N² ≤ 8π(N+1)`, using `Σ(1−|k|/(N+1))² = 1 + N(2N+1)/(3(N+1)) ≤ N+1`;
    - `q ≥ 8/(3π²) = 0.2702` on `|s| ≤ 2/(N+1)`, from `sin u ≥ 2u/π` and `(N+1)/(N+2) ≥ 2/3`;
    - Chebyshev's inequality outside `[−1,1]`.
  - The rest of the proof:
    - `0.27² = 0.0729 > 1/14`;
    - the split at `|W−x| = R` gives `2·16^N(β/R)^{2N} ≤ 2·16^{−N}` and `4·32^N N!σ^{2N}/R^{2N} ≤ 4·32^{−N}`;
    - `14·8π = 351.9 < 352`;
    - `2R/(N+1) ≥ 32σ/√N` exactly when `N ≥ 1`;
    - the sub-Gaussian bound `E|V|^{2N} ≤ 2·N!(2σ²)^N`, by the Gamma integral.
  - Numerics: the exact value is `(N+1)∫q² → 4π/3`, a factor `6` below the bound (A1). The Gauss–Hermite central
    weight satisfies `·√N → 1.2533` (A6). The threshold `N = 2.02·10⁷` is reproduced (A7).
- **Lemma D.**
  - Nicholson: `y(J_0²+Y_0²)` increases to `2/π` (order `0 < ½`).
  - Weierstrass product and Rayleigh's sum `Σj_{0,k}^{−2} = 1/4`.
  - The integrals `12/(πa)` and `2√π Σ^{−1/2}`.
  - In (b) both subcases give a factor `≤ θ`. In fact `sup_{y>2}|J_0| = 0.4028 < π^{−1/2}`.
  - Monte Carlo maxima respect both bounds. With three equal coefficients the maximum density is `0.290`, against
    the bound `6/π² = 0.608` (D4–D5).
- **Proposition B′.**
  - The moment expansion, and "a sequence contributes iff `k` is a relation" (by (I2)).
  - For `j ≤ L`, every nonzero relation of length `≤ j` is supported in `S_L`.
  - Hoeffding's constant `(2|c|)²/8 = |c|²/2`.
  - `W` has density `≤ D` by independence. `D < ∞` forces `σ_Y > 0`.
  - Exact check (E1), with a free part that carries a length-6 relation: the moments agree through order 5 and
    differ at order 6.
- **Proposition B.**
  - `S_L(Γ_T) ⊂ S_L` and `σ_Y² ≤ 2V_M`.
  - The inclusion of events, and the Chebyshev factor `V/(ρ/2)² = 4V/ρ²`.
  - The limit `T → ∞`, and the LI case.
- **Theorem C.**
  - `T_L < ∞` and `S_L ⊂ Γ_{T_L}`; `B ≪ M log²T_L`.
  - Lemma D(b) data with multiplicities:
    - `a ≥ M/(2T_L)`, which needs only `T_L ≥ M/√60`;
    - `#I ≫ T_L` distinct ordinates, from Riemann–von Mangoldt and `m_γ ≪ log γ`;
    - `Σ_I ≥ (M/(2T_L))²Σm_γ`, using `m_γ² ≥ m_γ`.
  - Small `T_L`: three ordinates above `max(t_χ, M) > T_L` lie outside `S_L`.
  - The limit:
    - `t/log(t+3)` is increasing;
    - `T_L/log T_L ≤ 2L/K` and `log T_L ≤ 2 log L`;
    - subsequences along which `T_L` stays bounded are handled separately.
  - Simplicity of the zeros is never assumed.
- **Corollary C′.** `C(t)` is finite, nondecreasing and tends to `∞`. `ψ` is nondecreasing. The inequality is
  strict, so GLI(ψ) holds and Theorem C applies.
- **Proposition E.**
  - The range, and `2T ≤ T*`.
  - `#I ≥ δT/(2C_χ)`.
  - `δ` cancels in the first term, so `C_2` depends on `c, M, χ` only.
  - The second term is `≪ L^{−5/2}`.
  - The printed truncation term with `(N+1)²` is weaker than the available factor `N`, but valid.
- **§3.**
  - Item 2 was constructed for `ψ(t) = 2^t`, with `e = (1, 4, 8, 13)` and `u = (3, 37, 769, 32771)`. All `10`
    relations of length `≤ 80` among `x_1..x_4` satisfy `|k_{n*}| ≥ 2^{e_{n*}−e_{n*−1}} ≥ ψ(x_{n*})`, and the gcd
    condition behind "index 1" holds (C2).
  - Item 3 is correct.
  - "For `ψ ≤ 3` GLI holds for every set" is correct: a nonzero relation among distinct positive reals has
    length `≥ 3`.
- **§4.**
  - The uniform bound on `Σ|c_γ|²` (it uses the counting law twice).
  - The dichotomy and its window of `L`.
  - `{β, 3β, 5β}`: the shortest relation has length `4` (C3).
  - `nβ ∈ S_L` for `L ≥ (n+n')/gcd(n,n')`.
  - The models of Theorem 9.22 have every `m_n ≥ 1`, so they are saturated at `L = 3`.
  - The independent-triples example is correct.

### 3.2 PRIMECORR.md, §§1–2 and §§4–5

- **§1.**
  - The window properties, and `g(−log n) = 0` for `y > 1/T`.
  - With my own window (X1), (1.1) holds to `2·10⁻⁸`. So the author's `5.4·10⁻⁴` residual is discretisation, as
    the file says.
- **Theorem P1.**
  - The expansion converges absolutely.
  - (a) `2m(m−1)` ordered sequences per ordinate.
  - (b) The symmetry `θ ↦ θ+π`, and padding.
  - (c) The exceptional `T` form a countable set: real zeros of the entire function `v̂`.
  - (d) The weight bound `c_0^ℓ`.
  - (e) Cumulants, and the explicit formula at `y = 0` with the even `g = Tw^{∗2j}(T·)`.
  - Exactly (E2): for `{1, 2, √2}` with unit weights, `Δ_3 = 6` and `Δ_5 = 220`, although the only odd relations
    of length `≤ 5` have length `3`. Padding is real.
- **H_corr and Theorem P2.**
  - P1(d) excludes every relation below `T`.
  - The three smallest ordinates give `D`.
  - `ρ ≥ 32(2V_M)^{1/2}N^{−1/2}` and `V_M^{>T} ≪ log²T/T`.
  - `C_1` and `C_2` do not depend on the window.
- **Corollary P3.** The arithmetic of `L'_j`, and both directions of (ii).
- **Proposition P5.**
  - The construction, and Kronecker's theorem.
  - The characteristic function is a product of `J_0` factors.
  - The counting law holds with `2A`.
  - The `σ²` factor in fact lies in `[(1+δ_∞/2)^{−2}, (1−δ_∞/2)^{−2}]`.
- **Corollary P6.**
  - The draft gives:
    - Theorem 9.22: integers `m_n ≥ 1` with `m_n ≤ C(M,θ)(log(n+2))^θ`, and `F = 0` a.e. on a set `S` with
      `λ(S) ≥ 3π/(4q)`;
    - Corollary 9.23(i): the counting law with `A_c = 2C(M,θ)`;
    - Corollary 9.23(ii): `σ² < ½` for `M ≤ M_0(θ)`;
    - Corollary 9.23(iii): rescaling to weight `½` on `Ωℕ` with the same `σ²` and the same atom;
    - §8 item 9: `A_c` is astronomically large.
  - P6 restates all of this correctly. It cites only (iii), which itself relies on (i) and (ii).
  - The proof (`δ ≤ ¼` keeps `σ² < 2`) is correct.

### 3.3 PRIMECORR.md §3 (band decoupling): the requested focus

- **Lemma P0′.**
  - *Existence of the limit.*
    - Under GRH, `ψ_χ(t) ≪ √t log²t`.
    - The boundary term tends to `0`.
    - In `u = log t` the integrand is `≪ (1+u)²(|f| + |f′|)(u)`.
  - *Range of validity.*
    - The inverse transform of a convolution is `2π` times a product. So `g_ε(u) = b(u−y)e^{−ε²u²/2}`, as
      stated.
    - `|h_ε(x+iτ)| ≤ e^{τ²/(2ε²)}(|β| ∗ φ_ε)(x)`. So for each fixed `ε`, `h_ε` is in the strip class with
      Gaussian decay, and the prime sum converges absolutely.
    - `h_ε` is not even. The form `χ(n)g(log n) + χ(n)g(−log n)` handles that. Alternatively, the odd part
      cancels on both sides for real `χ`.
  - *Majorants.*
    - Zero side: for `|γ| > max(2, 2S)`, `|h_ε(γ)| ≤ ‖β‖_1(2π)^{−1/2}ε^{−1}e^{−γ²/(8ε²)}`. Here
      `ε ↦ ε^{−1}e^{−γ²/(8ε²)}` increases on `(0, |γ|/2]`, so it is `≤ e^{−γ²/8}` for `ε ≤ 1`.
    - `Ω`-side: the same bound times `log(|r|+2)`.
    - Prime side: `C(1+u)²[|b′(u−y)| + (1+u)|b(u−y)|]`, where the factor `ε²u` from the Gaussian is `≤ u`.
      There is an analogous majorant with `u + y` for `Π^−`.
    - The `ε`-limit of the prime side is the conditional limit, by the same partial summation.
  - *The identity.*
    - `Ω_β` is the Fourier transform of the `C_c^∞` function `βΩ_χ`, so it is Schwartz.
    - `Π_β^−(y) ≪ ∫_0^∞(1+u)²(|b|+|b′|)(y+u)du → 0`.
  - *Numerics* (X2): for a window of P7 type (plateau `[9.8, 19.2]`, `g = 0.8`, zero-free margins, four band zeros
    with `β = 1`), the regularised identity holds to `≤ 3.3·10⁻⁷` for `ε = 0.4`–`0.5`.
- **Setting and H_dec.**
  - Such windows exist.
  - `β_I(s_i) = β_R(s_i) = 0` by continuity. Zeros at `s_1, s_2` are invisible to both windows, but (Z) excludes
    them.
  - The prime-side forms of (Z) and (N), as one-level sums, follow from P0′ at `y = 0`.
- **Lemma P1′.**
  - `(−1)^{a+b}M_{a,b} = E_{λ_A}[G_R^aG_I^b]`.
  - The comparison of the two indicators, and positivity.
  - The decoupled prediction is `E[G_R^a]·E[(G_I^{ind})^b]`. It also counts relations inside the band
    (`κ_R = 0`, `κ_I ≠ 0`), as it must. It is the right comparison object.
  - By (Z), every weight is `1` at height `≤ T*`.
  - Exact checks (E3):
    - `E ≡ 0` when the only relation is among the R-zeros, while the plain-moment excess is `6, 0, 460`;
    - a double band zero gives `E_{0,2} = 4`;
    - a length-3 relation through the band gives `E_{1,2} = 4`;
    - enumeration of sequences agrees with the definition.
- **Theorem P7: degenerate cases, each handled by the hypotheses or the proof.**
  1. *Empty band.* Excluded by (N).
  2. *Tiny bands.* The shape forces width `≥ s_1`. `s_1 → 0` is impossible, because (N) needs a zero in
     `(0, 3s_1)`, so `s_1 ≥ γ_1/3`.
  3. *`s_1 + g ≥ s_2 − g`.* The open margins cover `(s_1, s_2)`, so (Z) and (N) cannot both hold.
  4. *`s_1 < g`.* Then `[0, s_1−g]` is empty, and (Z) forbids every zero below `s_1+g`. This is consistent.
  5. *Bounded height (Case A).*
     - Finitely many zeros lie below `3 sup s_1`, so by pigeonhole some `γ_0` lies in `I_j` for infinitely many
       `j`.
     - (D) excludes every multiset relation through `γ_0`, since each has finite length and height.
     - So `γ_0` is simple and `γ_0 ∉ ⟨Γ∖{γ_0}⟩_ℚ`. TS holds with `T = γ_0`, and `DENSITY.md` Theorem 3.1
       applies.
  6. *Zeros exactly on the margins.* At `s_1−g` or `s_2+g`, `β_R = 1`; at `s_1+g` or `s_2−g`, `β_I = 1`. Weight
     `1` either way.
  7. *Multiple zeros.* In the band they are excluded by (D) at `(0,2)`. Outside the band they are harmless: their
     relation `e_ρ − e_ρ'` has `κ_I = 0`, and Proposition B′ uses distinct ordinates with `c_γ = m_γr_γ`.
  8. *Zeros in `(T*, T*+1]` (fractional weights) and above.* (D) does not constrain them, and the proof does not
     need them: Proposition B′ is applied to `Γ_{T*}`.
  9. *`s_2 < T*`.* Automatic for large `j`.
  10. *Cancellation in (D).* None is possible: `E` is a sum of nonnegative weights.
- **Can the hypotheses hold without the conclusion?** I found no way.
  - The deduction "(D) ⟹ `S_L(Γ_{T*}) ∩ I = ∅`" is exactly the input Proposition B′ needs.
  - Case B's rate is `D·R/(N+1) ≪ (s_1/log s_1)^{1/2}(log²L/L + L^{−1/2}) → 0`, and the truncation term is
    `≪ log²N/N³`.
  - On the true `χ_{−4}` bands `(T, 2T]` (C5), `#I/T = 0.45`–`0.90`, and the Lemma D(b) bound grows like
    `(T/log T)^{1/2}`, as used.
- **Comparison paragraph, apart from §2.1.**
  - "LI with simple zeros implies H_dec": correct.
  - "GLI + simple band zeros ⟹ H_dec", hence Theorem C is a special case of P7 for simple spectra: correct. A
    multiset relation with `κ_I ≠ 0` projects to an ordinate relation that meets a simple band ordinate.
  - The multiplicity-2 example: correct.
  - "5.4(1) does not imply H_dec": correct, up to the caveats in §5 f.
  - The closing remark on plain moments: correct, by (E3).

## 4. Required fixes

1. **`PRIMECORR.md` §3, comparison.** Repair the spectrum of the bullet "For simple spectra, H_dec is strictly
   weaker than GLI". Delete the multiples of `α` lying in `⋃_j[s_{1,j}−1, s_{2,j}+1]`, keeping `α` and `2α`. In
   the bullet "H_dec does not imply 5.4(1)", replace "`nα` via `2nα`" by "via any other multiple of `α` in
   `Γ_T`" (§2.1).
2. **`PRIMECORR.md` §2, "What is new".** Restrict "a double zero anywhere in `(T, πT]` violates H_corr(2,T)" to
   windows with `2ŵ ≥ c_0` on `(1, π]`, for example the standard bump. Alternatively, replace `πT` by `ξ_wT`, where
   `ξ_w := sup{ξ : 2ŵ ≥ c_0 on [1, ξ]}` (§2.2).
3. **Status lines.** Correct "Status: not yet refereed" in both files (§2.3).
4. **`INVLO.md` line 78.** Either fix `invlo_checks.py` (1) (multiprecision, or a log-scale evaluation for
   `|s| > 1`) or restrict the sentence to the values of `N` that the script actually confirms (§2.4).
5. **`PRIMECORR.md` §0 item 2.** State the band shape `2s_1 ≤ s_2 ≤ 3s_1`, the zero-free margins (Z), and that
   freeness is tested only against zeros up to `T* = ⌊L/2⌋⁴` (§2.5).
6. **Range of Proposition E.** In `INVLO.md` §4 (verdict) and in `PRIMECORR.md` Remark P4(b), write "from height
   `log³L` up to `C_2^{−1}α²δ·L log L`" (§2.6).
7. **`PRIMECORR.md` Remark P4(b).** Either prove the "cannot exclude the structure an atom needs" sentence by a
   model, or move it into the heuristic part (§2.7).

## 5. Optional improvements

- (a) `INVLO.md` Theorem C, line 236: "`N = ⌊L/2⌋ ≥ L/3`" fails for non-integer `L ∈ (3, 4)`; for example
  `L = 3.5` gives `N = 1 < 7/6`. The bound only uses `N + 1 ≥ L/2` and `N ≥ L/4`.
- (b) `INVLO.md` §0 overstates two scales.
  - Lemma A gives scale `N^{−1/2}` only when `|x−μ| + b ≪ √N σ`.
  - Proposition B's "mass `O(L^{−1/2})`" is really `≪ D(1 + |c−1| + B)L^{−1/2}`.
- (c) `PRIMECORR.md` Lemma P0′.
  - Give a precise reference for the strip-class explicit formula; the file quotes it from memory.
  - The "Alternatively, reduce to `C_c^∞`" route estimates only the prime-side cutoff error. The zero side and
    the `Ω` side need it too. Two integrations by parts give `|h_{ε,U} − h_ε|(r) ≪ (1+r²)^{−1}·o_{U→∞}(1)`.
- (d) Corollary P6: cite Corollary 9.23(i)–(iii), not only (iii). `σ² < ½` is (ii), and the counting constant is
  (i).
- (e) In the comparison, "LI ⟹ H_dec" and "GLI ⟹ H_dec" use the fact that admissible bands contain zeros
  (Riemann–von Mangoldt). For the "abstract spectra" of the heading, this needs a density assumption.
- (f) In the 5.4(1)-plus-triples example:
  - require `u_k, v_k ≥ (u_k+v_k)/3`, so that the set stays discrete, and require independence from the `b_n`;
  - the "double every element" variant reads Devin's 5.4(1) as a condition on the set of distinct ordinates, as
    `DENSITY.md` §3.3 does. Say so.
- (g) Theorem P1(b): state the weights behind "`Δ_5 = 161/90`". With unit weights the value is `220` (E2).
- (h) `p7_checks.py` uses a plateau-free window, so `β < 1` at the band zeros, and its `b` decays fast.
  - For windows of the type in Theorem P7 (a plateau, `g < 1`), `b(10) ≈ 1.4·10⁻²` and `Π^−(6) ≈ −0.06`, and the
    partial sums approach the limit only to `~10⁻²` (X3, `r5_piminus`).
  - Lemma P0′ is unaffected, since only `Π^− → 0` is used. But "negligible for `y ≥ 6`" holds only for the
    author's window.
- (i) The comment in `invlo_checks.py` (3), "the length-3 relation 2*(5/2) − 2 − 3 = 0", is wrong: that relation
  has length `4`. The order-3 discrepancy comes from `1 + 2 − 3 = 0` and `2·1 − 2 = 0`.
- (j) `INVLO.md` cites `DENSITY.md` Lemma 9.13, Remark 9.11 and §9. These are in the predecessor branch; the
  repository's `DENSITY.md` ends at §8. The facts used (`π_T` maps Haar measure to Haar measure; the Chebyshev
  truncation) are in `DENSITY.md` (3.1) and in the proof of Theorem 3.4. Cite those.
- (k) `INVLO.md` lines 417–418, "a sufficient hypothesis new here", sits oddly with "No novelty is claimed". Write
  "not considered earlier in this repository".
- (l) "What this shows": the qualifier "stable under small perturbations" is not needed for "cannot imply 'no
  atom'", because the models themselves are counterexamples. It is needed only for "characterise atoms".
- (m) `PRIMECORR.md` §0 item 1, "vanish iff there is no relation of odd length `≤ k`": add "among the zeros with
  `a_ρ > 0`", as P1(b) does. For `T` outside a countable set, that is all zeros.
- (n) The constants of Lemma A are crude. `∫q²` is a factor `6` below its bound, and `min q ≈ 0.708` against
  `0.27`. The `2·10⁷` threshold of Remark (2) could drop by about an order of magnitude.

## 6. Scripts and key outputs (`monotonicity/density2/referee5/`)

| Script | What it checks | Key output |
|---|---|---|
| `r5_lemmaA.py` | needle polynomial: exact `∫g²`, multiprecision outside `[−1,1]`; Gauss–Hermite; threshold of Remark (2) | for `N ≤ 400`: `g(0)−(N+1) ∈ [0,1]`; `(N+1)∫q² ≤ 4.27` (limit `4π/3`) against `8π`; `min q` on `|s| ≤ 2/(N+1)` `≥ 0.644` against `0.270`; `max |q|/(2(2|s|)^N) = 0.25`; central weight `·√N` → `1.2533`; Gaussian bound `< 1` iff `N > 2.02·10⁷` |
| `r5_lemmaD.py` | the `J_0` inequalities (mpmath); Lemma D(a), D(b) against Monte Carlo maxima | `max |J_0|/√(2/(πy)) = 0.99999998`; `e^{−y²/4} − J_0 ≥ 0` on `[0, j_{0,1})`; `sup_{y>2}|J_0| = 0.4028`; three equal coefficients: `0.290 ≤ 0.608`; D(b): all `6` cases within the bound (closest: `0.176` against `0.190`) |
| `r5_exact_moments.py` | exact group-ring moments: Proposition B′, P1 padding, Lemma P1′, the comparison spectrum | (E1) moments equal to order 5, differ at 6; (E2) `Δ_3 = 6`, `Δ_5 = 220`; (E3) `E ≡ 0`, plain excess `6, 0, 460`, `E_{0,2} = 4`; (E4) spectrum as written: `E_{1,2} = 240`, `E_{0,4} = 3960`, (D) fails; (E5) repaired: `E ≡ 0` |
| `r5_explicit_formula.py` | (1.1) with an independent window; regularised P0′ for a window of P7 type; the limit `ε → 0` | (1.1) to `≤ 2.1·10⁻⁸` (`T = 3, 6`; `y = 8`–`14`); regularised identity to `≤ 3.3·10⁻⁷` (`ε = 0.4`–`0.5`); partial sums within `10⁻²`–`10⁻¹` of the limit (slow conditional convergence) |
| `r5_piminus.py` | size of `Π^−_β` | author's window: `|Π^−| ≤ 2.2·10⁻⁵` for `y ≥ 6` (claim confirmed); window of P7 type: `Π^−(6) = −0.061` |
| `r5_examples.py` | window dependence (§2.2); INVLO §3 item 2; `{β,3β,5β}`; the P4(b) number; Lemma D(b) on the `χ_{−4}` bands | two-bump window: `2ŵ < c_0` on `[1.956, π]`; standard bump fine; GLI example: no violations, index `1`; shortest relation of length `4`; `L = 1.086·10⁵`; `#I/T = 0.45`–`0.90` |

Each script writes its output to `r5_<name>_results.txt` in the same directory.

## 7. What I could not check

- **The explicit formula.** I could not check the exact hypotheses of Iwaniec–Kowalski Theorem 5.12, which the
  file quotes from memory. My checks of the explicit formula are numerical: (1.1) to `2·10⁻⁸`, the regularised P0′
  to `3·10⁻⁷`.
- **The limit `ε → 0`.** For windows of P7 type, `b` decays slowly, so I could check the conditional limit
  numerically only to about `10⁻²`.
- **The draft.** I did not check Theorem 9.22 or Corollary 9.23 themselves, only that Corollary P6 restates them
  correctly.
- **The literature.**
  - Devin's Theorem 5.4(1) was checked only through `DENSITY.md` §3.3's paraphrase.
  - Not checked: `PRIMECORR.md` §6, and the Montgomery–Soundararajan range quoted in Remark P4(a).
- **Referee-4 computations cited in the files.** Per instructions, I read none of them.
  - `r4_lemmaA.py`: its number was reproduced independently.
  - `r4_P5_P4.py`: the `(T, πT]` claim turns out to depend on the window.
  - `r4_moments_exact.py`: `Δ_5 = 161/90` cannot be checked without the weights.
  - `r4b_P0prime.py`: not reproduced as such, but replaced by my own P0′ checks.
  - `r4b_P1prime.py`: its `6, 0, 460` was reproduced.
- **The hypotheses themselves.** Whether H_dec, GLI or H_corr hold for the true zeros is open, and neither file
  claims otherwise.
- **Constants.** I did not compute explicit values of `C_χ, c_χ, t_χ, C_2`; only their existence is used.
