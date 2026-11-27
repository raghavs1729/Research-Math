AI-produced (Anthropic's Claude), adversarial referee report, 2026-11-27.

# Referee report 4: `INVLO.md`, `PRIMECORR.md`, `LINE1.md` §4

## 0. Scope, method, summary

**Refereed (read in full):**
- `INVLO.md`, all of it;
- `PRIMECORR.md`, all of it (§6, the literature check, was read but not checked);
- `LINE1.md` §4.

**Context read:**
- `density/DENSITY.md` §1 and §3;
- `density2/HYP.md`;
- the predecessor's Lemma 9.13, the class `𝒞` and Conjecture U (scratchpad copy `DENSITY_pred2.md`);
- the statements of Theorem 9.22 and Corollary 9.23, and §6 and §8 of `THM922_DRAFT.md`, only to check what Corollary P6 assumes.

As instructed, I did not open `PROGRESS.md`, `RESULTS.md`, any `REFEREE-*.md`, or `referee1/`–`referee3/`.

I checked every step by hand. The numerics are in `referee4/` (§5). Choices I made:
- `L` is an integer wherever `⌊L/2⌋ ≥ L/3` is used. The claim fails for real `L ∈ (3,4)`, which is irrelevant.
- In `INVLO.md` §3, TS, FD and TS_ℤ are read verbatim as conditions on abstract frequency sets.

**Summary.**
- *No error in any labelled proof of `INVLO.md`.* All constants of Lemma A and Lemma D hold, as do the two Bessel facts, the moment identity of Proposition B′, the uniformity in `T` in Proposition B, the zero counting and the limits in Theorem C, the diagonal argument of Corollary C′, Proposition E, and both abstract examples of §3.
- *No error in any labelled proof of `PRIMECORR.md`*, apart from two small gaps (P1(b), P3(i)). The explicit-formula identity (1.1) and P1(e) hold numerically to `10⁻⁸` and `10⁻¹²` with a window different from the author's.
- *What is wrong:*
  1. a factor-2 notational error in `r_γ` (harmless; `PRIMECORR.md` imports it);
  2. a false sentence in `INVLO.md` §4: "cyclic models are saturated as soon as `L ≥ 3`";
  3. a sign error in `PRIMECORR.md` §0 item 1;
  4. a computation in Remark P4(b), and its §0 summary, that does not follow from Theorem P1;
  5. a logical overreach in the barrier conclusion (§0 item 4; the paragraph after Corollary P6).
- *Gaps (supplyable):* P1(b) needs the padding argument; P3(i) needs a line about the `M`-dependence of `C_1`; two range statements in `INVLO.md` §§3–4.
- *`LINE1.md` §4:*
  - Proposition 4.1 is correct.
  - Lemma 4.2 is honestly labelled a sketch, and its statement is plausible. The statement needs two precisions, and the sketch omits that the prime side converges only conditionally.
  - "The gap" uses counts of *distinct* zeros that are not known.

## 1. Verdict table

| # | Item (result and step) | Verdict |
|---|---|---|
| **INVLO.md** | | |
| I1 | Notation: `r_γ = 2M/(M+iγ)` (l. 36), cited as `DENSITY.md` §1 | WRONG AS STATED (harmless) |
| I2 | Lemma A: needle polynomial (degree `≤ N`, `g(0) ∈ [N+1,N+2]`, `0 ≤ q ≤ 2`, `∫q² ≤ 8π/(N+1)`, `q ≥ 0.27` on `\|s\| ≤ 2/(N+1)`, Chebyshev bound outside `[−1,1]`) | CORRECT |
| I3 | Lemma A: the bound (split at `R`, the two tail parts, `14·8π < 352`, `56·16^{−N}`); sub-Gaussian ⇒ moment condition | CORRECT |
| I4 | Lemma A: sharpness remark (Gauss–Hermite atom `≈ 1.24/√N`) | CORRECT |
| I5 | Lemma D: the two Bessel facts (Nicholson; product, Rayleigh) | CORRECT |
| I6 | Lemma D(a) and D(b) | CORRECT |
| I7 | Prop. B′: moment identity (frequency relations, `j ≤ L`), application of Lemma A (Hoeffding, density of `W`) | CORRECT |
| I8 | Prop. B: truncation, `S_L(Γ_T) ⊂ S_L`, `R` and `D` independent of `T`, limit `T → ∞`; the "in particular" | CORRECT |
| I9 | Hypothesis GLI(ψ) | CORRECT ("LI is GLI(∞)" is loose wording) |
| I10 | Theorem C (`S_L` low; mass `B`; zero-counting facts, which are unconditional; Lemma D(b)/(a); the limit in `K`) | CORRECT |
| I11 | Corollary C′ (diagonal argument) | CORRECT |
| I12 | Proposition E (proof) | CORRECT |
| I13 | Prop. E as summarised in the header (l. 26–28) and after the statement (l. 281–284) | WRONG AS STATED (minor overstatement: δ) |
| I14 | Remarks (1), (2) | CORRECT (one misattribution in (1)) |
| I15 | §3 item 1; item 2 (rank-one example: GLI(ψ) holds, TS, FD, TS_ℤ fail); item 3 (TS without GLI); item 5 | CORRECT |
| I16 | §3 item 4 ("incomparable"; "the disjunction … implies continuity") | CORRECT WITH GAP (hypothesis on ψ omitted) |
| I17 | §4, structured/unstructured split | CORRECT WITH GAP (range of `L`) |
| I18 | §4, "Cyclic models (`Γ' ⊂ βℕ`) are saturated as soon as `L ≥ 3`" | WRONG AS STATED |
| I19 | §4, verdict | CORRECT (modulo I17, I18) |
| **PRIMECORR.md** | | |
| P1 | §0 item 1: "moment = prediction **plus** a nonnegative sum" | WRONG AS STATED (sign for odd `k`) |
| P2 | §1: window `w = v∗v`, `ŵ ≥ 0`, `ŵ > 0` on `(−π,π)`, `c_0 > 0` | CORRECT |
| P3 | §1: explicit formula and (1.1) (normalisation, sign, admissible `g`, convergence, `Ω_T → 0`) | CORRECT |
| P4 | Lemma P0 | CORRECT |
| P5 | Theorem P1: identity `(−1)^kM_k = μ_k + Δ_k`, positivity, length `≤ k`, `≡ k (mod 2)` | CORRECT |
| P6 | Theorem P1 (a), (c), (d), (e) | CORRECT |
| P7 | Theorem P1 (b): the "iff" | CORRECT WITH GAP (padding not given) |
| P8 | Hypothesis H_corr; Theorem P2 (threshold `c_0^ℓ`, the three steps, `T ≥ C_1L^{3/2}log²L`) | CORRECT |
| P9 | Corollary P3 (i) | CORRECT WITH GAP (`C_1 = C_1(M)`) |
| P10 | Corollary P3 (ii) | CORRECT |
| P11 | "What is new and what is not" | CORRECT (wording, §4 of this report) |
| P12 | §3 (pairings amplify; mixed moments) | CORRECT (informal; padding bound checked) |
| P13 | Remark P4(a) (heuristic) | CORRECT WITH GAP (as a heuristic) |
| P14 | Remark P4(b), and §0 item 3 second bullet | WRONG AS STATED |
| P15 | Proposition P5 (construction; a.c. law of the infinite independent sum; counting law with `2A`; `σ²` factor) | CORRECT |
| P16 | Corollary P6: the statement, and what it assumes from the draft | CORRECT (conditional on the draft) |
| P17 | §0 item 4, last sentence; the paragraph "Hence …" after P6; the Remark after P6 | WRONG AS STATED (non sequitur) |
| **LINE1.md §4** | | |
| L1 | Proposition 4.1 | CORRECT |
| L2 | Lemma 4.2: statement and the "sketch" label | CORRECT WITH GAP |
| L3 | "The gap" and "What would revive the route" | CORRECT WITH GAP (two counting claims WRONG AS STATED) |

## 2. Findings for the items that are not CORRECT

### I1. Normalisation of `r_γ` (INVLO, l. 36; imported by PRIMECORR P2)

> "`r_γ = r_γ(M) = 2M/(M+iγ)` and `c_γ := m_γ r_γ`; `F_M = 1 + Re Σ_{γ∈Γ} c_γ e_γ` on the hull `A`, and `ν_M` is its law"

`DENSITY.md` §1 (l. 40) has `r_γ(M) = 4M/(M+iγ)`, and so do `HYP.md` and the predecessor. With INVLO's `r_γ`, `F_M − 1` is half of the true one. Its law is then the image of the true `ν_M` under `x ↦ (1+x)/2`. So the sentence "Under GRH + NCZ, `ν_M` is the law of `F_M` (`DENSITY.md` §1)" in the proof of Proposition B is false as written.

Nothing breaks:
- atoms correspond, `c ↔ 2c − 1`;
- every statement is for all `c`, with unspecified constants;
- the explicit bounds stay true with `4M`, only with other constants: `a ≥ M/(2T_L)` still holds, and `b ≤ 2MC_χ log(2T_L+3)/T_L` becomes `4MC_χ …`.

Fix: write `4M/(M+iγ)` and double the constant in `b`.

### I13. Proposition E's summaries drop `δ`

> Header: "For every large `L`, almost all ordinates in every dyadic block from height `log³L` up to `≍ α²L log L` occur in integer relations of length `≤ L`."

What is proved: for each `δ`, all `L ≥ L_0(δ)` and all `T ∈ [log³L, C_2^{−1}α²δ L log L]`, the free mass is `< δT log T`, i.e. a free proportion `< 2πδ(1+o(1))`.
- The upper end of the range shrinks linearly with `δ`. At heights `≍ α²L log L` with a fixed constant, only a fixed, non-small proportion is controlled.
- For `δ ≥ (1+o(1))/(2π)` the statement is trivial, since `Σ_{(T,2T]}m_γ ~ (T/2π)log T`.

The same applies to l. 281–284 ("This holds for every dyadic block between `log³L` and `≍ α²L log L`") and to §4 ("saturation of every dyadic block up to `≍ α²L log L`").

Fix: "all but a proportion `O(δ)` … up to `C_2^{−1}α²δ L log L`".

### I16. §3 item 4: the hypothesis on ψ is missing

> "So GLI is incomparable with each of TS, FD and TS_ℤ, and the disjunction "TS_ℤ or FD or GLI(ψ)" implies continuity of `ν_M`"

- Continuity holds only for ψ with `ψ(t)log(t+3)/t → ∞` (Theorem C).
- Incomparability needs `ψ > 3` eventually. For `ψ ≤ 3`, e.g. `ψ ≡ 0`, GLI(ψ) holds for every set, so it is implied by TS.

Fix: insert both hypotheses.

### I17. §4, the split: the range of `L` is missing

> "So once `N` is large compared with `(1 + Σ|c_γ|²)/(ηa)²`, the second case forces `B ≳ ηaN` … At every length `L` up to `≍ 1/ε`, either the short relations absorb all but two of the frequencies with coefficient `≥ a`, or they involve frequencies of `ℓ¹`-mass `≳ ηaL`."

The alternative "`B ≳ ηaL`" is derived only in the window
`C(1+Σ|c_γ|²)/(ηa)² ≤ L ≲ 64/ε`, with `56·16^{−⌊L/2⌋} ≤ η/2`.

Below that window, the second case of the dichotomy reads `1 + B + 2√Nσ_Y ≳ ηaN`. This can hold with `B = 0`, i.e. with `S_L = ∅`. So "at every length `L` up to `≍ 1/ε` … `≳ ηaL`" is not what is shown.

The window is nonempty uniformly over `𝒞`, which the text uses but does not state. Indeed `Σ|c_γ|² = Σ m_γ²|r_γ|² ≤ A_0 Σ m_γ log(γ+3)|r_γ|² ≤ C(M, A_0)` on `𝒞`, by the counting law.

Fix: state the window, and the uniform bound on `Σ|c_γ|²`.

### I18. §4: cyclic models are not saturated at `L = 3` in general

> "Cyclic models (`Γ' ⊂ βℕ`) are saturated as soon as `L ≥ 3` (`β + nβ − (n+1)β = 0`)."

The relation used needs `β, nβ, (n+1)β ∈ Γ'`.

**Counterexample.** `Γ' = {β, 3β, 5β}` has no relation of length `≤ 3`, so `S_3 = ∅`, and its shortest relations have length 4. `r4_moments_exact.py` (B2) checks this exhaustively, and confirms that the moments then agree up to order 3 only.

**A cyclic model in `𝒞`** (`r4_moments_exact.py` (B4)). Take `Γ' = 10⁻⁴·{60209, 102438, 129881}` (`≈ 6.02, 10.24, 12.99`) with `m = 1`.
- At `M = ½`, `σ² = 0.0856 ≤ 0.1556`, and the counting law holds trivially, so the model is in `𝒞`.
- Its shortest relation, `(−22,−53,52)`, has length 127. So `S_L = ∅` for `L ≤ 126`, and Proposition B′ applies to the model at those lengths.
- A smaller `β` pushes this length further.

(A Sidon set is not enough: `{1,3,4}` is Sidon but has `1 + 3 − 4 = 0`; see (B5).)

**What is true.**
- For `Γ' = β{1,…,N}` (all multiplicities positive, as in Theorem 9.22 of the draft), `S_3 = Γ'`; (B3) checks this.
- For any `Γ' ⊂ βℕ` with at least two elements, `nβ` lies in `S_L` as soon as `L ≥ (n+n')/gcd(n,n')` for some other `n'β ∈ Γ'`.
- So saturation holds at lengths `≳ 2·max{n : nβ` has coefficient `≥ a}`. This is uniform when `β` is bounded below, e.g. in the cyclic case of Question P″ with `Γ ⊂ ℕ`, but not uniformly over `𝒞`.

The conclusion that the moment method says nothing on the cyclic case of P″ (with `L ≍ 1/ε → ∞`) survives.

### P1. PRIMECORR §0 item 1: sign

> "The log-averaged `k`-th moment of the windowed prime sum … equals the *independent-phase prediction* plus a **nonnegative** sum over the integer relations of length `≤ k`"

Theorem P1 gives `(−1)^k M_k = μ_k + Δ_k`. For odd `k`, the prediction is `0` and `M_k = −Δ_k ≤ 0`, so the correction is non*positive*. The next bullet ("odd log-averaged moments are `≤ 0`") contradicts the sentence as written.

Fix: "`(−1)^k` times the moment equals …".

### P7. Theorem P1(b): the "iff" needs padding

> "*(b)* `μ_k = 0` by the symmetry `θ ↦ θ + π`. Then apply (c)."

(c) characterises "`Δ_ℓ = 0` for **all** `ℓ ≤ L`". The claim in (b) is about a single odd `k`: `Δ_k = 0` iff there is no odd relation of length `≤ k`.

For "⇒", a relation `κ` of odd length `ℓ < k` must force `Δ_k > 0`. This needs the padding argument of §3: append `(k−ℓ)/2` pairs `(ρ,+),(ρ,−)` with `a_ρ > 0`, which leaves the net vector unchanged.

**Example** (`r4_moments_exact.py`, model 2). Take ordinates `{1, 2, √2}`. The only odd relations of length `≤ 5` are `±(2,−1,0)`, of length 3, and yet `Δ_5 = 161/90 > 0`.

The statement of (b) is correct; the proof needs that line.

### P9. Corollary P3(i): `C_1` depends on `M`

> "(i) Suppose H_corr(`L_j, T_j`) holds along `L_j → ∞`, with `T_j ≥ C_1L_j^{3/2}log²(L_j+3)`. Then `ν_M` has no atoms, and the conclusions of `HYP.md` Theorem 1.3 hold."

In Theorem P2, `C_1 = C_1(M,χ)`. Its step 3 uses `V_M^{>T} ≤ C M² log²T/T`, valid for `T ≥ M`, and the ratio `M²/V_M`. The conclusions of HYP Theorem 1.3 need no atoms for **every** `M = m + ½`.

**Missing argument.** H_corr(`L,T`) implies H_corr(`L',T`) for `L' ≤ L`, and H_corr does not depend on `M`. Given `M'`, put `L'_j := min(L_j, ⌊(C_1(M)/C_1(M'))^{2/3}L_j⌋)`. Then `T_j ≥ C_1(M')L'^{3/2}_j log²(L'_j+3)`, so P2 at `M'` gives `ν_{M'}({c}) ≤ C_2(M')(1+|c−1|)L'^{−1/2}_j → 0`.

### P13. Remark P4(a) (heuristic): three caveats

The order-of-magnitude computation is right: `r4_P5_P4.py` (P4a) evaluates the ratio. The caveats:
- *`o` versus `O`.* "Errors `O_T(X^{1−k/2})` … bounded for `k = 2`, and `o(X^{−1/2})` for `k ≥ 3`." The needed size is `X^{1−k/2}(T log T)^{k/2}`, i.e. `O(1)`, `O(X^{−1/2})`, `O(X^{−1})`, …
- *An unargued identification.* "By Theorem P1 that average is the relation structure of the zeros below `≍ T`" presupposes that the shift-average of the Hardy–Littlewood main terms equals `μ_k(T)` to the precision of `Δ_k`. P1 says nothing about Hardy–Littlewood terms.
- *Different averages and windows.* P1 concerns logarithmic averages with smooth windows, while (a) uses natural averages with sharp windows. At fixed `T`, the natural-average moments of `S(x)/√x` do not converge: they are almost periodic in `log X`. So only orders of magnitude are compared.

The heuristic label is honest.

### P14. Remark P4(b) and §0 item 3: the count does not follow from Theorem P1

> "This needs only `≍ T log T/L` relations among the `N ≍ T log T` zeros up to `T`. Their contribution to `Δ_ℓ(T)` is `≍ T log T/L`, against `μ_ℓ(T) ≍ (T log T)^{ℓ/2}`: a relative excess `≍ (T log T)^{1−ℓ/2}/L`. … excludes the structure an atom needs only if `η ≪ (T log T)^{1−ℓ/2}/L`."

The count assumes that each relation contributes `≍ 1` to `Δ_ℓ`. By Theorem P1, a relation `κ` of length `ℓ'` contributes:
- `0` to `Δ_ℓ` if `ℓ < ℓ'` or `ℓ ≢ ℓ' (mod 2)`;
- `(ℓ'!/Π|κ_ρ|!)·Πa_ρ^{|κ_ρ|}` at `ℓ = ℓ'`, which lies between `c_0^{ℓ'}` and `ℓ'!` for zeros `≤ T`;
- a multiple `≍ (T log T)^{(ℓ−ℓ')/2}` of its weight for `ℓ > ℓ'`, by padding.

`r4_moments_exact.py`, model 3, has a single relation, of length 6: `Δ_2 = … = Δ_5 = 0` and `Δ_6 = 2·6!/(2!3!1!)·weight`.

Moreover, Proposition E speaks only for `L ≥ C_2T/(α²δ log L)`. For example, `α = 0.0891`, `δ = 0.1`, `T = 10³` need `L ≥ 1.1·10⁵` with `C_2 = 1` (`r4_P5_P4.py` (P4c)). So the saturating relations can all be far longer than any fixed moment order `K`, and then they change no moment of order `≤ K` at all.

§0 item 3 also differs from (b) by a factor `L`: "`≍ T log T` relations", "relative amount `≍ (T log T)^{1−ℓ/2}`".

The qualitative conclusion, that fixed-order moment hypotheses with any fixed relative precision cannot detect saturation, survives, and in a stronger form. The displayed relative excess and the threshold for `η` are not derived. Replace them by the per-length statement above.

### P17. The barrier conclusion does not follow (§0 item 4; after Corollary P6)

> §0 item 4: "Any proof must use information that is unstable under arbitrarily small perturbations of the zeros, such as the explicit formula as an exact identity"
>
> After P6: "Hence any property of spectra that distinguishes the two cases is unstable under perturbations … That excludes everything determined by the asymptotic statistics of the ordinates: counting functions and local correlations at any fixed scale relative to the mean spacing."

Corollary P6's content is correct, conditionally on the draft: the listed data do not decide. The inference to proofs is not.

**Why.** A proof uses some property `Q` of the true zeros with "`Q` + data ⇒ no atom".
- The cyclic model `C` of Theorem 9.22 only has to fail `Q`.
- If `Q` is stable, the perturbations `C'` of P5 fail `Q` too.
- Nothing is contradicted.

What P5 and Theorem 9.22 show is narrower: no property shared by `C` and stable under small perturbations can imply "no atom", or characterise atoms on a class that contains `C`. To force `Q` to be unstable, one would need atom-carrying spectra arbitrarily close to every spectrum with the true zeros' properties. That is the converse density statement, and it is not shown.

A stable property that `C` may lack: the zero-counting function `N(t) = (t/2π)log(qt/2πe) + O(log t)`. It is preserved, up to `O(log t)`, by the perturbations of P5, and Theorem 9.22 as quoted does not assert it.

The Remark after P6 ("What remains outside every such model is … the explicit formula, exactly") overstates in the same way.

Fix: "any proof must use information that the models of Theorem 9.22 lack. It may be stable (asymptotic statistics of the true zeros that these lattice spectra violate) or unstable (the explicit formula as an exact identity)."

### L2. LINE1 Lemma 4.2

> "Let `ĥ ≥ 0` be smooth and even, supported in `[−U,U]`, equal to `1` on `[−U/2,U/2]` … `|Σ_γ m_γ ĥ(γ) e^{iγy}| ≤ η N_ĥ`"
>
> Sketch: "The Weil explicit formula makes the left side equal to `D(y) − P(y)`. … `P(y)` is the prime sum over a window of relative width `≍ 1/U` around `e^y`. Brun–Titchmarsh bounds it …"

The label "sketch" is honest, and the statement is plausible. Three problems:
1. **One-sided versus two-sided.** With `Γ` the positive ordinates (notation of `DENSITY.md`), the left side is complex. The explicit formula gives only `Σ_{γ∈±Γ} = 2Σ_{γ>0}m_γĥ(γ)cos(γy)` (`ĥ` even, zeros symmetric), not the sine part.
   - Numerically, at `U = 600`, the sine part is as small as the cosine part: 0.143 versus 0.158. But the sketch proves at most the two-sided version.
   - The two-sided version is what "The gap" uses: near-recurrence gives a cosine sum `≈ N_ĥ`.
2. **Uniformity.** `η_1, η_2` and the `o(1)` depend on the shape of `ĥ`. Take `ĥ(r) = Φ(r/U)` with `Φ` fixed.
3. **The prime side is not a window sum.** `ĥ` has compact support. So `g = (2π)^{−1}∫ĥ(r)e^{−iru}dr = Uφ(Uu)` is Schwartz, and also entire of exponential type with bounded derivative on `ℝ`. Such a function cannot decay exponentially, not even as `u → +∞` alone: the log-integral theorem forces `∫log|g(u)|/(1+u²)du > −∞` unless `g ≡ 0`. Hence `∫_0^∞ e^{u/2}|g(u)|du = ∞`, and by the prime number theorem in intervals of fixed relative length, `Σ_nΛ(n)n^{−1/2}|g(log n − y)| = ∞`. So `P(y)` converges only conditionally, under GRH. This is the very issue `PRIMECORR.md` §3 raises for frequency-compact windows.

   **Missing argument.** Split at `|log n − y| = 1`.
   - Near part: apply Brun–Titchmarsh on unit `s`-intervals against `sup|φ|`. This gives `≪ e^{y/2}·log x/log(x/U) + U y e^{−y/2}`.
   - Far part: partial summation with `ψ(t,χ) ≪ √t log²t`. This gives `≪ U y²∫_{|s|>U}(|φ|+|φ'|)(1+|s|/U)²ds = O_A(U^{−A})`.

   The sketch's `≪ e^{y/2}` misses the factor `log x/log(x/U)` for `x ∈ [U, U^{1+δ}]`. That is harmless.

**Numerics** (`r4_line1.py`, `U = 600`). The maximum ratio is 0.158 on `[K/U, 2log U + 2loglog U − C]`. It is attained at prime spikes `y = log n`, and these match the explicit-formula prediction `|Λ(n)χ(n)n^{−1/2}Uφ(0)|/(2N)` to three digits (e.g. 0.158 at `n = 7`). They are the `o_U(1)` term, of size `O(1/log U)`.

### L3. LINE1 "The gap" and the closing paragraph

Two counting claims are wrong as stated:
- > "counting `#Γ_U ≫ U log U` distinct zeros"

  Riemann–von Mangoldt gives only `#Γ_U ≫ U`, from `m_γ ≪ log γ`. `≫ U log U` needs a positive proportion of distinct zeros, i.e. a pair-correlation-type multiplicity bound under GRH. That bound is not cited, and I did not check the literature for a fixed `L(s,χ)`. The conclusion (Dirichlet return time `U^{1−o(1)}`) does not change with `#Γ_U ≫ U`.
- > "In rank exactly 1, the return time is `2π/s ≥ (1+o(1)) log U` by counting."

  In rank one, `Γ_U ⊂ sℤ`, so `2π/s ≥ 2π#Γ_U/U`. The bound `(1+o(1))log U` needs `#Γ_U ≥ (1−o(1))N(U)`, i.e. almost all zeros simple, which is not known. `LINE1.md` §1 itself says the logarithmic bound comes from the explicit formula, "without any information on multiplicities". Lemma 4.2 covers return times up to `(2+o(1))log U` anyway, so only the phrase "(B) excludes only `s ∈ (π/log U, 2π/log U]`" depends on this.

Minor points:
- "a basis without exceptional rational dependence": a basis has no rational dependence. What is meant is a Diophantine condition, and the return-time bound is a heuristic; label it.
- > "Corollary 1.4 is the only structural consequence of an atom known here"

  Proposition 4.1 in the same section, and `INVLO.md` Corollary C′ and Proposition E, are structural consequences of an atom.
- > "(a) An unconditional theorem excluding rank-one `G` whose denominators grow like `log T`. That is Problem 1 of arXiv:2407.03530"

  By the description in `DENSITY.md` §3.3, that problem asks for two ordinates with irrational ratio, i.e. to exclude rank one altogether. (a) is a special case of it, not the problem itself. I have not read the cited paper.
- "(c) … resonant spectrum. Combined with Remark 7.6 …": this needs "resonant" to mean "cyclic", because Remark 7.6 excludes only cyclic `G`.

## 3. Required fixes

1. **INVLO l. 36:** `r_γ = 4M/(M+iγ)` as in `DENSITY.md` §1; double the constant in `b` in Theorem C and Prop. E. Or keep `2M` and drop the identification of the law with `ν_M`. Check the same in PRIMECORR P2. (I1)
2. **INVLO §4:** replace "Cyclic models (`Γ' ⊂ βℕ`) are saturated as soon as `L ≥ 3`" by the correct statement (I18), with `{β, 3β, 5β}` as a warning example.
3. **INVLO §4 split:** state the window `C(1+Σ|c_γ|²)/(ηa)² ≤ L ≲ 64/ε` for the alternative `B ≳ ηaL`, and the uniform bound `Σ|c_γ|² ≤ C(M,A_0)` on `𝒞`. (I17)
4. **INVLO header l. 26–28, l. 281–284 and §4:** restore `δ`: "all but a proportion `O(δ)` … up to `C_2^{−1}α²δ L log L`". (I13)
5. **INVLO §3 item 4:** add `ψ(t)log(t+3)/t → ∞` for the continuity claim, and `ψ > 3` eventually for incomparability. (I16)
6. **PRIMECORR §0 item 1:** fix the sign (`(−1)^k M_k = μ_k + Δ_k`). (P1)
7. **PRIMECORR P1(b):** add the padding step to the proof. (P7)
8. **PRIMECORR P3(i):** record `C_1 = C_1(M,χ)` and add the monotonicity-in-`L` argument for all `M`. (P9)
9. **PRIMECORR P4(b) and §0 item 3:** replace the count and the threshold by the per-length statement (zero below the relation's length; multinomial times weight at it; padding above). Note that Prop. E's lengths satisfy `L ≫ T/log T`, and make §0 agree with §4. (P14)
10. **PRIMECORR §0 item 4, the paragraph after P6, and the Remark after P6:** replace "any proof must use information that is unstable …" by the correct, narrower statement. (P17)
11. **LINE1 Lemma 4.2:**
    - state it for the two-sided sum (or the real part);
    - take `ĥ(r) = Φ(r/U)` with `Φ` fixed;
    - in the sketch, say that the prime side converges only conditionally, and treat the far tail by GRH partial summation. (L2)
12. **LINE1 "The gap":** weaken or cite the two counting claims (`#Γ_U ≫ U log U`; "`(1+o(1))log U` by counting"). Delete "Corollary 1.4 is the only structural consequence of an atom known here". (L3)

## 4. Optional improvements

- **INVLO Remark (1):** "which is sharp by Lemma A" should read "by the Gauss–Hermite example after Lemma A".
- **INVLO, GLI:** "LI is GLI(∞)" should read "LI implies GLI(ψ) for every ψ" (ψ takes values in `[0,∞)`). Also fix the grammar of the sentence after it.
- **INVLO Theorem C:**
  - say explicitly that for `t'_χ ≤ T_L ≤ T_0` (moderate `T_L`) the Lemma D(b) bound is `O_{M,χ}(1)`;
  - say that `t_χ` is chosen so that `#I ≥ 4`.
- **INVLO Lemma D(b):** `#I ≥ 3` suffices.
- **INVLO Remark (2):** quantify "far too large". With Gaussian `W` at `x = 0`, the bound of Lemma A drops below 1 only for `N > 2·10⁷` (`r4_lemmaA.py` (A5)).
- **INVLO §4:** "a new sufficient hypothesis" sits oddly with "No novelty is claimed" in the header. Say "new here".
- **INVLO §3 item 5:** FD is a finite-height condition, not a condition on "the tail".
- **PRIMECORR, "What is new":**
  - "equivalent, up to the tolerance `c_0^ℓ` …": the tolerance is effectively much smaller, because of padding. A double zero `γ_0 > T` gives `Δ_ℓ ≥ ŵ(γ_0/T)²(2s_2)^{(ℓ−2)/2}` for even `ℓ`. Already at `ℓ = 2`, a double zero at any `γ_0 ∈ (T, πT]` violates H_corr(`2,T`), since `4ŵ(γ_0/T)² ≥ c_0²` there (`r4_P5_P4.py` (P4d), both bumps). So H_corr also constrains zeros well above `T`.
  - "What it adds is the explicit bound atom `≤ C/√L`": Prop. B already has this bound. What P2 adds is the finite-height, prime-side hypothesis.
- **PRIMECORR §0 item 3:** "the error terms … exceed the main term" should read "the error *bounds*".
- **PRIMECORR §0 item 4:** "perturbing each ordinate … makes the spectrum `ℚ`-independent" should read "*some* perturbation …".
- **PRIMECORR P6:** "small `M` … on the frequencies `Ωℕ`" compresses Corollary 9.23(iii), which works at weight ½ with `Ω = 1/(2M)`. State it that way. Mention that the counting constant `A_c` is astronomically large (draft §8 item 9), and that the perturbed spectra need `2A_c`.
- **PRIMECORR P4(a):** the moment formula ignores the multiplicity of the position of the minimum among the `n_i`, and has `Σ_{n≤X}` where `n ∈ [X,2X]` is meant. Harmless for orders of magnitude.
- **LINE1 Prop. 4.1, prose:** "`T ≍ ε^{−2}/α`" drops the factor `log²(1/(αε))` of the statement.

## 5. Scripts and key outputs (`monotonicity/density2/referee4/`)

| Script | What it checks | Key outputs (`*_results.txt`) |
|---|---|---|
| `r4_lemmaA.py` | Needle polynomial in exact rationals, `N ≤ 60`; the constants; moment bound; both halves of `E p(W)`; Christoffel function | All properties hold for all `N` tested; `(N+1)∫q² ≈ 4.19 ≤ 8π`; `min q ≈ 0.708 ≥ 0.27`; `max\|q\|/(2\|T_N\|) < 1` on `1 ≤ \|s\| ≤ 10³`. `112π = 351.86 < 352`; `8/(3π²) = 0.27019`; `0.0729 > 1/14`. The maximal atom of a law with the first `2N` moments of `N(0,1)` is `1.2532/√N` at `N = 10⁴` (→ `√(π/2)`). The bound of Lemma A is `< 1` only for `N > 2.0·10⁷`. |
| `r4_lemmaD.py` | Nicholson monotonicity; `J_0 ≤ e^{−y²/4}`; `\|J_0\| < π^{−1/2}` on `y > 2`; the `12/(πa)` integral; sup densities against both bounds | `sup\|J_0\|/√(2/πy) = 0.999999997` on `(0,5000]`; `y(J_0²+Y_0²) − 2/π < 0` up to `10⁶`; `max\|J_0\|` on `[2,2000]` is `0.4028`. Every density is below its bound, e.g. `(1,1,1)`: `0.289 ≤ 0.608`. |
| `r4_moments_exact.py` | Prop. B′ identity and Thm P1 in exact arithmetic; saturation of cyclic models | B′: agreement for `j ≤ L` in all cases; first failure at the length of the shortest relation outside `S_L` (3, then 6). `{1,3,5}`: `S_3 = ∅`. Model in `𝒞` `10⁻⁴·{60209, 102438, 129881}`: `σ² = 0.0856`, shortest relation of length 127. P1: identity, (a), (e) exact. Model 2: `Δ_5 = 161/90 > 0` from padding only. Model 3: `Δ_2..Δ_5 = 0`, `Δ_6 = 15/8 = 2·6!/(2!3!1!)·weight`. |
| `r4_explicit_formula.py` | (1.1) with a different bump and the closed-form Gamma kernel of (I5); P1(e) at `y = 0`; zero counting | (1.1): `\|diff\| ≤ 1.3·10⁻⁸` (`T = 3`), `≤ 8.9·10⁻⁷` (`T = 5`). The two forms of `Ω_T` agree (e.g. `−1.250·10⁻⁴` at `y = 6`). P1(e): `\|diff\| ≤ 1.6·10⁻¹²` for 5 cases (prime sums over `n ≤ 13`). Distinct ordinates in `(t,2t]` ≈ `0.158·t log t`. |
| `r4_examples.py` | INVLO §3 item 2 with `ψ(t) = 2t`, by exhaustive relation search; TS_ℤ gcds; item 3 | Shortest relations with largest index `n* = 2,3,4,5`: 61, 35, 40, `>26`, against required `2^{e_{n*}−e_{n*−1}} = 8, 8, 16, 16 ≥ ψ(x_{n*})`. The divisibility claim holds in all relations found; all gcds are 1. |
| `r4_P5_P4.py` | P5 counting law and `σ²` factor; the factor 2 attained; P4 numbers; H_corr versus a double zero above `T` | `σ²` ratios within `[(1+δ)^{−2},(1−δ)^{−2}]` for `δ = 0.1, 0.25, 0.49`. Count `1 → 2` in a unit interval. `s_2(T)/(T log T) ≈ 0.44–0.46` (`T = 5..60`). Prop. E's `L`-range: `L ≥ 1.1·10⁵` for `T = 10³` (`α = 0.0891`, `δ = 0.1`). `c_0 = 0.961`; a double zero in `(T, πT]` violates H_corr(`2,T`). |
| `r4_line1.py` | Arcsine concentration (Prop 4.1); Lemma 4.2 on χ₋₄ zeros, `U = 600` | Concentration `= 1/3` exactly at `a = 8ε`. Lemma 4.2 ratio `≤ 0.158` on `[10/U, 12.5]`; the sine part `≤ 0.143`. Spikes at `y = log n` match the explicit formula to `10⁻³`. |

## 6. What I could not check

- **Quoted standard inputs, checked only numerically or by recall:**
  - the explicit formula as `PRIMECORR.md` states it from memory (the cited textbook theorem); (1.1) and P1(e) agree to `10⁻⁸`–`10⁻¹²`, which confirms sign and normalisation;
  - the monotonicity of `y(J_0²+Y_0²)` from Nicholson's formula (numerics only);
  - the Kolmogorov–Rogozin inequality used in Prop. 4.1;
  - the hypotheses and ranges of the short-interval moment theorem cited in P4(a).
- **Lemma 4.2 in full.** I checked its plausibility, the structure of a complete proof (§2, L2) and numerics at `U = 600`. I did not write out the full proof of the endpoint `2log U + 2loglog U − C`.
- **The draft.** I did not referee `THM922_DRAFT.md`; I only checked that P6 restates Theorem 9.22 and Corollary 9.23 correctly. It does, with the compression noted in §4.
- **The literature.** I did not check `PRIMECORR.md` §6 or the novelty statements. I did not read arXiv:2407.03530, so its problem list is taken from the repository's own description.
- **Known results on simple zeros of a fixed `L(s,χ)` under GRH** (pair-correlation-type multiplicity bounds), mentioned in L3: recalled, not checked.
- **Numerics cannot test the asymptotic statements.** Lemma A's bound is trivial below `N ≈ 2·10⁷`, and the `o_U(1)` of Lemma 4.2 decays like `1/log U`. Those statements were checked by hand only.
