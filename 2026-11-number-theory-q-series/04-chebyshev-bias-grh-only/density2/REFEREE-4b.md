AI-produced (Anthropic's Claude), adversarial referee report, 2026-11-27.

# Referee report 4b: `PRIMECORR.md` §3 (Theorem P7), and the fixes to report 4

## 0. Scope, method, summary

**Refereed:**
- the new §3 of `PRIMECORR.md`, "Beyond LI on the prime side: band decoupling (Theorem P7)", in full: Lemma P0′, the Setting, H_dec, Lemma P1′, Theorem P7, the comparison paragraph and the closing remark;
- the §0 summary of Theorem P7;
- `p7_checks.py` and its results;
- the fixes to required items 1–12 of `REFEREE-4.md`, via `git diff 0f309d0 c97d788`, in `INVLO.md`, `PRIMECORR.md` (§0, §2, §4, §5) and `LINE1.md` §4.

As before, I did not open `PROGRESS.md`, `RESULTS.md`, the other referee reports or `referee1/`–`referee3/`. New scripts: `referee4/r4b_*.py`.

**Summary.**
- *Lemma P0′ and Lemma P1′ are correct.* P0′ needs two small precisions (NCZ; the admissible class of test functions).
  - The regularised identity used in P0′'s proof holds numerically to `10⁻⁷`–`10⁻⁹`.
  - The decoupled prediction `(−1)^aM_{a,0}μ_b(n_I)` is the right object, with the right signs.
  - (Z) and the window shapes do give weight exactly 1 at heights `≤ T*`.
- *Theorem P7 is wrong as stated.* Its hypotheses hold vacuously for a band containing no zero.
  - Example: for `χ_{−4}` take `(s_1, s_2, g) = (1, 2.5, 0.4)`. The first zero is at `6.02`, so (Z) holds and (D) holds for every `L`.
  - The theorem would then give "no atoms" unconditionally under GRH, i.e. the open problem.
  - The proof breaks at step 2: with fewer than three free ordinates there is no density bound.
  - With the hypothesis "each band contains a zero" (or `s_{1,j} → ∞`), the proof is correct, including its use of Prop. B′, Lemma D and the truncation.
- *"GLI implies H_dec" and "Theorem C is the special case" are wrong.* Condition (D) at `(a,b) = (0,2)` forces the band zeros to be simple (a double zero gives `E_{0,2} = 4`). GLI and Theorem C allow multiple zeros.
  - Consequently "strictly weaker than GLI" holds only for simple spectra, and §0's "contains … Theorem C as a special case" is wrong.
  - The incomparability with the cited Theorem 5.4(1) is asserted without examples. I give two below.
- *Fixes 1–12 of report 4:* all applied and correct. Three of them introduced minor new inaccuracies:
  - INVLO §4 says the cyclic case of Question P″ consists of the full-segment models; in fact it allows missing frequencies;
  - LINE1 Lemma 4.2's `P(y)` drops one prime term (negligible);
  - LINE1 "would need a positive proportion of simple zeros" should read "would follow from".

  Separately, PRIMECORR P4(b) now claims to "follow from" P1 inside a remark labelled heuristic (a labelling point).

## 1. Verdict table

### Part A: the new §3 of `PRIMECORR.md`

| # | Item | Verdict |
|---|---|---|
| A1 | Lemma P0′: existence of `Π_β(y)` (partial summation under GRH, boundary term) | CORRECT |
| A2 | Lemma P0′: the regularisation (admissibility of `h_ε`; dominations uniform in `ε ≤ 1` on zero, `Ω` and prime sides; `Ω_β → 0`; `Π_β^− → 0`) | CORRECT WITH GAP |
| A3 | Lemma P0′ statement: label "(GRH)", with `G_β` summed over positive ordinates | WRONG AS STATED (minor: NCZ needed) |
| A4 | Setting (windows `β_I, β_R`) | CORRECT |
| A5 | H_dec (Z), and "On the prime side" (`η` "positive exactly on the margins") | CORRECT WITH GAP (closed versus open margins) |
| A6 | H_dec (D): the decoupled prediction `(−1)^aM_{a,0}μ_b(n_I)` (signs, multiplicities, threshold 1) | CORRECT (it forces simple band zeros) |
| A7 | Lemma P1′ (positivity; the sequences counted; weight exactly 1 at heights `≤ T*`) | CORRECT |
| A8 | Theorem P7 as stated (no condition that the band contain zeros) | WRONG AS STATED |
| A9 | Theorem P7, steps 1–4 for nonempty bands with `s_1` large (Prop. B′, Lemma D(b), truncation with `T* = N⁴`, limit) | CORRECT |
| A10 | Comparison: "LI implies H_dec" | CORRECT (needs simple zeros, as the text's own body says) |
| A11 | Comparison: "GLI implies H_dec … Theorem C of `INVLO.md` is the special case"; §0 item 2 "contains … (Theorem C) as a special case" | WRONG AS STATED |
| A12 | Comparison: "H_dec is strictly weaker than GLI", with the example | WRONG AS STATED (true for simple spectra; the example is under-specified) |
| A13 | Comparison: incomparable with the cited Theorem 5.4(1) | CORRECT WITH GAP (no examples given) |
| A14 | Remark (why plain moments do not suffice) | CORRECT |
| A15 | `p7_checks.py` and its description in the §3 intro | CORRECT (minor wording) |
| A16 | §0 item 2, second bullet (summary of P7) | WRONG AS STATED (inherits A8, A11) |

### Part B: the fixes to report 4

| # | Required fix of `REFEREE-4.md` | Verdict |
|---|---|---|
| B1 | `r_γ = 4M/(M+iγ)` in INVLO; constants in Thm C, Prop. E; PRIMECORR P2 | DONE, correct (no stale `2M` left) |
| B2 | INVLO §4, cyclic saturation | DONE; new minor inaccuracy (see §2) |
| B3 | INVLO §4, the window for `B ≳ ηaL` and the uniform bound on `Σ\|c_γ\|²` | DONE, correct |
| B4 | INVLO Prop. E summaries (header, after the statement, §4) | DONE, correct |
| B5 | INVLO §3 item 4 (hypotheses on ψ) | DONE, correct |
| B6 | PRIMECORR §0 item 1 (sign) | DONE, correct |
| B7 | PRIMECORR P1(b) (padding) | DONE, correct |
| B8 | PRIMECORR P3(i) (`C_1(M)`) | DONE, correct |
| B9 | PRIMECORR P4(b) and §0 item 3 | DONE, correct (a labelling point, §4) |
| B10 | PRIMECORR barrier (§0 item 4, after P6, the Remark) | DONE, correct |
| B11 | LINE1 Lemma 4.2 | DONE, correct (one prime term omitted, negligible) |
| B12 | LINE1 "The gap" and the closing paragraph | DONE, correct (one wording slip) |

The optional items of report 4 were also applied correctly: Remark (1), the GLI wording, the moderate range in Thm C, `#I ≥ 3`, Remark (2), "new here", §3 item 5, the "What is new" paragraph, "error bounds", "some perturbation", the P6 restatement, the P4(a) caveats, and the LINE1 prose.

## 2. Findings

### A8. Theorem P7 holds vacuously for an empty band

> "Suppose that along a sequence `L_j → ∞` the hypothesis H_dec(`L_j; s_{1,j}, s_{2,j}, g_j, T*_j`) holds, with: `T*_j := ⌊L_j/2⌋⁴`; `2s_{1,j} ≤ s_{2,j} ≤ 3s_{1,j}`; `s_{1,j}/log(s_{1,j}+3) = o(L_j)`. Then `ν_M` has no atoms"

Nothing requires the band to contain a zero.

**Counterexample to the statement.** Take `χ_{−4}`, `(s_1, s_2, g) = (1, 2.5, 0.4)` for every `j`.
- (Z): there is no zero in `[0.6, 1.4] ∪ [2.1, 2.9]`; the first ordinate is `6.02`.
- The band is empty, so `G_I = 0` and `n_I = 0`. By Lemma P0′, `Π_I(y) → 0`, so `M_{a,b} = 0` for `b ≥ 1`. Since `μ_b(0) = 0` for `b ≥ 1`, (D) reads `0 < 1` for every `L`.
- `2·1 ≤ 2.5 ≤ 3·1`, and `s_1/log(s_1+3) = 1/log 4 = o(L_j)`.
- NCZ is known for `q = 4`.

So under GRH alone every hypothesis of Theorem P7 holds, and the theorem would prove `ν_{1/2}({0}) = 0`, which this repository treats as open. `r4b_P0prime.py` (R4) confirms that (Z) holds, that `Π_I(y) ≈ Ω_β(y)` (within `4·10⁻⁴`), and that the one-level count gives `n_I ≈ 10⁻⁴`, i.e. 0. `r4b_P1prime.py` (Q5) has `E_{a,b} ≡ 0` for an empty band.

**Where the proof breaks.** Step 2 says: "For bounded `s_1`, use Lemma D(a) instead." Lemma D(a) needs three free ordinates with coefficient `≥ a`. With `n_I ≤ 2` the free part from the band has no bounded density (one cosine has the arcsine law, two have a logarithmic singularity).

**Fix.** Add "each band contains at least one zero" (`n_I ≥ 1`), or require `s_{1,j} → ∞`.
- If `s_{1,j} → ∞` along a subsequence, the band contains `≍ s_1 log s_1` zeros, and the proof as written works.
- If `s_{1,j}` stays bounded along a subsequence, then some zero `γ_0` lies in infinitely many bands. By Lemma P1′ it is then in no relation at all: any relation has a finite length and a finite height. So `γ_0 ∉ ⟨Γ∖{γ_0}⟩_ℚ`, TS holds, and `DENSITY.md` Theorem 3.1 gives absolute continuity.

The same correction is needed in §0 item 2 ("If, for infinitely many `L`, some band … is free in this sense, then `ν_M` has no atoms").

### A11, A12, A16. (D) forbids multiple zeros in the band; GLI does not

> "GLI implies H_dec. … So Theorem C of `INVLO.md` is the special case of Theorem P7 obtained from GLI, via the duality."
>
> "H_dec is strictly weaker than GLI, as conditions on abstract spectra."
>
> §0: "it contains the graded LI of `INVLO.md` (Theorem C) as a special case."

The decoupled prediction gives each copy of a band zero an independent phase. That is what keeps it one-level, since only `n_I` enters. As a consequence, two copies `ρ ≠ ρ'` of a multiple band zero form the multiset relation `e_ρ − e_{ρ'}`, with `κ_I ≠ 0`. Hence
`E_{0,2} = Σ_{γ∈I} 2m_γ(m_γ − 1)`,
which is `≥ 4` as soon as one band zero is double. `r4b_P1prime.py` (Q4) gives `E_{0,2} = 4` exactly. So (D) at `(a,b) = (0,2)` forces every band zero to be simple.

GLI(ψ) is a condition on relations among *distinct* ordinates, and `INVLO.md` Theorem C states "Zeros are not assumed simple". Therefore:
- GLI does not imply H_dec. Take an abstract spectrum that is `ℚ`-independent (so GLI(ψ) holds for every ψ) with every multiplicity equal to 2. Then (D) fails for every nonempty band.
- Theorem C is not a special case of Theorem P7. Theorem C covers spectra with multiple zeros, for which H_dec fails.
- "Strictly weaker than GLI" is true only for simple spectra (abstract sets).

The strictness example, "an LI band of zeros at height `T_j` for each `j`, embedded in a spectrum with short relations everywhere else", must also say:
- that the band zeros are simple;
- that each band has `≥ 1` zero;
- that the band zeros are free of relations of length `≤ L_j` with **all** zeros up to `T*_j` (e.g. `ℚ`-independent of everything else);
- that `s_{1,j}/log s_{1,j} = o(L_j)`.

**Correct statements.** GLI(ψ), together with simplicity of the zeros in one band per scale (e.g. in `(s_1, s_2)` with `s_1` just above `T_{L}`, as in the text), implies H_dec. Theorem C restricted to such spectra is then a special case of Theorem P7. For simple spectra, H_dec is strictly weaker than GLI.

Similarly, "LI implies H_dec" (A10) needs "LI with simple zeros", as the bullet's own body says. §0's "strictly weaker than LI" should read "than LI with simple zeros".

### A3, A2. Lemma P0′: NCZ, and the admissible class of test functions

**NCZ (A3).** The lemma is labelled "(GRH)", and `G_β` sums over the positive ordinates `Γ`. If `L(½,χ) = 0` with multiplicity `n_0`, the zero side of the explicit formula contains `n_0β(0)`. So `Π_β(y) = −G_β(ι(y)) − n_0β(0) + o(1)`, and the statement fails for every `β` with `β(0) ≠ 0` (e.g. `β_R`). Lemma P1′ and Theorem P7 assume NCZ, so nothing downstream breaks.

Fix: label the lemma "GRH + NCZ", or add the term `−n_0β(0)`.

**Admissibility (A2).** "this test function is admissible in the standard sense": the explicit formula of §1 is stated only for `g ∈ C_c^∞`, and `g_ε(u) = b(u−y)e^{−ε²u²/2}` is not compactly supported. The version needed is standard: `h` holomorphic in `|Im r| ≤ ½ + δ` with `|h(r)| ≪ (1+|r|)^{−1−δ}`. `h_ε` satisfies it, since it is entire with Gaussian decay in every horizontal strip.

Alternatively, reduce to `C_c^∞` by a smooth cutoff at `|u| ≤ U`. The cutoff error on the prime side is `≪ ∫_U^∞ e^{u/2 − ε²u²/2}du → 0`.

`r4b_P0prime.py` (R1) confirms the regularised identity for `ε = 0.5, 0.7, 1.0` to `|diff| ≤ 1.4·10⁻⁷`.

**The dominations.** I checked each one:
- *Zero side.* `|h_ε(γ)| ≤ ‖β‖_∞` for `|γ| ≤ max(2, 2S)`, where `supp β ⊂ [−S,S]`. Otherwise `|h_ε(γ)| ≤ C e^{−γ²/8}` uniformly in `ε ≤ 1`, because `ε ↦ ε^{−1}e^{−γ²/(8ε²)}` increases on `(0,1]` when `|γ| ≥ 2`.
- *`Ω` side.* `|Ω_χ| ≪ log(|r|+2)`.
- *Prime side.* After partial summation the integrand is dominated by `C(1+u)²[|b'(u−y)| + (1+|u|)|b(u−y)|]`, an integrable majorant independent of `ε ≤ 1`.

The text states only that "the bound depends on `f` only through `∫u²(|f|+|f'|)`". That uniform bound is not enough for dominated convergence; the majorant above should be written. This is a gap that can be filled.

The boundary term, `Ω_β → 0`, and `Π_β^− ≪ ∫_0^∞(1+u)²(|b|+|b'|)(y+u)du → 0` are all correct.

### A13. Incomparability with the cited Theorem 5.4(1): examples

> "Neither is a special case of the other, as conditions on frequency sets."

As written, this contrasts the two conditions but proves nothing. Examples (abstract spectra):
- **5.4(1) does not imply H_dec.**
  - Easiest: double every element of a set satisfying 5.4(1). That condition concerns `ℚ`-spans of ordinates, so it survives, while (D) fails at `(0,2)` for every nonempty band (A11).
  - With simple elements: start from a set satisfying 5.4(1) in which every element lies in some relation, such as `DENSITY.md` Example 3.3a (primary elements linked to far higher ones). For each `k`, add new `ℚ`-independent elements `u_k, v_k, u_k+v_k`, with some `u_k+v_k` in every interval of length `½` above height 1. Then:
    - every admissible band with `s_1 > 2.5` has an interior `(s_1+g, s_2−g)` of length `> ½`. So it either violates (Z) or contains some `u_k+v_k`, which lies in a length-3 relation with lower elements, and `E ≥ 1`;
    - every band with bounded `s_1` fails (D) for large `L_j`, because each of its finitely many elements lies in a relation of finite length and height;
    - the original free elements stay free, so 5.4(1) still holds.
- **H_dec does not imply 5.4(1).** Let `Γ` consist of:
  - `R = {nα}`;
  - bands `B_j` of `≥ 1` simple elements in `(s_{1,j}+g, s_{2,j}−g)`;
  - for each `x ∈ B_j`, two elements `z_x ∈ (T*_j+1, T*_j+2)` and `x + z_x`;

  with all of `⋃B_j ∪ {z_x} ∪ {α}` `ℚ`-independent. Then:
  - Choose the scales sparse: `s_{1,j+1} > T_j := T*_j + 3 + s_{2,j}`, with the margins of each band avoiding all other elements (take `α = 1` and `g < ¼`). This keeps the linking elements of earlier bands below the next band.
  - Every relation through `x ∈ B_j` uses `z_x` and `x + z_x`, which lie above `T*_j + 1`, outside `supp β_R`. Relations among elements of other bands and of `R` have `κ_I = 0`. So `E_{a,b} = 0`, and H_dec holds along `L_j` (with `s_{1,j}/log s_{1,j} = o(L_j)`, e.g. `L_j ≥ j·s_{1,j}`).
  - At `T = T_j`, every element `≤ T^{1/2}` lies in the `ℚ`-span of the other elements of `Γ_T`:
    - `nα` via `2nα`;
    - `x ∈ B_i` (`i ≤ j`) via `(x+z_x) − z_x`;
    - `z_x` via `(x+z_x) − x`.

    So 5.4(1) fails.

  The same spectrum fails GLI(ψ) for every `ψ > 3`, and it satisfies TS. So it can also serve as the strictness example of A12.

### A5. (Z) as a one-level condition

> "`η ≥ 0` is a bump that is positive exactly on the margins"

The margins in (Z) are closed, `[s_1−g, s_1+g]`. A continuous `η` cannot be positive exactly on a closed interval, because its positivity set is open.

Lemma P1′ needs only that no zero lies in the *open* margins:
- at `γ = s_1 − g`, `β_R = 1` and `β_I = 0`;
- at `γ = s_1 + g`, `β_I = 1` and `β_R = 0`;
- and likewise at `s_2 ± g`.

Fix: state (Z) with open margins; then take `η` an even `C_c^∞` bump positive exactly on them.

### A6. The decoupled prediction is the right object

In P1′ the hull expectation `(−1)^{a+b}M_{a,b} = E[G_R^aG_I^b]` is compared with `E[G_R^a]·E[(G_I^{ind})^b]`. This keeps every relation among the R-zeros, so allowed relations cancel, and it gives each band copy its own phase.
- Signs: `M_{a,0} = (−1)^aE[G_R^a]`, so `(−1)^aM_{a,0} = E[G_R^a]`. This is correct.
- Multiplicities: `n_I` counts copies, which is what makes the prediction one-level. The consequence is A11.

`r4b_P1prime.py` checks the identity in exact arithmetic, together with the brute-force count of the sequences with `κ_I ≠ 0` and total frequency `0`. Two cases:
- (Q2) Relations among R-zeros contribute nothing to `E`. The plain moment excess of the combined window is `6, 0, 460` for orders 3–5, so a single window does see them.
- (Q3) A band zero in a length-3 relation with R-zeros gives `E_{2,1} = 4 ≥ 1`.

### A15. `p7_checks.py`

The script checks Lemma P0′ for a *plateau-free* bump, not a window of the Theorem P7 type. Its docstring includes `Π_β^−`, but the code does not compute it; at `y ≥ 6` it is negligible. The §3 intro says the sums "match … to `10^{−3}`–`10^{−4}`, improving with the truncation". In fact:
- at `V = 2` the differences reach `7·10⁻³`;
- at `y = 9` the agreement is not monotone in `V` (`6.7·10⁻⁴`, `1.07·10⁻³`, `3.6·10⁻⁴`).

My check (`r4b_P0prime.py` (R2)) uses a plateau window `(8.13, 19.87)` with `g = 1.4`, both prime terms and smooth cutoffs. It gives agreement to `≤ 7·10⁻³`, limited by slow conditional convergence. The one-level count (R3) gives `n_I = 4.0008` for 4 zeros. Correct, with a wording fix.

### B2. A new inaccuracy in INVLO §4

> "Cyclic models with **all** frequencies `β, 2β, …, Nβ` present … These are the models of the cyclic case of Question P″ (`DENSITY.md` Remark 9.11), the predecessor's line (2), and those of its draft Theorem 9.22."

Question P″ allows any `Γ ⊂ G ∩ (0,∞)`, and Remark 9.11 states its cyclic case with `0 ≤ m_n ≤ A log(n+2)`. So that case contains models with missing frequencies (e.g. `m_n = 1_{n odd}`), not only the full segments. The next bullet ("uniform when `β` is bounded below, as in the cyclic case of Question P″") treats it correctly. The two bullets are inconsistent.

Fix: "These include the models of the draft Theorem 9.22. The cyclic case of Question P″ allows missing frequencies and is covered by the next bullet."

### B9, B11, B12. Minor points in the fixes

- **PRIMECORR P4(b)** now opens with "This part follows from Theorem P1 and `INVLO.md` Proposition E", inside a remark labelled heuristic. The per-length contributions and the length range are rigorous consequences. The final sentence ("It would take moments of order `≍ T/log T` …") is not. Either split (b) into a labelled observation plus a heuristic sentence, or drop "follows from".
- **LINE1 Lemma 4.2 sketch:** "`P(y) = Σ_nχ(n)Λ(n)n^{−1/2}g(log n − y)`" omits the second prime term `g(−log n − y)`. It is `O_A(U^{−A})` by the far-part bound (and by Schwartz decay for `n = 2` when `y` is tiny), so this is harmless, but it should be written.
- **LINE1 "The gap":** "`≫ U log U` distinct zeros would need a positive proportion of simple zeros" should read "would follow from …". Simplicity is sufficient, not necessary: if all zeros were double, `#Γ_U ≈ N(U)/2`. A bound `Σ_{γ≤U}m_γ² ≪ N(U)` is what is needed.

## 3. Required fixes

1. **Theorem P7:** add "each band contains at least one zero" (or `s_{1,j} → ∞`). Treat bounded `s_1` via TS and `DENSITY.md` Theorem 3.1 (A8). Make the same change in §0 item 2 and in the description of the band in the §3 intro.
2. **Comparison, GLI bullet and §0 item 2:**
   - replace "GLI implies H_dec" by "GLI plus simple zeros in one admissible band per scale implies H_dec";
   - delete "Theorem C … is the special case of Theorem P7", or restrict it to simple spectra;
   - delete "contains … (Theorem C) as a special case" in §0. (A11, A16)
3. **"Strictly weaker than GLI":** state it for simple spectra. Make the example precise: simple band zeros, `≥ 1` per band, free of relations with all zeros up to `T*_j`, `s_{1,j}/log s_{1,j} = o(L_j)`. The spectrum in A13 works. (A12)
4. **Lemma P0′:** label it GRH + NCZ (or add `−n_0β(0)`). State the admissible test-function class of the explicit formula used (or the cutoff reduction to `C_c^∞`), and write the `ε`-uniform majorant on the prime side. (A3, A2)
5. **Theorem 5.4(1) comparison:** add the two examples (A13), or weaken to "not shown to be comparable".
6. **INVLO §4:** correct "These are the models of the cyclic case of Question P″" as in B2.

## 4. Optional improvements

- "LI implies H_dec" should read "LI with simple zeros implies H_dec"; likewise "strictly weaker than LI" in §0 item 2.
- **(Z):** use open margins (A5).
- "free of **short** relations up to `T* = L⁴`" should read "`T* = ⌊L/2⌋⁴`".
- **§3 intro, the description of `p7_checks.py`:** say that the window has no plateau, that `Π^−` is omitted (negligible for `y ≥ 6`), and give the actual range: `≤ 7·10⁻³` at `V = 2`, `≤ 4·10⁻⁴` at `V = 4`, not monotone at `y = 9`. Optionally cite the regularised-identity check (R1).
- **Theorem P7:** note that for bounded `s_1` the hypothesis already implies TS, so the theorem's content lies in bands with `s_1 → ∞`.
- **PRIMECORR P4(b):** label (B9).
- **LINE1:** add the second prime term in Lemma 4.2's sketch, and change "would need" to "would follow from" (B11, B12).

## 5. Scripts and key outputs (`monotonicity/density2/referee4/`)

| Script | What it checks | Key outputs (`*_results.txt`) |
|---|---|---|
| `r4b_P1prime.py` | Lemma P1′ in exact arithmetic on finite multiset models, against a brute-force count of the sequences with `κ_I ≠ 0` and total frequency `0` | Identity and `E ≥ 0` hold in all models. (Q2) allowed R-relations give `E = 0`, while the plain moment excess is `6, 0, 460` (orders 3–5). (Q3) a band zero in a length-3 relation gives `E_{2,1} = 4`. (Q4) a double band zero gives `E_{0,2} = 4`, so (D) forces simplicity. (Q5) an empty band gives `E ≡ 0`, so (D) is vacuous. (Q6) a relation through a zero in `(T*,T*+1]` with weight `1/3` gives `E_{2,1} = 4/3`. |
| `r4b_P0prime.py` | Lemma P0′ for `χ_{−4}` with P7-type plateau windows, both prime terms; the regularised identity; `n_I` at `y = 0`; the empty band | (R1) regularised identity, `ε = 0.5, 0.7, 1.0`: `\|diff\| ≤ 1.4·10⁻⁷`. (R2) limit identity, smooth cutoffs at `u ≈ 12–16`: `\|diff\| ≤ 7·10⁻³` (slow conditional convergence). (R3) `n_I = 4.0008` for the band `(8.13, 19.87)` with 4 zeros. (R4) empty band `(1, 2.5, 0.4)`: (Z) holds, `Π_I(y) ≈ Ω_β(y)` within `4·10⁻⁴`, `n_I ≈ 10⁻⁴`. |

## 6. What I could not check

- **The explicit formula:** I did not check the precise admissible class of test functions against a source. I checked it only numerically, for the Gaussian-regularised test functions.
- **Lemma P0′ in the limit `ε → 0`:** this cannot be reached numerically. The limit identity (R2) converges slowly, and its agreement (`10⁻³`) reflects conditional convergence, not an error.
- **The cited Theorem 5.4(1):** I took its hypothesis from the `DENSITY.md` §3.3 table and did not reread the source.
- **The remaining text:** I did not reread `PRIMECORR.md` §§1–2 and §§5–6 beyond the diffs, nor `INVLO.md` or `LINE1.md` §§1–3 beyond the changed lines.
