# Second independent referee report: proofs/KRTW-general.md (A) and proofs/T1-A-even.md (B)

Date: 2026-11-26. I formed this verdict before reading `REFEREE-KRTW-general.md`. Note that the header of (A)
already summarises that report, so I was not blind to its conclusions. All code is my own and written from the
statements alone, not from `code/`:
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/krtwref/`
(`lib.py`, `test_i.py`, `test_i_big.py`, `test_refine.py`, `test_hl.py`, `test_B.py`, with a `.log` file for each;
the downloaded sources are in `dl/`, `gow/` and `pa/`).

## Verdict

* **(A) PROOF COMPLETE.** The main line (Lemma 2 → Lemma 3 → Lemma 4 → Theorem 1 → R → ∞ → Proposition 6 →
  x = 1 → Jacobi triple product) has no gap. I re-derived every step by hand and found no error. Every
  numbered statement also passed an independent exact or truncated check.
* **Statement matches the printed conjecture: YES.** I compared it with the KRTW arXiv TeX source, v2, which
  is the version accepted by Selecta Math. (DOI 10.1007/s00029-026-01198-y).
* **Novelty: NO PRIOR PROOF FOUND** for identity (i) when k ≥ 3, nor for Conjecture 2.2 at z = 1 when k ≥ 3.
  **PRIOR PROOF FOUND for the second half.** The evaluation of the reduced k-fold sum as the product is known
  for every k:
  - Chern–Li–Stanton–Xue–Yee 2024;
  - Kanade–Lovejoy, arXiv:2410.15408 (SIGMA 21 (2025) 021), Theorem 1.1 with a = 0. They prove it with
    Bailey pairs;
  - Andrews 2010 and Kim–Yee 2013, as attributed in Warnaar arXiv:2511.17034, display (Eq_AKY) after
    Corollary 5, case i = 0;
  - Dousse–Jang arXiv:2603.05300, Theorem "thm:AK", a = 0.

  So the new mathematical content is the first equality of Theorem 2 (the 2k-fold KRTW sum equals the k-fold
  sum), together with its R-refinement (Theorem 1 and Lemma 4).

## 1. Statement check against KRTW

**Identity (i).** It is printed twice in the KRTW source:
- in the Introduction, as one of "a number of new conjectures of Rogers–Ramanujan type";
- again right after `\begin{conjecture}\label{Con_Shun2}`, which is Conjecture 4.9 by counter count. There it
  is the first of the "conjectural Andrews–Gordon-type identities":

    Σ_{r_1..r_k ≥ 0, s_1..s_k ≥ 0} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}})
        = (q,q^{k+1},q^{k+2};q^{k+2})_∞ / ((q;q²)_∞ (q;q)_∞).

How this matches the document:
- **Ranges.** KRTW set r_{k+1} = s_0 := 0. With their convention (a;q)_n = (a)_∞/(aq^n)_∞, 1/(q)_n vanishes for
  negative n. (Their text says "nonnegative", a typo.) The vanishing forces r_1 ≥ … ≥ r_k and 0 ≤ s_1 ≤ … ≤ s_k,
  exactly as the document has it.
- **Product side.** Identical.
- **Identities (ii)–(iv)** of the document also match the printed ones, including Ω.

**Conjecture 2.2.** This is `Con_Shun`, Section 2, counter 2.2, eq. (Eq_q2):

    Σ_{λ_1≤k} (zq)^{|λ|} P_{2λ}(1,q,q²,…;q²) = Σ_{r,s} ∏ z^{r_i+s_i} q^{(r_i+s_i)²+s_i²+s_i} / (…).

At z = 1 its right side is the left side of (i).

GOW arXiv:1401.7718, Theorem 1.2 (C_n^{(1)}), eq. (1.9) by my count of the TeX, with n = 1, m = k and κ = 2k+4,
gives

    Σ_{λ_1≤k} q^{|λ|} P_{2λ}(…;q²) = (q²;q²)_∞(q^{k+2};q^{k+2})_∞ θ(q;q^{k+2})/(q)_∞²,

which equals the product of (i). So for k ≥ 1, Conjecture 2.2 at z = 1 is equivalent to (i), and Corollary 4
follows. The theorem needs m ≥ 1; k = 0 is trivial.

I also computed the Hall–Littlewood side directly from Macdonald's tableau formula for P_λ. It agrees with the
multisum and with the product for k = 1, 2, 3, 4 mod q^16. As a control, t = q, k = 1 reproduces the
Rogers–Ramanujan series G(q).

**"Conjecture 4.9, first line, at z = 1/q" (§6).** I count (Eq_CDn2) as KRTW (4.7). With the C_1^{(1)} product
theorem that KRTW attribute to Meurman–Primc, this remark is consistent. I did not read Meurman–Primc.

## 2. Line-by-line check of the proof (all confirmed)

- **(0.1) and the exponent in Lemma 4.** (r+2b)(r+2b+1)/2 + r(r−1)/2 = r² + 2rb + 2b² + b.
- **Lemma 2.** It reads off the coefficient of z^L in (−z;q)_R/(z²;q²)_∞ = 1/((z;q)_∞(−zq^R;q)_∞). Correct.
- **Lemma 3.** Correct. The exponent bookkeeping is t² + i(i−1)/2 + i(d+1−t), and (q^{d+1−t})_t vanishes for t > d.
- **Lemma 4 (∗).** Both sides reduce to (2.3). The left side uses Lemma 3 with y = xq^{2b} and d = M−2b ≥ 0, and
  [M−2b;t](q)_t/(q)_{M−2b} = 1/(q)_{M−2b−t}. (∗∗) and (∗∞) are correct.
- **Theorem 1, base case.** A_1(R;x) is exactly the left side of (∗∞) and B_1 its right side.
- **Theorem 1, inductive step.** The shifts s_i = b + s'_i and N_i = N'_i + 2b give
  N_i(N_i+1)/2 = N'_i(N'_i+1)/2 + 2bN'_i + 2b² + b. With
  (q)_R/((q)_{R−r_1}(q)_{r_1−r_2}) = [R;r_1](q)_{r_1}/(q)_{r_1−r_2}, this yields (3.1) with A_{j−1}(r_1; xq^{2b}).
  The shift M' = M + 2b reproduces the weights and moves the bound to M'_{j−1} ≥ 2b. Applying (∗∗) at
  M = M'_{j−1} gives B_j exactly. The interchange of summation is legitimate in Z[x][[q]]: both the q-order and
  the x-degree tend to infinity.
- **R → ∞.** (q)_R/(q)_{R−r_1} ≡ 1 mod q^{R−r_1+1}, a term has q-order ≥ r_1², and r_1² + R − r_1 + 1 ≥ R + 1. So
  A_k(R) ≡ KRTW_k mod q^{R+1}, and likewise B_k(R) ≡ Mel_k mod q^{R+1}. Correct.
- **Proposition 6, the limiting Bailey lemma.** I re-derived the limit ρ_1 → ∞, ρ_2 = −√(aq) from the general
  Bailey lemma: α' = a^{r/2}q^{r²/2}α and β' as stated.
- **Proposition 6, the specialisation a = x²q.** √(aq) = xq and a^{r/2}q^{r²/2} = x^r q^{r(r+1)/2}.
- **Proposition 6, the iteration.** β^{(1)} = 1/((q)_n(−xq)_n). Each later step cancels the (−xq)_k, and there
  are k+1 applications. The n → ∞ limit (orders ≥ M_1 and ≥ r) is correct, and
  (x²q²)_∞(1−x²q) = (x²q)_∞ gives (5.1).
- **Theorem 3.** At x = 1, (x²q)_n/(q)_n = 1. The exponent is ((k+2)n² + kn)/2. The map n ↦ −n−1 accounts for
  the −q^{2n+1} terms. JTP with Q = q^{k+2} and z = q^{k+1} gives (q,q^{k+1},q^{k+2};q^{k+2})_∞, and
  (−q)_∞ = 1/(q;q²)_∞.

## 3. Independent recomputation

| object | range | result |
|---|---|---|
| (i): the 2k-fold sum by direct enumeration of (r,s), vs the product | k = 1, 2 mod q^400; k = 3, 4 mod q^300; k = 5 mod q^250; k = 6 mod q^200; k = 7 mod q^160; k = 8 mod q^140 | all equal |
| Mel_k(1) (k-fold sum) vs the product | same ranges | all equal |
| negative control: k = 3 sum vs k = 4 product | mod q^60 | detected unequal |
| Lemma 3 | R, d ≤ 7, exact polynomial identity in (y, q) | true |
| Lemma 4 (∗), after multiplying by (q)_M²(q²;q²)_M | R, M ≤ 6, exact in (x, q) | true |
| Theorem 1, A_j(R;x) = B_j(R;x), A summed from its definition | (j, R_max) = (1,7), (2,6), (3,4), (4,3); mod q^40, exact in x | all true; control A_2(2) ≠ B_2(3) |
| Theorem 2, both equalities (KRTW_k(x) = Mel_k(x) = right side of (5.1)) | k ≤ 5, mod q^45, exact in x | true |
| Conjecture 2.2 at z = 1, Hall–Littlewood side | k ≤ 4, mod q^16 | equals the multisum and the product |

"Exact in x" means x was specialised to 2^256. A polynomial in x whose integer coefficients have absolute value
below 2^255 is determined by its value there, so these are exact identities in Z[x][[q]] mod q^N.

The finite refinement is not a polynomial identity: s is unbounded in A_j, and B_j carries (xq^{…})_∞. Its
polynomial core is (∗), which I checked exactly.

## 4. Findings, by severity

1. **[Moderate; literature, in the body of (A) §5]** The sentence "For odd k the statement is not in the texts
   I read … I claim no novelty for §5" is factually wrong as a statement about the literature.
   Mel_k(1) = product is known for all k:
   - CLSXY 2024 (J. Algebraic Combin.);
   - Kanade–Lovejoy arXiv:2410.15408 = SIGMA 21 (2025) 021, Theorem 1.1, a = 0. Here
     ∏[n_{i+1};n_i]/(q)_{n_m} = 1/∏(q)_{M_i−M_{i+1}}, a Bailey-pair proof valid for all m ≥ 1;
   - Warnaar arXiv:2511.17034 (SIGMA 22 (2026) 062), (Eq_AKY) case i = 0, citing Andrews 2010 Thm 3,
     CLSXY Thm 3.5, Kim–Yee 2013 Thm 1.4 and Kanade–Russell;
   - Dousse–Jang arXiv:2603.05300, Theorem thm:AK.

   The header of (A) corrects this for Warnaar only. Kanade–Lovejoy and Dousse–Jang should be added. The x-refined
   form (5.1) itself I did not find stated there; Kanade–Lovejoy use other pairs and weights. This is not
   exhaustive.
2. **[Minor; header of (A)]** The header of (A) already flags three statements. I agree with all three:
   - "known before: k = 1 (all four)" is supported in KRTW only for (i), (ii) and (iv);
   - the justification of J_k = C_k/C_1 in §8 by "x → x/q" should be replaced by rerunning Proposition 6 in
     base q²;
   - Remark (a) "both forms" is unproved for the second form.

   None of them touches the main line.
3. **[Minor; §5]** The identification with Bressoud–Ismail–Stanton Theorem 5.2 for even k was not checked by me.
   It is now irrelevant: the document's own Bailey proof is complete, and the identity has other sources (item 1).
4. **[Cosmetic]** KRTW write "1/(q;q)_n vanishes for n a nonnegative integer" (meaning negative). The document
   correctly states the vanishing convention and the explicit ranges.
5. **[Not checked]** I did not check:
   - Corollary 5 (y-refinement), beyond its proof sketch, which is correct as read;
   - §8 (F_k), which is labelled conjectural;
   - the "x = 1/q" side remark;
   - Meurman–Primc, Primc–Trupčević, and Bressoud's papers.

## 5. Prior art

I read the abstract and searched the full TeX of each paper for "KRTW", "2404.03851", "Capparelli",
"Tsuchioka" and the multisum.

| source | content relevant to (i) / Conjecture 2.2 / Conjecture 4.9 |
|---|---|
| KRTW arXiv:2404.03851 v2 (TeX, 20 Apr 2026; Selecta Math., DOI 10.1007/s00029-026-01198-y) | (i) is still conjectural in v2 (Introduction, and after Conjecture 4.9). Conjecture 2.2 is proved by them only for k = 0, 1. k = 2 is via Theorem 4.10. I did not see the Springer typeset version; v2 is the accepted version. |
| Kanade–Lovejoy arXiv:2410.15408 (SIGMA 2025) | Does not cite KRTW. Proves the CLSXY identity, whose a = 0 case is Mel_k(1) = product for all k. Nothing on the 2k-fold sum. |
| Dousse–Jang arXiv:2603.05300 | Does not cite KRTW. Particle-motion proofs of Andrews / Kim–Yee type identities with parity restrictions, plus a CLSXY reproof. Nothing on the 2k-fold sum. |
| Warnaar arXiv:2511.17034 (SIGMA 2026) | Cites KRTW for the A.2 appendix and the CMPP outlook. Contains the k-fold evaluation (all k). Nothing on (i). |
| Ding arXiv:2611.10654 | Proves type C_n^{(1)} CMPP for all boundary data. Quotes the GOW Hall–Littlewood formula. No base-q² multisum, nothing on Conjecture 2.2/4.9. |
| Bartholdi arXiv:2611.04020 | Three-row level-two Primc case. Cites KRTW Thm 4.10. Nothing for k ≥ 3. |
| Capparelli–Meurman–Primc arXiv:2611.32425 | Three A_2^{(2)} level-5 partition identities. Does not cite KRTW. |
| Kurşungöz arXiv:2605.13096 | k = 1 only. Mentions "a conjecture on such a generating series in [KRTW]" as open. |
| Lau–Ono 2608.05480; Šikić 2606.19994; Bhattacharya–Dai–Hou–Tong 2610.41155; Russell 2306.16251; Tsuchioka 2606.31821 | Cite KRTW in passing, or not at all; nothing on (i). |
| Semantic Scholar citation list of 2404.03851 (8 entries); arXiv API searches (abs:Capparelli, abs:Primc, "Hall-Littlewood" AND "Rogers-Ramanujan", au:Tsuchioka/Kanade/Warnaar) through 2026-11-24; two web searches | No proof of (i) for k ≥ 3 found. |

**Not available to me:** MathSciNet/zbMATH, Google Scholar, and full-text search. I did not read the original
papers of Andrews 2010, Kim–Yee, CLSXY, Kanade–Russell, Meurman–Primc, Bressoud or Melzer, nor the Selecta
typeset version.

## 6. Document (B): proofs/T1-A-even.md (lower priority, partial)

**Check of the statement.** I computed the following sum from Kac's Weyl–Kac formula, which is independent of
Lepowsky–Primc:

    Σ_N q^{−2N²/(2k+1)} mult_N(q),    mult_N = weight-(kΛ_0 + Nα − nδ) multiplicities of the level-k vacuum module.

Inserting c^0_{2N} = q^{−c/24−N²/k} mult_N (Kac–Peterson normalisation), this sum is exactly
q^{c/24}·(right side of (0.3)). I compared it with a brute-force evaluation of the Nahm sum
F_{2k} = Σ q^{½nᵀC(A_{2k})⁻¹n}/∏(q)_{n_i}:
- they are **equal** for k = 1 up to q^14, k = 2 up to q^10 and k = 3 up to q^6, in fractional powers q^{1/(2k+1)};
- a deliberately wrong boson radius is detected.

So (0.3) is numerically correct for k = 1, 2, 3, at low order.

**Citations.**
- The form quoted as (LP), P_{k,Q} = q^{c_k/24} η c^0_{2Q} with c_k = 2(k−1)/(k+2) and m restricted by
  Σ j m_j ≡ Q mod k, is the standard Lepowsky–Primc / KKMM parafermion fermionic formula. It is consistent with
  the check above: (0.2) followed by (LP) reproduces my Weyl–Kac computation. I could not read Lepowsky–Primc
  (Contemp. Math. 46) itself, so the quotation is verified only up to that consistency.
- The Kac–Peterson facts as used are stated correctly: string functions are modular forms of weight −1/2, and
  c^0_{2N} depends on N mod k.
- I did not check the fusion lemma of `E8-E6.md`, nor Proposition 1's matrix bookkeeping, beyond this end-to-end
  numerical agreement.
