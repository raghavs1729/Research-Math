# Referee report 2: `proofs/GKPS-signed.md` (rank-three signed Nahm sums of Gang–Kim–Park–Stubbs)

Referee: independent adversarial check, 2026-11-26. Files reviewed: `proofs/GKPS-signed.md`, `code/gkps_certify.py`,
`certificates/gkps_certify.txt`, `certificates/gkps_source_transcription.txt`. Source: arXiv:2411.06081 e-print
(`Nahm.tex`, v3 of 1 Oct 2025 = SciPost Phys. 19, 128 (2025)); v1 (9 Nov 2024) was also fetched and its rank-three
data (K, μ0, a, Q, claimed characters) are identical to v3 (the diffs are comments and punctuation).

## Verdict

* **Mathematics: PROOFS COMPLETE** for every example the document claims (I checked all 17, not just 8). I found no
  mathematical error. Each change of variables is unimodular and a bijection of the summation domains, and each
  classical identity is applied with the right parameters. No step is only checked numerically, provided (L4) gets a
  proper citation (Finding 6).
* **Claims about scope and literature: ERRORS** (Findings 1–3). The headline count is wrong, "all their versions
  and Wilson lines" is false, and the statement that no earlier proof exists is false for 5 of the examples.
* **Statements match the paper: YES.** For all 17 examples, K, v = μ0 − a, the sign convention
  (−q^{1/2})^{−(μ0−a)·m}, Q, the product side and the q-power prefactor all agree with `Nahm.tex`. The exception is
  three Wilson-line claims of the paper that the document omits (Finding 2).
* **Prior proof found: YES, for a subset.** 3-7 version 3, 3-9, 3-13 (both claims), 3-14 version 1 and 3-18 were
  already proved by Zagier (2007) and L. Wang (Adv. Math. 2024, arXiv:2211.04375). None of these is credited except
  3-18 (Andrews–Gordon) and 3-13 [1] (Li–Wang). For the genuinely signed examples (3-2, 3-3, 3-4, 3-5, 3-6, 3-7 v1–2,
  3-8, 3-10, 3-11, 3-12, 3-14 v2, 3-16, 3-17, 3-23) I found no earlier proof. The search was not exhaustive.

## Independent verification performed

1. **Statements against the source.** I read `Nahm.tex` directly: §3.4 lines 737–759 (formula), §4.1.3 lines
   947–1572 (examples), Appendix "Characters of rational VOAs" lines 2209–2262, and χ_F at line 837. I compared each
   table row of the document with it.
2. **Series check, written from the paper and not from the document.** Script:
   `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/gkps/ref_check.py`.
   Output:
   `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/gkps/ref_check_q150.txt`.
   The scratchpad is session-local, so copy these files if they should be kept.
   * The left side is the paper's sum with t = q^{1/2}, computed in exact integer FLINT arithmetic. The box is
     rigorous: |m| ≤ R comes from λ_min(K)|m|² − |v+2Q||m| ≤ T, plus exact pruning in m_3.
   * The right side is built from the appendix formulas: χ^{M(P,Q)}_{(r,s)} with its q^{h−c/24}, χ_F, and
     χ^{U(1)_k}_μ. The paper's prefactor q^P is applied as written. The total exponent P+h−c/24 must lie in ½Z, which
     tests the prefactors independently.
   * **Result: all 110 claims match through q^150.** The 110 claims are every claim of the 17 examples, all versions,
     the 3-12 family for a = ±1,…,±6, and 3-11 entered separately.
   * **Deliberately wrong variants, all detected (9 of 9):**
     * Wilson-line labels swapped (3-17);
     * U(1)_2 label 0↔1 swapped (3-14 v1);
     * 3-23 W_(2,−1,1) set to 0;
     * K_11 changed (3-16);
     * μ0+a used instead of μ0−a (3-5);
     * missing minus sign (3-14 v2);
     * wrong prefactor q^{−11/60} for 3-13;
     * the sign factor dropped for 3-4 v1 and for 3-23.
3. **Line-by-line check of the proofs, with a symbolic aid.** Script:
   `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/gkps/algebra_check.py`.
   Output: `.../scratchpad/gkps/algebra_check.txt`.
   * 91 SymPy identities confirm every quadratic-form rewrite and exponent bookkeeping step stated in §2:
     * the splitting off of the m_3 part as m_3(m_3−1)/2 + c·m_3 (3-2, 3-4, 3-16, 3-23);
     * the substituted remainders (3-2, 3-16);
     * the K_a form ½(aN+m_3)²+½m_1²+½m_3² and the parities;
     * the reduction of 3-3 to 3-4 v2 in (k, m_2, s);
     * the κ values and n²+i(3i∓1)/2 in 3-23;
     * the forms used in 3-5, 3-13, 3-14 (both versions, all α, β), 3-17 and 3-18.
   * (L4) was checked exactly as a polynomial identity (multiplied by (q)_{2n}) for n = 0,…,40.

## Numbered findings

1. **The count "18 of the 24 rank-three examples" is wrong; it is 17.**
   * The examples treated are 3-2 to 3-14, 3-16, 3-17, 3-18 and 3-23, which is 17 examples (3-11 included).
   * The untreated examples listed in §4 are 3-1, 3-15, 3-19 to 3-22 and 3-24, which is 7.
   * 17 + 7 = 24. Moreover, only 16 examples are treated in all their versions, because 3-2 version 3 is untreated.
   * Note also that the source itself says "we find 23 candidate rank-zero theories" while labelling 24.
2. **"With all their versions and Wilson-line insertions" is false.**
   * (a) 3-2 version 3 (osp(1|2)) is not treated. §4 admits this, which contradicts the Result paragraph.
   * (b) The paper also claims the following, which are absent from the document's table, its proofs and its
     certificate:
     * 3-2 v2: I_A[W_(1,−1,−1)] = q^{1/60}χ^{M(2,5)}_{(1,2)} (Nahm.tex l. 1002–1004);
     * 3-4 v2 and v3: I_A[W_(1,1,1)] = q^{1/60}χ^{M(2,5)}_{(1,2)} (l. 1052, 1060).
   * The transcription file has them (line 32, "same claims"); they were dropped when the table was compiled.
   * They are true (verified to q^150) and follow at once from the existing arguments:
     * 3-2 case 2: Q·m = m_1 − M does not depend on m_2 at fixed M, so (L1) still forces M = 0, which gives G;
     * 3-4 case 2: Q·m = m_1 + M, so after M = 0 the sum is Σ q^{m_1²}/(q)_{m_1} = G; case 3 follows by symmetry.
   * This is a coverage error, not a mathematical gap. Add three lines.
3. **Prior art missed: "no proof of these rank-three identities was found" (§4) is false for 5 examples (14 of the
   87 certified sums).** The examples with even v carry no sign, so they are ordinary Nahm sums, and they are on
   Zagier's rank-three list. GKPS §4.3 (l. 1737) itself says that Zagier's integral-B triples appear in their list.
   * **3-7 v3** (both claims) and **3-14 v1** (all four claims):
     * These are Zagier's Example 2 of Table 3. Its matrix is
       [[αh²+2, αh, −αh],[αh, α, 1−α],[−αh, 1−α, α]], with α = 1, h = −1 for 3-7 and α = 2, h = −1 for 3-14.
     * The vector is B = (ανh+B_1, αν, −αν). This gives ν = 0 with B_1 ∈ {0, 1} for 3-7, and
       (ν, B_1) = (−½,1), (0,1), (0,0), (½,0) for the four 3-14 Wilson lines.
     * Zagier proved this example himself (D. Zagier, *The dilogarithm function*, 2007). The proof is reproduced as
       L. Wang, arXiv:2211.04375 (Adv. Math. 2024), Theorem 2 / eq. (exam1-proof).
     * The document's 3-14 argument (Durfee-diagonal sum, then a theta series) *is* Zagier's argument.
     * 3-14 v2 is the same argument with ν = −3/4 and q^{1/2} → −q^{1/2}. This works because diag K = (4,2,2) is
       even, so the sign factor equals (−1)^{E(m)}. Minor: Wang's printed Theorem 2 drops the ν-dependence of the
       theta factor, which is correct only for ν ∈ Z; his derivation is correct.
   * **3-9**, both claims: Wang Theorem 6, (exam6-4) = φH and (exam6-1) = φG. The proof (Euler, then the q-binomial
     theorem) is the same as the document's.
   * **3-13**: I_A[1] = H² is Wang (exam10-5) and I_A[W_(2,1,1)] = G² is Wang (exam10-2). The Wilson-line case is
     therefore not new either. The document's short (L5) proof may be simpler than Wang's integral method.
   * **3-18**: Wang Theorem 9 and Andrews–Gordon; the document already says it is classical.
   * Pointers for the "not treated" list of §4, checked only at the level of matrix and B vector:
     * 3-21's K = Wang Example 11 with B = (2,1,0), (1,0,0), (0,0,0), i.e. (exam11-5), (11-3), (11-1);
     * 3-22's K = Wang Example 8 with B = (2,1,1), (0,0,0), i.e. (exam8-6), (8-5);
     * 3-15's K = C(A_3) with B = (0,0,0), (−1,1,0), (−1,1,1) is Cao–Wang, arXiv:2412.15767, Theorem 4.1, which cites
       GKPS §4.1.3.
     * For these, the product identities are proved; matching them to the SM(2,12) and osp(1|4)_1 characters is a
       separate step.
   * The truly signed cases appear to have no earlier proof (WebSearch over arXiv and the web, 2024–2026, including
     2412.15767, 2607.23257 Shi–Wang, and 2608.08780). INSPIRE / Google Scholar citation lists were not accessed.
4. **Changes of variables are correct and bijective.** I checked each one:
   * M = m_2+m_3 (3-2 c2, 3-4 c2) with 0 ≤ m_2 ≤ M;
   * N = m_1+m_2 (K_a, 3-9); N = m_1+m_3 (3-5) with 0 ≤ m_3 ≤ N;
   * (m_1, N, M) with N, M ≥ m_1 (3-13, 3-17);
   * (m_1, u, m_2) with u ∈ Z and d = u − m_1 (3-14): every d ∈ Z is reached and each diagonal fibre is infinite;
     the (L5) limit is applied correctly, including the sign (−1)^{m_2+m_3} = (−1)^d;
   * the nested sets N_3 = m_1, N_2 = m_1+m_2, N_1 = m_1+m_2+m_3 (3-18);
   * m_2 = m_1+1−ε+j (3-2, 3-16), with the complementary range killed by (E) at z = q^c, c ≤ 0.

   All are unimodular with the stated ranges.
5. **The classical identities are applied with the right parameters.** I checked each application:
   * (E) at z = q^c (vanishing for c ≤ 0; (q)_∞/(q)_{c−1} for c ≥ 1) and at z = −q^{1/2±m_1};
   * the factorisation (−q^{1/2−m};q)_∞ = φq^{−m²/2}(−q^{1/2};q)_m;
   * (L1) at z = 1 and at z = q^{1−N}, which vanishes by its factor k = N−1;
   * (L2) in 3-17, with a proof that is correct (generating function in x);
   * (C) with c ≥ 0 in every use;
   * (L3) with a = m_1 and R ∈ {0,1} in 3-16. Its imported proof (E7.md, Lemma 4) is correct; I rechecked
     j²+Rj+(m−j)(m+R+a) = m²+Rm+(m−j)(a−j).
   * (L5), which holds by the generating-function argument Σ x^Ay^B·RHS = 1/((x)_∞(y)_∞);
   * the Jacobi triple product: Σ(−1)^uq^{u²+u/2} = (q^{1/2};q)_∞(q²;q²)_∞ = (q)_∞/φ, with the shift u ↦ u−1;
   * Andrews–Gordon with k = 3 (mod 7) and k = 4 (mod 9), with the right linear terms;
   * the 3-23 antisymmetry (M odd and κ = 0). The m_1 = 0 case of Q = (4,−3,2) has M = −1 and an empty range; this
     is unstated but trivial.
   * The formal-series setting is fine: every intermediate family is bounded below by a positive-definite form, so
     regrouping and ungrouping are legitimate. The blanket remark in the document covers only the original family.
     Add one sentence saying that each re-expanded family (3-3, 3-14, 3-16, 3-17) is again summable.
6. **(L4), the finite pentagonal number theorem, is true but poorly cited.**
   * The citation "Andrews, Theory of Partitions, Example 3.10 type" is vague, and the certificate checks only n ≤ 6.
   * The identity is the standard Bailey pair β_n = 1/(q)_n, α_n = (−1)^n q^{n(3n−1)/2}(1+q^n) relative to a = 1.
     It is one step of Bailey's lemma applied to the unit Bailey pair, the same pair Rogers used for RR. E7.md
     already uses this lemma.
   * Cite it that way. I verified (L4) exactly for n ≤ 40.
7. **Smaller inaccuracies.**
   * The half-index formula is in GKPS §3.4, not "§3.3" (§3 has the subsections 3.1 matter theories,
     3.2 conditions, 3.3 simple lines, 3.4 half-index).
   * Everything else in the source summary (the quotation from §4.1, the search range ±17, the SciPost reference,
     §4.1.3) is accurate.
8. **The certificate does what the text says, with caveats.** Arithmetic is exact (Python ints and `Fraction`
   exponents). The caveats:
   * (a) The recorded `certificates/gkps_certify.txt` was run with N = 23, so it compares only up to q^22. The
     document does not state the precision. My rerun at N = 30 passes all 87 sums and the lemma checks.
   * (b) The truncation box is heuristic. B grows until all *integer* points on the faces m_i = B have exponent > N,
     which does not rigorously exclude points outside the box. In practice it agrees with my rigorous enumeration.
   * (c) The right sides are product forms (G, H, AG_i, φ, ϑ). The paper's q-prefactors are not certified. They are
     correct: I checked them independently through q^150.
   * (d) The lemma checks are finite ((L1) M ≤ 7, (L2) N,M ≤ 6, (L3) a ≤ 4, (L4) n ≤ 6, (L5) A,B ≤ 5). This is
     consistent with the file being labelled "checks, NOT proofs".
   * (e) The three Wilson lines of Finding 2 are missing.

## Examples checked / not checked

* **Checked in full** (statement against the paper, my own series to q^150, line-by-line proof):
  3-2 (v1, v2, including v2 W), 3-3, 3-4 (v1–v3, including the v2 and v3 W), 3-5, 3-6 (three versions), 3-7 (three
  versions), 3-8, 3-9, 3-10 (v1, v2), 3-11 (v1, v2), 3-12 (the proof for all a; series for a = ±1,…,±6), 3-13, 3-14
  (v1, v2), 3-16, 3-17, 3-18, 3-23. This covers every proof type: E+C, L1-delta, L1 with z = q^{1−N}, Euler+L1,
  E+L4+Durfee, L5, the L5 limit with JTP, E+L3+AG, L2+AG, and the vanishing cases.
* **Not checked:** 3-1, 3-2 v3, 3-15, 3-19, 3-20, 3-21, 3-22, 3-24. The document does not claim these; I only gave
  literature pointers for 3-15, 3-21 and 3-22.
* Also not checked: Zagier's 2007 text itself (I relied on Wang's reproduction of his proof); the full texts of
  Shi–Wang 2607.23257 and Milas–Wang (the dual of Example 6, possibly relevant to 3-1 = C(T_3)); INSPIRE citation
  lists.

## Required corrections (summary)

1. Change "18" to "17". Say that 3-2 v3 is excluded.
2. Add the three missing Wilson-line claims (3-2 v2, 3-4 v2, 3-4 v3), with the one-line proofs from Finding 2.
3. Rewrite the literature paragraph:
   * credit Zagier 2007 / Wang 2024 for 3-7 v3 and 3-14 v1 (Zagier's Example 2), for 3-9 (Wang Theorem 6) and for
     3-13 including the Wilson line (Wang Theorem 10);
   * note that 3-14 v2 is Zagier's argument with signs;
   * update §4 with the Wang and Cao–Wang pointers for 3-21, 3-22 and 3-15;
   * restrict the novelty claim to the signed cases.
4. Give a proper citation (Bailey pair) for (L4). Add the summability remark for the re-expanded families. State the
   certificate precision, or rerun it to a higher order.
