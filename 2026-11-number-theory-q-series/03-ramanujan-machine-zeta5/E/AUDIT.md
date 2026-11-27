# Referee audit of E/THEOREM.md (E1, E2, E3: three "new and unproven" CFs for e)

Audited 2026-11-26. THEOREM.md was not edited. Scripts and outputs are in `audit/`:
`audit_checks.py` (+ `audit_output.txt`) and `euler_family.py` (+ `euler_family_output.txt`), written from scratch.
`verify/e_checks.py` was not used.

## Verdict

**The mathematics is correct.** The Lemma and all three proofs check out step by step, both by hand and exactly in
sympy/rationals. Numerically, the recurrence values agree with the closed forms to 60 digits. **The novelty and
status framing is wrong and must be rewritten.** None of the three is a new result in substance:

* **All three are instances of a published general theorem with proof**, O. David, *On Euler polynomial continued
  fractions*, arXiv:2308.02567 (v1 Aug 2023, v2 1 Dec 2023), Theorem 1 ("Euler continued fractions", with a
  polynomial `f`). THEOREM.md's `Π` is exactly David's `f`, shifted.
* **E2 was already proved before v4 appeared.** It is an equivalence transform of `PCF[n+2, n] = 1/(2e−5)`. That
  formula is the k=1, a=1 case of Dougherty-Bliss and Zeilberger, arXiv:2004.00090 (31 Mar 2020), Theorem 1, and the
  same Table 3 lists it as "new and proven".
* **E3 is in the trivial Euler family (f = 1).** It is therefore a direct case of Euler's classical series-to-CF
  formula. It is also a Möbius/tail transform of Cohen's entry 5.10.8 (k = 1), which is stated without proof.

I found no explicit statement or proof of E1 anywhere. It follows from David's Thm 1 by a short computation
(h1 = x+3, h2 = x^2, f = x^3+3x^2+2x+2). THEOREM.md is acceptable as a self-contained explicit write-up, but it
must not be presented as settling open problems.

## Findings

1. **[Major: novelty/status] The Status, Standing and Openness paragraph (§0) is materially wrong.** It says no proof
   of E1 or E3 was located, and it presents only E3 as Euler-family. In fact, with b = −h1·h2 and
   f(x)a(x) = f(x−1)h1(x) + f(x+1)h2(x+1), all three fit David's Thm 1. `euler_family.py` solves the linear system
   for f over all factorisations of b, with deg f ≤ 6:
   * E1: (h1, h2, f) = (x+3, x², x³+3x²+2x+2), and f(n+1) = Π_E1(n);
   * E2: (−x−3, x², (x+2)(x²+5x+5)), and f(n+1) = Π_E2(n);
   * E3: (x+2, x², 1).

   Evaluating David's closed formula, `K = (f(1)h2(1)/f(0))·(1/S − 1)`, reproduces 6e/(2e−3), 3/(4e−10) and
   4e/(2e−1) to 40 digits. THEOREM.md's `U_n = (n+1)!² Π(n)` is David's `F_n = f(n)∏h2` with the index shifted, so
   the method is the same as well as the result.
   **Fix:** rewrite §0's status. Say "explicit proofs of three instances of David's Thm 1 (arXiv:2308.02567); no
   prior explicit statement found for E1." Cite the items in findings 2–3. Delete "located no proof of E1 and E3"
   and "Not independently audited".
2. **[Major: novelty] E2 is equivalent to a proven row of the same table.** The C-transform b_n/(a_n a_{n−1}) of E2
   equals n²(n+3)/((n+1)(n+3)·n(n+2)), which is n/((n+1)(n+2)). That is the C-transform of `PCF[n+2, n]`.
   Explicitly, with r_n = (n+1)(n+3)/(n+2) we get r_n(n+2) = a_n and r_n r_{n−1}·n = n²(n+3). Since r_0 = 3/2, the
   value is `3/(4e−10) = (3/2)·PCF[n+2,n] = (3/2)/(2e−5)`, checked to 40 digits.
   * `PCF[n+2,n]` is the row "1/(−5+2e) = 2 + 1/(3+2/(4+…)), a_n = 2+n, b_n = n, new and proven" of Table 3
     (n.txt lines ~730–736).
   * It is the a=1, k=1 case of Dougherty-Bliss–Zeilberger Thm 1, `[n+k : an] = a^D/((D−1)!(e^a − Σ_{s<D} a^s/s!))`
     with D = 3, which gives 1/(2e−5).
   * Cohen records the same link implicitly: 5.10.2 (= [5/2, n+1],[1/2, n]) and 5.10.5 have identical asymptotic
     constants and the same series for 1/e.

   **Fix:** say E2 is an equivalence transform of a CF already proved in 2020, and keep the explicit proof only as
   a cross-check.
3. **[Moderate: prior art] E3 and E2 appear in Cohen's dictionary.**
   * E3: Cohen 5.10.8 at k = 1 reads 2e = 3 + 6/(4 − 16/(13 − 45/(21 − …))), and E3 shares its tail from 13
     onward. Then 4 − 16/x = 6/(2e−3) gives 3 − 3/(7 − 16/x) = 4e/(2e−1), checked to 40 digits.
   * E2: 5.10.5 reads e = 5/2 + (3/4)/E2.

   Both are stated with series, but Cohen explicitly disclaims proofs ("I do not claim that I have proved the
   correctness of the CFs given here"; "Note that I do not claim to have proved the correctness of these series").
   **Fix:** cite the source concretely. The current text gives only "label 5.10.5, per ../STATUS-MAP-RM.md". The
   labels are in H. Cohen, *Continued Fractions of Polynomial Type: Theory and Encyclopedic Dictionary*,
   arXiv:2607.06581 (July 2026), section "Constants Related to Powers of e = exp(1)"; the GP database is described in
   arXiv:2409.06086. Add the 5.10.8 link for E3.
4. **[Minor: citation] Wrong title for v4.** §0 cites arXiv:1907.00205v4 under the Nature title. The v4 title (30 Apr
   2020) is *The Ramanujan Machine: Automatically Generated Conjectures on Fundamental Constants*. "Generating
   conjectures on fundamental constants with the Ramanujan Machine" is the Nature 590 (2021) title. The Nature
   Supplementary Information, SI Table 1, still labels all three rows "new and unproven" and still prints the
   E2 sign misprint.
   **Fix:** use the correct title for v4. Optionally add the SI status.
5. **[Minor: precision] The Lemma needs q_n ≠ 0 for p_n/q_n to be defined.** Here q_n = U_n U_{−1} T_n, so this
   amounts to T_n ≠ 0 for every n, not only T ≠ 0 in the limit. That holds in all three cases:
   * E1 and E3 have positive terms;
   * E2 is alternating, and its first term 1/330 dominates, so every partial sum is positive.

   **Fix:** add one line.
6. **[OK] Lemma.** W_n = −b_n W_{n−1} holds for n ≥ 1, and the j = 0 term is the empty product 1/(U_0U_{−1}).
   Both initial data are right: (s_{−1}, W_0) = (0, U_{−1}) for q and (1/U_{−1}, a_0U_{−1} − U_0) for p. The formula
   p_n/q_n = a_0 − U_0/U_{−1} + 1/(U_{−1}² T_n) follows, and it is exact for n ≤ 40 in all three cases.
7. **[OK] Recurrence for U.** The polynomial identity (n+1)²n²Π(n) = a_n n²Π(n−1) + b_nΠ(n−2) holds identically.
   U solves the recurrence directly for n = 1..60, including n = 1 with U_{−1} = Π(−1). U_n > 0 for n ≥ −1.
   (U_{−1}, U_0) is (2, 8) for E1, (10, 33) for E2 and (1, 1) for E3.
8. **[OK] Remaining steps.**
   * The products ∏(−b_i) are j!²(j+3)!/6 for E1, (−1)^j j!²(j+3)!/6 for E2 and j!²(j+2)!/2 for E3.
   * E1's telescoping identity holds with N = −(j²+3j+1)/3.
   * E2: the factorisations Π(j) = (j+3)c_{j+2} and Π(j−1) = (j+2)c_{j+1} hold, and the telescoping identity holds
     with N = −(j+3)/2.
   * The boundary terms are −1/6 for E1 and −3/10 for E2.
   * The limits of T are (2e−3)/36, (2e−5)/(60e) and e − 1/2. The T series match these to 50 digits.
   * The final algebra simplifies to the targets in sympy.
9. **[OK] Source quotation.**
   * Table 3 rows match n.txt at lines 716–722, 737–743 and 770–776. The numerators are −4, −20, −54, −112; +4, +20,
     +54, +112; and −3, −16, −45, −96.
   * The "new and unproven" definition is quoted correctly from lines 670–671. The section is "Appendix Section A".
   * The sign misprint in row 2 is confirmed. With the printed b_n = −(n+3)n², the CF equals 2.37296506…, which is
     below 3, while 3/(4e−10) = 3.43592503….

## Numerical table

Values come from the forward recurrence at n = 400, using mpmath with 60 digits.

| | PCF (forward recurrence, n = 400) | closed form | difference |
|---|---|---|---|
| E1 | 6.69372660322108358827775481645106774520919875 | 6e/(2e−3) | 5.6e−60 |
| E2 | 3.43592503917804363316130012184276753867310311 | 3/(4e−10) | 4.0e−60 |
| E3 | 2.45079934712112815793238233857772467712173806 | 4e/(2e−1) | 1.6e−60 |
| E2 as printed (b_n = −(n+3)n²) | 2.37296506030398963657750300766 | n/a | misprint confirmed |

## Prior art, by formula

| formula | status before this write-up |
|---|---|
| **E1** 6e/(2e−3) = PCF[n²+6n+7, −n²(n+3)] | **No explicit statement found outside RM sources.** It is an instance of David's Thm 1 (2023), with h1 = x+3, h2 = x², f = x³+3x²+2x+2; the instance is not written out there. It is listed as "Known" in the RM group's Ramanujan Library table, arXiv:2412.12361v2 App. F, as C[(−n³−3n²)/(n⁴+10n³+33n²+40n+14)] = 6e/(14e−21). There "Known" means it appears in the literature, here Raayoni et al. itself, not that it has been proved. |
| **E2** 3/(4e−10) = PCF[(n+1)(n+3), +n²(n+3)] | **Proved in substance before v4.** It is an equivalence transform of PCF[n+2, n] = 1/(2e−5), which is Dougherty-Bliss & Zeilberger, arXiv:2004.00090, Thm 1, posted 31 Mar 2020, and the same table's "new and proven" row. It is also an instance of David's Thm 1. It is stated without proof as Cohen 5.10.5, in arXiv:2607.06581 and the GP database of arXiv:2409.06086, together with the series 1/e = 2Σ(−1)^n/((n²+3n+1)(n²+5n+5)n!). It is listed as "Known" in arXiv:2412.12361 App. F as C[n/(n²+3n+2)] = 1/(4e−10). |
| **E3** 4e/(2e−1) = PCF[n²+3n+3, −n²(n+2)] | **A trivial-Euler-family case.** The case f = 1 is Euler's 1748 series-to-CF identity, written as David's Example 14 (2023); applying it is one line. The instance itself is not written out there. It is a tail/Möbius transform of Cohen 5.10.8 at k = 1, which is stated with a series but without proof. It is listed as "Known" in arXiv:2412.12361 App. F as C[(−n³−2n²)/(n⁴+4n³+7n²+6n+3)] = 4e/(6e−3). |

Sources I checked that do **not** treat these three formulas:
* Lu, arXiv:1907.05563, which covers two linear CFs;
* Kadyrov–Mashurov, arXiv:1912.03214, which covers linear e CFs and π;
* Lynch, arXiv:2012.12692, an expository paper on derangements;
* Yamamoto, arXiv:2403.09729, which only covers conjectures with deg b = 2 deg a;
* Brier–Naccache–Yifrach-Stav, arXiv:2211.01058, which covers ζ/π examples;
* the ramanujanmachine.com "prove our conjectures" page, which has no entry for these forms;
* web searches for the three closed forms.

None of these contains the three formulas. The Nature published version (Feb 2021, SI Table 1) still labels all three
"new and unproven".

URLs:
- https://arxiv.org/abs/2308.02567 (David, Euler polynomial CFs, Thm 1, Example 14)
- https://arxiv.org/abs/2004.00090 (Dougherty-Bliss–Zeilberger, Thm 1)
- https://arxiv.org/abs/2607.06581 (Cohen, dictionary: 5.10.2, 5.10.5, 5.10.8), https://arxiv.org/abs/2409.06086
- https://arxiv.org/abs/2412.12361 (Ramanujan Library, App. F)
- https://kaminer.technion.ac.il/files/2021/04/41586_2021_3229_MOESM1_ESM.pdf (Nature SI), https://www.nature.com/articles/s41586-021-03229-4
- https://arxiv.org/abs/1907.05563, https://arxiv.org/abs/1912.03214, https://arxiv.org/abs/2012.12692, https://arxiv.org/abs/2403.09729, https://arxiv.org/abs/2211.01058, https://ramanujanmachine.com/prove-our-conjectures/

**Caveat:** this audit used arXiv full texts (pdftotext/grep), Cohen's LaTeX source and web search. I did not search
MathOverflow or Wikipedia exhaustively, or the ramanujanmachine.com results PDFs item by item.
