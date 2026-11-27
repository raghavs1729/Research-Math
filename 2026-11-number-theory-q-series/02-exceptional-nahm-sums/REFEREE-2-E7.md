# Second independent referee report on `proofs/E7.md` (Theorem 4: C(E7)⁻¹ Nahm sum; Theorem 5: the k-family)

Refereed version: `proofs/E7.md` at commit 5eac89d (192 lines, including §6 "Addendum"). Lemmas 1–2 of
`proofs/E8-E6.md` were read (statements and proofs) and taken as given. `REFEREE-E7-KRTW.md` was opened only
after the verdicts below were fixed. My scripts are in the scratchpad (`e7ref/`: `fusion.py`, `nahm.c`, `run.py`,
`inter.py`, `series.py`). They were written from the statements, not from the project's code.

## Verdicts

| | Verdict |
|---|---|
| **Theorem 4 (E7)** | **PROOF COMPLETE.** Every step was re-derived by hand and re-computed. I found no mathematical error and no gap that affects validity. The remaining findings are editorial, attribution or reproducibility issues. |
| **Theorem 5 (family, all k ≥ 1)** | **PROOF COMPLETE, with one gap in the write-up (F1).** The proof is a sketch, but every k-dependent step was checked symbolically for general k and holds. The general-k Bailey iteration (Lemma 3 with k−1 inner variables) is used but not stated. It follows by the same three-line induction as Lemma 3. Nothing else depends on k: c_n does not depend on k at all. |
| **Novelty (Theorem 4)** | **NO PRIOR PROOF FOUND** (see the limits in §C). The source of the problem lists (T1,E7) with CFT "?" (Sun–Wang, Table 1). Cherednik–Feigin 2013 give only the dilogarithm value. |
| **Novelty (Theorem 5)** | **k = 2 is known; for k ≠ 2, NO PRIOR PROOF FOUND.** k = 2 is Melzer's (T1,T₄) identity (Melzer hep-th/9412154, conjecture), proved by Bressoud–Ismail–Stanton, *Change of base in Bailey pairs*, Thm 5.1, and Warnaar 2003, Thm 4.4 (as cited in Sun–Wang §2(2)). The document says "a known case of the Melzer identities" but gives no reference. For k = 3 the family is Theorem 4 itself. For other k I found no source. |

## Findings, by severity

None of the findings affects validity.

**F1 (gap in the write-up, Theorem 5; §6, lines 187–189).** The proof of Theorem 5 uses "iterating Bailey's lemma
k times from the unit pair (Lemma 3 with k−1 inner variables)". Lemma 3 is stated and proved only for three
iterations. The general statement is:

Σ_{N≥r_2≥…≥r_k≥0} x^{r_2+…+r_k} q^{r_2²+…+r_k²}/((q)_{N−r_2}(q)_{r_2−r_3}⋯(q)_{r_k}) = Σ_{r≤N} x^{kr}q^{kr²}α_r/((q)_{N−r}(xq)_{N+r}).

It follows verbatim from (i) and (ii) of the proof of Lemma 3 by induction. For k = 1 there is no inner variable and T_1(a) does not depend on a, but the formula still holds. It is also exactly the finite form of Andrews 1984, Theorem 1 / (3.1), with the unit pair.

*What closes it:* state Lemma 3 for general k (one line) and note the k = 1 and k = 2 degenerate cases.

The remaining k-dependent algebra is correct for all k. I checked: (k/4)a² + kar + (k+1)r² + r(r−1)/2 + j² + 2rj = (k/4)(a+2r)² + r(r−1)/2 + (r+j)². Hence c_n is the k-independent polynomial (3.1), Proposition 6 applies unchanged, and (q)_∞·LHS = 1 + Σ_{n≥1}(q^{kn²/4} + q^{(k+4)n²/4}).

**F2 (internal inconsistency, editorial; §5, second bullet, lines 172–173 vs §6).** §5 says the extension to longer
chains "is not pursued here". §6 then pursues it and claims a theorem.

*What closes it:* delete or update the bullet.

**F3 (missing attribution, Theorem 5 k = 2; line 184).** The claim "a known case of the Melzer identities" is correct (I verified numerically to q^150 and checked the defusion step myself) but has no citation.

*What closes it:* cite Melzer, hep-th/9412154, (1.7)/(2.1)–(2.7) (conjecture); Bressoud–Ismail–Stanton,
Ramanujan J. 4 (2000), Thm 5.1; Warnaar, Discrete Math. 272 (2003), Thm 4.4. With these citations, Theorem 5 at k = 2 is a new proof of a known identity, not a new identity.

**F4 (reproducibility; §0 line 6, §4 line 163, §6 line 192).** The logs `certificates/certify_E7.log` and
`certificates/test_family_k.log` are cited, but no `certificates/` directory exists under `monotonicity/rr-identities/`. Neither log is tracked anywhere in the repository (checked with `git ls-files`). `code/certify_E7.py` and `code/test_family_k.py` do exist.

*What closes it:* commit the logs, or change the references to the scripts.

**F5 (cosmetic, normalisation not explained; lines 13–14).** The document says "c♭_eff = 1/3 for E7, i.e. c = 1 in the present normalisation" but does not explain the conversion. In Cherednik–Feigin, c♭_eff = L♭/h with L♭(E7) = 6, h = 18. The Nahm effective central charge is rank − L♭ = 7 − 6 = 1, which equals 3c♭_eff. Both formulas give Sun–Wang's values 6/5 and 8/11 for E6 and E8. The statement is correct.

*What closes it:* one clause stating the conversion.

**F6 (cosmetic; Lemma 3 proof, line 93).** The parenthesis "(the referee compared the statements … verified (2.1)
exactly for N ≤ 6)" refers to an earlier review. I independently confirmed the citation locators in the current text:
- Andrews, Pacific J. Math. 114 (1984): (2.3) is the Bailey-pair definition, (2.4) is Bailey's lemma, (2.5)–(2.6) are the new pair, Theorem 1 is (2.7), the unit pair is (2.12)–(2.13), and Theorem 2 (3.1) is the limit N, b_i, c_i → ∞. I read these from the scanned pages.
- Warnaar, arXiv:0910.2062: (1.3) is the definition and Theorem 2.1/(2.2) is the ρ₁, ρ₂ → ∞ case α'_L = a^L q^{L²}α_L, β'_L = Σ a^r q^{r²}β_r/(q)_{L−r}.

The document's (ii) matches both. I verified (2.1) symbolically in x for N ≤ 4.

**No error found in:**
- §1: the fusion chain. I computed C(E7)⁻¹ from the diagram (det C = 2) in exact arithmetic; all three stage matrices match the printed ones entry by entry, and the final reordered matrix is C(E7)⁻¹. A⁽⁰⁾ gives ¾a² + a(N₂+N₃) + ΣN_j².
- The diagonal of C⁻¹ and the parity claim ½nᵀC⁻¹n ≡ ¾(n₁+n₃+n₇)² mod 1 (checked by hand and for all n ∈ {0..3}⁷), and both refinements.
- The unit pair and the x = 1 convention (removable singularity, (x)_r/(1−x) = (xq)_{r−1}), and the ρ₁, ρ₂ → ∞ limits.
- The interchange in Proposition 5: the (N,r) family is q-adically summable, since the exponent is ≥ N² and α_r ∈ q^{r(r−1)/2}Z[[q]].
- Lemma 4: the Cauchy identity derived from Lemma 1 with A → ∞ coefficientwise, the exponent identity j²+Rj+(m−j)(m+R+a) = m²+Rm+(m−j)(a−j), and the inner sum by Lemma 1.
- Proposition 5: (q^{a+1})_{m+2r} = (q)_{m+2r+a}/(q)_a; α_r(q)_a, including the a = 0 case and the a = r = 0 case; the exponent (9r²−r)/2.
- §3: the regrouping n = a+2r (summable in Z[[q^{1/4}]], exponent ≥ ¾n², finitely many (r,j) per n), (q)_{a+r−1} = (q)_{n−r−1}, and j ≤ a ⇔ 2r+j ≤ n.
- Proposition 6 for all n:
  - the substitution u = r+j;
  - the three reversal formulas, whose q-powers combine to q^r and whose signs combine to +1 (re-derived);
  - the ₂φ₁ with parameters (q^{−u}, q^{u−n}; q^{1−n}; q, q). It terminates at u, the terms with r > n−u vanish because (q^{u−n})_r = 0, and the denominators (q^{1−n})_r ≠ 0 for r ≤ u ≤ n−1;
  - q-Chu–Vandermonde in the form G–R (1.5.3) ₂φ₁(q^{−u},b;c;q,q) = (c/b)_u b^u/(c)_u. This is a polynomial identity after multiplying by (c)_u, so specialising c = q^{1−n} is legitimate;
  - c/b = q^{1−u}, and (q^{1−u})_u = 0 for u ≥ 1;
  - the boundary terms S₀ = S_n = 1/(1−q^n).

  So c_n = 1 + q^{n²} for every n ≥ 1. This is a proof, not a check.
- The final assembly and the factor ½: Σ_Z q^{3k²/4} + Σ_Z q^{7k²/4} = 2 + 2Σ_{k≥1}(…), so the k = 0 terms give exactly the "1".
- The modularity paragraph: q^{3k²/4} = e^{πik²(3τ/2)}, and the q^{−1/24} with c = 1.
- c(T1,E7) = 7·3/21 = 1 (Sun–Wang's formula with h(T1) = 3).
- The remark ½[(A1,D₃)+(A1,D₇)]: checked numerically to q^100.

## What I recomputed (own code, written from the statements)

1. **C(E7)⁻¹ and the three fusion moves**, exact (sympy). Every stage matrix agrees, and the final matrix equals C(E7)⁻¹.
2. **Theorem 4 by direct enumeration of N⁷**, in the variable t = q^{1/4}, exact integer arithmetic (C program; each 1/(q)_n is built incrementally, so the arithmetic is exact): **agreement through q^200** (t^800; 22 058 lattice points; the coefficient of q^200 is 11 365 684 153 082).
   - Negative controls (one matrix entry changed; a wrong right side) fail at once.
   - The reduced 4-fold sum F_{A⁽⁰⁾} (1.1) also agrees through q^200.
3. **Theorem 5 by direct enumeration** for k = 1, …, 8: agreement through q^150 for every k (the document checks k ≤ 6 mod q^80).
   - The k = 2 matrix is one fusion away from C(T₄)⁻¹ = (min(i,j)). The (T1,T₄) sum agrees with the k = 2 right side through q^150.
4. **Intermediate identities:**
   - Lemma 3 (2.1): exact as rational functions of x, q for N ≤ 4. The unit pair is exact for n ≤ 4.
   - Lemma 4 for a, R ≤ 4 mod q^71.
   - Proposition 5 against the direct triple sum T(a) for a ≤ 6 mod q^71.
   - c_n = 1 + q^{n²} exactly (sympy `cancel`) for n ≤ 12, and mod q^{n²+n+21} (mod 2³¹−1) for n ≤ 30.
   - S_u = 0 exactly for 1 ≤ u ≤ n−1, n ≤ 8.
5. **KKMM D_r remark:** Σ q^{nᵀC(D_r)⁻¹n}/∏(q)_{n_i} = (q)_∞⁻¹Σ_Z q^{rk²/4} for r = 3, 7 to q^100.

## C. Prior art

Method: the PDFs were downloaded and searched as full text (pdftotext; grep for E7/E_7/"(T1 ,"), with targeted
reading. I also listed all arXiv abstracts matching "Nahm sum" (API, sorted by date, 55 hits through 2611.12868) and ran web searches for '"E7" Nahm sum inverse Cartan', '"(T1, E7)"', and 'q^{3n²/4} q^{7n²/4} fermionic c = 1'.

| Source | Read? | What it contains relevant here |
|---|---|---|
| Sun–Wang, arXiv:2604.00847 | full text | Table 1: "(T1,E7) 1 ? −" (central charge 1, CFT "?", no reference or formula). Elsewhere it says only "for Nahm sums of type (E6,T1) and (E7,T1), we have not found good CFT interpretations". (2.2): the (T1,T_{2r}) identity (Melzer; proofs [6] = BIS Thm 5.1 and [51] = Warnaar Thm 4.4). That covers Theorem 5 at k = 2. A(X,Y) = C(X)⊗C(Y)⁻¹, so A(T1,E7) = C(E7)⁻¹. |
| Cherednik–Feigin, arXiv:1209.1978 (Adv. Math. 2013) | full text §3.6 | A♭ = A/2 with entries (ω_i,ω_j) = C⁻¹ for simply-laced types; Table 1 gives E7: L♭ = 6, c♭_eff = 1/3. Quote: "Our approach … cannot be directly extended to the ♭-case". No q-series formula for E7 anywhere (E7 occurs only in the table and in "E7 ↔ A2"). The document's description is accurate. |
| KKMM, hep-th/9211102; DKKMM, hep-th/9303013 | full text, E7 sections | They use B = C(E7)⁻¹ in Σ q^{mBm} (that is, matrix 2C⁻¹) with the parity restriction m₁+m₃+m₆ ≡ Q, which gives c = 7/10 (tricritical Ising). They make no remark on the halved matrix C⁻¹. |
| hep-th/9301046 (KKMM) | full text | No E7. |
| Wang–Wang arXiv:2607.08606; Wang–Wang 2512.07790 | full text | Only types D_r and (T1,T_r); 2512.07790 mentions E7 only when citing KKMM's 2C⁻¹. |
| Cao–Wang arXiv:2606.13590; 2508.12468; 2412.15767 | full text | Ranks 3–4 and lift operators; no E7, no (T1,·). |
| Mizuno arXiv:2609.09653 | full text | Kanade–Russell mod 9; no E7. |
| Xia arXiv:2609.08816 | full text | Lists the open rank ≤ 3 cases of Sun–Wang ((T1,G2), (G2,T1) are now proved); no E7. |
| Yoshida arXiv:2608.08780 | full text | (T1,E8) and (E8,T1); E7 appears only as (E7)₁×SU(2)₁ ⊂ (E8)₁. |
| 2611.12868, 2610.15587, 2610.00967, 2609.24129, 2609.28509, 2609.14940, 2607.23257, 2607.22904, 2606.25866, 2602.14483, 2608.08347 | full text, grep only | No occurrence of E7 / (T1,E7). |
| Warnaar, "50 years of Bailey's lemma", arXiv:0910.2062 | §§1–2 read | Statements used in Lemma 3 confirmed; no identity of this shape found by grep (not read in full). |
| Andrews, Pacific J. Math. 114 (1984) | pp. 268–273 read (scan) | Equation numbers confirmed (F6). Theorem 2 (3.1) is the Andrews–Gordon chain relative to a; it does not contain the coupling to Σ_a q^{γa²}/(q)_a. |
| Bressoud–Ismail–Stanton, arXiv:math/9909053 | Thm 5.1 and its proof read | The (T1,T_{2r}) family with an extra parameter a (base-change lemmas); this covers k = 2 of Theorem 5. No E7, no general-k family of this form. |
| Agarwal–Andrews–Bressoud, *The Bailey lattice*; Warnaar's A2-Bailey papers; Berkovich–McCoy–Schilling (N=2 / c = 1 fermionic) | **not read** (no full text obtained) | Unknown. In my judgement the Bailey lattice produces linear-term variants of Andrews–Gordon, not a quadratic extra variable, but this is not verified. |
| Melzer hep-th/9412154 | not read directly (search hit only) | Cited via Sun–Wang as the source of the (T1,T_{2r}) conjectures. |

Assessment:
- **Theorem 4:** no prior statement or proof found in the sources read. The source paper of the problem (2026) still marks it open.
- **Theorem 5:** k = 2 is known (Melzer, BIS, Warnaar). k = 3 is Theorem 4. k = 1 is the product (Σ q^{a²/4}/(q)_a)·G(q) of two known modular series, and I did not search for its theta-function form specifically. No source found for general k.
- The search did not cover the unread items in the table above or the non-arXiv literature (journals not on arXiv, theses).

Earlier report: after fixing my verdicts I read `REFEREE-E7-KRTW.md`. It agrees on the correctness of Theorem 4. Its E1 (wrong locators in Andrews 1984) has been fixed in the current text. It does not cover §6 (Theorem 5). F1–F4 above are new.
