# Referee report: `proofs/T1-BC.md` and `proofs/T1-A-even.md`

Referee 5 (hostile), 2026-11-26. Every check below was done with code written from scratch in
`/tmp/claude-0/-home-user-weighted-prime-races/3c4df68b-d770-5604-8028-1966b4f5b18d/scratchpad/referee5/code/`
(python-flint for truncated series, sympy for exact rational-function identities). The authors' scripts
and certificates were not opened. Downloads (untrusted, only grepped/read as text) are in
`.../scratchpad/referee5/dl/`.

## 0. Verdicts

| | correctness | completeness of proof | novelty |
|---|---|---|---|
| **Theorem B** ((T1,B_r), both parity parts, all r ≥ 1) | **correct** (re-derived by hand; brute force r ≤ 5 to q^120, r ≤ 8 to q^40, parts separately) | complete; Lemmas 1–3 are right as stated | **essentially none as a q-series result.** (0.2) is *verbatim* Warnaar 2003 (formula after (5.6)), reproduced as Cherednik–Feigin (2.3). The full theorem, for **every** r (not only even r), is Du–Ji–Shen–Xu Thm 1.5 (1),(3) transported by elementary rank-changing moves. The author's "for odd r … no direct transformation" is **false** (D1). Corollary B is DJSX Thms 1.1/1.2 (j = 0), with literally the same constant C = −c/24. What remains: the observation that the DJSX sums *are* the (T1,B_r) sums of Sun–Wang's table. |
| **Theorem C** ((T1,C_r) re-proof) | correct | complete | none claimed; fair. The proof is the L,M→∞ shadow of Warnaar's own. |
| **Corollary C** (SM(4,4r+6) dictionary, all i) | correct (r ≤ 5, all i, to q^40) | one gap (r = 1, i = 1 not covered by the quoted theorem), one misplaced index-shift remark in the proof (D6, D7) | plausible as a statement, but it is a four-line dissection of a theta series on top of two quoted theorems |
| **Theorem A / Corollary A** ((T1,A_{2k}) = Σ_N c^0_{2N} q^{N²/(k(2k+1))}) | **correct** (k ≤ 3 to q^40, k = 4 to q^16, k = 5 to q^8, against string functions computed from Weyl–Kac; fusion chain exact for k ≤ 8) | complete modulo two cited theorems (Lepowsky–Primc, Kac–Peterson); (LP) as used is correctly stated (checked k ≤ 7) | **not found in the literature I could read**; risk medium (see §3.4). Important context is missing from the file (D10, D11). |

**I found no mathematical error in either document.** Both survive every derivation and every numerical
test I could devise. The defects are (i) a wrong claim about the relation to prior work and two missed
prior sources, which together remove the novelty of Theorem B, and (ii) a wrong Slater label, a wrong
index shift, a small coverage gap, and assorted citation slips. The list is in §4.

---

## 1. What was checked, and how (scripts)

| script | content | result |
|---|---|---|
| `t1_matrices.py` | Sun–Wang appendix conventions: C(B_r) (rows …(2,−2),(−1,2)), D(B_r) = diag(2,…,2,1); C(C_r), D(C_r) = diag(1,…,1,2); A·D = C(Y)⁻¹D(Y) symmetric; ½nᵀADn = ¼ΣM_i² (B), ½ΣN_j² (C); C(A_{2k})⁻¹_{xy} = min(x,y) − xy/(2k+1); mᵀC(A_{k−1})⁻¹m = ΣM_j² − L²/k; central charges from (1.6) | all exact, r ≤ 8, k ≤ 8 |
| `t2_thmB.py`, `t13_extra.py` | (T1,B_r) by brute force with a generic pruned enumerator (validated against a naive box enumeration in `t2b_naive.py`); n_r even and n_r odd summed **separately** and compared with (0.2), (0.3), and the product form of (0.3) | r = 1..5 mod q^121; r = 6,7,8 mod q^41: all True |
| `t3_bailey.py` | exact rational-function identities in p = q^{1/2}: Gasper–Rahman (II.21) (n ≤ 3, generic A,b,c), its b→∞ limit (2.1) (n ≤ 5, generic A,c), (2.2), the general-A pair α_r = V_r(A;p)(−A)^r p^{r(3r−1)/2}, β_n = 1/((q)_n(−Ap;p)_{2n}) (n ≤ 6), (P1),(P2),(P3) (n ≤ 12), the pair the file calls "G(3)" (n ≤ 10); Lemma 1 applied to (P1) (once, twice), (P2), (P3), general-A pair; negative control (wrong y) fails | all True |
| `t10_misc.py` | Lemma 2, eq. (1.2), for (P1),(P2),(P3), s = 1..4, mod q^40 ((1.3) is s = 1); Slater-label test | True |
| `t4_thmC.py` | (T1,C_r) with all r+1 linear terms, r = 1..5: theta form (Thm C, C′) and Corollary C (character difference, with the exact relative power of q); positivity of the character series | all True mod q^41 |
| `strings.py` | level-k A_1^{(1)} vacuum-module weight multiplicities from the **Weyl–Kac** formula (no Lepowsky–Primc input); sanity: k = 1 gives q^{N²}/(q)_∞; N → −N and N → N+k symmetries; level-2 vacuum string 1,1,3,5,10,16,28,… | OK |
| `t6_LP.py` | (LP) as used: Σ_{m∈N^{k−1}, Σ j m_j ≡ Q (k)} q^{mᵀC(A_{k−1})⁻¹m}/∏(q)_{m_j} = (q)_∞ q^{−N²/k} mult_N(q), N ≡ Q, two representatives N per Q; leading exponent Q(k−Q)/k | k = 2..5 mod q^31, k = 6,7 mod q^21: True |
| `t5_thmA.py`, `t13_extra.py` | rank-2k Nahm sum of C(A_{2k})⁻¹ by brute force vs (0.2) and vs (0.3) with Weyl–Kac strings | k = 1,2,3 mod q^41; k = 4 mod q^17; k = 5 mod q^9: True |
| `t7_fusion.py` | (1.1) and the chain G_k → C(A_{2k})⁻¹ exactly (matrix times Cartan matrix = identity), k ≤ 8; G_k positive definite (det = 2^{k−1}/(2k+1)); §5 chain for odd rank, k ≤ 8 | exact |
| `t8_refined.py` | charge-refined Theorem A; Feigin–Stoyanovsky (A_2,T_k) comparison (see §3.3) | True, k ≤ 3 |
| `t9_djsx.py` | DJSX (1.7), (1.9) as printed, parity parts separately, ρ ≤ 5; exact matrix-level transformations (T1,B_r) ↔ DJSX for even **and odd** r | True / exact |
| `t12_ww.py` | Warnaar's product after (5.6) = R_r^{even}, r ≤ 5; Wang–Wang 2407.21725 (3.67) = R_3^{even}, (3.70) = q^{−3/4}R_3^{odd} | True mod q^61 |
| `t11_odd_and_cft.py` | T1-A-even §5: reduced rank-k form for C(A_{2k−1})⁻¹ (k ≤ 3); the (q)_∞²F_5 observation; the k = 2 "Ising × U(1)_10" example | True; see D16 |

---

## 2. `T1-BC.md`

### 2.1 Statement vs Sun–Wang conventions — correct

With Sun–Wang's appendix, A(T1,B_r)D(B_r) = C(B_r)⁻¹D(B_r) has entries 2min(i,j) (i,j < r), i (for (i,r)),
r/2 (for (r,r)), i.e. ½nᵀADn = ¼ΣM_i²; the long nodes 1..r−1 carry (q²;q²), the short node r carries (q;q).
c(T1,B_r) = (2r−1)·3/(2r+3). For C_r: AD = (min(i,j)), node r carries (q²;q²), c = 3(r+1)/(2r+3). All as in the
file. (Sun–Wang define B_r only for r ≥ 3; "r = 1" is (T1,A_1) and "r = 2" is (T1,C_2) with the nodes swapped.
Harmless, and the file says so.)

### 2.2 Lemmas 1–3 — correct

* Lemma 1 is Bailey's lemma with ρ_1 → ∞, ρ_2 = −(aq)^{1/2}; it is exactly (S2) of
  Bressoud–Ismail–Stanton §3 (citation correct). The three-line proof is right: the inner sum is
  ₂φ₁(q^{−(n−r)}, −yq^r; aq^{2r+1}; q, −yq^n), and −yq^n = c q^N/b because y² = aq, so (1.5.2) of Gasper–Rahman
  applies. Checked exactly.
* Lemma 2: the limit n → ∞ and the regrouping are legitimate (orders grow like k²/2); (1.2) checked
  numerically for s ≤ 4 on all three pairs.
* Lemma 3: the derivation from the b → ∞ limit of (II.21) is right, including the choice c = −p^{−n} in base
  p = q^{1/2} ((c;p)_r(p^{−n};p)_r = (q^{−n};q)_r, (Ap/c;p)_r(Ap^{n+1};p)_r = (A²q^{n+1};q)_r), the exponent
  bookkeeping p^{3C(r,2)+r} = p^{r(3r−1)/2}, the substitution p → −p for (P3), and the two-case choice of c for
  (P2) in base q². All three pairs hold exactly for n ≤ 12. **The label is wrong: (P3) is Slater's G(2), not
  G(3)** (D4).

### 2.3 Proofs of Theorem C and Theorem B — correct

Every step re-derived: the substitutions L_j = N′_j + μ (even), K_j = N′_j + μ (odd), M_j = 2L_j resp. 2K_j+1,
¼ΣM_j² = ΣL_j² resp. r/4 + ΣK_j(K_j+1); β̃ and β; (q;q)_{2L} = (q;q²)_L(q²;q²)_L and
(q;q)_{2K+1} = (q;q²)_{K+1}(q²;q²)_K; the prefactors (−q;q²)_∞/(q²;q²)_∞ and (−q²;q²)_∞/(q⁴;q²)_∞·1/(1−q²); the
bilateralisation (index −n−1); the 2-dissection n = 2m, 2m−1; the triple products. Corollary B's exponents
(1/8 − 1/(4r+6) = c/24 = r/4 − (2r+1)²/(8(2r+3))) are right.

Cosmetic: the bracket in (0.3) is a *single* triple product in base −q^{2r+3},
(−q, q^{2r+2}, −q^{2r+3}; −q^{2r+3})_∞ (this is how DJSX print it); for r = 3 it is even a pure product
(Wang–Wang (3.70), verified).

### 2.4 Corollary C — correct, with a wrong shift and a small gap

All sixteen-odd bookkeeping identities of §4 (T_0..T_3, N_{1,s}, exponent differences, weights, Kac labels,
agreement with Sun–Wang (3.6)–(3.13) including their B-vectors and constants C) re-derived and tested.
Two defects: D6 (the parenthetical "k → −k−1 in T_2" in the *i even* case is misplaced — no shift of T_2 is
needed there, only k → −k in the first sum of N_{1,2r+3+i}; after the shift as written the linear coefficient
of T_2 is E+j and no longer matches the second sum of N_{1,2r+3+i} term by term; the asserted identity itself
is true) and D7 (r = 1, i = 1 is asserted but Theorem C′ is quoted only for r ≥ 2).

### 2.5 Novelty of Theorem B — the file's account is wrong in three places

**(a) (0.2) is Warnaar's identity, verbatim.** arXiv:math/0110307, the unnumbered display between (5.6) and
(5.7) (L, M → ∞ limit of (5.6)):

    Σ_{n_1..n_{k−1}≥0} q^{N_1²+…+N_{k−1}²} / ( (q²;q²)_{n_1}⋯(q²;q²)_{n_{k−1}} (q;q²)_{n_{k−1}} )
        = ∏_{j≢2 (4), j≢0,±4k (8k+4)} (1−q^j)⁻¹ .

Since (q;q)_{2μ} = (q²;q²)_μ(q;q²)_μ, the left side with k = r+1 is exactly R_r^{even} in the variables
(n_1,…,n_{r−1},μ), and the product is ((−q;q²)_∞/(q²;q²)_∞)(q^{4r+4},q^{4r+8},q^{8r+12};q^{8r+12})_∞. Verified
numerically r ≤ 5. Warnaar obtains it from "the Bailey pair C(1)" — the same pair (P2). The author cites this
very paper for Theorem C and did not see the identity two pages earlier. It is moreover reproduced as
eq. (2.3) of Cherednik–Feigin, arXiv:1209.1978 (a paper Sun–Wang cite as [12]), next to (2.4) = Warnaar (5.14)
= Theorem C, with the comment that the two are the "♭" (= A/2, i.e. T1-type) counterparts of their B/C
series and the question "Are identities … from [War] somehow associated with root systems B, C?". The file
attributes (0.2) only for r = 1 (Rogers) and r = 3 (Wang–Wang).

**(b) For odd r the (T1,B_r) sum *is* directly transformable into the DJSX sum (1.9).** The file says:
"For odd r … follows from Theorem B and (1.9) (checked numerically for ρ ≤ 5), but I have no direct
transformation." There is one, and it uses the author's own Lemma S of `KRTW-k2.md` (= Wang–Wang 2402.06253
Lemma 2.1, eq. (2.4): 1/(q)_n = Σ_{m+2c=n} q^{m(m−1)/2}/((q)_m(q²;q²)_c)). Apply it to the short node of
(T1,B_r), r = 2ρ−3: n_r = s + 2c. The new short node s has diagonal entry r/2 + 1 = ρ − ½ and linear term −½;
the new long node c has row 2·row_r and diagonal 2r, i.e. it is "node r" of a long chain of length r. The
result is exactly the data obtained from DJSX's (BD, b_0 = (0,…,0,−½)) of (1.5) by the ρ−2 fusions
"node 2j+1 := node 1 + node 2j" in base q². Checked **exactly** (matrices, symmetrizers, linear terms) for
ρ ≤ 9 in `t9_djsx.py`, and numerically. The parity of n_r is preserved by all these moves, and DJSX's own
proof splits (1.7), (1.9) by the parity of n_ρ (their F_0, F_1, eq. (3.3)–(3.4)). Hence:

> **Theorem B, both parity parts, for every r ≥ 1, is a corollary of DJSX Theorem 1.5 (1), (3)**: for
> r = 2ρ−2 by ρ−2 fusions (as the file says), for r = 2ρ−3 by one application of Lemma S and ρ−2 fusions.

**(c) Corollary B is DJSX Theorems 1.1 and 1.2 with j = 0.** Their constants c_0 = (5−4ρ)/(32ρ−8) and
(7−4ρ)/(32ρ−24) are exactly −c(T1,B_r)/24 for r = 2ρ−2 and 2ρ−3, and they prove modularity on
Γ_1(128(4ρ−1)²), Γ_1(64(4ρ−3)²) with Robins' criterion — more precise than the file's appeal to Shimura.

**(d) "the proof given here is independent of theirs".** Only in the sense that half steps (S2) replace full
steps (S1) plus fusion. The ingredients coincide: parity split of n_r; C(1) for the even part; a G-pair in
base q² with q → −q for the odd part (DJSX use Wang–Wang's (1.16) = (S1)-chain on a G-pair, itself proved
with Slater's G pairs; Wang–Wang 2407.21725 do r = 3 with (G1),(G2), q → q², q → −q). BIS state
(S2)(S2) = (S1). So the half-step chain is the square root of the DJSX/Wang–Wang chain, as the file half-says.

**Not in:** Bressoud–Ismail–Stanton (their (S2)-chains start from the unit pair: Melzer's (T1,T_{2r}) sums,
Göllnitz–Gordon; no C(1)/G chain); Wang–Wang 2402.06253 (index (1,…,1,2) only: Thm 4.1 = Theorem C′);
Wang–Wang 2407.21725 (rank 3 only: (3.67) = R_3^{even}, (3.70) = q^{−3/4}R_3^{odd}, both verified).
I did not find the **odd part (0.3) as a single multisum for general r** anywhere; for even r it has the same
evaluation as Wang–Wang's (1.16) (q → q², q → −q) with half as many summations. So the honestly new items
are: the multisum (0.3) for general r (a routine (S2)-chain), and the identification of the DJSX families with
Sun–Wang's (T1,B_r).

**CFT.** Theorem B does not identify the "?" of Table 1; it gives modularity only. The file does not claim
more, but the task description ("hence modularity for (T1,B_r)") should not be read as resolving the "?".

### 2.6 Attributions for (T1,C_r) — accurate

Warnaar (5.14) published / (5.15) arXiv, "entry (20) in Slater's list", "(5.13) is a refinement of the Bailey
pair G(1)": all confirmed in the arXiv text. Wang–Wang (3.98) and Theorem 4.1 ("discovered by Warnaar
(unpublished)", k ≥ 2, 1 ≤ i ≤ k): confirmed. Their proof of Theorem 4.1 is a Bailey-chain proof ((S3),(S5) on
(G1*)), so a chain proof on a G(1)-type pair *is* written out for the companions; for B = 0 Warnaar's
iteration of (5.13) via his (2.19) is the polynomial refinement of exactly §3.

---

## 3. `T1-A-even.md`

### 3.1 Statement, G_k, fusion chain — correct

C(A_{2k})⁻¹_{xy} = min(x,y) − xy/(2k+1); c = 2k·3/(3+2k+1) = 3k/(k+2). The six entry checks of (1.1) are right.
Proposition 1: the induction is right in all five cases of x and on the diagonal; direction of the fusion
lemma (rank r → r+1, entry A_{ab} raised by 1) matches Lemma 2 of `E8-E6.md`; admissibility only needs to be
checked for G_k and is (G_k is positive definite). Exact check k ≤ 8. The use of Cauchy's identity and the
regrouping are fine. So (0.1), (0.2) are proved.

### 3.2 Use of Lepowsky–Primc; normalisation; modularity — correct as used

* Matrix and constraint: the sum is Σ q^{mᵀC(A_{k−1})⁻¹m}/∏(q)_{m_j}, i.e. Nahm matrix **2C(A_{k−1})⁻¹**, with
  Σ_j j·m_j ≡ Q (mod k). This is KKMM hep-th/9211102 eqs. (2.1)–(2.3) with n = k−1 (S^Q_{A_n} = q^{c_n/24}b^0_{2Q},
  attributed there to Lepowsky–Primc) and DKKMM hep-th/9303013 §5.1. Correct.
* **Normalisation that makes (0.3) true:** the Kac–Peterson one,
  c^0_{2N}(τ) = q^{−c/24 − N²/k} Σ_{n≥0} mult_{kΛ_0}(2N; grade n) q^n, c = 3k/(k+2), where 2N is the sl_2 weight
  (N ∈ Z, c^0_{2N} = c^0_{−2N} = c^0_{2N+2k}). Then (LP) reads P_{k,Q} = (q)_∞ q^{−N²/k}Σ mult q^n, verified from
  Weyl–Kac for k ≤ 7, and Corollary A is equivalent to the normalisation-free statement
  **F_{2k}(q) = Σ_{N∈Z} q^{−2N²/(2k+1)} · Σ_n mult_{kΛ_0}(2N;n) q^n**, verified k ≤ 5. The file never states the
  normalisation (D13).
* Modularity: ϑ_{k,Q}(τ) = Θ_{2Q,k}(τ/(2k+1)) (weight ½), c^0_{2Q} weight −½; each product has weight 0 with a
  finite-order multiplier trivial on a congruence subgroup; a finite sum of such is modular on the
  intersection. So, granting (LP) and Kac–Peterson, q^{−c/24}F_{2k} is a modular function and Conjecture 1.1
  holds for (T1,A_{2k}). The one-line justification in the file is thin but the claim is right.
* The k = 2 example and the U(1)_3 form for k = 1 (= Sun–Wang (3.24)) are right (k = 2 verified to q^40).

### 3.3 Things the file should say and does not

* **Charge-refined statement.** The proof gives more than (0.2): with the Z-valued charge
  N(n) = Σ_{x even}(x/2)n_x + Σ_{x odd}((x−2k−1)/2)n_x (2N ≡ Σ x n_x mod 2k+1),

      Σ_{n∈N^{2k}, N(n)=N} q^{½nᵀC(A_{2k})⁻¹n}/∏(q)_{n_x} = q^{−2N²/(2k+1)} Σ_n mult_{kΛ_0}(2N;n) q^n      (verified k ≤ 3).

* **Feigin–Stoyanovsky.** The pair (A_2,T_k) has the same central charge 3k/(k+2) (h(T1) = h(A_2),
  h(T_k) = h(A_{2k}), ranks 2·k = 1·2k), and its Nahm sum is the vacuum character of ŝl_2 at level k
  (Feigin–Stoyanovsky hep-th/9308079 §4.4; I verified the charge-refined version numerically, k ≤ 3). So,
  sector by sector, **(T1,A_{2k})[N] = q^{−2N²/(2k+1)} · (A_2,T_k)[N]**: the Sun–Wang "?" is the Gaussian twist
  of a 1993 character formula. This is the natural explanation of why ŝl_2 level-k string functions appear,
  and it is absent from the file.
* **Cherednik–Feigin.** arXiv:1209.1978 §3.6 poses exactly this series (their "♭-case" A♭ = ((ω_i,ω_j)) =
  C(A_n)⁻¹), lists c♭_eff = n/(n+4) "(maybe up to a common coefficient of proportionality)", and states that
  their method "cannot be directly extended to the ♭-case"; modularity there is conditional on Nahm's
  conjecture. The project's own `PRIOR-ART.md` knows this passage, but §4 of `T1-A-even.md` does not cite it.

### 3.4 Novelty of Theorem A

Not found in: KKMM hep-th/9211102, 9301046 and DKKMM hep-th/9303013 (only 2C⁻¹ and C_G ⊗ C_A⁻¹ forms; read
around all hits); Cherednik–Feigin (poses it, no formula); Feigin–Stoyanovsky hep-th/9308079 (dual partner
only); Keegan–Nahm 1103.4986; Warnaar–Zudilin 1001.1571; Genish–Gepner 1704.02467 / Hatayama et al.
math/0102113 (twisted string functions: 2C⁻¹-type again); Hatayama et al. math/9802042, Georgiev
q-alg/9504024, Kuniba–Nakanishi–Suzuki hep-th/9301018 (grep only); Sun–Wang and the 2026 papers Wang–Wang
2607.08606, Guo 2610.15587, Chen–Yin 2610.00967, Yoshida, Mizuno, Xia (grep); arXiv API queries
("Nahm sum(s)" × parafermion / "string functions" / Dynkin / symmetrizable / "inverse Cartan matrix"; listing
of all "Nahm sum" papers to 2026-11-09).

**Not checked:** A. Dunn's thesis *Nahm's conjecture in the Cartan case* (UQ 2015, Sun–Wang's [15] — the
most likely place for a prior statement; the repository page is a script-only app and I could not get the
PDF); C.-h. Lee's thesis; Zagier's *The dilogarithm function* (text not available here);
Bouwknegt–Ludwig–Schoutens and the spinon literature; Lepowsky–Primc's and Slater's original papers; any
full-text search engine. Risk: **medium**. Every ingredient (LP sum, free boson via Cauchy, fusion) is
standard and the result is a one-page corollary once G_k is written down; the physics literature on
PF_k × U(1) theories is large.

---

## 4. Numbered list of defects

Severity: **A** = invalidates a claim of the document; **B** = wrong statement or citation that must be
corrected, mathematics unaffected; **C** = gap / imprecision; **D** = typo / cosmetic.

**T1-BC.md**

1. **[A — novelty claim]** §5 Remark (2): "For odd r … I have no direct transformation." False. One
   application of Lemma S (self-fusion, `KRTW-k2.md`; = Wang–Wang 2402.06253 (2.4)) to the short node, followed
   by ρ−2 defusions, turns the (T1,B_{2ρ−3}) sum into the DJSX sum (1.9) exactly (short node: diagonal ρ−½,
   linear term −½). Exact check ρ ≤ 9. Consequently Theorem B is, for all r, a corollary of DJSX Thm 1.5
   (1),(3), and the opening summary ("for even r the two sums are related by the fusion lemma", "new at most in
   the sense that this link was not noticed") understates the dependence.
2. **[A — missed prior art]** (0.2) for every r is the unnumbered identity after (5.6) in Warnaar,
   math/0110307 (k = r+1), obtained there from C(1); reproduced as (2.3) in Cherednik–Feigin 1209.1978, where it
   is explicitly tied to the ♭-type B/C series. Not mentioned; the file attributes (0.2) only for r = 1, 3.
3. **[B]** Corollary B is presented as a consequence of Theorem B; it is DJSX Theorems 1.1 (j = 0) and 1.2
   (j = 0) (same C = −c/24, explicit level) plus the moves of item 1. "The proof given here is independent of
   theirs" should be weakened: same parity split, same Bailey pairs (C(1); G-pair in base q² with q → −q), half
   steps instead of full steps.
4. **[B — wrong citation, three occurrences]** "(P3) is Slater's G(3)", "(this is G(3))". The pair
   α_n = (−1)^n q^{n(3n+1)/4}(1−q^{n+1/2})/(1−q^{1/2}), β_n = 1/((q²;q²)_n(−q^{3/2};q)_n) (a = q) is Slater's **G(2)**
   in the transcriptions of Wang–Wang 2402.06253 (p. 469 of Slater), Wang–Wang 2407.21725 and DJSX. G(3) is
   a = 1, β_n = q^n/((q²;q²)_n(−q^{1/2};q)_n), α_n = (−1)^n q^{3n(n−1)/4}(1+q^{3n/2}) (verified to be a Bailey pair;
   consistent with Warnaar's "generalization of G(3)", which carries q^M). The mathematics is unaffected (the
   pair is derived in the file and holds exactly). I could not see Slater's paper either.
5. **[C]** The file says the half-step proof of Theorem C "I did not find written out". Wang–Wang's proof of
   their Theorem 4.1 is a written-out Bailey-chain proof on a G(1)-type pair (for i ≤ k), and Warnaar's
   derivation of (5.14)/(5.15) by iterating (5.13) is the polynomial version of §3. Not an error, but the
   sentence oversells.
6. **[D — misplaced shift]** §4, case i even: "T_1 + T_2 = −q^{(2r+3−j)/4} N_{1,2r+3+i} (k → −k−1 in T_2)". No
   shift of T_2 is needed (T_2 already has linear coefficient E−j = 2r+3+2s and constant 2r+3−j/2 =
   (2r+3−j)/4 + s/2); the substitution needed is k → −k in the first sum of N_{1,2r+3+i} against T_1. After
   k → −k−1 the summand of T_2 no longer matches the second sum of N term by term, so the hint does not
   exhibit the identity (which is nevertheless true; the shift belongs to the i odd case only).
7. **[C]** Corollary C is asserted for all r ≥ 1 and 1 ≤ i ≤ r+1, but for i ≤ r it rests on Theorem C′,
   quoted for r ≥ 2 only. The case r = 1, i = 1 (Σ q^{n²/2+n}/(q²;q²)_n = (−q^{3/2};q²)_∞, Euler) is true
   (verified) but not covered by the text.
8. **[C]** The character formula for the non-unitary SM(4,4r+6) is adopted ("the convention under which all
   the N = 1 formulas of Sun–Wang were found to hold") with a vague attribution (GKO treat the unitary series).
   Corollary C is rigorously a statement about explicit theta quotients; that these are the irreducible NS
   characters is a citation that should be made precise.
9. **[D]** (0.3): the bracket is the single triple product (−q, q^{2r+2}, −q^{2r+3}; −q^{2r+3})_∞. Corollary B:
   "Shimura Prop. 2.1" gives the transformation of the theta series; the triviality of the combined
   multiplier on a congruence subgroup is not argued (it is standard; DJSX do it properly). The "CFT ?" of
   Table 1 remains "?": no CFT is identified for (T1,B_r).

**T1-A-even.md**

10. **[C — missing context]** §4 (novelty) omits Cherednik–Feigin §3.6, where exactly this series (♭-case of
    A_n) is posed with c♭_eff = n/(n+4) and declared out of reach of their method — although the project's own
    `PRIOR-ART.md` records that passage.
11. **[C — missing context]** No mention that (A_2,T_k) (same central charge) is the Feigin–Stoyanovsky
    character formula for the level-k vacuum module of ŝl_2, nor of the charge-refined identity
    (T1,A_{2k})[N] = q^{−2N²/(2k+1)}·mult_N(q). Without this, "the same parafermions as in SU(2)_k … with
    another radius of the boson" looks like a coincidence.
12. **[B — citation]** "(in this form e.g. Kedem–Klassen–McCoy–Melzer, hep-th/9211102, §3, 'the
    Z_N-parafermion characters b^0_{2Q}')": the statement is in **§2**, eqs. (2.1)–(2.3) (with n = k−1,
    "Z_{n+1}-parafermionic", "branching functions"); the quoted phrase does not occur.
13. **[C]** The normalisation of c^0_{2Q} is not stated. The theorem is true for the Kac–Peterson
    normalisation c^0_{2N} = q^{−c/24−N²/k}Σ_n mult_{kΛ_0}(2N;n)q^n, c = 3k/(k+2) (and for no other power of q).
14. **[C]** "the '?' of Table 1 is an extension of (Z_k parafermions) ⊗ U(1)_{k(2k+1)}". F_{2k} is a Z-linear
    combination of k(2k+1) characters of the product theory with total weights s − 2s²/(2k+1) (not integers);
    it is not the character of an extended chiral algebra. "Combination of characters of", as in Sun–Wang
    (3.24), is what is proved.
15. **[C]** Corollary A is conditional on two quoted theorems (Lepowsky–Primc; Kac–Peterson modularity); the
    author did not see Lepowsky–Primc. I confirm the statement numerically (k ≤ 7) and against KKMM/DKKMM, but
    "Sun–Wang's Conjecture 1.1 holds for (T1,A_{2k})" should be flagged as "given (LP)". The multiplier
    argument is one sentence.
16. **[D]** §4 lists Chen–Yin 2610.00967 among "papers citing" Sun–Wang; it does not cite it (Guo 2610.15587
    does). §5: "between a quarter and a half of them vanish" — I get 18/61, 31/61, 31/61, 19/61 zeros in the
    four sectors (two sectors just over one half); max |coefficient| ≤ 4 confirmed. The observation was "not
    archived" by the author; it is reproducible with `t11_odd_and_cft.py`.
17. **[not checked]** "Andrews, *The Theory of Partitions*, (2.2.8)" for Cauchy's identity; "Shimura,
    Prop. 2.1"; Slater's and Lepowsky–Primc's originals.

## 5. Recommendation

Both documents are correct and can be kept as proofs. `T1-BC.md` must be rewritten as an *identification*
note: "the (T1,B_r) Nahm sum of Sun–Wang is, by fusion (r even) or self-fusion + fusion (r odd), the
Du–Ji–Shen–Xu sum (1.7)/(1.9); hence it is modular by their Theorems 1.1/1.2; its n_r-even part is Warnaar's
identity after (5.6); a direct (S2)-chain proof is given for convenience." The claim of an open case of
Conjecture 1.1 being settled should cite DJSX as the source of the modularity. `T1-A-even.md` stands as the
only item with a plausible claim to novelty; it needs items 10–15 fixed and Dunn's thesis consulted before
any claim is made.
