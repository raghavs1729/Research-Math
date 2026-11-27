# Referee report (second, independent) on `proofs/E8-E6.md`

Date: 2026-11-26. Formed without reading `REFEREE-E8-E6.md` (compared only at the end, see §5) and without using the
author's code or logs as evidence. My scripts: `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/e8ref/`
(`chain.py`, `series.py`/`series2.py`, `sanity.py`, `thm1.py`, `thm2.py`, `e6.py`, `rc.py`, `sec6.py`, `inter.py`).
Downloaded sources (text extracted with pdftotext) are in `.../scratchpad/dl/`.

## Verdict

- **Correctness: PROOF COMPLETE.** Theorem 1, Corollary 1 (all five mod-11 identities, including Sun–Wang (3.20)) and Theorem 2 are proved
  unconditionally from Lemma 1/Lemma 2, Andrews–Gordon and Rogers–Ramanujan. Corollary 2 and (0.6) (Sun–Wang (3.19)) are proved
  given Cao–Rosengren–Wang (CRW) Theorem 1.4, a published theorem that the document quotes correctly and whose status it describes correctly.
  I found no mathematical error. The only defect in a mathematical statement is a misleading "termwise" sentence in §4. It does not
  affect any result (Finding 1).
- **Novelty: NO PRIOR PROOF FOUND.** I found no proof of (3.19) or (3.20) in any source I could read (table in §4). The fusion
  lemma itself is known, and the document says so. What is new is the specific E8/E6 fusion chains. My search was not exhaustive (no full-text search of
  journals or of non-arXiv literature). The residual risk is that the argument is short and uses only known tools.

## 1. Findings, by severity

**F1 (minor; wording in an auxiliary argument). §4, lines 311–315 (the M(6,5) character translation of S₄).** The text says that the
four Jacobi-triple-product series for J_{7,15}, −qJ_{2,15}, J_{3,15} and the Rocha-Caridi numerators for (r,s) = (1,1), (4,1)·q³, (1,3)
"are termwise the same series". Read naturally (J_{7,15} ↔ N_{1,1}, −qJ_{2,15} ↔ q³N_{4,1}), this is **false**:
J_{7,15} = 1 − q⁷ − q⁸ + q²⁹ + …, while N_{1,1} = 1 − q − q²⁰ + q²⁹ + … The pairing is crossed:
N_{1,1} = Σ_k q^{30k²−k} (half of J_{7,15}) − Σ_k q^{30k²−11k+1} (half of −qJ_{2,15}), and
q³N_{4,1} = Σ_k q^{30k²+19k+3} (other half of −qJ_{2,15}) − Σ_k q^{30k²+29k+7} (other half of J_{7,15}).
Only the sums agree: J_{7,15} − qJ_{2,15} = N_{1,1} + q³N_{4,1}. That, together with J_{3,15} = N_{1,3}, which does hold individually, is all that is needed.
I checked this exactly to q⁴⁰⁰ (`rc.py`). The conclusion q^{−1/30}S₄ = χ₀+χ₃+2χ_{2/3} is correct. The sentence should say that the
four half-series of J_{7,15} − qJ_{2,15} coincide with the four half-series of N_{1,1} + q³N_{4,1}. The same imprecise sentence is in
the earlier referee report D4, from which the document took it. This affects only the restatement in CFT language, not Corollary 2.

**F2 (cosmetic). Line 102.** The remark that "(1.1) … is the c→0 limit of the q-Chu–Vandermonde sum" uses c both as the ₂φ₁ lower
parameter and as the summation index. The statement is correct (₂φ₁(q^{−a},q^{−b};0;q,q) = q^{−ab}), but it is confusing. It is not used.

**F3 (cosmetic). Lines 33, 40 and 61.** The letter i is reused as a bound index inside statements that are indexed by i.

**F4 (presentation, not a gap). Lemma 2, line 112.** Lemma 2 is stated for admissible (A,B) and does not assume positive definiteness, which is the right
generality. For the record, I checked that all intermediate matrices in both chains are positive definite anyway (`inter.py`). So even a
reader who only accepts the positive-definite case of the lemma (Cao–Wang's formulation) gets a complete proof.

**F5 (side material, correctly labelled).** §6 items 1–2 are sketches and item 3 is a conjecture, as the document says. I checked numerically
only (`sec6.py`): item 1 for r = 4, 5, 6 mod q⁴⁰, and the E7 conjecture mod q⁴⁰. Both agree. Nothing in §6 is used for §§0–4.

No error was found in any printed matrix, vector, quotation or citation locator that I could check.

## 2. Line-by-line check of the proof

- **§0.** The admissibility and summability set-up is correct. Q-values lie in (1/D)Z, because ½nᵀAn = Σ(A_ii/2)n_i² + Σ_{i<j}A_ij n_i n_j.
  Fusion does not change D, because A°_{00}/2 = A_aa/2 + A_bb/2 + A_ab.
  The Dynkin labelling as drawn (8 attached to 5 in E8, 6 attached to 3 in E6) matches the printed inverse matrices.
  The modularity line is correct: q^{1/88}(q⁵,q⁶,q¹¹;q¹¹)_∞ = Σ(−1)^n q^{11(n+1/22)²/2}, and 1/88 − 1/24 = −1/33.
  The M(2,11) data are correct: c = −232/11, h_i = −(i−1)(10−i)/22, h₅ − c/24 = −1/33, and the effective weights h_i − h₅ are 10/11, 6/11, 3/11, 1/11, 0 for i = 1…5.
  These match Yoshida's (3.79)–(3.80) pairing with B^{(1)}…B^{(4)}.
- **Lemma 1.** I checked each step. Euler gives the sum over c. Summing over j gives 1/(yq^i;q)_∞ = (y;q)_i/(y;q)_∞. The q-binomial theorem gives (xy;q)_∞/(x;q)_∞.
  The product is 1/((x;q)_∞(y;q)_∞). Extracting coefficients is legitimate because the coefficient of x^a y^b is a finite sum.
  Cao–Wang (1.5) cite the identity from Andrews, *Theory of Partitions*, p. 20.
- **Lemma 2.** (1.2) is right: the +1 at (a,b) produces m_a m_b, the new row produces c(Av)ᵀm, and A°_{00} = vᵀAv. The map (m,c) ↦ m + cv is a bijection
  from {(m,c) : m + cv = n} onto c ∈ [0, min(n_a,n_b)]. The other coordinates are fixed, and m_a = n_a − c, m_b = n_b − c ≥ 0.
  So the fibres are exactly the summation range in (1.1). Q° ≥ Q∘(fibre map) gives admissibility, and a summable family may be regrouped along finite fibres in Z((q^{1/D})).
  The formula uses the **old** A_ab in A°_{a0} = A_aa + A_ab, consistent with A° = PᵀAP + E_ab + E_ba. The displayed tables use this convention (checked).
  The formal-variable version is fine: each y-monomial has finitely many contributing n.
- **§2 (E8).** I recomputed the chain (`chain.py`, exact sympy). Starting from 2min(i,j) on nodes (1,7,2,8) and fusing 4:=7+8, 6:=1+8, 5:=1+4, 3:=1+2 gives
  exactly C(E8)⁻¹, which I computed independently from the Dynkin diagram. All four intermediate tables in the text were checked by hand.
  B(b) = (b₁, b₃, b₁+b₃, b₂+b₄, b₁+b₂+b₄, b₁+b₄, b₂, b₄) ✓, and B(b^{(i)}) = B^{(i)} for i = 1…5 ✓.
  The refined-form exponents and w = C⁻¹ε = (1,3,4,6,7,5,2,4) ✓.
- **§3.** The Andrews–Gordon identities are quoted correctly (k ≥ 2, 1 ≤ i ≤ k, modulus 2k+1, linear term N_i+…+N_{k−1}). Andrews PNAS 1974 Thm 1 and TofP Thm 7.8 are the right locators.
  k = 5 gives modulus 11. M_i+…+M₄ = Σ_j max(0, j−i+1)m_j ✓.
- **§4 (E6).** Chain 2:=5+6, 4:=1+6, 3:=1+2 from (4/3 2/3; 2/3 4/3) ⊕ (2) gives exactly C(E6)⁻¹ ✓ (table entries checked).
  B = (b₁, b₂+b₃, b₁+b₂+b₃, b₁+b₃, b₂, b₃) ✓, and the six B vectors of Corollary 2 are correct.
  RR is quoted correctly. The (Ex10) transcription from CRW (1.6a), (1.6b) under q³ → q into S₃, S₄ is correct (checked symbol by symbol against the arXiv text).
  The account of the literature status is accurate (see §4).

## 3. What I recomputed (independent code, written from the statements)

| Check | Method | Result |
|---|---|---|
| E8 matrix chain → C(E8)⁻¹, B(b), B^{(i)}, w | exact sympy | all match |
| E6 matrix chain → C(E6)⁻¹, B, six vectors | exact sympy | all match |
| Intermediate matrices positive definite | sympy | yes (all 7) |
| (0.3), i = 1…5 (E8 8-fold sum vs mod-11 product) | lattice enumeration with pruning, python-flint | **equal mod q⁴⁰⁰** (29 740–37 969 points per case); for B = 0, products i = 1,2,3,4 first differ at q¹,q²,q³,q⁴ |
| Controls | same code | AG 4×4 sums = products (q⁸⁰); a wrong B fails at q¹–q³; 2C(E8)⁻¹ fails |
| Theorem 1 for 13 vectors b (negative and fractional, e.g. (3,−2,1,−1), (1/2,−1/3,0,5/6)) | E8 sum vs AG-matrix sum | equal to q⁵⁰ |
| (Ex10) S₄, S₃ (both orders of B) | rank-2 sum vs J-quotients | equal to q¹⁵⁰ |
| Corollary 2, all six (incl. (0.6)) | E6 sum vs G·S, H·S | **equal to q¹⁴⁷** (14 624–17 970 points); control G·S₃, H·S₄ fail |
| Theorem 2 for 4 random rational b | E6 sum vs block sum | equal to q⁵⁰ |
| Rocha-Caridi translation of S₄ | exact series to q⁴⁰⁰ | sums equal, individual pieces cross-paired (F1) |
| §6: D_r (r = 4,5,6), E7 conjecture | series to q⁴⁰ | agree (numerical only) |

All numerical checks are only numerical. The proof does not rely on them.

## 4. Prior art

| Source | Read status | What it contains (relevant to E8/E6 with C⁻¹) |
|---|---|---|
| Sun–Wang, arXiv:2604.00847v1 | full text read (relevant sections) | A(T1,Y) = C(T1)⊗C(Y)⁻¹ = C(Y)⁻¹. (3.19), (3.20) are introduced with "We find these relations by checking the q-series to very high orders" and "We checked that…". **No proof.** Table 1 cites only their own equations. |
| Cherednik–Feigin, arXiv:1209.1978v4 | §§0.5, 2.5, 3.6 read | The ♭ case A♭ = A/2 (= C⁻¹ for simply-laced) appears in §3.6 only: dilogarithm values "calculated (numerically)", Table 1: L♭(E8) = 80/11, c♭_eff = 8/33; L♭(E6) = 24/5. (With L♭ = r − c_eff these give c_eff = 8/11 and 6/5, consistent with M(2,11) and with Sun–Wang.) There is **no q-series or product identity for E8♭/E6♭**, proved or conjectured. Separately, around (2.20) they prove modularity of the B = 0 sum of Zagier's Example 10 and report "a simple computer verification" of the coefficients (no bound stated). |
| Cao–Rosengren–Wang, arXiv:2303.01333v1 (JCTA 202 (2024) 105819) | intro and theorems read | Thm 1.2: S₁ = (−q³;q³)S₃, S₂ = (−q³;q³)S₄. Thm 1.3: product formulas for S₁, S₂. Thm 1.4: "(1.5) and the previously known identities (1.6) hold", via a Maple modular-identity package. They also say (1.6) "follow from the work of Cherednik and Feigin". Quoted correctly in the document. Nothing on E6/E8. |
| Warnaar–Zudilin, arXiv:1001.1571 | read | Matrices C(A_{N−1}) ⊗ T⁻¹ (Feigin–Stoyanovsky). Nothing on C(E8)⁻¹ or C(E6)⁻¹. |
| KKMM, hep-th/9301046 | §4 read | E8 with q^{mC⁻¹mᵀ} = 2C(E8)⁻¹ (Ising, M(3,4)), and C_G⁻¹ ⊗ C_n cosets (e.g. M(11,12)). Not (T1,E8). |
| Warnaar–Pearce, hep-th/9408136 | abstract and §2 grep | Proof of the E8 (2C⁻¹, Ising) identity. Different matrix. |
| DKKMM, hep-th/9303013 | grep | E8 ↔ Ising (2C⁻¹). Nothing on C⁻¹. |
| Nahm–Recknagel–Terhoeven, hep-th/9211034 | read | Dilogarithm/asymptotics. No E8 identity. |
| Kuniba–Nakanishi–Suzuki, hep-th/9301018 | grep | Parafermion string-function formulas ((E8,2) etc.). Not C(E8)⁻¹ Nahm sums. |
| Bytsko–Fring, hep-th/9809001 | grep | Nothing relevant. |
| Zagier, *The dilogarithm function* | §3 read | Rank 1–3 tables. No Dynkin-pair E8/E6 examples. |
| Vlasenko–Zwegers, arXiv:1104.4008 | grep and Table 1 context | Rank 2. Example 10 identities verified numerically ("unable to find a proof"). |
| C.-h. Lee, arXiv:1109.3667 | abstract and intro | Torsion in the Bloch group for C(X)⊗C(X′)⁻¹. No q-series identities. (1212.xxxx and the thesis were not read.) |
| L. Wang, arXiv:2210.10748 | grep | Rank 2 (Ex10 listed as conjecture). Nothing on E8/E6. |
| Cao–Wang, arXiv:2606.13590 | intro read | (1.5)–(1.6) = Lemma 2 for general rank (lift operators). Ranks 3–4 only. No E8/E6. |
| Cao–Wang, arXiv:2412.15767, 2508.12468 | grep | Lift-dual, ranks 3–4. No E8/E6. |
| Wang–Wang, arXiv:2607.08606 | intro read | Proves (T1,D_r). Conjecture 1.1 = ADET folklore. No E-type. |
| Wang–Wang, arXiv:2512.07790 | intro read | D_k (2C⁻¹), Bailey pairs. Mentions KKMM's E6/E7/E8 (2C⁻¹) only. |
| Milas–Wang, arXiv:2301.04532 | grep | Tadpole ranks 2–3. No E8/E6. |
| Yoshida, arXiv:2608.08780v1 | §§2.3, 2.5, 3.4 read | Physics (3d mirror symmetry). Uses (T1,E8) via Sun–Wang [29]. The vacuum identity is treated as "exact" but is not proved there. The four companion identifications are "verified through relative order q⁶" and "will not be promoted … to all-order". Charges Q̃∨ = −B^{(i)} as the document says. |
| Recent arXiv (abstract search "Nahm" Apr–Nov 2026, ~25 papers; the ones citing Sun–Wang grepped: 2607.23257, 2609.08816, 2610.15587, plus 2610.00967, 2607.22904, 2611.12868) | grep | None treats (T1,E8) or (T1,E6). 2609.08816 lists the remaining open Sun–Wang cases as rank ≤ 3 cases. |
| Mizuno's papers; Ekholm–Kucharski–Longhi; KLNS 2312.14905 | **not read** by me | Mizuno: symmetrizable rank 2–3 (from Sun–Wang's summary). EKL/KLNS: the linking lemma (as cited in the document; not re-verified). |
| Web searches (5 queries: E8 Nahm modulus 11; (T1,E8); inverse Cartan E6/E8 Sun–Wang; M(2,11) fermionic E8; arXiv API abs:Nahm∧E8/E6/Dynkin/"inverse Cartan") | done | No proof of either identity found. |

## 5. Comparison with the earlier referee report (read after the above)

My verdict agrees with `REFEREE-E8-E6.md`: no mathematical error, the E8 part is complete, the E6 part is complete relative to CRW Thm 1.4, and no prior proof was found.
The document has already absorbed that report's D1–D9: prior art for Lemma 2, the CRW status, the M(6,5) argument, the modularity line, the
formal-variable version and the attribution in §6.

New in this report:
- F1: the "termwise" sentence (present in both the earlier D4 and the document) is literally wrong. The identity is true as a statement about the sum.
- Higher-order confirmation: E8 to q⁴⁰⁰ and E6 to q¹⁴⁷.
- Theorem 1 checked with fractional and negative b.
- The intermediate matrices are positive definite.
- Precise statement of what Cherednik–Feigin contain for E8♭/E6♭: dilogarithm values only.
- Yoshida does not prove the vacuum identity.
