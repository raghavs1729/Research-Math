# Referee report: proofs/X-A-parafermion.md

Referee X-A, 2026-11-26. Independent scripts in `code/xa_referee/` (no project imports). Work log appended as I go;
final verdicts and summary at the end.

## Work log

### Step 1 (quotations), first pass
- Sun–Wang main.tex l.145–163 (generalized Nahm sum, "modular quadruple", Conjecture 1.1): note §0 quotation agrees
  verbatim (modulo the elided "where ⊗ is the Kronecker product"). Appendix l.699–836: SW's C(X) is the *transpose*
  of Kac's convention (C_{ij} = 2(α_i,α_j)/(α_j,α_j)): e.g. C(B_r)_{r−1,r} = −2 with α_r short, C(G_2)_{21} = −3 with
  α_2 long. With this convention C(X)D(X) is symmetric and equals the Gram matrix with short roots of square 2,
  so G = C(X)D(X) in (0.1) is correct and only G, D enter. Checked by hand for B_r, C_r, F_4, G_2, T_r.
- BKP main.tex: `thm38` l.1866–1872, K_P l.1859, D'_P and G^{mn} (charGij2) l.866–872, rect l.515, parafermion2 l.1400,
  ω_Ω l.1717–1722, ch = q^{−c_Λ} tr q^{L_Ω(0)} l.1806, [L_Ω(0), e_{β∨}] = 0 l.1775: all quotations accurate.
  c_{kΛ_0} = 0 is not stated by BKP but is immediate (L(0) and L_h(0) kill the vacuum).
- OT main.tex: last theorem l.1343–1350 quoted accurately; ρ_i l.371–373; ĝ[ν] l.487–497; D l.1215–1219;
  coinvariants l.1059–1071. Minor: OT never define (q^{ρ_i})_n; the reading (q^{ρ_i};q^{ρ_i})_n is forced by (C1)
  m_{p,i} ∈ ρ_i Z (l.685). OT's character definition (l.1275–1277) says "conformal energy −m" with weight q^m — a sign
  slip in OT (the proof uses +energy); harmless.
- Takenaka main.tex: `character2` l.1320–1325 quoted accurately; the misprint ∏(1−q^{1/2}) (l.1258–1261) is real,
  and the note's reading via l.1264 is right; D l.1279; `ST ch` l.1331–1342.

### Step 1 (dictionary), by hand
- Projections (OT labelling, l.150–177): A_{2l−1}: orbits {i,2l−i} (i<l, orthogonal pair) → square 1, fixed α_l → 2,
  ⟨(α_{l−1})_0, α_l⟩ = −1: r·Gram = C_l Gram with α_l long, D = (1,…,1,2) = SW's D(C_l). D_{l+1}: α_1..α_{l−1} fixed
  (square 2), ½(α_l+α_{l+1}) square 1 → B_l, D = (2,…,2,1) = SW's D(B_l). E_6: {1,6},{2,5} short, 3,4 long → F_4,
  D = two 2's, two 1's. D_4: ⅓(α_1+α_3+α_4) square 2/3, α_2 square 2, inner product −1 → r·Gram = [[2,−3],[−3,6]]
  = C(G_2)D(G_2) of SW. Hence r·ρ_i = d_i and q_OT = q^r turn OT's theorem literally into (0.1). No transposition
  issue arises because only G = CD and D enter the sum (SW's C is the transpose of Kac's, see above).
- Fixed-point algebras C_l, B_l, F_4, G_2 and horizontal subalgebras of A_{2l−1}^{(2)}, D_{l+1}^{(2)}, E_6^{(2)}, D_4^{(3)}
  agree (Kac Table Aff 2/3, node 0 with a_0 = 1). The B ↔ D_{l+1}^{(2)}, C ↔ A_{2l−1}^{(2)} assignment is confirmed
  numerically: cross-wiring (A_5^{(2)} strings against the (B_3,A_1) sum and vice versa) fails in every sector.
- Imaginary roots of OT's ĝ[ν]: h_{(j mod r)}⊗t^{j/r}; δ ↔ energy 1/r, so e^{−δ} = q_OT^{1/r}, q_OT = q^r. Central element:
  with ⟨x_θ0,x_{−θ0}⟩ = 2/⟨θ_0,θ_0⟩ one finds α_0^∨ = c − θ_0^∨, so c is Kac's canonical K (level k is Kac's level).
- Takenaka: order-4 lift; real-root modes in ¼+½Z (α_l, α_{l+1}), Heisenberg modes in ½Z; δ ↔ energy ½, q_Tak = q²;
  (α_i)_0 = ½(α_i+α_{2l+1−i}) has Gram ½C(T_l). Consistent with Kac's normalised form ((α_0|α_0) = 1, (α_i|α_i) = 2,
  (α_l|α_l) = 4) and with α_0 = ½δ − ε_1.
- Minor (§0, l.53–54): "the exponents of F_{X,k,c} lie in −(μ|μ)/(2k)+Z … G is even" is false for X = T_l (G_{ll} = 1):
  for T_l the exponents lie in −(μ|μ)/(2k) + ½Z and (μ|μ)/(2k) is only defined mod ½Z on a class when k is odd. Harmless
  (§5 sums over n ∈ ½Z), but should be restricted to the a_0 = 1 columns.
- Minor (§5, Theorem 3): for μ ∈ Q̄_T whose ½α_l-coordinate is odd, kΛ_0+μ−nδ is a weight only for n ∈ ½+Z, so
  "c_μ = c^{kΛ_0}_{kΛ_0+μ}" with Kac's definition (n ∈ Z in (0.3)) is identically 0; one must take λ = kΛ_0+μ−½δ.
  The sector formula of Theorem 3 (sum over n ∈ ½Z) is correct as written; only the notation c_μ needs this fix.

### Independent numerics (code/xa_referee/)
Engine: Freudenthal's formula with explicit positive roots (real roots per Kac Prop. 6.3 incl. the three root lengths of
A_{2l}^{(2)}, imaginary multiplicities per Kac Cor. 8.3), W̄-invariance, exact integer arithmetic; ρ̄, h^∨ derived from the
root data (not typed in). Different from the project's Racah recursion. Validated at level 1 (one-dimensional vacuum
spaces, all 11 algebras, `xr_level1_out.txt`); integrality/nonnegativity asserted at every step (wrong root data is
detected immediately, `xr_falsify_out.txt`).

## Table of independent computations (all outputs in `code/xa_referee/*.txt`)

| What | Script / output | Range | Result |
|---|---|---|---|
| Freudenthal engine sanity: level 1 (vacuum spaces 1-dim) for A_1,A_2,D_4,E_6 (untw), C_3,B_3,F_4,G_2 (tw), T_1,T_2,T_3 (A_{2l}^{(2)}); A_1^{(1)} mult(Λ_0−nδ) = p(n) | `xr_level1.py` → `xr_level1_out.txt` | depth 6–40, coords in [−2,2]^l | all pass |
| Sector identity Σ_n mult(kΛ_0+μ−nδ)q^n·∏(1−q^j)^{mult jδ} = q^{(μ|μ)/2k}F_{X,k,[μ]}, every class, 3 representatives each | `xr_sectors.py` → `xr_sectors_main.txt`, `xr_sectors_out.txt` | A_2^{(1)} k=2 (q^30), k=3 (q^22); D_4^{(1)} k=2 (q^16); A_5^{(2)} [C_3] k=2 (q^18); D_4^{(2)} [B_3] k=2 (q^18); D_4^{(3)} [G_2] k=2 (q^26), k=3 (q^20), k=4 (q^8); A_2^{(2)} [T_1] k=2..8 (k=4: q^30, k=6: q^22); A_4^{(2)} [T_2] k=2 (q^20), k=3 (q^7); A_6^{(2)} [T_3] k=2 (q^7); A_3^{(2)} [C_2] k=2 (q^30), k=3; C_4,B_4,F_4 k=2 (q^8–q^12); B_3 k=3 (q^6); A_3,D_5 k=2; E_6^{(1)} k=2 (q^8) | 0 failures |
| Sensitivity: wrong root data / wrong Heisenberg multiplicities / C_3^{(1)} instead of A_5^{(2)} / D = 2 for T_1 | `xr_falsify.py` | k=2,4 | all detected (non-integral Freudenthal or mismatch) |
| Sensitivity: cross-wiring A_5^{(2)} ↔ (B_3,A_1), D_4^{(2)} ↔ (C_3,A_1), C_3 Gram with D=1, A_2^{(2)} k=4 ↔ (A_1,A_3) | `xr_falsify2.py` → `xr_falsify_out.txt` | | fail in 8/8, 8/8, 8/8, 4/4 sectors; control 0/8 |
| Lemma 4 / (3.3) exactly, h^∨(g(X)) = h(X) with h^∨ computed from ρ, ∏(1−q^j)^{mult jδ} = ∏(q^{d_i};q^{d_i})_∞ | `xr_prefactor.py` → `xr_prefactor_out.txt` | A_1..A_7, D_4..D_8, E_6,7,8, B_3..B_6, C_2..C_6, F_4, G_2, T_1..T_6; k = 2..12; products to q^60 | all exact |
| Modularity near q = 1: free fit log F(e^{−t}) = a/t + b_0 + b_1t + b_2t² + b_3t³ | `xr_modular.py` → `xr_modular_out.txt` | (C_3,A_1), (B_3,A_1): t∈[0.10,0.30]; (D_5,A_1), (T_1,A_5): t∈[0.22,0.45] | b_1 = −c_SW/24 to 10 digits, b_2,b_3 ≈ 1e−11; constrained fit with c_SW: rms 2e−15…7e−15, with c_SW ± 1/10: rms ≈ 7e−5 |
| Same, controls | `xr_modular_controls.py` → `xr_modular_controls_out.txt` | rank-1 A=3; B_3 Gram with D(C_3) | power corrections b_2,b_3 ~ 1e−4…3e−2, no c' fits; RR control fits (rms 7e−13) |
| Consistency with `T1-A-even.md` (0.3) and with Theorem 3 total, (T_1,A_r) | `xr_t1even.py` → `xr_t1even_out.txt` | r = 2 (q^14), 4 (q^12), 6 (q^10) | Nahm sum = Σ_N q^{−2N²/(2k′+1)}·(ŝl_2 level-k′ strings) = Σ_m (A_2^{(2)} strings): 0 mismatches |

Observation (not an error): for the twisted pairs the q → 1 growth is governed by c_eff ≠ c:
c_eff(C_3,A_1) = 15/8 (c = 3), c_eff(B_3,A_1) = 3/2 (c = 15/4), while c_eff = c for (D_5,A_1) and (T_1,A_5). So Sun–Wang's
"central charge" c(X,Y) is *not* the effective central charge of the Nahm sum in the non-simply-laced cases (the S-image of a
twisted string function lives on the adjacent algebra g′, whose minimal weight is positive). The modular prefactor
q^{−c/24} is nevertheless exactly right (linear coefficient). This is worth a remark in the note, because the
asymptotic formula of Mizuno/Sun–Wang ties c_eff, not c, to the dilogarithm sum.

## Task 4: the modularity step (KP)

Sources I could obtain: Kac–Wakimoto 1990 (Kac's web page; already quoted in PRIOR-ART); Iohara–Saito arXiv:2502.03316 §1:
"Kac and Peterson … showed … (2) for a twisted affine Lie algebra, a congruence subgroup of SL_2(Z) acts on the C-span of
the normalized characters of each fixed central charge … For BC_l^{(2)} (or A_{2l}^{(2)}), this congruence subgroup is the
so-called θ-subgroup Γ_θ." Kac–Peterson (Adv. Math. 53) and Kac's book: not obtainable (ScienceDirect 403).
No accessible source states modularity of *string functions* of twisted algebras verbatim.

Assessment of the note's argument (§2 (KP)):
1. Γ-invariance of span{χ_Λ} + linear independence of the level-k theta functions gives a finite-dimensional representation
   of Γ on the vector of string functions (times the inverse theta representation). This alone does **not** give modularity
   for a congruence subgroup: a finite-dimensional representation of Γ need not have finite image, let alone congruence
   kernel. The note's sentence "so the vector … transforms under Γ by a finite-dimensional representation" therefore
   proves nothing by itself; the weight of the argument rests on the second half.
2. For r > 1, X_N^{(r)} ≠ A_{2l}^{(2)}: if S maps twisted characters (hence, via thetas, twisted string functions at τ) to
   string functions of the transposed algebra g′ = B_l^{(1)}, C_l^{(1)}, F_4^{(1)}, G_2^{(1)} (possibly with τ rescaled by r),
   then each twisted string function is an untwisted one composed with S (and a rescaling), hence modular for a
   (conjugate of a) principal congruence subgroup by KP Theorem A for untwisted algebras. This is a complete argument
   **given** Kac §13.9 (stated in KW Remark 4.3) and KP Theorem A for the non-simply-laced untwisted algebras; the note
   should say "transposed (untwisted, non-simply-laced) algebra" and mention the rescaling of τ. Acceptable, but the
   attribution "Kac–Peterson Theorem A" for the twisted statement is secondary (KW: "implicitly contained there").
3. For A_{2l}^{(2)} (Theorem 3, all T_l pairs): S preserves the A_{2l}^{(2)} characters, so step 2 is unavailable; one needs
   KP's (or Kac's Thm 13.8's) explicit transformation law *together with* the congruence property (finite image with
   congruence kernel). The note gives no argument for this; it is a **citation gap**: plausible and almost certainly in KP
   (who treat A_{2l}^{(2)} together with untwisted types), but not verified. "Congruence subgroup" for T_l is not justified
   from the texts at hand. A self-contained alternative the note could cite: the A_{2l}^{(2)} level-k modules are ν-twisted
   modules of the rational C_2-cofinite VOA L_k(sl_{2l+1}); twisted modular invariance (Dong–Li–Mason 2000) plus the
   congruence property for twisted modules (Dong–Ren, "Congruence property in orbifold theory", 2018) give a congruence
   subgroup for the characters, and the theta-function argument then passes to string functions. (I have not checked these
   papers' hypotheses in detail.) My numerics for (T_1,A_5) show modularity of weight 0 with the right prefactor (no power
   corrections), but cannot distinguish congruence from non-congruence.

## Task 5: novelty / attribution

- **KNS already state the total formula.** hep-th/9301018 (TeX l.375–425): sector formula ch(Ω^0_λ) = q^{−c_PF/24} Σ_{λ(n)≡λ}
  q^{K(n)}/(q)_n, the identification ch(Ω^Λ_λ) = η^r c^Λ_λ (l.351), and "The character of the space Ω̄^0 is now given by the
  same formula but *without* any restriction on the n-sum". For X simply laced, K(n) is exactly ½ nᵀ(C(X)⊗C(A_{ℓ−1})^{−1})n and
  c_PF = ℓ dim X/(ℓ+h) − r = c(X,A_{ℓ−1}). So for simply laced X, Theorem 1 *including the prefactor* is KNS's conjecture
  (proved by Georgiev for A, BKP in general) + Kac–Peterson; nothing beyond bookkeeping is added. §4's claim that "the
  statement that the unrestricted sums … are therefore modular with the prefactor …" was not found is technically true
  only in the sense that KNS do not write the word "modular"; the note should say plainly that KNS state the unrestricted
  sum as the total parafermion character.
- **Nahm (hep-th/0404120, §"Dynkin diagrams", around eq. (32))** identifies c_eff(Y,A_{k−1}) with the central charge of the
  parafermion coset F̌^k_Y (complement of the Heisenberg algebra in level k), and says the character identification is
  "known or conjectured … probably true in general", citing Kedem et al. and Kuniba–Nakanishi. The note does not cite this;
  it should (it is the physics origin of the "folklore conjecture" for Y = A_{k−1}).
- Is the folklore conjecture for (X simply laced, A_{k−1}) regarded as proved? I found no paper saying so explicitly; Sun–Wang
  treat (A_{n−1},A_{k−1}) as known only under Dunn's gcd condition, which suggests the combination KNS + Georgiev/BKP + KP had
  not been noticed in that community. So "known to experts in a different language, not recorded for these conjectures"
  (§4) is fair for the *conjecture literature*, but the simply-laced case should be presented as an observation, not a
  theorem of the note.
- **Genuinely new content**: the twisted identification (B, C, F_4, G_2 via OT; T_l via Takenaka) with Sun–Wang's D(X)
  convention, Lemma 4 for the twisted algebras, and the resulting cases (C_3,A_1), (B_3,A_1), (T_l,A_{k−1}). Mizuno's prior
  remark covers (C_2,A_1), (G_2,A_1) at k = 2 — correctly credited.
- Mizuno citation details: the remark is in **§4** ("Experiments on modular Nahm sums"), in the case-by-case comments on
  **Table 1 (rank 2)**, not "§5, on his Table 2" (arXiv PDF text, lines around "Then the modularity follows from the result by
  Kac and Peterson [8]"). Also the quoted second matrix [[1,−3/2],[−1,1]] is not symmetrizable with a positive definite AD for
  any positive diagonal D (Mizuno's normalisation A_{21} = (d_2/d_1)A_{12} forces d_2/d_1 = 2/3 and then det(AD) < 0); it is a
  typo in Mizuno (reproduced in the note's "verbatim" quote), presumably for [[1,−3/2],[−1/2,1]] with d = (3,1). Substance
  ((G_2,A_1)) is unaffected; the note should flag it.
- Du–Ji–Shen–Xu arXiv:2602.14483 (indices (2,…,2,1), (1,…,1,2)): their matrices are dense staircase matrices (duals-type), not
  C(B_r)/2 or C(C_r)/2; no overlap with (B_3,A_1), (C_3,A_1).
- Okado–Takenaka's published Correction (Algebr. Represent. Theory 26 (2023) 1691–1692): a search snippet says it concerns
  "Table 2, Twisted affine Dynkin diagrams" (not read; probably cosmetic, would not affect the character formula).
- Georgiev, Gepner, Cherednik–Feigin: not re-read in detail beyond what the note quotes (see "not checked").

## Verdicts

| Statement | Verdict |
|---|---|
| Lemma 4, (0.3)–(0.4), (3.1)–(3.3), table of §0 (h^∨, mult(jδ), Heisenberg product), §5 prefactor for T_l | **No error found** (re-derived by hand; exact machine check for all types, k ≤ 12) |
| Theorem 1 (X simply laced), sector form | **No error found.** Correct and fully supported by KNS + Georgiev/BKP + KP (untwisted). Not new beyond bookkeeping (see novelty). |
| Theorem 2 (X = B, C, F_4, G_2), sector form | **No error found** in the identity; dictionary (Gram, r·ρ_i = d_i, q_OT = q^r, B ↔ D_{l+1}^{(2)}, C ↔ A_{2l−1}^{(2)}) correct and numerically confirmed sector by sector. |
| Theorem 3 (T_l), sector form | **No error found** in the identity (minor notation fix for c_μ with odd ½α_l-coordinate); conditional on Takenaka (unrefereed preprint) as the note says. |
| Corollary 3 (modularity, congruence subgroup) | **Gap (citation)** for the twisted cases: complete only via the S-to-transposed-algebra route, which the note sketches but mis-weights; the finite-dimensional-representation sentence proves nothing by itself. |
| Theorem 3 / Corollary 5, modularity part (A_{2l}^{(2)}) | **Gap**: congruence property of A_{2l}^{(2)} string functions not established from accessible sources; argument in the note does not cover it. Numerically modular (weight 0, right prefactor). |
| Corollary 5 consistency with `T1-A-even.md` | **No error found** (r = 2, 4, 6 exact to q^10–q^14). |
| Remarks (a)–(c) | (a) correct; (b) explicitly unchecked, fine; (c) central charges do add to trD(X)(k−1). |

## Findings, ordered by severity

1. (Gap, §2 (KP), §3.5, §5 end) The modularity of twisted string functions is not proved by the argument given. The
   "finite-dimensional representation of Γ" step does not yield congruence (or even finite-index) invariance. For
   r > 1, ≠ A_{2l}^{(2)}, the S ↦ transposed-algebra step does yield it (given Kac §13.9 and KP for B,C,F,G^{(1)}); the note
   should make that the argument and mention the τ-rescaling. For A_{2l}^{(2)} (all of Theorem 3, Corollary 5, (T_3,A_1)) the
   congruence property is an unverified citation to KP / Kac 13.8; "congruence subgroup" is not justified from the texts
   available. Suggested alternative citation: Dong–Li–Mason twisted modular invariance + Dong–Ren congruence property for
   twisted modules of L_k(sl_{2l+1}).
2. (Attribution/novelty, §4, status paragraph) For X simply laced the full statement — unrestricted sum = total
   parafermion character Σ η^r c_λ with c_PF = c(X,A_{k−1}) — is written in KNS (l.351, l.423–425 of the TeX), and Nahm
   (hep-th/0404120, around (32)) identifies (Y,A_{k−1}) with the parafermion coset. Theorem 1 is therefore KNS's conjecture
   + its proof (Georgiev, BKP) + KP; the note should say this outright and cite Nahm 2004. The (A_{n−1},A_{k−1}) "without gcd"
   and (D_k,A_1) (Wang–Wang's conjecture) consequences are then consequences of 1990s–2021 literature, not of this note.
3. (Citation error, §0 status, §2, §4) Mizuno's remark is in §4, comments on Table 1 (rank 2), not "§5, on his Table 2";
   and the "verbatim" matrix [[1,−3/2],[−1,1]] is not admissible (typo in Mizuno, reproduced without comment).
4. (Overstatement, §0 l.53–54) Exponent claim "−(μ|μ)/(2k) + Z … G is even" false for T_l (½Z); harmless.
5. (Notation, §5 Theorem 3) c_μ := c^{kΛ_0}_{kΛ_0+μ} vanishes identically when μ has odd ½α_l-coordinate under definition
   (0.3) with n ∈ Z; use λ = kΛ_0 + μ − ½δ. The displayed sector identity (n ∈ ½Z) is correct.
6. (Remark worth adding) For twisted pairs c_eff ≠ c (C_3: 15/8 vs 3; B_3: 3/2 vs 15/4); Sun–Wang's c is the modular
   anomaly, not the growth constant. Nothing in the note is wrong, but the phrase "central charge" invites confusion.
7. (Minor) OT never define (q^{ρ_i})_n; the reading (q^{ρ_i};q^{ρ_i})_n follows from OT's (C1). OT's definition of ch has a
   sign slip ("conformal energy −m" ↦ q^m). OT's published Correction not read (snippet: concerns Table 2).
8. (Minor) c_{kΛ_0} = 0 in BKP is not stated by BKP (immediate).

## Not checked / could not check
- Kac–Peterson 1984 and Kac's book (Theorems 13.8–13.10): inaccessible (403). Dong–Ren / Dong–Li–Mason hypotheses not checked.
- The proofs inside BKP, OT, Takenaka (only the statements and conventions were checked); Takenaka's realization details
  (order-4 lift, normalisation of c) were checked only for consistency, and the combined statement numerically.
- Georgiev (5.8)–(5.11), Gepner 1410.1778, Cherednik–Feigin Theorem 1.2/Cor. 1.3 and §2 (remark (c)): not re-read.
- Dunn's thesis (inaccessible); physics papers Genish–Gepner, Nakanishi's dilogarithm/cluster survey, Lee 1109.3667: only
  searched, not read; no statement of the (X,A_{k−1}) modularity found there by search.
- E_7, E_8 and twisted cases with k ≥ 4 beyond G_2, T_1: no numerics (only exact prefactor checks).
- Congruence (vs. non-congruence) of the T_l cases cannot be decided numerically.
