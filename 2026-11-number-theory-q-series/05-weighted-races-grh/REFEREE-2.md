# Referee report 2: Theorem A, its corollaries, Lemma 6.3, Theorem C

Manuscript: `monotonicity/monotone-coupling.tex` (Nov 2026 draft).
Scope: Theorem A (i)-(iii); Cor 6.2, Lemma 6.3, Cor 6.4, Lemma 6.5/6.6,
Cor 6.7, Cor 6.8, Remarks 6.9-6.11; Theorem C and Lemma 7.1;
positioning and bibliography. Theorem B, the coupling and the certificate
are outside this report.

Classification: FATAL / REPAIRABLE (with fix) / MINOR.

**Summary verdict.** I found no FATAL error in this scope. I re-derived
every displayed identity of Section 5, and every sign and constant agrees.
An independent numerical test (my own sieve and script, not the
repository's) confirms the formula of Theorem A(i) to the claimed order,
and all numerical constants quoted in Remark 6.1 and Cor 6.8 reproduce.
The problems are in citations, in the wording of the GRH-status remark,
and in some imprecise hypotheses in the abstract and remarks.

---

## 1. Theorem A(i) / Proposition 5.3: re-derivation

1. **Exact identity.** I checked
   D_m(x) = ∫_{3/2}^x t^{-κ} dD_{m'}(t) = x^{-κ}D_{m'}(x) + κ∫_{3/2}^x t^{-κ-1}D_{m'}(t)dt.
   - The lower limit can be moved to 2 because D_{m'} = 0 on [3/2, 2).
   - The prime 2 and the primes dividing q are handled exactly: χ*(p) = 0
     for p | q, and the jump of D_{m'} at p is χ*(p)p^{m'}, which
     t^{-κ} turns into χ*(p)p^m.
   - Correct.
2. **Step 1 (prime powers).**
   - The identity h(p^k) log p = g(p^k)/k is correct.
   - The square term carries χ(p)^2 = 1 exactly for p ∤ q.
   - The k ≥ 3 bound O(x^{m+1/3} log x + log^2 x) also holds for
     m ∈ (-1/2, -1/3). There the terms with km < -1 are bounded, and
     there are O(log x) values of k.
   - This error is O(x^M/log^2 x) because M > 0 and M > m + 1/3. Correct.
3. **Step 2 (bias).**
   - Partial summation with f(y) = 0 and π(3/2) = ℓ(3/2) = 0 is correct.
   - The bound |f'| ≪ t^{2m-1} holds because (t/y)^{2κ} ≤ 1.
   - The unconditional PNT error O(t/log^2 t) is enough. It produces
     ∫ t^{2m}/log^2 t = ∫ t^{2M-1}/log^2 t ≪ y^{2M}/log^2 y = O(x^M/log^2 x).
   - Main term: (1/2)·y^{2M}/log y·(1/(2M) − 1/(2M')) = κx^M/(2MM' log x).
     This is correct, and the sign is positive.
4. **Steps 3-4.** Two integrations by parts against ψ, then ψ_1.
   - h'(x) = −κx^{m−1}/log x, which gives the sign in eq. (SLpsi1).
   - The bounds |h''| ≪ t^{m−2}/log t and |k| ≪ t^{m−2}/log^2 t are right.
   - R(x) ≪ (1+x^m) log x = o(x^M/log^2 x).
   - Lemma 5.2 (the explicit formula for ψ_1) follows correctly from the
     all-truncations formula, WR Lemma 1. Integrate it over [2, t] at
     fixed T, then let T → ∞. The 2^{ρ+1} constant is bounded and absorbed.
   - The trivial-zero term for even χ (−log x) integrates to O(t log t).
5. **Step 5.**
   - The identity m(m−1) = (m+ρ)(m−1−ρ) + ρ(ρ+1) is correct.
   - Hence κ + m(m−1)/s − m'(m'−1)/s' = κρ(ρ+1)/(ss').
   - The error O(|ρ|^{-2} x^M/log^2 x) is uniform in γ. Lemma 5.1 gives
     an O-constant that depends only on Re s, and |ρ(ρ+1)| ≥ |ρ|^2.
6. **Step 6.**
   - Each central zero contributes −κx^M/(MM' log x) to S_Λ, so it adds
     +n_0/(MM') to the bracket.
   - Total bracket: (1+2n_0)/(2MM') + 2Re Σ_{γ>0} x^{iγ}/((M+iγ)(M'+iγ)).
   - The sign of the oscillating sum is +, and the sign of the central-zero
     term is +. Both are correct.
   - The error is uniform in x with a constant depending on (χ, m, m').
     The constant degrades as M → 0, which the statement allows.
7. **Numerical test (independent).**
   - Setup: sieve to 10^7, q = 4, the 511 zeros in
     `prototype/explore/chi4_zeros.txt`, 161 log-spaced points in
     [10^3, 10^7]. I computed I(x) directly in extended precision.
   - Residual: res(x) = I(x) log x/(κx^M) − bracket(x).
   - Results for the size of res·log x:

     | (m, m') | max \|res·log x\| | behaviour for x ≥ 10^5 |
     |---|---|---|
     | (0, 1) | 2.50 | essentially constant, ≈ 2.46 |
     | (1, 2) | 0.25 | |
     | (0, 7) | 0.40 | |

   - In every case the quantity is bounded with no drift in log x. A wrong
     sign of the bias or of the n_0 term would make it grow like log x.
   - For (0, 1) the constant ≈ 2.46 is the expected second-order term.
     For example, the bias alone contributes (1/(2κ))(M^{-2} − M'^{-2}) ≈ 1.78.
   - My value 0.902 at x = 10^3 agrees with the paper's table in §9(2).

**Verdict:** Theorem A(i) is correct as stated, including uniformity.

## 2. Theorem A(ii), (iii); Cor 6.2

- (ii) is correct.
  - The bound |2Re Σ| ≤ B/2 gives bracket ≥ λ/2.
  - The threshold log x ≥ 4C_0/(κλ) is correct.
  - The two implications are correct; ties are handled correctly
    (D_{m'} ≥ 0 ⇒ D_m > 0, and D_m ≤ 0 ⇒ D_{m'} < 0).
- (iii) is correct.
  - Cauchy–Schwarz: Σ 4/(M^2+γ^2) = σ_M^2/(2M^2).
  - σ_M is strictly increasing in M; this is trivial termwise and needs
    no NCZ.
  - Hence σ_{M'}^2 ≤ 2 ⇒ σ_Mσ_{M'} < 2 ≤ 2 + 4n_0.
- Cor 6.2 is correct for all four densities, natural and logarithmic.
  - Every density is monotone under inclusion and ignores bounded sets.
  - Pairwise positivity holds for **all** pairs M < M' ≤ M_σ, because only
    σ_{M'}^2 ≤ 2 is needed.
  - Strictness and existence are correctly disclaimed.

## 3. Lemma 6.3 and Cor 6.4

### ANS14 statements, checked against arXiv:1306.1657

- **(1.6)-(1.7).** φ(y) = c + Re Σ_{λ_n ≤ X} r_n e^{iλ_n y} + E(y, X) for
  y ≥ y_0 and all X ≥ X_0. Here λ_n is a non-decreasing sequence of
  positive numbers, multiplicity allowed. The error satisfies
  (1/Y)∫|E(y, e^Y)|^2 → 0.
- **Cor 1.3(a).** If r_n ≪ λ_n^{-β} with β > 1/2 and
  Σ_{T<λ_n≤T+1} 1 ≪ log T, then φ is B^2-almost periodic and has a
  limiting distribution.
- **Thm 1.4.** This is the vector-valued version.
- **Thm 1.9.** If {λ_m} is linearly independent over Q, then the Fourier
  transform of the vector limit is exp(−iΣc_kξ_k)·Π_m J_0(Σ_k r_k(λ_m)ξ_k).
  For complex r the argument has to be read as the modulus.

**Verdict on Lemma 6.3:**
- The hypotheses are verified: r_γ = O(1/γ), the window count, and the
  error term from Lemma 2.2 at T = e^Y.
- The Cramér–Wold / Lévy argument is correct.
- The identification with |ξ_1 r_γ(M) + ξ_2 r_γ(M')| is the right reading
  of Thm 1.9.
- (MINOR) ANS Thm 1.4 and Thm 1.9 are already stated for ℝ^ℓ-valued φ, so
  the Cramér–Wold step can be replaced by a direct citation.
- (MINOR) Gaps item 2 says the ANS statements were taken from WR. I have
  now checked them against the source, and they match.

### Cor 6.4(a)

- The Kronecker argument is correct.
- The inhomogeneous Kronecker theorem needs γ_1, …, γ_J linearly
  independent over Q, which LI gives.
- The arithmetic −cos ε Σ_{j≤J} b_j + tail ≤ −B + ε(B+2) is correct.
- The "with equality" clause needs LI. The inequality needs GRH only.

### Cor 6.4(b)

- The boundary of {y ≤ 0 < y'} is null because both marginals are
  absolutely continuous under LI + NCZ (WR Thm 1).
- Null violation set for λ ≥ 0: the a.s. bound Y_M/M − Y_{M'}/M' ≥ κλ
  handles λ = 0 too.
- Positive density for λ < 0: I re-checked each step.
  - (Ydiff) holds, and on A_τ the difference is ≤ κλ/2.
  - Y_{M'} ∈ (0, 3ε_1) on A_τ ∩ B_1.
  - 3ε_1/M' ≤ 3κ|λ|/8 < κ|λ|/2.
  - W' is independent of A_τ and has full support by Lemma 2.6. The tail
    of c_γ(M') still has Σ c = ∞, Σ c^2 < ∞ and c → 0.

**Verdict:** Cor 6.4(b) is correct.

## 4. Lemma 6.6, Cor 6.7 (Landau), Cor 6.8

### Lemma 6.6

- The Mellin transform is correct. The inner integral gives
  p^{-s}(1/s − 1/(s+κ)) = κp^{-s}/(s(s+κ)).
- Fubini is justified for σ > m+1.
- Sign: F = −κP(s−m)/(s(s+κ)).
- The decomposition (Pstructure) and H_χ are correct.
  - The k ≥ 2 part of log ζ(2w) needs Re w > 1/4.
  - The k ≥ 3 part needs Re w > 1/3.

### The singularity at w = 1/2

- −(1/2) log ζ(2w) = (1/2) log(2w−1) + O(1).
- So P(σ) → −∞ as σ → 1/2+. Equivalently Σχ*(p)p^{-σ} → +∞, which is
  the positive bias, consistent with the bias term of Prop 5.3.
- Landau is used only on the real axis, so this logarithmic singularity is
  harmless. The argument never needs F to be singular at m + 1/2. It only
  needs "σ_c > m + β_0 ⇒ contradiction".

### Landau's lemma

- I checked the statement against MV07 Lemma 15.1, in the shifted form
  used in the paper.
- Hypotheses: A bounded and Riemann integrable on finite intervals, and
  one-signed for x > X_0. Both hold for I, and also for −I.
- σ_c ≤ m+1 because |I(x)| ≤ Σ_{p≤x} p^m ≪ x^{m+1}.

### Cor 6.7

- The proof is correct.
  - The tube U is simply connected, and U ∩ {Re w > c} is connected.
  - The continuation across σ_c is legitimate.
  - The identity e^G = L on Re w > β_0 gives zero-freeness.
- Real zeros in (1/2, 1) are correctly identified as the obstruction
  (Remark 6.9).

### Cor 6.8(a)

- There are no real zeros in (0, 1) for χ_3 and χ_4. The proof groups
  consecutive terms of the convergent alternating Dirichlet series:
  (4k+1)^{-σ} − (4k+3)^{-σ} > 0, and (3k+1)^{-σ} − (3k+2)^{-σ} > 0.
- This is elementary and correct. So β_0 = 1/2, and Cor 6.8(a) is fully
  proved, unconditionally.

### Cor 6.8(c)

- This part is stated correctly, with the unconditional quantity
  V_χ(M') = 4M'(Λ'/Λ)(M'+1/2).
- Under GRH and NCZ, V = σ^2. In general V = σ^2 + 4n_0, which does not
  matter here.
- I reproduced V_{χ4}(3/2) = 1.3685551 and V_{χ3}(3/2) = 1.0070099. The
  latter is also enclosed in `cert/aux_results.txt`.
- The identity I_{0,1}(x) = x^{-1}∫_2^x D_0 is correct.

## 5. Theorem C and Lemma 7.1

### Lemma 7.1

- The Fejér kernel argument is correct.
  - ∫ F(v) e^{iξv} dv = 2π(1−|ξ|)_+.
  - F ≥ 0.919 on |v| ≤ 1; only |v| ≤ 1/2 is needed.
  - The diagonal-band bound uses #{γ' : |γ'−γ| < 1} ≪ log(γ+3).
- The constant is uniform in U and T.
- The Fejér step itself does not use NCZ. NCZ enters only through the
  quoted expansion, WR Lemma 2.2: there |s| ≍ |γ| uses γ ≥ γ_1 > 0.
- The quoted pointwise error bound at T = e^{U+1} is O(1). I checked it
  against the companion paper's display (WR §3, the ANS verification):
  that bound holds for every T pointwise, not only in mean square.

### Theorem C(i)

- Partial summation is correct for every real b, including b < 0.
  - The lower limit 3/2 is fine because (log t)^{b−1} is smooth on
    [3/2, ∞).
  - R_m = 0 on [3/2, 2).
- The decomposition (Cdecomp) is correct.
- **The claim is genuinely pointwise.** The kernel e^{−M(u−v)} v^{b−2} is
  integrated against |E|. On [u/2, u] the integral is bounded by
  Σ_k e^{−Mk} (∫_{unit interval} |E|^2)^{1/2} ≤ K_m^{1/2}/(1−e^{−M}), using
  only the *uniform local* L^2 bound. On [log 2, u/2] the factor
  e^{−Mu/2} kills everything.
- So |R_b(e^u)| ≪ 1/u for every u ≥ 2, with no averaging. The mean-square
  control enters only through the sup over U in Lemma 7.1, which is a
  pointwise-in-U statement.
- Independent check: q = 4, m ∈ {0, 1}, b ∈ {−1, 0, 2, 3}, 41 points in
  [10^3, 10^7]. (E_{m,b} − E_{m,1}) log x stays in bounded ranges that
  agree with §9(4). For example, m = 0, b = 0 gives [2.79, 3.15].

### Theorem C(ii), (iii)

- Correct by Lemma 2.5(a), (b) and WR Thm 1.
- The densities are logarithmic.

## 6. Positioning and bibliography

### REPAIRABLE

**R1. Fujii citation is wrong.**
- The manuscript has "A. Fujii, Proc. Japan Acad. (1991)" (`\bibitem{Fuj91}`,
  cited at lines 374 and 1789).
- Platt–Trudgian's reference list, which is the manuscript's stated source,
  gives: A. Fujii, *Some generalizations of Chebyshev's conjecture*,
  Proc. Japan Acad. Ser. A Math. Sci. **64** (1988), no. 7, 260–263.
  Fujii's Theorem 1 is the α < 4.19 statement.
- Fix: correct the entry and rename the key to Fuj88. The same error is in
  `PRIOR-ART.md`.
- In the same list, Landau's paper is in two parts: *Über einige ältere
  Vermutungen und Behauptungen in der Primzahltheorie I & II*,
  Math. Z. **1** (1918), 1–24 and 213–219.
- Hardy–Littlewood is Acta Math. **41**, 119–196. PT gives the year 1917;
  it is often cited as 1916.

**R2. Remark 6.11 (status of GRH) needs more cautious and more accurate wording.**
- "has recently been announced and machine-checked in Lean" asserts the
  machine-checking as a fact. In the next clause the authors say they
  examined neither the argument nor the formalization.
- The bibliography marks both sources as second-hand. However,
  `monotonicity/GRH-PROOF-DIGEST.md` records that the Alpöge note and the
  Lean FRO statement were read in primary form, so the manuscript and the
  repository disagree. An unread PhilArchive paper claiming a vacuity
  defect is also logged there.
- Suggested wording: "A proof of GRH for Dirichlet L-functions was announced
  on 20 November 2026 [Alp26], accompanied by a Lean formalization that
  [Lean26] reports to have checked; no human-readable proof has been
  refereed. We have not examined either, and do not use them."
- The consequences listed need two corrections:
  - Cor 6.4(a): its equality clause needs LI. Only the inequality needs
    GRH alone (the remark says "Cor 6.4(a) needs GRH only").
  - "Theorems B, C(iii) ... would still require LI" should read "would
    still require LI, and (for χ other than χ_3, χ_4) NCZ".

**R3. The abstract's equivalence uses a hypothesis defined only under GRH.**
- The abstract says: "For q ∈ {3,4} and pairs with σ_{M'}^2 ≤ 2, the
  eventual positivity of this mean is equivalent to GRH(χ)".
- σ_{M'} is defined as a sum over γ > 0 under GRH, so it cannot be the
  hypothesis of an equivalence with GRH.
- Fix: state the hypothesis as V_χ(M') ≤ 2, as in Cor 6.8(c), or say
  "pairs with 4M'(Λ'/Λ)(M'+1/2, χ) ≤ 2".

**R4. Remark 6.10 overgeneralises the cited literature.**
- The remark says: "In all of these ... the result is established only for
  a limited range of the smoothing parameter (α, k)".
- This is true for Fujii, Platt–Trudgian and Washington.
- It is not true for Johnston (no parameter) or Zhao (no parameter; this
  is about Mertens error terms).
- Suzuki (arXiv:2411.07436, abstract) proves one direction, sign constancy
  ⇒ RH/GRH, plus sign constancy *on average*. That is not an equivalence
  of the "positive bias beats an absolutely convergent zero sum" type.
- Fix: describe Suzuki and Zhao accurately and restrict the "limited range"
  sentence to the parametric results.

### Verification status of every bibliography entry

**Verified online in this review:**
- ANS14: Cor 1.3, Thm 1.4 and Thm 1.9 read in the arXiv text. The QJM
  volume and pages were not checked.
- PT19: arXiv:1802.00089 read; α < 20.40442 and Fujii's α < 4.19 confirmed.
  The IJNT pages were not checked.
- Was24: arXiv:2209.12845. Theorem 3 is confirmed: for 0 < k ≤ 10.32,
  RH ⇔ ∫_1^x (π_k(t) − π(t^{k+1})) dt < 0 for all large x. The arXiv page
  lists no journal reference, so "Ramanujan J. 65 (2024) 783–795" is
  unverified.
- Joh23: arXiv:2201.06184. The statement is confirmed (RH ⇔ ∫_2^x (π − li) < 0
  for all x > 2). The arXiv page says only "to appear in Canad. Math. Bull.",
  so vol. 66 (2023) 185–195 is unverified.
- Zhao25: confirmed. Res. Number Theory 11 (2025) Paper 62; quadratic-character
  analogues mentioned.
- Suz25: confirmed. Ramanujan J. 68 (2025) art. 95.
- Hay25: arXiv:2512.23302 confirmed (v1 Dec 2025, v3 Jun 2026). Title and
  description (natural density 1, weight p^{-1/2}, GRH only) are accurate.
- BBLMR24: confirmed (PUMP J. 7 (2024) 29–50). The description is accurate:
  each class is compared with π(x^{k+1}; m, n), not class against class.

**Unverified or marked TODO** (must be checked before submission):
- HL18 (TODO)
- Lan18 (TODO; it has two parts, see R1)
- Fuj91 (TODO and wrong, see R1)
- LWS26 (TODO; abstract only, authors' initials missing)
- Alp26 (TODO)
- Lean26 (TODO)
- GW18
- AK23, Sheth24, SK25, Aym22, BMS26: not checked here
- ANS14, PT19, Was24, Joh23: journal volume and pages
- MV07 Lemma 15.1: content agrees with the standard statement, but the page
  was not checked
- RS94, Dav00, HW, Fel71, Bil99: standard

### Novelty claims

- They are suitably hedged.
- (MINOR) The identity in Prop 5.3(a) is one line of partial summation, and
  I_{m,m'} is literally a Riesz typical mean of order one. The claimed
  novelty is the *interpretation* as a link between two races, plus Cor 6.2.
  The introduction already says roughly this, but the sentence "What
  appears to be new is the observation ..." could be toned down further.

## 7. MINOR items

- M1. The abstract has "for instance when σ_Mσ_{M'} < 2". This is fine
  under GRH; note that it is 2 + 4n_0 in general.
- M2. Gaps item 2 says the ANS statements were taken only as reproduced in
  WR. They now match the source; update the item.
- M3. Lemma 2.2 widens WR's e^{u/2} to u·e^{u/2} and e^{−c√u} to u·e^{−c√u}.
  This is harmless, and the text says so.
- M4. Prop 2.4 cites "[WR, Proposition 2, eq. (8), and §2]". WR's
  Proposition 2 itself is stated for χ_4 only; the general-χ identity is
  eq. (8) in the text after it. The citation is acceptable, but say this
  explicitly.
- M5. Remark 7.2(3) (central zeros in Theorem C) is labelled as an
  expectation; the 1+2n_0 constant it predicts is correct.
- M6. Theorem A(ii): x_0 is effective in principle only. The paper says so
  (Gaps 5).

## 8. Certification

### Certified as correctly proved as written

All of the following are under the stated hypotheses, with the quoted
inputs (WR Lemma 1, WR Lemma 2.2 / Thm 1, ANS14, MV07 Lemma 15.1,
Kronecker) granted:

- Lemma 5.1, Lemma 5.2, Proposition 5.3(a), (b) = Theorem A(i), including
  uniformity and all signs and constants (n_0 term included).
- Theorem A(ii), Theorem A(iii).
- Remark 6.1(1)-(3); I reproduced the numerical values.
- Corollary 6.2.
- Lemma 6.3.
- Corollary 6.4(a), (b).
- Lemma 6.5 (quoted) and Lemma 6.6.
- Corollary 6.7.
- Corollary 6.8(a), (b), (c), and the displayed (0,1) equivalences for
  χ_4 and χ_3.
- Lemma 7.1, Theorem C(i) (genuinely pointwise), C(ii), C(iii).

### Not certified as written

- Abstract sentence on the GRH equivalence (R3). Corrected statement: "For
  q ∈ {3,4} and pairs with 4M'(Λ'/Λ)(M'+½, χ) ≤ 2, the eventual positivity
  of this mean is equivalent to GRH(χ)."
- Remark 6.11 (R2). Corrected statement: "… was announced [Alp26],
  accompanied by a Lean formalization reported in [Lean26] to have been
  checked; not refereed, not examined by us, not used. If confirmed,
  Theorem A, Cor 6.2 and the inequality in Cor 6.4(a) become unconditional,
  and Cor 6.8 becomes an unconditional positivity statement. Theorems B and
  C and the Main Theorem would still need NCZ (automatic for q = 3,4) and,
  for statements about races, LI."
- Remark 6.10's sentence "In all of these … limited range" (R4).
- Bibliography entry Fuj91 (R1), and the TODO entries listed above.
