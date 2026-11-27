# Referee report: `proofs/iv-all-k.md` (KRTW identity (iv) for all k; Theorems 4, 5, 6)

Referee 2. This is an adversarial read with independent computations. All scripts are in `referee2/`. They are my own
code: nothing is imported from `code/`, `referee/` or `rr-identities/`. The only libraries used are python-flint
(integer polynomials), sympy and mpmath. Date: 2026-11-27. Status: **DONE.**

## Overall verdict

**None: no fatal error and no gap.** The proof is complete and correct for every k ≥ 1. Its inputs are:

- `F-all-k.md` §§1–7, which was refereed separately;
- (B2) of `KRTW-F-moves.md` Theorem B;
- the (B1)-pair of `KRTW-F-moves.md` §4 for general a;
- standard facts: Bailey's lemma (L∞), Bailey inversion, q-Chu–Vandermonde, Cauchy and the Jacobi triple product.

Hence KRTW's fourth conjectural Andrews–Gordon-type identity (iv), as printed in arXiv:2404.03851v2, holds for all k ≥ 1.
Theorem 4, (xq)_∞K^Ω_k = J_k·Σλ_n^{k+1}b^Ω_n, and Theorem 6, K^Ω_k = B^{(k)}_{k+1}J_k, also hold for all k.

Every step was re-derived by hand. Every intermediate fact that iv-all-k extracts from a *proof* in F-all-k was checked
against the source, and independently in position space. The problems listed below are exposition only (E1–E8);
none affects correctness.

## Statement check (§0 against LITERATURE §1b, KRTW TeX lines 2435–2501)

Letter by letter, the §0 statement matches the print:

- the summation over all r_1,…,r_k, s_1,…,s_k ≥ 0;
- r_{k+1} = s_0 := 0;
- the summand q^{(r_i+s_i)²+s_i²}/((q;q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}});
- Ω = Σ_{i=1}^{k−1} q^{r_i+2s_i}(1−q^{2s_{i+1}−2s_i}) + q^{r_k+2s_k};
- the product (q^k,q^{k+1},q^{k+3},q^{k+4},q^{2k+4},q^{2k+4};q^{2k+4})_∞;
- the denominator (q²;q²)_∞(q;q)_∞.

The convention 1/(q;q)_m = 0 for m < 0 is the one `REFEREE-F-all-k.md` found in KRTW's TeX:
(a;q)_n := (a;q)_∞/(aq^n;q)_∞. I did not re-read the TeX for this.

Restricting to r decreasing and s increasing loses nothing. The non-monotone terms carry 1/(q)_{<0} or 1/(Q;Q)_{<0},
and multiplying the Laurent polynomial (1 − q^{2s_{i+1}−2s_i}) by 0 gives 0. So the left side of (iv) is K^Ω_k(1).

The extra grading x^{r_i+2s_i} differs from KRTW's z^{r_i+s_i}, but it is harmless at x = 1. The definitions of
B^{(k)}_{k+1} and J_k agree with `KRTW-general.md` §8. The factored form of b^Ω_n is correct (by hand, and S4 below).

## Section-by-section verdicts

| § | item | verdict | how checked |
|---|------|---------|-------------|
| 0 | statement; K^Ω_k; b^Ω_n; factored form, b^Ω_0 = 1 | correct | by hand; S4; computation 1 |
| 1 | Lemma 1.1 | correct; see note 1.1 | by hand; `L11op`; DP = naive enumeration |
| 2 | Lemma 2.1 | correct; see note 2.1 | by hand; S1; computation 3 |
| 2 | Lemma 2.2 | correct; see note 2.2 | by hand; computation 3 |
| 2 | Prop 2.3 | correct; see note 2.3 | by hand; `P23`, brute force for k ≤ 5 |
| 3 | Lemma 3.1 | correct; see note 3.1 | by hand; computation 4 |
| 3 | Theorem 3.2 | correct; see note 3.2 | by hand; `T32a`, `T32b`, `T32c` |
| 4 | Lemma 4.1 | correct for all N, M; see note 4.1 | by hand; S2, S3 |
| 4 | proof of Theorem 4 | correct; see note 4.2 and E1 | by hand; `TH4`, `L41s`, `B2` |
| 5 | Theorem 5 | correct; see note 5 | by hand; S4; `r2_x1.py`; computation 1 |
| 6 | Theorem 6 | correct; see note 6 | by hand; S5, S6, S7; `TH6a`, `TH6b` |

### Notes to the table

**1.1 Lemma 1.1.**
- The bond between levels i and i+1 is 1/((q)_{r_i−r_{i+1}}(Q;Q)_{s_{i+1}−s_i}), where (q)_{r_i−r_{i+1}} comes from index i and
  (Q;Q)_{s_{i+1}−s_i} from index i+1. In T, (R,S) = level i and (r,s) = level i+1.
- The i-th Ω-term (i < k) multiplies level i by q^{r_i+2s_i} and turns the bond into 1/(Q;Q)_{s_{i+1}−s_i−1}. This is
  exactly OTF_{i+1}, since the s = S term drops out.
- The i = k term is q^{r_k+2s_k}·u = Ou, because u does not depend on S. For k = 1, K^Ω_1 = ω(Ou).
- Splitting Ω into its 2k−1 monomials is legitimate because the x-degree bounds the indices (E6).

**2.1 Lemma 2.1.** I redid the normal-ordered coefficients by hand:
- WhW⁻¹ = P⁻¹ − xq⁻¹PY𝖲 + xPY, from (1.3);
- the coefficients of 𝖲^j𝔲^b, times (q)_j(Q;Q)_b, are q^jP⁻¹ + xPY in hE, q^jP⁻¹ + xQ^bPY in Eh̃, and xPY(1−Q^b) in xPY𝔲E.

Restriction to N² is harmless. For F supported in N², hF is supported in N². On N², h reads only (R,S) and (R−1,S),
and EWF = 0 for R < 0. The term xPY𝔲 reads EWF at (R,S+1) ∈ N² (E7).

**2.2 Lemma 2.2.** The summation by parts uses w(r+1,s)q^{−r−1} = xq^{r+2s}w(r,s) and F(−1,s) = 0. Each sum converges
x-adically for F ∈ 𝒱.

**2.3 Proposition 2.3.**
- hu = u + xOu holds, including R = 0, where 1/(q)_{−1} = 0.
- The telescoping Σ_{i=1}^{k−1}(a_i − a_{i−1}) = a_{k−1} − a_0 is right.
- a_{k−1} = K_k + x·(k-th term), and a_0 = K^♭_k by Lemma 2.2.
- The k = 1 case is right.
- I also checked k = 1 by hand as a termwise identity: x⁻¹(K^♭_1 − K_1) = K_1(xq).

**3.1 Lemma 3.1.**
- (a) h does not change S. On the slice S = σ it is h_{xq^{2σ}}, and (B-iii) holds, including at R ≤ e. So h is
  diagonal on Ξ_{e,σ}.
- (b) The identity is linear. It is valid for every finite combination of Ξ's, in particular for hφ. Here
  (TΞ_{e,s})(·,0) = w(e,s)φ_e(·;xq^{2s})/(Q;Q)_s is a finite combination of φ_e(·;x) by (B-iv), so v_x is defined (E3).
- (c) T(hφ) = λ(hφ − xOφ), and its slice S = 0 is λ(h_xF_0 − xq^RF_1). The two facts are correctly taken from the proof
  of F-all-k Prop 4.1:
  - F_1(R) = (xq)_2G(R+1), with G built from φ_n(·;x), n ≥ 1, gives v_x(q^RF_1) = 0 by (B-v);
  - the slice-0 equation of H′φ = εφ gives (1+x)μ_1 = εμ_0.
- The M = 0 case is right: F_1 = 0, φ_{L,0} = Ξ_{L,0}, and ε_{L,0} = (1+x)(q^{−L}+xq^L).

**3.2 Theorem 3.2.**
- The summability family is the one listed in the proof of F-all-k Thm 6.4: ord ≥ r + 2M_b + (k−1)(N_b+M_b), with
  r ≥ L_b. The factor q^{−r} leaves ords unchanged, so ω∘P⁻¹ can be applied termwise.
- The bookkeeping ũλW/(xq)_∞ = λt/((xq)_∞(x²Q;Q)_∞) is right, since t contains λ.

**4.1 Lemma 4.1.** Proof for all N, M:
- β^B_n/β^J_n = (1+x)q^n/(1+xq^{2n}), because the factors (x²;Q)_n/(Q;Q)_n cancel;
- ε − μ_N − 2x = (u−v)(1+xu²)(1−xuv)/(u²v), with u = q^N, v = q^M;
- the claim is then equivalent to a rational identity in u, v, x, which holds.

**4.2 Proof of Theorem 4.**
- The symmetrisation Σ_{N≥M}(f(N,M)+f(M,N))/(1+δ) = Σ_{n,m} f(n,m) is right.
- Convergence: ord λ_n^{k+1}b^Ω_n = (k+1)n−1 and ord λ_m^{k+1}β^J_m = (k+1)m.
- Passage from 𝔽 to Z[x][[q]]: every x^D-coefficient is a finite sum of rational functions regular at q = 0, and all four
  series lie in Z[x][[q]] (J_k = C_k/C_1 lies there because C_1 ≡ 1 mod q is invertible). The coefficientwise maps into
  Q[[x,q]], from the x-series with such coefficients and from Z[x][[q]], are injective ring homomorphisms.

**5 Theorem 5.**
- Each series is a q-adically summable family in Z[x][[q]], so x ↦ 1 commutes with the sums. This needs the factored
  form of b^Ω_n, because the definition contains 1/(1−x).
- b^Ω_n(1) = (−1)^n(q^{n²−2n}+q^{n²+2n}).
- Triple product with A = k+2, B = −2.
- J_k(1) from (B2).
- The product equals the right side of (iv).

**6 Theorem 6.**
- Step 1: the (B1)-pair holds for general a (re-read: a is generic throughout that proof).
- Step 2: (L∞) twice relative to xq.
- Step 3: Γ^{(1)}_{d,d} = (1−aq)/(1−aq^{2d+1}) and Γ^{(1)}_{d,d+1} = −aq^{2d+1}(1−aq)/(1−aq^{2d+1}), re-derived by partial fractions.
- Re-derived α*_e = (−1)^e x^{2e−1}q^{3e²−2e}(x²Q;Q)_{e−1}(1−w+xw²−x³w³)/(Q;Q)_e, and the factorisation
  1−w+xw²−x³w³ = (1−xw)(1+(x−1)w+x²w²). So α*_e = λ_e²b^Ω_e.
- Step 4: k−1 (L∞) steps relative to x. For k = 1 there is no step.
- Step 5: the limit n → ∞ holds coefficientwise in x (finitely many tuples and r per x^D, with q-adic convergence).

## Problems found

All are exposition-level or cosmetic. There is no fatal problem and no gap. Line numbers refer to `proofs/iv-all-k.md`
as of 07:18.

1. **E1 (exposition; Theorem 4 proof, line 166).** "Σ_mλ_m^{k+1}β^J_m = (x²Q;Q)_∞J_k(x) by (B2) … This is Theorem 4 in 𝔽."
   - (B2) is proved in 𝓡 = Z[[q]][[x]], not in 𝔽 = Q(q)((x)).
   - The transfer is true but is not stated. It holds because both sides of (B2) have x^D-coefficients that are finite
     sums of rational functions regular at q = 0, and the coefficientwise Taylor map from such series to Q[[q]][[x]] is
     injective.
   - *Fix:* one sentence. Alternatively, pass Theorem 3.2 together with Lemma 4.1 to Q[[q]][[x]] first, then use (B2)
     there. The same tacit step occurs in F-all-k's proof of Theorem 1.
2. **E2 (cosmetic; §0 "What is new", lines 46–47).** "for k = 1 … so Theorem 6 is the termwise identity (F_1) at xq".
   - Only the consequence K^Ω_1 = B^{(1)}_2J_1 is termwise.
   - The first assertion of Theorem 6 at k = 1, (xq)_∞B_2(xq) = Σλ_n²b^Ω_n, is not termwise. It is (B1) at k = 1 with
     x ↦ xq, rewritten by the change of a: (xq²;q)_∞B_2(xq) = Σα″_n, and Σα*_e = (1−xq)Σα″_d because Σ_eΓ^{(1)}_{d,e} = 1−xq.
   - *Fix:* "for k = 1, Theorems 4 and 6 follow from (F_1) and from (B1) at x ↦ xq".
3. **E3 (exposition; Lemma 3.1(b), lines 119–121).** Two facts are needed and only implicit:
   - (TF)(·,0) is a finite combination of the φ_e(·;x), by (B-i) and (B-iv);
   - the φ_e(·;x) are linearly independent (triangular: φ_e(R) = 0 for R < e and ≠ 0 at R = e).
   These make v_x, and comparisons of coefficient sums, well defined. "Uses only linearity" covers the rest.
   *Fix:* half a sentence.
4. **E4 (cosmetic; §6 item 3, line 199).** The auxiliary symbol w := q^{2e} clashes with the weight w(r,s) of §1.
   *Fix:* rename it, e.g. to y.
5. **E5 (cosmetic; §6 item 4, line 202).** "(for k = 1 there is no step and N_0 := n in item 2)" is hard to parse; nothing
   called N_0 occurs in item 2. *Fix:* "for k = 1 there is no step, and β^{(1)} = β″ of item 2 (whose sum is over
   n ≥ N_1 ≥ N_2)".
6. **E6 (minor; Lemma 1.1, lines 71–72).** "All sums converge x-adically, so the k pieces can be summed separately" is
   stated loosely. What is used is that only finitely many index tuples contribute to each power of x. That is what
   makes splitting Ω into its 2k−1 monomials, and regrouping into nested sums, legitimate. *Fix:* say so.
7. **E7 (minor; Lemma 2.1, lines 80–82).** The restriction to N² is delegated to the F-all-k §2 paragraph. That
   paragraph does mention "the operators of iv-all-k §2", so nothing is missing. For a self-contained reading, add:
   - h maps functions on N² to functions on N², and reads only (R,S) and (R−1,S);
   - xPY𝔲 reads EWF only at (R,S+1) ∈ N².
8. **E8 (note; lines 44–45).** The k = 2 attribution (KRTW Thm 4.10 together with Bartholdi 2611.04020) is not
   re-derived, as the file itself says. The project record (`rr-identities/PRIOR-ART.md`) names Bartholdi's Theorem A.
   I did not check this either. The proof covers k = 2 anyway, so this matters only for the novelty statement. Also,
   the status line (line 3, "adversarial read pending") should be updated.

## Independent computations

All of these are in `referee2/`. They use exact integer or rational arithmetic unless marked "num".

1. **`r2_iv_series.py 400 14` → `r2_iv_series_q400_k14.txt`.** Identity (iv) itself, **k = 1…14, mod q^401**.
   - The multisum is evaluated literally, level by level. Each Ω-term is inserted at its own bond as
     q^{r_i+2s_i}(1−q^{2(s_{i+1}−s_i)}), and the last term as q^{r_k+2s_k} at level k.
   - The level-by-level evaluation equals naive enumeration of all tuples for k ≤ 3, mod q^61.
   - The check is LHS·(q;q)_∞(q²;q²)_∞ = the numerator product. All pass.
   - Positive controls: (ii) and (iii) for k ≤ 4 pass.
   - Negative controls (all fail, as they must): LHS(k) against RHS(k+1); last Ω-term q^{r_k+s_k}; Ω without its
     (1−q^{…}) factors; (1−q^{s_{i+1}−s_i}) in place of (1−q^{2s_{i+1}−2s_i}).
2. **`r2_pointmass.py` → `r2_pointmass.txt`.** Lemma 2.1 (Th − hT = −xOT), Lemma 2.2 and hu = u + xOu.
   - Operators implemented from their definitions, exact Fractions.
   - Points (q,x) = (2/7, 3/11), (−3/5, 5/13), (1/3, −2/9).
   - Inputs: all point masses in [0,5]×[0,4] plus 6 random combinations.
   - Result: 3 × 1176 comparisons, 0 failures. The commutator Th − hT is nonzero at 354 of the 1080
     (test function, box point) pairs, so the test is not vacuous.
   - Negative control (O without its S-shift) detected on all 36 test functions.
3. **`r2_bivariate.py`** → `r2_bivariate_q41_k6.txt` (L = 41, k ≤ 6, brute force k ≤ 5), `r2_bivariate_q81_k8.txt`
   (L = 81, k ≤ 8, brute force k ≤ 5) and `r2_bivariate_q61_k12.txt` (L = 61, k ≤ 12, brute force k ≤ 3).
   - Arithmetic is in Z[[x,q]] mod (x^L, q^L). Every final series here has q-exponent ≥ x-exponent termwise, so the check
     covers all powers of x mod q^L. The exceptions are K^♭ and the first formula of Theorem 3.2, where
     x-exponent ≤ q-exponent + 1; those checks are complete mod q^{L−1}.
   - Built from definitions only: K_k, K^Ω_k (literal Ω), K^♭_k, C_k, J_k = C_k/C_1, B_{k+1}, B^{(k)}_{k+1}; b^Ω_n from its
     definition β^B_n(q^{−2n}+x−1+x²q^{2n})/(x(1+x)); t, ε from F-all-k §7 and Lemma 3.2.
   - Checks, all OK (87, 109 and 145 OK lines in the three runs):
     - `K=bf`, `KO=bf`, `Kf=bf`: DP = brute force;
     - `L11op`: **Lemma 1.1 in operator form**, Σ_iω(T^{i−1}OT^{k−i}u) with O, T, ω implemented as operators, equals the
       literal K^Ω_k;
     - `P23`, `P23bf`: Prop 2.3, DP and brute force;
     - `T32a`, `T32b`: both formulas of Theorem 3.2; `T32c`: the K_k analogue;
     - `L41s`: Lemma 4.1 summed;
     - `B2`;
     - `TH4`: **Theorem 4**;
     - `TH6a`, `TH6b`: **Theorem 6**, both forms;
     - `Fk`: (F_k).
   - Negative controls, all detected: perturbed b^Ω; Ω without its factors; perturbed ε for one mode; q^{−r_k} in place of
     q^{−r_1}; misplaced O-insertion.
   - My first draft had a truncation bug in my own code. K^♭ carries q^{−r_1}, so tuples whose weight exponent is ≥ L
     still contribute. I found and fixed it; it was not an error of the proof.
4. **`r2_eigen.py 5 22` → `r2_eigen_N5_R22.txt`.** Eigenvectors φ_{N,M}, **N ≤ 5**, all 21 modes, at the 3 points above.
   They are computed in **position space** by the pointwise triangular solve of Tφ = λφ, with R ≤ 22 and no Bailey
   recursion, in exact Fractions.
   - Exact checks, all pass for all 63 mode/point pairs:
     - the top slice equals φ_L(·;xq^{2M}), and H′φ = εφ;
     - Lemma 3.1(a): h is diagonal on Ξ;
     - Lemma 2.1 on φ: T(hφ) = λ(hφ − xOφ);
     - the Bailey expansions of F_0, h_xF_0, q^RF_1 and (Thφ)(·,0), solved on R ≤ L+2M and verified on R ≤ 22 (so they
       are finite, as claimed);
     - **v_x(q^RF_1) = 0**, **(1+x)v_x(h_xF_0) = εv_x(F_0)** and **v_x((Thφ)(·,0)) = λεW/(1+x)**;
     - v_x(F_0) = W_{N,M} of F-all-k Cor 4.2.
   - Numerical checks (num), with ω truncated at R ≤ 22. **Lemma 3.1** ω(hφ) = εω(φ)/(1+x), its P⁻¹ form, Lemma 3.1(b)
     and Lemma 4.3 all hold to relative error at most 5.9e−256, 2.5e−107 and 1.6e−226 at the three points (truncation
     limited). The negative control (ε of another mode) gives errors of order 1.
5. **`r2_symbolic.py` → `r2_symbolic.txt`** (sympy and exact rationals). All OK:
   - S1: the Lemma 2.1 coefficients and WhW⁻¹ = h̃;
   - S2: **Lemma 4.1 for all N, M**, via the identity in u = q^N, v = q^M;
   - S3: Lemma 4.1 symbolically for N, M ≤ 6, and exactly at 924 random rational points for N, M ≤ 20;
   - S4: b^Ω_n definition = factored form, and the x = 1 value, n ≤ 10;
   - S5: the step-1 Bailey pair for general a, n ≤ 7;
   - S6: the Γ^{(1)} re-expansion, and α*_e = λ_e²b^Ω_e for e ≤ 10;
   - S7: **Theorem 6, steps 1–4, at finite n**. β^{(k)}_n equals Σ_rλ_r^{k+1}b^Ω_r/((q)_{n−r}(xq)_{n+r}) exactly at two
     rational points for k ≤ 6, n ≤ 7. The (B1) analogue passes as a positive control. The negative control is
     detected 12/12.
6. **`r2_x1.py 400 14` → `r2_x1_q400_k14.txt`.** The x = 1 evaluations for k ≤ 14, mod q^401, all OK:
   - the theta series Σλ_n^{k+1}b^Ω_n(1) = (q^k,q^{k+4},q^{2k+4};q^{2k+4})_∞;
   - J_k(1) from the multisums C_k(1)/C_1(1);
   - Bressoud's B^{(k)}_{k+1}(1) = (q^k,q^{k+4},q^{2k+4};q^{2k+4})_∞/(q)_∞, which is Theorem 6 at x = 1.
   - Negative control detected.

## What I did not check

- The theorems of `F-all-k.md` (taken as known, per the brief), except:
  - the intermediate facts iv-all-k extracts from their proofs (computation 4);
  - Lemma 4.3 (num, computation 4);
  - the combination in Theorem 6.4/Lemma 7.1, through `T32c`.
- (B1), (B2), Bailey's lemma, Bailey inversion, q-Chu–Vandermonde, Cauchy and the triple product were taken as known.
  (B2), the general-a (B1)-pair, and the J_k(1) and theta evaluations were nevertheless confirmed by computation (`B2`,
  S5, `r2_x1.py`).
- The author's own scripts and certificates. I only read their final status lines ("ALL CHECKS PASSED").
- Novelty claims, and the k = 2 attribution (KRTW Thm 4.10 + Bartholdi); see E8.
- Whether the journal version (Selecta 2026) prints (iv) exactly as arXiv v2. Only the TeX quoted in LITERATURE §1b
  was compared.
