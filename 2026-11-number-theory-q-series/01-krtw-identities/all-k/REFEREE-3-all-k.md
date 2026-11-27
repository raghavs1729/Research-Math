# Referee 3 report on proofs/F-all-k.md (A) and proofs/iv-all-k.md (B)

Status: DONE (2026-11-27). Adversarial read by hand plus independent computer checks.

Independence: I did not open any REFEREE-*.md file other than this one, and did not open, run or import
anything from krtw4/referee/, krtw4/referee2/ or krtw4/code/ (author's code). All scripts are in
krtw4/referee3/, written from scratch for this report; outputs are the .txt files there. Nothing was
committed or pushed; the proof files were not edited.

## 1. Overall verdicts

| claim | verdict |
|---|---|
| (A) F-all-k.md: (F_k) K_k C_1 = B_{k+1} C_k in Z[x][[q]] for all k >= 1 (Theorem 1), hence KRTW (ii), (iii) for all k (Theorem 2) | **No error, no gap. Exposition only** (items 1-6 below). |
| (B) iv-all-k.md: KRTW (iv) for all k (Theorem 5), Theorem 4, Theorem 6 | **No error, no gap. Exposition only** (items 1-3, 7-8 below). |

I checked both proofs line by line. I redid every computation by hand (list in §4) and tested every
intermediate claim with my own code (§5). Both proofs are complete and correct, given the inputs that
the task allows me to take as known: (B1), (B2) of KRTW-F-moves.md, Bailey's lemma, Bailey inversion,
q-Chu-Vandermonde, Cauchy's identity and the Jacobi triple product. Each of these is used correctly and
in the right ring (see items 2 and 3 for two places where the text should say a little more). Nothing I
found changes a statement or a step.

**Statements against KRTW (LITERATURE §1b), letter by letter.** (ii) has the exponent
(r_i+s_i)^2+s_i^2+r_i+2s_i, the product (q,q^2,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4};q^{2k+4})_inf and the
denominator (q^2;q^2)_inf (q;q)_inf. (iii) has the exponent (r_i+s_i)^2+s_i^2, the product
(q^{k+1},q^{k+2},q^{k+2},q^{k+3},q^{2k+4},q^{2k+4};q^{2k+4})_inf and the same denominator. (iv) has the
same exponent, the product (q^k,q^{k+1},q^{k+3},q^{k+4},q^{2k+4},q^{2k+4};q^{2k+4})_inf, the denominator
(q^2;q^2)_inf (q;q)_inf, and Omega = sum_{i=1}^{k-1} q^{r_i+2s_i}(1-q^{2s_{i+1}-2s_i}) + q^{r_k+2s_k}. The
conventions are r_{k+1} = s_0 = 0 with the sum over all r_i, s_i >= 0. All of these match the print
exactly; the only difference is item 1. The claims "(iii) = K_k(1)", "(ii) = K_k(q)" and "(iv) = K^Om_k(1)"
are correct: with r increasing or s decreasing, some factor 1/(q)_m or 1/(Q;Q)_m with m < 0 vanishes.

## 2. Section-by-section verdict table

| file, place | content | verdict | how checked |
|---|---|---|---|
| A §0 | statements of (ii), (iii); K_k, B_{k+1}, C_k, J_k; Thm 1, Thm 2 | correct (item 1) | letter-by-letter vs KRTW; series_q |
| A §0 proof of Thm 2 | via KRTW-F-moves §3 | correct | the four theta evaluations at x = 1, q redone by hand |
| A §0 Cor. 3 | conditional HL statements | not checked | (cited external theorems) |
| A §1 | ring F = Q(q)((x)), passage to Q[[q]][[x]] and Z[x][[q]] | correct (item 2) | by hand |
| A §1 | V, T, u, omega, convergence, (1.1), (1.2), (1.3) | correct | by hand; (1.1) is what series_q/bivar compute |
| A §2 | H̃′, G̃; normal-ordered coefficients of H′E, EH̃′, GE, EG̃; restriction to N^2 | correct (item 5) | all four coefficient computations by hand; commutators.py (exact) |
| A (2.1) | H′, G in terms of h_{xY} | correct | by hand (h^2 expanded) |
| A §3 (B-i)-(B-iii) | Bailey's lemma (L∞), Cauchy, eigen-relation of h_a | correct | by hand; bailey_facts.py |
| A §3 (B-iv) | change of a: support, sum, φ_d(R;aq^2) = (aq)_2 φ_{d+1}(R+1;a); explicit Γ | correct (item 3) | partial-fraction argument checked; Γ symbolic (algebra.py), pointwise (bailey_facts.py) |
| A §3 (B-v), v_a rules | three-term relations, coefficient sums, h_a τ F identity | correct | by hand; symbolic; n >= 1 / κ_0 = 0 needed exactly as stated (see §3 below) |
| A §3 basis, (3.1) | Ξ-basis, triangular action of T | correct | by hand; (3.1) pointwise (bailey_facts.py) |
| A Prop 3.1 | existence, uniqueness, supports | correct (item 4) | by hand; eigen.py (position-space solve agrees) |
| A (3.2) | shift x -> xq^2 | correct | by hand; eigen.py exact |
| A Lemma 3.2 | H′φ = εφ, Gφ = gφ | correct | by hand; eigen.py exact, N <= 5 |
| A Prop 4.1 | three v_x equations, elimination, factorisation | correct | by hand; algebra.py (generic symbols); eigen.py (recursion, exact) |
| A Cor 4.2 | closed form of W_{N,M} | correct | algebra.py (iteration, N <= 10); eigen.py (against eigenvectors, N <= 5) |
| A Lemma 4.3 | ω(φ) = λW/(xq)_inf | correct | by hand; eigen.py (150-digit numerics) |
| A §5 | T♯ = W^{-1}T^tW, (5.1), transposed commutation, base-Q facts | correct | (5.1) derived by hand from the kernel of H′, G; commutators.py; algebra.py |
| A Prop 5.1, Cor 5.2 | recursion of v_ρ, Fricke factorisation, closed form of V♯ | correct | by hand; algebra.py; eigen.py |
| A Lemma 5.3 | ψ(u) = λV♯/(x^2Q;Q)_inf | correct | by hand; eigen.py numerics |
| A Lemma 6.1 | biorthogonality, diagonal value | correct | by hand; eigen.py exact, all 21^2 pairs N <= 5 |
| A Lemma 6.2 | valuations (a), (b), (c) | correct (item 2 for c1, c2) | by hand; valuations.py (x indeterminate) |
| A Lemma 6.3 | completeness | correct | by hand (min-ord argument re-derived) |
| A Thm 6.4 | expansion of u, all interchanges, spectral formula for K_k | correct (item 6) | by hand; expansion.py (exact mod x^14); bivar.py (spectral formula mod q^100) |
| A §7 | t_{N,M} closed form, Lemma 7.1, proof of Thm 1 | correct | own proof of Lemma 7.1 (§4); algebra.py N <= 10; bivar.py |
| B §0 | statement of (iv), K^Ω, B^(k), b^Ω and its closed form, k = 1 remarks | correct (item 1) | by hand; algebra.py; series_q |
| B Lemma 1.1 | operator form of the Ω-insertion | correct | by hand; lemma11_theta.py (bivariate, Ω literal) |
| B Lemma 2.1 | [T,h] = -xOT | correct | coefficients by hand; commutators.py exact |
| B Lemma 2.2 | ω(hF) = ω(P^{-1}F) | correct | by hand; commutators.py |
| B Prop 2.3 | telescoping K^Ω = x^{-1}(K♭ - K) | correct | by hand (indices checked); bivar.py |
| B Lemma 3.1 | ω∘h = ε/(1+x) ω on eigenvectors, incl. M = 0 | correct | by hand; eigen.py (both the slice-0 statement and the lemma itself) |
| B Thm 3.2 | spectral weights of K♭, K^Ω | correct | by hand; bivar.py (both formulas) |
| B Lemma 4.1 | factorisation of the weights | correct | own derivation (§4); algebra.py N <= 10 and generic u, v |
| B Thm 4 | identity, symmetrisation, ring passage | correct | by hand; bivar.py |
| B Thm 5 | x = 1, b^Ω_n(1), triple product, J_k(1) | correct | by hand; lemma11_theta.py; series_q (iv) mod q^400 |
| B Thm 6 | Bailey pair for general a, two (L∞) at xq, change of a, k-1 (L∞) at x, limit | correct (items 3, 7) | every step by hand; algebra.py (α* = λ^2 b^Ω, general-a pair); bivar.py |

## 3. Problems (all of severity "exposition"; none affects a statement or a step)

1. **Convention for the base-q^2 Pochhammer symbols (exposition).** A §0 line 15 and B §0 lines 16-17 state
   "1/(q;q)_m := 0 for m < 0" only. The sums also contain 1/(q^2;q^2)_{s_i-s_{i-1}}, and the same
   convention is needed there. It is the standard extension (q^2;q^2)_{-m} = 1/(q^{2-2m};q^2)_m, and KRTW's
   display is silent on both. *Fix:* "1/(q;q)_m = 1/(q^2;q^2)_m := 0 for m < 0".

2. **Infinite products in F = Q(q)((x)) (exposition).** A §1 (lines 55-59) says that rational functions
   are elements of F, but infinite products are not. (xq;q)_inf, (x^2Q;Q)_inf and 1/(aq;q)_inf with
   a = xq^{2s} appear in (B-ii) line 121, Lemma 4.3, Lemma 5.3, Thm 6.4 and §7, and in B Thm 3.2 and Thm 4.
   They do **not** converge x-adically when Q(q) carries the discrete topology: the x-coefficient of
   prod_{i<=n}(1-xq^i) is -(q+...+q^n), which never stabilises. They have to be defined as elements of F,
   for example by their Euler expansions sum (-1)^n q^{n(n+1)/2} x^n/(q)_n, which have coefficients in
   Q(q) ∩ Q[[q]]. Equivalently, define them as the preimages under the injective coefficientwise Taylor
   map, which A §7 (lines 322-324) already uses for (B1), (B2). With that definition, (B-ii) and the
   product manipulations hold in F by exactly the transfer argument of §7. *Fix:* one sentence in §1.
   Related: A Lemma 6.2(b) line 280 uses "ord c1, c2 >= 0" without writing c1, c2 down. They are
   c1 = (1-aQ)/(1-aQ^{2e+1}) and c2 = -aQ^{2e+1}c1 (a = x^2q^{2ρ}), from the c = 1 partial fractions;
   quoting them would make (b) self-contained.

3. **The explicit Γ^{(c)}_{d,e} formula is cited only to code (exposition).** A (B-iv) lines 123-125 gives
   the general Γ formula with the justification "Bailey inversion and q-Chu-Vandermonde; code/tower.py".
   The text does prove the properties that are used: support, coefficient sum (aq)_c, and
   Γ^{(0)}_{d,d} = 1. But B §6 item 3 (lines 205-206) takes the c = 1 coefficients "from the formula
   there". These follow in two lines from the numerators over (q)_{R-d}(aq)_{R+d+1}:
   1 - aq = A(1 - aq^{d+1}z) + B(1 - q^{-d}z) gives A = (1-aq)/(1-aq^{2d+1}) and B = -aq^{2d+1}A.
   *Fix:* put this derivation in the text, or drop the general formula. I verified the general formula
   anyway: symbolically for c <= 4, d <= 3, and pointwise and exactly for c <= 5, d <= 4.

4. **Last clause of Prop. 3.1 (exposition).** In A lines 155-156 the proof of "every T-eigenvector with
   finitely many Ξ-coefficients is a multiple of some φ_{N,M}" stops at "(e,σ) = (L,M)". It should add:
   "hence F/f_{L,M} = Ξ_{L,M} + (lower slices), and by the uniqueness part F = f_{L,M}φ_{N,M}". The same
   applies to the T♯ analogue in §5.

5. **Domain of the normal-ordering argument (exposition).** A Prop 2.1 proof, lines 87-93. "Equal
   normal-ordered coefficients give equal operators on W𝒱" is true. Strictly, though, it is applied to
   functions on Z^2 such as W H′F, which is supported in R >= 0, S >= -1 and has ord >= R + 2S + c. It
   would help to name that class and say that the needed rearrangements are finite recombinations of
   x-adically convergent sums. The restriction-to-N^2 argument that follows is complete; I re-derived it.

6. **Indexing in the summability families (exposition).** A Thm 6.4 proof, lines 299-301. The families
   are written with λ_b^{k-1}. The induction step T(sum ũ λ^{k-2} φ) uses the same family with exponent
   k-2, which has the same order bound because ord λ >= 0. This is harmless; one could say "for every
   j >= 0, with λ_b^j".

7. **Ring in Theorem 6 (exposition).** B line 41 does not name the ring. The proof works in 𝓡 = Z[[q]][[x]]
   (or F). The consequence K^Ω_k = B^(k)_{k+1} J_k then holds in Z[x][[q]] by Theorem 4, because
   (xq)_inf is a unit there. *Fix:* add "in Z[x][[q]]".

8. **Scope (not a problem).** The novelty and literature statements (A line 41; B lines 45-51, including
   the attribution of k = 2 to Bartholdi) and Corollary 3 of A are outside what I checked.

**Where n >= 1 is needed (asked in the brief).** It is needed exactly in (B-v), that is, in the v_a
rules for F with κ_0 = 0. For κ_0 ≠ 0 the rules really fail. τφ_0 is not a finite combination of the
φ_e, so v_a is undefined: bailey_facts.txt finds 15 nonzero coefficients among e <= 14. And the
pointwise identity h_a(τF) = F(·+1) + aq^{2R}F(·+1) - q^{-1}F fails at R = 0 when F(0) ≠ 0. Every use
is legitimate:
- Prop 4.1 applies the rules to G = sum κ′_e φ_{e+1}, so κ_0 = 0.
- Lemma 3.2 and Lemma 6.2(a) use τφ_{e+1} with n = e+1 >= 1.
- B Lemma 3.1 uses only Prop 4.1's equations.
- F_0, which may contain φ_0 (when L = 0), is only ever acted on by h_x, through (B-iii), which holds
  for all e.

## 4. Computations redone by hand (all agree with the text)

- The coefficients of S^j U^b (times (q)_j(Q;Q)_b) in H′E, EH̃′, GE and EG̃ (A Prop 2.1), and in hE, Eh̃
  and xPYUE (B Lemma 2.1). This includes the identities 1 - (1+q)(1-t) + q(1-t)(1-t/q) = t^2 and
  -(1+1/q)(q/t^2)(1-t) + (q/t^2)(1-t)(1-t/q) + 1/t^2 = 1 (t = q^j) in the G computation.
- H̃′ and G̃ from (1.3) (including WS^2 = x^2q^{-4}P^4Y^2S^2W), and (2.1) by expanding h_{xY}^2.
- (B-i) to (B-v), including α_n + β_n + γ_n = 0 and the boundary cases R < n.
- Lemma 3.2: the slice-M values ε, g.
- The three v_x equations of Prop 4.1, the elimination, and both factorisations.
- (5.1): H′♯ and G♯ derived from the kernel H′(b,a)w(b)/w(a).
- Prop 5.1: the slice equations, the elimination with c = q^ρ + q^{-ρ}, and the Fricke factorisation
  (c - A)(c - B) with A = xq^{N+M} + x^{-1}q^{-N-M} and B = q^L + q^{-L}.
- Cor 5.2: (q;q^2)_L/((q)_L(q^L;q)_L) = (1+q^L)/(Q;Q)_L.
- Lemma 6.2: which index is reached only through γ, and the x-coefficient of the bracket.
- Lemma 6.3: the min-ord argument.
- Theorem 6.4: each summability bound.
- **Lemma 7.1, my own proof for all N >= M >= 0.**
  - Monomials: both sides give (-1)^{N+M} x^{N+M} q^{2L^2+4LM+4M^2-L-M}.
  - q-only factors: (Q;Q)_M (Q^{L+1};Q)_M (Q;Q)_L = (Q;Q)_M (Q;Q)_N. The factor (1+q^L) matches
    q^{-M}(q^N+q^M).
  - Base-Q factors in x^2: (x^2Q;Q)_L (x^2q^{2L+2};Q)_{2M} / (x^2q^{2N};Q)_M = (x^2Q;Q)_{N-1}(1-x^2q^{2N+2M}).
  - What remains is (x;q)_{2N+1}(1 - xq^{N+M}) = (1 - xq^{2N}) (x;q)_{N+M+1} (xq^{N+M};q)_L. This holds
    because (x;q)_{N+M+1}(xq^{N+M};q)_L = (x;q)_{2N}(1 - xq^{N+M}).
  - For L = 0, (1+q^0)/(1+δ) = 1, and the identity reduces to (x;q)_{2N+1}(1-xq^{2N}) = (1-xq^{2N})^2 (x;q)_{2N}.
- B Lemma 4.1: ε - μ_N - 2x = (u-v)(1+xu^2)(1-xuv)/(u^2v) by expansion, and
  β^B_n(1+xq^{2n})q^{-n} = (1+x)β^J_n.
- B Lemma 1.1: the bond identity (1-Q^m)/(Q;Q)_m = 1/(Q;Q)_{m-1} for m >= 0.
- B Prop 2.3: the telescoping indices (a_{i-1} = ω(T^{i-1} h T^{k-i} u)).
- B Thm 5: b^Ω_n(1) = (-1)^n(q^{n^2-2n} + q^{n^2+2n}); the triple product with A = k+2, B = -2; J_k(1).
- B Thm 6:
  - the general-a Bailey pair, via inversion and q-Chu-Vandermonde, with (-q^{1-n}/a)_n = (-a)_n a^{-n} q^{-n(n-1)/2};
  - two (L∞) steps at a = xq;
  - the c = 1 change of a (item 3);
  - α*_e, including the exponent 3e^2 - 2e and 1 - y + xy^2 - x^3y^3 = (1-xy)(1+(x-1)y+x^2y^2);
  - the index bookkeeping of the k-1 further steps and the limit;
  - the k = 1 reduction to (B1) at x -> xq.
- KRTW-F-moves §3: the evaluations (q)_inf B_{k+1}(1), (q)_inf B_{k+1}(q), (Q;Q)_inf J_k(1) and
  (Q;Q)_inf J_k(q) as theta functions (n -> -n-1 pairings), which give the right sides of (iii), (ii).

## 5. Independent computations (scripts in krtw4/referee3/, outputs .txt there)

All arithmetic is exact (integers or rationals) unless marked numeric. Negative controls are in
parentheses; all were detected.

1. **series_q.py 400 10 -> series_q_400.txt**, and **series_q.py 300 12 -> series_q_300_k12.txt**.
   - (ii), (iii) and (iv) are evaluated directly from the printed multisums, by nested summation over the
     levels. Ω is expanded literally into its 2k-1 monomials. The results are compared with the printed
     products.
   - **PASS for k = 1..10 mod q^400 and k = 1..12 mod q^300.** KRTW's (i), which is not claimed, also agrees.
   - Negative controls: sum (iv) against product (iii); Ω without its "-q^{2s_{i+1}-2s_i}" part; product
     (iii) at k+1; sum (i) against product (iii).
   - brute_check.py -> brute_check.txt: the evaluator agrees with brute-force enumeration over all tuples
     (k <= 3, mod q^34).
2. **bivar.py 100 8 -> bivar_100.txt**, and **bivar.py 60 10 -> bivar_60_k10.txt**. Identities in
   Z[x][[q]] mod q^100 for k <= 8 (mod q^60 for k <= 10), all powers of x. Every monomial x^a q^b has
   a <= b, except in K♭, whose larger state space is accounted for, so the truncation is exact.
   - (F_k): K_k C_1 = B_{k+1} C_k. **PASS.**
   - K^Ω_k C_1 = B^(k)_{k+1} C_k (Ω literal). **PASS.**
   - (B1), (B2); B Thm 4; B Thm 6; B Prop 2.3. **PASS.**
   - A Thm 6.4's spectral formula from the closed forms of V and V♯, **without using Lemma 7.1**. **PASS.**
   - Both formulas of B Thm 3.2. **PASS.**
   - Negative controls: B_{k+2} in place of B_{k+1}; one spectral weight perturbed; Ω's last term dropped.
3. **commutators.py -> commutators.txt.** Exact on point masses at (q,x) = (2/7,3/5), (-3/4,5/11), (5/3,-2/9).
   - [T,H′] = [T,G] = 0 and Th - hT + xOT = 0. **PASS.**
   - (5.1) equals W^{-1}H′^tW and W^{-1}G^tW entrywise, and T♯ = W^{-1}T^tW. **PASS.**
   - [T♯,H′♯] = [T♯,G♯] = 0, and B Lemma 2.2. **PASS.**
   - Negative controls: H′ with its U-coefficient doubled; h with xP in place of xPY.
4. **eigen.py 5 -> eigen_N5.txt** (and eigen_N4.txt). My own eigenvectors of T and T♯, computed by a
   pointwise triangular solve in position space at two rational points, for all 21 modes N <= 5, exact.
   **ALL PASS:**
   - the supports, and the supports of the Bailey coefficients (L <= e <= L+2(M-S); M <= e <= M+L-ρ);
   - Lemma 3.2 and (3.2);
   - Prop 4.1 (the recursion and all three slice-0 statements at infinity) and Cor 4.2;
   - Prop 5.1 and Cor 5.2, and the eigen-relations of H′♯, G♯;
   - Lemma 6.1 for all 21^2 pairs;
   - Lemma 4.3, Lemma 5.3 and B Lemma 3.1, numeric at 150 digits from exact partial sums, relative
     error < 1e-140.
5. **valuations.py 5 -> valuations_N5.txt.** Lemma 6.2 (a), (b), (c) with x an indeterminate (exact in Q(x)),
   at q = 2/7 and -3/5, N <= 5. **PASS.** The bound 2(M-S) in (a) is attained on every slice. In (b) the
   observed orders are L - ρ, which is stronger than the >= 0 that is needed.
6. **expansion.py 14 -> expansion_x14.txt.** u = sum ũ_{N,M} φ_{N,M} (Thm 6.4), exactly in Q[[x]] mod x^14
   at q = 2/7 and -3/5, for (R,S) in [0,2]^2. **PASS.** (Dropping any one contributing mode breaks it.)
   This is a direct test of the completeness argument.
7. **algebra.py 10 -> algebra_N10.txt.** **ALL PASS:**
   - in fmpz_mpoly (exact cross-multiplication), for all N >= M >= 0 with N <= 10: Lemma 7.1 against the
     *definitions* of β^B, β^J; the closed form of t; Cor 4.2 from iterating Prop 4.1; Cor 5.2 from
     Prop 5.1; B Lemma 4.1;
   - the closed form of b^Ω and α*_e = λ_e^2 b^Ω_e for n <= 24;
   - the general-a Bailey pair for a in {xq, x}, n <= 6;
   - in sympy with generic symbols: the Γ formula and its sum; (B-v); the base-Q facts; the eliminations
     and factorisations of Prop 4.1 and Prop 5.1; the factorisations of B Lemma 4.1 and (B4).
   - Negative control: L7.1 with (xq;q)_3 in place of (xq;q)_4.
8. **bailey_facts.py -> bailey_facts.txt.** **PASS:**
   - (B-i), (B-iii), (B-iv) and (B-v), exactly; (B-ii) to 100 digits;
   - the v_a rules, and their failure when κ_0 ≠ 0 (see §3);
   - (3.1), and the displayed T♯Ξ♯ expansion of §5, pointwise.
9. **lemma11_theta.py 60 6 300 -> lemma11_theta.txt.** **PASS:**
   - B Lemma 1.1 on its own: K^Ω_k with Ω literal equals sum_i ω(T^{i-1}OT^{k-i}u), bivariate mod q^60,
     k <= 6;
   - the x = 1 evaluations of B Thm 5: sum λ^{k+1}b^Ω = (q^k,q^{k+4},q^{2k+4};q^{2k+4})_inf and
     (Q;Q)_inf J_k(1) = (q^{k+1},q^{k+3},q^{2k+4};q^{2k+4})_inf, mod q^300, k <= 10.

Reproduce: `cd krtw4/referee3 && python3 series_q.py 400 10 && python3 bivar.py 100 8 && python3 commutators.py
&& python3 eigen.py 5 && python3 valuations.py 5 && python3 expansion.py 14 && python3 algebra.py 10 &&
python3 bailey_facts.py && python3 lemma11_theta.py 60 6 300` (python-flint 0.9, sympy 1.14, mpmath 1.3;
a few minutes in total on 4 cores).

## 6. What I did not check

- The proofs of (B1) and (B2) in KRTW-F-moves.md. The task allows them as known. I did re-derive the
  general-a Bailey pair behind (B1), and verified (B1) and (B2) as series (§5, item 2). I did not read
  KRTW-general.md §8, where (B2) is said to be proved.
- Corollary 3 of A, which depends on Griffin–Ono–Warnaar Thm 1.3 and the C_2 Rogers–Selberg identity of
  Bartlett–Warnaar. The novelty and literature claims (§3, item 8) are not checked either.
- The author's certificates and code (by instruction), and the other referees' reports.
- KRTW's own derivation of (i)-(iv) from Conjecture 4.9 and Primc's conjecture. It is irrelevant here,
  since (ii)-(iv) are proved directly as q-series identities.
- No formal (proof-assistant) verification. The computer checks cover finite ranges (k <= 12,
  N <= 10, and so on); the all-k statements rest on the hand verification of §4.
