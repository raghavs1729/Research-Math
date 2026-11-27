# Referee report 2: `proofs/KRTW-F-moves.md` (KRTW identities (ii), (iii) for k = 3)

Second independent adversarial referee, 2026-11-26. I did not read `REFEREE-KRTW-F-moves.md`. I did not reuse
`code/krtw_moves_certify.py` or `code/krtw_referee/`. All my scripts are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/krtwF/`
(`chain.py`, `run_chain.py`, `series_test.py`, `xcheck.py`, `evals.py`; arXiv sources in `src/` and `gow/`).

## Verdict

* **Theorem A ((F_2), (F_3), hence (ii) and (iii) for k = 3): PROOF COMPLETE.** I found no mathematical error.
  The findings below are about exposition only.
* **Do the statements match the printed conjecture: YES.** They match arXiv:2404.03851v2 exactly: summation
  ranges, quadratic and linear exponents, and both product sides.
* **Consequence (C), KRTW Conjecture 2.2 at z = q for k = 3: follows as claimed.** It uses GOW Theorem 1.3 with
  n = 2, m = 3, which is a proved theorem.
* Theorem B and Corollary C (lower priority): the labelling is correct. Corollary C(a) is explicitly conditional on
  the Bartlett–Warnaar C₂ Rogers–Selberg identity, and (b) on GOW Theorem 1.3. I tested (B1) numerically for
  k = 3 and it passed.

## 1. Comparison with the printed statements (arXiv:2404.03851v2 source, `KRTW.tex`)

I checked the numbering by counting the shared theorem counter in the source.
* Conjecture 2.2 is `Con_Shun` (Eq_q2), the only numbered environment in §2 after Proposition 2.1.
* Conjecture 4.9 is `Con_Shun2` (§4 runs 4.1 Conj, 4.2 Conj, 4.3 Conj, 4.4 Prop, 4.5 Prop, 4.6 Lemma, 4.7 Conj,
  4.8 Prop, 4.9 Conj).

The displays printed after Conjecture 4.9 read as follows, summed over r_1..r_k ≥ 0 and s_1..s_k ≥ 0 with r_{k+1} = s_0 = 0:

    (ii)  Σ ∏ q^{(r_i+s_i)^2+s_i^2+r_i+2s_i}/((q;q)_{r_i-r_{i+1}}(q^2;q^2)_{s_i-s_{i-1}})
            = (q,q^2,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4};q^{2k+4})_∞ / ((q^2;q^2)_∞ (q;q)_∞)
    (iii) Σ ∏ q^{(r_i+s_i)^2+s_i^2}/(same)
            = (q^{k+1},q^{k+2},q^{k+2},q^{k+3},q^{2k+4},q^{2k+4};q^{2k+4})_∞ / ((q^2;q^2)_∞ (q;q)_∞)

At k = 3 these give (q,q²,q⁸,q⁹,q¹⁰,q¹⁰;q¹⁰) and (q⁴,q⁵,q⁵,q⁶,q¹⁰,q¹⁰;q¹⁰). Both match the document (line 41–42)
symbol for symbol, linear terms included.

KRTW sums over all r_i, s_i ≥ 0 with the convention 1/(q)_{n<0} = 0. That is the same as the document's ordered
ranges r_1 ≥ … ≥ r_k, s_1 ≤ … ≤ s_k. My series test below uses KRTW's form literally and agrees with the products.

## 2. Findings (all minor or informational)

1. **(minor, exposition) Theorem A display, lines 41–42.** The display gives no summation range. It is implicit
   from the definition of K_k on line 29 and agrees with KRTW, as noted above.
2. **(minor, the dependency is only sketched) §2, lines 114–116: J_3 = AG_1(x²;Q).** The chain proves
   K_3 = B_4 · AG_1(x²;Q). To restate this as (F_3), i.e. K_3 = B_4 J_3, one needs J_3 = C_3/C_1 = AG_1(x²;Q).
   The document justifies this in one sentence that points to (B2), which is proved in `KRTW-general.md` §8.
   I re-derived it:
   * Apply (L∞) twice in base Q, relative to a = x², to the unit pair. The α-side is
     (−1)^m x^{4m} q^{5m²−m} (x²;Q)_m (1−x²q^{4m})/((Q;Q)_m(1−x²)), which is (B2) for k = 3.
   * The β-limit is AG_1(x²;Q)/(Q;Q)_∞.

   So the sketch is correct, and I also confirmed C_3 = C_1 · AG_1(x²;Q) numerically. For (ii) and (iii)
   themselves this step is not even needed: K_3 = B_4 · AG_1(x²;Q) together with the two Rogers–Ramanujan
   identities in base q² is enough.
3. **(minor, wording) Two sentences are garbled but the mathematics is right.**
   * Lemma M2's proof (lines 95–97, "…without the entry M′_ab = M_ab+d changed…").
   * The parenthetical on line 33 about powers of x ("at least the number of its nonzero summation indices'
     total").

   The intended statements are correct: the substitution (A,B,c) ↦ (A−c, B−c, c) is a bijection
   {c ≤ min(A,B)} → N³; the factor q^{d·n′_a·n′_b} adds d to M_ab; and w ≥ 1 for every species, so each power of
   x has finitely many terms.
4. **(informational) The specialisation step, lines 171–186, is self-contained and correct.**
   * B_{k+1}(1) and B_{k+1}(q) come from (B1), proved in §4 by Bailey's lemma.
   * J_k(1) and J_k(q) come from (B2), proved in `KRTW-general.md` §8.
   * Bressoud's identities are not needed as an input.
   * I checked the bilateral recombinations by hand: at x = q, the terms n and −n−1 combine to
     Σ_Z (−1)^n q^{(k+2)n²+(k+1)n} = θ(q;q^{2k+4}) (Jacobi triple product).
   * Specialising termwise in Z[x][[q]] is legitimate, because the q-degree of each term is at least (k+2)n².
5. **(informational, not checked) Two claims I did not verify.**
   * The claim that k = 2 follows from KRTW Theorem 4.10 plus "Bartholdi's theorem (arXiv:2611.04020)".
   * The claim that Conjecture 2.2 holds at z = 1 for all k (from `KRTW-general.md`).

   Neither is needed for Theorem A. The document's own chain also covers k = 2.
6. **(informational) Corollary C(b) and consequence (C).** GOW (arXiv:1401.7718) Theorem 1.3 (`Thm_Main3`,
   D^{(2)}_{n+1}, n ≥ 2, κ = 2m+2n) states:

       Σ_{λ_1≤m} q^{2|λ|} P_{2λ}(1,q,…;q^{2n−2}) = (q^κ;q^κ)^n_∞ / ((q²;q²)_∞(q)^{n−1}_∞) · ∏_{i<j} θ(q^{j−i},q^{i+j−1};q^κ)

   With n = 2 and m = k, the right side is (q,q²,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4};q^{2k+4})_∞/((q²;q²)_∞(q)_∞),
   which is the RHS of (ii). Conjecture 2.2 reads:

       Σ_{λ_1≤k} (zq)^{|λ|} P_{2λ}(1,q,…;q²) = Σ ∏ z^{r_i+s_i} q^{(r_i+s_i)²+s_i²+s_i}/(…)

   At z = q, the left side is GOW's left side and the right side is the LHS of (ii). Hence (ii) for k = 3 implies
   Conjecture 2.2 at z = q for k = 3, as the document claims.

## 3. What I recomputed (independent code, exact arithmetic)

* **The moves.**
  * Euler splitting, 1/(q;q)_n = Σ_{u+v=n} q^v/((q²;q²)_u(q²;q²)_v): this is the coefficient of z^n in
    1/(z;q)_∞ = 1/((z;q²)_∞(zq;q²)_∞), and the q^v sits on the second factor. This matches Lemma M1 with
    B′_t = B_a + d.
  * Durfee rectangle identity: standard.
  * The transformation rules of M, B and w in Lemmas M1 and M2: I re-derived them by substitution and they agree
    with the document.
* **The chains.** I implemented split, fuse and defuse from the lemma statements, with exact Fraction
  preconditions on every DEFUSE, and built the data of K_2, K_3, B_3, B_4, C_2 and AG_1 from their definitions.
  * All 9 K-side and 4 B-side moves for k = 3 apply. Every DEFUSE precondition holds.
  * The two rank-8 end data coincide under a species matching, the same one as in the certificate.
  * For k = 2, all moves apply. The species ((n1′|n2″)|n3) is exactly free (M = 2, B = 0, w = 1, base Q, which is
    C_1). After removing it, tensoring with C_2 and one FUSE, the end data coincide.
  * I also checked that F(start datum) equals K_3, respectively B_4, computed directly from the definitions.
* **Every intermediate multisum, as a series.** For k = 3: all 21 data, x⁰…x¹², mod q⁶⁰. For k = 2: x⁰…x¹⁴,
  mod q⁶⁰. Every consecutive pair is equal, except the two tensor steps, which multiply by AG_1 or by C_2/C_1 and
  are not supposed to be equalities. The two ends are equal.
  * No intermediate exponent was negative in the enumerated range.
  * Negative controls fail, as they should: a split without the +d shift, and a fuse without the +d on M_ab.
* **(F) directly from the definitions** (x⁰…x¹⁰, mod q⁴⁰):
  * K_3 = B_4 · AG_1(x²;Q): true.
  * K_2 · C_1 = B_3 · C_2: true.
  * C_3 = C_1 · AG_1(x²;Q): true.
  * (B1) for k = 3: true.
  * Control K_3 = B_3 · AG_1: false, as it should be.
* **(ii) and (iii) directly from KRTW's printed form**, summing over all r_i, s_i ≥ 0 with 1/(q)_{neg} = 0, mod
  **q¹⁵⁰**: true for k = 3, and for k = 2 as a control. Identity (i) is also true, as a further control.
  Deliberately wrong variants fail at q³ and q² respectively: linear term r_i+s_i instead of r_i+2s_i, and product
  (q,q³,…) instead of (q,q²,…).
* **Specialisations:** (q)_∞B_4(1) = θ(q⁵;q¹⁰) and (q)_∞B_4(q) = θ(q;q¹⁰), mod q¹²⁰. The control θ(q⁴;q¹⁰) fails.

## 4. Not checked

* Theorem B for general k beyond the (B1) k = 3 test.
* The Δ_C computation of §5, and the normalisation caveat about the GOW display of the C_n Rogers–Selberg identity.
* The Bartlett–Warnaar citation details.
* §6 (commuting operators, spectral coefficients), and the novelty/prior-art claims.
* `certificates/krtw_four_test.txt` (mod q⁴⁰⁰). I did not rerun it; my own independent check goes to q¹⁵⁰.
