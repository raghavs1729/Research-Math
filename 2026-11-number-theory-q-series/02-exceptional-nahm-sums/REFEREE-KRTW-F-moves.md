# Referee report: proofs/KRTW-F-moves.md

Referee pass of 2026-11-26. I did not edit the document or its code. All my code is in
/tmp/claude-0/-home-user-weighted-prime-races/6df0c2f4-62d2-56bc-82b9-a48bb6ace891/scratchpad/referee_F/
and imports nothing from the project. Python 3 + python-flint, using exact integer and rational arithmetic throughout.

## Verdicts

* **Theorem A** (F_2), (F_3), hence KRTW (ii), (iii) for k = 3: **proof complete.** I re-derived Lemmas M1 and M2
  by hand. I replayed both chains with my own implementation of the moves, and every intermediate multisum has
  the right series. The passage from (F_k) to (ii), (iii) is correct. I found only presentational points (F5, F6).
* **Theorem B** (B1)–(B5), all k: **proof complete.** I re-derived every Bailey-pair step and the limiting forms
  of Bailey's lemma, and checked each item numerically for k ≤ 5.
* **Corollary C**: **proof complete, given the two cited theorems.** The caveat about the GOW normalisation is
  correct: the printed GOW display is false without the products, and true with them for generic (y_1, y_2).
  The 2k-fold multisum is right. One side formula in §5 (the k = 1 paragraph) is wrong as written (F2), and one
  citation is imprecise (F4).
* **§6 item 2** (not used in any proof): **error.** The operator H as printed does not commute with T. A corrected
  H′ does, and the "no other operator" claim is false (F1).

## Findings, by severity

**F1 (error, §6 item 2, lines 266–276).** With T, P, Y, 𝖲, 𝖴 exactly as defined there, the printed operator

    H = P⁻¹Y⁻¹(1 − q⁻¹𝖲) + q³x P⁻¹(1 − q⁻¹𝖲) + q⁻¹x P + q²x² PY(1 − q⁻²𝖴)

does **not** commute with T. Test: apply both sides to point masses F = δ_(r0,s0) with r0, s0 ∈ [−3, 4], evaluate at
(R, S) ∈ [−4, 6]², at two rational points (q, x) = (3/7, 5/11) and (−2/9, 13/4), in exact arithmetic. Result:
[H, T]F ≠ 0 at 2224 of 7744 sample points. Already for F = δ_(1,1) it fails at every (R, S) with R ≥ 2, S ≤ 0
(hamtest.py, hamtest2.py).

The commutant of T in the stated span {P^aY^b𝖲^c𝖴^d : |a|,|b| ≤ 2, 0 ≤ c,d ≤ 2} (225 monomials, 15876 exact
equations, commutant.py) has dimension **4**, not 1. It contains

    H′ = P⁻¹Y⁻¹(1 − 𝖲) + x P⁻¹(1 − 𝖲) + x P + x² PY(1 − 𝖴)

which commutes with T at both rational points (0 failures of 7744). It also contains H′² and a further,
independent operator

    G = P⁻²Y⁻¹(1 − 𝖲)(1 − q𝖲) + Y⁻¹ − x(1 + q⁻¹)𝖲 + x² Y(1 − 𝖴) + x² P²Y(1 − 𝖴),

which also gives 0 failures at both points (hamtest3.py). So:
(a) the displayed H should be H′, or else a different convention for T, 𝖲, 𝖴 must be stated;
(b) "with no other operator in the span … (except constants)" is false. H′² lies in that span trivially, and G is
a second, independent commuting operator. Two commuting operators fit the "rank two" interpretation better than one.

The script `ham.py` cited for this item is not in the repository, so the claim cannot be reproduced as cited.

**F2 (error in a side remark, §5 last paragraph, lines 254–256).** The k = 1 formula
"Σ_{a+b=N} x^N q^{a²+b²−N} ω(a−b) 2^{[a≠b]}/((Q)_a(Q)_b) — the case k = 1 of the display above" is false as written.
* Over ordered pairs (a, b), the factor 2^{[a≠b]} double counts. At N = 1 the formula gives 2/(1−q), but
  P′_{(2)}(√x, −√x) = x/(1−q).
* Restricted to a ≥ b, it equals P′_{(2^N)}(√x, −√x; q) = P_{(2^N)}(X_x; q) (checked by Kostka–Foulkes for
  N ≤ 7, mod q³⁰).
* The x^N-coefficient of the k = 1 case of the display is q^N times that. It is Σ_{a+b=N, ordered} q^{a²+b²} ω(a−b)/((Q)_a(Q)_b),
  with no q^{−N} and no 2^{[a≠b]}.

Fix: write "Σ_{a+b=N, a≥b}" and "= P_{(2^N)}(X_x;q)", or drop both the 2^{[a≠b]} and the q^{−N}.

**F3 (minor, reproducibility, §5 lines 235–236).** The comparison of the BW multisum with K_k is quoted from a
scratch script `bw2.py` that is not in the repository. I re-ran it independently (bwmulti.py). For k = 1, 2, 3
the multisum equals K_k, x⁰…x¹⁰ mod q³⁰. The claim holds, but the document should cite a script that is in
the repository.

**F4 (minor, citation, §5 "The cited theorems").** According to BW_nl.tex, (RS) comes from Theorem Thm_CnAndrews
(nonterminating limit, right-hand side in §"The right-hand side of the C_n Andrews transformation") combined with
Corollary Cor_conjectural, the Littlewood-type identity (Eq_sumM), at w, z → 0. In that limit only even λ survive,
with h = 1, and the prefactor is D → ∏(qx_i²)_∞ ∏_{i<j}(qx_ix_j)_∞. Theorem 3.1 enters only through
Cor_conjectural. The cleaner citation is "Thm_CnAndrews + (Eq_Rk) with b, c → ∞".
The normalisation caveat itself is right: BW's (Eq_Rk) carries the factor D, and GOW's display omits it.

**F5 (minor, presentation, §1 and §3).** Lemma M1 does not say which of a′, a″ is the copy with B shifted by d.
The chains only work, and the final data only match, if a′ keeps B_a and a″ gets B_a + d (this matches
`code/krtw_moves_certify.py`, function split). Say so in §3.

**F6 (minor, §3 "From (F_k) to (ii), (iii)", lines 160–170).**
* The ring: the text works in 𝓡 = Z[[q]][[x]] and then "lets x → 1". Substituting x = 1 or x = q is not defined
  on Z[[q]][[x]]. It is legitimate here because every series involved lies in the subring Z[x][[q]], where each
  q-coefficient is a polynomial in x (all exponents grow with the indices). It should be said once.
* The J factor: the remark "→ 2 for n ≥ 1 … similarly for J" is inaccurate for J. At x = 1 the factor
  (x²;Q)_m(1−x²q^{4m})/((Q;Q)_m(1−x²)) becomes 1 + q^{2m}, not 2. The theta-function identities stated are
  nevertheless correct; I re-derived all four by hand and checked them through (ii), (iii) mod q²⁰⁰.

**F7 (remark, Theorem A).** Correct as stated. The DEFUSE preconditions are equalities of rationals, and I checked
them exactly in my own code. Every intermediate datum has nonnegative exponents on the enumerated range. The
regrouping is legitimate in Z((q))[[x]], since w ≥ 1 leaves finitely many terms per power of x.

## Hand re-derivations (all agree with the text)

* **Lemma M1.** 1/(z;p)_∞ = 1/((z;p²)_∞(zp;p²)_∞) gives 1/(p;p)_n = Σ_{u+v=n} p^v/((p²;p²)_u(p²;p²)_v). The
  quadratic form in n_a = u+v gives M′_aa = M′_tt = M′_at = M_aa and B′_t = B_a + d.
* **Lemma M2.** The Durfee identity holds (checked for A, B ≤ 8 and d = 1, 2, mod q⁸⁰). Substituting A = a′+c,
  B = b′+c gives exactly the stated M′, B′, w′. The factor p^{a′b′} gives M′_ab = M_ab + d. The DEFUSE
  precondition is the inverse map.
* **Limiting cases of Bailey's lemma.** Take Bailey's lemma with ρ₁ → ∞:
  α′_n = (ρ₂)_n(−1)^n q^{n(n−1)/2}(aq/ρ₂)^n α_n/(aq/ρ₂)_n, and
  β′_n = Σ_j (ρ₂)_j(−1)^j q^{j(j−1)/2}(aq/ρ₂)^j β_j/((q)_{n−j}(aq/ρ₂)_n).
  Then ρ₂ = −√(aq) gives (S2), ρ₂ = −√a gives (S2′), and ρ₂ → ∞ gives (L∞). The unit pair (U) is standard.
  The limit (lim) holds coefficientwise in x.
* **(B1).** The inversion formula is correct. q-Chu–Vandermonde with c = −q, b = aq^n gives the stated α_n.
  (L∞)^{k+1} then gives (B1).
* **(B2).** In (S2)^{k+1}(U), base Q, a = x², the factors (−xq;Q)_j cancel step by step. This gives
  C_k/((Q)_∞(−xq;Q)_∞) on the β side, hence J_k = C_k/C_1.
* **(B3).** The α of (S2′)(U) simplifies via (1+x)(1−x²q^{4n})/((1+xq^{2n})(1−x²)) = (1−xq^{2n})/(1−x). The
  prefactor (x²Q;Q)_∞/((xq;q)_∞(−xq;Q)_∞) equals (−xQ;Q)_∞.
* **(B5).** Both statements follow as written.
* **(B4).** The identity for E checks symbolically (sympy).
* **J_3 = AG_1(x²;Q).** (L∞)²(U) in base Q with a = x² gives α_n = x^{4n}q^{5n²−n}(…), the (B2) term for k = 3.
* **§5(a).** L^{(0)}_k(√x, −√x; q) = L_k(x): checked by hand, and symbolically (sympy) for k = 1, 2 and
  n, m ≤ 2.
* **§5(d).** GOW Thm_Main3 (Theorem 1.3) with n = 2 and κ = 2m+4 gives exactly the product side of (ii).
  KRTW Con_Shun (Conjecture 2.2) at z = q has exactly the left side of (ii) as its right-hand side. So
  "(ii) ⇔ Conjecture 2.2 at z = q", given GOW, is correct (m = k ≥ 1, n = 2 ≥ 2).

## Table of independent computations

| # | What | Range | Result |
|---|------|-------|--------|
| 1 | All 30 intermediate data of both chains, k = 2 and k = 3 (chains.py): start, after SPLIT, each FUSE/DEFUSE, the k = 2 datum minus its free species (vs B_3/C_1), the direct sums with C_2 and AG_1. Each is evaluated by brute force and compared with a direct evaluation of the target | x⁰…x¹⁶ mod q¹⁰⁰ (also x¹² mod q⁶⁰) | all equal; every DEFUSE precondition holds exactly; final data equal up to permutation (rank 5 and rank 8) |
| 2 | K_k C_1 = B_{k+1} C_k, brute force without any moves (allx.py) | k = 2, 3, all powers of x, mod q¹²⁰ | true |
| 3 | J_3 = C_3/C_1 = AG_1(x²;Q) | x⁰…x¹⁶ mod q¹⁰⁰ | true |
| 4 | (ii), (iii) at x = q, x = 1 against the KRTW products, by direct summation (xq.py) | k = 1, 2, 3, 4, mod q²⁰⁰ | true (k = 4 numerical only) |
| 5 | Products of Theorem A against KRTW.tex (after \label{Con_Shun2}) | k = 3 | agree |
| 6 | Bailey pair β_n = 1/(Q;Q)_n relative to x; (B1)–(B4) for k = 1..5; (B5) for k = 2, 4; J_k = AG_{(k−1)/2}(x²;Q) for odd k (thmB.py), with a negative control (thmB_neg.py) | x⁰…x¹⁴ mod q⁹⁰ | all true; negative control false |
| 7 | E identity of (B4); L^{(0)} specialisation (sym.py) | symbolic | true |
| 8 | Generic-y C₂ Rogers–Selberg, y = (2t, −3t/5), left side from Kostka–Foulkes via charge (rs_generic.py) | m = 1, 2, 3, t-degree ≤ 12, mod q²⁴ | true with the products; false as printed in GOW |
| 9 | 2k-fold BW multisum at (√x, −√x) against K_k (bwmulti.py) | k = 1, 2, 3, x⁰…x¹⁰ mod q³⁰ | true |
| 10 | (HL_k) directly from Kostka–Foulkes (hl_direct.py) | k = 1, 2, 3, x⁰…x⁷ mod q³⁰ | true |
| 11 | k = 1 formula of the last paragraph of §5 | N ≤ 7 | false as written; true restricted to a ≥ b, as P′_{(2^N)} (F2) |
| 12 | Durfee and Euler identities; closed form (3.1); K_2 = (−xq;q)_∞ 𝔍^{(1)}_1 𝔍^{(2)}_1 (misc.py) | A, B ≤ 8; x⁰…x¹² mod q⁸⁰ | true |
| 13 | §6 item 2: [H, T] on finite-support functions on Z², exact arithmetic (hamtest*.py, commutant.py) | 2 rational points, 7744 evaluations each | printed H fails; H′ and G commute; commutant dimension 4 (F1) |

## Novelty (item 5)

* KRTW.tex states (ii), (iii) as conjectural for all k, and says the k = 2 case can be completed (Thm_kis2).
* Two web searches (standard engine) found no proof of (ii), (iii) or of Con_Shun / Con_Shun2 for k ≥ 3. They
  turned up only KRTW itself and unrelated 2026 papers (2609.07822 Kanade–Russell mod 12; 2610.18666
  Takigiku–Tsuchioka level 2).
* I could not do a full arXiv listing or full-text search, so novelty is unconfirmed beyond the project's own
  PRIOR-ART search. I did not read Bartholdi 2611.04020 beyond the brief's description (level 2 only).
* Theorem B: I could not decide whether (B3)/(B5) appear in Bressoud–Ismail–Stanton or in Warnaar's papers. They
  follow from (S2′) and (L∞) in three lines, so they are at most a "known-type" lemma; I would not claim them as
  new without a check of BIS §§3–4 and of Warnaar's papers on the A₂/C₂ Bailey lattice.

## Not checked

* §6 items 1, 3 and 4. These are numerical observations with no proof attached and are not used anywhere;
  I did not verify them.
* The certificates in the repository. I did not run them: by design, my checks are independent re-implementations.
* The correctness of the cited theorems of GOW and BW. Corollary C is conditional on them; my evidence for (RS)
  is numerical only (row 8).
