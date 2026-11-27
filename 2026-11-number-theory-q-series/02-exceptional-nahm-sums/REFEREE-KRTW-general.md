# Referee report on proofs/KRTW-general.md (hostile referee, round 4)

Date: 2026-11-26. Object: `proofs/KRTW-general.md` (Theorems 1–3, Lemmas 1–4, Proposition 6, Corollaries 4–5,
conjecture (F_k) of §8). Claim under review: KRTW identity (i),

    Σ_{r,s} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}}) = (q,q^{k+1},q^{k+2};q^{k+2})_∞ / ((q;q²)_∞ (q;q)_∞),

for all k ≥ 1.

Referee code (all written from scratch; the authors' certificate was read but not used):
`/tmp/claude-0/-home-user-weighted-prime-races/3c4df68b-d770-5604-8028-1966b4f5b18d/scratchpad/referee4/code/`
(`rser.py` own truncated-series library; `lemmas_exact.py`, `thm1.py`, `thm2.py`, `thm3.py`, `cor5.py`, `hl.py`,
`sec8.py`, `ek2.py`, `cross.gp`, `cross2.gp`, `eqnum.py`, `axparse.py`; `lemmas.gp` is an abandoned first attempt
and is not used). Downloads: `…/scratchpad/referee4/dl/`.

## Verdict

**Correctness: I was instructed to find the error. I did not find one.** Every step of Lemmas 1–4, Theorem 1,
Theorem 2 (both equalities), Proposition 6, Theorem 3 and Corollaries 4, 5 was re-derived by hand and each
statement was confirmed by independent exact/truncated computation (details in the log below). The proof of
identity (i) for all k ≥ 1 is, as far as I can tell, **correct and complete**: it uses only the q-binomial theorem,
Euler's two series, one limiting case of Bailey's lemma with the unit Bailey pair, and the Jacobi triple product,
and all interchanges of summation are legitimate in Z[x][[q]].

**Completeness:** the main line (Lemma 2 → Lemma 3 → Lemma 4 → Theorem 1 → R→∞ → Proposition 6 → x=1) has no gap.
The gaps I found are all in side remarks (defects 3, 4 below) or in the literature discussion (defects 1, 2).

**Novelty risk: moderate-to-low, but not zero, and the document under-cites.** In KRTW (arXiv v2 = Selecta 2026)
identity (i) is explicitly conjectural for general k (displayed in the Introduction as one of "a number of new
conjectures of Rogers–Ramanujan type" and again after Conjecture 4.9); proved there: k = 1 (Conjecture 2.2 for
k = 1) and k = 2 (Theorem 4.10 + (4.7) + Meurman–Primc). I found no later proof (arXiv API, Semantic Scholar
citations, full text of the most likely candidates; see §7 of the log). What is genuinely new is only the *first*
equality of Theorem 2 (KRTW-type 2k-fold sum = k-fold "Melzer-type" sum, with the grading x^{Σ(r_i+2s_i)}) together
with its finite form Theorem 1 / Lemma 4. The second half (k-fold sum = product) is **in the literature for all
k, odd and even** — contrary to what §5 of the document says (defect 1). Not checked: MathSciNet/zbMATH, Google
Scholar, non-arXiv literature, full texts of most 2024–2026 papers.

## Numbered list of defects (none affects the truth of Theorems 1–3)

1. **[Wrong/incomplete literature statement — moderate]** §5, "For **odd** k the statement is not in the texts I
   read … so it is certainly 'known to the experts'". The evaluation
   Σ q^{Σ M_i(M_i+1)/2}/∏(q)_{M_i−M_{i+1}} = (q,q^{k+1},q^{k+2};q^{k+2})_∞/((q;q²)_∞(q)_∞) for **every** k is stated in
   Warnaar, arXiv:2511.17034 (SIGMA 22 (2026) 062), display following Corollary 5 (the i = 0 case of the identity with
   exponent Σ binom(n_j+1,2) − n_2 − n_4 − … − n_{2i}, product (q^{i+1},q^{k−i+1},q^{k+2};q^{k+2})_∞/((q;q²)_∞(q)_∞)), and is
   attributed there to Andrews, *Parity in partition identities*, Ramanujan J. 23 (2010), Thm 3; Kim–Yee, JCTA 120
   (2013), Thm 1.4; Chern–Li–Stanton–Xue–Yee, J. Algebraic Combin. 60 (2024), Thm 3.5; Kanade–Russell,
   arXiv:2508.15113, Thm 13. KRTW themselves cite Andrews10 and KY13 in the very sentence of their Introduction
   that displays (i) ("very different from those in [Andrews10, CSXY22, Feigin09, KY13, Stembridge90, Tsuchioka23]").
   The document lists arXiv:2511.17034 in §9 as "checked" without noticing this. Consequence: Theorem 3 =
   (first equality of Theorem 2 at x = 1) + a known identity; Proposition 6 is a correct re-proof (with an
   x-refinement), not a new evaluation, and §5 must cite these sources. (I read Warnaar's TeX source; I did not
   read the original Andrews / Kim–Yee / CLSXY / Kanade–Russell papers. The x-refined form (5.1) may well be in
   Kanade–Lovejoy, arXiv:2410.15408, or Dousse–Jang, arXiv:2603.05300 — not checked.)

2. **[Unsupported statement — minor]** Header, "Known before: k = 1 (all four)". No reference is given. From
   KRTW.tex one gets (i) for k = 1 (Conjecture 2.2 holds for k = 1) and (ii) = (iv) for k = 1 (Proposition 4.8,
   eq. (4.12b) for k = 1, with Conjecture 2.2 for k = 1 and the level-one case of Primc's conjecture); for (iii) at k = 1 (the weight Λ_1 line of
   Conjecture 4.9) I found no statement in KRTW that it is proved. Either give a reference or weaken the sentence.
   (Numerically all four hold for k = 1; my check: mod q^40.)

3. **[Gap in a side statement — minor]** §8, "*Proved:* J_k(x) = C_k(x)/C_1(x). (Proposition 6 with (q, x) replaced
   by (q², x/q) …)". The substitution x ↦ x/q is not among the continuous substitutions set up in §1, and under it
   the prefactor 1/(x²q;q)_∞ of (5.1) becomes 1/(x²;q²)_∞ ∉ Z[x][[q]] (it is only after cancelling 1 − x² against
   the summand that one lands in the ring). The statement is true (my check: k ≤ 4 mod q^40, all powers of x) and
   has a one-line rigorous proof — rerun Proposition 6 verbatim in base q² with a = x² (√(aq²) = xq,
   a^{r/2}q^{r²} = x^r q^{r²}, (a;q²)_n/(1−a) = (x²q²;q²)_{n−1}) — but that is not what is written.

4. **[Gap in a remark — minor]** Remark (a) after Theorem 2: "This proves the conjecture (E_k) … in both of its
   forms (… equivalent … as explained there)". `KRTW-k2.md` §4 only has a one-phrase parenthesis. The direction
   needed here (first form ⇒ second form) requires the injectivity of the operator
   F(x) ↦ Σ_b x^{2bk} q^{k(2b²+b)} F(xq^{2b})/(q²;q²)_b on Z[x][[q]] (true: it is F + terms of higher q-order), which is
   stated in neither file. The second form is true numerically (`ek2.py`, k ≤ 5 mod q^40).

5. **[Imprecision — trivial]** §5: "(UBP), (S2) … are identities of rational functions of q and √a". (S2) contains
   a^{k/2} q^{k²/2}, so the field is Q(q^{1/2}, a^{1/2}) (or Q(q, √(aq)) after writing a^{k/2}q^{k²/2} = (√(aq))^k q^{k(k−1)/2}).
   With a = x²q, √(aq) = xq everything is in Z[x][[q]], as the document then says. Also not said: 1/(−xq;q)_n ∈ 𝓡
   (needed for β^{(i)}); obvious.

6. **[Citation precision — trivial]** §5: "Melzer …, (2.6)/(2.3), the character χ̂^{(2,4κ)}_{1,2}". The relevant
   formula is (2.6) with s = 2 (second fermionic form; M̃ = S of BIS). (2.3) is a different fermionic form of the
   same character (with (−q)_{N_1} / a q-binomial). BIS Theorem 5.2 (i = 1, a = x², k_BIS = k/2+1, third line) does
   coincide with (5.1) for even k — I checked every exponent against `papersF/math_9909053.txt`.

7. **[Hypothesis — trivial]** Corollary 4 is stated for k ≥ 0, but Griffin–Ono–Warnaar Theorem 1.2 requires m ≥ 1;
   k = 0 is the trivial identity 1 = 1 and should be said separately.

8. **[Notation — trivial]** §8 uses K_k for two different functions (K_k(X,c) and K_k(x) = K_k(xq^{−1/2}, q^{1/2})),
   and B_j(R;x) (Theorem 1) vs B_{k+1}(x), B^{(a)}_{k+1}(x) (Bressoud-type sums). "Theorem 2 evaluates K_k(X,1), the
   line c = 1 (z = 1/q)": correct only in the sense that this line meets KRTW's two curves X = c and X = c/q at
   z = 1/q; say so.

9. **[Wording — trivial]** Corollary 4, caveat paragraph: "The Hall–Littlewood side was not computed numerically
   here." It can be (and now has been, by me): see log §5. The "caveat" about reading GOW (1.9) from a text
   extraction can be removed: my reading of the same text agrees, and the n = 1 case is confirmed numerically.

10. **[Not checked by me]** Andrews' equation numbers (3.3.6), (2.2.5), (2.2.6) (no copy of the book; they agree
    with my memory); Bressoud's original papers (the document copies KRTW's citation Bressoud79/Bressoud80; the
    three instances used are confirmed numerically, k ≤ 4 mod q^40); the Meurman–Primc and Bartholdi inputs for
    the "known before" k = 2 statements (Bartholdi's Theorem B, boundary (2,0,0), does give the product of (ii)
    at k = 2: 1/(q²,q³,q⁴,q⁴,q⁵,q⁶;q⁸)_∞ — checked by hand).

No defect of severity "error in a proof" or "false statement" was found. In particular nothing conjectural in §8
is presented as proved: (F_k), the product formula in remark X1 and identities (ii)–(iv) are all explicitly
labelled unproved.

## Log of checks

### 1. Lemmas 1–4 (re-derived by hand; exact machine check)

* Hand re-derivation of every line: identity (0.1); Lemma 1 (q-binomial theorem
  (z;q)_n = Σ_i [n;i](−1)^i q^{i(i−1)/2} z^i, valid as a polynomial identity so z may be q^{d+1−t} with negative exponent;
  Euler: 1/(z;q)_∞ = Σ z^n/(q)_n, (z;q)_∞ = Σ(−1)^n q^{n(n−1)/2} z^n/(q)_n); Lemma 2 (coefficient of z^L in
  (−z;q)_R/(z²;q²)_∞ = 1/((z;q)_∞(−zq^R;q)_∞) — both expansions and the cancellation (−z;q)_R/(−z;q)_∞ = 1/(−zq^R;q)_∞);
  Lemma 3 (t = r+i, [R;r][R−r;i] = [R;t][t;i], exponent t² + i(i−1)/2 + i(d+1−t), (q^{d+1−t};q)_t = 0 for t > d and
  = [d;t](q)_t for t ≤ d); Lemma 4 (∗) (left side via Lemma 3 with y = xq^{2b}, d = M−2b, then Lemma 2; right side via
  Lemma 1 and (L+i, i) → (L, j)); (∗∗) (multiply by (xq^{R+M+1};q)_∞); (∗∞) (Lemma 2 and Euler). No error.
* `lemmas_exact.py` (python-flint, denominators cleared, exact in Z[q][x]): Lemma 2 for R ≤ 8, L ≤ 12 (117 cases),
  Lemma 3 for R ≤ 7, d ≤ 8 (72 cases), Lemma 4 (∗) for R, M ≤ 8 (81 cases): all hold exactly. A deliberately wrong
  variant of (∗) (shift R+L instead of R+L+1) is rejected in all 72 cases with M ≥ 1.

### 2. Theorem 1 (finite R-refinement)

* Statement precise (ranges R ≥ r_1 ≥ … ≥ r_j ≥ 0, 0 ≤ s_1 ≤ … ≤ s_j; r_{j+1} = s_0 = M_{j+1} = 0; factor (q)_R/(q)_{R−r_1}).
* Induction step re-derived in full: (q)_R/((q)_{R−r_1}(q)_{r_1−r_2}) = [R;r_1](q)_{r_1}/(q)_{r_1−r_2}; shift s_i = b + s'_i giving
  x^{2b(j−1)} q^{(j−1)(2b²+b)} A_{j−1}(r; xq^{2b}) (the factor (q)_r/(q)_{r−r_2} and the constraint r_2 ≤ r are exactly those of
  A_{j−1}(r;·)); induction hypothesis under the continuous substitution x → xq^{2b}; shift M'_i = M_i + 2b, range
  0 ≤ b ≤ ⌊M'_{j−1}/2⌋, denominator (q)_{M'_{j−1}−2b}; (∗∗) applied with M = M'_{j−1} and the same R (hypotheses R, M ≥ 0
  met). Interchange of Σ_{r,b} with Σ_{M'}: the family indexed by (r, b, M') has q-order ≥ (r+2b)(r+2b+1)/2 +
  Σ M'_i(M'_i+1)/2 with r ≤ R, hence is summable in Z[x][[q]]. Base case j = 1 is (∗∞). No error.
* `thm1.py` (brute force over (r,s) straight from the definition): A_j(R;x) = B_j(R;x) for j ≤ 4, R ≤ 5 mod q^40 and
  for j ≤ 5, R ≤ 8 mod q^60, in all powers of x: all cases OK. Negative control (A_j(R) against B_j(R+1)) fails at
  q^{R+1} as it must.

### 3. Theorem 2 and Proposition 6

* R → ∞: a term has q-order ≥ r_1², (q)_R/(q)_{R−r_1} ≡ 1 mod q^{R−r_1+1}, r_1² + R − r_1 + 1 ≥ R+1, r_1² ≥ R+1 for r_1 > R;
  (xq^{R+M_k+1})_∞ ≡ 1 mod q^{R+1}. Fine. Identification of the limits with KRTW_k(x), Mel_k(x): fine.
* Bailey step: (S2) re-derived from the general Bailey lemma (ρ_1 → ∞, ρ_2 = −√(aq)) and compared with BIS
  (`papersF/math_9909053.txt`: definition of Bailey pair "with parameters (a,q)", Bailey's Lemma at the end of §1,
  (UBP) and (S2) in §3) — identical. a = x²q, √(aq) = xq, a^{r/2}q^{r²/2} = x^r q^{r(r+1)/2}; β^{(1)}_n = 1/((q)_n(−xq)_n);
  induction for β^{(i)}; n → ∞ with the stated error orders; (x²q²;q)_∞(1−x²q) = (x²q;q)_∞; final exponent
  (k+1)n(n+1)/2 + n(n−1)/2. Fine.
* `thm2.py`: KRTW_k(x) = Mel_k(x) = single sum (5.1) for k ≤ 5 mod q^60 and k ≤ 6 mod q^90, all powers of x (x-degree
  up to 29): OK. Negative control (KRTW's grading z^{r+s}, i.e. x^{2(r+s)}) fails at q^1. Independent PARI/GP
  cross-check `cross2.gp`: first equality for k = 2, 3 and (5.1) for k = 3, mod q^45: OK.

### 4. Theorem 3 (identity (i))

* x = 1: (x²q;q)_n/(q)_n = 1; exponent ((k+2)n²+kn)/2; n ↦ −n−1 gives the terms −q^{2n+1}; Jacobi triple product with
  Q = q^{k+2}, z = q^{k+1} gives (q,q^{k+1},q^{k+2};q^{k+2})_∞; (−q)_∞ = 1/(q;q²)_∞. Correct for every k ≥ 0; parity of k
  plays no role.
* Left side compared with KRTW.tex (lines 2469–2475 and 632–638): exponent, denominators, r_{k+1} = s_0 := 0, sum
  over all r_i, s_i ≥ 0 with 1/(q)_{n<0} = 0. Identical.
* `thm3.py`: literal KRTW display (all r_i, s_i ≥ 0, vanishing reciprocal Pochhammers) k = 1, 2, 3 mod q^20;
  chain-enumerated brute force: k = 1, 2, 3, 4 mod q^500; k = 5 mod q^300. All OK. PARI cross-check `cross.gp`: k = 1, 2, 3 mod q^45: OK.

### 5. Corollaries 4, 5

* Corollary 4. GOW Theorem 1.2, eq. (1.9) (`papersH/1401.7718.txt`): for positive integers m, n, κ = 2m+2n+2,
  Σ_{λ_1≤m} q^{|λ|}P_{2λ}(1,q,q²,…;q^{2n}) = (q²;q²)_∞(q^{κ/2};q^{κ/2})_∞(q^κ;q^κ)^{n−1}_∞/(q)^{n+1}_∞ · ∏_{i=1}^n θ(q^i;q^{κ/2})
  ∏_{i<j} θ(q^{j−i},q^{i+j};q^κ), θ(a;q) = (a;q)_∞(q/a;q)_∞. At n = 1, m = k this is the product in the document and equals
  the product of (i). The deduction is airtight given GOW's published theorem; it needs nothing else from GOW. It
  is also exactly how KRTW say (i) arises (their Proposition 4.8: their (4.12a) "holds for z = 1").
  Independent numerical confirmation, `hl.py`: P_{2λ}(1,q,q²,…;q²) computed from Macdonald's tableau formula
  (III.5.11′, ψ_{λ/μ}); for k = 2, 3, 4 mod q^26: HL sum = GOW product = multisum (i). (The z-graded Conjecture 2.2
  also holds to that order; the document correctly does not claim it.)
* "Conjecture 4.9 (first line) holds at z = 1/q": uses (4.7) (C^{(1)}_{kΛ_0}(z) = D^{(2)}_{kΛ_0}(z/q)) and Conjecture 4.1 for
  n = 1, which KRTW attribute to Meurman–Primc [MP99]; the product of Conjecture 4.1 at n = 1, k_0 = k is the
  product of (i) (checked). Valid.
* Corollary 5: proof re-derived (Σ_{R≥r_1} y^R/(q)_{R−r_1} = y^{r_1}/(y)_∞; Euler expansion of (xq^{R+M_k+1})_∞;
  Σ_R y^R q^{Rl}/(q)_R = 1/(yq^l)_∞; (y)_∞/(yq^l)_∞ = (y)_l; polynomiality in y so that y = 1 is allowed). Correct.
  `cor5.py`: k ≤ 4, y ∈ {0, 1, −1, 2, −3, q, q², −q, 5q³}, mod q^40, all powers of x: OK.

### 6. Section 8

* Status labels are honest: (F_k), X1's product formula, (ii)–(iv) are marked unproved/numerical.
* Hand checks: the dictionary K_k(X,c) ↔ KRTW's lines (c² = zq, X = c resp. X = c/q; Conjecture 2.2 ↔ X = c = √z);
  (ii), (iii) at c = q^{1/2}, X = q^{±1/2}; both factorisations of the product of (iii) (the second uses the display
  before Theorem 4.1 of Warnaar math/0110307, which I read: (q^{k/2},q^{k/2+1},q^{k+1};q^{k+1}) for k−1 variables);
  J_k(1), J_k(q) by the triple product; (ii) = B_{k+1}(q)J_k(q), (iii) = B_{k+1}(1)J_k(1), (iv) = B^{(k)}_{k+1}(1)J_k(1) at
  the level of products; KRTW (3.6) is F^{(n)}_{a,δ} as quoted. All correct. (J_k = C_k/C_1: true, but see defect 3.)
* `sec8.py`, k ≤ 4 mod q^40: C_1 = (−xq;q²)_∞, J_kC_1 = C_k, J_k(1), J_k(q), three Bressoud instances, both statements
  of (F_k) (all powers of x), (ii), (iii), (iv), second form of (iii): all OK. Remark X1 (k ≤ 6; after x → 1/q a
  truncation mod q^120 is valid only mod q^85…q^105, which is what was compared): OK.

### 7. Novelty / openness

* (a) KRTW.tex (v2, journal-ref Selecta Math. 32 (2026) no. 5, Paper No. 100, confirmed via the arXiv API).
  Numbering reconstructed from the TeX source: Conjecture 2.2 = `Con_Shun`, eq. (2.12); Proposition 4.8; Conjecture
  4.9 = `Con_Shun2`; Theorem 4.10 = `Thm_kis2`; (4.7) = `Eq_CDn2`; (3.6) = `Eq_Fdef` — all as cited in the document.
  KRTW prove: Conjecture 2.2 for k = 0, 1; Proposition 4.8 ((4.12a) for n = 0, (4.12b) for n = 1, both for k = 1,
  (4.12a) at z = 1 via GOW + Primc–Trupčević; equation numbers reconstructed by counting in the TeX source); Theorem 4.10 (k = 2, four weights, bivariate, by functional
  equations). The four displayed identities are called "conjectural Andrews–Gordon-type identities"; (i) and
  (iii) are displayed as new conjectures in the Introduction. So (i) is conjectural there for k ≥ 3. Confirmed.
* (b) Searches (2026-11-26). arXiv API, sorted by date: all:Capparelli AND Meurman AND Primc (10 hits);
  abs:"Hall-Littlewood" AND abs:"Rogers-Ramanujan" (12); abs:"Andrews-Gordon" 2024-04…2026-11 (22);
  abs:"Hall-Littlewood" 2025–2026 in math.CO/NT/QA/RT (34); abs:CMPP / "Capparelli, Meurman, Primc" (9); Primc
  conjecture variants (5); multisum AND Rogers-Ramanujan (4); (Rogers-Ramanujan OR Andrews-Gordon OR Bailey pair
  OR Bailey lemma OR Nahm sum) 2024-04…2026-11 (130 titles scanned); authors Warnaar, Kanade, M. C. Russell,
  Tsuchioka, Bartholdi, Primc, Trupčević, Kurşungöz, Liuquan Wang, Xiangyu Ding. Semantic Scholar citations of
  arXiv:2404.03851: 8 papers (Ding 2611.10654; Bartholdi 2611.04020; Bhattacharya–Dai–Hou–Tong 2610.41155; Lau–Ono
  2608.05480; Šikić 2606.19994, 2404.01791; Warnaar 2511.17034; Russell 2306.16251). Full text/source read or
  grepped: Warnaar 2511.17034 (cites KRTW three times, none about Conjecture 2.2/4.9; contains defect 1's
  identity; its mixed-base identities (A.1) are of a different shape), Bhattacharya et al. (uses only KRTW's
  formula for P_{(2^r,1^s)}), Bartholdi (three-row level-two Primc), Ding (C_n^{(1)} products, no multisum),
  Kurşungöz 2605.13096 (k = 1). Abstracts read: Fischer–Gangl–Gutiérrez–Kumari 2610.08500, Takigiku–Tsuchioka
  2606.31821, Kılıç–Kurşungöz 2501.03189, Butorac–Kožić–Primc 2603.23120, CMP 2611.32425, Haijun Li 2610.35234,
  Berkovich–Dhar 2609.07490, Hakimoglu-Brown 2609.14940, Kanade–Lovejoy 2410.15408, Kanade–Russell 2508.15113,
  Dousse–Jang 2603.05300 and others. **No proof of (i) for k ≥ 3, and no statement of Theorem 2/Theorem 1/Lemma 4,
  was found.**
* (c) Is Theorem 2 or Lemma 4 already in BIS / Warnaar / Stembridge-type papers? I read BIS §§1–5, 7 and Warnaar
  math/0110307 §4 in the text files: the multisums there are single chains (possibly with changing base), not the
  coupled two-chain shape (q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}}, and nothing resembling (∗)/(∗∗) with a boundary
  parameter R appears. The KRTW multisum itself dates from 2024, so an earlier occurrence of Theorem 2 in this
  form is unlikely; Lemma 3 (a terminating q-Chu–Vandermonde-type polynomial identity) and Lemma 2 are surely
  classical in some form. **Could not check:** Stembridge 1990, Jouhet–Zeng, Ishikawa–Jouhet–Zeng, Andrews 2010,
  Kim–Yee 2013, CLSXY 2024, Kanade–Lovejoy, Dousse–Jang (not available here), MathSciNet/zbMATH, the printed
  Selecta version of KRTW (a note added in proof is not in the arXiv v2 source).

### 8. Things the authors' certificate does that I did not repeat, and vice versa

Not repeated: k = 6, 7, 8 of (i) by brute force; (F_k) for k = 5, 6, 7; Mel_k(1) for k up to 20. Added by me: exact
Lemma 4 up to R, M = 8; Theorem 1 up to j = 5, R = 8; the Hall–Littlewood side of Corollary 4; the second form of
(E_k); negative controls for each test; a PARI cross-check of my own Python library.
