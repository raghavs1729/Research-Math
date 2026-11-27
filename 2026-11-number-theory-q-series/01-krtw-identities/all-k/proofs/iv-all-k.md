# The KRTW identity (iv) for all k

Status (2026-11-27, krtw4 session): **complete proof written; certificate passes; two independent adversarial reads
done** (`../REFEREE-iv-all-k.md`, own code in `../referee2/`; `../REFEREE-3-all-k.md`, own code in `../referee3/`, written
without seeing the first report or any code): no error, no gap; their exposition items are incorporated below.
Certificate: `code/certify_iv_all_k.py` → `certificates/certify_iv_all_k.txt` (labels [D1]–[D8] below; ALL CHECKS
PASSED with q^40, k ≤ 8). Independent series check of (iv) itself mod q^300 for k ≤ 8: `code/series_check_q300.py`
→ `certificates/series_check_q300.txt` (separate code, Ω expanded literally).

This file builds on `F-all-k.md` (same directory; notation, §§1–7 and Theorem 6.4 used without change) and on
`../../rr-identities/proofs/KRTW-F-moves.md` (Theorem B (B2) and the evaluation of J_k(1) in its §3).

## 0. Statement

Source: Kanade–Russell–Tsuchioka–Warnaar, arXiv:2404.03851v2, the fourth display after Conjecture 4.9 (verbatim LaTeX
in `../LITERATURE-2026-11-27.md` §1b). With r_{k+1} = s_0 := 0, the sum over all r_1, …, r_k, s_1, …, s_k ≥ 0 and
1/(q;q)_m = 1/(q²;q²)_m := 0 for m < 0:

    (iv)   Σ Ω^{s_1,…,s_k}_{r_1,…,r_k}(q) ∏_{i=1}^k q^{(r_i+s_i)²+s_i²} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}})
             = (q^k, q^{k+1}, q^{k+3}, q^{k+4}, q^{2k+4}, q^{2k+4}; q^{2k+4})_∞ / ((q²;q²)_∞ (q;q)_∞),

    Ω^{s_1,…,s_k}_{r_1,…,r_k}(q) := Σ_{i=1}^{k−1} q^{r_i+2s_i}(1 − q^{2s_{i+1}−2s_i}) + q^{r_k+2s_k}.

Notation as in `F-all-k.md`: Q := q², λ_{N,M} = x^{N+M}q^{N²+M²}, λ_n := x^nq^{n²}, ε_{N,M} = q^{−N−M}(1+xq^{2N})(1+xq^{2M}),
β^B_n, β^J_m, t_{N,M} as in its §7, K_k, J_k as in its §0. Put

    K^Ω_k(x) := Σ_{r_1≥…≥r_k≥0, 0≤s_1≤…≤s_k} Ω^{s}_{r}(q) ∏_{i=1}^k x^{r_i+2s_i} q^{(r_i+s_i)²+s_i²} / ((q)_{r_i−r_{i+1}} (Q;Q)_{s_i−s_{i−1}}),
    B^{(k)}_{k+1}(x) := Σ_{N_1≥…≥N_{k+1}≥0} x^{ΣN_i} q^{ΣN_i² + N_k + N_{k+1}} / ((q)_{N_1−N_2}⋯(q)_{N_k−N_{k+1}} (Q;Q)_{N_{k+1}}),
    b^Ω_n := β^B_n (q^{−2n} + x − 1 + x²q^{2n}) / (x(1+x)),

so b^Ω_0 = 1 and, for n ≥ 1, b^Ω_n = (−1)^n q^{n²−2n} (x²Q;Q)_{n−1} (1−xq^{2n}) (1 − (1−x)q^{2n} + x²q^{4n}) / (x (Q;Q)_n)
(use (x²;Q)_n/(1−x) = (1+x)(x²Q;Q)_{n−1}). The left side of (iv) is K^Ω_k(1) (restricting to r decreasing and s
increasing loses nothing: every other term contains a factor 1/(q)_m or 1/(Q;Q)_m with m < 0).

**Theorem 4.** For every k ≥ 1, in Z[x][[q]],

    (xq;q)_∞ K^Ω_k(x) = J_k(x) · Σ_{n≥0} λ_n^{k+1} b^Ω_n(x).

**Theorem 5.** The KRTW identity (iv) holds for every k ≥ 1.

**Theorem 6 (x-refinement).** For every k ≥ 1, in Z[x][[q]], (xq;q)_∞ B^{(k)}_{k+1}(x) = Σ_n λ_n^{k+1} b^Ω_n(x); hence
K^Ω_k(x) = B^{(k)}_{k+1}(x) J_k(x), the second identity of the conjectural factorisation (F_k) of
`../../rr-identities/proofs/KRTW-general.md` §8 (its first identity is Theorem 1 of `F-all-k.md`).

What is new: Theorem 5 for k ≥ 3 (for k = 1 the left side and the right side of (iv) coincide with those of (ii) for
k = 1, a known case; k = 2 is covered, like (ii) and (iii), by KRTW's Theorem 4.10 together with the level-two
three-row case of Primc's conjecture, Bartholdi arXiv:2611.04020, Theorem A — this last attribution is taken from the
project's earlier literature record `../../rr-identities/PRIOR-ART.md` and was not re-derived; the proof below covers
k = 2 anyway), and Theorems 4 and 6 for k ≥ 2 (for k = 1, K^Ω_1(x) = K_1(xq), B^{(1)}_2(x) = B_2(xq) and J_1 = 1, so
Theorems 4 and 6 follow from (F_1) and from (B1) at x ↦ xq; no x-refinement of (iv) is in the literature as far as our
searches go: `../LITERATURE-2026-11-27.md`, `../NOVELTY-2026-11-27.md`).

**Idea.** K^Ω_k inserts the operator O = PY𝔲 once into the transfer-operator chain, at every level. O is a commutator
with T up to a factor: [T, h] = −xOT for the rank-one operator h = h_{xY} that already appears in H′ (Lemma 2.1). So the
k insertions telescope to a boundary term (Proposition 2.3): K^Ω_k = x^{−1}(K^♭_k − K_k), where K^♭_k is K_k with one extra
factor q^{−r_1}. On each eigenvector, ω∘h is ε/(1+x) times ω (Lemma 3.1) — this is the first slice equation of H′φ = εφ
read at infinity, already used in `F-all-k.md` §4. So the spectral weights of K^Ω_k are t_{N,M}(ε_{N,M}−1−x)/(x(1+x)), and
an elementary factorisation (Lemma 4.1) splits them into a product of a B-type and a J-type one-dimensional sum.

## 1. Operator form

Setting of `F-all-k.md` §1: 𝔽 = Q(q)((x)), the space 𝒱, T, u, ω, P, Y, 𝖲, 𝔲, W, E. Put O := PY𝔲, i.e.
(OF)(R,S) = q^{R+2S}F(R,S+1); O𝒱 ⊂ 𝒱.

**Lemma 1.1.** K^Ω_k(x) = Σ_{i=1}^{k} ω(T^{i−1} O T^{k−i} u)  (k ≥ 1).

*Proof.* As for (1.1) of `F-all-k.md`: with F_k := u and F_i := TF_{i+1}, F_i is the function of (r_i, s_i) obtained by
summing levels i+1, …, k (the weights w(r_j,s_j), j > i, and the bonds between levels i, …, k); the bond between levels
i and i+1 is 1/((q)_{r_i−r_{i+1}}(Q;Q)_{s_{i+1}−s_i}), with (R,S) = (r_i,s_i), (r,s) = (r_{i+1},s_{i+1}) in the definition of T, and
ω supplies w(r_1,s_1)/(Q;Q)_{s_1}. Expand Ω into its k terms. For i < k the i-th term multiplies the summand by
q^{r_i+2s_i}(1−Q^{s_{i+1}−s_i}); since (1−Q^m)/(Q;Q)_m = 1/(Q;Q)_{m−1} for m ≥ 0 (both sides 0 for m = 0, and both bond
factors vanish for m < 0), this replaces the bond by 1/((q)_{r_i−r_{i+1}}(Q;Q)_{s_{i+1}−s_i−1}) and multiplies by q^{r_i+2s_i}, i.e.
replaces F_i = TF_{i+1} by (PY𝔲TF_{i+1}) = OTF_{i+1}, because Σ_{r≤R, s≥S+1} w(r,s)F(r,s)/((q)_{R−r}(Q;Q)_{s−S−1}) = (TF)(R,S+1).
The k-th term multiplies F_k = u by q^{r_k+2s_k}, giving PYu = Ou (u does not depend on S). For each power of x only
finitely many index tuples (r, s) contribute (the exponent of x is Σ(r_i+2s_i)), so splitting Ω into its terms and
regrouping the 2k-fold sum into the nested sums T^{i−1}OT^{k−i}u are finite rearrangements, coefficient by coefficient. ∎

## 2. A commutator, and the telescoping

Let h := h_{xY} = P^{−1}(1−𝖲) + xPY, i.e. (hF)(R,S) = q^{−R}(F(R,S) − F(R−1,S)) + xq^{R+2S}F(R,S) (with F(−1,S) = 0); h𝒱 ⊂ 𝒱.

**Lemma 2.1.** Th − hT = −xOT on 𝒱.

*Proof.* Normal ordering exactly as in the proof of Proposition 2.1 of `F-all-k.md`. By (1.3) there, W𝖲W^{−1} = xq^{−1}P²Y𝖲,
so Wh = h̃W with h̃ := P^{−1} − xq^{−1}PY𝖲 + xPY. It suffices to show hE − Eh̃ = xPY𝔲E: then
hT − Th = hEW − Eh̃W = xPY𝔲EW = xOT. Coefficients of 𝖲^j𝔲^b, multiplied by (q)_j(Q;Q)_b (using 𝖲^jP^{−1} = q^jP^{−1}𝖲^j,
𝖲^j𝔲^b·PY = q^{−j}Q^bPY𝖲^j𝔲^b, 1/(q)_j − 1/(q)_{j−1} = q^j/(q)_j, (q)_j/(q)_{j−1} = 1−q^j, (Q;Q)_b/(Q;Q)_{b−1} = 1−Q^b):

    in hE:          q^jP^{−1} + xPY,
    in Eh̃:          q^jP^{−1} − xq^{−j}Q^b(1−q^j)PY + xq^{−j}Q^bPY = q^jP^{−1} + xQ^bPY,
    in xPY𝔲E:       xPY(1−Q^b),

and (q^jP^{−1} + xPY) − (q^jP^{−1} + xQ^bPY) = xPY(1−Q^b). Restriction to N² (as in `F-all-k.md` §2): h maps functions
on N² (extended by 0) to such functions and reads only (R,S) and (R−1,S); xPY𝔲 reads EWF only at (R,S+1) ∈ N² when
(R,S) ∈ N²; so the identity on Z², restricted to N², is the claim. ∎ [D1: exact, on point masses, two rational points.]

**Lemma 2.2.** ω(hF) = ω(P^{−1}F) for F ∈ 𝒱.

*Proof.* Fix s. Summation by parts (F(−1,s) = 0; all sums converge x-adically since ord w(r,s) = r+2s):
Σ_r w(r,s)q^{−r}(F(r,s) − F(r−1,s)) = Σ_r F(r,s)(w(r,s)q^{−r} − w(r+1,s)q^{−r−1}) = Σ_r w(r,s)F(r,s)(q^{−r} − xq^{r+2s}), because
w(r+1,s) = xq^{2r+2s+1}w(r,s). Adding Σ_r w(r,s)·xq^{r+2s}F(r,s) gives Σ_r w(r,s)q^{−r}F(r,s); divide by (Q;Q)_s and sum over s. ∎ [D2]

**Proposition 2.3.** For k ≥ 1, with K^♭_k(x) := ω(P^{−1}T^{k−1}u) (= K_k(x) with the extra factor q^{−r_1} in the summand),

    K^Ω_k(x) = x^{−1} (K^♭_k(x) − K_k(x)).

*Proof.* hu = u + xOu: indeed q^{−R}(1/(q)_R − 1/(q)_{R−1}) = 1/(q)_R (also at R = 0), and xq^{R+2S}/(q)_R = x(Ou)(R,S).
For i < k, Lemma 2.1 gives OT^{k−i}u = −x^{−1}(Th − hT)T^{k−i−1}u. With a_i := ω(T^ihT^{k−1−i}u) (0 ≤ i ≤ k−1) the first k−1 terms of
Lemma 1.1 add up to −x^{−1}Σ_{i=1}^{k−1}(a_i − a_{i−1}) = −x^{−1}(a_{k−1} − a_0), and
−x^{−1}a_{k−1} = −x^{−1}ω(T^{k−1}u) − ω(T^{k−1}Ou) = −x^{−1}K_k − (the k-th term of Lemma 1.1). Hence K^Ω_k = x^{−1}(a_0 − K_k), and
a_0 = ω(hT^{k−1}u) = ω(P^{−1}T^{k−1}u) by Lemma 2.2. For k = 1 the statement is ω(Ou) = x^{−1}(ω(hu) − ω(u)), the same computation. ∎
[D3: both sides by direct summation, k ≤ 8, all powers of x, mod q^40; K^Ω_k summed with the 2k−1 signed insertions of Ω
taken literally; a negative control is detected.]

(Proposition 2.3 is an identity between multisums; with Theorem 1 of `F-all-k.md` and Theorem 6 below it says
K^♭_k = (B_{k+1} + xB^{(k)}_{k+1})J_k.)

## 3. The spectral weights

**Lemma 3.1.** ω(hφ_{N,M}) = (ε_{N,M}/(1+x)) · ω(φ_{N,M}) for all N ≥ M ≥ 0.

*Proof.* (a) h is diagonal on the basis Ξ_{e,σ}: on the slice S = σ it acts as h_{xq^{2σ}} in R, so by (B-iii)
hΞ_{e,σ} = (q^{−e} + xq^{2σ+e})Ξ_{e,σ}. Hence hφ_{N,M} is again a finite combination of Ξ's.

(b) For every finite combination F of Ξ's, (TF)(·,0) is a finite combination of the φ_e(·;x) (by (B-i) and (B-iv)), whose
coefficients are unique since the φ_e(·;x) are linearly independent (φ_e(R;x) = 0 for R < e, ≠ 0 for R = e); and
(xq)_∞ ω(F) = v_x((TF)(·,0)). This is the computation in the proof of Lemma 4.3 of `F-all-k.md`, which uses only
linearity: for F = Ξ_{e,s}, (B-ii) gives (xq)_∞ω(Ξ_{e,s}) = w(e,s)(xq)_{2s}/(Q;Q)_s, while (B-i) gives
(TΞ_{e,s})(·,0) = w(e,s)φ_e(·;xq^{2s})/(Q;Q)_s, whose v_x is the same by (B-iv) (v_x(φ_e(·;xq^{2s})) = (xq)_{2s}).

(c) By Lemma 2.1 and Tφ = λφ (λ := λ_{N,M}): T(hφ) = hTφ − xOTφ = λ(hφ − xOφ). Its slice S = 0 is λ(h_xF_0 − x q^R F_1(R)),
with F_0 := φ(·,0), F_1 := φ(·,1) as in the proof of Proposition 4.1 of `F-all-k.md`. For M ≥ 1 that proof shows
v_x(q^RF_1) = 0 and (1+x)v_x(h_xF_0) = ε v_x(F_0) (the first slice equation of H′φ = εφ, read at infinity). For M = 0,
F_1 = 0 and φ_{L,0} = Ξ_{L,0}, so h_xF_0 = (q^{−L}+xq^L)F_0 and ε_{L,0} = (1+x)(q^{−L}+xq^L). In both cases
v_x((Thφ)(·,0)) = λ ε v_x(F_0)/(1+x) = λ ε W_{N,M}/(1+x), and by (b) and Lemma 4.3 of `F-all-k.md`,
ω(hφ) = λεW_{N,M}/((1+x)(xq)_∞) = (ε/(1+x)) ω(φ). ∎
[D4: (1+x)v_x(h_xF_0) = εv_x(F_0) on eigenvectors computed from the triangular system of T, N ≤ 5; D4b: the statement of
the lemma itself, N ≤ 5, by the finite Cauchy sums of (b) applied to hφ and φ, without v_x or the commutator.]

**Theorem 3.2.** For k ≥ 1, in 𝔽:

    (xq)_∞ (x²Q;Q)_∞ K^♭_k(x) = Σ_{N≥M≥0} λ_{N,M}^k t_{N,M} ε_{N,M}/(1+x),
    (xq)_∞ (x²Q;Q)_∞ K^Ω_k(x) = Σ_{N≥M≥0} λ_{N,M}^k t_{N,M} (ε_{N,M} − 1 − x)/(x(1+x)).

*Proof.* By Theorem 6.4 of `F-all-k.md`, T^{k−1}u = Σ ũ_{N,M}λ^{k−1}_{N,M}φ_{N,M}, and the family
(w(r,s)ũ_bλ_b^{k−1}φ_b(r,s)/(Q;Q)_s)_{b,r,s} is x-adically summable (that is how K_k = Σũλ^kW/(xq)_∞ is obtained there).
Multiplying the (r,s)-entry by q^{−r} does not change ords, so K^♭_k = ω(P^{−1}T^{k−1}u) = Σ_b ũ_bλ_b^{k−1}ω(P^{−1}φ_b), and
ω(P^{−1}φ_b) = ω(hφ_b) = ε_bω(φ_b)/(1+x) by Lemmas 2.2 and 3.1. With Lemma 4.3 and Lemma 7.1 of `F-all-k.md`,
ũ_{N,M}λ_{N,M}W_{N,M}/(xq)_∞ = λ_{N,M}t_{N,M}/((xq)_∞(x²Q;Q)_∞) (the same step as in its Theorem 6.4 and proof of Theorem 1), which gives the
first formula; the second follows with Proposition 2.3 and (xq)_∞(x²Q;Q)_∞K_k = Σλ^kt. ∎ [D6: the second formula against
direct summation, k ≤ 8, mod q^40, all powers of x; a perturbed ε is detected.]

## 4. The product side; proof of Theorem 4

**Lemma 4.1.** For all N ≥ M ≥ 0,
t_{N,M}(ε_{N,M} − 1 − x)/(x(1+x)) = λ_Nλ_M (b^Ω_Nβ^J_M + b^Ω_Mβ^J_N)/(1+δ_{NM}).

*Proof.* Put μ_n := q^{−2n} + x²q^{2n}, so b^Ω_n = β^B_n(μ_n + x − 1)/(x(1+x)). Since t_{N,M} = λ_Nλ_M(β^B_Nβ^J_M + β^B_Mβ^J_N)/(1+δ_{NM}), the claim
is equivalent to

    β^B_Nβ^J_M (ε − μ_N − 2x) + β^B_Mβ^J_N (ε − μ_M − 2x) = 0      (ε := ε_{N,M}).

With u := q^N, v := q^M, ε = 1/(uv) + xu/v + xv/u + x²uv, and expanding the right side shows
ε − μ_N − 2x = (u−v)(1+xu²)(1−xuv)/(u²v), ε − μ_M − 2x = (v−u)(1+xv²)(1−xuv)/(v²u). So the left side is
((u−v)(1−xuv)/(uv)) · [β^B_N(1+xu²)u^{−1}·β^J_M − β^B_M(1+xv²)v^{−1}·β^J_N], and β^B_n(1+xq^{2n})q^{−n} = (1+x)β^J_n for all n ≥ 0
(both are (−1)^nq^{n²−n}(x²;Q)_n(1−x²q^{4n})/((Q;Q)_n(1−x))), so the bracket is (1+x)(β^J_Nβ^J_M − β^J_Mβ^J_N) = 0. ∎
[D5: symbolic, N ≤ 12.]

*Proof of Theorem 4.* By Theorem 3.2 and Lemma 4.1, with λ_{N,M} = λ_Nλ_M and symmetrisation (the double series converge
x-adically: ord λ_n^{k+1}b^Ω_n ≥ (k+1)n − 1 for n ≥ 1, and ord λ_m^{k+1}β^J_m ≥ (k+1)m),

    (xq)_∞(x²Q;Q)_∞K^Ω_k = Σ_{N≥M} (λ_Nλ_M)^{k+1}(b^Ω_Nβ^J_M + b^Ω_Mβ^J_N)/(1+δ_{NM}) = (Σ_n λ_n^{k+1}b^Ω_n)(Σ_m λ_m^{k+1}β^J_m),

and Σ_mλ_m^{k+1}β^J_m = (x²Q;Q)_∞J_k(x) by (B2) of `KRTW-F-moves.md`. ((B2) is proved there in Z[[q]][[x]]; it holds in 𝔽
because the x^D-coefficients of both of its sides are finite sums of rational functions of q regular at q = 0, and the
coefficientwise Taylor map from such series to Q[[q]][[x]] is injective.) This is Theorem 4 in 𝔽. Passage to Z[x][[q]]: the coefficient of
each power of x on either side is a finite sum of rational functions of q regular at q = 0 (for K^Ω_k the exponent of x
bounds the indices; J_k = C_k/C_1 with C_1 ≡ 1 mod x; for the n-sum, ord ≥ n), so the identity holds in Q[[q]][[x]];
both sides lie in Z[x][[q]] (for K^Ω_k: each of the 2k−1 monomial insertions of Ω has exponent ≥ 0, and q-exponents grow
with the indices; for the n-sum: the n-th term is x^{(k+1)n−1}q^{(k+2)n²−2n} times a polynomial in x with coefficients in Z[[q]];
J_k ∈ Z[x][[q]] by `KRTW-F-moves.md` §3), and Z[x][[q]] → Q[[q]][[x]] is an injective ring homomorphism. ∎

## 5. Proof of Theorem 5 ((iv) for all k)

x ↦ 1 is a ring homomorphism Z[x][[q]] → Z[[q]]; on K^Ω_k it gives the left side of (iv) (for each power of q only
finitely many terms contribute), and (xq)_∞ ↦ (q)_∞. At x = 1, b^Ω_0 = 1 and, for n ≥ 1 (from the factored form in §0),
b^Ω_n = (−1)^nq^{n²−2n}(Q;Q)_{n−1}(1−q^{2n})(1+q^{4n})/(Q;Q)_n = (−1)^n(q^{n²−2n} + q^{n²+2n}), so

    Σ_n λ_n^{k+1}b^Ω_n |_{x=1} = 1 + Σ_{n≥1}(−1)^n(q^{(k+2)n²−2n} + q^{(k+2)n²+2n}) = Σ_{n∈Z}(−1)^n q^{(k+2)n²−2n} = (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_∞

by the Jacobi triple product Σ_{n∈Z}(−1)^nq^{An²+Bn} = (q^{A−B}, q^{A+B}, q^{2A}; q^{2A})_∞ with A = k+2, B = −2 [D8: mod q^200, k ≤ 10].
By `KRTW-F-moves.md` §3 (from (B2) and the triple product), (q²;q²)_∞J_k(1) = (q^{k+1}, q^{k+3}, q^{2k+4}; q^{2k+4})_∞. Theorem 4 at x = 1 gives

    K^Ω_k(1) = (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_∞ (q^{k+1}, q^{k+3}, q^{2k+4}; q^{2k+4})_∞ / ((q;q)_∞ (q²;q²)_∞),

which is the right side of (iv). ∎

## 6. Proof of Theorem 6 (the x-refinement)

Bailey pairs relative to a (base q): β_n = Σ_{r≤n} α_r φ_r(n;a), φ_r(n;a) = 1/((q)_{n−r}(aq)_{n+r}) as in `F-all-k.md` §3.

1. For any a, β_n = 1/(Q;Q)_n, α_n = (−1)^nq^{n²}(a²;Q)_n(1−aq^{2n})/((Q;Q)_n(1−a)) is a Bailey pair relative to a (the proof of (B1)
   in `KRTW-F-moves.md` §4 is for general a). Take a = xq.
2. Two steps (L∞) relative to xq ((B-i)): α″_n = (xq)^{2n}q^{2n²}α_n, β″_n = Σ_{n≥N_k≥N_{k+1}} (xq)^{N_k+N_{k+1}}q^{N_k²+N_{k+1}²}/((q)_{n−N_k}(q)_{N_k−N_{k+1}}(Q;Q)_{N_{k+1}}).
3. Change of a ((B-iv) with c = 1, a = x): φ_d(·;xq) = Γ^{(1)}_{d,d}(x)φ_d(·;x) + Γ^{(1)}_{d,d+1}(x)φ_{d+1}(·;x), Γ^{(1)}_{d,d}(a) = (1−aq)/(1−aq^{2d+1}),
   Γ^{(1)}_{d,d+1}(a) = −aq^{2d+1}(1−aq)/(1−aq^{2d+1}) (the two-line partial-fraction computation of the case c = 1 in
   `F-all-k.md` (B-iv)). So (α*, β″) is a Bailey pair relative to x with
   α*_e = (1−xq)α″_e/(1−xq^{2e+1}) − xq^{2e−1}(1−xq)α″_{e−1}/(1−xq^{2e−1}), i.e. α*_0 = 1 and, for e ≥ 1,
   α*_e = (−1)^e[x^{2e}q^{3e²+2e}(x²Q;Q)_e/(Q;Q)_e + x^{2e−1}q^{3e²−2e}(x²Q;Q)_{e−1}/(Q;Q)_{e−1}].
   With y := q^{2e}: α*_e = (−1)^e x^{2e−1}q^{3e²−2e}(x²Q;Q)_{e−1}(1 − y + xy² − x³y³)/(Q;Q)_e, and 1 − y + xy² − x³y³ = (1−xy)(1+(x−1)y+x²y²),
   so α*_e = λ_e² b^Ω_e (also for e = 0). [D5]
4. k−1 steps (L∞) relative to x: α^{(k)}_n = λ_n^{k−1}α*_n = λ_n^{k+1}b^Ω_n, and β^{(k)}_n is the sum defining B^{(k)}_{k+1} with the extra
   factor 1/(q)_{n−N_1} and N_1 ≤ n (for k = 1 there is no step, and β^{(1)} = β″ of item 2, whose sum is over
   n ≥ N_1 ≥ N_2).
5. Let n → ∞ (`KRTW-F-moves.md` §4, (lim)): for each power x^D only finitely many index tuples and finitely many r contribute
   (ord α^{(k)}_r ≥ r), and 1/(q)_{n−j}, 1/(xq)_{n+r} converge q-adically, coefficientwise in x; so B^{(k)}_{k+1}(x)/(q)_∞ = Σ_rα^{(k)}_r/((q)_∞(xq)_∞).

Hence (xq)_∞B^{(k)}_{k+1}(x) = Σ_nλ_n^{k+1}b^Ω_n, and Theorem 4 gives K^Ω_k = B^{(k)}_{k+1}J_k. ∎ [D7: k ≤ 8, mod q^40, all powers of x.]

(At x = 1 Theorem 6 is Bressoud's even-modulus identity B^{(k)}_{k+1}(1) = (q^k,q^{k+4},q^{2k+4};q^{2k+4})_∞/(q)_∞, which §5 does not need.)

## 7. Verification summary

`code/certify_iv_all_k.py 40 8` (ALL CHECKS PASSED, `certificates/certify_iv_all_k.txt`, 13 s): D1 Lemma 2.1, D2 Lemma 2.2
(exact, point masses); D3 Proposition 2.3 (direct summation, k ≤ 8, mod q^40, all powers of x; negative control); D4/D4b
Lemma 3.1 (N ≤ 5, symbolic); D5 Lemma 4.1 (N ≤ 12) and α* = λ²b^Ω (n ≤ 12); D6 Theorem 3.2 (k ≤ 8 mod q^40; negative
control); D7 Theorem 6 (k ≤ 8 mod q^40); D8 the triple-product evaluation (k ≤ 10 mod q^200). Independent check of (iv)
itself to q^300 for k ≤ 8: `certificates/series_check_q300.txt`.
