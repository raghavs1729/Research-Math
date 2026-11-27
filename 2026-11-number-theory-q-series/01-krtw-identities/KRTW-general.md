# The first Kanade–Russell–Tsuchioka–Warnaar Andrews–Gordon-type identity for all k

> **Referee notes (2026-11-26, `REFEREE-KRTW-general.md`; they supersede any contrary statement below).**
> No error found; every lemma re-derived and checked exactly; Theorem 3 confirmed by independent brute force
> (k ≤ 4 mod q^500, k = 5 mod q^300); Corollary 4 confirmed including a direct computation of the Hall–Littlewood
> side (k = 2,3,4 mod q^26), so its caveats can be dropped (k = 0 is the trivial case; GOW needs m ≥ 1).
> **Correction to §5:** the evaluation of the k-fold Melzer-type sum as the product is in the literature for
> ALL k (Warnaar arXiv:2511.17034, display after Corollary 5, case i = 0, attributed to Andrews 2010, Kim–Yee
> 2013, Chern–Li–Stanton–Xue–Yee 2024, Kanade–Russell 2025), so Proposition 6 is a re-proof. The new content is
> the first equality of Theorem 2 (the reduction of the 2k-fold KRTW sum to the k-fold sum, with its
> x-refinement), Theorem 1 and Lemma 4 — and hence KRTW identity (i) for k ≥ 3 and Conjecture 2.2 at z = 1.
> The minor points of both referee reports have been applied in the text (2026-11-26, third session): the
> "known before" sentence, the literature paragraph of §5 (Kanade–Lovejoy, Warnaar, Dousse–Jang,
> Chern–Li–Stanton–Xue–Yee added; the new part is the reduction), the proof of J_k = C_k/C_1 in §8, Remark (a)
> (only the first form of (E_k) is proved), and "Melzer (2.6), s = 2".


Status: **complete proof of identity (i) for every k ≥ 1** (written 2026-11-26), via two refinements that are
proved on the way: an x-refinement (Theorem 2, the former "conjecture (E_k)" of `KRTW-k2.md` §4) and a finite
R-refinement (Theorem 1). Identities (ii), (iii), (iv) are **not** proved: they are tested (§8) and reduced to a
new x-refined factorisation conjecture (F_k) (§8, numerical; `code/krtw_Fk_test.py`). Ingredients: the
q-binomial theorem, Euler's two series, Bailey's lemma (one limiting case) with the unit Bailey pair, and the
Jacobi triple product. Exact certificate: `code/krtw_general_certify.py` (log
`certificates/krtw_general_certify.txt`); high-order test of all four identities: `code/krtw_four_test.py` (log
`certificates/krtw_four_test.txt`). Labels like [L4], [T1] below refer to the checks of the certificate.

Source. S. Kanade, M. C. Russell, S. Tsuchioka, S. O. Warnaar, *Remarks on the conjectures of Capparelli,
Meurman, Primc and Primc*, arXiv:2404.03851v2 (Selecta Math. 2026). After Conjecture 4.9 they display four
"conjectural Andrews–Gordon-type identities". With r_{k+1} = s_0 := 0, 1/(q)_m := 0 for m < 0 and

    W(r,s) := ∏_{i=1}^k q^{(r_i+s_i)²+s_i²} / ( (q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} ),

(the sums run over r_1 ≥ … ≥ r_k ≥ 0 and 0 ≤ s_1 ≤ … ≤ s_k) they read

    (i)    Σ q^{Σ s_i} W          = (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ / ( (q;q²)_∞ (q;q)_∞ ),
    (ii)   Σ q^{Σ (r_i+2s_i)} W   = (q,q²,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4}; q^{2k+4})_∞ / ( (q²;q²)_∞ (q;q)_∞ ),
    (iii)  Σ W                    = (q^{k+1},q^{k+2},q^{k+2},q^{k+3},q^{2k+4},q^{2k+4}; q^{2k+4})_∞ / ( (q²;q²)_∞ (q;q)_∞ ),
    (iv)   Σ Ω W                  = (q^k,q^{k+1},q^{k+3},q^{k+4},q^{2k+4},q^{2k+4}; q^{2k+4})_∞ / ( (q²;q²)_∞ (q;q)_∞ ),
           Ω = Σ_{i<k} q^{r_i+2s_i}(1−q^{2s_{i+1}−2s_i}) + q^{r_k+2s_k}.

Known before: k = 1 for (i), (ii), (iv) (KRTW: Conjecture 4.9 holds for k = 1 and the level-one products are
classical; for (iii) at k = 1 the sum is Σ q^{(r+s)²+s²}/((q)_r(q²;q²)_s) = Bressoud's (q³,q³,q⁶;q⁶)_∞/(q)_∞, which is
the right side, but KRTW do not discuss this case separately); k = 2: (i) by KRTW Theorem 4.10 + their (4.7) +
Meurman–Primc, (ii)–(iv) by Theorem 4.10 + Bartholdi, arXiv:2611.04020 (three-row, level-two case of Primc's
conjecture; this deduction is ours, it is not stated in either paper). For k ≥ 3 I found no proof in the
literature (searches: §9). What is proved here:

**Theorem 3.** Identity (i) holds for every integer k ≥ 1 (and trivially for k = 0).

**Theorem 2 (x-refinement).** For every k ≥ 0, in Z[x][[q]],

    Σ_{r_1≥…≥r_k≥0, 0≤s_1≤…≤s_k} x^{Σ_i (r_i+2s_i)} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ( (q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} )
        = Σ_{M_1≥…≥M_k≥0} x^{M_1+…+M_k} q^{Σ_i M_i(M_i+1)/2} / ∏_{i=1}^k (q)_{M_i−M_{i+1}}                      (M_{k+1} := 0)
        = ( (−xq;q)_∞ / (x²q;q)_∞ ) Σ_{n≥0} (−1)^n x^{(k+1)n} q^{(k+1)n(n+1)/2 + n(n−1)/2} (x²q;q)_n (1−x²q^{2n+1}) / (q;q)_n .

Corollaries 4 (the z = 1 case of KRTW Conjecture 2.2) and 5 (a second parameter y marking r_1) are in §6.

Notation: (a;q)_n = ∏_{i=0}^{n−1}(1−aq^i), (q)_n = (q;q)_n, [n;k] = (q)_n/((q)_k(q)_{n−k}) (0 ≤ k ≤ n; := 0 otherwise),
N_i := r_i + 2s_i. The elementary identity used throughout:

    (r+s)² + s² + s = N(N+1)/2 + r(r−1)/2     for N = r + 2s.                                                   (0.1)

## 1. Formal setting

All identities are identities in the ring 𝓡 = Z[x][[q]] (formal power series in q whose coefficients are
polynomials in x) with its q-adic topology, or in Q(q)[x] when they are finite. A family (f_ι) in 𝓡 is
*summable* if for every n only finitely many f_ι are ≢ 0 mod q^n; summable families may be regrouped and
their order of summation interchanged. Every multiple series below is a summable family: its general term is
x^a q^e·u with u ∈ 𝓡, where e tends to infinity along the index set (e.g. e = Σ_i ((r_i+s_i)² + s_i² + s_i) ≥ Σ_i (r_i² + s_i²),
or e = Σ_i M_i(M_i+1)/2, or e ≥ (r+2b)(r+2b+1)/2 + Σ_i M_i(M_i+1)/2 in §3), so only finitely many index tuples have e < n. The substitutions x ↦ xq^{c} (c ≥ 0) and x ↦ 1 are continuous
ring homomorphisms 𝓡 → 𝓡, resp. 𝓡 → Z[[q]]. 1/(q)_n, 1/(q²;q²)_n ∈ Z[[q]], [n;k] ∈ Z[q], (xq^c;q)_∞ ∈ 𝓡 for c ≥ 1.

## 2. Four lemmas

**Lemma 1 (q-binomial theorem and Euler's series).** For n ≥ 0 and an indeterminate z,

    (z;q)_n = Σ_{i=0}^n [n;i] (−1)^i q^{i(i−1)/2} z^i ,      1/(z;q)_∞ = Σ_{n≥0} z^n/(q)_n ,      (z;q)_∞ = Σ_{n≥0} (−1)^n q^{n(n−1)/2} z^n/(q)_n .

(G. E. Andrews, *The Theory of Partitions*, (3.3.6), (2.2.5), (2.2.6). The first is a polynomial identity in z
over Z[q], so z may be specialised to any element of Z[q,q⁻¹][x].) [L1]

**Lemma 2.** For integers R, L ≥ 0:

    Σ_{t,b≥0, t+2b=L} q^{t(t−1)/2} [R;t] / (q²;q²)_b  =  Σ_{j=0}^{L} (−1)^j q^{Rj} / ( (q)_j (q)_{L−j} ).              (2.1)

*Proof.* By Lemma 1 (with z ↦ −z, and Euler's first series in base q²) the left side is the coefficient of z^L in
(−z;q)_R/(z²;q²)_∞. Since (z²;q²)_∞ = (z;q)_∞(−z;q)_∞ and (−z;q)_R/(−z;q)_∞ = 1/(−zq^R;q)_∞,

    (−z;q)_R/(z²;q²)_∞ = 1/( (z;q)_∞ (−zq^R;q)_∞ ) = Σ_{l≥0} z^l/(q)_l · Σ_{j≥0} (−q^R z)^j/(q)_j ,

whose coefficient of z^L is the right side. ∎ [L2: exact for R ≤ 8, L ≤ 10]

(For R → ∞ this is the self-fusion Lemma S of `KRTW-k2.md`; the point of (2.1) is that the "boundary" R enters the
right side only through the exponentials q^{Rj}.)

**Lemma 3.** For integers R, d ≥ 0 and an indeterminate y:

    Σ_{r=0}^{R} [R;r] y^r q^{r²} (yq^{d+r+1};q)_{R−r}  =  Σ_{t=0}^{min(R,d)} [R;t] [d;t] (q)_t y^t q^{t²}.                (2.2)

*Proof.* Expand (yq^{d+r+1};q)_{R−r} = Σ_i [R−r;i](−1)^i q^{i(i−1)/2} y^i q^{(d+r+1)i} (Lemma 1), put t = r + i and use
[R;r][R−r;i] = [R;t][t;i]. With r = t − i the exponent of q is
(t−i)² + i(i−1)/2 + (d+t−i+1)i = t² + i(i−1)/2 + i(d+1−t), so the left side equals

    Σ_{t=0}^{R} [R;t] y^t q^{t²} Σ_{i=0}^{t} [t;i](−1)^i q^{i(i−1)/2} (q^{d+1−t})^i = Σ_{t=0}^{R} [R;t] y^t q^{t²} (q^{d+1−t};q)_t

(Lemma 1 again). For t > d the product (q^{d+1−t};q)_t contains the factor 1 − q⁰ = 0; for t ≤ d it equals
(q)_d/(q)_{d−t} = [d;t](q)_t. ∎ [L3: exact for R, d ≤ 7] (For d = 0 this is the classical Σ_r [R;r]y^r q^{r²}/(yq)_r = 1/(yq)_R.)

**Lemma 4 (the key identity).** For integers R, M ≥ 0, in Q(q)[x]:

    Σ_{r=0}^{R} Σ_{b=0}^{⌊M/2⌋} [R;r] x^{r+2b} q^{(r+2b)(r+2b+1)/2 + r(r−1)/2} (xq^{r+M+1};q)_{R−r} / ( (q²;q²)_b (q)_{M−2b} )
        =  Σ_{L=0}^{M} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_{M−L} / ( (q)_L (q)_{M−L} ).                                       (∗)

Consequently, in 𝓡 (multiply (∗) by (xq^{R+M+1};q)_∞):

    Σ_{r=0}^{R} Σ_{b=0}^{⌊M/2⌋} [R;r] x^{r+2b} q^{(r+2b)(r+2b+1)/2 + r(r−1)/2} (xq^{r+M+1};q)_∞ / ( (q²;q²)_b (q)_{M−2b} )
        =  Σ_{L=0}^{M} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_∞ / ( (q)_L (q)_{M−L} ),                                           (∗∗)

and ("M = ∞")

    Σ_{r=0}^{R} Σ_{b≥0} [R;r] x^{r+2b} q^{(r+2b)(r+2b+1)/2 + r(r−1)/2} / (q²;q²)_b  =  Σ_{L≥0} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_∞ / (q)_L .   (∗∞)

*Proof of (∗).* Both sides will be brought to the same form

    Σ_{L=0}^{M} ( x^L q^{L(L+1)/2} / (q)_{M−L} ) · Σ_{j=0}^{L} (−1)^j q^{Rj} / ( (q)_j (q)_{L−j} ).                          (2.3)

Left side. By (0.1) the exponent is r² + 2br + 2b² + b, so the left side is

    Σ_b ( x^{2b} q^{2b²+b} / ((q²;q²)_b (q)_{M−2b}) ) Σ_r [R;r] (xq^{2b})^r q^{r²} ( (xq^{2b}) q^{(M−2b)+r+1} ; q )_{R−r} .

By Lemma 3 with y = xq^{2b}, d = M − 2b ≥ 0, the inner sum is Σ_t [R;t][M−2b;t](q)_t x^t q^{2bt+t²}, and
[M−2b;t](q)_t/(q)_{M−2b} = 1/(q)_{M−2b−t}. Hence the left side equals

    Σ_{t,b≥0, t+2b≤M} x^{t+2b} q^{t²+2bt+2b²+b} [R;t] / ( (q²;q²)_b (q)_{M−2b−t} )
        = Σ_{L=0}^{M} ( x^L q^{L(L+1)/2} / (q)_{M−L} ) Σ_{t+2b=L} q^{t(t−1)/2} [R;t] / (q²;q²)_b

(L = t + 2b and (0.1) once more), which is (2.3) by Lemma 2.
Right side. Expand (xq^{R+L+1};q)_{M−L} = Σ_i [M−L;i](−1)^i q^{i(i−1)/2} x^i q^{(R+L+1)i} (Lemma 1); since
L(L+1)/2 + i(i−1)/2 + (L+1)i = (L+i)(L+i+1)/2, the right side is

    Σ_{L,i≥0, L+i≤M} (−1)^i q^{Ri} x^{L+i} q^{(L+i)(L+i+1)/2} / ( (q)_L (q)_i (q)_{M−L−i} ),

which is (2.3) after renaming (L+i, i) → (L, j). ∎
*Proof of (∗∗).* (xq^{r+M+1})_{R−r}(xq^{R+M+1})_∞ = (xq^{r+M+1})_∞ and (xq^{R+L+1})_{M−L}(xq^{R+M+1})_∞ = (xq^{R+L+1})_∞. ∎
*Proof of (∗∞).* By (0.1) and Lemma 2 the left side is Σ_L x^L q^{L(L+1)/2} Σ_{j+l=L} (−1)^j q^{Rj}/((q)_j(q)_l)
= Σ_l (x^l q^{l(l+1)/2}/(q)_l) Σ_j (−1)^j q^{j(j−1)/2}(xq^{R+l+1})^j/(q)_j (as (j+l)(j+l+1)/2 = l(l+1)/2 + j(j−1)/2 + j(l+1)),
and the inner sum is (xq^{R+l+1};q)_∞ by Lemma 1. ∎
[L4: (∗) exactly in Z[x,q] after clearing denominators, for 0 ≤ R, M ≤ 7; L4s: (∗∗) for R, M ≤ 8 and (∗∞) for R ≤ 8 mod q^36.]

## 3. The R-refined identity

**Theorem 1.** For integers j ≥ 1 and R ≥ 0, in 𝓡:

    A_j(R;x) := Σ_{R≥r_1≥…≥r_j≥0, 0≤s_1≤…≤s_j} ( (q)_R / (q)_{R−r_1} ) ∏_{i=1}^{j} x^{r_i+2s_i} q^{(r_i+s_i)²+s_i²+s_i} / ( (q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} )
      = B_j(R;x) := Σ_{M_1≥…≥M_j≥0} x^{M_1+…+M_j} q^{Σ_i M_i(M_i+1)/2} (xq^{R+M_j+1};q)_∞ / ∏_{i=1}^{j} (q)_{M_i−M_{i+1}}

(r_{j+1} = s_0 = M_{j+1} := 0). Only the innermost variable M_j feels the bound R.

*Proof.* Induction on j. For j = 1, (q)_R/((q)_{R−r}(q)_r) = [R;r] and (0.1) show that A_1(R;x) is the left side of
(∗∞), and B_1(R;x) is its right side.

Let j ≥ 2. In A_j write r = r_1, b = s_1 and s_i = b + s'_i for i ≥ 2 (so 0 = s'_1 ≤ s'_2 ≤ … ≤ s'_j, and
s_i − s_{i−1} = s'_i − s'_{i−1} for i ≥ 2), N'_i = r_i + 2s'_i. By (0.1) the exponent of level i ≥ 2 is
N_i(N_i+1)/2 + r_i(r_i−1)/2 with N_i = N'_i + 2b, and N_i(N_i+1)/2 = N'_i(N'_i+1)/2 + 2bN'_i + (2b²+b). Since
(q)_R/((q)_{R−r_1}(q)_{r_1−r_2}) = [R;r_1]·(q)_{r_1}/(q)_{r_1−r_2}, regrouping gives

    A_j(R;x) = Σ_{r=0}^{R} Σ_{b≥0} [R;r] ( x^{r+2b} q^{(r+2b)(r+2b+1)/2 + r(r−1)/2} / (q²;q²)_b ) · x^{2b(j−1)} q^{(j−1)(2b²+b)} · A_{j−1}(r; xq^{2b}).     (3.1)

By the induction hypothesis A_{j−1}(r;xq^{2b}) = B_{j−1}(r;xq^{2b}). In B_{j−1}(r;xq^{2b}) put M'_i = M_i + 2b (1 ≤ i ≤ j−1):

    x^{2b(j−1)} q^{(j−1)(2b²+b)} · x^{ΣM_i} q^{2bΣM_i + ΣM_i(M_i+1)/2} = x^{ΣM'_i} q^{ΣM'_i(M'_i+1)/2},      (xq^{2b}·q^{r+M_{j−1}+1};q)_∞ = (xq^{r+M'_{j−1}+1};q)_∞ ,

the differences M_i − M_{i+1} = M'_i − M'_{i+1} (i ≤ j−2) are unchanged and (q)_{M_{j−1}} = (q)_{M'_{j−1}−2b}. Hence, summing over b
inside (allowed: summable family),

    A_j(R;x) = Σ_{M'_1≥…≥M'_{j−1}≥0} ( x^{ΣM'_i} q^{ΣM'_i(M'_i+1)/2} / ∏_{i=1}^{j−2} (q)_{M'_i−M'_{i+1}} )
                 × Σ_{r=0}^{R} Σ_{b=0}^{⌊M'_{j−1}/2⌋} [R;r] x^{r+2b} q^{(r+2b)(r+2b+1)/2 + r(r−1)/2} (xq^{r+M'_{j−1}+1};q)_∞ / ( (q²;q²)_b (q)_{M'_{j−1}−2b} ).

The second line is the left side of (∗∗) with M = M'_{j−1}; replacing it by the right side of (∗∗) and calling the
summation index L =: M_j (so M'_{j−1} ≥ M_j ≥ 0, with denominators (q)_{M'_{j−1}−M_j}(q)_{M_j}) gives exactly B_j(R;x). ∎
[T1: both sides, the left one summed directly from its definition, for j ≤ 4 and R ≤ 6, 5, 4, 3, mod q^36 in all powers of x.]

## 4. The x-refinement

Let KRTW_k(x) and Mel_k(x) denote the first and the second expression in Theorem 2.

*Proof of the first equality of Theorem 2.* k = 0 is 1 = 1. Let k ≥ 1 and R ≥ 0. The terms of KRTW_k(x) and of A_k(R;x)
are indexed by the same (r,s) (those with r_1 ≤ R in A_k), a term has q-order ≥ r_1², and
(q)_R/(q)_{R−r_1} = (q^{R−r_1+1};q)_{r_1} ≡ 1 mod q^{R−r_1+1}. As r_1² + R − r_1 + 1 ≥ R + 1 and r_1² ≥ R + 1 for r_1 > R,
A_k(R;x) ≡ KRTW_k(x) mod q^{R+1}. Likewise (xq^{R+M_k+1};q)_∞ ≡ 1 mod q^{R+1}, so B_k(R;x) ≡ Mel_k(x) mod q^{R+1}.
By Theorem 1, KRTW_k(x) ≡ Mel_k(x) mod q^{R+1} for every R. ∎
[T2a: k ≤ 3 by brute force; T2b: k ≤ 7 with the left side summed level by level; all powers of x, mod q^36.]

Remarks. (a) This proves the conjecture (E_k) of `KRTW-k2.md` §4 in its first form (the one displayed in
Theorem 2). The second form mentioned there (s_1 = 0 on the left, last weight q^{M_k²} on the right) is NOT
proved here: deducing it from the first form would need an injectivity argument that is not given; it has only
been checked numerically (k ≤ 5). (b) For k = 2, Theorem 2 is
equivalent to the identity (E) of `KRTW-k2.md` (proved there by a 14-move chain found by computer search); it now
has a short proof: (∗∞), one application of (∗∗), and R → ∞.

## 5. The Melzer-type sum: Bailey's lemma

A pair of sequences (α_n, β_n)_{n≥0} is a *Bailey pair relative to a* if β_n = Σ_{r=0}^{n} α_r/((q)_{n−r}(aq;q)_{n+r}) for all
n ≥ 0. Two classical facts (D. Bressoud, M. Ismail, D. Stanton, *Change of base in Bailey pairs*, Ramanujan J. 4
(2000), arXiv:math/9909053: Bailey's Lemma as stated in §1 there, its limiting case (S2) of §3 (ρ_1 → ∞, ρ_2 = −√(aq)),
and the unit Bailey pair (UBP); originally W. N. Bailey 1949, G. E. Andrews, Pacific J. Math. 114 (1984)):

* (UBP) β_n = δ_{n,0}, α_n = (−1)^n q^{n(n−1)/2} (a;q)_n (1−aq^{2n}) / ((q)_n (1−a)) is a Bailey pair relative to a.
* (S2) If (α,β) is a Bailey pair relative to a, then so is
  α'_r = a^{r/2} q^{r²/2} α_r,   β'_n = Σ_{k=0}^{n} (−√(aq);q)_k a^{k/2} q^{k²/2} β_k / ( (q)_{n−k} (−√(aq);q)_n ).

Both are identities of rational functions of q and √a (a generic). We use them with **a = x²q**, √(aq) = xq,
a^{r/2}q^{r²/2} = x^r q^{r(r+1)/2}; then 1/(aq;q)_m = 1/(x²q²;q)_m ∈ 𝓡 and (a;q)_n/(1−a) = (x²q²;q)_{n−1} (n ≥ 1), so everything
lies in 𝓡.

**Proposition 6.** For every k ≥ 0, in 𝓡,

    Mel_k(x) = ( (−xq;q)_∞ / (x²q;q)_∞ ) Σ_{n≥0} (−1)^n x^{(k+1)n} q^{(k+1)n(n+1)/2 + n(n−1)/2} (x²q;q)_n (1−x²q^{2n+1}) / (q)_n .       (5.1)

*Proof.* Let (α^{(0)},β^{(0)}) be (UBP) with a = x²q and let (α^{(i)},β^{(i)}) be obtained from it by i applications of
(S2). Then α^{(i)}_r = x^{ir} q^{ir(r+1)/2} α^{(0)}_r, and by induction on i ≥ 1

    β^{(i)}_n = (1/(−xq;q)_n) Σ_{n≥M_2≥…≥M_i≥0} ∏_{l=2}^{i} x^{M_l} q^{M_l(M_l+1)/2} / ( (q)_{n−M_2} (q)_{M_2−M_3} ⋯ (q)_{M_i} )

(for i = 1: β^{(1)}_n = 1/((q)_n(−xq)_n); in the step the factor (−xq)_k of (S2) cancels the prefactor 1/(−xq)_k of
β^{(i)}_k). The Bailey pair relation for i = k + 1 reads

    (1/(−xq)_n) Σ_{n≥M_1≥…≥M_k≥0} ∏_{l=1}^{k} x^{M_l} q^{M_l(M_l+1)/2} / ( (q)_{n−M_1} ⋯ (q)_{M_k} ) = Σ_{r=0}^{n} x^{(k+1)r} q^{(k+1)r(r+1)/2} α^{(0)}_r / ( (q)_{n−r} (x²q²;q)_{n+r} ).

Let n → ∞ in 𝓡: on the left 1/((−xq)_n (q)_{n−M_1}) ≡ 1/((−xq)_∞(q)_∞) mod q^{n−M_1+1} and the term has q-order ≥ M_1, on the
right 1/((q)_{n−r}(x²q²)_{n+r}) ≡ 1/((q)_∞(x²q²)_∞) mod q^{n−r+1} and the term has q-order ≥ r; the omitted terms (M_1 > n,
r > n) have q-order > n. Hence Mel_k(x)/((−xq)_∞(q)_∞) = (1/((q)_∞(x²q²;q)_∞)) Σ_{r≥0} x^{(k+1)r} q^{(k+1)r(r+1)/2} α^{(0)}_r, which is
(5.1) because (x²q²;q)_∞(1−x²q) = (x²q;q)_∞. ∎
[B1: the Bailey pair relation for (α^{(i)},β^{(i)}), i ≤ 3, n ≤ 4, mod q^36; B2: (5.1) for k ≤ 7 in all powers of x, mod q^36.]

This also proves the second equality of Theorem 2.

Identification in the literature (corrected after the second referee report, `REFEREE-2-KRTW.md`). At x = 1
the evaluation Mel_k(1) = ((−q)_∞/(q)_∞)(q, q^{k+1}, q^{k+2}; q^{k+2})_∞ of the reduced k-fold sum is KNOWN for every k:
Chern–Li–Stanton–Xue–Yee (J. Algebraic Combin. 2024, Theorem 3.5); Kanade–Lovejoy, arXiv:2410.15408 (SIGMA 21
(2025) 021), Theorem 1.1 with a = 0 (a Bailey-pair proof, all k); Warnaar, arXiv:2511.17034 (SIGMA 22 (2026) 062),
the display after Corollary 5, case i = 0, where it is attributed to Andrews (2010, Theorem 3), Kim–Yee (2013,
Theorem 1.4), Chern–Li–Stanton–Xue–Yee and Kanade–Russell; Dousse–Jang, arXiv:2603.05300, Theorem "thm:AK"
(a = 0). For even k it is also Bressoud–Ismail–Stanton, Theorem 5.2 (i = 1; this identification was not
re-checked by the referee), i.e. the case s = 2 of Melzer's conjecture (hep-th/9412154, (2.6)). So Proposition 6
is a re-proof (I did not find the x-refined form (5.1) itself stated in these papers, but it is the same Bailey
chain (S2)^{k+1}(UBP)), and **no novelty is claimed for §5**. The new part of Theorem 2 is its FIRST equality — the
reduction of the 2k-fold KRTW sum to the k-fold sum, with its x-refinement — together with Theorem 1 and Lemma 4.

**Proof of Theorem 3.** Put x = 1 in Theorem 2. The left side becomes the left side of (i) (the summand is
q^{Σ s_i}W). In (5.1), (q)_n/(q)_n cancels and

    Mel_k(1) = ((−q)_∞/(q)_∞) Σ_{n≥0} (−1)^n q^{((k+2)n²+kn)/2} (1 − q^{2n+1}) = ((−q)_∞/(q)_∞) Σ_{n∈Z} (−1)^n q^{((k+2)n²+kn)/2}

(the terms with −q^{2n+1} are the terms n ↦ −n−1 of the bilateral sum:
((k+2)(n+1)² − k(n+1))/2 = ((k+2)n²+kn)/2 + 2n + 1). By the Jacobi triple product
Σ_{n∈Z} (−1)^n z^n Q^{n(n−1)/2} = (z, Q/z, Q; Q)_∞ with Q = q^{k+2}, z = q^{k+1}, the sum is (q, q^{k+1}, q^{k+2}; q^{k+2})_∞, and
(−q;q)_∞ = 1/(q;q²)_∞ (Euler). ∎
[T3a: Mel_k(1) = product for k ≤ 20 mod q^300; T3b: the 2k-fold sum (i) = product for k ≤ 8 mod q^300; also
`krtw_four_test.txt`: k ≤ 8 mod q^400.]

## 6. Corollaries

**Corollary 4 (KRTW Conjecture 2.2 at z = 1).** For every k ≥ 0,

    Σ_{λ: λ_1≤k} q^{|λ|} P_{2λ}(1,q,q²,…;q²) = Σ_{r,s} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ( (q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} ),

i.e. KRTW's Conjecture 2.2 (their (2.12)) holds at z = 1 (P_λ(·;t) the Hall–Littlewood polynomial). *Proof.* By
Griffin–Ono–Warnaar, arXiv:1401.7718, Theorem 1.2, eq. (1.9) with n = 1, m = k (κ = 2k+4), the left side equals
(q²;q²)_∞ (q^{k+2};q^{k+2})_∞ θ(q;q^{k+2}) / (q)_∞² = (q,q^{k+1},q^{k+2};q^{k+2})_∞/((q;q²)_∞(q)_∞); by Theorem 3 so does the right side. ∎
(Caveats: I read (1.9) in a text extraction of the arXiv PDF; the reading agrees with KRTW's own statement that
(i) follows from their conjectures and with the fact that the product is the Meurman–Primc / Lepowsky–Milne
principally specialised A_1^{(1)} character of level k. The Hall–Littlewood side was not computed numerically
here. For general z Conjecture 2.2 remains open: the grading z^{Σ(r_i+s_i)} is not the grading x^{Σ(r_i+2s_i)} of Theorem 2.)

In KRTW's notation Theorem 3 also says: the first multisum of Conjecture 4.9, at z = 1/q, equals
𝓓^{(2)}_{kΛ_0}(1/q,q) = 𝓒^{(1)}_{kΛ_0}(1,q) (their (4.7) and the Meurman–Primc theorem); i.e. Conjecture 4.9 (first line) holds at z = 1/q.

**Corollary 5 (a second parameter).** For k ≥ 1, in Z[x,y][[q]],

    Σ_{r,s} y^{r_1} x^{Σ(r_i+2s_i)} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ( (q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} )
       = Σ_{M_1≥…≥M_k≥0} Σ_{l≥0} (−1)^l (y;q)_l x^{l+ΣM_i} q^{l(l+1)/2 + lM_k + ΣM_i(M_i+1)/2} / ( (q)_l ∏_{i=1}^k (q)_{M_i−M_{i+1}} ).

*Proof.* Work in Z[x][[q]][[y]]. Multiply Theorem 1 (j = k) by y^R/(q)_R and sum over R ≥ 0. Left: Σ_{R≥r_1} y^R/(q)_{R−r_1} =
y^{r_1}/(y;q)_∞. Right: expand (xq^{R+M_k+1};q)_∞ = Σ_l (−x)^l q^{l(l+1)/2+lM_k} q^{Rl}/(q)_l (Lemma 1) and use Σ_R y^R q^{Rl}/(q)_R = 1/(yq^l;q)_∞.
Multiply by (y;q)_∞ and note (y)_∞/(yq^l)_∞ = (y)_l. Both sides have coefficients of q^n that are polynomials in x, y. ∎
y = 1 is Theorem 2; y = 0 (all r_i = 0) is Theorem 1 with R = 0:
Σ_{0≤s_1≤…≤s_k} x^{2Σs_i} q^{Σ(2s_i²+s_i)}/∏(q²;q²)_{s_i−s_{i−1}} = Σ_M x^{ΣM_i} q^{ΣM_i(M_i+1)/2} (xq^{M_k+1};q)_∞/∏(q)_{M_i−M_{i+1}}.
[C5: k ≤ 3 at y = 0, −1, q, q², q³, −q, all powers of x, mod q^36.]

## 7. Verification

`certificates/krtw_general_certify.txt` (all "OK"; exact arithmetic with python-flint): L1–L3 and (∗) as polynomial
identities in Z[x,q]; (∗∗), (∗∞), Theorem 1, Theorem 2, Corollary 5, the Bailey pairs and (5.1) modulo q^36 in all
powers of x; Theorem 3 modulo q^300 (k ≤ 8 for the 2k-fold sum, k ≤ 20 for Mel_k(1)); X1, X2: the two numerical
remarks of §8. The 2k-fold sums are evaluated by two independent routines: brute-force enumeration of (r,s) for
small k, and a level-by-level transfer recursion (the definition summed one level at a time), which is compared
with the brute force in T2a/T2b and in the self-test printed at the top of `krtw_four_test.txt` (k ≤ 3, mod q^30,
all four sums; a deliberately wrong product is detected).
`certificates/krtw_four_test.txt`: all four identities (i)–(iv) hold modulo q^400 for k = 1,…,8.

## 8. The other three identities — not proved

Tested: (ii), (iii), (iv) hold modulo q^400 for 1 ≤ k ≤ 8 (`krtw_four_test.py`; no counterexample).

Why the method stops. Give level i the weight X^{N_i} q^{N_i(N_i+1)/2} · c^{r_i} q^{r_i(r_i−1)/2} and call the resulting sum
K_k(X,c). In KRTW's two-variable normalisation z^{Σ(r_i+s_i)}: c² = zq; 𝓓^{(2)}_{kΛ_0}(z) ↔ X = c, 𝓓^{(2)}_{kΛ_1}(z) ↔ X = c/q
(Conjecture 4.9), Conjecture 2.2 ↔ K_k(√z, √z). Theorem 2 evaluates K_k(X,1), the line c = 1 (z = 1/q); (i) is the
point X = c = 1. Identities (iii), (ii) are the points c = q^{1/2}, X = q^{∓1/2} (z = 1), and (iv) is a combination at
c = q^{1/2}. With a general c, Lemma 2 becomes the coefficient extraction from

    (−cz;q)_R / (z²;q²)_∞ = (−cz;q)_R / ( (z;q)_∞ (−z;q)_∞ ),

and the cancellation (−z;q)_R/(−z;q)_∞ = 1/(−zq^R;q)_∞ that drives Lemmas 2–4 happens only for c ∈ q^Z (c = 1 here;
c = q^{±1},… give variants on the lines z = q, q^{−3}, …, never z = 1). At c = q^{1/2} I found no closed form for the
R-refined sums. The product sides confirm that more is needed: (iii) is a product of two theta quotients,

    (iii) = [ (q^{k+2},q^{k+2},q^{2k+4};q^{2k+4})_∞/(q)_∞ ] · [ (q^{k+1},q^{k+3},q^{2k+4};q^{2k+4})_∞/(q²;q²)_∞ ]
          = [ Σ_{0≤s_1≤…≤s_k} q^{Σ s_i²}/∏(q²;q²)_{s_i−s_{i−1}} ] · (−q²;q²)_∞ (q^{k+2},q^{k+2},q^{2k+4};q^{2k+4})_∞/(q²;q²)_∞

(second form: the evaluation Σ q^{(M_1²+…+M_k²)/2}/∏(q)_{m_i} = ((−q^{1/2})_∞/(q)_∞)(q^{(k+1)/2},q^{(k+3)/2},q^{k+2};q^{k+2})_∞ displayed before
Theorem 4.1 of Warnaar, math/0110307, in base q²; the first bracket is the s-chain of (iii) alone; [X2]: k ≤ 8 mod
q^200), whereas (i) is a single theta quotient times (−q)_∞. A proof of (ii)–(iv) along these lines would need a
two-theta ("rank two") analogue of Theorem 1; none was found. Note also that the product sides of (ii)–(iv) come
from Primc's conjecture for three-row arrays, open for k ≥ 3, so (ii)–(iv) would not by themselves prove anything
about 𝓓^{(2)} partitions without Conjecture 4.9.

**A reduction of (ii), (iii), (iv) to an x-refined factorisation (conjectural step labelled).** Put

    K_k(x)      := Σ_{r,s} x^{Σ(r_i+2s_i)} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²} / ( (q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} )          ((iii) is K_k(1), (ii) is K_k(q)),
    K^Ω_k(x)    := the same sum with the factor Ω inserted                                                         ((iv) is K^Ω_k(1)),
    B^{(a)}_{k+1}(x) := Σ_{N_1≥…≥N_{k+1}≥0} x^{ΣN_i} q^{ΣN_i² + N_a+…+N_{k+1}} / ( (q)_{N_1−N_2}⋯(q)_{N_k−N_{k+1}} (q²;q²)_{N_{k+1}} ),   B_{k+1} := B^{(k+2)}_{k+1} (no linear term),
    C_k(x)      := Σ_{0≤s_1≤…≤s_k} x^{Σs_i} q^{Σs_i²} / ∏ (q²;q²)_{s_i−s_{i−1}}      (the s-chain of K_k alone; C_1(x) = (−xq;q²)_∞ by Euler),
    J_k(x)      := (1/(x²q²;q²)_∞) Σ_{n≥0} (−1)^n x^{(k+1)n} q^{(k+1)n²+n(n−1)} (x²;q²)_n (1−x²q^{4n}) / ( (q²;q²)_n (1−x²) ).

*Proved:* J_k(x) = C_k(x)/C_1(x). (Proof: rerun the proof of Proposition 6 in base Q = q² with a = x². The step
(S2) in base Q with ρ_2 = −√(aQ) = −xq reads α'_r = x^r q^{r²} α_r, β'_n = Σ_{m≤n} (−xq;q²)_m x^m q^{m²} β_m/((q²;q²)_{n−m}(−xq;q²)_n);
applying it k+1 times to the unit Bailey pair relative to a = x² in base q² gives β^{(k+1)}_n = (1/(−xq;q²)_n) Σ_{n≥M_1≥…≥M_k≥0} ∏_l x^{M_l} q^{M_l²}
/((q²;q²)_{n−M_1}⋯(q²;q²)_{M_k}) and α^{(k+1)}_r = x^{(k+1)r} q^{(k+1)r²}·(−1)^r q^{r(r−1)} (x²;q²)_r (1−x²q^{4r})/((q²;q²)_r(1−x²)); letting n → ∞ in
the Bailey pair relation, as in Proposition 6, gives C_k(x)/(−xq;q²)_∞ = J_k(x) (the sum over M is C_k(x) with
M_l = s_{k+1−l}), and C_1(x) = (−xq;q²)_∞ by Euler.) *Proved:* J_k(1) = (q^{k+1},q^{k+3},q^{2k+4};q^{2k+4})_∞/(q²;q²)_∞ and
J_k(q) = (q²,q^{2k+2},q^{2k+4};q^{2k+4})_∞/(q²;q²)_∞ (Jacobi triple product, exactly as in the proof of Theorem 3).
*Known:* Bressoud's even-modulus identities (D. Bressoud, J. Combin. Theory A 27 (1979); Mem. AMS 227 (1980); in the
notation of KRTW (3.6), F^{(n)}_{a,0}(1,q) = (q^{a+1},q^{2n−a+1},q^{2n+2};q^{2n+2})_∞/(q)_∞ for 0 ≤ a ≤ n, used here with n = k+1 and
a = k+1, 0, k−1; I did not re-derive them, the three instances are checked numerically in the log) give B_{k+1}(1), B_{k+1}(q) = B^{(1)}_{k+1}(1), B^{(k)}_{k+1}(1) =
(q^{k+2},q^{k+2},q^{2k+4};q^{2k+4})_∞/(q)_∞, (q,q^{2k+3},q^{2k+4};q^{2k+4})_∞/(q)_∞, (q^k,q^{k+4},q^{2k+4};q^{2k+4})_∞/(q)_∞.

**Conjecture (F_k) (numerical).** For all k ≥ 1, in Z[x][[q]]:

    K_k(x) = B_{k+1}(x) · J_k(x),   equivalently   K_k(x) · C_1(x) = B_{k+1}(x) · C_k(x),            and      K^Ω_k(x) = B^{(k)}_{k+1}(x) · J_k(x).

By the three "proved/known" facts, the first statement at x = 1 and x = q is exactly (iii) and (ii), and the second
at x = 1 is exactly (iv). So **(ii), (iii), (iv) for all k follow from (F_k)**; (F_k) itself is NOT proved. Evidence:
`code/krtw_Fk_test.py` → `certificates/krtw_Fk_test.txt`: both statements for k ≤ 7 modulo q^44 in all powers of x (the
log also checks J_k C_1 = C_k and the three specialisations against the products). In the form K_k·C_1 = B_{k+1}·C_k both
sides are multisums with positive terms of the same rank 2k+1 and the same species (k variables of base q, k+1 of
base q²): on the left the s-chain is coupled level by level to the r-chain and one base-q² variable is free, on
the right a single base-q² variable is coupled to all r-levels and the whole s-chain is free. The identity is not
termwise for k ≥ 2. Two natural R-refinements (factor [R;n_1], or (q)_R/(q)_{R−(N_1−N_{k+1})}, on the Bressoud side) are false for
j ≥ 2, so the analogue of Theorem 1 for (F_k) is still to be found.

Numerical side remark (not proved here beyond Theorem 2 itself): Theorem 2 may be specialised at x = 1/q (both
sides are summable after the substitution, the q-order minus the x-degree of a term being ≥ Σ N_i(N_i−1)/2), and
then KRTW_k(1/q) = Mel_k(1/q) = ((−q)_∞/(q)_∞) Σ_{i=1}^{k+1} (q^i,q^{k+2−i},q^{k+2};q^{k+2})_∞ was observed for k ≤ 6 mod q^80 ([X1]; the first
equality is Theorem 2, the second is only numerical). KRTW_k(1/q) is the value at z = 1/q of the second multisum of
Conjecture 4.9; the right side is the sum of the k+1 principally specialised level-k characters of A_1^{(1)}.

## 9. Openness and prior art (searches of 2026-11-26)

arXiv API queries (all:Capparelli AND Meurman AND Primc; authors Warnaar, Tsuchioka, Kanade, Russell, Primc;
abs: "Hall-Littlewood" AND "Rogers-Ramanujan"; abs: "Andrews-Gordon"; abs: Primc AND conjecture; abs: "Meurman-Primc";
"Hall-Littlewood" 2026; "Andrews-Gordon type" 2024–2026; "coloured partitions" AND conjecture; Melzer AND identities),
listings up to 2026-11-24. Papers found that cite or continue KRTW, with what I checked in the full text:

* X. Ding, arXiv:2611.10654 (7 Nov 2026): proves the CMPP conjecture of type C_n^{(1)} (products for 𝓒); quotes the
  Hall–Littlewood form of Griffin–Ono–Warnaar; no multisum, nothing on Conjecture 2.2/4.9.
* L. Bartholdi, arXiv:2611.04020 (2 Nov 2026): three-row level-two case of Primc's conjecture (the k = 2 products).
* Capparelli–Meurman–Primc, arXiv:2611.32425 (24 Nov 2026): three A_2^{(2)} identities (abstract only).
* Bhattacharya–Dai–Hou–Tong, arXiv:2610.41155: cites KRTW p. 11 for a Hall–Littlewood formula; unrelated.
* Warnaar, arXiv:2511.17034; Kurşungöz, arXiv:2605.13096 (k = 1 only); Lau–Ono, arXiv:2608.05480; Russell,
  arXiv:2306.16251: no proof of the displayed identities for k ≥ 3.

So, as far as these searches go, (i) for k ≥ 3 (equivalently Conjecture 2.2 at z = 1 for k ≥ 3) was open; the
usual caveat applies (no full-text search engine, no MathSciNet/zbMATH). The tools used are classical; the new
input is the grading x^{Σ(r_i+2s_i)} (which is not KRTW's z^{Σ(r_i+s_i)}) and the boundary parameter R of Theorem 1.
