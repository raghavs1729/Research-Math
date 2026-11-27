# The KRTW factorisation (F_k) for all k, and the KRTW identities (ii), (iii) for all k

Status (2026-11-27, krtw4 session): **complete proof written; certificate passes; two independent adversarial reads
done** (`../REFEREE-F-all-k.md`, own code in `../referee/`; `../REFEREE-3-all-k.md`, own code in `../referee3/`, written
without seeing the first report or any code): no error, no gap; their exposition items are incorporated below. Identity (iv) is proved in the companion file `iv-all-k.md`, which uses §§1–7.
Certificate: `code/certify_F_all_k.py` → `certificates/certify_F_all_k.txt` (labels [C1]–[C9] below).
Earlier work used here: `../rr-identities/proofs/KRTW-F-moves.md` (Theorem B (B1), (B2) and its §3 paragraph
"From (F_k) to (ii), (iii)"); `../rr-identities/proofs/KRTW-general.md` §8 (the reduction to (F_k)).

## 0. Statements

Source: S. Kanade, M. C. Russell, S. Tsuchioka, S. O. Warnaar, *Remarks on the conjectures of Capparelli, Meurman,
Primc and Primc*, arXiv:2404.03851v2 (Selecta Math. 32 (2026), Paper 100). After Conjecture 4.9 (`Con_Shun2`) they
display four "conjectural Andrews–Gordon-type identities" (verbatim LaTeX in `../LITERATURE-2026-11-27.md` §1b). With
r_{k+1} = s_0 := 0, the sums over all r_1, …, r_k, s_1, …, s_k ≥ 0 and 1/(q;q)_m = 1/(q²;q²)_m := 0 for m < 0, the second and third read

    (ii)   Σ ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+r_i+2s_i} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}})
             = (q,q²,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4}; q^{2k+4})_∞ / ((q²;q²)_∞ (q;q)_∞),
    (iii)  Σ ∏_{i=1}^k q^{(r_i+s_i)²+s_i²} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}})
             = (q^{k+1},q^{k+2},q^{k+2},q^{k+3},q^{2k+4},q^{2k+4}; q^{2k+4})_∞ / ((q²;q²)_∞ (q;q)_∞).

Notation: (a;p)_n as usual, (q)_n := (q;q)_n, Q := q². For k ≥ 1 put

    K_k(x)  := Σ_{r_1≥…≥r_k≥0, 0≤s_1≤…≤s_k} ∏_{i=1}^k x^{r_i+2s_i} q^{(r_i+s_i)²+s_i²} / ((q)_{r_i−r_{i+1}} (Q;Q)_{s_i−s_{i−1}}),
    B_{k+1}(x) := Σ_{N_1≥…≥N_{k+1}≥0} x^{ΣN_i} q^{ΣN_i²} / ((q)_{N_1−N_2}⋯(q)_{N_k−N_{k+1}} (Q;Q)_{N_{k+1}}),
    C_k(x)  := Σ_{0≤s_1≤…≤s_k} x^{Σs_i} q^{Σs_i²} / ∏_i (Q;Q)_{s_i−s_{i−1}},   J_k := C_k/C_1,   C_1(x) = (−xq;Q)_∞.

**Theorem 1 ((F_k) for all k).** For every k ≥ 1, in Z[x][[q]],  K_k(x)·C_1(x) = B_{k+1}(x)·C_k(x), i.e. K_k = B_{k+1}J_k.

**Theorem 2.** The KRTW identities (ii) and (iii) hold for every k ≥ 1.

*Proof of Theorem 2 from Theorem 1.* Exactly `KRTW-F-moves.md` §3, paragraph "From (F_k) to (ii), (iii)": (iii) is
K_k(1), (ii) is K_k(q), and (B1), (B2) of Theorem B there evaluate B_{k+1}(1), B_{k+1}(q), J_k(1), J_k(q) as theta
quotients (Jacobi triple product). ∎

**Corollary 3.** (a) KRTW's Conjecture 2.2 holds at z = q for every k (given Griffin–Ono–Warnaar, arXiv:1401.7718,
Theorem 1.3 with n = 2; `KRTW-F-moves.md` Corollary C(b)). (b) The Hall–Littlewood identity
K_k(x) = Σ_{λ_1≤k} q^{|λ|} P′_{2λ}(√x, −√x; q) holds for every k (given the C₂ Rogers–Selberg identity of
Bartlett–Warnaar; `KRTW-F-moves.md` Corollary C(a)). Both are conditional only on the cited published theorems.

What is new is Theorem 1 for k ≥ 4 (k ≤ 3 was proved by move chains in `KRTW-F-moves.md`, Theorem A, and k = 1 is a
termwise identity) and hence Theorem 2 for k ≥ 4. Identity (iv) is proved in `iv-all-k.md` (§9).

**Idea of the proof.** K_k is a matrix element ⟨ω, T^{k−1}u⟩ of a transfer operator T acting on functions of two
variables (R, S). T is triangular in a Bailey-type basis with eigenvalues λ_{N,M} = x^{N+M}q^{N²+M²} (N ≥ M ≥ 0).
Instead of computing its eigenfunctions (which are not hypergeometric terms), we use two *local* operators H′, G that
commute with T (they were found numerically in `KRTW-F-moves.md` §6; here they are proved to commute). Restricted to the
lowest S-slice, the two eigen-equations, read "at R = ∞", give a closed one-step recursion for the value of each
eigenfunction at infinity (§4). The same works for the reversed chain, whose eigenfunctions give the left eigenvectors
(§5). The left and right eigenvectors of the same mode meet in a single point, so the spectral expansion of K_k is
explicit (§6) and equals the known double-sum form of B_{k+1}J_k (§7).

## 1. Setting

**Ring.** 𝔽 := Q(q)((x)), the field of formal Laurent series in x over Q(q), with the x-adic valuation ord. Rational
functions of x, q are elements of 𝔽 (expand in x). For each D the coefficient of x^D in K_k, B_{k+1}, C_k is a finite sum
of rational functions of q regular at q = 0; so Theorem 1 in 𝔽 implies Theorem 1 in Q[[q]][[x]], hence
in Z[x][[q]] (all series involved lie in Z[x][[q]], `KRTW-F-moves.md` §3). (C_1 ≡ 1 mod x is used only to define
J_k = C_k/C_1.) **Infinite products** such as (xq;q)_∞, (x²Q;Q)_∞, 1/(aq;q)_∞ (a = xq^{2s}) do not converge x-adically in
𝔽 (Q(q) carries the discrete topology); they denote the elements of 𝔽 given by their Euler expansions, e.g.
(xq;q)_∞ = Σ_n (−1)^n q^{n(n+1)/2} x^n/(q)_n and 1/(zq;q)_∞ = Σ_n z^n q^n/(q)_n, whose x-coefficients are rational functions
of q regular at q = 0. Every identity used between such products and the series above (e.g. (B-ii), the product
manipulations of §§4–7, (B1), (B2)) holds in 𝔽 by the same transfer: the coefficientwise Taylor map from series with
such coefficients to Q[[q]][[x]] is an injective ring homomorphism, and the identity holds there.

**Functions.** 𝒱 := {F : N² → 𝔽 with inf ord F > −∞}; F(R,S) := 0 when R < 0 or S < 0. Put w(r,s) := x^{r+2s}q^{(r+s)²+s²}
(note ord w(r,s) = r + 2s) and

    (TF)(R,S) := Σ_{r=0}^{R} Σ_{s≥S} w(r,s) F(r,s) / ((q)_{R−r} (Q;Q)_{s−S}),       u(R,S) := 1/(q)_R,
    ω(F)      := Σ_{r,s≥0} w(r,s) F(r,s) / (Q;Q)_s.

The sums converge x-adically for F ∈ 𝒱, and T𝒱 ⊂ 𝒱. Unwinding the definitions (level i of K_k is the pair (r_i, s_i);
the factor 1/(q)_{r_k} is u, the factor 1/(Q;Q)_{s_1} is in ω):

    K_k(x) = ω(T^{k−1}u)     (k ≥ 1).                                                                         (1.1)

**Operators.** P, Y: multiplication by q^R, q^{2S}; 𝖲F(R,S) := F(R−1,S); 𝔲F(R,S) := F(R,S+1); W: multiplication by w.
Then 𝖲P = q^{−1}P𝖲, 𝔲Y = QY𝔲, all other pairs among P, Y, 𝖲, 𝔲 commute, and

    T = E·W,   E := Σ_{j,b≥0} 𝖲^j 𝔲^b / ((q)_j (Q;Q)_b),                                                   (1.2)
    W𝖲 = xq^{−1}P²Y 𝖲 W,   W𝔲 = x^{−2}q^{−2}P^{−2}Y^{−2} 𝔲 W                                                 (1.3)

((1.3): w(R,S)/w(R−1,S) = xq^{2R+2S−1}, w(R,S)/w(R,S+1) = x^{−2}q^{−2R−4S−2}).

## 2. Two commuting local operators

    H′ := (Y^{−1}+x)P^{−1}(1−𝖲) + xP(1+xY) − x²PY𝔲,
    G  := P^{−2}Y^{−1}(1−𝖲)(1−q𝖲) + Y^{−1} − x(1+q^{−1})𝖲 + x²(1+P²)Y(1−𝔲).

**Proposition 2.1.** TH′ = H′T and TG = GT on 𝒱.

*Proof.* Every operator below is a series Σ_{j,b} c_{j,b}(P,Y)𝖲^j𝔲^b with Laurent-polynomial coefficients written to the
left ("normal order"). They are applied to the class 𝒞 of functions on Z² supported in {R ≥ 0, S ≥ −1} with
ord F(R,S) ≥ R + 2S + c for some constant c (it contains WF and H′WF for F ∈ 𝒱); on 𝒞 such a series converges, because
F(R,S+b) has ord ≥ 2b + R + 2S + c, and all rearrangements below are finite recombinations of x-adically convergent sums.
Products are normal-ordered with 𝖲^j P^m = q^{−jm}P^m𝖲^j, 𝔲^b Y^n = Q^{bn}Y^n𝔲^b, and equal normal-ordered coefficients
give equal operators on W𝒱 (this is all that is used). Restriction to N²: we compute with functions on Z² vanishing
outside N²; 𝔲, E and H′ can produce values off N² (e.g. (H′F)(R,−1) = −x²q^{R−2}F(R,0)), but the restriction of EW to N² is
T, it reads only R′ ≤ R and S′ ≥ S, and EWF vanishes for R < 0; so the identity EWH′ = H′EW on Z², restricted to N², is
TH′ = H′T on 𝒱 (likewise for G and for the operators of `iv-all-k.md` §2).
By (1.3), W H′ = H̃′ W with

    H̃′ = (Y^{−1}+x)P^{−1} − xq^{−1}P(1+xY)𝖲 + xP(1+xY) − q^{−2}P^{−1}Y^{−1}𝔲.

It suffices to show H′E = EH̃′ (then TH′ = EWH′ = EH̃′W = H′EW = H′T). The coefficient of 𝖲^j𝔲^b, multiplied by
(q)_j(Q;Q)_b, is, using 1/(q)_j − 1/(q)_{j−1} = q^j/(q)_j and 1/(Q;Q)_{b−1} = (1−Q^b)/(Q;Q)_b,

    in H′E:  (Y^{−1}+x)q^jP^{−1} + xP(1+xY) − x²PY(1−Q^b),
    in EH̃′:  (Q^{−b}Y^{−1}+x)q^jP^{−1} − xq^{−j}(1−q^j)P(1+xQ^bY) + xq^{−j}P(1+xQ^bY) − q^jQ^{−b}(1−Q^b)P^{−1}Y^{−1},

and both equal q^jP^{−1}Y^{−1} + xq^jP^{−1} + xP + x²Q^bPY. For G: WG = G̃W with

    G̃ = P^{−2}Y^{−1} + Y^{−1} − x(1+q^{−1})𝖲 − x²q^{−1}(1+q^{−1})P²Y𝖲 + x²q^{−3}P²Y𝖲² + x²(1+P²)Y − q^{−2}(1+P^{−2})Y^{−1}𝔲,

and the coefficients of 𝖲^j𝔲^b (times (q)_j(Q;Q)_b) of GE and EG̃ both equal
q^{2j}P^{−2}Y^{−1} + Y^{−1} − x(1+q^{−1})(1−q^j) + x²Q^b(1+P²)Y (for GE use 1 − (1+q)(1−q^j) + q(1−q^j)(1−q^{j−1}) = q^{2j}). ∎
[C1a: commutators vanish on point masses, exact arithmetic at two rational points.]

Rewriting (algebra): with the rank-one operators h_a := P^{−1}(1−𝖲) + aP (acting in R, a possibly depending on S),

    H′ = (Y^{−1}+x) h_{xY} − x² PY𝔲,       G = Y^{−1}h_{xY}² + Y^{−1}(1−xY)² − x²(1+P²)Y𝔲.                      (2.1)

## 3. The Bailey basis and the eigenvectors of T

For a with all (aq)_n invertible, φ_e(R;a) := 1/((q)_{R−e}(aq)_{R+e}) (= 0 for R < e). Standard facts:

* (B-i) Σ_{r≤R} a^rq^{r²}φ_e(r;a)/(q)_{R−r} = a^eq^{e²}φ_e(R;a) (Bailey's lemma, limiting case (L∞)).
* (B-ii) Σ_{r≥0} a^rq^{r²}φ_e(r;a) = a^eq^{e²}/(aq;q)_∞ (Cauchy: Σ_j z^jq^{j²}/((q)_j(zq)_j) = 1/(zq)_∞, z = aq^{2e}).
* (B-iii) h_aφ_e(·;a) = (q^{−e}+aq^e)φ_e(·;a) (from φ_e(R−1;a)/φ_e(R;a) = (1−q^{R−e})(1−aq^{R+e})).
* (B-iv) Change of a: φ_d(·;aq^c) = Σ_{e=d}^{d+c} Γ^{(c)}_{d,e}(a)φ_e(·;a) with
  Γ^{(c)}_{d,e} = ((1−aq^{2e})/(1−a))((a)_{d+e}/(aq^{c+1})_{d+e})(−1)^{e−d}q^{binom(e−d,2)}a^{e−d}q^{(e+d)(e−d)}[c; e−d]_q
  (Bailey inversion and q-Chu–Vandermonde; `code/tower.py`), and Σ_e Γ^{(c)}_{d,e} = (aq;q)_c. Support and sum are
  algebraic (partial fractions in z = q^R): over the common denominator (q)_{R−d}(aq)_{R+d+c}, the numerator of φ_d(·;aq^c)
  is the constant (aq)_c, and the numerators p_e(z) = (q^{R−e+1})_{e−d}(aq^{R+e+1})_{d+c−e} of φ_e(·;a), d ≤ e ≤ d+c, are c+1
  polynomials of degree ≤ c, linearly independent (p_e(q^m) = 0 for d ≤ m < e, p_e(q^e) ≠ 0), hence a basis; all p_e(0) = 1, so
  the coefficient sum is the value (aq)_c at z = 0. Only the support, the coefficient sum and the case c = 1 are used
  (the general formula was also checked independently by the referees). For c = 1: with z = q^R,
  1 − aq = A(1 − aq^{d+1}z) + B(1 − q^{−d}z) gives Γ^{(1)}_{d,d} = A = (1−aq)/(1−aq^{2d+1}) and Γ^{(1)}_{d,d+1} = B = −aq^{2d+1}A.
  Also φ_d(R;aq²) = (aq)_2 φ_{d+1}(R+1;a).
* (B-v) For n ≥ 1, with τ := multiplication by q^R after the shift R → R+1, i.e. (τF)(R) = q^RF(R+1):
  τφ_n = α_nφ_{n−1} + β_nφ_n + γ_nφ_{n+1},
  α_n = q^{n−1}/((1−aq^{2n−1})(1−aq^{2n})), β_n = −q^{n−1}(1+aq^{2n})/((1−aq^{2n−1})(1−aq^{2n+1})), γ_n = aq^{3n}/((1−aq^{2n})(1−aq^{2n+1})),
  and α_n+β_n+γ_n = 0; likewise φ_n(R+1) and q^{2R}φ_n(R+1) are combinations of φ_{n−1}, φ_n, φ_{n+1} with coefficient
  sums 1 and 0. [C2: symbolic in q^n, q^R.] (These are partial fractions in z = q^R; n ≥ 1 is needed.)

For a finite combination F = Σ_e κ_e φ_e(·;a) put v_a(F) := Σ_e κ_e (for |q|, |a| < 1 this is lim_{R→∞}(q)_∞(aq)_∞F(R)).
Then: v_a(h_aF) = Σ_e(q^{−e}+aq^e)κ_e; by (B-v), if κ_0 = 0 then v_a(τF) = 0, v_a(F(·+1)) = v_a(F), v_a(q^{2R}F(·+1)) = 0, and
h_a(τF) = F(·+1) + aq^{2R}F(·+1) − q^{−1}F (a pointwise identity for R ≥ 0, using F(0) = 0), so v_a(h_aτF) = (1−q^{−1})v_a(F);
by (B-iv), v_a(φ_d(·;aq^{c})) = (aq;q)_c.

**Basis of 𝒱.** Ξ_{e,σ}(R,S) := δ_{S,σ}φ_e(R;xq^{2σ}). Every F ∈ 𝒱 is uniquely Σ_{e,σ} f_{e,σ}Ξ_{e,σ} (Bailey inversion in each
slice; locally finite). By (B-i), w(r,σ) = x^{2σ}q^{2σ²}(xq^{2σ})^rq^{r²} and (B-iv),

    T Ξ_{e,σ} = w(e,σ) Σ_{S=0}^{σ} (Q;Q)_{σ−S}^{−1} Σ_{e′=e}^{e+2(σ−S)} Γ^{(2σ−2S)}_{e,e′}(xq^{2S}) Ξ_{e′,S}.                    (3.1)

With (N, M) := (e+σ, σ) the diagonal entry is w(e,σ) = λ_{N,M} := x^{N+M}q^{N²+M²}; (e,σ) ↦ (N,M) is a bijection onto
{N ≥ M ≥ 0}, and the λ_{N,M} are pairwise distinct monomials. Write L := N − M.

**Proposition 3.1.** For N ≥ M ≥ 0 there is a unique φ_{N,M} ∈ 𝒱 with Tφ_{N,M} = λ_{N,M}φ_{N,M} of the form
Ξ_{L,M} + Σ_{S<M} Σ_e c^{(S)}_e Ξ_{e,S}; in it c^{(S)}_e ≠ 0 only for L ≤ e ≤ L + 2(M−S). Every T-eigenvector with finitely
many nonzero Ξ-coefficients is a multiple of some φ_{N,M}. φ_{N,M}(R,S) = 0 unless R ≥ L and S ≤ M.

*Proof.* By (3.1) the eigen-equation is triangular in S (downwards): for S < M,
(λ_{N,M} − λ_{e+S,S}) c^{(S)}_e = Σ_{σ>S} [known], and λ_{e+S,S} ≠ λ_{N,M} for S < M. The support follows from that of Γ.
For an eigenvector F with finitely many Ξ-coefficients look at its highest slice σ: there (λ − w(e,σ))f_{e,σ} = 0, so
(e,σ) = (L,M); hence F/f_{L,M} = Ξ_{L,M} + (lower slices), and by the uniqueness part F = f_{L,M}φ_{N,M}. ∎  [C3 confirms, from an independent computation of φ_{N,M} by this triangular system (`code/tower.py`).]

**Shift.** Since w(r,s+1) = x²q²·w|_{x→xq²}(r,s), the function (R,S) ↦ φ_{N,M}(R,S+1;x) is a T|_{x→xq²}-eigenvector with
eigenvalue λ_{N,M}(x)/(x²q²) = λ_{N−1,M−1}(xq²) and top Ξ_{L,M−1}|_{x→xq²}; by uniqueness

    φ_{N,M}(·, S+1; x) = φ_{N−1,M−1}(·, S; xq²)      (M ≥ 1).                                               (3.2)

**Lemma 3.2.** H′φ_{N,M} = ε_{N,M}φ_{N,M} and Gφ_{N,M} = g_{N,M}φ_{N,M}, where

    ε_{N,M} = q^{−N−M}(1+xq^{2N})(1+xq^{2M}),     g_{N,M} = q^{−2N}+x²q^{2N}+q^{−2M}+x²q^{2M}.

*Proof.* By (2.1), (B-iii), (B-iv) and (B-v) (note q^Rφ_e(R;xq^{2S+2}) = (xq^{2S+1})_2 (τφ_{e+1})(R) relative to xq^{2S}, and
e + 1 ≥ 1), H′ and G map functions with finitely many Ξ-coefficients to such functions and do not raise the top slice.
By Proposition 2.1, H′φ_{N,M} is a T-eigenvector for λ_{N,M}, hence a multiple of φ_{N,M}; the multiple is read off on the
slice S = M, where the 𝔲-terms vanish: (q^{−2M}+x)(q^{−L}+xq^{2M+L}) = ε_{N,M}, and for G
q^{−2M}(q^{−L}+xq^{2M+L})² + q^{−2M}(1−xq^{2M})² = g_{N,M}. ∎ [C3]

## 4. Values at infinity of the eigenvectors of T

Put W_{N,M}(x) := v_x(φ_{N,M}(·,0)).

**Proposition 4.1.** W_{L,0} = 1 and, for M ≥ 1,

    W_{N,M}(x) = x²(1+x)(x+q^{−1}) q^{2N+2M} (xq;q)_2 / ((1−q^{2N})(1−q^{2M})(1−x²q^{2N})(1−x²q^{2M})) · W_{N−1,M−1}(xq²).

*Proof.* φ_{L,0} = Ξ_{L,0}. Let M ≥ 1, F_0 := φ_{N,M}(·,0) (a combination of φ_e(·;x), e ≥ L) and F_1 := φ_{N,M}(·,1). By (3.2)
F_1 = Σ_e κ′_e φ_e(·;xq²) with Σκ′_e = W_{N−1,M−1}(xq²). By (B-iv), F_1(R) = (xq)_2 G(R+1) with G := Σ_e κ′_e φ_{e+1}(·;x), a
combination of φ_n(·;x), n ≥ 1, with v_x(G) = Σκ′_e. Put v_1 := (xq)_2 W_{N−1,M−1}(xq²). By the rules of §3 (applied to G):
q^RF_1(R) = (xq)_2(τG)(R) has v_x = 0, h_x(q^RF_1) has v_x = (1−q^{−1})v_1, and (1+q^{2R})F_1(R) = (xq)_2(1+q^{2R})G(R+1) has
v_x = v_1. The slices S = 0 of H′φ = εφ and Gφ = gφ read, by (2.1) (Y = 1 on S = 0),

    (1+x) h_x F_0 − x² q^R F_1(R) = ε F_0,       h_x² F_0 + (1−x)² F_0 − x²(1+q^{2R}) F_1(R) = g F_0.

Apply v_x to both equations and to h_x applied to the first; with μ_j := v_x(h_x^jF_0) (h_x is diagonal on F_0):

    (1+x)μ_1 − εμ_0 = 0,      (1+x)μ_2 − εμ_1 = x²(1−q^{−1}) v_1,      μ_2 + (1−x)²μ_0 − x² v_1 = g μ_0.

Eliminating μ_1, μ_2: μ_0 [ε² + (1+x)²((1−x)² − g)] = x²(1+x)(x+q^{−1}) v_1. With e_n := q^{−n}+xq^n one has ε = e_Ne_M,
g = e_N²+e_M²−4x, so the bracket is (e_N²−(1+x)²)(e_M²−(1+x)²), and e_n²−(1+x)² = q^{−2n}(1−q^{2n})(1−x²q^{2n}). ∎ [C4b]

**Corollary 4.2.** W_{N,M}(x) = (xq;q)_{2M} V_{N,M}(x) with

    V_{N,M}(x) := x^{2M} q^{2MN+M(2M−1)} (−x;q)_M / ((Q;Q)_M (Q^{N−M+1};Q)_M (xq^M;q)_M (x²q^{2N};Q)_M).

(Iterate Proposition 4.1 and use (x²q^{2M};Q)_M = (xq^M;q)_M(−xq^M;q)_M.) [C4a: compared with the eigenvectors of §3, N ≤ 7.]

**Lemma 4.3.** ω(φ_{N,M}) = λ_{N,M} W_{N,M} / (xq;q)_∞.
*Proof.* For the slice S = s, φ_{N,M}(·,s) = Σ_e c^{(s)}_e φ_e(·;xq^{2s}). By (B-ii), Σ_r w(r,s)φ_e(r;xq^{2s}) = w(e,s)/(xq^{2s+1})_∞,
so ω(φ) = (xq)_∞^{−1} Σ_s (Q;Q)_s^{−1} Σ_e c^{(s)}_e w(e,s)(xq)_{2s}. By (B-i), (Tφ)(·,0) = Σ_s (Q;Q)_s^{−1} Σ_e c^{(s)}_e w(e,s)φ_e(·;xq^{2s}),
whose v_x is the same double sum by (B-iv). As (Tφ)(·,0) = λφ(·,0), ω(φ) = λ v_x(φ(·,0))/(xq)_∞. ∎

## 5. The reversed chain and the left eigenvectors

    (T♯G)(r,s) := Σ_{R≥r} Σ_{S≤s} w(R,S) G(R,S) / ((q)_{R−r} (Q;Q)_{s−S})     (T♯ = W^{−1}T^tW as kernels [C1c]).

If T♯G = λG then ψ_G(F) := Σ_{r,s} w(r,s)G(r,s)F(r,s) satisfies ψ_G(TF) = λψ_G(F) whenever the double sum is finite (the
case below). Let φ^Q_e(S;a) := 1/((Q;Q)_{S−e}(aQ;Q)_{S+e}) and Ξ♯_{e,ρ}(R,S) := δ_{R,ρ}φ^Q_e(S;x²q^{2ρ}). Since
w(ρ,S) = x^ρq^{ρ²}(x²q^{2ρ})^SQ^{S²}, the base-Q versions of (B-i), (B-iv) give

    T♯Ξ♯_{e,ρ} = λ_{ρ+e,e} Σ_{r=0}^{ρ} (q)_{ρ−r}^{−1} Σ_{e′=e}^{e+ρ−r} Γ^{Q,(ρ−r)}_{e,e′}(x²q^{2r}) Ξ♯_{e′,r},

where Γ^{Q,(c)} is Γ^{(c)} with q → Q (change a → aQ^c in base Q; Σ_{e′}Γ^{Q,(c)}_{e,e′}(a) = (aQ;Q)_c). Exactly as in
Proposition 3.1 there is a unique T♯-eigenvector φ♯_{N,M} = Ξ♯_{M,L} + Σ_{ρ<L} Σ_e c_{e,ρ}Ξ♯_{e,ρ} for λ_{N,M}, with c_{e,ρ} ≠ 0 only for
M ≤ e ≤ M + L − ρ; φ♯_{N,M}(R,S) = 0 unless R ≤ L and S ≥ M; and, by the same top-slice argument and uniqueness, every
T♯-eigenvector with finitely many Ξ♯-coefficients is a multiple of some φ♯_{N,M}. [`code/sharp.py`, C5a]

Transposing Proposition 2.1, H′♯ := W^{−1}H′^tW and G♯ := W^{−1}G^tW commute with T♯ (the kernel identity
(T♯H′♯)(a,b) = w(a)^{−1}w(b)(H′T)(b,a) involves only finite sums, and applied to φ♯ all sums are finite: H′♯ is local, T♯
sums over finitely many S, and φ♯ vanishes for R > L); explicitly, with
𝖲♮F(R,S) := F(R+1,S), 𝔲♮F(R,S) := F(R,S−1) and h^Q_a := Y^{−1}(1−𝔲♮) + aY (acting in S),

    H′♯ = P^{−1}Y^{−1}(1−𝔲♮) + xP^{−1} + xP(1+xY)(1−𝖲♮),
    G♯  = (P^{−2}+1) h^Q_{x²P²} − x(1+q^{−1})𝖲♮ − x²(1+q)P²Y𝖲♮ + q x²P²Y(𝖲♮)²                                (5.1)

(e.g. H′♯(a,b) = w(b)H′(b,a)/w(a) entrywise; [C1b]). Base-Q facts used: h^Q_aφ^Q_e(·;a) = (Q^{−e}+aQ^e)φ^Q_e(·;a);
φ^Q_e(S;aQ) = c₁φ^Q_e(S;a) + c₂φ^Q_{e+1}(S;a) with c₁ = (1−aQ)/(1−aQ^{2e+1}), c₂ = −aQ^{2e+1}c₁, c₁+c₂ = 1−aQ (the case c = 1 of
(B-iv) in base Q); Q^Sφ^Q_e(S;aQ) = b_e(φ^Q_e(S;a) − φ^Q_{e+1}(S;a)), b_e = (1−aQ)/(Q^{−e}−aQ^{e+1});
Q^Sφ^Q_e(S;aQ²) is a combination of φ^Q_e, φ^Q_{e+1}, φ^Q_{e+2} (relative to a) with coefficient sum 0 (partial fractions in
Q^S; for |q|, |a| small, coefficient sums are values at S → ∞ of (Q;Q)_∞(aQ;Q)_∞·function). With v^Q_a defined as v_a
(coefficient sum relative to a), these give as in Lemma 3.2: H′♯φ♯_{N,M} = ε_{N,M}φ♯_{N,M}, G♯φ♯_{N,M} = g_{N,M}φ♯_{N,M} [C5a]
(on the top slice ρ = L the 𝖲♮-terms vanish and (5.1) reduces to q^{−L}h^Q + x(q^{−L}+q^L), resp. (q^{−2L}+1)h^Q).

**Proposition 5.1.** Let v_ρ := Σ_e c_{e,ρ} (coefficient sum of φ♯_{N,M}(ρ,·) relative to a_ρ := x²q^{2ρ}); v_L = 1. For ρ < L,

    v_ρ = −x q^{L+N+M−1} (1−x²q^{2ρ+2})(1−q^{2ρ+1}) / ((1−q^{L−ρ})(1−q^{L+ρ})(1−xq^{N+M−ρ})(1−xq^{N+M+ρ})) · v_{ρ+1}.

*Proof.* Write G_ρ := φ♯(ρ,·), a := a_ρ (so a_{ρ+1} = aQ), μ_j := v^Q_a((h^Q_a)^jG_ρ). The slice ρ of (5.1) gives

    [q^{−ρ}h^Q_a + x(q^{−ρ}+q^ρ) − ε] G_ρ = xq^ρ(1+xY) G_{ρ+1},
    (q^{−2ρ}+1) h^Q_a G_ρ − x(1+q^{−1})G_{ρ+1} − x²(1+q)q^{2ρ}Y G_{ρ+1} + q x²q^{2ρ} Y G_{ρ+2} = g G_ρ.

Applying v^Q_a (the Y-terms have coefficient sum 0, and v^Q_a(G_{ρ+1}) = (1−aQ)v_{ρ+1}):
q^{−ρ}μ_1 + (x(q^{−ρ}+q^ρ) − ε)μ_0 = xq^ρ(1−aQ)v_{ρ+1} and (q^{−2ρ}+1)μ_1 − x(1+q^{−1})(1−aQ)v_{ρ+1} = gμ_0. Eliminating μ_1 with
c := q^ρ+q^{−ρ} (note q^{−ρ}/(q^{−2ρ}+1) = 1/c and q^ρc = 1+q^{2ρ}):

    μ_0 [g + xc² − εc] = x(1−aQ)(q^{2ρ} − q^{−1}) v_{ρ+1}.

With e_n = √x ζ_n, ζ_n := z_n + z_n^{−1}, z_n := √x q^n, the bracket is x(c² − ζ_Nζ_Mc + ζ_N² + ζ_M² − 4) =
x(c − z_Nz_M − (z_Nz_M)^{−1})(c − z_N/z_M − z_M/z_N) (the Fricke relation), i.e.
q^{−L−N−M}(1−q^{L−ρ})(1−q^{L+ρ})(1−xq^{N+M−ρ})(1−xq^{N+M+ρ}). ∎ [C5c]

**Corollary 5.2.** V♯_{N,M} := v_0 = (−x)^L q^{L(2L−1)+2ML} · ((q;q²)_L/((q)_L(q^L;q)_L)) · (x²Q;Q)_L / ((xq^{2M+1};q)_L (xq^{2M+L};q)_L);
for L ≥ 1, (q;q²)_L/((q)_L(q^L;q)_L) = (1+q^L)/(Q;Q)_L; V♯_{M,M} = 1. [C5b, N ≤ 7]

**Lemma 5.3.** ψ_{N,M}(u) := Σ_{r,s} w(r,s)φ♯_{N,M}(r,s)/(q)_r = λ_{N,M} V♯_{N,M} / (x²Q;Q)_∞.
*Proof.* As Lemma 4.3, with the base-Q forms of (B-ii) (Σ_s w(r,s)φ^Q_e(s;a_r) = w(r,e)/(a_rQ;Q)_∞ = w(r,e)(x²Q;Q)_r/(x²Q;Q)_∞),
(B-i) for (T♯φ♯)(0,·) = λφ♯(0,·), and v^Q_{x²}(φ^Q_e(·;x²Q^r)) = (x²Q;Q)_r. ∎

## 6. Biorthogonality, completeness and the spectral expansion

Let ψ_{N,M}(F) := Σ_{r,s} w(r,s)φ♯_{N,M}(r,s)F(r,s) (r ≤ L, s ≥ M; converges for F ∈ 𝒱).

**Lemma 6.1.** ψ_a(φ_b) = 0 for modes a ≠ b, and ψ_{N,M}(φ_{N,M}) = λ_{N,M} / ((x²q^{2L+2};Q)_{2M} (xq^{2M+1};q)_{2L}).
*Proof.* For a ≠ b the sum defining ψ_a(Tφ_b) is finite (r ≤ L_a, s ≤ M_b), so it can be computed as λ_bψ_a(φ_b) and as
λ_aψ_a(φ_b); the λ's differ. By the supports (φ_{N,M}: R ≥ L, S ≤ M; φ♯_{N,M}: R ≤ L, S ≥ M) only (r,s) = (L,M) contributes to
ψ_{N,M}(φ_{N,M}), where w = λ, φ_{N,M}(L,M) = φ_L(L;xq^{2M}) = 1/(xq^{2M+1})_{2L}, φ♯_{N,M}(L,M) = 1/(x²q^{2L}Q;Q)_{2M}. ∎ [C6]

**Lemma 6.2 (valuations).** (a) ord c^{(S)}_e ≥ 2(M−S) for the coefficients of φ_{N,M}; hence ord φ_{N,M}(R,S) ≥ 2(M−S).
(b) ord c_{e,ρ} ≥ 0 for φ♯_{N,M}; hence ord φ♯_{N,M}(R,S) ≥ 0. (c) ord ψ_{L′+M′,M′}(Ξ_{e,σ}) ≥ e + 2σ, with equality and
ψ ≠ 0 for (e,σ) = (L′,M′); and ψ_{L′+M′,M′}(Ξ_{e,σ}) = 0 unless e ≤ L′, σ ≥ M′.
*Proof.* (a) By Lemma 3.2 and (2.1), slice S: [ε_{e+S,S} − ε_{N,M}]c^{(S)}_e = x²q^{2S}(xq^{2S+1})_2·(the φ_e-coefficient of
Σ_{e′}c^{(S+1)}_{e′}τφ_{e′+1}, relative to xq^{2S}). The bracket has ord 0 unless e + 2S = N+M (the largest index on slice S),
where it has ord 1 (its x-coefficient is q^{N+M−2S}+q^{2S−N−M}−q^{L}−q^{−L} ≠ 0); that index is reached only through γ, and
ord γ ≥ 1 (γ contains a = xq^{2S}). Induction downwards from S = M. (b) Same with the first equation of the proof of
Proposition 5.1: the right side carries a factor x, the bracket ε_{ρ+e,e} − ε_{N,M} is nonzero for ρ < L with ord ≤ 1 (if
its constant terms q^{−ρ−2e}, q^{−N−M} agree, its x-coefficient is q^ρ+q^{−ρ}−q^L−q^{−L} ≠ 0), and the re-expansion
coefficients c₁, c₂, b_e have ord ≥ 0. (c) ψ(Ξ_{e,σ}) = Σ_{e≤r≤L′} w(r,σ)φ♯(r,σ)φ_e(r;xq^{2σ}) with ord w(r,σ) = r+2σ; at
(e,σ) = (L′,M′) only r = L′ occurs and the term is the explicit monomial-times-unit of Lemma 6.1. ∎

**Lemma 6.3 (completeness).** If Δ ∈ 𝒱 and ψ_{N,M}(Δ) = 0 for all N ≥ M ≥ 0, then Δ = 0.
*Proof.* Expand Δ = Σ δ_{e,σ}Ξ_{e,σ}; the δ's have ord ≥ inf ord Δ (Bailey inversion has x-integral coefficients). If Δ ≠ 0
let m := min{ord δ_{e,σ} + e + 2σ}, attained on a finite nonempty set Z. Choose (e*,σ*) ∈ Z with e* minimal and then σ*
maximal, and evaluate ψ_{e*+σ*,σ*}(Δ) = Σ_{e≤e*,σ≥σ*} ψ(Ξ_{e,σ})δ_{e,σ} (an x-adically convergent sum by Lemma 6.2(c)). The term
(e*,σ*) has ord exactly m; a term with e < e*, or e = e* and σ > σ*, has ord ≥ ord δ_{e,σ} + e + 2σ > m (it is not in Z). So
ψ_{e*+σ*,σ*}(Δ) ≠ 0. ∎

**Theorem 6.4.** In 𝒱 (pointwise, x-adically convergent),

    u = Σ_{N≥M≥0} ũ_{N,M} φ_{N,M},     ũ_{N,M} := ψ_{N,M}(u)/ψ_{N,M}(φ_{N,M}) = V♯_{N,M}(x²q^{2L+2};Q)_{2M}(xq^{2M+1};q)_{2L}/(x²Q;Q)_∞,

and for every k ≥ 1:   K_k(x) = Σ_{N≥M≥0} λ_{N,M}^k · V_{N,M}V♯_{N,M}(x²q^{2L+2};Q)_{2M}(xq;q)_{2N} / ((xq;q)_∞(x²Q;Q)_∞).

*Proof.* ord ũ ≥ 0 (Corollary 5.2), so by Lemma 6.2(a) the (N,M)-term of Ũ := Σũφ has ord ≥ 2(M−S) at (R,S), and only
L ≤ R, M ≥ S contribute: Ũ ∈ 𝒱 converges pointwise. The family (w(r,s)φ♯_a(r,s)ũ_bφ_b(r,s)) is x-adically summable, so
ψ_a(Ũ) = Σ_bũ_bψ_a(φ_b) = ψ_a(u) by Lemma 6.1; Lemma 6.3 gives u = Ũ (ũ_{N,M} by Lemmas 5.3, 6.1). Summability, explicitly:
in the families (r,s,b) ↦ w(r,s)φ♯_a(r,s)ũ_bφ_b(r,s), w(r,s)ũ_bλ_b^jφ_b(r,s)/((q)_{R−r}(Q;Q)_{s−S}) and
w(r,s)ũ_bλ_b^jφ_b(r,s)/(Q;Q)_s (any j ≥ 0), a nonzero term has r ≥ L_b, s ≤ M_b and ord ≥ r + 2M_b + j(N_b+M_b) (Lemma 6.2(a), ord ũ ≥ 0,
ord φ♯ ≥ 0), so only finitely many terms lie below any given order. This gives T^{k−1}u = Σũλ^{k−1}φ (inductively in k) and,
by (1.1) and Lemma 4.3, K_k = Σ ũ_{N,M}λ^k W_{N,M}/(xq)_∞; insert W = (xq)_{2M}V and
(xq)_{2M}(xq^{2M+1})_{2L} = (xq)_{2N}. ∎ [C9: the expansion of u checked at (R,S) ≤ (2,2) mod x^7.]

## 7. The product identity; proof of Theorem 1

Let β^B_n := (−1)^nq^{n²}(x²;Q)_n(1−xq^{2n})/((Q;Q)_n(1−x)), β^J_m := (−1)^mq^{m²−m}(x²;Q)_m(1−x²q^{4m})/((Q;Q)_m(1−x²)) and
t_{N,M} := λ_{N,M}(β^B_Nβ^J_M + β^B_Mβ^J_N)/(1+δ_{NM}). By `KRTW-F-moves.md` (B4) (the identity for E there),

    t_{N,M} = λ_{N,M}(−1)^{N+M} q^{N²+M²−N−M}(q^N+q^M)(1+xq^{N+M})(1−xq^{2N})(1−xq^{2M})(x²;Q)_N(x²;Q)_M / ((Q;Q)_N(Q;Q)_M(1−x)(1−x²)(1+δ_{NM})).

**Lemma 7.1.** V_{N,M}V♯_{N,M}(x²q^{2L+2};Q)_{2M}(xq;q)_{2N} = t_{N,M}.
*Proof.* Collect factors (N = L+M): Q-factors (Q;Q)_M(Q^{L+1};Q)_M(Q;Q)_L = (Q;Q)_M(Q;Q)_N; x²-factors
(x²Q;Q)_L(x²q^{2L+2};Q)_{2M}/(x²q^{2N};Q)_M = (x²Q;Q)_{N−1}(1−x²q^{2N+2M}); the remaining x-factors combine, via
(x;q)_M(xq^M;q)_M = (x;q)_{2M}, (x;q)_{2M}(1−xq^{2M})(xq^{2M+1};q)_L = (x;q)_{N+M+1}, (x²;Q)_M = (x;q)_M(−x;q)_M and
(x;q)_{2N+1} = (x;q)_{N+M}(1−xq^{N+M})(xq^{N+M+1};q)_{L−1}(1−xq^{2N}) (L ≥ 1), to the x-part of t_{N,M}; the monomials agree
(x^{N+M}(−1)^{N+M}q^{2L²+4LM+4M²−L−M} on both sides, using q^N+q^M = q^M(1+q^L)); for L = 0 the factor (q;q²)_0/… = 1 matches
(q^N+q^M)/(1+δ_{NM}) = q^M. ∎ [C7: symbolic for N ≤ 8.]

*Proof of Theorem 1.* By Theorem 6.4 and Lemma 7.1, (xq)_∞(x²Q;Q)_∞K_k = Σ_{N≥M}λ_{N,M}^k t_{N,M} = Σ_{n,m≥0}(λ_nλ_m)^{k+1}β^B_nβ^J_m
with λ_n := x^nq^{n²} (symmetrisation; λ_{N,M} = λ_Nλ_M). By (B1) and (B2) of `KRTW-F-moves.md` Theorem B,
Σ_nλ_n^{k+1}β^B_n = (xq)_∞B_{k+1}(x) and Σ_mλ_m^{k+1}β^J_m = (x²Q;Q)_∞J_k(x). ((B1), (B2) are proved there in Z[[q]][[x]];
they hold in 𝔽 because the x^D-coefficients of all series involved are finite sums of rational functions of q regular
at q = 0, and the coefficientwise Taylor map from such series to Q[[q]][[x]] is injective.) ∎ [C8: K_k = Σλ^k t/(…) for k ≤ 6 mod q^30.]

## 8. Verification summary

`code/certify_F_all_k.py` (ALL CHECKS PASSED, `certificates/certify_F_all_k.txt`): C1 commutators and transposes (exact);
C2 (B-v) symbolically; C3/C5a eigen-relations of H′, G, H′♯, G♯ on eigenvectors computed independently from the
triangular systems of T and T♯; C4/C5b product formulas against those eigenvectors (N ≤ 7); C4b/C5c the factorisations
and recursions (N ≤ 8); C6 biorthogonality (all pairs N ≤ 3); C7 Lemma 7.1 (N ≤ 8); C8 (F_k) for k ≤ 6 mod q^30 in all
powers of x; C9 the expansion of u. Negative controls (perturbed t_{1,1}, perturbed V♯_{1,0}) are detected.
Independent series checks of (ii), (iii) to q^300 for k ≤ 8: `code/series_check_q300.py` (separate code).

## 9. Identity (iv); what is not addressed

Identity (iv) is proved in `iv-all-k.md`: in the language above K^Ω_k = Σ_{i=1}^k ω(T^{i−1}OT^{k−i}u) with O := PY𝔲, and
[T, h_{xY}] = −xOT makes the insertions telescope (no off-diagonal matrix elements of O are needed). KRTW Conjecture 2.2
for general z and Conjecture 4.9 are not addressed.
