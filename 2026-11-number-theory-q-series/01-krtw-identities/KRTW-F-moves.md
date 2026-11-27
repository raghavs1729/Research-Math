# The KRTW identities (ii), (iii): proof for k = 3, and the structure of the general case

> **Referee note (2026-11-26, `REFEREE-KRTW-F-moves.md`, fresh adversarial read with independent code in
> `code/krtw_referee/`).** Theorem A, Theorem B and Corollary C: no error found; both move chains replayed
> independently (all 30 intermediate multisums, x⁰…x¹⁶ mod q¹⁰⁰), Bailey steps re-derived by hand, (HL_k) confirmed
> for k ≤ 3 directly from Kostka–Foulkes polynomials. Two errors outside the proofs were found and are corrected
> in this version: the operator of §6 item 2 (now H′, G) and the k = 1 side formula at the end of §5. Open points
> the referee could not settle: novelty beyond the project's searches (no full-text search available), and whether
> (B3)/(B5) are in Bressoud–Ismail–Stanton or Warnaar (they are short consequences of the Bailey lemma and are
> not claimed as new). Corollary C is conditional on the cited theorems of Bartlett–Warnaar and Griffin–Ono–Warnaar.
> **Second, independent referee (parent session, `REFEREE-2-KRTW-F.md` on the parent branch):** Theorem A — proof
> complete, statements match the printed identities; its exposition requests (summation range in Theorem A, the
> step J_3 = Rogers–Ramanujan series in x², two unclear sentences) are incorporated in this version.

Status (2026-11-26, third session). **Proved here:** the x-refined factorisation (F_k) of `KRTW-general.md` §8 for
k = 2 and k = 3 (Theorem A), hence the second and third conjectural Andrews–Gordon-type identities of
Kanade–Russell–Tsuchioka–Warnaar (arXiv:2404.03851v2, after Conjecture 4.9) for **k = 3** (new; k = 1, 2 were
known). **Proved for all k:** closed forms of the two factors and of their product (Theorem B), which identify the
conjectural right-hand side of (F_k) with the C₂ Rogers–Selberg series at the point (√x, −√x); consequently (F_k)
is equivalent to a Hall–Littlewood statement (HL_k) of the same type as KRTW's Conjecture 2.2 (Corollary C,
conditional on two cited theorems), and identity (ii) is equivalent to Conjecture 2.2 at z = q. **Not proved:** (F_k),
(ii), (iii) for k ≥ 4, and identity (iv) for k ≥ 3. The chains of §3 were found by computer search and I do not
know the general pattern; what is known about the general case is collected in §7.

Certificates: `code/krtw_moves_certify.py` → `certificates/krtw_moves_certify.txt` (exact: every move of §3 with its
precondition, the final matching, and an independent series evaluation of all four ends, x⁰…x¹⁶ mod q⁶⁰);
`code/krtw_struct_certify.py` → `certificates/krtw_struct_certify.txt` (series checks of Theorem B and of (F_k), k ≤ 7,
all powers of x, mod q⁹⁰).

Notation. (a;q)_n as usual, (q)_n = (q;q)_n, Q := q², 1/(q)_n := 0 for n < 0. For k ≥ 1 (r_{k+1} = s_0 := 0):

    K_k(x)   := Σ_{r_1≥…≥r_k≥0, 0≤s_1≤…≤s_k} x^{Σ(r_i+2s_i)} ∏_i q^{(r_i+s_i)²+s_i²} / ( (q)_{r_i−r_{i+1}} (Q;Q)_{s_i−s_{i−1}} ),
    B_m(x)   := Σ_{N_1≥…≥N_m≥0} x^{ΣN_i} q^{ΣN_i²} / ( (q)_{N_1−N_2}⋯(q)_{N_{m−1}−N_m} (Q;Q)_{N_m} ),
    C_k(x)   := Σ_{0≤s_1≤…≤s_k} x^{Σs_i} q^{Σs_i²} / ∏_i (Q;Q)_{s_i−s_{i−1}},        J_k(x) := C_k(x)/C_1(x),    C_1(x) = (−xq;Q)_∞.

All series are elements of 𝓡 := Z[[q]][[x]]: in each of these sums the exponent of x of a term is at least the sum of
the differences r_i − r_{i+1}, s_i − s_{i−1} (resp. N_i − N_{i+1}, N_m) that index it, so only finitely many terms contribute to a
given power of x (see §1). KRTW's
identities (iii) and (ii) are the evaluations of K_k(1) and K_k(q) (`KRTW-general.md` §0).

**Conjecture (F_k)** (`KRTW-general.md` §8): K_k(x) = B_{k+1}(x) J_k(x), equivalently K_k C_1 = B_{k+1} C_k.

**Theorem A.** (F_2) and (F_3) hold. Consequently, for k = 3, with both sums running over all integers
r_1 ≥ r_2 ≥ r_3 ≥ 0 and 0 ≤ s_1 ≤ s_2 ≤ s_3 (and r_4 = s_0 := 0),

    Σ ∏_{i=1}^3 q^{(r_i+s_i)²+s_i²+r_i+2s_i}/((q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}}) = (q,q²,q⁸,q⁹,q¹⁰,q¹⁰;q¹⁰)_∞ / ((q²;q²)_∞ (q)_∞),        (ii, k=3)
    Σ ∏_{i=1}^3 q^{(r_i+s_i)²+s_i²}/((q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}})          = (q⁴,q⁵,q⁵,q⁶,q¹⁰,q¹⁰;q¹⁰)_∞ / ((q²;q²)_∞ (q)_∞).        (iii, k=3)

**Theorem B (all k ≥ 1).** In 𝓡:

    (B1)  (xq;q)_∞ B_{k+1}(x)   = Σ_{n≥0} (−1)^n x^{(k+1)n} q^{(k+2)n²} (x²;Q)_n (1−xq^{2n}) / ( (Q;Q)_n (1−x) ),
    (B2)  (x²Q;Q)_∞ J_k(x)      = Σ_{m≥0} (−1)^m x^{(k+1)m} q^{(k+2)m²−m} (x²;Q)_m (1−x²q^{4m}) / ( (Q;Q)_m (1−x²) ),
    (B3)  B_{k+1}(x) = (−xQ;Q)_∞ · Σ_{t_1≥…≥t_k≥0} x^{Σt_i} q^{Σt_i²} (−xq;Q)_{t_k} / ( ∏_{i<k}(Q;Q)_{t_i−t_{i+1}} · (Q;Q)_{t_k} (−xQ;Q)_{t_k} ),
    (B4)  (xq;q)_∞ (x²Q;Q)_∞ B_{k+1}(x) J_k(x) = L_k(x), where

          L_k(x) := ½ Σ_{n,m≥0} (−1)^{n+m} x^{(k+1)(n+m)} q^{(k+2)(n²+m²)−n−m}
                        × (x²;Q)_n (x²;Q)_m (q^n+q^m)(1+xq^{n+m})(1−xq^{2n})(1−xq^{2m}) / ( (Q;Q)_n (Q;Q)_m (1−x)(1−x²) ),

    (B5)  for k = 2m even:  B_{2m+1}(x) = (−xq;q)_∞ 𝔍^{(2)}_m(x)  and  J_{2m}(x) = 𝔍^{(1)}_m(x),  where
          𝔍^{(e)}_m(x) := Σ_{σ_1≥…≥σ_m≥0} x^{2Σσ_i} q^{2Σσ_i²} / ( ∏_{i<m}(Q;Q)_{σ_i−σ_{i+1}} · (Q;Q)_{σ_m} (−xq^e;Q)_{σ_m} ).

(In (B1), (B2), (B4) the quotients (x²;Q)_n/(1−x), (x²;Q)_m/(1−x²) are polynomials for n, m ≥ 1; for n = 0 the factor
(1−xq^{2n})/(1−x) is 1, and similarly in the other cases; so every term lies in 𝓡.)

L_k(x) is exactly the series L^{(0)}_k(y_1,y_2;q) of the level-k C_n Rogers–Selberg identity (n = 2) of
Bartlett–Warnaar, as displayed in Griffin–Ono–Warnaar, arXiv:1401.7718, §3.1, at (y_1,y_2) = (√x, −√x) (§5).

**Corollary C** (conditional on the two cited theorems; §5). (a) Given the C₂ Rogers–Selberg identity,
(F_k) is equivalent to

    (HL_k)   K_k(x) = Σ_{λ: λ_1≤k} q^{|λ|} P'_{2λ}(√x, −√x; q) = Σ_{λ: λ_1≤k} q^{|λ|} P_{2λ}(X_x; q),   X_x := {±√x·q^i : i ≥ 0},

(P_λ the Hall–Littlewood symmetric function with parameter q, P' the modified one), and by Theorem A (HL_k)
holds for k ≤ 3. (b) Given Griffin–Ono–Warnaar's Theorem 1.3 (n = 2), identity (ii) for a given k is equivalent to
KRTW's Conjecture 2.2 at z = q for that k; hence, by Theorem A, **Conjecture 2.2 holds at z = q for k = 3** (and at
z = 1 for all k, `KRTW-general.md`).

## 1. Generalized Nahm data and two moves

A *datum* is D = (M, B, w, d): M a symmetric r×r matrix over Q, B ∈ Q^r, w ∈ Z^r_{≥1}, d ∈ {1,2,4,…}^r, such that
½nᵀMn + B·n ∈ Z for all n ∈ N^r. Put

    F(D) := Σ_{n∈N^r} x^{w·n} q^{½nᵀMn+B·n} / ∏_i (q^{d_i};q^{d_i})_{n_i}  ∈ Z((q))[[x]].

(For each power x^D only the finitely many n with w·n = D contribute, so the coefficient of x^D is a finite sum
of Laurent series in q; no positivity of M is needed. In all data below the exponents are ≥ 0, so F(D) ∈ 𝓡.)

**Lemma M1 (SPLIT).** Let D′ be obtained from D by replacing the species a (base d = d_a) by two species a, t of
base 2d, both with the row of a (M′_{aa} = M′_{tt} = M′_{at} = M_{aa}, M′_{ak} = M′_{tk} = M_{ak}), with B′_a = B_a, B′_t = B_a + d,
w′_a = w′_t = w_a. Then F(D′) = F(D).
*Proof.* 1/(z;q^d)_∞ = 1/((z;q^{2d})_∞ (zq^d;q^{2d})_∞); comparing coefficients of z^n (Euler),
1/(q^d;q^d)_n = Σ_{u+v=n} q^{dv}/((q^{2d};q^{2d})_u (q^{2d};q^{2d})_v). Substitute this for the factor 1/(q^d;q^d)_{n_a} and note that
the summand of F(D) depends on (u,v) only through n_a = u+v, apart from the new factor q^{dv}. ∎

**Lemma M2 (FUSE / DEFUSE).** Let a ≠ b be species with d_a = d_b = d, put Q_d := q^d, and let D′ have one more
species c with d_c = d, w_c = w_a + w_b, B_c = B_a + B_b, M′_{cc} = M_{aa}+M_{bb}+2M_{ab}, M′_{ac} = M_{aa}+M_{ab}, M′_{bc} = M_{bb}+M_{ab},
M′_{ck} = M_{ak}+M_{bk} (k ≠ a,b,c), M′_{ab} = M_{ab} + d, all other entries as in D. Then F(D′) = F(D).
*Proof.* The Durfee rectangle identity 1/((Q_d)_A (Q_d)_B) = Σ_{c≥0} Q_d^{(A−c)(B−c)}/((Q_d)_c (Q_d)_{A−c} (Q_d)_{B−c}) (Andrews,
*The Theory of Partitions*, (3.3.10) in the limit, or Cauchy; it is `proofs/E8-E6.md`, Lemma 1). In F(D) write
n_a = A, n_b = B, replace 1/((Q_d)_A(Q_d)_B) by the right side of the identity, and substitute A = n′_a + c, B = n′_b + c (so that
A − c = n′_a, B − c = n′_b). Substituting n_a = n′_a + c, n_b = n′_b + c in ½nᵀMn + B·n gives a quadratic form in (n′_a, n′_b, c, others) whose
coefficients are exactly the entries of M′ and B′ listed above, except that the coefficient of n′_a n′_b is still M_{ab};
the factor Q_d^{(A−c)(B−c)} = q^{d n′_a n′_b} of the identity adds d to it, which gives M′_{ab} = M_{ab} + d. The power of x becomes
w_a(n′_a+c)+w_b(n′_b+c)+… = w_a n′_a + w_b n′_b + (w_a+w_b)c + … . The regrouping
is legitimate because for each power of x both sides are finite sums. ∎
DEFUSE(a,b,c) is the same identity read from right to left: it applies to a datum D′ whenever c has exactly
the row, B and w prescribed above in terms of a and b (with M_{ab} := M′_{ab} − d).

These are the only tools. The moves are implemented in `code/krtw_moves_certify.py` (functions `split`, `fuse`,
`defuse`; the precondition of every DEFUSE is checked in exact rational arithmetic and the program stops if it
fails), independently of the search code that found the chains.

## 2. The data of the four sums

* K_k: species a_j := r_j − r_{j+1} (base q, j = 1..k) and b_l := s_l − s_{l−1} (base Q, l = 1..k); r_i = Σ_{j≥i} a_j, s_i = Σ_{l≤i} b_l;
  exponent Σ_i (r_i+s_i)² + s_i², weight Σ_i (r_i + 2s_i).
* B_m: species n_j := N_j − N_{j+1} (base q, j < m), n_m := N_m (base Q); exponent ΣN_i², weight ΣN_i.
* C_k: species e_l := s_l − s_{l−1} (base Q); exponent Σs_i², weight Σs_i.
* AG_j(x²;Q) := Σ_{σ_1≥…≥σ_j≥0} x^{2Σσ} q^{2Σσ²}/∏(Q;Q)_{σ_i−σ_{i+1}}: species g_l (base Q).

By (B2) and the Rogers–Selberg/Andrews–Gordon series in base Q (the Bailey chain with j+1 = (k+1)/2 full steps
applied to the unit pair relative to x², see §4) J_k(x) = AG_{(k−1)/2}(x²;Q) for odd k; in particular **J_3(x) =
AG_1(x²;Q) = Σ_σ x^{2σ} q^{2σ²}/(Q;Q)_σ** and J_1 = 1. Written out for k = 3: the Rogers–Ramanujan series with a parameter
satisfies the classical Rogers–Selberg (Watson) identity

    Σ_{σ≥0} a^σ p^{σ²}/(p;p)_σ = (ap;p)_∞^{−1} Σ_{m≥0} (−1)^m a^{2m} p^{m(5m−1)/2} (a;p)_m (1−ap^{2m}) / ( (p;p)_m (1−a) )

(two steps (L∞) of §4 applied to the unit Bailey pair relative to a, in base p). For p = Q = q², a = x² the right
side is (x²Q;Q)_∞^{−1}Σ_m(−1)^m x^{4m}q^{5m²−m}(x²;Q)_m(1−x²q^{4m})/((Q;Q)_m(1−x²)), which is J_3(x) by (B2) with k = 3 (there (k+1)m = 4m and
(k+2)m²−m = 5m²−m). Hence J_3(x) = Σ_σ x^{2σ}q^{2σ²}/(Q;Q)_σ: the Rogers–Ramanujan series in base q² with parameter x².

## 3. Proof of Theorem A

SPLIT every base-q species (a_j → a_j′, a_j″ and n_j → n_j′, n_j″). Convention (Lemma M1 with a = a_j′, t = a_j″):
the primed copy keeps the linear term B of the original species, the doubly primed copy has B + 1 (the factor q^{v}
of the Euler splitting); the chains below are only valid with this reading. After that all species have base Q.
Notation in the tables: (u|c) is the species u after a DEFUSE in which c was removed (so the new summation
index is n_u + n_c); [u+v] is the species created by FUSE(u,v).

**k = 2.** K-side (start: a_1′, a_2′, b_1, b_2, a_1″, a_2″):

    FUSE   a_1′ , b_2
    DEFUSE a_2″        = fusion of  a_1′ , a_1″
    DEFUSE [a_1′+b_2]  = fusion of  (a_1′|a_2″) , a_2′                       → rank 5.

B-side (start: n_1′, n_2′, n_3, n_1″, n_2″):

    DEFUSE n_2″ = fusion of n_1′ , n_1″ ;      DEFUSE n_3 = fusion of (n_1′|n_2″) , n_2′ .

After these two moves the species ((n_1′|n_2″)|n_3) is *free* (zero coupling to the others) with exponent n²,
weight 1, base Q: it contributes the factor Σ_n x^n q^{n²}/(Q;Q)_n = C_1(x). Removing it leaves a datum for
B_3/C_1; its product with the datum of C_2 is a datum for B_3 C_2/C_1 = B_3 J_2, and one more move

    FUSE (n_2′|n_3) , e_1

gives the K-side result up to the permutation of species recorded in the certificate. Hence K_2 = B_3 J_2. In
closed form the common value is (p, u, v, α, γ ≥ 0; all Pochhammer symbols in base Q)

    K_2(x) = Σ x^{4p+2u+2v+α+γ} q^{p²+(p+u+v)²+(p+u+α)²+(p+v+γ)²+γ} / ( (Q)_p (Q)_u (Q)_v (Q)_α (Q)_γ ).                       (3.1)

(Direct reading of (3.1): sum α and γ by Euler, then p by Durfee with U = p+u, V = p+v:
K_2(x) = (−xq;q)_∞ Σ_{U,V} x^{2U+2V} q^{2U²+2V²}/((Q)_U (Q)_V (−xq;Q)_U (−xQ;Q)_V) = (−xq;q)_∞ 𝔍^{(1)}_1 𝔍^{(2)}_1, in accordance with (B5).)

**k = 3.** K-side (start: a_1′, a_2′, a_3′, b_1, b_2, b_3, a_1″, a_2″, a_3″):

    FUSE   a_1′ , b_2
    DEFUSE a_2″                  = fusion of  a_1′ , a_1″
    FUSE   a_2′ , b_3
    DEFUSE a_3′                  = fusion of  (a_1′|a_2″) , a_2′
    FUSE   (a_2′|a_3′) , b_2
    DEFUSE [(a_2′|a_3′)+b_2]     = fusion of  b_3 , [a_2′+b_3]
    FUSE   b_2 , (a_1″|a_2″)
    DEFUSE [a_1′+b_2]            = fusion of  ((a_1′|a_2″)|a_3′) , b_2
    DEFUSE [b_2+(a_1″|a_2″)]     = fusion of  (b_3|[(a_2′|a_3′)+b_2]) , a_3″              → rank 8.

B-side (start: n_1′, n_2′, n_3′, n_4, n_1″, n_2″, n_3″; then the species g_1 of J_3 = AG_1(x²;Q)):

    DEFUSE n_2″ = fusion of n_1′ , n_1″ ;    DEFUSE n_3′ = fusion of (n_1′|n_2″) , n_2′ ;    FUSE n_4 , g_1 ;    FUSE (n_2′|n_3′) , g_1    → rank 8.

The two final data coincide after the permutation of species listed in `certificates/krtw_moves_certify.txt`.
Each move is an identity by Lemmas M1, M2 (the preconditions of the DEFUSE moves are finitely many equalities
between rational numbers, verified by the certificate), so F(K_3-datum) = F(B_4-datum ⊗ AG_1-datum), i.e.
K_3 = B_4 J_3. ∎ (F_2), (F_3).

*From (F_k) to (ii), (iii)* (as in `KRTW-general.md` §8). All series involved (K_k, B_{k+1}, C_k, J_k and both sides
of (B1), (B2)) lie in the subring Z[x][[q]] of 𝓡 (every coefficient of q^N is a polynomial in x, because the
q-exponents grow with the summation indices), on which the substitutions x = 1 and x = q are defined and are
ring homomorphisms to Z[[q]]. At x = 1 and x = q the sums in (B1), (B2) are
bilateral theta series: with θ(z;p) := (z, p/z, p; p)_∞,

    (q)_∞ B_{k+1}(1) = Σ_{n∈Z} (−1)^n q^{(k+2)n²} = θ(q^{k+2};q^{2k+4}),           (q)_∞ B_{k+1}(q) = θ(q;q^{2k+4}),
    (Q;Q)_∞ J_k(1) = Σ_{m∈Z} (−1)^m q^{(k+2)m²−m} = θ(q^{k+1};q^{2k+4}),       (Q;Q)_∞ J_k(q) = θ(q²;q^{2k+4})

((x²;Q)_n(1−xq^{2n})/((Q;Q)_n(1−x)) → 2 for n ≥ 1, → 1 for n = 0 as x → 1; at x = q, (q²;Q)_n(1−q^{2n+1})/((Q;Q)_n(1−q)) = (1−q^{2n+1})/(1−q)
and the terms n and −n−1 of the bilateral sum combine; (xq)_∞|_{x=q} = (q)_∞/(1−q). For J: at x = 1 the factor
(x²;Q)_m(1−x²q^{4m})/((Q;Q)_m(1−x²)) is 1 + q^{2m} for m ≥ 1 (and 1 for m = 0), which turns the unilateral sum into Σ_{m∈Z}(−1)^m q^{(k+2)m²−m};
at x = q it is (1−q^{4m+2})/(1−q²), and m, −m−1 combine. These are
Bressoud's even-modulus identities, resp. the Jacobi triple product.) Hence K_k(1) = θ(q^{k+2};q^{2k+4})θ(q^{k+1};q^{2k+4})/((q)_∞(Q;Q)_∞)
and K_k(q) = θ(q;q^{2k+4})θ(q²;q^{2k+4})/((q)_∞(Q;Q)_∞), which are the right sides of (iii) and (ii). For k = 3 this is
Theorem A. (Series check of the two k = 3 identities: `certificates/krtw_four_test.txt`, mod q⁴⁰⁰.) ∎

Remark (what is new). For k = 2, (ii) and (iii) follow from KRTW Theorem 4.10 and Bartholdi's theorem
(arXiv:2611.04020); the proof above is independent of both and elementary. For k = 3 no proof was found in the
literature (searches of 26 Nov 2026: `PRIOR-ART.md`, §"KRTW second round").

## 4. Proof of Theorem B

Bailey pairs (α_n, β_n) relative to a in base p: β_n = Σ_{r≤n} α_r/((p;p)_{n−r}(ap;p)_{n+r}). Facts used (Bailey's lemma in the
two limiting cases; Bressoud–Ismail–Stanton, *Change of base in Bailey pairs*, §1, §3; rational-function
identities in p and √a):

 (L∞) α′_n = a^n p^{n²} α_n,  β′_n = Σ_{j≤n} a^j p^{j²} β_j/(p;p)_{n−j};
 (S2) α′_n = a^{n/2} p^{n²/2} α_n,  β′_n = Σ_{j≤n} (−√(ap);p)_j a^{j/2} p^{j²/2} β_j/((p;p)_{n−j}(−√(ap);p)_n);
 (S2′) α′_n = a^{n/2} p^{n(n+1)/2} ((1+√a)/(1+√a p^n)) α_n,  β′_n = Σ_{j≤n} (−√a;p)_j a^{j/2} p^{j(j+1)/2} β_j/((p;p)_{n−j}(−√a p;p)_n)
       (Bailey's lemma with ρ_1 → ∞, ρ_2 = −√a);
 (U)  the unit pair β_n = δ_{n0}, α_n = (−1)^n p^{n(n−1)/2}(a;p)_n(1−ap^{2n})/((p;p)_n(1−a));
 (lim) if (α,β) is a Bailey pair with coefficients in 𝓡 (for the specialisations below) then
       lim_{n→∞} β_n = (1/((p;p)_∞(ap;p)_∞)) Σ_r α_r (as in `KRTW-general.md`, Prop. 6).

*(B1).* Base p = q, a = x. Claim: β_n = 1/(Q;Q)_n, α_n = (−1)^n q^{n²}(x²;Q)_n(1−xq^{2n})/((Q;Q)_n(1−x)) is a Bailey pair. Indeed, by
the inversion of Bailey pairs α_n = ((1−aq^{2n})/(1−a)) Σ_{j≤n} (a;q)_{n+j}(−1)^{n−j}q^{(n−j)(n−j−1)/2}β_j/(q)_{n−j}; with β_j = 1/((q)_j(−q)_j) the sum
is ((a)_n(−1)^n q^{n(n−1)/2}/(q)_n)·₂φ₁(q^{−n}, aq^n; −q; q, q) and by the q-Chu–Vandermonde sum
₂φ₁(q^{−n}, b; c; q, q) = (c/b)_n b^n/(c)_n this is ((a)_n(−1)^n q^{n(n−1)/2}/(q)_n)·(−q^{1−n}/a)_n (aq^n)^n/(−q)_n; since
(−q^{1−n}/a)_n = (−a)_n a^{−n}q^{−n(n−1)/2}, α_n = (−1)^n q^{n²}(a)_n(−a)_n(1−aq^{2n})/((Q;Q)_n(1−a)), as claimed. Apply (L∞) k+1 times and (lim):
the β-side is B_{k+1}(x)/(q)_∞ (indices N_1 ≥ … ≥ N_{k+1}), the α-side is (1/((q)_∞(xq)_∞)) Σ x^{(k+1)n}q^{(k+1)n²}α_n. ∎
*(B2).* This is the statement J_k = C_k/C_1 of `KRTW-general.md` §8 (proved there: (S2)^{k+1}(U) in base Q, a = x²).
*(B3).* Base Q, a = x²: (S2′)(U) gives β_n = 1/((Q;Q)_n(−xQ;Q)_n), α_n = x^n q^{n²+n}((1+x)/(1+xq^{2n}))·(−1)^n q^{n²−n}(x²;Q)_n(1−x²q^{4n})/((Q;Q)_n(1−x²))
= (−1)^n x^n q^{2n²}(x²;Q)_n(1−xq^{2n})/((Q;Q)_n(1−x)). Apply (S2) k times (√(aQ) = xq, a^{n/2}Q^{n²/2} = x^n q^{n²}): α^{(k)}_n = x^{kn}q^{kn²}α_n is, by (B1),
the n-th term of (xq;q)_∞ B_{k+1}(x), and β^{(k)}_N = (1/(−xq;Q)_N) Σ_{N≥t_1≥…≥t_k} x^{Σt}q^{Σt²}(−xq;Q)_{t_k}/((Q;Q)_{N−t_1}⋯(Q;Q)_{t_k}(−xQ;Q)_{t_k}). By (lim),
(xq;q)_∞B_{k+1}(x) = (x²Q;Q)_∞ (−xq;Q)_∞^{−1} Σ_t (…), and (x²Q;Q)_∞/((xq;q)_∞(−xq;Q)_∞) = (−xQ;Q)_∞ because
(x²Q;Q)_∞ = (xq;q)_∞(−xq;q)_∞ and (−xq;q)_∞ = (−xq;Q)_∞(−xQ;Q)_∞. ∎
*(B5).* J_{2m} = 𝔍^{(1)}_m: (S2)(U) then (L∞)^m in base Q, a = x² has α_n = x^{(2m+1)n}q^{(2m+1)n²}α^U_n, the n-th term of (B2), and
β-limit 𝔍^{(1)}_m/(Q;Q)_∞. B_{2m+1} = (−xq;q)_∞𝔍^{(2)}_m: (S2′)(U) then (L∞)^m has α_n = x^{2mn}q^{2mn²}·(the α of (B3)) = the n-th term of (B1)
with k = 2m, and β-limit 𝔍^{(2)}_m/(Q;Q)_∞; (x²Q;Q)_∞/(xq;q)_∞ = (−xq;q)_∞. ∎
*(B4).* Write b_n, j_m for the n-th and m-th terms of (B1), (B2). Then (xq)_∞(x²Q;Q)_∞B_{k+1}J_k = Σ_{n,m} b_n j_m = ½ Σ_{n,m} (b_n j_m + b_m j_n),
and b_n j_m + b_m j_n = (−1)^{n+m} x^{(k+1)(n+m)} q^{(k+2)(n²+m²)} ((x²;Q)_n(x²;Q)_m/((Q;Q)_n(Q;Q)_m(1−x)(1−x²))) · E, with

    E = q^{−m}(1−xq^{2n})(1−x²q^{4m}) + q^{−n}(1−xq^{2m})(1−x²q^{4n}).

Put u = q^n, v = q^m. Then uv·E = u(1−xu²)(1−x²v⁴) + v(1−xv²)(1−x²u⁴) = (u+v)(1+xuv)(1−xu²)(1−xv²): expanding the right
side, (u+v)(1+xuv)(1 − x(u²+v²) + x²u²v²) = (u+v) − x(u³+v³) − x²uv(u³+v³) + x³u³v³(u+v), which is the expansion of the left
side term by term. Hence E = q^{−n−m}(q^n+q^m)(1+xq^{n+m})(1−xq^{2n})(1−xq^{2m}), which gives L_k(x). ∎

## 5. The Hall–Littlewood form (Corollary C)

*The series L_k.* Griffin–Ono–Warnaar (arXiv:1401.7718, §3.1, the display defining L^{(0)}_m after their equation
labelled Eq_Cn-RS-m) put, for y = (y_1,…,y_n),

    L^{(0)}_m(y;q) := Σ_{r∈N^n} (Δ_C(yq^r)/Δ_C(y)) ∏_i y_i^{2(m+1)r_i} q^{(m+1)r_i² + n·r_i(r_i−1)/2} ∏_{i,j} (−y_i/y_j)^{r_i} (y_iy_j)_{r_i}/(qy_i/y_j)_{r_i},
    Δ_C(y) := ∏_i (1−y_i²) ∏_{i<j} (y_i−y_j)(y_iy_j−1).

For n = 2, y = (α, −α), α² = x, (r_1,r_2) = (n,m): ∏y_i^{2(k+1)r_i}q^{(k+1)r_i²+r_i²−r_i} = x^{(k+1)(n+m)}q^{(k+2)(n²+m²)−n−m}; ∏_{i,j}(−y_i/y_j)^{r_i} = (−1)^{n+m};
∏_j (y_1y_j)_n/(qy_1/y_j)_n = (x)_n(−x)_n/((q)_n(−q)_n) = (x²;Q)_n/(Q;Q)_n, likewise for m; and
Δ_C(yq^r)/Δ_C(y) = (1−xq^{2n})(1−xq^{2m})(q^n+q^m)(1+xq^{n+m})/(2(1−x)²(1+x)). So **L^{(0)}_k(√x,−√x;q) = L_k(x)** (the function of (B4)).

*The cited theorems.* (RS) [Bartlett–Warnaar, *Hall–Littlewood polynomials and characters of affine Lie
algebras*, Adv. Math. 285 (2015), arXiv:1304.1602 — the nonterminating limit of their C_n Andrews transformation
(Theorem 4.2; the right-hand side is evaluated in their section "The right-hand side of the C_n Andrews
transformation") combined with their Littlewood-type identity for the functions R_k (the corollary of their
Theorem 3.1), in the limit b, c → ∞, where the prefactor tends to ∏(qy_i²)_∞∏_{i<j}(qy_iy_j)_∞; restated as "the level-m C_n Rogers–Selberg identity" in GOW §3.1]:

    Σ_{λ: λ_1≤m} q^{|λ|} P′_{2λ}(y;q) = L^{(0)}_m(y;q) / ( ∏_i (qy_i²;q)_∞ ∏_{i<j} (qy_iy_j;q)_∞ ).

Caveat on normalisation: the display in the arXiv version of GOW omits the infinite products on the right. With
them the identity is the classical Rogers–Selberg identity for n = 1 (GOW's own (3.1)), tends to the Littlewood
sum Σ_λ q^{|λ|}P′_{2λ} = ∏(qy_i²)_∞^{−1}∏(qy_iy_j)_∞^{−1} for m → ∞, and for n = 2, y = (√x,−√x) agrees with the explicit multisum for the left
side (below) for k = 1, 2, 3 (`code/krtw_struct_certify.py` checks (B4); the referee's independent scripts, archived in
`code/krtw_referee/`, compare the multisum with K_k (`bwmulti.py`, x⁰…x¹⁰ mod q³⁰), compute Σ_{λ_1≤k}q^{|λ|}P′_{2λ}(√x,−√x;q)
directly from Kostka–Foulkes polynomials (`hl_direct.py`, k ≤ 3, x⁰…x⁷ mod q³⁰), and test (RS) at a generic
point (y_1,y_2) for m ≤ 3 (`rs_generic.py`): true with the products, false without them). For
y = (√x, −√x) the product is (xq;q)_∞²(−xq;q)_∞ = (xq;q)_∞(x²Q;Q)_∞, the prefactor in (B4).
(BW 3.1) P′_μ(y;q) = Σ ∏_{l≥1} f^{(1)}_{r^{(l)},r^{(l+1)}}(y;q), summed over r^{(1)}, r^{(2)}, … ∈ N^n with |r^{(l)}| = μ′_l, where
f^{(τ)}_{r,s}(y;q) = ∏_i (y_i^{r_i}q^{r_i(r_i−1)/2})^τ ∏_{i,j} (qy_i/y_j)_{r_i−r_j}/(qy_i/y_j)_{r_i−s_j}. For μ = 2λ and y = (√x,−√x) this gives the 2k-fold sum

    Σ_{λ_1≤k} q^{|λ|}P′_{2λ}(√x,−√x;q) = Σ_{(a_l,b_l)_{l=1..k}} ∏_l x^{a_l+b_l} q^{a_l²+b_l²} ω(a_l−b_l) / ( (q)_{a_l−a_{l+1}} (q)_{b_l−b_{l+1}} (−q)_{a_l−b_{l+1}} (−q)_{b_l−a_{l+1}} ),

ω(0) = 1, ω(d) = ½q^{|d|(|d|−1)/2}(1+q^{|d|}), (a_{k+1},b_{k+1}) = (0,0), with 1/(−q)_{−j} := (−q^{1−j};q)_j.

*Proof of Corollary C.* (a) By (B4) and (RS), B_{k+1}(x)J_k(x) = Σ_{λ_1≤k} q^{|λ|}P′_{2λ}(√x,−√x;q); and
P′_μ(y_1,…,y_n;q) = P_μ(Y;q) for the alphabet Y = {y_i q^j : j ≥ 0} (definition of P′: the plethystic substitution p_r ↦ p_r/(1−q^r)).
(b) GOW Theorem 1.3 (D^{(2)}_{n+1}, n ≥ 2; quoted verbatim in `PRIOR-ART.md`) for n = 2, m = k reads
Σ_{λ_1≤k} q^{2|λ|}P_{2λ}(1,q,q²,…;q²) = (q,q²,q^{2k+2},q^{2k+3},q^{2k+4},q^{2k+4};q^{2k+4})_∞/((q²;q²)_∞(q)_∞), the right side of (ii); and the right side
of KRTW's Conjecture 2.2 at z = q is the left side of (ii) (z^{r_i+s_i}q^{…+s_i} at z = q gives the exponent (r_i+s_i)²+s_i²+r_i+2s_i). ∎

So the same 2k-fold sum K_k carries two different Hall–Littlewood interpretations: along KRTW's grading
z^{Σ(r_i+s_i)} (base q², alphabet {q^i}; conjectural) and along the grading x^{Σ(r_i+2s_i)} used here (base q, alphabet
{±√x q^i}; equivalent to (F_k)); they meet at x = q, z = q. For k = 1 the display above reads K_1(x) = Σ_{a,b≥0} x^{a+b}q^{a²+b²}ω(a−b)/((Q)_a(Q)_b) (ordered pairs), i.e.
P′_{(2^N)}(√x,−√x;q) = x^N Σ_{a+b=N} q^{a²+b²−N}ω(a−b)/((Q)_a(Q)_b) (checked for N ≤ 7 by the referee, `code/krtw_referee/hl_direct.py`;
an earlier version of this paragraph had a spurious factor 2^{[a≠b]}).

## 6. Further exact facts about the general case (no proof of (F_k) for k ≥ 4)

1. **The x-grading is the only grading that survives.** The chains of §3 preserve the weight w (the DEFUSE
   preconditions include w_c = w_a + w_b), but no second linear grading: in particular neither KRTW's z^{Σ(r_i+s_i)}
   nor a marking of r_1 is preserved (already for k = 2 the precondition for a_2″ fails for the grading y^{r_1}). This
   is why the finite R-refinement that proves identity (i) has no analogue here; the two natural
   two-parameter versions of (F_2) (weights a^{Σ(r_i+s_i)}b^{Σs_i} against a^{N_1+N_2}b^{N_3}·J_2(a or b)) are false (checked
   mod q²²).
2. **Commuting difference operators.** Let T be the transfer operator of the K-chain on functions of
   (R,S): (TF)(R,S) = Σ_{r≤R, s≥S} x^{r+2s}q^{(r+s)²+s²}F(r,s)/((q)_{R−r}(Q;Q)_{s−S}), and P, Y the multiplications by q^R, q^{2S}, 𝖲F(R,S) = F(R−1,S),
   𝖴F(R,S) = F(R,S+1). Then T commutes with

       H′ = P^{−1}Y^{−1}(1 − 𝖲) + x P^{−1}(1 − 𝖲) + x P + x² PY(1 − 𝖴),
       G  = P^{−2}Y^{−1}(1 − 𝖲)(1 − q𝖲) + Y^{−1} − x(1 + q^{−1})𝖲 + x² Y(1 − 𝖴) + x² P²Y(1 − 𝖴)

   (exact arithmetic at two rational points (q,x), on point masses; `code/krtw_referee/hamtest3.py`, `commutant.py`;
   numerical verification, not a proof). In the span of the monomials P^aY^b𝖲^c𝖴^d, |a|,|b| ≤ 2, 0 ≤ c,d ≤ 2, the commutant of T is
   4-dimensional: constants, H′, H′², G. [Correction: an earlier version of this item displayed an operator H with
   extra powers of q, which does not commute with T as defined here, and claimed uniqueness; both were wrong
   and were corrected by the referee, `REFEREE-KRTW-F-moves.md` F1.] H′ has the four "weights" q^{±R}, q^{±(R+2S)}: it is
   of relativistic (q-)Toda type of rank two with unequal root lengths, and two independent commuting operators
   are what a rank-two integrable system should have. This is consistent with (B4): the spectrum of T read off
   from L_k is {x^{n+m}q^{n²+m²}}, n, m ≥ 0, symmetric in (n,m), i.e. indexed by orbits of the Weyl group of C₂.
3. **Finite refinement.** With A_j(R;x) := Σ_{R≥r_1} ((q)_R/(q)_{R−r_1})·(summand of K_j) one has A_j(R;x) = Σ_{i≥0} (−xq^{R+1})^i q^{i(i−1)/2} 𝒦_j(i;x)/(q)_i,
   𝒦_j(0) = K_j, and for j = 2 (numerically, mod q³⁰): 𝒦_2(0;x) = (−xq;q)_∞ 𝔍^{(1)}_1(x)𝔍^{(2)}_1(x), 𝒦_2(1;x) = (−xq;q)_∞ 𝔍^{(1)}_1(x)𝔍^{(1)}_1(xq),
   but 𝒦_2(2;x) is not of this product form (its quotient by (−xq;q)_∞𝔍^{(1)}_1(x) has negative coefficients).
4. **Spectral coefficients.** Writing A_j(R;x) = (q)_R Σ_d α^{(j)}_d(x)/((q)_{R−d}(xq)_{R+d}) (Bailey pair relative to x), the α's
   are α^{(j)}_d = (x²Q;Q)_∞^{−1} Σ_{n≥m≥0, n−m≤d≤n+m} c^{n,m}_d(x)(x^{n+m}q^{n²+m²})^j, where the c^{n,m}_d are determined recursively by the level
   recursion; for m ≤ 1 (computed symbolically, n ≤ 3) they are products of cyclotomic-type factors, e.g.
   c^{1,1}_0 = −x²q²(1−x²q⁴)/(1−q²), c^{1,1}_1 = x²q²(1+xq)(1−xq²)²(1+xq²)/((1−q)(1−xq³)), c^{1,1}_2 = x²q⁵(1−xq)²(1+xq)(1−xq²)(1+xq²)/((1−q²)²(1−xq³)),
   with c^{1,1}_0 + c^{1,1}_1 + c^{1,1}_2 = the (1,1)-term of L_1-type (B4). A closed formula for all c^{n,m}_d, together with the
   normalisation Σ_{n,m} c^{n,m}_d = (x²Q;Q)_∞·(−1)^d x^d q^{d²+d(d−1)/2}(x)_d(1−xq^{2d})/((q)_d(1−x)) (which expresses A_0 = 1), would prove (F_k)
   for all k; I have not found it.
5. **Identity (iv)** (the Ω-weighted sum) is equivalent, given Bressoud's identities, to
   K^Ω_k(x) = B^{(k)}_{k+1}(x)J_k(x) at x = 1 (`KRTW-general.md` §8). Ω is a sum of k terms with different linear exponents, none
   of which is a function of the surviving grading, so the chains of §3 do not apply termwise. Open for k ≥ 3.

## 7. How the chains were found, and openness

Breadth-first search in the class of a datum under SPLIT/FUSE/DEFUSE with bounded rank
(`code/classsearch.py`, extended in this session to the x-weight; for k = 3: all species split first, then
only FUSE/DEFUSE in base Q, rank ≤ 11 on the K-side (7·10⁵ data) and ≤ 9 on the B-side, meeting in the middle).
The search is not part of the proof. For k = 4 the same search is out of reach (rank 12 after splitting).

Openness (search of 26 Nov 2026, details in `PRIOR-ART.md`): no proof of (ii)–(iv) for k ≥ 3, of KRTW
Conjecture 2.2 for general z, of Conjecture 4.9, or of Primc's D^{(2)}_3 conjecture at level k ≥ 3 was found
(arXiv listings to 24 Nov 2026, the papers citing KRTW, the trackers named there). The usual caveat applies:
no full-text search.
