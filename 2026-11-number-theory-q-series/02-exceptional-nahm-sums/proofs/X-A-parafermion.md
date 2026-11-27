# Sun–Wang's Conjecture 1.1 for the pairs (X, A_{k−1}): Nahm sums as total parafermion characters

> **Referee note (2026-11-26, `REFEREE-X-A-parafermion.md`; fresh adversarial read with independent code in
> `code/xa_referee/`, Freudenthal's formula instead of the Racah recursion).** No mathematical error found in the
> character identities, the dictionary to the cited theorems (all quotations checked against the TeX sources) or
> the prefactor (Lemma 4 re-derived by hand and checked exactly for all types, k ≤ 12); every sector agrees in all
> cases tested, and deliberately wrong inputs are detected. **Two findings changed this note:** (1) the modularity
> step: for B, C, F_4, G_2 it holds via the S-transformation to the transposed untwisted algebra, given Kac §13.9
> and Kac–Peterson's Theorem A; for A_{2l}^{(2)} (all tadpole pairs, including (T_1,A_r)) the congruence-subgroup property
> of string functions is an **unverified citation** — the primary sources could not be accessed; (2) for X simply
> laced the statement, prefactor included, is exactly Kuniba–Nakanishi–Suzuki's 1993 formula (proved by Georgiev
> and Butorac–Kožić–Primc) plus Kac–Peterson, so Theorem 1 is an observation, not a new theorem; the new content is
> the twisted part (Theorem 2, Lemma 4 for twisted algebras, Theorem 3) and the consequences for the conjectures
> of Sun–Wang and Wang–Wang. Smaller corrections (Mizuno locator, exponent lattice for T_l, definition of c_μ for
> A_{2l}^{(2)}, remark on c versus c_eff) are incorporated.

Written 2026-11-26 (third session). Certificate (series check, independent of the cited papers):
`code/xa_parafermion_certify.py` → `certificates/xa_parafermion_certify.txt`.

**Status in one paragraph.** For every Dynkin type X ∈ {A_l, D_l, E_6, E_7, E_8, B_l, C_l, F_4, G_2, T_l} and every k ≥ 2,
the generalized Nahm sum of the Sun–Wang pair (X, A_{k−1}) is the character of the whole parafermionic space
of the level-k vacuum module of an affine Lie algebra g(X), hence a finite sum of Kac–Peterson string functions
times an eta product, hence modular, with exactly the prefactor q^{−c(X,A_{k−1})/24} of Sun–Wang's Conjecture 1.1
(Theorems 1–3 below). **This is a corollary of published theorems**: Butorac–Kožić–Primc (X simply laced),
Okado–Takenaka (X = B_l, C_l, F_4, G_2), Takenaka (X = T_l; arXiv preprint), Lepowsky–Wilson, Kac–Peterson and
Kac's formula for the modular anomaly. What is done here is the identification with Sun–Wang's pairs, the
computation of the prefactor (Lemma 4: Sun–Wang's central charge is the Kac–Peterson anomaly, uniformly), and an
independent numerical check. **For X simply laced nothing beyond bookkeeping is added**: Kuniba–Nakanishi–Suzuki
(1993) state exactly this sum, with the prefactor q^{−c_PF/24} and "without any restriction on the n-sum", as the
character of the total parafermion space, together with ch(Ω^Λ_λ) = η^r c^Λ_λ; the formula was proved by Georgiev (A_n)
and Butorac–Kožić–Primc (all types), and Nahm (hep-th/0404120) names the parafermion coset as the conformal field
theory of the pairs (Y, A_{k−1}). The twisted identification (B, C, F_4, G_2, T) and the uniform prefactor are the new
bookkeeping. The mechanism was used before in two rank-two cases: Mizuno (arXiv:2305.02267, §4, comments on his
Table 1) notes that his sums for (C_2,A_1) and (G_2,A_1) are characters of L(2Λ_0) for D_3^{(2)}, D_4^{(3)} by Okado–Takenaka and
are modular by Kac–Peterson. Neither Sun–Wang (arXiv:2604.00847) nor Wang–Wang (arXiv:2512.07790) nor the later
papers on the conjecture that I could check mention this route. Consequences: the family (T_1, A_r) for **all** r
(odd r ≥ 5 was open in `T1-A-even.md`); Wang–Wang's concluding conjecture (½C(D_k), all k ≥ 3); Sun–Wang's rank-three
cases (C_3,A_1) and (B_3,A_1), listed as unresolved as late as September 2026 (Xia, arXiv:2609.08816), and
(T_3,A_1) (proved independently by Shi–Wang, arXiv:2607.23257); the pairs (A_{n−1},A_{k−1}) without the gcd restriction;
all (E_n, A_{k−1}). **Conditional on** the cited theorems exactly as quoted in §2; the weakest links are flagged
there: Takenaka's paper is an unrefereed preprint (needed for T_l only), and the Kac–Peterson statement for
twisted algebras could not be checked against the primary source — for B, C, F_4, G_2 it is reduced in §2 to the
untwisted statement by the S-transformation to the transposed algebra; **for A_{2l}^{(2)} (all pairs (T_l, A_{k−1}),
including (T_1, A_r)) the congruence-subgroup property of the string functions is a citation gap** (the identity
(5.1) with string functions does not depend on it; the word "modular" does). Proposition M in §2 (added after the
referee's report, not yet refereed) gives the classical argument for all affine types, which removes the
dependence on that citation.

## 0. The conjecture (verbatim) and notation

Sun–Wang, arXiv:2604.00847v1, §1 (generalized Nahm sum, after Mizuno): for D = diag(d_1,…,d_r) with positive
integers d_j and A ∈ Q^{r×r} with AD symmetric positive definite,

    f_{A,B,C,D}(τ) := Σ_{n∈N^r} q^{½ nᵀADn + nᵀB + C} / ( (q^{d_1};q^{d_1})_{n_1} ⋯ (q^{d_r};q^{d_r})_{n_r} ),

"The quadruple (A,B,C,D) is called modular if the associated generalized Nahm sum is modular", where modular
means (their §1) a modular form — necessarily of weight 0 — for a congruence subgroup of SL_2(Z).

> **Conjecture 1.1 (Sun–Wang).** "Let X and Y be Dynkin diagrams of type ABCDEFGT with Cartan matrices C(X) and
> C(Y), respectively. Let D(X) denote the diagonal matrix whose j-th diagonal entry is defined as the quotient
> α_j²/β² of the square lengths of the simple roots, where α_j is the j-th simple root and β is a short root.
> We set A(X,Y) = C(X)⊗C(Y)^{−1} and D(X,Y) = D(X)⊗D(Y) […]. We also set C(X,Y) = −c(X,Y)/24, where the central
> charge c(X,Y) is given by c(X,Y) = tr(D(X,Y))h(X)/(h(X)+h(Y)), and h(−) is the Coxeter number. Then
> (A(X,Y), 0, C(X,Y), D(X,Y)) is a modular quadruple."

For Y = A_{k−1} (D(Y) = 1, h(Y) = k, (C(Y)^{−1})_{st} = min(s,t) − st/k for 1 ≤ s,t ≤ k−1) this concerns

    F_{X,k}(q) := Σ_{p = (p_i^{(s)}) ∈ N^{l(k−1)}} q^{½ Σ_{i,j=1}^{l} Σ_{s,t=1}^{k−1} G_{ij} (min(s,t) − st/k) p_i^{(s)} p_j^{(t)}} / ∏_{i,s} (q^{d_i};q^{d_i})_{p_i^{(s)}},        (0.1)

where l = rank X, d_i = D(X)_{ii} and **G := C(X)D(X)** is the Gram matrix of the simple roots of X normalised so
that short roots have square length 2 (G_{ii} = 2d_i). The conjecture for the pair (X, A_{k−1}) says:

    q^{−c/24} F_{X,k}(q) is modular,        c = c(X, A_{k−1}) = tr(D(X))·(k−1)·h(X)/(h(X)+k).                       (0.2)

Sectors. Let Q̄ = ⊕_i Zα_i be the root lattice of X with the form (α_i|α_j) = G_{ij}. For p as in (0.1) put
μ(p) := Σ_i (Σ_s s p_i^{(s)}) α_i ∈ Q̄, and for a class c ∈ Q̄/kQ̄ let F_{X,k,c}(q) be the part of (0.1) with μ(p) ∈ c. Then
F_{X,k} = Σ_{c∈Q̄/kQ̄} F_{X,k,c}, and for X ≠ T_l the exponents of F_{X,k,c} lie in −(μ|μ)/(2k) + Z for μ ∈ c (because
½Σ G_{ij}min(s,t)p_i^{(s)}p_j^{(t)} ∈ Z and (μ|μ)/(2k) is well defined mod Z on c: G is even). For X = T_l (G_{ll} = 1) they lie in
−(μ|μ)/(2k) + ½Z; correspondingly the weights of L(kΛ_0) for A_{2l}^{(2)} have depths n ∈ ½Z (§5).

The affine algebra g(X). In Kac's notation (Table Aff; node 0 is the node with a_0 = 1 whose removal leaves X):

| X | A_l, D_l, E_l | C_l (l ≥ 2) | B_l (l ≥ 2) | F_4 | G_2 | T_l (l ≥ 1) |
|---|---|---|---|---|---|---|
| g(X) | X^{(1)} | A_{2l−1}^{(2)} (l ≥ 3), D_3^{(2)} (l = 2) | D_{l+1}^{(2)} | E_6^{(2)} | D_4^{(3)} | A_{2l}^{(2)} |
| r | 1 | 2 | 2 | 2 | 3 | 2 |
| h^∨(g(X)) = h(X) | h | 2l | 2l | 12 | 6 | 2l+1 |
| mult(jδ) | l | l (j even), l−1 (j odd) | l (j even), 1 (j odd) | 4 (j even), 2 (j odd) | 2 (3∣j), 1 (3∤j) | l |
| ∏_{j≥1}(1−q^j)^{mult(jδ)} | (q)_∞^l | (q)_∞^{l−1}(q²;q²)_∞ | (q)_∞(q²;q²)_∞^{l−1} | (q)_∞²(q²;q²)_∞² | (q)_∞(q³;q³)_∞ | (q)_∞^l |

In every column except the last: (i) the horizontal subalgebra (generated by e_i, f_i, i ≥ 1) has type X; (ii) Kac's normalised
invariant form restricted to the roots of X is the form G above (short roots of square length 2); (iii) the
dual Coxeter number of g(X) equals the Coxeter number of X; (iv) the lattice M of translations of the affine
Weyl group is Q̄ (Kac, *Infinite dimensional Lie algebras*, 3rd ed., §6.5: M = Q̄^∨ for r = 1, which is Q̄ for X simply
laced; M = Q̄ for r > 1, a_0 = 1); (v) the last row is **∏_{i=1}^{l} (q^{d_i};q^{d_i})_∞** — one factor for each simple root of X.
(B_2 = C_2 and D_3^{(2)} = A_3^{(2)} as algebras.) Here q = e^{−δ}, δ the primitive imaginary root. The last column (a_0 = 2)
is described in §5: there Q̄ is replaced by the lattice Q̄_T of finite parts of the weights of L(kΛ_0), whose Gram matrix
in a suitable basis is C(T_l) = C(T_l)D(T_l), and (iii), (v) hold as stated (h^∨(A_{2l}^{(2)}) = 2l+1 = h(T_l), all d_i = 1).

String functions (Kac §12.7; Kac–Peterson). For Λ dominant integral of level k and λ a weight of L(Λ),

    c^Λ_λ(τ) := q^{m_Λ(λ)} Σ_{n∈Z} mult_Λ(λ − nδ) q^n,     m_Λ(λ) = m_Λ − (λ|λ)/(2k),     m_Λ = (Λ+ρ|Λ+ρ)/(2(k+h^∨)) − (ρ|ρ)/(2h^∨).       (0.3)

c^Λ_λ depends only on the class of λ modulo kM + Cδ and on the W̄-orbit of that class. For Λ = kΛ_0 and
λ = kΛ_0 + μ, μ ∈ Q̄, write c_μ := c^{kΛ_0}_{kΛ_0+μ}; then m_{kΛ_0}(kΛ_0+μ) = m_{kΛ_0} − (μ|μ)/(2k) and (ρ = ρ̄ + h^∨Λ_0, (Λ_0|ρ̄) = (Λ_0|Λ_0) = 0)

    m_{kΛ_0} = (ρ̄|ρ̄)·( 1/(2(k+h^∨)) − 1/(2h^∨) ) = − k (ρ̄|ρ̄) / (2h^∨(k+h^∨)).                                                  (0.4)

## 1. Results

**Theorem 1 (X simply laced).** Let X ∈ {A_l (l ≥ 1), D_l (l ≥ 4), E_6, E_7, E_8}, k ≥ 2, g = X^{(1)}. Then

    q^{−c(X,A_{k−1})/24} F_{X,k}(q) = η(τ)^l Σ_{μ ∈ Q̄/kQ̄} c_μ(τ),                                                              (1.1)

the sum of all k^l string functions c_μ = c^{kΛ_0}_{kΛ_0+μ} of the level-k vacuum module (with repetitions). More
precisely, for every class c ∈ Q̄/kQ̄ and μ ∈ c:  q^{−c(X,A_{k−1})/24} F_{X,k,c}(q) = η(τ)^l c_μ(τ).

**Theorem 2 (X not simply laced).** Let X ∈ {C_l (l ≥ 2), B_l (l ≥ 2), F_4, G_2}, k ≥ 2, g = g(X) the twisted affine
algebra of the table. Then, with d_i = D(X)_{ii},

    q^{−c(X,A_{k−1})/24} F_{X,k}(q) = ∏_{i=1}^{l} η(d_iτ) · Σ_{μ ∈ Q̄/kQ̄} c_μ(τ),                                                    (1.2)

and sector by sector as in Theorem 1.

**Corollary 3.** Conjecture 1.1 of Sun–Wang holds for all pairs (X, A_{k−1}), X of type A, B, C, D, E, F, G, k ≥ 2:
q^{−c(X,A_{k−1})/24}F_{X,k} is a modular function for a congruence subgroup — given (KP-u) for X simply laced, and given (KP-t)
(i.e. Kac §13.9 and (KP-u) for the transposed algebra) for X = B_l, C_l, F_4, G_2; see §2. (The sectors F_{X,k,c} are modular as well, with
the same prefactor.) In particular:

* (D_k, A_1): f_{½C(D_k),0,(1−k)/24} is modular for every k ≥ 4 (and k = 3, D_3 = A_3) — the conjecture in the concluding
  remarks of Wang–Wang arXiv:2512.07790 ("The Nahm sum f_{½𝒞(D_k),0,(1−k)/24}(q) is modular for any k ≥ 3", proved there
  only for k = 3, 4 by Cao–Wang). Here c(D_k,A_1) = k(2k−2)/(2k) = k−1.
* (C_3, A_1) and (B_3, A_1): two of the seven rank-three cases listed as unproved by Sun–Wang (§1, list after
  item (h): (T_1,G_2), (G_2,T_1), (T_3,A_1), (C_3,A_1), (C_3,T_1), (B_3,A_1), (B_3,T_1)). According to Xia (arXiv:2609.08816,
  September 2026), after Shi–Wang ((T_3,A_1)) and Xia ((T_1,G_2), (G_2,T_1)) "the remaining unresolved cases are (C_3,A_1),
  (C_3,T_1), (B_3,A_1), (B_3,T_1)"; so after this note only (C_3,T_1) and (B_3,T_1) remain among Sun–Wang's rank ≤ 3 cases.
* (A_{n−1}, A_{k−1}) for all n, k ≥ 2: Sun–Wang quote Dunn's thesis for this under gcd((n−1)!, k) = 1.
* (E_6, A_{k−1}), (E_7, A_{k−1}), (E_8, A_{k−1}), (F_4, A_{k−1}), (G_2, A_{k−1}) for all k ≥ 2.

Each result is conditional on the cited theorems of §2 exactly as stated there (they are not re-derived here).

## 2. The cited theorems

**(LW) Lepowsky–Wilson** (*The structure of standard modules I*, Invent. Math. 77 (1984), Theorem A5.3; used in
this form by BKP §3.1 and OT (4.2)). Let ŝ ⊂ g be the homogeneous Heisenberg subalgebra (the sum of the
imaginary root spaces g_{jδ}, j ≠ 0, and the centre), ŝ^± its positive and negative parts. For a standard module
L(Λ) of level k > 0 the vacuum space L(Λ)^{ŝ⁺} := {v : ŝ⁺v = 0} satisfies U(ŝ⁻) ⊗ L(Λ)^{ŝ⁺} ≅ L(Λ) (isomorphism of graded
h̄-modules, x ⊗ v ↦ xv).

**(BKP) Butorac–Kožić–Primc**, *Parafermionic bases of standard modules for affine Lie algebras*, Math. Z. 298
(2021) 1003–1032, arXiv:2002.00435v2, §3.3, theorem with label `thm38` (Theorem 3.9 by the counters of the arXiv
source). Verbatim: "For any highest weight Λ as in (rect) we have ch L(Λ)^{ĥ⁺}_{Q(k)} = Σ_𝒫 D′_𝒫(q) C_𝒫(q) K_𝒫(q), where the sum
goes over all finite sequences 𝒫 = (𝒫_l,…,𝒫_1) of k_{α_1}+…+k_{α_l}−l nonnegative integers and D′_𝒫(q) is given by (charGij2)",
with K_𝒫(q) = q^{½Σ_{i,r}Σ_{m,n} K^{mn}_{ir} p_i^{(m)}p_r^{(n)}}, K^{mn}_{ir} = G^{mn}_{ir} − (mn/k)⟨α_i,α_r⟩, G^{mn}_{ir} = min{ν_r m, ν_i n}⟨α_i,α_r⟩, ν_i = k_{α_i}/k,
D′_𝒫 = 1/∏_i∏_{r=1}^{k_{α_i}−1}(q;q)_{p_i^{(r)}}, and C_𝒫 = 1 for Λ = kΛ_0. Here g is any untwisted affine algebra, Λ = k_0Λ_0 + k_jΛ_j is
"rectangular" (k = k_0 + k_j > 0; this includes kΛ_0), the parafermionic space is L(Λ)^{ĥ⁺}_{Q(k)} = L(Λ)^{ĥ⁺}/span{(ρ(kα)−1)v}
≅ ⊕_{μ∈Λ+Q/kQ} L(Λ)^{ĥ⁺}_μ in the simply-laced case (their (parafermion1)–(parafermion2)), and the character is
ch = q^{−c_Λ} tr q^{L_Ω(0)} with "ω_Ω = ω − ω_𝐡 […] the difference between the grading operators of the vertex operator
algebras L(kΛ_0) and M(k)" (their §3.3), c_{kΛ_0} = 0. Their Lemma in §3.3 shows [L_Ω(0), e_{β^∨}] = 0, i.e. the grading is
compatible with the identification of the sectors μ and μ + kβ. For simply laced g: k_{α_i} = k, ν_i = 1, ⟨α_i,α_r⟩ = C(X)_{ir},
so K^{mn}_{ir} = C(X)_{ir}(min(m,n) − mn/k) and the right-hand side is exactly (0.1): **ch L(kΛ_0)^{ĥ⁺}_{kQ} = F_{X,k}(q)**.
(For g = ŝl_{n+1} this is Georgiev, q-alg/9504024, (0.5); for D, E the input is Butorac–Kožić, arXiv:1902.10794.
The formula was conjectured by Kuniba–Nakanishi–Suzuki, hep-th/9301018.) The basis vectors of BKP are
h̄-weight vectors: a monomial with p_i^{(s)} quasi-particles of colour i and charge s has weight kΛ_0 + μ(p) modulo kQ;
so the part of the sum with μ(p) ∈ c is the character of the sector L(kΛ_0)^{ĥ⁺}_{μ}, μ ∈ c (Georgiev (5.8) for ŝl_{n+1}).

**(OT) Okado–Takenaka**, *Parafermionic bases of standard modules for twisted affine Lie algebras of type
A_{2l−1}^{(2)}, D_{l+1}^{(2)}, E_6^{(2)} and D_4^{(3)}*, arXiv:2109.08892v1; Algebr. Represent. Theory (2022), DOI 10.1007/s10468-022-10145-2
(Theorem 16 of the journal version; a published Correction, DOI 10.1007/s10468-022-10159-w, was not accessible to me —
the formula as quoted here is the arXiv one and is confirmed numerically in every family, §6). arXiv version: last
theorem of §5.2 (it has no label; it follows their Theorem `th:BS2`). Verbatim: "For each of the affine Lie algebras A_{2l−1}^{(2)}, D_{l+1}^{(2)}, E_6^{(2)}, D_4^{(3)}, we have

    ch L(kΛ_0)_{kQ}^{ĥ[ν]⁺} = Σ_{𝒫^{(k−1)}} q^{½Σ_{i,j=1}^{l}⟨(α_i)_{(0)},(α_j)_{(0)}⟩ Σ_{s,t=1}^{k−1} D^{(k)}_{s,t} p_i^{(s)}p_j^{(t)}} / ∏_{i=1}^{l}∏_{s=1}^{k−1}(q^{ρ_i})_{p_i^{(s)}} · ∏_{i=1}^{l} y_i^{Σ_{s=1}^{k−1} s p_i^{(s)}}

where the sum runs over all sequences 𝒫^{(k−1)} of l(k−1) nonnegative integers and D^{(k)}_{s,t} = min{s,t} − st/k."
Their conventions: g is simply laced of type A_{2l−1}, D_{l+1}, E_6, D_4 with form ⟨α,α⟩ = 2, ν the diagram automorphism of
order r, x_{(0)} the projection to the fixed space h_{(0)}, ρ_i = ½⟨(α_i)_{(0)},(α_i)_{(0)}⟩ (their (2.x), label `rho`), the twisted
affine algebra is ĝ[ν] = ⊕_{m∈(1/r)Z} g_{(rm)}⊗t^m ⊕ Cc with [d, x⊗t^m] = m x⊗t^m; the vacuum space is taken for
ĥ[ν]⁺ (the positive part of the full twisted Heisenberg algebra); the parafermionic space is the space of
kQ-coinvariants, "L(kΛ_0)_{kQ}^{ĥ[ν]⁺} ≃ ⊕_{μ∈kΛ_0+Q/kQ} L(kΛ_0)_μ^{ĥ[ν]⁺}"; and the grading operator is (their §5.1)
"D = −d − D^{ĥ[ν]⁺}, D^{ĥ[ν]⁺}|_{L(kΛ_0)^{ĥ[ν]⁺}_μ} = ⟨μ_{(0)},μ_{(0)}⟩/(2k)"; the y_i record the h_{(0)}-weight. The principal-subspace
input is Butorac–Sadowski, New York J. Math. 25 (2019) 71–106. The formula was conjectured by
Hatayama–Kuniba–Okado–Takagi–Tsuboi (math/0102113, Conjecture 5.3).

**(KP) Kac–Peterson**, *Infinite-dimensional Lie algebras, theta functions and modular forms*, Adv. Math. 53
(1984) 125–264, Theorem A and §4; Kac's book, §12.7 and §13.8–13.10. Statement used: for an affine Lie algebra g and
a dominant integral Λ of level k > 0, every string function c^Λ_λ of (0.3) is a modular form of weight −½ dim h̄ for a
congruence subgroup of SL_2(Z) (holomorphic on the upper half-plane, meromorphic at the cusps, with a multiplier
system of finite order). **The primary sources were not accessible in this session**, so the three cases have
different status:

* *(KP-u) g untwisted* (needed for X simply laced): this is Kac–Peterson's Theorem A as it is universally quoted
  (e.g. Kac–Wakimoto, Acta Appl. Math. 21 (1990), §0.6: the string function "is a modular form of weight −ℓ/2",
  transformation law "[Kac–Peterson, Theorem A]"; Kac's book, Corollary 13.10). Standard; taken as cited.
* *(KP-t) g twisted, not A_{2l}^{(2)}* (needed for X = B_l, C_l, F_4, G_2): Kac–Wakimoto (loc. cit.) Theorem 4.3: "Let Γ = Γ_0(r) if g
  is of type X_N^{(r)} ≠ A_{2l}^{(2)} and Γ = Γ_θ if it is of type A_{2l}^{(2)}. The following statement is proved in [Kac–Peterson] for r = 1
  and is implicitly contained there in the general case. Theorem 4.3. The C-span of the set of normalized characters
  {χ_Λ} is Γ-invariant", and their Remark 4.3 / Kac §13.9: the transformation S maps the characters of g of level k to
  those of the adjacent (transposed) untwisted algebra g′ = B_l^{(1)}, C_l^{(1)}, F_4^{(1)}, G_2^{(1)} for g = A_{2l−1}^{(2)}, D_{l+1}^{(2)}, E_6^{(2)}, D_4^{(3)}
  (with a rescaling of τ by r). Since χ_Λ = Σ_λ c^Λ_λ Θ_λ with classical theta functions Θ_λ of an even lattice (linearly
  independent over functions of τ, and transforming among theta functions under S), each string function of g at τ is
  a linear combination of string functions of g′ at −1/τ (suitably rescaled), which are modular for a principal
  congruence subgroup by (KP-u) for g′; hence so is c^Λ_λ. This argument is complete **given** Kac §13.9 and (KP-u) for
  the non-simply-laced untwisted algebras; Γ-invariance of the span alone would not suffice (a finite-dimensional
  representation need not have finite image). Mizuno (arXiv:2305.02267) uses Kac–Peterson for D_3^{(2)}, D_4^{(3)} in the same way:
  "Then the modularity follows from the result by Kac and Peterson".
* *(KP-A) g = A_{2l}^{(2)}* (needed for X = T_l): here S preserves the span of the characters (Kac's Theorem 13.8 treats
  A_{2l}^{(2)} together with the untwisted types; Γ_θ-invariance: Kac–Wakimoto Theorem 4.3; Iohara–Saito arXiv:2502.03316 §1:
  "For BC_l^{(2)} (or A_{2l}^{(2)}), this congruence subgroup is the so-called θ-subgroup Γ_θ"), and the argument of (KP-t) is not
  available. That the string functions are modular **for a congruence subgroup** then requires Kac–Peterson's explicit
  matrices for S and T² (a Weil-type representation with congruence kernel). I believe this is in Kac–Peterson §4,
  but I could not verify it: **citation gap, bypassed by Proposition M below**. An alternative route (hypotheses not checked by me): the level-k
  A_{2l}^{(2)}-modules are ν-twisted modules of the rational, C_2-cofinite vertex operator algebra L_k(sl_{2l+1}); modular
  invariance of twisted trace functions (Dong–Li–Mason, Comm. Math. Phys. 214 (2000)) and the congruence property in
  orbifold theory (Dong–Ren, 2018) give a congruence subgroup for the characters, and the theta-function argument
  passes to string functions. Numerically (`code/xa_modularity_numeric.py`, and the referee's free fits) the sums for
  (T_1,A_3), (T_1,A_5), (T_3,A_1) behave like weight-0 modular functions with exactly the predicted prefactor; this
  cannot distinguish congruence from non-congruence subgroups.

**Proposition M (modularity of string functions; classical argument, valid for every affine type — added after
the referee's report, not itself refereed).** Let g be an affine Lie algebra, l = dim h̄, Λ dominant integral of level
k > 0, λ a weight of L(Λ). Then c^Λ_λ is holomorphic on the upper half-plane, and there are N ≥ 1 and a character ε of
Γ(N), trivial on a principal congruence subgroup, such that c^Λ_λ(γτ) = ε(γ) v_θ(γ)^{−l} (cτ+d)^{−l/2} c^Λ_λ(τ) for γ ∈ Γ(N) (v_θ the multiplier of
the Jacobi theta function θ(τ) = Σq^{n²/2}, N divisible by 4); c^Λ_λ is meromorphic at the cusps. Consequently ∏_iη(d_iτ)·c^Λ_λ(τ) is
a modular function for a congruence subgroup.
*Proof (this is the argument of Kac–Peterson, reduced to classical facts).* (1) W = W̄ ⋉ t(M) with a lattice M ⊂ h̄*, and for
μ of level m and γ ∈ M: t_γ(μ) = μ + mγ − ((μ̄|γ) + ½m|γ|²)δ (Kac (6.5.2)); (·|·) is positive definite and rational on M, with
(P|M) ⊂ Z. (2) For μ of level m > 0 put Θ_μ := e^{−|μ|²δ/(2m)}Σ_{γ∈M}e^{t_γ(μ)} = e^{mΛ_0}Σ_{γ∈M+μ̄/m}e^{mγ−½m|γ|²δ}; in the coordinates e^{−δ} = e^{2πiτ}, e^{β} = e^{2πi(β|z)}
(β ∈ h̄*, z ∈ h̄), e^{Λ_0} = e^{2πit} this is the classical theta series with characteristic of the lattice (M, m(·|·)):
Θ_μ(τ,z,t) = e^{2πimt}Σ_{γ∈M+μ̄/m}e^{πimτ|γ|²+2πim(γ|z)}. (3) Splitting the sum over W in the Weyl–Kac character formula and in the
denominator identity along W̄ ⋉ t(M) gives e^{−m_Λδ}ch L(Λ) = A_{Λ+ρ}/A_ρ, A_μ := Σ_{w∈W̄}ε(w)Θ_{w(μ)} (Kac (12.7.?); this is where m_Λ of (0.3)
comes from), and, since mult_Λ is W-invariant, e^{−m_Λδ}ch L(Λ) = Σ_{λ mod kM+Cδ} c^Λ_λ(τ)Θ_λ(τ,z,t) (Kac (12.7.12)). (4) Classical
transformation formula for theta series of a positive definite rational lattice with rational characteristics
(Schoeneberg, *Elliptic Modular Functions*, Ch. IX; for the action of the full modular group on the span, Kac
Theorem 13.5): there is N_m such that for all γ = (a b; c d) ∈ Γ(N_m) and all μ of level m with μ̄ ∈ P̄,
Θ_μ(γτ, z/(cτ+d), t − c(z|z)/(2(cτ+d))) = ψ_m(d) v_θ(γ)^{l} (cτ+d)^{l/2} Θ_μ(τ,z,t) with a Dirichlet character ψ_m (the same for all μ). (5) By (3), (4)
the normalised character χ_Λ = A_{Λ+ρ}/A_ρ satisfies χ_Λ(γ·(τ,z,t)) = (ψ_{k+h^∨}/ψ_{h^∨})(d) χ_Λ(τ,z,t) on Γ(N_{k+h^∨}) ∩ Γ(N_{h^∨}), the action on (z,t) being
the one of (4) (the weights and the theta multipliers cancel; the exponential factors combine to level k). (6) ch L(Λ)
is holomorphic in z (Kac, Ch. 11), and the Θ_λ, λ mod kM + Cδ, are linearly independent over functions of τ (they
have disjoint supports in the characters e^{2πi(β|z)}). Comparing χ_Λ = Σ_λ c^Λ_λΘ_λ at γ·(τ,z,t) and at (τ,z,t) gives the
transformation law of c^Λ_λ with ε = ψ_{k+h^∨}ψ_{h^∨}^{−1}ψ_k^{−1}. (7) Holomorphy: mult_Λ(λ−nδ) is bounded by the coefficients of ∏(1−q^n)^{−dim g̊}.
Cusps: for γ ∈ SL_2(Z) (or Γ_θ) the same computation with the transformation of the spans of theta functions shows that
c^Λ_λ|γ is again a theta-decomposition coefficient of a quotient of theta series, hence has a q-expansion with finitely
many negative exponents. Finally η(d_iτ) has the multiplier of a theta function of weight ½ (η(24τ) is a theta series),
so the product ∏η(d_iτ)c^Λ_λ has weight 0 and a multiplier that is a character trivial on a principal congruence
subgroup. ∎
For A_{2l}^{(2)} the lattice M is the odd lattice Q̄_T of §5 (Gram matrix C(T_l)); the argument is the same (this is why the
full modular group is replaced by Γ_θ in Kac–Wakimoto's Theorem 4.3), and the conclusion — invariance under some
principal congruence subgroup — is unaffected. With Proposition M, the modularity statements of Corollary 3,
Theorem 3 and Corollary 5 depend only on the Weyl–Kac formula and classical theta functions, not on the wording
of Kac–Peterson's Theorem A.

**(SF) The strange formula** (Freudenthal–de Vries; Kac (12.1.8)): for a simple Lie algebra of type X with the
form in which long roots have square length 2, (ρ̄|ρ̄) = h^∨(X) dim(X)/12.

## 3. Proofs

### 3.1 Sectors of the vacuum space are strings

Let g be affine with a_0 = 1, L = L(kΛ_0), ŝ the homogeneous Heisenberg subalgebra. ŝ commutes with h̄ and ŝ⁻
is graded by −d with dim (ŝ⁻)_j = mult(jδ), so by (LW), for μ ∈ Q̄,

    Σ_n mult_{kΛ_0}(kΛ_0 + μ − nδ) q^n = ( Σ_n dim (L^{ŝ⁺})_{kΛ_0+μ−nδ} q^n ) / ∏_{j≥1}(1−q^j)^{mult(jδ)}.                              (3.1)

Hence, by (0.3) and property (v) of the table (∏_j(1−q^j)^{mult(jδ)} = ∏_i(q^{d_i};q^{d_i})_∞ = q^{−tr D(X)/24}∏_iη(d_iτ)):

    ∏_i η(d_iτ) · c_μ(τ) = q^{tr D(X)/24 + m_{kΛ_0}} · q^{−(μ|μ)/(2k)} Σ_n dim (L^{ŝ⁺})_{kΛ_0+μ−nδ} q^n.                                        (3.2)

### 3.2 The prefactor

**Lemma 4.** For every X of the table, with g = g(X): (ρ̄|ρ̄)/(2h^∨) = tr D(X)·(h(X)+1)/24, where (·|·) is Kac's normalised form
of g (short roots of X of square length 2), and consequently

    tr D(X)/24 + m_{kΛ_0} = − c(X, A_{k−1})/24.                                                                                 (3.3)

*Proof.* For X simply laced the form is the standard one, h^∨ = h, tr D = l and dim X = l(h+1), so (SF) gives
(ρ̄|ρ̄)/(2h) = l(h+1)/24. For X not simply laced the form is r times the standard form of X (long roots have square
length 2r), so (ρ̄|ρ̄) = r·h^∨(X)dim(X)/12 by (SF), and h^∨(g(X)) = h(X):

| X | r | h^∨(X) | dim X | (ρ̄|ρ̄) = r h^∨(X) dim X/12 | h(X) = h^∨(g) | (ρ̄|ρ̄)/(2h^∨(g)) | tr D(X) (h+1)/24 |
|---|---|---|---|---|---|---|---|
| C_l | 2 | l+1 | l(2l+1) | l(l+1)(2l+1)/6 | 2l | (l+1)(2l+1)/24 | (l+1)(2l+1)/24 |
| B_l | 2 | 2l−1 | l(2l+1) | l(2l−1)(2l+1)/6 | 2l | (2l−1)(2l+1)/24 | (2l−1)(2l+1)/24 |
| F_4 | 2 | 9 | 52 | 78 | 12 | 78/24 | 6·13/24 |
| G_2 | 3 | 4 | 14 | 14 | 6 | 28/24 | 4·7/24 |

(tr D(C_l) = (l−1)+2, tr D(B_l) = 2(l−1)+1, tr D(F_4) = 6, tr D(G_2) = 4.) Then by (0.4), with h = h(X) = h^∨(g):
tr D/24 + m_{kΛ_0} = (tr D/24)(1 − k(h+1)/(k+h)) = −tr D·h(k−1)/(24(h+k)) = −c(X,A_{k−1})/24. ∎

So the central charge of Sun–Wang's conjecture for (X, A_{k−1}) is exactly the Kac–Peterson modular anomaly of the
vacuum module of g(X) corrected by the Heisenberg factor.

### 3.3 Proof of Theorem 1

By (BKP) with Λ = kΛ_0 (c_Λ = 0), F_{X,k}(q) = tr q^{L_Ω(0)} on L(kΛ_0)^{ĥ⁺}_{kQ} ≅ ⊕_{μ∈Q/kQ} L(kΛ_0)^{ĥ⁺}_{kΛ_0+μ}, and F_{X,k,c} is the trace on the
summand μ ∈ c. On L(kΛ_0), L(0) = −d (the Sugawara operator and −d agree on the vacuum module, both vanish on
the highest weight vector and have the same commutators with g), and L_𝐡(0) acts on a vector of the vacuum
space of h̄-weight μ as (μ|μ)/(2k) (the modes h(n), n > 0, kill it, and h(0) acts by μ). So on (L^{ŝ⁺})_{kΛ_0+μ−nδ}:
L_Ω(0) = n − (μ|μ)/(2k), i.e.

    F_{X,k,c}(q) = q^{−(μ|μ)/(2k)} Σ_n dim (L^{ŝ⁺})_{kΛ_0+μ−nδ} q^n        (μ ∈ c).                                                         (3.4)

Insert (3.4) into (3.2) and use (3.3): η^l c_μ = q^{−c/24}F_{X,k,c}. Summing over Q̄/kQ̄ gives (1.1). ∎

### 3.4 Proof of Theorem 2

Let g = ĝ[ν] ⊕ Cd as in (OT), of type g(X) with X = C_l, B_l, F_4, G_2 for A_{2l−1}^{(2)}, D_{l+1}^{(2)}, E_6^{(2)}, D_4^{(3)}. Dictionary:

* The imaginary root spaces are g_{jδ} = h_{(j)}⊗t^{j/r}, so the homogeneous Heisenberg algebra ŝ is ĥ[ν], the primitive
  imaginary root δ corresponds to the energy 1/r, and e^{−δ} = q_OT^{1/r}: **q_OT = q^r**.
* The simple roots of the horizontal algebra g_{(0)} (type X) are the projections (α_i)_{(0)}, i = 1,…,l (one for each
  ν-orbit of simple roots), and r⟨(α_i)_{(0)},(α_j)_{(0)}⟩ = G_{ij}: an orbit of size r of mutually orthogonal roots projects to a
  vector of square length 2/r (short root, d_i = 1), a fixed root has square length 2 (long, d_i = r), and the
  off-diagonal entries are r·⟨(α_i)_{(0)},α_j⟩ = −(number of bonds). Hence r·ρ_i = G_{ii}/2 = d_i, and (·|·) = r⟨·,·⟩ on h_{(0)}^*.

Substituting q_OT = q^r, y_i = 1 in (OT): the exponent becomes ½ΣG_{ij}D^{(k)}_{st}p_i^{(s)}p_j^{(t)} and (q_OT^{ρ_i};q_OT^{ρ_i}) becomes (q^{d_i};q^{d_i}), so the
right-hand side is F_{X,k}(q) of (0.1), and the part with μ(p) ∈ c is the trace over the sector μ ∈ c. On the
other side, q_OT^{D} = q^{−rd − r⟨μ_{(0)},μ_{(0)}⟩/(2k)} = q^{n − (μ|μ)/(2k)} on (L^{ŝ⁺})_{kΛ_0+μ−nδ}. So (3.4) holds verbatim, and (3.2), (3.3) give
∏_iη(d_iτ)c_μ = q^{−c/24}F_{X,k,c}; summing over Q̄/kQ̄ gives (1.2). ∎

### 3.5 Proof of Corollary 3

By (KP-u), resp. (KP-t), each c_μ is a modular form of weight −l/2 on a congruence subgroup; ∏_iη(d_iτ) is a modular form
of weight l/2 on a congruence subgroup (with a multiplier system of finite order); the right-hand sides of (1.1), (1.2)
are finite sums of such products. ∎

Remarks. (0) For X not simply laced Sun–Wang's c(X,A_{k−1}) is the exponent of the modular prefactor but **not** the
effective central charge governing the growth as q → 1: numerically c_eff(C_2,A_1) = 1 (c = 2), c_eff(G_2,A_1) = 1 (c = 3),
c_eff(C_3,A_1) = 15/8 (c = 3), c_eff(B_3,A_1) = 3/2 (c = 15/4) (`certificates/xa_modularity_numeric.txt`; referee's fits). This
is consistent with (KP-t): the S-image of a twisted string function lives on the adjacent algebra g′. For X simply
laced and X = T_l, c_eff = c. (a) The number of distinct functions among the c_μ is the number of orbits of W̄ ⋉ (translations by kQ̄)
on Q̄/kQ̄; e.g. for (A_1,A_5) Sun–Wang's observed χ_0+χ_{3/2}+2χ_{5/6}+2χ_{4/3} is Σ_{μ∈Z/6} with the orbits {0},{3},{±1},{±2}.
(b) For X simply laced the sector functions η^l c_μ are characters of modules of the parafermion vertex operator
algebra K(g,k) (the commutant of the Heisenberg algebra in L_k(g)). For X not simply laced I expect (not checked,
not used) that they are characters of ν-twisted modules of K(g,k) for the simply laced g of (OT).
(c) Zagier/Mizuno duality sends (X, A_{k−1}) to (A_{k−1}, X) (central charges add up to tr D(X)(k−1)); the dual pairs
(A_{p−1}, X) are the coset sums treated by Cherednik–Feigin, arXiv:1209.1978 (Theorem 1.2/Corollary 1.3 and the
modularity statement of their §2: sums with exponent ½(b_1² + (b_1−b_2)² + … + b_{p−1}²), b_i ∈ P_−), for all reduced root systems,
up to an unspecified power of q. So both members of each dual pair {(X,A_{k−1}), (A_{k−1},X)} are covered by
representation theory, by two different mechanisms (parafermions; level-one cosets).

## 4. What exactly is new, and attribution

* The character formulas are due to Kuniba–Nakanishi–Suzuki (conjecture), Lepowsky–Primc (ŝl_2), Georgiev (ŝl_n),
  Butorac–Kožić–Primc (all untwisted types), Hatayama–Kuniba–Okado–Takagi–Tsuboi (twisted conjecture),
  Okado–Takenaka (twisted types except A_{2l}^{(2)}). The identification of sector characters with η^l·(string
  function) is stated in KNS (their formula for ch(Ω^Λ_λ)), Georgiev (5.9)–(5.11) and Gepner (arXiv:1410.1778).
  Modularity of string functions is Kac–Peterson.
* Not found in the literature (searches: `PRIOR-ART.md`): the statement that the *unrestricted* sums
  C(X)⊗C(A_{k−1})^{−1} with (q^{d_i};q^{d_i}) denominators are therefore modular with the prefactor q^{−tr D·h(k−1)/(24(h+k))}; the
  observation that Sun–Wang's symmetrizable pairs (B_l, C_l, F_4, G_2; A_{k−1}) are the *twisted* parafermion
  characters (Lemma 4: the uniform formula for the anomaly, and tr D(X) as the weight of the Heisenberg factor);
  the consequences listed in Corollary 3. Sun–Wang know the A_1 case ((A_1,A_r) = Z_{r+1} parafermions, KKMM) and
  quote Dunn for (A_{n−1},A_{k−1}) under a gcd condition; Wang–Wang motivate their conjecture by duality only.
* Prior use of the mechanism: Mizuno, *Remarks on Nahm sums for symmetrizable matrices* (arXiv:2305.02267; Ramanujan J.
  66 (2025)), §4 ("Experiments on modular Nahm sums"), comments on Table 1 (rank 2): "A = [[1,−1/2],[−1,1]] or [[1,−3/2],[−1,1]]. If
  b = 0, the Nahm sum coincides, up to a multiplication by a product of the dedekind eta function, with a certain
  character of the integrable highest weight module L(2Λ_0) of the affine Lie algebra of type D_3^{(2)} or D_4^{(3)}. Such a
  character formula was conjectured in [HKOTT] and recently proved in [Okado–Takenaka]. Then the modularity follows
  from the result by Kac and Peterson". These are the pairs (C_2,A_1), (G_2,A_1). [The second matrix as printed is not
  symmetrizable with positive definite AD; presumably [[1,−3/2],[−1/2,1]] with d = (3,1) is meant (referee's remark).]
  So the idea is Mizuno's for two rank-two cases at k = 2.
* For X simply laced the statement is literally Kuniba–Nakanishi–Suzuki's (hep-th/9301018): sector formula
  "ch(Ω^0_λ) = q^{−c_PF/24} Σ_{λ(n)≡λ mod ℓQ^∨} q^{𝒦(n)}/(q)_n", "ch(Ω^Λ_λ) = η(q)^r c^Λ_λ(q)", and "The character of the space Ω̄^0 is now given by
  the same formula but *without* any restriction on the n-sum"; with 𝒦(n) = ½nᵀ(C(X)⊗C(A_{ℓ−1})^{−1})n and c_PF = ℓ dim X/(ℓ+h) − r = c(X,A_{ℓ−1}).
  Nahm, *Conformal field theory and torsion elements of the Bloch group* (hep-th/0404120), section on Dynkin diagrams,
  identifies the pairs (Y, A_{k−1}) with the parafermion coset (the complement of the Heisenberg algebra at level k)
  and calls the character identification "known or conjectured … probably true in general". So Theorem 1 is the
  observation that the conjecture of KNS is now a theorem (Georgiev; Butorac–Kožić–Primc) and that this settles
  Sun–Wang's Conjecture 1.1 for (X, A_{k−1}), X simply laced — in particular Wang–Wang's conjecture. I found no paper
  that says so; Sun–Wang quote (A_{n−1},A_{k−1}) only under Dunn's gcd condition.
* New bookkeeping here: the identification of the pairs (B_l, C_l, F_4, G_2, T_l; A_{k−1}) in Sun–Wang's conventions with
  the twisted parafermion characters of Okado–Takenaka and Takenaka (the conjecture of
  Hatayama–Kuniba–Okado–Takagi–Tsuboi gives the string functions only "up to a power of q"), Lemma 4 and §5 (the power
  of q is exactly −c(X,A_{k−1})/24, uniformly), and the consequences (C_3,A_1), (B_3,A_1), (T_l,A_{k−1}).
* I regard Theorems 1–3 as "known to experts in a different language, not recorded for these conjectures".

## 5. The tadpole pairs (T_l, A_{k−1}) and A_{2l}^{(2)}

Sun–Wang's T_l has C(T_l) = the A_l Cartan matrix with last diagonal entry 1, D = 1, h(T_l) = 2l+1; C(T_l) is the Gram
matrix of the vectors e_1−e_2, …, e_{l−1}−e_l, e_l. For A_{2l}^{(2)} in Kac's labelling (a = (2,…,2,1), a_0^∨ = 1; α_0 is the
short end node, the normalised form has (α_0|α_0) = 1, (α_i|α_i) = 2 (0 < i < l), (α_l|α_l) = 4; h^∨ = 2l+1; mult(jδ) = l
for all j ≥ 1): the weights of L(kΛ_0) are kΛ_0 + μ − nδ with n ∈ ½Z and μ in the lattice Q̄_T with basis
α_1,…,α_{l−1}, ½α_l (indeed kΛ_0 − Σβ_iα_i = kΛ_0 − (β_0/2)δ + Σ_{0<i<l}(β_0−β_i)α_i + (β_0−2β_l)·½α_l), whose Gram matrix is C(T_l).

**(Tak) Takenaka**, *Vertex algebraic construction of modules for twisted affine Lie algebras of type A_{2l}^{(2)}*,
arXiv:2205.05271v1 (11 May 2022; I found no journal reference), theorem with label `character2`. Verbatim: "For
affine Lie algebras A_{2l}^{(2)}, we obtain

    ch L(kΛ_0)_{kQ}^{ĥ[ν]⁺} = Σ_{𝒫^{(k−1)}} q^{½Σ_{i,j=1}^{l}⟨(α_i)_{(0)},(α_j)_{(0)}⟩ Σ_{s,t=1}^{k−1} D^{(k)}_{s,t} p_i^{(s)}p_j^{(t)}} / ∏_{i=1}^{l}∏_{s=1}^{k−1}(q^{1/2})_{p_i^{(s)}} · ∏_{i=1}^{l} y_i^{Σ_{s=1}^{k−1} s p_i^{(s)}}

where the sum runs over all sequences 𝒫^{(k−1)} of l(k−1) nonnegative integers and D^{(k)}_{s,t} = min{s,t} − st/k", with
(q^{1/2})_n = ∏_{i=1}^{n}(1−q^{i/2}) [the paper prints ∏(1−q^{1/2}), an evident misprint: it is defined through 1/(q^{1/2})_n = Σ_j p_n(j)q^{j/2}],
g = sl_{2l+1}, ν the diagram automorphism α_i ↦ α_{2l+1−i} (lifted to an automorphism of order 4 of the lattice vertex
algebra), the same definitions of vacuum space, parafermionic space and grading operator D = −d − ⟨μ_{(0)},μ_{(0)}⟩/(2k) as in
(OT), and (their theorem `ST ch`) ch L(kΛ_0) = ∏_i(q^{1/2})_∞^{−1} · ch L(kΛ_0)^{ĥ[ν]⁺} by Lepowsky–Wilson. The abstract states that
this settles the conjecture of Hatayama–Kuniba–Okado–Takagi–Tsuboi for the vacuum modules of A_{2l}^{(2)}.

**Theorem 3 (X = T_l; conditional on (Tak) and (KP) for A_{2l}^{(2)}).** For l ≥ 1 and k ≥ 2,

    q^{−c(T_l,A_{k−1})/24} F_{T_l,k}(q) = η(τ)^l Σ_{μ ∈ Q̄_T/kQ̄_T} c_μ(τ),        c(T_l,A_{k−1}) = l(k−1)(2l+1)/(2l+1+k),                         (5.1)

where c_μ := c^{kΛ_0}_λ for any weight λ of L(kΛ_0) with λ ≡ kΛ_0 + μ mod Cδ are the string functions of A_{2l}^{(2)} at level k
(explicitly c_μ = q^{m_{kΛ_0} − (μ|μ)/(2k)}Σ_{n∈½Z} mult_{kΛ_0}(kΛ_0+μ−nδ)q^n; if the ½α_l-coordinate of μ is odd the depths n lie in ½+Z and one
takes λ = kΛ_0 + μ − ½δ in Kac's definition (0.3)); and sector by sector, for μ ∈ c ∈ Q̄_T/kQ̄_T:
q^{−(μ|μ)/(2k)} Σ_{n∈½Z} mult_{kΛ_0}(kΛ_0 + μ − nδ) q^n · (q)_∞^l = F_{T_l,k,c}(q). Consequently Sun–Wang's Conjecture 1.1 holds for (T_l, A_{k−1})
**provided** the string functions of A_{2l}^{(2)} are modular for a congruence subgroup ((KP-A) in §2: a citation I could not
verify; Proposition M in §2 gives the classical argument). The identity (5.1) itself depends only on (Tak).

*Proof.* Dictionary as in §3.4. The imaginary root spaces are h_{(j)}⊗t^{j/2} (dim h_{(0)} = dim h_{(1)} = l), so δ ↔ energy ½ and
q_Tak = q². The projections (α_i)_{(0)} = ½(α_i + α_{2l+1−i}), i = 1,…,l, have Gram matrix ½C(T_l) (square lengths 1,…,1,½;
neighbours −½), they are identified with the basis α_1,…,α_{l−1},½α_l of Q̄_T, and (·|·) = 2⟨·,·⟩. Substituting q_Tak = q², y_i = 1 in
(Tak) gives F_{T_l,k}(q) (exponent ½ΣC(T_l)_{ij}D^{(k)}_{st}p_i^{(s)}p_j^{(t)}, denominators (q;q)), and q_Tak^{D} = q^{n − (μ|μ)/(2k)} on the vectors of the
vacuum space of weight kΛ_0 + μ − nδ; so (3.4) holds. (3.1), (3.2) hold with ∏_j(1−q^j)^{mult(jδ)} = (q)_∞^l = q^{−l/24}η^l. For the
prefactor: by (0.4) with h^∨ = 2l+1 and (ρ̄|ρ̄) = l(l+1)(2l+1)/6 (ρ̄ for the horizontal algebra C_l with long roots of square
length 4, as in Lemma 4), m_{kΛ_0} = −k l(l+1)/(12(k+2l+1)), and l/24 + m_{kΛ_0} = (l/24)(1 − 2k(l+1)/(k+2l+1)) = −l(2l+1)(k−1)/(24(k+2l+1))
= −c(T_l,A_{k−1})/24. Modularity: (KP-A), see §2. ∎

**Corollary 5 (the family (T_1, A_r)).** For every r ≥ 1, with k = r+1,

    q^{−r/(8(r+4))} Σ_{n∈N^r} q^{½ nᵀC(A_r)^{−1}n} / ∏_i (q)_{n_i} = η(τ) Σ_{m mod r+1} c^{(r+1)Λ_0}_{(r+1)Λ_0 + mε}(τ)      (ε = ½α_1, string functions of A_2^{(2)}),

a modular function given (KP-A) (c = 3r/(r+4)); the summand m corresponds to the part of the Nahm sum with Σ_s s·n_s ≡ m mod r+1.
For even r this is consistent with Theorem A of `T1-A-even.md` (proved there without (Tak): the same function
equals Σ_N c^0_{2N}(τ)q^{N²/(k′(2k′+1))} with ŝl_2 string functions of level k′ = r/2); for odd r ≥ 5 it is the missing case of the
(T_1, X) family, now conditional on (Tak) and (KP). With Wang–Wang (arXiv:2607.08606, X = D_r) and `E8-E6.md`, `E7.md`
(X = E_6, E_7, E_8) this gives Sun–Wang's conjecture for (T_1, X) for every simply laced X; attribution table in
`RESULTS.md`.

Status of the inputs for X = T_l: (Tak) is an arXiv preprint (v1 only); its statement is confirmed numerically
by `certificates/xa_parafermion_certify.txt` for l = 1, k ≤ 8; l = 2, k ≤ 4; l = 3, k ≤ 3; l = 4, k = 2 (all sectors). For
l = 1 and r = k−1 even the statement also follows from `T1-A-even.md` only as far as the total sum is concerned.
(T_3,A_1) was proved independently by Shi–Wang (arXiv:2607.23257), by a different method.

## 6. Numerical certificate

`code/xa_parafermion_certify.py` computes weight multiplicities of L(kΛ_0) from the generalized Cartan matrix
alone (Racah recursion Σ_{w∈W} ε(w) mult(λ+ρ−wρ) = δ_{λ,Λ} for dominant λ, the affine Weyl group generated by
breadth-first search), multiplies the string of each μ by ∏_j(1−q^j)^{mult(jδ)} and compares with q^{(μ|μ)/(2k)}F_{X,k,c}
computed by direct summation of (0.1); for every class c ∈ Q̄/kQ̄ and 2^l representatives μ ∈ c. Cases and depths are
in the output file; all agree. This checks (3.4), the dictionary of §3.4 and (T); it is not a proof.

Small cases against known results: (A_1,A_1): F = (−q^{1/2};q)_∞, c = 1/2; (A_2,A_1): Zagier's rank-two example
A = ((1,−½),(−½,1)), c = 6/5; (C_2,A_1): q^{−1/12}F = η-quotient (Euler exponents of F in q^{1/2} have period 12:
2,−1,0,1,2,−2,2,1,0,−1,2,0; total weight 0 and leading exponent −1/12, in agreement with c = 2; this case is in
Mizuno's Table 1, proved by Wang–Wang).
