# The Nahm sum of C(A_{2k})⁻¹ (the pair (T1,A_{2k})): Z_k parafermions times a boson

> **Referee notes (2026-11-26, `REFEREE-T1-ABC.md`).** No mathematical error found; fusion chain exact for
> k ≤ 8; theorem checked against string functions computed from Weyl–Kac for k ≤ 5. The result is
> **conditional on two cited theorems not re-derived here: Lepowsky–Primc (matrix 2C(A_{k−1})⁻¹, constraint
> Σ j·m_j ≡ Q mod k) and Kac–Peterson (modularity of string functions).** Normalisation used: Kac–Peterson,
> c⁰_{2N} = q^{−c/24−N²/k} Σ_n mult(2N; n) q^n, c = 3k/(k+2); equivalently
> F_{2k} = Σ_N q^{−2N²/(2k+1)}·(weight-2N multiplicity series of the level-k vacuum ŝl_2-module), i.e. sector by
> sector (T1,A_{2k})[N] = q^{−2N²/(2k+1)}·(A_2,T_k)[N] with (A_2,T_k) the Feigin–Stoyanovsky vacuum character.
> Context to add: this is the A_{2k} member of the "♭-case" of Cherednik–Feigin arXiv:1209.1978 §3.6
> (c_eff = n/(n+4) in their normalisation), which they could not reach. "An extension of Z_k ⊗ U(1)" is loose:
> F_{2k} is a Z-linear combination of product-theory characters. The KKMM citation should be to their §2,
> (2.1)–(2.3). Chen–Yin 2610.00967 does not cite Sun–Wang. Not checked by anyone here: Dunn's thesis
> (*Nahm's conjecture in the Cartan case*), the likeliest place for a prior statement.


Written 2026-11-26. Certificate and numerical checks: `code/sw_certify_T1.py`, part (e), log
`certificates/sw_certify_T1.txt`.

Source of the problem. K. Sun and H. Wang, arXiv:2604.00847, Table 1: the pair (T1,A_r) has central charge
3r/(r+4), CFT "?", reference "−"; Conjecture 1.1 there predicts that q^{−c/24}Σ_n q^{½nᵀC(A_r)⁻¹n}/∏(q)_{n_i} is modular.
Known cases: r = 1 (Rogers–Ramanujan type, M(5,3)), r = 2 (Zagier, U(1)_3), r = 3 (= (T1,D_3); Zagier's Example 9,
proved by L. Wang; SM(28,2)). **This file treats all even r = 2k.** The odd case r ≥ 5 is not solved (§5).

## 0. Statement

Notation as in `E8-E6.md`. For k ≥ 1 and m = (m_1,…,m_{k−1}) ∈ N^{k−1} put
M_j = m_j + … + m_{k−1}, L = L(m) = Σ_j j·m_j = M_1 + … + M_{k−1}; then mᵀC(A_{k−1})⁻¹m = M_1²+…+M_{k−1}² − L²/k
(because C(A_{k−1})⁻¹_{ij} = min(i,j) − ij/k). For Q ∈ Z/kZ let

    P_{k,Q}(q) = Σ_{m∈N^{k−1}, L(m) ≡ Q (mod k)}  q^{mᵀC(A_{k−1})⁻¹m} / ∏_j (q)_{m_j}          (P_{1,0} = 1),
    ϑ_{k,Q}(q) = Σ_{N∈Z, N ≡ Q (mod k)}  q^{N²/(k(2k+1))}.

**Theorem A.** For every k ≥ 1,

    F_{2k}(q) := Σ_{n∈N^{2k}} q^{½ nᵀ C(A_{2k})⁻¹ n} / ∏_{i=1}^{2k} (q)_{n_i}
       = (1/(q;q)_∞) Σ_{t∈Z} Σ_{m∈N^{k−1}} q^{mᵀC(A_{k−1})⁻¹m + (L(m)+kt)²/(k(2k+1))} / ∏_j (q)_{m_j}                    (0.1)
       = (1/(q;q)_∞) Σ_{Q mod k} P_{k,Q}(q) ϑ_{k,Q}(q).                                                              (0.2)

(0.1), (0.2) are proved below with the Durfee rectangle identity only. Now P_{k,Q} is the Lepowsky–Primc
fermionic form of the Z_k-parafermion characters:

**(LP)** (J. Lepowsky, M. Primc, *Structure of the standard modules for the affine Lie algebra A_1^{(1)}*,
Contemp. Math. 46 (1985); in this form e.g. Kedem–Klassen–McCoy–Melzer, hep-th/9211102, §3, "the Z_N-parafermion
characters b^0_{2Q}"). For k ≥ 2 and all Q: P_{k,Q}(q) = q^{c_k/24} η(τ) c^0_{2Q}(τ), where c_k = 2(k−1)/(k+2) and
c^ℓ_m are the string functions of A_1^{(1)} at level k (Kac–Peterson, Adv. Math. 53 (1984)).
[I did not have the Lepowsky–Primc paper at hand; (LP) was checked for 2 ≤ k ≤ 7 and every Q against the
Kac–Peterson Hecke-type formula for c^0_{2Q}, mod q^{60} (k ≤ 5), q^{30} (k = 6, 7).]

**Corollary A** (modularity; explicit CFT). With c = 3k/(k+2) = 3r/(r+4) (r = 2k),

    q^{−c/24} F_{2k}(q) = Σ_{Q mod k} c^0_{2Q}(τ) · ϑ_{k,Q}(τ) = Σ_{N∈Z} c^0_{2N}(τ) q^{N²/(k(2k+1))}.                    (0.3)

Each c^0_{2Q} is a modular form of weight −½ and each ϑ_{k,Q} a theta series of weight ½ (both on congruence
subgroups, Kac–Peterson), so q^{−c/24}F_{2k} is a modular function: Sun–Wang's Conjecture 1.1 holds for (T1,A_{2k}).
In CFT terms (η c^0_{2Q} = Z_k parafermion character of weight Q(k−Q)/k, ϑ/η = characters of a free boson):

    q^{−c/24} F_{2k} = Σ_{s mod k(2k+1)} χ^{Z_k}_{0,2s}(τ) · χ^{U(1)_{k(2k+1)}}_{2s}(τ),

where χ^{U(1)_K}_j = η⁻¹Σ_n q^{K(n+j/(2K))²} (weight j²/(4K)); i.e. the "?" of Table 1 is an extension of
(Z_k parafermions) ⊗ U(1)_{k(2k+1)} — the same parafermions as in SU(2)_k = Z_k ⊗ U(1)_k, with another radius of
the boson. Examples: k = 1: F_2 = (q)_∞⁻¹Σ_N q^{N²/3} (Sun–Wang (3.24)); k = 2 (c = 3/2):

    q^{−1/16} F_4 = χ^{Ising}_0 (χ_0 + 2χ_4 + 2χ_8)^{U(1)_{10}} + χ^{Ising}_{1/2} (2χ_2 + 2χ_6 + χ_{10})^{U(1)_{10}}.

## 1. The rank k+1 form

Let G_k be the symmetric matrix of the quadratic form (so that Q_G(v) = ½vᵀG_kv)

    Q_G(a,b,m) = ab + mᵀC(A_{k−1})⁻¹m + (L + k(a−b))²/(k(2k+1))
              = ab + ( k(a−b)² + 2L(a−b) − 2L² )/(2k+1) + M_1²+…+M_{k−1}²

in the k+1 variables (b, m_1, …, m_{k−1}, a). Name these variables after the nodes 1, 2, 4, …, 2k−2, 2k of A_{2k}
(b ↔ 1, m_j ↔ 2j, a ↔ 2k). With n := 2k+1 and C⁻¹ := C(A_{2k})⁻¹, C⁻¹_{xy} = min(x,y) − xy/n, one reads off

    (G_k)_{xy} = C⁻¹_{xy} − [ {x,y} = {1, 2j} for some 1 ≤ j ≤ k−1 ]        (x,y ∈ {1,2,4,…,2k}).                    (1.1)

[Check: (b,b): 2k/n = 1 − 1/n. (a,a): 2k/n = 2k − 4k²/n. (a,b): 1 − 2k/n = 1/n = C⁻¹_{1,2k}. (m_i,m_j):
2min(i,j) − 4ij/n = C⁻¹_{2i,2j}. (a,m_j): 2j/n = 2j − 4jk/n = C⁻¹_{2j,2k}. (b,m_j): −2j/n = C⁻¹_{1,2j} − 1. Exact machine
check for k ≤ 10 in the certificate.]

**Proposition 1.** Starting from G_k, the k−1 fusions (Lemma 2 of `E8-E6.md`)

    node 2j+1 := node 1 + node 2j,        in the order j = k−1, k−2, …, 1,

produce C(A_{2k})⁻¹. Consequently F_{2k} = F_{G_k,0}.

*Proof.* Let A be the matrix before the step with index j: it lives on the nodes X_j = {1,2,4,…,2k} ∪ {2i+1 : i > j},
and we claim A_{xy} = C⁻¹_{xy} − [{x,y} = {1,2i}, i ≤ j] (true for j = k−1 by (1.1)). The fusion adjoins the node
c = 2j+1 with A°_{cx} = A_{1x} + A_{2j,x} for x ∈ X_j, A°_{cc} = A_{11} + A_{2j,2j} + 2A_{1,2j}, and raises A_{1,2j} by 1. Since
−xy/n is additive in x, we have to compare the "min" parts:
x = 1: min(1,1) + (min(2j,1) − 1) = 1 = min(2j+1,1). x = 2j: (1−1) + 2j = min(2j+1,2j). x = 2i, i < j (entry (1,2i)
still lowered): (1−1) + 2i = min(2j+1,2i). x = 2i with j < i ≤ k (entry (1,2i) already raised, or i = k, never
lowered): 1 + 2j = min(2j+1,2i). x = 2i+1, i > j: 1 + 2j = min(2j+1,2i+1). Diagonal: 1 + 2j + 2(1−1) = 2j+1 for
the min part and −(1 + 4j² + 4j)/n = −(2j+1)²/n for the rest. So row c of A° is row 2j+1 of C⁻¹, and after the
raising the claim holds for j−1. After the last step (j = 1) no entry is lowered and A = C(A_{2k})⁻¹ on all 2k
nodes. Admissibility of (G_k,0): Q_G ≥ ab + mᵀC(A_{k−1})⁻¹m and Q_G ≥ (L+k(a−b))²/(k(2k+1)), so Q_G ≤ E bounds m, ab
and |a−b|, hence a and b. Lemma 2 of `E8-E6.md` (applied k−1 times) gives F_{C(A_{2k})⁻¹,0} = F_{G_k,0}. ∎

## 2. Proof of Theorem A

By Proposition 1,

    F_{2k} = Σ_{m∈N^{k−1}} ( q^{mᵀC(A_{k−1})⁻¹m} / ∏(q)_{m_j} ) Σ_{t∈Z} q^{(L+kt)²/(k(2k+1))} Σ_{a,b≥0, a−b=t} q^{ab}/((q)_a(q)_b)

(regrouping of a summable family). The inner sum is Cauchy's identity Σ_{b≥0} q^{b(b+|t|)}/((q)_b(q)_{b+|t|}) = 1/(q)_∞
(proved from the Durfee rectangle identity in `E7.md`, proof of Lemma 4; Andrews, *The Theory of Partitions*,
(2.2.8)), for every t ∈ Z. This is (0.1). For fixed m with L(m) ≡ Q (mod k), N = L + kt runs over all integers
≡ Q (mod k) exactly once, so Σ_t q^{(L+kt)²/(k(2k+1))} = ϑ_{k,Q}; summing over m gives (0.2). ∎

*Proof of Corollary A.* Insert (LP) in (0.2): F_{2k} = q^{c_k/24} (η(τ)/(q)_∞) Σ_Q c^0_{2Q} ϑ_{k,Q} = q^{(c_k+1)/24} Σ_Q c^0_{2Q}ϑ_{k,Q},
and c_k + 1 = 3k/(k+2). The second form of (0.3) uses c^0_{2N} = c^0_{2Q} for N ≡ Q (mod k). For the CFT form write
N = s + k(2k+1)n′, s mod k(2k+1): q^{N²/K} = q^{K(n′+s/K)²}, K = k(2k+1), and η c^0_{2s} = χ^{Z_k}_{0,2s}. For k = 2:
ηc^0_0 = χ^{Ising}_0, ηc^0_2 = χ^{Ising}_{1/2}, and s runs mod 10 with s even ↔ Q = 0. ∎

## 3. Verification

Certificate part (e): (1.1) and Proposition 1 exactly for k ≤ 10; F_{2k} = F_{G_k,0} as power series (k ≤ 3 mod q^{60},
k = 4 mod q^{30}, k = 5 mod q^{20}); (0.3) with the Kac–Peterson formula for the string functions, same orders;
(LP) for k ≤ 7, sector by sector. (An earlier exploratory run, `sw_audit`-independent, had confirmed (0.3) for
k = 4 mod q^{24} before the proof was found.)

## 4. Novelty

I found no statement of (0.1)–(0.3) (for k ≥ 2) in: Sun–Wang (CFT "?"), the papers citing it that I could read
(Wang–Wang 2607.08606, Yoshida 2608.08780, Mizuno 2609.09653, Xia 2609.08816, Guo 2610.15587, Chen–Yin
2610.00967), the Wang-school papers listed in `PRIOR-ART.md`, Kedem–Klassen–McCoy–Melzer and Dasmahapatra et al.
(which treat 2C⁻¹-type forms). arXiv listing searches ("Nahm sum" with "parafermion", "string function",
"inverse Cartan") returned nothing relevant. The proof is a (k−1)-fold application of the known fusion lemma to
the known Lepowsky–Primc sums, so it is certainly within reach of the experts; the 1990s physics literature on
fermionic characters (e.g. of SU(2)_k and its orbifolds/"lens space" cosets, where products of the form
(parafermion) × (boson at another radius) are standard) may well contain equivalent formulas that the keyword
searches available to me cannot find. No full-text search engine was available.

## 5. What is missing: odd r

For r = 2k−1 the same fusions (node 2j+1 := node 1 + node 2j, j = k−1, …, 1, on the nodes 1, 2, 4, …, 2k−2) reduce
C(A_{2k−1})⁻¹ to the rank k form (exact check for k ≤ 8 in the certificate; same computation as in Proposition 1)

    F_{2k−1} = Σ_{b≥0} Σ_{m∈N^{k−1}} q^{b²/2 + M_1²+…+M_{k−1}² − (L + b/2)²/k} / ( (q)_b ∏(q)_{m_j} ),

i.e. the Lepowsky–Primc sum coupled to one "fermionic" variable b carrying half a unit of Z_k charge (k = 1:
Σ q^{b²/4}/(q)_b; k = 2: Zagier's Example 9). There is no Cauchy identity to remove b, and I have no closed formula
for r ≥ 5. Numerical observation only (r = 5, mod q^{60}, exploratory script, not archived): the four sectors of
F_5 (exponents ≡ 0, 5/12, 2/3, 3/4 mod 1) have rapidly growing coefficients, but after multiplication by (q)_∞² all
coefficients up to q^{60} have absolute value ≤ 4 and between a quarter and a half of them vanish — the behaviour of a weight-one (Hecke-type)
double theta series, as for string functions, not that of a product of two unary theta series. This fits the
heuristic "Z_κ parafermions at level κ = k − ½" (central charge 3κ/(κ+2)), but I did not identify the forms.
