> **Independent referee (coordinating session, `../REFEREE-2-GKPS.md`, 26 Nov): proofs complete for
> all 17 examples treated; statements match arXiv:2411.06081. Corrections to this file: 17 examples are
> treated, not 18 (16 in all versions; 3-2 v3 not treated); three claims of the paper are missing (3-2 v2,
> 3-4 v2 and v3 with Wilson lines; true, one line from the existing arguments); PRIOR PROOFS exist for 5
> examples: 3-7 v3 and 3-14 v1 (Zagier 2007, Example 2), 3-9, 3-13 (with Wilson line), 3-18 (L. Wang,
> Adv. Math. 2024, arXiv:2211.04375, Thms 6, 10, 9); for the untreated 3-15, 3-21, 3-22 the underlying
> q-series identities are in Wang and Cao–Wang arXiv:2412.15767. The statement in section 4 that no earlier
> proof was found is therefore false for those 5. No earlier proof was found for the genuinely signed
> cases (search not exhaustive).**

# The rank-three signed Nahm sum identities of Gang–Kim–Park–Stubbs: elementary proofs

Written 2026-11-26 (third session). Certificate (series checks): `code/gkps_certify.py` → `certificates/gkps_certify.txt`.

**Source.** D. Gang, H. Kim, B. Park, S. Stubbs, *Three dimensional topological field theories and Nahm sum
formulas*, arXiv:2411.06081, §4.1.3. For a positive definite integer 3×3 matrix K and v = μ_0 − a, Q ∈ Z³ put

    I(K, v, Q) := Σ_{m∈N³} q^{½ mᵀKm} (−q^{1/2})^{−v·m} q^{−Q·m} / ( (q)_{m_1}(q)_{m_2}(q)_{m_3} )            (their I_A[W_Q]; see §4 for the verbatim formula).

In the source (§4.1, first paragraph): "The identification of the indices with known RCFT characters are
conjectures based on the q-series expansion, except for a few cases with known Nahm sum representations."
The only rank-three case for which the source gives a proof reference is I_A[1] of example 3-13.

**Result.** All the identities of the list in §2 hold: 18 of the 24 rank-three examples of the source, with all their
versions and Wilson lines (87 sums in the certificate, counting the family 3-12 for |a| ≤ 4). One of them (3-18) is
the classical Andrews–Gordon identity and one (I_A[1] of 3-13) is attributed by the source to Li–Wang; for the
others the source gives no proof. The remaining examples are listed in §4. Every one of them follows from the q-binomial theorem
(L1) or the q-Chu–Vandermonde sum (L2) or Euler's identity (E) — which "integrate out" the summation variable
that carries the sign — followed by the Rogers–Ramanujan or Andrews–Gordon identities. In most cases the sum
collapses to a one-fold sum: the signed variable forces another variable to vanish.

Notation: G(q) = Σ_n q^{n²}/(q)_n = 1/(q,q⁴;q⁵)_∞, H(q) = Σ_n q^{n²+n}/(q)_n = 1/(q²,q³;q⁵)_∞ (Rogers–Ramanujan), and for i = 1, 2, 3

    AG_i(q) := Σ_{N_1≥N_2≥0} q^{N_1²+N_2²+N_i+…+N_2} / ( (q)_{N_1−N_2}(q)_{N_2} ) = ∏_{n≢0,±i (7)} (1−q^n)^{−1}        (Andrews–Gordon, k = 3),

i.e. AG_1 = 1/(q²,q³,q⁴,q⁵;q⁷)_∞ (linear term N_1+N_2), AG_2 = 1/(q,q³,q⁴,q⁶;q⁷)_∞ (linear term N_2), AG_3 = 1/(q,q²,q⁵,q⁶;q⁷)_∞ (no linear term).
All sums are over finitely many terms for each power of q^{1/2} (K is positive definite), so every rearrangement
below is legitimate in Z((q^{1/2})).

## 1. Lemmas

**(L1) q-binomial theorem.** Σ_{j=0}^{M} (−1)^j q^{j(j−1)/2} z^j / ((q)_j(q)_{M−j}) = (z;q)_M/(q)_M. In particular the sum is δ_{M,0} for
z = 1, and it is 0 for z = q^{−c} with an integer 0 ≤ c < M.

**(L2) q-Chu–Vandermonde.** For integers N, M ≥ 0: Σ_{j≥0} (−1)^j q^{j(j−1)/2} / ((q)_j(q)_{N−j}(q)_{M−j}) = q^{NM}/((q)_N(q)_M).
*Proof.* Multiply by x^M and sum over M: Σ_M x^M/(q)_{M−j} = x^j/(x)_∞, so the left side becomes
(x)_∞^{−1}Σ_j(−1)^j q^{j(j−1)/2}x^j/((q)_j(q)_{N−j}) = (x;q)_N/((q)_N(x)_∞) = 1/((q)_N(xq^N)_∞) = (q)_N^{−1}Σ_M x^Mq^{NM}/(q)_M by (L1) and Euler. ∎

**(E) Euler.** Σ_{n≥0} (−1)^n q^{n(n−1)/2} z^n/(q)_n = (z;q)_∞; for z = q^{c}, c ∈ Z: this is 0 if c ≤ 0 and (q)_∞/(q)_{c−1} if c ≥ 1.

**(C) Cauchy/Durfee.** Σ_{m≥0} q^{m(m+c)}/((q)_m(q)_{m+c}) = 1/(q)_∞ for every integer c ≥ 0.

**(L4) Finite pentagonal number theorem** (Rogers; e.g. Andrews, *The Theory of Partitions*, Example 3.10 type; checked
by the certificate): Σ_{i=−n}^{n} (−1)^i q^{i(3i−1)/2}/((q)_{n+i}(q)_{n−i}) = 1/(q)_n, and the same with i(3i+1)/2 (i ↦ −i).

**(L5) Durfee rectangle** (`E7.md`, Lemma 1): Σ_c q^{(A−c)(B−c)}/((q)_c(q)_{A−c}(q)_{B−c}) = 1/((q)_A(q)_B); and its limit (C) in the form
Σ_{m_2,m_3≥0, m_3−m_2=d} q^{m_2m_3}/((q)_{m_2}(q)_{m_3}) = 1/(q)_∞ for every d ∈ Z.

**(L3)** (`E7.md`, Lemma 4). For integers a, R ≥ 0: Σ_{m≥0} q^{m²+Rm}/((q)_m(q)_{m+R+a}) = (q)_∞^{−1} Σ_{j=0}^{a} [a;j]_q q^{j²+Rj}.

## 2. The identities and their proofs

Each line: label in the source, K, v = μ_0 − a, the Wilson-line charge Q, the value. φ := (−q^{1/2};q)_∞ = q^{1/48}χ_F;
ϑ_ε := (q)_∞^{−1}Σ_{u∈Z}q^{u²+εu} (ε = 0, 1; q^{−1/24}ϑ_0 = χ^{U(1)_2}_0, q^{5/24}ϑ_1 = χ^{U(1)_2}_1).
(In the source the values are written as q^{P}·(characters), e.g. q^{−11/60}χ^{M(2,5)}_{(1,1)} = H, q^{1/60}χ^{M(2,5)}_{(1,2)} = G, q^{−17/42}χ^{M(2,7)}_{(1,1)} = AG_1,
q^{−5/42}χ^{M(2,7)}_{(1,2)} = AG_2, q^{1/42}χ^{M(2,7)}_{(1,3)} = AG_3, q^{−13/80}χ_Fχ^{M(2,5)}_{(1,1)} = φH, q^{3/80}χ_Fχ^{M(2,5)}_{(1,2)} = φG, q^{−17/120}χ_F²χ^{M(2,5)}_{(1,1)} = φ²H,
q^{7/120}χ_F²χ^{M(2,5)}_{(1,2)} = φ²G; the prefactors were checked.)

| label | K | v = μ_0 − a | Q | I(K,v,Q) |
|---|---|---|---|---|
| 3-2, case 1 | [[2,−1,−1],[−1,2,1],[−1,1,1]] | (−2,2,1) | 0 ; (1,−1,−1) | H ; G |
| 3-2, case 2 | same | (−2,1,0) | 0 | H |
| 3-4, case 1 | [[2,1,1],[1,2,1],[1,1,1]] | (−2,−2,−3) | 0 ; (1,1,1) | H ; G |
| 3-4, cases 2, 3 | same | (−2,−1,−2) ; (−1,−2,−2) | 0 | H |
| 3-3 | [[2,1,−1],[1,2,0],[−1,0,1]] | (−2,−1,0) | 0 ; (1,1,0) | φH ; φG |
| 3-6 | [[3,−1,−1],[−1,1,0],[−1,0,1]] | (−5,1,1) ; (−3,0,1) ; (−3,1,0) | 0 ; (1,0,0) | 0 ; 0 |
| 3-7, versions 1, 2 | [[3,−1,1],[−1,1,0],[1,0,1]] | (−3,1,−1) ; (−5,1,−2) | 0 ; (1,0,0) | 0 ; 0 |
| 3-7, version 3 | same | (−2,0,0) | 0 ; (1,0,0) | φ²H ; φ²G |
| 3-8 | [[3,2,−1],[2,2,−1],[−1,−1,1]] | (−3,−3,1) | 0 ; (1,1,0) | 0 ; 0 |
| 3-9 | [[3,2,1],[2,2,1],[1,1,1]] | (−2,−2,0) | 0 ; (1,1,0) | φH ; φG |
| 3-10 | [[5,4,−2],[4,4,−2],[−2,−2,2]] | (−1,−2,0) ; (−2,−3,0) | 0 ; (2,2,−1) | G ; H |
| 3-23 | [[5,−3,2],[−3,3,−1],[2,−1,1]] | (−5,3,−3) | 0 ; (2,−1,1) ; (3,−2,1) ; (4,−3,2) | 0 ; 1 ; 1 ; 0 |
| 3-5 | [[3,−2,2],[−2,2,−1],[2,−1,2]] | (0,0,−1) | 0 ; (1,−1,1) | φG ; φH |
| 3-13 | [[4,2,2],[2,2,1],[2,1,2]] | (−4,−2,−2) | 0 ; (2,1,1) | H² ; G² |
| 3-14, version 1 | [[4,−2,2],[−2,2,−1],[2,−1,2]] | (−4,2,−2) | 0 ; (1,−1,1) ; (2,−1,1) ; (3,−2,2) | ϑ_1H ; ϑ_0H ; ϑ_0G ; ϑ_1G |
| 3-14, version 2 | same | (−5,3,−3) | 0 ; (1,−1,1) ; (2,−1,1) ; (3,−2,2) | −q^{−1/2}H/φ ; H/φ ; G/φ ; G/φ |
| 3-18 | [[6,4,2],[4,4,2],[2,2,2]] | (−6,−4,−2) | 0 ; (1,1,1) ; (2,2,1) ; (3,2,1) | ∏_{n≢0,±i (9)}(1−q^n)^{−1}, i = 1, 2, 3, 4 |
| 3-12 (family, a ∈ Z∖{0}; 3-11 is a = 2), case 1 | K_a = [[a²+1,a²,a],[a²,a²,a],[a,a,2]] | (−1−a,−2−a,−2) for a > 0, (1−a,−a,−2) for a < 0 | 0 ; (a,a,1) | H ; G |
| 3-12, case 2 | K_a | (−a,−1−a,−2) | 0 ; (a,a,1) | H ; G |
| 3-16 | [[4,−2,−1],[−2,2,1],[−1,1,1]] | (−4,2,1) | 0 ; (1,0,0) ; (2,−1,−1) | AG_1 ; AG_2 ; AG_3 |
| 3-17 | [[9,5,3],[5,4,1],[3,1,2]] | (−5,−4,−2) | 0 ; (2,1,1) ; (3,2,1) | AG_1 ; AG_2 ; AG_3 |

**3-4** (K = [[2,1,1],[1,2,1],[1,1,1]]; ½mᵀKm = m_1²+m_2²+m_1m_2 + ½m_3² + m_3(m_1+m_2)).
*Case 1.* (−q^{1/2})^{2m_1+2m_2+3m_3} = (−1)^{m_3}q^{m_1+m_2+3m_3/2}, so the m_3-sum is Σ(−1)^{m_3}q^{m_3(m_3−1)/2}(q^{m_1+m_2+2−ε})^{m_3}/(q)_{m_3}
with ε = 0 for Q = 0 and ε = 1 for Q = (1,1,1); by (E) it equals (q)_∞/(q)_{m_1+m_2+1−ε}. Hence

    I = (q)_∞ Σ_{m_1,m_2} q^{m_1(m_1+m_2+1−ε) + m_2²+(1−ε)m_2} / ( (q)_{m_1}(q)_{m_2}(q)_{m_1+m_2+1−ε} ) = Σ_{m_2} q^{m_2²+(1−ε)m_2}/(q)_{m_2}

by (C) with c = m_2+1−ε: this is H for ε = 0 and G for ε = 1.
*Case 2.* Put M := m_2+m_3. Then ½mᵀKm = m_1² + m_1M + ½M² + ½m_2², and (−q^{1/2})^{2m_1+m_2+2m_3} = (−1)^{m_2}q^{m_1+M−m_2/2}. The
m_2-sum is Σ_{m_2}(−1)^{m_2}q^{m_2(m_2−1)/2}/((q)_{m_2}(q)_{M−m_2}) = δ_{M,0} by (L1). So I = Σ_{m_1}q^{m_1²+m_1}/(q)_{m_1} = H. *Case 3* is case 2
with m_1 and m_2 interchanged (K is symmetric under this).

**3-2** (K = [[2,−1,−1],[−1,2,1],[−1,1,1]]). *Case 1.* (−q^{1/2})^{2m_1−2m_2−m_3} = (−1)^{m_3}q^{m_1−m_2−m_3/2}; the m_3-sum is (E) with z = q^{m_2−m_1+ε}
(ε = 0 for Q = 0, ε = 1 for Q = (1,−1,−1), which contributes q^{−m_1+m_2+m_3}), so only m_2 = m_1+1−ε+j, j ≥ 0, survive, with value
(q)_∞/(q)_j. The remaining exponent m_1²+m_2²−m_1m_2 + (1−ε)(m_1−m_2) equals m_1(m_1+j+1−ε) + j²+(1−ε)j, and

    I = (q)_∞ Σ_{m_1,j} q^{m_1(m_1+j+1−ε) + j²+(1−ε)j} / ( (q)_{m_1}(q)_{m_1+j+1−ε}(q)_j ) = Σ_j q^{j²+(1−ε)j}/(q)_j        by (C):

H for Q = 0, G for Q = (1,−1,−1).
*Case 2.* With M := m_2+m_3: ½mᵀKm = m_1² − m_1M + ½M² + ½m_2², and (−q^{1/2})^{2m_1−m_2} = (−1)^{m_2}q^{m_1−m_2/2}; by (L1) (z = 1) only M = 0
survives: I = Σ_{m_1}q^{m_1²+m_1}/(q)_{m_1} = H.

**3-10, 3-11, 3-12** (K_a = [[a²+1,a²,a],[a²,a²,a],[a,a,2]], a ∈ Z∖{0}; 3-10 is a = −2 with a different pair (μ_0, a-vector), 3-11
is a = 2). With N := m_1+m_2: **½mᵀK_am = ½(aN+m_3)² + ½m_1² + ½m_3²** = ½(aN+m_3)² + ½m_3² + ½N² − Nm_2 + ½m_2².
*Case 1, and 3-10.* In all three situations v·m ≡ aN + m_1 (mod 2) and −½v·m = (linear form in N, m_3) − ½m_1
[a > 0: −½v·m = ½((2+a)N − m_1) + m_3; a < 0: ½(aN − m_1) + m_3; 3-10: N − ½m_1]. So the m_1-sum at fixed N, m_3 is
Σ_{m_1}(−1)^{m_1}q^{m_1(m_1−1)/2}/((q)_{m_1}(q)_{N−m_1}) = δ_{N,0}, and I = Σ_{m_3} q^{m_3² + ℓ m_3}/(q)_{m_3}, where ℓ is the coefficient of m_3 in −½v·m − Q·m:
case 1: ℓ = 1 for Q = 0 (I = H), ℓ = 0 for Q = (a,a,1) (I = G); 3-10: ℓ = 0 for Q = 0 (I = G), ℓ = 1 for Q = (2,2,−1) (I = H).
*Case 2.* v·m ≡ aN + m_2 and −½v·m = ½aN + ½m_2 + m_3. In the variables (N, m_2, m_3) the m_2-sum is
Σ_{m_2}(−1)^{m_2}q^{m_2(m_2−1)/2}(q^{1−N})^{m_2}/((q)_{m_2}(q)_{N−m_2}) = (q^{1−N};q)_N/(q)_N = δ_{N,0} by (L1). Again I = Σ_{m_3}q^{m_3²+ℓm_3}/(q)_{m_3} with ℓ = 1 (Q = 0, I = H)
and ℓ = 0 (Q = (a,a,1), I = G).

*3-10, second version* (v = (−2,−3,0)): v·m ≡ m_2 and −½v·m = N + ½m_2, so as in case 2 the m_2-sum is (q^{1−N};q)_N/(q)_N = δ_{N,0}.

**3-3** (K = [[2,1,−1],[1,2,0],[−1,0,1]], v = (−2,−1,0): sign (−1)^{m_2}, linear term m_1+½m_2). The m_3-sum is
Σ_{m_3}q^{m_3(m_3−1)/2}(q^{1/2−m_1})^{m_3}/(q)_{m_3} = (−q^{1/2−m_1};q)_∞ = φ·q^{−m_1²/2}(−q^{1/2};q)_{m_1} (Euler), and (−q^{1/2};q)_{m_1} = Σ_{k+s=m_1}(q)_{m_1}q^{k²/2}/((q)_k(q)_s) (L1). Hence

    I = φ Σ_{k,s,m_2≥0} (−1)^{m_2} q^{k²+ks+½s²+m_2²+(k+s)m_2 + k+s+½m_2 − Q_1(k+s) − Q_2m_2} / ((q)_k(q)_s(q)_{m_2}),

which is φ times the sum of example 3-4, version 2, in the variables (m_1,m_2,m_3) = (k,m_2,s) (for Q = 0), resp. the same
sum with linear term −½m_2 (for Q = (1,1,0)); by the computation of 3-4 (M = m_2+s, (L1) with z = 1) they are
Σ_kq^{k²+k}/(q)_k = H and Σ_kq^{k²}/(q)_k = G. So I = φH, φG.

**3-9** (K = [[3,2,1],[2,2,1],[1,1,1]], v = (−2,−2,0): no sign). ½mᵀKm = ½(N+m_3)² + ½N² + ½m_1² with N = m_1+m_2, linear term N − Q·m.
Euler: Σ_{m_3}q^{(N+m_3)²/2}/(q)_{m_3} = q^{N²/2}(−q^{N+1/2};q)_∞ = q^{N²/2}φ/(−q^{1/2};q)_N; (L1): Σ_{m_1+m_2=N}q^{m_1²/2}/((q)_{m_1}(q)_{m_2}) = (−q^{1/2};q)_N/(q)_N. So
I = φΣ_Nq^{N²+N−Q_1N}/(q)_N = φH (Q = 0), φG (Q = (1,1,0)). [This is the pair (T_1,T_3) with a linear term; cf. Warnaar.]

**3-7** (K = [[3,−1,1],[−1,1,0],[1,0,1]]). *Versions 1, 2* (v = (−3,1,−1), (−5,1,−2)): the m_2-dependence is (−1)^{m_2}q^{m_2(m_2−1)/2}(q^{−m_1})^{m_2}/(q)_{m_2},
whose sum is (q^{−m_1};q)_∞ = 0 for every m_1 ≥ 0 (E). So I = 0, for every Q with Q_2 = 0.
*Version 3* (v = (−2,0,0): no sign, linear term m_1): by Euler Σ_{m_2}q^{m_2²/2−m_1m_2}/(q)_{m_2} = φ q^{−m_1²/2}(−q^{1/2};q)_{m_1} and
Σ_{m_3}q^{m_3²/2+m_1m_3}/(q)_{m_3} = φ/(−q^{1/2};q)_{m_1}; so I = φ²Σ_{m_1}q^{m_1²+m_1−Q_1m_1}/(q)_{m_1} = φ²H, φ²G.

**3-6** (K = [[3,−1,−1],[−1,1,0],[−1,0,1]]): for v = (−5,1,1) and (−3,1,0) the m_2-sum is (q^{−m_1};q)_∞ = 0 as in 3-7; for v = (−3,0,1) the
m_3-sum is. **3-8** (K = [[3,2,−1],[2,2,−1],[−1,−1,1]], v = (−3,−3,1)): the m_3-sum is (q^{−m_1−m_2};q)_∞ = 0. In both examples I = 0 for the
Wilson lines of the table (they do not involve the variable that is summed).

**3-23** (K = [[5,−3,2],[−3,3,−1],[2,−1,1]], v = (−5,3,−3); sign (−1)^{m_1+m_2+m_3}). The m_3-sum is (E) with z = q^{2m_1−m_2+2−Q_3}: it is
(q)_∞/(q)_{M−m_2} for m_2 ≤ M := 2m_1+1−Q_3 and 0 otherwise. What remains depends on m_2 through
(−1)^{m_2}q^{(3/2)m_2(m_2−M) + κm_2}/((q)_{m_2}(q)_{M−m_2}) with κ = 0, −½, +½, 0 for Q = 0, (2,−1,1), (3,−2,1), (4,−3,2).
*Q = 0 and Q = (4,−3,2)*: M is odd (2m_1+1, 2m_1−1) and κ = 0, so the summand is antisymmetric under m_2 ↦ M−m_2: I = 0.
*Q = (2,−1,1), (3,−2,1)*: M = 2m_1 =: 2n; with m_2 = n+i the m_2-sum is (−1)^n q^{−(3/2)n² ∓ n/2} Σ_i(−1)^iq^{i(3i∓1)/2}/((q)_{n+i}(q)_{n−i}) = (−1)^nq^{−(3/2)n²∓n/2}/(q)_n by (L4).
The remaining factor is (−1)^{n}q^{(5/2)n² ± n/2}(q)_∞/(q)_n, so I = (q)_∞Σ_nq^{n²}/(q)_n² = 1 (Durfee square, (C) with c = 0).

**3-5** (K = [[3,−2,2],[−2,2,−1],[2,−1,2]], v = (0,0,−1): sign (−1)^{m_3}, linear term ½m_3). With N := m_1+m_3:
½mᵀKm = (3/2)N² + m_2² − 2Nm_2 + ½m_3² − Nm_3 + m_2m_3, so the m_3-sum is Σ(−1)^{m_3}q^{m_3(m_3−1)/2}(q^{1+m_2−N})^{m_3}/((q)_{m_3}(q)_{N−m_3}) = (q^{1+m_2−N};q)_N/(q)_N by (L1); it vanishes
unless m_2 = N+j, j ≥ 0, where it is (q)_{N+j}/((q)_j(q)_N²)·(q)_N = (q)_{m_2}/((q)_j(q)_N). The remaining exponent is (3/2)N²+m_2²−2Nm_2 = ½N²+j², and
Q = (1,−1,1) contributes q^{−N+m_2} = q^{j}. Hence I = Σ_{N,j}q^{½N²+j²+εj}/((q)_N(q)_j) = φG (Q = 0), φH (Q = (1,−1,1)).

**3-13** (K = [[4,2,2],[2,2,1],[2,1,2]], v = (−4,−2,−2): no sign). With N = m_1+m_2, M = m_1+m_3: ½mᵀKm = N²+M² + (N−m_1)(M−m_1), and
−½v·m = N+M. By (L5) the m_1-sum is 1/((q)_N(q)_M), so I = Σ_{N,M}q^{N²+M²+N+M−Q·m}/((q)_N(q)_M) with Q·m = 0 resp. N+M: I = H², G².

**3-14** (K = [[4,−2,2],[−2,2,−1],[2,−1,2]]). ½mᵀKm = m_1² + u² + m_2m_3 with u := m_1+m_3−m_2 ∈ Z. For fixed m_1 and u the difference
d = m_3−m_2 = u−m_1 is fixed, and by (L5) Σ_{m_3−m_2=d}(±1)^{m_2+m_3}q^{m_2m_3}/((q)_{m_2}(q)_{m_3}) = (±1)^d/(q)_∞.
*Version 1* (no sign; −½v·m = 2m_1−m_2+m_3 = m_1+u): I = (q)_∞^{−1}Σ_{u∈Z}q^{u²+u−αu}·Σ_{m_1}q^{m_1²+m_1−βm_1}/(q)_{m_1}, where Q·m = αu + βm_1 with
(α,β) = (0,0), (1,0), (1,1), (2,1) for the four Wilson lines: ϑ_1H, ϑ_0H, ϑ_0G, ϑ_1G (Σ_uq^{u²−u} = Σ_uq^{u²+u}).
*Version 2* (sign (−1)^{m_1+m_2+m_3} = (−1)^{m_1+d}·(−1)^{2m_2}, i.e. (−1)^u after the (m_2,m_3)-sum; −½v·m = m_1 + (3/2)u):
I = (q)_∞^{−1}Σ_u(−1)^uq^{u²+(3/2−α)u}·Σ_{m_1}q^{m_1²+m_1−βm_1}/(q)_{m_1}. By the Jacobi triple product Σ_u(−1)^uq^{u²+u/2} = (q^{1/2};q)_∞(q²;q²)_∞ = (q)_∞/φ
(Euler: (q^{1/2};q)_∞(−q^{1/2};q)_∞(−q;q)_∞ = 1), Σ_u(−1)^uq^{u²−u/2} is the same, and Σ_u(−1)^uq^{u²+3u/2} = −q^{−1/2}(q)_∞/φ (shift u ↦ u−1). Hence
I = −q^{−1/2}H/φ, H/φ, G/φ, G/φ. (The source flags version 2 as physically suspicious; as q-series identities the claims hold.)

**3-18** (K = 2(min(4−i,4−j))_{i,j}, v = (−6,−4,−2): no sign): with N_3 = m_1, N_2 = m_1+m_2, N_1 = m_1+m_2+m_3 the sum is
Σq^{N_1²+N_2²+N_3²+N_i+…+N_3}/((q)_{N_1−N_2}(q)_{N_2−N_3}(q)_{N_3}), the Andrews–Gordon identities for k = 4 (classical; listed for completeness).

**3-16** (K = [[4,−2,−1],[−2,2,1],[−1,1,1]], v = (−4,2,1)). (−q^{1/2})^{4m_1−2m_2−m_3} = (−1)^{m_3}q^{2m_1−m_2−m_3/2}; the m_3-sum is (E) with z = q^{m_2−m_1+ε},
ε the coefficient of m_3 in −Q·m (ε = 0, 0, 1 for Q = 0, (1,0,0), (2,−1,−1)). So m_2 = m_1+1−ε+j, the sum equals (q)_∞/(q)_j, and the
remaining exponent 2m_1²+m_2²−2m_1m_2+2m_1−m_2 − Q_1m_1 − Q_2m_2 is

    Q = 0: m_1²+m_1+j²+j;      Q = (1,0,0): m_1²+j²+j;      Q = (2,−1,−1): m_1²+j².

In each case apply (L3) to the sum over j (R = 1, 1, 0 and a = m_1; the denominators are (q)_j(q)_{j+R+m_1}):

    I = (q)_∞ Σ_{m_1,j} q^{m_1²+δm_1+j²+Rj}/((q)_{m_1}(q)_j(q)_{j+R+m_1}) = Σ_{m_1≥i≥0} q^{m_1²+δm_1+i²+Ri}/((q)_i(q)_{m_1−i})        (δ = 1, 0, 0),

which is AG_1, AG_2, AG_3 respectively (N_1 = m_1, N_2 = i).

**3-17** (K = [[9,5,3],[5,4,1],[3,1,2]], v = (−5,−4,−2)). With N := m_1+m_2, M := m_1+m_3: **½mᵀKm = ½m_1² + 2N² + NM + M²**, and
(−q^{1/2})^{5m_1+4m_2+2m_3} = (−1)^{m_1}q^{2N+M−m_1/2}. By (L2) the m_1-sum is q^{NM}/((q)_N(q)_M), so

    I = Σ_{N,M≥0} q^{2N²+2NM+M²+2N+M − Q·m}/((q)_N(q)_M),        Q·m = 0, N+M, 2N+M for Q = 0, (2,1,1), (3,2,1).

With N_1 = N+M, N_2 = N the exponent is N_1²+N_2²+N_1+N_2 − Q·m, i.e. N_1²+N_2²+N_1+N_2, N_1²+N_2²+N_2, N_1²+N_2²: AG_1, AG_2, AG_3. ∎

## 3. Remarks

* The mechanism is uniform: after a unimodular change of variables the signed variable j has quadratic part ½j²
  and enters the other summation indices as M_i − j; then (L1)/(L2) remove it. For 3-2, 3-4 (cases 2, 3), 3-10, 3-11,
  3-12 the result is a Kronecker delta: the three-fold sum *is* a one-fold Rogers–Ramanujan sum. This is the
  q-series shadow of "integrating out" a pair of chiral multiplets in the 3d theory.
* The proofs do not use positivity of K beyond convergence, and give the refined statements with the Wilson lines.

## 4. Source statements and literature

Source text (latest arXiv version; published as SciPost Phys. 19, 128), §3.3, formula labelled `main formula` (I
did not see the compiled equation number): "q^Δ I^{W_Q}_{half}(q, ν=−1, η=1) = Σ_{m∈N^r} q^{½mᵗKm+Δ}(−q^{1/2})^{−(μ_0−a)ᵗm}q^{−Qᵗm}/((q)_{m_1}…(q)_{m_r}), which we
denote by q^Δ I_A[W_Q]"; K is the Chern–Simons level matrix, μ_0 and a integer vectors (R-symmetry mixing and
axial charge), Q the Wilson-line charge. The examples were found by a search over "all positive definite
integer K for r = 1, 2, 3 with entries between −17 and 17". The data (K, μ_0, a, Q, claimed characters with their
q-prefactors) of all rank-three examples were transcribed from the TeX source by a reader agent, which also
verified all 105 sums of the section (ranks 1–3, all versions) to q^25; the transcription is reproduced in
`certificates/gkps_source_transcription.txt`.

Rank-three examples of the source NOT treated in §2: 3-1 and 3-2 version 3 (osp(1|2)_1 characters), 3-15 (osp(1|4)_1;
the characters are not defined in the source), 3-19 – 3-22 (super-Virasoro SM(2,8), SM(2,12)), 3-24 (W_3(3,7)).
(For 3-13 the source says that I_A[1] follows from an identity proved by Li–Wang; the Wilson-line case and the
short proof above via (L5) are not in the source.) All of them hold to q^25
(reader agent's check); I expect the same lemmas to reduce most of them to known rank ≤ 2 identities but have not
done it. Example: for 3-19 (K = [[4,2,1],[2,2,0],[1,0,1]], v = (−3,−2,−2)) one has ½mᵀKm = ½m_1² + N² + ½M² with N = m_1+m_2, M = m_1+m_3 and
(−q^{1/2})^{−v·m} = (−1)^{m_1}q^{N+M−m_1/2}, so by (L2) I = Σ_{N,M}q^{N²+NM+½M²+N+M}/((q)_N(q)_M), which is the rank-two sum of the source's example 2-2
(matrix [[2,1],[1,1]], a Slater-type identity for SM(2,8)).

Later literature (searches by the reader agent: arXiv listings and web search for papers citing 2411.06081,
"signed Nahm sums", Gang–Park–Sohn arXiv:2512.23122 (abstract only), Cao–Wang arXiv:2412.15767 (abstract only)):
no proof of these rank-three identities was found; the citation lists (INSPIRE / Google Scholar) could not be
accessed, so this is not conclusive. Given how short the proofs are, the identities may well be known to the
authors or to experts; what is recorded here is that the source states them as conjectural and that they are
elementary.
