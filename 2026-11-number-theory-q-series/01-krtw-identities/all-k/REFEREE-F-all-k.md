# Referee report: `proofs/F-all-k.md` ((F_k) for all k; KRTW (ii), (iii))

Referee: adversarial read with independent computations. Scripts are in `referee/` (own code; nothing is imported
from `code/` or from `rr-identities`). Date: 2026-11-27. Status: **DONE.**

## Overall verdict

**No fatal error and no substantive gap found.** The proof of (F_k), K_k C_1 = B_{k+1} C_k in Z[x][[q]], is
complete and correct for every k ≥ 1. Its only external inputs are (B1) and (B2) of `KRTW-F-moves.md` Theorem B
(taken as known per the brief; they are used correctly) and standard facts: Bailey's lemma (L∞), Bailey inversion,
Cauchy's identity, and partial fractions. Theorem 2 follows: KRTW (ii) and (iii) hold for all k. I re-derived every
displayed identity by hand, and independent computations confirm each step (list below). The problems found are
exposition only (E1–E6). None affects correctness. Identity (iv) is not claimed, and nothing here bears on it.

## Statements (Theorem 1/2 vs the print)

- (ii) and (iii) in §0 match KRTW's TeX (lines 2470–2500 of `KRTW.tex`, also LITERATURE §1b) letter for letter: summation
  ranges, exponents (r_i+s_i)²+s_i²(+r_i+2s_i), and both products. KRTW's convention (a;q)_n := (a;q)_∞/(aq^n;q)_∞ gives
  1/(q;q)_m = 0 for m < 0 (KRTW's own sentence says "nonnegative", a typo in KRTW), as the author uses.
- (iii) = K_k(1) and (ii) = K_k(q): correct (x^{r_i+2s_i} → 1, resp. q^{r_i+2s_i}).
- The reduction from Theorem 1 (KRTW-F-moves §3) was re-derived. The theta series (q)_∞B_{k+1}(1) = θ(q^{k+2};q^{2k+4}),
  (q)_∞B_{k+1}(q) = θ(q;q^{2k+4}), (Q;Q)_∞J_k(1) = θ(q^{k+1};q^{2k+4}) and (Q;Q)_∞J_k(q) = θ(q²;q^{2k+4}) are correct (the
  n ↦ −n−1 pairings are checked), and their products are exactly the right sides of (iii) and (ii). x ↦ 1, q is a ring
  homomorphism on Z[x][[q]], and all series involved lie in Z[x][[q]].

## Section verdicts

| § | item | verdict |
|---|------|---------|
| 1 | ring 𝔽 = Q(q)((x)); (1.1) K_k = ω(T^{k−1}u); (1.2), (1.3) | correct. Each x^D-coefficient is a rational function of q, regular at q = 0; the Taylor map is an injective ring map, and Z[x][[q]] ⊂ Q[[x,q]] = Q[[q]][[x]], so the transfer 𝔽 → Z[x][[q]] is justified. (1.1) also checked by brute force. |
| 2 | Prop 2.1; H̃′ = WH′W⁻¹, G̃; (2.1) | correct. The normal-ordered coefficients of H′E, EH̃′, GE and EG̃ were recomputed by hand and all equal the printed common values (the P²Y-coefficient of EG̃ uses the same identity 1−(1+q)(1−q^j)+q(1−q^j)(1−q^{j−1}) = q^{2j}). Restriction to N² is harmless (E1). |
| 3 | (B-i)–(B-v), the v_a rules | correct. (B-iii) holds at R = 0 (φ_e(−1) = 0). (B-v) holds for R ≥ 0, n ≥ 1, and genuinely fails for n = 0 (τφ_0 has infinite support; `ref_tau0.py`). Every use has n ≥ 1: Lemma 3.2 (τφ_{e+1}), Prop 4.1 (G built from φ_{e+1}), Lemma 6.2(a). h_a(τF) = F(·+1)+aq^{2R}F(·+1)−q⁻¹F needs F(0) = 0, which holds there. |
| 3 | (3.1), Prop 3.1, (3.2), Lemma 3.2 | correct. T is triangular (slice σ maps to slices S ≤ σ; within a slice it is diagonal, since Γ^{(0)} = id). The λ_{N,M} are distinct monomials. The support [L, L+2(M−S)] follows from the support of Γ. Uniqueness holds, using the top slice of a finite eigenvector. H′ and G send Ξ_{e,σ} to slices σ and σ−1 with finitely many terms (via φ_e(R;aq²) = (aq)_2φ_{e+1}(R+1;a) and (B-v) with n = e+1 ≥ 1). The top-slice eigenvalues ε_{N,M} and g_{N,M} were recomputed. |
| 4 | Prop 4.1, Cor 4.2, Lemma 4.3 | correct. The three v_x relations, the elimination (the bracket equals (e_N²−(1+x)²)(e_M²−(1+x)²)) and the normalisation v_1 = (xq)_2 W_{N−1,M−1}(xq²) were re-derived. Cor 4.2 was proved for all N, M via the ratio V_{N,M}(x)(xq)_{2M}/(V_{N−1,M−1}(xq²)(xq³)_{2M−2}) = the factor of Prop 4.1. |
| 5 | T♯ = W⁻¹TᵗW, (5.1), base-Q facts, Prop 5.1, Cor 5.2 (incl. L = 0), Lemma 5.3 | correct. (5.1) was derived independently: W⁻¹𝖲♮W = xqP²Y𝖲♮ and W⁻¹𝔲♮W = x⁻²q²P⁻²Y⁻²𝔲♮. The base-Q coefficients are c₁ = (1−aQ)/(1−aQ^{2e+1}), c₂ = −(1−aQ)aQ^{2e+1}/(1−aQ^{2e+1}) and b_e as printed. The Fricke factorisation and the product of Cor 5.2 were re-derived. At L = 0 the factor is an empty product, so V♯_{M,M} = 1. The transposed commutation T♯H′♯φ♯ = H′♯T♯φ♯ involves only finite sums (E3). |
| 6 | Lemma 6.1 | correct. Both sums are finite. The supports meet only at (L, M). |
| 6 | Lemma 6.2 | correct and proved. (a) The top index e = L+2(M−S) on slice S can arise only as e′+2 with e′ the top index of slice S+1 (through γ): reaching it through β or α would need e′ = e−1 or e, outside the support of slice S+1. There the bracket has order exactly 1, with x-coefficient ≠ 0 since S < M ≤ N. (b) The bracket ε_{ρ+e,e}−ε_{N,M} is nonzero for ρ < L, with order ≤ 1 (E4); c₁, c₂, b_e are x-integral. (c) correct. Numerically, the bound in (a) is even attained at every nonzero entry. |
| 6 | Lemma 6.3 | correct. m is attained on a finite nonempty set Z. With e* minimal and then σ* maximal, every other term (e ≤ e*, σ ≥ σ*) lies outside Z and so has order > m. The interchange of sums is x-adically legitimate: the order of the (r,s,e) term is ≥ r + 2s + inf ord Δ, with r ≤ L*. |
| 6 | Thm 6.4 | correct. ord ũ = L ≥ 0 (the author only needs ≥ 0). The (r,s,b)-families for ψ_a(Ũ), TŨ and ω(·) have order ≥ r + 2M_b with L_b ≤ r, so they are summable (E5). |
| 7 | Lemma 7.1, final assembly | correct. Lemma 7.1 was proved by hand for all N ≥ M ≥ 0 (N ≥ 1; N = M = 0 is trivial); the printed chain of Pochhammer identities is accurate. Symmetrisation and the use of (B1), (B2) are correct. |
| 0 | Corollary 3 | the deduction from Theorem 1 plus KRTW-F-moves Corollary C is correct. The cited GOW and BW theorems were not re-checked (outside the brief). |

## Problems found (all exposition or minor; no fatal problem, no substantive gap)

- **E1 (exposition, Prop 2.1).** Functions live on N² extended by 0. 𝔲 and E read and produce values off N²
  (e.g. (H′F)(R,−1) = −x²q^{R−2}F(R,0) ≠ 0). TH′ = H′T holds for the restricted operators because E reads only S′ ≥ S
  and R′ ≤ R, and EWF vanishes for R < 0. *Fix:* add this sentence. "Normal order is unique" is not needed; only
  "equal coefficients ⟹ equal operators on W𝒱" is used.
- **E2 (exposition, (B-iv)).** The support [d, d+c] and the sum Σ_eΓ = (aq)_c follow algebraically from partial fractions
  in z = q^R. Take the common denominator (q)_{R−d}(aq)_{R+d+c}: the numerator of φ_d(·;aq^c) is the constant (aq)_c, the
  numerators of φ_e are c+1 independent polynomials of degree c, and coefficient sums are values at z = 0. *Fix:* state this
  instead of the analytic "R → ∞" argument and the pointer to `code/tower.py`. (The printed Γ formula is correct; it was
  checked symbolically for c ≤ 3, d ≤ 2.)
- **E3 (exposition, §5).** "Transposing Prop 2.1" needs one line, because φ♯ has infinite S-support. The kernel identity
  (T♯H′♯)(a,b) = w(a)⁻¹w(b)(H′T)(b,a) holds with finite sums, and on φ♯ all sums are finite (H′♯ is local, T♯ is finite
  in S, φ♯ = 0 for R > L).
- **E4 (minor, Lemma 6.2(b)).** State that the bracket is *nonzero* for ρ < L, because its x-coefficient
  q^ρ+q^{−ρ}−q^L−q^{−L} ≠ 0. This is what lets the H′♯ equation determine c_{e,ρ}. "ord ≤ 1" implies it but is easy to misread.
- **E5 (minor, Thm 6.4).** "The same summability gives T^{k−1}u = Σũλ^{k−1}φ" and the interchange with ω are asserted
  only. *Fix:* add the bound: the (r,s,b) term has order ≥ r + 2M_b + (k−1)(N_b+M_b), with r ≥ L_b, so only finitely many
  terms lie below any order.
- **E6 (cosmetic).** In §1, "C_1 ≡ 1 mod x" plays no role in the transfer argument. Only the invertibility of C_1,
  needed to define J_k, is used.

## Independent computations (all in `referee/`, own code, exact unless noted)

1. `ref_commute.py`: [T,H′], [T,G], [T♯,H′♯], [T♯,G♯] on all point masses (r₀,s₀) ≤ (4,3), evaluated on boxes; and
   (5.1) vs W⁻¹H′ᵗW and W⁻¹GᵗW entrywise on [0,4]×[0,3]. Run at (q,x) = (2/7, 3/11) and (−3/5, 5/13): 3480 exact
   comparisons per point, 0 mismatches. A perturbed H′ is detected.
2. `ref_eigen.py 7` → `ref_eigen_N7.out`: φ_{N,M} and φ♯_{N,M} for all N ≤ 7, computed in **position space** by a
   pointwise triangular solve of Tφ = λφ (no Bailey recursion), at the same two points. Checked: H′/G eigen-relations
   (ε, g) and H′♯/G♯ eigen-relations; slice-0 Bailey coefficients supported in [L, L+2M] with sum (xq)_{2M}V_{N,M}
   (Cor 4.2); slice-R=0 base-Q coefficients supported in [M, M+L] with sum V♯_{N,M} (Cor 5.2); shift (3.2); all 1296
   pairings ψ_a(φ_b) (Lemma 6.1). Total 11322 exact checks per point, 0 failures. Numerically (summed exactly, then
   converted): ω(φ) = λW/(xq)_∞ and ψ(u) = λV♯/(x²Q;Q)_∞ to relative error 1e−72 (q = 2/7) and 3e−32 (q = −3/5),
   both truncation-limited. (q)_∞(xq)_∞φ(R,0) → W, consistent with an O(|q|^{R−L−2M}) error.
3. `ref_valuation.py 5`: ord_x estimated from exact values at x = 10⁻²⁵ and 10⁻⁵⁰, N ≤ 5, both q. ord φ_{N,M}(R,S) =
   2(M−S) exactly at all 560 nonzero entries (Lemma 6.2(a) is sharp), and ord φ♯ ≥ 0 (Lemma 6.2(b)).
4. `ref_expansion.py`: u(R,S) = Σũ_{N,M}φ_{N,M}(R,S) (Thm 6.4) for R ≤ 3, S ≤ 2, L ≤ 3, M ≤ 10, at (2/7, 1/50) and
   (−3/5, 1/40): relative error 1.5e−81. The error decreases with Mmax (1e−9, 1e−21, 1e−37, 1e−77 for Mmax = 2, 3, 4, 6),
   so the large-M terms do contribute.
5. `ref_series.py 24 120 8` → `ref_series_x24_q120.out`: bivariate series, exact, mod x²⁴, q¹²⁰, k = 1…8. For every k,
   K_k (via (1.1)) = the spectral formula of Thm 6.4 (with V, V♯) = Σλ^k t_{N,M}/((xq)_∞(x²Q;Q)_∞) = B_{k+1}C_k/C_1
   (from the multisum definitions of B and C). A perturbed spectral term is detected.
6. `ref_bruteK.py`: brute-force 2k-fold sum = ω(T^{k−1}u) mod x¹⁰, q⁴⁵ for k ≤ 5, confirming (1.1).
7. `ref_krtw.py 300 12` → `ref_krtw_q300.out`: KRTW (ii) and (iii) as q-series mod q³⁰⁰ for **k = 1…12**: all hold.
   The negative control ((iii) for k = 4 against the product for k = 5) is detected.
8. `ref_symbolic.py` → `ref_symbolic.out` (sympy), 24 checks, all OK: the Prop 2.1 normal-ordered coefficients
   (symbols q^j, Q^b, P, Y) and WH′ = H̃′W, WG = G̃W; the (B-iv) Γ formula, its coefficient sum, and
   φ_d(R;aq²) = (aq)_2φ_{d+1}(R+1;a); (B-v) α, β, γ symbolically and coefficient sums 1, 0; the Prop 4.1 elimination and
   factorisation; Cor 4.2 from the recursion (N ≤ 10; also proved for all N by hand); the Prop 5.1 elimination and the
   Fricke factorisation; Cor 5.2 from the recursion (N ≤ 10); Lemma 7.1 (N ≤ 10; also proved for all N by hand);
   t_{N,M} = λ(β^B_Nβ^J_M + β^B_Mβ^J_N)/(1+δ).
9. `ref_tau0.py`: τφ_n has Bailey support {n−1, n, n+1} for n = 1, 2, 3, while τφ_0 has infinite support (boundary case of (B-v)).

## Not checked / outside the brief

- (B1) and (B2) themselves (previously refereed; taken as known). Only their use was checked, including the
  evaluations at x = 1 and x = q.
- Novelty and the literature claim "k ≤ 3 previously known": not re-searched.
- The theorems cited in Corollary 3 (GOW Thm 1.3, Bartlett–Warnaar C₂ Rogers–Selberg).
