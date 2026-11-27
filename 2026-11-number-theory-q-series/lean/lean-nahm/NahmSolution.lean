/-
NahmSolution.lean — the UNTRUSTED solution module: it proves, from the library `Nahm/`, exactly the
statements of `NahmStatement.lean` (same names, same docstrings, same statements; checked by
`scripts/compare_statements.sh`).  No `sorry`, no `axiom`, no `native_decide`.
The definitions used in the statements are those of `Nahm/Defs.lean` (Mathlib only), which are
inlined verbatim in `NahmStatement.lean`.
-/
import Nahm
import Nahm.Final

open PowerSeries Finset

namespace Nahm

/-- **The Andrews–Gordon identities** (G. E. Andrews, PNAS 71 (1974) 4082–4085, Theorem 1;
*The Theory of Partitions*, Theorem 7.8): for `k ≥ 1` and `1 ≤ i ≤ k`, in `ℚ[[q]]`,

  `Σ_{n_1, …, n_{k-1} ≥ 0} q^{N_1² + ⋯ + N_{k-1}² + N_i + ⋯ + N_{k-1}} / ((q;q)_{n_1} ⋯ (q;q)_{n_{k-1}})
      = ∏_{n ≥ 1, n ≢ 0, ±i (mod 2k+1)} (1 - q^n)⁻¹`,

`N_j = n_j + ⋯ + n_{k-1}`.  `agSum k i` is the left side (defined by its coefficients), `agProduct k i`
the right side, written as `1/∏_a (q^a; q^{2k+1})_∞` over the residues `1 ≤ a ≤ 2k`, `a ≠ i, 2k+1-i`. -/
theorem andrews_gordon (k i : ℕ) (hi : 1 ≤ i) (hik : i ≤ k) :
    (agSum k i : ℚ⟦X⟧) = agProduct k i := by
  exact andrews_gordon_top' k i hi hik

/-- The Andrews–Gordon identities over an arbitrary field `K` (in place of `ℚ`). -/
theorem andrews_gordon_field {K : Type*} [Field K] (k i : ℕ) (hi : 1 ≤ i) (hik : i ≤ k) :
    (agSum k i : K⟦X⟧) = agProduct k i := by
  exact andrews_gordon_top' k i hi hik

/-- The summands of the Andrews–Gordon sum with an index outside `{0, …, N}` do not contribute to the
coefficient of `q^N`; so the finite sum in the definition of `agSum` is the coefficient of `q^N` of
the whole (infinite) sum. -/
theorem coeff_agTerm_eq_zero {K : Type*} [Field K] (k i : ℕ) (n : Fin (k - 1) → ℕ) (N : ℕ)
    (h : n ∉ box (k - 1) N) : coeff N ((X : K⟦X⟧) ^ agExp k i n * ∏ j, (qPoch (n j))⁻¹) = 0 := by
  exact coeff_agTerm_eq_zero' k i n N h

/-- The inverse of the Cartan matrix of E8 (Bourbaki labelling, Mathlib's `CartanMatrix.E₈`). -/
theorem invCartanE8_eq_explicit : invCartanE8 =
    !![(4 : ℚ), 5, 7, 10, 8, 6, 4, 2;
       5, 8, 10, 15, 12, 9, 6, 3;
       7, 10, 14, 20, 16, 12, 8, 4;
       10, 15, 20, 30, 24, 18, 12, 6;
       8, 12, 16, 24, 20, 15, 10, 5;
       6, 9, 12, 18, 15, 12, 8, 4;
       4, 6, 8, 12, 10, 8, 6, 3;
       2, 3, 4, 6, 5, 4, 3, 2] := by
  exact invCartanE8_explicit

/-- **The E8 exponents are integers**: for every `B ∈ ℕ^8` and `n ∈ ℕ^8`,
`½ nᵀ C(E8)⁻¹ n + B·n` is the natural number `e8Exp B n`. -/
theorem e8Form_eq_e8Exp (B n : Fin 8 → ℕ) : (e8Exp B n : ℚ) = e8Form B n := by
  exact e8Form_eq_e8Exp' B n

/-- The summands of the E8 Nahm sum with an index outside `{0, …, N}` do not contribute to the
coefficient of `q^N`. -/
theorem coeff_e8Term_eq_zero {K : Type*} [Field K] (B n : Fin 8 → ℕ) (N : ℕ) (h : n ∉ box 8 N) :
    coeff N ((X : K⟦X⟧) ^ e8Exp B n * ∏ i, (qPoch (n i))⁻¹) = 0 := by
  exact coeff_e8Term_eq_zero' B n N h

/-- **The E8 vacuum identity** (Sun–Wang, arXiv:2604.00847, eq. (3.20); `proofs/E8-E6.md`, (0.4)):
in `ℚ[[q]]`,

  `Σ_{n ∈ ℕ^8} q^{½ nᵀ C(E8)⁻¹ n} / ∏_{i=1}^{8} (q;q)_{n_i} = 1/(q, q², q³, q⁴, q⁷, q⁸, q⁹, q¹⁰; q¹¹)_∞`. -/
theorem e8_vacuum :
    (e8Sum 0 : ℚ⟦X⟧) = (∏ a ∈ ({1, 2, 3, 4, 7, 8, 9, 10} : Finset ℕ), qPochInf a 11)⁻¹ := by
  exact e8_vacuum_top'

/-- The E8 vacuum identity over an arbitrary field `K`. -/
theorem e8_vacuum_field {K : Type*} [Field K] :
    (e8Sum 0 : K⟦X⟧) = (∏ a ∈ ({1, 2, 3, 4, 7, 8, 9, 10} : Finset ℕ), qPochInf a 11)⁻¹ := by
  exact e8_vacuum_top'

/-- **E8 companion `i = 4`** (`proofs/E8-E6.md`, (0.3); `B^{(4)} = (0,0,0,1,1,1,0,1)` in the paper's
labelling = `(0,1,1,1,1,0,0,0)` in Bourbaki's):
`Σ_n q^{½ nᵀ C(E8)⁻¹ n + B^{(4)}·n} / ∏ (q;q)_{n_i} = ∏_{n ≢ 0, ±4 (mod 11)} (1 - q^n)⁻¹`. -/
theorem e8_companion_4 : (e8Sum ![0, 1, 1, 1, 1, 0, 0, 0] : ℚ⟦X⟧) = agProduct 5 4 := by
  rw [← Bb_agb4]; exact e8_companion_top' 4 (by norm_num) (by norm_num)

/-- **E8 companion `i = 3`** (`B^{(3)} = (0,1,1,2,2,2,0,2)` in the paper's labelling
= `(0,2,2,2,2,1,1,0)` in Bourbaki's): `… = ∏_{n ≢ 0, ±3 (mod 11)} (1 - q^n)⁻¹`. -/
theorem e8_companion_3 : (e8Sum ![0, 2, 2, 2, 2, 1, 1, 0] : ℚ⟦X⟧) = agProduct 5 3 := by
  rw [← Bb_agb3]; exact e8_companion_top' 3 (by norm_num) (by norm_num)

/-- **E8 companion `i = 2`** (`B^{(2)} = (0,2,2,4,4,3,1,3)` in the paper's labelling
= `(1,3,3,4,4,2,2,0)` in Bourbaki's): `… = ∏_{n ≢ 0, ±2 (mod 11)} (1 - q^n)⁻¹`. -/
theorem e8_companion_2 : (e8Sum ![1, 3, 3, 4, 4, 2, 2, 0] : ℚ⟦X⟧) = agProduct 5 2 := by
  rw [← Bb_agb2]; exact e8_companion_top' 2 (by norm_num) (by norm_num)

/-- **E8 companion `i = 1`** (`B^{(1)} = (1,3,4,6,7,5,2,4)` in the paper's labelling
= `(2,4,5,7,6,4,3,1)` in Bourbaki's): `… = ∏_{n ≢ 0, ±1 (mod 11)} (1 - q^n)⁻¹`. -/
theorem e8_companion_1 : (e8Sum ![2, 4, 5, 7, 6, 4, 3, 1] : ℚ⟦X⟧) = agProduct 5 1 := by
  rw [← Bb_agb1]; exact e8_companion_top' 1 (by norm_num) (by norm_num)

/-- The inverse of the Cartan matrix of E7 (Bourbaki labelling, Mathlib's `CartanMatrix.E₇`). -/
theorem invCartanE7_eq_explicit : invCartanE7 =
    !![(2 : ℚ), 2, 3, 4, 3, 2, 1;
       2, 7/2, 4, 6, 9/2, 3, 3/2;
       3, 4, 6, 8, 6, 4, 2;
       4, 6, 8, 12, 9, 6, 3;
       3, 9/2, 6, 9, 15/2, 5, 5/2;
       2, 3, 4, 6, 5, 4, 2;
       1, 3/2, 2, 3, 5/2, 2, 3/2] := by
  exact invCartanE7_explicit

/-- **The E7 exponents are quarter-integers**: for every `n ∈ ℕ^7`, `4 · ½ nᵀ C(E7)⁻¹ n` is the
natural number `e7Exp n`. -/
theorem e7Form_eq_e7Exp (n : Fin 7 → ℕ) : (e7Exp n : ℚ) = 4 * e7Form n := by
  exact e7Form_eq_e7Exp' n

/-- The summands of the E7 Nahm sum with an index outside `{0, …, N}` do not contribute to the
coefficient of `X^N`. -/
theorem coeff_e7Term_eq_zero {K : Type*} [Field K] (n : Fin 7 → ℕ) (N : ℕ) (h : n ∉ box 7 N) :
    coeff N ((X : K⟦X⟧) ^ e7Exp n * ∏ i, (qPoch4 (n i))⁻¹) = 0 := by
  exact coeff_e7Term_eq_zero' n N h

/-- The terms `k ∉ [-N, N]` of the theta series do not contribute to the coefficient of `X^N`. -/
theorem coeff_e7ThetaTerm_eq_zero {K : Type*} [Field K] (k : ℤ) (N : ℕ) (h : k ∉ Icc (-(N : ℤ)) N) :
    coeff N ((X : K⟦X⟧) ^ (3 * k.natAbs ^ 2) + X ^ (7 * k.natAbs ^ 2)) = 0 := by
  exact coeff_e7ThetaTerm_eq_zero' k N h

/-- **The E7 Nahm-sum identity** (`proofs/E7.md`, Theorem 4), one-sided form: as power series in
`X = q^{1/4}` with rational coefficients,

  `Σ_{n ∈ ℕ^7} q^{½ nᵀ C(E7)⁻¹ n} / ∏_{i=1}^{7} (q;q)_{n_i}
      = (1/(q;q)_∞) (Σ_{k ≥ 0} q^{3k²/4} + Σ_{k ≥ 1} q^{7k²/4})`

(`qPochInf 4 4 = (X^4; X^4)_∞ = (q;q)_∞`). -/
theorem e7_identity : (e7Sum : ℚ⟦X⟧) = (qPochInf 4 4)⁻¹ * e7ThetaOneSided := by
  exact e7_identity_top'

/-- **The E7 Nahm-sum identity**, bilateral form (as printed):

  `Σ_{n ∈ ℕ^7} q^{½ nᵀ C(E7)⁻¹ n} / ∏ (q;q)_{n_i} = (1/(2 (q;q)_∞)) Σ_{k ∈ ℤ} (q^{3k²/4} + q^{7k²/4})`,

stated without division as `2 · (left side) = (q;q)_∞⁻¹ · Σ_{k ∈ ℤ} (q^{3k²/4} + q^{7k²/4})`. -/
theorem e7_identity_bilateral : 2 * (e7Sum : ℚ⟦X⟧) = (qPochInf 4 4)⁻¹ * e7Theta := by
  exact e7_identity_bilateral_top'

/-- The E7 identity (one-sided form) over an arbitrary field `K`. -/
theorem e7_identity_field {K : Type*} [Field K] :
    (e7Sum : K⟦X⟧) = (qPochInf 4 4)⁻¹ * e7ThetaOneSided := by
  exact e7_identity_top'

/-- The E7 identity (bilateral form) over an arbitrary field `K`. -/
theorem e7_identity_bilateral_field {K : Type*} [Field K] :
    2 * (e7Sum : K⟦X⟧) = (qPochInf 4 4)⁻¹ * e7Theta := by
  exact e7_identity_bilateral_top'

/-- Stability of the truncated products: `qPochInf a d` has the coefficients of the infinite product
`∏_{i ≥ 0} (1 - q^{a+di})` (every finite product over `i ≤ M` with `M ≥ N` has the same coefficient
of `q^N`). -/
theorem coeff_qPochInf {K : Type*} [Field K] (a d : ℕ) (ha : 1 ≤ a) (hd : 1 ≤ d) (N M : ℕ)
    (h : N ≤ M) :
    coeff N (qPochInf a d : K⟦X⟧)
      = coeff N (∏ i ∈ range (M + 1), (1 - X ^ (a + d * i) : K⟦X⟧)) := by
  exact KRTW.coeff_qPochInf' a d ha hd N M h

end Nahm
