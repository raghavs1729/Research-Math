/-
Nahm/Defs.lean — the definitions that occur in the final statements.  Imports only Mathlib.
(These definitions are inlined verbatim in `NahmStatement.lean`.)
-/
import Mathlib

namespace Nahm

open PowerSeries Finset

-- BEGIN DEFINITIONS
variable {K : Type*} [Field K]

/-- `(q;q)_n = (1 - q)(1 - q²)⋯(1 - q^n)`, as a formal power series in `q = X`. -/
noncomputable def qPoch (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (i + 1))

/-- `(q^a; q^d)_∞ = ∏_{i ≥ 0} (1 - q^{a+di})` for `a, d ≥ 1`, defined by its coefficients: the
coefficient of `q^N` is that of the finite product over `i ≤ N` (the later factors are
`≡ 1 mod q^{N+1}`). -/
noncomputable def qPochInf (a d : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => coeff N (∏ i ∈ range (N + 1), (1 - X ^ (a + d * i) : K⟦X⟧))

/-- All tuples of length `k` with entries in `{0, 1, …, N}`. -/
def box (k N : ℕ) : Finset (Fin k → ℕ) := Fintype.piFinset fun _ => range (N + 1)

/-- For a tuple `n = (n_1, …, n_r)` (stored as `n 0, …, n (r-1)`): the partial sum
`N_j = n_j + n_{j+1} + ⋯ + n_r` (0-based: `psum n j = Σ_{l ≥ j} n l`). -/
def psum {r : ℕ} (n : Fin r → ℕ) (j : Fin r) : ℕ := ∑ l ∈ univ.filter (fun l => j ≤ l), n l

/-- The Andrews–Gordon exponent `N_1² + ⋯ + N_{k-1}² + N_i + ⋯ + N_{k-1}` of a tuple
`n = (n_1, …, n_{k-1})` (0-based index `j` stands for `j + 1`; the linear part runs over `j + 1 ≥ i`). -/
def agExp (k i : ℕ) (n : Fin (k - 1) → ℕ) : ℕ :=
  ∑ j, psum n j ^ 2 + ∑ j ∈ univ.filter (fun j : Fin (k - 1) => i ≤ (j : ℕ) + 1), psum n j

/-- The Andrews–Gordon multiple sum
`Σ_{n_1, …, n_{k-1} ≥ 0} q^{N_1² + ⋯ + N_{k-1}² + N_i + ⋯ + N_{k-1}} / ((q;q)_{n_1} ⋯ (q;q)_{n_{k-1}})`,
defined by its coefficients: the coefficient of `q^N` is the finite sum, over the tuples with all
entries `≤ N`, of the coefficients of `q^N` of the summands (a summand with an entry `n_j > N` is
divisible by `q^{N_1²}` with `N_1² ≥ N_1 ≥ n_j > N`). -/
noncomputable def agSum (k i : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ n ∈ box (k - 1) N, coeff N (X ^ agExp k i n * ∏ j, (qPoch (n j))⁻¹)

/-- The Andrews–Gordon product `∏_{n ≥ 1, n ≢ 0, ±i (mod 2k+1)} (1 - q^n)⁻¹`, written as
`1 / ∏_{a} (q^a; q^{2k+1})_∞` over the residues `1 ≤ a ≤ 2k` with `a ≠ i`, `a ≠ 2k+1-i`. -/
noncomputable def agProduct (k i : ℕ) : K⟦X⟧ :=
  (∏ a ∈ (Icc 1 (2 * k)).filter (fun a => a ≠ i ∧ a ≠ 2 * k + 1 - i), qPochInf a (2 * k + 1))⁻¹
/-- The inverse `C(E8)⁻¹` (over `ℚ`) of the Cartan matrix of E8; `CartanMatrix.E₈` is Mathlib's Cartan
matrix of E8 in Bourbaki's labelling (chain `1 — 3 — 4 — 5 — 6 — 7 — 8`, node `2` attached to node `4`;
0-based in Lean). -/
noncomputable def invCartanE8 : Matrix (Fin 8) (Fin 8) ℚ := (Matrix.map CartanMatrix.E₈ (Int.cast : ℤ → ℚ))⁻¹

/-- The quadratic form `½ nᵀ C(E8)⁻¹ n + B·n` (a rational number). -/
noncomputable def e8Form (B : Fin 8 → ℕ) (n : Fin 8 → ℕ) : ℚ :=
  (1 / 2 : ℚ) * dotProduct (fun i => (n i : ℚ)) (Matrix.mulVec invCartanE8 fun i => (n i : ℚ))
    + dotProduct (fun i => (B i : ℚ)) fun i => (n i : ℚ)

/-- The exponent `½ nᵀ C(E8)⁻¹ n + B·n` as a natural number (`⌊·⌋₊`; it IS a natural number, see
`e8Form_eq_e8Exp`). -/
noncomputable def e8Exp (B : Fin 8 → ℕ) (n : Fin 8 → ℕ) : ℕ := ⌊e8Form B n⌋₊

/-- The E8 Nahm sum `Σ_{n ∈ ℕ^8} q^{½ nᵀ C(E8)⁻¹ n + B·n} / ∏_{i=1}^{8} (q;q)_{n_i}`, defined by its
coefficients (the coefficient of `q^N` is the finite sum over the tuples with all entries `≤ N`; a
summand with an entry `n_i > N` has exponent `≥ n_i² ≥ n_i > N`). -/
noncomputable def e8Sum (B : Fin 8 → ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ n ∈ box 8 N, coeff N (X ^ e8Exp B n * ∏ i, (qPoch (n i))⁻¹)
/-- `(q;q)_n` for `q = X^4`, i.e. as a power series in `X = q^{1/4}`: `∏_{i<n} (1 - X^{4(i+1)})`. -/
noncomputable def qPoch4 (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (4 * (i + 1)))

/-- The inverse `C(E7)⁻¹` (over `ℚ`) of the Cartan matrix of E7; `CartanMatrix.E₇` is Mathlib's Cartan
matrix of E7 in Bourbaki's labelling (chain `1 — 3 — 4 — 5 — 6 — 7`, node `2` attached to node `4`;
0-based in Lean). -/
noncomputable def invCartanE7 : Matrix (Fin 7) (Fin 7) ℚ := (Matrix.map CartanMatrix.E₇ (Int.cast : ℤ → ℚ))⁻¹

/-- The quadratic form `½ nᵀ C(E7)⁻¹ n` (a rational number with denominator dividing 4). -/
noncomputable def e7Form (n : Fin 7 → ℕ) : ℚ :=
  (1 / 2 : ℚ) * dotProduct (fun i => (n i : ℚ)) (Matrix.mulVec invCartanE7 fun i => (n i : ℚ))

/-- The exponent `4 · ½ nᵀ C(E7)⁻¹ n` of `X = q^{1/4}`, as a natural number (`⌊·⌋₊`; it IS a natural
number, see `e7Form_eq_e7Exp`). -/
noncomputable def e7Exp (n : Fin 7 → ℕ) : ℕ := ⌊4 * e7Form n⌋₊

/-- The E7 Nahm sum `Σ_{n ∈ ℕ^7} q^{½ nᵀ C(E7)⁻¹ n} / ∏_{i=1}^{7} (q;q)_{n_i}` as a power series in
`X = q^{1/4}` (`q = X^4`), defined by its coefficients (the coefficient of `X^N` is the finite sum over
the tuples with all entries `≤ N`; a summand with an entry `n_i > N` has exponent `≥ 3 n_i² > N`). -/
noncomputable def e7Sum : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ n ∈ box 7 N, coeff N (X ^ e7Exp n * ∏ i, (qPoch4 (n i))⁻¹)

/-- The bilateral theta series `Σ_{k ∈ ℤ} (q^{3k²/4} + q^{7k²/4})` as a power series in `X = q^{1/4}`:
`Σ_{k ∈ ℤ} (X^{3k²} + X^{7k²})`, defined by its coefficients (only `|k| ≤ N` contributes to `X^N`). -/
noncomputable def e7Theta : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ k ∈ Icc (-(N : ℤ)) N, coeff N (X ^ (3 * k.natAbs ^ 2) + X ^ (7 * k.natAbs ^ 2))

/-- The one-sided form `Σ_{k ≥ 0} q^{3k²/4} + Σ_{k ≥ 1} q^{7k²/4}` (in powers of `X = q^{1/4}`). -/
noncomputable def e7ThetaOneSided : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ k ∈ range (N + 1), coeff N (X ^ (3 * k ^ 2) + if k = 0 then 0 else X ^ (7 * k ^ 2))
-- END DEFINITIONS

end Nahm
