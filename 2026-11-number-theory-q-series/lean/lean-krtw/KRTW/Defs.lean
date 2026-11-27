/-
KRTW/Defs.lean — the definitions that occur in the final statements.  Imports only Mathlib.
(These definitions are inlined verbatim in `KRTWStatement.lean`.)
-/
import Mathlib

namespace KRTW

open PowerSeries Finset

-- BEGIN DEFINITIONS
variable {K : Type*} [Field K]

/-- `(q;q)_n = (1 - q)(1 - q²)⋯(1 - q^n)`, as a formal power series in `q = X`. -/
noncomputable def qPoch (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (i + 1))

/-- `(q²;q²)_n = (1 - q²)(1 - q⁴)⋯(1 - q^{2n})`. -/
noncomputable def q2Poch (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (2 * (i + 1)))

/-- For a tuple `r = (r_1, …, r_k)` (stored as `r 0, …, r (k-1)`): the entry following position `i`,
with the convention `r_{k+1} := 0`. -/
def nxt {k : ℕ} (r : Fin k → ℕ) (i : Fin k) : ℕ :=
  if h : (i : ℕ) + 1 < k then r ⟨i + 1, h⟩ else 0

/-- For a tuple `s = (s_1, …, s_k)`: the entry preceding position `i`, with the convention
`s_0 := 0`. -/
def prv {k : ℕ} (s : Fin k → ℕ) (i : Fin k) : ℕ :=
  if h : (i : ℕ) = 0 then 0 else s ⟨i - 1, by omega⟩

/-- All tuples of length `k` with entries in `{0, 1, …, N}`. -/
def box (k N : ℕ) : Finset (Fin k → ℕ) := Fintype.piFinset fun _ => range (N + 1)

/-- The summand of the (`x`-refined) Kanade–Russell–Tsuchioka–Warnaar sum:
`∏_{i=1}^k x^{r_i+2s_i} q^{(r_i+s_i)²+s_i²+s_i} / ((q;q)_{r_i-r_{i+1}} (q²;q²)_{s_i-s_{i-1}})`. -/
noncomputable def krtwTerm (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) : K⟦X⟧ :=
  ∏ i : Fin k, x ^ (r i + 2 * s i) * X ^ ((r i + s i) ^ 2 + s i ^ 2 + s i)
    * (qPoch (r i - nxt r i))⁻¹ * (q2Poch (s i - prv s i))⁻¹

/-- The `2k`-fold KRTW sum over `r_1 ≥ … ≥ r_k ≥ 0` and `0 ≤ s_1 ≤ … ≤ s_k`, defined by its
coefficients: the coefficient of `q^N` is the finite sum of the coefficients of `q^N` of the summands
with all `r_i, s_i ≤ N` (the other summands are divisible by `q^{N+1}`, see
`coeff_krtwTerm_eq_zero`). -/
noncomputable def krtwSum (x : K⟦X⟧) (k : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ r ∈ box k N, ∑ s ∈ box k N,
    if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prv s i ≤ s i) then coeff N (krtwTerm x r s) else 0

/-- The summand of the (`x`-refined) `k`-fold "Melzer-type" sum:
`∏_{i=1}^k x^{M_i} q^{M_i(M_i+1)/2} / (q;q)_{M_i-M_{i+1}}`. -/
noncomputable def melzerTerm (x : K⟦X⟧) {k : ℕ} (M : Fin k → ℕ) : K⟦X⟧ :=
  ∏ i : Fin k, x ^ (M i) * X ^ (M i * (M i + 1) / 2) * (qPoch (M i - nxt M i))⁻¹

/-- The `k`-fold sum over `M_1 ≥ … ≥ M_k ≥ 0`, defined by its coefficients (as `krtwSum`). -/
noncomputable def melzerSum (x : K⟦X⟧) (k : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ M ∈ box k N,
    if ∀ i, nxt M i ≤ M i then coeff N (melzerTerm x M) else 0

/-- `(q^a; q^d)_∞ = ∏_{i ≥ 0} (1 - q^{a+di})` for `a, d ≥ 1`, defined by its coefficients: the
coefficient of `q^N` is that of the finite product over `i ≤ N` (the later factors are
`≡ 1 mod q^{N+1}`). -/
noncomputable def qPochInf (a d : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => coeff N (∏ i ∈ range (N + 1), (1 - X ^ (a + d * i) : K⟦X⟧))

/-- The theta series `Σ_{n ∈ ℤ} (-1)^n q^{((k+2)n² + kn)/2}`, written as a sum over `n ≥ 0` by pairing
the terms `n` and `-n-1` (whose exponents differ by `2n+1`):
`Σ_{n ≥ 0} (-1)^n q^{((k+2)n² + kn)/2} (1 - q^{2n+1})`.  Defined by its coefficients: the terms with
`n > N` do not contribute to the coefficient of `q^N`. -/
noncomputable def thetaSum (k : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ n ∈ range (N + 1),
    coeff N ((-1 : K⟦X⟧) ^ n * X ^ (((k + 2) * n ^ 2 + k * n) / 2) * (1 - X ^ (2 * n + 1)))
-- END DEFINITIONS

end KRTW
