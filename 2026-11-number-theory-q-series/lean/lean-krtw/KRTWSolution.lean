/-
KRTWSolution.lean — the UNTRUSTED solution module: it proves, from the library `KRTW/`, exactly the
statements of `KRTWStatement.lean` (same names, same docstrings, same statements; checked by
`scripts/compare_statements.sh`).  No `sorry`, no `axiom`, no `native_decide`.
The definitions used in the statements are those of `KRTW/Defs.lean` (Mathlib only), which are
inlined verbatim in `KRTWStatement.lean`.
-/
import KRTW

open PowerSeries Finset

namespace KRTW

/-- **The reduction of the `2k`-fold Kanade–Russell–Tsuchioka–Warnaar sum to a `k`-fold sum,
with the refinement by `x`** (Theorem 2, first equality, of `proofs/KRTW-general.md`; this is the
new part of the proof of the conjecture).  For every `k ≥ 0`,

  `Σ_{r_1 ≥ … ≥ r_k ≥ 0, 0 ≤ s_1 ≤ … ≤ s_k}  x^{Σ_i (r_i + 2 s_i)}
        ∏_{i=1}^k q^{(r_i+s_i)² + s_i² + s_i} / ( (q;q)_{r_i - r_{i+1}} (q²;q²)_{s_i - s_{i-1}} )`
  `= Σ_{M_1 ≥ … ≥ M_k ≥ 0}  x^{M_1 + ⋯ + M_k}  ∏_{i=1}^k q^{M_i(M_i+1)/2} / (q;q)_{M_i - M_{i+1}}`

(`r_{k+1} = s_0 = M_{k+1} := 0`), as formal power series in `q`.

In Lean: `K` is an arbitrary field, `q = X` is the variable of `K⟦X⟧`, and `x` is an ARBITRARY element
of `K⟦X⟧`.  Taking `K = ℚ(x)` (or the fraction field of any ring containing `ℤ[x]`) and `x` the
constant series `x` gives the identity in `ℤ[x][[q]]`; taking `x = 1` gives the unrefined identity
(`krtw_sum_eq_melzer_sum`).  `krtwSum x k` and `melzerSum x k` are defined by their coefficients:
the coefficient of `q^N` is the finite sum, over the index tuples with all entries `≤ N`, of the
coefficient of `q^N` of the summand (`krtwTerm`, `melzerTerm`); the tuples with an entry `> N` do
not contribute (`coeff_krtwTerm_eq_zero`, `coeff_melzerTerm_eq_zero`). -/
theorem krtw_sum_eq_melzer_sum_refined {K : Type*} [Field K] (x : K⟦X⟧) (k : ℕ) :
    krtwSum x k = melzerSum x k := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact krtwSum_eq_melzerSum' x k

/-- **The `2k`-fold KRTW sum equals a `k`-fold sum** (`x = 1` in `krtw_sum_eq_melzer_sum_refined`):
for every `k ≥ 0`, in `ℚ[[q]]`,

  `Σ_{r_1 ≥ … ≥ r_k ≥ 0, 0 ≤ s_1 ≤ … ≤ s_k}
        ∏_{i=1}^k q^{(r_i+s_i)² + s_i² + s_i} / ( (q;q)_{r_i - r_{i+1}} (q²;q²)_{s_i - s_{i-1}} )`
  `= Σ_{M_1 ≥ … ≥ M_k ≥ 0}  ∏_{i=1}^k q^{M_i(M_i+1)/2} / (q;q)_{M_i - M_{i+1}}`

(`r_{k+1} = s_0 = M_{k+1} := 0`).  The left side is the left side of the first conjectural
Andrews–Gordon-type identity of S. Kanade, M. C. Russell, S. Tsuchioka and S. O. Warnaar,
*Remarks on the conjectures of Capparelli, Meurman, Primc and Primc*, arXiv:2404.03851v2 (displayed
after Conjecture 4.9; there the sum runs over all `r_i, s_i ≥ 0` with `1/(q)_n := 0` for `n < 0`,
which forces `r_1 ≥ … ≥ r_k` and `s_1 ≤ … ≤ s_k`).  The right side is a sum whose product formula is
known (Andrews; Kim–Yee; Chern–Li–Stanton–Xue–Yee; Kanade–Lovejoy, arXiv:2410.15408, Theorem 1.1
at `a = 0`; Warnaar, arXiv:2511.17034). -/
theorem krtw_sum_eq_melzer_sum (k : ℕ) :
    krtwSum (1 : ℚ⟦X⟧) k = melzerSum (1 : ℚ⟦X⟧) k := by
  exact krtw_sum_eq_melzer_sum_refined 1 k

/-- The summands of the KRTW sum with an index outside `{0, …, N}` do not contribute to the
coefficient of `q^N` (they are divisible by `q^{N+1}`); so the finite sum in the definition of
`krtwSum` is the coefficient of `q^N` of the whole (infinite) sum. -/
theorem coeff_krtwTerm_eq_zero {K : Type*} [Field K] (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) (N : ℕ)
    (h : r ∉ box k N ∨ s ∉ box k N) : coeff N (krtwTerm x r s) = 0 := by
  exact coeff_krtwTerm_eq_zero' x r s N h

/-- The summands of the `k`-fold sum with an index outside `{0, …, N}` do not contribute to the
coefficient of `q^N`. -/
theorem coeff_melzerTerm_eq_zero {K : Type*} [Field K] (x : K⟦X⟧) {k : ℕ} (M : Fin k → ℕ) (N : ℕ)
    (h : M ∉ box k N) : coeff N (melzerTerm x M) = 0 := by
  exact coeff_melzerTerm_eq_zero' x M N h

lemma krtw_product_field {K : Type*} [Field K] (k : ℕ) :
    krtwSum (1 : K⟦X⟧) k = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
      * (qPochInf 1 2 * qPochInf 1 1)⁻¹ := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact krtw_product' k

lemma krtw_theta_field {K : Type*} [Field K] (k : ℕ) :
    krtwSum (1 : K⟦X⟧) k * (qPochInf 1 2 * qPochInf 1 1) = thetaSum k := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact krtw_theta' k

/-- **The first Andrews–Gordon-type identity of Kanade, Russell, Tsuchioka and Warnaar, for every
`k`.**  In `ℚ[[q]]`, for every integer `k ≥ 0` (the case of interest is `k ≥ 1`; `k = 0` is `1 = 1`):

  `Σ_{r_1 ≥ … ≥ r_k ≥ 0, 0 ≤ s_1 ≤ … ≤ s_k}
        ∏_{i=1}^k q^{(r_i+s_i)² + s_i² + s_i} / ( (q;q)_{r_i - r_{i+1}} (q²;q²)_{s_i - s_{i-1}} )`
  `= (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ / ( (q;q²)_∞ (q;q)_∞ )`          (`r_{k+1} = s_0 := 0`).

This is the first of the four "conjectural Andrews–Gordon-type identities" displayed after
Conjecture 4.9 of S. Kanade, M. C. Russell, S. Tsuchioka, S. O. Warnaar, *Remarks on the conjectures
of Capparelli, Meurman, Primc and Primc*, arXiv:2404.03851v2 (Selecta Math., 2026); it is
conjectural there for `k ≥ 3`.  For `k ≤ 2` it follows from results there (for `k = 2`: Theorem 4.10
and identity (4.7), with the theorem of Meurman and Primc).  (There the sum runs over all `r_i, s_i ≥ 0`
with the convention `1/(q)_n := 0` for `n < 0`; that is the same sum.)

In Lean: `q = X`; `krtwSum 1 k` is the left side (defined by its coefficients, see `krtwSum`);
`qPochInf a d` is `(q^a; q^d)_∞ = ∏_{i ≥ 0} (1 - q^{a+di})` (defined by its coefficients, see
`qPochInf` and `coeff_qPochInf`), so the right side is
`(q;q^{k+2})_∞ (q^{k+1};q^{k+2})_∞ (q^{k+2};q^{k+2})_∞ · ((q;q²)_∞ (q;q)_∞)⁻¹`, where `⁻¹` is the
inverse in `ℚ⟦X⟧` of a series with constant term `1`. -/
theorem krtw_identity (k : ℕ) :
    krtwSum (1 : ℚ⟦X⟧) k
      = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
        * (qPochInf 1 2 * qPochInf 1 1)⁻¹ := by
  exact krtw_product_field k

/-- The KRTW identity in theta form (the output of the Bailey chain, before the Jacobi triple
product): for every `k ≥ 0`, in `ℚ[[q]]`,

  `(2k-fold KRTW sum) · (q;q²)_∞ (q;q)_∞ = Σ_{n ∈ ℤ} (-1)^n q^{((k+2)n² + kn)/2}`,

the right side being written as `Σ_{n ≥ 0} (-1)^n q^{((k+2)n²+kn)/2} (1 - q^{2n+1})` (`thetaSum k`:
the terms `n` and `-n-1` of the bilateral sum are paired). -/
theorem krtw_identity_theta (k : ℕ) :
    krtwSum (1 : ℚ⟦X⟧) k * (qPochInf 1 2 * qPochInf 1 1) = thetaSum k := by
  exact krtw_theta_field k

/-- The instance of the Jacobi triple product that is used:
`Σ_{n ∈ ℤ} (-1)^n q^{((k+2)n² + kn)/2} = (q, q^{k+1}, q^{k+2}; q^{k+2})_∞` (left side: `thetaSum k`). -/
theorem theta_eq_triple_product (k : ℕ) :
    (thetaSum k : ℚ⟦X⟧) = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2) := by
  exact (jtp_theta k).symm

/-- The (known) evaluation of the `k`-fold sum (Andrews; Kim–Yee; Chern–Li–Stanton–Xue–Yee;
Kanade–Lovejoy, arXiv:2410.15408, Theorem 1.1 at `a = 0`; Warnaar, arXiv:2511.17034), re-proved here
by Bailey's lemma and the Jacobi triple product: for every `k ≥ 0`,

  `Σ_{M_1 ≥ … ≥ M_k ≥ 0} ∏_{i=1}^k q^{M_i(M_i+1)/2} / (q;q)_{M_i - M_{i+1}}
      = (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ / ( (q;q²)_∞ (q;q)_∞ )`. -/
theorem melzer_identity (k : ℕ) :
    melzerSum (1 : ℚ⟦X⟧) k
      = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
        * (qPochInf 1 2 * qPochInf 1 1)⁻¹ := by
  rw [← krtw_sum_eq_melzer_sum]
  exact krtw_product_field k

/-- `qPochInf a d` is the infinite product `∏_{i ≥ 0} (1 - q^{a+di})` (`a, d ≥ 1`): its coefficient
of `q^N` is the coefficient of `q^N` of every finite product `∏_{i ≤ M} (1 - q^{a+di})` with
`M ≥ N`. -/
theorem coeff_qPochInf {K : Type*} [Field K] (a d : ℕ) (ha : 1 ≤ a) (hd : 1 ≤ d) (N M : ℕ)
    (h : N ≤ M) :
    coeff N (qPochInf a d : K⟦X⟧)
      = coeff N (∏ i ∈ range (M + 1), (1 - X ^ (a + d * i) : K⟦X⟧)) := by
  exact coeff_qPochInf' a d ha hd N M h

end KRTW
