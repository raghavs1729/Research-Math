/-
Nahm/Cauchy.lean — `(Q;Q)_∞` for a base `Q` with zero constant term, and Cauchy's identity
`Σ_{m ≥ 0} Q^{m(m+c)} / ((Q)_m (Q)_{m+c}) = 1/(Q;Q)_∞` (`c ≥ 0`), from the Durfee rectangle identity
(`RRFusion.durfee_inv`) by letting the size of the rectangle tend to infinity.
-/
import Nahm.Lattice
import RRFusion.Durfee
import KRTW.Summation

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K]

/-- `(Q;Q)_∞`, defined by its coefficients: the coefficient of `X^N` is that of `(Q;Q)_N`. -/
noncomputable def pinf (Q : K⟦X⟧) : K⟦X⟧ := PowerSeries.mk fun N => coeff N (poch Q N)

section base

variable (Q : K⟦X⟧) [Nil Q]

omit [Nil Q] in
lemma poch_add (N r : ℕ) : poch Q (N + r) = poch Q N * ∏ i ∈ range r, (1 - Q ^ (N + i + 1)) := by
  simp only [poch]; rw [prod_range_add]

lemma poch_trunc {N n : ℕ} (h : N ≤ n) : (X : K⟦X⟧) ^ (N + 1) ∣ poch Q n - poch Q N := by
  obtain ⟨r, rfl⟩ := Nat.exists_eq_add_of_le h
  rw [poch_add, show poch Q N * ∏ i ∈ range r, (1 - Q ^ (N + i + 1)) - poch Q N
      = -(poch Q N) * (1 - ∏ i ∈ range r, (1 - Q ^ (N + i + 1))) by ring]
  refine Dvd.dvd.mul_left (X_pow_dvd_one_sub_prod _ _ _ (fun i _ => ?_)) _
  rw [sub_sub_cancel]
  exact X_pow_dvd_pow_of_le Q (by omega)

lemma pinf_trunc (n : ℕ) : (X : K⟦X⟧) ^ (n + 1) ∣ pinf Q - poch Q n := by
  rw [X_pow_dvd_iff]
  intro m hm
  rw [map_sub, pinf, coeff_mk, sub_eq_zero]
  have := X_pow_dvd_iff.mp (poch_trunc Q (show m ≤ n by omega)) m (Nat.lt_succ_self m)
  rw [map_sub, sub_eq_zero] at this
  exact this.symm

lemma constantCoeff_pinf : constantCoeff (pinf Q) = 1 := by
  rw [← coeff_zero_eq_constantCoeff_apply, pinf, coeff_mk, coeff_zero_eq_constantCoeff_apply,
    constantCoeff_poch]

lemma pinf_mul_inv : pinf Q * (pinf Q)⁻¹ = 1 :=
  PowerSeries.mul_inv_cancel _ (by rw [constantCoeff_pinf]; exact one_ne_zero)

lemma pinf_inv_ne_zero : (pinf Q)⁻¹ ≠ 0 := by
  intro h
  have := pinf_mul_inv Q
  rw [h, mul_zero] at this
  exact zero_ne_one this

lemma ipinf_trunc (n : ℕ) : (X : K⟦X⟧) ^ (n + 1) ∣ (pinf Q)⁻¹ - iq Q n := by
  have e : (pinf Q)⁻¹ - iq Q n = (pinf Q)⁻¹ * iq Q n * (poch Q n - pinf Q) := by
    have h1 := pinf_mul_inv Q
    have h2 := poch_mul_iq Q n
    linear_combination (-(pinf Q)⁻¹) * h2 + (iq Q n) * h1
  rw [e]
  exact Dvd.dvd.mul_left ((dvd_sub_comm).mp (pinf_trunc Q n)) _

lemma ipinf_trunc_of_le {n N : ℕ} (h : n ≤ N) : (X : K⟦X⟧) ^ (n + 1) ∣ (pinf Q)⁻¹ - iq Q N :=
  (pow_dvd_pow X (by omega)).trans (ipinf_trunc Q N)

/-- The Durfee rectangle identity in `K⟦X⟧` (from `RRFusion.durfee_inv`). -/
lemma durfee (a b : ℕ) :
    iq Q a * iq Q b = ∑ c ∈ range (min a b + 1),
      Q ^ ((a - c) * (b - c)) * (iq Q (a - c) * iq Q (b - c) * iq Q c) :=
  RRFusion.durfee_inv Q (iq Q) (fun n => iq_mul_poch Q n) a b

lemma X_pow_dvd_Q_pow (m e : ℕ) (h : m ≤ e) : (X : K⟦X⟧) ^ m ∣ Q ^ e := X_pow_dvd_pow_of_le Q h

omit [Nil Q] in
lemma le_mul_add_self (m c : ℕ) : m ≤ m * (m + c) := by
  rcases Nat.eq_zero_or_pos m with h | h
  · simp [h]
  · exact Nat.le_mul_of_pos_right m (by omega)

variable [TopologicalSpace K] [DiscreteTopology K]

lemma summable_cauchy (c : ℕ) : Summable fun m => Q ^ (m * (m + c)) * (iq Q m * iq Q (m + c)) := by
  apply summable_nat
  intro m
  exact Dvd.dvd.mul_right (X_pow_dvd_Q_pow Q _ _ (le_mul_add_self _ _)) _

/-- **Cauchy's identity** `Σ_m Q^{m(m+c)}/((Q)_m (Q)_{m+c}) = 1/(Q;Q)_∞`. -/
theorem cauchy (c : ℕ) : ∑' m, Q ^ (m * (m + c)) * (iq Q m * iq Q (m + c)) = (pinf Q)⁻¹ := by
  set C := ∑' m, Q ^ (m * (m + c)) * (iq Q m * iq Q (m + c)) with hC
  set P := (pinf Q)⁻¹ with hP
  have key : ∀ A : ℕ, (X : K⟦X⟧) ^ (A + 1) ∣ P * P - C * P := by
    intro A
    have hd := durfee Q A (A + c)
    rw [min_eq_left (Nat.le_add_right A c)] at hd
    rw [← sum_range_reflect] at hd
    have hd' : iq Q A * iq Q (A + c) = ∑ m ∈ range (A + 1),
        Q ^ (m * (m + c)) * (iq Q m * iq Q (m + c)) * iq Q (A - m) := by
      rw [hd]
      refine sum_congr rfl (fun m hm => ?_)
      have hm' : m ≤ A := Nat.lt_succ_iff.mp (mem_range.mp hm)
      rw [show A - (A + 1 - 1 - m) = m by omega, show A + c - (A + 1 - 1 - m) = m + c by omega,
        show A + 1 - 1 - m = A - m by omega]
      ring
    -- (a)
    have ha : (X : K⟦X⟧) ^ (A + 1) ∣ P * P - iq Q A * iq Q (A + c) :=
      KRTW.cong_mul (ipinf_trunc Q A) (ipinf_trunc_of_le Q (Nat.le_add_right A c))
    -- (b)
    have hb : (X : K⟦X⟧) ^ (A + 1) ∣ iq Q A * iq Q (A + c)
        - ∑ m ∈ range (A + 1), Q ^ (m * (m + c)) * (iq Q m * iq Q (m + c)) * P := by
      rw [hd', ← sum_sub_distrib]
      refine dvd_sum (fun m hm => ?_)
      have hm' : m ≤ A := Nat.lt_succ_iff.mp (mem_range.mp hm)
      rw [← mul_sub, show (X : K⟦X⟧) ^ (A + 1) = X ^ m * X ^ (A - m + 1) by
        rw [← pow_add]; congr 1; omega]
      exact mul_dvd_mul (Dvd.dvd.mul_right (X_pow_dvd_Q_pow Q _ _ (le_mul_add_self _ _)) _)
        ((dvd_sub_comm).mp (ipinf_trunc Q (A - m)))
    -- (c)
    have hc : (X : K⟦X⟧) ^ (A + 1) ∣ ∑ m ∈ range (A + 1), Q ^ (m * (m + c)) * (iq Q m * iq Q (m + c)) * P
        - C * P := by
      rw [← sum_mul, ← sub_mul]
      refine Dvd.dvd.mul_right ?_ _
      rw [dvd_sub_comm]
      apply X_pow_dvd_tsum_sub_sum (summable_cauchy Q c)
      intro m hm
      have : A + 1 ≤ m := by simpa using hm
      exact Dvd.dvd.mul_right (X_pow_dvd_Q_pow Q _ _ (le_trans this (le_mul_add_self _ _))) _
    exact dvd_sub_trans (dvd_sub_trans ha hb) hc
  have h0 : P * P - C * P = 0 := by
    have := eq_of_forall_X_pow_dvd_sub (f := P * P - C * P) (g := 0)
      (fun A => (pow_dvd_pow X (Nat.le_succ A)).trans (by simpa using key A))
    exact this
  rw [← sub_mul, mul_eq_zero] at h0
  rcases h0 with h0 | h0
  · exact (sub_eq_zero.mp h0).symm
  · exact absurd h0 (pinf_inv_ne_zero Q)

end base

end Nahm
