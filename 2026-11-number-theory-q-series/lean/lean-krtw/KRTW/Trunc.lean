/-
KRTW/Trunc.lean — the summands with an index `> N` do not contribute to the coefficient of `q^N`
(this justifies the finite sums in the definitions of `krtwSum`, `melzerSum`).
-/
import KRTW.Defs

namespace KRTW

open PowerSeries Finset

variable {K : Type*} [Field K]

lemma X_pow_dvd_krtwTerm (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) (i : Fin k) :
    (X : K⟦X⟧) ^ ((r i + s i) ^ 2 + s i ^ 2 + s i) ∣ krtwTerm x r s := by
  refine Dvd.dvd.trans ?_ (Finset.dvd_prod_of_mem _ (Finset.mem_univ i))
  exact Dvd.dvd.mul_right (Dvd.dvd.mul_right (Dvd.dvd.mul_left (dvd_refl _) _) _) _

theorem coeff_krtwTerm_eq_zero' (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) (N : ℕ)
    (h : r ∉ box k N ∨ s ∉ box k N) : coeff N (krtwTerm x r s) = 0 := by
  have key : ∃ i, N + 1 ≤ (r i + s i) ^ 2 + s i ^ 2 + s i := by
    rcases h with h | h
    · simp only [box, Fintype.mem_piFinset, mem_range, not_forall, not_lt] at h
      obtain ⟨i, hi⟩ := h
      refine ⟨i, ?_⟩
      have : r i + s i ≤ (r i + s i) ^ 2 := Nat.le_self_pow (by norm_num) _
      omega
    · simp only [box, Fintype.mem_piFinset, mem_range, not_forall, not_lt] at h
      obtain ⟨i, hi⟩ := h
      exact ⟨i, by omega⟩
  obtain ⟨i, hi⟩ := key
  exact X_pow_dvd_iff.mp ((pow_dvd_pow X hi).trans (X_pow_dvd_krtwTerm x r s i)) N
    (Nat.lt_succ_self N)

theorem coeff_melzerTerm_eq_zero' (x : K⟦X⟧) {k : ℕ} (M : Fin k → ℕ) (N : ℕ)
    (h : M ∉ box k N) : coeff N (melzerTerm x M) = 0 := by
  simp only [box, Fintype.mem_piFinset, mem_range, not_forall, not_lt] at h
  obtain ⟨i, hi⟩ := h
  have h1 : (X : K⟦X⟧) ^ (M i * (M i + 1) / 2) ∣ melzerTerm x M := by
    refine Dvd.dvd.trans ?_ (Finset.dvd_prod_of_mem _ (Finset.mem_univ i))
    exact Dvd.dvd.mul_right (Dvd.dvd.mul_left (dvd_refl _) _) _
  have h2 : N + 1 ≤ M i * (M i + 1) / 2 := by
    have : M i ≤ M i * (M i + 1) / 2 := by
      rw [Nat.le_div_iff_mul_le (by norm_num)]
      nlinarith
    omega
  exact X_pow_dvd_iff.mp ((pow_dvd_pow X h2).trans h1) N (Nat.lt_succ_self N)

end KRTW
