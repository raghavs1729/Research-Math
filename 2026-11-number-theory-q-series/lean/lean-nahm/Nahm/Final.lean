/-
Nahm/Final.lean — auxiliary facts for the statement layer: the finite sums in the coefficient
definitions lose nothing, the explicit inverse Cartan matrices, the general-field versions of the
headline theorems without topology hypotheses.
-/
import Nahm.E7Final

namespace Nahm

open PowerSeries Finset KRTW

variable {K : Type*} [Field K]

/-! ### The truncations in the coefficient definitions -/

lemma coeff_X_pow_mul_eq_zero (e N : ℕ) (f : K⟦X⟧) (h : N < e) : coeff N ((X : K⟦X⟧) ^ e * f) = 0 :=
  X_pow_dvd_iff.mp (Dvd.dvd.mul_right (dvd_refl _) _) N h

lemma exists_gt_of_not_mem_box {r N : ℕ} {n : Fin r → ℕ} (h : n ∉ box r N) : ∃ j, N < n j := by
  simp only [box, Fintype.mem_piFinset, mem_range, not_forall, not_lt] at h
  obtain ⟨j, hj⟩ := h
  exact ⟨j, by omega⟩

lemma le_agExp (k i : ℕ) (n : Fin (k - 1) → ℕ) (j : Fin (k - 1)) : n j ≤ agExp k i n := by
  have := le_trans (le_tot n j) (tot_le_qexp (agb i) n)
  rwa [qexp_agb] at this

lemma coeff_agTerm_eq_zero' (k i : ℕ) (n : Fin (k - 1) → ℕ) (N : ℕ) (h : n ∉ box (k - 1) N) :
    coeff N ((X : K⟦X⟧) ^ agExp k i n * ∏ j, (qPoch (n j))⁻¹) = 0 := by
  obtain ⟨j, hj⟩ := exists_gt_of_not_mem_box h
  exact coeff_X_pow_mul_eq_zero _ _ _ (lt_of_lt_of_le hj (le_agExp k i n j))

lemma e8Exp_eq_add (B n : Fin 8 → ℕ) : e8Exp B n = e8Exp 0 n + ∑ i, B i * n i := by
  have h1 := e8Form_eq_e8Exp' B n
  have h0 := e8Form_eq_e8Exp' 0 n
  have : e8Form B n = e8Form 0 n + ∑ i, (B i : ℚ) * n i := by
    simp only [e8Form, dotProduct, Pi.zero_apply, Nat.cast_zero, zero_mul, sum_const_zero, add_zero]
  have h2 : ((e8Exp B n : ℕ) : ℚ) = ((e8Exp 0 n + ∑ i, B i * n i : ℕ) : ℚ) := by
    rw [h1, this, ← h0]; push_cast; rfl
  exact_mod_cast h2

lemma le_e8Exp (B n : Fin 8 → ℕ) (i : Fin 8) : n i ≤ e8Exp B n := by
  rw [e8Exp_eq_add]
  have := le_e8Exp_Bb 0 n i
  rw [Bb_zero] at this
  omega

lemma coeff_e8Term_eq_zero' (B n : Fin 8 → ℕ) (N : ℕ) (h : n ∉ box 8 N) :
    coeff N ((X : K⟦X⟧) ^ e8Exp B n * ∏ i, (qPoch (n i))⁻¹) = 0 := by
  obtain ⟨j, hj⟩ := exists_gt_of_not_mem_box h
  exact coeff_X_pow_mul_eq_zero _ _ _ (lt_of_lt_of_le hj (le_e8Exp B n j))

lemma coeff_e7Term_eq_zero' (n : Fin 7 → ℕ) (N : ℕ) (h : n ∉ box 7 N) :
    coeff N ((X : K⟦X⟧) ^ e7Exp n * ∏ i, (qPoch4 (n i))⁻¹) = 0 := by
  obtain ⟨j, hj⟩ := exists_gt_of_not_mem_box h
  exact coeff_X_pow_mul_eq_zero _ _ _ (lt_of_lt_of_le hj (le_e7Exp n j))

lemma coeff_e7ThetaTerm_eq_zero' (k : ℤ) (N : ℕ) (h : k ∉ Icc (-(N : ℤ)) N) :
    coeff N ((X : K⟦X⟧) ^ (3 * k.natAbs ^ 2) + X ^ (7 * k.natAbs ^ 2)) = 0 := by
  have hk : N < k.natAbs := by simp only [mem_Icc, not_and_or, not_le] at h; omega
  rw [map_add, coeff_X_pow, coeff_X_pow, if_neg (by nlinarith), if_neg (by nlinarith), add_zero]

/-! ### The explicit inverse Cartan matrices -/

theorem invCartanE8_explicit : invCartanE8 =
    !![(4 : ℚ), 5, 7, 10, 8, 6, 4, 2;
       5, 8, 10, 15, 12, 9, 6, 3;
       7, 10, 14, 20, 16, 12, 8, 4;
       10, 15, 20, 30, 24, 18, 12, 6;
       8, 12, 16, 24, 20, 15, 10, 5;
       6, 9, 12, 18, 15, 12, 8, 4;
       4, 6, 8, 12, 10, 8, 6, 3;
       2, 3, 4, 6, 5, 4, 3, 2] := by
  rw [invCartanE8_eq]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [A8B]

theorem invCartanE7_explicit : invCartanE7 =
    !![(2 : ℚ), 2, 3, 4, 3, 2, 1;
       2, 7/2, 4, 6, 9/2, 3, 3/2;
       3, 4, 6, 8, 6, 4, 2;
       4, 6, 8, 12, 9, 6, 3;
       3, 9/2, 6, 9, 15/2, 5, 5/2;
       2, 3, 4, 6, 5, 4, 2;
       1, 3/2, 2, 3, 5/2, 2, 3/2] := by
  rw [invCartanE7_eq]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [A7Q, A7B2] <;> norm_num

/-! ### The headline theorems over an arbitrary field -/

theorem andrews_gordon_top' (k i : ℕ) (hi : 1 ≤ i) (hik : i ≤ k) :
    (agSum k i : K⟦X⟧) = agProduct k i := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact andrews_gordon_top k i hi hik

theorem e8_vacuum_top' :
    (e8Sum 0 : K⟦X⟧) = (∏ a ∈ ({1, 2, 3, 4, 7, 8, 9, 10} : Finset ℕ), qPochInf a 11)⁻¹ := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact e8_vacuum_top

theorem e8_companion_top' (i : ℕ) (hi : 1 ≤ i) (hi5 : i ≤ 5) :
    (e8Sum (Bb (agb i)) : K⟦X⟧) = agProduct 5 i := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact e8_companion_top i hi hi5

theorem e7_identity_top' : (e7Sum : K⟦X⟧) = (qPochInf 4 4)⁻¹ * e7ThetaOneSided := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact e7_identity_top

theorem e7_identity_bilateral_top' : 2 * (e7Sum : K⟦X⟧) = (qPochInf 4 4)⁻¹ * e7Theta := by
  let _ : TopologicalSpace K := ⊥
  have _ : DiscreteTopology K := ⟨rfl⟩
  exact e7_identity_bilateral_top

end Nahm
