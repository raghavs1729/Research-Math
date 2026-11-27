/-
Nahm/E8.lean — Theorem 1 and Corollary 1 of `proofs/E8-E6.md`: the Nahm sums of `C(E8)⁻¹`.

`C(E8)⁻¹` is an iterated fusion (Lemma 2, `Nahm.fusion`) of the rank 4 Andrews–Gordon matrix
`(2 min(i,j))`.  In the paper's labelling (chain 1—2—3—4—5—6—7, node 8 on node 5) the four
Andrews–Gordon variables are the nodes 1, 7, 2, 8 and the fusions are
`4 := 7 + 8`, `6 := 1 + 8`, `5 := 1 + 4`, `3 := 1 + 2`.  The statements use Mathlib's
`CartanMatrix.E₈` (Bourbaki labelling); paper node `p` is Bourbaki node `σ p` with
`σ = (8, 7, 6, 5, 4, 3, 1, 2)` (1-based), i.e. `![7, 6, 5, 4, 3, 2, 0, 1]` 0-based
(`RRFusion.C8_eq_mathlib`).
-/
import Nahm.Fusion
import Nahm.AG
import RRFusion.Cartan

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

/-! ### The inverse Cartan matrix in Bourbaki labelling -/

/-- `C(E8)⁻¹` for Mathlib's `CartanMatrix.E₈`. -/
def A8B : Matrix (Fin 8) (Fin 8) ℤ :=
  !![ 4,  5,  7, 10,  8,  6,  4, 2;
      5,  8, 10, 15, 12,  9,  6, 3;
      7, 10, 14, 20, 16, 12,  8, 4;
     10, 15, 20, 30, 24, 18, 12, 6;
      8, 12, 16, 24, 20, 15, 10, 5;
      6,  9, 12, 18, 15, 12,  8, 4;
      4,  6,  8, 12, 10,  8,  6, 3;
      2,  3,  4,  6,  5,  4,  3, 2]

theorem E8_mul_A8B : CartanMatrix.E₈ * A8B = 1 := by decide

theorem A8B_mul_E8 : A8B * CartanMatrix.E₈ = 1 := by decide

theorem invCartanE8_eq : invCartanE8 = A8B.map (Int.cast : ℤ → ℚ) := by
  rw [invCartanE8]
  apply Matrix.inv_eq_right_inv
  have h := congrArg (fun M => M.map (Int.castRingHom ℚ)) E8_mul_A8B
  simp only [Matrix.map_mul] at h
  simpa using h

/-! ### The fusion chain -/

/-- The index set after four fusions. -/
abbrev I4 := Option (Option (Option (Option (Fin 4))))

/-- The four successive exponents. -/
def E8e0 (b : Fin 4 → ℕ) : (Fin 4 → ℕ) → ℕ := qexp b
def E8e1 (b : Fin 4 → ℕ) : (Option (Fin 4) → ℕ) → ℕ := fuseE 1 (E8e0 b) 1 3
def E8e2 (b : Fin 4 → ℕ) : (Option (Option (Fin 4)) → ℕ) → ℕ := fuseE 1 (E8e1 b) (some 0) (some 3)
def E8e3 (b : Fin 4 → ℕ) : (Option (Option (Option (Fin 4))) → ℕ) → ℕ :=
  fuseE 1 (E8e2 b) (some (some 0)) (some none)
def E8e4 (b : Fin 4 → ℕ) : (I4 → ℕ) → ℕ :=
  fuseE 1 (E8e3 b) (some (some (some 0))) (some (some (some 2)))

/-- The position (Bourbaki, 0-based) of each node of the fused index set. -/
def toFin8 : I4 → Fin 8
  | none => 5
  | some none => 3
  | some (some none) => 2
  | some (some (some none)) => 4
  | some (some (some (some j))) => ![7, 0, 6, 1] j

theorem toFin8_bijective : Function.Bijective toFin8 := by decide

/-- The linear term `B(b')` of Theorem 1, in Bourbaki labelling, where `b'` is the linear term on the
variables `m_j` (`b'_j = b_0 + ⋯ + b_j`) of the sum `Σ q^{Σ_l (M_l² + b_l M_l)}/∏(q)_{m_l}`. -/
def Bb (b : Fin 4 → ℕ) : Fin 8 → ℕ :=
  let b' : Fin 4 → ℕ := ![b 0, b 0 + b 1, b 0 + b 1 + b 2, b 0 + b 1 + b 2 + b 3]
  ![b' 1, b' 3, b' 0 + b' 3, b' 0 + b' 1 + b' 3, b' 1 + b' 3, b' 0 + b' 2, b' 2, b' 0]

lemma psum4_0 (m : Fin 4 → ℕ) : psum m 0 = m 0 + m 1 + m 2 + m 3 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_four]; simp
lemma psum4_1 (m : Fin 4 → ℕ) : psum m 1 = m 1 + m 2 + m 3 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_four]; simp
lemma psum4_2 (m : Fin 4 → ℕ) : psum m 2 = m 2 + m 3 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_four]; simp
lemma psum4_3 (m : Fin 4 → ℕ) : psum m 3 = m 3 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_four]; simp

/-- The exponent identity behind the chain: the fourfold fusion of `M_1² + ⋯ + M_4² + b·M` is
`½ nᵀ C(E8)⁻¹ n + B(b)·n`. -/
theorem E8e4_eq (b : Fin 4 → ℕ) (n : Fin 8 → ℕ) :
    ((E8e4 b (fun j => n (toFin8 j)) : ℕ) : ℚ) = e8Form (Bb b) n := by
  simp only [E8e4, E8e3, E8e2, E8e1, E8e0, fuseE, qexp, Fin.sum_univ_four, psum4_0, psum4_1, psum4_2,
    psum4_3]
  simp only [shift, toFin8]
  simp
  rw [e8Form, invCartanE8_eq]
  simp [Matrix.mulVec, dotProduct, Fin.sum_univ_eight, A8B, Bb]
  ring

lemma le_E8e4 (b : Fin 4 → ℕ) : ∀ n i, n i ≤ E8e4 b n := by
  have h0 : ∀ n i, n i ≤ E8e0 b n := fun n i => le_trans (le_tot n i) (tot_le_qexp b n)
  exact le_fuseE 1 _ (le_fuseE 1 _ (le_fuseE 1 _ (le_fuseE 1 _ h0 _ _) _ _) _ _) _ _

/-- The Bourbaki-labelled tuple of a tuple on the fused index set. -/
noncomputable def e8Equiv : I4 ≃ Fin 8 := Equiv.ofBijective toFin8 toFin8_bijective

lemma e8Exp_Bb (b : Fin 4 → ℕ) (n : Fin 8 → ℕ) : e8Exp (Bb b) n = E8e4 b (fun j => n (toFin8 j)) := by
  rw [e8Exp, ← E8e4_eq, Nat.floor_natCast]

lemma le_e8Exp_Bb (b : Fin 4 → ℕ) (n : Fin 8 → ℕ) (i : Fin 8) : n i ≤ e8Exp (Bb b) n := by
  rw [e8Exp_Bb]
  obtain ⟨j, rfl⟩ := toFin8_bijective.2 i
  exact le_E8e4 b (fun j => n (toFin8 j)) j

variable {K : Type*} [Field K]

lemma Bb_zero : Bb 0 = 0 := by decide

/-- **Integrality**: `½ nᵀ C(E8)⁻¹ n + B·n` is a natural number (namely `e8Exp B n`). -/
theorem e8Form_eq_e8Exp' (B n : Fin 8 → ℕ) : (e8Exp B n : ℚ) = e8Form B n := by
  have h0 := E8e4_eq 0 n
  rw [Bb_zero] at h0
  have h1 : e8Form B n = ((E8e4 0 (fun j => n (toFin8 j)) + ∑ i, B i * n i : ℕ) : ℚ) := by
    rw [Nat.cast_add, h0]
    simp only [e8Form, dotProduct, Pi.zero_apply, Nat.cast_zero, zero_mul, Finset.sum_const_zero, add_zero,
      Nat.cast_sum, Nat.cast_mul]
  rw [e8Exp, h1, Nat.floor_natCast]

variable [TopologicalSpace K] [DiscreteTopology K]

/-- A sum defined by its coefficients over boxes is the infinite sum. -/
lemma mk_box_eq_tsum {r : ℕ} (f : (Fin r → ℕ) → K⟦X⟧) (w : (Fin r → ℕ) → ℕ) (hw : ∀ n i, n i ≤ w n)
    (hf : ∀ n, (X : K⟦X⟧) ^ w n ∣ f n) :
    PowerSeries.mk (fun N => ∑ n ∈ box r N, coeff N (f n)) = ∑' n, f n := by
  have hsum : Summable f := by
    apply summable_of_dvd w _ hf
    intro N
    apply (box r N).finite_toSet.subset
    intro n hn
    simp only [Set.mem_ofPred_eq] at hn
    simp only [box, mem_coe, Fintype.mem_piFinset, mem_range]
    exact fun i => Nat.lt_succ_of_le ((hw n i).trans hn)
  ext N
  rw [coeff_mk, coeff_tsum hsum]
  symm
  apply tsum_eq_sum
  intro n hn
  simp only [box, Fintype.mem_piFinset, mem_range, not_forall, not_lt] at hn
  obtain ⟨i, hi⟩ := hn
  exact X_pow_dvd_iff.mp (hf n) N (by have := hw n i; omega)

/-- **Theorem 1 of `E8-E6.md`**: `F_{C(E8)⁻¹, B(b')} = Σ_{m ∈ ℕ⁴} q^{Σ_l (M_l² + b_l M_l)}/∏(q)_{m_l}`. -/
theorem e8_theorem1 (b : Fin 4 → ℕ) : (e8Sum (Bb b) : K⟦X⟧) = ∑' m, mterm X b m := by
  have hX : (X : K⟦X⟧) = X ^ 1 := (pow_one X).symm
  have h0 : ∀ n i, n i ≤ E8e0 b n := fun n i => le_trans (le_tot n i) (tot_le_qexp b n)
  have h1 := le_fuseE 1 _ h0 (1 : Fin 4) 3
  have h2 := le_fuseE 1 _ h1 (some (0 : Fin 4)) (some 3)
  have h3 := le_fuseE 1 _ h2 (some (some (0 : Fin 4))) (some none)
  -- the chain of fusions
  have f4 := fusion (Q := (X : K⟦X⟧)) 1 hX (E8e3 b) h3 (some (some (some 0))) (some (some (some 2))) (by decide)
  have f3 := fusion (Q := (X : K⟦X⟧)) 1 hX (E8e2 b) h2 (some (some 0)) (some none) (by decide)
  have f2 := fusion (Q := (X : K⟦X⟧)) 1 hX (E8e1 b) h1 (some 0) (some 3) (by decide)
  have f1 := fusion (Q := (X : K⟦X⟧)) 1 hX (E8e0 b) h0 1 3 (by decide)
  have hchain : ∑' n : I4 → ℕ, nterm (X : K⟦X⟧) (E8e4 b) n = ∑' m, mterm (X : K⟦X⟧) b m := by
    rw [E8e4, f4, E8e3, f3, E8e2, f2, E8e1, f1]
    rfl
  rw [← hchain, e8Sum]
  rw [mk_box_eq_tsum (fun n => X ^ e8Exp (Bb b) n * ∏ i, (qPoch (n i))⁻¹) (e8Exp (Bb b))
    (le_e8Exp_Bb b) (fun n => Dvd.dvd.mul_right (dvd_refl _) _)]
  let φ : (Fin 8 → ℕ) ≃ (I4 → ℕ) := Equiv.arrowCongr e8Equiv.symm (Equiv.refl ℕ)
  rw [← φ.tsum_eq]
  refine tsum_congr (fun n => ?_)
  simp only [nterm, φ, Equiv.arrowCongr_apply, Equiv.symm_symm, Equiv.coe_refl, Function.comp_apply,
    id]
  rw [e8Exp_Bb]
  congr 1
  exact (Fintype.prod_equiv e8Equiv (fun j => iq X (n (e8Equiv j))) (fun i => (qPoch (n i))⁻¹)
    (fun j => rfl)).symm

lemma Bb_agb5 : Bb (agb 5) = 0 := by decide
lemma Bb_agb4 : Bb (agb 4) = ![0, 1, 1, 1, 1, 0, 0, 0] := by decide
lemma Bb_agb3 : Bb (agb 3) = ![0, 2, 2, 2, 2, 1, 1, 0] := by decide
lemma Bb_agb2 : Bb (agb 2) = ![1, 3, 3, 4, 4, 2, 2, 0] := by decide
lemma Bb_agb1 : Bb (agb 1) = ![2, 4, 5, 7, 6, 4, 3, 1] := by decide

/-- Corollary 1 of `E8-E6.md`, all five cases: `F_{C(E8)⁻¹, B^{(i)}} = ∏_{n ≢ 0, ±i (mod 11)} (1-q^n)⁻¹`. -/
theorem e8_companion_top (i : ℕ) (hi : 1 ≤ i) (hi5 : i ≤ 5) :
    (e8Sum (Bb (agb i)) : K⟦X⟧) = agProduct 5 i := by
  rw [e8_theorem1]
  have h := agSum_eq_tsum (K := K) 5 i
  rw [← andrews_gordon_top 5 i hi hi5, h]

/-- The E8 vacuum identity `F_{C(E8)⁻¹, 0} = 1/(q, q², q³, q⁴, q⁷, q⁸, q⁹, q¹⁰; q¹¹)_∞`. -/
theorem e8_vacuum_top :
    (e8Sum 0 : K⟦X⟧) = (∏ a ∈ ({1, 2, 3, 4, 7, 8, 9, 10} : Finset ℕ), qPochInf a 11)⁻¹ := by
  rw [← Bb_agb5, e8_companion_top 5 (by norm_num) le_rfl, agProduct]
  congr 2

end Nahm
