/-
Nahm/E7.lean — §1 of `proofs/E7.md`: the Nahm sum of `C(E7)⁻¹` is the rank 4 sum
`Σ_{a, n_1, n_2, n_3} q^{¾a² + a(N_2+N_3) + N_1²+N_2²+N_3²} / ((q)_a (q)_{n_1} (q)_{n_2} (q)_{n_3})`
(three fusions, Lemma 2 of `E8-E6.md`), everything in powers of `X = q^{1/4}`.

In the paper's labelling (chain 1—2—3—4—5—6, node 7 on node 4) the rank 4 variables `(a, n_1, n_2, n_3)`
are the nodes `(1, 6, 2, 5)` and the fusions are `7 := 1 + 6`, `4 := 5 + 6`, `3 := 1 + 2`.  The
statements use Mathlib's `CartanMatrix.E₇` (Bourbaki labelling); paper node `p` (1-based) is Bourbaki
node `(7, 6, 5, 4, 3, 1, 2)_p`.
-/
import Nahm.Fusion
import Nahm.E8
import Nahm.E7Lemmas

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

/-- `C(E7)⁻¹` for Mathlib's `CartanMatrix.E₇`, times 2 (an integer matrix). -/
def A7B2 : Matrix (Fin 7) (Fin 7) ℤ :=
  !![4, 4,  6,  8,  6, 4, 2;
     4, 7,  8, 12,  9, 6, 3;
     6, 8, 12, 16, 12, 8, 4;
     8, 12, 16, 24, 18, 12, 6;
     6, 9, 12, 18, 15, 10, 5;
     4, 6,  8, 12, 10, 8, 4;
     2, 3,  4,  6,  5, 4, 3]

theorem E7_mul_A7B2 : CartanMatrix.E₇ * A7B2 = 2 := by decide

/-- `C(E7)⁻¹` (Bourbaki labelling) over `ℚ`. -/
def A7Q : Matrix (Fin 7) (Fin 7) ℚ := Matrix.of fun i j => (A7B2 i j : ℚ) / 2

theorem E7_mul_A7Q : (Matrix.map CartanMatrix.E₇ (Int.cast : ℤ → ℚ)) * A7Q = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_seven, CartanMatrix.E₇, A7Q, A7B2] <;> norm_num

theorem invCartanE7_eq : invCartanE7 = A7Q := by
  rw [invCartanE7]
  exact Matrix.inv_eq_right_inv E7_mul_A7Q

/-! ### The fusion chain -/

/-- The rank 4 exponent (in powers of `X = q^{1/4}`): `3a² + 4(a(N_2+N_3) + N_1²+N_2²+N_3²)`
for `v = (a, n_1, n_2, n_3)`. -/
def E7e0 (v : Fin 4 → ℕ) : ℕ := 3 * v 0 ^ 2 + 4 * qexp ![0, v 0, v 0] (Fin.tail v)
def E7e1 : (Option (Fin 4) → ℕ) → ℕ := fuseE 4 E7e0 0 1
def E7e2 : (Option (Option (Fin 4)) → ℕ) → ℕ := fuseE 4 E7e1 (some 3) (some 1)
def E7e3 : (Option (Option (Option (Fin 4))) → ℕ) → ℕ := fuseE 4 E7e2 (some (some 0)) (some (some 2))

abbrev I3 := Option (Option (Option (Fin 4)))

/-- The position (Bourbaki, 0-based) of each node of the fused index set. -/
def toFin7 : I3 → Fin 7
  | none => 4
  | some none => 3
  | some (some none) => 1
  | some (some (some j)) => ![6, 0, 5, 2] j

theorem toFin7_bijective : Function.Bijective toFin7 := by decide

lemma psum3_0 (m : Fin 3 → ℕ) : psum m 0 = m 0 + m 1 + m 2 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_three]; simp
lemma psum3_1 (m : Fin 3 → ℕ) : psum m 1 = m 1 + m 2 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_three]; simp
lemma psum3_2 (m : Fin 3 → ℕ) : psum m 2 = m 2 := by
  simp only [psum, Finset.sum_filter, Fin.sum_univ_three]; simp

/-- The exponent identity behind the chain: the threefold fusion of the rank 4 exponent is
`4 · ½ nᵀ C(E7)⁻¹ n`. -/
theorem E7e3_eq (n : Fin 7 → ℕ) : ((E7e3 (fun j => n (toFin7 j)) : ℕ) : ℚ) = 4 * e7Form n := by
  simp only [E7e3, E7e2, E7e1, E7e0, fuseE, qexp, Fin.sum_univ_three, psum3_0, psum3_1, psum3_2]
  simp only [shift, toFin7, Fin.tail]
  simp
  rw [e7Form, invCartanE7_eq]
  simp [Matrix.mulVec, dotProduct, Fin.sum_univ_seven, A7Q, A7B2]
  ring

end Nahm
