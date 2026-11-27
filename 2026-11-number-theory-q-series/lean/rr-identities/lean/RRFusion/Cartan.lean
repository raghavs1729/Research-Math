/-
Finite certificates: the Cartan matrices of E8 and E6 and their inverses.

Labelling (as in `proofs/E8-E6.md`, 1-based there, 0-based `Fin` indices here):

    E8:  1 — 2 — 3 — 4 — 5 — 6 — 7          E6:  1 — 2 — 3 — 4 — 5
                         |                                |
                         8                                6
-/
import Mathlib

namespace RRFusion

/-- Cartan matrix of E8: chain 1-2-3-4-5-6-7, node 8 attached to node 5. -/
def C8 : Matrix (Fin 8) (Fin 8) ℤ :=
  !![ 2, -1,  0,  0,  0,  0,  0,  0;
     -1,  2, -1,  0,  0,  0,  0,  0;
      0, -1,  2, -1,  0,  0,  0,  0;
      0,  0, -1,  2, -1,  0,  0,  0;
      0,  0,  0, -1,  2, -1,  0, -1;
      0,  0,  0,  0, -1,  2, -1,  0;
      0,  0,  0,  0,  0, -1,  2,  0;
      0,  0,  0,  0, -1,  0,  0,  2]

/-- The matrix printed in §2 of the document (claimed to be `C(E8)⁻¹`). -/
def A8 : Matrix (Fin 8) (Fin 8) ℤ :=
  !![2,  3,  4,  5,  6,  4,  2,  3;
     3,  6,  8, 10, 12,  8,  4,  6;
     4,  8, 12, 15, 18, 12,  6,  9;
     5, 10, 15, 20, 24, 16,  8, 12;
     6, 12, 18, 24, 30, 20, 10, 15;
     4,  8, 12, 16, 20, 14,  7, 10;
     2,  4,  6,  8, 10,  7,  4,  5;
     3,  6,  9, 12, 15, 10,  5,  8]

/-- `C8` is the matrix "2 on the diagonal, −1 for adjacent nodes, 0 otherwise" of the E8 diagram. -/
def adjE8 (i j : Fin 8) : Bool :=
  (i.val + 1 = j.val ∧ j.val ≤ 6) || (j.val + 1 = i.val ∧ i.val ≤ 6) ||
    (i.val = 4 ∧ j.val = 7) || (i.val = 7 ∧ j.val = 4)

theorem C8_eq_diagram :
    C8 = Matrix.of fun i j => if i = j then 2 else if adjE8 i j then -1 else 0 := by
  decide

/-- `C8` is Mathlib's `CartanMatrix.E₈` (Bourbaki labelling) after relabelling the nodes. -/
theorem C8_eq_mathlib :
    C8 = Matrix.of fun i j =>
      CartanMatrix.E₈ (![7, 6, 5, 4, 3, 2, 0, 1] i) (![7, 6, 5, 4, 3, 2, 0, 1] j) := by
  decide

theorem C8_mul_A8 : C8 * A8 = 1 := by decide

theorem A8_mul_C8 : A8 * C8 = 1 := by decide

theorem A8_symm : A8.transpose = A8 := by decide

/-- Cartan matrix of E6 (over ℚ): chain 1-2-3-4-5, node 6 attached to node 3. -/
def C6 : Matrix (Fin 6) (Fin 6) ℚ :=
  !![ 2, -1,  0,  0,  0,  0;
     -1,  2, -1,  0,  0,  0;
      0, -1,  2, -1,  0, -1;
      0,  0, -1,  2, -1,  0;
      0,  0,  0, -1,  2,  0;
      0,  0, -1,  0,  0,  2]

/-- The matrix printed in §4 of the document (claimed to be `C(E6)⁻¹`). -/
def A6 : Matrix (Fin 6) (Fin 6) ℚ :=
  !![4/3,  5/3, 2,  4/3, 2/3, 1;
     5/3, 10/3, 4,  8/3, 4/3, 2;
       2,    4, 6,    4,   2, 3;
     4/3,  8/3, 4, 10/3, 5/3, 2;
     2/3,  4/3, 2,  5/3, 4/3, 1;
       1,    2, 3,    2,   1, 2]

def adjE6 (i j : Fin 6) : Bool :=
  (i.val + 1 = j.val ∧ j.val ≤ 4) || (j.val + 1 = i.val ∧ i.val ≤ 4) ||
    (i.val = 2 ∧ j.val = 5) || (i.val = 5 ∧ j.val = 2)

theorem C6_eq_diagram :
    C6 = Matrix.of fun i j => if i = j then 2 else if adjE6 i j then -1 else 0 := by
  decide

/-- `C6` is Mathlib's `CartanMatrix.E₆` (Bourbaki labelling) after relabelling the nodes. -/
theorem C6_eq_mathlib :
    C6 = Matrix.of fun i j =>
      ((CartanMatrix.E₆ (![0, 2, 3, 4, 5, 1] i) (![0, 2, 3, 4, 5, 1] j) : ℤ) : ℚ) := by
  decide

theorem C6_mul_A6 : C6 * A6 = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [C6, A6, Matrix.mul_apply, Fin.sum_univ_six] <;> norm_num

theorem A6_mul_C6 : A6 * C6 = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [C6, A6, Matrix.mul_apply, Fin.sum_univ_six] <;> norm_num

theorem A6_symm : A6.transpose = A6 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [A6]

end RRFusion
