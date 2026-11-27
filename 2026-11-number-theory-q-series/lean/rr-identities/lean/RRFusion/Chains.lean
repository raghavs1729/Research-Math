/-
The fusion chains of §2 (E8) and §4 (E6) of `proofs/E8-E6.md`:
iterating `fuse` on the Andrews–Gordon matrix `2 min(i,j)` (rank 4) gives `C(E8)⁻¹`, and on
`(4/3 2/3 0; 2/3 4/3 0; 0 0 2)` gives `C(E6)⁻¹`, after the stated relabelling of the nodes.
-/
import RRFusion.Cartan
import RRFusion.Fusion

namespace RRFusion

/-! ### E8 -/

/-- `A⁽⁰⁾ = (2 min(i,j))_{1 ≤ i,j ≤ 4}` (0-based indices here).
The four variables are named after the E8 nodes 1, 7, 2, 8. -/
def AG4 : Fin 4 → Fin 4 → ℚ := fun i j => 2 * ((min i.val j.val + 1 : ℕ) : ℚ)

/-- After `node 4 := node 7 + node 8`. -/
def E8step1 : Option (Fin 4) → Option (Fin 4) → ℚ := fuse AG4 1 3
/-- After `node 6 := node 1 + node 8`. -/
def E8step2 := fuse E8step1 (some 0) (some 3)
/-- After `node 5 := node 1 + node 4`. -/
def E8step3 := fuse E8step2 (some (some 0)) (some none)
/-- After `node 3 := node 1 + node 2`. -/
def E8step4 := fuse E8step3 (some (some (some 0))) (some (some (some 2)))

/-- Where the E8 node `i + 1` sits in the fused index set `Option⁴ (Fin 4)`. -/
def nodeE8 : Fin 8 → Option (Option (Option (Option (Fin 4)))) :=
  ![some (some (some (some 0))),   -- node 1 = variable m₁
    some (some (some (some 2))),   -- node 2 = variable m₃
    none,                          -- node 3, created 4th
    some (some (some none)),       -- node 4, created 1st
    some none,                     -- node 5, created 3rd
    some (some none),              -- node 6, created 2nd
    some (some (some (some 1))),   -- node 7 = variable m₂
    some (some (some (some 3)))]   -- node 8 = variable m₄

theorem nodeE8_bijective : Function.Bijective nodeE8 := by decide

/-- **E8 chain certificate.**  The 4-fold fusion of `2 min(i,j)` is the matrix `A8`
(which is `C(E8)⁻¹` by `C8_mul_A8`, `A8_mul_C8`). -/
theorem E8_chain : ∀ i j : Fin 8, E8step4 (nodeE8 i) (nodeE8 j) = (A8 i j : ℚ) := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    simp [E8step4, E8step3, E8step2, E8step1, fuse, nodeE8, AG4, A8] <;> norm_num

/-- The linear terms: the 4-fold fusion of `b = (b₁,b₂,b₃,b₄)` is the vector `B(b)` of Theorem 1. -/
theorem E8_chain_linear (b : Fin 4 → ℚ) : ∀ i : Fin 8,
    fuseB (fuseB (fuseB (fuseB b 1 3) (some 0) (some 3)) (some (some 0)) (some none))
        (some (some (some 0))) (some (some (some 2))) (nodeE8 i) =
      ![b 0, b 2, b 0 + b 2, b 1 + b 3, b 0 + b 1 + b 3, b 0 + b 3, b 1, b 3] i := by
  intro i
  fin_cases i <;> simp [fuseB, nodeE8, add_assoc]

/-! ### E6 -/

/-- The rank 3 start `(4/3 2/3 0; 2/3 4/3 0; 0 0 2)`; variables named after the E6 nodes 1, 5, 6. -/
def E6start : Fin 3 → Fin 3 → ℚ :=
  ![![4/3, 2/3, 0], ![2/3, 4/3, 0], ![0, 0, 2]]

/-- After `node 2 := node 5 + node 6`. -/
def E6step1 : Option (Fin 3) → Option (Fin 3) → ℚ := fuse E6start 1 2
/-- After `node 4 := node 1 + node 6`. -/
def E6step2 := fuse E6step1 (some 0) (some 2)
/-- After `node 3 := node 1 + node 2`. -/
def E6step3 := fuse E6step2 (some (some 0)) (some none)

/-- Where the E6 node `i + 1` sits in the fused index set `Option³ (Fin 3)`. -/
def nodeE6 : Fin 6 → Option (Option (Option (Fin 3))) :=
  ![some (some (some 0)),   -- node 1
    some (some none),       -- node 2, created 1st
    none,                   -- node 3, created 3rd
    some none,              -- node 4, created 2nd
    some (some (some 1)),   -- node 5
    some (some (some 2))]   -- node 6

theorem nodeE6_bijective : Function.Bijective nodeE6 := by decide

/-- **E6 chain certificate.**  The 3-fold fusion of `E6start` is the matrix `A6`
(which is `C(E6)⁻¹` by `C6_mul_A6`, `A6_mul_C6`). -/
theorem E6_chain : ∀ i j : Fin 6, E6step3 (nodeE6 i) (nodeE6 j) = A6 i j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    simp [E6step3, E6step2, E6step1, fuse, nodeE6, E6start, A6] <;> norm_num

/-- The linear terms of Theorem 2: `B = (b₁, b₂+b₃, b₁+b₂+b₃, b₁+b₃, b₂, b₃)`. -/
theorem E6_chain_linear (b : Fin 3 → ℚ) : ∀ i : Fin 6,
    fuseB (fuseB (fuseB b 1 2) (some 0) (some 2)) (some (some 0)) (some none) (nodeE6 i) =
      ![b 0, b 1 + b 2, b 0 + b 1 + b 2, b 0 + b 2, b 1, b 2] i := by
  intro i
  fin_cases i <;> simp [fuseB, nodeE6, add_assoc]

end RRFusion
