/-
The fusion operation of Lemma 2 of `proofs/E8-E6.md` on the data `(A, B)` of a Nahm sum,
and the algebraic identity (1.2).

Matrices are functions `ι → ι → ℚ`; the enlarged index set is `Option ι`, the new species is `none`.
-/
import Mathlib

namespace RRFusion

open Finset

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- `Q_{A,B}(n) = ½ nᵀAn + Bᵀn`. -/
def Q (A : ι → ι → ℚ) (B : ι → ℚ) (n : ι → ℚ) : ℚ :=
  (1 / 2) * ∑ i, ∑ j, A i j * n i * n j + ∑ i, B i * n i

/-- Fusion of the species `a` and `b`: the matrix `A°` of Lemma 2 (new species = `none`). -/
def fuse (A : ι → ι → ℚ) (a b : ι) : Option ι → Option ι → ℚ
  | some i, some j => A i j + if (i = a ∧ j = b) ∨ (i = b ∧ j = a) then 1 else 0
  | some i, none => A i a + A i b
  | none, some j => A j a + A j b
  | none, none => A a a + A b b + 2 * A a b

/-- The linear term `B°` of Lemma 2. -/
def fuseB (B : ι → ℚ) (a b : ι) : Option ι → ℚ
  | some i => B i
  | none => B a + B b

/-- The vector `(m, c)` indexed by `Option ι`. -/
def extend (m : ι → ℚ) (c : ℚ) : Option ι → ℚ
  | some i => m i
  | none => c

/-- The vector `v = e_a + e_b`. -/
def pairVec (a b : ι) (i : ι) : ℚ := (if i = a then 1 else 0) + (if i = b then 1 else 0)

omit [Fintype ι] in
theorem fuse_symm {A : ι → ι → ℚ} (hA : ∀ i j, A i j = A j i) (a b : ι) :
    ∀ x y, fuse A a b x y = fuse A a b y x := by
  rintro (_ | i) (_ | j)
  · rfl
  · rfl
  · rfl
  · simp only [fuse, hA i j, and_comm, or_comm]

omit [Fintype ι] in
private lemma bump_eq {a b : ι} (hab : a ≠ b) (i j : ι) :
    (if (i = a ∧ j = b) ∨ (i = b ∧ j = a) then (1 : ℚ) else 0) =
      (if i = a then 1 else 0) * (if j = b then 1 else 0) +
        (if i = b then 1 else 0) * (if j = a then 1 else 0) := by
  by_cases h1 : i = a <;> by_cases h2 : j = b <;> by_cases h3 : i = b <;> by_cases h4 : j = a <;>
    simp_all

private lemma sum_mul_pairVec (a b : ι) (f : ι → ℚ) :
    ∑ i, f i * pairVec a b i = f a + f b := by
  simp [pairVec, mul_add, Finset.sum_add_distrib]

/-- **Identity (1.2).**  For symmetric `A`, `a ≠ b`, and `n = m + c (e_a + e_b)`:
`Q_{A°,B°}(m, c) = Q_{A,B}(n) + m_a m_b`. -/
theorem Q_fuse {A : ι → ι → ℚ} (hA : ∀ i j, A i j = A j i) (B : ι → ℚ) {a b : ι} (hab : a ≠ b)
    (m : ι → ℚ) (c : ℚ) :
    Q (fuse A a b) (fuseB B a b) (extend m c) =
      Q A B (fun i => m i + c * pairVec a b i) + m a * m b := by
  -- the three sums in which everything is expressed
  set S := ∑ i, ∑ j, A i j * m i * m j with hS
  set Ra := ∑ i, A i a * m i with hRa
  set Rb := ∑ i, A i b * m i with hRb
  have hT : ∀ x, ∑ j, A x j * m j = ∑ i, A i x * m i := fun x =>
    Finset.sum_congr rfl fun j _ => by rw [hA x j]
  -- left side
  have h1 : ∑ j, (A j a + A j b) * c * m j = c * (Ra + Rb) := by
    rw [hRa, hRb, ← Finset.sum_add_distrib, Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hbump : ∑ i, ∑ j, (if (i = a ∧ j = b) ∨ (i = b ∧ j = a) then (1 : ℚ) else 0) * m i * m j
      = 2 * (m a * m b) := by
    simp only [bump_eq hab, add_mul, Finset.sum_add_distrib]
    simp [Finset.sum_ite_eq']
    ring
  have h2 : ∑ i, ((A i a + A i b) * m i * c +
        ∑ j, (A i j + if (i = a ∧ j = b) ∨ (i = b ∧ j = a) then 1 else 0) * m i * m j)
      = c * (Ra + Rb) + S + 2 * (m a * m b) := by
    rw [← hbump, hS, hRa, hRb, mul_add, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib,
      ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have : ∑ j, (A i j + if (i = a ∧ j = b) ∨ (i = b ∧ j = a) then 1 else 0) * m i * m j
        = ∑ j, A i j * m i * m j
          + ∑ j, (if (i = a ∧ j = b) ∨ (i = b ∧ j = a) then (1 : ℚ) else 0) * m i * m j := by
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun j _ => by ring
    rw [this]
    ring
  -- right side
  have h3 : ∑ i, ∑ j, A i j * (m i + c * pairVec a b i) * (m j + c * pairVec a b j)
      = S + 2 * c * (Ra + Rb) + c ^ 2 * (A a a + A b b + 2 * A a b) := by
    have inner : ∀ i, ∑ j, A i j * (m i + c * pairVec a b i) * (m j + c * pairVec a b j)
        = (m i + c * pairVec a b i) * (∑ j, A i j * m j + c * (A i a + A i b)) := by
      intro i
      rw [← sum_mul_pairVec a b (fun j => A i j), Finset.mul_sum, ← Finset.sum_add_distrib,
        Finset.mul_sum]
      exact Finset.sum_congr rfl fun j _ => by ring
    simp only [inner]
    have e1 : ∑ i, (m i + c * pairVec a b i) * (∑ j, A i j * m j + c * (A i a + A i b))
        = ∑ i, ∑ j, A i j * m i * m j + c * ∑ i, (A i a * m i + A i b * m i)
          + c * ∑ i, (∑ j, A i j * m j + c * (A i a + A i b)) * pairVec a b i := by
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun i _ => ?_
      have : ∑ j, A i j * m i * m j = m i * ∑ j, A i j * m j := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun j _ => by ring
      rw [this]
      ring
    rw [e1, sum_mul_pairVec, hT a, hT b, Finset.sum_add_distrib, hA b a]
    ring
  have h4 : ∑ i, B i * (m i + c * pairVec a b i) = ∑ i, B i * m i + c * (B a + B b) := by
    rw [← sum_mul_pairVec a b B, Finset.mul_sum, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by ring
  simp only [Q, Fintype.sum_option, fuse, fuseB, extend]
  rw [h1, h2, h3, h4]
  ring

end RRFusion
