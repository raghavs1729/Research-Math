/-
KRTW/Tuples.lean — bookkeeping for index tuples: peeling off the first (`Fin.cons`) or the last
(`Fin.snoc`) entry of a tuple, for the conventions `nxt`, `prv` of `KRTW/Defs.lean`.
-/
import KRTW.Defs

namespace KRTW

open Finset

/-- First entry of a tuple (`0` for the empty tuple). -/
def hd {j : ℕ} (r : Fin j → ℕ) : ℕ := if h : 0 < j then r ⟨0, h⟩ else 0

/-- `prv` with the convention `s_0 := B`. -/
def prvB (B : ℕ) {j : ℕ} (s : Fin j → ℕ) (i : Fin j) : ℕ :=
  if h : (i : ℕ) = 0 then B else s ⟨i - 1, by omega⟩

/-- `nxt` with the convention `M_{k+1} := B`. -/
def nxtB (B : ℕ) {j : ℕ} (M : Fin j → ℕ) (i : Fin j) : ℕ :=
  if h : (i : ℕ) + 1 < j then M ⟨i + 1, h⟩ else B

lemma prv_eq_prvB {j : ℕ} (s : Fin j → ℕ) (i : Fin j) : prv s i = prvB 0 s i := rfl

lemma nxt_eq_nxtB {j : ℕ} (r : Fin j → ℕ) (i : Fin j) : nxt r i = nxtB 0 r i := rfl

lemma hd_cons {j : ℕ} (a : ℕ) (r' : Fin j → ℕ) : hd (Fin.cons a r' : Fin (j + 1) → ℕ) = a := by
  simp [hd]

lemma nxt_cons_zero {j : ℕ} (a : ℕ) (r' : Fin j → ℕ) :
    nxt (Fin.cons a r' : Fin (j + 1) → ℕ) 0 = hd r' := by
  unfold nxt hd
  by_cases h : 0 < j
  · rw [dif_pos (by simpa using h), dif_pos h]
    exact Fin.cons_succ (α := fun _ => ℕ) a r' ⟨0, h⟩
  · rw [dif_neg (by simpa using h), dif_neg h]

lemma nxt_cons_succ {j : ℕ} (a : ℕ) (r' : Fin j → ℕ) (i : Fin j) :
    nxt (Fin.cons a r' : Fin (j + 1) → ℕ) i.succ = nxt r' i := by
  unfold nxt
  by_cases h : (i : ℕ) + 1 < j
  · rw [dif_pos (by simpa using h), dif_pos h]
    exact Fin.cons_succ (α := fun _ => ℕ) a r' ⟨i + 1, h⟩
  · rw [dif_neg (by simpa using h), dif_neg h]

lemma prvB_cons_zero {j : ℕ} (B b : ℕ) (s' : Fin j → ℕ) :
    prvB B (Fin.cons b s' : Fin (j + 1) → ℕ) 0 = B := by
  simp [prvB]

lemma prvB_cons_succ {j : ℕ} (B b : ℕ) (s' : Fin j → ℕ) (i : Fin j) :
    prvB B (Fin.cons b s' : Fin (j + 1) → ℕ) i.succ = prvB b s' i := by
  unfold prvB
  rw [dif_neg (by simp)]
  by_cases h : (i : ℕ) = 0
  · rw [dif_pos h]
    have : (⟨(i.succ : ℕ) - 1, by omega⟩ : Fin (j + 1)) = 0 := Fin.ext (by simp [h])
    rw [this, Fin.cons_zero]
  · rw [dif_neg h]
    have : (⟨(i.succ : ℕ) - 1, by omega⟩ : Fin (j + 1))
        = Fin.succ ⟨(i : ℕ) - 1, by omega⟩ := Fin.ext (by simp; omega)
    rw [this, Fin.cons_succ]

lemma chain_cons {j : ℕ} (a : ℕ) (r' : Fin j → ℕ) :
    (∀ i, nxt (Fin.cons a r' : Fin (j + 1) → ℕ) i ≤ (Fin.cons a r' : Fin (j + 1) → ℕ) i)
      ↔ hd r' ≤ a ∧ ∀ i, nxt r' i ≤ r' i := by
  rw [Fin.forall_fin_succ]
  simp only [nxt_cons_zero, nxt_cons_succ, Fin.cons_zero, Fin.cons_succ]

lemma chainS_cons {j : ℕ} (B b : ℕ) (s' : Fin j → ℕ) :
    (∀ i, prvB B (Fin.cons b s' : Fin (j + 1) → ℕ) i ≤ (Fin.cons b s' : Fin (j + 1) → ℕ) i)
      ↔ B ≤ b ∧ ∀ i, prvB b s' i ≤ s' i := by
  rw [Fin.forall_fin_succ]
  simp only [prvB_cons_zero, prvB_cons_succ, Fin.cons_zero, Fin.cons_succ]

lemma nxtB_snoc_last {j : ℕ} (B n : ℕ) (M' : Fin j → ℕ) :
    nxtB B (Fin.snoc M' n : Fin (j + 1) → ℕ) (Fin.last j) = B := by
  simp [nxtB]

lemma nxtB_snoc_castSucc {j : ℕ} (B n : ℕ) (M' : Fin j → ℕ) (i : Fin j) :
    nxtB B (Fin.snoc M' n : Fin (j + 1) → ℕ) i.castSucc = nxtB n M' i := by
  unfold nxtB
  rw [dif_pos (by simp)]
  by_cases h : (i : ℕ) + 1 < j
  · rw [dif_pos h]
    have : (⟨(i.castSucc : ℕ) + 1, by simp⟩ : Fin (j + 1))
        = Fin.castSucc ⟨(i : ℕ) + 1, h⟩ := Fin.ext (by simp)
    rw [this, Fin.snoc_castSucc]
  · rw [dif_neg h]
    have : (⟨(i.castSucc : ℕ) + 1, by simp⟩ : Fin (j + 1)) = Fin.last j :=
      Fin.ext (by simp; omega)
    rw [this, Fin.snoc_last]

lemma chainM_snoc {j : ℕ} (B n : ℕ) (M' : Fin j → ℕ) :
    (∀ i, nxtB B (Fin.snoc M' n : Fin (j + 1) → ℕ) i ≤ (Fin.snoc M' n : Fin (j + 1) → ℕ) i)
      ↔ B ≤ n ∧ ∀ i, nxtB n M' i ≤ M' i := by
  rw [Fin.forall_fin_succ']
  simp only [nxtB_snoc_last, nxtB_snoc_castSucc, Fin.snoc_last, Fin.snoc_castSucc]
  exact and_comm

/-! ### Sums over boxes -/

lemma sum_box_zero {M : Type*} [AddCommMonoid M] (N : ℕ) (f : (Fin 0 → ℕ) → M) :
    ∑ r ∈ box 0 N, f r = f Fin.elim0 := by
  have : box 0 N = {Fin.elim0} := by
    ext r
    simp only [box, Fintype.mem_piFinset, mem_singleton]
    exact ⟨fun _ => Subsingleton.elim _ _, fun _ i => i.elim0⟩
  rw [this, sum_singleton]

lemma sum_box_succ {M : Type*} [AddCommMonoid M] (j N : ℕ) (f : (Fin (j + 1) → ℕ) → M) :
    ∑ r ∈ box (j + 1) N, f r
      = ∑ a ∈ range (N + 1), ∑ r' ∈ box j N, f (Fin.cons a r') := by
  rw [← Finset.sum_product']
  refine Finset.sum_bij' (fun r _ => (r 0, Fin.tail r)) (fun p _ => Fin.cons p.1 p.2)
    ?_ ?_ ?_ ?_ ?_
  · intro r hr
    simp only [box, Fintype.mem_piFinset, mem_product] at hr ⊢
    exact ⟨hr 0, fun i => hr i.succ⟩
  · rintro ⟨a, r'⟩ hp
    simp only [box, Fintype.mem_piFinset, mem_product] at hp ⊢
    intro i
    refine Fin.cases ?_ (fun i => ?_) i
    · simpa using hp.1
    · simpa using hp.2 i
  · intro r _
    simp [Fin.cons_self_tail]
  · rintro ⟨a, r'⟩ _
    simp
  · intro r _
    simp [Fin.cons_self_tail]

lemma sum_box_snoc {M : Type*} [AddCommMonoid M] (j N : ℕ) (f : (Fin (j + 1) → ℕ) → M) :
    ∑ r ∈ box (j + 1) N, f r
      = ∑ a ∈ range (N + 1), ∑ r' ∈ box j N, f (Fin.snoc r' a) := by
  rw [← Finset.sum_product']
  refine Finset.sum_bij' (fun r _ => (r (Fin.last j), Fin.init r)) (fun p _ => Fin.snoc p.2 p.1)
    ?_ ?_ ?_ ?_ ?_
  · intro r hr
    simp only [box, Fintype.mem_piFinset, mem_product] at hr ⊢
    exact ⟨hr (Fin.last j), fun i => hr i.castSucc⟩
  · rintro ⟨a, r'⟩ hp
    simp only [box, Fintype.mem_piFinset, mem_product] at hp ⊢
    intro i
    refine Fin.lastCases ?_ (fun i => ?_) i
    · simpa using hp.1
    · simpa using hp.2 i
  · intro r _
    simp [Fin.snoc_init_self]
  · rintro ⟨a, r'⟩ _
    simp
  · intro r _
    simp [Fin.snoc_init_self]

end KRTW
