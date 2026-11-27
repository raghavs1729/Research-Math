/-
Nahm/Multisum.lean — multiple sums over tuples `m = (m_1, …, m_r) ∈ ℕ^r` with partial sums
`M_l = m_l + ⋯ + m_r`:

  `Σ_{m ∈ ℕ^r} Q^{Σ_l (M_l² + b_l M_l)} / ((Q)_{m_1} ⋯ (Q)_{m_r})`,

and their expression through the nested sums of the Bailey chain (`chainG`, built with `bstep`).
-/
import Nahm.Cauchy
import Nahm.Defs
import KRTW.Tuples

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K]

/-! ### Tuples -/

/-- `m_1 + ⋯ + m_r`. -/
def tot {r : ℕ} (m : Fin r → ℕ) : ℕ := ∑ l, m l

/-- `Σ_l (M_l² + b_l M_l)`. -/
def qexp {r : ℕ} (b m : Fin r → ℕ) : ℕ := ∑ l, (psum m l ^ 2 + b l * psum m l)

lemma tot_cons {r : ℕ} (a : ℕ) (m : Fin r → ℕ) : tot (Fin.cons a m : Fin (r + 1) → ℕ) = a + tot m := by
  simp [tot, Fin.sum_univ_succ]

lemma psum_zero {r : ℕ} (m : Fin (r + 1) → ℕ) : psum m 0 = tot m := by
  simp [psum, tot]

lemma psum_cons_zero {r : ℕ} (a : ℕ) (m : Fin r → ℕ) :
    psum (Fin.cons a m : Fin (r + 1) → ℕ) 0 = a + tot m := by
  rw [psum_zero, tot_cons]

lemma psum_cons_succ {r : ℕ} (a : ℕ) (m : Fin r → ℕ) (l : Fin r) :
    psum (Fin.cons a m : Fin (r + 1) → ℕ) l.succ = psum m l := by
  simp only [psum, Finset.sum_filter]
  rw [Fin.sum_univ_succ]
  simp [Fin.succ_le_succ_iff]

lemma qexp_cons {r : ℕ} (b0 a : ℕ) (b m : Fin r → ℕ) :
    qexp (Fin.cons b0 b : Fin (r + 1) → ℕ) (Fin.cons a m) = (a + tot m) ^ 2 + b0 * (a + tot m) + qexp b m := by
  simp only [qexp, Fin.sum_univ_succ, psum_cons_zero, psum_cons_succ, Fin.cons_zero, Fin.cons_succ]

lemma tot_le_qexp {r : ℕ} (b m : Fin r → ℕ) : tot m ≤ qexp b m := by
  cases r with
  | zero => simp [tot]
  | succ r =>
    have h : psum m 0 ^ 2 + b 0 * psum m 0 ≤ qexp b m :=
      single_le_sum (f := fun l => psum m l ^ 2 + b l * psum m l) (fun _ _ => Nat.zero_le _) (mem_univ 0)
    rw [psum_zero] at h
    nlinarith

lemma le_tot {r : ℕ} (m : Fin r → ℕ) (l : Fin r) : m l ≤ tot m :=
  single_le_sum (f := m) (fun _ _ => Nat.zero_le _) (mem_univ l)

lemma mem_box_of_tot_le {r N : ℕ} {m : Fin r → ℕ} (h : tot m ≤ N) : m ∈ KRTW.box r N := by
  simp only [KRTW.box, Fintype.mem_piFinset, mem_range]
  intro l
  exact Nat.lt_succ_of_le ((le_tot m l).trans h)

/-! ### The summands -/

variable (Q : K⟦X⟧)

/-- The summand `Q^{Σ_l (M_l² + b_l M_l)} / ∏_l (Q)_{m_l}`. -/
noncomputable def mterm {r : ℕ} (b m : Fin r → ℕ) : K⟦X⟧ := Q ^ qexp b m * ∏ l, iq Q (m l)

lemma mterm_cons {r : ℕ} (b0 a : ℕ) (b m : Fin r → ℕ) :
    mterm Q (Fin.cons b0 b : Fin (r + 1) → ℕ) (Fin.cons a m)
      = (Q ^ b0) ^ (a + tot m) * Q ^ ((a + tot m) ^ 2) * iq Q a * mterm Q b m := by
  simp only [mterm, qexp_cons, Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ]
  rw [← pow_mul, pow_add, pow_add]
  ring

/-- The nested sums of the Bailey chain: `chainG 0 c n = 1/(Q)_n`,
`chainG (r+1) c n = Σ_{k ≤ n} Q^{c_0 k + k²} chainG r (c_1, …) k / (Q)_{n-k}`. -/
noncomputable def chainG : (r : ℕ) → (Fin r → ℕ) → ℕ → K⟦X⟧
  | 0, _ => iq Q
  | r + 1, c => bstep Q (Q ^ c 0) (chainG r (Fin.tail c))

/-- Reindexing `n = a + t`. -/
lemma sum_shift {M : Type*} [AddCommMonoid M] (f : ℕ → M) (t N : ℕ) :
    ∑ n ∈ range (N + 1), (if t ≤ n then f n else 0)
      = ∑ a ∈ range (N + 1), (if a + t ≤ N then f (a + t) else 0) := by
  rw [← sum_filter, ← sum_filter]
  have : (range (N + 1)).filter (fun n => t ≤ n)
      = ((range (N + 1)).filter (fun a => a + t ≤ N)).map (addRightEmbedding t) := by
    ext n
    simp only [mem_filter, mem_range, mem_map, addRightEmbedding_apply]
    constructor
    · rintro ⟨h1, h2⟩
      exact ⟨n - t, ⟨by omega, by omega⟩, by omega⟩
    · rintro ⟨a, ⟨h1, h2⟩, rfl⟩
      exact ⟨by omega, by omega⟩
  rw [this, sum_map]
  rfl

/-- The nested sums as finite sums over tuples:
`chainG r c n = Σ_{m ∈ ℕ^r, |m| ≤ n} Q^{Σ_l (M_l² + c_l M_l)} / ((Q)_{n-|m|} ∏_l (Q)_{m_l})`. -/
theorem chainG_box : ∀ (r : ℕ) (c : Fin r → ℕ) (n N : ℕ), n ≤ N →
    chainG Q r c n = ∑ m ∈ KRTW.box r N, if tot m ≤ n then iq Q (n - tot m) * mterm Q c m else 0
  | 0, c, n, N, _ => by
    rw [sum_box_zero]
    simp [chainG, tot, mterm, qexp]
  | r + 1, c, n, N, h => by
    rw [chainG, bstep, sum_box_succ]
    have hc : c = Fin.cons (c 0) (Fin.tail c) := (Fin.cons_self_tail c).symm
    -- left side
    have hL : ∑ k ∈ range (n + 1), (Q ^ c 0) ^ k * Q ^ (k ^ 2) * iq Q (n - k) * chainG Q r (Fin.tail c) k
        = ∑ m ∈ KRTW.box r N, ∑ k ∈ range (n + 1), if tot m ≤ k then
            (Q ^ c 0) ^ k * Q ^ (k ^ 2) * iq Q (n - k) * (iq Q (k - tot m) * mterm Q (Fin.tail c) m)
            else 0 := by
      rw [sum_comm]
      refine sum_congr rfl (fun k hk => ?_)
      have hk' : k ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hk)
      rw [chainG_box r (Fin.tail c) k N (le_trans hk' h), mul_sum]
      refine sum_congr rfl (fun m _ => ?_)
      split_ifs <;> simp
    rw [hL]
    -- right side
    have hR : ∑ a ∈ range (N + 1), ∑ m ∈ KRTW.box r N,
        (if tot (Fin.cons a m : Fin (r + 1) → ℕ) ≤ n then
          iq Q (n - tot (Fin.cons a m : Fin (r + 1) → ℕ)) * mterm Q c (Fin.cons a m) else 0)
        = ∑ m ∈ KRTW.box r N, ∑ a ∈ range (n + 1), if a + tot m ≤ n then
            iq Q (n - (a + tot m)) * ((Q ^ c 0) ^ (a + tot m) * Q ^ ((a + tot m) ^ 2) * iq Q a *
              mterm Q (Fin.tail c) m) else 0 := by
      rw [sum_comm]
      refine sum_congr rfl (fun m _ => ?_)
      rw [← sum_subset (range_subset_range.mpr (Nat.succ_le_succ h))]
      · refine sum_congr rfl (fun a _ => ?_)
        rw [tot_cons]
        conv_lhs => rw [hc]
        rw [mterm_cons]
      · intro a _ ha
        have : n < a := by simpa using ha
        rw [tot_cons, if_neg (by omega)]
    rw [hR]
    refine sum_congr rfl (fun m _ => ?_)
    rw [sum_shift (fun k => (Q ^ c 0) ^ k * Q ^ (k ^ 2) * iq Q (n - k) *
      (iq Q (k - tot m) * mterm Q (Fin.tail c) m)) (tot m) n]
    refine sum_congr rfl (fun a _ => ?_)
    split_ifs
    · rw [Nat.add_sub_cancel]; ring
    · rfl

/-- The truncated multiple sum: `Σ_{m ∈ ℕ^{r+1}, |m| ≤ N} (summand) = Σ_{n ≤ N} Q^{b_0 n + n²} chainG r (b_1, …) n`. -/
theorem trunc_eq (r : ℕ) (b0 : ℕ) (c : Fin r → ℕ) (N : ℕ) :
    ∑ m ∈ KRTW.box (r + 1) N, (if tot m ≤ N then mterm Q (Fin.cons b0 c) m else 0)
      = ∑ n ∈ range (N + 1), (Q ^ b0) ^ n * Q ^ (n ^ 2) * chainG Q r c n := by
  rw [sum_box_succ]
  have hR : ∑ n ∈ range (N + 1), (Q ^ b0) ^ n * Q ^ (n ^ 2) * chainG Q r c n
      = ∑ m ∈ KRTW.box r N, ∑ n ∈ range (N + 1), if tot m ≤ n then
          (Q ^ b0) ^ n * Q ^ (n ^ 2) * (iq Q (n - tot m) * mterm Q c m) else 0 := by
    rw [sum_comm]
    refine sum_congr rfl (fun n hn => ?_)
    have hn' : n ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hn)
    rw [chainG_box Q r c n N hn', mul_sum]
    refine sum_congr rfl (fun m _ => ?_)
    split_ifs <;> simp
  rw [hR, sum_comm]
  refine sum_congr rfl (fun m _ => ?_)
  rw [sum_shift (fun n => (Q ^ b0) ^ n * Q ^ (n ^ 2) * (iq Q (n - tot m) * mterm Q c m)) (tot m) N]
  refine sum_congr rfl (fun a _ => ?_)
  rw [tot_cons, mterm_cons]
  split_ifs
  · rw [Nat.add_sub_cancel]; ring
  · rfl

variable [Nil Q]

lemma X_pow_tot_dvd_mterm {r : ℕ} (b m : Fin r → ℕ) : (X : K⟦X⟧) ^ tot m ∣ mterm Q b m :=
  Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (tot_le_qexp b m)) _

variable [TopologicalSpace K] [DiscreteTopology K]

omit [DiscreteTopology K] in
lemma summable_mterm {r : ℕ} (b : Fin r → ℕ) : Summable (mterm Q b) := by
  apply summable_of_dvd tot _ (X_pow_tot_dvd_mterm Q b)
  intro n
  exact (KRTW.box r n).finite_toSet.subset (fun m hm => mem_box_of_tot_le hm)

lemma coeff_tsum_mterm {r : ℕ} (b : Fin r → ℕ) (N : ℕ) :
    coeff N (∑' m, mterm Q b m)
      = coeff N (∑ m ∈ KRTW.box r N, if tot m ≤ N then mterm Q b m else 0) := by
  rw [coeff_tsum (summable_mterm Q b), map_sum]
  rw [tsum_eq_sum (s := KRTW.box r N)]
  · refine sum_congr rfl (fun m _ => ?_)
    split_ifs with h
    · rfl
    · rw [map_zero]
      exact X_pow_dvd_iff.mp (X_pow_tot_dvd_mterm Q b m) N (by omega)
  · intro m hm
    have : ¬ tot m ≤ N := fun h => hm (mem_box_of_tot_le h)
    exact X_pow_dvd_iff.mp (X_pow_tot_dvd_mterm Q b m) N (by omega)

/-- **The multiple sum as a sum over the outermost variable of the chain**:
`Σ_{m ∈ ℕ^{r+1}} Q^{Σ_l (M_l² + b_l M_l)}/∏(Q)_{m_l} = Σ_{n ≥ 0} Q^{b_0 n + n²} chainG r (b_1, …, b_r) n`. -/
theorem tsum_mterm_eq (r : ℕ) (b0 : ℕ) (c : Fin r → ℕ) :
    ∑' m, mterm Q (Fin.cons b0 c : Fin (r + 1) → ℕ) m
      = ∑' n, (Q ^ b0) ^ n * Q ^ (n ^ 2) * chainG Q r c n := by
  have hs : Summable fun n => (Q ^ b0) ^ n * Q ^ (n ^ 2) * chainG Q r c n := by
    apply summable_nat
    intro n
    exact Dvd.dvd.mul_right (Dvd.dvd.mul_left (X_pow_dvd_pow_of_le Q (by nlinarith)) _) _
  ext N
  rw [coeff_tsum_mterm, trunc_eq, coeff_tsum hs, map_sum]
  symm
  apply tsum_eq_sum
  intro n hn
  have : N < n := by simpa using hn
  exact X_pow_dvd_iff.mp (Dvd.dvd.mul_right (Dvd.dvd.mul_left (X_pow_dvd_pow_of_le Q
    (show n ≤ n ^ 2 by nlinarith)) _) _) N this

end Nahm
