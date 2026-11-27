/-
KRTW/Bridge.lean — the recursively defined sums `Kf`, `Mel` of `KRTW/Limit.lean` are the explicit
sums over index tuples `krtwSum`, `melzerSum` of `KRTW/Defs.lean`; hence the first equality of
Theorem 2 for the latter.
-/
import KRTW.Limit
import KRTW.Tuples

namespace KRTW

open PowerSeries Finset

open scoped PowerSeries.WithPiTopology

set_option linter.unusedSectionVars false

variable {K : Type*} [Field K] [TopologicalSpace K] [DiscreteTopology K]

lemma qPoch_eq (n : ℕ) : (qPoch n : K⟦X⟧) = poch X n := rfl

lemma q2Poch_eq (n : ℕ) : (q2Poch n : K⟦X⟧) = poch (X ^ 2) n := by
  simp [q2Poch, poch, pow_mul]

/-- The weight of one level. -/
noncomputable def lvl (x : K⟦X⟧) (r s nr ps : ℕ) : K⟦X⟧ :=
  x ^ (r + 2 * s) * X ^ ((r + s) ^ 2 + s ^ 2 + s) * iq X (r - nr) * iq (X ^ 2) (s - ps)

/-- `krtwTerm` with the convention `s_0 := B`. -/
noncomputable def gterm (B : ℕ) (x : K⟦X⟧) {j : ℕ} (r s : Fin j → ℕ) : K⟦X⟧ :=
  ∏ i : Fin j, lvl x (r i) (s i) (nxt r i) (prvB B s i)

lemma krtwTerm_eq (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) : krtwTerm x r s = gterm 0 x r s := by
  simp only [krtwTerm, gterm, lvl, iq, qPoch_eq, q2Poch_eq, prv_eq_prvB]

lemma gterm_cons (B : ℕ) (x : K⟦X⟧) {j : ℕ} (a b : ℕ) (r' s' : Fin j → ℕ) :
    gterm B x (Fin.cons a r' : Fin (j + 1) → ℕ) (Fin.cons b s' : Fin (j + 1) → ℕ)
      = lvl x a b (hd r') B * gterm b x r' s' := by
  simp only [gterm, Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, nxt_cons_zero,
    nxt_cons_succ, prvB_cons_zero, prvB_cons_succ]

lemma tsum_shift_ite (B : ℕ) (u : ℕ → K⟦X⟧) (hu : ∀ n, (X : K⟦X⟧) ^ n ∣ u n) :
    ∑' m, u (m + B) = ∑' n, if B ≤ n then u n else 0 := by
  have hs : Summable (fun n => if B ≤ n then u n else 0) := by
    apply summable_nat
    intro n
    split_ifs
    · exact hu n
    · exact dvd_zero _
  rw [← hs.sum_add_tsum_nat_add B, sum_eq_zero, zero_add]
  · refine tsum_congr (fun m => ?_)
    rw [if_pos (Nat.le_add_left B m)]
  · intro i hi
    rw [if_neg (by simpa using hi)]

/-- Truncation of one level to `r, s ≤ N`, modulo `X^{N+1}`. -/
lemma Lv_trunc (ω : ℕ → K⟦X⟧) (j B : ℕ) (x : K⟦X⟧) (N : ℕ) :
    (X : K⟦X⟧) ^ (N + 1) ∣ Lv ω j B x - ∑ a ∈ range (N + 1), ∑ b ∈ range (N + 1),
      if B ≤ b then X ^ (a ^ 2 + 2 * a * b + tri (2 * b)) *
        (ω a * x ^ (a + 2 * b) * iq (X ^ 2) (b - B) * Tf j a b x) else 0 := by
  set u : ℕ → ℕ → K⟦X⟧ := fun a b => X ^ (a ^ 2 + 2 * a * b + tri (2 * b)) *
    (ω a * x ^ (a + 2 * b) * iq (X ^ 2) (b - B) * Tf j a b x) with hu
  have hub : ∀ a b, (X : K⟦X⟧) ^ b ∣ u a b := by
    intro a b
    apply dvd_X_pow_mul
    have := le_tri (2 * b)
    omega
  have hua : ∀ a b, (X : K⟦X⟧) ^ a ∣ u a b := by
    intro a b
    apply dvd_X_pow_mul
    have : a ≤ a ^ 2 := Nat.le_self_pow (by norm_num) a
    omega
  have hLv : Lv ω j B x = ∑' a, ∑' b, if B ≤ b then u a b else 0 := by
    unfold Lv
    refine tsum_congr (fun a => ?_)
    rw [← tsum_shift_ite B (u a) (hub a)]
    refine tsum_congr (fun b => ?_)
    simp only [hu, Nat.add_sub_cancel]
  have h1 : (X : K⟦X⟧) ^ (N + 1) ∣ (∑' a, ∑' b, if B ≤ b then u a b else 0)
      - ∑ a ∈ range (N + 1), ∑' b, if B ≤ b then u a b else 0 := by
    refine X_pow_dvd_tsum_sub_sum (summable_nat (fun a => X_pow_dvd_tsum (fun b => ?_))) _
      (fun a ha => X_pow_dvd_tsum (fun b => ?_))
    · split_ifs
      · exact hua a b
      · exact dvd_zero _
    · have ha' : N + 1 ≤ a := by simpa using ha
      split_ifs
      · exact (pow_dvd_pow X ha').trans (hua a b)
      · exact dvd_zero _
  have h2 : ∀ a, (X : K⟦X⟧) ^ (N + 1) ∣ (∑' b, if B ≤ b then u a b else 0)
      - ∑ b ∈ range (N + 1), if B ≤ b then u a b else 0 := by
    intro a
    refine X_pow_dvd_tsum_sub_sum (summable_nat (fun b => ?_)) _ (fun b hb => ?_)
    · split_ifs
      · exact hub a b
      · exact dvd_zero _
    · have hb' : N + 1 ≤ b := by simpa using hb
      split_ifs
      · exact (pow_dvd_pow X hb').trans (hub a b)
      · exact dvd_zero _
  rw [hLv]
  refine dvd_sub_trans h1 ?_
  rw [← sum_sub_distrib]
  exact dvd_sum (fun a _ => h2 a)

/-- The finite Fubini step: peeling off the first entries of `r` and `s`. -/
lemma Lv_cons (ω : ℕ → K⟦X⟧) (j B : ℕ) (x : K⟦X⟧) (N : ℕ) :
    (∑ a ∈ range (N + 1), ∑ b ∈ range (N + 1),
      if B ≤ b then X ^ (a ^ 2 + 2 * a * b + tri (2 * b)) *
        (ω a * x ^ (a + 2 * b) * iq (X ^ 2) (b - B) *
          ∑ r ∈ box j N, ∑ s ∈ box j N,
            if (hd r ≤ a ∧ ∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB b s i ≤ s i)
              then poch X a * iq X (a - hd r) * gterm b x r s else 0) else 0)
      = ∑ r ∈ box (j + 1) N, ∑ s ∈ box (j + 1) N,
          if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i)
            then ω (hd r) * poch X (hd r) * gterm B x r s else 0 := by
  rw [sum_box_succ]
  refine sum_congr rfl (fun a _ => ?_)
  rw [sum_congr rfl (fun r' _ => sum_box_succ j N _), sum_comm]
  refine sum_congr rfl (fun b _ => ?_)
  by_cases hb : B ≤ b
  · rw [if_pos hb, mul_sum, mul_sum]
    refine sum_congr rfl (fun r' _ => ?_)
    rw [mul_sum, mul_sum]
    refine sum_congr rfl (fun s' _ => ?_)
    simp only [chain_cons, chainS_cons, hd_cons, gterm_cons]
    by_cases h1 : (hd r' ≤ a ∧ ∀ i, nxt r' i ≤ r' i) ∧ (∀ i, prvB b s' i ≤ s' i)
    · rw [if_pos h1, if_pos ⟨h1.1, hb, h1.2⟩]
      unfold lvl
      have he : a ^ 2 + 2 * a * b + tri (2 * b) = (a + b) ^ 2 + b ^ 2 + b := by
        rw [tri_two_mul]; ring
      rw [he]
      ring
    · rw [if_neg h1, if_neg (by tauto)]
      simp
  · rw [if_neg hb]
    symm
    refine sum_eq_zero (fun r' _ => sum_eq_zero (fun s' _ => ?_))
    simp only [chainS_cons]
    rw [if_neg (by tauto)]

/-- One level, modulo `X^{N+1}`, given the tuple form of the inner sums. -/
lemma Lv_box (ω : ℕ → K⟦X⟧) (j B : ℕ) (x : K⟦X⟧) (N : ℕ)
    (ih : ∀ a b : ℕ, (X : K⟦X⟧) ^ (N + 1) ∣ Tf j a b x - ∑ r ∈ box j N, ∑ s ∈ box j N,
      if (hd r ≤ a ∧ ∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB b s i ≤ s i)
        then poch X a * iq X (a - hd r) * gterm b x r s else 0) :
    (X : K⟦X⟧) ^ (N + 1) ∣ Lv ω j B x - ∑ r ∈ box (j + 1) N, ∑ s ∈ box (j + 1) N,
      if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i)
        then ω (hd r) * poch X (hd r) * gterm B x r s else 0 := by
  refine dvd_sub_trans (Lv_trunc ω j B x N) ?_
  rw [← Lv_cons, ← sum_sub_distrib]
  refine dvd_sum (fun a _ => ?_)
  rw [← sum_sub_distrib]
  refine dvd_sum (fun b _ => ?_)
  split_ifs
  · rw [← mul_sub, ← mul_sub]
    exact Dvd.dvd.mul_left (Dvd.dvd.mul_left (ih a b) _) _
  · simp

/-- The recursively defined `T_j(R,B;x)` is the sum over tuples, modulo `X^{N+1}`. -/
theorem Tf_box : ∀ (j N R B : ℕ) (x : K⟦X⟧),
    (X : K⟦X⟧) ^ (N + 1) ∣ Tf j R B x - ∑ r ∈ box j N, ∑ s ∈ box j N,
      if (hd r ≤ R ∧ ∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i)
        then poch X R * iq X (R - hd r) * gterm B x r s else 0
  | 0, N, R, B, x => by
    rw [sum_box_zero, sum_box_zero, Tf_zero]
    have := poch_mul_iq (X : K⟦X⟧) R
    simp [hd, gterm, this]
  | j + 1, N, R, B, x => by
    have h := Lv_box (qbinom X R) j B x N (fun a b => Tf_box j N a b x)
    rw [← Tf_succ_eq_Lv] at h
    have e : ∀ r s : Fin (j + 1) → ℕ,
        (if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i)
          then qbinom X R (hd r) * poch X (hd r) * gterm B x r s else 0)
        = (if (hd r ≤ R ∧ ∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i)
          then poch X R * iq X (R - hd r) * gterm B x r s else 0) := by
      intro r s
      by_cases hR : hd r ≤ R
      · have hq : (qbinom X R (hd r) * poch X (hd r) : K⟦X⟧) = poch X R * iq X (R - hd r) := by
          simp only [qbinom, hR, if_true]
          calc (poch X R * iq X (hd r) * iq X (R - hd r) * poch X (hd r) : K⟦X⟧)
              = poch X R * iq X (R - hd r) * (iq X (hd r) * poch X (hd r)) := by ring
            _ = _ := by rw [iq_mul_poch, mul_one]
        rw [hq]
        by_cases hc : (∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i)
        · rw [if_pos hc, if_pos ⟨⟨hR, hc.1⟩, hc.2⟩]
        · rw [if_neg hc, if_neg (by tauto)]
      · rw [qbinom_of_lt X (not_le.mp hR),
          if_neg (fun h : (hd r ≤ R ∧ ∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB B s i ≤ s i) => hR h.1.1)]
        simp
    rw [sum_congr rfl (fun r _ => sum_congr rfl (fun s _ => e r s))] at h
    exact h

/-- The recursively defined `KRTW_k(x)` is the sum over tuples, modulo `X^{N+1}`. -/
theorem Kf_box (k N : ℕ) (x : K⟦X⟧) :
    (X : K⟦X⟧) ^ (N + 1) ∣ Kf k x - ∑ r ∈ box k N, ∑ s ∈ box k N,
      if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prv s i ≤ s i) then krtwTerm x r s else 0 := by
  cases k with
  | zero =>
    rw [sum_box_zero, sum_box_zero]
    simp [Kf, krtwTerm]
  | succ j =>
    have h := Lv_box (iq X) j 0 x N (fun a b => Tf_box j N a b x)
    have e : ∀ r s : Fin (j + 1) → ℕ,
        (if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prvB 0 s i ≤ s i)
          then iq X (hd r) * poch X (hd r) * gterm 0 x r s else 0)
        = (if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prv s i ≤ s i) then krtwTerm x r s else 0) := by
      intro r s
      simp only [krtwTerm_eq, prv_eq_prvB, iq_mul_poch, one_mul]
    rw [sum_congr rfl (fun r _ => sum_congr rfl (fun s _ => e r s))] at h
    exact h

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma coeff_eq_of_dvd_sub {f g : K⟦X⟧} {N : ℕ} (h : (X : K⟦X⟧) ^ (N + 1) ∣ f - g) :
    coeff N f = coeff N g := by
  have := X_pow_dvd_iff.mp h N (Nat.lt_succ_self N)
  rwa [map_sub, sub_eq_zero] at this

theorem krtwSum_eq_Kf (x : K⟦X⟧) (k : ℕ) : krtwSum x k = Kf k x := by
  ext N
  rw [coeff_eq_of_dvd_sub (Kf_box k N x), krtwSum, coeff_mk, map_sum]
  refine sum_congr rfl (fun r _ => ?_)
  rw [map_sum]
  refine sum_congr rfl (fun s _ => ?_)
  split_ifs <;> simp

/-! ### The Melzer side -/

/-- `melzerTerm` with the convention `M_{k+1} := B`. -/
noncomputable def mterm (B : ℕ) (x : K⟦X⟧) {j : ℕ} (M : Fin j → ℕ) : K⟦X⟧ :=
  ∏ i : Fin j, x ^ (M i) * X ^ tri (M i) * iq X (M i - nxtB B M i)

lemma melzerTerm_eq (x : K⟦X⟧) {k : ℕ} (M : Fin k → ℕ) : melzerTerm x M = mterm 0 x M := by
  simp only [melzerTerm, mterm, ← tri_eq, iq, qPoch_eq, nxt_eq_nxtB]

lemma mterm_snoc (B : ℕ) (x : K⟦X⟧) {j : ℕ} (n : ℕ) (M' : Fin j → ℕ) :
    mterm B x (Fin.snoc M' n : Fin (j + 1) → ℕ)
      = mterm n x M' * (x ^ n * X ^ tri n * iq X (n - B)) := by
  simp only [mterm, Fin.prod_univ_castSucc, Fin.snoc_castSucc, Fin.snoc_last,
    nxtB_snoc_castSucc, nxtB_snoc_last]

/-- The recursively defined `H_j(B;x)` is the sum over tuples, modulo `X^{N+1}`. -/
theorem Hf_box : ∀ (j N B : ℕ) (x : K⟦X⟧),
    (X : K⟦X⟧) ^ (N + 1) ∣ Hf j B x - ∑ M ∈ box j N,
      if ∀ i, nxtB B M i ≤ M i then mterm B x M else 0
  | 0, N, B, x => by
    rw [sum_box_zero]
    simp [Hf_zero, mterm]
  | j + 1, N, B, x => by
    set u : ℕ → K⟦X⟧ := fun n => X ^ tri n * (x ^ n * iq X (n - B) * Hf j n x) with hu
    have hud : ∀ n, (X : K⟦X⟧) ^ n ∣ u n := fun n => dvd_X_pow_mul (le_tri n) _
    have h0 : Hf (j + 1) B x = ∑' n, if B ≤ n then u n else 0 := by
      rw [Hf_succ, ← tsum_shift_ite B u hud]
      refine tsum_congr (fun m => ?_)
      simp only [hu, Nat.add_sub_cancel]
    have h1 : (X : K⟦X⟧) ^ (N + 1) ∣ (∑' n, if B ≤ n then u n else 0)
        - ∑ n ∈ range (N + 1), if B ≤ n then u n else 0 := by
      refine X_pow_dvd_tsum_sub_sum (summable_nat (fun n => ?_)) _ (fun n hn => ?_)
      · split_ifs
        · exact hud n
        · exact dvd_zero _
      · have hn' : N + 1 ≤ n := by simpa using hn
        split_ifs
        · exact (pow_dvd_pow X hn').trans (hud n)
        · exact dvd_zero _
    rw [h0]
    refine dvd_sub_trans h1 ?_
    rw [sum_box_snoc, ← sum_sub_distrib]
    refine dvd_sum (fun n _ => ?_)
    by_cases hB : B ≤ n
    · rw [if_pos hB]
      have e : (∑ M' ∈ box j N,
          if ∀ i, nxtB B (Fin.snoc M' n : Fin (j + 1) → ℕ) i ≤ (Fin.snoc M' n : Fin (j + 1) → ℕ) i
            then mterm B x (Fin.snoc M' n : Fin (j + 1) → ℕ) else 0)
          = X ^ tri n * (x ^ n * iq X (n - B) *
            ∑ M' ∈ box j N, if ∀ i, nxtB n M' i ≤ M' i then mterm n x M' else 0) := by
        rw [mul_sum, mul_sum]
        refine sum_congr rfl (fun M' _ => ?_)
        simp only [chainM_snoc, mterm_snoc]
        by_cases h : ∀ i, nxtB n M' i ≤ M' i
        · rw [if_pos ⟨hB, h⟩, if_pos h]
          ring
        · rw [if_neg (by tauto), if_neg h]
          simp
      rw [e, hu, ← mul_sub, ← mul_sub]
      exact Dvd.dvd.mul_left (Dvd.dvd.mul_left (Hf_box j N n x) _) _
    · rw [if_neg hB]
      have e : (∑ M' ∈ box j N,
          if ∀ i, nxtB B (Fin.snoc M' n : Fin (j + 1) → ℕ) i ≤ (Fin.snoc M' n : Fin (j + 1) → ℕ) i
            then mterm B x (Fin.snoc M' n : Fin (j + 1) → ℕ) else 0) = 0 := by
        refine sum_eq_zero (fun M' _ => ?_)
        simp only [chainM_snoc]
        rw [if_neg (by tauto)]
      rw [e, sub_zero]
      exact dvd_zero _

theorem melzerSum_eq_Mel (x : K⟦X⟧) (k : ℕ) : melzerSum x k = Mel k x := by
  ext N
  rw [Mel, coeff_eq_of_dvd_sub (Hf_box k N 0 x), melzerSum, coeff_mk, map_sum]
  refine sum_congr rfl (fun M _ => ?_)
  simp only [melzerTerm_eq, nxt_eq_nxtB]
  split_ifs <;> simp

/-- **Theorem 2, first equality**: the `2k`-fold KRTW sum equals the `k`-fold Melzer-type sum, for
every `k` and every power series `x`. -/
theorem krtwSum_eq_melzerSum' (x : K⟦X⟧) (k : ℕ) : krtwSum x k = melzerSum x k := by
  rw [krtwSum_eq_Kf, melzerSum_eq_Mel, Kf_eq_Mel]

end KRTW
