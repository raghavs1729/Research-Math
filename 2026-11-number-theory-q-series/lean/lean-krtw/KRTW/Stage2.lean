/-
KRTW/Stage2.lean — the evaluation of the `k`-fold sum at `x = 1` (Proposition 6 and the first half of
the proof of Theorem 3 of `proofs/KRTW-general.md`):

  `(q;q²)_∞ (q;q)_∞ · Mel_k(1) = Σ_{n ≥ 0} (-1)^n q^{((k+2)n²+kn)/2} (1 - q^{2n+1})`,

from the Bailey chain of `KRTW/Bailey.lean` by letting `n → ∞` (congruences modulo `X^{n+1}`).
-/
import KRTW.Bridge
import KRTW.Bailey
import KRTW.JTP

namespace KRTW

open PowerSeries Finset

open scoped PowerSeries.WithPiTopology

set_option linter.unusedSectionVars false

variable {K : Type*} [Field K] [TopologicalSpace K] [DiscreteTopology K]

/-! ### Congruences modulo a power of `X` -/

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma cong_mul {n : ℕ} {a b c d : K⟦X⟧} (h1 : (X : K⟦X⟧) ^ n ∣ a - b)
    (h2 : (X : K⟦X⟧) ^ n ∣ c - d) : (X : K⟦X⟧) ^ n ∣ a * c - b * d := by
  have : a * c - b * d = (a - b) * c + b * (c - d) := by ring
  rw [this]
  exact dvd_add (h1.mul_right _) (h2.mul_left _)

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma X_pow_dvd_prod_sub_one {ι : Type*} (s : Finset ι) (f : ι → K⟦X⟧) (n : ℕ)
    (h : ∀ i ∈ s, (X : K⟦X⟧) ^ n ∣ 1 - f i) : (X : K⟦X⟧) ^ n ∣ ∏ i ∈ s, f i - 1 :=
  (dvd_sub_comm).mp (X_pow_dvd_one_sub_prod s f n h)

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma X_pow_dvd_gp_sub_one (y : K⟦X⟧) (m e : ℕ) (h : (X : K⟦X⟧) ^ e ∣ y) :
    (X : K⟦X⟧) ^ e ∣ gp y m - 1 := by
  apply X_pow_dvd_prod_sub_one
  intro i _
  have : (1 : K⟦X⟧) - (1 - y * X ^ i) = y * X ^ i := by ring
  rw [this]
  exact h.mul_right _

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma poch_mul_iq_sub (q : K⟦X⟧) [Nil q] {n m : ℕ} (h : m ≤ n) :
    poch q n * iq q (n - m) = ∏ i ∈ range m, (1 - q ^ (n - m + i + 1)) := by
  have hp : poch q n = poch q (n - m) * ∏ i ∈ range m, (1 - q ^ (n - m + i + 1)) := by
    have hn : n = (n - m) + m := by omega
    conv_lhs => rw [hn]
    simp only [poch]
    rw [prod_range_add]
  rw [hp]
  calc (poch q (n - m) * ∏ i ∈ range m, (1 - q ^ (n - m + i + 1))) * iq q (n - m)
      = (poch q (n - m) * iq q (n - m)) * ∏ i ∈ range m, (1 - q ^ (n - m + i + 1)) := by ring
    _ = _ := by rw [poch_mul_iq, one_mul]

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma X_pow_dvd_poch_mul_iq_sub_one {n m : ℕ} (h : m ≤ n) :
    (X : K⟦X⟧) ^ (n - m + 1) ∣ poch X n * iq X (n - m) - 1 := by
  rw [poch_mul_iq_sub X h]
  apply X_pow_dvd_prod_sub_one
  intro i _
  have : (1 : K⟦X⟧) - (1 - X ^ (n - m + i + 1)) = X ^ (n - m + i + 1) := by ring
  rw [this]
  exact pow_dvd_pow X (by omega)

/-! ### The chain sums `Gc` as sums over tuples, and `(q)_n G_k(n) ≡ Mel_k mod q^{n+1}` -/

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma mterm_cons (x : K⟦X⟧) {j : ℕ} (m : ℕ) (M' : Fin j → ℕ) :
    mterm 0 x (Fin.cons m M' : Fin (j + 1) → ℕ)
      = (x ^ m * X ^ tri m * iq X (m - hd M')) * mterm 0 x M' := by
  simp only [mterm, Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, ← nxt_eq_nxtB,
    nxt_cons_zero, nxt_cons_succ]

omit [TopologicalSpace K] [DiscreteTopology K] in
theorem Gc_box (x : K⟦X⟧) : ∀ (j n N : ℕ), n ≤ N →
    Gc x j n = ∑ M ∈ box j N,
      if hd M ≤ n ∧ ∀ i, nxt M i ≤ M i then iq X (n - hd M) * mterm 0 x M else 0
  | 0, n, N, _ => by
    rw [sum_box_zero]
    simp [Gc_zero, hd, mterm]
  | j + 1, n, N, h => by
    rw [Gc_succ, sum_box_succ, ← sum_subset (range_subset_range.mpr (Nat.succ_le_succ h))]
    · refine sum_congr rfl (fun m hm => ?_)
      have hm' : m ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hm)
      rw [Gc_box x j m N (le_trans hm' h), mul_sum, mul_sum]
      refine sum_congr rfl (fun M' _ => ?_)
      simp only [hd_cons, chain_cons, mterm_cons]
      by_cases hc : hd M' ≤ m ∧ ∀ i, nxt M' i ≤ M' i
      · rw [if_pos hc, if_pos ⟨hm', hc⟩]
        ring
      · rw [if_neg hc, if_neg (by tauto)]
        simp
    · intro m _ hm
      refine sum_eq_zero (fun M' _ => ?_)
      simp only [hd_cons]
      rw [if_neg]
      intro hh
      exact hm (mem_range.mpr (Nat.lt_succ_of_le hh.1))

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma X_pow_hd_dvd_mterm (x : K⟦X⟧) : ∀ {j : ℕ} (M : Fin j → ℕ),
    (X : K⟦X⟧) ^ hd M ∣ mterm 0 x M
  | 0, M => by simp [hd]
  | j + 1, M => by
    have h0 : hd M = M 0 := by simp [hd]
    rw [h0]
    refine Dvd.dvd.trans ?_ (Finset.dvd_prod_of_mem _ (Finset.mem_univ (0 : Fin (j + 1))))
    exact Dvd.dvd.mul_right (Dvd.dvd.mul_left (pow_dvd_pow X (le_tri _)) _) _

/-- `(q)_n G_k(n) ≡ Mel_k(x) mod q^{n+1}`. -/
theorem poch_Gc_cong (x : K⟦X⟧) (k n : ℕ) :
    (X : K⟦X⟧) ^ (n + 1) ∣ poch X n * Gc x k n - Mel k x := by
  have hH := Hf_box k n 0 x
  have h1 : (X : K⟦X⟧) ^ (n + 1) ∣ poch X n * Gc x k n
      - ∑ M ∈ box k n, if ∀ i, nxtB 0 M i ≤ M i then mterm 0 x M else 0 := by
    rw [Gc_box x k n n le_rfl, mul_sum, ← sum_sub_distrib]
    refine dvd_sum (fun M hM => ?_)
    have hhd : hd M ≤ n := by
      unfold hd
      split_ifs with hk
      · exact Nat.lt_succ_iff.mp (mem_range.mp (Fintype.mem_piFinset.mp hM _))
      · exact Nat.zero_le _
    by_cases hc : ∀ i, nxt M i ≤ M i
    · have hc' : ∀ i, nxtB 0 M i ≤ M i := hc
      rw [if_pos ⟨hhd, hc⟩, if_pos hc']
      have e : poch X n * (iq X (n - hd M) * mterm 0 x M) - mterm 0 x M
          = (poch X n * iq X (n - hd M) - 1) * mterm 0 x M := by ring
      have e2 : (X : K⟦X⟧) ^ (n + 1) = X ^ (n - hd M + 1) * X ^ hd M := by
        rw [← pow_add]; congr 1; omega
      rw [e, e2]
      exact mul_dvd_mul (X_pow_dvd_poch_mul_iq_sub_one hhd) (X_pow_hd_dvd_mterm x M)
    · have hc' : ¬ ∀ i, nxtB 0 M i ≤ M i := hc
      rw [if_neg (fun h => hc h.2), if_neg hc']
      simp
  exact dvd_sub_trans h1 ((dvd_sub_comm).mp hH)

/-! ### `x = 1`: the limit `n → ∞` of the Bailey chain -/

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma gp_XX (r : ℕ) : (1 - X) * gp ((X : K⟦X⟧) * X) r = poch X (r + 1) := by
  simp only [gp, poch]
  rw [prod_range_succ']
  simp only [zero_add, pow_one]
  rw [mul_comm]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  ring

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- At `a = q`: `α_r (1 - q) = (-1)^r q^{r(r-1)/2} (1 - q^{2r+1})`. -/
lemma alpha_one (r : ℕ) :
    alpha (X : K⟦X⟧) r * (1 - X) = (-1) ^ r * X ^ tri (r - 1) * (1 - X ^ (2 * r + 1)) := by
  cases r with
  | zero => simp [alpha]
  | succ r =>
    have h := gp_XX (K := K) r
    have h2 := iq_mul_poch (X : K⟦X⟧) (r + 1)
    simp only [alpha, Nat.add_sub_cancel]
    linear_combination ((-1) ^ (r + 1) * X ^ tri r * (1 - X * X ^ (2 * r + 2)) * iq X (r + 1)) * h
      + ((-1) ^ (r + 1) * X ^ tri r * (1 - X * X ^ (2 * r + 2))) * h2

/-- The partial sums of the theta series: `Σ_{r ≤ n} (-1)^r q^{(k+1)r(r+1)/2 + r(r-1)/2} (1 - q^{2r+1})`. -/
noncomputable def ThetaN (k n : ℕ) : K⟦X⟧ :=
  ∑ r ∈ range (n + 1), X ^ ((k + 1) * tri r) * ((-1) ^ r * X ^ tri (r - 1) * (1 - X ^ (2 * r + 1)))

/-- `(q)_n Mel_k(1) ≡ (-q;q)_n Σ_{r ≤ n} (-1)^r q^{((k+2)r²+kr)/2} (1 - q^{2r+1}) mod q^{n+1}`. -/
theorem limit_one (k n : ℕ) :
    (X : K⟦X⟧) ^ (n + 1) ∣ poch X n * Mel k (1 : K⟦X⟧) - gp (-X) n * ThetaN k n := by
  have hc := bailey_chain (1 : K⟦X⟧) k n
  simp only [one_pow, one_mul, mul_one] at hc
  have hG := poch_Gc_cong (1 : K⟦X⟧) k n
  have hp : poch (X : K⟦X⟧) (2 * n + 1) = poch X n * ∏ i ∈ range (n + 1), (1 - X ^ (n + i + 1)) := by
    have : 2 * n + 1 = n + (n + 1) := by omega
    rw [this]
    simp only [poch]
    rw [prod_range_add]
  have h1 : (X : K⟦X⟧) ^ (n + 1) ∣ poch X (2 * n + 1) * Gc 1 k n - Mel k 1 := by
    have e : poch X (2 * n + 1) * Gc 1 k n - Mel k (1 : K⟦X⟧)
        = (∏ i ∈ range (n + 1), (1 - X ^ (n + i + 1))) * (poch X n * Gc 1 k n) - 1 * Mel k 1 := by
      rw [hp]; ring
    rw [e]
    refine cong_mul (X_pow_dvd_prod_sub_one _ _ _ (fun i _ => ?_)) hG
    have : (1 : K⟦X⟧) - (1 - X ^ (n + i + 1)) = X ^ (n + i + 1) := by ring
    rw [this]
    exact pow_dvd_pow X (by omega)
  have h2 : poch X n * (poch X (2 * n + 1) * Gc (1 : K⟦X⟧) k n)
      = gp (-X) n * ∑ r ∈ range (n + 1), X ^ ((k + 1) * tri r) *
          ((-1) ^ r * X ^ tri (r - 1) * (1 - X ^ (2 * r + 1))) *
            ((poch X n * iq X (n - r)) * gp (X * X ^ (n + r + 1)) (n - r)) := by
    rw [← gp_XX (2 * n), mul_assoc (1 - X), hc, BPsum]
    simp only [mul_sum]
    refine sum_congr rfl (fun r _ => ?_)
    rw [← alpha_one r]
    ring
  have h3 : (X : K⟦X⟧) ^ (n + 1) ∣ poch X n * (poch X (2 * n + 1) * Gc 1 k n)
      - gp (-X) n * ThetaN k n := by
    rw [h2, ThetaN, ← mul_sub, ← sum_sub_distrib]
    refine Dvd.dvd.mul_left (dvd_sum (fun r hr => ?_)) _
    have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
    have e : X ^ ((k + 1) * tri r) * ((-1) ^ r * X ^ tri (r - 1) * (1 - X ^ (2 * r + 1))) *
          ((poch X n * iq X (n - r)) * gp (X * X ^ (n + r + 1)) (n - r))
        - X ^ ((k + 1) * tri r) * ((-1 : K⟦X⟧) ^ r * X ^ tri (r - 1) * (1 - X ^ (2 * r + 1)))
        = X ^ ((k + 1) * tri r) * (((poch X n * iq X (n - r)) * gp (X * X ^ (n + r + 1)) (n - r)
            - 1 * 1) * ((-1) ^ r * X ^ tri (r - 1) * (1 - X ^ (2 * r + 1)))) := by ring
    have e2 : (X : K⟦X⟧) ^ (n + 1) = X ^ r * X ^ (n - r + 1) := by
      rw [← pow_add]; congr 1; omega
    rw [e, e2]
    refine mul_dvd_mul (pow_dvd_pow X ?_) (Dvd.dvd.mul_right (cong_mul
      (X_pow_dvd_poch_mul_iq_sub_one hr') (X_pow_dvd_gp_sub_one _ _ _ ?_)) _)
    · have := le_tri r
      nlinarith
    · exact Dvd.dvd.mul_left (pow_dvd_pow X (by omega)) _
  have e : poch X n * Mel k (1 : K⟦X⟧) - gp (-X) n * ThetaN k n
      = (poch X n * (poch X (2 * n + 1) * Gc 1 k n) - gp (-X) n * ThetaN k n)
        - poch X n * (poch X (2 * n + 1) * Gc 1 k n - Mel k 1) := by ring
  rw [e]
  exact dvd_sub h3 (h1.mul_left _)

/-! ### Assembly with the infinite products -/

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma gpneg_mul_poch (n : ℕ) : gp (-(X : K⟦X⟧)) n * poch X n = poch (X ^ 2) n := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [gp_succ, poch_succ, poch_succ, ← ih]
    ring

/-- `(q;q²)_n`. -/
noncomputable def oddP (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (1 + 2 * i))

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma oddP_mul_poch2 (n : ℕ) : (oddP n : K⟦X⟧) * poch (X ^ 2) n = poch X (2 * n) := by
  induction n with
  | zero => simp [oddP]
  | succ n ih =>
    have e : 2 * (n + 1) = 2 * n + 1 + 1 := by ring
    rw [e, poch_succ, poch_succ, poch_succ, ← ih, oddP, prod_range_succ, ← oddP]
    ring

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- `(-q;q)_n (q;q²)_n = (q^{n+1};q)_n`. -/
lemma gpneg_mul_oddP (n : ℕ) :
    gp (-(X : K⟦X⟧)) n * oddP n = ∏ i ∈ range n, (1 - X ^ (n + i + 1)) := by
  have hne : poch ((X : K⟦X⟧) ^ 2) n ≠ 0 := by
    intro h
    have := constantCoeff_poch ((X : K⟦X⟧) ^ 2) n
    rw [h] at this
    simp at this
  have hp : poch (X : K⟦X⟧) (2 * n) = poch X n * ∏ i ∈ range n, (1 - X ^ (n + i + 1)) := by
    have : 2 * n = n + n := by omega
    rw [this]
    simp only [poch]
    rw [prod_range_add]
  apply mul_right_cancel₀ hne
  calc gp (-(X : K⟦X⟧)) n * oddP n * poch (X ^ 2) n
      = gp (-X) n * (oddP n * poch (X ^ 2) n) := by ring
    _ = gp (-X) n * (poch X n * ∏ i ∈ range n, (1 - X ^ (n + i + 1))) := by
        rw [oddP_mul_poch2, hp]
    _ = (∏ i ∈ range n, (1 - X ^ (n + i + 1))) * poch (X ^ 2) n := by
        rw [← gpneg_mul_poch]; ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma theta_exp (k r : ℕ) : ((k + 2) * r ^ 2 + k * r) / 2 = (k + 1) * tri r + tri (r - 1) := by
  apply Nat.div_eq_of_eq_mul_right (by norm_num)
  cases r with
  | zero => simp
  | succ r =>
    have h1 := two_mul_tri (r + 1)
    have h2 := two_mul_tri r
    simp only [Nat.add_sub_cancel]
    nlinarith

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma thetaSum_cong (k n : ℕ) : (X : K⟦X⟧) ^ (n + 1) ∣ (thetaSum k : K⟦X⟧) - ThetaN k n := by
  rw [X_pow_dvd_iff]
  intro m hm
  rw [map_sub, sub_eq_zero, thetaSum, coeff_mk, ThetaN, map_sum]
  have hsub : range (m + 1) ⊆ range (n + 1) := range_subset_range.mpr hm
  rw [← sum_subset hsub]
  · refine sum_congr rfl (fun r _ => ?_)
    congr 1
    rw [theta_exp]
    ring
  · intro r _ hr
    have hr' : m + 1 ≤ r := by simpa using hr
    refine X_pow_dvd_iff.mp (Dvd.dvd.mul_right (pow_dvd_pow X ?_) _) m (Nat.lt_succ_self m)
    have := le_tri r
    nlinarith

/-- **Proposition 6 at `x = 1`** (with the bilateral sum folded):
`(q;q²)_∞ (q;q)_∞ Mel_k(1) = Σ_{n ≥ 0} (-1)^n q^{((k+2)n²+kn)/2} (1 - q^{2n+1})`. -/
theorem mel_theta (k : ℕ) :
    (qPochInf 1 2 * qPochInf 1 1 * Mel k 1 : K⟦X⟧) = thetaSum k := by
  apply eq_of_forall_X_pow_dvd_sub
  intro n
  have hR := limit_one (K := K) k n
  have hq1 : (X : K⟦X⟧) ^ (n + 1) ∣ qPochInf 1 1 - poch X n := by
    have h := qPochInf_trunc (K := K) 1 1 n le_rfl le_rfl
    have e : ∏ i ∈ range (n + 1), (1 - X ^ (1 + 1 * i) : K⟦X⟧) = poch X n * (1 - X ^ (n + 1)) := by
      rw [← poch_succ]
      simp only [poch]
      refine prod_congr rfl (fun i _ => ?_)
      ring_nf
    rw [e] at h
    have e2 : (qPochInf 1 1 : K⟦X⟧) - poch X n
        = (qPochInf 1 1 - poch X n * (1 - X ^ (n + 1))) - X ^ (n + 1) * poch X n := by ring
    rw [e2]
    exact dvd_sub h (dvd_mul_right _ _)
  have hq2 : (X : K⟦X⟧) ^ (n + 1) ∣ qPochInf 1 2 - oddP n := by
    have h := qPochInf_trunc (K := K) 1 2 n le_rfl (by norm_num)
    have e : ∏ i ∈ range (n + 1), (1 - X ^ (1 + 2 * i) : K⟦X⟧) = oddP n * (1 - X ^ (1 + 2 * n)) := by
      rw [oddP, prod_range_succ]
    rw [e] at h
    have e2 : (qPochInf 1 2 : K⟦X⟧) - oddP n
        = (qPochInf 1 2 - oddP n * (1 - X ^ (1 + 2 * n))) - X ^ (1 + 2 * n) * oddP n := by ring
    rw [e2]
    exact dvd_sub h (Dvd.dvd.mul_right (pow_dvd_pow X (by omega)) _)
  have s1 : (X : K⟦X⟧) ^ (n + 1) ∣ qPochInf 1 2 * qPochInf 1 1 * Mel k 1
      - oddP n * poch X n * Mel k 1 := by
    rw [← sub_mul]
    exact (cong_mul hq2 hq1).mul_right _
  have s2 : (X : K⟦X⟧) ^ (n + 1) ∣ oddP n * poch X n * Mel k 1
      - oddP n * (gp (-X) n * ThetaN k n) := by
    rw [mul_assoc, ← mul_sub]
    exact hR.mul_left _
  have s3 : (X : K⟦X⟧) ^ (n + 1) ∣ oddP n * (gp (-X) n * ThetaN k n) - ThetaN k n := by
    have e : (oddP n : K⟦X⟧) * (gp (-X) n * ThetaN k n) - ThetaN k n
        = (gp (-X) n * oddP n - 1) * ThetaN k n := by ring
    rw [e, gpneg_mul_oddP]
    refine Dvd.dvd.mul_right (X_pow_dvd_prod_sub_one _ _ _ (fun i _ => ?_)) _
    have : (1 : K⟦X⟧) - (1 - X ^ (n + i + 1)) = X ^ (n + i + 1) := by ring
    rw [this]
    exact pow_dvd_pow X (by omega)
  have s4 : (X : K⟦X⟧) ^ (n + 1) ∣ ThetaN k n - thetaSum k :=
    (dvd_sub_comm).mp (thetaSum_cong k n)
  exact (pow_dvd_pow X (Nat.le_succ n)).trans
    (dvd_sub_trans (dvd_sub_trans (dvd_sub_trans s1 s2) s3) s4)

/-- The KRTW sum in theta form:
`(2k-fold KRTW sum) · (q;q²)_∞ (q;q)_∞ = Σ_{n ∈ ℤ} (-1)^n q^{((k+2)n²+kn)/2}`. -/
theorem krtw_theta' (k : ℕ) :
    krtwSum (1 : K⟦X⟧) k * (qPochInf 1 2 * qPochInf 1 1) = thetaSum k := by
  rw [krtwSum_eq_Kf, Kf_eq_Mel, ← mel_theta]
  ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma constantCoeff_qPochInf (a d : ℕ) (ha : 1 ≤ a) :
    constantCoeff (qPochInf a d : K⟦X⟧) = 1 := by
  obtain ⟨a', rfl⟩ := Nat.exists_eq_add_of_le ha
  rw [← coeff_zero_eq_constantCoeff_apply, qPochInf, coeff_mk]
  simp [add_comm 1 a']

/-- **KRTW identity (i)** for every `k`:
`(2k-fold KRTW sum) = (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ / ((q;q²)_∞ (q;q)_∞)`. -/
theorem krtw_product' (k : ℕ) :
    krtwSum (1 : K⟦X⟧) k = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
      * (qPochInf 1 2 * qPochInf 1 1)⁻¹ := by
  have h := krtw_theta' (K := K) k
  rw [← jtp_theta] at h
  have hne : constantCoeff (qPochInf 1 2 * qPochInf 1 1 : K⟦X⟧) ≠ 0 := by
    rw [map_mul, constantCoeff_qPochInf 1 2 le_rfl, constantCoeff_qPochInf 1 1 le_rfl]
    simp
  rw [← h, mul_assoc, PowerSeries.mul_inv_cancel _ hne, mul_one]

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- Stability of the truncated products: `qPochInf a d` has the coefficients of the infinite product. -/
theorem coeff_qPochInf' (a d : ℕ) (ha : 1 ≤ a) (hd : 1 ≤ d) (N M : ℕ) (h : N ≤ M) :
    coeff N (qPochInf a d : K⟦X⟧) = coeff N (∏ i ∈ range (M + 1), (1 - X ^ (a + d * i) : K⟦X⟧)) := by
  rw [qPochInf, coeff_mk]
  have hM : M + 1 = (N + 1) + (M - N) := by omega
  rw [hM]
  symm
  apply coeff_eq_of_dvd_sub
  apply prod_trunc
  intro i
  have : (1 : K⟦X⟧) - (1 - X ^ (a + d * (N + 1 + i))) = X ^ (a + d * (N + 1 + i)) := by ring
  rw [this]
  apply pow_dvd_pow
  have : N + 1 + i ≤ d * (N + 1 + i) := Nat.le_mul_of_pos_left _ hd
  omega

end KRTW
