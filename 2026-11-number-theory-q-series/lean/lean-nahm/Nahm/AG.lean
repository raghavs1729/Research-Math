/-
Nahm/AG.lean — the Andrews–Gordon identities in the printed product form, for the statement-level
definitions `agSum`, `agProduct` of `Nahm/Defs.lean`.
-/
import Nahm.JTP

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K]

lemma qPochInf_eq (a d : ℕ) : (Nahm.qPochInf a d : K⟦X⟧) = KRTW.qPochInf a d := rfl

lemma pinf_X : pinf (X : K⟦X⟧) = KRTW.qPochInf 1 1 := by
  ext N
  rw [pinf, coeff_mk, KRTW.qPochInf, coeff_mk]
  have e : ∏ i ∈ range (N + 1), (1 - X ^ (1 + 1 * i) : K⟦X⟧) = poch X (N + 1) := by
    simp only [poch]; exact prod_congr rfl (fun i _ => by ring_nf)
  rw [e]
  have := X_pow_dvd_iff.mp (poch_trunc (X : K⟦X⟧) (Nat.le_succ N)) N (Nat.lt_succ_self N)
  rw [map_sub, sub_eq_zero] at this
  exact this.symm

lemma constantCoeff_qPochInf' (a d : ℕ) (ha : 1 ≤ a) : constantCoeff (KRTW.qPochInf a d : K⟦X⟧) = 1 :=
  constantCoeff_qPochInf a d ha

lemma qexp_agb {r : ℕ} (i : ℕ) (n : Fin r → ℕ) :
    qexp (agb i) n = ∑ j, psum n j ^ 2 + ∑ j ∈ univ.filter (fun j : Fin r => i ≤ (j : ℕ) + 1), psum n j := by
  rw [qexp, sum_add_distrib, sum_filter]
  congr 1
  refine sum_congr rfl (fun j _ => ?_)
  simp only [agb]
  split_ifs <;> simp

lemma prod_range_succ_eq_Icc (f : ℕ → K⟦X⟧) (D : ℕ) : ∏ a ∈ range D, f (a + 1) = ∏ a ∈ Icc 1 D, f a := by
  have : Icc 1 D = (range D).map (addRightEmbedding 1) := by
    ext a
    simp only [mem_Icc, mem_map, mem_range, addRightEmbedding_apply]
    constructor
    · intro h; exact ⟨a - 1, by omega, by omega⟩
    · rintro ⟨b, hb, rfl⟩; omega
  rw [this, prod_map]
  rfl

/-- Assembling the product side: `(q^c, q^i, q^D; q^D)_∞ / (q;q)_∞ = 1/∏_{a ≠ c, i} (q^a;q^D)_∞`
for `D = 2k+1`, `c = D - i`. -/
lemma product_side (k i : ℕ) (hi : 1 ≤ i) (hik : i ≤ k) :
    (KRTW.qPochInf 1 1 : K⟦X⟧)⁻¹ * (KRTW.qPochInf (2 * k + 1 - i) (2 * k + 1) *
      KRTW.qPochInf i (2 * k + 1) * KRTW.qPochInf (2 * k + 1) (2 * k + 1)) = agProduct k i := by
  set D := 2 * k + 1 with hD
  set c := 2 * k + 1 - i with hc
  set R := ∏ a ∈ (Icc 1 (2 * k)).filter (fun a => a ≠ i ∧ a ≠ 2 * k + 1 - i), (KRTW.qPochInf a D : K⟦X⟧)
    with hR
  have hRc : constantCoeff R ≠ 0 := by
    rw [hR, map_prod]
    rw [prod_eq_one (fun a ha => by
      simp only [mem_filter, mem_Icc] at ha
      exact constantCoeff_qPochInf' a D ha.1.1)]
    exact one_ne_zero
  have hP : (KRTW.qPochInf 1 1 : K⟦X⟧)
      = R * (KRTW.qPochInf c D * KRTW.qPochInf i D * KRTW.qPochInf D D) := by
    have hI : Icc 1 D = insert D (Icc 1 (2 * k)) := by
      ext a; simp only [mem_Icc, mem_insert, hD]; omega
    rw [qPochInf_dissect D (by omega), prod_range_succ_eq_Icc (fun a => (KRTW.qPochInf a D : K⟦X⟧)) D,
      hI, prod_insert (by simp [hD])]
    rw [← prod_filter_mul_prod_filter_not (Icc 1 (2 * k)) (fun a => a ≠ i ∧ a ≠ 2 * k + 1 - i)]
    have hs : (Icc 1 (2 * k)).filter (fun a => ¬(a ≠ i ∧ a ≠ 2 * k + 1 - i)) = {i, c} := by
      ext a
      simp only [mem_filter, mem_Icc, mem_insert, mem_singleton, hc]
      constructor
      · intro h; by_contra h'; push Not at h'; exact h.2 ⟨h'.1, h'.2⟩
      · rintro (h | h) <;> refine ⟨⟨by omega, by omega⟩, ?_⟩ <;> omega
    rw [hs, prod_pair (by omega), ← hR]
    ring
  show (KRTW.qPochInf 1 1 : K⟦X⟧)⁻¹ * _ = R⁻¹
  rw [hP]
  have hT : constantCoeff (KRTW.qPochInf c D * KRTW.qPochInf i D * KRTW.qPochInf D D : K⟦X⟧) ≠ 0 := by
    rw [map_mul, map_mul, constantCoeff_qPochInf' c D (by omega), constantCoeff_qPochInf' i D hi,
      constantCoeff_qPochInf' D D (by omega)]
    simp
  generalize (KRTW.qPochInf c D * KRTW.qPochInf i D * KRTW.qPochInf D D : K⟦X⟧) = T at hT ⊢
  rw [PowerSeries.mul_inv_rev, mul_comm T⁻¹, mul_assoc, PowerSeries.inv_mul_cancel _ hT, mul_one]

variable [TopologicalSpace K] [DiscreteTopology K]

lemma agSum_eq_tsum (k i : ℕ) :
    (agSum k i : K⟦X⟧) = ∑' n : Fin (k - 1) → ℕ, mterm X (agb i) n := by
  ext N
  rw [agSum, coeff_mk, coeff_tsum_mterm, map_sum]
  refine sum_congr rfl (fun n _ => ?_)
  rw [mterm, qexp_agb]
  split_ifs with h
  · rfl
  · rw [map_zero]
    have hlt : N < qexp (agb i) n := lt_of_lt_of_le (by omega) (tot_le_qexp _ _)
    rw [qexp_agb] at hlt
    exact X_pow_dvd_iff.mp (Dvd.dvd.mul_right (dvd_refl _) _) N hlt

/-- **The Andrews–Gordon identities** (over an arbitrary field, `1 ≤ i ≤ k`). -/
theorem andrews_gordon_top (k i : ℕ) (hi : 1 ≤ i) (hik : i ≤ k) :
    (agSum k i : K⟦X⟧) = agProduct k i := by
  rcases Nat.lt_or_ge k 2 with hk | hk
  · -- k = 1, i = 1: both sides are 1
    obtain rfl : k = 1 := by omega
    obtain rfl : i = 1 := by omega
    have h1 : (agSum 1 1 : K⟦X⟧) = 1 := by
      ext N
      rw [agSum, coeff_mk]
      have : box (1 - 1) N = {Fin.elim0} := by
        ext r
        simp only [box, Fintype.mem_piFinset, mem_singleton]
        exact ⟨fun _ => funext fun j => absurd j.2 (by omega), fun _ j => absurd j.2 (by omega)⟩
      rw [this, sum_singleton]
      simp [agExp]
    have h2 : (agProduct 1 1 : K⟦X⟧) = 1 := by
      rw [agProduct]
      have : (Icc 1 (2 * 1)).filter (fun a => a ≠ 1 ∧ a ≠ 2 * 1 + 1 - 1) = ∅ := by decide
      rw [this, prod_empty, inv_one]
    rw [h1, h2]
  obtain ⟨s, rfl⟩ : ∃ s, k = s + 2 := ⟨k - 2, by omega⟩
  rw [agSum_eq_tsum]
  change ∑' n : Fin (s + 1) → ℕ, mterm X (agb i) n = _
  rw [← product_side (s + 2) i hi hik, ← pinf_X]
  rcases Nat.lt_or_ge i 2 with hi2 | hi2
  · obtain rfl : i = 1 := by omega
    rw [ag_theta_one, ← jtp_gen (2 * (s + 2)) 1 (by omega) le_rfl,
      show 2 * (s + 2) + 1 - 1 = 2 * (s + 2) by omega]
  · obtain ⟨p, rfl⟩ : ∃ p, i = p + 2 := ⟨i - 2, by omega⟩
    obtain ⟨q, rfl⟩ : ∃ q, s = q + p := ⟨s - p, by omega⟩
    rw [ag_theta_two, ← jtp_gen (p + 2 * q + 3) (p + 2) (by omega) (by omega),
      show 2 * (q + p + 2) + 1 - (p + 2) = p + 2 * q + 3 by omega,
      show p + 2 * q + 3 + (p + 2) = 2 * (q + p + 2) + 1 by ring]

/-- **The Andrews–Gordon identities** for modulus `2k+1`, `1 ≤ i ≤ k`:
`Σ_{n_1,…,n_{k-1} ≥ 0} q^{N_1²+⋯+N_{k-1}² + N_i+⋯+N_{k-1}} / ((q)_{n_1}⋯(q)_{n_{k-1}})
   = ∏_{n ≥ 1, n ≢ 0, ±i (mod 2k+1)} (1 - q^n)⁻¹`. -/
theorem andrews_gordon' (k i : ℕ) (hi : 1 ≤ i) (hik : i ≤ k) :
    (agSum k i : ℚ⟦X⟧) = agProduct k i := by
  let _ : TopologicalSpace ℚ := ⊥
  have _ : DiscreteTopology ℚ := ⟨rfl⟩
  exact andrews_gordon_top k i hi hik

end Nahm
