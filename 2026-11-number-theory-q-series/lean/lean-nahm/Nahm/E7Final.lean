/-
Nahm/E7Final.lean — Theorem 4 of `proofs/E7.md`:

  `Σ_{n ∈ ℕ^7} q^{½ nᵀ C(E7)⁻¹ n} / ∏ (q)_{n_i} = (1/(q;q)_∞) (Σ_{k ≥ 0} q^{3k²/4} + Σ_{k ≥ 1} q^{7k²/4})
                                              = (1/(2 (q;q)_∞)) Σ_{k ∈ ℤ} (q^{3k²/4} + q^{7k²/4})`,

as power series in `X = q^{1/4}`.  §1 (fusion) is `Nahm/E7.lean`, §2 (Lemma 4, Proposition 5) is
`Nahm/E7Lemmas.lean`, Proposition 6 is `Nahm/Prop6.lean`; here: §3 (collecting `n = a + 2r`).
-/
import Nahm.E7
import Nahm.Prop6

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K]

local notation "Q4" => ((X : K⟦X⟧) ^ 4)

lemma iq4_eq (n : ℕ) : iq Q4 n = (qPoch4 n : K⟦X⟧)⁻¹ := by
  rw [iq, poch, qPoch4]
  congr 1
  exact prod_congr rfl (fun i _ => by rw [← pow_mul])

lemma pinf_X4 : pinf Q4 = (qPochInf 4 4 : K⟦X⟧) := by
  ext N
  rw [pinf, coeff_mk, qPochInf, coeff_mk]
  have e : ∏ i ∈ range (N + 1), (1 - X ^ (4 + 4 * i) : K⟦X⟧) = poch Q4 (N + 1) := by
    simp only [poch]; exact prod_congr rfl (fun i _ => by rw [← pow_mul]; ring_nf)
  rw [e]
  have := X_pow_dvd_iff.mp (poch_trunc Q4 (Nat.le_succ N)) N (Nat.lt_succ_self N)
  rw [map_sub, sub_eq_zero] at this
  exact this.symm

/-! ### Integrality and the exponent -/

lemma le_E7e0 : ∀ (v : Fin 4 → ℕ) (i : Fin 4), v i ≤ E7e0 v := by
  intro v i
  refine Fin.cases ?_ (fun i => ?_) i
  · simp only [E7e0]; nlinarith
  · have h1 : (Fin.tail v) i ≤ qexp ![0, v 0, v 0] (Fin.tail v) := le_trans (le_tot _ i) (tot_le_qexp _ _)
    simp only [Fin.tail] at h1
    simp only [E7e0]; omega

lemma le_E7e3 : ∀ n i, n i ≤ E7e3 n :=
  le_fuseE 4 _ (le_fuseE 4 _ (le_fuseE 4 _ le_E7e0 _ _) _ _) _ _

noncomputable def e7Equiv : I3 ≃ Fin 7 := Equiv.ofBijective toFin7 toFin7_bijective

lemma e7Exp_eq (n : Fin 7 → ℕ) : e7Exp n = E7e3 (fun j => n (toFin7 j)) := by
  rw [e7Exp, ← E7e3_eq, Nat.floor_natCast]

/-- **Integrality**: `4 · ½ nᵀ C(E7)⁻¹ n` is a natural number (namely `e7Exp n`). -/
theorem e7Form_eq_e7Exp' (n : Fin 7 → ℕ) : (e7Exp n : ℚ) = 4 * e7Form n := by
  rw [e7Exp_eq, E7e3_eq]

lemma le_e7Exp (n : Fin 7 → ℕ) (i : Fin 7) : n i ≤ e7Exp n := by
  rw [e7Exp_eq]
  obtain ⟨j, rfl⟩ := toFin7_bijective.2 i
  exact le_E7e3 (fun j => n (toFin7 j)) j

variable [TopologicalSpace K] [DiscreteTopology K]

/-- The E7 sum is the rank 4 sum (three fusions). -/
theorem e7Sum_eq_rank4 : (e7Sum : K⟦X⟧) = ∑' v : Fin 4 → ℕ, nterm Q4 E7e0 v := by
  have hQ : Q4 = X ^ 4 := rfl
  have h1 := le_fuseE 4 _ le_E7e0 (0 : Fin 4) 1
  have h2 := le_fuseE 4 _ h1 (some (3 : Fin 4)) (some 1)
  have f3 := fusion (Q := Q4) 4 hQ E7e2 h2 (some (some 0)) (some (some 2)) (by decide)
  have f2 := fusion (Q := Q4) 4 hQ E7e1 h1 (some 3) (some 1) (by decide)
  have f1 := fusion (Q := Q4) 4 hQ E7e0 le_E7e0 0 1 (by decide)
  have hchain : ∑' n : I3 → ℕ, nterm Q4 E7e3 n = ∑' v : Fin 4 → ℕ, nterm Q4 E7e0 v := by
    rw [E7e3, f3, E7e2, f2, E7e1, f1]
  rw [← hchain, e7Sum]
  rw [mk_box_eq_tsum (fun n => X ^ e7Exp n * ∏ i, (qPoch4 (n i))⁻¹) e7Exp le_e7Exp
    (fun n => Dvd.dvd.mul_right (dvd_refl _) _)]
  let φ : (Fin 7 → ℕ) ≃ (I3 → ℕ) := Equiv.arrowCongr e7Equiv.symm (Equiv.refl ℕ)
  rw [← φ.tsum_eq]
  refine tsum_congr (fun n => ?_)
  simp only [nterm, φ, Equiv.arrowCongr_apply, Equiv.symm_symm, Equiv.coe_refl, Function.comp_apply, id]
  rw [e7Exp_eq]
  congr 1
  exact (Fintype.prod_equiv e7Equiv (fun j => iq Q4 (n (e7Equiv j))) (fun i => (qPoch4 (n i))⁻¹)
    (fun j => iq4_eq _)).symm

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma r4_dvd (a : ℕ) (m : Fin 3 → ℕ) :
    (X : K⟦X⟧) ^ (a + tot m) ∣ X ^ (3 * a ^ 2) * iq Q4 a * mterm Q4 ![0, a, a] m := by
  simp only [mterm]
  have h1 := tot_le_qexp ![0, a, a] m
  rw [← pow_mul, show X ^ (3 * a ^ 2) * iq Q4 a * (X ^ (4 * qexp ![0, a, a] m) * ∏ l, iq Q4 (m l))
    = X ^ (3 * a ^ 2 + 4 * qexp ![0, a, a] m) * (iq Q4 a * ∏ l, iq Q4 (m l)) by rw [pow_add]; ring]
  exact Dvd.dvd.mul_right (pow_dvd_pow X (by nlinarith)) _

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma r4_nterm (a : ℕ) (m : Fin 3 → ℕ) :
    nterm Q4 E7e0 (Fin.cons a m) = X ^ (3 * a ^ 2) * iq Q4 a * mterm Q4 ![0, a, a] m := by
  simp only [nterm, E7e0, Fin.cons_zero, Fin.tail_cons, mterm, Fin.prod_univ_succ, Fin.cons_succ]
  rw [pow_add, ← pow_mul]
  ring

lemma r4_inner (a : ℕ) : ∑' m : Fin 3 → ℕ, X ^ (3 * a ^ 2) * iq Q4 a * mterm Q4 ![0, a, a] m
    = X ^ (3 * a ^ 2) * iq Q4 a * ∑' m : Fin 3 → ℕ, mterm Q4 ![0, a, a] m :=
  Summable.tsum_mul_left _ (summable_mterm _ _)

/-- Peeling off the variable `a`: the rank 4 sum is `Σ_a q^{3a²/4} T(a) / (q)_a`. -/
theorem rank4_eq : ∑' v : Fin 4 → ℕ, nterm Q4 E7e0 v
    = ∑' a, X ^ (3 * a ^ 2) * iq Q4 a * ∑' m : Fin 3 → ℕ, mterm Q4 ![0, a, a] m := by
  have hfin : ∀ N, {p : ℕ × (Fin 3 → ℕ) | p.1 + tot p.2 ≤ N}.Finite := by
    intro N
    apply ((range (N + 1)) ×ˢ (KRTW.box 3 N)).finite_toSet.subset
    rintro ⟨a, m⟩ h
    simp only [Set.mem_ofPred_eq] at h
    simp only [coe_product, Set.mem_prod, mem_coe, mem_range]
    exact ⟨by omega, mem_box_of_tot_le (by omega)⟩
  have hs : Summable fun p : ℕ × (Fin 3 → ℕ) => X ^ (3 * p.1 ^ 2) * iq Q4 p.1 * mterm Q4 ![0, p.1, p.1] p.2 :=
    summable_of_dvd _ hfin (fun p => r4_dvd p.1 p.2)
  have hs1 : ∀ a, Summable fun m => X ^ (3 * a ^ 2) * iq Q4 a * mterm Q4 ![0, a, a] m := by
    intro a
    apply summable_of_dvd (fun m => tot m) _ (fun m => (pow_dvd_pow X (Nat.le_add_left _ _)).trans (r4_dvd a m))
    intro N
    exact (KRTW.box 3 N).finite_toSet.subset (fun m hm => mem_box_of_tot_le hm)
  rw [← (Fin.consEquiv (fun _ : Fin 4 => ℕ)).tsum_eq (nterm Q4 E7e0)]
  have e : ∀ p : ℕ × (Fin 3 → ℕ), nterm Q4 E7e0 (Fin.consEquiv (fun _ : Fin 4 => ℕ) p)
      = X ^ (3 * p.1 ^ 2) * iq Q4 p.1 * mterm Q4 ![0, p.1, p.1] p.2 := fun p => r4_nterm p.1 p.2
  rw [tsum_congr e, hs.tsum_prod' hs1]
  exact tsum_congr fun a => r4_inner a

/-! ### §3: collecting the terms with `a + 2r = n` -/

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma ua_qbinom (Q : K⟦X⟧) [Nil Q] (a r j : ℕ) (hj : j ≤ a) (h : 1 ≤ a + 2 * r) :
    ua Q (Q ^ a) r * qbinom Q a j
      = (-1) ^ r * Q ^ tri (r - 1) * (1 - Q ^ (a + 2 * r)) * poch Q (a + r - 1)
          * iq Q r * iq Q j * iq Q (a - j) := by
  have hb : qbinom Q a j = poch Q a * iq Q j * iq Q (a - j) := by simp [qbinom, hj]
  rw [hb]
  cases r with
  | zero =>
    obtain ⟨b, rfl⟩ : ∃ b, a = b + 1 := ⟨a - 1, by omega⟩
    rw [show b + 1 + 2 * 0 = b + 1 by ring, show b + 1 + 0 - 1 = b by omega, show (0 : ℕ) - 1 = 0 from rfl]
    simp only [ua, tri_zero, pow_zero, mul_one, one_mul, iq_zero]
    rw [poch_succ]
    ring
  | succ s =>
    have hg := gq_pow Q a 0 s
    rw [pow_one, add_zero] at hg
    simp only [ua, Nat.add_sub_cancel]
    have hi := poch_mul_iq Q a
    have e1 : Q ^ a * Q ^ (2 * s + 2) = Q ^ (a + 2 * (s + 1)) := by rw [← pow_add]; ring_nf
    rw [e1, show a + (s + 1) - 1 = a + s by omega, ← hg]
    ring

/-- The summand of `(q;q)_∞ · F` after inserting Proposition 5. -/
noncomputable def H7 (a r : ℕ) : K⟦X⟧ :=
  X ^ (3 * a ^ 2) * (Q4 ^ (r ^ 2) * alE7 Q4 a r *
    ∑ j ∈ range (a + 1), qbinom Q4 a j * Q4 ^ (j ^ 2 + 2 * r * j))

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma X_pow_dvd_H7 (a r : ℕ) : (X : K⟦X⟧) ^ (3 * a ^ 2 + 4 * r ^ 2) ∣ H7 a r := by
  rw [H7, pow_add]
  refine mul_dvd_mul_left _ (Dvd.dvd.mul_right (Dvd.dvd.mul_right ?_ _) _)
  rw [← pow_mul]

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- The coefficient `c_n` of `q^{3n²/4}`: for `n ≥ 1`, by Proposition 6. -/
lemma sum_H7 (n : ℕ) :
    ∑ r ∈ range (n + 1), (if 2 * r ≤ n then H7 (n - 2 * r) r else 0)
      = X ^ (3 * n ^ 2) * (1 + if n = 0 then 0 else Q4 ^ (n ^ 2)) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp [H7, alE7, ua]
  rw [if_neg (by omega), ← prop6 Q4 n hn, mul_sum]
  refine sum_congr rfl (fun r _ => ?_)
  split_ifs with h
  · obtain ⟨a, rfl⟩ : ∃ a, n = a + 2 * r := ⟨n - 2 * r, by omega⟩
    rw [show a + 2 * r - 2 * r = a by omega, H7, mul_sum, mul_sum, mul_sum]
    rw [← sum_subset (range_subset_range.mpr (by omega : a + 1 ≤ a + 2 * r + 1))]
    · refine sum_congr rfl (fun j hj => ?_)
      have hj' : j ≤ a := Nat.lt_succ_iff.mp (mem_range.mp hj)
      rw [if_pos (by omega), alE7]
      have hm : ((Q4 ^ a) ^ r * Q4 ^ (r ^ 2)) ^ 3 * Q4 ^ (r ^ 2) * Q4 ^ (j ^ 2 + 2 * r * j) *
          X ^ (3 * a ^ 2) * Q4 ^ tri (r - 1)
          = X ^ (3 * (a + 2 * r) ^ 2) * Q4 ^ (r ^ 2 + tri (r - 1) + j ^ 2 + 2 * r * j) := by
        simp only [← pow_mul, ← pow_add]
        congr 1
        ring
      rw [show X ^ (3 * a ^ 2) * (Q4 ^ (r ^ 2) * (((Q4 ^ a) ^ r * Q4 ^ (r ^ 2)) ^ 3 * ua Q4 (Q4 ^ a) r) *
          (qbinom Q4 a j * Q4 ^ (j ^ 2 + 2 * r * j)))
          = ((Q4 ^ a) ^ r * Q4 ^ (r ^ 2)) ^ 3 * Q4 ^ (r ^ 2) * Q4 ^ (j ^ 2 + 2 * r * j) * X ^ (3 * a ^ 2) *
            (ua Q4 (Q4 ^ a) r * qbinom Q4 a j) by ring,
        ua_qbinom Q4 a r j hj' (by omega), show a + 2 * r - r - 1 = a + r - 1 by omega]
      rw [show ((Q4 ^ a) ^ r * Q4 ^ (r ^ 2)) ^ 3 * Q4 ^ (r ^ 2) * Q4 ^ (j ^ 2 + 2 * r * j) * X ^ (3 * a ^ 2) *
          ((-1) ^ r * Q4 ^ tri (r - 1) * (1 - Q4 ^ (a + 2 * r)) * poch Q4 (a + r - 1) * iq Q4 r * iq Q4 j *
            iq Q4 (a - j))
          = (((Q4 ^ a) ^ r * Q4 ^ (r ^ 2)) ^ 3 * Q4 ^ (r ^ 2) * Q4 ^ (j ^ 2 + 2 * r * j) * X ^ (3 * a ^ 2) *
            Q4 ^ tri (r - 1)) * ((-1) ^ r * (1 - Q4 ^ (a + 2 * r)) * poch Q4 (a + r - 1) * iq Q4 r * iq Q4 j *
            iq Q4 (a - j)) by ring, hm]
      ring
    · intro j _ hj
      rw [if_neg (by simp only [mem_range, not_lt] at hj; omega), mul_zero]
  · rw [eq_comm]
    refine mul_eq_zero_of_right _ (Finset.sum_eq_zero (fun j _ => ?_))
    rw [if_neg (by omega)]

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma X_pow_dvd_alE7_term (a r : ℕ) : (X : K⟦X⟧) ^ (4 * r ^ 2) ∣ Q4 ^ (r ^ 2) * alE7 Q4 a r *
    ∑ j ∈ range (a + 1), qbinom Q4 a j * Q4 ^ (j ^ 2 + 2 * r * j) := by
  refine Dvd.dvd.mul_right (Dvd.dvd.mul_right ?_ _) _
  rw [← pow_mul]

/-- `(q;q)_∞ F = Σ_n q^{3n²/4} c_n` with `c_0 = 1`, `c_n = 1 + q^{n²}`. -/
theorem e7Sum_tsum : (e7Sum : K⟦X⟧)
    = (pinf Q4)⁻¹ * ∑' n, X ^ (3 * n ^ 2) * (1 + if n = 0 then 0 else Q4 ^ (n ^ 2)) := by
  rw [e7Sum_eq_rank4, rank4_eq]
  have hg : ∀ a, Summable fun r => Q4 ^ (r ^ 2) * alE7 Q4 a r *
      ∑ j ∈ range (a + 1), qbinom Q4 a j * Q4 ^ (j ^ 2 + 2 * r * j) :=
    fun a => summable_nat (fun r => (pow_dvd_pow X (by nlinarith)).trans (X_pow_dvd_alE7_term a r))
  have e1 : ∀ a, X ^ (3 * a ^ 2) * iq Q4 a * ∑' m : Fin 3 → ℕ, mterm Q4 ![0, a, a] m
      = (pinf Q4)⁻¹ * ∑' r, H7 a r := by
    intro a
    rw [prop5, show ∑' r, (H7 a r : K⟦X⟧) = ∑' r, X ^ (3 * a ^ 2) * (Q4 ^ (r ^ 2) * alE7 Q4 a r *
      ∑ j ∈ range (a + 1), qbinom Q4 a j * Q4 ^ (j ^ 2 + 2 * r * j)) from rfl, (hg a).tsum_mul_left]
    have := poch_mul_iq Q4 a
    linear_combination (X ^ (3 * a ^ 2) * ((pinf Q4)⁻¹ * ∑' r, Q4 ^ (r ^ 2) * alE7 Q4 a r *
      ∑ j ∈ range (a + 1), qbinom Q4 a j * Q4 ^ (j ^ 2 + 2 * r * j))) * this
  rw [tsum_congr e1]
  have hdivH : ∀ a r, (X : K⟦X⟧) ^ (max a r) ∣ H7 a r := fun a r =>
    (pow_dvd_pow X (by
      rcases le_total a r with h | h
      · rw [max_eq_right h]; nlinarith
      · rw [max_eq_left h]; nlinarith)).trans (X_pow_dvd_H7 a r)
  have hsH : Summable (Function.uncurry (H7 : ℕ → ℕ → K⟦X⟧)) :=
    summable_nat2 (fun p => hdivH p.1 p.2)
  have hs1 : ∀ a, Summable ((H7 : ℕ → ℕ → K⟦X⟧) a) :=
    fun a => summable_nat (fun r => (pow_dvd_pow X (le_max_right a r)).trans (hdivH a r))
  have hs2 : ∀ r, Summable fun a => (H7 a r : K⟦X⟧) :=
    fun r => summable_nat (fun a => (pow_dvd_pow X (le_max_left a r)).trans (hdivH a r))
  have hsA : Summable fun a => ∑' r, (H7 a r : K⟦X⟧) :=
    summable_nat (fun a => X_pow_dvd_tsum (fun r => (pow_dvd_pow X (le_max_left a r)).trans (hdivH a r)))
  rw [hsA.tsum_mul_left, ← hsH.tsum_comm' hs1 hs2]
  congr 1
  set F : ℕ → ℕ → K⟦X⟧ := fun r M => if 2 * r ≤ M then H7 (M - 2 * r) r else 0 with hF
  have hFd : ∀ r M, (X : K⟦X⟧) ^ M ∣ F r M := by
    intro r M
    simp only [hF]
    split_ifs with h
    · exact (pow_dvd_pow X (by
        have : M - 2 * r ≤ 3 * (M - 2 * r) ^ 2 := by nlinarith
        have : 2 * r ≤ 4 * r ^ 2 := by nlinarith
        omega)).trans (X_pow_dvd_H7 (M - 2 * r) r)
    · exact dvd_zero _
  have hsh := tsum_tsum_shift (fun r => 2 * r) (fun r => by omega) F hFd
  have e2 : ∀ r, ∑' a, (H7 a r : K⟦X⟧) = ∑' m, F r (m + 2 * r) := by
    intro r
    refine tsum_congr (fun m => ?_)
    simp only [hF, if_pos (Nat.le_add_left _ _), Nat.add_sub_cancel]
  rw [tsum_congr e2, hsh]
  refine tsum_congr (fun M => ?_)
  rw [← sum_H7]
  refine sum_congr rfl (fun r _ => ?_)
  simp only [hF]
  split_ifs <;> rfl

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma e7Theta_term (n : ℕ) : (X : K⟦X⟧) ^ (3 * n ^ 2) * (1 + if n = 0 then 0 else Q4 ^ (n ^ 2))
    = X ^ (3 * n ^ 2) + if n = 0 then 0 else X ^ (7 * n ^ 2) := by
  split_ifs
  · ring
  · rw [← pow_mul, mul_add, mul_one, ← pow_add]
    congr 2
    ring

lemma e7ThetaOneSided_eq :
    (e7ThetaOneSided : K⟦X⟧) = ∑' n, X ^ (3 * n ^ 2) * (1 + if n = 0 then 0 else Q4 ^ (n ^ 2)) := by
  have hd : ∀ n, (X : K⟦X⟧) ^ n ∣ X ^ (3 * n ^ 2) * (1 + if n = 0 then 0 else Q4 ^ (n ^ 2)) :=
    fun n => Dvd.dvd.mul_right (pow_dvd_pow X (by nlinarith)) _
  ext N
  rw [e7ThetaOneSided, coeff_mk, coeff_tsum (summable_nat hd)]
  symm
  rw [tsum_eq_sum (s := range (N + 1))]
  · exact sum_congr rfl (fun n _ => by rw [e7Theta_term])
  · intro n hn
    have : N < n := by simpa using hn
    exact X_pow_dvd_iff.mp (hd n) N this

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma sum_Icc_natAbs (f : ℕ → K) : ∀ M : ℕ,
    ∑ k ∈ Icc (-(M : ℤ)) M, f k.natAbs = f 0 + 2 * ∑ m ∈ range M, f (m + 1)
  | 0 => by simp
  | M + 1 => by
    have hI : Icc (-((M + 1 : ℕ) : ℤ)) ((M + 1 : ℕ) : ℤ)
        = insert (-((M : ℤ) + 1)) (insert ((M : ℤ) + 1) (Icc (-(M : ℤ)) M)) := by
      ext k; simp only [mem_Icc, mem_insert]; push_cast; omega
    rw [hI, sum_insert (by simp only [mem_insert, mem_Icc]; omega),
      sum_insert (by simp only [mem_Icc]; omega), sum_Icc_natAbs f M, sum_range_succ]
    have h1 : (-((M : ℤ) + 1)).natAbs = M + 1 := by omega
    have h2 : ((M : ℤ) + 1).natAbs = M + 1 := by omega
    rw [h1, h2]
    ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma e7Theta_eq_two_mul : (e7Theta : K⟦X⟧) = 2 * e7ThetaOneSided := by
  rw [two_mul]
  ext N
  rw [map_add, e7Theta, e7ThetaOneSided, coeff_mk, coeff_mk,
    sum_Icc_natAbs (fun m => coeff N ((X : K⟦X⟧) ^ (3 * m ^ 2) + X ^ (7 * m ^ 2))) N,
    sum_range_succ']
  simp only [if_neg (Nat.succ_ne_zero _), if_true]
  simp only [pow_two, mul_zero, pow_zero, add_zero]
  rw [map_add]
  ring

/-- **Theorem 4 of `E7.md`** (one-sided form), over an arbitrary field. -/
theorem e7_identity_top : (e7Sum : K⟦X⟧) = (qPochInf 4 4)⁻¹ * e7ThetaOneSided := by
  rw [e7Sum_tsum, e7ThetaOneSided_eq, pinf_X4]

/-- **Theorem 4 of `E7.md`** (bilateral form), over an arbitrary field:
`2 F = (1/(q;q)_∞) Σ_{k ∈ ℤ} (q^{3k²/4} + q^{7k²/4})`. -/
theorem e7_identity_bilateral_top : 2 * (e7Sum : K⟦X⟧) = (qPochInf 4 4)⁻¹ * e7Theta := by
  rw [e7_identity_top, e7Theta_eq_two_mul]
  ring

end Nahm
