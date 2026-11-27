/-
Nahm/JTP.lean — the Jacobi triple product in the general form needed for the Andrews–Gordon
products, and the dissection of `(q;q)_∞` into residue classes.

  `(q^c; q^{c+i})_∞ (q^i; q^{c+i})_∞ (q^{c+i}; q^{c+i})_∞
      = Σ_{n ≥ 0} (-1)^n q^{c n(n+1)/2 + i n(n-1)/2} (1 - q^{i(2n+1)})`      (`c, i ≥ 1`)

(from the finite triple product `KRTW.limit_tp`), and `(q;q)_∞ = ∏_{a=1}^{D} (q^a; q^D)_∞`.
-/
import Nahm.AGTheta
import Nahm.Defs
import KRTW.JTP
import KRTW.Stage2

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K]

lemma term_gen (c i m : ℕ) :
    ((-X ^ c) ^ (tri m) * (-X ^ i) ^ (tri m - m) * (1 + (-X ^ i) ^ (2 * m + 1)) : K⟦X⟧)
      = (-1) ^ m * X ^ (c * tri m + i * tri (m - 1)) * (1 - X ^ (i * (2 * m + 1))) := by
  have hs := tri_pred m
  rw [show tri m - m = tri (m - 1) by omega, hs]
  generalize tri (m - 1) = s
  rw [neg_pow, neg_pow (X ^ i : K⟦X⟧) s, Odd.neg_pow ⟨m, rfl⟩ (X ^ i : K⟦X⟧)]
  have hsg : ((-1 : K⟦X⟧) ^ s) ^ 2 = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  rw [← pow_mul, ← pow_mul, ← pow_mul]
  linear_combination ((-1 : K⟦X⟧) ^ m * X ^ (c * (s + m) + i * s) * (1 - X ^ (i * (2 * m + 1)))) * hsg

lemma summable_thetaP_term [TopologicalSpace K] [DiscreteTopology K] (c i : ℕ) (hc : 1 ≤ c) :
    Summable fun n : ℕ => ((-1) ^ n * X ^ (c * tri n + i * tri (n - 1)) * (1 - X ^ (i * (2 * n + 1))) : K⟦X⟧) := by
  apply summable_nat
  intro n
  refine Dvd.dvd.mul_right (Dvd.dvd.mul_left (pow_dvd_pow X ?_) _) _
  have := le_tri n
  nlinarith

variable [TopologicalSpace K] [DiscreteTopology K]

/-- **The Jacobi triple product** (the instances `(q^c, q^i, q^{c+i}; q^{c+i})_∞`). -/
theorem jtp_gen (c i : ℕ) (hc : 1 ≤ c) (hi : 1 ≤ i) :
    (KRTW.qPochInf c (c + i) * KRTW.qPochInf i (c + i) * KRTW.qPochInf (c + i) (c + i) : K⟦X⟧)
      = thetaP X c i := by
  have hD : Nil (X ^ (c + i) : K⟦X⟧) := ⟨by rw [map_pow, constantCoeff_X, zero_pow (by omega)]⟩
  ext N
  have hlim := limit_tp (-X ^ c) (-X ^ i) (X ^ (c + i) : K⟦X⟧) (by ring)
    (dvd_neg.mpr (dvd_pow_self _ (by omega))) (neg_ne_zero.mpr (pow_ne_zero _ X_ne_zero))
    (pow_ne_zero _ X_ne_zero) N
  have h1 := qPochInf_trunc (K := K) i (c + i) N hi (by omega)
  have h2 := qPochInf_trunc (K := K) c (c + i) N hc (by omega)
  have h3 := qPochInf_trunc (K := K) (c + i) (c + i) N (by omega) (by omega)
  have e1 : ∏ j ∈ range (N + 1), (1 - X ^ (i + (c + i) * j) : K⟦X⟧)
      = ∏ j ∈ range (N + 1), (1 + (-X ^ i) * (X ^ (c + i)) ^ j) :=
    prod_congr rfl fun j _ => by ring
  have e3 : ∏ j ∈ range (N + 1), (1 - X ^ (c + i + (c + i) * j) : K⟦X⟧) = poch (X ^ (c + i)) (N + 1) :=
    prod_congr rfl fun j _ => by ring
  have h2' : (X : K⟦X⟧) ^ (N + 1) ∣ ∏ j ∈ range (N + 1), (1 - X ^ (c + (c + i) * j) : K⟦X⟧)
      - ∏ j ∈ range N, (1 + (-X ^ c) * (X ^ (c + i)) ^ j) := by
    have := prod_trunc (fun j => (1 - X ^ (c + (c + i) * j) : K⟦X⟧)) (N + 1) N 1
      (fun j => by simp only [sub_sub_cancel]; exact pow_dvd_pow X (by nlinarith))
    rw [show ∏ j ∈ range N, (1 + (-X ^ c) * (X ^ (c + i)) ^ j : K⟦X⟧)
        = ∏ j ∈ range N, (1 - X ^ (c + (c + i) * j)) from prod_congr rfl fun j _ => by ring]
    exact this
  rw [e1] at h1
  rw [e3] at h3
  have h2'' := dvd_add h2 h2'
  rw [sub_add_sub_cancel] at h2''
  have h12 := dvd_mul_sub_mul' (dvd_mul_sub_mul' h2'' h1) h3
  have hfin := dvd_add h12 hlim
  rw [show ∀ A B C P1 P2 P3 S : K⟦X⟧, A * B * C - P1 * P2 * P3 + (P3 * (P2 * P1) - S)
      = A * B * C - S from fun _ _ _ _ _ _ _ => by ring] at hfin
  have hc' := X_pow_dvd_iff.mp hfin N (Nat.lt_succ_self N)
  rw [map_sub, sub_eq_zero] at hc'
  rw [hc', thetaP, coeff_tsum (summable_thetaP_term c i hc), map_sum]
  symm
  rw [tsum_eq_sum (s := range (N + 1))]
  · exact sum_congr rfl (fun m _ => by rw [term_gen])
  · intro m hm
    have : N < m := by simpa using hm
    refine X_pow_dvd_iff.mp (Dvd.dvd.mul_right (Dvd.dvd.mul_left (pow_dvd_pow X ?_) _) _) N
      (Nat.lt_succ_self N)
    have := le_tri m
    nlinarith

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma prod_cong {ι : Type*} (s : Finset ι) (f g : ι → K⟦X⟧) (n : ℕ)
    (h : ∀ a ∈ s, (X : K⟦X⟧) ^ n ∣ f a - g a) : (X : K⟦X⟧) ^ n ∣ ∏ a ∈ s, f a - ∏ a ∈ s, g a := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [prod_insert ha, prod_insert ha]
    exact cong_mul (h a (mem_insert_self a s)) (ih (fun b hb => h b (mem_insert_of_mem hb)))

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma prod_double (D : ℕ) : ∀ n : ℕ,
    ∏ a ∈ range D, ∏ j ∈ range n, (1 - X ^ (a + 1 + D * j) : K⟦X⟧)
      = ∏ t ∈ range (D * n), (1 - X ^ (t + 1))
  | 0 => by simp
  | n + 1 => by
    simp only [prod_range_succ (fun j => (1 - X ^ (_ + 1 + D * j) : K⟦X⟧)) n, prod_mul_distrib,
      prod_double D n, show D * (n + 1) = D * n + D by ring, prod_range_add]
    congr 1
    refine prod_congr rfl (fun a _ => ?_)
    ring_nf

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- `(q;q)_∞ = ∏_{a=1}^{D} (q^a; q^D)_∞`. -/
theorem qPochInf_dissect (D : ℕ) (hD : 1 ≤ D) :
    (KRTW.qPochInf 1 1 : K⟦X⟧) = ∏ a ∈ range D, KRTW.qPochInf (a + 1) D := by
  apply eq_of_forall_X_pow_dvd_sub
  intro n
  rcases n with _ | n
  · simp
  have h1 := qPochInf_trunc (K := K) 1 1 n le_rfl le_rfl
  have h2 := prod_cong (range D) (fun a => (KRTW.qPochInf (a + 1) D : K⟦X⟧))
    (fun a => ∏ j ∈ range (n + 1), (1 - X ^ (a + 1 + D * j))) (n + 1)
    (fun a _ => qPochInf_trunc (a + 1) D n (by omega) hD)
  rw [prod_double] at h2
  have h3 := prod_trunc (fun t => (1 - X ^ (t + 1) : K⟦X⟧)) (n + 1) (n + 1) (D * (n + 1) - (n + 1))
    (fun t => by simp only [sub_sub_cancel]; exact pow_dvd_pow X (by omega))
  rw [show n + 1 + (D * (n + 1) - (n + 1)) = D * (n + 1) by
    have : n + 1 ≤ D * (n + 1) := Nat.le_mul_of_pos_left _ hD
    omega] at h3
  have e : ∏ j ∈ range (n + 1), (1 - X ^ (1 + 1 * j) : K⟦X⟧) = ∏ t ∈ range (n + 1), (1 - X ^ (t + 1)) :=
    prod_congr rfl (fun t _ => by ring_nf)
  rw [e] at h1
  have := dvd_sub_trans (dvd_sub_trans h1 ((dvd_sub_comm).mp h3)) ((dvd_sub_comm).mp h2)
  exact this

end Nahm
