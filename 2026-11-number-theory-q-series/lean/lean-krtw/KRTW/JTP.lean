import KRTW.Basic
import KRTW.Defs

namespace KRTW
open PowerSeries Finset
variable {K : Type*} [Field K]

/-! ### The homogeneous `q`-binomial theorem -/

/-- The `t`-th term of the homogeneous `q`-binomial expansion of `∏_{i<N} (w + z Q^{i+1})`. -/
noncomputable def qbt (Q w z : K⟦X⟧) (N t : ℕ) : K⟦X⟧ :=
  qbinom Q N t * z ^ t * Q ^ (tri t) * w ^ (N - t)

lemma qbt_zero (Q w z : K⟦X⟧) [Nil Q] (N : ℕ) : qbt Q w z (N + 1) 0 = qbt Q w z N 0 * w := by
  simp [qbt, pow_succ]

lemma qbt_succ (Q w z : K⟦X⟧) [Nil Q] (N t : ℕ) :
    qbt Q w z (N + 1) (t + 1) = qbt Q w z N (t + 1) * w + qbt Q w z N t * (z * Q ^ (N + 1)) := by
  rcases Nat.lt_or_ge N t with h | h
  · simp [qbt, qbinom_of_lt Q h, qbinom_of_lt Q (show N + 1 < t + 1 by omega),
      qbinom_of_lt Q (show N < t + 1 by omega)]
  rcases Nat.eq_or_lt_of_le h with h' | h'
  · subst h'
    simp only [qbt, qbinom_self, qbinom_of_lt Q (Nat.lt_succ_self t), Nat.sub_self, tri_succ]
    ring
  · obtain ⟨r, rfl⟩ : ∃ r, N = t + 1 + r := ⟨N - t - 1, by omega⟩
    simp only [qbt]
    rw [qbinom_succ_succ, show t + 1 + r - t = r + 1 by omega,
      show t + 1 + r + 1 - (t + 1) = r + 1 by omega, show t + 1 + r - (t + 1) = r by omega,
      tri_succ]
    ring

lemma qbt_of_lt (Q w z : K⟦X⟧) {N t : ℕ} (h : N < t) : qbt Q w z N t = 0 := by
  simp [qbt, qbinom_of_lt Q h]

/-- Homogeneous `q`-binomial theorem: `∏_{i<N} (w + z Q^{i+1}) = Σ_t [N;t]_Q z^t Q^{tri t} w^{N-t}`. -/
lemma qbinomial (Q w z : K⟦X⟧) [Nil Q] (N : ℕ) :
    ∏ i ∈ range N, (w + z * Q ^ (i + 1)) = ∑ t ∈ range (N + 1), qbt Q w z N t := by
  induction N with
  | zero => simp [qbt]
  | succ N ih =>
    rw [prod_range_succ, ih, sum_range_succ' (qbt Q w z (N + 1))]
    simp_rw [qbt_succ, qbt_zero]
    rw [sum_add_distrib, ← sum_mul, ← sum_mul]
    have h1 : ∑ t ∈ range (N + 1), qbt Q w z N t
        = ∑ t ∈ range N, qbt Q w z N (t + 1) + qbt Q w z N 0 := sum_range_succ' _ _
    have h2 : ∑ t ∈ range (N + 1), qbt Q w z N (t + 1)
        = ∑ t ∈ range N, qbt Q w z N (t + 1) + qbt Q w z N (N + 1) := sum_range_succ _ _
    rw [h2, qbt_of_lt Q w z (Nat.lt_succ_self N), h1]
    ring

/-! ### The finite triple product -/

lemma qbinom_symm (Q : K⟦X⟧) {R r : ℕ} (h : r ≤ R) : qbinom Q R (R - r) = qbinom Q R r := by
  rw [qbinom, qbinom, if_pos (Nat.sub_le R r), if_pos h, Nat.sub_sub_self h]; ring

lemma prod_pow_tri (Q : K⟦X⟧) (n : ℕ) : ∏ i ∈ range n, Q ^ (i + 1) = Q ^ (tri n) := by
  induction n with
  | zero => simp
  | succ n ih => rw [prod_range_succ, ih, tri_succ, ← pow_add]

/-- Factoring the left side of the `q`-binomial theorem at `N = 2n+1`, `w = Q^{n+1}`, `z = b`. -/
lemma lhs_split (a b Q : K⟦X⟧) (hQ : Q = a * b) (n : ℕ) :
    ∏ i ∈ range (2 * n + 1), (Q ^ (n + 1) + b * Q ^ (i + 1))
      = (Q ^ (tri n) * b ^ n * (Q ^ (n + 1)) ^ (n + 1))
        * ((∏ i ∈ range (n + 1), (1 + b * Q ^ i)) * ∏ i ∈ range n, (1 + a * Q ^ i)) := by
  rw [show 2 * n + 1 = n + (n + 1) by ring, prod_range_add]
  have h1 : ∏ i ∈ range n, (Q ^ (n + 1) + b * Q ^ (i + 1))
      = ∏ i ∈ range n, (Q ^ (i + 1) * b * (1 + a * Q ^ (n - 1 - i))) := by
    refine prod_congr rfl (fun i hi => ?_)
    have hi' := mem_range.mp hi
    obtain ⟨r, hr⟩ : ∃ r, n = i + 1 + r := ⟨n - i - 1, by omega⟩
    rw [show n - 1 - i = r by omega, hr, hQ]; ring
  have h2 : ∏ i ∈ range (n + 1), (Q ^ (n + 1) + b * Q ^ (n + i + 1))
      = ∏ i ∈ range (n + 1), (Q ^ (n + 1) * (1 + b * Q ^ i)) :=
    prod_congr rfl (fun i _ => by ring)
  simp only [h1, h2, prod_mul_distrib]
  simp only [prod_pow_tri, prod_const, card_range]
  rw [ prod_range_reflect (fun i => 1 + a * Q ^ i) n]
  ring

lemma pair_eq (a b Q : K⟦X⟧) (hQ : Q = a * b) (j m : ℕ) :
    qbt Q (Q ^ (j + m + 1)) b (2 * (j + m) + 1) j
      + qbt Q (Q ^ (j + m + 1)) b (2 * (j + m) + 1) (j + m + 1 + m)
    = (Q ^ (tri (j + m)) * b ^ (j + m) * (Q ^ (j + m + 1)) ^ (j + m + 1))
      * (qbinom Q (2 * (j + m) + 1) j * a ^ (tri m) * b ^ (tri m - m) * (1 + b ^ (2 * m + 1))) := by
  unfold qbt
  rw [show 2 * (j + m) + 1 - j = j + 2 * m + 1 by omega,
    show 2 * (j + m) + 1 - (j + m + 1 + m) = j by omega]
  have hs : qbinom Q (2 * (j + m) + 1) (j + m + 1 + m) = qbinom Q (2 * (j + m) + 1) j := by
    rw [show j + m + 1 + m = 2 * (j + m) + 1 - j by omega, qbinom_symm Q (by omega)]
  rw [hs]
  obtain ⟨s, hs'⟩ : ∃ s, tri m = s + m := ⟨tri m - m, by have := le_tri m; omega⟩
  rw [show tri m - m = s by omega]
  generalize qbinom Q (2 * (j + m) + 1) j = qb
  subst hQ
  have e1 : qb * b ^ j * (a * b) ^ tri j * ((a * b) ^ (j + m + 1)) ^ (j + 2 * m + 1)
      = qb * (a ^ (tri j + (j + m + 1) * (j + 2 * m + 1))
          * b ^ (j + tri j + (j + m + 1) * (j + 2 * m + 1))) := by ring
  have e2 : qb * b ^ (j + m + 1 + m) * (a * b) ^ tri (j + m + 1 + m) * ((a * b) ^ (j + m + 1)) ^ j
      = qb * (a ^ (tri (j + m + 1 + m) + (j + m + 1) * j)
          * b ^ (j + m + 1 + m + tri (j + m + 1 + m) + (j + m + 1) * j)) := by ring
  have e3 : ((a * b) ^ tri (j + m) * b ^ (j + m) * ((a * b) ^ (j + m + 1)) ^ (j + m + 1))
      * (qb * a ^ (tri m) * b ^ s * (1 + b ^ (2 * m + 1)))
      = qb * (a ^ (tri (j + m) + (j + m + 1) * (j + m + 1) + tri m)
          * b ^ (tri (j + m) + j + m + (j + m + 1) * (j + m + 1) + s))
        + qb * (a ^ (tri (j + m) + (j + m + 1) * (j + m + 1) + tri m)
          * b ^ (tri (j + m) + j + m + (j + m + 1) * (j + m + 1) + s + (2 * m + 1))) := by ring
  have t1 := tri_add j m
  have t2 := tri_add (j + m + 1) m
  have t3 := tri_succ (j + m)
  have t4 := two_mul_tri m
  rw [e1, e2, e3,
    show tri j + (j + m + 1) * (j + 2 * m + 1) = tri (j + m) + (j + m + 1) * (j + m + 1) + tri m by
      nlinarith,
    show j + tri j + (j + m + 1) * (j + 2 * m + 1)
      = tri (j + m) + j + m + (j + m + 1) * (j + m + 1) + s by nlinarith,
    show tri (j + m + 1 + m) + (j + m + 1) * j
      = tri (j + m) + (j + m + 1) * (j + m + 1) + tri m by nlinarith,
    show j + m + 1 + m + tri (j + m + 1 + m) + (j + m + 1) * j
      = tri (j + m) + j + m + (j + m + 1) * (j + m + 1) + s + (2 * m + 1) by nlinarith]

/-- The finite triple product (with `+` signs), for `Q = a b`:
`(-b;Q)_{n+1} (-a;Q)_n = Σ_{m ≤ n} [2n+1; n-m]_Q a^{tri m} b^{tri m - m} (1 + b^{2m+1})`. -/
lemma finite_tp (a b Q : K⟦X⟧) [Nil Q] (hQ : Q = a * b) (hb : b ≠ 0) (hQ0 : Q ≠ 0) (n : ℕ) :
    (∏ i ∈ range (n + 1), (1 + b * Q ^ i)) * ∏ i ∈ range n, (1 + a * Q ^ i)
      = ∑ m ∈ range (n + 1),
          qbinom Q (2 * n + 1) (n - m) * a ^ (tri m) * b ^ (tri m - m) * (1 + b ^ (2 * m + 1)) := by
  have hC : (Q ^ (tri n) * b ^ n * (Q ^ (n + 1)) ^ (n + 1) : K⟦X⟧) ≠ 0 :=
    mul_ne_zero (mul_ne_zero (pow_ne_zero _ hQ0) (pow_ne_zero _ hb))
      (pow_ne_zero _ (pow_ne_zero _ hQ0))
  apply mul_left_cancel₀ hC
  rw [← lhs_split a b Q hQ n, qbinomial, mul_sum,
    show 2 * n + 1 + 1 = (n + 1) + (n + 1) by ring, sum_range_add,
    ← sum_range_reflect (fun t => qbt Q (Q ^ (n + 1)) b (2 * n + 1) t) (n + 1), ← sum_add_distrib]
  refine sum_congr rfl (fun m hm => ?_)
  have hm' := mem_range.mp hm
  obtain ⟨j, rfl⟩ : ∃ j, n = j + m := ⟨n - m, by omega⟩
  rw [show j + m + 1 - 1 - m = j by omega, pair_eq a b Q hQ j m, Nat.add_sub_cancel]

/-! ### Passing to the limit -/

lemma poch_add_mul_iq (Q : K⟦X⟧) [Nil Q] (L r : ℕ) :
    (X : K⟦X⟧) ^ (L + 1) ∣ poch Q (L + r) * iq Q L - 1 := by
  have hp : poch Q (L + r) = poch Q L * ∏ i ∈ range r, (1 - Q ^ (L + i + 1)) := by
    simp only [poch]; rw [prod_range_add]
  rw [hp, show poch Q L * (∏ i ∈ range r, (1 - Q ^ (L + i + 1))) * iq Q L - 1
      = (∏ i ∈ range r, (1 - Q ^ (L + i + 1))) * (poch Q L * iq Q L) - 1 by ring,
    poch_mul_iq, mul_one, ← dvd_neg, neg_sub]
  exact (X_pow_dvd_one_sub_prod _ _ _ (fun i _ => by
    rw [sub_sub_cancel]; exact X_pow_dvd_pow_of_le Q (by omega)))

lemma poch_mul_qbinom (Q : K⟦X⟧) [Nil Q] {n m : ℕ} (h : m ≤ n) :
    (X : K⟦X⟧) ^ (n - m + 1) ∣ poch Q (n + 1) * qbinom Q (2 * n + 1) (n - m) - 1 := by
  have h1 := X_pow_dvd_iq_sub_qbinom Q (show n - m ≤ 2 * n + 1 by omega)
  have h2 := poch_add_mul_iq Q (n - m) (m + 1)
  rw [show n - m + (m + 1) = n + 1 by omega] at h2
  rw [show poch Q (n + 1) * qbinom Q (2 * n + 1) (n - m) - 1
      = (poch Q (n + 1) * iq Q (n - m) - 1)
        - poch Q (n + 1) * (iq Q (n - m) - qbinom Q (2 * n + 1) (n - m)) by ring]
  exact dvd_sub h2 (((pow_dvd_pow X (by omega)).trans h1).mul_left _)

lemma limit_tp (a b Q : K⟦X⟧) [Nil Q] (hQ : Q = a * b) (ha : (X : K⟦X⟧) ∣ a) (hb : b ≠ 0)
    (hQ0 : Q ≠ 0) (n : ℕ) :
    (X : K⟦X⟧) ^ (n + 1) ∣ poch Q (n + 1) *
        ((∏ i ∈ range (n + 1), (1 + b * Q ^ i)) * ∏ i ∈ range n, (1 + a * Q ^ i))
      - ∑ m ∈ range (n + 1), a ^ (tri m) * b ^ (tri m - m) * (1 + b ^ (2 * m + 1)) := by
  rw [finite_tp a b Q hQ hb hQ0 n, mul_sum, ← sum_sub_distrib]
  refine dvd_sum (fun m hm => ?_)
  have hm' := mem_range.mp hm
  have e : ∀ p q A B C : K⟦X⟧, p * (q * A * B * C) - A * B * C = (p * q - 1) * (A * (B * C)) := by
    intros; ring
  have hX : (X : K⟦X⟧) ^ (n + 1) = X ^ (n - m + 1) * X ^ m := by
    rw [← pow_add]; congr 1; omega
  rw [e, hX]
  exact mul_dvd_mul (poch_mul_qbinom Q (by omega))
    (Dvd.dvd.mul_right ((pow_dvd_pow X (le_tri m)).trans (pow_dvd_pow_of_dvd ha _)) _)

/-! ### Truncations and the final assembly -/

lemma prod_trunc (f : ℕ → K⟦X⟧) (n L r : ℕ) (h : ∀ i, (X : K⟦X⟧) ^ n ∣ 1 - f (L + i)) :
    (X : K⟦X⟧) ^ n ∣ ∏ i ∈ range (L + r), f i - ∏ i ∈ range L, f i := by
  rw [prod_range_add]
  have := X_pow_dvd_one_sub_prod (range r) (fun i => f (L + i)) n (fun i _ => h i)
  rw [show (∏ i ∈ range L, f i) * (∏ i ∈ range r, f (L + i)) - ∏ i ∈ range L, f i
      = -(∏ i ∈ range L, f i) * (1 - ∏ i ∈ range r, f (L + i)) by ring]
  exact this.mul_left _

lemma dvd_mul_sub_mul' {n : ℕ} {f f' g g' : K⟦X⟧} (hf : (X : K⟦X⟧) ^ n ∣ f - f')
    (hg : (X : K⟦X⟧) ^ n ∣ g - g') : (X : K⟦X⟧) ^ n ∣ f * g - f' * g' := by
  rw [show f * g - f' * g' = f * (g - g') + (f - f') * g' by ring]
  exact dvd_add (hg.mul_left _) (hf.mul_right _)

lemma qPochInf_trunc (a d N : ℕ) (ha : 1 ≤ a) (hd : 1 ≤ d) :
    (X : K⟦X⟧) ^ (N + 1) ∣ (qPochInf a d : K⟦X⟧) - ∏ i ∈ range (N + 1), (1 - X ^ (a + d * i)) := by
  rw [X_pow_dvd_iff]
  intro m hm
  rw [map_sub, qPochInf, coeff_mk, sub_eq_zero]
  have := prod_trunc (fun i => (1 - X ^ (a + d * i) : K⟦X⟧)) (m + 1) (m + 1) (N - m)
    (fun i => by simp only [sub_sub_cancel]; exact pow_dvd_pow X (by nlinarith))
  rw [show m + 1 + (N - m) = N + 1 by omega] at this
  have hc := X_pow_dvd_iff.mp this m (Nat.lt_succ_self m)
  rw [map_sub, sub_eq_zero] at hc
  exact hc.symm

lemma term_eq (k m : ℕ) :
    ((-X ^ (k + 1)) ^ (tri m) * (-X) ^ (tri m - m) * (1 + (-X) ^ (2 * m + 1)) : K⟦X⟧)
      = (-1) ^ m * X ^ (((k + 2) * m ^ 2 + k * m) / 2) * (1 - X ^ (2 * m + 1)) := by
  obtain ⟨s, hs'⟩ : ∃ s, tri m = s + m := ⟨tri m - m, by have := le_tri m; omega⟩
  have t := two_mul_tri m
  rw [show tri m - m = s by omega, hs',
    Nat.div_eq_of_eq_mul_left two_pos (show (k + 2) * m ^ 2 + k * m = ((k + 1) * (s + m) + s) * 2
      by nlinarith),
    Odd.neg_pow ⟨m, rfl⟩, neg_pow (X ^ (k + 1) : K⟦X⟧), neg_pow (X : K⟦X⟧) s]
  have hsg : ((-1 : K⟦X⟧) ^ s) ^ 2 = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
  linear_combination ((-1 : K⟦X⟧) ^ m * X ^ ((k + 1) * (s + m) + s) * (1 - X ^ (2 * m + 1))) * hsg

/-- Jacobi triple product, the case needed: (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ = Σ_{n∈ℤ} (-1)^n q^{((k+2)n²+kn)/2}. -/
theorem jtp_theta (k : ℕ) :
    (qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2) : K⟦X⟧) = thetaSum k := by
  have : Nil (X ^ (k + 2) : K⟦X⟧) := ⟨by simp⟩
  ext N
  have hlim := limit_tp (-X ^ (k + 1)) (-X) (X ^ (k + 2) : K⟦X⟧) (by ring)
    (dvd_neg.mpr (dvd_pow_self _ (Nat.succ_ne_zero k))) (neg_ne_zero.mpr X_ne_zero)
    (pow_ne_zero _ X_ne_zero) N
  have h1 := qPochInf_trunc (K := K) 1 (k + 2) N le_rfl (by omega)
  have h2 := qPochInf_trunc (K := K) (k + 1) (k + 2) N (by omega) (by omega)
  have h3 := qPochInf_trunc (K := K) (k + 2) (k + 2) N (by omega) (by omega)
  have e1 : ∏ i ∈ range (N + 1), (1 - X ^ (1 + (k + 2) * i) : K⟦X⟧)
      = ∏ i ∈ range (N + 1), (1 + (-X) * (X ^ (k + 2)) ^ i) :=
    prod_congr rfl fun i _ => by ring
  have e3 : ∏ i ∈ range (N + 1), (1 - X ^ (k + 2 + (k + 2) * i) : K⟦X⟧) = poch (X ^ (k + 2)) (N + 1) :=
    prod_congr rfl fun i _ => by ring
  have h2' : (X : K⟦X⟧) ^ (N + 1) ∣ ∏ i ∈ range (N + 1), (1 - X ^ (k + 1 + (k + 2) * i) : K⟦X⟧)
      - ∏ i ∈ range N, (1 + (-X ^ (k + 1)) * (X ^ (k + 2)) ^ i) := by
    have := prod_trunc (fun i => (1 - X ^ (k + 1 + (k + 2) * i) : K⟦X⟧)) (N + 1) N 1
      (fun i => by simp only [sub_sub_cancel]; exact pow_dvd_pow X (by nlinarith))
    rw [show ∏ i ∈ range N, (1 + (-X ^ (k + 1)) * (X ^ (k + 2)) ^ i : K⟦X⟧)
        = ∏ i ∈ range N, (1 - X ^ (k + 1 + (k + 2) * i)) from prod_congr rfl fun i _ => by ring]
    exact this
  rw [e1] at h1
  rw [e3] at h3
  have h2'' := dvd_add h2 h2'
  rw [sub_add_sub_cancel] at h2''
  have h12 := dvd_mul_sub_mul' (dvd_mul_sub_mul' h1 h2'') h3
  have hfin := dvd_add h12 hlim
  rw [show ∀ A B C P1 P2 P3 S : K⟦X⟧, A * B * C - P1 * P2 * P3 + (P3 * (P1 * P2) - S)
      = A * B * C - S from fun _ _ _ _ _ _ _ => by ring] at hfin
  have hc := X_pow_dvd_iff.mp hfin N (Nat.lt_succ_self N)
  rw [map_sub, sub_eq_zero] at hc
  rw [hc, thetaSum, coeff_mk, map_sum]
  exact sum_congr rfl (fun m _ => by rw [term_eq])

end KRTW
