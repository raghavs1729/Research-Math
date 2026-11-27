/-
Nahm/E7Lemmas.lean — Lemma 4 and Proposition 5 of `proofs/E7.md`, over an arbitrary base `Q`.

* Lemma 4: `Σ_m Q^{m²+Rm} / ((Q)_m (Q)_{m+R+a}) = (1/(Q;Q)_∞) Σ_{j ≤ a} [a;j] Q^{j²+Rj}`
  (from Cauchy's identity and the Durfee identity).
* Proposition 5: `T(a) = Σ_{n ∈ ℕ³} Q^{N_1²+N_2²+N_3² + a(N_2+N_3)}/∏(Q)_{n_i}
     = ((Q)_a/(Q;Q)_∞) Σ_r Q^{r²} α_r Σ_{j ≤ a} [a;j] Q^{j²+2rj}`, with `α_r` the third iterate of the
  unit Bailey pair relative to `Q^a` (Lemma 3 of E7.md is `IsBP.step` applied to `isBP_iq`).
-/
import Nahm.AGTheta

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K] (Q : K⟦X⟧) [Nil Q]

omit [Nil Q] in
/-- `(Q^{a+k+1};Q)_m = (Q)_{a+k+m} / (Q)_{a+k}`. -/
lemma gq_pow (a k m : ℕ) : gq Q (Q ^ a * Q ^ (k + 1)) m * poch Q (a + k) = poch Q (a + k + m) := by
  rw [poch_add Q (a + k) m, gq, mul_comm]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  rw [← pow_add, ← pow_add]
  congr 2
  omega

/-- A Bailey pair relative to `Q^a`: `β_n = (Q)_a Σ_{r ≤ n} α_r / ((Q)_{n-r} (Q)_{a+n+r})`. -/
lemma isBP_pow_eq (a : ℕ) {α β : ℕ → K⟦X⟧} (h : IsBP Q (Q ^ a) α β) (n : ℕ) :
    β n = poch Q a * ∑ r ∈ range (n + 1), α r * iq Q (n - r) * iq Q (n + r + a) := by
  have hb := h n
  rw [bps] at hb
  apply mul_left_cancel₀ (poch_ne_zero' Q (a + 2 * n))
  have e1 : poch Q (a + 2 * n) = poch Q a * gq Q (Q ^ a * Q) (2 * n) := by
    have := gq_pow Q a 0 (2 * n)
    rw [pow_one, add_zero] at this
    rw [← this]; ring
  rw [e1, mul_assoc, hb, mul_sum, mul_sum, mul_sum]
  refine sum_congr rfl (fun r hr => ?_)
  have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
  have e2 : gq Q (Q ^ a * Q ^ (n + r + 1)) (n - r) = poch Q (a + 2 * n) * iq Q (n + r + a) := by
    have := gq_pow Q a (n + r) (n - r)
    rw [show a + (n + r) + (n - r) = a + 2 * n by omega] at this
    rw [← this, show n + r + a = a + (n + r) by ring, mul_assoc, poch_mul_iq, mul_one]
  rw [e2, e1]
  ring

variable [TopologicalSpace K] [DiscreteTopology K]

/-- **Lemma 4 of `E7.md`.** -/
theorem lemma4 (R a : ℕ) :
    ∑' m, Q ^ (m ^ 2 + R * m) * (iq Q m * iq Q (m + R + a))
      = (pinf Q)⁻¹ * ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + R * j) := by
  set F : ℕ → ℕ → K⟦X⟧ := fun j m => qbinom Q a j * Q ^ (j ^ 2 + R * j + (m - j) * (m + R + a)) *
    (iq Q (m - j) * iq Q (m + R + a)) with hF
  have hdvd : ∀ j m, (X : K⟦X⟧) ^ m ∣ F j m := by
    intro j m
    refine Dvd.dvd.mul_right (Dvd.dvd.mul_left (X_pow_dvd_pow_of_le Q ?_) _) _
    rcases Nat.lt_or_ge m j with h | h
    · nlinarith
    · rcases Nat.eq_zero_or_pos m with h0 | h0
      · omega
      · have : m - j ≤ (m - j) * (m + R + a) := Nat.le_mul_of_pos_right _ (by omega)
        have h2 : j ≤ j ^ 2 := by nlinarith
        omega
  have hsh := tsum_tsum_shift (fun j => j) (fun j => le_refl j) F hdvd
  -- right side as a double sum
  have hR : (pinf Q)⁻¹ * ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + R * j)
      = ∑' j, ∑' m', F j (m' + j) := by
    rw [mul_sum, tsum_eq_sum (s := range (a + 1))]
    · refine sum_congr rfl (fun j _ => ?_)
      rw [← cauchy Q (R + a + j), ← (summable_cauchy Q (R + a + j)).tsum_mul_right]
      refine tsum_congr (fun m' => ?_)
      simp only [hF, Nat.add_sub_cancel]
      rw [show m' + j + R + a = m' + (R + a + j) by ring, pow_add]
      ring
    · intro j hj
      have : a < j := by simpa using hj
      simp [hF, qbinom_of_lt Q this]
  rw [hR, hsh]
  refine tsum_congr (fun m => ?_)
  -- the inner sum by the Durfee identity
  have hd := durfee Q a m
  have hsub : ∑ j ∈ range (m + 1), (if j ≤ m then F j m else 0)
      = ∑ j ∈ range (min a m + 1), F j m := by
    symm
    rw [← sum_subset (range_subset_range.mpr (Nat.succ_le_succ (min_le_right a m)))]
    · refine sum_congr rfl (fun j hj => ?_)
      rw [if_pos (by have := mem_range.mp hj; omega)]
    · intro j hj hj'
      have h1 : j ≤ m := by have := mem_range.mp hj; omega
      have h2 : a < j := by
        simp only [mem_range, not_lt] at hj'
        omega
      rw [if_pos h1]
      simp [hF, qbinom_of_lt Q h2]
  rw [hsub]
  have e : ∀ j ∈ range (min a m + 1), F j m
      = Q ^ (m ^ 2 + R * m) * iq Q (m + R + a) * poch Q a *
          (Q ^ ((a - j) * (m - j)) * (iq Q (a - j) * iq Q (m - j) * iq Q j)) := by
    intro j hj
    have hj' : j ≤ min a m := Nat.lt_succ_iff.mp (mem_range.mp hj)
    have hja : j ≤ a := le_trans hj' (min_le_left a m)
    have hjm : j ≤ m := le_trans hj' (min_le_right a m)
    have hexp : j ^ 2 + R * j + (m - j) * (m + R + a) = m ^ 2 + R * m + (a - j) * (m - j) := by
      obtain ⟨s, rfl⟩ := Nat.exists_eq_add_of_le hjm
      obtain ⟨t, rfl⟩ := Nat.exists_eq_add_of_le hja
      rw [show j + s - j = s by omega, show j + t - j = t by omega]
      ring
    simp only [hF, qbinom, hja, if_true]
    rw [hexp, pow_add]
    ring
  rw [sum_congr rfl e, ← mul_sum, ← hd]
  have := poch_mul_iq Q a
  linear_combination (-(Q ^ (m ^ 2 + R * m) * iq Q (m + R + a) * iq Q m)) * this

/-- The outermost sum for a Bailey pair relative to `Q^a`. -/
theorem top_pow (a : ℕ) (α : ℕ → K⟦X⟧) :
    ∑' n, Q ^ (n ^ 2) * ∑ r ∈ range (n + 1), α r * iq Q (n - r) * iq Q (n + r + a)
      = (pinf Q)⁻¹ * ∑' r, Q ^ (r ^ 2) * α r *
          ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + 2 * r * j) := by
  set F : ℕ → ℕ → K⟦X⟧ := fun r M => Q ^ (M ^ 2) * (α r * iq Q (M - r) * iq Q (M + r + a)) with hF
  have hdvd : ∀ r M, (X : K⟦X⟧) ^ M ∣ F r M := fun r M =>
    Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _
  have hsh := tsum_tsum_shift (fun r => r) (fun r => le_refl r) F hdvd
  have e1 : ∑' n, Q ^ (n ^ 2) * ∑ r ∈ range (n + 1), α r * iq Q (n - r) * iq Q (n + r + a)
      = ∑' M, ∑ r ∈ range (M + 1), if r ≤ M then F r M else 0 := by
    refine tsum_congr (fun M => ?_)
    rw [mul_sum]
    refine sum_congr rfl (fun r hr => ?_)
    rw [if_pos (Nat.lt_succ_iff.mp (mem_range.mp hr))]
  rw [e1, ← hsh]
  have e2 : ∀ r, ∑' m, F r (m + r) = (pinf Q)⁻¹ * (Q ^ (r ^ 2) * α r *
      ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + 2 * r * j)) := by
    intro r
    have hs : Summable fun m => Q ^ (m ^ 2 + 2 * r * m) * (iq Q m * iq Q (m + 2 * r + a)) :=
      summable_nat (fun m => Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _)
    rw [show (pinf Q)⁻¹ * (Q ^ (r ^ 2) * α r * ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + 2 * r * j))
        = Q ^ (r ^ 2) * α r * ((pinf Q)⁻¹ * ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + 2 * r * j))
        by ring, ← lemma4, ← hs.tsum_mul_left]
    refine tsum_congr (fun m => ?_)
    simp only [hF, Nat.add_sub_cancel, show m + r + r + a = m + 2 * r + a by ring]
    rw [show (m + r) ^ 2 = r ^ 2 + (m ^ 2 + 2 * r * m) by ring, pow_add]
    ring
  rw [tsum_congr e2]
  have hs : Summable fun r => Q ^ (r ^ 2) * α r *
      ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + 2 * r * j) :=
    summable_nat (fun r => Dvd.dvd.mul_right (Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _) _)
  rw [hs.tsum_mul_left]

/-- `α_r(a)`: the unit Bailey pair relative to `Q^a`, iterated three times. -/
noncomputable def alE7 (a r : ℕ) : K⟦X⟧ := ((Q ^ a) ^ r * Q ^ (r ^ 2)) ^ 3 * ua Q (Q ^ a) r

/-- **Proposition 5 of `E7.md`** (`T(a)`, before inserting the explicit `α_r`). -/
theorem prop5 (a : ℕ) :
    ∑' m : Fin 3 → ℕ, mterm Q ![0, a, a] m
      = poch Q a * ((pinf Q)⁻¹ * ∑' r, Q ^ (r ^ 2) * alE7 Q a r *
          ∑ j ∈ range (a + 1), qbinom Q a j * Q ^ (j ^ 2 + 2 * r * j)) := by
  have hb : (![0, a, a] : Fin 3 → ℕ) = Fin.cons 0 ![a, a] := rfl
  rw [hb, tsum_mterm_eq]
  have hc : IsBP Q (Q ^ a) (alE7 Q a) (chainG Q 2 ![a, a]) := by
    have h := ((isBP_iq Q (Q ^ a)).step).step
    have e : alE7 Q a = fun r => (Q ^ a) ^ r * Q ^ (r ^ 2) * ((Q ^ a) ^ r * Q ^ (r ^ 2) *
        ((Q ^ a) ^ r * Q ^ (r ^ 2) * ua Q (Q ^ a) r)) := by
      funext r; simp only [alE7]; ring
    rw [e]
    exact h
  rw [tsum_congr (fun n => by rw [pow_zero, one_pow, one_mul, isBP_pow_eq Q a hc n])]
  have hs : Summable fun n => Q ^ (n ^ 2) *
      ∑ r ∈ range (n + 1), alE7 Q a r * iq Q (n - r) * iq Q (n + r + a) :=
    summable_nat (fun n => Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _)
  rw [show (fun n => Q ^ (n ^ 2) * (poch Q a * ∑ r ∈ range (n + 1), alE7 Q a r * iq Q (n - r) *
      iq Q (n + r + a))) = fun n => poch Q a * (Q ^ (n ^ 2) *
        ∑ r ∈ range (n + 1), alE7 Q a r * iq Q (n - r) * iq Q (n + r + a)) from
    funext fun n => by ring, hs.tsum_mul_left, top_pow]

end Nahm
