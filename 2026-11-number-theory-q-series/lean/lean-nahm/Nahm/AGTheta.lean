/-
Nahm/AGTheta.lean — the Andrews–Gordon identities in theta form, for every modulus `2k+1` and every
`1 ≤ i ≤ k`, over an arbitrary base `Q`:

  `Σ_{m ∈ ℕ^{k-1}} Q^{M_1²+⋯+M_{k-1}² + M_i+⋯+M_{k-1}} / ((Q)_{m_1}⋯(Q)_{m_{k-1}})
      = (1/(Q;Q)_∞) Σ_{n ≥ 0} (-1)^n Q^{(2k+1-i) n(n+1)/2 + i n(n-1)/2} (1 - Q^{i(2n+1)})`

(the right side is `(1/(Q;Q)_∞) Σ_{n ∈ ℤ} (-1)^n Q^{((2k+1)n² + (2k+1-2i)n)/2}` with the terms `n` and
`-n-1` paired).  Proof (Andrews; Agarwal–Andrews–Bressoud): the unit Bailey pair relative to `Q`,
`k - i` steps of the Bailey chain relative to `Q`, the change of `a` (`Q → 1`), `i - 2` steps relative
to `1`, and Cauchy's identity for the outermost sum (for `i = 1`: `k - 1` steps relative to `Q`).
-/
import Nahm.Multisum

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K] (Q : K⟦X⟧) [Nil Q]

lemma ua_Q_mul (r : ℕ) : ua Q Q r * (1 - Q) = (-1) ^ r * Q ^ tri (r - 1) * (1 - Q ^ (2 * r + 1)) := by
  cases r with
  | zero => simp [ua]
  | succ r =>
    have h := gq_QQ Q r
    have h2 := iq_mul_poch Q (r + 1)
    simp only [ua, Nat.add_sub_cancel]
    have e : Q * Q ^ (2 * r + 2) = Q ^ (2 * (r + 1) + 1) := by ring
    rw [← e]
    linear_combination ((-1) ^ (r + 1) * Q ^ tri r * (1 - Q * Q ^ (2 * r + 2)) * iq Q (r + 1)) * h
      + ((-1) ^ (r + 1) * Q ^ tri r * (1 - Q * Q ^ (2 * r + 2))) * h2

/-- The chain with all levels relative to `Q`. -/
theorem chainG_isBP_Q : ∀ (s : ℕ) (c : Fin s → ℕ), (∀ l, c l = 1) →
    IsBP Q Q (fun r => (Q ^ r * Q ^ (r ^ 2)) ^ (s + 1) * ua Q Q r) (chainG Q s c)
  | 0, c, _ => by
    have h := isBP_iq Q Q
    simp only [zero_add, pow_one]
    exact h
  | s + 1, c, hc => by
    have ih := chainG_isBP_Q s (Fin.tail c) (fun l => hc l.succ)
    have h := ih.step
    have e1 : chainG Q (s + 1) c = bstep Q Q (chainG Q s (Fin.tail c)) := by
      simp only [chainG, hc 0, pow_one]
    have e2 : (fun r => (Q ^ r * Q ^ (r ^ 2)) ^ (s + 1 + 1) * ua Q Q r)
        = fun r => Q ^ r * Q ^ (r ^ 2) * ((Q ^ r * Q ^ (r ^ 2)) ^ (s + 1) * ua Q Q r) := by
      funext r; ring
    rw [e1, e2]
    exact h

/-- `γ_r = (-1)^r Q^{(q+1)(r²+r) + r(r-1)/2}`. -/
noncomputable def gam (q r : ℕ) : K⟦X⟧ := (-1) ^ r * Q ^ ((q + 1) * (r ^ 2 + r) + tri (r - 1))

/-- `α'_r = γ_r - Q^{2r-1} γ_{r-1}`. -/
noncomputable def alph1 (q r : ℕ) : K⟦X⟧ :=
  gam Q q r - if r = 0 then 0 else Q ^ (2 * r - 1) * gam Q q (r - 1)

theorem chainG_isBP_one (q : ℕ) : ∀ (p : ℕ) (c : Fin (q + p) → ℕ),
    (∀ l : Fin (q + p), c l = if (l : ℕ) < p then 0 else 1) →
    IsBP Q 1 (fun r => (Q ^ (r ^ 2)) ^ p * alph1 Q q r) (chainG Q (q + p) c)
  | 0, c, hc => by
    have h := (chainG_isBP_Q Q q c (fun l => by simp [hc l])).toOne Q (gam Q q) (fun r => by
      rw [mul_assoc, ua_Q_mul, gam]
      rw [show (Q ^ r * Q ^ (r ^ 2)) ^ (q + 1) = Q ^ ((q + 1) * (r ^ 2 + r)) by ring]
      ring)
    simp only [pow_zero, one_mul]
    exact h
  | p + 1, c, hc => by
    have ih := chainG_isBP_one q p (Fin.tail c) (fun l => by
      rw [Fin.tail, hc l.succ]
      simp only [Fin.val_succ]
      by_cases h : (l : ℕ) < p
      · rw [if_pos (by omega), if_pos h]
      · rw [if_neg (by omega), if_neg h])
    have h := ih.step
    have hc0 : c 0 = 0 := by rw [hc 0]; simp
    have e1 : chainG Q (q + (p + 1)) c = bstep Q 1 (chainG Q (q + p) (Fin.tail c)) := by
      change bstep Q (Q ^ c 0) (chainG Q (q + p) (Fin.tail c)) = _
      rw [hc0, pow_zero]
    have e2 : (fun r => (Q ^ (r ^ 2)) ^ (p + 1) * alph1 Q q r)
        = fun r => (1 : K⟦X⟧) ^ r * Q ^ (r ^ 2) * ((Q ^ (r ^ 2)) ^ p * alph1 Q q r) := by
      funext r; ring
    rw [e1, e2]
    exact h

lemma poch_ne_zero' (n : ℕ) : poch Q n ≠ 0 := by
  intro h
  have := constantCoeff_poch Q n
  rw [h] at this
  simp at this

/-- A Bailey pair relative to `1`: `β_n = Σ_{r ≤ n} α_r / ((Q)_{n-r} (Q)_{n+r})`. -/
lemma isBP_one_eq {α β : ℕ → K⟦X⟧} (h : IsBP Q 1 α β) (n : ℕ) :
    β n = ∑ r ∈ range (n + 1), α r * iq Q (n - r) * iq Q (n + r) := by
  have hb := h n
  rw [gq_one_eq_poch, bps] at hb
  apply mul_left_cancel₀ (poch_ne_zero' Q (2 * n))
  rw [hb, mul_sum]
  refine sum_congr rfl (fun r hr => ?_)
  rw [gq_top Q n r (Nat.lt_succ_iff.mp (mem_range.mp hr))]
  ring

/-- A Bailey pair relative to `Q`: `β_n = (1 - Q) Σ_{r ≤ n} α_r / ((Q)_{n-r} (Q)_{n+r+1})`. -/
lemma isBP_Q_eq {α β : ℕ → K⟦X⟧} (h : IsBP Q Q α β) (n : ℕ) :
    β n = ∑ r ∈ range (n + 1), (1 - Q) * α r * iq Q (n - r) * iq Q (n + r + 1) := by
  have hb := h n
  rw [bps] at hb
  apply mul_left_cancel₀ (poch_ne_zero' Q (2 * n + 1))
  rw [← gq_QQ, mul_assoc, hb, mul_sum, mul_sum]
  refine sum_congr rfl (fun r hr => ?_)
  rw [gq_top' Q n r (Nat.lt_succ_iff.mp (mem_range.mp hr)), ← gq_QQ]
  ring

variable [TopologicalSpace K] [DiscreteTopology K]

/-- The outermost sum, relative to `1`, by Cauchy's identity. -/
theorem top_one (α : ℕ → K⟦X⟧) :
    ∑' n, Q ^ (n ^ 2) * ∑ r ∈ range (n + 1), α r * iq Q (n - r) * iq Q (n + r)
      = (pinf Q)⁻¹ * ∑' r, Q ^ (r ^ 2) * α r := by
  set F : ℕ → ℕ → K⟦X⟧ := fun r M => Q ^ (M ^ 2) * (α r * iq Q (M - r) * iq Q (M + r)) with hF
  have hdvd : ∀ r M, (X : K⟦X⟧) ^ M ∣ F r M := fun r M =>
    Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _
  have hsh := tsum_tsum_shift (fun r => r) (fun r => le_refl r) F hdvd
  have e1 : ∑' n, Q ^ (n ^ 2) * ∑ r ∈ range (n + 1), α r * iq Q (n - r) * iq Q (n + r)
      = ∑' M, ∑ r ∈ range (M + 1), if r ≤ M then F r M else 0 := by
    refine tsum_congr (fun M => ?_)
    rw [mul_sum]
    refine sum_congr rfl (fun r hr => ?_)
    rw [if_pos (Nat.lt_succ_iff.mp (mem_range.mp hr))]
  rw [e1, ← hsh]
  have e2 : ∀ r, ∑' m, F r (m + r) = Q ^ (r ^ 2) * α r * (pinf Q)⁻¹ := by
    intro r
    have hs := summable_cauchy Q (2 * r)
    rw [← cauchy Q (2 * r), ← hs.tsum_mul_left]
    refine tsum_congr (fun m => ?_)
    simp only [hF, Nat.add_sub_cancel, show m + r + r = m + 2 * r by ring]
    rw [show (m + r) ^ 2 = r ^ 2 + m * (m + 2 * r) by ring, pow_add]
    ring
  rw [tsum_congr e2]
  have hs : Summable fun r => Q ^ (r ^ 2) * α r :=
    summable_nat (fun r => Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _)
  rw [hs.tsum_mul_right, mul_comm]

/-- The outermost sum, relative to `Q`, by Cauchy's identity. -/
theorem top_Q (α : ℕ → K⟦X⟧) :
    ∑' n, Q ^ n * Q ^ (n ^ 2) * ∑ r ∈ range (n + 1), (1 - Q) * α r * iq Q (n - r) * iq Q (n + r + 1)
      = (pinf Q)⁻¹ * ∑' r, (1 - Q) * Q ^ (r ^ 2 + r) * α r := by
  set F : ℕ → ℕ → K⟦X⟧ := fun r M => Q ^ M * Q ^ (M ^ 2) *
    ((1 - Q) * α r * iq Q (M - r) * iq Q (M + r + 1)) with hF
  have hdvd : ∀ r M, (X : K⟦X⟧) ^ M ∣ F r M := fun r M =>
    Dvd.dvd.mul_right (Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q le_rfl) _) _
  have hsh := tsum_tsum_shift (fun r => r) (fun r => le_refl r) F hdvd
  have e1 : ∑' n, Q ^ n * Q ^ (n ^ 2) * ∑ r ∈ range (n + 1), (1 - Q) * α r * iq Q (n - r) * iq Q (n + r + 1)
      = ∑' M, ∑ r ∈ range (M + 1), if r ≤ M then F r M else 0 := by
    refine tsum_congr (fun M => ?_)
    rw [mul_sum]
    refine sum_congr rfl (fun r hr => ?_)
    rw [if_pos (Nat.lt_succ_iff.mp (mem_range.mp hr))]
  rw [e1, ← hsh]
  have e2 : ∀ r, ∑' m, F r (m + r) = (1 - Q) * Q ^ (r ^ 2 + r) * α r * (pinf Q)⁻¹ := by
    intro r
    have hs := summable_cauchy Q (2 * r + 1)
    rw [← cauchy Q (2 * r + 1), ← hs.tsum_mul_left]
    refine tsum_congr (fun m => ?_)
    simp only [hF, Nat.add_sub_cancel, show m + r + r + 1 = m + (2 * r + 1) by ring]
    rw [← pow_add, show m + r + (m + r) ^ 2 = (r ^ 2 + r) + m * (m + (2 * r + 1)) by ring, pow_add]
    ring
  rw [tsum_congr e2]
  have hs : Summable fun r => (1 - Q) * Q ^ (r ^ 2 + r) * α r :=
    summable_nat (fun r => Dvd.dvd.mul_right (Dvd.dvd.mul_left (X_pow_dvd_pow_of_le Q (by nlinarith)) _) _)
  rw [hs.tsum_mul_right, mul_comm]

/-- The theta series `Σ_{n ≥ 0} (-1)^n Q^{c n(n+1)/2 + i n(n-1)/2} (1 - Q^{i(2n+1)})`
(`= Σ_{n ∈ ℤ} (-1)^n Q^{((c+i)n² + (c-i)n)/2}`, the terms `n` and `-n-1` paired). -/
noncomputable def thetaP (c i : ℕ) : K⟦X⟧ :=
  ∑' n, (-1) ^ n * Q ^ (c * tri n + i * tri (n - 1)) * (1 - Q ^ (i * (2 * n + 1)))

omit [Nil Q] [TopologicalSpace K] [DiscreteTopology K] in
lemma two_tri_pred (t : ℕ) : 2 * tri (t - 1) + t = t ^ 2 := by
  cases t with
  | zero => simp
  | succ t => have := two_mul_tri t; simp only [Nat.add_sub_cancel]; nlinarith

/-- The Andrews–Gordon multiple sum, linear coefficients: `b_l = 1` iff `l + 1 ≥ i` (0-based `l`). -/
def agb (i : ℕ) {r : ℕ} : Fin r → ℕ := fun l => if i ≤ (l : ℕ) + 1 then 1 else 0

/-- **Andrews–Gordon, theta form, `i = 1`**. -/
theorem ag_theta_one (s : ℕ) :
    ∑' m, mterm Q (agb 1 : Fin (s + 1) → ℕ) m = (pinf Q)⁻¹ * thetaP Q (2 * (s + 2)) 1 := by
  have hb : (agb 1 : Fin (s + 1) → ℕ) = Fin.cons 1 (fun _ => 1) := by
    funext l
    refine Fin.cases ?_ (fun l => ?_) l <;> simp [agb]
  rw [hb, tsum_mterm_eq, pow_one]
  have hc := chainG_isBP_Q Q s (fun _ => 1) (fun _ => rfl)
  rw [tsum_congr (fun n => by rw [isBP_Q_eq Q hc n]), top_Q]
  congr 1
  refine tsum_congr (fun r => ?_)
  rw [show (1 - Q) * Q ^ (r ^ 2 + r) * ((Q ^ r * Q ^ (r ^ 2)) ^ (s + 1) * ua Q Q r)
      = Q ^ ((s + 2) * (r ^ 2 + r)) * (ua Q Q r * (1 - Q)) by ring, ua_Q_mul]
  have e : (s + 2) * (r ^ 2 + r) + tri (r - 1) = 2 * (s + 2) * tri r + 1 * tri (r - 1) := by
    have := two_mul_tri r
    nlinarith
  rw [← e, pow_add]
  ring

/-- **Andrews–Gordon, theta form, `i = p + 2 ≥ 2`** (modulus `2k+1`, `k = q + p + 2`). -/
theorem ag_theta_two (p q : ℕ) :
    ∑' m, mterm Q (agb (p + 2) : Fin (q + p + 1) → ℕ) m
      = (pinf Q)⁻¹ * thetaP Q (p + 2 * q + 3) (p + 2) := by
  have hb : (agb (p + 2) : Fin (q + p + 1) → ℕ)
      = Fin.cons 0 (fun l : Fin (q + p) => if (l : ℕ) < p then 0 else 1) := by
    funext l
    refine Fin.cases ?_ (fun l => ?_) l
    · simp [agb]
    · simp only [agb, Fin.cons_succ, Fin.val_succ]
      by_cases h : (l : ℕ) < p
      · rw [if_neg (by omega), if_pos h]
      · rw [if_pos (by omega), if_neg h]
  rw [hb, tsum_mterm_eq]
  have hc := chainG_isBP_one Q q p _ (fun _ => rfl)
  rw [tsum_congr (fun n => by rw [pow_zero, one_pow, one_mul, isBP_one_eq Q hc n]), top_one]
  congr 1
  -- the theta series
  set f : ℕ → K⟦X⟧ := fun r => Q ^ (r ^ 2) * ((Q ^ (r ^ 2)) ^ p * gam Q q r) with hf
  set g : ℕ → K⟦X⟧ := fun r => Q ^ (r ^ 2) * ((Q ^ (r ^ 2)) ^ p *
    (if r = 0 then 0 else Q ^ (2 * r - 1) * gam Q q (r - 1))) with hg
  have hsf : Summable f := summable_nat (fun r => Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _)
  have hsg : Summable g := summable_nat (fun r => Dvd.dvd.mul_right (X_pow_dvd_pow_of_le Q (by nlinarith)) _)
  have e1 : ∑' r, Q ^ (r ^ 2) * ((Q ^ (r ^ 2)) ^ p * alph1 Q q r) = ∑' r, f r - ∑' r, g r := by
    rw [← hsf.tsum_sub hsg]
    refine tsum_congr (fun r => ?_)
    simp only [hf, hg, alph1]
    ring
  have e2 : ∑' r, g r = ∑' t, g (t + 1) := by
    rw [hsg.tsum_eq_zero_add]
    simp [hg]
  have hsg' : Summable fun t => g (t + 1) := (summable_nat_add_iff 1).mpr hsg
  rw [e1, e2, ← hsf.tsum_sub hsg', thetaP]
  refine tsum_congr (fun t => ?_)
  simp only [hf, hg, if_neg (Nat.succ_ne_zero t), Nat.add_sub_cancel, gam]
  have ht := two_tri_pred t
  have ea : t ^ 2 + t ^ 2 * p + ((q + 1) * (t ^ 2 + t) + tri (t - 1))
      = (p + 2 * q + 3) * tri t + (p + 2) * tri (t - 1) := by
    have := tri_pred t
    zify at ht this ⊢
    linear_combination (-(↑p + ↑q + 2 : ℤ)) * ht + (-(↑p + 2 * ↑q + 3 : ℤ)) * this
  have eb : (t + 1) ^ 2 + (t + 1) ^ 2 * p + (2 * (t + 1) - 1 + ((q + 1) * (t ^ 2 + t) + tri (t - 1)))
      = (p + 2 * q + 3) * tri t + (p + 2) * tri (t - 1) + (p + 2) * (2 * t + 1) := by
    rw [show 2 * (t + 1) - 1 = 2 * t + 1 by omega, ← ea]
    ring
  calc Q ^ t ^ 2 * ((Q ^ t ^ 2) ^ p * ((-1) ^ t * Q ^ ((q + 1) * (t ^ 2 + t) + tri (t - 1)))) -
      Q ^ (t + 1) ^ 2 * ((Q ^ (t + 1) ^ 2) ^ p * (Q ^ (2 * (t + 1) - 1) *
        ((-1) ^ t * Q ^ ((q + 1) * (t ^ 2 + t) + tri (t - 1)))))
      = (-1) ^ t * Q ^ (t ^ 2 + t ^ 2 * p + ((q + 1) * (t ^ 2 + t) + tri (t - 1)))
        - (-1) ^ t * Q ^ ((t + 1) ^ 2 + (t + 1) ^ 2 * p + (2 * (t + 1) - 1 +
          ((q + 1) * (t ^ 2 + t) + tri (t - 1)))) := by
        rw [← pow_mul, ← pow_mul]; simp only [pow_add]; ring
    _ = _ := by rw [ea, eb, pow_add]; ring

end Nahm
