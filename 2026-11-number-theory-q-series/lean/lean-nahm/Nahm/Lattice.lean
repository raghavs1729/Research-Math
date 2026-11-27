/-
Nahm/Lattice.lean — the change of `a` of the Bailey lattice (Agarwal–Andrews–Bressoud), in the case
`a = q → a = 1`: a Bailey pair relative to `Q` gives a Bailey pair relative to `1` with the same `β`.

If `(α, β)` is a Bailey pair relative to `Q` and `α_r (1 - Q) = γ_r (1 - Q^{2r+1})`, then `(α', β)` with
`α'_r = γ_r - Q^{2r-1} γ_{r-1}` (`α'_0 = γ_0`) is a Bailey pair relative to `1`.
-/
import Nahm.Bailey
import KRTW.Stage2

namespace Nahm

open PowerSeries Finset KRTW

variable {K : Type*} [Field K] (Q : K⟦X⟧) [Nil Q]

omit [Nil Q] in
lemma gq_one_eq_poch (m : ℕ) : gq Q (1 * Q) m = poch Q m := by
  simp only [gq, poch, one_mul]
  refine prod_congr rfl (fun i _ => ?_)
  rw [pow_succ']

omit [Nil Q] in
lemma gq_QQ (m : ℕ) : (1 - Q) * gq Q (Q * Q) m = poch Q (m + 1) := by
  simp only [gq, poch]
  rw [prod_range_succ']
  simp only [zero_add, pow_one]
  rw [mul_comm]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  ring

/-- `(Q^{n+r+1};Q)_{n-r} = (Q)_{2n} / (Q)_{n+r}`. -/
lemma gq_top (n r : ℕ) (h : r ≤ n) :
    gq Q (1 * Q ^ (n + r + 1)) (n - r) = poch Q (2 * n) * iq Q (n + r) := by
  have := poch_mul_iq_sub Q (show n - r ≤ 2 * n by omega)
  rw [show 2 * n - (n - r) = n + r by omega] at this
  rw [this, gq]
  refine prod_congr rfl (fun i _ => ?_)
  rw [one_mul, ← pow_add]
  congr 2
  omega

/-- `(Q^{n+r+2};Q)_{n-r} = (Q)_{2n+1} / (Q)_{n+r+1}`. -/
lemma gq_top' (n r : ℕ) (h : r ≤ n) :
    gq Q (Q * Q ^ (n + r + 1)) (n - r) = poch Q (2 * n + 1) * iq Q (n + r + 1) := by
  have := poch_mul_iq_sub Q (show n - r ≤ 2 * n + 1 by omega)
  rw [show 2 * n + 1 - (n - r) = n + r + 1 by omega] at this
  rw [this, gq]
  refine prod_congr rfl (fun i _ => ?_)
  rw [← pow_succ', ← pow_add]
  congr 2
  omega

/-- The change of `a` (`Q → 1`) of the Bailey lattice. -/
theorem IsBP.toOne {α β : ℕ → K⟦X⟧} (h : IsBP Q Q α β) (γ : ℕ → K⟦X⟧)
    (hγ : ∀ r, α r * (1 - Q) = γ r * (1 - Q ^ (2 * r + 1))) :
    IsBP Q 1 (fun r => γ r - if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) β := by
  intro n
  have hb := h n
  -- `(Q)_{2n+1} β_n = Σ_r γ_r (1 - Q^{2r+1}) (Q)_{2n+1} / ((Q)_{n-r} (Q)_{n+r+1})`
  have h1 : poch Q (2 * n + 1) * β n
      = poch Q (2 * n + 1) * ∑ r ∈ range (n + 1),
          γ r * (1 - Q ^ (2 * r + 1)) * iq Q (n - r) * iq Q (n + r + 1) := by
    conv_lhs => rw [← gq_QQ, mul_assoc, hb, bps, mul_sum]
    rw [mul_sum]
    refine sum_congr rfl (fun r hr => ?_)
    have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
    rw [gq_top' Q n r hr', ← hγ r]
    ring
  have hne : (1 - Q ^ (2 * n + 1) : K⟦X⟧) ≠ 0 := one_sub_pow_ne_zero Q (2 * n)
  apply mul_left_cancel₀ hne
  rw [gq_one_eq_poch, ← mul_assoc, mul_comm (1 - Q ^ (2 * n + 1)), ← poch_succ, h1, bps]
  simp only [mul_sum]
  -- right side
  have h2 : ∀ r ∈ range (n + 1), (1 - Q ^ (2 * n + 1)) *
      ((γ r - if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) * iq Q (n - r) *
        gq Q (1 * Q ^ (n + r + 1)) (n - r))
      = poch Q (2 * n + 1) * ((γ r - if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) *
          iq Q (n - r) * iq Q (n + r)) := by
    intro r hr
    have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
    rw [gq_top Q n r hr', poch_succ]
    ring
  rw [sum_congr rfl h2, ← mul_sum, ← mul_sum]
  congr 1
  -- the combinatorial identity
  have hsplit : ∑ r ∈ range (n + 1),
      (γ r - if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) * iq Q (n - r) * iq Q (n + r)
      = ∑ r ∈ range (n + 1), γ r * iq Q (n - r) * iq Q (n + r)
        - ∑ r ∈ range n, Q ^ (2 * r + 1) * γ r * iq Q (n - (r + 1)) * iq Q (n + (r + 1)) := by
    rw [show ∑ r ∈ range (n + 1),
        (γ r - if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) * iq Q (n - r) * iq Q (n + r)
        = ∑ r ∈ range (n + 1), (γ r * iq Q (n - r) * iq Q (n + r)
          - (if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) * iq Q (n - r) * iq Q (n + r))
        from sum_congr rfl (fun r _ => by ring), sum_sub_distrib]
    rw [sum_range_succ' (fun r => (if r = 0 then 0 else Q ^ (2 * r - 1) * γ (r - 1)) *
        iq Q (n - r) * iq Q (n + r))]
    simp only [if_neg (Nat.succ_ne_zero _), if_true, zero_mul, add_zero, Nat.add_sub_cancel]
    congr 1
  rw [hsplit, sum_range_succ, sum_range_succ (fun r => γ r * iq Q (n - r) * iq Q (n + r))]
  have hlast : γ n * (1 - Q ^ (2 * n + 1)) * iq Q (n - n) * iq Q (n + n + 1)
      = γ n * iq Q (n - n) * iq Q (n + n) := by
    rw [Nat.sub_self, ← iq_succ_mul Q (n + n), show 2 * n + 1 = n + n + 1 by ring]
    ring
  rw [hlast]
  have hterm : ∀ r ∈ range n, γ r * (1 - Q ^ (2 * r + 1)) * iq Q (n - r) * iq Q (n + r + 1)
      = γ r * iq Q (n - r) * iq Q (n + r)
        - Q ^ (2 * r + 1) * γ r * iq Q (n - (r + 1)) * iq Q (n + (r + 1)) := by
    intro r hr
    have hr' : r < n := mem_range.mp hr
    have e1 : iq Q (n - (r + 1)) = iq Q (n - r) * (1 - Q ^ (n - r)) := by
      rw [show n - r = (n - (r + 1)) + 1 by omega, iq_succ_mul]
    have e2 : iq Q (n + r) = iq Q (n + r + 1) * (1 - Q ^ (n + r + 1)) := (iq_succ_mul Q _).symm
    have e3 : Q ^ (n + r + 1) = Q ^ (2 * r + 1) * Q ^ (n - r) := by
      rw [← pow_add]; congr 1; omega
    rw [e1, e2, show n + (r + 1) = n + r + 1 by ring, e3]
    ring
  rw [sum_congr rfl hterm, sum_sub_distrib]
  ring

end Nahm
