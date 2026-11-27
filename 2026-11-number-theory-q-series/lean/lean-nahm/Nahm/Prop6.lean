/-
Nahm/Prop6.lean — Proposition 6 of `proofs/E7.md`: the finite identity `c_n = 1 + q^{n²}` (`n ≥ 1`),
where (E7.md (3.1))

  `c_n = (1 - q^n) Σ_{r, j ≥ 0, 2r + j ≤ n} (-1)^r q^{(3r²-r)/2 + j² + 2rj} (q)_{n-r-1} / ((q)_r (q)_j (q)_{n-2r-j})`,

over an arbitrary base `Q` with zero constant term (`(3r²-r)/2 = r² + r(r-1)/2`).
-/
import Nahm.Bailey
import KRTW.JTP
import KRTW.Stage2

namespace Nahm

open PowerSeries Finset KRTW

variable {K : Type*} [Field K] (Q : K⟦X⟧) [Nil Q]

omit [Nil Q] in
/-- The inner sum `S_u` of Proposition 6 (after substituting `u = r + j`). -/
noncomputable def p6S (n u : ℕ) : K⟦X⟧ :=
  ∑ r ∈ range (n + 1),
    (if r ≤ u ∧ u + r ≤ n then
      (-1) ^ r * Q ^ (tri (r - 1)) * poch Q (n - r - 1) * iq Q r * iq Q (u - r) * iq Q (n - u - r)
    else 0)

omit [Nil Q] in
lemma p6_step1 (n : ℕ) :
    ∑ r ∈ range (n + 1), ∑ j ∈ range (n + 1),
      (if 2 * r + j ≤ n then
        (-1) ^ r * Q ^ (r ^ 2 + tri (r - 1) + j ^ 2 + 2 * r * j) * (1 - Q ^ n) * poch Q (n - r - 1)
          * iq Q r * iq Q j * iq Q (n - 2 * r - j)
      else 0)
      = ∑ u ∈ range (n + 1), Q ^ (u ^ 2) * (1 - Q ^ n) * p6S Q n u := by
  simp only [p6S, mul_sum]
  conv_rhs => rw [sum_comm]
  refine sum_congr rfl (fun r _ => ?_)
  simp only [mul_ite, mul_zero]
  rw [← sum_filter, ← sum_filter]
  refine sum_nbij' (· + r) (· - r) ?_ ?_ ?_ ?_ ?_
  · intro a ha; simp only [mem_filter, mem_range] at ha ⊢; omega
  · intro a ha; simp only [mem_filter, mem_range] at ha ⊢; omega
  · intro a _; simp
  · intro a ha; simp only [mem_filter, mem_range] at ha ⊢; omega
  · intro a ha
    simp only [mem_filter, mem_range] at ha
    rw [Nat.add_sub_cancel, show n - (a + r) - r = n - 2 * r - a by omega,
      show r ^ 2 + tri (r - 1) + a ^ 2 + 2 * r * a = (a + r) ^ 2 + tri (r - 1) by ring, pow_add]
    ring

omit [Nil Q] in
lemma p6S_zero (n : ℕ) : p6S Q n 0 = poch Q (n - 1) * iq Q n := by
  unfold p6S
  rw [sum_eq_single 0]
  · simp
  · intro b _ hb; rw [if_neg (by omega)]
  · intro h; simp at h

omit [Nil Q] in
lemma p6S_self (n : ℕ) : p6S Q n n = poch Q (n - 1) * iq Q n := by
  unfold p6S
  rw [sum_eq_single 0]
  · simp
  · intro b _ hb; rw [if_neg (by omega)]
  · intro h; simp at h

omit [Nil Q] in
lemma p6S_symm {n u : ℕ} (hu : u ≤ n) : p6S Q n u = p6S Q n (n - u) := by
  unfold p6S
  refine sum_congr rfl (fun r _ => ?_)
  by_cases h : r ≤ u ∧ u + r ≤ n
  · rw [if_pos h, if_pos (by omega), Nat.sub_sub_self hu]; ring
  · rw [if_neg h, if_neg (by omega)]

lemma p6_tri {r v : ℕ} (h : r ≤ v) : tri (r - 1) + v = tri r + (v - r) := by
  cases r with
  | zero => simp
  | succ k => rw [tri_succ]; simp only [Nat.add_sub_cancel]; omega

lemma p6_key (v k : ℕ) (hk1 : 1 ≤ k) (hk : k ≤ v) :
    ∑ r ∈ range (v + 1), (-1) ^ r * Q ^ (tri r) * qbinom Q v r * (Q ^ k) ^ (v - r) = 0 := by
  have h := qbinomial Q (Q ^ k) (-1) v
  have h0 : ∏ i ∈ range v, (Q ^ k + (-1) * Q ^ (i + 1)) = 0 := by
    apply prod_eq_zero (i := k - 1) (mem_range.mpr (by omega))
    rw [Nat.sub_add_cancel hk1]; ring
  rw [h0] at h
  rw [h]
  refine sum_congr rfl (fun r _ => ?_)
  simp only [qbt]; ring

lemma p6_A {n v : ℕ} (hv1 : 1 ≤ v) (hv : 2 * v ≤ n) : Q ^ v * poch Q v * p6S Q n v = 0 := by
  have hS : p6S Q n v = ∑ r ∈ range (v + 1),
      (-1) ^ r * Q ^ (tri (r - 1)) * poch Q (n - r - 1) * iq Q r * iq Q (v - r) * iq Q (n - v - r) := by
    unfold p6S
    rw [← sum_subset (s₁ := range (v + 1)) (fun x hx => by simp at hx ⊢; omega)]
    · refine sum_congr rfl (fun r hr => ?_)
      simp only [mem_range] at hr
      rw [if_pos (by omega)]
    · intro x _ hx
      simp only [mem_range] at hx
      rw [if_neg (by omega)]
  rw [hS, mul_sum]
  have hterm : ∀ r ∈ range (v + 1),
      Q ^ v * poch Q v * ((-1) ^ r * Q ^ (tri (r - 1)) * poch Q (n - r - 1) * iq Q r * iq Q (v - r)
        * iq Q (n - v - r))
      = ∑ t ∈ (range (v - 1)).powerset, (∏ s ∈ t, -Q ^ (n - 2 * v + s + 1))
          * ((-1) ^ r * Q ^ (tri r) * qbinom Q v r * (Q ^ (t.card + 1)) ^ (v - r)) := by
    intro r hr
    simp only [mem_range] at hr
    have e1 : poch Q (n - r - 1) * iq Q (n - v - r)
        = ∏ s ∈ range (v - 1), (1 - Q ^ (n - 2 * v + s + 1) * Q ^ (v - r)) := by
      have := poch_mul_iq_sub Q (n := n - r - 1) (m := v - 1) (by omega)
      rw [show n - r - 1 - (v - 1) = n - v - r by omega] at this
      rw [this]
      refine prod_congr rfl (fun s _ => ?_)
      rw [← pow_add]; congr 2; omega
    have e2 : qbinom Q v r = poch Q v * iq Q r * iq Q (v - r) := by
      simp [qbinom, show r ≤ v by omega]
    have e3 : Q ^ v * Q ^ (tri (r - 1)) = Q ^ (tri r) * Q ^ (v - r) := by
      rw [← pow_add, ← pow_add, add_comm v, p6_tri (by omega)]
    have e4 : ∏ s ∈ range (v - 1), (1 - Q ^ (n - 2 * v + s + 1) * Q ^ (v - r))
        = ∑ t ∈ (range (v - 1)).powerset,
            (∏ s ∈ t, -Q ^ (n - 2 * v + s + 1)) * (Q ^ (v - r)) ^ t.card := by
      rw [show ∏ s ∈ range (v - 1), (1 - Q ^ (n - 2 * v + s + 1) * Q ^ (v - r))
          = ∏ s ∈ range (v - 1), ((-Q ^ (n - 2 * v + s + 1)) * Q ^ (v - r) + 1) from
          prod_congr rfl (fun _ _ => by ring), prod_add]
      refine sum_congr rfl (fun t _ => ?_)
      rw [prod_mul_distrib, prod_const, prod_const_one, mul_one]
    calc Q ^ v * poch Q v * ((-1) ^ r * Q ^ (tri (r - 1)) * poch Q (n - r - 1) * iq Q r
          * iq Q (v - r) * iq Q (n - v - r))
        = (-1) ^ r * (Q ^ v * Q ^ (tri (r - 1))) * (poch Q v * iq Q r * iq Q (v - r))
          * (poch Q (n - r - 1) * iq Q (n - v - r)) := by ring
      _ = _ := by
        rw [e3, ← e2, e1, e4, mul_sum]
        refine sum_congr rfl (fun t _ => ?_)
        have hp : (Q ^ (t.card + 1)) ^ (v - r) = Q ^ (v - r) * (Q ^ (v - r)) ^ t.card := by
          rw [← pow_mul, ← pow_mul, ← pow_add]; congr 1; ring
        rw [hp]
        ring
  rw [sum_congr rfl hterm, sum_comm]
  apply sum_eq_zero
  intro t ht
  have hc : t.card ≤ v - 1 := by
    have := card_le_card (mem_powerset.mp ht)
    rwa [card_range] at this
  rw [← mul_sum, p6_key Q v (t.card + 1) (by omega) (by omega), mul_zero]

lemma p6_middle {n u : ℕ} (h1 : 1 ≤ u) (h2 : u + 1 ≤ n) :
    Q ^ (u ^ 2) * (1 - Q ^ n) * p6S Q n u = 0 := by
  obtain ⟨v, hv1, hv2, hvu, hS⟩ : ∃ v, 1 ≤ v ∧ 2 * v ≤ n ∧ v ≤ u ∧ p6S Q n u = p6S Q n v := by
    by_cases h : 2 * u ≤ n
    · exact ⟨u, h1, h, le_rfl, rfl⟩
    · exact ⟨n - u, by omega, by omega, by omega, p6S_symm Q (by omega)⟩
  have hA := p6_A Q hv1 hv2
  have e : Q ^ (u ^ 2) = Q ^ (u ^ 2 - v) * Q ^ v := by
    rw [← pow_add, Nat.sub_add_cancel]; nlinarith
  calc Q ^ (u ^ 2) * (1 - Q ^ n) * p6S Q n u
      = Q ^ (u ^ 2 - v) * (1 - Q ^ n) * (Q ^ v * poch Q v * p6S Q n v) * iq Q v := by
        rw [hS, e]
        linear_combination (-(Q ^ (u ^ 2 - v) * Q ^ v * (1 - Q ^ n) * p6S Q n v)) * poch_mul_iq Q v
    _ = 0 := by rw [hA]; ring

/-- **Proposition 6 of `E7.md`.** -/
theorem prop6 (n : ℕ) (hn : 1 ≤ n) :
    ∑ r ∈ range (n + 1), ∑ j ∈ range (n + 1),
      (if 2 * r + j ≤ n then
        (-1) ^ r * Q ^ (r ^ 2 + tri (r - 1) + j ^ 2 + 2 * r * j) * (1 - Q ^ n) * poch Q (n - r - 1)
          * iq Q r * iq Q j * iq Q (n - 2 * r - j)
      else 0)
      = 1 + Q ^ (n ^ 2) := by
  rw [p6_step1]
  have hend : (1 - Q ^ n) * (poch Q (n - 1) * iq Q n) = 1 := by
    obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
    have h := iq_mul_poch Q (m + 1)
    rw [poch_succ] at h
    rw [Nat.add_sub_cancel]
    linear_combination h
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [sum_range_succ, sum_range_succ', sum_eq_zero (fun u hu => p6_middle Q (by omega)
      (by simp only [mem_range] at hu; omega)), p6S_zero, p6S_self]
  linear_combination (1 + Q ^ ((m + 1) ^ 2)) * hend

end Nahm
