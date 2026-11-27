/-
# Ramanujan's series  `∑_{n≥0} C(2n,n)³ (6n+1) / 2^{8n} = 4/π`

S. Ramanujan, *Modular equations and approximations to π*, Quart. J. Math. 45 (1914), eq. (28).

WZ proof.  With `h(n,k) = (1/2)_n³ (1/2)_k² / (4^n n! (n+k)!²)`, i.e. `h(0,0) = 1` and

  `h(n,k+1)/h(n,k) = σ(n,k) = (2k+1)² / (4 (k+n+1)²)`,
  `h(n+1,k)/h(n,k) = ρ(n,k) = (2n+1)³ / (32 (n+1) (k+n+1)²)`,

and `Q(n,k) = (4k+6n+1)/8`, the pair `F = n h`, `G = Q h` satisfies
`F(n+1,k) − F(n,k) = G(n,k+1) − G(n,k)`.  Here `G(n,0) = (6n+1) C(2n,n)³ / (8 · 2^{8n})`,
`F(0,k) = 0`, and `G(0,K) = (4K+1)/8 · (C(2K,K)/4^K)² → 1/(2π)` by Wallis' product.
The pair is the term of Gauss' `₂F₁(1/2 − n, 1/2 − n; 1; 1)` at index `k + n`; in binomial form
`8 G(n,k) = C(2n,n)³ C(2k,k)² (6n+4k+1) / (2^{8n} 2^{4k} C(n+k,n)²)`, which is the family printed in
J. Guillera, *Generators of some Ramanujan formulas*, Ramanujan J. 11 (2006) (arXiv:1104.0392),
section "Some more series" (without certificate).
-/
import PiSeries.Term
import PiSeries.Wallis

open Finset Filter Topology

namespace PiSeries.Ramanujan6

/-- `h(n,k+1)/h(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ := (2 * (k : ℝ) + 1) ^ 2 / (4 * ((k : ℝ) + n + 1) ^ 2)

/-- `h(n+1,k)/h(n,k)` -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  (2 * (n : ℝ) + 1) ^ 3 / (32 * ((n : ℝ) + 1) * ((k : ℝ) + n + 1) ^ 2)

/-- the WZ certificate: `G = Q h` -/
noncomputable def Q (n k : ℕ) : ℝ := (4 * (k : ℝ) + 6 * n + 1) / 8

lemma compat (n k : ℕ) : ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1) := by
  unfold ρ σ
  push_cast
  have h1 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + (n + 1) + 1) ≠ 0 := by positivity
  have h3 : ((k : ℝ) + 1 + n + 1) ≠ 0 := by positivity
  have h4 : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

/-- **The WZ equation** (divided by `h(n,k)`), an identity between rational functions. -/
lemma wz_rat (n k : ℕ) : ((n : ℝ) + 1) * ρ n k - n = Q n (k + 1) * σ n k - Q n k := by
  unfold ρ σ Q
  push_cast
  have h1 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h4 : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma σ_pos (n k : ℕ) : 0 < σ n k := by unfold σ; positivity

lemma ρ_pos (n k : ℕ) : 0 < ρ n k := by unfold ρ; positivity

lemma σ_le (n : ℕ) (hn : 1 ≤ n) (k : ℕ) :
    σ n k ≤ (((k : ℝ) + 1 * n) / ((k : ℝ) + 1 * n + 1)) ^ 2 := by
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  unfold σ
  rw [one_mul, div_pow, div_le_div_iff₀ (by positivity) (by positivity)]
  have e : (2 * (k : ℝ) + 1) ^ 2 ≤ 4 * ((k : ℝ) + n) ^ 2 := by nlinarith
  calc (2 * (k : ℝ) + 1) ^ 2 * ((k : ℝ) + n + 1) ^ 2
      ≤ 4 * ((k : ℝ) + n) ^ 2 * ((k : ℝ) + n + 1) ^ 2 :=
        mul_le_mul_of_nonneg_right e (by positivity)
    _ = ((k : ℝ) + n) ^ 2 * (4 * ((k : ℝ) + n + 1) ^ 2) := by ring

lemma ρ0_le (n : ℕ) : ρ n 0 ≤ 1 / 4 := by
  unfold ρ
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  push_cast
  nlinarith [mul_nonneg hn hn, mul_nonneg (mul_nonneg hn hn) hn]

lemma Q_abs_le (n : ℕ) (_hn : 1 ≤ n) (k : ℕ) : |Q n k| ≤ (6 * (n : ℝ) + 1) * ((k : ℝ) + 1) := by
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  unfold Q
  rw [abs_of_nonneg (by positivity)]
  nlinarith [mul_nonneg hn hk]

lemma Q0_nonneg (n : ℕ) : 0 ≤ Q n 0 := by unfold Q; positivity

/-- the boundary row: `h(0,K) = ∏_{j<K} (2j+1)²/(4(j+1)²)` -/
lemma hterm_zero_left_eq (K : ℕ) : hterm σ ρ 0 K = wallisTerm K := by
  rw [hterm_zero_left]
  unfold wallisTerm σ
  simp

/-- the boundary limit `G(0,K) → 1/(2π)` -/
lemma tendsto_boundary :
    Tendsto (fun K => Q 0 K * hterm σ ρ 0 K) atTop (𝓝 (1 / (2 * Real.pi))) := by
  have hx : Tendsto (fun K : ℕ => 1 / ((K : ℝ) + 1)) atTop (𝓝 0) :=
    tendsto_one_div_add_atTop_nhds_zero_nat
  have hq : Tendsto (fun K : ℕ => Q 0 K / (2 * (K : ℝ) + 1)) atTop (𝓝 (1 / 4)) := by
    have h1 : Tendsto (fun K : ℕ => (4 - 3 * (1 / ((K : ℝ) + 1))) / (8 * (2 - 1 / ((K : ℝ) + 1))))
        atTop (𝓝 ((4 - 3 * 0) / (8 * (2 - 0)))) :=
      (tendsto_const_nhds.sub (hx.const_mul 3)).div
        ((tendsto_const_nhds.sub hx).const_mul 8) (by norm_num)
    have h2 : ((4 - 3 * 0) / (8 * (2 - 0)) : ℝ) = 1 / 4 := by norm_num
    rw [h2] at h1
    refine h1.congr (fun K => ?_)
    unfold Q
    have hK : ((K : ℝ) + 1) ≠ 0 := by positivity
    have hK2 : (2 * (K : ℝ) + 1) ≠ 0 := by positivity
    have hK3 : (2 - 1 / ((K : ℝ) + 1)) ≠ 0 := by
      have : 2 - 1 / ((K : ℝ) + 1) = (2 * (K : ℝ) + 1) / ((K : ℝ) + 1) := by field_simp; ring
      rw [this]; positivity
    rw [div_eq_div_iff (mul_ne_zero (by norm_num) hK3) hK2]
    push_cast
    field_simp
    ring
  have h := tendsto_mul_wallisTerm hq
  have e : (1 / 4 : ℝ) * (2 / Real.pi) = 1 / (2 * Real.pi) := by
    field_simp; norm_num
  rw [e] at h
  simpa [hterm_zero_left_eq] using h

/-- `h(n,0) = C(2n,n)³ / 2^{8n}` -/
lemma hterm_zero_right_eq (n : ℕ) :
    hterm σ ρ n 0 = (Nat.choose (2 * n) n : ℝ) ^ 3 / 2 ^ (8 * n) := by
  induction n with
  | zero => simp [hterm]
  | succ n ih =>
    rw [hterm_succ_zero, ih]
    unfold ρ
    have hstep := choose_two_mul_succ n
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have hc : (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)
        = 2 * (2 * (n : ℝ) + 1) * (Nat.choose (2 * n) n : ℝ) / ((n : ℝ) + 1) := by
      rw [eq_div_iff h1, mul_comm]; exact hstep
    rw [hc, show 8 * (n + 1) = 8 * n + 8 by ring, pow_add]
    push_cast
    field_simp
    ring

/-- the WZ sum: `∑_n G(n,0) = 1/(2π)` -/
theorem hasSum_G0 : HasSum (fun n => Q n 0 * hterm σ ρ n 0) (1 / (2 * Real.pi)) :=
  hasSum_of_ratios σ ρ Q (fun n => 6 * (n : ℝ) + 1) 1 (1 / 4) _ compat wz_rat σ_pos
    (fun n => ρ_pos n 0) le_rfl σ_le (by norm_num) ρ0_le Q_abs_le Q0_nonneg tendsto_boundary

/-- **Ramanujan (1914):** `∑_{n≥0} C(2n,n)³ (6n+1) / 2^{8n} = 4/π`. -/
theorem ramanujan_six_n_plus_one :
    HasSum (fun n : ℕ => (Nat.choose (2 * n) n : ℝ) ^ 3 * (6 * n + 1) / 2 ^ (8 * n))
      (4 / Real.pi) := by
  have h := hasSum_G0.mul_left 8
  have e : (8 : ℝ) * (1 / (2 * Real.pi)) = 4 / Real.pi := by
    field_simp; norm_num
  rw [e] at h
  refine h.congr_fun ?_  
  intro n
  rw [hterm_zero_right_eq]
  unfold Q
  push_cast
  ring

end PiSeries.Ramanujan6
