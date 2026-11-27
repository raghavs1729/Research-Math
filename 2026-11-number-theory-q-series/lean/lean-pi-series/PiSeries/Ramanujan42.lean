/-
# Ramanujan's series  `∑_{n≥0} C(2n,n)³ (42n+5) / 2^{12n} = 16/π`

S. Ramanujan, *Modular equations and approximations to π*, Quart. J. Math. 45 (1914), eq. (29).
First WZ proof: J. Guillera, *Some binomial series obtained by the WZ-method*, Adv. Appl. Math. 29
(2002) (arXiv:math/0503345), by a different pair and Carlson's theorem.

WZ proof used here.  With `h(0,0) = 1` and

  `h(n,k+1)/h(n,k) = σ(n,k) = (2k+2n+1)² / (4 (k+2n+1)²)`,
  `h(n+1,k)/h(n,k) = ρ(n,k) = (2n+1)³ (2k+2n+1)² / (128 (n+1) (k+2n+1)² (k+2n+2)²)`,
  `Q(n,k) = (16k³ + (104n+36)k² + (224n²+160n+24)k + 168n³+188n²+62n+5) / (32 (k+2n+1)²)`,

the pair `F = n h`, `G = Q h` satisfies `F(n+1,k) − F(n,k) = G(n,k+1) − G(n,k)`.
Here `Q(n,0) = (42n+5)/32`, `h(n,0) = C(2n,n)³/2^{12n}`, `F(0,k) = 0`,
`h(0,K) = (C(2K,K)/4^K)²` and `G(0,K) → 1/(2π)` by Wallis' product.
(`h(n,k) = (1/2)_n (1/2)_{n+k}² / (4^{2n} n! (2n+k)!²)` up to normalisation: the term of Gauss'
`₂F₁(1/2 − n, 1/2 − n; 1; 1)` at index `k + 2n`.)
-/
import PiSeries.Term
import PiSeries.Wallis

open Finset Filter Topology

namespace PiSeries.Ramanujan42

/-- `h(n,k+1)/h(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ :=
  (2 * (k : ℝ) + 2 * n + 1) ^ 2 / (4 * ((k : ℝ) + 2 * n + 1) ^ 2)

/-- `h(n+1,k)/h(n,k)` -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  (2 * (n : ℝ) + 1) ^ 3 * (2 * (k : ℝ) + 2 * n + 1) ^ 2 /
    (128 * ((n : ℝ) + 1) * ((k : ℝ) + 2 * n + 1) ^ 2 * ((k : ℝ) + 2 * n + 2) ^ 2)

/-- the WZ certificate: `G = Q h` -/
noncomputable def Q (n k : ℕ) : ℝ :=
  (16 * (k : ℝ) ^ 3 + (104 * n + 36) * (k : ℝ) ^ 2 + (224 * (n : ℝ) ^ 2 + 160 * n + 24) * k
      + 168 * (n : ℝ) ^ 3 + 188 * (n : ℝ) ^ 2 + 62 * n + 5) / (32 * ((k : ℝ) + 2 * n + 1) ^ 2)

lemma compat (n k : ℕ) : ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1) := by
  unfold ρ σ
  push_cast
  have h1 : ((k : ℝ) + 2 * n + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + 2 * n + 2) ≠ 0 := by positivity
  have h3 : ((k : ℝ) + 2 * (n + 1) + 1) ≠ 0 := by positivity
  have h4 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h5 : ((k : ℝ) + 1 + 2 * n + 1) ≠ 0 := by positivity
  have h6 : ((k : ℝ) + 1 + 2 * n + 2) ≠ 0 := by positivity
  field_simp
  ring

/-- **The WZ equation** (divided by `h(n,k)`), an identity between rational functions. -/
lemma wz_rat (n k : ℕ) : ((n : ℝ) + 1) * ρ n k - n = Q n (k + 1) * σ n k - Q n k := by
  unfold ρ σ Q
  push_cast
  have h1 : ((k : ℝ) + 2 * n + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + 2 * n + 2) ≠ 0 := by positivity
  have h4 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h5 : ((k : ℝ) + 1 + 2 * n + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma σ_pos (n k : ℕ) : 0 < σ n k := by unfold σ; positivity

lemma ρ_pos (n k : ℕ) : 0 < ρ n k := by unfold ρ; positivity

lemma σ_le (n : ℕ) (hn : 1 ≤ n) (k : ℕ) :
    σ n k ≤ (((k : ℝ) + 2 * n) / ((k : ℝ) + 2 * n + 1)) ^ 2 := by
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  unfold σ
  rw [div_pow, div_le_div_iff₀ (by positivity) (by positivity)]
  have e : (2 * (k : ℝ) + 2 * n + 1) ^ 2 ≤ 4 * ((k : ℝ) + 2 * n) ^ 2 := by nlinarith
  calc (2 * (k : ℝ) + 2 * n + 1) ^ 2 * ((k : ℝ) + 2 * n + 1) ^ 2
      ≤ 4 * ((k : ℝ) + 2 * n) ^ 2 * ((k : ℝ) + 2 * n + 1) ^ 2 :=
        mul_le_mul_of_nonneg_right e (by positivity)
    _ = ((k : ℝ) + 2 * n) ^ 2 * (4 * ((k : ℝ) + 2 * n + 1) ^ 2) := by ring

lemma ρ0_eq (n : ℕ) : ρ n 0 = (2 * (n : ℝ) + 1) ^ 3 / (512 * ((n : ℝ) + 1) ^ 3) := by
  unfold ρ
  push_cast
  have h1 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h3 : (2 * (n : ℝ) + 2) ≠ 0 := by positivity
  field_simp
  ring

lemma ρ0_le (n : ℕ) : ρ n 0 ≤ 1 / 64 := by
  rw [ρ0_eq]
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith [mul_nonneg hn hn, mul_nonneg (mul_nonneg hn hn) hn]

lemma Q_nonneg (n k : ℕ) : 0 ≤ Q n k := by unfold Q; positivity

lemma Q_abs_le (n : ℕ) (_hn : 1 ≤ n) (k : ℕ) : |Q n k| ≤ (2 * (n : ℝ) + 1) * ((k : ℝ) + 1) := by
  rw [abs_of_nonneg (Q_nonneg n k)]
  unfold Q
  rw [div_le_iff₀ (by positivity)]
  have key : (2 * (n : ℝ) + 1) * ((k : ℝ) + 1) * (32 * ((k : ℝ) + 2 * n + 1) ^ 2)
      - (16 * (k : ℝ) ^ 3 + (104 * n + 36) * (k : ℝ) ^ 2 + (224 * (n : ℝ) ^ 2 + 160 * n + 24) * k
        + 168 * (n : ℝ) ^ 3 + 188 * (n : ℝ) ^ 2 + 62 * n + 5)
      = 256 * (n : ℝ) ^ 3 * (k : ℝ) + 88 * (n : ℝ) ^ 3 + 256 * (n : ℝ) ^ 2 * (k : ℝ) ^ 2
        + 416 * (n : ℝ) ^ 2 * (k : ℝ) + 196 * (n : ℝ) ^ 2 + 64 * (n : ℝ) * (k : ℝ) ^ 3
        + 216 * (n : ℝ) * (k : ℝ) ^ 2 + 288 * (n : ℝ) * (k : ℝ) + 130 * (n : ℝ) + 16 * (k : ℝ) ^ 3
        + 60 * (k : ℝ) ^ 2 + 72 * (k : ℝ) + 27 := by ring
  have hpos : (0 : ℝ) ≤ 256 * (n : ℝ) ^ 3 * (k : ℝ) + 88 * (n : ℝ) ^ 3 + 256 * (n : ℝ) ^ 2 * (k : ℝ) ^ 2
        + 416 * (n : ℝ) ^ 2 * (k : ℝ) + 196 * (n : ℝ) ^ 2 + 64 * (n : ℝ) * (k : ℝ) ^ 3
        + 216 * (n : ℝ) * (k : ℝ) ^ 2 + 288 * (n : ℝ) * (k : ℝ) + 130 * (n : ℝ) + 16 * (k : ℝ) ^ 3
        + 60 * (k : ℝ) ^ 2 + 72 * (k : ℝ) + 27 := by positivity
  linarith

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
    have h1 : Tendsto (fun K : ℕ => 1 / 4 - (2 + 1 / ((K : ℝ) + 1)) * (1 / ((K : ℝ) + 1)) / 32)
        atTop (𝓝 (1 / 4 - (2 + 0) * 0 / 32)) :=
      tendsto_const_nhds.sub (((tendsto_const_nhds.add hx).mul hx).div_const 32)
    have h2 : ((1 / 4 - (2 + 0) * 0 / 32 : ℝ)) = 1 / 4 := by norm_num
    rw [h2] at h1
    refine h1.congr (fun K => ?_)
    unfold Q
    have hK : ((K : ℝ) + 1) ≠ 0 := by positivity
    have hK2 : (2 * (K : ℝ) + 1) ≠ 0 := by positivity
    push_cast
    field_simp
    ring
  have h := tendsto_mul_wallisTerm hq
  have e : (1 / 4 : ℝ) * (2 / Real.pi) = 1 / (2 * Real.pi) := by
    field_simp; norm_num
  rw [e] at h
  simpa [hterm_zero_left_eq] using h

/-- `h(n,0) = C(2n,n)³ / 2^{12n}` -/
lemma hterm_zero_right_eq (n : ℕ) :
    hterm σ ρ n 0 = (Nat.choose (2 * n) n : ℝ) ^ 3 / 2 ^ (12 * n) := by
  induction n with
  | zero => simp [hterm]
  | succ n ih =>
    rw [hterm_succ_zero, ih, ρ0_eq]
    have hstep := choose_two_mul_succ n
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have hc : (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)
        = 2 * (2 * (n : ℝ) + 1) * (Nat.choose (2 * n) n : ℝ) / ((n : ℝ) + 1) := by
      rw [eq_div_iff h1, mul_comm]; exact hstep
    rw [hc, show 12 * (n + 1) = 12 * n + 12 by ring, pow_add]
    field_simp
    ring

lemma Q0_eq (n : ℕ) : Q n 0 = (42 * (n : ℝ) + 5) / 32 := by
  unfold Q
  push_cast
  have h1 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

/-- the WZ sum: `∑_n G(n,0) = 1/(2π)` -/
theorem hasSum_G0 : HasSum (fun n => Q n 0 * hterm σ ρ n 0) (1 / (2 * Real.pi)) :=
  hasSum_of_ratios σ ρ Q (fun n => 2 * (n : ℝ) + 1) 2 (1 / 64) _ compat wz_rat σ_pos
    (fun n => ρ_pos n 0) (by norm_num) σ_le (by norm_num) ρ0_le Q_abs_le (fun n => Q_nonneg n 0)
    tendsto_boundary

/-- **Ramanujan (1914):** `∑_{n≥0} C(2n,n)³ (42n+5) / 2^{12n} = 16/π`. -/
theorem ramanujan_forty_two_n_plus_five :
    HasSum (fun n : ℕ => (Nat.choose (2 * n) n : ℝ) ^ 3 * (42 * n + 5) / 2 ^ (12 * n))
      (16 / Real.pi) := by
  have h := hasSum_G0.mul_left 32
  have e : (32 : ℝ) * (1 / (2 * Real.pi)) = 16 / Real.pi := by
    field_simp; norm_num
  rw [e] at h
  refine h.congr_fun ?_
  intro n
  rw [hterm_zero_right_eq, Q0_eq]
  ring

end PiSeries.Ramanujan42
