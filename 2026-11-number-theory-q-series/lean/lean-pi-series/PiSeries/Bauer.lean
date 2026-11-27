/-
# Bauer's series  `∑_{n≥0} (−1)ⁿ C(2n,n)³ (4n+1) / 2^{6n} = 2/π`

G. Bauer, *Von den Coefficienten der Reihen von Kugelfunctionen einer Variablen*, J. reine angew.
Math. 56 (1859), 101–121.  First WZ proof: S. B. Ekhad and D. Zeilberger, *A WZ proof of
Ramanujan's formula for π* (1994, arXiv:math/9306213), via a terminating identity and Carlson's
theorem.

The series converges only conditionally (its terms have size `≍ n^{-1/2}`), so the statement is the
convergence of the partial sums `∑_{n<N}`.

WZ proof used here.  With `h(0,0) = 1` and

  `h(n,k+1)/h(n,k) = σ(n,k) = (2k+2n+1)(2k+1) / (4 (k+n+1)²)`,
  `h(n+1,k)/h(n,k) = ρ(n,k) = −(2n+1)² (2k+2n+1) / (8 (n+1) (k+n+1)²)`,
  `Q(n,k) = (4n+2k+1)/2`,

the pair `F = n h`, `G = Q h` satisfies `F(n+1,k) − F(n,k) = G(n,k+1) − G(n,k)`.  In binomial form
`h(n,k) = (−1)ⁿ C(2n,n)² C(2n+2k,n+k) C(2k,k) / (2^{6n} 2^{4k} C(n+k,n))`: the family printed
(without certificate) in J. Guillera, *Generators of some Ramanujan formulas*, Ramanujan J. 11
(2006), section "Some more series".  Here `h(n,0) = (−1)ⁿ C(2n,n)³/2^{6n}`, `F(0,k) = 0`,
`h(0,K) = (C(2K,K)/4^K)²` and `G(0,K) = (2K+1)/2 · h(0,K) → 1/π` by Wallis' product.
Since `|h(N,0)| = (C(2N,N)/4^N)³ ≤ (2/(2N+1))^{3/2}`, also `N ∑_k h(N,k) → 0`.
-/
import PiSeries.WZSigned
import PiSeries.Wallis

open Finset Filter Topology

namespace PiSeries.Bauer

/-- `h(n,k+1)/h(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ :=
  (2 * (k : ℝ) + 2 * n + 1) * (2 * (k : ℝ) + 1) / (4 * ((k : ℝ) + n + 1) ^ 2)

/-- `h(n+1,k)/h(n,k)` -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  -((2 * (n : ℝ) + 1) ^ 2 * (2 * (k : ℝ) + 2 * n + 1)) /
    (8 * ((n : ℝ) + 1) * ((k : ℝ) + n + 1) ^ 2)

/-- the WZ certificate: `G = Q h` -/
noncomputable def Q (n k : ℕ) : ℝ := (4 * (n : ℝ) + 2 * k + 1) / 2

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

lemma σ_le (n : ℕ) (hn : 1 ≤ n) (k : ℕ) :
    σ n k ≤ (((k : ℝ) + 1) / ((k : ℝ) + 1 + 1)) ^ 2 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  unfold σ
  rw [div_pow, div_le_div_iff₀ (by positivity) (by positivity)]
  push_cast
  have key : ((k : ℝ) + 1) ^ 2 * (4 * ((k : ℝ) + ((m : ℝ) + 1) + 1) ^ 2)
      - (2 * (k : ℝ) + 2 * ((m : ℝ) + 1) + 1) * (2 * (k : ℝ) + 1) * ((k : ℝ) + 1 + 1) ^ 2
      = 4 * (m : ℝ) ^ 2 * (k : ℝ) ^ 2 + 8 * (m : ℝ) ^ 2 * (k : ℝ) + 4 * (m : ℝ) ^ 2
        + 4 * (m : ℝ) * (k : ℝ) ^ 3 + 14 * (m : ℝ) * (k : ℝ) ^ 2 + 16 * (m : ℝ) * (k : ℝ)
        + 8 * (m : ℝ) + (k : ℝ) ^ 2 + 4 * (k : ℝ) + 4 := by ring
  have hpos : (0 : ℝ) ≤ 4 * (m : ℝ) ^ 2 * (k : ℝ) ^ 2 + 8 * (m : ℝ) ^ 2 * (k : ℝ) + 4 * (m : ℝ) ^ 2
        + 4 * (m : ℝ) * (k : ℝ) ^ 3 + 14 * (m : ℝ) * (k : ℝ) ^ 2 + 16 * (m : ℝ) * (k : ℝ)
        + 8 * (m : ℝ) + (k : ℝ) ^ 2 + 4 * (k : ℝ) + 4 := by positivity
  linarith

lemma ρ0_eq (n : ℕ) : ρ n 0 = -((2 * (n : ℝ) + 1) ^ 3 / (8 * ((n : ℝ) + 1) ^ 3)) := by
  unfold ρ
  push_cast
  have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
  field_simp
  ring

lemma Q_abs_le (n : ℕ) (_hn : 1 ≤ n) (k : ℕ) : |Q n k| ≤ (2 * (n : ℝ) + 1) * ((k : ℝ) + 1) := by
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  unfold Q
  rw [abs_of_nonneg (by positivity)]
  nlinarith [mul_nonneg hn hk]

/-- the boundary row: `h(0,K) = ∏_{j<K} (2j+1)²/(4(j+1)²)` -/
lemma hterm_zero_left_eq (K : ℕ) : hterm σ ρ 0 K = wallisTerm K := by
  rw [hterm_zero_left]
  unfold wallisTerm σ
  refine Finset.prod_congr rfl (fun j _ => ?_)
  push_cast
  ring

/-- the boundary limit `G(0,K) → 1/π` -/
lemma tendsto_boundary :
    Tendsto (fun K => Q 0 K * hterm σ ρ 0 K) atTop (𝓝 (1 / Real.pi)) := by
  have h := tendsto_wallisTerm.const_mul (1 / 2 : ℝ)
  have e : (1 / 2 : ℝ) * (2 / Real.pi) = 1 / Real.pi := by ring
  rw [e] at h
  refine h.congr (fun K => ?_)
  rw [hterm_zero_left_eq]
  unfold Q
  push_cast
  ring

/-- `h(n,0)² = ((C(2n,n)/4ⁿ)²)³` -/
lemma hterm_zero_right_sq (n : ℕ) : (hterm σ ρ n 0) ^ 2 = wallisTerm n ^ 3 := by
  induction n with
  | zero => simp [hterm, wallisTerm]
  | succ n ih =>
    rw [hterm_succ_zero, mul_pow, ih, wallisTerm_succ, mul_pow, ρ0_eq]
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    congr 1
    field_simp
    ring

/-- `(C(2n,n)/4ⁿ)² ≤ 2/(2n+1)` (from Wallis' inequalities) -/
lemma wallisTerm_le (n : ℕ) : wallisTerm n ≤ 2 / (2 * (n : ℝ) + 1) := by
  have hW := Real.Wallis.le_W n
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hW2 : (1 / 2 : ℝ) ≤ Real.Wallis.W n := by
    have h1 : (1 / 2 : ℝ) ≤ (2 * (n : ℝ) + 1) / (2 * n + 2) := by
      rw [div_le_div_iff₀ (by norm_num) (by positivity)]; linarith
    have h2 : (1 : ℝ) ≤ Real.pi / 2 := by linarith [Real.pi_gt_three]
    calc (1 / 2 : ℝ) = 1 / 2 * 1 := by ring
      _ ≤ (2 * (n : ℝ) + 1) / (2 * n + 2) * (Real.pi / 2) :=
          mul_le_mul h1 h2 (by norm_num) (by positivity)
      _ ≤ Real.Wallis.W n := hW
  have hmul := wallisTerm_mul_W n
  have hpos : (0 : ℝ) < 2 * (n : ℝ) + 1 := by positivity
  have hwt : 0 ≤ wallisTerm n := by
    unfold wallisTerm; exact prod_nonneg fun j _ => by positivity
  rw [le_div_iff₀ hpos]
  have h3 : wallisTerm n * (2 * (n : ℝ) + 1) * (1 / 2) ≤ (2 * (n : ℝ) + 1) * wallisTerm n * Real.Wallis.W n := by
    rw [mul_comm (wallisTerm n)]
    exact mul_le_mul_of_nonneg_left hW2 (by positivity)
  linarith

/-- `N · |h(N,0)| → 0` -/
lemma tendsto_mul_abs_hterm :
    Tendsto (fun n : ℕ => (n : ℝ) * |hterm σ ρ n 0|) atTop (𝓝 0) := by
  have hx : Tendsto (fun n : ℕ => Real.sqrt (2 * (1 / ((n : ℝ) + 1)))) atTop (𝓝 0) := by
    have h1 : Tendsto (fun n : ℕ => 2 * (1 / ((n : ℝ) + 1))) atTop (𝓝 (2 * 0)) :=
      (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul 2
    rw [mul_zero] at h1
    have h2 := (Real.continuous_sqrt.tendsto 0).comp h1
    rw [Real.sqrt_zero] at h2
    exact h2
  refine squeeze_zero (fun n => by positivity) (fun n => ?_) hx
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hwt0 : 0 ≤ wallisTerm n := by
    unfold wallisTerm; exact prod_nonneg fun j _ => by positivity
  have hwt := wallisTerm_le n
  have hpos : (0 : ℝ) < 2 * (n : ℝ) + 1 := by positivity
  have hsq : ((n : ℝ) * |hterm σ ρ n 0|) ^ 2 ≤ 2 * (1 / ((n : ℝ) + 1)) := by
    rw [mul_pow, sq_abs, hterm_zero_right_sq]
    have h3 : wallisTerm n ^ 3 ≤ (2 / (2 * (n : ℝ) + 1)) ^ 3 := pow_le_pow_left₀ hwt0 hwt 3
    have h4 : (n : ℝ) ^ 2 * (2 / (2 * (n : ℝ) + 1)) ^ 3 ≤ 2 * (1 / ((n : ℝ) + 1)) := by
      rw [div_pow, mul_one_div, mul_div_assoc', div_le_div_iff₀ (by positivity) (by positivity)]
      nlinarith [mul_nonneg hn hn, mul_nonneg (mul_nonneg hn hn) hn,
        mul_nonneg (mul_nonneg hn hn) (mul_nonneg hn hn)]
    exact (mul_le_mul_of_nonneg_left h3 (by positivity)).trans h4
  have := Real.abs_le_sqrt hsq
  rwa [abs_of_nonneg (by positivity)] at this

/-- `h(n,0) = (−1)ⁿ C(2n,n)³ / 2^{6n}` -/
lemma hterm_zero_right_eq (n : ℕ) :
    hterm σ ρ n 0 = (-1 : ℝ) ^ n * (Nat.choose (2 * n) n : ℝ) ^ 3 / 2 ^ (6 * n) := by
  induction n with
  | zero => simp [hterm]
  | succ n ih =>
    rw [hterm_succ_zero, ih, ρ0_eq]
    have hstep := choose_two_mul_succ n
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have hc : (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)
        = 2 * (2 * (n : ℝ) + 1) * (Nat.choose (2 * n) n : ℝ) / ((n : ℝ) + 1) := by
      rw [eq_div_iff h1, mul_comm]; exact hstep
    rw [hc, show 6 * (n + 1) = 6 * n + 6 by ring, pow_add (2 : ℝ), pow_succ (-1 : ℝ) n]
    field_simp
    ring

/-- the WZ limit: `∑_{n<N} G(n,0) → 1/π` -/
theorem tendsto_sum_G0 :
    Tendsto (fun N => ∑ n ∈ range N, Q n 0 * hterm σ ρ n 0) atTop (𝓝 (1 / Real.pi)) := by
  refine tendsto_of_wz_signed (hterm σ ρ) Q (fun n => |hterm σ ρ n 0|)
    (fun n => 2 * (n : ℝ) + 1) _ ?_ ?_ Q_abs_le tendsto_mul_abs_hterm tendsto_boundary
  · intro n k
    rw [hterm_succ_left compat, hterm_succ_right]
    linear_combination (hterm σ ρ n k) * wz_rat n k
  · intro n hn k
    exact abs_hterm_le_of_sigma_le n (σ_pos n) (σ_le n hn) k

/-- **Bauer (1859):** `∑_{n<N} (−1)ⁿ C(2n,n)³ (4n+1) / 2^{6n} → 2/π` as `N → ∞`. -/
theorem bauer_alternating :
    Tendsto (fun N : ℕ => ∑ n ∈ range N,
        (-1 : ℝ) ^ n * (Nat.choose (2 * n) n : ℝ) ^ 3 * (4 * n + 1) / 2 ^ (6 * n))
      atTop (𝓝 (2 / Real.pi)) := by
  have h := tendsto_sum_G0.const_mul 2
  have e : (2 : ℝ) * (1 / Real.pi) = 2 / Real.pi := by ring
  rw [e] at h
  refine h.congr (fun N => ?_)
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  rw [hterm_zero_right_eq]
  unfold Q
  push_cast
  ring

end PiSeries.Bauer
