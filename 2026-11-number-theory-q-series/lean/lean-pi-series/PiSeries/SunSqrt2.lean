/-
# Sun's series  `∑_{k≥0} (360k³+612k²+230k+15) C(3k,k) C(6k,3k)² / ((3k+1)(3k+2) 2^{15k}) = 32√2/π`

Conjectured by Zhi-Wei Sun, *New series involving binomial coefficients* (arXiv:2307.03086),
Conjecture 5.25(i).  The proof is Theorem 4 of `monotonicity/ramanujan-series/RESULTS.md`.

WZ pair: `F = n h`, `G = Q h` with `h(0,0) = 1` and

  `h(n,k+1)/h(n,k) = σ(n,k) = (4k+6n+1)(4k+6n+3) / (16 (k+n+1)(k+3n+1))`,
  `h(n+1,k)/h(n,k) = ρ(n,k) = 3(6n+1)(6n+5)(4k+6n+1)(4k+6n+3)(4k+6n+5)
                               / (1024 (n+1)(k+n+1)(k+3n+1)(k+3n+2)(k+3n+3))`,
  `Q(n,k) = (64k³+(160n+144)k²+(240n²+432n+92)k+360n³+612n²+230n+15) / (64 (k+3n+1)(k+3n+2))`.

`h(n,0) = C(3n,n) C(6n,3n)²/2^{15n}`, `64 Q(n,0) h(n,0)` is Sun's summand, `F(0,k) = 0`, and
`h(0,K) = (1/4)_K (3/4)_K / K!²`, so that `G(0,K) → 1/(Γ(1/4)Γ(3/4)) = sin(π/4)/π = √2/(2π)`.
In this file that limit is obtained from Euler's product for the sine at `x = 1/4`
(Mathlib: `Real.tendsto_euler_sin_prod`), which is the reflection formula in product form.
-/
import PiSeries.Term
import PiSeries.Wallis

open Finset Filter Topology

namespace PiSeries.SunSqrt2

/-- `h(n,k+1)/h(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ :=
  (4 * (k : ℝ) + 6 * n + 1) * (4 * (k : ℝ) + 6 * n + 3) /
    (16 * ((k : ℝ) + n + 1) * ((k : ℝ) + 3 * n + 1))

/-- `h(n+1,k)/h(n,k)` -/
noncomputable def ρ (n k : ℕ) : ℝ :=
  3 * (6 * (n : ℝ) + 1) * (6 * (n : ℝ) + 5) * (4 * (k : ℝ) + 6 * n + 1) * (4 * (k : ℝ) + 6 * n + 3)
      * (4 * (k : ℝ) + 6 * n + 5) /
    (1024 * ((n : ℝ) + 1) * ((k : ℝ) + n + 1) * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2)
      * ((k : ℝ) + 3 * n + 3))

/-- the WZ certificate: `G = Q h` -/
noncomputable def Q (n k : ℕ) : ℝ :=
  (64 * (k : ℝ) ^ 3 + (160 * n + 144) * (k : ℝ) ^ 2 + (240 * (n : ℝ) ^ 2 + 432 * n + 92) * k
      + 360 * (n : ℝ) ^ 3 + 612 * (n : ℝ) ^ 2 + 230 * n + 15) /
    (64 * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2))

lemma compat (n k : ℕ) : ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1) := by
  unfold ρ σ
  push_cast
  have h1 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + 3 * n + 1) ≠ 0 := by positivity
  have h3 : ((k : ℝ) + 3 * n + 2) ≠ 0 := by positivity
  have h4 : ((k : ℝ) + 3 * n + 3) ≠ 0 := by positivity
  have h5 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h6 : ((k : ℝ) + (n + 1) + 1) ≠ 0 := by positivity
  have h7 : ((k : ℝ) + 3 * (n + 1) + 1) ≠ 0 := by positivity
  have h8 : ((k : ℝ) + 1 + n + 1) ≠ 0 := by positivity
  have h9 : ((k : ℝ) + 1 + 3 * n + 1) ≠ 0 := by positivity
  have h10 : ((k : ℝ) + 1 + 3 * n + 2) ≠ 0 := by positivity
  have h11 : ((k : ℝ) + 1 + 3 * n + 3) ≠ 0 := by positivity
  field_simp
  ring

/-- **The WZ equation** (divided by `h(n,k)`), an identity between rational functions. -/
lemma wz_rat (n k : ℕ) : ((n : ℝ) + 1) * ρ n k - n = Q n (k + 1) * σ n k - Q n k := by
  unfold ρ σ Q
  push_cast
  have h1 : ((k : ℝ) + n + 1) ≠ 0 := by positivity
  have h2 : ((k : ℝ) + 3 * n + 1) ≠ 0 := by positivity
  have h3 : ((k : ℝ) + 3 * n + 2) ≠ 0 := by positivity
  have h4 : ((k : ℝ) + 3 * n + 3) ≠ 0 := by positivity
  have h5 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h9 : ((k : ℝ) + 1 + 3 * n + 1) ≠ 0 := by positivity
  have h10 : ((k : ℝ) + 1 + 3 * n + 2) ≠ 0 := by positivity
  field_simp
  ring

lemma σ_pos (n k : ℕ) : 0 < σ n k := by unfold σ; positivity

lemma ρ_pos (n k : ℕ) : 0 < ρ n k := by unfold ρ; positivity

lemma σ_le (n : ℕ) (hn : 1 ≤ n) (k : ℕ) :
    σ n k ≤ (((k : ℝ) + 3 * n) / ((k : ℝ) + 3 * n + 1)) ^ 2 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  unfold σ
  rw [div_pow, div_le_div_iff₀ (by positivity) (by positivity)]
  push_cast
  have key : ((k : ℝ) + 3 * ((m : ℝ) + 1)) ^ 2
        * (16 * ((k : ℝ) + ((m : ℝ) + 1) + 1) * ((k : ℝ) + 3 * ((m : ℝ) + 1) + 1))
      - (4 * (k : ℝ) + 6 * ((m : ℝ) + 1) + 1) * (4 * (k : ℝ) + 6 * ((m : ℝ) + 1) + 3)
        * ((k : ℝ) + 3 * ((m : ℝ) + 1) + 1) ^ 2
      = 108 * (m : ℝ) ^ 4 + 216 * (m : ℝ) ^ 3 * (k : ℝ) + 576 * (m : ℝ) ^ 3
        + 108 * (m : ℝ) ^ 2 * (k : ℝ) ^ 2 + 672 * (m : ℝ) ^ 2 * (k : ℝ) + 1017 * (m : ℝ) ^ 2
        + 16 * (m : ℝ) * (k : ℝ) ^ 3 + 160 * (m : ℝ) * (k : ℝ) ^ 2 + 582 * (m : ℝ) * (k : ℝ)
        + 696 * (m : ℝ) + 17 * (k : ℝ) ^ 2 + 104 * (k : ℝ) + 144 := by ring
  have hpos : (0 : ℝ) ≤ 108 * (m : ℝ) ^ 4 + 216 * (m : ℝ) ^ 3 * (k : ℝ) + 576 * (m : ℝ) ^ 3
        + 108 * (m : ℝ) ^ 2 * (k : ℝ) ^ 2 + 672 * (m : ℝ) ^ 2 * (k : ℝ) + 1017 * (m : ℝ) ^ 2
        + 16 * (m : ℝ) * (k : ℝ) ^ 3 + 160 * (m : ℝ) * (k : ℝ) ^ 2 + 582 * (m : ℝ) * (k : ℝ)
        + 696 * (m : ℝ) + 17 * (k : ℝ) ^ 2 + 104 * (k : ℝ) + 144 := by positivity
  linarith

lemma ρ0_eq (n : ℕ) : ρ n 0 = 3 * (6 * (n : ℝ) + 1) ^ 2 * (6 * (n : ℝ) + 5) ^ 2 * (2 * (n : ℝ) + 1) /
    (1024 * ((n : ℝ) + 1) ^ 3 * (3 * (n : ℝ) + 1) * (3 * (n : ℝ) + 2)) := by
  unfold ρ
  push_cast
  have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h2 : (3 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h3 : (3 * (n : ℝ) + 2) ≠ 0 := by positivity
  have h4 : (3 * (n : ℝ) + 3) ≠ 0 := by positivity
  field_simp
  ring

lemma ρ0_le (n : ℕ) : ρ n 0 ≤ 27 / 32 := by
  rw [ρ0_eq, div_le_div_iff₀ (by positivity) (by positivity)]
  have key : 27 * (1024 * ((n : ℝ) + 1) ^ 3 * (3 * (n : ℝ) + 1) * (3 * (n : ℝ) + 2))
      - 3 * (6 * (n : ℝ) + 1) ^ 2 * (6 * (n : ℝ) + 5) ^ 2 * (2 * (n : ℝ) + 1) * 32
      = 373248 * (n : ℝ) ^ 4 + 981504 * (n : ℝ) ^ 3 + 933120 * (n : ℝ) ^ 2 + 375360 * (n : ℝ)
        + 52896 := by ring
  have hpos : (0 : ℝ) ≤ 373248 * (n : ℝ) ^ 4 + 981504 * (n : ℝ) ^ 3 + 933120 * (n : ℝ) ^ 2
      + 375360 * (n : ℝ) + 52896 := by positivity
  linarith

lemma Q_nonneg (n k : ℕ) : 0 ≤ Q n k := by unfold Q; positivity

lemma Q_abs_le (n : ℕ) (_hn : 1 ≤ n) (k : ℕ) : |Q n k| ≤ (3 * (n : ℝ) + 3) * ((k : ℝ) + 1) := by
  rw [abs_of_nonneg (Q_nonneg n k)]
  unfold Q
  rw [div_le_iff₀ (by positivity)]
  have key : (3 * (n : ℝ) + 3) * ((k : ℝ) + 1) * (64 * ((k : ℝ) + 3 * n + 1) * ((k : ℝ) + 3 * n + 2))
      - (64 * (k : ℝ) ^ 3 + (160 * n + 144) * (k : ℝ) ^ 2 + (240 * (n : ℝ) ^ 2 + 432 * n + 92) * k
        + 360 * (n : ℝ) ^ 3 + 612 * (n : ℝ) ^ 2 + 230 * n + 15)
      = 1728 * (n : ℝ) ^ 3 * (k : ℝ) + 1368 * (n : ℝ) ^ 3 + 1152 * (n : ℝ) ^ 2 * (k : ℝ) ^ 2
        + 4368 * (n : ℝ) ^ 2 * (k : ℝ) + 2844 * (n : ℝ) ^ 2 + 192 * (n : ℝ) * (k : ℝ) ^ 3
        + 1760 * (n : ℝ) * (k : ℝ) ^ 2 + 3408 * (n : ℝ) * (k : ℝ) + 1882 * (n : ℝ)
        + 128 * (k : ℝ) ^ 3 + 624 * (k : ℝ) ^ 2 + 868 * (k : ℝ) + 369 := by ring
  have hpos : (0 : ℝ) ≤ 1728 * (n : ℝ) ^ 3 * (k : ℝ) + 1368 * (n : ℝ) ^ 3 + 1152 * (n : ℝ) ^ 2 * (k : ℝ) ^ 2
        + 4368 * (n : ℝ) ^ 2 * (k : ℝ) + 2844 * (n : ℝ) ^ 2 + 192 * (n : ℝ) * (k : ℝ) ^ 3
        + 1760 * (n : ℝ) * (k : ℝ) ^ 2 + 3408 * (n : ℝ) * (k : ℝ) + 1882 * (n : ℝ)
        + 128 * (k : ℝ) ^ 3 + 624 * (k : ℝ) ^ 2 + 868 * (k : ℝ) + 369 := by positivity
  linarith

/-! ### The boundary row `h(0,K) = (1/4)_K (3/4)_K / K!²` and Euler's sine product -/

/-- `(4K+1) h(0,K) = ∏_{j<K} (1 − (1/4)²/(j+1)²)` -/
lemma boundary_prod (K : ℕ) :
    (4 * (K : ℝ) + 1) * hterm σ ρ 0 K = ∏ j ∈ range K, (1 - (1 / 4 : ℝ) ^ 2 / ((j : ℝ) + 1) ^ 2) := by
  induction K with
  | zero => simp [hterm]
  | succ K ih =>
    rw [hterm_succ_right, prod_range_succ, ← ih]
    unfold σ
    push_cast
    have h1 : ((K : ℝ) + 1) ≠ 0 := by positivity
    field_simp
    ring_nf

/-- **`(4K+1) (1/4)_K (3/4)_K / K!² → 2√2/π`** (Euler's product for `sin (π/4)`). -/
lemma tendsto_boundary_prod :
    Tendsto (fun K : ℕ => (4 * (K : ℝ) + 1) * hterm σ ρ 0 K) atTop (𝓝 (2 * Real.sqrt 2 / Real.pi)) := by
  have h := (Real.tendsto_euler_sin_prod (1 / 4)).const_mul (4 / Real.pi)
  have e : Real.sin (Real.pi * (1 / 4)) = Real.sqrt 2 / 2 := by
    rw [show Real.pi * (1 / 4) = Real.pi / 4 by ring, Real.sin_pi_div_four]
  rw [e] at h
  have e2 : 4 / Real.pi * (Real.sqrt 2 / 2) = 2 * Real.sqrt 2 / Real.pi := by ring
  rw [e2] at h
  refine h.congr (fun K => ?_)
  rw [boundary_prod]
  have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
  field_simp

/-- the boundary limit `G(0,K) → √2/(2π) = 1/(Γ(1/4)Γ(3/4))` -/
lemma tendsto_boundary :
    Tendsto (fun K => Q 0 K * hterm σ ρ 0 K) atTop (𝓝 (Real.sqrt 2 / (2 * Real.pi))) := by
  have hx : Tendsto (fun K : ℕ => 1 / ((K : ℝ) + 1)) atTop (𝓝 0) :=
    tendsto_one_div_add_atTop_nhds_zero_nat
  have hq : Tendsto (fun K : ℕ => Q 0 K / (4 * (K : ℝ) + 1)) atTop (𝓝 (1 / 4)) := by
    have h1 : Tendsto (fun K : ℕ => 1 / 4 - (16 + 1 / ((K : ℝ) + 1)) * (1 / ((K : ℝ) + 1))
          / (64 * (1 + 1 / ((K : ℝ) + 1))))
        atTop (𝓝 (1 / 4 - (16 + 0) * 0 / (64 * (1 + 0)))) :=
      tendsto_const_nhds.sub (((tendsto_const_nhds.add hx).mul hx).div
        ((tendsto_const_nhds.add hx).const_mul 64) (by norm_num))
    have h2 : ((1 / 4 - (16 + 0) * 0 / (64 * (1 + 0)) : ℝ)) = 1 / 4 := by norm_num
    rw [h2] at h1
    refine h1.congr (fun K => ?_)
    unfold Q
    have hK : ((K : ℝ) + 1) ≠ 0 := by positivity
    have hK2 : ((K : ℝ) + 2) ≠ 0 := by positivity
    have hK4 : (4 * (K : ℝ) + 1) ≠ 0 := by positivity
    have hK3 : (1 + 1 / ((K : ℝ) + 1)) ≠ 0 := by positivity
    push_cast
    field_simp
    ring
  have h := hq.mul tendsto_boundary_prod
  have e : (1 / 4 : ℝ) * (2 * Real.sqrt 2 / Real.pi) = Real.sqrt 2 / (2 * Real.pi) := by
    have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
    field_simp; ring
  rw [e] at h
  refine h.congr (fun K => ?_)
  have hK4 : (4 * (K : ℝ) + 1) ≠ 0 := by positivity
  field_simp

/-! ### Identification of `h(n,0)` with the binomial coefficients -/

lemma cast_fact_succ (m : ℕ) : (((m + 1).factorial : ℕ) : ℝ) = ((m : ℝ) + 1) * (m.factorial : ℝ) := by
  push_cast [Nat.factorial_succ]; ring

/-- `(n+1)(2n+1)(2n+2) C(3n+3,n+1) = (3n+1)(3n+2)(3n+3) C(3n,n)` -/
lemma choose_three_step (n : ℕ) :
    (Nat.choose (3 * (n + 1)) (n + 1) : ℝ) * (((n : ℝ) + 1) * (2 * n + 1) * (2 * n + 2))
      = (3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3) * (Nat.choose (3 * n) n : ℝ) := by
  rw [Nat.cast_choose ℝ (by omega : n + 1 ≤ 3 * (n + 1)), Nat.cast_choose ℝ (by omega : n ≤ 3 * n)]
  rw [show 3 * (n + 1) - (n + 1) = 2 * n + 1 + 1 by omega, show 3 * n - n = 2 * n by omega,
    show 3 * (n + 1) = 3 * n + 1 + 1 + 1 by ring]
  simp only [cast_fact_succ]
  have h1 : ((n.factorial : ℕ) : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero _)
  have h2 : (((2 * n).factorial : ℕ) : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero _)
  have h3 : ((n : ℝ) + 1) ≠ 0 := by positivity
  have h4 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h5 : (2 * (n : ℝ) + 1 + 1) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

/-- `(3n+1)(3n+2)(3n+3) C(6n+6,3n+3) = 8 (6n+1)(6n+3)(6n+5) C(6n,3n)` -/
lemma choose_six_step (n : ℕ) :
    (Nat.choose (6 * (n + 1)) (3 * (n + 1)) : ℝ) * ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3))
      = 8 * ((6 * (n : ℝ) + 1) * (6 * n + 3) * (6 * n + 5)) * (Nat.choose (6 * n) (3 * n) : ℝ) := by
  have e1 := choose_two_mul_succ (3 * n)
  have e2 := choose_two_mul_succ (3 * n + 1)
  have e3 := choose_two_mul_succ (3 * n + 2)
  rw [show 6 * (n + 1) = 2 * (3 * n + 2 + 1) by ring, show 3 * (n + 1) = 3 * n + 2 + 1 by ring,
    show 6 * n = 2 * (3 * n) by ring]
  push_cast at e1 e2 e3
  linear_combination (3 * (n : ℝ) + 1) * (3 * n + 2) * e3
    + (3 * (n : ℝ) + 1) * 2 * (2 * (3 * (n : ℝ) + 2) + 1) * e2
    + 2 * (2 * (3 * (n : ℝ) + 2) + 1) * 2 * (2 * (3 * (n : ℝ) + 1) + 1) * e1

/-- `h(n,0) = C(3n,n) C(6n,3n)² / 2^{15n}` -/
lemma hterm_zero_right_eq (n : ℕ) :
    hterm σ ρ n 0
      = (Nat.choose (3 * n) n : ℝ) * (Nat.choose (6 * n) (3 * n) : ℝ) ^ 2 / 2 ^ (15 * n) := by
  induction n with
  | zero => simp [hterm]
  | succ n ih =>
    rw [hterm_succ_zero, ih, ρ0_eq]
    have h1 : ((n : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (2 * (n : ℝ) + 1) ≠ 0 := by positivity
    have h3 : (2 * (n : ℝ) + 2) ≠ 0 := by positivity
    have h4 : (3 * (n : ℝ) + 1) ≠ 0 := by positivity
    have h5 : (3 * (n : ℝ) + 2) ≠ 0 := by positivity
    have h6 : (3 * (n : ℝ) + 3) ≠ 0 := by positivity
    have hc3 : (Nat.choose (3 * (n + 1)) (n + 1) : ℝ)
        = (3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3) * (Nat.choose (3 * n) n : ℝ)
          / (((n : ℝ) + 1) * (2 * n + 1) * (2 * n + 2)) := by
      rw [eq_div_iff (by positivity)]; exact choose_three_step n
    have hc6 : (Nat.choose (6 * (n + 1)) (3 * (n + 1)) : ℝ)
        = 8 * ((6 * (n : ℝ) + 1) * (6 * n + 3) * (6 * n + 5)) * (Nat.choose (6 * n) (3 * n) : ℝ)
          / ((3 * (n : ℝ) + 1) * (3 * n + 2) * (3 * n + 3)) := by
      rw [eq_div_iff (by positivity)]; exact choose_six_step n
    rw [hc3, hc6, show 15 * (n + 1) = 15 * n + 15 by ring, pow_add]
    field_simp
    ring

/-- the WZ sum: `∑_n G(n,0) = √2/(2π)` -/
theorem hasSum_G0 : HasSum (fun n => Q n 0 * hterm σ ρ n 0) (Real.sqrt 2 / (2 * Real.pi)) :=
  hasSum_of_ratios σ ρ Q (fun n => 3 * (n : ℝ) + 3) 3 (27 / 32) _ compat wz_rat σ_pos
    (fun n => ρ_pos n 0) (by norm_num) σ_le (by norm_num) ρ0_le Q_abs_le (fun n => Q_nonneg n 0)
    tendsto_boundary

/-- **Sun's conjecture 5.25(i) of arXiv:2307.03086:**
`∑_{k≥0} (360k³+612k²+230k+15) C(3k,k) C(6k,3k)² / ((3k+1)(3k+2) 2^{15k}) = 32√2/π`. -/
theorem sun_series_32_sqrt2_div_pi :
    HasSum (fun k : ℕ => (360 * (k : ℝ) ^ 3 + 612 * (k : ℝ) ^ 2 + 230 * k + 15)
        * (Nat.choose (3 * k) k : ℝ) * (Nat.choose (6 * k) (3 * k) : ℝ) ^ 2
        / ((3 * (k : ℝ) + 1) * (3 * k + 2) * 2 ^ (15 * k)))
      (32 * Real.sqrt 2 / Real.pi) := by
  have h := hasSum_G0.mul_left 64
  have e : (64 : ℝ) * (Real.sqrt 2 / (2 * Real.pi)) = 32 * Real.sqrt 2 / Real.pi := by
    have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
    field_simp; ring
  rw [e] at h
  refine h.congr_fun ?_
  intro n
  rw [hterm_zero_right_eq]
  unfold Q
  have h4 : (3 * (n : ℝ) + 1) ≠ 0 := by positivity
  have h5 : (3 * (n : ℝ) + 2) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

end PiSeries.SunSqrt2
