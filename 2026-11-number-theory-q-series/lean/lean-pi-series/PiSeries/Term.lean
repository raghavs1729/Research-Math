/-
# Bivariate hypergeometric terms given by their two ratios, and the packaged WZ theorem

`hterm σ ρ n k` is the term `h(n,k)` with `h(0,0) = 1`, `h(n+1,0)/h(n,0) = ρ(n,0)` and
`h(n,k+1)/h(n,k) = σ(n,k)`.  When the two ratio functions are compatible
(`ρ(n,k) σ(n+1,k) = σ(n,k) ρ(n,k+1)`), also `h(n+1,k)/h(n,k) = ρ(n,k)` for every `k`.

`hasSum_of_ratios` reduces a series `∑_n Q(n,0) h(n,0) = L` to
* two identities between rational functions (compatibility and the WZ equation),
* three inequalities between rational functions,
* one limit of a single hypergeometric term, `Q(0,K) h(0,K) → L`.
-/
import PiSeries.WZ

open Finset Filter Topology

namespace PiSeries

/-- The bivariate hypergeometric term with `h(0,0) = 1`, `h(n+1,0) = ρ(n,0) h(n,0)`,
`h(n,k+1) = σ(n,k) h(n,k)`. -/
noncomputable def hterm (σ ρ : ℕ → ℕ → ℝ) (n k : ℕ) : ℝ :=
  (∏ j ∈ range n, ρ j 0) * ∏ j ∈ range k, σ n j

section
variable {σ ρ : ℕ → ℕ → ℝ}

lemma hterm_zero_right (n : ℕ) : hterm σ ρ n 0 = ∏ j ∈ range n, ρ j 0 := by simp [hterm]

lemma hterm_zero_left (k : ℕ) : hterm σ ρ 0 k = ∏ j ∈ range k, σ 0 j := by simp [hterm]

lemma hterm_succ_right (n k : ℕ) : hterm σ ρ n (k + 1) = hterm σ ρ n k * σ n k := by
  simp [hterm, prod_range_succ, mul_assoc]

lemma hterm_succ_zero (n : ℕ) : hterm σ ρ (n + 1) 0 = hterm σ ρ n 0 * ρ n 0 := by
  simp [hterm, prod_range_succ]

lemma hterm_succ_left (compat : ∀ n k, ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1)) (n k : ℕ) :
    hterm σ ρ (n + 1) k = hterm σ ρ n k * ρ n k := by
  induction k with
  | zero => exact hterm_succ_zero n
  | succ k ih => rw [hterm_succ_right, ih, hterm_succ_right, mul_assoc, compat, mul_assoc]

lemma hterm_pos (hσ : ∀ n k, 0 < σ n k) (hρ : ∀ n, 0 < ρ n 0) (n k : ℕ) : 0 < hterm σ ρ n k :=
  mul_pos (prod_pos fun j _ => hρ j) (prod_pos fun j _ => hσ n j)

lemma hterm_zero_le_pow {r : ℝ} (hρ : ∀ n, 0 < ρ n 0) (hr : ∀ n, ρ n 0 ≤ r) (n : ℕ) :
    hterm σ ρ n 0 ≤ r ^ n := by
  rw [hterm_zero_right]
  calc ∏ j ∈ range n, ρ j 0 ≤ ∏ _j ∈ range n, r :=
        prod_le_prod (fun j _ => (hρ j).le) (fun j _ => hr j)
    _ = r ^ n := by simp

/-- If `σ(n,k) ≤ ((k+a)/(k+a+1))²` for all `k` (with `a ≥ 1`), then
`h(n,k) ≤ h(n,0) a²/(k+1)²`. -/
lemma hterm_le_of_sigma_le {a : ℝ} (ha : 1 ≤ a) (n : ℕ) (hσ0 : ∀ k, 0 < σ n k)
    (hc : 0 ≤ hterm σ ρ n 0)
    (hσ : ∀ k : ℕ, σ n k ≤ (((k : ℝ) + a) / ((k : ℝ) + a + 1)) ^ 2) (k : ℕ) :
    hterm σ ρ n k ≤ hterm σ ρ n 0 * a ^ 2 / ((k : ℝ) + 1) ^ 2 := by
  have ha0 : 0 < a := by linarith
  have key : hterm σ ρ n k ≤ hterm σ ρ n 0 * (a / ((k : ℝ) + a)) ^ 2 := by
    induction k with
    | zero => simp [div_self ha0.ne']
    | succ k ih =>
      rw [hterm_succ_right]
      have hk : (0 : ℝ) < (k : ℝ) + a := by positivity
      have hk1 : (0 : ℝ) < (k : ℝ) + a + 1 := by positivity
      calc hterm σ ρ n k * σ n k
          ≤ (hterm σ ρ n 0 * (a / ((k : ℝ) + a)) ^ 2) * (((k : ℝ) + a) / ((k : ℝ) + a + 1)) ^ 2 :=
            mul_le_mul ih (hσ k) (hσ0 k).le (by positivity)
        _ = hterm σ ρ n 0 * (a / (((k + 1 : ℕ) : ℝ) + a)) ^ 2 := by
            push_cast
            field_simp
            ring
  have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have hka : (0 : ℝ) < (k : ℝ) + a := by positivity
  have h1 : a / ((k : ℝ) + a) ≤ a / ((k : ℝ) + 1) :=
    div_le_div_of_nonneg_left ha0.le hk (by linarith)
  have h2 : (a / ((k : ℝ) + a)) ^ 2 ≤ (a / ((k : ℝ) + 1)) ^ 2 :=
    pow_le_pow_left₀ (by positivity) h1 2
  calc hterm σ ρ n k ≤ hterm σ ρ n 0 * (a / ((k : ℝ) + a)) ^ 2 := key
    _ ≤ hterm σ ρ n 0 * (a / ((k : ℝ) + 1)) ^ 2 := mul_le_mul_of_nonneg_left h2 hc
    _ = hterm σ ρ n 0 * a ^ 2 / ((k : ℝ) + 1) ^ 2 := by rw [div_pow]; ring

end

/-- **The WZ theorem for a pair given by its ratios.**  Let `σ, ρ, Q` be functions `ℕ → ℕ → ℝ`
(in the applications: explicit rational functions) such that
* `ρ(n,k) σ(n+1,k) = σ(n,k) ρ(n,k+1)`   (compatibility),
* `(n+1) ρ(n,k) − n = Q(n,k+1) σ(n,k) − Q(n,k)`   (the WZ equation divided by `h(n,k)`),
* `σ > 0`, `0 < ρ(n,0) ≤ r < 1`,
* `σ(n,k) ≤ ((k+cn)/(k+cn+1))²` for `n ≥ 1` (some `c ≥ 1`),
* `|Q(n,k)| ≤ A(n)(k+1)` for `n ≥ 1`, `Q(n,0) ≥ 0`,
* `Q(0,K) h(0,K) → L`, where `h = hterm σ ρ`.
Then `∑_n Q(n,0) h(n,0) = L`. -/
theorem hasSum_of_ratios (σ ρ Q : ℕ → ℕ → ℝ) (A : ℕ → ℝ) (c r L : ℝ)
    (compat : ∀ n k, ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1))
    (hrat : ∀ n k : ℕ, ((n : ℝ) + 1) * ρ n k - n = Q n (k + 1) * σ n k - Q n k)
    (hσpos : ∀ n k, 0 < σ n k) (hρpos : ∀ n, 0 < ρ n 0)
    (hc : 1 ≤ c)
    (hσ : ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ,
      σ n k ≤ (((k : ℝ) + c * n) / ((k : ℝ) + c * n + 1)) ^ 2)
    (hr : r < 1) (hρ : ∀ n, ρ n 0 ≤ r)
    (hA : ∀ n : ℕ, 1 ≤ n → ∀ k : ℕ, |Q n k| ≤ A n * ((k : ℝ) + 1))
    (hQ0 : ∀ n, 0 ≤ Q n 0)
    (hL : Tendsto (fun K => Q 0 K * hterm σ ρ 0 K) atTop (𝓝 L)) :
    HasSum (fun n => Q n 0 * hterm σ ρ n 0) L := by
  have hr0 : 0 < r := (hρpos 0).trans_le (hρ 0)
  have hpos := hterm_pos (ρ := ρ) hσpos hρpos
  refine hasSum_of_wz (hterm σ ρ) Q (fun n => r ^ n * (c * n) ^ 2) A L ?_ (fun n k => (hpos n k).le)
    ?_ hA ?_ hL (fun n => mul_nonneg (hQ0 n) (hpos n 0).le)
  · intro n k
    rw [hterm_succ_left compat, hterm_succ_right]
    linear_combination (hterm σ ρ n k) * hrat n k
  · intro n hn k
    have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
    have hcn : 1 ≤ c * n := by nlinarith
    have h1 := hterm_le_of_sigma_le (ρ := ρ) hcn n (hσpos n) (hpos n 0).le (hσ n hn) k
    have h2 : hterm σ ρ n 0 * (c * n) ^ 2 ≤ r ^ n * (c * n) ^ 2 :=
      mul_le_mul_of_nonneg_right (hterm_zero_le_pow hρpos hρ n) (by positivity)
    exact h1.trans (div_le_div_of_nonneg_right h2 (by positivity))
  · have h3 := (tendsto_pow_const_mul_const_pow_of_abs_lt_one 3
      (by rw [abs_of_pos hr0]; exact hr : |r| < 1)).const_mul (c ^ 2)
    rw [mul_zero] at h3
    refine h3.congr (fun n => ?_)
    ring

end PiSeries
