/-
# The Wallis limit `(2K+1) (C(2K,K)/4^K)² → 2/π`, and the recurrence of `C(2n,n)`

`wallisTerm K = ∏_{j<K} (2j+1)²/(4(j+1)²) = ((1/2)_K / K!)² = (C(2K,K)/4^K)²` is the boundary row
`h(0,K)` of the WZ pairs for Ramanujan's two series.  Its asymptotic is Wallis' product
(Mathlib: `Real.Wallis.tendsto_W_nhds_pi_div_two`).
-/
import Mathlib

open Finset Filter Topology

namespace PiSeries

/-- `∏_{j<K} (2j+1)²/(4(j+1)²)`. -/
noncomputable def wallisTerm (K : ℕ) : ℝ :=
  ∏ j ∈ range K, (2 * (j : ℝ) + 1) ^ 2 / (4 * ((j : ℝ) + 1) ^ 2)

lemma wallisTerm_succ (K : ℕ) :
    wallisTerm (K + 1) = wallisTerm K * ((2 * (K : ℝ) + 1) ^ 2 / (4 * ((K : ℝ) + 1) ^ 2)) :=
  prod_range_succ _ _

/-- `(2K+1) · wallisTerm K · W(K) = 1`, where `W(K)` is the Wallis partial product. -/
lemma wallisTerm_mul_W (K : ℕ) : (2 * (K : ℝ) + 1) * wallisTerm K * Real.Wallis.W K = 1 := by
  induction K with
  | zero => simp [wallisTerm, Real.Wallis.W]
  | succ K ih =>
    rw [wallisTerm_succ, Real.Wallis.W_succ]
    push_cast
    have h1 : (2 * (K : ℝ) + 1) ≠ 0 := by positivity
    have h2 : ((K : ℝ) + 1) ≠ 0 := by positivity
    have h3 : (2 * (K : ℝ) + 3) ≠ 0 := by positivity
    have e : (2 * ((K : ℝ) + 1) + 1) * (wallisTerm K * ((2 * (K : ℝ) + 1) ^ 2 / (4 * ((K : ℝ) + 1) ^ 2)))
        * (Real.Wallis.W K * ((2 * (K : ℝ) + 2) / (2 * (K : ℝ) + 1) * ((2 * (K : ℝ) + 2) / (2 * (K : ℝ) + 3))))
        = (2 * (K : ℝ) + 1) * wallisTerm K * Real.Wallis.W K := by
      field_simp
      ring
    rw [e, ih]

/-- **Wallis:** `(2K+1) ∏_{j<K} (2j+1)²/(4(j+1)²) → 2/π`. -/
theorem tendsto_wallisTerm :
    Tendsto (fun K : ℕ => (2 * (K : ℝ) + 1) * wallisTerm K) atTop (𝓝 (2 / Real.pi)) := by
  have h := (Real.Wallis.tendsto_W_nhds_pi_div_two).inv₀ (by positivity : Real.pi / 2 ≠ 0)
  rw [inv_div] at h
  refine h.congr (fun K => ?_)
  have hW := Real.Wallis.W_pos K
  rw [inv_eq_one_div, eq_comm, eq_div_iff hW.ne']
  exact wallisTerm_mul_W K

/-- If `q(K)/(2K+1) → ℓ` then `q(K) · wallisTerm K → ℓ · 2/π`. -/
theorem tendsto_mul_wallisTerm {q : ℕ → ℝ} {ℓ : ℝ}
    (hq : Tendsto (fun K : ℕ => q K / (2 * (K : ℝ) + 1)) atTop (𝓝 ℓ)) :
    Tendsto (fun K : ℕ => q K * wallisTerm K) atTop (𝓝 (ℓ * (2 / Real.pi))) := by
  refine (hq.mul tendsto_wallisTerm).congr (fun K => ?_)
  have h1 : (2 * (K : ℝ) + 1) ≠ 0 := by positivity
  field_simp

/-- `(n+1) C(2n+2,n+1) = 2(2n+1) C(2n,n)`, over `ℝ`. -/
lemma choose_two_mul_succ (n : ℕ) :
    ((n : ℝ) + 1) * (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)
      = 2 * (2 * (n : ℝ) + 1) * (Nat.choose (2 * n) n : ℝ) := by
  have h := Nat.succ_mul_centralBinom_succ n
  rw [Nat.centralBinom_eq_two_mul_choose, Nat.centralBinom_eq_two_mul_choose] at h
  exact_mod_cast h

end PiSeries
