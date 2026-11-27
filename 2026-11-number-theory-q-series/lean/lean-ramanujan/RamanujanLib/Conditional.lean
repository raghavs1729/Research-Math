/-
RamanujanLib/Conditional.lean — consequences of the explicit criterion (Mathlib only).

* `ramanujan_of_schoenfeld`: a bound of Schoenfeld's shape
  `|π(x) − (c + ∫_2^x dt/log t)| < √x log x / (8π)` for `x ≥ 2657` gives Ramanujan's inequality
  for every real `x ≥ e⁴⁰`.
* `ramanujan_real_of_nat`: the inequality at `⌊x⌋` implies it at `x` (`x ≥ 3`).
* `ramanujan_sharp_reduction`: that bound plus the inequality at the integers of
  `[38 358 837 683, e⁴⁰)` give it for every real `x ≥ 38 358 837 683`.
-/
import RamanujanLib.Explicit

open Real

noncomputable section

namespace RamanujanLib

theorem ramanujan_of_schoenfeld (c : ℝ) (hc : |c| ≤ 2)
    (hS : ∀ x : ℝ, 2657 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..x, 1 / log t)|
        < √x * log x / (8 * π))
    (x : ℝ) (hx : exp 40 ≤ x) :
    (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
  have hpi : (3.14 : ℝ) < π := Real.pi_gt_d2
  refine ramanujan_of_sqrt_log_bound (1 / (8 * π)) c (by positivity) ?_ hc ?_ x hx
  · rw [div_le_div_iff₀ (by positivity) (by norm_num)]
    linarith
  · intro y hy
    have h2657 : (2657 : ℝ) ≤ y := by
      have h1 : (2657 : ℝ) ≤ exp 20 := le_trans (by norm_num) exp_twenty_ge
      have h2 : exp 20 ≤ exp 39 := exp_le_exp.mpr (by norm_num)
      linarith
    calc |(Nat.primeCounting ⌊y⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..y, 1 / log t)|
        ≤ √y * log y / (8 * π) := (hS y h2657).le
      _ = 1 / (8 * π) * (√y * log y) := by ring

/-- The inequality at the integer `⌊x⌋` implies the inequality at the real `x ≥ 3`: the left side
does not change, and `t / log t` and `π(t / e)` are nondecreasing for `t ≥ e`. -/
theorem ramanujan_real_of_nat (x : ℝ) (hx : 3 ≤ x)
    (hn : (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * (⌊x⌋₊ : ℝ) / log (⌊x⌋₊ : ℝ)
        * (Nat.primeCounting ⌊(⌊x⌋₊ : ℝ) / exp 1⌋₊ : ℝ)) :
    (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
  have hE0 : 0 < exp 1 := exp_pos 1
  have hE3 : exp 1 ≤ 3 := by
    have := Real.exp_one_lt_d9
    linarith
  have hx0 : 0 ≤ x := by linarith
  have hnx : (⌊x⌋₊ : ℝ) ≤ x := Nat.floor_le hx0
  have h3n : (3 : ℝ) ≤ (⌊x⌋₊ : ℝ) := by
    have : (3 : ℕ) ≤ ⌊x⌋₊ := Nat.le_floor (by exact_mod_cast hx)
    exact_mod_cast this
  have hlogn : 0 < log (⌊x⌋₊ : ℝ) := log_pos (by linarith)
  have hlogx : 0 < log x := log_pos (by linarith)
  -- `n / log n ≤ x / log x`
  have hanti := Real.log_div_self_antitoneOn
    (show (⌊x⌋₊ : ℝ) ∈ {t : ℝ | exp 1 ≤ t} by simp only [Set.mem_ofPred_eq]; linarith)
    (show x ∈ {t : ℝ | exp 1 ≤ t} by simp only [Set.mem_ofPred_eq]; linarith) hnx
  have hdiv : (⌊x⌋₊ : ℝ) / log (⌊x⌋₊ : ℝ) ≤ x / log x := by
    have hpos : 0 < log x / x := div_pos hlogx (by linarith)
    have := inv_anti₀ hpos hanti
    rwa [inv_div, inv_div] at this
  -- `π(n / e) ≤ π(x / e)`
  have hπ : (Nat.primeCounting ⌊(⌊x⌋₊ : ℝ) / exp 1⌋₊ : ℝ)
      ≤ (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
    have : ⌊(⌊x⌋₊ : ℝ) / exp 1⌋₊ ≤ ⌊x / exp 1⌋₊ :=
      Nat.floor_le_floor (div_le_div_of_nonneg_right hnx hE0.le)
    exact_mod_cast Nat.monotone_primeCounting this
  have hcoef : exp 1 * (⌊x⌋₊ : ℝ) / log (⌊x⌋₊ : ℝ) ≤ exp 1 * x / log x := by
    rw [mul_div_assoc, mul_div_assoc]
    exact mul_le_mul_of_nonneg_left hdiv hE0.le
  calc (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * (⌊x⌋₊ : ℝ) / log (⌊x⌋₊ : ℝ)
        * (Nat.primeCounting ⌊(⌊x⌋₊ : ℝ) / exp 1⌋₊ : ℝ) := hn
    _ ≤ exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) :=
        mul_le_mul hcoef hπ (Nat.cast_nonneg _) (by positivity)

/-- Reduction of the sharp threshold `38 358 837 683` to a Schoenfeld-type bound and a statement
about the integers below `e⁴⁰`. -/
theorem ramanujan_sharp_reduction (c : ℝ) (hc : |c| ≤ 2)
    (hS : ∀ x : ℝ, 2657 ≤ x →
      |(Nat.primeCounting ⌊x⌋₊ : ℝ) - (c + ∫ t in (2 : ℝ)..x, 1 / log t)|
        < √x * log x / (8 * π))
    (hFinite : ∀ n : ℕ, 38358837683 ≤ n → (n : ℝ) < exp 40 →
      (Nat.primeCounting n : ℝ) ^ 2
        < exp 1 * (n : ℝ) / log (n : ℝ) * (Nat.primeCounting ⌊(n : ℝ) / exp 1⌋₊ : ℝ))
    (x : ℝ) (hx : 38358837683 ≤ x) :
    (Nat.primeCounting ⌊x⌋₊ : ℝ) ^ 2
      < exp 1 * x / log x * (Nat.primeCounting ⌊x / exp 1⌋₊ : ℝ) := by
  rcases le_or_gt (exp 40) x with h | h
  · exact ramanujan_of_schoenfeld c hc hS x h
  · refine ramanujan_real_of_nat x (by linarith) ?_
    have h1 : 38358837683 ≤ ⌊x⌋₊ := Nat.le_floor (by exact_mod_cast hx)
    have h2 : (⌊x⌋₊ : ℝ) < exp 40 := lt_of_le_of_lt (Nat.floor_le (by linarith)) h
    exact hFinite ⌊x⌋₊ h1 h2

end RamanujanLib
