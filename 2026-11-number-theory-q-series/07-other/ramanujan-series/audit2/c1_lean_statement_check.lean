import RamanujanSeries
open SunSeries16
-- the headline statement, restated from scratch (no project definitions) and closed by the project's theorem
example : HasSum (fun k : ℕ => (3 * (k : ℝ) + 2) ^ 2 * 16 ^ k /
      ((2 * k + 1) * (6 * k + 1) * (6 * k + 5) * ((Nat.choose (6 * k) (3 * k) : ℕ) : ℝ)))
      (Real.pi / (2 * Real.sqrt 3)) := by
  have h := sun_series_16
  have e : ∀ k : ℕ, Nat.centralBinom (3 * k) = Nat.choose (6 * k) (3 * k) := by
    intro k; rw [Nat.centralBinom_eq_two_mul_choose]; congr 1; ring
  simpa [e] using h
-- summability and the tsum form follow from HasSum
example : Summable (fun k : ℕ => (3 * (k : ℝ) + 2) ^ 2 * 16 ^ k /
      ((2 * k + 1) * (6 * k + 1) * (6 * k + 5) * (Nat.centralBinom (3 * k) : ℝ))) := sun_series_16.summable
example : (∑' k : ℕ, (3 * (k : ℝ) + 2) ^ 2 * 16 ^ k /
      ((2 * k + 1) * (6 * k + 1) * (6 * k + 5) * (Nat.centralBinom (3 * k) : ℝ))) = Real.pi / (2 * Real.sqrt 3) :=
  sun_series_16.tsum_eq
#print axioms sun_series_16
