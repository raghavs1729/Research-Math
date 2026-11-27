/-
Tau/Main.lean — assembly of Theorem 1 of PROOF.md and its corollaries.
-/
import Tau.Defs
import Tau.Descent
import Tau.Sieve

namespace Tau

lemma cert_461 : certOK seeds461 [3, 7, 43, 89, 199, 967] 88 = true := by
  decide +kernel

lemma cert_M599 : certOK seedsM599 [89] 22 = true := by
  decide +kernel

lemma form_eq (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 = Qf (x ^ 11, y ^ 2) := by
  simp only [Qf]; ring

theorem no_solution_461' (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ 461 := by
  intro h
  rw [form_eq] at h
  obtain ⟨s, hs, n, hn⟩ := orbit_461 _ h
  exact sieve (by norm_num) cert_461 hs n x y hn

theorem no_solution_M599' (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ -599 := by
  intro h
  rw [form_eq] at h
  obtain ⟨s, hs, n, hn⟩ := orbit_M599 _ h
  exact sieve (by norm_num) cert_M599 hs n x y hn

lemma heckeSeq_four (t p : ℤ) : heckeSeq t p 4 = p ^ 22 - 3 * p ^ 11 * t ^ 2 + t ^ 4 := by
  simp only [heckeSeq]; ring

theorem heckeSeq_four_ne' (t p : ℤ) : heckeSeq t p 4 ≠ 461 ∧ heckeSeq t p 4 ≠ -599 := by
  rw [heckeSeq_four]
  exact ⟨no_solution_461' p t, no_solution_M599' p t⟩

/-- Any `τ : ℕ → ℤ` with `τ 1 = 1` and the Hecke recurrence at the powers of a prime `p` satisfies
`τ (p ^ k) = heckeSeq (τ p) p k`. -/
lemma eq_heckeSeq (τ : ℕ → ℤ) (p : ℕ) (h1 : τ 1 = 1)
    (hrec : ∀ k : ℕ, τ (p ^ (k + 2)) = τ p * τ (p ^ (k + 1)) - (p : ℤ) ^ 11 * τ (p ^ k)) :
    ∀ k, τ (p ^ k) = heckeSeq (τ p) p k
  | 0 => by simp [heckeSeq, h1]
  | 1 => by simp [heckeSeq]
  | k + 2 => by
    rw [hrec k, eq_heckeSeq τ p h1 hrec (k + 1), eq_heckeSeq τ p h1 hrec k]
    rfl

end Tau
