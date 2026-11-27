/-
Tau/Descent.lean — Lemma A (seed bound) of PROOF.md: every integer solution of `Q (U, V) = α` lies in
the `T`-orbit of a solution `(z, b)` with `|z|` minimal along the orbit, and for such a pair
`2 z² ≤ |z² − α|`.  For `α = 461` this gives `|z| ≤ 12`, for `α = −599` it gives `|z| ≤ 24`; the seeds
are then found by a finite search (`decide +kernel`).
-/
import Tau.Basic

namespace Tau

/-- Every point of `ℤ × ℤ` has, in its `T`-orbit, a point whose first coordinate has minimal absolute
value among all points of the orbit. -/
lemma exists_min_orbit (p : ℤ × ℤ) :
    ∃ m : ℤ, ∀ n : ℤ, ((T ℤ ^ m) p).1.natAbs ≤ ((T ℤ ^ n) p).1.natAbs := by
  classical
  have hex : ∃ k : ℕ, ∃ n : ℤ, ((T ℤ ^ n) p).1.natAbs = k := ⟨_, 0, rfl⟩
  obtain ⟨m, hm⟩ := Nat.find_spec hex
  refine ⟨m, fun n => ?_⟩
  rw [hm]
  exact Nat.find_min' hex ⟨n, rfl⟩

/-- The key inequality at an orbit minimum: if `(z, b)` solves `Q = α`, `a = 3z − b`, and
`|z| ≤ |a|`, `|z| ≤ |b|`, then `2 z² ≤ |z² − α|`. -/
lemma key_ineq (z a b α : ℤ) (hQ : z ^ 2 - 3 * z * b + b ^ 2 = α) (hab : a + b = 3 * z)
    (ha : |z| ≤ |a|) (hb : |z| ≤ |b|) : 2 * z ^ 2 ≤ |z ^ 2 - α| := by
  have hprod : a * b = z ^ 2 - α := by
    have : a = 3 * z - b := by linarith
    subst this; linarith
  have h1 : 0 ≤ (|a| - |z|) * (|b| - |z|) := mul_nonneg (by linarith) (by linarith)
  have h2 : |a + b| ≤ |a| + |b| := abs_add_le a b
  have h3 : |a + b| = 3 * |z| := by rw [hab, abs_mul]; norm_num
  have h4 : |a| * |b| = |z ^ 2 - α| := by rw [← abs_mul, hprod]
  have h5 : |z| ^ 2 = z ^ 2 := sq_abs z
  have hz : 0 ≤ |z| := abs_nonneg z
  nlinarith

/-- Generic descent: every solution of `Q = α` is `T ^ n` of a solution `(z, b)` with
`2 z² ≤ |z² − α|`. -/
lemma descent (p : ℤ × ℤ) :
    ∃ s : ℤ × ℤ, Qf s = Qf p ∧ 2 * s.1 ^ 2 ≤ |s.1 ^ 2 - Qf p| ∧ ∃ n : ℤ, (T ℤ ^ n) s = p := by
  obtain ⟨m, hm⟩ := exists_min_orbit p
  refine ⟨(T ℤ ^ m) p, Qf_zpow m p, ?_, -m, ?_⟩
  · set z := ((T ℤ ^ m) p).1
    set b := ((T ℤ ^ m) p).2
    set a := ((T ℤ ^ (m - 1)) p).1
    have hb : b = ((T ℤ ^ (m + 1)) p).1 := (fst_zpow_succ m p).symm
    have hab : a + b = 3 * z := by
      have := snd_zpow m p
      linarith
    have hQ : z ^ 2 - 3 * z * b + b ^ 2 = Qf p := Qf_zpow m p
    apply key_ineq z a b (Qf p) hQ hab
    · have := hm (m - 1)
      rw [Int.abs_eq_natAbs, Int.abs_eq_natAbs]; exact_mod_cast this
    · have := hm (m + 1)
      rw [Int.abs_eq_natAbs, Int.abs_eq_natAbs, hb]; exact_mod_cast this
  · rw [← Equiv.Perm.mul_apply, ← zpow_add, neg_add_cancel, zpow_zero, Equiv.Perm.one_apply]

/-! ### Finite seed search -/

/-- All pairs `(z, b)` with `|z| ≤ B`, `|b| ≤ C` and `z² − 3zb + b² = α` belong to the list `S`. -/
def seedsOK (α : ℤ) (B C : ℕ) (S : List (ℤ × ℤ)) : Bool :=
  (List.range (2 * B + 1)).all fun i => (List.range (2 * C + 1)).all fun j =>
    let z : ℤ := (i : ℤ) - B
    let b : ℤ := (j : ℤ) - C
    !(z * z - 3 * z * b + b * b == α) || S.contains (z, b)

lemma seedsOK_sound {α : ℤ} {B C : ℕ} {S : List (ℤ × ℤ)} (h : seedsOK α B C S = true)
    (z b : ℤ) (hz : |z| ≤ B) (hb : |b| ≤ C) (hQ : z ^ 2 - 3 * z * b + b ^ 2 = α) : (z, b) ∈ S := by
  rw [abs_le] at hz hb
  simp only [seedsOK, List.all_eq_true, List.mem_range] at h
  have hi : (z + B).toNat < 2 * B + 1 := by omega
  have hj : (b + C).toNat < 2 * C + 1 := by omega
  have := h _ hi _ hj
  have ez : (((z + B).toNat : ℕ) : ℤ) - B = z := by omega
  have eb : (((b + C).toNat : ℕ) : ℤ) - C = b := by omega
  rw [ez, eb] at this
  have hQ' : z * z - 3 * z * b + b * b = α := by rw [← hQ]; ring
  simpa [hQ'] using this

/-- The seeds for `α = 461`. -/
def seeds461 : List (ℤ × ℤ) := [(-1, -23), (-1, 20), (1, -20), (1, 23)]

/-- The seeds for `α = −599`. -/
def seedsM599 : List (ℤ × ℤ) := [(-24, -47), (-24, -25), (24, 25), (24, 47)]

lemma seedsOK_461 : seedsOK 461 12 50 seeds461 = true := by decide +kernel

lemma seedsOK_M599 : seedsOK (-599) 24 50 seedsM599 = true := by decide +kernel

/-- Every integer solution of `U² − 3UV + V² = 461` is `T ^ n` of one of the four seeds. -/
theorem orbit_461 (p : ℤ × ℤ) (hp : Qf p = 461) : ∃ s ∈ seeds461, ∃ n : ℤ, (T ℤ ^ n) s = p := by
  obtain ⟨s, hs, hineq, n, hn⟩ := descent p
  refine ⟨s, ?_, n, hn⟩
  rw [hp] at hs hineq
  obtain ⟨z, b⟩ := s
  simp only [Qf] at hs hineq
  have hz2 : 3 * z ^ 2 ≤ 461 := by
    rcases abs_cases (z ^ 2 - 461) with ⟨h, _⟩ | ⟨h, _⟩ <;> nlinarith [sq_nonneg z]
  have hz' : -12 ≤ z ∧ z ≤ 12 := by constructor <;> nlinarith
  have hsq : (2 * b - 3 * z) ^ 2 = 5 * z ^ 2 + 4 * (461) := by linear_combination 4 * hs
  have hb' : -50 ≤ b ∧ b ≤ 50 := by
    constructor
    · by_contra hc; replace hc := lt_of_not_ge hc
      have : 3 * z - 2 * b ≥ 64 := by omega
      nlinarith
    · by_contra hc; replace hc := lt_of_not_ge hc
      have : 2 * b - 3 * z ≥ 64 := by omega
      nlinarith
  have hz : |z| ≤ ((12 : ℕ) : ℤ) := by rw [abs_le]; push_cast; omega
  have hb : |b| ≤ ((50 : ℕ) : ℤ) := by rw [abs_le]; push_cast; omega
  exact seedsOK_sound seedsOK_461 z b hz hb hs

/-- Every integer solution of `U² − 3UV + V² = −599` is `T ^ n` of one of the four seeds. -/
theorem orbit_M599 (p : ℤ × ℤ) (hp : Qf p = -599) :
    ∃ s ∈ seedsM599, ∃ n : ℤ, (T ℤ ^ n) s = p := by
  obtain ⟨s, hs, hineq, n, hn⟩ := descent p
  refine ⟨s, ?_, n, hn⟩
  rw [hp] at hs hineq
  obtain ⟨z, b⟩ := s
  simp only [Qf] at hs hineq
  have hz2 : z ^ 2 ≤ 599 := by
    rcases abs_cases (z ^ 2 - -599) with ⟨h, _⟩ | ⟨h, _⟩ <;> nlinarith [sq_nonneg z]
  have hz' : -24 ≤ z ∧ z ≤ 24 := by constructor <;> nlinarith
  have hsq : (2 * b - 3 * z) ^ 2 = 5 * z ^ 2 + 4 * (-599) := by linear_combination 4 * hs
  have hb' : -50 ≤ b ∧ b ≤ 50 := by
    constructor
    · by_contra hc; replace hc := lt_of_not_ge hc
      have : 3 * z - 2 * b ≥ 30 := by omega
      nlinarith
    · by_contra hc; replace hc := lt_of_not_ge hc
      have : 2 * b - 3 * z ≥ 30 := by omega
      nlinarith
  have hz : |z| ≤ ((24 : ℕ) : ℤ) := by rw [abs_le]; push_cast; omega
  have hb : |b| ≤ ((50 : ℕ) : ℤ) := by rw [abs_le]; push_cast; omega
  exact seedsOK_sound seedsOK_M599 z b hz hb hs

end Tau
