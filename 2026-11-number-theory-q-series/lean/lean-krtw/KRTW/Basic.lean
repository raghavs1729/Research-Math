/-
KRTW/Basic.lean — q-Pochhammer symbols, their inverses and Gaussian binomials as formal power series.

Everything lives in `K⟦X⟧` for a field `K`.  The "base" `q` is an arbitrary power series with zero
constant term (class `Nil q`); in the application `q = X` or `q = X ^ 2`.
-/
import Mathlib

namespace KRTW

open PowerSeries Finset

/-! ### Triangular numbers -/

/-- The triangular number `n (n + 1) / 2`. -/
def tri : ℕ → ℕ
  | 0 => 0
  | n + 1 => tri n + (n + 1)

@[simp] lemma tri_zero : tri 0 = 0 := rfl

lemma tri_succ (n : ℕ) : tri (n + 1) = tri n + (n + 1) := rfl

@[simp] lemma tri_one : tri 1 = 1 := rfl

lemma two_mul_tri (n : ℕ) : 2 * tri n = n * (n + 1) := by
  induction n with
  | zero => rfl
  | succ n ih => rw [tri_succ, mul_add, ih]; ring

lemma tri_eq (n : ℕ) : tri n = n * (n + 1) / 2 := by
  rw [← two_mul_tri]; simp

lemma tri_add (a b : ℕ) : tri (a + b) = tri a + tri b + a * b := by
  have h1 := two_mul_tri (a + b)
  have h2 := two_mul_tri a
  have h3 := two_mul_tri b
  nlinarith

lemma tri_two_mul (b : ℕ) : tri (2 * b) = 2 * b ^ 2 + b := by
  have h1 := two_mul_tri (2 * b)
  nlinarith

lemma le_tri (n : ℕ) : n ≤ tri n := by
  cases n with
  | zero => simp
  | succ n => rw [tri_succ]; omega

variable {K : Type*} [Field K]

/-- The base `q` has zero constant term. -/
class Nil (q : K⟦X⟧) : Prop where
  cc : constantCoeff q = 0

instance : Nil (X : K⟦X⟧) := ⟨by simp⟩

instance (n : ℕ) : Nil ((X : K⟦X⟧) ^ (n + 1)) := ⟨by simp⟩

lemma X_dvd_of_nil (q : K⟦X⟧) [Nil q] : (X : K⟦X⟧) ∣ q :=
  X_dvd_iff.mpr Nil.cc

lemma X_pow_dvd_pow (q : K⟦X⟧) [Nil q] (n : ℕ) : (X : K⟦X⟧) ^ n ∣ q ^ n :=
  pow_dvd_pow_of_dvd (X_dvd_of_nil q) n

lemma X_pow_dvd_pow_of_le (q : K⟦X⟧) [Nil q] {m n : ℕ} (h : m ≤ n) : (X : K⟦X⟧) ^ m ∣ q ^ n :=
  (pow_dvd_pow X h).trans (X_pow_dvd_pow q n)

section base

variable (q : K⟦X⟧)

/-- `(q;q)_n = ∏_{i=1}^n (1 - q^i)`. -/
noncomputable def poch (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - q ^ (i + 1))

/-- `1/(q;q)_n`. -/
noncomputable def iq (n : ℕ) : K⟦X⟧ := (poch q n)⁻¹

@[simp] lemma poch_zero : poch q 0 = 1 := by simp [poch]

lemma poch_succ (n : ℕ) : poch q (n + 1) = poch q n * (1 - q ^ (n + 1)) := by
  simp [poch, prod_range_succ]

@[simp] lemma iq_zero : iq q 0 = 1 := by simp [iq]

variable [Nil q]

lemma constantCoeff_one_sub_pow (n : ℕ) : constantCoeff (1 - q ^ (n + 1)) = 1 := by
  simp [Nil.cc]

lemma one_sub_pow_ne_zero (n : ℕ) : (1 - q ^ (n + 1) : K⟦X⟧) ≠ 0 := by
  intro h
  have := constantCoeff_one_sub_pow q n
  rw [h] at this
  simp at this

lemma constantCoeff_poch (n : ℕ) : constantCoeff (poch q n) = 1 := by
  induction n with
  | zero => simp
  | succ n ih => rw [poch_succ, map_mul, ih, constantCoeff_one_sub_pow, mul_one]

lemma poch_mul_iq (n : ℕ) : poch q n * iq q n = 1 :=
  PowerSeries.mul_inv_cancel _ (by rw [constantCoeff_poch]; exact one_ne_zero)

lemma iq_mul_poch (n : ℕ) : iq q n * poch q n = 1 := by
  rw [mul_comm, poch_mul_iq]

lemma iq_succ_mul (n : ℕ) : iq q (n + 1) * (1 - q ^ (n + 1)) = iq q n := by
  calc iq q (n + 1) * (1 - q ^ (n + 1))
      = iq q (n + 1) * (poch q n * (1 - q ^ (n + 1))) * iq q n := by
        rw [show iq q (n + 1) * (poch q n * (1 - q ^ (n + 1))) * iq q n
          = iq q (n + 1) * (1 - q ^ (n + 1)) * (poch q n * iq q n) by ring, poch_mul_iq, mul_one]
    _ = iq q n := by rw [← poch_succ, iq_mul_poch, one_mul]

lemma iq_eq (n : ℕ) : iq q n = iq q (n + 1) - q ^ (n + 1) * iq q (n + 1) := by
  rw [← iq_succ_mul]; ring

lemma constantCoeff_iq (n : ℕ) : constantCoeff (iq q n) = 1 := by
  have h := congrArg constantCoeff (poch_mul_iq q n)
  simpa [constantCoeff_poch] using h

/-- The Gaussian binomial coefficient `[n;k]_q`, as a power series (zero for `k > n`). -/
noncomputable def qbinom (n k : ℕ) : K⟦X⟧ :=
  if k ≤ n then poch q n * iq q k * iq q (n - k) else 0

omit [Nil q] in
lemma qbinom_of_lt {n k : ℕ} (h : n < k) : qbinom q n k = 0 := by
  simp [qbinom, Nat.not_le.mpr h]

@[simp] lemma qbinom_zero_right (n : ℕ) : qbinom q n 0 = 1 := by
  simp [qbinom, poch_mul_iq]

@[simp] lemma qbinom_self (n : ℕ) : qbinom q n n = 1 := by
  simp [qbinom, poch_mul_iq]

/-- The `q`-Pascal rule `[n+1;k+1] = [n;k+1] + q^{n-k} [n;k]` (valid for all `n`, `k`). -/
lemma qbinom_succ_succ (n k : ℕ) :
    qbinom q (n + 1) (k + 1) = qbinom q n (k + 1) + q ^ (n - k) * qbinom q n k := by
  rcases Nat.lt_or_ge n k with h | h
  · rw [qbinom_of_lt q (show n + 1 < k + 1 by omega), qbinom_of_lt q (show n < k + 1 by omega),
      qbinom_of_lt q h]; simp
  rcases Nat.eq_or_lt_of_le h with h' | h'
  · subst h'
    rw [qbinom_self, qbinom_self, qbinom_of_lt q (Nat.lt_succ_self k)]; simp
  · have e1 : qbinom q (n + 1) (k + 1) = poch q (n + 1) * iq q (k + 1) * iq q (n - k) := by
      simp [qbinom, h]
    have e2 : qbinom q n (k + 1) = poch q n * iq q (k + 1) * iq q (n - (k + 1)) := by
      simp [qbinom, Nat.succ_le_of_lt h']
    have e3 : qbinom q n k = poch q n * iq q k * iq q (n - k) := by
      simp [qbinom, h]
    have e4 : iq q (n - (k + 1)) = iq q (n - k) * (1 - q ^ (n - k)) := by
      have : n - k = (n - (k + 1)) + 1 := by omega
      rw [this, iq_succ_mul]
    have e5 : iq q k = iq q (k + 1) * (1 - q ^ (k + 1)) := (iq_succ_mul q k).symm
    have e6 : q ^ (n + 1) = q ^ (n - k) * q ^ (k + 1) := by
      rw [← pow_add]; congr 1; omega
    rw [e1, e2, e3, e4, e5, poch_succ, e6]
    ring

/-- `[R;r] = (1/(q)_r) ∏_{i<r} (1 - q^{R-r+i+1})` for `r ≤ R`. -/
lemma qbinom_eq_iq_mul {R r : ℕ} (h : r ≤ R) :
    qbinom q R r = iq q r * ∏ i ∈ range r, (1 - q ^ (R - r + i + 1)) := by
  have hp : poch q R = poch q (R - r) * ∏ i ∈ range r, (1 - q ^ (R - r + i + 1)) := by
    have hR : R = (R - r) + r := by omega
    conv_lhs => rw [hR]
    simp only [poch]
    rw [prod_range_add]
  simp only [qbinom, h, if_true]
  rw [hp]
  calc (poch q (R - r) * ∏ i ∈ range r, (1 - q ^ (R - r + i + 1))) * iq q r * iq q (R - r)
      = (poch q (R - r) * iq q (R - r)) * (iq q r * ∏ i ∈ range r, (1 - q ^ (R - r + i + 1))) := by
        ring
    _ = _ := by rw [poch_mul_iq, one_mul]

/-- A product of factors `≡ 1 mod X^n` is `≡ 1 mod X^n`. -/
lemma X_pow_dvd_one_sub_prod {ι : Type*} (s : Finset ι) (f : ι → K⟦X⟧) (n : ℕ)
    (h : ∀ i ∈ s, (X : K⟦X⟧) ^ n ∣ 1 - f i) : (X : K⟦X⟧) ^ n ∣ 1 - ∏ i ∈ s, f i := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [prod_insert ha]
    have h1 := h a (mem_insert_self a s)
    have h2 := ih (fun i hi => h i (mem_insert_of_mem hi))
    have : 1 - f a * ∏ i ∈ s, f i = (1 - f a) + f a * (1 - ∏ i ∈ s, f i) := by ring
    rw [this]
    exact dvd_add h1 (Dvd.dvd.mul_left h2 _)

/-- `[R;r] ≡ 1/(q)_r mod X^{R-r+1}`. -/
lemma X_pow_dvd_iq_sub_qbinom {R r : ℕ} (h : r ≤ R) :
    (X : K⟦X⟧) ^ (R - r + 1) ∣ iq q r - qbinom q R r := by
  rw [qbinom_eq_iq_mul q h]
  have : iq q r - iq q r * ∏ i ∈ range r, (1 - q ^ (R - r + i + 1))
      = iq q r * (1 - ∏ i ∈ range r, (1 - q ^ (R - r + i + 1))) := by ring
  rw [this]
  apply Dvd.dvd.mul_left
  apply X_pow_dvd_one_sub_prod
  intro i _
  have : (1 : K⟦X⟧) - (1 - q ^ (R - r + i + 1)) = q ^ (R - r + i + 1) := by ring
  rw [this]
  exact X_pow_dvd_pow_of_le q (by omega)

end base

end KRTW
