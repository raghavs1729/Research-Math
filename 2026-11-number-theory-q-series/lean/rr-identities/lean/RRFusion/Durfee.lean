/-
Lemma 1 of `proofs/E8-E6.md` (Durfee rectangle identity)

    1/((q)_a (q)_b) = Σ_{c=0}^{min(a,b)} q^{(a−c)(b−c)} / ((q)_{a−c} (q)_{b−c} (q)_c)

for all natural numbers `a`, `b`, in three forms:

* `durfee_gauss`  : denominators cleared, in an arbitrary commutative ring `R` with `q : R`
                    (so in particular in `ℤ[X]` with `q = X`):
                    `Σ_c q^{(a−c)(b−c)} [a;c]_q [b;c]_q (q;q)_c = 1`;
* `durfee_inv`    : for any family `u n` of inverses of `(q;q)_n` in a commutative ring;
* `durfee_field`  : in a field in which `q` is not a root of unity, with honest divisions;
                    `durfee_ratFunc` is the case of the rational function field `ℚ(q)`.

The proof is by induction (q-Pascal recurrence), not via the q-binomial theorem.
-/
import Mathlib

namespace RRFusion

open Finset

variable {R : Type*} [CommRing R]

/-- `(q;q)_n = ∏_{j=1}^{n} (1 − q^j)`. -/
def qPoch (q : R) (n : ℕ) : R := ∏ j ∈ range n, (1 - q ^ (j + 1))

/-- The Gaussian binomial coefficient `[n;k]_q`, defined by the q-Pascal recurrence
`[n+1;k+1] = [n;k] + q^{k+1} [n;k+1]`.  See `gauss_mul_qPoch_mul_qPoch` for the product formula. -/
def gauss (q : R) : ℕ → ℕ → R
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => gauss q n k + q ^ (k + 1) * gauss q n (k + 1)

/-- The falling q-factorial `(1 − q^a)(1 − q^{a−1}) ⋯ (1 − q^{a−c+1}) = (q)_a / (q)_{a−c}`. -/
def qFalling (q : R) (a c : ℕ) : R := ∏ j ∈ range c, (1 - q ^ (a - j))

@[simp] lemma qPoch_zero (q : R) : qPoch q 0 = 1 := by simp [qPoch]

lemma qPoch_succ (q : R) (n : ℕ) : qPoch q (n + 1) = qPoch q n * (1 - q ^ (n + 1)) := by
  simp [qPoch, prod_range_succ]

@[simp] lemma gauss_zero_right (q : R) (n : ℕ) : gauss q n 0 = 1 := by cases n <;> rfl

@[simp] lemma gauss_zero_succ (q : R) (k : ℕ) : gauss q 0 (k + 1) = 0 := rfl

lemma gauss_succ_succ (q : R) (n k : ℕ) :
    gauss q (n + 1) (k + 1) = gauss q n k + q ^ (k + 1) * gauss q n (k + 1) := rfl

lemma gauss_eq_zero_of_lt (q : R) : ∀ {n k : ℕ}, n < k → gauss q n k = 0
  | 0, k + 1, _ => rfl
  | n + 1, k + 1, h => by
    rw [gauss_succ_succ, gauss_eq_zero_of_lt q (n := n) (k := k) (by omega),
      gauss_eq_zero_of_lt q (n := n) (k := k + 1) (by omega)]
    simp

@[simp] lemma qFalling_zero_right (q : R) (a : ℕ) : qFalling q a 0 = 1 := by simp [qFalling]

lemma qFalling_succ (q : R) (a c : ℕ) :
    qFalling q a (c + 1) = qFalling q a c * (1 - q ^ (a - c)) := by
  simp [qFalling, prod_range_succ]

lemma qFalling_succ_succ (q : R) (a c : ℕ) :
    qFalling q (a + 1) (c + 1) = (1 - q ^ (a + 1)) * qFalling q a c := by
  simp [qFalling, prod_range_succ', mul_comm]

lemma qFalling_eq_zero_of_lt (q : R) {a c : ℕ} (h : a < c) : qFalling q a c = 0 :=
  prod_eq_zero (mem_range.2 h) (by simp)

/-- `(q)_a / (q)_{a−c}` as a product. -/
lemma qFalling_mul_qPoch (q : R) {a c : ℕ} (h : c ≤ a) :
    qFalling q a c * qPoch q (a - c) = qPoch q a := by
  induction c with
  | zero => simp
  | succ c ih =>
    have e : a - c = (a - (c + 1)) + 1 := by omega
    rw [← ih (by omega), qFalling_succ]
    conv_rhs => rw [e, qPoch_succ, ← e]
    ring

/-- `[a;c]_q (q)_c = (q)_a / (q)_{a−c}` (valid for all `a`, `c`; both sides vanish for `c > a`). -/
lemma gauss_mul_qPoch (q : R) : ∀ a c : ℕ, gauss q a c * qPoch q c = qFalling q a c
  | _, 0 => by simp
  | 0, c + 1 => by simp [qFalling_eq_zero_of_lt q (Nat.succ_pos c)]
  | a + 1, c + 1 => by
    rw [gauss_succ_succ, qPoch_succ, add_mul, ← mul_assoc, gauss_mul_qPoch q a c, mul_assoc,
      ← qPoch_succ, gauss_mul_qPoch q a (c + 1), qFalling_succ_succ, qFalling_succ]
    rcases le_or_gt c a with h | h
    · have : q ^ (c + 1) * q ^ (a - c) = q ^ (a + 1) := by
        rw [← pow_add]; congr 1; omega
      linear_combination (-qFalling q a c) * this
    · simp [qFalling_eq_zero_of_lt q h]

/-- The product formula `[n;k]_q (q)_k (q)_{n−k} = (q)_n`, identifying `gauss` with the usual
Gaussian binomial coefficient. -/
theorem gauss_mul_qPoch_mul_qPoch (q : R) {n k : ℕ} (h : k ≤ n) :
    gauss q n k * qPoch q k * qPoch q (n - k) = qPoch q n := by
  rw [gauss_mul_qPoch, qFalling_mul_qPoch q h]

/-- The sum of the Durfee rectangle identity, with denominators cleared. -/
def durfeeSum (q : R) (a b : ℕ) : R :=
  ∑ c ∈ range (a + 1), q ^ ((a - c) * (b - c)) * qFalling q a c * gauss q b c

lemma durfeeSum_zero_left (q : R) (b : ℕ) : durfeeSum q 0 b = 1 := by
  simp [durfeeSum]

lemma durfeeSum_zero_right (q : R) (a : ℕ) : durfeeSum q a 0 = 1 := by
  simp [durfeeSum, sum_range_succ']

lemma durfeeSum_succ_succ (q : R) (a b : ℕ) :
    durfeeSum q (a + 1) (b + 1) =
      (1 - q ^ (a + 1)) * durfeeSum q a b + q ^ (a + 1) * durfeeSum q (a + 1) b := by
  unfold durfeeSum
  rw [sum_range_succ' _ (a + 1), sum_range_succ' _ (a + 1), mul_add, mul_sum, mul_sum, ← add_assoc,
    ← sum_add_distrib]
  congr 1
  · refine sum_congr rfl fun c hc => ?_
    have hc' : c ≤ a := Nat.lt_succ_iff.1 (mem_range.1 hc)
    simp only [Nat.add_sub_add_right, gauss_succ_succ, qFalling_succ_succ]
    rcases lt_or_ge c b with h | h
    · have e : q ^ (c + 1) * q ^ ((a - c) * (b - c)) = q ^ (a + 1) * q ^ ((a - c) * (b - (c + 1))) := by
        rw [← pow_add, ← pow_add]
        congr 1
        obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_lt h
        obtain ⟨e, rfl⟩ := Nat.exists_eq_add_of_le hc'
        have h1 : c + e - c = e := by omega
        have h2 : c + d + 1 - c = d + 1 := by omega
        have h3 : c + d + 1 - (c + 1) = d := by omega
        rw [h1, h2, h3]
        ring
      linear_combination ((1 - q ^ (a + 1)) * qFalling q a c * gauss q b (c + 1)) * e
    · rw [gauss_eq_zero_of_lt q (Nat.lt_succ_of_le h)]
      ring
  · simp only [Nat.sub_zero, qFalling_zero_right, gauss_zero_right, mul_one]
    ring

/-- **Durfee rectangle identity, cleared of denominators** (first form). -/
theorem durfeeSum_eq_one (q : R) : ∀ a b : ℕ, durfeeSum q a b = 1 := by
  intro a b
  induction b generalizing a with
  | zero => exact durfeeSum_zero_right q a
  | succ b ih =>
    cases a with
    | zero => exact durfeeSum_zero_left q _
    | succ a => rw [durfeeSum_succ_succ, ih, ih]; ring

/-- **Lemma 1 (Durfee rectangle identity), polynomial form.**  In any commutative ring, for all
`a b : ℕ`:  `Σ_{c=0}^{min(a,b)} q^{(a−c)(b−c)} [a;c]_q [b;c]_q (q;q)_c = 1`. -/
theorem durfee_gauss (q : R) (a b : ℕ) :
    ∑ c ∈ range (min a b + 1), q ^ ((a - c) * (b - c)) * (gauss q a c * gauss q b c * qPoch q c)
      = 1 := by
  rw [← durfeeSum_eq_one q a b, durfeeSum]
  have hsub : range (min a b + 1) ⊆ range (a + 1) :=
    range_mono (Nat.succ_le_succ (min_le_left a b))
  rw [← sum_subset hsub]
  · refine sum_congr rfl fun c _ => ?_
    rw [← gauss_mul_qPoch q a c]
    ring
  · intro c _ hc
    have : b < c := by
      simp only [mem_range, not_lt] at hc ⊢
      rename_i hc'
      have := mem_range.1 hc'
      omega
    rw [gauss_eq_zero_of_lt q this, mul_zero]

/-- **Lemma 1, inverse form.**  If `u n` is an inverse of `(q;q)_n` for every `n`, then
`u_a u_b = Σ_{c=0}^{min(a,b)} q^{(a−c)(b−c)} u_{a−c} u_{b−c} u_c`. -/
theorem durfee_inv (q : R) (u : ℕ → R) (hu : ∀ n, u n * qPoch q n = 1) (a b : ℕ) :
    u a * u b = ∑ c ∈ range (min a b + 1), q ^ ((a - c) * (b - c)) * (u (a - c) * u (b - c) * u c) := by
  have key : ∀ {n k : ℕ}, k ≤ n → gauss q n k * u n = u k * u (n - k) := by
    intro n k h
    have h1 := gauss_mul_qPoch_mul_qPoch q h
    linear_combination (u n * u k * u (n - k)) * h1 + (u k * u (n - k)) * hu n
      - (gauss q n k * u n) * hu k
      - (gauss q n k * u n * u k * qPoch q k) * hu (n - k)
  calc u a * u b = u a * u b * ∑ c ∈ range (min a b + 1),
        q ^ ((a - c) * (b - c)) * (gauss q a c * gauss q b c * qPoch q c) := by
        rw [durfee_gauss, mul_one]
    _ = _ := by
        rw [mul_sum]
        refine sum_congr rfl fun c hc => ?_
        have hc' : c ≤ min a b := Nat.lt_succ_iff.1 (mem_range.1 hc)
        have ha := key (le_trans hc' (min_le_left a b))
        have hb := key (le_trans hc' (min_le_right a b))
        linear_combination (q ^ ((a - c) * (b - c)) * gauss q b c * qPoch q c * u b) * ha
          + (q ^ ((a - c) * (b - c)) * u c * u (a - c) * qPoch q c) * hb
          + (q ^ ((a - c) * (b - c)) * u (a - c) * u (b - c) * u c) * hu c

/-- **Lemma 1 in a field** in which `q` is not a root of unity: identity (1.1) as printed. -/
theorem durfee_field {K : Type*} [Field K] (q : K) (hq : ∀ n : ℕ, 0 < n → q ^ n ≠ 1) (a b : ℕ) :
    1 / (qPoch q a * qPoch q b) =
      ∑ c ∈ range (min a b + 1),
        q ^ ((a - c) * (b - c)) / (qPoch q (a - c) * qPoch q (b - c) * qPoch q c) := by
  have hne : ∀ n, qPoch q n ≠ 0 := fun n =>
    prod_ne_zero_iff.2 fun j _ => sub_ne_zero.2 (hq (j + 1) j.succ_pos).symm
  have := durfee_inv q (fun n => (qPoch q n)⁻¹) (fun n => inv_mul_cancel₀ (hne n)) a b
  simp only [mul_inv, div_eq_mul_inv]
  simpa using this

open RatFunc in
/-- **Lemma 1 as an identity of rational functions** in `ℚ(q)`. -/
theorem durfee_ratFunc (a b : ℕ) :
    1 / (qPoch (X : RatFunc ℚ) a * qPoch X b) =
      ∑ c ∈ range (min a b + 1),
        X ^ ((a - c) * (b - c)) / (qPoch X (a - c) * qPoch X (b - c) * qPoch X c) := by
  refine durfee_field _ (fun n hn h => ?_) a b
  have h' : algebraMap (Polynomial ℚ) (RatFunc ℚ) (Polynomial.X ^ n) = algebraMap _ _ 1 := by
    simpa using h
  have := congrArg Polynomial.natDegree (algebraMap_injective ℚ h')
  simp at this
  omega

end RRFusion
