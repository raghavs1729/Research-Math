/-
# Correctness of the bit-mask sieve of `AllxLib.Finite.Defs`

Main result (`Allx.sieve_bit_iff`): for `n0 ≤ n < n0 + K`, `2 ≤ n` and `n0 + K ≤ B * B`,
bit `n - n0` of `sieve n0 K B` is clear iff `n` is prime.
-/
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.Nat.Bitwise
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import AllxLib.Finite.Defs

namespace Allx

/-! ### The geometric mask -/

theorem geo_eq (q J : ℕ) : geo q J = (2 ^ (q * J) - 1) / (2 ^ q - 1) := rfl

/-- The recursive form of `geo`. -/
def geoR (q : ℕ) : ℕ → ℕ
  | 0 => 0
  | J + 1 => 2 ^ q * geoR q J + 1

lemma geoR_mul (q J : ℕ) : (2 ^ q - 1) * geoR q J = 2 ^ (q * J) - 1 := by
  induction J with
  | zero => simp [geoR]
  | succ J ih =>
    obtain ⟨a, ha⟩ : ∃ a, 2 ^ q = a + 1 := ⟨2 ^ q - 1, by have := Nat.one_le_two_pow (n := q); omega⟩
    obtain ⟨b, hb⟩ : ∃ b, 2 ^ (q * J) = b + 1 :=
      ⟨2 ^ (q * J) - 1, by have := Nat.one_le_two_pow (n := q * J); omega⟩
    have e : 2 ^ (q * (J + 1)) = (b + 1) * (a + 1) := by rw [Nat.mul_succ, Nat.pow_add, ha, hb]
    rw [e]
    rw [ha, hb] at ih
    simp only [Nat.add_sub_cancel] at ih
    rw [geoR, ha, Nat.add_sub_cancel, Nat.mul_add, Nat.mul_left_comm, ih]
    have : (b + 1) * (a + 1) = ((a + 1) * b + a * 1) + 1 := by ring
    rw [this, Nat.add_sub_cancel]

lemma geo_eq_geoR {q : ℕ} (hq : 1 ≤ q) (J : ℕ) : geo q J = geoR q J := by
  have h2 : 2 ≤ 2 ^ q := by
    calc 2 = 2 ^ 1 := rfl
      _ ≤ 2 ^ q := Nat.pow_le_pow_right (by norm_num) hq
  rw [geo_eq, ← geoR_mul]
  exact Nat.mul_div_cancel_left _ (by omega)

lemma testBit_geoR {q : ℕ} (hq : 1 ≤ q) :
    ∀ J i, (geoR q J).testBit i = true ↔ q ∣ i ∧ i / q < J := by
  intro J
  induction J with
  | zero => intro i; simp [geoR]
  | succ J ih =>
    intro i
    have h1 : 1 < 2 ^ q := Nat.one_lt_two_pow (by omega)
    rw [geoR, Nat.testBit_two_pow_mul_add _ h1]
    split_ifs with h
    · rw [Nat.testBit_one_eq_true_iff_self_eq_zero]
      constructor
      · rintro rfl; simp
      · rintro ⟨hd, -⟩; exact Nat.eq_zero_of_dvd_of_lt hd h
    · have hi : i = (i - q) + q := by omega
      rw [ih (i - q)]
      constructor
      · rintro ⟨hd, hlt⟩
        refine ⟨?_, ?_⟩
        · rw [hi]; exact Nat.dvd_add hd (Nat.dvd_refl q)
        · rw [hi, Nat.add_div_right _ (by omega)]; omega
      · rintro ⟨hd, hlt⟩
        refine ⟨(Nat.dvd_sub hd (Nat.dvd_refl q)), ?_⟩
        rw [hi, Nat.add_div_right _ (by omega)] at hlt; omega

lemma testBit_geo {q : ℕ} (hq : 1 ≤ q) (J i : ℕ) :
    (geo q J).testBit i = true ↔ q ∣ i ∧ i / q < J := by
  rw [geo_eq_geoR hq]; exact testBit_geoR hq J i

/-! ### One mask -/

theorem firstMul_eq (n0 q : ℕ) : firstMul n0 q = q * max 2 ((n0 + q - 1) / q) := by
  show q * (bif Nat.ble 2 ((n0 + q - 1) / q) then (n0 + q - 1) / q else 2) = _
  congr 1
  cases h : Nat.ble 2 ((n0 + q - 1) / q) with
  | true => have := Nat.le_of_ble_eq_true h; simp only [cond_true]; omega
  | false =>
    have : ¬ 2 ≤ (n0 + q - 1) / q := fun h' => by simp [Nat.ble_eq_true_of_le h'] at h
    simp only [cond_false]; omega

theorem mask_eq (n0 K q : ℕ) :
    mask n0 K q = geo q (K / q + 1) <<< (firstMul n0 q - n0) := rfl

lemma le_firstMul {n0 q : ℕ} (hq : 1 ≤ q) : n0 ≤ firstMul n0 q := by
  rw [firstMul_eq]
  have h1 : n0 + q - 1 < q * ((n0 + q - 1) / q + 1) := Nat.lt_mul_div_succ _ (by omega)
  have h2 : q * ((n0 + q - 1) / q) ≤ q * max 2 ((n0 + q - 1) / q) :=
    Nat.mul_le_mul_left _ (le_max_right _ _)
  rw [Nat.mul_succ] at h1
  omega

/-- Bit `n - n0` of `mask n0 K q` is set iff `q` is a proper divisor of `n`. -/
lemma testBit_mask {n0 K q n : ℕ} (hq : 1 ≤ q) (h0 : n0 ≤ n) (hK : n < n0 + K) :
    (mask n0 K q).testBit (n - n0) = true ↔ q ∣ n ∧ q < n := by
  have hs := le_firstMul (n0 := n0) hq
  have hsd : q ∣ firstMul n0 q := by rw [firstMul_eq]; exact Nat.dvd_mul_right _ _
  have hs2 : 2 * q ≤ firstMul n0 q := by
    rw [firstMul_eq, Nat.mul_comm 2 q]; exact Nat.mul_le_mul_left _ (le_max_left _ _)
  rw [mask_eq, Nat.testBit_shiftLeft, Bool.and_eq_true, decide_eq_true_eq, testBit_geo hq]
  have e : n - n0 - (firstMul n0 q - n0) = n - firstMul n0 q := by omega
  rw [e]
  constructor
  · rintro ⟨hle, hd, -⟩
    have hsn : firstMul n0 q ≤ n := by omega
    refine ⟨?_, by omega⟩
    have := Nat.dvd_add hd hsd
    rwa [Nat.sub_add_cancel hsn] at this
  · rintro ⟨hd, hlt⟩
    obtain ⟨c, rfl⟩ := hd
    have hc2 : 2 ≤ c := by
      by_contra hc
      have : c ≤ 1 := by omega
      have := Nat.mul_le_mul_left q this
      omega
    have hcm : (n0 + q - 1) / q ≤ c := by
      have : n0 + q - 1 < q * (c + 1) := by rw [Nat.mul_succ]; omega
      have := (Nat.div_lt_iff_lt_mul (by omega : 0 < q)).2 (by rwa [Nat.mul_comm] at this)
      omega
    have hsn : firstMul n0 q ≤ q * c := by
      rw [firstMul_eq]; exact Nat.mul_le_mul_left _ (max_le hc2 hcm)
    refine ⟨by omega, Nat.dvd_sub (Nat.dvd_mul_right _ _) hsd, ?_⟩
    have h1 : (q * c - firstMul n0 q) / q ≤ K / q := Nat.div_le_div_right (by omega)
    omega

/-! ### The OR-loop -/

theorem sieveLoop_zero (n0 K q acc : ℕ) : sieveLoop n0 K 0 q acc = acc := rfl

theorem sieveLoop_succ (n0 K fuel q acc : ℕ) :
    sieveLoop n0 K (fuel + 1) q acc = sieveLoop n0 K fuel (q + 1) (acc ||| mask n0 K q) := by
  show (bif Nat.beq (acc ||| mask n0 K q) 0 then sieveLoop n0 K fuel (q + 1) 0
    else sieveLoop n0 K fuel (q + 1) (acc ||| mask n0 K q)) = _
  cases h : Nat.beq (acc ||| mask n0 K q) 0 with
  | true => rw [cond_true, Nat.eq_of_beq_eq_true h]
  | false => rfl

lemma testBit_sieveLoop (n0 K i : ℕ) : ∀ fuel q acc,
    (sieveLoop n0 K fuel q acc).testBit i = true ↔
      (acc.testBit i = true ∨ ∃ q', q ≤ q' ∧ q' < q + fuel ∧ (mask n0 K q').testBit i = true) := by
  intro fuel
  induction fuel with
  | zero =>
    intro q acc
    rw [sieveLoop_zero]
    constructor
    · exact Or.inl
    · rintro (h | ⟨q', h1, h2, -⟩)
      · exact h
      · omega
  | succ fuel ih =>
    intro q acc
    rw [sieveLoop_succ, ih, Nat.testBit_lor, Bool.or_eq_true]
    constructor
    · rintro ((h | h) | ⟨q', h1, h2, h3⟩)
      · exact Or.inl h
      · exact Or.inr ⟨q, le_rfl, by omega, h⟩
      · exact Or.inr ⟨q', by omega, by omega, h3⟩
    · rintro (h | ⟨q', h1, h2, h3⟩)
      · exact Or.inl (Or.inl h)
      · rcases Nat.eq_or_lt_of_le h1 with rfl | hlt
        · exact Or.inl (Or.inr h3)
        · exact Or.inr ⟨q', by omega, by omega, h3⟩

/-- Bit `n - n0` of the sieve is set iff `n` has a proper divisor `2 ≤ q < B`. -/
lemma testBit_sieve {n0 K B n : ℕ} (h0 : n0 ≤ n) (hK : n < n0 + K) :
    (sieve n0 K B).testBit (n - n0) = true ↔ ∃ q, 2 ≤ q ∧ q < B ∧ q ∣ n ∧ q < n := by
  show (sieveLoop n0 K (B - 2) 2 0).testBit (n - n0) = true ↔ _
  rw [testBit_sieveLoop]
  constructor
  · rintro (h | ⟨q, h1, h2, h3⟩)
    · simp at h
    · exact ⟨q, h1, by omega, (testBit_mask (by omega) h0 hK).1 h3⟩
  · rintro ⟨q, h1, h2, h3⟩
    exact Or.inr ⟨q, h1, by omega, (testBit_mask (by omega) h0 hK).2 h3⟩

/-! ### Primality -/

lemma prime_iff_no_small_divisor {B n : ℕ} (h2 : 2 ≤ n) (hB : n < B * B) :
    n.Prime ↔ ¬ ∃ q, 2 ≤ q ∧ q < B ∧ q ∣ n ∧ q < n := by
  constructor
  · rintro hp ⟨q, hq2, -, hd, hlt⟩
    rcases (Nat.dvd_prime hp).1 hd with h | h <;> omega
  · intro H
    by_contra hnp
    have hq : n.minFac.Prime := Nat.minFac_prime (by omega)
    have hsq : n.minFac ^ 2 ≤ n := Nat.minFac_sq_le_self (by omega) hnp
    rw [pow_two] at hsq
    have hq2 := hq.two_le
    refine H ⟨n.minFac, hq2, ?_, Nat.minFac_dvd n, ?_⟩
    · by_contra hge
      have : B * B ≤ n.minFac * n.minFac := Nat.mul_le_mul (by omega) (by omega)
      omega
    · have : n.minFac * 2 ≤ n.minFac * n.minFac := Nat.mul_le_mul_left _ hq2
      omega

/-- The raw bit test used by the scan. -/
lemma bit_test_iff (C i : ℕ) :
    Nat.beq (Nat.land (Nat.shiftRight C i) 1) 0 = true ↔ C.testBit i = false := by
  have e : Nat.land (Nat.shiftRight C i) 1 = C / 2 ^ i % 2 := by
    show (C >>> i) &&& 1 = _
    rw [Nat.and_one_is_mod, Nat.shiftRight_eq_div_pow]
  rw [e, Nat.testBit_eq_decide_div_mod_eq, decide_eq_false_iff_not]
  constructor
  · intro h; have := Nat.eq_of_beq_eq_true h; omega
  · intro h
    have : C / 2 ^ i % 2 = 0 := by omega
    rw [this]; rfl

/-- **Correctness of the sieve.** -/
theorem sieve_bit_iff {n0 K B n : ℕ} (h0 : n0 ≤ n) (hK : n < n0 + K) (h2 : 2 ≤ n)
    (hB : n0 + K ≤ B * B) :
    Nat.beq (Nat.land (Nat.shiftRight (sieve n0 K B) (n - n0)) 1) 0 = true ↔ n.Prime := by
  rw [bit_test_iff, prime_iff_no_small_divisor (B := B) h2 (by omega),
    ← testBit_sieve (B := B) h0 hK]
  simp

end Allx
