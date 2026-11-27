/-
Tau/Basic.lean — the recurrence map `T (U, V) = (V, 3V − U)`, the invariant binary quadratic form
`Q (U, V) = U² − 3UV + V²`, and compatibility of `T` with ring homomorphisms (in particular reduction
modulo `q`), for integer powers `T ^ n`, `n : ℤ`.
-/
import Mathlib

namespace Tau

/-- The quadratic form `U² − 3UV + V²` (discriminant 5). -/
def Qf {R : Type*} [CommRing R] (p : R × R) : R := p.1 ^ 2 - 3 * p.1 * p.2 + p.2 ^ 2

/-- The map `(U, V) ↦ (V, 3V − U)` as a permutation of `R × R` (inverse `(U, V) ↦ (3U − V, U)`). -/
def T (R : Type*) [CommRing R] : Equiv.Perm (R × R) where
  toFun p := (p.2, 3 * p.2 - p.1)
  invFun p := (3 * p.1 - p.2, p.1)
  left_inv p := by ext <;> simp
  right_inv p := by ext <;> simp

variable {R S : Type*} [CommRing R] [CommRing S]

@[simp] lemma T_apply (p : R × R) : T R p = (p.2, 3 * p.2 - p.1) := rfl

@[simp] lemma T_inv_apply (p : R × R) : (T R)⁻¹ p = (3 * p.1 - p.2, p.1) := rfl

lemma Qf_T (p : R × R) : Qf (T R p) = Qf p := by
  simp only [Qf, T_apply]; ring

lemma Qf_T_inv (p : R × R) : Qf ((T R)⁻¹ p) = Qf p := by
  simp only [Qf, T_inv_apply]; ring

lemma zpow_succ_apply (n : ℤ) (p : R × R) : (T R ^ (n + 1)) p = T R ((T R ^ n) p) := by
  rw [add_comm, zpow_one_add, Equiv.Perm.mul_apply]

lemma zpow_pred_apply (n : ℤ) (p : R × R) : (T R ^ (n - 1)) p = (T R)⁻¹ ((T R ^ n) p) := by
  have h := zpow_succ_apply (R := R) (n - 1) p
  rw [sub_add_cancel] at h
  rw [h]; simp

lemma Qf_zpow (n : ℤ) (p : R × R) : Qf ((T R ^ n) p) = Qf p := by
  induction n using Int.induction_on with
  | zero => simp
  | succ k ih => rw [zpow_succ_apply, Qf_T, ih]
  | pred k ih => rw [zpow_pred_apply, Qf_T_inv, ih]

/-- Componentwise application of a ring homomorphism. -/
def pmap (f : R →+* S) (p : R × R) : S × S := (f p.1, f p.2)

lemma pmap_T (f : R →+* S) (p : R × R) : pmap f (T R p) = T S (pmap f p) := by
  simp [pmap, map_sub, map_mul, map_ofNat]

lemma pmap_T_inv (f : R →+* S) (p : R × R) : pmap f ((T R)⁻¹ p) = (T S)⁻¹ (pmap f p) := by
  simp [pmap, map_sub, map_mul, map_ofNat]

lemma pmap_zpow (f : R →+* S) (n : ℤ) (p : R × R) :
    pmap f ((T R ^ n) p) = (T S ^ n) (pmap f p) := by
  induction n using Int.induction_on with
  | zero => simp
  | succ k ih => rw [zpow_succ_apply, zpow_succ_apply, pmap_T, ih]
  | pred k ih => rw [zpow_pred_apply, zpow_pred_apply, pmap_T_inv, ih]

/-- The first coordinate of `T ^ (n + 1)` is the second coordinate of `T ^ n`. -/
lemma fst_zpow_succ (n : ℤ) (p : R × R) : ((T R ^ (n + 1)) p).1 = ((T R ^ n) p).2 := by
  rw [zpow_succ_apply]; rfl

/-- `Z_{n+1} = 3 Z_n − Z_{n−1}` in pair form: `T ^ n p = (Z_n, 3 Z_n − Z_{n−1})`. -/
lemma snd_zpow (n : ℤ) (p : R × R) :
    ((T R ^ n) p).2 = 3 * ((T R ^ n) p).1 - ((T R ^ (n - 1)) p).1 := by
  have h := zpow_succ_apply (R := R) (n - 1) p
  rw [sub_add_cancel] at h
  rw [h]; simp

end Tau
