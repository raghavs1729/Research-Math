/-
KRTW/Lemma4.lean — the key identity (Lemma 4 of `proofs/KRTW-general.md`) in its two infinite forms
(**) and (*∞).

The proof given here differs from the write-up: instead of the q-binomial theorem and Lemmas 2, 3 it
uses three two-term recurrences (in `R`, with `x ↦ xq`, and in `M`), each checked termwise.
Throughout `q = X`, and `x : K⟦X⟧` is an arbitrary power series.
-/
import KRTW.Summation

namespace KRTW

open PowerSeries Finset

open scoped PowerSeries.WithPiTopology

set_option linter.unusedSectionVars false

variable {K : Type*} [Field K] [TopologicalSpace K] [DiscreteTopology K]

/-! ### The coefficients `c(R,n) = Σ_j (-1)^j q^{Rj} / ((q)_j (q)_{n-j})` -/

/-- `c(R,n) = Σ_{j=0}^n (-1)^j q^{Rj} / ((q)_j (q)_{n-j})` (right side of Lemma 2). -/
noncomputable def cR (R n : ℕ) : K⟦X⟧ :=
  ∑ j ∈ range (n + 1), (-1) ^ j * X ^ (R * j) * iq X j * iq X (n - j)

omit [TopologicalSpace K] [DiscreteTopology K] in
@[simp] lemma cR_zero_right (R : ℕ) : cR (K := K) R 0 = 1 := by
  simp [cR]

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma cR_succ_succ (R n : ℕ) :
    cR (K := K) (R + 1) (n + 1) = cR R (n + 1) + X ^ R * cR R n := by
  simp only [cR]
  rw [sum_range_succ' _ (n + 1), sum_range_succ' _ (n + 1), mul_sum, add_right_comm,
    ← sum_add_distrib]
  congr 1
  apply sum_congr rfl
  intro j _
  rw [Nat.add_sub_add_right, ← iq_succ_mul X j]
  ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma one_sub_mul_cR (R n : ℕ) :
    (1 - X ^ (n + 1)) * cR (K := K) R (n + 1) = cR (R + 1) n - X ^ R * cR R n := by
  set g1 : ℕ → K⟦X⟧ := fun j =>
    (-1) ^ j * X ^ (R * j) * (iq X j * (1 - X ^ j)) * iq X (n + 1 - j) with hg1
  set g2 : ℕ → K⟦X⟧ := fun j =>
    (-1) ^ j * X ^ (R * j) * X ^ j * iq X j * (iq X (n + 1 - j) * (1 - X ^ (n + 1 - j))) with hg2
  have key : ∀ j ∈ range (n + 1 + 1),
      (1 - X ^ (n + 1)) * ((-1) ^ j * X ^ (R * j) * iq X j * iq X (n + 1 - j)) = g1 j + g2 j := by
    intro j hj
    have hj' : j ≤ n + 1 := Nat.lt_succ_iff.mp (mem_range.mp hj)
    have : (X : K⟦X⟧) ^ (n + 1) = X ^ j * X ^ (n + 1 - j) := by
      rw [← pow_add]; congr 1; omega
    simp only [hg1, hg2]
    rw [this]; ring
  simp only [cR]
  rw [mul_sum, sum_congr rfl key, sum_add_distrib, sum_range_succ' g1, sum_range_succ g2]
  have e1 : g1 0 = 0 := by simp [hg1]
  have e2 : g2 (n + 1) = 0 := by simp [hg2]
  rw [e1, e2, add_zero, add_zero, mul_sum, sub_eq_add_neg, add_comm, ← sum_neg_distrib]
  congr 1
  · apply sum_congr rfl
    intro j hj
    have hj' : j ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hj)
    simp only [hg2]
    have : n + 1 - j = (n - j) + 1 := by omega
    rw [this, iq_succ_mul]
    ring
  · apply sum_congr rfl
    intro j _
    simp only [hg1]
    rw [iq_succ_mul, Nat.add_sub_add_right]
    ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma one_sub_mul_cR_zero (n : ℕ) :
    (1 - X ^ (n + 2)) * cR (K := K) 0 (n + 2) = cR 0 n := by
  have h1 := one_sub_mul_cR (K := K) 0 (n + 1)
  have h2 := cR_succ_succ (K := K) 0 n
  rw [h1, h2]; simp

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- Lemma 2 at `R = 0`: `Σ_j (-1)^j/((q)_j (q)_{n-j})` is `1/(q²;q²)_{n/2}` for even `n`, `0` for odd `n`. -/
lemma cR_zero_even_odd (b : ℕ) :
    cR (K := K) 0 (2 * b) = iq (X ^ 2) b ∧ cR (K := K) 0 (2 * b + 1) = 0 := by
  induction b with
  | zero =>
    refine ⟨by simp, ?_⟩
    have h := one_sub_mul_cR (K := K) 0 0
    simp only [cR_zero_right, pow_zero, mul_one, sub_self, zero_add] at h
    rcases mul_eq_zero.mp h with h | h
    · exact absurd h (one_sub_pow_ne_zero X 0)
    · simpa using h
  | succ b ih =>
    constructor
    · have h := one_sub_mul_cR_zero (K := K) (2 * b)
      rw [ih.1, ← iq_succ_mul (X ^ 2) b, ← pow_mul] at h
      have e : 2 * (b + 1) = 2 * b + 2 := by ring
      rw [e]
      apply mul_left_cancel₀ (one_sub_pow_ne_zero X (2 * b + 1))
      rw [h]
      have e' : 2 * (b + 1) = 2 * b + 1 + 1 := by ring
      rw [e']; ring
    · have h := one_sub_mul_cR_zero (K := K) (2 * b + 1)
      rw [ih.2] at h
      have e : 2 * (b + 1) + 1 = 2 * b + 1 + 2 := by ring
      rw [e]
      rcases mul_eq_zero.mp h with h | h
      · exact absurd h (one_sub_pow_ne_zero X (2 * b + 1 + 1))
      · exact h

/-! ### Euler's series for `(xq^{a+1};q)_∞` -/

/-- `Einf a x = (x q^{a+1}; q)_∞`, defined by Euler's series
`Σ_l (-1)^l q^{l(l+1)/2 + a l} x^l / (q)_l`. -/
noncomputable def Einf (a : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑' l, X ^ (tri l + a * l) * ((-1) ^ l * x ^ l * iq X l)

lemma Einf_summable (a : ℕ) (x : K⟦X⟧) :
    Summable (fun l => X ^ (tri l + a * l) * ((-1) ^ l * x ^ l * iq X l)) :=
  summable_nat (fun l =>
    Dvd.dvd.mul_right (pow_dvd_pow X (le_trans (le_tri l) (Nat.le_add_right _ _))) _)

lemma Einf_mul_pow (a n : ℕ) (x : K⟦X⟧) : Einf a (x * X ^ n) = Einf (a + n) x := by
  simp only [Einf]
  apply tsum_congr
  intro l
  ring

lemma Einf_mul_X (a : ℕ) (x : K⟦X⟧) : Einf a (x * X) = Einf (a + 1) x := by
  have := Einf_mul_pow a 1 x
  rwa [pow_one] at this

/-- `(xq^{a+1};q)_∞ = (1 - xq^{a+1}) (xq^{a+2};q)_∞`. -/
lemma Einf_eq (a : ℕ) (x : K⟦X⟧) : Einf a x = (1 - x * X ^ (a + 1)) * Einf (a + 1) x := by
  have hs := Einf_summable (K := K)
  have h0 : ∀ a : ℕ, Einf a x = 1 + ∑' l, X ^ (tri (l + 1) + a * (l + 1)) *
      ((-1) ^ (l + 1) * x ^ (l + 1) * iq X (l + 1)) := by
    intro a
    rw [Einf, (hs a x).tsum_eq_zero_add]
    simp
  have hs1 : ∀ a : ℕ, Summable (fun l => X ^ (tri (l + 1) + a * (l + 1)) *
      ((-1) ^ (l + 1) * x ^ (l + 1) * iq X (l + 1))) := by
    intro a
    exact summable_nat (fun l => Dvd.dvd.mul_right (pow_dvd_pow X
      (le_trans (Nat.le_succ l) (le_trans (le_tri (l + 1)) (Nat.le_add_right _ _)))) _)
  have hs2 : Summable (fun l => x * X ^ (a + 1) * (X ^ (tri l + (a + 1) * l) *
      ((-1) ^ l * x ^ l * iq X l))) := by
    exact summable_nat (fun l => Dvd.dvd.mul_left (Dvd.dvd.mul_right (pow_dvd_pow X
      (le_trans (le_tri l) (Nat.le_add_right _ _))) _) _)
  have e2 : x * X ^ (a + 1) * Einf (a + 1) x = ∑' l, x * X ^ (a + 1) *
      (X ^ (tri l + (a + 1) * l) * ((-1) ^ l * x ^ l * iq X l)) := by
    rw [Einf]
    exact ((hs (a + 1) x).tsum_mul_left (x * X ^ (a + 1))).symm
  rw [sub_mul, one_mul, e2, h0 a, h0 (a + 1), add_sub_assoc, ← (hs1 (a + 1)).tsum_sub hs2]
  congr 1
  apply tsum_congr
  intro l
  rw [← iq_succ_mul X l, tri_succ]
  ring

lemma X_pow_dvd_Einf_sub_one (a : ℕ) (x : K⟦X⟧) : (X : K⟦X⟧) ^ (a + 1) ∣ Einf a x - 1 := by
  rw [Einf, (Einf_summable a x).tsum_eq_zero_add]
  simp only [tri_zero, mul_zero, add_zero, pow_zero, mul_one, iq_zero, add_sub_cancel_left]
  apply X_pow_dvd_tsum
  intro l
  apply Dvd.dvd.mul_right
  apply pow_dvd_pow
  have h1 := le_tri (l + 1)
  have h2 : a ≤ a * (l + 1) := Nat.le_mul_of_pos_right a (Nat.succ_pos l)
  omega

/-! ### The common form `Z(R,M;x)` and the right side `U(R,M;x)` of (**) -/

/-- `Z(R,M;x) = Σ_{n ≤ M} x^n q^{n(n+1)/2} c(R,n) / (q)_{M-n}` (the form (2.3) of the write-up). -/
noncomputable def Zf (R M : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑ n ∈ range (M + 1), X ^ tri n * (x ^ n * cR R n * iq X (M - n))

/-- Right side of (**): `Σ_{L ≤ M} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_∞ / ((q)_L (q)_{M-L})`. -/
noncomputable def Uss (R M : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑ L ∈ range (M + 1), X ^ tri L * (x ^ L * iq X L * iq X (M - L) * Einf (R + L) x)

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma Zf_succ (R M : ℕ) (x : K⟦X⟧) :
    (1 - X ^ (M + 1)) * Zf R (M + 1) x
      = (1 - x * X ^ (R + M + 1)) * Zf R M x + x * X ^ (M + 1) * Zf (R + 1) M x := by
  set h1 : ℕ → K⟦X⟧ := fun n =>
    X ^ tri n * (x ^ n * cR R n * (iq X (M + 1 - n) * (1 - X ^ (M + 1 - n)))) with hh1
  set h2 : ℕ → K⟦X⟧ := fun n =>
    X ^ tri n * (x ^ n * X ^ (M + 1 - n) * iq X (M + 1 - n) * ((1 - X ^ n) * cR R n)) with hh2
  have key : ∀ n ∈ range (M + 1 + 1),
      (1 - X ^ (M + 1)) * (X ^ tri n * (x ^ n * cR R n * iq X (M + 1 - n))) = h1 n + h2 n := by
    intro n hn
    have hn' : n ≤ M + 1 := Nat.lt_succ_iff.mp (mem_range.mp hn)
    have : (X : K⟦X⟧) ^ (M + 1) = X ^ (M + 1 - n) * X ^ n := by
      rw [← pow_add]; congr 1; omega
    simp only [hh1, hh2]
    rw [this]; ring
  simp only [Zf]
  rw [mul_sum, sum_congr rfl key, sum_add_distrib, sum_range_succ h1, sum_range_succ' h2]
  have e1 : h1 (M + 1) = 0 := by simp [hh1]
  have e2 : h2 0 = 0 := by simp [hh2]
  rw [e1, e2, add_zero, add_zero]
  have s1 : ∑ n ∈ range (M + 1), h1 n
      = ∑ n ∈ range (M + 1), X ^ tri n * (x ^ n * cR R n * iq X (M - n)) := by
    apply sum_congr rfl
    intro n hn
    have hn' : n ≤ M := Nat.lt_succ_iff.mp (mem_range.mp hn)
    simp only [hh1]
    have : M + 1 - n = (M - n) + 1 := by omega
    rw [this, iq_succ_mul]
  have s2 : ∑ n ∈ range (M + 1), h2 (n + 1)
      = x * X ^ (M + 1) * ∑ n ∈ range (M + 1), X ^ tri n * (x ^ n * cR (R + 1) n * iq X (M - n))
        - x * X ^ (R + M + 1) * ∑ n ∈ range (M + 1), X ^ tri n * (x ^ n * cR R n * iq X (M - n)) := by
    rw [mul_sum, mul_sum, ← sum_sub_distrib]
    apply sum_congr rfl
    intro n hn
    have hn' : n ≤ M := Nat.lt_succ_iff.mp (mem_range.mp hn)
    simp only [hh2]
    rw [one_sub_mul_cR, Nat.add_sub_add_right, tri_succ]
    have e : (X : K⟦X⟧) ^ (M + 1) = X ^ (M - n) * X ^ (n + 1) := by
      rw [← pow_add]; congr 1; omega
    have e' : (X : K⟦X⟧) ^ (R + M + 1) = X ^ R * (X ^ (M - n) * X ^ (n + 1)) := by
      rw [← e, ← pow_add]; congr 1
    rw [e', e]
    ring
  rw [s1, s2]
  ring

lemma Uss_succ_M (R M : ℕ) (x : K⟦X⟧) :
    (1 - X ^ (M + 1)) * Uss R (M + 1) x = Uss R M x + x * X ^ (M + 1) * Uss (R + 1) M x := by
  set h1 : ℕ → K⟦X⟧ := fun L =>
    X ^ tri L * (x ^ L * iq X L * (iq X (M + 1 - L) * (1 - X ^ (M + 1 - L))) * Einf (R + L) x)
    with hh1
  set h2 : ℕ → K⟦X⟧ := fun L =>
    X ^ tri L * (x ^ L * X ^ (M + 1 - L) * iq X (M + 1 - L) * (iq X L * (1 - X ^ L))
      * Einf (R + L) x) with hh2
  have key : ∀ L ∈ range (M + 1 + 1),
      (1 - X ^ (M + 1)) * (X ^ tri L * (x ^ L * iq X L * iq X (M + 1 - L) * Einf (R + L) x))
        = h1 L + h2 L := by
    intro L hL
    have hL' : L ≤ M + 1 := Nat.lt_succ_iff.mp (mem_range.mp hL)
    have : (X : K⟦X⟧) ^ (M + 1) = X ^ (M + 1 - L) * X ^ L := by
      rw [← pow_add]; congr 1; omega
    simp only [hh1, hh2]
    rw [this]; ring
  simp only [Uss]
  rw [mul_sum, sum_congr rfl key, sum_add_distrib, sum_range_succ h1, sum_range_succ' h2]
  have e1 : h1 (M + 1) = 0 := by simp [hh1]
  have e2 : h2 0 = 0 := by simp [hh2]
  rw [e1, e2, add_zero, add_zero]
  congr 1
  · apply sum_congr rfl
    intro L hL
    have hL' : L ≤ M := Nat.lt_succ_iff.mp (mem_range.mp hL)
    simp only [hh1]
    have : M + 1 - L = (M - L) + 1 := by omega
    rw [this, iq_succ_mul]
  · rw [mul_sum]
    apply sum_congr rfl
    intro L hL
    have hL' : L ≤ M := Nat.lt_succ_iff.mp (mem_range.mp hL)
    simp only [hh2]
    rw [iq_succ_mul, Nat.add_sub_add_right, tri_succ]
    have e : (X : K⟦X⟧) ^ (M + 1) = X ^ (M - L) * X ^ (L + 1) := by
      rw [← pow_add]; congr 1; omega
    have e' : R + (L + 1) = R + 1 + L := by omega
    rw [e, e']
    ring

/-- `U(R,M;x) = (xq^{R+M+1};q)_∞ Z(R,M;x)`. -/
lemma Uss_eq (M : ℕ) : ∀ (R : ℕ) (x : K⟦X⟧), Uss R M x = Einf (R + M) x * Zf R M x := by
  induction M with
  | zero => intro R x; simp [Uss, Zf]
  | succ M ih =>
    intro R x
    apply mul_left_cancel₀ (one_sub_pow_ne_zero (X : K⟦X⟧) M)
    rw [Uss_succ_M, ih R, ih (R + 1)]
    have hZ := Zf_succ R M x
    have hE := Einf_eq (R + M) x
    have e1 : R + 1 + M = R + M + 1 := by omega
    have e2 : R + (M + 1) = R + M + 1 := by omega
    rw [e1, e2, hE]
    linear_combination (-Einf (R + M + 1) x) * hZ

/-! ### The left side of (**) and (*∞): the recurrence in `R` -/

/-- `W_G(R;x) = Σ_{r ≤ R} [R;r] x^r q^{r²} G(xq^r)` for an arbitrary function `G`. -/
noncomputable def Wf (G : K⟦X⟧ → K⟦X⟧) (R : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑ r ∈ range (R + 1), X ^ (r ^ 2) * (qbinom X R r * x ^ r * G (x * X ^ r))

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma Wf_succ (G : K⟦X⟧ → K⟦X⟧) (R : ℕ) (x : K⟦X⟧) :
    Wf G (R + 1) x = Wf G R x + x * X ^ (R + 1) * Wf G R (x * X) := by
  have h1 : Wf G (R + 1) x = ∑ r ∈ range (R + 1),
      X ^ ((r + 1) ^ 2) * (qbinom X (R + 1) (r + 1) * x ^ (r + 1) * G (x * X ^ (r + 1))) + G x := by
    rw [Wf, sum_range_succ']; simp
  have h2 : Wf G R x = ∑ r ∈ range (R + 1),
      X ^ ((r + 1) ^ 2) * (qbinom X R (r + 1) * x ^ (r + 1) * G (x * X ^ (r + 1))) + G x := by
    have : Wf G R x = ∑ r ∈ range (R + 1 + 1),
        X ^ (r ^ 2) * (qbinom X R r * x ^ r * G (x * X ^ r)) := by
      rw [sum_range_succ _ (R + 1), qbinom_of_lt X (Nat.lt_succ_self R)]; simp [Wf]
    rw [this, sum_range_succ']; simp
  have h3 : x * X ^ (R + 1) * Wf G R (x * X) = ∑ r ∈ range (R + 1),
      X ^ ((r + 1) ^ 2) * (X ^ (R - r) * qbinom X R r * x ^ (r + 1) * G (x * X ^ (r + 1))) := by
    rw [Wf, mul_sum]
    apply sum_congr rfl
    intro r hr
    have hr' : r ≤ R := Nat.lt_succ_iff.mp (mem_range.mp hr)
    have e : (X : K⟦X⟧) ^ ((r + 1) ^ 2) * X ^ (R - r) = X ^ (R + 1) * X ^ (r ^ 2) * X ^ r := by
      rw [← pow_add, ← pow_add, ← pow_add]; congr 1
      have : (r + 1) ^ 2 = r ^ 2 + 2 * r + 1 := by ring
      omega
    have a : x * X * X ^ r = x * X ^ (r + 1) := by ring
    rw [a]
    linear_combination (-(qbinom X R r * x ^ (r + 1) * G (x * X ^ (r + 1)))) * e
  rw [h1, h2, h3, add_right_comm, ← sum_add_distrib]
  congr 1
  apply sum_congr rfl
  intro r _
  rw [qbinom_succ_succ]
  ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma Wf_zero (G : K⟦X⟧ → K⟦X⟧) (x : K⟦X⟧) : Wf G 0 x = G x := by
  simp [Wf]

omit [TopologicalSpace K] [DiscreteTopology K] in
/-- Two solutions of `F(R+1;x) = F(R;x) + x q^{R+1} F(R;xq)` with the same value at `R = 0` agree. -/
lemma rec_unique (F1 F2 : ℕ → K⟦X⟧ → K⟦X⟧)
    (h1 : ∀ R x, F1 (R + 1) x = F1 R x + x * X ^ (R + 1) * F1 R (x * X))
    (h2 : ∀ R x, F2 (R + 1) x = F2 R x + x * X ^ (R + 1) * F2 R (x * X))
    (h0 : ∀ x, F1 0 x = F2 0 x) : ∀ R x, F1 R x = F2 R x := by
  intro R
  induction R with
  | zero => exact h0
  | succ R ih => intro x; rw [h1, h2, ih x, ih (x * X)]

lemma Einf_succ' (a : ℕ) (x : K⟦X⟧) :
    Einf (a + 1) x = Einf a x + x * X ^ (a + 1) * Einf a (x * X) := by
  rw [Einf_mul_X, Einf_eq a x]; ring

lemma Uss_succ_R (R M : ℕ) (x : K⟦X⟧) :
    Uss (R + 1) M x = Uss R M x + x * X ^ (R + 1) * Uss R M (x * X) := by
  simp only [Uss]
  rw [mul_sum, ← sum_add_distrib]
  apply sum_congr rfl
  intro L _
  have e : R + 1 + L = R + L + 1 := by omega
  rw [e, Einf_succ']
  ring

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma sum_range_even (f : ℕ → K⟦X⟧) (hodd : ∀ b, f (2 * b + 1) = 0) (M : ℕ) :
    ∑ n ∈ range (M + 1), f n = ∑ b ∈ range (M + 1), if 2 * b ≤ M then f (2 * b) else 0 := by
  rw [← sum_filter, ← sum_filter_of_ne (p := fun n => Even n)]
  · symm
    apply sum_nbij' (fun b => 2 * b) (fun n => n / 2)
    · intro b hb
      simp only [mem_filter, mem_range] at hb ⊢
      exact ⟨by omega, even_two_mul b⟩
    · intro n hn
      simp only [mem_filter, mem_range] at hn ⊢
      obtain ⟨hn1, hn2⟩ := hn
      have := Nat.two_mul_div_two_of_even hn2
      exact ⟨by omega, by omega⟩
    · intro b _
      simp
    · intro n hn
      simp only [mem_filter, mem_range] at hn
      exact Nat.two_mul_div_two_of_even hn.2
    · intro b _
      rfl
  · intro n _ hn
    by_contra hodd'
    rw [Nat.not_even_iff_odd] at hodd'
    obtain ⟨b, rfl⟩ := hodd'
    exact hn (hodd b)

/-- The inner `s`-sum of the left side of (**):
`G_M(y) = Σ_{b ≤ M/2} y^{2b} q^{2b²+b} / ((q²;q²)_b (q)_{M-2b})`. -/
noncomputable def GM (M : ℕ) (y : K⟦X⟧) : K⟦X⟧ :=
  ∑ b ∈ range (M + 1),
    if 2 * b ≤ M then X ^ tri (2 * b) * (y ^ (2 * b) * iq (X ^ 2) b * iq X (M - 2 * b)) else 0

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma Zf_zero (M : ℕ) (x : K⟦X⟧) : Zf 0 M x = GM M x := by
  simp only [Zf, GM]
  rw [sum_range_even]
  · apply sum_congr rfl
    intro b _
    rw [(cR_zero_even_odd b).1]
  · intro b
    rw [(cR_zero_even_odd b).2]; simp

/-- **Lemma 4, form (**).**  For all `R, M ≥ 0` and every power series `x`:
`Σ_{r ≤ R} Σ_{b ≤ M/2} [R;r] x^{r+2b} q^{(r+b)²+b²+b} (xq^{r+M+1};q)_∞ / ((q²;q²)_b (q)_{M-2b})
   = Σ_{L ≤ M} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_∞ / ((q)_L (q)_{M-L})`. -/
theorem star2 (M R : ℕ) (x : K⟦X⟧) :
    Wf (fun y => Einf M y * GM M y) R x = Uss R M x := by
  refine rec_unique (Wf (fun y => Einf M y * GM M y)) (fun R x => Uss R M x)
    (Wf_succ _) (fun R x => Uss_succ_R R M x) ?_ R x
  intro x
  rw [Wf_zero, Uss_eq, Zf_zero, zero_add]

/-! ### The form (*∞) -/

/-- `G_∞(y) = Σ_{b ≥ 0} y^{2b} q^{2b²+b} / (q²;q²)_b`. -/
noncomputable def Ginf (y : K⟦X⟧) : K⟦X⟧ :=
  ∑' b, X ^ tri (2 * b) * (y ^ (2 * b) * iq (X ^ 2) b)

/-- Right side of (*∞): `Σ_{L ≥ 0} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_∞ / (q)_L`. -/
noncomputable def Uinf (R : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑' L, X ^ tri L * (x ^ L * iq X L * Einf (R + L) x)

lemma Uinf_summable (R : ℕ) (x : K⟦X⟧) :
    Summable (fun L => X ^ tri L * (x ^ L * iq X L * Einf (R + L) x)) :=
  summable_nat (fun L => Dvd.dvd.mul_right (pow_dvd_pow X (le_tri L)) _)

lemma Uinf_add_mul (R : ℕ) (x y a : K⟦X⟧) :
    Uinf R x + a * Uinf R y = ∑' L, (X ^ tri L * (x ^ L * iq X L * Einf (R + L) x)
      + a * (X ^ tri L * (y ^ L * iq X L * Einf (R + L) y))) := by
  have hs2 : Summable (fun L => a * (X ^ tri L * (y ^ L * iq X L * Einf (R + L) y))) :=
    summable_nat (fun L => Dvd.dvd.mul_left (Dvd.dvd.mul_right (pow_dvd_pow X (le_tri L)) _) _)
  have e2 : a * Uinf R y = ∑' L, a * (X ^ tri L * (y ^ L * iq X L * Einf (R + L) y)) := by
    unfold Uinf
    exact ((Uinf_summable R y).tsum_mul_left a).symm
  rw [e2]
  unfold Uinf
  exact ((Uinf_summable R x).tsum_add hs2).symm

lemma Uinf_succ (R : ℕ) (x : K⟦X⟧) :
    Uinf (R + 1) x = Uinf R x + x * X ^ (R + 1) * Uinf R (x * X) := by
  rw [Uinf_add_mul]
  unfold Uinf
  refine tsum_congr (fun L => ?_)
  have e : R + 1 + L = R + L + 1 := by omega
  have hE := Einf_succ' (R + L) x
  rw [e, hE]
  ring

lemma Uinf_zero (x : K⟦X⟧) : Uinf 0 x = Ginf x := by
  set F : ℕ → ℕ → K⟦X⟧ := fun b M =>
    X ^ tri M * (x ^ M * ((-1) ^ (M - b) * iq X b * iq X (M - b))) with hF
  have e1 : Uinf 0 x = ∑' b, ∑' m, F b (m + b) := by
    simp only [Uinf]
    apply tsum_congr
    intro L
    have hm : X ^ tri L * (x ^ L * iq X L * Einf L x)
        = ∑' l, (X ^ tri L * (x ^ L * iq X L)) *
            (X ^ (tri l + L * l) * ((-1) ^ l * x ^ l * iq X l)) := by
      rw [(Einf_summable L x).tsum_mul_left (X ^ tri L * (x ^ L * iq X L)), Einf]; ring
    rw [zero_add, hm]
    apply tsum_congr
    intro l
    simp only [hF]
    rw [Nat.add_sub_cancel, tri_add]
    ring
  rw [e1, tsum_tsum_shift (fun b => b) (fun b => le_refl b) F
    (fun b M => Dvd.dvd.mul_right (pow_dvd_pow X (le_tri M)) _)]
  have e2 : ∀ M : ℕ, (∑ b ∈ range (M + 1), if b ≤ M then F b M else 0)
      = X ^ tri M * (x ^ M * cR 0 M) := by
    intro M
    have : ∀ b ∈ range (M + 1), (if b ≤ M then F b M else 0) = F b M := by
      intro b hb
      rw [if_pos (Nat.lt_succ_iff.mp (mem_range.mp hb))]
    rw [sum_congr rfl this]
    simp only [hF, cR]
    rw [← mul_sum, ← mul_sum, ← sum_range_reflect]
    congr 2
    apply sum_congr rfl
    intro j hj
    have hj' : j ≤ M := Nat.lt_succ_iff.mp (mem_range.mp hj)
    have a1 : M + 1 - 1 - j = M - j := by omega
    have a2 : M - (M - j) = j := by omega
    rw [a1, a2]
    ring
  rw [tsum_congr e2]
  set g : ℕ → K⟦X⟧ := fun M => X ^ tri M * (x ^ M * cR 0 M) with hg
  have hinj : Function.Injective (fun b : ℕ => 2 * b) := fun a b h => by
    simp only at h; omega
  have hsupp : Function.support g ⊆ Set.range (fun b : ℕ => 2 * b) := by
    intro M hM
    rcases Nat.even_or_odd' M with ⟨b, hb | hb⟩
    · exact ⟨b, hb.symm⟩
    · exfalso
      apply hM
      simp only [hg]
      rw [hb, (cR_zero_even_odd b).2]; simp
  rw [← hinj.tsum_eq hsupp]
  simp only [Ginf, hg]
  apply tsum_congr
  intro b
  rw [(cR_zero_even_odd b).1]

/-- **Lemma 4, form (*∞).**  For all `R ≥ 0` and every power series `x`:
`Σ_{r ≤ R} Σ_{b ≥ 0} [R;r] x^{r+2b} q^{(r+b)²+b²+b} / (q²;q²)_b
   = Σ_{L ≥ 0} x^L q^{L(L+1)/2} (xq^{R+L+1};q)_∞ / (q)_L`. -/
theorem starInf (R : ℕ) (x : K⟦X⟧) : Wf Ginf R x = Uinf R x := by
  refine rec_unique (Wf Ginf) Uinf (Wf_succ _) Uinf_succ ?_ R x
  intro x
  rw [Wf_zero, Uinf_zero]

end KRTW
