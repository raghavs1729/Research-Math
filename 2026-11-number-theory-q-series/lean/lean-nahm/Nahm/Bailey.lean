/-
Nahm/Bailey.lean — Bailey's lemma in the limit `ρ₁, ρ₂ → ∞`, relative to an arbitrary `A`, over an
arbitrary base `Q` (a power series with zero constant term), cleared of the denominators
`(AQ;Q)_{2n}`; the unit Bailey pair; the Bailey chain started from an arbitrary Bailey pair.

A pair of sequences `(α, β)` is a Bailey pair relative to `A` if
`β_n = Σ_{r ≤ n} α_r / ((Q;Q)_{n-r} (AQ;Q)_{n+r})`; multiplied by `(AQ;Q)_{2n}` this reads
`(AQ;Q)_{2n} β_n = Σ_{r ≤ n} α_r (AQ^{n+r+1};Q)_{n-r} / (Q;Q)_{n-r}` (`IsBP`).
Bailey's lemma (`ρ₁, ρ₂ → ∞`): then `(A^r Q^{r²} α_r, Σ_k A^k Q^{k²} β_k / (Q;Q)_{n-k})` is again a
Bailey pair relative to `A` (`IsBP.step`).
-/
import KRTW.Basic
import KRTW.Bailey

namespace Nahm

open PowerSeries Finset KRTW

variable {K : Type*} [Field K] (Q : K⟦X⟧)

/-- `(y;Q)_n = ∏_{i<n} (1 - y Q^i)`. -/
noncomputable def gq (y : K⟦X⟧) (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - y * Q ^ i)

@[simp] lemma gq_zero (y : K⟦X⟧) : gq Q y 0 = 1 := by simp [gq]

lemma gq_succ (y : K⟦X⟧) (n : ℕ) : gq Q y (n + 1) = gq Q y n * (1 - y * Q ^ n) := by
  simp [gq, prod_range_succ]

lemma gq_succ' (y : K⟦X⟧) (n : ℕ) : gq Q y (n + 1) = (1 - y) * gq Q (y * Q) n := by
  simp only [gq]
  rw [prod_range_succ', pow_zero, mul_one, mul_comm]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  rw [pow_succ]; ring

lemma gq_add (y : K⟦X⟧) (m n : ℕ) : gq Q y (m + n) = gq Q y m * gq Q (y * Q ^ m) n := by
  simp only [gq]
  rw [prod_range_add]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  rw [pow_add]; ring

lemma gq_split (A : K⟦X⟧) (r i N : ℕ) (hi : i ≤ N) :
    gq Q (A * Q ^ (r + i + r + 1)) i * gq Q (A * Q ^ (2 * (r + i) + 1)) (2 * (r + N) - 2 * (r + i))
      = gq Q (A * Q ^ (2 * r + i + 1)) (N - i) * gq Q (A * Q ^ (r + N + r + 1)) N := by
  have a1 : r + i + r + 1 = 2 * r + i + 1 := by omega
  have a2 : A * Q ^ (2 * (r + i) + 1) = A * Q ^ (2 * r + i + 1) * Q ^ i := by
    rw [mul_assoc, ← pow_add]; congr 2; omega
  have a3 : A * Q ^ (r + N + r + 1) = A * Q ^ (2 * r + i + 1) * Q ^ (N - i) := by
    rw [mul_assoc, ← pow_add]; congr 2; omega
  have a4 : i + (2 * (r + N) - 2 * (r + i)) = (N - i) + N := by omega
  rw [a1, a2, a3, ← gq_add, ← gq_add, a4]

/-! ### The unit Bailey pair relative to `A` -/

/-- `α_r = (-1)^r Q^{r(r-1)/2} (A;Q)_r (1 - A Q^{2r}) / ((Q;Q)_r (1 - A))`, written without
division: `α_0 = 1`, `α_{r+1} = (-1)^{r+1} Q^{r(r+1)/2} (AQ;Q)_r (1 - A Q^{2r+2}) / (Q;Q)_{r+1}`. -/
noncomputable def ua (A : K⟦X⟧) : ℕ → K⟦X⟧
  | 0 => 1
  | r + 1 => (-1) ^ (r + 1) * Q ^ tri r * gq Q (A * Q) r * (1 - A * Q ^ (2 * r + 2)) * iq Q (r + 1)

/-- `Σ_{r ≤ n} α_r (AQ^{n+r+1};Q)_{n-r} / (Q;Q)_{n-r}` (`= (AQ;Q)_{2n} β_n` for a Bailey pair). -/
noncomputable def bps (A : K⟦X⟧) (α : ℕ → K⟦X⟧) (n : ℕ) : K⟦X⟧ :=
  ∑ r ∈ range (n + 1), α r * iq Q (n - r) * gq Q (A * Q ^ (n + r + 1)) (n - r)

/-- `(α, β)` is a Bailey pair relative to `A` (cleared of the denominator `(AQ;Q)_{2n}`). -/
def IsBP (A : K⟦X⟧) (α β : ℕ → K⟦X⟧) : Prop :=
  ∀ n, gq Q (A * Q) (2 * n) * β n = bps Q A α n

variable [Nil Q]

/-- **The unit Bailey pair**: `Σ_{r ≤ n} α_r / ((Q)_{n-r} (AQ;Q)_{n+r}) = δ_{n,0}` (cleared). -/
theorem ubp (A : K⟦X⟧) (n : ℕ) : bps Q A (ua Q A) n = if n = 0 then 1 else 0 := by
  cases n with
  | zero => simp [bps, ua]
  | succ N =>
    rw [if_neg (Nat.succ_ne_zero N)]
    apply (mul_eq_zero.mp _).resolve_left (one_sub_pow_ne_zero Q N)
    set w : ℕ → K⟦X⟧ := fun r =>
      (-1) ^ (r + 1) * Q ^ tri r * gq Q (A * Q) r * gq Q (A * Q ^ (N + r + 2)) (N + 1 - r) * iq Q r *
        (iq Q (N + 1 - r) * (1 - Q ^ (N + 1 - r))) with hw
    have hw0 : w 0 = -(iq Q N * gq Q (A * Q ^ (N + 2)) (N + 1)) := by
      simp only [hw]
      have : iq Q (N + 1 - 0) * (1 - Q ^ (N + 1 - 0)) = iq Q N := iq_succ_mul Q N
      rw [this]
      simp
      ring
    have hwN : w (N + 1) = 0 := by simp [hw]
    have key : ∀ r ∈ range (N + 1),
        (1 - Q ^ (N + 1)) * (ua Q A (r + 1) * iq Q (N + 1 - (r + 1)) *
          gq Q (A * Q ^ (N + 1 + (r + 1) + 1)) (N + 1 - (r + 1))) = w r - w (r + 1) := by
      intro r hr
      have hr' : r ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hr)
      have e1 : N + 1 - r = (N - r) + 1 := by omega
      have e2 : Q ^ (N + 1) = Q ^ (r + 1) * Q ^ (N - r) := by
        rw [← pow_add]; congr 1; omega
      have e3 : iq Q r = iq Q (r + 1) * (1 - Q ^ (r + 1)) := (iq_succ_mul Q r).symm
      have e4 : A * Q ^ (N + r + 2) * Q = A * Q ^ (N + 1 + (r + 1) + 1) := by
        rw [mul_assoc, ← pow_succ]; congr 2; omega
      have e5 : A * Q * Q ^ r = A * Q ^ (r + 1) := by rw [mul_assoc, ← pow_succ']
      have e6 : A * Q ^ (N + (r + 1) + 2) = A * Q ^ (N + 1 + (r + 1) + 1) := by
        congr 2; omega
      simp only [hw, ua]
      rw [Nat.add_sub_add_right, e1, gq_succ' Q (A * Q ^ (N + r + 2)) (N - r), e4,
        iq_succ_mul Q (N - r), e3, gq_succ Q (A * Q) r, e5, e6, tri_succ, e2]
      have e7 : Q ^ (2 * r + 2) = Q ^ (r + 1) * Q ^ (r + 1) := by
        rw [← pow_add]; congr 1; omega
      have e8 : Q ^ (N + r + 2) = Q ^ (r + 1) * Q ^ (N - r) * Q ^ (r + 1) := by
        rw [← pow_add, ← pow_add]; congr 1; omega
      rw [e7, e8]
      ring
    simp only [bps]
    rw [sum_range_succ', mul_add, mul_sum, sum_congr rfl key, sum_range_sub', hwN, hw0]
    have e0 : iq Q N = (1 - Q ^ (N + 1)) * iq Q (N + 1) := by
      rw [mul_comm, iq_succ_mul]
    simp only [ua, Nat.sub_zero, add_zero, one_mul]
    have e9 : N + 1 + 1 = N + 2 := rfl
    rw [e9, e0]
    ring

/-! ### The inner sum of Bailey's lemma -/

/-- The `q`-Pascal rule `[n+1;k+1] = Q^{k+1} [n;k+1] + [n;k]`. -/
lemma qbinom_succ_succ' (n k : ℕ) :
    qbinom Q (n + 1) (k + 1) = Q ^ (k + 1) * qbinom Q n (k + 1) + qbinom Q n k := by
  rcases Nat.lt_or_ge n k with h | h
  · rw [qbinom_of_lt Q (show n + 1 < k + 1 by omega), qbinom_of_lt Q (show n < k + 1 by omega),
      qbinom_of_lt Q h]; simp
  rcases Nat.eq_or_lt_of_le h with h' | h'
  · subst h'
    rw [qbinom_self, qbinom_self, qbinom_of_lt Q (Nat.lt_succ_self k)]; simp
  · have e1 : qbinom Q (n + 1) (k + 1) = poch Q (n + 1) * iq Q (k + 1) * iq Q (n - k) := by
      simp [qbinom, h]
    have e2 : qbinom Q n (k + 1) = poch Q n * iq Q (k + 1) * iq Q (n - (k + 1)) := by
      simp [qbinom, Nat.succ_le_of_lt h']
    have e3 : qbinom Q n k = poch Q n * iq Q k * iq Q (n - k) := by
      simp [qbinom, h]
    have e4 : iq Q (n - (k + 1)) = iq Q (n - k) * (1 - Q ^ (n - k)) := by
      have : n - k = (n - (k + 1)) + 1 := by omega
      rw [this, iq_succ_mul]
    have e5 : iq Q k = iq Q (k + 1) * (1 - Q ^ (k + 1)) := (iq_succ_mul Q k).symm
    have e6 : Q ^ (n + 1) = Q ^ (k + 1) * Q ^ (n - k) := by
      rw [← pow_add]; congr 1; omega
    rw [e1, e2, e3, e4, e5, poch_succ, e6]
    ring

/-- `Σ_i [N;i] B^i Q^{i²} (BQ^{i+1};Q)_{N-i} = 1` (the limit `b → ∞`, `d → 0`, `bd` fixed, of
`q`-Chu–Vandermonde). -/
theorem inner_one : ∀ (N : ℕ) (B : K⟦X⟧),
    ∑ i ∈ range (N + 1), qbinom Q N i * (B ^ i * Q ^ (i ^ 2) * gq Q (B * Q ^ (i + 1)) (N - i)) = 1
  | 0, B => by simp
  | N + 1, B => by
    have ih := inner_one N (B * Q)
    set c : ℕ → ℕ → K⟦X⟧ → K⟦X⟧ := fun N i B =>
      qbinom Q N i * (B ^ i * Q ^ (i ^ 2) * gq Q (B * Q ^ (i + 1)) (N - i)) with hc
    set D : ℕ → K⟦X⟧ := fun i =>
      qbinom Q N i * (B ^ (i + 1) * Q ^ (i ^ 2 + 2 * i + 1) * gq Q (B * Q ^ (i + 2)) (N - i)) with hD
    have h0 : c (N + 1) 0 B = c N 0 (B * Q) - D 0 := by
      simp only [hc, hD, qbinom_zero_right, pow_zero, one_mul, zero_add, Nat.sub_zero, pow_one,
        mul_one, zero_pow (two_ne_zero)]
      rw [gq_succ']
      have : B * Q * Q = B * Q ^ 2 := by ring
      rw [this]
      ring
    have hs : ∀ i, c (N + 1) (i + 1) B = c N (i + 1) (B * Q) - D (i + 1) + D i := by
      intro i
      simp only [hc, hD]
      rw [qbinom_succ_succ', Nat.add_sub_add_right]
      rcases Nat.lt_or_ge i N with h | h
      · obtain ⟨t, rfl⟩ : ∃ t, N = i + 1 + t := ⟨N - i - 1, by omega⟩
        rw [show i + 1 + t - i = t + 1 by omega, show i + 1 + t - (i + 1) = t by omega,
          gq_succ']
        have a1 : B * Q ^ (i + 1 + 1) * Q = B * Q * Q ^ (i + 1 + 1) := by ring
        have a2 : B * Q ^ (i + 1 + 2) = B * Q * Q ^ (i + 1 + 1) := by ring
        rw [a1, a2, show i + 2 = i + 1 + 1 by ring]
        ring
      · rw [qbinom_of_lt Q (show N < i + 1 by omega)]
        ring
    show ∑ i ∈ range (N + 1 + 1), c (N + 1) i B = 1
    rw [sum_range_succ', h0, sum_congr rfl (fun i _ => hs i)]
    have tel : ∑ i ∈ range (N + 1), (c N (i + 1) (B * Q) - D (i + 1) + D i)
        = ∑ i ∈ range (N + 1), c N (i + 1) (B * Q) - (D (N + 1) - D 0) := by
      have := sum_range_sub (fun i => D i) (N + 1)
      rw [sum_add_distrib, sum_sub_distrib, ← this, sum_sub_distrib]
      ring
    have hDN : D (N + 1) = 0 := by simp [hD, qbinom_of_lt Q (Nat.lt_succ_self N)]
    have hcN : c N (N + 1) (B * Q) = 0 := by simp [hc, qbinom_of_lt Q (Nat.lt_succ_self N)]
    rw [tel, hDN]
    have : ∑ i ∈ range (N + 1 + 1), c N i (B * Q) = 1 := by
      rw [sum_range_succ, hcN, add_zero]; exact ih
    rw [sum_range_succ'] at this
    linear_combination this

/-- The inner sum in the form used: `Σ_i B^i Q^{i²} (BQ^{i+1};Q)_{N-i} / ((Q)_i (Q)_{N-i}) = 1/(Q)_N`. -/
theorem inner_iq (N : ℕ) (B : K⟦X⟧) :
    ∑ i ∈ range (N + 1), B ^ i * Q ^ (i ^ 2) * iq Q i * iq Q (N - i) * gq Q (B * Q ^ (i + 1)) (N - i)
      = iq Q N := by
  have h := inner_one Q N B
  have h1 : ∑ i ∈ range (N + 1), qbinom Q N i * (B ^ i * Q ^ (i ^ 2) * gq Q (B * Q ^ (i + 1)) (N - i))
      = poch Q N * ∑ i ∈ range (N + 1),
          B ^ i * Q ^ (i ^ 2) * iq Q i * iq Q (N - i) * gq Q (B * Q ^ (i + 1)) (N - i) := by
    rw [mul_sum]
    refine sum_congr rfl (fun i hi => ?_)
    have hi' : i ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hi)
    simp only [qbinom, hi', if_true]
    ring
  rw [h1] at h
  have h2 := iq_mul_poch Q N
  linear_combination (iq Q N) * h - (∑ i ∈ range (N + 1),
    B ^ i * Q ^ (i ^ 2) * iq Q i * iq Q (N - i) * gq Q (B * Q ^ (i + 1)) (N - i)) * h2

/-! ### Bailey's lemma and the chain -/

/-- **Bailey's lemma** (`ρ₁, ρ₂ → ∞`), cleared of denominators: if `Δ_k = Σ_{r≤k} α_r
(AQ^{k+r+1};Q)_{k-r}/(Q)_{k-r}` for all `k`, then
`Σ_{k ≤ n} A^k Q^{k²} Δ_k (AQ^{2k+1};Q)_{2n-2k}/(Q)_{n-k} = Σ_{r ≤ n} A^r Q^{r²} α_r (AQ^{n+r+1};Q)_{n-r}/(Q)_{n-r}`. -/
theorem bailey_step (A : K⟦X⟧) (α Δ : ℕ → K⟦X⟧) (hΔ : ∀ k, Δ k = bps Q A α k) (n : ℕ) :
    ∑ k ∈ range (n + 1), A ^ k * Q ^ (k ^ 2) * iq Q (n - k) * Δ k * gq Q (A * Q ^ (2 * k + 1)) (2 * n - 2 * k)
      = bps Q A (fun r => A ^ r * Q ^ (r ^ 2) * α r) n := by
  set F : ℕ → ℕ → K⟦X⟧ := fun k r => A ^ k * Q ^ (k ^ 2) * iq Q (n - k) *
    (α r * iq Q (k - r) * gq Q (A * Q ^ (k + r + 1)) (k - r)) *
      gq Q (A * Q ^ (2 * k + 1)) (2 * n - 2 * k) with hF
  have h1 : ∑ k ∈ range (n + 1), A ^ k * Q ^ (k ^ 2) * iq Q (n - k) * Δ k *
      gq Q (A * Q ^ (2 * k + 1)) (2 * n - 2 * k)
      = ∑ k ∈ range (n + 1), ∑ r ∈ range (k + 1), F k r := by
    refine sum_congr rfl (fun k _ => ?_)
    rw [hΔ, bps]
    simp only [mul_sum, sum_mul, hF]
  rw [h1, sum_range_triangle, bps]
  refine sum_congr rfl (fun r hr => ?_)
  have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
  obtain ⟨N, rfl⟩ := Nat.exists_eq_add_of_le hr'
  rw [Nat.add_sub_cancel_left]
  have hterm : ∀ i ∈ range (N + 1), F (r + i) r
      = (α r * gq Q (A * Q ^ (r + N + r + 1)) N) * (A ^ r * Q ^ (r ^ 2)) *
        ((A * Q ^ (2 * r)) ^ i * Q ^ (i ^ 2) * iq Q i * iq Q (N - i) *
          gq Q (A * Q ^ (2 * r) * Q ^ (i + 1)) (N - i)) := by
    intro i hi
    have hi' : i ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hi)
    have hs := gq_split Q A r i N hi'
    have e1 : A * Q ^ (2 * r) * Q ^ (i + 1) = A * Q ^ (2 * r + i + 1) := by ring
    have e2 : A ^ (r + i) * Q ^ ((r + i) ^ 2) = A ^ r * Q ^ (r ^ 2) * ((A * Q ^ (2 * r)) ^ i * Q ^ (i ^ 2)) := by
      ring
    simp only [hF]
    rw [Nat.add_sub_add_left, Nat.add_sub_cancel_left, e1, e2]
    linear_combination (A ^ r * Q ^ (r ^ 2) * ((A * Q ^ (2 * r)) ^ i * Q ^ (i ^ 2)) * iq Q (N - i) *
      α r * iq Q i) * hs
  rw [sum_congr rfl hterm, ← mul_sum, inner_iq]
  ring

/-- One step of the Bailey chain on `β`: `β'_n = Σ_{k ≤ n} A^k Q^{k²} β_k / (Q)_{n-k}`. -/
noncomputable def bstep (A : K⟦X⟧) (β : ℕ → K⟦X⟧) (n : ℕ) : K⟦X⟧ :=
  ∑ k ∈ range (n + 1), A ^ k * Q ^ (k ^ 2) * iq Q (n - k) * β k

/-- Bailey's lemma for Bailey pairs. -/
theorem IsBP.step {A : K⟦X⟧} {α β : ℕ → K⟦X⟧} (h : IsBP Q A α β) :
    IsBP Q A (fun r => A ^ r * Q ^ (r ^ 2) * α r) (bstep Q A β) := by
  intro n
  rw [← bailey_step Q A α (fun k => gq Q (A * Q) (2 * k) * β k) (fun k => h k) n, bstep, mul_sum]
  refine sum_congr rfl (fun k hk => ?_)
  have hk' : k ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hk)
  have e : gq Q (A * Q) (2 * n) = gq Q (A * Q) (2 * k) * gq Q (A * Q ^ (2 * k + 1)) (2 * n - 2 * k) := by
    have a1 : A * Q ^ (2 * k + 1) = A * Q * Q ^ (2 * k) := by ring
    rw [a1, ← gq_add]; congr 1; omega
  rw [e]; ring

/-- The unit Bailey pair `(α, δ)`. -/
theorem isBP_unit (A : K⟦X⟧) : IsBP Q A (ua Q A) (fun n => if n = 0 then 1 else 0) := by
  intro n
  rw [ubp]
  split_ifs with h
  · subst h; simp
  · simp [h]

omit [Nil Q] in
lemma bstep_delta (A : K⟦X⟧) : bstep Q A (fun n => if n = 0 then 1 else 0) = iq Q := by
  funext n
  rw [bstep, sum_range_succ', sum_eq_zero (fun k _ => by simp)]
  simp

/-- The first step of the chain: `(A^r Q^{r²} α_r, 1/(Q)_n)` is a Bailey pair relative to `A`. -/
theorem isBP_iq (A : K⟦X⟧) : IsBP Q A (fun r => A ^ r * Q ^ (r ^ 2) * ua Q A r) (iq Q) := by
  have h := (isBP_unit Q A).step
  rwa [bstep_delta] at h

end Nahm
