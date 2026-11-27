/-
KRTW/Bailey.lean — Bailey's lemma (the limiting case (S2) of `proofs/KRTW-general.md` §5) and the
unit Bailey pair, relative to `a = x²q`, in a form cleared of the denominators `(aq;q)_m`,
`(-xq;q)_m`; the Bailey chain.  All sums are finite; `x : K⟦X⟧` is arbitrary, `q = X`.
-/
import KRTW.Basic

namespace KRTW

open PowerSeries Finset

variable {K : Type*} [Field K]

/-- `(y;q)_n = ∏_{i<n} (1 - y q^i)` for an arbitrary power series `y`. -/
noncomputable def gp (y : K⟦X⟧) (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - y * X ^ i)

@[simp] lemma gp_zero (y : K⟦X⟧) : gp y 0 = 1 := by simp [gp]

lemma gp_succ (y : K⟦X⟧) (n : ℕ) : gp y (n + 1) = gp y n * (1 - y * X ^ n) := by
  simp [gp, prod_range_succ]

lemma gp_succ' (y : K⟦X⟧) (n : ℕ) : gp y (n + 1) = (1 - y) * gp (y * X) n := by
  simp only [gp]
  rw [prod_range_succ', pow_zero, mul_one, mul_comm]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  rw [pow_succ]; ring

lemma gp_add (y : K⟦X⟧) (m n : ℕ) : gp y (m + n) = gp y m * gp (y * X ^ m) n := by
  simp only [gp]
  rw [prod_range_add]
  congr 1
  refine prod_congr rfl (fun i _ => ?_)
  rw [pow_add]; ring

lemma tri_pred (i : ℕ) : tri i = tri (i - 1) + i := by
  cases i with
  | zero => simp
  | succ i => simp [tri_succ]

/-- **q-Chu–Vandermonde** (polynomial form):
`Σ_i [N;i] (-1)^i q^{i(i-1)/2} d^i (b;q)_i (bdq^i;q)_{N-i} = (d;q)_N`. -/
theorem cv (d : K⟦X⟧) : ∀ (N : ℕ) (b : K⟦X⟧),
    ∑ i ∈ range (N + 1), qbinom X N i *
      ((-1) ^ i * X ^ tri (i - 1) * d ^ i * gp b i * gp (b * d * X ^ i) (N - i)) = gp d N
  | 0, b => by simp
  | N + 1, b => by
    have ih1 := cv d N b
    have ih2 := cv d N (b * X)
    set T : ℕ → K⟦X⟧ := fun i =>
      (-1) ^ i * X ^ tri (i - 1) * d ^ i * gp b i * gp (b * d * X ^ i) (N + 1 - i) with hT
    have hA : ∑ i ∈ range (N + 1), qbinom X N i * T i = (1 - b * d * X ^ N) * gp d N := by
      rw [← ih1, mul_sum]
      refine sum_congr rfl (fun i hi => ?_)
      have hi' : i ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hi)
      have e : N + 1 - i = (N - i) + 1 := by omega
      have e2 : b * d * X ^ i * X ^ (N - i) = b * d * X ^ N := by
        rw [mul_assoc, ← pow_add]; congr 2; omega
      simp only [hT]
      rw [e, gp_succ, e2]; ring
    have hB : ∑ i ∈ range (N + 1), X ^ (N - i) * qbinom X N i * T (i + 1)
        = -(d * X ^ N * (1 - b)) * gp d N := by
      rw [← ih2, mul_sum]
      refine sum_congr rfl (fun i hi => ?_)
      have hi' : i ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hi)
      have e : (X : K⟦X⟧) ^ (N - i) * X ^ i = X ^ N := by
        rw [← pow_add]; congr 1; omega
      have e3 : b * d * X ^ (i + 1) = b * X * d * X ^ i := by ring
      simp only [hT]
      rw [Nat.add_sub_cancel, Nat.add_sub_add_right, gp_succ' b i, e3, tri_pred i]
      linear_combination (-(qbinom X N i * (-1) ^ i * X ^ tri (i - 1) * d ^ (i + 1) * (1 - b) *
        gp (b * X) i * gp (b * X * d * X ^ i) (N - i))) * e
    have h1 : ∑ i ∈ range (N + 1 + 1), qbinom X (N + 1) i * T i
        = ∑ i ∈ range (N + 1), qbinom X (N + 1) (i + 1) * T (i + 1) + T 0 := by
      rw [sum_range_succ']; simp
    have h2 : ∑ i ∈ range (N + 1), qbinom X N i * T i
        = ∑ i ∈ range (N + 1), qbinom X N (i + 1) * T (i + 1) + T 0 := by
      have : ∑ i ∈ range (N + 1), qbinom X N i * T i
          = ∑ i ∈ range (N + 1 + 1), qbinom X N i * T i := by
        rw [sum_range_succ _ (N + 1), qbinom_of_lt X (Nat.lt_succ_self N)]; simp
      rw [this, sum_range_succ']; simp
    have h3 : ∑ i ∈ range (N + 1 + 1), qbinom X (N + 1) i * T i
        = ∑ i ∈ range (N + 1), qbinom X N i * T i
          + ∑ i ∈ range (N + 1), X ^ (N - i) * qbinom X N i * T (i + 1) := by
      rw [h1, h2, add_right_comm, ← sum_add_distrib]
      congr 1
      refine sum_congr rfl (fun i _ => ?_)
      rw [qbinom_succ_succ]; ring
    rw [h3, hA, hB, gp_succ]
    ring

/-! ### The unit Bailey pair relative to `A` -/

/-- `α_r = (-1)^r q^{r(r-1)/2} (A;q)_r (1 - A q^{2r}) / ((q)_r (1 - A))`, written without division:
`α_0 = 1`, `α_{r+1} = (-1)^{r+1} q^{r(r+1)/2} (Aq;q)_r (1 - A q^{2r+2}) / (q)_{r+1}`. -/
noncomputable def alpha (A : K⟦X⟧) : ℕ → K⟦X⟧
  | 0 => 1
  | r + 1 => (-1) ^ (r + 1) * X ^ tri r * gp (A * X) r * (1 - A * X ^ (2 * r + 2)) * iq X (r + 1)

/-- The right side of the Bailey pair relation, multiplied by `(Aq;q)_{2n}`:
`Σ_{r ≤ n} α_r (Aq^{n+r+1};q)_{n-r} / (q)_{n-r}`  (`= β_n (Aq;q)_{2n}`). -/
noncomputable def BPsum (A : K⟦X⟧) (α : ℕ → K⟦X⟧) (n : ℕ) : K⟦X⟧ :=
  ∑ r ∈ range (n + 1), α r * iq X (n - r) * gp (A * X ^ (n + r + 1)) (n - r)

/-- **The unit Bailey pair**: `Σ_{r ≤ n} α_r / ((q)_{n-r} (Aq;q)_{n+r}) = δ_{n,0}` (cleared of the
denominator `(Aq;q)_{2n}`). -/
theorem ubp (A : K⟦X⟧) (n : ℕ) : BPsum A (alpha A) n = if n = 0 then 1 else 0 := by
  cases n with
  | zero => simp [BPsum, alpha]
  | succ N =>
    rw [if_neg (Nat.succ_ne_zero N)]
    apply (mul_eq_zero.mp _).resolve_left (one_sub_pow_ne_zero (X : K⟦X⟧) N)
    set w : ℕ → K⟦X⟧ := fun r =>
      (-1) ^ (r + 1) * X ^ tri r * gp (A * X) r * gp (A * X ^ (N + r + 2)) (N + 1 - r) * iq X r *
        (iq X (N + 1 - r) * (1 - X ^ (N + 1 - r))) with hw
    have hw0 : w 0 = -(iq X N * gp (A * X ^ (N + 2)) (N + 1)) := by
      simp only [hw]
      have : iq (X : K⟦X⟧) (N + 1 - 0) * (1 - X ^ (N + 1 - 0)) = iq X N := iq_succ_mul X N
      rw [this]
      simp
      ring
    have hwN : w (N + 1) = 0 := by simp [hw]
    have key : ∀ r ∈ range (N + 1),
        (1 - X ^ (N + 1)) * (alpha A (r + 1) * iq X (N + 1 - (r + 1)) *
          gp (A * X ^ (N + 1 + (r + 1) + 1)) (N + 1 - (r + 1))) = w r - w (r + 1) := by
      intro r hr
      have hr' : r ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hr)
      have e1 : N + 1 - r = (N - r) + 1 := by omega
      have e2 : (X : K⟦X⟧) ^ (N + 1) = X ^ (r + 1) * X ^ (N - r) := by
        rw [← pow_add]; congr 1; omega
      have e3 : iq (X : K⟦X⟧) r = iq X (r + 1) * (1 - X ^ (r + 1)) := (iq_succ_mul X r).symm
      have e4 : A * X ^ (N + r + 2) * X = A * X ^ (N + 1 + (r + 1) + 1) := by
        rw [mul_assoc, ← pow_succ]; congr 2; omega
      have e5 : A * X * X ^ r = A * X ^ (r + 1) := by rw [mul_assoc, ← pow_succ']
      have e6 : A * X ^ (N + (r + 1) + 2) = A * X ^ (N + 1 + (r + 1) + 1) := by
        congr 2; omega
      simp only [hw, alpha]
      rw [Nat.add_sub_add_right, e1, gp_succ' (A * X ^ (N + r + 2)) (N - r), e4,
        iq_succ_mul X (N - r), e3, gp_succ (A * X) r, e5, e6, tri_succ, e2]
      have e7 : (X : K⟦X⟧) ^ (2 * r + 2) = X ^ (r + 1) * X ^ (r + 1) := by
        rw [← pow_add]; congr 1; omega
      have e8 : (X : K⟦X⟧) ^ (N + r + 2) = X ^ (r + 1) * X ^ (N - r) * X ^ (r + 1) := by
        rw [← pow_add, ← pow_add]; congr 1; omega
      rw [e7, e8]
      ring
    simp only [BPsum]
    rw [sum_range_succ', mul_add, mul_sum, sum_congr rfl key, sum_range_sub', hwN, hw0]
    have e0 : iq (X : K⟦X⟧) N = (1 - X ^ (N + 1)) * iq X (N + 1) := by
      rw [mul_comm, iq_succ_mul]
    simp only [alpha, Nat.sub_zero, add_zero, one_mul]
    have e9 : N + 1 + 1 = N + 2 := rfl
    rw [e9, e0]
    ring

/-! ### Bailey's lemma (S2) -/

lemma sum_range_triangle {M : Type*} [AddCommMonoid M] (f : ℕ → ℕ → M) (n : ℕ) :
    ∑ k ∈ range (n + 1), ∑ r ∈ range (k + 1), f k r
      = ∑ r ∈ range (n + 1), ∑ i ∈ range (n - r + 1), f (r + i) r := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ, ih, sum_range_succ (fun r => f (n + 1) r) (n + 1)]
    conv_rhs => rw [sum_range_succ]
    rw [← add_assoc, ← sum_add_distrib]
    congr 1
    · refine sum_congr rfl (fun r hr => ?_)
      have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
      have e : n + 1 - r + 1 = (n - r + 1) + 1 := by omega
      have e2 : r + (n - r + 1) = n + 1 := by omega
      rw [e]
      conv_rhs => rw [sum_range_succ, e2]
    · simp

/-- q-Chu–Vandermonde in the form used in Bailey's lemma (`b = d = -yq`). -/
lemma cv' (y : K⟦X⟧) (N : ℕ) :
    ∑ i ∈ range (N + 1), iq X i * iq X (N - i) *
      (X ^ tri i * y ^ i * gp (-(y * X)) i * gp (y ^ 2 * X ^ (i + 2)) (N - i))
      = iq X N * gp (-(y * X)) N := by
  rw [← cv (-(y * X)) N (-(y * X)), mul_sum]
  refine sum_congr rfl (fun i hi => ?_)
  have hi' : i ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hi)
  have e : -(y * X) * -(y * X) * X ^ i = y ^ 2 * X ^ (i + 2) := by ring
  have h1 : iq (X : K⟦X⟧) N * poch X N = 1 := iq_mul_poch X N
  have hm : (-1 : K⟦X⟧) ^ i * (-1) ^ i = 1 := by rw [← mul_pow]; simp
  simp only [qbinom, hi', if_true]
  rw [e, neg_pow (y * X) i, tri_pred i]
  linear_combination (-(iq X i * iq X (N - i) * X ^ tri (i - 1) * (y * X) ^ i * gp (-(y * X)) i *
      gp (y ^ 2 * X ^ (i + 2)) (N - i) * (-1) ^ i * (-1) ^ i)) * h1
    + (-(iq X i * iq X (N - i) * X ^ tri (i - 1) * (y * X) ^ i * gp (-(y * X)) i *
      gp (y ^ 2 * X ^ (i + 2)) (N - i))) * hm

lemma bailey_inner (x : K⟦X⟧) (r N : ℕ) :
    ∑ i ∈ range (N + 1), X ^ tri (r + i) * (x ^ (r + i) * iq X (N - i) * iq X i *
      gp (-(x * X)) (r + i) * gp (x ^ 2 * X * X ^ (2 * r + i + 1)) (N - i))
      = X ^ tri r * x ^ r * iq X N * gp (-(x * X)) (r + N) := by
  have h := cv' (x * X ^ r) N
  have e1 : -(x * X) * X ^ r = -(x * X ^ r * X) := by ring
  rw [gp_add (-(x * X)) r N, e1]
  have e3 : X ^ tri r * x ^ r * iq X N * (gp (-(x * X)) r * gp (-(x * X ^ r * X)) N)
      = (X ^ tri r * x ^ r * gp (-(x * X)) r) * (iq X N * gp (-(x * X ^ r * X)) N) := by ring
  rw [e3, ← h, mul_sum]
  refine sum_congr rfl (fun i _ => ?_)
  have e2 : x ^ 2 * X * X ^ (2 * r + i + 1) = (x * X ^ r) ^ 2 * X ^ (i + 2) := by ring
  rw [tri_add, gp_add (-(x * X)) r i, e1, e2]
  ring

lemma gp_split (A : K⟦X⟧) (r i N : ℕ) (hi : i ≤ N) :
    gp (A * X ^ (r + i + r + 1)) i * gp (A * X ^ (2 * (r + i) + 1)) (2 * (r + N) - 2 * (r + i))
      = gp (A * X ^ (2 * r + i + 1)) (N - i) * gp (A * X ^ (r + N + r + 1)) N := by
  have a1 : r + i + r + 1 = 2 * r + i + 1 := by omega
  have a2 : A * X ^ (2 * (r + i) + 1) = A * X ^ (2 * r + i + 1) * X ^ i := by
    rw [mul_assoc, ← pow_add]; congr 2; omega
  have a3 : A * X ^ (r + N + r + 1) = A * X ^ (2 * r + i + 1) * X ^ (N - i) := by
    rw [mul_assoc, ← pow_add]; congr 2; omega
  have a4 : i + (2 * (r + N) - 2 * (r + i)) = (N - i) + N := by omega
  rw [a1, a2, a3, ← gp_add, ← gp_add, a4]

/-- **Bailey's lemma, case (S2)**, relative to `a = x²q`, cleared of denominators.  If
`Δ_k = (-xq;q)_k · Σ_{r ≤ k} α_r (aq^{k+r+1};q)_{k-r}/(q)_{k-r}` for all `k`, then
`Σ_{k ≤ n} x^k q^{k(k+1)/2} Δ_k (aq^{2k+1};q)_{2n-2k}/(q)_{n-k}
   = (-xq;q)_n · Σ_{r ≤ n} x^r q^{r(r+1)/2} α_r (aq^{n+r+1};q)_{n-r}/(q)_{n-r}`. -/
theorem bailey_step (x : K⟦X⟧) (α Δ : ℕ → K⟦X⟧)
    (hΔ : ∀ k, Δ k = gp (-(x * X)) k * BPsum (x ^ 2 * X) α k) (n : ℕ) :
    ∑ k ∈ range (n + 1), X ^ tri k * (x ^ k * iq X (n - k) * Δ k *
      gp (x ^ 2 * X * X ^ (2 * k + 1)) (2 * n - 2 * k))
      = gp (-(x * X)) n * BPsum (x ^ 2 * X) (fun r => X ^ tri r * x ^ r * α r) n := by
  set F : ℕ → ℕ → K⟦X⟧ := fun k r => X ^ tri k * (x ^ k * iq X (n - k) *
    (gp (-(x * X)) k * (α r * iq X (k - r) * gp (x ^ 2 * X * X ^ (k + r + 1)) (k - r))) *
      gp (x ^ 2 * X * X ^ (2 * k + 1)) (2 * n - 2 * k)) with hF
  have h1 : ∑ k ∈ range (n + 1), X ^ tri k * (x ^ k * iq X (n - k) * Δ k *
      gp (x ^ 2 * X * X ^ (2 * k + 1)) (2 * n - 2 * k))
      = ∑ k ∈ range (n + 1), ∑ r ∈ range (k + 1), F k r := by
    refine sum_congr rfl (fun k _ => ?_)
    rw [hΔ, BPsum]
    simp only [mul_sum, sum_mul, hF]
  rw [h1, sum_range_triangle, BPsum, mul_sum]
  refine sum_congr rfl (fun r hr => ?_)
  have hr' : r ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hr)
  obtain ⟨N, rfl⟩ := Nat.exists_eq_add_of_le hr'
  rw [Nat.add_sub_cancel_left]
  have hterm : ∀ i ∈ range (N + 1), F (r + i) r
      = (α r * gp (x ^ 2 * X * X ^ (r + N + r + 1)) N) *
        (X ^ tri (r + i) * (x ^ (r + i) * iq X (N - i) * iq X i *
          gp (-(x * X)) (r + i) * gp (x ^ 2 * X * X ^ (2 * r + i + 1)) (N - i))) := by
    intro i hi
    have hi' : i ≤ N := Nat.lt_succ_iff.mp (mem_range.mp hi)
    have hs := gp_split (x ^ 2 * X) r i N hi'
    simp only [hF]
    rw [Nat.add_sub_add_left, Nat.add_sub_cancel_left]
    linear_combination (X ^ tri (r + i) * x ^ (r + i) * iq X (N - i) * gp (-(x * X)) (r + i) *
      α r * iq X i) * hs
  rw [sum_congr rfl hterm, ← mul_sum, bailey_inner]
  ring

/-! ### The Bailey chain -/

/-- `G_i(n) = Σ_{n ≥ M_1 ≥ … ≥ M_i ≥ 0} ∏_l x^{M_l} q^{M_l(M_l+1)/2} / ((q)_{n-M_1} (q)_{M_1-M_2} ⋯ (q)_{M_i})`,
defined by peeling off the largest variable. -/
noncomputable def Gc (x : K⟦X⟧) : ℕ → ℕ → K⟦X⟧
  | 0, n => iq X n
  | i + 1, n => ∑ k ∈ range (n + 1), X ^ tri k * (x ^ k * iq X (n - k) * Gc x i k)

lemma Gc_zero (x : K⟦X⟧) (n : ℕ) : Gc x 0 n = iq X n := rfl

lemma Gc_succ (x : K⟦X⟧) (i n : ℕ) :
    Gc x (i + 1) n = ∑ k ∈ range (n + 1), X ^ tri k * (x ^ k * iq X (n - k) * Gc x i k) := rfl

/-- **The Bailey chain** `(S2)^{i+1}(UBP)` relative to `a = x²q`:
`(x²q²;q)_{2n} G_i(n) = (-xq;q)_n Σ_{r ≤ n} x^{(i+1)r} q^{(i+1)r(r+1)/2} α_r (x²q^{n+r+2};q)_{n-r}/(q)_{n-r}`. -/
theorem bailey_chain (x : K⟦X⟧) : ∀ (i n : ℕ),
    gp (x ^ 2 * X * X) (2 * n) * Gc x i n = gp (-(x * X)) n *
      BPsum (x ^ 2 * X) (fun r => X ^ ((i + 1) * tri r) * x ^ ((i + 1) * r) * alpha (x ^ 2 * X) r) n
  | 0, n => by
    have h := bailey_step x (alpha (x ^ 2 * X)) (fun k => if k = 0 then 1 else 0)
      (fun k => by
        rw [ubp]
        split_ifs with hk
        · subst hk; simp
        · simp) n
    rw [sum_range_succ', sum_eq_zero (fun k _ => by simp)] at h
    simp only [zero_add, tri_zero, pow_zero, Nat.sub_zero, mul_zero, if_true, mul_one, one_mul,
      pow_one] at h
    simp only [Gc_zero, zero_add, one_mul]
    rw [← h]; ring
  | i + 1, n => by
    have h := bailey_step x
      (fun r => X ^ ((i + 1) * tri r) * x ^ ((i + 1) * r) * alpha (x ^ 2 * X) r)
      (fun k => gp (x ^ 2 * X * X) (2 * k) * Gc x i k) (fun k => bailey_chain x i k) n
    have e : ∀ k ∈ range (n + 1), gp (x ^ 2 * X * X) (2 * k) *
        gp (x ^ 2 * X * X ^ (2 * k + 1)) (2 * n - 2 * k) = gp (x ^ 2 * X * X) (2 * n) := by
      intro k hk
      have hk' : k ≤ n := Nat.lt_succ_iff.mp (mem_range.mp hk)
      have a1 : x ^ 2 * X * X ^ (2 * k + 1) = x ^ 2 * X * X * X ^ (2 * k) := by ring
      have a2 : 2 * k + (2 * n - 2 * k) = 2 * n := by omega
      rw [a1, ← gp_add, a2]
    have hL : gp (x ^ 2 * X * X) (2 * n) * Gc x (i + 1) n
        = ∑ k ∈ range (n + 1), X ^ tri k * (x ^ k * iq X (n - k) *
          (gp (x ^ 2 * X * X) (2 * k) * Gc x i k) *
            gp (x ^ 2 * X * X ^ (2 * k + 1)) (2 * n - 2 * k)) := by
      rw [Gc_succ, mul_sum]
      refine sum_congr rfl (fun k hk => ?_)
      rw [← e k hk]; ring
    rw [hL, h]
    congr 1
    simp only [BPsum]
    refine sum_congr rfl (fun r _ => ?_)
    ring

end KRTW
