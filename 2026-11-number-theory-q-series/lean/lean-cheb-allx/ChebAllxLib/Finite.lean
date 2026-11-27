/-
ChebAllxLib/Finite.lean — an UNCONDITIONAL finite range: `Σ_p χ₋₄(p) e^{−p/x} < 0` for every
`0 < x ≤ 250000`, with no hypothesis on zeros.

Input: `Allx.theta1_wlog_nonneg_of_le` (lean-allx; kernel-checked): `Θ₁(t) ≥ 0` for `2 ≤ t ≤ 10⁷`.

Method.  Let `f(n) = e^{−n/x}/log n`, `M = 10⁷`, `δ = f(M) − f(M+1)`.  Split `f = h + r` with

    h(n) = f(n) − f(M) − δ (M − n)   (n ≤ M),      h(n) = 0   (n > M)

(`f` minus its secant through `M`, `M+1`, cut off at `M`).  `h` is convex on `n ≥ 3`, vanishes
beyond `M`, and has the same second difference as `f` at `n = 3`.  By the series identity
`hasSum_theta1_mul_secondDiff` applied to `h` — only `Θ₁(4), …, Θ₁(M)` occur —

    −Σ_p χ₋₄(p) log p · h(p) ≥ Θ₁(4) (f(3) − 2 f(4) + f(5)),      Θ₁(4) = log 3,

and the remainder is bounded trivially: `Σ_p log p · r(p) ≤ ((M+1)² + x) e^{−M/x}`.  For
`x ≤ 250000` the first quantity is `≥ 0.03 e^{−4/x}` and the second `≤ 0.002 e^{−4/x}`.
-/
import ChebAllxLib.Main

open Finset Filter Topology

noncomputable section

namespace ChebAllx

/-- `f` minus its secant through `M`, `M+1`, cut off at `M`. -/
def headPart (M : ℕ) (f : ℕ → ℝ) (n : ℕ) : ℝ :=
  if n ≤ M then f n - f M - (f M - f (M + 1)) * ((M : ℝ) - n) else 0

/-- The remainder `f − headPart M f`: the secant up to `M`, and `f` beyond. -/
def tailPart (M : ℕ) (f : ℕ → ℝ) (n : ℕ) : ℝ :=
  if n ≤ M then f M + (f M - f (M + 1)) * ((M : ℝ) - n) else f n

lemma head_add_tail (M : ℕ) (f : ℕ → ℝ) (n : ℕ) : headPart M f n + tailPart M f n = f n := by
  unfold headPart tailPart
  split_ifs <;> ring

lemma headPart_convex {M : ℕ} {f : ℕ → ℝ}
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2)) (n : ℕ) (hn : 3 ≤ n) :
    2 * headPart M f (n + 1) ≤ headPart M f n + headPart M f (n + 2) := by
  have h := hconv n hn
  unfold headPart
  rcases lt_trichotomy (n + 1) M with h1 | h1 | h1
  · rw [if_pos (by omega : n + 1 ≤ M), if_pos (by omega : n ≤ M), if_pos (by omega : n + 2 ≤ M)]
    push_cast
    linarith
  · rw [if_pos (by omega : n + 1 ≤ M), if_pos (by omega : n ≤ M),
      if_neg (by omega : ¬ n + 2 ≤ M)]
    subst h1
    push_cast
    rw [show n + 1 + 1 = n + 2 from rfl]
    linarith
  · rw [if_neg (by omega : ¬ n + 1 ≤ M), if_neg (by omega : ¬ n + 2 ≤ M)]
    by_cases h2 : n ≤ M
    · have hnM : n = M := by omega
      subst hnM
      rw [if_pos le_rfl]
      simp
    · rw [if_neg h2]
      simp

lemma summable_headPart (M : ℕ) (f : ℕ → ℝ) :
    Summable fun p : Nat.Primes => Real.log (p : ℕ) * headPart M f p := by
  refine summable_of_hasFiniteSupport
    (((Set.finite_Iic M).preimage (Subtype.val_injective.injOn)).subset ?_)
  intro p hp
  have hp' : Real.log ((p : ℕ) : ℝ) * headPart M f p ≠ 0 := hp
  by_contra h
  apply hp'
  have hpM : ¬ (p : ℕ) ≤ M := h
  unfold headPart
  rw [if_neg hpM, mul_zero]

lemma theta1_four :
    Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) ((0 + 4 : ℕ) : ℝ) = Real.log 3 := by
  rw [theta1_natCast]
  have h : (Finset.Icc 0 (0 + 4)).filter Nat.Prime = {2, 3} := by decide
  have h2 : ZMod.χ₄ ((2 : ℕ) : ZMod 4) = 0 := by decide
  have h3 : ZMod.χ₄ ((3 : ℕ) : ZMod 4) = -1 := by decide
  rw [h, Finset.sum_pair (by norm_num), h2, h3]
  push_cast
  ring

/-- **The head (no hypothesis).**  For `f` convex on `n ≥ 3` and `5 ≤ M ≤ 10⁷`:
`−Σ_p χ₋₄(p) log p · h(p) ≥ log 3 · (f(3) − 2 f(4) + f(5))`, `h = headPart M f`. -/
lemma log_three_mul_le_neg_tsum_head {M : ℕ} (hM5 : 5 ≤ M) (hMle : M ≤ 10000000) {f : ℕ → ℝ}
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2)) :
    Real.log 3 * (f 3 - 2 * f 4 + f 5)
      ≤ -∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
          * headPart M f p := by
  have hH := hasSum_theta1_mul_secondDiff (headPart_convex (M := M) hconv) (summable_headPart M f)
  refine le_trans (le_of_eq ?_) (sum_le_hasSum {0} (fun k _ => ?_) hH)
  · rw [Finset.sum_singleton, theta1_four]
    unfold headPart
    rw [if_pos (by omega : 3 ≤ M), if_pos (by omega : 4 ≤ M), if_pos (by omega : 5 ≤ M)]
    push_cast
    ring
  · by_cases hk : k + 4 ≤ M
    · refine mul_nonneg (Allx.theta1_wlog_nonneg_of_le ?_ ?_) ?_
      · push_cast
        linarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)]
      · exact_mod_cast hk.trans hMle
      · have := headPart_convex (M := M) hconv (k + 3) (by omega)
        rw [show k + 3 + 1 = k + 4 from rfl, show k + 3 + 2 = k + 5 from rfl] at this
        linarith
    · have h0 : headPart M f (k + 3) - 2 * headPart M f (k + 4) + headPart M f (k + 5) = 0 := by
        unfold headPart
        rw [if_neg (by omega : ¬ k + 4 ≤ M), if_neg (by omega : ¬ k + 5 ≤ M)]
        by_cases h3 : k + 3 ≤ M
        · have hkM : k + 3 = M := by omega
          rw [if_pos h3, ← hkM]
          simp
        · rw [if_neg h3]
          simp
      rw [h0, mul_zero]

/-! ### The weight `e^{−n/x}/log n` -/

/-- `f(n) = e^{−n/x}/log n`. -/
def chebSeq (x : ℝ) (n : ℕ) : ℝ := Real.exp (-(n : ℝ) / x) / Real.log n

lemma chebSeq_nonneg (x : ℝ) (n : ℕ) : 0 ≤ chebSeq x n :=
  div_nonneg (Real.exp_pos _).le (Real.log_natCast_nonneg n)

lemma chebSeq_convex {x : ℝ} (hx : 0 < x) (n : ℕ) (hn : 3 ≤ n) :
    2 * chebSeq x (n + 1) ≤ chebSeq x n + chebSeq x (n + 2) := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 3 := ⟨n - 3, by omega⟩
  set u : ℕ → ℝ := fun j => Real.exp (-((j + 3 : ℕ) : ℝ) / x) with hu
  set v : ℕ → ℝ := fun j => 1 / Real.log ((j + 3 : ℕ) : ℝ) with hv
  have hlogpos : ∀ j : ℕ, 0 < Real.log ((j + 3 : ℕ) : ℝ) := fun j =>
    Real.log_pos (by push_cast; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)])
  have hu0 : ∀ j, 0 ≤ u j := fun j => (Real.exp_pos _).le
  have hstep : ∀ j : ℕ, u (j + 1) = u j * Real.exp (-1 / x) := by
    intro j
    simp only [hu]
    rw [← Real.exp_add]
    congr 1
    push_cast
    ring
  have hr0 : 0 < Real.exp (-1 / x) := Real.exp_pos _
  have hr1 : Real.exp (-1 / x) ≤ 1 := by
    rw [Real.exp_le_one_iff]
    have : 0 < 1 / x := one_div_pos.2 hx
    rw [neg_div]
    linarith
  have hu1 : ∀ j, u (j + 1) ≤ u j := by
    intro j
    rw [hstep j]
    nlinarith [hu0 j]
  have hu2 : ∀ j, 0 ≤ secondDiff u j := by
    intro j
    unfold secondDiff
    rw [show j + 2 = j + 1 + 1 from rfl, hstep (j + 1), hstep j]
    nlinarith [mul_nonneg (hu0 j) (sq_nonneg (1 - Real.exp (-1 / x)))]
  have hv0 : ∀ j, 0 ≤ v j := fun j => (one_div_pos.2 (hlogpos j)).le
  have hv1 : ∀ j, v (j + 1) ≤ v j := by
    intro j
    simp only [hv]
    refine one_div_le_one_div_of_le (hlogpos j) (Real.log_le_log ?_ ?_)
    · push_cast; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
    · push_cast; linarith
  have hv2 : ∀ j, 0 ≤ secondDiff v j := by
    intro j
    have h := inv_log_convex j
    unfold secondDiff
    simp only [hv]
    rw [show j + 1 + 3 = j + 3 + 1 by omega, show j + 2 + 3 = j + 3 + 2 by omega]
    linarith
  have h := secondDiff_mul_nonneg hu0 hv0 hu1 hv1 hu2 hv2 k
  unfold secondDiff at h
  simp only [hu, hv] at h
  rw [show k + 1 + 3 = k + 3 + 1 by omega, show k + 2 + 3 = k + 3 + 2 by omega] at h
  unfold chebSeq
  simp only [div_eq_mul_one_div (Real.exp _)]
  linarith

lemma chebSeq_antitone {x : ℝ} (hx : 0 < x) {m n : ℕ} (hm : 2 ≤ m) (hmn : m ≤ n) :
    chebSeq x n ≤ chebSeq x m := by
  unfold chebSeq
  have hm' : (1 : ℝ) < m := by exact_mod_cast (by omega : 1 < m)
  have hmn' : (m : ℝ) ≤ n := by exact_mod_cast hmn
  have hlm : 0 < Real.log m := Real.log_pos hm'
  have hln : Real.log m ≤ Real.log n := Real.log_le_log (by linarith) hmn'
  have he : Real.exp (-(n : ℝ) / x) ≤ Real.exp (-(m : ℝ) / x) := by
    rw [Real.exp_le_exp]
    exact div_le_div_of_nonneg_right (by linarith) hx.le
  exact div_le_div₀ (Real.exp_pos _).le he hlm hln

/-- The trivial bound for the remainder: `Σ_p log p · r(p) ≤ ((M+1)² + x) e^{−M/x}`. -/
lemma tsum_log_mul_tail_le {x : ℝ} (hx : 0 < x) {M : ℕ} (hM : 2 ≤ M)
    (hs : Summable fun p : Nat.Primes => Real.log (p : ℕ) * tailPart M (chebSeq x) p) :
    ∑' p : Nat.Primes, Real.log (p : ℕ) * tailPart M (chebSeq x) p
      ≤ (((M : ℝ) + 1) ^ 2 + x) * Real.exp (-(M : ℝ) / x) := by
  set B : ℕ → ℝ := fun n =>
    if n ≤ M then ((M : ℝ) + 1) * Real.exp (-(M : ℝ) / x) else Real.exp (-(n : ℝ) / x) with hB
  set r : ℝ := Real.exp (-1 / x) with hr
  have hr0 : 0 < r := Real.exp_pos _
  have hr1 : r < 1 := by
    rw [hr, Real.exp_lt_one_iff]
    have : 0 < 1 / x := one_div_pos.2 hx
    rw [neg_div]
    linarith
  have hBshift : ∀ n : ℕ, B (n + (M + 1)) = Real.exp (-((M : ℝ) + 1) / x) * r ^ n := by
    intro n
    simp only [hB]
    rw [if_neg (by omega), hr, ← Real.exp_nat_mul, ← Real.exp_add]
    congr 1
    push_cast
    ring
  have hgeo : Summable fun n : ℕ => Real.exp (-((M : ℝ) + 1) / x) * r ^ n :=
    (summable_geometric_of_lt_one hr0.le hr1).mul_left _
  have hBs : Summable B := by
    rw [← summable_nat_add_iff (M + 1)]
    exact hgeo.congr fun n => (hBshift n).symm
  have hBsum : ∑' n, B n
      = ((M : ℝ) + 1) * (((M : ℝ) + 1) * Real.exp (-(M : ℝ) / x))
        + Real.exp (-((M : ℝ) + 1) / x) * (1 - r)⁻¹ := by
    rw [← hBs.sum_add_tsum_nat_add (M + 1)]
    congr 1
    · rw [Finset.sum_congr rfl (fun n hn => by
        simp only [hB]
        rw [if_pos (by have := Finset.mem_range.1 hn; omega)]), Finset.sum_const,
        Finset.card_range, nsmul_eq_mul]
      push_cast
      ring
    · rw [tsum_congr hBshift, tsum_mul_left, tsum_geometric_of_lt_one hr0.le hr1]
  -- the geometric tail: `e^{−(M+1)/x}/(1 − r) ≤ x e^{−M/x}`
  have htail : Real.exp (-((M : ℝ) + 1) / x) * (1 - r)⁻¹ ≤ x * Real.exp (-(M : ℝ) / x) := by
    have e1 : Real.exp (-((M : ℝ) + 1) / x) = Real.exp (-(M : ℝ) / x) * r := by
      rw [hr, ← Real.exp_add]
      congr 1
      ring
    have h1r : 0 < 1 - r := by linarith
    have hkey : r ≤ x * (1 - r) := by
      -- `e^{1/x} ≥ 1 + 1/x`, i.e. `1 ≥ r (1 + 1/x)`
      have h := Real.add_one_le_exp (1 / x)
      have hrinv : Real.exp (1 / x) * r = 1 := by
        rw [hr, ← Real.exp_add]
        simp [neg_div]
      have hx1 : 1 / x * x = 1 := by field_simp
      nlinarith [mul_le_mul_of_nonneg_right h hr0.le, mul_pos hr0 hx]
    rw [e1, mul_assoc, mul_comm x]
    refine mul_le_mul_of_nonneg_left ?_ (Real.exp_pos _).le
    rw [← div_eq_mul_inv, div_le_iff₀ h1r]
    exact hkey
  -- termwise comparison
  have hterm : ∀ n : ℕ, (if n.Prime then Real.log n * tailPart M (chebSeq x) n else 0) ≤ B n := by
    intro n
    have hBnn : 0 ≤ B n := by
      simp only [hB]
      split_ifs <;> positivity
    by_cases hp : n.Prime
    · rw [if_pos hp]
      have hn2 : 2 ≤ n := hp.two_le
      have hlogn : 0 < Real.log n := Real.log_pos (by exact_mod_cast (by omega : 1 < n))
      unfold tailPart
      simp only [hB]
      by_cases hnM : n ≤ M
      · rw [if_pos hnM, if_pos hnM]
        have hfM : 0 ≤ chebSeq x M := chebSeq_nonneg x M
        have hfM1 : 0 ≤ chebSeq x (M + 1) := chebSeq_nonneg x (M + 1)
        have hMn : (0 : ℝ) ≤ (M : ℝ) - n := sub_nonneg.2 (by exact_mod_cast hnM)
        have hMn' : (M : ℝ) - n ≤ M := by linarith [(Nat.cast_nonneg n : (0 : ℝ) ≤ n)]
        have hlogM : 0 < Real.log M := Real.log_pos (by exact_mod_cast (by omega : 1 < M))
        have hlognM : Real.log n ≤ Real.log M :=
          Real.log_le_log (by exact_mod_cast (by omega : 0 < n)) (by exact_mod_cast hnM)
        have hsec : chebSeq x M + (chebSeq x M - chebSeq x (M + 1)) * ((M : ℝ) - n)
            ≤ chebSeq x M * ((M : ℝ) + 1) := by
          have h1 : (chebSeq x M - chebSeq x (M + 1)) * ((M : ℝ) - n) ≤ chebSeq x M * M :=
            mul_le_mul (by linarith) hMn' hMn hfM
          linarith
        have hsec0 : 0 ≤ chebSeq x M + (chebSeq x M - chebSeq x (M + 1)) * ((M : ℝ) - n) :=
          add_nonneg hfM (mul_nonneg
            (sub_nonneg.2 (chebSeq_antitone hx hM (Nat.le_succ M))) hMn)
        calc Real.log n * (chebSeq x M + (chebSeq x M - chebSeq x (M + 1)) * ((M : ℝ) - n))
            ≤ Real.log M * (chebSeq x M * ((M : ℝ) + 1)) :=
              mul_le_mul hlognM hsec hsec0 hlogM.le
          _ = ((M : ℝ) + 1) * Real.exp (-(M : ℝ) / x) := by
              unfold chebSeq
              field_simp
      · rw [if_neg hnM, if_neg hnM]
        unfold chebSeq
        rw [mul_div_cancel₀ _ hlogn.ne']
    · rw [if_neg hp]
      exact hBnn
  have hsN : Summable fun n : ℕ => if n.Prime then Real.log n * tailPart M (chebSeq x) n else 0 := by
    have h1 := (summable_subtype_iff_indicator (s := {n : ℕ | n.Prime})
      (f := fun n : ℕ => Real.log n * tailPart M (chebSeq x) n)).1 hs
    refine h1.congr fun n => ?_
    by_cases hp : n.Prime
    · rw [Set.indicator_of_mem (show n ∈ {n : ℕ | n.Prime} from hp), if_pos hp]
    · rw [Set.indicator_of_notMem (show n ∉ {n : ℕ | n.Prime} from hp), if_neg hp]
  rw [ChebLib.tsum_primes_eq (fun n => Real.log n * tailPart M (chebSeq x) n)]
  calc ∑' n : ℕ, (if n.Prime then Real.log n * tailPart M (chebSeq x) n else 0)
      ≤ ∑' n, B n := hsN.tsum_le_tsum hterm hBs
    _ = ((M : ℝ) + 1) * (((M : ℝ) + 1) * Real.exp (-(M : ℝ) / x))
        + Real.exp (-((M : ℝ) + 1) / x) * (1 - r)⁻¹ := hBsum
    _ ≤ ((M : ℝ) + 1) * (((M : ℝ) + 1) * Real.exp (-(M : ℝ) / x))
        + x * Real.exp (-(M : ℝ) / x) := by linarith
    _ = (((M : ℝ) + 1) ^ 2 + x) * Real.exp (-(M : ℝ) / x) := by ring

lemma tailPart_chebSeq_nonneg {x : ℝ} (hx : 0 < x) {M : ℕ} (hM : 2 ≤ M) (n : ℕ) :
    0 ≤ tailPart M (chebSeq x) n := by
  unfold tailPart
  split_ifs with h
  · exact add_nonneg (chebSeq_nonneg x M) (mul_nonneg
      (sub_nonneg.2 (chebSeq_antitone hx hM (Nat.le_succ M)))
      (sub_nonneg.2 (by exact_mod_cast h)))
  · exact chebSeq_nonneg x n

lemma two_le_add_one_div {t : ℝ} (ht : 0 < t) : 2 ≤ t + 1 / t := by
  have h : t + 1 / t - 2 = (t - 1) ^ 2 / t := by
    field_simp
    ring
  have h2 : 0 ≤ (t - 1) ^ 2 / t := by positivity
  linarith

/-- The numerical comparison: for `0 < x ≤ 250000` and `M = 10⁷`,
`((M+1)² + x) e^{−M/x} < log 3 · (f(3) − 2 f(4) + f(5))`, `f(n) = e^{−n/x}/log n`. -/
lemma finite_numeric {x : ℝ} (hx : 0 < x) (hxN : x ≤ 250000) :
    ((((10000000 : ℕ) : ℝ) + 1) ^ 2 + x) * Real.exp (-((10000000 : ℕ) : ℝ) / x)
      < Real.log 3 * (chebSeq x 3 - 2 * chebSeq x 4 + chebSeq x 5) := by
  set E4 := Real.exp (-4 / x) with hE4
  set y := Real.exp (1 / x) with hy
  have hE4pos : 0 < E4 := Real.exp_pos _
  have hypos : 0 < y := Real.exp_pos _
  have a_lo := Real.log_three_gt_d9
  have a_hi := Real.log_three_lt_d9
  have b_lo := Real.log_five_gt_d9
  have b_hi := Real.log_five_lt_d9
  have l2_lo := Real.log_two_gt_d9
  have l2_hi := Real.log_two_lt_d9
  have ha : 0 < Real.log 3 := by linarith
  have hb : 0 < Real.log 5 := by linarith
  have hl2 : 0 < Real.log 2 := by linarith
  have e3 : chebSeq x 3 = E4 * y / Real.log 3 := by
    unfold chebSeq
    rw [hE4, hy, ← Real.exp_add]
    push_cast
    congr 2
    ring
  have e4 : chebSeq x 4 = E4 / (2 * Real.log 2) := by
    unfold chebSeq
    push_cast
    rw [Real.log_four_eq]
  have e5 : chebSeq x 5 = E4 / y / Real.log 5 := by
    unfold chebSeq
    rw [hE4, hy, ← Real.exp_sub]
    push_cast
    congr 2
    ring
  -- `y/a + 1/(y b) ≥ 4/(a+b)` (twice AM–GM), and `4/(a+b) − 1/log 2 ≥ 0.03`
  have hcore : (0.03 : ℝ) ≤ y / Real.log 3 - 2 / (2 * Real.log 2) + 1 / y / Real.log 5 := by
    have h1 := two_le_add_one_div hypos
    have h2 := two_le_add_one_div (t := y * Real.log 5 / Real.log 3) (by positivity)
    have e : (y / Real.log 3 + 1 / y / Real.log 5) * (Real.log 3 + Real.log 5)
        = (y + 1 / y) + (y * Real.log 5 / Real.log 3 + 1 / (y * Real.log 5 / Real.log 3)) := by
      field_simp
      ring
    have key : 4 / (Real.log 3 + Real.log 5) ≤ y / Real.log 3 + 1 / y / Real.log 5 := by
      rw [div_le_iff₀ (by positivity)]
      linarith
    have h3 : (1.4770 : ℝ) ≤ 4 / (Real.log 3 + Real.log 5) := by
      rw [le_div_iff₀ (by positivity)]
      linarith
    have h4 : 2 / (2 * Real.log 2) ≤ (1.4428 : ℝ) := by
      rw [div_le_iff₀ (by positivity)]
      linarith
    linarith
  have hR : 0.03 * E4 ≤ Real.log 3 * (chebSeq x 3 - 2 * chebSeq x 4 + chebSeq x 5) := by
    rw [e3, e4, e5]
    have e : Real.log 3 * (E4 * y / Real.log 3 - 2 * (E4 / (2 * Real.log 2))
        + E4 / y / Real.log 5)
        = Real.log 3 * E4 * (y / Real.log 3 - 2 / (2 * Real.log 2) + 1 / y / Real.log 5) := by
      ring
    rw [e]
    have h1 : 0 ≤ Real.log 3 * E4 * (y / Real.log 3 - 2 / (2 * Real.log 2)
        + 1 / y / Real.log 5 - 0.03) :=
      mul_nonneg (mul_pos ha hE4pos).le (by linarith)
    have h2 : 0 ≤ (Real.log 3 - 1) * E4 := mul_nonneg (by linarith) hE4pos.le
    nlinarith
  -- the left side
  have hM : (((10000000 : ℕ) : ℝ)) = 10000000 := by norm_num
  have hexp : Real.exp (-((10000000 : ℕ) : ℝ) / x) ≤ E4 * Real.exp (-39) := by
    rw [hE4, ← Real.exp_add, Real.exp_le_exp, hM]
    have h39 : (39 : ℝ) ≤ (10000000 - 4) / x := by
      rw [le_div_iff₀ hx]
      linarith
    have e : -(10000000 : ℝ) / x = -4 / x + -((10000000 - 4) / x) := by ring
    rw [e]
    linarith
  have h39 : Real.exp (-39) < 1 / 66000000000000000 := by
    have hpow : Real.exp 39 = Real.exp 1 ^ 39 := by
      rw [← Real.exp_nat_mul]
      norm_num
    have hbig : (66000000000000000 : ℝ) < Real.exp 39 := by
      rw [hpow]
      calc (66000000000000000 : ℝ) < (2.7 : ℝ) ^ 39 := by norm_num
        _ ≤ Real.exp 1 ^ 39 :=
          pow_le_pow_left₀ (by norm_num) (by linarith [Real.exp_one_gt_d9]) 39
    rw [Real.exp_neg, one_div]
    exact inv_strictAnti₀ (by norm_num) hbig
  have hfac : (((10000000 : ℕ) : ℝ) + 1) ^ 2 + x ≤ 110000000000000 := by
    rw [hM]
    have : ((10000000 : ℝ) + 1) ^ 2 = 100000020000001 := by norm_num
    linarith
  have hexp39 : 0 < Real.exp (-39) := Real.exp_pos _
  calc ((((10000000 : ℕ) : ℝ) + 1) ^ 2 + x) * Real.exp (-((10000000 : ℕ) : ℝ) / x)
      ≤ 110000000000000 * (E4 * Real.exp (-39)) :=
        mul_le_mul hfac hexp (Real.exp_pos _).le (by norm_num)
    _ < 0.03 * E4 := by nlinarith [mul_pos hE4pos hexp39]
    _ ≤ Real.log 3 * (chebSeq x 3 - 2 * chebSeq x 4 + chebSeq x 5) := hR

/-- **Chebyshev's sum is negative on a finite range, with NO hypothesis.**  For every real
`0 < x ≤ 250000`:  `Σ_p χ₋₄(p) e^{−p/x} < 0`. -/
theorem chebyshevSum_neg_of_le {x : ℝ} (hx : 0 < x) (hxN : x ≤ 250000) : chebyshevSum x < 0 := by
  have hconv := chebSeq_convex hx
  have hlog : ∀ p : Nat.Primes, Real.log (p : ℕ) ≠ 0 := fun p =>
    (Real.log_pos (by exact_mod_cast p.2.one_lt)).ne'
  have hsf : Summable fun p : Nat.Primes => Real.log (p : ℕ) * chebSeq x p := by
    refine (summable_primes_exp_neg_div hx).congr fun p => ?_
    unfold chebSeq
    have := hlog p
    field_simp
  have hsh := summable_headPart 10000000 (chebSeq x)
  have hst : Summable fun p : Nat.Primes =>
      Real.log (p : ℕ) * tailPart 10000000 (chebSeq x) p := by
    refine (hsf.sub hsh).congr fun p => ?_
    rw [← head_add_tail 10000000 (chebSeq x) p]
    ring
  have hχ : ∀ p : Nat.Primes, |((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)| ≤ 1 := by
    intro p
    have h := abs_w4_le_one (p : ℕ)
    unfold Coupling.w4 at h
    rwa [abs_neg] at h
  have hbound : ∀ (g : ℕ → ℝ) (p : Nat.Primes),
      ‖((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * g p‖
        ≤ ‖Real.log (p : ℕ) * g p‖ := by
    intro g p
    rw [mul_assoc, norm_mul, Real.norm_eq_abs]
    exact (mul_le_mul_of_nonneg_right (hχ p) (norm_nonneg _)).trans (by rw [one_mul])
  have hsχh : Summable fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.log (p : ℕ) * headPart 10000000 (chebSeq x) p :=
    Summable.of_norm_bounded hsh.norm (hbound _)
  have hsχt : Summable fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.log (p : ℕ) * tailPart 10000000 (chebSeq x) p :=
    Summable.of_norm_bounded hst.norm (hbound _)
  have hC : chebyshevSum x
      = ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
          * headPart 10000000 (chebSeq x) p
        + ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
          * tailPart 10000000 (chebSeq x) p := by
    rw [← hsχh.tsum_add hsχt]
    unfold chebyshevSum
    refine tsum_congr fun p => ?_
    rw [← mul_add, head_add_tail]
    unfold chebSeq
    have := hlog p
    field_simp
  have hhead := log_three_mul_le_neg_tsum_head (M := 10000000) (by norm_num) le_rfl hconv
  have htail1 : ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
        * tailPart 10000000 (chebSeq x) p
      ≤ ∑' p : Nat.Primes, Real.log (p : ℕ) * tailPart 10000000 (chebSeq x) p := by
    refine hsχt.tsum_le_tsum (fun p => ?_) hst
    have h0 : 0 ≤ Real.log (p : ℕ) * tailPart 10000000 (chebSeq x) p :=
      mul_nonneg (Real.log_natCast_nonneg _) (tailPart_chebSeq_nonneg hx (by norm_num) _)
    have h1 := (abs_le.1 (hχ p)).2
    rw [mul_assoc]
    nlinarith
  have htail2 := tsum_log_mul_tail_le hx (M := 10000000) (by norm_num) hst
  have hnum := finite_numeric hx hxN
  rw [hC]
  linarith

end ChebAllx
