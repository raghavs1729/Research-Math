/-
ChebAllxLib/Discrete.lean — second-order summation by parts over the INFINITE range, in discrete
form (Mathlib only; no number theory).

For a real sequence `f : ℕ → ℝ` put `Δ²f(k) = f(k) − 2 f(k+1) + f(k+2)` (`secondDiff`).  If
`Δ²f ≥ 0` (f is a convex sequence) and `f(k) → 0`, then

* `f` is nonincreasing and nonnegative (`antitone_of_secondDiff_nonneg`, `nonneg_of_secondDiff_nonneg`);
* `f(n) = Σ_{j ≥ 0} (j+1) Δ²f(n+j)` for every `n` (`hasSum_mul_secondDiff_shift`): `f` is a
  mixture, with the nonnegative weights `Δ²f(k)`, of the kernels `n ↦ (k+1−n)₊`;
* for every sequence `a` with `Σ |a(n)| f(n) < ∞`:

      Σ_n a(n) f(n) = Σ_k A₁(k+1) Δ²f(k),      A₁(m) = Σ_{n < m} a(n) (m − n)

  (`hasSum_secondSum_mul_secondDiff`).  There are no boundary terms and no growth hypothesis on
  `a`: the exchange of the two sums is Tonelli's theorem for the absolutely summable double family
  `a(n) (j+1) Δ²f(n+j)`, whose rows sum to `a(n) f(n)`.
* consequently, if `A₁(k+1) > 0` for every `k` and `f(0) > 0`, then `Σ_n a(n) f(n) > 0`
  (`tsum_mul_pos_of_secondSum_pos`).

This is the discrete counterpart of `Σ_p c(p) F(p) = ∫ Θ₁(t) F''(t) dt` (Lemma 1 of
`monotonicity/chebyshev1853/CHEBYSHEV1853.md`); it needs the second sum `A₁ = Θ₁` only at the
integers and no differentiability of the weight.
-/
import Mathlib

open Finset Filter Topology

namespace ChebAllx

/-- The second difference `Δ²f(k) = f(k) − 2 f(k+1) + f(k+2)` of a real sequence. -/
def secondDiff (f : ℕ → ℝ) (k : ℕ) : ℝ := f k - 2 * f (k + 1) + f (k + 2)

/-- The second sum `A₁(m) = Σ_{n < m} a(n) (m − n)` of a real sequence (the discrete `Θ₁`). -/
def secondSum (a : ℕ → ℝ) (m : ℕ) : ℝ := ∑ n ∈ range m, a n * ((m : ℝ) - n)

section Convex

variable {f : ℕ → ℝ}

lemma secondDiff_shift (f : ℕ → ℝ) (n j : ℕ) :
    secondDiff (fun i => f (n + i)) j = secondDiff f (n + j) := rfl

/-- For a convex sequence the first differences `f(k) − f(k+1)` are nonincreasing. -/
lemma diff_antitone (hD : ∀ k, 0 ≤ secondDiff f k) : Antitone (fun k => f k - f (k + 1)) := by
  refine antitone_nat_of_succ_le fun k => ?_
  have h := hD k
  unfold secondDiff at h
  show f (k + 1) - f (k + 1 + 1) ≤ f k - f (k + 1)
  rw [show k + 1 + 1 = k + 2 from rfl]
  linarith

lemma diff_nonneg (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0)) (k : ℕ) :
    0 ≤ f k - f (k + 1) := by
  have hG : Tendsto (fun k => f k - f (k + 1)) atTop (𝓝 0) := by
    simpa using h0.sub ((tendsto_add_atTop_iff_nat 1).2 h0)
  exact (diff_antitone hD).le_of_tendsto hG k

/-- A convex sequence tending to `0` is nonincreasing. -/
lemma antitone_of_secondDiff_nonneg (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0)) :
    Antitone f :=
  antitone_nat_of_succ_le fun k => by linarith [diff_nonneg hD h0 k]

/-- A convex sequence tending to `0` is nonnegative. -/
lemma nonneg_of_secondDiff_nonneg (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0))
    (k : ℕ) : 0 ≤ f k :=
  (antitone_of_secondDiff_nonneg hD h0).le_of_tendsto h0 k

/-- Summation by parts twice, finite form. -/
lemma sum_range_mul_secondDiff (f : ℕ → ℝ) (K : ℕ) :
    ∑ j ∈ range K, ((j : ℝ) + 1) * secondDiff f j = f 0 - f K - K * (f K - f (K + 1)) := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [sum_range_succ, ih]
    unfold secondDiff
    push_cast
    ring

lemma mul_diff_le (hD : ∀ k, 0 ≤ secondDiff f k) (m K : ℕ) (hmK : m ≤ K) :
    ((K : ℝ) - m) * (f K - f (K + 1)) ≤ f m - f K := by
  induction K, hmK using Nat.le_induction with
  | base => simp
  | succ K hK ih =>
    have h1 := hD K
    unfold secondDiff at h1
    have h2 : (0 : ℝ) ≤ (K : ℝ) - m := sub_nonneg.2 (Nat.cast_le.2 hK)
    have h3 : (0 : ℝ) ≤ ((K : ℝ) + 1 - m) * ((f K - f (K + 1)) - (f (K + 1) - f (K + 2))) :=
      mul_nonneg (by linarith) (by linarith)
    push_cast
    rw [show K + 1 + 1 = K + 2 from rfl]
    nlinarith

/-- For a convex sequence tending to `0`: `K (f(K) − f(K+1)) → 0`. -/
lemma tendsto_mul_diff (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0)) :
    Tendsto (fun K : ℕ => (K : ℝ) * (f K - f (K + 1))) atTop (𝓝 0) := by
  have hhalf : Tendsto (fun K : ℕ => K / 2) atTop atTop :=
    tendsto_atTop_atTop.2 fun b => ⟨2 * b, fun a ha => by omega⟩
  have hup : Tendsto (fun K : ℕ => 2 * f (K / 2)) atTop (𝓝 0) := by
    simpa using (h0.comp hhalf).const_mul 2
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hup ?_ ?_
  · intro K
    exact mul_nonneg (Nat.cast_nonneg K) (diff_nonneg hD h0 K)
  · intro K
    have h1 := mul_diff_le hD (K / 2) K (Nat.div_le_self K 2)
    have h2 := nonneg_of_secondDiff_nonneg hD h0 K
    have h3 := diff_nonneg hD h0 K
    have h4 : ((K / 2 : ℕ) : ℝ) * 2 ≤ K := by exact_mod_cast Nat.div_mul_le_self K 2
    show (K : ℝ) * (f K - f (K + 1)) ≤ 2 * f (K / 2)
    nlinarith

/-- **Mixture representation at `0`.**  A convex sequence tending to `0` satisfies
`f(0) = Σ_{j ≥ 0} (j+1) Δ²f(j)`. -/
theorem hasSum_mul_secondDiff (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0)) :
    HasSum (fun j : ℕ => ((j : ℝ) + 1) * secondDiff f j) (f 0) := by
  rw [hasSum_iff_tendsto_nat_of_nonneg (fun j => mul_nonneg (by positivity) (hD j))]
  simp_rw [sum_range_mul_secondDiff]
  have h := ((tendsto_const_nhds (x := f 0)).sub h0).sub (tendsto_mul_diff hD h0)
  simpa using h

/-- **Mixture representation.**  A convex sequence tending to `0` satisfies
`f(n) = Σ_{j ≥ 0} (j+1) Δ²f(n+j)` for every `n`. -/
theorem hasSum_mul_secondDiff_shift (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0))
    (n : ℕ) : HasSum (fun j : ℕ => ((j : ℝ) + 1) * secondDiff f (n + j)) (f n) := by
  have hD' : ∀ k, 0 ≤ secondDiff (fun i => f (n + i)) k := fun k => hD (n + k)
  have h0' : Tendsto (fun i => f (n + i)) atTop (𝓝 0) := by
    have := (tendsto_add_atTop_iff_nat n).2 h0
    simpa [add_comm] using this
  exact hasSum_mul_secondDiff hD' h0'

/-- If a convex sequence tending to `0` has `f(0) > 0`, some second difference is positive. -/
lemma exists_secondDiff_pos (hD : ∀ k, 0 ≤ secondDiff f k) (h0 : Tendsto f atTop (𝓝 0))
    (hf0 : 0 < f 0) : ∃ k, 0 < secondDiff f k := by
  by_contra hcon
  push Not at hcon
  have hz : ∀ k, secondDiff f k = 0 := fun k => le_antisymm (hcon k) (hD k)
  have h := hasSum_mul_secondDiff hD h0
  simp only [hz, mul_zero] at h
  exact hf0.ne' (h.unique hasSum_zero)

end Convex

/-- **Second-order summation by parts over the infinite range.**  Let `f` be a convex real sequence
tending to `0` and `a` a real sequence with `Σ |a(n)| f(n) < ∞`.  Then

    Σ_n a(n) f(n) = Σ_k A₁(k+1) Δ²f(k),    A₁(m) = Σ_{n < m} a(n) (m − n),

the series on the right being convergent. -/
theorem hasSum_secondSum_mul_secondDiff {a f : ℕ → ℝ} (hD : ∀ k, 0 ≤ secondDiff f k)
    (h0 : Tendsto f atTop (𝓝 0)) (hs : Summable fun n => |a n| * f n) :
    HasSum (fun k : ℕ => secondSum a (k + 1) * secondDiff f k) (∑' n, a n * f n) := by
  -- the double family `g(n, j) = a(n) (j+1) Δ²f(n+j)`; its rows sum to `a(n) f(n)`
  set g : ℕ × ℕ → ℝ := fun x => a x.1 * (((x.2 : ℝ) + 1) * secondDiff f (x.1 + x.2)) with hg
  have hrow : ∀ n, HasSum (fun j => g (n, j)) (a n * f n) := fun n =>
    (hasSum_mul_secondDiff_shift hD h0 n).mul_left (a n)
  have habs : ∀ x, |g x| = |a x.1| * (((x.2 : ℝ) + 1) * secondDiff f (x.1 + x.2)) := by
    intro x
    simp only [hg, abs_mul]
    rw [abs_of_nonneg (by positivity : (0 : ℝ) ≤ (x.2 : ℝ) + 1), abs_of_nonneg (hD _)]
  have hrowabs : ∀ n, HasSum (fun j => |g (n, j)|) (|a n| * f n) := by
    intro n
    simp_rw [habs]
    exact (hasSum_mul_secondDiff_shift hD h0 n).mul_left |a n|
  have hgabs : Summable (fun x => |g x|) := by
    rw [summable_prod_of_nonneg (fun x => abs_nonneg _)]
    refine ⟨fun n => (hrowabs n).summable, ?_⟩
    have : (fun n => ∑' j, |g (n, j)|) = fun n => |a n| * f n :=
      funext fun n => (hrowabs n).tsum_eq
    rw [this]
    exact hs
  have hgs : Summable g := Summable.of_abs hgabs
  have htot : HasSum g (∑' n, a n * f n) := by
    have h := hgs.hasSum
    rw [hgs.tsum_prod' (fun n => (hrow n).summable)] at h
    have : (fun n => ∑' j, g (n, j)) = fun n => a n * f n :=
      funext fun n => (hrow n).tsum_eq
    rwa [this] at h
  -- regroup along the antidiagonals `n + j = k`
  have h1 := (Finset.HasAntidiagonal.sigmaAntidiagonalEquivProd (A := ℕ)).hasSum_iff.2 htot
  have h2 := h1.sigma (fun k => hasSum_fintype _)
  convert h2 using 1
  funext k
  simp only [Function.comp_apply, Finset.HasAntidiagonal.sigmaAntidiagonalEquivProd_apply]
  rw [Finset.sum_coe_sort (antidiagonal k) g, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk,
    secondSum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun n hn => ?_
  have hnk : n ≤ k := Nat.lt_succ_iff.1 (Finset.mem_range.1 hn)
  simp only [hg, Nat.add_sub_cancel' hnk, Nat.cast_sub hnk]
  push_cast
  ring

/-- The same, written out (no auxiliary definitions). -/
theorem hasSum_secondSum_mul_secondDiff' {a f : ℕ → ℝ}
    (hconv : ∀ k, 2 * f (k + 1) ≤ f k + f (k + 2)) (h0 : Tendsto f atTop (𝓝 0))
    (hs : Summable fun n => |a n| * f n) :
    HasSum (fun k : ℕ => (∑ n ∈ range (k + 1), a n * ((k : ℝ) + 1 - n))
        * (f k - 2 * f (k + 1) + f (k + 2))) (∑' n, a n * f n) := by
  have hD : ∀ k, 0 ≤ secondDiff f k := fun k => by
    unfold secondDiff
    linarith [hconv k]
  have h := hasSum_secondSum_mul_secondDiff hD h0 hs
  simpa only [secondSum, secondDiff, Nat.cast_add, Nat.cast_one] using h

/-- **Positivity.**  Let `f` be a convex real sequence tending to `0` with `f(0) > 0`, and `a` a real
sequence with `Σ |a(n)| f(n) < ∞` all of whose second sums `A₁(k+1) = Σ_{n ≤ k} a(n) (k+1−n)` are
positive.  Then `Σ_n a(n) f(n) > 0`. -/
theorem tsum_mul_pos_of_secondSum_pos {a f : ℕ → ℝ} (hD : ∀ k, 0 ≤ secondDiff f k)
    (h0 : Tendsto f atTop (𝓝 0)) (hs : Summable fun n => |a n| * f n)
    (hA : ∀ k, 0 < secondSum a (k + 1)) (hf0 : 0 < f 0) : 0 < ∑' n, a n * f n := by
  have h := hasSum_secondSum_mul_secondDiff hD h0 hs
  obtain ⟨k, hk⟩ := exists_secondDiff_pos hD h0 hf0
  rw [← h.tsum_eq]
  exact h.summable.tsum_pos (fun i => mul_nonneg (hA i).le (hD i)) k (mul_pos (hA k) hk)

/-- The product of two nonnegative, nonincreasing, convex sequences is convex:
`Δ²(uv)(k) = u(k+2) Δ²v(k) + v(k+1) Δ²u(k) + (u(k+2) − u(k+1) + u(k+1) − u(k)) (v(k+1) − v(k))`. -/
lemma secondDiff_mul_nonneg {u v : ℕ → ℝ} (hu : ∀ k, 0 ≤ u k) (hv : ∀ k, 0 ≤ v k)
    (hu1 : ∀ k, u (k + 1) ≤ u k) (hv1 : ∀ k, v (k + 1) ≤ v k)
    (hu2 : ∀ k, 0 ≤ secondDiff u k) (hv2 : ∀ k, 0 ≤ secondDiff v k) (k : ℕ) :
    0 ≤ secondDiff (fun n => u n * v n) k := by
  have key : secondDiff (fun n => u n * v n) k
      = u (k + 2) * secondDiff v k + v (k + 1) * secondDiff u k
        + (u (k + 2) - u k) * (v (k + 1) - v k) := by
    unfold secondDiff
    ring
  rw [key]
  have h1 := mul_nonneg (hu (k + 2)) (hv2 k)
  have h2 := mul_nonneg (hv (k + 1)) (hu2 k)
  have h3 : 0 ≤ (u (k + 2) - u k) * (v (k + 1) - v k) :=
    mul_nonneg_of_nonpos_of_nonpos (by linarith [hu1 k, hu1 (k + 1)]) (by linarith [hv1 k])
  linarith

/-- A convex sequence which is arbitrarily small infinitely often tends to `0` (it is then
nonincreasing and nonnegative). -/
lemma tendsto_zero_of_secondDiff_nonneg {g : ℕ → ℝ} (hD : ∀ k, 0 ≤ secondDiff g k)
    (hfreq : ∀ ε > 0, ∀ N, ∃ n ≥ N, |g n| < ε) : Tendsto g atTop (𝓝 0) := by
  -- the first differences are nondecreasing
  have hmono : ∀ k m, k ≤ m → g (k + 1) - g k ≤ g (m + 1) - g m := by
    intro k m hkm
    induction m, hkm using Nat.le_induction with
    | base => exact le_rfl
    | succ m _ ih =>
      have := hD m
      unfold secondDiff at this
      rw [show m + 1 + 1 = m + 2 from rfl]
      linarith
  -- hence `g` is nonincreasing: one positive difference would force `g → +∞`
  have hanti : ∀ k, g (k + 1) ≤ g k := by
    intro k
    by_contra hk
    push Not at hk
    set δ := g (k + 1) - g k with hδ
    have hδ0 : 0 < δ := by linarith
    have hlin : ∀ n, k ≤ n → g k + ((n : ℝ) - k) * δ ≤ g n := by
      intro n hn
      induction n, hn using Nat.le_induction with
      | base => simp
      | succ n hn ih =>
        have := hmono k n hn
        push_cast
        nlinarith
    obtain ⟨N, hN⟩ := exists_nat_gt ((1 + |g k|) / δ + k)
    obtain ⟨n, hn, hgn⟩ := hfreq 1 one_pos (max N k)
    have h1 := hlin n (le_trans (le_max_right _ _) hn)
    have h2 : (N : ℝ) ≤ n := by exact_mod_cast le_trans (le_max_left _ _) hn
    have h3 : (1 + |g k|) / δ < (n : ℝ) - k := by linarith
    have h4 : 1 + |g k| < ((n : ℝ) - k) * δ := by rwa [div_lt_iff₀ hδ0] at h3
    have h5 := (abs_lt.1 hgn).2
    have h6 := neg_abs_le (g k)
    linarith
  have hA : Antitone g := antitone_nat_of_succ_le hanti
  have hnn : ∀ n, 0 ≤ g n := by
    intro n
    by_contra hneg
    push Not at hneg
    obtain ⟨m, hm, hgm⟩ := hfreq (-g n) (by linarith) n
    have h1 := hA hm
    have h2 := (abs_lt.1 hgm).1
    linarith
  rw [Metric.tendsto_atTop]
  intro ε hε
  obtain ⟨m, _, hgm⟩ := hfreq ε hε 0
  refine ⟨m, fun n hn => ?_⟩
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (hnn n)]
  exact lt_of_le_of_lt (hA hn) (lt_of_le_of_lt (le_abs_self _) hgm)

end ChebAllx
