/-
ChebAllxLib/Main.lean — Chebyshev's 1853 assertion at EVERY x > 0, and for every convex weight,
under GRH for `L(s, χ₋₄)`.

Inputs:
* `Allx.theta1_wlog_pos` (lean-allx): under `GRHFor chi4`,
  `Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p) > 0` for every real `t > 3`; used at the integers
  `t = 4, 5, 6, …` only;
* `ChebAllx.tsum_mul_pos_of_secondSum_pos` (`Discrete.lean`): second-order summation by parts over
  the infinite range, `Σ_n a(n) f(n) = Σ_k A₁(k+1) Δ²f(k)`.

With `a(j) = [j+3 prime] (−χ₋₄(j+3)) log(j+3)` one has `A₁(k+1) = Θ₁(k+4)`, so for every convex
sequence `f` on `n ≥ 3` with `Σ_p log p · |f(p)| < ∞`:

    −Σ_p χ₋₄(p) log p · f(p) = Σ_{k ≥ 0} Θ₁(k+4) · (f(k+3) − 2 f(k+4) + f(k+5)),

a series of nonnegative terms.  For the sum without logarithms take `f(n) = w(n)/log n`, which is
convex when `w` is (the product of two nonnegative nonincreasing convex sequences is convex).
-/
import ChebAllxLib.Discrete
import AllxLib.Final
import ChebLib

open Finset Filter Topology

noncomputable section

namespace ChebAllx

/-- `a(j) = [j+3 prime] (−χ₋₄(j+3)) log(j+3)`: the coefficients of the θ-race, shifted by 3. -/
def thetaCoeff (j : ℕ) : ℝ :=
  if (j + 3).Prime then Coupling.w4 (j + 3) * Real.log ((j + 3 : ℕ) : ℝ) else 0

lemma w4_two : Coupling.w4 2 = 0 := by
  unfold Coupling.w4
  have : ZMod.χ₄ ((2 : ℕ) : ZMod 4) = 0 := by decide
  rw [this]; simp

lemma abs_w4_le_one (n : ℕ) : |Coupling.w4 n| ≤ 1 := by
  unfold Coupling.w4
  rw [abs_neg]
  have h : ∀ a : ZMod 4, ZMod.χ₄ a = 0 ∨ ZMod.χ₄ a = 1 ∨ ZMod.χ₄ a = -1 := by decide
  rcases h (n : ZMod 4) with h | h | h <;> rw [h] <;> norm_num

lemma sum_range_shift3 (h : ℕ → ℝ) (h0 : h 0 = 0) (h1 : h 1 = 0) (h2 : h 2 = 0) (m : ℕ)
    (hm : h (m + 3) = 0) :
    ∑ p ∈ Finset.range (m + 3 + 1), h p = ∑ j ∈ Finset.range m, h (j + 3) := by
  rw [Finset.sum_range_succ, hm, add_zero, Finset.sum_range_succ' h (m + 2),
    Finset.sum_range_succ' (fun i => h (i + 1)) (m + 1),
    Finset.sum_range_succ' (fun i => h (i + 1 + 1)) m]
  simp only [zero_add, h0, h1, add_zero]
  rw [show (1 : ℕ) + 1 = 2 from rfl, h2, add_zero]

/-- The second sums of `thetaCoeff` are the values of `Θ₁` at the integers: `A₁(k+1) = Θ₁(k+4)`. -/
lemma secondSum_coeff (k : ℕ) :
    secondSum thetaCoeff (k + 1)
      = Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) ((k + 4 : ℕ) : ℝ) := by
  unfold Coupling.Theta1 secondSum
  rw [Nat.floor_natCast, Finset.sum_filter, ← Nat.range_succ_eq_Icc_zero,
    show k + 4 + 1 = (k + 1) + 3 + 1 by omega, sum_range_shift3 _ _ _ _ (k + 1)]
  · refine Finset.sum_congr rfl fun j _ => ?_
    unfold thetaCoeff
    split_ifs with hp
    · push_cast; ring
    · rw [zero_mul]
  · split_ifs
    · push_cast; ring
    · rfl
  · simp [Nat.not_prime_zero]
  · simp [Nat.not_prime_one]
  · simp [w4_two]

lemma secondSum_coeff_pos (hGRH : GRHFor chi4) (k : ℕ) : 0 < secondSum thetaCoeff (k + 1) := by
  rw [secondSum_coeff]
  refine Allx.theta1_wlog_pos hGRH ?_
  push_cast
  linarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)]

/-- From a convergent series over the primes: the terms are small at arbitrarily large primes. -/
lemma exists_prime_small {g : ℕ → ℝ} (hsum : Summable fun p : Nat.Primes => g p) {ε : ℝ}
    (hε : 0 < ε) (N : ℕ) : ∃ p, N ≤ p ∧ p.Prime ∧ |g p| < ε := by
  have h1 : Summable ({n : ℕ | n.Prime}.indicator g) :=
    (summable_subtype_iff_indicator (s := {n : ℕ | n.Prime}) (f := g)).1 hsum
  obtain ⟨N0, hN0⟩ := Metric.tendsto_atTop.1 h1.tendsto_atTop_zero ε hε
  obtain ⟨p, hpN, hp⟩ := Nat.exists_infinite_primes (max N N0)
  refine ⟨p, le_trans (le_max_left _ _) hpN, hp, ?_⟩
  have := hN0 p (le_trans (le_max_right _ _) hpN)
  rwa [Real.dist_eq, sub_zero, Set.indicator_of_mem (show p ∈ {n : ℕ | n.Prime} from hp)] at this

/-- Rewriting a sum over the primes weighted by `χ₋₄` as a sum over `j = p − 3 ≥ 0`. -/
lemma tsum_primes_chi4_eq (F : ℕ → ℝ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * F p
      = -∑' j : ℕ, (if (j + 3).Prime then Coupling.w4 (j + 3) * F (j + 3) else 0) := by
  rw [ChebLib.tsum_primes_eq (fun n => ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * F n), ← tsum_neg]
  have hinj : Function.Injective (fun j : ℕ => j + 3) := fun a b h => by simpa using h
  rw [← hinj.tsum_eq (f := fun n : ℕ => if n.Prime then ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) * F n
    else 0)]
  · refine tsum_congr fun j => ?_
    unfold Coupling.w4
    split_ifs
    · push_cast; ring
    · simp
  · intro n hn
    rw [Function.mem_support] at hn
    have h3 : 3 ≤ n := by
      by_contra hlt
      push Not at hlt
      apply hn
      interval_cases n
      · simp [Nat.not_prime_zero]
      · simp [Nat.not_prime_one]
      · have : ZMod.χ₄ ((2 : ℕ) : ZMod 4) = 0 := by decide
        simp
    exact ⟨n - 3, by simp; omega⟩

/-- The common set-up: for `f` convex on `n ≥ 3` with `Σ_p log p · |f(p)| < ∞`, the shifted
sequence `f(j+3)` is convex and tends to `0`, `Σ_j |a(j)| f(j+3) < ∞`, and
`Σ_p χ₋₄(p) log p · f(p) = −Σ_j a(j) f(j+3)`. -/
lemma log_weight_setup {f : ℕ → ℝ}
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) :
    (∀ k, 0 ≤ secondDiff (fun j => f (j + 3)) k) ∧ Tendsto (fun j => f (j + 3)) atTop (𝓝 0) ∧
      (Summable fun j => |thetaCoeff j| * f (j + 3)) ∧
      ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p
        = -∑' j : ℕ, thetaCoeff j * f (j + 3) := by
  set fs : ℕ → ℝ := fun j => f (j + 3) with hfs
  have hD : ∀ k, 0 ≤ secondDiff fs k := by
    intro k
    have h := hconv (k + 3) (by omega)
    unfold secondDiff
    simp only [hfs]
    rw [show k + 1 + 3 = k + 3 + 1 by omega, show k + 2 + 3 = k + 3 + 2 by omega]
    linarith
  have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have h0 : Tendsto fs atTop (𝓝 0) := by
    refine tendsto_zero_of_secondDiff_nonneg hD fun ε hε N => ?_
    obtain ⟨p, hpN, hp, hsmall⟩ := exists_prime_small (g := fun n => Real.log n * f n) hsum
      (mul_pos hε hlog2) (N + 3)
    refine ⟨p - 3, by omega, ?_⟩
    simp only [hfs]
    rw [show p - 3 + 3 = p by omega]
    have hlogp : Real.log 2 ≤ Real.log p :=
      Real.log_le_log (by norm_num) (by exact_mod_cast hp.two_le)
    rw [abs_mul, abs_of_nonneg (by linarith : (0 : ℝ) ≤ Real.log p)] at hsmall
    by_contra hge
    push Not at hge
    nlinarith [abs_nonneg (f p)]
  have hfs0 : ∀ j, 0 ≤ fs j := nonneg_of_secondDiff_nonneg hD h0
  -- absolute summability of `a(j) f(j+3)`
  have hsN : Summable (fun j : ℕ => ({n : ℕ | n.Prime}.indicator (fun n : ℕ => Real.log n * f n))
      (j + 3)) :=
    (summable_nat_add_iff 3).2
      ((summable_subtype_iff_indicator (s := {n : ℕ | n.Prime})
        (f := fun n : ℕ => Real.log n * f n)).1 hsum)
  have hs : Summable fun j => |thetaCoeff j| * fs j := by
    refine Summable.of_norm_bounded hsN.abs fun j => ?_
    rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (abs_nonneg _) (hfs0 j))]
    unfold thetaCoeff
    by_cases hp : (j + 3).Prime
    · rw [if_pos hp, Set.indicator_of_mem (show j + 3 ∈ {n : ℕ | n.Prime} from hp)]
      have hl : 0 ≤ Real.log ((j + 3 : ℕ) : ℝ) := Real.log_natCast_nonneg _
      rw [abs_mul, abs_mul, abs_of_nonneg hl, abs_of_nonneg (hfs0 j)]
      have := abs_w4_le_one (j + 3)
      have h2 : 0 ≤ Real.log ((j + 3 : ℕ) : ℝ) * fs j := mul_nonneg hl (hfs0 j)
      simp only [hfs] at h2 ⊢
      nlinarith [abs_nonneg (Coupling.w4 (j + 3))]
    · rw [if_neg hp, abs_zero, zero_mul]
      exact abs_nonneg _
  have e : ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p
      = -∑' j : ℕ, thetaCoeff j * fs j := by
    have e1 := tsum_primes_chi4_eq (fun n => Real.log n * f n)
    beta_reduce at e1
    have e2 : ∀ j : ℕ, (if (j + 3).Prime then
        Coupling.w4 (j + 3) * (Real.log ((j + 3 : ℕ) : ℝ) * f (j + 3)) else 0) = thetaCoeff j * fs j := by
      intro j
      unfold thetaCoeff
      split_ifs
      · simp only [hfs]; ring
      · rw [zero_mul]
    rw [← tsum_congr e2, ← e1]
    exact tsum_congr fun p => by ring
  exact ⟨hD, h0, hs, e⟩

/-- **The series identity (discrete form of Lemma 1 of the note; no hypothesis on zeros).**  For
every real sequence `f` convex on `n ≥ 3` with `Σ_p log p · |f(p)| < ∞`:

    −Σ_p χ₋₄(p) log p · f(p) = Σ_{k ≥ 0} Θ₁(k+4) · (f(k+3) − 2 f(k+4) + f(k+5)),

where `Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p)`; all the second differences are `≥ 0`. -/
theorem hasSum_theta1_mul_secondDiff {f : ℕ → ℝ}
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) :
    HasSum (fun k : ℕ => Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) ((k + 4 : ℕ) : ℝ)
        * (f (k + 3) - 2 * f (k + 4) + f (k + 5)))
      (-∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p) := by
  obtain ⟨hD, h0, hs, e⟩ := log_weight_setup hconv hsum
  have h := hasSum_secondSum_mul_secondDiff hD h0 hs
  rw [e, neg_neg]
  convert h using 1
  funext k
  rw [secondSum_coeff]
  unfold secondDiff
  rw [show k + 1 + 3 = k + 4 by omega, show k + 2 + 3 = k + 5 by omega]

/-- `Θ₁` at an integer, written out. -/
lemma theta1_natCast (m : ℕ) :
    Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) (m : ℝ)
      = ∑ p ∈ (Finset.Icc 0 m).filter Nat.Prime,
          -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * ((m : ℝ) - p) := by
  unfold Coupling.Theta1 Coupling.w4
  rw [Nat.floor_natCast]

/-- **Core theorem (logarithmic weights).**  Under GRH for `L(s, χ₋₄)`: for every real sequence `f`
which is convex on `n ≥ 3`, with `Σ_p log p · f(p)` absolutely convergent and `f(3) > 0`,
`Σ_p χ₋₄(p) log p · f(p) < 0`. -/
theorem tsum_log_mul_neg (hGRH : GRHFor chi4) {f : ℕ → ℝ}
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) (hf3 : 0 < f 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p < 0 := by
  obtain ⟨hD, h0, hs, e⟩ := log_weight_setup hconv hsum
  have hpos := tsum_mul_pos_of_secondSum_pos hD h0 hs (secondSum_coeff_pos hGRH) hf3
  rw [e]
  linarith

/-- The sequence `1/log(j+3)` is convex. -/
lemma inv_log_convex (k : ℕ) :
    2 * (1 / Real.log ((k + 3 + 1 : ℕ) : ℝ))
      ≤ 1 / Real.log ((k + 3 : ℕ) : ℝ) + 1 / Real.log ((k + 3 + 2 : ℕ) : ℝ) := by
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have ha : 0 < Real.log ((k + 3 : ℕ) : ℝ) := Real.log_pos (by push_cast; linarith)
  have hb : 0 < Real.log ((k + 3 + 2 : ℕ) : ℝ) := Real.log_pos (by push_cast; linarith)
  have hc : 0 < Real.log ((k + 3 + 1 : ℕ) : ℝ) := Real.log_pos (by push_cast; linarith)
  have habc : Real.log ((k + 3 : ℕ) : ℝ) + Real.log ((k + 3 + 2 : ℕ) : ℝ)
      ≤ 2 * Real.log ((k + 3 + 1 : ℕ) : ℝ) := by
    rw [← Real.log_mul (by push_cast; linarith) (by push_cast; linarith), ← Real.log_rpow
      (by push_cast; linarith)]
    refine Real.log_le_log (by push_cast; nlinarith) ?_
    rw [Real.rpow_two]
    push_cast
    nlinarith
  rw [show 2 * (1 / Real.log ((k + 3 + 1 : ℕ) : ℝ)) = 2 / Real.log ((k + 3 + 1 : ℕ) : ℝ) by ring,
    div_add_div _ _ ha.ne' hb.ne', div_le_div_iff₀ hc (mul_pos ha hb)]
  nlinarith [sq_nonneg (Real.log ((k + 3 : ℕ) : ℝ) - Real.log ((k + 3 + 2 : ℕ) : ℝ)),
    mul_pos ha hb, mul_nonneg (add_pos ha hb).le (sub_nonneg.2 habc)]

/-- **General theorem (convex weights).**  Under GRH for `L(s, χ₋₄)`: for every real sequence `w`
which is convex on `n ≥ 3` (`2 w(n+1) ≤ w(n) + w(n+2)`), with `Σ_p w(p)` (absolutely) convergent
and `w(3) > 0`:  `Σ_p χ₋₄(p) w(p) < 0`. -/
theorem tsum_mul_neg (hGRH : GRHFor chi4) {w : ℕ → ℝ}
    (hconv : ∀ n, 3 ≤ n → 2 * w (n + 1) ≤ w n + w (n + 2))
    (hsum : Summable fun p : Nat.Primes => w p) (hw3 : 0 < w 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w p < 0 := by
  -- `w(j+3)` is convex and tends to `0`, hence nonincreasing and nonnegative
  set ws : ℕ → ℝ := fun j => w (j + 3) with hws
  have hD : ∀ k, 0 ≤ secondDiff ws k := by
    intro k
    have h := hconv (k + 3) (by omega)
    unfold secondDiff
    simp only [hws]
    rw [show k + 1 + 3 = k + 3 + 1 by omega, show k + 2 + 3 = k + 3 + 2 by omega]
    linarith
  have h0 : Tendsto ws atTop (𝓝 0) := by
    refine tendsto_zero_of_secondDiff_nonneg hD fun ε hε N => ?_
    obtain ⟨p, hpN, hp, hsmall⟩ := exists_prime_small hsum hε (N + 3)
    refine ⟨p - 3, by omega, ?_⟩
    simp only [hws]
    rwa [show p - 3 + 3 = p by omega]
  have hws0 : ∀ j, 0 ≤ ws j := nonneg_of_secondDiff_nonneg hD h0
  have hws1 : ∀ j, ws (j + 1) ≤ ws j := fun j =>
    antitone_of_secondDiff_nonneg hD h0 (Nat.le_succ j)
  -- `v(j) = 1/log(j+3)` is nonnegative, nonincreasing and convex
  set v : ℕ → ℝ := fun j => 1 / Real.log ((j + 3 : ℕ) : ℝ) with hv
  have hlogpos : ∀ j : ℕ, 0 < Real.log ((j + 3 : ℕ) : ℝ) := fun j =>
    Real.log_pos (by push_cast; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)])
  have hv0 : ∀ j, 0 ≤ v j := fun j => (one_div_pos.2 (hlogpos j)).le
  have hv1 : ∀ j, v (j + 1) ≤ v j := by
    intro j
    simp only [hv]
    refine one_div_le_one_div_of_le (hlogpos j) (Real.log_le_log ?_ ?_)
    · push_cast; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
    · push_cast; linarith
  have hv2 : ∀ k, 0 ≤ secondDiff v k := by
    intro k
    have h := inv_log_convex k
    unfold secondDiff
    simp only [hv]
    rw [show k + 1 + 3 = k + 3 + 1 by omega, show k + 2 + 3 = k + 3 + 2 by omega]
    linarith
  -- apply the core theorem to `f(n) = w(n)/log n`
  have key := tsum_log_mul_neg hGRH (f := fun n => w n / Real.log n) ?_ ?_ ?_
  · have e : ∀ p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
        * (w p / Real.log (p : ℕ)) = ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w p := by
      intro p
      have : Real.log (p : ℕ) ≠ 0 :=
        (Real.log_pos (by exact_mod_cast p.2.one_lt)).ne'
      field_simp
    rwa [tsum_congr e] at key
  · intro n hn
    obtain ⟨k, rfl⟩ : ∃ k, n = k + 3 := ⟨n - 3, by omega⟩
    have h := secondDiff_mul_nonneg hws0 hv0 hws1 hv1 hD hv2 k
    unfold secondDiff at h
    simp only [hws, hv] at h
    rw [show k + 1 + 3 = k + 3 + 1 by omega, show k + 2 + 3 = k + 3 + 2 by omega] at h
    simp only [div_eq_mul_one_div (w _)]
    linarith
  · refine hsum.congr fun p => ?_
    have : Real.log (p : ℕ) ≠ 0 :=
      (Real.log_pos (by exact_mod_cast p.2.one_lt)).ne'
    field_simp
  · exact div_pos hw3 (Real.log_pos (by norm_num))

/-- The same for a weight given as a real function, convex on `[3, ∞)`. -/
theorem tsum_mul_neg_of_convexOn (hGRH : GRHFor chi4) {w : ℝ → ℝ}
    (hconv : ConvexOn ℝ (Set.Ici 3) w)
    (hsum : Summable fun p : Nat.Primes => w (p : ℕ)) (hw3 : 0 < w 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w (p : ℕ) < 0 := by
  refine tsum_mul_neg hGRH (w := fun n => w n) ?_ hsum (by simpa using hw3)
  intro n hn
  have hn3 : (3 : ℝ) ≤ n := by exact_mod_cast hn
  have key := hconv.2 (Set.mem_Ici.2 hn3) (Set.mem_Ici.2 (by linarith : (3 : ℝ) ≤ (n : ℝ) + 2))
    (by norm_num : (0 : ℝ) ≤ 1 / 2) (by norm_num : (0 : ℝ) ≤ 1 / 2) (by norm_num)
  have e : (1 / 2 : ℝ) • (n : ℝ) + (1 / 2 : ℝ) • ((n : ℝ) + 2) = n + 1 := by
    simp only [smul_eq_mul]; ring
  rw [e] at key
  simp only [smul_eq_mul] at key
  push_cast
  linarith

/-- The same for the logarithmically weighted sum and a real function convex on `[3, ∞)`. -/
theorem tsum_log_mul_neg_of_convexOn (hGRH : GRHFor chi4) {f : ℝ → ℝ}
    (hconv : ConvexOn ℝ (Set.Ici 3) f)
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f (p : ℕ)) (hf3 : 0 < f 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f (p : ℕ)
      < 0 := by
  refine tsum_log_mul_neg hGRH (f := fun n => f n) ?_ hsum (by simpa using hf3)
  intro n hn
  have hn3 : (3 : ℝ) ≤ n := by exact_mod_cast hn
  have key := hconv.2 (Set.mem_Ici.2 hn3) (Set.mem_Ici.2 (by linarith : (3 : ℝ) ≤ (n : ℝ) + 2))
    (by norm_num : (0 : ℝ) ≤ 1 / 2) (by norm_num : (0 : ℝ) ≤ 1 / 2) (by norm_num)
  have e : (1 / 2 : ℝ) • (n : ℝ) + (1 / 2 : ℝ) • ((n : ℝ) + 2) = n + 1 := by
    simp only [smul_eq_mul]; ring
  rw [e] at key
  simp only [smul_eq_mul] at key
  push_cast
  linarith

/-! ### Chebyshev's sum -/

lemma convexOn_exp_neg_div (x : ℝ) : ConvexOn ℝ (Set.Ici 3) (fun t : ℝ => Real.exp (-t / x)) := by
  have h : ConvexOn ℝ Set.univ (fun t : ℝ => Real.exp (-t / x)) := by
    have := convexOn_exp.comp_linearMap ((-1 / x) • (LinearMap.id : ℝ →ₗ[ℝ] ℝ))
    simp only [Set.preimage_univ] at this
    refine this.congr ?_ |>.subset (Set.subset_univ _) convex_univ
    intro t _
    simp only [Function.comp_apply, LinearMap.smul_apply, LinearMap.id_coe, id_eq, smul_eq_mul]
    congr 1; ring
  exact h.subset (Set.subset_univ _) (convex_Ici 3)

lemma summable_primes_exp_neg_div {x : ℝ} (hx : 0 < x) :
    Summable fun p : Nat.Primes => Real.exp (-((p : ℕ) : ℝ) / x) := by
  have h : Summable fun n : ℕ => Real.exp (n * (-1 / x)) :=
    Real.summable_exp_nat_mul_iff.2 (by
      have : 0 < 1 / x := one_div_pos.2 hx
      rw [neg_div]; linarith)
  have h' : Summable fun n : ℕ => Real.exp (-(n : ℝ) / x) :=
    h.congr fun n => by congr 1; ring
  exact h'.subtype {n : ℕ | n.Prime}

lemma summable_primes_log_mul_exp_neg_div {x : ℝ} (hx : 0 < x) :
    Summable fun p : Nat.Primes => Real.log (p : ℕ) * Real.exp (-((p : ℕ) : ℝ) / x) := by
  have h := ChebLib.summable_mul_exp_neg_div (fun n : ℕ => Real.log n) 1 (fun n => by
    rw [one_mul, abs_of_nonneg (Real.log_natCast_nonneg n)]
    exact (Real.log_le_self (Nat.cast_nonneg n))) hx
  have h' : Summable fun n : ℕ => Real.log n * Real.exp (-(n : ℝ) / x) :=
    h.congr fun n => by rw [neg_div]
  exact h'.subtype {n : ℕ | n.Prime}

/-- **Chebyshev's 1853 assertion, pointwise.**  Under GRH for `L(s, χ₋₄)`:
`Σ_p χ₋₄(p) e^{−p/x} < 0` for EVERY real `x > 0`. -/
theorem chebyshevSum_neg (hGRH : GRHFor chi4) {x : ℝ} (hx : 0 < x) : chebyshevSum x < 0 :=
  tsum_mul_neg_of_convexOn hGRH (convexOn_exp_neg_div x) (summable_primes_exp_neg_div hx)
    (Real.exp_pos _)

/-- Under GRH for `L(s, χ₋₄)`: `Σ_p χ₋₄(p) log p · e^{−p/x} < 0` for EVERY real `x > 0`. -/
theorem chebyshevLogSum_neg (hGRH : GRHFor chi4) {x : ℝ} (hx : 0 < x) : chebyshevLogSum x < 0 :=
  tsum_log_mul_neg_of_convexOn hGRH (convexOn_exp_neg_div x)
    (summable_primes_log_mul_exp_neg_div hx) (Real.exp_pos _)

end ChebAllx
