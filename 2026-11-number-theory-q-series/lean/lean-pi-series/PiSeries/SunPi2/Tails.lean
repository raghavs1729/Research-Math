/-
# WP-T: tails of the rows `a(N,·)` (see `N2-PLAN.md`).

`ρ n k → ρ̄(x) = (4x²−1)²/(64(1−x²)²)` for `k = xn`, and `ρ̄(x) < 1 ⟺ x < x* = √3/2`.
With `R(y) = (4y²−1)/(1−y²)` (so `ρ̄ = R²/64`, `R(y) < 8 ⟺ y < x*`) we prove the explicit bounds

* `k ≤ yn`, `n ≥ 3`, `2/5 ≤ y² < 3/4`:  `ρ n k ≤ (1+4/n)(16R+408/n)²/16384`,
* `yn ≤ k < n`, `1/4 ≤ y² < 1`:          `ρ n k ≥ (16R)²/(1+1/n)²/16384`,

and deduce the geometric decay of any positive bounded `f` with `f (n+1) k = f n k * ρ n k`
outside the window `|k − x*N| ≤ δN`.
-/
import PiSeries.SunPi2.Defs

open Finset Filter Topology

namespace PiSeries.SunPi2

lemma xstar_sq : xstar ^ 2 = 3 / 4 := by
  unfold xstar
  rw [div_pow, Real.sq_sqrt (by norm_num)]
  norm_num

lemma xstar_nonneg : 0 ≤ xstar := by unfold xstar; positivity

lemma xstar_gt : (86 : ℝ) / 100 < xstar := by
  have h := xstar_sq
  have h0 := xstar_nonneg
  by_contra hc
  rw [not_lt] at hc
  nlinarith

lemma xstar_lt : xstar < 87 / 100 := by
  have h := xstar_sq
  have h0 := xstar_nonneg
  by_contra hc
  rw [not_lt] at hc
  nlinarith

/-- `ρ` as a product of the `n`-prefactor and the `k`-dependent part -/
lemma rho_split {m κ : ℝ} {n k : ℕ} (hmn : (n : ℝ) = m) (hkk : (k : ℝ) = κ) :
    ρ n k = (4 * m + 1) * (4 * m + 3) ^ 2 / (8 * m * (2 * m - 1) * (4 * m - 1))
      * ((8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * ((8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7))
          / ((m - κ) * (m + κ + 1)) ^ 2) / 16384 := by
  unfold ρ
  rw [hmn, hkk, div_mul_div_comm, div_div]
  congr 1 <;> ring

lemma pref_den_pos {m : ℝ} (hm : 1 ≤ m) : 0 < 8 * m * (2 * m - 1) * (4 * m - 1) := by
  have h1 : 0 < m := by linarith
  have h2 : 0 < 2 * m - 1 := by linarith
  have h3 : 0 < 4 * m - 1 := by linarith
  positivity

lemma pref_le {m : ℝ} (hm : 3 ≤ m) :
    (4 * m + 1) * (4 * m + 3) ^ 2 / (8 * m * (2 * m - 1) * (4 * m - 1)) ≤ 1 + 4 / m := by
  have hmpos : 0 < m := by linarith
  rw [div_le_iff₀ (pref_den_pos (by linarith))]
  have e : (1 + 4 / m) * (8 * m * (2 * m - 1) * (4 * m - 1))
      = (m + 4) * (8 * (2 * m - 1) * (4 * m - 1)) := by
    field_simp
  rw [e]
  nlinarith [mul_nonneg (sub_nonneg.2 hm) (sub_nonneg.2 hm)]

lemma pref_ge {m : ℝ} (hm : 1 ≤ m) :
    1 ≤ (4 * m + 1) * (4 * m + 3) ^ 2 / (8 * m * (2 * m - 1) * (4 * m - 1)) := by
  rw [le_div_iff₀ (pref_den_pos hm)]
  nlinarith [sq_nonneg m]

/-- the `k`-part, upper bound for `0 ≤ κ ≤ y m` -/
lemma kpart_le {y R m κ : ℝ} (hy1 : y ^ 2 < 3 / 4)
    (hR : R * (1 - y ^ 2) = 4 * y ^ 2 - 1) (hR1 : 1 ≤ R)
    (hm : 1 ≤ m) (hκ0 : 0 ≤ κ) (hk : κ ≤ y * m) :
    (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * ((8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7))
      / ((m - κ) * (m + κ + 1)) ^ 2 ≤ (16 * R + 408 / m) ^ 2 := by
  have hmpos : 0 < m := by linarith
  have hκsq : κ ^ 2 ≤ y ^ 2 * m ^ 2 := by
    have := pow_le_pow_left₀ hκ0 hk 2
    rwa [mul_pow] at this
  have hm2 : 0 < m ^ 2 := by positivity
  have hym : y ^ 2 * m ^ 2 ≤ 3 / 4 * m ^ 2 := mul_le_mul_of_nonneg_right hy1.le hm2.le
  have hD4 : m ^ 2 / 4 ≤ m ^ 2 - κ ^ 2 := by linarith
  have hD0 : 0 < m ^ 2 - κ ^ 2 := by linarith
  have hκm : κ < m := by
    by_contra hc
    rw [not_lt] at hc
    have := pow_le_pow_left₀ hmpos.le hc 2
    linarith
  have h4R : (4 + R) * y ^ 2 = R + 1 := by linear_combination -hR
  -- `|4κ² − m²| ≤ R (m² − κ²)`
  have hA1 : 4 * κ ^ 2 - m ^ 2 ≤ R * (m ^ 2 - κ ^ 2) := by
    have h := mul_le_mul_of_nonneg_left hκsq (show (0 : ℝ) ≤ 4 + R by linarith)
    have h' : (4 + R) * (y ^ 2 * m ^ 2) = (R + 1) * m ^ 2 := by rw [← mul_assoc, h4R]
    rw [h'] at h
    linarith
  have hA2 : m ^ 2 - 4 * κ ^ 2 ≤ R * (m ^ 2 - κ ^ 2) := by
    have := mul_nonneg (sub_nonneg.2 hR1) hD0.le
    have := sq_nonneg κ
    linarith
  set E : ℝ := (16 * R + 408 / m) * (m ^ 2 - κ ^ 2) with hE
  have hDm : m / 4 ≤ (m ^ 2 - κ ^ 2) / m := by
    rw [le_div_iff₀ hmpos]
    linarith
  have hE' : 16 * R * (m ^ 2 - κ ^ 2) + 102 * m ≤ E := by
    have : E = 16 * R * (m ^ 2 - κ ^ 2) + 408 * ((m ^ 2 - κ ^ 2) / m) := by
      rw [hE]; ring
    rw [this]; linarith
  have habs : |(8 * κ + 4 * m + 3) * (8 * κ - 4 * m + 6)| ≤ E := by
    rw [abs_le]
    constructor
    · have e : (8 * κ + 4 * m + 3) * (8 * κ - 4 * m + 6)
          = -(16 * (m ^ 2 - 4 * κ ^ 2)) + (72 * κ + 12 * m + 18) := by ring
      rw [e]; linarith
    · have e : (8 * κ + 4 * m + 3) * (8 * κ - 4 * m + 6)
          = 16 * (4 * κ ^ 2 - m ^ 2) + (72 * κ + 12 * m + 18) := by ring
      rw [e]; linarith
  have hN12 : (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * ((8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7))
      ≤ E ^ 2 := by
    have h1 : (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) ≤ (8 * κ + 4 * m + 3) ^ 2 := by
      have e : (8 * κ + 4 * m + 3) ^ 2 - (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3)
          = 2 * (8 * κ + 4 * m + 3) := by ring
      linarith
    have h1' : 0 ≤ (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) := by positivity
    have h2 : (8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7) ≤ (8 * κ - 4 * m + 6) ^ 2 := by
      have e : (8 * κ - 4 * m + 6) ^ 2 - (8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7) = 1 := by ring
      linarith
    calc (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * ((8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7))
        ≤ (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * (8 * κ - 4 * m + 6) ^ 2 :=
          mul_le_mul_of_nonneg_left h2 h1'
      _ ≤ (8 * κ + 4 * m + 3) ^ 2 * (8 * κ - 4 * m + 6) ^ 2 :=
          mul_le_mul_of_nonneg_right h1 (sq_nonneg _)
      _ = ((8 * κ + 4 * m + 3) * (8 * κ - 4 * m + 6)) ^ 2 := by ring
      _ ≤ E ^ 2 := sq_le_sq' (abs_le.1 habs).1 (abs_le.1 habs).2
  have hden2 : m ^ 2 - κ ^ 2 ≤ (m - κ) * (m + κ + 1) := by
    have e : (m - κ) * (m + κ + 1) = m ^ 2 - κ ^ 2 + (m - κ) := by ring
    linarith
  have hden2' : 0 < ((m - κ) * (m + κ + 1)) ^ 2 := by
    have : 0 < (m - κ) * (m + κ + 1) := lt_of_lt_of_le hD0 hden2
    positivity
  rw [div_le_iff₀ hden2']
  calc _ ≤ E ^ 2 := hN12
    _ = (16 * R + 408 / m) ^ 2 * (m ^ 2 - κ ^ 2) ^ 2 := by rw [hE]; ring
    _ ≤ (16 * R + 408 / m) ^ 2 * ((m - κ) * (m + κ + 1)) ^ 2 :=
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hD0.le hden2 2) (sq_nonneg _)

/-- explicit upper bound for `ρ n k`, `k ≤ y n` -/
lemma rho_below_bound {y R : ℝ} (hy1 : y ^ 2 < 3 / 4)
    (hR : R * (1 - y ^ 2) = 4 * y ^ 2 - 1) (hR1 : 1 ≤ R)
    {n k : ℕ} (hn : 3 ≤ n) (hk : (k : ℝ) ≤ y * n) :
    ρ n k ≤ (1 + 4 / (n : ℝ)) * (16 * R + 408 / (n : ℝ)) ^ 2 / 16384 := by
  have hm : (3 : ℝ) ≤ n := by exact_mod_cast hn
  rw [rho_split rfl rfl]
  have hX1 := pref_le hm
  have hX2 := kpart_le hy1 hR hR1 (by linarith) (Nat.cast_nonneg k) hk
  have hX1pos : 0 ≤ (4 * (n : ℝ) + 1) * (4 * (n : ℝ) + 3) ^ 2
      / (8 * (n : ℝ) * (2 * (n : ℝ) - 1) * (4 * (n : ℝ) - 1)) :=
    le_trans zero_le_one (pref_ge (by linarith))
  exact div_le_div_of_nonneg_right
    ((mul_le_mul_of_nonneg_left hX2 hX1pos).trans (mul_le_mul_of_nonneg_right hX1 (sq_nonneg _)))
    (by norm_num)

/-- the `k`-part, lower bound for `y m ≤ κ < m` -/
lemma kpart_ge {y R m κ : ℝ} (hy0 : 1 / 2 ≤ y)
    (hR : R * (1 - y ^ 2) = 4 * y ^ 2 - 1) (hR1 : 0 ≤ R)
    (hm : 1 ≤ m) (hk : y * m ≤ κ) (hκm : κ < m) :
    (16 * R) ^ 2 / (1 + 1 / m) ^ 2
      ≤ (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * ((8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7))
        / ((m - κ) * (m + κ + 1)) ^ 2 := by
  have hmpos : 0 < m := by linarith
  have hym : 0 ≤ y * m := by positivity
  have hκ0 : 0 ≤ κ := hym.trans hk
  have hκsq : y ^ 2 * m ^ 2 ≤ κ ^ 2 := by
    have := pow_le_pow_left₀ hym hk 2
    rwa [mul_pow] at this
  have hD0 : 0 < m ^ 2 - κ ^ 2 := by
    have := pow_lt_pow_left₀ hκm hκ0 (two_ne_zero)
    linarith
  have h4R : (4 + R) * y ^ 2 = R + 1 := by linear_combination -hR
  have hA1 : R * (m ^ 2 - κ ^ 2) ≤ 4 * κ ^ 2 - m ^ 2 := by
    have h := mul_le_mul_of_nonneg_left hκsq (show (0 : ℝ) ≤ 4 + R by linarith)
    have h' : (4 + R) * (y ^ 2 * m ^ 2) = (R + 1) * m ^ 2 := by rw [← mul_assoc, h4R]
    rw [h'] at h
    linarith
  have hRD : 0 ≤ R * (m ^ 2 - κ ^ 2) := mul_nonneg hR1 hD0.le
  have hb : 0 ≤ 8 * κ - 4 * m := by
    have := mul_le_mul_of_nonneg_right hy0 hmpos.le
    linarith
  have hN12 : (16 * (R * (m ^ 2 - κ ^ 2))) ^ 2
      ≤ (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * ((8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7)) := by
    have h1 : (8 * κ + 4 * m) ^ 2 ≤ (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) := by
      have e : (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) - (8 * κ + 4 * m) ^ 2
          = 4 * (8 * κ + 4 * m) + 3 := by ring
      linarith
    have h2 : (8 * κ - 4 * m) ^ 2 ≤ (8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7) := by
      have e : (8 * κ - 4 * m + 5) * (8 * κ - 4 * m + 7) - (8 * κ - 4 * m) ^ 2
          = 12 * (8 * κ - 4 * m) + 35 := by ring
      linarith
    calc (16 * (R * (m ^ 2 - κ ^ 2))) ^ 2 ≤ (16 * (4 * κ ^ 2 - m ^ 2)) ^ 2 :=
          pow_le_pow_left₀ (by positivity) (by linarith) 2
      _ = (8 * κ + 4 * m) ^ 2 * (8 * κ - 4 * m) ^ 2 := by ring
      _ ≤ (8 * κ + 4 * m + 1) * (8 * κ + 4 * m + 3) * (8 * κ - 4 * m) ^ 2 :=
          mul_le_mul_of_nonneg_right h1 (sq_nonneg _)
      _ ≤ _ := mul_le_mul_of_nonneg_left h2 (by positivity)
  have hmκ : 0 < m - κ := by linarith
  have hden1 : 0 < (m - κ) * (m + κ + 1) := by positivity
  have hden2 : (m - κ) * (m + κ + 1) ≤ (m ^ 2 - κ ^ 2) * (1 + 1 / m) := by
    have h : m - κ ≤ (m ^ 2 - κ ^ 2) * (1 / m) := by
      rw [mul_one_div, le_div_iff₀ hmpos]
      have e : m ^ 2 - κ ^ 2 - (m - κ) * m = (m - κ) * κ := by ring
      have := mul_nonneg hmκ.le hκ0
      linarith
    have e : (m - κ) * (m + κ + 1) = m ^ 2 - κ ^ 2 + (m - κ) := by ring
    rw [e, mul_add, mul_one]
    linarith
  have h11 : 0 < (1 + 1 / m) ^ 2 := by positivity
  rw [div_le_div_iff₀ h11 (by positivity)]
  calc (16 * R) ^ 2 * ((m - κ) * (m + κ + 1)) ^ 2
      ≤ (16 * R) ^ 2 * ((m ^ 2 - κ ^ 2) * (1 + 1 / m)) ^ 2 :=
        mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hden1.le hden2 2) (sq_nonneg _)
    _ = (16 * (R * (m ^ 2 - κ ^ 2))) ^ 2 * (1 + 1 / m) ^ 2 := by ring
    _ ≤ _ := mul_le_mul_of_nonneg_right hN12 h11.le

/-- explicit lower bound for `ρ n k`, `y n ≤ k < n` -/
lemma rho_above_bound {y R : ℝ} (hy0 : 1 / 2 ≤ y)
    (hR : R * (1 - y ^ 2) = 4 * y ^ 2 - 1) (hR1 : 0 ≤ R)
    {n k : ℕ} (hn : 1 ≤ n) (hk : y * n ≤ (k : ℝ)) (hkn : k < n) :
    (16 * R) ^ 2 / (1 + 1 / (n : ℝ)) ^ 2 / 16384 ≤ ρ n k := by
  have hm : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hκm : (k : ℝ) < n := by exact_mod_cast hkn
  rw [rho_split rfl rfl]
  have hX1 := pref_ge hm
  have hX2 := kpart_ge hy0 hR hR1 hm hk hκm
  have hX2pos : 0 ≤ (16 * R) ^ 2 / (1 + 1 / (n : ℝ)) ^ 2 := by positivity
  have := mul_le_mul hX1 hX2 hX2pos (le_trans zero_le_one hX1)
  rw [one_mul] at this
  exact div_le_div_of_nonneg_right this (by norm_num)

theorem rho_le_below (η : ℝ) (hη : 0 < η) (hη' : η < 1 / 10) :
    ∃ q : ℝ, 0 < q ∧ q < 1 ∧ ∃ n0 : ℕ, ∀ n ≥ n0, ∀ k : ℕ, (k : ℝ) ≤ (xstar - η) * n → ρ n k ≤ q := by
  have hx1 := xstar_gt
  have hx2 := xstar_lt
  have hxs := xstar_sq
  obtain ⟨y, hy⟩ : ∃ y : ℝ, y = xstar - η := ⟨_, rfl⟩
  rw [← hy]
  have hy0 : 0 < y := by linarith
  have hy76 : (76 : ℝ) / 100 < y := by linarith
  have hylt : y < xstar := by linarith
  have hysq : y ^ 2 < 3 / 4 := by
    have := pow_lt_pow_left₀ hylt hy0.le two_ne_zero
    linarith
  have hysq' : (2 : ℝ) / 5 ≤ y ^ 2 := by
    have := pow_le_pow_left₀ (by norm_num) hy76.le 2
    linarith
  have h1y : 0 < 1 - y ^ 2 := by linarith
  obtain ⟨R, hRdef⟩ : ∃ R : ℝ, R = (4 * y ^ 2 - 1) / (1 - y ^ 2) := ⟨_, rfl⟩
  have hR : R * (1 - y ^ 2) = 4 * y ^ 2 - 1 := by rw [hRdef]; exact div_mul_cancel₀ _ h1y.ne'
  have hR1 : 1 ≤ R := by rw [hRdef, le_div_iff₀ h1y]; linarith
  have hR8 : R < 8 := by rw [hRdef, div_lt_iff₀ h1y]; linarith
  have hlim : Tendsto (fun n : ℕ => (1 + 4 / (n : ℝ)) * (16 * R + 408 / (n : ℝ)) ^ 2 / 16384) atTop
      (𝓝 ((1 + 0) * (16 * R + 0) ^ 2 / 16384)) :=
    ((tendsto_const_nhds.add (tendsto_const_div_atTop_nhds_zero_nat 4)).mul
      ((tendsto_const_nhds.add (tendsto_const_div_atTop_nhds_zero_nat 408)).pow 2)).div_const _
  have hval : (1 + 0) * (16 * R + 0) ^ 2 / 16384 = R ^ 2 / 64 := by ring
  rw [hval] at hlim
  have hR2 : R ^ 2 / 64 < 1 := by
    have := pow_lt_pow_left₀ hR8 (by linarith) two_ne_zero
    linarith
  have hq : R ^ 2 / 64 < (1 + R ^ 2 / 64) / 2 := by linarith
  obtain ⟨n0, hn0⟩ := eventually_atTop.1 (hlim.eventually_lt_const hq)
  refine ⟨(1 + R ^ 2 / 64) / 2, by positivity, by linarith, max n0 3, ?_⟩
  intro n hn k hk
  have h3 : 3 ≤ n := (le_max_right _ _).trans hn
  have hn0' : n0 ≤ n := (le_max_left _ _).trans hn
  exact (rho_below_bound hysq hR hR1 h3 hk).trans (hn0 n hn0').le

theorem rho_ge_above (η : ℝ) (hη : 0 < η) (hη' : η < 1 / 10) :
    ∃ q : ℝ, 1 < q ∧ ∃ n0 : ℕ, ∀ n ≥ n0, ∀ k : ℕ, (xstar + η) * n ≤ k → k < n → q ≤ ρ n k := by
  have hx1 := xstar_gt
  have hx2 := xstar_lt
  have hxs := xstar_sq
  obtain ⟨y, hy⟩ : ∃ y : ℝ, y = xstar + η := ⟨_, rfl⟩
  rw [← hy]
  have hy0 : 1 / 2 ≤ y := by linarith
  have hy1 : y < 1 := by linarith
  have hylt : xstar < y := by linarith
  have hysq : 3 / 4 < y ^ 2 := by
    have := pow_lt_pow_left₀ hylt xstar_nonneg two_ne_zero
    linarith
  have hysq' : y ^ 2 < 1 := by
    have := pow_lt_pow_left₀ hy1 (by linarith) two_ne_zero
    linarith
  have h1y : 0 < 1 - y ^ 2 := by linarith
  obtain ⟨R, hRdef⟩ : ∃ R : ℝ, R = (4 * y ^ 2 - 1) / (1 - y ^ 2) := ⟨_, rfl⟩
  have hR : R * (1 - y ^ 2) = 4 * y ^ 2 - 1 := by rw [hRdef]; exact div_mul_cancel₀ _ h1y.ne'
  have hR8 : 8 < R := by rw [hRdef, lt_div_iff₀ h1y]; linarith
  have hlim : Tendsto (fun n : ℕ => (16 * R) ^ 2 / (1 + 1 / (n : ℝ)) ^ 2 / 16384) atTop
      (𝓝 ((16 * R) ^ 2 / (1 + 0) ^ 2 / 16384)) :=
    (Tendsto.div tendsto_const_nhds
      ((tendsto_const_nhds.add (tendsto_const_div_atTop_nhds_zero_nat 1)).pow 2)
      (by norm_num)).div_const _
  have hval : (16 * R) ^ 2 / (1 + 0) ^ 2 / 16384 = R ^ 2 / 64 := by ring
  rw [hval] at hlim
  have hR2 : 1 < R ^ 2 / 64 := by
    have := pow_lt_pow_left₀ hR8 (by norm_num) two_ne_zero
    linarith
  have hq : (1 + R ^ 2 / 64) / 2 < R ^ 2 / 64 := by linarith
  obtain ⟨n0, hn0⟩ := eventually_atTop.1 (hlim.eventually_const_lt hq)
  refine ⟨(1 + R ^ 2 / 64) / 2, by linarith, max n0 1, ?_⟩
  intro n hn k hk hkn
  have h1 : 1 ≤ n := (le_max_right _ _).trans hn
  have hn0' : n0 ≤ n := (le_max_left _ _).trans hn
  exact (hn0 n hn0').le.trans (rho_above_bound hy0 hR (by linarith) h1 hk hkn)

/-- the column recurrence, multiplied up -/
lemma f_prod (f : ℕ → ℕ → ℝ)
    (hρ : ∀ n k : ℕ, 1 ≤ n → k < n → f (n + 1) k = f n k * ρ n k)
    {n1 N k : ℕ} (h1 : 1 ≤ n1) (hk : k < n1) (hN : n1 ≤ N) :
    f N k = f n1 k * ∏ n ∈ Ico n1 N, ρ n k := by
  induction N, hN using Nat.le_induction with
  | base => simp
  | succ N hN ih =>
    rw [hρ N k (h1.trans hN) (hk.trans_le hN), ih, Finset.prod_Ico_succ_top hN, mul_assoc]

lemma pow_tail {r δ : ℝ} (hr0 : 0 < r) (hr1 : r < 1) {m N : ℕ} (h : δ / 4 * N ≤ m) :
    r ^ m ≤ (r ^ (δ / 4)) ^ N := by
  calc r ^ m = r ^ (m : ℝ) := (Real.rpow_natCast r m).symm
    _ ≤ r ^ (δ / 4 * N) := Real.rpow_le_rpow_of_exponent_ge hr0 hr1.le h
    _ = (r ^ (δ / 4)) ^ (N : ℝ) := Real.rpow_mul hr0.le _ _
    _ = (r ^ (δ / 4)) ^ N := Real.rpow_natCast _ N

/-- upper tail: push the column up to `N' = ⌊k/y⌋₊`, where `ρ ≥ q1 > 1` on the way -/
lemma tail_upper (f : ℕ → ℕ → ℝ) (B : ℝ)
    (hpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < f n k) (hB : ∀ n k : ℕ, 1 ≤ n → k < n → f n k ≤ B)
    (hρ : ∀ n k : ℕ, 1 ≤ n → k < n → f (n + 1) k = f n k * ρ n k)
    {y q1 δ : ℝ} {m0 : ℕ} (hq1 : 1 < q1) (hy : 0 < y)
    (h1 : ∀ n ≥ m0, ∀ k : ℕ, y * n ≤ k → k < n → q1 ≤ ρ n k)
    {N k : ℕ} (hN2 : m0 ≤ N) (hkN : k < N) (hkN' : (1 + δ / 2) * N ≤ (k : ℝ) / y) (hδ : 0 ≤ δ)
    (hδN : 1 ≤ δ / 4 * N) :
    ∃ m : ℕ, δ / 4 * N ≤ m ∧ f N k ≤ B * (q1⁻¹) ^ m := by
  have hq1pos : 0 < q1 := by linarith
  have hNpos : 1 ≤ N := by omega
  have hNr : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  have hdiv0 : 0 ≤ (k : ℝ) / y := by positivity
  obtain ⟨N', hN'def⟩ : ∃ N' : ℕ, N' = ⌊(k : ℝ) / y⌋₊ := ⟨_, rfl⟩
  have hδN0 : 0 ≤ δ / 2 * N := by positivity
  have hNN' : N ≤ N' := by
    rw [hN'def]
    exact Nat.le_floor (by linarith)
  have hN'le : (N' : ℝ) ≤ k / y := by rw [hN'def]; exact Nat.floor_le hdiv0
  have hN'gt : (k : ℝ) / y < N' + 1 := by rw [hN'def]; exact Nat.lt_floor_add_one _
  have hprod : q1 ^ (N' - N) ≤ ∏ n ∈ Ico N N', ρ n k := by
    have : ∏ _n ∈ Ico N N', q1 ≤ ∏ n ∈ Ico N N', ρ n k := by
      apply Finset.prod_le_prod
      · intro _ _; exact hq1pos.le
      · intro n hn
        rw [Finset.mem_Ico] at hn
        apply h1 n (hN2.trans hn.1) k _ (lt_of_lt_of_le hkN hn.1)
        have hnN' : (n : ℝ) ≤ N' := by exact_mod_cast hn.2.le
        have := hnN'.trans hN'le
        rw [le_div_iff₀ hy] at this
        linarith
    simpa using this
  have hfN' : f N' k = f N k * ∏ n ∈ Ico N N', ρ n k := f_prod f hρ hNpos hkN hNN'
  have hfB : f N' k ≤ B := hB N' k (hNpos.trans hNN') (lt_of_lt_of_le hkN hNN')
  have hfpos := hpos N k hNpos hkN
  have hle : f N k * q1 ^ (N' - N) ≤ B :=
    (mul_le_mul_of_nonneg_left hprod hfpos.le).trans (hfN' ▸ hfB)
  refine ⟨N' - N, ?_, ?_⟩
  · rw [Nat.cast_sub hNN']
    linarith
  · have hqpow : 0 < q1 ^ (N' - N) := pow_pos hq1pos _
    calc f N k ≤ B / q1 ^ (N' - N) := by rwa [le_div_iff₀ hqpow]
      _ = B * (q1⁻¹) ^ (N' - N) := by rw [inv_pow, div_eq_mul_inv]

/-- lower tail: start the column at `n1 ≈ k/y`, from where on `ρ ≤ q0 < 1` -/
lemma tail_lower (f : ℕ → ℕ → ℝ) (B : ℝ) (hB0 : 0 ≤ B)
    (hB : ∀ n k : ℕ, 1 ≤ n → k < n → f n k ≤ B)
    (hρ : ∀ n k : ℕ, 1 ≤ n → k < n → f (n + 1) k = f n k * ρ n k)
    (hρpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < ρ n k)
    {y q0 δ : ℝ} {n0 : ℕ} (hy : 0 < y) (hy1 : y < 1)
    (h0 : ∀ n ≥ n0, ∀ k : ℕ, (k : ℝ) ≤ y * n → ρ n k ≤ q0)
    {N k : ℕ} (hN1 : 2 * (n0 + 1) ≤ N) (hkN' : (k : ℝ) / y ≤ (1 - δ / 2) * N) (hδ : 0 ≤ δ)
    (hδ1 : δ ≤ 1) (hδN : 1 ≤ δ / 4 * N) :
    ∃ m : ℕ, δ / 4 * N ≤ m ∧ f N k ≤ B * q0 ^ m := by
  have hNr : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  have hN1' : (2 : ℝ) * (n0 + 1) ≤ N := by exact_mod_cast hN1
  have hdiv0 : 0 ≤ (k : ℝ) / y := by positivity
  have hδN0 : 0 ≤ δ / 2 * N := by positivity
  have hδN1 : δ * N ≤ 1 * N := mul_le_mul_of_nonneg_right hδ1 hNr
  obtain ⟨n1, hn1def⟩ : ∃ n1 : ℕ, n1 = max (n0 + 1) ⌈(k : ℝ) / y⌉₊ := ⟨_, rfl⟩
  have hceil : ⌈(k : ℝ) / y⌉₊ ≤ N := Nat.ceil_le.2 (by linarith)
  have hn1N : n1 ≤ N := by rw [hn1def]; exact max_le (by omega) hceil
  have hn1pos : 1 ≤ n1 := by rw [hn1def]; exact le_trans (by omega) (le_max_left _ _)
  have hn0n1 : n0 ≤ n1 := by rw [hn1def]; exact le_trans (by omega) (le_max_left _ _)
  have hn1r : (1 : ℝ) ≤ n1 := by exact_mod_cast hn1pos
  have hkn1 : (k : ℝ) ≤ y * n1 := by
    have h : (k : ℝ) / y ≤ n1 := by
      rw [hn1def]
      exact (Nat.le_ceil _).trans (by exact_mod_cast le_max_right (n0 + 1) ⌈(k : ℝ) / y⌉₊)
    rw [div_le_iff₀ hy] at h
    linarith
  have hkltn1 : k < n1 := by
    have h : y * n1 < 1 * n1 := mul_lt_mul_of_pos_right hy1 (by linarith)
    have : (k : ℝ) < n1 := by linarith
    exact_mod_cast this
  have hρn : ∀ n ∈ Ico n1 N, 0 ≤ ρ n k := by
    intro n hn
    rw [Finset.mem_Ico] at hn
    exact (hρpos n k (hn1pos.trans hn.1) (lt_of_lt_of_le hkltn1 hn.1)).le
  have hprod : ∏ n ∈ Ico n1 N, ρ n k ≤ q0 ^ (N - n1) := by
    have : ∏ n ∈ Ico n1 N, ρ n k ≤ ∏ _n ∈ Ico n1 N, q0 := by
      apply Finset.prod_le_prod hρn
      intro n hn
      rw [Finset.mem_Ico] at hn
      apply h0 n (hn0n1.trans hn.1) k
      have hnr : (n1 : ℝ) ≤ n := by exact_mod_cast hn.1
      have := mul_le_mul_of_nonneg_left hnr hy.le
      linarith
    simpa using this
  refine ⟨N - n1, ?_, ?_⟩
  · rw [Nat.cast_sub hn1N]
    rcases max_choice (n0 + 1) ⌈(k : ℝ) / y⌉₊ with h | h
    · have : (n1 : ℝ) = n0 + 1 := by rw [hn1def, h]; push_cast; ring
      linarith
    · have h' : (n1 : ℝ) < k / y + 1 := by
        rw [hn1def, h]
        exact Nat.ceil_lt_add_one hdiv0
      linarith
  · rw [f_prod f hρ hn1pos hkltn1 hn1N]
    exact mul_le_mul (hB n1 k hn1pos hkltn1) hprod (Finset.prod_nonneg hρn) hB0

theorem tail_bound (f : ℕ → ℕ → ℝ) (B : ℝ)
    (hpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < f n k) (hB : ∀ n k : ℕ, 1 ≤ n → k < n → f n k ≤ B)
    (hρ : ∀ n k : ℕ, 1 ≤ n → k < n → f (n + 1) k = f n k * ρ n k)
    (δ : ℝ) (hδ : 0 < δ) (hδ' : δ < 1 / 10) :
    ∃ q : ℝ, 0 < q ∧ q < 1 ∧ ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ, k < N → δ * N ≤ |(k : ℝ) - xstar * N| →
      f N k ≤ B * q ^ N := by
  obtain ⟨q0, hq0, hq0', n0, h0⟩ := rho_le_below (δ / 2) (by linarith) (by linarith)
  obtain ⟨q1, hq1, m0, h1⟩ := rho_ge_above (δ / 2) (by linarith) (by linarith)
  have hx1 := xstar_gt
  have hx2 := xstar_lt
  have hBpos : 0 < B := (hpos 1 0 le_rfl Nat.one_pos).trans_le (hB 1 0 le_rfl Nat.one_pos)
  have hρpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < ρ n k := by
    intro n k hn hk
    have h := hpos (n + 1) k (by omega) (by omega)
    rw [hρ n k hn hk] at h
    exact (mul_pos_iff_of_pos_left (hpos n k hn hk)).1 h
  have hq1pos : 0 < q1 := by linarith
  obtain ⟨r, hr⟩ : ∃ r : ℝ, r = max q0 q1⁻¹ := ⟨_, rfl⟩
  have hr0 : 0 < r := by rw [hr]; exact lt_max_of_lt_left hq0
  have hr1 : r < 1 := by rw [hr]; exact max_lt hq0' (inv_lt_one_of_one_lt₀ hq1)
  have hrq0 : q0 ≤ r := by rw [hr]; exact le_max_left _ _
  have hrq1 : q1⁻¹ ≤ r := by rw [hr]; exact le_max_right _ _
  refine ⟨r ^ (δ / 4), Real.rpow_pos_of_pos hr0 _, Real.rpow_lt_one hr0.le hr1 (by positivity),
    max (max (2 * (n0 + 1)) m0) ⌈4 / δ⌉₊, ?_⟩
  intro N hN k hkN hk
  have hN1 : 2 * (n0 + 1) ≤ N := ((le_max_left _ _).trans (le_max_left _ _)).trans hN
  have hN2 : m0 ≤ N := ((le_max_right _ _).trans (le_max_left _ _)).trans hN
  have hN3 : 4 / δ ≤ (N : ℝ) :=
    (Nat.le_ceil _).trans (by exact_mod_cast (le_max_right _ _).trans hN)
  have hNr : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  have hδN : 1 ≤ δ / 4 * N := by
    rw [div_le_iff₀ hδ] at hN3
    linarith
  suffices ∃ m : ℕ, δ / 4 * N ≤ m ∧ f N k ≤ B * r ^ m by
    obtain ⟨m, hm, hf⟩ := this
    exact hf.trans (mul_le_mul_of_nonneg_left (pow_tail hr0 hr1 hm) hBpos.le)
  rcases le_abs.1 hk with hk | hk
  · -- upper tail: `(x* + δ) N ≤ k`
    have hy1 : 0 < xstar + δ / 2 := by linarith
    have hc : (1 + δ / 2) * (xstar + δ / 2) ≤ xstar + δ := by
      have := mul_nonneg hδ.le (show (0 : ℝ) ≤ 1 - xstar - δ / 2 by linarith)
      linarith
    have hkN' : (1 + δ / 2) * N ≤ (k : ℝ) / (xstar + δ / 2) := by
      rw [le_div_iff₀ hy1]
      have := mul_le_mul_of_nonneg_right hc hNr
      linarith
    obtain ⟨m, hm, hf⟩ := tail_upper f B hpos hB hρ hq1 hy1 h1 hN2 hkN hkN' hδ.le hδN
    exact ⟨m, hm, hf.trans (mul_le_mul_of_nonneg_left
      (pow_le_pow_left₀ (inv_nonneg.2 hq1pos.le) hrq1 _) hBpos.le)⟩
  · -- lower tail: `k ≤ (x* − δ) N`
    have hy0 : 0 < xstar - δ / 2 := by linarith
    have hy0' : xstar - δ / 2 < 1 := by linarith
    have hc : xstar - δ ≤ (1 - δ / 2) * (xstar - δ / 2) := by
      have := mul_nonneg hδ.le (show (0 : ℝ) ≤ 1 - xstar + δ / 2 by linarith)
      linarith
    have hkN' : (k : ℝ) / (xstar - δ / 2) ≤ (1 - δ / 2) * N := by
      rw [div_le_iff₀ hy0]
      have := mul_le_mul_of_nonneg_right hc hNr
      linarith
    obtain ⟨m, hm, hf⟩ := tail_lower f B hBpos.le hB hρ hρpos hy0 hy0' h0 hN1 hkN'
      hδ.le (by linarith) hδN
    exact ⟨m, hm, hf.trans (mul_le_mul_of_nonneg_left
      (pow_le_pow_left₀ hq0.le hrq0 _) hBpos.le)⟩

end PiSeries.SunPi2
