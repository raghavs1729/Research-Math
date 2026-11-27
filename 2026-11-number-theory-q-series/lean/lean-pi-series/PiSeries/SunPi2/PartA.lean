/-
# Sun's `64√2/π²` series, Part A (exact): the WZ pair on the triangle `k < n`

Work package WP-A of `N2-PLAN.md`.  With the definitions of `PiSeries/SunPi2/Defs.lean`:

* `a_pos`, `a_succ_left` : positivity and the column relation `a (n+1) k = a n k * ρ n k` (`k < n`);
* `T_one`, `T_succ`, `sum_t_eq` : `∑_{n<N} t n = 256 T N` (WZ telescoping on the triangle);
* `t_pos`, `summable_t`, `a_le_T`, `T_le` : `a n k ≤ T n ≤ (∑' t)/256`;
* `a_eq_cf` : the factorial closed form of `a N k` for `k < N ≤ 2k+2`;
* `hasSum_t_of_tendsto_T` : `T N → L` implies `HasSum t (256 L)`.

Auxiliary declarations live in the namespace `PiSeries.SunPi2.PartA`.
-/
import PiSeries.SunPi2.Defs

open Finset Filter Topology

namespace PiSeries.SunPi2

namespace PartA

/-! ### Signs of the odd factors -/

/-- an odd integer shifted by multiples of `8k`, `4n` is non-zero -/
lemma odd_ne (n k : ℕ) (c : ℤ) (hc : c % 2 = 1) : (8 * (k : ℝ) - 4 * n + c) ≠ 0 := by
  have h : (8 * (k : ℤ) - 4 * n + c) ≠ 0 := by omega
  exact_mod_cast h

/-- two consecutive odd integers `v`, `v + 2` with `v ≡ 1 (mod 4)` have a positive product -/
lemma odd2_pos (v : ℤ) (hv : v % 4 = 1) : (0 : ℝ) < (v : ℝ) * ((v : ℝ) + 2) := by
  rcases (by omega : 1 ≤ v ∨ v ≤ -3) with h | h
  · have h' : (1 : ℝ) ≤ v := by exact_mod_cast h
    exact mul_pos (by linarith) (by linarith)
  · have h' : (v : ℝ) ≤ -3 := by exact_mod_cast h
    exact mul_pos_of_neg_of_neg (by linarith) (by linarith)

/-- four consecutive odd integers centred at a multiple of 4 have a positive product -/
lemma odd4_pos (u : ℤ) (hu : u % 4 = 0) :
    (0 : ℝ) < ((u : ℝ) - 3) * ((u : ℝ) - 1) * ((u : ℝ) + 1) * ((u : ℝ) + 3) := by
  have e : ((u : ℝ) - 3) * ((u : ℝ) - 1) * ((u : ℝ) + 1) * ((u : ℝ) + 3)
      = ((u : ℝ) ^ 2 - 9) * ((u : ℝ) ^ 2 - 1) := by ring
  rw [e]
  rcases (by omega : u = 0 ∨ 4 ≤ u ∨ u ≤ -4) with h | h | h
  · subst h; norm_num
  · have h' : (4 : ℝ) ≤ u := by exact_mod_cast h
    have h2 : (16 : ℝ) ≤ (u : ℝ) ^ 2 := by nlinarith
    exact mul_pos (by linarith) (by linarith)
  · have h' : (u : ℝ) ≤ -4 := by exact_mod_cast h
    have h2 : (16 : ℝ) ≤ (u : ℝ) ^ 2 := by nlinarith
    exact mul_pos (by linarith) (by linarith)

lemma sigma_odd_pos (n k : ℕ) :
    (0 : ℝ) < (8 * (k : ℝ) - 4 * n + 9) * (8 * (k : ℝ) - 4 * n + 11) * (8 * (k : ℝ) - 4 * n + 13)
      * (8 * (k : ℝ) - 4 * n + 15) := by
  have h := odd4_pos (8 * (k : ℤ) - 4 * n + 12) (by omega)
  push_cast at h
  have e : (8 * (k : ℝ) - 4 * n + 9) * (8 * (k : ℝ) - 4 * n + 11) * (8 * (k : ℝ) - 4 * n + 13)
      * (8 * (k : ℝ) - 4 * n + 15)
      = (8 * (k : ℝ) - 4 * n + 12 - 3) * (8 * (k : ℝ) - 4 * n + 12 - 1) * (8 * (k : ℝ) - 4 * n + 12 + 1)
        * (8 * (k : ℝ) - 4 * n + 12 + 3) := by ring
  rw [e]; exact h

lemma rho_odd_pos (n k : ℕ) :
    (0 : ℝ) < (8 * (k : ℝ) - 4 * n + 5) * (8 * (k : ℝ) - 4 * n + 7) := by
  have h := odd2_pos (8 * (k : ℤ) - 4 * n + 5) (by omega)
  push_cast at h
  have e : (8 * (k : ℝ) - 4 * n + 5) * (8 * (k : ℝ) - 4 * n + 7)
      = (8 * (k : ℝ) - 4 * n + 5) * (8 * (k : ℝ) - 4 * n + 5 + 2) := by ring
  rw [e]; exact h

lemma σ_pos {n k : ℕ} (h : k + 1 < n) : 0 < σ n k := by
  have h1 : (0 : ℝ) < (n : ℝ) - k - 1 := by
    have : ((k : ℝ) + 1) < n := by exact_mod_cast h
    linarith
  have h2 := sigma_odd_pos n k
  unfold σ
  apply div_pos
  · have : (0 : ℝ) < ((n : ℝ) - k - 1) ^ 2 := by positivity
    positivity
  · exact mul_pos (by positivity) h2

lemma ρ_pos {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : 0 < ρ n k := by
  have h0 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have h1 : (0 : ℝ) < (n : ℝ) - k := by
    have : (k : ℝ) < n := by exact_mod_cast hk
    linarith
  have h2 := rho_odd_pos n k
  have h3 : (0 : ℝ) < 2 * (n : ℝ) - 1 := by linarith
  have h4 : (0 : ℝ) < 4 * (n : ℝ) - 1 := by linarith
  unfold ρ
  apply div_pos
  · have e : (4 * (n : ℝ) + 1) * (4 * (n : ℝ) + 3) ^ 2 * (8 * (k : ℝ) + 4 * n + 1)
        * (8 * (k : ℝ) + 4 * n + 3) * (8 * (k : ℝ) - 4 * n + 5) * (8 * (k : ℝ) - 4 * n + 7)
        = (4 * (n : ℝ) + 1) * (4 * (n : ℝ) + 3) ^ 2 * (8 * (k : ℝ) + 4 * n + 1)
          * (8 * (k : ℝ) + 4 * n + 3) * ((8 * (k : ℝ) - 4 * n + 5) * (8 * (k : ℝ) - 4 * n + 7)) := by
      ring
    rw [e]
    exact mul_pos (by positivity) h2
  · have : (0 : ℝ) < (n : ℝ) := by linarith
    positivity

lemma A_pos (n : ℕ) : 0 < A n := by
  unfold A
  refine mul_pos (by norm_num) (Finset.prod_pos fun m hm => ?_)
  have hm' := Finset.mem_Ico.mp hm
  exact ρ_pos hm'.1 (by omega)

/-- compatibility of the two ratios (all denominators are non-zero for `k + 1 < n`) -/
lemma compat {n k : ℕ} (h : k + 1 < n) : ρ n k * σ (n + 1) k = σ n k * ρ n (k + 1) := by
  have hn : ((k : ℝ) + 2) ≤ n := by exact_mod_cast h
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have d1 : (n : ℝ) ≠ 0 := by linarith
  have d2 : (2 * (n : ℝ) - 1) ≠ 0 := by linarith
  have d3 : (4 * (n : ℝ) - 1) ≠ 0 := by linarith
  have d4 : ((n : ℝ) - k) ≠ 0 := by linarith
  have d5 : ((n : ℝ) + k + 1) ≠ 0 := by linarith
  have d6 : (2 * (k : ℝ) + 1) ≠ 0 := by linarith
  have d7 : ((n : ℝ) + 1 + k + 1) ≠ 0 := by linarith
  have d8 : ((n : ℝ) - (k + 1)) ≠ 0 := by linarith
  have d9 : ((n : ℝ) + (k + 1) + 1) ≠ 0 := by linarith
  have d10 : (2 * ((k : ℝ) + 1) + 1) ≠ 0 := by linarith
  have o5 : (8 * (k : ℝ) - 4 * n + 5) ≠ 0 := by exact_mod_cast odd_ne n k 5 (by norm_num)
  have o7 : (8 * (k : ℝ) - 4 * n + 7) ≠ 0 := by exact_mod_cast odd_ne n k 7 (by norm_num)
  have o9 : (8 * (k : ℝ) - 4 * n + 9) ≠ 0 := by exact_mod_cast odd_ne n k 9 (by norm_num)
  have o11 : (8 * (k : ℝ) - 4 * n + 11) ≠ 0 := by exact_mod_cast odd_ne n k 11 (by norm_num)
  have o13 : (8 * (k : ℝ) - 4 * n + 13) ≠ 0 := by exact_mod_cast odd_ne n k 13 (by norm_num)
  have o15 : (8 * (k : ℝ) - 4 * n + 15) ≠ 0 := by exact_mod_cast odd_ne n k 15 (by norm_num)
  have p5 : (8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 9) ≠ 0 := by
    rw [show 8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 9 = 8 * (k : ℝ) - 4 * n + 5 by ring]; exact o5
  have p7 : (8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 11) ≠ 0 := by
    rw [show 8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 11 = 8 * (k : ℝ) - 4 * n + 7 by ring]; exact o7
  have p9 : (8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 13) ≠ 0 := by
    rw [show 8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 13 = 8 * (k : ℝ) - 4 * n + 9 by ring]; exact o9
  have p11 : (8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 15) ≠ 0 := by
    rw [show 8 * (k : ℝ) - 4 * ((n : ℝ) + 1) + 15 = 8 * (k : ℝ) - 4 * n + 11 by ring]; exact o11
  unfold ρ σ
  push_cast
  field_simp
  ring

end PartA

open PartA

/-! ### Positivity and the column relation -/

theorem a_pos {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : 0 < a n k := by
  have _ := hn
  unfold a
  refine mul_pos (A_pos n) (Finset.prod_pos fun j hj => ?_)
  have hj' := Finset.mem_range.mp hj
  exact σ_pos (by omega)

theorem a_succ_left {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : a (n + 1) k = a n k * ρ n k := by
  induction k with
  | zero => rw [a_zero_right, a_zero_right, A_succ hn]
  | succ k ih =>
    rw [a_succ_right, a_succ_right, ih (by omega)]
    have h := compat (n := n) (k := k) hk
    linear_combination a n k * h

theorem T_one : T 1 = t 0 / 256 := by
  simp [T, t, a_zero_right, A_one]

theorem t_pos (n : ℕ) : 0 < t n := by
  unfold t
  have h1 : (0 : ℝ) < (Nat.choose (2 * n) n : ℝ) := by
    exact_mod_cast Nat.choose_pos (by omega)
  have h2 : (0 : ℝ) < (Nat.choose (4 * n) (2 * n) : ℝ) := by
    exact_mod_cast Nat.choose_pos (by omega)
  positivity

theorem a_le_T {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : a n k ≤ T n := by
  unfold T
  exact Finset.single_le_sum (f := fun j => a n j)
    (fun j hj => (a_pos hn (Finset.mem_range.mp hj)).le) (Finset.mem_range.mpr hk)

namespace PartA

/-! ### The WZ certificate (copied from `ramanujan-series/lean/RamanujanSeries/Sun1pi2.lean`,
with numerators and denominators of `σ`, `ρ` normalised as in `Defs.lean`) -/

/-- numerator of the certificate -/
def NUM (n k : ℝ) : ℝ :=
  (3145728*n + 1310720)*k^6 + (-6291456*n^2 + 3670016*n + 2621440)*k^5
  + (-786432*n^3 - 11337728*n^2 - 344064*n + 1724416)*k^4
  + (9437184*n^4 + 786432*n^3 - 8716288*n^2 - 1884160*n + 376832)*k^3
  + (-8060928*n^5 + 10010624*n^4 + 3973120*n^3 - 3434496*n^2 - 892672*n - 15936)*k^2
  + (2752512*n^6 - 6324224*n^5 + 1916928*n^4 + 2363392*n^3 - 525824*n^2 - 170112*n - 12672)*k
  + (-344064*n^7 + 1134592*n^6 - 922624*n^5 - 195328*n^4 + 318016*n^3 - 1680*n^2 - 10044*n - 945)

/-- denominator of the certificate -/
def RD (n k : ℝ) : ℝ := 131072 * n * (2*n - 1) * (2*n + 1) * (4*n - 1) * (k - n)^2 * (2*k + 1)

/-- the WZ certificate `R(n,k)` -/
noncomputable def Rc (n k : ℝ) : ℝ := NUM n k / RD n k

/-- numerator and denominator of `σ`, as functions of real variables -/
def σN (n k : ℝ) : ℝ :=
  (2 * k + 3) * (n - k - 1) ^ 2
      * ((8 * k + 4 * n + 1) * (8 * k + 4 * n + 3) * (8 * k + 4 * n + 5) * (8 * k + 4 * n + 7))
def σD (n k : ℝ) : ℝ :=
  (2 * k + 1) * (n + k + 1) ^ 2
      * ((8 * k - 4 * n + 9) * (8 * k - 4 * n + 11) * (8 * k - 4 * n + 13) * (8 * k - 4 * n + 15))

/-- numerator and denominator of `ρ`, as functions of real variables -/
def ρN (n k : ℝ) : ℝ :=
  (4 * n + 1) * (4 * n + 3) ^ 2 * (8 * k + 4 * n + 1) * (8 * k + 4 * n + 3)
      * (8 * k - 4 * n + 5) * (8 * k - 4 * n + 7)
def ρD (n k : ℝ) : ℝ := 131072 * n * (2 * n - 1) * (4 * n - 1) * (n - k) ^ 2 * (n + k + 1) ^ 2

lemma σ_eq (n k : ℕ) : σ n k = σN n k / σD n k := rfl

lemma ρ_eq (n k : ℕ) : ρ n k = ρN n k / ρD n k := rfl

set_option maxRecDepth 100000 in
set_option maxHeartbeats 4000000 in
/-- the WZ equation, cleared of denominators (a polynomial identity) -/
theorem wz_poly (n k : ℝ) :
    (ρN n k - ρD n k) * (RD n (k + 1) * σD n k * RD n k)
      = ρD n k * (NUM n (k + 1) * σN n k * RD n k - NUM n k * (RD n (k + 1) * σD n k)) := by
  unfold ρN ρD RD σN σD NUM
  ring

set_option maxRecDepth 100000 in
set_option maxHeartbeats 4000000 in
/-- the boundary row `k = n - 1` (here `n = m + 1`, `k = m`), cleared of denominators -/
theorem wz_boundary_poly (m : ℝ) :
    (ρN (m + 1) m - ρD (m + 1) m) * (σD (m + 1 + 1) m * RD (m + 1) m)
      = -(ρN (m + 1) m * σN (m + 1 + 1) m * RD (m + 1) m)
        - NUM (m + 1) m * (ρD (m + 1) m * σD (m + 1 + 1) m) := by
  unfold ρN ρD RD σN σD NUM
  ring

/-- the WZ equation `ρ - 1 = R(n,k+1) σ - R(n,k)` for `k + 1 < n` -/
lemma wz_nat {n k : ℕ} (h : k + 1 < n) :
    ρ n k - 1 = Rc n ((k : ℝ) + 1) * σ n k - Rc n k := by
  have hn : ((k : ℝ) + 2) ≤ n := by exact_mod_cast h
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have p1 : (0 : ℝ) < n := by linarith
  have p2 : (0 : ℝ) < 2 * (n : ℝ) - 1 := by linarith
  have p3 : (0 : ℝ) < 4 * (n : ℝ) - 1 := by linarith
  have p4 : (0 : ℝ) < (n : ℝ) - k := by linarith
  have p5 : (0 : ℝ) < ((k : ℝ) - n) ^ 2 := by nlinarith
  have p6 : (0 : ℝ) < ((k : ℝ) + 1 - n) ^ 2 := by nlinarith
  have hρ : ρD n k ≠ 0 := by
    have : 0 < ρD n k := by unfold ρD; positivity
    exact this.ne'
  have hσ : σD n k ≠ 0 := by
    have : 0 < σD n k := by unfold σD; exact mul_pos (by positivity) (sigma_odd_pos n k)
    exact this.ne'
  have hR : RD n k ≠ 0 := by
    have : 0 < RD n k := by unfold RD; positivity
    exact this.ne'
  have hR' : RD n ((k : ℝ) + 1) ≠ 0 := by
    have : 0 < RD n ((k : ℝ) + 1) := by unfold RD; positivity
    exact this.ne'
  rw [ρ_eq, σ_eq]
  unfold Rc
  rw [div_sub_one hρ, div_mul_div_comm, div_sub_div _ _ (mul_ne_zero hR' hσ) hR,
    div_eq_div_iff hρ (mul_ne_zero (mul_ne_zero hR' hσ) hR)]
  linear_combination wz_poly n k

/-- the boundary row: `ρ(n,n-1) - 1 = -ρ(n,n-1) σ(n+1,n-1) - R(n,n-1)` with `n = m + 1` -/
lemma wz_bdry (m : ℕ) :
    ρ (m + 1) m - 1 = -(ρ (m + 1) m * σ (m + 1 + 1) m) - Rc ((m : ℝ) + 1) m := by
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have p2 : (0 : ℝ) < 2 * ((m : ℝ) + 1) - 1 := by linarith
  have p2' : (0 : ℝ) < 2 * ((m : ℝ) + 1) + 1 := by linarith
  have p3 : (0 : ℝ) < 4 * ((m : ℝ) + 1) - 1 := by linarith
  have p4 : (0 : ℝ) < (m : ℝ) + 1 - m := by linarith
  have p5 : (0 : ℝ) < ((m : ℝ) - ((m : ℝ) + 1)) ^ 2 := by
    rw [show ((m : ℝ) - ((m : ℝ) + 1)) ^ 2 = 1 by ring]; norm_num
  have q1 : (0 : ℝ) < 8 * (m : ℝ) - 4 * ((m : ℝ) + 1 + 1) + 9 := by linarith
  have q2 : (0 : ℝ) < 8 * (m : ℝ) - 4 * ((m : ℝ) + 1 + 1) + 11 := by linarith
  have q3 : (0 : ℝ) < 8 * (m : ℝ) - 4 * ((m : ℝ) + 1 + 1) + 13 := by linarith
  have q4 : (0 : ℝ) < 8 * (m : ℝ) - 4 * ((m : ℝ) + 1 + 1) + 15 := by linarith
  have hρ : ρD ((m : ℝ) + 1) m ≠ 0 := by
    have : 0 < ρD ((m : ℝ) + 1) m := by unfold ρD; positivity
    exact this.ne'
  have hσ : σD ((m : ℝ) + 1 + 1) m ≠ 0 := by
    have : 0 < σD ((m : ℝ) + 1 + 1) m := by unfold σD; positivity
    exact this.ne'
  have hR : RD ((m : ℝ) + 1) m ≠ 0 := by
    have : 0 < RD ((m : ℝ) + 1) m := by unfold RD; positivity
    exact this.ne'
  rw [ρ_eq, σ_eq]
  unfold Rc
  simp only [Nat.cast_add, Nat.cast_one]
  rw [div_sub_one hρ, div_mul_div_comm, ← neg_div, div_sub_div _ _ (mul_ne_zero hρ hσ) hR,
    div_eq_div_iff hρ (mul_ne_zero (mul_ne_zero hρ hσ) hR)]
  linear_combination (ρD ((m : ℝ) + 1) m) * wz_boundary_poly m

/-- telescoping of the WZ equation along the row `N`, `k < m < N` -/
lemma tele (N : ℕ) : ∀ m : ℕ, m < N →
    ∑ k ∈ range m, (a (N + 1) k - a N k) = Rc N m * a N m - Rc N 0 * a N 0 := by
  intro m
  induction m with
  | zero => intro _; simp
  | succ m ih =>
    intro hm
    rw [sum_range_succ, ih (by omega), a_succ_left (by omega) (by omega), a_succ_right]
    have h := wz_nat (n := N) (k := m) hm
    push_cast
    linear_combination a N m * h

/-! ### The closed form of `A n = a(n,0)` -/

/-- `(n+1) C(2n+2,n+1) = 2(2n+1) C(2n,n)`, over `ℝ` -/
lemma choose_two_mul_succ (n : ℕ) :
    ((n : ℝ) + 1) * (Nat.choose (2 * (n + 1)) (n + 1) : ℝ)
      = 2 * (2 * (n : ℝ) + 1) * (Nat.choose (2 * n) n : ℝ) := by
  have h := Nat.succ_mul_centralBinom_succ n
  rw [Nat.centralBinom_eq_two_mul_choose, Nat.centralBinom_eq_two_mul_choose] at h
  exact_mod_cast h

/-- `(2n+1)(2n+2) C(4n+4,2n+2) = 4(4n+1)(4n+3) C(4n,2n)` -/
lemma choose_four_step (n : ℕ) :
    (2 * (n : ℝ) + 1) * (2 * n + 2) * (Nat.choose (4 * (n + 1)) (2 * (n + 1)) : ℝ)
      = 4 * (4 * (n : ℝ) + 1) * (4 * n + 3) * (Nat.choose (4 * n) (2 * n) : ℝ) := by
  have e1 := choose_two_mul_succ (2 * n)
  have e2 := choose_two_mul_succ (2 * n + 1)
  rw [show 4 * (n + 1) = 2 * (2 * n + 1 + 1) by ring, show 2 * (n + 1) = 2 * n + 1 + 1 by ring,
    show 4 * n = 2 * (2 * n) by ring]
  push_cast at e1 e2
  linear_combination (2 * (n : ℝ) + 1) * e2 + 2 * (2 * (2 * (n : ℝ) + 1) + 1) * e1

lemma four_sub_ne (m : ℕ) (c : ℤ) (hc : c % 4 ≠ 0) : (4 * (m : ℝ) - c) ≠ 0 := by
  have h : (4 * (m : ℤ) - c) ≠ 0 := by omega
  exact_mod_cast h

/-- `ρ(n,0)` for `n = m + 1` -/
lemma ρ_zero (m : ℕ) :
    ρ (m + 1) 0 = (4 * (m : ℝ) + 5) ^ 2 * (4 * (m : ℝ) + 7) ^ 3 * ((4 * (m : ℝ) - 1) * (4 * (m : ℝ) - 3))
      / (131072 * ((m : ℝ) + 1) ^ 3 * (2 * (m : ℝ) + 1) * (4 * (m : ℝ) + 3) * ((m : ℝ) + 2) ^ 2) := by
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have p2 : (0 : ℝ) < 2 * ((m : ℝ) + 1) - 1 := by linarith
  have p3 : (0 : ℝ) < 4 * ((m : ℝ) + 1) - 1 := by linarith
  have p4 : (0 : ℝ) < (m : ℝ) + 1 - 0 := by linarith
  have hρ : ρD ((m : ℝ) + 1) 0 ≠ 0 := by
    have : 0 < ρD ((m : ℝ) + 1) 0 := by unfold ρD; positivity
    exact this.ne'
  rw [ρ_eq]
  simp only [Nat.cast_add, Nat.cast_one, Nat.cast_zero]
  rw [div_eq_div_iff hρ (by positivity)]
  unfold ρN ρD
  ring

/-- **closed form of `A n`**, `n = m + 1`:
`A n = 512 n³ (2n-1) C(2n,n)² C(4n,2n)³ / (2^{22n} (4n-7)(4n-5)(4n-3))` -/
lemma A_closed (m : ℕ) :
    A (m + 1) = 512 * ((m : ℝ) + 1) ^ 3 * (2 * (m : ℝ) + 1) * (Nat.choose (2 * (m + 1)) (m + 1) : ℝ) ^ 2
        * (Nat.choose (4 * (m + 1)) (2 * (m + 1)) : ℝ) ^ 3
      / (2 ^ (22 * (m + 1)) * ((4 * (m : ℝ) - 3) * (4 * (m : ℝ) - 1) * (4 * (m : ℝ) + 1))) := by
  induction m with
  | zero =>
    rw [A_one]
    norm_num [Nat.choose]
  | succ m ih =>
    have o3 : (4 * (m : ℝ) - 3) ≠ 0 := by exact_mod_cast four_sub_ne m 3 (by norm_num)
    have o1 : (4 * (m : ℝ) - 1) ≠ 0 := by exact_mod_cast four_sub_ne m 1 (by norm_num)
    have e2 := choose_two_mul_succ (m + 1)
    have e4 := choose_four_step (m + 1)
    push_cast at e2 e4
    have c2 : (Nat.choose (2 * (m + 1 + 1)) (m + 1 + 1) : ℝ)
        = 2 * (2 * ((m : ℝ) + 1) + 1) * (Nat.choose (2 * (m + 1)) (m + 1) : ℝ) / ((m : ℝ) + 1 + 1) := by
      rw [eq_div_iff (by positivity)]; linear_combination e2
    have c4 : (Nat.choose (4 * (m + 1 + 1)) (2 * (m + 1 + 1)) : ℝ)
        = 4 * (4 * ((m : ℝ) + 1) + 1) * (4 * ((m : ℝ) + 1) + 3) * (Nat.choose (4 * (m + 1)) (2 * (m + 1)) : ℝ)
          / ((2 * ((m : ℝ) + 1) + 1) * (2 * ((m : ℝ) + 1) + 2)) := by
      rw [eq_div_iff (by positivity)]; linear_combination e4
    rw [A_succ (by omega), ih, ρ_zero, c2, c4, show 22 * (m + 1 + 1) = 22 * (m + 1) + 22 by ring, pow_add]
    push_cast
    have o5 : (4 * ((m : ℝ) + 1) - 3) ≠ 0 := by
      have : (0 : ℝ) ≤ m := Nat.cast_nonneg m
      have : (0 : ℝ) < 4 * ((m : ℝ) + 1) - 3 := by linarith
      exact this.ne'
    have o6 : (4 * ((m : ℝ) + 1) - 1) ≠ 0 := by
      have : (0 : ℝ) ≤ m := Nat.cast_nonneg m
      have : (0 : ℝ) < 4 * ((m : ℝ) + 1) - 1 := by linarith
      exact this.ne'
    -- `field_simp` looks for the side conditions in its own normal form
    have o3' : ((m : ℝ) * 4 - 3) ≠ 0 := by rwa [mul_comm]
    have o1' : ((m : ℝ) * 4 - 1) ≠ 0 := by rwa [mul_comm]
    have o5' : (((m : ℝ) + 1) * 4 - 3) ≠ 0 := by rwa [mul_comm]
    have o6' : (((m : ℝ) + 1) * 4 - 1) ≠ 0 := by rwa [mul_comm]
    field_simp
    ring

/-- `g(n,0) = R(n,0) a(n,0) = -t(n)/256`, `n = m + 1` -/
lemma g_zero (m : ℕ) : -(Rc ((m : ℝ) + 1) 0 * A (m + 1)) = t (m + 1) / 256 := by
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have o3 : (4 * (m : ℝ) - 3) ≠ 0 := by exact_mod_cast four_sub_ne m 3 (by norm_num)
  have o1 : (4 * (m : ℝ) - 1) ≠ 0 := by exact_mod_cast four_sub_ne m 1 (by norm_num)
  have eN : NUM ((m : ℝ) + 1) 0 = -((4 * (m : ℝ) - 3) * (4 * (m : ℝ) - 1) * (4 * (m : ℝ) + 1)
      * (4 * (m : ℝ) + 3) * (1344 * ((m : ℝ) + 1) ^ 3 + 944 * ((m : ℝ) + 1) ^ 2 + 156 * ((m : ℝ) + 1) + 9)) := by
    unfold NUM; ring
  have eD : RD ((m : ℝ) + 1) 0 = 131072 * ((m : ℝ) + 1) ^ 3 * (2 * (m : ℝ) + 1) * (2 * (m : ℝ) + 3)
      * (4 * (m : ℝ) + 3) := by
    unfold RD; ring
  rw [A_closed]
  unfold Rc t
  rw [eN, eD]
  push_cast
  field_simp
  ring

end PartA

/-! ### The telescoped identity `∑_{n<N} t n = 256 T N` -/

theorem T_succ {N : ℕ} (hN : 1 ≤ N) : T (N + 1) = T N + t N / 256 := by
  obtain ⟨m, rfl⟩ : ∃ m, N = m + 1 := ⟨N - 1, by omega⟩
  have h1 := tele (m + 1) m (by omega)
  have h2 := wz_bdry m
  have h3 := g_zero m
  rw [sum_sub_distrib, a_zero_right] at h1
  push_cast at h1
  have e1 : T (m + 1 + 1)
      = ∑ k ∈ range m, a (m + 1 + 1) k + a (m + 1) m * ρ (m + 1) m
        + a (m + 1) m * ρ (m + 1) m * σ (m + 1 + 1) m := by
    unfold T
    rw [sum_range_succ, sum_range_succ, a_succ_right, a_succ_left (by omega) (by omega)]
  have e2 : T (m + 1) = ∑ k ∈ range m, a (m + 1) k + a (m + 1) m := by
    unfold T
    rw [sum_range_succ]
  rw [e1, e2]
  linear_combination h1 + (a (m + 1) m) * h2 + h3

theorem sum_t_eq {N : ℕ} (hN : 1 ≤ N) : ∑ n ∈ range N, t n = 256 * T N := by
  induction N with
  | zero => omega
  | succ N ih =>
    rcases Nat.eq_zero_or_pos N with h0 | h0
    · subst h0
      rw [zero_add, sum_range_one, T_one]; ring
    · rw [sum_range_succ, ih h0, T_succ h0]; ring

/-! ### Summability and the bound `a n k ≤ T n ≤ (∑' t)/256` -/

theorem summable_t : Summable t := by
  have hg (k : ℕ) : Summable (fun n : ℕ => (n : ℝ) ^ k * (1 / 64 : ℝ) ^ n) :=
    summable_pow_mul_geometric_of_norm_lt_one k
      (by rw [Real.norm_eq_abs, abs_of_pos (by norm_num)]; norm_num)
  have hs : Summable (fun n : ℕ => 1344 * ((n : ℝ) ^ 3 * (1 / 64 : ℝ) ^ n)
      + 944 * ((n : ℝ) ^ 2 * (1 / 64 : ℝ) ^ n) + 156 * ((n : ℝ) ^ 1 * (1 / 64 : ℝ) ^ n)
      + 9 * ((n : ℝ) ^ 0 * (1 / 64 : ℝ) ^ n)) :=
    ((((hg 3).mul_left 1344).add ((hg 2).mul_left 944)).add ((hg 1).mul_left 156)).add
      ((hg 0).mul_left 9)
  refine Summable.of_nonneg_of_le (fun n => (t_pos n).le) (fun n => ?_) hs
  have c2 : (Nat.choose (2 * n) n : ℝ) ≤ 2 ^ (2 * n) := by
    exact_mod_cast Nat.choose_le_two_pow (2 * n) n
  have c4 : (Nat.choose (4 * n) (2 * n) : ℝ) ≤ 2 ^ (4 * n) := by
    exact_mod_cast Nat.choose_le_two_pow (4 * n) (2 * n)
  have h1 : (Nat.choose (2 * n) n : ℝ) ^ 2 * (Nat.choose (4 * n) (2 * n) : ℝ) ^ 3
      ≤ ((2 : ℝ) ^ (2 * n)) ^ 2 * ((2 : ℝ) ^ (4 * n)) ^ 3 :=
    mul_le_mul (pow_le_pow_left₀ (by positivity) c2 2) (pow_le_pow_left₀ (by positivity) c4 3)
      (by positivity) (by positivity)
  have h2 : ((2 : ℝ) ^ (2 * n)) ^ 2 * ((2 : ℝ) ^ (4 * n)) ^ 3 = 2 ^ (16 * n) := by
    rw [← pow_mul, ← pow_mul, ← pow_add]; congr 1; ring
  have h3 : (2 : ℝ) ^ (22 * n) = 2 ^ (16 * n) * 64 ^ n := by
    rw [show (64 : ℝ) = 2 ^ 6 by norm_num, ← pow_mul, ← pow_add]; congr 1; ring
  rw [h2] at h1
  have hP : (0 : ℝ) < 1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9 := by positivity
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have e : 1344 * ((n : ℝ) ^ 3 * (1 / 64 : ℝ) ^ n) + 944 * ((n : ℝ) ^ 2 * (1 / 64 : ℝ) ^ n)
      + 156 * ((n : ℝ) ^ 1 * (1 / 64 : ℝ) ^ n) + 9 * ((n : ℝ) ^ 0 * (1 / 64 : ℝ) ^ n)
      = (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * 2 ^ (16 * n)
        / (2 ^ (16 * n) * 64 ^ n) := by
    rw [one_div_pow]
    field_simp
  rw [e]
  unfold t
  rw [h3, div_le_div_iff₀ (by positivity) (by positivity)]
  have h4 : (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9)
      * ((Nat.choose (2 * n) n : ℝ) ^ 2 * (Nat.choose (4 * n) (2 * n) : ℝ) ^ 3)
      ≤ (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * 2 ^ (16 * n) :=
    mul_le_mul_of_nonneg_left h1 hP.le
  have h5 : (0 : ℝ) ≤ (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * 2 ^ (16 * n)
      * (2 ^ (16 * n) * 64 ^ n) := by positivity
  have h6 : (0 : ℝ) < 2 ^ (16 * n) * 64 ^ n := by positivity
  calc (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * (Nat.choose (2 * n) n : ℝ) ^ 2
        * (Nat.choose (4 * n) (2 * n) : ℝ) ^ 3 * (2 ^ (16 * n) * 64 ^ n)
      = (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9)
          * ((Nat.choose (2 * n) n : ℝ) ^ 2 * (Nat.choose (4 * n) (2 * n) : ℝ) ^ 3)
          * (2 ^ (16 * n) * 64 ^ n) := by ring
    _ ≤ (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * 2 ^ (16 * n)
          * (2 ^ (16 * n) * 64 ^ n) := mul_le_mul_of_nonneg_right h4 h6.le
    _ ≤ (1344 * (n : ℝ) ^ 3 + 944 * (n : ℝ) ^ 2 + 156 * n + 9) * 2 ^ (16 * n)
          * ((2 * (n : ℝ) + 1) * (2 ^ (16 * n) * 64 ^ n)) := by nlinarith

theorem T_le (N : ℕ) (hN : 1 ≤ N) : T N ≤ (∑' n, t n) / 256 := by
  have h := summable_t.sum_le_tsum (range N) (fun n _ => (t_pos n).le)
  rw [sum_t_eq hN] at h
  linarith

theorem hasSum_t_of_tendsto_T {L : ℝ} (h : Tendsto T atTop (𝓝 L)) : HasSum t (256 * L) := by
  rw [hasSum_iff_tendsto_nat_of_nonneg (fun n => (t_pos n).le)]
  have h' : Tendsto (fun N => 256 * T N) atTop (𝓝 (256 * L)) := h.const_mul 256
  refine h'.congr' ?_
  filter_upwards [eventually_ge_atTop 1] with N hN
  exact (sum_t_eq hN).symm

/-! ### The factorial closed form `a N k = cf N k` for `k < N ≤ 2k + 2` -/

namespace PartA

lemma fact1 (m : ℕ) : (((m + 1).factorial : ℕ) : ℝ) = ((m : ℝ) + 1) * (m.factorial : ℝ) := by
  push_cast [Nat.factorial_succ]; ring

lemma fact2 (m : ℕ) :
    (((m + 2).factorial : ℕ) : ℝ) = ((m : ℝ) + 2) * ((m : ℝ) + 1) * (m.factorial : ℝ) := by
  rw [show m + 2 = (m + 1) + 1 from rfl, fact1, fact1]; push_cast; ring

lemma fact4 (m : ℕ) :
    (((m + 4).factorial : ℕ) : ℝ)
      = ((m : ℝ) + 4) * ((m : ℝ) + 3) * ((m : ℝ) + 2) * ((m : ℝ) + 1) * (m.factorial : ℝ) := by
  rw [show m + 4 = (m + 2) + 2 from rfl, fact2, fact2]; push_cast; ring

lemma fact6 (m : ℕ) :
    (((m + 6).factorial : ℕ) : ℝ)
      = ((m : ℝ) + 6) * ((m : ℝ) + 5) * ((m : ℝ) + 4) * ((m : ℝ) + 3) * ((m : ℝ) + 2) * ((m : ℝ) + 1)
        * (m.factorial : ℝ) := by
  rw [show m + 6 = (m + 4) + 2 from rfl, fact2, fact4]; push_cast; ring

lemma fact12 (m : ℕ) :
    (((m + 12).factorial : ℕ) : ℝ)
      = ((m : ℝ) + 12) * ((m : ℝ) + 11) * ((m : ℝ) + 10) * ((m : ℝ) + 9) * ((m : ℝ) + 8) * ((m : ℝ) + 7)
        * (((m : ℝ) + 6) * ((m : ℝ) + 5) * ((m : ℝ) + 4) * ((m : ℝ) + 3) * ((m : ℝ) + 2) * ((m : ℝ) + 1))
        * (m.factorial : ℝ) := by
  rw [show m + 12 = (m + 6) + 6 from rfl, fact6, fact6]; push_cast; ring

/-- the diagonal `N = k + 1` of the closed form, after cancellation -/
lemma cf_diag (k : ℕ) :
    cf (k + 1) k = 2 ^ 12 * (4 * (k : ℝ) + 3) * (2 * (k : ℝ) + 1) * ((12 * k + 4).factorial : ℝ)
      / (2 ^ (22 * (k + 1)) * ((2 * k).factorial : ℝ) * ((6 * k + 2).factorial : ℝ)
        * ((2 * k + 1).factorial : ℝ) ^ 2) := by
  unfold cf
  rw [show 8 * k + 4 * (k + 1) = 12 * k + 4 by ring,
    show 4 * k + 4 - 2 * (k + 1) = 2 * (k + 1) by omega,
    show 8 * k + 8 - 4 * (k + 1) = 4 * (k + 1) by omega,
    show 4 * k + 2 * (k + 1) = 6 * k + 2 by ring,
    show 2 * (k + 1) - 2 = 2 * k by omega,
    show k + 1 - k - 1 = 0 by omega,
    show k + 1 + k = 2 * k + 1 by ring]
  have f1 : (((4 * (k + 1)).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f2 : (((2 * (k + 1)).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f3 : (((2 * k).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f4 : (((6 * k + 2).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f5 : (((2 * k + 1).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

/-- `ρ(k+1,k)` -/
lemma ρ_diag (k : ℕ) :
    ρ (k + 1) k = (4 * (k : ℝ) + 5) * (4 * (k : ℝ) + 7) ^ 2 * (12 * (k : ℝ) + 5) * (12 * (k : ℝ) + 7)
        * (4 * (k : ℝ) + 1) * (4 * (k : ℝ) + 3)
      / (131072 * ((k : ℝ) + 1) * (2 * (k : ℝ) + 1) * (4 * (k : ℝ) + 3) * (2 * (k : ℝ) + 2) ^ 2) := by
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have p2 : (0 : ℝ) < 2 * ((k : ℝ) + 1) - 1 := by linarith
  have p3 : (0 : ℝ) < 4 * ((k : ℝ) + 1) - 1 := by linarith
  have p4 : (0 : ℝ) < (k : ℝ) + 1 - k := by linarith
  have hρ : ρD ((k : ℝ) + 1) k ≠ 0 := by
    have : 0 < ρD ((k : ℝ) + 1) k := by unfold ρD; positivity
    exact this.ne'
  rw [ρ_eq]
  simp only [Nat.cast_add, Nat.cast_one]
  rw [div_eq_div_iff hρ (by positivity)]
  unfold ρN ρD
  ring

/-- `σ(k+2,k)` -/
lemma σ_diag (k : ℕ) :
    σ (k + 1 + 1) k = (2 * (k : ℝ) + 3) * ((12 * (k : ℝ) + 9) * (12 * (k : ℝ) + 11)
        * (12 * (k : ℝ) + 13) * (12 * (k : ℝ) + 15))
      / ((2 * (k : ℝ) + 1) * (2 * (k : ℝ) + 3) ^ 2 * ((4 * (k : ℝ) + 1) * (4 * (k : ℝ) + 3)
        * (4 * (k : ℝ) + 5) * (4 * (k : ℝ) + 7))) := by
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have q1 : (0 : ℝ) < 8 * (k : ℝ) - 4 * ((k : ℝ) + 1 + 1) + 9 := by linarith
  have q2 : (0 : ℝ) < 8 * (k : ℝ) - 4 * ((k : ℝ) + 1 + 1) + 11 := by linarith
  have q3 : (0 : ℝ) < 8 * (k : ℝ) - 4 * ((k : ℝ) + 1 + 1) + 13 := by linarith
  have q4 : (0 : ℝ) < 8 * (k : ℝ) - 4 * ((k : ℝ) + 1 + 1) + 15 := by linarith
  have hσ : σD ((k : ℝ) + 1 + 1) k ≠ 0 := by
    have : 0 < σD ((k : ℝ) + 1 + 1) k := by unfold σD; positivity
    exact this.ne'
  rw [σ_eq]
  simp only [Nat.cast_add, Nat.cast_one]
  rw [div_eq_div_iff hσ (by positivity)]
  unfold σN σD
  ring

/-- the diagonal: `a (k+1) k = cf (k+1) k` -/
lemma a_diag (k : ℕ) : a (k + 1) k = cf (k + 1) k := by
  induction k with
  | zero =>
    rw [a_zero_right, A_one, cf_diag]
    norm_num [Nat.factorial]
  | succ k ih =>
    rw [a_succ_right, a_succ_left (by omega) (by omega), ih, cf_diag, cf_diag, ρ_diag, σ_diag]
    rw [show 12 * (k + 1) + 4 = (12 * k + 4) + 12 by ring,
      show 2 * (k + 1) = 2 * k + 2 by ring,
      show 6 * (k + 1) + 2 = (6 * k + 2) + 6 by ring,
      show 2 * k + 2 + 1 = (2 * k + 1) + 2 by ring,
      fact12 (12 * k + 4), fact2 (2 * k), fact6 (6 * k + 2), fact2 (2 * k + 1),
      show 22 * (k + 1 + 1) = 22 * (k + 1) + 22 by ring, pow_add]
    have f3 : (((2 * k).factorial : ℕ) : ℝ) ≠ 0 := by positivity
    have f4 : (((6 * k + 2).factorial : ℕ) : ℝ) ≠ 0 := by positivity
    have f5 : (((2 * k + 1).factorial : ℕ) : ℝ) ≠ 0 := by positivity
    have f6 : (((12 * k + 4).factorial : ℕ) : ℝ) ≠ 0 := by positivity
    push_cast
    field_simp
    ring

end PartA

namespace PartA

/-- `ρ(N,k)` for `N = 2d + e + 1`, `k = d + e` -/
lemma ρ_col (d e : ℕ) :
    ρ (2 * d + e + 1) (d + e)
      = (8 * (d : ℝ) + 4 * e + 5) * (8 * (d : ℝ) + 4 * e + 7) ^ 2 * (16 * (d : ℝ) + 12 * e + 5)
          * (16 * (d : ℝ) + 12 * e + 7) * (4 * (e : ℝ) + 1) * (4 * (e : ℝ) + 3)
        / (131072 * (2 * (d : ℝ) + e + 1) * (4 * (d : ℝ) + 2 * e + 1) * (8 * (d : ℝ) + 4 * e + 3)
          * ((d : ℝ) + 1) ^ 2 * (3 * (d : ℝ) + 2 * e + 2) ^ 2) := by
  have hd : (0 : ℝ) ≤ d := Nat.cast_nonneg d
  have he : (0 : ℝ) ≤ e := Nat.cast_nonneg e
  have p2 : (0 : ℝ) < 2 * (2 * (d : ℝ) + e + 1) - 1 := by linarith
  have p3 : (0 : ℝ) < 4 * (2 * (d : ℝ) + e + 1) - 1 := by linarith
  have p4 : (0 : ℝ) < 2 * (d : ℝ) + e + 1 - ((d : ℝ) + e) := by linarith
  have hρ : ρD (2 * (d : ℝ) + e + 1) ((d : ℝ) + e) ≠ 0 := by
    have : 0 < ρD (2 * (d : ℝ) + e + 1) ((d : ℝ) + e) := by unfold ρD; positivity
    exact this.ne'
  rw [ρ_eq]
  simp only [Nat.cast_add, Nat.cast_mul, Nat.cast_one, Nat.cast_ofNat]
  rw [div_eq_div_iff hρ (by positivity)]
  unfold ρN ρD
  ring

/-- the column step of the closed form, `N = 2d + e + 1`, `k = d + e` -/
lemma cf_succ (d e : ℕ) :
    cf (2 * d + e + 1 + 1) (d + e) = cf (2 * d + e + 1) (d + e) * ρ (2 * d + e + 1) (d + e) := by
  rw [ρ_col]
  unfold cf
  rw [show 8 * (d + e) + 4 * (2 * d + e + 1 + 1) = (16 * d + 12 * e + 4) + 4 by ring,
    show 4 * (d + e) + 4 - 2 * (2 * d + e + 1 + 1) = 2 * e by omega,
    show 8 * (d + e) + 8 - 4 * (2 * d + e + 1 + 1) = 4 * e by omega,
    show 4 * (d + e) + 2 * (2 * d + e + 1 + 1) = (8 * d + 6 * e + 2) + 2 by ring,
    show 2 * (2 * d + e + 1 + 1) - 2 = (4 * d + 2 * e) + 2 by omega,
    show 2 * (2 * d + e + 1 + 1) = (4 * d + 2 * e) + 4 by ring,
    show 4 * (2 * d + e + 1 + 1) = (8 * d + 4 * e + 4) + 4 by ring,
    show 2 * d + e + 1 + 1 - (d + e) - 1 = d + 1 by omega,
    show 2 * d + e + 1 + 1 + (d + e) = (3 * d + 2 * e + 1) + 1 by ring,
    show 8 * (d + e) + 4 * (2 * d + e + 1) = 16 * d + 12 * e + 4 by ring,
    show 4 * (d + e) + 4 - 2 * (2 * d + e + 1) = 2 * e + 2 by omega,
    show 8 * (d + e) + 8 - 4 * (2 * d + e + 1) = 4 * e + 4 by omega,
    show 4 * (d + e) + 2 * (2 * d + e + 1) = 8 * d + 6 * e + 2 by ring,
    show 2 * (2 * d + e + 1) - 2 = 4 * d + 2 * e by omega,
    show 2 * (2 * d + e + 1) = (4 * d + 2 * e) + 2 by ring,
    show 4 * (2 * d + e + 1) = 8 * d + 4 * e + 4 by ring,
    show 2 * d + e + 1 - (d + e) - 1 = d by omega,
    show 2 * d + e + 1 + (d + e) = 3 * d + 2 * e + 1 by ring,
    fact4 (16 * d + 12 * e + 4), fact2 (8 * d + 6 * e + 2), fact2 (4 * d + 2 * e),
    fact4 (4 * d + 2 * e), fact4 (8 * d + 4 * e + 4), fact1 d, fact1 (3 * d + 2 * e + 1),
    fact2 (2 * e), fact4 (4 * e),
    show 22 * (2 * d + e + 1 + 1) = 22 * (2 * d + e + 1) + 22 by ring, pow_add]
  have f1 : (((8 * d + 4 * e + 4).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f2 : (((16 * d + 12 * e + 4).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f3 : (((2 * e).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f4 : (((4 * d + 2 * e).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f5 : (((8 * d + 6 * e + 2).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f6 : (((4 * e).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f7 : ((d.factorial : ℕ) : ℝ) ≠ 0 := by positivity
  have f8 : (((3 * d + 2 * e + 1).factorial : ℕ) : ℝ) ≠ 0 := by positivity
  push_cast
  field_simp
  ring

/-- the column step of the closed form, for `k < N`, `N + 1 ≤ 2k + 2` -/
lemma cf_succ' {N k : ℕ} (hk : k < N) (hN : N + 1 ≤ 2 * k + 2) : cf (N + 1) k = cf N k * ρ N k := by
  obtain ⟨d, e, rfl, rfl⟩ : ∃ d e, N = 2 * d + e + 1 ∧ k = d + e :=
    ⟨N - k - 1, 2 * k + 1 - N, by omega, by omega⟩
  exact cf_succ d e

end PartA

theorem a_eq_cf {N k : ℕ} (hk : k < N) (hN : N ≤ 2 * k + 2) : a N k = cf N k := by
  have key : ∀ d : ℕ, d ≤ k + 1 → a (k + 1 + d) k = cf (k + 1 + d) k := by
    intro d
    induction d with
    | zero => intro _; exact a_diag k
    | succ d ih =>
      intro hd
      rw [show k + 1 + (d + 1) = (k + 1 + d) + 1 by ring, a_succ_left (by omega) (by omega),
        ih (by omega), cf_succ' (by omega) (by omega)]
  have h := key (N - k - 1) (by omega)
  rwa [show k + 1 + (N - k - 1) = N by omega] at h

end PiSeries.SunPi2
