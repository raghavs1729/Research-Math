/-
KRTW/Limit.lean — the limit `R → ∞` of Theorem 1: the first equality of Theorem 2 of
`proofs/KRTW-general.md`, `KRTW_k(x) = Mel_k(x)`, for the recursively defined sums.
-/
import KRTW.Theorem1

namespace KRTW

open PowerSeries Finset

open scoped PowerSeries.WithPiTopology

set_option linter.unusedSectionVars false

variable {K : Type*} [Field K] [TopologicalSpace K] [DiscreteTopology K]

/-- One level of the KRTW sum with an arbitrary weight `ω r` in place of `[R;r]`:
`Σ_{r ≥ 0} Σ_{s ≥ B} ω(r) x^{r+2s} q^{(r+s)²+s²+s} / (q²;q²)_{s-B} · T_j(r,s;x)`. -/
noncomputable def Lv (ω : ℕ → K⟦X⟧) (j B : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑' r, ∑' b, X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
    (ω r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x)

/-- The `2k`-fold KRTW sum with the `x`-grading, `KRTW_k(x)`, defined one level at a time. -/
noncomputable def Kf : ℕ → K⟦X⟧ → K⟦X⟧
  | 0, _ => 1
  | j + 1, x => Lv (iq X) j 0 x

lemma Lv_exp_ge_right (r b B : ℕ) : b ≤ r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B)) := by
  have := le_tri (2 * (b + B))
  omega

lemma Lv_exp_ge_left (r b B : ℕ) : r ≤ r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B)) := by
  have : r ≤ r ^ 2 := Nat.le_self_pow (by norm_num) r
  omega

lemma Lv_inner_summable (ω : ℕ → K⟦X⟧) (j B : ℕ) (x : K⟦X⟧) (r : ℕ) :
    Summable (fun b => X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
      (ω r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x)) :=
  summable_nat (fun b => dvd_X_pow_mul (Lv_exp_ge_right r b B) _)

lemma Lv_outer_summable (ω : ℕ → K⟦X⟧) (j B : ℕ) (x : K⟦X⟧) :
    Summable (fun r => ∑' b, X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
      (ω r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x)) :=
  summable_nat (fun r => X_pow_dvd_tsum (fun b => dvd_X_pow_mul (Lv_exp_ge_left r b B) _))

lemma Tf_succ_eq_Lv (j R B : ℕ) (x : K⟦X⟧) : Tf (j + 1) R B x = Lv (qbinom X R) j B x := by
  rw [Tf_succ, Lv]
  symm
  apply tsum_eq_sum
  intro r hr
  have hr' : R < r := by
    simp only [mem_range, not_lt] at hr
    omega
  simp [qbinom_of_lt X hr']

/-- `A_{j+1}(R;x) ≡ KRTW_{j+1}(x) mod q^{R+1}` (in the generality of an arbitrary lower bound `B`). -/
lemma X_pow_dvd_Lv_sub (R j B : ℕ) (x : K⟦X⟧) :
    (X : K⟦X⟧) ^ (R + 1) ∣ Lv (iq X) j B x - Lv (qbinom X R) j B x := by
  refine X_pow_dvd_tsum_sub_tsum (Lv_outer_summable _ j B x) (Lv_outer_summable _ j B x)
    (fun r => ?_)
  refine X_pow_dvd_tsum_sub_tsum (Lv_inner_summable _ j B x r) (Lv_inner_summable _ j B x r)
    (fun b => ?_)
  have e : X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
        (iq X r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x)
      - X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
        (qbinom X R r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x)
      = (X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) * (iq X r - qbinom X R r)) *
        (x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x) := by ring
  rw [e]
  apply Dvd.dvd.mul_right
  rcases Nat.lt_or_ge R r with h | h
  · exact dvd_X_pow_mul (le_trans (Nat.succ_le_of_lt h) (Lv_exp_ge_left r b B)) _
  · have h1 : (X : K⟦X⟧) ^ (R + 1) = X ^ r * X ^ (R - r + 1) := by
      rw [← pow_add]; congr 1; omega
    rw [h1]
    exact mul_dvd_mul (pow_dvd_pow X (Lv_exp_ge_left r b B)) (X_pow_dvd_iq_sub_qbinom X h)

lemma Mel_succ (j : ℕ) (x : K⟦X⟧) :
    Mel (j + 1) x = ∑' m, X ^ tri m * (x ^ m * iq X m * Hf j m x) := by
  simp only [Mel, Hf_succ, add_zero]

/-- `B_{j+1}(R;x) ≡ Mel_{j+1}(x) mod q^{R+1}`. -/
lemma X_pow_dvd_Bf_sub_Mel (j R : ℕ) (x : K⟦X⟧) :
    (X : K⟦X⟧) ^ (R + 1) ∣ Bf j R x - Mel (j + 1) x := by
  rw [Mel_succ, Bf]
  refine X_pow_dvd_tsum_sub_tsum (summable_nat (fun m => dvd_X_pow_mul (le_tri m) _))
    (summable_nat (fun m => dvd_X_pow_mul (le_tri m) _)) (fun m => ?_)
  have e : X ^ tri m * (x ^ m * iq X m * Einf (R + m) x * Hf j m x)
      - X ^ tri m * (x ^ m * iq X m * Hf j m x)
      = (X ^ tri m * (x ^ m * iq X m * Hf j m x)) * (Einf (R + m) x - 1) := by ring
  rw [e]
  apply Dvd.dvd.mul_left
  exact (pow_dvd_pow X (by omega)).trans (X_pow_dvd_Einf_sub_one (R + m) x)

/-- **Theorem 2, first equality** (recursive form): `KRTW_k(x) = Mel_k(x)` for every `k` and every
power series `x`. -/
theorem Kf_eq_Mel : ∀ (k : ℕ) (x : K⟦X⟧), Kf k x = Mel k x
  | 0, x => rfl
  | j + 1, x => by
    apply eq_of_forall_X_pow_dvd_sub
    intro n
    have h1 := X_pow_dvd_Lv_sub n j 0 x
    have h2 := X_pow_dvd_Bf_sub_Mel j n x
    rw [← Tf_succ_eq_Lv, thm1] at h1
    have h3 : (X : K⟦X⟧) ^ (n + 1) ∣ Kf (j + 1) x - Mel (j + 1) x := dvd_sub_trans h1 h2
    exact (pow_dvd_pow X (Nat.le_succ n)).trans h3

end KRTW
