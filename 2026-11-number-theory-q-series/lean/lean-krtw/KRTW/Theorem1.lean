/-
KRTW/Theorem1.lean — Theorem 1 of `proofs/KRTW-general.md` (the `R`-refined identity
`A_j(R;x) = B_j(R;x)`), by induction on `j`.

The multiple sums are defined recursively, one level at a time; `KRTW/Bridge.lean` identifies them with
the explicit sums over index tuples.
-/
import KRTW.Lemma4

namespace KRTW

open PowerSeries Finset

open scoped PowerSeries.WithPiTopology

set_option linter.unusedSectionVars false

variable {K : Type*} [Field K] [TopologicalSpace K] [DiscreteTopology K]

lemma tsum_mul_left' (a : K⟦X⟧) (f : ℕ → K⟦X⟧) (hf : ∀ m, (X : K⟦X⟧) ^ m ∣ f m) :
    a * ∑' m, f m = ∑' m, a * f m :=
  ((summable_nat hf).tsum_mul_left a).symm

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma dvd_X_pow_mul {n e : ℕ} (h : n ≤ e) (r : K⟦X⟧) : (X : K⟦X⟧) ^ n ∣ X ^ e * r :=
  Dvd.dvd.mul_right (pow_dvd_pow X h) r

/-! ### The Melzer-type chain sums -/

/-- `H_j(B;x) = Σ_{B ≤ N_1 ≤ … ≤ N_j} ∏_i x^{N_i} q^{N_i(N_i+1)/2} / (q)_{N_i - N_{i-1}}` (`N_0 := B`),
defined by peeling off the smallest variable `N_1 = m + B`. -/
noncomputable def Hf : ℕ → ℕ → K⟦X⟧ → K⟦X⟧
  | 0, _, _ => 1
  | j + 1, B, x => ∑' m, X ^ tri (m + B) * (x ^ (m + B) * iq X m * Hf j (m + B) x)

lemma Hf_zero (B : ℕ) (x : K⟦X⟧) : Hf 0 B x = 1 := rfl

lemma Hf_succ (j B : ℕ) (x : K⟦X⟧) :
    Hf (j + 1) B x = ∑' m, X ^ tri (m + B) * (x ^ (m + B) * iq X m * Hf j (m + B) x) := rfl

/-- Shifting all summation variables by `c`:
`H_j(B+c;x) = x^{jc} q^{j c(c+1)/2} H_j(B; xq^c)`. -/
lemma Hf_shift : ∀ (j B c : ℕ) (x : K⟦X⟧),
    Hf j (B + c) x = X ^ (j * tri c) * (x ^ (j * c) * Hf j B (x * X ^ c))
  | 0, B, c, x => by simp [Hf_zero]
  | j + 1, B, c, x => by
    rw [Hf_succ, Hf_succ]
    have hd : ∀ m, (X : K⟦X⟧) ^ m ∣
        X ^ tri (m + B) * ((x * X ^ c) ^ (m + B) * iq X m * Hf j (m + B) (x * X ^ c)) :=
      fun m => dvd_X_pow_mul (le_trans (Nat.le_add_right m B) (le_tri _)) _
    have e : X ^ ((j + 1) * tri c) * (x ^ ((j + 1) * c) *
        ∑' m, X ^ tri (m + B) * ((x * X ^ c) ^ (m + B) * iq X m * Hf j (m + B) (x * X ^ c)))
        = (X ^ ((j + 1) * tri c) * x ^ ((j + 1) * c)) *
        ∑' m, X ^ tri (m + B) * ((x * X ^ c) ^ (m + B) * iq X m * Hf j (m + B) (x * X ^ c)) := by
      ring
    rw [e, tsum_mul_left' _ _ hd]
    refine tsum_congr (fun m => ?_)
    have e1 : m + (B + c) = (m + B) + c := by ring
    rw [e1, Hf_shift j (m + B) c x, tri_add]
    ring

/-- `Mel_j(x) = Σ_{M_1 ≥ … ≥ M_j ≥ 0} x^{ΣM_i} q^{Σ M_i(M_i+1)/2} / ∏ (q)_{M_i - M_{i+1}}`. -/
noncomputable def Mel (j : ℕ) (x : K⟦X⟧) : K⟦X⟧ := Hf j 0 x

/-- `B_{j+1}(R;x)` of Theorem 1 (the Melzer-type sum with the factor `(xq^{R+M+1};q)_∞`, `M` the
smallest summation variable). -/
noncomputable def Bf (j R : ℕ) (x : K⟦X⟧) : K⟦X⟧ :=
  ∑' m, X ^ tri m * (x ^ m * iq X m * Einf (R + m) x * Hf j m x)

/-! ### The KRTW chain sums -/

/-- `T_j(R,B;x) = Σ_{R ≥ r_1 ≥ … ≥ r_j ≥ 0, B ≤ s_1 ≤ … ≤ s_j} ((q)_R/(q)_{R-r_1})
  ∏_i x^{r_i+2s_i} q^{(r_i+s_i)²+s_i²+s_i} / ((q)_{r_i-r_{i+1}} (q²;q²)_{s_i-s_{i-1}})` (`s_0 := B`),
defined by peeling off the first level `(r_1, s_1) = (r, b + B)`.  Note
`(r+s)²+s²+s = r² + 2rs + tri(2s)`. -/
noncomputable def Tf : ℕ → ℕ → ℕ → K⟦X⟧ → K⟦X⟧
  | 0, _, _, _ => 1
  | j + 1, R, B, x => ∑ r ∈ range (R + 1), ∑' b,
      X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
        (qbinom X R r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x)

lemma Tf_zero (R B : ℕ) (x : K⟦X⟧) : Tf 0 R B x = 1 := rfl

lemma Tf_succ (j R B : ℕ) (x : K⟦X⟧) :
    Tf (j + 1) R B x = ∑ r ∈ range (R + 1), ∑' b,
      X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
        (qbinom X R r * x ^ (r + 2 * (b + B)) * iq (X ^ 2) b * Tf j r (b + B) x) := rfl

/-- Shifting all `s`-variables by `c`:
`T_j(R,B+c;x) = x^{2cj} q^{j(2c²+c)} T_j(R,B; xq^{2c})`. -/
lemma Tf_shift : ∀ (j R B c : ℕ) (x : K⟦X⟧),
    Tf j R (B + c) x = X ^ (j * tri (2 * c)) * (x ^ (2 * c * j) * Tf j R B (x * X ^ (2 * c)))
  | 0, R, B, c, x => by simp [Tf_zero]
  | j + 1, R, B, c, x => by
    rw [Tf_succ, Tf_succ, mul_sum, mul_sum]
    refine sum_congr rfl (fun r _ => ?_)
    have hd : ∀ b, (X : K⟦X⟧) ^ b ∣ X ^ (r ^ 2 + 2 * r * (b + B) + tri (2 * (b + B))) *
        (qbinom X R r * (x * X ^ (2 * c)) ^ (r + 2 * (b + B)) * iq (X ^ 2) b *
          Tf j r (b + B) (x * X ^ (2 * c))) := by
      intro b
      apply dvd_X_pow_mul
      have := le_tri (2 * (b + B))
      omega
    rw [tsum_mul_left' _ _ hd, tsum_mul_left' _ _ (fun b => Dvd.dvd.mul_left (hd b) _)]
    refine tsum_congr (fun b => ?_)
    have e1 : b + (B + c) = (b + B) + c := by ring
    have e2 : 2 * (b + B + c) = 2 * (b + B) + 2 * c := by ring
    rw [e1, Tf_shift j r (b + B) c x, e2, tri_add]
    ring

/-! ### Theorem 1 -/

/-- **Theorem 1** (`A_j(R;x) = B_j(R;x)` for `j ≥ 1`, here with `j` replaced by `j + 1`). -/
theorem thm1 : ∀ (j R : ℕ) (x : K⟦X⟧), Tf (j + 1) R 0 x = Bf j R x
  | 0, R, x => by
    have h := starInf R x
    have e1 : Tf 1 R 0 x = Wf Ginf R x := by
      rw [Tf_succ]
      simp only [Wf, Ginf, add_zero, Tf_zero, mul_one]
      refine sum_congr rfl (fun r _ => ?_)
      have hd : ∀ b, (X : K⟦X⟧) ^ b ∣
          X ^ tri (2 * b) * ((x * X ^ r) ^ (2 * b) * iq (X ^ 2) b) := by
        intro b
        apply dvd_X_pow_mul
        have := le_tri (2 * b)
        omega
      have e : X ^ r ^ 2 * (qbinom X R r * x ^ r *
          ∑' b, X ^ tri (2 * b) * ((x * X ^ r) ^ (2 * b) * iq (X ^ 2) b))
          = (X ^ r ^ 2 * (qbinom X R r * x ^ r)) *
          ∑' b, X ^ tri (2 * b) * ((x * X ^ r) ^ (2 * b) * iq (X ^ 2) b) := by ring
      rw [e, tsum_mul_left' _ _ hd]
      refine tsum_congr (fun b => ?_)
      ring
    have e2 : Bf 0 R x = Uinf R x := by
      simp only [Bf, Uinf, Hf_zero, mul_one]
    rw [e1, e2, h]
  | j + 1, R, x => by
    have ih := thm1 j
    set Ψ : ℕ → K⟦X⟧ := fun M => X ^ tri M * (x ^ M * Hf j M x) with hΨ
    have hΨdvd : ∀ M, (X : K⟦X⟧) ^ M ∣ Ψ M := fun M => dvd_X_pow_mul (le_tri M) _
    -- the right side
    have hB : Bf (j + 1) R x = ∑' M, Ψ M * Uss R M x := by
      set F' : ℕ → ℕ → K⟦X⟧ := fun L M =>
        Ψ M * (X ^ tri L * (x ^ L * iq X L * iq X (M - L) * Einf (R + L) x)) with hF'
      have e1 : Bf (j + 1) R x = ∑' L, ∑' m, F' L (m + L) := by
        unfold Bf
        refine tsum_congr (fun L => ?_)
        rw [Hf_succ]
        have hd : ∀ m, (X : K⟦X⟧) ^ m ∣
            X ^ tri (m + L) * (x ^ (m + L) * iq X m * Hf j (m + L) x) :=
          fun m => dvd_X_pow_mul (le_trans (Nat.le_add_right m L) (le_tri _)) _
        have e : X ^ tri L * (x ^ L * iq X L * Einf (R + L) x *
            ∑' m, X ^ tri (m + L) * (x ^ (m + L) * iq X m * Hf j (m + L) x))
            = (X ^ tri L * (x ^ L * iq X L * Einf (R + L) x)) *
            ∑' m, X ^ tri (m + L) * (x ^ (m + L) * iq X m * Hf j (m + L) x) := by ring
        rw [e, tsum_mul_left' _ _ hd]
        refine tsum_congr (fun m => ?_)
        simp only [hF', hΨ, Nat.add_sub_cancel]
        ring
      rw [e1, tsum_tsum_shift (fun b => b) (fun b => le_refl b) F'
        (fun b M => Dvd.dvd.mul_right (hΨdvd M) _)]
      refine tsum_congr (fun M => ?_)
      unfold Uss
      rw [mul_sum]
      refine sum_congr rfl (fun L hL => ?_)
      rw [if_pos (Nat.lt_succ_iff.mp (mem_range.mp hL))]
    -- the left side
    have hA : Tf (j + 1 + 1) R 0 x
        = ∑' M, Ψ M * Wf (fun y => Einf M y * GM M y) R x := by
      rw [Tf_succ]
      simp only [add_zero]
      have hr : ∀ r : ℕ, ∑' b, X ^ (r ^ 2 + 2 * r * b + tri (2 * b)) *
            (qbinom X R r * x ^ (r + 2 * b) * iq (X ^ 2) b * Tf (j + 1) r b x)
          = ∑' M, Ψ M * (X ^ r ^ 2 * (qbinom X R r * x ^ r *
              (Einf M (x * X ^ r) * GM M (x * X ^ r)))) := by
        intro r
        set F : ℕ → ℕ → K⟦X⟧ := fun b M =>
          Ψ M * (X ^ r ^ 2 * (qbinom X R r * x ^ r * (Einf M (x * X ^ r) *
            (X ^ tri (2 * b) * ((x * X ^ r) ^ (2 * b) * iq (X ^ 2) b * iq X (M - 2 * b))))))
          with hF
        have e1 : ∀ b : ℕ, X ^ (r ^ 2 + 2 * r * b + tri (2 * b)) *
            (qbinom X R r * x ^ (r + 2 * b) * iq (X ^ 2) b * Tf (j + 1) r b x)
            = ∑' m, F b (m + 2 * b) := by
          intro b
          have hT : Tf (j + 1) r b x
              = X ^ ((j + 1) * tri (2 * b)) * (x ^ (2 * b * (j + 1)) * Bf j r (x * X ^ (2 * b))) := by
            have h := Tf_shift (j + 1) r 0 b x
            rw [zero_add] at h
            rw [h, ih]
          rw [hT]
          unfold Bf
          have hd : ∀ m, (X : K⟦X⟧) ^ m ∣ X ^ tri m * ((x * X ^ (2 * b)) ^ m * iq X m *
              Einf (r + m) (x * X ^ (2 * b)) * Hf j m (x * X ^ (2 * b))) :=
            fun m => dvd_X_pow_mul (le_tri m) _
          have e : X ^ (r ^ 2 + 2 * r * b + tri (2 * b)) *
              (qbinom X R r * x ^ (r + 2 * b) * iq (X ^ 2) b *
                (X ^ ((j + 1) * tri (2 * b)) * (x ^ (2 * b * (j + 1)) *
                  ∑' m, X ^ tri m * ((x * X ^ (2 * b)) ^ m * iq X m *
                    Einf (r + m) (x * X ^ (2 * b)) * Hf j m (x * X ^ (2 * b))))))
              = (X ^ (r ^ 2 + 2 * r * b + tri (2 * b)) *
                (qbinom X R r * x ^ (r + 2 * b) * iq (X ^ 2) b *
                  (X ^ ((j + 1) * tri (2 * b)) * x ^ (2 * b * (j + 1))))) *
                  ∑' m, X ^ tri m * ((x * X ^ (2 * b)) ^ m * iq X m *
                    Einf (r + m) (x * X ^ (2 * b)) * Hf j m (x * X ^ (2 * b))) := by ring
          rw [e, tsum_mul_left' _ _ hd]
          refine tsum_congr (fun m => ?_)
          simp only [hF, hΨ, Nat.add_sub_cancel]
          have e3 : r + m + 2 * b = m + 2 * b + r := by ring
          rw [Hf_shift j m (2 * b) x, Einf_mul_pow, Einf_mul_pow, e3, tri_add]
          ring
        rw [tsum_congr e1, tsum_tsum_shift (fun b => 2 * b) (fun b => by omega) F
          (fun b M => Dvd.dvd.mul_right (hΨdvd M) _)]
        refine tsum_congr (fun M => ?_)
        unfold GM
        simp only [mul_sum, mul_ite, mul_zero, hF]
      rw [sum_congr rfl (fun r _ => hr r),
        ← Summable.tsum_finsetSum (fun r _ => summable_nat
          (fun M => Dvd.dvd.mul_right (hΨdvd M) _))]
      refine tsum_congr (fun M => ?_)
      simp only [Wf]
      rw [mul_sum]
    rw [hA, hB]
    refine tsum_congr (fun M => ?_)
    rw [star2]

end KRTW
