/-
KRTW/Summation.lean — infinite sums of formal power series (`q`-adic, i.e. coefficientwise) and the
few rearrangement lemmas that the proof needs.

`K⟦X⟧` carries the product topology of the discrete topology on the coefficient field `K`
(`[TopologicalSpace K] [DiscreteTopology K]` are explicit hypotheses, so that nothing here depends on a
choice of instance; the final statements do not mention any topology).
-/
import KRTW.Basic

namespace KRTW

open PowerSeries Finset

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K] [TopologicalSpace K] [DiscreteTopology K]

lemma coeff_tsum {ι : Type*} {f : ι → K⟦X⟧} (hf : Summable f) (n : ℕ) :
    coeff n (∑' i, f i) = ∑' i, coeff n (f i) :=
  hf.map_tsum (coeff n) (WithPiTopology.continuous_coeff K n)

omit [DiscreteTopology K] in
/-- A family whose terms are divisible by higher and higher powers of `X` is summable. -/
lemma summable_of_dvd {ι : Type*} {f : ι → K⟦X⟧} (g : ι → ℕ) (hg : ∀ n, {i | g i ≤ n}.Finite)
    (h : ∀ i, (X : K⟦X⟧) ^ g i ∣ f i) : Summable f := by
  rw [WithPiTopology.summable_iff_summable_coeff]
  intro n
  apply summable_of_hasFiniteSupport
  apply (hg n).subset
  intro i hi
  simp only [Function.mem_support] at hi
  by_contra hlt
  simp only [Set.mem_ofPred_eq, not_le] at hlt
  exact hi (X_pow_dvd_iff.mp (h i) n hlt)

lemma summable_nat {f : ℕ → K⟦X⟧} (h : ∀ i, (X : K⟦X⟧) ^ i ∣ f i) : Summable f :=
  summable_of_dvd (fun i => i) (fun n => Set.finite_Iic n) h

lemma summable_nat2 {f : ℕ × ℕ → K⟦X⟧} (h : ∀ p, (X : K⟦X⟧) ^ (max p.1 p.2) ∣ f p) :
    Summable f := by
  apply summable_of_dvd (fun p => max p.1 p.2) _ h
  intro n
  apply (Set.finite_Iic ((n, n) : ℕ × ℕ)).subset
  intro p hp
  simp only [Set.mem_ofPred_eq] at hp
  exact ⟨le_trans (le_max_left _ _) hp, le_trans (le_max_right _ _) hp⟩

lemma X_pow_dvd_tsum {ι : Type*} {f : ι → K⟦X⟧} {n : ℕ} (h : ∀ i, (X : K⟦X⟧) ^ n ∣ f i) :
    (X : K⟦X⟧) ^ n ∣ ∑' i, f i := by
  by_cases hf : Summable f
  · rw [X_pow_dvd_iff]
    intro m hm
    rw [coeff_tsum hf]
    have : ∀ i, coeff m (f i) = 0 := fun i => X_pow_dvd_iff.mp (h i) m hm
    simp [this]
  · rw [tsum_eq_zero_of_not_summable hf]
    exact dvd_zero _

/-- Truncation of a sum: modulo `X^n`, only the terms not divisible by `X^n` matter. -/
lemma X_pow_dvd_tsum_sub_sum {ι : Type*} {f : ι → K⟦X⟧} (hf : Summable f) (s : Finset ι) {n : ℕ}
    (h : ∀ i ∉ s, (X : K⟦X⟧) ^ n ∣ f i) : (X : K⟦X⟧) ^ n ∣ ∑' i, f i - ∑ i ∈ s, f i := by
  rw [X_pow_dvd_iff]
  intro m hm
  rw [map_sub, coeff_tsum hf, map_sum, sub_eq_zero]
  apply tsum_eq_sum
  intro i hi
  exact X_pow_dvd_iff.mp (h i hi) m hm

lemma X_pow_dvd_tsum_sub_tsum {ι : Type*} {f g : ι → K⟦X⟧} (hf : Summable f) (hg : Summable g)
    {n : ℕ} (h : ∀ i, (X : K⟦X⟧) ^ n ∣ f i - g i) :
    (X : K⟦X⟧) ^ n ∣ ∑' i, f i - ∑' i, g i := by
  rw [← hf.tsum_sub hg]
  exact X_pow_dvd_tsum h

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma eq_of_forall_X_pow_dvd_sub {f g : K⟦X⟧} (h : ∀ n, (X : K⟦X⟧) ^ n ∣ f - g) : f = g := by
  ext m
  have := X_pow_dvd_iff.mp (h (m + 1)) m (Nat.lt_succ_self m)
  rwa [map_sub, sub_eq_zero] at this

omit [TopologicalSpace K] [DiscreteTopology K] in
lemma dvd_sub_trans {a b c d : K⟦X⟧} (h1 : d ∣ a - b) (h2 : d ∣ b - c) : d ∣ a - c := by
  have : a - c = (a - b) + (b - c) := by ring
  rw [this]; exact dvd_add h1 h2

/-- Triangular rearrangement: `Σ_b Σ_m F(b, m + σ b) = Σ_M Σ_{b : σ b ≤ M} F(b, M)`. -/
lemma tsum_tsum_shift (σ : ℕ → ℕ) (hσ : ∀ b, b ≤ σ b) (F : ℕ → ℕ → K⟦X⟧)
    (hF : ∀ b M, (X : K⟦X⟧) ^ M ∣ F b M) :
    ∑' b, ∑' m, F b (m + σ b)
      = ∑' M, ∑ b ∈ range (M + 1), if σ b ≤ M then F b M else 0 := by
  set f : ℕ → ℕ → K⟦X⟧ := fun b M => if σ b ≤ M then F b M else 0 with hf
  have hdvd : ∀ b M, (X : K⟦X⟧) ^ M ∣ f b M := by
    intro b M
    simp only [hf]
    split_ifs
    · exact hF b M
    · exact dvd_zero _
  have hsum : Summable (Function.uncurry f) := by
    apply summable_nat2
    rintro ⟨b, M⟩
    simp only [Function.uncurry_apply_pair, hf]
    split_ifs with h
    · have : max b M = M := max_eq_right (le_trans (hσ b) h)
      rw [this]; exact hF b M
    · exact dvd_zero _
  have h1 : ∀ b, Summable (f b) := fun b => summable_nat (hdvd b)
  have h2 : ∀ M, Summable (fun b => f b M) := by
    intro M
    apply summable_of_ne_finset_zero (s := range (M + 1))
    intro b hb
    simp only [mem_range, not_lt] at hb
    simp only [hf]
    rw [if_neg]
    have := hσ b
    omega
  have hcomm := Summable.tsum_comm' hsum h1 h2
  have e1 : ∑' b, ∑' m, F b (m + σ b) = ∑' b, ∑' M, f b M := by
    apply tsum_congr
    intro b
    have h := (h1 b).sum_add_tsum_nat_add (σ b)
    rw [← h, sum_eq_zero, zero_add]
    · apply tsum_congr
      intro m
      simp [hf]
    · intro i hi
      simp only [mem_range] at hi
      simp only [hf]
      rw [if_neg (by omega)]
  rw [e1, ← hcomm]
  apply tsum_congr
  intro M
  apply tsum_eq_sum
  intro b hb
  simp only [mem_range, not_lt] at hb
  simp only [hf]
  rw [if_neg]
  have := hσ b
  omega

end KRTW
