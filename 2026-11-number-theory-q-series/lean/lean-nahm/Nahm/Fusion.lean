/-
Nahm/Fusion.lean — Lemma 2 of `proofs/E8-E6.md` (fusion of two species) for infinite sums.

For a finite index type `ι`, an exponent function `E : (ι → ℕ) → ℕ` (in powers of `X`), a base
`Q = X^d` and two distinct indices `a ≠ b`, the Nahm-type sum
`F_E = Σ_{n : ι → ℕ} X^{E n} / ∏_i (Q;Q)_{n_i}` does not change if one adjoins a new index `0 = a + b`
with exponent `E°(m, c) = E(m + c(e_a + e_b)) + d m_a m_b` (`fuseE`).  The proof regroups the sum
according to `n = m + c(e_a + e_b)` and uses the Durfee rectangle identity.
-/
import Nahm.Multisum

namespace Nahm

open PowerSeries Finset KRTW

open scoped PowerSeries.WithPiTopology

variable {K : Type*} [Field K] (Q : K⟦X⟧) {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The summand `X^{E n} / ∏_i (Q;Q)_{n_i}`. -/
noncomputable def nterm (E : (ι → ℕ) → ℕ) (n : ι → ℕ) : K⟦X⟧ := X ^ E n * ∏ i, iq Q (n i)

/-- `m + c (e_a + e_b)`. -/
def shift (a b : ι) (c : ℕ) (m : ι → ℕ) : ι → ℕ := fun i => m i + if i = a ∨ i = b then c else 0

/-- The exponent of the fused sum, on `Option ι` (the new index is `none`):
`E°(n) = E(m + c(e_a + e_b)) + d m_a m_b` with `m = n ∘ some`, `c = n none`. -/
def fuseE (d : ℕ) (E : (ι → ℕ) → ℕ) (a b : ι) (n : Option ι → ℕ) : ℕ :=
  E (shift a b (n none) (fun i => n (some i))) + d * (n (some a) * n (some b))

omit [Fintype ι] in
lemma shift_unshift (a b : ι) (c : ℕ) (n : ι → ℕ) (ha : c ≤ n a) (hb : c ≤ n b) :
    shift a b c (fun i => n i - if i = a ∨ i = b then c else 0) = n := by
  funext i
  simp only [shift]
  split_ifs with h
  · rcases h with rfl | rfl <;> omega
  · simp

omit [Fintype ι] in
lemma unshift_shift (a b : ι) (c : ℕ) (m : ι → ℕ) :
    (fun i => shift a b c m i - if i = a ∨ i = b then c else 0) = m := by
  funext i; simp [shift]

omit [Fintype ι] in
lemma shift_injective (a b : ι) (c : ℕ) : Function.Injective (shift a b c) := by
  intro m m' h
  rw [← unshift_shift a b c m, ← unshift_shift a b c m', h]

omit [Fintype ι] in
lemma le_fuseE (d : ℕ) (E : (ι → ℕ) → ℕ) (hE : ∀ n i, n i ≤ E n) (a b : ι) :
    ∀ n i, n i ≤ fuseE d E a b n := by
  intro n i
  have key : ∀ j, n (some j) ≤ E (shift a b (n none) (fun i => n (some i))) := fun j =>
    le_trans (by simp [shift]) (hE _ j)
  rcases i with _ | j
  · exact le_trans (le_trans (by simp [shift]) (hE _ a)) (Nat.le_add_right _ _)
  · exact le_trans (key j) (Nat.le_add_right _ _)

variable {Q}

/-- The fusion lemma. -/
theorem fusion [Nil Q] [TopologicalSpace K] [DiscreteTopology K] (d : ℕ) (hQ : Q = X ^ d)
    (E : (ι → ℕ) → ℕ) (hE : ∀ n i, n i ≤ E n) (a b : ι) (hab : a ≠ b) :
    ∑' n : Option ι → ℕ, nterm Q (fuseE d E a b) n = ∑' n : ι → ℕ, nterm Q E n := by
  -- the summand as a function of `(c, m)`
  set F : ℕ → (ι → ℕ) → K⟦X⟧ := fun c m =>
    X ^ (E (shift a b c m) + d * (m a * m b)) * (iq Q c * ∏ i, iq Q (m i)) with hF
  -- after reindexing `m ↦ n = m + c(e_a + e_b)`
  set G : ℕ → (ι → ℕ) → K⟦X⟧ := fun c n =>
    if c ≤ n a ∧ c ≤ n b then F c (fun i => n i - if i = a ∨ i = b then c else 0) else 0 with hG
  have hGF : ∀ c m, G c (shift a b c m) = F c m := by
    intro c m
    simp only [hG]
    rw [if_pos ⟨by simp [shift], by simp [shift]⟩, unshift_shift]
  -- divisibility
  have hEsup : ∀ n : ι → ℕ, univ.sup n ≤ E n := fun n => Finset.sup_le (fun i _ => hE n i)
  have hdivG : ∀ c n, (X : K⟦X⟧) ^ (c ⊔ univ.sup n) ∣ G c n := by
    intro c n
    simp only [hG]
    split_ifs with h
    · simp only [hF]
      rw [shift_unshift a b c n h.1 h.2]
      refine Dvd.dvd.mul_right (pow_dvd_pow X ?_) _
      have h1 := hE n a
      have h2 := hEsup n
      omega
    · exact dvd_zero _
  have hfin : ∀ N, {p : ℕ × (ι → ℕ) | p.1 ⊔ univ.sup p.2 ≤ N}.Finite := by
    intro N
    apply ((range (N + 1)) ×ˢ (Fintype.piFinset fun _ : ι => range (N + 1))).finite_toSet.subset
    rintro ⟨c, n⟩ h
    simp only [Set.mem_ofPred_eq, sup_le_iff] at h
    simp only [coe_product, Set.mem_prod, mem_coe, mem_range, Fintype.mem_piFinset]
    exact ⟨by omega, fun i => Nat.lt_succ_of_le (le_trans (Finset.le_sup (mem_univ i)) h.2)⟩
  have hGs : Summable (Function.uncurry G) :=
    summable_of_dvd (fun p => p.1 ⊔ univ.sup p.2) hfin (fun p => hdivG p.1 p.2)
  have hGc : ∀ c, Summable (G c) := by
    intro c
    apply summable_of_dvd (fun n => univ.sup n) _ (fun n => (pow_dvd_pow X le_sup_right).trans (hdivG c n))
    intro N
    apply (Fintype.piFinset fun _ : ι => range (N + 1)).finite_toSet.subset
    intro n h
    simp only [Set.mem_ofPred_eq] at h
    simp only [mem_coe, Fintype.mem_piFinset, mem_range]
    exact fun i => Nat.lt_succ_of_le (le_trans (Finset.le_sup (mem_univ i)) h)
  have hGn : ∀ n, Summable fun c => G c n := by
    intro n
    apply summable_of_ne_finset_zero (s := range (n a + 1))
    intro c hc
    simp only [mem_range, not_lt] at hc
    simp only [hG]
    rw [if_neg (by omega)]
  -- Step 1: the left side as a double sum over `(c, m)`
  have e1 : ∑' n : Option ι → ℕ, nterm Q (fuseE d E a b) n = ∑' p : ℕ × (ι → ℕ), F p.1 p.2 := by
    rw [← (Equiv.piOptionEquivProd (α := ι) (β := fun _ => ℕ)).symm.tsum_eq]
    refine tsum_congr (fun p => ?_)
    simp only [nterm, fuseE, Equiv.piOptionEquivProd_symm_apply, hF, Fintype.prod_option]
  -- Step 2: reindex each slice
  have e2 : ∀ c, ∑' m, F c m = ∑' n, G c n := by
    intro c
    have hsupp : Function.support (G c) ⊆ Set.range (shift a b c) := by
      intro n hn
      simp only [Function.mem_support, hG] at hn
      split_ifs at hn with h
      · exact ⟨_, shift_unshift a b c n h.1 h.2⟩
      · exact absurd rfl hn
    rw [← (shift_injective a b c).tsum_eq hsupp]
    exact tsum_congr (fun m => (hGF c m).symm)
  have hdivF : ∀ c m, (X : K⟦X⟧) ^ (c ⊔ univ.sup m) ∣ F c m := by
    intro c m
    simp only [hF]
    refine Dvd.dvd.mul_right (pow_dvd_pow X ?_) _
    have h1 := hE (shift a b c m) a
    have h2 := hEsup (shift a b c m)
    have h3 : univ.sup m ≤ univ.sup (shift a b c m) :=
      Finset.sup_mono_fun (fun i _ => by simp [shift])
    simp only [shift, true_or, if_true] at h1
    omega
  have hFs : Summable (Function.uncurry F) :=
    summable_of_dvd (fun p => p.1 ⊔ univ.sup p.2) hfin (fun p => hdivF p.1 p.2)
  have hFc : ∀ c, Summable (F c) := by
    intro c
    apply summable_of_dvd (fun m => univ.sup m) _ (fun m => (pow_dvd_pow X le_sup_right).trans (hdivF c m))
    intro N
    apply (Fintype.piFinset fun _ : ι => range (N + 1)).finite_toSet.subset
    intro n h
    simp only [Set.mem_ofPred_eq] at h
    simp only [mem_coe, Fintype.mem_piFinset, mem_range]
    exact fun i => Nat.lt_succ_of_le (le_trans (Finset.le_sup (mem_univ i)) h)
  have e3 : ∑' p : ℕ × (ι → ℕ), F p.1 p.2 = ∑' n, ∑' c, G c n := by
    rw [show ∑' p : ℕ × (ι → ℕ), F p.1 p.2 = ∑' p : ℕ × (ι → ℕ), Function.uncurry F p from rfl,
      hFs.tsum_prod' hFc]
    simp only [Function.uncurry_apply_pair]
    rw [tsum_congr e2]
    exact (hGs.tsum_comm' hGc hGn).symm
  rw [e1, e3]
  -- Step 3: the inner sums by the Durfee identity
  refine tsum_congr (fun n => ?_)
  rw [tsum_eq_sum (s := range (min (n a) (n b) + 1))]
  · have hd := durfee Q (n a) (n b)
    have hprod : ∀ f : ι → K⟦X⟧, ∏ i, f i = f a * f b * ∏ i ∈ (univ.erase a).erase b, f i := by
      intro f
      rw [mul_assoc, mul_prod_erase _ _ (mem_erase.mpr ⟨hab.symm, mem_univ b⟩), mul_prod_erase _ _ (mem_univ a)]
    rw [nterm, hprod, show X ^ E n * (iq Q (n a) * iq Q (n b) * ∏ i ∈ (univ.erase a).erase b, iq Q (n i))
      = X ^ E n * (∏ i ∈ (univ.erase a).erase b, iq Q (n i)) * (iq Q (n a) * iq Q (n b)) by ring, hd, mul_sum]
    refine sum_congr rfl (fun c hc => ?_)
    have hc' : c ≤ n a ∧ c ≤ n b := by
      have := Nat.lt_succ_iff.mp (mem_range.mp hc); exact ⟨by omega, by omega⟩
    simp only [hG, hF]
    rw [if_pos hc', shift_unshift a b c n hc'.1 hc'.2, hprod]
    have hrest : ∏ i ∈ (univ.erase a).erase b, iq Q (n i - if i = a ∨ i = b then c else 0)
        = ∏ i ∈ (univ.erase a).erase b, iq Q (n i) := by
      refine prod_congr rfl (fun i hi => ?_)
      simp only [mem_erase] at hi
      rw [if_neg (by tauto), Nat.sub_zero]
    rw [hrest]
    simp only [true_or, or_true, if_true]
    rw [hQ, ← pow_mul, pow_add]
    ring
  · intro c hc
    simp only [mem_range, not_lt] at hc
    simp only [hG]
    rw [if_neg (by omega)]

end Nahm
