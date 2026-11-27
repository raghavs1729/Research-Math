/-
Tau/Sieve.lean — Steps 2 and 3 of PROOF.md: periodicity of the orbit modulo `q` and the congruence
sieve on the index `n` modulo `M`.

The finite computations are done on natural numbers (`stepN`, `iterN`, `isPow`) so that the kernel
evaluates them with its built-in (GMP) arithmetic; `cast_iterN` links them to `T (ZMod q)`.  The
witness primes need not be prime: the argument works for every modulus `q > 1`.
-/
import Tau.Basic

namespace Tau

/-! ### Computable layer -/

/-- One step `(a, b) ↦ (b, 3b − a mod q)` on residues `0 ≤ a, b < q`. -/
def stepN (q : ℕ) (p : ℕ × ℕ) : ℕ × ℕ := (p.2, (3 * p.2 + (q - p.1)) % q)

/-- `r` steps. -/
def iterN (q r : ℕ) (p : ℕ × ℕ) : ℕ × ℕ := (stepN q)^[r] p

/-- Reduction of an integer pair modulo `q`, as a pair of naturals in `[0, q)`. -/
def redN (q : ℕ) (s : ℤ × ℤ) : ℕ × ℕ := ((s.1 % (q : ℤ)).toNat, (s.2 % (q : ℤ)).toNat)

/-- `a` is an `e`-th power modulo `q` (brute force over all residues). -/
def isPow (e q a : ℕ) : Bool := (List.range q).any fun x => x ^ e % q == a

/-- The class `n ≡ r (mod M)` of the orbit of `s` is killed by `q`: the orbit of `s mod q` has period
dividing `M`, and at step `r` the first coordinate is not an 11th power or the second is not a
square modulo `q`. -/
def classOK (q M r : ℕ) (s : ℤ × ℤ) : Bool :=
  1 < q && iterN q M (redN q s) == redN q s &&
    (!(isPow 11 q (iterN q r (redN q s)).1) || !(isPow 2 q (iterN q r (redN q s)).2))

/-- For every seed and every class `r mod M`, some modulus in `Qs` kills the class. -/
def certOK (S : List (ℤ × ℤ)) (Qs : List ℕ) (M : ℕ) : Bool :=
  S.all fun s => (List.range M).all fun r => Qs.any fun q => classOK q M r s

/-! ### Link with `ZMod q` -/

/-- A pair of naturals viewed in `ZMod q × ZMod q`. -/
def castP (q : ℕ) (p : ℕ × ℕ) : ZMod q × ZMod q := ((p.1 : ZMod q), (p.2 : ZMod q))

lemma stepN_lt {q : ℕ} (hq : 0 < q) {p : ℕ × ℕ} (h : p.1 < q ∧ p.2 < q) :
    (stepN q p).1 < q ∧ (stepN q p).2 < q :=
  ⟨h.2, Nat.mod_lt _ hq⟩

lemma castP_stepN {q : ℕ} {p : ℕ × ℕ} (h : p.1 < q) :
    castP q (stepN q p) = T (ZMod q) (castP q p) := by
  simp only [castP, stepN, T_apply, ZMod.natCast_mod, Nat.cast_add, Nat.cast_mul,
    Nat.cast_sub h.le, ZMod.natCast_self, Prod.mk.injEq, true_and]
  push_cast; ring

lemma iterN_succ (q r : ℕ) (p : ℕ × ℕ) : iterN q (r + 1) p = stepN q (iterN q r p) :=
  Function.iterate_succ_apply' _ _ _

lemma iterN_lt {q : ℕ} (hq : 0 < q) (r : ℕ) {p : ℕ × ℕ} (h : p.1 < q ∧ p.2 < q) :
    (iterN q r p).1 < q ∧ (iterN q r p).2 < q := by
  induction r with
  | zero => exact h
  | succ k ih => rw [iterN_succ]; exact stepN_lt hq ih

lemma castP_iterN {q : ℕ} (hq : 0 < q) (r : ℕ) {p : ℕ × ℕ} (h : p.1 < q ∧ p.2 < q) :
    castP q (iterN q r p) = (T (ZMod q) ^ r) (castP q p) := by
  induction r with
  | zero => rfl
  | succ k ih =>
    rw [iterN_succ, castP_stepN (iterN_lt hq k h).1, ih, pow_succ', Equiv.Perm.mul_apply]

lemma redN_lt {q : ℕ} (hq : 0 < q) (s : ℤ × ℤ) : (redN q s).1 < q ∧ (redN q s).2 < q := by
  have hq' : (0 : ℤ) < q := by exact_mod_cast hq
  have h1 := Int.emod_nonneg s.1 hq'.ne'
  have h2 := Int.emod_lt_of_pos s.1 hq'
  have h3 := Int.emod_nonneg s.2 hq'.ne'
  have h4 := Int.emod_lt_of_pos s.2 hq'
  simp only [redN]
  omega

lemma castP_redN {q : ℕ} (hq : 0 < q) (s : ℤ × ℤ) :
    castP q (redN q s) = pmap (Int.castRingHom (ZMod q)) s := by
  have hq' : (q : ℤ) ≠ 0 := by exact_mod_cast hq.ne'
  have key : ∀ a : ℤ, (((a % (q : ℤ)).toNat : ℕ) : ZMod q) = (a : ZMod q) := fun a => by
    rw [← Int.cast_natCast, Int.toNat_of_nonneg (Int.emod_nonneg _ hq'), ZMod.intCast_mod]
  simp only [castP, redN, pmap, eq_intCast, key]

lemma isPow_of {e q a : ℕ} [NeZero q] (ha : a < q) (x : ZMod q) (h : (a : ZMod q) = x ^ e) :
    isPow e q a = true := by
  simp only [isPow, List.any_eq_true, List.mem_range, beq_iff_eq]
  refine ⟨x.val, ZMod.val_lt x, ?_⟩
  have h1 : ((x.val ^ e % q : ℕ) : ZMod q) = (a : ZMod q) := by
    rw [ZMod.natCast_mod, Nat.cast_pow, ZMod.natCast_zmod_val, h]
  have h2 := congrArg ZMod.val h1
  rw [ZMod.val_natCast, ZMod.val_natCast, Nat.mod_mod, Nat.mod_eq_of_lt ha] at h2
  exact h2

/-- Periodicity: if the orbit of `x₀` under `T (ZMod q)` returns after `M` steps, then
`T ^ n x₀ = T ^ (n mod M) x₀` for every `n : ℤ`. -/
lemma zpow_apply_emod {q M : ℕ} (hM : 0 < M) (x₀ : ZMod q × ZMod q)
    (hper : (T (ZMod q) ^ M) x₀ = x₀) (n : ℤ) :
    (T (ZMod q) ^ n) x₀ = (T (ZMod q) ^ (n % (M : ℤ)).toNat) x₀ := by
  have hM' : (M : ℤ) ≠ 0 := by exact_mod_cast hM.ne'
  have hnn : 0 ≤ n % (M : ℤ) := Int.emod_nonneg _ hM'
  have hper' : (T (ZMod q) ^ (M : ℤ)) x₀ = x₀ := by rwa [zpow_natCast]
  conv_lhs => rw [← Int.emod_add_mul_ediv n (M : ℤ)]
  rw [zpow_add, zpow_mul, Equiv.Perm.mul_apply,
    Equiv.Perm.zpow_apply_eq_self_of_apply_eq_self hper', ← zpow_natCast,
    Int.toNat_of_nonneg hnn]

/-- **The sieve.** If the certificate checks, then no point `T ^ n s` (for a seed `s ∈ S`, `n : ℤ`)
has the form `(x ^ 11, y ^ 2)`. -/
theorem sieve {S : List (ℤ × ℤ)} {Qs : List ℕ} {M : ℕ} (hM : 0 < M) (h : certOK S Qs M = true)
    {s : ℤ × ℤ} (hs : s ∈ S) (n : ℤ) (x y : ℤ) : (T ℤ ^ n) s ≠ (x ^ 11, y ^ 2) := by
  intro heq
  set r := (n % (M : ℤ)).toNat with hr
  have hrM : r < M := by
    have hM' : (0 : ℤ) < M := by exact_mod_cast hM
    have h1 := Int.emod_nonneg n hM'.ne'
    have h2 := Int.emod_lt_of_pos n hM'
    omega
  simp only [certOK, List.all_eq_true, List.mem_range, List.any_eq_true] at h
  obtain ⟨q, -, hq⟩ := h s hs r hrM
  simp only [classOK, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq, Bool.or_eq_true,
    Bool.not_eq_true'] at hq
  obtain ⟨⟨hq1, hper⟩, hkill⟩ := hq
  have hq0 : 0 < q := by omega
  have : NeZero q := ⟨hq0.ne'⟩
  set f := Int.castRingHom (ZMod q)
  have hbd := redN_lt hq0 s
  -- reduce `heq` modulo `q`
  have hmod := congrArg (pmap f) heq
  rw [pmap_zpow, ← castP_redN hq0] at hmod
  have hper' : (T (ZMod q) ^ M) (castP q (redN q s)) = castP q (redN q s) := by
    rw [← castP_iterN hq0 M hbd, hper]
  rw [zpow_apply_emod hM _ hper' n, ← hr, ← castP_iterN hq0 r hbd] at hmod
  have hlt := iterN_lt hq0 r hbd
  simp only [castP, pmap, map_pow, f, eq_intCast, Prod.mk.injEq] at hmod
  rcases hkill with hk | hk
  · rw [isPow_of hlt.1 _ hmod.1] at hk; exact Bool.noConfusion hk
  · rw [isPow_of hlt.2 _ hmod.2] at hk; exact Bool.noConfusion hk

end Tau
