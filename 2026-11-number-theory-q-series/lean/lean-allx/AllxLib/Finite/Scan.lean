/-
# Soundness of the integer scan for the log-weighted sum `Θ₁`

With `cw p = w4 p * log p` (`w4 p = −χ₄(p)`) and `Θ₁(t) = ∑_{p ≤ t} cw(p) (t − p)`
(`Coupling.Theta1 cw`), the scan of `AllxLib.Finite.Defs` keeps integers `up, um, vp, vm` with

  `(up − um) log 2 ≤ 256 ∑_{p < n} cw p`,   `(vp − vm) log 2 ≤ 256 ∑_{p < n} cw(p) (n − 1 − p)`

(`Allx.Inv`).  These are *inequalities between real numbers*; the integers themselves are
never claimed to be "correct".  They are maintained because at each prime `n` the scan has
checked `2^a ≤ n^256 < 2^(a+1)`, i.e. `a log 2 ≤ 256 log n < (a+1) log 2`, and primality is
read off the sieve mask (`Allx.sieve_bit_iff`).  At each `n ≥ 4` the scan checks `vm < vp`
for the updated `V`, which gives `Θ₁(n) > 0`.

Main results: `Allx.scan_sound`, `Allx.Good.step` (one chunk), `Allx.good_zero`, and the
passage from integers to reals `Allx.theta1_pos_real`, `Allx.theta1_nonneg_real`.
-/
import AllxLib.Coupling.Convexity
import AllxLib.Finite.Sieve

open Real Finset

namespace Allx

/-- The log-weighted coefficient `cw p = −χ₄(p) log p`. -/
noncomputable def cw (p : ℕ) : ℝ := Coupling.w4 p * Real.log p

/-- `∑_{p < n} cw p`. -/
noncomputable def ThLt (n : ℕ) : ℝ := ∑ p ∈ range n with p.Prime, cw p

/-- `∑_{p < n} cw(p) (n − 1 − p)`. -/
noncomputable def Vsem (n : ℕ) : ℝ := ∑ p ∈ range n with p.Prime, cw p * ((n : ℝ) - 1 - p)

lemma ThLt_zero : ThLt 0 = 0 := by simp [ThLt]

lemma Vsem_zero : Vsem 0 = 0 := by simp [Vsem]

lemma ThLt_succ (n : ℕ) : ThLt (n + 1) = ThLt n + if n.Prime then cw n else 0 := by
  unfold ThLt
  rw [sum_filter, sum_range_succ, ← sum_filter]

lemma Vsem_succ (n : ℕ) : Vsem (n + 1) = Vsem n + ThLt n := by
  unfold Vsem ThLt
  rw [sum_filter, sum_range_succ, ← sum_filter, ← sum_add_distrib]
  have : (if n.Prime then cw n * (((n + 1 : ℕ) : ℝ) - 1 - n) else 0) = 0 := by
    split_ifs <;> simp
  rw [this, add_zero]
  refine sum_congr rfl fun p _ => ?_
  push_cast; ring

lemma theta1_nat (n : ℕ) : Coupling.Theta1 cw n = Vsem (n + 1) := by
  unfold Coupling.Theta1 Vsem
  rw [Nat.floor_natCast, Nat.range_succ_eq_Icc_zero]
  refine sum_congr rfl fun p _ => ?_
  push_cast; ring

/-! ### The weight `w4` -/

lemma chi4_nat (p : ℕ) : (ZMod.χ₄ (p : ZMod 4) : ℤ) =
    if p % 4 = 1 then 1 else if p % 4 = 3 then -1 else 0 := by
  rw [ZMod.χ₄_nat_mod_four]
  have : p % 4 < 4 := Nat.mod_lt _ (by norm_num)
  interval_cases h : p % 4 <;> rfl

lemma w4_eq (p : ℕ) :
    Coupling.w4 p = (if p % 4 = 3 then 1 else 0) - (if p % 4 = 1 then 1 else 0) := by
  unfold Coupling.w4
  rw [chi4_nat]
  split_ifs <;> simp_all

lemma cw_three {n : ℕ} (h : n % 4 = 3) : cw n = Real.log n := by
  rw [cw, w4_eq]; simp [h]

lemma cw_one {n : ℕ} (h : n % 4 = 1) : cw n = -Real.log n := by
  rw [cw, w4_eq]; simp [h]

lemma cw_other {n : ℕ} (h3 : n % 4 ≠ 3) (h1 : n % 4 ≠ 1) : cw n = 0 := by
  rw [cw, w4_eq]; simp [h3, h1]

/-! ### Logarithm bounds from integer comparisons -/

lemma log_lower {n a : ℕ} (h : 2 ^ a ≤ n ^ 256) :
    (a : ℝ) * Real.log 2 ≤ 256 * Real.log n := by
  have h' : ((2 : ℝ)) ^ a ≤ (n : ℝ) ^ 256 := by exact_mod_cast h
  have := Real.log_le_log (by positivity) h'
  rw [Real.log_pow, Real.log_pow] at this
  simpa using this

lemma log_upper {n a : ℕ} (hn : 0 < n) (h : n ^ 256 + 1 ≤ 2 ^ (a + 1)) :
    256 * Real.log n < ((a : ℝ) + 1) * Real.log 2 := by
  have h' : (n : ℝ) ^ 256 < ((2 : ℝ)) ^ (a + 1) := by exact_mod_cast Nat.lt_of_succ_le h
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have := Real.log_lt_log (by positivity) h'
  rw [Real.log_pow, Real.log_pow] at this
  simpa using this

/-! ### The invariant and the soundness of the scan -/

/-- The semantic invariant of the scan at position `n` (before `n` is processed). -/
def Inv (n up um vp vm : ℕ) : Prop :=
  ((up : ℝ) - um) * Real.log 2 ≤ 256 * ThLt n ∧ ((vp : ℝ) - vm) * Real.log 2 ≤ 256 * Vsem n

theorem scan_zero (C n0 n a up um vp vm : ℕ) :
    scan C n0 0 n a up um vp vm = (true, a, up, um, vp, vm) := rfl

theorem scan_succ (C n0 fuel n a up um vp vm : ℕ) :
    scan C n0 (fuel + 1) n a up um vp vm =
      bif (bif Nat.ble 4 n then Nat.ble (vm + um + 1) (vp + up) else true) then
        bif (bif Nat.ble 2 n then Nat.beq (Nat.land (Nat.shiftRight C (n - n0)) 1) 0
            else false) then
          bif (bif Nat.ble (2 ^ liftAll (n ^ 256) a) (n ^ 256) then
              Nat.ble (n ^ 256 + 1) (2 ^ (liftAll (n ^ 256) a + 1)) else false) then
            bif Nat.beq (n % 4) 3 then
              scan C n0 fuel (n + 1) (liftAll (n ^ 256) a) (up + liftAll (n ^ 256) a) um
                (vp + up) (vm + um)
            else bif Nat.beq (n % 4) 1 then
              scan C n0 fuel (n + 1) (liftAll (n ^ 256) a) up (um + liftAll (n ^ 256) a + 1)
                (vp + up) (vm + um)
            else scan C n0 fuel (n + 1) (liftAll (n ^ 256) a) up um (vp + up) (vm + um)
          else scanFail
        else scan C n0 fuel (n + 1) a up um (vp + up) (vm + um)
      else scanFail := rfl

lemma scanFail_ne (s : ℕ × ℕ × ℕ × ℕ × ℕ) : scanFail ≠ (true, s) := by
  intro h
  have := congrArg Prod.fst h
  simp [scanFail] at this

/-- **Soundness of the scan.** `C` is any mask whose bits decide primality on `[n0, n0+K)`. -/
theorem scan_sound (C n0 K : ℕ)
    (hC : ∀ n, n0 ≤ n → n < n0 + K → 2 ≤ n →
      (Nat.beq (Nat.land (Nat.shiftRight C (n - n0)) 1) 0 = true ↔ n.Prime)) :
    ∀ (fuel n a up um vp vm : ℕ) (s : ℕ × ℕ × ℕ × ℕ × ℕ), n0 ≤ n → n + fuel ≤ n0 + K →
      Inv n up um vp vm → scan C n0 fuel n a up um vp vm = (true, s) →
      Inv (n + fuel) s.2.1 s.2.2.1 s.2.2.2.1 s.2.2.2.2 ∧
        ∀ k, n ≤ k → k < n + fuel → 4 ≤ k → 0 < Coupling.Theta1 cw k := by
  intro fuel
  induction fuel with
  | zero =>
    intro n a up um vp vm s _ _ hI h
    rw [scan_zero] at h
    obtain rfl := (Prod.mk.inj h).2
    exact ⟨hI, fun k h1 h2 => by omega⟩
  | succ fuel ih =>
    intro n a up um vp vm s hn0 hnK hI h
    have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
    rw [scan_succ] at h
    -- the positivity check
    cases hc1 : (bif Nat.ble 4 n then Nat.ble (vm + um + 1) (vp + up) else true) with
    | false => rw [hc1, cond_false] at h; exact absurd h (scanFail_ne s)
    | true =>
    rw [hc1, cond_true] at h
    have hV : (((vp + up : ℕ) : ℝ) - ((vm + um : ℕ) : ℝ)) * Real.log 2 ≤ 256 * Vsem (n + 1) := by
      rw [Vsem_succ]; push_cast; nlinarith [hI.1, hI.2]
    have hpos : 4 ≤ n → 0 < Coupling.Theta1 cw n := by
      intro h4
      rw [Nat.ble_eq_true_of_le h4, cond_true] at hc1
      have hlt : vm + um < vp + up := Nat.lt_of_succ_le (Nat.le_of_ble_eq_true hc1)
      have hlt' : ((vm + um : ℕ) : ℝ) < ((vp + up : ℕ) : ℝ) := by exact_mod_cast hlt
      rw [theta1_nat]
      have : 0 < (((vp + up : ℕ) : ℝ) - ((vm + um : ℕ) : ℝ)) * Real.log 2 :=
        mul_pos (by linarith) hlog2
      linarith
    -- how to conclude from the recursive call
    have key : ∀ a' up' um' : ℕ,
        ((up' : ℝ) - um') * Real.log 2 ≤ 256 * ThLt (n + 1) →
        scan C n0 fuel (n + 1) a' up' um' (vp + up) (vm + um) = (true, s) →
        Inv (n + (fuel + 1)) s.2.1 s.2.2.1 s.2.2.2.1 s.2.2.2.2 ∧
          ∀ k, n ≤ k → k < n + (fuel + 1) → 4 ≤ k → 0 < Coupling.Theta1 cw k := by
      intro a' up' um' hU hs
      obtain ⟨h1, h2⟩ := ih (n + 1) a' up' um' (vp + up) (vm + um) s (by omega) (by omega)
        ⟨hU, hV⟩ hs
      refine ⟨by rwa [show n + (fuel + 1) = n + 1 + fuel by omega], fun k hk1 hk2 hk4 => ?_⟩
      rcases Nat.eq_or_lt_of_le hk1 with rfl | hlt
      · exact hpos hk4
      · exact h2 k hlt (by omega) hk4
    -- the primality test
    cases hc2 : (bif Nat.ble 2 n then Nat.beq (Nat.land (Nat.shiftRight C (n - n0)) 1) 0
        else false) with
    | false =>
      rw [hc2, cond_false] at h
      have hnp : ¬ n.Prime := by
        intro hp
        have h2 := hp.two_le
        rw [Nat.ble_eq_true_of_le h2, cond_true] at hc2
        rw [(hC n hn0 (by omega) h2).2 hp] at hc2
        exact absurd hc2 (by simp)
      refine key a up um ?_ h
      rw [ThLt_succ, if_neg hnp, add_zero]; exact hI.1
    | true =>
    rw [hc2, cond_true] at h
    have h2 : 2 ≤ n := by
      by_contra hlt
      have : Nat.ble 2 n = false := by
        cases hb : Nat.ble 2 n with
        | false => rfl
        | true => exact absurd (Nat.le_of_ble_eq_true hb) hlt
      rw [this, cond_false] at hc2
      exact absurd hc2 (by simp)
    rw [Nat.ble_eq_true_of_le h2, cond_true] at hc2
    have hp : n.Prime := (hC n hn0 (by omega) h2).1 hc2
    -- the check of the exponent
    cases hc3 : (bif Nat.ble (2 ^ liftAll (n ^ 256) a) (n ^ 256) then
        Nat.ble (n ^ 256 + 1) (2 ^ (liftAll (n ^ 256) a + 1)) else false) with
    | false => rw [hc3, cond_false] at h; exact absurd h (scanFail_ne s)
    | true =>
    rw [hc3, cond_true] at h
    generalize liftAll (n ^ 256) a = a' at h hc3
    have hlo : 2 ^ a' ≤ n ^ 256 := by
      cases hb : Nat.ble (2 ^ a') (n ^ 256) with
      | false => rw [hb, cond_false] at hc3; exact absurd hc3 (by simp)
      | true => exact Nat.le_of_ble_eq_true hb
    rw [Nat.ble_eq_true_of_le hlo, cond_true] at hc3
    have hhi : n ^ 256 + 1 ≤ 2 ^ (a' + 1) := Nat.le_of_ble_eq_true hc3
    have hL := log_lower hlo
    have hU := log_upper (by omega : 0 < n) hhi
    have hTh : ThLt (n + 1) = ThLt n + cw n := by rw [ThLt_succ, if_pos hp]
    cases hr3 : Nat.beq (n % 4) 3 with
    | true =>
      rw [hr3, cond_true] at h
      refine key a' (up + a') um ?_ h
      rw [hTh, cw_three (Nat.eq_of_beq_eq_true hr3)]
      push_cast; nlinarith [hI.1]
    | false =>
    rw [hr3, cond_false] at h
    have hn3 : n % 4 ≠ 3 := fun h' => by rw [h'] at hr3; exact absurd hr3 (by decide)
    cases hr1 : Nat.beq (n % 4) 1 with
    | true =>
      rw [hr1, cond_true] at h
      refine key a' up (um + a' + 1) ?_ h
      rw [hTh, cw_one (Nat.eq_of_beq_eq_true hr1)]
      push_cast; nlinarith [hI.1]
    | false =>
      rw [hr1, cond_false] at h
      have hn1 : n % 4 ≠ 1 := fun h' => by rw [h'] at hr1; exact absurd hr1 (by decide)
      refine key a' up um ?_ h
      rw [hTh, cw_other hn3 hn1, add_zero]; exact hI.1

/-! ### Chunks -/

/-- The state `s = (a, up, um, vp, vm)` at position `n` satisfies the invariant, and
`Θ₁(k) > 0` for all integers `4 ≤ k < n`. -/
def Good (n : ℕ) (s : ℕ × ℕ × ℕ × ℕ × ℕ) : Prop :=
  Inv n s.2.1 s.2.2.1 s.2.2.2.1 s.2.2.2.2 ∧ ∀ k, 4 ≤ k → k < n → 0 < Coupling.Theta1 cw k

theorem good_zero : Good 0 (0, 0, 0, 0, 0) := by
  refine ⟨⟨?_, ?_⟩, fun k _ h => by omega⟩
  · simp [ThLt_zero]
  · simp [Vsem_zero]

theorem runChunk_eq (n0 K B a up um vp vm : ℕ) :
    runChunk n0 K B a up um vp vm =
      bif Nat.ble (n0 + K) (B * B) then scan (sieve n0 K B) n0 K n0 a up um vp vm
      else scanFail := rfl

/-- One kernel-checked chunk transports `Good` from `n0` to `n0 + K`. -/
theorem Good.step {n0 K B a up um vp vm : ℕ} {s : ℕ × ℕ × ℕ × ℕ × ℕ}
    (hG : Good n0 (a, up, um, vp, vm)) (hc : runChunk n0 K B a up um vp vm = (true, s)) :
    Good (n0 + K) s := by
  rw [runChunk_eq] at hc
  cases hb : Nat.ble (n0 + K) (B * B) with
  | false => rw [hb, cond_false] at hc; exact absurd hc (scanFail_ne s)
  | true =>
    rw [hb, cond_true] at hc
    have hB := Nat.le_of_ble_eq_true hb
    obtain ⟨h1, h2⟩ := scan_sound (sieve n0 K B) n0 K
      (fun n h0 hK h2 => sieve_bit_iff h0 hK h2 hB) K n0 a up um vp vm s le_rfl le_rfl hG.1 hc
    refine ⟨h1, fun k hk4 hk => ?_⟩
    rcases Nat.lt_or_ge k n0 with hlt | hge
    · exact hG.2 k hk4 hlt
    · exact h2 k hge hk hk4

/-! ### From integers to reals -/

/-- On `[n, n+1]` the function `Θ₁` is affine (its kinks are at integers). -/
lemma theta1_affine (c : ℕ → ℝ) (n : ℕ) {x : ℝ} (h1 : (n : ℝ) ≤ x) (h2 : x ≤ n + 1) :
    Coupling.Theta1 c x =
      (n + 1 - x) * Coupling.Theta1 c n + (x - n) * Coupling.Theta1 c (n + 1 : ℕ) := by
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hx0 : 0 ≤ x := hn0.trans h1
  rw [Coupling.Theta1_eq_sum_max c hx0 (x := ((n + 1 : ℕ) : ℝ)) (by push_cast; linarith),
    Coupling.Theta1_eq_sum_max c hn0 (x := ((n + 1 : ℕ) : ℝ)) (by push_cast; linarith),
    Coupling.Theta1_eq_sum_max c (Nat.cast_nonneg _) (x := ((n + 1 : ℕ) : ℝ)) le_rfl,
    mul_sum, mul_sum, ← sum_add_distrib]
  refine sum_congr rfl fun p hp => ?_
  have hp' : p ≤ n + 1 := by
    have := (mem_Icc.1 (mem_filter.1 hp).1).2
    rwa [Nat.floor_natCast] at this
  push_cast
  rcases Nat.lt_or_ge n p with hlt | hle
  · have hpn : p = n + 1 := by omega
    subst hpn
    push_cast
    rw [max_eq_right (by linarith), max_eq_right (by linarith), max_eq_right (by linarith)]
    ring
  · have : (p : ℝ) ≤ n := by exact_mod_cast hle
    rw [max_eq_left (by linarith), max_eq_left (by linarith), max_eq_left (by linarith)]
    ring

lemma cw_two : cw 2 = 0 := cw_other (by decide) (by decide)

lemma ThLt_three : ThLt 3 = 0 := by
  have h0 : ThLt 1 = 0 := by rw [ThLt_succ, ThLt_zero]; simp [Nat.not_prime_zero]
  have h1 : ThLt 2 = 0 := by rw [ThLt_succ, h0]; simp [Nat.not_prime_one]
  rw [ThLt_succ, h1, cw_two]; simp

lemma theta1_two : Coupling.Theta1 cw (2 : ℕ) = 0 := by
  have h0 : ThLt 1 = 0 := by rw [ThLt_succ, ThLt_zero]; simp [Nat.not_prime_zero]
  have h1 : ThLt 2 = 0 := by rw [ThLt_succ, h0]; simp [Nat.not_prime_one]
  rw [theta1_nat, Vsem_succ, Vsem_succ, Vsem_succ, Vsem_zero, ThLt_zero, h0, h1]; simp

lemma theta1_three : Coupling.Theta1 cw (3 : ℕ) = 0 := by
  have h := theta1_two
  rw [theta1_nat] at h
  rw [theta1_nat, Vsem_succ, h, ThLt_three, add_zero]

/-- If `Θ₁(k) > 0` at the integers `4 ≤ k ≤ N`, then `Θ₁(t) > 0` for all real `3 < t ≤ N`. -/
theorem theta1_pos_real (N : ℕ) (hN : ∀ k : ℕ, 4 ≤ k → k ≤ N → 0 < Coupling.Theta1 cw k)
    {x : ℝ} (hx3 : 3 < x) (hxN : x ≤ N) : 0 < Coupling.Theta1 cw x := by
  set n := ⌊x⌋₊ with hn
  have hx0 : 0 ≤ x := by linarith
  have hnx : (n : ℝ) ≤ x := Nat.floor_le hx0
  have hxn : x < n + 1 := Nat.lt_floor_add_one x
  have hn3 : 3 ≤ n := Nat.le_floor (by norm_num; linarith)
  have hnN : n ≤ N := Nat.floor_le_of_le hxN
  rcases eq_or_lt_of_le hnx with heq | hlt
  · rw [← heq]
    have : 4 ≤ n := by
      by_contra h
      have : n = 3 := by omega
      rw [this] at heq; push_cast at heq; linarith
    exact hN n this hnN
  · have hn1 : n + 1 ≤ N := by
      by_contra h
      have : n = N := by omega
      rw [this] at hlt; linarith
    have hpos : 0 < Coupling.Theta1 cw (n + 1 : ℕ) := hN (n + 1) (by omega) hn1
    have hnn : 0 ≤ Coupling.Theta1 cw n := by
      rcases Nat.eq_or_lt_of_le hn3 with h | h
      · rw [← h, theta1_three]
      · exact (hN n (by omega) hnN).le
    rw [theta1_affine cw n hnx hxn.le]
    have : 0 < x - n := by linarith
    have : 0 ≤ n + 1 - x := by linarith
    positivity

/-- Under the same hypothesis, `Θ₁(t) ≥ 0` for all real `2 ≤ t ≤ N`. -/
theorem theta1_nonneg_real (N : ℕ) (hN : ∀ k : ℕ, 4 ≤ k → k ≤ N → 0 < Coupling.Theta1 cw k)
    {x : ℝ} (hx2 : 2 ≤ x) (hxN : x ≤ N) : 0 ≤ Coupling.Theta1 cw x := by
  rcases le_or_gt x 3 with h3 | h3
  · have := theta1_affine cw 2 (x := x) (by norm_num; linarith) (by norm_num; linarith)
    rw [this, theta1_two, show ((2 + 1 : ℕ)) = 3 from rfl, theta1_three]; simp
  · exact (theta1_pos_real N hN h3 hxN).le

end Allx
