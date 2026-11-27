/-
# The finite range of the log-weighted race mod 4, kernel-checked

With `c(p) = −χ₄(p) log p` (`Coupling.w4 p * Real.log p`) and
`Θ₁(t) = ∑_{p ≤ t} c(p) (t − p)` (`Coupling.Theta1`):

  `Θ₁(t) > 0` for every real `3 < t ≤ 10⁷`,   `Θ₁(t) ≥ 0` for every real `2 ≤ t ≤ 10⁷`.

Proof: `Θ₁` is affine on each `[n, n+1]`, so it suffices that `Θ₁(k) > 0` at the integers
`4 ≤ k ≤ 10⁷` (`Allx.theta1_cw_nat_pos`, from `AllxLib.Finite.Chain`).  That is established
by an integer scan evaluated by the Lean kernel (`decide +kernel`, no `native_decide`) in
chunks of `10⁴` integers; see `AllxLib.Finite.Defs` (programs), `AllxLib.Finite.Sieve`
(the bit-mask sieve decides primality), `AllxLib.Finite.Scan` (soundness of the scan), and
`scripts/gen_chunks.py` (generator of the chunk and chain files).
-/
import AllxLib.Finite.Chain

namespace Allx

/-- `Θ₁(t) = ∑_{p ≤ t} (−χ₄(p)) log p · (t − p) > 0` for every real `3 < t ≤ 10⁷`. -/
theorem theta1_wlog_pos_of_le {t : ℝ} (h3 : 3 < t) (hN : t ≤ 10000000) :
    0 < Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t :=
  theta1_pos_real 10000000 (fun _ h4 hk => theta1_cw_nat_pos h4 hk) h3 (by exact_mod_cast hN)

/-- `Θ₁(t) = ∑_{p ≤ t} (−χ₄(p)) log p · (t − p) ≥ 0` for every real `2 ≤ t ≤ 10⁷`. -/
theorem theta1_wlog_nonneg_of_le {t : ℝ} (h2 : 2 ≤ t) (hN : t ≤ 10000000) :
    0 ≤ Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t :=
  theta1_nonneg_real 10000000 (fun _ h4 hk => theta1_cw_nat_pos h4 hk) h2 (by exact_mod_cast hN)

end Allx
