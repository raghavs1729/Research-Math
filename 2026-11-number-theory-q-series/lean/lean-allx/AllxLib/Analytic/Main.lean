/-
AllxLib/Analytic/Main.lean — the analytic range: under GRH for `L(s, χ₋₄)`,
`Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p) > 0` for every `t ≥ 10⁷`.

Assembly of `theta1_wlog_ge` (prime-power bookkeeping, `PrimePowers.lean`), `psi1chi_le` (explicit
formula, `Psi1.lean`) and `final_ineq` (`Elementary.lean`).
-/
import AllxLib.Analytic.Psi1
import AllxLib.Analytic.Elementary

open Real

noncomputable section

namespace Allx

/-- **The analytic range.**  Under GRH for `L(s, χ₋₄)`:  `Θ₁(t) > 0` for every `t ≥ 10⁷`. -/
theorem theta1_wlog_pos_of_ge (hGRH : GRHFor chi4) {t : ℝ} (ht : 10000000 ≤ t) :
    0 < Coupling.Theta1 (fun p => Coupling.w4 p * Real.log p) t := by
  have h1 := theta1_wlog_ge (t := t) (by linarith)
  have h2 := psi1chi_le hGRH (t := t) (by linarith)
  have h3 := final_ineq ht
  linarith

end Allx
