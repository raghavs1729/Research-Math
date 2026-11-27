/-
RaceLib/PNT.lean — the prime number theorem `θ(x) ~ x` from the Zeta23 port of PNT+.

`MediumPNT` (Zeta23/FromPNTPlus/MediumPNT.lean, ported from the PrimeNumberTheoremAnd project of
Kontorovich, Tao et al.) states  ψ(x) − x = O(x exp(−c (log x)^{1/10})).  Together with Mathlib's
`Chebyshev.isBigO_psi_sub_theta_sqrt` (ψ − θ = O(√x)) this gives θ(x)/x → 1.
-/
import ChebLib.Race.Defs
import Zeta23.FromPNTPlus.MediumPNT

open Filter Topology Asymptotics

noncomputable section

namespace RaceLib

theorem tendsto_psi_div_self : Tendsto (fun x : ℝ => Chebyshev.psi x / x) atTop (𝓝 1) := by
  obtain ⟨c, hc, hO⟩ := MediumPNT
  have hexp : Tendsto (fun x : ℝ => Real.exp (-c * (Real.log x) ^ ((1 : ℝ) / 10))) atTop (𝓝 0) := by
    apply Real.tendsto_exp_atBot.comp
    have h1 : Tendsto (fun x : ℝ => (Real.log x) ^ ((1 : ℝ) / 10)) atTop atTop :=
      (tendsto_rpow_atTop (by norm_num)).comp Real.tendsto_log_atTop
    exact h1.const_mul_atTop_of_neg (show -c < 0 by linarith)
  have ho : (fun x : ℝ => x * Real.exp (-c * (Real.log x) ^ ((1 : ℝ) / 10)))
      =o[atTop] (fun x : ℝ => x) := by
    have h1 : (fun x : ℝ => Real.exp (-c * (Real.log x) ^ ((1 : ℝ) / 10)))
        =o[atTop] (fun _ : ℝ => (1 : ℝ)) := (isLittleO_one_iff ℝ).mpr hexp
    have h2 := (isBigO_refl (fun x : ℝ => x) atTop).mul_isLittleO h1
    simpa using h2
  have hψ : (Chebyshev.psi - id) =o[atTop] (fun x : ℝ => x) := hO.trans_isLittleO ho
  have h2 := hψ.tendsto_div_nhds_zero.add_const 1
  rw [zero_add] at h2
  refine h2.congr' ?_
  filter_upwards [eventually_gt_atTop 0] with x hx
  simp only [Pi.sub_apply, id]
  field_simp
  ring

/-- **The prime number theorem**, `θ(x)/x → 1`, proved from `MediumPNT`. -/
theorem pntTheta : PNTTheta := by
  unfold PNTTheta
  have hsqrt : (fun x : ℝ => Real.sqrt x) =o[atTop] (fun x : ℝ => x) := by
    rw [isLittleO_iff_tendsto' (by
      filter_upwards [eventually_gt_atTop 0] with x hx h0
      exact absurd h0 hx.ne')]
    have h : Tendsto (fun x : ℝ => (Real.sqrt x)⁻¹) atTop (𝓝 0) :=
      tendsto_inv_atTop_zero.comp Real.tendsto_sqrt_atTop
    refine h.congr fun x => ?_
    rw [Real.sqrt_div_self]
  have h2 : (Chebyshev.psi - Chebyshev.theta) =o[atTop] (fun x : ℝ => x) :=
    Chebyshev.isBigO_psi_sub_theta_sqrt.trans_isLittleO hsqrt
  have h4 := tendsto_psi_div_self.sub h2.tendsto_div_nhds_zero
  rw [sub_zero] at h4
  refine h4.congr fun x => ?_
  simp only [Pi.sub_apply]
  ring

end RaceLib
