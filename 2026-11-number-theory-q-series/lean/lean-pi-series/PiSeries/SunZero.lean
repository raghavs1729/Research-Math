/-
# Sun's zero-sum series (MathOverflow 506695)

  `∑_{k≥1} C(2k,k) 16^k (1656k⁴ − 2604k³ + 1158k² − 229k + 15)
            / (k² (2k−1)² (6k−1)(6k−5) C(6k,3k) C(3k,k)) = 0`.

Conjectured by Zhi-Wei Sun; the proof is Theorem 7 of `monotonicity/ramanujan-series/RESULTS.md`.
There is a WZ pair `(Φ, R Φ)` supported on the triangle `0 ≤ k < n`.  Here it is normalised as
`a(n,k) = 6π² Φ(n,k)`, which is rational: `a(1,0) = t(1) = −32/75` (`t` = Sun's summand) and

  `a(n,k+1)/a(n,k) = σ(n,k) = −4(4k+4n+7)(k+1)(n−k−1)(2k+2n+1)(2k+2n+2)
                               / ((4k+4n+3)(2k+2n+3)(2k+2n+4)(2k+4n+3)(2k+6n+3))`,
  `a(n+1,k)/a(n,k) = ρ(n,k) = (4k+4n+7)(4n+1)(4n+2)(4n+3)(4n+4)(2n+1)(2n+2)(2k+2n+1)(2k+2n+2)
       / (2(4k+4n+3)(n−k)(2k+2n+3)(2k+2n+4)(2k+4n+3)(2k+4n+5)(2k+6n+3)(2k+6n+5)(2k+6n+7))`,
  `R(n,k) = P(n,k) / (3(6n+1)(6n+5)(k−n)(2k+4n+3)(2k+6n+3)(2k+6n+5)(4k+4n+3))`

(`certificates/T7.gp`).  Then `−R(n,0) a(n,0) = t(n+1)` for `n ≥ 1`, so by `PiSeries.tri_sum`
`∑_{m=1}^{N} t(m) = T_N := ∑_{k<N} a(N,k)`; and `|σ| ≤ 1/24` on the triangle, `0 < ρ(n,0) ≤ 1/6`, so
`|T_N| ≤ (24/23)(32/75) 6^{1−N} → 0`.  Exact checks of all identities and of the certificates of the
inequalities: `derivation/t7_check.py`.
-/
import PiSeries.Triangle
import PiSeries.SunSqrt2

open Finset Filter Topology

namespace PiSeries.SunZero

/-! ### The WZ data -/

/-- numerator of the certificate -/
def PN (n k : ℝ) : ℝ :=
  (576 * n + 480) * k ^ 6 + (6912 * n ^ 2 + 10080 * n + 3600) * k ^ 5
    + (31104 * n ^ 3 + 66720 * n ^ 2 + 46576 * n + 10480) * k ^ 4
    + (59904 * n ^ 4 + 180000 * n ^ 3 + 195472 * n ^ 2 + 90464 * n + 14920) * k ^ 3
    + (25920 * n ^ 5 + 145152 * n ^ 4 + 263360 * n ^ 3 + 211992 * n ^ 2 + 77796 * n + 10454) * k ^ 2
    + (-62208 * n ^ 6 - 165120 * n ^ 5 - 137152 * n ^ 4 - 7880 * n ^ 3 + 44364 * n ^ 2 + 21138 * n
        + 2861) * k
    + (-52992 * n ^ 7 - 208128 * n ^ 6 - 334416 * n ^ 5 - 279712 * n ^ 4 - 127288 * n ^ 3
        - 29432 * n ^ 2 - 2609 * n + 12)

/-- denominator of the certificate -/
def RD (n k : ℝ) : ℝ :=
  3 * (6 * n + 1) * (6 * n + 5) * (k - n) * (2 * k + 4 * n + 3) * (2 * k + 6 * n + 3)
    * (2 * k + 6 * n + 5) * (4 * k + 4 * n + 3)

/-- `−σ = σM/σD` -/
def σM (n k : ℝ) : ℝ :=
  4 * (4 * k + 4 * n + 7) * (k + 1) * (n - k - 1) * (2 * k + 2 * n + 1) * (2 * k + 2 * n + 2)

def σD (n k : ℝ) : ℝ :=
  (4 * k + 4 * n + 3) * (2 * k + 2 * n + 3) * (2 * k + 2 * n + 4) * (2 * k + 4 * n + 3)
    * (2 * k + 6 * n + 3)

def ρN (n k : ℝ) : ℝ :=
  (4 * k + 4 * n + 7) * (4 * n + 1) * (4 * n + 2) * (4 * n + 3) * (4 * n + 4) * (2 * n + 1)
    * (2 * n + 2) * (2 * k + 2 * n + 1) * (2 * k + 2 * n + 2)

def ρD (n k : ℝ) : ℝ :=
  2 * (4 * k + 4 * n + 3) * (n - k) * (2 * k + 2 * n + 3) * (2 * k + 2 * n + 4)
    * (2 * k + 4 * n + 3) * (2 * k + 4 * n + 5) * (2 * k + 6 * n + 3) * (2 * k + 6 * n + 5)
    * (2 * k + 6 * n + 7)

/-- `a(n,k+1)/a(n,k)` -/
noncomputable def σ (n k : ℕ) : ℝ := -σM n k / σD n k

/-- `a(n+1,k)/a(n,k)` (for `k < n`) -/
noncomputable def ρ (n k : ℕ) : ℝ := ρN n k / ρD n k

/-- the WZ certificate -/
noncomputable def R (n k : ℕ) : ℝ := PN n k / RD n k

/-! ### The three polynomial identities -/

theorem compat_poly (n k : ℝ) :
    ρN n k * σM (n + 1) k * (σD n k * ρD n (k + 1))
      = σM n k * ρN n (k + 1) * (ρD n k * σD (n + 1) k) := by
  unfold ρN ρD σM σD
  ring

set_option maxRecDepth 100000 in
set_option maxHeartbeats 4000000 in
/-- the WZ equation, cleared of denominators -/
theorem wz_poly (n k : ℝ) :
    (ρN n k - ρD n k) * (RD n (k + 1) * σD n k * RD n k)
      = ρD n k * (PN n (k + 1) * (-σM n k) * RD n k - PN n k * (RD n (k + 1) * σD n k)) := by
  unfold ρN ρD RD σM σD PN
  ring

set_option maxRecDepth 100000 in
set_option maxHeartbeats 4000000 in
/-- the boundary row `k = n − 1` (`n = m + 1`, `k = m`), cleared of denominators -/
theorem wz_boundary_poly (m : ℝ) :
    (ρN (m + 1) m - ρD (m + 1) m) * (σD (m + 1 + 1) m * RD (m + 1) m)
      = -(ρN (m + 1) m * (-σM (m + 1 + 1) m) * RD (m + 1) m)
        - PN (m + 1) m * (ρD (m + 1) m * σD (m + 1 + 1) m) := by
  unfold ρN ρD RD σM σD PN
  ring

end PiSeries.SunZero
