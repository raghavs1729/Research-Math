/-
# Algebraic core of Theorem 11 (Sun's series for 64 √2 / π²)

  ∑_{k ≥ 0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^(22k)) = 64 √2 / π².

Formalized here: the two rational-function identities behind the WZ pair on the triangle `k < n`
(RESULTS.md, Theorem 11, Part A), for the explicit ratios `σ`, `ρ` and the certificate `R`:

* `wz_rat`      : `ρ n k - 1 = R n (k+1) * σ n k - R n k`           (WZ equation, `k ≤ n - 2`)
* `wz_boundary` : `ρ n (n-1) - 1 = -ρ n (n-1) * σ (n+1) (n-1) - R n (n-1)`   (the row `k = n - 1`)

as identities over `ℝ` under the non-vanishing of the denominators.  NOT formalized: the telescoping on
the triangle, the identification with the binomial summand, and the Laplace asymptotics (Part B).
-/
import Mathlib

namespace SunSeries1pi2

set_option maxRecDepth 100000
set_option maxHeartbeats 4000000

/-- numerator of the certificate -/
def NUM (n k : ℝ) : ℝ :=
  (3145728*n + 1310720)*k^6 + (-6291456*n^2 + 3670016*n + 2621440)*k^5
  + (-786432*n^3 - 11337728*n^2 - 344064*n + 1724416)*k^4
  + (9437184*n^4 + 786432*n^3 - 8716288*n^2 - 1884160*n + 376832)*k^3
  + (-8060928*n^5 + 10010624*n^4 + 3973120*n^3 - 3434496*n^2 - 892672*n - 15936)*k^2
  + (2752512*n^6 - 6324224*n^5 + 1916928*n^4 + 2363392*n^3 - 525824*n^2 - 170112*n - 12672)*k
  + (-344064*n^7 + 1134592*n^6 - 922624*n^5 - 195328*n^4 + 318016*n^3 - 1680*n^2 - 10044*n - 945)

/-- denominator of the certificate -/
def RD (n k : ℝ) : ℝ := 131072 * n * (2*n - 1) * (2*n + 1) * (4*n - 1) * (k - n)^2 * (2*k + 1)

/-- the WZ certificate -/
noncomputable def R (n k : ℝ) : ℝ := NUM n k / RD n k

/-- numerator and denominator of `σ(n,k) = a(n,k+1)/a(n,k)` -/
noncomputable def σN (n k : ℝ) : ℝ :=
  (2*k + 3) * (n - k - 1)^2 *
      ((4*k + 2*n + 1/2) * (4*k + 2*n + 3/2) * (4*k + 2*n + 5/2) * (4*k + 2*n + 7/2))
noncomputable def σD (n k : ℝ) : ℝ :=
  (2*k + 1) * (n + k + 1)^2 *
      ((4*k - 2*n + 9/2) * (4*k - 2*n + 11/2) * (4*k - 2*n + 13/2) * (4*k - 2*n + 15/2))
noncomputable def σ (n k : ℝ) : ℝ := σN n k / σD n k

/-- numerator and denominator of `ρ(n,k) = a(n+1,k)/a(n,k)` -/
noncomputable def ρN (n k : ℝ) : ℝ :=
  (n + 3/4)^2 * (n + 1/4) * (4*k + 2*n + 1/2) * (4*k + 2*n + 3/2) * (4*k - 2*n + 5/2) * (4*k - 2*n + 7/2)
noncomputable def ρD (n k : ℝ) : ℝ := 512 * n * (n - 1/4) * (2*n - 1) * (n - k)^2 * (n + k + 1)^2
noncomputable def ρ (n k : ℝ) : ℝ := ρN n k / ρD n k

/-- **The WZ equation of Theorem 11, cleared of denominators** (a polynomial identity in `n`, `k`):
`(ρN - ρD) · RD(n,k+1) · σD · RD(n,k) = ρD · (NUM(n,k+1) · σN · RD(n,k) - NUM(n,k) · RD(n,k+1) · σD)`. -/
theorem wz_poly (n k : ℝ) :
    (ρN n k - ρD n k) * (RD n (k + 1) * σD n k * RD n k)
      = ρD n k * (NUM n (k + 1) * σN n k * RD n k - NUM n k * (RD n (k + 1) * σD n k)) := by
  unfold ρN ρD RD σN σD NUM
  ring

/-- **WZ equation** `ρ - 1 = R(n,k+1) σ - R(n,k)` as a rational-function identity. -/
theorem wz_rat (n k : ℝ) (hρ : ρD n k ≠ 0) (hσ : σD n k ≠ 0) (hR : RD n k ≠ 0) (hR' : RD n (k + 1) ≠ 0) :
    ρ n k - 1 = R n (k + 1) * σ n k - R n k := by
  unfold ρ R σ
  rw [div_sub_one hρ, div_mul_div_comm, div_sub_div _ _ (mul_ne_zero hR' hσ) hR,
    div_eq_div_iff hρ (mul_ne_zero (mul_ne_zero hR' hσ) hR)]
  linear_combination wz_poly n k

/-- **Boundary row `k = n - 1`, cleared of denominators.** -/
theorem wz_boundary_poly (n : ℝ) :
    (ρN n (n - 1) - ρD n (n - 1)) * (σD (n + 1) (n - 1) * RD n (n - 1))
      = -(ρN n (n - 1) * σN (n + 1) (n - 1) * RD n (n - 1))
        - NUM n (n - 1) * (ρD n (n - 1) * σD (n + 1) (n - 1)) := by
  unfold ρN ρD RD σN σD NUM
  ring

/-- **Boundary row `k = n - 1`**: `ρ(n,n-1) - 1 = -ρ(n,n-1) σ(n+1,n-1) - R(n,n-1)`. -/
theorem wz_boundary (n : ℝ) (hρ : ρD n (n - 1) ≠ 0) (hσ : σD (n + 1) (n - 1) ≠ 0) (hR : RD n (n - 1) ≠ 0) :
    ρ n (n - 1) - 1 = -(ρ n (n - 1) * σ (n + 1) (n - 1)) - R n (n - 1) := by
  unfold ρ R σ
  rw [div_sub_one hρ, div_mul_div_comm, ← neg_div, div_sub_div _ _ (mul_ne_zero hρ hσ) hR,
    div_eq_div_iff hρ (mul_ne_zero (mul_ne_zero hρ hσ) hR)]
  linear_combination (ρD n (n - 1)) * wz_boundary_poly n

end SunSeries1pi2
