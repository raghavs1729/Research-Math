# Certified negativity of δ′(M) on the compact middle range

This note states exactly which inequalities `certify_derivative.py` evaluates,
with a proof or a precise citation for each, so that the code can be checked
against it line by line.  Lemma labels (`lem:…`, `prop:…`, `eq:…`) refer to
`paper/weighted-races.tex`.  Section 9 lists everything that is *not*
certified.

## 0. Statement

`χ` is the real primitive character mod `q ∈ {3, 4}` (both odd), `M = m + 1/2 > 0`,
`γ` runs over the positive ordinates of `L(s, χ)` (with multiplicity),

    c_γ(M) = 4M / sqrt(M² + γ²)          (= 2 a_γ(m) in the paper),
    c_γ'(M) = dc_γ/dM = 4γ² / (M² + γ²)^{3/2} > 0,
    Φ_M(t) = Π_γ J_0(c_γ(M) t),
    δ(M)   = 1/2 + (1/π) ∫_0^∞ (sin t / t) Φ_M(t) dt          (lem:invert(b), eq:gilpelaez).

**Claim certified per cell.**  For a cell `[M_a, M_b]` the program outputs
numbers `lo ≤ hi` and the word `PASS` only if it has verified in ball
arithmetic that `δ` is differentiable on `[M_a, M_b]` and

    lo ≤ δ′(M) ≤ hi < 0     for every real M ∈ [M_a, M_b].

Adjacent cells share their endpoint exactly (endpoints are exact rationals),
so a contiguous run of `PASS` cells proves that `δ` is strictly decreasing on
the union.

## 1. Hypotheses

* **(H1) GRH(χ)**, so that the ordinates are real.  **LI(χ)** is needed only
  to identify `δ(M)` with the logarithmic density of the weighted race
  (thm:interp(ii)); the computation is about the function `δ(M)` displayed
  above.
* **(H2) The variance identity and its derivative** (prop:varid, eq:varidq,
  lem:Dprop(i); theorems of the paper under GRH and `L(1/2,χ) ≠ 0`):

      σ²(M) := Σ_γ c_γ²/2 = Σ_γ 8M²/(M²+γ²)       = 4M (log Λ)′(M + 1/2),
      D(M)  := Σ_γ c_γ c_γ′ = Σ_γ 16Mγ²/(M²+γ²)²  = 4 (log Λ)′(M+1/2) + 4M (log Λ)″(M+1/2),

  `Λ(s) = (q/π)^{(s+1)/2} Γ((s+1)/2) L(s,χ)`.
* **(H3) The zero lists.**  `chi4_zeros.txt` (511 ordinates) and
  `chi3_zeros.txt` (537): each printed decimal is within `10⁻²⁰` of an
  ordinate, the zeros are simple, and there is no other zero with
  `0 < γ ≤ γ_N` (`certify_arb.py --K 0`, sign changes of Hardy's `Z` plus a
  winding-number count).  The prover loads each ordinate as the ball
  `decimal ± 10⁻²⁰` and checks that the balls are strictly ordered.
* **(H4)** Correctness of FLINT/Arb (python-flint 0.9.0), including
  `arb.bessel_j`, `digamma`, Hurwitz `zeta`, `sinc`, `exp`, `log`, `sqrt`.

## 2. Head and tail; the exact tail quantities

Head `H` = the first `N` listed ordinates (`N` = the whole list in all runs:
511 resp. 537).  Tail = every other zero; by (H3) every tail zero has
`γ > γ_N`.  Write `Φ = P·Ψ`, `P = Π_H`, `Ψ = Π_tail`.

**(2.1) Exact tail sums.**  From (H2),

    V_T(M) := Σ_tail c_γ²/2  = σ²(M) − Σ_H 8M²/(M²+γ²),
    D_T(M) := Σ_tail c_γ c_γ′ = D(M)  − Σ_H 16Mγ²/(M²+γ²)².

(`tail_VD`, `head_sums`.)  The program requires both to be certified `> 0` at
the two endpoints of the cell.

**(2.2) Monotonicity in M.**  Each term `8M²/(M²+γ²)` is increasing in `M`, so
`V_T(M_a) ≤ V_T(M) ≤ V_T(M_b)` on the cell.  For `D_T`:
`∂_M [16Mγ²/(M²+γ²)²] = 16γ²(γ²−3M²)/(M²+γ²)³ > 0` when `γ² > 3M²`; the
program checks `γ_N² > 3M_b²` (true up to `M_b = 369` resp. `404`), so
`D_T(M_a) ≤ D_T(M) ≤ D_T(M_b)`.  Hence the balls
`V := hull(V_T(M_a), V_T(M_b))`, `DT := hull(D_T(M_a), D_T(M_b))` contain
`V_T(M)`, `D_T(M)` for every `M` in the cell.  (This avoids the cancellation
that evaluating `σ²(M) − Σ_H` on a wide `M`-ball would cause.)

**(2.3) Closed forms.**  For real `s ≥ 3` (`loglam_derivs_series`; only
`s = M + 1/2 ≥ 4.1` occurs in the δ′ runs):

    (log Λ)′(s) = ½ log(q/π) + ½ ψ((s+1)/2) + (L′/L)(s),
    (log Λ)″(s) = ¼ ψ′((s+1)/2) + (L′/L)′(s),          ψ′(x) = ζ(2, x),
    (L′/L)(s)  = − Σ_n χ(n) Λ(n) n^{−s},     (L′/L)′(s) = Σ_n χ(n) Λ(n) log n · n^{−s}.

The series are summed over prime powers `n ≤ N₀ = 20000` in ball arithmetic;
the omitted part is bounded by the integral test (`log x · x^{−s}` and
`log² x · x^{−s}` are decreasing on `[N₀, ∞)`):

    |Σ_{n>N₀} χ(n)Λ(n) n^{−s}|       ≤ N₀^{1−s} ( L/a + 1/a² ),
    |Σ_{n>N₀} χ(n)Λ(n) log n n^{−s}| ≤ N₀^{1−s} ( L²/a + 2L/a² + 2/a³ ),      L = log N₀, a = s − 1,

and added as a radius.  `selftest` checks these balls against an independent
mpmath evaluation and against the complex-step evaluation used by
`race_density_arb.py`.

## 3. Bessel inequalities used

All from the product `J_0(z) = Π_k (1 − z²/j_{0,k}²)` (lem:phi(i)):
for `|z| < j_{0,1} = 2.4048…`,

    −log J_0(z) = Σ_{n≥1} (λ_n/n) z^{2n},     φ(z) := J_1(z)/J_0(z) = Σ_{n≥1} 2λ_n z^{2n−1},
    λ_n = Σ_k j_{0,k}^{−2n} > 0,   λ_1 = 1/4.

**(B1)** `h(z) := (−log J_0(z) − z²/4)/z⁴ = Σ_{n≥2}(λ_n/n) z^{2n−4}` is
nondecreasing on `[0, j_{0,1})`.  Hence for `0 ≤ z ≤ z₀ < j_{0,1}`:

    0 ≤ −log J_0(z) − z²/4 ≤ C₄ z⁴,     C₄ := h(z₀).

(With `z₀ = 1` this is eq:quartic, `C₄ = c₄ = 0.01762…`.)

**(B2)** `(φ(z) − z/2)/z³ = Σ_{n≥2} 2λ_n z^{2n−4}` is nondecreasing, hence
for `0 ≤ z ≤ z₀`:

    z/2 ≤ φ(z) ≤ z/2 + κ z³,     κ := (φ(z₀) − z₀/2)/z₀³      (κ = 0.07508… at z₀ = 1).

The program uses `z₀ = 1` and computes `C₄`, `κ` in Arb from `J_0(1)`, `J_1(1)`.

**(B3)** `|J_n(x)| ≤ 1` for real `x`, `n ≥ 0` (Bessel's integral
`J_n(x) = (1/π)∫_0^π cos(nτ − x sin τ) dτ`); `J_0′ = −J_1`,
`J_n′ = (J_{n−1} − J_{n+1})/2`, `J_1′(x) = J_0(x) − J_1(x)/x`.  Consequently
`|J_0″| = |J_2 − J_0|/2 ≤ 1` and `|J_1″| = |J_3 − 3J_1|/4 ≤ 1`.

**(B4)** `|J_1(z)| ≤ z/2` for `z ≥ 0` (lem:J0env(i)).

**(B5)** `|J_0(x)| ≤ x^{−1/2}` for `x > 0` (eq:watson; Watson §13.74).

## 4. Differentiation under the integral sign

lem:diffint(c) gives `δ′(M) = (1/π)∫_0^∞ ∂_MΦ_M(t) (sin t/t) dt`, but its
hypothesis "at least 8 ordinates `≤ 2M_a`" **fails on most of the priority
range**: `γ_8(χ₄) = 25.73`, `γ_8(χ₃) = 28.22`, so the lemma as stated applies
only for `M_a ≥ 12.87` resp. `14.11` (for `χ₄` a single ordinate lies below
`2·3.6`).  That hypothesis is used in the paper only to produce an integrable
dominant for `t ≥ 16`; we replace it by the dominant of Section 6, which the
program verifies for every cell.  The argument:

1. *Pointwise.*  For every `t ≥ 0` and every `M > 0`, `M ↦ Φ_M(t)` is
   differentiable with (eq:dPhi)

       ∂_MΦ_M(t) = − Σ_γ c_γ′ t · J_1(c_γ t) · Π_{γ′≠γ} J_0(c_{γ′} t),

   absolutely convergent.  This is lem:diffint(a); its proof does not use the
   8-ordinate hypothesis (it only uses `a_γ(M_b) → 0` as `γ → ∞`).
2. *Dominant.*  By (B4) and `|J_0| ≤ 1`: `|∂_MΦ_M(t)| ≤ Σ_γ c_γ′t · c_γt/2 = (t²/2) D(M) ≤ (t²/2) D_max`
   for all `t ≥ 0`; and for `t ≥ T₁` the bound (6.1) below,
   `|∂_MΦ_M(t)| ≤ (t²/2) D_max · A · t^{−(K−1)/2}` with `K ≥ 8` (enforced by the
   program).  Times `|sin t|/t ≤ min(1, 1/t)` this is integrable on `[0, ∞)`,
   uniformly for `M ∈ [M_a, M_b]`.
3. *Inversion.*  By (B5), `|Φ_M(t)| ≤ Π_{γ∈S} (c_γ(M_a) t)^{−1/2}` for
   `t ≥ T₁`, so `Φ_M ∈ L¹` and lem:invert(b) applies for each `M` in the cell.

The standard theorem (pointwise derivative + integrable uniform dominant) then
gives differentiability of `δ` on `[M_a, M_b]` and the formula for `δ′`.

## 5. Zone 1: `0 ≤ t ≤ T₁`

`T₁` is a rational chosen by the program subject to the **verified** condition

    c_{γ_N}(M_b) · T₁ ≤ z₀ = 1.                                                    (5.0)

Since `c_γ(M)` is increasing in `M` and decreasing in `γ`, (5.0) gives
`c_γ(M) t ≤ z₀` for every tail zero, every `M` in the cell and `t ≤ T₁`.
Let `c_Γ := c_{γ_N}(M_b)` (upper bound taken).

**(5.1) Tail product (eq:psi-enclose).**  By (B1) summed over the tail,

    Ψ = exp(−V_T t²/2 − ε),   0 ≤ ε ≤ C₄ t⁴ Σ_tail c_γ⁴ ≤ C₄ t⁴ · c_Γ² · 2V_T(M_b),

so `G·e^{−C₄S₄t⁴} ≤ Ψ ≤ G`, with `G = exp(−V_T(M) t²/2)` and
`S₄ := 2 c_Γ² V_T(M_b)`.

**(5.2) Tail log-derivative.**  `Ψ_M(t) = exp(Σ_tail log J_0(c_γ t))`; the
termwise `M`-derivatives `−φ(c_γ t) c_γ′ t` are bounded by
`(φ(z₀)/z₀) c_γ c_γ′ t²`, and `Σ_tail c_γ c_γ′ = D_T(M) ≤ D_T(M_b)` uniformly
near the cell, so (exactly as in the proof of lem:diffint(a), with `0.517`
replaced by `φ(1) < 0.5751`) `∂_MΨ = −Ψ·B_T`, `B_T(t) := Σ_tail φ(c_γ t) c_γ′ t`.
By (B2), using `Σ_tail c_γ³ c_γ′ ≤ c_Γ² D_T`:

    (t²/2) D_T ≤ B_T(t) ≤ (t²/2) D_T · (1 + 2κ c_Γ² t²).

**(5.3) The integrand.**  With `dP := ∂_M P = −t Σ_{γ∈H} c_γ′ J_1(c_γ t) Π_{γ′∈H∖γ} J_0(c_{γ′} t)`,

    F(t, M) := (sin t/t) ∂_MΦ_M(t) = (sin t/t) · Ψ · (dP − P·B_T).

`Cell.integrand` evaluates, for a box `[t_a, t_b] × [M_a, M_b]`, with `t` and
`M` the corresponding balls,

    F ∈ sinc(t) · exp(−V t²/2) · hull(e^{−C₄S₄t⁴}, 1) · ( dP − P · (t²/2)·DT·hull(1, 1 + 2κc_Γ²t²) ),

which by (2.2), (5.1), (5.2) contains `F(t, M)` for every point of the box.
`P` and `dP` are formed with prefix/suffix products (no division by `J_0`).
The coefficient balls are computed from expressions in which `M` occurs once,
`c_γ = 4/sqrt(1+(γ/M)²)`, `c_γ′ = (4/γ)(1+(M/γ)²)^{−3/2}`.

**(5.4) Bessel values on wide balls (`bessel01`).**  For the ball
`x = c_γ t` with midpoint `x₀ > 0` and radius `r`, Taylor's theorem with (B3):

    J_0(x) ∈ J_0(x₀) ± ( |J_1(x₀)| r + r²/2 ),
    J_1(x) ∈ J_1(x₀) ± ( |J_0(x₀) − J_1(x₀)/x₀| r + r²/2 ),

intersected with `[−1, 1]` (B3) and, for `J_0`, with `±x_lo^{−1/2}` when
`x_lo > 1` (B5).  `J_0(x₀)`, `J_1(x₀)` are Arb evaluations at an exact point.

**(5.5) Quadrature = interval Riemann sum.**  For any partition
`0 = t_0 < … < t_n = T₁` (exact dyadic rationals; the program verifies that
the final leaves tile `[0, T₁]`),

    ∫_0^{T₁} F(t, M) dt ∈ Σ_j (t_{j+1} − t_j) · F([t_j, t_{j+1}] × [M_a, M_b])      for every M in the cell,

because an integral over an interval lies in (length) × (range of the
integrand).  No quadrature rule, extrapolation or sampled derivative is
involved.  The partition is refined adaptively (largest radius first) until
the radius is below `tol × |midpoint|` or refinement stops helping; this
steering uses floats and affects only the width of the result.

## 6. Zone 2: `t ≥ T₁`

Let `S` = the first `K` head zeros, where `K` is the largest index with
`c_γ(M_a) T₁ ≥ 1` certified for all of them (the program requires `K ≥ 8`).
For `M` in the cell and `t ≥ T₁`, `c_γ(M) t ≥ c_γ(M_a) T₁ ≥ 1` for `γ ∈ S`.
In eq:dPhi bound `|J_1(c_γ t)| ≤ c_γ t/2` (B4), every factor outside
`S∖{γ}` by 1, and every factor in `S∖{γ}` by (B5):

    |∂_MΦ_M(t)| ≤ (t²/2) D(M) · max_γ Π_{γ′∈S∖{γ}} (c_{γ′}(M_a) t)^{−1/2}
               ≤ (t²/2) D_max · A · t^{−(K−1)/2},      A := c_{γ_1}(M_b)^{1/2} · Π_{γ′∈S} c_{γ′}(M_a)^{−1/2}.      (6.1)

(If `γ ∈ S`, removing its factor multiplies the full product over `S` by
`(c_γ(M_a) t)^{1/2} ≤ (c_{γ_1}(M_b) t)^{1/2}`; if `γ ∉ S` the full product is
at most itself times `(c_{γ_1}(M_b) t)^{1/2} ≥ 1`.)  Here
`D_max := upper(Σ_H c_γ c_γ′ on the M-ball) + D_T(M_b) ≥ D(M)`.  With
`|sin t|/t ≤ 1/t`:

    | ∫_{T₁}^∞ F(t, M) dt | ≤ (D_max/2) · A · ∫_{T₁}^∞ t^{1−(K−1)/2} dt
                           = (D_max/2) · A · T₁^{(5−K)/2} · 2/(K−5)  =: E₂.

(`Cell.remainder`.)  `T₁` is taken as small as possible subject to
`E₂ ≤ 10⁻¹⁰` (a steering target; whatever `E₂` results is added rigorously).

## 7. Assembly and output columns

    δ′(M) ∈ (1/π) · ( Σ_j (t_{j+1}−t_j) F_j + [−E₂, E₂] )     for every M ∈ [M_a, M_b].

Each line of `results_q*.txt` is

    M_a M_b hi=… lo=… E2=… Emodel=… rad=… T1=… K=… N=… tcells=… VT=[…] DT=[…] secs=… PASS

* `lo ≤ δ′(M) ≤ hi` on the whole cell; both printed decimals are rounded
  outward (`hi` up, `lo` down) from the Arb endpoints.  `PASS` iff `hi < 0`.
* **Error budget.**  `rad` = radius of the final enclosure = everything:
  `E2` (large-`t` remainder / π, Section 6) + `Emodel` (a rigorous upper bound
  for the part of the width caused by the tail model, i.e. by the two `hull`
  factors in (5.3), / π) + the rest, which is the interval Riemann-sum width
  caused by the `t`-subdivision and by the width of the `M`-cell.  In all runs
  `E2 < 2·10⁻¹⁰` and `Emodel < 10⁻⁷`, so the radius (about 2 % of `|δ′|`) is
  entirely subdivision/cell width.
* `VT`, `DT` = enclosures of `V_T`, `D_T` over the cell; `tcells` = number of
  `t`-boxes; `T1`, `K` as above.

A cell that is not certified is bisected (up to 6 times) and both halves are
retried; a line `… FAIL` is written if that does not help.  The `# SUMMARY`
block at the end of each file gives the union of certified cells and the
worst (closest to 0) certified upper bound.

## 8. Auxiliary certificates (`aux_results.txt`)

**(a) `Q_χ(M) = Σ_γ 4M³/(M²+γ²)^{3/2}`.**  Head sum over all listed zeros in
ball arithmetic.  For the unlisted zeros (`γ > γ_N`) put
`a_γ := 4M²/(M²+γ²)`, so that the term of `Q` is `a_γ^{3/2}/2`.  Two tail
moments are exact:

    A₁ := Σ_tail a_γ  = V_T(M)/2,
    A₂ := Σ_tail a_γ² = Σ_tail 16M⁴/(M²+γ²)² = 2V_T(M) − M·D_T(M)

(the second because `16Mγ²/(M²+γ²)² = 16M/(M²+γ²) − 16M³/(M²+γ²)²`, i.e.
`D_T = (2/M)V_T − A₂/M`), and `a_γ ≤ a_max := 4M²/(M²+γ_N²)`.  Hence

    A₂ / (2·sqrt(a_max))  ≤  Σ_tail a_γ^{3/2}/2  ≤  min( sqrt(a_max)·A₁ , sqrt(A₁A₂) ) / 2,

using `a² ≤ sqrt(a_max)·a^{3/2}`, `a^{3/2} ≤ sqrt(a_max)·a`, and Cauchy–Schwarz
`Σ a^{1/2}·a ≤ sqrt(Σa · Σa²)`.  (The bound `sqrt(a_max)·A₁/2 ≤ (M/γ_N)·V_T/2`
is the one suggested in the task; the Cauchy–Schwarz bound is about twice as
sharp and the lower bound is new.)  `Q` is strictly increasing in `M` (each
term `4(1+γ²/M²)^{−3/2}` is), so `Q(M) < 1` at one point proves `M < M_0(χ)`,
and `Q(M) > 1` proves `M > M_0(χ)`.

**(b1) `M*` with `σ²_{M*} = 2`.**  `σ²(M) = 4M (log Λ)′(M+1/2)` is strictly
increasing in `M` (`D > 0`).  For `s = M + 1/2 < 4` the derivative `L′(s)` is
taken by the complex step of `race_density_arb.sigma2_arb`
(`Im L(s+ih)/h`, `h = 2⁻⁶⁰`, remainder by Cauchy's estimate on `|z−s| ≤ 1/4`
with `|L(z)| ≤ |z|/Re z`), run at 384 bits.

**(b2) `Λ(M,M′) = 1/(MM′) − Σ_γ 4/(|M+iγ||M′+iγ|)`, `M < M′`.**  Head sum in
ball arithmetic.  With `W_M := Σ_{γ>T} 4/(M²+γ²) = V_T(M)/(2M²)` (exact),

    W_{M′} ≤ Σ_{γ>T} 4/(|M+iγ||M′+iγ|) ≤ sqrt(W_M W_{M′}),

the lower bound because `|M+iγ| ≤ |M′+iγ|`, the upper by Cauchy–Schwarz.

## 9. What is and is not certified — caveats

1. Everything is conditional on (H1)–(H4).  The zero-list certificate (H3) is
   the output of `prototype/explore/certify_arb.py --K 0`, a separate
   program; this directory uses its result and does not re-prove it.  (It was
   re-run in the same environment on 2026-11-24: all certifications passed,
   511/511 and 537/537 sign changes, winding numbers exactly 511 and 537,
   106.8 s.)
2. The variance identity and its derivative (H2) are theorems of the paper
   (under GRH), not machine-checked here; what is machine-checked is the
   evaluation of their right-hand sides.
3. The paper's lem:diffint(c) is used with a different (weaker, verified)
   decay hypothesis; the modified argument is Section 4 above and should be
   added to the paper if these results are cited.  Likewise (5.2) uses the
   tail log-derivative for arguments up to `z₀ = 1` rather than the paper's
   `1/2`; the proof is the same with `φ(1)` in place of `0.517`.
4. Floats are used (i) to steer the adaptive subdivision and the choice of
   `T₁` — every choice is re-verified in Arb before use — and (ii) to print.
   Printed `hi`/`lo` are outward-rounded decimals of the Arb endpoints.
5. The precision of the integrand is 64 bits; this is immaterial for
   validity (balls only get wider).
6. `crosscheck.py` and `selftest` are **not** part of the proof; they are
   independent float/mpmath sanity checks (Section 10).
7. Soundness, not sharpness, is certified: the enclosures are about 4 % wide
   relative to `|δ′|`; no claim is made about the exact value of `δ′`.
8. Stretch files (`results_q*_stretch.txt`) use the same method and code
   (coarser cells, ratio 1.01, tolerance 0.04); on those ranges monotonicity
   is *also* covered by the paper's prop:criterion-eval, whose endpoint→cell
   step and GRH-conditional count `N⁺` are listed as uncertified in
   rem:monotone-rigor.  The direct certificate here uses neither: the only
   cell-wide analytic inputs are the termwise monotonicity statements of (2.2)
   and the ball evaluation itself.  It does not touch the range `M > 320`
   resp. `M > 350`, i.e. lem:analytic-leg is still needed beyond.
9. The method is restricted to `M_b < γ_N/√3` (369 resp. 404) by (2.2) and to
   `s = M + 1/2 ≥ 3` by (2.3); both are checked at run time.
10. The value `M_0(χ)` of item (a) is enclosed only to the width allowed by
   the two-moment tail bound (about `2·10⁻⁴` resp. `3.5·10⁻⁴`).

## 10. Cross-checks performed (non-rigorous, `crosscheck.py`)

* **B.** A float64 Gil–Pelaez evaluation of `δ(M)` (numpy/scipy, Gaussian tail
  with the exact tail variance from mpmath) lies inside every certified
  density interval of `race_mono_grid_certified_results.txt` for
  `m ∈ [1.5, 100]`.
* **C.** At several `M` per modulus, the float `δ′(M)` from the derivative
  formula and a central finite difference (`h = 10⁻³`) of the float `δ` agree
  to `~10⁻⁹` and both lie inside the certified enclosure of the cell
  containing `M`.
* **D.** *Certified vs certified:* for consecutive grid points `m₁ < m₂` of
  the certified density grid inside the covered range, the interval
  `Σ_cells [lo, hi]·|cell ∩ [M₁, M₂]|` overlaps the certified interval for
  `δ(m₂) − δ(m₁)`.

Outputs: `crosscheck_q4.txt`, `crosscheck_q3.txt`.

## 11. Results (run of 2026-11-24, python-flint 0.9.0, 3 worker processes)

| file | range of M | cells | cell ratio / tol | worst certified upper bound for δ′ | enclosure width / abs(δ′) | max E2 | max Emodel | wall time (CPU) |
|---|---|---|---|---|---|---|---|---|
| `results_q4.txt` | [3.6, 17.4] | 317 | 1.005 / 0.02 | −1.675291e-3 on [17.3701, 17.4] | 4.0 % | 3.2e-11 | 4.9e-8 | 10.4 min (1806 s) |
| `results_q3.txt` | [4.4, 19.2] | 296 | 1.005 / 0.02 | −1.549494e-3 on [19.1282, 19.2] | 4.0 % | 3.2e-11 | 5.4e-8 | 12.7 min (2276 s) |
| `results_q4_stretch.txt` | [17.4, 320] | 293 | 1.01 / 0.04 | −1.224122e-5 on [317.8839, 320] | 8.0 % | 3.1e-11 | 1.4e-8 | 5.0 min (898 s) |
| `results_q3_stretch.txt` | [19.2, 350] | 292 | 1.01 / 0.04 | −1.098446e-5 on [347.3088, 350] | 8.0 % | 3.2e-11 | 1.4e-8 | 8.4 min (1514 s) |

Every cell passed at the first attempt (no bisection, no `FAIL` line), the
cells of each file are contiguous, and in every cell `hi/lo ≥ 0.9608`
(priority) resp. `≥ 0.9231` (stretch): the certified upper bound is never
closer to 0 than 92 % of the certified lower bound, so the sign is nowhere
near marginal.  "Worst" in the table is the absolute value closest to 0; it
sits at the top of each range simply because `|δ′|` decreases like
`M^{-3/2}`.  (The machine was shared; wall times are indicative.)

`aux_results.txt`: `3.688854 < M_0(χ₄) < 3.689059`,
`4.483812 < M_0(χ₃) < 4.484161`; `Q_{χ₄}(3.68) ≤ 0.994137`,
`Q_{χ₃}(4.47) ≤ 0.992173`; `1.82391477 ≤ M*(χ₄) < 1.82391478`,
`2.12775095 ≤ M*(χ₃) < 2.12775096`; `Λ(M,M′)` for the six pairs per
character (positive for `(0,1), (0,2), (0,3)` for both characters and for
`(1,2)` when `q = 3`; **negative** for `(1,2)` when `q = 4` and for
`(1,3), (2,3)` for both).

Cross-checks (`crosscheck_q4.txt`, `crosscheck_q3.txt`): B — the float
`δ(M)` lies inside all 20 certified density intervals per modulus; C — 11
points per modulus, float formula and finite difference agree to `≤ 2.4·10⁻⁹`
and all lie inside the certified cell enclosures; D — all 15 grid intervals
per modulus (`m` from 4 to 100) consistent.

## 12. Reproduction

    sh run_all.sh          # selftest, aux, both priority ranges, both stretch ranges, cross-checks
                           # (resumable: delete a results file to recompute it)

Single cell: `python3 certify_derivative.py cell --q 4 --Ma 3.6 --Mb 3.618 --tol 0.02`.
The computation is deterministic up to the order of lines in the result files
(cells are farmed out to 3 worker processes).
