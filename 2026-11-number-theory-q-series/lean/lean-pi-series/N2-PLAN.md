# N2 plan: Sun's `64√2/π²` series in Lean (work packages)

Target (Sun, arXiv:2307.03086, Conj. 5.33(i); Theorem 11 of `../ramanujan-series/RESULTS.md`):

    ∑_{k≥0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^{22k}) = 64√2/π².

All definitions are in `PiSeries/SunPi2/Defs.lean` (namespace `PiSeries.SunPi2`): `σ ρ : ℕ → ℕ → ℝ` (explicit
rational functions), `A n = 9/256 · ∏_{1≤m<n} ρ m 0`, `a n k = A n · ∏_{j<k} σ n j`, `T N = ∑_{k<N} a N k`,
`t` (Sun's summand), `xstar = √3/2`, `cf N k` (factorial closed form), lemmas `A_one`, `A_succ`,
`a_zero_right`, `a_succ_right`. **Do not change `Defs.lean`.**

## Mathematical structure (reorganised to be formalizable; differs from RESULTS.md Part B)

Part A (exact). For `n ≥ 1`, `k < n`: `a n k > 0`, `a (n+1) k = a n k * ρ n k`, and `T (n+1) = T n + t n / 256`,
`T 1 = t 0 / 256`; hence `∑_{n<N} t n = 256 T N`, `T` is increasing and bounded by `B := (∑' t)/256`, and
`a n k ≤ T n ≤ B`. For `k < N ≤ 2k+2`: `a N k = cf N k`.

Part B (`T N → √2/(4π²)`). Let `k0 = ⌊x* N⌋`, `x* = √3/2`.
* (T) Tails. `ρ̄(x) = (4x²−1)²/(64(1−x²)²)` is the limit of `ρ(n,k)`, `k = xn`; `ρ̄ < 1` for `x < x*`, `> 1` for
  `x > x*`. For `k ≤ (x*−δ)N`: `a(N,k) = a(n1,k) ∏_{n1≤n<N} ρ(n,k) ≤ B q^{N−n1}` with `n1 ≈ k/(x*−δ/2)`, `N−n1 ≥ cN`.
  For `k ≥ (x*+δ)N`: `a(N,k) = a(N',k)/∏_{N≤n<N'} ρ(n,k) ≤ B q1^{-(N'−N)}`, `N' ≈ k/(x*+δ/2)`. So everything outside
  `|k − x*N| ≤ δN` is `≤ N·B·q^N → 0`. (Uses only positivity, `a ≤ B`, and bounds on `ρ`. No log-concavity.)
* (R) Window. `σ(N,k) = Λ(k/N)(1+O(1/N))`, `Λ(x) = (1−x)²(2x+1)⁴/((1+x)²(2x−1)⁴)`, `Λ(x*) = 1`,
  `(log Λ)'(x*) = −24`. Hence for every `ε` there are `δ, C`: `|log σ(N,k) + 24(k−x*N)/N| ≤ ε|k−x*N|/N + C/N` for
  `|k−x*N| ≤ δN`; multiplying up, `a(N,k)/a(N,k0)` lies between `exp(−(12±ε)(k−k0)²/N ∓ C'/N)`.
* (G) `(1/√N) ∑_{|j|≤δN} exp(−c j²/N) → √(π/c)`.
* (C) Centre. `a(N,k0) √N = cf N k0 · √N → √3/(√2 π^{5/2})` by Stirling for the 11 factorials.
* Assembly: `T N = a(N,k0)√N · (1/√N)∑_{window} a(N,k)/a(N,k0) + tails`, squeeze with `ε → 0`:
  `T N → √3/(√2π^{5/2}) · √(π/12) = √2/(4π²)`. With Part A: `∑ t = 256 · √2/(4π²) = 64√2/π²`.

## Work packages (one file each, in `PiSeries/SunPi2/`, namespace `PiSeries.SunPi2`)

Rules for every package: Lean v4.33.0-rc2 + Mathlib (pinned); `import Mathlib` and (if stated) `import PiSeries.SunPi2.Defs`;
no `sorry`, no `axiom`, no `native_decide`; the theorem names and statements below must be proved **exactly as written**
(auxiliary lemmas are free); check with `export PATH="$HOME/.elan/bin:$PATH"; cd monotonicity/lean-pi-series;
lake env lean PiSeries/SunPi2/<File>.lean` (this compiles one file against the already-built dependencies;
do not run a bare `lake build`, other agents are compiling in the same project; never `lake clean`, never touch `.lake/`);
`open Finset Filter Topology` as needed; do not edit files of other packages; no git commit/push.
If a statement below is false or unprovable as written, say so precisely and prove the nearest correct statement.

### WP-A `PartA.lean` (imports Defs)
```lean
theorem a_pos {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : 0 < a n k
theorem a_succ_left {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : a (n + 1) k = a n k * ρ n k
theorem T_one : T 1 = t 0 / 256
theorem T_succ {N : ℕ} (hN : 1 ≤ N) : T (N + 1) = T N + t N / 256
theorem sum_t_eq {N : ℕ} (hN : 1 ≤ N) : ∑ n ∈ range N, t n = 256 * T N
theorem t_pos (n : ℕ) : 0 < t n
theorem summable_t : Summable t
theorem a_le_T {n k : ℕ} (hn : 1 ≤ n) (hk : k < n) : a n k ≤ T n
theorem T_le (N : ℕ) (hN : 1 ≤ N) : T N ≤ (∑' n, t n) / 256
theorem a_eq_cf {N k : ℕ} (hk : k < N) (hN : N ≤ 2 * k + 2) : a N k = cf N k
theorem hasSum_t_of_tendsto_T {L : ℝ} (h : Tendsto T atTop (𝓝 L)) : HasSum t (256 * L)
```
Sketch. Positivity: `σ n j > 0` for `j + 1 < n` (with `u = 8j−4n+12 ∈ 4ℤ`, the four odd numbers are `u∓3, u∓1`, product
`(u²−9)(u²−1) > 0` since `u = 0` or `u² ≥ 16`); `ρ m 0 > 0` for `m ≥ 1`. Column relation by induction on `k` from the
compatibility `ρ n k * σ (n+1) k = σ n k * ρ n (k+1)` (valid for `k + 1 < n`, i.e. all denominators non-zero).
`T_succ`: with the certificate `R n k = NUM n k / (131072 n (2n−1)(2n+1)(4n−1)(k−n)²(2k+1))` (`NUM`: see
`../ramanujan-series/lean/RamanujanSeries/Sun1pi2.lean`, which already proves the two polynomial identities
`wz_poly`, `wz_boundary_poly` for real `n k` — copy them), for `k ≤ n−2`: `a(n+1,k) − a(n,k) = g(n,k+1) − g(n,k)`,
`g = R a`; the row `k = n−1`: `a(n+1,n−1) − a(n,n−1) + a(n+1,n) + g(n,n−1) = 0`; so
`T(n+1) − T(n) = −g(n,0) = −R(n,0) A(n) = t(n)/256`. The last equality needs `A n` in closed form:
`A n = −t n / (256 R n 0)`; prove it by induction with `ρ n 0` and the recurrences of `C(2n,n)`, `C(4n,2n)`
(cf. `PiSeries/Wallis.lean: choose_two_mul_succ`, `PiSeries/SunSqrt2.lean: choose_six_step` for the technique;
the case `n = 1` of `T_succ` (no `k ≤ n−2`) must be handled). `summable_t`: ratio `t(n+1)/t(n) ≤ 1/2` (true limit 1/64) or
compare with a geometric series. `a_eq_cf`: first the diagonal `a (k+1) k = cf (k+1) k` by induction on `k`
(`a(k+2,k+1) = a(k+1,k) ρ(k+1,k) σ(k+2,k)`), then induction on `N` from `k+1` up to `2k+2` with `a_succ_left` and
`cf (N+1) k = cf N k * ρ N k` (factorial recurrences; valid when `N + 1 ≤ 2k+2`). Checked numerically in
`derivation/n2_check.py`.

### WP-G `Gaussian.lean` (Mathlib only)
```lean
theorem tendsto_gaussian_sum (c δ : ℝ) (hc : 0 < c) (hδ : 0 < δ) :
    Tendsto (fun N : ℕ => (1 / Real.sqrt N) *
        ∑ i ∈ range (2 * ⌊δ * N⌋₊ + 1), Real.exp (-c * ((i : ℝ) - ⌊δ * N⌋₊) ^ 2 / N))
      atTop (𝓝 (Real.sqrt (Real.pi / c)))
```
Sketch. `f(s) = exp(−c s²/N)` is even and decreasing on `[0,∞)`, `∫_ℝ f = √(πN/c)` (`integral_gaussian`). With
`M = ⌊δN⌋`: `∑_{|j|≤M} f(j) = 1 + 2∑_{j=1}^M f(j)`, and `∫_0^{M+1} f − 1 ≤ ∑_{j=1}^M f(j) ≤ ∫_0^M f ≤ ½√(πN/c)`
(`AntitoneOn.sum_le_integral`, `AntitoneOn.integral_le_sum` or by hand); `∫_{M+1}^∞ f ≤ ∫_{δN}^∞ exp(−cδ s) ds`-type
bound `→ 0` after dividing by `√N` (e.g. `f(s) ≤ exp(−cδ s)` for `s ≥ δN`, integral `= exp(−cδ²N)/(cδ)`).

### WP-S `Stirling.lean` (Mathlib only), then WP-C `Centre.lean` (imports Defs, Stirling) — same agent
```lean
theorem tendsto_log_factorial_linear (c : ℝ) (hc : 0 < c) (m : ℕ → ℕ) (u : ℕ → ℝ) (B : ℝ)
    (hu : ∀ N, |u N| ≤ B) (hm : ∀ᶠ N in atTop, (m N : ℝ) = c * N + u N) :
    Tendsto (fun N : ℕ => Real.log ((m N).factorial : ℝ)
        - (c * N * Real.log N + c * N * (Real.log c - 1) + u N * Real.log N + u N * Real.log c
           + 1 / 2 * Real.log (2 * Real.pi * c * N)))
      atTop (𝓝 0)

theorem tendsto_cf_centre :
    Tendsto (fun N : ℕ => cf N ⌊xstar * N⌋₊ * Real.sqrt N) atTop
      (𝓝 (Real.sqrt 3 / (Real.sqrt 2 * Real.pi ^ (5 / 2 : ℝ))))
```
Sketch (S). Mathlib: `Stirling.tendsto_stirlingSeq_sqrt_pi : Tendsto stirlingSeq atTop (𝓝 √π)` with
`stirlingSeq n = n! / (√(2n) (n/e)^n)`. So `log m! − (m log m − m + ½ log(2π m)) → 0` along `m = m N → ∞`; then
`m log m = (cN+u)(log N + log c + log(1+u/(cN)))` and `(cN+u) log(1+u/(cN)) − u → 0`, `log(2πm) − log(2πcN) → 0`.
Sketch (C). `k = k0 = x*N − θ`, `θ ∈ [0,1)`. The eleven factorials `m!` of `cf N k0` have `m = c N + u` with
(numerator) `(c,u) = (4,0), (8x*+4, −8θ), (4x*−2, 4−4θ)`; (denominator) `(2,0), (2,−2), (4x*+2, −4θ), (8x*−4, 8−8θ),
(1−x*, θ−1)` twice, `(1+x*, −θ)` twice. Signed sums: `∑±c = 0`, `∑±u = 0`. The `N`-linear exponent is
`N(∑±c log c − 22 log 2) = N(x*·I2 + I1 − 22 log 2) = 0` with
`I2 = 8 log(8x+4) + 4 log(4x−2) − 4 log(4x+2) − 8 log(8x−4) + 2 log(1−x) − 2 log(1+x) = 0` at `x = x*`
(equivalent to `(2x+1)²(1−x) = (2x−1)²(1+x)`, i.e. `8x³ = 6x`), and
`I1 = 8 log 2 + 4 log(8x+4) − 2 log(4x−2) − 4 log 2 − 2 log(4x+2) + 4 log(8x−4) − 2 log(1−x) − 2 log(1+x) = 22 log 2`
(uses `4x²−1 = 2`, `1−x² = 1/4`). The `θ`-dependent part of `∑±u log c` is `−θ·I2 = 0`; the rest is
`4 log(4x−2) + 2 log 2 − 8 log(8x−4) + 2 log(1−x) = log((1−x)²/(1024(2x−1)⁴)) = log(1/16384)`. The square-root factors give
`(2πN)^{-5/2} · 1/(1−x²) = 4 (2πN)^{-5/2}`. Prefactor `2^{12}(4N−1)(2k0+1) ~ 2^{15} x* N²`. Result:
`cf N k0 · √N → 2^{15} x* · 4 · (2π)^{-5/2}/16384 = √3/(√2 π^{5/2})`.
Useful: `(2−√3)`, `(√3−1)` positivity, `√3² = 3`; `⌊x*N⌋₊ ≤ x*N < ⌊x*N⌋₊ + 1`; for `N ≥ 8` the natural subtractions in
`cf N k0` are genuine (`2N ≤ 4k0+4`, `k0+1 ≤ N`). Check numerically first (python, mpmath).

### WP-R `Ratio.lean` (imports Defs)
```lean
theorem log_sigma_window (ε : ℝ) (hε : 0 < ε) :
    ∃ δ > 0, ∃ C : ℝ, ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ, |(k : ℝ) - xstar * N| ≤ δ * N →
      0 < σ N k ∧
      |Real.log (σ N k) + 24 * ((k : ℝ) - xstar * N) / N| ≤ ε * |(k : ℝ) - xstar * N| / N + C / N

theorem prod_sigma_window (ε : ℝ) (hε : 0 < ε) :
    ∃ δ > 0, ∃ C : ℝ, ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ, |(k : ℝ) - ⌊xstar * N⌋₊| ≤ δ * N →
      (⌊xstar * N⌋₊ ≤ k →
        Real.exp (-(12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N) ≤ ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i ∧
        ∏ i ∈ Ico ⌊xstar * N⌋₊ k, σ N i ≤ Real.exp (-(12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N)) ∧
      (k ≤ ⌊xstar * N⌋₊ →
        Real.exp ((12 - ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N - C / N) ≤ ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i ∧
        ∏ i ∈ Ico k ⌊xstar * N⌋₊, σ N i ≤ Real.exp ((12 + ε) * ((k : ℝ) - ⌊xstar * N⌋₊) ^ 2 / N + C / N))
```
Sketch. `σ N k = Λ(k/N) · E`, where `E` is a product of twelve factors `1 + O(1/N)`, uniformly for `k/N ∈ [0.8, 0.93]`
(there `N − k ≥ 0.07N`, `4k − 2N ≥ 1.2N`): `(2k+3)/(2k+1)`, `((N−k−1)/(N−k))²`, `((N+k)/(N+k+1))²`,
`(8k+4N+1+2i)/(8k+4N)`, `(8k−4N)/(8k−4N+9+2i)` (`i<4`). `g(x) = log Λ(x) = 2 log(1−x) + 4 log(2x+1) − 2 log(1+x) − 4 log(2x−1)`,
`g(x*) = 0`, `HasDerivAt g (−24) x*` (`g'(x) = −4/(1−x²) − 16/(4x²−1)`), so
`|g(x) + 24(x−x*)| ≤ ε|x−x*|` near `x*`. Second statement: sum the first over `i` between `k0` and `k`
(`∑_{i<j}(i − θ) = j²/2 + O(j)`, absorb `C·j/N ≤ ε j²/N + C²/(4εN)`); the `ε, C, δ` of the two statements need not agree.

### WP-T `Tails.lean` (imports Defs)
```lean
theorem rho_le_below (η : ℝ) (hη : 0 < η) (hη' : η < 1 / 10) :
    ∃ q : ℝ, 0 < q ∧ q < 1 ∧ ∃ n0 : ℕ, ∀ n ≥ n0, ∀ k : ℕ, (k : ℝ) ≤ (xstar - η) * n → ρ n k ≤ q
theorem rho_ge_above (η : ℝ) (hη : 0 < η) (hη' : η < 1 / 10) :
    ∃ q : ℝ, 1 < q ∧ ∃ n0 : ℕ, ∀ n ≥ n0, ∀ k : ℕ, (xstar + η) * n ≤ k → k < n → q ≤ ρ n k
theorem tail_bound (f : ℕ → ℕ → ℝ) (B : ℝ)
    (hpos : ∀ n k : ℕ, 1 ≤ n → k < n → 0 < f n k) (hB : ∀ n k : ℕ, 1 ≤ n → k < n → f n k ≤ B)
    (hρ : ∀ n k : ℕ, 1 ≤ n → k < n → f (n + 1) k = f n k * ρ n k)
    (δ : ℝ) (hδ : 0 < δ) (hδ' : δ < 1 / 10) :
    ∃ q : ℝ, 0 < q ∧ q < 1 ∧ ∃ N0 : ℕ, ∀ N ≥ N0, ∀ k : ℕ, k < N → δ * N ≤ |(k : ℝ) - xstar * N| →
      f N k ≤ B * q ^ N
```
Sketch. `ρ n k = P(n) · (8k+4n+1)(8k+4n+3)(8k−4n+5)(8k−4n+7) / (2048 (n−k)²(n+k+1)²)`,
`P(n) = (4n+1)(4n+3)²/(64 n(2n−1)(4n−1)) → 1/8`; note `(8k−4n+5)(8k−4n+7) > 0` always (consecutive odd numbers).
Below (`k ≤ yn`, `y = x*−η`): if `2k ≥ n`, every factor is monotone in `k`, so `ρ ≤` its value at `k = yn`, which tends to
`ρ̄(y) = (4y²−1)²/(64(1−y²)²) < 1`; if `2k < n`, crude bounds give `ρ ≤ 1/2` for large `n`.
Above (`yn ≤ k < n`, `y = x*+η`): `ρ ≥` value at `k = yn` `→ ρ̄(y) > 1`. (`ρ̄(x) < 1 ⟺ 12x² < 9`.)
`tail_bound`: see "Tails" above; `f N k = f n1 k · ∏_{n ∈ Ico n1 N} ρ n k` by induction.

### Assembly `Main.lean` (main session)
`T N → √2/(4π²)` from the five packages, then `HasSum t (64√2/π²)`.
