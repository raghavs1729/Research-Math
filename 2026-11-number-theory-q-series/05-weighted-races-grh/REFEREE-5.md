# Referee report 5: Theorem E, Proposition C₀, Corollary C′ (expanded proof), Lemma C″, and the envelope certificate

Date 2026-11-24. Scope: `BRIDGE-ATTEMPT.md` (Proposition C₀, expanded proof of Corollary C′,
Lemma C″, §3, §3A), and `envelope-cert/` (`ANALYTIC-LEG.md`, `analytic_leg.py`, `METHOD.md`,
`certify_envelope.py`, `results_q4.txt`, `results_q3.txt`, `SUMMARY.md`). I relied on
`REFEREE-3.md` for Lemma D, Lemma J, Lemma L and Theorem C's inequality, and re-tested them
numerically where they enter. Nothing under review was edited. My scripts and outputs are in
`monotonicity/referee5/` (list in §8). I made no commit; a concurrent process in this
repository swept some `referee5/` files into its own commits while I worked.

Reviewed versions (SHA-256, first 12 hex digits): `BRIDGE-ATTEMPT.md` 42d109ac0cd5,
`ANALYTIC-LEG.md` 0ea5cde930ac, `analytic_leg.py` 06d9b32cd2d6, `METHOD.md` 7e688a148a20,
`certify_envelope.py` 032caf29180a, `results_q4.txt` eefd98da6472, `results_q3.txt`
bc06e9e7ea82, `SUMMARY.md` a4b476f3c62d, shared module `cert/certify_derivative.py`
f4674b09e182. python-flint 0.9.0, FLINT 3.6.0.

## 0. Verdict

**No FATAL error. One REPAIRABLE item (R1), which does not touch any theorem. Thirteen MINOR items.**

* **Theorem E is proved as written.** For every real primitive `χ` mod `q ≥ 3` (both parities)
  with GRH(χ) and `L(1/2,χ) ≠ 0`, and every `M ≥ 40`:
  `−Mδ′(M) ≥ G(σ²(M))/(sqrt(2π)σ(M)) ≥ 0.3219/(sqrt(2π)σ(M)) > 0`. I re-derived every step
  (§1), reran the Arb script (output identical), and tested the conclusion numerically for 15
  characters and on 400 + 8 abstract amplitude sets (no violation; the inequality is
  asymptotically sharp in the abstract class, so there is no hidden slack to hide an error).
* **Proposition C₀, the expanded proof of Corollary C′, and Lemma C″ are correct** (§2).
* **The certificate `certify_envelope.py` is sound** (§3): it evaluates exactly the hypothesis
  of Corollary C′, case (b), `n = 12`, with Theorem C's own `NUM`; every rounding goes the right
  way; no float is on a certified path. A full fresh rerun reproduces both results files
  line for line (765 + 765 cells, all `PASS`, exact tiling of `[1, 2000]`). My independent
  mpmath evaluation of the exact cell quantity at 13 cells per modulus satisfies
  `LB ≤ E` with `LB/E ∈ [0.983, 0.997]`.
* **The chain for `q = 3, 4` is complete** (§4), and the general statement of §5 is supported.
* **R1 (REPAIRABLE, presentation of the proof of Theorem E):** `analytic_leg.py` prints
  numbers labelled `>=` / `<=` with round-to-nearest, so several displayed inequalities in
  `analytic_leg.out` and `ANALYTIC-LEG.md` are false in their last digit, including the
  displayed step (E5) "`G(236.9711) ≥ 0.321972`" (true lower endpoint `0.3219715`) and
  "`T_c ≤ 8.2662`" (`T_c = 8.266219…`). The certified decision uses the Arb balls, and the
  theorem's constant `0.3219` is correct. Fix and corrected table in §6.

## 1. Theorem E

### 1.1 Inputs, re-derived

* **Variance identity (V).** `ξ(s,χ) = (q/π)^{(s+a)/2}Γ((s+a)/2)L(s,χ)` is entire of order 1.
  Hadamard: `ξ′/ξ(s) = B(χ) + Σ_ρ(1/(s−ρ) + 1/ρ)`, `Re B(χ) = −Σ Re(1/ρ)`. For real `χ`, `B` is
  real and the zeros are symmetric under conjugation, so `ξ′/ξ(s) = Σ_ρ 1/(s−ρ)` (paired). Under
  GRH and `L(1/2,χ) ≠ 0`, with `s = M + 1/2`: `1/(M−iγ) + 1/(M+iγ) = 2M/(M²+γ²)`. Hence
  `h(M) := (log ξ)′(M+1/2) = Σ_{γ>0} 2M/(M²+γ²)` (with multiplicity) and
  `σ²(M) = Σ 8M²/(M²+γ²) = 4M·h(M)`,
  `h = ½log(q/π) + ½ψ((s+a)/2) + (L′/L)(s)`. This is the paper's `eq:varidq`, valid for both
  parities. Correct.
* **(V′).** `Σ(M²+γ²)^{−1} = h/(2M)`; termwise `d/dM` (locally uniformly convergent) gives
  `Σ(M²+γ²)^{−2} = h/(4M³) − h′/(4M²)`, so `S_4 = 256M⁴Σ(M²+γ²)^{−2} = 64Mh − 64M²h′ = 16σ² − 64M²h′`,
  with `h′ = ¼ψ′((s+a)/2) + (L′/L)′(s)`. Correct. Note `h′ = Σ 2(γ²−M²)/(M²+γ²)²`, so `h′ > 0`
  is not a triviality (the trivial bound from `c_j ≤ 4` is only `S_4 ≤ 32σ²`); it fails at
  small `M` for characters with a low zero.
* **Dirichlet-series bounds.** `|(L′/L)(s)| ≤ ΣΛ(n)n^{−s} ≤ Σ_{n≥2} log n·n^{−s}`,
  `|(L′/L)′(s)| ≤ ΣΛ(n)log n·n^{−s} ≤ Σ_{n≥2}log²n·n^{−s}`, for any character. First term plus
  `∫_2^∞`: `∫_2^∞ log x·x^{−s}dx = 2^{1−s}(L/(s−1) + 1/(s−1)²)`,
  `∫_2^∞ log²x·x^{−s}dx = 2^{1−s}(L²/(s−1) + 2L/(s−1)² + 2/(s−1)³)`; the integrands decrease on
  `[2,∞)` when `log 2 > 1/s` resp. `log 2 > 2/s`, i.e. for `s > 2.89`. `Z_1`, `Z_2` as stated.
  Correct.

### 1.2 The steps

| step | claim | check |
|---|---|---|
| 0 | `c_j ≤ 4` | `c_j = 4/sqrt(1+γ_j²/M²) < 4`. Trivial. |
| 1 | `h′ > 0`, hence `S_4 ≤ 16σ²`, all `χ`, `M ≥ 40`, both parities | `ψ′(x) = Σ_{k≥0}(x+k)^{−2} ≥ 1/x`; `ψ′` decreasing, so `¼ψ′((s+a)/2) ≥ ¼ψ′((s+1)/2) ≥ 1/(2(s+1))` for `a ∈ {0,1}` (odd parity is the worst case). `(s+1)Z_2(s) = (s+1)2^{−s}·[L² + 2(L²/(s−1) + …)]` is a product of positive decreasing functions for `s > 1/log 2 − 1`, so `1/(2(s+1)) > Z_2(s)` at `s = 40.5` gives it for all `s ≥ 40.5`. Margin at `40.5`: `0.0120481`. The bound does not depend on `q`. Proved. Numerically `h′(40) = 0.01234` (odd), `0.01266` (even) for every character tested. |
| 2 | `σ²(M) ≥ 236.9711` for all `χ`, `M ≥ 40` | `σ²` increasing termwise. `σ²(40) = 160h(40) ≥ 160[½log(q/π) + ½ψ((40.5+a)/2) − Z_1(40.5)]`, increasing in `q`. Real primitive characters exist for conductors 3 (odd), 4 (odd), 5 (even), 7 (odd), 8 (one odd, one even), 11, 12, 13, …: smallest odd conductor 3, smallest even conductor 5. Values: `(a,q) = (1,3)`: `236.971101087`; `(0,5)`: `275.837463`. Minimum is `χ_3`, and the bound is tight there (`σ²_{χ_3}(40) = 236.9711010873`). My scan of all 39 fundamental discriminants `|d| ≤ 60` confirms the ordering (`referee5/thmE_even.out`). Proved. |
| 3 | counting lemma `N_χ(Γ) ≥ (Y/2)[h(Y) − (2Y/Γ)h(Γ)]` | In `σ²(Y) = Σ 8Y²/(Y²+γ²)`: each zero `γ ≤ Γ` contributes `< 8`; for `γ > Γ`, `8Y²/(Y²+γ²) ≤ (2Y²/Γ²)·8Γ²/(Γ²+γ²)` is equivalent to `Γ² ≤ 2Y² + γ²`, true. The high zeros are then bounded by `(2Y²/Γ²)σ²(Γ)` (all terms positive, so adding the low zeros back only weakens it). With `σ² = 4Yh(Y)`, `4Γh(Γ)`: the lemma. Valid for all `Y, Γ > 0`; rigorous; counts with multiplicity (the sums are with multiplicity). |
| 3 | `N_χ(200) ≥ 11.375` for all `χ` | `h(Y) ≥ ½log(q/π) + ½ψ(·) − Z_1(Y+½)`, `h(Γ) ≤ ½log(q/π) + ½ψ(·) + Z_1(Γ+½)`; coefficient of `½log(q/π)` is `1 − w = 3/5 ≥ 0`, so the bound increases with `q`: uniform in `q`, worst case smallest conductor of each parity. Values: `11.37517` (`a=1, q=3`), `14.21017` (`a=0, q=5`). So `N_χ(200) ≥ 12` (an integer), in particular `≥ 9`. Proved. |
| 3 | nine ordinates with multiplicity suffice | Yes. Lemma D and Proposition C₀ are stated for a nondecreasing sequence; a repeated ordinate gives a repeated factor `J_0(c t)`, and the decay argument only counts factors. Nine distinct zeros are not needed. |
| 3 | `c_i(M) ≥ 4/sqrt(26)` for `i ≤ 9`, `M ≥ 40` | `γ_i ≤ 200 ≤ 5M`. Correct. |
| 4 | envelopes | `t ≤ τ`: `c_jt ≤ 4τ = z_* < j_{0,1}`, `0 < J_0(z) ≤ e^{−z²/4}`, so `|Φ| ≤ e^{−σ²t²/2}`, `|Φ^{(j)}| ≤ e^{−(σ²−8)t²/2}`. `t ≥ τ`, `j ∉ 𝓛`: `|J_0(c_jt)| ≤ b(c_jt) ≤ b(c_jτ) = e^{−c_j²τ²/4}`; product `= e^{−(τ²/2)(σ² − Σ_𝓛c_i²/2)} ≤ ρ^{(σ²−72)/8}` (`τ² = log(1/ρ)/4`; monotone in the exponent, no sign condition needed). `i ∈ 𝓛`: `|J_0(c_it)| ≤ min(1,(8/π)(c_κt)^{−1/2}) = p(t)`, `T_c = (8/π)²/c_κ`. Leave-one-out: `ρ^{(σ²−80)/8}p^8` in both cases `j ∈ 𝓛`, `j ∉ 𝓛`. Correct. |
| 5 | Proposition C₀ applies | `z_1 = z_* < 2.4 < j_{0,1}`; `t_1 = τ = 0.449 ≤ π/2`; `τ = z_*/4 ≤ z_*/c_1`; `Σ² = σ²`; `Q_4 = 16σ² ≥ S_4`; ≥ 9 ordinates (step 3), `Σγ^{−2} < ∞`. Correct. |
| 5 | first integral | `e^{−x}cos t ≥ 1 − x − t²/2` on `[0,τ]`, `x ≥ 0` (two cases, both checked). Gaussian moments and the Mills bound `∫_τ^∞e^{−σ²t²/2} ≤ e^{−σ²τ²/2}/(σ²τ)` are right; `16β_*σ²·3σ^{−4} = 48β_*/σ²`. |
| 5 | second, third integrals | `∫_0^∞p^9 = T_c + (2/7)T_c = (9/7)T_c`; `∫_0^∞t²p^8 = T_c³/3 + T_c³ = (4/3)T_c³`; `(σ²/2π)·sqrt(π/2)(σ²−8)^{−3/2} = (sqrt(2π)σ)^{−1}·½(1−8/σ²)^{−3/2}` (needs `σ² > 8`). Collecting gives exactly `G(σ²)`. |
| 6 | `G` increasing on `[15,∞)` | each subtracted term decreasing: `1/x`; `(1−8/x)^{−3/2}` (`x > 8`); `x^{−1/2}ρ^{x/8}`; `x^{1/2}ρ^{x/8}` (`x > 4/log(1/ρ) = 4.954`); `x^{3/2}ρ^{x/8}` (`x > 12/log(1/ρ) = 14.862`). Proved, not assumed. |

**Uniformity in `M` on `[40, ∞)`.** Obtained from four evaluations at one point by four
monotonicity statements, each proved above: `(s+1)Z_2(s)` decreasing; `σ²(M)` increasing;
`c_i(M)` increasing (so `γ_i ≤ 200` suffices for all `M ≥ 40`); `G` increasing on `[15,∞)`
with `σ_0² = 236.97 > 15`. Nothing is interpolated.

**What Proposition C₀ needs, exactly:** Lemma D (≥ 9 positive ordinates, `Σγ^{−2} < ∞`);
Lemma L on `[0, t_1]` with `t_1 ≤ min(π/2, z_1/sup_jc_j)`; `|J_1(z)| ≤ z/2`; `|sin t| ≤ t`;
`Σε_jc_j² = S_4/16`; numbers `Σ² ≥ σ²`, `Q_4 ≥ S_4`; measurable `W ≥ |Φ|` on `[t_1,∞)` and
`W¹ ≥ |Φ^{(j)}|` on `[0,∞)` for every `j`. It does not need Lemma J, any ordinate, or any
ordering of the zeros beyond the value of `sup_jc_j`.

**High-precision recomputation (mpmath, 30 digits, `referee5/` inline run):**
`z_* = 1.79714921691`, `τ = 0.449287304228`, `J_0(z_*) = 0.341644119`, `β_* = 0.0255529512`,
`T_c = 8.2662190805`, (E2) margin `0.0120481928`, `σ_0² = 236.971101087`, count `11.3751732`,
`G(σ_0²) = 0.3219715139`, `G(236.9711) = 0.3219714989`. The exact value of the three C₀ integrals
at `σ² = σ_0²`, times `sqrt(2π)σ`, is `0.3221796` (pieces `0.992916 − 0.0000073 − 0.670729`),
so the closed-form bounding of step 5 loses only `2·10⁻⁴`, in the right direction.

### 1.3 The script `analytic_leg.py`

Read line by line. `Z1`, `Z2`, `beta_star`, `G`, the two families, (E2)–(E5) and the side
conditions (`s ≥ 4`, `Y + ½ ≥ 4`, `2Y ≤ κM_E`, `σ_0² > 80`) are the formulas of the proof.
All decisions use lower endpoints of balls (`e2.lower() > 0`, `cnt.lower() >= 9`,
`g.lower() > 0`, `x0 = sig0.lower()` with `G` increasing). Rerun: output identical to
`analytic_leg.out`. The only defect is the printing (R1).

### 1.4 Numerical tests of the conclusion (non-rigorous; `referee5/thmE_numeric.py`, `thmE_even.py`)

`−Mδ′` computed two ways: Lemma D's formula `f(1) + R`, and a central difference of the
Gil–Pelaez integral; zeros above the list through the cumulant expansion of `log J_0` to order
`t⁶` with tail sums from my own closed forms. The two agree to 6 digits in all cases, which also
re-tests Lemma D for 15 characters. Column "bound/truth" must be `≤ 1`.

| character | `γ_1` | `M=40`: truth | bound | ratio | `M=60` ratio | `M=100` ratio |
|---|---|---|---|---|---|---|
| `q=3` | 8.040 | 1.71390e−2 | 8.3441e−3 | 0.487 | 0.746 | 0.782 |
| `q=4` | 6.021 | 1.60069e−2 | 1.12109e−2 | 0.700 | 0.762 | 0.794 |
| `d=5` (even) | 6.649 | 1.54242e−2 | 1.12317e−2 | 0.728 | 0.769 | 0.801 |
| `d=−7` | 4.476 | 1.42886e−2 | 1.08265e−2 | 0.758 | 0.788 | 0.814 |
| `d=8` (even) | 4.900 | 1.40833e−2 | 1.06950e−2 | 0.759 | 0.790 | 0.817 |
| `d=−8` | 3.576 | 1.39489e−2 | 1.06652e−2 | 0.765 | 0.793 | 0.818 |
| `d=12`, `13` (even) | 3.80, 3.12 | | | 0.778, 0.781 | 0.805, 0.808 | 0.829, 0.831 |
| `d=−43` | 0.836 | 1.09970e−2 | 9.0748e−3 | 0.825 | 0.843 | 0.859 |
| `d=−67` | 0.604 | 1.04752e−2 | 8.7604e−3 | 0.836 | 0.852 | 0.867 |
| `d=−148` | 0.505 | 9.6970e−3 | 8.2711e−3 | 0.853 | 0.867 | 0.879 |
| `d=−267` | 0.349 | 9.2152e−3 | 7.9556e−3 | 0.863 | 0.876 | 0.887 |
| `d=−232` | 0.311 | 9.3236e−3 | 8.0275e−3 | 0.861 | 0.873 | 0.885 |
| `d=−403` | 0.288 | 8.9174e−3 | 7.7556e−3 | 0.870 | 0.881 | 0.891 |
| **`d=−163`** (own zeros) | 0.2029 | 9.6132e−3 | 8.2170e−3 | 0.855 | 0.868 | 0.880 |

* **Against the direct certificate** (`cert/results_q{3,4}.txt`, `_stretch.txt`): the cell
  containing `M` gives `−Mδ′ ∈ [M|hi|, M|lo|]`. `q=3`: `M=40` `[1.6532e−2, 1.7908e−2]`,
  `M=60` `[1.2287e−2, 1.3310e−2]`, `M=100` `[8.681e−3, 9.404e−3]`; `q=4`: `[1.5412e−2, 1.6695e−2]`,
  `[1.1565e−2, 1.2526e−2]`, `[8.244e−3, 8.930e−3]`. Theorem E's bound (`8.344e−3`, `9.541e−3`,
  `7.024e−3`; `1.1211e−2`, `9.187e−3`, `6.787e−3`) lies below the certified lower end every time,
  and my truth values lie inside the certified intervals.
* **Conductor 163 at `M = 40`:** own Hardy-`Z` scan (109 zeros below 100, smooth count 109.32;
  first 84 agree with `bridge/zeros_q163.txt` to `2·10⁻¹³`): `σ² = 556.58`, `S_4/(16σ²) = 0.858`,
  `−Mδ′ = 9.6132e−3`, bound `8.2170e−3`. No conflict; low first zeros do not threaten the
  theorem (the worst case is `χ_3`, the smallest variance).
* **Abstract fuzz test** (`thmE_fuzz.py`): the proof uses only `c_j ≤ 4`, `S_4 ≤ 16σ²`, nine
  amplitudes `≥ 4/sqrt(26)`, and `σ²` in the range where `G > 0`. On 400 random amplitude sets
  with these properties (plus a Gaussian component) `f(1) + R ≥ G(σ²)/(sqrt(2π)σ)` always;
  smallest ratio `1.0187`. On the extremal family (`K` amplitudes equal to 4 plus a Gaussian of
  variance `8K`, so `S_4 = 16σ²`) the ratio is `1.370, 1.104, 1.038, 1.012, 1.003, 1.0007` for
  `σ² = 240, 256, 320, 960, 4000, 16000`: the inequality is asymptotically sharp in this class.
  An algebra slip of relative size `10⁻²` anywhere in step 5 would have shown up.

## 2. Proposition C₀, Corollary C′ (expanded proof), Lemma C″

**Proposition C₀: correct.** `|R| ≤ (1/2π)Σε_jc_j²∫t²|Φ^{(j)}| ≤ (Q_4/(32π))∫t²W¹`. On `[0,t_1]`
every `c_jt ≤ c_1t_1 ≤ z_1`, so Lemma L gives factorwise
`Φ ≥ exp(−σ²t²/2 − β(z_1)S_4t⁴) ≥ exp(−Σ²t²/2 − β(z_1)Q_4t⁴)` (`β > 0`), `cos t ≥ 0`; on
`[t_1,∞)`, `Φcos t ≥ −W`. `c_1` must be `sup_jc_j` (it is, by the notation of the document).

**Corollary C′: correct and uniform on `[M_a, M_b]`.** For fixed `M` in the cell:
(i) `c_i`, `σ²`, `S_4` increase termwise in `M`, so `Σ² = σ²(M_b)`, `Q_4 = S_4(M_b)` are
admissible and `t_1 ≤ z_1/c_1(M_b) ≤ z_1/c_1(M)`. (ii) `b(c_i(M)t) ≤ b(c_i(M_a)t)`. (iii) with
`m = min(t, z_*/c_n(M_a))` and `i > n`: `b(c_i(M)t) ≤ b(c_i(M_a)t) ≤ b(c_i(M_a)m) = exp(−c_i(M_a)²m²/4)`
because `c_i(M_a)m ≤ c_n(M_a)m ≤ z_*`; the product over `i > n` is `exp(−v_n(M_a)m²/2)`, and
with `j` removed it is that divided by `b(c_j(M_a)m)`. (iv) `b(c_1(M_a)t) ≤ b(c_j(M_a)t)`
(`j ≤ n`) and `b(c_1(M_a)t) ≤ b(c_j(M_a)m)` (`j > n`, since `c_j(M_a)m ≤ c_1(M_a)t`). So
`W = B_cell`, `W¹ = B¹_cell` satisfy the hypotheses of C₀ for every `M` in the cell, with the
same `MAIN_cell`. Case (b) needs: the `n` listed ordinates are all the ordinates in `(0, γ_n]`
(with multiplicity), and `γ_1` is the lowest ordinate. Both are stated.

**Lemma C″: correct** (single `M`). `|Φ^{(j)}| ≤ B/b(c_jt)` for `j ≤ n`; for `j > n`,
`≤ B/b(c_j min(t, z_*/c_n)) ≤ B/b(c_nt)` in case (b), `≤ B ≤ B/b(c_nt)` in case (a); the
weights `ε_jc_j²` are nonnegative and sum to `S_4/16`; `NUM″ ≤ NUM` because
`b(c_1t) ≤ b(c_jt)`. The lemma has no cell form and none is claimed; the certificate does not
use it.

## 3. Audit of `certify_envelope.py` and the results files

### 3.1 What it evaluates

Exactly the hypothesis of Corollary C′, case (b), `n = 12`, with **Theorem C's own `NUM`**
(`B¹ = B/b(c_1t)`; the product in the code runs over `ca_lo[1:]`), not the sharper `NUM″`.
My independent `NUM` (computed with `B¹`) is 0.6–0.9 % below the certified `NUM` in all 26
cells compared, which also excludes the `NUM″` variant.

### 3.2 Line-by-line findings

| item | finding |
|---|---|
| amplitudes | `c = 4/sqrt(1+(γ/M)²)`, each variable once. Used: lower endpoints of `c_i(M_a)` in `b` (correct, `b` nonincreasing), upper endpoint of `c_n(M_a)` in `r_lo = lo(lo(z_*)/hi(c_n))`, upper endpoint of `c_1(M_b)` for `t_1`. |
| `log_b_upper` | `z = lo(c_lo·t)` is an exact lower bound of the argument. Four exact comparisons: `z < lo(z_*)` → `−z²/4`; `z > hi(z_far)` → far piece; `z > hi(z_*)` → `log ρ`; otherwise `−lo(z_*)²/4`. Each is `≥ log b(z) ≥ log b(argument)`. Correct. |
| tail factor | `exp(−v_lo·min(t, r_lo)²/2) ≥ exp(−(v/2)min(t, z_*/c_n)²)` since `v ≥ v_lo > 0`, `z_*/c_n ≥ r_lo`. Correct. `v_lo > 0` is required. |
| `NUM` | `Σ t_{k+1}²·B1_hi[k]·(t_{k+1}−t_k)`, an upper sum because `B¹` is nonincreasing; plus `rem2`; then `hi(hi(S_4(M_b))/(32π)·…)`. Correct. |
| far remainder | requires (verified in Arb) `lo(c_10(M_a))·T_end > hi(z_far)`; nine factors `i = 2..10` in the far piece, factors 11, 12 bounded by 1, tail factor frozen at `min(T_end, r_lo)`: `rem2 = (8/π)^9Π_{i=2}^{10}c_i^{−1/2}·e^{…}·T_end^{−3/2}/(3/2)`, `rem0 = …T_end^{−7/2}/(7/2)` (`B ≤ B¹`). Exponents and index range (`ca_lo[1:10]`) checked. Correct. |
| `TAILF` | suffix sums of `B_hi[k]·Δ_k` plus `rem0`; upper. Correct. |
| `MAIN` | lower Riemann sum `Σ_{k=1}^{k_1} lo(G(t_k))(t_k − t_{k−1})` of a positive decreasing integrand, with `hi(σ²(M_b))`, `hi(β(z_1))`, `hi(S_4(M_b))`; early stop only drops positive terms. `t_1 < lo(π/2)` and `hi(hi(c_1(M_b))t_1) ≤ z_1` verified in Arb after the float guess. `β` from `bessel_j` at an exact rational, `lo(J_0) > 0` required. Correct. |
| decision | `LB = lo(lo(MAIN − TAILF) − NUM)`, `PASS` iff `LB > 0`. `K` is display only. Correct. |
| floats | only: grid placement, `T_end`, the guess of `k_1`, the early stop, printing. Every float-chosen quantity is an exact dyadic rational (exactness asserted) and is re-verified in Arb where it matters. Mutant m19 (float cap ×1.2) is refused by the Arb check. No float on a certified path. |
| zero balls | `arb(decimal, "1e-20")` at 128 bits: I checked for all 1048 lines that the ball contains `[d − 10⁻²⁰, d + 10⁻²⁰]` exactly (rational comparison at 400 bits); radius `1.0000000036e−20`. Strict ordering checked on load. |
| `σ²`, `S_4` for `s ≥ 4` | `loglam_derivs_series`: signs of `(L′/L)`, `(L′/L)′`, the integral-test tails, `ψ′ = ζ(2,·)`, `S_4 = 64Mh − 64M²h′` all correct. My mpmath values (via `dirichlet`) agree to all printed digits and lie on the correct side of the printed bounds. |
| `σ²`, `S_4` for `s < 4` | `h` by complex step with Cauchy remainder `M_L h²/(R³(1−(h/R)²))`, `M_L = (s+R)/(s−R)` from Abel summation with partial sums bounded by 1 (true for `χ_3`, `χ_4`). `S_4 ≤ Σ_{list}c⁴ + c_N²·2V_T`, `V_T = σ² − Σ_{list}c²/2 = Σ_{unlisted}c²/2`. **Rigorous given the unseen tail**: under (H1)–(H3) every unlisted zero has `γ > γ_N`, so `c ≤ c_N` and `Σ_{unlisted}c⁴ ≤ c_N²Σ_{unlisted}c²`; the tail variance is exact from the variance identity. It needs completeness **and simplicity** of the whole list (a double zero would hide a large amplitude in `V_T`); both are in (H3). The enclosure is tight (`S_4` upper bound exceeds my value by `< 10⁻⁵` relative). |
| results files | 765 lines each, all `PASS`, no `ERROR`/`FAIL`, in order, no duplicate; adjacent cells share endpoints exactly; first endpoint 1, last 2000; `n=12` on every line; 127 list-mode cells (126 `list/list`, 1 `list/series`); printed values mutually consistent (`referee5/check_results.out`). |
| reproducibility | full fresh rerun of both ranges (75 s + 91 s): identical to the stored files except `secs` and the timestamp. `selftest`, `points_q*.txt`, `consistency_q*.txt`, `analytic_leg.out` also reproduce identically. |

### 3.3 Independent recomputation (`referee5/ref_lib.py`, `compare_cells.py`)

Own mpmath code (30 digits, piecewise Gauss quadrature between all breakpoints, closed-form
far tail, closed forms through `mpmath.dirichlet`), evaluating the exact cell quantity
`E = MAIN_cell − TAILF_cell − NUM_cell` with the certificate's `(z_1, t_1)`. 13 cells per modulus:

| cell | `q=4`: `LB` (cert) | `E` (mine) | `LB/E` | `q=3`: `LB` | `E` | `LB/E` |
|---|---|---|---|---|---|---|
| `[1, 1.01]` (list) | 1.467219e−1 | 1.482706e−1 | 0.9896 | 1.032459e−1 | 1.050308e−1 | 0.9830 |
| `[1.2942, 1.3071]` | 1.906912e−1 | 1.917813e−1 | 0.9943 | 1.794049e−1 | 1.806523e−1 | 0.9931 |
| `[1.9827, 2.0025]` | 1.782553e−1 | 1.789278e−1 | 0.9962 | 2.023387e−1 | 2.030377e−1 | 0.9966 |
| `[3.4578, 3.4923]` (list) | 1.083973e−1 | 1.088750e−1 | 0.9956 | 1.375163e−1 | 1.379533e−1 | 0.9968 |
| `[3.4923, 3.5272]` (list/series) | 1.073323e−1 | 1.078072e−1 | 0.9956 | 1.362480e−1 | 1.366823e−1 | 0.9968 |
| `[4.4773, 4.522]` | 8.341092e−2 | 8.382594e−2 | 0.9950 | 1.070335e−1 | 1.074115e−1 | 0.9965 |
| `[7.9698, 8.0494]` | 4.896672e−2 | 4.925764e−2 | 0.9941 | 5.974664e−2 | 6.002489e−2 | 0.9954 |
| `[17.141, 17.3124]` | 2.750565e−2 | 2.767645e−2 | 0.9938 | 3.059705e−2 | 3.077051e−2 | 0.9944 |
| `[39.9288, 40.328]` | 1.525471e−2 | 1.535345e−2 | 0.9936 | 1.633693e−2 | 1.643806e−2 | 0.9938 |
| `[99.7281, 100.7253]` | 8.279664e−3 | 8.336565e−3 | 0.9932 | 8.704450e−3 | 8.762634e−3 | 0.9934 |
| `[397.6236, 401.5998]` | 3.449757e−3 | 3.475381e−3 | 0.9926 | 3.573136e−3 | 3.599247e−3 | 0.9927 |
| `[993.1969, 1003.1288]` | 1.980509e−3 | 1.995811e−3 | 0.9923 | 2.039159e−3 | 2.054713e−3 | 0.9924 |
| `[1973.3722, 1993.1059]` | 1.318804e−3 | 1.329219e−3 | 0.9922 | 1.353296e−3 | 1.363900e−3 | 0.9922 |

In every cell also: certified `MAIN ≤` mine, `TAILF ≥` mine, `NUM ≥` mine, `sig2b`, `S4b ≥` mine,
`v ≤` mine. The 0.4–1.7 % gap is the Riemann-sum slack.

**Against the truth where no direct certificate exists** (`[1, 3.6]`, `lowM_truth.py`, full
lists, non-rigorous): certified `LB`, `MAIN − TAILF`, `NUM` versus my `−Mδ′`, `f(1)`, `|R|`:
`LB/truth = 0.625, 0.690, 0.766, 0.842, 0.899, 0.890, 0.882, 0.869` at
`M = 1, 1.1, 1.25, 1.5, 2, 2.5, 3, 3.5` (`q = 4`) and `0.508, 0.611, 0.717, 0.821, 0.917, 0.940, 0.941, 0.930`
(`q = 3`); always `MAIN − TAILF ≤ f(1)` and `NUM ≥ |R|`. My truth values reproduce the column
"true `−Mδ′`" of the table in `BRIDGE-ATTEMPT.md` §3 to within one unit of the last printed digit.

### 3.4 Mutation tests (`referee5/mutation_tests.py`, on copies in `referee5/mut/`)

13 cells, `q = 4`. "Detected by `E`" = the mutant prints `LB` above my independent value.

| mutant | outcome | detected by |
|---|---|---|
| unmutated | 13 `PASS`, max `LB/E = 0.9962` | — |
| sign of `NUM` flipped | 13 `PASS`, `LB/E` up to 2.56 | `E` (13/13); direct certificate (5 cells) |
| `TAILF` dropped | 13 `PASS`, `LB/E` up to 1.79 | `E` (6/13, the cells with `M < 5`); **not** by `consistency.py` |
| far remainder dropped | 13 `PASS`, unchanged to 5 digits | nothing (remainder `≤ 3·10⁻⁵`); covered by code reading |
| `NUM` as a lower sum | `LB/E` up to 1.0021 | `E` (7/13) |
| `MAIN` as an upper sum | `LB/E` up to 1.0031 | `E` (6/13) |
| tail factor ×3 | `LB/E` up to 1.62 | `E` (13/13); direct certificate (5) |
| `σ²`, `S_4` at `M_a` instead of `M_b` | `LB/E` up to 1.015 | `E` (12/13) |
| `ρ = 0.30` (false envelope) | `LB/E` up to 1.026 | `E` (6/13) |
| zero radius 0 / `10⁻³` | 13 `PASS`, `LB ≤ E` | (sound either way) |
| zero radius 0.2 | 13 refused: balls not ordered | program (fails closed) |
| zeros 4, 5 swapped | 13 refused | program |
| `γ_3 + 10⁻³` | 13 `PASS` | `selftest` refuses (`S_4` list vs series) |
| `γ_3 + 10⁻⁶` | 13 `PASS` | nothing (below every tolerance; this is hypothesis (H3)) |
| zero 5 deleted | 13 `PASS` | `selftest` refuses |
| `γ_1 := 0.2029` | 6 refused (`v` or `V_T` not positive), 7 `PASS` | `selftest` refuses |
| precision 40 bits | 13 `PASS`, `LB ≤ E` | (sound) |
| NaN injected in `σ²` | 13 refused | program |
| `t_1` beyond admissible | 13 refused | program |

Conclusions: every refusal path fails closed. Logic errors that raise `LB` by more than the
quadrature slack (0.4–1.7 %) are caught by my independent evaluation. The cell program cannot detect a wrong zero list by
itself (in series mode it never compares the list with the closed forms); that is exactly
hypothesis (H3), and `selftest` covers gross errors only. The comparison with the direct
certificate does not exercise `TAILF`.

## 4. The chain for `q = 3, 4`

Claim: Theorem B on `(0, M_0]`, the envelope certificate on `[1, 2000]`, Theorem E on
`[40, ∞)` give `δ_q` strictly decreasing on `(0, ∞)` without `thm:monotone-large`.

* **Differentiability:** `δ ∈ C¹(0,∞)` by Lemma D. Needed only on `[1, ∞)`.
* **`(0, M_0]`:** Theorem B (certified in `REFEREE-1.md` §4.3; I re-derived it: for fixed
  phases `X_M/M − X_{M′}/M′ ≥ ∫_M^{M′}(1 − Q(u))u^{−2}du > 0` when `M < M′ ≤ M_0`; strictness from
  the full support of `S`, since `Σc_j = ∞`). It is a statement about the model, so LI is not
  used.
* **`[1, 2000]`:** `δ′ ≤ −LB/M_b < 0` on each of 765 contiguous closed cells.
* **`[40, ∞)`:** Theorem E, `δ′ < 0`.
* **Overlaps:** `M_0 > 3.68` resp. `4.47` (`cert/aux_results.txt`; my float values
  `Q_{χ_4}(3.68) = 0.99409`, `Q_{χ_3}(4.47) = 0.99211` agree with the certified enclosures).
  Only `M_0 > 1` is needed, and that follows from numbers already inside the envelope
  certificate: `Q = Σc_j³/16 ≤ c_1σ²/8`, and at `M = 1.01`: `c_1 ≤ 0.671`, `σ² ≤ 0.6292`
  (`q = 4`), `c_1 ≤ 0.503`, `σ² ≤ 0.4599` (`q = 3`), so `Q(1.01) ≤ 0.053` resp. `0.029`. The
  certificate and Theorem E overlap on `[40, 2000]`.
* **Strictness:** strictly decreasing on `(0, M_0]` and `δ′ < 0` on `[1, ∞)`, with
  `1 < M_0`: strictly decreasing on `(0, ∞)`.

**The chain is complete.** Theorem E is needed only for `M > 2000`. `thm:monotone-large`,
`prop:criterion-eval`, `lem:analytic-leg` are not used.

**Exact trusted inputs.**

1. GRH(`χ_q`), `q ∈ {3,4}`. (`L(1/2,χ_q) ≠ 0` is unconditional and also follows from the zero
   count below. LI only to read `δ` as the density of the race.)
2. The variance identity `σ² = 4M(log ξ)′(M+1/2)` and its termwise `M`-derivative
   (`S_4 = 64Mh − 64M²h′`): paper `prop:varid`, `eq:varidq`, `lem:Dprop`(i); re-derived in §1.1.
3. Zero lists `prototype/explore/chi4_zeros.txt` (511), `chi3_zeros.txt` (537); SHA-256
   `0e73e00a…ddfbc`, `e3075235…8aaa`, matching `cert-audit/SHA256SUMS`. Certification:
   `prototype/explore/certify_arb.py --K 0` (sign change of Hardy's `Z` in each `±10⁻²⁰`
   window; winding number of `L` around `(−1/2, 3/2) × (−1/2, T)` equal to the list length,
   hence complete and simple); log `cert-audit/certify_arb_K0.log`, rerun by Referee 4. I did
   not rerun it; I spot-checked the 12 lowest ordinates of each list (12/12 certified sign
   changes; exactly 12 sign changes of `Z` below `γ_12 + 0.5` on a 0.01 grid). What the
   envelope certificate uses of it: the 12 lowest ordinates, the absence of other zeros with
   `0 < γ ≤ γ_12`, simplicity; and for the 127 cells with `M_a < 3.5` the whole list, complete
   and simple.
4. FLINT 3.6.0 through python-flint 0.9.0: ball arithmetic, `exp`, `log`, `sqrt`, `cos`, `π`,
   `digamma`, Hurwitz `ζ(2,·)`, `bessel_j(0)` at 20 exact rationals, and for `M < 3.5`
   `dirichlet_char.l` at complex argument.
5. The shared module `cert/certify_derivative.py` (`load_zeros`, `loglam_derivs_series`,
   `loglam_d1_complex_step`, `make_grid`, printing); read in this audit.
6. Pencil proofs: Lemma D, Lemma J, Lemma L, Theorem C (Referee 3), with the paper's
   `lem:phi`(i),(iii), `lem:J0env`(i),(ii), (J2), `lem:invert`(b); Proposition C₀ and
   Corollary C′ (this report); Theorem B (Referee 1; re-derived here).
7. For `M > 2000`: Theorem E (this report) and four Arb evaluations of elementary functions.

## 5. The general theorem that follows

**Theorem (any real primitive `χ`).** Let `χ` be a real primitive character mod `q ≥ 3`,
assume GRH(χ) and `L(1/2,χ) ≠ 0`, let `γ_j` be the positive ordinates with multiplicity,
`δ_χ(M) = P(1 + Σ_jc_j(M)cos θ_j > 0)` with independent uniform phases, and let `M_0(χ)` be the
unique root of `Q_χ(M) = Σ_j4M³/(M²+γ_j²)^{3/2} = 1`. Then

1. `δ_χ ∈ C¹(0, ∞)`;
2. `δ_χ` is strictly decreasing on `(0, M_0(χ)]`;
3. for every `M ≥ 40`, `−Mδ_χ′(M) ≥ G(σ²(M))/(sqrt(2π)σ(M)) ≥ 0.3219/(sqrt(2π)σ(M)) > 0`; in
   particular `δ_χ` is strictly decreasing on `[40, ∞)`;
4. `M_0(χ) < 40` always, so the middle range `(M_0(χ), 40)` is a nonempty interval on which
   nothing is asserted.

Under LI(χ), `δ_χ(M)` is the logarithmic density of the weighted race.

**Supported by what is proved:** yes. (1) is Lemma D; (2) is Theorem B; (3) is Theorem E.
(4) is not in the documents but is immediate: `c ≤ 4` gives `Q = Σc³/16 ≥ S_4/64`, and at
`M = 40`, `S_4 = 16σ² − 64·1600·h′ ≥ 16·236.97 − 102400·0.0127 > 2490`, so `Q(40) > 38`.

**The increase in the middle range is not a theorem of the documents under review.** There it
is numerical evidence: conductor 163, `M_0 = 0.164`, `δ(0.6) = 0.59033 < δ(0.9) = 0.59414 < δ(1.2) = 0.59603`,
`δ(2) = 0.58664`. I reproduced these four values with my own zeros and closed forms
(`referee5/check_163_lowM.py`). A rigorous enclosure is claimed in
`monotonicity/nonmonotone/NONMONOTONE.md`, which appeared during this review; I have not
refereed it (my float values lie inside its enclosures). The accurate wording is therefore:
"monotone on `(0, M_0(χ)]` and on `[40, ∞)`; open on `(M_0(χ), 40)` in general; numerically
(and, subject to an unrefereed certificate, rigorously) not monotone there for conductor 163".

## 6. Findings

### FATAL
None.

### REPAIRABLE

**R1. Displayed inequalities of Theorem E's proof that are false in the last digit.**
`analytic_leg.py` prints lower endpoints with `:.6f`/`:.4f`/`:.3f` (round to nearest) under the
labels `>=`, `<=`; `ANALYTIC-LEG.md` copies them, some with further rounding. With directed
rounding (`referee5/analytic_leg_digits.out`) the certified values are:

| `M_E` | `κ` | `Y` | (E2) ≥ | `σ_0²` ≥ | `N` ≥ | `T_c` ≤ | `G(σ_0²)` ≥ |
|---|---|---|---|---|---|---|---|
| 300 | 4 | 100 | 0.0016583 | 2979.71089 | 70.3092 | 6.684128 | 0.4974001 |
| 100 | 4 | 40 | 0.0049261 | 774.18188 | 19.1144 | 6.684128 | 0.4899184 |
| 60 | 4 | 48 | 0.0081300 | 403.61089 | 14.9398 | 6.684128 | 0.4804794 |
| 50 | 4 | 40 | 0.0097087 | 318.27738 | 11.3751 | 6.684128 | 0.4750812 |
| 45 | 5 | 40 | 0.0107526 | 277.06752 | 12.9847 | 8.266220 | 0.4681158 |
| **40** | **5** | **40** | **0.0120481** | **236.97110** | **11.3751** | **8.266220** | **0.3219715** |
| 38.5 | 5 | 28 | 0.0124999 | 225.17929 | 9.0645 | 8.266220 | 0.0249902 |

False as written in `ANALYTIC-LEG.md`:

* Step 6: "(E5): `G(236.9711) ≥ 0.321972`" (certified lower endpoint `0.3219715`) and
  "`T_c ≤ 8.2662`" (`T_c = 8.266219…`); Constants: "`T_c = 8.2661…`" (it is `8.2662…`).
* Table, row 40 (the row of the theorem): `G ≥ 0.3220`, `T_c ≤ 8.2662`.
* Table, other rows: `N ≥ 70.31` (300); `N ≥ 14.94`, `G ≥ 0.4805` (60); (E2) `≥ 0.009709`,
  `σ_0² ≥ 318.28`, `G ≥ 0.4751` (50); (E2) `≥ 0.010753`, `σ_0² ≥ 277.07` (45); (E2) `≥ 0.012500`,
  `σ_0² ≥ 225.18`, `N ≥ 9.065`, `G ≥ 0.0250` (38.5); `T_c ≤ 6.6841` (four rows), `T_c ≤ 8.2662`
  (three rows).

`analytic_leg.out` has the same defect in 20 of its 42 printed numbers (e.g. `0.321972`,
`2979.7109`, `14.940`, `0.468116`, and `T_c` on every line). Each excess is at most one unit of the last printed digit.

**Not affected:** the script's decisions (they use the balls); in the proof text
`σ_0² ≥ 236.9711`, `N ≥ 11.375`, `β_* ≤ 0.025553`, (E2) margin `≥ 0.01204`; the theorem's
constant `0.3219`; the constant `0.475` quoted for `M ≥ 50` in `BRIDGE-ATTEMPT.md` §3A.
*Fix:* print lower endpoints rounded down and upper endpoints rounded up (as
`certify_envelope.py` does with `e_dn`/`e_up`), regenerate `analytic_leg.out`, and replace the
table and the two sentences by the values above (e.g. "(E5): `G(236.9711) ≥ 0.32197`, with
`β_* ≤ 0.025553`, `T_c ≤ 8.26622`").

### MINOR

1. **Sub-range maxima of `K` omit the cell that straddles the right endpoint.** To cover
   `[3.6, 20]` (`q=4`) the bound is `0.3614` (cell `[19.8995, 20.0984]`), not `0.3612`;
   `[4.4, 20]` (`q=3`): `0.3364`, not `0.3359`; up to `M = 400`: `0.4219` / `0.4169`, not
   `0.4217` / `0.4168`. (`BRIDGE-ATTEMPT.md` §3, `METHOD.md` §6, `SUMMARY.md`.)
2. **"The bound loses 4–17 %" (`METHOD.md` §8.4; "at most 17 %" in `BRIDGE-ATTEMPT.md` §3)**
   holds on the common range with the direct certificate only. On `[1, 3.6]` the loss is up to
   38 % (`q=4`, `M=1`) and 49 % (`q=3`, `M=1`) against my numerical truth. The bound is still
   positive by a wide margin there.
3. **`suf[len−1]` is left at 0** in `certify_cell`; it should be `rem0`. Unreachable
   (`t_1 < π/2 < 8.15 < T_end`), so no certified value is affected. One-line fix:
   initialise `suf[-1] = rem0`.
4. **"The smallest threshold this argument reaches is `M ≥ 38.5`"** (`BRIDGE-ATTEMPT.md` §3A):
   it is the smallest on the half-integer grid, for `κ ∈ {3,…,7}`, even `Y`, and the fixed
   constants. `ANALYTIC-LEG.md` words it correctly.
5. **Multiplicities.** Theorem E counts a zero of multiplicity `m` as `m` independent
   amplitudes `c ≤ 4`; that is the stated model and the proof is correct for it. If a multiple
   zero were modelled as one amplitude `mc`, the step `c_jτ ≤ z_*` would fail. Under LI zeros are
   simple, so the race reading is unaffected. A sentence would prevent a misreading.
6. **The count is 12, not 11.** `N_χ(200) ≥ 11.375` and `N` is an integer; the "why 40" paragraph
   says "certifies 11 zeros".
7. **`consistency.py` does not exercise `TAILF`** (mutant "`TAILF` dropped" passes it) nor the
   far remainder. It validates `MAIN − NUM` against the direct certificate. The independent
   evaluation of §3.3 covers `TAILF`.
8. **Resume logic.** `done_cells` accepts any existing `PASS` line without recording the hash
   of the script or of the zero files. No stale line is present (fresh rerun identical).
   Suggest writing the SHA-256 of `certify_envelope.py`, `certify_derivative.py` and the zero
   file into the header.
9. **Series-mode cells never compare the zero list with the closed forms** beyond `v > 0`;
   `selftest` does so at five values of `M` and tolerates an ordinate error of about `10⁻⁴`.
   This is hypothesis (H3), correctly listed, but "fail-closed" should not be read as covering
   a wrong zero file.
10. **"Share … nothing else"** (`BRIDGE-ATTEMPT.md` §0 item 2: "the zero files and the
    closed-form evaluation of `σ²`"). They also share the evaluation of `h′`, `load_zeros`,
    `make_grid` and the helpers; `METHOD.md` §7 and `SUMMARY.md` caveat 4 say so.
11. **`SUMMARY.md` "under GRH"** for the consequence: the computational hypotheses (H2)–(H5)
    belong in the same sentence, as in `BRIDGE-ATTEMPT.md` §3.
12. **`Z_1`, `Z_2` "for `s ≥ 4`"** is stronger than needed (`s > 2/log 2`); harmless.
13. **Lemma C″** is proved for a single `M` only and is not used by any certificate; say so
    next to the statement.

## 7. What I certify, and what I do not

**Certified as proved as written**

* Proposition C₀.
* Corollary C′ with its expanded proof, cases (a) and (b), uniformly on `[M_a, M_b]`.
* Lemma C″ (single `M`).
* **Theorem E**: statement, proof, and the four Arb inequalities, with the constant `0.3219`,
  for every real primitive `χ` mod `q ≥ 3`, both parities, GRH(χ), `L(1/2,χ) ≠ 0`, every
  `M ≥ 40`, zeros counted with multiplicity as independent amplitudes. Also the rows
  `M_E = 300, 100, 60, 50, 45, 38.5` of the table as decided by the script (with the corrected
  numbers of R1).
* The certificate: under (H1)–(H5), `δ_q′(M) ≤ −LB/M_b < 0` for every `M ∈ [1, 2000]`,
  `q ∈ {3, 4}`; the program evaluates the hypothesis of Corollary C′ soundly; the results files
  are complete, contiguous, all `PASS`, and reproducible.
* The chain of §4 and the general theorem of §5 (items 1–4).

**Not certified**

* The displayed numbers listed in R1 as printed (corrected values given).
* The sub-range maxima of minor item 1 and the "4–17 %" of minor item 2 as printed.
* Anything that is an input: the zero-list certificate (not rerun by me), FLINT/Arb, the
  paper's lemmas, `cert/aux_results.txt` (not needed: `M_0 > 1` suffices and follows from the
  certificate's own numbers).
* Non-monotonicity for conductor 163 or 232 as a theorem: in the documents under review it is
  numerical; the certificate in `monotonicity/nonmonotone/` is outside this review.

**Corrected statements**

* `ANALYTIC-LEG.md`, Step 6: "(E5): `G(236.9711) ≥ 0.32197` (with `β_* ≤ 0.025553`,
  `T_c ≤ 8.26622`)"; Constants: "`T_c = 8.2662…`"; table as in R1.
* `BRIDGE-ATTEMPT.md` §3 / `METHOD.md` §6 / `SUMMARY.md`: "`K_cell ≤ 0.3614` on the cells
  covering `[3.6, 20]` (`q = 4`), `≤ 0.3364` on those covering `[4.4, 20]` (`q = 3`); `≤ 0.4219`
  resp. `0.4169` on those covering `M ≤ 400`."
* `METHOD.md` §8.4: "on the common range with the direct certificate the bound loses 4–17 % of
  `|δ′|`; near `M = 1` it loses up to one half."
* `BRIDGE-ATTEMPT.md` §3A: "on the searched grid the smallest threshold is `M ≥ 38.5`."
* Consequence (`SUMMARY.md`): "Under GRH(`χ_q`), the variance identity, the certified zero
  lists and FLINT/Arb: `δ_q` is strictly decreasing on `(0, ∞)` for `q = 3, 4` (Theorem B on
  `(0, 1.01]`, the envelope certificate on `[1, 2000]`, Theorem E on `[40, ∞)`); LI is needed
  only to read `δ_q` as the density of the race."

## 8. Referee files (`monotonicity/referee5/`)

| file | purpose |
|---|---|
| `ref_lib.py`, `compare_cells.py`, `compare_cells.out` | independent mpmath evaluation of the cell quantity; 13 cells per modulus |
| `check_results.py`, `check_results.out` | structure of the results files (all `PASS`, tiling, internal consistency) |
| `rerun_q4.txt`, `rerun_q3.txt`, `*.rerun.*` | full fresh reruns (identical to the stored outputs up to `secs`/timestamp) |
| `mutation_tests.py`, `mutation_tests.out`, `mut/` | 19 mutants of the prover and one unmutated copy |
| `thmE_numeric.py`, `.out` | Theorem E vs numerical `−Mδ′` for `q = 3, 4`, conductor 163 (own zeros, `myzeros_163.txt`) and six low-zero characters; comparison with the direct certificate |
| `thmE_even.py`, `.out` | even characters; minimum of `σ²(40)` and of the count over `|d| ≤ 60` |
| `thmE_fuzz.py`, `.out`, `thmE_fuzz_extreme.out` | fuzz test of the abstract inequality; extremal family |
| `analytic_leg_digits.out` | Theorem E constants with directed rounding (R1) |
| `lowM_truth.py`, `.out` | certified `LB` vs numerical truth on `[1, 10]`; `Q(1)`, `Q(M_0)` |
| `check_163_lowM.py` | conductor 163, low `M`, own zeros |
| `first12_check.out` | sign changes of `Z` at the 12 lowest ordinates |

All numerical "truth" values in this report are floating point (Simpson, double precision);
they are tests, not certificates. Total CPU used: about 8 minutes, one core, `nice -n 10`.
