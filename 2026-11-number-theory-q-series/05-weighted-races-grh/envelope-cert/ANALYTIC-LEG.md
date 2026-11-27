# Theorem E: δ′(M) < 0 for every real primitive character and every M ≥ 40, without zero data

Status: complete pencil proof below; the four numerical inequalities it reduces
to are evaluated in Arb by `analytic_leg.py` (output `analytic_leg.out`).
**New on 2026-11-24; refereed in `../REFEREE-5.md` (internal, AI-conducted: proved
as written; the printing defect R1 and the minor items are repaired here, see
"Review history" at the end).** Notation and the results quoted
(Lemma D, Lemma J, Lemma L, Proposition C₀) are those of
`monotonicity/BRIDGE-ATTEMPT.md` §1–§2.

## Statement

Let `χ` be a real primitive Dirichlet character mod `q ≥ 3`, of parity
`a ∈ {0, 1}` (`χ(−1) = (−1)^a`); assume GRH(χ) and `L(1/2, χ) ≠ 0`. Let
`0 < γ_1 ≤ γ_2 ≤ …` be the positive ordinates of `L(s,χ)` (with multiplicity: a
zero of multiplicity `m` is listed `m` times and counts as `m` independent
amplitudes, each `≤ 4`; see "What is and is not claimed"),
`M > 0`, `s = M + 1/2`,

    c_j = 4M/sqrt(M²+γ_j²),  σ² = Σ c_j²/2,  S_4 = Σ c_j⁴,  Φ = Π J_0(c_j t),  δ(M) = 1/2 + (1/π)∫_0^∞ Φ sin t/t dt.

Constants: `ρ = 0.446`, `z_* = 2 sqrt(log(1/ρ)) = 1.79714…`, `τ := z_*/4 = 0.44928…`
(so `ρ = e^{−z_*²/4}` and `e^{−xτ²/2} = ρ^{x/8}`), `κ = 5`,
`c_κ := 4/sqrt(1+κ²)`, `T_c := (8/π)²/c_κ = 8.2662…`,
`β_* := β(z_*) = log(ρ/J_0(z_*))/z_*⁴ ≤ 0.025553`, and

    G(x) := 1 − (1/2 + 48β_*)/x − (1/2)(1 − 8/x)^{−3/2}
              − sqrt(2πx)·[ ρ^{x/8}/(π x τ) + (9/(7π))·T_c·ρ^{(x−72)/8} ]
              − (x^{3/2}/sqrt(2π))·(4/3)·T_c³·ρ^{(x−80)/8}.

**Theorem E.** For every such `χ` and every `M ≥ 40`:

    −M·δ′(M) ≥ G(σ²(M)) / (sqrt(2π)·σ(M)) ≥ 0.3219 / (sqrt(2π)·σ(M)) > 0.

In particular `δ_χ` is strictly decreasing on `[40, ∞)`.

The theorem uses no ordinate of any zero. Inputs: the variance identity
(`prop:varid`, `eq:varidq` of the paper, valid for either parity),

    σ²(M) = 4M·h(M),   h(M) := (log Λ)′(s) = ½ log(q/π) + ½ ψ((s+a)/2) + (L′/L)(s),            (V)

its termwise `M`-derivative (as in `METHOD.md` (H2)),

    S_4(M) = 64M·h(M) − 64M²·h′(M),   h′(M) = ¼ ψ′((s+a)/2) + (L′/L)′(s),                      (V′)

and, for real `s ≥ 4` (more than is needed: `s > 2/log 2 = 2.885…` suffices for
the monotonicity of the integrands below), the Dirichlet-series bounds

    |(L′/L)(s)|  ≤ Σ_{n≥2} log n · n^{−s}  ≤ Z_1(s) := L 2^{−s} + 2^{1−s}(L/(s−1) + 1/(s−1)²),
    |(L′/L)′(s)| ≤ Σ_{n≥2} log² n · n^{−s} ≤ Z_2(s) := L² 2^{−s} + 2^{1−s}(L²/(s−1) + 2L/(s−1)² + 2/(s−1)³),

`L = log 2` (first term plus the integral over `[2, ∞)`; the integrands
`log x·x^{−s}`, `log²x·x^{−s}` are decreasing on `[2, ∞)` because
`log 2 > 2/s`). `Z_1`, `Z_2` are decreasing in `s`.

## Step 1. `S_4 ≤ 16σ²` for `M ≥ 40`

By (V), (V′): `S_4 = 16σ² − 64M²h′`. `ψ′(x) = Σ_{k≥0}(x+k)^{−2} ≥ ∫_x^∞ u^{−2}du = 1/x`
and `ψ′` is decreasing, so `h′ ≥ ¼ψ′((s+1)/2) − Z_2(s) ≥ 1/(2(s+1)) − Z_2(s)`.
`(s+1)·Z_2(s)` is decreasing for `s ≥ 4`, so it suffices that
`1/(2(s+1)) > Z_2(s)` at `s = 40.5`: **(E2)**, margin `≥ 0.01204`. Hence `h′ > 0`
and `S_4 ≤ 16σ²` for all `M ≥ 40`.

## Step 2. `σ²(M) ≥ σ_0² := 236.97` for `M ≥ 40`

`σ²(M) = Σ 8M²/(M²+γ²)` is increasing in `M`, so `σ²(M) ≥ σ²(40) = 160·h(40)`, and
by (V), `h(40) ≥ ½log(q/π) + ½ψ((40.5+a)/2) − Z_1(40.5)`. Odd characters have
`q ≥ 3`, even ones `q ≥ 5`. **(E3)**: the minimum over the two families
`(a, q_min) ∈ {(1, 3), (0, 5)}` of `160·[½log(q_min/π) + ½ψ((40.5+a)/2) − Z_1(40.5)]`
is `≥ 236.9711`.

## Step 3. At least nine ordinates are `≤ 5M`, for `M ≥ 40`

**Counting lemma.** For `0 < Y` and `Γ > 0`:
`N_χ(Γ) := #{j : γ_j ≤ Γ} ≥ (Y/2)·[ h(Y) − (2Y/Γ)·h(Γ) ]`.

*Proof.* In `σ²(Y) = Σ_j 8Y²/(Y²+γ_j²)` each zero with `γ_j ≤ Γ` contributes at
most 8. For `γ > Γ`: `8Y²/(Y²+γ²) ≤ (2Y²/Γ²)·8Γ²/(Γ²+γ²)`, because this is
equivalent to `Γ² + γ² ≤ 2Y² + 2γ²`, true since `γ > Γ`. Summing,
`σ²(Y) ≤ 8N_χ(Γ) + (2Y²/Γ²)σ²(Γ)`; insert (V): `σ²(Y) = 4Y h(Y)`,
`σ²(Γ) = 4Γ h(Γ)`. ∎

Take `Y = 40`, `Γ = 200`, `w := 2Y/Γ = 2/5`. By (V),

    h(Y) − w·h(Γ) ≥ (1−w)·½log(q/π) + ½ψ((Y+½+a)/2) − (w/2)ψ((Γ+½+a)/2) − Z_1(Y+½) − w·Z_1(Γ+½),

and `(1−w) ≥ 0`, so the right side is smallest at `q = q_min`. **(E4)**: for both
families the lower bound of `(Y/2)[…]` is `≥ 11.375`. Hence `N_χ(200) ≥ 12` (it
is an integer), in particular `≥ 9`, for every `χ`, and for `M ≥ 40` the nine lowest ordinates satisfy
`γ_i ≤ 200 ≤ 5M`, i.e.

    c_i(M) ≥ c_κ = 4/sqrt(26)        (i = 1, …, 9).

## Step 4. Envelopes

Let `𝓛 = {1,…,9}`. Every `c_j ≤ 4`, so `c_j²/2 ≤ 8`.

*`0 ≤ t ≤ τ`.* Then `c_j t ≤ 4τ = z_* < j_{0,1}` for all `j`, and
`0 < J_0(c_j t) ≤ e^{−c_j²t²/4}` (positivity of the Taylor coefficients of
`−log J_0`). Hence

    |Φ(t)| ≤ e^{−σ²t²/2},        |Φ^{(j)}(t)| ≤ e^{−(σ² − c_j²/2)t²/2} ≤ e^{−(σ²−8)t²/2}   (every j).

*`t ≥ τ`.* For `j ∉ 𝓛`: `|J_0(c_j t)| ≤ b(c_j t) ≤ b(c_j τ) = e^{−c_j²τ²/4}` (Lemma J;
`b` nonincreasing; `c_j τ ≤ z_*`). For `i ∈ 𝓛`: `|J_0(c_i t)| ≤ p(t) :=
min(1, (8/π)(c_κ t)^{−1/2})` ((J1), (J2), Step 3); `p(t) = 1` for `t ≤ T_c` and
`p(t) = (T_c/t)^{1/2}` for `t ≥ T_c`. The zeros outside `𝓛` carry variance
`σ² − Σ_{𝓛} c_i²/2 ≥ σ² − 72`. Therefore

    |Φ(t)| ≤ ρ^{(σ²−72)/8}·p(t)^9,        |Φ^{(j)}(t)| ≤ ρ^{(σ²−80)/8}·p(t)^8   (every j):

for `j ∈ 𝓛` eight factors `p` remain and the other zeros give `ρ^{(σ²−72)/8}`;
for `j ∉ 𝓛` nine factors `p ≤ 1` remain and the other zeros minus `j` have
variance `≥ σ² − 72 − 8`.

Define `W(t) := ρ^{(σ²−72)/8}p(t)^9` (`t ≥ τ`) and `W¹(t) := e^{−(σ²−8)t²/2}`
(`t ≤ τ`), `W¹(t) := ρ^{(σ²−80)/8}p(t)^8` (`t > τ`).

## Step 5. Proposition C₀

Apply Proposition C₀ with `z_1 = z_*` (`< 2.4 < j_{0,1}`), `t_1 = τ`
(`≤ π/2`, and `c_1 t_1 ≤ 4τ = z_1`), `Σ² = σ²`, `Q_4 = 16σ²` (Step 1), `W`, `W¹`:

    −Mδ′ ≥ (1/π)∫_0^τ e^{−σ²t²/2 − 16β_*σ²t⁴} cos t dt − (1/π)∫_τ^∞ W dt − (σ²/(2π))∫_0^∞ t² W¹ dt.

*First integral.* `e^{−x} ≥ 1 − x` and `cos t ≥ 1 − t²/2 > 0` on `[0, τ]` give
`e^{−x}cos t ≥ 1 − x − t²/2` for `x ≥ 0` (if `1 − x ≥ 0` multiply the two lower
bounds; otherwise the left side is `≥ 0 ≥ 1 − x`). With the Gaussian moments
`∫_0^∞ t²e^{−σ²t²/2} = sqrt(π/2)σ^{−3}`, `∫_0^∞ t⁴e^{−σ²t²/2} = 3 sqrt(π/2)σ^{−5}` and
`∫_τ^∞ e^{−σ²t²/2}dt ≤ e^{−σ²τ²/2}/(σ²τ) = ρ^{σ²/8}/(σ²τ)`:

    (1/π)∫_0^τ … ≥ (1/(sqrt(2π)σ))·[1 − 1/(2σ²) − 48β_*/σ²] − ρ^{σ²/8}/(πσ²τ).

*Second integral.* `∫_τ^∞ p^9 ≤ T_c + ∫_{T_c}^∞ (T_c/t)^{9/2}dt = (9/7)T_c`, so
`(1/π)∫_τ^∞ W ≤ (9/(7π))·T_c·ρ^{(σ²−72)/8}`.

*Third integral.* `∫_0^τ t² e^{−(σ²−8)t²/2} ≤ sqrt(π/2)(σ²−8)^{−3/2}` and
`∫_τ^∞ t² p^8 ≤ T_c³/3 + ∫_{T_c}^∞ T_c⁴t^{−2}dt = (4/3)T_c³`, so

    (σ²/(2π))∫ t²W¹ ≤ (1/(sqrt(2π)σ))·(1/2)(1 − 8/σ²)^{−3/2} + (σ²/(2π))·(4/3)T_c³·ρ^{(σ²−80)/8}.

Multiplying by `sqrt(2π)σ` and collecting: `sqrt(2π)·σ·(−Mδ′) ≥ G(σ²)`.

## Step 6. `G(σ²) ≥ G(σ_0²) > 0`

For `x ≥ 15` every subtracted term of `G` is decreasing in `x`: `1/x`,
`(1−8/x)^{−3/2}`, `x^{−1/2}ρ^{x/8}`, `x^{1/2}ρ^{x/8}` (for `x > 4/log(1/ρ) = 4.954…`) and
`x^{3/2}ρ^{x/8}` (for `x > 12/log(1/ρ) = 14.862…`). So `G` is increasing on
`[15, ∞)`, and by Step 2 `G(σ²(M)) ≥ G(236.9711)`. **(E5)**: `G(236.9711) ≥ 0.32197`
(with `β_* ≤ 0.025553`, `T_c ≤ 8.26622`). ∎

## The Arb evaluation (`analytic_leg.py`, `analytic_leg.out`)

For a triple `(M_E, κ, Y)` the script evaluates (E2)–(E5) in ball arithmetic
(128 bits), for both parity families, taking lower endpoints, and checks the
side conditions `M_E + ½ ≥ 4`, `Y + ½ ≥ 4`, `2Y ≤ κM_E`, `σ_0² > 80`.

| `M_E` | `κ` | `Y` | (E2) margin ≥ | `σ_0²` ≥ | `N_χ(κM_E)` ≥ | `T_c` ≤ | `G(σ_0²)` ≥ |
|---|---|---|---|---|---|---|---|
| 300 | 4 | 100 | 0.0016583 | 2979.71089 | 70.3092 | 6.684128 | 0.4974001 |
| 100 | 4 | 40 | 0.0049261 | 774.18188 | 19.1144 | 6.684128 | 0.4899184 |
| 60 | 4 | 48 | 0.0081300 | 403.61089 | 14.9398 | 6.684128 | 0.4804794 |
| 50 | 4 | 40 | 0.0097087 | 318.27738 | 11.3751 | 6.684128 | 0.4750812 |
| 45 | 5 | 40 | 0.0107526 | 277.06752 | 12.9847 | 8.266220 | 0.4681158 |
| **40** | **5** | **40** | **0.0120481** | **236.97110** | **11.3751** | **8.266220** | **0.3219715** |
| 38.5 | 5 | 28 | 0.0124999 | 225.17929 | 9.0645 | 8.266220 | 0.0249902 |

(Every entry is an outward-rounded bound: lower endpoints of the Arb balls
rounded down, upper endpoints rounded up, as printed by `analytic_leg.py` into
`analytic_leg.out`; identical to the corrected table in §6 of `REFEREE-5.md`.
Also printed there: `β(z_*) ≤ 0.0255530` on every line.)

(The proof above is written for the row `M_E = 40`; for another row replace
`40, 200, 5` by `M_E, κM_E, κ` and `Y` accordingly — the argument is unchanged.
`38.5` is the smallest half-integer `M_E` for which some `(κ, Y)` with
`κ ∈ {3,…,7}`, `Y` even in `[16, 78]`, passes with the fixed constants `ρ`,
`z_1 = z_*`: the best on the searched grid, not a proved limit of the method.) `β_*` is evaluated from `arb.bessel_j` at the ball `z_*`; a by-hand
alternative is `β(z) ≤ (1/64)/(1 − z²/5.76)`, i.e. `β_* ≤ 0.0356`, from
`λ_n ≤ λ_2·j_{0,1}^{−2(n−2)}` and `j_{0,1} > 2.4`, which changes `G` in the third
decimal.

## What is and is not claimed

* Claimed: the inequality of the theorem for all real primitive `χ` under GRH
  and `L(1/2,χ) ≠ 0`, all `M ≥ 40`. Uniform in `q`; the worst case of every
  estimate is the smallest conductor of each parity.
* Dependence: (V) and its derivative (paper; (V′) is `lem:Dprop`(i) rewritten),
  Lemma D and Proposition C₀ (`BRIDGE-ATTEMPT.md`), Lemma J (which rests on the
  paper's `lem:phi`, `lem:J0env`, (J2)), and FLINT/Arb for the four numbers.
  None of `lem:Nplus`, `lem:Nlow`, `lem:window`, `lem:negstruct`, `lem:tailexp`,
  `lem:analytic-leg`, `prop:criterion-eval` is used.
* Not claimed: the two-sided asymptotic `θ(m) ∈ [0.90, 1.14]` of
  `thm:monotone-large`; anything below `M = 38.5`; any statement for
  non-real or imprimitive characters.
* Multiplicities: the theorem is about the model in which a zero of multiplicity
  `m` contributes `m` independent amplitudes `c ≤ 4` (the counting lemma, the
  variance identity and Lemma D all count with multiplicity, and nine ordinates
  with multiplicity suffice in Step 3). If a multiple zero were instead modelled
  as one amplitude `mc`, the step `c_jτ ≤ z_*` of Step 4 would fail. Under LI
  the zeros are simple, so the race interpretation is unaffected.
* The range below 40 is never empty (Referee 5): `M_0(χ) < 40` for every real
  primitive `χ`, where `M_0` is the root of `Q_χ(M) = Σ_j 4M³/(M²+γ_j²)^{3/2} = 1`
  up to which the coupling theorem (Theorem B) applies. Indeed `c_j ≤ 4` gives
  `Q = Σc_j³/16 ≥ S_4/64`, and at `M = 40`, by (V′), Step 2 and
  `h′ ≤ ¼ψ′(20.25) + Z_2(40.5) ≤ 0.0127`:
  `S_4 = 16σ² − 64·1600·h′ ≥ 16·236.97 − 102400·0.0127 > 2490`, so `Q(40) > 38`.
  So for a general `χ` there is always an interval `(M_0(χ), 40)` on which
  neither Theorem B nor Theorem E says anything.
* Consistency checks (Arb; not part of the proof): for `q ∈ {3, 4}` and
  `M ∈ {40, 50, 100, 400, 2000}` the bound `G(σ²)/(sqrt(2π)σ)` is `0.50`–`0.88` of
  the certified envelope bound `LB` of `results_q*.txt` at the same `M`, never
  above it (e.g. `q = 4`, `M = 50`: `≈ 0.010220` vs `≈ 0.013341`; true value
  `≈ 0.01368`; these three are midpoints, quoted for orientation, not bounds).
  Referee 5 tested the conclusion numerically for 15 characters and on 408
  abstract amplitude sets without a violation, and found the closed-form
  bounding of Step 5 to lose only `2·10⁻⁴` against the exact value of the three
  integrals at `σ² = σ_0²`. Non-monotone characters (conductors 163, 232; `BRIDGE-ATTEMPT.md`
  §4.3) fail only for `M < 1.2`, far below the range of the theorem.
* Why 40 and not lower: for `τ ≤ t ≤ T_c` the envelope is the constant
  `ρ^{(σ²−80)/8}`, which must beat `σ³T_c³`; and the variance-identity zero count
  is weak (it certifies 12 zeros below 200 for `χ_3` — the bound is `11.375` and
  the count is an integer — where there are over 100).
  Either a better count or a decaying envelope on `[τ, T_c]` would lower the
  threshold; I have not pursued it because the Arb certificate covers
  `[1, 2000]` for `q = 3, 4`.

## Review history

* 2026-11-24: first version (statement, proof, `analytic_leg.py`).
* 2026-11-24, `../REFEREE-5.md` (internal, AI-conducted). Verdict: Theorem E proved as
  written, for every real primitive `χ` of either parity and every `M ≥ 40`, with
  the constant `0.3219`; every step re-derived; script decisions correct. One
  repairable item and minor items concerning this file:
  * **R1** (repaired). `analytic_leg.py` printed lower/upper endpoints with
    round-to-nearest under the labels `>=` / `<=`, so 20 of the 42 numbers in
    `analytic_leg.out`, and the table and Step 6 here, were false in the last
    digit. The script now prints with directed rounding (`p_dn`, `p_up`), the
    output was regenerated, and the table, the displayed (E5) (`G(236.9711) ≥
    0.32197`, `T_c ≤ 8.26622`) and the constant `T_c = 8.2662…` were replaced.
    The thresholds `4/log(1/ρ)`, `12/log(1/ρ)` in Step 6 are now given by
    truncation (`4.954…`, `14.862…`) instead of a rounded value. No decision, and
    not the constant `0.3219`, was affected. (The first half of this repair was
    made by another agent that was stopped mid-task; its edits were checked
    line by line against a fresh run of the script and kept.)
  * Minor 4: "smallest threshold 38.5" qualified as the best on the searched grid.
  * Minor 5: multiplicities — sentence in the statement and an item above.
  * Minor 6: the count is 12, not 11.
  * Minor 12: `Z_1`, `Z_2` need only `s > 2/log 2`; noted where they are introduced.
  * Added fact: `M_0(χ) < 40` always (item above).
