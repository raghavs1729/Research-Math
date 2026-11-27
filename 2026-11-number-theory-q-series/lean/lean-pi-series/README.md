# lean-pi-series: Ramanujan-type series for 1/π, machine-checked in Lean 4

Status (2026-11-24). **Five series are proved completely** (no `sorry`, no `axiom`, no `native_decide`;
`#print axioms` gives `[propext, Classical.choice, Quot.sound]` for each). Two of them were open conjectures of
Z.-W. Sun, proved in this repository: the `32√2/π` series and the `64√2/π²` series (a Ramanujan-type series for
`1/π²`, whose proof needs a Laplace-method asymptotic; section 6).

## 1. What is proved

Each statement is in its own Mathlib-only file `PiSeriesStatement/*.lean` (with a prose docstring giving the
formula as its author wrote it), and is proved verbatim in `PiSeriesSolution.lean`.

```lean
-- Ramanujan 1914, eq. (28)                                   PiSeriesStatement/Ramanujan4.lean
theorem ramanujan_series_four_div_pi :
    HasSum (fun n : ℕ => (Nat.choose (2 * n) n : ℝ) ^ 3 * (6 * n + 1) / 2 ^ (8 * n))
      (4 / Real.pi)

-- Ramanujan 1914, eq. (29)                                   PiSeriesStatement/Ramanujan16.lean
theorem ramanujan_series_sixteen_div_pi :
    HasSum (fun n : ℕ => (Nat.choose (2 * n) n : ℝ) ^ 3 * (42 * n + 5) / 2 ^ (12 * n))
      (16 / Real.pi)

-- Z.-W. Sun, arXiv:2307.03086, Conjecture 5.25(i)  (NEW)      PiSeriesStatement/SunSqrt2.lean
theorem sun_series_thirty_two_sqrt_two_div_pi :
    HasSum (fun k : ℕ => (360 * (k : ℝ) ^ 3 + 612 * (k : ℝ) ^ 2 + 230 * k + 15)
        * (Nat.choose (3 * k) k : ℝ) * (Nat.choose (6 * k) (3 * k) : ℝ) ^ 2
        / ((3 * (k : ℝ) + 1) * (3 * k + 2) * 2 ^ (15 * k)))
      (32 * Real.sqrt 2 / Real.pi)

-- Z.-W. Sun, arXiv:2307.03086, Conjecture 5.33(i)  (NEW)      PiSeriesStatement/SunPiSq.lean
theorem sun_series_sixty_four_sqrt_two_div_pi_sq :
    HasSum (fun k : ℕ => (1344 * (k : ℝ) ^ 3 + 944 * (k : ℝ) ^ 2 + 156 * k + 9)
        * (Nat.choose (2 * k) k : ℝ) ^ 2 * (Nat.choose (4 * k) (2 * k) : ℝ) ^ 3
        / ((2 * (k : ℝ) + 1) * 2 ^ (22 * k)))
      (64 * Real.sqrt 2 / Real.pi ^ 2)

-- Bauer 1859 (conditionally convergent: limit of partial sums) PiSeriesStatement/Bauer.lean
theorem bauer_series_two_div_pi :
    Filter.Tendsto (fun N : ℕ => ∑ n ∈ Finset.range N,
        (-1 : ℝ) ^ n * (Nat.choose (2 * n) n : ℝ) ^ 3 * (4 * n + 1) / 2 ^ (6 * n))
      Filter.atTop (nhds (2 / Real.pi))
```

In words:

* `∑_{n≥0} C(2n,n)³ (6n+1) / 2^{8n} = 4/π` and `∑_{n≥0} C(2n,n)³ (42n+5) / 2^{12n} = 16/π`
  (S. Ramanujan, *Modular equations and approximations to π*, Quart. J. Math. 45 (1914));
* `∑_{k≥0} (360k³+612k²+230k+15) C(3k,k) C(6k,3k)² / ((3k+1)(3k+2) 2^{15k}) = 32√2/π`
  (conjectured by Z.-W. Sun; the proof is Theorem 4 of `../ramanujan-series/RESULTS.md`, found in this repository);
* `∑_{k≥0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^{22k}) = 64√2/π²`
  (conjectured by Z.-W. Sun; the proof is Theorem 11 of `../ramanujan-series/RESULTS.md`, reorganised in `N2-PLAN.md`);
* `∑_{n<N} (-1)^n C(2n,n)³ (4n+1) / 2^{6n} → 2/π` (G. Bauer 1859). This series converges only conditionally
  (its terms have size `≈ n^{-1/2}`), so Mathlib's `HasSum`, which means unconditional summability, would be
  false; the statement is the convergence of the partial sums.

`HasSum f a` for the four `HasSum` statements: the family `f : ℕ → ℝ` is unconditionally summable with sum `a`; all terms
are positive, so this is ordinary convergence of the series.

Numerical sanity check of the values (60 digits, `mpmath`): all differences `< 10^{-59}`.

Remark. The `32√2/π` series was also formalized, independently and in parallel, in
`../ramanujan-series/lean/RamanujanSeries/Sun1pi.lean` (`SunSeries1pi.sun_series_1pi`, same statement up to notation).

## 2. Method

Four of the five proofs (all but the `64√2/π²` series, section 6) are Wilf–Zeilberger proofs of one shape (the shape of Theorem 4 in `RESULTS.md`):
a hypergeometric term `h(n,k)` given by its two ratios `σ = h(n,k+1)/h(n,k)`, `ρ = h(n+1,k)/h(n,k)`
(explicit rational functions) and a rational certificate `Q(n,k)` with

    (n+1) h(n+1,k) − n h(n,k) = Q(n,k+1) h(n,k+1) − Q(n,k) h(n,k),            F = n·h,  G = Q·h.

`G(n,0)` is the summand of the series and `F(0,k) = 0`. Summing over the finite rectangle `n < N`, `k < K`
and letting `K → ∞`, then `N → ∞`, gives `∑_n G(n,0) = lim_K G(0,K)`: the whole value of the series is the
limit of **one** boundary hypergeometric term,

| series | `h(0,K)` | boundary limit | Mathlib input |
|---|---|---|---|
| Ramanujan `4/π`, `16/π`, Bauer | `(C(2K,K)/4^K)²` | `(2K+1) h(0,K) → 2/π` | `Real.Wallis.tendsto_W_nhds_pi_div_two` |
| Sun `32√2/π` | `(1/4)_K (3/4)_K / K!²` | `(4K+1) h(0,K) → 2√2/π = 4 sin(π/4)/π` | `Real.tendsto_euler_sin_prod` at `x = 1/4` |

Only integers `n, k` occur. **Carlson's theorem, analytic continuation in `k`, and limits of parameter-dependent
infinite sums are not used** (the published WZ proofs of the Ramanujan series use Carlson's theorem; see `PRIOR-ART.md`).

The pairs:

| series | `σ(n,k)` | `ρ(n,k)` | `Q(n,k)` |
|---|---|---|---|
| `4/π` | `(2k+1)²/(4(k+n+1)²)` | `(2n+1)³/(32(n+1)(k+n+1)²)` | `(4k+6n+1)/8` |
| `16/π` | `(2k+2n+1)²/(4(k+2n+1)²)` | `(2n+1)³(2k+2n+1)²/(128(n+1)(k+2n+1)²(k+2n+2)²)` | `(16k³+(104n+36)k²+(224n²+160n+24)k+168n³+188n²+62n+5)/(32(k+2n+1)²)` |
| `32√2/π` | `(4k+6n+1)(4k+6n+3)/(16(k+n+1)(k+3n+1))` | `3(6n+1)(6n+5)(4k+6n+1)(4k+6n+3)(4k+6n+5)/(1024(n+1)(k+n+1)(k+3n+1)(k+3n+2)(k+3n+3))` | `(64k³+(160n+144)k²+(240n²+432n+92)k+360n³+612n²+230n+15)/(64(k+3n+1)(k+3n+2))` |
| Bauer | `(2k+2n+1)(2k+1)/(4(k+n+1)²)` | `−(2n+1)²(2k+2n+1)/(8(n+1)(k+n+1)²)` | `(4n+2k+1)/2` |

The first is the family printed (without certificate) in Guillera's *Generators of some Ramanujan formulas* (2006),
section "Some more series"; the second and fourth were recomputed with the repository's engine
(`../ramanujan-series/code/wz.py`, exploration in `derivation/`); the third is Theorem 4 of `RESULTS.md`.

## 3. Layout

| file | content |
|---|---|
| `PiSeries/WZ.lean` | finite WZ identity `wz_finite_identity`; the limit argument `hasSum_of_wz` |
| `PiSeries/Term.lean` | `hterm σ ρ` (term defined by its ratios); `hasSum_of_ratios`: reduces a series to two rational identities, three rational inequalities and one boundary limit |
| `PiSeries/WZSigned.lean` | signed variant `tendsto_of_wz_signed` (for Bauer's alternating series) |
| `PiSeries/Wallis.lean` | `(2K+1)(C(2K,K)/4^K)² → 2/π`; recurrence of `C(2n,n)` |
| `PiSeries/Ramanujan6.lean`, `Ramanujan42.lean`, `SunSqrt2.lean`, `Bauer.lean` | the four proofs |
| `PiSeries/SunPi2/*.lean` | Sun's `64√2/π²` series (section 6) |
| `PiSeriesStatement/*.lean` | trusted statements: Mathlib only, one theorem per file, proof `sorry` (deliberate) |
| `PiSeriesSolution.lean` | the same statements, proved from `PiSeries/` |
| `scripts/PrintAxioms.lean`, `scripts/compare_statements.sh` | axiom audit; comparison of Statement and Solution (text and elaborated types) |
| `derivation/` | scripts that found and checked the WZ pairs (not needed to check the proofs) |
| `PRIOR-ART.md` | what was already formalized elsewhere (nothing for these series) |
| `N2-PLAN.md` | the reorganised proof of the `64√2/π²` series and its work packages |

## 4. Pins and how to check

Lean `leanprover/lean4:v4.33.0-rc2`, Mathlib `51e6992efd06126df61a496bebf8f49482a4e129` (the pins of
`../ramanujan-series/lean`).

    curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
    export PATH="$HOME/.elan/bin:$PATH"
    cd lean/lean-pi-series
    lake exe cache get          # Mathlib binaries; never build Mathlib from source
    lake build                  # the warnings "declaration uses sorry" come from PiSeriesStatement/*.lean only (deliberate)
    lake env lean scripts/PrintAxioms.lean
    ./scripts/compare_statements.sh

Output of `lake env lean scripts/PrintAxioms.lean`:

    'ramanujan_series_four_div_pi' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ramanujan_series_sixteen_div_pi' depends on axioms: [propext, Classical.choice, Quot.sound]
    'sun_series_thirty_two_sqrt_two_div_pi' depends on axioms: [propext, Classical.choice, Quot.sound]
    'bauer_series_two_div_pi' depends on axioms: [propext, Classical.choice, Quot.sound]
    'sun_series_sixty_four_sqrt_two_div_pi_sq' depends on axioms: [propext, Classical.choice, Quot.sound]

`scripts/compare_statements.sh` checks that every statement file imports only Mathlib and defines nothing, that each
statement (docstring included) occurs byte for byte in `PiSeriesSolution.lean`, that the solution side contains no
`sorry` / `axiom` / `native_decide`, and that `#check @name` prints the same type from both sides. It is a weak
substitute for `leanprover/comparator`, which is not installed here.

## 5. What is and is not formalized

Formalized, for each of the four series of section 2: the WZ equation (as an identity between rational functions, by
`field_simp; ring`), the compatibility of the two ratios, the finite telescoping identity for all `N, K`, positivity
and decay bounds, summability, both limits, the boundary limit (from Wallis' product, resp. Euler's sine product, both
in Mathlib), and the identification of `G(n,0)` with the binomial summand.

Not formalized: how the pairs were found (Gosper's algorithm / the seed search); this is not needed to check the proofs.

## 6. Sun's `64√2/π²` series (Conjecture 5.33(i)) — complete

    ∑_{k≥0} (1344k³+944k²+156k+9) C(2k,k)² C(4k,2k)³ / ((2k+1) 2^{22k}) = 64√2/π²

This is Theorem 11 of `RESULTS.md`. Its WZ pair `(a, R·a)` lives on the triangle `0 ≤ k < n`, and the value of the
series is not a single boundary term but the limit of the row sums `T_N = ∑_{k<N} a(N,k)`, a Laplace-type asymptotic.
Everything is formalized (`PiSeries/SunPi2/`, about 2700 lines), following the reorganised proof of `N2-PLAN.md`:

| file | content |
|---|---|
| `Defs.lean` | `σ`, `ρ` (the two ratios of `a`), `a`, `T`, the summand `t`, `x* = √3/2`, the closed form `cf` |
| `PartA.lean` | exact part: `a > 0` on the triangle, `a(n+1,k) = ρ a(n,k)`, the WZ equation with the degree-6 certificate, the boundary row `k = n−1`, `T(N+1) = T(N) + t(N)/256`, hence `∑_{n<N} t(n) = 256 T_N`; summability of `t`; the factorial closed form of `a(N,k)` for `k < N ≤ 2k+2` |
| `Tails.lean` | `ρ(n,k) ≤ q < 1` for `k ≤ (x*−η)n` and `≥ q' > 1` for `k ≥ (x*+η)n` (large `n`), hence `a(N,k) ≤ B qᴺ` for `|k − x*N| ≥ δN`, using only positivity and `a ≤ B := ∑t/256` |
| `Ratio.lean` | `|log σ(N,k) + 24(k−x*N)/N| ≤ ε|k−x*N|/N + C/N` near `x*N` (from `HasDerivAt (log Λ) (−24) x*`), and the resulting Gaussian profile `∏σ = exp(−(12±ε)j²/N ± C/N)` |
| `Gaussian.lean` | `(1/√N) ∑_{|j|≤δN} exp(−c j²/N) → √(π/c)` (sum–integral comparison, `integral_gaussian`) |
| `Stirling.lean`, `Centre.lean` | `log m! = cN log N + cN(log c − 1) + u log N + u log c + ½ log(2πcN) + o(1)` for `m = cN + u`; `a(N,⌊x*N⌋)·√N → √3/(√2 π^{5/2})` (eleven factorials) |
| `Assembly.lean` | squeeze: these four inputs give `T_N → √3/(√2π^{5/2})·√(π/12)` |
| `Main.lean` | `T_N → √2/(4π²)` and the theorem |

Differences from the written proof in `RESULTS.md` (Part B): the tails are controlled through the ratio in `n` and
the a-priori bound `a(n,k) ≤ T_n ≤ ∑t/256` instead of log-concavity in `k` (no case analysis for `k < N/2`, no odd-`N`
junction); the window is `|k − x*N| ≤ δN` with `ε`-`δ` bookkeeping instead of `N^{3/5}`; near the peak all eleven
Gamma values in the closed form are factorials of integers, so Mathlib's Stirling formula for `n!` suffices.
The work packages were proved by parallel agents against fixed statements; the only trust needed is in Lean's kernel
(`#print axioms` of the final theorem lists the three standard axioms).
