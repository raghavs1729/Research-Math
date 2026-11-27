# Envelope certificate of δ′(M) < 0: method

This note states exactly what `certify_envelope.py` evaluates and why each step
is a valid bound, so that the code can be checked against it line by line.
Mathematical statements (Lemma D, Lemma J, Lemma L, Proposition C₀, Theorem C,
Corollary C′) refer to `monotonicity/BRIDGE-ATTEMPT.md` §1–§2; `lem:…`, `prop:…`
refer to `paper/weighted-races.tex`. Section 8 lists what is *not* certified.

## 0. Statement certified per cell

`χ` = the real primitive character mod `q ∈ {3, 4}`; `M = m + 1/2 > 0`;
`γ_1 < γ_2 < …` the positive ordinates of `L(s,χ)`;

    c_j(M) = 4M/sqrt(M² + γ_j²),   Φ_M(t) = Π_j J_0(c_j(M) t),
    δ(M) = 1/2 + (1/π) ∫_0^∞ Φ_M(t) sin t / t dt,
    σ²(M) = Σ_j c_j²/2,   S_4(M) = Σ_j c_j⁴.

For a cell `[M_a, M_b]` (exact rationals) the program prints `LB` and `PASS`
only if it has verified in ball arithmetic that `LB > 0` and

    −M·δ′(M) ≥ MAIN − TAILF − NUM ≥ LB        for every real M ∈ [M_a, M_b],

hence `δ′(M) ≤ −LB/M_b < 0` on the cell (column `dprime_hi`). Adjacent cells
share their endpoint exactly, so a contiguous run of `PASS` cells proves that
`δ` is strictly decreasing on the union.

## 1. Hypotheses

* **(H1) GRH(χ)** (ordinates real). LI is needed only to read `δ(M)` as the
  logarithmic density of the weighted race.
* **(H2) Variance identity and its `M`-derivative** (`prop:varid`,
  `eq:varidq`, `lem:Dprop`(i)): with `h(M) := (log Λ)′(M + 1/2)`,
  `Λ(s) = (q/π)^{(s+1)/2} Γ((s+1)/2) L(s,χ)`,

      σ²(M) = Σ_j 8M²/(M²+γ_j²) = 4M·h(M),
      S_4(M) = 256 M⁴ Σ_j (M²+γ_j²)^{−2} = 64M·h(M) − 64M²·h′(M),   h′(M) = (log Λ)″(M + 1/2).

  The second line is the first divided by `8M²` and differentiated termwise in
  `M` (`Σ_j (M²+γ_j²)^{−1} = h/(2M)`, so `Σ_j (M²+γ_j²)^{−2} = h/(4M³) − h′/(4M²)`);
  equivalently `S_4/16 = 2σ² − M·D` with `D = (σ²)′` of `lem:Dprop`(i).
* **(H3) Zero lists** `prototype/explore/chi4_zeros.txt` (511 ordinates),
  `chi3_zeros.txt` (537): each decimal within `10⁻²⁰` of an ordinate, zeros
  simple, no other zero with `0 < γ ≤ γ_N`
  (`prototype/explore/certify_arb.py --K 0`). Loaded as balls `decimal ± 10⁻²⁰`
  by `cert/certify_derivative.load_zeros`, which checks strict ordering.
* **(H4)** Correctness of FLINT/Arb (python-flint 0.9.0): ball arithmetic,
  `exp`, `log`, `sqrt`, `cos`, `pi`, `digamma`, Hurwitz `zeta`, `bessel_j` at
  exact points, and (for `M < 3.5`) `dirichlet_char.l`.
* **(H5) Pencil proofs** of Lemma D, Lemma J, Lemma L, Proposition C₀ and
  Corollary C′ (`BRIDGE-ATTEMPT.md`; Lemma D, J, L and Theorem C checked in
  `REFEREE-3.md`, Proposition C₀ and Corollary C′ in `REFEREE-5.md`; both
  reviews internal and AI-conducted).

## 2. The inequality (Corollary C′, case (b), n lowest ordinates; n = 12)

Constants `ρ = 446/1000` (exact), `z_* = 2·sqrt(log(1/ρ))`, `z_far = (8/(πρ))²`,

    b(z) = exp(−z²/4) (0 ≤ z ≤ z_*),   ρ (z_* ≤ z ≤ z_far),   (8/π) z^{−1/2} (z ≥ z_far),
    β(z_1) = (−log J_0(z_1) − z_1²/4)/z_1⁴          (0 < z_1 < j_{0,1}).

`b` is continuous and nonincreasing, and `|J_0| ≤ b` (Lemma J). With

    v := σ²(M_a) − Σ_{i≤n} c_i(M_a)²/2        (= Σ_{i>n} c_i(M_a)²/2 by (H2), (H3)),
    B(t)  := Π_{i≤n} b(c_i(M_a) t) · exp(−(v/2)·min(t, z_*/c_n(M_a))²),
    B¹(t) := B(t)/b(c_1(M_a) t),

and any `z_1 ∈ (0, j_{0,1})`, `0 < t_1 ≤ min(π/2, z_1/c_1(M_b))`:

    MAIN  := (1/π) ∫_0^{t_1} exp(−σ²(M_b)t²/2 − β(z_1)S_4(M_b)t⁴) cos t dt,
    TAILF := (1/π) ∫_{t_1}^∞ B(t) dt,
    NUM   := (S_4(M_b)/(32π)) ∫_0^∞ t² B¹(t) dt,

Corollary C′ gives `−Mδ′(M) ≥ MAIN − TAILF − NUM` for every `M ∈ [M_a, M_b]`.
(The monotonicities used inside the cell — `c_i`, `σ²`, `S_4` increasing in `M`,
`b` nonincreasing — are in the proof of Corollary C′; the program does not
re-derive them, it evaluates the three numbers.)

## 3. What the program computes, with the direction of every rounding

All quantities are Arb balls at 128 bits; `lo(x)`, `hi(x)` denote the exact
lower/upper endpoints of a ball, themselves exact Arb numbers.

**(3.1) Amplitudes.** `c_i(M) = 4/sqrt(1 + (γ_i/M)²)` (each variable occurs
once). Used: `lo(c_i(M_a))` for `i ≤ n` (smaller amplitude ⇒ larger `b`),
`hi(c_n(M_a))` (for the breakpoint `z_*/c_n`), `hi(c_1(M_b))` (for `t_1`).

**(3.2) Closed forms** (`sig2_S4`), at the exact rationals `M_a` and `M_b`:

* `s = M + 1/2 ≥ 4`: `h`, `h′` by `cert/certify_derivative.loglam_derivs_series`
  (Dirichlet series over prime powers `≤ 20000` with integral-test tails,
  `ψ`, `ψ′ = ζ(2,·)` in Arb; see `cert/METHOD.md` (2.3)); then `σ² = 4Mh`,
  `S_4 = 64Mh − 64M²h′`.
* `1 ≤ s < 4` (cells with `M < 3.5`): `h` by the complex step
  `cert/certify_derivative.loglam_d1_complex_step` (384 bits, Cauchy remainder),
  `σ² = 4Mh`; and `S_4` is enclosed **without `h′`**, from the full list:

      Σ_{i≤N} c_i⁴ ≤ S_4 ≤ Σ_{i≤N} c_i⁴ + c_N² · 2V_T,      V_T = σ² − Σ_{i≤N} c_i²/2,

  since by (H3) every unlisted zero has `c ≤ c_N`, so `Σ_{tail} c⁴ ≤ c_N² Σ_{tail} c²`.

  Column `cf=` records which form was used at `M_a`/`M_b`.
* `selftest` checks, at `M ∈ {3.5, 5, 10, 40, 150}` and both moduli, that the
  series value of `h` overlaps the complex-step value and that the series value
  of `S_4` lies in the list enclosure (two evaluations of `S_4` that share
  nothing but the zero list and `σ²`).

Used: `hi(σ²(M_b))`, `hi(S_4(M_b))`, `v_lo := lo(σ²(M_a) − Σ_{i≤n} c_i(M_a)²/2)`
(required `> 0`), `r_lo := lo(lo(z_*)/hi(c_n(M_a)))`.

**(3.3) Grid.** Exact dyadic rationals `0 = t_0 < t_1 < … < t_K = T_end`
(multiples of `2⁻²⁴`; exactness of every Arb conversion is checked): 3000
uniform steps on `[0, ≈12/σ_1]` (`σ_1² ≈ σ² − c_1²/2`, a float used only to place
points), then a geometric grid of ratio `1.01` up to `T_end`. Verified in Arb:
`lo(c_10(M_a))·T_end > hi(z_far)`.

**(3.4) Envelope values (upper).** For a grid point `t` (exact) and `i ≤ n`, let
`z := lo(lo(c_i(M_a))·t)` (an exact lower bound of every argument `c_i(M)t`,
`M` in the cell). Since `b` is nonincreasing, `b(c_i(M)t) ≤ b(z)`, and

    z < lo(z_*)                    ⇒  b(z) = exp(−z²/4);
    z > hi(z_far)                  ⇒  b(z) = (8/π) z^{−1/2};
    hi(z_*) < z ≤ hi(z_far)        ⇒  b(z) ≤ ρ          (b ≤ ρ on [z_*, ∞));
    lo(z_*) ≤ z ≤ hi(z_*)          ⇒  b(z) ≤ exp(−lo(z_*)²/4)   (both candidate values are ≤ this).

All comparisons are between exact numbers, hence decided. The tail factor uses
`m := min(t, r_lo)` (exact) and `v_lo`: `exp(−(v/2)min(t, z_*/c_n)²) ≤ exp(−v_lo m²/2)`.
The program accumulates logarithms in balls and takes `hi(exp(·))`, giving
`B1_hi[k] ≥ B¹(t_k)`, `B_hi[k] ≥ B(t_k)`.

**(3.5) NUM (upper).** `t² B¹(t) ≤ t_{k+1}² B¹(t_k)` on `[t_k, t_{k+1}]` (`B¹`
nonincreasing), so

    ∫_0^{T_end} t² B¹ dt ≤ Σ_k t_{k+1}² · B1_hi[k] · (t_{k+1} − t_k).

Beyond `T_end`: for `t ≥ T_end` and `2 ≤ i ≤ 10`, `c_i(M_a) t ≥ c_10(M_a)T_end > z_far`,
so `b(c_i t) = (8/π)(c_i t)^{−1/2}`; the factors `i = 11..n` are `≤ 1`; the tail
factor is nonincreasing in `t`. Hence

    ∫_{T_end}^∞ t² B¹ dt ≤ (8/π)^9 Π_{i=2}^{10} lo(c_i(M_a))^{−1/2} · exp(−v_lo·min(T_end, r_lo)²/2) · T_end^{−3/2}/(3/2)  =: rem2,
    ∫_{T_end}^∞    B  dt ≤ (same prefactor) · T_end^{−7/2}/(7/2)                                              =: rem0

(`B ≤ B¹`). `NUM := hi( hi(S_4(M_b))/(32π) · (sum + rem2) )`.

**(3.6) TAILF (upper).** `B` nonincreasing: `∫_{t_k}^{t_{k+1}} B ≤ B_hi[k](t_{k+1} − t_k)`.
Suffix sums plus `rem0` give an upper bound of `∫_{t_k}^∞ B` for every grid
index `k`.

**(3.7) MAIN (lower).** For each candidate `z_1 ∈ {0.40, 0.50, …, 2.30}` (exact
rationals, all `< 2.4 < j_{0,1}`, `lem:phi`(iii)) the program takes `t_1` := the
largest grid point with (float guess, then **verified in Arb**)
`t_1 < lo(π/2)` and `hi(hi(c_1(M_b))·t_1) ≤ z_1`. `β(z_1)` is evaluated from
`arb.bessel_j(0)` at the exact point `z_1` (required: `lo(J_0(z_1)) > 0`), upper
endpoint taken. The integrand

    G(t) = exp(−hi(σ²(M_b))t²/2 − hi(β)hi(S_4(M_b))t⁴)·cos t

is positive and decreasing on `[0, π/2)`, so `∫_0^{t_1} G ≥ Σ_{k=1}^{k_1} lo(G(t_k))(t_k − t_{k−1})`
(lower Riemann sum; the sum may be stopped early — dropping positive terms
only lowers it). `MAIN := lo(sum/π)`.

**(3.8) Decision.** For each candidate, `val := lo(MAIN − hi(suffix[k_1]/π))`;
the candidate with the largest `val` is kept (all candidates are valid, so
choosing among them is free). `LB := lo(val − NUM)`. **`PASS` iff `LB > 0`**
(a certified Arb comparison). `K := hi(NUM/val)` is printed as a margin
indicator only; the decision never uses it.

## 4. Fail-closed behaviour

"Fail-closed" means: every check listed below is a certified Arb comparison,
and a cell that cannot be decided is not `PASS`. It does **not** mean that a
wrong zero file would be detected: that is hypothesis (H3). In series mode a
cell compares the zero list with the closed forms only through `v > 0`;
`selftest` compares `S_4`(series) with the list enclosure at five values of `M`
and so detects gross errors only (Referee 5: an ordinate shifted by `10⁻³`, a
deleted zero or `γ_1 := 0.2029` is refused by `selftest`; a shift of `10⁻⁶` is
below every tolerance).

The program raises (the cell is written as `ERROR … FAIL`) if: a ball is not
finite; `v` is not certified positive; a grid point is not exactly
representable; `T_end` is not certified large enough; `t_1` is not certified
admissible; `J_0(z_1)` or the `MAIN` integrand is not certified positive; the
closed forms are not certified positive. Floats occur only (i) in the placement
of grid points, `T_end` and the guess of `t_1` — each then an exact rational
verified in Arb — (ii) in the early stop of the `MAIN` sum, and (iii) in
printing. Printed `LB`, `MAIN` are rounded down, `TAILF`, `NUM`, `K`,
`dprime_hi` rounded up, from the Arb endpoints.

## 5. Output columns

    M_a M_b LB= MAIN= TAILF= NUM= K= dprime_hi= sig2b= S4b= v= z1= t1= Tend= rem= npts= n= cf= secs= PASS

`LB ≤ −Mδ′(M)` on the cell; `dprime_hi = −LB/M_b ≥ δ′(M)`; `sig2b`, `S4b` upper
bounds of `σ²(M_b)`, `S_4(M_b)`; `v` lower bound of `v`; `rem` = `rem2`; `npts`
grid size; `cf` closed-form mode at `M_a`/`M_b`. `# SUMMARY` lines give the
union of certified cells and the worst values.

The second header line, `# sha256 prover=… shared=… zeros=…`, records the
SHA-256 of `certify_envelope.py`, of `cert/certify_derivative.py` and of the
zero file used. `range` resumes onto an existing results file only if that
line equals the current hashes; otherwise it stops and asks for the file to be
deleted (so a stale `PASS` line produced by different code or data cannot be
kept).

## 6. Results (python-flint 0.9.0, one niced process; rerun of 2026-11-24 after `REFEREE-5.md`)

| file | certified union | cells | FAIL/ERROR | worst `K` (upper bound) | smallest `LB` | weakest bound on `δ′` | wall time |
|---|---|---|---|---|---|---|---|
| `results_q4.txt` | `[1, 2000]` | 765 | 0 | 0.44115 on `[1973.3722, 1993.1059]` | 1.318804e−3 | `δ′ ≤ −6.616829e−7` | 75 s (first run: 118 s) |
| `results_q3.txt` | `[1, 2000]` | 765 | 0 | 0.43819 on `[1973.3722, 1993.1059]` | 1.353296e−3 | `δ′ ≤ −6.789885e−7` | 76 s (first run: 103 s) |

Sub-ranges, counting every cell that meets the range (so including the cell that
straddles its right endpoint): `K ≤ 0.3614` on the cells covering `[3.6, 20]`
(`q = 4`, attained on `[19.8995, 20.0984]`), `≤ 0.3364` on those covering
`[4.4, 20]` (`q = 3`); `K ≤ 0.4219` resp. `0.4169` on those covering `M ≤ 400`
(attained on `[397.6236, 401.5998]`). "Worst" sits at the top of the range
because `K` increases slowly towards `1/2` and `LB ∼ M^{−1/2}`. On the 127 cells
with `M_a < 3.5` the closed forms are in `list` mode at `M_a` (126 `list/list`,
1 `list/series`).

## 7. Consistency with the direct certificate (`consistency.py`; not part of the proof)

`monotonicity/cert/results_q*.txt` and `results_q*_stretch.txt` certify
`lo ≤ δ′(M) ≤ hi < 0` on cells covering `[3.6, 320]` (`q = 4`) and `[4.4, 350]`
(`q = 3`), by an interval Riemann sum of the oscillatory integrand itself. If an
envelope cell and a direct cell intersect in `[x, y]`, then for `M = x`:
`LB ≤ x|δ′(x)| ≤ x·|lo|`. Checked for every intersecting pair:

| modulus | pairs | violations | max `LB/(x·lo)` (must be ≤ 1) | range of the enclosures `[LB/(y·lo), LB/(x·hi)]` of `LB/(M·δ′)` |
|---|---|---|---|---|
| 4 | 1061 | 0 | 0.9453 | `[0.8307, 1.0239]` |
| 3 | 1028 | 0 | 0.9450 | `[0.8798, 1.0237]` |

(absolute values of `lo`, `hi` understood). The upper end exceeds 1 only through
the width of the direct enclosures (`hi/lo ≥ 0.92`). `consistency_q*.txt` list
about every 26th pair. This comparison uses the printed decimals as floats; it
is a cross-check of two certificates, not a third one.

**Scope of this comparison** (Referee 5, mutation tests). On the common range
`TAILF ≤ 1.9·10⁻³` and the far remainder `≤ 3·10⁻⁵` are far below the 4–8 % width
of the direct enclosures, so the comparison validates `MAIN − NUM` only: a
prover with `TAILF` dropped, or with the far remainder dropped, passes it.
`TAILF` and the remainder are covered by code reading and by the referee's
independent mpmath evaluation of the exact cell quantity
`E = MAIN_cell − TAILF_cell − NUM_cell` at 13 cells per modulus, including the
cells near `M = 1` where `TAILF` is 43 % (`q = 4`) resp. 62 % (`q = 3`) of `MAIN`:
`LB ≤ E` with
`LB/E ∈ [0.983, 0.997]` (`../REFEREE-5.md` §3.3).

What the two certificates share: (H1)–(H4), the zero files, and the module
`cert/certify_derivative.py`, from which this prover imports `load_zeros`,
`loglam_derivs_series` (`h`, `h′`), `loglam_d1_complex_step`, `make_grid` and the
exact-rational and printing helpers. What they do not share: the integrand, the quadrature, the
tail model, the use of Bessel functions on wide balls (the envelope certificate
evaluates `J_0` only at the 20 exact points `z_1`), and the analytic input
(Corollary C′ here; differentiation under the integral with an explicit
dominant there).

## 8. Caveats — what is not certified

1. Everything is conditional on (H1)–(H5). (H5) is a pencil proof: the program
   certifies the *hypothesis* of Corollary C′, not the corollary.
2. (H3) is the output of a separate program and is not re-proved here. Only the
   first 12 ordinates enter the envelope; completeness of the list is used for
   `v` (case (b)) and, in `list` mode, for the enclosure of `S_4`.
3. (H2) is a theorem of the paper; what is machine-checked is the evaluation of
   its right-hand side. The formula for `S_4` is derived in (H2) above and
   cross-checked against the list in `selftest`.
4. Soundness, not sharpness. **On the range common with the direct certificate**
   (`[3.6, 320]` resp. `[4.4, 350]`) the bound loses 4–17 % of `|δ′|` (Section 7).
   Near `M = 1` it loses more: against a floating-point evaluation of `−Mδ′`
   (Referee 5), `LB/truth = 0.625, 0.690, 0.766, 0.842, 0.899` at
   `M = 1, 1.1, 1.25, 1.5, 2` (`q = 4`) and `0.508, 0.611, 0.717, 0.821, 0.917`
   (`q = 3`), i.e. up to 38 % resp. 49 % at `M = 1`. The bound is positive by a
   wide margin there (`K ≤ 0.08`); what grows is `TAILF/MAIN`.
5. The range stops at `M = 2000` by choice (cost is linear in the number of
   cells; nothing in the method degrades). Beyond, and already from `M = 40`,
   Theorem E (`ANALYTIC-LEG.md`) applies. Below `M = 1` the criterion itself
   weakens (`TAILF` approaches `MAIN`); there Theorem B applies
   (`M_0 > 3.68` resp. `4.47`, `cert/aux_results.txt`). For the chain only
   `M_0 > 1` is needed, and it follows from this certificate's own numbers
   (Referee 5): `Q = Σc_j³/16 ≤ c_1σ²/8` is increasing in `M`, and at `M = 1.01`,
   `c_1 ≤ 4.04/γ_1 ≤ 0.671` and `σ² ≤ 0.6292` (`sig2b` of the first cell) for
   `q = 4`, `c_1 ≤ 0.503` and `σ² ≤ 0.4599` for `q = 3`, so `Q(1.01) ≤ 0.053` resp.
   `0.029`, hence `M_0 > 1.01`.
6. `K` and `dprime_hi` are derived quantities printed for convenience.
7. The certificate concerns `q = 3, 4` only.

## 9. Reproduction

    sh run_all.sh
    python3 certify_envelope.py cell --q 4 --Ma 3.6 --Mb 3.636

About 3 minutes of CPU in total; deterministic (identical output on every
run, apart from `secs` and the timestamp).

## 10. Review history

* 2026-11-24: prover, results and this note written.
* 2026-11-24, `../REFEREE-5.md` (internal, AI-conducted). Audit of the prover line
  by line, a full fresh rerun (identical), an independent mpmath evaluation of
  the cell quantity at 26 cells, 19 mutants. Verdict: the certificate is sound;
  it evaluates exactly the hypothesis of Corollary C′, case (b), `n = 12`, with
  Theorem C's own `NUM`; every rounding goes the right way; no float on a
  certified path; every refusal path fails closed. Items concerning this
  directory, and what was done:
  * Minor 1 (sub-range maxima of `K` omitted the straddling cell): Section 6
    corrected (`0.3614` / `0.3364` up to 20, `0.4219` / `0.4169` up to 400).
  * Minor 2 ("loses 4–17 %" without its range): Section 8.4 now states the range
    and the loss near `M = 1`.
  * Minor 3 (`suf[len−1]` left at 0 in `certify_cell`): initialised to `rem0`.
    The index is unreachable (`t_1 < π/2 < T_end`), so no certified value could
    change, and none did.
  * Minor 7 (`consistency.py` does not exercise `TAILF`): stated in Section 7
    and in the docstring of `consistency.py`.
  * Minor 8 (resume logic accepted any existing `PASS` line): results files now
    carry the SHA-256 of the prover, the shared module and the zero file, and
    `range` refuses to resume onto a file with other hashes (Section 5).
  * Minor 9 ("fail-closed" and a wrong zero file): Section 4 and the docstring
    of the prover.
  * Minor 10 (what is shared with `../cert/`): Section 7 lists the imported
    functions.
  * Added fact: only `M_0 > 1` is needed, and it follows from the certificate
    (Section 8.5).
* **Rerun after these code changes (2026-11-24).** `results_q4.txt` and
  `results_q3.txt` were deleted and recomputed with the modified prover
  (SHA-256 `cf52de54…e1cd35`): 765 + 765 cells, all `PASS`. Compared with the
  previous files after removing the `secs=` field and the `#` lines: **all 1530
  cell lines identical** (every printed digit of `LB`, `MAIN`, `TAILF`, `NUM`,
  `K`, `dprime_hi`, `sig2b`, `S4b`, `v`, `z1`, `t1`, `Tend`, `rem`). The only
  differences are the new `# sha256` header line, the run times and the
  timestamp of the `# SUMMARY` lines. `points_q*.txt` (apart from `secs`),
  `consistency_q*.txt`, `selftest.out` and `analytic_leg.out` were regenerated
  and are byte-identical to the previous ones (`analytic_leg.out`: identical to
  the R1-corrected output).
