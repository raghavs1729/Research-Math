# Independent audit of `monotonicity/cert/` — certified negativity of δ′(M)

Auditor's working directory: `monotonicity/cert-audit/` (nothing was written
anywhere else in the repository; every interpreter that touched the audited tree was run with
`python3 -B`).  Date of the audit: 2026-11-24.  python-flint 0.9.0,
mpmath 1.4.1, Python 3.13.

## 0. Verdict

* **No RIGOR BREAK was found.**  Every inequality of `METHOD.md` was re-derived
  and is correct, is used in the correct direction and inside its range of
  validity, and is implemented as stated by `certify_derivative.py`.
* **Three GAPs** (section 3): G1 the differentiation-under-the-integral lemma
  that the certificate needs is *not* the paper's `lem:diffint` (the needed
  argument is supplied in full below); G2 the certificate is bound to its
  hypothesis (H3) — the zero lists — only by convention, and the checker
  silently accepts a wrong list (closed here by `SHA256SUMS` + `replay.sh`);
  G3 the result files were produced by earlier versions of the code and the
  shipped `run_all.sh` does not recompute them (closed here by a from-scratch
  replay that reproduces all 613 lines field for field).
* **Eleven MINOR** items (section 3.3).
* The claim was **re-proved by a second, independently written rigorous
  program** (`indep_enclose.py`: different quadrature, different treatment of
  the zeros, different closed-form evaluation) on **every one of the 317 + 296
  audited cells**: every independent enclosure is negative and lies strictly
  inside the audited one.

The exact statement that is certified, and its hypotheses, is in section 7.

## 1. What was audited, and when

The audited directory changed during the audit (the author was still
working).  State at the start (05:45 UTC): `certify_derivative.py` of commit
`ff715a2` (sha256 `c18eb4aa…`), `run_all.log` and `run_stage2.log` present.
At 05:58–05:59 the author changed `certify_derivative.py` (one line: the text
"section 6" → "section 8" in the header written to `aux_results.txt`),
`METHOD.md` (§4 wording, §9 items 8–10, new §11 "Results") and
`aux_results.txt` (header line), added `crosscheck_q{3,4}.txt` and
`results_q3_stretch.txt`, and **removed `run_all.log` and `run_stage2.log`**
(commit `40edd9d`, 06:00).  Everything below refers to the final state, pinned
in `SHA256SUMS`:

| file | sha256 (first 16 hex) |
|---|---|
| `monotonicity/cert/certify_derivative.py` | `f4674b09e1825744` |
| `monotonicity/cert/METHOD.md` | `09c149122d515936` |
| `monotonicity/cert/results_q4.txt` | `17310a535d769973` |
| `monotonicity/cert/results_q3.txt` | `5d61ee7cd2ceae18` |
| `monotonicity/cert/aux_results.txt` | `bd7d1a0ea8d226a1` |
| `monotonicity/cert/crosscheck.py` | `6ad14160cb0183b3` |
| `prototype/explore/chi4_zeros.txt` | `0e73e00a6c1609bd` |
| `prototype/explore/chi3_zeros.txt` | `e3075235a847f43e` |

`run_all.log` (listed in the audit request) no longer exists and cannot be
hashed.  As read at 05:45 it contained the `selftest` output, an **older**
`aux` block (produced with the cruder tail bound `(M/γ_N)·V_T/2` for `Q`:
`Q(3.68) ∈ [0.993839, 0.994382]`, `3.688482 < M_0(χ₄) < 3.689307`, …) and the
two `# SUMMARY` blocks of the priority ranges.  The older aux intervals contain
the final ones of `aux_results.txt`, so the two are consistent; the log was
stale with respect to the code but not wrong.

A `__pycache__/certify_derivative.cpython-313.pyc` appeared in the audited
directory at 06:25:58; it was not created by the audit (all audit runs use
`-B`; verified on a scratch copy that `python3 -B certify_derivative.py …`
writes no bytecode, and the multiprocessing start method is `fork`).

## 2. Part 1 — the mathematics of METHOD.md

Notation as in METHOD.md.  "Code" gives the place in
`certify_derivative.py` (line numbers of the final version).

### 2.1 Summary table

| # | statement of METHOD.md | re-derivation | direction / range | code implements it? |
|---|---|---|---|---|
| 1 | §0: `c′ = 4γ²/(M²+γ²)^{3/2}`, `δ = 1/2 + (1/π)∫ sinc·Φ` | ✔ (Gil-Pelaez at x = −1 for the symmetric variable S = Σ c cos θ) | — | ✔ l.439–440 (`c = 4/sqrt(1+(γ/M)²)`, `c′ = (4/γ)(1+(M/γ)²)^{-3/2}`, algebraically identical, M occurs once) |
| 2 | (H2) `σ² = Σ 8M²/(M²+γ²) = 4M(log Λ)′(M+½)`, `D = dσ²/dM = Σ 16Mγ²/(M²+γ²)² = 4(log Λ)′ + 4M(log Λ)″` | ✔ from `Λ′/Λ(s) = Σ_ρ Re 1/(s−ρ)` (Hadamard, real χ, `Re B = −Σ Re 1/ρ`); needs **every** non-trivial zero on the line and `L(½,χ) ≠ 0` | — | ✔ `tail_VD`, `head_sums` l.330–356 |
| 3 | (2.1) exact tail sums `V_T = σ² − Σ_H`, `D_T = D − Σ_H` | ✔ — exact **given** #2 and (H3) | — | ✔ l.350–351; both required `> 0` |
| 4 | (2.2) `V_T`, `D_T` nondecreasing in M on the cell | ✔ termwise: `∂_M 8M²/(M²+γ²) = 16Mγ²/(…)² > 0`; `∂_M 16Mγ²/(M²+γ²)² = 16γ²(γ²−3M²)/(M²+γ²)³ > 0` iff `γ² > 3M²` (sum of nondecreasing functions; no termwise differentiation of a series is needed) | needs `γ > √3 M_b` for every tail zero | ✔ run-time check l.423 `γ_N.lower()² > 3·M_b.upper()²`; hulls l.429–430 |
| 5 | (2.3) closed forms and Dirichlet-series tails | ✔ `∫_N^∞ log x·x^{-s}dx = N^{1−s}(L/a+1/a²)`, `∫_N^∞ log²x·x^{-s}dx = N^{1−s}(L²/a+2L/a²+2/a³)`; integrands decreasing for `x ≥ e^{2/s}` | `s ≥ 3` enforced l.254 | ✔ l.259–284; `arb(2).zeta(x)` is the Hurwitz `ζ(2,x) = ψ′(x)` (checked) |
| 6 | (B1) `0 ≤ −log J₀(z) − z²/4 ≤ C₄z⁴` on `[0,z₀]` | ✔ positive Taylor coefficients `λ_n/n` (re-computed as exact rationals: `λ₁..λ₄ = 1/4, 1/32, 1/192, 11/12288`) | only for `z₀ < j₀,₁`; used with `z₀ = 1` | ✔ l.459 |
| 7 | (B2) `z/2 ≤ J₁/J₀(z) ≤ z/2 + κz³` on `[0,z₀]` | ✔ same argument | same | ✔ l.460 |
| 8 | (B3) `\|J_n\| ≤ 1`, `\|J₀″\| ≤ 1`, `\|J₁″\| ≤ 1`; (B4) `\|J₁(z)\| ≤ z/2`; (B5) `\|J₀(x)\| ≤ x^{-1/2}` | ✔ (Bessel's integral; `J₀″ = (J₂−J₀)/2`, `J₁″ = (J₃−3J₁)/4`; Watson 13.74 gives even `√(2/(πx))`) | all valid for all real `x` (`x > 0` for B5) | ✔ `bessel01` l.365–397 |
| 9 | (5.0) tail arguments `c_γ(M)t ≤ z₀` on zone 1 | ✔ `c_γ(M)` increasing in M, decreasing in γ | this is the "below the first zero of J₀" condition | ✔ **checked in Arb** l.474 for the `T₁` actually used |
| 10 | (5.1) `G e^{−C₄S₄t⁴} ≤ Ψ ≤ G`, `S₄ = 2c_Γ²V_T(M_b)` | ✔ `Σ_tail c⁴ ≤ c_Γ² Σ_tail c² = 2c_Γ²V_T(M) ≤ 2c_Γ²V_T(M_b)` | lower bound uses an **upper** bound of `S₄` ✔ | ✔ l.463–464, 541–544 |
| 11 | (5.2) `∂_MΨ = −ΨB_T`, `(t²/2)D_T ≤ B_T ≤ (t²/2)D_T(1+2κc_Γ²t²)` | ✔ `Σ_tail c³c′ ≤ c_Γ² D_T` | ✔ | ✔ l.465, 545–546 |
| 12 | (5.3) `F = sinc·Ψ·(dP − P·B_T)` enclosed on a `(t,M)` box | ✔ product rule; every factor is enclosed for each point of the box, so the interval product contains F (dependency only costs width) | — | ✔ l.515–551; prefix/suffix products, no division by `J₀` |
| 13 | (5.4) Taylor enclosures of `J₀, J₁` on wide balls | ✔ Lagrange remainder with #8; cap by `[−1,1]` and by `x_lo^{-1/2}` only when `x_lo > 1` is certified | ✔ | ✔; property-tested, 0 misses in 7500 samples (§4.3) |
| 14 | (5.5) interval Riemann sum | ✔ mean value property; **there is no quadrature error term**, hence none to audit | — | ✔ tiling of `[0,T₁]` verified with exact rationals l.616–623 |
| 15 | (6.1) and `E₂` | ✔ including both cases of the "max over γ" (γ ∈ S: one factor `(c_γ(M_a)t)^{1/2} ≤ (c_{γ₁}(M_b)t)^{1/2}` is returned; γ ∉ S: multiply by `(c_{γ₁}(M_b)t)^{1/2} ≥ 1`); `∫_{T₁}^∞ t^{1−(K−1)/2}dt = T₁^{(5−K)/2}·2/(K−5)` | needs `c_γ(M_a)T₁ ≥ 1` on S and `c_{γ₁}(M_b)T₁ ≥ 1`, `K > 5` | ✔ l.471–492 (`K ≥ 8` enforced, both conditions checked in Arb) |
| 16 | §4 differentiation under the integral sign | ✔ but see **G1**: the argument is METHOD.md's, not the paper's lemma | — | the dominant is verified per cell (`K ≥ 8`) |
| 17 | §7 assembly, outward rounding of printed numbers | ✔ `up/dn` = `nextafter` of `float(arb.upper()/lower())` (tested on 20000 random balls), decimals rounded with `ROUND_CEILING/FLOOR` from the exact binary value | — | ✔ l.137–172, 642–645; PASS is decided on the Arb ball (`enc.upper() < 0`) |
| 18 | (8a) `Q` tail: `A₂/(2√a_max) ≤ Σ_tail a^{3/2}/2 ≤ min(√a_max·A₁, √(A₁A₂))/2`, `A₁ = V_T/2`, `A₂ = 2V_T − M·D_T` | ✔ (`16Mγ²/(M²+γ²)² = 16M/(M²+γ²) − 16M³/(M²+γ²)²`; `a² ≤ √a_max·a^{3/2}`, `a^{3/2} ≤ √a_max·a`, Cauchy–Schwarz) | ✔ | ✔ `Q_enclosure` l.817–847 |
| 19 | (8b1) root of `σ² = 2` by complex step | ✔ `Im L(s+ih)/h − L′(s) = Σ_{k≥1} (−1)^k h^{2k} L^{(2k+1)}(s)/(2k+1)!`, Cauchy with `\|L(z)\| ≤ \|z\|/Re z` (partial sums of χ bounded by 1 for q = 3, 4) gives `≤ M_L h²/(R³(1−(h/R)²))` | `Re z ≥ s − R ≥ 3/4 > 0` ✔ | ✔ l.287–314 |
| 20 | (8b2) `W_{M′} ≤ Σ_{γ>T} 4/(\|M+iγ\|\|M′+iγ\|) ≤ √(W_M W_{M′})`, `W_M = V_T(M)/(2M²)` | ✔ (termwise `\|M+iγ\| ≤ \|M′+iγ\|`; Cauchy–Schwarz) | ✔ | ✔ `Lambda_enclosure` l.850–878 |

### 2.2 The specific traps the audit was asked to look for

* **Inequality in the wrong direction.**  None.  Every place where a one-sided
  bound enters was checked: `S₄` and `η₂ = 2κc_Γ²` are upper bounds taken with
  `.upper()` (l.463–465), the Gaussian uses the hull of `V_T` over the cell,
  the remainder uses lower bounds `c_γ(M_a)` for the decaying factors and the
  upper bound `c_{γ₁}(M_b)` for the returned factor, `D_max ≥ D(M)` on the cell.
* **(B1)/(B2) outside `c_γt < j₀,₁`.**  They are applied only to unlisted
  zeros and only for `t ≤ T₁`, and `c_{γ_N}(M_b)·T₁ ≤ z₀ = 1` is verified in
  Arb for the `T₁` used (l.474).  Listed zeros never go through (B1)/(B2);
  their Bessel functions are evaluated.  (See MINOR m1, m2 for two guards that
  are weaker than they look but are not exercised.)
* **Tail variance exact only given completeness and GRH — is that stated?**
  Yes: (H1)–(H3) of METHOD.md and the docstring.  Precisely: the identity #2
  sums over *all* zeros, so the "tail" equals `Σ_{γ>γ_N}` only if (i) every
  non-trivial zero is on the critical line (GRH for the whole zero set, not
  just up to T), (ii) `L(½,χ) ≠ 0`, (iii) the list is complete and simple up
  to `γ_N`.  (ii) and (iii) are certified by `certify_arb.py --K 0` (re-run by
  the audit: 511/511 and 537/537 sign changes, winding numbers exactly 511 and
  537, `L(½,χ) > 0`).
* **Zero ordinates read at finite precision.**  Each ordinate is loaded as
  `arb(decimal, "1e-20")` at 128 bits (l.189); the audit checked that this
  ball contains `[dec − 10⁻²⁰, dec + 10⁻²⁰]` exactly (the parse error is
  *added*, radius − 10⁻²⁰ = 3.6·10⁻²⁹ > 0).  The radius is exactly what
  `certify_arb.py` proves: Hardy's Z has certified opposite signs at the two
  exact decimals `dec ± 10⁻²⁰`, the 511 (537) windows are disjoint, and the
  winding number equals the number of windows, so each window holds exactly
  one simple zero and there is no other zero with `−1/2 < Im s < γ_N + 1/4`.
  Mutation tests show the result is moreover insensitive to the radius
  (perturbing `γ₁` by 0.05 moves the enclosure by 0.13 %).
* **Floats on a certified path.**  None.  Floats are used for (a) steering
  (candidate `T₁`, initial breakpoints, refinement order, grid generation) —
  every choice is an exact `Fraction` that is then verified in Arb or by the
  exact tiling test — and (b) printing.  `fr()` rejects float inputs.
* **Monotonicity in M asserted but not proved.**  All uses were located:
  `c_γ(M)` increasing in M / decreasing in γ (trivial), `V_T`, `D_T`
  nondecreasing (#4, with its hypothesis checked at run time), `Q` and `σ²`
  strictly increasing (termwise; `σ²′ = D > 0`).  The sanity test l.427
  (`VTa.upper() ≤ VTb.upper()`) is not a proof and is not needed.  No
  monotonicity of `δ′` in M is used anywhere.
* **Quadrature error term.**  There is none: zone 1 is an interval Riemann
  sum over an exact tiling, zone 2 is a closed-form bound.

## 3. Findings

### 3.1 RIGOR BREAK

None.

### 3.2 GAP

**G1 — The lemma that turns the computed integral into δ′ is not in the paper.**
The program encloses `I(M) := (1/π)∫₀^∞ sinc(t)·∂_MΦ_M(t) dt`.  That
`I(M) = δ′(M)` is the paper's `lem:diffint(c)`, whose hypothesis ("at least 8
ordinates ≤ 2M_a") holds only for `M_a ≥ 12.87` (χ₄) resp. `14.11` (χ₃), i.e.
on a small part of the priority ranges.  METHOD.md §4 replaces it by a
dominated-convergence argument whose dominant the program verifies per cell;
that argument is correct, but (a) it exists only in METHOD.md, (b) it refers
to the proof of `lem:diffint(a)` for the pointwise formula, and (c) it is
silent about cell endpoints.  The complete argument, to be added to the paper
if the certificate is cited:

1. *Pointwise formula without any hypothesis on the ordinates.*  Fix `t ≥ 0`
   and a compact `[M_a, M_b] ⊂ (0,∞)`.  Each `M ↦ c_γ(M) = 4M/√(M²+γ²)` is
   holomorphic on the bounded complex neighbourhood
   `U = {M : dist(M, [M_a, M_b]) < r}`, `r = min(M_a, γ₁)/2` (there
   `Re(M²+γ²) ≥ γ² − (Im M)² ≥ 3γ²/4`), with `|c_γ(M)| ≤ C_U/γ` on U.  Since `|1 − J₀(z)| ≤ (|z|²/4)e^{|Im z|}` and
   `Σ_γ γ^{-2} < ∞`, the product `Φ_M(t) = Π_γ J₀(c_γ(M)t)` converges locally
   uniformly on U, so it is holomorphic in M and (Weierstrass) its
   M-derivative is the limit of the derivatives of the partial products.  The
   derivative of the partial product over `γ ≤ G` is
   `−Σ_{γ≤G} c_γ′t·J₁(c_γt)·Π_{γ′≤G, γ′≠γ} J₀(c_{γ′}t)`; for real M its terms
   are bounded by `c_γc_γ′t²/2` (B4, `|J₀| ≤ 1`), which is summable
   (`Σ c_γc_γ′ = D(M) < ∞`), and each truncated product converges to the full
   one, so by dominated convergence for series the limit is `eq:dPhi`.
2. *Dominant.*  For `M ∈ [M_a, M_b]`: `|∂_MΦ_M(t)| ≤ (t²/2)D_max` for all
   `t ≥ 0`, and `≤ (t²/2)D_max·A·t^{-(K−1)/2}` for `t ≥ T₁` (6.1).  Times
   `|sin t/t| ≤ min(1, 1/t)` this is integrable on `[0,∞)` iff `K > 5`; the
   program enforces `K ≥ 8` and the two inequalities `c_γ(M_a)T₁ ≥ 1` (γ ∈ S),
   `c_{γ₁}(M_b)T₁ ≥ 1` used in (6.1).
3. *Conclusion.*  `Φ_M ∈ L¹` for each M in the cell (B5 on the K factors), so
   `δ(M) = 1/2 + (1/π)∫ sinc·Φ_M` holds on the cell, and the standard theorem
   gives that `δ|_{[M_a,M_b]}` is differentiable with derivative `I(M)`
   (one-sided at `M_a`, `M_b`).  At a breakpoint shared by two certified cells
   the two one-sided derivatives are the same integral, so δ is differentiable
   there.  At the outer endpoints the certificate yields one-sided derivatives
   (right at 3.6 / 4.4; at 17.4 / 19.2 the stretch cells supply the other
   side).  For the conclusion "δ strictly decreasing on the closed interval"
   one-sided derivatives suffice (mean value theorem on each closed cell); for
   a two-sided statement at 3.6 / 4.4 note that the dominant persists on
   `[M_a − ε, M_b]` because `K = 86` resp. `95` there, far above 8.

*Status:* the certificate itself is sound; what is missing is the lemma in the
paper.  Hypothesis to document next to any citation: "δ′(M) is given by
`lem:diffint(c)` in the form of METHOD.md §4 / AUDIT.md G1".

**G2 — The certificate is tied to hypothesis (H3) only by convention.**
`certify_derivative.py` accepts any zero file that has ≥ 100 strictly
increasing lines.  It cannot know whether the file is the one certified by
`certify_arb.py --K 0`, and nothing in a result line records which file was
used.  The mutation tests make the consequence concrete: with `γ₁` deleted
from `chi4_zeros.txt` the checker prints `PASS` for `[3.6, 3.618]` with the
enclosure `[−3.5040e-2, −3.3666e-2]`, which does **not** contain the true
values (`δ′(3.6) = −3.3423e-2`, `δ′(3.618) = −3.3100e-2`); deleting `γ₃`, `γ₁₀₀` or `γ₄₀₀`, or moving
`γ₁` by 0.5, also prints `PASS` (there the enclosure happens to still contain
the truth).  Only gross corruption is caught (disorder, duplicate, a spurious
low zero → "tail variance not certified positive", truncated file).  So the
result files prove the claim only *together with* a run of `certify_arb.py
--K 0` on byte-identical zero files.
*Added hypothesis / fix:* "(H3′) the zero files are those with sha256
`0e73e00a6c1609bdb4eec62c13292435d0beee55750735a51f0fab84024ddfbc` (χ₄) and
`e3075235a847f43eb35057d3998a853f9383e27e958cc5fd3c4115dad7eb8aaa` (χ₃), for
which `certify_arb.py --K 0` passes."  `SHA256SUMS` and stages 0 and 2 of
`replay.sh` implement it; both were run (all certifications PASSED).  It would
be better still if `load_zeros` verified the hash and each result file
recorded it.

**G3 — Provenance: the result files predate the code, and the shipped
reproduction script does not recompute them.**  `results_q4.txt` was completed
at 05:29:21 and `results_q3.txt` at 05:42:01; `certify_derivative.py` was
modified at 05:38 (between the start and the end of the q = 3 run) and again
at 05:58.  `git diff` shows that the changes after commit `f25da0c` touch only
`Q_enclosure` (sharper tail bound), docstrings, an early return when there is
nothing to do, and one header string — not the δ′ path — but a reader cannot
see that from the files.  Moreover `run_all.sh` is resumable: with the result
files present, `sh run_all.sh` skips every cell ("q=4: 0 cells to do") and
re-verifies nothing, and `summarize()` trusts the token `PASS`.
*Closed by the audit:* the priority ranges were recomputed from scratch with
the final code (`replay.sh` stage 4, output in `replay/`).  Result: the same
317 + 296 cells, and **every field of every line (hi, lo, E2, Emodel, rad, T1,
K, N, tcells, VT, DT) is identical** to the audited files (only `secs` and the
line order differ); `aux_results.txt` is reproduced byte for byte.
(q = 4: 1560 s, q = 3: 1322 s with one worker; `replay/results_q*.txt`,
`replay/verify_results.log`.)

### 3.3 MINOR

* **m1** (`Cell.__init__` l.456–462).  The guard "`J₀(z₀) > 0`" does not prove
  `z₀ < j₀,₁` (J₀ is positive again on `(5.52, 8.65)`).  Harmless for the
  hard-coded `z₀ = 1`; the mutant `z₀ = 6` was stopped only because `C₄`
  happened to come out negative.  Should be `z₀ < 2.4048`.
* **m2** (`Cell.integrate`, `Cell.integrand`).  The range condition (5.0) is
  checked in `remainder()`, not where the tail model is used.  In the shipped
  flow `T₁` always comes from `choose_T1 → remainder`, so the check is always
  executed; a caller passing another `T₁` to `integrate` would bypass it.
* **m3** One-sided derivatives at the outer endpoints of a range (see G1.3).
* **m4** (`parse_results`/`summarize`).  The summary re-parses printed bounds
  with `float()` and prints them again rounded up: the `# SUMMARY` of
  `results_q4.txt` says `δ′ ≤ −1.675290e-3` while the cell line (and
  METHOD.md §11) says `−1.675291e-3`.  Both are valid (the round trip can only
  move the printed upper bound up by one unit of the last digit); a consumer
  that uses `float(hi)` directly may be off by half an ulp in either direction
  (immaterial here: the nearest certified upper bound is 1.5·10⁻³ from 0).
* **m5** (`summarize`).  A line marked `PASS` with a positive `hi` is counted
  as certified ("fully certified: True" on a hand-made file).  The driver
  never writes such a line, but the summary is not a verifier;
  `verify_results.py` of this audit is (exact rationals: tiling, `lo ≤ hi < 0`).
* **m6** Black-box insensitivity.  The enclosures are 4 % wide, so many
  rigor-destroying mutations do not change the verdict or even leave the
  truth inside the output (Bessel radii ÷1000, midpoint sampling in t, tail
  variance ×3, dropping `E₂`, …; table in §4.2).  `crosscheck.py` (C) cannot
  detect them either.  Rigor therefore rests on the derivation (section 2),
  the enclosure-property tests (§4.3) and the independent implementation
  (section 5), not on the cross-checks of METHOD.md §10.
* **m7** `run_all.log` was stale (older aux bound) and has since been deleted;
  METHOD.md §12 still describes `run_all.sh` as the reproduction, which only
  recomputes after the result files are deleted (now stated there).
* **m8** `Emodel`, `rad`, `E2` columns are diagnostics; `E2`/`Emodel` are
  printed as upper bounds, fine.  The statement "E2 < 2·10⁻¹⁰ in all runs"
  (§7) is true (max 3.2·10⁻¹¹).
* **m9** For anyone extending `crosscheck.py`: `mpmath.dirichlet(s, chi,
  derivative=2)` returns a wrong value in mpmath 1.4.1 (0.990 instead of
  −0.0107 at s = 4.1, χ₄).  `crosscheck.py`/`selftest` use `mp.diff` and are
  not affected.
* **m10** (H2) is quoted from the paper; the audit re-derived it (table #2)
  and checked it numerically in a way that would expose a missing low zero or
  a wrong closed form: `A(M) := V_T(M)/(8M²)` must be a sum `Σ_{γ>T} 1/(M²+γ²)`
  over zeros above `T`; `tail_sanity.py` verifies in ball arithmetic the
  necessary inequalities `0 < [A(M₁)−A(M₂)]/(M₂²−M₁²) ≤ A(M₁)/(M₂²+T²)` for
  eleven M per modulus, and finds `A(2.5)` within 0.15 % and the implied
  `Σγ⁻⁴` within 0.5 % of the Riemann–von Mangoldt main terms.
* **m11** The stretch files (`results_q*_stretch.txt`) were not part of the
  claim under audit; they were hashed but not replayed.  (They use the same
  code path.)

## 4. Part 2 — attacks on the code

All mutation work was done on the auditor's byte-identical copy
`mutants/certify_derivative_copy.py`; zero-file mutants live in a scratch
directory.  Raw outputs: `mutants/mutation_results.json`,
`mutants/mutation_all.log`, `mutants/proptest_results.json`,
`mutants/artificial_family.log`, `mutants/driver_checks.log`.

### 4.1 python-flint semantics (`flint_semantics.py`, all 69 checks pass)

* Comparisons of balls are `True` only when certain; overlapping balls give
  `False` for `<, >, <=, >=, ==`; NaN/∞ balls compare `False` and are not
  finite; division by a ball containing 0, `sqrt` of a negative ball, `log`
  of a ball containing 0 are non-finite (so `need_finite` / `not (x > 0)`
  guards fail closed).
* `arb(1, 1) >= 0` is `False`: radii are 30-bit `mag`s rounded **up**, the
  ball is `[−ε, 2+ε]`.  Conservative; it means `rad_ball(u)` has radius ≥ u.
* `upper()`, `lower()`, `abs_upper()`, `abs_lower()` return exact (zero-radius)
  bounds on the correct side at precisions 30, 53, 64 (checked with exact
  rational arithmetic); `mid()`, `rad()` are exact.
* `float(arb)` is within one ulp, so `nextafter(float(x.upper()), +inf)` is an
  upper bound (20000 random balls, exact rational comparison).
* `arb("6.02…", "1e-20")` contains `[dec − 10⁻²⁰, dec + 10⁻²⁰]`.
* `union` contains both arguments and everything between; `intersection`
  contains the true intersection and raises `ValueError` when it is empty
  (this is what stops the mutants `J1 → x/2`, `J0 → cos`).
* `x.bessel_j(n)` is `J_n(x)` (argument is the object), `arb(2).zeta(x)` is
  the Hurwitz `ζ(2,x)`, `sinc` is `sin x/x`; checked against mpmath.
* Printing: `Decimal(float)` is exact and `ROUND_CEILING/FLOOR` give outward
  decimals; re-parsing with `float()` is the only non-monotone step (m4).

### 4.2 Mutations — outcome of every one

Test cells: `q=4 [3.6, 3.618]` (first cell), `q=4 [17.3701, 17.4]` (worst
margin), and for zero-file mutants also `q=3 [19.1282, 19.2]`.  "truth" =
the auditor's rigorous point enclosures of δ′ at `M_a`, the midpoint and `M_b`
(relative width ≤ 7·10⁻⁶), computed with the **unmutated** zero lists.
`RAISES`/`FAIL` = the mutated checker refuses (fail-closed); `PASS, EXCLUDES`
= it prints PASS with an enclosure that misses the truth (a wrong certificate
that only an independent evaluation exposes); `PASS, CONTAINS` = it prints
PASS and the truth is still inside.  Last column: enclosure-property tests of
§4.3 (number of sample points *not* contained).

| mutation | what is changed | outcome per test cell | enclosure-property test (L1 Bessel / L2 integrand misses) |
|---|---|---|---|
| `baseline` | unmutated copy | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `tail_sign_flip` | flip the sign of the tail term: dP - P*BT -> dP + P*BT | q=4 [3.6,3.618]: PASS, EXCLUDES truth [-3.2189e-02, -3.0928e-02]<br>q=4 [17.3701,17.4]: PASS, EXCLUDES truth [-1.4986e-03, -1.4400e-03] | L1: 0+0 of 7500; L2: 58/678, 232/678 |
| `dP_sign_flip` | flip the sign of the head derivative dP | q=4 [3.6,3.618]: FAIL (fail-closed) [3.0928e-02, 3.2189e-02]<br>q=4 [17.3701,17.4]: FAIL (fail-closed) [1.4400e-03, 1.4986e-03] |  |
| `sinc_sign_flip` | truth known POSITIVE: delta -> 1 - delta (sinc -> -sinc) | q=4 [3.6,3.618]: FAIL (fail-closed) [3.2596e-02, 3.3926e-02]<br>q=4 [17.3701,17.4]: FAIL (fail-closed) [1.6753e-03, 1.7435e-03] |  |
| `J1_to_J2` | replace J_1 by J_2 in bessel01 | q=4 [3.6,3.618]: PASS, EXCLUDES truth [-6.5410e-03, -6.1730e-03]<br>q=4 [17.3701,17.4]: PASS, EXCLUDES truth [-2.2556e-04, -2.1672e-04] | L1: 1855+2951 of 7500; L2: 368/678, 510/678 |
| `J1_to_x_over_2` | replace J_1(x) by its first Taylor term x/2 | q=4 [3.6,3.618]: RAISES (fail-closed): ValueError: empty intersection<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ValueError: empty intersection |  |
| `J1_to_minus_J1` | replace J_1 by -J_1 | q=4 [3.6,3.618]: FAIL (fail-closed) [3.0928e-02, 3.2190e-02]<br>q=4 [17.3701,17.4]: FAIL (fail-closed) [1.4399e-03, 1.4987e-03] |  |
| `J0_to_cos` | replace J_0 by cos in bessel01 | q=4 [3.6,3.618]: RAISES (fail-closed): ValueError: empty intersection<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ValueError: empty intersection |  |
| `bessel_radius_div1000` | shrink the Bessel Taylor radii by 1000 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3900e-02, -3.2584e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7426e-03, -1.6745e-03] | L1: 5776+5722 of 7500; L2: 470/678, 487/678 |
| `bessel_drop_second_order` | drop the r^2/2 Lagrange term of the Bessel enclosures | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 1526+1630 of 7500; L2: 0/678, 0/678 |
| `t_midpoint_sampling` | evaluate the integrand at the t-midpoint instead of on the t-box | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3565e-02, -3.2957e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7145e-03, -1.7045e-03] | L1: 0+0 of 7500; L2: 119/678, 280/678 |
| `M_left_endpoint_only` | use M = M_a instead of the M-ball for all head coefficients | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.4088e-02, -3.2751e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7461e-03, -1.6777e-03] | L1: 0+0 of 7500; L2: 74/678, 31/678 |
| `V_left_endpoint_only` | tail variance/derivative taken at M_a only | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3929e-02, -3.2599e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7436e-03, -1.6754e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `zero_radius_removed` | zero balls without the 1e-20 radius | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `E2_dropped` | drop the large-t remainder E2 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `tail_model_radius_dropped` | drop both hull factors of the tail model (Psi = G, BT = BT0) | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `tail_model_wrong_direction` | quartic correction with the wrong sign: Psi in G*[1, e^{+C4S4t^4}] | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `gaussian_tail_dropped` | drop the unlisted zeros altogether (G = 1) | q=4 [3.6,3.618]: PASS, EXCLUDES truth [-3.5036e-02, -3.3663e-02]<br>q=4 [17.3701,17.4]: PASS, EXCLUDES truth [-1.8774e-03, -1.8039e-03] | L1: 0+0 of 7500; L2: 108/678, 360/678 |
| `tail_variance_x3` | tail variance and its derivative multiplied by 3 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3494e-02, -3.2182e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7321e-03, -1.6642e-03] | L1: 0+0 of 7500; L2: 498/678, 597/678 |
| `digamma_dropped` | closed form without the digamma term | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: tail variance not certified positive: [-5.05657366190<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: tail variance not certified positive: [-71.8950883656 |  |
| `chi_sign_error` | chi_4(3) = +1 and chi_3(2) = +1 in the Dirichlet series | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: tail variance not certified positive: [-0.19320922179<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: ?/?, 0/678 |
| `series_tail_dropped` | drop the Dirichlet-series tail radii | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `cprime_wrong_exponent` | c' = (4/g)(1+(M/g)^2)^{-1} instead of ^{-3/2} | q=4 [3.6,3.618]: PASS, EXCLUDES truth [-3.6081e-02, -3.4666e-02]<br>q=4 [17.3701,17.4]: PASS, EXCLUDES truth [-2.0259e-03, -1.9465e-03] | L1: 0+0 of 7500; L2: 66/678, 233/678 |
| `c_large_gamma_approx` | c = 4M/g instead of 4M/sqrt(M^2+g^2) in the head | q=4 [3.6,3.618]: PASS, EXCLUDES truth [-2.8962e-02, -2.7827e-02]<br>q=4 [17.3701,17.4]: PASS, EXCLUDES truth [-5.0392e-04, -4.8418e-04] | L1: 0+0 of 7500; L2: 286/678, 589/678 |
| `watson_too_strong` | large-t remainder with abs(J0(x)) <= 1/x instead of x^{-1/2} | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `z0_2p3` | tail-model range z0 = 2.3 (< j_{0,1}: still a valid certificate) | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7436e-03, -1.6753e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `z0_2p5` | tail-model range z0 = 2.5 (> j_{0,1} = 2.4048: INVALID) | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: J0(z0) not positive<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: J0(z0) not positive |  |
| `z0_6` | tail-model range z0 = 6 (INVALID, but J0(6) > 0) | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: tail-model constants not positive<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: tail-model constants not positive |  |
| `z0_range_check_removed_T1x4` | remove the check c_{g_N} T1 <= z0 and integrate the tail model up to 4*T1max | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03] |  |
| `pass_test_inverted` | PASS test inverted (hi > 0) | q=4 [3.6,3.618]: FAIL (fail-closed) [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: FAIL (fail-closed) [-1.7435e-03, -1.6753e-03] |  |
| `zero1_plus_1e-3` | perturb gamma_1 by +1e-3 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3927e-02, -3.2597e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `zero1_plus_0.05` | perturb gamma_1 by +0.05 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3969e-02, -3.2637e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7436e-03, -1.6754e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] | L1: 0+0 of 7500; L2: 0/678, 4/678 |
| `zero1_plus_0.5` | perturb gamma_1 by +0.5 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.4282e-02, -3.2939e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7442e-03, -1.6759e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6132e-03, -1.5500e-03] |  |
| `zero250_plus_0.05` | perturb gamma_250 by +0.05 | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `zero10_out_of_order` | gamma_10 := gamma_11 + 0.01 (order broken) | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: zero balls not strictly ordered<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: zero balls not strictly ordered<br>q=3 [19.1282,19.2]: RAISES (fail-closed): ArithmeticError: zero balls not strictly ordered |  |
| `zero5_duplicated` | gamma_5 listed twice | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: zero balls not strictly ordered<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: zero balls not strictly ordered<br>q=3 [19.1282,19.2]: RAISES (fail-closed): ArithmeticError: zero balls not strictly ordered |  |
| `drop_zero1` | drop gamma_1 from the list | q=4 [3.6,3.618]: PASS, EXCLUDES truth [-3.5040e-02, -3.3666e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7545e-03, -1.6858e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6211e-03, -1.5576e-03] | L1: 0+0 of 7500; L2: 564/678, 441/678 |
| `drop_zero3` | drop gamma_3 from the list | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3936e-02, -3.2605e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7471e-03, -1.6787e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6153e-03, -1.5520e-03] | L1: 0+0 of 7500; L2: 427/678, 368/678 |
| `drop_zero100` | drop gamma_100 from the list | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `drop_zero400` | drop gamma_400 from the list | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `drop_last_zero` | drop the LAST listed zero (still a complete list: valid) | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] | L1: 0+0 of 7500; L2: 0/678, 0/678 |
| `drop_last_100` | drop the last 100 listed zeros (still complete: valid) | q=4 [3.6,3.618]: PASS, CONTAINS truth [-3.3926e-02, -3.2596e-02]<br>q=4 [17.3701,17.4]: PASS, CONTAINS truth [-1.7435e-03, -1.6753e-03]<br>q=3 [19.1282,19.2]: PASS, CONTAINS truth [-1.6127e-03, -1.5495e-03] |  |
| `spurious_zero` | insert a spurious zero at gamma_1 + 1 | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: tail variance not certified positive: [-1.48485243449<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: tail variance not certified positive: [-2.67313427753<br>q=3 [19.1282,19.2]: RAISES (fail-closed): ArithmeticError: tail variance not certified positive: [-2.01187163129 |  |
| `file_truncated_50` | zero file truncated to 50 lines | q=4 [3.6,3.618]: RAISES (fail-closed): ArithmeticError: zero file too short<br>q=4 [17.3701,17.4]: RAISES (fail-closed): ArithmeticError: zero file too short<br>q=3 [19.1282,19.2]: RAISES (fail-closed): ArithmeticError: zero file too short |  |

Reading of the table:

* **Sign and "truth positive" tests fail closed.**  `sinc → −sinc` (i.e.
  δ → 1 − δ, whose derivative is positive), `dP → −dP`, `J₁ → −J₁` and the
  inverted PASS test all end in `FAIL` with a certified *positive* enclosure.
* **Artificial amplitude family with known positive derivative**
  (`mutants/mutate.py art`): for the dilation family `c_γ(u) = (M₀/u)c_γ(M₀)`
  one has `δ̃(u) = P(S > −u/M₀)`, `δ̃′(u) = f_S(−u/M₀)/M₀ > 0`.  Feeding this
  family to the unmutated `Cell.choose_T1 / integrate` and the verbatim PASS
  test gives

  | family | checker enclosure of δ̃′ | float reference `f_S(−u_a/M₀)/M₀` | PASS |
  |---|---|---|---|
  | q=4, M₀=5, u ∈ [5, 5.02] | [2.075922e-02, 2.160608e-02] | 2.118569e-02 | **False** |
  | q=4, M₀=17.4, u ∈ [17.35, 17.4] | [2.414135e-03, 2.512623e-03] | 2.463466e-03 | **False** |
  | q=3, M₀=8, u ∈ [8, 8.03] | [9.956467e-03, 1.036269e-02] | 1.016047e-02 | **False** |

  The checker fails closed and its positive enclosure contains the
  independently computed value.
* **Wrong formulas are either refused or produce an enclosure that the
  independent values contradict.**  Refused: `J₁ → x/2` and `J₀ → cos` (empty
  intersection with the a-priori caps of `bessel01`).  Accepted with an
  enclosure that excludes the truth: tail sign flip, `J₁ → J₂`, wrong `c`,
  wrong `c′`, tail dropped.  A checker cannot in general detect that it
  computes the wrong quantity; this is what section 5 is for.
* **The range check (5.0).**  The mutant that removes the check
  `c_{γ_N}T₁ ≤ z₀` and integrates the tail model up to `4·T₁max` (where
  `c_γt` reaches 4 > j₀,₁) changes no printed digit, because the head product
  is astronomically small there.  That the check is present and is applied to
  the `T₁` actually used was therefore established by reading (table #9), not
  by testing.
* **Closed-form errors are caught by the positivity of the tail variance**
  (digamma dropped, wrong χ at small M, spurious zero) — a useful built-in
  consistency check.
* **Zero-list mutations:** disorder/duplicate/truncation are refused; a
  deleted or displaced zero is accepted (G2).  Deleting a *low* zero (γ₁, γ₃)
  is exposed by the integrand-level test L2 (427–564 of 678 sample points
  missed) even where the final enclosure still contains the truth; deleting
  γ₁₀₀ or γ₄₀₀ is invisible at every level, because such a zero still
  satisfies `c_γT₁ ≤ z₀` and the exact tail variance absorbs it (the
  certificate is then wrong only in the 10⁻⁷-level constants `S₄`, `η₂`).
  Deleting the *last* zeros is a legitimate certificate (shorter but still
  complete list) and is sound.
* **Radius-shrinking mutations are invisible at the output** (m6) and are
  exposed only by the property tests below; the three that concern bounds on
  the *unlisted* zeros (`tail_model_radius_dropped`,
  `tail_model_wrong_direction`, `series_tail_dropped`) act at the 10⁻⁷ level,
  below anything a numerical experiment can resolve, and were verified by
  derivation only (table #5, #10, #11).

### 4.3 Enclosure-property tests of the unmutated code (`mutants/proptest.py`)

* L1 — `bessel01(ball)` against Arb at 200 bits (spot-checked against mpmath)
  on 1500 random balls (including balls `[0, 2x₀]`, radii up to 1.5, and
  balls with `x_lo > 1` where the Watson cap acts), 5 sample points each:
  **0 misses out of 7500** for J₀ and for J₁.
* L2 — `Cell.integrand(t-box)` against the auditor's independent rigorous
  point values of `F(t, M)` at the corners and midpoints of 76 `(t, M)` boxes
  of three widths, for the cells `[3.6, 3.618]` and `[17.3701, 17.4]`:
  **0 misses out of 678 + 678**.

  This is a direct, quantitative confirmation that the audited integrand
  enclosure contains the true integrand (the two programs share no code and
  treat zeros 41…511 completely differently).

### 4.4 Driver

* Range driver on a failing mutant: writes `# not certified …` and `… FAIL`
  lines and `fully certified: False`; exceptions in a worker are caught and
  reported as "not certified".
* `summarize()` detects a gap and a `FAIL` line; accepts overlapping cells
  (correct) and a `PASS` line with positive `hi` (m5).
* `certify_cell(4, 3.6, 3.618)` with float endpoints raises `TypeError`.

## 5. Part 3 — independent second implementation

### 5.1 `indep_enclose.py` (rigorous)

Written from the definition, sharing with the audited code only the zero
files, the hypotheses and FLINT.  Differences of route:

| | audited | auditor |
|---|---|---|
| quadrature on `[0,T₁]` | adaptive interval Riemann sums, ~550 t-boxes | Arb's rigorous Gauss–Legendre integrator `acb.integral` (Petras); the integrand is enclosed on **complex** balls, non-finite outside its domain of validity |
| zeros with explicit Bessel factors | all 511 / 537 | the first N′ = 36…43 only (chosen per cell) |
| all other zeros | unlisted only: Gaussian × `[e^{−C₄S₄t⁴}, 1]`, first term of φ | listed zeros N′+1…N **and** unlisted: `−log J₀(z) = Σ (λ_n/n) z^{2n}` through `z¹²` with exact rational Rayleigh sums; listed power sums `S_{2n}`, `U_{2n}` evaluated exactly, unlisted parts in `[0, c_Γ^{2n−2}·2V_T]`; monotone-ratio remainder at `z₀ = 3/2`, valid for complex t |
| closed forms `σ²`, `D` | Dirichlet series over prime powers ≤ 20000 + digamma + Hurwitz ζ(2,·) | Hurwitz-zeta **power series** (`arb_series.zeta`) and `lgamma` series: `log Λ(s+x)` to order 2 |
| Bessel on balls | real second-order Taylor | complex Lipschitz/Taylor with `\|J_n(z)\| ≤ e^{\|Im z\|}`, `\|J₁(z)\| ≤ \|z\|e^{\|Im z\|}/2` |
| coefficients over a cell | ball arithmetic in M | monotone hulls of endpoint values |
| large t | `\|J₀(x)\| ≤ x^{-1/2}`, `T₁ ≈ 10.5 … 1.9` | `\|J₀(x)\| ≤ √(2/(πx))`, own `T₁ ≈ 6.7 … 1.3` |
| zero parsing | `arb(str, "1e-20")` | exact `Fraction` → ball, + 10⁻²⁰ |
| M-cell | one ball | 4 sub-cells, union |

For each fixed M the true integrand `F(t; M) = sinc(t)·∂_MΦ_M(t)` is an entire
function of t; `Cell.integrand` returns, for every complex ball inside
`|t| ≤ 1.3·T₁`, a ball containing `F(ball; M)` for all M in the (sub-)cell, so
`acb.integral` returns a ball containing `∫₀^{T₁} F(t; M)dt` for every such M.
The derivation of every bound is in the docstring of `indep_enclose.py`.

**Result (all cells, not a sample):** `indep_q4.txt`, `indep_q3.txt`.

| modulus | cells | meet the audited enclosure | independent enclosure < 0 | independent ⊂ audited | worst independent upper bound | max rel. width |
|---|---|---|---|---|---|---|
| q = 4, [3.6, 17.4] | 317 | 317 | 317 | 317 | −1.706259e-03 on [17.3701, 17.4] | 1.18 % |
| q = 3, [4.4, 19.2] | 296 | 296 | 296 | 296 | −1.574568e-03 on [19.1282, 19.2] | 1.19 % |

Fourteen cells per modulus (both endpoints of the range, evenly spread cells,
and the three cells with the worst — closest to 0 — reported margin):

**q = 4**

| cell [M_a, M_b] | audited lo | audited hi | independent lo | independent hi | meet | indep. < 0 | indep. inside audited | T1 / K / N′ (indep.) |
|---|---|---|---|---|---|---|---|---|
| [3.6, 3.618] | -3.392577e-02 | -3.259580e-02 | -3.345853e-02 | -3.306509e-02 | True | True | True | 6.6875 / 86 / 40 |
| [4.1583, 4.179] | -2.558579e-02 | -2.458343e-02 | -2.523378e-02 | -2.493704e-02 | True | True | True | 5.8906 / 88 / 41 |
| [4.78, 4.8038] | -1.945721e-02 | -1.869432e-02 | -1.918895e-02 | -1.896378e-02 | True | True | True | 5.2188 / 90 / 42 |
| [5.5223, 5.5499] | -1.466194e-02 | -1.408704e-02 | -1.445952e-02 | -1.429039e-02 | True | True | True | 4.2031 / 82 / 38 |
| [6.3804, 6.4123] | -1.107292e-02 | -1.063890e-02 | -1.091951e-02 | -1.079308e-02 | True | True | True | 3.7344 / 85 / 39 |
| [7.3717, 7.4085] | -8.391459e-03 | -8.062575e-03 | -8.274582e-03 | -8.180037e-03 | True | True | True | 3.3281 / 88 / 41 |
| [8.4754, 8.5177] | -6.442840e-03 | -6.190366e-03 | -6.352712e-03 | -6.280962e-03 | True | True | True | 2.7188 / 82 / 37 |
| [9.7927, 9.8416] | -4.918902e-03 | -4.726061e-03 | -4.849783e-03 | -4.795568e-03 | True | True | True | 2.4375 / 85 / 39 |
| [11.3152, 11.3717] | -3.768226e-03 | -3.620474e-03 | -3.715065e-03 | -3.673939e-03 | True | True | True | 2.1875 / 90 / 41 |
| [13.0744, 13.1397] | -2.895594e-03 | -2.782169e-03 | -2.854647e-03 | -2.823325e-03 | True | True | True | 1.7969 / 83 / 38 |
| [15.0323, 15.1074] | -2.251132e-03 | -2.162924e-03 | -2.219191e-03 | -2.195028e-03 | True | True | True | 1.6250 / 87 / 40 |
| [17.1978, 17.2837] | -1.769672e-03 | -1.700301e-03 | -1.744478e-03 | -1.725631e-03 | True | True | True | 1.3438 / 81 / 36 |
| [17.2837, 17.3701] | -1.754026e-03 | -1.685286e-03 | -1.729061e-03 | -1.710378e-03 | True | True | True | 1.3438 / 82 / 37 |
| [17.3701, 17.4] | -1.743515e-03 | -1.675291e-03 | -1.712669e-03 | -1.706259e-03 | True | True | True | 1.3438 / 82 / 37 |

**q = 3**

| cell [M_a, M_b] | audited lo | audited hi | independent lo | independent hi | meet | indep. < 0 | indep. inside audited | T1 / K / N′ (indep.) |
|---|---|---|---|---|---|---|---|---|
| [4.4, 4.422] | -2.739566e-02 | -2.632233e-02 | -2.701847e-02 | -2.670149e-02 | True | True | True | 5.7656 / 85 / 38 |
| [5.033, 5.0581] | -2.104574e-02 | -2.022077e-02 | -2.075649e-02 | -2.051158e-02 | True | True | True | 5.1094 / 86 / 39 |
| [5.7569, 5.7856] | -1.613017e-02 | -1.549857e-02 | -1.590915e-02 | -1.572078e-02 | True | True | True | 4.5312 / 88 / 40 |
| [6.5528, 6.5855] | -1.247272e-02 | -1.198380e-02 | -1.230163e-02 | -1.215588e-02 | True | True | True | 4.0469 / 90 / 41 |
| [7.4959, 7.5333] | -9.554015e-03 | -9.179632e-03 | -9.422793e-03 | -9.311624e-03 | True | True | True | 3.2969 / 82 / 37 |
| [8.5748, 8.6176] | -7.330599e-03 | -7.043190e-03 | -7.229540e-03 | -7.144832e-03 | True | True | True | 2.9531 / 84 / 38 |
| [9.8095, 9.8585] | -5.638790e-03 | -5.417795e-03 | -5.560826e-03 | -5.496243e-03 | True | True | True | 2.6406 / 87 / 40 |
| [11.222, 11.2781] | -4.351560e-03 | -4.181159e-03 | -4.291206e-03 | -4.241872e-03 | True | True | True | 2.3750 / 90 / 41 |
| [12.838, 12.9021] | -3.370150e-03 | -3.237987e-03 | -3.323055e-03 | -3.285338e-03 | True | True | True | 1.9531 / 84 / 38 |
| [14.6143, 14.6873] | -2.643236e-03 | -2.539669e-03 | -2.606192e-03 | -2.576916e-03 | True | True | True | 1.7656 / 86 / 39 |
| [16.7196, 16.8031] | -2.060445e-03 | -1.979686e-03 | -2.031431e-03 | -2.008856e-03 | True | True | True | 1.6094 / 91 / 41 |
| [18.9385, 19.0331] | -1.640366e-03 | -1.576067e-03 | -1.617178e-03 | -1.599379e-03 | True | True | True | 1.3281 / 84 / 37 |
| [19.0331, 19.1282] | -1.625556e-03 | -1.561841e-03 | -1.602577e-03 | -1.584941e-03 | True | True | True | 1.3281 / 84 / 37 |
| [19.1282, 19.2] | -1.612692e-03 | -1.549494e-03 | -1.587717e-03 | -1.574568e-03 | True | True | True | 1.3281 / 84 / 37 |

Robustness of the independent program (`indep_variants.log`): for 7 cells /
points, 9 parameter sets (series order 3, 6, 9; `z₀` = 1, 1.5, 2.2; remainder
goal 10⁻⁷ … 10⁻²⁰, i.e. `T₁` from 4.6 to 9.8 and N′ from 21 to 74; domain
margin; 110 vs 200 bits) give mutually intersecting enclosures in every case.

### 5.2 Plain high-precision evaluation (`plain_eval.py`, mpmath only, 30 digits)

No FLINT.  All listed zeros explicit, unlisted zeros as a Gaussian with the
exact variance, 24-point Gauss–Legendre on 18 panels up to `T = 18/σ`
(doubling T and quadrupling the panels changes δ′ by < 10⁻³¹).  20 points per
modulus, endpoints included.  **Every value lies inside the audited `[lo, hi]`
of its cell and inside the auditor's rigorous point enclosure** (whose relative
width is 10⁻⁶ … 7·10⁻⁶, entirely due to the unlisted zeros):

**q = 4**

| M | plain mpmath δ′(M) | audited cell | audited [lo, hi] | in cell interval | auditor's rigorous point enclosure | in |
|---|---|---|---|---|---|---|
| 3.60000 | -3.3423211414e-02 | [3.6,3.618] | [-3.392577e-02,-3.259580e-02] | INSIDE | [-3.3423231391e-02,-3.3423191169e-02] | INSIDE |
| 4.32632 | -2.3320497361e-02 | [4.3057,4.3272] | [-2.389345e-02,-2.295689e-02] | INSIDE | [-2.3320513312e-02,-2.3320481647e-02] | INSIDE |
| 5.05263 | -1.7191153667e-02 | [5.049,5.0742] | [-1.747354e-02,-1.678906e-02] | INSIDE | [-1.7191167048e-02,-1.7191140767e-02] | INSIDE |
| 5.77895 | -1.3219760438e-02 | [5.7754,5.8042] | [-1.343467e-02,-1.290814e-02] | INSIDE | [-1.3219772140e-02,-1.3219749380e-02] | INSIDE |
| 6.50526 | -1.0507761601e-02 | [6.4765,6.5088] | [-1.075794e-02,-1.033617e-02] | INSIDE | [-1.0507772138e-02,-1.0507751795e-02] | INSIDE |
| 7.23158 | -8.5748242447e-03 | [7.2263,7.2624] | [-8.716711e-03,-8.374988e-03] | INSIDE | [-8.5748339625e-03,-8.5748153360e-03] | INSIDE |
| 7.95789 | -7.1476937722e-03 | [7.9437,7.9834] | [-7.280699e-03,-6.995239e-03] | INSIDE | [-7.1477028323e-03,-7.1476855446e-03] | INSIDE |
| 8.68421 | -6.0626517700e-03 | [8.646,8.6892] | [-6.206047e-03,-5.962724e-03] | INSIDE | [-6.0626603445e-03,-6.0626440536e-03] | INSIDE |
| 9.41053 | -5.2171724752e-03 | [9.4101,9.4571] | [-5.297108e-03,-5.089587e-03] | INSIDE | [-5.2171806310e-03,-5.2171651769e-03] | INSIDE |
| 10.13684 | -4.5445427049e-03 | [10.0899,10.1403] | [-4.653738e-03,-4.471356e-03] | INSIDE | [-4.5445505239e-03,-4.5445357682e-03] | INSIDE |
| 10.86316 | -3.9998748001e-03 | [10.819,10.873] | [-4.091630e-03,-3.931254e-03] | INSIDE | [-3.9998823268e-03,-3.9998681605e-03] | INSIDE |
| 11.58947 | -3.5520669629e-03 | [11.543,11.6007] | [-3.633098e-03,-3.490704e-03] | INSIDE | [-3.5520742276e-03,-3.5520605779e-03] | INSIDE |
| 12.31579 | -3.1790066456e-03 | [12.3154,12.3769] | [-3.227985e-03,-3.101487e-03] | INSIDE | [-3.1790136964e-03,-3.1790004832e-03] | INSIDE |
| 13.04211 | -2.8646115634e-03 | [13.0094,13.0744] | [-2.921895e-03,-2.807399e-03] | INSIDE | [-2.8646184182e-03,-2.8646055980e-03] | INSIDE |
| 13.76842 | -2.5969484718e-03 | [13.7424,13.8111] | [-2.645945e-03,-2.542196e-03] | INSIDE | [-2.5969551498e-03,-2.5969426824e-03] | INSIDE |
| 14.49474 | -2.3670041300e-03 | [14.4447,14.5169] | [-2.418458e-03,-2.323677e-03] | INSIDE | [-2.3670106566e-03,-2.3669985021e-03] | INSIDE |
| 15.22105 | -2.1678623357e-03 | [15.1829,15.2588] | [-2.211220e-03,-2.124537e-03] | INSIDE | [-2.1678687190e-03,-2.1678568482e-03] | INSIDE |
| 15.94737 | -1.9941404696e-03 | [15.8796,15.9589] | [-2.040391e-03,-1.960438e-03] | INSIDE | [-1.9941467253e-03,-1.9941351093e-03] | INSIDE |
| 16.67368 | -1.8415958761e-03 | [16.6081,16.6911] | [-1.883245e-03,-1.809397e-03] | INSIDE | [-1.8416020139e-03,-1.8415906361e-03] | INSIDE |
| 17.40000 | -1.7068458353e-03 | [17.3701,17.4] | [-1.743515e-03,-1.675291e-03] | INSIDE | [-1.7068518664e-03,-1.7068407014e-03] | INSIDE |

**q = 3**

| M | plain mpmath δ′(M) | audited cell | audited [lo, hi] | in cell interval | auditor's rigorous point enclosure | in |
|---|---|---|---|---|---|---|
| 4.40000 | -2.6990439556e-02 | [4.4,4.422] | [-2.739566e-02,-2.632233e-02] | INSIDE | [-2.6990460930e-02,-2.6990416663e-02] | INSIDE |
| 5.17895 | -1.9597572064e-02 | [5.1598,5.1855] | [-2.003695e-02,-1.925192e-02] | INSIDE | [-1.9597589329e-02,-1.9597553881e-02] | INSIDE |
| 5.95789 | -1.4846323397e-02 | [5.9315,5.9611] | [-1.520144e-02,-1.460585e-02] | INSIDE | [-1.4846337988e-02,-1.4846308361e-02] | INSIDE |
| 6.73684 | -1.1631536625e-02 | [6.718,6.7515] | [-1.187128e-02,-1.140573e-02] | INSIDE | [-1.1631549321e-02,-1.1631523797e-02] | INSIDE |
| 7.51579 | -9.3640575824e-03 | [7.4959,7.5333] | [-9.554015e-03,-9.179632e-03] | INSIDE | [-9.3640689302e-03,-9.3640463165e-03] | INSIDE |
| 8.29474 | -7.7089442767e-03 | [8.2808,8.3222] | [-7.850355e-03,-7.542580e-03] | INSIDE | [-7.7089546030e-03,-7.7089341939e-03] | INSIDE |
| 9.07368 | -6.4655092159e-03 | [9.0578,9.103] | [-6.585167e-03,-6.326932e-03] | INSIDE | [-6.4655187452e-03,-6.4655000460e-03] | INSIDE |
| 9.85263 | -5.5081586118e-03 | [9.8095,9.8585] | [-5.638790e-03,-5.417795e-03] | INSIDE | [-5.5081675164e-03,-5.5081501559e-03] | INSIDE |
| 10.63158 | -4.7553543767e-03 | [10.6234,10.6765] | [-4.834208e-03,-4.644764e-03] | INSIDE | [-4.7553627664e-03,-4.7553465047e-03] | INSIDE |
| 11.41053 | -4.1524647570e-03 | [11.391,11.4479] | [-4.229121e-03,-4.063280e-03] | INSIDE | [-4.1524727130e-03,-4.1524573608e-03] | INSIDE |
| 12.18947 | -3.6618655803e-03 | [12.1533,12.214] | [-3.738362e-03,-3.591891e-03] | INSIDE | [-3.6618731719e-03,-3.6618585800e-03] | INSIDE |
| 12.96842 | -3.2570093931e-03 | [12.9666,13.0314] | [-3.307325e-03,-3.177728e-03] | INSIDE | [-3.2570166804e-03,-3.2570027361e-03] | INSIDE |
| 13.74737 | -2.9187553135e-03 | [13.6974,13.7658] | [-2.983480e-03,-2.866504e-03] | INSIDE | [-2.9187623295e-03,-2.9187489417e-03] | INSIDE |
| 14.52632 | -2.6330325703e-03 | [14.4693,14.5416] | [-2.692814e-03,-2.587291e-03] | INSIDE | [-2.6330393484e-03,-2.6330264554e-03] | INSIDE |
| 15.30526 | -2.3893151524e-03 | [15.2848,15.3612] | [-2.431790e-03,-2.336504e-03] | INSIDE | [-2.3893217258e-03,-2.3893092621e-03] | INSIDE |
| 16.08421 | -2.1796032182e-03 | [16.066,16.1463] | [-2.217524e-03,-2.130653e-03] | INSIDE | [-2.1796096044e-03,-2.1795975263e-03] | INSIDE |
| 16.86316 | -1.9977288299e-03 | [16.8031,16.8871] | [-2.041666e-03,-1.961602e-03] | INSIDE | [-1.9977350494e-03,-1.9977233169e-03] | INSIDE |
| 17.64211 | -1.8388738219e-03 | [17.5741,17.6619] | [-1.880311e-03,-1.806594e-03] | INSIDE | [-1.8388798890e-03,-1.8388684657e-03] | INSIDE |
| 18.42105 | -1.6992291950e-03 | [18.3804,18.4723] | [-1.732254e-03,-1.664367e-03] | INSIDE | [-1.6992351248e-03,-1.6992239798e-03] | INSIDE |
| 19.20000 | -1.5757506626e-03 | [19.1282,19.2] | [-1.612692e-03,-1.549494e-03] | INSIDE | [-1.5757564724e-03,-1.5757455767e-03] | INSIDE |

Finite differences of the plain δ(M) (central, h = 10⁻³, Richardson) against
the plain derivative formula: 

| q | M | formula δ′(M) | Richardson finite difference of δ | abs. difference |
|---|---|---|---|---|
| 4 | 3.7 | -0.03168473022627092 | -0.0316847302262704 | 5.18e-16 |
| 4 | 9 | -0.005670165781535199 | -0.005670165781535195 | 3.29e-18 |
| 4 | 17.3 | -0.001724448381748536 | -0.001724448381748536 | 6.23e-20 |
| 3 | 4.5 | -0.02583239095525719 | -0.025832390955257 | 1.9e-16 |
| 3 | 10 | -0.005352004568448188 | -0.005352004568448186 | 2.14e-18 |
| 3 | 19.1 | -0.001590788904080494 | -0.001590788904080494 | 4.2e-20 |

### 5.3 Certified densities (`race_mono_grid_certified_results.txt`)

For the six consecutive pairs of grid points inside each range
(`M = 4.5, 5.5, 6.5, 8.5, 10.5, 12.5, 15.5`): (a) by the mean value theorem
the certified difference quotient `[δ(M₂) − δ(M₁)]/(M₂ − M₁)` must meet the
hull of the δ′ enclosures over `[M₁, M₂]`; (b) `Σ_cells [lo, hi]·|cell ∩ [M₁,M₂]|`
must meet the certified interval for `δ(M₂) − δ(M₁)`.  Both hold for all pairs,
with the audited enclosures (4 % wide) and with the auditor's (1.2 % wide):

| q | [M₁, M₂] | certified δ(M₂) − δ(M₁) | ∫δ′ from audited cells | ∫δ′ from independent enclosures |
|---|---|---|---|---|
| 4 | [4.5, 5.5] | [−0.017721605, −0.017721545] | [−0.018075507, −0.017367074] | [−0.017826499, −0.017617202] |
| 4 | [5.5, 6.5] | [−0.012372761, −0.012372701] | [−0.012619796, −0.012125193] | [−0.012445240, −0.012300579] |
| 4 | [6.5, 8.5] | [−0.016273884, −0.016273821] | [−0.016598822, −0.015948271] | [−0.016367678, −0.016180600] |
| 4 | [8.5, 10.5] | [−0.010357394, −0.010357326] | [−0.010564450, −0.010150396] | [−0.010416238, −0.010299407] |
| 4 | [10.5, 12.5] | [−0.007254428, −0.007254354] | [−0.007399207, −0.007109221] | [−0.007294868, −0.007214129] |
| 4 | [12.5, 15.5] | [−0.007635591, −0.007635510] | [−0.007787922, −0.007482685] | [−0.007677585, −0.007593588] |
| 3 | [4.5, 5.5] | [−0.021207978, −0.021207895] | [−0.021631732, −0.020783984] | [−0.021334429, −0.021082887] |
| 3 | [5.5, 6.5] | [−0.014740864, −0.014740785] | [−0.015034965, −0.014445780] | [−0.014828844, −0.014653034] |
| 3 | [6.5, 8.5] | [−0.019145422, −0.019145346] | [−0.019527818, −0.018762426] | [−0.019259610, −0.019032194] |
| 3 | [8.5, 10.5] | [−0.011954593, −0.011954516] | [−0.012193202, −0.011715296] | [−0.012024796, −0.011884718] |
| 3 | [10.5, 12.5] | [−0.008241761, −0.008241681] | [−0.008406293, −0.008076833] | [−0.008289445, −0.008194355] |
| 3 | [12.5, 15.5] | [−0.008553804, −0.008553718] | [−0.008724629, −0.008382711] | [−0.008602555, −0.008505450] |

(The midpoints of the independent integral intervals agree with the certified
differences to 8·10⁻⁷ or better.)

### 5.4 Auxiliary certificates (`aux_check.py`)

Every line of `aux_results.txt` was re-certified with the auditor's closed
forms (deflated Hurwitz series, valid down to s = 1) and zero balls:
`Q_{χ₄}(3.68) < 1` and `Q_{χ₃}(4.47) < 1` (also with the cruder tail bound
`(M/γ_N)V_T/2` alone); `Q(3.688854) < 1 < Q(3.689059)` (χ₄),
`Q(4.483812) < 1 < Q(4.484161)` (χ₃); `σ²(1.82391477) ≤ 2 < σ²(1.82391478)`
(χ₄), `σ²(2.12775095) ≤ 2 < σ²(2.12775096)` (χ₃); all twelve `Λ(M, M′)`
enclosures meet the auditor's and have the same certified sign; the eight
`σ_m²` intervals contain the auditor's 20-digit balls.

## 6. Part 4 — packaging

* `SHA256SUMS` — sha256 of the two zero lists, `certify_arb.py`, the certified
  density grid, and every file of `monotonicity/cert/` as it stood at the end
  of the audit (paths relative to the repository root; check with
  `sha256sum -c monotonicity/cert-audit/SHA256SUMS` from the root).
* `replay.sh` — re-verifies from scratch; stages, what each one proves and the
  measured run times are in its header.  Stage 4 (the two ranges) refuses to
  resume an existing output file, because the audited driver would skip cells
  that are already marked PASS.
* What was actually run for this report: stages 1, 2, 3, 6, 7 by
  `OUT=replay STAGES="1 2 3 6 7" sh replay.sh`; stage 4 by the same two
  commands with one worker (`replay/_queue_ranges.sh`, because of the 2-core
  limit of the audit); stages 0 and 5 afterwards; stage 8 by the same commands
  from the top-level directory (outputs in `mutants/`).  Logs: `replay/`.
* Measured run times (wall clock, one process, shared and heavily loaded
  machine): stage 1: 3 s; stage 2 (`certify_arb.py --K 0`): 138 s; stage 3: 5 s;
  stage 4: 1560 s (q = 4) + 1322 s (q = 3) with one worker — about 25 min with
  the default two workers; stage 5: < 1 s; stage 6 (independent enclosure of
  all 613 cells): 361 s; stage 7 (plain mpmath values, grid, finite
  differences): 375 s; stage 8 (mutations, optional): about 12 min.  Total
  without stage 8: about 63 min on one core, about 40 min with `WORKERS=2`.
* The two runs of the independent program (top level and `replay/`) are
  identical in every digit; so are the two runs of `plain_eval.py`.
* Files of this audit: `AUDIT.md` (this report), `SHA256SUMS`, `replay.sh`;
  `indep_enclose.py` (independent rigorous enclosure), `plain_eval.py`
  (mpmath), `compare.py`, `verify_results.py`, `aux_check.py`,
  `tail_sanity.py`, `flint_semantics.py`, `indep_variants.py`,
  `make_tables.py`; `mutants/` (copy of the audited code, `mutate.py`,
  `proptest.py`, results); raw outputs `indep_q{3,4}.txt`, `plain_q{3,4}.txt`,
  `compare_*.log`, `certify_arb_K0.log`, `indep_variants.log`, `tables.md`;
  `replay/` (outputs of `replay.sh`).

## 7. Final statement: what is certified, and under which hypotheses

Let χ be the real primitive character mod q ∈ {3, 4}.

**Hypotheses.**
(A) GRH(χ): every non-trivial zero of `L(s, χ)` has real part 1/2.  [Used to
define `c_γ`, and — for the whole zero set, not only up to T — to identify
`σ²(M) − Σ_{listed}` with the sum over the unlisted ordinates.]
(B) The classical facts used by `certify_arb.py` (no zeros of `L(s,χ)` with
`Re s ≥ 1` or `Re s ≤ 0` other than the trivial ones) and by METHOD.md §3
(Hadamard product of `J₀`, Bessel's integral, Watson's inequality).
(C) Correctness of FLINT/Arb as shipped in python-flint 0.9.0 (ball
arithmetic, `bessel_j`, `digamma`, Hurwitz zeta, `sinc`, `exp`, `log`, `sqrt`,
Dirichlet `L` and Hardy `Z`), and of CPython's `int`, `Fraction`, `Decimal`.
(D) The zero files are the ones pinned in `SHA256SUMS`.
Not hypotheses any more: completeness and simplicity of the zero lists up to
`T = 640.165707` (χ₄, 511 ordinates) and `700.024747` (χ₃, 537 ordinates), the
`10⁻²⁰` accuracy of each listed ordinate, and `L(1/2, χ) ≠ 0` are *proved* by
`certify_arb.py --K 0` under (B), (C), (D) (re-run by the audit).  LI is not
used at all.

**Certified statement.**  Under (A)–(D), let γ run over the positive ordinates
of `L(s,χ)`, `c_γ(M) = 4M/√(M²+γ²)`, and
`δ(M) = 1/2 + (1/π)∫₀^∞ (sin t/t)·Π_γ J₀(c_γ(M)t) dt`.  Then δ is
differentiable on `[3.6, 17.4]` (q = 4) resp. `[4.4, 19.2]` (q = 3) (one-sided
at the two endpoints of each range, which is all that monotonicity needs), `δ′(M) = (1/π)∫₀^∞ (sin t/t)·∂_M Π_γ J₀(c_γ(M)t) dt`,
and for every cell `[M_a, M_b]` listed in `results_q4.txt` / `results_q3.txt`
and every real `M ∈ [M_a, M_b]`

    lo ≤ δ′(M) ≤ hi < 0,

with the printed `lo`, `hi` of that cell.  The cells tile the two ranges
exactly, so `δ′(M) < 0` for every M in `[3.6, 17.4]` resp. `[4.4, 19.2]`, with
`δ′(M) ≤ −1.675291·10⁻³` (q = 4) and `δ′(M) ≤ −1.549494·10⁻³` (q = 3)
throughout, and δ is strictly decreasing on each closed range.  (The step
"enclosed integral = δ′" is G1: proved in METHOD.md §4 / AUDIT.md G1, not by
the paper's `lem:diffint`.)  The same statement, with sharper bounds
(`δ′ ≤ −1.706259·10⁻³` resp. `−1.574568·10⁻³`), is proved independently by
`indep_enclose.py` on the same cells.

Also certified under the same hypotheses (`aux_results.txt`, re-certified
independently): `Q_{χ₄}(3.68) < 1`, `Q_{χ₃}(4.47) < 1`;
`3.688854 < M₀(χ₄) < 3.689059`, `4.483812 < M₀(χ₃) < 4.484161`;
`1.82391477 ≤ M*(χ₄) < 1.82391478`, `2.12775095 ≤ M*(χ₃) < 2.12775096` for the
root of `σ²_M = 2`; and the enclosures and signs of `Λ(M, M′)` for
`M < M′` in `{1/2, 3/2, 5/2, 7/2}` listed there (positive for `(m, m′) = (0,1),
(0,2), (0,3)` for both characters and for `(1,2)` when q = 3; negative for
`(1,2)` when q = 4 and for `(1,3)`, `(2,3)` for both).

**Not certified here:** that `δ(M)` is the logarithmic density of the weighted
race (needs LI and the paper's `thm:interp`); anything outside the two ranges
(the stretch files were not audited); the paper's `thm:monotone-large` and the
coupling theorem, which the certificate is meant to connect.
