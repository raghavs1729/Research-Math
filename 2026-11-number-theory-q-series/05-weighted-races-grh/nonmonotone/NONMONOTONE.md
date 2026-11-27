# The weighted bias density is NOT monotone for the character of conductor 163

Computer-assisted proof, 2026-11-24. Code, raw outputs and the assembled
`results.txt` are in this directory. Nothing here is committed to the paper.
Independently audited (`../nonmonotone-audit/AUDIT.md`: no rigor break; one gap,
G1, closed in this revision — see §7, Review history).

**Verdict: the floating-point finding of `BRIDGE-ATTEMPT.md` §4.3 is TRUE under
rigorous evaluation.** The certified enclosures contain the float values
(`0.59033`, `0.59603`), and they are separated by `5.7·10⁻³` while each is
`4·10⁻⁸` wide.

---------------------------------------------------------------------------

## 1. Statement

Notation as in `paper/weighted-races.tex` and `CORE-NOTES.md`. `χ = χ_{−163}` is
the Kronecker symbol `(−163/·)`: the real primitive character of conductor 163,
odd (`χ(−1) = −1`, `a = 1`). `M = m + 1/2 > 0`. The positive ordinates `γ` of
`L(s,χ)` are listed with multiplicity,

    c_γ(M) = 4M / sqrt(M² + γ²),      X_M = 1 + Σ_{γ>0} c_γ(M) cos θ_γ     (θ_γ i.i.d. uniform),
    δ(M) := P(X_M > 0).

**Theorem 1 (model form).** Assume GRH for `L(s, χ_{−163})`. Then

    0.590328985850 ≤ δ(3/5) ≤ 0.590329024804   <   0.596034510065 ≤ δ(6/5) ≤ 0.596034556476 ,

so `δ(6/5) − δ(3/5) ∈ [0.005705485, 0.005705571]`. More completely,

| `M` | `m = M − 1/2` | certified enclosure of `δ(M)` | width |
|---|---|---|---|
| 3/10 | −0.20 | [0.599269990566, 0.599270019978] | 2.9e−8 |
| 57/100 | 0.07 | [0.590250631768, 0.590250669413] | 3.8e−8 |
| 3/5 | 0.10 | [0.590328985850, 0.590329024804] | 3.9e−8 |
| 117/100 | 0.67 | [0.596052042574, 0.596052089770] | 4.7e−8 |
| 6/5 | 0.70 | [0.596034510065, 0.596034556476] | 4.6e−8 |
| 2 | 1.50 | [0.586642408419, 0.586642448867] | 4.0e−8 |

hence

    δ(0.3)  >  δ(0.57)  <  δ(0.6)  <  δ(1.17)  >  δ(1.2)  >  δ(2) ,

the largest certified rise being `δ(1.17) − δ(0.57) ∈ [0.005801373, 0.005801459]`.
(`0.57` and `1.17` were chosen by a float scan as close to the local minimum and
maximum; nothing is claimed about the exact location of the extrema.)
The proof is the certificate of §3. Besides GRH and the trusted software its
only premise is hypothesis (H-Z) of §2 on the two zero files, which the
programs enforce by sha256.

**Theorem 1′ (race form).** Assume GRH and LI for `L(s, χ_{−163})`. Let
`δ_χ(m)` be the logarithmic density of
`{x ≥ 2 : −Σ_{p≤x} χ(p) p^m log p > 0}` (it exists and equals `δ(m + 1/2)` by
Theorem `thm:interp`(ii) of the paper, whose third hypothesis
`L(1/2,χ) ≠ 0` is certified here: `L(1/2, χ_{−163}) = 0.0685317582960750993…`).
Then

    δ_χ(m_1) < δ_χ(m_2)    for  m_1 = 1/10 < m_2 = 7/10        (M_1 = 3/5, M_2 = 6/5),

hence `m ↦ δ_χ(m)` is **not monotone** on `(−1/2, ∞)`: it is neither
non-increasing (`δ_χ(0.1) < δ_χ(0.7)`) nor non-decreasing
(`δ_χ(−0.2) > δ_χ(0.1)`, `δ_χ(0.7) > δ_χ(1.5)`). The monotonicity proved /
conjectured for `q = 3, 4` is therefore not a general fact about real
primitive characters.

**Theorem 2 (conductor 232, same hypotheses for `χ_{−232} = (−232/·)`, odd).**

| `M` | `m` | certified enclosure of `δ(M)` | width |
|---|---|---|---|
| 2/5 | −0.10 | [0.606642556775, 0.606642613806] | 5.7e−8 |
| 71/100 | 0.21 | [0.597435124908, 0.597435193395] | 6.9e−8 |
| 18/25 | 0.22 | [0.597438040240, 0.597438108755] | 6.9e−8 |
| 1 | 0.50 | [0.598231098965, 0.598231174771] | 7.6e−8 |
| 101/100 | 0.51 | [0.598233617676, 0.598233693573] | 7.6e−8 |
| 3/2 | 1.00 | [0.594070789210, 0.594070852885] | 6.4e−8 |

so `δ(0.4) > δ(0.71) < δ(0.72) < δ(1) < δ(1.01) > δ(1.5)`, with
`δ(1) − δ(18/25) ∈ [0.000792990, 0.000793135]` and
`δ(1.01) − δ(0.71) ∈ [0.000798424, 0.000798569]`.
(`L(1/2, χ_{−232}) = 0.24629913501197558…` certified `> 0`.)

**Proposition 3 (the coupling threshold).** Under GRH,

    0.16441017 < M_0(χ_{−163}) < 0.16441035 ,        0.2511058 < M_0(χ_{−232}) < 0.2511076 ,

where `M_0` is the root of `Q(u) = Σ_γ 4u³/(u²+γ²)^{3/2} = 1`. So both
increases occur above `M_0`, consistent with Theorem B of `CORE-NOTES.md`
(strict decrease on `(0, M_0]`), and Theorem B cannot be extended beyond
`M = 0.6 ≈ 3.65·M_0` for conductor 163 (already `δ(0.57) < δ(0.6)`).

---------------------------------------------------------------------------

## 2. Hypotheses, exactly

| used for | hypothesis | status |
|---|---|---|
| ordinates real; variance identity `Σ_γ c_γ²/2 = 4M (ξ'/ξ)(M+1/2)` | GRH(χ) | **assumed** |
| the race has a limiting distribution equal to the law of `X_M` (and so `δ_χ(m) = δ(M)`) | LI(χ) | **assumed** (Theorem 1′ only; Theorem 1 does not use it) |
| Hadamard product without a central zero; `γ_1 > 0` | `L(1/2,χ) ≠ 0` | **certified** (Arb ball `> 0`; also follows from the zero count) |
| the listed zeros are all the zeros with `0 < γ < T_cert`, each simple | **(H-Z)** below: the zero files are the certified ones | **certified** (§3.1), unconditionally, *for the files with the pinned hashes*; enforced by the programs |
| root number `+1`, zeros of `ξ` in the critical strip, `ξ` entire of order 1 | classical | cited, not re-proved |
| `|J_0(x)| ≤ x^{−1/2}` (Watson §13.74); `0 < J_0(x) ≤ e^{−x²/4}` and `−log J_0(x) − x²/4 ≤ c_4 x⁴` on `[0,1]`; Gauss–Legendre remainder; Fourier inversion (paper Lemma `lem:invert`) | classical / paper | cited |
| correctness of FLINT/Arb 3.6.0 through python-flint 0.9.0 (`acb_dirichlet_l`, `acb_dirichlet_hardy_z`, the L-series jet, `lgamma`, `digamma`, `J_0`, ball arithmetic) | software | trusted |

**Hypothesis (H-Z) (zero files).** The enclosures of §1 are statements about the
two zero files

    zeros_q163.txt   sha256 b9d38aba6df821dea401035ea6ae24aa045305e757f81ea1a6233b632ad25b9d   (133 windows, T_cert = 7617/64)
    zeros_q232.txt   sha256 42f40b29d79cc5228478acfd0480f23939c6f3b5d0479acaf178c661fa0af378   (114 windows, T_cert = 6351/64)

and they are theorems only because *for exactly these files* the certificate of
§3.1 passes: every window `[d − 10⁻²⁰, d + 10⁻²⁰]` holds exactly one simple zero
ordinate and there is no other zero with `0 < γ < T_cert`. A zero file is just a
list of decimals; `density.py` cannot see from the list that nothing is missing
(with `γ_1` deleted the arithmetic goes through and produces a false "enclosure"
`δ(3/5) ≈ 0.63504`: audit mutant m1). (H-Z) is discharged three times over —
`zeros_cert.py` in generate mode, `zeros_cert.py --verify` on the stored files,
and the audit's independent re-certification of the windows and of the counts
266 / 228 with a different contour — and it is **enforced mechanically** (§3.1,
"Binding"): `density.py` and `m0.py` refuse any zero file whose sha256 is not
recorded in a passing certificate record, every `DELTA` line prints the hash,
and `make_results.py` checks that the printed hashes are the two above.

**The variance identity for `χ ≠ χ_4`.** The paper proves `prop:varid` for `χ_4`
and states (`eq:varidq`) that the proof applies verbatim to a real primitive
`χ`. What the extension needs is exactly: `ξ(s,χ)` entire of order 1 and real on
ℝ; root number `+1`, i.e. `τ(χ) = i√q` for odd real primitive `χ` (Gauss;
verified numerically for `q = 163, 232` in the audit); GRH; and
`L(1/2,χ) ≠ 0` (certified here). The audit re-derived the extension from the
Hadamard product (`AUDIT.md` §2.3) and found it correct; it deserves one
sentence in the paper if this result is cited.

LI forces the zeros to be simple; the tail bounds below never use simplicity of
the *unlisted* zeros (a zero of multiplicity `k` is `k` factors `J_0` and
contributes `k` times to the variance identity), and the listed ones are
certified simple, so the multiplicity convention of `prop:varid` is matched
exactly.

---------------------------------------------------------------------------

## 3. The certificate

### 3.1 Zero data (`zeros_cert.py`, outputs `zeros_q*.txt`, `zeros_q*.cert.txt`)

* **Character.** `flint.dirichlet_char(163, 162)` (resp. `(232, 115)`) is checked
  to be real, primitive, of conductor `q`, odd, and equal to the Kronecker
  symbol `(−q/n)` for all `n ≤ 2q` (exact integer arithmetic).
* **Localisation.** Each listed decimal `d` has certified opposite signs of the
  Hardy function `Z(t,χ)` at `d ± 10⁻²⁰` (192-bit balls). 133 windows for
  conductor 163 (`γ_1 = 0.20290133749879861298…`, `γ_2 = 2.36853394678810445…`,
  `γ_133 = 118.218287414…`), 114 for conductor 232
  (`γ_1 = 0.31067616665127619…`).
* **Completeness and simplicity.** Argument principle for the completed
  function `Λ(s) = (q/π)^{(s+1)/2} Γ((s+1)/2) L(s,χ)` on
  `R_T = (−1,2) × (−T,T)`. By `Λ(s) = Λ(1−s)` and `Λ(s̄) = \overline{Λ(s)}` the
  four quarters of `∂R_T` contribute equally, so
  `#{zeros of Λ in R_T} = (2/π)·Δ_C arg Λ`, `C: 2 → 2+iT → 1/2+iT`, and
  `Δ_C arg Λ = θ(T) + Arg L(2+iT) + Δ_{horizontal} arg L` with
  `θ(T) = (T/2) log(q/π) + Im logΓ(3/4 + iT/2)` (principal `logΓ` = continuous
  argument), the *principal* `Arg L(2+iT)` (because `|L(2+it) − 1| ≤ ζ(2) − 1 < 1`
  keeps `Re L > 0` on the whole vertical side), and the horizontal part obtained
  by covering `[1/2, 2] + iT` with segments on each of which the Arb image of a
  ball containing the segment excludes `0`. Result (every block):

  | conductor | `T_cert` | windows `K` | certified count in `R_T` |
  |---|---|---|---|
  | 163 | 7617/64 = 119.015625 | 133 | `266 = 2K` (ball `266.000000000 ± 3e−14`) |
  | 232 | 6351/64 = 99.234375 | 114 | `228 = 2K` |

  (intermediate checkpoints near 20, 40, 60, 80, 100 agree: 34, 76, 120, 168,
  218 for 163). Since each window and its mirror image hold at least one zero,
  `count = 2K` proves: one **simple** zero per window, **no other zero** of `Λ`
  with `|Im s| < T_cert` — none off the line, none real (so `L(1/2,χ) ≠ 0`, no
  exceptional real zero), none missed by the scan. This part is unconditional.
* **Binding of the file to the certificate (audit G1).** Whenever windows and
  count have passed for a zero file, `zeros_cert.py` appends to the certificate
  log a line
  `HASH-RECORD sha256=<hash of the zero file> file=… q=… K=… T_cert=p/d W=2K mode=generate|verify status=PASS`.
  `common.load_certified_zeros` — the only loader used by `density.py`, `m0.py`
  and `crosscheck.py` — computes the sha256 of the file it is given and refuses
  (exit status 2, no output line) unless that hash has a passing record whose
  `K`, `T_cert`, `W = 2K` are consistent with the file.
  `zeros_cert.py q --verify [--file F] [--T p/d]` re-certifies an *existing*
  file (all windows at 192 bits, positivity, disjointness, order, `< T_cert`,
  `L(1/2) > 0`, and the count at `T_cert`) without scanning or modifying it, and
  writes the record only if everything passes (15 s for conductor 163, 6 s for
  232). Generate mode reproduces `zeros_q232.txt` byte for byte (same hash; also
  observed by the audit). The record is a consistency lock against stale,
  truncated or edited zero files, not a cryptographic signature — the log is a
  text file — which is why the two hashes are also pinned in (H-Z), in
  `common.HZ_SHA256` and in `results.txt`.
* **Second, heuristic guard.** `density.py` and `m0.py` also refuse when the
  exact tail variance `V` differs from the Riemann–von Mangoldt estimate of the
  unlisted zeros by more than three top-zero contributions plus 1%
  (`rvm_guard`; actual agreement 0.02–0.6% in all runs). With `γ_1` deleted `V = 7.21`
  against `0.0349` and the run is refused even when the hash requirement is
  waived. The guard is float arithmetic, can only refuse, and does **not** see a
  single missing high zero (mutant m2, `γ_67`); that is what the hash binding is
  for.
* **The `χ_4`, `χ_3` lists** of `prototype/explore` used in the sanity runs were
  re-certified by the same `--verify` mode at `T = 640.165707` and `700.024747`
  (`zeros_q4.cert.txt`, `zeros_q3.cert.txt`: `W = 1022 = 2·511`, `1074 = 2·537`,
  all windows with radius `10⁻²⁰`), so the radius `10⁻²⁰` that the loader gives
  every file is certified for them too (audit M3).

### 3.2 The density (`density.py`)

For rational `M`, with `b_j = c_{γ_j}(M)` for the listed zeros and
`c_top = b_last` (every unlisted zero has `c_γ < c_top`):

1. **Gil–Pelaez** (paper Lemma `lem:invert`):
   `δ − 1/2 = (1/π) ∫_0^∞ Φ(t) (sin t / t) dt`, `Φ = Π_γ J_0(c_γ t)`.
2. **Exact tail variance.** `V := σ²(M) − Σ_listed b_j²/2`, with
   `σ²(M) = 4M [ ½ log(q/π) + ½ ψ((M + 3/2)/2) + (L'/L)(M + 1/2) ]` (GRH;
   `prop:varid`/`eq:varidq`, `a = 1`). `L'` is enclosed by the complex step
   `Im L(s + i2⁻⁶⁰)/2⁻⁶⁰` with the Cauchy remainder from
   `|L(z)| ≤ S_max |z| / Re z` (`S_max = max_x |Σ_{n≤x} χ(n)|`, computed exactly:
   9 for conductor 163, 8 for 232), and must overlap Arb's independent power-series
   jet of `L`. `V > 0` is required (it is `0.0349` and `0.1395` at `M = 0.6, 1.2`).
   `V` is exactly `Σ_unlisted c_γ²/2`, so `Σ_unlisted c_γ⁴ ≤ c_top²·2V`.
3. **Two-sided tail model.** For `0 ≤ t ≤ 1/c_top`:
   `e^{−Vt²/2 − c_b t⁴} ≤ Π_unlisted J_0(c_γ t) ≤ e^{−Vt²/2}`, `c_b = 0.018·c_top²·2V`
   (`0.018 ≥ −log J_0(1) − 1/4`, re-checked in Arb). No statistical assumption
   on the unlisted zeros is made.
4. **Main integral** on `[0, T]`, `T ≤ 1/c_top`: integrand
   `P(t) e^{−Vt²/2} (1 + e^{−c_b t⁴})/2 · sin t/t` with `P` the product of the 133
   certified factors; composite 16-point Gauss–Legendre (nodes and weights
   isolated in Arb), 192-bit balls, with the remainder
   `N h (h/r)^{32} Q_16 · sup` where the sup over the complex `r`-neighbourhood is
   bounded by `exp(r Σ b_j + V r²/2 + r)(1 + e^{8 c_b r⁴})/2`.
5. **Model error** `E1 = ∫_0^T |P| e^{−Vt²/2} (1 − e^{−c_b t⁴})/2 · |sin t / t| dt`:
   upper Riemann sums of Arb range enclosures on adaptively bisected cells
   (first 40 factors as balls; the remaining factors bounded by 1, or by
   `e^{−b_j²t²/4}` on cells where `b_j t ≤ 1`; `|sin t/t| ≤ min(1, 1/t)`).
6. **Truncation** `E2 = |∫_T^∞ Φ sin t/t| ≤ Π_{j≤K} b_j^{−1/2} · T^{−K/2}/(K/2)`
   from `|J_0(x)| ≤ x^{−1/2}`; valid for every `K`, optimised over `K`.

   *The regime here differs from `q = 3, 4`*: `c_1 = 3.79` (`M = 0.6`) and
   `3.94` (`M = 1.2`) are close to the supremum 4, so `Φ` is dominated at
   moderate `t` by the slowly decaying, oscillating factor `J_0(c_1 t)`
   (frequency `c_1 ≈ 3.9`, resolved by cells of length `h ≤ 0.094`). With only
   the few large amplitudes the envelope would be useless
   (`K = 3`: a multiple of `∫_T^∞ t^{−5/2} dt`, about `4·10⁻³` at `T = 14`); the bound actually used
   takes `K = 72` envelopes at `T = 28.1` (`M = 0.6`) and `T = 14.2` (`M = 1.2`),
   both below `1/c_top = 49.3` resp. `24.6`, giving `E2 ≤ 3.2·10⁻¹⁵`. Of these 72
   factors, 70 have `b_j T > 1` (`b_72 T = 0.976` resp. `0.988`; an earlier
   version of this note wrongly said "all"). This is harmless:
   `|J_0(x)| ≤ x^{−1/2}` holds for every `x > 0`, the last two factors are merely
   bounded by slightly more than 1, and `K` is chosen to minimise the resulting
   bound (the count of factors above 1 is printed in each `guard` line). Convergence of the Gil–Pelaez integral itself needs
   only `K ≥ 1`.

Budget at `M = 3/5` / `6/5` (raw integral units; the density error is the sum
divided by `π`): `E1 = 6.1e−8 / 7.3e−8`, `E2 = 3.1e−15`, Gauss remainder
`1.4e−21 / 2.2e−21`, ball radius of the quadrature sum `2.4e−19 / 4.8e−21`.
The width of every enclosure is thus the model error `E1`, i.e. the price of
knowing the zeros only up to height 119; it is `10⁵` times smaller than the gap.

The final comparison `upper(δ(3/5)) < lower(δ(6/5))` is made on the printed,
outward-rounded 12-decimal endpoints in exact rational arithmetic
(`make_results.py`).

### 3.3 `M_0` (`m0.py`)

`Q_list(u) ≤ Q(u) ≤ Q_list(u) + V(u)·u / (2 sqrt(u² + γ_last²))` (each unlisted
term is `(c_γ²/2)·u/(2 sqrt(u²+γ²))`), `Q` strictly increasing; the bracket is
certified by two Arb evaluations.

---------------------------------------------------------------------------

## 4. Cross-checks (none of them is part of the proof)

| check | result |
|---|---|
| same code on `χ_4`, `m = 0, 1` (511 certified zeros of `prototype/explore`) | `[0.995927912125, 0.995927938565]`, `[0.797630385848, 0.797630412326]`: inside the published enclosures `[0.995927866782, 0.995927983908]`, `[0.797630378395, 0.797630419780]`; round to `0.995928`, `0.797630` |
| stability: conductor 163 with only the first 60 zeros | `[0.590328859673, 0.590329153124]`, `[0.596034403752, 0.596034753220]`: contain the 133-zero enclosures |
| mpmath (own `L`, `L'`, tanh–sinh quadrature, Gaussian surrogate tail) | `δ(0.6) = 0.590329004992`, `δ(1.2) = 0.5960345189881`, `δ(0.57) = 0.5902506512188`, `δ(1.17) = 0.5960520510988`: all inside |
| numpy/scipy composite Simpson, `h = 5·10⁻⁴` | same values to 12 digits: all inside |
| mpmath re-location of the first 12 ordinates | agree with the Arb windows to `1e−25` |
| mpmath closed-form `σ²` | agrees with the Arb ball to 15 digits; unlisted variance agrees with the Riemann–von Mangoldt estimate to 0.1% (`0.034875` vs `0.034903`) |
| Monte Carlo, 8·10⁶ samples, common phases, `N(0,V)` surrogate tail | `δ(0.6) = 0.590234 ± 0.000187`, `δ(1.2) = 0.596057 ± 0.000190`, paired difference `+0.005823 ± 0.000092` (certified: `0.0057055`) |
| conductor 232: mpmath / Simpson | `0.597438059916`, `0.598231110576`, `0.597435145332`, `0.598233629374`: all inside |
| conductor 232: Monte Carlo, 1.6·10⁷ samples | `δ(1) − δ(0.72) = +0.000821 ± 0.000044` (certified: `0.000793`) |
| earlier floats (`bridge/lowzero_real_q163.out`) | `0.59032886`, `0.59603454`: within `1.5·10⁻⁷` of the enclosures (84 zeros + surrogate there) |
| the zero-count routine of `zeros_cert.py` run on `χ_4` at `T = 640.165707` and `χ_3` at `T = 700.024747` (`count_crosscheck.py`) | `1022 = 2·511` and `1074 = 2·537`: exactly the counts certified by the different, full-rectangle winding method of `certify_arb.py`; also `100 = 2·50` between `γ_50` and `γ_51` for both. (The full-rectangle routine itself is impractical for conductor 163: it did not finish in 13 minutes even for the box up to `T = 3.2`; hence the symmetric count.) |

| mutated zero files (`raw_mutant_refusal.txt`): `γ_1` deleted, `γ_67` deleted | `density.py`: refused, exit 2, no `DELTA` line (hash has no certificate record); with the hash requirement waived, the `γ_1` mutant is refused by the RvM guard (`V = 7.2139` vs `0.0349`); `zeros_cert.py --verify`: `VERIFY FAILED: W = 266 but 2K = 264`, no record written |
| independent audit (`../nonmonotone-audit/`) | own rigorous enclosures (Arb `acb.integral`, one-sided tail model, 133/100 and 114/80 zeros, own zero count on `(−2,3)×(−T,T)`) contain all twelve enclosures of §1 and prove the inequalities on their own |

Limits of the float checks (audit M4). The mpmath value integrates only over
`[0,16]`, and mpmath, Simpson and Monte Carlo all replace the unlisted zeros by
the *Gaussian* `e^{−Vt²/2}`. The test "value lies in the certified enclosure" is
meaningful only because that Gaussian is the upper end of the certified
two-sided model on `0 ≤ t ≤ 1/c_top` (beyond it the integrand is below `10⁻¹⁵`),
so these values differ from the truth by at most `E1/π ≤ 4·10⁻⁸`. All three read
the same zero file as `density.py` and therefore **cannot detect a wrong zero
file** (on the `γ_1` mutant Simpson agrees with the false value): only the hash
binding and the count can.

Not evidence (audit M5): the remark that the full-rectangle routine is
impractical for conductor 163, and all timings in §6, are operational notes. The
symmetric count is justified by the argument of §3.1 (re-derived in `AUDIT.md`
§1.1), not by the other method being slow.

---------------------------------------------------------------------------

## 5. What is and is not rigorous

**Rigorous** (modulo §2: GRH, (H-Z), the cited classical facts, FLINT/Arb): the zero
windows, their completeness and simplicity up to `T_cert`, `L(1/2,χ) > 0`, the
values of `σ²`, every `δ` enclosure in §1, the inequalities between them, the
brackets for `M_0`.

**Assumed, not proved:** GRH for `L(s,χ_{−163})` (resp. `χ_{−232}`); LI for the
race interpretation. Without LI, Theorem 1 is still a theorem about the random
series `X_M` with independent phases.

**Imported:** the limiting-distribution theorem `thm:interp` of the paper (not
re-verified here or in the audit); the variance identity `prop:varid`/`eq:varidq`
(proved in the paper for `χ_4` only; the extension to these characters was
re-derived in the audit, see §2); Theorem B of `CORE-NOTES.md` (only used for the
remark on `M_0`).

**Not claimed:** the location of the local minimum/maximum, the exact interval
on which `δ' > 0` (floats: about `[0.57, 1.17]` for 163, `[0.71, 1.01]` for 232),
or the sign of `δ'` at any point — only inequalities between `δ` at the listed
points. Nothing is claimed for other characters (conductor 403 was "too
marginal to call" in floats and has not been examined here).

**Caveats.**
* The code was written in one session. It has since had one independent audit
  (§7), which re-derived every inequality and re-computed every enclosure with
  different code; the revision that closes G1 (hash binding, `--verify`, RvM
  guard) has **not** itself been re-audited. Its internal
  protection against bugs is: the `χ_4` reproduction, the 60-zero stability
  run, the agreement with two independent float quadratures to `10⁻¹¹`, the
  Monte Carlo, the integer zero counts (`W = 2K` to `3·10⁻¹⁴` at eleven
  checkpoints, with `K` equal to the count found independently in
  `bridge/zeros_q163.txt`: 84 below height 80), and the reproduction by the
  same count routine of the certified counts 511 and 537 for `χ_4`, `χ_3`.
* The enclosures are fail-closed: any uncertified sign, non-finite ball,
  non-positive tail variance, out-of-range truncation, zero file without a
  passing certificate record, or tail variance far from the RvM estimate ends
  the run with exit status 2 and no `DELTA` line. `--uncertified-zeros` waives
  the hash requirement only; its output lines start with `UNCERTIFIED`, which
  `make_results.py` ignores.
* The Gauss nodes are seeded by mpmath and then isolated in Arb; no other
  floating-point number enters a certified quantity (floats choose the scan
  grid, the truncation point `T`, the cells to bisect and the points `M`).
* Runs are deterministic. Outputs printed with `%.2e` (`E1`, widths) are
  informational; the certified objects are the `lower=`/`upper=` endpoints.

---------------------------------------------------------------------------

## 6. Files and reproduction

All commands single-core, from this directory; times are wall-clock on a loaded
4-core machine under `nice -n 10`.

| command | output | time |
|---|---|---|
| `python3 zeros_cert.py 163 120` | `zeros_q163.txt`, `zeros_q163.cert.txt` (resumable, `zeros_q163.state.json`) | 5 min |
| `python3 zeros_cert.py 232 100` | `zeros_q232.txt`, `zeros_q232.cert.txt` | 2 min |
| `python3 zeros_cert.py 163 --verify` ; `python3 zeros_cert.py 232 --verify` | re-certify the stored files; append `HASH-RECORD` to the certificate logs | 15 s, 6 s |
| `python3 zeros_cert.py 4 --verify --file ../../prototype/explore/chi4_zeros.txt --T 640165707/1000000` (and `3 … chi3_zeros.txt --T 700024747/1000000`) | `zeros_q4.cert.txt`, `zeros_q3.cert.txt` (needed by the sanity run) | 1 s each |
| `python3 density.py --q 163 --M 3/10 57/100 3/5 117/100 6/5 2 --out raw_density_q163.txt` | six enclosures | 1.5 min (3–26 s each) |
| `python3 density.py --q 232 --M 2/5 71/100 18/25 1 101/100 3/2 --out raw_density_q232.txt` | six enclosures | 1.6 min |
| `python3 m0.py 163 --out raw_m0.txt` ; `python3 m0.py 232 --out raw_m0.txt` | `M_0` brackets | 2 s each |
| `python3 density.py --q 4 --M 1/2 3/2 --out raw_sanity_chi4.txt` | sanity | 1.1 min |
| `python3 density.py --q 163 --M 3/5 6/5 --nzeros 60 --out raw_stability_q163.txt` | stability | 0.5 min |
| `python3 crosscheck.py 163 3/5 6/5 57/100 117/100 --mc 8e6 --out raw_crosscheck_q163.txt` | mpmath, Simpson, Monte Carlo | 1.8 min |
| `python3 crosscheck.py 232 18/25 1 71/100 101/100 --mc 1.6e7 --out raw_crosscheck_q232.txt` | idem | 2.6 min |
| `python3 count_crosscheck.py --out raw_count_crosscheck.txt` | test of the count routine on `χ_4`, `χ_3` | 1 s |
| `python3 make_results.py` | `results.txt` (exact comparison of endpoints; exit status 0 iff all checks pass) | 1 s |

`common.py`: character selection, exact Kronecker symbol, exact partial sums,
the certified zero-file loader and the pinned hashes of (H-Z).
`raw_mutant_refusal.txt`: transcript of the refusals on the two mutated zero
files (the mutants themselves are not kept: delete `γ_1`, resp. the 67th
ordinate, from a copy of `zeros_q163.txt` and pass it with `--zeros`). The
`--out` files are appended to, not overwritten: delete the `raw_*.txt` files
before a clean reproduction. The first run of `zeros_cert.py` in generate mode
predates the hash records; the records in the stored logs come from `--verify`
(generate mode now writes one per block).

Total for the certificate proper (zeros + densities + `M_0`, both conductors):
about 10 minutes of single-core time.

---------------------------------------------------------------------------

## 7. Review history

**v1 (2026-11-24).** Zero certificate, twelve enclosures, `M_0` brackets,
cross-checks, as described above without the binding of §3.1.

**Independent audit (2026-11-24, `../nonmonotone-audit/AUDIT.md`).** Verdict: no
rigor break. The auditor re-derived the symmetric zero count, the tail model,
the variance identity for these characters, the Gauss–Legendre remainder and
the truncation bound; re-certified all 133 + 114 windows; reproduced the counts
266 and 228 on a different rectangle with a different criterion (plus `W = 2`
at `T = 1/2`: no real zero, one zero pair); and computed independent rigorous
enclosures (Arb's `acb.integral`, a one-sided tail model, different numbers of
zeros and truncation points) that contain all twelve enclosures of §1 and on
their own prove every inequality stated there. Findings and their disposition:

| item | finding | disposition in v2 |
|---|---|---|
| **G1** (gap) | completeness of the zero list tied to the density certificate only by convention: any zero file was accepted; with `γ_1` deleted `density.py` printed false "certified" values, exit status 0 | **fixed**: `HASH-RECORD` in the certificate log; `load_certified_zeros` refuses unrecorded files (used by `density.py`, `m0.py`, `crosscheck.py`); hash in every `DELTA` line, in the `M_0` lines and in `results.txt`; `zeros_cert.py --verify`; RvM guard on the tail variance; (H-Z) stated in §2. The `γ_1`- and `γ_67`-deleted mutants are now refused (`raw_mutant_refusal.txt`) |
| M1 | variance identity listed as "imported" although the paper proves it only for `χ_4` | sentence added in §2 and §5; extension re-derived in the audit (§2.3 there) |
| M2 | "all `K = 72` envelopes have `b_j T > 1`" is false: 70 do | corrected in §3.2 item 6; harmless (`|J_0(x)| ≤ x^{−1/2}` for all `x > 0`); the count is now printed per run |
| M3 | every zero file is given radius `10⁻²⁰`, including the `χ_4`/`χ_3` lists (sanity runs only) | those lists re-certified with `10⁻²⁰` windows by `--verify`; the loader now requires that record |
| M4 | float cross-checks integrate to `t = 16` / use the Gaussian surrogate; containment test valid only because the Gaussian is the upper end of the model | stated in §4 and in `crosscheck.py`; also stated that they cannot detect a wrong zero file |
| M5 | "full-rectangle routine impractical" and the timings are not evidence | stated in §4 |
| M6 | (auditor's own runs) numpy used several threads | `crosscheck.py` now pins the BLAS thread count to 1 |

**v2 (2026-11-24, this revision).** After the changes all twelve densities were
recomputed: the twelve `DELTA` lines are identical to those of v1 character for
character, apart from the `time=` field and the added `zeros_sha256=` field
(the hashes verified here equal those pinned by the audit). The `M_0`
brackets, the `χ_4` sanity enclosures, the 60-zero stability enclosures and all
float cross-check values are likewise unchanged. `make_results.py`: 64 checks,
all passed. No theorem, enclosure or inequality of §1 changed. The v2 changes
themselves (hash binding, `--verify`, guard) have not been re-audited.
