# Audit of `monotonicity/nonmonotone/` (non-monotonicity for conductors 163 and 232)

Auditor: independent re-derivation and re-computation, 2026-11-24.
Nothing in `monotonicity/nonmonotone/` was modified; all audit code, logs and
mutants are in this directory. `SHA256SUMS` pins the audited files (zero files,
certificate logs, state files, code, raw outputs, `results.txt`,
`NONMONOTONE.md`) exactly as audited; they were unchanged at the end of the audit.

## 0. Verdict

* **No RIGOR BREAK found.** Every inequality used was re-derived, is used in the
  right direction and inside its range of validity, and is implemented as stated.
* **One GAP** (G1, the same pitfall as in `cert-audit` G2): `density.py` and
  `m0.py` accept any zero file; completeness is tied to the certificate only by
  convention. Closed for the audited files by hypothesis (H-Z) below plus my
  independent re-certification of the windows and of the counts.
* MINOR items in §5.

**What is certified (with (H-Z) and trusting FLINT/Arb 3.6 / python-flint 0.9):**

* *Unconditionally:* for χ = (−163/·) (the flint character `dirichlet_char(163,162)`)
  every listed window holds exactly one zero, the zero is simple, and there are
  no other zeros of Λ(s,χ) with |Im s| < 7617/64 (none off the line, none real;
  so L(1/2,χ) ≠ 0). The same for (−232/·) (`dirichlet_char(232,115)`) below 6351/64.
* *Under GRH(χ) (model form, Theorem 1 / Theorem 2):* the twelve δ(M) enclosures
  of `NONMONOTONE.md` §1 hold, in particular
  δ(3/5) ≤ 0.590329024804 < 0.596034510065 ≤ δ(6/5) for 163 and
  δ(18/25) ≤ 0.597438108755 < 0.598231098965 ≤ δ(1) for 232. Here
  δ(M) = P(1 + Σ_{γ>0} c_γ(M) cos θ_γ > 0), zeros counted with multiplicity
  (each zero of multiplicity k contributes k independent terms).
  My independent enclosures (different code, integrator, tail model, number of
  zeros, truncation point, envelope constant) contain all twelve audited
  enclosures and on their own already prove both inequalities and all the
  increase/decrease comparisons of §1.
* *Under GRH(χ) + LI(χ) (race form, Theorem 1′):* m ↦ δ_χ(m), the logarithmic
  density of {x ≥ 2 : −Σ_{p≤x} χ(p) p^m log p > 0}, is not monotone on (−1/2, ∞).
  This additionally imports the paper's `thm:interp`(ii) (stated there for every
  real primitive χ under GRH + NCZ + LI; NCZ = L(1/2,χ) ≠ 0 is certified) and the
  variance identity `prop:varid`/`eq:varidq` (proved in the paper for χ₄,
  extended in prose; I re-derived the extension, §2.3). No hypothesis on real
  zeros in (1/2,1) is needed (GRH excludes them; in fact the count excludes all
  real zeros unconditionally).

## 1. Zero data (item 1)

### 1.1 The symmetric count, re-derived

Λ(s) = (q/π)^{(s+1)/2} Γ((s+1)/2) L(s,χ), χ odd real primitive. Λ is entire
(the poles of Γ((s+1)/2) at s = −1, −3, … cancel the trivial zeros of L), its
zeros are exactly the non-trivial zeros of L, and
Λ(1−s) = (τ(χ)/(i√q)) Λ(s) with τ(χ) = i√q (Gauss; checked numerically:
τ = 12.76714533…i = i√163 and 15.23154621…i = i√232, `a6_misc.log`), so
Λ(s) = Λ(1−s); also Λ(s̄) = conj Λ(s) (χ real). Hence Λ(1−s̄) = conj Λ(s): |Λ| is
symmetric about Re s = 1/2 and about the real axis.

Rectangle R_T = (−1,2) × (−T,T), boundary traversed counterclockwise from 2:
C₁: 2 → 2+iT → 1/2+iT, C₂: 1/2+iT → −1+iT → −1, C₃: −1 → −1−iT → 1/2−iT,
C₄: 1/2−iT → 2−iT → 2.
* C₂ is the image of reversed C₁ under s ↦ 1−s̄, and Λ(1−s̄) = conj Λ(s):
  Δ_{C₂} arg Λ = −Δ_{−C₁} arg Λ = Δ_{C₁} arg Λ.
* C₃ is the image of C₁ under s ↦ 1−s, Λ(1−s) = Λ(s): same change.
* C₄ is the conjugate of reversed C₂: same change.
So #zeros(R_T) = 4Δ_{C₁}/(2π) = (2/π)Δ_{C₁} arg Λ, provided Λ ≠ 0 on ∂R_T:
on Re s = 2 by the Euler product, on Re s = −1 by reflection, on the top edge
because the certified balls on [1/2,2]+iT exclude 0 (and reflection gives
[−1,1/2]+iT), on the bottom edge by conjugation.

Along C₁: the factor (q/π)^{(s+1)/2}Γ((s+1)/2) has continuous argument
Im[((s+1)/2) log(q/π) + logΓ((s+1)/2)] because (s+1)/2 stays in Re ≥ 3/4 where
the principal logΓ is analytic; it is 0 at s = 2 and θ(T) = (T/2)log(q/π) +
Im logΓ(3/4 + iT/2) at 1/2+iT. I checked that python-flint's `acb.lgamma` is
the continuous log-gamma: Im lgamma(3/4+60i) = 186.0535464… (= mpmath
`loggamma`), while Arg Γ(3/4+60i) = −2.442…. On Re s = 2, |L − 1| ≤ ζ(2) − 1 < 1
for all t, so Re L > 0 on the whole edge, L(2) > 0, and the continuous change
equals the principal Arg L(2+iT). On the horizontal piece the change over a
sub-segment whose Arb image (a rectangle, convex) excludes 0 is the principal
Arg of the ratio of endpoint values (a closed convex set missing 0 lies in an
open half-plane through 0). The code (`Cert.count`, `horizontal_arg`) does
exactly this, at 128 bits, with exact dyadic subdivision points, and requires
the Arg balls to have radius < 1e−6 and W to contain a unique integer.

Conclusion logic: each window [d ∓ 1e−20] has certified opposite signs of Arb's
Hardy Z, hence an odd number (≥ 1) of zeros counted with multiplicity on the
line; the mirror windows likewise. W = 2K forces exactly one simple zero per
window and nothing else in R_T. Correct.

### 1.2 Independent checks (`a1_character_zeros.py`, `a2_count.py`)

* **Character.** flint's `dirichlet_char(163,162)` / `(232,115)` is real,
  primitive, conductor q, odd, and χ(n) (evaluated by flint) equals my own
  Kronecker symbol (−q/n) for all 1 ≤ n ≤ 400 (both q). Arb's `chi.l` agrees with
  mpmath's Hurwitz-based `dirichlet` built from my Kronecker list at s = 2
  (to 20 digits) and at 1/2, 1/2+3.7i, 1/2+15.25i (15 digits). So Arb's L is
  called with the right character.
* **Ordinates.** My own mpmath Hardy-Z scan (step 0.02, own θ) re-locates the
  first 20 ordinates: max deviation from the file 1.7e−25 (163) / 1.6e−25 (232).
  γ₁ = 0.2029013374987986…, 0.3106761666512761….
* **Sign changes.** All 133 (resp. 114) windows re-certified with my own loop
  (192 bits): opposite signs, windows disjoint, increasing, positive, below T_cert.
* **Count, different contour and criterion.** `a2_count.py` uses the rectangle
  (−2,3)×(−T,T) (right edge σ = 3, |L−1| ≤ ζ(3)−1), a disc criterion
  sup|L(ball) − L(c)| < |L(c)| on each sub-segment, and its own bookkeeping:

  | q | T | W (mine) | listed below T |
  |---|---|---|---|
  | 163 | 1/2 | 2 | 1 |
  | 163 | 5914/64 = 92.406 | 200 | 100 |
  | 163 | 7617/64 = 119.016 | 266 | 133 |
  | 232 | 1/2 | 2 | 1 |
  | 232 | 4675/64 = 73.047 | 160 | 80 |
  | 232 | 6351/64 = 99.234 | 228 | 114 |

  The top counts equal the audited ones (266, 228). The T = 1/2 counts show
  directly that there is no real zero and exactly one zero pair below 1/2.
* **Riemann–von Mangoldt.** (T/2π)log(qT/2πe) = 133.26 at T_cert (163) and 113.81
  (232); θ(T)/π = 133.38 / 113.94; listed K = 133 / 114.
* `zeros_q*.txt` equals the zero list in `zeros_q*.state.json`; the block log is
  consistent (W = 2K at every checkpoint).

## 2. Tail and quadrature (item 2), re-derived

### 2.1 Large-t remainder E2
|J₀(x)| ≤ √(2/(πx)) < x^{−1/2} for **all** x > 0 (Watson §13.74: x(J₀²+Y₀²)
increases to 2/π; numerically sup √x|J₀(x)| on (0,2000] = 0.797885 = √(2/π)).
So the envelope is valid for every factor, including those with b_jT < 1 (where
it is merely weaker than |J₀| ≤ 1), and |∫_T^∞ Φ sin t/t| ≤ Π_{j≤K} b_j^{−1/2}
T^{−K/2}/(K/2) for every K. Recomputed: K = 72 is the optimum and E2 = 3.118e−15
(M = 3/5, T = 28.125) and 3.092e−15 (M = 6/5, T = 14.234): matches the printed
3.12e−15 / 3.09e−15. (The text's "all of them with b_jT > 1" is false — see M2 —
but irrelevant to validity.)

### 2.2 Tail model for unlisted zeros
For 0 ≤ x ≤ 1: J₀(x) = Π_k(1 − x²/j_k²) gives −log J₀(x) = Σ_n σ_n x^{2n}/n with
σ_n = Σ_k j_k^{−2n} > 0, σ₁ = 1/4. Hence J₀(x) ≤ e^{−x²/4}, and
(−log J₀(x) − x²/4)/x⁴ is increasing, so ≤ −log J₀(1) − 1/4 = 0.0176211 < 0.018.
Every unlisted γ is ≥ T_cert > γ_last, so c_γ < c_top = b_last, and for
t ≤ T ≤ 1/c_top all arguments c_γt ≤ 1. With V = Σ_unlisted c²/2 (multiplicity
included), Σ c⁴ ≤ c_top²·2V. So e^{−Vt²/2 − 0.018 c_top² 2V t⁴} ≤ Ψ(t) ≤ e^{−Vt²/2}:
correct, and the code enforces T·c_top ≤ 1 in Arb. Unlisted multiple zeros are
harmless (k factors). The midpoint/half-width split and the E1 integrand are
correct; the E1 cell bound uses |J₀| ≤ 1 or e^{−b²t²/4} only where b_j·t ≤ 1 is
certified, and |sin t/t| ≤ 1/a only when a ≥ 1 is certified.

### 2.3 Tail variance
ξ(s) = Λ(s) is entire of order 1, real, ξ(s) = ξ(1−s); Hadamard with
B = −Σ Re 1/ρ gives ξ'/ξ(s) = Σ_ρ 1/(s−ρ) (symmetric summation). Under GRH and
with no central zero, at s = M + 1/2 the pair ρ = 1/2 ± iγ contributes
2M/(M²+γ²), and c_γ²/2 = 8M²/(M²+γ²) = 4M·2M/(M²+γ²); so
Σ_{γ>0} c_γ²/2 = 4M ξ'/ξ(M+1/2) = 4M[½log(q/π) + ½ψ((M+3/2)/2) + L'/L(M+1/2)],
with multiplicity. This is `eq:varidq` with a = 1, m+1 = M+1/2; the code
evaluates exactly that (`digamma` at (s+1)/2 with s = M+1/2). L' from the complex
step: Im L(s+ih)/h − L'(s) involves only odd k ≥ 3, bounded by
M_L h²/(R³(1−(h/R)²)) with M_L = S_max (s+R)/(s−R) from Abel summation — correct;
and it must overlap Arb's jet. My σ² from the jet alone agrees to all printed
digits (8.37716380533145, 12.2068736380107, …), and float mpmath agrees too.

### 2.4 Gauss–Legendre
Remainder h^{2n+1}(n!)⁴/((2n+1)((2n)!)³)·|f^{(2n)}| with Cauchy
|f^{(2n)}| ≤ (2n)! M_b/r^{2n} gives N·M_b·h·Q_n·(h/r)^{2n}, Q_n =
(n!)⁴/((2n+1)((2n)!)²): matches the code. The disc bounds |J₀(bz)| ≤ e^{b|y|},
|e^{−Vz²/2}| ≤ e^{Vy²/2}, |sin z/z| ≤ e^{|y|}, Re z⁴ ≥ −8y⁴ are correct. Nodes are
isolated (sign change of P₁₆ at both ends of a 1e−45 ball, P₁₆' ≠ 0 on it);
weights 2/((1−x²)P'²) enclosed. Gil–Pelaez in the form
δ = 1/2 + (1/π)∫₀^∞ Φ sin t/t holds since Φ is real (symmetric law), integrable
(K ≥ 3 envelopes), and the law of X_M has no atom.

## 3. Independent rigorous re-evaluation (item 3)

`a3_indep_density.py` (my code; nothing imported from the audited directory):
* zero windows: first N lines of the file, sign changes re-certified (§1.2),
  completeness from **my** count at T0 (unlisted γ ≥ T0, c_top = 4M/√(M²+T0²));
* σ² from Arb's L-series jet only;
* main term: `acb.integral` (Arb's rigorous adaptive Gauss–Legendre with analytic
  error bounds) of Π J₀(b_j t)·e^{−Vt²/2}·sinc t on [0,T] (entire integrand);
* my tail model: |Ψ − e^{−Vt²/2}| ≤ e^{−Vt²/2}·C·S4·t⁴ (one-sided around the
  Gaussian, not their midpoint model), C = −log J₀(1) − 1/4 enclosed in Arb;
  error ≤ C·S4·∫₀^T |P| e^{−Vt²/2} t³|sin t| dt by an upper Riemann sum of Arb
  interval enclosures (cells 1/64);
* t-tail with √(2/(πx)) envelopes (different constant), optimised over K;
* T = min(0.95/c_top, 40) on a 1/64 grid (different from theirs).

Results (`a3_runs.log`, `a6_misc.log`):

| q | M | N | mine | audited | inside |
|---|---|---|---|---|---|
| 163 | 3/5 | 133 | [0.590328717617, 0.590329292367] | [0.590328985850, 0.590329024804] | yes |
| 163 | 3/5 | 100 | [0.590328409967, 0.590329599644] | same | yes |
| 163 | 6/5 | 133 | [0.596034368907, 0.596034669069] | [0.596034510065, 0.596034556476] | yes |
| 163 | 6/5 | 100 | [0.596034200377, 0.596034821648] | same | yes |
| 232 | 18/25 | 114 | [0.597437751927, 0.597438367905] | [0.597438040240, 0.597438108755] | yes |
| 232 | 18/25 | 80 | [0.597437299501, 0.597438797993] | same | yes |
| 232 | 1 | 114 | [0.598230879642, 0.598231341511] | [0.598231098965, 0.598231174771] | yes |
| 232 | 1 | 80 | [0.598230528663, 0.598231652204] | same | yes |

and likewise for M = 3/10, 57/100, 117/100, 2 (163) and 2/5, 71/100, 101/100, 3/2
(232): all twelve audited enclosures lie inside mine. My enclosures alone give
δ(6/5) − δ(3/5) ≥ 5.705e−3 (N = 133 and N = 100), δ(1) − δ(18/25) ≥ 7.92e−4
(N = 114 and N = 80), and reproduce every increase/decrease of §1.

Non-rigorous (`a4_plain_mc.py`): plain composite 20-point Gauss on [0,60] with the
Gaussian surrogate: δ(0.6) = 0.590329004992, δ(1.2) = 0.596034518988,
δ(0.72) = 0.597438059916, δ(1) = 0.598231110576 (inside). Monte Carlo, 4·10⁶
common phases: paired δ(1.2) − δ(0.6) = +0.005738 ± 0.000151,
δ(1) − δ(0.72) = +0.000808 ± 0.000085 (certified 0.0057055 / 0.000793).

## 4. Mutation tests (item 4; `a5_mutants.sh`, `a5_mutants.log`, copies in `mutants/`)

| mutant | effect | detected? |
|---|---|---|
| m1: delete γ₁ from `zeros_q163.txt` | density.py prints "certified" δ(3/5) ∈ [0.635042760165, 0.635043389430], δ(6/5) ∈ [0.612136128413, 0.612137162595] — false (truth ≈ 0.590329 / 0.596035), and the order reverses | **No.** exit 0, no error. crosscheck.py reads the same file and would agree with the wrong value (Simpson on the mutant, §6) |
| m2: delete γ₆₇ | [0.590328985432, 0.590329025159] / [0.596034508239, 0.596034555565]: still contains the truth, but the tail model's premise c < c_top is false, so it is not a proof | **No** |
| m3: tail-model sign flip e^{−Vt²/2} → e^{+Vt²/2} | δ(3/5) shifts to 0.589898… (outside the true enclosure) | not by density.py; yes by make_results.py (Simpson/mpmath containment would FAIL) |
| m4: drop E1 from the radius | width 2e−15 around the model midpoint 0.5903290053 | not by density.py; yes by make_results.py (the Gaussian Simpson value 0.590329004992 falls outside) |
| m5: zeros_cert.py "misses" sign changes | `*** COUNT MISMATCH: W = 34 but 2K = 10 … Block NOT accepted.` | **yes** |
| m6: zeros_cert.py θ(T) sign error | `count ball [17.028…] does not determine an integer` | **yes** |

So the count is an effective guard **when zeros_cert.py runs**, but nothing that
consumes the zero file re-checks it (G1).

## 5. Findings

### 5.1 RIGOR BREAK
None.

### 5.2 GAP

**G1 — Completeness of the zero list is tied to the density certificate only by
convention.** `common.load_zeros` (used by `density.py`, `m0.py`) checks only
ordering and positivity. The completeness/simplicity certificate exists only as
the log `zeros_q*.cert.txt`, which records neither a hash of the zero file nor
the file contents, and there is no verify-only mode of `zeros_cert.py`.
Mutation m1 shows that an incomplete file yields a silently wrong "certified"
enclosure (and reverses the claimed inequality); the non-rigorous cross-checks
use the same file and do not notice. *Needed hypothesis:*
**(H-Z)** the zero files are those with
sha256 `b9d38aba6df821dea401035ea6ae24aa045305e757f81ea1a6233b632ad25b9d`
(`zeros_q163.txt`) and `42f40b29d79cc5228478acfd0480f23939c6f3b5d0479acaf178c661fa0af378`
(`zeros_q232.txt`), on which the completeness certificate passes.
For these exact files the audit discharges (H-Z): all windows re-certified
(`a1`), counts 266 / 228 reproduced with a different contour and criterion
(`a2`), and `zeros_cert.py 232 100` re-run from scratch reproduces the file
byte for byte (§6). *Fix:* have `zeros_cert.py` write the sha256 of the zero
file into the cert log, have `load_zeros` require a matching hash (or a
`--verify` mode that re-certifies windows and re-counts at T_cert, ~5 min), and
print the hash in every `DELTA` line.

### 5.3 MINOR

* **M1 (hypothesis table).** `NONMONOTONE.md` §2 lists the variance identity as
  "imported"; the paper proves `prop:varid` for χ₄ only and says it "applies
  verbatim". The extension needs exactly: Λ entire of order 1, root number +1
  (τ(χ) = i√q, Gauss, here also verified numerically), real on ℝ, GRH,
  L(1/2,χ) ≠ 0 (certified). Re-derived in §2.3; correct. Worth one sentence in
  the paper if cited.
* **M2.** §3.2 item 6 says the K = 72 envelopes are "all of them with b_jT > 1":
  only 70 are (b₇₂T = 0.976 at M = 3/5, 0.988 at 6/5). Harmless: |J₀(x)| ≤ x^{−1/2}
  for all x > 0.
* **M3.** `load_zeros` gives every file radius 1e−20, also the χ₄/χ₃ files of
  `prototype/explore` (sanity runs only); the windows certified there are not
  necessarily 1e−20. Not used in the theorems.
* **M4.** Cross-check `(A)` of crosscheck.py integrates only to t = 16 and both
  float checks use the Gaussian surrogate; that is fine for a sanity check, and
  the containment test in make_results.py is valid only because the Gaussian is
  the *upper* end of the two-sided model (true for t ≤ 1/c_top; beyond, the
  integrand is ~1e−15).
* **M5.** The note "(The full-rectangle routine itself is impractical …)" and the
  timing tables are not evidence; the replacement method is sound (§1.1).
* **M6.** numpy used several threads in my float checks despite the one-core
  request (OMP threads not pinned in the first run); results unaffected.

## 6. Provenance (`a7_replay.sh`, `a7_replay.log`)

* `zeros_cert.py 232 100` re-run from scratch in a clean copy (`replay/`):
  `zeros_q232.txt` **byte-identical** to the audited file; the BLOCK lines of
  the certificate log identical up to timings (W = 36, 80, 128, 178, 228).
* `density.py --q 163 --M 3/5 6/5` re-run in the clean copy: the two `DELTA`
  lines identical to `raw_density_q163.txt` up to the `time=` field.
* `crosscheck.py` (Simpson part) on mutant m1 (γ₁ deleted): δ = 0.635043357516
  and 0.612137120882, i.e. it agrees with the wrong m1 "enclosures", so the
  float cross-checks of `make_results.py` would pass on an incomplete file. The
  same run prints V_tail = 7.2139 next to the Riemann–von Mangoldt tail estimate
  0.034903: a cheap automatic guard (require V_tail within, say, a factor 2 of
  the RvM estimate) would have caught m1, but nothing checks it.
* `sha256sum -c SHA256SUMS` in `../nonmonotone`: all files unchanged after the audit.

## 7. Files

* `SHA256SUMS` — audited files as audited.
* `a1_character_zeros.py` / `.log` — character, L-values, lgamma semantics,
  mpmath ordinates, RvM, window re-certification.
* `a2_count.py` / `a2_count.log` — independent zero count (σ = 3 contour).
* `a3_indep_density.py` / `a3_runs.log` — independent rigorous enclosures.
* `a4_plain_mc.py` / `.log` — float evaluation and Monte Carlo.
* `a5_mutants.sh` / `.log`, `mutants/` — mutation tests.
* `a6_misc.py` / `.log` — containment table, Gauss sums, functional equation,
  J₀ envelope, C4, file/state consistency.
* `a7_replay.sh` / `.log`, `replay/` — re-run from a clean copy.
