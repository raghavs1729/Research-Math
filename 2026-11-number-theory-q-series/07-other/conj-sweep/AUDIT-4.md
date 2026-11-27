# AUDIT-4 — adversarial audit of `R10-traffic-496.md` (arXiv:2609.01562, Conjecture `conj:496`)

Date: 26 Nov 2026. Auditor's own programs and logs: `audit4/`. Nothing outside `audit4/` and this file was touched.

## Verdict

**PROOF COMPLETE** (no error found) for what R10 actually states: `ρ(n) < 1` for all `n ≥ 496`, the explicit
("Equivalently, …") form of the conjecture. The verbal form (maximum of `f_B` at `(1,1)` and `(n−1,n−1)` for every
antidiagonal obstruction) follows from it through results that are proved in the source paper; R10 says so explicitly and
does not re-prove them. I checked that the reduction covers every antidiagonal obstruction (§1).

**NOT NOVEL.** A public Lean 4 / Mathlib formalisation of exactly the same inequality was merged on 20 Sep 2026
(the-omega-institute/trureturing, issue #8249, PR #8934) — two months before R10. R10's status line
("GitHub issue search for the arXiv id: no hits") is wrong in effect and must be corrected (§4).

Only presentational remarks remain (§2.7); none is a gap.

---

## 1. Source check

Source: `https://arxiv.org/e-print/2609.01562` → gzip tar, sha256 `08473a86…5ec57e`; members `max_traffic.tex`
(61321 bytes, sha256 `50320d41…c55ac7`), `Location300.png`, `00README.json`. Read only (python3 -I + tarfile), in
`/tmp/arx260901562/`. arXiv API: a single version **v1**, 1 Sep 2026, authors J. Gil, Z. Liang, A. Odetola, M. Weiner,
math.CO, "19 pages, 7 figures"; no v2, no journal ref, no withdrawal.

| Item | Source (`max_traffic.tex`) | R10 | Match |
|---|---|---|---|
| `D(n)` | Prop. `prop:antidiag-criterion` (l. 907): `D(n) = (1/n)·C(2n−2,n−1)` | `2C(2n−2,n−1) − C(2n−1,n−1)` | same number; R10 quotes the bracket from the *proof* (l. 912: `f(1,1) − f(1,0)`), not the displayed definition |
| `G(n,a)` | l. 907: `G(n,a) = (n−2a+1)/(n−a)·C(n,a)·C(n−2,a−1)` | `[2 − (n−1)/(n−a)]·C(n,a)·C(n−2,a−1)` | same number; again the proof's bracket form (l. 912) |
| `R_n(a)` | l. 916: `G(n,a)/D(n)` | same | exact |
| `ρ(n)` | l. 917: `max_{1≤a<n/2} R_n(a)` | same | exact |
| Conjecture | l. 940–942, verbatim as quoted in R10 | | exact |
| data | `ρ(n) > 1` for `5 ≤ n ≤ 375`; anomaly for `8 ≤ n ≤ 375`; 82 of 120 values in `376..495`; verified `ρ<1` for `496..2000` | R10: "`ρ(n) > 1` for `8 ≤ n ≤ 375`" | harmless imprecision (source: from 5) |
| `lem:monotone` | l. 919–924: quotient `a(a−1)(n−2a+3)/((n−a+1)²(n−2a+1))`, `<1` when `n−2a ≥ √n` | R10 uses only the quotient, with its own window | exact |

Both forms of `G`, `D` agree exactly for all `2 ≤ n ≤ 400`, `1 ≤ a ≤ n−1` (`audit4/a1.log`).

**What the source assumes for the equivalence.**
* `prop:antidiag-ties` (`n ≥ 2`, `1 ≤ a ≤ n−1`): reflection in the antidiagonal fixes `B`, so `f_B(1,1)=f_B(n−1,n−1)`,
  `f_B(1,0)=f_B(n,n−1)`, `f_B(0,1)=f_B(n−1,n)`; `f_B(0,1)<f_B(1,0)` if `2a<n`, reversed if `2a>n`.
* "By symmetry, we may assume `2a<n`" (l. 902) — the case `2a>n` is the mirror image in `y=x`, which fixes `(1,1)`, `(n−1,n−1)`.
* `prop:antidiag-criterion`: `f_B(1,0) − f_B(1,1) = G(n,a) − D(n)`.
* l. 917: "recalling that for `n ≥ 9` the maximum is attained on `𝒞`" (the six corner-adjacent points) — this is the
  classification of Sections 3–6: unnamed theorem for `B∈R_2` (`n ≥ 8`), `thm:R1` (`n ≥ 9`, `3<c+d<2n−3`), `thm:R0`
  (`n ≥ 5`), `lem:axes` for boundary `B`; the exceptional set `𝒮={(1,2),(2,1),(n−1,n−2),(n−2,n−1)}` meets `x+y=n` only for `n=3`.

**Does `ρ(n) < 1` (max over `1 ≤ a < n/2`) cover every antidiagonal obstruction?** Yes, with the source's results:
* `1 ≤ a < n/2`: `R_n(a)<1` ⇒ `f_B(1,0)<f_B(1,1)`, and `f_B(0,1)<f_B(1,0)`; with the classification the maximum is exactly at `(1,1)`, `(n−1,n−1)`.
* `n/2 < a ≤ n−1`: mirror image.
* `2a = n` (even `n`, `B=(n/2,n/2)` on the diagonal): **not** covered by `ρ`, and not needed — `thm:R0` with `c<n≤2c`
  (true for `n=2c`) gives `f_B(A) ≤ f_B(1,1)`, and `f_B(1,1)=f_B(n−1,n−1)`.
* `a ∈ {0,n}` (`B=(0,n)` or `(n,0)`, also "on the antidiagonal"): `lem:axes`(ii)/(iii); `f_B(1,1)=f(1,1)=f(n−1,n−1)=f_B(n−1,n−1)`.

**Mismatches between R10 and the conjecture as printed.**
1. R10 proves the second sentence of the conjecture; the first sentence then rests on the source's classification
   (Sections 3–6) and Props. 7.1/7.2, which R10 does not re-prove and I did not referee in full. R10 states this
   ("the equivalence itself is the source's"). Not an error, but the claim "PROVED" for the verbal form is conditional on
   the published reductions being correct.
2. The source's word "Equivalently" is slightly loose: form 1 ⇐ `ρ(n)<1`; conversely form 1 only gives `ρ(n) ≤ 1` (a tie
   `ρ(n)=1` would still have the maximum attained at `(1,1)`). Irrelevant here: R10 proves the stronger strict statement,
   and `R_n(a)=1` never occurs for `3 ≤ n ≤ 16100` (exact).

**Independent brute-force support for the reduction** (my own DP over `B`-avoiding paths, no binomial formulas, no
theorem of the source; `audit4/a1.log`, `audit4/a6.log`):
* `f_B(1,0) − f_B(1,1) = G − D` exactly for all `2 ≤ n ≤ 70`, `1 ≤ a ≤ n−1`.
* For `9 ≤ n ≤ 70`, all `n+1` antidiagonal obstructions, all grid points `A`: the maximum lies in `𝒞`, and form 1 holds
  iff `ρ(n) ≤ 1` (it never holds there).
* Full grid, all `0 ≤ a ≤ n`: form 1 **fails** for `n = 9, 20, 375, 377, 495` (two obstructions each, e.g. `n=495`:
  `B=(240,255)` → maximisers `(1,0),(495,494)`; `B=(255,240)` → `(0,1),(494,495)`), and **holds with maximiser set exactly
  `{(1,1),(n−1,n−1)}` for every antidiagonal `B`** for `n = 376, 423, 496, 497, 498`. So the verbal conjecture is
  confirmed directly at `n = 496, 497, 498`, and the equivalence is confirmed in both directions on these samples.

---

## 2. Proof audit (line by line)

No inequality fails. Every displayed inequality was re-derived by hand, checked symbolically (`audit4/a5_symbolic.py`,
sympy) and tested with rigorous ball arithmetic (python-flint `arb`, 128 bits) for **all** `3 ≤ n ≤ 3000`,
`1 ≤ a < n/2` — 2 248 500 pairs, 0 failures of any step (`audit4/a2_bounds.py`, `a2.log`).

### 2.1 Step 0, formula (1) — correct
`2C(2n−2,n−1) − C(2n−1,n−1) = (2 − (2n−1)/n)·C(2n−2,n−1) = C(2n−2,n−1)/n`;
`C(n−2,a−1) = C(n,a)·a(n−a)/(n(n−1))`; `C(2n−2,n−1) = C(2n,n)·n/(2(2n−1))`. Hence
`G = (u+1)·a·C(n,a)²/(n(n−1))`, `D = C(2n,n)/(2(2n−1))`, and (1) follows. Exact equality `(1) = G/D` verified for all
`2 ≤ n ≤ 400`, `1 ≤ a ≤ n−1` (it holds for all `a`, not only `2a<n`).

### 2.2 (B1) — correct, both inequalities strict
`w_{m+1}/w_m = (2m+1)/(2m+2)` ✓. `(2m+1)²(m+5/4) − (2m+2)²(m+1/4) = 1/4 > 0` ⇒ `w_m²(m+¼)` strictly increasing;
`(2m+1)²(m+3/2) − (2m+2)²(m+½) = −(m+½) < 0` ⇒ `w_m²(m+½)` strictly decreasing; `m·w_m² → 1/π` (Wallis; sympy limit
`= 1/π`), so both tend to `1/π`, and strict monotonicity gives strict inequalities for every `m ≥ 0`. Ball-arithmetic
check of the two-sided bound and of both monotonicities for `1 ≤ m ≤ 6002`: OK.

### 2.3 (B2), (B3) — correct
Even `n=2m`, `a=m−t`, `1 ≤ t ≤ m−1`: `C(2m,m−t)/C(2m,m) = m!²/((m−t)!(m+t)!) = ∏_{i=1}^{t}(m−t+i)/(m+i)` ✓;
`1−x ≤ e^{−x}` and `1/(m+i) ≥ 1/(m+t)` give `≤ exp(−t²/(m+t))` ✓.
Odd `n=2m+1`, `a=m−j`, `0 ≤ j ≤ m−1`, `u=2j+1`: `C(n,m−j)/C(n,m) = m!(m+1)!/((m−j)!(m+1+j)!) =
∏_{i=1}^{j}(m−j+i)/(m+1+i)` ✓ (empty product for `j=0`); `Σ ≥ j(j+1)/(m+1+j)` ✓; `C(2m+1,m) = ½C(2m+2,m+1)` ✓.
Index ranges match `1 ≤ a < n/2` exactly (`t ≥ 1`; `j ≥ 0`).

### 2.4 Inequality (2) — correct for both parities, strict
Even: `√(π(2m+½))/(π(m+¼)) = 2/√(π(2m+½)) = 2/√(π(n+½))` (identity, sympy `= 0` difference); exponent
`2t²/(m+t) = u²/(n+u) ≥ (u²−1)/(n+u)` ✓ (R10 silently weakens `u²` to `u²−1`; legitimate).
Odd: `C(2m+1,m)² < 4·16^m/(π(m+5/4))` (B1 at `m+1`), `C(4m+2,2m+1) > 4·16^m/√(π(2m+3/2))` (B1 at `2m+1`);
`m+5/4 = (n+3/2)/2`, `2m+3/2 = n+½`, ratio `= 2√(n+½)/(√π(n+3/2)) < 2/√(π(n+½))` ⇔ `n+½ < n+3/2` ✓;
exponent `2j(j+1)/(m+1+j) = (u²−1)/(n+u)` ✓. B1 is only used with index `≥ 1`.

### 2.5 Step 2 — correct
* `a/(n−1)·2(2n−1)/n = (n−u)(2n−1)/(n(n−1))`, and `2n(n−1) − (n−u)(2n−1) = 2un − n − u = n(2u−1) − u > 0` for
  `1 ≤ u < n` ✓. So (3) holds strictly.
* `u²−1 = v²−2v`, `n+u = n+v−1` ✓.
  `(v²−2v)/(n+v−1) − [v²/(n+v) − 2v/n] = v(3nv − 2n + 2v² − 2v)/(n(n+v)(n+v−1)) ≥ 0` for every real `v ≥ 1` ✓.
  `v²/(n+v) − [v²/n − v³/n²] = v⁴/(n²(n+v)) ≥ 0` ✓. (Exact rational check on all pairs, `n ≤ 3000`.) Neither step needs
  `v ≥ 2` or `u ≤ n−2`; both hold for all real `v ≥ 1`.
* Case `v ≤ 1.2√n`: `max_v v e^{−v²/n} = √(n/(2e))` at `v=√(n/2)` ✓; `v³/n² ≤ 1.728/√n`, `2v/n ≤ 2.4/√n` ✓.
* Case `u ≥ 1.2√n − 1`: `d/du[(u²−1)/(n+u)] = (u²+2un+1)/(n+u)²` ✓, so
  `(log φ)′ < 0 ⇔ P(u) := (u+1)(u²+2un+1) − (n+u)² > 0`. Since `(u+1)(u²+2un+1) − 2u²n = u³+u²+u+2un+1 > 0`,
  `(n+u)² ≤ 2u²n` is sufficient ✓, i.e. `u ≥ n/(√(2n)−1)` ✓. With `s=√n`:
  `(1.2s−1)(√2 s−1) − s² = 0.697056…s² − 2.614213…s + 1`, larger root `s = 3.317990…`, i.e. `n ≥ 11.00906…`; so the
  threshold `n ≥ 12` is right (at `n=11` the sufficient condition fails by `0.0027`). R10's printed `0.697`, `2.615` are
  rounded in the safe direction. The condition is a lower bound on `u`, so it holds on the whole half-line.
  I also verified `P(1.2√n−1) > 0` by ball arithmetic for all `3 ≤ n ≤ 20000` and `P′(u)=3u²+4un+1>0`.
* **Real-variable monotonicity is legitimate here**: `φ` is a smooth function of real `u>−1`; for an integer
  `u ≥ u₀ := 1.2√n−1`, `φ(u) ≤ φ(u₀)`, and the first-case chain applies at the real point `v₀ = 1.2√n` because the chain
  is valid for all real `v ≥ 1` (see above) and `v₀ ≤ 1.2√n` trivially. The two cases are exhaustive
  (`v ≤ 1.2√n` or `v ≥ 1.2√n`).
* Final line: `√(n/(n+½)) < 1` ✓. `φ(u) ≤ √(n/(2e))·e^{4.128/√n}` holds for all `12 ≤ n ≤ 3000`, all `u` (ball arithmetic).

### 2.6 Steps 3 and 4 — correct
`4/√(2πe) = 0.967882898076573…`, `ln(√(2πe)/4) = 0.0326441720847821…`, `(4.128/that)² = 15990.72152…`;
bound at `n=15990`: `1.00000074` (not `<1`), at `n=15991`: `0.99999972 < 1` ✓. The bound is decreasing in `n`, so all
`n ≥ 15991` are covered; Step 4's range `496..16100` overlaps. The strictness chain is `ρ(n) < bound(n) < 1`; no
`<`/`≤` slip matters anywhere.
Step 4: `R_n(a−1) < R_n(a) ⇔ a(a−1)(n−2a+3) < (n−a+1)²(n−2a+1)` follows from `lem:monotone`'s quotient (re-derived:
`G(n,a−1)/G(n,a) = (n−2a+3)/(n−2a+1) · (n−a)/(n−a+1) · a/(n−a+1) · (a−1)/(n−a)`; exact check `n ≤ 400`) and `R_n>0`.
If `a_min` is the least violator, `R_n` is strictly increasing on `[1, a_min−1]`, so the maximum is on
`[a_min−1, n/2)` ✓ (this window argument does not need unimodality). The session program implements this correctly.

### 2.7 Remarks (presentation only — no gap)
1. R10 l. 12–14: `G`, `D` are quoted in the form appearing inside the source's proof of Prop. 7.2, not as displayed in its
   statement. Identical values; say so.
2. R10 l. 17: the source reports `ρ(n) > 1` for `5 ≤ n ≤ 375` (anomaly from `n = 8`).
3. R10 l. 56–65: the chain is announced for `v = u+1 ≥ 2` (integers) but is then used at the real point `v₀ = 1.2√n`;
   add "valid for all real `v ≥ 1`".
4. R10 l. 75: the text does not say what the window is when no `a` violates the inequality (the code takes `{a_max}`).
   For `n ≥ 28` a violator always exists (`a = a_max`), so this never matters for `n ≥ 496`.
5. R10 l. 87 ("Novelty: … no hits; no further search yet") is superseded by §4 below.
6. The theorem's by-product `ρ(n) ≤ (4/√(2πe))e^{4.128/√n}` (`n ≥ 12`) is sharp in the constant: `ρ(3000)=0.98061`,
   `ρ(16100)=0.97320`, limit `0.96788`. The proof works only because `4/√(2πe) < 1`; there is no hidden slack to worry about.

---

## 3. Independent recomputation (own programs, in `audit4/`)

| File | What | Result |
|---|---|---|
| `a1_formula_and_source.py` → `a1.log` | exact: forms of `G`,`D`; (1); `lem:monotone` quotient (`n ≤ 400`); brute-force DP check of Prop. 7.2 and of form 1 ⇔ `ρ ≤ 1` (`n ≤ 70`) | all OK |
| `a2_bounds.py` → `a2.log` | **3(a)**: `arb` balls, every inequality of Steps 1–2 and the final bound vs exact `ρ(n)`, all `3 ≤ n ≤ 3000`, all `a` | 0 failures; bound dominates for every `n` (also for `3 ≤ n ≤ 11`); `max ρ/bound = 0.93959` at `n=3000`; e.g. `ρ(496)=0.996147 < 1.164985`, `ρ(2000)=0.982773 < 1.061476` |
| `a3_all_a.py` → `a3_full.log`, `a3_all_a_3_16100.tsv` | **3(b), no monotonicity used**: exact integers (FLINT `fmpz`), `(n−2a+1)·a·2(2n−1)·C(n,a)² vs n(n−1)C(2n,n)` for **every** `1 ≤ a < n/2`, every `3 ≤ n ≤ 16100` — 64 794 450 pairs, 95 s | `R_n(a) < 1` for all `a`, all `496 ≤ n ≤ 16100`; `R_n(a)=1` never; `ρ(n)>1` for exactly 454 values: `n=3`, `5..375`, and 82 values in `376..495`; largest `495` |
| `a4_step4.gp` → `a4_step4.log` | **3(b), second implementation**: PARI/GP, direct `binomial()`, window `[a_min−1, a_max]` with the criterion evaluated at every `a` (479 554 window pairs, largest window 45) | same conclusions; additionally `{a : R_n(a−1) ≥ R_n(a)}` is always an interval ending at `a_max` (unimodal) for `3 ≤ n ≤ 16100`; `ρ(494..498) = 0.9962124, 1.0000241, 0.9961471, 0.9999553, 0.9963327`; closest call above 495 is `ρ(497) = 0.99995528` |
| `a5_symbolic.py` → `a5.log` | sympy/arb checks of all algebra in §2 and the constants of Step 3 | all OK |
| `a6_bruteforce_form1.py` → `a6.log` | full-grid brute force of the verbal conjecture, `n ∈ {9,20,375,376,377,423,495,496,497,498}` | see §1 |

**3(c).** `n ∈ 376..495` with `ρ(n) ≤ 1` (two independent programs, identical, 38 values, all with `ρ<1`):
`376, 378, 423, 425, …, 463 (odd), 466, 468, …, 494 (even)` — exactly the source's list and R10's. `ρ(495) = 1.0000240709 > 1`
(argmax `a = 240`, the only `a` with `R_495(a) > 1`); `ρ(496) = 0.9961471034` (argmax `a = 240`).
These agree with `cert/traffic496_cert.txt`.

---

## 4. Novelty

**The conjecture has already been (reportedly) proved, publicly, in Lean 4.**

* GitHub issue https://github.com/the-omega-institute/trureturing/issues/8249 — "research: antidiagonal grid-traffic
  threshold 496 (Conjecture 7.4)", opened 16 Sep 2026, closed as *completed* 20 Sep 2026. Target: Conjecture 7.4 of
  arXiv:2609.01562v1, in the strict arithmetic form
  `∀ n a, 496 ≤ n → 1 ≤ a → 2a < n → n·(n−2a+1)/(n−a)·C(n,a)·C(n−2,a−1)/C(2n−2,n−1) < 1` — this is literally
  `G(n,a)/D(n) < 1`, the statement of R10. The grid consequence is taken from the paper as "literature input"
  (same split as R10, same domain accounting as §1: Section 6, Props. 7.1/7.2, Thm 5.1, Lemma 2.2).
* PR https://github.com/the-omega-institute/trureturing/pull/8934 — "math: prove the uniform antidiagonal grid traffic
  inequality", merged 20 Sep 2026.
* Lean file (read at commit `4143c05`):
  https://github.com/the-omega-institute/trureturing/blob/4143c05313b9c81c08e0f5e895eef889e14b97f5/D5/S3/Arith/FactorialRatio/GridAntidiagonalTrafficBound.lean
  — 652 lines, `def claim`, `theorem result : claim`, a private `uniform_peak_bound` (Stirling/entropy envelope at
  same-parity peaks); no `sorry`, `axiom` or `native_decide` in the text. Problem note:
  https://github.com/the-omega-institute/trureturing/blob/4143c05313b9c81c08e0f5e895eef889e14b97f5/Problems/gil-liang-odetola-weiner-antidiagonal-traffic.md
  **I did not compile the Lean file**, and the issue page itself was not reachable from this session (content taken from
  the GitHub search API and raw.githubusercontent.com). The route differs from R10's (peak location + Stirling envelope
  + local certificates vs. R10's uniform Gaussian bound + finite exact check), so R10 is at most an independent,
  shorter human-readable proof.
* PR https://github.com/the-omega-institute/trureturing/pull/8676 (18 Sep 2026) records an earlier finite check of the
  same conjecture "`n ≤ 10^5`, no counterexample".
* Why R10's search missed it: `"2609.01562" is:issue` returns 0 because the issue writes `2609.01562v1`; the keyword
  query `"maximal traffic" obstruction is:issue` finds it.

Other searches (all negative):
* arXiv API `id_list=2609.01562`: only v1; arXiv full-text query `"maximal traffic" AND obstruction`: only the paper itself.
* Semantic Scholar citations endpoint: HTTP 200, `{"data": []}` (no citing papers); the paper-record endpoint stayed at 429.
* OpenAlex title search: one record (the paper itself).
* Zenodo (`"2609.01562"`, `"maximal traffic" AND obstruction AND lattice`, title search): 0 records.
* GitHub repositories `"maximal traffic" lattice paths obstruction`: 0.
* Two web searches (title + "conjecture 496 antidiagonal proof"; authors + "proved OR proof OR Lean"): only
  https://arxiv.org/pdf/2609.01562 itself; no follow-up paper, no refutation.

Nobody has refuted the conjecture; one prior public (formal) proof claim exists.

---

## 5. Required corrections to R10

1. Status/novelty line: cite the prior Lean proof (trureturing #8249 / PR #8934, 20 Sep 2026); do not present the result
   as new. Any priority or "first proof" wording must go.
2. Optional wording fixes listed in §2.7 (items 1–4).
