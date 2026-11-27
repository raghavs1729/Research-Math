# Referee report on `CHEBYSHEV1853.md`

Date: 2026-11-24. Scope: §§0–10 of `CHEBYSHEV1853.md` and the scripts and outputs in this directory. I rely on
`../explicit/EXPLICIT.md` ([E]) Theorem 1 (`Θ_1(t) > 0` for `t > p_0` under GRH(χ)) and its Lemma 3(a). I re-derived
everything else from first principles or recomputed it with my own code. My scripts are in `referee_checks/`:
`nc.py` (sieve, C, W, N_C), `ef.py` (explicit formula against direct sums), `beta.py` (an independent
implementation of κ^♭ and β(a)), and `fujii.py` (C_α counterexamples at 40 digits). Each runs in under 3 minutes
on one core at nice 10.

**Verdict.** I found no FATAL error. The mathematics of Theorems 1–6 and Corollary 4a holds as written. Every
derivation I redid came out the same, and every number I recomputed agrees with the document to the printed
digits. There are four REPAIRABLE items. Each is a gap in labelling or in a certificate, not in the mathematics.
I have closed the numerical gaps myself, below. The rest are MINOR.

---

## 1. Checks performed (what was verified, and how)

### 1.1 Kernel identity, convexity, convex weights (Lemma 1, Theorem 1)

* **Partial summation.** Partial summation over `[p_0^-, T]` gives `Σ_{p≤T} s(p) log p F(p) = F(T)Θ(T) − F'(T)Θ_1(T) + ∫_{p_0}^T Θ_1 F''`.
  * There are no lower boundary terms: `Θ(p_0^-) = 0` and `Θ_1(p_0) = 0`.
  * The prime 2 (q = 4) and the prime 3 (q = 3) have `s = 0`.
  * At infinity the terms vanish when `tF → 0` and `t²F' → 0`, because unconditionally `|Θ| ≤ c_θ t` and `|Θ_1| ≤ c_θ t²/2`. Correct.
* **Convexity of the kernel.** With `F = e^{−t/x}/log t` I get symbolically
  `F'' = e^{−t/x}[1/(x² log t) + 2/(x t log² t) + 1/(t² log² t) + 2/(t² log³ t)]`.
  Every term is positive for every `t > 1` and every `x > 0`. This includes `t → 1+` (the terms blow up, positively)
  and `t ≫ x`. The range of integration is `[3,∞)` (resp. `[2,∞)`), so it is covered.
* **General weights.**
  * `F = w/log` and `F'' = w''/log + 2w'(1/log)' + w(1/log)''`, with `(1/log)' < 0` and
    `(1/log)'' = (log t + 2)/(t² log³ t) > 0` on `t > 1`. So "w ≥ 0, nonincreasing, convex" does imply that
    `w/log` is convex on `(1,∞)`, which is what the argument needs.
  * The decay conditions in the statement (`C²`, `Σ_p w(p) < ∞`, `t w(t) → 0`) are sufficient. `T²F'(T) → 0` follows from
    `|w'(T)| ≤ 2(w(T/2) − w(T))/T`, which is convexity, so `T²|w'(T)| ≤ 4 (T/2) w(T/2) → 0`.
  * `Σ w(p) < ∞` is not implied by the other conditions. Example: `w = 1/(t log log t)`.
  * Strict positivity holds because a nonincreasing w that is not identically 0 on `(3,∞)` is positive on some `(3, 3+δ)`.
* **Fujii's weights.** `w = e^{−(t/x)^α}` gives `w'' = α t^{α−2} x^{−α} e^{−(t/x)^α}[α (t/x)^α − (α − 1)] ≥ 0` for `α ≤ 1`, every
  `x`, every `t > 0`. The decay conditions hold, and so does `p^{−σ}`, `σ > 1`.
* **Verdict on Lemma 1 and Theorem 1:** correct.

### 1.2 Explicit formula with Γ(ρ) (Lemma 3)

* **Pole structure.** I derived it from `ψ_e = (2πi)^{−1}∫ (−L'/L)(s) Γ(s) y^s ds` for odd primitive χ:
  * each zero ρ gives `−Γ(ρ) y^ρ`;
  * `s = 0` gives `−b`, with `b = (L'/L)(0)`; here `L(0,χ) ≠ 0` and Γ has a simple pole;
  * `s = −2n` gives simple Γ-poles;
  * `s = −(2n+1)` gives double poles, where a trivial zero meets a Γ-pole.
* **The double pole at `s = −1`.** Its residue is `(log y + 1 − γ + c)/y`, where c is the constant term of `L'/L` at −1.
  This is exactly the `1/y` part of `[(1+b_0)e^{−1/y} + E_1(1/y)]/y` with `b_0 = c`. I computed c independently in mpmath:
  `b_0 = −0.0602669141769` (χ_{−4}) and `0.1438316445629` (χ_{−3}), in agreement with [E].
* **Term-by-term integration.** I redid the integration of [E] Lemma 3(a):
  * `Γ(ρ+2)/(ρ(ρ+1)) = Γ(ρ)`;
  * `y^{−2}∫_1^∞ t e^{−t/y} = (1+1/y)e^{−1/y}`;
  * `y^{−2}∫_1^∞ log t e^{−t/y} = E_1(1/y)/y`;
  * `Σ 1/((2n−1)2n(2n+1)) = log 2 − 1/2`;
  * `|R_1| ≤ (2/5)σ_0²/y²`.

  All correct.
* **Numerical test.** I compared the formula with the direct sum `Σ χ(n)Λ(n)e^{−n/y}` at y = 1, 2, 5, 20, 100, 1000, using all 511 (537) zeros
  of `prototype/explore/chi{4,3}_zeros.txt` and the closed form
  `Σ_n u^{2n}/(2n(2n+1)) = 1 − ½ log(1−u²) − atanh(u)/u` for R_4.
  * χ_{−4} differences: `−3.3e−9, −1.55e−9, −3.6e−10, −2.7e−11, −1.1e−12, ~1e−13`. This is the document's table, and the residual is the R_1 truncation.
  * A first run, with a series for R_4 truncated at 400 terms, gave residuals near `1.5e−7/y²`. The cause was my R_4, not the formula.
  * χ_{−3}: same pattern.
  * A wrong b, b_0 or sign would show as O(1) or O(1/y). It does not.
* **The two inequalities used for `y ≥ 4`.** I checked them term by term, including `(1+u)e^{−u} ≥ 1 − u²/2`, A&S 5.1.20, and
  `b − r^+(4) = 0.11753 / 0.22742`.

### 1.3 S_Γ (Lemma 4)

* `|Γ(1/2+iγ)|² = π/cosh πγ`. Recomputed from the list: `2Σ_j |Γ(1/2+iγ_j)| = 3.91993271294e−4` (χ_{−4}) and
  `1.65313139264e−5` (χ_{−3}), identical to `cert_abel.out`.
* **Closed form.** I re-derived `Σ_ρ 1/(ρ(1−ρ)) = 2 Σ Re(1/ρ) = log(q/π) − γ + 2(L'/L)(1,χ)` from the Hadamard product at
  s = 1 for odd χ (`ψ(1) = −γ`). It is unconditional, and it equals `Σ|ρ|^{−2}` under GRH.
  Value for χ_{−4}: `0.155567979924`.
* **The completeness argument is sound.**
  * Disjoint certified sign-change intervals give distinct zeros. So the unlisted zeros, and extra multiplicities, carry
    `Σ'|ρ|^{−2} = d = 0.0034836`.
  * d is almost exactly the expected tail of the genuine zeros above the end of the list (height ≈ 640):
    `(1/π)(log(qT/2π)+1)/T ≈ 3.5e−3`.
  * Any missing pair therefore has `|γ| ≥ √(2/d − 1/4) = 23.955`. A missing zero at γ = 0 would contribute 4 > d, so
    `L(1/2) ≠ 0` is not even needed.
* **The tail bound is rigorous and covers all zeros beyond the list.** `|Γ(ρ)| = |ρ|^{−2}·h(γ)` with h decreasing, so
  `Σ'|Γ(ρ)| ≤ d·h(g) = 2.3e−16`. The arithmetic is right.
* Overstatement (MINOR M4): "all the zeros up to height 24" should read 23.955.

### 1.4 Volterra mixture (Lemma 5)

* **Kernel.** The kernel is `μ(λ) = ∫_0^∞ λ^{u−1}/Γ(u) du`. It is positive because `1/Γ > 0` on `u > 0`, and finite because
  `1/Γ(u) ~ u` at 0 and decays faster than exponentially.
* **Laplace identity.** `∫ e^{−λt} μ = ∫_0^∞ t^{−u} du = 1/log t` for `t > 1`. Correct.
* **Mixture identity.** `−C(x) = ∫ W((λ+1/x)^{−1}) μ dλ` is correct. The interchange is justified by **Fubini**: absolute
  convergence of `Σ_p log p e^{−p/x}∫e^{−λp}μ = Σ_p e^{−p/x}`. "Tonelli" is the wrong word because s(p) changes sign (MINOR).
* **General α.** The version with `(λ + x^{−α})^{−1/α}` is correct.
* **Lemma 5(ii).** `∫_0^1 u λ^{u−1} du = (1 − λ + λ log λ)/(λ log² λ)` is correct.
* **Expansion.** `μ λ log²λ = 1 + 2γ/a + 6(γ²/2 − π²/12)/a² + …` with `a = log(1/λ)` gives 1.0511 at `λ = 10^{−8}`, as printed.
* **Theorem 4.** I re-derived the substitution `λ = v/x`, the `c_a` integral (`= 2 log((1+√2)/2)`), `A(V)`, the m-loss term,
  and the asymptotic constant `2 log 2`. I also re-derived the true constant `γ + 2 log 2` from
  `∫_0^∞ e^{−s²} log s ds = −(√π/4)(γ + 2 log 2)`. All correct.

### 1.5 Corollary 4a: values, constant, junction

My own sieve, in double precision with fsum (`nc.py`), gives `N_C(x) = −C(x) log x/√x` for χ_{−4} as follows:

| x | 68.5 | 69 | 69.79 | 70 | 71 | 100 | 10³ | 10⁴ | 10⁵ | 10⁶ |
|---|---|---|---|---|---|---|---|---|---|---|
| N_C | 0.879172 | 0.880335 | 0.882150 | 0.882628 | 0.884877 | 0.934571 | 1.059554 | 1.106609 | 1.090514 | 1.062529 |

* The crossing of 0.88 is at `x = 68.8554`.
* The peak is `1.1087617`, at x ≈ 6360. The document's grid maximum is 1.108761 at 6309.6. Same.
* χ_{−3}: `N_C(13.43) = 0.88005` and `N_C(14) = 0.88777`. Same.
* Every C value equals the document's table to all printed digits.

**The junction.** The Arb certificates cover `[69.783, 2001.01]` and `[1995.26, 5.01·10^6]`. The analytic range starts at
`X_A = 10^{6.6} = 3.98·10^6`.

* I re-implemented κ^♭, κ^⋆, `K_E`, `K_I` and β(a) independently (`beta.py`). Results:
  * `K_E = 0.8837362`, `K_I = 0.8830301`;
  * min β over `[10^{6.6}, 10^{14}]` is `0.8960482`, at the left end, for χ_{−4};
  * min β is `0.8860773`, at `a = 10^{6.63}`, for χ_{−3};
  * `A(64) − A(1) − c_a = 1.13694`;
  * `Ē = ε_a F(2^{22}/10^6) = 8.489e−11 ≤ 8.5e−11`.

  These agree with the document.
* **Where the margins come from.** 0.016 and 0.006 are `min β − 0.88`, and 0.003 is `K_I − 0.88` for `x ≥ 10^{14}`.
  These are closed forms. Thirty-digit evaluation has errors near `10^{−28}`, so floating point cannot flip any of these
  decisions. What matters is whether the formulas and the monotonicity argument are right, and I checked both by hand:
  * the max in κ switches at `y = 2.1537^8 = 463.6`;
  * every fixed-k term is monotone;
  * newly present odd k contribute at most `c_3 √p_0 p_0^{−k/2}` resp. `c_θ p_0 p_0^{−k/2}`, summing geometrically to the J term;
  * at most 26 odd k occur below `Y_I`, each with `c_θ τ_k ≤ 10^{−3200}`, which is covered by `10^{−3000}`;
  * the regime constants `K_E` and `K_I` are lower bounds on `(Y_E, Y_I]` and `(Y_I, ∞)`.
* **The grid argument is sound.**
  * `log x / log²(x/V)` is decreasing in x, since its derivative has the sign of `−(L−c)(L+c)`.
  * The code's use of `log a` in the numerator is conservative.
  * The m-loss term is monotone as claimed.
* **The Arb Taylor certificates (`certify_grid.py arb`).**
  * The Stirling-polynomial recursion is correct: `d/dh[P_k(w)e^{−w}] = w(P_k − P_k')e^{−w}`.
  * The remainder bound is correct: `w ∈ [v, v e^{2δ}]` and `Σ_n f(n/X) ≤ X∫f + max f`.
  * The tail beyond `60x` is correct (`f` is decreasing past m ≤ 8).
  * Cells tile.
  * `log x/√x` is evaluated at the right end on cells with left end `> e²`.
  * The certified lower bound `0.88079` on `[69.78, 70.19]` is consistent with my values.
* The full chain `[70, ∞)` is covered.

### 1.6 Theorem 5

The infinite range is handled rigorously.

* `y²W(y) = ∫_{p_0}^∞ Θ_1 e^{−t/y}`. On `[p_0, T]` this is bounded below by the part over `[5,7]` (resp. `[3,5]`), with
  `A = 2(4 log 3 − 2 log 5) = 2.35115` (resp. `2 log 2`).
* The tail `t > T` is bounded by `|Θ_1(t)| ≤ c_θ t²/2`, which needs `ϑ_o(u) ≤ c_θ u` for all u. That holds by the
  certificate up to 2^22 and by BKLNW above `3 594 641`. It gives `(c_θ/2) y e^{−T/y}(T² + 2Ty + 2y²)`.
* The ratio is increasing in y and equals `2.37e−9` at `y = T/60`. I recomputed it by hand.
* Then `C < 0` follows through Lemma 5(iii), because all arguments `(λ+1/x)^{−1} < x`.
* **BKLNW.** I extracted Corollary 11.2 of arXiv:2002.11068v2 from the PDF. The rows `(k, M_k, X_0) = (2, 0.2, 3 594 641)` and
  `(3, 1, 89 967 803)` are there exactly as cited. There is one gap here, which is REP-1.

### 1.7 Landau converse (Theorem 2)

* Step 1, the Mellin transforms, is correct. Near `x = 0`, `|C(x)| = O(e^{−3/x})`, so the `[0,1]` piece is entire.
* The Euler-product identities for P and Z are correct, including the regrouping into H and H_1, the identity
  `Σ_{p∤q} log p p^{−2s} = −ζ'/ζ(2s) − Σ log p/(p^{2s}(p^{2s}−1)) − log ℓ ℓ^{−2s}`, and analyticity in `Re s > 1/2`.
* The real axis is handled correctly:
  * `L(σ,χ_{−4})` and `L(σ,χ_{−3})` are positive for `σ > 0`, as alternating series with decreasing terms;
  * `ζ(2σ) > 0` for `σ > 1/2`;
  * the pole of `ζ(2s)` is at the boundary point 1/2.
* Landau's theorem applies to the one-signed tail. Continuation through a thin neighbourhood of the real ray is legitimate,
  because the overlap with `{Re s > σ_c}` is connected.
* `exp(continued log L) = L`, so there are no zeros in `Re s > 1/2`. The functional equation gives the same for `Re s < 1/2`.
* The W version, through `L'/L`, is correct.
* The cycle (a)⇒(b),(d); (b)⇒(b'); (c)⇒(b'); (d)⇒(d'); (b'),(d')⇒(a); (a)⇒(c) is complete.
* The unconditional remark ("`C < 0` for arbitrarily large x") is correct.

### 1.8 Fujii exponents

* **Definition.** I read Platt–Trudgian (arXiv:1802.00089v1, full text extracted from the PDF). Their sum is
  `Σ_{p>2}(−1)^{(p−1)/2} e^{−(xp)^α}` with `x → 0`. Under `x ↦ 1/x` it is the document's `C_α`.
* **The counterexamples.** I recomputed them independently at 40 digits, summing over all primes `p ≤ 1.2x`:
  * `C_2048(623 658) = +1.6335030788763`;
  * `C_4096(26 862.745) = +0.0088978249978`;
  * `π(x;4,1) − π(x;4,3)` first becomes positive at 26 861.

  Both counterexamples are confirmed, and they are unconditional.
* **The criterion.** From the 511 listed zeros: `Σ_ρ|Γ(ρ/α)| = 20.12766214` against `Γ(1/(2α))/2 = 20.12766492` at
  α = 20.40442. The root is `α_* = 20.40442354`. This reproduces Platt–Trudgian's numbers (`< 20.1276643`, `> 20.1276649`)
  and the document's table.
* **Normalisation.** The document's version, `(1/α)Σ|Γ(ρ/α)|` against `Γ(1+1/(2α))`, is the same inequality, because
  `αΓ(1+1/(2α)) = Γ(1/(2α))/2`.
* **What is and is not proved for `1 < α < 20.4`.** The document claims proofs only for `α ≤ 1` (F1) and for the identity
  (F2). Everything else is labelled exploration. That is accurate.

### 1.9 Attribution and novelty

* **Platt–Trudgian** are described accurately:
  * Fujii's `0 < α < 4.19` using only `γ_1 > 6` and `Σγ^{−2} < 1/5`;
  * their Theorem 1 for `α < 20.40442`;
  * 1000 lcalc zeros checked by Arb and by Turing's method up to height 1127;
  * the S_1…S_4 decomposition with `o(·)` terms;
  * the remark that 20.40443 is out of reach of the method.

  The document also correctly notes that their `S_1 = −½Γ(1/(2α))x^{−1/2} + o(x^{−1/2})`, as printed for the unweighted sum,
  is missing the logarithm. This does not affect their criterion. There is no all-x statement and no explicit o-term there.
* **Suzuki (arXiv:2411.07436v3).** I checked it via the HTML text. It contains sign-for-`x ≥ x_0` ⇔ GRH(χ_4) equivalences for
  other weights (Theorem 3(2), the weight `log(x/n)/√n`; the Knapowski–Turán weight is cited), but no statement for all
  `x > 0`. The document's claim that the all-x form is "apparently new" is plausible and properly hedged. Suzuki's
  Theorem 3(2) and Knapowski–Turán should be cited as close relatives of Theorem 2(b') (MINOR M6).
* **Corollary 7** is properly conditional ("*If* that proof is correct …; GRH(χ) remains a hypothesis of every theorem").
  Small issues are under M7.

---

## 2. Findings

### FATAL

None.

### REPAIRABLE

**REP-1 (certificates stop at the last prime, not at `T = 2^22`).** Several statements are asserted on `(0, 2^22]`:

* Lemma 7: "`Θ_1 > 0` on `(p_0, T]` (piecewise linear)";
* (P1): `u − c_1√u ≤ ϑ_o(u)` on `(0, U_1]`;
* (P3): `|Θ(u)| < c_3√u` on `(0, U_1]`.

`cert_abel.py` A1/A2 test only at primes `p ≤ T`, so the final segment `[4 194 301, 4 194 304]` is not covered. For A1,
`Θ_1` is linear on that segment with slope `Θ(4194301)`, and positivity at the left end does not imply positivity at T.
A3 is fine, because `|Θ|/√u` decreases between primes. A2 is not: `(u − ϑ_o(u))/√u` increases on the segment.

*Fix:* add the endpoint `u = T` to A1 and A2. I evaluated it: `Θ_1(T) = 5.89·10^9` (χ_{−4}) and `6.05·10^9` (χ_{−3}), and
`(T − ϑ_o(T))/√T = 1.0144 / 1.0146 < c_1`. Both are double precision with margins of nine orders, resp. a factor 2, so
the statements are true. The certificate should still include the endpoint.

**REP-2 (Taylor certificates labelled "certified" use floating-point pieces).** §0 defines "certified" as "no
floating-point decision". In `certify_grid.py arb`, three pieces are evaluated in IEEE double and then multiplied by
`1.000001`:

* the Taylor remainder R;
* the prime tails `xΓ(m+1,60) + 60^m e^{−60}`;
* the tail argument `xf`.

Section 5 discloses this. The safety factor exceeds double rounding by about 10^9, and R is below `10^{−7}` of `−C`, so no
decision can flip.

*Fix:* evaluate R and the tails in Arb (a few lines), or relabel the two certificates "Arb with
double-precision remainder (inflated by 10^{−6})". Corollary 4a is unaffected.

**REP-3 (LI statements presented as facts).** Two places state as fact an argument that §9.7 itself says is only a sketch:

* the status paragraph ("the constant cannot exceed `√π/2 − Σ|Γ(ρ)| = 0.88583…` if the ordinates are linearly independent");
* the parenthesis in Corollary 4a ("Under GRH and LI one has `liminf N_C = √π/2 − S_Γ^{true}`, so 0.88 cannot be replaced by 0.8859").

The claim is very likely true. The log-singularity of `log L` at ρ gives a zero term `Σ Γ(ρ)x^{iγ}(1 + O(1/log x))` in `N_C`.
Kronecker under LI, which needs simple zeros and linear independence, then makes the liminf `√π/2 − S_Γ`. But this is not proved in the note.

*Fix:* write "expected (sketch, §9.7)" in both places, or prove it. The required proof is a standard
Mellin computation with `log L` and uniform control of the zero sum.

**REP-4 (statement of Theorem 1's weight class and of (F1)).** Three small corrections:

* Theorem 1 is stated for `w` on `[3,∞)` (χ_{−4}). For the χ_{−3} version claimed in Theorem 6, the domain is `[2,∞)`.
  Say so explicitly.
* The Remark after Theorem 1 says that "*every* nonnegative, nonincreasing, convex weight inherits the sign". It omits the
  decay conditions (`Σ_p w(p) < ∞`, `t w(t) → 0`), which are needed. Without them the sum need not converge, as with
  `w = 1/(t log log t)`.
* (F1) says that `e^{−(t/x)^α}` "is convex exactly when α ≤ 1". As a statement on `[3,∞)` for fixed x this is false: for
  α = 2 the weight is convex on `[3,∞)` whenever `x ≤ 3√2`. The correct statement is "convex on `[3,∞)` for every
  `x > 0` exactly when α ≤ 1". What the method needs is that `w/log` be convex on `[3,∞)`. For α > 1 and large x this fails
  near `t ≍ x`, because the term `−α(α−1)/(e x² log x)` beats the positive terms of order `1/(x² log² x)`. So the
  conclusion that the method stops at α = 1 stands.

### MINOR

* **M1.** Theorem 5's statement contains floating-point claims in parentheses: the rounding-budget scan
  `[0.103, 1.02·10^8]` and the grid scan. Move them out of the theorem.
* **M2.** Lemma 5(iii) says "Tonelli"; it should be "Fubini (absolute convergence)".
* **M3.** In "`E(y)/√y ≤ ε_a √π/2 = 7.63·10^{−4}`", the value is `7.622·10^{−4}`. The inequality as printed is still true.
* **M4.** In Lemma 4, "complete up to height 24" should read "up to height 23.955". One could also cite Platt–Trudgian's
  Turing-method verification of all zeros of `L(s,χ_4)` up to height 1127.
* **M5.** Corollary 7 cites "Theorems 1 and 4" for the bound `0.88√x/log x`; it should cite Corollary 4a.
* **M6.** For Theorem 2 and §7, cite Suzuki's Theorem 3(2) and Knapowski–Turán as nearby sign-for-large-x ⇔ GRH results.
* **M7.** Corollary 7 rests on `../GRH-PROOF-DIGEST.md`. The digest also records two things the corollary does not mention:
  * the default-kernel stage of `lake comparator` was still running;
  * an unread PhilArchive paper alleges vacuity.

  The conditional wording is adequate, but one sentence noting this would be more careful.
* **M8.** In §1, the "maximum 1.10876 at x ≈ 6310" is the grid maximum. The true maximum is `1.1087617` near x ≈ 6360,
  which makes no difference.
* **M9.** Theorem 3: the threshold "κ > 0 for y ≥ 413 / 750" is a grid statement. It is labelled as such, and is not used.

---

## 3. Certification list

**Proved as written** (given [E] Theorem 1 and [E] Lemma 3(a)):

* Lemma 1 and both kernel instances, including `F'' > 0` for all `t > 1`, `x > 0`.
* Theorem 1 (GRH(χ_{−4})): `C(x) < 0` and `W(x) > 0` for every `x > 0`, and the convex-weight statement under the hypotheses
  as written (C², nonincreasing, convex, ≥ 0, not ≡ 0, `Σ_p w(p) < ∞`, `t w(t) → 0`). This includes Fujii's weights for
  `0 < α ≤ 1`, all x, and `p^{−σ}`, σ > 1.
* Theorem 2, all equivalences.
* Lemma 3, Lemma 4 (`S_Γ ≤ 3.9199328·10^{−4} / 1.6531314·10^{−5}` under GRH), Lemma 5, Theorem 3, and Theorem 4.
* Corollary 4a: `−C(x) ≥ 0.88 √x/log x` for every `x ≥ 70` under GRH(χ_{−4}), and `liminf N_C ≥ √π/2 − S_Γ`. This is subject
  to REP-2, which is labelling only.
* Theorem 5, unconditional: `C(x) < 0` and `W(x) > 0` for `0 < x ≤ 69 905`, both characters. This is subject to REP-1;
  I verified the missing endpoint, and the statement is true.
* Theorem 6, the χ_{−3} analogues, including `−C_3(x) ≥ 0.88√x/log x` for `x ≥ 14`.
* §6 (F1), (F2), and the two unconditional counterexamples `C_2048(623 658) > 0` and `C_4096(26 862.745) > 0`.
* Corollary 7, as a conditional statement.

**Not certified as written, with corrected statements:**

* The status-paragraph and Corollary 4a claims "under GRH and LI, `liminf N_C = √π/2 − S_Γ`; the constant cannot exceed
  0.88583…; 0.88 cannot be replaced by 0.8859". *Corrected:* "Under GRH and LI one expects (sketch, not proved here)
  `liminf N_C = √π/2 − S_Γ`".
* (F1) "the weight is convex exactly when α ≤ 1". *Corrected:* "`e^{−(t/x)^α}` is convex on `[3,∞)` for every `x > 0`
  if and only if `α ≤ 1`".
* Remark after Theorem 1: "every nonnegative, nonincreasing, convex weight". *Corrected:* add "with `Σ_p w(p) < ∞` and
  `t w(t) → 0`".
* Lemma 4 remark "all zeros up to height 24". *Corrected:* "up to height 23.955".
* Theorem 5's parenthetical claims for `[0.103, 1.02·10^8]` and `N_C ≥ 0.88` on `[74.2, 1.02·10^8]`. These are floating
  point with an assumed libm error bound, as the document says, and are not part of the theorem.
* Everything in §6 beyond (F1), (F2) and the two counterexamples (F3–F5, `α_W ≈ 6.76`, monotonicity of S/β). These are
  exploration, correctly labelled.
