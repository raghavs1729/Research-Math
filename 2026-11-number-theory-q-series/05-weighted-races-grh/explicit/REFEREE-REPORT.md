# Referee report on `monotonicity/explicit/EXPLICIT.md`

Date: 2026-11-24. Scope: every claim in §§1–7 of EXPLICIT.md (§8 was read for
Lemma 6 and the cone remark only). All checks below were done from scratch with my own
code, in `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/`
(`ef.py`, `scan.py`, `cert.py`, `Lb.py`, `ident.py`, `nsmall.py`). The author's
scripts were used only to read the raw output files. One core, each run under 2 minutes.

**Verdict in one line.** I found no FATAL error. Every identity, sign and constant I
re-derived agrees with the note. There are four REPAIRABLE items: two concern sources and
rigour of inputs, one is a proof sketch to write out, and one is about novelty claims.
There are also some MINOR points.

---------------------------------------------------------------------------

## 1. Identities (re-derived)

**Lemma 2 / (2.1).** Correct. Θ vanishes on [0, p_0) because s(ℓ) = 0, and Θ_1(p_0) = 0.
The Stieltjes integral starts at p_0^- and the jump at p_0 is included. The second
integration by parts has no boundary term at p_0. For h(t) = (x−t)/log t, I get
h'' = 2/(t log²t) + (x−t)(1/(t² log²t) + 2/(t² log³t)), which matches the note. It is
strictly positive on all of (1, x], so it is positive on [p_0, x] for both p_0 = 3 (χ_4)
and p_0 = 2 (χ_3). p = 2 (χ_4) and p = 3 (χ_3) have weight 0 and play no role. Starting
the integral below p_0 would change nothing, because Θ_1 ≡ 0 there. h(x) = 0 and
h'(x) = −1/log x, so the boundary term is Θ_1(x)/log x with a + sign. I checked the
identity numerically to 30 digits, with piecewise quadrature, at x = 200.5.

**Theorem 5 (cone).** Correct for all real m ≤ 0 and m < m' ≤ 1, including m' ≤ 0.
f = t^m(1−(t/x)^κ) ≥ 0. f' = t^{m−1}[m − m'(t/x)^κ] ≤ 0 in both cases m' ≥ 0 and
m' < 0; the second uses |m'| < |m|. f'' = t^{m−2}[m(m−1) − m'(m'−1)(t/x)^κ] ≥ 0. Here
m(m−1) ≥ 0. When m' ∈ [0,1], m'(m'−1) ≤ 0. When m' < 0, m'(m'−1) < m(m−1), because
u(u−1) is decreasing on u < 1/2. 1/log t is positive, decreasing and convex, and
(fg)'' = f''g + 2f'g' + fg'' ≥ 0. F(x) = 0 and F'(x) = −κx^{m−1}/log x. I verified the
identity and the lower bound to ~30 digits at (x,m,m') = (200.5,0,1), (100.3,−5,−4.9),
(57.7,−30,1), (40.2,−1,−0.999), (90.1,0,0.01) and (33.3,−0.5,1). These cover m very
negative, m' very close to m, and m' < 0.

**Lemma 6 (§8).** I re-derived it. With u = (t^{−κ} − x^{−κ})/log t, one has
−u' = κt^{−κ−1}g, dΦ = κt^{−κ−1}A dt and Φ(p_0) = 0, so the stated w_x is correct and
positive for κ > 0. The remark "Φ_{m,m'} ≥ κ t^{m−1}Θ_1" is also correct. It is Lemma 2
with F = f_t, without the 1/log factor.

**Theorem 6, the Mellin identity.** ∫_1^∞ A x^{−s−2} dx = P(s)/(s(s+1)) is correct.
P = −log L(s,χ) + ½ log L(2s,χ_0) + H, with H analytic in Re s > 1/3, is also correct.

**Corollary 4.** D_0 = D_1/x + R. Both implications follow from R > 0, and the stated
reading of the result is correct.

**§4 consequence.** For q = 4, 2π(t;4,1) = π(t) − 1 − D_0(t) for t ≥ 2, so
∫_2^x(2π(t;4,1) − li) = ∫_2^x(π − li) − (x−2) − xR(x). That is correct. For q = 3 the
extra term is −(x−3)_+, and the conclusion is still correct.

## 2. Explicit formula for ψ_1(x,χ) (Lemma 3(a))

I re-derived the residues of (−L'/L)(s) x^{s+1}/(s(s+1)):

- at ρ: −x^{ρ+1}/(ρ(ρ+1));
- at 0: −bx, since L(0,χ) ≠ 0 for odd χ;
- at −1 (double pole): log x + 1 + b_0;
- at −(2n+1), n ≥ 1: −x^{−2n}/(2n(2n+1)).

The formula agrees with the note. At x = 1 the constant is consistent: the Hadamard
product gives Σ1/(ρ(ρ+1)) = b_0 − b + log 2, and Σ_{n≥1} 1/(2n(2n+1)) = 1 − log 2. I
also re-derived b and b_0 from the functional equation, and they agree with the closed
forms in the note.

Numerical test (`ef.py`). I computed b and b_0 independently with mpmath, used the 511
(χ_4) and 537 (χ_3) zeros, and computed ψ_1 directly. "diff" is ψ_1 direct minus the
formula:

| x | diff χ_4 | diff χ_3 | crude omitted-zero size x^{3/2} log(qT/2π)/(πT) |
|---|---|---|---|
| 1 | 3.48e-3 | 3.09e-3 | 3.0e-3 / 2.6e-3 |
| 3 | 1.6e-3 | −6.5e-7 | 1.6e-2 |
| 5.5 | −1.2e-6 | 1.6e-5 | 3.9e-2 |
| 10 | 7.5e-5 | 1.7e-4 | 9.5e-2 |
| 20 | −1.6e-3 | −8.9e-4 | 0.27 |
| 50 | −3.2e-4 | 3.1e-3 | 1.06 |
| 200 | −8.3e-2 | 3.6e-3 | 8.5 |

At x = 1 the difference is exactly the non-oscillating tail Σ_{γ>T} 1/(ρ(ρ+1)), which is
≈ 0.00349 for χ_4. At every x the error is consistent with truncating the zero sum, and
there is no sign of a wrong constant or a missing term.

A note on the author's check (iii): at x = 10³ and 10⁴ the allowed error (110 and 3487)
is so large that the test is vacuous there. The χ_3 difference of 69 at 10⁴ is not quoted
in §2. The tests that actually say something are those at small x.

## 3. Constants

- **σ_0².** Σ_ρ 1/|ρ|² = log(q/π) − γ + 2(L'/L)(1,χ) under GRH. I re-derived it from
  the Hadamard formula at s = 1, with Re B = −ΣRe 1/ρ and Re 1/(1−ρ) = ½|ρ|^{−2} on the
  critical line. mpmath gives 0.15556797990 (χ_4) and 0.1132299704 (χ_3); these differ
  from the note only in the 11th digit, at the level of numerical differentiation. The
  inequality is used in the right direction: |ρ(ρ+1)| = |ρ||ρ+1| > |ρ|² because
  Re ρ > 0, so S = Σ1/|ρ(ρ+1)| < σ_0². S ≈ 0.15380 / 0.11256 from the zeros plus tail.
- **Step 2 bracket.** bt − log t − 1 − b_0 equals 0.31122 at t = 3 (χ_4) and 0.05942 at
  t = 2 (χ_3). It is increasing for t > 1/b. The dropped sum Σ_n t^{−2n}/(2n(2n+1))
  enters −ψ_1 with a + sign, so dropping it is legitimate.
- **Prime powers.** χ(p^k) = χ(p) for odd k and 1 for even k (p ∤ q). The signs in Step 1
  are correct, and the even-k terms (k ≥ 4) are ≥ 0. B_2 uses the substitution u = y².
  The integrals (2/3, 4/5, 3/4, 5/6) are correct, and so is the count of odd k ≥ 7,
  ≤ max(0, (K−5)/2) with K = log t/log p_0. The split at Y_1 in Step 3 is correct: every
  dropped term has the right sign.
- **c_1.** I recomputed it independently up to 10⁷: 2.0711926180 (q = 4) and
  2.0819411967 (q = 3), both attained as y → 1423^-. sup(ϑ_o − y)/√y < 0. This matches.
- **L(t).** My own code (`Lb.py`), to 30 digits:

  | t | L_{χ_4} | L_{χ_3} |
  |---|---|---|
  | 10⁴ | 0.068154 | 0.018223 |
  | 10⁶ | 0.343225 | 0.358818 |
  | 10⁸ | 0.448221 | 0.483658 |
  | 10¹⁰ | 0.486725 | 0.527391 |

  - L(10²⁰+) = 0.510677 / 0.553014.
  - Positivity thresholds: 4891.69 / 8525.58.
  - L ≥ 1/2 from 5.4768·10¹¹ / 3.398·10⁸.
  - Monotonicity point e^{5 log p_0 + 14/5}: 3996 / 527.

  Every value matches the note and bounds.out. Each one lies below the exact values. For
  χ_4, N is 0.8305, 0.8785, 0.7774, 0.7460 and G is 0.5234, 0.6849, 0.6447, 0.6485 at
  10⁴, 10⁶, 10⁸, 10¹⁰. For χ_3, N is 0.8899, 0.8538, 0.7960, 0.7988.
- **Corollary 3 grid.** I recomputed it independently. β(a) ≥ 1/2 holds on the grid
  from a = 1.5136·10⁸ (χ_4) and 2.512·10⁷ (χ_3) up to the point where L ≥ 1/2. The
  smallest margin from 10¹⁰ on is 0.03463. This matches.
- **Sources** (I read the actual texts this time):
  - Dusart, arXiv:1002.0442, Thm 5.2. It does contain the row k = 3, η_3 = 1,
    x_3 = 89 967 803, i.e. |ϑ(x) − x| < x/ln³x. That is correct as quoted. See item
    REP-1 for its input.
  - Büthe, arXiv:1511.02032, Thm 2. (1.6) says x − ϑ(x) ≤ 1.95√x for
    1423 ≤ x ≤ 10¹⁹; (1.7) says x − ϑ(x) > 0.05√x for 1 ≤ x ≤ 10¹⁹. Correct as quoted.
  - Johnston, arXiv:2201.06184v2. Thm 1.1 and Thm 1.2 are as described, and Lemma 2.6
  has 0.04621. Lemma 2.2 quotes RS (2.2) as `x − x/log x < θ(x) < x + 0.5x/log x` for
  x ≥ 41. The note attributes the upper bound "for x > 1" to this quotation, which is
  wrong; see MINOR-2.
- **No RH for ζ.** Confirmed. The sieve, Büthe and Dusart use only finite
  verifications of RH plus zero-free regions. GRH(χ) is used only in Step 2 and
  Lemma 3(b).

## 4. Junction and finite computation

- **R.** The test p·D_0(p) − D_1(p) > 0 is done in exact int128 arithmetic. R is
  continuous and equals A − B/x on each gap, so it is monotone there. Positivity at the
  primes therefore gives positivity everywhere, and the note's argument is right. The
  last gap [9 999 999 967, 10¹⁰] is not put into the bins, but R(10¹⁰) > 0 covers it.
- **Θ_1.** Θ_1 is continuous and piecewise linear with nodes at the race primes, so
  checking it at the nodes is enough. In the C code, however, the sign of Θ_1 is decided
  in `long double`, not exactly (REP-3). The proof needs only p ≤ 10⁴. I certified
  that range:
  - exact big-integer comparison of Π_{s=+1} r^{p−r} with Π_{s=−1} r^{p−r} for p ≤ 300;
  - 200-bit arithmetic up to 10⁴.

  For χ_4 the minimum is Θ_1(7) = 1.1756, and for χ_3 it is Θ_1(5) = 2.0794 (`cert.py`).
- **Independent scan to 10⁶ and 10⁷** (`scan.py`, numpy):
  - on (3, 10⁶] (q = 4) and (2, 10⁶] (q = 3), R ≤ 0 and Θ_1 ≤ 0 never occur, either at
    primes or at left limits;
  - min_{x>p_0, at primes} R = 2/7 at 7 (q = 4) and 0.6 at 5 (q = 3);
  - min N on [10⁴, 10⁷]: 0.7534644 at the prime 721 291 (q = 4), the same as the
    table's k = 19 entry. The rigorous gap bound is 0.75345, and the true minimum on the
    nearest competing gap is larger;
  - for q = 3: 0.767957 at 6 015 047;
  - checkpoint values (D_0, D_1, R, N, G at 10⁴…10⁷) agree digit for digit with
    out_q4/out_q3.
- **By hand.**
  - R(7) = 2/7: on [7, 11), D_0 = 1 and D_1 = 3 − 5 + 7 = 5, so R = 1 − 5/x. On [5, 7),
    R = 2/x. So 2/7 is the minimum over x ≥ 5.
  - Corollary 3 at 19: on [19, 23), D_0 = 1, D_1 = 5 and R(19) = 14/19. N(19) =
    0.7368 · 0.6755 = 0.4977 < 1/2, so the corollary fails there, as the note states.
    N = 1/2 is crossed at x ≈ 19.46.
  - On [23, 29), R = 2 − 28/x (D_0 = 2, D_1 = 28), which confirms 0.5038 and 0.550.
  - A fine-grid scan gives min N = 0.6134 on [32, 10⁵] (q = 4) and 0.72299 at 11 on
    [7, 10⁵] (q = 3).

## 5. Landau converse (Theorem 6)

The argument is complete once the steps below are written out.

1. A ≥ 0 eventually, A = O(x²), and A is locally bounded. Landau's theorem then puts a
   real singularity at the abscissa σ_c.
2. P(s)/(s(s+1)) continues analytically along the real segment (1/2, ∞). This needs
   L(σ,χ) ≠ 0 for real σ > 1/2, which holds for χ_4 and χ_3 because the alternating
   series is > 0. It also needs ζ(2σ)(1 − ℓ^{−2σ}) ≠ 0 for σ > 1/2.
3. Hence σ_c ≤ 1/2. log L(s,χ) then has a single-valued continuation to Re s > 1/2,
   so L(s,χ) ≠ 0 there.
4. The functional equation for real χ gives the rest of the strip.

**On real zeros.** The non-vanishing of L on the real segment is genuinely needed and
cannot be dropped. A real zero β > 1/2 would give a real singularity of P at β, with
P → +∞ as s → β+. That is *consistent* with A ≥ 0, so Landau alone would only give
σ_c = β and no contradiction. For χ_4 and χ_3 the elementary positivity of L(σ) settles
it, and the note says so. Write this dependence out explicitly.

---------------------------------------------------------------------------

## REPAIRABLE items

**REP-1 (input of Lemma 4(b)).** Dusart arXiv:1002.0442 (dated 2007/2010) derives
Thm 5.2 from Gourdon's 10¹³-zero RH verification. That computation is unpublished and
not formally rigorous; the paper's introduction cites it as [9]. The note says "no RH
for ζ" and all inputs are unconditional. That is true in spirit, but this input is not
of certified quality.

Only Corollary 3 for χ_4 at x > 10²⁰ depends on it; the RS fallback (b') covers
everything else, as §6 says. Corollary 3 needs only |θ(y) − y| ≤ εy for y ≥ 10¹⁰, with
ε ≤ (3/2)(c_0 − 1/2) ≈ 0.0166 and with the jump at 10²⁰ kept above 1/2.

*Fix:* cite a bound built on a rigorous verification instead. Options are the published
Dusart, Ramanujan J. 45 (2018), checking which verification height it uses, or
Broadbent–Kadiri–Lumley–Ng–Wilk, "Sharper bounds for the Chebyshev function θ(x)",
Math. Comp. 90 (2021), which relies on Platt–Trudgian's rigorous RH to 3·10¹². Any
ε ≤ 0.016 for y ≥ 10¹⁰ suffices, and the constants of Cor. 3 do not change.

**REP-2 (Lemma 3(a), contour estimate).** The formula is correct (§2 above), but the
proof cites only a "standard bound". The horizontal segments in the critical strip need
heights T_k → ∞ chosen with |T_k − γ| ≫ 1/log T_k. On those heights
L'/L(σ + iT_k) ≪ log² T_k, uniformly in −1 ≤ σ ≤ 2 (Davenport Ch. 16/19). With the
kernel |s|^{−2} the horizontal integrals are then ≪ x^{c+1} log² T_k / T_k² → 0.

On Re s = −U (U even), the functional equation gives L'/L ≪ log(q|s|). The integral is
then ≪ x^{1−U} ∫ log(|s|)/|s|² → 0 for x > 1.

*Fix:* insert these two lines. The absolute convergence of Σ x^{ρ+1}/(ρ(ρ+1)) makes
the limit over T_k harmless.

**REP-3 (exactness of the decisive Θ_1 check).** Theorem 1 for t ≤ 10⁴ rests on the
sign of Θ_1(p), which is decided in long double. The rounding analysis is plausible, but
the decision is not exact, and the stated margin "0.14 t^{3/2}" is wrong at the bottom
of the range: G(7) = 0.0635 for q = 4.

*Fix:* replace the check with an exact or interval check for p ≤ 10⁴. I have done this
(cert.py: exact big integers for p ≤ 300 and 200-bit arithmetic to 10⁴, all positive,
min Θ_1 = 1.1756 at p = 7 for χ_4 and 2.0794 at p = 5 for χ_3). Put a script of that
kind in the directory and state the true margin.

Similarly, the N_low values used in Corollary 3 on [32, 1.52·10⁸] are long-double
values. Their margin is ≥ 0.1, so the result is safe, but say plainly that this is not
certified, or certify that range in interval arithmetic.

**REP-4 (priority and novelty claims, §2 "Consequence" and §4).**

(a) Lemma 2 with the kernel h'' is the second-order partial summation of Ingham
(pp. 103–104). Johnston uses the same identity as his Lemma 2.7, with the same positive
kernel 1/(t²log²t) + 2/(t²log³t) acting on the integrated error ∫(ψ − t). What is new
here is only that the integrated function Θ_1 is one-signed, so the positive kernel
loses nothing. Claim (c) in §4, "Lemma 2 … appears to be new", should be restated in
that form. Lemma 6 and the cone of Theorem 5 are a genuine use of this observation.

(b) Theorem 6 is the Landau mechanism. Its classical prototype for χ_4 is Hardy–
Littlewood and Landau (1918): GRH(χ_4) is equivalent to Σ_p χ_4(p)e^{−p/x} → −∞. That
work should be cited together with Pintz (1991), whom Johnston credits for the
π − li case.

(c) Liu–Wang–Shi and Platt–Trudgian/Suzuki/Fujii remain unread, so "new" must stay
hedged, as §9 already says.

## MINOR items

1. The cone is written as "m ≤ 0 < m' ≤ 1" in the status paragraph, in §4(c) and at
   the end of §8. Theorem 5 itself, correctly, has m ≤ 0, m < m' ≤ 1, which includes
   m' ≤ 0. Use one form everywhere.
2. In Lemma 4(b'), "θ(x) < x + x/(2 log x) for x > 1 … as quoted in Johnston's
   Lemma 2.2" is inaccurate: Johnston states it for x ≥ 41. This is harmless, since the
   bound is used only for y ≥ 10¹⁰. The RS constant 2.05282 is also given without a
   precise equation number. Moreover, the cited bounds (RS 2.05282, Büthe 1.95) plus
   the log ℓ shift do **not** give c_1 for small y: at y = 4, log 2/√4 ≈ 0.35. Below
   1423, Lemma 4(a) rests on the sieve alone. That is trivially checkable, and I
   checked it, but say so.
3. §5: "inf of N over [11, 10^10] between 0.3880 and 0.3944". The value is exactly
   N(11) = (6/11)·log 11/√11 = 0.394360. On [11, 13), R = 2 − 16/x and N is increasing.
4. `out_q4_1e10.txt` prints the k = 10 bin line twice. This is cosmetic.
5. §2 check (iii): the agreement quoted at 10³ and 10⁴ is uninformative, because the
   allowed error is of order 10²–10³. The χ_3 difference of 69 at 10⁴ is omitted.
   Replace these with small-x checks, which are informative (table in §2 above).
6. Corollary 3 parenthetical: "linearly independent" should be "linearly independent
   and simple" (or "LI in the usual sense"). Also say that liminf N = liminf G, because
   the ∫Θ_1 h'' term is ≥ 0 and is O(1/log x) relative to the main term.
7. Theorem 7 says Corollary 3 holds for χ_3 for every x ≥ 7. That is true, but it would
   also hold from about x ≈ 5.6, since N(5) = 0.432 and R = 2 − 7/x on [5, 7). Optional.

---------------------------------------------------------------------------

## Statement-by-statement verdict

These are proved as written. REP-2 and REP-3 are routine write-ups, and I have already
carried out the computation REP-3 needs.

- **Lemma 2** and identity (2.1) (no hypothesis).
- **Lemma 3(a)**, the explicit formula. The formula is correct and the proof is
  standard; it needs the two-line estimate in REP-2.
- **Lemma 3(b)** (σ_0² majorant, right direction).
- **Lemma 6** (§8).
- **Theorem 1** (χ_4), given the exact check of Θ_1 for p ≤ 10⁴ (REP-3, done here).
- **Theorem 2.** Positivity for every x > 3, and the refined bound for x ≥ 10⁴.
- **Corollary 4** at every x > 3, with the quantitative form for x ≥ 23.
- **Theorem 5**, the full cone m ≤ 0, m < m' ≤ 1, every x > 3.
- **Theorem 6**, the GRH equivalence. Complete, using the elementary positivity of
  L(σ,χ) on the real axis.
- **Theorem 7** (χ_3): Theorems 1, 2, 5, 6 and Corollary 4 for x > 2.

The dependence of Theorems 1, 2, 5, 6 and 7 on Dusart is avoidable via (b'), so REP-1
does not affect them.

These are proved only up to REP-1, which concerns the source certificate for x > 10²⁰:

- **Corollary 3** (χ_4), R(x) > √x/(2 log x) for x ≥ 23. It is fully proved on
  [23, 10²⁰], given the (uncertified but wide-margin) long-double N_low on
  [32, 1.52·10⁸], my check to 10⁷, and the analytic grid. Beyond 10²⁰ it needs a
  rigorous replacement for Dusart's Thm 5.2 (REP-1).
- **Corollary 3** for χ_3 (x ≥ 7) does not need (b), because the fallback limit 0.5245
  is above 1/2. It is proved, with the same caveat about long-double N_low below
  2.52·10⁷.

Not established as stated:

- The novelty claims of §4(c) for Lemma 2, and the implicit priority for the mechanism
  of Theorem 6 (REP-4). These are wording and citation issues, not mathematical ones.
- §8 beyond Lemma 6 is, as the author says, a programme and not a proof. Its Λ values
  are uncertified.
