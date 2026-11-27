# REFEREE-4: adversarial report on `monotonicity/GRH-ONLY.md`

Date: 2026-11-24. Referee: AI (Anthropic's Claude), acting as adversarial referee in analytic
number theory and probability. Object: `GRH-ONLY.md` (md5 `e314991e579760844b2c86d0c31b8574`,
60865 bytes) with `monotonicity/grhonly/`. The document was not edited. Nothing was committed.

## 0. Verdict

**No FATAL error found.** Every Theorem / Proposition / Lemma of the document was re-derived
from first principles and I found the proofs correct; every certified constant I could
recompute was reproduced. The two points I expected to break, the all-even-moments lemma
(Lemma 1.1) and the non-resonance bound (Theorems 6.1, 6.2), both hold as stated; for
Theorem 6.2 I upgraded the evidence (exact-integer verification of the non-resonance
conditions, Arb evaluation of `J`).

There are **7 REPAIRABLE** items (all at the level of statements around the theorems: one
over-claim in Corollary 3.5, one wrong factor 2 in a literature comparison, hypothesis lines,
status tags, two inaccurate citations) and a list of MINOR items. None changes a number in
the summary table.

What I did not re-prove: the quoted inputs (Q1), (Q4) of the companion paper (I checked that
the quotes agree with `lem:WRexp`, `lem:trunc`, `lem:squares`, `prop:unlogged`, and that the
sieve agrees); the correctness of FLINT/Arb; Platt's computation (tier B).

## 1. What was checked, and how

All computations on one core at `nice -n 10`. My scripts lived in the session scratchpad
(ephemeral); the two that carry rigor are reproduced in the Appendix.

| check | result |
|---|---|
| rerun `certify_arb.py --K 0` (zero lists: sign changes + winding count) | PASSED, 124 s: 511 and 537 zeros, winding balls `[511 ± 3e-12]`, `[537 ± 3e-11]` |
| rerun `certified_bounds.py`, `nonresonance.py` in a scratch copy | outputs identical to the committed `*_results.txt` |
| independent recomputation of §4–5 constants in mpmath (no FLINT): `C(χ)`, `tail_2`, `J_A^+`, head log-sum, `b_δ`, `μ_δ(T)`, `c'`, excess, `C^+`, `σ_M^2`, `τ_M`, `V_M^+`, tier B `E_X` | all agree to 12+ digits, except the two `m = 0` values of `σ^2`, where mpmath's numerical differentiation of `L` at `s = 1` is only good to `3e-10`; for those I used the closed forms `L'/L(1,χ_4) = γ + 2log 2 + 3log π − 4log Γ(1/4)`, `L'/L(1,χ_3) = γ + log 2π − 3log(Γ(1/3)/Γ(2/3))`, which give `σ_{1/2}^2 = 0.15556797992358592886…`, `0.11322996985747234191…`, as in the document |
| Lemma 4.3 against the truth: `J_A` by direct quadrature of `(Λ'/Λ)(s+½)s^{-2}` | `J_30(χ_4) = 0.0659660381`, `J_30^+ = 0.0659679497`, difference `1.91e-6 ≤` stated slack `2.06e-6`; same for `χ_3` |
| Lemma 4.2 at `t = 0, γ_1, 50, 200, 400, 600` (`χ_4`, `δ = 0.9`) | `δ·Re(Λ'/Λ) = 0.125, 1.136, 1.602, 2.421, 2.759, 3.056`; zero sums over the list agree to `≤ 0.021` (the omitted zeros above `T`); bound `μ_δ(t) = 2.01 … 4.49` holds |
| independent sieve to `10^7` (q = 4, 3; four weight pairs; non-prime `x`) | θ-identity residual `≤ 6.5e-14`, exact unlogged identity `≤ 9.1e-14`; brute-force `u`-quadrature at `x = 3000.5` agrees to `8e-7`; `y(E_M − ∫E_{M'}dμ)` bounded |
| (NR1)–(NR3) in **exact integer arithmetic** (ordinates × 10^24) | all margins reproduced; minimum `7.9412e-7` (`χ_3`, NR3, `N ≥ 10`) |
| `J` of Theorem 6.2 in **Arb** at a rational feasible point | `1 − J ∈ [0.884276434474771 ± 5e-16]`, `[0.922718034063003 ± 4e-16]`, `[0.911240824904685 ± 2e-16]`, `[0.944555174897096 ± 5e-16]` |
| closed form of `∫Q^2`, `∫κ^2` on an exact grid, `N = 1,2,3` | agreement `≤ 1.1e-16` |
| Theorem 6.1 on toy models (below, item 7) | never violated; violated when a hypothesis is dropped |
| Lemma 1.1 on the worst resonant configuration (below, item 2) | moments finite, consistent with the bound |
| literature quotes read in the arXiv PDFs: Fiorilli 1210.6946 (Rem. 1.3, Lem. 2.5, 2.7), Devin 1706.06394 (Thm 2.1, 2.2, 2.3, 2.5, Cor. 5.8, Def. 5.3), Bailleul 2007.05763 (Thm 1.7), Martin–Ng 1710.00088 (Def. 1.3, Thm 1.5), ANS 1306.1657 (Thm 1.14), BMOR 2005.02989 (abstract), Platt 1305.3087 (Thm 7.1), Platt–Trudgian 2610.29018 (abstract, §2.1), Hayani 2512.23302 (abstract) | quotes accurate except items R2, R7 below |

## 2. Item-by-item

### Item 1. Resolvent identity (Theorem 2.1, Corollary 2.4, Remark 2.3): correct

Re-derived. `R_m(x) = ∫ t^{-κ} dR_{m'}(t) = x^{-κ}R_{m'}(x) + κ∫_2^x t^{-κ-1}R_{m'}(t)dt`;
multiplying by `2Mx^{-M}` and substituting `t = e^{y-u}` gives exactly
`(M/M')E^θ_{M'} + (1−M/M')∫_0^∞ Me^{-Mu}E^θ_{M'}(e^{y-u})du`. Total mass
`M/M' + (1 − M/M') = 1`; `\hat μ(ξ) = (M/M')(M'+iξ)/(M+iξ)`; `|\hat μ| ≤ 1 ⇔ M ≤ M'`;
`\hat μ(γ) r_γ(M') = r_γ(M)`; multipliers multiply, so the semigroup law holds (also checked
on the density of the convolution, agreement `1e-12`). The sign convention `e^{-iξu}` is the
right one for a shift `y ↦ y − u`. (b): from (1.1) and a split at `u = y/2` (the split needs
`y ≥ 4`; for `2 ≤ y < 4` the bound is trivial; MINOR). (c): signs and the constant `1`
check (`∫_2^x t^{-κ-1}t^{M'}/(2M')dt·2MM'x^{-M} → 1`); the interchange is legitimate because
the truncated formula is uniform and `Σ|γ|^{-2} < ∞`. Corollary 2.4 is consistent with
Theorem A(i) of `monotone-coupling.tex` (same bracket after normalisation).
The exact unlogged form of Remark 2.3 is correct (factor `y/(y−u) = log x/log t`).
Sieve: see table. For `(m,m') = (−0.4, 0)` the quantities `y·Δ` are still growing at `10^7`
(10.4, 13.5 for q = 4, 3); this is the expected `1/M` blow-up of the constant
(`li(x^M)` is not yet asymptotic at `My = 1.6`), not a contradiction.

### Item 2. Limiting distribution, Lemma 1.1, Theorem 1.3: correct

**Lemma 1.1 is proved as written; nothing is assumed.** The mechanism, which answers the
question in the brief:

- `S` has only positive frequencies, so `|S|^{2k} = |S^k|^2` and `S^k` is a generalised
  Dirichlet polynomial with frequencies `ω ≥ kT`.
- The mean value of `e^{i(ω−ω')y}` over `[0,Y]`, `Y ≥ 1`, is bounded in absolute value by
  `min(1, 2/|ω−ω'|)`. Exact coincidences `ω = ω'` (i.e. every additive relation
  `γ_1+…+γ_k = γ'_1+…+γ'_k`, whatever its origin, including multiple zeros) are counted with
  weight 1, the same weight as near-coincidences `|ω−ω'| ≤ 3`. No orthogonality is used.
- The total weight of `k`-fold sums in a unit window is `b_k(n) ≪ log^{2k−1}(n+3)/(n+1)`,
  by induction using only the local count `≪ log` and `|c_j| ≤ A_1/λ_j`. I re-derived the
  induction (the convolution sum is `2H_{n+1}/(n+2)`; the shifted term costs a factor 2).
- Hence the weighted number of solutions is `≤ Σ_n b_k(n)·Σ_{n'} b_k(n')h(n−n')
  ≪ Σ_{n ≥ kT−1} log^{4k−1}n/n^2 ≪ log^{4k−1}T/T`. This is the "absolute convergence after
  pairing" the brief asks about, and it is explicit in the proof.

Stress test: frequencies = integers `n ∈ (T, X]` with multiplicity `⌈log n⌉`, `c = 1/n`, all
phases aligned (maximal resonance compatible with the hypotheses). Then
`S(y) ≈ ½(log^2(1/|y|) − log^2 T)` near `y = 0`, and the moments converge as `X → ∞` but
slowly (the integrand peaks at `1/|y| ≈ e^{4k}`): for `T = 100`,
`T·mean|S|^4/log^7 T = 0.20, 0.38, 0.47, 0.50` at `X = 2^{12}, 2^{15}, 2^{18}, 2^{21}`.
Two consequences worth recording in the document:
(i) `C_k` is of factorial size (`≈ (4k)!`); irrelevant for the theorems.
(ii) In this model `P(S > A) ≍ e^{-\sqrt{2A}}`, so *exponential* moments are infinite: the
class of functions covered by Lemma 1.1 does not have exponential moments in general. This
matches the Rubinstein–Sarnak / Devin tail `exp(−c√A)` and bears on Corollary 3.5 (R1).

Theorem 1.3(a): correct given (Q1) (`4Y^3e^{-Y/2}` checked). (b): the `W_p` argument is
correct (coupling through the same `y`; Cauchy in `T`; identification of the limit in `W_2`;
Villani Thm 6.9). Lemma 1.2 is correct (characters of `A` are indexed by their frequency;
ergodicity; full support of Haar measure). The Hausdorff–Young route of Remark 1.5 is
correct and gives `F_M ∈ L^p` independently of Lemma 1.1.

**Dependency map** (answers "what survives if Lemma 1.1 failed"): Theorem 3.1 for *every*
convex `φ` does **not** depend on Lemma 1.1. It needs only `F_M ∈ L^1(A)`, the multiplier
identity, and the identification of `ν_M` as the weak limit (B^2 theory, `k = 1`). Lemma 1.1
with `k ≥ 2` is used only for the assertion that *time averages* of `φ(E_M(e^y))` converge
for `φ` of polynomial growth of degree `> 2`.

### Item 3. Convex order (Theorems 3.1, 3.4, Corollary 3.5): theorems correct, one over-claim

Theorem 3.1: the Fourier computation is right (`∫Φ\bar ψ = \hat μ(ω)\hat F_{M'}(ψ)`, both
functions supported on `{1, e_γ^{±1}}`), Fubini is justified by `‖F_{M'}‖_1 < ∞`, Jensen is
applied to a probability measure with `g ∈ L^1(μ)` for a.e. `θ`, and Tonelli after
subtracting a supporting line is legitimate. Because the proof is on the group, there is no
`φ(E + o(1))` issue at all; the `O(1/y)` error only enters through Theorem 1.3(b), where it
is harmless for polynomial growth.

Proposition 3.2: correct. At `θ = ι(∓1/T)` the cosine part is `≤ σ_M^2/2` (resp. `≥ 0`),
all sine terms have one sign for `γ < T`, and the window `[T/4, T/2]` contributes
`≫ (N(T/2) − N(T/4))/T ≫ log T`. Lemma 3.3 and Theorem 3.4: correct (equality in Jensen
forces `φ = ℓ` on the essential range of `u ↦ F_{M'}(τ_{-u}θ)`, which is unbounded on both
sides for a.e. `θ`; zero set of a nonnegative convex function is an interval).

**R1 (REPAIRABLE). Corollary 3.5 over-claims for the exponential functional, and its last
sentence is vacuous as worded.**
(a) The preamble says each listed functional "is the limit of the corresponding logarithmic
time average". For `∫e^{tx}dν_M` this is not proved: Theorem 1.3(b) covers polynomial
growth only. Moreover, under GRH alone finiteness of `∫e^{tx}dν_M` is not known for any
`t ≠ 0` (the available tail bound is `exp(−c√A)`, and item 2(ii) shows that this order of
tail does occur in the abstract class), so the bullet may be empty.
(b) "any bound of the form `ν((−∞,0]) ≤ ∫Ψdν` with `Ψ ≥ 1_{(−∞,0]}` convex that holds at
weight `M'` holds at every `M < M'`": such an inequality holds for every measure, so nothing
is asserted.
*Fix.* Replace the preamble by: "Each of the following is a strictly increasing function of
`M`; the first three are the limits of the corresponding time averages (Theorem 1.3(b))."
Add to the last bullet: "(a statement about `ν_M` only; under GRH alone we do not know that
this integral is finite for any `t ≠ 0`)". Replace the last sentence by: "Moreover, for
convex `Ψ ≥ 1_{(−∞,0]}` and `M < M'`: `ν_M((−∞,0]) ≤ ∫Ψdν_M ≤ ∫Ψdν_{M'}`; so a bound for
`ν_M((−∞,0])` obtained by transporting a convex test function from a heavier weight is never
better than the same test function at `M`." (This is what Remark 3.7 uses.)

### Item 4. Both lead sets (Proposition 3.2, Theorem 3.6): correct

`{D_m < 0} = {E_M < 0}` for `x > 1`; `(−∞,0)` is open; portmanteau gives
`liminf ν_{M,Y}((−∞,0)) ≥ ν_M((−∞,0)) > 0`; the change from `(1/Y)∫_0^Y dy` to
`(1/log X)∫_2^X dt/t` costs `O(1/Y)`. The quantitative part is correct
(`m_γc_γ/2 = |\hat F_M(e_γ)| ≤ ‖F_M‖_1 = 1 + 2∫F_M^-`; Cauchy–Schwarz with
`(F^-)^2 ≤ (F−1)^2` on `{F<0}`); `γ_1/√3 − ½ = 2.976` for q = 4.

### Item 5. Variance with multiplicities, Lemmas 4.2, 4.3, Propositions 4.4, 4.5: correct

Theorem 4.1: correct (`‖F_M − 1‖_2^2 = ½Σ m_γ^2|r_γ|^2`).
Lemma 4.2: derived from the Hadamard identity, Binet's formula
(`0 ≤ ½ − 1/t + 1/(e^t−1) ≤ t/12`, `Re(1/2w) > 0`) and `|L'/L(s)| ≤ −ζ'/ζ(Re s)`; I
re-derived `b_δ` (the `−½log 2` from `log|w|` turns `log(q/π)` into `log(q/2π)`) and
checked the bound numerically (table). `μ_{0.9}(T_4) = 4.5219`, `μ_{0.9}(T_3) = 4.4327`. So
"multiplicity ≤ 4 just above `T`" is a consequence of Lemma 4.2 under GRH + NCZ and does
not use any external zero-count.
Lemma 4.3: `∫_A^∞ 2ds/(s(s^2+γ^2)) = γ^{-2}log(1+γ^2/A^2)` and both elementary integrals
are right; verified against quadrature (table) and, through the head sums, against the zero
lists: the exact log-tail `0.014645819` (`χ_4`) agrees with the Riemann–von Mangoldt
estimate `0.014652` to `4e-4` relative.
Proposition 4.4: correct (`A ≥ 2x_0` gives `log(γ^2+4x_0^2) ≤ log A^2 + log(1+γ^2/A^2)`;
the `c'τ_M` term is an identity whatever the sign of `c'`). Proposition 4.5: correct
(`δb_δ − 1 = 0.614, 0.485 > 0`). Recomputed: `V_{1/2}^+(χ_4) = 0.16963083`,
`V_{1/2}^+(χ_3) = 0.12545141`, `C^+ = 0.1700925058, 0.1256235992`, tier B
`0.15557363, 0.11323421`, `E_X = 2.825e-6, 2.119e-6`. Thresholds `m < 0.7220`,
`0.9184`, `0.7778`, `0.9947` and `M_σ = 1.82391, 2.12775` confirmed.

**R2 (REPAIRABLE). The comparison with Bennett–Martin–O'Bryant–Rechnitzer is wrong by a
factor 2.** The quoted inequality is correct (arXiv:2005.02989, abstract:
`|N(T,χ) − (T/π·log(qT/2πe) − χ(−1)/4)| ≤ 0.22737ℓ + 2log(1+ℓ) − 0.5`,
`ℓ = log(q(T+2)/2π) > 1.567`, `T ≥ 5/7`). But `N(T,χ)` counts zeros with `|γ| ≤ T`, and for
real `χ` an ordinate `γ` of multiplicity `m_γ` makes `N` jump by `2m_γ`. At `T_4` the bound
is `4.763`, so `2m_γ ≤ 9.53`, i.e. **`m_γ ≤ 4`, not `≤ 9`** (`T_3`: `4.660`, again
`m_γ ≤ 4`). Hence the sentence "under GRH the Poisson-kernel bound above is sharper" is not
supported: the integer bounds coincide and BMOR is unconditional (as real numbers `4.52`
vs `4.76`).
*Fix.* "(… gives `2m_γ ≤ 9.53`, i.e. `m_γ ≤ 4`, unconditionally, at the same height; Lemma
4.2 gives the same integer bound, is self-contained under GRH, and has the smooth form
needed in Proposition 4.4.)" Nothing downstream changes.

Tier B is correctly labelled as conditional. Platt's Theorem 7.1 is quoted correctly
(`max(10^8/q, 7.5·10^7/q + 200)` for even `q`, `3.75·10^7/q + 200` for odd `q`); the
Platt–Trudgian phrase "(and are simple)" is in their §2.1.

### Item 6. Cantelli (Theorems 5.1, 5.2, 5.3, Proposition 5.4): correct

Direction of Cantelli (`P(Z − EZ ≤ −t) ≤ Var/(Var + t^2)`, `t = 1`), portmanteau for the
closed set `(−∞,0]`, ties counted with the loser, `E_M ⋚ 0 ⇔ D_m ⋚ 0` for `x > 1`, the
stretch `y < log 2` contributing `O(1/Y)`: all right. `1/1.1696308345 = 0.854970620`,
`1/1.1254514113 = 0.888532361`. Theorem 5.3: correct (orthogonality gives `V'`;
`{F ≤ 0} ⊂ {U ≤ −(1−B)}`); bracketed table values reproduced (`0.999464…`, `0.996494…`,
`0.978593…`, `0.963911…`, `0.903337…`). Proposition 5.4: correct
(`B ≤ c_1M log^2(G+3)`, `V' ≤ c_2M^2 log^2(G+3)/G`, the choice of `G`).

**R3 (REPAIRABLE, wording). "Under LI the truth is `exp(−c/M^2)`-small."** The cited
Theorem 1(iii) of the companion paper is an *upper bound* (sub-Gaussian Chernoff). The true
size under LI is far smaller: since `Σ c_γ(M) = ∞` with `c_γ ≈ 4M/γ`, a Montgomery-type
large-deviation computation gives `log(1/P) ≍ exp(c M^{-1/2})` up to powers of `M`
(heuristic, not proved here). *Fix:* "Under LI one has the upper bound
`exp(−1/(8C(χ)M^2))` (`weighted-races.tex`, Thm 1(iii)); the true order is expected to be
doubly exponential in `M^{-1/2}`."

MINOR: "At `m = 0` it gains nothing" is true in tier A only (`c_{γ_1} = 0.3310 < 2V^+ =
0.3393`); with simple zeros `2V = 0.3111 < 0.3310` and soft thresholding gains slightly.
The first-order criterion "gain iff `a > 2V`" is exact (`dR/ds|_0 = (2V − a)/(1+V)^2`).

### Item 7. Non-resonance bound (Theorems 6.1, 6.2): correct; evidence upgraded

**What is claimed.** A deterministic inequality on the group `A`:
`1_{F ≤ 0} ≤ (1 − κF)^2` with `κ = κ_0(1 − Σb_n cos ϑ_n) ≥ 0`, expanded into six integrals.
Characters of `A` are indexed by their real frequency, so "the integral of a monomial is 1
iff its total frequency is 0" is exact. I re-derived every item:

1. `∫Q^2`: formal degree 4 in the head; (NR1) makes it equal to the independent-phase value.
   `E G = −½Σab`, `E G^2 = ½Σ(a−b)^2 + E[L^2H^2]` (the cross term has odd degree),
   `E[L^2H^2] = ¼ΣaaΣbb + ½(Σab)^2 − ⅜Σa^2b^2` (re-derived by pairings; `¼+½−⅜ = ⅜`).
2. `∫QκW_𝒦 = ∫QκW_t = 0`: `Qκ` has degree 3; its frequencies `|n·γ| ≤ 3γ_N ≤ T` cannot be
   an unseen ordinate, and (NR2) excludes the listed ones.
3. `∫κ^2W_𝒦^2 = (∫κ^2)V_𝒦`: (NR1) and (NR3), including the diagonal `γ' = γ''`
   (the script includes `i = j` in the sums).
4. `∫κ^2W_𝒦W_t = 0`: for `γ'' > T` the size argument `γ'' ≤ 2γ_N + T' = T` is right; for
   `γ'' ∈ (T', T]` (NR3) applies.
5. `∫κ^2W_t^2 ≤ κ_0^2(1+Σb)^2V_t`: sup-norm only.

**How the unseen zeros are handled.** They are never tested. They enter (i) through
frequencies that are too large to match any head monomial of degree `≤ 3` or any
`n·γ ± γ'` with `|n|_1 ≤ 2`, `γ' ≤ T'`, and (ii) through the single number `V_t`, bounded by
`V_M^+ − V_𝓗 − V_𝒦` (Proposition 4.4), in a term where no orthogonality is used. Unseen zeros
may resonate with the head at order `≥ 4`, with `𝒦`, or among themselves, with any
multiplicity; the proof does not care. This is sound.

**How a numerical check excludes an exact relation.** The relations to exclude form a
finite list of linear forms in the *certified* ordinates, each with coefficient mass
`≤ 4` (NR1: `|n|_1 ≤ 4`; NR2: `≤ 3 + 1`; NR3: `≤ 2 + 2`). Each ordinate is enclosed with
radius `10^{-20}` and the list is complete up to `T` (winding count, rerun by me), so each
form is known to within `4·10^{-20}`. I recomputed all forms in exact integer arithmetic on
the decimal midpoints (Appendix A): the smallest absolute values are

| | NR1 | NR2 | NR3 |
|---|---|---|---|
| `χ_4`, `N = 8` | 3.1380e-2 | 8.6878e-3 | 2.5389e-5 |
| `χ_4`, `N = 16` | 6.4631e-3 | 6.9545e-4 | 3.6214e-6 |
| `χ_4`, `N = 20` | 2.7660e-3 | 2.3755e-4 | 1.4769e-6 |
| `χ_3`, `N = 10` | 8.4219e-3 | 7.1150e-4 | 7.9412e-7 |
| `χ_3`, `N = 20` | 8.5570e-4 | 9.5385e-5 | 7.9412e-7 |

identical to the double-precision output. `7.9·10^{-7} − 4·10^{-20} > 0`: the conditions
are proved, given the certification of the lists.

**Are 0.8842 / 0.9112 rigorous?** Yes. The inequality is a theorem, the hypotheses are
verified as above, any feasible `(κ_0, b)` gives a bound, and I evaluated `J` in Arb at a
rational feasible point with `V_M^+` rounded up (Appendix B):

    χ_4  tier A  N =  8:  J ∈ [0.115723565525229 ± 5e-16]   density ≥ 0.8842764
    χ_4  tier B  N = 16:  J ∈ [0.0772819659369966 ± 4e-17]  density ≥ 0.9227180
    χ_3  tier A  N = 10:  J ∈ [0.0887591750953149 ± 5e-17]  density ≥ 0.9112408
    χ_3  tier B  N = 20:  J ∈ [0.0554448251029045 ± 4e-17]  density ≥ 0.9445551

**R5 (REPAIRABLE, status).** The tag "[rig-fp]" for Theorem 6.2 can be replaced by
[cert-A] / [cert-B] once `nonresonance.py` evaluates `J` in Arb and performs the comparisons
of (NR1)–(NR3) on balls or exact integers (30 lines; Appendix). As it stands the tag is
honest and the numbers are correct.

**Toy-model tests** (the theorem is abstract, so it can be tested on any frequency set
satisfying the hypotheses):

- *LI-like*: `χ_4` amplitudes, `N = 4`, 196 independent `𝒦`-phases, independent tail of
  variance `0.0177`, `2·10^6` samples: `P(F ≤ 0) = 0.0057`, `E(1−κF)^2 = 0.0986`,
  `J = 0.1139`; `E Q^2` and `E κ^2` match the closed forms to `2e-5`, `5e-4`.
- *Adversarial resonant tail, exact quadrature on a `1024^2` grid*: head `γ = (1, √2)`,
  `T = 3γ_2`, tail = an arbitrary function of the two head phases with spectrum in
  `{n : |n·γ| > T}` (so the tail is a deterministic function of the head, the opposite of
  LI), obtained by projecting indicators of `{1 + H < t}` on the allowed spectrum; 4
  amplitude pairs × 4 tail variances × 6 thresholds × 2 shapes × 2 signs.
  `E(1−κF)^2 ≤ J` in every case; `max P(F ≤ 0)/J = 0.51`.
- *Adversarial with nonempty `𝒦`* (6-torus Monte Carlo, `10^6` samples, 60 random
  configurations; tail characters such as `e_{K_3}e_{K_4}`, `e_{γ_2}^9`,
  `e_{γ_1}e_{K_3}e_{K_4}`, all of frequency `> T`): `E(1−κF)^2 ≤ J` within `0.64` standard
  errors; `max P/J = 0.35`.
- *Independent tail with the Cantelli-extremal two-point law*: `max P/J = 0.56`.
- *Controls (hypothesis dropped)*: with an order-3 head resonance in the tail,
  `E(1−κF)^2` exceeds `J` (by `9.9e-4` in the worst case found); with a tail frequency
  `γ_1 + K_1 ≤ T`, the cross term `E κ^2W_𝒦W_t` is nonzero at 94 standard errors. So the
  tests can see a failure, and (NR2), (NR3) are genuinely used.

MINOR: "decidable by comparing enclosures" should be "verifiable when true" (an exact
relation could not be detected by enclosures). "can resonate with nothing" should be "with
nothing that enters items 1–4". The check "closed formula for `J` against exact grid
quadrature, `N = 1,2,3`" announced in §10 is not in any script of `grhonly/`; I did it
(agreement `≤ 1.1e-16`). Remark 6.3.3: I get `0.509` for the LI value of `‖F−1‖_4`, not
`≈ 0.52`. Remark 6.3.5: the sign criterion `σ_M^2/2 > M^2Σ(M^2+γ^2)^{-2}/Σ(M^2+γ^2)^{-1}`
is right (re-derived).

### Item 8. Theorem 7.1, Proposition 7.3: correct, hypothesis line to fix

Theorem 7.1(a): inclusion from Corollary 2.4; `σ_Mσ_{M'} < 2 ⇒ λ_χ > 0` by Cauchy–Schwarz
(`Σ4/(|M+iγ||M'+iγ|) ≤ σ_Mσ_{M'}/(2MM')`). (b): the interlacing argument is right
(`δ^≥(m+) ≤ δ^>(m) ≤ δ^≥(m)`, equality at continuity points; superadditivity of lower
density on eventually disjoint tie sets). Proposition 7.3: correct (trivial annihilator
gives `A = 𝕋 × A'`; arcsine factor; weak continuity in `M` plus no atom at 0).

**R4 (REPAIRABLE).** (i) Theorem 7.1 is headed "(GRH; no LI)", but its proof cites Corollary
2.4, proved here under GRH + NCZ, and part (c) uses Theorem 5.1 (GRH + NCZ). Either head it
"GRH + NCZ" or cite Theorem A of `monotone-coupling.tex` for (a), (b) (no NCZ, with the
`n_0`-modified margin) and keep NCZ for (c). (ii) Summary row 7.1 says "ties are immaterial
outside a countable set of weights"; the theorem proves this only for weights in an
interval `I` on which `λ_χ > 0` (e.g. `M ≤ M_σ(χ)`). Add "in the Theorem-A range".
(iii) `I` is an interval of `M` in the hypothesis and of `m` in the conclusion.

**R7 (REPAIRABLE, citations).** (i) Devin Thm 2.2 (hypothesis
`⟨λ_1..λ_N⟩_ℚ ∩ ⟨Z^* − {λ_i}⟩_ℚ = {0}`; conclusion: `μ_S` has no atoms and `δ(S)` exists)
and Bailleul Thm 1.7 are quoted fairly; Devin proves continuity, the absolute continuity in
Proposition 7.3 is the document's own (correct) addition. (ii) "Martin–Ng need three
self-sufficient zeros for their (many-contestant) statement" is inaccurate: Martin–Ng
Theorem 1.5(a) is a *two-way* statement (three self-sufficient ordinates ⇒ weakly
inclusive). (iii) The definition is attributed to Martin–Ng, but their Definition 1.3 is on
the *multiset* of ordinates and forces the zero to be simple; the document's definition (on
the set `Γ` of distinct ordinates, multiplicity allowed) is Devin's Definition 5.3.

### Item 9. Erratum for the companion paper: correct

The companion remark applies ANS14 Thm 1.14 with the ordinates listed with multiplicity and
gets `1 + ½Σ_γ|r_γ|^2`; the correct value is `1 + ½Σ_{distinct} m_γ^2|r_γ|^2 = 1 + V_m`.
The replacement text is accurate (`σ_m^2 = 4M(Λ'/Λ)(m+1)`, equality iff all zeros simple,
`C_2 ≤ 0.170093, 0.125624` rounded up correctly, Fiorilli Lem. 2.5 and Devin Thm 2.1 do
have the squared multiplicities; I read both). Suggested additions (MINOR):

- ANS14 state their theorems for a **non-decreasing** sequence `(λ_n)` and Thm 1.14 as
  `c^2 + ½Σ|r_n|^2` (arXiv v1 eq. (1.22); the companion cites "(1.24)", check the journal
  numbering). So the companion's "literal" application was faithful to the source; the
  formula is valid only after merging equal frequencies. Say so, to make clear where the
  slip comes from.
- Existence of `ν_m` (their Cor. 1.3) needs no change of convention; the last sentence of
  the erratum is harmless but not needed.
- `RESEARCH-CONTEXT.md` §2 also says the dissolution law and monotonicity "need only GRH +
  NCZ for the *limiting distribution*"; they concern the model variable `X_m`, which is the
  limiting distribution only under LI. Same correction.

### Item 10. Novelty and positioning: accurate and modest

The table credits Fiorilli and Devin for the variance and for the Chebyshev-type density
bound, RS94/Kaczorowski/Devin for positivity at `m = 0`, and claims as new the convex
order, the explicit multiplicity control, the evaluation, and the non-resonance bound. I
found nothing in the papers I read that contradicts these claims; Devin's Cor. 5.8 is
indeed `δ_lower ≥ 1 − Var/m_S^2` and would give `0.844` for `χ_4` if the zeros are simple.
The header claim on NCZ (Platt–Trudgian: even `χ`, `q ≤ 10^9`; odd `χ`, `q ≤ 10^{10}`) is
accurate; "if odd" should read "if `χ` is odd". Hayani's result is a natural-density
statement, as the document says.

## 3. Findings by class

**FATAL:** none.

**REPAIRABLE** (fixes given in §2):

- **R1** Corollary 3.5: time-average claim not proved for `e^{tx}` (and finiteness unknown
  under GRH alone); last sentence vacuous as worded.
- **R2** Lemma 4.2, parenthesis on BMOR: the bound gives `m_γ ≤ 4`, not `≤ 9` (`N(T,χ)`
  counts `±γ`); "sharper" is unsupported.
- **R3** After Proposition 5.4: `exp(−c/M^2)` is an upper bound under LI, not "the truth".
- **R4** Theorem 7.1: hypothesis line (NCZ), and the summary row on ties must be restricted
  to the Theorem-A range.
- **R5** Theorem 6.2: status tag can and should be upgraded to ball arithmetic (done here).
- **R6** Summary table row "1.4" refers to a theorem that does not exist (it is Theorem
  1.3); the scripts' docstrings use an older numbering (`sieve_checks.py`: "Theorem 1.1",
  "Remark 1.3", "Theorem 2.1" for the convex order; `certified_bounds.py`: "Theorem 5.2"
  for the soft threshold). A reader cannot match script to statement.
- **R7** Citations: Martin–Ng (two-way, multiset definition).

**MINOR:**

1. Proof of Theorem 2.1(b): the split at `u = y/2` needs `y ≥ 4`.
2. `c' := δ(log A/2 + b_δ) − 1` is ambiguous; it is `δ((log A)/2 + b_δ) − 1` (`= 2.1447`,
   `2.0152`).
3. "At `m = 0` it gains nothing": tier A only.
4. `certified_bounds.py` restricts the soft threshold to `t > μ_δ(T)·4M/T`; this is not
   needed and costs a little (`m = −0.45`, `χ_4`: `0.999498` instead of `0.999464`).
5. "decidable", "can resonate with nothing": wording (item 7).
6. Grid-quadrature check of `J` not in the scripts; `0.52` should be `0.51`.
7. Theorem 6.1 should say explicitly that the list of ordinates in `(0, T]` is complete
   (it is used in (NR2), (NR3)); this is what part [D] of `certify_arb.py` provides.
8. Remark 1.5 / Theorem 3.6: RS94 Thm 1.2 bounds the mass outside a ball; cite Devin
   Thm 2.5(1) for the two one-sided statements.
9. Numbering gaps (no 1.4, 7.2).
10. Erratum: the three additions of item 9.

## 4. Statements certified as proved as written

Lemma 1.1; Lemma 1.2; Theorem 1.3 (given (Q1), (Q3)); Theorem 2.1(a)–(d); Remark 2.3
(exact unlogged identity); Corollary 2.4; Theorem 3.1; Proposition 3.2; Lemma 3.3; Theorem
3.4; Theorem 3.6; Theorem 4.1; Lemma 4.2; Lemma 4.3; Proposition 4.4; Proposition 4.5;
the numerical table of §4 [cert-A] and the tier-B values (conditional as stated); Theorem
5.1; Theorem 5.2(a)–(c) (`0.8549`, `0.8885` under GRH alone; `0.8653`, `0.8982` under the
stated simplicity hypothesis); the table of `1/(1+V_M^+)` and soft-threshold values; Theorem
5.3; Proposition 5.4; Theorem 6.1; Theorem 6.2 (`0.8842`, `0.9112` under GRH alone;
`0.9227`, `0.9445` under the tier-B hypothesis), now with exact verification of
(NR1)–(NR3) and Arb enclosures of `J`; the §6 table of other weights (reproduced);
Theorem 7.1(a), (b) under GRH and (c) under GRH + NCZ; Proposition 7.3; the erratum of §9.

## 5. Statements not certified as written, with corrected statements

- **Corollary 3.5.** Corrected: "Each of `V_M`, `∫|x−c|^p dν_M` (`p ≥ 1`),
  `∫(c−x)^+dν_M`, `∫(x−c)^+dν_M` is strictly increasing in `M` and is the limit of the
  corresponding time average. For real `t ≠ 0`, `M ↦ ∫e^{tx}dν_M ∈ (0,∞]` is nondecreasing,
  and strictly increasing on `{M < M'}` whenever it is finite at `M'` (finiteness is not
  known under GRH alone; no time-average statement is made). For convex
  `Ψ ≥ 1_{(−∞,0]}` and `M < M'`: `ν_M((−∞,0]) ≤ ∫Ψdν_M ≤ ∫Ψdν_{M'}`."
- **Parenthesis after Lemma 4.2.** Corrected: BMOR gives `2m_γ ≤ 9.53`, hence `m_γ ≤ 4`
  unconditionally at `T_4` (and at `T_3`); Lemma 4.2 gives the same integer bound.
- **Sentence after Proposition 5.4.** Corrected: upper bound `exp(−1/(8C(χ)M^2))` under LI;
  true order expected doubly exponential in `M^{-1/2}`.
- **Theorem 7.1, heading and summary row.** Corrected: (a), (b) under GRH (via Theorem A of
  `monotone-coupling.tex`) or GRH + NCZ (via Corollary 2.4); (c) under GRH + NCZ; the tie
  statement holds for all but countably many weights *in an interval on which
  `λ_χ(M,M') > 0`*.
- **Sentence after Proposition 7.3.** Corrected: Martin–Ng, Theorem 1.5(a): three
  self-sufficient ordinates (multiset sense) imply that the two-way race is weakly
  inclusive.
- **Theorem 6.2 status.** Corrected: [cert-A] / [cert-B] with the enclosures of item 7.

## Appendix A. Exact-integer verification of (NR1)–(NR3)

```python
from decimal import Decimal
import itertools, bisect
Z = 'prototype/explore/'; S = 10**24
def load(fn): return [int(Decimal(l.strip())*S) for l in open(Z+fn) if l.strip()]
def vecs(N, order):
    seen = set()
    for k in range(1, order+1):
        for idx in itertools.combinations_with_replacement(range(N), k):
            for sg in itertools.product((1,-1), repeat=k):
                v = [0]*N
                for i, s in zip(idx, sg): v[i] += s
                if any(v): seen.add(tuple(v))
    return seen
def nearest(L, s):
    j = bisect.bisect_left(L, s)
    return min(abs(L[jj]-s) for jj in (j-1, j) if 0 <= jj < len(L))
for fn in ('chi4_zeros.txt', 'chi3_zeros.txt'):
    g = load(fn); T = g[-1] + S//4
    for N in (1,2,3,4,5,6,8,10,12,16,20):
        head, rest = g[:N], g[N:]; assert 3*head[-1] < T
        dot = lambda v: abs(sum(a*b for a, b in zip(v, head)))
        V4 = vecs(N, 4)
        m1 = min(dot(v) for v in V4)
        m2 = min(nearest(rest, dot(v)) for v in V4 if sum(map(abs, v)) <= 3)
        D = sorted([rest[i]+rest[j] for i in range(len(rest)) for j in range(i, len(rest))]
                 + [rest[j]-rest[i] for i in range(len(rest)) for j in range(i+1, len(rest))])
        m3 = min(nearest(D, dot(v)) for v in V4 if sum(map(abs, v)) <= 2)
        print(fn, N, m1/S, m2/S, m3/S)      # each must exceed 4e-20
```

Runtime 3 s. Every linear form has coefficient mass `≤ 4`, so the enclosure radius
`10^{-20}` moves it by at most `4·10^{-20}`.

## Appendix B. Arb evaluation of `J`

With `g` = ordinates as `arb(mid) ± 1e-20`, `M = 1/2`, `a_n = 4M/sqrt(M^2+g_n^2)`,
`𝒦 = {i ≥ N : g_i < T − 2g_{N−1}}` decided on balls, `V_t = V^+ + 10^{-10} − V_𝓗 − V_𝒦`
with `V^+` the decimal upper bounds `0.1696308345`, `0.1555736302`, `0.1254514113`,
`0.1132342076`; `(κ_0, b)` = the float minimiser with each `b_n` rounded **down** to 6
decimals (so `Σb_n ≤ 1` exactly) and `κ_0` rounded to 6 decimals, converted to exact
rationals; then

    J = (1−κ_0)^2 + (1−κ_0)κ_0·Σab + κ_0^2 [ ½Σ(a−b)^2 + ¼ΣaaΣbb + ½(Σab)^2 − ⅜Σa^2b^2
                                             + (1+½Σbb)V_𝒦 + (1+Σb)^2 V_t ]

evaluated at 200 bits. Results in item 7.
