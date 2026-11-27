# Referee report 1 — the coupling half and the assembled Main Theorem

Manuscript: `monotonicity/monotone-coupling.tex` (2435 lines, file time
2026-11-24 05:42:40 UTC; the PDF alongside is the same build).
Companion: `paper/weighted-races.tex` ("WR").

Scope: Lemma 3.1 (the function Q), Lemma 3.2 (deterministic coupling
inequality), Section 4 (Proposition 4.1, Theorem B, Remarks 4.2–4.3,
Proposition 4.4, Proposition 4.5, Remarks 4.6–4.7), Lemma 2.6 as far as
Theorem B uses it, and Section 8 (Lemma 8.1, certificate statements
(C1)–(C3), Theorem 8.2, Corollary 8.3, Remark 8.4) including every use of
WR. Not in scope: Theorems A and C, and the code of
`monotonicity/cert/certify_derivative.py`.

Everything below was re-derived by hand; numerical statements were re-run
with scripts of my own (Section 5). No file other than this report was
written to the repository.

---

## 1. Verdict

**No FATAL error.** Every lemma, proposition and theorem in scope is
correctly proved as written, with the two qualifications below. The
mathematics of the coupling (Lemmas 3.1–3.2, Proposition 4.1, Theorem B,
Propositions 4.4–4.5) is elementary and I could not break it; the new
differentiability Lemma 8.1 is correct on the whole range M > 0; the three
ranges of the Main Theorem overlap with room to spare.

**Two REPAIRABLE findings**, both about the logical status of the Main
Theorem and not about a proof step:

* **R1.** The Main Theorem is a correct *implication*, not a theorem about
  δ_q under GRH. The abstract, the title and two sentences of the
  introduction present it (and WR Corollary 7) more strongly than that.
* **R2.** Remark 8.4 misdescribes the rigor status of WR Corollary 7 in
  three places, and understates what a direct certificate would remove.
  The repository now contains exactly such a certificate (complete, both
  moduli, up to M = 320 resp. 350), which the manuscript does not use.

**Seven MINOR findings** (stale numbers, a wrong technical term, a proof
by reference to an unstated part of another proof, hypothesis wording).

---

## 2. REPAIRABLE findings

### R1. The Main Theorem is a conditional reduction; the framing overstates it

*What is stated.* Main Theorem (intro) and Theorem 8.2: "assume GRH(χ_q);
suppose (C1)–(C3) hold … and that [WR, Corollary 7] holds as stated there.
Then δ_q is strictly decreasing on (−1/2, ∞)." Corollary 8.3: "Grant the
certificate … and grant [WR, Corollary 7]."

*Status.* As an implication this is rigorous and its proof (three intervals
sharing endpoints; mean value theorem on the middle one) is complete. But
two of its hypotheses are not mathematical hypotheses in the usual sense:
they are "a computation that we did not audit is right" and "a corollary of
another unrefereed preprint, itself stated in a 'computational rigor
class', is right". So what the manuscript proves is

> Theorem B + (certificate) + (WR Cor. 7) ⟹ full monotonicity,

and nothing stronger. Against that:

* the title says "Chebyshev's bias decreases with the weight";
* the abstract says the combination "gives full monotonicity of δ_4 and
  δ_3 on (−1/2, ∞)";
* §1.1 says WR "proves that δ_4 is strictly decreasing for m ≥ 16.82",
  whereas WR itself says "in the computational rigor class described in
  Remark 13";
* the Main Theorem ends "that is, [WR, Conjecture 1] holds", attached to
  the LI sentence. Conjecture 1 of WR is the statement about
  δ_q(m) = P(X > 0) and needs no LI; it is the *first* conclusion.

*Fix.* State the result in two layers.

(a) A reduction with purely mathematical hypotheses (this is Theorem 8.2
with (C3) and the WR citation merged into one hypothesis):

> **Theorem 8.2′.** Let q ∈ {3,4} and assume GRH(χ_q). Let 0 < a ≤ b.
> If Q_{χ_q}(a) ≤ 1, if 𝒟_{χ_q}(M) < 0 for every M ∈ [a, b], and if δ_q is
> strictly decreasing on {m : M ≥ b}, then δ_q is strictly decreasing on
> (−1/2, ∞). Under LI(χ_q) the same holds for the logarithmic density of
> {x : D_{m,b}(x) > 0}, for every real b.

(b) A statement labelled *computer-assisted*, listing what is taken on
trust, in the form "Assume GRH(χ_q). Then δ_q is strictly decreasing on
(−1/2, ∞), **provided** (i) the zero lists (511 / 537 ordinates, ±10⁻²⁰,
complete to height 640.165707 / 700.024747) are correct, (ii) the runs
`results_q4.txt`, `results_q3.txt`, `aux_results.txt` are correct, and
(iii) WR Corollary 7 holds." Replace "gives full monotonicity" in the
abstract by "reduces full monotonicity to", or by "gives, conditionally on
a computer-assisted certificate and on [WR, Cor. 7]," and move "that is,
Conjecture 1 holds" to the first conclusion.

*What would make it a theorem under GRH(χ_q) alone* (LI is needed only to
read δ_q as the density of the race; NCZ is unconditional for q = 3, 4):

1. Use a direct certificate of (C2) up to M ≥ 300. Then WR Proposition 5
   (the evaluated criterion on [17.317, 320] resp. [19.119, 350]) is not
   needed at all, and with it WR Lemmas 22–23 on that range.
   **This certificate now exists in the repository**:
   `cert/results_q4_stretch.txt` (293 cells, [17.4, 320], 0 failed) and
   `cert/results_q3_stretch.txt` (292 cells, [19.2, 350], 0 failed),
   as of commit `40edd9d`. I re-parsed both: cells are contiguous with exact
   4-decimal rational endpoints, all `PASS`, largest upper bound
   −1.224×10⁻⁵ resp. −1.098×10⁻⁵ (consistent with WR's asymptotic
   −D/(2√(2π)σ³) ≈ −1.28×10⁻⁵ at M ≈ 319, q = 4). The manuscript still says
   "we make no claim about such a run".
2. For M ≥ 300 the only remaining input is WR Lemma 24 (analytic leg).
   It needs a referee-level check of its chain — WR Lemmas 14, 15, 16, 17,
   19, 20(iii), 21, 22, 23 — and the logarithmic-derivative constants in its
   proof (−0.15, −0.09, −0.007, with the auxiliary 0.001) must be verified
   in exact or interval arithmetic; WR Remark 13 says they are "still checked in float64 only".
   (I reproduced the M = 300, q = 4 evaluation by hand:
   TAIL_an/MAIN_an ≈ 6.2×10⁻⁹ / 1.36×10⁻⁵ ≈ 4.6×10⁻⁴, matching WR's
   "< 5×10⁻⁴". I did not audit the chain.)
3. An independent audit of `certify_derivative.py` and of the zero-list
   certification (`prototype/explore/certify_arb.py`).

After 1–3 the honest statement is "Theorem (computer-assisted). Assume
GRH(χ_q). Then δ_q is strictly decreasing on (−1/2, ∞)."

### R2. Remark 8.4 misdescribes WR Remark 13

Remark 8.4(3) says: "for M ≥ 300 it is a theorem with hand-checkable
constants, while on 17.32 ≤ M ≤ 320 … it rests on an evaluated criterion,
re-evaluated in ball arithmetic, with three items listed in that remark as
not certified (an endpoint-to-cell step, the GRH-conditional zero count,
and floating-point checks of four constants) … a direct certificate there
would remove the dependence on the first of those items."

I read WR Remark 13, Proposition 5, Lemmas 20–24 and the script
`paper/v3/verify_monotonicity_arb.py`. Three corrections:

1. **The float64 item belongs to M ≥ 300, not to [17.32, 320].** WR:
   "on M ≥ 300 the ball arithmetic verifies the inequality at the single
   point M = 300; the extension to all M ≥ 300 is Lemma 24, whose
   logarithmic-derivative constants are still checked in float64 only."
   So the leg the manuscript calls "a theorem with hand-checkable
   constants" is precisely the one carrying the unverified constants.
2. **The endpoint-to-cell step is not a computational gap and nothing is
   sampled.** WR Proposition 5 is a cell-uniform criterion: D_min =
   max(D(M_a), D(M_b))·(M_a/M_b)³ (WR Lemma 20(ii)), β at M_b, Bessel
   suprema over the argument range [2a_γ(M_a)t_j, 2a_γ(M_b)t_{j+1}], zero
   counts at M_b. The script evaluates exactly that per cell
   (`interval_bounds(Ma, Mb, …)`). WR's remark only says the lemma is a
   human proof that the script does not re-prove. The status on
   [17.317, 320] is therefore "computer-assisted, resting on human-proved
   WR Lemmas 14–23", the same kind of status as (C2).
3. **A direct certificate removes more than "the first item".** On the
   range it covers it removes WR Proposition 5 entirely (items one and
   two). What survives is Lemma 24 and item three.

*Fix.* Replace Remark 8.4(3) by: "[WR, Cor. 7] has two legs. On
17.317 ≤ M ≤ 320 (resp. 19.119 ≤ M ≤ 350) it is computer-assisted: the
cell-uniform criterion of [WR, Prop. 5], whose analytic part is [WR,
Lemmas 14–23], evaluated in ball arithmetic on 148 (resp. 147) cells. On
M ≥ 300 it is [WR, Lemma 24]: an inequality certified in ball arithmetic
at M = 300, extended by a monotonicity argument whose four constants were
checked in floating point only. We audited neither leg." Then either adopt
the stretch certificates (set b_4 = 320, b_3 = 350, replace (C3) by
b_q ≥ 300, cite only [WR, Lemma 24]) or say explicitly that they exist and
are not used.

---

## 3. MINOR findings

**M1. Display (8.2) and the description of (C1) are stale.** The text says
`aux_results.txt` "reports" 3.688482 < M_0(χ_4) < 3.689307 and
4.483177 < M_0(χ_3) < 4.484582, with Q_{χ_4}(3.68) ≤ 0.99439 and
Q_{χ_3}(4.47) ≤ 0.99253. The file now reports
3.688854 < M_0(χ_4) < 3.689059, 4.483812 < M_0(χ_3) < 4.484161,
Q(3.68) ∈ [0.994002, 0.994137], Q(4.47) ∈ [0.991977, 0.992173]. The
printed numbers are the *earlier* run (git history confirms) and coincide
to all digits with the floating-point enclosure from the crude tail bound
(4.3), which I reproduced: [3.6884825, 3.6893068] and
[4.4831777, 4.4845818]. Consequences: (i) no logical harm, the printed
enclosures contain the certified ones and Q(3.68), Q(4.47) < 1 either way;
(ii) "agree with our independent floating-point evaluation" is circular as
printed; (iii) the parenthesis in (C1) describing the method ("in the
sharper form 0 < Q_χ − Q_Γ ≤ (M/Γ)·½Σ…") no longer describes the
certificate, which uses two exact tail moments and Cauchy–Schwarz
(METHOD.md §8; I checked that inequality chain, it is correct, and it does
need the derivative identity [WR, Lemma 20(i)], which the manuscript
lists). Update the numbers and the parenthesis.

**M2. "Comonotone coupling" is the wrong term** (after Proposition 4.4,
and Remark 6.1(3)). The coupling used is the *equal-phase* coupling. The
comonotone (quantile) coupling X_M = F_M⁻¹(U) nests {X_{M'} > 0} ⊂
{X_M > 0} if and only if δ(m) ≥ δ(m'), so for it Proposition 4.4 would be
false whenever δ is monotone beyond M_0 (which the certificate asserts).
Proposition 4.4 is a statement about the equal-phase coupling only; say
so.

**M3. Lemma 8.1 proves pointwise differentiability by citing the *proof*
of [WR, Lemma 21(a)]**, a lemma whose stated hypothesis (eight ordinates
below 2M_a) is the one being removed. The claim that part (a) does not use
that hypothesis is true (I read the proof), but a referee should not have
to. A self-contained argument is five lines: with
P_N(M) = ∏_{j≤N} J_0(c_j(M)t), one has
P_N′(M) = −Σ_{j≤N} c_j′ t J_1(c_j t) ∏_{j'≤N, j'≠j} J_0(c_{j'} t), each term
bounded by c_j c_j′ t²/2, summable uniformly on [M_a, M_b] since
Σ c_j c_j′ ≤ 16 M_b Σγ⁻²; and
|∏_{j'≤N, j'≠j} − ∏_{j'≠j}| ≤ Σ_{j'>N} (c_{j'}(M_b)t)²/4 → 0 uniformly in
M; so P_N′ converges uniformly on [M_a, M_b] and the limit of P_N is
differentiable with the termwise derivative.

**M4. Abstract, first sentence.** "A companion paper showed that under GRH
and LI the logarithmic density … exists": WR Theorem 1 also assumes
L(1/2, χ) ≠ 0 for a general real primitive χ. Add it.

**M5. What δ_χ means without LI.** Theorem B(i)–(ii), Lemma 8.1 and the
first conclusion of Theorem 8.2 are stated under GRH + NCZ only, where
δ_χ(m) := P(X_M > 0) with one independent phase per zero *counted with
multiplicity*. That is a definition of a model quantity; it is the density
of the race only under LI (which forces simple zeros). No proof in scope
uses simplicity, but one sentence should say that without LI these are
statements about the model law. (Lemma 3.2 itself would survive repeated
ordinates sharing a phase, since it allows arbitrary x_j ∈ [−1, 1]; the
formula of Lemma 8.1 and the certificate are specific to independent
phases.)

**M6. Theorem B(iii), "(i) and (ii) are statements about the races
themselves".** True for the densities. It should not be read as nesting of
the actual sets {x : D_m(x) > 0}: the joint law of two races at the same x
is the rotated-phase coupling of Lemma 6.3, not the equal-phase one. The
pair (m, m') = (1, 2) for χ_4 shows the difference: λ < 0, so by Corollary
6.4(b) the violation set has positive density, while M' = 2.5 < M_0 and
the equal-phase coupling nests surely. Remark 6.1(3) says this; a forward reference at Theorem B(iii) would
prevent the misreading.

**M7. Small points.**
(a) Remark 8.4 and the table are consistent with the result files (317 and
296 cells, worst bounds −1.675290×10⁻³ and −1.549494×10⁻³; rounding in the
macros is toward 0, i.e. in the safe direction).
(b) Remark 4.6: 3.298 and 4.103 are valid only with the tail of Σγ⁻³
included (head alone gives 3.29864 / 4.10429, which are upper bounds for
the quantity); with the tail bound I get [3.29821, 3.29864] and
[4.10345, 4.10429], so the printed digits are right.
(c) "nested surely" (§1.3) is literally true only because of the
convention X_M := 1 off Ω_0; "almost surely" is the natural word.

---

## 4. Item-by-item answers

### 4.1 Lemma 3.1 (the function Q) — correct

* Finiteness and continuity: g(t) = 4t³(1+t²)^{−3/2} = 4(1+t⁻²)^{−3/2}
  satisfies 0 < g(t) < min(4, 4t³), so on (0, U] the series is dominated by
  4U³Σγ_j⁻³ ≤ 4U³γ_1⁻¹Σγ_j⁻² < ∞ (Weierstrass). ✔
* Strict monotonicity: each term is strictly increasing and the series
  converges. ✔
* Q(u) → 0 as u → 0⁺ from Q(u) < 4u³Σγ_j⁻³; Q(u) → ∞ because the sequence
  is infinite and each term tends to 4. ✔ (The infinitude is in the
  standing hypothesis of §3 and is supplied in §4 by (2.2).)
* Unique root M_0. ✔
* Σ_j Δ_j = ∫_M^{M'} Q(u)u⁻² du: Δ_j = ∫_M^{M'} 4u(u²+γ_j²)^{−3/2} du and
  Σ_j 4u(u²+γ_j²)^{−3/2} = Q(u)/u²; the interchange is Tonelli for
  nonnegative terms; finiteness from continuity of Q on [M, M']. ✔
  Numerically (χ_4, 511 zeros, (M, M') = (1, 5)): both sides
  0.24771181707259189959… to 30 digits.

### 4.2 Lemma 3.2 (deterministic coupling inequality) — correct

* (i): 0 ≤ 1/γ − 1/√(M²+γ²) ≤ M²/(2γ³), so
  Σ x_j(4/γ_j − 4/√(M²+γ_j²)) converges absolutely and φ(M) converges in
  the given order whenever Σ x_j/γ_j does. ✔
* The hypothesis is the right one and is exactly what the application
  supplies: Ω_0 = {Σ cos θ_j/γ_j converges} has probability 1 by
  Kolmogorov's one-series theorem (independent, centred, variances
  1/(2γ_j²), summable), in the fixed order of increasing γ. On Ω_0 the
  series for X_M converges for **all** M > 0 at once — more than the
  "simultaneously for M and M'" that is needed.
* (ii): the identity and both bounds follow from |x_j| ≤ 1 and Lemma
  3.1(iii). ✔
* (iii), including η = 0: φ(M) ≥ φ(M') + η gives both implications; for
  M < M' ≤ M_0 the integrand 1 − Q is positive on [M, M'), so η > 0,
  **including M' = M_0**. ✔

No independence is used anywhere in §3. The only probabilistic inputs to
Theorem B are P(Ω_0) = 1 and full support of X_{M'}.

### 4.3 Theorem B and Proposition 4.1 — correct

* Probability space: [0, 2π)^ℕ with product Lebesgue measure; X_M = Mφ(M)
  on Ω_0 with x_j = cos θ_j. Each X_M has the law of (1.4). ✔
* (i): on Ω_0, X_M/M ≥ X_{M'}/M' + η; the two events on the left are
  disjoint and inside {X_M > 0}. ✔
* Strictness: needs only P(X_{M'} ∈ (−ηM', 0)) > 0 for the nonempty open
  interval, i.e. full support. This is Lemma 2.6, proved in the manuscript
  itself (I checked every constant: Chebyshev gives 1/4, head within ε/4
  of S, total 3ε/4), and its hypotheses hold for c_γ(M') by (2.2). **No
  absolute continuity, no atom-freeness and no LI are used here.** ✔
* (ii): strictly decreasing on (−1/2, M_0 − 1/2], right endpoint
  included. ✔
* (iii): WR Theorem 1 states — under GRH(χ) and NCZ(χ) — existence of
  ν_m; under LI(χ) in addition, ν_m = law(X_m) (via [ANS14, Thm 1.9]), ν_m
  absolutely continuous with support ℝ, the logarithmic density of
  {R_m > 0} exists and equals ν_m((0, ∞)) ∈ (1/2, 1), continuity in m. WR
  Proposition 1 transfers this to D_m. The manuscript's Theorem 2.3 quotes
  this accurately; the normalisations agree (2a_γ = c_γ = 4M/√(M²+γ²);
  constant 1 from the prime squares).
  - Real zeros: under GRH the only possible real zero in the strip is
    s = 1/2, excluded by NCZ (unconditional for q = 3, 4).
  - Zero at γ = 0: excluded by NCZ; with n_0 > 0 the constant would be
    1 + 2n_0 (Remark 3.3), which only helps.
  - Multiplicities: LI of the *multiset* forces simple zeros.
  - Ties: ν_m({0}) = 0 by absolute continuity, so {D_m = 0} has
    logarithmic density 0 and {D_m > 0}, {D_m ≥ 0} have the same density.
  The extension to all real b is Theorem C(iii) (other referee).
* (iv): see 4.5.

### 4.4 Proposition 4.4 (sharpness) — correct

With cos θ_j < −1 + ε' for j ≤ J:
difference ≤ η + ε'Σ_{j≤J}Δ_j + 2Σ_{j>J}Δ_j ≤ η + |η|/4 + |η|/4 = η/2;
head H' ∈ [h_0, h_0 + M'|η|/4]; tail W' independent of A with full support;
on A ∩ B, 0 < X_{M'} < M'|η|/2 and X_M/M < 0. ✔ The interval statement
follows (an interval not inside (0, M_0] contains two points
M_0 ≤ M < M', where Q > 1 on (M, M')). ✔

The examples η_{χ_4}(1, 5) ≈ 0.55, η_{χ_4}(2, 5) ≈ 0.09 are right: with
the rigorous tail bound, η(1, 5) ∈ [0.552158, 0.552288] and
η(2, 5) ∈ [0.091074, 0.091188].

Caveat on reading numerics: the failure for η < 0 uses arbitrarily high
zeros (the tail must reach −h_0). In a model truncated to 511 zeros it can
be *impossible* for η slightly negative: for (3.6, 3.8), η = −9.1×10⁻⁵ and
a linear program over all phases shows min X_M given X_{M'} = 0⁺ is
+1.8×10⁻⁴. This is a truncation artefact, not a counterexample (with a
free far tail the minimum is −3.3×10⁻⁴). The manuscript's own examples
(3.68, 4.5) and (3.6, 5) do fail in the truncated model (LP minima −0.046
and −0.107), so its sentence "the probability of a violation being
positive but small" is correct for what was simulated.

### 4.5 Proposition 4.5 and Theorem B(iv) — correct

* (a) lower: Q(M) < 4M³Σγ⁻³. Upper: Q(M) > g(M/γ_1), and g(t) = 1 ⟺
  4^{2/3}t² = 1 + t² ⟺ t = κ_0 = (4^{2/3} − 1)^{−1/2} =
  0.811149238535981…; g(κ_0) = 1 to 30 digits. ✔
* (b): Σγ⁻³ ≤ γ_1⁻¹Σγ⁻² = C/(2γ_1). C = (log Λ)″(1/2) = Σ_ρ γ⁻² =
  2Σ_{γ>0}γ⁻². ✔
* (c): termwise 4M³(M²+γ²)^{−3/2} < 4M²/(M²+γ²), and
  Σ 4M²/(M²+γ²) = σ_M²/2. ✔
* (d): partial summation gives Σ_{γ>1}γ⁻³ ≤ C_R′(3/2 log q + 3/4);
  3/4 ≤ 0.69 log q for q ≥ 3 (0.758 ≥ 0.75); Σ_{γ≤1}γ⁻³ ≤ C_R′γ_1⁻³ log q;
  4·(1 + 2.19) = 12.76 ≤ 12.8 in both cases γ_1 ≤ 1 and γ_1 > 1. WR Lemma
  9(b) does state N_χ(T) ≤ C_R′T log(qT) for T ≥ 1 with C_R′ = 2 + 2C_R,
  and C_R = 2.8 follows from [BMOR21, Cor. 1.2] as WR says
  (0.298 + 4.358/log 6 = 2.730). c_1 = 12.8 × 7.6 = 97.28 ≤ 98. ✔

Numerical values (my evaluation, mpmath 30 digits, closed-form C from
ψ′(3/4)/4 + (log L)″(1/2), tail bound (4.3)):

| | χ_4 | χ_3 |
|---|---|---|
| γ_1 | 6.020948904… | 8.039737155… |
| C(χ) | 0.156029649889644884… | 0.113402156852758571… |
| (γ_1/2C)^{1/3} | 2.682106 | 3.284961 |
| (4Σγ⁻³)^{−1/3} | [3.298207, 3.298637] | [4.103454, 4.104290] |
| M_0 | [3.688482, 3.689307] | [4.483178, 4.484582] |
| κ_0 γ_1 | 4.883888 | 6.521427 |
| M_σ | ≥ 1.8239147 (certificate: [1.82391477, 1.82391478)) | ≥ 2.1277508 (certificate: [2.12775095, 2.12775096)) |

All agree with Remark 4.6 and §9(1). Q < σ²/2 checked at M = 0.1, 1, 2,
3.6, 3.68, 4.4, 4.47, 5, 10, 50 for both characters.

### 4.6 Lemma 8.1 (differentiability on the whole range) — correct

Formula: ∂_M c_γ = 4γ²(M²+γ²)^{−3/2} = c_γ′, J_0′ = −J_1, so (8.1) is WR's
(dPhi) with 2a_γ = c_γ. ✔

Dominant: in the γ-th term, |J_1(c_γ t)| ≤ c_γ t/2; at least seven of the
first eight ordinates survive in ∏_{γ'≠γ}, each with
|J_0(c_{γ'}t)| ≤ min(1, (c_{γ'}t)^{−1/2}) ≤ min(1, (c_* t)^{−1/2}) because
c_{γ'}(M) ≥ c_* = c_{γ_8}(M_a) on [M_a, M_b]; the rest are ≤ 1. Hence
|∂_MΦ_M(t)| ≤ (t²/2) D(M) min(1, (c_* t)^{−7/2}), and with
|sin t|/t ≤ min(1, 1/t) the dominant is
(D_max/2)·min(t², c_*^{−7/2} t^{−5/2}), integrable on (0, ∞), uniform on
[M_a, M_b]. D(M) = Σ 16Mγ²/(M²+γ²)² ≤ 16MΣγ⁻². |Φ_M(t)| ≤
min(1, (c_* t)⁻⁴) ∈ L¹, so WR Lemma 7(b) (which requires only Φ ∈ L¹)
gives the inversion formula. ✔ Bessel facts used: |J_1(z)| ≤ |z|/2,
|J_0| ≤ 1, √x|J_0(x)| ≤ √(2/π) (Watson §13.74; WR eq. (37) is indeed that
display). Numerically sup √x|J_0| = 0.797884 and sup|J_1|/(x/2) = 1.

The eight-ordinate hypothesis of WR Lemma 21 is used there only to get
explicit constants (c_γ ≥ 4/√5 for γ ≤ 2M_a). Replacing it by "the first
eight ordinates, with c_* depending on M_a" costs nothing for
differentiability. The modification is legitimate. (See M3 for the
exposition.)

Independent numerical check (511/537 zeros, Gaussian model of the unseen
zeros, Simpson on (0, 60]): the formula (8.1) against a central finite
difference of the inversion integral, and against the certificate's cell
enclosure:

| q | M | formula | finite difference | certificate cell [lo, hi] |
|---|---|---|---|---|
| 4 | 1.0 | −2.349141e−1 | −2.349140e−1 | (below a_4) |
| 4 | 3.6 | −3.342322e−2 | −3.342322e−2 | [−3.3926e−2, −3.2596e−2] |
| 4 | 5.0 | −1.754806e−2 | −1.754806e−2 | [−1.7819e−2, −1.7120e−2] |
| 4 | 10.0 | −4.660368e−3 | −4.660368e−3 | [−4.7404e−3, −4.5546e−3] |
| 4 | 17.4 | −1.706872e−3 | −1.706872e−3 | [−1.7435e−3, −1.6753e−3] |
| 3 | 4.4 | −2.699045e−2 | −2.699045e−2 | [−2.7396e−2, −2.6322e−2] |
| 3 | 8.0 | −8.278043e−3 | −8.278043e−3 | [−8.4076e−3, −8.0780e−3] |
| 3 | 19.2 | −1.575775e−3 | −1.575775e−3 | [−1.6127e−3, −1.5495e−3] |

### 4.7 Theorem 8.2, Corollary 8.3 — correct as implications

* Cell boundaries: Lemma 8.1 gives differentiability on all of (0, ∞);
  (C2) asserts 𝒟 < 0 at every point of the closed interval [a_q, b_q]
  (closed cells sharing endpoints). The mean value theorem gives strict
  decrease on [a_q, b_q]. Continuity of δ′ is not needed. ✔
* Gluing three intervals sharing endpoints. ✔
* Overlaps with the certified numbers:

| | coupling range | certified range | WR range |
|---|---|---|---|
| q = 4 | M ≤ M_0, M_0 > 3.688854 | [3.6, 17.4] | M ≥ 17.317 (stated as ≥ 17.32) |
| q = 3 | M ≤ M_0, M_0 > 4.483812 | [4.4, 19.2] | M ≥ 19.119 (stated as ≥ 19.12) |

  Left overlap 0.0888 resp. 0.0838; right overlap 0.083 resp. 0.081.
  (C1) holds with margin: Q_{χ_4}(3.6) ≤ 0.94213, Q_{χ_3}(4.4) ≤ 0.95333
  already with the crude tail bound. I re-parsed `results_q4.txt` and
  `results_q3.txt`: 317 and 296 cells, contiguous, exact endpoints, all
  `PASS`, all upper bounds negative.
* The definition of δ_q is the same in the manuscript, in the certificate
  and in WR §6 (P(X > 0), GRH only). ✔
* WR Corollary 7 needs 12 certified ordinates below 2M_a at its left edge:
  γ_12(χ_4) = 34.19996 < 34.634, γ_12(χ_3) = 37.55180 < 38.238 (and γ_13
  is above in both cases), which is why its thresholds sit where they do.

Status: see R1 and R2.

### 4.8 Hypothesis hygiene

| hypothesis | where used |
|---|---|
| GRH(χ) | ordinates real, so Q, σ, η, X_M are defined; variance identity and its derivative (certificate tails, (C1)); WR Theorem 1; WR's zero counts |
| NCZ(χ) | bias constant equal to 1 (n_0 = 0); genus-zero Hadamard product behind the variance identity; N_χ = N_f/2 in WR Lemma 9(b) (Prop. 4.5(d)); unconditional for q = 3, 4 |
| LI(χ) | only Theorem B(iii) and the last sentence of Theorem 8.2 / Cor. 8.3: ν_m = law(X_M), hence δ = density of the race |
| simple zeros | never used on its own; implied by LI; see M5 |

Theorem B(i), (ii), (iv), Propositions 4.1, 4.4, 4.5, Lemma 8.1 and the
first conclusion of Theorem 8.2 do not use LI. ✔

---

## 5. Independent Monte Carlo and adversarial check (χ_4)

Own script (numpy, seed 20261124), 511 listed zeros, no tail, 3×10⁵ phase
vectors, all weights evaluated on the same phases. "min" and "max" are the
sample extremes of X_M/M − X_{M'}/M'; the bounds are η and
∫(1+Q)/u² of Lemma 3.2 for the truncated list.

| (M, M') | η | sample min | sample max | upper bound | X_{M'} > 0 ≥ X_M | samples in (−ηM', 0] | δ̂(M) − δ̂(M') |
|---|---|---|---|---|---|---|---|
| (0.3, 3.68) | 2.904231 | 2.915828 | 3.209132 | 3.218957 | 0 | 107808 | 0.35936 |
| (2.0, 3.6) | 0.122074 | 0.129769 | 0.315857 | 0.322370 | 0 | 17747 | 0.09373 |
| (3.0, 3.68) | 0.014049 | 0.017926 | 0.105863 | 0.109140 | 0 | 2234 | 0.02705 |
| (3.5, 3.68) | 0.000915 | 0.002015 | 0.026107 | 0.027035 | 0 | 152 | 0.00640 |
| (1.0, 5.0) | 0.552288 | 0.572713 | 1.030460 | 1.047712 | 0 | 73974 | 0.28921 |
| (3.7, 4.2) | −0.005604 | −0.002263 | 0.067136 | 0.069954 | 0 | — | 0.01410 |
| (4.0, 6.0) | −0.077978 | −0.061388 | 0.230712 | 0.244644 | 101 | — | 0.03663 |
| (5.0, 8.0) | −0.173626 | −0.142284 | 0.297448 | 0.323626 | 1076 | — | 0.03324 |
| (8.0, 12.0) | −0.272765 | −0.214669 | 0.307125 | 0.356099 | 2834 | — | 0.01919 |
| (12.0, 17.4) | −0.347699 | −0.250299 | 0.314020 | 0.399423 | 3448 | — | 0.01296 |

* The two-sided bound of Lemma 3.2 holds in every sample for every pair.
* Sure nesting below M_0: no violation for any pair with η ≥ 0, including
  (1, 5), which straddles M_0.
* Failure above M_0: violations at rates 3.4×10⁻⁴ to 1.15×10⁻² (the
  manuscript's rates for (5, 8) and (8, 12), 3.5×10⁻³ and 9.1×10⁻³, are
  reproduced: 3.6×10⁻³ and 9.4×10⁻³).
* Gap inequality of Theorem B(i): column 8 ≥ column 7 / 3×10⁵ in all five
  cases with η > 0.
* δ̂ keeps decreasing through the pairs with η < 0, as Proposition 4.4
  allows.
* Deterministic check by linear programming over x ∈ [−1, 1]⁵¹¹ (minimise
  X_M subject to X_{M'} = 10⁻⁹): the minimum is negative for (3.7, 3.75),
  (3.69, 3.8), (5, 8) and positive for (3.6, 3.68), (3, 3.689), (3.2, 4),
  (2, 5), (2.5, 5), in agreement with the sign of η except for the
  truncation artefact discussed in 4.4.

The scripts were run from a scratch directory and are not in the
repository.

---

## 6. Certification lists

### Certified as correctly proved as written

* Lemma 2.6 (full support), as used by Theorem B.
* Lemma 3.1 (i), (ii), (iii).
* Lemma 3.2 (i), (ii), (iii); Remark 3.3.
* Proposition 4.1 (pathwise nesting).
* Theorem B (i), (ii), (iv). Theorem B (iii) for b ∈ {0, 1} given WR
  Theorem 1 and Proposition 1 as quoted (which I read and whose statements
  match); for general b it rests on Theorem C(iii), outside my scope.
* Remark 4.2; Remark 4.3 (formal identity, sign and algebra right, flagged
  by the authors as not used).
* Proposition 4.4, including the interval statement.
* Proposition 4.5 (a)–(d); the numbers of Remark 4.6; the exact identity
  and the integral 2λ/π in Remark 4.7 (the rest of 4.7 is labelled
  heuristic).
* Lemma 8.1 (the result; see M3 on the write-up).
* Theorem 8.2 and Corollary 8.3 **as implications**.

### Not certified

* The Main Theorem as an assertion about δ_q under GRH (R1). It is
  conditional on (C1)–(C2) (code not audited by the authors or by me) and
  on WR Corollary 7.
* WR Corollary 7 itself: I read its proof and Remark 13, read the Arb
  script, and reproduced one evaluation; I did not referee WR Lemmas
  14–24.
* Remark 8.4(3) as a description of WR's rigor status (R2).
* Display (8.2) and the (C1) parenthesis as a description of the current
  `aux_results.txt` (M1).
* The term "comonotone coupling" (M2).

### Corrected statements

1. **Main Theorem**: Theorem 8.2′ of R1(a), plus the computer-assisted
   statement of R1(b). With the stretch certificates adopted:
   > Assume GRH(χ_q), q ∈ {3,4}. Grant (C1) at a_4 = 3.6, a_3 = 4.4, grant
   > (C2) on [3.6, 320] (q = 4) resp. [4.4, 350] (q = 3), and grant
   > [WR, Lemma 24] (δ_q′ < 0 for M ≥ 300). Then δ_q is strictly
   > decreasing on (−1/2, ∞); under LI(χ_q) so is the logarithmic density
   > of {x : D_{m,b}(x) > 0} for every real b.
2. **Abstract**: "… reduces full monotonicity of δ_4 and δ_3 on
   (−1/2, ∞) to a ball-arithmetic certificate on a compact range and to
   the large-weight result of the companion paper; this statement is
   conditional on both." First sentence: add "and L(1/2, χ) ≠ 0".
3. **§1.1**: "… and shows, in a computer-assisted sense described there,
   that δ_4 is strictly decreasing for m ≥ 16.82 …".
4. **Remark 8.4(3)**: the replacement text in R2.
5. **(8.2)**: 3.688854 < M_0(χ_4) < 3.689059, 4.483812 < M_0(χ_3) <
   4.484161, Q_{χ_4}(3.68) ≤ 0.994137, Q_{χ_3}(4.47) ≤ 0.992173; delete or
   rewrite the method parenthesis in (C1).
6. **After Proposition 4.4 and in Remark 6.1(3)**: "equal-phase coupling"
   in place of "comonotone coupling", with the sentence "other couplings
   of the same marginals may nest the events on a larger range; the
   quantile coupling does so exactly where δ is monotone."
