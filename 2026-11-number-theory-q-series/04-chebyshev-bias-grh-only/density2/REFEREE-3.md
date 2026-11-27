AI-produced (Anthropic's Claude), adversarial referee report, 2026-11-27.

# Referee report 3: `TOY.md` §1 and §2.2, `HYP.md`

**Scope.** Refereed in full: `TOY.md` §1 (Lemma 1.1, Theorem 1.2, Lemma 1.3, Proposition 1.4
and the remarks after them) and §2.2 (setting, Theorem 2.3 with steps (i)–(v) and Remarks
(1)–(4), Corollary 2.4 with Steps 1–5, Remarks on Corollary 2.4). Also refereed: what the
introduction (lines 8–20), the end of §2.1 ("What is proved") and the "Proved instances"
paragraph of §5 say about these results; and `HYP.md` in full. Inputs read for definitions:
`DENSITY.md` §1, §3, Theorem 2.1, Remark 7.5 and the statements of §7; `LITERATURE.md`; the
author's `toy_toeplitz.py`, `toy_toeplitz2.py` and their outputs.

**Not read**, as instructed: `REFEREE-1-TOY.md`, `REFEREE-2.md`, `referee1/`, `referee2/`,
`PROGRESS.md`. Every citation of these in `TOY.md` is treated as unverified. Where such a
citation carries mathematics or a number, it was redone here from scratch.

**Summary.** I found **no mathematical error** in any labelled result. I checked every step of
Lemma 1.1, Theorem 1.2, Lemma 1.3, Proposition 1.4, Theorem 2.3 (i)–(v), Corollary 2.4 (Steps
1–5), and HYP Proposition 1.1, Lemma 1.2, Theorem 1.3 and Corollary 1.4 by hand. I also checked
the identities of Theorem 2.3 (i) numerically with an independent method, and the Galois step
(iv) in exact arithmetic. Every attempt to build a counterexample failed, for the reasons the
proofs give (§3 below). The problems found are of three kinds:

- **Scope.** Three places say "hierarchical (Toeplitz)" or "Toeplitz weights" where only the
  one-hole-at-0 class of §2.2 is covered: the end of §2.1, the heading of §2.2, and the last
  sentence of §5.
- **One numerical remark.** Remark (1) on Corollary 2.4 says that periodic symmetric patterns are
  flat for special `M`. This is false for some patterns, e.g. `m_n = 1_{n odd}`. The remark
  also quotes 108,205 flat triples; an independent, validated count gives **108,164**.
- **Two small gaps and some wording.** The gaps are in Theorem 1.2 (i) (passing from the coset
  ring to an integer combination of coset indicators) and in HYP §2 item 4 (the word "all").
  The wording points are misattributions and terminology.

## 1. Verdict table

| # | Item | Verdict |
|---|---|---|
| 1 | TOY Lemma 1.1 (mod-2 lemma): statement and proof | CORRECT |
| 2 | TOY remark after Lemma 1.1 (the identity in `(ℤ/3)²`) | CORRECT |
| 3 | TOY Theorem 1.2: statement | CORRECT |
| 4 | TOY Theorem 1.2, step (i) (structure of `2w`, parity of `P̂`) | CORRECT WITH GAP |
| 5 | TOY Theorem 1.2, step (ii) (orbits, BV representative, derivative) | CORRECT |
| 6 | TOY Theorem 1.2, step (iii) (flat set ⇒ `P ≡ c`, `c` odd) | CORRECT |
| 7 | TOY Lemma 1.3 (a) (BV along a.e. orbit, invariant set `A_0`) | CORRECT |
| 8 | TOY Lemma 1.3 (b) (disintegration (1.1), monotone class) | CORRECT |
| 9 | TOY Lemma 1.3 (c) (densities, `λ_A(Z) = 0`, `ρ̃ = ρ`) | CORRECT |
| 10 | TOY Proposition 1.4: decomposition, no dense cosets, (a) | CORRECT |
| 11 | TOY Proposition 1.4 (b) (true race) | CORRECT |
| 12 | TOY Remarks (1)–(3) after Proposition 1.4 | CORRECT (two attributions to tidy; optional fixes 5, 6) |
| 13 | TOY §2.2 setting: definition and the equivalence with "`P_k`-periodic off `P_kℤ`" | CORRECT |
| 14 | TOY §2.2 setting: `ν_2+ν_3` example of a regular Toeplitz sequence outside the class | CORRECT |
| 15 | TOY §2.2 setting: examples (valuation family = Proposition 2.1; periodic sequences) | CORRECT |
| 16 | TOY Theorem 2.3: statement (0/1 patterns and the integer-pattern extension) | CORRECT |
| 17 | Theorem 2.3 step (i): both forms of `w'`, the `p_a` criterion, atom masses, (R) | CORRECT |
| 18 | Theorem 2.3 step (ii): periodicity of `Φ` and the `L²` bound | CORRECT |
| 19 | Theorem 2.3 step (iii): density points and the Markov count | CORRECT |
| 20 | Theorem 2.3 step (iii): symmetry (`A_r = \bar A_{P−r} = −A_r`) | CORRECT |
| 21 | Theorem 2.3 step (iii): `B(s) = B_*` and `V_{i+1} = β_*` on `≥ 0.6(q−1)` residues | CORRECT |
| 22 | Theorem 2.3 step (iv): Galois, `β_* ∈ ℤ` | CORRECT |
| 23 | Theorem 2.3 step (v): parity, `κ` odd | CORRECT |
| 24 | Theorem 2.3 Remarks (1) and (3) (interpretation, scope) | CORRECT |
| 25 | Theorem 2.3 Remark (2) (numerical checks) | CORRECT as a description (referee-2 citations unverified; reproduced independently here; optional fix 8) |
| 26 | Theorem 2.3 Remark (4) (integer patterns) | CORRECT |
| 27 | Corollary 2.4: statement and "In particular" | CORRECT |
| 28 | Corollary 2.4 Step 1 (cases (a), (b), (c)) | CORRECT |
| 29 | Corollary 2.4 Step 2 (`v_{i+1} = −β_*1_{≠0}`, `β_* ∈ {0,−1}`) | CORRECT |
| 30 | Corollary 2.4 Step 3 (the bound) | CORRECT |
| 31 | Corollary 2.4 Step 4 (the level `±1`) | CORRECT |
| 32 | Corollary 2.4 Step 5 (periodicity; `κ = 1 − 2m_P` for every period) | CORRECT |
| 33 | Remarks on Corollary 2.4, (1) (periodic flat pieces; 108,205 triples) | **WRONG AS STATED** (minor) |
| 34 | Remarks on Corollary 2.4, (2) (proportion → 1) | CORRECT |
| 35 | Remarks on Corollary 2.4, (3) (open question; Ramanujan sums) | CORRECT (optional fix 9) |
| 36 | TOY introduction, lines 8–20 | CORRECT (wording; optional fixes 1, 3, 4) |
| 37 | TOY end of §2.1, "What is proved" | **WRONG AS STATED** (scope) |
| 38 | TOY heading of §2.2, "Toeplitz weights obey the same parity law" | **WRONG AS STATED** (scope) |
| 39 | TOY §5 "Proved instances": first two sentences | CORRECT |
| 40 | TOY §5 "Proved instances": the "ghost mechanism … hierarchically (Toeplitz)" sentence | **WRONG AS STATED** (scope) |
| 41 | HYP Proposition 1.1 (the bound `1/N`) | CORRECT |
| 42 | HYP Lemma 1.2 | CORRECT |
| 43 | HYP Hypothesis TS_ℤ and the monotonicity remark | CORRECT |
| 44 | HYP Theorem 1.3 (both cases) and its consequences | CORRECT |
| 45 | HYP Corollary 1.4 | CORRECT |
| 46 | HYP remark after Corollary 1.4 | CORRECT |
| 47 | HYP §2 item 1 (TS ⇒ TS_ℤ) | CORRECT |
| 48 | HYP §2 item 2 (TS_ℤ ⇏ TS: `X = {2^kα}`) | CORRECT |
| 49 | HYP §2 item 3 (TS_ℤ and FD incomparable: `{2^kα}`; Example 3.3a) | CORRECT |
| 50 | HYP §2 item 4: "the hypotheses LITERATURE.md lists … all imply FD" | CORRECT WITH GAP |
| 51 | HYP §2 item 4: parity exception (`X = {(2k+1)α}`), the GRH-only items, "strictly weaker than FD", the TS-implying list | CORRECT |
| 52 | HYP §2 item 5 (proof mechanism) | CORRECT (interpretation) |
| 53 | HYP "Status of the hypotheses" | not checkable here (literature status); consistent with `LITERATURE.md` |

No item is WRONG (an error in the mathematics).

## 2. Items that are not CORRECT

### 2.1 Scope of Theorem 2.3: rows 37, 38, 40

Quoted text:

> *What is proved:* … For hierarchical (Toeplitz) weights the parity law of Theorem 1.2 holds,
> by a different argument (Theorem 2.3). (end of §2.1)

> ### 2.2 Toeplitz weights obey the same parity law

> a "ghost" mechanism … is excluded by Theorem 2.3 when the blocks are organised
> hierarchically (Toeplitz), but not for non-hierarchical concatenations of periodic blocks (§5)

*Problem.* Theorem 2.3 is proved only for the class of §2.2: one-sided sequences that are
`P_k`-periodic off `P_kℤ` along one chain `P_k | P_{k+1}`, i.e. one hole per period, located at
the residue `0`. Step (i) uses that there is a single hole class (so `m'` is `P`-periodic).
Step (ii) uses that the hole is at `0`: the remainder `F''` has its frequencies in `Pℤ`, so
`Φ` is the same function on every sub-arc, and steps (iii)–(v) rest on this. TOY's own Remark
(3) after Theorem 2.3 lists hierarchical (regular Toeplitz) sequences that are not covered:

- holes at `P_k/2`, where the remainder is antiperiodic;
- holes in several classes, e.g. `m_n = [ν_2(n)+ν_3(n) odd]`, whose level-`k` remainder has
  phases that differ between adjacent sub-arcs. By the setting paragraph of §2.2 this sequence
  is a regular Toeplitz sequence whose only non-periodic position is `0`.

So "hierarchical (Toeplitz) weights" and "Toeplitz weights" claim more than is proved. The
introduction (lines 13–14) and the first sentence of §5 "Proved instances" state the class
correctly. The three places above should use the same wording.

### 2.2 Remarks on Corollary 2.4, item (1): row 33

Quoted text:

> (1) The periodic case does occur. For periodic symmetric patterns, `F_M` is flat at level
> `1 − 2m_Q` on some sub-arcs for special values of `M`. Referee 2 found 108,205 such (pattern,
> sub-arc, `M`) triples among the 11,956 primitive symmetric patterns with `Q ≤ 22`, all at
> levels `±1` …

*The formula.* Let `m` be symmetric with minimal period `Q`, `p(r) := m_r` (so `p(0) = m_Q`),
`W(s) := Σ_{r mod Q} p(r) cos(2πrs/Q)` and `ε := 2πM/Q`. On the sub-arc `J_s` of length `2π/Q`,

    F_M = (1 − 2m_Q) + 2B(s)e^{−εu},   B(s) = ε f_s(e^{−ε})/(1 − e^{−Qε}),   f_s(x) := Σ_{j=0}^{Q−1} x^j W(s−j)

(this is (R) of TOY (i) with the full mass `W`, solved periodically). So `F_M` is flat on `J_s`
iff `f_s` has a root `x = e^{−ε} ∈ (0,1)`, and the level is then `1 − 2m_Q`.

*(a) "For periodic symmetric patterns" is false as a universal statement.* Counterexample:
`m_n = 1_{n odd}` (`Q = 2`). Here `W = (1, −1)`, `f_0 = 1 − x` and `f_1 = x − 1`. The only root
is `x = 1`, i.e. `M = 0`, so `F_M` is never flat for `M > 0`. The same holds for
`m_n = 1_{n even}` (`f_s = 1 + x`) and for both primitive symmetric patterns with `Q = 3`. In
all, 78 of the 11,956 primitive symmetric patterns with `2 ≤ Q ≤ 22` have no flat sub-arc for
any `M > 0`, and 11,878 have at least one (`r3_periodic_flat.py`).

*(b) The count 108,205 is not reproduced.*
- The pattern count 11,956 is reproduced exactly (`2 ≤ Q ≤ 22`; `Q = 1` adds the two constant
  sequences).
- I counted, for each (pattern, sub-arc), the distinct `M > 0` with `B(s) = 0`. This gives
  **108,164** triples: 56,285 at level `−1` and 51,879 at level `+1`.
- The classifier was validated in 60-digit arithmetic, with 0 disagreements, on 1,600 random
  (pattern, sub-arc) pairs and on all 48 pairs where a naive count disagrees with it
  (`r3_periodic_flat_hpcheck.py`).
- A naive double-precision count (plain `numpy.roots`, `|Im| < 1e−7`) gives 108,212. It is
  inflated by two artefacts:
  - when `W(s) = 0` exactly (an algebraic zero), rounding to `~1e−16` creates spurious roots
    `x ≈ 10^{−8}…10^{−5}`;
  - exact double roots are split in two. Example: `m = 001100001000011000` (`Q = 18`), `s = 8`
    has the double root `x = 2cos(4π/9)`, where `f_s` and `f_s'` vanish to 120 digits. It
    gives one value of `M`.
- The quoted 108,205 lies between these two counts. I did not read the script that produced
  it, so I cannot say which artefact, if either, it contains.

The level claim (only `±1`, indeed exactly `1 − 2m_Q`) is confirmed, and it is automatic from
the formula. I also confirmed eight flat cases with an independent evaluation of `F_M`: the
deviation from the level on `J_s` is `≤ 1.2e−14`.

*(c)* "`Q ≤ 22`" should be "`2 ≤ Q ≤ 22`". For `Q = 1`, `m ≡ 0` gives `F ≡ 1`, which is flat
everywhere for every `M`.

### 2.3 Theorem 1.2, step (i): row 4 (gap, minor)

Quoted text:

> By Cohen's idempotent theorem (…), each level set of `2ŵ` lies in the coset ring of `G`.
> Hence `2ŵ = Σ_i n_i 1_{g_i+L_i}` (a finite sum, `n_i ∈ ℤ`).

"Hence" needs one more fact: the indicator of every member of the coset ring is a
`ℤ`-combination of indicators of cosets. *Missing argument.* Let `𝓜` be the `ℤ`-span of coset
indicators. It is closed under products, because two cosets meet in a coset or in `∅`. So the
family `{S : 1_S ∈ 𝓜}`:
- contains all cosets;
- is closed under `∩` (products);
- is closed under `∪`, by `1_{A∪B} = 1_A + 1_B − 1_A1_B`;
- is closed under `∖`, by `1_{A∖B} = 1_A − 1_A1_B`.

It is therefore a ring of sets containing the cosets, hence contains the coset ring. Then
`2ŵ = Σ_{v∈V} v·1_{2ŵ=v}` has the stated form. The same sentence is used in Proposition 1.4.

### 2.4 HYP §2 item 4: row 50 (gap in a universal claim)

Quoted text:

> The hypotheses that `LITERATURE.md` lists for continuity of the law or for existence of the
> density all imply FD, and so imply the disjunction.

`LITERATURE.md` also lists Akbary–Ng–Shahabi 1306.1657 ("continuity only under LI-type
hypotheses", tagged MEMORY). It does not say what those hypotheses are, and HYP does not
mention them, so "all" is not verified for this item.

A second point: Martin–Ng Theorem 1.7(a) is about zeros of `ζ` and the race `π(x)` versus
`li(x)`. FD or TS_ℤ for the ordinates of `ζ` makes sense as a condition on a set of
frequencies. But Theorem 1.3 here and Theorem 3.4 of `DENSITY.md` are proved for the race
`D_m` of a real primitive `χ`. So "and so imply the disjunction. It implies continuity of
`ν_M`" does not apply to that race as proved.

HYP already makes this kind of caveat for hypotheses stated for several `L`-functions. It
should extend it to these two items, or replace "all" by the explicit list.

## 3. Notes on the items judged CORRECT, including the counterexample attempts

**Lemma 1.1.** `Δ_ℓ` corresponds to multiplication by `[−ℓ]+[0]` in `𝔽_2[G]`. This ring has no
zero divisors because `G ⊂ ℝ` is totally ordered: the top coefficient of a product is the
product of the top coefficients. `Q` is finitely supported, and the congruence is between
functions on `G`, possibly not finitely supported. The `Δ`'s are applied before reducing, so
this is fine. The `(ℤ/3)²` identity also checks: `0` lies on all 4 lines, and every other point
on exactly one.

**Theorem 1.2.**
- (i) The Lagrange polynomials are taken in `M(A)` with unit `δ_e`. The measure with transform
  `1_{g+L}` is `e_g m_{L^⊥}`, and `λ_A(L^⊥) = 0` when `L ≠ 0`, because `L^⊥` has infinite index.
  Lemma 1.1 is applied to `Q = 1_{0} − P̂`.
- (ii) `h = F − P̃` satisfies `(D+M)h = Mσ`. Lemma 1.3 then applies with `ρ = −Mh`. The
  right-continuous representative is unique, hence covariant.
- (iii) At a.e. `t ∈ Z_θ`, the point is a density point of `Z_θ ∩ {F̃_θ = c}` and `F̃_θ` is
  differentiable there, so the derivative is `0`. Then `P = c` a.e. on `S`, and `DENSITY.md`
  Lemma 1.1 applies on any `A` of §0: the image of `π` is closed and has trivial annihilator.

**Lemma 1.3.**
- `Du_ε = φ_ε⋆μ` holds as distributions, and both sides are finite measures that agree on
  trigonometric polynomials.
- The variation is a countable supremum of `L¹_loc`-continuous functionals, so Fatou applies.
- In (b), the domination by `|ν_θ|(supp ψ + [−1,1])` is integrable by (a).
- In (c), the singular part of `ν_θ` is carried by `{D|ν_θ| = ∞}`, a Lebesgue-null set. The
  `limsup` over rational `δ` equals the full `limsup`, because `δ ↦ |ν|((−δ,δ))` is monotone.
- The uniqueness of the Lebesgue decomposition gives `ρ̃ = ρ`.

**Proposition 1.4.**
- `supp D` is locally finite.
- The "superincreasing" choice of the `ℓ_i` makes `ℓ_T ≠ 0` for `T ≠ ∅`. That is all the
  contradiction needs; distinctness is not used.
- (b): `m` bounded and supported in finitely many APs plus a finite set gives
  `Σ_{0<γ≤T} m_γ = O(T)`, against `N(T) ~ (T/2π)log T`.

**§2.2 setting.**
- The equivalence holds as stated. The converse needs only well-definedness, which is
  `P_{k+1}`-periodicity off `P_{k+1}ℤ` applied to `P_k r'`.
- The `ν_2+ν_3` sequence is `P`-periodic off `Pℤ` for no `P ≥ 2`. I checked this for
  `P = 2^a`, `P = 3^b`, `P` with a prime factor `≥ 5`, and `P = 2^a3^bc`.
- Period doubling is the valuation family with `q_k = 2`, `c_k = 1_{k even}`.

**Theorem 2.3, step (i).**
- *First form.* From `1/(1−e^{iθ}) = (1+i cot(θ/2))/2` we get
  `w' = 1/2 + Σ_r p(r) sin((P/2−r)t)/sin(Pt/2)`.
- *Second form.* Substitute `t = 2π(s+u)/P`. The factors `(−1)^s` cancel, and the pair
  `r, P−r` of the symmetric part contributes `z + \bar z`, which is real.
- *Atom masses.* The masses are `(2π/P)Σ_k (P/P_k)1_{(P/P_k)|s}V_k(sP_k/P mod q_k)`. For
  `q ∤ s`, only the term `k = i+1` survives.
- All of this was verified against an independent evaluation of `F'`
  (`r3_toeplitz_forms.py`; method below): error `≤ 2.4e−13` for both forms; `(R)` to
  `2.2e−16`; the general mass recursion, at all `s` including `q | s`, to `4.9e−15`.

**Theorem 2.3, step (ii).** `Σ_{n∈Pℤ∖0} 4M²/n² = (4π²/3)M²/P²`, and `m² = m` for 0/1 weights.

**Theorem 2.3, step (iii).**
- *Markov count.* `#bad · 0.1λ(J) < λ(I∖E) ≤ 0.01qλ(J)`, so `#bad < 0.1q`.
- *Adjacent good pair.* `2⌊0.1q⌋ < q−1` for every `q ≥ 2`.
- *Agreement of `F'`.* Real-analyticity of `F'` on open sub-arcs (a distributional ODE with an
  analytic right side) turns agreement on a set of measure `≥ 0.8` into agreement on all of
  `(0,1)`.
- *Symmetry.* `A_{P−r} = −\bar A_r` uses that `p_a` is odd and real, and
  `1 − e^{2πir/P} ≠ 0` gives `p_a ≡ 0`.
- *Count of residues.* `(q−1) − 2⌊0.1q⌋ ≥ 0.6(q−1)` holds for all `q ≥ 2`.

**Theorem 2.3, steps (iv)–(v).**
- (iv) The fibres of `a ↦ ae_0` are cosets of the stabiliser, so on an orbit every Galois
  conjugate of `β` occurs equally often, and a value taken on `> 1/2` of an orbit is
  rational.
- (v) `ρ ∈ [ε/(e^ε−1), ε/(1−e^{−ε})] ⊂ [1−ε, 1+ε]`, and both cases `|n| ≤ d/(8ε)` and
  `|n| > d/(8ε)` check.
- Remark (4): `∫|Φ|² = 8M²Σ_{j≥1} m_{jP}²/(M²+j²P²)`, which is a tail of a convergent
  series. Nothing else in the proof uses `m ∈ {0,1}`.

**Theorem 2.3: counterexample attempts.**
- *Non-symmetric patterns.* The decisive identity is
  `w'_{s−1} − w'_s = −Im[Σ_r A_re^{iω_ru}]/sin πu`. I verified it numerically against the
  independent evaluation, in 6 non-symmetric cases with `P ≤ 84`, to `≤ 3.2e−14`. In those
  cases, the smallest value over all adjacent pairs of `sup_u|w'_{s−1} − w'_s|` is `0.162`.
  Over 3,374 adjacent pairs from random non-symmetric periodic patterns with `Q ≤ 22` it is
  `0.158`. Adjacent sub-arcs never agree, so no level set of positive measure exists. The
  attempt fails for the reason the proof gives.
- *Levels that are not odd integers.* (v) excludes them once `β_* ∈ ℤ`, and (iv) needs only
  that a non-rational value cannot fill 60% of the residues. An exact check
  (`r3_galois.py`) covered 124,522 patterns: all symmetric 0/1 patterns with `q ≤ 30`, all 0/1
  patterns with `q ≤ 14`, and all `{0,1,2}`-valued patterns with `q ≤ 9`. It found 0
  violations of (iv) and of Corollary 2.4 Steps 1(b) and 2. A non-rational value never
  occupies more than 50% of `{1,…,q−1}`, which is the orbit bound.
- The theorem does not exclude odd `|κ| ≥ 3`, and does not claim to. At a single level the
  arithmetic already allows `κ = −3` at `q = 20`: a symmetric 0/1 pattern has `V = −2` on 12 of
  the 19 nonzero residues. This is smaller than the `q = 21` example of Remark (3).

**Corollary 2.4.**
- *Atoms when the periods are bounded but `m` is not periodic.*
  - Step 1(a) forces `B(s) = B_*` on all `q` sub-arcs of a dense arc, hence `V_{i+1}` constant
    on the nonzero residues.
  - Step 2 (exact Fourier inversion; checked exhaustively as (C2) in `r3_galois.py`) then
    gives `v_{i+1} = −β_*1_{≠0}`.
  - So a bounded-`q` sequence with an atom has eventually constant patterns, and Step 5 makes
    it periodic. There is no room for a counterexample.
  - As an illustration only, not a proof, small-ball scans (`r3_smallball.py`) bracket the
    true `F_M` between the depth-`K` approximant and the bound (ii). They give bounded values
    of `λ{|F_M − κ| < η}/(2η)` down to `η = 10^{−4}`, at `κ ∈ {0, ±1, ±3}`, for:
    - period doubling (`P = 2^22`, `M = 0.5, 1, 2`);
    - three random bounded-`q` sequences;
    - a near-periodic bounded-`q` sequence (`1_{≠0}` except a non-constant `q = 4` pattern
      every fifth level).
  - The same scan detects the known atom of a periodic control: the ratio grows like
    `0.2/(2η)`, i.e. mass `1/5`.
- *Level `±1`.* Step 4 is right: `|κ − 1 − 2β^{(i)}| → 0` with `β^{(i)} ∈ {0, −1}` forces
  `κ ∈ {±1}` and `β^{(i)}` eventually constant. Step 5 is right as well: the identity
  `m_P = m_{PP'} = m_{P'}` uses both periods.
- *"In particular".* For a valuation family, `m` is periodic iff `(c_k)` is eventually
  constant: if `Q` is a period, `k(P_K + Q) = max{j : P_j | Q}` for large `K`. So the
  statement is equivalent to "no atoms unless `m` is periodic".
- *Remark (2).* Checked. The two cases `|n| ≤ 1/(2ε)` and `|n| > 1/(2ε)` give `|Φ| ≥ 1` on a
  good `U_s`.
- *Remark (3).* `c_q(e) = μ(q)` on the units, and `c_21 = 1` on 12 of 20 residues. Both
  checked exactly.

**Proposition 2.1 (outside the assigned items, used by the introduction).** The jumps of `F_M`
for period doubling, `4πM/3` at `0` and `π` and `±(4πM/3)2^{−a}` at `π/2, π/4, π/8`, are
reproduced by the depth-22 recursion to within `2·10^{−6}` (`r3_smallball.py`). That is the
size of the jump the truncation leaves in the remainder, `(4πM/3)4^{−11}`.

**HYP Proposition 1.1.**
- `K = H^⊥` has order `N`.
- `G = H + ⟨Λ⟩` makes the characters `e_γ`, `γ ∈ Λ`, separate the points of `K`.
- Each `D_{k,k'}` is a nonzero trigonometric polynomial in the distinct characters `e_{±γ}`.
- So `|S(θ)| ≤ 1` a.e. and `Nν_M({c}) ≤ 1`.
- *Counterexample attempt for the bound `1/N`.*
  - The bound cannot be improved in this generality. In an abstract cyclic model (`A = 𝕋`,
    `H = Nℤ`, `Λ = {1,2}`, `F_H := c − F_Λ` on `[0,2π/N)`, extended periodically), the level
    set has measure exactly `1/N`, one point per `K`-orbit (`r3_hyp.py` (a)).
  - The bound fails if `H + ⟨Λ⟩ ≠ G`: with `H = 4ℤ` and `Λ = {2}` one gets `1/2 > 1/4`
    (`r3_hyp.py` (b)).
  - In HYP, `G = H + ⟨Λ⟩` holds because `G = ⟨Γ⟩` and `Γ ⊂ H ∪ Λ`. So the attempt fails.

**HYP Lemma 1.2, Theorem 1.3, Corollary 1.4.**
- Lemma 1.2: if TS fails, `Λ ⊂ ⟨Γ_{>T}⟩_ℚ ⊂ ⟨Γ∖Λ⟩_ℚ` for `T ≥ max Λ`, and a finitely generated
  torsion abelian group is finite.
- Theorem 1.3, case 1 is the contrapositive of Lemma 1.2 with `Λ = Γ_T`. Case 2 uses
  `Γ∖H_T ⊂ Γ_T`, which is finite.
- Theorem 1.3's consequences: these need only `ν_M({c}) = 0` for every `c`, through (1.1) and
  (2.3) of `DENSITY.md`, for every `m > −1/2`. The interval `I` of `DENSITY.md` §2 is not
  needed.
- Corollary 1.4: the index is a nondecreasing integer-valued function of `T`, bounded by
  `1/α`, so it is eventually constant; nested subgroups of equal finite index are equal; and
  `ord(γ) | |G/H_∞|`.

**HYP §2.**
- Item 1: a nonzero form vanishing on `Γ_{>T}` maps `G/⟨Γ_{>T}⟩` onto a nonzero subgroup of
  `ℚ`.
- Item 2: `⟨X_{>T}⟩ = 2^{k(T)}ℤα`, and `⟨X⟩_ℚ = ℚα`.
- Item 3:
  - `{2^kα}` fails FD for every `N_0`, `ε`: once `#X_T > N_0`, the set `X_T∖Λ_T` is
    nonempty and spans `ℚα`.
  - Example 3.3a satisfies Devin's 5.4(1): the only elements of `X ∩ (0,T]` that involve a
    primary `b_p ∈ (log T, T^{1/2−ε}]` are `b_p` itself, since `b_p + b_{N(p,j)} ≥ e^{b_p} > T`
    and primary indices are never linking.
  - Example 3.3a fails TS_ℤ: `b_n = (b_n + b_{N(n,j)}) − b_{N(n,j)}`.
- Item 4, parity exception: the basis `e_k − (2k+1)e_0` of the relation lattice spans the
  integer kernel. The minor is `±1`; checked with `K = 12` in `r3_hyp.py`. Its coordinate sums
  are `−2k`, which are even, and `gcd(2k+1, 2k+3) = 1`.

## 4. Required fixes

1. **Scope of Theorem 2.3** (rows 37, 38, 40). Replace "hierarchical (Toeplitz) weights" (end
   of §2.1) and "Toeplitz weights" (heading of §2.2) by "Toeplitz weights with one hole per
   period, located at 0 (the class of §2.2)". In §5, say that the ghost mechanism is
   excluded "when the blocks are organised hierarchically with one hole per period at `0`
   (the class of §2.2)". Also say that it is not excluded for the other Toeplitz hierarchies
   of Remark (3) after Theorem 2.3, nor for non-hierarchical concatenations.
2. **Remarks on Corollary 2.4, (1)** (row 33).
   - (a) Replace "For periodic symmetric patterns" by "For most periodic symmetric patterns
     (11,878 of the 11,956 primitive ones with `2 ≤ Q ≤ 22`; none with `Q ≤ 3`, e.g.
     `m_n = 1_{n odd}` is never flat)".
   - (b) Replace "108,205" by a count whose method is stated. For example: 108,164 distinct
     `M` over all (pattern, sub-arc), exact double roots counted once (`r3_periodic_flat.py`,
     validated by `r3_periodic_flat_hpcheck.py`). Alternatively, drop the exact figure.
   - (c) Write `2 ≤ Q ≤ 22`.
3. **Theorem 1.2 (i)** (row 4). After "lies in the coset ring of `G`", add the two-line
   argument of §2.3 above: the `ℤ`-span of coset indicators is closed under products, so it
   contains the indicators of all members of the coset ring.
4. **HYP §2 item 4** (row 50). Either replace "all" by the explicit list, or add that
   Akbary–Ng–Shahabi's "LI-type hypotheses" are not stated in `LITERATURE.md`, and that
   Martin–Ng Theorem 1.7(a) concerns `ζ` and `π(x)` versus `li(x)`, outside the setting in
   which Theorem 1.3 and `DENSITY.md` Theorem 3.4 are proved.

## 5. Optional improvements

1. *Terminology.* The `q_k` are ratios of consecutive periods; the periods are the `P_k`, which
   are never bounded and (for `k ≥ 2`) never prime. "Bounded or eventually prime periods"
   (introduction line 16, title of Corollary 2.4, §5) should read "bounded or eventually prime
   ratios `q_k`".
2. End of §2.1: "every valuation family that is not eventually periodic". This can be read as a
   condition on `(c_k)`, which would exclude period doubling itself (`c_k = 1_{k even}` is
   periodic). Write "with `(c_k)` not eventually constant (equivalently, `m` not periodic)",
   as in Corollary 2.4.
3. Introduction, lines 8–9. "A level set of positive measure can only occur at an odd integer
   level … in two cases" can be misread. Suggested wording: "In two cases we prove that a
   level set of positive measure can only occur at an odd integer level".
4. Introduction, line 11: "finitely many arithmetic progressions" should read "finitely many
   arithmetic progressions plus a finite set", as in Proposition 1.4(a).
5. §1 Remark (2). Theorem 1.2 says nothing about the models of `DENSITY.md` §7 because their
   weights are not integers. Proposition 1.4(b) is not the reason.
6. §1 Remark (1). `DENSITY.md` Remark 7.5 is labelled "Sketch only". Write "sketches" instead
   of "treats".
7. Theorem 2.3 (i). "The second formula below depends on this choice": the first formula
   depends on it too, since it comes from `Σ_{r=1}^{P−1}`.
8. Theorem 2.3 Remark (2). In `toy_toeplitz_results.txt` (1), the closed form is compared with
   partial Fourier sums. Their error is `~1e−4`, and for random patterns it does not decrease
   from `N = 2e5` to `4e5` (`2.13e−4`, `2.12e−4`). The closed form is nevertheless exact: it
   agrees with an independent radial-integral evaluation to `2.3e−15` on both cases of that
   file (`r3_toeplitz_forms.py` (E)). Cite an exact check.
9. Remarks on Corollary 2.4, (3). A smaller instance than `q = 21` exists. At `q = 20` a
   symmetric 0/1 pattern has `V = −2` on 12 of the 19 nonzero residues, so at a single level
   the arithmetic allows `κ = −3` (`r3_galois.py`).
10. The label list at the top of `TOY.md` omits "Corollary". `HYP.md`'s list includes it.
11. HYP §2 item 4, last bullet. For completeness, also list Devin's 5.4(3)–(4) and Devin
    1909.03975, Corollary 2, among the hypotheses that imply TS.

## 6. Scripts (`monotonicity/density2/referee3/`)

All six scripts write their own `*_results.txt`.

**Method shared by `r3_toeplitz_forms.py` and `r3_periodic_flat.py`.** `F'` is evaluated
through the radial-integral form of the resolvent:

    S(t) = Σ_{n≥1} p(n)e^{int}/(M+in) = −i∫_0^∞ e^{iMx} R(e^{it−x}) dx,   R(z) = Σ_{n≥1} p(n)z^n.

This is valid off `(2π/P)ℤ` by Abel summation followed by dominated convergence. It uses none
of the trigonometric manipulations of the proof.

| Script | What it checks | Key output |
|---|---|---|
| `r3_toeplitz_forms.py` | Theorem 2.3 (i): both forms of `w'`, `B(s)` constant in `u`, (R), and the general mass recursion at all `s`. The symmetry-step difference formula. The author's closed form (`toy_toeplitz.py`, imported) | Sanity check against an exact case: `4.4e−16`. Forms of `w'`: `≤ 2.4e−13`. (R): `≤ 2.2e−16`. Mass recursion: `≤ 4.9e−15`. Difference formula: `≤ 3.2e−14`. Smallest `sup_u|w'_{s−1}−w'_s|` over adjacent pairs in non-symmetric cases: `0.162`. Author's closed form: `≤ 2.3e−15` |
| `r3_galois.py` | Step (iv); Corollary 2.4 Steps 1(b) and 2; Remark (3) on Corollary 2.4. Exact arithmetic in `ℤ[ζ_q]` | 124,522 patterns, 0 violations. A non-rational value fills at most 50% of the residues (the proof needs `< 60%`). Smallest `q` allowing `|κ| ≥ 3` at one level: `q = 20` (`κ = −3`). `c_21 = 1` on 12/20 |
| `r3_periodic_flat.py` | Remark (1) on Corollary 2.4: census of flat sub-arcs for periodic symmetric patterns, `2 ≤ Q ≤ 22`. Independent re-check of eight flat cases. Non-symmetric periodic patterns | 11,956 patterns (reproduced). 11,878 with a flat sub-arc. **108,164 triples** (−1: 56,285; +1: 51,879). 9 exact double roots, counted once. Flat cases confirmed to `≤ 1.2e−14`. Non-symmetric: smallest `sup` difference `0.158` over 3,374 pairs |
| `r3_periodic_flat_hpcheck.py` | Validates the root classifier of `r3_periodic_flat.py` against 60-digit solves | 1,600 random pairs and 48 flagged pairs: 0 disagreements. A naive double-precision count gives 108,212 |
| `r3_smallball.py` | Corollary 2.4, as an illustration: small-ball masses with a two-sided bracket from (ii). Proposition 2.1 jumps | Bounded ratios (period doubling; random and near-periodic bounded `q`). Control atom of mass 1/5 detected (ratio `1001` at `η = 10^{−4}`). Jumps `4πM/3` and `±(4πM/3)2^{−a}` confirmed |
| `r3_hyp.py` | HYP Proposition 1.1: sharpness, and the role of `G = H + ⟨Λ⟩`. The lattice facts of HYP §2 | Level set of measure exactly `1/N` (`N = 2, 3, 5, 8`), one point per orbit. `1/2 > 1/4` when `H + ⟨Λ⟩ ≠ G`. Indices `2^{k(T)}` and `1`. Relation-lattice basis spans the kernel (minor `1`), coordinate sums even |

## 7. What I could not check

- **Citations of earlier referees.** I did not read them, by instruction. This covers the
  attributions in Lemma 1.1, Proposition 1.4 and Remarks (2)–(3); referee 2's 12-case and
  `q ≤ 30` checks; and the 108,205 count. I redid the mathematics they support. I cannot
  say how the 108,205 figure was obtained.
- **Quoted theorems.** Cohen's idempotent theorem and Helson's theorem, including their
  numbering in Rudin, and the Riemann–von Mangoldt formula for `L(s,χ)` were used as standard
  inputs and not re-derived.
- **Literature.** I took the statements of Devin, Martin–Ng, Bailleul, Akbary–Ng–Shahabi,
  Hayani and Kaczorowski–Ramaré only from `LITERATURE.md`; I read none of the papers. HYP's
  "Status of the hypotheses" (that none of TS, TS_ℤ, FD is known) is a literature claim I
  could not check.
- **`DENSITY.md` inputs.** Theorem 3.1, Theorem 3.4 and Lemma 1.1 were read but not refereed
  in full. The certified bound `ν_{1/2}({0}) < 0.0891` (Theorem 5.4) was not checked.
- **Flat-piece census.** The census was validated on a sample of 1,648 pairs in 60-digit
  arithmetic, not solved exhaustively at high precision. So 108,164 is a validated count, not
  a proved one.
- **Not in the assigned items.** I did not check Remark (3)'s case of holes at `P_k/2`, nor
  the argument for it that TOY attributes to referee 2. I did not try to settle the open
  question of Remark (3) on Corollary 2.4 (`|κ| ≥ 3` for unbounded composite `q_k`). I did not
  referee §§3–5 of `TOY.md` beyond the "Proved instances" paragraph, nor Proposition 2.1
  beyond its jump values.
