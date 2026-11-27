AI-produced (Anthropic's Claude), adversarial referee report, 2026-11-27.

# Referee report 2: `TOY.md` §0 and §2.2 (Theorem 2.3), and `HYP.md`

**Scope.**
- Read in full: `TOY.md` §0 and §2.2.
- Read because they cite Theorem 2.3: the introduction of `TOY.md` (lines 8–17), the end of
  §2.1 ("What is proved"), and the sentence of §5 that cites the theorem.
- Read in full: `HYP.md`.
- Context: `DENSITY.md` §1 and §3, and `LITERATURE.md`.
- Not read, as instructed: `toy_toeplitz*.py` and their outputs. So the claims of TOY
  Remark (2) about those scripts are not checked here; my own code checks the same identities
  independently.

Scripts and outputs are in `monotonicity/density2/referee2/` (`r2_*.py`, `r2_*_results.txt`),
written from scratch for this report.

**Bottom line.** I found no error in the mathematics of Theorem 2.3 or of `HYP.md` §1. I
re-derived every step by hand, and checked numerically where possible (agreement to `1e−11`
or better). I could not build a counterexample: an exhaustive periodic search gives 108,205
flat pieces, all at levels `±1`.

The problems are statements of **scope**:
- The class of Theorem 2.3 is mis-described ("the regular Toeplitz sequences whose only hole
  is `n = 0`" is false as an identification), and the introduction's "single-hole … any
  patterns" overstates it.
- `TOY.md` §2.1 credits Theorem 2.3 with "no atoms" for the Grandi family. The statement does
  not give this; the proof does (Corollary O1 below).
- `HYP.md` §2.4's "Every hypothesis listed in `LITERATURE.md` implies FD" is false as written.

## Verdicts

| # | Item | Verdict |
|---|---|---|
| T0 | TOY §0 setting | CORRECT |
| T1 | §2.2: definition of `m`; `P_k`-periodicity off `P_kℤ`; valuation family; period doubling | CORRECT |
| T2 | §2.2: "These are the regular Toeplitz sequences whose only hole is `n = 0`" | WRONG as stated (counterexample; harmless for the proof) |
| T3 | §2.2: every periodic sequence is included | CORRECT |
| T4 | Intro: "single-hole Toeplitz 0/1 sequence, with any periods and any patterns" | WRONG as stated (scope overstated) |
| T5 | Proof (i): both forms of `w'`, atom masses, level split, (R), real-analyticity | CORRECT WITH GAP (one line) |
| T6 | Proof (ii): periodicity of `F''`, `L²` bound | CORRECT |
| T7 | Proof (iii): martingale, Markov count, analytic continuation | CORRECT |
| T8 | Proof (iii): symmetry (`A_r`, independence of exponentials, `p_a ≡ 0`, every `v_k` symmetric) | CORRECT |
| T9 | Proof (iii): `B(s) = B_*`, the count `≥ 0.6(q−1)` | CORRECT |
| T10 | Proof (iv): proportion `1/[ℚ(β):ℚ]`, `β_* ∈ ℤ` | CORRECT |
| T11 | Proof (v): bounds on `ρ` and on the modulus of `κ−1−2nρ`, the contradiction | CORRECT |
| T12 | Quantifiers in `i`; circularity; statement = what the proof gives | CORRECT |
| T13 | §2.1: "Theorem 2.3 below shows that the Grandi mechanism does not by itself create atoms" | CORRECT WITH GAP (needs Corollary O1, which the proof gives) |
| H1 | Proposition 1.1 | CORRECT (one garbled sentence) |
| H2 | Lemma 1.2 | CORRECT |
| H3 | TS_ℤ: definition, monotonicity | CORRECT |
| H4 | Theorem 1.3 (both cases) and its consequences | CORRECT |
| H5 | Corollary 1.4 | CORRECT ("strengthens Cor. 3.2" holds only for the atom case) |
| H6 | §2 item 1: TS ⇒ TS_ℤ | CORRECT |
| H7 | §2 item 2: `{2^kα}`: TS fails, TS_ℤ holds | CORRECT |
| H8 | §2 item 3: `{2^kα}` fails FD; Example 3.3a fails TS_ℤ | CORRECT |
| H9 | §2 item 4: "Every hypothesis listed in `LITERATURE.md` implies FD" | WRONG as stated (counterexample); true after scoping |
| H10 | §2 item 4, other bullets | CORRECT |
| H11 | §2 item 5 (proof mechanism) | heuristic; one misleading phrase |

---

## Part A — `TOY.md`

### T0. §0 — CORRECT
- **Dual group.** The dual of `A` is a countable subgroup `G ⊂ ℝ`:
  - `A` is connected, so the dual is torsion-free;
  - `ι` has dense image, so `χ ↦ g` is injective;
  - `A` is metrizable, so the dual is countable.
- **Identities.** `De_g = ige_g`. Coefficientwise, `(D+M)F = M + 2MΣ m(g)e_g = 2Mw`.
  `F ∈ L²` iff `Σm²/(1+g²) < ∞`.
- **True race.** `Re(m_γ r_γ e_γ) = m_γ(2M/(M+iγ))e_γ + m_γ(2M/(M−iγ))e_{−γ}`, as in
  `DENSITY.md` §1.

### T1, T3. The class and its examples — CORRECT
- **Definition.**
  - `k(n)` is finite, and `n/P_{k(n)} ≢ 0 (mod q_{k(n)+1})`. So `v_k(0)` is never used by `m`;
    it only normalises `V_k`.
  - For `P_k ∤ n`, `m_n` depends only on `n mod P_k`.
- **Converse.** Every one-sided sequence that is `P_k`-periodic off `P_kℤ` along a chain
  `P_k | P_{k+1}` arises this way: put `v_{k+1}(r) := m_{P_k r'}` for any `r' ≡ r` with
  `q_{k+1} ∤ r'`.
- **Valuation family.** With `q_k = 2` and `c_k = 1_{k even}` it gives `ε_n` of Proposition 2.1.
- **Periodic sequences** with period `Q ≥ 2`: take `q_1 = Q`, `q_k = 2` and `v_k(1) = m_Q` for
  `k ≥ 2`.
  - Constant sequences need `Q = 2`, since `q_1 = 1` is not allowed.
  - `r2_toeplitz_class.py` (c): 300 random periodic sequences are rebuilt exactly.

### T2. The identification with Toeplitz sequences — WRONG as stated
Quote (lines 249–250): "These are the regular Toeplitz sequences whose only hole is `n = 0` (for
every `k`, `m` is `P_k`-periodic off `P_kℤ`)."

Take the usual meaning: every position `n ≠ 0` is periodic, and the skeleta have density `→ 1`.
Then the class is strictly larger than the parenthetical.

**Counterexample: `m_n := 1` iff `ν_2(n) + ν_3(n)` is odd.**
- Every `n ≠ 0` is periodic, with period `2^{ν_2(n)+1}3^{ν_3(n)+1}`. This holds also for the
  two-sided sequence `m_{|n|}`.
- The `6^k`-skeleton `{ν_2 < k, ν_3 < k}` has density `(1−2^{−k})(1−3^{−k}) → 1`.
- `0` is the only non-periodic position.
- But `m` is `P`-periodic off `Pℤ` for **no** `P ≥ 2`. Periodicity would make `m` constant on
  each class `n + Pℤ` with `P ∤ n`.
  - *`P` has an odd prime factor.* Take the class of `2^j`, `j > a := ν_2(P)`. Its members
    `2^j + kP` with `k` odd have `ν_2 = a`.
    - If `3 | P`, they also have `ν_3 = 0`; choose `j ≢ a (mod 2)`.
    - If `3 ∤ P`, choose `k` odd (Chinese remainder theorem) so that `ν_3` has the wrong parity.
  - *`P = 2^a`.* The class of `3` contains `3 + 2^a k` with `3 ∤ k`, where `m = 0 ≠ 1 = m_3`.
- `r2_toeplitz_class.py` (a) confirms this for every `P ≤ 3000`.

The level-`k` skeleton of this sequence has two hole classes, `2^kℤ` and `3^kℤ`. Its remainder
is the sum of a `2π/2^k`-periodic and a `2π/3^k`-periodic function. So step (iii) cannot
compare adjacent sub-arcs.

The proof uses only the parenthetical property, so nothing in it is affected. What Theorem 2.3
covers is the Toeplitz sequences with **one hole per period at every level, located at the
residue 0** (one-hole patterns with the hole in the last position).

*Side remark.* For non-symmetric patterns only the one-sided sequence `(m_n)_{n≥1}` has this
structure. The toy weights `n ↦ m_{|n|}` are then `P_k`-periodic only on each half-line. The
proof handles this correctly: `p` is the periodic extension of the one-sided sequence, and
`p_a` measures the defect. But the sentence suggests that the toy weights themselves are
one-hole Toeplitz.

### T4. The introduction overstates the scope — WRONG as stated
Quote (lines 13–14): "Theorem 2.3: on the circle, the weights form a single-hole Toeplitz 0/1
sequence, with any periods and any patterns."

"Single-hole" in the usual sense (one hole per period, as in Jacobs–Keane-type constructions)
allows hole classes `h_k + P_kℤ` with `h_k ≢ 0`. These are not covered.

**Example.** `P_1 = 2`, `q_k = 3` for `k ≥ 2`, `h_k = P_k/2 = 3^{k−1}`. Put `m_n = a` for `n` even
and `m_n = c_{ν_3(n)}` for `n` odd, with `c_0 ≠ a`.
- It is `P_k`-periodic off `h_k + P_kℤ` for every `k`.
- It has no hole in `ℤ` at all.
- It is `P`-periodic off `Pℤ` for no `P ≤ 3000` (`r2_toeplitz_class.py` (b)).
- Its level-`k` remainder has frequencies in `h_k + P_kℤ`, so it is **anti**-periodic:
  `F''(t+2π/P) = −F''(t)`.

Remark (3) (lines 331–333) should list both non-covered families: one hole per period not
at `0`, and hole sets made of several classes (e.g. `2^kℤ ∪ 3^kℤ`). See O4 for the
antiperiodic case.

### T5. Step (i) — CORRECT WITH GAP (one line)
**Checked by hand.**
- **First form of `w'`.** From `1/(1−e^{iθ}) = (1+i cot(θ/2))/2`,
  `w' = 1/2 + Σ_r p(r)[cos rt − cot(Pt/2) sin rt] = 1/2 + Σ_r p(r) sin((P/2−r)t)/sin(Pt/2)`.
- **Second form.** With `t = 2π(s+u)/P`, the factors `(−1)^s` coming from `e^{iπ(s+u)}` and from
  `sin(π(s+u))` cancel. The symmetric part of `p` drops out, since the terms `r` and `P−r`
  cancel. This gives the `p_a`-form.
- **Atom masses.** `Σ_k e^{iPkt} = (2π/P)Σ_j δ_{2πj/P}`, with the normalisation
  `δ̂_0 = 1/2π` of §2.1.
- **Level split.** Level `k` contributes `Σ_a v_k(a)e^{2πias/(P/P_{k−1})}·(P/P_k)·1_{(P/P_k)|s}`.
- **Recursion (R).** `R_Mδ_0` jumps by `1`, so `F'` jumps by `2M(2π/P)V = 2εV` at `2πs/P`. This
  gives (R) with exactly the stated signs and factors. `V` is real for symmetric patterns.

**Checked numerically** (`r2_subarcs.py`, 12 cases):
- *Cases:* symmetric patterns, non-symmetric ones, asymmetry at level 1 only and at level 2
  only, and period doubling.
- *Method:* smoothed direct Fourier sums with `2^20` coefficients and a flat-top `C^∞` window,
  evaluated at `u ∈ [1/4, 3/4]`.
- *Agreement with step (i):*
  - both forms of `w'`, to `≤ 1.7e−11`;
  - `w' = 1/2` on open sub-arcs in the symmetric case, to `≤ 7.3e−12`;
  - `F' = 1 + 2B(s)e^{−εu}`, to `≤ 8.6e−14`. Here `B` is the periodic solution of the **full**
    recursion, with atom masses from a direct DFT of `p`, not from `V`;
  - `C(s) = V_{i+1}(s mod q)` for `q ∤ s`, to `≤ 2.4e−13`;
  - (R) itself, to `≤ 6.5e−14`.
- *Non-symmetric cases:* `F'` is visibly not of the form `1 + 2Be^{−εu}`.
- *Also:* `p_a ≡ 0` iff `v_1, …, v_{i+1}` are symmetric, in all cases.

**The gap.** "On each open sub-arc `w' = …`" identifies a distribution with a pointwise sum.
One line should say why:
- the Abel means converge in `D'(𝕋)`, because the coefficients are bounded;
- off `(2π/P)ℤ` they converge locally uniformly, because `Σ_{n≥1}p(n)z^n = R(z)/(1−z^P)` is
  rational with poles only at `P`-th roots of unity;
- `F'` is then the real-analytic solution of `(D+M)F' = 2Mw'` on each open sub-arc.

In the non-symmetric case `F'` also has logarithmic singularities at the sub-arc endpoints.
This is harmless, but should be said, since only `F' ∈ L²` is known a priori.

**Clarity.** The second form of `w'` depends on the representatives `r`: under `r ↦ r−P`,
`e^{−2πiru/P}` gains a factor `e^{2πiu}`. It is correct for `r = 1, …, P−1`, the range of the
generating function. State this range on line 270, and again in the symmetry step (line 297),
which also uses it (distinct `ω_r`, pairing `r ↔ P−r`).

### T6. Step (ii) — CORRECT
- **Periodicity.** The frequencies of `F''` lie in `Pℤ`, so `Φ` is the same `L²` function for
  every `s`.
- **Bound.** `∫_0^1|Φ|² = Σ_{n∈Pℤ∖0} m_{|n|}·4M²/(M²+n²) ≤ 2Σ_{k≥1} 4M²/(P²k²) =
  (4π²/3)M²/P²`, using `m² = m`.
- **Numerically**, computed as `F − F'` from two independent direct sums:
  - it is `2π/P`-periodic to `≤ 1.1e−14`;
  - its grid `L²` norm equals the coefficient sum, and is below the bound.

### T7–T9. Step (iii) — CORRECT
- **Martingale.** The level-`i` arcs refine and generate the Borel σ-algebra, so
  `λ(E∩I_i(t))/λ(I_i(t)) → 1_E(t)` a.e.
- **Markov count.** With `b` bad sub-arcs, `0.1b < 0.01q`, so `b < 0.1q` and `b ≤ ⌊0.1q⌋`. In
  particular `b = 0` when `q < 10`.
- **Analytic continuation.** `U_s ∩ U_{s'}` has measure `≥ 0.8`. Two functions real-analytic on
  `(0,1)` that agree on a set of positive measure agree everywhere. Applying
  `(P/2π)d/du + M` passes the agreement to `w'`.
- **Symmetry.**
  - I recomputed `w'_{s−1} − w'_s = −Im[Σ_r A_r e^{iω_r u}]/sin(πu)`. The text writes "reads
    `Im Σ_r A_r e^{iω_r u} = 0`", which is the same; the sign is immaterial.
  - With `ω_{P−r} = −ω_r`, writing `Im Z = (Z−Z̄)/2i` gives `Σ(A_r − Ā_{P−r})e^{iω_r u}`.
  - Oddness and reality of `p_a` give `Ā_{P−r} = −A_r`. Since `1 − e^{2πir/P} ≠ 0`, `p_a ≡ 0`.
    One pair of adjacent good sub-arcs suffices.
  - *Numerically:* the difference formula holds to `2.5e−11`.
  - *Numerically:* whenever `p_a ≢ 0`, the difference is nonzero on **every** adjacent pair
    (minimum `1.3e−2`). This includes a case with `v_1` symmetric and only `v_2` not.
- **Counts.** `2⌊0.1q⌋ < q−1` and `(q−1) − 2⌊0.1q⌋ ≥ 0.6(q−1)` hold for every `q ≥ 2`. I also
  checked them for all `q ≤ 10^6` (`r2_misc.py`). The pairs lie inside `I`, so `q ∤ s` and
  (R) applies.

### T10. Step (iv) — CORRECT
- **Galois action.** `σ_a(V(e)) = V(ae)`, because the `v(r)` are rational.
- **Proportion.**
  - `{a : V(ae_0) = β}` is a coset of `Gal(ℚ(ζ_q)/ℚ(β))`.
  - The map `a ↦ ae_0` has equal fibres, and that set is a union of fibres.
  - Together these give the proportion `1/[ℚ(β):ℚ]`.

**Numerically** (`r2_galois.py`, exact arithmetic in `ℤ[x]/(Φ_q)`):
- *Range:* all patterns for `q ≤ 13`; 300 random patterns (half symmetric) for `14 ≤ q ≤ 30`.
- *Galois action:* 2,274,586 checks of `V(ae) = σ_a V(e)`, with 0 mismatches.
- *Proportion law:* 0 mismatches over all orbits, values and patterns. The degree is computed
  **without** the Galois action, as the number of distinct eigenvalues of multiplication
  by `β`. 182 exact sympy minimal polynomials agree with it.
- *Majority values:* every value taken on `≥ 60%` of `{1,…,q−1}` is a rational integer.

**Observation, not an error.** The majority value need not be `0` or `−1`.
- For `q = 21` and `v = 1_{(ℤ/21)^×}`, `V(e) = c_{21}(e) = μ(21) = 1` on the 12 units, i.e. on
  60% of `{1,…,20}`.
- For `q = 13·17` the units are 87% of the residues, above the 80% needed.

So for unbounded composite `q_k` the method allows `n = 1`, i.e. `κ = 3`. Theorem 2.3 claims
only "odd", so this is consistent; see O3.

### T11. Step (v) — CORRECT
- **`ρ` bounds.** `ρ(0) = ε/(1−e^{−ε}) ≤ 1+ε/2+ε²/12 ≤ 1+ε` and `ρ(1) = ε/(e^ε−1) ≥ 1−ε/2`. On
  a grid, `(ρ−1)/ε ∈ [−0.5, 0.582]` for `ε ≤ 1`.
- **Lower bound.** Both cases of the bound on `|κ−1−2nρ|` are right. `ε_0(κ)` exists; grid
  values are in `r2_misc.py` (3).
- **Contradiction.** It follows because `P_i → ∞`.

### T12. Quantifiers, circularity, statement — CORRECT
The order of the argument is sound:
1. One density point `t` is fixed.
2. Symmetry of `v_1, …, v_{i+1}` is proved for every large `i`. This uses only
   real-analyticity, in both cases, not (R).
3. (R) is then available at every level.
4. Steps (iv)–(v) are run at one large `i`.

The proof gives exactly the statement: every `v_k` is symmetric, and `κ` is odd.

**Counterexample search: none found.**
- **Periodic, symmetric** (`r2_periodic_flat.py`; my own derivation).
  - Explicit form on sub-arcs: `F = (1−2m_Q) + [2ε/(1−e^{−εQ})]e^{−εu}g_s(e^{−ε})`.
  - All 11,956 primitive nonzero patterns with `Q ≤ 22`, all `M > 0`: 108,205
    (pattern, sub-arc, `M`) flat pieces, in 11,880 patterns.
  - They all lie at levels **`+1` and `−1` only**, namely `1 − 2m_Q`.
  - Ten random hits, re-checked by direct Fourier sums, are flat to `≤ 3e−14`.
- **Periodic, non-symmetric.** For every pattern with `Q ≤ 12`, `w` is non-constant on every
  sub-arc. `F` is flat on a sub-arc iff `w` is constant there, whatever `M`. So there is no
  flat piece for any `M`.
- **Period doubling.**
  - A crude smoothed-sum small-ball test (`r2_misc.py` (4)) suggested growth at `κ = ±1` for
    `M = 2`.
  - A deep-level evaluation resolves it (`r2_pd_smallball.py`: level `2^{22}`, sup error
    `< 1e−6`).
  - Near `−1`: no mass within `3e−4`. Near `+1`: the ratio stays bounded (`≈ 0.8` at
    `δ = 1e−4`).
  - The law is lumpy, not atomic, consistent with O1.

### T13. §2.1, "What is proved" — CORRECT WITH GAP
Quote (lines 238–239): "Theorem 2.3 below shows that the Grandi mechanism does not by itself
create atoms."

As stated, Theorem 2.3 allows atoms at odd levels, and periodic weights do have atoms at `±1`
(T12). So the sentence does not follow from the statement. It does follow from the **proof**:
see O1, with `q_k = 2`. Either add that corollary, or weaken the sentence to "creates no atom
at a level that is not an odd integer".

Related: §5 (lines 429–430) says "Theorem 2.3 (symmetric single-hole Toeplitz weights, any
periods)". This understates the theorem: non-symmetric patterns are covered too, and they have
no level set of positive measure at all.

## Part B — `HYP.md`

### H1. Proposition 1.1 — CORRECT
- **Finite `K`.** `K = H^⊥` is isomorphic to the dual of `G/H`, by annihilator duality for
  the discrete group `G` (via (I2)). So `|K| = N`.
- **Separation.** `G = H + ⟨Λ⟩` implies that the `e_λ` (`λ ∈ Λ`) separate the points of `K`:
  if they agree at `k` and `k'`, so does every `e_g`.
- **Invariant version.** `F̃_H := N^{−1}Σ_{k∈K}F_H(·k)` is `K`-invariant and equals `F_H` a.e.
  (compare Fourier coefficients). The proof should say this.
- **Counting.** The counting identity holds.
- **`D_{k,k'}`.** It is a nonzero trigonometric polynomial in the pairwise distinct characters
  `e_{±γ}`, `γ ∈ Λ ⊂ (0,∞)`.
- **`Λ = ∅`.** Then `N = 1`, and the bound is trivial.

**Wording slip.** Line 26 says "The characters `g ↦ e_g(k)` of `G/H`, `k ∈ K`, separate the
points of `K`". Characters of `G/H` are functions on `G/H`. The intended statement: the
functions `k ↦ e_g(k)` (`g ∈ G`; they depend only on `g + H`) separate the points of `K`.

**Numerical check** (`r2_hyp_orbits.py`), on two finite models:
- *Model A:* on `𝕋`, frequencies `2^k`, `N = 2, 4, 8, 16`.
- *Model B:* on `𝕋²`, with `G/H ≅ ℤ/8 × ℤ/2` and `N = 16`, computed by Smith normal form.

Results:
- `F_H` is `K`-invariant to `2.5e−14`.
- `F_Λ` separates every sampled `K`-orbit (minimum gap `> 0` over `10^5`–`2·10^5` samples).
- `E|S_δ| = Nλ(E_δ)` within Monte Carlo error.
- In Model B, `P(|S_δ| ≥ 2)/P(|S_δ| ≥ 1)` is `7e−2`, `4e−3` and `0` for `δ = 3e−2, 3e−3,
  3e−4`, as `|S| ≤ 1` a.e. predicts.

### H2–H4. Lemma 1.2, TS_ℤ, Theorem 1.3 — CORRECT
- **Lemma 1.2.**
  - For `T ≥ max Λ`, `Γ_{>T} ⊂ Γ∖Λ`.
  - If TS fails, every `λ ∈ Λ` has a positive multiple in `⟨Γ∖Λ⟩`.
  - A finitely generated torsion abelian group is finite.
- **TS_ℤ.** The index is nondecreasing in `T`, because the tails shrink.
- **Theorem 1.3, first case.** Lemma 1.2 with `Λ = Γ_T` (finite) gives TS. Then `DENSITY.md`
  Theorem 3.1 applies.
- **Theorem 1.3, second case.** `Γ∖H_T ⊂ Γ_T`, and Proposition 1.1 applies for every `T`.
- **Consequences.** They follow from (1.1) as in `DENSITY.md` 2.1(d) and 3.4. NCZ is
  automatic for `q = 4`.

### H5. Corollary 1.4 — CORRECT
- **Eventual constancy.** A bounded nondecreasing integer function of `T` is eventually
  constant. Nested subgroups of equal finite index coincide.
- **Multiples `n ≤ 1/α`.** The order of `γ mod H_∞` divides `|G/H_∞| ≤ 1/α`, and
  `H_∞ ⊂ ⟨Γ_{>T}⟩` for every `T`.

Minor (line 71): "strengthens `DENSITY.md` Corollary 3.2" is true only for the **atom** case.
Corollary 3.2 also covers laws that are not absolutely continuous but have no atoms.

### H6–H8. §2, items 1–3 — CORRECT
- **Item 1.** A nonzero ℚ-linear form vanishing on the tail produces an element of infinite
  order in `G/⟨Γ_{>T}⟩`, so the index is infinite.
- **`{2^kα}`.**
  - Index: `2^{k(T)}`, computed exactly by gcd; it is `2, 16, …, 2^{20}` for `T = 1, …, 10^6`.
  - TS fails.
  - `G_T = ℤα`, so every nonzero homomorphism `G_T → ℤ` is nonzero on all `≈ log₂T`
    elements: FD fails.
- **Example 3.3a.**
  - `b_n = (b_n + b_{N(n,j)}) − b_{N(n,j)}`, with both terms in `X`.
  - `b_{N(n,j)}` is in `X`, because `{b_n : n ≥ 1}` includes the linking indices.
  - So the index is `1`. A truncated version was checked by Smith normal form.

### H9. §2 item 4: "Every hypothesis listed in `LITERATURE.md` implies FD" — WRONG as stated
Quote (line 105): "Every hypothesis listed in `LITERATURE.md` implies FD, and so implies the
disjunction."

`LITERATURE.md` (third check) also lists Devin's Theorem 2.3 hypothesis,
`Σk_γγ = 0 ⇒ Σk_γ ≡ 0 (mod 2)`.

**Counterexample: `X = {(2k+1)α : k ≥ 0}`.** Checked in `r2_hyp_orbits.py` (5).
- **Parity holds.** The relation lattice has the basis `e_k − x_k e_0`, and each coordinate
  sum `1 − x_k` is even.
- **FD fails**, since `G_T = ℤα`.
- **TS_ℤ fails.** The index is `1`, since two consecutive odd numbers generate `ℤ`.
- Hence `X` fails the disjunction "TS_ℤ or FD" as well.

`LITERATURE.md` also lists results that assume GRH only (Hayani; Kaczorowski–Ramaré, which was
not read). GRH is not known to imply FD.

**Fix.** Restrict the claim to the hypotheses that `LITERATURE.md` lists for continuity or for
existence of the density:
- Devin, Theorems 2.2 and 5.4(1)–(4);
- Devin 1909.03975, Corollary 2;
- Bailleul, Theorem 1.7;
- Martin–Ng, Theorems 1.5(a) and 1.7(a);
- LI.

Also add that the hypotheses stated for several `L`-functions (Martin–Ng, Bailleul,
Devin 1909.03975) must be read for the union of the ordinates. FD, TS and TS_ℤ are defined in
`DENSITY.md` for a single real `χ`.

### H10. §2 item 4, other bullets — CORRECT
- "Strictly weaker than FD": shown by `{2^kα}`.
- The hypotheses that imply TS also imply TS_ℤ, since TS ⇒ TS_ℤ.

### H11. §2 item 5 — heuristic
Quote (lines 112–113): "identical remainders on the translates force the finite part to
separate them."

The two mechanisms differ:
- In Proposition 1.1 the separation is not forced by the remainders. It is proved from
  `G = H + ⟨Λ⟩`, and it is what bounds the level set to one point per orbit.
- In `TOY.md` step (iii), the invariant remainder forces the non-invariant parts to
  **coincide** on good translates.

Suggested wording: "in both, the remainder is invariant under a finite group of translations,
so on each orbit the level set is governed by the non-invariant part."

**Status and novelty.**
- `HYP.md` claims no novelty with respect to the literature.
- Its comparisons concern abstract frequency sets, and are correctly hedged ("no claim that
  ordinates look like this").
- `LITERATURE.md` searched only for ℚ-span conditions. If TS_ℤ is to be presented as new, a
  search for finite-index ("subgroup generated by the large ordinates") conditions is still
  needed.

## Required fixes
1. **`TOY.md` lines 249–250.** Replace "These are the regular Toeplitz sequences whose only hole
   is `n = 0`" by the definition actually used. For example: "Equivalently, for some chain
   `P_k | P_{k+1}`, `m` is `P_k`-periodic off `P_kℤ` for every `k` (one hole per period at each
   level, located at `0`). This is a proper subclass of the regular Toeplitz sequences whose
   only non-periodic position is `0`; e.g. `m_n = [ν_2(n)+ν_3(n) odd]` is not in it."
2. **`TOY.md` lines 13–14.** Replace "single-hole Toeplitz 0/1 sequence" by "Toeplitz 0/1
   sequence with one hole per period, located at `0`".
   - In Remark (3), add: "one hole per period not at `0` (e.g. `h_k = P_k/2`, whose remainder
     is antiperiodic)".
   - Also add: "holes in several classes (e.g. `2^kℤ ∪ 3^kℤ`)".
3. **`TOY.md` lines 238–239.** Either add Corollary O1, which gives "no atoms" for the Grandi
   family, or weaken the sentence to "creates no atom at a level that is not an odd integer".
4. **`HYP.md` line 105.** Scope the claim as in H9, and name Devin's parity hypothesis as an
   exception, with the example `{(2k+1)α}`.
5. **`TOY.md` step (i).** Add the one-line justification of T5, and the range `r = 1, …, P−1`
   (lines 270 and 297).
6. **`HYP.md` line 26.** Fix the wording of the separation bullet, and write out the `K`-average
   version of `F_H`.

## Optional improvements

**O1. A stronger corollary, from the given proof.** Suppose `sup_k q_k < ∞`, or `q_k` is prime
for all large `k`. If `λ{F_M = κ} > 0`, then `m` is periodic and `κ = 1 − 2m_P ∈ {±1}`, where
`P` is a period. In particular, the period-doubling weights and every non-periodic valuation
family give an atomless `F_M`, for every `M`.

*Proof sketch* (bounded case). Let `Q := sup q_k`, and use density `→ 1` instead of `0.99`.
1. **All sub-arcs are good.** For large `i`,
   `λ(J_s∖E) ≤ λ(I∖E) ≤ λ(I)/(4Q) ≤ λ(J_s)/4` for **every** sub-arc of `I`.
2. **`V` is constant.** So `B(s) = B_*` on all `q` sub-arcs, and (R) gives `V_{i+1}(e) = β_*`
   for **every** `e ≠ 0`. In the prime case, step (iv) plus transitivity of `(ℤ/q)^×` on
   `{1, …, q−1}` gives the same.
3. **Constant pattern.** Invert the DFT, using `Σ_e V(e) = q v(0) = 0`: `v_{i+1}(r) = −β_*`
   for all `r ≠ 0`. So `β_* ∈ {0, −1}`.
4. **`β_*` is the same at all large levels.** On `U_s`, `|Φ| ≥ |κ−1−2β_*| − 2ε`, and
   `∫|Φ|² → 0`. So `κ = 1 + 2β_*`, and `β_* = (κ−1)/2` for all large `i`.
5. **Periodicity.** For `i ≥ i_1`, `v_{i+1} = ((1−κ)/2)·1_{≠0}`. So `m = (1−κ)/2` on
   `P_{i_1+1}ℤ∖0`, and `m` is `P_{i_1+1}`-periodic. ∎

Without either assumption, the same argument shows that, along a density point, the
proportion of `e` with `V_{i+1}(e) = (κ−1)/2` tends to `1`.

**O2. Integer patterns.** The proof works verbatim for integer-valued patterns with
`Σm²/(1+n²) < ∞`. Only the constant in (ii) changes; alternatively, use that the tail tends
to `0`. The authors' notes already say this; it belongs in the statement.

**O3. Open question.** Can `|κ| ≥ 3` occur when the `q_k` are unbounded and composite? The
Galois step does not exclude it (T10). By O1 it cannot occur when the `q_k` are bounded.

**O4. The antiperiodic case** (hole at `P_k/2`) looks amenable.
- On adjacent good sub-arcs, `F'_s + F'_{s+1} = 2κ`.
- In the symmetric case the a.c. part of `w'` is the constant `1/2 − m_P`.
- That forces `κ = 1 − 2m_P`, which is odd.

I checked this only as a sketch.

**O5.** In Remark (2), mention that the identities were reproduced independently (this
report).

## Scripts (`monotonicity/density2/referee2/`)

| Script | What it checks | Key output |
|---|---|---|
| `r2_subarcs.py` | step (i)/(ii) against direct smoothed Fourier sums: both forms of `w'`, the exponential form, (R), atom masses, the symmetry-lemma identity, remainder periodicity, the `L²` bound | all agreements `≤ 2.5e−11`; 12 cases, 0 failures |
| `r2_periodic_flat.py` | exhaustive flat-piece search, periodic weights | 108,205 flat pieces, levels `{−1, +1}` only; none for non-symmetric patterns |
| `r2_galois.py` | step (iv) for `q ≤ 30`, exact arithmetic | 0 failures; a majority value `+1` occurs (`q = 21`) |
| `r2_toeplitz_class.py` | the setting claims | the identification is false (two examples); periodic inclusion holds |
| `r2_misc.py` | counting inequalities, `ρ` bounds, `ε_0(κ)`, a small-ball exploration | all true; the exploration is superseded by the next script |
| `r2_pd_smallball.py` | period doubling: the law near `±1`, at a deep level | no atom (consistent with O1) |
| `r2_hyp_orbits.py` | HYP Proposition 1.1 on finite models; the examples of HYP §2 | consistent; the parity hypothesis implies neither FD nor TS_ℤ |
