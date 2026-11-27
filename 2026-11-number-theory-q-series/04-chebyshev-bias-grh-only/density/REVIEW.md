# Review history of DENSITY.md

(Section 5 of the working log, published separately; "I" is the session that wrote DENSITY.md.)

**First adversarial read (fresh `opus` subagent, 2026-11-26), of the version without
Propositions 4.4, 7.2, 7.3 and Remarks 2.2a, 7.4.** Verdict: no fatal error. It re-derived
every Theorem/Proposition of §§1–7 present at the time, checked (I5) numerically against
the stored zeros of `χ_{−4}` (agreement `10^{−12}`), tested Theorem 5.2 on 232 cyclic models
and Theorem 5.3 by exact quadrature on 206 toy models and Monte Carlo on 319 adversarial
models with resonant tails (no violation), audited `reflection_certified.py` line by line
and reran it (byte-identical output). Items found and applied:

| item | change |
|---|---|
| 1 (REPAIRABLE) | Theorem 3.5(b) claimed invariance for an arbitrary subset `S` of the level set; false (take `S ∩` a small ball). Restated for the full level set `{F_M = c}`; proof text cleaned. Corollary 3.6 unaffected. |
| 2 | "no atom" was equated with "density exists" in §8 and Remark 2.5; reworded (no atom ⇒ density exists; converse not claimed). |
| 3 | "incomparable with TS" for Devin Thm 5.4(1) was only half shown; Example 3.3a (supplied by the referee) added. |
| 4 | prior art: the level-parameter analogue of Theorem 2.1 is classical (Kaczorowski–Ramaré per Devin; not read); now said in §0 and after Corollary 2.2. |
| 5 | "total mass `< 1`": one-line justification added. |
| 6 | §7 intro and §8 item 5: Theorem 5.3 also uses (NR1)–(NR3); option "arithmetic information" added. |
| 7 | Remark 7.6: "no nonzero tame direction". |
| 8 | Chebyshev's 1853 sum is a mollification "in spirit" only. |
| 9 | §0 row C: hypotheses of Theorem 4.1 stated. |
| 10 | Remark 6.3: mean term `(1−e^{−Mh})²` added. |

**Second adversarial read (another fresh `opus` subagent), of the parts added after the
first read** (Theorem 3.5(b) restated, Example 3.3a, the additions to §2, Proposition 4.4,
Propositions 7.2, 7.3 with `model_certified.py`, the model test). Verdict: no fatal error.
It re-derived Propositions 4.4 and 7.2 (and tested 7.2 numerically at `M = 0.01`), audited
`model_certified.py` line by line, reran it (byte-identical), and checked Proposition 7.3
independently at 50 digits (`max_J|F_N| = 1.24·10^{−15}`, `Σc = 0.27168058…`); for the model
family it verified `σ²_M < 2` near `M = 1/2` (`0.355, 0.543, 0.764` at `M = 0.4, 0.5, 0.6`),
consistency with Theorems 5.2 and 2.1 (`V = 0.586`, bound `0.485` vs atom `0.048`), the
reflection identity (`10^{−14}`), and that on `J` the functions at `M' = 0.4, 0.6` are
`±0.2·G` with `G e^{M't}` constant. Items found and applied:

| item | change |
|---|---|
| 1 (REPAIRABLE) | Over-claim in the model test: the model lacks, besides integrality and the explicit formula, weights `≥ 1`, a non-cyclic hull, and the growth of `N(T)` (unbounded `F`). Remark 7.4, §8 item 5 and §0 row O now say "a proof must use something the model lacks" and list all five. Proposition 7.2 was strengthened to give all weights `≥ 1` (so that feature is not the essential one). |
| 2 (REPAIRABLE) | Example 3.3a needed `ε < 1/6` and disjoint primary/linking indices; rewritten with windows `(log T, T^{1/2−ε}]` and `b_{N(n,1)} ≥ e^{b_n}` (works for every `ε < 1/2`). |
| 3 | `model_certified.py` now asserts `Σc_n < 0.2717`; "weights `≤ 1.10`" labelled [num]. |
| 4–10 | Remark 2.2a(iv) (only resolvents and their products are covered by Prop 4.4); Corollary 3.6 wording; total atom mass over `(0, M_σ]`; the "nothing to transfer" claim extended from `J` to the whole tie set; simplified constants in Prop 7.2; numbering (Prop 4.4 moved to §4.3, remarks of §7 renumbered 7.4–7.7); §8 lists which results need GRH; "continuous version" in Prop 4.4. |

**Not independently refereed after these last edits:** the strengthened form of
Proposition 7.2 (weights `≥ 1`: the extra inequalities are one line each) and the rewritten
Example 3.3a (it follows the referee's own suggestion).

**Third and fourth adversarial reads (fresh `opus` subagents, continued session).**
- Read of §9.1–9.5: no fatal error. Applied: Cor 9.4 restated under NCZ (the central-zero
  case is conditional); Remark 9.3 now uses codimension-one subgroups and is labelled a
  sketch; Thm 9.1 strengthened ("no atoms at all if infinitely many `b_k ≠ 0`"); Helson,
  Cohen and Davenport/Jaffard cited; the exact identity `A_1/A_0 = (1+e^u+e^{2u})/3` in
  Ex 9.2; the table row fixed.
- Read of Thm 9.6: no fatal error. The theorem holds and was re-derived. The construction was
  reimplemented at `M = 0.05, ½, 3` and checked on the Fourier side through `j = 6`.
  Applied: growth sharpened to `(log n)^{5/2}`; the `σ²` addendum (an integer counterexample
  at `M = ½` on `25ℕ` with `σ² < 0.5`); the `M`-dependence in Lemma 9.7(iv) corrected; the
  over-claim at the end of Remark 9.8 removed; the refined question unified as P″, with the
  local counting law; prior art for the technique noted (Menshov–Olevskii correction).
- Not refereed: Remark 9.9 (the Bernstein observation and Question P″), written after the
  fourth read.

**Fifth adversarial read (fresh `opus` subagent): Theorem 9.16 to Remark 9.19.** No fatal error.
- Applied:
  - Pólya–Carlson *does* apply after continuation across one arc. My Remarks 7.5 and 9.16(iii)
    had claimed the opposite; both fixed. Theorem 9.16 is now labelled a folklore-level
    consequence of Pólya–Carlson.
  - The case `I = 𝕋`.
  - The weaker hypothesis `F ∈ c + ℂe^{−Mt}` on `I`.
  - Remark (ii) restricted to bounded weights, with a direct proof.
  - Cor 9.16′ uses a generator `h` of `⟨g_i⟩` (`γ_0/q` need not lie in `G`).
  - Prop 9.17: `θ_0` even, `ℓ_1 ≤ |J|/8`, explicit candidate count.
  - Lemma 9.18 levels restricted to `C ⊂ ℝ∖(1+2ℤ)`, with `σ² ≤ 2−δ` in "P″ = U".
  - Remark 9.19: constant in `Φ`; the later arrows labelled a sketch; the `H^{1/2}` claim
    weakened.
- Referee's scripts added: `arc_rigidity_tests.py` (flat arcs at `F ≡ 1, −1, −3` with
  `c = 1 + 2g(∞)`) and `right_limit_real_weights_check.py` (Prop 9.17: `|F| ≤ 1.5·10^{−15}` on `S`;
  right-limit pattern matched to `2·10^{−5}`).
- Lemma 9.20 was written after this read and is not yet refereed.

**Sixth adversarial read (fresh `opus` subagent): Lemma 9.20.** The lemma is correct. The class
was widened (`φ_γ + φ_{−γ} ∈ ℤ`, the odd part free). Remark (i) was false as stated: Fekete's
theorem needs a conjugation-symmetric set, and certificates for an arc in a half-circle need the
widening. Remark (ii) was empty: Theorem 9.6's set sits in a half-circle arc, where real
coefficients already give the bound. Remark (iii) needed `j ≥ 1` and the `j > J` step; its
caveat was wrong (the threshold is `V/(1+V)`). Remark (iv) needed a fix to its hypothesis, and
the narrow-class duality question is negative. Status: standard and of low use. Rewritten as
Remark 9.20.

**Seventh adversarial read (fresh `opus` subagent): Theorem 9.22 / Corollary 9.23.** Verdict:
correct, with fixable gaps; nothing fatal.
- Re-derived: Lemmas A, B, C, O, F and E; the sigma–delta bookkeeping with spill; positivity;
  integrality; convergence; vanishing; measure.
- Numerics:
  - Lemma B against worst-case errors: ratio `≤ 0.076` of the bound.
  - Toy cascade (`ref922_multiband.py`): multiplicities integral to `3.6·10^{−12}`, `min m = 1`,
    vanishing on `S` to `3.6·10^{−17}`.
  - Starting the offset too early (violating (J3)) gives negative multiplicities, so (J3) is
    needed.
- Applied:
  - Lemma 4.2's `H¹` reason fixed (`|f(n)| ≤ 2M²A(log n)^{θ′}n^{−2}`);
  - `b_0` depends on `M` directly;
  - `θ ≥ 0.92 ⇒ q = 60` flagged as tight; `M_0` labelled numerical;
  - **scope**: P″ fails only with a free counting constant (`A_c ≈ 10^{14}`, larger at small
    `M`), so the fixed-`A_0` form (Conjecture U) and the bounded case are untouched.
  - The text calling P″ open or likely true was updated (P″ paragraph, Remarks 9.11, 9.14, 9.21,
    summary table).
- Not applied (cosmetic or crude, noted only): Lemma F could give `θ′ = 3/(1+β) + ε`; the Lemma O
  test in `thm922_check.py` uses `A = 1`, `q = 20`.
