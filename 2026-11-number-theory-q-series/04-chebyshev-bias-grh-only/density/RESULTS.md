# Results: logarithmic density of Chebyshev's bias under GRH alone

AI-produced (Anthropic's Claude), 2026-11-26. Proofs: `DENSITY.md`. Log: `PROGRESS.md`.
Literature: `LITERATURE.md`.

## Summary (15 lines)

1. **Main target NOT achieved.** It is not proved that under GRH alone the limiting
   distribution of `(log x/√x)(π(x;4,3) − π(x;4,1))` has no atom at `0`; the existence of the
   logarithmic density of `{π(x;4,3) > π(x;4,1)}` under GRH alone stays open. (Literature
   check: still open; it is Problem 13 of arXiv:2407.03530 in the form "ties have density 0".
   Nothing here concerns Problem 12.)
2. **Fallback (2), achieved — improved certified GRH-only bounds.** New LI-free identity
   (time reversal on the hull): `F(θ) + F(θ⁻¹) = 2 + R(θ)`, `|R| ≤ σ²`. Consequences: the
   limiting law is symmetric about `1` up to `σ²` without LI; under GRH(`χ₋₄`) the lower
   logarithmic density of `{π(x;4,3) > π(x;4,1)}` is `> 0.9109` (repository: `0.8842`), under
   GRH(`χ₋₃`) that of `{π(x;3,2) > π(x;3,1)}` is `> 0.9347` (was `0.9112`) [cert-A, Arb]; the
   possible atom at `0`, hence upper minus lower density, is `< 0.0891` resp. `< 0.0653`
   (was `0.1157`). Tier-B values of the repository (`0.9227`, `0.9445`) are not improved.
3. **Fallback (4)-type reduction.** The atom at `0` of the weighted race `D_m` equals the
   jump at `m` of the (monotone) bias as a function of the weight. So under GRH alone the
   density exists for all but countably many weights `m ∈ (−1/2, 1.32]` (`q = 4`); the
   classical race has no atom at `0` iff the bias is continuous in the weight at `m = 0`
   (and then its density exists). No explicit non-exceptional weight is known.
4. **Fallback (1).** Existence (indeed absolute continuity, resp. continuity) under GRH +
   hypotheses TS / FD on the ordinates, each weaker than the hypotheses found in the
   literature (Devin Thm 2.2 and 5.4(1), Bailleul Thm 1.7, Martin–Ng). Contrapositive: an
   atom forces every ordinate to be a rational combination of ordinates above any height.
5. **Fallback (3) NOT achieved**: ties / near-ties of bounded size are not shown to have
   density zero; the obstruction is identified (a local limit statement for primes in short
   intervals).
6. Also proved: every analytic mollification of the race has a density under GRH alone
   (finite-order smoothing creates no new atoms); a level set of the hull function has
   empty interior; and — the obstruction — positive-weight models with an atom at `0`,
   including one **at the classical weight `M = 1/2` with `σ² < 0.55`** (Prop 7.3, certified
   in Arb). Such a model has the hull, the positive weights, the whole family in the weight
   with its identities, and time reversal; so no argument using only these (in particular
   none of the three strategies suggested after note 2, `PROGRESS.md` §4a) can exclude the
   atom. A proof must use some feature the models lack; `DENSITY.md` §7.4 and §8 list five
   (among them integrality of the multiplicities and arithmetic of the dense orbit).

Independent referee by the coordinating session: `REFEREE-2-DENSITY.md` (26 Nov): R and A
**PROOF COMPLETE**, certified numbers reproduced (0.91092486193, 0.93476139844); minor findings
1–2 applied here, cosmetic items listed there.

7. **Integrality and parity (§9, added in the continued session).** The multiplicities of
   zeros are integers, while the bias constant is `½`.
   - In nested model classes this mismatch provably forbids an atom at `0` (Theorem 9.1).
     Real weights in the same class do not forbid it (Example 9.2).
   - **But integrality does not suffice in general.** Theorem 9.6 constructs positive
     integer multiplicities `m_n ≪ (log n)^{5/2}` on the frequencies `n ∈ ℕ` whose model vanishes
     on a set of measure `> 0.1·2π`. So the natural conjecture "integer multiplicities ⇒ no
     atom" (Question P) is false as stated.
   - What the counterexample violates is the counting law `N(T+1) − N(T) ≪ log T`, a
     theorem for the true zeros. Within the construction this is forced by Bernstein's
     inequality. The refined Question P″ (integer multiplicities + counting law + `σ² < 2`)
     is the sharpest open spectral form of the problem: "no" would complete a spectral no-go
     theorem; "yes" on every hull would prove the density under GRH alone.
   - **Arcs versus fat sets (Theorem 9.16, Corollary 9.16′).** In the cyclic model, integer
     multiplicities of any polynomial growth forbid a flat arc at every level outside `1 + 2ℤ`.
     This is a short, folklore-level consequence of Pólya–Carlson once one shows that `g` continues
     across the arc. The same holds in every rank-one hull, with no GRH. Since Theorem 9.6 has an atom on a nowhere-dense set, the open problem is exactly a
     "fat-set Pólya–Szegő theorem": false at growth `(log n)^{5/2}`, open for bounded
     multiplicities. Proposition 9.17 shows that the right-limit (Breuer–Simon) route to it
     must use integrality from the start. Lemma 9.18: in compact model classes, "no atom" is
     equivalent to a uniform small-ball bound.
   - **Question P″ fails (Theorem 9.22, Corollary 9.23).** Noise-shaped rounding (sigma–delta
     with `(1−z)^k`, a cutoff equal to 1 on the target arc, the comb correction of Theorem 9.6
     applied to the *shaped* error) gives positive integer multiplicities `m_n ≪ (log n)^θ`, for
     any `θ > 0`, whose model vanishes on a set of positive measure. For small `M` the example
     also has `σ² < ½` and obeys the local counting law, but with constant `≈ 10^{14}`. So
     integrality, positivity, the counting law with a free constant and `σ² < 2` do not exclude
     the atom. What is left on this line: the counting law with the true constant `A_0`
     (Conjecture U), and multiplicities bounded by a constant, which GRH does not supply.

## What is proved (all under GRH for the character concerned; no LI)

| Result | Where | Status |
|---|---|---|
| atom at `0` = jump of the bias in the weight; countably many exceptional weights; two-sided modulus | `DENSITY.md` Thm 2.1, 2.3, Prop 2.4 | proof; extends `GRH-ONLY.md` Thm 7.1 (which did not give existence) |
| TS ⇒ `ν_M` absolutely continuous; FD ⇒ `ν_M` continuous (Hölder); tame directions; necessary conditions for an atom | Thm 3.1, 3.4, 3.5, Cor 3.2, 3.6 | proof; FD quotes the Turán–Nazarov inequality |
| analytic smoothing ⇒ absolutely continuous law | Thm 4.1 | proof |
| no level set has interior | Prop 4.3 | proof (Weil explicit formula) |
| time reversal, approximate symmetry, `ν_M((−∞,0]) ≤ V_M/(1+(1−σ_M²)²)` (σ_M² < 2) | Thm 5.1, 5.2 | proof |
| refined bound with finite non-resonance | Thm 5.3 | proof |
| certified constants | Thm 5.4, table | `reflection_certified.py` (Arb) |
| band / short-interval reformulations | Prop 6.1, 6.2 | proof |
| finite-order smoothing creates no atoms | Prop 4.4 | proof |
| positive-weight models with an atom at `0`: `σ² > 2` (7.1); `σ²` arbitrarily small, small `M` (7.2); `M = 1/2`, `σ² < 0.55` (7.3) | Prop 7.1, 7.2, 7.3 | proof; 7.3 computer-assisted (`model_certified.py`, Arb) |
| parity rigidity in nested comb classes; real-weight sharpness | Thm 9.1, Ex 9.2 | proof; refereed (no fatal error) |
| positive **integer** multiplicities (`≪ (log n)^{5/2}`; also at `M = ½` on `25ℕ` with `σ² < 0.5`) with an atom at `0` | Thm 9.6 | proof; refereed (no fatal error) |
| reduction to an integer-weight harmonic-analysis question | Cor 9.4 (+ Remark 9.8) | proof (conditional; Question P false as stated, P″ open) |
| cyclic model, integer multiplicities of polynomial growth: no flat arc at levels `∉ 1+2ℤ`; same in rank-one hulls | Thm 9.16, Cor 9.16′ | proof (Pólya–Carlson; folklore-level); refereed |
| bounded real weights: atom on a fat nowhere-dense set whose right limits are carried by it | Prop 9.17 | proof; refereed; numerically checked |
| compact classes: no atom ⇔ uniform small-ball bound | Lemma 9.18 | proof; refereed |
| integer-certificate duality (standard; low use: `φ = 1` gives only `λ(S) ≤ V/(1+V)`) | Remark 9.20 | proof; refereed |
| **positive integer multiplicities `m_n ≪ (log n)^θ` (any `θ > 0`) with an atom at `0`; for small `M` also `σ² < ½` and the counting law, so Question P″ (free constant `A`) has a negative answer. The counting constant is `≈ 10^{14}`, so Conjecture U and the bounded case are untouched** | Thm 9.22, Cor 9.23 (heuristic: Remark 9.21) | proof; refereed (seventh read, no fatal error); toy cascade checked numerically |

## Certified constants (lower bounds for the lower logarithmic density; GRH only)

| race | tier A (in-repo certified zeros) | previous | tier B (simple zeros to `10^8/q`) | previous |
|---|---|---|---|---|
| `π(x;4,3) > π(x;4,1)` | **0.910924861** (Thm 5.3, `N = 20`); 0.901150240 (Thm 5.2) | 0.884276434 | 0.919189115 | **0.922718033** |
| `π(x;3,2) > π(x;3,1)` | **0.934761398**; 0.929813746 | 0.911240824 | 0.941625186 | **0.944555790** |

## Conjectural / not proved

- `ν_{1/2}({0}) = 0` under GRH alone (the main problem).
- Ties and near-ties have density zero under GRH alone.
- That the growth of the zero-counting function and a non-cyclic hull can also be built
  into the obstruction model (expected, not carried out; Remark 7.4).
- Remarks 2.5, 7.5 are sketches.

## Novelty evidence and caveats

- Searches (three subagent passes; search engines and the arXiv PDFs of Devin 1706.06394 and
  1909.03975, Martin–Ng 1710.00088, Bailleul 2007.05763, Fiorilli 1210.6946,
  Akbary–Ng–Shahabi 1306.1657, Hayani 2512.23302, HKMN 2407.03530): the time-reversal
  identity, an LI-free approximate symmetry, a bound "lower density `≥ 1/2` under GRH alone"
  or any explicit GRH-only number for `q = 4` were **not found** there; the nearest prior
  result is Devin's Theorem 2.3 (exact symmetry under a parity hypothesis on the integer
  relations, via a half-turn translation of the torus). The tail criterion (Cor 3.2) was
  not found. Theorem 2.1 extends this repository's own Theorem 7.1; its level-parameter
  analogue is classical. Rubinstein–Sarnak 1994 and the wider literature were not read:
  "not found in the papers listed" is all that is claimed.
- The gain of §3 over the literature is modest (same mechanism as Devin's; independence
  replaced by conditioning), and none of TS, FD is known to hold.
- Rubinstein–Sarnak, Kaczorowski–Ramaré and Nazarov were not re-read; (I5), (I6) are quoted
  from memory of standard sources (see `DENSITY.md` §1, §3.2).

## Files

`DENSITY.md` (theorems and proofs) · `LITERATURE.md` · `PROGRESS.md` ·
`reflection_certified.py` → `reflection_certified_results.txt` [cert] ·
`model_certified.py` → `model_certified_results.txt` (+ `model_M0.5_coeffs.txt`) [cert] ·
`parity_nested_check.py`, `integer_counterexample_demo.py` (+ results) [num] ·
`explore_reflection.py`, `explore_general.py`, `explore_model_arc.py` [num, searches] ·
`sanity_primes.py` → `sanity_primes_results.txt` [num].
