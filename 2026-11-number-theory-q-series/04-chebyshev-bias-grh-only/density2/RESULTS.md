# Results of the second attack (density2): no atom at 0 under GRH alone

AI-produced (Anthropic's Claude), 2026-11-27. Write-ups: `TOY.md`, `HYP.md`, `INVLO.md`,
`LINE1.md`. Log: `PROGRESS.md`. Referee reports: `REFEREE-1-TOY.md`, `REFEREE-2.md` (and
`REFEREE-3.md`). Nothing here concerns Problem 12 of arXiv:2407.03530; the
relevant item is Problem 13 (ties). No result is verified in Lean.

## Summary

1. **Main statement NOT proved.** Under GRH(`χ_{−4}`) alone it is still open whether the
   limiting law of `(log x/√x)(π(x;4,3) − π(x;4,1))` has an atom at `0`. The same is open for
   `χ_{−3}` and the weighted races `D_m`.
2. **New sufficient hypotheses** (fallback (1)-type: existence of the densities under GRH +
   a weakening of LI):
   - `HYP.md` Theorem 1.3: GRH + NCZ + **TS_ℤ** (the large ordinates generate subgroups of
     unbounded index).
   - `INVLO.md` Theorem C: GRH + NCZ + **graded LI** GLI(ψ) (every integer relation among the
     ordinates has length `≥ ψ(largest ordinate involved)`, with `ψ(t) log t/t → ∞`).

   TS_ℤ is strictly weaker than TS and incomparable with FD. GLI is incomparable with TS, FD
   and TS_ℤ. Together with `DENSITY.md` (TS, FD), the disjunction "TS_ℤ or FD or GLI" suffices.
3. **What an atom forces** (structure, all GRH + NCZ):
   - `HYP.md` Corollary 1.4: an atom of mass `α` makes `⟨Γ_{>T}⟩` eventually constant, of
     index `≤ 1/α`. Each ordinate has a multiple `nγ`, `n ≤ 1/α`, that is an integer
     combination of ordinates above any height.
   - `INVLO.md` Corollary C′ and Proposition E: infinitely many short relations (length
     `≤ Ct/log t`). More strongly, for every large `L`, almost all ordinates of every dyadic
     block up to height `≍ α²L log L` occur in integer relations of length `≤ L`.
   - `LINE1.md` Proposition 4.1: an atom bounds the size of free sets of large-coefficient
     ordinates by `C/α²`.

   None of these is known to be contradictory, even together. They are statements about
   exact integer relations, to which the statistics of the zeros are blind (`LINE1.md` §2).
   The explicit formula excludes only cyclic spectra and slowly growing rank-one
   denominators (`LINE1.md` §1).
4. **The toy question** (integer weights on the circle; line (2), handed to the predecessor
   session):
   - `TOY.md` Theorem 1.2 + Proposition 1.4: if `w` is a finite measure (weights on finitely
     many APs), a positive-measure level set occurs only at odd integers. This never applies
     to the true race.
   - `TOY.md` Theorem 2.3: the same parity law for Toeplitz 0/1 (or integer) weights with one
     hole per period at `0`, any periods and patterns. Non-symmetric patterns admit no
     positive-measure level set at all.
   - `TOY.md` Corollary 2.4: bounded or eventually prime period ratios `q_k`, or valuation families: an
     atom forces periodic weights and level `±1`. So period doubling is atomless.
5. **Lines assessed and killed**:
   - line (1), transport of the `ℚ`-independence of `{log p}` through the explicit formula:
     killed in naive form (`LINE1.md` §§1–3). The structure/arithmetic combination of §4
     leaves an exact gap: return time `≤ (2+o(1))log U` is needed, `U^{1−o(1)}` available;
   - line (3), the prime side: killed (correlations of `χ(p)` beyond GRH needed);
   - Conjecture U by inverse Littlewood–Offord: done as far as it goes (`INVLO.md` §4). It
     proves the unsaturated half. The saturated half contains the open cyclic toy question.
6. **Prime-side hypotheses** (`PRIMECORR.md`, GRH + NCZ):
   - Theorem P1 (exact duality): log-averaged moments of windowed prime sums at relative
     width `1/T` equal the independent-phase prediction plus a nonnegative sum over integer
     relations among the zeros.
   - Theorem P2: if these moments match up to order `L` at one scale `T ≫ L^{3/2}log²L`, every
     atom is `≤ C/√L`. The infinite version is LI with simple zeros. So prime-side
     correlation hypotheses of this kind are zero-side ones in disguise.
   - Remark P4 (heuristic): pointwise Hardy–Littlewood asymptotics cannot supply them.
   - Theorem P7 (band decoupling): a prime-side hypothesis strictly weaker than LI with simple
     zeros. Suppose mixed moments of two windowed sums, with disjoint frequency bands, show a
     free **nonempty** band at height `s` with `s/log s = o(L)`, for infinitely many `L`. Then
     there are no atoms. For spectra with simple zeros it contains GLI (Theorem C). Corrected
     after referee 4b: the first version held vacuously for empty bands.
   - Proposition P5 / Corollary P6: LI spectra are dense. Hence, given the predecessor's draft
     Theorem 9.22, the spectral data cannot decide the question. No stable property shared by
     its cyclic models can imply "no atom". A proof must use a property of the true zeros that
     those models lack.

## Status of the reads

| Result | File | Reads |
|---|---|---|
| Lemma 1.1, Theorem 1.2, Lemma 1.3, Proposition 1.4 | `TOY.md` §1 | **two reads**: referees 1 and 3 (fixes applied; no mathematical error) |
| Propositions 2.1, 3.1, 3.2 | `TOY.md` §§2–3 | referee 1 (fixes applied) |
| Theorem 2.3 | `TOY.md` §2.2 | **two reads**: referees 2 and 3 (no mathematical error; scope fixes applied) |
| Corollary 2.4 | `TOY.md` §2.2 | suggested by referee 2 (O1), re-derived; read by referee 3 (correct; Remark (1) corrected) |
| Proposition 1.1, Lemma 1.2, Theorem 1.3, Corollary 1.4 | `HYP.md` | **two reads**: referees 2 and 3 (no mathematical error; fixes applied) |
| Proposition 4.1; Lemma 4.2 (proof sketched only) | `LINE1.md` §4 | referee 4 (Prop. 4.1 correct; Lemma 4.2 sketch honest; fixes applied) |
| Lemmas A, D; Propositions B′, B, E; Theorem C; Corollary C′ | `INVLO.md` | referee 4 (all labelled proofs correct; summaries and §4 fixed) |
| Lemma P0; Theorems P1, P2; Corollary P3; Proposition P5; Corollary P6 (conditional) | `PRIMECORR.md` | referee 4 (correct; two small gaps filled; barrier conclusion narrowed) |
| Lemmas P0′, P1′; Theorem P7 (band decoupling) | `PRIMECORR.md` §3 | referee 4b (P7 wrong as stated for empty bands; corrected with condition (N); GLI comparison restricted to simple zeros; P0′/P1′ correct with gaps filled) |

## Disclosure

All write-ups, proofs, scripts and referee reports in `monotonicity/density2/` were produced
by an AI system (Anthropic's Claude), with independent AI referee passes as listed. They have
not been checked by a human mathematician.
