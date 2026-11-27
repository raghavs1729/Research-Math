# Envelope certificate: summary (2026-11-24)

## What is certified

Under GRH(χ), the paper's variance identity and its `M`-derivative, the
certified zero lists, FLINT/Arb, and the pencil proofs of Lemma D / Proposition
C₀ / Corollary C′ (`../BRIDGE-ATTEMPT.md`; hypotheses (H1)–(H5) of `METHOD.md`):

| modulus | statement | range of `M = m + 1/2` | cells | not certified | worst margin |
|---|---|---|---|---|---|
| `q = 4` | `δ′(M) < 0`, with `−Mδ′(M) ≥ LB > 0` per cell | `[1, 2000]` | 765 (ratio 1.01) | 0 | `NUM/(MAIN−TAILF) ≤ 0.44115` (top cell); `LB ≥ 1.3188e−3` |
| `q = 3` | same | `[1, 2000]` | 765 (ratio 1.01) | 0 | `≤ 0.43819` (top cell); `LB ≥ 1.3532e−3` |

* The decision per cell is `lower endpoint of (MAIN − TAILF − NUM) > 0` in Arb
  (128 bits). `K = NUM/(MAIN − TAILF)` is a margin indicator: on the cells
  covering `[3.6, 20]` (`q = 4`) it is at most `0.3614`, on those covering
  `[4.4, 20]` (`q = 3`) at most `0.3364`, on those covering `M ≤ 400` at most
  `0.4219` (`q = 4`) resp. `0.4169` (`q = 3`) (these maxima include the cell that
  straddles the right endpoint), and it creeps towards `1/2` as `M → ∞`. The
  criterion never has less than a factor 2.26 to spare.
* Sharpness (not soundness): on the range common with the direct certificate
  the bound loses 4–17 % of `M|δ′|`; near `M = 1` it loses more, up to 38 %
  (`q = 4`) and 49 % (`q = 3`) at `M = 1`, against a floating-point evaluation
  of `−Mδ′` (`../REFEREE-5.md` §3.3). The bound stays positive there by a wide
  margin (`K ≤ 0.08`).
* Zero data used: the 12 lowest ordinates (as balls `± 10⁻²⁰`) and the
  completeness of the lists; plus, for the 127 cells with `M_a < 3.5`, the full
  list for an enclosure of `S_4`.
* Runtime: 75 s (`q = 4`) + 76 s (`q = 3`) wall in the rerun of 2026-11-24 (118 s
  + 103 s in the first run, on a more loaded machine), one niced process; about
  0.1 s per cell.
* Required ranges `[3.6, 400]` (`q = 4`) and `[4.4, 400]` (`q = 3`) are covered,
  extended down to `1` and up to `2000`.

## Second, independent certificate on `[3.6, 320]` / `[4.4, 350]`

`monotonicity/cert/` certifies `lo ≤ δ′ ≤ hi < 0` there by enclosing the
oscillatory integral. `consistency.py` compares all intersecting cell pairs
(1061 for `q = 4`, 1028 for `q = 3`): the necessary relation `LB ≤ x·|lo|` holds
in every pair (largest ratio 0.945); the envelope bound captures between 83 %
and (up to enclosure width) 100 % of `M|δ′|`. Sample lines: `consistency_q4.txt`,
`consistency_q3.txt`. The certificates share the hypotheses (H1)–(H4), the zero
files and the module `../cert/certify_derivative.py` (zero loading, evaluation of
`(log Λ)′` and `(log Λ)″`, `M`-grid, printing); they share no integrand, no
quadrature and no tail model. The comparison validates `MAIN − NUM` only: on the
common range `TAILF` is below the width of the direct enclosures, so it would
not notice an error in `TAILF` (Referee 5's mutant "`TAILF` dropped" passes it);
`TAILF` is covered by the referee's independent evaluation of the cell quantity
at 26 cells (`LB/E ∈ [0.983, 0.997]`).

## Zero-data-free analytic leg (Theorem E, `ANALYTIC-LEG.md`)

For **every** real primitive `χ` (GRH, `L(1/2,χ) ≠ 0`) and every `M ≥ 40`:
`−Mδ′(M) ≥ G(σ²)/(sqrt(2π)σ) ≥ 0.3219/(sqrt(2π)σ(M)) > 0`. Pencil proof from
Proposition C₀ with `c_j ≤ 4`, `S_4 ≤ 16σ²`, a uniform lower bound for `σ²`, and
nine zeros below `5M` obtained from the variance identity at two weights; four
elementary inequalities evaluated in Arb (`analytic_leg.py`). It is independent
of the paper's `lem:analytic-leg` (`M ≥ 300`, `q = 3, 4`) and of
`thm:monotone-large`, and uniform in `q`. Refereed in `../REFEREE-5.md` (proved
as written; constant `0.3219` correct). For every real primitive `χ` one also
has `M_0(χ) < 40` (`Q(40) ≥ S_4(40)/64 > 38`), so for a general character there is
always a nonempty range `(M_0(χ), 40)` on which nothing is proved.

## Consequence for the conjecture (`q = 3, 4`)

Theorem B (coupling; `M ≤ M_0`, `M_0 > 3.68` resp. `4.47` certified in
`../cert/aux_results.txt`) + this certificate (`[1, 2000]`) + Theorem E
(`[40, ∞)`) give `δ_q` strictly decreasing on `(0, ∞)`, under GRH(`χ_q`), the
variance identity, the certified zero lists, FLINT/Arb and the pencil proofs
(hypotheses (H1)–(H5));
LI is needed only to read `δ_q` as the density of the race. (Only `M_0 > 1` is
needed for the overlap, and it follows from the certificate's own numbers:
`Q ≤ c_1σ²/8 ≤ 0.053` resp. `0.029` at `M = 1.01`, `METHOD.md` §8.5. Theorem E is
needed only for `M > 2000`.) The direct certificate of `../cert/` re-proves
`[3.6, 320]` resp. `[4.4, 350]` by a different method.

## Caveats

1. Conditional on (H1)–(H5). The program certifies the hypothesis of
   Corollary C′; the corollary, Lemma D and Proposition C₀ are pencil proofs
   (Lemma D, Theorem C's inequality and Corollary C′ were checked by Referee 3;
   Proposition C₀, the expanded proof of C′ and Theorem E by Referee 5).
   All refereeing is internal and AI-conducted.
2. The zero-list certificate and the variance identity are inputs, not
   re-proved here. Theorem B is an input, not re-derived here (Referee 5
   re-derived it). "Fail-closed" does not cover a wrong zero file: that is
   hypothesis (H3); `selftest` detects gross errors only.
3. The comparison with `../cert/` uses printed decimals in floating point: it
   is a cross-check, not a proof.
4. The closed-form code for `(log Λ)′`, `(log Λ)″` is shared with `../cert/`; an
   error there would affect both certificates. Mitigation: `selftest` compares
   the series with the complex-step evaluation and compares `S_4` with an
   enclosure computed from the zero list alone.
5. `q = 3, 4` only. For other characters nothing is certified below `M = 40`,
   and monotonicity is numerically false for conductors 163 and 232 near
   `M ≈ 0.6–1.2` (`../BRIDGE-ATTEMPT.md` §4.3).

## Files

`certify_envelope.py` (prover), `METHOD.md`, `results_q4.txt`, `results_q3.txt`
(one line per cell), `points_q*.txt` (degenerate cells for the table in
`../BRIDGE-ATTEMPT.md` §3), `consistency.py`, `consistency_q*.txt`,
`analytic_leg.py`, `analytic_leg.out`, `ANALYTIC-LEG.md`, `selftest.out`,
`run_all.sh`.

## Review history

* 2026-11-24: certificate, Theorem E and this summary written.
* 2026-11-24, `../REFEREE-5.md` (internal, AI-conducted): no fatal error. Theorem E,
  Proposition C₀, Corollary C′, Lemma C″ and the envelope certificate certified;
  the chain for `q = 3, 4` complete. One repairable item, thirteen minor; all
  applied:
  * R1 (displayed `≥`/`≤` numbers in the proof of Theorem E printed with
    round-to-nearest): `analytic_leg.py` prints lower endpoints rounded down
    and upper endpoints rounded up, `analytic_leg.out` regenerated,
    `ANALYTIC-LEG.md` table and Step 6 replaced by the corrected values.
  * In this file: sub-range maxima of `K` corrected (minor 1); the range of
    "loses 4–17 %" stated (minor 2); scope of the comparison with `../cert/`
    (minor 7) and the full list of what is shared (minor 10); the consequence
    stated with its computational hypotheses (minor 11); `M_0 > 1` and
    `M_0(χ) < 40` recorded.
  * In the code: `suf[-1] = rem0` (minor 3); SHA-256 header and guarded resume
    (minor 8); docstrings on "fail-closed" (minor 9) and on `consistency.py`
    (minor 7). In `METHOD.md`: minors 1, 2, 7, 8, 9, 10. In `ANALYTIC-LEG.md`:
    minors 4, 5, 6, 12. In `../BRIDGE-ATTEMPT.md`: minors 1, 2, 4, 5, 6, 7, 10,
    11, 13.
  * After the code changes both ranges were recomputed from scratch: all 1530
    cell lines identical to the previous files apart from the `secs` field; the
    files differ only by the new `# sha256` header line and timestamps
    (`METHOD.md` §10).
  * The first part of this repair (R1 and the first edits of this file) was
    begun by another agent that was stopped mid-task; those edits were checked
    against the referee's report, the results files and a fresh run of
    `analytic_leg.py`, and kept where correct.
