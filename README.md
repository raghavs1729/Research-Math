# Research-Math

Mathematics research by Raghav Sharma, carried out with AI (Anthropic's Claude models).
Every project states on its first page what was produced by AI and how it was checked.
Corrections are welcome: raghavsharma1729@gmail.com.

| Project | What it contains |
|---|---|
| [`2026-11-number-theory-q-series/`](2026-11-number-theory-q-series/) | Proofs of published conjectures in q-series and number theory (Kanade–Russell–Tsuchioka–Warnaar; Nahm sums for E6/E7/E8; the Ramanujan Machine's unidentified ζ(5) constants; Zudilin's relations; Ramanujan's τ), results on Chebyshev's bias under GRH alone, with Lean 4 formalizations and every referee report |
| [weighted-prime-races](https://github.com/raghavs1729/weighted-prime-races) (separate repository) | Weighted Chebyshev prime races: preprint [doi:10.5281/zenodo.21947101](https://doi.org/10.5281/zenodo.21947101), code and data |

## Verify it yourself

The Lean proofs are in `2026-11-number-theory-q-series/lean/`. With
[elan](https://github.com/leanprover/elan) installed:

```bash
cd 2026-11-number-theory-q-series/lean/lean-krtw
lake exe cache get && lake build && lake env lean scripts/PrintAxioms.lean
```

Details in [`2026-11-number-theory-q-series/VERIFY.md`](2026-11-number-theory-q-series/VERIFY.md).

## License

Code: MIT. Text, proofs and data: CC BY 4.0.
