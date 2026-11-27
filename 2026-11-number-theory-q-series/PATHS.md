# Path map

The files were written in a working repository whose folder names differ from this public layout.
When a file refers to a path on the left, look on the right.

| Path in the files | Path in this folder |
|---|---|
| `monotonicity/lean-krtw`, `../lean-krtw` | `lean/lean-krtw` |
| `monotonicity/lean-nahm`, `monotonicity/lean-tau`, `monotonicity/lean-allx`, … | `lean/lean-nahm`, `lean/lean-tau`, `lean/lean-allx`, … |
| `monotonicity/rr-identities/lean`, `../rr-identities/lean` | `lean/rr-identities/lean` |
| `monotonicity/rr-identities/…`, `../rr-identities/…` (proofs, code, certificates, referee reports) | `02-exceptional-nahm-sums/…` (the KRTW proofs and their referee reports are also copied to `01-krtw-identities/`) |
| `code/krtw_referee/` (cited in `01-krtw-identities/KRTW-F-moves.md`) | `02-exceptional-nahm-sums/code/krtw_referee/` |
| `monotonicity/rm-open/…` | `03-ramanujan-machine-zeta5/…` |
| `monotonicity/density/`, `density2/`, `density3/` | `04-chebyshev-bias-grh-only/density*/` |
| `monotonicity/chebyshev1853`, `explicit`, `grhonly`, `cert*`, `nonmonotone*`, `hkmn23`, manuscripts | `05-weighted-races-grh/…` |
| `monotonicity/tau-race`, `monotonicity/open-misc/tau`, `monotonicity/open-misc/tau2` | `06-ramanujan-tau/tau-race`, `tau-p4`, `tau-p2` |
| `../CORE-NOTES.md`, `../explicit/EXPLICIT.md` (cited in `06-ramanujan-tau/tau-race/`) | `05-weighted-races-grh/CORE-NOTES.md`, `05-weighted-races-grh/explicit/EXPLICIT.md` |
| `monotonicity/conj-sweep`, `ramanujan-series`, `akbary-francis`, `ramanujan-ap` | `07-other/…` |

To build a Lean project, `cd` into its folder under `lean/` (not `monotonicity/…`).
Paths under `/tmp/...` in referee reports are the referees' scratch locations; those scripts were not archived unless stated.
`PROGRESS.md` and `HANDOFF.md` files are working logs and are not published; references to them, and to
"the parent branch" or "the parent session", can be ignored. The review history of
`04-chebyshev-bias-grh-only/density/DENSITY.md`, which those logs held, is published as `density/REVIEW.md`.
