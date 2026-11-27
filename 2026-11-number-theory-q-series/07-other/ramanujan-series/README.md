# Ramanujan-type series: WZ proofs with certificates

Start with `RESULTS.md` (summary on top; statements, sources and full proofs of Theorems 1–11).

| path | content |
|---|---|
| `RESULTS.md` | what was proved (11 conjectured identities from Z.-W. Sun's lists: 10 apparently new proofs, 1 — T5 — already in Chu 2020), what is machine-checked, what resisted, how novelty was checked |
| `AUDIT-INDEPENDENT.md`, `audit2/` | independent audit of all eleven proofs and of the novelty claims (by another session); `proofs/AUDIT-T3-T4-T11.md`, `proofs/audit1/` = earlier referee pass |
| `OPEN-LIST.md` | index of ~95 conjectural identities open on 2026-11-24, with sources, digits verified, scan results, status |
| `survey/` | the three literature surveys behind `OPEN-LIST.md` and their numerical verification scripts/outputs; `verify_G.py` = Guillera's 32 conjectured half-shifted Ramanujan series (arXiv:1802.02023), all verified to 75+ digits |
| `certificates/` | plain-text WZ certificates (PARI/GP); `check_all.out.txt` = output of `gp -q < T*.gp` |
| `proofs/verify_T*.py` (+ `.out.txt`) | per-theorem verification: WZ equation and auxiliary identities/inequalities exactly (sympy), numerics >= 104 digits; `verify_R4_uniform_bounds.py` = the parameter-uniform bounds of Theorems 6, 8, 9, 10 |
| `lean/` | Lean 4 + Mathlib: `Sun16.lean` (complete formal proof of Theorem 1), `Sun1pi.lean` (complete formal proof of Theorem 4, the 1/pi series), `Sun1pi2.lean` (WZ identities of Theorem 11); `Axioms.out.txt` |
| `code/` | the search machinery: `wz.py` (Gosper/Zeilberger on PARI, seeds of Au), `seeds_extra.py`, `ebisu_parse.py`, `ebisu_seeds.json` (149 seeds from Ebisu's tables), `scan.py`, `run_scan*.py`, `match_hits.py`, `targets.py`, exploration scripts and logs; `canonical_data.py` (canonical hypergeometric data of the eleven series), `chu_family_test.py` (are our identities, plain / bisected / trisected, in the Omega, Phi, 2H2 families of Chu–Zhang / Chu?), `multisect_scan.py` + `match_multisect.py` (signature scan and exact matching for targets grouped in blocks of 2 or 3 terms), `gscan/`, `msect/` (their logs for Guillera's and Sun's targets), `g13_reduction.py` |
| `PROGRESS.md` | chronological log |

Reproduce: `pip install sympy mpmath cypari2`; `apt install pari-gp`; `python3 proofs/verify_T1_sun_16k.py` etc.;
`cd certificates && for f in T*.gp; do gp -q < $f; done`; `cd lean && lake exe cache get && lake build`.
