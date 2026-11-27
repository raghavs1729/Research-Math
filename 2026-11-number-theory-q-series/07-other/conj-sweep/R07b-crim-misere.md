# R07b. CRIM under misère play: complete solution of the losing positions (addendum to R07)

Not a conjecture of the source (arXiv:2606.16828), but the source studies the misère game through
"Conway pairs" and asks for characterisations; the method of R07 settles the misère outcome as well.
Misère play: the player who cannot move **wins**; so the empty partition is a winning (`𝒩`) position
and a nonempty position is losing iff none of its options is losing. Encoding as in R07 (boundary
words over `{E,N}`, a move deletes one letter and re-canonicalises).

## Theorem 2
Let `B⁻` be the deterministic automaton with states `0,…,15`, start state `0`, and the table

| state | a word reaching it | accepting | on `E` | on `N` |
|---|---|---|---|---|
| 0 | `ε` | no | 1 | 0 |
| 1 | `E` | no | 2 | 3 |
| 2 | `EE` | no | 4 | 4 |
| 3 | `EN` | yes | 5 | 1 |
| 4 | `EEE` | no | 2 | 6 |
| 5 | `ENE` | yes | 7 | 8 |
| 6 | `EEEN` | yes | 9 | 4 |
| 7 | `ENEE` | yes | 10 | 11 |
| 8 | `ENEN` | no | 1 | 12 |
| 9 | `EEENE` | yes | 6 | 2 |
| 10 | `ENEEE` | yes | 7 | 13 |
| 11 | `ENEEN` | no | 14 | 13 |
| 12 | `ENENN` | no | 14 | 8 |
| 13 | `ENEEEN` | no | 15 | 11 |
| 14 | `ENEENE` | no | 11 | 11 |
| 15 | `ENEEENE` | no | 13 | 7 |

A partition `λ` is a losing position of misère CRIM if and only if `B⁻` accepts `w(λ)`.

*Proof.* As for Theorem 1 of R07, with "the empty word is rejected" in place of "accepted":
(i) `δ(0,N) = 0`; (ii) `δ(q,E)` is accepting iff `q` is, for every state `q`; (iii) state 0 is not
accepting — so `B⁻(v) ⇔ B⁻(canon(v))` for every word and words with empty `canon` are rejected;
(iv) for every word `w` starting with `E` and ending with `N`: `B⁻(w)` iff no one-letter deletion of
`w` is accepted (exhaustive search over the 162 reachable configurations). Induction on `|w|`. ∎
(`code/crim_misere_dfa.py`, output `cert/crim_misere_cert.txt`; the automaton was learned from
misère Grundy values of words of length `≤ 18`.)

## Corollaries
1. Every misère losing position has even rank (no accepted canonical word of odd length; see the
   table of counts in `cert/crim_misere_cert.txt`). `[1]` (word `EN`) is misère-losing.
2. Joint outcome. With Theorem 1 of R07, for each `L = rows + columns` the numbers of partitions that
   are losing in both, only in normal, only in misère, in neither play are, for
   `L = 4, 6, 8, 10, 12, 14`: `(1,1,2,0)`, `(7,4,3,2)`, `(32,8,8,16)`, `(142,18,18,78)`, `(600,40,40,344)`,
   `(2472,88,88,1448)`; for odd `L` every position is winning in both.

## What this does not give
`conj:hookswap` of the source ("`λ` has Conway pair `(0,1)` iff `λ` is an even-by-even hook") needs the
set `{G⁻ = 1}`, whose Nerode index does not stabilise in our experiments (67 classes and growing with
suffixes of length 9), so the closure method does not apply. Experimentally the set
`{G⁺ = 0, G⁻ = 1}` has a stable 8-state description and `{G⁺ = 1, G⁻ = 0}` (Conway pair `(1,0)`, asked
for in the source) a stable 14-state description (words of length up to 20), consistent with
`conj:hookswap`; these are observations, not proofs.

## Novelty
Probably NOT NEW: a Lean statement `misereP_iff_isP` (a misère P-position classification) is claimed
in the open pull request #76 of `umaia1234/agentic-conjectures` (25 Oct 2026; `NOVELTY-2.md`, item 1).
It was not read here; Theorem 2 is an independent proof.
