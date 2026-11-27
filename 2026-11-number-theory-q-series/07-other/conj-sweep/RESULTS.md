# conj-sweep: final results (session closed 26 Nov 2026)

Sweep of conjectures printed in 2024–2026 arXiv papers (math.NT / math.CO / math.CA): each was
attacked by falsification first, then proved where the structure allowed.

* **PROVED** = complete proof in the file (possibly a reduction, by cited theorems, to a finite exact
  computation that was carried out). **DISPROVED** = exact counterexample to the statement as printed.
  **tested** = survived computation only; nothing is claimed.
* **Audit** = report of an independent adversarial referee (fresh agent, own programs).
* **Novelty** = what the searches of **26 Nov 2026** found (`NOVELTY-1.md`, `NOVELTY-2.md`, the novelty
  sections of the audits, GitHub issue/PR searches, `1337isnot1337/solved-conjectures`). "No prior
  found" means only that; none of the third-party proofs mentioned was read or checked here.
  Several automated projects are working the same papers, so every status is as of that date.

## A. Settled here, no prior proof or correction found (checked 26 Nov 2026)

| # | statement (source) | outcome | file | audit | novelty |
|---|---|---|---|---|---|
| R06 | 3-sumfree sequence `S_{1,d−1,2d−1}`: "preperiod `d+3`, period `d` modulo `8d−1`, except `d = 3, 5`" (van Berkel–Bosma, arXiv:2609.16843, `conj:except`) | **DISPROVED as printed** (modulus is `8d−3`; `d = 4` is a further exception). Corrected statement **PROVED** for all `d ≥ 6`. **Theorem 2:** `S_{f,d−f,2d−f}` is periodic mod `8d−2f−1`, preperiod `d+3`, period `d`, for all `f ≥ 1`, `d ≥ max(6, 2f+1)` (the boundary family `f+g = d` not covered by the source's theorem) | `R06-3-sumfree.md` | `AUDIT-2.md` (`f = 1`), `AUDIT-3.md` (general `f`): proof complete | no prior correction or proof found; source v2 (11 Nov 2026) still prints `8d−1`. Caveat: the source's tables already establish the 120 instances `f ≥ 2`, `6 ≤ d ≤ 25` by computation. A public issue (`ndesmo/leaning` #24, 25 Nov) targets the printed form |
| R08 | 2-colour overpartitions: `p̄_ℓ(4n+3) ≡ 0 (8)` for all `ℓ ≡ 0,1 (4)`; `p̄₉(2^α(8n+7)) ≡ 0 (16)`; `p̄₈(32n+28) ≡ 0 (32)`; `p̄₆(24n+17), p̄₆(24n+23) ≡ 0 (64)`; `p̄₉(36n+33) ≡ 0 (128)` (arXiv:2607.16608, three closing conjectures) | **PROVED** (first two by sums-of-squares arguments, the others by Radu's finite check) | `R08-two-colour-overpartitions.md` | `AUDIT-2.md`: proof complete | no citing paper, issue, PR or record found (id and keyword searches) |
| R03a | `A₆(12n+10) ≡ 0 (3⁵)` (exact), `A₃(12n+7), A₃(12n+11) ≡ 0 (3⁴)`, `A₄(6n+5), A₇(12n+11) ≡ 0 (3³)`, `A₈(18n+15) ≡ 0 (3⁵)` (arXiv:2611.05623, Conj. 6.3, 6.4) | **PROVED** (Radu certificates) | `R03-radu-congruences.md` | `AUDIT-2.md`: complete (certificates re-run with the referee's own PARI implementation) | no proof found; open public issues with numerics only. High risk: the paper's Conj. 6.1, 6.2 were proved by others on 23 and 26 Nov |
| R05 | Hankel determinants of `M(x)^r`, Motzkin (arXiv:2502.21050, Conj. 1, 2) | **DISPROVED as printed**: sign error, `H₉(M⁹) = H₁₀(M⁹) = +256` (also `r = 12, 21, 24`); `H₇ + H₁₀ = 0` at `r = 4`. Repaired forms tested `r ≤ 33` | `R05-motzkin-hankel.md` | `AUDIT-2.md`: confirmed | the `r = 9, 12, 21, 24` sign error was not found elsewhere; the `r = 4` vanishing sum was already noted (trureturing PR #6686, 9 Sep 2026). Minor |

## B. Settled here independently, but NOT NEW (prior public claim exists)

| # | statement (source) | outcome | file | audit | prior work found 26 Nov 2026 |
|---|---|---|---|---|---|
| R02 | four overpartition congruences mod 7, 17, 19, 23 (X. Wei, arXiv:2603.08510, Conj. 1, 2 and closing remarks) | **PROVED** (Hecke `T(p²)` on `M_k(Γ₀(4))` + Sturm; mod 19 also Θ-operator, weight 719/2; mod 7 by Radu) | `R02-overpartition-congruences.md` | `AUDIT-R02.md`: A, B, C complete; D complete after one citation was replaced by a direct proof | S. K. Jha, Research Square preprint 10.21203/rs.3.rs-11336379/v1 (26 Oct 2026, unrefereed) claims all four |
| R07 | CRIM: "Every 𝒫-position is of even rank"; "a meld of losing thick hooks is losing"; `conj:bloby`; "if" half of `conj:pself` (arXiv:2606.16828) | **PROVED** via a 17-state automaton for all losing positions, exact closure check | `R07-crim.md` | `AUDIT-2.md`: complete | Z. Zeng's boundary-word classification (29 Aug 2026, GitHub archive `ArtificialZeng/AI-Has-Taste`); Lean PRs #73, #76 of `umaia1234/agentic-conjectures` (25 Oct 2026). Our automaton agrees with Zeng's criterion on all words of length ≤ 22 |
| R07b | CRIM, misère play: all losing positions (addendum, not a printed conjecture) | **PROVED** (16-state automaton) | `R07b-crim-misere.md` | not refereed separately | a misère classification is claimed in the same Lean PR #76 |
| R10 | maximal traffic with an antidiagonal obstruction: `ρ(n) < 1` for all `n ≥ 496` (arXiv:2609.01562, `conj:496`) | **PROVED**: `ρ(n) < (4/√(2πe))e^{4.128/√n}` for `n ≥ 12`, plus exact check `496 ≤ n ≤ 16100` | `R10-traffic-496.md` | `AUDIT-4.md`: complete | Lean proof in `the-omega-institute/trureturing` (issue #8249, PR #8934, 20 Sep 2026) |
| R09 | P-NIM: partitions between the staircase and a heavy rectangle are heavy (arXiv:2506.04991, `conj:heavychoppedrect`) | **DISPROVED**: `a = 8`, `b = 7`, `λ = [9,9,7,7,6,5,4,2]`, `SG = 3 ≠ 16` (two independent solvers) | `R09-pnim-heavy.md` | not refereed separately | same failure `(8,7)` in trureturing #8467 (17 Sep 2026) |
| R01 | `[a_1]_q⋯[a_k]_q[b]_{q^r}`: necessity clause "if `k ≤ 3` or `r ≤ 3`" (arXiv:2605.12822) | **DISPROVED**: `(1+q)⁶(1+q³)` | `R01-qint-unimodality.md` | hand check + two programs + Lean | trureturing PR #6805 (10 Sep 2026) |
| R04 | `det[lcm(i,j)+gcd(i,j)]_{n×n} = 0` iff binary expansion of `n` begins with `10` (arXiv:2601.16626) | "if" direction **PROVED** (row relation); "only if" open | `R04-lcm-gcd-eigenvalue.md` | not refereed | "if": trureturing issue #9180 (20 Sep 2026) |
| R03b | `ped(225n+t) ≡ 0 (mod 192)` (arXiv:2503.06228) | re-proved (Radu); **was never open** — proved in its own source (misreading of an extracted `Conjecture` environment) | `R03-radu-congruences.md` | — | the source's Theorem 1.2 |

Also noted in `R07-crim.md` (addendum): the source's printed Sprague–Grundy conjecture for the rectairs
`R^k_{r,r−1}` has the wrong parity at `k = r−2` (value 3 iff `r` is even), checked for `7 ≤ r ≤ 15`;
novelty not searched.

## C. Tested only (all survived; nothing proved)

| # | conjecture (source) | range reached | program / output | note (26 Nov 2026) |
|---|---|---|---|---|
| T1 | `13/66·π(s) ≤ π(a,b) ≤ ½π(s)` (Chen–Zhu, arXiv:2506.03620) | 676 729 coprime pairs `3 ≤ a ≤ 400`, `b ≤ 3000` | `code/axby.py` | upper bound proved for `a ≤ 8` in arXiv:2506.03625 |
| T2 | least pure period of 3-move subtraction games ∈ `{a+b, c+a, c+b}` (arXiv:2609.05358) | 14 671 617 triples `a ≤ 120`, `b ≤ 300`, `c ≤ 800` | `code/subgame3.c`, `cert/subgame3_test_big.txt` | partial results by others (Zenodo 24536124, 22661435) |
| T3 | unimodality of `rmax` on separable permutations (arXiv:2404.18517, Conj. 1–3) | `n ≤ 700` | `code/separable.py` | since proved by others in Lean (Nov 2026) |
| T4 | `COS, CES(p^{2α}n + (7p^{2α}+1)/8) ≡ 0 (p^α)`, `p ∈ {3,5,7}` (arXiv:2611.05196) | `(p,α)` up to `(3,4), (5,3), (7,2)` | `code/cubic_partsum.py` | no proof found |
| T5 | Grundy values of KING on generalised staircases (arXiv:2501.14640, two conjectures) | odd `r ≤ 31`, even `r ≠ c ≤ 30`, `k ≤ 130` | `code/king_gs.c` | not searched beyond issue search (0 hits) |
| T6 | Divisor Nim: `g(P) ≤ 2h_min` (arXiv:2610.06925) | 2 heaps ≤ 1200, 3 ≤ 180, 4 ≤ 60, 5 ≤ 30 | `code/divisor_nim.py` | a public preregistration exists (trureturing) |
| T7 | wall/sink nimber duality for `{m, m+δ, 2m+δ}` (arXiv:2601.18715) | 4800 pairs `m ≤ 40`, `δ ≤ 120` | `code/sink_duality.py` | not searched |
| T8 | closed form of `∫₀^∞ x⁹/(cos x + cosh x)⁶ dx` (arXiv:2401.01385) | agreement to 68 digits | `code/berndt_integral.py` | issue search: 0 hits |
| T9 | Berkovich–Dhar positivity (arXiv:2503.17571, Conj. 4.1, 2.6) | all `L ≤ 160` | `code/berkovich_dhar.py` | issue search: 0 hits |
| T10 | `N_k = 26` for all `k ≥ 30`, log-concavity of `k`-regular partitions (arXiv:2406.10987) | `30 ≤ k ≤ 400`, `n ≤ 4000` | `code/kregular_logconcave.py` | issue search: 0 hits |
| — | 2-adic Conj. 6.1, 6.2 of arXiv:2611.05623 | `k ≤ 2100`, `n < 160` | `code/ds_2adic.py` | proved by others 23 and 26 Nov 2026 |
| — | repaired forms in R01 and R05; "only if" halves of R04 and of CRIM `conj:pself` | see the files | | |

## D. Other files

* `CANDIDATES.md`: the harvest (verbatim statements, status of each candidate, why some were dropped).
* `PROGRESS.md`: chronological log, including the search lesson: a GitHub search for a bare arXiv id
  misses issues that write `<id>v1`; search `"<id>v1"`, author names and keywords, in issues and pull
  requests, and grep `CONTENTS.md` of `1337isnot1337/solved-conjectures`, before investing in a target.
* `lean/ConjSweep.lean` (Lean 4 core, no Mathlib): the R01 counterexample (kernel `decide` for the
  failed hypothesis), and by `native_decide` the R05 determinants, the R06 sequences for `d = 4, 6`
  and the closure check of the R07 automaton.
* `code/`, `cert/`: programs and their outputs; `audit/`, `audit2/`, `audit3/`, `audit4/`: the referees'
  independent programs.
* A second harvest restricted to conjectures of standing was started and cancelled on the
  coordinator's instruction (no `HARVEST-2.md` was produced).

No human has read any of this.
