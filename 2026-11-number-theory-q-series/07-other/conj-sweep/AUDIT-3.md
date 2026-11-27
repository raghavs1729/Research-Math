# AUDIT-3. Adversarial referee report on R06, Theorem 2 (boundary family `f + g = d`)

Date: 26 Nov 2026. Object: `R06-3-sumfree.md`, section "Theorem 2 (the whole boundary family `f + g = d`)"
(lines 83–127). Everything below was recomputed with my own programs (`audit3/`); nothing was taken
from `code/`, `cert/` or `audit2/` (the author's script was only re-run at the end, read-only).

## Verdict

**PROOF COMPLETE.** For all integers `f ≥ 1`, `d ≥ max(6, 2f+1)`, `m = 8d−2f−1`:

    S_{f, d−f, 2d−f} = {f, d−f} ∪ [2d−f, 3d−f−1] ∪ {4d−f} ∪ ⋃_{j≥0} [10d−3f−2+jm, 11d−3f−3+jm],

and `(k, p, m) = (d+3, d, 8d−2f−1)` are the minimal values in the sense of the source (`k`, `p` =
number of 1s in the minimal preperiod / period of the characteristic sequence; equivalently
preperiod and period length of the difference sequence). I found **no error** in the statement or in
the proof: every endpoint is correct, every inequality holds on the whole range (no failing pair
`(f, d)`, including `d = 2f+1` with large `f` and `f = 1, 2` with `d = 6`), the eight types are
exhaustive, the coverings are attained by distinct smaller elements, and the cases exhaust the
complement. The four pairs with `2f+1 ≤ d < 6` all fail, as stated.

One defect of wording outside the proof must be fixed (**D1**, definition of "preperiod length": read
literally, R06's paraphrase would give `k = 4`, not `d+3`); the remaining items (**D2–D5**) are
cosmetic. A novelty caveat is in §3 (the 120 instances with `2 ≤ f`, `d ≤ 25` are already tabulated
in the source).

## 1. Proof, line by line, for general `(f, d)`

Notation: `g = d−f`, `h = 2d−f`, `C = [2d−f, 3d−f−1]`, `E = 4d−f`, `B_j = [10d−3f−2+jm, 11d−3f−3+jm]`,
`T` = right-hand side. All inequalities are linear in `(f, d)`; each was proved on the region
`R = {f ≥ 1, d ≥ 6, d ≥ 2f+1}` exactly (values at the vertices `(1,6)`, `(5/2,6)` and along the
recession directions `(0,1)`, `(1,2)`) and its exact solution set was determined
(`audit3/proof_check.py`, Part A; 46 inequalities, all hold on `R`).

**Induction principle.** For `x > h`: `x ∈ S` iff `x` is not a sum of three distinct elements of
`S ∩ [1, x)` (if `s_n` is the largest entry `< x` and `x` is not such a sum, then `s_{n+1} ≤ x` by
leastness and `s_{n+1} ≥ x` by choice of `n`; the converse is the definition; summands of a 3-sum of
positive integers are `< x`). `S ∩ [1, h] = {f, g, h} = T ∩ [1, h]`. Hence (I)+(II) ⇒ `S = T` by
induction on `x`. Needs `f < g` (⇔ `d ≥ 2f+1`), `g < h` (`d ≥ 1`), `f ≥ 1` (positivity). Correct.

**(I).**
* `x ∈ C`, `x > h`: least 3-sum is `f+g+h = 3d−f > 3d−f−1`. Correct (identity).
* `x = E`: summands lie in `{f, g} ∪ C`. One element of `C`: `d + c ∈ [3d−f, 4d−f−1]`. At least two:
  `≥ f + (2d−f) + (2d−f+1) = 4d−f+1` (uses `f ≤ g`, i.e. `d ≥ 2f`). Correct.
* `x ∈ B_j`: `10d−3f−2 − m = 2d−f−1`, so `B_j ≡ [2d−f−1, 3d−f−2]`; `f, g, E` and `C` are `< m`
  (`E < m` ⇔ `4d ≥ f+2`), so they are their own residues; `U = [2d−f−1, 3d−f−1]`; the residue
  `3d−f−1` belongs to the single element `3d−f−1` of `T` (`f, g < 2d−f−1`, `E > 3d−f−1`). The types
  are the 8 subsets of `{a, b, e}` completed by `U`-elements: exhaustive. `σ < 2m`: the largest
  residue sum is `10d−3f−3` (type `eUU`) and `10d−3f−3 < 16d−4f−2` ⇔ `6d ≥ f`; `σ > 0` trivially;
  hence `σ ∈ [2d−f−1, 3d−f−2] ∪ [10d−3f−2, 11d−3f−3]`.

  | type | printed interval | exact range (distinct elements) | `> 3d−f−2` needs | `< 10d−3f−2` needs |
  |---|---|---|---|---|
  | `UUU` | `[6d−3f−3, 9d−3f−3]` | `[6d−3f−3, 9d−3f−5]` | `3d ≥ 2f+2` | `d ≥ 0` |
  | `aUU` | `[4d−f−2, 6d−f−2]` | `[4d−f−2, 6d−f−3]` | `d ≥ 1` | `4d ≥ 2f+1` |
  | `bUU` | `[5d−3f−2, 7d−3f−2]` | `[5d−3f−2, 7d−3f−3]` | `2d ≥ 2f+1` | `3d ≥ 1` |
  | `eUU` | `[8d−3f−2, 10d−3f−3]` | same | `5d ≥ 2f+1` | identity (**tight**; needs distinctness) |
  | `abU` | `[3d−f−1, 4d−f−1]` | same | identity (**tight**) | `6d ≥ 2f+2` |
  | `aeU` | `[6d−f−1, 7d−f−1]` | same | `d ≥ 0` | `3d ≥ 2f+2` |
  | `beU` | `[7d−3f−1, 8d−3f−1]` | same | `4d ≥ 2f` | `d ≥ 1` |
  | `abe` | `5d−f` | same | `2d ≥ −1` | `5d ≥ 2f+3` |

  Every condition in the table follows from `d ≥ f+1`, `f ≥ 1`; R06's "(using `d ≥ 2f+1`)" is
  sufficient and more than needed. `d ≥ 6` is not used in (I). The exact ranges were computed from
  the multiplicity model (residues `[2d−f−1, 3d−f−2]` unlimited, `3d−f−1` once) and the model was
  validated against brute force over triples of distinct actual elements of `T` (three blocks) for
  `d ≤ 22`. No type meets a target interval for any `f ≤ 30`, `max(6,2f+1) ≤ d ≤ 150`. Correct.

**(II).** Two-sums of distinct elements of `C` are exactly `[4d−2f+1, 6d−2f−3]` (`d ≥ 2`), three-sums
exactly `[6d−3f+3, 9d−3f−6]` (`d ≥ 3`): all integers attained, checked as sets.
* `[3d−f, 4d−f−1]`: `x = f + g + (x−d)`, `x−d ∈ C`, `f < g < x−d`. Correct.
* `[4d−f+1, 10d−3f−3]`: `I1 = f+c_1+c_2 = [4d−f+1, 6d−f−3]`, `I2 = g+c_1+c_2 = [5d−3f+1, 7d−3f−3]`,
  `I3 = c_1+c_2+c_3 = [6d−3f+3, 9d−3f−6]`, `I4 = E+c_1+c_2 = [8d−3f+1, 10d−3f−3]`: all four pairs of
  endpoints correct. Covering: `5d−3f+1 ≤ 6d−f−2` ⇔ `d+2f ≥ 3`; `6d−3f+3 ≤ 7d−3f−2` ⇔ **`d ≥ 5`**;
  `8d−3f+1 ≤ 9d−3f−5` ⇔ **`d ≥ 6`** (the only place where `d ≥ 6` is used; independent of `f`);
  implicitly also `end(I1) ≤ end(I2)` ⇔ `d ≥ 2f` (see D3). Summands are `≤ E < x`; `E ∉ C`,
  `f, g ∉ C`, so the three summands are distinct. Correct.
* gap after `B_j`: `11d−3f−2 − m = 3d−f−1`, so `y ∈ [3d−f−1, 10d−3f−3]`; `ρ ∈ [2d−f−1, 3d−f−2]`.
  `J1 = ρ+f+g = [3d−f−1, 4d−f−2]`, `J2 = ρ+f+c = [4d−f−1, 6d−f−3]`, `J3 = ρ+g+c = [5d−3f−1, 7d−3f−3]`,
  `J4 = ρ+c_1+c_2 = [6d−3f, 9d−3f−5]`, `J5 = ρ+E+c = [8d−3f−1, 10d−3f−3]`: endpoints correct, every
  integer attained (sum of integer intervals). Covering needs `d+2f ≥ 1`, `d ≥ 2f` (ends ordered),
  `d ≥ 2`, `d ≥ 3`. `β ∈ B_j` is `< x` and `> E ≥ s_1, s_2` (`6d ≥ 2f+3`); distinctness clear. Correct.
* Exhaustion: `(h, ∞) ∖ T = [3d−f, 4d−f−1] ∪ [4d−f+1, 10d−3f−3] ∪ ⋃_j [11d−3f−2+jm, 10d−3f−3+(j+1)m]`
  (pieces consecutive and disjoint; needs `E < min B_0`, `d−1 < m`). Checked as a partition. Correct.

**Minimality.** Periodic part of the characteristic sequence: `1^d 0^{7d−2f−1}`, least period `m`.
`E ∈ S`, `E + m = 12d−3f−1 ∉ S` (`11d−3f−3 < 12d−3f−1 < 18d−5f−3` ⇔ `d ≥ −1`, `6d ≥ 2f+3`), and
`c(n) = c(n+m)` for all `n > E`; so the minimal preperiod is `[1, E]` with `d+3` ones. Equivalently
`D_{d+3} = 6d−2f−2 ≠ 7d−2f = D_{2d+3}` and `(1^{d−1}, 7d−2f)` is primitive. Difference sequence
`(d−2f, d, 1^{d−1}, d+1, 6d−2f−2, \overline{1^{d−1}, 7d−2f})`: recomputed, correct.

**Search for a failing pair.** None exists: Part A is exact on `R`; in addition all set-level
claims were evaluated for all 3566 pairs `f ≤ 30`, `max(6,2f+1) ≤ d ≤ 150` (Part B: 0 failures).
Below the range the only claim that fails is the covering of `[4d−f+1, 10d−3f−3]`
(pairs `(1,3), (1,4), (1,5), (2,5)`); for `d = 5` the uncovered point is `8d−3f` (`37`, resp. `34`),
which is exactly the first deviation of `S` from `T`.

## 2. Independent verification

`audit3/sf3.c` (own C program, "pull" algorithm with pair counts and inclusion–exclusion for
distinctness; reproduces the source's `S_{1,2,3}`, `S_{1,2,4}`, `S_{3,7,13}`, `S_{4,7,8}`; agrees with a
naive triple-loop Python implementation on 53 sequences, `naive_cross.py`).

| range | pairs with `d ≥ 6` | window | (A) `S ≠ T` | (B) (I)/(II) violated for `T` |
|---|---|---|---|---|
| `f ≤ 30`, `2f+1 ≤ d ≤ 150` (requested) | 3566 | `[1, 10d−3f−3+12m]` (12 periods) | 0 | 0 |
| `31 ≤ f ≤ 74`, `d ≤ 150` (all remaining `f`) | 1980 | 12 periods | 0 | 0 |
| `f ≤ 149`, `2f+1 ≤ d ≤ 300` (all pairs) | 22346 | 10 periods | 0 | 0 |

(A) = greedy simulation compared with `T`; (B) = direct test of "`x ∈ T` ⇔ `x` is not a sum of three
distinct elements of `T`" for all `h < x ≤ X`. A third route (`proof_check.py`, Part C): number of
representations by FFT inclusion–exclusion (`A_3 − 3A_{21} + 2A_{111})/6`, 4 periods, 3566 pairs:
0 violations; the least number of representations of a non-element is 1. Part D: minimal `(k,p,m)`
of the characteristic sequence equals `(d+3, d, 8d−2f−1)` for all 3566 pairs.
Remark: each simulated pair is also proved independently of R06 by the source's Lemma `lemma`
(`t = 3`, `K = 10d−3f−3`, `M = m`, `Q = 12` or `10`: `(Q−6)M ≥ 2K` ⇔ e.g. `28d ≥ 6f`).

**Pairs with `2f+1 ≤ d < 6`** (exactly four; the pattern fails for all four; `small_cases.py`):

| `(f,d)` | sequence | first deviation | observed `(k,p,m)` | pattern would be |
|---|---|---|---|---|
| `(1,3)` | `S_{1,2,5} = 1,2,5,6,7,11,16,21,26,…` | `16 ∈ S` | `(5,1,5)` | `(6,3,21)` |
| `(1,4)` | `S_{1,3,7} = 1,3,7,8,9,10,15,29,35,36,56,…` | `29 ∈ S` | `(7,3,27)` | `(7,4,29)` |
| `(1,5)` | `S_{1,4,9} = 1,4,9,…,13,19,37,45,46,72,…` | `37 ∈ S` | `(13,10,104)` | `(8,5,37)` |
| `(2,5)` | `S_{2,3,8} = 2,3,8,…,12,18,34,42,43,50,74,…` | `34 ∈ S` | `(208,237,2215)` | `(8,5,35)` |

(`(k,p,m)` observed on `[1, 3000]`, `[1, 6000]`, `[1, 20000]`, `[1, 400000]`; the last agrees with the
source's table row `(2,3,8)`.) R06's sentence "The pattern fails for `(1,3), (1,4), (1,5), (2,5)`, the
only pairs with `2f+1 ≤ d < 6`" is correct.

## 3. Comparison with the source

arXiv:2609.16843, e-print downloaded 26 Nov 2026 to `/tmp/arxiv_2609_16843_audit3/src.bin`
(sha256 `a1002c5d…b1b58c`, gzip of 13 Nov 2026; members `main.tex`, `dlistcor`, `00README.json`), read
only through `python3 -I` + `tarfile` (`audit3/src_*.py`), nothing extracted or executed.

(a) **Confirmed.** `main.tex` l. 997–1004, Theorem `thm:gend>1`: "Let `d ≥ 2`, and if `d = 2` assume
`(f,g) ≠ (1,2)`. Then for every `f ≥ 1` and every `g ≥ f+1` such that `f+g ≥ d+1` the greedy 3-sumfree
sequence `S_{f,g,g+d}` is ultimately periodic modulo `5g+2d+3(f−1)`, with preperiod length `g+f+2` and
period length `g+f`." Its proof (l. 1015, 1021) says the condition `f+g ≥ d+1` is used and that
otherwise "there may be elements in this block that do belong to `S`" (for `f+g = d` that element is
`f+2h = 4d−f = E`). At `f+g = d` its formulas would give modulus `7d−2f−3` and preperiod `d+2`, not
`8d−2f−1` and `d+3`. So `f+g = d` is not covered. (The short form `thm:k=3d>1short` requires `g ≥ max(f+1, d)`.)

(b) **No theorem, conjecture or remark of the source covers the family.** The only statements on
3-sumfree families are `thm:k=1d=1`, `thm:d=1` (`d = 1`), `thm:f=1d>1` (`f = 1`, `g ≥ d+1`),
`thm:gend>1`, and Conjecture `conj:except` (`f = 1` only, with the modulus `8d−1`). The string "f+g"
occurs only in `thm:gend>1`/its proof and in the caption of Table `dat` ("`f+g ≤ d`"). The
introduction only says (l. 160–162) that the `f = 1` family "is just one example of various other
parameter families for which computations of small values readily provide conjectural period and
preperiod lengths". No formula `8d−2f−1` and no statement about `g = d−f` appears.

**Tables** (`src_table.py`, `src_table.log`). Table `short` (`main.tex`, `d ≤ 10`, 84 rows) contains
the 11 boundary triples with `f ≥ 2`, `6 ≤ d ≤ 10`; the full list `dlistcor` (`\input` at l. 1411,
2024 rows, `d ≤ 25`) contains all 120 boundary triples with `f ≥ 2`, `6 ≤ d ≤ 25`. **All 120 rows
agree with `(d+3, d, 8d−2f−1)`**; in particular `(2,4,10) ↦ (9,6,43)`, `(2,8,18) ↦ (13,10,75)`,
`(3,7,17) ↦ (13,10,73)`, `(4,6,16) ↦ (13,10,71)`. The only other boundary row is `(2,3,8)`
(`d = 5`): `(208,237,2215)`, outside the range of Theorem 2 and consistent with §2. Rows with `f = 1`,
`g = d−1` are absent (excluded as "covered" by the conjecture).

**Novelty caveat (N1).** According to the source (§ "Computational results", Theorem
`thm:smallshort`), its tabulated `(k,p,m)` are proved by Lemma `lemma` plus a finite computation.
Hence for the 120 individual triples with `2 ≤ f`, `6 ≤ d ≤ 25` the conclusion of Theorem 2 is already
established in the source; new in R06 are the general formula, the cases `d > 25` or `f = 1`, and a
computation-free proof. R06's wording "the source's tables of unexplained small cases" is accurate
but should mention this.

## 4. Defects and remarks (none affects `S = T`)

* **D1 (R06 l. 6–8, definition; affects the clause "(minimal) preperiod length `d+3`" of both
  Theorems as worded).** R06 defines: "`S` consists of `k` entries `a_1 < … < a_k` and all integers
  `b_i + jm`; lengths are minimal". Read literally (no condition `a_k < b_1`) the least `k` is **4**,
  not `d+3`: since `C ∖ {3d−f−1} = [2d−f, 3d−f−2]` lies in the residue classes of the blocks,
  `T = {f, d−f, 3d−f−1, 4d−f} ∪ {b + jm : j ≥ 0, b ∈ [2d−f, 3d−f−2] ∪ {10d−3f−2}}`
  (`audit3/decomp_check.py`, 3566 pairs). The source means the sequential notion (`k` = number of 1s
  in the minimal preperiod of the characteristic sequence = preperiod length of `D`; e.g. it counts
  `k = 6` for `S_{1,2,3}` although only `1, 5, 13` lie outside the periodic classes), and under that
  definition `d+3` is correct and proved (R06 argues via `D`). Fix: add "`a_k < b_1`" or define
  `k, p` through the characteristic/difference sequence, as `AUDIT-2.md` B1 does.
* **D2 (l. 107, types `UUU`, `aUU`, `bUU`).** The printed upper ends `9d−3f−3`, `6d−f−2`, `7d−3f−2` are
  not attained (true maxima `9d−3f−5`, `6d−f−3`, `7d−3f−3`, since only one element has residue
  `3d−f−1`). As enclosures of `σ` they are valid and still `< 10d−3f−2`; the argument is unaffected.
* **D3 (l. 113–115 and 117–119).** "They cover the interval because …" lists only the conditions
  `start(I_{i+1}) ≤ end(I_i)+1`; one also needs the right ends to be increasing, i.e.
  `6d−f−3 ≤ 7d−3f−3` ⇔ `d ≥ 2f` for `I1, I2` (and `J2, J3`). True under the hypothesis, unstated.
* **D4 (l. 110).** "(using `d ≥ 2f+1`)": correct, but (I) for the blocks needs only `d ≥ f+1`; the
  hypothesis `d ≥ 2f+1` enters through `f < g` (seeds, least sums) and D3, and `d ≥ 6` only through
  `8d−3f+1 ≤ 9d−3f−5`. For `d = 5` exactly one integer, `8d−3f`, is uncovered.
* **D5 (l. 98, 120).** The induction principle, the two-/three-sum facts for `C`, the exhaustion of
  the complement and the minimality argument are inherited from the proof of the Theorem by
  reference ("As for the Theorem", "Minimality as before"); all four transfer verbatim (checked
  above). `σ < 2m` is asserted without the one-line bound `σ ≤ 10d−3f−3`.
* The author's check `code/sumfree3_family_verify.py 14 100` was re-run (`python3 -B`, read-only) and
  reproduces `cert/sumfree3_family_verify.txt`; its count 1186 is correct.

## 5. Programs and logs (all in `audit3/`)

| file | purpose |
|---|---|
| `sf3.c` (`gcc -O3 -march=native`) | simulation `seq`, representation counts `rep`, family check (A)+(B) |
| `family_f30_d150.log`, `family_f31-74_d150.log`, `family_f*_d300.log` | results of §2 |
| `naive_cross.py`, `.log` | naive reference implementation vs `sf3` |
| `proof_check.py`, `proof_check_f30_d150.log` | Part A inequalities, B set claims, C FFT counts, D minimal `(k,p,m)` |
| `small_cases.py`, `.log` | the four pairs with `d < 6`; spot checks of 10 source table rows (all agree) |
| `decomp_check.py`, `.log` | defect D1 |
| `src_list.py`, `src_grep.py`, `src_table.py`, `src_table.log` | read-only inspection of the arXiv source |
