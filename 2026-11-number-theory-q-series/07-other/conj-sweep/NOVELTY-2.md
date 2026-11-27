# NOVELTY-2: prior-art search, batch 2 (eight items)

Search date: 2026-11-26. Tools: WebSearch, WebFetch, `curl` (arXiv API and arXiv HTML, Semantic Scholar
API, Crossref, Zenodo API, GitHub public issue-search API, `raw.githubusercontent.com`, whataifound.org
`/data/entries.json`). No papers were stored in the repository.

Limits of this search:
* The GitHub REST API for repository contents is blocked in this session (HTTP 403 from the proxy). Issue
  and PR pages were therefore read through WebFetch, which returns a summary written by a small model. Text
  quoted from GitHub pages below is as that summary gave it. Exceptions, read directly as raw files:
  Zeng's LaTeX source (item 1), the CKZ result README (item 6), and the README/CONTENTS files of the
  projects listed under "Automated-conjecture projects checked".
* I read the proofs found only far enough to see what they claim. I did not check any of them.

## Summary table

| # | Source (latest version) | Verdict |
|---|---|---|
| 1 | arXiv:2606.16828 **v1**, 15 Jun 2026 | **Already settled (unrefereed)**: main conjecture and a complete P-position classification by Z. Zeng, preprint dated 29 Aug 2026 (GitHub, 16 Sep 2026). Stack, meld and parity-self-conjugate conjectures: Lean proofs claimed in GitHub PRs, 25 Oct 2026. |
| 2 | arXiv:2609.16843 **v2**, 11 Nov 2026 | **No proof/disproof found.** v2 still prints modulus 8d-1 and exceptions d = 3, 5. A public issue opened 25 Nov 2026 targets the conjecture as printed. |
| 3 | arXiv:2502.21050 **v1**, 28 Feb 2025 | **No proof/disproof found** for Conjecture 3. No corrected version. A GitHub PR (9 Sep 2026) refutes the printed form of Conjecture 7 only and says it did not check Conjecture 3. |
| 4 | arXiv:2607.16608 **v1**, 18 Jul 2026 | **No proof/disproof found.** No citing paper, no public issue. |
| 5 | arXiv:2609.05358 **v2**, 16 Sep 2026 (Conjecture 56 in v2; 54 in v1) | **No proof/disproof found; partial results exist**: S. Levinstein (Zenodo, 25 Nov 2026) proves Conj. 57, Conj. 52 <=> Conj. 56, and "two partial period exclusions toward Conjecture 56"; Z. Zeng (Zenodo, 8 Sep 2026) settles the family {2,5,c}. |
| 6 | arXiv:2404.18517 **v1**, 29 Apr 2024 (Discrete Appl. Math. 355 (2024), 169-179) | **Already settled (unrefereed)**: Conjecture 2 by `the-omega-institute/trureturing` (Lean, merged 17 Nov 2026); Conjectures 1 and 3 by `1337isnot1337/solved-conjectures` (Lean, 21-24 Nov 2026). |
| 7 | arXiv:2611.05196 **v1**, 3 Nov 2026 | **No proof/disproof found.** trureturing issue #26479 (4 Nov 2026) closed "unresolved". |
| 8 | arXiv:2611.05623 **v1**, 3 Nov 2026 | **No proof/disproof found** for Conj. 6.3 and 6.4 (unchanged since NOVELTY-1). 6.1 and 6.2 have claimed proofs. |

Version check: `https://export.arxiv.org/api/query?id_list=...` for all eight. Only 2609.16843 and
2609.05358 have a v2.

---

## 1. arXiv:2606.16828 (Bašić, Gottlieb, Krnc, CRIM)

**Latest version.** v1, 15 Jun 2026 (24 pages). No later version.

**Searched.**
* Semantic Scholar citations: empty.
* arXiv API: `all:CRIM AND (partitions | "boundary word" | losing)`, `au:Zeng_Zijian`, `au:Gottlieb_Eric`,
  `"Young diagram" AND game AND (row|column) AND (impartial|Grundy|"P-positions")`.
* WebSearch: title of Zeng's preprint; `"CRIM" game integer partitions "boundary word" ... "even rank" conjecture proved`.
* Zenodo API: `CRIM AND losing`, `"Losing Positions in CRIM"`, `"2606.16828"`, `title:CRIM AND (game OR partition*)`.
* GitHub issue search: `2606.16828`, `CRIM Bašić Gottlieb Krnc`, `CRIM "even rank"`, and the same inside
  `the-omega-institute/trureturing`, `umaia1234/agentic-conjectures`, `1337isnot1337/solved-conjectures`.

**Found.**
* **Zijian Zeng, "A Boundary-Word Classification of Losing Positions in CRIM"**, dated 29 Aug 2026.
  Not on arXiv, Zenodo or Crossref. It is item BW2-B02-P10 of the public GitHub archive
  `ArtificialZeng/AI-Has-Taste` ("450+ open problems in mathematics investigated by AI agents").
  * PDF: https://github.com/ArtificialZeng/AI-Has-Taste/blob/main/pdfs/BigWIN2/BW2-B02-P10_A_Boundary_Word_Classification_of_Losing_Positions_in_CRIM.pdf
  * Source (read in full): https://raw.githubusercontent.com/ArtificialZeng/AI-Has-Taste/main/reproducible_source_code/BigWIN2/BW2-B02-P10/source/paper/main.tex
  * The only commit shown for the PDF is `f3ae144`, 16 Sep 2026 ("Archive 223 current-workflow manuscripts as
    BigWIN2 with paired sources").
  * Content. A partition is encoded by its boundary word `w = 0^{a_1} 1^{b_1} ... 0^{a_k} 1^{b_k}`
    (`|w|` = rows + columns). A move deletes one letter, then strips leading 1s and trailing 0s. Let `L` =
    length of the maximal prefix of pairs with `a_i = 1`, `b_i` odd; `R` = length of the maximal suffix of
    pairs with `a_i` odd, `b_i = 1`; `M = min(k, L+R)`.
    **Theorem 1.1: `w` is a P-position iff `|w|` is even and `M(w)` is even.**
    Corollary 1.2: every P-position has even rank (the main conjecture). Proof by induction on `|w|`, about
    3 pages; computer checks (all partitions of size <= 50, all words of length <= 22) are stated to be
    checks only.
  * This criterion is a regular condition on the boundary word, i.e. the same kind of statement as our
    17-state automaton. The two descriptions can be compared mechanically.
  * The preprint does not mention the meld, stack or parity-self-conjugate conjectures.
* **`umaia1234/agentic-conjectures` PR #73**, "Add crim-hook-stacks: Lean proofs of CRIM conj:bloby and
  conj:main", opened 25 Oct 2026, still open (not merged).
  https://github.com/umaia1234/agentic-conjectures/pull/73
  Lean formalization of Zeng's classification (`isP_ofPairs_iff`), the even-rank conjecture
  (`mainConjecture_holds`), and the losing-hook stacking conjecture ("conj:bloby", = the stack conjecture)
  by a second-player strategy that reduces it to the even-rank statement. Stated sorry-free, standard
  axioms only. "No human has reviewed this work, no novelty is claimed." It also notes an erratum: the
  paper's three-part theorem as printed includes `[1,1,1]`, which is an N-position.
* **`umaia1234/agentic-conjectures` PR #76**, "Add crim-pself-meld: Lean proofs of the CRIM meld and pself
  conjectures", opened 25 Oct 2026, still open, stacked on #73.
  https://github.com/umaia1234/agentic-conjectures/pull/76
  Claims: every meld of losing thick hooks, in any bracketing, is a P-position (`meldConjecture_holds`);
  the parity self-conjugate conjecture (`pselfConjecture_holds`); a misère P-position classification
  (`misereP_iff_isP`); only one direction of the hook-swap (Conway pair) conjecture.
* Partial results that do not reach the full conjecture:
  * S. S. Mehat, "Losing positions of CRIM for partitions with at most six parts", artifact
    10.5281/zenodo.22070744, 23 Aug 2026: all P-positions with at most six parts (or largest part <= 6).
  * H. Ma, W. Zhang, 10.5281/zenodo.24406584, 20 Nov 2026: alternative proof of the four-part
    classification; even rank and hook-stacking only "with at most four rows or columns".
  * trureturing PR #26523 (merged 5 Nov 2026): Lean proof of the four-part classification (Question 1);
    says "Positions with five or more rows ... remain open".
    https://github.com/the-omega-institute/trureturing/pull/26523
  * trureturing issue #8292 (16 Sep 2026, open): meld conjecture, a bounded search that produced nothing.
  * trureturing PRs #8914 and #6828: refute printed Sprague-Grundy conjectures of the paper (Conjecture 2
    at r = 1). Not our statements.
* Earlier literature (LCTR arXiv:2207.04990, 2207.05599; Impartial Chess on partitions 2501.14640; Nim on
  partitions 2506.04991; Row Impartial Terminus 2509.03390): different games. No complete solution of
  row/column removal found before 2026.

**Verdict.** Already settled, in unrefereed form. The main conjecture and a complete characterisation of
P-positions as a regular condition on the boundary word: Zeng, 29 Aug 2026 (public on GitHub by
16 Sep 2026). Stack and meld conjectures, and the parity self-conjugate conjecture: Lean proofs claimed in
PRs #73 and #76, 25 Oct 2026 (open PRs, no human review). Our automaton proof is an independent
confirmation, not a first proof. Nothing is on arXiv or in a journal.

---

## 2. arXiv:2609.16843 (van Berkel, Bosma, t-sumfree sequences)

**Latest version.** v2, 11 Nov 2026 (41 pages; v1 of 15 Sep 2026 was much shorter). The conjecture on
`S_{1,d-1,2d-1}` appears only in v2, as Conjecture 4. In v1, "Conjecture 4" is the folklore conjecture that
every 3-sumfree sequence is ultimately periodic. v2 text, from the arXiv HTML:

> Conjecture 4. For every d >= 3 the 3-sumfree sequence S_{1,d-1,2d-1} is ultimately periodic, and, with the
> exception of d = 3 or d = 5, the preperiod length is d+3 and the period length is d modulo 8d-1.

The authors add that they are "convinced the following result is correct" but "have not written down a
full proof". So the modulus 8d-1 and the exception list {3, 5} are still as printed.

**Searched.**
* Semantic Scholar citations: empty (after one 429).
* arXiv API: `sumfree AND (periodic|periodicity) AND (sequence|sequences)`; `au:Bosma_W`.
* GitHub issue search: `2609.16843`, `sumfree Berkel Bosma`, `repo:the-omega-institute/trureturing sumfree`.
* Zenodo API: `"2609.16843"`, `sumfree AND periodic*`.

**Found.**
* `ndesmo/leaning` issue #24, "P4: Periodicity of 3-sumfree sequences S_{1,d-1,2d-1}", opened 25 Nov 2026,
  open, labelled `approved`, `priority:high`. https://github.com/ndesmo/leaning/issues/24
  It targets the conjecture as printed (8d-1, exceptions 3 and 5), plans "a block-induction proof modelled
  on the paper's Theorem 19", and has no result, comment or PR. Its `PROBLEMS.md` entry does not notice
  that the modulus is wrong.
* trureturing: PR #26092 / #28953 (Nov 2026) formalize Theorem 19 of the paper. Issue #34312 and PR #34359
  (25-26 Nov 2026) refute a printed preperiod formula of the companion paper arXiv:2609.18522
  (Conjecture 9, 2-sumfree). None mentions Conjecture 4, 8d-3, or d = 4.
* arXiv:2606.17447 (O. Shtrezi, "The greedy 3-sumfree sequence S_{1,g,g+1}"): a different family.
* Zenodo 22758978 (N. Wierenga, "A nonperiodic greedy sum-free set"): different object.

**Verdict.** No proof or disproof found, and no erratum: the corrected modulus 8d-3 and the extra
exception d = 4 appear nowhere. Time-sensitive: an automated project opened an issue on this exact
conjecture yesterday, and the authors say they could write a proof.

---

## 3. arXiv:2502.21050 (Ying Wang, Yingrui Zhang, Hankel determinants of Motzkin convolution powers)

**Latest version.** v1, 28 Feb 2025. No journal reference or DOI in the arXiv record; Crossref has no
published version under this title. Conjecture 3 as printed (arXiv HTML):

> For r = 0 (mod 3): H_{rn}(F(x,r)) = H_{rn+1}(F(x,r)) = (-1)^n (n+1)^{r-1}, and
> H_{rn+2} + H_{rn-1} = 3^2 α (n+1)^{r-1} with |α| = r(r-3)/18.

**Searched.**
* Semantic Scholar citations: 3 papers (below).
* arXiv API: `Hankel AND Motzkin AND (convolution | "shifted periodic")`.
* WebSearch: title + `journal OR conjecture OR proof OR "shifted periodic"`.
* Crossref: bibliographic query for the title.
* GitHub issue search: `2502.21050`, `Motzkin Hankel convolution`, same inside trureturing.
* Zenodo API: `"2502.21050"`, `Motzkin AND Hankel AND convolution`.

**Found.**
* trureturing PR #6686, "math: refute the printed universal form of Conjecture 7, arXiv:2502.21050v1",
  merged 9 Sep 2026. https://github.com/the-omega-institute/trureturing/pull/6686
  Counterexample to Conjecture 7 (r = 1, 2 mod 3) as printed: r = 4, n = 0, `H_0 = 1`, `H_4 = -1`. The PR
  says the authors probably meant r >= 8, and that it "does not refute Conjecture 3 or Conjecture 8" and
  checked none of their values.
* Citing papers, none of which treats Conjecture 3:
  * F. Liu, Y. Wang, Y. Zhang, Z. Zhang, arXiv:2503.17187 (21 Mar 2025): extends Cigler's shifted-Hankel
    identities to power series satisfying a quadratic equation, with Motzkin numbers as an application.
  * H. Sun, arXiv:2611.27755 (20 Nov 2026; also Zenodo 24388865): Cigler's recurrence for Motzkin-path
    Hankel determinants via the Hirota-Miwa equation. Its conjecture references are to Cigler.
  * Y. Wang, Z. Zhang, arXiv:2409.18609: an earlier paper (weighted lattice paths).
* No paper by G. Xin and coauthors on these conjectures was found.

**Verdict.** No proof or disproof found, and no corrected statement. The sign error at r = 9, 12, 21, 24
(for example `H_9(M^9) = +256`) has not been noticed publicly. Only the printed form of Conjecture 7 has
been refuted (GitHub PR, 9 Sep 2026).

---

## 4. arXiv:2607.16608 (Bharadwaj, Sujatha, Chandankumar, 2-color overpartitions)

**Latest version.** v1, 18 Jul 2026. No later version, no journal reference.

**Searched.**
* Semantic Scholar citations: empty.
* arXiv API: `"2-color overpartition" | "two-color overpartition"`, `au:Chandankumar`, the 30 newest papers
  with "overpartition(s)" in the title (Nov 2025 to Nov 2026).
* WebSearch: `"2-color overpartition" ... Bharadwaj Sujatha Chandankumar conjecture proof`; a second query
  with `"4n+3" "modulo 8"` and Sellers / Saikia / Tang.
* Crossref: title query; `query.author=Sumit Kumar Jha`, Research Square prefix, since Aug 2026.
* GitHub issue search: `2607.16608`, `"2-color overpartition"`, `Bharadwaj Sujatha Chandankumar`,
  `repo:the-omega-institute/trureturing overpartition`.
* Zenodo API: `"2607.16608"`, `"2-color overpartition" OR "two-color overpartition"`.

**Found.**
* Nothing that cites the paper or states its closing conjectures.
* Near misses, all on other functions:
  * Anakha V, C. Ray, arXiv:2609.39697 (30 Sep 2026), colored overpartitions `C_{m,r,s}(n)`. Their mod-8
    family is for `C_{2^α,1,2}` (two colors on parts *not* divisible by m); `p̄_ℓ` would be `C_{ℓ,2,1}`. The
    paper does not cite 2607.16608 and has no `4n+3` congruence.
  * L. Zhang, Q.-M. Zhang, arXiv:2611.16145 (restricted overpartition function `t̄^{(3)}`).
  * J. Math. Sci. DOI 10.1007/s10958-026-08664-0 (Oct 2026), "2-color overpartition pairs".
  * Saikia, Talukdar, Tang, arXiv:2609.11806 (k-colored overpartitions into odd parts).
* trureturing has 14 overpartition issues; none is about this paper.
* Jha's only relevant Research Square preprint is the one already recorded in NOVELTY-1 (Wei's conjectures).

**Verdict.** No proof or disproof found for any of the five closing conjectures.

---

## 5. arXiv:2609.05358 (H. Manabe, purely periodic three-move subtraction games)

**Latest version.** v2, 16 Sep 2026 (38 pages, with Lean ancillary files). Our statement is
**Conjecture 56 (period restriction)** in v2 (Conjecture 54 in v1):

> Let S be non-additive with a >= 2 and gcd(a,b,c) = 1. If G_S is purely periodic then its least period
> satisfies p in {a+b, c+a, c+b}.

v2 also has Conjecture 52 (pure periodicity iff the angle `ρ = c mod (a+b)` is admissible; "Only necessity
remains open") and Conjecture 57 (admissible angle); 56 and 57 together imply 52 for `c >= 2(a+b)`.

**Searched.**
* Semantic Scholar citations: empty.
* arXiv API: `"subtraction games" | "subtraction game"` (12 newest).
* WebSearch: title + `conjecture least period "a+b" proof OR counterexample`.
* GitHub issue search: `2609.05358`, `subtraction Manabe`, inside trureturing and agentic-conjectures.
* Zenodo API: `"2609.05358"`, `"subtraction game" AND "three-move"`.
* README of `s66lem/conjectures`; manifest of `ArtificialZeng/AI-Has-Taste`.

**Found.**
* **S. Levinstein, "Admissible Angles for Purely Periodic Three-Move Subtraction Games"**,
  10.5281/zenodo.24536124, 25 Nov 2026, with Lean (`github.com/s66lem/conjectures`, result 15). Abstract:
  proves Conjecture 57 "in a strictly stronger form", removes the hypothesis `c >= 2(a+b)` "so Conjecture 52
  is equivalent to Conjecture 56 for every c > b", and establishes "two partial period exclusions toward
  Conjecture 56". Conjecture 56 itself is not claimed.
* Z. Zeng, "Pure Periodicity for the Three-Move Subtraction Games {2,5,c}", 10.5281/zenodo.22661435,
  8 Sep 2026: for `{2,5,c}`, purely periodic iff `c mod 7` in {1,2,3,5,6}, with least period `c+2`, 7 or
  `c+5`. So the conjecture holds in this one-parameter family. His archive also lists "Exact
  Sprague-Grundy Words for ... {2,5,c}" and "Necessity of admissible angles for the periods c+a and c+b"
  (Conjecture 57).
* M. Jo, 10.5281/zenodo.22669492 (9 Sep 2026): exact Grundy sequences for `{2,b,2b-6}`, b odd.
* trureturing issue #8384 (16-17 Sep 2026), exactly this conjecture: closed "not planned", "design
  abstention; no mathematical attempt admitted", "Conjecture 54 remains OPEN".
  https://github.com/the-omega-institute/trureturing/issues/8384
  Issue #13306 (5 Oct 2026, open) is about Conjecture 57 only.
* arXiv:2611.31336 (Y. Pan) and arXiv:2611.28065 (Zhao, Zhang, Shi): quadratic periods, but for four-move
  sets (a fixed triple plus one move). Not counterexamples to a three-move statement.
* Ward's conjecture (arXiv:1606.04029), which Conjecture 56 sharpens, is still a conjecture.

**Verdict.** No proof or counterexample of Conjecture 56 found. Partially settled: Levinstein (25 Nov 2026)
gives two partial period exclusions and the equivalence with Conjecture 52; Zeng (8 Sep 2026) the family
`{2,5,c}`. Time-sensitive: Levinstein's package appeared yesterday and stops one step short.

---

## 6. arXiv:2404.18517 (Chen, Kitaev, Zhang, separable permutations)

**Latest version.** v1, 29 Apr 2024. Published: Discrete Appl. Math. 355 (2024), 169-179,
DOI 10.1016/j.dam.2024.05.004. Conjecture 1: the four record statistics on separable permutations are
unimodal with peak at 2 (n >= 3). Conjecture 2: `rmax`/`lmin` on irreducible ones, peak at 3 (n >= 5).
Conjecture 3: `rmax`/`lmin` on reducible ones, peak at 1.

**Searched.**
* Semantic Scholar citations: 6 papers (Gao-Kitaev-Li-Li arXiv:2610.08865; Liao-Yang-Yu 2510.12046;
  Gao-Kitaev-Li-Ruan 2509.14341; Han-Kitaev-Zhang 2408.12865; two without identifiers). By title, none is
  about these conjectures.
* arXiv API: `"separable permutations"` (12 newest).
* WebSearch: `"separable permutations" "right-to-left maxima" unimodal conjecture Chen Kitaev Zhang proof`.
* GitHub issue search: `2404.18517`, `"separable permutations" unimodal`, `"separable permutations" Kitaev`,
  `repo:the-omega-institute/trureturing CKZ`.
* Zenodo API: `"2404.18517"`, `"separable permutations" AND unimodal*`.

**Found.**
* **Conjecture 2: `the-omega-institute/trureturing`.** PRs #12720, #12759, #12795, #12836, #12912
  (3-4 Oct 2026) and **PR #31245, "Formalize the full CKZ record peak-three theorem"**, merged
  17 Nov 2026. https://github.com/the-omega-institute/trureturing/pull/31245
  "Proves the source-v1 weak formulation of Chen-Kitaev-Zhang, arXiv:2404.18517v1 §3 Conjecture 2": for
  every n >= 5 the four distributions rise weakly below 3 and decline weakly from 3. Method: Newton-kernel
  positivity from a weighted quadratic; Lean module `RecordAllIndex.lean`. Tracking issue #12664 is open.
* **Conjectures 1 and 3: `1337isnot1337/solved-conjectures` issue #157**, "Chen-Kitaev-Zhang's record
  conjectures", opened 21 Nov 2026, closed as completed 24 Nov 2026.
  https://github.com/1337isnot1337/solved-conjectures/issues/157
  Result package: https://github.com/1337isnot1337/solved-conjectures/tree/main/results/chen-kitaev-zhang-records
  README (read directly): "The original Conjectures 1 and 3 ... hold in full." Method: Schröder-type
  generating-function relation, three low-record comparisons reduced to positive-coefficient inequalities
  with 16 Pell-type certificates, plus the trureturing tail theorem for Conjecture 2. "Lean status: complete
  for the stated Conjectures 1 and 3" (`CKZBridge.lean`; "the combined repository gate is pending"). "The
  notes remain unrefereed and priority is subject to correction."
* trureturing issue #20489 (22 Oct 2026, open) only preregisters Conjecture 1; no result there.

**Verdict.** Already settled, in unrefereed form: all three conjectures have public Lean-backed proofs
(Conjecture 2 by 17 Nov 2026, Conjectures 1 and 3 by 24 Nov 2026). Nothing on arXiv or in a journal.

---

## 7. arXiv:2611.05196 (Bardhan, Saikia, COS/CES in cubic partitions)

**Latest version.** v1, 3 Nov 2026. No later version.

**Searched.**
* Semantic Scholar citations: empty (after two 429s).
* arXiv API: `"cubic partitions" AND ("odd part" | "even part" | "sum of parts")`, `au:Bardhan_Gaurab`.
* WebSearch: `"cubic partitions" "odd part" "even part" sums divisibility Bardhan Saikia conjecture COS CES proof`.
* GitHub issue search: `2611.05196`, `"cubic partitions" Bardhan Saikia`, inside trureturing.
* Zenodo API: `"2611.05196"`, `"cubic partitions" AND (COS OR "odd part")`.

**Found.**
* trureturing issue #26479, "Bardhan-Saikia full COS/CES prime-power conjecture (arXiv:2611.05196v1,
  Section 4 Eq. 4.1)", opened and closed "not planned" on 4 Nov 2026.
  https://github.com/the-omega-institute/trureturing/issues/26479
  "no Lean proof or refutation of the complete joint conjecture was produced". A numerical check of 3,748
  parameter pairs with m <= 20,000 found no counterexample. The α = 1 case is already in the source paper.
* No follow-up by Bardhan and Saikia: their later arXiv paper of the same day (2611.05160) is on another
  topic. No other paper on COS/CES found.

**Verdict.** No proof or disproof found for α >= 2.

---

## 8. arXiv:2611.05623 (Yaqubi, Mirzavaziri), re-check of Conjectures 6.3 and 6.4

**Latest version.** v1, 3 Nov 2026 (unchanged).

**Searched.** trureturing issues #26878, #26887, #30694 re-read; GitHub issue search `2611.05623`,
`"designated summands"`, `repo:the-omega-institute/trureturing Yaqubi`; `facebookresearch/atlas-lean` PR #8061;
Zenodo API `"designated summands"`, `"2611.05623"`; arXiv API `"designated summands"`; README of
`s66lem/conjectures`.

**Found.**
* #26878 (Conj. 6.3): still open. Last comment 9 Nov 2026 (a prior-art note: another repository,
  "aftd-platform", proves the instances n = 0, 1 and the sharpness witness, "so Conjecture 6.3 is not
  settled there"). "the required coefficient vanishing has not yet been proved." No linked PR.
  https://github.com/the-omega-institute/trureturing/issues/26878
* #26887 (all-index eta-quotient congruence for 6.3): open, all activity 5 Nov 2026. "The exact event remains
  unresolved." PR #27019 merged a classical partial theta identity only.
* #30694 (15 Nov 2026): the probe timed out; "remains unresolved".
* New since NOVELTY-1: trureturing issue #34438 (26 Nov 2026) preregisters Conjecture **6.2**; issues #33469
  and #33480 concern 6.1. Levinstein's Zenodo record 10.5281/zenodo.24552279 (26 Nov 2026) is 6.2 only.
* `facebookresearch/atlas-lean` PR #8061 (25 Nov 2026): its only designated-summands item is Sellers'
  mod-4 congruences for `PD_2`; nothing from this paper.
* Conjecture 6.4 is mentioned nowhere. No use of Radu's algorithm on this paper was found.
* arXiv: no new designated-summands paper since arXiv:2611.28840 (21 Nov 2026).

**Verdict.** No proof or disproof found for Conjectures 6.3 and 6.4, as of 26 Nov 2026. The paper is still
being worked on by automated groups (6.1 proved 23 Nov, 6.2 on 26 Nov), so priority remains time-sensitive.

---

## Automated-conjecture projects checked (all items)

| Project | How checked | Relevant hits |
|---|---|---|
| `the-omega-institute/trureturing` (issues, PRs) | GitHub issue search by arXiv id and keyword, repo-scoped; pages via WebFetch | items 1, 2 (other conjectures), 3 (Conj. 7), 5, 6, 7, 8 |
| `umaia1234/agentic-conjectures` | issue search, README | item 1 (PRs #73, #76) |
| `1337isnot1337/solved-conjectures` | issue search, README, CONTENTS.md | item 6 (issue #157) |
| `ArtificialZeng/AI-Has-Taste` | README and `archive-manifest.json`, keyword grep | items 1 and 5 |
| `s66lem/conjectures` (S. Levinstein) | README (39 results) | item 5 (Conj. 57), item 8 (Conj. 6.2) |
| `ndesmo/leaning` | issue #24, PROBLEMS.md | item 2 (open target) |
| `facebookresearch/atlas-lean` | issue search, PR #8061 | none |
| `openai/math` | README, CONTENTS.md keyword grep | none |
| whataifound.org | `/data/entries.json` (290 entries), keyword grep | none |
| Zenodo | API, by arXiv id and keyword | items 1, 5, 8 |
| Research Square (S. K. Jha) | Crossref author query | none new |

Not reachable: GitHub REST API for repository contents (403); Google Scholar was not tried.
