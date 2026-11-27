# NOVELTY-1: prior-art search for six settled conjectures

Search date: 2026-11-26. Tools: WebSearch, WebFetch, `curl` (arXiv API, arXiv abs/html pages,
Semantic Scholar API, Crossref, OpenAlex, Zenodo API, OEIS, GitHub public issue-search API).
No papers were stored in the repository.

## Summary table

| # | Source (latest version) | Verdict |
|---|---|---|
| 1 | arXiv:2603.08510 **v1**, 9 Mar 2026 | **Already settled (claimed, unrefereed)**: all four families, by S. K. Jha, Research Square preprint, posted 26 Oct 2026 |
| 2 | arXiv:2611.05623 **v1**, 3 Nov 2026 | **No proof/disproof found** for Conj. 6.3 and 6.4. (Conj. 6.1 and 6.2, not on our list, have claimed proofs.) |
| 3 | arXiv:2605.12822 **v1**, 12 May 2026 | **Already noticed**: same counterexample `(1+q)^6(1+q^3)` in a public GitHub PR, 10 Sep 2026. Not in any arXiv version or citing paper. |
| 4 | arXiv:2601.16626 **v1**, 23 Jan 2026 (= J. Integer Seq. 28 (2025), Art. 25.7.1) | **"If" direction already public** (same row relation, GitHub issue, 20 Sep 2026, unrefereed). "Only if": no proof found anywhere. |
| 5 | arXiv:2506.03620 **v1**, 4 Jun 2025 (publ. Front. Combin. Number Theory 1 (2026), 16-25) | Lower bound 13/66: **no proof/disproof found**. Upper bound 1/2: **partially settled (3 <= a <= 8) by Chen-Zhu's companion paper** arXiv:2506.03625. |
| 6 | arXiv:2503.06228 **v1**, 8 Mar 2025 (publ. Acta Math. Hungar., DOI 10.1007/s10474-025-01573-9) | **Already settled by the source paper itself**: its Theorem 1.2 proves its Conjecture 1.1. |

Version check for all six: `https://export.arxiv.org/api/query?id_list=...` and the
"Submission history" block of each `https://arxiv.org/abs/ID` page. Every paper has only v1.

---

## 1. arXiv:2603.08510 (X. Wei, overpartitions mod 7, 17, 19, 23)

**Latest version.** v1, 9 Mar 2026. No later version; the four statements are still conjectural there
(two unnumbered in "Closing remarks", Conjectures 4.1 and 4.2).

**Searched.**
* Semantic Scholar citations of arXiv:2603.08510: empty list (HTTP 200 after two 429s).
* arXiv API: `all:overpartition AND all:congruences` (40 newest), `all:overpartitions AND (modulo 7|11|13|17|19|23)`,
  `au:Wei_XuanLing`, `all:overpartition AND all:Hecke AND all:congruences`,
  `all:overpartition AND all:"half-integral weight"`, `au:Jha_Sumit`.
* WebSearch: title of the paper + "conjecture proof"; `overpartition "modulo 7" Ramanujan-type "(mod 7)"`.
* Zenodo API (`overpartition AND congruence ...`, `"2603.08510"`): nothing.
* GitHub issue search (`2603.08510`, `overpartition Wei`): only an arXiv-listing bot.
* Older literature read in full text (arXiv HTML): arXiv:2309.01792v4, arXiv:1911.05799v3, arXiv:2303.15895v2.

**Found.**
* **Sumit Kumar Jha, "On Four Conjectural Congruence Families for the Overpartition Function"**,
  Research Square preprint, DOI 10.21203/rs.3.rs-11336379/v1,
  https://www.researchsquare.com/article/rs-11336379/v1. Manuscript dated 11 Oct 2026, posted
  26 Oct 2026, 14 pages, status "under review" at The Ramanujan Journal (reviewer agreed 28 Oct 2026).
  Not on arXiv. States an "AI-assisted preparation" note. Abstract: "We establish four Ramanujan-type
  congruence families for the overpartition function proposed by Xuanling Wei. Two are the numbered
  conjectures modulo 7 and 19; the other two are the proposed square-prime families modulo 17 and 23."
  Its theorems match ours exactly:
  * Thm 1.1 = Wei Conj. 4.1: `p(896n+176) = p(896n+688) = p(896n+816) = 0 (mod 7)` (Radu-Sellers criterion).
  * Thm 1.2 (mod 17): `p(17m) = -2H(7,m) (mod 17)` for `m = 3 (mod 8)` (Cohen Eisenstein series), plus a
    Hecke relation for every odd prime Q; Q = 11 gives `p(17*11^2 m) = 0` for `(m/11) = -1`.
  * Thm 1.3 = Wei Conj. 4.2 (mod 19), with explicit Hecke identities at 13^2 and 17^2.
  * Thm 1.4 (mod 23): Hecke relation for every odd Q; Q = 13 gives `p(23*13^2 m) = 0` for `(m/13) = 1`.
  It goes beyond the four statements (relations for all odd primes Q mod 17 and 23, prime-power towers).
  I read the first pages and the reference list only; I did not check its proofs.
* Older literature does **not** contain these congruences:
  * Ryan, Sirolli, Villegas-Morales, Zheng, "Explicit families of congruences for the overpartition
    function", Ramanujan J. 65 (2024), arXiv:2309.01792. Thm 1.5 gives `p(m l^2 n) = 0 (mod m)` for
    **all** n prime to l with a Legendre condition, for m = 17 only with
    l in {167, 541, 911, 1013, 1153, 1867, 1931, 2543, 2683, 2887, 3019, 3023, 3329, 4243, 4651}, for
    m = 19 only l = 2207, and nothing for m = 23. l = 11 (mod 17) and l = 13 (mod 23) are absent; Wei's
    families need the extra restriction n = 3 or 5 (mod 8). For m = 7 their families are `p(7 l^3 n)`,
    `p(7 l^2 n)`, a different shape from `p(16n)`. They state: "the results in this article give the
    first known congruences for m > 7."
  * Ryan, Scherr, Sirolli, Treneer, "Congruences satisfied by eta-quotients" (arXiv:1911.05799) and
    Barquero-Sanchez et al. (arXiv:2303.15895): only `p(Q^3 n) = 0 (mod l^j)`, l in {3,5,7}, for large
    primes `Q = -1 (mod 16 l^j)`.
  * Treneer 2006 is existence-only. No paper found with a mod-7 congruence of the form `p(16(56n+k))`.

**Verdict.** Already settled (as a claim in an unrefereed preprint) by Jha, 26 Oct 2026, one month
before our write-up. Our proofs are an independent confirmation, not a first proof. Not contained in
the pre-2026 literature.

---

## 2. arXiv:2611.05623 (Yaqubi, Mirzavaziri, k-tuples with designated summands)

**Latest version.** v1, 3 Nov 2026. In the paper our statements are Conjecture 6.3
(`A_6(12n+10) = 0 (mod 3^5)`, and 3^5 exact) and Conjecture 6.4(a)-(c). The authors write that 6.3 is
"a fact we have not yet been able to prove".

**Searched.**
* Semantic Scholar citations: empty.
* arXiv API: `all:"designated summands"` (13 hits, all read by title; abstracts of arXiv:2611.15164 and
  arXiv:2611.28840 read), `au:Yaqubi`.
* WebSearch: `"designated summands" k-tuples congruences Yaqubi Mirzavaziri conjecture modulo 243 proof`;
  `"bipartitions with designated summands" OR "tripartitions ..." congruences modulo 9 27 81`;
  `"Conjecture 6.4" OR "Conjecture 6.3" ... "A_6(12n+10)" OR "A_3(12n+7)" OR "A_8(18n+15)"`;
  `"partition triples with designated summands" OR "PD_{-3}" ...`.
* Zenodo API `"designated summands"`; GitHub issue search `Yaqubi`, `2611.05623`.

**Found.**
* arXiv:2611.15164v2 (Li, Tang, 11/24 Nov 2026) and arXiv:2611.28840 (Saikia, 21 Nov 2026): both about
  Lin's function `PDO_t(n)` mod powers of 3 and 2. Different function; they do not cite the A_k conjectures.
* GitHub `the-omega-institute/trureturing` (public, AI-driven conjecture tracker):
  * issues #26878 and #26887 (5 Nov 2026) and #30694 (15 Nov 2026) target Conjecture 6.3. Status
    there: "The exact event remains unresolved"; only numerical checks for n <= 250 and the witness
    `A_6(22) = 486 (mod 729)`. https://github.com/the-omega-institute/trureturing/issues/26878
  * PR #33580 (merged 23 Nov 2026): Lean proof of **Conjecture 6.1** (2-adic, `A_k(4n+3)`). Not on our list.
* Zenodo 10.5281/zenodo.24552279 (S. Levinstein, 26 Nov 2026): claimed Lean-verified proof of
  **Conjecture 6.2** (2-adic, `B_k(4n+3)`). Not on our list. Its repository `github.com/s66lem/conjectures`
  lists only 6.2 from this paper.
* Nothing at all mentions Conjecture 6.4.
* Earlier literature covers only k = 1 and k = 2: Hao-Shen, "Congruences modulo 9 for bipartitions
  with designated summands" (Turkish J. Math. 42 (2018)); Naika-Shivashankar (Bol. Soc. Mat. Mex. 24
  (2018)); Xia (JNT 159 (2016)) for PD(n) mod 9, 27. Naika-Harishkumar (Palestine J. Math. 2022) treat
  l-regular triples, a different generating function. No paper found with `(f_6/(f_1 f_2 f_3))^k` for
  k in {3,4,6,7,8} on these progressions.

**Verdict.** No proof/disproof found for Conjectures 6.3 and 6.4. This paper is being actively worked
on by automated groups (6.1 and 6.2 fell on 23 and 26 Nov 2026; 6.3 has open public issues), so
priority is time-sensitive.

---

## 3. arXiv:2605.12822 (Connelly, Ito, Martinez, Shevchenko, Yang, q-Fibonomial unimodality)

**Latest version.** v1, 12 May 2026. Conjecture 5.4 still contains "Moreover, if k <= 3 or r <= 3, this
condition is also necessary". The authors tested k <= 5, r <= 6, max{a_i,b} <= 15.

**Searched.**
* Semantic Scholar citations: empty.
* arXiv API: `all:Fibonomial`, `all:unimodal AND all:"q-integers"`, `all:unimodal AND all:Fibonomial`.
* WebSearch: title + `conjecture counterexample unimodal product of q-integers "[b]_{q^r}"`.
* Zenodo API: `Fibonomial AND unimodal*`, `"q-integers" AND unimodal*`, `"2605.12822"`.
* GitHub issue search: `2605.12822`, `Fibonomial`, `Connelly Ito Martinez`.

**Found.**
* **GitHub `the-omega-institute/trureturing` PR #6805, "feat: refute Conjecture 5.4 necessity at r=3
  and k=6"**, opened and merged 10 Sep 2026.
  https://github.com/the-omega-institute/trureturing/pull/6805
  Same counterexample: `[2]_q^6 [2]_{q^3} = (1+q)^6(1+q^3)`, coefficients 1, 6, 15, 21, 21, 21, 21, 15, 6, 1,
  k = 6, a_i = 2, b = 2, r = 3. Formalized in Lean (`QProductNecessityRefutation`). The PR refutes only
  the `r <= 3` branch; it says the `k <= 3` branch is not refuted and sufficiency is not addressed. It
  makes "no claim of first discovery" and does not say the authors were told.
* Zenodo preprints on the main Bergeron-Ceballos-Kustner conjecture (not on Conjecture 5.4, judged from
  their abstracts; PDFs not read): hainzhao (n = 4, 6 Aug 2026, 10.5281/zenodo.21826970); K. Fathi
  (columns <= 14, 30 Aug 2026, zenodo.22177919; full conjecture claimed, 14 Oct 2026, zenodo.23356229);
  M. Stefan-Donovan (n = 4 with Lean, 12 Sep 2026, zenodo.22719863, cites 2605.12822); S. Jiang (full
  conjecture claimed, 9 Oct 2026, zenodo.23260785); Anonymous (25 Nov 2026, zenodo.24541015).

**Verdict.** The failure of the necessity clause at r = 3, k = 6 has already been noticed publicly
(GitHub PR, 10 Sep 2026, same polynomial). It has not been noticed in a later arXiv version or in any
citing paper. Our larger families of counterexamples (k = 6, 7) were not found elsewhere.

---

## 4. arXiv:2601.16626 (Merikoski, Haukkanen, Sasaki, Tossavainen, LCM/GCD generalized eigenvalue -1)

**Latest version.** v1, 23 Jan 2026. Published as J. Integer Seq. 28 (2025), Article 25.7.1 (received
5 May 2025, published 2 Dec 2025), https://cs.uwaterloo.ca/journals/JIS/VOL28/Tossavainen/tossa12.html.
Conjecture 5 is stated in both; the abstract still calls it a conjecture. OEIS A004754 links the paper
and carries no comment about a proof.

**Searched.**
* Semantic Scholar citations: empty.
* arXiv API: `all:"LCM matrix"`, `all:"LCM matrices"`, `all:"generalized eigenvalue" AND all:GCD`.
* WebSearch: `"generalized eigenvalues" "LCM matrices" "GCD matrices" ... binary expansion begins with 10 proof`;
  `matrix lcm(i,j)+gcd(i,j) singular ... A004754`.
* Zenodo API (`"LCM matrix"`, `A004754`, `"2601.16626"`); GitHub issue search (`2601.16626`,
  `Merikoski Haukkanen Sasaki Tossavainen`, `A004754 LCM GCD`); OEIS A004754 JSON.

**Found.**
* **GitHub `the-omega-institute/trureturing` issue #9180**, "Research Conjecture 5: binary
  classification of the LCM-GCD eigenvalue -1", opened 20 Sep 2026, still open.
  https://github.com/the-omega-institute/trureturing/issues/9180
  * "If" direction: "The binary-10 direction has a uniform paper-level proof", by the kernel vector
    `x_r = -2, x_{2r} = 3, x_{4r} = -1` (r = 2^s) for `4r <= n < 6r`. This is our row relation
    `R_{4r} - 3R_{2r} + 2R_r = 0` up to sign. Not Lean-checked, not refereed.
  * "Only if": reduced to non-vanishing of a Schur-complement determinant `det T_n`. "The full
    conjecture remains OPEN."
* `facebookresearch/atlas-lean` PR #7724 (25 Nov 2026) formalizes only Theorem 6 of the paper (Cauchy
  interlacing), not the conjecture.
* arXiv:2609.07352, arXiv:2607.01265, Zenodo 21212231: LCM-matrix papers on other questions.

**Verdict.** "If" direction: already in the public record (unrefereed GitHub issue, 20 Sep 2026) with
the same explicit relation; no journal or arXiv proof found. "Only if" direction: no proof found.

---

## 5. arXiv:2506.03620 (Chen, Zhu, primes of the form ax+by)

**Latest version.** v1, 4 Jun 2025. Published: Front. Combin. Number Theory 1 (2026), 16-25,
DOI 10.3934/fcnt.2026002 (journal page not checked). Conjecture 1.4: `pi(a,b) >= (13/66) pi(s)`,
equality iff (a,b) = (3,166). Conjecture 1.5: `pi(a,b) <= (1/2) pi(s)`, equality iff (a,b) = (3,5).

**Searched.**
* Semantic Scholar citations of arXiv:2506.03620 and of arXiv:2506.03625: the same three papers (below).
* arXiv API: `all:"primes of the form" AND all:Frobenius`, `all:"ax+by" AND all:primes AND cat:math.NT`,
  `au:"Honghu Liu"`, `all:"Ramirez Alfonsin" AND all:primes`, `all:"numerical semigroup" AND all:primes AND all:Frobenius`.
* WebSearch: `Chen Zhu "Primes of the form ax+by" conjecture "13/66"`; `"Note on primes of the form ax + by"
  Lithuanian Mathematical Journal 2026`; `... "Chen and Zhu" conjecture "13/66" proved OR disproved`.
* Crossref for DOI 10.1007/s10986-026-09723-4 (reference list); OpenAlex and Semantic Scholar for its
  abstract; Zenodo API; GitHub issue search.
* Full texts read: arXiv:2506.03620, arXiv:2506.03625, arXiv:2510.01781, arXiv:2509.08599.

**Found.**
* **Y.-G. Chen, H. Zhu, "The number of primes not in a numerical semigroup"**, arXiv:2506.03625 (same
  day; Int. J. Number Theory 22 (2026), 307-317). `pi*(a,b)` = number of primes not of the form
  `au+bv`. Since every such prime is `<= s(a,b)`, `pi*(a,b) = pi(s) - pi(a,b)`. Its Conjecture 1.4,
  `pi*(a,b) >= (1/2) pi(s)` with equality iff a = 1, (2,3), (2,5), (3,5), is therefore the same as
  the upper bound `pi(a,b) <= (1/2) pi(s)` for a >= 3. **Theorem 1.5 there: "Conjecture 1.3 and
  Conjecture 1.4 are both true for a <= 8."** So the upper bound is proved for 3 <= a <= 8, all b.
  The source paper does not point this out.
* Asymptotics consistent with both bounds: Chen-Zhu Thm 1.1 (fixed a, b -> infinity:
  `pi(a,b) = (1/2 - 1/(2(a-1)) + o(1)) pi(s)`, so 1/4 for a = 3); Ding-Zhai-Zhao (`~ 1/2` as a -> infinity).
* Citing papers, none of which settles either bound:
  * Ding, Komatsu, Liu, arXiv:2510.01781 (RACSAM, DOI 10.1007/s13398-026-01888-4): l-Frobenius
    generalization; quotes the `pi* >= 1/2 pi` conjecture as open.
  * Ding, Wang, Zhang, arXiv:2509.08599: coprime-coefficient variant.
  * H. Liu, "Note on primes of the form ax + by", Lith. Math. J. 66 (2026), 414-423,
    DOI 10.1007/s10986-026-09723-4. Abstract **not reachable** (Springer cookie wall; publisher withholds
    it from Semantic Scholar and OpenAlex). A search snippet describes it as a GRH-conditional result on
    primes with exactly l+1 representations.
* No public tracker issue on this paper.

**Verdict.** Lower bound (13/66): no proof/disproof found. Upper bound (1/2): partially settled, for
3 <= a <= 8, by Chen-Zhu arXiv:2506.03625 Theorem 1.5; open in general. Our numerical verification
(a <= 400, b <= 3000) is evidence only.

---

## 6. arXiv:2503.06228 (Nath, Sarma, ped(n) mod 192)

**Latest version.** v1, 8 Mar 2025. Published in Acta Math. Hungar., DOI 10.1007/s10474-025-01573-9.

**Searched.**
* Full text (arXiv HTML) of arXiv:2503.06228 and of arXiv:2403.16735v2.
* Semantic Scholar citations of arXiv:2503.06228 (empty) and of arXiv:2403.16735 (two, both by Nath et al.).
* arXiv API: `all:"PED partitions"`, `abs:ped AND abs:congruences AND cat:math.NT`,
  `all:partitions AND all:"modulo 192"`, `all:"even parts distinct"`, `all:"distinct even parts"`,
  `all:Nath AND all:Sarma AND all:conjecture`.
* WebSearch: `Nath Sarma ... conjecture ped(225n+43) modulo 192 proof`; `... "modulo 192" "9t^2" ... Conjecture 7.1`.
* Zenodo API, GitHub issue search.

**Found.**
* The conjecture `ped(225n+43), ped(225n+88), ped(225n+133), ped(225n+223) = 0 (mod 192)` was posed by
  H. Nath, "New congruences for partitions where the even parts are distinct", arXiv:2403.16735,
  INTEGERS 24 (2024), #A79, Conjecture 3.1.
* arXiv:2503.06228 restates it as "Conjecture 1.1" and then proves it: **"Theorem 1.2. Conjecture 1.1
  is true."** Abstract: "Using an algorithm developed by Radu, we prove congruences modulo 192 which
  were conjectured by Nath". Section 3 is the proof: reduce to `a(25n+r) = 0 (mod 16)`, r in
  {4, 9, 14, 24}, for `f_2^4 f_3^6 f_4 / f_1^11`; Radu data (m, M, N) = (25, 12, 60), bound 53; checked
  in Mathematica.
* The statement the paper leaves open is **Conjecture 7.1**: for a prime p >= 5 with (-2/p) = -1,
  (t,6) = 1, p | t and 1 <= j <= p-1,
  `ped(9 t^2 n + 9 t^2 j/p + (57 t^2 - 1)/8) = 0 (mod 192)`. The case t = p = 5 is exactly Conjecture 1.1.
  No proof of Conjecture 7.1 found.
* arXiv:2510.22006 (Dockery, Jameson) and arXiv:2503.12583 (Paudel, Sellers, Wang): unrelated to mod 192.

**Verdict.** Already settled by Nath and Sarma in the source paper itself (Theorem 1.2, March 2025).
Our Radu certificate is a re-proof. A new target in this paper would be Conjecture 7.1.

---

## Sources not reachable or not consulted

* Springer page of H. Liu's Lith. Math. J. note (abstract): not reachable.
* GitHub REST API for repository contents (issue bodies, files): refused by the session proxy (HTTP 403,
  would need repository attachment). Issue and PR pages were read through WebFetch summaries instead, so
  quotations from them are second-hand.
* MathSciNet, zbMATH, Google Scholar: not consulted (no access with the allowed tools).
* PDFs of the six Zenodo q-Fibonomial records: only abstracts read.
* Jha's Research Square PDF: first two pages and references read; proofs not checked.
