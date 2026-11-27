# Prior-art check for results (A)-(D), 2026-11-24

Method: web search plus reading primary text where it could be fetched (arXiv abstract pages, arXiv API listings, and pdftotext of Washington II, Boran et al., Platt-Trudgian, and the Hamieh-Kadiri-Martin-Ng list). "Read" = I read the text. "Abstract only" = I saw only the abstract or a search snippet. Coverage is limited; "not found" means not found by these searches, not a proof of absence.

Not read in full: Rubinstein-Sarnak 1994, Hardy-Littlewood 1918, Landau, Knapowski-Turan, Kaczorowski's papers, Fujii 1988, the Ford-Konyagin survey 1910.08983, the second half of Bailleul-Hayani-Untrau, Shimada-Koyama (MDPI returned 403; I saw only the Preprints.org and IDEAS summaries).

## Verdicts

| Q | Verdict |
|---|---|
| 1. Monotonicity of a race density in a continuous weight; coupling, peakedness, Anderson | Not found. (B) looks novel. |
| 2. Gerard-Washington and Washington II | A neighbouring problem that must be cited. It does not contradict the claim, which has to be phrased carefully. |
| 3. Smoothed or Riesz-mean results valid for all large x | The method of (A) is anticipated: Fujii/Platt-Trudgian, Johnston, Zhao and Washington II Thm 3 are all "positive bias beats an absolutely convergent zero sum, for all large x". The statement of (A) (nesting across weights) and the identity D_m - x^{-κ}D_{m'} = Riesz mean were not found. |
| 4. 2025-26 papers | Nothing found that overlaps (A)-(D). Nothing states b-independence for general real b. |
| 5. Unimodality of Rubinstein-Sarnak limits | Not found. |
| 6. User's own arXiv submission | Not visible as of 2026-11-24. |

## Q1. Monotonicity in a continuous weight; coupling; peakedness; Anderson

Not found, from several angles:
- Searches for "monotone density", "weight", "prime race", "Anderson's inequality", "peakedness", "coupling" returned only the standard Rubinstein-Sarnak/Devin/Lamzouri/Fiorilli material.
- Bailleul-Hayani-Untrau, arXiv:2603.20093 (read the first ~100k of 160k characters). "Coupling" appears only as the definition of a transport plan in the Wasserstein setup. There is no use of Anderson, no peakedness, no unimodality, and no monotonicity in a parameter. The paper builds general maps t_q and studies convergence rates. I did not read the last 60k characters (Sections 5.4-5.5).
- Closest "parameter dependence" results (none is monotonicity in a continuous weight):
  - Fiorilli, "Highly biased prime number races", arXiv:1210.6946: densities can be pushed arbitrarily close to 1, which is dependence on the modulus.
  - Meng / Ford-Konyagin-type results on products of k primes, arXiv:1606.04877 (abstract only): the bias oscillates with the parity of k. That is non-monotone in a discrete parameter.
  - Dihedral/quaternion papers (arXiv:2001.06671, abstract only): densities tend to 1/2 as a family parameter grows.
  - Rubinstein-Sarnak's own monotonicity in q (bias dissolves as q grows with r fixed) is a limit statement, not a monotonicity statement.
- Harper-Lamzouri, arXiv:1509.07188: normal comparison and orderings of weakly correlated variables for many contestants. This is a comparison technique but not a monotone-in-weight or coupling argument.
- No result was found that varies a continuous parameter inside the Rubinstein-Sarnak random variable and proves a density ordering.

Verdict: (B) appears novel. Caveat: I could not read Rubinstein-Sarnak or Devin in full. A proof by sample-path coupling of a one-parameter family of explicit sums is easy to miss in a keyword search.

## Q2. Gerard-Washington and Washington II

Read: Washington II (arXiv:2209.12845, full PDF text) and Boran-Byun-Li-Miller-Reyes (arXiv:2309.16007, Theorems 1.1 and 1.2 and Remark 1). Gerard-Washington itself: abstract only, via Washington II's description.

**Gerard-Washington, "Sums of powers of primes", Ramanujan J. 45 (2018) 171-180.** Described in Washington II: it showed π_k(x) ~ π(x^{k+1}), proved an Ω± result, and gave a heuristic for why π_k(x) - π(x^{k+1}) is usually negative when k>0 and usually positive when -1<k<0. Washington II says the Ω± proof was based on an incorrect formula and is valid only if RH is false.

**Washington, "Sums of powers of primes II", arXiv:2209.12845, published Ramanujan J. 65 (2024) 783-795.** Abstract: "For a real number k, define π_k(x) = Σ_{p≤x} p^k. When k > 0, we prove that π_k(x) − π(x^{k+1}) = Ω±( x^{1/2+k} / log x · log log log x ) as x → ∞, and we prove a similar result when −1 < k < 0. This strengthens a result in a paper by J. Gerard and the author and it corrects a flaw in a proof in that paper. We also quantify the observation from that paper that π_k(x) − π(x^{k+1}) is usually negative when k > 0 and usually positive when −1 < k < 0."

The theorems read in the text:
- Thm 1 (RH): Ω± with x^{k+1/2}·logloglog x/log x for k>0, and x^{(k+1)/2}·logloglog x/log x for k<0.
- Thm 2: ∫_1^∞ (π_k(t) − π(t^{k+1}))/t^{k+2} dt = −log(k+1)/(k+1), which is negative for k>0 and positive for -1<k<0.
- Thm 3: for 0<k≤10.32, RH is equivalent to ∫_1^x (π_k(t) − π(t^{k+1})) dt < 0 for all sufficiently large x. This is an equivalence with RH for an integrated (Riesz-type) mean valid for all large x. Washington says the methods of Sections 2-3 were "inspired by [5]" = Johnston (see Q3).
- There is no logarithmic density, no Rubinstein-Sarnak distribution, no limiting distribution, and no dependence of a density on k. "Usually negative" is made precise only through Thm 2 (an average) and Thm 3 (an integrated sign).

**Boran-Byun-Li-Miller-Reyes, "Sums of powers of primes in arithmetic progression", arXiv:2309.16007, PUMP J. Undergrad. Res. 7 (2024) 29-50.** Thm 1.1: π_k(x;m,n) − π(x^{k+1};m,n) is O(x^{k+1} exp(−c√log x)) for k>-1. Thm 1.2: ∫_1^∞ (π_k(t;m,n) − π(t^{k+1};m,n))/t^{k+2} dt = −log(k+1)/((k+1)φ(m)). Remark 1 cites Chebyshev and Rubinstein-Sarnak and says that a finer analysis "necessitates ... the Riemann Hypothesis". The paper has numerics for x in [10^4, 10^8]. It compares the k-th power sum in a class with the unweighted count in the same class at x^{k+1}. It does not race two residue classes against each other and gives no densities.

**Consequence for the project.** The statement "the growing-weight regime m>0 appears untouched for races in arithmetic progressions" is not contradicted. Washington and Boran et al. never compare two classes with the same weight and never compute a Rubinstein-Sarnak density. They are a neighbouring problem (single weighted sum versus π(x^{k+1})) and must be cited. The paper should use language like: "sums of k-th powers of primes, k>0, have been studied as approximations to π(x^{k+1}) [GW, Wa, BBLMR], but the race between residue classes weighted by p^k, and its logarithmic density, has not been treated." Do not write "untouched" without the qualification.

There is also a useful structural link. The race D_k for a character χ is the χ-twisted analogue of π_k(x) − π(x^{k+1}), and the identity relating weights, integrated against t^{-k-1}, is the same Stieltjes manipulation Washington uses.

## Q3. Smoothed, Riesz and Cesaro versions valid for all large x

Closest anticipations, each read at the level stated:

1. **Hardy-Littlewood 1918 / Landau.** Per Platt-Trudgian's introduction (read): Chebyshev's conjecture Σ_{p>2} (−1)^{(p−1)/2} e^{−xp} → −∞ as x→0 is equivalent to GRH for χ_4. This is an Abel-mean "for all small x" statement.
2. **Fujii ("Some generalizations of Chebyshev's conjecture", Proc. Japan Acad. Ser. A Math. Sci. 64 (1988), no. 7, 260-263; corrected 2026-11-24 from "1991", following Platt-Trudgian's reference list) and Platt-Trudgian, "Fujii's development on Chebyshev's conjecture", arXiv:1802.00089, IJNT 15 (2019) 639-644 (read).** Fujii: for 0<α<4.19, Σ (−1)^{(p−1)/2} e^{−(xp)^α} → −∞ as x→0 is equivalent to GRH for χ_4. Platt-Trudgian: the same for 0<α<20.40442. The proof is a bias term of size x^{-1/2}Γ(1/(2α)) against Σ_ρ |Γ(1/(2α)+iγ/α)/ρ| x^{-1/2}, which converges absolutely because of the Γ decay. Positive bias dominates an absolutely convergent zero sum, for all small x. This is the same mechanism as (A) with an exponential kernel instead of the power kernel t^{-κ-1}, and for a single weight.
3. **Pintz / Johnston, "On the average value of π(t) − li(t)", Canad. Math. Bull. 66 (2023) 185-195, arXiv:2201.06184 (abstract and search summary only).** RH is equivalent to ∫_2^x (π(t)−li(t)) dt < 0 for all x>2. Also, unconditionally, integrals of (π(t)−li(t))/t^2 and related weighted versions are negative for all x>2.
4. **Zhao, "On the mean values of the error terms in Mertens' theorems", Res. Number Theory 11 (2025) paper 62, arXiv:2411.18903 (abstract only).** RH is equivalent to ∫_2^X E_i(x) dx > 0 for all X>2 (i=1,2), with analogues for Dirichlet characters and APs.
5. **Washington II Thm 3 (read; see Q2).** This is the only item found that mixes a power weight p^k with such an integrated sign statement.
6. **Liu-Wang-Shi, "Two problems in comparative prime number theory: the five-way race modulo 11, and the sign of the integrated prime counting error", Zenodo 10.5281/zenodo.22208026, 31 Aug 2026 (abstract only).** They "classify the residue-class pairs for which the integrated prime counting error is eventually negative", answering a problem from Hamieh-Kadiri-Martin-Ng (arXiv:2407.03530). Methods use Turing-certified zero lists.
7. **Suzuki, "On variants of Chebyshev's conjecture", arXiv:2411.07436, Ramanujan J. 68 (2025) 95 (abstract only).** Sign constancy of weighted summatory functions of Λ implies RH/GRH, and the paper proves sign constancy on average. Per Hayani's paper, Suzuki uses log p·log(x/p) smoothing weights.
8. **Alamadhi-Planat-Solé, arXiv:1112.2398 (abstract only).** Claims the inequality π'(x;q,N) > π'(x;q,R), with π' = π − ψ/log x, holding for all x is equivalent to GRH for q. I have not checked the proof, and the claim looks strong. Treat with caution.
9. **Kaczorowski.** Search returned only descriptions: under GRH, {x: π(x;4,3) > π(x;4,1)} has no natural density and its upper natural density is below 1. No Abel/Cesaro statement attributed to him was found. Not read.
10. **Knapowski-Turan.** Nothing relevant beyond the citation lists. Not read.

What I did NOT find: (i) the specific identity D_m(x) − x^{−κ}D_{m'}(x) = κ∫ t^{−κ−1}D_{m'}(t)dt used as a Riesz mean in a prime race with an explicit formula; (ii) the statement that lead sets nest across weights, or that D_{m'}(x) ≥ 0 implies D_m(x) > 0 for large x; (iii) the weights p^m for m>0 in any of the above. Items 2-6 use the same positive-bias-versus-absolutely-convergent-sum idea, so the proof technique of (A) is standard. The paper should say so, and present (A) as a new application: the kernel is the weight shift, and the conclusion is monotone inheritance across weights.

## Q4. Recent (2025-2026) papers

Read or checked:
- **Shimada-Koyama, "Weighted Prime Number Theorem on Arithmetic Progressions with Refinements", Mathematics 13(16):2564 (Aug 2025).** From IDEAS/Preprints.org summaries only: the weight is p^{-w}, with w in [0,1); proves π_w(x;q,a) ~ Li(x^{1−w})/φ(q), extending Shimada's q=1 result. No densities, no m>0, no monotonicity visible. Not read in full.
- **Koyama, arXiv:2607.28931** (v1 31 Jul 2026, v2 29 Aug 2026; revises the April 2026 manuscript "A Hidden Hierarchy of Chebyshev's Bias and the Dominance of −1 (mod N)"). Abstract: Gaussian-mollified Weil explicit formula; under the Deep Riemann Hypothesis, a secondary bias governed by log L(1,χ); ranking such as 7>3>5>1 mod 8 with error (log x)/√x. A deterministic regularized approach, not a density in a weight. No overlap with (A)-(D) seen (abstract only).
- **Hayani, arXiv:2512.23302** (v1 29 Dec 2025, v3 25 Jun 2026). Under GRH only, natural density 1 for weighted counts with weight p^{-1/2} (the Aoki-Koyama weight); explicit error terms; says Sheth independently has similar results in the automorphic setting. Abstract only, plus the search report that it reduces to higher moments. This is the closest paper to the "GRH-only, for most x" flavour but concerns the single weight 1/√p, and natural not logarithmic density.
- **Hayani, arXiv:2607.23150** (25 Jul 2026), prime ideal races and 2p-th roots. Unrelated to weights.
- **Bailleul-Hayani, arXiv:2508.04087** (Aug 2025), several competitors in prime ideal races. No weights.
- **Bailleul-Hayani-Untrau, arXiv:2603.20093** (20 Mar 2026). Wasserstein rates for races; see Q1. No weight-dependence.
- **Bhattacharya-Martin-Simpson, arXiv:2507.13504, "Correlations of error terms for weighted prime counting functions"** (v1 17 Jul 2025; v2 6 May 2026, to appear in Acta Arith.). Weights named in the abstract: 1/p and log(1−1/p)^{-1}. Computes densities of joint sign patterns of pairs of error terms under RH+LI (example: ψ(x)−x and a Λ(n)/n-weighted sum share a sign on a set of log density about 0.9865). It does not discuss p^a(log p)^b and does not state b-independence for general b. Joint sign-agreement densities are the nearest published object to a "coupling of two weights", but on the zeta side and for different error terms.
- **Sheth, arXiv:2405.01512** (v1 2 May 2024; MPCPS 179 (2025) 331-349). GRH; Euler product at the centre; Chebyshev's bias in the Aoki-Koyama framework. No p^m weights.
- **Koyama-Sheth, arXiv:2509.04187** (Sep 2025), newforms; bias determined by the order of vanishing at the centre. No weights.
- **Lay, arXiv:1505.03589, "Sign changes in Mertens' first and second theorems"** and **Lamzouri, IJNT 12 (2016)**: Mertens weights (p^{-1}); already in the project's citations via Martin.
- **Mossinghoff-Trudgian, "Oscillations in weighted arithmetic sums", IJNT 17 (2021) 1697-1716, arXiv:2007.14537** (abstract only): sums involving Ω(n), ω(n); not prime races. Cite only as Martin did.
- **Devin, Lamzouri, Fiorilli:** no weighted-race paper from them turned up in 2025-26. Devin arXiv:1706.06394 and 1909.03975 give existence of logarithmic limiting distributions for general L-function sums; no b-dependence.

arXiv listing searches (export.arxiv.org API, run 2026-11-24):
- Query "weighted prime race" returned 0 results.
- "prime race" AND weight returned 0 results.
- A search over abstracts for "prime number race", "prime races", "Chebyshev's bias", "Chebyshev bias" returned 44 results. The newest are Koyama 2607.28931, Hayani 2607.23150, Bates et al. "Ties in Function Field Prime Races" 2603.21005, and Hayani 2512.23302. No paper on weights p^m with m>0, and no paper on b-independence.
- The "logarithmic density" + weight/monotone + prime query returned only Hayani 2512.23302, BMS 2507.13504, and one unrelated paper.

**b-independence.** Not found for general real b. Your own Proposition (b ∈ {0,1}) and the BMS abstract on Mertens weights remain the known cases. Note the mechanism: multiplying by (log p)^b changes the explicit-formula coefficients by factors (log x)^b times a term O(1/log x), which suggests the claim for all real b is expected, but I found no source stating it.

## Q5. Unimodality / log-concavity of Rubinstein-Sarnak distributions

Not found. Searches on unimodal, log-concave, "symmetric unimodal", "Pearson random walk" and Rubinstein-Sarnak returned only the existence, smoothness and absolute-continuity statements (Rubinstein-Sarnak; Devin) and tail estimates (Lamzouri). Smoothness under LI is stated by Rubinstein-Sarnak, but I found no statement about shape.

Cautions from general probability (my own knowledge, not verified in the literature search):
- Wintner (1938): the convolution of symmetric unimodal distributions is symmetric unimodal. This does not apply directly, because each summand 2a_γ cos θ_γ is an arcsine law, which is U-shaped, not unimodal.
- Hence the project's remark (Anderson's inequality when the head sum is symmetric unimodal) is the right way to proceed. A general proof would need a different tool (e.g. the Pearson random walk projection results). Cited names for that line of work, unchecked: Pearson random walks, Borwein-Straub-Wan-Zudilin.

## Q6. User's own submission

Not visible. Searches for "weighted prime race", Sharma, "dissolution law" and the arXiv API returned nothing. The arXiv author-listing query for "Sharma" and prime/Chebyshev/logarithmic density returned only unrelated Sharma papers (e.g. 2610.01422, 2608.27318). No identifier to report.

## Must-cite list (beyond the existing core citations)

Weighted-sum neighbours (cite for the m>0 regime):
1. J. Gerard, L. C. Washington, "Sums of powers of primes", Ramanujan J. 45 (2018) 171-180.
2. L. C. Washington, "Sums of powers of primes II", arXiv:2209.12845; Ramanujan J. 65 (2024) 783-795. (Thm 3 is the closest published relative of (A).)
3. M. Boran, J. Byun, Z. Li, S. J. Miller, S. Reyes, "Sums of powers of primes in arithmetic progression", arXiv:2309.16007; PUMP J. Undergrad. Res. 7 (2024) 29-50.

Method precedents for (A):
4. Hardy-Littlewood (Acta Math. 41, 119-196; dated 1917 by Platt-Trudgian), Landau (Math. Z. 1 (1918), 1-24 and 213-219, two parts), as cited by Platt-Trudgian.
5. A. Fujii, "Some generalizations of Chebyshev's conjecture", Proc. Japan Acad. Ser. A Math. Sci. 64 (1988), no. 7, 260-263 (corrected from "1991"), and D. Platt, T. Trudgian, arXiv:1802.00089; IJNT 15 (2019).
6. D. R. Johnston, arXiv:2201.06184; Canad. Math. Bull. 66 (2023) 185-195.
7. T. Zhao, arXiv:2411.18903; Res. Number Theory 11 (2025) paper 62.
8. M. Suzuki, arXiv:2411.07436; Ramanujan J. 68 (2025) 95.
9. Liu-Wang-Shi, Zenodo 10.5281/zenodo.22208026 (Aug 2026), if the integrated-sign claim is used.

Recent weighted-race literature:
10. Hayani, arXiv:2512.23302 (GRH only, weight p^{-1/2}).
11. Bhattacharya-Martin-Simpson, arXiv:2507.13504 (already cited).
12. Bailleul-Hayani-Untrau, arXiv:2603.20093 (if rates are discussed).
13. Koyama, arXiv:2607.28931 (to show awareness; no overlap seen).

Suggested wording change: in the discussion of (A), say the identity and mechanism are "in the spirit of Hardy-Littlewood, Fujii, Platt-Trudgian, Johnston and Washington II", and that what is new is the application to inheritance of leads across weights. Keep "(B) and (C-general-b) not found in the literature searched" as an honest claim.
