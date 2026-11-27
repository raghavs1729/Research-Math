# SURVEY-C — recent-activity sweep, open identities outside Sun's lists, Au's WZ-seed method

Survey date: 2026-11-24. Scope: part (C) of the literature survey (everything except (A) the
Guillera/Gourevitch/Cullen 1/pi^k families and (B) Z.-W. Sun's own conjecture lists).

How this was produced
- arXiv API sweep (export.arxiv.org), about 130 queries, sorted by submission date, restricted to papers
  submitted or updated 2023-01 .. 2026-11 (raw results: `api/sweep.json`, `api/sweep2.json`; ~800 records,
  ~250 relevant after keyword filtering). Only title/abstract were read unless marked **[read]**
  (PDF downloaded, `pdftotext`, relevant parts read) or **[TeX]** (arXiv source read).
- MathOverflow via api.stackexchange.com (15 search strings, 2022-01 .. 2026-11; `mo/mo.json`, answers and
  comments for the relevant questions in `mo/ans1.json`).
- Every identity in Section 2 was checked numerically with `verify_C.py` (mpmath, 100-digit working
  precision, agreement reported up to 80 digits; output in `verify_C.out.txt`). 63 entries, 0 mismatches.
- Nothing here is a proof. "OPEN" means: stated as conjectural/unproved in the cited source and I did not see
  a proof in the sweep. "OPEN?" means the status is genuinely unclear (reason given). The arXiv API does not
  search full text, and journals-only papers (e.g. many of W. Chu's) are not covered, so a proof may exist
  that I did not see.

Notation: (a)_n = a(a+1)...(a+n-1); H_n^{(s)} = sum_{k=1}^n k^-s, H_n = H_n^{(1)};
U_n^{(s)} = sum_{i=0}^{n-1} (i+1/6)^-s; T_n^{(s)} = 2^s H_{2n}^{(s)} - H_n^{(s)};
C(n,k) binomial; G = Catalan; L_d(s) = sum_{n>=1} (d/n) n^-s (Kronecker symbol), so L_{-4}(2) = G,
L_{-4}(4) = beta(4); zeta(bar5,1) = sum_{n>m>=1} (-1)^n/(n^5 m);
Li_{s1,s2}(a,b) = sum_{n1>n2>=1} a^n1 b^n2/(n1^s1 n2^s2).

---------------------------------------------------------------------------------------------------

## 1. Chronological sweep 2023-2026

Format: arXiv id — authors — title — what it proves / poses. (Sun-list items are only summarised; they
belong to survey B. Items about pure supercongruences, partitions, Nahm sums etc. were dropped.)

### Background items updated in the window (pre-2023 submissions)
- 2101.12592 — Cohen, Guillera — Rational hypergeometric Ramanujan identities for 1/pi^c: survey and
  generalizations **[read]** — the reference table of the 13 convergent + 8 "divergent" rational formulas;
  conjectural Taylor expansions with L(D,1), L(D,2) coefficients (Sect. 6). Status table: see Section 4.
- 2201.01676 (v4 Aug 2026) — Au — Iterated integrals and multiple polylogarithm at algebraic arguments
  **[read parts]** — CMZV machinery; proves many of Sun's Apery-type conjectures (Props 7.10-7.13) and
  Borwein's sum (-1)^n/(n^4 C(2n,n)), sum (-1)^n/(n^5 C(2n,n)). Poses Galois-descent Conjectures 1.5, 1.6
  for multiple polylogarithms (not hypergeometric-series identities).
- 2212.02986 — Au — Multiple zeta values, WZ-pairs and infinite sums computations (Ramanujan J. 2025) — the
  first "seed" paper: many Sun conjectures via Gauss2F1/Dixon3F2/Watson3F2/Dougall5F4.
- 2210.07238 (v9) — Sun — Series with summands involving harmonic numbers — source of "Conjecture 14/29/38"
  quoted below (survey B).
- 2211.11484 — Wei — On two conjectural series for pi and their q-analogues — proves two Sun series.

### 2023
- 2302.00757 — Cooper — Apery-like sequences defined by four-term recurrence relations — new sequences and
  series for 1/pi (modular; title/abstract only).
- 2303.05402 — Wei — On a conjectural series of Sun for beta(4) — proves Sun's beta(4) series.
- 2303.07887 — Wei — Some fast convergent series for zeta(4) and zeta(5) — proves 4 Sun conjectures.
- 2304.00360 — Campbell — On a conjecture on a series of convergence rate 1/2 **[read parts]** — proves a
  Sun harmonic-number conjecture; mentions as open the symbolic evaluation of
  sum (1/32)^k C(2k,k)^2 H_{2k}^{(2)} and sum (1/32)^k C(2k,k)^2 H_{2k}^2 (no conjectured value given).
- 2304.09753 — Wei — On some conjectures of Z.-W. Sun involving harmonic numbers — proves ten.
- 2305.00498 — Hou, He, Wang — Ramanujan-inspired series for 1/pi involving harmonic numbers — proves five
  Sun conjectures by differentiating WZ pairs.
- 2306.02641 — Wei — eight Sun conjectures (binomials and harmonic numbers).
- 2306.04638 — Zhou — Sun's series via cyclotomic multiple zeta values — proves/generalises several Sun
  conjectures (levels 4, 8, 12, 16, 24).
- 2308.06440 — Wei, Xu — harmonic numbers of 3-order — proves three Sun conjectures.
- 2310.03699 — Hou, Sun — Taylor coefficients and series involving harmonic numbers — 58 identities, 8 of
  them Sun conjectures (coefficient extraction from parametric WZ identities).
- 2310.04642 — Wei, Xu — two conjectural series (zeta(7), zeta(3)^2) of Sun proved.
- 2310.05112 — Campbell — Proofs of conjectures on Ramanujan-type series of level 3 **[read parts]** — proves
  three of Aldawoud's quadratic-irrational level-3 conjectures; leaves open the N=13 and N=59 cases of
  Aldawoud's Table 3.7 (needs lambda*(39), alpha(39), lambda*(177), alpha(177)). Modular, not WZ. See 2602.09352.
- 2312.14051 (v4 Oct 2024; J. Symb. Comput. 130 (2025) 102421) — Au — Wilf-Zeilberger seeds and non-trivial
  hypergeometric identities **[read]** — defines WZ seeds; proves Gourevitch 1/pi^3, Cullen and Zhao 1/pi^4,
  Zhao's zeta(5), Guillera's 1/pi^2 (1930n^2+549n+45), Guillera's zeta(3) (483n^2-245n+30), several Sun
  conjectures. Names as out of reach: 1/pi^2 with 7^-4 (CG-10 below). See Section 3.
- 2312.16827 (v14) — Guillera — Heuristic derivation of Zudilin's supercongruences for rational Ramanujan
  series **[read parts]** — bilateral series; conjectures that the x^{2m+1} Taylor coefficient of the
  extended series is r*L(chi,m+1) (see E1, E2 below).

### 2024
- 2401.08997 — Broadhurst — Five families of rapidly convergent evaluations of zeta values — Lambert-series
  families (not hypergeometric); proofs included.
- 2401.14197 — Sun, Zhou — proof of Sun's conjectures with 8^k C(2k,k)/(C(3k,k)C(6k,3k)) (level-8 CMZV).
- 2402.08485 — Campbell — Ramanujan-type series associated with Heegner number 163 — proves Bagis-Glasser
  conjecture. 2403.07291 — Campbell — An extension of the Chudnovsky algorithm.
- 2402.08693 — Hakimoglu — beta-integral acceleration of 2F1; compendium of series for pi (proved).
- 2403.04555 (v4 Jul 2026; Adv. Math. 502 (2026) 111181) — Au — q-analogues of WZ seeds and Ramanujan
  1/pi^k formulas **[read parts]** — q-analogues of all seeds and of the 1/pi^k formulas; states the p-adic
  companions (e.g. sum_{n<p} (1/2)_n^7/(1)_n^7 2^-6n (6n+1)(28n^2+8n+1) == (-1/p) p^3 mod p^7) as open.
- 2403.09729 — Yamamoto — proof and generalisation of 38 Ramanujan Machine continued-fraction conjectures.
- 2403.16945 — Campbell, Glasser, Zhou — inverse binomial series via CMZVs (solves open problems of Chen).
- 2404.00711 / 2411.15116 / 2607.25173 — Allen, Grove, Long, Tu (+Rosen) — Explicit hypergeometric-
  modularity method I-III — 3F2(1), 4F3(1) values as L-values of explicit eta-quotients (proved); 2412.07054,
  2502.08760 — Rosen — L-values of weight-3 / weight-2 CM forms as 3F2(1), 2F1(1) (one auxiliary irrationality
  Conjecture 4.6 in 2412.07054; no open series identities).
- 2404.16199 — Au, Charlton — creative telescoping and generating functions of MZV variants.
- 2405.02430 — Chen, Koutschan, Wang — How to generate all possible rational WZ forms? (structure theory).
- 2405.02776, 2505.14392, 2502.15249, 2407.00621, 2509.05897, 2602.14314 — Campbell — hypergeometric
  accelerations; three-parameter generalisations of Guillera's 1/pi^2 formulas; q-analogues (all proved).
- 2410.19289 — van Hoeij, Tsai, Ye — modularity of Sun's conjectural 1/pi formulas (open since 2011).
- 2410.19186 — Cooper, Huber, Opoku — Ramanujan-Fine integrals for level 10. 2408.08844 — Babei, Roy,
  Swisher et al. — supercongruences from Ramanujan-Sato series.
- 2410.07107 — Candelas, de la Ossa, McGovern — conjectural identities: ratios of weight-4 L-values as sums of
  Gromov-Witten invariants (not hypergeometric-term series; out of WZ scope).

### 2025
- 2501.10090 — Cohen, Zudilin — Variations on a theme of Apery (continued fractions; proved).
- 2503.00570 (Bull. Aust. Math. Soc. 113 (2026) 152-160) — Guillera — The WZ method and flawless WZ pairs
  **[read parts]** — passing to G~(n,k) = G(n-k,k) etc. removes the boundary terms; simpler proofs of Cullen
  1/pi^4, the (1930n^2+549n+45) 1/pi^2 formula, Gourevitch 1/pi^3; Maple code.
- 2503.05599 (Ramanujan J. 69:83, 2026) — Zhou — Notes on certain binomial harmonic sums of Sun's type —
  proves Sun's Conjecture 29 type ((log m)/pi) identities by modular parametrisation.
- 2504.01975 — Guillera — Fast formulas for zeta(2,a) and zeta(3,a) (WZ; proved).
- 2506.08245, 2506.10321 — Zuniga — Fast Ramanujan-type series for logarithms I, II **[read parts]** —
  WZ/beta-integral proofs; one log 2 series (eq. (23)) explicitly left unproved (C1 below).
- 2506.01865 — Sun, Zhou — fast irrational series for L_d(2) (Guillera-Rogers theory).
- 2506.19221 — Feng, Hou — Finding congruences with the WZ method.
- 2507.19502 — Chen, Nguyen — product of four consecutive harmonic numbers (elementary; Furdui conjecture).
- 2508.00343 — Valloud — Streamlined WZ proofs of Van Hamme supercongruences (procedure to find WZ pairs).
- 2508.05770 — Shalev — log-sine integrals and Apery-like sums (proved).
- 2509.19960 — Au — Linear relations between L-values of newforms and moments of elliptic K **[read parts]**
  — proves e.g. L(eta(2tau)^12, 5) = (pi^6/1024) sum_{n>=0} (4n+1) (1/2)_n^8/(1)_n^8; conjectures the
  supercongruence sum_{n<p} (1/2)_n^8/(1)_n^8 (1+4n) == p a_p(f) (mod p^6), p >= 7.
- 2510.15008, 2601.16232, 2601.20027, 2603.21079 — Gonzalez Layja — Apery-like series of weight 4, 5 (proved).
- 2502.17533 — Raz et al. — From Euler to AI: unifying formulas for constants (conservative matrix fields).
- 2502.17840 — Xiong et al. — LeanComb: combinatorial identities benchmark in Lean.

### 2026
- 2601.12495 — Zuniga — Fast computing formulas for some Dirichlet L-series (uses one Sun conjecture for
  L_{-3}(2) as "fast conjectural series").
- 2601.18645 — Jiang, Sun — over ten Sun conjectures with C(4k,k) and harmonic numbers proved.
- 2602.06679 — Zudilin — A supercongruence fantasy on Fibonacci and Lucas (after Guillera).
- 2602.08721 (v4, 9 Nov 2026) — Au — Discovering hypergeometric series with harmonic numbers via WZ seeds
  **[TeX read]** — ten families; fast series for zeta(5,3); Hilbert-series Conjectures 1.1, 1.3, 3.x, 4.1-4.5;
  about 40 explicit identities marked "?=" (found from WZ pairs, proofs not completed because g(n) and
  lim_n sum_k F(n,k) were not evaluated). **Main source of Section 2.** Mathematica package WZSeedTools.
- 2602.09352 — Campbell, Cooper, Ye — Quadratic irrational analogues of Ramanujan's series for 1/pi
  **[read parts]** — large classification, with proofs by modular methods, of quadratic-irrational series.
  I did not check whether Aldawoud's N=13, N=59 cases are included.
- 2602.12091 — Sun, Zhou — Series involving central binomial coefficients and higher-order harmonic numbers
  **[read parts]** — modular proofs of e.g. sum C(2k,k)^3 (H_{2k}^{(3)} - (43/352) H_k^{(3)}) (42k+5)/4096^k
  = 555 zeta(3)/(77 pi) - 32G/11 (Sun's Conjecture 38) and the L_{-7}(2) identity.
- 2603.29973 — Sun — Various conjectural series identities (150+ new conjectures; survey B).
- 2604.03327, 2609.22304 — Le Lan — proofs of Sun's series (Cauchy product; beta function); leaves open Sun's
  supercongruence and Sun's (4.25) [sum (U(k)-V(k)/(4k))/(k(3k-1)(3k-2)(-8)^k C(4k,k)) = -pi^2/4].
- 2604.04006 — Qi — proof of two supercongruences of Guillera and Zudilin.
- 2604.06239, 2604.09723, 2604.14219 — Shvets — Ramanujan-Machine conjecture Z2 and Apery limits (modular).
- 2604.15172 — Hou, Sun — Evaluations of some series via the WZ method — proves
  sum (28k^2+10k+1) C(2k,k)^5/((6k+1)(-64)^k C(3k,k)C(6k,3k)) = 3/pi and the derivative sums
  sum_k d^m/dk^m (...) conjectured on MO 508743.
- 2604.22377 — Hou, Mu — A construction method for WZ seeds **[TeX read]** — seven seeds with one more
  parameter than Au's (bilateral Dougall 2H2, 5H5, 7H7 type; listed in Section 3.4); "transpose" theorem
  (G(k,n) is again a seed); q-analogues.
- 2605.04472 — Xiong, Lv, Liu, Wang — WZ-LLM: WZ proof plans -> Lean 4 sketches (34% on 100 classic
  identities). 2607.01544 — Lau, Lee, Ono — Formalized q-series in Lean (Rogers-Ramanujan).
- 2606.28966 — Mendusic — cubic Ramanujan-Sato series at D=-163 on Gamma_0(2)^+.
- 2607.09721 — Shalyt, Kalisch, Schneider, et al. — The Ramanujan Challenge for AI **[read]** — 8 problems
  with withheld proofs (continued fractions, Apery limits, one harmonic double sum) and 2 open conjectures
  (a knot A-polynomial integral = 4 pi^2/85; gcd(d_n a_n, d_n b_n) = e^{o(n)} for Apery's sequences); no open
  hypergeometric-term series. 2610.40551 — results: one of the two open conjectures reportedly proved by a
  submission (abstract; I did not determine which).
- 2607.10135 — Bhandari — recurrences for hypergeometric moments and binomial-harmonic sums.
- 2607.11924, 2607.24867 — Fernandez Refolio — analytic evaluation of classical Ramanujan series.
- 2607.14646 — Sun, Zhou — Multiple Clausen values and deformed Apery-like series (proved).
- 2607.15303, 2607.27229 — Patel — Sun's derivative-sum Conjectures 4.1 and 4.2-4.6 (the latter via exact WZ
  certificates and Au's seeds; four misprints in Sun's statements corrected; 4.6 "through weight eleven").
- 2609.09209 — Adegoke; 2611.28059 — Adegoke, Frontczak, Goy — families generalising Ramanujan/Forsyth-
  Glaisher series with harmonic numbers (proved).
- 2610.00551 — Campbell — solves the level-11 Ramanujan-Sato coefficient problem of Cooper et al.
- 2610.04770 — Shalev — Borwein-Bradley-type identities for Ramanujan-type series for 1/pi^4 **[read
  parts]** — parametric extensions of Cullen/Zhao via Au's Dougall7F6 seed on the "Gamma reflection
  hyperplane", proved by Guillera's periodicity argument; proves ten Sun conjectures (Conj. 70(i), 69 of
  [Sun 2210.07238]; 5.47(i) and the last two of 5.45 of [Sun's 2021 book]); the first two of 5.45 (zeta(5))
  are not obtained. Announces further numerically found identities.
- 2610.16844 — Chen, Liu, Lu — The WZ pair and finite congruences — meromorphic WZ pair for the
  Guillera-Rogers series for 3 ln 2 (completes the WZ proof left open there).
- 2610.17240 — Au — Holonomic sequences arising from meromorphic modular forms **[read parts]** — modular
  creative telescoping; Conjectures 1.6/1.7/6.1/6.6 on denominators (not series identities).
- 2610.27799 — Leibtag et al. — New directions in irrationality (CMFs; mu(zeta(3)) < 5.4197 claimed).
- 2610.34482 — Shalev — A symmetric extension of Zhao's series for zeta(5) **[read]** — four-parameter
  extension (S_5-symmetric summand from Seed7 on a hyperplane); new companion
  sum_{n>=0} (-1)^n (2n)!^10 (5n)!/(n!^15 (10n)!) * 3(5532n^4+5464n^3+2173n^2+416n+32)/((10n+1)(10n+3)) = 240/pi^2;
  proves five of the six "?=" identities of Au 2602.08721 Ex. X; the sixth (zeta(bar5,1)) stays open.
- 2611.16156 — Fischler, Lai, Sprang, Zudilin — Many more odd zeta values are irrational.

### MathOverflow 2022-2026 (non-Sun authors; Sun's own questions are listed at the end of Section 2)
- 281009 (Y. Zhao, 2017): source of the Zhao zeta(5) and 1/pi^4 series — proved (Au).
- 466853 (Zuniga, 2024-03, 0 answers): unproved log 2 series (C1). 487156 (Zuniga, 2025-02): log p series.
- 473746 (Guillera, 2024-06, 0 answers): L-value Taylor-coefficient conjecture (E1, E2).
  473848 (Guillera, 2024-06, 0 answers): new supercongruence pattern. 480598 (Guillera): primitive multiplier.
- 479718 (H. Cohen, 2024-09): Gamma(1/4)^4/pi^3 series — answered (Gosper-summable / Dixon).
- 486961 (Deyi Chen, 2025-01/02): five series; (1),(2),(4) known/proved per Guillera and Sun; (3),(5) no proof
  posted (D1, D2). 507693 (Deyi Chen, 2026-02, 0 answers): D3-D5.
- 501176 (2025-10): sum (2n)!^4/(2^8n n!^8) (8n^2+10n+3)/(n+1)^3 = 32/pi^2 — proved in the question itself.
- 508743, 508825 (Deyi Chen, 2026-03): derivative sums — proved (Hou-Sun 2604.15172; Beenakker's answer).
- 508879 (2026-03): S(k) = sum C(2n,n)C(4n,2n)/((kn+1)64^n); S(1), S(2), S(3) known (Yabu); k >= 4 no
  conjectured closed form.

---------------------------------------------------------------------------------------------------

## 2. Open identities found (exact ASCII statements)

Columns of information for each item: statement; source; status; id in `verify_C.py`; digits verified.
All items below agree to >= 76 digits (80 = the cap) in `verify_C.out.txt`. "Digits" for identities with
right-hand side 0 means |lhs| < 10^-80.

### 2.A  K. C. Au, arXiv:2602.08721v4 (9 Nov 2026), Section 4 — identities carrying "?=" in the paper

Au's own description of their status: each is the coefficient extraction of sum_n G(n,0) for an explicit WZ
pair built from a seed, but "explicit expressions for g(n) and lim_{n->oo} sum_k F(n,k) are unavailable", so
"strictly speaking, conjectural". These are therefore the most promising WZ targets: the WZ pair is known,
what is missing is the evaluation of the boundary terms (Au's asymptotic technique or Guillera's flawless
pairs). I found no later paper proving them except where noted.

**Example VI.** a_n = (1/4)^n (1)_n^3/(1/2)_n^3, sums over n >= 1. WZ pairs (with d the shift of k):
F1 = Seed2(a-n+1, b+1/2, c-n+2, k+d+n-1/2), F2 = Gauss2F1(a-n+1, b-n+1, c+1, k+d+n-1/2). Base identity
sum a_n (3n-1)/n^3 = pi^2/2. (Equations (4.2)-(4.5), (4.9) of the paper are proved: Sun's conjectures,
Li-Chu, Electron. Res. Arch. 31 (2023).)

- (4.6)  sum a_n [ -(17n-6)/(6n^5) - (3n-1) H_n (n H_n - 2)/n^4 ]
         = -(128/9) Li_4(1/2) + 533 pi^4/6480 - (16/27) log^4 2 - (8/27) pi^2 log^2 2.
         [Au-VI-(4.6); OPEN; 80 digits]
- (4.7)  sum a_n [ -(7n-2)/(4n^5) + (4n-1) H_n/(2n^4) - (3n-1)(n H_n - 1) H_{2n}/n^4 ]
         = -(112/9) Li_4(1/2) + 629 pi^4/12960 - (14/27) log^4 2 + (2/27) pi^2 log^2 2.   [OPEN; 80]
- (4.8)  sum a_n [ -(19n-6)/(24n^5) + (4n-1) H_{2n}/n^4 - (3n-1) H_{2n}^2/n^3 ]
         = -(80/9) Li_4(1/2) - 17 pi^4/5184 - (10/27) log^4 2 + (4/27) pi^2 log^2 2.      [OPEN; 80]
- (4.10) sum a_n [ -(13n-5)/(3n^6) + (5n-2) H_n/(2n^5)
                   - (3n-1)(-3 H_n^{(2)} + 3n H_n H_n^{(2)} - 2n H_n^{(3)})/(3n^4) ]
         = -17 pi^2 zeta(3)/12 + 341 zeta(5)/24 + pi^4 log 2/36.                           [OPEN; 80]
- (4.11) sum a_n [ -(29n-13)/(6n^6) + (n-2) H_n/(4n^5)
                   + (3n-1)(5n H_n^{(3)} + 6 H_{2n}^{(2)} - 6n H_n H_{2n}^{(2)})/(3n^4) ]
         = -19 pi^2 zeta(3)/8 + 713 zeta(5)/48 + 13 pi^4 log 2/72.                         [OPEN; 80]
- (4.12) sum a_n [ (19n-7)/(3n^6) - (17n-6) H_n/(2n^5)
                   - (3n-1)(-9 H_n^2 + 3n H_n^3 + 4n H_n^{(3)})/(3n^4) ]
         = -(512/9) Li_5(1/2) + 67 pi^2 zeta(3)/36 + 1705 zeta(5)/72 + 64 log^5 2/135
           + (32/81) pi^2 log^3 2 - 533 pi^4 log 2/1620.                                   [OPEN; 80]
- (4.13) sum a_n [ (200n-71)/(12n^6) - (19n-6) H_{2n}/(4n^5) + 3(4n-1) H_{2n}^2/n^4
                   - (3n-1)(6 H_{2n}^3 + 17 H_n^{(3)})/(3n^3) ]
         = -(320/9) Li_5(1/2) + 545 pi^2 zeta(3)/72 - 24149 zeta(5)/288 + 8 log^5 2/27
           - (16/81) pi^2 log^3 2 + 17 pi^4 log 2/1296.                                    [OPEN; 80]
- (4.14) sum a_n [ -(7n-3)/(2n^7) + (5n-2) H_n^{(2)}/n^5 - (3n-1)(2 (H_n^{(2)})^2 - H_n^{(4)})/(2n^3) ]
         = -pi^6/1440.                                                                     [OPEN; 80]
- (4.15) sum a_n [ (9n-1)/(4n^7) + (n-2) H_n^{(2)}/(4n^5) + (5n-2) H_{2n}^{(2)}/n^5
                   - (3n-1)(3 H_n^{(4)} + 8 H_n^{(2)} H_{2n}^{(2)})/(4n^3) ] = -17 pi^6/2880.  [OPEN; 80]
- (4.16) sum a_n [ (73n-21)/(8n^7) + (n-2) H_{2n}^{(2)}/n^5
                   - (3n-1)(21 H_n^{(4)} + 32 (H_{2n}^{(2)})^2 + 32 H_{2n}^{(4)})/(8n^3) ]
         = -337 pi^6/5760.                                                                 [OPEN; 80]
- (4.17) sum a_n [ -(3n-1)(n H_n - 1)(-8 + 7n^3 H_n^{(3)} + 8n^3 H_{2n}^{(3)})/(8n^7) ]
         = -7 zeta(3)^2/3 + (2/3) pi^2 zeta(3) log 2 - pi^6/224.                           [OPEN; 80]
- (4.18) sum a_n [ (299n-129)/(36n^7) - (29n-13) H_n/(3n^6) + (n-2) H_n^2/(4n^5)
                   - (17n-6) H_{2n}^{(2)}/(3n^5)
                   + (3n-1)(-40 H_n^{(3)} + 40n H_n H_n^{(3)} + 3n H_n^{(4)} + 48 H_n H_{2n}^{(2)}
                            - 24n H_n^2 H_{2n}^{(2)})/(12n^4) ]
         = (184/9) zeta(bar5,1) - 179 zeta(3)^2/6 + (19/3) pi^2 zeta(3) log 2 - (713/18) zeta(5) log 2
           + 827 pi^6/25920 - (13/54) pi^4 log^2 2.                                        [OPEN; 80]
  (The bar on the 5 is in the TeX source; it is lost in the PDF text layer.)

**Example VII.** a_n = 2^-6n (1/2)_n^3/(1)_n^3, p0 = 42n+5, p2 = 35+162n, p5 = 77+366n, sums over n >= 0.
WZ pairs: F1 = Seed2(a-n+1/2, b+1/2, c+n+1, d+k+n), F2 = Dixon6F5(a-n+1/2, b+1, c+1, d+k+n+1/2),
F3 = Seed9(a-2n+1, b-n+1/2, c-n+1/2, k+d+n-1/6). Base: sum a_n (42n+5) = 16/pi.
((4.19) is a case of Sun's Conjecture 29, proved by Zhou; (4.22), (4.23) proved by Li-Chu 2024; (4.25) not
marked.)

- (4.20) sum a_n [ (31+162n)/(3+18n) + (-H_n/2 - H_{3n} + 2 H_{6n}) p0 ] = 48 log 2/pi.
         [Au-VII-(4.20); OPEN?: marked "?=", but Au's introduction says it "can be proved by both the WZ and
         modular approaches [Zhou]"; 80 digits]
- (4.21) sum a_n [ (35+162n)/(1+6n) + (-H_n/2 + 2 U_n^{(1)}) p0 ] = 32/sqrt(3) + 80 log 2/pi.  [OPEN; 80]
- (4.24) sum a_n [ (43+396n+972n^2)/(9(1+2n)(1+6n)^2) + (-H_{3n}^{(2)}/4 + H_{6n}^{(2)}) p0 ] = 14 pi/9.
         [OPEN; 80]
- (4.26) with X = H_n + 2 H_{3n} - 4 H_{6n}:
         sum a_n [ (6n+5)(78n+17)/(9(2n+1)(6n+1)^2) + (p0/4) X^2 - (162n+31) X/(18n+3) ]
         = 144 log^2 2/pi - 34 pi/9.                                                       [OPEN; 80]
- (4.27) sum a_n [ 14 U_n^{(1)} - 8/(6n+1)^2 - 23/(6n+3)
                   + (1/6)(3 H_n^2 + 12 H_{2n} U_n^{(1)} - 3 H_n (H_{2n} + 4 U_n^{(1)}) - 2 U_n^{(2)}) p0
                   + H_{2n} p2/(1+6n) - H_n p5/(2+12n) ]
         = -120 L_{-3}(2)/pi - 104 pi/9 + 160 log^2 2/pi + 64 log 2/sqrt(3).               [OPEN; 80]
- (4.28) sum a_n [ -9/(1+2n)^2 + (1/3)(-2 H_n^{(3)} + 17 H_{2n}^{(3)}) p0 ] = 80 zeta(3)/pi - 128 G/3.
         [OPEN?: Au: "would follow from two conjectures of Sun"; one of them (Conj. 38) is now proved by
         Sun-Zhou 2602.12091, the status of the other (in Sun 2307.03086) I did not determine; 80]
- (4.29) sum a_n [ (-479 - 1275(1+2n) H_n + 1275(1+2n) H_{2n} + (357 + 1428n + 1428n^2) H_{2n}^{(2)})/(1+2n)^2
                   + (-4 H_n^{(3)} + 51 (H_{2n} - H_n) H_{2n}^{(2)}) p0 ]
         = -2432 G - 64 zeta(3)/pi + 816 pi log 2.                                         [OPEN; 80]
- (4.30) sum a_n [ 25/(9(1+2n)^3) + (-H_n^{(4)}/16 + H_{2n}^{(4)}) p0 ] = 5 pi^3/54.
         [OPEN; also Sun 2307.03086 (5.14); 80]
- (4.31) sum a_n [ 4(-5620 + 621(1+2n)^2 H_n^{(2)})/(27(1+2n)^3) + (1/2)((H_n^{(2)})^2 + H_n^{(4)}) p0 ]
         = -10874 pi^3/405.                                                                [OPEN; 80]
- (4.32) sum a_n [ 25 H_n^{(2)}/(1+2n) + 4(-3103 + 621(1+2n)^2 H_{2n}^{(2)})/(27(1+2n)^3)
                   + (H_n^{(4)}/4 + H_n^{(2)} H_{2n}^{(2)}) p0 ] = -6002 pi^3/405.          [OPEN; 80]
- (4.33) sum a_n [ 5(-673 + 270(1+2n)^2 H_{2n}^{(2)})/(27(1+2n)^3) + (H_n^{(4)}/16 + (H_{2n}^{(2)})^2) p0 ]
         = -3251 pi^3/810.                                                                 [OPEN; 80]
- (4.34) sum a_n [ 9 p0 (16 (H_n - H_{2n})(2 H_n^{(3)} - 17 H_{2n}^{(3)}) - H_n^{(4)})
                   + 3888 (H_n - H_{2n})/(2n+1)^2 - 2016 H_n^{(3)} + 17136 H_{2n}^{(3)} + 9392/(2n+1)^3 ]
         = 73728 Im Li_3(1/2 + i/2) + 69120 zeta(3) log 2/pi - 22904 pi^3/15 - 2304 pi log^2 2.  [OPEN; 80]
- (4.35) sum a_n [ H_n^{(2)} ((17 H_{2n}^{(3)} - 2 H_n^{(3)}) p0 - 27/(2n+1)^2)
                   + 4 (5 p0 H_n^{(5)} + (391 H_{2n}^{(3)} - 46 H_n^{(3)})/(2n+1) - 159 p0 H_{2n}^{(5)}
                        - 4/(2n+1)^4) ]
         = 1536 L_{-4}(4) + 440 pi zeta(3) - 9856 zeta(5)/pi.                              [OPEN; 80]
- (4.36) sum a_n [ (H_{2n} - H_n)((144 H_{2n}^{(4)} - 9 H_n^{(4)}) p0 + 400/(2n+1)^3)
                   - 63 H_n^{(4)} + 6(42n+5) H_n^{(5)} - 8064 n H_{2n}^{(5)}
                   + 16 (63 H_{2n}^{(4)} - 60 H_{2n}^{(5)} - 25/(2n+1)^4) ]
         = (80/3) pi^3 log 2 - 2976 zeta(5)/pi.                                            [OPEN; 80]
- (4.37) sum a_n [ -50 H_n^{(3)}/(2n+1) + 425 H_{2n}^{(3)}/(2n+1)
                   + p0 (6 H_n^{(5)} - 2 H_n^{(3)} H_{2n}^{(2)} + 17 H_{2n}^{(2)} H_{2n}^{(3)} - 191 H_{2n}^{(5)})
                   - 27 H_{2n}^{(2)}/(2n+1)^2 - 2/(2n+1)^4 ]
         = 512 L_{-4}(4) + 120 pi zeta(3) - 2960 zeta(5)/pi.                               [OPEN; 80]
- (4.38) sum a_n [ p0 (-27 H_n^{(2)} H_n^{(4)} - 172 H_n^{(6)} + 432 H_n^{(2)} H_{2n}^{(4)} + 11008 H_{2n}^{(6)})
                   + 1200 H_n^{(2)}/(2n+1)^3 - 2484 (H_n^{(4)} - 16 H_{2n}^{(4)})/(2n+1) - 192/(2n+1)^5 ]
         = 44 pi^5/15.                                                                     [OPEN; 80]
- (4.39) sum a_n [ p0 (27 (H_n^{(2)})^3 + 81 H_n^{(4)} H_n^{(2)} + 13022 H_n^{(6)} - 829952 H_{2n}^{(6)})
                   + 7452 (H_n^{(2)})^2/(2n+1) - 134880 H_n^{(2)}/(2n+1)^3 + 7452 H_n^{(4)}/(2n+1)
                   + 428928/(2n+1)^5 ] = 120706 pi^5/105.                                   [OPEN; 80]
- (4.40) sum a_n [ p0 (-155 H_n^{(6)} - 108 H_n^{(4)} H_{2n}^{(2)} + 1728 H_{2n}^{(2)} H_{2n}^{(4)} + 9920 H_{2n}^{(6)})
                   - 2700 H_n^{(4)}/(2n+1) + 4800 H_{2n}^{(2)}/(2n+1)^3 + 43200 H_{2n}^{(4)}/(2n+1)
                   - 1152/(2n+1)^5 ] = 0.                                                   [OPEN; 80]
- (4.41) sum a_n [ p0 (4 (H_n^{(3)})^2 - 68 H_{2n}^{(3)} H_n^{(3)} + 289 (H_{2n}^{(3)})^2 + 2 H_n^{(6)} - 127 H_{2n}^{(6)})
                   + 108 H_n^{(3)}/(2n+1)^2 - 918 H_{2n}^{(3)}/(2n+1)^2 - 54/(2n+1)^5 ]
         = 3072 Im Li_{4,1}(i,1) + 1024 Im Li_{4,1}(i,-1) + 2048 L_{-4}(4) log 2 + 3600 zeta(3)^2/pi
           - 9434 pi^5/945.                                                                [OPEN; 80]

**Example VIII.** a_n = 9^-n (1/2)_n (1/4)_n (3/4)_n/(1)_n^3, p0 = 1+8n, n >= 0.
WZ pair: F = Seed10(a+3b-n+1/2, 2b+1, c+k+n). Base: sum a_n (1+8n) = 2 sqrt(3)/pi.

- VIII.1 sum a_n [ 1 + (1/6)(-4 H_n + 3 H_{2n}) p0 ] = 8 log 2/(sqrt(3) pi).
         [Au-VIII-1; = Sun's Conjecture 14 of 2210.07238, so a survey-B item; OPEN?; 80 digits.
         **Misprint warning**: the paper prints a_n p0 (1 + (1/6)(-4H_n+3H_{2n})); that version is false
         numerically (1.1108... vs 1.0190...). The form above is the one that holds.]
- VIII.4 sum a_n [ -4 H_n + 3 H_{2n} + (1/24)(32 H_n^2 - 48 H_n H_{2n} + 18 H_{2n}^2 + 5 H_n^{(2)}) p0 ]
         = 32 log^2 2/(sqrt(3) pi) - 37 pi/(24 sqrt(3)).                                   [OPEN; 80]
- VIII.5 sum a_n [ -H_n + H_{4n} + (1/16)(4 H_n^2 - 8 H_n H_{4n} + 4 H_{4n}^2 + H_n^{(2)} - 4 H_{4n}^{(2)}) p0 ]
         = -3 sqrt(3) pi/16 + 2 sqrt(3) log^2 2/pi + sqrt(3) log^2 3/(8 pi) + sqrt(3) log 3 log 2/pi.
         [OPEN; 80]
- VIII.6 sum a_n [ -5 H_n^{(2)} + 18 H_{2n}^{(2)} + (1/6)(20 H_n H_n^{(2)} - 15 H_{2n} H_n^{(2)} + 14 H_n^{(3)}
                   - 72 H_n H_{2n}^{(2)} + 54 H_{2n} H_{2n}^{(2)} - 108 H_{2n}^{(3)}) p0 ]
         = 45 L_{-3}(2)/2 - 94 zeta(3)/(sqrt(3) pi) + 4 pi log 2/(3 sqrt(3)).              [OPEN; 80]
- VIII.7 sum a_n p0 [ 50 (H_n^{(2)})^2 - 360 H_{2n}^{(2)} H_n^{(2)} + 648 (H_{2n}^{(2)})^2 + 41 H_n^{(4)} - 648 H_{2n}^{(4)} ]
         = -2 pi^3/(15 sqrt(3)).                                                           [OPEN; 80]
  (VIII.2 [(2 + p0 (H_{4n} - H_n)) -> 4 sqrt(3) log 2/pi + sqrt(3) log 3/pi] is proved (Zhou);
   VIII.3 [p0(-5H_n^{(2)} + 18 H_{2n}^{(2)}) -> pi/sqrt(3)] is proved (Wei-Ruan 2024).)

**Example IX.** a_n = (-27/512)^n (1/2)_n (1/6)_n (5/6)_n/(1)_n^3, p0 = 15+154n, p1 = 61+110n, p2 = 1+6n,
n >= 0. WZ pair: F = Seed1(-a+b-c-n+1/2, b-c+1, c+k+2n). Base: sum a_n (15+154n) = 32 sqrt(2)/pi.

- IX.1 sum a_n [ (H_{2n} - H_n) p0 + p1/(3(1+2n)) ] = 64 sqrt(2) log 2/pi.                 [OPEN; 80]
- IX.3 sum a_n [ (1/15)(-4 H_n^{(2)} + 15 H_{2n}^{(2)}) p0 - p2/(1+2n)^2 ] = -16 sqrt(2) pi/45.  [OPEN; 80]
- IX.4 sum a_n [ (1/18)(9 H_n^2 - 18 H_n H_{2n} + 9 H_{2n}^2 + H_n^{(2)}) p0 - (H_n - H_{2n}) p1/(3(1+2n))
                 - p2/(2(1+2n)^2) ] = 64 sqrt(2) log^2 2/pi - 88 sqrt(2) pi/27.            [OPEN; 80]
- IX.5 sum a_n [ -2 (H_n - T_{3n}^{(1)}) + (1/154)(3 H_n^2 + H_n^{(2)} - 6 H_n T_{3n}^{(1)} + 3 (T_{3n}^{(1)})^2
                 - 3 T_{3n}^{(2)}) p0 ] = 1200 sqrt(2) log^2 2/(77 pi) - 16 sqrt(2) pi/21.  [OPEN; 80]
- IX.6 sum a_n [ (1/3)(77 H_{2n}^{(2)} + (6n-1)/(2n+1)^3) + (H_n - T_{3n}^{(1)}) p2/(2(1+2n)^2)
                 - 308 H_n^{(2)}/45
                 + (1/90)(12 H_n H_n^{(2)} + 8 H_n^{(3)} - 45 H_n H_{2n}^{(2)} - 60 H_{2n}^{(3)}
                          - 12 H_n^{(2)} T_{3n}^{(1)} + 45 H_{2n}^{(2)} T_{3n}^{(1)}) p0 ]
       = 512 L_{-8}(2)/45 - 832 sqrt(2) zeta(3)/(45 pi) - (8/9) sqrt(2) pi log 2.          [OPEN; 80]
  (IX.2 [1/3 + (1/154)(-H_n + T_{3n}^{(1)}) p0 -> 80 sqrt(2) log 2/(77 pi)] is proved (Zhou).)

**Example X.** a_n = (-2^10/5^5)^n (1)_n^9/((1/2)_n^5 (1/5)_n (2/5)_n (3/5)_n (4/5)_n), n >= 1,
p0 = 30 - 425n + 2275n^2 - 5600n^3 + 5532n^4, p1 = 54 - 680n + 3185n^2 - 6720n^3 + 5532n^4.
WZ pair: F = Seed7(a+n+1/2, b-n, c+n+1/2, d+n+1/2, e+k+n+1/2).

- X.1 sum a_n [ p1/n^10 + p0 (-3 H_n + 2 H_{2n} + H_{5n})/n^9 ]
      = 786432 zeta(bar5,1) - 294912 zeta(3)^2 - 1523712 zeta(5) log 2 + 32768 pi^6/105.
      [Au-X-1; OPEN: the only one of the six Ex. X identities not proved by Shalev 2610.34482, who writes
      that it "presumably requires the full five-parameter family"; 80 digits]

**Sun's (4.1), singled out by Au** (remark after (3.32) in Example IV; WZ pairs of Example IV are
F1 = Dougall5F4(a+n+1/2, b+1/2, c, d-n, k+d+n+1/2), F2 = Dougall5F4(a, b, c, d-n, e+k+2n+1),
F3 = Seed3(a, b-n, c, k+d+n+1), all with hypergeometric part (-1/27)^n (1)_n^2/((1/3)_n (2/3)_n)):

- S41  sum_{n>=1} (-1)^n/(n^2 C(2n,n) C(3n,n)) * [ (7n-2)/(2n-1) * H_{n-1} - 1/(6n) ]
       = pi^2 log 2/6 - 11 zeta(3)/12.
       [Sun-(4.1)/Au-IV; OPEN; Au: "evidently does not follow from the formulas above because of the presence
       of log 2. We encourage interested readers to seek other WZ-pairs that might prove it." 80 digits.
       **Misprint warning**: Au displays it as sum a_n ((7n-2)/(2n-1) H_{n-1} - 1/(6n)) with
       a_n = (-1/27)^n (1)_n^2/((1/3)_n(2/3)_n) = (-1)^n/(C(2n,n)C(3n,n)); numerically the factor 1/n^2 is
       needed (0.0681... vs 0.0383... without it).]

Structural conjectures in the same paper (not single identities; listed for completeness): Conjectures 1.1,
1.3 = 3.5, 3.1, 3.2, 3.4, 3.7, 3.8, 4.1-4.5: the Hilbert series sum_N dim V_N(S) t^N of the space of
coefficient-extracted summands is rational, e.g. 1/((1-t)(1-t^2)(1-t^3)(1-t^4)(1-t^5)) for the
(-1/2^10)^n (1)_n^5/(1/2)_n^5 family and for the Zhao zeta(5) family; Problem 1.2 (relation between WZ and
modular proofs).

### 2.B  Rational 1/pi^2 series of the Cohen-Guillera table not seen proved (cross-check with survey A)

Au, arXiv:2312.14051 (Sect. 1 and 3.4) explains that with his seed list only Seed1, Seed4, Seed10 can produce
a ratio z outside the multiplicative group generated by {n^n}, and names B1 as not provable by the technique.

- B1  sum_{n>=0} 7^-4n (1/2)_n (1/8)_n (3/8)_n (5/8)_n (7/8)_n/(1)_n^5 * (1920n^2 + 304n + 15) = 56 sqrt(7)/pi^2.
      [CG-10; Guillera; OPEN (printed with "?=" by Au, Oct 2024); 80 digits]
- B2  sum_{n>=0} (-1/80^3)^n (1/2)_n (1/3)_n (2/3)_n (1/6)_n (5/6)_n/(1)_n^5 * (5418n^2 + 693n + 29) = 128 sqrt(5)/pi^2.
      [CG-7; OPEN? — conjectural in the 2021 survey, no proof seen in my sweep; 80]
- B3  sum_{n>=0} (3/5)^6n (1/2)_n (1/3)_n (2/3)_n (1/6)_n (5/6)_n/(1)_n^5 * (532n^2 + 126n + 9) = 375/(4 pi^2).
      [CG-8; OPEN? (same remark); 80]
- B4  sum_{n>=0} (-1/2^10)^n (1/2)_n (1/4)_n (3/4)_n (1/6)_n (5/6)_n/(1)_n^5 * (1640n^2 + 278n + 15) = 256/(sqrt(3) pi^2).
      [CG-9; OPEN? (same remark); 80]
- B5  sum_{n>=0} (-1/48)^n (1/2)_n (1/3)_n (2/3)_n (1/4)_n (3/4)_n/(1)_n^5 * (252n^2 + 63n + 5) = 48/pi^2.
      [CG-5; OPEN? (same remark); 80]
  For B2-B5 I only verified the numerics and that no 2023-2026 abstract in the sweep claims a proof; whether
  older WZ proofs exist is for survey A to settle. (#6, #11, #12, #13 are proved by Au, see Section 4; that
  #1-#4 are Guillera's WZ-proved formulas is from my memory of the older literature, not re-checked here.)

### 2.C  Zuniga's unproved fast series for log 2

- C1  log 2 = (1/576) sum_{n>=0} 8^-n * (686430n^3 + 1317033n^2 + 798173n + 153712)
              /((12n+1)(12n+5)(12n+7)(12n+11)) * (3n)! (4n)! (6n)!/(n! (12n)!)
      equivalently (arXiv:2506.08245 eq. (23))
      log 2 = sum_{n>=1} (686430n^3 - 742257n^2 + 223397n - 13858)/(3n(2n-1)(3n-1)(3n-2)) * (2^13 3^3)^-n
              * (1)_n (1/2)_n (1/3)_n (2/3)_n/((1/12)_n (5/12)_n (7/12)_n (11/12)_n).
      [Zuniga-log2, Zuniga-log2-b; J. Zuniga, MO 466853 (Mar 2024: "remains un-proven. After trying to prove this
      series using Au WZ seeds with no success ...") and arXiv:2506.08245 (June 2025: "just conjectured since
      it has not been possible to prove it so far"); used to verify the 3*10^12-digit record of log 2.
      OPEN as of June 2025, nothing later seen; 80 digits. Ratio 1/221184 = 2^-13 3^-3.]

### 2.D  Deyi Chen (MathOverflow)

- D1  sum_{k>=0} (-1)^k C(6k,k) C(5k,k) C(4k,k)^2 C(3k,k) (4k+1)(6k+1)(2264k^4 + 6118k^3 + 6267k^2 + 2884k + 503)
              /(384 (2k+1)^9 C(2k,k)^10) = zeta(3).
      [Chen-486961-(3); MO 486961 eq. (3), Feb 2025, "related to Example II of arXiv:2312.14051"; OPEN? — no
      proof posted; by Guillera's answer for eq. (2) one expects it from Au's WZ pair for Zhao's 1/pi^4 series
      by n -> n+1/2; ratio -27/256; 80 digits]
- D2  sum_{k>=1} (-1)^k (28k^2 + 10k + 1) C(2k,k)^2/((6k+1)(2k+1)^2 k^3 C(3k,k) C(6k,3k)) = 24 - 20 zeta(3).
      [Chen-486961-(5); MO 486961 eq. (5); OPEN? — no proof posted; "upside-down"/half-shift companion of
      Sun's 3/pi series proved by Hou-Sun 2604.15172; ratio -1/27; 80]
- D3  sum_{k>=0} (H_{2k} - H_k) C(2k,k)^3/256^k
      = ((log 256 - sqrt(3) pi)/6) * sqrt(3) Gamma(1/3)^6/(4 * 2^(2/3) * pi^4).
      [Chen-507693-(2); MO 507693 eq. (2), Feb 2026, unanswered; OPEN; 80]
- D4  sum_{k>=0} [ 4(H_{4k} - H_k) - 4 log 396 + sqrt(58) pi ] C(4k,2k) C(2k,k)^2/396^(4k) = 0.
- D5  sum_{k>=0} [ 3(2 H_{6k} - H_{3k} - H_k) - 3 log 255 + 2 sqrt(7) pi ] C(6k,3k) C(3k,k) C(2k,k)/255^(3k) = 0.
      [Chen-507693-C1a, -C1b; instances of Chen's "Conjecture 1": for a Ramanujan series sum (ak+b) f(k),
      Re sum f'(k) = -2 Im(tau) pi sum f(k), tau as in the table on p. 7 of 2101.12592v1; unanswered; OPEN
      (probably accessible by modular methods a la Zhou 2503.05599, not obviously by WZ); |lhs| < 10^-99]

### 2.E  Guillera's L-value conjecture for the x-extended series (two explicit instances)

For a rational Ramanujan-like series sum_n prod_{i=0}^{2m} (s_i)_n/(1)_n * P(n) z0^n = sqrt((-1)^m chi)/pi^m
put S(x) = sum_{n>=0} prod_i (s_i)_{n+x}/(1)_{n+x} * P(n+x) * (-1)^n |z0|^{n+x} (for z0 < 0), with
(s)_y := Gamma(s+y)/Gamma(s). Conjecture (Guillera, arXiv:2312.16827 Sect. 6; MO 473746, June 2024,
unanswered): the coefficient of x^{2m+1} in S(x) is r * L(chi, m+1) with r rational. The rationals below
were identified here numerically (discrete Cauchy integral, 128 nodes; 76 digits):

- E1  s = (1/2, 1/4, 3/4), P(n) = 65n + 8, z0 = -(16/63)^2, chi = -7 (S(0) = 9 sqrt(7)/pi):
      [x^3] S(x) = (2205/2) L_{-7}(2).        [Guillera-L(-7,2); OPEN; 76 digits]
- E2  s = (1/2, 1/3, 2/3, 1/6, 5/6), P(n) = 5418n^2 + 693n + 29, z0 = -1/80^3, chi = 5 (S(0) = 128 sqrt(5)/pi^2):
      [x^5] S(x) = -5376000 L_5(3).           [Guillera-L(5,3); OPEN; 76 digits]
      (also observed: [x^1] = [x^3] = 0 and [x^2] S/S(0) = -(15/2) pi^2 for E2.)
  Written out, [x^{2m+1}] S is a series with the same hypergeometric part twisted by polynomials in
  psi^{(j)}(s_i + n) - i.e. a harmonic-number identity of the Sun/Au type; E1 is of the same nature as the
  Sun-Zhou theorem with L_{-7}(2) (2602.12091) and may follow from modular methods; E2 (a 1/pi^2 series)
  has no known approach.

### 2.F  Further precisely stated open items (not verified here or not single hypergeometric identities)

- Sun's Conjecture 5.45, first two identities (zeta(5); upside-down analogue of Cullen's series) — not
  reached by Shalev 2610.04770 (survey B).
- Sun's (4.25) [Le Lan 2609.22304, Conjecture 1] and Sun's supercongruence [Le Lan 2604.03327, Conjecture 1].
- Aldawoud's level-3 quadratic-irrational series N = 13, N = 59 (Campbell 2310.05112, Sect. 4; modular).
- p-adic companions of the 1/pi^3, 1/pi^4 series (Au 2403.04555, Sect. 1), Guillera's "alternative
  supercongruences" (MO 473848), Au's weight-6 supercongruence (2509.19960).
- Unanswered MathOverflow questions by Z.-W. Sun containing conjectural series (for survey B):
  436279 (zeta(5), second-order harmonic numbers), 456433 (sqrt(3)/pi), 456444 (L(2,(-3/.))), 456779
  (sqrt(3) pi, rate 1/9), 486518 (pi^10, rate 1/64), 487145 (pi^2), 487361 (rate -1/27 bunch), 492147,
  495096, 495403 (irrational series), 507600 (three harmonic-number series).

Count: 12 (Ex. VI) + 19 (Ex. VII) + 5 (Ex. VIII) + 5 (Ex. IX) + 1 (Ex. X) + 1 (S41) + 5 (B) + 1 (C) + 5 (D)
+ 2 (E) = 56 verified items. Nine carry "OPEN?" ((4.20), (4.28), VIII.1, B2-B5, D1, D2) and two more are
Sun-list items (S41, (4.30)); the remaining 45 are, as far as this sweep can tell, open identities with
hypergeometric summands (mostly twisted by harmonic numbers) outside Sun's lists: 38 from Au's paper, B1,
C1, D3-D5, E1-E2.

---------------------------------------------------------------------------------------------------

## 3. Au's WZ-seed method: technical summary

Source: K. C. Au, "Wilf-Zeilberger seeds and non-trivial hypergeometric identities", arXiv:2312.14051v4
(J. Symbolic Comput. 130 (2025) 102421), and the recap in arXiv:2602.08721v4 Sect. 2. Formulas are copied
from the text layer / TeX source. (The coordinator is reading 2312.14051 in full; this is kept short.)

### 3.1 Definition and criterion

Definition 3.1. Let f(a_1,...,a_m,k) be a hypergeometric term in a_i and k. It is a WZ-seed in the variables
a_i if for all A_i in Z and K in Z,
    F(n,k) = f(a_1 + A_1 n, ..., a_m + A_m n, k + K n)
has a hypergeometric WZ-mate G(n,k), i.e. F(n+1,k) - F(n,k) = G(n,k+1) - G(n,k).

Theorem 3.3 (criterion). f is a WZ-seed iff for each i there is a rational function r_i with
Delta_{a_i}(f) = Delta_k(r_i f) (forward differences). (Proof: (S_a1^A1 ... S_am^Am S_k^K - 1) lies in the
ideal (S_a1 - 1, ..., S_am - 1, S_k - 1) of the Laurent ring of shift operators.) This also gives a faster
algorithm for the mate than Gosper when A_i, K are large (implemented in WZSeedTools).

Heuristic source ("Problem 3.2"): if sum_{k>=0} f(a_1,...,a_m,k) exists and is independent of the a_i (i.e. f
is a classical summation formula divided by its Gamma-product value), then f is a seed — true in every
computed case (Gessel 1995 used this implicitly).

Example. Gauss2F1(a,b,c,k) has r_1 = k(c+k-1)/(a(a-c+1)), r_2 = k(c+k-1)/(b(b-c+1)), r_3 = -k/(a+b-c)   (3.1).

### 3.2 The central proposition (Proposition 2.1 of 2312.14051 = Proposition 2.2 of 2602.08721)

Suppose F, G: N^2 -> C satisfy F(n+1,k) - F(n,k) = G(n,k+1) - G(n,k) and
 * sum_{k>=0} F(0,k) converges,
 * sum_{n>=0} G(n,0) converges,
 * lim_{k->oo} G(n,k) = g(n) exists for each n in N and sum_{n>=0} g(n) converges.
Then lim_{n->oo} sum_{k>=0} F(n,k) exists and is finite, and
    sum_{k>=0} F(0,k) + sum_{n>=0} g(n) = sum_{n>=0} G(n,0) + lim_{n->oo} sum_{k>=0} F(n,k).
(Proof: sum the WZ equation over 0 <= k < K, 0 <= n < N, let K -> oo, then N -> oo.)

Proposition 3.4 (vanishing of the limit). For F(n,k) = A^n B^k prod_i Gamma(a_i + b_i n + c_i k)
/ prod_i Gamma(d_i + e_i n + f_i k) with sum b_i = sum e_i, sum c_i = sum f_i, put p(x) = |x|^x (p(0) = 1) and
    E_F(x) = |A| |B|^x prod_i p(b_i + c_i x) / prod_i p(e_i + f_i x),   x >= 0.
Then |F(n, xn)| = E_F(x)^n O(n^C), and lim_n sum_k F(n,k) = 0 if sup_{x>=0} E_F(x) < 1.

Search procedure (Sect. 3.4). For a target sum_n z^n (a_1)_n...(a_m)_n/((b_1)_n...(b_m)_n) R(n): for each seed and
integer shifts (A,B,...,K) in [-N,N], compute the "exponent" of F(n,0) (e.g. for Seed3:
B^-B 2^2D D^-D K^-K (A-B)^(A-B) (A+K)^(A+K) (B+2K)^(B+2K) (D+K)^(D+K) (2A-B-2D)^(2A-B-2D) (A-B-D)^(-A+B+D)
(2A-B+2K)^(-2A+B-2K) (A-D+K)^(-A+D-K), with 0^0 = 1), match it with z, then choose the offsets
(a,b,...,k0) in (1/l)Z/Z so that the Pochhammer part matches. R(n) need not be matched: it "pops up" from
Gosper. Free parameters a, b, ... are kept; coefficient extraction at a = b = ... = 0 (Laurent/Taylor) gives
the target and its harmonic-number twists. Limitation: only Seed1, Seed4, Seed10 give z outside the group
generated by {n^n}; hence e.g. (-1/80^3)^n (5418n+263) and the 7^-4n series B1 are not reachable.

Loose ends listed by Au: (i) offsets are only determined mod Z (trial and error; sometimes G(n,0) differs
from the target by a telescoping series, Ex. VII, VIII); (ii) g(n) may be non-zero (then recurse on g);
(iii) sum_k F(0,k) may need hypergeometric transformations (Ex. V, VI); (iv) lim_n sum_k F(n,k) may be
non-zero (Sect. 5, Ex. XV-XVII, asymptotic analysis). Guillera's "flawless WZ pairs" (2503.00570) remove
(ii)/(iv) for 1/pi^k formulas by replacing G(n,k) with G(n-k,k) (or F(n,k) with F(n,k-n)).

### 3.3 List of seeds (each is "summand / Gamma-product value" of a classical summation formula)

Calibration (my numerical check, `scripts/seedsum.py`, generic real parameters, non-terminating sums over
k >= 0): Gauss2F1 -> 1; Dixon3F2 -> 1/2; Watson3F2 -> 1/2; Dougall5F4 -> 1; Seed1 -> 1/(4 sqrt(pi));
Seed3 -> 1; Seed4 -> 1. So the normalisation is "independent of the parameters", not always "= 1".
Dougall7F6, Balanced3F2, Seed2, Seed7, Seed9, Seed10 come from terminating formulas (the factors (-1)^(a+e)
etc. are the analytic continuation in the terminating parameter) and were not tested this way.

(1) Gauss2F1(a,b,c,k) = Gamma(c-a) Gamma(a+k) Gamma(c-b) Gamma(b+k)
        / [Gamma(a) Gamma(b) Gamma(k+1) Gamma(c+k) Gamma(-a-b+c)].                    (Gauss 2F1)
(2) Dixon3F2(a,b,c,k) = Gamma(a-b+1) Gamma(a-c+1) Gamma(2a+k) Gamma(b+k) Gamma(c+k) Gamma(2a-b-c+1)
        / [Gamma(a) Gamma(b) Gamma(c) Gamma(k+1) Gamma(a-b-c+1) Gamma(2a-b+k+1) Gamma(2a-c+k+1)].   (Dixon 3F2)
(3) Watson3F2(a,b,c,k) = 2^(-2a-2b+2c) Gamma(-a+c+1/2) Gamma(2a+k) Gamma(-b+c+1/2) Gamma(2b+k) Gamma(c+k)
        / [Gamma(a) Gamma(b) Gamma(k+1) Gamma(2c+k) Gamma(-a-b+c+1/2) Gamma(a+b+k+1/2)].      (Watson 3F2)
(4) Dougall5F4(a,b,c,d,k) = (a+2k) Gamma(a+k) Gamma(b+k) Gamma(c+k) Gamma(d+k) Gamma(a-b-c+1) Gamma(a-b-d+1) Gamma(a-c-d+1)
        / [Gamma(b) Gamma(c) Gamma(d) Gamma(k+1) Gamma(a-b+k+1) Gamma(a-c+k+1) Gamma(a-d+k+1) Gamma(a-b-c-d+1)].
                                                                              (very-well-poised 5F4)
(5) Dougall7F6(a,b,c,d,e,k) = (-1)^(a+e) (a+2k) Gamma(a+k) Gamma(b+k) Gamma(c+k) Gamma(d+k) Gamma(e+k)
        Gamma(a-b-c+1) Gamma(a-b-d+1) Gamma(a-c-d+1) Gamma(-a+b+c+e) Gamma(-a+b+d+e) Gamma(-a+c+d+e)
        Gamma(2a-b-c-d-e+k+1)
        / [Gamma(b) Gamma(c) Gamma(d) Gamma(e) Gamma(k+1) Gamma(-a+b+e) Gamma(a-b+k+1) Gamma(-a+c+e) Gamma(a-c+k+1)
           Gamma(-a+d+e) Gamma(a-d+k+1) Gamma(a-e+k+1) Gamma(a-b-c-d+1) Gamma(2a-b-c-d-e+1) Gamma(-a+b+c+d+e+k)].
                                                  (very-well-poised 2-balanced terminating 7F6)
(6) Balanced3F2(a,b,c,d,k) = Gamma(d-a) Gamma(a+k) Gamma(d-b) Gamma(b+k) Gamma(d-c) Gamma(c+k) (-1)^(a+b+c-d)
        / [Gamma(a) Gamma(b) Gamma(c) Gamma(k+1) Gamma(d+k) Gamma(-a-b+d) Gamma(-a-c+d) Gamma(-b-c+d)
           Gamma(a+b+c-d+k+1)].                          (1-balanced terminating 3F2, Saalschutz)
(7) Dixon6F5(a,b,c,k) [only in 2602.08721] = 2^(4a-2b-2c) (-1)^(b+c+k) Gamma(a+1/2) Gamma(a-b+1) Gamma(a-c+1)
        Gamma(-b+k+1) Gamma(b+k) Gamma(-c+k+1) Gamma(c+k) (2a-b-c+2k+1) Gamma(a-b-c+3/2) Gamma(2a-b-c+1)
        Gamma(2a-b-c+k+1)
        / [Gamma(k+1) Gamma(2a-b+k+1) Gamma(2a-c+k+1) Gamma(2a-b-2c+k+2) Gamma(2a-2b-c+k+2)].
        (Origin: Dixon 3F2 + Whipple's 3F2(1) -> 6F5(-1). Not normalised to a constant: numerically
        sum_{k>=0} Dixon6F5(a,b,c,k)/(-1)^(b+c) = pi^3/(2 sin(pi b) sin(pi c)) (three test points, 25 digits),
        i.e. a 1-periodic function of the parameters, which is all the seed property needs.)
(8) Seed1(a,c,k) = 2^(a-2c-2k) Gamma(-a+2c-1) Gamma(a+2k) / [Gamma(a) Gamma(k+1) Gamma(-a+c-1/2) Gamma(c+k)].
        (a -> a/2, b -> a/2+1/2 in Gauss2F1)
(9) Seed2(a,c,d,k) = 2^(-2c-2k) (-1)^(a+c-d) Gamma(-a+2d-1) Gamma(a+2k) Gamma(d-c) Gamma(c+k)
        / [Gamma(a) Gamma(c) Gamma(k+1) Gamma(-a+d-1/2) Gamma(d+k) Gamma(-a-2c+2d-1) Gamma(a+c-d+k+3/2)].
        (a -> a/2, b -> a/2+1/2 in Balanced3F2; the last argument is printed "+2/3" in the PDF text layer,
        which is the fraction 3/2 with digits swapped by pdftotext)
(10) Seed3(a,b,d,k) = 2^(2d) (a+2k) Gamma(a-b+1/2) Gamma(a+k) Gamma(b+2k) Gamma(d+k) Gamma(2a-b-2d+1)
        / [Gamma(b) Gamma(d) Gamma(k+1) Gamma(a-b-d+1/2) Gamma(2a-b+2k+1) Gamma(a-d+k+1)].
        (b -> b/2, c -> b/2+1/2 in Dougall5F4)
(11) Seed4(a,b,k) = 3^(b+1) (a+2k) Gamma(3a-2b) Gamma(a+k) Gamma(b+3k)
        / [Gamma(b) Gamma(k+1) Gamma(a-b) Gamma(3a-b+3k+1)].
        (b -> b/3, c -> b/3+1/3, d -> b/3+2/3 in Dougall5F4)
(12) Seed7(a,b,d,e,k) = (-1)^(a+e) (a+2k) Gamma(a-b+1/2) Gamma(a+k) Gamma(b+2k) Gamma(d+k) Gamma(e+k)
        Gamma(2a-b-2d+1) Gamma(-a+b+e+1/2) Gamma(-2a+b+2d+2e) Gamma(2a-b-d-e+k+1/2)
        / [Gamma(b) Gamma(d) Gamma(e) Gamma(k+1) Gamma(a-b-d+1/2) Gamma(-2a+b+2e) Gamma(2a-b+2k+1) Gamma(-a+d+e)
           Gamma(a-d+k+1) Gamma(a-e+k+1) Gamma(2a-b-d-e+1/2) Gamma(-a+b+d+e+k+1/2)].
        (b -> b/2, c -> b/2+1/2 in Dougall7F6)
(13) Seed9(a,b,d,k) = (a+2k) (-1)^d 2^(4a-2b-2d) Gamma(a-b+1/2) Gamma(a+k) Gamma(b+2k) Gamma(d+2k)
        Gamma(-2a+2b+d+1) Gamma(-2a+b+2d+1) Gamma(2a-b-d+k)
        / [Gamma(b) Gamma(d) Gamma(k+1) Gamma(-a+d+1/2) Gamma(-2a+b+d) Gamma(2a-b+2k+1) Gamma(2a-d+2k+1)
           Gamma(-a+b+d+k+1)].
        (b -> b/2, c -> b/2+1/2, d -> d/2, e -> d/2+1/2 in Dougall7F6)
(14) Seed10(a,d,k) = (-1)^(a-d) 3^(-a-3k) Gamma(-a+3d-2) Gamma(a+3k)
        / [Gamma(a) Gamma(k+1) Gamma(-2a+3d-3) Gamma(d+k) Gamma(a-d+k+2)].
        (a -> a/3, b -> a/3+1/3, c -> a/3+2/3 in Balanced3F2)
Kummer, Whipple and Bailey 2F1/3F2 sums do not appear as separate seeds in Au's list (Watson is (3);
Seed5, Seed6, Seed8 are not listed in either paper — "many seeds which do not generate useful examples are
not included"). All seeds have q-analogues (2403.04555).

### 3.4 Hou-Mu extension (arXiv:2604.22377, Theorem 2.1) — seven seeds with one extra parameter

Notation Gamma[u_1,...; v_1,...] = prod Gamma(u_i)/prod Gamma(v_j).
(1) F0(k,a,b) = z^(b+k) (z-1)^(a-b) Gamma[a+k; b+k] Gamma(-a+b)                       [from (a)_k/(b)_k z^k]
(2) F0(k,a,b,c,d) = (-1)^(c+d) Gamma[a+k, b+k; c+k, d+k] Gamma[a+b-c-d+2; a-c+1, a-d+1, b-c+1, b-d+1]
    [bilateral Dougall 2H2; d = 1 gives Gauss2F1]
(3) F0(k,a,b,c) = (-1)^k Gamma[2a+k, 2b+k; 2a+2c+k, 2b+2c+k] Gamma[2c, a-b+c+1/2; a-b-c+1/2]
(4) F0(k,a,b,c,d,e) = (-1)^(d+e) Gamma[a+k, b+k, c+k; d+k, e+k, a+b+c-d-e+2+k]
      * Gamma[a+b-d-e+2, a+c-d-e+2, b+c-d-e+2; a-d+1, a-e+1, b-d+1, b-e+1, c-d+1, c-e+1]
    [e = a+b+c+1-d gives Balanced3F2]
(5) F0(k,a,b,c,d) = (-1)^(a+b+c+d) Gamma[a+k, b+k, c+k; 2d-a+k, 2d-b+k, 2d-c+k]
      * Gamma[a+b+c-3d+3/2; a+b-2d+1, a+c-2d+1, b+c-2d+1] / Gamma[a-d+1/2, b-d+1/2, c-d+1/2]
    [d = a+1 gives Dixon3F2]
(6) F0(k,a,b,c,d,e) = (-1)^(a+b+c+d+e) (e+2k) Gamma[a+k, b+k, c+k, d+k; e-a+1+k, e-b+1+k, e-c+1+k, e-d+1+k]
      * Gamma[a+b+c+d-2e; a+b-e, a+c-e, a+d-e, b+c-e, b+d-e, c+d-e]
    [e = a gives Dougall5F4]
(7) F0(k,a,b,c,d,e,f) = (-1)^f (f+2k) Gamma[a+k, b+k, c+k; f-a+1+k, f-b+1+k, f-c+1+k]
      * Gamma[d+k, e+k, 3f-a-b-c-d-e+1+k; f-d+1+k, f-e+1+k, a+b+c+d+e-2f+k]
      * Gamma[a+b+c+d-2f, a+b+c+e-2f, a+b+d+e-2f; a+b-f, a+c-f, a+d-f, a+e-f, b+c-f, b+d-f]
      * Gamma[a+c+d+e-2f, b+c+d+e-2f; b+e-f, c+d-f, c+e-f, d+e-f]
    [f = a gives Dougall7F6]
Theorem 2.2 there: if F(n,k) = F0(Kn+k, An+a, Bn+b, ...) has mate G (not rational in k), then G(k,n) is
again a WZ seed in the parameters n, a, b, ....

### 3.5 Worked examples (seed, shifts, resulting identity)

From 2312.14051 (the seed/shift data for I, II are stated on p. 10; the identification for III is mine, by
matching Gamma factors with the displayed F):

1. Example I (Cullen, 1/pi^4). Dougall7F6 with A = B = C = D = E = -1, K = 2 and a = b = c = d = e = 1/2,
   k0 = 1; displayed as
   F(n,k) = (2k+3n+1/2) Gamma(1-2n)^3 Gamma(n+1/2)^3 Gamma(k+n+1/2)^5 Gamma(k+4n)
            / [Gamma(1/2-n)^7 Gamma(2n)^2 Gamma(k-n+3/2) Gamma(k+2n+1)^5],
   G(n,k) = F(n,k) (-2k+2n-1)(68728320 k^2 n^8 + ... 61 monomials)/(131072 n^5 (k+2n+1)^5 (4k+6n+1)).
   lim_n sum_k F = 0, g != 0. Shift n -> n+a, divide by the Gamma prefactor, compare the a^-5 Laurent
   coefficient at a = 0:
   (1/(4096 pi^2)) sum_n (-1)^n (1/2)_n^10/(1)_n^10 ... reduces to Guillera's
   sum (1/2)_n^5/(1)_n^5 (-1)^n 2^-10n (820n^2+180n+13) = 128/pi^2, giving
   sum_{n>=0} 2^-12n (1/2)_n^7 (1/4)_n (3/4)_n/(1)_n^9 (43680n^4+20632n^3+4340n^2+466n+21) = 2048/pi^4.
2. Example II (Zhao, 1/pi^4). Dougall7F6 with A = 1, B = C = D = E = 0, K = 1, a = ... = e = 1/2, k0 = 1:
   F(n,k) = (-1)^n (1/2+2k+3n) Gamma(1-n)^3 Gamma(n+1/2)^3 Gamma(k+n+1/2)^4 Gamma(k+2n+1/2) Gamma(k+3n)
            / [Gamma(k+3/2) Gamma(1/2-n)^3 Gamma(n) Gamma(2n) Gamma(k+n+1) Gamma(k+2n+1)^4];
   a^-5 coefficient + Guillera's sum (1/2)_n^5/(1)_n^5 (-1)^n 2^-2n (20n^2+8n+1) = 8/pi^2 give
   sum_{n>=0} (-27/256)^n (1/2)_n^5 (1/4)_n (3/4)_n (1/3)_n (2/3)_n/(1)_n^9 (4528n^4+3180n^3+972n^2+147n+9) = 768/pi^4.
3. Example III (Sun; MO 456443). F(n,k) = Gauss2F1(a-n, b-n, c+1, k+2n+1/2), i.e.
   F = Gamma(-a+c+n+1) Gamma(a+k+n+1/2) Gamma(-b+c+n+1) Gamma(b+k+n+1/2)
       / [Gamma(a-n) Gamma(b-n) Gamma(k+2n+3/2) Gamma(c+k+2n+3/2) Gamma(-a-b+c+2n+1)];
   g = 0, limit 0. With a_n = 2^-6n (1/2)_n (1)_n^3/((1/4)_n^2 (3/4)_n^2), p_n = 42n^2-23n+3:
   sum_{n>=1} a_n p_n/(n^3 (2n-1)) = pi^2/2, plus four harmonic twists (coefficients of a, b, c).
4. Example IV (Guillera). Seed7; displayed
   F(n,k) = n (2k+2n+3/2) Gamma(n+1)^3 Gamma(k+n+1)^3 Gamma(2k+n+1) Gamma(k+2n+1)
            / [Gamma(k+3/2) Gamma(n+1/2)^3 Gamma(k+n+3/2)^3 Gamma(2k+3n+3)]:
   sum_{n>=1} (4/27)^n (1)_n^7/((1/2)_n^5 (1/3)_n (2/3)_n) (92n^3-84n^2+27n-3)/n^7 = 8 pi^4.
   (The paper's display has "21n"; numerically the coefficient must be 27, in agreement with table #20 of
   2101.12592 and with Sun's comment on MO 486961. Checked to 30 digits.)
5. Example VII (Guillera, 1/pi^2). Seed2;
   F(n,k) = 2^(-2k-6n) Gamma(-n-1) Gamma(n-1/2) Gamma(k+3n+2) Gamma(2k+3n+3/2)
            / [Gamma(-n-5/2) Gamma(-n-1/2) Gamma(n) Gamma(n+1) Gamma(k+2n+1) Gamma(k+2n+2) Gamma(k+2n+3)];
   g = 0, limit 0, sum_k F(0,k+a) = sum_{n>=1} G(n-1,a):
   sum_{n>=0} (-3^6/2^12)^n (1/2)_n (1/3)_n (2/3)_n (1/6)_n (5/6)_n/(1)_n^5 (1930n^2+549n+45) = 384/pi^2.
6. Example IX (Guillera, zeta(3)). F(n,k) = (-4)^n (4k+6n+3) Gamma(n+1)^2 Gamma(2n+1) Gamma(2k+n+1) Gamma(k+2n+1)^2
            / [Gamma(n+1/2)^2 Gamma(k+n+3/2)^2 Gamma(2k+5n+3)] (g = 0, limit 0):
   sum_{n>=1} (-2^8/5^5)^n (1)_n^5/((1/2)_n (1/5)_n (2/5)_n (3/5)_n (4/5)_n) (-483n^2+245n-30)/n^5 = 896 zeta(3).
7. Example X (Zhao, zeta(5)). Seed7; F(n,k) = (-1)^n (4k+6n+3) Gamma(n+1)^4 Gamma(2n+1) Gamma(2k+n+1) Gamma(k+2n+1)^4
            / [Gamma(n+1/2)^6 Gamma(k+n+3/2)^4 Gamma(2k+5n+3)]; g != 0 (needs two further evaluations):
   sum_{n>=1} (-2^10/5^5)^n (1)_n^9/((1/2)_n^5 (1/5)_n (2/5)_n (3/5)_n (4/5)_n) (5532n^4-5600n^3+2275n^2-425n+30)/n^9
   = -380928 zeta(5).   In 2602.08721 the same pair is written
   F(n,k) = Seed7(a+n+1/2, b-n, c+n+1/2, d+n+1/2, e+k+n+1/2).

From 2602.08721 (explicit seed calls):

8. Example I. F(n,k) = Dougall5F4(a, b, c, 1+d-n, k+e+n); sum_{n>=1} G(n-1,0) = sum_{k>=0} F(0,k) reads
   sum_{n>=1} (-1)^(n-1) (1-d)_{n-1} (a+e+1)_{n-1} (b+e+1)_{n-1} (c+e+1)_{n-1} (a-b-d+1)_{n-1} (a-c-d+1)_{n-1} P(n)
              / [(e+1)_n (a-b+e+1)_n (a-c+e+1)_n (a-d+e+1)_{2n} (a-b-c-d+1)_n]
   = sum_{k>=0} (a+2e+2k+2) (a+e+1)_k (b+e+1)_k (c+e+1)_k (d+e+1)_k
              / [(e+1)_{k+1} (a-b+e+1)_{k+1} (a-c+e+1)_{k+1} (a-d+e+1)_{k+1}],
   P(n) = a^3 - a^2(b+c+2d-3e) + a(b(c+d-2e) + cd - 2ce + d^2 - 5de + 3e^2)
          + n(5a^2 - a(3b+3c+7d-11e) + b(c+2d-4e) + 2cd - 4ce + 2d^2 - 8de + 5e^2)
          + n^2(9a - 3(b+c+2d-3e)) + e(b(c+2d-e) + 2cd - ce + 2d^2 - 2de + e^2) + 5n^3.
   Constant term: sum_{n>=1} (-1/4)^n (1)_n/(1/2)_n * 5/n^3 = -2 zeta(3) (the Apery series
   (5/2) sum (-1)^(n-1)/(n^3 C(2n,n)) = zeta(3)).
9. Example II. F(n,k) = Dougall5F4(a, b-n+1, c, d+1/2, k+n).
10. Example III ((-1/2^10)^n (1)_n^5/(1/2)_n^5, the Amdeberhan-Zeilberger series
    sum_{n>=1} (-1/2^10)^n (1)_n^5/(1/2)_n^5 (205n^2-160n+32)/n^5 = -2 zeta(3)): uses the known five-parameter
    extension (Au 2212.02986 Ex. VI; Chu-Zhang 2014), whose k-side is
    sum_{k>=0} (a+e+2k+2) (a+1)_k (b+1)_k (c+1)_k (d+1)_k/[(e+1)_{k+1} (a-b+e+1)_{k+1} (a-c+e+1)_{k+1} (a-d+e+1)_{k+1}]
    (the seed call is not printed; P(n) is in the Mathematica file). The "generators" of degrees 1-5 evaluate
    to -pi^4/60, -2 zeta(5), -2 zeta(3)^2, -zeta(7) and the zeta(5,3) series
    sum a_n [ (351-1755n+2240n^2)/n^10 - (11 H_n^{(5)} - H_{2n}^{(5)}) (205n^2-160n+32)/n^5 ]
    = -4 zeta(5,3) - pi^8/4725 (all proved there).
11. Example IV ((-1/27)^n (1)_n^2/((1/3)_n (2/3)_n)): F1 = Dougall5F4(a+n+1/2, b+1/2, c, d-n, k+d+n+1/2),
    F2 = Dougall5F4(a, b, c, d-n, e+k+2n+1), F3 = Seed3(a, b-n, c, k+d+n+1). Base:
    sum_{n>=1} a_n (5-32n+56n^2)/(n^3 (2n-1)^2) = -4 zeta(3); all three needed for Sun's
    sum_{n>=1} a_n (17-136n+408n^2-640n^3+560n^4)/(n^5 (2n-1)^4) = 180 zeta(5) - 56 pi^2 zeta(3)/3.
12. Example V (Gourevitch family). F(n,k) = Seed9(1+a-2n, 1/2+b-n, 1/2+c-n, k+n) with Guillera's flawless-
    pair technique gives a three-parameter formula whose a = b = c = 0 case is
    sum_{n>=0} 2^-6n (1/2)_n^7/(1)_n^7 (168n^3+76n^2+14n+1) = 32/pi^3.
13. Examples VI-IX: seeds and shifts as quoted in Section 2.A.

---------------------------------------------------------------------------------------------------

## 4. Things that were open in older lists but are now PROVED — do not re-prove

Rational 1/pi^k and "upside-down" formulas (numbering of the Cohen-Guillera table, arXiv:2101.12592):
- #11 Gourevitch: sum 2^-6n (1/2)_n^7/(1)_n^7 (168n^3+76n^2+14n+1) = 32/pi^3 — Au 2312.14051 Ex. XI/App. B;
  simpler: Guillera 2503.00570 Sect. 4.5; parametric: Au 2602.08721 Ex. V.
- #13 Cullen 1/pi^4 (43680n^4+...) = 2048/pi^4 — Au Ex. I; Guillera 2503.00570 Sect. 4.1; Shalev 2610.04770.
- #12 Zhao 1/pi^4 (4528n^4+...) = 768/pi^4 — Au Ex. II; Shalev 2610.04770.
- #6 sum (-3^6/2^12)^n ... (1930n^2+549n+45) = 384/pi^2 — Au Ex. VII; Guillera 2503.00570 Sect. 4.2.
- #17 (upside-down) sum (-2^8/5^5)^n (1)_n^5/((1/2)_n(1/5)_n(2/5)_n(3/5)_n(4/5)_n) (483n^2-245n+30)/n^5
  = -896 zeta(3) — Au Ex. IX.
- #19 (upside-down) sum 2^-6n (1)_n^7/(1/2)_n^7 (21n^3-22n^2+8n-1)/n^7 = pi^4/8 — Au Ex. XI (4.2).
- #20 (upside-down) sum (4/27)^n (1)_n^7/((1/2)_n^5(1/3)_n(2/3)_n) (92n^3-84n^2+27n-3)/n^7 = 8 pi^4 — Au Ex. IV;
  and Sun's sum (4/27)^n (1/2)_n^7/((1/6)_n(5/6)_n(1)_n^5) (92n^3+54n^2+12n+1)/(6n+1) = 12/pi^2 — Au Ex. VIII.
- #21 (upside-down) Zhao's zeta(5) series = -380928 zeta(5) — Au Ex. X; self-contained with four parameters:
  Shalev 2610.34482 (also the companion 240/pi^2 series).

Harmonic-number / zeta-value identities:
- Au 2602.08721 Ex. X, identities 2-6 (zeta(7), pi^8, zeta(9), zeta(9), zeta(5)^2 + pi^10) — Shalev 2610.34482
  (still printed with "?=" in Au's v4).
- Sun's ten 1/pi^4- and pi^4-type harmonic conjectures (Conj. 69, 70(i) of 2210.07238; 5.47(i) and last two of
  5.45 of the 2021 book) — Shalev 2610.04770.
- Sun, MO 486353: sum (-1/27)^n (1)_n^2/((1/3)_n(2/3)_n) (17-136n+408n^2-640n^3+560n^4)/(n^5(2n-1)^4)
  = 180 zeta(5) - 56 pi^2 zeta(3)/3 — Au 2602.08721 (1.5).
- Sun's seven conjectures (4.9)-(4.15) of "New series involving binomial coefficients (I)" — Au 2602.08721
  ((3.15)-(3.21), (3.26), (3.32)). Sun's (181), (183), (185) of 2210.07238 — Au 2602.08721 Ex. V.
- Sun, MO 456443: sum 2^-6n (1/2)_n(1)_n^3/((1/4)_n^2(3/4)_n^2) (42n^2-23n+3)/(n^3(2n-1)) = pi^2/2 and four
  harmonic twists — Au 2312.14051 Ex. III. Further Sun conjectures are proved in Au 2312.14051 Ex. V, VI,
  XII-XIV, XVI, XVII (see the paper for the statements).
- Rate-1/4 family (1/4)^n (1)_n^3/(1/2)_n^3 with H_n, H_{2n}, H^{(2)}, H^{(3)} ((4.2)-(4.5), (4.9) of Au) —
  Li-Chu, Electron. Res. Arch. 31 (2023) 4611-4636.
- Rate-1/64 family: (4.22), (4.23) of Au — Li-Chu, Bull. Malays. Math. Sci. Soc. 47 (2024) 94; Sun's
  Conjecture 38 (H^{(3)}, 555 zeta(3)/(77 pi) - 32G/11) and the second-order L_{-7}(2), L_{-8}(2), L_{-3}(2),
  L_{-4}(2) identities — Sun-Zhou 2602.12091 (modular; "inaccessible to our WZ approach" per Au).
- Sun's Conjecture 29 ((log m)/pi twists of every classical Ramanujan series; includes Au (4.19), VIII.2,
  IX.2) — Zhou 2503.05599; Hou-He-Wang 2305.00498.
- Sun's conjectures proved by Wei / Wei-Xu / Wei-Ruan (2211.11484, 2303.05402, 2303.07887, 2304.09753,
  2306.02641, 2308.06440, 2310.04642), Hou-Sun (2310.03699, 2604.15172), Zhou (2306.04638), Sun-Zhou
  (2401.14197, 2506.01865, 2607.14646), Jiang-Sun (2601.18645), Campbell (2304.00360), Le Lan (2604.03327,
  2609.22304), Patel (2607.15303, 2607.27229), Au (2201.01676 Props 7.10-7.13, 2212.02986).
- Sun's 2011 1/pi conjectures with T_k(b,c)-type sums — van Hoeij-Tsai-Ye 2410.19289 (modular
  parametrisation).
- Deyi Chen's derivative sums sum_k f^{(m)}(k), f = 1/(k^2 C(2k,k)) (MO 508743) — Hou-Sun 2604.15172;
  MO 508825 — closed form by Beenakker. Deyi Chen MO 486961 (1), (2), (4): (1) is Guillera 2003/Au Ex. XI,
  (2) follows from Au's Cullen pair by n -> n+1/2 (Guillera's answer), (4) is in 2101.12592 p. 16 (Sun's comment).
- Guillera-Rogers 3 ln 2 series: WZ proof completed — Chen-Liu-Lu 2610.16844.
- Bagis-Glasser d = 163 series — Campbell 2402.08485; three Aldawoud level-3 conjectures — Campbell
  2310.05112; level-11 coefficient formula (Cooper et al. 2015) — Campbell 2610.00551.
- Borwein's sum (-1)^n/(n^4 C(2n,n)) and sum (-1)^n/(n^5 C(2n,n)) evaluations — Au 2201.01676 Ex. 7.6.
- L(eta(2tau)^12, 5) = (pi^6/1024) sum (4n+1) (1/2)_n^8/(1)_n^8 — Au 2509.19960 (Fourier-Legendre + L-value
  relations).
- Guillera-Zudilin supercongruences (172n^2+75n+9; 11n+3) — Qi 2604.04006. Ramanujan-Machine conjectures —
  Yamamoto 2403.09729, Shvets 2604.06239, Wang 2601.08461.

---------------------------------------------------------------------------------------------------

## Appendix: files

- `verify_C.py`, `verify_C.out.txt` — numerical verification (63 entries; 7 are calibration entries for proved
  identities).
- `api/sweep.json`, `api/sweep2.json`, `api/filtered.txt` — raw arXiv sweep; `mo/mo.json`, `mo/ans1.json`.
- `pdf/` — downloaded PDFs and text layers; `src/2602.08721/`, `src/2604.22377/` — TeX sources.
- `scripts/` — sweep, MO and test scripts (`guillera_L.py` finds the rationals of E1, E2).
