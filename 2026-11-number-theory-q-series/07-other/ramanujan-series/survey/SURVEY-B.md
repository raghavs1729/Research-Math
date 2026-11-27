# SURVEY-B — Zhi-Wei Sun's conjectural hypergeometric series (with/without harmonic numbers): what is still open

Survey date: 2026-11-24.  Companion files (same directory):
`verify_B.py` (mpmath re-computation of every identity listed as open; importable list `ENTRIES`),
`verify_B.out.txt` (its output: 89 entries, 0 failures, worst agreement 87 digits).

## 0. Scope, notation, method, caveats

**Notation (plain ASCII).** `C(n,k)` = binomial coefficient. `H_n^{(m)} = sum_{0<j<=n} 1/j^m`, `H_n = H_n^{(1)}`.
`zeta(s)` Riemann zeta; `G` = Catalan's constant = L_{-4}(2); `beta(4) = sum_{k>=0} (-1)^k/(2k+1)^4`;
`K = L_{-3}(2) = sum_{n>=0} (1/(3n+1)^2 - 1/(3n+2)^2)`; `L_{-3}(4)` likewise with 4th powers.
All sums run to infinity; the lower index is stated in each line (`k>=0` or `k>=1`).

**Sources read (TeX source of the latest arXiv version, downloaded 2026-11-24 unless noted).**

| tag | paper | version read |
|---|---|---|
| H | Sun, "Series with summands involving harmonic numbers", arXiv:2210.07238 | v9 (2025-03-03; = published version, Springer PROMS 464) |
| I | Sun, "New series involving binomial coefficients (I)" (Nanjing Univ. J. Math. Biquarterly 41 (2024)) | = arXiv:2307.03086**v6** (2023-11-16) |
| II | Sun, "New series involving binomial coefficients (II)", arXiv:2307.03086 | v10 (2026-02-06) |
| III | Sun, "New series involving binomial coefficients (III)", arXiv:2506.01870 | v6 (2026-02-10) |
| V | Sun, "Various conjectural series identities", arXiv:2603.29973 | v3 (2026-04-13) |
| HS1 | Hou–Sun, "Taylor coefficients and series involving harmonic numbers", arXiv:2310.03699 | v3 (2026-02-10) |
| HS2 | Hou–Sun, "Evaluations of some series via the WZ method", arXiv:2604.15172 | v1 (2026-04-16) |
| JS | Jiang–Sun, "On series identities involving C(4k,k) and harmonic numbers", arXiv:2601.18645 | v1 |
| SZ1 | Sun–Zhou, "Series involving central binomial coefficients and higher-order harmonic numbers", arXiv:2602.12091 | v2 (2026-03-02) |
| SZ2 | Sun–Zhou, "Multiple Clausen values and deformed Apery-like series", arXiv:2607.14646 | v2 (abstract/intro only) |
| SZ3 | Sun–Zhou, "Fast converging irrational series for L(2,(d/.))", arXiv:2506.01865 | v3 (abstract only) |
| Au-MZV | Au, "Multiple zeta values, WZ-pairs and infinite sums computations", arXiv:2212.02986 | v4 (grep only) |
| Au-seed | Au, "Wilf-Zeilberger seeds and non-trivial hypergeometric identities", arXiv:2312.14051 | v4 (grep only) |
| Au-H | Au, "Discovering hypergeometric series with harmonic numbers via Wilf-Zeilberger seeds", arXiv:2602.08721 | **v4 (2026-11-09)** — intro, Sec. 4 examples VI–IX and all remarks about Sun read |
| Z1, Z2 | Y. Zhou, "Sun's series via cyclotomic multiple zeta values" arXiv:2306.04638v3; "Notes on certain binomial harmonic sums of Sun's type" arXiv:2503.05599v2 | abstract + list of the equations of H they cite |
| W* | C. Wei (and Wei–Xu, Wei–Ruan): arXiv:2211.11484, 2303.05402, 2303.07887, 2304.09753, 2306.02641, 2308.06440, 2310.04642 | grep only |
| other | Cohen–Guillera arXiv:2101.12592; Zuniga arXiv:2601.12495, 2506.08245, 2506.10321; Campbell–Cooper–Ye arXiv:2602.09352; Campbell 2505.14392, 2502.15249; Hou–Mu 2604.22377; Shalev 2610.04770, 2610.34482; Au 2403.04555, 2610.17240, 2201.01676; Xu–Zhao 2203.04184; Le Lan 2604.03327 | grep only (abstracts read for Shalev, Hou–Mu, Campbell–Cooper–Ye) |
| MO | Sun's MathOverflow questions 2023-06 … 2026-03 (user 124654) via the StackExchange API, with all answers of questions 456443, 456444, 456528, 456779, 486353, 486518, 487033, 487145, 487361, 488056, 506695, 507592, 507600, 507679, 507709, 508663, 508768, 508769, 508774, 508786, 508844, 508845, 508850 | fetched 2026-11-24 |

**How status was determined.** (1) Remarks in Sun's own latest versions (II v10 and V v3 explicitly say which earlier
conjectures "remain open" / were "confirmed by ..."). (2) Au-H v4 (two weeks old) says for each family what is proved,
what is only obtained "?=" (i.e. without rigorous proof) and what is "inaccessible". (3) A local full-text corpus of
the ~45 TeX sources above was searched for the defining polynomial of each kernel in both orderings
(scripts `sg.py`, `ng.py` in this directory: e.g. `35k-8`, or the coefficient set `{145,104,18}` within a 90-character
window, which also catches Au's ascending Pochhammer-style notation). (4) arXiv API listings 2023-06 … 2026-11 for
authors Sun, Wei, Au, Y. Zhou, Campbell, Chu and for abstract keywords (harmonic numbers / series / conjecture / WZ /
Ramanujan-type). (5) Answers to Sun's MathOverflow questions.

**Caveats (honest limits).**
* "OPEN" below always means "open as far as I could determine on 2026-11-24 with the above corpus".
  Journal-only papers were *not* accessible (e.g. Li–Chu, Electron. Res. Arch. 31 (2023) 4611–4636; Li–Chu, Bull. Malays.
  Math. Sci. Soc. 47 (2024); Li–Chu, Axioms 2025; Wei–Ruan, JMAA 537 (2024); Hou–He–Wang, J. Difference Equ. Appl. 30
  (2024); C. Wei's 2024–2026 journal papers). I only know of them through citations in Au-H, SZ1 and Sun's remarks.
* Sun's book "New Conjectures in Number Theory and Combinatorics" (2021) and arXiv:1102.5649 / 1911.05456 / 1010.4298
  were not mined for hypergeometric candidates: their hypergeometric items are old and, as far as the remarks in H, I, II
  show, mostly settled (Guillera, Guillera–Rogers, Hessami Pilehrood, Zudilin, Au, Zhou, Ablinger, Xu–Zhao …). The open
  list below is therefore drawn from 2022–2026 material (H, I, II, III, V, MO).
* github.com was blocked, so Au's Mathematica package / notebook could not be consulted.
* The arXiv id 2209.01731 in the brief is unrelated; Au's papers are 2212.02986, 2312.14051, 2403.04555, 2602.08721.
* Two statements in V v3 are **numerically false as printed** (see B4, B5) and one has a misprinted "=" (B72); the
  corrected forms (found by PSLQ, verified to 95 digits) are given.

**Numerics.** Every identity in Sections 1–2 was summed directly (working precision 95 digits, harmonic numbers by
incremental accumulation) until six consecutive terms were below 1e-88, and compared with the closed form:
all 89 agree to >= 87 digits (see table at the end, produced from `verify_B.out.txt`).

---------------------------------------------------------------------------------------------------

## 1. OPEN — hypergeometric, no harmonic numbers

| id | identity | source | status / how checked |
|---|---|---|---|
| A1 | `sum_{k>=1} C(2k,k) 16^k (1656k^4-2604k^3+1158k^2-229k+15) / (k^2 (2k-1)^2 (6k-1)(6k-5) C(6k,3k) C(3k,k)) = 0` | MO 506695 (Sun, 2026-01-06) | OPEN. Only answer (C. Leibovici, 2026-01-17) is a partial-fraction remark, no proof. Coefficients `1656,2604` occur nowhere in the corpus. Rate 4/27. |
| A2 | `sum_{k>=1} (-1)^(k-1) (448k^4-104k^3+24k^2-8k+1) C(2k,k)^2 / ((2k-1)^4 k^3 C(6k,3k) C(3k,k)) = 20 zeta(3)` | MO 487361 eq.(3) (Sun, 2025-02-07); note 448k^4-…=(16k^2+2k+1)(28k^2-10k+1) | OPEN. MO question has no answer; HS2 (its first theorem, eq. `id-1`) proves only the companion `… (28k^2-8k+1) … = 2 log 2` (eq.(2) of the same MO post). `448,104,24` not in corpus. Rate -1/27. |
| A3 | `sum_{k>=1} (10k-3) / ((2k-1) k^2 3^k C(2k,k) C(3k,k)) = K/2` | II v10 Conj 5.1(i), eq.(5.1) "(81K)" [correction 2026-11-24: not in MO 456444, which is the `(42k^2-27k+4)(-9)^k` series proved by Au] | OPEN. MO: no answer. Zuniga arXiv:2601.12495v2 (Feb 2026) still quotes it as "conjectured identity" (and used it to compute L_{-3}(2)). Rate 1/81. |
| A4 | `sum_{k>=1} (4k-1) 3^k / ((2k-1) k^2 C(2k,k) C(3k,k)) = 2K` | II v10 Conj 5.1(i), 2nd eq. | OPEN (same check as A3; rate 1/9). |
| A5 | `sum_{k>=1} (8k-3) C(4k,2k) / (k (4k-1) 9^k C(2k,k)^2) = sqrt(3) pi/18` | II v10 Conj 5.9(i) eq.(8k-3); MO 456779 (2023-10-20) | OPEN. MO: no answer; no third-party hit. Rate 1/9. |
| A6 | `sum_{k>=1} (13k^2-15k+2) 9^(k-1) / (k (3k-1)(3k-2) C(4k,k)) = pi/sqrt(3)` | II v10 Conj 4.17(i) | OPEN. Sun (Remark after it): equivalent to `sum (221k^2-6411k+844) 9^(k-1)/(k C(4k,k)) = 243 - 20pi/sqrt3`, "unable to prove by the integration method". Single binomial coefficient; rate 243/256 (slow). |
| A7 | `sum_{k>=0} (360k^3+612k^2+230k+15) C(3k,k) C(6k,3k)^2 / ((3k+1)(3k+2) 2^(15k)) = 32 sqrt(2)/pi` | II v10 Conj 5.25(i), eq.(5.20) | OPEN. V v3 Remark 5.11 (Apr 2026) still says "was conjectured to have the value". No third-party hit. Rate 27/32. |
| A8 | `sum_{k>=0} (1344k^3+944k^2+156k+9) C(2k,k)^2 C(4k,2k)^3 / ((2k+1) 2^(22k)) = 64 sqrt(2)/pi^2` | II v10 Conj 5.33(i) (dated 2025-06-26), eq.(5.31) | OPEN. V v3 Remark 5.23 still "conjectured". No third-party hit. Rate 1/64. |
| A9 | `sum_{k>=1} (77k^3-86k^2+29k-3) (-32)^k / (k^3 (3k-1)(3k-2) C(3k,k) C(4k,k)^2) = -6G` | II v10 Conj 4.33(i) (2025-02-28) | OPEN. No third-party hit for `77,86,29`. Rate -27/512. |
| A10 | `sum_{k>=0} (3k+2)^2 16^k / ((2k+1)(6k+1)(6k+5) C(6k,3k)) = pi/(2 sqrt(3))` | III v6 Conj 3.2, eq.(16^k) | OPEN. Sun: equivalent to `sum 16^k(54k^2-36k-12)/((6k+1)C(6k,3k)) = -4 - pi/sqrt3`, "which we cannot prove via integration". Rate 1/4. |

Remarks.
* A3/A4/A5/A9 are of the "Zeilberger/Guillera" type `sum R(k) z^k / prod C(..)`; the natural attack is a WZ pair
  from one of Au's seeds (cf. the Chu–Zhang list, Math. Comp. 83 (2014), whose Examples Sun uses as starting points).
* III v6 Conj 3.6–3.9 (`sum (57k^2+2k+1) C(6k,3k)/216^k = 9 sqrt3`, `sum (6k-1) C(6k,3k)/256^k = -2/(3 sqrt3)`,
  `sum (84k^2-4k+1) C(6k,3k)/512^k = 4 sqrt2`, `sum (6k-1)^2 C(6k,3k)/(-512)^k = 4 sqrt6/27`) are stated as conjectures
  by Sun but are **not listed here**: the generating function `sum C(6k,3k) y^(3k) = (1/3) sum_{w^3=1} (1-4wy)^(-1/2)` is
  algebraic, so these are routine; not WZ targets. (Not verified in `verify_B.py`.)
* III v6 Conj 5.1–5.3 and 6.1–6.2 (26 irrational Ramanujan-type series for 1/pi and irrational "Zeilberger-type" series for
  L_d(2), via Sun's "Duality Principle") are hypergeometric without harmonic numbers, but have quadratic-irrational
  arguments; they belong to the modular world (Campbell–Cooper–Ye arXiv:2602.09352 classify conjugate quadratic-irrational
  1/pi series; SZ3 treats L_d(2) series by Guillera–Rogers theory). Status not determined; not WZ targets; not listed.
* Still-conjectural hypergeometric *base* identities quoted by Sun but **not due to Sun** (Cohen–Guillera, arXiv:2101.12592),
  as described in V v3 (Remarks 5.20–5.22, 5.26, 5.27: "was conjectured to have the value"):
  `sum (5418k^2+693k+29) C(2k,k)^2C(3k,k)^2C(6k,3k)/(-2^18 3^6 5^3)^k = 128 sqrt5/pi^2`,
  `sum (1640k^2+278k+15) C(2k,k)^2C(3k,k)C(4k,2k)C(6k,3k)/(-2^22 3^3)^k = 256/(sqrt3 pi^2)`,
  `sum (1920k^2+304k+15) C(2k,k)^3C(4k,2k)C(8k,4k)/(2^18 7^4)^k = 56 sqrt7/pi^2`,
  and the two 1/pi^4 series with kernels `C(2k,k)^7C(3k,k)C(4k,2k)/(-2^24)^k` and `C(2k,k)^8C(4k,2k)/2^(32k)`.
  (Au-seed has hits for `1920,304,15` and `43680,20632`, so some of these may already be proved there — I did not
  read those passages; treat as "status not determined".)

---------------------------------------------------------------------------------------------------

## 2. OPEN — hypergeometric with harmonic numbers (grouped by hypergeometric kernel)

Convention: "no hit" = the kernel polynomial (both orderings) does not occur in any third-party source of the corpus.

### K1. `C(2k,k)^2 C(3k,k) / m^k`, classical Ramanujan 1/pi kernels — V v3 Conj 3.1–3.3 = MO 507592 (2026-01-30), eqs (1)–(6)

Base series: `sum (5k+1) C(2k,k)^2C(3k,k)/(-192)^k = 4/(sqrt3 pi)`, `sum (51k+7)…/(-12)^(3k) = 12 sqrt3/pi`,
`sum (9k+1)…/(-8640)^k = 4 sqrt15/(5 pi)`.

* **B1** `sum_{k>=0} C(2k,k)^2C(3k,k)/(-192)^k * ((5k+1)(H_{2k}-H_k) + (3k+1)/(2k+1)) = 4 log3/(sqrt3 pi)`
* **B2** `sum_{k>=0} C(2k,k)^2C(3k,k)/(-192)^k * ((5k+1)(H_{3k}-H_k) + (4k+3)/(6k+3)) = (4/(3 sqrt3 pi)) log(64/3)`
* **B3** `sum_{k>=0} C(2k,k)^2C(3k,k)/(-12)^(3k) * ((51k+7)(H_{2k}+H_k) - 3(13k+7)/(2k+1)) = -36 sqrt3 log3/pi`
* **B4** `sum_{k>=0} C(2k,k)^2C(3k,k)/(-12)^(3k) * ((51k+7)(3H_{3k}-7H_k) + (180k+93)/(2k+1)) = 36 sqrt3 log(108)/pi`
  — **corrected**: V v3 and MO 507592 eq.(4) print `- (180k+93)/(2k+1)`; that version is false
  (lhs = -91.8375…, rhs = 92.9301…). PSLQ on {S, sum k/(2k+1)·kernel, sum 1/(2k+1)·kernel, sqrt3 log2/pi, sqrt3 log3/pi}
  gives the "+" sign; verified to 95 digits.
* **B5** `sum_{k>=0} C(2k,k)^2C(3k,k)/(-8640)^k * ((9k+1)(H_{2k}-2H_k) + (27k+13)/(6k+3)) = 48 log3/(sqrt15 pi)`
  — **corrected** in the same way (printed sign `-` is false: lhs = -4.3203…, rhs = 4.3340…).
* **B6** `sum_{k>=0} C(2k,k)^2C(3k,k)/(-8640)^k * ((9k+1)(3H_{3k}-H_k) + 1/(6k+3)) = (12/(sqrt15 pi)) log(320/243)`

Status: OPEN. MO 507592 has two answers, neither a proof (Nolord: heuristic "differentiate the binomial coefficient";
C. Beenakker: "I expect Cohen–Guillera to cover it, although I have not checked"). These are *not* instances of Sun's
General Conjecture (H Conj 29, proved by Zhou Z2 §3), which fixes the combination `3H_{3k}+2H_{2k}-5H_k`; here one needs a
second, independent relation per kernel. No hit for `51k+7`, `63k+8`, `342k+25` outside Cohen–Guillera's table.
WZ remark: Au-H Examples VIII/IX do exactly this for the kernels 48^2 and (-32)^3 (one seed, coefficient extraction in
2–3 auxiliary parameters) and list their analogues with "?=" — the same strategy should apply.

### K2. `C(2k,k)^2 C(4k,2k) / (-12288)^k` — V v3 Conj 3.4, 3.5 (MO 507592 (7),(8); MO 507600 (3))

* **B7** `sum_{k>=0} C(2k,k)^2C(4k,2k)/(-12288)^k * ((28k+3)(H_{2k}-H_k) + (20k+9)/(4k+2)) = 4 log(432)/(sqrt3 pi)`
* **B8** `sum_{k>=0} C(2k,k)^2C(4k,2k)/(-12288)^k * ((28k+3)(H_{4k}-H_{2k}) + (8k+5)/(4k+2)) = (8/(sqrt3 pi)) log(16/3)`
* **B9** `sum_{k>=0} C(2k,k)^2C(4k,2k)/(-12288)^k * ((28k+3)(36H_{4k}^{(2)} - 9H_{2k}^{(2)} - H_k^{(2)}) + 12/(4k+1)) = 16 pi/(3 sqrt3)`

Status: OPEN. MO 507600 (2026-01-30) has zero answers. B9 is the analogue of H Conj 41 (eq.(119), kernel (-1024)^k,
whose own status I could not determine).

### K3. `C(2k,k) C(3k,k) C(6k,3k) / m^k` — V v3 Conj 3.6, 3.7 (MO 507592 (9),(10))

* **B10** `sum_{k>=0} C(2k,k)C(3k,k)C(6k,3k)/(-15)^(3k) * ((63k+8)(3H_{2k}-5H_k) + (198k+83)/(4k+2)) = 30 sqrt15 log3/pi`
* **B11** `sum_{k>=0} C(2k,k)C(3k,k)C(6k,3k)/(-96)^(3k) * ((342k+25)(3H_{2k}-7H_k) + (1122k+571)/(2k+1)) = 64 sqrt6 log(2^7 3^6)/pi`

Status: OPEN (same check as K1).

### K4. `C(2k,k)^3 / 4096^k` (Ramanujan's 42k+5)

* **B12** `sum_{k>=0} C(2k,k)^3/4096^k * (9(42k+5)(H_{2k}^{(4)} - H_k^{(4)}/16) + 25/(2k+1)^3) = 5 pi^3/6`
  — I eq.(5.14) = II v10 Conj 5.13(i) (2023-06-18).

Status: OPEN. Au-H v4 Example VII lists exactly this identity (in the form `… = 5 pi^3/54`) with "?=", i.e. found by WZ
coefficient extraction but **without** rigorous evaluation of the boundary terms, and says "(4.30) was also conjectured
by Sun". SZ1 treats orders 2 and 3 only. (All other Sun conjectures on this kernel are proved, see Section 3.)
Au-H Example VII contains ~15 further "?=" identities on this kernel (e.g.
`sum a_n(-9/(2n+1)^2 + (1/3)(17H_{2n}^{(3)} - 2H_n^{(3)})(42n+5)) ?= 80 zeta(3)/pi - 128G/3`, which Au says would follow
from two Sun conjectures — both of those are now covered by SZ1, see Section 3).

### K5. `1/(k^3 C(2k,k)^3)` (Zeilberger's 21k-8)

* **B13** `sum_{k>=1} ((21k-8)(2H_{3k-1} - H_{k-1}) - 4(27k-11)/(3(3k-1))) / (k^3 C(2k,k)^3) = 3 zeta(3)`
  — V v3 Conj 3.8 = MO 507600 (1).

Status: OPEN (MO: no answers). Note the unusual `H_{3k-1}`: the known WZ results on this kernel (Au-MZV, after Thm 2.2 (per H Remark):
`sum ((21k-8)(H_{2k-1}-H_{k-1}) - 7/2)/(k^3C(2k,k)^3) = zeta(3)`; HS2 Thm on f_7: all derivative sums up to order 7)
only produce H_{k}, H_{2k}.

### K6. `z^k / (k^3 C(2k,k)^2 C(3k,k))`, z = 8, 64, -27

Base: `sum (10k-3)8^k/(…) = pi^2/2` (Guillera–Rogers), `sum (11k-3)64^k/(…) = 8pi^2` (Guillera), `sum (15k-4)(-27)^(k-1)/(…) = K` (Hessami Pilehrood).

* **B14** `sum_{k>=1} 8^k ((10k-3)(H_{2k-1}-H_{k-1}) - 1) / (k^3 C(2k,k)^2 C(3k,k)) = 7 zeta(3)/2` — H Conj 23(i), eq.(63)
* **B15** `sum_{k>=1} 8^k ((10k-3)(H_{3k-1}-H_{k-1}) - 8/3) / (k^3 C(2k,k)^2 C(3k,k)) = (2 pi^2 log2 + 7 zeta(3))/4` — H Conj 23(i)
* **B16** `sum_{k>=1} 8^k ((50k-15) H_{k-1}^{(2)} + 4/k) / (k^3 C(2k,k)^2 C(3k,k)) = pi^4/24` — II v10 Conj 5.15(i)
* **B17** `sum_{k>=1} 64^(k-1) ((11k-3)(2H_{2k-1}+H_{k-1}) - 4) / (k^3 C(2k,k)^2 C(3k,k)) = 7 zeta(3)/2` — H Conj 26(i), eq.(71)
* **B18** `sum_{k>=1} 64^(k-1) ((11k-3)(3H_{3k-1}-6H_{k-1}) - 7) / (k^3 C(2k,k)^2 C(3k,k)) = (6 pi^2 log2 - 21 zeta(3))/8` — H Conj 26(i)
* **B19** `sum_{k>=1} 64^k ((55k-15) H_{k-1}^{(2)} + 8/k) / (k^3 C(2k,k)^2 C(3k,k)) = 8 pi^4/3` — II v10 Conj 5.17(i)
* **B20** `sum_{k>=1} (-27)^k ((15k-4)(H_{2k-1}-H_{k-1}) + 1/(2k-1)) / (k^3 C(2k,k)^2 C(3k,k)) = -4 pi^3/(3 sqrt3)` — H Conj 25(i) (2023-08-26)
* **B21** `sum_{k>=1} (-27)^(k-1) ((15k-4)(9H_{3k-1}^{(2)} - 5H_{k-1}^{(2)}) - 18/(3k-1)) / (k^3 C(2k,k)^2 C(3k,k)) = 9 L_{-3}(4)` — V v3 Conj 3.9 = MO 507600 (2)

Status: OPEN. II v10 (Feb 2026), Remarks after Conj 5.15 and 5.17, say explicitly that (63) and (71) "remain open".
Au-H's kernels are (-1/27)^n(1)_n^2/((1/3)_n(2/3)_n), 2^-6(1/2)_n^3, (-2^-10)(1)_n^5/(1/2)_n^5, (-1/4)(1)_n/(1/2)_n,
2^-6(1/2)_n^7, 4^-n(1)_n^3/(1/2)_n^3, the 48^2 and the (-32)^3 kernel — none of the K6 kernels. HS1 (first theorem of its introduction) proves only the
companion `sum (-27)^k((15k-4)(3H_{3k-1}-H_{k-1}) - 9)/(…) = -4pi^3/sqrt3` (H Conj 24(i)) by the Pilehroods' WZ pair, which is
the obvious starting point for B20/B21. `50k-15`, `55k-15`, `11k-3`+harmonic: no hit.

### K7. `81^k / (k^3 C(2k,k)^2 C(4k,2k))`

Base: `sum (35k-8)81^k/(…) = 12 pi^2` (Guillera–Rogers).

* **B22** `sum_{k>=1} 81^k ((35k-8)(H_{4k-1}-H_{k-1}) - 35/4) / (k^3 C(2k,k)^2 C(4k,2k)) = 12 pi^2 log3 + 39 zeta(3)` — H Conj 27(i), eq.(75)
* **B23** `sum_{k>=1} 81^k ((350k-80) H_{k-1}^{(2)} + 27/k) / (k^3 C(2k,k)^2 C(4k,2k)) = 4 pi^4` — II v10 Conj 5.19(i)

Status: OPEN. II v10 Remark after Conj 5.19: (75) "remains open". `35k-8`, `350k-80`: no hit.
(H Remark after Conj 28: Sun could not even guess the value of the analogous (-144)^k, 5k-1 series.)

### K8. `1/(k^3 (2k-1) C(2k,k) C(3k,k)^2)` and `(-1)^k/(k^3 (2k-1) C(2k,k) C(4k,2k)^2)` — third-order harmonic numbers

* **B24** `sum_{k>=1} (17(145k^2-104k+18) H_{k-1}^{(3)} + 28(2k-1)/k^2) / (k^3 (2k-1) C(2k,k) C(3k,k)^2) = 528 zeta(5) - 46 pi^2 zeta(3)` — II v10 Conj 5.21(i)
* **B25** `sum_{k>=1} (-1)^k (59(410k^2-197k+24) H_{k-1}^{(3)} - 62(2k-1)/k^2) / (k^3 (2k-1) C(2k,k) C(4k,2k)^2) = 6(57 pi^2 zeta(3) - 652 zeta(5))` — II v10 Conj 5.23(i)

Status: OPEN as far as determined. Proved on these kernels: the base series and the four first/second-order identities
(Au, WZ; MO 456443 answer and II Remark 5.22), and all derivative sums `sum f^{(m)}(k)`, m<=7 (HS2, Thms for f_8, f_9), and
for 410k^2-… four first/second-order identities (HS1). The derivative sums of order 3 are *polynomials* in harmonic
numbers, so B24/B25 (pure H^{(3)}) do not follow directly; but they should be reachable by coefficient extraction from the
same seeds with more parameters.

### K9. `3^(-k)` resp. `3^k` over `(2k-1) k^2 C(2k,k) C(3k,k)` — II v10 Conj 5.1(i) (2023-10-15)

* **B26** `sum_{k>=1} ((10k-3)(3H_{3k-1}+4H_{2k-1}-6H_{k-1}) - 12) / ((2k-1) k^2 3^k C(2k,k) C(3k,k)) = 4 pi^3/(27 sqrt3)`
* **B27** `sum_{k>=1} 3^k ((4k-1)(3H_{3k-1}-3H_{2k-1}-H_{k-1}) + 6k/(2k-1)) / ((2k-1) k^2 C(2k,k) C(3k,k)) = 8 pi^3/(27 sqrt3)`

Status: OPEN (base identities A3, A4 are themselves open).

### K10. `(-16)^k / ((2k-1) k^2 C(2k,k) C(4k,2k))` — II v10 Conj 5.3(i), 5.5(i) (2023-08/10)

Base (Chu–Zhang Thm 9): `sum (-16)^k (5k-1)/((2k-1)k^2C(2k,k)C(4k,2k)) = -pi^2/2`.

* **B28** `sum_{k>=1} (-16)^k/(k^2 C(2k,k) C(4k,2k)) * ((5k-1)/(2k-1) H_{k-1} - 1/(8k)) = pi^2 log2 - 21 zeta(3)/4`
* **B29** `sum_{k>=1} (-16)^k/(k^2 C(2k,k) C(4k,2k)) * ((5k-1)/(2k-1) H_{2k-1} + 1/(12k)) = (pi^2/3) log2 - 35 zeta(3)/6`
* **B30** `sum_{k>=1} (-16)^k/((2k-1) k^2 C(2k,k) C(4k,2k)) * ((5k-1) H_{4k-1} - (46k-5)/(24k)) = (pi^2/6) log2 - 77 zeta(3)/12`
* **B31** `sum_{k>=1} (-16)^k ((5k-1) H_{k-1}^{(2)} - 4(4k-1)/(2k-1)^2) / ((2k-1) k^2 C(2k,k) C(4k,2k)) = pi^4/6`
* **B32** `sum_{k>=1} (-16)^k ((5k-1) H_{k-1}^{(3)} + 12(4k-1)/(2k-1)^3) / ((2k-1) k^2 C(2k,k) C(4k,2k)) = -4 pi^2 zeta(3)`

Status: OPEN as far as determined: still conjectures in II v10 (Feb 2026); `(5k-1)H` occurs only in Sun's own lists.
HS1 proves the analogues for the neighbouring kernels `(-16)^k(20k^2-8k+1)/((2k-1)^2k^3C(2k,k)C(4k,2k))` and
`16^k(6k-1)/((2k-1)k^2C(2k,k)C(4k,2k))` by exactly this WZ/Taylor-coefficient method, so these look very reachable.
Compare B71 (the same phenomenon, "log 2 present", which Au could not reach).

### K11. `(-256)^k / (k^2 (2k-1)(4k-1)(4k-3) C(3k,k) C(6k,3k))` — II v10 Conj 5.7 (2023-11-18)

P(k) = 344k^3-386k^2+115k-9, Q(k) = 1952k^3-1732k^2+315k+12; base (Chu–Zhang Ex. 44): `sum (-256)^k P(k)/(…) = -8pi^2`.

* **B33** `sum_{k>=1} (-256)^k (P(k)(H_{2k-1}-2H_{k-1}) - (4k-1)(58k^2+181k-66)/(26k)) / (k^2(2k-1)(4k-1)(4k-3) C(3k,k) C(6k,3k)) = (1680/13) zeta(3) - 32 pi^2 log2`
* **B34** `sum_{k>=1} (-256)^k (P(k)(2H_{6k-1}-H_{3k-1}-3H_{k-1}) - Q(k)/(13k)) / (k^2(2k-1)(4k-1)(4k-3) C(3k,k) C(6k,3k)) = (2464/13) zeta(3) - 64 pi^2 log2`
  (II v10 prints the lower limit k=0, an obvious misprint for k=1.)

Status: OPEN as far as determined (`344,386,115`: no hit). HS1 proves the analogues for (-16)^k, P=112k^3-116k^2+35k-3 and for
64^k, 368k^3-… on the same denominator shape.

### K12. `C(4k,2k) / (k (4k-1) 9^k C(2k,k)^2)` — II v10 Conj 5.9(i) (2023-10-16); MO 456779

* **B35** `sum_{k>=1} C(4k,2k) ((8k-3)(5H_{2k-1}-4H_{k-1}) - 6) / (k (4k-1) 9^k C(2k,k)^2) = 3K/2`
* **B36** `sum_{k>=1} C(4k,2k) (8k-3)(2H_{2k-1}^{(2)} - 5H_{k-1}^{(2)}) / (k (4k-1) 9^k C(2k,k)^2) = pi^3/(36 sqrt3)`

Status: OPEN (base identity A5 open; MO no answer).

### K13. `C(3k,k) C(4k,k)^2 / ((3k+1)(3k+2) m^k)` and the reciprocal kernel — II v10 Conj 4.25–4.33 (Feb 2025)

P1 = 88k^3+108k^2+36k+3 (m=1024, base = 8/pi, Chu–Zhang Ex. 110); P2 = 368k^3+400k^2+118k+9 (m=4096, 16/pi, Ex. 106);
P3 = 896k^3+992k^2+296k+21 (m=-2^14, 32/pi, Ex. 97); P4 = 74k^3-84k^2+29k-3 (reciprocal, 256^k, 6pi^2, Ex. 51);
P5 = 77k^3-86k^2+29k-3 (reciprocal, (-32)^k, base A9 open).

* **B37** `sum_{k>=0} C(3k,k)C(4k,k)^2 (P1(k)(2H_{2k}-3H_k) + 40k^2+56k+14) / ((3k+1)(3k+2) 1024^k) = 32 log2/pi`
* **B38** `sum_{k>=0} C(3k,k)C(4k,k)^2 (P1(k)(2H_{4k}-3H_{2k}) + 40k^2+32k+8) / ((3k+1)(3k+2) 1024^k) = 16 log2/pi`
* **B39** `sum_{k>=0} C(3k,k)C(4k,k)^2 (P2(k)(H_{2k}+H_k) - 128k^2-136k-31) / ((3k+1)(3k+2) 4096^k) = -64 log2/pi`
* **B40** `sum_{k>=0} C(3k,k)C(4k,k)^2 (P2(k)(2H_{4k}-5H_{2k}) + 192k^2+212k+46) / ((3k+1)(3k+2) 4096^k) = 96 log2/pi`
* **B41** `sum_{k>=0} C(3k,k)C(4k,k)^2 (P3(k)(4H_{4k}-4H_{2k}+3H_k) - 576k^2-544k-110) / ((3k+1)(3k+2) (-2^14)^k) = -256 log2/pi`
* **B42** `sum_{k>=1} 256^k (P4(k)(3H_{2k-1}-4H_{k-1}) - 57k^2+33k-5) / (k^3(3k-1)(3k-2) C(3k,k)C(4k,k)^2) = 24 pi^2 log2 - 84 zeta(3)`
* **B43** `sum_{k>=1} 256^k (P4(k)(4H_{4k-1}-7H_{2k-1}) + 51k^2-15k-1) / (k^3(3k-1)(3k-2) C(3k,k)C(4k,k)^2) = 24 pi^2 log2 - 84 zeta(3)`
* **B44** `sum_{k>=1} 256^k (P4(k)(4H_{4k-1}^{(2)} - H_{2k-1}^{(2)} - 4H_{k-1}^{(2)}) + 4k-4) / (k^3(3k-1)(3k-2) C(3k,k)C(4k,k)^2) = 2 pi^4`
* **B45** `sum_{k>=1} (-32)^(k-1) (P5(k)(4H_{4k-1}^{(2)} - H_{2k-1}^{(2)} - 4H_{k-1}^{(2)}) - 8k+4) / (k^3(3k-1)(3k-2) C(3k,k)C(4k,k)^2) = (3/4) beta(4)`

Status: OPEN as far as determined: `88,108,36`, `368,400,118`, `896,992,296`, `74,84,29`, `77,86,29` occur only in II and V
(V v3 Conj 5.8–5.10 re-states the first three as "sum g'(k) = 0 / 32i", i.e. still conjectural in Apr 2026).
(B42 and B43 have the same right-hand side, so their difference is a rational-function-free relation worth checking first.)

### K15. `C(3k,k) C(6k,3k)^2 / ((3k+1)(3k+2) 2^(15k))` — II v10 Conj 5.25(i)

* **B46** `sum_{k>=0} C(3k,k)C(6k,3k)^2 ((360k^3+612k^2+230k+15)(H_{2k}-H_k) - 180k^2+36k+23) / ((3k+1)(3k+2) 2^(15k)) = 48 sqrt2 log2/pi`

Status: OPEN (base A7 open). (A second identity with `4H_{6k}-3H_{3k}-H_k` and an auxiliary f(k) is in II; not transcribed.)

### K16. `C(2k,k)^2 C(4k,2k)^3 / ((2k+1) 2^(22k))` — II v10 Conj 5.33(i) (2025-06-26)

P(k) = 1344k^3+944k^2+156k+9.

* **B47** `sum_{k>=0} C(2k,k)^2C(4k,2k)^3 ((3H_{4k}-2H_{2k}-H_k)P(k) + 672k^2+404k+34 + 1/(4k+2)) / ((2k+1) 2^(22k)) = 352 sqrt2 log2/pi^2`
* **B48** `sum_{k>=0} C(2k,k)^2C(4k,2k)^3 ((32H_{4k}^{(2)} - 8H_{2k}^{(2)} - H_k^{(2)})P(k) + 32(12k+5)) / ((2k+1) 2^(22k)) = 352 sqrt2/3`

Status: OPEN (base A8 open; no hit).

### K17. `C(2k,k)^5 / ((6k+1) (-64)^k C(3k,k) C(6k,3k))` — II v10 Conj 5.37(i), 5.39

Base `sum (28k^2+10k+1)·kernel = 3/pi` and the two first-order identities (II Conj 5.35) are now proved (HS2, theorem with eqs (2-1)–(2-3)).

* **B49** `sum_{k>=0} C(2k,k)^5 ((28k^2+10k+1)(10H_{2k}^{(2)} - 3H_k^{(2)}) + 2) / ((6k+1)(-64)^k C(3k,k) C(6k,3k)) = pi/2`
* **B50** `sum_{k>=0} C(2k,k)^5 ((28k^2+10k+1)(28H_{2k}^{(3)} - 3H_k^{(3)}) + 8/(2k+1)) / ((6k+1)(-64)^k C(3k,k) C(6k,3k)) = 75 zeta(3)/pi - 24G`

Status: OPEN as far as determined. HS2 (which has the WZ pair for this kernel) states only the order-1 results;
V v3 Conj 5.12 (g_12'', g_12''' sums) was still conjectural three days earlier. Prime WZ target: the pair exists.

### K18. `256^k / (k^7 C(2k,k)^7)` (Guillera's pi^4/8; base proved by Au-seed)

P(k) = 21k^3-22k^2+8k-1.

* **B51** `sum_{k>=1} 256^k P(k) (64H_{2k-1}^{(6)} - 65H_{k-1}^{(6)}) / (k^7 C(2k,k)^7) = 31 pi^10/3780` — II v10 Conj 5.41(i); MO 486518 (2025-01-24, no answer)
* **B52** `sum_{k>=1} 256^k (14P(k)(H_{2k-1}-H_{k-1}) - 63k^2+44k-8) / (k^7 C(2k,k)^7) = (93/2) zeta(5) + pi^4 log2` — V v3 Conj 3.11 (2026-04-12)

Status: OPEN as far as determined (no third-party hit with harmonic numbers on this kernel). Not determined: the two
earlier companions quoted in II Remark 5.42, `… (P(k)(4H_{2k-1}^{(2)}-5H_{k-1}^{(2)}) - 6k+2) = pi^6/24` (H eq.(178), Conj 65)
and `… (P(k)(16H_{2k-1}^{(4)}+7H_{k-1}^{(4)}) + 4/k) = 31pi^8/1440` (I eq.(5.20)) — II still calls them "conjectured".

### K19. `4096^k / (k^7 C(2k,k)^6 C(3k,k))` (Guillera's 8pi^4; base proved by Au-seed Example IV) — II v10 Conj 5.45; MO 487145 (no answer)

Q(k) = 92k^3-84k^2+27k-3.

* **B53** `sum_{k>=1} 4096^k (Q(k)(5H_{2k-1}-2H_{k-1}) - 94k^2+57k-9) / (k^7 C(2k,k)^6 C(3k,k)) = 2976 zeta(5)`
* **B54** `sum_{k>=1} 4096^k (Q(k)(H_{3k-1}-3H_{k-1}) - 88k^2/3+18k-3) / (k^7 C(2k,k)^6 C(3k,k)) = 16(2 pi^4 log2 - 93 zeta(5))`
* **B55** `sum_{k>=1} 4096^k (Q(k)(8H_{2k-1}^{(2)} - 7H_{k-1}^{(2)}) - 44k+12) / (k^7 C(2k,k)^6 C(3k,k)) = (16/3) pi^6`
* **B56** `sum_{k>=1} 4096^k (Q(k)(64H_{2k-1}^{(4)} + 13H_{k-1}^{(4)}) + 32/k) / (k^7 C(2k,k)^6 C(3k,k)) = (248/45) pi^8`

Status: OPEN as far as determined (`92,84,27` only in Au-seed for the base identity and in Cohen–Guillera).

### K20. `C(2k,k)^8 C(4k,2k) / 2^(32k)` (Cohen–Guillera 1/pi^4 kernel) — II v10 Conj 5.47 (2025-03-01)

R(k) = 43680k^4+20632k^3+4340k^2+466k+21 (base `sum R(k)·kernel = 2048/pi^4`: "conjectural" according to Sun; hits in Au-seed not read).

* **B57** `sum_{k>=0} C(2k,k)^8C(4k,2k)/2^(32k) * (R(k)(17H_{2k}^{(4)} - H_k^{(4)}) + (1796k+193)/(2k+1)) = 8704/45`
* **B58** `sum_{k>=0} C(2k,k)^8C(4k,2k)/2^(32k) * (R(k)(127H_{2k}^{(6)} - 2H_k^{(6)}) + 1323(4k+1)/(2k+1)^3) = (126976/945) pi^2`
* **B59** `sum_{k>=0} C(2k,k)^8C(4k,2k)/2^(32k) * (R(k)(4354H_{2k}^{(8)} - 17H_k^{(8)}) - 4182(4k+1)/(2k+1)^5) = -(197632/4725) pi^4`

Status: OPEN as far as determined. Caveat: these are the "upside-down" twins of II Conj 5.43 (Chen's kernel), which
N. Shalev claims to have proved in Oct 2026 by a 4-parameter Dougall-7F6 seed identity (see Section 3); the same
parametric identity summed "over the half-integers" may well give B57–B59. Check arXiv:2610.04770 before attacking.

### K21. 1/pi^2 kernels with first-order harmonic numbers ("log m / pi^2") — H Conj 59–62 (2022-12 / 2023-01)

* **B60** `sum_{k>=0} C(2k,k)^3C(3k,k)C(4k,2k)/(-24^4)^k * ((252k^2+63k+5)(4H_{4k}+3H_{3k}-7H_k) + 504k+63) = 192 log(24)/pi^2`
* **B61** `sum_{k>=0} C(2k,k)^2C(3k,k)^2C(6k,3k)/10^(6k) * (3(532k^2+126k+9)(H_{6k}-H_k) + 532k+63) = 1125 log(10)/(4 pi^2)`
* **B62** `sum_{k>=0} C(2k,k)^2C(3k,k)^2C(6k,3k)/(-2^18)^k * (6(1930k^2+549k+45)(H_{6k}-H_k) + 3860k+549) = 6912 log2/pi^2`
* **B63** `sum_{k>=0} C(2k,k)^2C(3k,k)^2C(6k,3k)/(-2^18 3^6 5^3)^k * (2(5418k^2+693k+29)(H_{6k}-H_k) + 3612k+231) = (128 sqrt5/pi^2) log(2^6 3^2 5)`

Status: OPEN as far as determined: no third-party hit beyond the base identities (base of B62 proved by Au-seed; base of
B63 itself still conjectural, see end of Section 1, so B63 is certainly open). V v3 Conj 5.17–5.20 (Apr 2026) restate the
first-derivative sums `sum g'(k)` for these kernels as conjectures (the real/imaginary parts encode B60–B63).
H Conj 63, 64 (kernels (-2^22 3^3)^k and (2^18 7^4)^k) are of the same kind (not transcribed: they need Sun's auxiliary H(k)).

### K22. `z^k / ((2k+1)(6k+1)(6k+5) C(6k,3k))`, z = 8, 16, 1, 27 — III v6 Conj 3.1–3.4 (2025)

Base: `sum (63k^2+78k+22)8^k/(…) = 3pi/2` (III Thm 1.1), A10 (open), `sum (21k^2+27k+8)/(…) = 8pi/(9 sqrt3)` (Chu–Zhang Ex. 115),
`sum 27^k(111k^2+133k+36)/(…) = 16pi/(3 sqrt3)` (Chu–Zhang Ex. 58). Write D(k) = (2k+1)(6k+1)(6k+5) C(6k,3k).

* **B64** `sum_{k>=0} 8^k ((63k^2+78k+22) H_{3k} + 12k+8) / D(k) = 6G - (3/2) pi log2`
* **B65** `sum_{k>=0} (3k+2)^2 16^k H_k / D(k) = (5/2)K - pi log(432)/(6 sqrt3)`
* **B66** `sum_{k>=0} 16^k ((3k+2)^2 (3H_{3k}-2H_{2k}) + 2k(3k+2)/(2k+1)) / D(k) = (5/4)K + (pi/(6 sqrt3)) log(27/256)`
* **B67** `sum_{k>=0} ((21k^2+27k+8)(H_{2k}+2H_k) + (3k+1)(5k+4)/(6k+3)) / D(k) = 14K - 16 log3 pi/(3 sqrt3)`
* **B68** `sum_{k>=0} ((21k^2+27k+8) H_{3k} + 2k + 13/9) / D(k) = (8/3)K - 8 log3 pi/(9 sqrt3)`
* **B69** `sum_{k>=0} 27^k ((111k^2+133k+36)(H_{2k}-3H_k) - (3k+1)(3k+4)/(2k+1)) / D(k) = 16 log3 pi/sqrt3 - 48K`
* **B70** `sum_{k>=0} 27^k ((111k^2+133k+36) H_{3k} + 30k+19) / D(k) = 12K`

Status: OPEN as far as determined (`63,78,22`, `111,133,36`: no hit). HS2 (Thm for f_3) proves the *derivative* sums for the
21k^2+27k+8 kernel (`sum f_3'(k) = -12K`, `sum f_3'' = 32pi^3/(9 sqrt3)`, `sum f_3''' = 54(2pi^2 K - 27 L_{-3}(4))`), whose
first-order content is the single combination `H_{3k}-H_{6k}`; B67, B68 are different combinations.
Single-binomial kernel: also accessible by beta-integral methods (as in JS for C(4k,k)).

### K23. `(-1)^k / (k^2 C(2k,k) C(3k,k))` — Sun (I) eq.(4.1)

* **B71** `sum_{k>=1} (-1)^k/(k^2 C(2k,k) C(3k,k)) * ((7k-2)/(2k-1) H_{k-1} - 1/(6k)) = (pi^2 log2)/6 - (11/12) zeta(3)`

Status: OPEN — the best-documented open case. Au-H v4 (2026-11-09), end of Example V: the seven neighbouring conjectures
(I) (4.9)–(4.15) follow from his formulas, but (4.1) "evidently does not follow from the formulas above because of the
presence of log 2. We encourage interested readers to seek other WZ-pairs that might prove it." (Statement transcribed
from Au-H with a_n = (-1/27)^n (1)_n^2/((1/3)_n(2/3)_n) = (-1)^n/(C(2n,n)C(3n,n)) and the factor 1/n^2 restored;
numerically confirmed to 96 digits in exactly this form. HS1 proves the companion
`sum (-1)^k/(k^2C(2k,k)C(3k,k)) ((7k-2)/(2k-1)(2H_{2k-1}-H_{k-1}) + 1/k) = -(3/2) zeta(3)`.)

### K24. `1/(k^2 C(2k,k))`, quadratic in harmonic numbers — V v3 Conj 3.10 (2026-03-09)

* **B72** `sum_{k>=1} H_{k-1}^{(2)}/(k^2 C(2k,k)) * (H_{2k} + 26H_k - 8/k) = (34 pi^2 zeta(3) - 287 zeta(5))/54`
  (V v3 prints "-" instead of "="; the value with this sign was confirmed by PSLQ and to 92 digits.)

Status: OPEN as far as determined, but probably **not hard**: inverse central binomial sums of this type are reducible to
level-3/level-6 cyclotomic MZVs (Zhou Z1; Au's MZV package; SZ2). Low priority for a WZ project.

### K25. `z^k / (k (3k-1)(3k-2) C(4k,k))`, z = 9/8, 9, -1/2, -1/8 — II v10 Conj 4.15, 4.19, 4.21, 4.23 (Feb 2025)

Base identities: II Thm 1.1 (z=9/8), A6 (z=9, open), II Thm 1.2 (z=-1/2, -1/8). D(k) = k(3k-1)(3k-2) C(4k,k).

* **B73** `sum_{k>=1} (9/8)^(k-1) ((95k^2-84k+16)(H_{4k-1}-H_{k-1}) - 25k+12) / D(k) = (log3/sqrt3) pi + (15/4)K`
* **B74** `sum_{k>=1} 9^(k-1) ((13k^2-15k+2)(3H_{2k-1}-4H_{k-1}) + 12(4k-1)) / D(k) = (27/4)K`
* **B75** `sum_{k>=1} 9^(k-1) ((13k^2-15k+2)(8H_{4k-1}-9H_{2k-1}) - 160k+84) / D(k) = (4 log3/sqrt3) pi - (93/4)K`
* **B76** `sum_{k>=1} ((77k^2-53k+10)(H_{2k-1}-H_{k-1}) - (209k^2-113k+18)/(2k)) / (D(k) (-2)^k) = 3 log^2 2`
* **B77** `sum_{k>=1} ((77k^2-53k+10)(2H_{4k-1}-H_{2k-1}-H_{k-1}) + 11k-27+6/k) / (D(k) (-2)^k) = -pi^2/2`
* **B78** `sum_{k>=1} ((415k^2-343k+62)(H_{4k-1}-H_{k-1}) - (581k^2-229k-6)/(4k)) / (D(k) (-8)^k) = -pi^2/4`

Status: OPEN as far as determined (JS proves the C(4k,k)-in-the-numerator conjectures II 4.1–4.13, not these).
Single-binomial kernel; integration methods (as in II §2–3, JS) are the natural tool rather than WZ.

### K26. `C(4k,k) / ((3k+1)(3k+2) 16^k)`, second-order harmonic numbers — II v10 Conj 4.3(i), 4th identity

* **B79** `sum_{k>=0} C(4k,k) ((11k^2+8k+1)(H_{2k}^{(2)} + (3/4)H_k^{(2)}) + (4k+1)/(2k+1)) / ((3k+1)(3k+2) 16^k) = pi^2/6`

Status: OPEN as far as determined: JS (its second theorem) proves the three first-order identities of Conj 4.3 but not this one
(no `pi^2/6` in JS).

### Not individually listed (but open as far as determined): Sun's "derivative" conjectures, V v3 Sections 4–5

V v3 Conj 4.2, 4.4–4.6 and 5.1–5.27 assert values of `sum_k f^{(m)}(k)` for 30 hypergeometric terms f extended by Gamma
functions (m up to 9); these are hypergeometric sums with *polynomials* in harmonic numbers and are exactly what
"differentiate the WZ pair in n" produces (HS2, second observation). Status seen: Conj 4.1 (f_1) proved (Section 3);
Conj 4.3 (f_3 = 1/(x^4 C(2x,x))) proved by SZ2 (its abstract gives A_{4,4} = 352 zeta(5,3)/15 + 752537 pi^8/10206000, which is
Sun's value); HS2 proves the analogous statements for 15 other kernels (its f_1 … f_15) — *not* the g_1 … g_27 of V Section 5
(MO 508786: H. Cohen notes that such observations were made without proof in Cohen–Guillera §6). I did not verify these
numerically; they are a large reservoir of further WZ targets (e.g. V Conj 5.12: `sum g_12'' = -12 pi`, `sum g_12''' = 288G - 30 pi^2 i`
for the kernel of K17, equivalent in content to B49/B50-type statements).

---------------------------------------------------------------------------------------------------

## 3. PROVED recently — do not re-prove

(Only what I actually saw; "per X" = I rely on a statement in X and did not read the proof.)

**From H (arXiv:2210.07238v9; conjecture numbers of v9).**
* Conj 1(i) `sum (-1)^(k-1)/(k^3C(2k,k)) (H_{2k-1}^{(2)} - 123/16 H_{k-1}^{(2)}) = 451 zeta(5)/40 - 14 pi^2 zeta(3)/15` — Au-MZV Cor. 2.9 (per H Remark).
* Conj 2–10, 12–16 (Section 2: one or two binomial coefficients; C(2k,k)/8^k, /16^k with H^{(3)}; 2^k C(3k,k) families; C(4k,2k)C(2k,k)/m^k
  with 2H_{4k}-H_{2k}; C(6k,3k)C(3k,k)/864^k …) — Y. Zhou, Z1 (cites H Conj 2.2–2.8) and Z2 (cites H Conj 9, 10, 12–16, eqs (29),(33)–(38)).
* Conj 11 `sum 48^k/(k(2k-1)C(2k,k)C(4k,2k)) (6H_{4k-1}-9H_{2k-1}+2H_{k-1}+6/(2k-1)) = 2pi^3/sqrt3` — Au-seed (formula found verbatim there, derived from a three-variable identity).
* Conj 17(i) `… (21k-8)(H_{2k-1}^{(2)} - 25/8 H_{k-1}^{(2)}) … = 47pi^4/2880`, Conj 20(i) `… (3k-1)16^k (H_{2k-1}^{(2)} - 5/4 H_{k-1}^{(2)}) … = pi^4/24`,
  Conj 22(i) `… (4k-1)(-64)^k … = -16 beta(4)` — Au-MZV Cor. 2.3, 2.11 (per H Remarks).
* Conj 18(i) `sum (21k-8)/(k^3C(2k,k)^3) (H_{2k-1}^{(3)} + 43/8 H_{k-1}^{(3)}) = 711 zeta(5)/28 - 29 pi^2 zeta(3)/14` — Wei–Xu arXiv:2310.04642 (formula present in the source).
* Conj 19(i) `sum (3k-1)(-8)^k/(k^3C(2k,k)^3) (H_{2k-1}^{(2)} - 5/4 H_{k-1}^{(2)}) = -2 beta(4)` — C. Wei arXiv:2303.05402 (by its title "On a conjectural series of Sun for beta(4)"; source only grepped).
* Conj 21(i) `sum (3k-1)16^k/(k^3C(2k,k)^3) (H_{2k-1}^{(3)} + 7/8 H_{k-1}^{(3)}) = pi^2 zeta(3)/2` — Li–Chu, Electron. Res. Arch. 31 (2023) (per Au-H, remark closing Example VI: the non-"?=" formulas there "were conjectured by Sun and proved in" Li–Chu).
* Conj 24(i) `sum (-27)^k((15k-4)(3H_{3k-1}-H_{k-1}) - 9)/(k^3C(2k,k)^2C(3k,k)) = -4pi^3/sqrt3` — HS1, first theorem of the introduction.
* Conj 29 (General Conjecture: `(log|m|)/pi` companions of all 36 rational Ramanujan series, parts (i)–(iv)) — Y. Zhou Z2 §3 via modular equations (per Au-H); see also Hou–He–Wang 2024. Includes e.g. `sum C(2k,k)^3/4096^k ((42k+5)(H_{2k}-H_k)+7) = 32 log2/pi`.
* Conj 30 `sum (4k+1)C(2k,k)^3/(-64)^k H_{2k}^{(3)} = 15zeta(3)/(4pi) - 2G` — Wei–Xu (per SZ1); Conj 33 `(6k+1)/256^k, H_{2k}^{(3)} - 7/64 H_k^{(3)}` — Li–Chu, Axioms 2025 (per SZ1);
  Conj 35 `(6k+1)/(-512)^k` — Li–Chu 2024 (per SZ1); Conj 37 `sum (42k+5)C(2k,k)^3/4096^k (H_{2k}^{(2)} - 25/92 H_k^{(2)}) = 2pi/69` — C. Wei arXiv:2211.11484;
  **Conj 38** `sum (42k+5)C(2k,k)^3/4096^k (H_{2k}^{(3)} - 43/352 H_k^{(3)}) = 555 zeta(3)/(77 pi) - 32G/11` — **Sun–Zhou SZ1 (Feb 2026), modular (Eichler integrals)**;
  Au-H stresses that it is "inaccessible to our WZ approach despite the availability of three WZ-pairs" — a *WZ* proof would be new.
* Conj 51 (eqs (141),(142): C(2k,k)^5/(-4096)^k with H_k, H_{2k}) and Conj 53 (eq (148): C(2k,k)^5/(-2^20)^k with H_{2k}-H_k) — C. Wei (per V v3 Remarks 5.13, 5.14; the second-order companions Conj 52, 54 were not checked); H^{(3)} companions — Wei–Xu arXiv:2308.06440 (per HS1); an H^{(4)} companion — HS1.
* Conj 66–68 (eqs (181),(183),(185): Gourevich kernel C(2k,k)^7/2^(20k) with H, H^{(2)}, H^{(4)}) — Au-H Example on the 1/pi^3 kernel.
* Status **not determined** (neither counted as open nor as proved): H Conj 31 (256^k with H_k, H_{2k} and Gamma(5/6)), 39 (216^k), 40, 41 ((-1024)^k),
  42 (48^2; Au-H Ex. VIII still lists `sum a_n(8n+1)(1 + (3H_{2n}-4H_n)/6) ?= 8 log2/(sqrt3 pi)` with "?="), 45–50, 55–58, 63–65, 69, 70.

**From I / II (arXiv:2307.03086).**
* (I) (4.9)–(4.15) (seven identities on `(-1)^k/(k^2C(2k,k)C(3k,k))`) — Au-H Example V; (I) (4.12) `sum (-1)^k(28k^2-18k+3)/((2k-1)^3k^4C(2k,k)C(3k,k)) = -pi^4/45` — Au-H; (I) (5.21), (5.30) — Au (per V Remarks 5.4, 5.24).
* `sum (145k^2-104k+18)/((2k-1)k^3C(2k,k)C(3k,k)^2) = pi^2/3` and four harmonic companions — Au (WZ), MO 456443 answer (2023-10-16), per II Remark 5.22.
* `sum (35k^2-29k+6)3^k/(k(3k-1)(3k-2)C(4k,k)) = sqrt3 pi` — Au-seed (per II Remark 1.2); `sum (10k-1)C(6k,3k)C(3k,k)/((2k+1)512^k) = 0` — Au (per V Remark 5.1).
* II v10 Conj 4.1–4.13 (C(4k,k)/16^k, 24^k, (-25)^k, (-72)^k, 128^k, (-256)^k with H(k)=2H_{4k}-3H_{2k}+H_k etc.) — Jiang–Sun JS (integration), except the H^{(2)} identity B79.
* II v10 Conj 5.11 `sum C(2k,k)^3/4096^k ((42k+5)H_k^{(3)} - 352/(2k+1)^2) = (32/7)(335 zeta(3)/pi - 224G)` — appears as eq. (SunH3e) among the results of SZ1 (I saw the statement in SZ1's introduction; proof not read).
* II v10 Conj 5.27 `sum (-1)^k(28k^2-8k+1)C(2k,k)^2/(k(2k-1)^2C(6k,3k)C(3k,k)) = -2 log2` — HS2, first theorem (eq. `id-1`), by the WZ method.
* II v10 Conj 5.29 (60k^2-26k+3 kernel, four identities) and Conj 5.31 (364k^2-227k+36 kernel, three identities) — HS1 (stated there as confirming "Conjecture 5.29 / 5.31").
* II v10 Conj 5.35 (`sum (28k^2+10k+1)C(2k,k)^5/((6k+1)(-64)^kC(3k,k)C(6k,3k)) = 3/pi` and the two log2/pi companions) — HS2, theorem with eqs (2-1)–(2-3), by the WZ method.
* II v10 Conj 5.43 (D. Chen's kernel `C(4k,2k)/((4k-1)k^7C(2k,k)^8)`, P(k)=5460k^4-…: values pi^6/189, pi^8/1350, 2pi^10/31185, and Cohen's pi^12 one) —
  **claimed** by Noam Shalev, MO 487033 answer of 2026-10-27 (4-parameter identity from Au's Dougall7F6 seed + Guillera's periodicity argument; the author
  states the results were obtained with AI tools), cf. arXiv:2610.04770 "Borwein–Bradley-type identities for Ramanujan-type series for 1/pi^4"
  (abstract not extractable from the source; body not read). Treat as proved-pending-refereeing.

**From V / MO 2026.**
* MO 486353 `sum (-1)^k(560k^4-640k^3+408k^2-136k+17)/(k^5(2k-1)^4C(2k,k)C(3k,k)) = 180 zeta(5) - (56/3)pi^2 zeta(3)` — Au (MO answer; Au-H Example III).
* V v3 Conj 4.1 (`sum f_1^{(m)}(k)`, f_1 = Gamma(x)^2/(2 Gamma(2x)), m=1,2,3: -3K/2, 4pi^3/(27 sqrt3), (3/4)(2pi^2K - 27L_{-3}(4))) — Qiping Zhou, MO 508769 answer (2026-07-23, beta-integral proof); also within SZ2's A_{1,n}.
* V v3 Conj 4.3 (1/(x^4C(2x,x)), derivatives 2–5) = MO 508768 Conj II — SZ2 (per abstract). MO 508768 Conj I, MO 508774 (2)–(7), MO 508844 Conj 1–2, MO 508845, MO 508850 (all but its formula (9)) — HS2 (per Sun's own MO answers of 2026-04-17).
* MO 508663 ("Few conjectural series for zeta(5) and zeta(6)", 2026-03-02): identity (1) — M. Wildon (from Li–Chu 2023 formulas, found with an LLM and checked by him); (2) = Chu 2020 Example 3.6, (3), (4) — N. Kalinin (MO answers 2026-03-03).
* MO 507709 (four vanishing series `sum C(2k,k)^3/m^k [H_{2k}^{(2)} - c H_k^{(2)} + const] = 0`, m = -64, 256, -512, 4096) — SZ1 (first one also E. Naslund on MO).
* MO 488044 `sum (22k^2-92k+11)C(4k,k)/16^k = -5` — M. Alekseyev (MO); MO 488056 (2)–(5) — JS.

---------------------------------------------------------------------------------------------------

## 4. Non-hypergeometric conjectures (short; secondary interest)

Not WZ-reachable in the single-sum sense; listed only for orientation. I did **not** check their status beyond the
fact that they appear as conjectures in the versions cited, and did not verify them numerically.
T_n(b,c) = coefficient of x^n in (x^2+bx+c)^n; f_k = Franel numbers sum_j C(k,j)^3; g_k = sum_j C(k,j)^2 C(2j,j).

1. `sum_{k>=0} (24k+5)/76^(2k) C(2k,k) T_k(1,-12) T_{2k}(8,-3) = sqrt(38(120+73 sqrt3))/(6 pi)` — V v3 Conj 2.1 (X1), 2026-03-18.
2. `sum_{k>=0} (585k+136)/(-85^2)^k C(2k,k) T_k(8,-9) T_{2k}(7,-9) = 85 sqrt255/(6 pi)` — V v3 Conj 2.1 (X2).
3. `sum_{k>=0} (16k+3)/(-202^2)^k C(2k,k) T_k(19,-20) T_{2k}(9,-5) = 43 sqrt101/(75 pi)` — V v3 Conj 2.1 (X3), 2026-03-20 (the abstract's showcase).
4. `sum_{k>=0} (3k+1)/(-48)^k f_k T_k(4,-2) = 4 sqrt2/(3 pi)` — arXiv:1911.05456v10 (2020), family "f_k T_k(b,c)".
5. `sum_{k>=0} (99k+23)/(-288)^k f_k T_k(8,-2) = 39 sqrt2/pi` — ibid.
6. `sum_{k>=0} (8k+3)/(-81)^k g_k T_k(7,-8) = 9 sqrt3/(4 pi)` — ibid., family "g_k T_k(b,c)".
7. `sum_{k>=0} (4k+1)/(-1089)^k g_k T_k(31,-32) = 33/(16 pi)` — ibid.
8. `sum_{k>=0} (39600310408k+39624469807)/39216^(2k) C(2k,k+1) C(4k,2k) T_k(39202,1) = 1334161 sqrt817/pi` — ibid., (III12').

(Sun's best-known item of this kind, the "$520 series" `sum (1054k+233)/480^k C(2k,k) sum_j C(k,j)^2C(2j,j)(-1)^j 8^(2j-k) = 520/pi`
from arXiv:1102.5649 eq.(3.24), is — to my recollection, not re-checked here — long settled (Rogers–Straub), so it is not listed as open.)

---------------------------------------------------------------------------------------------------

## 5. Suggested priorities for a WZ project (my reading of the evidence)

1. **K17 (B49, B50)**: a WZ pair for the kernel was published in April 2026 (HS2) and only order-1 consequences were drawn.
2. **K4 (B12)** and the other "?=" identities of Au-H Example VII: the WZ pairs are known; what is missing is the rigorous
   evaluation of `sum g(n)` and `lim_n sum_k F(n,k)` (Au says so explicitly). Same for Au-H Examples VIII, IX (48^2 and (-32)^3 kernels).
3. **K1–K3 (B1–B11)**: ten(+one) fresh 2026 conjectures on six classical Ramanujan kernels, no proof attempts on record;
   mirror Au-H Examples VIII/IX (find a seed proving the base series, extract first/second-order coefficients).
4. **K6, K7 (B14–B23)**: flagged "remains open" by Sun in Feb 2026; base identities have WZ proofs (Guillera, Guillera–Rogers,
   Hessami Pilehrood), and HS1 already did one member (H Conj 24(i)).
5. **B71** ((I) (4.1)): Au explicitly asks for "other WZ-pairs"; **K10 (B28–B32)** has the same log 2 feature on a kernel HS1's method handles.
6. **A1, A2, A3–A5, A9**: parameter-free hypergeometric identities with rates 4/27, -1/27, 1/81, 1/9, 1/9, -27/512 — classical WZ-seed search problems.
7. K13, K16, K18, K19, K20, K21: higher-degree kernels; B57–B59 should be compared with Shalev's parametric identity first.

---------------------------------------------------------------------------------------------------

## 6. Numerical verification table (from `verify_B.out.txt`)

`terms` = number of terms summed; `digits` = floor(-log10(|lhs-rhs| / max(1,|rhs|))) at working precision 95.

| id | terms | abs(lhs-rhs) | digits | source |
|---|---|---|---|---|
| A1 | 110 | 3.6305e-94 | 93 | MO 506695 (2026-01-06) |
| A2 | 64 | 8.9888e-95 | 95 | MO 487361 eq.(3) (2025-02-07) |
| A3 | 51 | 9.3634e-97 | 96 | II v10 Conj 5.1(i) eq.1 |
| A4 | 97 | 2.0599e-95 | 94 | II v10 Conj 5.1(i) eq.2 |
| A5 | 97 | 8.5207e-95 | 94 | II v10 Conj 5.9(i) eq.1; MO 456779 |
| A6 | 3794 | 1.4153e-87 | 87 | II v10 Conj 4.17(i) |
| A7 | 1183 | 1.9795e-88 | 88 | II v10 Conj 5.25(i) eq.1 (=(5.20)) |
| A8 | 55 | 0.0 | 95 | II v10 Conj 5.33(i) eq.1 (=(5.31)) |
| A9 | 75 | 7.4907e-96 | 95 | II v10 Conj 4.33(i) eq.1 |
| A10 | 150 | 2.0469e-92 | 91 | III v6 Conj 3.2 |
| B1 | 353 | 1.141e-90 | 89 | 2603.29973v3 Conj 3.1; MO 507592 |
| B2 | 353 | 1.8071e-90 | 89 | 2603.29973v3 Conj 3.1 |
| B3 | 80 | 0.0 | 95 | 2603.29973v3 Conj 3.2 |
| B4 | 81 | 2.397e-94 | 95 | 2603.29973v3 Conj 3.2 (sign of rational term corrected) |
| B5 | 53 | 1.4981e-95 | 95 | 2603.29973v3 Conj 3.3 (sign of rational term corrected) |
| B6 | 53 | 1.8727e-96 | 95 | 2603.29973v3 Conj 3.3 |
| B7 | 59 | 2.9963e-95 | 95 | 2603.29973v3 Conj 3.4 |
| B8 | 59 | 3.7453e-96 | 95 | 2603.29973v3 Conj 3.4 |
| B9 | 60 | 0.0 | 95 | 2603.29973v3 Conj 3.5; MO 507600 (3) |
| B10 | 311 | 8.5747e-91 | 91 | 2603.29973v3 Conj 3.6 |
| B11 | 40 | 9.5881e-94 | 95 | 2603.29973v3 Conj 3.7 |
| B12 | 56 | 0.0 | 95 | II v10 Conj 5.13(i) (=I (5.14)); Au 2602.08721v4 '?=' |
| B13 | 55 | 1.4981e-95 | 95 | 2603.29973v3 Conj 3.8; MO 507600 (1) |
| B14 | 84 | 2.2472e-95 | 95 | 2210.07238v9 Conj 23(i) eq.1 (=(63)) |
| B15 | 84 | 1.4981e-95 | 95 | 2210.07238v9 Conj 23(i) eq.2 |
| B16 | 85 | 0.0 | 95 | II v10 Conj 5.15(i) |
| B17 | 393 | 8.3987e-90 | 89 | 2210.07238v9 Conj 26(i) eq.1 (=(71)) |
| B18 | 393 | 6.5305e-90 | 89 | 2210.07238v9 Conj 26(i) eq.2 |
| B19 | 399 | 9.0142e-90 | 91 | II v10 Conj 5.17(i) |
| B20 | 153 | 8.3296e-93 | 93 | 2210.07238v9 Conj 25(i) |
| B21 | 152 | 1.173e-92 | 92 | 2603.29973v3 Conj 3.9; MO 507600 (2) |
| B22 | 185 | 4.6982e-92 | 93 | 2210.07238v9 Conj 27(i) (=(75)) |
| B23 | 187 | 5.4173e-92 | 93 | II v10 Conj 5.19(i) |
| B24 | 46 | 3.4083e-94 | 93 | II v10 Conj 5.21(i) |
| B25 | 36 | 7.9486e-93 | 92 | II v10 Conj 5.23(i) |
| B26 | 52 | 3.7453e-96 | 95 | II v10 Conj 5.1(i) eq.3 |
| B27 | 97 | 6.7416e-95 | 94 | II v10 Conj 5.1(i) eq.4 |
| B28 | 151 | 1.0077e-92 | 91 | II v10 Conj 5.3(i) first |
| B29 | 151 | 1.1281e-92 | 92 | II v10 Conj 5.3(i) second |
| B30 | 151 | 1.2622e-92 | 92 | II v10 Conj 5.3(i) third |
| B31 | 150 | 1.1835e-92 | 93 | II v10 Conj 5.5(i) first |
| B32 | 150 | 8.5693e-93 | 93 | II v10 Conj 5.5(i) second |
| B33 | 393 | 1.8471e-90 | 91 | II v10 Conj 5.7 first |
| B34 | 394 | 1.9797e-90 | 92 | II v10 Conj 5.7 second |
| B35 | 98 | 8.0525e-95 | 94 | II v10 Conj 5.9(i) eq.2 |
| B36 | 98 | 4.9158e-95 | 94 | II v10 Conj 5.9(i) eq.3 |
| B37 | 391 | 9.1139e-90 | 89 | II v10 Conj 4.25 first |
| B38 | 391 | 1.0335e-89 | 89 | II v10 Conj 4.25 second |
| B39 | 113 | 1.0337e-93 | 94 | II v10 Conj 4.27 first |
| B40 | 114 | 2.9963e-94 | 94 | II v10 Conj 4.27 second |
| B41 | 69 | 5.9925e-95 | 95 | II v10 Conj 4.29 |
| B42 | 244 | 5.1812e-91 | 92 | II v10 Conj 4.31(i) first |
| B43 | 245 | 9.5677e-91 | 91 | II v10 Conj 4.31(i) second |
| B44 | 243 | 5.0218e-91 | 92 | II v10 Conj 4.31(i) third |
| B45 | 74 | 9.3634e-97 | 96 | II v10 Conj 4.33(i) eq.2 |
| B46 | 1180 | 2.285e-88 | 88 | II v10 Conj 5.25(i) eq.2 |
| B47 | 56 | 1.1985e-94 | 95 | II v10 Conj 5.33(i) eq.2 |
| B48 | 56 | 1.1985e-93 | 95 | II v10 Conj 5.33(i) eq.3 |
| B49 | 68 | 1.8727e-96 | 95 | II v10 Conj 5.37(i) |
| B50 | 69 | 0.0 | 95 | II v10 Conj 5.39 |
| B51 | 55 | 9.5881e-94 | 95 | II v10 Conj 5.41(i); MO 486518 |
| B52 | 56 | 3.5955e-94 | 95 | 2603.29973v3 Conj 3.11 |
| B53 | 116 | 7.6705e-93 | 95 | II v10 Conj 5.45 first; MO 487145 |
| B54 | 116 | 1.9176e-93 | 95 | II v10 Conj 5.45 second |
| B55 | 115 | 7.6705e-93 | 95 | II v10 Conj 5.45 third |
| B56 | 117 | 0.0 | 95 | II v10 Conj 5.45 fourth |
| B57 | 32 | 4.794e-94 | 95 | II v10 Conj 5.47 first |
| B58 | 32 | 3.8352e-93 | 95 | II v10 Conj 5.47 second |
| B59 | 32 | 0.0 | 95 | II v10 Conj 5.47 third |
| B60 | 59 | 5.9925e-95 | 96 | 2210.07238v9 Conj 59(i) |
| B61 | 73 | 3.5955e-94 | 95 | 2210.07238v9 Conj 60(i) |
| B62 | 126 | 9.5881e-94 | 95 | 2210.07238v9 Conj 61(i) |
| B63 | 22 | 4.794e-94 | 95 | 2210.07238v9 Conj 62(i) |
| B64 | 104 | 2.3596e-94 | 93 | III v6 Conj 3.1 |
| B65 | 151 | 2.856e-92 | 91 | III v6 Conj 3.2 second |
| B66 | 152 | 9.5502e-93 | 92 | III v6 Conj 3.2 third |
| B67 | 55 | 7.9588e-96 | 95 | III v6 Conj 3.3 first |
| B68 | 55 | 4.6817e-97 | 96 | III v6 Conj 3.3 second |
| B69 | 243 | 4.9889e-91 | 91 | III v6 Conj 3.4 first |
| B70 | 242 | 7.4146e-91 | 91 | III v6 Conj 3.4 second |
| B71 | 67 | 7.6077e-97 | 96 | Sun (I) (4.1) as quoted in Au 2602.08721v4 |
| B72 | 151 | 9.6311e-93 | 92 | 2603.29973v3 Conj 3.10 (TeX has '-' misprinted for '='; value fixed by PSLQ) |
| B73 | 101 | 1.2734e-94 | 94 | II v10 Conj 4.15(ii) |
| B74 | 3831 | 1.3828e-87 | 87 | II v10 Conj 4.19 eq.1 |
| B75 | 3821 | 1.3748e-87 | 87 | II v10 Conj 4.19 eq.2 |
| B76 | 75 | 1.1236e-95 | 95 | II v10 Conj 4.21 eq.1 |
| B77 | 75 | 7.4907e-96 | 95 | II v10 Conj 4.21 eq.2 |
| B78 | 53 | 1.1236e-95 | 95 | II v10 Conj 4.23 |
| B79 | 389 | 8.0388e-90 | 89 | II v10 Conj 4.3(i) eq.4 (the only part of Conj 4.3 not in Jiang-Sun 2601.18645) |

Total: 89 identities, all agreeing to >= 87 digits.
