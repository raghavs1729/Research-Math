# AUDIT-INDEPENDENT — adversarial referee report and novelty check on `ramanujan-series/`

Date: 2026-11-24. Nothing under `RESULTS.md`, `OPEN-LIST.md`, `PROGRESS.md`, `proofs/`, `certificates/`, `lean/`,
`code/`, `survey/` was modified. My own scripts and their outputs are in `audit2/` (outputs are named `*.out.txt`,
because the repository `.gitignore` ignores `*.out`). Nothing was committed.

## 0. Verdict in one paragraph

All eleven identities are **true** (each re-summed independently to >= 228 digits) and all nine certificate files are
correct WZ certificates (re-derived by two independent paths). The proofs of T1, T2, T3, T4, T5, T7, T11 are complete
as written, up to small textual points; the proofs of T6, T8, T9, T10 are correct in substance but one step
("the bounds are uniform in the parameter, differentiate termwise") is asserted, not written out, and is justified by
a principle that is false as stated. The Lean claim for T1 is accurate. **The novelty claim is wrong for one identity
and needs correcting for a second:**

* **T5 is not new.** It is, term by term, equation (7.10) of W. Chu, *Infinite series formulae related to Gauss and
  Bailey 2F1(1/2)-sums*, Acta Math. Sci. 40B (2020) 293–315 — a theorem proved five years before Sun conjectured it.
* **T2 was not first conjectured by Sun.** It is J. Guillera's 2018 conjecture `G_2(1/2, e^{-pi sqrt6}) = pi/6`
  (the Ramanujan series `sum (8n+1)(1/2)_n(1/4)_n(3/4)_n/(n!^3 9^n) = 2 sqrt3/pi` "with the shift 1/2"). I found no
  proof of it, so the proof here may still be new, but the attribution and the "why nobody tried" remark are wrong.

For the other nine (in particular T4 and T11) I found no prior proof; the residual risk is named in A.4.

Findings by severity:

| # | severity | item |
|---|---|---|
| F1 | **FATAL for the novelty claim of T5** (the proof itself is fine) | T5 = Chu 2020, eq. (7.10) |
| R1 | REPAIRABLE | T2: statement is Guillera's (2018), also Cohen–Guillera (2021); attribution and context must be fixed |
| R2 | REPAIRABLE | the novelty check searched for the "defining polynomial" as a string; it cannot see index shifts or Pochhammer forms (this is how F1 and R1 were missed) |
| R3 | REPAIRABLE | T3: wrong source citations (MO 456444 and Zuniga do not contain this identity) |
| R4 | REPAIRABLE | T6, T8, T9, T10: uniformity in the parameter and termwise differentiation are asserted, with an invalid one-line justification; "every analytic step is written out" (Summary, item 2) is overstated for these four |
| M1–M6 | MINOR | see Part D |

---

## Part A — Novelty

### A.1 What was consulted (all on 2026-11-24)

* Sun's lists, TeX sources of the latest arXiv versions: 2307.03086**v10** (6 Feb 2026), 2506.01870**v6** (10 Feb 2026),
  2603.29973**v3** (13 Apr 2026); version histories from the abstract pages. Sun's arXiv listing (all 2024–2026
  entries), his later papers 2604.15172 (Hou–Sun, *Evaluations of some series via the WZ method*), 2607.14646
  (Sun–Zhou), 2602.12091 (Sun–Zhou), 2601.18645 (Jiang–Sun), 2609.04176 (Sun, 3 Sep 2026); his home page
  (`maths.nju.edu.cn/~zwsun`, updated 4 Aug 2026).
* MathOverflow through the StackExchange API: questions, **all answers, comments and revision histories** of
  456444, 456779, 506695; the bodies of all 77 questions Sun asked since 2023; content search on MathOverflow and
  Math.SE for each of the eleven series; the follow-up thread 506724.
* K. C. Au: arXiv:2312.14051v4 (every displayed series), 2212.02986 (his first WZ paper; the number 2209.01731 in the
  brief is an unrelated statistics paper), 2403.04555v4 (q-analogues, 21 examples), 2602.08721v4 (9 Nov 2026; every
  kernel `a_n`), his public formula repository.
* W. Chu: Acta Math. Sci. 40B (2020) 293–315 (full text), Math. Comp. 80 (2011) 2223–2251 (full text, 100 examples),
  arXiv:2103.07872, 2104.01599, ERA 31 (2023), Li–Chu Mathematics 12 (2024) 1381 (abstract and text).
* Guillera arXiv:1802.02023v4, 2503.00570; Cohen–Guillera arXiv:2101.12592; Hou–Mu 2604.22377; Le Lan 2609.22304,
  2604.03327; Patel 2607.27229, 2607.15303; Zhou 2503.05599, 2608.16735; Zuniga 2601.12495v2; Shalev 2610.04770,
  2610.34482; Zudilin 2506.20289; Campbell 2405.02776, 2502.15249, 2509.05897; Xu 2512.18121; Wei–Xu 2308.06440,
  2310.04642; Fernández Refolio 2607.24867, 2607.11924; Shalyt et al. 2607.09721; and others (about 60 sources).
* Citation graphs: Semantic Scholar and OpenAlex "cited by" lists for 2307.03086, 2506.01870, 2603.29973, 2312.14051,
  2604.15172 and for Sun's "Series with summands involving harmonic numbers"; arXiv abstract sweeps Apr–Nov 2026
  ("Wilf-Zeilberger", "WZ"+series, Sun+series+conjecture, Ramanujan-type, hypergeometric+conjecture) and author
  listings (Au, Sun, Hou, Zhou, Wei, Chu, Zudilin, Guillera, Campbell, Mao).
* Full-text search of all downloaded sources for each defining polynomial **and for its index shifts `P(k+j)`,
  `|j| <= 3`**, in both coefficient orders, and for the binomial/Pochhammer patterns.

Could **not** be read (publisher blocks): Chu–Zhang, Math. Comp. 83 (2014) 475–512; Chu, Ramanujan J. 55 (2021)
239–270 (very well-poised Omega-sum) and Ramanujan J. 26 (2011); Chu, Ramanujan J. (2021) Gould–Hsu; Chu, Ann.
Polon. Math. (2022). Google Scholar itself is not reachable from here. The arXiv API was rate-limited for part of
the session (author sweeps were redone through the HTML search).

### A.2 Statements in the sources (all confirmed verbatim)

| id | source, exact location | status in the source |
|---|---|---|
| T1, T6 | arXiv:2506.01870v6, Section 3, second conjecture (= Conj. 3.2), eq. `\label{16^k}` and the next display. Remark: "\eqref{16^k} is equivalent to the identity `sum 16^k(54k^2-36k-12)/((6k+1) C(6k,3k)) = -4 - pi/sqrt3`, which we cannot prove via integration." | conjecture; no "confirmed by" remark |
| T2, T9 | arXiv:2307.03086v10, Conj. 5.9 [2023-10-16] (i), eq. `\label{8k-3}` and the next display; MO 456779 eqs. (1), (2) | conjecture; MO: 0 answers, comments only |
| T3 | arXiv:2307.03086v10, Conj. 5.1 [2023-10-15] (i), second display: `\sum_{k=1}^\infty\f{(4k-1)3^k}{(2k-1)k^2\bi{2k}k\bi{3k}k}=2K` | conjecture; restated as conjectured in arXiv:2609.04176 (3 Sep 2026): "The author [S26] conjectured further that …" |
| T4, T8, T10 | arXiv:2307.03086v10, Conj. 5.25 [2023-11-13], (i) eq. `\label{2^15}`, (ii) both displays; Remark "rate of convergence is 27/32" | conjecture; arXiv:2603.29973v3: "was conjectured to have the value `32\sqrt2/\pi` by the author [(5.20)]{II}" (other items in that paper carry "confirmed by …"; this one does not) |
| T5 | arXiv:2307.03086v10, Conj. 4.33 [2025-02-28] (i), eq. `\label{-32}`; Remark "rate of convergence is -27/512" | stated as conjecture — **but see F1** |
| T7 | MO 506695 (6 Jan 2026) | 1 answer (C. Leibovici, 17 Jan 2026: partial fractions and 4F3 forms of the partial sum, "only related to the partial sum"), 1 comment; no proof. Not in Sun's arXiv lists |
| T11 | arXiv:2307.03086v10, Conj. 5.33 [2025-06-26] (i), eq. `\label{22k}`; Remark "rate of convergence is 1/64" | conjecture; arXiv:2603.29973v3: "was conjectured to have the value `64\sqrt2/\pi^2` by the author [(5.31)]{II}" |

All statements in RESULTS.md agree with the sources (polynomials, binomials, powers, right-hand sides, conjecture
numbers and dates).

### A.3 Findings

**F1 (FATAL for the novelty claim of T5). T5 is equation (7.10) of Chu (2020).**
W. Chu, Acta Math. Sci. 40B (2020), no. 2, 293–315, doi:10.1007/s10473-020-0201-y, Theorem 7.2 with `a = c = 1`
(`W_k(1,1,1) = (17+88k+145k^2+77k^3)/((1+2k)(3+4k)^2)`), Example 7.3:

    2G = sum_{n>=0} (-1/32)^n (1)_n (1)_n (1)_{3n} / ((3/2)_n (5/2)_{2n} (5/2)_{2n}) * (17 + 88n + 145n^2 + 77n^3)/9      (7.10)

with a proof (limit `e -> 1` of the transformation for the recurrence pattern [31-1] of the Gauss/Bailey Omega-sum,
`Theta(1,1) = 2G` by Whipple's quadratic transformation). Chu comments: "Two other identities may have appeared in
the vast literature on Catalan's constant, even though the author fails to locate them".
`audit2/a6_T5_is_Chu2020.py`: Chu's n-th term equals `-(1/3) x` (Sun's term at `k = n+1`) exactly for `n < 400`, the
term ratios agree identically (flint), `P(n+1) = 77n^3+145n^2+88n+17` for `P(k) = 77k^3-86k^2+29k-3`, and the sum
is `2G` to 60 digits. The polynomial `77n^3+145n^2+88n+17` is also the numerator of the certificate `R(n,0)` printed
in RESULTS.md, Theorem 5: the series `sum_n gamma(n,0) = 2G` of that proof *is* (7.10).
Consequences: T5 must be removed from the list of new results (or kept as "new proof of Chu 2020, (7.10); answers
Sun's Conj. 4.33(i), which Sun was apparently unaware was known"); `OPEN-LIST.md` row A9 and the sentence "A1, A5,
A7, A9 are within reach of Au's seeds … at shifts/base points nobody seems to have tried" are wrong for A9.
(Side remark, outside the eleven: Chu's (7.12), `2 ln 2 = sum (-1/2)^n n!(3n)!/(3+4n)! (34+101n+77n^2)/4`, is the base
series of Sun's `77k^2-53k+10` family in the same way.)

**R1 (REPAIRABLE). T2 is Guillera's conjecture (2018), not originally Sun's.**
With `k = n+1`, Sun's series is `(1/18) sum_{n>=0} (5+8n) [1, 3/4, 5/4; 3/2, 3/2, 3/2]_n 9^{-n}`, i.e. Ramanujan's
`sum (1+8n)(1/2)_n(1/4)_n(3/4)_n/(n!^3 9^n) = 2 sqrt3/pi` with all Pochhammer indices shifted by 1/2.
* J. Guillera, *Ramanujan series with a shift*, arXiv:1802.02023v4 (J. Aust. Math. Soc. 107 (2019) 367–380), Table 1
  (level 2: `q = e^{-pi sqrt6}`, `a = 1/(2 sqrt3)`, `b = 8/(2 sqrt3)`, `z = 1/9`) and Table 2 "Some conjectured
  values of `G_2(1/2,q)`": the entry is `pi/6`. That is exactly T2. The paper proves only level 4; levels 1, 2, 3 are
  "patterns observed by the author".
* H. Cohen, J. Guillera, arXiv:2101.12592, Section 6, table of `S(H,a,P;1/2)`, row 19 (`D = -3`, `a_0 = 2/3`,
  value `pi/3`): the same identity, stated as an observation.
* K. C. Au, arXiv:2602.08721v4, Example 8: the three-parameter summand `S(a,b,c;n)` obtained from
  `Seed10(a+3b-n+1/2, 2b+1, c+k+n)` equals `-(5+8n)[1,3/4,5/4;3/2,3/2,3/2]_n 9^{-n}` at `(a,b,c) = (-1/2,0,1/2)`.
  Au does not state T2 and works with the Taylor coefficients at `(0,0,0)` only, so T2 is not proved there; but the
  WZ family is the same one (the repository's own scan lists "exact match: … Seed10, Seed4" for A5).
All three facts are verified in `audit2/a7_T2_is_shifted_Ramanujan.py` (exactly, and numerically to 60 digits).
I found no proof of Guillera's level-2 shifted values anywhere, so the *proof* of T2 (and of T9, its `x`-derivative
companion) may well be new — but RESULTS.md must (a) attribute the conjecture to Guillera 2018 / Cohen–Guillera
2021 as well as Sun 2023, (b) say that it is the 1/2-shift of Ramanujan's `1/9` series, (c) cite Au's Example 8.
Fix: three sentences in the *Source* line of Theorem 2 and in the Summary; change "conjectured by Z.-W. Sun" to
"stated as open by Z.-W. Sun (T2 goes back to Guillera 2018)".

**R2 (REPAIRABLE, method). The novelty check was blind to equivalent forms.** "Full-text search … for the defining
polynomials: no third-party occurrence" cannot find (7.10) (index shifted by one, Pochhammer notation, polynomial
written in ascending order) or Guillera's table (the series is given by `(a, b, z)` and a shift). Fix: for every
claimed-new identity compute the canonical hypergeometric data (argument `z`, multiset of Pochhammer parameters mod
1, polynomial up to `k -> k+j`) and search for those; and read the Chu school's Abel-lemma/Omega-sum papers for the
seed in question (Gauss/Bailey 2F1(1/2): Chu 2020; Gauss 2F1(1)/Dougall 2H2: Chu 2011; Dougall 5F4: Chu–Zhang 2014;
very well-poised 7F6: Chu 2021), since "Bailey2F1 / Gauss2nd / Dougall" hits are exactly their territory.

**R3 (REPAIRABLE). Wrong citations for T3.**
* "MathOverflow 456444 (2023, no answer)": the question 456444 ("A curious series for `L(2,(-3/.))`", all
  revisions checked) contains `sum (42k^2-27k+4)(-9)^k/(k^3(2k-1) C(2k,k) C(3k,k) C(4k,2k)) = -6K` and its
  harmonic companion, **not** T3 (nor A3). That series was moreover *proved* by Au (arXiv:2312.14051, the example
  with `(-1/48)^n`). T3 does not occur in any of Sun's MathOverflow questions since 2023. The same wrong reference
  is in `OPEN-LIST.md` (row A3) and in "How novelty was checked (2)".
* "quoted as conjectural by Zuniga, arXiv:2601.12495v2": Zuniga quotes the *other* identity of Conj. 5.1(i),
  `sum (10k-3)/((2k-1)k^2 3^k C(2k,k)C(3k,k)) = K/2` (rate 1/81; row A3), not T3.
Fix: delete both references from Theorem 3; a correct supporting reference is Sun, arXiv:2609.04176 (Sep 2026),
which restates T3 as conjectured.

### A.4 Classification

| id | verdict | basis |
|---|---|---|
| T1 | OPEN as far as a serious search can tell | Sun III v6 (Feb 2026) conjecture with "cannot prove via integration"; no citing paper treats it; not in Le Lan's Beta-integral paper |
| T2 | UNCLEAR: statement known since Guillera 2018 (conjecture); no proof found | R1 |
| T3 | OPEN as far as a serious search can tell | still "conjectured" by Sun in Sep 2026; citations to fix (R3) |
| T4 | OPEN as far as a serious search can tell | not in Au 2023/2024/2026 (checked series by series; Au's only `27/32` example is the vanishing sum `sum (27/32)^n (1/6)_n(5/6)_n/n!^2 (10n-1)/(2n+1) = 0`, a different kernel); not in Chu 2020 (its `27/32` results are 4F3 evaluations and the same vanishing sum (5.14)); not in Chu 2011 (in that framework T4 is the half-integer pattern `[3/2,3/2,1,3]` at `(1/4,3/4,1,1)`; Chu's Example 92 at the same base point is a different series for `1/(pi sqrt2)` with rate `(3/8)^6`); Sun Apr 2026 "conjectured". Residual risk: Chu–Zhang 2014 unread (the repository's scan found no Dougall-5F4 hit for T4, which makes an overlap there unlikely) |
| T5 | **ALREADY PROVED**: W. Chu, Acta Math. Sci. 40B (2020), Example 7.3, eq. (7.10) | F1 |
| T6 | OPEN as far as a serious search can tell | companion of T1 |
| T7 | OPEN as far as a serious search can tell | MO 506695: partial-sum remarks only; no paper found |
| T8, T10 | OPEN as far as a serious search can tell | companions of T4; not in Au 2026 (no kernel with `27/32`) |
| T9 | OPEN as far as a serious search can tell (companion of T2, see R1) | MO 456779 eq. (2): no answer |
| T11 | OPEN as far as a serious search can tell | kernel `64^{-k}(1/4)_k^3(3/4)_k^3/((3/2)_k k!^5)` occurs in none of Au's lists (2312.14051, 2403.04555, 2602.08721v4), not in Chu's Dougall-7F6 dual-series papers (2103.07872, 2104.01599), not in Hou–Sun 2604.15172 / Hou–Mu 2604.22377; Sun Apr 2026 "conjectured". Residual risk: Chu, Ramanujan J. 55 (2021) (very well-poised Omega-sum = Dougall 7F6, the seed family of T11) could not be read; after F1 this risk is real and should be closed by someone with library access before T11 is announced as new |

---

## Part B — Proofs

### B.1 The WZ equations, re-derived independently (`audit2/a1_wz_certificates.py`)

The `.gp` files are parsed as text (only the assignments; the authors' `print` lines are ignored) and the identities
are re-stated in my script: `rho - 1 = R(n,k+1) sg - R` (T1, T2/T9, T3, T5, T6, T7, T11),
`(n+1) rho - n = Q(n,k+1) sg - Q` (T4/T10, T8), compatibility `rho(n,k) sg(n+1,k) = sg(n,k) rho(n,k+1)` (all), and
the triangle row `rho(n,n-1) - 1 = -rho(n,n-1) sg(n+1,n-1) - R(n,n-1)` (T7, T11). Two paths:
exact `Fraction` arithmetic at 300 random rational points per identity, and python-flint `fmpq_mpoly` with rational
functions kept as uncancelled numerator/denominator pairs (final numerator `== 0`). **All nine files pass on both
paths**, including the parameters `a, b` (T6) and `s` (T8). Negative controls (one coefficient changed in T4, T5,
T7, T11) are rejected by both paths. The polynomial `NUM` of T11 is textually identical in `certificates/T11.gp`,
`lean/RamanujanSeries/Sun1pi2.lean` and `proofs/verify_T11…py`.

### B.2 The pairs as numbers (`audit2/a2_exact_pairs.py`, exact rationals; 36 checks, all pass)

For T4, T11, T5, T7 the function `F` is computed from a closed form typed in from the statement in RESULTS.md
(not from the ratios), and on grids: the closed form has the printed ratios; the WZ equation holds **term by term**
(for T7, T11 in all four cases `k <= n-2`, `k = n-1`, `k = n`, `k > n`, up to `n = 44`, `k = 49`); the boundary
data (`h(n,0)`, `h(0,k)`, `Q(n,0)h(n,0) = t(n)/64`; `a(1,0) = t(0)/256`, `gam(n,0) = -t(n)/256`;
`psi(0,k)`, `gamma(n,0) = -t(n+1)/3`; `pi^2 Gam(n,0) = -t(n+1)/6`); the telescoped finite identities
(T4 for `N, K <= 12`; **identity (A) of T11 for `N = 1..44`**; the finite identity of T7 for `N <= 39`); and every
rational inequality used in a limit step, as an exact polynomial identity or coefficient positivity.

### B.3 T4 in full — correct and complete

Line by line:
1. `h > 0`, closed form `h(n,k) = 2^{-2k-5n}(1/2)_{3n}(1/2)_{2k+3n}/((1/2)_n (k+n)! n! (k+3n)!)`. Checked.
2. `1 - sigma - (3/4)(n+1)/(k+3n+1) = (4nk+4k+16n+1)/(16(k+n+1)(k+3n+1)) >= 0`. Exact.
3. `h(n,k) <= h(n,0)(m/(k+m))^c`, `c = 3(n+1)/4 <= m = 3n+1` (so each factor `1-c/(j+m)` is in `[0,1)`), using
   `sum_{j<k} 1/(j+m) >= log((k+m)/m)`. Correct.
4. `sum_k (m/(k+m))^c <= 1 + m/(c-1) <= 33/5` for `n >= 2` (it is 9 for `n = 1`; the true values of
   `sum_k h(n,k)/h(n,0)` are 2.946, 3.063, …, 3.835 (n = 40), 3.992 (n = 1000), tending to 4). Correct.
5. `0 <= Q <= k+3n+3` (coefficients), so `Q(n,K)h(n,K) -> 0` for `n >= 1` (exponent `c >= 3/2`). Correct. In fact
   `h(n,k) ~ C k^{-n-1}`; the weaker exponent suffices.
6. `n = 0`: `h(0,k) = (1/4)_k(3/4)_k/k!^2`, `Q(0,K)/K -> 1`, `(a)_K/K! ~ K^{a-1}/Gamma(a)`, limit
   `1/(Gamma(1/4)Gamma(3/4)) = 1/(sqrt2 pi)`. Correct; numerically `K (Q(0,K)h(0,K) - 1/(sqrt2 pi)) -> -0.2110`.
7. The finite identity with `u = n h` (so `u(0,k) = 0`; this is what makes the divergent row `n = 0` harmless),
   then `K -> oo` for fixed `N >= 1` (legitimate: finitely many `n`, each limit exists; `sum_k h(N,k)` converges
   for `N >= 1`), then `N -> oo` with `N sum_k h(N,k) <= 7N(27/32)^N` (`rho(n,0) <= 27/32` exactly: the cross
   difference is `1119744n^4+2944512n^3+2799360n^2+1126080n+158688`). There is no interchange of limits.
   Numerically `sum_{n<N} t(n)/64 + N sum_k h(N,k) = 1/(sqrt2 pi)` to 33 digits for `N = 1..8`.
Only remark (MINOR): "`<= 7N(27/32)^N`" is proved for `N >= 2`; irrelevant for the limit.

### B.4 T11 in full — correct; Part B is a complete Laplace argument at textbook level, not machine-checked

**Part A.** Positivity of `a(n,k)` (`4k-2n+9/2 = 1/2 - j`, `j = 2n-4k-4` even); `C_n` is the rational number
`-2^{4-8n}(3/4)_n^2(1/4)_n/((-1/4)_n(2n-2)!)` (using `Gamma(3/4)^2Gamma(1/4)/Gamma(-1/4) = -pi sqrt2/4`), so the
whole of Part A is rational arithmetic. Telescoping over `1 <= n <= N-1`, `0 <= k <= N`: left side
`T_N - a(1,0)`, right side `-sum gam(n,0)` because `gam(n,N+1) = 0`. All of this is confirmed exactly (B.2).
Since (A) is exact and `t(n) ~ c 64^{-n}`, `T_N` certainly converges; Part B only has to identify the limit.

**Part B**, every line (script `audit2/a4_T11_partB.py`; exact for (B1), 50 digits otherwise):
* (B1) The product formula for `r(k)` on `I` is an identity; each factor is nonnegative and decreasing, so `r`
  decreases on `I` and `a` is log-concave there. Below `I`: all four numbers `4k-2N+9/2+i` are `<= -1/2`; the lower
  bound `((N-k-1)(N+2k)^2/((N+k+1)(N-2k)^2))^2`, the polynomial identity
  `(N-k-1)(N+2k)^2-(N+k+1)(N-2k)^2 = 6N^2k-8k^3-2N^2-8k^2 >= 4N^2(k-1)`, the case `k = 0` (`N >= 4`), and the single
  remaining index `k = (N-3)/2` for odd `N` (numbers `-3/2,-1/2,1/2,3/2`) are all correct; the index sets do cover
  `0 <= k < min I` for both parities. Exact test for **every** index and every `3 <= N <= 600`: no failure
  (smallest `r-1` below `I`: 0.328; it decays like `N^{-1/2}`: 0.168, 0.103, 0.050 at `N` = 2000, 5001, 20000 — so
  the written argument, not numerics, is what proves it).
* (B2) `L(x*) = 0` iff `8x^3 = 6x`; `L'(x) = -4/(1-x^2) - 16/(4x^2-1)`, `L'(x*) = -24`; `sup|L''| = 388` on the
  window. `N |log r(k) - L(k/N)| <= 41.3, 41.11, 41.09, 41.09` on `[0.8,0.93]` for `N = 10^3..10^6` (uniform, as
  claimed); the constant in `-24j/N + O(1/N + j^2/N^2)` for `|j| <= J+1` is `65, 47, 42.3, 41.5`; the summed form
  has `sup|…| N^{1/5} = 25, 17.8, 15.1, 14.7`. Correct. (The `O(N^{-1/5})` is honest but slow — about `14.7 N^{-1/5}`
  at the edge of the window; harmless because the Gaussian weight there is `e^{-12N^{1/5}}`.)
* (B3) Correct (unimodal sum versus integral, error at most the maximum 1).
* (B4) Directions of both inequalities are right (they use only that `r` decreases on `I`); `r(k0+J) <= e^{-12J/N}`
  and `r(k0-J-1) >= e^{12J/N}` hold for every `N` in `(164, 30000]` and at `10^5, 10^6`; `1-e^{-u} >= u/2` is used
  for `u = 12J/N <= 1.59`; `a(min I)/a(k0) <= e^{-4N^{3/5}}` holds with a wide margin (the true size is
  `(4/27)^{N(1+o(1))}`); measured tails are far below the claimed `N^{2/5}e^{-11N^{1/5}}`.
* (B5) I re-derived `phi` by hand: the `N log N` terms cancel (`(4x+2)-2(1-x)-2(1+x)-(4x-2) = 0`), the linear
  terms cancel, `log C_N` contributes `-10 N log 2`, and the `log 2` coefficient is `2(2x+1)-2(2x-1)-10 = -6`.
  `phi(x*) = 0` exactly by the argument in the text; `x(1-x)/((1+x)(2x-1)^4) = sqrt3/8` at `x*` (by hand), so
  `psi0(x*) = sqrt3/(sqrt2 pi^{5/2})` and `psi0(x*) sqrt(pi/12) = sqrt2/(4pi^2)`. Uniform remainder:
  `N |log a - (N phi - (1/2)log N + log psi0)| <= 7.76` at five points of the window for `N = 10^2..10^6`.
* Combination: no circularity ((B4) uses (B1), (B2); (B5) is independent; the value of the series is never used).
  `T_N/(sqrt2/(4pi^2)) - 1` is below `10^{-43}` for `N = 10^2..10^6` (as it must be, by (A)).

I have no objection to "proved" for T11, with the existing caveat "written out in full, not machine-checked".
What is and is not formal: Lean has only the two rational identities (`wz_rat`, `wz_boundary`) under
non-vanishing hypotheses — not the triangle telescoping, not `gam(n,0) = -t(n)/256`, not Part B.

### B.5 Two further proofs in full: T5 and T7

**T5** (correct; not new, F1). `psi(0,k) = 4^k/((2k+1)^2 C(2k,k))`, `sum = 2G`
(`int_0^{pi/2} theta/sin theta = 2G`; both to 30 digits). `1-sigma-(n+3/4)/(k+2n+3/2) =
(8nk+6k+6n^2+12n+5)/(2(2k+2n+3)(2k+4n+3))`; `c = n+3/4 <= m = 2n+3/2`; `|rho(n,0)| <= 1/8`; `0 <= R <= 2k+2n+2`;
`1 + m/(c-1) <= 6` for `N >= 1`; `gamma(0,K) sqrt(K) -> 0.8862` so `gamma(0,K) = O(K^{-1/2})`. The row `n = 0`
converges although `c = 3/4 < 1` there (the bound is not needed for `n = 0`). After `K -> oo`:
`sum_{n<N} gamma(n,0) = 2G - sum_k psi(N,k)` to 30 digits, `N = 1..6`. The remark about the "first-order coefficient
in nu" is not used in the proof, which is good, because it is not justified in the text.

**T7** (correct). All four cases of the WZ equation and the finite identity
`sum_{m<=N} t(m) = 6 pi^2 sum_{k<N} Phi(N,k)` hold exactly. The two bounds of the limit step are only asserted in
RESULTS.md; both are true for all indices: with `n = k+2+v`, `den - 24 num` of `|Phi(n,k+1)/Phi(n,k)|` has
nonnegative coefficients in `(k,v)`, and with `n = m+1`, `den - 6 num` of `rho(n,0)` has nonnegative coefficients
(`a2`, symbolic). Hence `|sum_{k<N} Phi(N,k)| <= (24/23)|Phi(N,0)| <= C 6^{-N}`.

T2 (inequalities `sigma <= (2k+1)/(2k+5)` for `n >= 1`, `rho(n,0) <= 1/8`, `R <= k+2`: symbolic; boundary data)
and T3 (WZ for real `x`; `G'(n,0) = -(3/(4pi))t(n+1)` for `n <= 6`; `D_0 = -3K/(2pi)`; `I = -4K`;
`sum_{n<N} G'(n,0) = D_0 - D_N` and the bound on `D_N` for `N <= 4`) were read and spot-checked: no error.

### B.6 T6, T8, T9, T10 — correct in substance, one step is a sketch (R4)

The final identities are right (B.7), the parametrised WZ equations hold identically in the parameter (B.1), the
digamma/derivative formulas are right (`a3`), and the constancy claims hold numerically to 30+ digits
(`sum_n Q_s(n,0)h_s(n,0) = 1/(4 sqrt pi)` for `s = 0, ±1/20, ±1/10`; `sum_n Q(n,b)h(n,b) = 1/(sqrt2 pi)`;
`sum_n g(n,c) = (psi(c+2/3)-psi(c+1/3))/3`; `sum_n G(n,0;a,b) = sum_k F(0,k;a,b)`; `sum_k F(0,k;a,0) = 2^{-1/2}`).
But:
* "Exactly as in Theorem 4 (the bounds are uniform for `|s| <= 1/10`)", "the proof of Theorem 2 applies verbatim to
  the rows `x = k+c`", "the bounds of Theorem 1 persist … uniformly" are assertions. The inequalities of Theorems
  1, 2, 4 were verified as "nonnegative coefficients in `(n,k)`", which says nothing about `k` replaced by `k + c`
  with `c < 0`, or about the deformed ratios `sigma_s`, `rho_s`, `Q_s`. They are true and easy, but they are not in
  the text or in the scripts.
* "the series converge uniformly with all terms analytic in `c`, [hence] can be differentiated termwise" is not a
  valid principle on a real interval (uniform limits of real-analytic functions need not be differentiable).
  **Fix (one line each):** the differentiated series converges uniformly, because
  `d/dc log(term_n) = O(log(n+2))` uniformly on `|c| <= 1/10` (these are exactly the digamma expressions the text
  computes) while the terms are `O(n^A rho^n)` with `rho < 1`; together with convergence at one point this gives
  termwise differentiability. (Alternatively: the same bounds on a complex disc and Weierstrass.)
Until this is written, "proved" for T6, T8, T9, T10 means "proof with one routine step sketched".

### B.7 The constants (`audit2/a5_series_200digits.py`)

Direct summation of the series as printed by Sun (binomials by `math.comb`, harmonic numbers by running sums,
mpmath at 260 digits; 3100 terms for the three series of rate 27/32):

| T1 | T2 | T3 | T4 | T5 | T6 | T7 | T8 | T9 | T10 | T11 |
|---|---|---|---|---|---|---|---|---|---|---|
| 230 | 234 | 235 | 229 | 229 | 232 | 230 | 229 | 233 | 228 | 234 |

(digits of agreement with the closed form; `K = (zeta(2,1/3)-zeta(2,2/3))/9`).

---

## Part C — Lean

* Pins: `lean-toolchain` = `leanprover/lean4:v4.33.0-rc2`, Mathlib `51e6992efd06126df61a496bebf8f49482a4e129`
  (as stated in PROGRESS.md; **not** v4.36.0, so `monotonicity/lean/.lake/packages` does not apply). A build with
  exactly these pins exists in `monotonicity/lean-race/.lake/packages`; I copied it and built a byte-identical copy
  of the `lean/` sources in a scratch directory (nothing written under `lean/`, no download, Mathlib not rebuilt,
  2 cores, `nice`).
* `lake build`: "Build completed successfully (8700 jobs)", no warnings, no `sorry`/`axiom`/`native_decide` in the
  sources. `#print axioms SunSeries16.sun_series_16`: `[propext, Classical.choice, Quot.sound]`.
  `leanchecker RamanujanSeries.Sun16`: exit 0. Log: `audit2/c1_lean_rebuild.log.txt`.
* The statement, with nothing to unfold (the theorem does not mention the auxiliary definitions):

      SunSeries16.sun_series_16 : HasSum
        (fun k => (3 * ↑k + 2) ^ 2 * 16 ^ k / ((2 * ↑k + 1) * (6 * ↑k + 1) * (6 * ↑k + 5) * ↑(3 * k).centralBinom))
        (Real.pi / (2 * √3))

  `k : ℕ`, casts to `ℝ`, `Nat.centralBinom (3k) = C(6k,3k)`, right-hand side `π/(2√3)`. It is a `HasSum`, so
  summability is part of the statement (no `tsum`-of-a-non-summable-series loophole). In
  `audit2/c1_lean_statement_check.lean` I restate it with `Nat.choose (6*k) (3*k)` written out and derive the
  `Summable` and `∑' … = π/(2√3)` forms from the project's theorem; it compiles.
  **This is identity T1.** The claim "T1 is proved completely in Lean" is accurate.
* `Sun1pi2.lean` proves what PROGRESS.md says and no more: `wz_poly`, `wz_rat`, `wz_boundary_poly`, `wz_boundary`
  for real `n, k` under non-vanishing of the denominators.

---

## Part D — Wording of RESULTS.md

* **Summary 1, "Eleven identities conjectured by Z.-W. Sun (2023–2026), all still stated as open in the latest
  versions of his lists"**: literally true of Sun's lists, but misleading: T5 was a theorem since 2020 (F1) and T2
  is Guillera's conjecture (R1). Should read "ten" (or "eleven, one of which — T5 — turned out to be Chu 2020,
  (7.10)"), with the T2 attribution corrected.
* **Summary 2, "every analytic step is written out below with its reference"**: not for T6, T8, T9, T10 (R4).
* **T3 source line and "How novelty was checked (2)"**: R3.
* M1 (MINOR). "T11 … a new Ramanujan-type series for 1/π²": the *series* is Sun's (2025); what is new is the proof.
  Say "a proof of Sun's series".
* M2 (MINOR). T7: the constants `1/24` and `1/6` are asserted; add "(coefficientwise, after `n = k+2+v`, resp.
  `n = m+1`)" or weaken to any constant `< 1` (`5/11` is immediate).
* M3 (MINOR). T4: `N sum_k h(N,k) <= 7N(27/32)^N` is shown for `N >= 2`.
* M4 (MINOR). README/PROGRESS refer to `proofs/verify_*.out`, `certificates/check_all.out`, `lean/Axioms.out`,
  `survey/verify_*.out`; none of these is in the repository (`*.out` is git-ignored). Likewise the earlier audit
  says "all scratch scripts and their outputs are in this directory" (`t11_partA.py`, `t11_B1.py`, …): they are
  not there, so that audit cannot be re-run. (Its numbers that I recomputed — 41.09, 164, 2.95/3.06/3.84/3.96, the
  values of `256 T_N` — all agree with mine.)
* M5 (MINOR). `OPEN-LIST.md`, note "Why existing methods had not reached the rows we prove": false for A9 (= T5),
  and for A5 (= T2) the family is the one in Au's Example 8.
* M6 (MINOR). The brief's "arXiv:2209.01731" is indeed not Au; Au's first WZ paper is arXiv:2212.02986 (worth
  saying, since PROGRESS only says the number is wrong).

---

## Final table

| id | identity (value) | novelty | proof | machine-checked |
|---|---|---|---|---|
| T1 | `(3k+2)^2 16^k/((2k+1)(6k+1)(6k+5)C(6k,3k))`, `pi/(2 sqrt3)` | open as far as found | correct, complete | **all of it**: Lean 4 + Mathlib, `HasSum`, rebuilt here, axioms `propext, Classical.choice, Quot.sound` |
| T2 | `(8k-3)C(4k,2k)/(k(4k-1)9^k C(2k,k)^2)`, `sqrt3 pi/18` | unclear: Guillera's 2018 conjecture (shifted Ramanujan series), no proof found | correct, complete | WZ equation + inequalities exactly (sympy by authors; flint/Fraction here) |
| T3 | `(4k-1)3^k/((2k-1)k^2 C(2k,k)C(3k,k))`, `2 L_{-3}(2)` | open as far as found (citations wrong) | correct, complete | WZ equation exactly; analytic steps numerically |
| T4 | `(360k^3+612k^2+230k+15)C(3k,k)C(6k,3k)^2/((3k+1)(3k+2)2^{15k})`, `32 sqrt2/pi` | open as far as found (residual risk: Chu–Zhang 2014 unread) | correct, complete | WZ equation, finite identity, all inequalities exactly |
| T5 | `(77k^3-86k^2+29k-3)(-32)^k/(k^3(3k-1)(3k-2)C(3k,k)C(4k,k)^2)`, `-6G` | **already proved: Chu, Acta Math. Sci. 40B (2020), (7.10)** | correct, complete | as T4 |
| T6 | T1 with `H_k`, `(5/2)K - pi log432/(6 sqrt3)` | open as far as found | correct; uniformity/termwise differentiation sketched | WZ with parameters exactly; rest numerically |
| T7 | Sun's vanishing series (MO 506695), `0` | open as far as found | correct, complete | WZ on the triangle + finite identity + both bounds exactly |
| T8 | T4 with `H_{2k}-H_k`, `48 sqrt2 log2/pi` | open as far as found | correct; uniformity sketched | WZ with parameter exactly; rest numerically |
| T9 | T2 with `5H_{2k-1}-4H_{k-1}`, `3K/2` | open as far as found (companion of Guillera's conjecture) | correct; uniformity sketched | WZ exactly; rest numerically |
| T10 | T4 with `4H_{6k}-3H_{3k}-H_k`, `192 sqrt2 log2/pi` | open as far as found | correct; uniformity sketched | WZ exactly; rest numerically |
| T11 | `(1344k^3+944k^2+156k+9)C(2k,k)^2C(4k,2k)^3/((2k+1)2^{22k})`, `64 sqrt2/pi^2` | open as far as found (residual risk: Chu, Ramanujan J. 55 (2021) unread) | correct; Part A exact, Part B Laplace argument refereed line by line, numerics to `N = 10^6` | Lean: the two rational WZ identities only; Part A exactly by CAS (identity (A) for `N <= 44`); Part B not machine-checked |

All eleven: series re-summed to >= 228 digits against the closed form.

## Files in `audit2/`

| file | content |
|---|---|
| `a1_wz_certificates.py` | independent parser + two-path check of all certificates, with negative controls |
| `a2_exact_pairs.py` | T4, T11, T5, T7 (and T2 inequalities) in exact rational arithmetic: closed forms, WZ equation term by term, boundary data, finite identities, inequalities |
| `a3_limit_steps.py` | numerical tests of every limit step (T4, T10, T8, T5, T7, T3, T2/T9, T6) |
| `a4_T11_partB.py` | Theorem 11 Part B: (B1) exactly for all `N <= 600`, (B2)–(B5) with measured O-constants up to `N = 10^6` |
| `a5_series_200digits.py` | all eleven series to >= 228 digits |
| `a6_T5_is_Chu2020.py` | T5 = Chu 2020 (7.10), exactly |
| `a7_T2_is_shifted_Ramanujan.py` | T2 = Guillera 2018 Table 2 = Cohen–Guillera #19 = Au 2026 Ex. 8 at `(-1/2,0,1/2)` |
| `c1_lean_statement_check.lean`, `c1_lean_rebuild.log.txt` | Lean rebuild, axioms, restated theorem |
