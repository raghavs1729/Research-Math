# RESULTS — WZ proofs, with certificates, of conjectured Ramanujan-type series

**Summary (2026-11-24; revised after the independent audit `AUDIT-INDEPENDENT.md`).**
1. Eleven conjectured identities from Z.-W. Sun's lists (2023–2026), one more from Guillera's tables (T12) and one from Cohen–Guillera's (T13), are **proved** here by the Wilf–Zeilberger method, each with an explicit certificate `R(n,k)` and a script checking the WZ equation as an exact rational-function identity. **Ten** of the thirteen were not found in any source we could read and are, with the caveats of item 3, new *proofs of conjectured series* (the series themselves are Sun's, resp. Guillera's — not ours; for two of these, T1 and T7, the WZ pair lies in classical families treated in papers we could not read, so priority is not established); **two (T2 and T5) turned out to be already in the literature**, in papers of W. Chu that neither Sun nor Guillera nor our first survey had noticed, **and one (T13) is a known statement with a proof indicated on MathOverflow**:
   - T1 `sum_{k>=0} (3k+2)^2 16^k/((2k+1)(6k+1)(6k+5) C(6k,3k)) = pi/(2 sqrt3)` (Sun, arXiv:2506.01870v6, Conj. 3.2);
   - T2 `sum_{k>=1} (8k-3) C(4k,2k)/(k(4k-1) 9^k C(2k,k)^2) = sqrt3 pi/18` (Sun, arXiv:2307.03086v10, Conj. 5.9(i); Guillera's conjecture of 2018, arXiv:1802.02023, Table 2) — **NOT NEW: this is Example 19 of W. Chu, Kodai Math. J. 41 (2018) 512–530** (`pi sqrt3 = sum [1,3/4,5/4;3/2,3/2,3/2]_k (5+8k)/9^k`, by trisection of a cubic transformation of Dougall's sum; found by us on 2026-11-24 after the audit); ours is a second proof;
   - T3 `sum_{k>=1} (4k-1) 3^k/((2k-1) k^2 C(2k,k) C(3k,k)) = 2 L_{-3}(2)` (Sun, arXiv:2307.03086v10, Conj. 5.1(i); restated as open in arXiv:2609.04176);
   - T4 `sum_{k>=0} (360k^3+612k^2+230k+15) C(3k,k) C(6k,3k)^2/((3k+1)(3k+2) 2^(15k)) = 32 sqrt2/pi` (ibid., Conj. 5.25(i)) — a Ramanujan-type series for 1/pi;
   - T5 `sum_{k>=1} (77k^3-86k^2+29k-3)(-32)^k/(k^3(3k-1)(3k-2) C(3k,k) C(4k,k)^2) = -6G` (ibid., Conj. 4.33(i)) — **NOT NEW: this is eq. (7.10) of W. Chu, Acta Math. Sci. 40B (2020) 293–315** (index shifted by one); ours is only a second proof;
   - T6 `sum_{k>=0} (3k+2)^2 16^k H_k/((2k+1)(6k+1)(6k+5) C(6k,3k)) = (5/2) L_{-3}(2) - pi log(432)/(6 sqrt3)` (arXiv:2506.01870v6, Conj. 3.2);
   - T7 `sum_{k>=1} C(2k,k) 16^k (1656k^4-2604k^3+1158k^2-229k+15)/(k^2(2k-1)^2(6k-1)(6k-5) C(6k,3k) C(3k,k)) = 0` (Sun, MO 506695);
   - **T11 `sum_{k>=0} (1344k^3+944k^2+156k+9) C(2k,k)^2 C(4k,2k)^3/((2k+1) 2^(22k)) = 64 sqrt2/pi^2`** (arXiv:2307.03086v10, Conj. 5.33(i)) — Sun's conjectured Ramanujan-type series for 1/pi^2; its proof needs, besides the WZ pair, a Laplace-method asymptotic for an explicit log-concave sum (written out in full, not machine-checked);
   - **T12 `sum_{n>=0} (28n+17) (1)_n(3/4)_n(5/4)_n/((3/2)_n^3 (-48)^n) = 32 log(27/16)`**, i.e. `sum_{k>=1} (28k-11) C(4k,2k)/(k(4k-1)(-48)^k C(2k,k)^2) = (4/3) log 2 - log 3` — another **conjecture of Guillera (2018)** (half-shift of Ramanujan's `-1/48` series; arXiv:1802.02023, Tables 1–2), proved after the audit; not in Sun's lists, not in Chu's Kodai paper;
   - T13 `sum_{n>=0} P(n+1/2) n!^16 (4n+2)!/(2n+1)!^10 = 2048 pi^4/15` (`P` = Cullen's quartic; the half-shift of Cullen's `1/pi^4` series, Cohen–Guillera 2021, Sect. 6, item 13) — **NOT A NEW STATEMENT**: it is eq. (2) of D. Chen's MathOverflow question 486961 (1 Feb 2025), `sum 30(4k+1)(5460k^4+13499k^3+12601k^2+5265k+831) C(4k,2k)/((2k+1)^9 C(2k,k)^8) = (4 pi)^4`, and J. Guillera answered there that it follows from Au's WZ pair for Cullen's series "replacing n with n+1/2"; no written proof was found, and ours shows that the substitution alone gives `0 = 0` and that the identity is the next order, which needs the Sun–Au–Wei harmonic identity. We found the MO post only after writing the proof (it was in our own SURVEY-C in binomial form);
   - T8, T10 the two harmonic-number companions of T4 (values `48 sqrt2 log2/pi`, `192 sqrt2 log2/pi`; Conj. 5.25(ii)), and T9 the companion of T2 with `5H_{2k-1}-4H_{k-1}` (value `3 L_{-3}(2)/2`; Conj. 5.9(i)).
2. **Machine-checked:** T1 **and T4** are proved completely in Lean 4 + Mathlib, no `sorry`, axioms `propext, Classical.choice, Quot.sound` only: `lean/RamanujanSeries/Sun16.lean` (theorem `SunSeries16.sun_series_16`) and `lean/RamanujanSeries/Sun1pi.lean` (theorem `SunSeries1pi.sun_series_1pi : HasSum (fun k => (360k^3+612k^2+230k+15) * (3k).choose k * (3k).centralBinom^2 / ((3k+1)(3k+2) 2^(15k))) (32 * sqrt 2 / pi)`; WZ equation, finite identity, decay of the rows `n >= 1`, the limit `1/(Gamma(1/4)Gamma(3/4)) = sqrt2/(2 pi)` via Euler's limit formula and the reflection formula, identification with Sun's summand). For T2, T3, T5–T11 the WZ equations and all rational-function inequalities used (including, for T6, T8, T9, T10, the bounds uniform in the deformation parameter, `proofs/verify_R4_uniform_bounds.py`) are checked exactly by scripts; the analytic steps are written out below with references and are **not** formalized (for T11 only the certificate identities are in Lean, `Sun1pi2.lean`).
3. **Novelty caveat.** Two independent audits found all proofs correct, but the second one found that our first novelty check was blind to index shifts and Pochhammer notation (T5), and our own follow-up found that it was also blind to *multisection* (T2 is a trisected Dougall acceleration in Chu 2018). The check was redone on canonical hypergeometric data including bisections and trisections (see the last section). Two papers of the Chu school could still not be read (Chu–Zhang, Math. Comp. 83 (2014); Chu, Ramanujan J. 55 (2021)); their method is known to us from Guillera, arXiv:1611.04385, and Chu's Kodai paper, and their seed families (Omega, Phi, 2H2) were enumerated independently (`code/chu_family_test.py`) on canonical data, plain, bisected and trisected: **T11 occurs in none of them in any form; T4 only in bisected form (Phi and 2H2 families); T12 and T3 only in trisected form (Omega); T1 and T7 lie inside these classical families at small patterns**, so for T1 and T7 priority is *not established* (see the table in the last section). "New" therefore means "not found after this search", not a guarantee; given that two of twelve were already known, the reader should treat the remaining ten with corresponding caution.
4. **Not proved (honest failures).** The famous open rational 1/pi^2 series (7^-4, 80^-3, (3/5)^6, 2^-10/1640, -1/48/252, golden) have no signature hit in any of our 169 seeds; Gourevitch's 1/pi^3 and Cullen's 1/pi^4 series were *already proved* by Au in 2023 (the brief was out of date). Sun's A2 (20 zeta(3)), A3 (K/2) and ~70 harmonic-number companions remain open.
5. Method: Au's WZ seeds (arXiv:2212.02986, 2312.14051, 2602.08721), extended by 149 one-parameter seeds generated automatically from Ebisu's tables of strange 2F1 evaluations and 5 classical summations at z != 1; own Gosper/Zeilberger engine on PARI (`code/wz.py`). Details: `PROGRESS.md`, `OPEN-LIST.md`.

---

## 0. Conventions and the two general facts used

`(a)_m` rising factorial; `C(n,k)` binomial; `K = L_{-3}(2) = sum_{m>=0} (1/(3m+1)^2 - 1/(3m+2)^2)`;
`G` Catalan's constant; `Cl_2(t) = sum_{m>=1} sin(mt)/m^2`.

**Telescoping lemma.** If `u, v : N x N -> R` satisfy `u(n+1,k) - u(n,k) = v(n,k+1) - v(n,k)` for all
`n,k >= 0`, then for all `N, K >= 0`

    sum_{k<K} u(N,k) - sum_{k<K} u(0,k) = sum_{n<N} v(n,K) - sum_{n<N} v(n,0).

(Sum the equation over `0 <= k < K`, `0 <= n < N`; both sides telescope. Formalized as
`SunSeries16.finite_identity`.) This is the finite form of Au's Proposition 2.1 (arXiv:2312.14051).

**How the pairs were found.** Each pair is `F(n,k) = seed(a_0 + A n, ..., k + k_0 + K n)` for a hypergeometric
term `seed(a,...,k)` whose sum over `k` does not depend on the parameters (a "WZ seed" in the sense of
Au), and `G = R F` with `R` computed by Gosper's algorithm. The search (`code/run_scan.py`, `code/match_hits.py`)
is exhaustive over a box of shifts and base points; see `PROGRESS.md`. None of this is needed to *check* the
proofs: each proof below only uses the explicit `R` and the displayed ratios.

In each proof `sigma(n,k) = f(n,k+1)/f(n,k)` and `rho(n,k) = f(n+1,k)/f(n,k)`; the WZ equation
`f(n+1,k) - f(n,k) = g(n,k+1) - g(n,k)` with `g = R f` is equivalent to the rational identity

    rho(n,k) - 1 = R(n,k+1) sigma(n,k) - R(n,k),                                   (WZ)

which is what the scripts verify (`sympy.cancel(...) == 0`), together with the compatibility
`rho(n,k) sigma(n+1,k) = sigma(n,k) rho(n,k+1)` that makes `f` well defined from its ratios.

---

## Theorem 1.  `sum_{k>=0} (3k+2)^2 16^k / ((2k+1)(6k+1)(6k+5) C(6k,3k)) = pi/(2 sqrt 3)`

*Source.* Sun, "New series involving binomial coefficients (III)", arXiv:2506.01870v6, Conjecture 3.2 (first
identity); Sun remarks that an equivalent form "cannot [be proved] via integration".

*Pair.* `f(0,0) = 1`,

    sigma(n,k) = -(k+3n+1)(2k+1) / (3 (2k+6n+3)(k+n+1)),
    rho(n,k)   = 2 (n+1)(k+3n+1)(k+3n+2)(k+3n+3) / ((2k+6n+3)(2k+6n+5)(2k+6n+7)(k+n+1)),
    R(n,k)     = 3 (k+3n+2)^2 / ((2k+6n+3)(2k+6n+5)).

Closed form: `f(n,k) = (-1/3)^k 4^-n (1/2)_n (3/2)_{2n} (1)_n^2 (1)_{k+3n} (1/2)_k / ((1)_{2n} (3/4)_n (5/4)_n (3/2)_{k+3n} (1)_{k+n})`.
Origin: Ebisu's evaluation `2F1(2a, 1-a; 2a+1/2; -1/3) = sqrt(pi) Gamma(2a+1/2)/(3^a Gamma(a+1/2)^2)`
(A. Ebisu, *Special values of the hypergeometric series*, Mem. AMS 248 (2017), formula (1,3,3-1)(xiv)) as a seed,
`a -> 1/2 + n`, `k -> k + n`.

*Proof.* (WZ) holds. `|sigma| <= 1/3`, `0 < rho(n,0) <= 1/3`, `0 <= R <= 1`, hence
`|f(n,k)| <= 3^{-n-k}` and `|g(n,k)| <= 3^{-n-k}`. In the telescoping lemma let `K -> oo`:
`sum_{n<N} g(n,0) = sum_k f(0,k) - sum_k f(N,k)`, and `|sum_k f(N,k)| <= (3/2) 3^{-N} -> 0`. Now
`f(0,k) = (-1/3)^k/(2k+1)`, so `sum_k f(0,k) = sqrt3 arctan(1/sqrt3) = pi/(2 sqrt3)`. Finally
`f(n,0) = 16^n/((6n+1) C(6n,3n))`, so `g(n,0) = R(n,0) f(n,0)` is the summand. ∎

*Machine check.* Complete: `lean/RamanujanSeries/Sun16.lean` (every step above). Script: `proofs/verify_T1_sun_16k.py`.

---

## Theorem 2.  `sum_{k>=1} (8k-3) C(4k,2k) / (k (4k-1) 9^k C(2k,k)^2) = sqrt(3) pi/18`

**Not a new result.** With `k = n+1` this is Example 19 of W. Chu, *q-series reciprocities and further pi-formulae*,
Kodai Math. J. 41 (2018) 512–530 (received October 2017): `pi sqrt3 = sum_{k>=0} [1, 3/4, 5/4; 3/2, 3/2, 3/2]_k (5+8k)/9^k`,
obtained there from the cubic transformation with iteration pattern [30111] of the very well-poised sum
`Omega(a;b,c,d,e)` (Theorem 20 of that paper, the `q -> 1` case of Chu–Zhang, J. Difference Equ. Appl. 18 (2012)) at
`a = c = 1, b = 1/2, d = 1/3, e = 2/3`, followed by trisection. In our language: with `n = 3m` the pair below is the
integer pattern `[3;1,1,1,0]` of Dougall's sum at `(1; 1/3, 2/3, 1, 1/2)`, where `Omega = sum 2/((3k+1)(3k+2))`. We found
this only after the second audit, with `code/chu_family_test.py` (which reports the *trisected* canonical data of
T2 inside the Dougall family) and then by reading the Kodai paper (open access on J-STAGE). Guillera's paper
(February 2018) lists the value as conjectured, Sun (2023) as an open conjecture; neither cites Chu's example.
The WZ proof below is kept as a second proof; Theorem 9 (the harmonic-number companion) is not in Chu's paper.

*Source and attribution.* With `k = n+1` the series is `(1/18) sum_{n>=0} (5+8n) (1)_n (3/4)_n (5/4)_n / ((3/2)_n^3 9^n)`,
i.e. Ramanujan's `sum (1+8n) (1/2)_n (1/4)_n (3/4)_n/(n!^3 9^n) = 2 sqrt3/pi` with every Pochhammer index shifted by 1/2.
In this form it was **conjectured by J. Guillera**, *Ramanujan series with a shift*, arXiv:1802.02023 (J. Aust. Math. Soc.
107 (2019) 367–380): Table 1 (level 2, `q = e^{-pi sqrt6}`, `z = 1/9`) and Table 2 ("Some conjectured values of
`G_2(1/2,q)`": the value `pi/6`); that paper proves only the level-4 analogues. It reappears as an observation in
H. Cohen and J. Guillera, arXiv:2101.12592, Section 6, and as Conjecture 5.9(i) (first identity) of Z.-W. Sun,
arXiv:2307.03086v10, and MathOverflow 456779 (2023, no answer). *Relation to Au's work:* the three-parameter summand
`S(a,b,c;n)` of Au's Example VIII (arXiv:2602.08721, from `Seed10(a+3b-n+1/2, 2b+1, c+k+n)`) specialises at
`(a,b,c) = (-1/2,0,1/2)` to a constant multiple of our summand (checked in `audit2/a7_T2_is_shifted_Ramanujan.py`).
Au uses that family only through its Taylor coefficients at `(0,0,0)` (harmonic-number extensions of Ramanujan's
series) and neither states nor proves the identity at the half-integer point, where one needs the convergence
analysis below and the evaluation `sum_k 1/((3k+1)(3k+2)) = pi/(3 sqrt3)` of the other boundary. (Superseded: the identity is in Chu 2018, see above.)

*Pair.* `f(0,0) = 1/2`,

    sigma(n,k) = (2k+1)(2k+n+3)(3k+n+1)(3k+n+2)(3k+n+3) / ((2k+n+1)(2k+2n+3)(3k+2n+3)(3k+2n+4)(3k+2n+5)),
    rho(n,k)   = (n+1)(4n+3)(4n+5)(2k+n+2)(3k+n+1) / (9(2n+1)(2k+n+1)(2k+2n+3)(3k+2n+3)(3k+2n+4)),
    R(n,k)     = (3k+2n+2)(6k+8n+5) / (9(2n+1)(2k+n+1)).

Closed form: `f = -(pi/sqrt3) F`, `F(n,k) = (2k+n+1) 3^(1/2-2n) n! Gamma(k+1/2) Gamma(3k+n+1) / (Gamma(-2n-1/2) Gamma(k+n+3/2) Gamma(n+1/2) Gamma(3k+2n+3))`
`= Seed4(-n, -1/2-2n, k+1/2+n)` (Au's Seed4: Dougall's 5F4 with `b -> b/3, c -> b/3+1/3, d -> b/3+2/3`).

*Proof.* (WZ) holds; `f > 0`, `g > 0`. `f(0,k) = 1/((3k+1)(3k+2))`, so
`sum_k f(0,k) = sum_k (1/(3k+1) - 1/(3k+2)) = int_0^1 dx/(1+x+x^2) = pi/(3 sqrt3)`.
For `n >= 1`: `sigma(n,k) <= (2k+1)/(2k+5)`, hence `f(n,k) <= 3 f(n,0)/((2k+1)(2k+3))` and
`sum_k f(n,k) <= (3/2) f(n,0)`; `rho(n,0) <= 1/8` gives `f(N,0) <= 8^{-N}/2 -> 0`.
`R(n,k) <= k+2`, so `g(n,K) <= (K+2) 3 f(n,0)/((2K+1)(2K+3)) -> 0` for `n >= 1`, and
`g(0,K) = (6K+5)/(9(2K+1)(3K+1)) -> 0`. The telescoping lemma with `K -> oo`, then `N -> oo`, gives
`sum_n g(n,0) = pi/(3 sqrt3)`. Finally `g(n,0) = 2 t(n+1)` with `t` the summand (equal at `n=0`, same ratio). ∎

*Check.* `proofs/verify_T2_sun_9k.py` (exact: compatibility, (WZ), boundary data, the three inequalities as
"numerator and denominator have nonnegative coefficients"; numerics 110 digits).

---

## Theorem 3.  `sum_{k>=1} (4k-1) 3^k / ((2k-1) k^2 C(2k,k) C(3k,k)) = 2K`

*Source.* Sun, arXiv:2307.03086v10, Conjecture 5.1(i) (second identity); restated as a conjecture in Sun,
arXiv:2609.04176 (3 Sep 2026). (An earlier version of this file also cited MathOverflow 456444 and Zuniga,
arXiv:2601.12495; both concern *other* series — MO 456444 is the `(42k^2-27k+4)(-9)^k` series, proved by Au in
arXiv:2312.14051, and Zuniga quotes the first identity of Conj. 5.1(i), our row A3 — and have been removed.)

*Pair with a continuous second variable* (`n >= 0` integer, `x > -1` real):

    F(n,x) = c_n (3/4)^x (-x)_{2n} Gamma(x+n+1)/Gamma(x+n+3/2),
    c_n    = 9^-n Gamma(3n+3/2) n!^2 / (4 (2n)! (3n)! Gamma(n+5/6) Gamma(n+7/6)),
    R(n,x) = (8n^2 + 6n - (2n+3)x - x^2) / ((3n+1)(3n+2)),      G = R F.

Origin: Ebisu's `2F1(2a, 3a-1/2; 3a; 3/4)` (loc. cit., family (1,3,3-1)) as a seed, `a -> 1/2+n`, `k -> x - 2n`.
For integer `x = k` the pair is trivial (`G(n,0) = 0`); the identity is its first-order term in `x`.

*Proof.*
1. `F(n+1,x) - F(n,x) = G(n,x+1) - G(n,x)` for all real `x > -1`: (WZ) holds as a rational identity for the
   ratios `sigma = (3/4)(x+1)(x+n+1)/((x+1-2n)(x+n+3/2))`, `rho` (script), so the equation holds wherever
   `F(n,x) != 0`, and both sides are analytic in `x`.
2. Differentiate in `x` and put `x = k`: with `'` = `d/dx`,
   `F'(n+1,k) - F'(n,k) = G'(n,k+1) - G'(n,k)` for integers `n,k >= 0`. `F'(n,x)` and `G'(n,x)` are
   `O(x^{2n+3} (3/4)^x)`, so `D_n := sum_{k>=0} F'(n,k)` converges absolutely and `G'(n,K) -> 0`. Summing over `k`:
   `D_{n+1} - D_n = -G'(n,0)`, hence `sum_{n<N} G'(n,0) = D_0 - D_N`.
3. `G'(n,0) = -(3/(4 pi)) t(n+1)`, `t` the summand: for `n >= 1`, `(-x)_{2n}` has a simple zero at 0 with
   derivative `-(2n-1)!`, so `G'(n,0) = -R(n,0) c_n (2n-1)! n!/Gamma(n+3/2)`; `G'(0,0) = -(3/2) c_0/Gamma(3/2)
   = -9/(8 pi)` (using `Gamma(5/6)Gamma(7/6) = pi/3`) and `G'(1,0) = -7/(160 pi) = -(3/(4 pi)) t(2)`; for `n >= 1` the ratio
   `G'(n+1,0)/G'(n,0) = t(n+2)/t(n+1)` is checked exactly.
4. `D_0 = -3K/(2 pi)`. By the Beta integral, `F(0,x) = (c_0/sqrt(pi)) int_0^1 (3t/4)^x (1-t)^{-1/2} dt`,
   `c_0 = 3/(8 sqrt(pi))`; differentiating under the integral and summing the geometric series (dominated
   convergence), `D_0 = (3/(8 pi)) I`, `I = int_0^1 log(3t/4) / ((1-3t/4) sqrt(1-t)) dt`.
   With `t = 1-u^2`, `u = tan(theta)/sqrt3`:
   `I = (8/sqrt3) int_0^{pi/3} [log cos(theta+pi/6) + log cos(theta-pi/6) - 2 log cos(theta)] d theta`
   (because `3 cos^2 - sin^2 = 4 cos(theta+pi/6) cos(theta-pi/6)`). By `int_0^alpha log(2 cos phi) d phi = (1/2) Cl_2(pi - 2 alpha)`:
   `I = (8/sqrt3) [ (1/2) Cl_2(2pi/3) - Cl_2(pi/3) ] = -(8/sqrt3) Cl_2(2pi/3) = -4K`,
   using `Cl_2(pi/3) = (3/2) Cl_2(2pi/3)` (duplication) and `Cl_2(2pi/3) = (sqrt3/2) K`.
5. `D_N -> 0`. Again by the Beta integral, `F(N,x) = (c_N/sqrt(pi)) (-x)_{2N} int_0^1 (3t/4)^x t^N (1-t)^{-1/2} dt`.
   Since `(-x)_{2N} y^x = y^{2N} (d/dy)^{2N} y^x`, we get `sum_{k>=0} d/dx[(-x)_{2N} y^x]_{x=k} = Psi_N(y) := y^{2N} (d/dy)^{2N} [log y/(1-y)]`
   for `0 < y < 1`, and (dominated convergence; the terms are bounded by `C_N (1+k)^{2N} y^k (1+|log y|)`, whose sum over `k`
   is `<= C'_N (1+|log y|)` for `0 < y <= 3/4`, integrable against `t^N (1-t)^{-1/2}`) `D_N = (c_N/sqrt(pi)) int_0^1 t^N (1-t)^{-1/2} Psi_N(3t/4) dt`.
   By Leibniz' rule and `log(1+u) - sum_{i<=M} (-1)^{i-1} u^i/i = (-1)^M int_0^u s^M/(1+s) ds` with `u = (1-y)/y`:
   `Psi_N(y) = -(2N)! y^{2N} (1-y)^{-2N-1} int_0^{(1-y)/y} s^{2N}/(1+s) ds`, so `|Psi_N(y)| <= (2N)!/((2N+1) y)` and
   `|D_N| <= c_N (2N)! 4 B(N,1/2) / (3 sqrt(pi) (2N+1)) = O(9^{-N})`.
6. Therefore `sum_{n>=0} G'(n,0) = D_0 = -3K/(2 pi)`, i.e. `sum_{k>=1} t(k) = 2K`. ∎

*Check.* `proofs/verify_T3_sun_3k.py`: (WZ) and step 3 exactly; steps 2, 4, 5 numerically (D_N by direct summation
against the integral formula and the bound; `I = -4K` to 30 digits); final identity to 110 digits.

---

## Theorem 4.  `sum_{k>=0} (360k^3+612k^2+230k+15) C(3k,k) C(6k,3k)^2 / ((3k+1)(3k+2) 2^(15k)) = 32 sqrt2/pi`

*Source.* Sun, arXiv:2307.03086v10, Conjecture 5.25(i), eq. (5.20); still "conjectured" in arXiv:2603.29973v3 (Apr 2026).

*Pair.* `h(0,0) = 1`,

    sigma(n,k) = (4k+6n+1)(4k+6n+3) / (16 (k+n+1)(k+3n+1)),
    rho(n,k)   = 3(6n+1)(6n+5)(4k+6n+1)(4k+6n+3)(4k+6n+5) / (1024 (n+1)(k+n+1)(k+3n+1)(k+3n+2)(k+3n+3)),
    Q(n,k)     = (64k^3+(160n+144)k^2+(240n^2+432n+92)k+360n^3+612n^2+230n+15) / (64 (k+3n+1)(k+3n+2)),

and the WZ pair is `F(n,k) = n h(n,k)`, `G(n,k) = Q(n,k) h(n,k)`:
`(n+1) rho - n = Q(n,k+1) sigma - Q(n,k)`. (`n h(n,k)` is a constant times `Seed1(1/2+n, 1+2n, k+n)`, Au's Seed1.)
`h(n,0) = C(3n,n) C(6n,3n)^2/2^(15n)`, `h(0,k) = (1/4)_k (3/4)_k/k!^2`, `Q(n,0) h(n,0) = t(n)/64`.

*Proof.* `h > 0`. From `1 - sigma(n,k) >= (3/4)(n+1)/(k+3n+1)`:
`h(n,k) <= h(n,0) ((3n+1)/(k+3n+1))^{3(n+1)/4}` (as `prod_{j<k} (1 - c/(j+m)) <= (m/(k+m))^c`). Hence for `n >= 2`,
`sum_k h(n,k) <= 7 h(n,0)` (since `sum_{k>=0} (m/(k+m))^c <= 1 + m/(c-1) = 1 + 4(3n+1)/(3n-1) <= 33/5`), and for `n >= 1`,
`Q(n,K) h(n,K) <= (K+3n+3) h(n,K) -> 0` (`Q <= k+3n+3`, exponent `>= 3/2`). For `n = 0`:
`Q(0,K) h(0,K) ~ K (1/4)_K (3/4)_K/K!^2 -> 1/(Gamma(1/4) Gamma(3/4)) = sin(pi/4)/pi = 1/(sqrt2 pi)`.
The telescoping lemma (with `u = n h`, `v = Q h`, `u(0,k) = 0`) and `K -> oo` give, for `N >= 1`,
`sum_{n<N} Q(n,0) h(n,0) = 1/(sqrt2 pi) - N sum_k h(N,k)`, and for `N >= 2`, `N sum_k h(N,k) <= 7N (27/32)^N -> 0`
(`rho(n,0) <= 27/32`). So `sum_n t(n)/64 = 1/(sqrt2 pi)`. ∎

*Check.* `proofs/verify_T4_sun_1pi_sqrt2.py` (exact: compatibility, (WZ), boundary data, the three inequalities;
numerics: the limit and the final identity to 105 digits with 4200 terms).

---

## Theorem 5.  `sum_{k>=1} (77k^3-86k^2+29k-3)(-32)^k / (k^3 (3k-1)(3k-2) C(3k,k) C(4k,k)^2) = -6G`

**Not a new result.** With `k = n+1` this is equation (7.10) of W. Chu, *Infinite series formulae related to Gauss and
Bailey 2F1(1/2)-sums*, Acta Math. Sci. 40B (2020), no. 2, 293–315 (Theorem 7.2 with `a = c = 1`, Example 7.3):
`2G = sum_{n>=0} (-1/32)^n (1)_n^2 (1)_{3n} / ((3/2)_n (5/2)_{2n}^2) (17+88n+145n^2+77n^3)/9`, proved there by the modified
Abel lemma. Chu's summand is `-(1/3) t(n+1)`, exactly our `gamma(n,0)` (independently confirmed in
`audit2/a6_T5_is_Chu2020.py`). We missed it because we searched for the polynomial `77k^3-86k^2+29k-3` as a string.
The WZ proof below is kept as a second proof.

*Source.* Sun, arXiv:2307.03086v10, Conjecture 4.33(i) (first identity, dated 2025-02-28), stated there as a conjecture.

*Pair.* `psi(0,0) = 1`,

    sigma(n,k) = (2k+n+1)(2k+n+2) / ((2k+2n+3)(2k+4n+3)),
    rho(n,k)   = -(n+1)(3n+1)(3n+2)(3n+3)(2k+n+1) / ((4n+1)(4n+3)(2k+2n+3)(2k+4n+3)(2k+4n+5)),
    R(n,k)     = ((16n+12)k^2+(66n^2+86n+28)k+77n^3+145n^2+88n+17) / ((4n+1)(4n+3)(2k+4n+3)),   gamma = R psi.

Origin: `Seed1(-nu, 1+nu, k+1/2+nu)` (Au's Seed1) vanishes at integers `nu = n` through a factor `1/Gamma(-nu)`;
`psi` is its first-order coefficient (same ratios, same certificate), normalised by `psi(0,0) = 1`.

*Proof.* (WZ) holds. `psi(0,k) = 4^k/((2k+1)^2 C(2k,k))`, and
`sum_k psi(0,k) = int_0^1 arcsin(x)/(x sqrt(1-x^2)) dx = int_0^{pi/2} theta/sin(theta) d theta = 2G`
(from `arcsin(x)/sqrt(1-x^2) = sum_k 4^k x^{2k+1}/((2k+1) C(2k,k))`, positive terms).
`1 - sigma(n,k) >= (n+3/4)/(k+2n+3/2)` gives `|psi(n,k)| <= |psi(n,0)| ((2n+3/2)/(k+2n+3/2))^{n+3/4}`; with
`|rho(n,0)| <= 1/8` and `0 <= R <= 2k+2n+2` this yields `sum_k |psi(N,k)| <= 6 * 8^{-N}` for `N >= 1` and `gamma(n,K) -> 0`
(for `n = 0`: `gamma(0,K) = O(K^{-1/2})`). Telescoping: `sum_n gamma(n,0) = 2G`, and `gamma(n,0) = -(1/3) t(n+1)`. ∎

*Check.* `proofs/verify_T5_sun_catalan.py`.

---

## Theorem 6.  `sum_{k>=0} (3k+2)^2 16^k H_k / ((2k+1)(6k+1)(6k+5) C(6k,3k)) = (5/2) K - pi log(432)/(6 sqrt3)`

*Source.* Sun, arXiv:2506.01870v6, Conjecture 3.2 (second identity).

*Two-parameter extension of the pair of Theorem 1* (`(a,b)` real near `(0,0)`, `u = k+3n+2a+b`):

    F(n,k;a,b) = (-1)^{n+k} 3^{-(k+n+b)} (3/4)^{a+n+1/2} Gamma(u+1) Gamma(k+b-a+1/2) Gamma(2n+2a+3/2) Gamma(n+a+1)^2
                 / ( Gamma(2n+2a+1) Gamma(1/2-a-n) Gamma(u+3/2) Gamma(k+n+b+1) Gamma(n+a+3/4) Gamma(n+a+5/4) ),
    R(n,k;a,b) = 3 (u+2)^2 / (4 (u+2)^2 - 1),      G = R F.

*Proof.*
(i) (WZ) holds identically in `a, b` (script). (ii) Let `B = [-1/10,1/10]^2`. For `(a,b) in B`, `n, k >= 0`:
`sigma = F(n,k+1)/F(n,k) = -(1/3)(u+1)(k+b-a+1/2)/((u+3/2)(k+n+b+1))` satisfies `-1/3 <= sigma < 0`
(both quotients `(u+1)/(u+3/2)` and `(k+b-a+1/2)/(k+n+b+1)` lie in `(0,1]`), `0 < F(n+1,0)/F(n,0) <= 1/2`, and
`0 < R <= 1` (as `u+2 >= 1`); these six inequalities are checked exactly on the whole box in
`proofs/verify_R4_uniform_bounds.py`. Hence `|F(n,k)| <= |F(0,0;a,b)| 2^{-n} 3^{-k}`, `sum_k |F(N,k)| -> 0`,
`G(n,K) -> 0`, and the telescoping lemma gives `Phi(a,b) := sum_n G(n,0;a,b) = sum_k F(0,k;a,b)` for every
`(a,b) in B`. By the differentiation lemma of the section "Theorems 8, 9, 10" below (the logarithmic derivatives of
`F(n,0;a,b)` in `a` and `b` are sums of at most 12 digamma values at arguments in `[3/10, 3n+4]`, plus
`psi(1/2-a-n) = psi(1/2+a+n) - pi tan(pi a)`, hence `O(log(n+2))` uniformly on `B`; `R` and its partial derivatives
are bounded on `B`), `Phi` is `C^1` on `B` and its partial derivatives are obtained termwise from the left-hand
series. The right-hand series is differentiated in (iii), (iv) only along the two axes, where it is explicit.
(iii) `sum_k F(0,k;a,0) = 2^{-1/2}` for all small `a`: it equals a Gamma factor times `2F1(2A,1-A;2A+1/2;-1/3)`,
`A = a+1/2`, and Ebisu's formula quoted in Theorem 1 with the duplication formula gives the constant.
So `d/da Phi(0,0) = 0`.
(iv) `F(0,k;0,b) = (sqrt6/(2 pi)) (-1)^k 3^{-k-b}/(k+b+1/2)`, so
`d/db Phi(0,0) = -(sqrt6/(2 pi)) (2 log3 S_1 + 4 S_2)`, `S_1 = sum (-1/3)^k/(2k+1) = pi/(2 sqrt3)`,
`S_2 = sum (-1/3)^k/(2k+1)^2 = sqrt3 Ti_2(1/sqrt3)`, and `Ti_2(tan t) = t log tan t + (1/2) Cl_2(2t) + (1/2) Cl_2(pi-2t)`
(integrate `int_0^x arctan(s)/s ds` by parts; Lewin, *Polylogarithms and associated functions*, Ch. 2) gives
`S_2 = -(sqrt3 pi/12) log 3 + (15/8) K`. Hence `d/db Phi(0,0) = -(15/4)(sqrt6/pi) K`.
(v) `R(n,0;a,b)` depends on `2a+b` only, so `(d/da - 2 d/db)` kills it; and
`(d/da - 2 d/db) log F(n,0;a,b)|_0 = 4 psi(n+1) - 2 psi(2n+1) + psi(n+1/2) - 3 psi(1/2) + log 27 = 3 H_n + log 432`
(the terms `2 psi(2n+3/2) - psi(n+3/4) - psi(n+5/4) = 2 log 2` by duplication, `psi(1/2-n) = psi(1/2+n)`).
Since `G(n,0;0,0) = (sqrt6/pi) t(n)`: `sum_n t(n) (3 H_n + log 432) = (15/2) K`, and Theorem 1 gives the claim. ∎

*Check.* `proofs/verify_T6_sun_16k_harmonic.py` ((WZ) with parameters exactly; (ii)–(v) numerically to 30–40
digits; final identity to 110 digits). The derivative `d/db` alone gives one more (new) identity involving `H_{6n}, H_{3n}, H_n`.

---

## Theorem 7.  `sum_{k>=1} C(2k,k) 16^k (1656k^4-2604k^3+1158k^2-229k+15) / (k^2 (2k-1)^2 (6k-1)(6k-5) C(6k,3k) C(3k,k)) = 0`

*Source.* Sun, MathOverflow 506695 (6 Jan 2026); one answer (partial fractions, no proof).

*Pair supported on the triangle `0 <= k < n`:*

    Phi(n,k) = (-1)^{k+1} (4k+4n+3) 2^{1-6n} k! (4n)! (2n)! / ( pi (n-k-1)! (2k+2n+1)(2k+2n+2) Gamma(k+2n+3/2) Gamma(k+3n+3/2) )   (k < n),
    Phi(n,k) = 0 (k >= n);   Gam(n,k) = R(n,k) Phi(n,k) (k < n),  Gam(n,n) = -Phi(n+1,n),  Gam(n,k) = 0 (k > n),
    R(n,k) = P(n,k) / (3(6n+1)(6n+5)(k-n)(2k+4n+3)(2k+6n+3)(2k+6n+5)(4k+4n+3)),

`P` the polynomial of degree 6 in `k` printed in `proofs/verify_T7_sun_zero.py`. These are the values at integer
`nu = n` of `F(nu,k) = Seed3(1/2-2nu, -2nu, 1/2-3nu, k+1/2+2nu)` (Au's Seed3) and of `R F`; `F` contains
`Gamma(k-nu+1)/Gamma(-2nu)`, which is `0` for `k >= n` and finite for `k < n`.

*Proof.* The WZ equation `Phi(n+1,k) - Phi(n,k) = Gam(n,k+1) - Gam(n,k)` holds for all `n,k >= 0`: for `k <= n-2` it is
(WZ); for `k = n-1` it is a second rational identity in `n`; for `k = n` it is the definition of `Gam(n,n)`; for `k > n`
all terms vanish (script: both identities, and all cases exactly on a grid). The telescoping lemma with `K = N+1`:
`sum_{n<N} Gam(n,0) = -sum_{k<N} Phi(N,k)`. On the support (`0 <= k <= n-2`) `|Phi(n,k+1)/Phi(n,k)| <= 1/24`, and
`0 < Phi(n+1,0)/Phi(n,0) <= 1/6` for `n >= 1`, so `|sum_{k<N} Phi(N,k)| <= (24/23) |Phi(1,0)| 6^{1-N} -> 0`.
(Both constants are verified in the script by coefficient positivity after the substitutions that parametrise the
range: `n = k+2+v`, `v >= 0`, for the first — `1/24 - |sigma|` then has nonnegative coefficients in `(k,v)` — and
`n = m+1`, `m >= 0`, for the second.) `Gam(n,0) = -t(n+1)/(6 pi^2)`. ∎

Remark: the proof gives the exact finite identity `sum_{m<=N} t(m) = 6 pi^2 sum_{k<N} Phi(N,k)` (a rational number).

*Check.* `proofs/verify_T7_sun_zero.py`.


---

## Theorems 8, 9, 10: harmonic-number companions by differentiating in a continuous parameter

**Principle.** The WZ equation is a rational identity, so the pairs above remain WZ pairs when the second
variable is real (`x = k + c`) or when the seed parameters are moved. The boundary analysis is then redone for
every value of the parameter in `I = [-1/10, 1/10]` (the inequalities, with the parameter, are verified exactly in
`proofs/verify_R4_uniform_bounds.py`; they are listed in each proof), which gives an identity
`sum_n G(n; c) = Phi(c)` for every `c in I` with `Phi` explicit. It is differentiated by the following standard fact.

**Differentiation lemma** (e.g. Rudin, *Principles of Mathematical Analysis*, Thm. 7.17). Let `u_n in C^1(I)`,
`sum_n u_n(c_0)` convergent for one `c_0`, and `sum_n sup_I |u_n'| < oo`. Then `sum_n u_n` converges uniformly on
`I` to a `C^1` function and `(sum_n u_n)' = sum_n u_n'`.

It is applied to `u_n = G(n; c) = R F` (resp. `Q h`), with `u_n' = u_n (d/dc log F + d/dc log R)`. In each case
below (a) `|F(n; c)| <= C r^n` on `I` with `r < 1`, from the uniform bound on `rho`; (b) `d/dc log F(n; c)` is a
finite sum (at most 12 terms, with fixed integer coefficients) of values `psi(t)` at arguments
`3/10 <= t <= 3n + 4`, plus bounded terms, and `|psi(t)| <= 4 + log(1+t)` for `t >= 3/10`, so it is
`O(log(n+2))` uniformly on `I`; (c) `R` and `dR/dc` (resp. `Q`, `dQ/dc`) are rational functions whose denominators
are products of linear forms bounded below by `1/2` on `I` and whose numerators are polynomials, hence
`O((n+1)^A)` uniformly on `I` for a fixed `A`. So `sup_I |u_n'| = O((n+1)^A log(n+2) r^n)`, which is summable.
`d/dc log F` produces digamma values, i.e. harmonic numbers. Theorem 6 is the same principle with two parameters.

### Theorem 8.  `sum_{k>=0} C(3k,k) C(6k,3k)^2 [P(k)(H_{2k}-H_k) - 180k^2+36k+23] / ((3k+1)(3k+2) 2^(15k)) = 48 sqrt2 log2/pi`, `P(k) = 360k^3+612k^2+230k+15`

*Source.* Sun, arXiv:2307.03086v10, Conjecture 5.25(ii), first identity.

*Proof.* For `|s| <= 1/10` let `h_s(n,k) = 2^{-2k-5n+s-3/2} Gamma(3n+s+1/2) Gamma(2k+3n-s+1/2) / (Gamma(n+s+1/2) Gamma(k+n-s+1) n! Gamma(k+3n+1))`
(`n h_s(n,k) = Seed1(1/2+s+n, 1+s+2n, k-s+n)`), with certificate `Q_s(n,k) = N_s(n,k)/(64(2n+2s+1)(k+3n+1)(k+3n+2))`
(`N_s` in `proofs/verify_T8_sun_1pi_sqrt2_harmonic.py`); `(n+1) rho_s - n = Q_s(n,k+1) sigma_s - Q_s(n,k)` identically in `s`.
For every `|s| <= 1/10` and `n, k >= 0` (exact, `proofs/verify_R4_uniform_bounds.py`): `h_s > 0` (all Gamma
arguments are `>= 2/5`), `0 < sigma_s`, `1 - sigma_s(n,k) >= (7/10)(n+1)/(k+3n+1)`, `0 < rho_s(n,0) <= 9/10`,
`0 <= Q_s(n,k) <= 2(k+3n+3)`, `Q_s(0,k)/k -> 1`. Then the proof of Theorem 4 goes through with these constants:
with `c = 7(n+1)/10 <= m = 3n+1`, `h_s(n,k) <= h_s(n,0) (m/(k+m))^c`; for `n >= 1` (`c >= 7/5`)
`sum_k h_s(n,k) <= (1 + m/(c-1)) h_s(n,0) <= 11 h_s(n,0)` and `Q_s(n,K) h_s(n,K) <= 2(K+3n+3) h_s(n,0) (m/(K+m))^{7/5} -> 0`;
`h_s(N,0) <= h_s(0,0) (9/10)^N`, so `N sum_k h_s(N,k) -> 0`. For `n = 0`, by the duplication formula
`Gamma(2K-s+1/2) = 2^{2K-s-1/2} Gamma(K-s/2+1/4) Gamma(K-s/2+3/4)/sqrt(pi)` and `Gamma(K+p)/Gamma(K+q) ~ K^{p-q}`,
`Q_s(0,K) h_s(0,K) ~ K 2^{-2K+s-3/2} Gamma(2K-s+1/2)/(Gamma(K-s+1) K!) -> 1/(4 sqrt(pi))`. Hence
`sum_n Q_s(n,0) h_s(n,0) = 1/(4 sqrt(pi))` for every `|s| <= 1/10`, independent of `s`.
The differentiation lemma applies (`r = 9/10`; `d/ds log h_s(n,0) = log 2 + psi(3n+s+1/2) - psi(3n-s+1/2) - psi(n+s+1/2) + psi(n-s+1)`;
the denominator of `Q_s(n,0)` is `64(2n+2s+1)(3n+1)(3n+2)`). Differentiate at `s = 0`:
`d/ds log h_s(n,0) = log 2 - psi(n+1/2) + psi(n+1) = 3 log 2 - 2(H_{2n}-H_n)` and
`d/ds log Q_s(n,0) = (360n^2-72n-46)/P(n)` (exact). Hence
`sum_n t(n) [3 log 2 - 2(H_{2n}-H_n) + (360n^2-72n-46)/P(n)] = 0`, and Theorem 4 gives the claim. ∎

### Theorem 9.  `sum_{k>=1} C(4k,2k) [(8k-3)(5H_{2k-1} - 4H_{k-1}) - 6] / (k (4k-1) 9^k C(2k,k)^2) = 3K/2`

*Source.* Sun, arXiv:2307.03086v10, Conjecture 5.9(i), second identity.

*Proof.* Use the pair of Theorem 2 with real second variable: `f(n,x)` (Gamma closed form of Theorem 2, `x > -1/2`),
`g = R f`; `f(n,x) > 0` for `x >= -1/10` (the only sign-carrying factor, `1/Gamma(-2n-1/2)`, does not depend on
`x`), and `f(0,x) = 1/((3x+1)(3x+2))` identically in `x`. For all real `x >= -1/10` (exact,
`proofs/verify_R4_uniform_bounds.py`): `0 < sigma(n,x) <= (2x+1)/(2x+5)` for `n >= 1`, `0 < R(n,x) <= x+2`, and
`0 < rho(n,c) <= 1/5` for `|c| <= 1/10`, `n >= 0`. Fix `|c| <= 1/10` and apply the telescoping lemma to
`u(n,k) = f(n,k+c)`, `v(n,k) = g(n,k+c)`: for `n >= 1`,
`f(n,k+c) <= f(n,c) (2c+1)(2c+3)/((2k+2c+1)(2k+2c+3))` (the product of the bounds telescopes), so
`sum_k f(n,k+c) <= (2c+3) f(n,c)/2 <= (8/5) f(0,c) 5^{-n} -> 0` and
`g(n,K+c) <= (K+c+2) f(n,K+c) -> 0`; `g(0,x) = (6x+5)/(9(2x+1)(3x+1)) -> 0`. Hence, for every `|c| <= 1/10`,
`Phi(c) := sum_n g(n,c) = sum_k f(0,k+c) = sum_k 1/((3k+3c+1)(3k+3c+2)) = (psi(c+2/3) - psi(c+1/3))/3`.
The differentiation lemma applies (`r = 1/5`; `d/dx log f(n,x) = 2/(2x+n+1) + psi(x+1/2) + 3 psi(3x+n+1) - psi(x+n+3/2) - 3 psi(3x+2n+3)`;
the denominator of `R(n,x)` is `9(2n+1)(2x+n+1)`).
So `Phi'(0) = (psi'(2/3) - psi'(1/3))/3 = -3K`. On the other side
`d/dx log f(n,x)|_0 = 2/(n+1) + psi(1/2) + 3 psi(n+1) - psi(n+3/2) - 3 psi(2n+3) = 4H_n - 5H_{2n+1} + 1/(2n+2)` and
`d/dx log R(n,x)|_0 = 6/(8n+5) - 1/(2n+2)`; with `g(n,0) = 2 t(n+1)`:
`2 sum_m t(m) [4H_{m-1} - 5H_{2m-1} + 6/(8m-3)] = -3K`. ∎

### Theorem 10.  `sum_{k>=0} C(3k,k) C(6k,3k)^2 [P(k)(4H_{6k} - 3H_{3k} - H_k) + f(k)] / ((3k+1)(3k+2) 2^(15k)) = 192 sqrt2 log2/pi`, `f(k) = (1296k^3+1980k^2+912k+139)/((3k+1)(3k+2))`

*Source.* Sun, arXiv:2307.03086v10, Conjecture 5.25(ii), second identity.

*Proof.* Use the pair of Theorem 4 with real second variable: `h(n,x) = 2^{-2x-5n} (1/2)_{3n} Gamma(2x+3n+1/2) / (sqrt(pi) (1/2)_n n! Gamma(x+n+1) Gamma(x+3n+1))`
(positive for `x >= -1/10`; ratios `sigma(n,x)`, `rho(n,x)` of Theorem 4 with `k` replaced by `x`), on the rows
`x = k + b`, `|b| <= 1/10`. For all real `x >= -1/10`, `n >= 0` (exact, `proofs/verify_R4_uniform_bounds.py`):
`0 < sigma`, `1 - sigma(n,x) >= (3/4)(n+1)/(x+3n+1)`, `(3/4)(n+1) <= x+3n+1`, `0 <= Q(n,x) <= x+3n+3`,
`Q(0,x)/x -> 1`, and `0 < rho(n,b) <= 9/10` for `|b| <= 1/10`. So, with `c = 3(n+1)/4` and `m = 3n+1+b`,
`h(n,k+b) <= h(n,b) (m/(k+m))^c`; for `n >= 1`, `sum_k h(n,k+b) <= (1 + m/(c-1)) h(n,b) <= 10 h(n,b)` and
`Q(n,K+b) h(n,K+b) -> 0`; `h(N,b) <= h(0,b) (9/10)^N`; and
`Q(0,x) h(0,x) ~ x Gamma(x+1/4) Gamma(x+3/4)/(Gamma(1/4) Gamma(3/4) Gamma(x+1)^2) -> 1/(sqrt2 pi)` as `x -> oo`.
The telescoping lemma on the rows `x = k+b` (with `u = n h`, `v = Q h`) therefore gives
`sum_n Q(n,b) h(n,b) = 1/(sqrt2 pi)` for every `|b| <= 1/10`. The differentiation lemma applies (`r = 9/10`; the
denominator of `Q(n,x)` is `64(x+3n+1)(x+3n+2)`). Differentiating at `b = 0`:
`d/dx log h(n,x)|_0 = -2 log 2 + 2 psi(3n+1/2) - psi(n+1) - psi(3n+1) = (4H_{6n} - 3H_{3n} - H_n) - 6 log 2`,
`d/dx log Q(n,x)|_0 = f(n)/P(n)` (exact). Hence `sum_n t(n) [(4H_{6n}-3H_{3n}-H_n) + f(n)/P(n)] = 6 log 2 * 32 sqrt2/pi`. ∎

*Check.* `proofs/verify_T8_sun_1pi_sqrt2_harmonic.py`, `proofs/verify_T9_T10_harmonic_kshift.py` (WZ with the parameter and the
derivatives of the certificates exactly; constancy of `Phi` and the digamma identities numerically to 30–40 digits;
final identities to 104–110 digits). Together with Theorem 4 this settles all of Sun's Conjecture 5.25(i),(ii)
(series part; the supercongruences in that conjecture are not addressed).

---

## Theorem 11.  `sum_{k>=0} (1344k^3+944k^2+156k+9) C(2k,k)^2 C(4k,2k)^3 / ((2k+1) 2^(22k)) = 64 sqrt2/pi^2`

*Source.* Sun, arXiv:2307.03086v10, Conjecture 5.33(i), eq. (5.31) (dated 2025-06-26); still "conjectured" in
arXiv:2603.29973v3 (Apr 2026). A Ramanujan-type series for `1/pi^2` (rate 1/64). Write `t(n)` for the summand.

*Status of this proof.* Part A (algebra) is checked exactly by `proofs/verify_T11_sun_1pi2_sqrt2.py`. Part B is a
classical Laplace-method argument written out in full below; its algebraic ingredients are checked exactly and its
conclusions numerically (N up to 20000), but Part B is **not machine-checked** and is the longest analytic step in this file.

**Part A: a WZ pair on a triangle.** For integers `n >= 1`, `k >= 0` let

    a(n,k) = C_n (2k+1) Gamma(4k+2n+1/2) / ( (n-k-1)!^2 (n+k)!^2 Gamma(4k-2n+9/2) )   (k < n),     a(n,k) = 0  (k >= n),
    C_n    = 2^(6-8n) Gamma(n+3/4)^2 Gamma(n+1/4) / ( sqrt2 pi Gamma(n-1/4) Gamma(2n-1) ),
    gam(n,k) = R(n,k) a(n,k) (k < n),   gam(n,n) = -a(n+1,n),   gam(n,k) = 0 (k > n),
    R(n,k) = NUM(n,k) / ( 131072 n (2n-1)(2n+1)(4n-1) (k-n)^2 (2k+1) ),

`NUM` the polynomial of degree 6 in `k` printed in the script. (`a`, `gam` are, up to the constant `e^{-i pi/4}`, the
coefficients of `eps^-2` of `F(nu,k) = Seed9(1-2nu, 1/4-nu, 3/4-nu, k+nu)` and of its Gosper mate at `nu = n + eps`; Au's
Seed9 is Dougall's 7F6 with `b,c -> b/2, b/2+1/2` and `d,e -> d/2, d/2+1/2`. This remark is not used in the proof.)
`a(n,k) > 0` for `0 <= k < n` (if `4k-2n+9/2 < 0` it equals `1/2 - j` with `j` even).
With `sigma(n,k) = a(n,k+1)/a(n,k)`, `rho(n,k) = a(n+1,k)/a(n,k)`:

    sigma = (2k+3)(n-k-1)^2 (4k+2n+1/2)_4 / ( (2k+1)(n+k+1)^2 (4k-2n+9/2)_4 ),
    rho   = (n+3/4)^2 (n+1/4) (4k+2n+1/2)(4k+2n+3/2)(4k-2n+5/2)(4k-2n+7/2) / ( 512 n (n-1/4)(2n-1) (n-k)^2 (n+k+1)^2 ).

*WZ equation* `a(n+1,k) - a(n,k) = gam(n,k+1) - gam(n,k)` for all `n >= 1`, `k >= 0`: for `k <= n-2` it is the rational
identity (WZ); for `k = n-1` it is a second rational identity in `n`; for `k = n` it is the definition of `gam(n,n)`; for
`k > n` every term is 0. Moreover `gam(n,0) = -t(n)/256` for `n >= 1` and `a(1,0) = t(0)/256`. Telescoping over
`1 <= n < N`, `0 <= k <= N`:

    sum_{n<N} t(n) = 256 T_N,        T_N := sum_{k=0}^{N-1} a(N,k).                                   (A)

**Part B: `T_N -> sqrt2/(4 pi^2)`.** Fix `N` large, write `a(k) = a(N,k)`, `r(k) = a(k+1)/a(k) = sigma(N,k)`,
`x* = sqrt3/2`, `k0 = floor(x* N)`, `theta = x* N - k0`, `J = floor(N^{3/5})`.

(B1) *Shape.* Let `I = {k : N/2 - 1 <= k <= N-1}`. On `I`, `4k-2N+9/2 >= 1/2`, and
`r(k) = (1 + 2/(2k+1)) ((N-k-1)/(N+k+1))^2 prod_{i=0}^{3} (1 + (4N-4)/(4k-2N+9/2+i))` is a product of nonnegative decreasing
functions of `k`: `a` is log-concave on `I`. For `0 <= k <= N/2 - 2` all four numbers `4k-2N+9/2+i` are negative, so
`r(k) = ((2k+3)/(2k+1)) ((N-k-1)/(N+k+1))^2 prod_i (4k+2N+1/2+i)/(2N-4k-9/2-i)`; each factor of the product is
`>= (2N+4k)/(2N-4k)` (larger numerator, smaller positive denominator), hence
`r(k) >= ((N-k-1)(N+2k)^2 / ((N+k+1)(N-2k)^2))^2 >= 1`, because
`(N-k-1)(N+2k)^2 - (N+k+1)(N-2k)^2 = 6N^2 k - 8k^3 - 2N^2 - 8k^2 >= 4N^2 (k-1) >= 0` for `1 <= k <= N/2`
(use `8k^3 <= 2N^2 k`, `8k^2 <= 2N^2`); and `r(0) >= 3((N-1)/(N+1))^2 >= 1` for `N >= 4`.
For odd `N` the one remaining index `k = (N-3)/2` has the four numbers `-3/2,-1/2,1/2,3/2` (product `9/16`), so
`r(k) >= ((N+1)/(3N-1))^2 (16/9)(4N-11/2)^4 >= (16/81)(4N-11/2)^4 >= 1`.
Hence `a(k) <= a(min I)` for all `k < min I`.

(B2) *Ratio in the window.* For `k/N = x` in `[0.8, 0.93]`,
`log r(k) = L(x) + O(1/N)`, `L(x) = log[(1-x)^2 (2x+1)^4 / ((1+x)^2 (2x-1)^4)]`, uniformly (each factor of `r` differs from its
limit by a factor `1 + O(1/N)`). `L(x*) = 0` (`(1-x)(2x+1)^2 = (1+x)(2x-1)^2` iff `8x^3 = 6x`), `L'(x*) = -24`, `L''` bounded.
So for `|j| <= J+1`: `log r(k0+j) = -24 j/N + O(1/N + j^2/N^2)`, and summing (over `|j| <= J`),
`log a(k0+j) = log a(k0) - 12 j^2/N + O(|j|/N + |j|^3/N^2) = log a(k0) - 12 j^2/N + O(N^{-1/5})`.

(B3) *Centre.* `sum_{|j|<=J} a(k0+j) = a(k0) (1+O(N^{-1/5})) sum_{|j|<=J} e^{-12 j^2/N} = a(k0) sqrt(pi N/12) (1+o(1))`
(for the unimodal function `e^{-12 t^2/N}` the sum over `Z` differs from the integral `sqrt(pi N/12)` by at most 1, and the
part `|j| > J` is `O(N e^{-12 N^{1/5}})`).

(B4) *Tails.* By log-concavity on `I` and (B2): for `k >= k0+J` in `I`, `a(k) <= a(k0+J) r(k0+J)^{k-k0-J}` with
`r(k0+J) <= e^{-12J/N}`; for `k <= k0-J` in `I`, `a(k) <= a(k0-J) r(k0-J-1)^{-(k0-J-k)}` with `r(k0-J-1) >= e^{12J/N}`.
Both tails are therefore `<= a(k0 +- J) N/(6J) = a(k0) O(N^{2/5} e^{-11 N^{1/5}})`. In particular
`a(min I) <= a(k0) e^{-4 N^{3/5}}`, and by (B1) `sum_{k < min I} a(k) <= N a(min I)`. So everything outside the window is
`o(a(k0) sqrt N)`.

(B5) *Size at the centre.* By Stirling's formula `log Gamma(z) = (z-1/2) log z - z + (1/2) log 2pi + O(1/z)` applied to the
eight Gamma functions in `a(N,k)` (all arguments are `>= 0.13 N` for `k/N` near `x*`), for `k = xN`:
`log a(N,k) = N phi(x) - (1/2) log N + log psi0(x) + O(1/N)`, where

    phi(x)  = -6 log 2 + 2(2x+1) log(2x+1) - 2(1-x) log(1-x) - 2(1+x) log(1+x) - 2(2x-1) log(2x-1),     phi' = L,
    psi0(x) = 2^{5/2} pi^{-5/2} x (1-x) / ((1+x)(2x-1)^4).

`phi(x*) = 0` exactly (with `(1+sqrt3)^2 = 2(2+sqrt3)`, `(sqrt3-1)^2 = 2(2-sqrt3)` all logarithms cancel), `phi'(x*) = 0`, so
`N phi(k0/N) = O(1/N)` and `a(k0) sqrt N -> psi0(x*) = sqrt3/(sqrt2 pi^{5/2})`.

(B3)–(B5): `T_N -> psi0(x*) sqrt(pi/12) = sqrt2/(4 pi^2)`. With (A): `sum_n t(n) = 256 sqrt2/(4 pi^2) = 64 sqrt2/pi^2`. ∎

*Check.* `proofs/verify_T11_sun_1pi2_sqrt2.py`: Part A exactly (ratios, both rational identities, all cases of the WZ
equation on a grid, `gam(n,0) = -t(n)/256`, the finite identity (A) for N = 1,2,3,5); Part B: `phi' = L`, `phi''(x*) = -24`, `psi0(x*) sqrt(pi/12) = sqrt2/(4 pi^2)` and the polynomial identity of (B1) exactly (`phi(x*) = 0` to 55 digits; the exact argument is in (B5)); log-concavity on `I`,
monotonicity below (B1) in exact rational arithmetic for every N from 4 to 400 (all indices, including the junction and the
odd-N index), `a(k0) sqrt N/psi0(x*) -> 1` and `T_N/(a(k0) sqrt(pi N/12)) -> 1` numerically for N = 40 … 20000;
`T_N = sqrt2/(4 pi^2)` to 27 digits; the series to 110 digits.

---

## Theorem 12 (added after the audit).  `sum_{n>=0} (28n+17) (1)_n (3/4)_n (5/4)_n / ((3/2)_n^3 (-48)^n) = 32 log(27/16)`

Equivalently `sum_{k>=1} (28k-11) C(4k,2k) / (k (4k-1) (-48)^k C(2k,k)^2) = (4/3) log 2 - log 3`.

*Source.* **Conjecture of J. Guillera**, *Ramanujan series with a shift*, arXiv:1802.02023 (J. Aust. Math. Soc. 107 (2019)
367–380), Table 1 (level 2, `q = -e^{-pi sqrt 9}`, `a = 3 sqrt3/16`, `b = 28 sqrt3/16`, `z = -1/48`) and Table 2
(`-i G_2(1/2,q) = (3/2) ln 3 - 2 ln 2`): it is Ramanujan's series
`sum (28n+3)(1/2)_n(1/4)_n(3/4)_n/(n!^3 (-48)^n) = 16/(pi sqrt3)` with every Pochhammer index shifted by 1/2, in the
normalisation `sqrt|z| sum_n [(1)_n(3/4)_n(5/4)_n/(3/2)_n^3] (a + b/2 + b n) z^n = V`. It is the companion of
Theorem 2 (same table, `z = 1/9`). Guillera's paper proves only level 4; Cohen–Guillera (arXiv:2101.12592, Sect. 6)
list these values as observations. All 32 tabulated values (levels 1, 2, 3) are re-verified to 75+ digits in
`survey/verify_G.py` (one misprint: for `z = 1/81` the value is `pi/2 - 4 arctan(1/(2 sqrt2))`, not `pi/2 + ...`).

*Pair.* For integers `n, k >= 0`

    psi(n,k) = (2k+2n+3/2) 3^(1-n) (-1)^n n! Gamma(2n+3/2) Gamma(k+n+1) Gamma(3k+2n+3/2)
               / ( Gamma(k+n+3/2) Gamma(n+1/2) Gamma(3k+4n+4) ),
    sigma(n,k) = (k+n+1)(4k+4n+7)(6k+4n+3)(6k+4n+5)(6k+4n+7) / (4 (2k+2n+3)(3k+4n+4)(3k+4n+5)(3k+4n+6)(4k+4n+3)),
    rho(n,k)   = -(n+1)(4n+3)(4n+5)(k+n+1)(4k+4n+7)(6k+4n+3)(6k+4n+5)
                 / (12 (2n+1)(2k+2n+3)(3k+4n+4)(3k+4n+5)(3k+4n+6)(3k+4n+7)(4k+4n+3)),
    R(n,k)     = (36k^3+(156n+120)k^2+(232n^2+362n+138)k+112n^3+264n^2+203n+51) / (3 (2n+1)(3k+4n+4)(4k+4n+3)),   gamma = R psi.

Origin: Au's `Seed4(1/2, b, k+1/2+n)` (Dougall's 5F4 sum with `(b/3, b/3+1/3, b/3+2/3)`) vanishes at `b = -n` through the
factor `1/Gamma(b)`; `psi` is the coefficient of `b+n` (same ratios, same certificate), as in Theorem 5.

*Proof.*
1. (WZ) and the compatibility hold (script; `certificates/T12.gp`), so `psi(n+1,k) - psi(n,k) = gamma(n,k+1) - gamma(n,k)`.
2. Bounds (all coefficientwise, script): `0 < sigma(n,k) <= (k+3n+3)/(k+3n+5)`, `-1/48 <= rho(n,0) < 0`, `0 <= R(n,k) <= k+2`.
   Hence `|psi(n,k)| <= |psi(n,0)| (3n+3)(3n+4)/((k+3n+3)(k+3n+4))` (telescoping product),
   `sum_k |psi(n,k)| <= (3n+4)|psi(n,0)|`, `|psi(N,0)| <= (3/8) 48^{-N}`, and
   `|gamma(n,K)| <= (K+2)|psi(n,0)| (3n+3)(3n+4)/((K+3n+3)(K+3n+4)) -> 0` as `K -> oo`, for every `n >= 0`.
3. The telescoping lemma with `K -> oo`, then `N -> oo`: `sum_{n>=0} gamma(n,0) = B := sum_{k>=0} psi(0,k)`, where
   `psi(0,k) = (3/4)(4k+3) k! Gamma(3k+3/2) / (Gamma(k+3/2) (3k+3)!) > 0`.
4. `B = log(27/16)`.
   (a) `k!/Gamma(k+3/2) = pi^{-1/2} int_0^1 s^k (1-s)^{-1/2} ds` (Beta integral).
   (b) For `|w| <= 1`: `E(w) := sum_{k>=0} (4k+3) Gamma(3k+3/2)/(3k+3)! w^{3k+3} = (4 sqrt(pi)/3) [1 - (1/3) sum_{omega^3=1} (1+omega w)(1-omega w)^{1/2}]`.
       Indeed the coefficient of `w^j` in `(1+omega w)(1-omega w)^{1/2}` is `omega^j c_j`, `c_j` the coefficient in
       `(1+w)(1-w)^{1/2}`; summing over the three cube roots of unity leaves `3 c_j` for `3 | j` and `0` otherwise;
       `c_0 = 1`, and `c_{3k+3} = -(3/(4 sqrt(pi))) (4k+3) Gamma(3k+3/2)/(3k+3)!` from the binomial series
       (`(1-w)^{1/2} = sum_j e_j w^j`, `e_j = -Gamma(j-1/2)/(2 sqrt(pi) j!)`, so `c_j = e_j + e_{j-1} = -(4j-3) Gamma(j-3/2)/(4 sqrt(pi) j!)`
       for `j >= 1`; also checked exactly for `k < 40` in the script). The series
       converges absolutely for `|w| <= 1`.
   (c) All terms being positive, sum and integral may be interchanged (monotone convergence); with `s = w^3`:
       `B = (3/(4 sqrt(pi))) int_0^1 (1-s)^{-1/2} E(s^{1/3}) ds/s = int_0^1 sum_{omega^3=1} [ (1-w^3)^{-1/2} - (1+omega w)(1-omega w)^{1/2}(1-w^3)^{-1/2} ] dw/w`.
   (d) `int_0^1 3[(1-w^3)^{-1/2} - 1] dw/w = int_0^1 [(1-v)^{-1/2} - 1] dv/v = 2 log 2` (`v = w^3`, then `v = 1-y^2`).
   (e) For each cube root of unity `omega`: `(1-omega w)^{1/2}(1-w^3)^{-1/2} = q_omega(w)^{-1/2}`,
       `q_omega(w) = 1 + omega w + conj(omega) w^2` (the product of the other two factors of `1-w^3`), with the branch
       `s(w) = sqrt(q_omega(w))` continuous on `[0,1]`, `s(0) = 1`. The function
       `A_omega(w) = (2 + omega w + 2 s)(omega + 2 conj(omega) w - 2 omega s)` satisfies
       `A_omega'/A_omega = (1 - (1+omega w)/s)/w` (script: exact, using only `s^2 = q_omega`), and it does not vanish on
       `[0,1]` (if the first factor vanished then `4 q_omega = (2+omega w)^2`, i.e. `3 conj(omega) w^2 = 0`, but
       `A_omega(0) = -4 omega`; if the second vanished then `4 omega^2 q_omega = (omega + 2 conj(omega) w)^2`, i.e. `3 omega^2 = 0`).
       Hence `Re int_0^1 [1 - (1+omega w)/s] dw/w = log|A_omega(1)| - log|A_omega(0)| = log 3 - log 4`
       (`A_1(1) = (3+2 sqrt3)(3-2 sqrt3) = -3`; for `omega != 1`, `q_omega(1) = 0` and `|(2+omega)(omega + 2 conj(omega))| = 3`).
   Since `B` is real, (c), (d), (e) give `B = 2 log 2 + 3 log(3/4) = log(27/16)`.
5. `gamma(n,0) = t(n)/32`, `t` the summand: `psi(0,0) = 3/8`, `R(n,0) = (28n+17)/(12(2n+1))`, and
   `gamma(n+1,0)/gamma(n,0) = t(n+1)/t(n)` exactly (script). Therefore `sum_n t(n) = 32 log(27/16)`. ∎

*Check.* `proofs/verify_T12_guillera_shift_48.py` (exact: (WZ), compatibility, the bounds of step 2, step 5, the
coefficients in 4(b) for `k < 40`, the antiderivative in 4(e); numerics: every step of 4 to 18–38 digits, the finite
identity `sum_{n<N} gamma(n,0) = B - sum_k psi(N,k)` for `N = 1, 2, 4`, and the final identity in both forms to 110 digits).

*Novelty.* Canonical data: `z = -1/48`, parameters `(1/4)(3/4)(1) / (1/2)^3` mod 1. The identity is stated only as a
conjecture/observation in Guillera 2018 and Cohen–Guillera 2021; it is not in Sun's lists; the three papers citing
Guillera's (OpenAlex, 2026-11-24: Cohen–Guillera; Campbell, arXiv:2310.05112, on unshifted level-3 series; one
2026 paper on irrational series) do not prove it; full-text search of our ~110 downloaded sources for the polynomial
`28n+17` / `28k-11` and its shifts: no occurrence. Like Theorem 2 it comes from a Dougall-5F4 specialisation (Seed4)
with step 1/3 in the parameters; in the integer-pattern Dougall family of Chu–Zhang it can only appear in trisected
form (rate `-1/48^3`, blocks of three terms; `code/chu_family_test.py ... G` tests this). Since Theorem 2 turned out to be exactly such a trisected example (Chu, Kodai Math. J. 41 (2018), Example 19), we
read that paper in full: its 40 examples are series for `pi`, `pi^2`, `1/pi`, `zeta(3)`, `G` with rates `1/9, -1/27, 1/64,
-1/8, 1/4, 1/729`; none has rate `-1/48` or a logarithm as value. Structurally, Theorem 2 is the trisection of the
pattern `[3;1,1,1,0]` at a point where the Omega-sum is elementary, whereas Theorem 12 corresponds to `[6;2,2,2,3]` at
`(3/2; 1/2, 5/6, 7/6, 1)` taken at first order at a pole of the family, and its boundary sum `B` is not a Dougall
evaluation (step 4 is a separate argument). Chu–Zhang 2014 (unread; by its abstract: series for `pi`, `zeta(3)`, `G`) and
Chu 2021 (unread) remain a residual risk. "Proof apparently new", not a guarantee.

*Remark (what Theorem 12 is in Chu's language).* `code/match_multisect.py` confirms exactly that the series of Theorem 12,
grouped in blocks of three terms, is the `n`-series of the integer-pattern Dougall pairs
`Dougall5F4(1/2, -m, 1/3-m, 2/3-m; k+1/2+3m)` (iteration pattern `[6;2,2,2,3]` at `Omega(3/2; 1/2, 5/6, 7/6, 1)`; this is our
pair with `n = 3m`) and `[5;1,1,2,2]` at `Omega(5/3; 2/3, 1, 7/6, 7/6)`; i.e. it is a "trisection series" in the sense of
Chu's Kodai paper, for patterns larger than those used there, at base points whose Omega-sum is not given by Dougall's
theorem (step 4 above evaluates it).

*Further members of Guillera's tables.* `code/gscan/` (signature scan + exact matching for all 32 tabulated series):
exact WZ pairs exist for five of them: `z = 1/9` (Theorem 2), `z = -1/48` (Theorem 12), and
- `z = -1/4`, level 2 (value `ln 2`): matched by Gauss, Dixon, Dougall, ... seeds with integer patterns — it lies in the
  Dougall family and is therefore very probably already in the Chu–Zhang list; **not claimed**;
- `z = -27/512`, level 1 (value `(3 sqrt3/8) ln 2`): matched by Seed1 / Gauss's second theorem / Bailey, i.e. the
  seed family of Chu 2020 (the paper that contains Theorem 5); **not claimed** until that paper has been checked;
- `z = -9/16`, level 3 (value `(sqrt3/4)(3 ln 3 - 2 ln 2)`, i.e.
  `S := sum_{n>=0} (10n+7)(1)_n(5/6)_n(7/6)_n/(3/2)_n^3 (-9/16)^n ?= 8 log 3 - (16/3) log 2`): matched only by Ebisu-type
  seeds. **Partial result, not a proof** (`code/g13_reduction.py`): the pair from `2F1(a,1-2a;3/2-2a;-1/3)`, `a -> -n`,
  `k -> k-n`, is finite on the triangle `k < n` with poles on `n <= k <= 2n`; telescoping on the triangle gives the exact
  relation `S = 28/9 + (4/9) U`, where `U = sum_{n>=1} (6 + 2/(2n+1) - 1/(n+1)) (5/6)_n(7/6)_n/((5/4)_n(7/4)_n) (3/4)^n`
  is the series of diagonal terms (WZ equation exact; the triangle identity, the constants and `T_N -> 0` checked
  numerically to 38 digits, not written out as a proof). So Guillera's value is equivalent to
  `U = 18 log 3 - 12 log 2 - 7` (true to 45 digits), a series at `3/4` that we did not evaluate: a second pair for `U`
  (Ebisu's `2F1(3a,1-a;3/2;9/8)`) leads to a divergent diagonal series. **Open.**
The other 27 tabulated series (all those with larger discriminants) have no signature hit in our 169 seeds.

---

## Theorem 13 (added 2026-11-25).  `sum_{n>=0} P(n+1/2) n!^16 (4n+2)! / (2n+1)!^10 = 2048 pi^4/15`, `P(x) = 43680x^4+20632x^3+4340x^2+466x+21`

Equivalently `sum_{n>=0} P(n+1/2) (1)_n^7 (3/4)_n (5/4)_n / ((3/2)_n^9 4096^n) = 1024 pi^4/15`; equivalently, in the notation
of Cohen–Guillera, `S(H,a,P;1/2) = (2048/15)/pi^4` for Cullen's series
`sum (1/2)_n^7 (1/4)_n (3/4)_n/(1)_n^9 P(n)/4096^n = 2048/pi^4` (the value of the same series with every Pochhammer index
and the argument of `P` shifted by 1/2).

**Not a new statement.** In binomial form the identity is eq. (2) of Deyi Chen, MathOverflow question 486961 (update of
1 February 2025): `sum_{k>=0} 30(4k+1)(5460k^4+13499k^3+12601k^2+5265k+831) C(4k,2k)/((2k+1)^9 C(2k,k)^8) = (4 pi)^4`
(the two forms agree term by term: `P(k+1/2) = 8(5460k^4+...+831)` and `k!^16 (4k+2)!/(2k+1)!^10 = 2(4k+1) C(4k,2k)/((2k+1)^9 C(2k,k)^8)`).
J. Guillera answered there (1 Feb 2025): "the proof of your second formula can be obtained replacing n with n+1/2 in
the function G(n,k) of the WZ pair found by Kam [Au]". We know of no written proof. The proof below shows what that
remark hides: at `n -> n+1/2` Au's pair, suitably normalised, gives only the Amdeberhan–Zeilberger series (`T3 = 0`);
the identity in question is the *derivative* in the shift at 1/2 and needs, in addition, the harmonic-number identity of
Sun–Au–Wei. We found the MathOverflow post after writing the proof: it was listed in our own `survey/SURVEY-C.md` (entry
for MO 486961) in binomial form, and the check-first search had used the Pochhammer form and the value `(2048/15)/pi^4`.
So: **statement known (Chen 2025; Cohen–Guillera 2021 as an observation), proof route indicated (Guillera 2025); what is
ours is only the written-out proof.**

*Source.* H. Cohen and J. Guillera, *Rational hypergeometric Ramanujan identities for 1/pi^c: survey and
generalizations*, arXiv:2101.12592, Section 6, item 13: `S(H,a,P;1/2) = (2048/15)/pi^4`, stated as a numerical
observation (the section lists the values at `x = 1/2` as observations; for `1/pi` they go back to Guillera's
conjectures of 2018, cf. Theorems 2 and 12). Cullen's series itself (J. Cullen 2010, conjecture) was proved by
K. C. Au, arXiv:2312.14051, Example I.

*Pair.* For real `n > 0` and integers `k >= 0` let

    E(n,k) = (2k+3n+1/2) Gamma(n+1/2)^10 Gamma(k+n+1/2)^5 Gamma(k+4n) / ( Gamma(2n)^5 Gamma(k-n+3/2) Gamma(k+2n+1)^5 )

(with `1/Gamma` entire), and let `R(n,k)` be the rational function in `certificates/T13.gp` (numerator of degree 9 in `k`,
81 monomials, of degree 6 in `k` in the denominator). Then, as rational identities,

    -E(n+1,k) - E(n,k) = R(n,k+1) E(n,k+1) - R(n,k) E(n,k),                                   (WZ13)
    R(n,0) = (2n-1) P(n) / (16384 n^4 (6n+1)),      lim_{k->oo} R(n,k)/k^3 = -(820n^2+180n+13)/(8192 n^5).

So for every real `a > 0`, `U(m,k) = (-1)^m E(a+m,k)`, `V(m,k) = (-1)^m R(a+m,k) E(a+m,k)` is a WZ pair in `(m,k)`.
`E` is Au's function `F` of Example I (his Dougall-7F6 seed in the direction `(-1,-1,-1,-1,-1;2)` at the base point
`(1/2,...,1/2)`) divided by `Gamma(1-2n)^3 Gamma(2n)^3 / (Gamma(n+1/2)^7 Gamma(1/2-n)^7)`, which equals `cos^4(pi n)/(8 pi^4 sin^3(pi n))`
and changes sign under `n -> n+1`; `R` is his certificate.

*Proof.*
1. *The identity for generic `a`.* Put `v(n) = lim_k R(n,k) E(n,k) = -(820n^2+180n+13) Gamma(n+1/2)^10 / (4096 n^5 Gamma(2n)^5)` and

       T1(a) = sum_{k>=0} E(a,k),   T2(a) = sum_{m>=0} (-1)^m v(a+m),   T3(a) = sum_{m>=0} (-1)^m R(a+m,0) E(a+m,0).

   Au's identity `sum_k F(a,k) + sum_n g(n+a) = sum_n G(n+a,0)` (arXiv:2312.14051, Example I, "for a generic `a`";
   it is his Proposition 2.1, the vanishing of `lim_n sum_k F(n+a,k)` being his Proposition 3.4), divided by the
   antiperiodic factor at `n = a`, reads `T1(a) + T2(a) = T3(a)`. **This step is quoted from Au**; we re-checked it only
   numerically (30 digits at `a = 0.3, 0.45, 0.55, 0.7`).
2. *Analyticity near `a = 1/2`.* On the disc `D = {|a - 1/2| <= 1/10}` every term of the three series is holomorphic
   (`Gamma` at arguments with positive real part, `1/Gamma` entire, the denominators of `R(n,0)` and `v` do not vanish), and
   the series converge uniformly on `D`: the terms of `T1` are `O(k^{-3})` uniformly (Stirling: `E(a,k) ~ 2 c(a) k^{-3}`,
   the exponent does not depend on `a`), and the ratios of consecutive terms of `T2` and `T3` are fixed rational functions of
   `n = a+m` tending to `-1/1024` and `-1/4096` as `|n| -> oo`. Hence `T1, T2, T3` are holomorphic on `D` (Weierstrass),
   may be differentiated termwise, and by step 1 and continuity `T1 + T2 = T3` on `D`.
3. *At `a = 1/2`.* `E(1/2,k) = 2/(k+1)^3` and
   `d/dn log E(n,k)|_{n=1/2} = 3/(2k+2) + 6 psi(k+1) - 6 psi(k+2) = -9/(2(k+1))`, so `T1(1/2) = 2 zeta(3)`, `T1'(1/2) = -9 zeta(4)`.
   `(-1)^m v(m+1/2) = -(1/32) A_m (205m^2+250m+77)`, `A_m = (-1/1024)^m (1)_m^5/(3/2)_m^5`, and
   `d/dn log v(n)|_{n=m+1/2} = (410m+250)/(205m^2+250m+77) - 10 (H_{2m+1} - H_m)`; hence
   `T2'(1/2) = (10/32) sum_m A_m [ (205m^2+250m+77)(H_{2m+1}-H_m) - 41m - 25 ] = (10/32) * 48 zeta(4) = 15 zeta(4)`
   by the identity `sum_m A_m [(77+250m+205m^2)(H_{2m+1}-H_m) - 41m-25] = 48 zeta(4)` (conjectured by Z.-W. Sun; proved by
   K. C. Au, Ramanujan J. 66 (2025) = arXiv:2212.02986, and by C. Wei, arXiv:2303.07887).
   Finally `R(n,0) E(n,0) = -P(n) Gamma(n+1/2)^15 Gamma(4n) / (16384 n^4 Gamma(2n)^5 Gamma(2n+1)^5 Gamma(1/2-n))` vanishes at every
   `n = m+1/2` (so `T3(1/2) = 0`, and step 2 at `a = 1/2` is the Amdeberhan–Zeilberger series `sum A_m(205m^2+250m+77) = 64 zeta(3)`),
   and since `d/dn [1/Gamma(1/2-n)]_{n=m+1/2} = -(-1)^m m!`:
   `(-1)^m d/dn [R(n,0)E(n,0)]_{n=m+1/2} = P(m+1/2) m!^16 (4m+2)! / (2048 (2m+1)!^10)`.
4. Differentiating `T1 + T2 = T3` at `a = 1/2`:
   `(1/2048) sum_m P(m+1/2) m!^16 (4m+2)!/(2m+1)!^10 = -9 zeta(4) + 15 zeta(4) = 6 zeta(4) = pi^4/15`. ∎

*Check.* `proofs/verify_T13_cullen_half_shift.py`: (WZ13), its compatibility condition, `R(n,0)` and the leading coefficient
of `R` exactly (PARI; the certificate is recomputed by Gosper's algorithm and written to `certificates/T13.gp`); every
formula of step 3 numerically (25–40 digits); `T1+T2 = T3` at four values of `a`; `T3'(1/2) = pi^4/15` to 45 digits by
numerical differentiation; the two cited identities to 35 digits; the final identity to 120 digits.

*What is and is not proved here.* The theorem is a corollary of three published results (Au's parametric identity for
Cullen's series; the Amdeberhan–Zeilberger series; the Sun–Au–Wei harmonic-number identity) plus the expansion at the
half-integer point carried out above. Step 1 is quoted, not re-proved. No Lean.

*Novelty.* None claimed for the statement (see "Not a new statement" above). The same MathOverflow question contains, as
eq. (3), the half-shift of Zhao's series as a `zeta(3)` series,
`sum (-27/256)^m (1)_m^5 (3/4)_m (5/4)_m (5/6)_m (7/6)_m/(3/2)_m^9 (2264m^4+6118m^3+6267m^2+2884m+503) = 384 zeta(3)`
(= Cohen–Guillera item 12, `9216 zeta(3)/pi^7`), without proof. By the method above it follows at *leading* order from Au's
pair for Zhao's series (direction `(1,0,0,0,0;1)`): there `T1(1/2) = 2 zeta(3)/pi^4`, `T2(1/2)` is Guillera's proved series
`sum (-1/4)^m (1)_m^5/(3/2)_m^5 (10m^2+14m+5) = 7 zeta(3)/2` and `T3(1/2)` is the shifted Zhao series (checked numerically,
`code/cg12_explore.py`; not written up as a theorem, since statement and route are known).

---

## What resisted, and why

* **The famous rational 1/pi^2 series** (Guillera `7^-4`, `80^-3`; Almkvist–Guillera `(3/5)^6`; `1640/2^-10`; golden):
  no signature hit for any of the 169 seeds in the scanned boxes. For `7^-4`, `80^-3`, `(3/5)^6` the exponent is not in the
  group generated by the `m^m` and the bases of the seeds (Au's obstruction); new seeds with other bases would be needed.
* **Gourevitch (1/pi^3), Cullen and Zhao (1/pi^4)**: already theorems (Au, arXiv:2312.14051, 2023). Not ours.
* **A2** (`= 20 zeta(3)`) and **A3** (`= K/2`): signature hits but no exact summand match (A2), no hit (A3). This now
  includes the multisection search (`code/multisect_scan.py`, `code/match_multisect.py`): for A2, 631 signature hits
  with blocks of 1, 2 or 3 terms in the Dougall, Phi, Dixon, Watson, Whipple, Seed2, Seed3 families (boxes 3–6), none
  of which matches the summand exactly; for A3 no signature hit at all for q = 1, 2, 3.
* **Guillera's half-shifted series** (32 conjectured values, arXiv:1802.02023): two are proved here (one of them, T2, was
  already in Chu 2018), two more are probably known (Chu), one is reduced to another series (see the end of Theorem 12),
  and for the remaining 27 (all larger discriminants, e.g. `z = 1/7^4`, `1/99^2`, `-1/882^2`, the Chudnovsky-type
  `-1/640320^3` shift) no seed has the right exponent `z`.
* **Harmonic-number companions** (B1–B79): B65, B46, B35 and the second identity of Conj. 5.25(ii) are done (Theorems 6, 8, 9, 10).
  Second-order companions (B36, B45, B27, B48) need second derivatives, whose squares of first-order harmonic numbers do not
  cancel in the single pairs we have; B66 needs a second pair for the kernel of Theorem 1. Not done.

## Independent audits

1. A referee pass over Theorems 3, 4 and 11 (every step of Part B re-derived; (B1) tested in exact arithmetic for all
   N <= 1000; the Stirling expansion re-derived in sympy; identity (A) checked exactly for N <= 40) found no error;
   its small textual repairs are incorporated above. Report: `proofs/AUDIT-T3-T4-T11.md`.
2. A second, fully independent audit (`AUDIT-INDEPENDENT.md`, scripts `audit2/`) re-derived all nine certificates on two
   exact paths, checked the pairs term by term from the closed forms, re-read every limit argument, re-summed all
   eleven series to 228–235 digits and rebuilt the Lean files. Verdict: **all eleven proofs correct; T1 fully
   machine-checked**; but (F1) T5 is Chu 2020 (7.10); (R1) T2 is Guillera's conjecture; (R2) the novelty method was
   blind to equivalent forms; (R3) two citations for T3 were wrong; (R4) the uniformity/termwise-differentiation
   step of T6, T8, T9, T10 was only asserted. All five points are repaired in this revision (R4:
   "Differentiation lemma" and `proofs/verify_R4_uniform_bounds.py`).
These are readings by independent agents, not formal proofs; T1 and T4 are formalized end to end.

## How novelty was checked (revised 2026-11-24)

*Statements.* The latest arXiv versions of Sun's lists (2307.03086v10, 2506.01870v6, 2603.29973v3) state T1, T3, T4,
T6, T8–T11 as conjectures; T3 is restated as open in arXiv:2609.04176 (Sep 2026); T7 is MathOverflow 506695 (one
answer, on partial sums only); T2/T9 are MathOverflow 456779 (no answer).

*Canonical data (the method that replaces the string search).* For each identity: argument `z`, multiset of
Pochhammer parameters mod 1, and the polynomial up to `k -> k+j`:

| id | z | upper / lower parameters mod 1 | closest known results | verdict |
|---|---|---|---|---|
| T1, T6 | 1/4 | `(1/3)(2/3)(1) / (1/6)(1/2)(5/6)` | Sun's own companions; not in Le Lan arXiv:2609.22304 (Beta integrals) | not found |
| T2 | 1/9 | `(1/4)(3/4)(1) / (1/2)^3` | **Chu, Kodai Math. J. 41 (2018), Example 19**; Guillera 2018 (conjecture), Cohen–Guillera 2021 (observation), Au 2026 Ex. VIII (same WZ family, other point) | known |
| T9 | 1/9 | as T2, with `5H_{2k-1}-4H_{k-1}` | no harmonic-number version found (Chu 2018 has none; Au's Ex. VIII expands at another point) | not found |
| T12 | -1/48 | `(1/4)(3/4)(1) / (1/2)^3` | Guillera 2018 (conjecture), Cohen–Guillera 2021 (observation); not in Chu 2018 | not found |
| T3 | 1/9 | `(1)^2 / (1/3)(2/3)` | Sun arXiv:2609.04176 (still open) | not found |
| T4, T8, T10 | 27/32 | `(1/6)^2 (1/2) (5/6)^2 / (1/3)(2/3)(1)^3` | Au's and Chu 2020's only `27/32` kernel is `(1/6)(5/6)/(1)^2` (a vanishing sum); Chu 2011 Ex. 92 has the same base point but rate `(3/8)^6` | not found |
| T5 | -27/512 | `(1/3)(2/3)(1)^3 / (1/4)^2 (1/2) (3/4)^2` | **Chu 2020, (7.10)** | known |
| T7 | 4/27 | `(1/2)(1) / (1/6)(5/6)` | MO 506695 (partial sums only) | not found |
| T11 | 1/64 | `(1/4)^3 (3/4)^3 / (1/2)(1)^5` | not among the kernels of Au (2312.14051, 2403.04555, 2602.08721), Chu (2103.07872, 2104.01599), Hou–Sun 2604.15172, Hou–Mu 2604.22377 | not found |

(Computed exactly by `code/canonical_data.py`, output `code/canonical_data.out.txt`; classes common to numerator and
denominator are cancelled, so rational factors such as `1/(2k+1)` are absorbed.)

*Literature read for the seed families concerned* (by the independent auditor, full text): Chu, Math. Comp. 80 (2011)
(Gauss 2F1 / Dougall 2H2; 100 examples); Chu, Acta Math. Sci. 40B (2020) (Gauss/Bailey 2F1(1/2)); Chu,
arXiv:2103.07872, 2104.01599; all displayed series of Au's three papers; Guillera 1802.02023, 2503.00570;
Cohen–Guillera 2101.12592; Hou–Mu, Le Lan, Patel, Zhou, Zuniga, Shalev, Zudilin, Campbell, Xu, Wei–Xu and others
(about 60 sources; list in `AUDIT-INDEPENDENT.md`, A.1), with full-text search for each polynomial *and its index
shifts* `P(k+j)`, `|j| <= 3`, and for the Pochhammer patterns.

*What could not be read.* Chu–Zhang, *Accelerating Dougall's 5F4-sum and infinite series involving pi*, Math. Comp. 83
(2014) 475–512 (125 examples), Chu, *Infinite series identities derived from the very well-poised Omega-sum*,
Ramanujan J. 55 (2021) 239–270, and (by us; the independent auditor had it) Chu, Math. Comp. 80 (2011) 2223–2251
(Dougall's bilateral 2H2-series, 100 examples). Routes tried: AMS (bot check, also from a headless browser), Springer
(abstract and reference list only), Unpaywall / OpenAlex / CORE / Semantic Scholar (no open copy), zbMATH (a review of
the first only), ResearchGate (blocked); a Wayback Machine snapshot of the AMS PDF of Chu–Zhang 2014 exists
(`web.archive.org/web/20250906133047/https://www.ams.org/journals/mcom/2014-83-285/S0025-5718-2013-02701-9/S0025-5718-2013-02701-9.pdf`)
but that host is blocked by the egress policy of the machine this was written on.

*What their method is* (read in full: Guillera, arXiv:1611.04385, "we show how to prove the examples of a paper by Chu
and Zhang using the WZ-algorithm", and Chu, Kodai Math. J. 41 (2018) 512–530). They accelerate
`Omega(a;b,c,d,e) = sum_k (a+2k) [b,c,d,e]_k/[1+a-b,1+a-c,1+a-d,1+a-e]_k` (the partial sum of Dougall's 5H5) and its
limit `Phi(a;b,c,d)` (`e -> oo`, alternating) by "iteration patterns" `[n_a,n_b,n_c,n_d,n_e]` of integers
(`[3,0,1,2,2]`, `[3,0,2,2]`, `[3,1,1,1,3]`, `[2,0,1,1,1]`, `[2,0,1,0,1]`, `[3,1,1,1,1]`, `[3,0,1,1,1]`, `[4,1,1,1,1]`, ...)
at base points where the sum is known (Dougall's theorem, its limits, `zeta(3)`, `G`), sometimes followed by
**bisection or trisection** of the resulting series. In WZ language this is the family
`Dougall5F4(a+An, b+Bn, c+Cn, d+Dn; k+x+Kn)` with a free offset `x` (the offset turns the 5F4 term into the 5H5 term
with `(a';b',c',d',e') = (a+2x; b+x, c+x, d+x, a+x)`), resp. the well-poised `4F3(-1)` seed for `Phi`, resp. Gauss's
2F1(1) seed for the 2H2 paper; bisection/trisection correspond to taking blocks of 2 or 3 consecutive `n`.

*Independent enumeration* (`code/chu_family_test.py`, logs `code/chu_family_*.log.txt`): all integer patterns in a box
(Omega: |shifts| <= 4, also 5 for T4/T11; Phi and 2H2: <= 6), all base points and offsets with denominators dividing
24, compared on canonical data (`z`, Pochhammer classes mod 1) with each identity **and with its bisection and
trisection**. Positive controls that the enumeration must find, and does: Guillera's three `1/pi^2` series (rates `-1/4`,
`1/16`, `27/64`), Chu–Zhang's Examples 42 and 45 (as quoted by Guillera), Ramanujan's `1/4` and `1/64` series, and
T2 (trisected; = Chu 2018, Example 19). A hit is only a necessary condition (same `z` and parameter classes).

| id | Omega family | Phi family | 2H2 (Gauss) family | exact summand match known (`OPEN-LIST.md`) | our reading |
|---|---|---|---|---|---|
| T11 | none | none | none | Au's Seed9, Dixon6F5 only | not in the Chu families in any form: lowest risk |
| T4 | none | bisected only | bisected only (it *is* a half-integer Gauss pattern) | Seed1, Gauss2nd, Bailey | not in Omega; in Chu 2011 / Chu 2020 only if a large pattern with bisection was used; the auditor found nothing with rate 27/32 for this kernel in either paper |
| T12 | trisected | none | none | Seed4, Seed10, Ebisu | same structure as T2, but the value is a logarithm and the base sum is not a Dougall evaluation |
| T3 | trisected | none | none | Seed4, Ebisu | needs a derivative in the offset; Zuniga (2026), who works with Guillera's Dougall–WZ program, still treats the sister series as conjectural |
| T7 | plain, bisected | plain, bi-, trisected | plain, bi-, trisected | Dixon, Watson, Seed2, Seed3 | lies in classical families at small patterns; a vanishing sum is unlikely to be listed as a "pi-formula", but this is not a guarantee |
| T1 | bisected | plain, bi-, trisected | plain, bi-, trisected | Gauss2F1, Kummer, WP4F3(-1), Watson, Seed2, Seed3, Ebisu | **lies in the seed families of Chu 2011 and Chu–Zhang 2014** (at base points where the seed degenerates, i.e. as a limiting case); whether it is among their 225 examples we could not check; the auditor's full-text search of Chu 2011 found no occurrence |

*Exactly what the enumeration covers* (so that a reader with library access knows what remains to be checked by eye).
Coordinates: `F(n,k) = Seed(a+An, b+Bn, ...; k+x+Kn)`; for the Omega family Chu's iteration pattern is
`[d_a; d_b, d_c, d_d, d_e] = [A+2K; B+K, C+K, D+K, A+K]` at the base point `(a+2x; b+x, c+x, d+x, a+x)`.
- Omega (Dougall 5F4 / 5H5): all `(A,B,C,D,K)` with `|A|,|B|,|C|,|D| <= 4`, `0 <= K <= 4` (for T4, T11 plain and
  bisected also `<= 5`); this contains every pattern with entries in `0..4` (and many more, e.g. `[5;1,1,2,2]`,
  `[6;2,2,2,3]`), in particular all patterns quoted from Chu–Zhang and Chu: `[2,0,1,1,1]`, `[2,0,1,0,1]`, `[3,0,1,1,1]`,
  `[3,0,1,2,2]`, `[3,1,1,1,1]`, `[3,1,1,1,3]`, `[4,1,1,1,1]`. Phi (`e -> oo`): `|A|,|B|,|C| <= 6`, `0 <= K <= 6`.
  2H2 (Gauss): the same box 6. Dougall 7F6 seed: box 4, plain and bisected only.
- Base points and offsets: all values in `(1/24) Z` mod 1 (for 7F6: `(1/12) Z`), i.e. denominators 1, 2, 3, 4, 6, 8, 12, 24.
- Forms: the series itself, and the series grouped in blocks of 2 and of 3 consecutive terms.
- NOT covered: denominators 5, 7, 9, 16, ...; patterns outside the boxes; multisections of order >= 4; further
  confluent limits (two parameters to infinity other than through the 2H2 family); series obtained as linear
  combinations of two accelerations or from Chu's "dual"/inverse-series relations (Gould–Hsu); q-analogues.

*Check of the Chu papers that are available elsewhere* (done by the coordinating session's auditor, page by page, and
relayed to us; not our own reading): Chu, Math. Comp. 80 (2011) (Examples 1–100), Chu, Acta Math. Sci. 40B (2020),
Chu, Ukr. Math. J. 74 (2022) (= arXiv:2103.07872), Li–Chu, ERA 31 (2023), Chu, Contrib. Discrete Math. 15 (2020),
Li–Chu, Mathematics 12 (2024): **T1, T3, T4, T7, T11, T12 are not in these six papers**, plain or multisected. Near
misses recorded there: Chu 2011, Examples 47–48, have T1's parameters `[1,1/3,2/3; 3/2,5/6,7/6]` but at rate `1/64` with
a quadratic polynomial (`4 pi/sqrt3 = 7 + sum_{k>=1} [...]_k (4+39k+63k^2)/(64^k k)`, which we re-verified numerically);
this is a different series, not a multisection of T1 (a trisection of T1 would have parameters with denominators
9 and 18); T7's parameters occur in Chu 2011 at rates `-1/4`, `-1/27` and in Chu 2020 at `2/27`, not at `4/27`; the
only rate-`27/32` material is Chu 2020, Section 5 (kernels with at most 4 parameters). So for these six identities
**the residual risk is confined to Chu–Zhang, Math. Comp. 83 (2014) and Chu, Ramanujan J. 55 (2021)**, which neither
session could open, minus what the enumeration above excludes.

T4 itself comes from Gauss's 2F1(1) theorem with a half-integer pattern, and T11 from Au's Seed9 through a
second-order pole, a route outside the Abel-lemma framework. For T1 (the identity that is fully formalized in Lean)
the honest statement is: *proof new to us and machine-checked; not in Chu 2011 (checked by the auditor); priority not
established* until someone reads Chu–Zhang 2014 with the canonical data `z = 1/4`, `(1/3)(2/3)(1)/(1/6)(1/2)(5/6)` in hand.

So: T2 and T5 are known (Chu 2018, Chu 2020); for the other ten (T1, T3, T4, T6–T12), "open as far as a serious
search can tell" — not a guarantee.
