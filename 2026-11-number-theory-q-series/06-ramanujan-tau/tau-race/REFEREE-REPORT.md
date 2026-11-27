# Referee report on TAU-RACE.md (adversarial pass, 2026-11-26)

Referee: an AI agent (Claude), instructed to find the error. Scope: `TAU-RACE.md` (all 585 lines), the
programs `cert_constants.py`, `cert_small.py`, `bounds.py`, `sym2_phi.py`, `sym2_verify.py`, `tau_table.py`,
`arbutil.py`, `independent_check.gp`, and the files in `data/`. Scratch work (not part of the repository):
`/tmp/claude-0/-home-user-weighted-prime-races/559a4317-d376-57ac-9c95-b35c412a8522/scratchpad/referee/`.

## Verdict

**Repairable gaps. No fatal error found.**

I re-derived every lemma on paper and attacked each with an independent numerical test. Theorems 1, 2, 3, 4
survive: every identity, sign, residue, multiplicity and inequality I re-derived is correct, every program
computes what the text says, and all five certificate files reproduce byte for byte. What I found is two
small gaps (one uncertified sub-interval behind a constant in Theorem 2; one mis-transcribed citation
under Theorem 4) and a list of wrong or unsupported digits and statements in the prose. None of them
touches the positivity statements `Θ_1(t) > 0 (t > t*)`, `R(x) > 0 (x > x*)`, the equivalence, or `RH_S(60)`.

Limits of this verdict are listed at the end ("Not checked").

---------------------------------------------------------------------------

## Findings

### FATAL

None.

### REPAIRABLE

**R1. Theorem 2, last sentence ("`R(x) > √x/(3 log x)` for every `x ≥ 100`"); §1 "Numerical truth"
("`N(x) ≥ 0.7951` on `[100, 9 999 991]`", "`G(t) ≥ 0.5814` on `[100, 9 999 991]`", "`G(t) ≥ 0.6257` for all
real `t ∈ [10^4, 9 999 991]`"); §4 "The constants".**
The certificate does not cover the stated ranges. `cert_small.py` accumulates a gap `[q, p]` into
`infG[a]`, `infN[a]` only if `q >= a`. For `a = 100` the first gap is `[101, 103]`; for `a = 10^4` it is
`[10007, 10009]`. The pieces `[100, 101)` (inside the gap `[97, 101]`) and `[10^4, 10007)` (inside
`[9973, 10007]`) are certified by nothing. §4 says it correctly ("gaps with `q ≥ 100`"), but Theorem 2
and §1 then claim `x ≥ 100`, `t ≥ 10^4`.
The statements are true. Fix from numbers already in `small_range_10000000.json` (`kernel_constants`):
on `[97, 101)`, `R(x) = D_0(97) − D_1(97)/x` increases (`D_1(97) = 148.36 > 0`), so
`N(x) ≥ R(97)·log 101/√101 ≥ 2.2406·0.4592 = 1.0289`; `Θ_1` increases (`Θ(97) = 12.71 > 0`), so
`G(t) ≥ Θ_1(97)/101^{3/2} = 648.8646/1015.04 = 0.6392`; on `[9973, 10007)`,
`G(t) ≥ Θ_1(9973)/10007^{3/2} = 629075.32/1001050.2 = 0.6284`. Better: change the test in the program to
"gap intersects `[a, N]`" (`p > a`) and rerun. My PARI float check: min `N` on `[100,101]` = 1.0530,
min `G` on `[100,101]` = 0.6870, min `G` on `[10^4,10007]` = 0.6320.

**R2. (F3), "In Shimura's normalisation the completed function is `π^{-3s/2} Γ(s) Γ((s−k+2)/2) D(s)` …
with `k = 12` and `s ↦ s + 11` this is `γ_S` up to a constant."**
False as written. `π^{-3(s+11)/2} Γ(s+11) Γ((s+1)/2) = 2^{s+10} π^{-5} · γ_S(s)`: the quotient is `2^s` times a
constant, not a constant. If the quoted function were symmetric under `s ↦ 2k−1−s`, then
`Λ_S(1−s)/Λ_S(s) = 2^{2s−1}`, contradicting `Λ_S(s) = Λ_S(1−s)`. So the transcription has lost a factor
`2^{-s}` (the symmetric function is `2^{-s} π^{-3s/2} Γ(s) Γ((s−k+2)/2) D(s)`). The functional equation with
exactly `γ_S = Γ_R(s+1)Γ_C(s+11)` is what Lemma V0(c), the reality of `Z`, and Lemma V3 stand on, so the
citation must be right. `γ_S` itself is the correct, standard factor (Gelbart–Jacquet), and I confirmed
it numerically three ways (PARI `lfuncheckfeq` = −67 bits; `Λ_S` from the Poisson formula equals PARI's
`lfunlambda` at six points, including one off the line; the Mellin transform of `K` equals `γ_S` at three
points). Fix: correct the quotation (I could not consult Shimura's paper; check the `2^{-s}` against the
original), or cite the `Γ_R(s+1)Γ_C(s+k−1)` form directly.

### MINOR

**M1. Theorem 4, "(the first is `1/2 + i·7.0702241402869…`)".** Wrong digits. The certified isolating interval
is `[7.070224140286882, 7.070224140290521]` (width `2^{-38}`), which supports `7.0702241402` and no more. The
quoted string is the truncated lower endpoint, not the zero. PARI: `7.07022414028743499…`.

**M2. §5.2, "each is bisected to an interval `[l_j, r_j]` of length `2^{-38}`".** False for zeros #46–#55 (all
above 60): widths `2^{-31}, 2^{-31}, 2^{-30}, 2^{-27}, 2^{-19}, 2^{-22}, 2^{-16}, 2^{-13}, 2^{-11}, 2^{-11}`
(`sym2_rh_H60.json`). The bisection hits `sm == 0` and stops. Cause: §5.1 states
`|E| ≤ 3.2·10^{-100} + 4.4·10^{-95}` but omits the dominant error, the radii of the table entries. For
`k ≥ 447` (`y > 80`, `N = 1`) the term `tail_n` is `∫_1^∞` rather than `∫_2^∞` and gives radius up to
`8.3·10^{-64}`; the resulting radius of `Z(t)` is about `2·10^{-64}` for every `t`, against `|Z(72)| ≈ 2·10^{-61}`
and `|Z(71.458)| = 2.3·10^{-64}`. Nothing is invalid (the program carries the true radii, all 1153 grid signs
are certified, the 41 zeros below 60 do have width `2^{-38}`, and Lemma V3 uses the actual intervals), but
the text misdescribes the accuracy, and the scan to 72 is close to the resolution limit. Fix: state the
table radii; bound the n-tail from `N+1`, or raise `XMAX`.

**M3. §1, "`ε_1 = 0.2/log²(10^7) = 7.6985·10^{-4}`"; Remark 4.3, "`ε_1 = 1/log 10^7 = 0.06205`".** The values
are `7.69844·10^{-4}` and `0.062042`. Both printed numbers are upper bounds, so harmless, but "=" is wrong
in the last digit. `bounds.py` uses the exact Arb values.

**M4. §1, "Certified constants (provenance in §6)".** Provenance is in §7.

**M5. §1, "(Float: the minimum of `G` over primes in `[10^4, 10^7]` is `0.625775` at `p = 9 132 733` …)".** Not
supported by any file in `data/` (the string occurs nowhere). I recomputed it in floats:
`0.6257750353` at `p = 9132733`. Correct, but unsourced.

**M6. §7, environment claims.** "Python 3.11": the installed interpreter is 3.13.18 (python-flint 0.9.0 is
right). "PARI/GP (version printed at the top of `data/independent_check.out`)": the stored file has no
version line; it is stale relative to `independent_check.gp` (rerun differs from the stored file by exactly
the line `PARI/GP version [2, 15, 4]`). "3 min" for that script: 8 s here.

**M7. §3, "`L_A` from `t ≈ 1.12·10^5`".** `bounds.out` says `10^5.046 = 111173`. Safe side, last digit off.

**M8. §1, after Theorem 4: "without it they hold with 'GRH(Δ)' replaced by 'GRH(Δ) and RH_S(60)'".**
Not true for Theorem 3 as an equivalence: (c) or (e) gives GRH(Δ) by Landau, not `RH_S(60)`. Say: without
Theorem 4, (a) + `RH_S(60)` ⟹ (b), (d), and (c) or (e) ⟹ (a).

**M9. Lemmas 2 and 3, "exactly as in `../explicit/EXPLICIT.md`, Lemma 3".** That proof is written for a
Dirichlet character with simple trivial zeros and the variable `s`. Here one needs it for `L_S'/L_S(2s)`
(horizontal lines `Im 2s = ±T_k`, double trivial zeros at odd `w ≥ 11`, left line `Re 2s = −2U − 1/2`). The
adaptation is routine and I see no obstacle, and the resulting formulas pass every numerical test I made
(below), but the estimate is cited, not written.

**M10. Program comments that do not match the code (no effect on results).** `cert_constants.py`, `c3`
docstring: "check … `psi(y) < 1.04 y`" is not checked (`psi_up` is dead code, with a comment "fix below"
and no fix). `cert_small.py` docstring: output is `small_range_<N>.json`, not `small_range.json`.
`sym2_phi.py` usage line: default `XMAX=130`; the run used 160.

---------------------------------------------------------------------------

## Answers to the ten questions of the work order

1. **Lemma 1.** Correct. Even powers `k = 2j`: `c_{2j}(p) log p = Λ_S(p^j) − log p` with local roots
   `α², 1, β²`; signs in (2.1) are right. Numerical defect of the identity at `u = 10.5, 100.5, 5000.3, 65536.5`:
   `≤ 1.2·10^{-35}`; of (2.1) at `t = 10^3, 10^5, 10^6`: `≤ 1.9·10^{-29}`.
2. **Lemma 2.** Correct. Trivial zeros at `−11/2 − n`, simple; `s = 0`, `s = −1` regular; residues `−a_0 t`,
   `+a_{−1}`; telescoping sum `2/9`; `Σ 1/|ρ|² = 2Λ'/Λ(1)` under GRH needs only (F4), no simplicity; (2.2) holds
   for `t ≥ 1` because `a_{−1} − 2/9 > 0`.
3. **Lemma 3.** Correct. Kernel, residues at `ρ/2`, `0`, `−1/2`, `−1`, `−w/2`, multiplicities
   (`3,5,7,9`: 1; odd `≥ 11`: 2; even `≥ 12`: 1), all signs. `a_{S,0} ≥ 0` is proved correctly under `RH_S(2)`:
   on-line zeros give ratio `≤ 3` (equality at `γ = 0`), off-line zeros have `|γ| > 2 > 2 − β`, so monotonicity
   of `a/(a²+γ²)` applies. `K_S` is a valid upper bound.
4. **Lemmas 4, 5.** Correct: `56g² + 414g ≥ 0`; `a/(a²+u²) ≥ 1/(2+u²)` on `[1,2]`; the convexity step (needs
   `γ − d = H' > 0`, true); `c`; Binet bound and the factor bookkeeping `1/(24 x_j T)`; Fubini for the
   oscillatory part (absolutely convergent); the final use for zeros anywhere in the strip
   (`|ρ(ρ+2)| ≥ γ²`, `t^{1+β/2} ≤ t^{3/2}`). `T(60, 51.45) = 0.18460434` recomputed.
5. **Lemma 6, §3.** Correct: every coefficient of `L_A`, `L_B`, the count of odd `k ≥ 7`, `2c_θ`, the threshold
   526.2, the switch at `10^14` (indicator off at `t = 10^14` exactly, consistent with Lemma 6), the overlap.
   `(9 999 991, 10^7)` is irrelevant: the analytic bound takes over at `10^6`.
6. **§4.** Lemma 7, `h''`, `A_1(T)`, `A_2(T)` re-derived; signs right. Negative patch handled exactly (not
   bounded). `x ≥ 97`, small `x`, `0.15` correct. `1/3`: see R1.
7. **§5.** V0, V1 (normalisation right; `f` is Schwartz), V2 (both exponential factors, range of σ), `E_tail`,
   `E_alias` correct. Lemma V3 is airtight: multiple zeros, a zero at `1/2`, real zeros, off-line pairs
   (`1/(2−β) + 1/(1+β) ≥ 4/3`), interval-valued `γ_j` (lower bound of the subtracted terms in balls),
   `|Rem'| ≤ Rem` (`a ≥ 1`), `e^{1/16}`. Two sign-change intervals cannot hold the same zero (they lie in
   different grid cells and the sign at the shared grid point is certified nonzero). `Z = γ_S L_S` on the
   line with `γ_S ≠ 0`, real and continuous. Zeros up to 72 are needed only to make `Rem` small. See M2.
8. **§6.** Correct: the identity for `P(s)` (and `P_1`), `H` analytic in `Re s > 1/3`, `L(σ, Δ) > 0` on the real
   axis from `Δ(iy) > 0`, Landau's hypotheses, both converses.
9. **Numbers.** All numbers in §1–§5 match the data files except M1, M3, M5, M7 and the range claims of R1.
10. **Hidden hypotheses.** None found. Nothing on ζ beyond unconditional ψ, θ bounds and `ζ(2s) ≠ 0` for
    `Re s > 1/2`. `L(s, Sym²Δ)` is entire and in the Iwaniec–Kowalski class. `RH_S(2)` / `RH_S(60)` are used
    consistently.

---------------------------------------------------------------------------

## Verified by re-derivation

On paper (every line recomputed): Lemma 1; Lemma 2 (Mellin kernel, four residue types, `2/9`, Hadamard
argument, (2.2)); Lemma 3 (as in item 3 above, including the duplication-formula form of `γ_S'/γ_S`);
Lemma 4; Lemma 5 (all six steps and the formula for `T(H,H')`); Lemma 6 (both integrations); §3 (assembly,
`E_odd`, monotonicity, small `t`, `t*`); Lemma 7 (`h'`, `h''`, `A_1`, `A_2`); Theorem 2 for `x ≥ 97` and
`x < 97`; Lemma V0 (a), (b), (c); V1; V2; the two sum-to-integral bounds (`c_0` monotone in β, `M_θ`,
`E_tail`); the recurrence for `a_S(p^k)`; Lemma V3; §6 (both Dirichlet-series identities, Landau).
By hand arithmetic: `L_A(10^6) = 0.21205`, `L_A(10^14+) = 0.54583`, `T(60, 51.45) = 0.18461`, `b' = 7.2360`.

Programs read line by line: `gamma_upper` argument order is right (python-flint 0.9.0:
`z.gamma_upper(s) = Γ(s, z)`; tested: `arb(5).gamma_upper(2.5) = 0.10001 = Γ(2.5, 5)`). No float decides
anything (floats only choose `N`, `u_lo`, `u_hi`, each valid for any value). Every ball comparison that
could silently return False is on the conservative side. All tail bounds are upper bounds.

Numerical checks run (scripts in the scratch directory):

| check | result |
|---|---|
| Rerun `cert_constants.py 10000000 10000000`, `cert_small.py 10000000`, `sym2_phi.py 640 160 4`, `sym2_verify.py 640 60 72 200000`, `bounds.py` on a copy | `constants.json`, `small_range_10000000.json`, `sym2_phi_M640.json`, `sym2_rh_H60.json`, `bounds.out` byte-identical to `data/`; sha256 prefixes of §7 confirmed |
| Rerun `independent_check.gp` | identical to stored output except the missing version line (M6) |
| Lemma 1 at four `u` (PARI, 38 digits) | defects `5.9·10^{-39}` … `1.2·10^{-35}` |
| Identity (2.1) and the four bounds of §3 at `t = 10^3, 10^5+0.5, 10^6` (PARI) | defect `≤ 1.9·10^{-29}`; `ψ_1/t^{3/2} = −0.0064, 0.0243, −0.0239` (bound `−0.118, −0.113, −0.112`); `J_ζ/t^{3/2} = 0.6083, 0.6619, 0.6642` (bound `0.4655, 0.6030, 0.6309`); `J_S/t^{3/2} = −0.0164, −0.0119, −0.0009` (bound `+0.312, +0.224, +0.207`); odd part `0.0228, 0.0302, 0.0076` (bound `1.106, 0.463, 0.284`) |
| Lemma 3 formula against the direct sum, 161 zeros below 150 from PARI | `t = 2`: direct 0, formula `0.0019`; `t = 3.5`: `−0.0009`; `t = 7.3`: `−1.6441` vs `−1.6472`; `t = 20.5`: `−10.2301` vs `−10.2228`; `t = 101.7`: `−83.673` vs `−83.621`; `t = 1000.5`: `−520.184` vs `−520.151` |
| `(L_S'/L_S)(−2)` directly and by the functional equation; `K_S`; `a_{S,0}`; `(12/5)Λ_S'/Λ_S(2)` | `−0.2282526551` both ways; `0.3417553424`; `0.9410336579`; `0.7059630031` |
| True sums (PARI) | `Σ_{|γ|≤60} 2/|ρ(ρ+2)| = 0.28606 ≤ C_S`; `Σ_{|γ|>60} 2/γ² ≈ 0.10409 ≤ T_60` |
| Lemma 5, Binet bound vs `Re γ_S'/γ_S(2+iT)` | `T = 8`: `0.9944 ≤ 1.8283`; `51.45`: `3.18311 ≤ 3.19020`; `60`: `3.40619 ≤ 3.41134`; `200`: `5.19263 ≤ 5.19320` |
| `T(60, 51.45)` (PARI) | `0.184604339…`, `c = 1.0052159` |
| `Rem(T)` from PARI `lfun` and PARI zeros `≤ 72` | `Rem(0) = 0.06835`, `Rem(30) = 0.07976`, `Rem(59.9375) = 0.18411`, `Rem(60) = 0.18486` (program: `≤ 0.068551`, `≤ 0.079962`, max `≤ 0.185063`) |
| Zero counts (PARI `lfunzeros`, not rigorous) | 41 below 60, 55 below 72, 161 below 150 |
| `Λ_S` from the Poisson table vs PARI `lfunlambda` | `t = 0`: `0.0041368748687`; `10`: `−9.8823484·10^{-8}`; `30`: `−2.2927765·10^{-24}`; `50`: `−8.8545406·10^{-41}`; `Λ_S(2) = 0.00516555034690975608302261527257`; `Λ_S(−0.7+3i) = 0.00142882831118 − 0.00133336648751 i`; all digits agree |
| Mellin transform of `K` vs `γ_S` (mpmath quadrature) | `s = 2`: `0.0064098638` both; `0.5`: `0.00818292944567`; `1.3`: `0.00634871463291` |
| Lemma 7 with `T = 97` at `x = 150.5`, `1000` (PARI `intnum`) | `xR(x) = 409.50993367960245328653…` both sides; `5056.29638595761063235…` both sides (38 digits) |
| Gaps `[97,101]`, `[9973,10007]` (R1) | `N ≥ 1.0290`, `G ≥ 0.6393`; `N ≥ 1.0022`, `G ≥ 0.6284` |
| `ψ(y) ≥ y − √2√y`, `θ(y) < y` to `10^7` (numpy floats, independent sieve) | sup over `y ≥ 2` is `1.3319` (at `3^-`); `max(θ(p) − p) = −1.208`; `ψ(10^7) = 9998539.4033`, `θ(10^7) = 9995179.3179` |
| Table radii and `Z` radii (M2) | max table radius `8.3·10^{-64}` at `k = 447`; radius of `Z`: `1.0–2.2·10^{-64}` |
| (F6) against the text of arXiv:2002.11068v2 | Corollary 11.2, (3.39): row `k = 2`, `M_k = 0.2`, `X_0 = 3 594 641` is there |
| `D_0(17)` (Remark 4.2) | `−0.2273 < 0` |

## Not checked

* The originals of Deligne, Hecke, Shimura, Gelbart–Jacquet, Rosser–Schoenfeld; Iwaniec–Kowalski
  Theorem 5.6 and Proposition 5.7 are quoted from memory (the statements used are the standard ones).
* Correctness of FLINT/Arb (incomplete gamma, `acb.integral`, elementary functions) and of PARI.
* I did not write a second rigorous implementation of the table `Φ_k` or of the completeness test; the
  cross-checks above are floating point (PARI) at selected points.
* `selftest.py` (mutation tests, direct functional-equation test) was not rerun; `explore_float.py`,
  `LITERATURE.md`, `RESULTS.md`, `PROGRESS.md` were not reviewed.
* The contour estimates behind Lemmas 2 and 3 were checked for plausibility and through their numerical
  consequences, not written out (M9).
* Remarks 4.2 and 4.3 beyond the numbers quoted from `bounds.out`.
