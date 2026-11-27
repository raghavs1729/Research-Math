# Independent numerical audit of TAU-RACE.md (claims (1)–(6))

Date: 2026-11-26. Auditor: AI agent (Claude), working blind: all scripts below were written and run
**before** opening `code/`; only TAU-RACE.md §0, §1, §5 and the rational intervals in
`data/sym2_rh_H60.json` were read beforehand. `code/` and the other `data/` files were read afterwards
(§C, §D). Nothing in the repository was modified except this file.

**Versions.** The claims were taken from, and the code was first read at, commit `f2a58f1` (10:11). While the
audit was running the repository moved to commit `5d8ab00` (10:50: referee repairs; the Φ table and the zero
certificate were regenerated). I then read the diff `f2a58f1..5d8ab00` of `code/`, `data/` and TAU-RACE.md and
repeated every comparison against the new files. Where the two versions differ this report says so.

Tools: PARI/GP 2.15.4, Python 3.13 with mpmath 1.3.0, python-flint 0.9.0 (exact integers only), numpy 2.5.3.
**Nothing in this audit is interval arithmetic.** Everything is floating point (PARI 19–60 digits, mpmath
30–50 digits, IEEE double / x87 long double for the sums over primes), except τ(p), which is exact.
Agreement of independent floating-point routes is evidence, not proof.

Scratch directory (scripts, logs, zero lists):
`/tmp/claude-0/-home-user-weighted-prime-races/559a4317-d376-57ac-9c95-b35c412a8522/scratchpad/audit/`

## A. Summary

| item | verdict | discrepancy |
|---|---|---|
| (1) σ_0², a_0, a_{−1}, zero sum | PASS | none (30 digits, two independent routes; zero sum agrees to 8·10^{-11}) |
| (2) Λ_S'/Λ_S(2), L_S(1/2), L_S(1), L_S(2) | PASS | none (30 digits, two independent routes + Euler product) |
| (3) 41 zeros ≤ 60, completeness | PASS for the count/completeness; text errata | at `f2a58f1`: quoted first zero `7.0702241402869…` wrong in its last two digits (true 7.07022414028743…), and 10 of 55 intervals wider than the stated 2^{-38}; both repaired in `5d8ab00`. New in `5d8ab00`: §5.1 misstates the size of Z(60), Z(72) by 6 and 9 orders of magnitude (harmless) |
| (4) Θ_1, R at primes ≤ 10^7, t*, x*, min G, A_1, A_2 | PASS | none |
| (5) L_B values, monotonicity, T(60, 51.45), true tail | PASS | none |
| (6) explicit formulas for ψ_1 and J_S | PASS | discrepancies ≤ 2.5·10^{-4} (ψ_1, 3035 zeros) and ≤ 5.5·10^{-3} (J_S, 608 zeros), of the expected size; exact at t = 1 to 10^{-57} |

No disagreement with any certified enclosure in `data/*.json` was found, in either version.

## B. Details

### (1) Constants of L(s, Δ)

Methods (both different from the repository's "theta integral with incomplete gamma in Arb"):

* **Route A** (`a1_delta.gp`): PARI `lfunmf` + `lfun` with derivative, 60 digits (`lfuncheckfeq` = −260 bits).
  `σ_0² = 2(−log 2π + ψ(13/2) + L'/L(1))`. `L'/L(0)`, `L'/L(−1)` both directly and through the functional equation.
* **Route B** (`a1b_delta_mellin.py`, no PARI, no incomplete gamma): mpmath quadrature (50 digits) of
  `∫_1^∞ Δ(iy)[y^{s+9/2} ± y^{11/2−s}](log y) dy` with `Δ(iy)` evaluated from the η-product (`qp`).
* **Zero sum** (`a1_zeros.gp`, `a1c_zero_sum_corrected.py`): PARI `lfunzeros` up to T, plus the tail
  `∫_T^∞ 2/(1/4+t²) dN_0(t)`, `N_0(t) = (1/π) arg[(2π)^{-s}Γ(s+11/2)]_{s=1/2+it}`, plus the boundary term `−2S(T)/(1/4+T²)`.

| quantity | route A (PARI) | route B (mpmath) | claimed |
|---|---|---|---|
| σ_0² | 0.112252140471908526745583426253246256296648… | 0.1122521404719085267455834262532462562966 | 0.112252140471908526745583426253… |
| L'/L(0, Δ) | 0.170657847591640096454240918080193115593… | 0.170657847591640096454240918080193115593 | 0.170657847591640… |
| L'/L(−1, Δ) | 0.281799519641618645635973197303320048659… | 0.2817995196416186456359731973033200486591 | ∈ [0.281796, 0.281803] |
| L(1, Δ) | 0.839345512031942086488146670154135692299… | 0.8393455120319420864881466701541356922992 | (0.8393…) |
| L'/L(1), L'/L(2) | 0.101091806245366805018006162515257…, 0.058326202571456631904649951668199… | | |

Zero sum versus σ_0² (the number of zeros found equals `N_0(T) + S(T)` with `|S(T)| < 0.23` in all three runs):

| T | zeros 0 < γ ≤ T | Σ_{|γ|≤T} 1/|ρ|² | + smooth tail | diff to σ_0² | + boundary term, diff |
|---|---|---|---|---|---|
| 500 | 540 | 0.105404936138 | 0.112250826316 | −1.3·10^{-6} | −3.5·10^{-9} |
| 1000 | 1298 | 0.108387478697 | 0.112251686098 | −4.5·10^{-7} | +9.6·10^{-10} |
| 2000 | 3035 | 0.110099410482 | 0.112252148685 | +8.2·10^{-9} | +8.4·10^{-11} |

Verdict: **PASS**. Caveats: floating point; T = 2000 was run at 19 digits (the 38-digit run exceeded 30 min),
T = 500, 1000 at 38 digits. First zero γ_1 = 9.22237939992110252224.

### (2) Constants of L_S = L(s, Sym²Δ)

* **Route A** (`a2_sym2.gp`): PARI `lfuncreate` from `direuler` with the motivic local factor, gamma shifts
  `[-10,0,1]`, weight 23, N = 1, ε = 1; 57 digits. `lfuncheckfeq` = −109 bits at 57 digits (−67 bits at 38 digits):
  PARI's own check is weaker than the working precision; the values are nevertheless confirmed to 30 digits by route B.
  (`lfunsympow` exists in this PARI but answers "not yet implemented".)
* **Route B** (`a2b_sym2_afe.py`, no PARI, no theta table): smoothed approximate functional equation,
  `Λ_S(s) = Σ_{n≤80} a_S(n)[F(s,n) + F(1−s,n)]`, `F(s,x) = (1/2πi)∫_{Re w=3} γ_S(s+w) x^{-(s+w)} dw/w`, Mellin–Barnes
  integrals by mpmath quadrature (34 digits, |Im w| ≤ 76), `a_S(n)` from my exact τ.
* **Route C** (`a3_rem.py`): Euler product over p ≤ 10^7 with exact τ(p) (double precision).

| quantity | route A | route B | claimed |
|---|---|---|---|
| L_S(1/2) | 0.505549375222713515001575353290033064… | 0.50554937522271351500157535329 | 0.505549375222713515… |
| L_S(1) | 0.631792945727883203011116326434919338… | 0.631792945727883203011116326435 | 0.631792945727883203… |
| L_S(2) | 0.805875209448668759610154942699503675… | 0.8058752094486687596101549427 | 0.805875209448668759… |
| Λ_S'/Λ_S(2) | 0.294151251276053089971061019233191… | 0.294151251276053089971061019233 | ∈ [0.2941460, 0.2941565] |
| L_S'/L_S(2) | 0.160153260311665050252223746092… | route C: 0.16015326031800 (tail ≤ 3·10^{-7}) | |
| Λ_S'/Λ_S(1), Λ_S'/Λ_S(3) | 0.0992107753484266096056350582737…, 0.479608199517955135600428233816… | same 30 digits | |
| L_S'/L_S(0), L_S'/L_S(−2) | 0.941033657929609593099770092828, −0.228252655128807821783911971603 | (direct = via functional equation) | K_S ≤ 0.34176 |

`C_S = (12/5)Λ_S'/Λ_S(2) = 0.7059630030… ≤ 0.70598`. Verdict: **PASS**.
Remark: the true value of `(L_S'/L_S)(−2)` is −0.22825; `K_S = 0.34176` is a valid upper bound for it and for its
absolute value (it comes from `|L_S'/L_S(3)| ≤ −3ζ'/ζ(3)`; the true `L_S'/L_S(3)` is +0.07554).

### (3) Zeros of L_S up to 60 (and 72), completeness

* Zeros (`a3_sym2_zeros.gp`): PARI `lfunzeros`, 38 digits, heights 72, 200, 400: **41** zeros in (0, 60],
  **55** in (0, 72], 240 in (0, 200], 608 in (0, 400]. First zero **7.0702241402874349921024418838…**;
  γ_41 = 58.962868339014035…, γ_42 = 60.621217659536835…
* Comparison with `data/sym2_rh_H60.json` (`a3_compare_intervals.py`, exact rational comparison): all 55 of my
  zeros lie inside the corresponding stored interval; the intervals are disjoint and increasing; 41 have right
  end ≤ 60, 14 have left end > 60.
* Completeness test 1, argument principle (`a3_sym2_zeros.gp`): `N(T) = (1/π)[arg γ_S(1/2+iT) + arg L_S(1/2+iT)]`,
  the second argument followed continuously from its principal value at `2 + iT` (where `|log L_S| ≤ 3 log ζ(2) < π/2`)
  along the horizontal segment in 4000 and in 9000 steps (results agree to 10^{-35}):

  | T | 10 | 20 | 30 | 40 | 50 | 59.9 | 60 | 60.5 | 72 | 199.63 | 399.63 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | N(T) (integer to 10^{-33}) | 1 | 6 | 13 | 22 | 31 | 41 | 41 | 41 | 55 | 240 | 607 |
  | S(T) | −0.416 | −0.214 | −0.066 | 0.597 | 0.113 | −0.193 | −0.301 | −0.844 | 0.179 | 0.128 | −0.134 |

  equal to the number of zeros found in every case.
* Completeness test 2, the document's identity (`a3_rem.py`): `Re Λ_S'/Λ_S(2+iT)` from the Euler product (p ≤ 10^7,
  my τ table) and mpmath digamma; `Rem(T)` = that minus the 55 zeros' terms; then the zeros in (72, 400] and the
  smooth density beyond 400:

  | T | Re Λ_S'/Λ_S(2+iT) | Rem(T) (mine) | claimed | zeros in (72,400] | density tail > 400 | residual |
  |---|---|---|---|---|---|---|
  | 0 | 0.294151251282 | 0.0683503 | 0.0686 (repo upper bound 0.068551) | 0.0498976 | 0.0184552 | −2.5·10^{-6} |
  | 10 | 1.29202796106 | 0.0694285 | (≤ 0.069629) | 0.0509658 | 0.0184652 | −2.5·10^{-6} |
  | 30 | 2.44750070754 | 0.0797615 | (≤ 0.079962) | 0.0612179 | 0.0185460 | −2.5·10^{-6} |
  | 60 | 3.13376644903 | 0.1848615 | ≤ 0.185063 (certificate of `5d8ab00`: ≤ 0.185062) | 0.1660394 | 0.0188247 | −2.6·10^{-6} |

  The residual has the size of the boundary term `2w(400)·S(400)`. Maximum of `Rem(j/8)`, 0 ≤ j ≤ 480: 0.184862 at T = 60
  (claimed bound 0.185063; the difference 2.0·10^{-4} is the repository's Euler-tail allowance for P = 2·10^5).
  My own deletion test (zero #1, #8, #21, #41 removed, grid 1/8): max Rem = 0.7420, 0.7407, 0.7560, 0.8400 (> 2/3);
  the document's 0.743, 0.741, 0.758, 0.842 are upper bounds with P = 2·10^4, grid 1/4 (allowance 0.0017): consistent.
* Completeness test 3, sum rule (`a3c_sum_rule.py`): `Σ_ρ 1/|ρ|² = 2Λ_S'/Λ_S(1) = 0.198421550696853…`; zeros ≤ H + smooth
  tail + boundary term minus this: +1.2·10^{-6} (H = 72), −1.3·10^{-7} (H = 200), −9.8·10^{-9} (H = 400). A missing zero
  pair at height g would cost 2/(1/4+g²) ≥ 1.25·10^{-5} for g ≤ 400.

Verdict: **PASS** for "exactly 41 zeros with 0 < Im ≤ 60 on the line, list complete" (floating point; test 1 is not
rigorous in its step size, tests 2–3 use a smooth-density tail). **Errata found:**

1. (`f2a58f1`; repaired in `5d8ab00`.) Theorem 4 and my task statement quote the first zero as `7.0702241402869…`.
   The zero is 7.07022414028743…; `7.07022414028688` is the *left endpoint* of the stored isolating interval
   `[7.0702241402868822, 7.0702241402905202]`. `5d8ab00` prints `7.0702241402…`, which is correct.
2. (`f2a58f1`; repaired in `5d8ab00`.) §5.2 said every sign change is bisected to length 2^{-38}; in the certificate of
   `f2a58f1` (sha256 `e08b0806…`) zeros #46–#55 (all above 60) had widths 2^{-31}, 2^{-31}, 2^{-30}, 2^{-27}, 2^{-19},
   2^{-22}, 2^{-16}, 2^{-13}, 2^{-11}, 2^{-11}. Cause: §D.1. In the certificate of `5d8ab00` (sha256 `1e6a096b…`) all 55
   intervals have width 2^{-38}; all 55 still contain my zeros (re-checked).
3. (New in `5d8ab00`, not repaired.) The sentence added to §5.1, "`|Z(t)|` of size `10^{-44}` at `t = 60` and `10^{-52}` at
   `t = 72` away from zeros", is wrong: `|Z(60)| = 9.32·10^{-50}`, `|Z(59.5)| = 2.58·10^{-49}`, `|Z(71)| = 1.67·10^{-61}`,
   `|Z(72)| = 2.06·10^{-61}` (`a9_Z_size.gp`: PARI `lfunlambda`, equal to `|γ_S L_S|`; the same six digits come out of the
   repository's own Φ table evaluated in floating point). The conclusion drawn there is unaffected: the radius of
   the enclosure of Z is `≈ 10^{-83}`, still 22 orders of magnitude below `|Z(72)|`.
4. (New in `5d8ab00`, cosmetic.) The regenerated certificate has `rem_max = 0.185062`, `rem_bound = 0.196998`;
   §5.3 still prints `0.185063`, `0.196999` (valid, since larger; but "max Rem unchanged" in §10 is off by one digit).

### (4) Small range with exact τ(p)

* τ (`a4_tau_table.py`): `Δ = q(Σ(−1)^n(2n+1)q^{n(n+1)/2})^8` mod `q^{10^7}` with FLINT `fmpz_poly.mul_low` (21 s). Checks:
  known values, multiplicativity, `τ(p²) = τ(p)² − p^{11}`, `τ(p) ≡ 1 + p^{11}` (mod 691) at **all** 664 579 primes,
  Deligne at all primes; equality with PARI `ramanujantau` (class-number formula) at 607 primes: all p ≤ 1999,
  300 random p ≤ 10^7, and 1223, 9132733, 9132889, 9999991 (`a4_crosscheck_tau.py`): 0 mismatches.
  (This is the same construction as the repository's `tau_table.py`, written independently; the PARI comparison is
  the independent leg.)
* Sums (`a4_checks.py`): cumulative sums in x87 long double (64-bit mantissa) and in double; identical to the 12
  digits printed. Closed forms in mpmath, 50 digits.

| claim | my value | verdict |
|---|---|---|
| Θ_1(p) > 0 at all primes 5 ≤ p ≤ 10^7 | true; minimum 0.212761797245 at p = 5 | ok |
| R(p) > 0 at all primes 7 ≤ p ≤ 10^7 | true; minimum 0.160815812419 at p = 7 | ok |
| t* = 4.266791956347541… | 4.26679195634754108927818019391 | ok |
| x* = 5.518054631888274… | 5.51805463188827498152859434296 | ok |
| min_{p ∈ [10^4,10^7]} Θ_1(p)/p^{3/2} ≈ 0.625775 at 9132733 | 0.625775035267 at p = 9 132 733 | ok |
| A_1(97) = 0.671972273685… | 0.6719722736852640993407984 | ok |
| A_2(97) = 10.326446649… | 10.32644664916132390886086 | ok |

Also: Θ_1(3) = −0.3675968, Θ_1(5) = 0.2127618, R(5) = −0.0787046, R(7) = 0.1608158; final sums at 9 999 991:
Θ = 3801.734953658, Θ_1 = 20958813966.4690, D_0 = 271.623191660, D_1 = 1087226617.6091 (inside the repository's
enclosures). Over real arguments (`a4_real_minima.py`, float): inf G = 0.611712 on [100, 9999991] (at t = 1223) and
0.625775 on [10^4, 9999991]; min N = 0.174700 on [11, …] (at x = 11) and 0.795108 on [100, …] (at 9132889). The
document's 0.5814, 0.6257, 0.1718, 0.7951 are valid lower bounds; the first and third are not sharp (see §D.5).
Verdict: **PASS**.

### (5) L_B, monotonicity, T(60, 51.45), true tail (`a5_LB.py`, mpmath 40 digits + numpy scan)

With σ_0² = 0.11225214047…, a_0 = 0.170658, C_S = 0.70598, K_S = 0.34176, T_60 = 0.184605, ε_1 = 7.69843661·10^{-4}:

| t | L_A(t) | L_B(t) | claimed |
|---|---|---|---|
| 10^6 | 0.212045947 | 0.027440947 | ≥ 0.2120, ≥ 0.0274 |
| 10^7 | 0.348398586 | 0.163793586 | ≥ 0.3483, ≥ 0.1637 |
| 10^14 (without ε_1) | 0.546345115 | 0.361740115 | |
| 10^14 + (with ε_1) | 0.545831885 | 0.361226885 | ≥ 0.5458, ≥ 0.3612 |
| ∞ | 0.553901297 | 0.369296297 | ≥ 0.5539, ≥ 0.3692 |

* Monotonicity: exact derivative of L_B (checked against numerical differentiation) scanned at 2·10^6 points, uniform
  in log t on [527, 10^14]: `L_B' > 0` everywhere, `min t^{5/4} L_B'(t) = 1.283` (at t = 527); the sampled values are
  increasing. (In fact `L_B' > 0` on the whole grid of (32, 527) as well; 527 is where the single term
  `t^{-5/14}(log_2 t − 5)` starts to decrease, `exp(5 log 2 + 14/5) = 526.23`.) L_B = 0 at t = 7.051·10^5, L_A = 0 at t = 1.111·10^5.
* `T(60, 51.45) = 0.184604339768 ≤ 0.184605`. (Minimum over H' on a 0.01 grid: 0.184603086 at H' = 51.51.)
* True tail: `Σ_ρ 1/|ρ|² = 2Λ_S'/Λ_S(1) = 0.198421550697`; the 41 zeros and their conjugates give 0.146376831431; hence
  `Σ_{|γ|>60} 1/|ρ|² = 0.0520447193` and `Σ_{|γ|>60} 2/|ρ|² = 0.104089439`. Also
  `Σ_ρ 1/(ρ(ρ+2)) = (Λ_S'/Λ_S(1) − Λ_S'/Λ_S(3))/2 = −0.190198712` gives `|Σ_{|γ|>60} 2 Re 1/(ρ(ρ+2))| = 0.104031481`. So
  (zeros above 60 assumed on the line) `Σ_{|γ|>60} 2/|ρ(ρ+2)| ∈ [0.10403, 0.10409]`, a factor 1.77 below `T_60 = 0.184605`.

Verdict: **PASS**. Caveat: the monotonicity scan is a grid scan in double precision, not a proof (the paper argument
in §3 of TAU-RACE.md is elementary and correct: every other t-dependent term is monotone for all t > 0).

### (6) Explicit formulas (`a6_explicit.py`, `a6b_identities_t1.py`; mpmath 30 digits, my PARI zeros)

ψ_1(t) against `−Σ_ρ t^{ρ+1}/(ρ(ρ+1)) − a_0 t + a_{−1} − Σ_n t^{−9/2−n}/((n+11/2)(n+9/2))`, zeros of L(s,Δ) with |γ| ≤ T:

| t | ψ_1(t) direct | lhs − rhs, T = 500 | T = 1000 | T = 2000 (3035 zeros) | absolute bound t^{3/2}Σ_{|γ|>2000}|ρ|^{-2} | expected size |
|---|---|---|---|---|---|---|
| 10 | −1.92486552599732 | −5.2·10^{-4} | +1.3·10^{-4} | +2.9·10^{-5} | 0.068 | 3·10^{-5} |
| 30 | −4.13885560075730 | −4.0·10^{-3} | +1.9·10^{-4} | −2.5·10^{-4} | 0.354 | 1·10^{-4} |
| 100 | −14.3196460979777 | −6.0·10^{-2} | +1.0·10^{-3} | +4.1·10^{-5} | 2.153 | 7·10^{-4} |

J_S(t) against `−Σ_ρ 2t^{1+ρ/2}/(ρ(ρ+2)) − (L_S'/L_S)(0)t + 2√t + (L_S'/L_S)(−2) − (trivial zeros)`, zeros of L_S with |γ| ≤ T:

| t | J_S(t) direct | lhs − rhs, T = 60 | T = 100 | T = 200 | T = 400 (608 zeros) | absolute bound 2t^{5/4}Σ_{|γ|>400}|ρ|^{-2} | expected size |
|---|---|---|---|---|---|---|---|
| 10 | −3.69397684167979 | +2.5·10^{-4} | −3.3·10^{-3} | +3.6·10^{-3} | −1.4·10^{-3} | 0.44 | 1·10^{-3} |
| 30 | −12.9932767844395 | −4.7·10^{-2} | −1.6·10^{-2} | +1.1·10^{-3} | −1.1·10^{-3} | 1.73 | 4·10^{-3} |
| 100 | −84.0986025596883 | −4.1·10^{-2} | −5.0·10^{-2} | −1.4·10^{-2} | +5.5·10^{-3} | 7.78 | 1.5·10^{-2} |

"Absolute bound" uses the exact tails `σ_0² − Σ_{|γ|≤T}|ρ|^{-2}` resp. `2Λ_S'/Λ_S(1) − Σ_{|γ|≤T}|ρ|^{-2}`. "Expected size" is the
heuristic `≈ 2 t^{3/2} T^{-2}(1 + (density)/log t)` (resp. `4 t^{5/4} T^{-2}(…)`) for an oscillating tail. All discrepancies are
of the expected size, change sign, and decrease as T grows. At t = 1, where both zero sums have closed forms
(`Λ'/Λ(1) − Λ'/Λ(2)` and `Λ_S'/Λ_S(1) − Λ_S'/Λ_S(3)`) and the trivial sums telescope (2/9 and 8/9 + 2/9 + 1/10), the two
formulas hold to 4·10^{-61} and 1.4·10^{-57}. This confirms every residue (signs of `a_0 t`, `a_{−1}`, `+2√t`, the
multiplicities of the trivial zeros) to the full precision of the constants. Verdict: **PASS**.

## C. Comparison with the repository after unblinding

Every number of mine was compared with `data/constants.json`, `data/small_range_10000000.json`, `data/sym2_rh_H60.json`
(both versions), `data/bounds`, `data/sym2_verify_H60`, `data/selftest`, `data/independent_check` (`.out` at `f2a58f1`,
`.txt` at `5d8ab00`). `constants.json` and `small_range_10000000.json` are byte-identical in the two commits.

* All 17 numerical enclosures in `constants.json` contain my values (`a8_check_enclosures.py`; σ_0², a_0, L(1,Δ),
  L'/L(1), Λ(1): to the last of 30 digits; `LS_logderiv_2`, `LambdaS_logderiv_2`, `a_minus1`, `L_logderiv_2`: my value
  sits near the middle of each enclosure; `K_S`, `gammaS_logderiv_0/2`: 14 digits). `c_psi = √2` (sup attained as
  y → 2^-) and `θ(y) < y` up to 10^7 also hold in my float check; ψ(10^7) = 9998539.4033459, θ(10^7) = 9995179.3178562.
* `sym2_rh_H60.json`: `LS_half`, `LS_1`, `LS_2` enclosures (30 digits) contain my values; zeros: see (3).
* The repository's Φ table (`5d8ab00`), evaluated in floating point through its Poisson formula, gives
  Z(0), Z(30), Z(59.5), Z(60), Z(60.2), Z(66), Z(71), Z(71.8), Z(72) equal to PARI's `lfunlambda` to the 6 digits
  compared (signs included) — an independent check of the table itself.
* `small_range_10000000.json`: `t_star`, `x_star`, `A1`, `A2` at 97, the values at p = 5, 7 and the final sums agree in
  every digit; `infG`, `infN` are below my floating minima as they must be.
* `bounds.out`/`bounds.txt`: my L_A values agree to the 6 digits printed. `bounds.py` uses the exact `T(60, 51.45) = 0.1846043…`
  rather than 0.184605, so its `L_B(10^6) ≥ 0.027441` versus my 0.0274409 (with 0.184605): no conflict.
* `independent_check.out`/`.txt` (the repository's own PARI run) agrees with my PARI numbers digit for digit, as expected:
  it uses the same PARI constructions that my task statement suggested. The genuinely independent legs of this
  audit are therefore routes B (mpmath Mellin quadrature of the η-product; mpmath Mellin–Barnes approximate
  functional equation), the argument-principle count, the sum rules, the Euler product at P = 10^7 with my own τ
  table, and the explicit-formula tests with 3035 / 608 zeros (the repository used 277 zeros of L(s,Δ), no J_S test
  away from t = 1).

**Digit-level disagreements:** only the errata of item (3): at `f2a58f1` the quoted digits `…2869` of the first zero
and the interval widths (both repaired); at `5d8ab00` the orders of magnitude of `|Z(60)|`, `|Z(72)|` in §5.1 and the
last digit of `max Rem`. Everything else agrees in every digit compared.

## D. Places where the repository code or text is non-rigorous, loose, or wrong

None of the points below invalidates a certified statement. "Open" = still present at `5d8ab00`.

1. **Repaired in `5d8ab00` (found independently here and by the referee).** `sym2_verify.py` stops bisecting silently
   when the sign of Z cannot be decided (`if sm == 0: break`). At `f2a58f1` this happened for zeros #46–#55 because
   the Φ_k balls were far wider than §5.1 suggested: radii up to 8.3·10^{-64} (the n-tail bound `∫_N^∞` with N = 1,
   valid but loose), giving Z an uncertainty of 2.2·10^{-64} against `|Z| ≈ 2·10^{-61}` near t = 72. The new n-tail
   ("first omitted term + integral from N+1"; I checked it: the summand `n³e^{−βn^{2/3}}` decreases for `βn^{2/3} ≥ 9/2`,
   which the code asserts) gives radii ≤ 2.91·10^{-84} and an uncertainty of Z of 1.2·10^{-83}; I confirm both
   numbers from the stored table. The silent `break` is still in the code; the certificate now records
   `max_interval_width`, which makes it visible.
2. **Repaired in `5d8ab00`.** Quoted first zero `7.0702241402869…` (left end of the interval, not the zero).
3. **Open (text).** §5.1: `|Z(60)| ≈ 10^{-44}`, `|Z(72)| ≈ 10^{-52}` should read `≈ 10^{-49}`, `≈ 2·10^{-61}` (item (3), erratum 3).
4. **Open (cosmetic).** §5.3 prints `max Rem ≤ 0.185063`, bound `0.196999`; the current certificate says `0.185062`, `0.196998`.
   "Rem(0) is 0.0686" is the certified upper bound 0.068551 rounded; the value is 0.068350. The mutation figures
   0.743, 0.741, 0.758, 0.842 are upper bounds computed with P = 2·10^4 on a grid of 1/4 (`selftest.py`), not with
   the parameters of the certificate; my float values on the grid 1/8 are 0.7420, 0.7407, 0.7560, 0.8400.
5. **Repaired in `5d8ab00`; I missed it on first reading.** `cert_small.py` at `f2a58f1` took the gaps `[q, p]` with
   `q ≥ a`, which leaves `[a, next prime)` uncovered when a is not prime (a = 100, 1000, 10^4, …). Now `p > a`.
   The certificate is byte-identical (the extra gaps are not minimisers); my float minima over real arguments,
   with the partial gap at a included, are those quoted under (4). Still true in both versions:
   `min(Θ_1(q), Θ_1(p))/p^{3/2}` is a lower bound for G on `[q, p]` only because both endpoint values are positive
   (guaranteed by check (A), not asserted at that line), and two of the constants called "numerical truth" are
   certified lower bounds that are not sharp: `G ≥ 0.5814` on [100, N] (true infimum 0.6117), `N ≥ 0.1718` on [11, N]
   (true minimum 0.1747).
6. **Open (not checked by code).** Monotonicity of L_A, L_B on [527, 10^14] is a paper argument; `bounds.py` only
   prints the sentence. The argument is correct and my scan confirms it.
7. **Open (self-test only).** `selftest.py` (F): for y < 1 the n-tail of Φ is omitted (comment only). M3 also records
   that shifting a zero by 0.3 does *not* raise max Rem (0.1865 before and after): Rem detects missing zeros, not
   misplaced ones; the sign-change step is what pins the locations.
8. **Remark.** `K_S` is described as "an upper bound for (L_S'/L_S)(−2)"; the true value is negative (−0.22825). The
   bound holds for the absolute value too, so whichever sign the proof needs is covered.
9. Checked and found sound (both versions): the tail bounds in `cert_constants.py` (C1 geometric tail, Y-truncation,
   C2 Euler tails), the K-kernel tails in `sym2_phi.py`, `E_tail`, `E_alias`, `M_θ`, the Lipschitz factor `e^{1/16}`,
   the rigorous running maximum of Rem, and the outward rounding in `arbutil.py`. The `acb.integral` callback in
   `sym2_phi.K_arb` ignores the `analytic` flag, which is legitimate there (integrand built from `exp`, integer
   powers and division only); `cert_constants.py` passes the flag to `log` as it must.
10. Not audited: the proofs (Lemmas 0–7, V0–V3), the literature inputs (F1)–(F6), and Arb/FLINT/PARI themselves.

## E. Scripts (sha256)

All in the scratch directory named at the top. `.out`/`.log` files next to them hold the outputs quoted here
(`a7_inspect_repo_data_f2a58f1.out` and `…_5d8ab00.out` for the two versions of the Φ table and zero certificate;
the `f2a58f1` JSON files were extracted with `git show` into `repo_f2a58f1/` in the scratch directory).
`a1_zeros.gp` and `a3_sym2_zeros.gp` are templates: `sed "s/TT/1000/" a1_zeros.gp | gp -q` (T = 500, 1000 at 38 digits;
T = 2000 with `realprecision` changed to 19 and output prefix `delta_zeros_p19_`), `sed "s/HMAX/400/" a3_sym2_zeros.gp | gp -q`
(HMAX = 72, 200, 400).

    a7ea9740262d7ba3b1894cb58fe7c9d88d03feef8ee18fcce7db7490e99f6070  a1_delta.gp
    5b0c77d9361d67828fddf1f417dca735c34f9ff46f655a03a827c7a7fff0011c  a1_zeros.gp
    3d7f37e9cc1dacb80d06cc753f9a57cb5bc21f04352262f465412ea6cd5dccc8  a1b_delta_mellin.py
    d82508eade6db934e5eb7b9f88e35a7c6a72ace1f724d1cc7b4f7ed3406e4094  a1c_zero_sum_corrected.py
    e8a0f95763156a16156fec79575fe0fc8a6864a6fd23e69d0d9b085fb385fd05  a2_sym2.gp
    5f6a4852cd4f6da7e7a059819d3c2532a7b3950100e02b145dd45eb018f5f98b  a2b_sym2_afe.py
    dbd3a50c71948c218fab18008ff1da08c4a74bc5831641e0e35c4dd71bbd0046  a3_sym2_zeros.gp
    7e374555beea3535196d1a53c1b57d1f4d0fec014c9a7d1e2feadb12495ca642  a3_compare_intervals.py
    cc91d316bd2f3b27c1a33dd26e9af7ff2578cff21e341b5fe4a777d4ef2cb6be  a3_rem.py
    0c6d3bc14bdd069806822fba3810222e4fb034ed83fbff622df2f5d9a6386aaf  a3c_sum_rule.py
    911d870903ae30b934f6462da94647ff74e7ee3a080e474e45ce8ff70becc102  a4_tau_table.py
    573f04857d34b88ff1772d163a707fe7ed94332e651e3940113f9f3b35ddf331  a4_checks.py
    fe1d5353833281f6deb1f4b719f27f7ec45ccdcf1ae67c3f02f30aa1b268accd  a4_crosscheck_tau.py
    a66bea563b1986d043407c0ba73c4f503df62d49699b8479017ee8d41f34b5a8  a4_crosscheck_tau.gp
    4ab0e703293d732fca060c2c654a9fc9089a2336525a93d1c63ee9353dd8ab18  a4_real_minima.py
    77dd6cf11c1003840ed7614edf1b460b86b39b0e9dd975de7c0efb0b27a6a3ba  a5_LB.py
    512e753f37a6e71bf2e4dabbed9f82d6038352eab2cb0ac9548d2a29f4a70cd5  a6_explicit.py
    dcf8b2c58c428d1e91de58e49fdb432c92a92f9174feb4587808279fe362534a  a6b_identities_t1.py
    be4a99133a71f5c8e79800317e709904a7e1e76deb398d0c3a88f58def93e624  a7_inspect_repo_data.py
    7717cb26fbe2a2e82d1940f592ade3aaef6480e665b48f57081629c7c74b186c  a8_check_enclosures.py
    827c7868221a4bda1d4cb6d7dd366fb5a1d955b8d62a7c25642e39b180b246f7  a9_Z_size.gp

Run times (4 cores): τ table 21 s; `a4_checks.py` 1 min; Δ zeros to 1000: 1 min (38 digits), to 2000: 12 min (19 digits);
Sym² zeros to 72 / 200 / 400: 23 s / 50 s / 2.3 min; `a2b_sym2_afe.py` 5.5 min; `a4_crosscheck_tau.gp` 1.6 min; the rest seconds.
