# rm-open — results (Ramanujan Machine conjectures and open identifications)

**Summary (2026-11-26).**
1. **U1 — PROVED, twice audited.** The continued fraction `PCF[8(n^5+(n+1)^5) − 15(n^3+(n+1)^3) + 9(2n+1), −64n^10]`,
   listed only numerically ("1.20426...", "did not fit any combination of zetas") in Table 1, §3.6 of
   arXiv:2308.11829v2 (Elimelech et al., PNAS 2024), equals **96/(372 ζ(5) − π^5)**. Proof: `U1/THEOREM.md`
   (θ_3^8, Eichler integral of a weight-6 Eisenstein series, Mellin transform). 300–498 digit agreement.
2. **U2 — PROVED, audited twice (own audit + parent referee: "proof complete"); L-value form now unconditional.** Row 5 of the same table
   ("2.45174...", `a_n = 8σ_5 − 12σ_3 + 7(2n+1)`): `1/PCF = (93/4) ζ(5) + (π^6/64) S_0/S_1`,
   `S_0 = Σ C(2k,k)^4/256^k`, `S_1 = Σ C(2k,k)^4/256^k (H_{2k} − H_k − log 2)`; equivalently
   **PCF = 8 L(f,1)/(186 ζ(5) L(f,1) − π^4 L(f,2))**, `f = η(2τ)^4 η(4τ)^4` (299 digits). Proof: `U2/THEOREM.md`.
   **Z — Zudilin's relations (6) (SIGMA 2018) for `(r,t) = (1/2,1/2)`: PROVED, audited; novelty essentially nil.**
   `L(f,1)/F_1(1) = −1/16`, `L(f,2)/F_2(1) = 1/32`, `L(f,3)/F_3(1) = −3/256` (i.e. `L(f,1) = −S_1/2`,
   `L(f,2) = π^2 S_0/16`). Every ingredient is published (Borwein–Straub–Wan 2013; Rogers–Wan–Zucker 2015; Zudilin's
   own Remarks 2, 5); only the assembly with the constants was not found stated. `Z/THEOREM.md`, `Z/AUDIT.md`.
   **Z2, Z3 — Zudilin's relations (6) for `(1/2,1/4)` and `(1/2,1/3)`: PROVED (2026-11-26), refereed** (`Z/REFEREE-2-Z2Z3.md`: "proof complete" for both cases, minor fixes; ratios
   reproduced independently to 52–54 digits; an earlier audit had been stopped unfinished).
   `(1/2,1/4)`, `f = η_4^16/(η_2^4η_8^4)` (16.4.a.a): `L(f,k)/F_k(1) = −1/8, 1/32, −3/448`. `(1/2,1/3)`,
   `f = η_6^14/(η_2^3η_18^3) − 3η_2^3η_6^2η_18^3` (36.4.a.a): `−1/3, 1/18, −1/90`. Method: for `r = 1/2` the first three
   Frobenius coefficients are `(1/π)∫_0^1 U_j(x) dx/√(x(1−x))` with `U_j` the Frobenius solutions of the symmetric-square
   operator of `3F2(1/2,t,1−t;1,1;x)`; for `t = 1/4, 1/3` these are `(1, y, y^2)` times a weight-2 form in a
   Fricke-symmetric Hauptmodul, `F_1, F_2` become complete Mellin integrals of a twist of `f` (level 8, resp. 12),
   and Manin's trick converts to `L(f,k)`; `F_3 = (csc^2 πt + 1/3)π^2 F_1` from a bilateral sum. Literature survey
   (`Z/LITERATURE-13.md`): no relation with `F_k(1)` was found proved for any of the 13 cases other than `(1/2,1/2)`;
   for `(1/2,1/4)` the L-values were known (Rogers–Wan–Zucker), for `(1/2,1/3)` everything was numerical. The other
   11 cases are **not** reached by this method (in 8 of them `L(f,2) = 0` and the `k = 2` relation is trivial; `k = 3`
   is equivalent to `k = 1` in all). `Z/THEOREM-2.md`.
3. **U3 — PROVED (2026-11-26), audited ("complete after listed fixes", fixes applied).** Row 6 ("22.8410..."):
   `1/PCF[8σ5 + 20σ3 − 5(2n+1), −64n^10] = (93/4)ζ(5) − 4 + (π^6/64) S_0 K/(1 + S_1 K)`, `K = 4 + π^2(Z_0 − 2S_0)`, with one
   further hypergeometric constant `Z_0 = Σ C(2k,k)^4/(256^k (k+1)) = 1.0433551893...` (quasi-period-type; no
   relation to `L(f,1), L(f,2), π` found by PSLQ, nothing proved about independence). Proof: `U3/THEOREM.md`
   (gauge map U2 → U3; Legendre-type relation from the Casoratian; bilateral sum for `Z` and the contiguity
   relation `Σ A_k/(k+1) = 2Σ A_k/(1−2k) − 8/π^2`). 48-digit agreement; auditor's independent evaluation 35 digits.
4. **Factorial reduction (arXiv:2111.04468, Arnold Math. J. 2024) — PROVED, refereed by the parent session
   ("complete modulo listed gaps"; gaps closed):** Conjecture 1.3 (the Euler family has FR, with the exact value of
   `GCD_n`) **for `B, m ≠ 0` and `x_1..x_d` not positive integers**; outside that range the conjecture **as printed is
   false** (`d = 2, B = −3/4, m = 1, x = (−2,1,−2,2)`: `p_n = q_n = 0` for `n ≥ 1`). Consequences: Conjecture 1.2 and the
   "if" half of 1.1 only for even-degree `b` with at least `d` roots that are not positive integers; the "if" half
   of Conjecture 1.4 (not both roots positive integers). Prior partial result: O. David, arXiv:2308.02567v1,
   Example 35 (`ζ(d)` case). `FR/THEOREM.md`.
   **"Only if" direction (2026-11-26; refereed, `FR/AUDIT-2.md`, fixes applied):** for degree (1,2), a prime `p ∤ ω`
   modulo which `A^2+4B` is not a nonzero square never divides two consecutive `GCD_n` (Theorem C); hence the paper's
   two-step gcd `GCD2_n = gcd(p_n,q_n,p_{n−1},q_{n−1})` has factorial size **iff** `a_n` is in the families of
   Conjecture 1.4 (with `k = 0` allowed in (Eq. 6)), is `≤ C^n √(n!)` otherwise, and also for every `a_n` when `b` is an
   irreducible quadratic (Theorems D, E). So Conjecture 1.4 holds with `GCD2` in place of `GCD`, and the paper's
   Conjecture 3 (`GCD2_n ≐ GCD_n`) in degree (1,2) would imply the "only if" halves; these remain **open** for the
   `limsup` definition with `GCD_n` (transient common factors not excluded; exact search of 1.39M fractions: no
   counterexample). Printed statements that are false as stated: (Eq. 6) excludes `k = 0`, but
   `PCF[n, 2(n+1)(n+2)]` has FR (non-terminating; misprint-level); terminating `PCF[n+1, n(n−1)]`,
   `PCF[n^2+1, (n−1)^2(n^2−2)]` have FR ("only if" of 1.4, 1.1); `PCF[n−3, 2(n−1)(n−3)]` has `p_n = q_n = 0` ("if" of 1.4).
5. **Status map:** every itemised conjecture on the project's website lists (ζ(3), Catalan, π², mixed zeta, log 2)
   is already proved in the literature; see `STATUS-MAP.md`.
6. **Families of arXiv:2412.16818 (NeurIPS 2024) — PROVED (low standing, not independently audited):** the two
   "awaiting proof" families of §4.3 (`a_n = i+2j+1, b_n = n^2+(i+k)n`: explicit Möbius transforms of log 2 / π;
   `a_n = −2n+j−1, b_n = −n^2+jn+k`: the smaller root) and the Gauss/lemniscate formulas (8), (9). They are
   corollaries of the classical Stieltjes–Rogers fraction and of Euler's fraction + Gauss's theorem. `SR/THEOREM.md`.
   **Dropped:** FrontierMath part, on instruction (`epoch-small-diophantine/` holds only an interrupted search).
7. Nothing here is formalized in Lean (modular forms / Mellin transforms); exact computer-algebra and
   rational-arithmetic checks accompany every algebraic step, high-precision numerics every analytic step.

| id | statement | status | files |
|---|---|---|---|
| U1 | `PCF[8σ5 − 15σ3 + 9(2n+1), −64n^10] = 96/(372ζ(5) − π^5)` | proved; audits `U1/audit/`, `U1/AUDIT-2.md` | `U1/` |
| U2 | `1/PCF[8σ5 − 12σ3 + 7(2n+1), −64n^10] = (93/4)ζ(5) + (π^6/64)S_0/S_1` | proved; audit `U2/AUDIT.md` | `U2/` |
| U2' | `= (93/4)ζ(5) − (π^4/8)L(f,2)/L(f,1)` | proved (via Z); 299 digits | `U2/`, `Z/` |
| Z | `L(f,1) = −S_1/2`, `L(f,2) = π^2S_0/16`; Zudilin's (6) for `(1/2,1/2)` with constants `−1/16, 1/32, −3/256` | proved, audited `Z/AUDIT.md`; ingredients all published (assembly only) | `Z/` |
| Z2 | Zudilin's (6) for `(1/2,1/4)`: `F_1 = −8L(f,1)`, `F_2 = 32L(f,2)`, `F_3 = −(448/3)L(f,3)`, `f` = 16.4.a.a | proof written, self-checked (symbolic + 15–30 digit numerics); refereed (`Z/REFEREE-2-Z2Z3.md`: proof complete); relation not found in the literature (L-values known: RWZ) | `Z/THEOREM-2.md` |
| Z3 | Zudilin's (6) for `(1/2,1/3)`: `F_1 = −3L(f,1)`, `F_2 = 18L(f,2)`, `F_3 = −90L(f,3)`, `f` = 36.4.a.a | proof written, self-checked (symbolic + 15–30 digit numerics); refereed (`Z/REFEREE-2-Z2Z3.md`: proof complete); not found in the literature (numerical only before) | `Z/THEOREM-2.md` |
| U3 | `1/PCF[8σ5 + 20σ3 − 5(2n+1), −64n^10] = (93/4)ζ(5) − 4 + (π^6/64)S_0K/(1+S_1K)`, `K = 4 + π^2(Z_0 − 2S_0)` | proved; audit `U3/AUDIT.md` | `U3/` |
| FR-A | Conj. 1.3 of arXiv:2111.04468 for `x_1..x_d ∉ Z_{≥1}`, `B,m ≠ 0`; counterexample outside | proved | `FR/` |
| FR-B | "if" half of Conj. 1.4 | proved | `FR/` |
| FR-C/D/E | deg (1,2), non-terminating: the paper's `GCD2_n` has factorial size iff `a_n` in the families of Conj. 1.4 (`k = 0` allowed); `≤ C^n √(n!)` otherwise and for irreducible `b` | proved, refereed `FR/AUDIT-2.md`; "only if" for `GCD_n` open (= Conj. 3 in this range) | `FR/` §5–§7 |
| FR-6 | printed 1.4 false for `k = 0` in (Eq. 6) (non-terminating, misprint-level); printed 1.1/1.4 false on terminating fractions | proved / verified exactly | `FR/` §6 |
| SR | families F1, F2 and eqs (8), (9) of arXiv:2412.16818 | proved | `SR/` |
| E | three "new and unproven" fractions for e, Table 3 of arXiv:1907.00205v4: `6e/(2e−3)`, `3/(4e−10)`, `4e/(2e−1)` | proved, audited `E/AUDIT.md`; **novelty nil** (instances of David's Thm 1, arXiv:2308.02567; E2 ⇐ Dougherty-Bliss–Zeilberger 2020); no prior explicit statement found only for the first | `E/` |
| RL | arXiv:2412.12361 Table 3: `C[(n+k)/(2n(n+1))] = 1 + ((k+1)/5)·M(k+2;7/2;1/2)/M(k+1;5/2;1/2)`, Möbius in `R_2` for integer `k` (the four "still unproven" rows) | proved, audited `RL/AUDIT.md`; **novelty nil** (Perron's fraction for contiguous `1F1`) | `RL/` |
| RL-m | same table: printed `1/(2R_2−1)` and `(R_2+1)/(6R_2−1)` are false (misprints); correct `1/(2R_2−2)` (proved), `(R_2+1)/(6R_2−6)` (numerical) | refuted as printed | `RL/` §0 |
| RL-C2 | Conjecture 2 of arXiv:2412.12361 = Gauss's logarithm fraction | classical (shown by equivalence transformation) | `RL/` §5 |
| RL-C1 | Conjecture 1 of arXiv:2412.12361 (convergence rules) false as printed: constant `f < −1/4`; `f_n = −1/4 − c/(n(n+1))`, `c > 1/16` | disproved as printed (textbook facts) | `RL/` §6 |

(`σ_d = n^d + (n+1)^d`; convention `PCF[a,b] = a_0 + b_1/(a_1 + b_2/(a_2 + ...))`.)

## Errata found in the sources (printed statements that are false or inconsistent as printed)

"Referee" = confirmed by an independent adversarial audit (file given). "Proved" = proved in this directory;
"numerical" = high-precision evaluation / PSLQ only. All are stated for the conventions of the respective paper.

| # | source | printed statement | what is true | status | referee |
|---|---|---|---|---|---|
| 1 | arXiv:2412.12361v2 (ICLR 2025), Table 3 | `C[n/(2n(n+1))] = 1/(2R_2 − 1)` (`= 0.5490...`), `R_2 = √(πe/2) erf(1/√2)` | `= 1/(2R_2 − 2) = 1.2174747521...` | proved (`RL/`, case `k = 0`) | yes (`RL/AUDIT.md`; also in the ICLR PDF; no v3) |
| 2 | same, Table 3 | `C[(1−2n)/(4(n+2)(n+3))] = (R_2+1)/(6R_2 − 1)` (`= 0.3229...`) | value `0.9783165014...` `= (R_2+1)/(6R_2 − 6)` to 40 digits | printed value refuted numerically; corrected value numerical only | yes (same) |
| 3 | same, Conjecture 1, case "`1 + 4f_n = C + o(1)`, `C ≠ 1` ⇒ converges" | includes `C < 0` | `C[f]` diverges for every constant `f < −1/4` (e.g. `f = −1`: approximants `0, ∞, 1, 0, ...`; `f = −1/2`: period 4; `f = −1/3`: period 6) | proved (`RL/` §6 (A)); textbook | no (written after the audit) |
| 4 | same, Conjecture 1, case "`1 + 4f_n = O(n^k)`, `k < 0` ⇒ converges" | — | `C[−1/4 − c/(n(n+1))]` diverges for every `c > 1/16` (explicit solution `2^{−n}Γ(n+2+σ)/Γ(n+2)`, `σ = 1/2 + i√(4c − 1/4)`; example `c = 1/4`) | proved (`RL/` §6 (B)); sharpness of the parabola theorem | no |
| 5 | same, Conjecture 1, case "`f_n = O(n^k)`, `k < 0` ⇒ `ε_n = O((n!)^k)`" | literal rate | for `f_n = c/n`, `c > 1`: `ε_n + ε_{n+1} ≥ κ c^n n^{−2c}/n!`, not `O(1/n!)` | proved (`RL/` §6 (C)); pedantic (an exponential factor) | no |
| 6 | arXiv:2212.09470v2 (ICML 2023), Table 3 | fraction with signs `{−1,−1,−1,1,1}` (repeated), `a = [1,112,2,143,1,2,2,45,1,1,3,2,53,1,1, 3,240,4,271,...]` `= (−12047J_1(1) + 6928J_0(1))/(8(1777J_0(1) + 3090J_1(1)))` (`= 7.6·10^{−9}`) | value `0.9910312496...` `= (6928J_0(1) − 12047J_1(1))/(8(3090J_1(1) − 1777J_0(1)))` (denominator `8J_6(1)`) | numerical only (60 digits, PSLQ) | no |
| 7 | same, Table 1, row `(−3+5e)/(−6+6e) = 1 + 1/(36 + 1/(2 − 1/(4 − ...)))`, signs `{1,1,−1,−1}` | coefficient data `[2,2,18,3,3,4,3,5,2,2,6,2, 2,6,54,7,...]` | that list gives `2.4862094076...`, not the stated value (no sign pattern of period 12 repairs it); the displayed value with the printed signs has coefficients `[1,36,2,4,1,8,28,2,1,12,2,2, 5,180,6,4,1,24,64,2,1,28,2,2, ...]` | numerical only | no |
| 8 | arXiv:1907.00205v4, Table 3 (also the Nature 2021 supplement) | row `3/(−10+4e) = 3 + 4/(8 + 20/(15 + ...))`: `b_n = −(n+3)n^2` | the displayed fraction and the value need `b_n = +(n+3)n^2` (with the minus sign the value is `2.37297...`) | proved (`E/`, Theorem E2) | yes (`E/AUDIT.md`) |
| 9 | arXiv:2111.04468v2 (Arnold Math. J. 2024), Conjecture 1.4, family (Eq. 6) "for any `k ∈ Q∖{0}`" | "only if": FR ⇒ `a_n` in the printed families | `PCF[n, 2(n+1)(n+2)]` (`k = 0`, `m = 1`, `B = 2`; non-terminating) has FR but is in neither printed family; likewise the paper's own example `PCF[7n + k, 8n^2 − 2]` at `k = 0` | proved (`FR/` Theorem B (i), §6 item 0); misprint-level: `k = 0` must be allowed in (Eq. 6) | found by the referee (`FR/AUDIT-2.md`) |
| 10 | same, Conjecture 1.4 "only if", printed "for every `b_n`" | — | terminating `PCF[n+1, n(n−1)]`: `p_n = q_n ≥ (n+1)!`, FR, but `A^2 + 4B = 5` (not (Eq. 6)) and `a_n ≠ 2kn` (not (Eq. 7)) | proved (`FR/` §6 item 1) | yes (`FR/AUDIT-2.md`) |
| 11 | same, Conjecture 1.1 "only if" | FR ⇒ `b_n` splits over `Q` | terminating `PCF[n^2+1, (n−1)^2(n^2−2)]`: `p_n = q_n ≥ n!^2`, FR, `b` has the irreducible factor `n^2 − 2` | proved (`FR/` §6 item 2) | yes |
| 12 | same, Conjecture 1.4 "if" | families ⇒ FR | `PCF[n−3, 2(n−1)(n−3)]` (family (Eq. 6): `B = 2`, `m = 1`, `k = −3`): `p_n = q_n = 0` for `n ≥ 2`; `PCF[n−1, −(n−1)(n−2)]` (family (Eq. 7), `k = 1/2`): `p_n = q_n = 0` for `n ≥ 1` | proved (`FR/` §6 item 3) | yes |
| 13 | same, Conjecture 1.3 "for any ... `x_1..x_{2d} ∈ Q`" | Euler family has FR | `d = 2`, `B = −3/4`, `m = 1`, `x = (−2, 1, −2, 2)`: `p_n = q_n = 0` for `n ≥ 1` | proved (`FR/` §4 Remark 2) | found by the parent's referee |
| 14 | same, Conjecture 3 "for any PCF: `GCD2_n ≐ GCD_n`" | — | terminating `PCF[n+1, n(n−1)]`: `GCD_n ≥ (n+1)!` while `(ln GCD2_{n+1} − ln n!)/n = −2.09, −2.39, −2.80` at `n = 50, 100, 199` | numerical only | observed by the referee |

Items 3–5 and 10–14 concern ranges the authors presumably did not intend (signs, terminating fractions); items
1, 2, 6, 7, 8, 9 are misprints or inconsistent data in tables/statements. Not errata but relevant: Conjecture 2 of
arXiv:2412.12361 is Gauss's logarithm fraction (`RL/` §5).

**Zudilin's relations (6), per-case status (survey `Z/LITERATURE-13.md`, 2026-11-26; constants `L(f,k)/F_k(1)` recomputed
numerically to ≥ 27 digits by the survey).** No relation with `F_k(1)` was found proved in the literature for any
case except through the ingredients of `(1/2,1/2)`. `k = 3 ⇔ k = 1` in all cases (Zudilin's Remark 5 + functional
equation). "trivial" = `L(f,2) = 0` (root number −1).

| (r,t) | newform | k=1, k=2, k=3 constants | status |
|---|---|---|---|
| (1/2,1/2) | 8.4.a.a | −1/16, 1/32, −3/256 | proved, audited (`Z/THEOREM.md`; assembly of published ingredients) |
| (1/2,1/4) | 16.4.a.a | −1/8, 1/32, −3/448 | proved, refereed (`Z/THEOREM-2.md`, `Z/REFEREE-2-Z2Z3.md`); L-values known (RWZ), relation was numerical |
| (1/2,1/3) | 36.4.a.a | −1/3, 1/18, −1/90 | proved, refereed (`Z/THEOREM-2.md`, `Z/REFEREE-2-Z2Z3.md`); was numerical only |
| (1/2,1/6) | 72.4.a.b | 1/3, 0, −1/468 | numerical only; k=2 trivial |
| (1/3,1/3) | 27.4.a.a | 1/9, 0, −1/243 | numerical only; k=2 trivial |
| (1/3,1/4) | 9.4.a.a (CM) | −1/18, 1/54, −1/216 | L-values known as Γ-products; relations numerical only |
| (1/3,1/6) | 108.4.a.a | 1/2, 0, −1/504 | numerical only; k=2 trivial |
| (1/4,1/4) | 32.4.a.a | 1/8, 0, −3/1280 | numerical only; k=2 trivial |
| (1/4,1/6) | 144.4.a.f | −5/6, 1/36, −5/2304 | numerical only |
| (1/6,1/6) | 216.4.a.c | 5/6, 0, −5/4752 | numerical only; k=2 trivial |
| (1/5,2/5) | 25.4.a.b | −1/5, 4/125, −3/625 | closest to proved: Bönisch–Golyshev–Klemm arXiv:2510.03939 Thm 1 gives the conifold periods as periods of `f`; the dictionary to `F_k(1)` and the rational constants are not written anywhere |
| (1/8,3/8) | 128.4.a.b | 1/2, 0, −3/2816 | numerical only; k=2 trivial |
| (1/10,3/10) | 200.4.a.f | 7/10, 0, −21/34000 | numerical only; k=2 trivial |
| (1/12,5/12) | 864.4.a.a | −13/6, 0, −13/39744 | numerical only; k=2 trivial |

**Stopped here (2026-11-26) on instruction of the parent session** (afterwards one more target was added: Z2, Z3 above, later refereed in `Z/REFEREE-2-Z2Z3.md`; then stopped again on instruction). Not done: Zudilin's relations (6) for the remaining 11 cases (`Z/LITERATURE-13.md`; `(1/2,1/6)` and all cases without `r = 1/2` are outside the method of `Z/THEOREM-2.md`); a uniform proof of the 26 "new and
unproven" rows of arXiv:2212.09470 (24 agree numerically with their printed values; `ESMA/NOTES.md`); the γ
conjecture of arXiv:2610.13919; the CMF conjectures of arXiv:2507.08138; Conjecture 2 of arXiv:2111.04468; the
"only if" halves of Conjectures 1.1/1.4 there for `GCD_n` (= its Conjecture 3 in degree (1,2)).
