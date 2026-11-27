# Audit of Theorem U3 (hostile referee, 2026-11-26)

**Verdict: complete after listed fixes.** I found no gap in the proof of the main formula
`1/PCF = (93/4)ζ(5) − 4 + (π^6/64) S_0 K/(1 + S_1 K)`. Every step checks out by hand and numerically
(§2–§7). The fixes concern statements around the proof: one false "equivalent form" in the boxed theorem
(finding 1) and wrong decimals in Remark 1 (finding 2). The other findings are wording.

I checked U2's results where U3 uses them, and re-read the U2 proofs behind each "word for word" claim (Lemma 3,
Prop. 4, Prop. 5, Cor. 6, Prop. 7, Prop. 8). My own scripts are in `audit/`:
`a1_exact.py` (exact/symbolic, 3 s), `a2_numerics.py` (direct continued fraction and constants, 7 s),
`a3_analytic.py` (quadrature and Frobenius/local-series checks of §3–§6, about 1 min). None of them imports
`verify/` or the pickled constants in `../code/`. I did not run `verify/u3_numerics.py`. It was still running
in another process at the time of the audit.

## Findings

1. **Major (false statement in the theorem box; the main formula is unaffected).** The third form
   `1/U3 = 1/U2 − 4 + π^6 S_0/(64 S_1 (1 + S_1 K))` has the wrong sign. From the main formula and
   `1/U2 = (93/4)ζ(5) + π^6 S_0/(64 S_1)`:
   `(π^6 S_0/64)(K/(1+S_1K) − 1/S_1) = −π^6 S_0/(64 S_1 (1+S_1K))`. Numerically the stated right side is
   `−7.228…`, not `0.04378…` (`a2` §3, `a1` (G)). The sign-corrected form agrees to 60 digits and matches
   U2 §9's L-value form `1/U2 − 4 + π^4 L(f,2)/(8 L(f,1)(1 − 2L(f,1)K))`. None of the `verify/` scripts tests this line.
   **Fix:** write `1/U3 = 1/U2 − 4 − π^6 S_0/(64 S_1 (1 + S_1 K))`, and optionally add the L-value version
   `= 1/U2 − 4 + π^4 L(f,2)/(8 L(f,1)(1 − 2 L(f,1) K))`.

2. **Minor (wrong numbers in Remark 1).** The remark gives `K = −7.78328…`, `1 + S_1K = 6.51836…` and
   `1/PCF = 0.0437808…`. The correct values are `K = −7.7834942513952961922…`,
   `1 + S_1K = 6.5185080678713317792…` and `1/PCF = 0.0437809038563519708…`. The last digit shown for 1/PCF is
   wrong: the value truncates to `0.0437809`. The final 30-digit PCF in §0 is correct.
   **Fix:** replace the three numbers.

3. **Minor (overclaim in §0).** "Numerically (PSLQ, 60 digits, three bases) `Z_0` is not a combination of
   `L(f,1), L(f,2), π`" says more than PSLQ can show. PSLQ only rules out relations up to a coefficient-height
   bound. **Fix:** write "PSLQ at 60 digits finds no integer relation among {…} with coefficients of size below
   H (state H)". Remark 2 already words this correctly.

4. **Minor (status line).** §0 lists `verify/u3_numerics.py` as "passing". When I audited, it had been running
   for more than 20 minutes and had not finished. Its part (6) also compares against a pickled 48-digit value
   from `../code/u_consts.pkl`, not against a value computed inside the script.
   **Fix:** state this only after the run has finished, and record its output. Optionally cite
   `audit/a2_numerics.py`, which evaluates the fraction directly from the recurrence (N = 1.28·10^5 plus
   Richardson extrapolation) and agrees with the closed form to 34 digits.

5. **Cosmetic (uniqueness remark after Lemma 3.3).** "The choice `c_2 = −1/16` is the one for which this
   cancellation occurs" holds modulo `ν` and scaling. Adding `λν` to `r_*` leaves `μ_*` unchanged, because
   `νθν − νθν = 0`. Adding `λr` adds `λ G(g)_n = 128λn + 64λ + O(1/n)` (checked in `a1` (E)).
   **Fix:** add "(up to adding multiples of ν and scaling)".

6. **Cosmetic (Prop. 3.2 proof).** I checked the "word for word" claim, and it does transfer. In detail:
   - (i): the boundary terms at `x = 0` need only `u > 0` and the growth `O(log^4)`. The terms at `x = 1` are
     exactly `B_5(u)`, which vanishes for `μ_*` with `(θ^iμ_*(1)) = (1, −1/2, 1/4, −1/8, 7/128)`. I recomputed
     this with my own integration by parts, and it also vanishes for a general holomorphic `ρ`.
   - (ii): the expansion `μ_* = μ_4^* log^4 x + O(1 + |log x|^3)` near 0 does not care that `μ_*(1) ≠ 0`.
     Away from 0, `μ_*` is bounded, so its contribution to `u^5 g_*(u−1)` tends to 0.
   - (iii): it uses only the Frobenius coordinates of `r_*` and U2's identity `∫(νθy_j − y_jθν) = ξ_j`. That
     identity and its Fubini justification involve only `ν` and `y_j`, so they are unaffected.
   - U2's bound `g(n) = O(n^{−3})` and Cor. 6 fail for `g_*`, since `g_*(n) ~ 1/n`. The text correctly says
     that they are not used.

   **Fix (optional):** add one sentence listing these points explicitly.

7. **Cosmetic (cross-file consistency).** `../U2/THEOREM.md` §9 still says "U3 … not proved", and it calls the
   reduction of the Taylor coefficients of `Z` to `Z_0` "found by PSLQ and not proved". That reduction is now
   Prop. 6.1 plus the Legendre relation. **Fix:** add a pointer in U2 §9 to U3/THEOREM.md.

### Items checked and found correct (no fix needed)

- **Lemma 2.1:**
  - `a_n/8 = ã_3(n)`, and `ã_3 − (n+1)^5 − n^5 = 5(2n+1)(4n^2+4n+3)/8`.
  - The monotonicity induction holds.
  - The Casoratian is `1/(8(n+1)^5)`, checked exactly for n < 200.
  - The identification `p_{n−1} = 8^n n!^5 w3_n`, `q_{n−1} = 8^n n!^5 v3_n` was checked against the convergents
    computed directly from `a_n, b_n`, with `(p_{−1},p_0) = (1,a_0)` and `(q_{−1},q_0) = (0,1)`, exactly for n ≤ 40.
- **Lemma 2.2:**
  - I eliminated `X_{n−1}` and `X_{n+2}` myself: `R_3` applied to `G(X)` is ≡ 0.
  - `G(w)_0, G(w)_1 = 1, 23/8` and `G(v)_0, G(v)_1 = 4, 93/8`.
  - `G(v) = v3 + 4w3` holds exactly for n < 202.
- **Prop. 3.1:**
  - The local series of `r_*` and `μ_*` were recomputed with an independently derived coefficient recurrence
    (`L4 s^m = Σ_{j=−3}^{1} P_j(m) s^{m+j}`, pivot `m(m−1)^2(m−2)`).
  - `θ^i r_*(1) = (1, 0, −1/8)`, and the Lemma-3 right side is `−1/8 − u(u+1/2)`.
  - The relation for `R_*` holds, by hand and with sympy.
  - The growth estimates hold, using `0 < Z ≤ Y ≤ u+2` for `u > 0`. `1/G ~ u^2` gives a periodic function that
    is `O(1/u)`, hence 0.
  - The Mellin transform agrees numerically at u = 0.7 and 2.5 to about 1e−24.
  - The connection coefficients, solved from four conditions at x = 1/2, match `2(−1)^j η*_{3−j}` to 1e−31 and
    re-check at x = 0.3.
- **Lemma 3.3:**
  - The decomposition with `c(n) = 32n^4+112n^3+152n^2+92n+21` is correct.
  - `ρ_4` contributes `O(1/n)`.
  - `Σ_{j≤3} m_j G(β_j)_n` has numerator degree 4 and denominator degree 5, and behaves like `839/n`.
  - `G(β_0) ~ 16n^2` and `G(β_1) ~ −32n^2`.
  - Numerically, `G(g_*)_n` decays (about 5e−16 at n = 399), so it is in fact much smaller than `O(1/n)`.
- **Prop. 4.1 / 5.1:**
  - The signs, the factor −192 (from `ad − bc`), and `−24X = −4 ⇒ X = 1/6` are correct.
  - Numerically, `−(n+1)^5 C[g,g_*]_n/24 = 0.1666…667` (22 digits) for n = 0..3, using g and g_* from quadrature.
- **Prop. 6.1:**
  - The reflection `γ(z) = tan^4(πz) γ(−z−1/2)` is correct.
  - The decomposition (6.2), with `m = j+1` and `u − m + 1 = −(v + j − 1/2)`, is correct.
  - The correction term is `tan^4(πu)/(u N(1/2−u)) = 16π^2 u^3 + O(u^4)`, using `N(1/2) = π^2/16`.
  - At `±i∞`, `|tan πu| ≤ coth π` for `|Im u| ≥ 1`, giving `O(|Im u|^{−2})`.
  - Liouville: `sec^2 tan` has no `v^{−1}` term, since it is a derivative of `csc^2(πv)`, and the limits of
    `cπ tan` at `±i∞` are `±iπc`.
  - The `e_2` limit is legitimate: `Z/N` is analytic at 1/2, and `Z~(v)` converges uniformly near `v = 0`
    (terms `O(j^{−3})`). Numerically, `T_Z(1/2−δ)/t^4 → e_2` with an error of about `e_1πδ`.
  - The telescoping step holds: `H(0) = 0`, `H(k) → 1/π^2`.
- **§7:**
  - I redid the coefficient extraction by hand and with my own sympy bookkeeping.
  - `[ε^3](N Y_*(−ε)) = 0`, so the `ζ(3)` term cancels.
  - `192μ_4^* = E_1` and `g_*(0) = (93/4)ζ(5)E_1 + (π^2/16)(e_0+e_2−4τ_0)`.
  - The Legendre relation gives `E_1 = 8(1+S_1K)/(π^4S_0)`, and the by-product `e_1` is correct.
  - Signs: `2S_0 − Z_0 ≥ A_0 = 1 > 4/π^2`, so `K < 0`, `1 + S_1K > 1` and `E_1 > 0`.
- **Remark 2:** the elliptic-integral formula for `Z_0` holds numerically to 1e−31.
- **Source:** the untrusted text extraction of arXiv:2308.11829 shows Table 1, sixth row:
  `8(n^5+(n+1)^5) + 20(n^3+(n+1)^3) − 5(2n+1)`, `−64n^10`, `22.8410...`, and the quoted sentence
  "Three of the discovered formulas … did not fit any combination of zetas (in numerical tests using PSLQ)".
  Table 4 has the same row (first of "Lines 5 to 7 … did not find a closed-form expression"). The convention
  matches: my direct evaluation gives 22.8410….

## Numerical table (own code, independent of `verify/`)

| quantity | value / residual | source |
|---|---|---|
| PCF, direct (N = 128000, Richardson) | 22.84100856576798561590966746243729 (last Richardson increment 4e−27) | a2 §1 |
| 1/PCF direct − closed form | 7.7e−35 | a2 §3 |
| closed form 1/(…) | 22.84100856576798561590966746243732787 | a2 §3 |
| S_0 (own sum + Richardson; vs 4F3) | 1.1186363871641870683496192575256 (Δ 1.5e−39) | a2 §2 |
| S_1 (own sum; vs parameter-derivative of 5F4) | −0.70900136746192943753111978298986 (Δ 3e−38) | a2 §2 |
| Z_0 (own sum; vs 5F4) | 1.0433551893343198886076074043165 (Δ 1.8e−38) | a2 §2 |
| K, 1 + S_1K | −7.78349425139529619, 6.51850806787133178 | a2 §3 |
| third form as stated − formula | −7.27 (**wrong sign**); sign fixed: 6e−61 | a2 §3 |
| Σ A_k/(k−1/2) + 8/π^2 + Z_0 | 1.8e−38 | a2 §3 |
| τ_2 − 2π^2τ_0, τ_3 − (4π^2/3)τ_1 | 9e−30, −3e−30 | a3 |
| ζ_2 − π^2(e_0+e_2), ζ_3 − ((4π^3/3)e_1 − 16π^2) | 1.2e−29, 1e−29 | a3 |
| e_1 vs 8(1+4S_1+π^2S_1Z_0)/(πS_0) | −20.799475356133531439898850 (both) | a3 |
| connection coefficients ν, r, r_* (solved at x = 1/2 vs formula) | 6e−37, 5e−30, 1e−31 | a3 |
| Prop 3.1 Mellin at u = 0.7, 2.5 | 7.6e−24, −1.6e−24 | a3 |
| Prop 3.2(i) residual at u = 0.37, 1.6 | 5e−11 (quadrature of x^{−0.63}log^4), 6e−30 | a3 |
| μ_4^*: connection / formula / from g_*(1)−ã(0)g_*(0) | 0.0024925767106439338804 (all three) | a3 |
| g_*(0): quadrature / −2[ε^3](Ξ Y_*(−ε)) | 1.935251347948611770291 (both) | a3 |
| (3.1) at n = 2..5 | agreement to all 22 printed digits | a3 |
| Prop 5.1: −(n+1)^5 C[g,g_*]_n/24, n = 0..3 | 0.1666666666666666666667 | a3 |
| Prop 4.1: g_*(0)/(192μ_4^*) − 4 | 0.04378090385635197083893732 (= 1/PCF direct) | a3 |
| T_Z closed form at u = 0.2, 0.3+0.4i, −0.35 | 1.5e−26, 1.7e−25, 1.4e−25; direct bilateral sum at 0.2: 1.9e−24 | a3 |
| Remark 2 elliptic integral − Z_0 | 1.4e−31 | a2 |

## Prior art

I found no prior identification of 22.8410085657679856… (or of 1/PCF = 0.04378090385…). I also found none of
`Z_0 = 5F4(1/2,1/2,1/2,1/2,1;1,1,1,2;1) = Σ C(2k,k)^4/(256^k(k+1))` as a quasi-period of `η(2τ)^4η(4τ)^4`.

- arXiv:2308.11829 (source) lists the value only numerically: https://arxiv.org/pdf/2308.11829
- Yamamoto, "Proof and generalization of conjectures of Ramanujan Machine", arXiv:2403.09729
  (https://arxiv.org/abs/2403.09729). Its abstract does not mention this row. I did not check the full text
  line by line.
- Exact-value searches for "22.84100856", "1.0433551893" and "0.0437809038" returned nothing relevant.
- Related context (quasi-periods at conifold fibres of the 14 hypergeometric CY families, the framework cited in
  Remark 2):
  - Bönisch–Klemm–Scheidegger–Zagier, arXiv:2203.09426 (https://ar5iv.labs.arxiv.org/html/2203.09426)
  - The CMP version: https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1007/s00220-024-05006-6/ADEK-CMP24.pdf
  - Neither gives this 5F4 constant explicitly in what I read; a closer reading of their §4 tables is advisable
    before any novelty claim about `Z_0` itself.
- The level-8 L-value `L(f,2) = (π^2/16) 4F3(1/2^4;1^3;1)` (Zagier), as quoted in "The Explicit
  Hypergeometric-Modularity Method I", arXiv:2404.00711 (https://arxiv.org/pdf/2404.00711). This is only
  background for `S_0`.

"New" therefore still means "not found after this search", as the theorem file says.
