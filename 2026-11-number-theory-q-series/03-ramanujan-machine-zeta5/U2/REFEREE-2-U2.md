# Referee report 2 on U2/THEOREM.md (independent, adversarial)

Date: 2026-11-26. I formed this verdict before opening `AUDIT.md`. Afterwards I only skimmed it for its list of
issues, which agrees with mine. My scripts are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/u2ref/` (`r1_exact.py` ... `r5_u3.py`).

## Verdict

**PROOF COMPLETE** for Theorem U2, i.e. the hypergeometric form
`1/PCF = (93/4)ζ(5) + (π^6/64) S_0/S_1`, convergence of the fraction, `p_{n−1} = 8^n n!^5 w_n`, and the series
for `8/PCF`. I checked every step of §1–§7 by hand and found no fatal error and no real gap. The only places
that need a short extra sentence are routine (findings 5–7). Three sentences outside the proof state more than
is proved (findings 1–3). The L-value Corollary is honestly labelled conditional, and §1–§7 do not use it.
§9 (U3) is honestly labelled numerical.

## Reading of the statement

`S_0 = Σ C(2k,k)^4/256^k` converges absolutely, with terms `~ 1/(π^2 k^2)`. `S_1 = Σ A_k (H_{2k} − H_k − log 2)`
also converges absolutely: `H_{2k} − H_k − log 2 = −1/(4k) + O(k^{−2})`, so the terms are `~ −1/(4π^2 k^3)`. No
regularisation is involved anywhere. The doc's line 18 ("every term ≤ 0; terms O(k^{−3})") is correct. In the
proof, `8 S_1 = D = Y'(0) − 8 log 2 · Y(0)` is obtained by differentiating the series for `Y` term by term. That is
legitimate because `Y` is a locally uniformly convergent series of analytic functions.

## Findings

### Fatal
None.

### Gaps
None that affect the theorem. Findings 5–7 are routine estimates that should be written out.

### Minor

1. **Line 38, "So this constant is *not* a combination of zeta values".** This is not proved. No independence
   or irrationality result for `S_0/S_1` relative to zeta values is available or claimed.
   *Fix:* "is not visibly / not expected to be a combination of zeta values (PSLQ finds none)".
2. **§8 Remark 2, "by the proof this constant is −D/π^2".** The proof contains no asymptotic analysis of `w_n`.
   Only `w_n ≥ 1/(16 n^2)` and `g(n) = O(n^{−3})` are used.
   *Fix:* either label it numerical, or add a singularity analysis of `ψ_0` at `x = 1` (transfer theorem).
3. **§8 Remark 3, last sentence** (invariance of `r` under `y(x) ↦ x^{−1/2} y(1/x)`). This is not proved; it is a
   heuristic reading of Prop. 8. *Fix:* label it as such.
4. **Remark 2a / Prop. 8: attribution and a missed consequence.**
   - I checked Zudilin's source (arXiv:1805.00544v2, `sigma18-086.tex`). His normalised Frobenius generating
     function at `r = t = 1/2`, `z = 1` is `π^{−2}·Y(ε)/N(ε)`, so `F_j(1) = τ_j/π^2`. Prop. 8's first relation is
     therefore exactly his `F_0(1) = F_2(1)/(2π^2)`, and the proof of it is correct.
   - The second relation, `τ_3 = (4π^2/3) τ_1`, is the case `(1/2,1/2)` of the statement after Zudilin's Lemma 3
     ("$F_1(1)$ and $F_3(1)/\pi^2$ are rationally proportional", which he calls rigorously established).
   - The method (a bilateral series that is periodic in ε, Golyshev–Mellit) is the one in Zudilin §3.
   - *Fix:* acknowledge both points.
   - Also worth stating: with eq. (4) (`S_0 = 16L(f,2)/π^2`, proved by Zagier and by Rogers–Wan–Zucker), Prop. 8
     gives `F_2(1) = 32 L(f,2)`. This proves the second inclusion `L(f,2)/F_2(1) ∈ Q` of Zudilin's (6) for
     `(r,t) = (1/2,1/2)`. Only the first inclusion, `L(f,1)/F_1(1) = −1/2`, remains conditional; that is the
     Corollary's hypothesis.
5. **§3, convergence of the formal solutions at `x = 1`.** The fact used is true: every formal power-series
   solution at a regular singular point of a Fuchsian scalar equation converges in the largest disc free of
   other singularities. The citation is imprecise, though. Coddington–Levinson ch. 4 treats first-order
   systems, so the scalar-to-system reduction (or a majorant argument) should be named.
   *Fix:* cite precisely, or add the standard three-line majorant argument based on the pivot
   `(p+1)(p+2)^2(p+3)`, which I recomputed: 12, 72, 240, ... .
6. **§5 first paragraph, and Prop. 7's interchange justification.**
   - The bound is stated only on `[1/2,1]`. On `(0,1/2]` one also needs `|ν|, |θν| ≤ C|log x|^3`, which makes
     `x^{k−1/4}|log x|^3` integrable (also for `k = 0`).
   - The same uniform bound `|A_k(ε)| ≤ C(k+1)^{−2}` should be stated explicitly for `ŷ_j(u) = [ε^j]Σ A_k(ε)/(k+ε+u)`.
   - Trivial to add.
7. **Prop. 5(iii), `g(n) = O(n^{−3})`.** The doc gives the `[1/2,1]` part. The `[0,1/2]` part is
   `≤ C 2^{−n}` because `|μ| ≤ C log^4 x` there. Add one clause.
8. **Prop. 2, `A_n ≥ 1/(16 n^2)`.** This holds for `n ≥ 1` only (`n = 0` is trivial). Cosmetic.
9. **§9, "Exact (proved; v6)".**
   - The claim is correct. `v6` verifies symbolically that the gauge map sends solutions of R_2 to solutions
     of R_3 for `n ≥ 1`. A solution of R_3 is determined by its values at `n = 0, 1`, and `v6`'s finite check
     includes those.
   - The text should say that this is why a finite check proves it for all `n`.
   - The rest of §9 is correctly labelled "not proved, numerical".

## Step-by-step status (§1–§7)

| Step | Status |
|---|---|
| §1: Frobenius basis, `L4Φ = ε^4 x^ε`, `λ_k = 8(H_{2k}−H_k)`, `N'(0)/N(0) = 8 log 2`, `D = 8S_1 < 0` | correct |
| Lemma 1 (`L5` is the exterior square of `L4`) | correct. I recomputed it independently in sympy: `L5(fθh − hθf) = 0` modulo `L4` |
| Prop. 2: `w_n` formula, positivity, recurrence from `L5ψ_0 = 0` | correct. The coefficient identity is valid because `ψ_0` is a convergent power series |
| Prop. 2: identification of `p_n, q_n` for **all** `n` | correct and proved, not just checked. Substituting `U_{n−1} = 8^n n!^5 c_n` turns the fraction's recurrence into `(n+1)^5 c_{n+1} = ã(n)c_n − n^5 c_{n−1}` for `n ≥ 1`. The initial values `p_{−1} = w_0 = 1`, `p_0 = 8w_1 = a_0`, `q_{−1} = v_0 = 0`, `q_0 = 8v_1 = 1` match, so induction finishes it. I verified `ã = a/8` exactly |
| Casoratian `D_n = 1/(8(n+1)^5)` and monotonicity of `v_n/w_n` | correct |
| Convergence | correct. It is proved in Cor. 6 via `g(n)/w_n → 0`, so `ℓ < ∞`; `ℓ > 0` holds by monotonicity |
| §3: holomorphic solutions at `x = 1` | correct (exponents 0,1,1,2; three-dimensional holomorphic space). Citation: finding 5 |
| `ν`, `r` coefficients and `θ^iν(1)`, `θ^ir(1)` | correct (recomputed exactly: `ν = s + 3/4 s^2 + 29/48 s^3`, `r = s^2 + 5/4 s^3 + 377/288 s^4`) |
| Lemma 3 boundary polynomial `F_2 + F_1/2 − u(u+1/2)F_0` | correct (recomputed for a generic holomorphic solution) |
| Prop. 4(i) | correct. The periodic function lives on the real half-line and has a limit at `+∞`, so it is constant. No Carlson-type growth argument is needed |
| Prop. 4(ii) | correct. Checked the telescoping of `R_0` and the bound `Y(u) ≤ u+2` |
| Prop. 5(i): `B_5 ≡ 0` from `(θ^iμ(1)) = (0,0,−2,3,−9/4)` | correct (recomputed independently) |
| Prop. 5(ii): `u^5 g(u−1) → 24μ_4` | correct |
| Prop. 5(iii) and Cor. 6 | correct (two initial values plus a second-order recurrence with nonzero leading coefficient). Finding 7 |
| (5.1) by comparing principal parts at `u = 0` | correct. Both sides are meromorphic near 0 and agree on `(0,δ)` |
| Prop. 7: `μ_4 = N_0 D/6`, the `1/12` coefficient, `g(0) = −2[ε^3](ΞY(−ε))`, (5.2) | correct. Rederived by hand, including `A_k(ε)Γ(k+ε+1)^4/Γ(k+ε+3/2)^4 = N(ε)/(k+ε+1/2)^4`. Finding 6 |
| Prop. 8: bilateral periodicity argument | correct. `T` is meromorphic with period 1 and poles only on `1/2+Z`. The decay `O(1/|Im u|)` holds uniformly on the strip (Stirling for `k ≥ 0`; for `k < 0` the reflected form with `|tan πu| ≤ coth π`). Principal parts at `1/2` are correct (`cot^4 z = z^{−4}(1 − 4z^2/3 + …)`, `S_2, S_3, S_4` expansions checked). `E` is entire, periodic and bounded, hence constant by Liouville; its limits `∓iπc` force `c = 0`. Evaluating at `u = 0` gives `τ_2 = 2π^2τ_0`. Nothing is missing; no Carlson theorem is needed |
| §7 algebra | correct. Rederived by hand: `[ε^3]` coefficient `= −(2/π^4)(π^6τ_0 + 186ζ(5)τ_1)`, and the final formula |

**Use of the conditional identity.** §1–§7 use only the definitions of `S_0` and `S_1`. `L(f,·)` appears only in
the Corollary, §8 and §9. The Corollary is stated as conditional on `8S_1 = −16L(f,1)`, i.e. Zudilin's
`L(f,1)/F_1(1) = −1/2`. I confirmed in the source that Zudilin's remark following eq. (6) says none of the
relations (6) "seem to be proved".

## What I recomputed

- **Exact (sympy):**
  - `ã(n) = a_n/8`.
  - Lemma 1, by my own implementation in the differential ring.
  - The local pivot at `x = 1`, `ν`, `r`, `μ` and their θ-values.
  - The Lemma 3 boundary polynomial for a generic holomorphic solution.
  - `B_5(u) ≡ 0`.
- **The fraction itself:**
  - `p_n, q_n` computed in exact integers for `n ≤ 3000`.
  - Neville/Richardson extrapolation of `q_n/p_n` in `1/n` gives `PCF = 2.45174606462878291704788012916881`
    (12 nodes).
  - The formula, with `S_0` (via `hyper`) and `S_1` (via `nsum`, and independently via `D/8` from `mp.diff`),
    gives `2.45174606462878291704788012916881306`. They agree to about 31 digits.
- **Prop. 8 numerically (40 digits):**
  - `τ_2 − 2π^2τ_0 ≈ 3e−39` and `τ_3 − (4π^2/3)τ_1 = 0`.
  - The displayed identity holds at `u = 0.2+0.3i` to `1e−40`.
- **Corollary (conditional):** I computed `L(f,1) = 0.3545006837309647…` and `L(f,2) = 0.6900311631233975…`
  from the functional equation (level 8, sign +1). Both `S_0 − 16L(f,2)/π^2` and `8S_1 + 16L(f,1)` are below
  `3e−30`.
- **§9:** Richardson extrapolation of the U3 fraction agrees with both printed formulas for `1/U3` to about
  30 digits (`0.0437809038563519708389373215186033`).
- **Not checked:** the 95/299-digit claims, the BKSZ (arXiv:2203.09426) statement in §0, and the `../code/*`
  scripts.
