# Adversarial audit of Theorem U2 (THEOREM.md, 2026-11-26)

## Verdict

**COMPLETE AFTER LISTED FIXES.** I found no fatal error and no wrong step in §1 to §7. The proof of the
hypergeometric form of Theorem U2 is correct. Issues 1 to 4 are small justification gaps that one or two
sentences each would close. Issues 5 to 8 are wording or citation problems. The Corollary is correctly
labelled as conditional. My literature search found no proof of the identity `8 S_1 = −16 L(f,1)` with its
exact constant. However, the §0 description of BKSZ is inaccurate (issue 6).

## Issues

1. **§3, the sentence "convergent for |s| < 1"** (fillable gap). The pivot argument gives only a *formal*
   solution. Convergence needs the standard theorem that a formal power-series solution of a Fuchsian equation
   at a regular singular point converges up to the nearest other singularity. Fix: cite that theorem, for
   example Ince §16 or Coddington–Levinson Ch. 4, or give the majorant argument. The singularities in `s` are
   at 0, 1 and ∞, so the radius is 1, as the text says.
2. **§5, first paragraph; Prop. 7, "uniform convergence"** (fillable gap). Two interchanges need an explicit
   bound: `[ε^j]` with `∫ Σ_k`, and termwise integration of `∫ ν θΦ dx`.
   - The series `θΦ = Σ (k+ε) A_k(ε) x^{k+ε}` does not converge uniformly up to `x = 1`, because its terms
     are about `k^{−1}`.
   - Fix: state that `|A_k(ε)| ≤ C (k+1)^{−2}` uniformly for `|ε| ≤ 1/4`. This holds because
     `(1/2+ε)_k/(1+ε)_k ~ (Γ(1+ε)/Γ(1/2+ε)) k^{−1/2}`.
   - Then use `|ν(x)| ≤ C(1−x)` on `[1/2,1]` and `O(|log x|^3)` on `(0,1/2]`, which give
     `∫_0^1 |ν| x^{k} dx = O(k^{−2})`. Now `Σ_k k · k^{−2} · k^{−2} < ∞`, so Fubini applies.
   - Cauchy's estimates on `|ε| = 1/4` then justify taking `[ε^j]` termwise. The same bound makes
     `Σ_{k≥1} [ε^j] A_k(ε)/(k+ε+u)` analytic for `|u| < 1/2`.
3. **Prop. 5(ii), "μ = μ_4 log^4 x + O(|log x|^3)"** (fillable gap, one line). The full expansion of `μ` from
   (1.2) also contains terms such as `log^4 x (φ_0^2 − 1)` and `log^5 x · x(...)`, which come from `θφ_m`. They
   are `O(x |log x|^5) = O(1)` on `(0,1/2]`, so the claim holds, but the text should say so. With this the
   error term in `u^5 g(u−1) → 24 μ_4` is justified: `u^5 ∫_0^1 O(1+|log x|^3) x^{u−1} dx = O(u)`.
4. **Prop. 8, decay at ±i∞** (fillable gap, cosmetic). Say "for `|Im u| ≥ 1`". There `|z| ≥ 1` for every
   term, so the Stirling bound `|Γ(z+1/2)/Γ(z+1)| ≤ C|z|^{−1/2}` holds uniformly on `Re z ≥ 0`. Also
   `|tan πu| ≤ coth π < 2` there. Boundedness of `E` on the rest of the strip follows from continuity. The
   argument is otherwise correct, and the limits `tan(πu) → ±i` are right.
5. **Prop. 8, display** (cosmetic). The text says "for all `u ∉ 1/2 + Z`" and then "(−1/2 < Re u < 1)".
   - Inside that strip the left side, written with `Y/N`, equals `T`, and `T` is defined for all
     `u ∉ 1/2+Z`.
   - Fix: state it on the strip only, or state it for `T`.
   - In Prop. 7, the letter `k` serves as both the summation index and the `[ε^k]` index. Rename one of them.
6. **§0, status paragraph** (cosmetic, literature). "BKSZ … treat such conifold periods numerically" is
   inaccurate.
   - BKSZ (CMP 405 (2024) 134), §4, *prove* for exactly this family (N = 8, the mirror of four quadrics in
     P^7) that `Φ*Ω = 2(2πi)^3 f(τ) dτ∧dz_1∧dz_2` with `f = η(2τ)^4η(4τ)^4`. Hence the periods of
     Π-anti-invariant cycles at the conifold have the form `α^+ ω_f^+ + α^− ω_f^−`. They prove the same for
     quasiperiods, but only "up to a multiplicative constant" (their eq. (4.25)–(4.28)).
   - The exact rational entries (their Table 3, eq. (3.11)) are numerical.
   - Fix: rephrase. Optionally, check whether their §4 can be made to give the constant in `8S_1 = −16L(f,1)`.
     The Corollary should stay conditional until then.
7. **§9, first bullet labelled "proved; v6"** (cosmetic). `v6_u3_gauge.py` proves two things:
   - the gauge identity, symbolically in `n`;
   - `G(w) = w^{(3)}` and `G(v−4w) = v^{(3)}`. The check covers only `n < 40`, but both sides satisfy `R_3`
     and agree at `n = 0, 1`, so it is a proof.

   It does **not** prove "exchanges dominant and minimal behaviour" or the "Consequently `1/U3 + 4 = …`"
   sentence. That sentence needs the existence of asymptotic expansions of `w_n` and `v_n` in powers of
   `n+1/2`. Move those two statements out of the "proved" bullet, or prove them. §9 is otherwise clearly
   labelled as not proved.
8. **§1** (cosmetic). "`Y` is analytic for `Re ε > −1/2`" is true but undersells the result: the terms are
   `O(k^{−2})` for all `ε ∉ −N`. Nothing depends on this.

## Line-by-line results (all steps not listed above: correct)

- **§1.** `L4 Φ = ε^4 x^ε` is correct: the ratio `A_k/A_{k−1} = ((k−1/2+ε)/(k+ε))^4` cancels every term with
  `k ≥ 1`. So `y_0..y_3` are solutions, and their leading `log^j` terms make them a basis. (1.2) is correct.
  `λ_k = 8(H_{2k}−H_k)` is correct. `N'(0)/N(0) = 4(ψ(1)−ψ(1/2)) = 8 log 2` is correct. `D = 8S_1 < 0` is
  correct.
- **Lemma 1.** In `v1`, `top` encodes `(1−x)θ^4y = x(2y_3 + (3/2)y_2 + (1/2)y_1 + y_0/16)`, which is exactly
  `θ^4y = x(θ+1/2)^4y`. The functions `theta` (chain rule on the jet variables) and `shifted`
  (`(θ+c)^k m`) are correct, and the `L5` used matches §1. The computation is an exact rational-function
  identity in the jet variables, so it is a legitimate proof on `(0,1)`.
- **Prop. 2.** I re-derived `ψ_0 = y_0^2 + y_0θg_1 − g_1θy_0`; the logs cancel, giving the coefficient
  `Σ A_jA_m(1 + (m−j)λ_m)`. The symmetrisation is correct. The bounds `(1/2)_n/n! ≥ 1/(2√n)`, giving
  `A_n ≥ 1/(16n^2)`, are correct. The gauge `U_{n−1} = 8^n n!^5 c_n` is correct. The base cases are correct
  (`a_0 = 3 = 8·(3/8)`). The Casoratian satisfies `(n+1)^5 D_n = n^5 D_{n−1}` for `n ≥ 1`, so
  `D_n = 1/(8(n+1)^5)`.
- **§3 / Lemma 3.** I derived the boundary polynomial by hand.
  - The `θ^4` term contributes `F_3 − uF_2 + u^2F_1 − u^3F_0`. The `−x(θ+1/2)^4` term contributes
    `−F_3 + (u−1)F_2 + (1/2−u^2)F_1 + (u^3−u^2−u/2)F_0`.
  - Their sum is `−F_2 − F_1/2 + u(u+1/2)F_0`, which matches the lemma and `v2`(c).
  - The integration-by-parts rule and its hypotheses (`c > −1`, polylogarithmic growth at 0) are correct.
  - `ν`: with `θ = −(1−s)d/ds`, I get `θν(1) = −1` and `θ^2ν(1) = 2c_2 − c_1 = 1/2`. The condition
    `F_2 + F_1/2 = 0` forces `c_2 = (3/4)c_1`, as in the text.
  - `r`: `(0, 0, 2)`.
- **Prop. 4.** `π_1(u+1) = π_1(u)` on `(0,∞)` is correct. A periodic function with a limit at `+∞` is constant.
  `ν̂(u) = B(u,2) + O(u^{−3}) + O(2^{−u}/u)` is correct. The shift identity for `R_0` is correct. In the bound
  `Y(u) ≤ (u+1)^2[(u+1)^{−2} + ∫_0^∞(u+1+t)^{−2}dt] = u+2`, every inequality holds. So `π_2 = O(u^{−1}) → 0`.
- **Prop. 5.** The multipliers `−u^5`, `+ã(u)`, `−(u+1)^5` with arguments `g(u−1)`, `g(u)`, `g(u+1)` are
  correct; the sign `ã(u) = −2P(−(u+1))` is checked. Recomputing from `μ = −s^2 − (3/2)s^3 − (169/96)s^4`, I
  get `θ^iμ(1) = (0, 0, −2, 3, −9/4)`. The routine `mellin_boundary` in `v2` reproduces my hand-derived L4
  polynomial, so its `B_5` is trustworthy, and `B_5 ≡ 0`. `g` is continuous on `(−1,∞)` by dominated
  convergence. (iii) is correct: the second-order recurrence has a nonzero leading coefficient, and
  `X_1 = ã(0)g(0) − 24μ_4 = g(1)`.
- **Cor. 6.** The logic is correct: `ℓ = g(0)/(192μ_4)` is finite, `ℓ > 0` because `v_n/w_n` increases
  strictly from 0, so `g(0) ≠ 0`. `μ_4 ≠ 0` is supplied by Prop. 7.
- **§5 / Prop. 7.**
  - The principal parts are correct: `[ε^j](u+ε)^{−1} = (−1)^j u^{−j−1}`, and (5.1) follows, with signs
    correct.
  - Only `(2,3)` and `(3,2)` contribute to `log^4`, since `(3,3)` cancels by antisymmetry. The coefficient is
    `1/(2!2!) − 1/(3!1!) = 1/12`, giving `μ_4 = (N_0η_1 − N_1η_0)/6 = N(0)D/6`.
  - `∫(νθx^a − x^aθν) = (2a+1)ν̂(a+1)` with `a = k + ε`. The Gamma identity is correct, the sum equals
    `2ζ(3,1/2+ε)`, and `g(0) = −2[ε^3](ΞY(−ε))` and (5.2) are correct, signs included.
- **Prop. 8.**
  - The reflection formula is correct, and both halves identify with `Y/N` at `u` and at `1/2−u`.
  - In the strip, poles occur only at `u = 1/2`. The zeros of `1/Γ(u+k+1)` are zeros of the terms, not poles.
    Periodicity puts the remaining poles at `1/2+Z`, and `T` is analytic at the integers.
  - At `u = 0`, `tan^4 = O(u^4)`. I checked the principal part at `1/2` using `cot^4 z = z^{−4}(1 − 4z^2/3 + …)`.
  - I checked `S_2, S_3, S_4` (csc/cot expansions) and the Liouville argument.
  - As a consistency check, `T_0 = τ_0(1+2π^2u^2) + τ_1(u + (4π^2/3)u^3) + O(u^4)`.
- **§7.** `N(ε)N(−ε) = (ε cot πε)^4` is correct. The Hurwitz coefficients `7ζ(3)`, `−3·15ζ(4) = −π^4/2`,
  `6·31ζ(5)`, `−10·63ζ(6) = −2π^6/3` are correct. I redid the final algebra by hand and the τ_0-terms sum
  to `−π^6/3 − 2π^6/3`. The relations `τ_0 = π^2Y(0)` and `τ_1 = π^2D` are correct. The Theorem and the
  Corollary formulas follow (I checked the algebra between them). The series identity is `8ℓ`, by telescoping.

## Literature (item 11)

- **Zudilin, SIGMA 14 (2018) 086 (arXiv:1805.00544v2), read in full around eqs. (4)–(6).**
  - Eq. (4), `4F3(1/2^4;1^3;1) = 16L(f,2)/π^2`, is stated as **proved**. He writes "they prove [34]" (Zagier,
    *Arithmetic and topology of differential equations*, Proc. 7ECM 2018) and "independently established in
    [23]" (Rogers–Wan–Zucker, Ramanujan J. 37 (2015), arXiv:1303.2259).
  - Eq. (6), `L(f,1)/F_1(1) ∈ Q`, is introduced with "numerical calculations suggest conjectural inclusions".
    Remark 2 adds: "None of the relations in (6) seem to be proved."
  - **Normalisation.** For `r = t = 1/2`, his (5) is
    `π^{−2} Σ_{n≥0} (Γ(1/2+ε+n)/Γ(1+ε+n))^4 + O(ε^4) = π^{−2}(Y/N)(ε)`.
    So `F_0(1) = S_0`, `F_1(1) = τ_1/π^2 = D = 8S_1`, and `F_2(1) = τ_2/π^2`. Hence "8S_1 = −16L(f,1)" is
    exactly `L(f,1)/F_1(1) = −1/16`; the normalisation matches.
  - Remark 2 also mentions the "extra rational relation `F_0(1) = F_2(1)/(2π^2)`" for `r = t = 1/2`. This is
    exactly `τ_2 = 2π^2τ_0`, which Prop. 8 proves. The authors may want to point this out.
- **Bönisch–Klemm–Scheidegger–Zagier, arXiv:2203.09426 (CMP 2024).** See issue 6. §4 contains a proof, for
  this N = 8 family, that the conifold periods are combinations of `ω_f^±` (and the quasiperiods `η_F^±`),
  up to a multiplicative constant. The precise rationals are only numerical. I did not find a statement that
  proves `L(f,1)/F_1(1) = −1/16` itself.
- **Web searches** (extended) for later proofs found none. Long–Tu–Yui–Zudilin (arXiv:1705.01663) is about
  supercongruences. Bönisch–Golyshev–Klemm, arXiv:2510.03939 (2025), proves modularity of mixed periods for
  other families and does not mention this one. I found nothing from Rogers, Wan, Guo, Tauraso, Osburn–Straub,
  Dawsey–McCarthy or Allen that settles it. **Conclusion:** I found no proof of the `L(f,1)` identity with
  its constant, so the Corollary should remain conditional.

## Own computations

- I ran `v1`, `v2`, `v3`, `v5` and `v6`; all pass. I read `v1`, `v2`, `v3`, `v5` and `v6` critically.
- **Brute force.** I computed partial sums to K = 160000 and fitted the tail `a/K + b/K^2 + c/K^3`:
  - `S_0 = 1.118636387164187068…`
  - `S_1 = −0.709001367461929437…`
  - `64S_1/(1488ζ(5)S_1 + π^6S_0) = 2.451746064628782916…`, which agrees with the stated value
    2.4517460646287829170 to 18 digits.
- **PARI** (`lfunmf`, `mfinit([8,4])`): `L(f,1) = 0.35450068373096471877`, so
  `−2L(f,1) = −0.70900136746192943755`, which equals `S_1`. `L(f,2) = 0.69003116312339752512`, and
  `16L(f,2)/π^2` equals `S_0`.
- **By hand:** the Lemma 3 boundary polynomial, `θ^iν(1)`, `θ^ir(1)` and `θ^iμ(1)`, the μ_4 coefficient
  (1/12), the Casoratian, the Prop. 8 principal parts, and the §7 algebra.
