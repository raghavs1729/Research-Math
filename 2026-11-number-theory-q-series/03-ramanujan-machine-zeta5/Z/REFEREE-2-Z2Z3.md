# Referee report on `THEOREM-2.md` (Zudilin's relations (6) for (r,t) = (1/2,1/4) and (1/2,1/3))

Independent adversarial referee, 2026-11-26. I read `THEOREM-2.md` in full, took `THEOREM.md` as given (Lemma 1 (ii),
(Z1)–(Z2), and the values of `L(f_2,·)` in terms of `S_0`, `S_1`), and read Zudilin's TeX source (arXiv:1805.00544, `sigma18-086.tex`).
I formed this verdict before looking at `audit2/`.

## Verdicts

| | Theorem Z2, (1/2,1/4), f = 16.4.a.a | Theorem Z3, (1/2,1/3), f = 36.4.a.a |
|---|---|---|
| **Verdict** | **PROOF COMPLETE** (minor presentation points only, §Findings 2–7; none is a mathematical gap) | **PROOF COMPLETE** (same; Finding 1 concerns this case but does not affect it) |
| **Numerics reproduced** | **YES**: the three ratios −1/8, 1/32, −3/448 agree to 52–54 digits; my two independent summation runs agree to 1e-40, so at least 40 digits are certified | **YES**: the three ratios −1/3, 1/18, −1/90 agree to 52–54 digits; at least 40 digits certified |
| **Normalisation matches Zudilin** | **YES** | **YES** |
| **Novelty** | Not found proved in the literature. Ingredients classical/known (RWZ for the L-side; the twist trick is in RWZ, Bao and EHMM III §6) | Not found proved in the literature. The closest item is an unpublished Hsu–Tu–Yang identity, quoted in EHMM III Ex. 6.2. It concerns `L(12.4.a.a,2)` and different (non-Frobenius) 4F3 values |

I found no error in any displayed formula, constant or sign. Every step of §§1–6 was re-derived by hand, and the key
steps were re-checked numerically with my own code.

## 1. Normalisation against Zudilin (read from the TeX source)

* **Definition.** Zudilin's eq. (5) (`eq:comp`) divides by `Γ(r)Γ(1−r)Γ(t)Γ(1−t)` and has bilateral = unilateral + O(ε^4).
  It uses `z^{n+ε}`, which equals 1 at z = 1. For r = 1/2 the prefactor is `1/(πΓ(t)Γ(1−t))`. This is exactly the
  document's `F(ε)`. Eq. (6) (`LF4`) is `L(f,k)/F_k(1) ∈ Q`, k = 1, 2, 3, with the standard `L(f,s) = Σ a(n)n^{−s}`. **Match.**
* **Remark numbering.** Remark 2 (`rem1`) contains "None of the relations in (6) seem to be proved"; Remark 5 (`remB`)
  is the F_1/F_3 proportionality. The document's citations are correct.
* **Pairs and forms (Table 1).** Zudilin lists (1/2,1/3) ↔ `η_6^14/(η_2^3η_18^3) − 3η_2^3η_6^2η_18^3`, level 36, old label 36.4.1.a.
  He lists (1/2,1/4) ↔ `η_4^16/(η_2^4η_8^4)`, level 16, old label 16.4.1.a. Both agree with the document.
  PARI confirms that these eta expressions are the newforms 16.4.a.a and 36.4.a.a. Their coefficients agree for n ≤ 400,
  beyond the Sturm bounds 8 and 24, and both new spaces are 1-dimensional.
* **Rationals.** Zudilin prints **no** constants for any of the 14 weight-4 cases, so there is nothing to disagree with.
  *Curiosity, not an error:* he prints exactly `−1/8, 1/32, −3/448` (k = 1, 2, 3) in §5, but for a different example:
  the 6F5(1/2^6) family and the weight-6 level-8 form `η_2^12 + 32η_2^4η_8^8`. These coincide with Z2's constants.
  The document does not mention this. A reader who searches Zudilin's paper for "3/448" will find that example, not (1/2,1/4).

## 2. Numerics (recomputed independently; scripts in `scratchpad/zud2/`)

**L-values** (PARI/GP, `mfinit([N,4],0)`, `mfeigenbasis`, `lfunmf`; `lfuncheckfeq` = −265 and −263; root numbers +1; `\p 70`):

    L(16.4.a.a,1) = 0.8785749639883094440231224865716792749021310570418
    L(16.4.a.a,2) = 1.1136967437017574953710011825926461769631746453537
    L(16.4.a.a,3) = 1.0838984164082429558256288549336151979249837871566
    L(36.4.a.a,1) = 2.0607353967982738558069586184896744814811433976009
    L(36.4.a.a,2) = 1.4862973952817495895031380594590791018419553615920
    L(36.4.a.a,3) = 1.1299246189844927458509354836536322745171654311314

**F_j(1)** come straight from Zudilin's definition. For each n, the ε-Taylor coefficients of the summand are computed
exactly from `loggamma` and polygamma values. The partial sums are Richardson-extrapolated with two different
(block, order) choices, (16,60) and (25,48), at 120-digit working precision. The two runs agree to 1e-40. The
`hyper`/`taylor` route of the authors' scripts is not used.

    (1/2,1/4): F_0..F_3 = 1.08731032623967599238600937326066598591747691, −7.02859971190647555218497989257343419921704846,
                          35.6382957984562398518720378429646776628215887, −161.862163516964281403293909003419869556797579
    (1/2,1/3): F_0..F_3 = 1.10457125931086534557110485310724212983670894, −6.18220619039482156742087585546902344444343019,
                          26.7533531150714926110564850702634238331551965, −101.693215708604347126584193528826904706544889

(After forming my verdict I compared these with `audit2/F_values.txt`: they agree in every printed digit.)

**Ratios.**

* (1/2,1/4): −1/8, 1/32, −3/448, with |ratio − claim| = 5e-54, 2e-55, 6e-57.
* (1/2,1/3): −1/3, 1/18, −1/90, with |ratio − claim| = 2e-53, 6e-55, 2e-56.
* Control (1/2,1/2) with 8.4.a.a: −1/16, 1/32, −3/256, all to 1e-52 (THEOREM.md reproduced).
* PSLQ on `[ratio, 1]` (coefficients ≤ 1e6) returns exactly these rationals.

**Deliberately wrong variants.** In every one of these, PSLQ finds no relation with coefficients ≤ 1e6:

* the untwisted forms 8.4.a.a or 12.4.a.a in place of 16.4.a.a or 36.4.a.a;
* the two cases swapped (36.4.a.a with t = 1/4, and 16.4.a.a with t = 1/3);
* a normalisation off by π;
* an index shift `L(f,k)/F_{k−1}`.

The test therefore discriminates.

**Intermediate identities** (agreement in brackets):

* `F_1 + 2πS_0` and `F_2 + 16πS_1` [1e-51].
* `F_1 = −(32/π)L(8.4.a.a,2)`, `F_2 = 32πL(8.4.a.a,1)` [1e-52].
* `F_1 = −(12√3/π)L(12.4.a.a,2)`, `F_2 = 8√3πL(12.4.a.a,1)` [1e-52]. PARI confirms g = 12.4.a.a.
* `F_3 = (csc^2πt + 1/3)π^2F_1` for t = 1/4, 1/3, 1/2 [1e-51].
* All six relations of Prop. C, i.e. the twist identities between PARI L-values of levels 16/8 and 36/12 [1e-60].
* Coefficient identities 16.4.a.a = 8.4.a.a ⊗ χ_{−4} and 36.4.a.a = 12.4.a.a ⊗ χ_{−3} for n ≤ 400.
* My own evaluation of `U_0, U_1, U_2` from the defining series, against Prop. A(ii) at y = 1.5, 2, 3 and Prop. B(ii)
  at y = 0.7, 1, 1.6 [1e-29].
* Cor. A and the §5 integral formula for j = 0, 1, 2, including both F_0 "incomplete period" formulas [1e-23 … 1e-31].

## 3. Line-by-line check of the proof

* **Definition / Lemma 1.**
  * `c(s)u_t(s)` is the summand. `c(s) = (1/π)B(s+½,½)` holds for Re s > −½.
  * `|u_t(n+ε)| ≤ C(n+1)^{−3/2}` is uniform on |ε| ≤ 1/8, because each factor is Γ(s+a)/Γ(s+1) with exponents summing to −3/2, and 1/8 < t.
  * The dominating function `C'x^{−1/8}` is integrable against dx/√(x(1−x)).
  * Fubini and Cauchy then give the claim. **Correct.**
* **Lemma 2.**
  * The recurrence gives `D_tU = ε^3u_t(ε)x^ε`.
  * The ψ-values check out: ψ(¼)+ψ(¾) = 2ψ(1) − 6log2; ψ(⅓)+ψ(⅔) = 2ψ(1) − 3log3; ψ'(½) = 3ψ'(1); ψ'(t)+ψ'(1−t) = π²csc²πt.
  * This gives M = 256, 108 and c_t = π², 2π²/3.
  * The uniqueness claim follows from the triple exponent 0 at x = 0 (local basis 1, log, log²). **Correct.**
* **Lemma 3** (bilateral sum).
  * I re-derived the reflection identity for n = −m.
  * Decay O(n^{−2}) on both sides; the order-4 zero of the n < 0 part at u = 0; the pole structure (double pole at ½,
    simple poles at ±t); Liouville; K = 0 and B + C_1 + C_2 = 0.
  * I re-derived the cot Taylor coefficients by hand (−S, SC, −S(3S−2)/3).
  * F_1 = −π²(C_1+C_2)(S−1) and F_3 = −(π⁴/3)(C_1+C_2)(S−1)(3S+1), hence F_3/F_1 = π²(S + 1/3). **Correct.**
  * The uniform O(1/|Im u|) bound is only sketched (Finding 3).
* **§2 classical facts.**
  * (J1)–(J3): dk/dy agrees with THEOREM.md (Z2); duplication; Landen; k'(1/y) = k(y).
  * (B1)–(B2): I checked b(−1/(3τ)) = −i√3τc(τ) by hand, and h(1/(3y)) = 1 − h(y) follows.
  * (T1) Clausen, (T2) Pfaff, (T3) Whipple: I re-ran `z2_symbolic.py`. S1–S3 are genuine polynomial-identity proofs in
    Q(h,t)[w, w'] with w'' eliminated. Valid proofs.
* **Prop. A / Cor. A** (t = 1/4).
  * P, z and X; `1 − z = P²/θ_4^8`; dz/dy = −πzP; z(√2) = −1 via Landen + Jacobi; dX/dz.
  * √(X(1−X)) = 2√(−z)(1+z)/(1−z)²; (iii).
  * The chain λ = k² → λ_− = −θ_2^4/θ_4^4 (T2) → z = 4λ_−(1−λ_−) (T1) → X (T3): the domains are monotone and regular,
    and dimension 3 on both sides gives the bijection.
  * Asymptotic matching via Lemma 2.
  * θ_2²θ_3²θ_4^4(σ) = 4f_2(σ/4); y = 4t. **All correct.**
* **§4.**
  * The folding identities under t ↦ 1/(8t).
  * F_1 = −128πρ_1 = −(32/π)L(f_2,2) and F_2 = 64π²ρ_0 = 32πL(f_2,1). **Correct.**
* **Prop. B / §5** (t = 1/3).
  * Abel's Wronskian, C = −√3/(2π) (this matches the known Wronskian of 2F1(⅓,⅔;1;h) and 2F1(⅓,⅔;1;1−h)).
  * h(1/√3) = ½; dx/dy; √(h(1−h)) = 3√3η^3η_3^3/a³.
  * g_0 = η_1^6η_3^2 + 9η_1^3η_3^2η_9^3 via a(q) = b(q) + 3c(q³).
  * Fricke: g_0(−1/(3τ)) = 9τ⁴g_0(τ), checked by hand.
  * ρ_2 = ρ_0/3; F_1 = −(12√3/π)L(g,2), F_2 = 8√3πL(g,1). **Correct.**
* **Lemma 4.**
  * Conjugations: τ/(4τ+1) = W_8T^{−1/2}W_8 and τ/(3τ+1) = W_3T^{−1}W_3, both checked by hand. The signs −(4τ+1)^4 and
    −(3τ+1)^4 follow from f_2(τ+½) = −f_2 and g_0(τ+1) = −g_0.
  * Twist identities: Sturm bound 8 for Γ_0(16) and 24 for Γ_0(36), with 400 coefficients checked. That is enough even
    with Shimura's crude level Nm² (Γ_0(128), Γ_0(108), Sturm 64, 72).
  * Eta-quotient conditions: I checked the Ligozat congruences and the square condition by hand for η_4^16η_2^{−4}η_8^{−4}
    at level 16, and for η_6^14η_2^{−3}η_18^{−3}, η_2^3η_6^2η_18^3 and η_2^6η_6^2 at level 36. All hold, with trivial
    character. **Correct** (but see Finding 1).
* **Prop. C** (Manin's trick).
  * χ_{−3}(n) = (ζ^n − ζ^{−n})/(i√3) and χ_{−4}(n) = (i^n − i^{−n})/(2i).
  * f(it) = (2/√3)Im g_0(⅓+2it), resp. Im f_2(¼+it).
  * γ = (1 0; 3 1), resp. (1 0; 4 1): w − 1/N = −1/(N(Nτ+1)).
  * I recomputed all six contour sums (j = 0, 1, 2 for both cases):
    * (b): `−iρ_0 − 6ρ_1`, `−(2/3)iρ_0`, `−iρ_0/9 + (2/3)ρ_1`;
    * (a): `−8ρ_1`, `−iρ_0/2`, `ρ_1/2`.
    Hence Λ_f(1), Λ_f(2), Λ_f(3) = 2√3ρ_1, ρ_0/(3√3), ρ_1/(6√3), resp. 8ρ_1, ρ_0/2, ρ_1/2. **Correct.**
  * No Atkin–Lehner sign of f itself enters: only the Fricke eigenvalue +1 of f_2 (W_8) and of g_0 (W_3; equivalently
    g|W_12 = +g, since g(−1/(12τ)) = 144τ⁴g(τ)). These signs are verified analytically and numerically. They are
    consistent with root number +1 for 16.4.a.a and 36.4.a.a (PARI).
  * (a) reproduces RWZ (30)–(31).
* **End of proof.** The rational arithmetic is re-done and correct. Non-vanishing holds because f_2(it) > 0 and g_0(iy) > 0.
  The restatements `F_3 = −(448/3)L(f,3)`, `F_2 = 18L(f,2)`, `F_1 = −3L(f,1)`, `F_3 = −90L(f,3)` are correct.

## 4. Findings (by severity; no major or fatal findings)

1. **(Minor, justification error in a side check.)** Location: Lemma 4 proof and `verify/z2_symbolic.py` (S5).
   * S5 runs Ligozat's cusp-order formula for `η_2^6η_6^2` "at level 12". That quotient is not on Γ_0(12):
     Σ(12/δ)r_δ = 40 ≢ 0 mod 24. The output contains the non-integral "orders" 5/3, which signals the misuse.
   * S5 also never checks Ligozat's congruence and character conditions, only the orders.
   * The claim that g has level 12 is true (g = a(q²)η_2^3η_6^3 with a ∈ M_1(Γ_0(3),χ_{−3}) and η_2^3η_6^3 ∈
     S_3(Γ_0(12),χ_{−3}); PARI identifies g as 12.4.a.a). But the document never proves it.
   * The proof does not need it. Both summands of g are eta quotients on Γ_0(36) with trivial character, so g⊗χ_{−3} ∈
     M_4(Γ_0(36)). And even with the crude Nm² level, 400 coefficients exceed the Sturm bound.
   * Fix: replace "level 12" by "level 36 termwise" in the Lemma 4 argument, or add the a(q²)η_2^3η_6^3 argument; move
     that S5 line to level 36.
2. **(Minor, presentation.)** Prop. B(i) writes the Wronskian as `√3a² dy/dh` before y ↦ h is known to be invertible.
   This is not circular: the chain rule gives `W(h(y))·h'(y) = √3a²` directly, hence h' = −2πh(1−h)a² < 0. One line to fix.
3. **(Minor, sketched estimate.)** Lemma 3 asserts `T(u) = O(1/|Im u|)` uniformly in Re u ("each term is O(|u+n|^{−2})").
   This needs Stirling's estimate uniformly in |arg s| ≤ π − δ for the terms with Re(u+n) near 0, and the reflected form
   for n < 0. It is standard and true, and the conclusion is confirmed to 1e-51, but it is not written out.
4. **(Minor, sketched estimate.)** Prop. C applies Cauchy's theorem on the ideal triangle (0, 1/N, i∞) with "the integrand
   decays exponentially at the three cusps". This is correct because f_2 and g_0 are cusp forms for congruence subgroups:
   in the local coordinate at 0, the triangle becomes a vertical half-strip of width 4, resp. 3. The closing arcs are not
   discussed.
5. **(Wording.)** §0: "F_0 is *not* an L-value" is not proved. Only the incomplete-period formulas for F_0 are proved
   (both confirmed numerically). The sentence should be marked heuristic (cf. Zudilin, Remark 2).
6. **(Typo-level.)** Lemma 4 proof: "τ ↦ τ/(Nτ+1) is W_N T^{−1} W_N^{−1}, resp. W_8 T^{−1/2} W_8^{−1}". Here N = 3 in the
   first expression, and the "resp." order is reversed relative to (a), (b). Both conjugations are correct.
7. **(Information for readers.)** The coincidence of Z2's constants with Zudilin's printed weight-6 constants (§1 above)
   deserves a footnote, so that nobody "checks" Z2 against the wrong example in Zudilin's paper.

## 5. Novelty

* **Searches (2026-11-26).** Four web searches: Zudilin (6) proofs for 16.4.a.a and 36.4.a.a; conifold/4F3(½,½,¼,¾) L-values;
  4F3(½,½,⅓,⅔) and level 36; 2024–26 proofs of (6).
* **Sources read.**
  * EHMM III (Grove–Long–Rosen–Tu, arXiv:2607.25173), §6.
  * Dawsey–McCarthy (arXiv:2101.06303). Its statement that "all the relations in Table 1 should now be known" refers to
    the *finite-field* relations, via LTYZ supercongruences, not to (6).
  * The abstract of Grove (arXiv:2507.19971), which concerns 3F2 / weight 3.
* **Result.** No proof of (6) for (1/2,1/4) or (1/2,1/3) was found. This agrees with `LITERATURE-13.md`.
* **Closest prior work.**
  * EHMM III Ex. 6.1 proves `L(16.4.a.a,1) = (π/4)·4F3(½^4;1^3;1)` by an Euler integral of the symmetric square and the
    shift `f_8(τ+¼) = −i f_16(τ)`. This is the same twist mechanism as Prop. C(a), applied to an ε^0 value.
  * EHMM III Ex. 6.2 quotes unpublished Hsu–Tu–Yang: `L(12.4.a.a,2)` as a combination of two non-Frobenius 4F3(1)'s,
    obtained from the Euler integral of `F(κ_3(3),t)` against `g(τ/2)`. That is the same parametrisation idea as §5,
    applied to different hypergeometric values.
* **Assessment.** The method's ingredients are known, and the L-side for (1/2,1/4) is RWZ (30)–(31). The identities
  `F_1, F_2 ↔ L(f_2,·), L(g,·)` for the ε^1 and ε^2 Frobenius coefficients, and hence (6) in these two cases, appear
  unrecorded. Describe them as "new assembly of known techniques, not found in the literature".

## 6. What I could not check

* Bönisch's Bonn thesis (2023) and the Hsu–Tu–Yang preprint, which is "in preparation" and not available.
* BKSZ and BGK in full; I relied on `LITERATURE-13.md` for them.
* A complete citation list of SIGMA 14 (2018) 086 (MathSciNet / Google Scholar not accessible).
* The classical facts (B1)–(B2) (Borwein cubic inversion, a(−1/(3τ)) = −i√3τa(τ)) and Jacobi's theta identities are
  cited, not re-proved. I confirmed them numerically only, at sample points, to 1e-29.

**Reproducibility** (scratch, `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/zud2/`):

* `lv.gp`, `lv.out`: PARI L-values and newform identification.
* `fj_indep.py`, `fj_indep.out`, `Fvals.pkl`: F_j(1).
* `ratios.py`: ratios, intermediate identities and wrong variants.
* `steps.py`, `steps2.py`: Props A/B, Cor. A, §5, F_0.
