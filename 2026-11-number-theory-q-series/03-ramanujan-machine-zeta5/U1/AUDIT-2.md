# AUDIT-2 of THEOREM.md (U1): referee proof-read and novelty search

Date: 2026-11-26. Scope: (1) a line-by-line logical check of §1–§5 of `THEOREM.md`; (2) a search for any earlier
identification of the constant. Numerical agreement (about 498 digits, audit scripts a1–a5) was taken as given and
not rerun. The only new computations were the asymptotic constant in issue 1 and the first terms of `64^n w_n`.

## Verdict

**COMPLETE AFTER LISTED FIXES.** The main identity `PCF = 96/(372 ζ(5) − π^5)` is proved correctly. I found no
fatal error, and every sign, index and constant I recomputed by hand is right. One claim in the theorem box, the
rate `p_n/q_n − PCF ~ c/n`, is not proved in the text (issue 1). It can be proved in a few lines, and a proof is
sketched below. Everything else is cosmetic or a one-line gap.

## Part 1: proof-read

Results of the checks you asked for (everything here was recomputed by hand):

* **(a) Lemma A.** This is a legitimate proof step. `v1_operator_identity.py` works in the differential ring
  `Q(x)[s, u, u', h_0..h_5]` with `s' = −1/(2s)` and `s^2 = 1 − x`, uses the hypergeometric ODE to replace `u''`,
  and sets `D h_k = h_{k+1} T'`. Here `u`, `u'` and the `h_k` are free symbols, so the identity is proved formally,
  which is stronger than needed and therefore valid for the actual functions. The branch of `√(1−x)` is the same
  one used in `T'`, so either branch works. `h` is arbitrary because it enters only through `h^{(k)}(T)`.
  Independent check: the leading symbol of `L` is `(1−x)^2 θ^5`. Applying it to `u^4 h(T)` with
  `T' = 1/(x√(1−x)u^2)` gives `h^{(5)}/(√(1−x)u^6)`, which matches (this is Bol's identity for a symmetric fourth
  power). `T'` is `W/u^2`, where `W = 1/(x√(1−x))` is the hypergeometric Wronskian. The expansion
  `ã(n) = a_n/8 = 2T^5 + (5/4)T^3 + T/16` and the coefficient action of `L` are both correct.
* **(b) §3.** For `t > 1`: `λ < 1/2`, so (C6) applies and `√(1−x) = 1 − 2λ > 0`. `dx/dt = −π x√(1−x)u^2` and
  `T' = 1/(x√(1−x)u^2)` are correct. `F/((1−2λ)θ_3^{12}) = λ(1−λ)` holds because
  `(θ_4^4 − θ_2^4)/θ_3^4 = 1 − 2λ`. Hence `LQ̃ = x/8`.
  The nome expansion is correct: `λ = 16q − 128q^2 + …` and `x = 4λ(1−λ) = 64q + …`. `θ_3^8 Φ = 8q + …`, so
  `Q̃ = x/8 + O(x^2)`, giving `v_0 = 0` and `v_1 = 1/8`. The coefficient equations of `LQ̃ = x/8` are:
  `n = 0` gives `0 = 0`; `n = 1` gives `v_1 − ã(0)v_0 = 1/8`, which is consistent; `n ≥ 2` is homogeneous.
  For `P`: `w_1 = ã(0) = 1/4 = 4u_1` and `a_0 = 2 = 8w_1`. `LP = 0` on `(0,1)` gives the recurrence because `LP` is
  a convergent power series on `|x| < 1`.
* **(c) Proposition 1.** The base cases and the induction step `U_n = a_n U_{n−1} − 64n^{10} U_{n−2}` are correct.
  The step at index `n ≥ 1` uses the coefficient recurrence at `n + 1 ≥ 2`, so it is valid for both `w` and `v`.
* **(d) §4.**
  * The Casoratian recursion `(n+1)^5 D_n = n^5 D_{n−1}` is correct, with `D_0 = 1/8`. The sign of (4.1) is
    positive.
  * Lower bound for `w_n`: Bernoulli's inequality `(1 − 1/j)^{3/4} ≤ 1 − 3/(4j)` holds for `j ≥ 2`. For `n = 1`
    the product is empty, so the bound is `1/4` with equality. `w_n ≥ u_n` holds because `u^3` has non-negative
    coefficients and constant term 1. The bound `ℓ ≤ 1/2 + 32(ζ(2) − 1)` is a little sharper than the stated one,
    and the stated one is valid.
  * `0 < v_n ≤ ℓ w_n`, which gives radius of convergence ≥ 1.
  * `Q̃` is real-analytic on `(0,1)`: `dx/dt ≠ 0`, so `t(x)` is real-analytic by the inverse function theorem, and
    the identity theorem applies.
  * The Abelian lemma is applied with `α_n = (n+1)w_{n+1} > 0` and `β_n/α_n = v_{n+1}/w_{n+1} → ℓ`. It needs
    `Σα_n x^n → ∞`, which holds because `n w_n ≥ n^{−1/2}/16`.
  * The limit `x → 1−` is the same as `t → 1+`, by continuity of the inverse of `λ`. `g`, `φ`, `g_t` and `φ_t` are
    analytic at `t = 1`. `g(1/t) = t^4 g(t)` gives `g_t(1) = −2g(1)`. (4.3) is correct.
* **(e) §5.**
  * The Laplace integral evaluates to `12e^{−a}(2 + a)/a^5`.
  * `|a°_n| ≤ (16 + 1040 + 1024)σ_5(n) = 2080 σ_5(n)`, so termwise integration is justified by dominated
    convergence.
  * `(t−1)^3(t+1) = t^4 − 2t^3 + 2t − 1` is correct.
  * The substitution `t → 1/t` together with `F(i/t) = −t^6F(it)` is correct.
  * The Dirichlet series `16 ζ(s)ζ(s−5)(1 − 65·2^{−s} + 64·4^{−s})` is correct.
  * Poles of the right side. They do not matter, since `M` is entire and `s = 4, 5` are regular points. For the
    record they all cancel:
    * `s = 6`: the factor is `1 − 65/64 + 1/64 = 0`.
    * `s = 1`: `ζ(−4) = 0`.
    * `s = 0`: the factor is `1 − 65 + 64 = 0`.
    * `s = −k`: a trivial zero of `ζ(s)` or `ζ(s−5)`.
  * Arithmetic: `M(5) = 186ζ(5)/π^5`, `M(4) = 1/4` and `ℓ = (372ζ(5) − π^5)/96` are all correct.
* **(f) Lemma B.** I re-derived it from (C5):
  * `E_6(q^4)` uses the thetas at `2τ`. The only irrationality, `r = θ_3^2θ_4^2`, appears only as `r^2 = AB`.
  * `E_6(q)` uses the thetas at `τ/2`, where `r' = θ_2^2θ_3^2` appears only as `r'^2 = AC`.
  * Both closed forms in the text are correct.
  * The coefficients `63·(16, −1040, 1024) = (1008, −65520, 64512)` match `−2E_6(q) + 130E_6(q^2) − 128E_6(q^4)`,
    and the constant terms cancel (`−2 + 130 − 128 = 0`).
  * The transformation `F(i/t) = −t^6F(it)` is correct.
* **(g) What the statement claims beyond the proof.** See issues 1 and 2.

### Numbered issues

1. **Theorem box, "p_n/q_n − PCF ~ c/n (c ≈ 0.724)".** Severity: *gap, fillable in a short paragraph.* It is stated
   as part of the theorem, but only Remark 2 supports it, with a heuristic. Proof sketch:
   * `p_n/q_n − 1/ℓ = (ℓ − r_{n+1})/(ℓ r_{n+1})`, where `r_m = v_m/w_m`.
   * The tail is `ℓ − r_m = Σ_{k≥m} 1/(8(k+1)^5 w_k w_{k+1})`.
   * `u` continues analytically to `C \ [1,∞)` and equals `u(1) + B_0(1−x)^{1/2} + O(1−x)` at `x = 1`, with
     `u(1) = Γ(1/4)^2/(2π^{3/2})` and `B_0 = Γ(−1/2)/Γ(1/4)^2`.
   * Flajolet–Odlyzko transfer then gives `w_n ~ C n^{−3/2}` with `C = Γ(1/4)^4/(2π^{9/2}) = 0.50040…`. This
     agrees with Kotešovec's asymptotic in OEIS A397472.
   * Therefore `p_n/q_n − PCF ~ c/n` with `c = 1/(8C^2ℓ^2) = 4608 π^9 / (Γ(1/4)^8 (372ζ(5) − π^5)^2) = 0.723945…`.
   * Fix: add this argument, or move the rate claim into Remark 2.
2. **Theorem box, the binomial form of `w_n` via Clausen.** Severity: *cosmetic.* The text says this is "not used in
   the proof", but it is stated in the theorem. It is a direct consequence of Clausen's formula
   `2F1(a,b;a+b+1/2;x)^2 = 3F2(2a,2b,a+b;2a+2b,a+b+1/2;x)` with `a = b = 1/4`, together with
   `(1/2)_k/k! = C(2k,k)/4^k`. Fix: state this and give a precise citation for Clausen's formula.
3. **§3, "Q̃ is analytic at x = 0".** Severity: *one line.* Add that `λ(q)` is locally biholomorphic at `q = 0`
   (`λ'(0) = 16`), so the real nome `q = e^{−πt(x)}` on a small interval `(0,ε)` is the analytic inverse branch.
   Then `Q̃` on `(0,ε)` is the restriction of the analytic germ.
4. **§4, "Q'/P' = (gφ)_t/g_t".** Severity: *cosmetic.* The division also needs `g_t ≠ 0` on `(1,∞)`. This follows
   from `P'(x) > 0` and `dx/dt ≠ 0`. Say so.
5. **§4, "Q = Q̃ on (0,1) (both real-analytic)".** Severity: *cosmetic.* Add: `t(x)` is real-analytic by the inverse
   function theorem, since `dx/dt = −πx√(1−x)u^2 ≠ 0`.
6. **§5, "For Re s > 6 … −∫_1^∞ t^{5−s}F = ∫_0^1 t^{s−1}F".** Severity: *cosmetic.* `∫_0^1` converges for every `s`,
   because `F(it) = O(t^{−6}e^{−π/t})` as `t → 0`. `Re s > 6` is needed only for the Dirichlet-series
   interchange. Optionally add the pole-cancellation remark from (e).
7. **§1 (C5) and (C2)–(C3), references.** Severity: *cosmetic.* Conventions such as nome `e^{iπτ}` versus `e^{2πiτ}`
   and the normalization of `E_6` are easy to get wrong. Give an exact equation number for
   `E_6(τ) = (1/2)(θ_2^4+θ_3^4)(θ_3^4+θ_4^4)(θ_4^4−θ_2^4)` with nome `q = e^{iπτ}`. I checked it to `O(q^2)` by hand,
   and `v2` checks it to `O(q^300)`.
   * The finite q-series check in v2 is only a confirmation. The proof rests on (C5) plus the exact polynomial
     identity, and that logic is sound.
   * Lemma B is also standard: `F` is the weight-6 Eisenstein form for Γ_0(4).
8. **§0, "1.20426..." versus 1.2042534.** Severity: *cosmetic, but worth one sentence.* The printed value is not a
   truncation of the true limit. This is consistent with the slow `0.724/n` convergence: convergents print as
   1.20426 for depths of roughly 4·10^4 to 1.1·10^5, and `a5_source_value.py` tabulates this. State it, so a reader
   does not suspect a different row.
9. **§0, "Openness" paragraph.** Severity: *cosmetic.* Update the list of sources searched with the ones in Part 2,
   especially OEIS A397472 and Cohen arXiv:2607.06581.
10. **Remark 4(iii), the `K(z)` integral formula.** Severity: *cosmetic.* It is marked "verified to 60 digits", which
    is fine because it is a remark, not part of the theorem.

No fatal issues.

## Part 2: novelty

**Conclusion: I found no earlier identification of this constant, of the continued fraction, or of the Apéry limit
of the recurrence.** One near-miss narrows the novelty claim: OEIS **A397472** (A. Shvets, created June 2026) already
records the numerator sequence `64^n w_n`, namely 1, 16, 496, 19456, 860656, 40950016, …. It gives the sequence as
`64^n[z^n] 2F1(1/4,1/4;1;z)^4`, the binomial-convolution formula, the order-5 recurrence (the same one as here) and
the asymptotics. It does not mention the companion solution, the Apéry limit, ζ(5), π^5 or the Ramanujan Machine. So
"`p_{n−1} = 8^n n!^5 w_n` with the binomial form" is not new as a description of the sequence. What is new is the
identification of the continued fraction with it, and the evaluation of the limit.

Sources checked (R = read or grepped in full; S = abstract or snippet only):

| # | Source | How | Finding |
|---|---|---|---|
| 1 | arXiv:2308.11829 v1 and v2 (v3 does not exist; the API lists v2 as the last version, journal ref PNAS 121(25) e2321440121) | R (pdftotext) | The row appears in Table 1 and in Table 4, row 5, in both versions, only as "1.20426...". Table 4 says rows 5–7 are "only formula candidates as we did not find a closed-form expression". |
| 2 | PNAS version and SI (pnas.org, PMC11194572) | not reachable (403 / captcha / Europe PMC 503) | Not checked. arXiv v2 carries the PNAS journal reference. |
| 3 | Cohen, arXiv:2607.06581 (Oct 2026 version, encyclopedic dictionary of polynomial continued fractions), full TeX | R | 25 ζ(5) lines, all for `b_n = −n^10` or `−16n^10`, combining ζ(5), ζ(3) and log 2. No `π^5`/ζ(5) entry, no `−64n^10` entry, no 1.2042. |
| 4 | Cohen, arXiv:2409.06086 TeX (CFfile2) | R | Same ζ(5) entries. Nothing relevant. |
| 5 | Yamamoto, arXiv:2403.09729 TeX | R | Proves only the three `−n^10` rows of Table 1. |
| 6 | RM group: 2412.12361, 2412.16818, 2502.17533, 2507.08138, 2610.13919, 2610.27799, 2607.09721 (Ramanujan Challenge), 2610.40551 (Challenge results), 2607.27794 | R (grep) | No 1.2042, no `π^5`, no `−64n^10` ζ(5) formula. The Challenge problems do not include this row. |
| 7 | David, arXiv:2303.09318; Brier–Naccache–Yifrach-Stav, arXiv:2211.01058; 2212.09470; 2601.08461; 2602.03027 | R (grep) | Nothing. |
| 8 | Shvets: 2604.06239, 2604.14219 (Eichler integrals, level-8 Apéry limit), 2604.06238, 2604.09723, 2605.19773, 2606.15462, 2607.19427, 2611.01481 | R (grep) | The same method family (eta products and Eichler integrals for Apéry limits), but for other sequences. No ζ(5), no `π^5`. |
| 9 | OEIS A397472, A397369, A395287 | R (JSON entries) | Sequence, recurrence and asymptotics present; no Apéry limit (see above). |
| 10 | Shvets, Zenodo 10.5281/zenodo.20940394 (linked from A397472) | S (record page) | Supercongruences only. |
| 11 | L. Gao, Zenodo 20527727; ramanujanmachine.com/?p=1291 (forum proof, June 2026) | S | `b_n = −n^{2m}` families only. No `−64n^10`, no `π^5`. |
| 12 | Chamberland–Straub, arXiv:2011.03400; Golyshev, arXiv:0709.1968; arXiv:2109.05359 | R (grep) | Nothing. |
| 13 | Almkvist, math/0703261 (5th-order CY-related equations); AESZ tables math/0507430; Almkvist–Zudilin math/0402386; 0808.1480; 0902.4786 | R (grep for the operator `8θ^4+16θ^3+17θ^2`, ζ(5), π^5) | Our operator is not listed. |
| 14 | github.com/openai/math CONTENTS.md (raw) | R (grep) | Nothing on ζ(5) continued fractions or Ramanujan Machine constants. |
| 15 | github.com/RamanujanMachine (code search `64*n**10`, LIReC analyze_pcfs.py) | R (one file) | Only `−n^10`-type pairs. |
| 16 | Semantic Scholar list of papers citing 2308.11829 (10 entries) | R (list) | All relevant entries are covered in rows 6 and 13 above. |
| 17 | Web searches: "372 ζ(5)", "31ζ(5)/8 − π^5/96", "1488ζ(5)", "1.20426", "−64n^10 zeta(5)", 2F1(1/4,1/4;1)^4 Apéry limit, Eisenstein/Eichler with ζ(5) and π^5 | S | No hit. Zeilberger's oAperyLimits2 page has no ζ(5) or π^5 match. |

What could still make "new" false: the PNAS SI, which I could not reach (it most likely mirrors the arXiv appendix);
non-indexed preprints such as Zenodo or ResearchGate beyond those found; and private communications. Nothing found
weakens the claim beyond the A397472 point above.
