# Audit of `monotonicity/akbary-francis/RESULT.md` (Akbary–Francis, q = 13)

Auditor: independent referee, 2026-11-24. Nothing under review was edited. Audit code and logs are in
`audit/` (`audit/logs/`). `audit/SHA256SUMS_audited` pins the audited files (RESULT.md, the zero files,
certificates and the main scripts) as they were when audited. Nothing was committed. All runs used `nice -n 10`
and at most 2 processes.

Primary source, read by me: arXiv:1810.11524v1, fetched again from arXiv. It is byte-identical to `source/af.pdf`
and to the copy on the first author's web page. I read the **TeX source** (arXiv e-print), so every conjugation bar
and sign below comes from the LaTeX itself, not from pdftotext. I read §1 in full, §3, §4 (Thms 4.2, 4.3 =
`omegaJ`, `omegaJ2`), §5 (proofs of Thm 1.14 and 1.16), §6 (Thm 6.1, Lemma 6.2, Cor. 6.3, Lemma 6.4, proof of
Thm 1.18), §7 (computation of F_q) and §8 (proof of Thm 1.22).
Journal version: **Math. Comp. 89 (2020), no. 322, 993–1026, DOI 10.1090/mcom/3463** (title "…in a fixed
arithmetic progression"). Its abstract (Crossref) states the same results: q ≤ 10 and q = 12, 14; infinitely many
k for q ≤ 400 000; non-squares. **I could not read the full text** (the AMS returns 403). See G1.

## 0. Verdict

* **No RIGOR BREAK found.** Theorem A, Theorem A′, Corollary B and the explicit q = 17 failure all survive
  independent re-derivation and re-computation, with different code, different character constructions and a
  different zero-certification method.
* **Two GAPs, both about wording and novelty, not mathematics:**
  * **G1.** The journal version was not consulted. The literature search in RESULT.md was not done.
  * **G2.** The statement of Theorem A cites an "announced Lean-checked proof" of GRH. This must come out of
    every theorem statement.
* MINOR items are in §9. One of them is an *improvement*: for q = 17 and q = 19 the "infinitely many k" statement
  also holds unconditionally, by the argument of Corollary B.

## 1. Conventions, re-derived from the TeX (task item 1)

Everything RESULT.md §1 states agrees with the source:

* **(3)** (`GeneralIneq`, eq. (3)): n/(φ(n)(log(φ(q) log n))^{1/φ(q)}) > 1/C(q,a). The exponent 1/φ(q) applies to
  log(φ(q) log n), with a single log inside.
* **f(x;q,a)** = (log(φ(q)θ(x;q,a)))^{1/φ(q)} C(q,a)^{−1} ∏_{p̄≤x}(1−1/p̄). (3) holds for N̄_k ⇔ log f < 0 on
  [p̄_k, p̄_{k+1}) (eq. (4)).
* **Thm 1.8:** C(q,a)^{φ(q)} = e^{−γ} ∏_p (1−1/p)^{α(p;q,a)}, with α = φ(q)−1 on the class and −1 off it.
* **Thm 1.18** (`oscillation`, under SH_{q,a}, m = Ind_q(a)):
  φ(q) log f = (√x log x)^{−1} Σ_χ **χ̄(a)** Σ_{ρ∈Z(χ′)} x^{i Im ρ}/(ρ(ρ−1)) − m R_{q,a}/((m−1)x^{(m−1)/m} log x),
  plus O(1/(√x log² x)) (the O-term sits outside the bracket, so it is not divided by φ(q)). The bar is
  `\overline{\chi}(a)` in the TeX, and the denominator is ρ(ρ−1).
  For ρ = 1/2+iγ, ρ(ρ−1) = −(1/4+γ²). For m = 2 this gives
  φ√x log x · log f = −S_a(log x) − 2R + O(1/log x), with S_a(u) = Σ χ̄(a) Σ e^{iγu}/(1/4+γ²).
  So **(3) fails near x ⇔ S_a(log x) < −2R** (up to O(1/log x)). This matches RESULT (★).
* **Cor. 1.19** (`chapter6`): x₀ exists ⇔ lim sup of (7) < 2R_{q,a}. **Thm 1.22** is stated for q ≤ 10 and
  q = 12, 14 only, "for all positive integers k".
* **Derivation check.** I traced Lemma 6.2 → Cor. 6.3 → Prop. 3.6. With Nicolas's Lemma 3.4,
  F_ρ(x) = −x^{ρ−1}/((ρ−1)log x) + r_ρ, so −F_ρ/ρ = x^{ρ−1}/(ρ(ρ−1) log x). This is the sign in (6), and the
  −mR term comes from Prop. 3.5. McCurley's formula carries χ̄(a), as in the standard explicit formula.
* **Numerical test, independent of RESULT's sieve** (`audit/psieve.c`, `audit/compare_G.py`). I sieved all primes
  ≡ 1 mod 13 up to 10⁹, with my own C(13,1). At 88 points u = log x ∈ [12, 20.7] I compared
  G = 12√x log x · log f with −S_1 − 4 (zeros below 400):
  * u ≥ 18: rms difference 0.075 and correlation 0.99996, while G ranges over [−5.69, −1.40].
  * With the wrong sign (+S − 4) the rms is 2.60 and the correlation −1.000.
  * With threshold 0 instead of 4 the rms is 3.93.

  So sign, normalisation and the constant 2R = 4 are right.

## 2. The exact-tail argument and the zero list (task item 2)

* **Logic.** Under GRH_13 every zero is 1/2+iγ with weight 1/(1/4+γ²) > 0, and the weights of all zeros, counted
  with multiplicity, sum to F_13. Then S(u₀) ≤ S_T(u₀) + (F_13 − W), provided two things hold:
  * every listed window holds a zero of the labelled function;
  * listed zeros are distinct elements of the multiset. Windows of one label are disjoint, windows of different
    labels belong to different factors, and the partner zeros have negative ordinate.

  Completeness is **not** used. A zero at s = 1/2, allowed by GRH, would sit in the tail with |cos| ≤ 1, so it is
  covered. The argument is correct.
* **F_13 closed form.** I re-derived it from Λ′/Λ(1, χ): F(χ) = 2 Re Λ′/Λ(1,χ) = log(q/π) + 2 Re L′/L(1,χ) +
  ψ((1+κ)/2). This is the paper's (50) (−γ for odd χ, −γ − 2 log 2 for even χ), and F_1 = 2 + γ − log 4π.
  I evaluated it in mpmath with my own characters and L, L′ at 1 from Hurwitz ζ (`audit/mp_check.py`):
  **F_13 = 6.45484563118166436057782527964**, equal to RESULT in all 30 digits. Paper: 6.45484.
* **W and the tail.** In mpmath, with ordinates as exact decimals: **W = 6.38320081487160952753045515817**,
  **tail = 0.07164481631005483304737**. Both are identical to `certificate_q13_a1.json`.
* **Independent re-certification of every window** (`audit/recert_zeros.py`). It does not use `hardy_z`,
  `dirichlet_char`, or any audited script:
  * Characters are built by hand from the least primitive root.
  * L(s,χ) = q^{−s} Σ χ(m) ζ(s, m/q), using Arb Hurwitz ζ as power series.
  * A complex **Krawczyk test** on the box [1/2 ± 10⁻³⁰] + i[d ± 10⁻³⁰] proves existence and uniqueness of a zero
    in each window.
  * Windows of one label are checked to be disjoint.

  Results:

  | file | windows | failures | max \|d − γ\| (refined) | label convention |
  |---|---|---|---|---|
  | zeros_q13 (sha 98115f…) | 4208 | 0 | 2.6e−42 | Conrey, all 11 labels |
  | zeros_q11 (sha eea2f7…) | 3385 | 0 | 3.7e−42 | Conrey |
  | zeros_q17 (sha 217bf7…) | 1301 | 0 | 5.0e−43 | Conrey |
  | zeros_q19 (sha 221c69…) | 1476 | 0 | 7.6e−43 | Conrey |

  A deliberately wrong entry (`6 0.883961`) is rejected, which is the negative control. The ordinates are good to
  about 10⁻⁴², far better than the claimed 10⁻³⁰.
* **Hash binding.** `afcommon.load_zeros` refuses any file whose hash differs from `SHA256SUMS`, and the
  certificate JSON records the hash. My logs record the same hashes, so the certificate is bound to the file I
  re-certified. This fixes the lesson of `nonmonotone-audit` G1 for these consumers.
* **Completeness.** It is not needed for Lemma 2. I did not re-run the [N] count. I checked its contour bookkeeping
  in `verify_zeros.py`: #zeros of Λ(χ) with |Im| < T = (D(χ)+D(χ̄))/π, which matches the windows of χ plus the
  windows of χ̄. It is consistent. Independently, F − W = 0.071645 against the density estimate 0.071665 is
  consistent with no zero missing.

## 3. S_1(u₀), precision and sensitivity (task item 3)

* `audit/mp_check.py` works at 80 digits. It parses u₀ = 18133579.0743 as an exact decimal (Arb's `arb("…")` in
  `certify.py` also encloses the exact decimal).
  * **S_T(u₀) = −4.34330566187385447981405661174.**
  * A-priori error from |γ − d| < 10⁻³⁰: ≤ 1.2·10⁻²² (that is, Σ w·u₀·10⁻³⁰).
  * **S_T + tail = −4.2716608455637996467667**, margin **0.27166084556380** below −4.

  All of these equal the Arb certificate.
* **Sensitivity.** dS_T/du(u₀) = −0.0029 (u₀ sits at a local minimum). S_T(u₀ ± δ) for δ from 10⁻⁹ to 0.05
  varies smoothly:
  * linear behaviour at 10⁻⁹;
  * quadratic behaviour with S″ ≈ 98 at 10⁻⁴ to 10⁻³;
  * −4.338 and −4.341 at ±0.01.

  There is no precision artefact.
* **Theorem A′.** With my own phases arg χ_n(a) = 2π ind(n) ind(a)/12 (labels verified as Conrey by the Krawczyk
  run), the six certificates reproduce to all digits:

  | a | bound |
  |---|---|
  | 1 | −4.27166 |
  | 3 | −4.16187 |
  | 4 | −4.17302 |
  | 9 | −4.25216 |
  | 10 | −4.24720 |
  | 12 | −4.18781 |

  The window contribution 2w cos(γu − arg χ(a)) follows from χ̄(a)e^{iγu} + χ(a)e^{−iγu}, which I re-derived.
* q = 17 (u₀ = 11761.4153), q = 19 (165.7829) and q = 11 (16133306.6866): bounds −7.896, −8.895 and −2.681,
  identical to RESULT.

## 4. From one point to infinitely many x and k (task item 4)

* Under GRH_13, S is an absolutely convergent trigonometric series (Σ weights = F_13 < ∞).
* Lemma 3 (Dirichlet pigeonhole on finitely many ordinates, tail < ε/4) is correct as written: τ = n t₀ > U and
  ‖nα_j‖ < 1/Q.
* The O-term of Thm 1.18 is a genuine O(1/(√x log² x)) for fixed (q, a), uniform in x ≥ X₀. I checked its
  ingredients:
  * Prop. 3.1: O(1/x).
  * Prop. 3.5: θ(√x; q, b) by the PNT, plus x^{1/3}.
  * Lemma 3.4, with Σ|ρ|/|1−ρ|² convergent.
  * Lemma 6.4: O(1/x).

  So "−S(u) − 4 ≥ δ − ε" at u = u₀ + τ gives log f > 0 once u > C₀/(δ − ε). Since f is constant on
  [p̄_k, p̄_{k+1}), the failure is a failure at k = π(x; 13, 1) → ∞.
* The lim sup statements (≥ 4.2716 for (7), ≥ 0.2716 for 12√x log x · log f) follow.
* Remark 2 is honest: there is no claim at x₀ itself.

## 5. Corollary B (task item 5)

* Thm 1.16(a) requires only: SH_{q,a} false, 1 − Θ < b < 1/2, and no singularity of 𝓛 on (1−b, 1). Its proof
  (Thm 4.3 + Prop. 3.6) is the standard Landau/Ingham argument. I checked it line by line:
  * G(s) continues to W_{δ,b};
  * Landau gives σ₀ ≤ 1−b;
  * the log-branch singularity of H₁ at ρ, with Re ρ > 1−b, gives the contradiction;
  * the −1/(x^{1/2} log x) term is o(x^{−b}) because b < 1/2.

  It is sound. A cosmetic slip: d²H₁ and d²H₂ are interchanged in the displayed expansions; immaterial.
* For a = 1 the residues of 𝓛 are Σ_χ m_ρ(χ) > 0, so SH_{13,1} ⇔ GRH_13, and ¬GRH_13 gives Θ > 1/2 (functional
  equation).
* The hypothesis is a **real-zero condition**, and it is met:
  * certified in RESULT by the [N] count;
  * by Platt via the paper's Cor. 1.13;
  * **independently by me** (`audit/real_zeros.py`, Arb balls on [1/2, 0.95] plus a Cauchy/maximum-principle
    bound on [0.95, 1]): L(σ, χ) ≠ 0 for σ ∈ [1/2, 1] for all 11 χ mod 13.

  ζ has no zeros in (0, 1).
* So Corollary B follows from GRH_13 ∨ ¬GRH_13. "No hypothesis" is right in the usual sense: it depends on a
  peer-reviewed theorem (Math. Comp.) whose proof I checked, not on any conjecture. Recommended wording:
  "unconditionally (using Akbary–Francis Thm 1.16(a))".

## 6. The explicit q = 17, k = 5 failure (task item 6)

* **C(17,1).** I used my own formula: φ log C(q,1) = −γ − log(1−1/q) + Σ_{χ≠χ₀}[−log|L(1,χ)| +
  Re Σ_p Σ_{k≥2}(χ(p)^k − χ(p))/(k p^k)]. Only real parts enter for a = 1, so no branch issue arises. Truncating
  at p ≤ 10⁷ gives a rigorous tail error < 2·10⁻⁷ (`audit/mertens_C.py`):
  * log C(17,1) = 0.0854342802 ± 2·10⁻⁷, so 1/C(17,1) = 0.9181134793.
  * Cross-check against the direct product over primes ≤ 10⁹: E(10⁹) = 0.08543470, off by 4·10⁻⁷, which is the
    expected size of log f there.
  * The same method gives log C(13,1) = 0.1003542810, against RESULT's 0.1003542812. It also gives
    C(5,1) = 1.225238, which confirms RESULT's erratum: the paper's figure caption quotes C(5,1) as C^{−1}.
* **Primes.** 103, 137, 239, 307, 409 are the first five primes ≡ 1 (mod 17); checked by trial division.
* **(3) at N̄_5** = 423464858827 (`audit/explicit17.py`, 40 digits):
  * LHS = N̄_5/(φ(N̄_5)(log(16 log N̄_5))^{1/16}) = **0.9179425463**;
  * RHS = 1/C(17,1) = **0.9181134793**;
  * log(LHS/RHS) = **−1.86·10⁻⁴**, about 900 times the error bar.

  **(3) fails at k = 5**; also at k = 6 (−6.3·10⁻⁶) and k = 7. The count up to 10⁹ is 175 failing k, all below
  10⁶, as in RESULT.
* I also confirmed q = 23, k = 5 (−2.7·10⁻⁴) and q = 43, k = 4 (−1.1·10⁻⁴).
* **What the paper says.** Thm 1.22 (all k ≥ 1) is claimed only for q ≤ 10 and q = 12, 14. For q > 14 the paper
  claims only Cor. 1.15 ((3) holds for infinitely many k) and Cor. 1.17 (non-squares), and it notes F_p > 4 for
  7 < p ≤ 149. It nowhere asserts (3) for all k when q = 17, and it does not mention any failure for a = 1. (The only
  failures it shows are for non-square a mod 7 in its plots.)
* So the q = 17 failure is **not** a counterexample to anything the paper claims. It shows, by a five-prime
  computation, that the literal q = 17 analogue of Thm 1.22 ("GRH_17 ⇔ (3) for all k ≥ 1") has a false right-hand
  side. Correct wording: "the q = 17 analogue of Theorem 1.22 holds only if GRH_17 is false."
* The direction of the inequality and the power of the logarithm are correct.

## 7. F_q table spot checks (task item 7)

`audit/fq_spot.py` uses my own character tables for composite q, a brute-force conductor and primitive lift, and
mpmath:

| q | F_q | 2R_{q,1} | verdict |
|---|---|---|---|
| 8 | 0.7532641 | 8 | < (paper 0.75326) |
| 12 | 0.6451559 | 8 | < (paper 0.64516) |
| 15 | 1.9967427 | 8 | < |
| 36 | 6.2576883 | 8 | < |
| 60 | 8.4580905 | 16 | < |
| 120 | 32.473283 | 32 | > |
| 11 | 4.2609839 | 4 | > |

All agree with `fq_table_1000.txt` to about 10⁻⁶ (my numerical derivative limits the accuracy) and give the same
verdicts.

## 8. Novelty (task item 8)

* **The paper.** It computed F_q (Table 1, F_11 = 4.26098, F_13 = 6.45484), noted F_p > 4 for 7 < p ≤ 149, asked
  for the finite set {q : F_q < 2R_{q,1}} and posed Question 1.24. It does **not** decide the sign of lim sup − 2R
  for any q with F_q > 2R, and it does not state any failure of (3) for a = 1.
* **Later work.**
  * Crossref, OpenAlex and Semantic Scholar list one citing work, Aoudjit–Berkane–Dusart 2021 (NNTDM, on Robin's
    criterion; unrelated).
  * Francis's 2017 Lethbridge thesis (read; text-searched) has nothing on q = 11, 13 or 17.
  * Web searches found nothing on q = 11, 13, 17 or on Question 1.24.
* **Residual risk:** the full journal text (G1).
* **Novelty statement.** As far as can be determined, the following are new:
  * the q = 13 results, i.e. the first modulus where the analogue of Thm 1.22 provably fails (in the "infinitely
    often" sense) unconditionally;
  * the explicit finite failures for a = 1;
  * the 13 extra moduli.

  Question 1.24 is **not** answered. Theorem A gives only lim sup ≥ 4.2716 ≈ 0.662 F_13.

## 9. Findings

### RIGOR BREAK

None.

### GAP

* **G1 (novelty / citation).** RESULT says "Journal version: not checked" and "I have not searched the later
  literature". Both must be fixed before any claim of novelty is made public:
  * Cite Math. Comp. 89 (2020) 993–1026, DOI 10.1090/mcom/3463, and say that statement numbers refer to arXiv v1.
  * Obtain the journal PDF (institutional access) and check that it adds nothing on q = 11, 13, 17 (the abstract
    adds nothing).
  * Record the search done here (§8).

  Until then, word the claims as "not found in the arXiv version, the 2017 thesis, or the citing literature we
  could find".
* **G2 (statement of Theorem A).** The verdict says GRH_13 is "covered by the Lean-checked, unrefereed proof
  announced on 2026-11-20". This is unverifiable here and does not belong in a theorem statement; a reader could
  take Theorem A as unconditional. §8 of RESULT says it is not used as a fact.
  *Fix:* delete the parenthesis from §0. If needed, keep a separate caveat outside the statement.

### MINOR

1. **Strengthen §7.2.** "Under GRH_17 (resp. GRH_19), infinitely many failures" can be made unconditional by
   copying Corollary B. If GRH_q fails, Thm 1.16(a) applies, because Platt (the paper's Cor. 1.13, q ≤ 400 000)
   excludes real zeros. My own ball check (`audit/real_zeros.py`) also excludes real zeros in [1/2, 1] for all
   15 characters mod 17 (run with h = 0.03, r = 0.25, N = 3000; `audit/logs/real_zeros_q17.log`). For q = 19 I
   did not run the check; Platt covers it.
2. **§6, q = 11.** "If the ordinates are linearly independent, sup S_1 = F_11 > 4" should read
   **inf S_1 = −F_11 < −4** (Kronecker; sup S_1 = F_11 holds trivially, at u = 0). The certified
   "lim sup ≥ 2.6813 for q = 11" needs GRH_11; state the hypothesis.
3. **Paper, Cor. 1.19.** The "only if" direction is not justified at the boundary lim sup = 2R; there the
   lower-order terms decide. This does not affect Theorem A, which has strict inequality with margin 0.27.
   RESULT uses only the valid direction.
4. **Corollary B consequence.** The sentence "'GRH_13 ⇔ (3) for all k' is equivalent to 'GRH_13 is false'" is
   correct but convoluted. Clearer: "unconditionally, (3) fails for infinitely many k, so no criterion of the form
   'GRH_13 ⇔ (3) for all k ≥ k₀' holds unless GRH_13 is false."
5. **Lemma 2 proof text.** It says "ball arithmetic (256 bits) … ± 10⁻²²". Record that the ordinates are known to
   about 10⁻⁴² (§2 here). The 10⁻³⁰ window is a deliberately loose certificate.
6. **§3.1 claims.** "Exactly one simple zero per window" and "no zero off the line below 400" rest on the [N] count.
   I re-certified uniqueness and simplicity in each window by Krawczyk. I did not re-run the count, which is not
   needed for any stated theorem.
7. **Floating-point sweep.** The §7.1 sweep is double/long double. RESULT says so, and the five Arb-checked cases
   carry the claim. Of those I re-checked (17,5), (23,5) and (43,4) independently.

## 10. What is established

* **Theorem A (RH + GRH for the 11 non-principal L(s,χ) mod 13).** Inequality (3) of Akbary–Francis with q = 13,
  a = 1 fails for infinitely many k. Moreover lim sup_{x→∞} Σ_χ Σ_{ρ∈Z(χ′)} x^{iImρ}/(ρ(ρ−1)) ≥ 4.27166 > 4 = 2R_{13,1},
  so condition (8) of Cor. 1.19 fails. **Established.**
* **Theorem A′ (same hypothesis).** The same holds for every square class a ∈ {1, 3, 4, 9, 10, 12} mod 13.
  **Established.**
* **Corollary B (no unproved hypothesis; uses Akbary–Francis Thm 1.16(a), Math. Comp. 2020).** Inequality (3) with
  q = 13, a = 1 fails for infinitely many k. Hence no Nicolas-type criterion "GRH_13 ⇔ (3) for all k ≥ k₀" can hold
  unless GRH_13 is false. **Established.** The real-zero hypothesis is certified independently.
* **Explicit, unconditional.** (3) with q = 17, a = 1 fails at k = 5, N̄_5 = 103·137·239·307·409:
  0.9179425 < 1/C(17,1) = 0.9181135, margin 1.9·10⁻⁴ in log against an error bar of 2·10⁻⁷. It also fails at
  q = 23, k = 5 and q = 43, k = 4. **Established.**
  * This is not a counterexample to the paper (which claims nothing for q = 17). It shows that the literal q = 17
    analogue of Thm 1.22 has a false right-hand side.
  * With MINOR 1, q = 17 and 19 also fail for infinitely many k unconditionally.
* **Table q ≤ 1000.** It is consistent with every spot check.
* **q = 11.** Open. Under GRH_11 one has only lim sup ≥ 2.6813 < 4.
* **Question 1.24.** Not answered.
