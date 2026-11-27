# SOURCES-CHECK — verbatim source statements for the tau-race write-up

Date of check: 2026-11-26. Everything below marked **[READ]** was seen by me in the
source text itself (PDF text layer, or the rendered page image for the scanned Zagier
paper and the Balkanova–Frolenkov page). Everything marked **[NOT OBTAINED]** was not
seen and nothing about it is reconstructed from memory. **[READ — SECONDARY]** means I
read a *later paper's quotation* of the statement, not the original.

Conventions for the quotations: the wording, equation numbers, constants and inequality
signs are exactly as printed. Displayed formulas are transliterated to one-line
ASCII/Unicode (fractions written with "/", products/sums written with ∏/Σ, a bar written
as ᾱ); nothing else is changed. Page numbers are printed page numbers unless "PDF p." is
written.

Downloaded PDFs are in (outside the repository)
`/tmp/claude-0/-home-user-weighted-prime-races/559a4317-d376-57ac-9c95-b35c412a8522/scratchpad/sources/`.

---

## A. Functional equation of the symmetric square L-function, level 1, weight k

### A.1 Zagier 1977 — **[READ: full statement seen, page images of pp. 106, 108, 109, 111, 155]**

D. Zagier, "Modular forms whose Fourier coefficients involve zeta-functions of quadratic
fields", in *Modular Functions of One Variable VI*, Lecture Notes in Math. 627, Springer
(1977), pp. 105–169.
URL fetched: `https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1007/BFb0065299/ModFormsZetaValues.pdf`
(linked from `https://people.mpim-bonn.mpg.de/zagier/`; scanned, 65 PDF pages, PDF p. n = printed p. 104+n).

Setting, p. 108 (§1 "Identities for the Rankin zeta function"):

> "Throughout §§ 1 – 4 we restrict ourselves to modular forms for the full modular group
> Γ = SL_2(Z)/{±1}; the results could be generalized to arbitrary congruence subgroups,
> but this would involve considerable technical complication and no essentially new ideas.
> We denote by k an even integer >2, by S_k the space of cusp forms of weight k on Γ,
> equipped with the Petersson scalar product [...]"

Definition, p. 108, eq. (2):

> "For each normalized Hecke eigenform f(z) = Σ a(n) q^n we set
>
> (2)  D_f(s) = ∏_p (1 − α_p² p^{−s})^{−1} (1 − α_p ᾱ_p p^{−s})^{−1} (1 − ᾱ_p² p^{−s})^{−1}   (Re(s) > k),
>
> where the product is over all primes and α_p, ᾱ_p are defined by
>
>   α_p + ᾱ_p = a(p),  α_p ᾱ_p = p^{k−1}
>
> (by Deligne's theorem, previously the Ramanujan-Petersson conjecture, the numbers α_p
> and ᾱ_p are complex conjugates). The function D_f(s) is related to the Rankin zeta
> function by"

Functional equation, p. 109, eqs. (3)–(5):

> "(3)  D_f(s) = [ζ(2s−2k+2) / ζ(s−k+1)] Σ_{n=1}^∞ a(n)² / n^s
>
> and hence, by the results of Rankin [17], has a meromorphic continuation to the entire
> complex plane, satisfies the functional equation
>
> (4)  D_f^*(s) = 2^{−s} π^{−3s/2} Γ(s) Γ((s−k+2)/2) D_f(s) = D_f^*(2k−1−s),
>
> and is related to the norm of f in the Petersson metric by
>
> (5)  (f,f) = [(k−1)! / (2^{2k−1} π^{k+1})] D_f(k)."

Answers to the specific questions:
- **The factor 2^{−s} IS printed** in (4). The Gamma arguments are **Γ(s)** and
  **Γ((s−k+2)/2)**; the power of π is π^{−3s/2}; the reflection is s ↔ 2k−1−s; no sign
  (i.e. sign +1) and no conductor factor.
- **Attribution of (4):** to "the results of Rankin [17]" ([17] = R. A. Rankin,
  "Contributions to the theory of Ramanujan's function τ(n) and similar arithmetical
  functions, I", Proc. Camb. Phil. Soc. 35 (1939) 351–372), via identity (3). Note that at
  this point only "meromorphic continuation" is asserted.
- **Holomorphy ("entire"):** stated separately, p. 111, attributed to Shimura with
  Zagier's own second proof:

> "From equation (11) we deduce that D_f(s+k−1) is also holomorphic in this region. On the
> other hand, the Euler product defining D_f(s) is absolutely convergent for Re(s) > k, so
> D_f(s) is certainly holomorphic in this half-plane and, by the functional equation, also
> in the half-plane Re(s) < k−1. Theorem 1 therefore implies the following result, which
> was proved by Shimura [21] in 1975 by a different method.
>
> Corollary 1 (Shimura): The function D_f(s) defined by (2) has a holomorphic continuation
> to the whole complex plane."

  ([21] = G. Shimura, "On the holomorphy of certain Dirichlet series", Proc. London Math.
  Soc. 31 (1975) 79–98.) Also p. 106 (Introduction): "As another corollary of the main
  identity we obtain a new proof of a recent result of Shimura [21] on the holomorphy of
  the Rankin zeta function." And p. 155 (in the later part on forms on Γ_0(D) with
  Nebentypus, i.e. NOT the level-1 case): "First of all, the functions D_f(s) is entire
  (proved by Shimura, [21], Theorem 1) and satisfies a functional equation (proved by Asai
  [1], Theorem 3)."
- Caveat: what Zagier states to be holomorphic on all of C is **D_f(s)** (Corollary 1), not
  literally the completed D_f^*(s). The word "entire" for the *completed* function is in
  source A.2 below. (It also follows from (4) + Corollary 1 — see Bookkeeping.)
- Zagier's standing hypothesis is "k an even integer >2"; k = 12, f = Δ is covered
  (dim S_12 = 1, Δ is the normalized eigenform).

### A.2 Balkanova–Frolenkov 2018 (published version) — **[READ: full statement seen, page image of p. 36]**

O. Balkanova, D. Frolenkov, "The mean value of symmetric square L-functions", Algebra &
Number Theory 12:1 (2018), 35–59 (first page 35; doi 10.2140/ant.2018.12.35).
URL fetched: `https://msp.org/ant/2018/12-1/ant-v12-n1-p02-s.pdf`. Page 36, §2:

> "Let S_{2k}(1) denote the space of holomorphic cusp forms of weight 2k ≥ 2 with respect
> to the full modular group. Denote by H_{2k} the normalized Hecke basis for S_{2k}(1).
> Every f ∈ H_{2k} has a Fourier expansion of the form
>
>   f(z) = Σ_{n≥1} λ_f(n) n^{k−1/2} exp(2πinz),   (2-1)
>   λ_f(1) = 1.   (2-2)
>
> For ℜs > 1 the associated symmetric square L-function is given by
>
>   L(sym² f, s) = ζ(2s) Σ_{n=1}^∞ λ_f(n²) / n^s.   (2-3)
>
> Let Γ(s) be the Gamma function and define
>
>   L_∞(s) := π^{−3s/2} Γ((s+1)/2) Γ((s−1)/2 + k) Γ(s/2 + k).   (2-4)
>
> Shimura [1975] showed that the completed L-function
>
>   Λ(sym² f, s) := L_∞(s) L(sym² f, s)
>
> is entire and satisfies the functional equation
>
>   Λ(sym² f, s) = Λ(sym² f, 1 − s).   (2-5)"

**Warning: the weight here is 2k**, so for Δ one must put k = 6 in (2-4):
L_∞(s) = π^{−3s/2} Γ((s+1)/2) Γ((s+11)/2) Γ((s+12)/2).
The arXiv preprint of the same paper (arXiv:1610.06331v2, "On the mean value of symmetric
square L-functions", p. 2–3, eq. (2.4)) **[READ]** prints instead
"L_∞(s) := π^{−3s/2} Γ((s+1)/2) Γ((s+k−1)/2) Γ((s+k)/2)" with the same weight-2k
set-up — i.e. the preprint's (2.4) is inconsistent with its own weight convention and
was changed in the published version. Cite the published version.

### A.3 Liu 2016 (preprint) — **[READ: full statement seen, PDF text]**

S. Liu, "The first moment of central values of symmetric square L-functions in the weight
aspect", arXiv:1610.07652v2 (26 Oct 2016). URL fetched: `https://arxiv.org/pdf/1610.07652`.
Here the weight is k ("holomorphic cusp forms of weight k for Γ_0(1) = SL(2, Z)"). p. 1:

> "L(s, sym² f) = ζ(2s) Σ_{n≥1} λ_f(n²) / n^s   (Re(s) > 1),
> which is a GL(3) L-function by Gelbart−Jacquet [GoJa1978]. By Shimura [Sh1975], Zagier
> [Za1977], and more generally [GoJa1978], L(s, sym² f) has analytic continuation to the
> whole s-plane."

p. 5, §2.1:

> "For the symmetric square L-function L(s, sym² f) let
>   L_∞(s) = π^{−(3/2)s} Γ((s+1)/2) Γ((s+k−1)/2) Γ((s+k)/2).
> Then the complete L-function Λ(s, sym² f) = L_∞(s) L(s, sym² f) satisfies the functional
> equation
> (5)  Λ(s, sym² f) = Λ(1 − s, sym² f)."

### A.4 Other candidates

- **Shimura 1975**, Proc. London Math. Soc. (3) 31, 79–98 (Theorem 1, R(s, f, χ)):
  **[NOT OBTAINED]** — `https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/plms/s3-31.1.79`
  returned HTTP 403. Only the attributions quoted above (Zagier p. 111 and p. 155 "[21],
  Theorem 1"; Balkanova–Frolenkov "Shimura [1975] showed ...") were read.
- **Iwaniec–Kowalski**, *Analytic Number Theory*, §5.12, Theorem 5.6, Proposition 5.7:
  **[NOT OBTAINED]** — no legitimate open text found. Only secondary citations read:
  - **[READ — SECONDARY]** arXiv:2412.01802 ("A new bound on the relative error in the
    Chebotarev density theorem"), §4: "From [8, Proposition 5.7(3)], we have
    Re B(χ) = −Σ_ρ Re(1/ρ)", and "Lemma 4.2. If t ∈ R and χ ∈ Irr(G), then
    #{β + iγ : β ≥ 0, |γ − t| ≤ 1, L(β + iγ, χ) = 0} ≪ log(D_K^{χ(1)} N f_χ (|t| + 3)^{χ(1) n_K}).
    Proof. The proof is similar to that of [8, Proposition 5.7(1)]." ([8] = Iwaniec–Kowalski.)
  - **[READ — SECONDARY]** arXiv:0803.2821 (Suzuki–Ricotta–Fesenko, "Mean-periodicity and
    zeta functions"), Appendix A: their Proposition A.1 "is a slight extension of [20,
    Proposition 5.7]", with the remark "the proof of [20, Proposition 5.7] only requires
    that (s(s − 1))^r L(f, s) is an entire function of order one and L′/L(s) is bounded
    with respect to the conductor in every vertical strip contained in some right-half
    plane", and a displayed conclusion of the shape "L′(s)/L(s) = Σ_{|s−ρ|<1} 1/(s−ρ) + O(log |t|)".
  These confirm the *content* (items (1) and (3) of Prop. 5.7) but the exact IK wording,
  the O(log q(f,s)) form, equation numbers and page numbers were not seen. Do not quote IK
  verbatim on the basis of this file.
- **Cohen–Strömberg** (GSM 179): **[NOT OBTAINED]** (not attempted beyond search; no open text).
- **LMFDB** page for Sym²Δ: **[NOT OBTAINED]** (`www.lmfdb.org` returned a reCAPTCHA
  "Checking your browser" page for every URL tried).
- Kohnen–Zagier 1981, Dummigan 2001, Watkins: not fetched (A.1–A.3 sufficed).

---

## B. Explicit bound |θ(x) − x| < 0.2 x / log² x for x ≥ 3 594 641

### B.1 Dusart 2010 preprint — **[READ: full statement seen, PDF text]**

P. Dusart, "Estimates of Some Functions Over Primes without R.H.", arXiv:1002.0442v1
(2 Feb 2010; only version; 20 pp.). URL fetched: `https://arxiv.org/pdf/1002.0442`.

Introduction, p. 2:

> "We will prove the following results:
>   ϑ(x) − x < (1/36 260) x   for x > 0,
>   |ϑ(x) − x| ⩽ 0.2 x / ln² x   for x ⩾ 3 594 641."

§5 "On the differences between ϑ and identity function", p. 4:

> "Theorem 5.2. We have
>   |ϑ(x) − x| < η_k x / ln^k x   for x ⩾ x_k
> with
>   k   | 0 | 1      | 1           | 2     | 2         | 2           | 2
>   η_k | 1 | 1.2323 | 0.001       | 3.965 | 0.2       | 0.05        | 0.01
>   x_k | 1 | 2      | 908 994 923 | 2     | 3 594 641 | 122 568 683 | 7 713 133 853
> and
>   k   | 3     | 3      | 3          | 3           | 4
>   η_k | 20.83 | 10     | 1          | 0.78        | 1300
>   x_k | 2     | 32 321 | 89 967 803 | 158 822 621 | 2"

Notes: (i) the theorem has strict "<", the introduction has "⩽"; (ii) the proof (p. 4–5)
starts "We use the estimates of |ψ(x) − x| with Proposition 3.2"; the paper's inputs
include the unpublished Wedeniwski/Gourdon verification of RH for the first 10^13 zeros
(p. 2: "Recently, Wedeniwski [24] then Gourdon [9] manage to compute zeros in a parallel
way and prove that the Riemann Hypothesis is true at least for first 10^13 nontrivial
zeros."); (iii) this is an unrefereed preprint.

### B.2 Broadbent–Kadiri–Lumley–Ng–Wilk — **[READ: full statement seen, PDF text]**

S. Broadbent, H. Kadiri, A. Lumley, N. Ng, K. Wilk, "Sharper bounds for the Chebyshev
function θ(x)", arXiv:2002.11068v2 (27 Jan 2021; latest arXiv version; 33 pp.).
URL fetched: `https://arxiv.org/pdf/2002.11068`. Published as Math. Comp. 90 (2021),
no. 331, 2281–2315 (journal reference as printed in the bibliography of arXiv:2204.02588v3;
the journal PDF itself was not read — theorem numbering in the journal may differ).

§3.8 "Computational examples", arXiv p. 16:

> "In our next examples we provide valid proofs of Theorem 4.2 of [15]. The proof of the
> main theorem in [15] is incorrect as it applies an incorrect result of Ramaré [39]. The
> stated leading constant in Theorem 1.1 of [39] is off by a factor of over 100. The
> results given here are optimal, verifiable through computation.
>
> Corollary 11.2. [15, Theorem 4.2] We have
> (3.39)  |ϑ(x) − x| < M_k x / log^k x   for x ≥ X_0
> with
>   k   | 0 | 1      | 1           | 2     | 2         | 2           | 2
>   M_k | 1 | 1.2323 | 0.001       | 3.965 | 0.2       | 0.05        | 0.01
>   X_0 | 1 | 2      | 908 994 923 | 2     | 3 594 641 | 112 568 683 | 7 713 113 853
> and
>   k   | 3     | 3      | 3          | 3           | 3           | 4
>   M_k | 20.83 | 10     | 1          | 0.78        | 0.5         | 151.3
>   X_0 | 2     | 32 321 | 89 967 803 | 158 822 621 | 767 135 587 | 2"

Proof bullet for the row we need (p. 16):

> "For m_2, M_2 = 0.2, numerical confirmation gives validity for x ∈ [3 594 641, 10^8).
> Table 13 then gives validity for x ∈ [10^8, 7 · 10^11), which Table 12 extends to 10^19.
> From here, combining Table 10 and Table 9 gives validity for all x ≥ 10^19."

([15] = P. Dusart, Ramanujan J. 45 (2018) 227–251; [39] = O. Ramaré, "An explicit density
estimate for Dirichlet L-series", Math. Comp. 85 (2016), 325–356.)

Table 3 "Case k = 2", arXiv p. 5 (historical table), row:
"Dusart (1999) [12] | m_2 = 0.2, X_0 = 3 594 641 | M_2 = 0.2, X_1 = 3 594 641"
([12] = P. Dusart, "Inégalités explicites pour ψ(x), θ(x), π(x) et les nombres premiers",
C. R. Math. Acad. Sci. Soc. R. Can. 21 (1999), no. 2, 53–59) — i.e. BKLNW attribute the
0.2 / 3 594 641 bound originally to Dusart 1999 (that paper was not read).

Table 15 (arXiv p. 32–33), "For k ∈ {1, 2, 3, 4, 5}, −m_k x/(log x)^k < θ(x) − x for
x > X_0 and θ(x) − x < M_k x/(log x)^k for x > X_1": its rows are indexed by
log X_0 = log X_1 ∈ {0, log(10^5), log(5·10^5), log(10^6), log(5·10^6), log(10^7), ...};
there is **no row at 3 594 641**. Rows read in the k = 2 column (m_2 = M_2):
X_0 = 1 (log X_0 = 0): 3.9649 · 10^0;  10^6: 3.2309 · 10^{−1};  5·10^6: 1.5952 · 10^{−1};
10^7: 1.3675 · 10^{−1};  10^8: 5.0612 · 10^{−2};  10^9: 2.0087 · 10^{−2};  10^10: 8.7526 · 10^{−3}.
So the citable BKLNW statement for "0.2, x ≥ 3 594 641" is **Corollary 11.2**, not Table 15
(Table 15 gives e.g. 0.15952 for x > 5·10^6 and 0.32309 for x > 10^6).

### B.3 Dusart 2018 (Ramanujan J.) — **[NOT OBTAINED]** (original); **[READ — SECONDARY]** quotations

P. Dusart, "Explicit estimates of some functions over primes", Ramanujan J. 45 (2018),
227–251, doi 10.1007/s11139-016-9839-4. `link.springer.com` returned a JavaScript "Client
Challenge"/login redirect; no open copy found (HAL search empty). The published Theorem
4.2 was therefore NOT seen. Two later papers quote it:

- BKLNW Corollary 11.2 (quoted in B.2) is headed "[15, Theorem 4.2]".
- C. Axler, "New estimates for some functions defined over primes", arXiv:1703.08032v2,
  p. 3 **[READ]**: "In the direction of (2.5), Dusart found in [11, Theorem 5.2] and [13,
  Theorem 4.2] the following explicit estimates for the distance between x and ϑ(x).
  Lemma 2.3 (Dusart, [11], [13]). We have |ϑ(x) − x| < η_k x / log^k x for every
  x ≥ x_0(k) with
    k      | 1           | 2         | 2             | 3          | 3           | 4
    η_k    | 0.001       | 0.2       | 0.01          | 1          | 0.5         | 151.3
    x_0(k) | 908 994 923 | 3 594 641 | 7 713 133 853 | 89 967 803 | 767 135 587 | 2"
  ([11] = arXiv:1002.0442, [13] = the Ramanujan J. paper, then "to appear").

All three texts read (Dusart 2010 Thm 5.2; Axler's Lemma 2.3; BKLNW Cor. 11.2) agree on
the entry **k = 2, η_2 = 0.2, x ≥ 3 594 641, strict inequality**.

### B.4 Discrepancies, errata and caveats come across

1. **Dusart 2018's proof is flagged as incorrect by BKLNW** (quote in B.2, p. 16), and
   again BKLNW footnote 3, arXiv p. 20: "Unfortunately the main theorem of [39] is
   incorrect and thus bounds claimed in [15] are likely affected, in particular Theorems
   3.2, 3.3 and Table 1 for bounds for ψ(x), and consequently Theorem 4.2 and Table 2 for
   θ(x). In addition this unfortunately affects the main theorem [1, Theorem 1.1]."
   ([1] = Axler, Integers 18 (2018) A52.) Consequence: **cite BKLNW Corollary 11.2 (which
   re-proves the table) rather than Dusart 2018 Theorem 4.2 alone**; Dusart's 2010 preprint
   Theorem 5.2 has a different proof (via its Proposition 3.2) and is not the object of
   that remark, but is unrefereed.
2. **Digit discrepancies between the printed tables** (irrelevant to the 0.2 row):
   η_2 = 0.05: Dusart 2010 prints x_k = 122 568 683, BKLNW prints 112 568 683.
   η_2 = 0.01: Dusart 2010 and Axler print 7 713 133 853, BKLNW prints 7 713 113 853.
   η_3 = 0.78: BKLNW's table prints 158 822 621 (as Dusart 2010) but BKLNW's own proof
   bullet says "x ∈ [158 882 621, 10^10)". η_4: Dusart 2010 has 1300, Axler/BKLNW 151.3.
   Which digits are in the published Dusart 2018 was not determined.
3. **Fiori–Kadiri–Swidinsky**, "Sharper bounds for the Chebyshev function ψ(x)",
   arXiv:2204.02588v3 **[READ]**, p. 3, item (2): "It has recently been discovered that the
   Hiary [14] subconvexity bound for |ζ(1/2 + it)| ≤ a_1 t^{1/6} log t, with a_1 = 0.63
   relies on an erroneous explicit version of van der Corput second derivative test due to
   Cheng-Graham [7]. See [24, Footnote 3] for details. After accounting for this correction,
   the constant, a_1, changes to 0.77. [...] Because the value 0.63 was used in [18,
   Theorem 1.1], the error would render many results concerning the prime number theorem
   unreliable (for instance [36, 28, 30, 5])." Here [18] = Kadiri–Lumley–Ng zero-density
   (J. Math. Anal. Appl. 2018) and [30] = Platt–Trudgian, "The error term in the prime
   number theorem", Math. Comp. 90 (2021) 871–881. BKLNW state (arXiv p. 3) that they use
   "[3, Theorem 2], [4, Theorem 1] and [37, Theorem 1]" with [37] = that Platt–Trudgian
   paper, and that "Büthe's technique leads to better bounds when x < e^{2300} while Platt
   and Trudgian's works better for larger values of x". So BKLNW's ε(b) for very large b
   (and hence the "x ≥ 10^19" tail of Corollary 11.2) inherits a dependence on a result FKS
   call unreliable; FKS add that with 0.77 they "still recover, and actually improve on
   these previous results for E_ψ(x)". *My assessment (not a quotation):* the 0.2/log² x
   bound needs only ε ≲ 0.2/log² x in that range, far weaker than anything at stake, so the
   k = 2, 0.2 row is not plausibly endangered; but no published erratum to BKLNW itself was
   found, and none of the papers read says so explicitly.
4. Johnston–Yang, arXiv:2204.01980v2 **[READ]**, footnote 1: "As noted in [Bro+21, Theorem
   16] there is a small error in [Büt16, Theorem 1]. However, in the case when α = 0 there
   is no difference to the result." (An error in Büthe 2016 noted *by* BKLNW; no error *in*
   BKLNW is asserted there.) In Fiori–Kadiri–Swidinsky, arXiv:2206.12557v2, I found no
   remark correcting BKLNW (BKLNW is used as [2], e.g. "[2, Corollary 5.1]", "[2, Table 13
   and 14]").
5. Web searches for an erratum/corrigendum to BKLNW returned nothing.

### B.5 Rosser–Schoenfeld 1962 — **[NOT OBTAINED]**

J. B. Rosser, L. Schoenfeld, "Approximate formulas for some functions of prime numbers",
Illinois J. Math. 6 (1962), 64–94. Project Euclid
(`https://projecteuclid.org/journals/illinois-journal-of-mathematics/volume-6/issue-1/Approximate-formulas-for-some-functions-of-prime-numbers/10.1215/ijm/1255631807.pdf`
and `.../download/pdf_1/euclid.ijm/1255631807`) served only an Incapsula bot-challenge
page; `web.archive.org` is not reachable from this session (connection reset by the
proxy / fetch tool refuses the host). Statements (3.15), (3.16) and Theorem 9 were
therefore NOT seen. Only this secondary information was read (BKLNW, arXiv p. 4,
historical tables, where x(1 − m_k/(log x)^k) ≤ θ(x) for x ≥ X_0 and
θ(x) ≤ x(1 + M_k/(log x)^k) for x ≥ X_1):

- **[READ — SECONDARY]** Table 1 "Case k = 0", block labelled "Rosser & Schoenfeld (1962)
  [42]": rows (m_0, X_0) = (0.16, 101), (0.05, 1 427), (0.02, 7 481), and a row
  (M_0, X_1) = (0.01624, 1) — i.e. θ(x) ≤ 1.01624 x for x ≥ 1, consistent with "Theorem 9:
  θ(x) < 1.01624 x for x > 0" (the label is vertically centred against these four rows in
  the text layout; theorem number not shown there).
- **[READ — SECONDARY]** Table 2 "Case k = 1", block "Rosser & Schoenfeld (1962) [42]":
  rows (m_1, X_0, M_1, X_1) = (0.50, 563, 0.50, 1) and (0.47, 569, 0.47, 569) — i.e.
  θ(x) ≤ x(1 + 1/(2 log x)) for x ≥ 1 [x > 1] and θ(x) ≥ x(1 − 1/(2 log x)) for x ≥ 563.
  The form "θ(x) > x(1 − 1/log x) for x ≥ 41" suggested in the task was **not** confirmed
  by anything read.
- BKLNW p. 2 also quote "x/log x < π(x) < 1.25506 x/log x, for all x ≥ 17 [42, Corollary 1]".

Equation numbers (3.15)/(3.16)/Theorem 9 must be checked against the paper before being
cited.

---

## C. arXiv:2408.00187 — **[READ: full text, PDF]**

URL fetched: `https://arxiv.org/pdf/2408.00187` (v1, 31 Jul 2024, the only version;
"28 pages, 3 figures"; MSC 11M06, 11Y35).

- **Title:** "A method for verifying the generalized Riemann hypothesis"
- **Authors:** Ghaith Hiary, Summer Ireland, Megan Kyi
- **Abstract (verbatim):** "Riemann numerically approximated at least three zeta zeros.
  According to Edwards, Riemann even took steps to verify that the lowest zero he computed
  was indeed the first zeta zero. This approach to verification is developed, improved, and
  generalized to a large class of L-functions. Results of numerical calculations
  demonstrating the efficacy of the method are presented."

**Method in two sentences (my summary of pp. 2–3, 8–13, 18–19).** Using the Hadamard
product of the completed (entire, order 1) function ξ_L, the sum over all nontrivial zeros
w_{1,δ} = Σ_ρ 1/(ρ − δ) is evaluated exactly as
(1/2) log N − (r/2) log π + (1/2) Σ_j ψ_0((1 − δ + μ_j)/2) + L′/L(1 − δ) (Corollary 5),
which needs only one rigorous evaluation of L′/L at a point of absolute convergence
(δ ≤ 0, e.g. δ = −1; truncation bound Lemma 6; interval arithmetic). Because for δ ≤ 0
every zero (with its symmetric partners) contributes a nonnegative amount, the
contribution C(Z, δ) of the zeros already located on the critical line by sign changes of
ξ_L(1/2 + it), plus the minimal contribution of a hypothetical off-line (or missed, or
multiple) zero of height ≤ η, would exceed w_{1,δ}; this contradiction proves RH — and
completeness/simplicity of the zero list — up to height η (Theorem 7), with no use of the
argument principle or of Turing's method.

Supporting quotations:
- p. 3: "Specifically, rather than fall back on a generalized Backlund method, which would
  require using a numerically-involved application of the argument principle, we re-examine
  and develop the Riemann method in a more general setting. [...] Our generalization is
  also simple to derive and justify, requiring a single numerical evaluation of a
  logarithmic derivative of the L-function at a special point in the region of absolute
  convergence."
- p. 8 (§4, hypotheses on the class): "we use the analytic normalization of L-functions, so
  the critical line is σ = 1/2. We consider L-functions of order 1 only."; "We suppose that
  the Dirichlet coefficients a(n) are real"; Euler product of degree r, "We will further
  assume that |α_{p,j}| ≤ 1"; "(8) ξ_L(s) := γ(s)L(s), γ(s) := ϵ N^{s/2} π^{−sr/2}
  ∏_{j=1}^r Γ(s/2 + μ_j/2), extends to an entire function and satisfies the functional
  equation (9) ξ_L(s) = \overline{ξ_L(1 − \bar s)}" [the bars are lost in the PDF text
  layer; the surrounding text ("If ϵ = ±1, then the functional equation simplifies to
  ξ_L(s) = ξ_L(1 − s)") shows the general form is the conjugate one]; p. 9: "We assume
  L(1) ≠ 0".
- p. 12, Theorem 7 (v): "If F(η, δ) + C(Z, δ) > w_{1,δ}, then the list of zeros in Z is
  complete. This means that every ρ in the upper half-plane of height ≤ η is on the
  critical line, is simple, and belongs to some subinterval in the set Z."
- p. 19, §7.3 "The Ramanujan τ L-function": "In the notation of §4, we have r = 2, N = 1,
  ϵ = 1, μ_1 = 11/2, and μ_2 = 13/2. We applied Theorem 7 to verify the RH using δ = −1,
  τ = 9877.3. [...] Our set Z was obtained using [Rub] and [LMFDB]. With the aid of
  Corollary 5 and Lemma 6, applied with K = 10^5, we computed
  w_{1,δ} ∈ [0.1671717623, 0.1672414682], C(Z, δ) ∈ [0.1663983945, 0.1663983946]. Based on
  this input data, Theorem 7 succeeded in verifying the RH for L(s) for η = 84, a window
  which includes 46 zeros with nonnegative ordinates."
- Other worked examples (§7): ζ near height 10^28 + 501675.8 (Theorems 11, 13, using
  explicit bounds on S(u)), L(s, χ_d) with d = −1159523 (η = 32, 74 zeros), the elliptic
  curve 37.a1 (η = 10, 5 zeros). The paper does **not** treat a symmetric square
  L-function. (My remark, not from the paper: Sym²Δ would fit the paper's framework (8)
  with r = 3, N = 1, ϵ = 1, (μ_1, μ_2, μ_3) = (1, 11, 12), by A.2/A.3.)
- p. 18: input zeros "from [LMFDB] and [Zet]" are assumed "accurate within ±10^{−10}"
  ("Our working assumption is ..."), computations in mpmath interval arithmetic.

---

## Bookkeeping for (A)

Take Zagier's (4): D_f^*(s) = 2^{−s} π^{−3s/2} Γ(s) Γ((s−k+2)/2) D_f(s) = D_f^*(2k−1−s),
with α_p ᾱ_p = p^{k−1}. Put s = s′ + k − 1. Then α_p² p^{−s} = (α_p p^{−(k−1)/2})² p^{−s′}
etc., so D_f(s′ + k − 1) = L(s′, Sym² f) in the analytic normalisation (Euler product
absolutely convergent for Re s′ > 1), and 2k − 1 − s = (1 − s′) + k − 1, so the functional
equation is **s′ ↔ 1 − s′ with sign +1 and conductor 1**. For the gamma factor, with
Γ_R(s) = π^{−s/2} Γ(s/2) and Γ_C(s) = 2 (2π)^{−s} Γ(s):

  Γ_R(s′+1) Γ_C(s′+k−1) = 2^{−(s′+k−2)} π^{−(3s′+2k−1)/2} Γ((s′+1)/2) Γ(s′+k−1),
  2^{−s} π^{−3s/2} Γ(s) Γ((s−k+2)/2) |_{s=s′+k−1} = 2^{−(s′+k−1)} π^{−3(s′+k−1)/2} Γ(s′+k−1) Γ((s′+1)/2),

hence

  **D_f^*(s′+k−1) = (1/2) π^{−(k−2)/2} · Γ_R(s′+1) Γ_C(s′+k−1) L(s′, Sym² f).**

The quotient is a constant ((1/2)π^{−5} for k = 12) precisely because of the printed
factor 2^{−s}; without it one would get 2^{s′+k−2} π^{−(k−2)/2} Γ_R Γ_C, which is not a
constant multiple (this agrees with the reference computation in the task). Checked
numerically at s′ = 0.3, 1.7, 2.5 for k = 12 (ratio 0.0016338818215267 = (1/2)π^{−5}).

Cross-check with A.2/A.3: by the duplication formula
Γ((s+κ−1)/2) Γ((s+κ)/2) = 2^{2−s−κ} π^{1/2} Γ(s+κ−1), the Liu / Balkanova–Frolenkov factor
(weight κ; κ = 2k in B–F) is L_∞(s) = π^{−3s/2} Γ((s+1)/2) Γ((s+κ−1)/2) Γ((s+κ)/2)
= **π^{κ} · Γ_R(s+1) Γ_C(s+κ−1)** (numerically confirmed for κ = 12: ratio = π^{12}),
with Λ(s) = Λ(1 − s). So the three sources agree: for f = Δ,

  Λ(s, Sym²Δ) := Γ_R(s+1) Γ_C(s+11) L(s, Sym²Δ) = Λ(1 − s, Sym²Δ),

equivalently gamma factor π^{−3s/2} Γ((s+1)/2) Γ((s+11)/2) Γ((s+12)/2) up to the constant π^{12}.

Entirety of the completed function: stated verbatim in A.2 ("is entire", attributed to
Shimura 1975). From Zagier alone it is a one-line deduction (mine, not quoted): by
Corollary 1, D_f is entire; Γ(s)Γ((s−k+2)/2) is holomorphic for Re(s) > k − 2, so D_f^* is
holomorphic there, and by (4) also for Re(s) < k + 1; the two half-planes cover C.
Zagier's (3) and Liu's/B–F's Dirichlet series ζ(2s) Σ λ_f(n²) n^{−s} are the same function
under s = s′ + k − 1 (not separately quoted from Zagier beyond (2)–(3)).

Suggested citation chain: functional equation — Rankin 1939 via Zagier 1977, eq. (4),
p. 109; holomorphy — Shimura 1975, Theorem 1 (as cited by Zagier 1977, Corollary 1,
p. 111, who gives an independent proof); "completed function entire, s ↔ 1 − s" in the
analytic normalisation — Balkanova–Frolenkov, ANT 12 (2018), (2-4)–(2-5), p. 36 (weight 2k
there). Shimura's and Iwaniec–Kowalski's own wording remain unverified here.
