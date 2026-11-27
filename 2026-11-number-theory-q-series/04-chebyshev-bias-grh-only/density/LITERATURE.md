# Literature: what is known without LI about the existence of the density

AI-produced (Anthropic's Claude), 2026-11-26. Two web checks by a subagent (search engines,
arXiv abstract/HTML pages); each item is tagged **READ** (seen in fetched text) or
**MEMORY** (recalled, not re-read). Nothing here was checked against a journal version.

## Top findings

1. **Open.** No result was found giving, under GRH alone, the existence of the logarithmic
   density of `{x : π(x;4,3) > π(x;4,1)}`, nor the absence of an atom at `0` of the limiting
   distribution. No preprint posted since 2026-11-20 claiming it was found (coverage: search
   engines only, not a full scan of the arXiv listings).
2. **Problem number.** In Hamieh–Kadiri–Martin–Ng, *Comparative prime number theory problem
   list*, arXiv:2407.03530, the relevant item is **Problem 13** (attributed to Pomerance,
   posed by Martin; READ in the HTML version): "Prove unconditionally that the set of ties
   `π(x;q,a) = π(x;q,b)` has natural or logarithmic density 0." Its background says that
   under GRH and a very mild version of LI (Devin, Thm 2.2) the limiting distribution is
   continuous, so near-ties `|π(x;q,a) − π(x;q,b)| < f(x)`, `f = o(√x/log x)`, have log
   density 0. **Problem 12** is Koyama's `1/√p`-weighted formulation — a different question;
   nothing in this directory bears on it (RESEARCH-CONTEXT §12). **Problem 1** of the list
   asks to prove that not all ordinates are rational multiples of one another; this is open
   and is relevant to §6 of `DENSITY.md`. No problem of the list asks verbatim "how much of
   LI is needed for the density to exist" (the list was not read in full).
3. **Weakest hypothesis in print** (as far as found): Devin, *Chebyshev's bias for analytic
   L-functions*, arXiv:1706.06394 (Math. Proc. Cambridge Philos. Soc. 2020), **Theorem 2.2**:
   GRH-type hypothesis plus the existence of `λ_1,…,λ_N` among the ordinates with
   `⟨λ_1,…,λ_N⟩_ℚ ∩ ⟨Z*_S ∖ {λ_1,…,λ_N}⟩_ℚ = {0}`; conclusion: the limiting distribution is
   continuous and the density exists (READ in paraphrase; statement as recalled, see the
   verification block at the end). One self-sufficient ordinate (Definition 5.3: not in the
   ℚ-span of the others) is the case `N = 1`. Bailleul, arXiv:2007.05763, restates this
   (his Thm 2.15) and applies it to number-field races (Thm 1.7) (READ).
4. The reported Lean proof of GRH (announced 2026-11-20; READ only on a Wikipedia page with
   maintenance tags; not corroborated by us) does not touch the question: the atom problem
   is about linear relations among the zeros, which GRH does not control.

## Items

| Source | What it gives without LI | Status |
|---|---|---|
| Rubinstein–Sarnak, Experiment. Math. 3 (1994) | GRH ⇒ limiting distribution `μ` exists, is the push-forward of Haar measure on the closure torus; `δ = μ((0,∞))` needs `μ({0}) = 0`, obtained from LI (then `μ` has an analytic density) | MEMORY |
| Devin 1706.06394, Thm 2.1, 2.2, 2.5, Def. 5.3, Cor. 5.8 | existence of `μ` for analytic L-functions; continuity under the direct-sum hypothesis above (proof: product splitting of the torus, decay of `μ̂`, Wiener); positivity of both densities; Chebyshev-type lower bound | READ (HTML) |
| Martin–Ng, *Inclusive prime number races*, arXiv:1710.00088 | finitely many self-sufficient zeros ⇒ races are inclusive (every ordering has positive density) | abstract READ; Thm 1.5, Def. 1.3 MEMORY |
| Fiorilli, arXiv:1210.6946, §2, Remark 1.3 | variance with squared multiplicities; Chebyshev-type bound; existence of the density stated as open without LI | MEMORY (quoted from `GRH-ONLY.md` §8, which read it) |
| Bailleul 2007.05763 | explicit Kronecker–Weyl; Devin's criterion restated | READ |
| Fiorilli–Jouve 2001.05428 | nothing on atoms/continuity without LI | abstract READ |
| Hayani 2512.23302 (and Sheth) | GRH only: natural density 1 for the `p^{-1/2}`-weighted race (endpoint `m = −1/2`, outside our range) | abstract READ |
| Hayani 2607.23150; Koyama 2607.28931; Bailleul–Hayani–Untrau 2603.20093 | not about existence of the density under GRH alone | abstracts READ |
| Akbary–Ng–Shahabi 1306.1657 | `B²`-almost periodicity, existence of limiting distributions; continuity only under LI-type hypotheses | MEMORY |
| Heath-Brown, Acta Arith. 60 (1992); Bleher et al. | absolutely continuous limiting laws for lattice-point error terms **without** assuming independence: the frequencies `√n` have a *provable* arithmetic independence (square roots of squarefree integers), which is what replaces LI. No analogue is known for zeros. | MEMORY |
| Jessen–Wintner (1935) | pure-type and continuity theorems for infinite convolutions; require independent summands | MEMORY |

Ties / near-ties: no unconditional or GRH-only result found (Problem 13 is exactly this).

## Verification block (second check; statements READ in the arXiv PDFs via pdftotext)

- **Devin 1706.06394, Theorem 2.2** (verbatim up to symbol repair): "Suppose that there
  exists `N ≥ 1` and `λ_1,…,λ_N ∈ Z*_S` such that `⟨λ_1,…,λ_N⟩_ℚ ∩ ⟨Z*_S − {λ_1,…,λ_N}⟩_ℚ = {0}`
  where `⟨·⟩_ℚ` denotes the ℚ-span of a set of real numbers. Then the distribution `μ_S` is
  continuous (i.e. `μ_S` assigns zero mass to finite sets), and `δ(S)` exists." Proof via the
  splitting `A_T = A(λ_1,…,λ_N) × A'_T`.
- **Devin, Definition 5.3**: "(1) An ordinate `γ ∈ Z*_S` is self-sufficient if it is not in the
  ℚ-span of `Z*_S − {γ}`. (2) For `U,V > 0`, we say that an ordinate `γ ∈ Z*_S(U)` is
  `(U,V)`-self-sufficient if it is not in the ℚ-span of `Z*_S(V) − {γ}`."
- **Devin, Theorem 5.4**: "(1) Suppose that there exists `ε > 0` and `T_0 > 0` such that for
  every `T ≥ T_0` the set `Z*_S(T^{1/2−ε})` contains a `(T^{1/2−ε},T)`-self-sufficient zero
  `γ_T`. Then `δ(S)` exists. (2) Suppose `Z*_S` contains at least one self-sufficient element,
  then the distribution `μ_S` is continuous. (3) … three or more self-sufficient elements,
  then `μ_S` admits a density `φ ∈ L¹`. (4) … infinite, then `μ_S` admits a density which is
  in the Schwartz space." Remark 9: (1)–(4) "are ordered by increasing strength but are all
  weaker than LI". Devin also notes that Kaczorowski–Ramaré [KR03, Th. 3] obtain existence
  of a logarithmic density "under the Riemann Hypothesis only" for a comparable set in the
  Selberg class (**not read by us**; presumably for all but countably many levels).
- **Martin–Ng 1710.00088**: Definition 1.3 (self-sufficient: not a finite ℚ-linear
  combination of the other elements of `Γ(q)`); Theorem 1.5(a): if the union of the sets of
  self-sufficient ordinates over `χ(a) ≠ χ(b)` has at least three elements, the race between
  `a` and `b` is weakly inclusive; Theorem 1.7(a): three self-sufficient zeros of `ζ` ⇒ the
  densities of `{π(x) > li(x)}`, `{π(x) < li(x)}` exist.
- **Devin 1909.03975** ("Limiting properties of the distribution of primes in an arbitrarily
  large number of residue classes"), Corollary 2: same direct-sum hypothesis
  (`⟨γ_1,…,γ_M⟩_ℚ ∩ ⟨Z(q) ∖ {γ_1,…,γ_M}⟩_ℚ = {0}`) plus a separation condition ⇒ every race mod
  `q` is weakly inclusive and ties have density zero.
- **Bailleul 2007.05763, Theorem 1.7**: `Span_ℚ Θ = Span_ℚ T ⊕ Span_ℚ(Θ ∖ T)` for a non-empty
  finite `T` (plus a non-vanishing condition) ⇒ the density exists.
- No statement of the form "a ℚ-linear form on the span of the ordinates that is nonzero on
  only finitely many ordinates", "the ordinates above height `T` do not span all ordinates",
  or "an atom forces every ordinate to be a rational combination of ordinates above any
  height" was found in these four papers (the wider literature was not searched).

Comparison with `DENSITY.md` §3: hypothesis TS is implied by Devin's Theorem 2.2 hypothesis
and is strictly weaker as a condition on a set of frequencies (Example 3.3); it is
incomparable with Theorem 5.4(1). Hypothesis FD is implied by both (it is Theorem 5.4(1)
with up to `N_0` exceptional ordinates and "arbitrarily large `T`" instead of "all
`T ≥ T_0`").

## Third check (prior art for the time-reversal identity of `DENSITY.md` §5)

Subagent, 2026-11-26; pdftotext of Devin 1706.06394, Fiorilli 1210.6946, Martin–Ng
1710.00088, Akbary–Ng–Shahabi 1306.1657, Bailleul 2007.05763, Hayani 2512.23302 (grep for
"symmetr", "time reversal", "without LI") and one web search. Rubinstein–Sarnak 1994 not read.

- **Nearest prior result: Devin, Theorem 2.3** (READ, verbatim up to symbol repair): "Suppose
  that the set `S` satisfies the conditions of Theorem 2.1. Suppose that for every
  `(k_γ) ∈ ℤ^{(Z_S)}` one has `Σ k_γ γ = 0 ⇒ Σ k_γ ≡ 0 [mod 2]`. Then the distribution `μ_S` is
  symmetric with respect to `m_S`." Proof (§5.2, Lemma 5.6): under that parity hypothesis
  the point `(1/2,…,1/2)` lies in the closure torus, so the *translation* `θ_γ ↦ −θ_γ` is
  available. This is exact symmetry under an (unproved) arithmetic hypothesis; it does not
  use inversion of the group. Without the hypothesis the paper says nothing about symmetry.
- Fiorilli (§2, Lemma 2.7: `δ ≥ 1 − 2Var/E²` when `E/√Var` is large), Devin (Lemma 5.7,
  Cor. 5.8: plain Chebyshev), Martin–Ng, Akbary–Ng–Shahabi (symmetry of the density only
  under LI), Bailleul, Hayani: **no** statement of the identity `F(θ) + F(θ^{−1}) = 2 + R`, of
  an approximate symmetry without LI, of a bound of the type "lower density `≥ 1/2` under
  GRH alone", or of any explicit numerical lower bound for the bias of the `q = 4` race
  under GRH alone, was found.
- Caveat: a 15-minute check of six papers; Rubinstein–Sarnak, Sheth and the wider
  literature were not read. The claim of novelty for Theorem 5.1 is therefore "not found
  in the papers listed", no more.

## Fourth check (power series with integer or finitely many coefficient values; for `DENSITY.md` 9.16–9.19)

Search engines plus two subagent reads of ar5iv HTML, 2026-11-26.
- **Breuer–Simon**, *Natural boundaries and spectral theory*, Adv. Math. 226 (2011),
  arXiv:1002.0823 (READ: §§1–4 from ar5iv; §§5–8 only via a model summary).
  - Thm 1.3: analytic continuation across an arc `I` implies that every right limit of the
    bounded coefficients is reflectionless on `I`.
  - Thm 1.4: the same under `sup_r ∫_I |f(re^{iθ})| dθ < ∞`.
  - "Reflectionless" is defined only for arcs.
  - Thm 5.1, per the summary and not verified: finitely many values implies a strong natural
    boundary or eventual periodicity.
  - No statement for sets of positive measure was found.
- **Sauzin–Tiozzo**, *Generalised continuation by means of right limits*, arXiv:1301.1175
  (READ: about 40%). Renascent right limits; arcs only; nothing on sets of positive measure
  or on finitely many values.
- **Borwein–Erdélyi–Littmann**, *Polynomials with coefficients from a finite set*, Trans. AMS
  360 (2008) (READ: author's revised preprint `BEL_rev.pdf`, via pdftotext).
  - §2, Thm 2: polynomials with coefficients from a finite set cannot tend to zero uniformly on
    an arc.
  - Thm 3: Duffin–Schaeffer. A series with coefficients from a finite set that is bounded in a
    sector is rational.
  - Verbatim: "It is an open problem whether Theorem 2 remains true after replacing arcs on the
    unit circle with sets of positive measure"; this would follow from an analogue of Thm 3 for
    sets of positive measure.
  - §4, Question 4 (Konyagin's conjecture): for 0/1 polynomials with distinct exponents,
    `∫_E |Σ_{j≤n} z^{k_j}| |dz| ≥ c(E)` for every fixed `E` of positive measure. Konyagin showed
    that `E` must vary with `ε` in his counterexamples.
  - A search found no later proof of either (Erdélyi proved the subarc case of Konyagin's).
- **Duffin–Schaeffer**, Amer. J. Math. 67 (1945) (MEMORY, plus search snippets). Coefficients
  from a finite set and bounded in a sector imply rational. Their integer-coefficient example,
  bounded in a sector but with natural boundary, shows that finiteness matters.
  **Borwein–Erdélyi–Littmann**, Trans. AMS 360 (2008), give a shorter proof (abstract via
  search).
- **Pólya**, Math. Ann. 99 (1928) (MEMORY). A series `Σ a_n z^{−n}` with `a_n ∈ ℤ`, analytic in
  the complement of a compact set of transfinite diameter `< 1`, is rational. The
  transfinite diameter of a circular arc of length `ℓ` on the unit circle is `sin(ℓ/4)`
  (MEMORY). **Fatou** (1906): an integer-coefficient series that represents a rational
  function is `R/Q` with `R, Q ∈ ℤ[z]` and `Q(0) = 1` (MEMORY).
- **Not found:** any fat-set (positive-measure, empty-interior) analogue of the theorems of
  Pólya, Szegő, Duffin–Schaeffer or Breuer–Simon, and anything on pseudocontinuation of
  integer-coefficient series. Nor was a statement of Theorem 9.16 found as such. It is a short
  consequence of Pólya–Carlson (rational or natural boundary), once `g` is continued across the
  arc, and is folklore-level (referee, from memory: Carlson, Math. Z. 11 (1921); Pólya, Math.
  Ann. 99 (1928); Bieberbach, *Analytische Fortsetzung* (1955); R. M. Robinson, Trans. AMS 130
  (1968)).
- Possibly relevant, **not read**: multivariate Pólya–Carlson (Yu, arXiv:2312.00428) and
  Bell–Chen–Nguyen–Zannier (arXiv:2306.02590), for a higher-rank version of Corollary 9.16′.
- For Theorem 9.22 (noise shaping), MEMORY:
  - Daubechies–DeVore, *Approximating a bandlimited function using very coarsely quantized data:
    a family of stable sigma-delta modulators of arbitrary order*, Ann. of Math. 158 (2003);
  - Güntürk, *One-bit sigma-delta quantization with exponential accuracy*, Comm. Pure Appl. Math.
    56 (2003).

  They are used only as context. The theorem uses its own greedy multi-bit scheme with filter
  `(1−z)^k`, assigned by band of origin.
