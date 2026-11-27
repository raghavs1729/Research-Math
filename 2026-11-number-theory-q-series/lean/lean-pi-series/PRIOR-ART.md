# Prior art: are these series already formalized? (search of 2026-11-24)

Searched: Mathlib (GitHub code search and file reads), other Lean 4 repositories on GitHub, the Lean
Zulip archive (web search, weak coverage), the Isabelle AFP mirror, HOL Light, Metamath `set.mm`, Coq/Rocq
libraries, arXiv 2023–2026. Nothing below was built by us; claims about other repositories are quoted
from their READMEs.

**Result: no formal proof of R1 (`(42n+5)/2^{12n}`, `16/π`), R2 (`(6n+1)/2^{8n}`, `4/π`) or Bauer's series
(`(-1)^n (4n+1)/2^{6n}`, `2/π`) was found in any proof assistant.** So all three were done here.

## What exists
- **Mathlib**: `Mathlib/Analysis/Real/Pi/Chudnovsky.lean` only *defines* the Chudnovsky sum ("does not (yet!)
  contain a proof"); `proof_wanted chudnovskySum_eq_pi_inv`. `Mathlib/NumberTheory/ModularForms/RamanujanFormula.lean`
  is about derivatives of Eisenstein series. No Ramanujan-type `1/π` series, no WZ machinery.
  `docs/1000.yaml`: only Ramanujan–Skolem and Hardy–Ramanujan, without declarations.
- **google-deepmind/formal-conjectures**: states Gourevitch's `32/π³` series with `sorry`; none of R1, R2, Bauer, Sun.
- **NERV-WMY/pi2-lean**: README claims a sorry-free proof of `Σ (6n)!/(n!)^6 (532n²+126n+9)/10^{6n} = 375/(4π²)`
  (a `1/π²` series, different method). Not built by us.
- **Ripple** (arXiv:2607.13531): the 1103/26390 series, conditional on a Chowla–Selberg hypothesis.
- **xiangyazi24/Ramanujan_Challenge**: Chudnovsky, conditional on an unformalized CM evaluation, uses `native_decide`.
- **bjbraams/lean-DLMF**: Dixon's `₃F₂(1)` by a WZ pair (WZ telescoping infrastructure, not a `1/π` series).
- **Isabelle AFP / HOL Light / Metamath / Coq**: nothing for a Ramanujan-type `1/π` series (BBP-type formulas exist in
  HOL Light and Coq; Chyzak–Mahboubi–Sibut-Pinote–Tassi formalized Apéry's recurrences by creative telescoping in Coq).
- **Sun's Conjectures 5.25(i) and 5.33(i)** (arXiv:2307.03086v10, Feb 2026; restated as conjectured in
  arXiv:2603.29973, Mar 2026): no proof and no formalization found other than this repository.

## The WZ pairs in the literature (transcribed from the arXiv TeX sources)
- **Ekhad–Zeilberger**, *A WZ proof of Ramanujan's formula for π* (arXiv:math/9306213): proves **Bauer's** series
  `2/π = Σ (-1)^k (4k+1) (1/2)_k³/k!³` (not the `42n+5` series), from the terminating identity
  `Γ(3/2+n)/(Γ(3/2)Γ(n+1)) = Σ_k (-1)^k (4k+1) (1/2)_k² (-n)_k/(k!² (3/2+n)_k)`, certificate
  `G = (2k+1)²/((2n+2k+3)(4k+1)) F`, and Carlson's theorem at `n = -1/2`.
- **Guillera 2002** (arXiv:math/0503345): first WZ proof of R1, via
  `G(n,k) = (-1)^{n+k} (20n+2k+3) C(2k,k)² C(2n,n)² C(4n-2k,2n-k) / (2^{10n} 2^{2k} C(2n,k) C(n+k,n))`,
  `F = 64 n²/((4n-2k-1)(20n+2k+3)) G`, Zeilberger's `H(n,k) = F(n+1,n+k) + G(n,n+k)` and Carlson's theorem.
- **Guillera 2006** (arXiv:1104.0392), Table 1: generators for Bauer/R2 (`j=1`), R2/`20n+3` (`j=2`), `20n+3`/R1 (`j=4`),
  constants fixed by Carlson's theorem and the limit `k → 1/2`. Section "Some more series" prints, without certificates,
  `Σ_n C(2n,n)³ C(2k,k)² (6n+4k+1) / (2^{8n} 2^{4k} C(n+k,n)²) = 4/π` — this is exactly the family used here for R2.

## What is different here
The pairs used in this project have the shape `F = n h`, `G = Q h` with `F(0,k) = 0`, and the telescoping
is done on the finite rectangle `n < N`, `k < K` for integers only. The value of the series is the limit of the
single boundary term `G(0,K)` (Wallis' product, resp. Euler's sine product at `1/4`). Carlson's theorem,
analytic continuation in `k`, and limits of infinite sums depending on a parameter are not used.
