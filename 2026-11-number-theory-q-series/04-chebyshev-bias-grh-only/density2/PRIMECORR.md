# Prime-side correlation hypotheses: an exact duality, a conditional theorem, and the barrier

AI-produced (Anthropic's Claude), 2026-11-27. Labels as in `monotonicity/density/DENSITY.md`:
**Theorem / Proposition / Lemma / Corollary** means a complete proof is given here, modulo
the quoted standard inputs. **Remark (heuristic)** means not proved. The main statement
(`ν_{1/2}({0}) = 0` under GRH alone) is **not** proved here. Nothing here concerns Problem 12 of
arXiv:2407.03530. Status: **not yet refereed**.

This answers the parent session's request (2026-11-27): find a conditional theorem whose extra
input is a **prime-side** correlation hypothesis instead of LI, and state the barrier.

## 0. Summary

1. **Exact duality (Theorem P1, GRH + NCZ).** Take a window `w ≥ 0` with `ŵ ≥ 0`, and the
   windowed prime sum `Π_T(y) = Σ_n χ(n)Λ(n)n^{−1/2}·Tw(T(y − log n))` (window of relative width
   `≍ 1/T` around `x = e^y`). Let `M_k(T)` be its log-averaged `k`-th moment. Then
   `(−1)^kM_k(T)` equals the *independent-phase prediction* plus a **nonnegative** sum over the
   integer relations of length `≤ k` among the zeros. The prediction is computable from finite
   prime sums (one-level data).
   - In particular the odd log-averaged moments are `≤ 0`, and vanish iff there is no relation
     of odd length `≤ k`.
   - The moments up to order `L` match the prediction iff there is no relation of length `≤ L`.

   So at this level a prime-side correlation hypothesis is a zero-side relation hypothesis,
   rewritten through the explicit formula: **not a hypothesis of a different kind.**
2. **Conditional theorems.**
   - **Theorem P2 (GRH + NCZ + H_corr).** Suppose that at one scale `T ≥ C L^{3/2}log²L` the
     first `L` log-averaged moments exceed the prediction by less than `c_0^ℓ`. Then every atom
     of `ν_M` has mass `≤ C(1+|c−1|)L^{−1/2}`. The infinite version (Corollary P3) is
     equivalent to LI with simple zeros, so it gives no atoms but is not weaker than LI.
   - **Theorem P7 (GRH + NCZ + band decoupling), §3.** Take mixed moments of two windowed prime
     sums with disjoint frequency bands. Their excess over the "decoupled" prediction counts
     only the relations that involve the band. Suppose that for infinitely many `L` some
     **nonempty** band at height `s`, with `s/log s = o(L)`, is free in this sense. Then
     `ν_M` has no atoms. (For bounded `s` the hypothesis already implies TS.)
   - This prime-side hypothesis is strictly weaker than LI with simple zeros. For spectra
     with simple zeros, it contains the graded LI of `INVLO.md` (Theorem C) as a special
     case. Without simplicity it does not: (D) forbids multiple zeros in the band, while GLI
     allows them.
3. **Why Hardy–Littlewood cannot supply it (Remark P4, heuristic with a computation).** Two
   reasons:
   - The needed moments concern windows of relative width `1/T` with `T` fixed as `x → ∞`, with
     logarithmic averaging. There the error *bounds* of any pointwise Hardy–Littlewood
     asymptotic (`O(x^{1/2+ε})` per shift tuple, as in Montgomery–Soundararajan) exceed the
     main term by a factor `≍ (x/T)^{k/2}x^{−1/2+ε} → ∞`.
   - The relations that an atom forces (`INVLO.md` Prop. E) have lengths up to
     `L ≍ T/(α²δ log T)` at height `T`. A relation of length `ℓ'` changes no moment of order
     `< ℓ'`. So hypotheses on moments of any fixed order `K`, with any precision, cannot see
     them. One needs moments of order `≍ T/log T` at scale `T`, at the precision of single
     relations.

   The density question sits at the level of **exact** long-range prime correlations of very
   high order: those that detect individual relations among zeros.
4. **The spectral data cannot decide (Proposition P5, plus the sibling's Theorem 9.22 draft).**
   LI spectra are dense: some perturbation of each ordinate by less than any prescribed `δ(γ)`
   makes the spectrum `ℚ`-independent, with an absolutely continuous law. The predecessor's
   Theorem 9.22 (draft, not refereed) gives cyclic integer spectra with an atom, the counting
   law and `σ² < 1/2`. If it holds, the data "integer positive multiplicities, counting law,
   `σ² < 2`" do not decide the question.
   - No property shared by those models and stable under small perturbations can imply
     "no atom".
   - A proof must use information that the models of Theorem 9.22 lack. That information may be
     stable (asymptotic statistics of the true zeros that those lattice spectra violate, such
     as the zero-counting function) or unstable (the explicit formula as an exact identity,
     which excludes the cyclic models: `LINE1.md` §1).
5. **Literature** (§6): no derivation of the existence of race densities from
   prime-correlation hypotheses was found.

## 1. Windowed prime sums

*Window.* Fix `v ∈ C_c^∞(−½,½)`, even, `v ≥ 0`, `∫v = 1`, and put `w := v∗v`. Then:
- `w` is even, `w ≥ 0`, `supp w ⊂ [−1,1]` and `∫w = 1`;
- `ŵ(ξ) := ∫w(u)e^{iξu}du = v̂(ξ)²` is even, Schwartz, `≥ 0`, with `ŵ(0) = 1`;
- `ŵ > 0` on `(−π, π)`, since `v̂(ξ) = ∫v(u)cos(ξu)du > 0` for `|ξ| < π`.

Put `c_0 := min_{|ξ|≤1} ŵ(ξ) > 0`.

*Windowed sum.* For `T > 0` and `y > 1/T`,

    Π_T(y) := Σ_{n≥1} χ(n)Λ(n) n^{−1/2} · T w(T(y − log n)),

a finite sum over `e^{y−1/T} ≤ n ≤ e^{y+1/T}`. It is a weighted sum of `χ(n)Λ(n)` over a window of
relative width `≍ 1/T` around `x = e^y`.

*The explicit formula* (Weil; e.g. Iwaniec–Kowalski, *Analytic Number Theory*, Thm 5.12, quoted
from memory; only its structure is used). For `χ` primitive and real of parity `a`, and
`g ∈ C_c^∞(ℝ)` with `h(r) := ∫g(u)e^{iru}du`:

    Σ_ρ h(γ_ρ) = (1/2π)∫ h(r) Ω_χ(r) dr − Σ_n Λ(n)n^{−1/2}(χ(n)g(log n) + χ(n)g(−log n)),

with `ρ = ½ + iγ_ρ` the nontrivial zeros (with multiplicity) and
`Ω_χ(r) = log(q/π) + Re(Γ'/Γ)((½+a+ir)/2)`, smooth with derivatives `≪ log(|r|+2)`.

Take `g(u) := Tw(T(u − y))`. Then `h(r) = e^{iry}ŵ(r/T)`, and `g(−log n) = 0` for `y > 1/T`. Under
GRH + NCZ the zeros are `½ ± iγ`, `γ ∈ Γ`, each with multiplicity `m_γ`, since `χ` is real. Hence

    Π_T(y) = Ω_T(y) − G_T(ι(y)),     G_T := Σ_{γ∈Γ} m_γ ŵ(γ/T)(e_γ + e_{−γ}),     (1.1)

where `Ω_T(y) := (2π)^{−1}∫e^{iry}ŵ(r/T)Ω_χ(r)dr`.
- `G_T` is a continuous function on the hull `A`: the series converges absolutely.
- `Ω_T(y) = O_{T,A}(y^{−A})`, as the Fourier transform of a Schwartz function.

*Numerical check of (1.1)* (`primecorr_checks.py`, `primecorr_checks_results.txt`). For
`χ_{−4}`, with the 511 certified zeros up to height `640` and the primes up to `5.7·10⁵`, take
`T = 4, 8` and `y = 7, 9, 11, 13`. Both sides agree to `≤ 5.4·10^{−4}`, for values of size
`0.2`–`5`. This confirms the sign and normalisation used here. The residual comes from the
discretised window and the truncated zero sum.

**Lemma P0 (log-averaged moments; GRH + NCZ).** For all `T > 0` and `k ≥ 1`,

    M_k(T) := lim_{Y→∞} (1/Y)∫_{y_0}^Y Π_T(y)^k dy   exists and equals   (−1)^k ∫_A G_T^k dλ_A,

where `y_0 := 1 + 1/T`, so that (1.1) holds on the range of integration.

*Proof.* By (1.1), `Π_T(y)^k − (−G_T(ι(y)))^k → 0`, because `Ω_T → 0` and `G_T` is bounded.
Time averages of `f∘ι`, for `f ∈ C(A)`, converge to `∫_A f dλ_A`: the Kronecker flow on `A` is
uniquely ergodic, as in (I1) of `DENSITY.md`. ∎

`M_k(T)` is a log-averaged `k`-point correlation of `χ(n)Λ(n)` over windows of relative width
`≍ 1/T`.

## 2. The duality and the conditional theorem

*Multiset convention.* List the zeros with positive ordinate, with multiplicity, as a set
`Z⁺`, and put `a_ρ := ŵ(γ_ρ/T) ≥ 0`, so that `G_T = Σ_{ρ∈Z⁺}a_ρ(e_{γ_ρ} + e_{−γ_ρ})`.
- A **multiset relation** is a nonzero `κ ∈ ℤ^{(Z⁺)}` with `Σ_ρ κ_ργ_ρ = 0`; its length is
  `‖κ‖_1`. Two distinct copies of a multiple zero give the relation `e_ρ − e_{ρ'}`, of length
  `2`.
- The **independent-phase model** is `G_T^{ind} := Σ_{ρ∈Z⁺} 2a_ρ cos θ_ρ`, with `θ_ρ` independent
  and uniform. Put `μ_k(T) := E[(G_T^{ind})^k]`.

**Theorem P1 (duality; GRH + NCZ).** For all `T > 0` and `k ≥ 1`,

    (−1)^k M_k(T) = μ_k(T) + Δ_k(T),     Δ_k(T) := Σ' Π_{i≤k} a_{ρ_i} ≥ 0.

The sum `Σ'` runs over sequences `(ρ_i, σ_i)_{i≤k} ∈ (Z⁺ × {±1})^k` whose net vector
`κ_ρ := Σ_{i:ρ_i=ρ}σ_i` is a multiset relation; it has length `≤ k` and `≡ k (mod 2)`.
Consequently:
- (a) `Δ_1 = 0`, and `Δ_2(T) = Σ_γ 2m_γ(m_γ−1)ŵ(γ/T)²`: only multiple zeros contribute.
- (b) For odd `k`, `μ_k = 0`. So `M_k(T) = −Δ_k(T) ≤ 0`, with equality iff no multiset relation of
  odd length `≤ k` exists among the `ρ` with `a_ρ > 0`.
- (c) For each `L`: `Δ_ℓ(T) = 0` for all `ℓ ≤ L` iff no multiset relation of length `≤ L` exists among
  the `ρ` with `a_ρ > 0`. For all `T` outside a countable set, `a_ρ > 0` for every `ρ`.
- (d) Every multiset relation of length `ℓ` supported on zeros with `γ_ρ ≤ T` contributes at least
  `c_0^ℓ` to `Δ_ℓ(T)`.
- (e) The prediction `μ_k(T)` is prime-side. It is a polynomial in
  `s_{2j}(T) := Σ_{ρ∈Z⁺}ŵ(γ_ρ/T)^{2j}` (`2j ≤ k`), and by the explicit formula at `y = 0`,

      2s_{2j}(T) = (2π)^{−1}∫ŵ(r/T)^{2j}Ω_χ(r)dr − 2Σ_{n ≤ e^{2j/T}} Λ(n)χ(n)n^{−1/2}·T w^{∗2j}(T log n),

  a finite prime sum.

*Proof.*
- *Expansion.* `E_{λ_A}[G_T^k] = Σ Π_i a_{ρ_i}·E[e_{Σ_iσ_iγ_{ρ_i}}]`, summed over all sequences, and
  `E[e_g] = 1_{g=0}`. In the independent model the expectation of `Π_i e^{iσ_iθ_{ρ_i}}` is
  `1_{κ=0}`. Since `κ = 0` implies `Σσ_iγ_{ρ_i} = 0`, the difference is the sum over sequences
  with `κ ≠ 0` and `Σ_ρκ_ργ_ρ = 0`. All terms are `≥ 0`. The sum converges absolutely, since the
  total over all sequences is `(2Σ_ρa_ρ)^k < ∞`.
- *The net vector.* Each step changes `‖κ‖_1` by `±1`, so `‖κ‖_1 ≤ k` and `‖κ‖_1 ≡ k (mod 2)`.
- *(a)* No relation has length `1`. A length-`2` relation `σ_1γ_{ρ_1} + σ_2γ_{ρ_2} = 0` with
  `κ ≠ 0` needs opposite signs and `γ_{ρ_1} = γ_{ρ_2}` with `ρ_1 ≠ ρ_2`. This gives
  `2m_γ(m_γ−1)` ordered sequences per ordinate.
- *(b)* `μ_k = 0` by the symmetry `θ ↦ θ + π`.
  - If there is no odd relation of length `≤ k`, then every term of `Σ'` vanishes, since
    `‖κ‖_1 ≡ k` is odd.
  - Conversely, let `κ` be a relation of odd length `ℓ ≤ k`. Realise it by a sequence of length
    `ℓ` as in (c). Then append `(k−ℓ)/2` cancelling pairs `(ρ,+),(ρ,−)` with `a_ρ > 0`
    (*padding*); the net vector is still `κ`. This gives a term of `Σ'` of positive weight, so
    `Δ_k > 0`.
  - Padding really occurs. For the ordinates `{1, 2, √2}`, the only odd relations of length
    `≤ 5` have length `3`, yet `Δ_5 = 161/90 > 0` (referee 4, `referee4/r4_moments_exact.py`).
- *(c)* A relation `κ` of length `ℓ ≤ L` is realised by a sequence of length `ℓ`: list each `ρ`
  `|κ_ρ|` times, with sign `sgn κ_ρ`. Its weight is `Π a_ρ^{|κ_ρ|} > 0`. Conversely, without such
  relations every term of `Σ'` vanishes.
  - `a_ρ = 0` iff `γ_ρ/T` is a real zero of `v̂`, an entire function `≢ 0`. These zeros are
    discrete, so each `ρ` excludes countably many `T`.
- *(d)* The same sequence has weight `≥ c_0^ℓ`, since `ŵ ≥ c_0` on `[−1,1]`.
- *(e)* The cumulants of `G^{ind}` are `κ_r = Σ_ρ(2a_ρ)^r κ_r(cos θ)`, and they vanish for odd
  `r`. Moments are polynomials in cumulants. For `s_{2j}`, apply the explicit formula with
  `g(u) = Tw^{∗2j}(Tu)` and `h(r) = ŵ(r/T)^{2j}`. Here `g` is even with support in
  `|u| ≤ 2j/T`, and the sum over all zeros equals `2s_{2j}`, by NCZ and the symmetry. ∎

**Hypothesis H_corr(L, T)** (prime side). For `2 ≤ ℓ ≤ L`,
`(−1)^ℓ M_ℓ(T) − μ_ℓ(T) < c_0^ℓ`. Both terms are prime-side: `M_ℓ(T)` by definition, and `μ_ℓ(T)`
by P1(e). In words, the log-averaged `ℓ`-point correlations of `χ(n)Λ(n)` at relative scale
`1/T` exceed their one-level ("zero-diagonal") prediction by less than `c_0^ℓ`, for `ℓ ≤ L`.

**Theorem P2 (GRH + NCZ + H_corr; quantitative).** Fix `M > 0` and `χ`, and let
`r_γ(M) := 4M/(M+iγ)` as in `DENSITY.md` §1. There are constants `C_1 = C_1(M,χ)` and
`C_2 = C_2(M,χ)` such that: if H_corr(L,T) holds for some `L ≥ 2` and
`T ≥ C_1L^{3/2}log²(L+3)`, then `ν_M({c}) ≤ C_2(1+|c−1|)L^{−1/2}` for every `c ∈ ℝ`.

*Proof.*
1. *No short relations below `T`.* By P1(d) and H_corr, there is no multiset relation of length
   `ℓ ∈ [2, L]` supported on zeros with `γ_ρ ≤ T`, and none has length `1`.
   - Hence the zeros with `γ ≤ T` are simple: length `2`.
   - The ordinates in `Γ_T` satisfy no integer relation of length `≤ L`: a relation `k` on
     `Γ_T` gives a multiset relation of the same length, choosing one zero per ordinate.
   - So `S_L(Γ_T) = ∅` in the notation of `INVLO.md`.
2. *Small ball.* Apply `INVLO.md` Proposition B′ to `Γ' := Γ_T`, `c_γ := m_γr_γ(M)`, with
   `S_L = ∅`, `B = 0`, `σ_Y² ≤ 2V_M`, and `D := 6/(π²|r_{γ_3}(M)|)` (Lemma D(a), the three
   smallest ordinates; `T ≥ γ_3`). With `N := ⌊L/2⌋` and `R := 16(|c−1| + 2(2NV_M)^{1/2})`,

       λ{|F_{M,T} − c| ≤ ρ} ≤ 352·D·R/(N+1) + 56·16^{−N},     ρ := 2R/(N+1).

3. *Truncation* (as in `INVLO.md` Proposition B):
   `ν_M({c}) ≤ λ{|F_{M,T} − c| ≤ ρ} + 4V_M^{>T}/ρ²`.
   - `ρ ≥ 32(2V_M)^{1/2}N^{−1/2}`.
   - `V_M^{>T} ≤ C_χM² log²(T+3)/T`, since `m_γ ≪ log γ` and `N(t) ≪ t log t`.

   So `4V_M^{>T}/ρ² ≤ C N log²(T+3)/T ≤ C L^{−1/2}` once `T ≥ C_1 L^{3/2}log²(L+3)`.
4. Adding up gives `ν_M({c}) ≤ C_2(1+|c−1|)L^{−1/2}`. ∎

**Corollary P3 (GRH + NCZ).**
- (i) Suppose H_corr(L_j, T_j) holds along `L_j → ∞`, with `T_j ≥ C_1L_j^{3/2}log²(L_j+3)`. Then
  `ν_M` has no atoms, and the conclusions of `HYP.md` Theorem 1.3 hold. In particular, for
  `χ_{−4}` the logarithmic density of `{π(x;4,3) > π(x;4,1)}` exists.
- (ii) The hypothesis of (i) holds (for some such sequence) iff all zeros are simple and the
  ordinates are `ℚ`-linearly independent.

*Proof.* (i) H_corr does not depend on `M`, and H_corr(L,T) implies H_corr(L',T) for `L' ≤ L`.
The hypothesis of (i) is stated with `C_1 = C_1(M_0,χ)` for one `M_0`. Given any `M`, put

    L'_j := min(L_j, ⌊(C_1(M_0)/C_1(M))^{2/3} L_j⌋).

Then `T_j ≥ C_1(M)L'^{3/2}_j log²(L'_j+3)`. Theorem P2 at `M` gives
`ν_M({c}) ≤ C_2(M)(1+|c−1|)L'^{−1/2}_j → 0`. So `ν_M` has no atoms for every `M`, as `HYP.md`
Theorem 1.3 requires (`M = m + ½`, all `m > −½`). This argument is due to referee 4.

(ii):
- ⇐: under LI and simplicity, `Δ ≡ 0` by P1(c).
- ⇒: a multiset relation `κ` of length `ℓ` has finite support. For `j` with `L_j ≥ ℓ` and
  `T_j ≥ max_ρ γ_ρ` over the support, P1(d) contradicts H_corr. ∎

*What is new and what is not.* H_corr(L,T) is a finite list of prime-side numbers: moments of
order `≤ L` at one scale. But each `M_ℓ(T)` is a limit as `x → ∞`, so H_corr is still an
asymptotic statement about primes at all heights.
- By Theorem P1, H_corr implies two facts: the zeros below `T` are simple, and they satisfy no
  relation of length `≤ L`.
- It also constrains zeros well above `T`, because of padding. A double zero `γ_0 > T`
  contributes `Δ_ℓ ≥ ŵ(γ_0/T)²(2s_2)^{(ℓ−2)/2}` for even `ℓ`. Already at `ℓ = 2`, a double
  zero anywhere in `(T, πT]` violates H_corr(2,T) (referee 4, `referee4/r4_P5_P4.py`).
- So the effective tolerance is much smaller than `c_0^ℓ`.

Proposition B of `INVLO.md` already gives the bound `atom ≤ C/√L` from LI up to length `L`.
What Theorem P2 adds is that the hypothesis is prime-side and stated at one finite scale.

## 3. Beyond LI on the prime side: band decoupling (Theorem P7)

Plain moments of one windowed sum cannot allow low relations: pairings amplify them (see the
remark at the end of this section). Mixed moments of two windowed sums with **disjoint
frequency supports** can. This needs windows whose Fourier transforms have compact support.

*Status.* This section was read by referee 4 (`REFEREE-4b.md`), and the corrections from that
read are included. In particular, the first version of Theorem P7 held vacuously for an empty
band; condition (N) below fixes that.

*Numerical checks.*
- `p7_checks.py` uses a plateau-free bump window on `[9, 21]`, which contains four zeros, and
  omits the term `Π^−` (negligible for `y ≥ 6`). The truncated prime sums match the band's zero
  sum to `≤ 7·10^{−3}` at truncation `V = 2` and `≤ 4·10^{−4}` at `V = 4`. The agreement is not
  monotone in `V` at `y = 9`.
- Referee 4 checked the regularised identity of the proof to `≤ 1.4·10^{−7}`
  (`referee4/r4b_P0prime.py` (R1)).
- Referee 4 also checked a plateau window of Theorem P7 type to `≤ 7·10^{−3}`, limited by the
  slow conditional convergence, and the one-level count (R2, R3).

**Lemma P0′ (explicit formula for frequency-compact windows; GRH + NCZ).** Let `β ∈ C_c^∞(ℝ)` be
even and real, with `supp β ⊂ [−S,S]`. Let `b(v) := (2π)^{−1}∫β(r)e^{−irv}dr`, a real even
Schwartz function. Then the limit

    Π_β(y) := lim_{X→∞} Σ_{n≤X} χ(n)Λ(n)n^{−1/2} b(y − log n)

exists for every `y`. Moreover `Π_β(y) = −G_β(ι(y)) + o(1)` as `y → ∞`, where
`G_β := Σ_{γ∈Γ} m_γ β(γ)(e_γ + e_{−γ})` is a trigonometric polynomial. Without NCZ, a zero of
multiplicity `n_0` at the central point would add the term `−n_0β(0)`.

*Proof.* **Existence of the limit.** Under GRH, `ψ_χ(t) := Σ_{n≤t}χ(n)Λ(n) ≪ √t log²t`. For
Schwartz `f`, partial summation gives

    Σ_{n≤X}χ(n)Λ(n)n^{−1/2}f(log n) = ψ_χ(X)X^{−1/2}f(log X) − ∫_1^X ψ_χ(t)(d/dt)[t^{−1/2}f(log t)]dt.

The boundary term tends to `0`. In the variable `u = log t`, the integrand is
`≪ (1+u)²(|f| + |f'|)(u)`, which is integrable. So the limit exists.

**Regularisation.** Put `h_ε := (β e^{i·y}) ∗ φ_ε`, with `φ_ε` the centred Gaussian density of
standard deviation `ε ≤ 1`.
- `h_ε` is entire, with Gaussian decay in every horizontal strip. So it belongs to the
  standard admissible class of the explicit formula: `h` holomorphic in `|Im r| ≤ ½ + δ`, with
  `|h(r)| ≪ (1+|r|)^{−1−δ}` (e.g. Iwaniec–Kowalski Thm 5.12, quoted from memory).
  Alternatively, reduce to `C_c^∞` by a smooth cutoff at `|u| ≤ U`: the cutoff error on the prime
  side is `≪ ∫_U^∞ e^{u/2−ε²u²/2}du → 0`.
- Its inverse transform is `g_ε(u) = b(u − y)e^{−ε²u²/2}`.

Let `ε → 0` in the explicit formula for `h_ε`, using three dominations, all uniform in `ε ≤ 1`.
- *Zero side.* `|h_ε(γ)| ≤ ‖β‖_∞` for `|γ| ≤ max(2, 2S)`. Otherwise `|h_ε(γ)| ≤ Ce^{−γ²/8}`,
  because `ε ↦ ε^{−1}e^{−γ²/(8ε²)}` increases on `(0,1]` when `|γ| ≥ 2`.
- *`Ω`-side.* `|Ω_χ(r)| ≪ log(|r|+2)`.
- *Prime side.* After partial summation, the integrand is dominated by
  `C(1+u)²[|b'(u−y)| + (1+|u|)|b(u−y)|]`. This majorant is integrable and independent of `ε`.
  (Referee 4 supplied these majorants, `REFEREE-4b.md` A2.)

**The identity.** Under NCZ the limit is `Π_β(y) = Ω_β(y) − G_β(ι(y)) − Π_β^−(y)`, where:
- `Ω_β(y) := (2π)^{−1}∫β(r)e^{iry}Ω_χ(r)dr → 0`, being the Fourier transform of a `C_c^∞`
  function;
- `Π_β^−(y) := lim Σ_{n≤X}χ(n)Λ(n)n^{−1/2}b(y + log n)`. By the same bound it is
  `≪ ∫_0^∞ (1+u)²(|b|+|b'|)(y+u)du → 0`. ∎

**Setting for Theorem P7.** Fix:
- `L ≥ 2`, `N := ⌊L/2⌋`, and `T* ≥ 1`;
- a band `0 < s_1 < s_2 < T*` and a margin `g ∈ (0,1)`;
- even windows `β_I, β_R ∈ C_c^∞(ℝ)` with values in `[0,1]`, such that:
  - `β_I = 1` on `[s_1+g, s_2−g]`, and `supp β_I ∩ [0,∞) ⊂ [s_1, s_2]`;
  - `β_R = 1` on `[0, s_1−g] ∪ [s_2+g, T*]`, and `supp β_R ∩ [0,∞) ⊂ [0, s_1] ∪ [s_2, T*+1]`.

Write `Π_I := Π_{β_I}` and `Π_R := Π_{β_R}`.

**Hypothesis H_dec(L; s_1, s_2, g, T*)** (prime side).
- (Z) No zero has ordinate in the open margins `(s_1−g, s_1+g) ∪ (s_2−g, s_2+g)`.
- (N) The band contains a zero: `n_I ≥ 1`. Here `n_I := Σ_{s_1<γ<s_2} m_γ` is the number of zeros
  in the band, counted with multiplicity.
- (D) For all `a ≥ 0` and `b ≥ 1` with `a + b ≤ L`, the log-averaged joint moment
  `M_{a,b} := lim_{Y→∞}(1/Y)∫_{1}^{Y}Π_R^aΠ_I^b dy` satisfies

      (−1)^{a+b}M_{a,b} − (−1)^a M_{a,0}·μ_b(n_I) < 1,

  where `μ_b(n) := E[(Σ_{i≤n}2cos θ_i)^b]` for independent uniform `θ_i`.

*Why the margins can be open.* A zero at `s_1−g` or `s_2+g` has `β_R = 1` and `β_I = 0`. A zero at
`s_1+g` or `s_2−g` has `β_I = 1` and `β_R = 0`.

*On the "prime side".*
- (Z) says that the one-level sum `Σ_ρ η(γ_ρ)` vanishes, where `η ≥ 0` is an even `C_c^∞` bump
  that is positive exactly on the open margins. By the explicit formula this is an identity
  between a Gamma-factor integral and a prime sum.
- Given (Z), `n_I = ½Σ_ρβ_I(γ_ρ)` is a one-level sum, and (N) says that it exceeds `½`.

So every quantity in H_dec is a limit of prime sums.

**Lemma P1′ (band duality; GRH + NCZ).** Assume (Z), and put `I := Γ ∩ (s_1, s_2)`. For all `a, b`,

    E_{a,b} := (−1)^{a+b}M_{a,b} − (−1)^a M_{a,0}μ_b(n_I) ≥ 0.

Here `E_{a,b}` is a sum of nonnegative weights over the sequences of length `a + b` (`a` from the
zeros in `supp β_R`, `b` from the zeros in `I`) whose net multiset vector `κ = (κ_R, κ_I)` is a
relation with `κ_I ≠ 0`.
- A relation supported on zeros of height `≤ T*` has weight exactly `1`.
- In particular `E_{0,2} = Σ_{γ∈I}2m_γ(m_γ−1)`, so (D) at `(a,b) = (0,2)` forces the zeros in
  the band to be simple.

*Proof.* By Lemma P0′ and unique ergodicity, `(−1)^{a+b}M_{a,b} = E_{λ_A}[G_R^aG_I^b]`, where
`G_I = Σ_{γ∈I}m_γ(e_γ + e_{−γ})` by (Z) and the shape of `β_I`. Expand as in Theorem P1. In the
two terms, each sequence contributes the indicators

    1[total frequency sum = 0]   and   1[R-part frequency sum = 0]·1[κ_I = 0].

- If `κ_I = 0`, the band part has frequency sum `0`, so the two indicators coincide.
- If `κ_I ≠ 0`, the second vanishes.

Hence the difference is `Σ_{κ_I≠0, total=0}(weights) ≥ 0`. The weight of a sequence is
`Π β_R(γ)` over its R-terms. By (Z) and the shape of `β_R` it equals `1` when all R-terms have
height `≤ T*`. The term `(−1)^aM_{a,0}μ_b(n_I)` is `E[G_R^a]·E[(G_I^{ind})^b]`, where `G_I^{ind}` is
a sum of `n_I` independent `2cos θ`. Two copies of a multiple band zero form the multiset
relation `e_ρ − e_{ρ'}`, which gives the formula for `E_{0,2}`. ∎

**Theorem P7 (GRH + NCZ + band decoupling).** Suppose that along a sequence `L_j → ∞` the
hypothesis H_dec(L_j; s_{1,j}, s_{2,j}, g_j, T*_j) holds, including (N), with:
- `T*_j := ⌊L_j/2⌋⁴`;
- `2s_{1,j} ≤ s_{2,j} ≤ 3s_{1,j}`;
- `s_{1,j}/log(s_{1,j}+3) = o(L_j)`.

Then `ν_M` has no atoms, for every `M > 0`, and the conclusions of `HYP.md` Theorem 1.3 hold.

*Proof.* Pass to a subsequence along which either `s_{1,j}` stays bounded or `s_{1,j} → ∞`.

*Case A: `s_{1,j}` bounded.* The bands then lie in a bounded region and, by (N), each contains a
zero. So some zero `γ_0` lies in `I_j` for infinitely many `j`. For those `j`, Lemma P1′ and (D)
exclude every multiset relation of length `≤ L_j`, supported on zeros of height `≤ T*_j`, that
involves `γ_0`. Every relation has finite length and finite height, and `L_j, T*_j → ∞`. So `γ_0`
lies in no relation at all: it is simple, and `γ_0 ∉ ⟨Γ∖{γ_0}⟩_ℚ`. Hence TS holds
(`DENSITY.md` §3), and `DENSITY.md` Theorem 3.1 gives absolute continuity of every `ν_M`.

*Case B: `s_{1,j} → ∞`.* Fix `j` and drop it from the notation.
1. *The band is free of short relations.* By Lemma P1′ and (D), there is no multiset relation of
   length `≤ L` supported on zeros of height `≤ T*` that involves a zero of `I`. In terms of
   ordinates, `S_L(Γ_{T*}) ∩ I = ∅`.
2. *Small ball.* Apply `INVLO.md` Proposition B′ to `Γ' := Γ_{T*}` and `c_γ := m_γr_γ(M)`, with:
   - `B' := Σ_{γ≤T*}|c_γ| ≤ C_χM log²(T*+3)` and `σ' := (2V_M)^{1/2}`;
   - `D` from Lemma D(b), with the free set `I`. By (Z) and Riemann–von Mangoldt, for large `s_1`
     the band `[s_1+g, s_2−g]` contains `≍ s_1 log s_1` zeros with multiplicity, and `≫ s_1`
     distinct ordinates. As in the proof of `INVLO.md` Theorem C, this gives
     `D ≤ C(M,χ)(1 + (s_1/log(s_1+3))^{1/2})`.

   This yields `λ{|F_{M,T*} − c| ≤ ρ} ≤ 352·D·R/(N+1) + 56·16^{−N}`, with
   `R/(N+1) ≪ (1 + log²T*)/N + N^{−1/2}`.
3. *Truncation.* As in Proposition B, `ν_M({c}) ≤ that + 4V_M^{>T*}/ρ²`. With `T* = N⁴` the last
   term is `≪ N·log²N/N⁴ → 0`.
4. *Limit.* With `D ≤ C(1 + (s_1/log(s_1+3))^{1/2})` and `R/(N+1) ≪ log²L/L + L^{−1/2}` (as
   `log T* ≍ log L`),

       ν_M({c}) ≤ C(c,M,χ)·(1 + (s_1/log(s_1+3))^{1/2})·(log²L/L + L^{−1/2}) + o(1).

   Under `s_1/log(s_1+3) = o(L)` both products tend to `0`. Indeed,
   `(s_1/log s_1)^{1/2}L^{−1/2} → 0`, and `(s_1/log s_1)^{1/2}log²L/L ≤ o(L^{1/2})·log²L/L → 0`. ∎

So the content of Theorem P7 lies in bands with `s_1 → ∞`. For bounded `s_1` the hypothesis
already implies TS.

*Comparison* (abstract spectra; the examples are from referee 4, `REFEREE-4b.md` A11–A13).
- **LI with simple zeros implies H_dec.** Then `E ≡ 0`. For large `s_1`, bands with zero-free
  margins (choose the `s_i` in gaps between ordinates) contain zeros. So H_dec holds for every
  `L`, with bands at any large height.
- **GLI and H_dec.**
  - Assume GLI(ψ), with `ψ(t)log t/t → ∞`, and in addition that the zeros are simple in one
    admissible band per scale. Then H_dec holds. Indeed, all relations of length `≤ L` lie
    below `T_L` (`INVLO.md` Theorem C), a band just above `T_L` is then free, and
    `T_L/log T_L = o(L)`.
  - So for spectra with simple zeros, `INVLO.md` Theorem C is a special case of Theorem P7.
  - Without simplicity it is not. GLI concerns distinct ordinates, whereas a double zero in the
    band gives `E_{0,2} = 4`. Example: a `ℚ`-independent spectrum with all multiplicities `2`
    satisfies GLI(ψ) for every ψ, but fails (D) for every nonempty band.
- **For simple spectra, H_dec is strictly weaker than GLI.** Take the spectrum of `REFEREE-4b.md`
  A13:
  - `Γ = {nα} ∪ ⋃_jB_j ∪ {z_x, x+z_x : x ∈ ⋃_jB_j}`;
  - `B_j` is a band of simple elements in `(s_{1,j}+g, s_{2,j}−g)`, and `z_x ∈ (T*_j+1, T*_j+2)`;
  - the elements of `⋃B_j ∪ {z_x} ∪ {α}` are `ℚ`-independent;
  - the scales are sparse (`s_{1,j+1} > T*_j + 3 + s_{2,j}`) and `L_j ≥ j·s_{1,j}`.

  Every relation through `x ∈ B_j` uses `z_x` and `x + z_x`, which lie outside `supp β_R`. So
  `E_{a,b} = 0` and H_dec holds along `L_j`. But GLI(ψ) fails for every `ψ > 3`: the set contains
  `(n+1)α − nα − α = 0`.
- **Devin's Theorem 5.4(1) and H_dec are incomparable.**
  - *5.4(1) does not imply H_dec.* Double every element of a set satisfying 5.4(1). That
    condition concerns `ℚ`-spans, so it survives, while (D) fails at `(0,2)` for every nonempty
    band.
    - With simple elements: start from `DENSITY.md` Example 3.3a, in which every element lies
      in a relation. Add `ℚ`-independent triples `{u_k, v_k, u_k+v_k}` with some `u_k+v_k` in
      every interval of length `½` above height `1`.
    - Then every admissible band with `s_1 > 2.5` either violates (Z) or contains some
      `u_k+v_k`, which lies in a relation of length `3`, so that `E ≥ 1`. Every band with bounded
      `s_1` fails (D) for large `L_j`.
  - *H_dec does not imply 5.4(1).* The spectrum of the previous item satisfies H_dec. At
    `T = T*_j + 3 + s_{2,j}`, every element `≤ T^{1/2}` lies in the `ℚ`-span of the other elements
    of `Γ_T`:
    - `nα` via `2nα`;
    - `x` via `(x+z_x) − z_x`;
    - `z_x` via `(x+z_x) − x`.

    So 5.4(1) fails. That spectrum satisfies TS, so it is not a new case of "no atoms".

*Remark (why plain moments do not suffice).* With one positive window, a relation `κ` of
length `ℓ` among low zeros is allowed by H_dec. It contributes to the excess of every moment of
order `k ≡ ℓ (mod 2)`, `k > ℓ`, at least its weight times `(Σ_ρ a_ρ²)^{(k−ℓ)/2}`: append
cancelling pairs. So a threshold on plain moments forbids it. The decoupled prediction removes
exactly those contributions, since sequences with `κ_I = 0` cancel. Referee 4's exact check
(`referee4/r4b_P1prime.py` (Q2)) confirms this. Relations among R-zeros contribute nothing to
`E`, while the plain moment excess of the combined window is `6, 0, 460` at orders `3`–`5`.

## 4. The barrier on the prime side

**Remark P4 (heuristic; a computation, not a theorem).**

(a) *Hardy–Littlewood does not reach these moments.* Use a natural average over `x ∈ [X, 2X]` and
a sharp window of length `h := X/T`, for simplicity. The `k`-th moment of the window sum of
`χ(n)Λ(n)` is

    (1/X) Σ_{d ∈ [0,h]^{k−1}} (h − max d) Σ_{n ≤ X} Π_i χ(n+d_i)Λ(n+d_i).

Assume a twisted Hardy–Littlewood asymptotic `S_χ(d)X + E(d)`, uniform for shifts `≤ X`, with
`|E(d)| ≤ X^{1/2+ε}`. This is the form used by Montgomery–Soundararajan for `ψ`.
- The main term is of Gaussian size `≍ (h log(X/h))^{k/2}`.
- The error terms add up to `≪ X^{−1}·h·h^{k−1}·X^{1/2+ε} = h^k X^{−1/2+ε}`.
- Their ratio is `≍ (X/T)^{k/2}X^{−1/2+ε}(log T)^{−k/2} → ∞` as `X → ∞`, for every `k ≥ 2` and
  fixed `T`.

So pointwise Hardy–Littlewood asymptotics with power-saving error bounds say nothing about
`M_k(T)` at fixed `T`. Errors of size `X^{1−k/2}(T log T)^{k/2}` per shift tuple would be
needed, i.e. `O(1)`, `O(X^{−1/2})`, `O(X^{−1})`, … for `k = 2, 3, 4, …`. This is not expected for
individual tuples, whose errors are believed to be of size `≍ X^{1/2}`. What the moments see is
the *average* of the error terms `E(d)` over long ranges of shifts.

Three caveats (referee 4):
- The identification of that average with the relation structure of the zeros below `≍ T`
  presupposes that the shift-average of the Hardy–Littlewood main terms equals `μ_k(T)` to the
  precision of `Δ_k`. Theorem P1 says nothing about Hardy–Littlewood terms.
- Theorem P1 concerns logarithmic averages with smooth windows. The computation above uses
  natural averages with sharp windows. At fixed `T`, the natural-average moments of
  `S(x)/√x` do not converge: they are almost periodic in `log X`. So only orders of magnitude
  are compared.
- The displayed moment formula ignores the multiplicity of the position of the minimum among
  the `n_i`, and `Σ_{n≤X}` stands for `n ∈ [X,2X]`. Both are harmless for orders of magnitude.

The proved Gaussian-moment range of Montgomery–Soundararajan is short windows (`h ≤ X^{1/k}` for
the `k`-th moment). In zero terms that is zeros of height `≫ X^{1−1/k}`, the opposite end. The
classical second-moment duality (Goldston–Montgomery: pair correlation versus primes in short
intervals) is also of that kind. At our scale the second moment is trivial (P1(a)).

(b) *Fixed-order Gaussianity cannot see what an atom needs.* Unlike (a), this part is not
heuristic, except its last sentence. It is a rigorous consequence of Theorem P1 and
`INVLO.md` Proposition E. By Theorem P1, a relation `κ` of length `ℓ'` contributes:
- `0` to `Δ_ℓ` if `ℓ < ℓ'` or `ℓ ≢ ℓ' (mod 2)`;
- `(ℓ'!/Π_ρ|κ_ρ|!)·Π_ρa_ρ^{|κ_ρ|}` at `ℓ = ℓ'`, which lies between `c_0^{ℓ'}` and `ℓ'!` for zeros
  of height `≤ T`;
- by padding, a multiple `≍ (T log T)^{(ℓ−ℓ')/2}` of its weight for `ℓ > ℓ'`, `ℓ ≡ ℓ'`.

An atom of mass `α` forces saturation (Proposition E): for every `δ`, all but a proportion
`O(δ)` of the zeros in the dyadic blocks up to height `C_2^{−1}α²δ·L log L` lie in relations of
length `≤ L`. At height `T` this concerns lengths `L ≳ C_2T/(α²δ log T)`. For example,
`α = 0.0891`, `δ = 0.1`, `T = 10³` and `C_2 = 1` give `L ≥ 1.1·10⁵` (referee 4,
`referee4/r4_P5_P4.py`).

So the saturating relations may all be far longer than any fixed moment order `K`. Then they
change **no** moment of order `≤ K` at all. A hypothesis "the moments of order `≤ K` are
Gaussian up to a relative error `η`", with `K` fixed or `K = o(T/log T)`, cannot exclude the
structure an atom needs, however small `η` is. *(Heuristic:)* it would take moments of order
`≍ T/log T` at scale `T`, at the precision of single relations. That is the precise sense in
which the density question lies at the level of **exact** prime correlations of very high
order.

## 5. The spectral data cannot decide

**Proposition P5 (LI spectra are dense; no hypothesis).** Let `Γ ⊂ (0,∞)` be countable and
discrete, `m : Γ → ℤ_{≥1}`, and `δ : Γ → (0, ½)` arbitrary. Then there is an injection `γ ↦ γ'`
with `|γ' − γ| < δ(γ)·min(1, γ/2)` such that `Γ' := {γ'}` is linearly independent over `ℚ`.
Suppose moreover that `Σm_γ²/γ² < ∞`. Then for every `M > 0`, the model
`F' := 1 + Re Σ_γ m_γr_{γ'}(M)e_{γ'}` (with at least three terms) has an absolutely continuous
law with a bounded continuous density. Moreover:
- if `(Γ,m)` satisfies the counting law `Σ_{γ∈(t,t+1]}m_γ ≤ A log(t+3)`, then `(Γ',m)` satisfies
  it with `2A` in place of `A`, after adjusting the additive constant inside the logarithm;
- `σ_M²` changes by a factor in `[(1+δ_∞)^{−2}, (1−δ_∞)^{−2}]`, where `δ_∞ := sup δ`.

*Proof.*
- *Construction.* Enumerate `Γ = {γ_1, γ_2, …}`. Choose `γ'_j` in the open interval of
  allowed values, outside the countable set `⟨γ'_1, …, γ'_{j−1}⟩_ℚ`. Then each `γ'_j` lies
  outside the `ℚ`-span of the earlier ones, so `Γ'` is independent.
- *The law.* By Kronecker's theorem the closure torus is `𝕋^{Γ'}`, and its Haar measure is the
  product measure. So `F' − 1 = Σ m_γ|r_{γ'}|cos(θ_γ + arg r_{γ'})`, with independent uniform
  `θ_γ`, converging in `L²`.
  - Its characteristic function is `Π_γ J_0(ξm_γ|r_{γ'}|)`, a convergent infinite product.
  - It is bounded by any three of its factors, hence integrable (`INVLO.md` Lemma D(a)).
  - So the law has a bounded continuous density.
- *Counting law.* `|γ'−γ| < 1/2`, so `(t,t+1]` receives only zeros from `(t−½, t+3/2]`, which
  is covered by two unit intervals.
- *`σ²`.* `|γ'−γ| < δ_∞γ/2` gives `(M²+γ²)/(M²+γ'²) ∈ [(1+δ_∞)^{−2}, (1−δ_∞)^{−2}]`. ∎

**Corollary P6 (barrier; conditional on the predecessor's Theorem 9.22).** Assume Theorem 9.22
and Corollary 9.23(iii) of `THM922_DRAFT.md` (predecessor branch, **draft, not refereed**).
These give, for `0 < θ ≤ 1` and a small `M ≤ M_0(θ)`, a model at the classical weight `½`. Put
`Ω := 1/(2M)`. The frequencies are `Ωℕ`, with positive integer multiplicities
`m_n ≤ C(M,θ)(log(n+2))^θ`, and:
- the counting law holds with a constant `A_c` (astronomically large; draft §8, item 9);
- `σ² < ½`;
- `λ{F = 0} > 0`.

Then the data "positive integer multiplicities, the local counting law, `σ² < 2`" do not
decide whether the model's law has an atom at `0`.
- Proposition P5, applied to that spectrum, gives spectra with the same multiplicities, the
  counting law with constant `2A_c`, and `σ² < 2`, whose laws are absolutely continuous.
- They are arbitrarily close to it: `|γ'−γ| < δ(γ)` for any prescribed `δ`.

*Proof.* Combine P5 with the assumed theorem, taking `δ ≤ ¼` so that `σ²` stays below `2`. ∎

*What this shows, and what it does not* (corrected after referee 4).
- **It shows:** no property of spectra that is shared by the models of Theorem 9.22 and stable
  under small perturbations of the ordinates can imply "no atom", or characterise atoms on a
  class containing those models.
- **It does not show** that a proof must use unstable information. A proof uses some property
  `Q` of the true zeros with "`Q` + data ⇒ no atom". The models of Theorem 9.22 only have to
  fail `Q`. So `Q` may be stable:
  - asymptotic statistics of the true zeros that those lattice spectra violate, e.g. the
    zero-counting function `N(t) = (t/2π)log(qt/2πe) + O(log t)`, which the perturbations of
    P5 preserve up to `O(log t)` and the draft does not assert;
  - or local spacing statistics.

  `Q` may also be unstable: the explicit formula as an exact identity. It excludes the cyclic
  spectra outright (`DENSITY.md` Remark 7.6; `LINE1.md` §1).
- To force unstable information, one would need atom-carrying spectra arbitrarily close to
  every spectrum with the true zeros' known properties. That converse density statement is not
  shown.

*Remark.* The same barrier with real weights is `DENSITY.md` §7 (Proposition 7.3 and its
models). With Theorem 9.22 it extends to integer multiplicities obeying the counting law. Any
proof must use a property of the true zeros that these models lack.

## 6. Literature check (web search, 2026-11-27; abstracts and snippets only)

- Montgomery–Soundararajan, *Primes in short intervals*, Comm. Math. Phys. 2004
  ([arXiv math/0409258](https://arxiv.org/abs/math/0409258)): Gaussian moments of
  `ψ(x+h) − ψ(x) − h` from a uniform Hardy–Littlewood conjecture with power-saving error, for
  short `h`. This is the short-window, natural-average regime of Remark P4(a).
- Lemke Oliver–Soundararajan, *Unexpected biases in the distribution of consecutive primes*
  ([arXiv:1603.03720](https://arxiv.org/abs/1603.03720)): Hardy–Littlewood heuristics for
  consecutive-prime biases, not for the existence of race densities.
- Prime-race results without full LI all use **zero-side** hypotheses:
  - Devin ([1706.06394](https://arxiv.org/abs/1706.06394),
    [1909.03975](https://arxiv.org/abs/1909.03975));
  - Martin–Ng ([1710.00088](https://arxiv.org/abs/1710.00088));
  - Bailleul ([2007.05763](https://arxiv.org/abs/2007.05763));
  - Fiorilli ([1210.6946](https://arxiv.org/abs/1210.6946)).

  No result deriving the existence of the logarithmic density from prime-correlation
  hypotheses was found.
- *An effective Linear Independence conjecture for the zeros of the Riemann zeta function*
  ([arXiv:2311.04860](https://arxiv.org/abs/2311.04860)): a quantitative **strengthening** of
  LI (lower bounds for nonzero linear forms in the ordinates). It is not a weakening like GLI
  (`INVLO.md`).
- Best–Trudgian, *Linear relations of zeroes of the zeta-function*
  ([arXiv:1209.3843](https://arxiv.org/abs/1209.3843)): numerical exclusion of small relations
  among the first zeros of `ζ`. This is the kind of statement that GLI and H_corr(L,T)
  formalise.
- Not read: Rubinstein–Sarnak's discussion of dependence among zeros, Kaczorowski–Ramaré, and
  the full texts above. The duality P1 is presumably folklore in some form; **no novelty is
  claimed.**
