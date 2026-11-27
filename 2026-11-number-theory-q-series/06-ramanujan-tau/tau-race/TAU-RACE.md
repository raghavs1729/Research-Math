# The Cesàro mean of Ramanujan's τ over primes is positive for every x > 5.52 (GRH for L(s, Δ) only)

Working notes, 2026-11-24 to 26. Everything in `monotonicity/tau-race/` was produced for this note; nothing
outside it was modified. All review mentioned below was done by AI agents (Anthropic's Claude); no human
has refereed this. Status of every statement is given in §9.

---------------------------------------------------------------------------

## 0. Notation and conventions

    τ(n)            Ramanujan's function,  Δ(z) = Σ τ(n) e^{2πinz} = e^{2πiz} Π (1 − e^{2πinz})^{24}
    λ(n)            = τ(n) n^{-11/2};   λ(p) = α_p + β_p,  α_p β_p = 1,  |α_p| = |β_p| = 1  (Deligne)
    c_k(p)          = α_p^k + β_p^k ∈ [−2, 2];  c_1 = λ(p),  c_2 = λ(p)² − 2,  c_{k+1} = λ(p) c_k − c_{k−1}
    L(s, Δ)         = Σ λ(n) n^{-s} = Π_p (1 − α_p p^{-s})^{-1} (1 − β_p p^{-s})^{-1}     (analytic normalisation)
    Λ(s)            = (2π)^{-s} Γ(s + 11/2) L(s, Δ) = Λ(1 − s)      (entire; Hecke)
    L_S(s)          = L(s, Sym²Δ) = Π_p [(1 − α_p² p^{-s})(1 − p^{-s})(1 − β_p² p^{-s})]^{-1} = Σ a_S(n) n^{-s}
    γ_S(s)          = Γ_R(s+1) Γ_C(s+11) = π^{-(s+1)/2} Γ((s+1)/2) · 2 (2π)^{-(s+11)} Γ(s+11)
    Λ_S(s)          = γ_S(s) L_S(s) = Λ_S(1 − s)                    (entire; Shimura 1975)
    Θ(t)            = Σ_{p≤t} λ(p) log p
    Θ_1(t)          = ∫_0^t Θ(u) du = Σ_{p≤t} λ(p) log p · (t − p)       (continuous, piecewise linear)
    D_0(x), D_1(x)  = Σ_{p≤x} λ(p),  Σ_{p≤x} λ(p) p
    R(x)            = Σ_{p≤x} λ(p) (1 − p/x) = D_0(x) − D_1(x)/x = (1/x) ∫_0^x D_0(t) dt
    G(t), N(x)      = Θ_1(t)/t^{3/2},   R(x) log x/√x
    ψ_Δ, ψ_S, ψ     = Σ_{n≤t} Λ_Δ(n),  Σ_{n≤t} Λ_S(n),  Σ_{n≤t} Λ(n)   with
                      Λ_Δ(p^k) = c_k(p) log p,   Λ_S(p^k) = (c_{2k}(p) + 1) log p,   Λ = von Mangoldt
    θ(y)            = Σ_{p≤y} log p

"GRH(Δ)": every zero of `L(s, Δ)` in `0 < Re s < 1` has real part 1/2. "GRH(Sym²Δ)": the same for `L_S`.
"RH_S(H)": every zero `ρ` of `L_S` with `0 < Re ρ < 1` and `|Im ρ| ≤ H` has real part 1/2.
No hypothesis on ζ, no linear independence, no Deep Riemann Hypothesis.

"Certified" means: decided in Arb ball arithmetic (python-flint 0.9.0) or exact integers, every
comparison certain, correct given FLINT/Arb and the program named. "Numerical" means floating point.

Facts used without proof (all classical):

 (F1) Deligne: `|α_p| = 1`.
 (F2) Hecke: `Λ(s)` is entire, `Λ(s) = Λ(1−s)`, bounded in vertical strips.
 (F3) Functional equation of the symmetric square (Rankin 1939; Shimura 1975; Gelbart–Jacquet 1978):
      `Λ_S = γ_S L_S`, `γ_S(s) = Γ_R(s+1)Γ_C(s+k−1)`, `k = 12`, is entire, `Λ_S(s) = Λ_S(1−s)`, and `L_S` is of finite
      order. Sources actually read for this note (quotations, page numbers and the bookkeeping in
      `SOURCES-CHECK.md`):
      (i) D. Zagier, "Modular forms whose Fourier coefficients involve zeta-functions of quadratic fields",
      Modular Functions of One Variable VI, LNM 627 (1977): for a normalised eigenform f of weight k on
      `SL_2(ℤ)` and `D_f(s) = Π_p (1−α_p²p^{-s})^{-1}(1−α_pᾱ_p p^{-s})^{-1}(1−ᾱ_p²p^{-s})^{-1}` (`α_pᾱ_p = p^{k−1}`), eq. (4), p. 109:
      "`D_f^*(s) = 2^{−s} π^{−3s/2} Γ(s) Γ((s−k+2)/2) D_f(s) = D_f^*(2k−1−s)`" (attributed there to Rankin, with
      meromorphic continuation), and p. 111: "Corollary 1 (Shimura): The function `D_f(s)` defined by (2) has
      a holomorphic continuation to the whole complex plane." With `s = s' + k − 1`: `D_f(s) = L_S(s')` and
      `D_f^*(s'+k−1) = (1/2) π^{−(k−2)/2} Γ_R(s'+1)Γ_C(s'+k−1) L_S(s')`, a constant multiple of `Λ_S(s')`, and
      `2k−1−s ↔ 1−s'`. `Λ_S` is entire: `D_f` is entire and `Γ(s)Γ((s−k+2)/2)` is holomorphic for `Re s > k−2`, so
      `D_f^*` is holomorphic there, hence, by the functional equation, also for `Re s < k+1`.
      (ii) O. Balkanova, D. Frolenkov, Algebra Number Theory 12 (2018), (2-4)–(2-5), p. 36 (weight `2k`, so `k = 6`
      for Δ): "`L_∞(s) := π^{−3s/2} Γ((s+1)/2) Γ((s−1)/2 + k) Γ(s/2 + k)`. Shimura [1975] showed that the completed
      L-function `Λ(sym² f, s) := L_∞(s) L(sym² f, s)` is entire and satisfies the functional equation
      `Λ(sym² f, s) = Λ(sym² f, 1 − s)`"; this `L_∞` is a constant multiple of `γ_S` (duplication formula).
      Not read: Shimura (Proc. London Math. Soc. 31 (1975) 79–98), Gelbart–Jacquet, Rankin; in particular the
      finite order of `L_S` (it follows from Rankin's integral against the Eisenstein series, or from
      automorphy on GL(3)) is cited, not checked. Numerical evidence, labelled as such: PARI's `lfuncheckfeq`
      (−67 bits; −111 bits in the second referee's independent set-up, −2 bits with a wrong gamma shift), and
      `Φ(1/y) = yΦ(y)` (§5) checked directly to `10^{-83}` at three points in ball arithmetic (`selftest.py`, test F).
      Theorem 4 and everything that uses it depend on this functional equation with exactly this `γ_S`.
 (F4) Hadamard factorisation for `Λ` and `Λ_S` (entire of order 1, real on the real axis, symmetric):
      `(Λ'/Λ)(s) = B + Σ_ρ (1/(s−ρ) + 1/ρ)` with `Re B = −Σ_ρ Re(1/ρ)`, hence
      `Re (Λ'/Λ)(s) = Σ_ρ Re 1/(s−ρ)` (absolutely convergent), and `N(T+1) − N(T) ≪ log(|T|+2)`
      (Iwaniec–Kowalski, *Analytic Number Theory*, Theorem 5.6, Proposition 5.7; both L-functions are in
      their class). Neither function vanishes at `s = 0, 1`: `Λ(1)`, `Λ_S(1)` are certified nonzero
      (`data/constants.json`: `L(1,Δ) = 0.8393…`; `data/sym2_rh_H60.json`: `L_S(1) = 0.6317…`).
 (F5) The bounds for `L'/L` needed to move contours: for `−1/2 ≤ σ ≤ 2`,
      `(L'/L)(s) = Σ_{|s−ρ|<1} 1/(s−ρ) + O(log(|s|+2))` (IK Prop. 5.7(2)); to the left of the strip,
      `L'/L(s) = −L'/L(1−s) − (γ'/γ)(s) − (γ'/γ)(1−s) ≪ log(|s|+2)` at distance `≥ 1/4` from the trivial zeros.
 (F6) Unconditional bounds for θ: `|θ(y) − y| < 0.2 y/log² y` for `y ≥ 3 594 641`. Sources read (quotations in
      `SOURCES-CHECK.md`): Broadbent–Kadiri–Lumley–Ng–Wilk, Math. Comp. 90 (2021), arXiv:2002.11068v2,
      Corollary 11.2, eq. (3.39) ("`|ϑ(x) − x| < M_k x/log^k x` for `x ≥ X_0`", row `M_2 = 0.2`, `X_0 = 3 594 641`); the same
      inequality is in P. Dusart, arXiv:1002.0442, Theorem 5.2 (row `η_2 = 0.2`, `x_2 = 3 594 641`). BKLNW remark
      that the proof of the corresponding theorem in Dusart's published paper (Ramanujan J. 45 (2018),
      Theorem 4.2) is incorrect and re-prove the corollary; so BKLNW is the reference. (Fiori–Kadiri–Swidinsky,
      arXiv:2204.02588, call some zero-density inputs used by BKLNW for very large x unreliable; no erratum to
      BKLNW was found. This is why Remark 4.3 matters.) Fallback, used only in Remark 4.3:
      Rosser–Schoenfeld (Illinois J. Math. 6 (1962)), in the form quoted in BKLNW's Tables 1–2 (the original was
      not obtained): `θ(y) ≤ 1.01624 y` for `y ≥ 1` and `θ(y) ≥ y(1 − 1/(2 log y))` for `y ≥ 563`; Remark 4.3 uses only
      the weaker `θ(y) ≥ y(1 − 1/log y)` for `y ≥ 10^7`.

---------------------------------------------------------------------------

## 1. Results

Certified constants (provenance in §7):

    σ_0²    = 2 (Λ'/Λ)(1) = 0.112252140471908526745583426253…       (= Σ_ρ |ρ|^{-2} under GRH(Δ))
    a_0     = (L'/L)(0, Δ) = 0.170657847591640…,      a_{−1} = (L'/L)(−1, Δ) ∈ [0.281796, 0.281803]
    (Λ_S'/Λ_S)(2) ∈ [0.2941460, 0.2941565],           C_S := (12/5)(Λ_S'/Λ_S)(2) ≤ 0.70598
    K_S     ≤ 0.34176   (an upper bound for (L_S'/L_S)(−2))
    c_ψ     = √2:  ψ(y) ≥ y − √2·√y for 0 < y ≤ Y_1 = 10^7;   θ(y) < y for 0 < y ≤ 10^7
    ε_1     = 0.2/log²(10^7) = 7.69844…·10^{-4},  c_θ = 1 + ε_1
    T_60    ≤ 0.184605  (bound for Σ_{|Im ρ|>60} 2/|ρ(ρ+2)| over ALL zeros of L_S in the strip; Lemma 5)
    t*      = 4.266791956347541…   (zero of Θ_1 in (3,5):  t* = (2λ(2) log 2 + 3λ(3) log 3)/(λ(2) log 2 + λ(3) log 3))
    x*      = 5.518054631888274…   (zero of R in (5,7):    x* = D_1(5)/D_0(5))

Define for `t ≥ 1`

    E_odd(t) = (3/4) t^{-1/6} + (5/6) t^{-3/10} + t^{-5/14} · max(0, (log t/log 2 − 5)/2),
    L_A(t)   = (2/3)(1 − ε_1·1_{t > 10^14}) − σ_0² − a_0 t^{-1/2} − ((4/5)√2 + C_S) t^{-1/4} − 2 t^{-1} − K_S t^{-3/2} − 2 c_θ E_odd(t),
    L_B(t)   = L_A(t) − T_60.

`L_A`, `L_B` are nondecreasing on `[527, 10^14]` and on `(10^14, ∞)`. Values (`data/bounds.txt`):
`L_A(10^6) ≥ 0.2120`, `L_A(10^7) ≥ 0.3483`, `L_A(10^14+) ≥ 0.5458`, `lim L_A ≥ 0.5539`;
`L_B(10^6) ≥ 0.0274`, `L_B(10^7) ≥ 0.1637`, `L_B(10^14+) ≥ 0.3612`, `lim L_B ≥ 0.3692`.

**Theorem 1 (integrated θ-sum of τ).**
(a) Assume GRH(Δ). Then `Θ_1(t) ≥ L_B(t) t^{3/2}` for every `t ≥ 1`, and

    Θ_1(t) = Σ_{p≤t} τ(p) p^{-11/2} log p · (t − p) > 0     for every t > t* = 4.2667…

(and `Θ_1(t) ≤ 0` for `t ≤ t*`).
(b) Assume GRH(Δ) and GRH(Sym²Δ). Then moreover `Θ_1(t) ≥ L_A(t) t^{3/2}` for every `t ≥ 1`.

**Theorem 2 (Cesàro mean).** Assume GRH(Δ). Then

    R(x) = Σ_{p≤x} τ(p) p^{-11/2} (1 − p/x) = (1/x) ∫_0^x ( Σ_{p≤t} τ(p) p^{-11/2} ) dt > 0     for every x > x* = 5.5180…

(and `R(x) ≤ 0` for `x ≤ x*`). More precisely, for every `x ≥ 97`,

    R(x) ≥ Θ_1(x)/(x log x) + 0.6719 ≥ L_B(x) √x/log x + 0.6719,

and `R(x) > 0.15 √x/log x` for every `x ≥ 11`. If also GRH(Sym²Δ) holds, `L_B` may be replaced by `L_A`,
and `R(x) > √x/(3 log x)` for every `x ≥ 100`.

**Theorem 3 (equivalence).** The following are equivalent:
(a) GRH(Δ); (b) `R(x) > 0` for every `x > x*`; (c) `R(x) ≥ 0` for all sufficiently large x;
(d) `Θ_1(t) > 0` for every `t > t*`; (e) `Θ_1(t) ≥ 0` for all sufficiently large t.

**Theorem 4 (computer-assisted; no hypothesis).** `RH_S(60)` holds: `L(s, Sym²Δ)` has exactly 41 zeros
with `0 < Im s ≤ 60` and `0 < Re s < 1`; they are simple and have real part 1/2 (the first is
`1/2 + i·7.0702241402…`), and `L_S` has no zero on the real segment `(0, 1)`.

Theorems 1(a), 2 and the implications (a)⟹(b), (a)⟹(d) of Theorem 3 use Theorem 4; without it they hold
with the hypothesis "GRH(Δ)" replaced by "GRH(Δ) and RH_S(60)". The converse implications of Theorem 3 and
Theorem 1(b) do not use Theorem 4.

*Sign.* The sum is biased towards **positive** values: the prime-square term of `−L'/L(s,Δ)` is
`Σ_p (λ(p)² − 2) log p · p^{-2s}`, and `λ(p)² − 2 = λ(p²) − 1` has mean `−1` because `L_S` has no pole. This is the
bias of Sarnak's letter and of Koyama–Kurokawa, in a form valid at every x. Sizes: bias `2/3`, zeros of
`L(s,Δ)` at most `σ_0² = 0.1123`, zeros of `L_S` above height 60 at most `0.1846` (and of lower order
under GRH(Sym²Δ)).

*Numerical truth* (certified, `data/small_range_10000000.json`: minimum of the gap bounds over all gaps
`[q, p]` of consecutive primes with `p > a`, so that the whole of `[a, 9 999 991]` is covered): `G(t) ≥ 0.6257`
for all real `t ∈ [10^4, 9 999 991]`, `G(t) ≥ 0.5814` on `[100, 9 999 991]`; `N(x) ≥ 0.7951` on `[100, 9 999 991]`,
`N(x) ≥ 0.1718` on `[11, 9 999 991]`. (Float, `data/explore_float_1e7.txt`: the minimum of `G` over primes in
`[10^4, 10^7]` is `0.625775` at `p = 9 132 733`; `2/3 − σ_0² = 0.5544`.)

---------------------------------------------------------------------------

## 2. The decomposition and the three explicit formulas

### Lemma 1 (decomposition; uses only (F1))

For `u > 0`,

    Θ(u) = ψ_Δ(u) − [ψ_S(√u) − ψ(√u)] − Σ_{k≥3 odd} Θ^{(k)}(u^{1/k}),      Θ^{(k)}(y) := Σ_{p≤y} c_k(p) log p,

and `|Θ^{(k)}(y)| ≤ 2θ(y)`. Consequently, for `t > 0`,

    Θ_1(t) = ψ_{1,Δ}(t) − J_S(t) + J_ζ(t) − Σ_{k≥3 odd} ∫_0^t Θ^{(k)}(u^{1/k}) du,                         (2.1)

    ψ_{1,Δ}(t) = Σ_n Λ_Δ(n)(t − n)_+,    J_S(t) = Σ_n Λ_S(n)(t − n²)_+,    J_ζ(t) = Σ_n Λ(n)(t − n²)_+ = ∫_0^{√t} ψ(y) 2y dy.

*Proof.* `ψ_Δ(u) = Σ_{k≥1} Σ_{p^k≤u} c_k(p) log p`. The terms with `k = 2j` even sum to
`Σ_j Σ_{p^j ≤ √u} c_{2j}(p) log p = Σ_{p^j ≤ √u} (Λ_S(p^j) − log p) = ψ_S(√u) − ψ(√u)`, because the local factor
of `L_S` at p has the roots `α_p², 1, β_p²`. The term `k = 1` is `Θ(u)`. `|c_k| ≤ 2` by (F1). Integrate over
`[0, t]`; `∫_0^t ψ_S(√u) du = Σ Λ_S(n) (t − n²)_+`, similarly for ψ. ∎

The point of this form: **all** even prime powers are carried by `ψ_S` and `ψ`, so no sign-indefinite
even power has to be bounded trivially, and the main term comes from `ψ` (not θ).

### Lemma 0 (moving the contour)

Let `𝓛` be `L(·, Δ)` or `L_S`, `m ∈ {1, 2}`, `t > 1`, `c > 1/m`. Then

    (1/2πi) ∫_{(c)} (−𝓛'/𝓛)(ms) t^{s+1}/(s(s+1)) ds = Σ (residues of the integrand at its poles in Re s < c),

where the poles are `s = ρ/m` (ρ a nontrivial zero of `𝓛`), `s = w_0/m` (`w_0` a trivial zero), `s = 0`, `s = −1`, and
the sum over ρ converges absolutely.

*Proof.* Two bounds. (i) By (F4), (F5): there are `T_k ∈ [k, k+1]` with `(𝓛'/𝓛)(σ ± iT_k) ≪ log² T_k` uniformly for
`−1/2 ≤ σ ≤ 2` (there are `≪ log k` ordinates in `[k, k+1]`, so some `T_k` has distance `≫ 1/log k` from all of
them; then each of the `≪ log k` terms `1/(s−ρ)` in (F5) is `≪ log k`). (ii) For `Re w ≤ −1/2` and distance
`≥ 1/4` from the trivial zeros, `(𝓛'/𝓛)(w) ≪ log(|w| + 2)`: by the functional equation
`(𝓛'/𝓛)(w) = −(𝓛'/𝓛)(1−w) − (γ'/γ)(w) − (γ'/γ)(1−w)`, the first term is bounded (`Re(1−w) ≥ 3/2`), and `γ'/γ` is a
constant plus digamma functions of `w + 11/2` resp. `(w+1)/2`, `(w+11)/2`, `(w+12)/2`, with `ψ(z) ≪ log(|z|+2)` at
distance `≥ 1/8` from the poles of ψ (reflection formula and Stirling).
Integrate over the rectangle with vertices `c ± iT_k/m`, `−U' ± iT_k/m`, where `mU'` is a positive integer for
`L(·,Δ)` (trivial zeros at half-integers) and `mU' ∈ ℕ + 1/2` for `L_S` (trivial zeros at integers). On the
horizontal sides `|s(s+1)| ≥ T_k²/m²` and, by (i), (ii) and the boundedness of `𝓛'/𝓛` in `Re w ≥ 2`, the integrand is
`≪ (m²/T_k²) t^{σ+1} [log² T_k + log(|σ| + 2)]`; so these two integrals are `≪_t log² T_k/T_k²`, uniformly in `U'`
(`t > 1`). On the left side, by (ii), the integral is `≪ t^{1−U'} ∫_{−∞}^{∞} log(mU' + m|y| + 2)/(U'² + y²) dy ≪ t^{1−U'} log(U'+2)/U'`.
Let `U' → ∞`, then `k → ∞`. The residues at `s = ρ/m` are `−(μ_ρ/m) t^{ρ/m+1}/((ρ/m)(ρ/m+1))`, of modulus `≤ m μ_ρ t²/|ρ|²`,
and `Σ_ρ |ρ|^{-2} < ∞` (F4); so the sum over ρ needs no grouping. ∎

### Lemma 2 (explicit formula for ψ_{1,Δ})

For `t ≥ 1`,

    ψ_{1,Δ}(t) = − Σ_ρ t^{ρ+1}/(ρ(ρ+1)) − a_0 t + a_{−1} − Σ_{n≥0} t^{-9/2−n}/((n + 11/2)(n + 9/2)),

the sum over the nontrivial zeros of `L(s,Δ)` (with multiplicity) being absolutely convergent. Under
GRH(Δ), `Σ_ρ 1/|ρ(ρ+1)| < Σ_ρ 1/|ρ|² = σ_0² = 2(Λ'/Λ)(1)`, and therefore, for `t ≥ 1`,

    ψ_{1,Δ}(t) ≥ − σ_0² t^{3/2} − a_0 t.                                                                  (2.2)

*Proof.* For `t > 1`, `ψ_{1,Δ}(t) = (1/2πi) ∫_{(2)} (−L'/L)(s, Δ) t^{s+1}/(s(s+1)) ds`. By Lemma 0
(`m = 1`; the trivial zeros of `L(s,Δ)` are the poles of `Γ(s + 11/2)`, at `s = −11/2 − n`, `n ≥ 0`, simple) the
integral is the sum of the residues to the left. Residues: at `s = ρ`: `−t^{ρ+1}/(ρ(ρ+1))`; at `s = 0` (where `L'/L` is regular since
`L(0,Δ) = Λ(1)/Γ(11/2) ≠ 0`): `−a_0 t`; at `s = −1` (regular, `L(−1,Δ) = (2π)^{-1}Λ(2)/Γ(9/2) ≠ 0`):
`(−L'/L)(−1)·t^0/(−1) = a_{−1}`; at `s = −11/2 − n`: `−t^{-9/2−n}/((−11/2−n)(−9/2−n))`. Both sides are continuous
at `t = 1`.
For the second statement, by (F4) at `s = 1`: `(Λ'/Λ)(1) = Σ_ρ Re 1/(1−ρ)`, and under GRH
`Re 1/(1−ρ) = (1/2)/|ρ|²`; `|ρ+1| > |ρ|`. Finally `Σ_{n≥0} 1/((n+11/2)(n+9/2)) = 2/9` (telescoping), so for
`t ≥ 1` the last two terms are `≥ a_{−1} − 2/9 > 0.059` (certified). ∎

Checks (`code/independent_check.gp`, PARI, 38 digits). (i) At `t = 1` the formula says
`0 = −Σ_ρ 1/(ρ(ρ+1)) − a_0 + a_{−1} − 2/9`, and by (F4) `Σ_ρ 1/(ρ(ρ+1)) = (Λ'/Λ)(1) − (Λ'/Λ)(2)`; PARI gives
`−3.9·10^{-38}` for the combination. This is only a consistency check of the bookkeeping: by the functional
equation the combination reduces identically to `ψ(11/2) − ψ(9/2) − 2/9 = 0`, whatever the values of `L'/L`; it
tests the list of trivial zeros and the two residues at `0, −1`, not the analytic input. The independent test
is (ii). (ii) Direct sum against the formula with 277 pairs of zeros:
difference `9.9·10^{-5}` at `t = 2.5`, `−1.7·10^{-4}` at `t = 5.5`, `−6.6·10^{-4}` at `t = 10`, within the tail
bound `0.0103·t^{3/2}`. (iii) `a_0 = log 2π − ψ(11/2) − σ_0²/2` (functional equation) agrees with direct
evaluation. Values of the true sums: `Σ_{|γ|<300} 1/|ρ|² = 0.10194`, `Σ_{|γ|<300} 1/|ρ(ρ+1)| = 0.10155`.

### Lemma 3 (explicit formula for J_S)

For `t ≥ 1`,

    J_S(t) = − Σ_ρ 2 t^{1+ρ/2}/(ρ(ρ+2)) − a_{S,0} t + 2 t^{1/2} + (L_S'/L_S)(−2) − Σ_{w} (m_w/2) t^{1−w/2}/((w/2)(w/2 − 1)),

where ρ runs over the nontrivial zeros of `L_S`, `a_{S,0} = (L_S'/L_S)(0)`, and `w` runs over the integers
`w ≥ 3` such that `−w` is a trivial zero of `L_S`, with multiplicity `m_w`
(`w ∈ {3,5,7,9}`: `m = 1`; `w ≥ 11` odd: `m = 2`; `w ≥ 12` even: `m = 1`). All terms of the last sum are positive.
Consequently, with `Z_S(t) := |Σ_ρ 2 t^{1+ρ/2}/(ρ(ρ+2))|`:

    J_S(t) ≤ Z_S(t) + 2 t^{1/2} + K_S       for t ≥ 1,  provided RH_S(2) holds.                           (2.3)

*Proof.* For `t > 1` and `c > 1/2`, `J_S(t) = (1/2πi) ∫_{(c)} (−L_S'/L_S)(2s) t^{s+1}/(s(s+1)) ds`
(termwise: `(1/2πi)∫ (t/n²)^s t/(s(s+1)) ds = t(1 − n²/t)_+`). By Lemma 0 (`m = 2`) the integral is the sum of the
residues to the left. Residues, with `w = 2s`: at `s = ρ/2`: `−(1/2) t^{ρ/2+1}/((ρ/2)(ρ/2+1)) = −2t^{1+ρ/2}/(ρ(ρ+2))`;
at `s = 0` (`L_S(0) = Λ_S(1)/γ_S(0) ≠ 0`): `−a_{S,0} t`; at `s = −1/2` (the simple trivial zero `w = −1`):
`−(1/2) t^{1/2}/((−1/2)(1/2)) = 2t^{1/2}`; at `s = −1` (`L_S(−2) = Λ_S(3)/γ_S(−2) ≠ 0`, regular):
`(−L_S'/L_S)(−2)·t^0/(−1) = (L_S'/L_S)(−2)`; at `s = −w/2`, `w ≥ 3`: `−(m_w/2) t^{1−w/2}/((−w/2)(1−w/2))`, and
`(−w/2)(1 − w/2) = (w/2)(w/2 − 1) > 0`. The trivial zeros are the poles of `γ_S`: `Γ((s+1)/2)` gives
`s = −1, −3, −5, …`, `Γ(s+11)` gives `s = −11, −12, −13, …`.
For (2.3): the trivial-zero sum is subtracted and positive; `(L_S'/L_S)(−2) = −(γ_S'/γ_S)(−2) − (γ_S'/γ_S)(3) − (L_S'/L_S)(3)`
and `|(L_S'/L_S)(3)| ≤ −3ζ'/ζ(3)`, which gives `K_S ≤ 0.34176` (certified). It remains to see
`a_{S,0} ≥ 0`: by the functional equation `a_{S,0} = −(γ_S'/γ_S)(0) − (Λ_S'/Λ_S)(1)`, with
`−(γ_S'/γ_S)(0) = 1.04024…`, and `(Λ_S'/Λ_S)(1) = Σ_ρ Re 1/(1−ρ) ≤ 3 Σ_ρ Re 1/(2−ρ) = 3(Λ_S'/Λ_S)(2) ≤ 0.88247`: for a zero
on the line, `[ (1/2)/(1/4+γ²) ] / [ (3/2)/(9/4+γ²) ] = (9/4+γ²)/(3(1/4+γ²)) ≤ 3`; for a zero `β + iγ` with `|γ| > 2`,
`a ↦ a/(a²+γ²)` increases on `(0, |γ|)`, so `Re 1/(1−ρ) = (1−β)/((1−β)²+γ²) ≤ (2−β)/((2−β)²+γ²) = Re 1/(2−ρ)`.
Hence `a_{S,0} ≥ 0.1577 > 0` (certified, `aS0_lower`). ∎

Checks (PARI). (As in Lemma 2, the identity at `t = 1` reduces, by the functional equation, to the digamma
identity `(γ_S'/γ_S)(0) − (γ_S'/γ_S)(−2) + 2 = Σ_w (m_w/2)/((w/2)(w/2−1)) = 1 + 1/9 + 1/10`; it tests the list of trivial
zeros, their multiplicities and the residues at `0, −1/2, −1`, nothing more. The independent tests are the
direct sums against the formula in `AUDIT-NUMERICS.md`, item 6, and in `REFEREE-REPORT.md`.) At `t = 1`: `0 = −[(Λ_S'/Λ_S)(1) − (Λ_S'/Λ_S)(3)] − a_{S,0} + 2 + (L_S'/L_S)(−2) − Σ_w (m_w/2)/((w/2)(w/2−1))`;
with the last sum truncated at `w = 4000` PARI gives `7.5013·10^{-4}`, and the truncated tail is
`Σ_{w>4000} (3/2)·2/w² ≈ 7.5·10^{-4}`. Values: `a_{S,0} = 0.94103`, `(L_S'/L_S)(−2) = −0.22825`.

### Lemma 4 (zeros of L_S, low part)

If `RH_S(H)` holds for some `H ≥ 2`, then

    Σ_{|Im ρ| ≤ H} 2/|ρ(ρ+2)|  ≤  (12/5) (Λ_S'/Λ_S)(2)  =  C_S ≤ 0.70598,

and under GRH(Sym²Δ) the same bound holds for the sum over all zeros.

*Proof.* For `ρ = 1/2 + iγ`: `2/|ρ(ρ+2)| = 2/√((1/4+γ²)(25/4+γ²))` and `Re 1/(2−ρ) = (3/2)/(9/4+γ²)`. The quotient is
`(4/3)(9/4+g)/√((1/4+g)(25/4+g))` with `g = γ²`; it is `≤ 12/5` (its value at `g = 0`), since
`25(g+9/4)² ≤ 81(g+1/4)(g+25/4)` ⇔ `0 ≤ 56g² + 414g`. All terms of `Σ_ρ Re 1/(2−ρ) = (Λ_S'/Λ_S)(2)` are positive. ∎

`(Λ_S'/Λ_S)(2) = (γ_S'/γ_S)(2) + (L_S'/L_S)(2)` is computed from the Euler product
`−(L_S'/L_S)(s) = Σ_p log p · [ (c_2 x − 2x²)/(1 − c_2 x + x²) + x/(1−x) ]`, `x = p^{-s}`, `c_2 = τ(p)²/p^{11} − 2` (exact
rationals for `p ≤ 10^7`), with the tail `|…| ≤ 3 Σ_{p>P} log p/(p²−1) ≤ 3.03 (log P + 1)/P` for `P ≥ 10` and `Re s = 2`: the bracket is
`Σ_{k≥1} (α_p^{2k} + 1 + β_p^{2k}) x^k`, of modulus `≤ 3|x|/(1−|x|) = 3/(p²−1)`; `1/(n²−1) ≤ 1.01/n²` for `n ≥ 11`; and
`Σ_{n>P} log n/n² ≤ ∫_P^∞ log u/u² du = (log P + 1)/P` because `log u/u²` decreases for `u > √e`.

### Lemma 5 (zeros of L_S above height H; unconditional)

Let `H > H' > 0`. For the zeros `ρ = β + iγ` of `L_S` in the critical strip,

    Σ_{|γ|>H} 2/γ²  ≤  T(H, H') := 2c [ (3/2)(log(H'/2π) + 1)/H' + (Σ_j x_j²)/(3H'³) + (Σ_j 1/(24 x_j))/(2H'²) + 6 log ζ(2)/H'² ],

    c = √2 / arctan((H − H')/√2),      (x_1, x_2, x_3) = (3/2, 13/2, 7).

In particular `T_60 := T(60, 51.45) ≤ 0.184605`, and for every `t ≥ 1`

    | Σ_{|γ|>60} 2 t^{1+ρ/2}/(ρ(ρ+2)) |  ≤  T_60 · t^{3/2}.

*Proof.* Put `F(T) = Re (Λ_S'/Λ_S)(2 + iT) = Σ_ρ (2−β)/((2−β)² + (T−γ)²)` (F4). For `a ∈ [1, 2]`,
`a/(a² + u²) ≥ 1/(2 + u²)` (⇔ `2a − a² + u²(a−1) ≥ 0`), so `F(T) ≥ Σ_ρ 1/(2 + (T−γ)²)`, every term positive.
Let `γ > H` and `d = γ − H' ≥ H − H'`. Since `T ↦ T^{-2}` is convex on `(0,∞)` and the kernel is even,

    ∫_{H'}^∞ T^{-2}/(2 + (T−γ)²) dT ≥ ∫_{−d}^{d} (γ+u)^{-2}/(2+u²) du ≥ γ^{-2} ∫_{−d}^{d} du/(2+u²) = γ^{-2} √2 arctan(d/√2) ≥ γ^{-2}·2/c.

The zeros are symmetric under `γ ↦ −γ` (real coefficients) and `F` is even; so
`Σ_{|γ|>H} 2/γ² ≤ 2c ∫_{H'}^∞ F(T) T^{-2} dT`. Now `F = Re γ_S'/γ_S(2+iT) + Re L_S'/L_S(2+iT)`.
(i) `Re (L_S'/L_S)(2+iT) = −Σ_n Λ_S(n) n^{-2} cos(T log n)`, and `|∫_{H'}^∞ cos(aT) T^{-2} dT| ≤ 2/(aH'²)` (one
integration by parts), so the contribution is at most `(2/H'²) Σ_n |Λ_S(n)|/(n² log n) ≤ (6/H'²) Σ Λ(n)/(n² log n) = 6 log ζ(2)/H'²`.
(ii) `(γ_S'/γ_S)(s) = −(3/2) log π + (1/2)[ψ((s+1)/2) + ψ((s+11)/2) + ψ((s+12)/2)]` (duplication formula). For
`z = x + iy`, `x, y > 0`, Binet's formula `ψ(z) = log z − 1/(2z) − 2∫_0^∞ u du/((u²+z²)(e^{2πu}−1))` with
`|u² + z²| ≥ Im z² = 2xy` and `∫_0^∞ u/(e^{2πu}−1) du = 1/24` gives `Re ψ(z) ≤ log|z| + 1/(24xy)`. With `z_j = x_j + iT/2`:
`log|z_j| ≤ log(T/2) + 2x_j²/T²`. Hence `Re (γ_S'/γ_S)(2+iT) ≤ (3/2) log(T/2π) + Σ_j [x_j²/T² + 1/(24 x_j T)]`, and
`∫_{H'}^∞ log(T/2π) T^{-2} dT = (log(H'/2π) + 1)/H'`.
Finally `|ρ(ρ+2)| ≥ γ²` and `|t^{1+ρ/2}| = t^{1+β/2} ≤ t^{3/2}`. The number is evaluated in `bounds.py`. ∎

(For orientation: by the closed form `Σ_ρ 1/|ρ|² = 2(Λ_S'/Λ_S)(1) = 0.198422` (PARI) and the 41 zeros of
Theorem 4, the true value of `Σ_{|γ|>60} 2/γ²` is about `0.104`, if those zeros are on the line.)

### Lemma 6 (inputs on ψ and θ; unconditional)

(a) (Certified, `cert_constants.py` part C3.) `ψ(y) ≥ y − √2 √y` and `θ(y) < y` for `0 < y ≤ Y_1 = 10^7`. (The
supremum of `(y − ψ(y))/√y` on `(0, 10^7]` is `√2`, approached as `y → 2^-`; between consecutive prime
powers `(y − ψ)/√y` increases, so the supremum over each interval is its limit at the right end, which is
what the program evaluates.)
(b) (F6.) For `y ≥ 10^7`: `(1 − ε_1) y ≤ θ(y) ≤ ψ(y)` and `θ(y) ≤ (1 + ε_1) y`.

Hence `θ(y) ≤ c_θ y` for all `y > 0`, and

    J_ζ(t) = ∫_0^{√t} ψ(y) 2y dy ≥ (2/3)(1 − ε_1·1_{t > Y_1²}) t^{3/2} − (4/5)√2 t^{5/4}        (t > 0):

for `√t ≤ Y_1` integrate (a); for `√t > Y_1`,
`J_ζ ≥ (2/3)Y_1³ − (4/5)√2 Y_1^{5/2} + (2/3)(1−ε_1)(t^{3/2} − Y_1³) ≥ (2/3)(1−ε_1)t^{3/2} − (4/5)√2 t^{5/4}`.

---------------------------------------------------------------------------

## 3. Proof of Theorem 1

Let `t ≥ 1`. Use (2.1).

* `ψ_{1,Δ}(t) ≥ −σ_0² t^{3/2} − a_0 t` by (2.2) (GRH(Δ)).
* `J_ζ(t) ≥ (2/3)(1 − ε_1 1_{t>10^14}) t^{3/2} − (4/5)√2 t^{5/4}` (Lemma 6).
* `J_S(t) ≤ Z_S(t) + 2t^{1/2} + K_S` by (2.3). Under GRH(Sym²Δ), `|t^{1+ρ/2}| = t^{5/4}` and `Z_S(t) ≤ C_S t^{5/4}`
  (Lemma 4). Under `RH_S(60)` only: the zeros with `|γ| ≤ 60` contribute at most `C_S t^{5/4}` (Lemma 4) and the
  others at most `T_60 t^{3/2}` (Lemma 5), so `Z_S(t) ≤ C_S t^{5/4} + T_60 t^{3/2}`.
* Odd `k ≥ 3`: `|∫_0^t Θ^{(k)}(u^{1/k}) du| ≤ 2 ∫_0^t θ(u^{1/k}) du ≤ 2 c_θ (k/(k+1)) t^{1+1/k}`, and the term vanishes unless
  `2^k ≤ t`. For `k = 3, 5` keep `3/4`, `5/6`; for odd `k ≥ 7` use `(k/(k+1)) t^{1+1/k} ≤ t^{8/7}`; there are at most
  `max(0, (log t/log 2 − 5)/2)` such k. Total `≤ 2 c_θ E_odd(t) t^{3/2}`.

Adding up: `Θ_1(t) ≥ L_A(t) t^{3/2}` under GRH(Δ) + GRH(Sym²Δ), and `Θ_1(t) ≥ L_B(t) t^{3/2}` under
GRH(Δ) + RH_S(60), i.e. under GRH(Δ) by Theorem 4.

*Monotonicity.* All of `t^{-1/2}, t^{-1/4}, t^{-1}, t^{-3/2}, t^{-1/6}, t^{-3/10}` decrease, and
`t^{-5/14}(log t/log 2 − 5)` decreases for `log t > 5 log 2 + 14/5`, i.e. `t > 526.2`. So `L_A`, `L_B` are nondecreasing
on `[527, 10^14]` and on `(10^14, ∞)`; since `L_B(10^6) ≥ 0.0274 > 0` and `L_B(10^14 +) ≥ 0.3612 > 0` (`bounds.py`,
Arb), `Θ_1(t) > 0` for all `t ≥ 10^6`.

*Small t.* `Θ_1` is continuous and piecewise linear with nodes at the primes; `Θ_1 = 0` on `[0, 2]`,
`Θ_1(t) = λ(2) log 2 · (t − 2) < 0` on `(2, 3]` (`λ(2) = −24/2^{11/2}`), and on `[3, 5]`,
`Θ_1(t) = Θ(3) t − (2λ(2) log 2 + 3λ(3) log 3)` with `Θ(3) = 0.2901… > 0`, which vanishes exactly at `t*`
(certified enclosure `[4.266791956347541, 4.266791956347542]`) and is positive on `(t*, 5]`. At every prime
`5 ≤ p ≤ 9 999 991` (the largest prime below `10^7`), `Θ_1(p) > 0` is certified (`cert_small.py`, check (A):
664 577 primes, exact integers `τ(p)`, 128-bit balls, no failure). By linearity between nodes, `Θ_1 > 0` on
`(t*, 9 999 991] ⊃ (t*, 10^6]`. ∎

The ranges overlap from `10^6` to `10^7`; the analytic bound `L_B` is positive from `t ≈ 7.1·10^5`
(`L_A` from `t ≈ 1.112·10^5`; grid values in `data/bounds.txt`).

---------------------------------------------------------------------------

## 4. Proof of Theorem 2

### Lemma 7 (second-order partial summation; Ingham, Johnston Lemma 2.7; cf. EXPLICIT.md Lemma 2)

Let `x > 2`, `2 ≤ T ≤ x`, and `F ∈ C²[2, x]`. Then

    Σ_{p≤T} λ(p) log p · F(p) = F(T)Θ(T) − F'(T)Θ_1(T) + ∫_2^T Θ_1(t) F''(t) dt.

(Stieltjes integration by parts twice; `Θ_1(2) = 0`.) Take `F = h`, `h(t) = (x − t)/log t`: `h(x) = 0`,
`h'(t) = −1/log t − (x−t)/(t log² t)`, and

    h''(t) = (x + t)/(t² log² t) + 2(x − t)/(t² log³ t) > 0      (1 < t ≤ x).

With `T = x`: `x R(x) = Θ_1(x)/log x + ∫_2^x Θ_1 h''`. With a general `T ≤ x`:

    ∫_2^T Θ_1 h'' = x·A_1(T) + A_2(T),
    A_1(T) = D_0(T) − Θ(T)/log T − Θ_1(T)/(T log² T),
    A_2(T) = −D_1(T) + (TΘ(T) − Θ_1(T))/log T + Θ_1(T)/log² T.

Here, unlike for `χ_4`, `Θ_1` is negative on `(2, t*)`, so the kernel term cannot simply be dropped; but its
part up to `T = 97` is positive:

    A_1(97) ∈ [0.671972273685, 0.671972273686],     A_2(97) ∈ [10.326446649, 10.326446650]       (certified).

*Proof of Theorem 2 for x ≥ 97.* By Theorem 1, `Θ_1 ≥ 0` on `[97, x]`, and `h'' > 0`, so

    x R(x) = Θ_1(x)/log x + x A_1(97) + A_2(97) + ∫_{97}^x Θ_1 h'' ≥ Θ_1(x)/log x + 0.6719 x,

i.e. `R(x) ≥ Θ_1(x)/(x log x) + 0.6719 ≥ L_B(x)√x/log x + 0.6719 > 0` (the middle term is `> 0` by Theorem 1
even where `L_B < 0`).

*x < 97.* On `[p, p')` (consecutive primes), `R(x) = D_0(p) − D_1(p)/x` is monotone, and `R` is continuous. So
`R > 0` on `[7, 9 999 991]` because `R(p) > 0` at every prime `7 ≤ p ≤ 9 999 991` (certified, check (B)). On `[5, 7]`,
`R(x) = D_0(5) − D_1(5)/x` with `D_0(5) = 0.7596… > 0`, `D_1(5) = 4.1916… > 0`: increasing, zero at
`x* = D_1(5)/D_0(5)`. On `[3, 5]`: monotone with `R(3) = −0.1767…`, `R(5) = −0.0787…`, so negative. On `(2, 3]`:
`λ(2)(1 − 2/x) < 0`. Hence `R(x) > 0` ⇔ `x > x*`.

*The constants.* For real `x ∈ [q, p]` between consecutive primes `11 ≤ q < p ≤ 9 999 991`:
`N(x) ≥ min(R(q), R(p)) · log p/√p` (`log x/√x` decreases for `x > e²`). The program certifies that the minimum
of these gap bounds is `≥ 0.1718` over all gaps with `q ≥ 11`, and `≥ 0.7951` over all gaps with `p > 100` (these
cover `[100, 9 999 991]`; the first of them is `[97, 101]`). For `x ≥ 9·10^6`:
`N(x) ≥ L_B(x) ≥ L_B(9·10^6) ≥ 0.1589`, resp. `N(x) ≥ L_A(9·10^6) ≥ 0.3435 > 1/3` under GRH(Sym²Δ). ∎

**Remark 4.1 (lead inheritance).** `D_0(x) = D_1(x)/x + R(x)`. So under GRH(Δ), for every `x > x*`:
`Σ_{p≤x} τ(p) p^{-9/2} ≥ 0 ⟹ Σ_{p≤x} τ(p) p^{-11/2} > 0`. (This is the pair `m = 0, m' = 1` of Theorem A of
`../CORE-NOTES.md` for τ, at every x.)

**Remark 4.2 (what is not claimed).** The un-smoothed sum `D_0(x) = Σ_{p≤x} τ(p)p^{-11/2}` is not one-signed
(`D_0(17) < 0`; under GRH it changes sign infinitely often unless very strong linear dependences hold among
the zeros), and nothing pointwise is proved here for the Koyama–Kurokawa weight `τ(p)p^{-6}`. Other smoothings
(`e^{-p/x}`, higher Riesz means) follow from `Θ_1 > 0` on `(t*, ∞)` by Lemma 7 for large x, but the initial
negative patch `(2, t*)` makes "every x" a separate finite computation for each weight; not done.

**Remark 4.3 (robustness of the cited input).** Replace (F6) by the classical Rosser–Schoenfeld bounds, i.e.
`ε_1 = 1/log 10^7 = 0.062042…` and `c_θ = 1.01624` in the definition of `L_A`, `L_B`. What changes is the explicit
function in Theorem 1: `L_B(10^6) ≥ 0.0230`, `L_B(9·10^6) ≥ 0.1562`, `L_B(10^14+) ≥ 0.3202`, `lim L_B ≥ 0.3284`;
`L_A(10^6) ≥ 0.2076`, `L_A(9·10^6) ≥ 0.3408`, `L_A(10^14+) ≥ 0.5048`, `lim L_A ≥ 0.5130` (`data/bounds.txt`, last lines).
Everything else survives as stated: the sign statements of Theorems 1 and 2, Theorem 3, and also the two
constants of Theorem 2, since `0.1562 > 0.15` and `0.3408 > 1/3` (the finite range does not use (F6)).

---------------------------------------------------------------------------

## 5. Proof of Theorem 4 (RH for L(s, Sym²Δ) up to height 60)

### 5.1 Evaluating Λ_S rigorously

`a_S(n)` is multiplicative and rational: `a_S(p^k) = c·a_S(p^{k−1}) − c·a_S(p^{k−2}) + a_S(p^{k−3})`,
`c = α_p² + 1 + β_p² = τ(p)²/p^{11} − 1`; `|a_S(n)| ≤ d_3(n) ≤ n²` by (F1). Let

    K(z) = 4z ∫_0^∞ u^9 exp(−πz²/u² − 2πu) du      (|arg z| < π/4),        Φ(z) = Σ_{n≥1} a_S(n) K(nz).

**Lemma V0.** (a) `∫_0^∞ K(x) x^{s−1} dx = γ_S(s)` for `Re s > −1`. (b) For `|arg z| = θ < π/4` and `0 < η < 1`:
`|K(z)| ≤ C_η |z| exp(−b_{η,θ} |z|^{2/3})`, `C_η = 4·9!/(2πη)^{10}`, `b_{η,θ} = 3π (cos 2θ · (1−η)²)^{1/3}`.
(c) `Φ` is analytic in the sector, `Λ_S(s) = ∫_0^∞ Φ(y) y^{s−1} dy` for all `s ∈ ℂ`, and `Φ(1/z) = zΦ(z)`.

*Proof.* (a) Integrate in x first: `∫_0^∞ x^s e^{-πx²/u²} dx = (1/2) π^{-(s+1)/2} u^{s+1} Γ((s+1)/2)`, then
`∫_0^∞ u^{s+10} e^{-2πu} du = (2π)^{-(s+11)} Γ(s+11)`. (b) `|exp(−πz²/u²)| = exp(−π|z|² cos 2θ/u²)`; by AM–GM,
`π|z|² cos 2θ/u² + π(1−η)u + π(1−η)u ≥ 3π (cos 2θ (1−η)² |z|²)^{1/3}`; and `∫_0^∞ u^9 e^{-2πηu} du = 9!/(2πη)^{10}`.
(c) By (a), (b) and absolute convergence, `∫_0^∞ Φ(y)y^{s−1}dy = γ_S(s)L_S(s)` for `Re s > 1`. By Mellin inversion,
`Φ(y) = (1/2πi)∫_{(2)} Λ_S(s) y^{-s} ds`. The line may be moved to `Re s = −1`: `|L_S(2+iT)| ≤ ζ(2)³`; on `Re s = −1`,
`|L_S(−1+iT)| = |γ_S(2−iT)/γ_S(−1+iT)|·|L_S(2−iT)| ≪ (|T|+2)^{9/2}` by the functional equation (F3) and Stirling
(`|Γ(x+iy)| ≍_x |y|^{x−1/2} e^{−π|y|/2}`); `L_S` is entire of finite order (F3), so by Phragmén–Lindelöf
`L_S(σ+iT) ≪ (|T|+2)^{9/2}` uniformly for `−1 ≤ σ ≤ 2`; and by Stirling `γ_S(σ+iT) ≪ (|T|+2)^{A} e^{−3π|T|/4}` there
(`γ_S` has no pole in `Re s > −1`, and on `Re s = −1` the simple pole of `γ_S` at `s = −1` is cancelled by the zero
of `L_S`; the estimate is needed only for `|T| ≥ 1`). So `Λ_S(σ+iT) → 0` exponentially, uniformly in the strip,
the horizontal segments vanish, and `Φ(y) = (1/2πi)∫_{(−1)} Λ_S(s)y^{-s}ds`; substituting `s ↦ 1−s` and using (F3)
again gives `Φ(y) = y^{-1}Φ(1/y)` for `y > 0`, hence in the sector by analytic continuation. Then `Φ(y)` decays faster than any power as `y → 0` as well, the Mellin
integral converges for all s and represents the entire function `Λ_S`. ∎

**Lemma V1 (trapezoidal rule = Poisson summation).** For `h > 0` and every `s ∈ ℂ`,

    h Σ_{k∈ℤ} Φ(e^{kh}) e^{skh} = Σ_{m∈ℤ} Λ_S(s + 2πim/h),

and by `Φ(e^{−kh}) = e^{kh} Φ(e^{kh})` the left side is `h[Φ(1) + Σ_{k≥1} Φ(e^{kh})(e^{skh} + e^{(1−s)kh})]`.

*Proof.* `f(v) = Φ(e^v)e^{sv}` extends analytically to `|Im v| < π/4` with decay faster than any exponential,
uniformly in `|Im v| ≤ π/8` (Lemma V0(b), (c)); so f is a Schwartz function and Poisson summation applies:
`h Σ_k f(kh) = Σ_m f̂(m/h)`, `f̂(ξ) = ∫ Φ(e^v) e^{(s − 2πiξ)v} dv = Λ_S(s − 2πiξ)`. ∎

**Lemma V2 (decay of Λ_S in the strip, elementary).** For `0 < θ < π/4`, `−1 ≤ σ ≤ 2`, `T ∈ ℝ`:

    |Λ_S(σ + iT)| ≤ M_θ e^{−θ|T|},     M_θ := 2 ∫_1^∞ y · sup_± |Φ(y e^{±iθ})| dy.

*Proof.* Let `T > 0`. Rotate the ray in `Λ_S(s) = ∫_0^∞ Φ(z) z^{s−1} dz` to `arg z = θ` (Cauchy; the arcs at 0 and
∞ vanish by the decay of Φ at both ends), split at `|z| = 1`, and apply `Φ(z) = z^{-1}Φ(1/z)` on the inner part:

    Λ_S(s) = e^{iθs} ∫_1^∞ Φ(ye^{iθ}) y^{s−1} dy + e^{iθ(s−1)} ∫_1^∞ Φ(ye^{−iθ}) y^{−s} dy.

`|e^{iθs}| = |e^{iθ(s−1)}| = e^{−θT}`, and `y^{σ−1}, y^{−σ} ≤ y` for `y ≥ 1`, `−1 ≤ σ ≤ 2`. For `T < 0` use `−θ`. ∎

With `θ = π/8`, `η = 1/5`: `b' = b_{η,θ} = 7.236 > 9/2`, so `n ↦ n³ e^{−b'(ny)^{2/3}}` decreases for `n ≥ 1`, `y ≥ 1`,
and `Σ_{n≥1} n³ e^{−b'(ny)^{2/3}} ≤ e^{−β} + (3/2)Γ(6, β)/β^6`, `β = b'y^{2/3}`; with `c_0 = Γ(6,b')e^{b'}/b'^6`,
`|Φ(ye^{±iθ})| ≤ C_η y (1 + (3/2)c_0) e^{−b'y^{2/3}}` and `M_θ ≤ 2C_η(1 + (3/2)c_0)(3/2)Γ(9/2, b')/b'^{9/2} ≤ 109.74`.

**Evaluation.** `h = 2π/640`, `k ≤ k_max = 517` (`e^{k_max h} > 160`). For `−1 ≤ Re s ≤ 2`, `|Im s| ≤ 73`:

    Λ_S(s) = h [ Φ_0 + Σ_{k=1}^{517} Φ_k (e^{skh} + e^{(1−s)kh}) ] + E,     |E| ≤ E_tail + E_alias ≤ 3.2·10^{-100} + 4.4·10^{-95},

where the `Φ_k` are balls (radius at most `2.91·10^{-84}`, set by the tolerance of the integrator), so that the
computed enclosure of `Z(t)` has radius about `10^{-83}` for `0 ≤ t ≤ 72`, against `|Z(60)| ≈ 9·10^{-50}` and
`|Z(72)| ≈ 2·10^{-61}`:

* `Φ_k ∋ Φ(e^{kh})` (`sym2_phi.py`): for `n ≤ max(1, ⌊160/y⌋)`, `K(ny)` by Arb's rigorous integrator on an exact
  interval `[u_lo, u_hi]`, plus `∫_0^{u_lo} ≤ u_lo^{10} e^{−πx²/u_lo²}` and `∫_{u_hi}^∞ ≤ Γ(10, 2πu_hi)/(2π)^{10}`; the terms
  `n > N` are bounded by `C_η y [(N+1)³ e^{−β(N+1)^{2/3}} + (3/2)Γ(6, β(N+1)^{2/3})/β^6]` (`η = 1/10`, `θ = 0`, `β = b_η y^{2/3}`;
  first term plus integral, the summand `n³e^{−βn^{2/3}}` being decreasing);
* `E_tail`: `Σ_{k>k_max} h|Φ(e^{kh})|·2e^{2kh} ≤ 2C_η(1+(3/2)c_0)(3/2)Γ(9/2, b_η y_0^{2/3})/b_η^{9/2}`, `y_0 = e^{k_max h}`;
* `E_alias = Σ_{m≠0} |Λ_S(s + 640im)| ≤ 2M_θ e^{73θ} e^{−640θ}/(1 − e^{−640θ})` (Lemmas V1, V2).

On the critical line `Z(t) := Λ_S(1/2 + it) = h[Φ_0 + Σ_k 2Φ_k e^{kh/2} cos(tkh)] + E` is real.
Consistency: `L_S(2)` from this formula is `0.80587520944866875961…`, inside the Euler-product enclosure
`0.8059 ± 3.8·10^{-5}` (asserted by the program) and equal to PARI's value to the 30 digits printed;
`L_S(1/2) = 0.505549375222713515001575353290…`, `L_S(1) = 0.631792945727883203011116326435…` (PARI: same digits).

### 5.2 Zeros on the line

`sym2_verify.py` evaluates the sign of `Z` (certain) on the grid `t = j/16`, `0 ≤ j ≤ 1152`, and finds 55 sign
changes in `(0, 72]`, 41 of them in `(0, 60]`; each is bisected to an interval `[l_j, r_j]` of length `2^{-38}`
(all 55; `max_interval_width = 3.638·10^{-12}` in the certificate) with certified opposite signs at the ends. So there are zeros `1/2 + iγ_j`, `γ_j ∈ [l_j, r_j]`, `j = 1..55`.
(They agree with PARI's `lfunzeros` to all digits compared.)

### 5.3 Completeness without a contour integral

**Lemma V3.** Let `γ_1 < … < γ_J` be ordinates of zeros of `L_S` on the critical line (one zero counted for
each j), `w(u) = (3/2)/(9/4 + u²)`, and for real T

    Rem(T) := Re (Λ_S'/Λ_S)(2 + iT) − Σ_j [w(T − γ_j) + w(T + γ_j)].

Then (i) `Rem(T) ≥ 0` and `|Rem'(T)| ≤ Rem(T)`; (ii) if there is a zero `ρ_0 = β_0 + iγ_0` of `L_S` in the critical
strip which is not among the `1/2 ± iγ_j` (as a multiset: a multiple zero at some `1/2 ± iγ_j` counts), then
`Rem(γ_0) ≥ 2/3`. Hence: if `Rem(T) < 2/3` for all `T ∈ [0, H]`, the zeros of `L_S` with `0 ≤ Im ρ ≤ H` are exactly the
`1/2 + iγ_j` with `γ_j ≤ H`, each simple, and there is no zero with `Im ρ = 0`.

*Proof.* By (F4), `Re (Λ_S'/Λ_S)(2+iT) = Σ_ρ f_ρ(T)`, `f_ρ(T) = a/(a² + (T−γ)²)`, `a = 2 − β ∈ (1, 2)`; for `ρ = 1/2 ± iγ_j`,
`f_ρ = w(T ∓ γ_j)`. So `Rem` is the sum of `f_ρ` over the remaining zeros: (i) `Rem ≥ 0`, and
`|f_ρ'| = f_ρ · 2|u|/(a²+u²) ≤ f_ρ/a ≤ f_ρ`. (ii) If `β_0 = 1/2`: `f_{ρ_0}(γ_0) = 1/(3/2) = 2/3`. If `β_0 ≠ 1/2`: `1 − β_0 + iγ_0` is
another remaining zero (functional equation and reality), and `1/(2−β_0) + 1/(1+β_0) ≥ 4/3`. ∎

The program computes, for `T = j/8`, `0 ≤ j ≤ 480`, an upper bound for `Rem(T)`: the Euler product over
`p ≤ 2·10^5` with tail `≤ 3.03(log P+1)/P`, `γ_S'/γ_S` in Arb, and a lower bound for the found zeros' terms
(all 55, with `γ_j` replaced by the interval `[l_j, r_j]`). Result: `max_j Rem(j/8) ≤ 0.185062` (at `T = 60`), so by (i)
`Rem(T) ≤ 0.185062·e^{1/16} ≤ 0.196998 < 2/3` on `[0, 60]`. By Lemma V3 and the symmetry `ρ ↦ ρ̄`, Theorem 4
follows. ∎

(`Rem` is small because the zeros up to 72 are known; at `T = 0` it is `0.0686`, the contribution of all zeros
above 72. Mutation test: deleting any one of the zeros #1, #8, #21, #41 from the list raises `max Rem` to
`0.743, 0.741, 0.758, 0.842 > 2/3`: the test rejects, `selftest.py` M3.)

---------------------------------------------------------------------------

## 6. Proof of Theorem 3

(a)⟹(b), (a)⟹(d): Theorems 2, 1. (b)⟹(c), (d)⟹(e): trivial. It remains to prove (c)⟹(a) and (e)⟹(a).

(c)⟹(a) (Landau). Let `A(x) = xR(x) = Σ_{p≤x} λ(p)(x − p)`; `|A(x)| ≤ 2x²`. For `Re s > 1`,

    F(s) := ∫_1^∞ A(x) x^{-s-2} dx = P(s)/(s(s+1)),      P(s) = Σ_p λ(p) p^{-s}.

From the Euler products, `log L(s,Δ) = Σ_p Σ_k c_k(p)/(k p^{ks})` and `log L_S(w) − log ζ(w) = Σ_p Σ_j c_{2j}(p)/(j p^{jw})`, so

    P(s) = log L(s, Δ) − (1/2) log L_S(2s) + (1/2) log ζ(2s) + H(s),

`H(s) = (1/2) Σ_p Σ_{j≥2} c_{2j}(p)/(j p^{2js}) − Σ_p Σ_{k≥3} c_k(p)/(k p^{ks})` (bounded coefficients), absolutely convergent and
analytic in `Re s > 1/3`. On the real segment `σ > 1/2`: `ζ(2σ) > 0` and finite; `L_S(2σ) > 0` (Euler product,
`2σ > 1`); and `L(σ, Δ) ≠ 0`, because

    (2π)^{-w} Γ(w) L(w − 11/2, Δ) = ∫_0^∞ Δ(iy) y^{w−1} dy = ∫_1^∞ Δ(iy)(y^{w−1} + y^{11−w}) dy > 0      (w real),

`Δ(iy) = e^{-2πy} Π(1 − e^{-2πny})^{24} > 0`. Hence `F` continues analytically to a neighbourhood of every real
`σ > 1/2`. If `A(x) ≥ 0` for `x ≥ x_0`, Landau's theorem (a Mellin transform of an eventually one-signed
function has a singularity at the real point of its abscissa of convergence) shows that the integral
converges and `F` is analytic in `Re s > 1/2`. So `P`, hence `log L(s, Δ) = P + (1/2)log L_S(2s) − (1/2)log ζ(2s) − H`, is
analytic in `Re s > 1/2` (there `L_S(2s) ≠ 0`, `ζ(2s) ≠ 0, ∞`): `L(s,Δ) ≠ 0` for `Re s > 1/2`, which is GRH(Δ) by
the functional equation.

(e)⟹(a): the same with `A` replaced by `Θ_1`: `∫_1^∞ Θ_1(t) t^{-s-2} dt = P_1(s)/(s(s+1))`,
`P_1(s) = Σ_p λ(p) log p · p^{-s} = −(L'/L)(s,Δ) + (L_S'/L_S)(2s) − (ζ'/ζ)(2s) + H_1(s)` (Lemma 1), which is meromorphic in
`Re s > 1/3` with poles in `Re s > 1/2` exactly at the zeros of `L(s,Δ)` there (simple poles with residue
`−m_ρ ≠ 0`), and no pole on the real segment `σ > 1/2`. Landau's theorem gives analyticity in `Re s > 1/2`, so no
such zeros. ∎

---------------------------------------------------------------------------

## 7. Certificates, programs, reproducibility

All programs in `code/`; outputs in `data/`. Python 3.13.18, python-flint 0.9.0 (Arb), PARI/GP (version printed at the top of `data/independent_check.txt`) for the
independent checks. Run times on 4 cores.

| program | what | output (sha256, first 16 hex) | time |
|---|---|---|---|
| `tau_table.py` | exact `τ(n)`, `n ≤ N`, from `η^{24} = (Σ(−1)^k(2k+1)q^{k(k+1)/2})^8` (FLINT `fmpz_poly`); checks: first 12 values, `τ ≡ σ_11 mod 691`, multiplicativity, `τ(p²) = τ(p)² − p^{11}` | — | 30 s (`N = 10^7`) |
| `cert_constants.py 10000000 10000000` | `σ_0²`, `a_0`, `a_{−1}`, `L(1,Δ)`; Euler products at 2 (`p ≤ 10^7`); `C_S`, `K_S`, `a_{S,0} > 0`; `c_ψ`, `θ(y) < y` | `constants.json` `815cedbfa8a8ced7` | 36 s |
| `cert_small.py 10000000` | checks (A), (B); gap bounds for `G`, `N`; `A_1(T)`, `A_2(T)`; `t*`, `x*` | `small_range_10000000.json` `2b4e1f9ce309cf54` | 31 s |
| `bounds.py` | `L_A`, `L_B`, `T(H)` in Arb; asserts its rational constants dominate `constants.json` | `bounds.txt` | 5 s |
| `sym2_phi.py 640 160 4` | table `Φ_k`, `0 ≤ k ≤ 517`, 330 bits | `sym2_phi_M640.json` `ada15ff1be6ce72b` | 13 s |
| `sym2_verify.py 640 60 72 200000` | Theorem 4 | `sym2_rh_H60.json` `1e6a096b3c35965b`, `sym2_verify_H60.txt` | 46 s |
| `selftest.py` | functional equation of Φ checked directly; mutation tests M1–M4 | `selftest.json`, `selftest.txt` | 2 min |
| `independent_check.gp` | PARI: all constants, residue identities at `t = 1`, direct-versus-formula, zeros, `A_1`, `A_2` | `independent_check.txt` | 3 min |
| `explore_float.py` | long-double exploration (not used in proofs) | — | 1 min |

How `σ_0²` is certified (`cert_constants.py`, C1). With `w = s + 11/2`:
`Λ*(w) := (2π)^{-w}Γ(w)L(w − 11/2, Δ) = Σ_n τ(n)[E(w,n) + E(12−w,n)]`, `E(w,n) = ∫_1^∞ e^{-2πny}y^{w−1}dy = (2πn)^{-w}Γ(w, 2πn)`
(from `Δ(i/y) = y^{12}Δ(iy)`), and `Λ*'(w) = Σ_n τ(n) ∫_1^∞ e^{-2πny} log y (y^{w−1} − y^{11−w}) dy`. At `w = 13/2`: 40 terms;
incomplete gamma in Arb; the derivative integrals by Arb's rigorous integrator on `[1, Y_n]`, `Y_n` an exact
integer, plus the tail bound `e^{-aY}Y^{6.5}/(a − 6.5/Y)`, `a = 2πn`; the terms `n > 40` are bounded by a geometric
series using `|τ(n)| ≤ n^{6.5}`. `σ_0² = 2Λ*'/Λ*(13/2)`. Second implementation: PARI's `lfun` (different
algorithm) gives `0.11225214047190852674558342625324625632`, inside the enclosure.

Independent implementations, by leg:

| leg | certified implementation | second implementation | agreement |
|---|---|---|---|
| `τ(p)` | `η^{24}` with FLINT | PARI `mfcoefs` (`n ≤ 10^5`) and `ramanujantau` (Hurwitz class numbers, `p ≤ 2·10^4`) | identical |
| `σ_0², a_0, a_{−1}, L(1)` | theta integral in Arb | PARI `lfun` | 30 digits |
| `(L_S'/L_S)(2)`, `(Λ_S'/Λ_S)(2)` | Euler product, `p ≤ 10^7` | PARI `lfun` of a hand-built `lfuncreate` | inside enclosure |
| `Θ, Θ_1, D_0, D_1, A_1, A_2` at 97, 997, 99991 | Arb | PARI floats | all digits printed |
| zeros of `L_S` up to 60 | Poisson evaluation + sign changes | PARI `lfunzeros` | 41 zeros, all digits |
| functional equation of `Λ_S` | used via `Φ(1/y) = yΦ(y)` | `Φ` from its series at `y` and `1/y`, `y = 9/10, 3/4, 1/2` | `≤ 1.2·10^{-83}` |
| explicit formulas (residues) | proof | identities at `t = 1` and direct sums, PARI | `10^{-38}` resp. within truncation |
| `T(60)` | Lemma 5 | closed form + 41 zeros: true tail `≈ 0.104 ≤ 0.1846` | consistent |

Mutation tests (`selftest.py`; each mutant must be rejected, all are): M1 sign-flipped λ, and
`τ(3), τ(5), τ(11)` negated → check (A)/(B) fails; M2 `σ_0² := 0.7` → no positivity threshold, and `H = 12` →
`L_B` never positive; M3 one zero deleted → `Rem ≥ 2/3`; M4 `τ(2) := τ(2)+1` → `σ_0²` leaves the pinned enclosure.

---------------------------------------------------------------------------

## 8. Relation to the literature (details and quotations in `LITERATURE.md`)

* Koyama–Kurokawa (Proc. Japan Acad. A 98 (2022)): under the Deep Riemann Hypothesis,
  `Σ_{p≤x} τ(p)p^{-6} = (1/2) log log x + c + o(1)`. Sheth (Math. Proc. Camb. Phil. Soc. 179 (2025), Cor. B):
  the same under GRH alone outside a set of finite logarithmic measure; hence positivity on a set of
  natural density 1 under GRH alone is **known** (with Hayani's remark arXiv:2512.23302 (1.5)), and is not
  claimed here. Koyama–Sheth (arXiv:2509.04187): all newforms, under DRH.
* Sarnak's letter to Mazur (2007): under GRH, `(log x/√x) Σ_{p≤x} τ(p)p^{-11/2}` has a limiting logarithmic
  distribution with mean 1, the mean coming from the prime squares / `Sym²`. Theorems 1–2 are the all-x,
  explicit, smoothed counterpart; Theorem 3 is the converse.
* Hayani (arXiv:2512.23302), Theorem 1.4: Cesàro means of the mod-4 race at every x under GRH with a
  non-explicit constant (Dirichlet characters only). Suzuki (Ramanujan J. 68 (2025); IJNT 21 (2025)):
  equivalences "GRH ⇔ eventual sign" for Dirichlet characters, and for the Selberg class with a
  different (screw-function) sum with archimedean terms; τ is not treated explicitly.
* G. Hiary, S. Ireland, M. Kyi, "A method for verifying the generalized Riemann hypothesis"
  (arXiv:2408.00187, 2024): they evaluate the Hadamard sum `Σ_ρ 1/(ρ−δ)` at one real point `δ ≤ 0` from `L'/L` in the
  region of absolute convergence and use that every zero contributes nonnegatively, so that the zeros found
  on the line by sign changes plus the least possible contribution of a missed or off-line zero below a
  height η would exceed the total: RH and completeness up to η without the argument principle or Turing's
  method. Their examples are ζ, `L(s, χ_d)`, `L(s, Δ)` (η = 84, 46 zeros) and an elliptic curve; no symmetric
  square. **Lemma V3 is a variant of this method** (the positive partial-fraction sum at the points `2 + iT`,
  `0 ≤ T ≤ H`, with a Lipschitz bound, instead of one real point); it was written before that paper was read
  and is not claimed as a new method. What is claimed is the result, Theorem 4, for `L(s, Sym²Δ)`.
* Not found in the literature searched: any statement for τ valid for every x beyond an explicit `x_0`;
  the equivalence of Theorem 3; a rigorous verification of RH for `L(s, Sym²Δ)` to any height (Theorem 4).
  "Not found" is a statement about the searches recorded in `LITERATURE.md`, not a proof of absence.
* Method: `../explicit/EXPLICIT.md` (this repository) for `χ_4, χ_3`. New ingredients needed for τ:
  Lemma 1 (all even prime powers through `ψ_S − ψ`), the treatment of `Sym²` without GRH (Lemmas 4, 5 and
  Theorem 4), the closed form `σ_0² = 2Λ'/Λ(1)` certified from the theta integral, and the handling of the
  initial negative patch in Lemma 7.

---------------------------------------------------------------------------

## 9. Status and trust

* Theorems 1(b): proof on paper + certified constants + certified finite range. Conditional on GRH(Δ),
  GRH(Sym²Δ).
* Theorems 1(a), 2, 3: as above, and they use Theorem 4. Conditional on GRH(Δ) only.
* Theorem 4: computer-assisted (Arb), no hypothesis. Trusted: FLINT/Arb (incomplete gamma, rigorous
  integration, elementary functions), the programs `sym2_phi.py`, `sym2_verify.py`, and facts (F1), (F3), (F4).
* Cited inputs: (F1)–(F6). (F6) can be weakened to Rosser–Schoenfeld (Remark 4.3).
* Review: see `REFEREE-REPORT.md` (in-house adversarial referee, an AI agent) and §10.
* Not done: Lean formalisation; `e^{-p/x}` and other weights; the weights and density statements T3, T4 of
  the original work order (dropped by instruction); a larger H in Theorem 4 (H = 100 would give
  `T_100 ≤ 0.1215` and `lim L_B ≥ 0.4324`; only the table parameters need to change).

## 10. Review history

**Referee pass (2026-11-26, AI agent instructed to find the error; report in `REFEREE-REPORT.md`).** Verdict:
"repairable gaps; no fatal error found"; Lemmas 0–7 (Lemma 0 was added afterwards, see M9), the assembly of §3,
Lemmas V0–V3 and §6 were re-derived by the referee, all five certificates reproduced byte for byte, and
identity (2.1), the explicit formula for `J_S`, `T(60, 51.45)`, `Rem(60)` and the zero counts were
re-checked with independent PARI computations. Findings and what was done:

* R1 (the gap bounds for `G`, `N` did not cover `[100,101)` and `[10^4, 10007)`): `cert_small.py` now takes all gaps
  `[q,p]` with `p > a`; rerun; the minima are unchanged (the certificate file is byte-identical), and the text of
  §1 and §4 says which gaps are used.
* R2 (the quotation of Shimura's gamma factor was wrong by a factor `2^{-s}`): the quotation is removed; (F3)
  now states `γ_S` in the form used, says that the originals were not consulted, and points to the two
  numerical confirmations of the functional equation.
* M1 (digits of the first zero), M3 (last digits of `ε_1`), M4 (wrong section number), M7 (`1.112·10^5`),
  M8 (wording of the remark after Theorem 4), M6 (Python version; `independent_check.txt` regenerated with
  the PARI version line): corrected.
* M2 (table radii up to `8·10^{-64}` because of a loose n-tail at `N = 1`; ten isolating intervals above height 60
  wider than stated): the n-tail is now "first term + integral from `N+1`"; table radii `≤ 2.91·10^{-84}`; all 55
  intervals have length `2^{-38}`; §5.1 states the radii. `sym2_phi_M640.json` and `sym2_rh_H60.json` changed
  (new hashes in §7); the verdict and `max Rem` are unchanged.
* M5 (float minimum of `G` not in a data file): `data/explore_float_1e7.txt` added.
* M9 (contour shift only cited): written out as Lemma 0.
* M10 (stale docstrings, dead variable): removed.

The referee did not consult the originals of Deligne, Hecke, Shimura, Gelbart–Jacquet, Rosser–Schoenfeld,
quoted Iwaniec–Kowalski from memory, confirmed (F6) against the arXiv text of BKLNW, and did not write a
second rigorous implementation of the `Φ` table or of the completeness test. The changes made after the
report (Lemma 0, the new n-tail, the gap rule) were not re-refereed.

**Independent numerical audit (2026-11-26, a second AI agent; own implementations first, repository code read
afterwards; report in `AUDIT-NUMERICS.md`).** Floating point only (PARI and mpmath, exact `τ(p)`); six items,
all PASS, no digit of any certified enclosure contradicted: (1) `σ_0²`, `a_0`, `a_{−1}` by PARI `lfun` and by
mpmath quadrature of the η-product, and the zero sum with 3035 zeros to height 2000 plus tail, which matches
`σ_0²` to `8·10^{-11}`; (2) the `Sym²` constants by PARI and by a Mellin–Barnes computation in mpmath; (3) 41 zeros
in `(0, 60]`, 55 in `(0, 72]`, all inside the stored isolating intervals, and — by a different method from §5.3 —
an argument-principle count giving 41 at `T = 60` and 55 at `T = 72`; `Rem(0) = 0.06835`, `Rem(60) = 0.18486`;
(4) positivity of `Θ_1`, `R` at every prime up to `9 999 991`, `t*`, `x*`, `A_1(97)`, `A_2(97)`, the float minimum of `G`;
(5) the values of `L_B`, its monotonicity on `[527, 10^14]` (scan at `2·10^6` points; in this note it is a paper
argument), `T(60, 51.45) = 0.18460434`, true tail above 60 between `0.10403` and `0.10409`; (6) both explicit
formulas, with discrepancies of the size of the truncation and `10^{-57}` at `t = 1`. Its two remaining errata
(the size of `|Z|` in §5.1 and the last digit of `max Rem` after the table was regenerated) are corrected above.
The audit notes that the bounds `G ≥ 0.5814` on `[100, N]` and `N ≥ 0.1718` on `[11, N]` are valid but not sharp
(true infima `0.6117`, `0.1747`), and that `K_S` bounds a quantity whose true value is `−0.228`.

**Second referee (2026-11-26, commissioned by the parent session; independent agent, own code; report
`REFEREE-2-TAU.md` on the parent branch).** No fatal error and no wrong step in Theorems 1–4; each rated
"complete modulo listed gaps"; the three changes made after the first report (Lemma 0, the n-tail of the `Φ`
table, the gap rule) were checked and found correct; `τ(n)` for `n ≤ 10^7` was recomputed from `(E_4³ − E_6²)/1728`
(agreement everywhere), as were the small range, `t*`, `x*`, `σ_0²`, the `Sym²` constants and the 41 zeros. Repairs
made in response: G1 — (F3) now quotes two sources that were actually read (Zagier 1977, eq. (4) and
Corollary 1; Balkanova–Frolenkov 2018), with the translation to `γ_S`; what remains unread (Shimura,
Gelbart–Jacquet, Rankin; the finite order of `L_S`) is said there. G2 — the decay needed for the contour shift
in Lemma V0(c) is proved (functional equation, Stirling, Phragmén–Lindelöf). G3 — (F6) quotes BKLNW
Corollary 11.2 and Dusart's preprint from the texts; Remark 4.3 is restated with computed constants (it had
the roles of Theorems 1 and 2 reversed). Minor — the "t = 1" checks are labelled as bookkeeping tautologies;
the Euler-product tail is proved; `RESULTS.md` no longer says "T2 is known" ambiguously; Hiary–Ireland–Kyi is
cited, and Lemma V3 is described as a variant of their method. These repairs were not refereed again.

