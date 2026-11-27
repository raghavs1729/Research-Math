# Line (1), the dual side: what the explicit formula says about relations among the zeros

AI-produced (Anthropic's Claude), 2026-11-27. Status: **assessment of a route**, with one
standard lemma. Nothing here concerns Problem 12 of arXiv:2407.03530. The main statement is
not proved here.

Notation:
- `Γ`, `m_γ`, `G = ⟨Γ⟩` and the hull `A` are as in `monotonicity/density/DENSITY.md` §1;
- `Z := Σ_{γ∈±Γ} m_γ δ_γ` is the zero measure.

By the explicit formula (`DENSITY.md` (I5)), `Ẑ` on `ℝ ∖ {0}` is
`−Σ_n Λ(n)χ(n)n^{−1/2}(δ_{log n} + δ_{−log n})` plus a smooth density. The atoms sit on
`log ℕ`, whose `ℚ`-span has the `ℚ`-independent basis `{log p}`. The question (parent session,
2026-11-27) is whether this independence on the dual side can be transported into a
constraint on the `ℚ`-relations among the zeros, i.e. on the hull.

## 1. What the explicit formula excludes

- **Cyclic `G`** (`DENSITY.md` Remark 7.6). If all ordinates lie in `(2π/P)ℤ`, then `Ẑ` is
  `P`-periodic. So the atom `−χ(p)(log p)/√p` at `log p` would recur at `log p + kP`, whereas
  the atoms at `log n` have masses `Λ(n)/√n → 0`. This excludes **every** cyclic spectrum, in
  particular the spectra `Ωℕ` of all known integer counterexamples (`DENSITY.md` Ex. 9.2,
  Thm 9.6). They are maximally resonant on the dual side.
- **Rank one with slowly growing denominators.** Suppose `G ⊂ αℚ`. Let `D(T)` be the lcm of
  the denominators of `γ/α` over `γ ∈ Γ_T`, and put `P := 2πD(T)/α`. Then `e^{iγP} = 1` for all
  `γ ≤ T`.
  - Apply the Weil explicit formula at `y = 0` and at `y = P`, with an even test function
    whose Fourier transform is supported in `[−T,T]` and equal to `1` on `[−T/2,T/2]`.
  - The zero sides agree. At `y = 0` that side is `≍ T log T`: all phases are `1`, and this is
    the density term of the formula.
  - The prime side at `y = P` is a sum over `n` in a window of relative width `≍ 1/T` around
    `e^P`. Trivially it is `≪ P e^{P/2} + T P e^{−P/2}`.
  - Hence `P ≥ 2 log T − 3 log log T` for large `T`. If `P` stays bounded, the zeros up to `T`
    would lie in an AP of fixed spacing, which is impossible.

  So `D(T) ≥ (α/π − o(1)) log T`. Counting distinct ordinates, `#Γ_T ≤ T D(T)/α + 1`, already
  gives `D(T) ≫ #Γ_T/T`. The explicit formula gives the logarithmic bound without any
  information on multiplicities. Equivalently,
  by the Landau–Gonek formula, the numerators of `γD′/α` must equidistribute mod `1` for
  `D′ ≪ log T`. (This is the "`P_T ≳ log T`" of `PROGRESS.md` §2, written out. It is a
  standard argument, and we claim nothing new for it.)

## 2. Why the naive transport stops there

**Each relation among the zeros, and each recurrence of the hull orbit, becomes a coincidence
on the prime side, but only at specific, astronomically large points.**
- Suppose the hull orbit `ι(y) = (e^{iγy})` returns `δ`-close to a point `θ_0` of `A` (or to
  the identity) at time `y_0`. Then the smoothed prime sums
  `S_φ(y) = Σ_n Λ(n)χ(n)n^{−1/2} φ(y − log n)` satisfy `S_φ(y + y_0) ≈ S_φ(y)`, or the relation
  dictated by `θ_0`.
- In rank `1`, `y_0 = 2πD/α`. In finite rank `d`, `y_0 ≲ δ^{−d}` by Dirichlet. After an atom,
  `θ_0` is a point of the finite group `H_∞^⊥` of `HYP.md` Corollary 1.4.
- The requirement on the primes concerns windows of relative width `≍ δ` at `x ≈ e^{y_0}`. Its
  size is within what GRH allows: `|S_φ| ≤ Σ_γ m_γ|φ̂(γ)|`.
- No unconditional information on primes in short windows at **prescribed** points of that
  size is available. Results for almost all short intervals (Selberg, Huxley) say nothing at
  prescribed points.

**Statistical information on zeros cannot see exact relations.**
- Pair correlation (Montgomery), Landau–Gonek and moments of `S_φ` give asymptotics for
  averaged quantities.
- A `ℚ`-relation is a measure-zero event. Perturbing the ordinates slightly to rationals
  changes none of these statistics (`PROGRESS.md` §2).

**Relations among the zeros are multiplicative relations among Bohr points.**
- For a zero `t_j`, put `z^{(j)} := (p^{−it_j})_p ∈ 𝕋^{primes}`.
- Because `{log p}` is `ℚ`-independent, `Σ_j n_j t_j = 0` holds iff `Π_j (z^{(j)})^{n_j} = 1`
  coordinatewise. (`(Σ n_j t_j) log p ∈ 2πℤ` for `p = 2, 3` forces `Σ n_j t_j = 0`.)
- So LI is a statement about **unlikely intersections**. The Kronecker line meets the "zero
  set of `L`" on the Bohr torus, and the intersection points should satisfy no multiplicative
  relations.
- Zilber–Pink-type finiteness is known only for algebraic subvarieties of finite-dimensional
  tori (Bombieri–Masser–Zannier; **MEMORY**). Here the Euler product does not converge on the
  critical line, the torus is infinite-dimensional, and the points lie on a transcendental
  one-parameter subgroup. Even for finite Dirichlet polynomials, independence of real zeros
  would need Schanuel-type transcendence input. No unconditional route is visible.

**The crystalline-measure / Fourier-quasicrystal theory does not apply.**
- Kurasov–Sarnak, *Stable polynomials and crystalline measures*, arXiv:2004.05678;
  Olevskii–Ulanovskii, *Fourier quasicrystals with unit masses*, arXiv:2009.12810; Favorov,
  arXiv:2311.02728 (zeros of absolutely convergent Dirichlet series with bounded spectrum).
- These classify discrete measures whose Fourier transform is also discrete, under
  growth/density conditions: zero sets of exponential polynomials, linear density.
- The zero measure `Z` has density `(1/2π) log T`, and `Ẑ` has an absolutely continuous part.
  Kurasov–Sarnak remark that the zeta explicit formula alone does not give a Fourier
  quasicrystal.
- So no structure theorem from that theory constrains `Γ` (web check 2026-11-27, two
  searches; abstracts and secondary pages only).

## 3. Verdict

**Killed in its naive form.**
- The dual side, as encoded by the explicit formula, excludes cyclic `G`. In rank one it
  forces `D(T) ≳ log T`. Under GRH, every finer constraint on the relations becomes a
  statement about primes at prescribed large points, where nothing unconditional is known.
- An atom forces `⟨Γ_{>T}⟩` to be eventually constant, of index `≤ 1/α` (`HYP.md`
  Corollary 1.4). It does not force a resonant (low-rank, AP-like) spectrum. So the
  "maximal resonance" of the known counterexamples is not a feature that an atom provably
  shares.

## 4. Joining the structure side and the arithmetic side: the exact gap

The parent session (2026-11-27) proposed a closing strategy:
- (A) an atom forces low effective rank of the truncated zero sets (few independent phases);
- (B) the explicit formula excludes low-rank (resonant) configurations, robustly;
- if the bounds overlap, the problem is solved.

Here are both halves, in their strongest rigorous forms available here, and the gap between
them.

**(A) Proposition 4.1 (an atom bounds the free sets; GRH + NCZ).** Suppose `ν_M({c}) = α > 0`.
Let `ε ∈ (0,1)`, and let `T` be so large that `V_M^{>T} ≤ αε²/4`. Since
`V_M^{>T} ≪_M (log T)²/T` (multiplicities `m_γ ≪ log γ`), some
`T ≍_M α^{−1}ε^{−2}log²(1/(αε))` suffices. Let `Λ ⊂ Γ_T` be **independent modulo the rest**: the
images of `Λ` in `⟨Γ_T⟩_ℚ/⟨Γ_T∖Λ⟩_ℚ` are linearly independent. Then

    #{λ ∈ Λ : m_λ|r_λ(M)| ≥ 8ε} ≤ C/α²,     C absolute.

*Proof.* By `DENSITY.md` Lemma 9.13 (predecessor branch),
`λ_{A_T}{|F_{M,T} − c| < 2ε} ≥ α − V^{>T}/ε² ≥ 3α/4`.
- Put `H := ⟨Γ_T∖Λ⟩`. Then `G_T/H` is free with basis the images of `Λ`: a `ℤ`-dependence
  would be a `ℚ`-dependence.
- So `K := H^⊥ ⊂ A_T` is a torus, and `k ↦ (e_λ(k))_{λ∈Λ}` is an isomorphism `K → 𝕋^Λ`.
- For `θ ∈ A_T` and `k ∈ K`,
  `F_{M,T}(θk) = c_θ + Σ_{λ∈Λ} a_λ cos(φ_λ(k) + β_λ(θ))`, with `a_λ := m_λ|r_λ|` and the
  `φ_λ` i.i.d. uniform under Haar measure on `K`.
- Fubini (invariance of `λ_{A_T}` under `K`) bounds `3α/4` by the concentration function at
  scale `4ε` of a sum of independent terms `a_λcos φ_λ`.
- For `a_λ ≥ 8ε`, each term has concentration `≤ 1/3` at that scale (arcsine law).
- The Kolmogorov–Rogozin inequality (Petrov, *Sums of independent random variables*,
  ch. III; **MEMORY**) gives `3α/4 ≤ C′(n·2/3)^{−1/2}`, where `n` is the number of such `λ`. ∎

This is the quantitative, multi-element form of the direct-sum argument (Devin, Thm 2.2;
`DENSITY.md` §3). What it constrains is **independence of the low zeros modulo the higher
ones**. All but `C/α²` of the zeros `γ ≲ 1/ε` must be `ℚ`-combinations of the other zeros up to
height `T ≍ α^{−1}ε^{−2}log²(1/(αε))`.

It does **not** bound the rank of `⟨Γ_{1/ε}⟩` itself. The low zeros may be independent among
themselves, as in Example 3.3 of `DENSITY.md`, where every element is dependent on higher ones.
It also gives no bound on the length or coefficients of the relations, and no
arithmetic-progression structure. Corollary 1.4 of `HYP.md` is the infinite-height,
`ℤ`-version of the same phenomenon. (Further structural consequences of an atom are in
`INVLO.md`: Corollary C′ and Proposition E.)

**(B) Lemma 4.2 (no early recurrence; GRH; standard, proof only sketched).** Fix an even
`Φ ∈ C_c^∞(−1,1)` with `0 ≤ Φ ≤ 1` and `Φ = 1` on `[−½,½]`. Put `ĥ(r) := Φ(r/U)` and
`N_ĥ := Σ_{γ>0} m_γĥ(γ) ≍ U log U`. Then for `K/U ≤ y ≤ 2log U + 2log log U − C`,

    |Σ_{γ>0} m_γ ĥ(γ) cos(γy)| ≤ η N_ĥ,

where `η ≤ η_1(K) + η_2(C) + o_{U→∞}(1)`. Here `η_1(K) → 0` as `K → ∞` and `η_2(C) → 0` as
`C → ∞`, and all three depend on `Φ`.

The lemma is stated for the two-sided sum `Σ_{γ∈±Γ} = 2Σ_{γ>0}m_γĥ(γ)cos(γy)`, which is all
the explicit formula gives; the zeros are symmetric. The sine part is numerically of the same
size (referee 4, `referee4/r4_line1.py`, at `U = 600`), but that is not proved here. The
cosine version is what "The gap" uses: a near-recurrence makes the cosine sum `≈ N_ĥ`.

*Proof sketch* (standard; Landau–Gonek type; **not written out in full**, so by the labels
of `DENSITY.md` this is a sketch):
- The Weil explicit formula makes the left side, times `2`, equal to `D(y) − P(y)`.
- The density term `D(y)` is `N_ĥ` times a Schwartz function of `Uy`, up to `O(log U)`, plus
  the contribution of the logarithmic weight near `ξ = 0`, which is `O(1/y)`. So
  `|D(y)| ≤ η_1(K)N_ĥ + O(U/K)` for `y ≥ K/U`.
- `P(y) = Σ_nχ(n)Λ(n)n^{−1/2}[g(log n − y) + g(−log n − y)]` with `g(u) = Uφ(Uu)`, where `φ`
  is the inverse transform of `Φ`. The second term is `O_A(U^{−A})`, by the far-part bound
  below (and by Schwartz decay for `n = 2` when `y` is tiny). The function `g` is Schwartz but, being of exponential type, does not
  decay exponentially. So the prime side converges only **conditionally**, under GRH. Split
  it at `|log n − y| = 1`.
  - *Near part:* Brun–Titchmarsh on unit intervals in `Uu`, against `sup|φ|`. This gives
    `≪ e^{y/2}·log x/log(x/U) + U y e^{−y/2}`, with `x = e^y`. For `y ≤ log U` the window holds
    at most two integers.
  - *Far part:* partial summation with `ψ(t,χ) ≪ √t log²t` (GRH). This gives
    `≪ U y²∫_{|s|>U}(|φ|+|φ'|)(1+|s|/U)²ds = O_A(U^{−A})`.

  The factor `log x/log(x/U)` matters only for `x ∈ [U, U^{1+δ}]`, and there it is harmless. ∎

In words: the zeros up to `U` are not biased, in weighted density, towards any arithmetic
progression of spacing `≥ π/((1+o(1)) log U)`, i.e. more than half the mean spacing. No
unconditional information on primes reaches further. For `y ≥ (2+δ) log U` the window holds
`≫ √x` primes, and a bias of the required size `≍ √x log U` is compatible with every sieve
bound.

**The gap.**
- **What (B) can use:** a weighted near-recurrence of the zeros up to `U` at a time
  `y_0 ≤ (2+o(1)) log U`. That is an approximate AP structure of `Γ_U` with spacing
  `s = 2π/y_0 ≥ π/((1+o(1)) log U)`.
- **What (A) can supply:** dependence of the low zeros on higher zeros, with no size control.
  Even granting the strongest conceivable structural output, bounded rank `d` of `⟨Γ_U⟩`: there
  are `#Γ_U ≫ U` distinct zeros (Riemann–von Mangoldt with `m_γ ≪ log γ`), and they force a GAP
  of size `≫ U`. (`≫ U log U` distinct zeros would follow from a bound
  `Σ_{γ≤U}m_γ² ≪ N(U)`, e.g. from a positive proportion of simple zeros. This was not checked for
  a fixed `L(s,χ)`.)
  - For `d ≥ 2`, simultaneous return then takes time `≳ U^{1−o(1)}`, under a Diophantine
    condition on the basis. This is a heuristic (Dirichlet-type counting), not a lemma.
  - In rank exactly `1`, `Γ_U ⊂ sℤ`, so the return time is `2π/s ≥ 2π#Γ_U/U`. The bound
    `≥ (1+o(1)) log U` would need almost all zeros simple, which is not known. §1 above
    gives `≫ log T` for the denominators from the explicit formula, without information on
    multiplicities. If the bound `(1+o(1))log U` held, (B) would exclude only
    `s ∈ (π/log U, 2π/log U]`.
- **The pair of exponents:** (B) needs the return time to be `≤ (2+o(1)) log U`, i.e.
  `(log U)^1`. Bounded-rank structure (`d ≥ 2`) gives at best `U^{1−o(1)}` (heuristically).
  They meet only for rank-one, AP-like configurations at the top scale with spacing above
  half the mean spacing. Neither Proposition 4.1 nor any known consequence of an atom
  (`HYP.md` Corollary 1.4, `INVLO.md` Corollary C′ and Proposition E) produces those.
- **What would close the gap** (none available):
  - (A′) an atom forces approximate AP structure of `Γ_U` at spacing `≥ π/((1−δ) log U)` for
    infinitely many `U`;
  - (B′) unconditional information on primes in windows of relative width `1/U` at
    prescribed `x ≥ U^{2+δ}`, which is beyond sieve methods;
  - (B″) a version of (B) that sees dependence of low zeros on high zeros rather than
    recurrence. The explicit formula is linear in the zero measure, while relations are
    multiplicative (§2), so no such version is in sight.

**What would revive the route:**
- (a) An unconditional theorem excluding rank-one `G` whose denominators grow like `log T`.
  This is a special case of Problem 1 of arXiv:2407.03530. By the description in
  `DENSITY.md` §3.3, that problem asks for two ordinates with irrational ratio, i.e. it
  excludes rank one altogether. It is open.
- (b) A structure theorem for measures with discrete support, density `≍ log T`, and
  Fourier transform "atomic on `log ℕ` plus absolutely continuous".
- (c) A proof that an atom forces a **cyclic** spectrum. Combined with Remark 7.6, which
  excludes exactly the cyclic `G`, this would settle the problem. The known structural
  consequences of an atom are far from that:
  - `HYP.md` Corollary 1.4;
  - Proposition 4.1 above;
  - `INVLO.md` Corollary C′ and Proposition E.
