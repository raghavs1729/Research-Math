# Short integer relations bound the atoms: a moment method

AI-produced (Anthropic's Claude), 2026-11-27. Labels as in `monotonicity/density/DENSITY.md`:
**Theorem / Proposition / Lemma** means a complete proof is given here. The main statement
(`ν_{1/2}({0}) = 0` under GRH alone) is **not** proved here. Nothing here concerns Problem 12 of
arXiv:2407.03530. No novelty is claimed: the literature was not searched beyond
`monotonicity/density/LITERATURE.md`. Status: **not yet refereed**.

This note answers the request to attack Conjecture U (predecessor branch, `DENSITY.md` §9,
after Lemma 9.13) by inverse Littlewood–Offord methods. The inverse principle — concentration
forces additive structure in the frequencies — takes a precise form for Haar measure on a
subtorus: **near-atoms force short integer relations.**
- Lemma A: if two laws share their first `2N` moments, and one of them has a density and
  sub-Gaussian tails, the other has small-ball probability `O(N^{−1/2})` at scale `N^{−1/2}`.
- Lemma D: density bounds for sums of cosines with independent phases.
- Proposition B′ (no hypothesis): the moments of a hull function of order `≤ L` see only the
  integer relations of length `≤ L`.
- Proposition B (GRH + NCZ): if all relations of length `≤ L` among the ordinates involve
  only finitely many ordinates, every atom of `ν_M` has mass `O(L^{−1/2})`, with an explicit
  dependence on those ordinates.
- Theorem C (GRH + NCZ + GLI(ψ) with `ψ(t)·log t/t → ∞`): a "graded LI" hypothesis gives no
  atoms, hence the densities. GLI(ψ) says every integer relation among the ordinates has
  length `≥ ψ(largest ordinate involved)`. It is incomparable with TS, FD and TS_ℤ (§3).
- Corollary C′ (GRH + NCZ): an atom forces infinitely many integer relations among the
  ordinates of length `≤ C·t/log t`, where `t` is the largest ordinate involved.
- Proposition E (GRH + NCZ): an atom of mass `α` forces *saturation*. For every `δ > 0` and
  every large `L`, all but a proportion `O(δ)` of the ordinates in every dyadic block from
  height `log³L` up to `C_2^{−1}α²δ·L log L` occur in integer relations of length `≤ L`.
- §4 is the assessment of Conjecture U by this route: it handles exactly the models without
  short relations, and the rest contains the cyclic toy question.

Numerical checks: `invlo_checks.py` (`invlo_checks_results.txt`).

Notation (`DENSITY.md` §1):
- `Γ` is the set of positive ordinates of `L(s,χ)`, with multiplicities `m_γ`;
- `r_γ = r_γ(M) = 4M/(M+iγ)` and `c_γ := m_γ r_γ` (the normalisation of `DENSITY.md` §1);
- `F_M = 1 + Re Σ_{γ∈Γ} c_γ e_γ` on the hull `A`, and `ν_M` is its law;
- `V_M := ½Σ_γ|c_γ|² = ‖F_M − 1‖²_2` and `V_M^{>T} := ½Σ_{γ>T}|c_γ|²`.

For a finitely supported `k ∈ ℤ^{(Γ)}` (or `ℤ^{Γ'}`), the **length** is `‖k‖_1 = Σ|k_γ|`, and
`k` is a **relation** if `Σ_γ k_γγ = 0`. For `L ≥ 1`, `S_L` is the set of frequencies that
occur (with `k_γ ≠ 0`) in some nonzero relation of length `≤ L`.

## 1. Two lemmas on real random variables

**Lemma A (moment comparison; no hypothesis).** Let `N ≥ 1`, and let `W, Z` be real random
variables with finite moments of order `2N` and `E Z^j = E W^j` for `0 ≤ j ≤ 2N`. Suppose:
- `W` has a density `f_W`;
- `W = μ + U + V` with `μ ∈ ℝ`, `|U| ≤ b` a.s., and `E|V|^{2N} ≤ 2·N!·(2σ²)^N` (`b ≥ 0`, `σ > 0`).

For `x ∈ ℝ` put `R := 16(|x−μ| + b + 2√N σ)` and `D_x := ess sup_{|t−x|≤R} f_W(t)`. Then

    P(|Z − x| ≤ 2R/(N+1)) ≤ 352·D_x R/(N+1) + 56·16^{−N}.

Note that `2R/(N+1) ≥ 32σ/√N`. The moment condition on `V` holds if
`E e^{λV} ≤ e^{λ²σ²/2}` for all real `λ`. Indeed then `P(|V| ≥ t) ≤ 2e^{−t²/(2σ²)}`, and
`E|V|^{2N} = ∫_0^∞ 2N t^{2N−1} P(|V| ≥ t) dt ≤ 2·N!(2σ²)^N`.

*Proof.* **A needle polynomial.** Let `F_N(φ) := Σ_{|k|≤N}(1 − |k|/(N+1))e^{ikφ}
= (N+1)^{−1}(sin((N+1)φ/2)/sin(φ/2))²`, the Fejér kernel, so `0 ≤ F_N ≤ N+1`. Put
`G(φ) := F_N(φ − π/2) + F_N(φ + π/2)`.
- `G` is an even trigonometric polynomial of degree `≤ N`. So `G(φ) = g(cos φ)` for a real
  polynomial `g` of degree `≤ N`.
- `g(0) = G(π/2) = (N+1) + F_N(π) ∈ [N+1, N+2]`, and `0 ≤ g ≤ 2(N+1)` on `[−1,1]`.
- `∫_{−1}^1 g² ds = ∫_0^π G(φ)² sin φ dφ ≤ 4∫_{−π}^{π} F_N² dφ`, and
  `∫_{−π}^{π}F_N² = 2πΣ_{|k|≤N}(1−|k|/(N+1))² ≤ 2π(N+1)`.

Put `q := g/g(0)`. Then:
- `q(0) = 1` and `0 ≤ q ≤ 2` on `[−1,1]`;
- `∫_{−1}^1 q² ≤ 8π/(N+1)`;
- `q ≥ 0.27` on `|s| ≤ 2/(N+1)`. Indeed, `|s| ≤ 2/(N+1) ≤ sin(π/(N+1))` gives
  `|arccos s − π/2| ≤ π/(N+1)`. For `|δ| ≤ π/(N+1)` we have
  `F_N(δ) ≥ (N+1)^{−1}((2/π)(N+1))² = 4(N+1)/π²`, using `sin u ≥ 2u/π` on `[0, π/2]`. So
  `q(s) ≥ 4(N+1)/(π²(N+2)) ≥ 8/(3π²) > 0.27`;
- `|q(s)| ≤ 2|T_N(s)| ≤ 2(2|s|)^N` for `|s| ≥ 1`. This is Chebyshev's inequality for a
  polynomial of degree `≤ N` bounded by `2` on `[−1,1]`, with `T_N` the Chebyshev polynomial.

`invlo_checks.py` (1) confirms these properties numerically for `N ≤ 300`.

**The bound.** `p(t) := q((t−x)/R)²` is a polynomial of degree `≤ 2N` with `p ≥ 0` on `ℝ`, and
`p ≥ 0.27² > 1/14` on `|t−x| ≤ 2R/(N+1)`. So

    P(|Z−x| ≤ 2R/(N+1)) ≤ 14·E p(Z) = 14·E p(W).

The equality uses only the moments of order `≤ 2N`. Split `E p(W)` at `|W − x| = R`:
- `E[p(W); |W−x| ≤ R] ≤ D_x R ∫_{−1}^1 q² ≤ 8πD_xR/(N+1)`;
- `E[p(W); |W−x| > R] ≤ 4·4^N E|W−x|^{2N}/R^{2N}`. Here `|W − x| ≤ β + |V|` with
  `β := |x−μ| + b`, so `E|W−x|^{2N} ≤ 2^{2N−1}β^{2N} + 2^{2N}N!(2σ²)^N`;
- since `R ≥ 16β`, the first part contributes `≤ 4·(16^N/2)·16^{−2N} = 2·16^{−N}`;
- since `R ≥ 32√Nσ` and `N! ≤ N^N`, the second contributes `≤ 4·32^N N!σ^{2N}/R^{2N} ≤ 4·32^{−N}
  ≤ 2·16^{−N}`.

Hence `P(|Z−x| ≤ 2R/(N+1)) ≤ 14(8πD_xR/(N+1) + 4·16^{−N})`, and `14·8π < 352`. ∎

*Sharpness of the rate.* `N^{−1/2}` cannot be improved. Let `W` be standard Gaussian and `Z`
its `(N+1)`-point Gauss–Hermite quadrature, which matches the moments up to order `2N+1`. Then
`Z` has an atom of mass `≈ 1.24/√N` (`invlo_checks.py` (4)). The constants in Lemma A are
crude; only the rates matter below.

**Lemma D (densities of random cosine sums; no hypothesis).** Let `J` be finite, `c_γ ∈ ℂ`, and
`Y := Re Σ_{γ∈J} c_γ e^{iθ_γ}`, with `θ_γ` independent and uniform on `[0,2π)`.
- (a) If at least three `γ ∈ J` have `|c_γ| ≥ a > 0`, then `Y` has a density `≤ 6/(π²a)`.
- (b) Let `I ⊂ J` with `#I ≥ 3` and `a ≤ |c_γ| ≤ b` for `γ ∈ I`. Put
  `θ := max(e^{−a²/b²}, π^{−1/2})` and `Σ_I := Σ_{γ∈I}|c_γ|²`. Then `Y` has a density
  `≤ (2π)^{−1}[2√π Σ_I^{−1/2} + 12θ^{#I−3}/(πa)]`.

*Proof.* The characteristic function of `Y` is `φ_Y(ξ) = Π_{γ∈J} J_0(ξ|c_γ|)`, and a density
bound is `(2π)^{−1}∫|φ_Y|`. Two facts about `J_0` are used.
- `|J_0(y)| ≤ (2/(π|y|))^{1/2}` for `y ≠ 0`. By Nicholson's formula, `y(J_0² + Y_0²)(y)`
  increases to `2/π`.
- `0 < J_0(y) ≤ e^{−y²/4}` for `|y| < j_{0,1} ≈ 2.405`. Use the product
  `J_0(y) = Π_k(1 − y²/j_{0,k}²)`, the inequality `1 − u ≤ e^{−u}`, and Rayleigh's sum
  `Σ_k j_{0,k}^{−2} = 1/4`.

Both are checked numerically in `invlo_checks.py` (2) and (5).

(a) Bound the three factors by `min(1, (2/(π|ξ|a))^{1/2})` and the others by `1`. Then
`∫|φ_Y| ≤ ∫ min(1, (2/(π|ξ|a))^{3/2}) dξ = 12/(πa)`.

(b) Split the integral at `|ξ| = 2/b`.
- For `|ξ| ≤ 2/b`, every `γ ∈ I` has `ξ|c_γ| ≤ 2 < j_{0,1}`. So
  `|φ_Y| ≤ exp(−ξ²Σ_I/4)`, whose integral over `ℝ` is `2√π Σ_I^{−1/2}`.
- For `|ξ| > 2/b`, take `γ ∈ I`.
  - If `ξ|c_γ| ≤ 2`, then `|J_0(ξ|c_γ|)| ≤ e^{−ξ²a²/4} < e^{−a²/b²}`.
  - Otherwise `|J_0(ξ|c_γ|)| < π^{−1/2}`.

  So every factor from `I` is `≤ θ`. Keep three of them for decay, as in (a):
  `|φ_Y(ξ)| ≤ θ^{#I−3} min(1, (2/(π|ξ|a))^{3/2})`. Its integral is `≤ 12θ^{#I−3}/(πa)`. ∎

## 2. Hull functions

**Proposition B′ (finite models; no hypothesis).** Let `Γ' ⊂ (0,∞)` be finite and `c_γ ∈ ℂ`.
Put `F := 1 + Re Σ_{γ∈Γ'} c_γ e_γ` on the closure torus `A_{Γ'}` of `{(e^{iγy})_{γ∈Γ'} : y ∈ ℝ}`,
with Haar measure `λ`. Let `L ≥ 2`, `N := ⌊L/2⌋`, and `S_L ⊂ Γ'` as above. Put
`B := Σ_{γ∈S_L}|c_γ|` and `σ_Y² := Σ_{γ∈Γ'∖S_L}|c_γ|²`. Let `D` be a density bound for
`Y := Re Σ_{γ∈Γ'∖S_L} c_γe^{iθ_γ}` (independent uniform phases), e.g. from Lemma D. Then for
every `x ∈ ℝ`, with `R := 16(|x−1| + B + 2√Nσ_Y)`,

    λ{|F − x| ≤ 2R/(N+1)} ≤ 352·D·R/(N+1) + 56·16^{−N}.

The same holds with `σ_Y` replaced by any `σ' ≥ σ_Y` in `R`, and `B` by any `B' ≥ B`.

*Proof.* Let `X := Re Σ_{γ∈S_L} c_γ e_γ` under Haar measure on `A_{S_L}`, independent of `Y`,
and put `W := 1 + X + Y`. Let `Z` be `F` under `λ`.

*Moments.* `Re(c e_γ) = ½(c e_γ + \bar c e_{−γ})`. So `(Z−1)^j` is `2^{−j}` times a sum, over
sequences `(γ_i, σ_i)_{i≤j} ∈ (Γ' × {±1})^j`, of `Π_i c_{γ_i}^{(σ_i)} e_{Σσ_iγ_i}`, where
`c^{(+1)} = c` and `c^{(−1)} = \bar c`. Since `∫ e_g dλ = 1_{g=0}`, the sequence contributes iff
`Σ_γ k_γγ = 0`, where `k_γ := Σ_{i:γ_i=γ}σ_i` and `‖k‖_1 ≤ j`.

The same expansion holds for `W − 1`. There a sequence contributes iff `k_γ = 0` for all
`γ ∉ S_L` and `Σ_{γ∈S_L}k_γγ = 0`. For `j ≤ L` the two conditions agree: a nonzero relation of
length `≤ L` is supported in `S_L`, by definition. Hence `E Z^j = E W^j` for `j ≤ L`, in
particular for `j ≤ 2N`. (`invlo_checks.py` (3) checks this identity in exact arithmetic on
an example, and shows that it fails at order `L+1` when a longer relation is present.)

*Lemma A applies* with `μ = 1`, `U = X`, `V = Y`, `b = B` and `σ = σ_Y`:
- `|X| ≤ B`;
- each term of `Y` has mean `0` and modulus `≤ |c_γ|`, so by Hoeffding's lemma and
  independence `E e^{λY} ≤ e^{λ²σ_Y²/2}`;
- `W` has a density `≤ D`, because `Y` does and is independent of `X`. ∎

**Proposition B (GRH + NCZ).** Fix `M > 0`, `c ∈ ℝ`, `L ≥ 2`, and `N := ⌊L/2⌋`. Suppose `S_L` (for
the ordinates `Γ`) is finite. Put `B := Σ_{γ∈S_L}|c_γ|` and
`R := 16(|c−1| + B + 2(2NV_M)^{1/2})`. Let `D` bound the density of
`Re Σ_{γ∈J}c_γe^{iθ_γ}` (independent uniform phases) for every finite `J ⊃ J_0`, for some fixed
finite `J_0 ⊂ Γ∖S_L`. By Lemma D the bound depends only on `J_0`. Then

    ν_M({c}) ≤ 352·D·R/(N+1) + 56·16^{−N}.

In particular, suppose the ordinates satisfy **no** nonzero integer relation of length `≤ L`.
Then `ν_M({c}) ≤ C(M,χ)(1 + |c−1|)L^{−1/2}`: take `J_0` the first three ordinates and Lemma
D(a).

*Proof.*
- Let `T` exceed `S_L` and `J_0`.
- Relations among `Γ_T` of length `≤ L` are relations among `Γ`, so `S_L(Γ_T) ⊂ S_L`. Hence,
  for `Γ' := Γ_T`:
  - the quantity `B` of Proposition B′ is `≤ B`;
  - `σ_Y² ≤ Σ_{γ≤T}|c_γ|² ≤ 2V_M`;
  - `Γ_T∖S_L(Γ_T) ⊃ J_0`, so the free sum has density `≤ D`.
- Apply Proposition B′ to the truncation `F_{M,T}` on `A_T`, with `B' := B` and
  `σ' := (2V_M)^{1/2}`, so that `R` does not depend on `T`. The law of `F_{M,T}` equals that
  of `F_{M,T}∘π_T` under `λ_A` (`DENSITY.md` Lemma 9.13, predecessor branch).
- Put `ρ := 2R/(N+1)`. Then
  `{|F_M − c| < ρ/2} ⊂ {|F_{M,T}∘π_T − c| ≤ ρ} ∪ {|F_M − F_{M,T}∘π_T| ≥ ρ/2}`. Under GRH + NCZ,
  `ν_M` is the law of `F_M` (`DENSITY.md` §1), so by Chebyshev's inequality

      ν_M({c}) ≤ 352·D·R/(N+1) + 56·16^{−N} + 4V_M^{>T}/ρ².

- Let `T → ∞`. The bound is unchanged except that `V_M^{>T} → 0`. ∎

**Hypothesis GLI(ψ) ("graded linear independence").** `ψ : [0,∞) → [0,∞)` is nondecreasing.
Every nonzero integer relation among the ordinates has length at least `ψ` of the largest
ordinate it involves: `Σk_γγ = 0`, `k ≠ 0` ⇒ `‖k‖_1 ≥ ψ(max{γ : k_γ ≠ 0})`.

So relations are allowed, but a relation reaching height `t` must be long, of length
`≥ ψ(t)`. LI implies GLI(ψ) for every ψ. Searches for small integer relations among computed
zeros test statements of this kind, up to a height.

**Theorem C (GRH + NCZ + GLI(ψ) with `ψ(t)·log(t+3)/t → ∞`).** For every `M > 0`, `ν_M` has no
atoms. Consequently the conclusions of `HYP.md` Theorem 1.3 hold. That is, for every
`m > −1/2`:
- the logarithmic densities of `{D_m > 0}` and `{D_m ≥ 0}` exist;
- near-ties have logarithmic density `0`.

For `χ_{−4}`, under GRH + GLI(ψ) the logarithmic density of `{x : π(x;4,3) > π(x;4,1)}` exists.

*Proof.* Fix `c` and `L ≥ 2`. Let `T_L := sup({0} ∪ {t ≥ 0 : ψ(t) ≤ L})`, which is finite because
`ψ → ∞`.
- *`S_L` is low.* If `k ≠ 0` is a relation with `‖k‖_1 ≤ L`, then `ψ(max supp k) ≤ L`, so
  `supp k ⊂ (0, T_L]`. Hence `S_L ⊂ Γ_{T_L}`, a finite set.
- *Its mass.* `B ≤ Σ_{γ≤T_L}|c_γ| ≤ 4MΣ_{γ≤T_L}m_γ/γ ≤ C_χ M(log(T_L+3))²`. This is partial
  summation with `N(t,χ) ≪ t log(t+3)`.
- *The free part has a small density.* We need `D ≤ C(M,χ)(1 + (T_L/log(T_L+3))^{1/2})`. Two
  facts come from the Riemann–von Mangoldt formula for `L(s,χ)`:
  - multiplicities satisfy `m_γ ≤ C_χ log(γ+3)`;
  - for `t ≥ t_χ`, `Σ_{γ∈(t,2t]} m_γ ≥ c_χ t log t`. Hence `(t, 2t]` contains at least
    `c_χ t log t/(C_χ log(2t+3)) ≥ c′_χ t` distinct ordinates.

  Zeros are not assumed simple.

  For `T_L ≥ max(t_χ, M)`, apply Lemma D(b) with `J_0 = I := Γ ∩ (T_L, 2T_L]`, which lies
  outside `S_L`:
  - `a := min_I|c_γ| ≥ |r_{2T_L}| ≥ M/(2T_L)` and `b := max_I|c_γ| ≤ 4MC_χ log(2T_L+3)/T_L`.
    So `a/b ≥ (8C_χ log(2T_L+3))^{−1}`. Here `t_χ` is also chosen so that `#I ≥ 3`.
  - `#I ≥ c′_χT_L`, so `θ^{#I−3} ≤ exp(−c″_χ T_L/log²T_L)` for large `T_L`. Hence the second
    term of Lemma D(b), `12θ^{#I−3}/(πa) ≤ (24T_L/(πM))exp(−c″_χT_L/log²T_L)`, is
    `O_χ(M^{−1})`.
  - `Σ_I = Σ_I m_γ²|r_γ|² ≥ (M/(2T_L))²Σ_I m_γ ≥ c_χM² log T_L/(4T_L)`. So the first term of
    Lemma D(b) is `≤ C_χ M^{−1}(T_L/log T_L)^{1/2}`.

  For smaller `T_L`, use Lemma D(a) with three ordinates in `(t_χ', t_χ' + H_χ]`, where
  `t_χ' := max(t_χ, M)` and `H_χ` is such that every interval of length `H_χ` contains three
  ordinates. On any bounded range of `T_L`, either bound is `O_{M,χ}(1)`.
- *Proposition B*, with `N = ⌊L/2⌋ ≥ L/3`, gives

      ν_M({c}) ≤ C(c,M,χ)·(1 + (T_L/log(T_L+3))^{1/2})·[(1 + log²(T_L+3))/L + L^{−1/2}] + 56·16^{−⌊L/2⌋}.

- *The limit.* Given `K > 0`, `ψ(t) ≥ Kt/log(t+3)` for `t ≥ t_K`.
  - If `T_L > t_K + 1`, pick `t ∈ (T_L−1, T_L]` with `ψ(t) ≤ L`. Then
    `L ≥ K(T_L−1)/log(T_L+2)`.
  - Then `(T_L/log(T_L+3))^{1/2}L^{−1/2} ≤ C K^{−1/2}`. Also `log T_L ≤ 2 log L` for large
    `L`, so the term with `log²(T_L+3)/L` tends to `0`.
  - If `T_L` stays bounded along a sequence of `L`, the whole bound tends to `0` along it.

  `K` is arbitrary, so `ν_M({c}) = 0`.

The consequences follow as in `HYP.md` Theorem 1.3, i.e. as in `DENSITY.md` Theorems 2.1(d)
and 3.4. NCZ is automatic for `q = 4`. ∎

**Corollary C′ (what an atom forces; GRH + NCZ).** Suppose `ν_M({c}) > 0` for some `M > 0`,
`c ∈ ℝ`. For a relation `k`, write `t_k := max{γ : k_γ ≠ 0}`. Then there is `C < ∞` such that
infinitely many nonzero integer relations `k` among the ordinates satisfy

    ‖k‖_1 ≤ C·t_k/log(t_k+3).

In particular `t_k → ∞` along them, since only finitely many relations with bounded `t_k` have
bounded length.

*Proof.* Suppose not. Then for every integer `C ≥ 1` the set of relations with
`‖k‖_1 ≤ Ct_k/log(t_k+3)` is finite. Let `τ_C` be the largest `t_k` in it (`0` if it is empty);
`τ_C` is nondecreasing in `C`.
- Put `C(t) := max{C ∈ ℕ : 1 ≤ C ≤ t, τ_C < t}`, with `C(t) := 0` if there is no such `C`. It
  is finite and nondecreasing. It tends to `∞`, because every `τ_C` is finite.
- Put `ψ(t) := C(t)·t/log(t+3)`. It is nondecreasing, and `ψ(t)log(t+3)/t = C(t) → ∞`.
- GLI(ψ) holds. Take a relation `k` with `t_k = t`.
  - If `C(t) = 0`, then `ψ(t) = 0`.
  - Otherwise `C := C(t)` has `τ_C < t_k`, so `k` is not in the finite set for `C`. Hence
    `‖k‖_1 > Ct_k/log(t_k+3) = ψ(t_k)`.

Theorem C then gives `ν_M({c}) = 0`, a contradiction. ∎

**Proposition E (an atom saturates the spectrum with short relations; GRH + NCZ).** Fix
`M > 0`, `c ∈ ℝ` and `δ ∈ (0,1]`. There is a constant `C_2` (depending on `c, M, χ`) with the
following property. Suppose `ν_M({c}) = α > 0`. Then for every `L ≥ L_0(α, c, M, χ, δ)` and
every `T` with `log³L ≤ T ≤ C_2^{−1}α²δ·L log L`,

    Σ_{γ ∈ (T,2T] ∖ S_L} m_γ < δ·T log T.

Since `Σ_{γ∈(T,2T]} m_γ ~ (T/2π)log T`, this says: all but a proportion `< 2πδ(1+o(1))` of the
ordinates in `(T, 2T]` (with multiplicity) occur in integer relations of length `≤ L`. This
holds for every dyadic block between `log³L` and `C_2^{−1}α²δ·L log L`, and it is nontrivial
only for `δ < 1/(2π)`. The upper end of the range shrinks linearly with `δ`. Corollary C′ is
the weak form "`S_L ⊄` any bounded set".

*Proof.* Let `N := ⌊L/2⌋` and `T* := N⁴`. Suppose `I := (Γ∖S_L) ∩ (T,2T]` has
`Σ_I m_γ ≥ δT log T`, with `T` in the stated range (so `2T ≤ T*` for large `L`).

Apply Proposition B′ to `Γ' := Γ_{T*}`. As in Proposition B, `S_L(Γ_{T*}) ⊂ S_L`, and we take:
- `B' := Σ_{γ≤T*}|c_γ| ≤ C_χM log²(T*+3)`;
- `σ' := (2V_M)^{1/2}`;
- `D` from Lemma D(b) with this `I`.

The Lemma D(b) data, as in the proof of Theorem C:
- `a ≥ M/(2T)` and `b ≤ 4MC_χ log(2T+3)/T`;
- `#I ≥ δT log T/(C_χ log(2T+3)) ≥ δT/(2C_χ)`;
- `Σ_I ≥ (M/(2T))²δT log T`.

Hence `D ≤ C(δ^{−1}T/log T)^{1/2}M^{−1} + C(T/M)exp(−c_χδT/log²T)`. Also
`R/(N+1) ≤ C(1 + log²N + √N)/N ≤ C/√N`. Proposition B′ and the truncation step of Proposition
B give

    α ≤ C(T/(δ L log T))^{1/2} + C·T·exp(−c_χδT/log²T)·L^{−1/2} + 56·16^{−N} + 4V_M^{>N⁴}(N+1)²/(1024·2V_M).

The last two terms tend to `0` as `L → ∞`; here `ρ ≥ 32(2V_M)^{1/2}/√N`.
- *First term.* If `T ≤ √L`, then `T/log T ≤ √L`. If `√L < T ≤ C_2^{−1}α²δL log L`, then
  `log T ≥ ½log L` and `T/log T ≤ 2C_2^{−1}α²δL`. Either way the first term is `≤ α/3` for
  large `L`, once `C_2` is large.
- *Second term.* `T/log²T` increases for `T ≥ e²`. So `T ≥ log³L` gives
  `T/log²T ≥ log³L/(9(log log L)²) ≥ 4(c_χδ)^{−1}log L` for large `L`. With `T ≤ L²`, the second
  term is then `≤ C·L²·L^{−4}·L^{−1/2} ≤ α/3` for large `L`.

So `α < α` for large `L`, a contradiction. ∎

*Remarks.*
- (1) *Where `t/log t` comes from.* The moment method resolves scale `L^{−1/2}`, which is
  sharp by the Gauss–Hermite example after Lemma A. The free part of the spectrum above `T_L` spreads the law over scale
  `≍ (log T_L/T_L)^{1/2}`. The argument works when the first scale is finer:
  `T_L ≪ L log T_L`. The low part `S_L` enters only through its sup norm `B`, which costs
  nothing here.
- (2) Proposition B is effective. Suppose a computation certifies that the ordinates up to `T`
  satisfy no integer relation of length `≤ L`. Then Proposition B′ for `Γ_T`, plus the
  truncation step, bounds the atoms of `ν_M` by an explicit function of `T`, `L` and the
  computed ordinates. The constants of Lemma A are far too large for this to beat the
  certified bound `0.0891` of `DENSITY.md` Theorem 5.4. For a Gaussian `W` at `x = 0`, the
  bound of Lemma A drops below `1` only for `N > 2·10⁷` (referee 4, `referee4/r4_lemmaA.py`).
  Better constants were not pursued.

## 3. Comparison with the other hypotheses

All examples are abstract frequency sets, as in `HYP.md` §2; there is no claim that ordinates
look like this.
1. **LI implies GLI(ψ)** for every `ψ`.
2. **GLI(ψ) implies none of TS, FD, TS_ℤ**, for any `ψ`.
   - *The set.* Take `α > 0` and integers `0 ≤ e_1 < e_2 < …` with
     `e_n − e_{n−1} ≥ log₂ ψ(n+1)`. Take distinct odd primes `u_n` with
     `x_n := αu_n/2^{e_n} ∈ (n, n+1)`; for `e_n` large such primes exist, by the prime number
     theorem.
   - *GLI.* Let `k ≠ 0` be a relation, `n*` the largest index with `k_{n*} ≠ 0`, and `E'` the
     largest `e_n` over the other indices of the support (there are some, since `x_{n*} ≠ 0`).
     Multiplying `Σk_nu_n2^{−e_n} = 0` by `2^{e_{n*}}` gives
     `2^{e_{n*}−E'} | k_{n*}u_{n*}`. As `u_{n*}` is odd,
     `‖k‖_1 ≥ |k_{n*}| ≥ 2^{e_{n*}−e_{n*−1}} ≥ ψ(n*+1) ≥ ψ(x_{n*})`.
   - *TS and FD fail*, because `X ⊂ ℚα`, as for `{2^kα}` in `HYP.md` §2.
   - *TS_ℤ fails.* For `m < m'`, since `u_{m'}` is an odd prime `≠ u_m`,
     `⟨x_m, x_{m'}⟩ = α2^{−e_{m'}}gcd(u_m2^{e_{m'}−e_m}, u_{m'})ℤ = α2^{−e_{m'}}ℤ`. Hence
     `⟨X_{>T}⟩ = αℤ[1/2] = ⟨X⟩` for every `T`: the index is `1`.
3. **TS (hence TS_ℤ and FD) does not imply GLI(ψ)** once `ψ > 3` eventually. Take
   `X = {β} ∪ {nα : n ≥ 1}` with `β ∉ ℚα`.
   - TS holds: `β ∉ ⟨X_{>T}⟩_ℚ = ℚα` for `T ≥ β`.
   - GLI fails: `(n+1)α − nα − α = 0` has length `3` and reaches height `(n+1)α → ∞`.
4. So, for `ψ > 3` eventually, GLI(ψ) is incomparable with each of TS, FD and TS_ℤ. For
   `ψ ≤ 3` it holds for every set, so TS implies it. For `ψ(t)log(t+3)/t → ∞`, the disjunction
   "TS_ℤ or FD or GLI(ψ)" implies continuity of `ν_M` (Theorem C here, `HYP.md` Theorem 1.3,
   `DENSITY.md` Theorem 3.4).
5. **Mechanism.** TS and TS_ℤ ask for *one* direction in which the tail of the spectrum is
   free: a `ℚ`-linear form, or a finite quotient. FD asks for the same at finite heights. GLI
   asks for no free direction. It only asks that the relations be long compared with the
   heights they reach. In the example of item 2 the spectrum has rank one.

## 4. What this gives for Conjecture U

Recall the class `𝒞(M, A_0, s)` (`DENSITY.md` Lemma 9.13 ff., predecessor branch). It
consists of finite `Γ'` with integer multiplicities, the counting law, and
`Σm_γ8M²/(M²+γ²) ≤ s`, with **no condition on relations**. Conjecture U asks for
`ω(ε) := sup_𝒞 λ{|F_{Γ'}| < ε} → 0`.
- **A uniform bound on the class.** The counting law gives `m_γ ≤ A_0 log(γ+3)`. Hence, over
  `𝒞`,
  `Σ|c_γ|² = Σm_γ²|r_γ|² ≤ A_0Σm_γ log(γ+3)|r_γ|² ≤ C(M, A_0)`.
- **The structured/unstructured split.** Apply Proposition B′ at `x = 0`, with `D` from Lemma
  D(a), to a model in `𝒞` with `λ{|F_{Γ'}| < ε} ≥ η`.
  - Since `R ≥ 16`, the radius `2R/(N+1)` is `≥ ε` whenever `N + 1 ≤ 32/ε`.
  - Take such an `N` with `56·16^{−N} ≤ η/2`.

  Then one of the following holds:
  - fewer than three frequencies outside `S_L` have coefficient `≥ a`;
  - `352·(6/(π²a))·R/(N+1) ≥ η/2`, i.e. `1 + B + 2√Nσ_Y ≳ ηaN`, with
    `σ_Y² ≤ Σ|c_γ|² ≤ C(M,A_0)`.

  In the second case, `B ≳ ηaL`, i.e. the `ℓ¹`-mass of `S_L` grows linearly in `L`. This is
  derived only in the window

      C(M,A_0)·(ηa)^{−2} ≤ L ≲ 64/ε     (with 56·16^{−⌊L/2⌋} ≤ η/2).

  The window is nonempty once `ε` is small. Below it the second case may hold with `B = 0`.

  So a counterexample to U must be **relation-saturated** at every length `L` in that window.
  Either the short relations absorb all but two of the frequencies with coefficient `≥ a`, or
  they involve frequencies of `ℓ¹`-mass `≳ ηaL`. Models that are not saturated in this sense
  obey the small-ball bound of Proposition B′.
- **Saturation of cyclic models, and the toy.**
  - Cyclic models with **all** frequencies `β, 2β, …, Nβ` present (all multiplicities `≥ 1`)
    are saturated at `L = 3`, by `β + nβ − (n+1)β = 0`. These include the models of the
    predecessor's draft Theorem 9.22. The cyclic case of Question P″ (`DENSITY.md` Remark
    9.11, the predecessor's line (2)) also allows missing frequencies (`m_n = 0`); it is covered
    by the next bullet.
  - An arbitrary `Γ' ⊂ βℕ` need not be saturated at small `L`. For example, `{β, 3β, 5β}` has
    no relation of length `≤ 3` (referee 4, `REFEREE-4.md` I18). In general `nβ ∈ S_L` once
    `L ≥ (n+n')/gcd(n,n')` for some other `n'β ∈ Γ'`. So saturation holds at lengths
    `≳ 2·max{n : nβ` has coefficient `≥ a}`. That is uniform when `β` is bounded below, as in
    the cyclic case of Question P″, but not uniformly over `𝒞`.
  - On the models of the cyclic case of P″, `L ≍ 1/ε → ∞` lies in the saturated regime, and
    the moment method says nothing there.
- **Verdict.** The inverse Littlewood–Offord route proves the "unstructured half" of U
  (Proposition B′). Together with `LINE1.md` Proposition 4.1 (free sets of large-coefficient
  zeros are bounded), it gives the structure an atom needs: **many short relations among the
  large-coefficient ordinates, and few free directions.** For the true race this is
  Proposition E: for every `δ`, all but a proportion `O(δ)` of every dyadic block up to height
  `C_2^{−1}α²δ·L log L` is covered by relations of length `≤ L`.
  - Saturation alone is not contradictory. A spectrum made of independent triples
    `{a, b, a+b}` is saturated at `L = 3`. It has no atom, but for a different reason: TS
    holds there, and Jessen–Wintner applies to the independent blocks.
  - So a hypothetical atom needs saturation **and** the failure of TS, TS_ℤ and FD at once.

  The route cannot by itself prove U, because the saturated half contains the open cyclic
  toy question. For the true race the structure statement yields a sufficient hypothesis
  new here (Theorem C). The missing step is unchanged: excluding saturated integer spectra.
