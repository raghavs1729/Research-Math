# Mahler measure and Jensen's formula for the atom at `0`

Working paper, 2026-11-27, branch `claude/kind-dijkstra-c1a5gz-density3`. AI-produced (Anthropic's
Claude). **The main question (no atom at `0` under GRH alone, hence a logarithmic density for
Chebyshev's bias) is NOT answered here, and Conjecture U is NOT proved.** Nothing here concerns
Problem 12 of arXiv:2407.03530.

Labels as in `monotonicity/density/DENSITY.md` (cited below as **D**): **Theorem / Proposition /
Lemma / Corollary** = complete proof here from the quoted inputs; **Remark / Sketch / Heuristic** =
not a complete proof; **[num]** = numerical evidence (scripts in this directory). No LI anywhere.

## 0. Summary

The brief for this line: integer coefficients cannot be small on large sets (Mahler measure,
Jensen's formula, Kronecker, Lehmer, Erdős–Turán). Use this to attack Conjecture U of **D** §9.

**Final summary (honest).**
- **Proved, refereed:**
  - the Jensen–Newton bound (Lemma 2.1, Cor 2.3);
  - under GRH + NCZ + hypothesis V, the limiting law `ν_M` is absolutely continuous, with
    `∫log|x−c| dν_M ≥ log(2Mm_{γ_0}/|M+iγ_0|)` and
    `ν_M((c−ε,c+ε)) ≤ C/log(1/ε)` (`C = 2.055` for `γ_1(χ_{−4})`, `M = ½`, simple zeros) (Theorem 3.3).
    Hence the logarithmic density of `{π(x;4,3) > π(x;4,1)}` exists under GRH + V. V is strictly
    weaker than TS;
  - a finite-level version, LV (Theorem 3.4);
  - a special case of Conjecture U: models with a Newton vertex at an ordinate `≤ Y`, with rate
    `1/log(1/ε)` (Corollary 3.5);
  - the deficit criterion (Props 4.1–4.2, Cor 4.3).
- **Killed:** the Mahler/Jensen route to U, for arguments that use vertex data. In rank one with a
  dense spectrum the only vertex is the top ordinate, which loses `log T`. The loss is real
  (Theorem 9.6). The same vertex data occur in classes where U fails, real weights and integers
  without the counting law (Remark 5.3). The arithmetic theory adds `O(1)`; Erdős–Turán points the
  wrong way; `2w_T` is not small on near-atom sets (§5).
- **Not provable with present tools:** V for an actual zero. It would solve Problem 1 of
  arXiv:2407.03530, and the explicit formula at finite resolution cannot see dominated relations
  (Remark 3.12).
- **Not proved:** the main theorem, Conjecture U, Conjecture L.

| # | Statement | Hypotheses | Status |
|---|---|---|---|
| 2.1 | **Jensen–Newton bound.** For a trigonometric polynomial on a compact abelian group, `∫ log|P| ≥ log|c_{g_0}|` for every vertex `g_0` of the Newton polytope of `P`. For a finite model this gives `λ{|F − c| < δ} ≤ (log(|M+iγ_0|/(2M|w_{γ_0}|)) + ½log⁺(|1−c|²+V) + 1/(2e)) / log(1/δ)` for every vertex `γ_0` (Cor 2.3). | none | proved; checked [num] |
| 3.3 | **A low vertex excludes atoms.** Hypothesis V(`γ_0`): in every nontrivial integer relation `Σ k_γ γ = 0` among the ordinates, `|k_{γ_0}| < Σ_{γ≠γ_0}|k_γ|` (equivalently: `γ_0` is a vertex of `conv(±Γ_T)` for every `T`). Then `∫ log|x−c| dν_M(x) ≥ log(2M m_{γ_0}/|M+iγ_0|)` for every `c`, `ν_M((c−ε,c+ε)) ≤ C/log(1/ε)` with explicit `C`, and `ν_M` is **absolutely continuous**; all the densities exist. V is implied by TS and strictly weaker (Example 3.7); not implied by FD (Example 3.8). | GRH, NCZ, V | proved |
| 3.4 | Finite-level version: `ν_M({c}) ≤ 2 liminf_T log(|M+iμ_T|/(2M))/τ_T`, `μ_T` = smallest vertex ordinate at level `T` (for the race: `= 2 liminf log⁺μ_T/log T`). So "low vertices" (`log μ_T = o(τ_T)` along a sequence) exclude atoms; this covers also super-lacunary spectra of rank one (Remark 3.9). | GRH, NCZ (or any hull function) | proved |
| 3.5 | **A special case of Conjecture U**: on the subclass of models with a vertex at an ordinate `≤ Y`, `sup λ{|F| < ε} ≤ C(M,A_0,Y)/log(1/ε)`. | none | proved |
| 4.1–4.2 | **Deficits.** With the Mahler deficit `δ_T(c) := −∫log|F_{M,T} − c|` and `τ_T := log(1/V^{>T}_M) = log T + O(log log T)`: an atom of mass `η` at `c` forces `δ_T(c) ≥ (η/2)τ_T − O(1)` for **every** `T` (Prop 4.2), hence `ν_M({c}) ≤ 2 liminf δ_T(c)/τ_T`; Prop 4.1 is the small-ball form. The factor `2` is not known to be sharp. | GRH, NCZ (or any hull function) | proved |
| 5 | **Killed: the Mahler/Jensen route to Conjecture U** (for arguments that apply Lemma 2.1 to `F − c` or to polynomials with the same vertex data; Remark 5.3). Without a low vertex, in particular in rank one with a dense spectrum (`max Γ_T ≥ T^{1−o(1)}`), where every known atom model lives, Lemma 2.1 gives only `δ_T ≤ log T + O(1)`, i.e. "`ν({c}) ≤ 2`". The loss is real: the integer models of **D** Thm 9.6 have `δ_T ≥ 0.05 log T − O(log log T)` (Prop 5.2). The vertex data that Lemma 2.1 uses are the same for classes in which U fails (Remark 5.3): they see neither the counting law nor positivity, nor integrality beyond `|m| ≥ 1` at one point. The arithmetic theory (Kronecker, Lehmer, Smyth, Dobrowolski) adds `O(1)` at best and is beaten by the trivial vertex bound `M(2w) ≥ 2`. Erdős–Turán bounds zeros away from the circle from **above**, while U needs a bound from **below** (§5.4). The integer polynomial `2w_T` is not small on near-atom sets (§5.1). | — | proofs + [num] |

**In one sentence.** Of the theory of Mahler measure, only the analytic core (Jensen's formula
along a one-parameter subgroup) acts on the race. It acts through the Newton polytope of the
truncated hull function in `G_T ⊗ ℝ`, not through integrality. It settles every configuration with
a low vertex (Theorems 3.3–3.4). This includes configurations covered neither by the literature
hypotheses tabulated in **D** §3.3 nor by **D** Theorems 3.1, 3.4, 3.5 (Example 3.7). It is blind
when there is no low vertex. That happens in rank one with a dense spectrum, which is the case that
matters.

A **sufficient** statement in this language (§7; it implies U, the converse is not shown): in rank
one, for integer models obeying the counting law, the self-inversive polynomial `z^N F_N(z)` keeps
`log N − o(log N)` worth of zeros off the unit circle (`Σ_{|ζ|>1} log|ζ| ≥ log(1/|C_N|) − o(log N)`).
Lemma 2.1 says nothing about this. By a heuristic of the sibling session (branch `-density`,
`DENSITY.md` Remark 9.21), even this may fail, together with U (§7).

---------------------------------------------------------------------------

## 1. Setting and quoted inputs

Notation of **D** §1: `χ` real primitive, `Γ` the distinct positive ordinates of `L(s,χ)`, `m_γ ≥ 1`
their multiplicities, `r_γ(M) = 4M/(M+iγ)`, `A` the hull (closure of `ι(ℝ)` in `𝕋^Γ`), `λ_A` its Haar
measure, `e_γ` the coordinate characters, `G ⊂ ℝ` the group generated by `Γ` (the dual of `A`),
`F_M = 1 + Re Σ_γ m_γ r_γ e_γ = 1 + Σ_{g∈±Γ} m_g (2M/(M+ig)) e_g` (`m_{−γ} := m_γ`), `ν_M` its law.
Under GRH(`χ`) + NCZ(`χ`), `ν_M` is the limiting logarithmic distribution of
`E_M(x) = 2M x^{−M}(log x) Σ_{p≤x} χ*(p) p^m`, `M = m + ½` (**D** (I1)). `V_M := ½Σ_γ m_γ²|r_γ|² =
Σ_γ m_γ² 8M²/(M²+γ²) = ‖F_M − 1‖_2²`.

Truncations: `Γ_T := Γ ∩ (0,T]`, `G_T` the group generated by `Γ_T`, `A_T` its dual torus,
`π_T : A → A_T` the restriction (onto, `π_{T*}λ_A = λ_{A_T}`), `F_{M,T} := 1 + Re Σ_{γ∈Γ_T} m_γ r_γ e_γ`,
`V_M^{>T} := ½Σ_{γ>T} m_γ²|r_γ|² = ‖F_M − F_{M,T}∘π_T‖_2²`.

**Finite models.** A finite set `Γ' ⊂ (0,∞)`, real weights `w_γ ≠ 0` (`γ ∈ Γ'`), `G'` the group
generated by `Γ'`, `A'` its dual torus, `F_{Γ'} := 1 + Re Σ_{γ∈Γ'} w_γ r_γ e_γ`,
`V_{Γ'} := ½Σ w_γ²|r_γ|²`. The class `𝒞(M, A_0, s)` of **D** §9 consists of the finite models with
integer weights `m_γ ≥ 1`, the counting law `Σ_{γ∈Γ'∩(t,t+1]} m_γ ≤ A_0 log(t+3)` (`t ≥ 0`) and
`Σ m_γ 8M²/(M²+γ²) ≤ s`; Conjecture U says `ω(ε) := sup_{𝒞(½,A_0,0.1556)} λ_{A'}{|F_{Γ'}| < ε} → 0`.

**Mahler measure.** On a probability space, for `f ∈ L¹`, `m(f) := ∫ log|f| ∈ [−∞, ∞)` (well
defined, since `log⁺|f| ≤ |f|`). For a polynomial, `exp m(P)` is its Mahler measure `M(P)`.

**Newton polytope.** Let `G` be a torsion-free abelian group, so that `G → G ⊗_ℤ ℝ` is injective. For a
trigonometric polynomial `P = Σ_{g∈S} c_g e_g` (`S` finite, all `c_g ≠ 0`) on the dual group of `G`,
`New(P) := conv(S) ⊂ G ⊗ ℝ`. A point `g_0 ∈ S` is a **vertex** if `g_0 ∉ conv(S ∖ {g_0})`. Equivalently
there is a homomorphism `φ : G_S → ℤ`, `G_S` the subgroup generated by `S`, with `φ(g_0) > φ(g)` for
all `g ∈ S ∖ {g_0}`. (Indeed `G_S ≅ ℤ^d` and `G_S ⊗ ℝ ≅ ℝ^d`. A point outside the convex hull of a finite
set is strictly separated from it by a linear form. The separating forms form an open set, which
contains a form with rational coefficients in a basis of `G_S`, and a positive multiple of that form
is integral.)

**Deficits.** `δ_T(c) := −m(F_{M,T} − c)` (on `A_T`), `τ_T := log(1/V_M^{>T})`. For the race, `V_M^{>T} > 0`
for every `T`.

Quoted inputs.

- **(J) Jensen's formula** (e.g. W. Rudin, *Real and Complex Analysis*, 3rd ed., 15.17–15.18; numbering
  recalled). If `p(z) = a Π_{k=1}^n (z − ζ_k)` with `a ≠ 0`, then
  `(1/2π)∫_0^{2π} log|p(e^{is})| ds = log|a| + Σ_k log max(1, |ζ_k|)`. In particular
  `m(p) ≥ log|a|`. (This is Mahler's formula `M(p) = |a| Π max(1,|ζ_k|)`.)
- **(D1)** **D** (I1), (I2), Lemma 1.1 (a nonzero trigonometric polynomial on `A` vanishes only on a
  null set), Theorem 3.5 and its proof (tame directions), and the deductions "no atom ⇒ densities
  exist" of **D** Theorem 2.1(d) and Theorem 3.1.
- **(D2)** **D** Lemma 9.13, which holds for every level `c` with the same proof:
  `ν_M((c−ε, c+ε)) ≤ λ_{A_T}{|F_{M,T} − c| < 2ε} + V_M^{>T}/ε²`.
- **(D3)** **D** (3.1): `V_M^{>T} ≤ C_1 M² log²(T+3)/T`, `C_1 = C_1(χ)`.
- **(D4)** The models of **D** Example 9.2, Proposition 7.3 and Theorem 9.6, and the sharpness remarks
  after Lemma 9.13 (rescaled Example 9.2 at `M = 0.05` on `10ℕ`).
- **(RvM)** Riemann–von Mangoldt for `L(s,χ)`, `χ` primitive mod `q`: the number of zeros with
  `0 < β < 1`, `|γ| ≤ T` is `(T/π) log(qT/2π) − T/π + O(log qT)` (H. Davenport, *Multiplicative
  Number Theory*, ch. 16; statement recalled). For real `χ` the zeros are symmetric, so
  `N(T) := Σ_{0<γ≤T} m_γ = (T/2π) log(qT/2πe) + O(log qT)`.

---------------------------------------------------------------------------

## 2. Jensen's formula along a direction: the Jensen–Newton bound

**Lemma 2.1 (Jensen–Newton).** Let `A` be a compact abelian group whose dual is a torsion-free group
`G`, and `P = Σ_{g∈S} c_g e_g` a trigonometric polynomial on `A` (`S ⊂ G` finite, `c_g ≠ 0` for `g ∈ S`).
If `g_0` is a vertex of `conv(S)`, then

    m(P) = ∫_A log|P| dλ_A ≥ log|c_{g_0}|.

*Proof.* Let `G_S` be the subgroup generated by `S`, `A_S` its dual torus and `π : A → A_S` the
restriction map. It is a continuous surjective homomorphism, so `π_*λ_A = λ_{A_S}`. Moreover `P = Q∘π`
with `Q := Σ_{g∈S} c_g e_g` on `A_S`, hence `m(P) = m(Q)`. Choose `φ : G_S → ℤ` with
`φ(g_0) > φ(g)` for `g ∈ S ∖ {g_0}`. Then `j(s) := (g ↦ e^{iφ(g)s})` is a continuous homomorphism
`ℝ → A_S = Hom(G_S, 𝕋)`, and for `θ ∈ A_S`

    q_θ(s) := Q(θ·j(s)) = Σ_{g∈S} c_g e_g(θ) e^{iφ(g)s}.

Put `k_+ := φ(g_0) = max φ(S)` and `k_− := min φ(S)`. If `k_+ = k_−` then `S = {g_0}` and the claim is
trivial. Otherwise `|q_θ(s)| = |p_θ(e^{is})|` with `p_θ(z) := Σ_{g∈S} c_g e_g(θ) z^{φ(g)−k_−}`, a polynomial
of degree `k_+ − k_−` whose leading coefficient is `c_{g_0} e_{g_0}(θ)`, because `g_0` is the only
element of `S` with `φ(g) = k_+`. By (J),
`(1/2π)∫_0^{2π} log|q_θ(s)| ds ≥ log|c_{g_0}|` for **every** `θ`. The map `(θ,s) ↦ Q(θ j(s))` is
continuous and `λ_{A_S}` is invariant under translation by `j(s)`. Applying Tonelli to `log⁺|Q(θj(s))|`
(bounded) and to `log⁻|Q(θj(s))|` (nonnegative) separately gives

    m(Q) = ∫_{A_S} [ (1/2π)∫_0^{2π} log|Q(θ j(s))| ds ] dλ_{A_S}(θ) ≥ log|c_{g_0}|. ∎

**Lemma 2.2 (small ball from the Mahler measure).** Let `(X, μ)` be a probability space, `f ∈ L²(μ)`
with `m(f) > −∞`, `0 < δ < 1` and `η := μ{|f| < δ}`. Then

    η log(1/δ) ≤ −m(f) + ½ log⁺ ∫|f|² dμ + 1/(2e).

*Proof.* Let `B := {|f| < δ}`. Then `∫_B log|f| ≤ η log δ`. If `η < 1`, Jensen's inequality for the
concave `log`, applied to the probability measure `μ|_{B^c}/(1−η)`, gives
`∫_{B^c} log|f| = ((1−η)/2) avg_{B^c} log|f|² ≤ ((1−η)/2) log(∫|f|²/(1−η))`
`≤ ½log⁺∫|f|² + ((1−η)/2) log(1/(1−η))`, and `x log(1/x) ≤ 1/e` on `(0,1]`. Add the two pieces. If
`η = 1` then `m(f) ≤ log δ` directly. ∎

**Corollary 2.3 (finite models).** Let `F_{Γ'}` be a finite model, `γ_0 ∈ Γ'` a vertex of `conv(±Γ')`
in `G' ⊗ ℝ`, and `c ∈ ℂ`. Then

    m(F_{Γ'} − c) ≥ log(2M|w_{γ_0}|/|M+iγ_0|),
    λ_{A'}{|F_{Γ'} − c| < δ} ≤ [ log(|M+iγ_0|/(2M|w_{γ_0}|)) + ½log⁺(|1−c|² + V_{Γ'}) + 1/(2e) ] / log(1/δ)   (0 < δ < 1).

For integer weights, `|w_{γ_0}| ≥ 1`.

*Proof.* The support of `F_{Γ'} − c` is `±Γ'` together with `0` if `c ≠ 1`. Since `0` is the midpoint of
`±γ_0`, adding it does not change the convex hull, so `γ_0` is a vertex of the hull of the support.
Distinct ordinates are distinct elements of `G'`, so the coefficient of `e_{γ_0}` is
`w_{γ_0}·2M/(M+iγ_0)`. Lemma 2.1 gives the first bound. For the second, apply Lemma 2.2, using
`‖F_{Γ'} − c‖_2² = |1−c|² + V_{Γ'}`. ∎

**Remark 2.4 (what is used, and what is not).** (i) Lemma 2.1 sees `P` only through **one**
coefficient. Corollary 2.3 holds for real or signed weights, for any relations among the ordinates,
with or without the counting law. Integrality enters only through `|w_{γ_0}| ≥ 1`. (ii) For a fixed
exposing `φ`, equality in Lemma 2.1 holds iff for a.e. `θ` all roots of `p_θ` lie in the closed unit
disc. It does occur: in Example 3.7 below, at `c = 1` with all `m_n = 1`, `m(F_N − 1) = log|c_{γ_0}|`
exactly. (iii) For one variable the lemma is the classical `M(P) ≥ |leading coefficient|`. The
multivariable form (vertex coefficients) is also classical; the proof above, by restriction to a
one-parameter subgroup that exposes the vertex, is the standard one. (iv) The vertex restriction
matters. For `10 + 2cos x + 2cos y` the constant term is not a vertex, and indeed
`m = 2.2816 < log 10` (referee's check).

**Numerical check 2.5 [num]** (`jensen_newton_check.py` → `jensen_newton_check_results.txt`). For 24
random Laurent polynomials on `𝕋²`, `𝕋³` (3–11 terms, random complex coefficients; vertices by linear
programming; offset-grid quadrature), the minimum of `m(P) − max_v log|c_v|` is `−0.0000`
(quadrature error), and the maximum of the defect in Lemma 2.2 is `−0.19 < 0`.

---------------------------------------------------------------------------

## 3. A low vertex excludes atoms

**Definition 3.0.** For `γ_0 ∈ Γ`, **hypothesis V(`γ_0`)**: for every `T ≥ γ_0`, `γ_0` is a vertex of
`conv(±Γ_T) ⊂ G_T ⊗ ℝ`. (Since `G_T ⊗ ℝ → G_{T'} ⊗ ℝ → G ⊗ ℝ` are injective, this can be read in
`G ⊗ ℝ`.)

**Proposition 3.1 (equivalent forms of V).** For `γ_0 ∈ Γ` the following are equivalent:

(i) V(`γ_0`);

(ii) every nontrivial integer relation `Σ_{γ∈Γ} k_γ γ = 0` (finitely many `k_γ ≠ 0`) among the
distinct positive ordinates satisfies `|k_{γ_0}| < Σ_{γ≠γ_0} |k_γ|`. In words: `γ_0` never carries half
of the `ℓ¹`-mass of a linear relation;

(iii) some total order `≺` on the group `G` (`a ≺ b ⇒ a + c ≺ b + c`) has `γ_0 = max(±Γ)`.

*Proof.* (i)⇔(ii). Fix `T ≥ γ_0`. Let `S := ±(Γ_T ∖ {γ_0})`; it is symmetric, so `0 ∈ conv(S)` when
`S ≠ ∅`. If `γ_0 ∈ conv(S ∪ {−γ_0})`, say `γ_0 = −λγ_0 + (1−λ)s` with `s ∈ conv(S)` and `λ < 1`, then
`γ_0 = ((1−λ)/(1+λ))s ∈ conv(S ∪ {0}) = conv(S)`. So `γ_0` is not a vertex iff `γ_0 = Σ_j q_j γ_j` in
`G ⊗ ℝ` with distinct `γ_j ∈ Γ_T ∖ {γ_0}` and `Σ|q_j| ≤ 1`. The set of such `(q_j)` is a polyhedron
defined over `ℚ`, so if it is nonempty it has a rational point. Clearing denominators gives a
relation `Nγ_0 − Σ k_j γ_j = 0` with `N ≥ 1` and `Σ|k_j| ≤ N`. Conversely, such a relation gives
`γ_0 = Σ (k_j/N) γ_j`. A relation with `k_{γ_0} = 0` satisfies the inequality in (ii) unless it is
trivial, and a relation with `k_{γ_0} < 0` can be negated.

(iii)⇒(ii). Let `≺` be such an order and suppose `Nγ_0 = Σ_j k_j γ_j` with `N ≥ 1`, `Σ|k_j| =: K ≤ N`
and `γ_j ≠ γ_0`. Write `γ'_j := sign(k_j)γ_j ∈ ±Γ ∖ {γ_0}`, so that `γ'_j ≺ γ_0`. Then
`Nγ_0 = Σ|k_j|γ'_j ≺ Kγ_0` if `K ≥ 1`, i.e. `(N−K)γ_0 ≺ 0`; if `K = 0` then `Nγ_0 = 0`, impossible. But
`−γ_0 ≺ γ_0`, so `0 ≺ 2γ_0` and `0 ≺ γ_0`, and hence `0 ≼ (N−K)γ_0`. Contradiction.

(i)⇒(iii). For finite `Φ ⊂ ±Γ ∖ {γ_0}`, take `T ≥ max(γ_0, max|Φ|)` and a homomorphism
`u : G_T → ℝ` with `u(γ_0) > u(g)` for all `g ∈ Φ` (possible by (i)). Extend `u` to `ũ : G → ℝ`; this
is possible because `ℝ` is divisible, hence an injective `ℤ`-module. Define `x ≻ 0` iff `ũ(x) > 0`,
or `ũ(x) = 0` and `x > 0` as a real number (recall `G ⊂ ℝ`). This is a total group order with
`γ_0 ≻ g` for all `g ∈ Φ`. The positive cones of total orders on `G` form a closed subset `𝒪` of the
compact space `{0,1}^G`. The subsets `O_Φ := {P ∈ 𝒪 : γ_0 − g ∈ P ∀g ∈ Φ}` are closed and nonempty,
and they have the finite intersection property. So `⋂_Φ O_Φ ≠ ∅`. ∎

**Lemma 3.2 (a bounded direction).** Under V(`γ_0`) there is a homomorphism `φ : G → ℝ` with
`φ(γ_0) = 1` and `|φ(γ)| ≤ 1` for all `γ ∈ Γ`.

*Proof.* For `T ≥ γ_0` put `K_T := {φ ∈ Hom(G_T, ℝ) : φ(γ_0) = 1, |φ(γ)| ≤ 1 ∀γ ∈ Γ_T}`. `K_T` is nonempty
by V, after normalising a separating form. It is compact: `Γ_T` generates `G_T`, so
`φ ↦ (φ(γ))_{γ∈Γ_T}` is a linear injection with closed image, and the constraints bound it. Restriction
maps `K_{T'} → K_T` (`T < T'`) are continuous. An inverse limit of nonempty compact Hausdorff spaces
along `T = T_1 < T_2 < … → ∞` is nonempty. A point of it is a homomorphism on `⋃ G_T = G` with the
required properties. ∎

**Theorem 3.3 (GRH + NCZ + V(`γ_0`)).** For every `M > 0`:

(a) for every `c ∈ ℂ`, `∫_ℝ log|x − c| dν_M(x) = m(F_M − c) ≥ log(2M m_{γ_0}/|M+iγ_0|)`;

(b) for every `c ∈ ℝ` and `0 < ε < 1`,

    ν_M((c−ε, c+ε)) ≤ [ log(|M+iγ_0|/(2M m_{γ_0})) + ½log⁺((1−c)² + V_M) + 1/(2e) ] / log(1/ε);

(c) `ν_M` is absolutely continuous.

Consequently, as in **D** Theorem 3.1, for every `m > −½` the logarithmic densities of
`{D_m > 0}`, `{D_m ≥ 0}` and `{D_m(x) > c x^M/log x}` (every real `c`) exist, and near-ties have
logarithmic density `0`. In particular, under GRH(`χ_{−4}`) and V(`γ_0`) for one ordinate `γ_0` of
`L(s,χ_{−4})`, the logarithmic density of `{x : π(x;4,3) > π(x;4,1)}` exists. For `γ_0 = γ_1 = 6.0209…`,
`M = ½`, and if all zeros are simple (so that `V_{1/2} = σ²_{1/2} = 0.15557`), the constant in (b) at
`c = 0` is `1.7987 + 0.0723 + 0.1839 = 2.0549` [num]. So `ν_{1/2}((−ε,ε)) ≤ 2.055/log(1/ε)`, e.g.
`≤ 0.0892` at `ε = 10^{−10}`.

*Proof.* (a) For `T ≥ γ_0`, Corollary 2.3 with `Γ' = Γ_T` and the vertex `γ_0` gives
`m(F_{M,T} − c) ≥ log(2Mm_{γ_0}/|M+iγ_0|)`, and `m(F_{M,T}∘π_T − c) = m(F_{M,T} − c)` because
`π_{T*}λ_A = λ_{A_T}`. Now `f_T := F_{M,T}∘π_T − c → f := F_M − c` in `L²(A)`. Take a subsequence `T_k` with
`f_{T_k} → f` a.e. Then `log|f_{T_k}| ≤ |f_{T_k}|`, and `|f_{T_k}| → |f|` in `L¹` and a.e. Fatou's lemma
applied to `|f_{T_k}| − log|f_{T_k}| ≥ 0` gives `limsup_k m(f_{T_k}) ≤ ∫ limsup_k log|f_{T_k}| = m(f)`
(reverse Fatou). Hence `m(f) ≥ log(2Mm_{γ_0}/|M+iγ_0|)`. Finally `ν_M` is the law of `F_M`.

(b) Lemma 2.2 with `f = F_M − c`, `‖F_M − c‖_2² = (1−c)² + V_M`.

(c) Let `φ` be as in Lemma 3.2. It is tame in the sense of **D** §3.4 (`|φ(γ)| ≤ 1 ≤ log(γ+3)/log 3`).
Its fibre over `1` contains `γ_0`. Let `C^φ_1 := Σ_{g∈±Γ, φ(g)=1} m_g (2M/(M+ig)) e_g ∈ L²(A)`, and let
`C^φ_{1,T}` be the sum over `g ∈ ±Γ_T`. The support of `C^φ_{1,T}` is a subset of `±Γ_T` that contains
`γ_0`, and `γ_0` is a vertex of `conv(±Γ_T)`, hence of the hull of that subset. By Lemma 2.1,
`m(C^φ_{1,T}) ≥ log(2Mm_{γ_0}/|M+iγ_0|)`. As in (a), `C^φ_{1,T} → C^φ_1` in `L²` and reverse Fatou give
`m(C^φ_1) > −∞`. In particular **`C^φ_1 ≠ 0` a.e.**

Now follow the proof of **D** Theorem 3.5(c). Suppose `ν_M(B) > 0` for a Lebesgue-null Borel set `B`. The
set `S_1` of `θ` for which the real-analytic function `h_θ` of **D** Thm 3.5(a) is constant then has
positive measure, and the argument of **D** Thm 3.5(b) gives `C^φ_v = 0` a.e. on `S_1` for every
`v ≠ 0`. With `v = 1` this contradicts `C^φ_1 ≠ 0` a.e. The finiteness of a fibre is used in
**D** Thm 3.5(c) only to get "nonzero a.e." from **D** Lemma 1.1, and here Jensen supplies it. ∎

**Theorem 3.4 (finite-level version; GRH + NCZ).** Let `μ_T := min{γ ∈ Γ_T : γ is a vertex of
conv(±Γ_T)}`. (The largest element of `Γ_T` is always a vertex: the real embedding `G_T → ℝ` is a
homomorphism that it maximises on `±Γ_T`.) For every `M > 0` and `c ∈ ℝ`,

    ν_M({c}) ≤ 2 liminf_{T→∞} log(|M+iμ_T|/(2M)) / τ_T = 2 liminf_{T→∞} (log⁺ μ_T)/(log T).

In particular, if `log μ_T = o(log T)` along some sequence `T → ∞` (**hypothesis LV**, "low
vertices"), then `ν_M` has no atoms and the conclusions about densities hold. If `μ_T ≤ T^a` for
arbitrarily large `T`, then `ν_M({c}) ≤ 2a`. The first inequality holds for **every** hull function
`1 + Re Σ_{γ∈Γ} w_γ r_γ e_γ` (real `w_γ`, `|w_γ| ≥ 1`, `Σ w_γ²/γ² < ∞`, `V^{>T} > 0` for all `T`). Its
proof uses only Corollary 2.3 and Corollary 4.3, neither of which needs number theory.

*Proof.* Corollary 2.3 gives `δ_T(c) ≤ log(|M+iμ_T|/(2M))` (using `|w_{μ_T}| ≥ 1`), and
Corollary 4.3 below gives the first inequality. For the equality, use
`log(|M+iμ|/(2M)) = log⁺μ + O_M(1)` and `τ_T = log T + O(log log T)` (Lemma 4.0). ∎

**Corollary 3.5 (a special case of Conjecture U; no hypothesis).** For `Y ≥ 1` let
`𝒞_V(M, A_0, s; Y)` be the set of models in `𝒞(M, A_0, s)` in which some ordinate `γ ≤ Y` is a vertex
of `conv(±Γ')`. Put `V*(M, A_0) := 8M²A_0² Σ_{t≥0} log²(t+3)/(M²+t²)`. Then, for `0 < ε < 1` and every
real `c`,

    sup_{𝒞_V(M,A_0,s;Y)} λ{|F_{Γ'} − c| < ε} ≤ [ ½log((M²+Y²)/(4M²)) + ½log⁺((1−c)² + V*(M,A_0)) + 1/(2e) ] / log(1/ε).

So Conjecture U holds on `𝒞_V`, with rate `1/log(1/ε)`.

*Proof.* Corollary 2.3 with `|w_γ| = m_γ ≥ 1` and `γ ≤ Y`, together with `V_{Γ'} ≤ V*`. For the latter,
group the ordinates into the intervals `(t, t+1]`, `t = 0, 1, 2, …`. On each,
`Σ m_γ² ≤ (Σ m_γ)² ≤ A_0² log²(t+3)` and `M² + γ² ≥ M² + t²`. ∎

The class `𝒞_V` is the "easy" part of `𝒞`. It contains every model with one ordinate `≤ Y` outside the
`ℚ`-span of the others, and no rank-one model with top ordinate `> Y` (Remark 3.9). The cyclic models
with a dense spectrum are where the obstruction lives (§5).

**Proposition 3.6 (TS implies V).** If hypothesis TS of **D** §3 holds, then V(`γ_0`) holds for some
`γ_0 ∈ Γ`.

*Proof.* TS gives a nonzero `ℚ`-linear form `φ` on `G ⊗ ℚ` with finite support `Λ := {γ : φ(γ) ≠ 0}`.
Let `K := max_{±Λ} φ > 0` and `Λ_K := {g ∈ ±Λ : φ(g) = K}`, a finite nonempty set. Let `g_0` be a vertex
of `conv(Λ_K)`, exposed by a form `ψ`. For each `T` and small `ε_T > 0`, `φ + ε_Tψ` has the unique
maximiser `g_0` on the finite set `±Γ_T`. Every element outside `±Λ` has `φ = 0 < K`, and elements of
`±Λ ∖ Λ_K` have `φ < K`. Take `γ_0 := |g_0|`; the vertex property is symmetric under `g ↦ −g`. ∎

**Example 3.7 (V without TS, FD, or a finite tame fibre).** Let `b, d > 0` with `b/d ∉ ℚ`, and
`Γ_0 := {b + nd : n ≥ 0}` with any multiplicities `m_n ≥ 1` (bounded, say). Then:
- **V(`b`) holds.** At level `T` let `N_T := max{n : b + nd ≤ T}`. The homomorphism
  `φ_T := b^* − ε d^*` (coordinates in the basis `b, d` of `G = ℤb + ℤd`), with `0 < ε < 1/N_T`, has
  `φ_T(b) = 1` and `φ_T(b+nd) = 1 − εn ∈ (0,1)` for `1 ≤ n ≤ N_T`. The negatives lie in `(−1, 0)`.
- **TS fails.** Two elements `b + n'd`, `b + n''d` above any height span `G ⊗ ℚ`.
- **FD fails** (**D** §3, for any `N_0`). Once `#Γ_T ≥ N_0 + 2`, `Γ_T ∖ Λ_T` contains two elements,
  which span everything. Hence also Devin's Thm 5.4(1), the direct-sum hypothesis and
  self-sufficiency fail.
- **D Thm 3.5(c) does not apply.** A tame `φ` satisfies `|φ(b) + nφ(d)| ≤ C log(b+nd+3)` for all `n`,
  hence `φ(d) = 0`. So every nonzero fibre of a tame `φ` is empty or all of `±Γ_0` (infinite).

For the hull of `Γ_0` (any `M > 0`, any multiplicities `1 ≤ m_n ≤ B`), the proof of Theorem 3.3 applies
verbatim. It uses no number theory, and the convergence step in **D** Thm 3.5(a) needs only
`Σ_γ m_γ²(γ+3)^{1/2}/γ² < ∞`, true here. So the law of the hull
function is absolutely continuous, with the explicit modulus (b). (Directly,
`F − 1 = 2Re(z_1 h(z_2))` with `h(z) = Σ m_n (2M/(M+i(b+nd))) z^n ∈ H²`, `h(0) ≠ 0`. Jensen in `z_2` is
exactly the vertex at `b`.) Numerically [num, `jensen_newton_check.py` (3)]: for `b = γ_1(χ_{−4})`,
`d = 2√2`, `N ≤ 200`, `m(F_N − c) ∈ [−0.1418, −0.0492]` for `c = 0, 2`. With all `m_n = 1`, the computed
`m(F_N − 1)` is `−1.7987 ± 2·10^{−4}` against the bound `−1.798681`. The referee checked that `h_N` has
no zeros in the disc for `N = 4, 16, 64, 200`, so there the bound is attained exactly. The small-ball
bounds of Corollary 2.3 hold with room (e.g. `λ{|F_{200}| < 10^{−2}} = 6.6·10^{−4} ≤ 0.445`).

**Example 3.8 (FD without LV; destruction schedule).** Choose `ℚ`-independent positive reals `a_k`
(`k ≥ 1`) and, for every element `x` introduced, a fresh independent `y_x`. Introduce
`u_x := 2x + y_x`, which satisfies `x = (u_x − y_x)/2`, so `x` is not a vertex once
`u_x, y_x ≤ T`. Take `y_{a_k} := e_k` with `a_k ≈ T_k^{1/3}`, `e_k ≈ T_k/2`, `T_{k+1} = T_k^{3/2}`, and
`y_x ≈ x²` for every other element. Every element is eventually destroyed (each new element gets its
own destroyer), and the counting law holds trivially since the set is sparse. Then:
- **LV fails.** At level `T`, every vertex `x` is undestroyed, `2x + y_x > T`. So `x ≳ T^{1/3}` and
  `log μ_T/log T ≥ 1/3 − o(1)`.
- **FD holds with `N_0 = 2`.** At levels `T ∈ [T_k, T_{k+1})` take `γ_T = a_k ≤ T^{1/3}` and
  `Λ_T = {a_k, u_{a_k}}`. The only other elements whose coordinate on the primitive `a_k` is nonzero
  are destroyers of `u_{a_k}` (and of these), of size `≳ u_{a_k}² ≈ T_k²/4 > T_{k+1}`, hence absent
  at level `T`. So `a_k ∉ ⟨Γ_T ∖ Λ_T⟩_ℚ`.

So LV (and V) is not implied by FD. With Example 3.7, FD and LV are incomparable as conditions on
frequency sets, and V is strictly weaker than TS. (The construction is elementary, but this example
is stated with less detail than the rest of §3.)

**Remark 3.9 (test against the known models with atoms).** In rank one (`G_T` cyclic) the only
vertices of `conv(±Γ_T)` are `±max Γ_T`. So V fails whenever `Γ` is infinite, and LV fails whenever
`max Γ_T ≥ T^{1−o(1)}`. That holds for `ℕ`, `Ωℕ`, spectra with bounded gaps, and the race. Every atom model
of **D** §§7, 9 is of rank one with such a spectrum: Example 9.2 (also rescaled to `10ℕ`),
Propositions 7.1–7.3, Theorem 9.6 (on `ℤ`, `25ℤ`, `50ℤ`). The infinite-rank product model of
`density2/PROGRESS.md` §7 (pieces with full support on `β_kℕ`) also has no low vertex: at level `T` the
vertices are `±β_k⌊T/β_k⌋`, the tops of the pieces, so `μ_T ≥ T/2`. So all of them escape
Theorems 3.3–3.4. They escape through the geometry of `Γ`, **not** through their weights: the theorems
hold for real weights too (Remark 2.4).

*LV in rank one (observed by the referee).* LV can hold in rank one when the spectrum has huge gaps.
For `Γ = {n_k} ⊂ ℕ` with `n_{k+1} = 2^{n_k}` and weights `|w| ≥ 1`, at `T = n_{k+1} − 1` one has
`μ_T = n_k ≈ log₂ T` and `τ_T ≈ 2 log T`. The general form of Theorem 3.4 then gives no atoms, in
accordance with the classical theory of lacunary series (Zygmund; recalled, not used). Hadamard
lacunarity (`n_{k+1} = q n_k`) gives only `ν({c}) ≤ 1`.

**Remark 3.10 (what V asks of the zeros).** V(`γ_1`) fails iff some relation such as
`2γ_1 = γ_a − γ_b` or `3γ_1 = γ_a + γ_b − γ_c` (`γ_1` carrying half the `ℓ¹`-mass) holds. Like TS and FD,
V cannot be verified by a finite computation. It implies that `G` is not cyclic, and even that `G_T`
has rank `≥ 2` for large `T`. That is far beyond what is known (Problem 1 of arXiv:2407.03530 asks
whether two ordinates have irrational ratio). Under LI every ordinate satisfies V.

**Remark 3.11 (another proof of (a)).** Proposition 3.1(iii) puts `F_M − c` in the form `e_{γ_0}·Φ` with
the spectrum of `Φ` in `{g ≼ 0}`. Jensen's inequality for `H¹` of a compact group with ordered dual
(H. Helson, D. Lowdenslager, Acta Math. 99 (1958); new proof by N. Asmar, S. Montgomery-Smith, Bull.
Austral. Math. Soc. 55 (1997), arXiv:math/9503213; statement taken from the abstract, papers not
re-read) then gives (a) directly on `A`. The truncation proof above is self-contained.

**Remark 3.12 (can V be proved for an actual zero under GRH alone? No, with present tools).** Asked
by the parent session; time-boxed.
- *V for one zero already implies an open problem.* V(`γ_0`) implies `γ_0/γ ∉ ℚ` for every ordinate
  `γ > γ_0`, because `γ_0 = (p/q)γ` with `p < q` gives the dominated relation `qγ_0 − pγ = 0`. So V for
  any single zero, in particular for `γ_1`, implies that two ordinates of `L(s,χ)` have an irrational
  ratio. That is Problem 1 of arXiv:2407.03530, open even under GRH. LV implies the same, since a
  vertex below the top ordinate forces `rank G_T ≥ 2`.
- *The explicit formula at finite resolution cannot see dominated relations.* The sibling's Lemma 4.2
  (branch `-density2`, `LINE1.md` §4: no weighted near-recurrence of `Γ_U` before time
  `(2+o(1)) log U`) and the counting law are inequalities about the zero measure at finite
  resolution. A dominated relation with cancellation, e.g. `2γ_1 = γ_a − γ_b` with `γ_a, γ_b` high, can
  be created or destroyed by moving one high zero by an arbitrarily small amount, because differences
  of zeros are dense in `ℝ`. Such a move leaves the counting law intact and changes each sum in
  Lemma 4.2 by `o(1)`, or not at all if the zero lies above `U`. So no argument that uses only these
  inequalities proves V(`γ_0`), for `γ_0 = γ_1` or any other choice. Nor does a dominated relation
  force a recurrence of the whole zero measure, which is the only thing Lemma 4.2 detects.
- *The other suggested choices reduce to the same.* The top ordinate of a truncation is always a
  vertex, but of size `T`, which is useless (Prop 5.1). A zero chosen as the maximum of `±Γ` in a group
  order adapted to the explicit formula is V for that zero (Prop 3.1(iii)), so the first bullet
  applies.
- *What finite computation can do.* A dominated relation `Nγ_1 = Σ k_jγ_j` whose coefficients `k_j`
  all have the same sign involves only zeros below `Nγ_1`. With bounded coefficients these finitely
  many cases can be checked with certified zeros, and excluded unless an exact relation holds. The
  relations with cancellation reach arbitrarily high zeros and cannot be checked this way.

So V stands as a hypothesis of LI type. It is strictly weaker than TS, but unverifiable for the same
reason as all hypotheses of this kind: it contains Problem 1.

---------------------------------------------------------------------------

## 4. Mahler deficits and atoms

Everything in this section holds for every hull function `F = 1 + Re Σ_{γ∈Γ} w_γ r_γ e_γ` with real
weights, `Σ w_γ²/γ² < ∞` and `V^{>T} > 0` for all `T`, in particular for the race under GRH + NCZ (where
`ν_M` is its law); only Lemma 4.0 concerns the race specifically.

**Lemma 4.0 (size of `τ_T` for `L(s,χ)`; the counting estimates used are unconditional).** For fixed
`M`: `log T − 2 log log T − O(1) ≤ τ_T ≤ log T − log log T + O(1)`.

*Proof.* Upper bound on `V^{>T}`: (D3). Lower bound:
`V^{>T} ≥ Σ_{γ∈(T,2T]} m_γ 8M²/(M²+γ²) ≥ (N(2T) − N(T))·8M²/(M²+4T²) ≫_M log T/T` by (RvM). ∎

**Proposition 4.1 (small balls from deficits).** For every `M > 0`, `c ∈ ℝ`, `T`, and `0 < ε < ½`,

    ν_M((c−ε,c+ε)) ≤ [δ_T(c) + C_c] / log(1/(2ε)) + V_M^{>T}/ε²,    C_c := ½log⁺((1−c)²+V_M) + 1/(2e).

With `ε = ε_T := (V_M^{>T} τ_T)^{1/2}` (`→ 0`) this reads
`ν_M((c−ε_T,c+ε_T)) ≤ 2(δ_T(c) + C_c)/(τ_T − log τ_T − 2 log 2) + 1/τ_T`.

*Proof.* (D2), and Lemma 2.2 applied on `A_T` to `F_{M,T} − c`, with `‖F_{M,T} − c‖² ≤ (1−c)² + V_M`. For
`ε = ε_T`: `V^{>T}/ε_T² = 1/τ_T` and `log(1/(2ε_T)) = ½(τ_T − log τ_T) − log 2`. ∎

**Proposition 4.2 (an atom forces a deficit).** If `ν_M({c}) = η > 0`, then for **every** `T`

    δ_T(c) ≥ (η/2) τ_T − ½log⁺((1−c)² + V_M) − 1/e.

*Proof.* Let `S := {F_M = c}`, `λ_A(S) = η`, and `R_T := F_M − F_{M,T}∘π_T`. On `S`,
`F_{M,T}∘π_T − c = −R_T`. Concavity of `log` on `S` and on `S^c`, as in Lemma 2.2, gives
`m(F_{M,T} − c) ≤ (η/2)log(∫_S|R_T|²/η) + ½log⁺((1−c)²+V_M) + 1/(2e)`. Then use `∫_S|R_T|² ≤ V_M^{>T}`
and `(η/2)log(1/η) ≤ 1/(2e)`. ∎

**Corollary 4.3 (deficit criterion).** `ν_M({c}) ≤ 2κ(c)`, where `κ(c) := liminf_{T→∞} δ_T(c)/τ_T`.
In particular `δ_T(c) = o(τ_T)` along some sequence excludes an atom at `c`.

*Proof.* If `ν_M({c}) = η > 0`, divide Proposition 4.2 by `τ_T → ∞` and take the liminf. If `η = 0`
there is nothing to prove; note that `κ ≥ 0`, since `δ_T(c) ≥ −½log⁺((1−c)²+V_M)` by concavity. ∎
(Proposition 4.1 gives the same conclusion and, in addition, a small-ball bound at the scale `ε_T`.)

**Remark 4.4 (the exponent).** So Propositions 4.1–4.2 contain one inequality. Proposition 4.2 holds
for every `T`, Proposition 4.1 at the scale `ε_T`. The factor `2` is not known to be sharp: for
Example 9.2 the atom is `1/3` and `κ ≈ 0.32` [num], so there `η ≈ κ`. What is known about `κ`:

| situation | bound for `δ_T` | `κ` | conclusion |
|---|---|---|---|
| any model, vertex `μ_T` (Cor 2.3) | `≤ log(|M+iμ_T|/(2M))` | `≤ liminf log(|M+iμ_T|/(2M))/τ_T` | Thm 3.4 |
| V(`γ_0`) | `≤ log(|M+iγ_0|/(2M))`, bounded | `0` | Thm 3.3 |
| rank one with `max Γ_T ≥ T^{1−o(1)}` (Lemma 2.1 only) | `≤ log(|M+iT|/(2M))` | `≤ 1` | `ν({c}) ≤ 2`: void |
| Example 9.2 (real weights, atom `1/3`) [num] | `≈ 0.32 log N` | `≈ 1/3` | consistent: `1/3 ≤ 2/3` |
| Theorem 9.6 (integer, atom `> 0.1`) | `≥ 0.05 log T − O(log log T)` (Prop 5.2) | `≥ 0.05` | — |
| integer near-atoms, `N ≤ 512` [num] | `0.7 – 1.0` | — | — |
| LI model on the true zeros, `T ≤ 640` [num] | `0.08 – 0.105` (stable) | `0` | — |

For the race the target is `κ(0) = 0`. Lemma 2.1 alone gives `κ ≤ 1` in rank one with a dense
spectrum, and integer models without the counting law can have `κ ≥ 0.05`.

---------------------------------------------------------------------------

## 5. The rank-one wall: why vertex-type Jensen arguments cannot give Conjecture U

**Proposition 5.1 (rank one).** If `G_T` is cyclic, the vertices of `conv(±Γ_T)` are `±γ_max(T)`. So
Corollary 2.3 gives exactly `δ_T(c) ≤ log(|M+iγ_max|/(2M m_{γ_max}))`. If moreover
`γ_max(T) ≥ T^{1−o(1)}` (`ℕ`, `Ωℕ`, bounded gaps; for the race this is automatic), Corollary 4.3 with
Lemma 4.0 gives only `ν_M({c}) ≤ 2`. *Proof.* `conv(±Γ_T)` is the segment `[−γ_max, γ_max]` of
`G_T ⊗ ℝ ≅ ℝ`. ∎

**Proposition 5.2 (the loss is real for integer models; no hypothesis).** The models `F` of **D**
Theorem 9.6 (integer `m_n ≥ 1`, `m_n ≪ (log n)^{5/2}`, `G = ℤ`, `F = 0` a.e. on a set of normalised
measure `η > 0.1`) satisfy, for their truncations `F_N`,

    δ_N(0) = −m(F_N) ≥ (η/2)(log N − 5 log log N) − O(1) ≥ 0.05 log N − 0.25 log log N − O(1).

*Proof.* `V^{>N} = Σ_{n>N} m_n² 8M²/(M²+n²) ≪ (log N)^5/N`, so `log(1/V^{>N}) ≥ log N − 5 log log N − O(1)`.
The proof of Proposition 4.2 uses only that `F = 0` on a set of measure `η`, that
`‖F − F_N‖² = V^{>N}` and that `V := ‖F − 1‖² < ∞`; it applies verbatim on `𝕋`. ∎

(The implied constants are astronomically large, since the multiplicities of Theorem 9.6 reach
`≈ 3·10^5` below `4^7`. So this is an asymptotic statement, not one visible at accessible `N`.)

**Remark 5.3 (what Lemmas 2.1–2.2 cannot do).** The no-go below is limited to the following kind of
argument: Lemma 2.1 applied to `F − c`, or to any polynomial with the same vertex data, together with
Lemma 2.2. Lemma 2.1 could also be applied to auxiliary polynomials built from `F` (products of
translates, compositions with endomorphisms of `A`, resultants). Nothing here excludes such a use.
Such an argument sees a finite model only through `|w_{γ_0}|` at a vertex and through `V_{Γ'}`. The
latter enters via `½log⁺(|1−c|² + V_{Γ'})`. In rank one the vertex is the top ordinate, and two classes
supply the same data as the truncations of the true race:

- *Real weights* (`A_0 ≥ 2`): the class `𝒞^ℝ` obtained by replacing "integer `m_γ ≥ 1`" with "real
  `w_γ ≥ ½`, top weight `≥ 1`". It contains the truncations of the rescaled Example 9.2: spectrum
  `10ℕ`, weights `4.868` (index `≢ 0 mod 30`) and `½`, `σ² = 0.1439`, counting law with `A_0 ≈ 1.96`;
  truncate at an index `≢ 0 mod 30`. These truncations satisfy
  `λ{|F_T| < K‖F − F_T‖_2} ≥ 1/3 − 1/K²`, because the limit vanishes on a set of measure `1/3`. So
  `ω^ℝ(ε) ≥ 1/3` for every `ε` (sharpness remark after **D** Lemma 9.13), and U is **false** on `𝒞^ℝ`.
  Every model of `𝒞` and of `𝒞^ℝ` has `V ≤ V*(½, A_0)` (proof of Cor 3.5), so for rank-one members the
  bound of Corollary 2.3 is the same. In higher rank a non-top vertex of a model in `𝒞^ℝ` may carry
  weight `½`, which costs at most `log 2`.
- *Integers without the counting law*: the truncations of **D** Theorem 9.6 (top multiplicity `≥ 1`,
  `V ≤ ‖F − 1‖² < ∞`), with near-atoms of measure `> 0.1`. Corollary 2.3 differs from its form on `𝒞`
  only by the constant `½log(1 + V)`. Whether these models fit the `σ²`-bound `s = 0.1556` of U (after
  rescaling to `50ℕ`) rests on **D**'s uncertified estimate `σ² ≲ 0.125`. They fit `s = 0.5` on `25ℕ`.
- *Signed or nonnegative integers with the counting law* (heuristic only): Remark 9.21 of the sibling
  session's `DENSITY.md` (branch `-density`; heuristic, not a theorem) suggests that noise-shaped rounding
  gives atoms for integer multiplicities `≪ (log n)^θ`, any `θ > 0`, on `G = ℤ`. If that is right, U fails
  in rank one even with positivity and the counting law.

So an argument of the limited kind above holds verbatim on `𝒞^ℝ`, where U fails. It uses neither
integrality beyond `|w| ≥ 1` at one point, nor positivity, nor the counting law, so it cannot prove U.
This is the brief's test: the real-weight and `(log)^{5/2}` models do **not** escape such an argument,
so it cannot be a proof of U. Theorem 3.3 claims less, and those models do escape it, through the
geometry of their spectra (Remark 3.9).

**Remark 5.4 (faces do not help).** Let `τ` be a face of `conv(supp(F_{M,T} − c))`, exposed by
`φ : G_T → ℤ` with `φ ≡ K` on `τ` and `φ < K` on the rest of the support. Along `j_φ` the top
coefficient of `s ↦ F_{M,T}(θ j_φ(s)) − c` is `F_τ(θ)`, where `F_τ := Σ_{g∈τ} c_g e_g`. Jensen as in
Lemma 2.1 gives `m(F_{M,T} − c) ≥ m(F_τ)`. Every face contains a vertex, and Lemma 2.1 applied to `F_τ`
returns only the vertex bound. To beat it, one needs `m(F_τ) ≥ −o(log T)` for a face without low
vertices.

- For a *high* face (all ordinates of `τ` are `≥ T^{1−o(1)}`) this is impossible. By concavity,
  `m(F_τ) ≤ ½log‖F_τ‖_2² = ½log Σ_{g∈τ} m_g² 4M²/(M²+g²) ≤ −½log T + o(log T)`, since
  `Σ_{γ≤T} m_γ² ≤ A_0 log(T+3)·N(T) ≪ T log² T`. So high faces give at best `δ_T ≤ ½log T + o(log T)`,
  hence only `ν({c}) ≤ 1` by Corollary 4.3.
- A face with low ordinates in its relative interior but only high vertices (for instance the top
  edge of a `{(n,0), (n,±1)}`-type configuration in `G_T ⊗ ℝ ≅ ℝ²`) has a face polynomial `F_τ` that is
  itself a finite model of lower rank. It has no constant term, and its coefficients `≍ m_g/|g|` have a
  logarithmically divergent sum, so no single term dominates. Bounding `m(F_τ)` from below is then a
  problem of the same type as the original one.

### 5.1 Question (1) of the brief: small values of `F` versus the integer polynomial `2w`

On a finite model, `(D + M)F_{Γ'} = 2M w_{Γ'}` with `2w_{Γ'} = 1 + 2Σ_{g∈±Γ'} m_g e_g` an **integer**
Laurent polynomial. Its vertex coefficients are `2m_γ ≥ 2`, so `m(2w_{Γ'}) ≥ log 2` (Lemma 2.1), and by
Lemma 2.2 `λ{|2w_{Γ'}| < δ} ≤ (−log 2 + ½log‖2w_{Γ'}‖_2² + 1/(2e))/log(1/δ)`. This controls only the set
where `|w| < δ` with `log(1/δ) ≫ log‖2w‖_2 ≈ ½ log N(T)`. The near-atom sets of `F` do not lie there:

- **Exactly on the hull** the transfer "`F = 0` on `Z` ⇒ `DF = 0` at density points ⇒ `w = 0` on `Z`"
  needs `F` to be differentiable along the flow, i.e. `w` locally integrable along orbits. That is
  false for the true race: `w` is not even a measure (`density2/TOY.md` §1, Remark (2)).
- **On finite models**, `F` and `DF` are uncorrelated: `∫ F·DF dλ = ½∫ D(F²) dλ = 0`. Small values of
  `F` carry derivative `DF = 2Mw − MF` of typical size. [num, `rank_one_deficit.py`] In the LI model
  on the true zeros, the rms of `w_T` on `{|F_T| < 0.02}` agrees with `‖w_T‖_2 = √(2n)` to within `3%`
  for `T ≥ 200` (`32.3` vs `32.0` at `T = 640`; `8%` at `T = 50`). On the near-zero sets of the
  structured near-atoms it is `O(1)`, but not small: `2.3–3.2` (Example 9.2), `4.5`, `7.2` (beam-search
  0/1 near-atoms).
- **Averaging along the flow** does give smallness: `(1/h)∫_0^h 2Mw(τ_uθ)du = (F(τ_hθ) − F(θ))/h + M·avg F`,
  so `|avg w| ≤ ε/(Mh) + ε/2` on flow segments inside `{|F| < ε}`. But the averaged `w` has coefficients
  `m_γ(e^{iγh}−1)/(iγh)`, which are not integers. And it is small only for `h ≫ ε`, while `{|F| < ε}`
  need contain no flow segment of fixed length; for the hull, level sets have empty interior
  (**D** Prop 4.3).
- **The "positive kernel" representation** `F = ∫_0^∞ 2Me^{−Mu} w(τ_{−u}·) du` does not combine with Jensen:
  the translates `τ_{−u}w` have coefficients `m_γ e^{−iγu}`, not integers, and `log|∫Kf|` is not bounded
  below by `∫ K log|f|`. The near-zeros of `F` are produced by cancellation in this average.

So the correct "Mahler-type lower bound for `log|F_T|` with integer-times-analytic coefficients" is
Lemma 2.1 applied to `F_T` itself (Corollary 2.3). In rank one it costs exactly `log T` (Prop 5.1).

### 5.2 Question (2): Kronecker rigidity and a cyclotomic dichotomy

- `2w_{Γ'} ≡ 1 (mod 2)` has vertex coefficients `2m_γ ≥ 2`, so `M(2w_{Γ'}) ≥ 2 > 1` always. It is never of
  Kronecker type (`M = 1`: `±` monomial × cyclotomic factors `Φ_n(z^a)`; Kronecker 1857, and in several
  variables Boyd, Smyth 1981; recalled, not used). The Lehmer-type lower bounds (Smyth 1971:
  `M ≥ 1.3247…` for non-reciprocal; Dobrowolski 1979: `M(α) ≥ 1 + c(log log d/log d)³`; recalled, not
  used) are weaker than this trivial bound. The extremal value `M(2w) = 2` does occur: `m ≡ 1` on `25ℕ`
  gives `m(2w_N) = log 2` [num], all zeros on the circle (the referee confirmed this with certified
  roots for `N = 500, 5000`). But the polynomial is not monic, so Kronecker's theorem says nothing about
  roots of unity. Bounds in terms of the number of terms (Dobrowolski–Lawton–Schinzel: `M(P) ≥ 1 + c(k)`
  for a non-cyclotomic `P` with `k` terms; recalled, not used) are likewise constants beaten by the
  vertex bound.
- `F_{Γ'}` has non-integer coefficients, and the arithmetic theory does not apply to it.
- A dichotomy for `w` (flat, `m(2w) ≈ log‖2w‖_2`, versus concentrated) exists, but it does not transfer
  to `F`. [num] Concentrated `w` with an atom (Example 9.2: `m(2w) = 0.29` vs `log‖2w‖ = 5.7`); flat `w`
  without one (LI model); concentrated `w` without one (`25ℕ`, Toeplitz); moderately flat `w` with a
  near-atom (beam: `2.25` vs `3.4`). The real-weight Proposition 7.3 has `m(2w) = −2.55`, which is
  impossible for integers, yet its `F` is no smaller on fat sets than allowed for integer models. *Sketch
  (not proved):* even in the flat branch the transfer back to `F` is by an `L⁰` Bernstein inequality
  (Arestov's theorem, Izv. AN SSSR 45 (1981); recalled, not checked in this generality). It loses
  `log T`, which would leave `δ_T ≤ ½ log T + O(1)` and hence only `ν({c}) ≤ 1` from Corollary 4.3.

### 5.3 Question (3): the exponent

Answered by Propositions 4.1–4.2, Corollary 4.3 and Remark 4.4. Lemma 9.13 needs precision `ε_T ≈ T^{−1/2+o(1)}`
(tail `V^{>T}` between `≫ log T/T` and `≪ log² T/T`, Lemma 4.0), and Lemma 2.1 gives
`λ{|F_T| < ε} ≤ (log T + O(1))/log(1/ε)` in rank one: about `2` at `ε = T^{−1/2}`. Jensen is useful
only at precision `T^{−A}` with `A → ∞`. What is needed is `δ_T = o(log T)`, exponent `0` in place of `1`.
An atom of mass `η` forces exponent `≥ η/2`.

### 5.4 Question (4): Erdős–Turán, Bilu, zeros near the circle

In rank one put `P_N(z) := z^N F_N(z)`, a self-inversive polynomial of degree `2N` (`F_N` is real on the
circle). By (J), `δ_N = log(1/|C_N|) − Σ_{|ζ|>1} log|ζ|`, with `C_N` the top coefficient. A statement like
`δ_N = o(log N)` therefore needs a **lower** bound for `Σ_{|ζ|>1} log|ζ|`: zeros must stay away from the
circle. The Erdős–Turán theorem (P. Erdős, P. Turán, Ann. of Math. 51 (1950)) bounds the angular
discrepancy of the zeros by `O(√(n log(L/√|a_0a_n|)))`, `L` the sum of the moduli of the coefficients.
Here `log(L/|C_N|) ≈ log N`. The radial count, at most `O(log(L/√|a_0a_n|)/δ)` zeros with `||ζ| − 1| ≥ δ`,
follows from Jensen's formula itself (a source is Hughes–Nikeghbali, Compositio Math. 2008, as
recalled by the referee; not checked). Both bound from **above** how many zeros lie away from the
circle. That is the wrong direction: zeros near the circle make `δ_N` close to its maximum
`log(1/|C_N|)`. The discrepancy `√(N log N)` is in any case far coarser than the `o(log N)` needed.
Bilu's theorem (Yu. Bilu, Duke Math. J. 89 (1997)) is qualitative: Galois orbits of algebraic numbers
of small height equidistribute on the circle. It concerns integer polynomials. For `2w_N` it applies,
but by §5.1 `w_N` is not small on near-atom sets.

### 5.5 Siegel's lemma, Bombieri–Vaaler, counting

These produce small integer solutions, i.e. they work in the counterexample direction. The
corresponding count is the entropy heuristic of `density2/TOY.md` §5. With multiplicities bounded by
`W`, near-atoms of measure `η` at precision `τ` should need `η log(1/τ) ≲ log W`. *Observation:* with
growing alphabets this count is violated by **D** Theorem 9.6. There `W ≍ (log N)^{5/2}` and
`τ ≍ (log N)^{5/2}N^{−1/2}`, so the count predicts `η ≲ 5 log log N/log N → 0`, yet `η > 0.1`. This only
happens for astronomically large `N`, because of the constants. So counting cannot be the mechanism
behind U: structured cascades beat it. (Not a contradiction with `density2`'s Conjecture 5.1, which is
for fixed `W`.)

### 5.6 Two related constructions elsewhere (sibling branch `-density`, `DENSITY.md` §9)

- *Certificate duality* (sibling's Remark 9.20, refereed there as correct, standard, weak). If `F = c`
  a.e. on `S`, then `dist(κφ_0, ℤ) ≤ (2M)^{−1}‖F − c‖‖(M − D)φ‖_{L²(A∖S)}` for every trigonometric
  polynomial `φ` with `φ_0 ∈ ℤ` and `φ_γ + φ_{−γ} ∈ ℤ` (`κ = (c−1)/2`). This is the dual form of integrality
  (Fekete). With `φ = 1` it gives only `λ(S) ≤ V/(1+V)`; it recovers arcs through Fekete's theorem; in
  the narrow class (all `φ_g ∈ ℤ`) no certificate exists for an arc in a half-circle; for fat nowhere
  dense `S` it is open. From the side of this paper: a certificate needs `(M − D)φ` small off `S`, i.e.
  `φ ≈ C_j e^{Mt}` on the components of `A ∖ S`, not `φ` small. Lemmas 2.1–2.2 control where an integer
  polynomial is small, so they do not constrain certificates. Not pursued.
- *Noise shaping* (sibling's Remark 9.21, heuristic). Its integer noise-transfer filters `H` with
  `h_0 = 1` are limited exactly by Jensen. Since `m(H) ≥ log|h_0| = 0`, `|H| ≤ τ` on an arc of normalised
  length `2/q` forces `(2/q)log(1/τ) ≤ (1 − 2/q)log‖H‖_∞`. That is `log(1/τ)/log‖H‖_∞ ≤ (q−2)/2`, the
  bound `β ≤ (q−2)/2` quoted there. So in that construction Jensen's formula is the limiting
  constraint, not an obstruction: integer polynomials *can* be tiny on short arcs, at a polynomial
  cost in height. This is the same polynomial-precision phenomenon as in §5.3.

---------------------------------------------------------------------------

## 6. Numerical checks [num]

All at `M = ½`. Scripts in this directory; outputs `*_results.txt`.

**`jensen_newton_check.py`**: Lemma 2.1 (24 random multivariable polynomials, min slack `−0.0000`),
Lemma 2.2 (max defect `−0.19`), Example 3.7 (table below), the constant of Theorem 3.3(b)
(`2.0549`).

| Example 3.7, `N` | 4 | 16 | 64 | 200 |
|---|---|---|---|---|
| `m(F_N)` (`m_n = 1`) | `−0.084` | `−0.059` | `−0.058` | `−0.059` |
| `m(F_N − 1)` (`m_n = 1`; bound `−1.7987`) | `−1.7986` | `−1.7987` | `−1.7988` | `−1.7985` |
| `m(F_N − 1)` (`m_n ∈ {1,2,3}` random) | `−1.167` | `−1.681` | `−1.123` | `−1.066` |

**`rank_one_deficit.py`** (cyclic models; `δ_N = −m(F_N)`, `J_N = log(1/|C_top|)` = Jensen's upper bound
for `δ_N`):

| model | `N` | `δ_N` | `J_N` | `δ_N/log N` | `m(2w)` | `log‖2w‖₂` | rms `w` on near-zero set / `‖w‖₂` |
|---|---|---|---|---|---|---|---|
| `1` on `25ℕ` | 50000 | 0.003 | 10.82 | 0.000 | 0.693 | 4.84 | 74.5 / 63.3 |
| random 0/1 | 8192 | 0.435 | 9.01 | 0.048 | 4.23 | 5.20 | 76.6 / 90.9 |
| Toeplitz base 5 | 15625 | 0.125 | 9.66 | 0.013 | 2.05 | 5.52 | 105.5 / 125.0 |
| Ex 9.2 (real) | 30, 300, 3000, 30000 | 1.06, 1.81, 2.57, 3.34 | 4.09 … 11.00 | 0.310 → 0.324 | 0.288 | 2.26 → 5.71 | 2.3–3.2 / 4.8 → 151 |
| beam 0/1, `η = 0.15` | 512 | 0.962 | 6.21 | 0.154 | 2.25 | 3.42 | 4.49 / 15.3 |
| beam 0/1, `η = 0.30` | 512 | 0.696 | 6.24 | 0.111 | 2.22 | 3.28 | 7.22 / 13.3 |
| Prop 7.3 (real) | 100 | 2.607* | 8.98 | 0.566 | −2.549 | 2.11 | 0.000 / 4.13 |

*For Proposition 7.3 the grid value `2.617` is biased, because `|F| ≈ 10^{−15}` on the arc `J`. The value `2.607052` comes from certified roots (python-flint) with the exact coefficients `C_n = c_n(½ − in)`; the script prints it, and the referee obtained the same value by split mpmath quadrature. Rounding the weights `w_n` to double alone moves it to `2.6133`.

The beam near-atoms come from `density2/toy_beam2.py` (width 150, `N = 512`; log `beam_search_log.txt`,
sequences `beam_best_eta*.npy`). LI model on the 511 stored zeros of `L(s,χ_{−4})` below `640`
(`prototype/explore/chi4_zeros.txt`), Monte Carlo with `10^6` points: `m(F_T)` is
`−0.083, −0.093, −0.100, −0.104, −0.105` at `T = 50, 100, 200, 400, 640`, `m(F_T − 1) ≈ −1.55` (vertex
bound `−1.80`), and rms `w_T` on `{|F_T| < 0.02}` is `5.8, 9.7, 16.0, 23.8, 32.3` against
`‖w_T‖₂ = 6.3, 10.0, 15.6, 24.1, 32.0`.

Readings: for every integer model tried the deficit stays `≤ 1` while Jensen's bound grows like
`log N`. The real-weight atom model Example 9.2 has deficit growing like `(1/3) log N`. That is more
than the `(1/6) log N` that Proposition 4.2 requires for an atom of mass `1/3`, because its tail is
`O(1/N)` pointwise on the arc rather than `N^{−1/2}`. Its `w` is bounded on the arc. The real-weight
Proposition 7.3 visibly violates the integrality constraint `m(2w) ≥ log 2`.

---------------------------------------------------------------------------

## 7. What remains, in this language

**Conjecture L (Mahler form; a sufficient statement, open).** There is `ψ(T) = o(log T)` such that
every model in `𝒞(½, A_0, 0.1556)` with ordinates `≤ T` satisfies `−m(F_{Γ'}) ≤ ψ(T)`.

**Proposition 7.1.** Conjecture L implies Conjecture U. (The converse is not shown: U could hold while
some models have super-small values on thin sets.)

*Proof.* Let `F ∈ 𝒞` and `0 < ε < 1/2`. Put `T := ε^{−3}` and let `F_T` be the truncation of `F` to its
ordinates `≤ T`. It is again in `𝒞`, since a sub-multiset satisfies the counting law and the
`σ²`-bound. As in (D2),
`λ{|F| < ε} ≤ λ{|F_T| < 2ε} + V^{>T}(F)/ε²`. By the grouping in the proof of Corollary 3.5,
`V^{>T}(F) ≤ 8M²A_0² Σ_{t≥⌊T⌋} log²(t+3)/(M²+t²) ≪ A_0² log² T/T` uniformly over `𝒞`, so the second term
is `≪ ε log²(1/ε)`. By Lemma 2.2 the first is `≤ (ψ(ε^{−3}) + ½log(1+V*) + 1/(2e))/log(1/(2ε))`.
Both tend to `0` uniformly. ∎

**Status of L.**
- It holds on `𝒞_V(Y)` with `ψ = O(1)` (Corollary 3.5).
- Its analogue without the counting law fails. The integer models of **D** Theorem 9.6 have
  `ψ(T) ≥ 0.05 log T − O(log log T)` (Proposition 5.2): on `25ℕ` with `σ² < 0.5`, i.e. for `s = 0.5`, and
  on `50ℕ` with `σ² ≲ 0.125` if **D**'s uncertified estimate holds.
- In rank one with a dense spectrum, Lemma 2.1 gives `ψ(T) = log T + O(1)`. Remark 5.3 describes the
  arguments that cannot do better.
- In rank one, by §5.4, L says: **for integer cyclic models with the counting law,
  `Σ_{|ζ|>1} log|ζ| ≥ log(1/|C_N|) − ψ(N)` for the zeros of `z^N F_N`.** Not all zeros may hug the
  circle. That is a statement about the complex zeros of `F_N` in which integrality and the counting
  law would have to work together. No result of the theory of Mahler measure that I know of
  (Kronecker, Lehmer, Dobrowolski, Bilu, Erdős–Turán) addresses it.
- *Caveat.* The sibling's noise-shaping heuristic (§5.6, Remark 5.3) suggests that U, and hence L, may be
  **false** in rank one even with positivity and the counting law. The bounded-multiplicity case, its
  case (b), is untouched by that heuristic.

**Heuristic 7.2 (where the deficit sits in rank one; not a proof).** By (J),
`δ_N = log(1/|C_N|) − Σ_{|ζ|>1} log|ζ|` for `P_N = z^N F_N`, which has `N` outer zeros counted with
pairing. For integer models the outer zeros sit at `|ζ| − 1 ≈ (0.8–0.94)(log N)/N` (median), so
`Σ log|ζ| ≈ log N` and `δ_N = O(1)` [num, `zeros_heights.py`: random 0/1, all ones; `N ≤ 1000`]. An atom of
mass `η` forces the outer zeros towards the circle on average (Proposition 4.2): `Σ log|ζ|` drops by at
least `(η/2) log N − O(1)`. The beam near-atom of §6 has `82` of its `576` zeros **on** the circle and `δ = 0.95`.
Erdős–Turán-type results see zeros at the scale `(log N)/N`, but not how near they are.

**Numerical check 7.3 [num] (maximising the deficit).** Beam search (width 40) maximising
`δ_N = −m(F_N)` over 0/1 weights at `M = ½` (`deficit_beam.py` → `deficit_beam_results.txt`):

| `N` | 8 | 64 | 256 | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|---|---|---|
| `δ_N`, all frequencies (`σ² ≈ 0.74`) | 0.844 | 0.924 | 1.012 | 1.032 | 1.058 | 1.064 | 1.071 |
| Jensen bound `J_N` | 2.08 | 4.14 | 5.11 | 6.16 | 6.85 | 7.56 | 8.29 |
| `δ_N`, frequencies `≥ 13` (`σ² ≤ 0.1598`) | 0.000 | 0.081 | 0.088 | 0.090 | 0.091 | 0.092 | 0.093 |

The optimised deficits stay bounded or nearly so (`δ_N/log N` falls from `0.41` to `0.13`), while
Jensen's bound grows like `log N`. This is weak evidence for L in the bounded case. A greedy beam
cannot find structured cascades such as Theorem 9.6 or the sibling's noise shaping, which only
develop at enormous `N`.

**Status of this line.** Killed as a route to Conjecture U, for the reasons in §5. Kept:
- the Jensen–Newton bound, a small-ball tool for any configuration with a low vertex;
- Theorem 3.3: hypothesis V, a new LI-type sufficient condition, strictly weaker than TS, giving
  absolute continuity and an explicit `log`-modulus;
- Theorem 3.4: the finite-level version, which also covers super-lacunary spectra of rank one;
- Corollary 3.5: U on `𝒞_V`;
- the deficit criterion (Props 4.1–4.2, Cor 4.3).

**Novelty.** Searches by the referee and by me found no prior use of hypothesis V, or of a Jensen or
Newton-polytope criterion, for continuity of prime-race limiting distributions. The weakened-LI
literature uses self-sufficiency and direct-sum hypotheses (Devin, arXiv:1706.06394, 1909.03975;
Martin–Ng, arXiv:1710.00088). Related tools: the Helson–Lowdenslager Jensen inequality (Remark 3.11),
and ordered-lattice lower bounds for Mahler measure (S. Akhtari, J. Vaaler, arXiv:1810.12413; not
read). What is new here is the application, and V ⇒ absolute continuity via **D** Thm 3.5. The tools
are not new.

## 8. Files and review

`jensen_newton_check.py` → `jensen_newton_check_results.txt` [num] · `rank_one_deficit.py` →
`rank_one_deficit_results.txt` [num] · `deficit_beam.py` → `deficit_beam_results.txt`,
`deficit_best_N4096_W40_n1.npy` [num] · `zeros_heights.py` → `zeros_heights_results.txt` [num] ·
`beam_best_eta0.15_N512.npy`, `beam_best_eta0.3_N512.npy`, `beam_search_log.txt` (output of
`density2/toy_beam2.py`, used as integer near-atoms) · `PROGRESS.md` (log, review details).

*Review.* One adversarial read by a fresh `opus` subagent (2026-11-27), of the version without §5.4's
corrections, Remark 5.4, §5.6, Heuristic 7.2, Numerical check 7.3:
- Verdict: **repairable issues only**. Every Lemma, Proposition, Theorem and Corollary was checked
  line by line, including Theorem 3.3(c) for a bounded direction with an infinite fibre.
- No theorem as stated excludes an atom in a known atom model.
- The two scripts were reproduced byte for byte, with independent checks (certified roots, new random
  polynomials).
- Applied: §4 restructured (Prop 4.2 implies the atom form of 4.1; "exact exponent" and "conversely"
  dropped; the factor 2 is not known to be sharp); the "cannot prove U" claims limited to arguments that
  use vertex data only, with L called a sufficient statement; rank-one claims qualified by
  `max Γ_T ≥ T^{1−o(1)}`, with the super-lacunary example (LV in rank one) added; Bilu's theorem
  correctly attributed; Prop 7.3 deficit corrected to `2.607`; cosmetic fixes.
- Remark 5.4, §5.6, Heuristic 7.2 and Numerical check 7.3 were written after the read and have not
  been refereed.
